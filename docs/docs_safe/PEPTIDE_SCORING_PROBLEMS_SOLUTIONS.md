# Peptide Discovery Scoring: Problems and Solutions

Date: 2026-03-25

## Context

Our naive `discover_peptides.py` used weighted-sum scoring on single MS1 scans. The user identified 4 real problems by visually inspecting the evidence viewer, where the scored isotope envelopes clearly did not match what a human would accept. This document records those problems and the solutions researched from established proteomics tools.

---

## Problem 1: M-1/M-2 Peaks (False Monoisotopic Assignment)

### What the user saw

ADDTWEPFASGKT had peaks at M-1 and M-2 (one and two isotope spacings *below* the putative monoisotopic peak). This means the matched envelope likely belongs to a different, heavier ion whose higher isotope peaks happen to overlap with our candidate's expected positions.

### Why it's a problem

If significant intensity exists below the monoisotopic m/z, the envelope is not self-consistent. The "M+0" we scored is actually someone else's M+2 or M+3. Our naive scorer had no concept of this -- it only looked rightward from M+0 and never checked leftward. This means we could assign a high confidence score to a completely wrong ion.

### Solution: ms_deisotope left-search + averagine theoretical distribution

ms_deisotope addresses this via the `left_search_limit` parameter in its `charge_state_determination()` method. Rather than explicitly penalizing M-1 peaks, it takes a more principled approach:

1. For each experimental peak, the algorithm searches *leftward* (to lower m/z) by `left_search_limit` isotope spacings to find candidate monoisotopic peaks that might be below the current peak.
2. If a peak exists at M-1, it generates a candidate fit where that lower peak is the true monoisotopic and the current peak is its M+1.
3. Both hypotheses (current peak = M+0, and lower peak = M+0) are scored against averagine theoretical distributions.
4. The peak dependency graph resolves the conflict: each experimental peak is assigned to exactly one isotopic envelope, and the assignment with the best score wins.

### Implementation

```python
deconvoluted_peaks, _ = ms_deisotope.deconvolute_peaks(
    peaks,
    averagine=ms_deisotope.peptide,
    scorer=ms_deisotope.PenalizedMSDeconVFitter(10., 1.0),
    left_search_limit=1,   # check 1 position below each peak
)
```

The `_find_previous_putative_peak()` method in `base.py` handles the leftward search. The `AveraginePeakDependenceGraphDeconvoluter` builds a hypergraph of all competing fits and resolves them via greedy maximization over disjoint subgraphs.

Alternative approaches from the literature: THRASH (Horn et al., 2000) uses iterative subtraction -- it processes the most intense peak first, assigns its full isotopic envelope, then *deletes* those peaks from the spectrum. By the time it reaches our candidate region, the contaminating peaks from the heavier ion are already gone.

---

## Problem 2: Peaks Between Isotopes (Charge State Ambiguity)

### What the user saw

AEVVFTANDSGPRR had peaks between M+0 and M+1 at the assumed z=2 spacing (~0.5 Da). These interstitial peaks at ~0.25 Da spacing suggest the true charge state is z=4, not z=2. At z=2, isotope peaks are spaced at 1.003/2 = 0.5015 Da. At z=4, they are spaced at 1.003/4 = 0.2508 Da. A z=4 ion places peaks at positions that a z=2 model sees as "between" its expected isotope positions.

### Why it's a problem

Our naive scorer assumed a fixed charge state (from the database) without verifying it against the spectral evidence. If the true charge state differs from the assumed one, the isotope envelope scoring is meaningless -- we are comparing against the wrong theoretical distribution. Worse, the naive scorer had no mechanism to detect this: interstitial peaks were simply ignored.

### Solution: brute-force charge enumeration, best averagine fit wins

ms_deisotope determines charge by exhaustive fitting across all charge states in a configurable range (default z=1 to z=8):

1. For each experimental peak, call `_fit_all_charge_states()` to generate candidate isotope fits at every charge state.
2. For each charge, generate a theoretical averagine isotope pattern at the corresponding neutral mass (since neutral_mass = m/z * z - z * proton_mass, different charge assumptions imply different neutral masses and therefore different theoretical envelopes).
3. Score each fit using PenalizedMSDeconVFitter.
4. The scorer's `select` method picks the best charge state assignment.

For the AEVVFTANDSGPRR case: at z=4, the algorithm would find peaks at the expected 0.25 Da spacing, the averagine envelope for the z=4 neutral mass would match the observed intensities, and this fit would outscore the z=2 hypothesis (which has unexplained interstitial peaks).

### Implementation

```python
deconvoluted_peaks, _ = ms_deisotope.deconvolute_peaks(
    peaks,
    averagine=ms_deisotope.peptide,
    scorer=ms_deisotope.PenalizedMSDeconVFitter(10., 1.0),
    charge_range=(1, 8),         # try z=1 through z=8
    use_quick_charge=True,       # pre-screen by peak spacing
)

for peak in deconvoluted_peaks:
    print(f"neutral_mass={peak.neutral_mass:.4f}  charge={peak.charge}  "
          f"score={peak.score:.2f}")
```

The optional `use_quick_charge` optimization (from Hoopmann's QuickCharge algorithm) pre-screens charge states by examining peak spacing patterns: `1 / (m/z difference)` between neighboring peaks gives a candidate charge. Charge states missing a peak at M+0 or M+1 are immediately rejected, reducing computation.

Alternative approaches: THRASH uses the Patterson autocorrelation method (adapted from crystallography) to read the charge directly from spectral periodicity. UniDec/IsoDec uses a neural network with phase encoding -- multiplying m/z by candidate z and checking which z produces coherent phase alignment.

---

## Problem 3: Isotope Shape Mismatch (Co-eluting Interference)

### What the user saw

ALGISPFHEHAEVVFTANDSGP z=2+ had M+3 unreasonably high relative to the rest of the envelope. The overall shape did not match what theory predicts for a peptide of that mass. This indicates co-eluting interference: another ion's isotope peak lands on top of our M+3 position, inflating its intensity.

### Why it's a problem

Our naive scorer used a weighted sum that could still produce a passing score even when individual isotope peaks deviated wildly from theory. It had no quantitative measure of *shape* agreement between the observed and theoretical envelopes. A single contaminated peak (like the inflated M+3) could go undetected.

### Solution: G-test penalty, cosine similarity threshold >= 0.85

The PenalizedMSDeconVFitter in ms_deisotope combines two scoring dimensions:

1. **MSDeconV score** (per-peak): For each isotope peak, computes `sqrt(theoretical_intensity) * mass_accuracy * abundance_diff`. The abundance_diff term penalizes asymmetrically -- observed < theoretical gets a linear penalty, observed > theoretical gets a sqrt penalty (more tolerant of excess, which handles minor interference gracefully).

2. **G-test penalty** (whole-envelope shape): The scaled G-test measures divergence between the normalized observed and theoretical distributions:
   ```
   G = 2 * SUM(o_i * ln(o_i / e_i))
   ```
   where o_i and e_i are normalized to sum to 1. G = 0 means perfect shape match; large G means the envelope shape is wrong.

3. **Combined score**: `S = MSDeconV_score * (1 - |G_test| * penalty_factor)`. A badly shaped envelope (high G) heavily penalizes the final score regardless of how well individual peaks match in m/z.

For the ALGISPFHEHAEVVFTANDSGP case: the inflated M+3 would cause the G-test to detect a large divergence from the averagine prediction, driving the penalty term up and the final score down below the minimum_score threshold of 10.

Established thresholds from the literature:
- Hardklor (Hoopmann et al., 2007): cosine similarity >= 0.85-0.95 for acceptance
- MSFragger (Teo et al., 2021): Kullback-Leibler divergence with size-dependent thresholds (2-peak > 0.05, 3-peak > 0.1, 4-peak > 0.2), and ALL subsets of a cluster must also pass
- UniDec/IsoDec: cosine similarity >= 0.7

### Implementation

The default `deconvolute_peaks()` call already applies PenalizedMSDeconV scoring. Fits with a final score below 10 are rejected:

```python
scorer = ms_deisotope.PenalizedMSDeconVFitter(
    minimum_score=10.0,     # reject fits scoring below this
    penalty_factor=1.0,     # weight of G-test penalty
)
```

Additionally, the `InterferenceDetection` module measures how much of the total signal in the m/z region is NOT explained by the matched peaks: `interference = 1 - (matched_intensity / region_intensity)`. A score near 0 means clean signal; near 1 means the region is dominated by unmatched (interfering) peaks.

---

## Problem 4: Single Scan vs Integrated Signal

### Why single-scan is unreliable

Real peptides elute over many scans (typically 10-60 seconds of chromatographic time). A single MS1 scan may:

- Catch the peptide on the rising or falling edge of its elution peak, where intensities are low and noise dominates
- Miss it entirely due to stochastic noise fluctuation
- Sample it at a moment when a transient interfering ion is present
- Produce isotope ratios distorted by ion statistics at low abundance

All established feature detection tools (MaxQuant, OpenMS, Dinosaur) treat peptide features as three-dimensional objects in m/z, retention time, and intensity space. None of them score peptides from single scans.

### Solution: XIC extraction, chromatographic peak integration

The standard approach in the field:

1. **Extract XICs**: For each candidate monoisotopic m/z at charge z, collect intensity at m/z +/- tolerance across all scans. Do the same for M+1 (m/z + 1.00336/z), M+2, etc.
2. **Detect the elution peak**: Find the chromatographic peak in the M+0 XIC (Gaussian fit or local maximum detection).
3. **Define integration window**: Peak apex +/- 2*sigma (or a fixed number of scans).
4. **Integrate**: Sum (or fit-integrate) each isotope XIC within the window.
5. **Score the integrated envelope**: Compare the integrated isotope intensities against theory.

This converts noisy single-scan observations into robust, integrated measurements. MaxQuant fits Gaussian elution profiles and uses the area under the curve for quantification. Dinosaur builds "hills" (chromatographic traces) and uses cosine correlation against shifted averagine patterns. OpenMS FeatureFinderCentroided fits simultaneous Gaussians to all mass traces of a feature, enforcing consistent elution shape.

### Implementation

ms_deisotope provides two levels:

**Single-scan** (`deconvolute_peaks()`): Operates on one spectrum at a time. Does not integrate across scans, but its multi-iteration graph-based deconvolution resolves overlapping envelopes within each scan.

**Multi-scan** (`LCMSFeatureProcessor`): Works on pre-extracted LC-MS features. For each feature, it:
- Walks through time points where all isotopic features overlap (`FeatureSetIterator`)
- Scores the isotope fit at each time point using ScaledGTestFitter
- Aggregates scores across time using a thresholded mean (mean of scores above `threshold_scale * max(scores)`)
- Validates feature quality with three multiplicative dimensions: isotopic consistency (G-test), spacing regularity, and chromatographic peak shape (via `AdaptiveMultimodalChromatogramShapeFitter`)

For our pipeline, the practical approach is: run `deconvolute_peaks()` per scan for deisotoping, then build XICs from the deconvoluted peaks for cross-scan integration and final scoring.

---

## Chosen Solution: ms_deisotope

### Why ms_deisotope

- **Already installed** in our environment (`pip install ms-deisotope`)
- **Handles all 4 problems** in a single library: left-search for M-1 (Problem 1), exhaustive charge enumeration (Problem 2), PenalizedMSDeconV with G-test shape penalty (Problem 3), and LC-MS feature processing for multi-scan integration (Problem 4)
- **Open source** (Apache 2.0), pure Python with C extensions for performance
- **Actively maintained** by Joshua Klein (mobiusklein)
- **Well-architected**: peak dependency graph resolves overlapping envelopes globally, not greedily

### Key API: deconvolute_peaks()

```python
import ms_deisotope

deconvoluted_peaks, targeted = ms_deisotope.deconvolute_peaks(
    peaks,                                                    # centroided peak list
    averagine=ms_deisotope.peptide,                          # C4.9384 H7.7583 N1.3577 O1.4773 S0.0417
    scorer=ms_deisotope.PenalizedMSDeconVFitter(10., 1.0),   # MSDeconV * (1 - G-test penalty)
    charge_range=(1, 8),                                      # try all charge states
    left_search_limit=1,                                      # check for M-1
    use_quick_charge=True,                                    # pre-screen by spacing
)
```

### Scoring: PenalizedMSDeconVFitter

The default scorer combines:
- **MSDeconV** (per-peak): `sqrt(theo_intensity) * mass_accuracy * abundance_diff`
- **Scaled G-test** (whole-envelope): `G = 2 * SUM(o_i * ln(o_i / e_i))`
- **Combined**: `score = MSDeconV * (1 - |G| * penalty_factor)`
- **Threshold**: minimum_score = 10 (fits below this are rejected)

### Theoretical distributions: averagine + brainpy

1. Scale the averagine composition to the target neutral mass (divide mass by 111.1254 Da, multiply element counts, adjust hydrogen to match)
2. Compute the isotope distribution using **brainpy** (the BRAIN algorithm -- Baffling Recursive Algorithm for Isotopic distributioN calculations), a Cython-accelerated exact combinatorial method
3. Truncate trailing peaks to retain 95% of total signal
4. Cache patterns at 1 Da resolution, shift to exact m/z on retrieval

### Deconvolution pipeline

1. For each peak in the spectrum, for each charge state (8 down to 1):
   - Search left (1 position) for alternative monoisotopic peaks
   - Generate averagine theoretical pattern
   - Match to experimental peaks within 20 ppm
   - Scale theoretical to match experimental (sum method)
   - Score with PenalizedMSDeconV
2. Insert all fits into peak dependency graph (hypergraph: peaks are nodes, fits are hyperedges)
3. Find connected components, select best non-overlapping fits
4. Subtract assigned signal from peaks
5. Repeat (up to 10 iterations, checking 0.1% convergence)
6. Return `DeconvolutedPeakSet` with neutral masses, charge states, and scores

---

## References

### Tools

1. **ms_deisotope** -- Joshua Klein. Deisotoping and charge state deconvolution of complex mass spectra. Apache 2.0.
   - Source: https://github.com/mobiusklein/ms_deisotope
   - Docs: https://mobiusklein.github.io/ms_deisotope/docs/_build/html/

2. **brainpy** -- Joshua Klein. BRAIN isotopic distribution calculator (Cython).
   - Source: https://github.com/mobiusklein/brainpy
   - Install: `pip install brain-isotopic-distribution`

3. **MaxQuant** -- Cox & Mann. 3D feature detection with Gaussian elution profile fitting.
   - Not open source (free for academic use).

4. **OpenMS / pyopenms** -- FeatureFinderCentroided, seed-and-fit feature detection.
   - Docs: https://abibuilder.informatik.uni-tuebingen.de/archive/openms/Documentation/nightly/html/TOPP_FeatureFinderCentroided.html
   - Install: `pip install pyopenms`

5. **Dinosaur** -- Teleman et al. Hill-based feature detection with cosine correlation scoring.
   - Source: https://github.com/fickludd/dinosaur (JVM-based)

6. **MSFragger** -- Teo et al. Fast deisotoping with KL divergence scoring.

7. **Hardklor** -- Hoopmann et al. Multi-distribution deconvolution with dot product scoring.

8. **DeconTools / DeconEngineV2** -- PNNL. THRASH implementation with AreaFitter scoring.
   - Source: https://github.com/PNNL-Comp-Mass-Spec/DeconEngineV2

9. **UniDec / IsoDec** -- Michael Marty. Neural network charge detection with phase encoding.
   - Source: https://github.com/michaelmarty/UniDec

### Papers

1. Cox, J. & Mann, M. "MaxQuant enables high peptide identification rates, individualized p.p.b.-range mass accuracies and proteome-wide protein quantification." *Nature Biotechnology* 26, 1367-1372 (2008). https://www.nature.com/articles/nbt.1511

2. Teleman, J. et al. "Dinosaur: A Refined Open-Source Peptide MS Feature Detector." *J. Proteome Res.* 15, 2143-2151 (2016). https://pmc.ncbi.nlm.nih.gov/articles/PMC4933939/

3. Horn, D.M. et al. "Automated reduction and interpretation of high resolution electrospray mass spectra of large molecules." *J. Am. Soc. Mass Spectrom.* 11, 320-332 (2000).

4. Hoopmann, M.R. et al. "High-speed data reduction, feature detection, and MS/MS spectrum quality assessment of shotgun proteomics data sets using high-resolution mass spectrometry." *Anal. Chem.* 79, 5620-5632 (2007). https://pmc.ncbi.nlm.nih.gov/articles/PMC3891918/

5. Senko, M.W. et al. "Determination of monoisotopic masses and ion populations for large biomolecules from resolved isotopic distributions." *J. Am. Soc. Mass Spectrom.* 6, 229-233 (1995).

6. Teo, G.C. et al. "Fast Deisotoping Algorithm and Its Implementation in the MSFragger Search Engine." *J. Proteome Res.* 20, 498-505 (2021). https://pmc.ncbi.nlm.nih.gov/articles/PMC8864561/

7. Dittwald, P. et al. "BRAIN: A Universal Tool for High-Throughput Calculations of the Isotopic Distribution for Mass Spectrometry." *Anal. Chem.* 85, 1991-1998 (2013).

8. Loos, M. et al. "Prediction, Detection, and Validation of Isotope Clusters in Mass Spectrometry Data." *Anal. Chem.* 87, 5738-5744 (2015). https://pmc.ncbi.nlm.nih.gov/articles/PMC5192443/

9. Jaitly, N. et al. "Decon2LS: An open-source software package for automated processing and visualization of high resolution mass spectrometry data." *BMC Bioinformatics* 10, 87 (2009). https://pmc.ncbi.nlm.nih.gov/articles/PMC2666663/
