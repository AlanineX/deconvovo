# Isotope Scoring Research: Solving Four MS1 Peptide Scoring Problems

Research date: 2026-03-25

This document summarizes how established proteomics tools solve four specific problems
with naive single-scan, single-envelope MS1 peptide scoring.

---

## Problem 1: Single Scan vs Integrated Signal

**The problem:** Checking each MS1 scan independently discards information. Real peptides
elute over many scans (typically 10-60 seconds of chromatographic time), so a single scan
may catch the peptide at low abundance on the rising/falling edge, or miss it entirely due
to noise fluctuation.

### How Established Tools Solve This

#### MaxQuant: 3D Feature Detection (Cox & Mann, 2008)

MaxQuant treats peptide features as **three-dimensional objects** in m/z, retention time,
and intensity space. The algorithm:

1. **Detect centroids** in each MS1 scan
2. **Build mass traces** by connecting centroids with similar m/z across consecutive scans
3. **Cluster mass traces** into isotope patterns using graph-theoretical methods
4. **Fit Gaussian elution profiles** to each mass trace
5. **Integrate the area under the curve (AUC)** of the fitted elution profile

Key insight: MaxQuant requires sufficient intensity correlation over elution time between
isotope traces. The XIC (extracted ion chromatogram) area, not any single scan's intensity,
is used for quantification. For SILAC pairs, it allows small RT shifts between light/heavy
forms due to isotope effects.

Reference: Cox & Mann, "MaxQuant enables high peptide identification rates, individualized
p.p.b.-range mass accuracies and proteome-wide protein quantification," Nature Biotechnology
26, 1367-1372 (2008). https://www.nature.com/articles/nbt.1511

#### OpenMS FeatureFinderCentroided

OpenMS finds features through an iterative seed-and-fit approach:

1. **Seed identification**: Find local maxima in m/z/RT space. Seeds are scored by
   the geometric mean of intensity score, mass trace score, and isotope pattern score.
2. **Isotope trace assembly**: For each seed, extend mass traces in the RT dimension
   by following peaks with consistent m/z.
3. **Simultaneous Gaussian fitting**: Fit a Gaussian RT profile to ALL mass traces of a
   feature at the same time, enforcing consistent elution shape across isotope peaks.
4. **Truncation and scoring**: After fitting, traces below a quality threshold are removed.
   The reported feature intensity is based on the fitted model, not the raw (noisy) data.

The maximum slope parameter controls separation of overlapping elution peaks. Traces
are validated by checking that their intensity envelope across RT follows the expected
bell-shaped profile.

Reference: OpenMS FeatureFinderCentroided documentation.
https://abibuilder.informatik.uni-tuebingen.de/archive/openms/Documentation/nightly/html/TOPP_FeatureFinderCentroided.html

Python: `pip install pyopenms`

```python
import pyopenms as oms

# Load centroided MS1 data
exp = oms.MSExperiment()
opts = oms.PeakFileOptions()
opts.setMSLevels([1])
fh = oms.MzMLFile()
fh.setOptions(opts)
fh.load("input.mzML", exp)

# Run feature detection
ff = oms.FeatureFinder()
features = oms.FeatureMap()
seeds = oms.FeatureMap()
params = oms.FeatureFinder().getParameters("centroided")
ff.run("centroided", exp, features, params, seeds)

# Access integrated features
for f in features:
    print(f"m/z={f.getMZ():.4f}  RT={f.getRT():.1f}s  intensity={f.getIntensity():.0f}")
```

#### Dinosaur: Hill-Based Feature Detection (Teleman et al., 2016)

Dinosaur builds features through four steps:

1. **Centroiding**: Convert profile spectra to centroid peaks
2. **Hill construction**: Assemble centroids with similar m/z across consecutive scans
   into "hills" (= full chromatographic traces of one ion isotope). Uses a sliding
   average of the last 3 peaks in a hill for matching, reducing m/z fluctuation effects.
   Non-greedy matching: stores all candidates before selecting by minimal m/z difference.
3. **Hill clustering**: Group hills by theoretically possible m/z differences based on
   carbon/sulfur isotopes and expected charge states (typically 2-7).
4. **Isotopic deconvolution**: Extract charge-state-consistent features from clusters
   using intensity-based seeding and **cosine correlation** against shifted averagine
   peptide patterns.

Performance: 98% of MS/MS identifications had matching features in benchmark data
(competitors achieved 81-96%).

Reference: Teleman et al., "Dinosaur: A Refined Open-Source Peptide MS Feature Detector,"
J. Proteome Res. 15, 2143-2151 (2016). https://pmc.ncbi.nlm.nih.gov/articles/PMC4933939/

Source: https://github.com/fickludd/dinosaur (JVM-based, not pip-installable)

#### ms_deisotope

ms_deisotope operates on individual scans (not LC features), but its graph-based
deconvolution across multiple iterations effectively integrates signal by resolving
overlapping envelopes within each scan. For cross-scan integration, it would need to be
combined with a feature finder.

### Concrete Algorithm for Our Pipeline

**Recommended approach: Build XICs, then score the integrated envelope.**

```
1. For each candidate monoisotopic m/z at charge z:
   a. Extract XIC: collect intensity at m/z +/- tolerance across all scans
   b. Extract XICs for M+1, M+2, ... at m/z + i*1.00336/z
   c. Detect the elution peak in the M+0 XIC (Gaussian fit or local max)
   d. Define integration window: peak apex +/- 2*sigma (or +/- N scans)
   e. Sum (or fit-integrate) each isotope XIC within the window
   f. Score the INTEGRATED isotope envelope against theory
```

This converts noisy single-scan observations into robust, integrated measurements.

### Pip-Installable Libraries

| Library | Install | Use Case |
|---------|---------|----------|
| pyopenms | `pip install pyopenms` | Full feature detection (FeatureFinderCentroided) |
| ms_deisotope | `pip install ms-deisotope` | Per-scan deisotoping + charge deconvolution |
| pyteomics | `pip install pyteomics` | Mass calculations, mzML reading, XIC extraction |

---

## Problem 2: M-1 / M-2 Peaks as False Positive Indicator

**The problem:** If significant peaks exist BELOW the putative monoisotopic peak (at
m/z - 1.00336/z, m/z - 2*1.00336/z, etc.), the envelope likely belongs to a different,
heavier ion whose higher isotope peaks happen to overlap with our candidate.

### How Established Tools Solve This

#### ms_deisotope: left_search_limit Parameter

ms_deisotope explicitly addresses this via the `left_search_limit` parameter in its
`charge_state_determination()` method:

- **`left_search_limit`**: Maximum number of neutron shifts to search to the LEFT
  (decrease) from each query peak. If the algorithm finds peaks at M-1 that fit a
  larger isotope envelope better, it reassigns the monoisotopic peak.
- **`right_search_limit`**: Maximum number of neutron shifts to search to the RIGHT.

The graph-based deconvolver (`AveraginePeakDependenceGraphDeconvoluter`) builds a
dependency network where multiple isotopic fits may claim the same experimental peak.
It resolves conflicts by assigning each peak to its best-scoring solution via greedy
maximization over disjoint subgraphs. This means that if a peak at position X is better
explained as the M+2 of a heavier ion than as the M+0 of a lighter ion, the heavier
assignment wins.

Reference: ms_deisotope deconvolution documentation.
https://mobiusklein.github.io/ms_deisotope/docs/_build/html/deconvolution/deconvolution.html

#### MaxQuant / THRASH: Iterative Subtraction

The THRASH algorithm (Horn et al., 2000), used in Decon2LS and adapted in MaxQuant:

1. Start at the most intense peak in the spectrum
2. Determine charge state (Fourier/Patterson method)
3. Generate theoretical isotope distribution via Averagine + Mercury
4. Score the fit (least-squares, area fit, or chi-square)
5. **Delete the assigned isotopic peaks from the spectrum** (set intensities to zero)
6. Repeat from step 1 with remaining peaks

This subtraction approach inherently handles M-1 contamination: when a heavy ion is
processed first (it's more intense), its isotope peaks are removed from the spectrum.
When the algorithm later encounters the region where a lighter ion's M+0 would be,
the contaminating peaks are already gone.

Reference: Horn et al., "Automated reduction and interpretation of high resolution
electrospray mass spectra of large molecules," J. Am. Soc. Mass Spectrom. 11, 320-332
(2000).

#### Isotope Cluster Validation (Zubarev Group, 2016)

The approach from "Prediction, Detection, and Validation of Isotope Clusters in Mass
Spectrometry Data" uses mass-specific abundance ratio validation:

1. For each mass window (e.g., 10-250 Da), compute database-derived confidence
   intervals for the ratio M+0/M+1, M+0/M+2, etc.
2. If the observed ratio of the putative monoisotopic peak to its successor falls
   outside the 99% confidence interval, the cluster is SPLIT -- the leading peak is
   removed as likely belonging to a different ion (e.g., hydrogen loss artifact with
   mass difference ~1.008, which mimics a 13C isotope spacing of ~1.003).

Reference: Loos et al., "Prediction, Detection, and Validation of Isotope Clusters in
Mass Spectrometry Data," Anal. Chem. 87, 5738-5744 (2015).
https://pmc.ncbi.nlm.nih.gov/articles/PMC5192443/

### Concrete Algorithm for Our Pipeline

```
def check_below_monoisotopic(spectrum, mono_mz, charge, tolerance_da=0.02):
    """
    Check for peaks below the monoisotopic m/z that indicate interference.
    Returns a penalty factor (0.0 = severe contamination, 1.0 = clean).
    """
    neutron_spacing = 1.00335483 / charge
    penalty = 1.0

    for i in range(1, 3):  # Check M-1 and M-2
        check_mz = mono_mz - i * neutron_spacing
        # Find nearest peak within tolerance
        peak_intensity = find_peak(spectrum, check_mz, tolerance_da)
        mono_intensity = find_peak(spectrum, mono_mz, tolerance_da)

        if peak_intensity > 0 and mono_intensity > 0:
            ratio = peak_intensity / mono_intensity
            # For a clean monoisotopic peak, M-1 should be ~0
            # Any significant signal is penalizing
            if i == 1 and ratio > 0.1:
                penalty *= max(0.0, 1.0 - ratio)
            elif i == 2 and ratio > 0.05:
                penalty *= max(0.0, 1.0 - 2 * ratio)

    return penalty
```

### Pip-Installable Libraries

| Library | Install | Relevant Feature |
|---------|---------|-----------------|
| ms_deisotope | `pip install ms-deisotope` | Graph-based deconvolution with left_search_limit |
| pyopenms | `pip install pyopenms` | Deisotoping via Deconv2D, FeatureFinder |

---

## Problem 3: Peaks Between Isotopes = Charge State Ambiguity

**The problem:** If peaks appear at half-Da spacing between M+0 and M+1, the signal
likely comes from a higher charge state ion. For example, a z=2 ion has isotope peaks
spaced at 1.003/2 = 0.5017 Da. If we are trying to score a z=1 envelope at the same
position, the interstitial peaks at +0.5 Da indicate the true charge is 2+ (or higher).

### How Established Tools Solve This

#### THRASH: Fourier Transform / Patterson Charge Determination (Horn et al., 2000)

The original THRASH algorithm determines charge state using the **Patterson autocorrelation
method** (adapted from crystallography):

1. Take a local window of the spectrum around the peaks of interest
2. Compute the Fourier transform of the peak pattern
3. Examine the autocorrelation (Patterson function) -- peaks in the autocorrelation
   appear at spacings corresponding to 1/z Da
4. The dominant spacing reveals the charge state: spacing of 1.0 Da = z=1,
   spacing of 0.5 Da = z=2, spacing of 0.33 Da = z=3, etc.

This is elegant because it does not require prior knowledge of the monoisotopic peak --
it directly reads the charge from the spectral periodicity.

Reference: Horn et al., "Automated reduction and interpretation of high resolution
electrospray mass spectra of large molecules," J. Am. Soc. Mass Spectrom. 11, 320-332
(2000).

#### ms_deisotope: Exhaustive Charge State Fitting

ms_deisotope determines charge by brute-force fitting:

1. For each experimental peak, call `_fit_all_charge_states()` to generate candidate
   isotope fits across the specified charge range (default z=1 to z=8).
2. For each charge state, generate a theoretical averagine isotope pattern at that
   charge spacing.
3. Score each fit using the configured scorer (MSDeconVFitter, PenalizedMSDeconVFitter,
   or ScaledGTestFitter).
4. Use the scorer's `select` method to pick the best charge state assignment.

The **`use_quick_charge`** optimization pre-screens charge states: any charge state
missing a peak at the monoisotopic or A+1 position is immediately rejected, reducing
computation.

For overlapping patterns at different charge states, the graph-based deconvolver
(`AveraginePeakDependenceGraphDeconvoluter`) resolves conflicts:
- Peaks are nodes; isotopic fits are hyperedges
- Disjoint subgraphs are solved independently
- Each peak is assigned to exactly one isotopic envelope

```python
import ms_deisotope

deconvoluted_peaks, _ = ms_deisotope.deconvolute_peaks(
    peaks,
    averagine=ms_deisotope.peptide,
    scorer=ms_deisotope.PenalizedMSDeconVFitter(10., 1.0),
    charge_range=(1, 8),              # Try z=1 through z=8
    left_search_limit=1,              # Check 1 peak below monoisotopic
    right_search_limit=0,
    use_quick_charge=True,
)

for peak in deconvoluted_peaks:
    print(f"neutral_mass={peak.neutral_mass:.4f}  charge={peak.charge}  "
          f"score={peak.score:.2f}")
```

Reference: ms_deisotope documentation.
https://mobiusklein.github.io/ms_deisotope/docs/_build/html/deconvolution/deconvolution.html

#### Dinosaur: Charge-Consistent Feature Deconvolution

Dinosaur clusters hills (chromatographic traces) by checking all theoretically possible
m/z differences based on carbon and sulfur natural isotopes across charge states 2-7.
It then deconvolves clusters into charge-state-consistent features:

- Seeds isotope patterns from the 100 most intense hills in each cluster
- Matches intensity profiles against shifted averagine patterns using **cosine
  correlation scoring**
- Selects the highest-scoring charge state + shift combination

Reference: Teleman et al. (2016). https://pmc.ncbi.nlm.nih.gov/articles/PMC4933939/

#### MSFragger: Fast Deisotoping (Teo et al., 2021)

MSFragger's algorithm determines charge states by:

1. Generate all possible peak clusters based solely on m/z distances
2. Consider charge states from 1+ to (precursor charge - 1)
3. Validate each cluster against averagine using **Kullback-Leibler divergence**
4. Apply size-dependent score thresholds: 2 peaks > 0.05; 3 peaks > 0.1; 4 peaks > 0.2
5. Critically: **ALL subsets of a cluster** must also score above threshold

Reference: Teo et al., "Fast Deisotoping Algorithm and Its Implementation in the
MSFragger Search Engine," J. Proteome Res. 20, 498-505 (2021).
https://pmc.ncbi.nlm.nih.gov/articles/PMC8864561/

#### Hardklor: Multi-Distribution Deconvolution (Hoopmann et al., 2007)

Hardklor handles charge state ambiguity through iterative depth:

1. Divide spectrum into ~5-6 Th windows
2. Generate averagine models for candidate masses at each charge state
3. First try fitting a single isotope distribution per window
4. If correlation < threshold (0.85-0.95), try combinations of 2 distributions
5. Continue increasing deconvolution depth (typically up to 3, max 5)
6. Accept the first combination that explains the observed peaks above threshold

The scoring uses the **dot product** (cosine similarity) between observed and
theoretical peak vectors.

Reference: Hoopmann et al., "High-speed data reduction, feature detection, and MS/MS
spectrum quality assessment of shotgun proteomics data sets using high-resolution mass
spectrometry," Anal. Chem. 79, 5620-5632 (2007).
https://pmc.ncbi.nlm.nih.gov/articles/PMC3891918/

### Concrete Algorithm for Our Pipeline

```
def check_interstitial_peaks(spectrum, mono_mz, charge, tolerance_da=0.02):
    """
    Check for peaks between expected isotope positions that suggest
    the true charge state is higher than assumed.
    Returns True if interstitial peaks suggest wrong charge assignment.
    """
    neutron_spacing = 1.00335483 / charge

    # Check at half-spacing positions (would indicate charge = 2*charge)
    mono_intensity = find_peak(spectrum, mono_mz, tolerance_da)
    if mono_intensity == 0:
        return False

    interstitial_count = 0
    interstitial_intensity = 0

    for i in range(3):  # Check between M+i and M+i+1
        half_mz = mono_mz + (i + 0.5) * neutron_spacing
        peak_int = find_peak(spectrum, half_mz, tolerance_da)
        if peak_int > 0.1 * mono_intensity:  # >10% of mono = suspicious
            interstitial_count += 1
            interstitial_intensity += peak_int

    # If 2+ interstitial peaks found, likely wrong charge
    return interstitial_count >= 2
```

**Better approach: Try all charge states and pick the best fit**, as ms_deisotope does.

### Pip-Installable Libraries

| Library | Install | Charge Determination Method |
|---------|---------|---------------------------|
| ms_deisotope | `pip install ms-deisotope` | Exhaustive fitting + graph-based resolution |
| pyopenms | `pip install pyopenms` | FeatureFinder with charge range, Deisotoper |

---

## Problem 4: Isotope Shape Does Not Match Theory

**The problem:** If the observed isotope envelope (e.g., M+3 anomalously high) does not
match the theoretical distribution for a peptide of that mass, the signal likely contains
interference from a co-eluting ion. We need a quantitative score for how well the
observed pattern matches theory.

### Theoretical Distribution Generation

#### Averagine Model (Senko et al., 1995)

The "average amino acid" has composition:

```
C_4.9384  H_7.7583  N_1.3577  O_1.4773  S_0.0417
Average mass: 111.1254 Da
```

To estimate an elemental formula for a peptide of mass M:
1. Number of averagine units = M / 111.1254
2. Multiply by each element count, round C/N/O/S to nearest integer
3. Adjust H count to match the target mass exactly

Reference: Senko et al., "Determination of monoisotopic masses and ion populations for
large biomolecules from resolved isotopic distributions," J. Am. Soc. Mass Spectrom. 6,
229-233 (1995).

#### Mercury Algorithm (Rockwood, 1995)

The Mercury algorithm generates theoretical isotope distributions using the **Fast Fourier
Transform (FFT)** method:

1. For each element, represent the isotope distribution as a polynomial
2. Multiply the polynomials (convolution in mass domain = multiplication in Fourier domain)
3. Use FFT for efficient computation
4. Result: theoretical peak positions and relative intensities

Properties: Places isotope peaks within millidaltons of true centroids. The height of
each nominal peak equals the integrated area of the corresponding exact peak.

Reference: Rockwood, "Ultrahigh-speed calculation of isotope distributions," Rapid
Commun. Mass Spectrom. 9, 103 (1995).
https://pubs.acs.org/doi/10.1021/ac951158i

#### BRAIN Algorithm (Claesen et al., 2012)

The Baffling Recursive Algorithm for Isotopic distributioN calculations uses Newton's
identities to recursively compute center-masses and probabilities of isotopic peaks
without polynomial expansion. Implemented in the `brainpy` Python library.

Reference: Dittwald et al., "BRAIN: A Universal Tool for High-Throughput Calculations
of the Isotopic Distribution for Mass Spectrometry," Anal. Chem. 85, 1991-1998 (2013).

### Scoring Methods

#### 1. Cosine Similarity (Dot Product)

The most widely used method. Used by Dinosaur, Hardklor, and many others.

**Formula:**
```
                    SUM_i(obs_i * theo_i)
cos(theta) = -----------------------------------
             sqrt(SUM_i(obs_i^2)) * sqrt(SUM_i(theo_i^2))
```

Where `obs_i` and `theo_i` are the intensities of the i-th isotope peak in the
observed and theoretical distributions respectively.

- Range: 0.0 (no match) to 1.0 (perfect match)
- Insensitive to overall intensity scaling
- Hardklor uses threshold of 0.85-0.95 for acceptance

**Implementation:**
```python
import numpy as np

def cosine_similarity(observed, theoretical):
    """Both inputs are arrays of isotope peak intensities."""
    obs = np.array(observed, dtype=float)
    theo = np.array(theoretical, dtype=float)
    dot = np.dot(obs, theo)
    norm = np.linalg.norm(obs) * np.linalg.norm(theo)
    return dot / norm if norm > 0 else 0.0
```

#### 2. G-Test (Scaled) -- ms_deisotope ScaledGTestFitter

A likelihood-ratio statistic measuring divergence between observed and expected
distributions:

**Formula:**
```
G = 2 * SUM_i( o_i * (log(o_i) - log(e_i)) )
```

Where `o_i` and `e_i` are normalized (sum-to-1) observed and expected intensities.

- Minimizing fitter: lower = better fit
- Intensity-independent (both distributions normalized)
- Sensitive but less tolerant of detector noise

#### 3. Kullback-Leibler Divergence -- MSFragger

MSFragger validates isotope clusters using KL divergence against averagine:

**Formula:**
```
KL(P || Q) = SUM_i( P_i * log(P_i / Q_i) )
```

With size-dependent thresholds: 2-peak clusters > 0.05, 3-peak > 0.1, 4-peak > 0.2.
ALL subsets of a cluster must also pass.

#### 4. MSDeconV Score -- ms_deisotope MSDeconVFitter

A maximizing scorer combining m/z accuracy and intensity match:

**Formula:**
```
S(e, t) = sqrt(intensity(t)) * s_mz(e, t) * s_int(e, t)
```

Where:
- `s_mz`: Penalizes m/z deviation (0 if > tolerance)
- `s_int`: Asymmetric intensity mismatch penalty (different formulas for under-
  and over-prediction)
- Higher = better fit
- Default minimum_score threshold: 10.0

#### 5. PenalizedMSDeconV Score -- ms_deisotope PenalizedMSDeconVFitter

Combines MSDeconV with a G-test shape penalty:

**Formula:**
```
S(e, t) = M(e, t) * (1 - penalty_factor * G(e, t))
```

Where M is the MSDeconV score and G is the normalized G-test. This is preferred for
high-quality isotopic patterns because it weights both intensity magnitude AND
distribution shape.

#### 6. Least-Squares Fit -- Decon2LS / THRASH

Uses the R^2 coefficient of determination:

**Formula:**
```
R^2 = 1 - SUM_i(obs_i_hat - theo_i_hat)^2 / SUM_i(theo_i_hat^2)
```

Where `_hat` denotes normalization by respective maxima.

Decon2LS offers three selectable scoring functions: area fit, peak fit, and chi-square
fit. The user-specified intensity threshold prevents noise from affecting the score.

#### 7. Chi-Square Test

Classical goodness-of-fit:

**Formula:**
```
chi^2 = SUM_i( (obs_i - expected_i)^2 / expected_i )
```

Used as an option in Decon2LS/THRASH. Sensitive to total counts and requires appropriate
normalization for mass spectrometry data.

### Tool Implementations Summary

| Tool | Scoring Method | Threshold | Reference |
|------|---------------|-----------|-----------|
| Hardklor | Dot product (cosine) | 0.85-0.95 | Hoopmann et al. 2007 |
| Dinosaur | Cosine correlation | Not published | Teleman et al. 2016 |
| ms_deisotope | MSDeconV, PenalizedMSDeconV, G-test | 10.0 (MSDeconV) | mobiusklein/ms_deisotope |
| MSFragger | Kullback-Leibler divergence | Size-dependent | Teo et al. 2021 |
| Decon2LS/THRASH | Least-squares / chi-square / area fit | User-specified | Jaitly et al. 2009 |
| MaxQuant | Correlation-based (proprietary) | Not published | Cox & Mann 2008 |
| OpenMS | Isotope pattern score (geometric mean) | Configurable | Kohlbacher et al. |

### Concrete Algorithm for Our Pipeline

```python
import numpy as np
from brainpy import isotopic_variants  # pip install brain-isotopic-distribution

def score_isotope_envelope(observed_intensities, neutral_mass, charge,
                           n_peaks=5, method='cosine'):
    """
    Score observed isotope envelope against Averagine theoretical prediction.

    Parameters
    ----------
    observed_intensities : list of float
        Observed intensities at M+0, M+1, M+2, ...
    neutral_mass : float
        Neutral monoisotopic mass of the peptide
    charge : int
        Charge state
    n_peaks : int
        Number of isotope peaks to consider
    method : str
        'cosine', 'gtest', or 'kl'

    Returns
    -------
    float : score (higher = better for cosine; lower = better for gtest/kl)
    """
    # Generate theoretical distribution from averagine
    averagine = {"C": 4.9384, "H": 7.7583, "N": 1.3577, "O": 1.4773, "S": 0.0417}
    avg_mass = 111.1254
    n_residues = neutral_mass / avg_mass

    composition = {
        "C": round(averagine["C"] * n_residues),
        "N": round(averagine["N"] * n_residues),
        "O": round(averagine["O"] * n_residues),
        "S": round(averagine["S"] * n_residues),
    }
    # Adjust H to match mass
    current_mass = (composition["C"] * 12 + composition["N"] * 14.003074
                    + composition["O"] * 15.994915 + composition["S"] * 31.972071)
    composition["H"] = round((neutral_mass - current_mass) / 1.00794)

    # Compute theoretical isotope distribution
    theo_dist = isotopic_variants(composition, npeaks=n_peaks)
    theo = np.array([p.intensity for p in theo_dist[:n_peaks]])

    # Truncate/pad observed to same length
    obs = np.zeros(n_peaks)
    for i, val in enumerate(observed_intensities[:n_peaks]):
        obs[i] = max(0, val)

    # Normalize
    theo_norm = theo / theo.sum() if theo.sum() > 0 else theo
    obs_norm = obs / obs.sum() if obs.sum() > 0 else obs

    if method == 'cosine':
        dot = np.dot(obs, theo)
        norm = np.linalg.norm(obs) * np.linalg.norm(theo)
        return dot / norm if norm > 0 else 0.0

    elif method == 'gtest':
        # G-test (lower = better)
        mask = (obs_norm > 0) & (theo_norm > 0)
        if mask.sum() == 0:
            return float('inf')
        return 2.0 * np.sum(obs_norm[mask] * (np.log(obs_norm[mask])
                                                - np.log(theo_norm[mask])))

    elif method == 'kl':
        # KL divergence (lower = better)
        mask = (obs_norm > 0) & (theo_norm > 0)
        if mask.sum() == 0:
            return float('inf')
        return np.sum(obs_norm[mask] * np.log(obs_norm[mask] / theo_norm[mask]))
```

### Pip-Installable Libraries for Theoretical Distributions

| Library | Install | Method | Notes |
|---------|---------|--------|-------|
| brain-isotopic-distribution (brainpy) | `pip install brain-isotopic-distribution` | BRAIN algorithm | Used by ms_deisotope internally. Fast Cython implementation. |
| IsoSpecPy | `pip install IsoSpecPy` | Fine structure calculator | Hyperfast, supports arbitrary coverage thresholds |
| molmass | `pip install molmass` | Polynomial expansion | Also calculates average/monoisotopic mass, elemental composition |
| pyteomics | `pip install pyteomics` | mass.isotopic_composition_abundance | General proteomics toolkit with mass calculation |
| pyopenms | `pip install pyopenms` | CoarseIsotopePatternGenerator | Part of the full OpenMS framework |
| ms_deisotope | `pip install ms-deisotope` | Averagine + brainpy | Full deconvolution pipeline, multiple scorers |

---

## Combined Scoring Strategy: Putting It All Together

Based on this research, a robust MS1 peptide scorer should:

### Step 1: Build XICs and Integrate (solves Problem 1)
- Extract ion chromatograms for each isotope peak across all scans
- Fit Gaussian (or similar) elution profiles
- Integrate area under the curve for each isotope peak
- Use integrated intensities for all subsequent scoring

### Step 2: Check Below Monoisotopic (solves Problem 2)
- Extract XIC at M-1 and M-2 positions
- If M-1 intensity > 10-15% of M+0, apply penalty or reject
- Better: use graph-based deconvolution (ms_deisotope) to determine if the
  envelope is better explained as part of a heavier ion

### Step 3: Determine Charge State (solves Problem 3)
- Try all plausible charge states (z=1 through z=6 for typical peptides)
- For each z, check for interstitial peaks at half-spacing
- Score each charge hypothesis against averagine theory
- Select the charge state with the best isotope fit score

### Step 4: Score Isotope Shape (solves Problem 4)
- Generate theoretical distribution via averagine + BRAIN/brainpy
- Compare observed (integrated) envelope to theoretical using cosine similarity
- Apply PenalizedMSDeconV-style scoring: weight by both intensity and shape match
- Reject envelopes with cosine < 0.85 or G-test > threshold

### Composite Score Formula

```
final_score = isotope_shape_score * (1 - M_minus1_penalty) * charge_confidence * XIC_quality

where:
  isotope_shape_score = cosine_similarity(observed_envelope, theoretical_envelope)
  M_minus1_penalty    = clamp(M_minus1_intensity / M0_intensity, 0, 1)
  charge_confidence   = best_charge_score / second_best_charge_score  (or 1.0 if only one)
  XIC_quality         = R^2 of Gaussian fit to elution profile (0 to 1)
```

---

## References

1. Cox, J. & Mann, M. "MaxQuant enables high peptide identification rates, individualized
   p.p.b.-range mass accuracies and proteome-wide protein quantification." Nature
   Biotechnology 26, 1367-1372 (2008).
   https://www.nature.com/articles/nbt.1511

2. Teleman, J. et al. "Dinosaur: A Refined Open-Source Peptide MS Feature Detector."
   J. Proteome Res. 15, 2143-2151 (2016).
   https://pmc.ncbi.nlm.nih.gov/articles/PMC4933939/

3. Horn, D.M. et al. "Automated reduction and interpretation of high resolution
   electrospray mass spectra of large molecules." J. Am. Soc. Mass Spectrom. 11,
   320-332 (2000).

4. Hoopmann, M.R. et al. "High-speed data reduction, feature detection, and MS/MS
   spectrum quality assessment of shotgun proteomics data sets using high-resolution
   mass spectrometry." Anal. Chem. 79, 5620-5632 (2007).
   https://pmc.ncbi.nlm.nih.gov/articles/PMC3891918/

5. Senko, M.W. et al. "Determination of monoisotopic masses and ion populations for
   large biomolecules from resolved isotopic distributions." J. Am. Soc. Mass Spectrom.
   6, 229-233 (1995).

6. Rockwood, A.L. "Ultrahigh-speed calculation of isotope distributions." Rapid Commun.
   Mass Spectrom. 9, 103 (1995).
   https://pubs.acs.org/doi/10.1021/ac951158i

7. Teo, G.C. et al. "Fast Deisotoping Algorithm and Its Implementation in the MSFragger
   Search Engine." J. Proteome Res. 20, 498-505 (2021).
   https://pmc.ncbi.nlm.nih.gov/articles/PMC8864561/

8. Jaitly, N. et al. "Decon2LS: An open-source software package for automated processing
   and visualization of high resolution mass spectrometry data." BMC Bioinformatics 10,
   87 (2009).
   https://pmc.ncbi.nlm.nih.gov/articles/PMC2666663/

9. Loos, M. et al. "Prediction, Detection, and Validation of Isotope Clusters in Mass
   Spectrometry Data." Anal. Chem. 87, 5738-5744 (2015).
   https://pmc.ncbi.nlm.nih.gov/articles/PMC5192443/

10. Dittwald, P. et al. "BRAIN: A Universal Tool for High-Throughput Calculations of the
    Isotopic Distribution for Mass Spectrometry." Anal. Chem. 85, 1991-1998 (2013).

11. OpenMS FeatureFinderCentroided documentation.
    https://abibuilder.informatik.uni-tuebingen.de/archive/openms/Documentation/nightly/html/TOPP_FeatureFinderCentroided.html

12. ms_deisotope documentation and source code.
    https://github.com/mobiusklein/ms_deisotope
    https://mobiusklein.github.io/ms_deisotope/docs/_build/html/

13. brainpy (BRAIN isotopic distribution calculator).
    https://github.com/mobiusklein/brainpy

14. IsoSpecPy (fine isotopic structure calculator).
    https://pypi.org/project/IsoSpecPy/
    https://matteolacki.github.io/IsoSpec/

15. molmass (molecular mass and isotope distributions).
    https://pypi.org/project/molmass/

16. Pyteomics (Python proteomics framework).
    https://pyteomics.readthedocs.io/
