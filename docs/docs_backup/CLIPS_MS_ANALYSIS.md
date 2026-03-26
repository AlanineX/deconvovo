# CLIPS-MS and Alternative Tool Source Code Analysis

Research date: 2026-03-25

This document analyzes CLIPS-MS and alternative open-source tools for isotope pattern
scoring, charge state determination, false positive rejection, and feature detection in
mass spectrometry data. Since CLIPS-MS is proprietary (Cerno Bioscience), the analysis
focuses on its published algorithm description plus deep source code analysis of three
open-source alternatives: ms_deisotope, DeconTools/DeconEngineV2, and UniDec/IsoDec.

---

## 1. CLIPS-MS (Cerno Bioscience) -- PROPRIETARY, NO SOURCE CODE

**Product:** MassWorks / CLIPS (Calibrated Line-shape Isotope Profile Search)
**Website:** https://cernobioscience.com/clips/
**Source code:** NOT available (closed-source commercial product)

### Algorithm Overview

CLIPS is an isotope-pattern-based algorithm for chemical formula determination. Unlike
tools that rely solely on mass accuracy, CLIPS adds a second dimension: **Spectral
Accuracy** -- matching the full isotope envelope shape.

### Key Technical Details (from published descriptions)

**Step 1: TrueCal Instrument Line Shape Calibration**
- Calibrates the actual instrument line shape to a known mathematical function
- This allows generating theoretical isotope profiles using the same line shape
  as the calibrated instrument
- Works even on unit-mass-resolution instruments (single quadrupole)

**Step 2: Formula Candidate Generation**
- Conventional formula determination based on mass accuracy
- Applies search constraints: mass tolerance, allowed elements, DBE, electron state, charge

**Step 3: CLIPS Isotope Pattern Scoring**
- For each candidate formula, calculates the "true mass spectrum" using the calibrated
  line shape function
- Matches the calculated spectrum against the unknown ion over a defined mass range
- Calculates **Spectral Accuracy** as:

```
Spectral Accuracy = 100% - (RMSE * 100%)
```

Where RMSE is the root mean square error between normalized theoretical and observed
isotope profiles. An RMSE of 0.005 yields Spectral Accuracy of 99.5%.

**Step 4: Ranking**
- Candidate formulas ranked by Spectral Accuracy
- Claims <1% relative spectral error for reliable identifications
- Can discriminate candidates with spectral differences as small as 0.1%

**Variant: sCLIPS (self-Calibrating)**
- For high-resolution instruments (TOF, Orbitrap, FT-ICR)
- Does not require running calibration standards
- Uses self-calibration of the line shape from the data itself

### Limitations for Our Use Case
CLIPS is designed for **formula determination** (identifying unknowns), not for
**deconvolution of known peptide isotope patterns** in HDX-MS. It requires commercial
software (MassWorks) and is not applicable to our pipeline.

---

## 2. ms_deisotope (mobiusklein/ms_deisotope) -- PRIMARY ANALYSIS

**Repository:** https://github.com/mobiusklein/ms_deisotope
**Language:** Python (with C extensions)
**License:** Apache 2.0
**Author:** Joshua Klein (mobiusklein)
**Description:** A library for deisotoping and charge state deconvolution of complex mass spectra

This is the most relevant and thoroughly documented open-source tool for our purposes.

### 2.1 Isotope Pattern Scoring

**File:** `src/ms_deisotope/scoring.py`

The library provides **7 scoring functions**, all inheriting from `IsotopicFitterBase`.
Scorers are either minimizing (lower = better) or maximizing (higher = better).

#### MSDeconVFitter (lines 415-464) -- MAXIMIZING

Implementation of the scoring from Liu et al. (2010) "Deconvolution and database search
of complex tandem mass spectra of intact proteins."

For each peak pair (observed, theoretical), computes three factors:

```python
def score_peak(self, obs, theo, mass_error_tolerance=0.02,
               minimum_signal_to_noise=1):
    if obs.signal_to_noise < minimum_signal_to_noise:
        return 0.

    mass_error = np.abs(obs.mz - theo.mz)
    if mass_error <= mass_error_tolerance:
        mass_accuracy = 1 - mass_error / mass_error_tolerance
    else:
        mass_accuracy = 0

    if (obs.intensity < theo.intensity and
        (((theo.intensity - obs.intensity) / obs.intensity) <= 1)):
        abundance_diff = 1 - ((theo.intensity - obs.intensity) / obs.intensity)
    elif (obs.intensity >= theo.intensity and
          (((obs.intensity - theo.intensity) / obs.intensity) <= 1)):
        abundance_diff = np.sqrt(
            1 - ((obs.intensity - theo.intensity) / obs.intensity))
    else:
        abundance_diff = 0.

    score = np.sqrt(theo.intensity) * mass_accuracy * abundance_diff
    return score
```

The total score sums across all peaks: `score = sum(score_peak(obs, theo) for each pair)`

Key properties:
- **mass_accuracy**: Linear penalty from 0 to 1 based on m/z error vs tolerance
- **abundance_diff**: Asymmetric -- uses sqrt when obs > theo (tolerates excess intensity
  more than deficit), linear when obs < theo
- **weighting**: Each peak weighted by sqrt(theoretical intensity), so dominant isotope
  peaks contribute more

#### PenalizedMSDeconVFitter (lines 467-482) -- DEFAULT SCORER

Combines MSDeconV with a G-test penalty:

```python
def evaluate(self, peaklist, observed, expected, **kwargs):
    score = self.msdeconv.evaluate(peaklist, observed, expected)
    penalty = abs(self.penalizer.evaluate(peaklist, observed, expected))
    return score * (1 - penalty * self.penalty_factor)
```

Formula: `S(e, t) = MSDeconV(e, t) * (1 - |G_scaled(e, t)|)`

This penalizes fits where the overall intensity distribution diverges from theoretical,
even if individual peaks score well.

#### ScaledGTestFitter (lines 349-372) -- MINIMIZING

G-test (log-likelihood ratio) after normalizing both distributions to sum to 1:

```python
def evaluate(self, peaklist, observed, expected, **kwargs):
    total_observed = sum(p.intensity for p in observed)
    total_expected = sum(p.intensity for p in expected)
    total_expected += eps  # eps = 1e-4
    normalized_observed = [obs.intensity / total_observed for obs in observed]
    normalized_expected = [theo.intensity / total_expected for theo in expected]
    g_score = 2 * sum([obs * np.log(obs / theo)
        for obs, theo in zip(normalized_observed, normalized_expected)])
    return g_score
```

Formula: `G = 2 * sum(o_i * (log(o_i) - log(e_i)))` where both are normalized to sum=1.

#### LeastSquaresFitter (lines 385-409) -- MINIMIZING

Normalized least-squares coefficient of determination:

```python
def evaluate(self, peaklist, observed, expected, **kwargs):
    exp_max = max(p.intensity for p in observed)
    theo_max = max(p.intensity for p in expected)
    sum_of_squared_errors = 0
    sum_of_squared_theoreticals = 0
    for e, t in zip(observed, expected):
        normed_expr = e.intensity / exp_max
        normed_theo = t.intensity / theo_max
        sum_of_squared_errors += (normed_theo - normed_expr) ** 2
        sum_of_squared_theoreticals += normed_theo ** 2
    return sum_of_squared_errors / sum_of_squared_theoreticals
```

Both distributions are max-normalized before comparison. Score is ratio of squared-error
to squared-theoretical.

#### Other Scorers

- **ChiSquareFitter** (lines 375-382): `sum((obs - theo)^2 / theo)` -- standard chi-squared
- **GTestFitter** (lines 330-346): Raw G-test without normalization
- **DotProductFitter** (lines 527-529): `sum(e_i * t_i)` -- simple dot product
- **DistinctPatternFitter** (lines 507-515): Scaled G-test weighted by interference detection

#### InterferenceDetection (lines 490-504)

Measures how much of the m/z region's total signal is NOT explained by the matched peaks:

```python
def detect_interference(self, experimental_peaks, lower=None, upper=None):
    region = self.peaklist.between(lower, upper)
    included_intensity = sum(p.intensity for p in experimental_peaks)
    region_intensity = sum(p.intensity for p in region)
    score = 1 - (included_intensity / region_intensity)
    return score
```

Returns 0 when all signal is explained; approaches 1 when interfering peaks dominate.

### 2.2 Charge State Determination

**File:** `src/ms_deisotope/deconvolution/utils.py` (ChargeIterator)
**File:** `src/ms_deisotope/deconvolution/exhaustive.py` (charge iteration logic)

#### Approach: Exhaustive Charge State Search

The system does NOT determine charge state from peak spacing a priori. Instead, it
**tries all charge states** in a configured range and picks the best-scoring fit.

```python
# ChargeIterator iterates from HIGH charge to LOW (descending absolute value)
class ChargeIterator(object):
    def __init__(self, lo, hi):
        self.set_bounds(lo, hi)
        self.make_sequence()

    def make_sequence(self):
        self.index = 0
        self.size = self.upper - self.lower + 1
        # Iterates from upper to lower charge state
        self.values = [self.sign * (self.upper - i) for i in range(self.size)]
```

Default charge range is (1, 8), iterable from 8 down to 1.

#### QuickCharge Optimization

An optional `quick_charge()` function pre-filters likely charge states by examining
peak spacing patterns before exhaustive fitting. When enabled, ChargeIterator calls
`sequence_from_quickcharge()` to prioritize likely candidates.

#### How Charge States Are Evaluated

For each (peak, charge) pair from `_get_all_peak_charge_pairs()`:
1. Generate theoretical isotopic pattern using averagine at that charge
2. Match theoretical peaks to experimental peaks within PPM tolerance
3. Scale theoretical pattern to experimental intensity
4. Score the fit using the configured scorer
5. Reject if only 1 real peak matched (for charge > 1)
6. Reject if scorer.reject(fit) returns True

The best fit across all charge states is selected by the peak dependency graph solver.

### 2.3 False Positive Rejection

Multiple layers of false positive rejection:

#### Layer 1: Placeholder Detection (base.py `_check_fit`)

```python
def _check_fit(self, fit):
    if len(drop_placeholders(fit.experimental)) == 1 and fit.charge > 1:
        return False  # Reject: only 1 real peak for charge > 1
    if self.scorer.reject(fit):
        return False  # Reject: score below threshold
    return True
```

A "placeholder" peak (intensity = 1.0) is inserted when no experimental peak matches
a theoretical peak position. If all but one experimental peak is a placeholder, the
fit is rejected for charge states > 1.

#### Layer 2: Score Thresholds

- MinimizeFitSelector: reject if score > minimum_score
- MaximizeFitSelector: reject if score < minimum_score
- PenalizedMSDeconVFitter default minimum_score = 10

#### Layer 3: Incremental Truncation

`fit_incremental_truncation()` generates progressively truncated versions of the
theoretical pattern and re-scores. This catches cases where trailing theoretical peaks
have no experimental counterpart but the core pattern matches well.

#### Layer 4: Peak Dependency Graph (exhaustive.py)

For complex spectra with overlapping isotope patterns, a peak dependency graph
resolves conflicts:
- Each isotopic fit registers its peak dependencies
- Connected components identify overlapping fits
- Three solver strategies:
  - **'top'**: Select single best fit per connected component
  - **'disjoint'**: Greedy selection of non-overlapping fits
  - **'iterative'**: Subtraction-based iterative resolution

#### Layer 5: Signal Subtraction

After a fit is accepted, its signal is subtracted from experimental peaks:

```python
def subtraction(self, isotopic_cluster, error_tolerance):
    for peak in isotopic_cluster:
        match = self.peaklist.has_peak(peak.mz, error_tolerance)
        if match is not None:
            existing = match.intensity
            match.intensity -= peak.intensity
            if (match.intensity < 0) or (peak.intensity > (existing * 0.7)):
                match.intensity = 1.  # Set to placeholder level
```

This prevents the same signal from being assigned to multiple isotope patterns.

#### Layer 6: M-1 Peak Handling

The `_find_previous_putative_peak()` method in base.py searches for peaks one isotopic
spacing BELOW the current seed peak, allowing the algorithm to check whether the
current peak is actually an M+1 peak of a lower-mass species rather than the monoisotopic:

```python
def _find_previous_putative_peak(self, mz, charge, step=1, tolerance=ERROR_TOLERANCE):
    shift = isotopic_shift(charge)
    prev_peak = mz - shift
    # Search for peaks at prev_peak position
    # If found, create candidate with corrected monoisotopic m/z
```

This generates alternative (peak, charge) candidates where the monoisotopic is shifted
left, which are then scored against the original assignment.

### 2.4 Feature Detection (Single Scan vs Integrated)

**File:** `src/ms_deisotope/feature_map/lcms_feature.py`
**File:** `src/ms_deisotope/feature_map/feature_fit.py`

#### Single-Scan Deconvolution

The core `deconvolute_peaks()` operates on a **single scan**. It processes one centroided
peak list at a time and returns deconvoluted peaks for that scan.

#### Multi-Scan Feature Assembly

The `feature_map` module provides LC-MS feature detection:

**LCMSFeature** class tracks peaks across retention time:
- `LCMSFeatureTreeNode`: Stores peaks at a single time point
- `LCMSFeatureTreeList`: Binary-search-indexed list of nodes across RT
- `split_sparse()`: Segments features with large RT gaps

**Feature Quality Scoring** (`feature_fit.py`):

Three multiplicative quality dimensions:

```python
def profile_qc(eic, averagine, truncate_after, smooth):
    v = 1.0
    v *= isotopic_consistency(eic, averagine, truncate_after)  # G-test of isotope pattern
    v *= spacing_fit(eic)        # Regular retention time spacing
    v *= shape_fit(eic, smooth)  # Chromatographic peak morphology
    return v
```

- **isotopic_consistency**: Uses ScaledGTestFitter to verify isotope pattern holds across scans
- **spacing_fit**: Penalizes irregular temporal spacing (suggests poor chromatographic peak)
- **shape_fit**: Uses AdaptiveMultimodalChromatogramShapeFitter for peak shape validation

#### Peak Retention Strategy

`TopNRetentionStrategy` handles unmatched peaks after deconvolution:
- Treats remaining peaks as potential monoisotopic peaks of unobserved distributions
- Applies mass threshold (default 850 Da) and intensity threshold (5% of base peak)
- Retains top N peaks (default 50) sorted by intensity

### 2.5 Theoretical Isotope Pattern Generation

**File:** `src/ms_deisotope/averagine.py`

Uses the Senko averagine model to interpolate composition from neutral mass:

```python
# Pre-defined averagine models
peptide = Averagine({"C": 4.9384, "H": 7.7583, "N": 1.3577, "O": 1.4773, "S": 0.0417})
glycopeptide = Averagine({"C": 10.93, "H": 15.75, "N": 1.6577, "O": 6.4773, "S": 0.02054})
glycan = Averagine({'C': 7.0, 'H': 11.8333, 'N': 0.5, 'O': 5.16666})
```

**Algorithm for `Averagine.scale(mz, charge)`:**
1. Calculate neutral mass from m/z and charge
2. Scale averagine composition by (neutral_mass / averagine_base_mass)
3. Round element counts to integers
4. Adjust hydrogen count to match target neutral mass

**Caching:** `AveragineCache` rounds m/z values and caches theoretical patterns, shifting
the cached pattern to the exact m/z on retrieval.

**Truncation:** `TheoreticalIsotopicPattern.truncate_after(fraction)` drops trailing peaks
that comprise the last (1-fraction) of total intensity (default: keeps 80%).

---

## 3. DeconTools / DeconEngineV2 (PNNL)

**Repository (DeconTools):** https://github.com/PNNL-Comp-Mass-Spec/DeconTools
**Repository (DeconEngineV2):** https://github.com/PNNL-Comp-Mass-Spec/DeconEngineV2
**Language:** C# (DeconTools), C++ (DeconEngineV2)
**License:** Apache 2.0
**Author:** PNNL Computational Mass Spectrometry group

### 3.1 Isotope Pattern Scoring

**File:** `DeconTools.Backend/ProcessingTasks/FitScoreCalculators/AreaFitter.cs`

The AreaFitter computes a normalized least-squares fit between theoretical and observed
isotope profiles:

```
fitScore = sum((normalized_theoretical_i - normalized_observed_i)^2) /
           sum(normalized_theoretical_i^2)
```

Where:
- Theoretical and observed data are each normalized by their own maximum intensity
- Observed intensities at theoretical m/z positions are obtained via linear interpolation
- Score of 0.0 = perfect match; 1.0 = complete mismatch; -1 = insufficient data
- Values > 0.15 are considered low quality

**File:** `DeconTools.Backend/ProcessingTasks/FitScoreCalculators/IsotopicProfileFitScoreCalculator.cs`

Orchestrates the scoring pipeline:
1. Find most abundant peak in theoretical profile
2. Find corresponding peak in observed profile (within 0.1 Da)
3. Calculate m/z offset between them
4. Generate theoretical XY profile using observed FWHM
5. Shift theoretical profile by the offset
6. Call AreaFitter.GetFit() for the actual scoring

This alignment step is crucial: it corrects for mass calibration drift before scoring.

### 3.2 Charge State Determination (THRASH/Horn Transform)

**File:** `DeconEngineV2/clsHornTransform.cpp`

Uses the THRASH (Thorough High Resolution Analysis of Spectra by Horn) algorithm:
1. Sort peaks by intensity (process strongest first)
2. For each peak, call `FindTransform()` which tries multiple charge states
3. Filter by maximum charge state constraint
4. Validate monoisotopic m/z when `IsActualMonoMZUsed` is enabled:
   - Calculate expected monoisotopic: `mono_mw / charge + 1.00727638`
   - Acceptance threshold: 20% below expected isotope spacing (1.003/charge)
   - Correct mass estimates when calculated and observed align

### 3.3 False Positive Rejection

- **Fit score threshold**: Fits > 0.15 are typically rejected as low quality
- **Interference score**: Measures likelihood of overlapping isotopic distributions
  (0 = no interference)
- **Minimum peptide intensity**: Peaks below threshold are skipped entirely
- **Charge state bounds**: Filters transforms exceeding MaxCharge parameter

### 3.4 Feature Detection

**Related tool:** https://github.com/PNNL-Comp-Mass-Spec/LC-IMS-MS-Feature-Finder
- Takes deisotoped features from DeconTools as input
- Groups across LC dimension and IMS dimension
- Specifically designed for LC-IMS-MS data (Waters Synapt-type instruments)

---

## 4. UniDec / IsoDec (michaelmarty/UniDec)

**Repository:** https://github.com/michaelmarty/UniDec
**Language:** Python + C (compiled DLLs)
**License:** See repository
**Author:** Michael Marty
**Description:** Universal Deconvolution of Mass and Ion Mobility Spectra

### 4.1 Isotope Pattern Scoring

**File:** `unidec/IsoDec/match.py`

Primary scoring metric is **cosine similarity** (CSS):

```python
@njit(fastmath=True)
def calculate_cosinesimilarity(cent_intensities, iso_intensities, shift, max_shift,
                               minusoneaszero=True):
    # Returns: ab / (sqrt(a^2) * sqrt(b^2))
    # Normalized dot product between centroid and theoretical intensities
```

Key features:
- Accounts for missed monoisotopic peaks via the `shift` parameter
- `minusoneaszero=True`: Treats missing M-1 peaks as zero intensity
- Typical acceptance threshold: CSS >= 0.7

**Pattern matching validation:**

```python
@njit(fastmath=True)
def isodist_match(isodists1, isodist2, css_threshold=0.7, ppm_tol=20):
    css = calc_css_from_data(isodists1, isodist2, ppm_tol=ppm_tol)
    if css < css_threshold:
        return False
    return True
```

### 4.2 Charge State Determination -- NEURAL NETWORK APPROACH

**File:** `unidec/IsoDec/encoding.py`

UniDec/IsoDec uses a **novel phase-encoding** approach with a neural network:

```python
@njit(fastmath=True)
def encode_phase(centroids, maxz=50, phaseres=8):
    """
    Encode charge phases as a 2D histogram (maxz x phaseres).
    For each candidate charge z, compute:
        phase = (m/z * z / mass_diff_c) % 1
    Accumulate intensity into phase bins.
    """
    phases = np.zeros((maxz, phaseres))
    rescale = centroids[:, 0] / mass_diff_c
    for i in range(maxz):
        phase = (rescale * (i + 1)) % 1
        phaseindexes = np.floor(phase * phaseres)
        for j in range(len(centroids)):
            phases[i, int(phaseindexes[j])] += centroids[j, 1]
    phases /= np.amax(phases)
    return phases
```

**Concept:** If peaks are spaced at 1/z Da intervals (the isotope spacing for charge z),
then multiplying m/z by z and taking modulo (mass_diff_c = 1.0033) will align all peaks
to the same phase. The neural network learns to recognize which charge state produces
the most coherent phase alignment.

**Prediction modes:**
- Mode 0: Batch neural network prediction
- Mode 1: Per-peak neural network prediction
- Mode 2: THRASH algorithm (fallback)
- Mode 3: Ensemble of neural network + THRASH
- Mode 4: Full probability vector for per-charge scoring

**Pre-trained models:** `phase_model_1.pth`, `phase_model_4.pth`, `phase_model_8.pth`
(different phase resolutions).

### 4.3 False Positive Rejection

**File:** `unidec/IsoDec/match.py`

Multi-stage validation in `make_shifted_peak()`:

1. **Minimum peaks threshold:** `len(matchedindexes) >= minpeaks`
2. **CSS threshold:** `shiftscore >= css_thresh` (configurable, default ~0.7)
3. **Area coverage:** `areacovered > minareacovered` OR top-three peaks condition
4. **Charge-specific validation for z=1:**
   ```python
   if z == 1:
       if len(matchedindexes) == 2:
           if isomatches[0] == 0 and isomatches[1] == 1:
               ratio = int2 / int1
               if p1low < ratio < p1high:
                   minpeaks = 2  # Accept 2-peak pattern if ratio is valid
   ```
5. **Adaptive shift limits by charge:**
   ```python
   if z < 3:
       maxshift = 1    # Low charge: allow 1 missed monoisotopic
   elif z < 6:
       maxshift = 2    # Mid charge: allow 2
   else:
       maxshift = config.maxshift  # High charge: configurable
   ```
6. **Mass proximity check for merging:**
   ```python
   close_check = abs(pk.monoiso - nearest_mass) <= config.maxshift * config.mass_diff_c * 1.1
   ```

**Missed monoisotopic handling:**

```python
@njit(fastmath=True)
def within_ppm_plus_mm(mass1, mass2, ppm_tol=20, max_mm=1, mass_diff_c=1.0033):
    # Tests mass compatibility allowing for missed monoisotopic peaks
    # Adjusts mass2 by mm * 1.0033 for mm in range(-max_mm, max_mm+1)
```

### 4.4 Feature Detection

IsoDec processes spectra **scan-by-scan** with iterative knockdown:

1. **Knockdown rounds**: Process spectrum multiple times with decreasing sensitivity
2. **Peak removal**: Successfully matched isotopic patterns are subtracted
3. **Remaining signal**: Processed in subsequent rounds to find weaker features

For LC-MS data, `process_file()` iterates across scans and `pks_to_mass()` converts
matched peaks to a mass spectrum.

### 4.5 Isotope Distribution Generation

**File:** `unidec/modules/isotopetools.py`

Two methods:

1. **FFT convolution** (`isojim`): Uses pre-computed FFT transforms for each element
   (H, C, N, O, S, metals) and convolves them based on molecular formula

2. **Fast approximation** (`isomike`): Exponential + Gaussian model, Numba-accelerated:
   ```python
   @njit(fastmath=True)
   def isomike(mass, length=128):
       # Combined exponential decay and Gaussian envelope model
   ```

3. **Neural network** (IsoGen): Trained model for fast isotope distribution prediction

---

## 5. Comparative Summary

| Feature | ms_deisotope | DeconTools | UniDec/IsoDec | CLIPS-MS |
|---------|-------------|------------|---------------|----------|
| **Scoring metric** | PenalizedMSDeconV (MSDeconV * G-test penalty) | Normalized least-squares (AreaFitter) | Cosine similarity | Spectral Accuracy (RMSE) |
| **Score range** | Higher = better (maximizing) | 0-1, lower = better | 0-1, higher = better | 0-100%, higher = better |
| **Charge determination** | Exhaustive search over range | THRASH (Horn transform) | Neural network (phase encoding) | N/A (formula ID tool) |
| **False positive layers** | 6 layers (placeholder, score, truncation, dependency graph, subtraction, M-1) | Fit threshold + interference + charge bounds | CSS threshold + peak count + area coverage + charge-specific rules | Spectral Accuracy threshold |
| **Feature detection** | LCMSFeature (multi-scan, RT + shape + isotope QC) | LC-IMS-MS Feature Finder (separate tool) | Scan-by-scan with knockdown | N/A |
| **Isotope model** | Averagine (Senko) with caching | Averagine + Mercury algorithm | FFT convolution / IsoGen neural net | TrueCal calibrated line shape |
| **Language** | Python + C extensions | C# + C++ | Python + C DLLs + PyTorch | Proprietary |
| **Open source** | Yes (Apache 2.0) | Yes (Apache 2.0) | Yes | No |

---

## 6. Key Algorithms Worth Adopting

### From ms_deisotope:
- **PenalizedMSDeconV scoring**: The combination of per-peak scoring (mass accuracy *
  abundance match * sqrt(intensity)) with a G-test penalty for overall distribution shape
  is the gold standard in the field
- **Peak dependency graph**: Essential for complex spectra with overlapping isotope patterns
- **Incremental truncation**: Smart handling of partially-observed isotope patterns
- **Signal subtraction with 70% guard**: Prevents over-subtraction

### From UniDec/IsoDec:
- **Phase encoding for charge detection**: Novel and fast -- could replace exhaustive search
  for initial charge state estimation
- **Cosine similarity scoring**: Simpler and potentially more robust than G-test for
  isotope pattern matching, especially with missed peaks
- **Adaptive maxshift by charge**: Practical heuristic for missed monoisotopic handling

### From DeconTools:
- **Profile-mode AreaFitter with interpolation**: When working with profile (not centroided)
  data, interpolating observed intensities at theoretical positions gives more accurate scores
- **Alignment before scoring**: Correcting m/z offset using the most abundant peak before
  scoring prevents calibration drift from corrupting fit scores
