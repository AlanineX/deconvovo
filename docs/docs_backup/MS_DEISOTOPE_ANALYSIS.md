# ms_deisotope Source Code Analysis

Deep-dive analysis of the `ms_deisotope` Python library by Joshua Klein (mobiusklein),
based on reading the installed source at:
`.venv-ms/lib/python3.12/site-packages/ms_deisotope/`

---

## 1. How Does It Score Isotope Patterns?

**File:** `ms_deisotope/scoring.py`

ms_deisotope provides **7 different scoring functions**, all inheriting from `IsotopicFitterBase`.
Each takes a list of observed experimental peaks and a list of matching theoretical peaks, and
returns a numeric score. Scorers are either **minimizing** (lower = better, e.g. G-test, least
squares) or **maximizing** (higher = better, e.g. MSDeconV).

### 1a. ScaledGTestFitter (the G-test) -- DEFAULT COMPONENT

The G-test (log-likelihood ratio test) after normalizing both distributions to sum to 1.0:

```python
class ScaledGTestFitter(IsotopicFitterBase):
    def evaluate(self, peaklist, observed, expected, **kwargs):
        total_observed = sum(p.intensity for p in observed)
        total_expected = sum(p.intensity for p in expected)
        total_expected += eps  # eps = 1e-4, prevents division by zero
        normalized_observed = [obs.intensity / total_observed for obs in observed]
        normalized_expected = [theo.intensity / total_expected for theo in expected]
        g_score = 2 * sum([obs * np.log(obs / theo)
                           for obs, theo in zip(normalized_observed, normalized_expected)])
        return g_score
```

This is a **minimizing** scorer. G = 0 means perfect match. The formula is:

    G = 2 * SUM( o_i * ln(o_i / e_i) )

where o_i and e_i are the normalized intensities of the ith observed and expected peaks.

### 1b. MSDeconVFitter -- MAXIMIZING scorer

An implementation of the scoring from the MSDeconV paper (Liu et al., 2010). Scores each
peak individually based on mass accuracy and abundance agreement:

```python
def score_peak(self, obs, theo, mass_error_tolerance=0.02, minimum_signal_to_noise=1):
    if obs.signal_to_noise < minimum_signal_to_noise:
        return 0.

    mass_error = np.abs(obs.mz - theo.mz)
    if mass_error <= mass_error_tolerance:
        mass_accuracy = 1 - mass_error / mass_error_tolerance
    else:
        mass_accuracy = 0

    if obs.intensity < theo.intensity and (((theo.intensity - obs.intensity) / obs.intensity) <= 1):
        abundance_diff = 1 - ((theo.intensity - obs.intensity) / obs.intensity)
    elif obs.intensity >= theo.intensity and (((obs.intensity - theo.intensity) / obs.intensity) <= 1):
        abundance_diff = np.sqrt(1 - ((obs.intensity - theo.intensity) / obs.intensity))
    else:
        abundance_diff = 0.

    score = np.sqrt(theo.intensity) * mass_accuracy * abundance_diff
    return score
```

Key points:
- `mass_accuracy` is a linear penalty: 1.0 at perfect match, 0.0 at the tolerance boundary
- `abundance_diff` penalizes asymmetrically: if observed < theoretical, linear penalty; if observed > theoretical, square-root penalty (more tolerant of over-estimation)
- Score is weighted by `sqrt(theoretical_intensity)` -- heavier isotope peaks contribute more
- Total score = sum of per-peak scores

### 1c. PenalizedMSDeconVFitter -- THE ACTUAL DEFAULT SCORER

This is the **default scorer used by `AveraginePeakDependenceGraphDeconvoluter`** (and therefore by `deconvolute_peaks()`). It combines MSDeconV with a G-test penalty:

```python
class PenalizedMSDeconVFitter(IsotopicFitterBase):
    def __init__(self, minimum_score=10, penalty_factor=1., mass_error_tolerance=0.02):
        self.select = MaximizeFitSelector(minimum_score)
        self.msdeconv = MSDeconVFitter(mass_error_tolerance=mass_error_tolerance)
        self.penalizer = ScaledGTestFitter()
        self.penalty_factor = penalty_factor

    def evaluate(self, peaklist, observed, expected, **kwargs):
        score = self.msdeconv.evaluate(peaklist, observed, expected)
        penalty = abs(self.penalizer.evaluate(peaklist, observed, expected))
        return score * (1 - penalty * self.penalty_factor)
```

The formula is:

    S(e, t) = MSDeconV(e, t) * (1 - |G_test(e, t)| * penalty_factor)

This is a **maximizing** scorer with a default minimum_score of 10. It takes the MSDeconV
score and multiplies it by a shape-penalty derived from the G-test. If the shape is
perfect (G ~ 0), the penalty is ~0 and the full MSDeconV score is retained. If the
shape is poor (G is large), the score is heavily reduced.

### 1d. Other Scorers

| Scorer | Type | Formula |
|--------|------|---------|
| `GTestFitter` | minimizing | Raw G-test on unnormalized intensities |
| `ChiSquareFitter` | minimizing | SUM((obs - theo)^2 / theo) |
| `LeastSquaresFitter` | minimizing | SUM((norm_obs - norm_theo)^2) / SUM(norm_theo^2), where both are max-normalized |
| `DistinctPatternFitter` | minimizing | G_test * (interference + eps) / (npeaks * scale) * domain_scale |
| `DotProductFitter` | neither | SUM(obs_i * theo_i) |

### 1e. InterferenceDetection

Not a scorer per se, but used by `DistinctPatternFitter`. Measures how much of the total
signal in the m/z region of the fit is NOT accounted for by the fit's peaks:

```python
score = 1 - (included_intensity / region_intensity)
```

A score of 0 means no interference; a score near 1 means the fit's peaks account for
very little of the total signal in that region.

---

## 2. How Does It Handle M-1 / M-2 Peaks?

**Files:** `ms_deisotope/deconvolution/exhaustive.py`, `ms_deisotope/deconvolution/base.py`,
`ms_deisotope/deconvolution/utils.py`

### 2a. Left-Search: Looking for Peaks BELOW the Monoisotopic

ms_deisotope does NOT explicitly check for "M-1 peaks" as contamination. Instead, it uses
a **left_search_limit** parameter to consider whether the current peak might actually be
an M+1 or M+2 peak of a LOWER monoisotopic peak.

In `_get_all_peak_charge_pairs()` (exhaustive.py):

```python
for charge in charges:
    target_peaks.add((peak, charge))

    # Look Left
    for i in range(1, left_search_limit):
        prev_peak = has_previous_peak_at_charge(self, peak, charge, i)
        if prev_peak is None:
            continue
        target_peaks.add((prev_peak, charge))
        if recalculate_starting_peak:
            target_peaks.update(self._find_previous_putative_peak(
                peak.mz, charge, i, 2 * error_tolerance))

    # Look Right
    for i in range(1, right_search_limit):
        nxt_peak = has_successor_peak_at_charge(self, peak, charge, i)
        ...
```

The `has_previous_peak_at_charge` function (utils.py):

```python
def has_previous_peak_at_charge(peak_collection, peak, charge=2, step=1,
                                 error_tolerance=ERROR_TOLERANCE):
    prev = peak.mz - isotopic_shift(charge) * step
    return peak_collection.has_peak(prev, error_tolerance)
```

### 2b. What This Means

When examining peak P:
- **left_search_limit=3** (default for simple deconvolution): It checks if there are peaks
  at P - 1/z, P - 2/z, and P - 3/z. If found, those lower peaks become candidate
  monoisotopic peaks, and P is considered as potentially their M+1, M+2, or M+3.
- **left_search_limit=1** (default for PeakDependenceGraph): It checks one position to
  the left.

### 2c. NO Explicit M-1 Penalty

ms_deisotope does **NOT** explicitly penalize for the presence of an M-1 peak (a peak
one neutron BELOW the monoisotopic). Unlike THRASH which explicitly checks for and
rejects patterns where an M-1 peak is detected, ms_deisotope takes a different approach:

1. It generates a theoretical pattern starting from the putative monoisotopic peak
2. It matches experimental peaks to this theoretical pattern
3. The scoring function evaluates goodness-of-fit

If you assign the wrong peak as monoisotopic (e.g., the actual M+1 peak), the theoretical
pattern will have poor shape agreement with the observed peaks, producing a bad score. The
correct assignment (where the actual monoisotopic peak is one position to the left) will
score better and be selected.

The key defense mechanism is the **left search**: for each peak encountered, the algorithm
also tries assigning peaks to its left as the monoisotopic peak. The scoring then
naturally selects the correct assignment.

### 2d. Single-Peak Rejection

There IS a hard filter: if a fit has only one real (non-placeholder) peak AND the charge
state is > 1, the fit is rejected:

```python
def _check_fit(self, fit):
    if len(drop_placeholders(fit.experimental)) == 1 and fit.charge > 1:
        return False
    if self.scorer.reject(fit):
        return False
    return True
```

---

## 3. How Does It Determine Charge State?

**Files:** `ms_deisotope/deconvolution/exhaustive.py`, `ms_deisotope/deconvolution/utils.py`

### 3a. Brute-Force Enumeration + Best Score Selection

The primary strategy is straightforward: **try every charge state in the range and pick
the one with the best isotopic fit score**.

In `charge_state_determination()` (exhaustive.py):

```python
def charge_state_determination(self, peak, error_tolerance=ERROR_TOLERANCE,
                                charge_range=(1, 8), ...):
    results = self._fit_all_charge_states(
        peak, error_tolerance=error_tolerance, charge_range=charge_range, ...)
    try:
        result = self.scorer.select(results)
        return result
    except ValueError:
        return None
```

`_fit_all_charge_states` calls `_get_all_peak_charge_pairs` which iterates over all
charges from max down to min:

```python
charges = ChargeIterator(*charge_range)  # e.g., (1, 8) -> [8, 7, 6, 5, 4, 3, 2, 1]
for charge in charges:
    target_peaks.add((peak, charge))
    # also look left and right for each charge...
```

For each (peak, charge) candidate, `fit_theoretical_distribution()` generates the
averagine-based theoretical pattern at that charge, finds matching experimental peaks,
and scores the fit. Then `scorer.select()` picks the best.

### 3b. How z=2 vs z=4 Is Resolved

For a peak at m/z = 500:
- At z=2: isotope spacing = ~1.003/2 = ~0.5015 Da. It looks for peaks at 500.50, 501.00, 501.50...
- At z=4: isotope spacing = ~1.003/4 = ~0.2508 Da. It looks for peaks at 500.25, 500.50, 500.75...

Both z=2 and z=4 will find the peak at ~500.50 (it's the M+1 for z=2 and M+2 for z=4).
The distinction comes from:

1. **Additional peaks**: z=4 expects peaks at ~500.25 and ~500.75 that z=2 does not
2. **Theoretical intensity distribution**: The averagine model produces DIFFERENT isotopic
   envelopes for the neutral masses implied by z=2 vs z=4 (because neutral_mass = (mz * z) - z * proton_mass)
3. **Scoring**: The fit with the better shape match wins

### 3c. QuickCharge Algorithm (Optional)

An implementation of Hoopmann's QuickCharge algorithm (2007) for pre-screening feasible
charge states by looking at neighboring peaks:

```python
def quick_charge(peak_set, index, min_charge, max_charge):
    min_intensity = peak_set[index].intensity / 4.
    charges = np.zeros(max_charge, dtype=int)
    for j in range(index + 1, size):
        if peak_set[j].intensity < min_intensity:
            continue
        diff = peak_set[j].mz - peak_set[index].mz
        if diff > 1.1:
            break
        raw_charge = 1 / diff
        charge = int(raw_charge + 0.5)
        remain = raw_charge - int(raw_charge)
        if 0.2 < remain < 0.8:  # fractional part too far from integer
            continue
        if charge < min_charge or charge > max_charge:
            continue
        charges[charge] = 1
    # Also search backwards...
    return np.where(charges)[0]
```

Logic:
- Look at neighboring peaks within 1.1 Da
- Calculate `1 / (mz_difference)` as the candidate charge
- If the fractional part is between 0.2 and 0.8, reject (not close enough to integer)
- Mark feasible charges

This is **disabled by default** (`use_quick_charge=False`). When enabled, it pre-filters
the charge states to try, reducing computation.

### 3d. This is NOT THRASH

Despite conceptual similarity to THRASH (Horn et al., 2000), ms_deisotope differs:

- **THRASH** walks the spectrum finding the most abundant peak, then tries to assign it
  as part of an isotopic cluster by checking spacing patterns and using mercury-based
  theoretical patterns.
- **ms_deisotope** visits every peak, generates all candidate (monoisotopic_peak, charge)
  pairs via left/right search, fits averagine-based theoretical patterns, scores them,
  and resolves conflicts via a peak dependency graph.

The peak dependency graph approach is the key differentiator: instead of greedy
assignment, it builds a hypergraph of all possible fits and solves for the globally
optimal disjoint set.

---

## 4. How Does It Compute Theoretical Isotope Distributions?

**File:** `ms_deisotope/averagine.py`

### 4a. The Averagine Model

ms_deisotope uses the **averagine** approach from Senko et al. (1995). An "averagine" is
the average elemental composition of a monomer unit for a class of biomolecule.

Pre-defined averagines:

```python
peptide = Averagine({"C": 4.9384, "H": 7.7583, "N": 1.3577, "O": 1.4773, "S": 0.0417})
glycopeptide = Averagine({"C": 10.93, "H": 15.75, "N": 1.6577, "O": 6.4773, "S": 0.02054})
glycan = Averagine({'C': 7.0, 'H': 11.8333, 'N': 0.5, 'O': 5.16666})
permethylated_glycan = Averagine({'C': 12.0, 'H': 21.8333, 'N': 0.5, 'O': 5.16666})
heparin = Averagine({'H': 10.5, 'C': 6, 'S': 0.5, 'O': 5.5, 'N': 0.5})
heparan_sulfate = Averagine({'H': 10.667, 'C': 6.0, 'S': 1.333, 'O': 9.0, 'N': 0.667})
```

### 4b. Scaling to Target Mass

The `Averagine.scale()` method interpolates the composition for a given m/z and charge:

```python
def scale(self, mz, charge=1, charge_carrier=PROTON):
    neutral = neutral_mass(mz, charge, charge_carrier)
    scale = neutral / self.base_mass
    scaled = {}
    for elem, count in self.base_composition.items():
        scaled[elem] = round(count * scale)
    # Hydrogen correction to match target mass
    scaled_mass = calculate_mass(scaled)
    delta_hydrogen = round(scaled_mass - neutral)
    H = scaled["H"]
    if H > delta_hydrogen:
        scaled["H"] = H - delta_hydrogen
    else:
        scaled["H"] = 0
    return scaled
```

Steps:
1. Convert m/z + charge to neutral mass
2. Divide neutral mass by average monomer mass to get scaling factor
3. Multiply each element count by the scaling factor, round to integer
4. Calculate the mass of the rounded composition
5. Adjust hydrogen count to minimize mass error (hydrogen correction)

### 4c. Isotopic Distribution: brainpy

The actual isotope distribution computation is delegated to the **brainpy** library
(also by mobiusklein). This is a Python/Cython implementation of the Baffling Recursive
Algorithm for Isotopic distributioN calculations (BRAIN).

```python
from brainpy import isotopic_variants

def isotopic_cluster(self, mz, charge=1, ...):
    composition = self.scale(mz, charge, charge_carrier)
    peaklist = isotopic_variants(composition, charge=charge)  # <-- brainpy call
    tid = TheoreticalIsotopicPattern(peaklist, peaklist[0].mz, 0)
    tid.shift(mz)
    if truncate_after < 1.0:
        tid.truncate_after(truncate_after)
    if ignore_below > 0:
        tid.ignore_below(ignore_below)
    return tid
```

`brainpy.isotopic_variants()` computes the exact isotopic distribution for the given
elemental composition using polynomial expansion of isotopic probabilities. This is
NOT the Poisson approximation or the FFT-based Mercury algorithm. It is an exact
combinatorial method (the BRAIN algorithm) with optional Cython acceleration.

### 4d. AveragineCache

To avoid recomputing isotopic patterns for similar masses, `AveragineCache` rounds
the m/z to the nearest truncation unit (default: 1.0 Da) and caches the result:

```python
def has_mz_charge_pair(self, mz, charge=1, ...):
    key_mz = round(mz / self.cache_truncation) * self.cache_truncation
    if (key_mz, charge, charge_carrier) in self.backend:
        return self.backend[key_mz, charge, charge_carrier].clone().shift(mz)
    else:
        tid = self.averagine.isotopic_cluster(key_mz, charge, ...)
        self.backend[key_mz, charge, charge_carrier] = tid.clone()
        return tid
```

The cached pattern is cloned and shifted to the exact target m/z.

### 4e. Pattern Truncation and Filtering

Two post-processing steps are applied to theoretical patterns:

1. **truncate_after** (default 0.95): Drops trailing peaks such that the retained peaks
   account for 95% of the total signal. The pattern is then renormalized.

2. **ignore_below** (default 0.001): Drops peaks with intensity below 0.1% of the
   total. The pattern is renormalized after removal.

### 4f. Scaling Theoretical to Match Experimental

Before scoring, the theoretical pattern is scaled to match experimental intensities.
The `scale_method` (default "sum") determines how:

- **"sum"**: Scale so that SUM(theoretical) = SUM(experimental)
- **"max"**: Scale so that the most abundant theoretical peak matches its experimental counterpart
- **"meanscale"**: Weighted average scaling factor, weighted by theoretical * experimental^2
- **"top3"**: Average of scale factors from the 3 most abundant theoretical peaks

### 4g. Poisson Approximation (Fallback)

There is also a fast Poisson approximation available:

```python
def _poisson_approximate(mass, n_peaks, charge=1):
    lmbda = mass / 1800.0
    # ... computes Poisson distribution with lambda = mass/1800
```

This is not used in the main deconvolution pipeline but is available as a utility.

---

## 5. How Does It Handle Feature Detection (Multi-Scan)?

**Files:** `ms_deisotope/feature_map/feature_processor.py`, `ms_deisotope/feature_map/feature_fit.py`,
`ms_deisotope/feature_map/lcms_feature.py`, `ms_deisotope/feature_map/dependence_network.py`

### 5a. Two Separate Systems

ms_deisotope has **two distinct deconvolution systems**:

1. **Single-scan deconvolution** (`deconvolution/` package): Works on one spectrum at a time.
   This is `deconvolute_peaks()` and the `ScanProcessor` class.

2. **LC-MS feature-level deconvolution** (`feature_map/` package): Works on pre-extracted
   LC-MS features (chromatographic peaks that span multiple scans).

### 5b. Single-Scan (ScanProcessor)

`ScanProcessor` (processor.py) iterates scan-by-scan. Each scan is independently
peak-picked and deconvoluted. There is no cross-scan integration at this level.

The multi-pass iteration within a single scan is controlled by the `iterations` parameter
(default 10). Each iteration:
1. Builds the peak dependency graph
2. Solves it for disjoint best fits
3. Subtracts assigned signal
4. Checks convergence: if `(signal_before - signal_after) / signal_after < convergence`, stop

### 5c. LC-MS Feature Deconvolution (LCMSFeatureProcessor)

`LCMSFeatureProcessor` works on **pre-extracted LC-MS features** -- chromatographic traces
that have already been detected across multiple scans. The input is an `LCMSFeatureMap`
containing `LCMSFeature` objects (chromatographic peaks with m/z, time, and intensity).

The algorithm:

1. For each feature (chromatographic trace), try all charge states
2. For each charge, generate a theoretical isotopic pattern at the feature's m/z
3. Find matching features at the theoretical isotopic peaks' m/z values
4. Use a `FeatureSetIterator` to walk through time points where all features overlap
5. At each time point, extract the current intensity of each isotopic feature, score the
   fit using the same scorers (default: `g_test_scaled`)
6. Aggregate scores across time points using a thresholded mean

```python
def _fit_feature_set(self, mz, error_tolerance, charge, ...):
    base_tid = self.create_theoretical_distribution(mz, charge, ...)
    feature_groups = self.match_theoretical_isotopic_distribution(base_tid, error_tolerance)
    for features in product(*feature_groups):
        feat_iter = FeatureSetIterator(features)
        scores = []
        for eid in feat_iter:  # iterate through time points
            cleaned_eid, tid, n_missing = self.conform_envelopes(eid, base_tid)
            score = self.scorer.evaluate(None, cleaned_eid, tid)
            scores.append(score)
        final_score = self._find_thresholded_score(scores, threshold_scale)
        fit = LCMSFeatureSetFit(features, base_tid, final_score, charge, ...)
```

The `_find_thresholded_score` takes the mean of scores above `threshold_scale * max(scores)`:

```python
def _find_thresholded_score(self, scores, percentage):
    scores = np.array(scores)
    maximum = scores.max()
    threshold = maximum * percentage
    return scores[scores > threshold].mean()
```

### 5d. Feature Dependency Graph

Like single-scan deconvolution, the LC-MS feature processor uses a dependency graph
(`FeatureDependenceGraph`) to resolve overlapping isotopic patterns among features.
It finds disjoint best fits and finalizes them.

### 5e. Finalization: Scan-by-Scan Deconvoluted Peaks

During finalization, each accepted feature fit is iterated time-point by time-point,
creating individual `DeconvolutedPeak` instances for each time point wrapped in
`DeconvolutedLCMSFeatureTreeNode` objects. These are assembled into a
`DeconvolutedLCMSFeature` representing the full LC-MS feature.

### 5f. Ion Mobility Frame Processor

`IonMobilityFrameProcessor` (mobility_frame_processor.py) extends the same pattern for
ion mobility data. It processes each frame (which contains mobility-separated spectra)
and supports `IonMobilityDeconvolutedPeak` objects that carry a `drift_time` attribute.

---

## 6. Peak Dependency Network: The Core Innovation

**Files:** `ms_deisotope/peak_dependency_network/peak_network.py`,
`ms_deisotope/peak_dependency_network/subgraph.py`

### 6a. The Problem

In a complex mass spectrum, multiple isotopic pattern hypotheses can claim the same
experimental peak. A greedy approach (first-come-first-served) may not find the global
optimum.

### 6b. The Solution: Hypergraph

ms_deisotope builds a **peak dependency graph** (actually a hypergraph):
- **Nodes**: individual experimental peaks (`PeakNode`)
- **Hyperedges**: isotopic fits (`IsotopicFitRecord`), each connecting multiple peaks

Each `PeakNode` maintains a `links` dictionary mapping each fit that uses this peak to
that fit's score.

### 6c. Clustering and Solving

After all fits are inserted:

1. **Find connected components** (`_gather_independent_clusters`): Group all fits that
   share any peak into `DependenceCluster` objects
2. **Solve each cluster independently**: Three solver strategies are available:
   - **"top"** (default): Just take the single best fit from each cluster
   - **"disjoint"**: Find the greedy best set of non-overlapping fits (using `GreedySubgraphSelection`)
   - **"iterative"**: Iteratively pick the best fit, subtract its signal, re-score
     overlapping fits, and repeat

### 6d. Greedy Disjoint Subset Selection

The `GreedySubgraphSelection` class sorts fits by score and greedily assigns them to
non-overlapping layers:

```python
def layout_layers(envelopes, overlap_fn=peak_overlap, maximize=True):
    layers = [[]]
    envelopes.sort(key=lambda x: x.score, reverse=maximize)
    for envelope in envelopes:
        for layer in layers:
            collision = False
            for member in layer:
                if overlap_fn(envelope, member):
                    collision = True
                    break
            if not collision:
                layer.append(envelope)
                placed = True
                break
        if not placed:
            layers.append([envelope])
    return layers
```

The first layer contains the best non-overlapping subset.

Two fits "overlap" if they share any peak index (`peak_overlap` function):

```python
def peak_overlap(a, b):
    return len(a.peak_indices & b.peak_indices) > 0
```

---

## 7. Complete Deconvolution Pipeline Summary

The default call `deconvolute_peaks(peaklist)` executes:

1. **Prepare peaklist**: Clone and reindex the centroided peak list
2. **Construct deconvoluter**: `AveraginePeakDependenceGraphDeconvoluter` with:
   - Averagine model: `peptide` (C4.9384 H7.7583 N1.3577 O1.4773 S0.0417)
   - Scorer: `PenalizedMSDeconVFitter(minimum_score=10, penalty_factor=1.0)`
   - Charge range: (1, 8)
   - Error tolerance: 20 ppm
3. **For each iteration** (up to 10, checking convergence):
   a. **Populate graph**: For each peak in the spectrum:
      - For each charge in range [8, 7, ..., 1]:
        - Search left (1 position) and right (0 positions) for alternative monoisotopic peaks
        - For each (candidate_monoisotopic_peak, charge) pair:
          - Generate averagine theoretical pattern at that m/z and charge
          - Match theoretical peaks to experimental peaks within 20 ppm
          - Scale theoretical to match experimental (sum method)
          - Score with PenalizedMSDeconV
        - Insert best fits into peak dependency graph
   b. **Solve graph**: Find connected components, take best fit from each
   c. **Subtract**: Remove assigned signal from peaks
   d. **Check convergence**: If signal change < 0.1%, stop
4. **Merge isobaric peaks**: Combine peaks with same mass and charge
5. **Return** `DeconvolutedPeakSet`

---

## 8. Key Constants

```python
TRUNCATE_AFTER = 0.95    # Keep 95% of theoretical pattern signal
MAX_ITERATION = 10       # Maximum deconvolution passes
ERROR_TOLERANCE = 2e-5   # 20 ppm mass tolerance
IGNORE_BELOW = 0.001     # Ignore theoretical peaks below 0.1%
CONVERGENCE = 1e-3       # Convergence threshold (0.1% signal change)
SCALE_METHOD = "sum"     # Scale theoretical to match experimental by sum
```

The isotopic shift (distance between isotope peaks) is computed from the C13-C12 mass
difference:

```python
_neutron_shift = calculate_mass({"C[13]": 1}) - calculate_mass({"C[12]": 1})  # ~1.00335 Da

def isotopic_shift(charge=1):
    return _neutron_shift / float(charge)
```

---

## 9. Key Source Files Reference

| File | Purpose |
|------|---------|
| `scoring.py` | All scoring functions (GTest, MSDeconV, PenalizedMSDeconV, etc.) |
| `averagine.py` | Averagine models, isotopic pattern generation, caching |
| `deconvolution/base.py` | `DeconvoluterBase` -- matching, scaling, subtraction |
| `deconvolution/exhaustive.py` | `ExhaustivePeakSearchDeconvoluterBase`, `PeakDependenceGraphDeconvoluterBase` |
| `deconvolution/averagine_based.py` | `AveragineDeconvoluter`, `AveraginePeakDependenceGraphDeconvoluter` |
| `deconvolution/utils.py` | `quick_charge`, `has_previous_peak_at_charge`, `ChargeIterator` |
| `deconvolution/api.py` | `deconvolute_peaks()` -- the high-level API |
| `peak_dependency_network/peak_network.py` | `PeakDependenceGraph`, `DependenceCluster` |
| `peak_dependency_network/subgraph.py` | `GreedySubgraphSelection`, `ConnectedSubgraph` |
| `peak_set.py` | `DeconvolutedPeak`, `DeconvolutedPeakSet`, `Envelope` |
| `envelope_statistics.py` | `a_to_a2_ratio`, `most_abundant_mz`, `average_mz` |
| `feature_map/feature_processor.py` | `LCMSFeatureProcessor` -- multi-scan feature deconvolution |
| `feature_map/feature_fit.py` | `LCMSFeatureSetFit`, `DeconvolutedLCMSFeature` |
| `feature_map/lcms_feature.py` | `LCMSFeature` chromatographic feature representation |
| `processor.py` | `ScanProcessor` -- scan-by-scan pipeline |
| `constants.py` | Default parameter values |
