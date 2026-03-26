# HDX Reconstruction: Findings and Implementation Plan

Research consolidated: 2026-03-25

This document consolidates findings from our HDX-MS reconstruction research
and proposes a versioned implementation plan to fix the broken deuterium
uptake pipeline. Each version is designed to be independently testable.

---

## 1. Current State (Broken)

**Script:** `skills/public/hdx-process-extract/scripts/reconstruct_hdx_from_raw.py`

### The Bug

The function `best_feature_for_peptide()` computes a `neutral_mass_est` by
subtracting isotope index offsets from each matched peak before averaging:

```python
iso_idx = np.arange(matched_int.size, dtype=float)
used = matched_int > 0
neutral_series = (matched_mz[used] * z - z * PROTON) - iso_idx[used] * NEUTRON
neutral_mass_est = float(np.average(neutral_series, weights=matched_int[used]))
```

This collapses every isotopic peak back to an estimate of the monoisotopic
mass. The `- iso_idx[used] * NEUTRON` term strips out the deuterium mass
shift, which is the very signal we need to measure. The result is that
`neutral_mass_est` approximates the theoretical monoisotopic mass regardless
of how much deuterium has been incorporated.

### Why Exchange Ratios Are Nonsensical

Because the script measures monoisotopic mass (not centroid mass), the
`delta_mass_vs_baseline` between timepoints reflects only mass accuracy noise,
not deuterium incorporation. Small random fluctuations in mass accuracy produce
deltas of a few hundredths of a Da in both directions. When `finalize_metrics()`
divides these random deltas by the (also-random) delta of the longest timepoint,
the resulting `exchange_ratio_vs_max` values scatter wildly -- observed range
was [-36, +11] instead of the expected [0, 1].

### Additional Issues

| Issue | Description |
|-------|-------------|
| No centroid mass | Computes monoisotopic mass estimate, not envelope centroid |
| No undeuterated control | Uses `peptide_mapping` run implicitly, no explicit validation |
| No fully deuterated control | Normalizes to longest timepoint, not a maxD control |
| No back-exchange correction | All uptake values would be systematically underestimated |
| No exchangeable amide count | Does not compute N (prolines, N-terminal) |
| Single scan only | Picks best-scoring scan, ignores chromatographic peak |
| No RT constraint | Searches ALL scans, may match noise or wrong peptides |
| Arbitrary scoring weights | Hardcoded `[1.0, 0.7, 0.5, ...]` biases toward monoisotopic |
| Single charge state | No multi-charge combination or consistency check |

---

## 2. What Established Tools Do

### The Centroid Mass Method (Gold Standard)

Every established HDX-MS tool computes deuterium uptake via the centroid
(intensity-weighted average) of the isotopic envelope. This is the method
introduced by Zhang & Smith (1993) and codified in the community consensus
guidelines (Masson et al. 2019).

**Step 1 -- Centroid on the m/z scale:**

```
              sum_i( mz_i * I_i )
C_mz    = ---------------------------
               sum_i( I_i )
```

where `mz_i` and `I_i` are the observed m/z and intensity of each isotopic
peak in the envelope. This uses the OBSERVED positions, not theoretical ones.

**Step 2 -- Convert to neutral mass:**

```
M_centroid = C_mz * z  -  z * m_proton
```

where `m_proton = 1.00728 Da`.

**Step 3 -- Deuterium uptake (Da):**

```
uptake(t) = M_centroid(t)  -  M_centroid(t=0)
```

where `t=0` is the undeuterated control. Each deuteron adds ~1.006 Da, so
uptake in Da approximately equals the number of deuterons incorporated.

**Step 4 -- Back-exchange correction (requires fully deuterated control):**

```
                    m(t) - m(0)
D_corrected(t) = ---------------  *  N
                  m(FD) - m(0)
```

where:
- `m(FD)` = centroid mass of the fully deuterated control
- `N` = maximum exchangeable amide hydrogens = `len(seq) - count(Pro) - 1`

The **relative fractional uptake (RFU)** is the normalized version:

```
              m(t) - m(0)
RFU(t) = -------------------       (range: 0 to 1)
           m(FD) - m(0)
```

### Tool Landscape

| Tool | Input | Computes centroids? | Back-exchange? | Language |
|------|-------|---------------------|----------------|----------|
| DynamX (Waters) | Raw vendor data | Yes | No | Proprietary |
| HDExaminer (Trajan) | Raw vendor data | Yes (pattern matching) | No | Proprietary |
| TheDeuteriumCalculator | mzML + peptide list | Yes (weighted avg mass) | No (recovery rate param) | Python |
| pyHXExpress | mzML / HDExaminer export | Yes (centroid + binomial fit) | No | Python |
| ExMS2 | mzML | Yes (full pipeline) | Yes | MATLAB |
| PyHDX | DynamX/HDExaminer CSV | No (imports pre-computed uptake) | Yes (`apply_control()`) | Python |
| DECA | DynamX/HDX Workbench CSV | No (imports pre-computed uptake) | Yes (global + LEAP) | Python |
| HaDeX | DynamX cluster CSV | No (imports centroid column) | Yes | R |
| HDXBoxeR | HDExaminer CSV | No (imports pre-computed uptake) | Yes | R |

### Key Insight

No established Python tool provides a complete raw-mzML-to-uptake pipeline.
The field is split: vendor software handles raw data extraction, and
open-source tools handle downstream analysis. TheDeuteriumCalculator is the
closest Python reference for centroid calculation from mzML.

---

## 3. Available Reference Implementations

### TheDeuteriumCalculator (OUWuLab)

**The best reference for what our script tries to do.** Single Python script
(~1400 lines) that reads mzML and computes weighted average mass:

```
WeightedAveMass = [SUM(mz_n * I_n / I_total)] * z - z * 1.00627
```

Key design decisions to learn from:
- Uses Accurate Mass and Time (AMT) matching with RT window (+/- 30 sec)
- Gaussian fitting to evaluate isotope distribution quality (R^2 cutoff)
- Sliding window to sum intensities across retention time range
- Binary search for peak matching within PPM tolerance
- Reports `mass_deuterated - mass_undeuterated` as uptake

Limitations: requires pre-identified peptide list, no FD control correction,
single monolithic script.

**Citation:** Espada et al. (2023) *J Proteome Res* 22:585-592,
doi:10.1021/acs.jproteome.2c00558

### PyHDX (Jhsmit/PyHDX)

**Best for post-processing** once we have uptake values.

Key functions to integrate with:
- `correct_d_uptake(peptides, drop_first=1, d_percentage=100)` -- adjusts
  exchangeable residue count, marks prolines and N-terminal as non-exchanging
- `apply_control(experiment, fd_control, nd_control)` -- computes RFU with
  error propagation: `RFU = (u - n) / (f - n)`
- `DeltaGFit` -- PyTorch-based per-residue Gibbs free energy fitting using
  coverage matrix (peptide-to-residue mapping)

Input format (DynamX state data CSV):
```
start, end, sequence, state, exposure, uptake, uptake_sd, maxuptake
```

**Citation:** Smit et al. (2021) *Anal Chem* 93:9497-9504,
doi:10.1021/acs.analchem.1c02155

### HDXrate (Jhsmit/HDXrate)

Calculates intrinsic amide hydrogen exchange rates from amino acid sequence,
temperature, and pH. Implements the Englander group methodology.

```python
from hdxrate import k_int_from_sequence
k_int = k_int_from_sequence(sequence, temperature=293.15, pH_read=7.0)
# Returns per-residue intrinsic exchange rates (s^-1)
```

Needed for protection factor calculation: `PF = k_int / k_obs`.

**Citation:** Based on Bai et al. (1993) and Connelly et al. (1993)

### ExMS2 (MATLAB)

The most complete raw-to-uptake pipeline, but MATLAB-based. Compiled
standalone (no license needed). Reads mzML, identifies peptides, extracts
envelopes, computes uptake. Best algorithmic reference for a full pipeline.

**Source:** https://hx2.med.upenn.edu/

### pyHXExpress (tuttlelm/pyHXExpress)

Useful reference for:
- `hxex_mzml_utils.py` -- mzML reading via pyteomics
- `hxex.py` -- centroid calculation: `centroid = sum(mz * y) / sum(y)`
- `count_amides()` -- exchangeable amide calculation with proline handling
- Binomial fitting for EX1 multimodal detection

**Citation:** Tuttle et al. (2025) *JASMS*, PMC11952558

---

## 4. Compatibility with Our Pipeline

### Current Pipeline Steps

| Step | Script | Status |
|------|--------|--------|
| 1. mzML conversion | (msConvert, external) | Working |
| 2. Peptide discovery | `scripts/hdx/discover_peptides.py` (v3) | Working |
| 3. HDX reconstruction | `reconstruct_hdx_from_raw.py` (v1) | **BROKEN** |
| 4. Kinetics fitting | `fit_hdx_kinetics.py` | Working (but fed garbage) |

### What Step 4 Expects

The kinetics fitting script (`fit_hdx_kinetics.py`) reads a CSV with these
required columns:

| Column | Type | Description |
|--------|------|-------------|
| `fragment_id` | str | e.g. `AVRGSPAIN\|21-29` |
| `time_min` | float | Labeling time in minutes |
| `deuterium_uptake_pct` | float | Uptake as percentage (0-100) |

Or alternatively: `exchange_ratio_vs_max` (0-1), which gets multiplied by 100.

The fitting script groups by `fragment_id`, aggregates duplicate timepoints
by mean, requires at least `--min-points` (default 4) timepoints, and fits
mono-exponential, stretched-exponential, and bi-exponential models with AICc
selection.

**What v2 reconstruction must output:** A CSV with `fragment_id`, `time_min`,
and `deuterium_uptake_pct` (or at minimum `exchange_ratio_vs_max`). Values
must be in a sensible range -- roughly 0-100% for pct or 0-1 for ratio.

### What PyHDX Expects (for downstream residue-level analysis)

If we want to feed results to PyHDX for per-residue deltaG fitting, we need
a CSV in DynamX-like format:

| Column | Description |
|--------|-------------|
| `start` | Residue start position (1-based) |
| `end` (or `stop`) | Residue end position |
| `sequence` | Peptide amino acid sequence |
| `state` | Protein state label (e.g. "apo", "holo") |
| `exposure` | Labeling time in seconds |
| `uptake` | Absolute deuterium uptake (Da) |
| `uptake_sd` | Standard deviation of uptake |
| `maxuptake` | Maximum exchangeable amides (N) |

Optional for back-exchange correction: `fd_uptake`, `fd_uptake_sd`.

PyHDX uses half-open intervals: `stop = end + 1`.

---

## 5. Proposed Versions

### v1 (current, broken): `reconstruct_hdx_from_raw.py`

**What it does:**
- Naive isotope peak matching with hardcoded weights
- Computes monoisotopic mass estimate (subtracts `iso_idx * NEUTRON`)
- Single best scan per peptide per run (no XIC)
- No RT constraint (searches all scans)
- No back-exchange correction
- Normalizes to longest timepoint (not FD control)

**Result:** exchange ratios in [-36, +11] range -- scientifically meaningless.

**Disposition:** Superseded. Keep for reference, do not use.

---

### v2 (centroid fix): Fix the critical bug

**Goal:** Minimal change, maximum impact. Get exchange ratios into [0, 1].

**Changes from v1:**

1. **Replace monoisotopic mass with centroid mass.** In `best_feature_for_peptide()`,
   change the mass calculation:

   ```python
   # BEFORE (v1 -- wrong):
   neutral_series = (matched_mz[used] * z - z * PROTON) - iso_idx[used] * NEUTRON
   neutral_mass_est = float(np.average(neutral_series, weights=matched_int[used]))

   # AFTER (v2 -- correct):
   mass_series = matched_mz[used] * z - z * PROTON
   centroid_mass = float(np.average(mass_series, weights=matched_int[used]))
   ```

   No `iso_idx` subtraction. The centroid mass now includes the deuterium shift.

2. **Compute uptake as centroid difference.** In `finalize_metrics()`:

   ```python
   uptake_da = centroid_mass(t) - centroid_mass(mapping)
   ```

   This replaces `delta_mass_vs_baseline` with a physically meaningful quantity.

3. **Compute exchangeable amides (N).** For each peptide:

   ```python
   n_prolines = sequence[1:].count('P')   # exclude first residue
   N = len(sequence) - n_prolines - 1     # subtract N-terminal
   ```

4. **Compute fractional uptake without FD control (theoretical max).**

   ```python
   fractional_uptake = uptake_da / (N * 1.00628)
   ```

   This assumes zero back-exchange (underestimates true fractional uptake by
   10-40%), but gives values in approximately [0, 1] which is far better than
   the current [-36, +11].

5. **Report both absolute uptake (Da) and fractional uptake.**

**What stays the same:**
- Same isotope matching approach (hardcoded peak positions)
- Same scoring (weighted intensity sum)
- Same single-scan strategy (best scan per peptide per run)
- Same CLI interface and output structure
- No RT constraint yet
- No FD control correction yet

**Expected outcomes:**
- Exchange ratios in approximately [0, 1] range
- Known peptides show increasing uptake with time (monotonic for EX2)
- Correlation with HDExaminer reference values from Result.xlsx
- Kinetics fitting produces sensible half-lives

**Implementation effort:** ~30 minutes. Change ~10 lines of code.

**Test:** Run on 11 known peptides, compare uptake vs HDExaminer values.

---

### v3 (ms_deisotope envelopes): Proper envelope detection

**Goal:** Replace naive isotope matching with ms_deisotope deconvolution,
matching the approach used in peptide discovery v2/v3.

**Changes from v2:**

1. **Use ms_deisotope for per-scan deconvolution.** Replace the manual
   `match_targets()` + `best_feature_for_peptide()` with:

   ```python
   from ms_deisotope import deconvolute_peaks
   result = deconvolute_peaks(
       peaklist,
       charge_range=(charge_min, charge_max),
       error_tolerance=ppm * 1e-6,
       left_search_limit=3,
       truncate_after=0.95,
       iterations=10,
   )
   ```

2. **Match deconvoluted envelopes to target peptides** by neutral mass
   (same approach as `discover_peptides.py` Phase 2).

3. **Compute centroid from the deconvoluted envelope.** Use the observed
   m/z values from the ms_deisotope-identified envelope, not theoretical
   positions:

   ```python
   env_mz = np.array([p.mz for p in envelope.experimental])
   env_int = np.array([p.intensity for p in envelope.experimental])
   centroid_mz = np.average(env_mz, weights=env_int)
   centroid_mass = centroid_mz * z - z * PROTON
   ```

4. **Add cosine similarity quality score.** Score the envelope against
   the averagine theoretical distribution (same as discovery v3):

   ```python
   from ms_deisotope.averagine import peptide as averagine
   theo = averagine.isotopic_cluster(neutral_mass, charge=z)
   cosine = dot(theo_int, obs_int) / (norm(theo_int) * norm(obs_int))
   ```

   Reject envelopes with cosine < 0.7.

5. **Better charge state handling.** ms_deisotope resolves charge states
   properly via averagine pattern fitting (PenalizedMSDeconV). If the same
   peptide is detected at multiple charge states, compute centroid mass from
   each and combine using intensity-weighted average:

   ```python
   M_combined = sum(I_z * M_z) / sum(I_z)
   ```

**What stays the same:**
- Still single best scan per run
- Still no RT constraint
- Still no XIC integration
- Same output format as v2

**Expected improvements over v2:**
- More accurate centroid (envelope boundaries determined by averagine model)
- Better charge state assignment (pattern-based, not just highest intensity)
- Quality filtering via cosine similarity
- Handles overlapping envelopes via peak dependency graph

**Dependencies added:** `ms_deisotope >= 0.0.60` (already required by
peptide discovery).

**Implementation effort:** ~2 hours. Reuse patterns from `discover_peptides.py`.

**Test:** Same 11 known peptides. Expect tighter correlation with reference
values than v2.

---

### v4 (XIC integration): Full pipeline

**Goal:** Extract centroid from the full chromatographic peak, not just one
scan. Add RT constraints. Add optional back-exchange correction. Produce
output compatible with both kinetics fitting and PyHDX.

**Changes from v3:**

1. **Add RT constraint from peptide discovery.** Use the `rt_apex` and
   `rt_start`/`rt_end` columns from the v3 discovery output (or from a
   reference peptide list) to limit the scan search window:

   ```python
   rt_window = (rt_apex - 0.5, rt_apex + 0.5)  # +/- 30 sec
   ```

   For deuterated runs, allow a slightly wider window to account for RT shift
   from deuteration.

2. **Deconvolute all scans in the RT window.** Instead of searching all scans
   and keeping the best, process only scans within the RT window and keep ALL
   detections.

3. **Build XIC features.** Group consecutive detections into XIC features
   (same approach as discovery v3 Phase 3). Require minimum 2 consecutive
   scans.

4. **Compute centroid from averaged envelope.** For each XIC feature, compute
   the centroid from the summed (or averaged) isotope envelope across all scans
   in the elution window:

   ```python
   # Sum envelopes across scans (aligned by isotope index)
   summed_int = np.zeros(n_peaks)
   for scan in feature_scans:
       summed_int += scan.envelope_intensities
   centroid_mz = np.average(env_mz, weights=summed_int)
   ```

   This reduces noise from individual scans and gives a more robust centroid.

5. **Cosine similarity filter per scan.** Reject individual scan detections
   where the envelope shape disagrees with theory. Only include scans with
   cosine > 0.7 in the averaged envelope.

6. **Optional back-exchange correction.** If a fully deuterated (FD) control
   run is provided:

   ```python
   # In the FD run, compute centroid for each peptide
   centroid_fd = compute_centroid(fd_run, peptide)

   # Corrected fractional uptake
   RFU = (centroid_t - centroid_0) / (centroid_fd - centroid_0)
   D_corrected = RFU * N
   ```

   Propagate uncertainty using the PyHDX formula:

   ```python
   sigma_RFU = sqrt(
       (1/(f-n))**2 * sigma_u**2 +
       ((u-f)/(f-n)**2)**2 * sigma_n**2 +
       ((n-u)/(f-n)**2)**2 * sigma_f**2
   )
   ```

7. **D2O percentage correction.** If the labeling buffer is not 100% D2O,
   scale the theoretical maximum:

   ```python
   N_effective = N * d2o_fraction
   ```

8. **Dual output format.** Produce two output CSVs:

   **a. Kinetics-compatible output** (for `fit_hdx_kinetics.py`):

   | Column | Description |
   |--------|-------------|
   | `fragment_id` | `SEQUENCE\|start-end` |
   | `time_min` | Labeling time in minutes |
   | `deuterium_uptake_da` | Absolute uptake (Da) |
   | `deuterium_uptake_pct` | Uptake as % of max exchangeable |
   | `exchange_ratio_vs_max` | RFU (0-1), with FD correction if available |
   | `centroid_mass` | Observed centroid mass |
   | `n_scans` | Number of scans in XIC feature |
   | `cosine_sim` | Mean cosine similarity across scans |

   **b. PyHDX-compatible output** (for downstream residue-level analysis):

   | Column | Description |
   |--------|-------------|
   | `start` | 1-based residue start |
   | `end` | Residue end (inclusive) |
   | `stop` | `end + 1` (half-open, PyHDX convention) |
   | `sequence` | Peptide sequence |
   | `state` | Protein state label |
   | `exposure` | Labeling time in seconds |
   | `uptake` | Absolute uptake (Da) |
   | `uptake_sd` | Standard deviation (from replicate scans or runs) |
   | `maxuptake` | N (exchangeable amides) |
   | `fd_uptake` | FD control uptake (if available) |
   | `fd_uptake_sd` | FD control uptake SD (if available) |

**New CLI arguments:**

```
--fd-pattern       Regex for fully deuterated control run (e.g. "24h_FD")
--rt-source        CSV or discovery output with RT information
--rt-tolerance     RT window half-width in minutes (default 0.5)
--d2o-fraction     D2O fraction in labeling buffer (default 0.9)
--min-cosine       Minimum cosine similarity (default 0.7)
--min-scans        Minimum scans per XIC feature (default 2)
--output-pyhdx     Also write PyHDX-compatible CSV (flag)
```

**Implementation effort:** ~1 day. Largest change is the XIC integration logic.

**Test:** Run on 11 known peptides and top 200 v3-discovered peptides.
Compare against HDExaminer reference values.

---

## 6. Test Plan

### Ground Truth

- **11 known peptides** from `Result.xlsx` with HDExaminer uptake values
- **6 timepoints** per peptide (10 sec to 24 h)
- **Reference exchange ratios** from PEAKS/HDExaminer (the gold standard for
  this dataset)

### Per-Version Test Protocol

For each version (v2, v3, v4):

1. **Run on 11 known peptides** across all timepoints.

2. **Compare exchange ratios** against HDExaminer reference:
   - Scatter plot: our uptake (Da) vs HDExaminer uptake (Da)
   - Pearson correlation (expect r > 0.9 for v2, r > 0.95 for v3/v4)
   - RMSE in Da (expect < 0.5 Da for v2, < 0.3 Da for v3/v4)
   - % of values within 10% of reference (expect > 70% for v2, > 85% for v4)

3. **Check monotonicity:** For each peptide under EX2 kinetics, uptake should
   increase with time (or plateau). Flag non-monotonic peptides.

4. **Kinetics validation:** Feed output to `fit_hdx_kinetics.py`:
   - Do models converge? (v1 mostly fails to converge meaningfully)
   - Are R^2 values > 0.8 for known peptides?
   - Do half-lives fall in physiologically reasonable ranges?

5. **Broader test on v3-discovered peptides:** Run on top 200 peptides from
   peptide discovery v3 output. Check:
   - Distribution of exchange ratios (should cluster in [0, 1])
   - No peptides with uptake > N (impossible; indicates a bug)
   - Coverage heatmap makes biological sense (exposed loops exchange fast,
     beta sheets exchange slow)

### Test Data Locations

| File | Path | Description |
|------|------|-------------|
| Known peptides | `Result.xlsx` (or extracted CSV) | 11 peptides, HDExaminer values |
| Discovery output | `output/03_hdx_peptide_discovery_v3/` | v3-discovered peptides |
| mzML files | (raw dir, converted) | All timepoint runs |
| v1 (broken) output | `output/03_hdx_reconstruction_old/` | Baseline for comparison |

### Success Criteria

| Metric | v2 target | v3 target | v4 target |
|--------|-----------|-----------|-----------|
| Exchange ratio range | [~0, ~1] | [~0, ~1] | [0, 1] (corrected) |
| Pearson r vs reference | > 0.85 | > 0.92 | > 0.95 |
| RMSE (Da) | < 0.8 | < 0.5 | < 0.3 |
| Kinetics R^2 (median) | > 0.7 | > 0.8 | > 0.9 |
| Peptides with uptake > N | 0 | 0 | 0 |

---

## 7. Integration with Downstream Tools

### Kinetics Fitting (`fit_hdx_kinetics.py`)

**Required columns:** `fragment_id`, `time_min`, `deuterium_uptake_pct`
(or `exchange_ratio_vs_max`).

The fitting script applies three kinetic models:

- `mono_exp`: `U_inf * (1 - exp(-k*t))`
- `stretched_exp`: `U_inf * (1 - exp(-(k*t)^beta))`
- `bi_exp`: `U_inf * (w*(1-exp(-k_fast*t)) + (1-w)*(1-exp(-k_slow*t)))`

Selects by AICc. Reports half-life, k_app, and empirical early/late slopes.

**What reconstruction must provide:** Sensible uptake values (0-100% range)
at each timepoint. The fitting script handles missing data and duplicate
averaging internally.

### PyHDX (Residue-Level Analysis)

To use PyHDX downstream for per-residue Gibbs free energy fitting:

1. Export a CSV with DynamX-compatible columns (see v4 output format above).
2. PyHDX reads this via `read_dynamx()` or `hdxms-datasets` format registry.
3. Apply `correct_d_uptake(drop_first=1)` to mark prolines and N-terminal.
4. If FD control is available, apply `apply_control()` for RFU.
5. PyHDX's `DeltaGFit` fits per-residue deltaG using the peptide-residue
   coverage matrix.

**Key requirement:** `uptake` must be in Da (absolute mass shift), not
percentage. PyHDX expects absolute uptake and computes fractional internally.

### HDXrate (Intrinsic Rates)

If we want to compute protection factors:

```python
from hdxrate import k_int_from_sequence
k_int = k_int_from_sequence(
    sequence,
    temperature=293.15,  # K (20 C)
    pH_read=7.4,
    d_percentage=90.0,
)
# PF = k_int / k_obs (from kinetics fit)
```

This requires the kinetics fitting to produce reliable `k_obs` values,
which in turn requires reliable uptake values from reconstruction.

---

## 8. References

### Foundational HDX-MS Methods

| Paper | DOI | Relevance |
|-------|-----|-----------|
| Zhang & Smith (1993) *Protein Sci* 2:522-531 | 10.1002/pro.5560020404 | Foundational HDX-MS method, back-exchange correction formula |
| Engen & Smith (2001) *Anal Chem* 73:256A-265A | -- | Established centroid method as HDX standard |
| Engen (2009) *Anal Chem* 81:7870-7875 | 10.1021/ac802405a | Review of HDX-MS analysis methods |
| Konermann et al. (2011) *Mass Spec Rev* 30:268-298 | 10.1002/mas.20290 | Comprehensive review of HDX-MS theory |

### Community Guidelines

| Paper | DOI | Relevance |
|-------|-----|-----------|
| Masson et al. (2019) *Nature Methods* 16:595-602 | 10.1038/s41592-019-0459-y | Community consensus guidelines for HDX-MS |
| Hageman & Bhatt (2024) *Chem Rev* | 10.1021/acs.chemrev.4c00438 | Comprehensive review of HDX-MS computational tools |

### Practical HDX-MS Analysis

| Paper | DOI | Relevance |
|-------|-----|-----------|
| Houde et al. (2011) *Anal Chem* 83:5789-5794 | 10.1021/ac200258g | Multi-charge-state handling, biopharm HDX |
| Weis et al. (2006) *JASMS* 17:1700-1703 | 10.1016/j.jasms.2006.07.025 | Practical HDX data analysis considerations |
| Bai et al. (1993) *Proteins* 17:75-86 | -- | Intrinsic exchange rates (basis for HDXrate) |
| Senko et al. (1995) *JASMS* 6:229-233 | -- | Averagine model for isotope distributions |

### Software Papers

| Tool | Paper | DOI |
|------|-------|-----|
| TheDeuteriumCalculator | Espada et al. (2023) *J Proteome Res* 22:585-592 | 10.1021/acs.jproteome.2c00558 |
| PyHDX | Smit et al. (2021) *Anal Chem* 93:9497-9504 | 10.1021/acs.analchem.1c02155 |
| pyHXExpress | Tuttle et al. (2025) *JASMS* | PMC11952558 |
| HaDeX | Puchala et al. (2020) *Bioinformatics* 36:4516-4518 | 10.1093/bioinformatics/btaa587 |
| HDXBoxeR | Janowska et al. (2024) *Bioinformatics* 40:btae479 | 10.1093/bioinformatics/btae479 |
| DECA | Lumpkin & Komives (2019) *Mol Cell Proteomics* 18:2516-2523 | 10.1074/mcp.TIR119.001731 |

### GitHub Repositories

| Tool | Repository | Language |
|------|------------|----------|
| TheDeuteriumCalculator | github.com/OUWuLab/TheDeuteriumCalculator | Python |
| PyHDX | github.com/Jhsmit/PyHDX | Python |
| HDXrate | github.com/Jhsmit/HDXrate | Python |
| pyHXExpress | github.com/tuttlelm/pyHXExpress | Python |
| DECA | github.com/komiveslab/DECA | Python |
| HaDeX | github.com/hadexversum/HaDeX | R |
| HDXBoxeR | github.com/mkajano/HDXBoxeR | R |
| ExMS2 | hx2.med.upenn.edu | MATLAB |
| hdxms-datasets | github.com/Jhsmit/hdxms-datasets | Python |
