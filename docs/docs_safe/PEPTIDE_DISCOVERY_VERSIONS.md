# Peptide Discovery: Version Comparison and Algorithms

## Results Comparison

| Metric | v1 (naive) | v2 (ms_deisotope) | v3 (XIC + cosine + M-1) |
|---|---|---|---|
| **Peptides found** | 2687 | 1570 | 1068 |
| **Known peptides found** | 11/11 (100%) | 11/11 (100%) | 11/11 (100%) |
| **Sequence coverage** | 100% | 100% | 100% |
| **Estimated false positives** | ~2676 | ~1559 | ~1057 |
| **Scoring method** | Weighted intensity sum | PenalizedMSDeconV | Integrated × cosine × M-1 penalty |
| **Cosine similarity** | Not computed | Not computed | median 0.989, min 0.826 |
| **Mass error (ppm)** | Not reported | -0.61 ± 7.17 | -0.54 ± 6.38 |
| **Scans per peptide** | 1 (single best) | 1 (single best) | median 10, max 1232 |
| **Charge state method** | Best intensity across z=1-6 | ms_deisotope averagine fit | ms_deisotope averagine fit |
| **M-1/M-2 check** | None | None (but left-search) | Explicit penalty |
| **Shape validation** | None | G-test (within ms_deisotope) | G-test + cosine similarity |
| **Runtime (8 workers)** | ~5 min | ~15 min (sequential) | ~5 min (parallel) |
| **Script** | discover_peptides_v1 (deleted) | discover_peptides_v2_legacy.py | discover_peptides.py |
| **Output** | output/03_hdx_peptide_discovery/ | output/03_hdx_peptide_discovery_v2_legacy/ | output/03_hdx_peptide_discovery_v3/ |

---

## v1: Naive MS1 Isotope Matching

**Script:** (original, superseded)

### Algorithm

#### Step 1: In-silico pepsin digest
Generate every possible sub-peptide of the protein sequence.

```
For each start position (0 to N-1):
    For each end position (start+min_len to start+max_len):
        peptide = sequence[start:end]
        neutral_mass = pyteomics.mass.fast_mass(peptide)
```

For TTR (129 residues, min_len=5, max_len=40): ~3,870 candidate peptides.

#### Step 2: Compute theoretical m/z
For each peptide at each charge state z=1 to z=6:

```
mz_monoisotopic = (neutral_mass + z × 1.007276) / z
```

Keep only targets where 300 < m/z < 2000.

#### Step 3: Single-scan search
For each MS1 scan, for each target:

1. Compute tolerance: `tol = mz × ppm / 1e6`
2. Check 5 isotope peaks (M+0 through M+4):
   - `mz_M+i = mz_mono + i × 1.003355 / z`
   - Find strongest peak within ±tol
3. Require ≥3 of 5 peaks matched
4. Score: `Σ(weight_i × intensity_i)` where weights = [1.0, 0.8, 0.6, 0.4, 0.2]
5. Keep best score per (sequence, charge) across all scans

#### Step 4: Deduplicate
Keep best charge state per sequence (highest score).

### Limitations
- No isotope shape validation (M+3 could be 10× too high and still pass)
- No M-1/M-2 check (overlapping ions create false matches)
- No charge state resolution (z=2 vs z=4 distinguished only by intensity, not pattern)
- Single scan (no chromatographic integration)
- Score is raw intensity (no statistical meaning)

---

## v2: ms_deisotope Per-Scan Deconvolution

**Script:** `discover_peptides_v2_legacy.py`

### Algorithm

#### Step 1: In-silico digest
Same as v1.

#### Step 2: Per-scan deconvolution with ms_deisotope
For each MS1 scan:

```python
result = deconvolute_peaks(
    peaklist,                    # [(mz, intensity), ...]
    charge_range=(1, 6),         # try all charge states
    error_tolerance=15e-6,       # 15 ppm
    left_search_limit=3,         # look 3 Da left for true monoisotopic
    truncate_after=0.95,         # keep isotopes covering 95% of signal
    iterations=10,               # iterative subtraction rounds
)
```

**What ms_deisotope does internally:**

1. **Peak indexing:** Sort all peaks by m/z for fast lookup.

2. **For each peak, try all charge states z=1 to z=6:**
   - Compute expected neutral mass: `M = mz × z - z × 1.007276`
   - Generate theoretical isotope distribution from averagine model:
     - Interpolate elemental composition (C, H, N, O, S) from neutral mass using averagine coefficients (Senko et al. 1995)
     - Compute exact isotope distribution using brainpy (BRAIN algorithm)
     - Truncate at 95% cumulative intensity
   - Match observed peaks to theoretical positions (within ppm tolerance)
   - Score using PenalizedMSDeconVFitter

3. **PenalizedMSDeconVFitter scoring:**
   ```
   Per-peak score = mass_accuracy × abundance_agreement × sqrt(theoretical_intensity)

   where:
     mass_accuracy = 1 - |mz_observed - mz_theoretical| / tolerance
     abundance_agreement = 1 - |intensity_observed/total - intensity_theoretical/total|

   MSDeconV_score = Σ per_peak_scores

   G_test = 2 × Σ (observed_i × ln(observed_i / expected_i))

   Final_score = MSDeconV_score × (1 - |G_test|)
   ```
   The G-test penalizes envelopes whose intensity distribution doesn't match the theoretical shape.

4. **Left-search for monoisotopic:** For each candidate envelope, also try shifting the monoisotopic assignment 1-3 peaks to the left. The correct assignment produces the best score.

5. **Peak dependency graph:** When multiple envelope candidates share peaks, build a graph where peaks are nodes and envelopes are hyperedges. Solve by selecting the highest-scoring non-conflicting set.

6. **Iterative subtraction:** After assigning the best envelopes, subtract their peaks from the spectrum and repeat (up to 10 iterations) to find weaker overlapping envelopes.

#### Step 3: Match to candidates
For each deconvoluted envelope with score ≥ 10:
- Compare its neutral mass to all candidate peptide masses
- Match if within ±0.02 Da AND ±15 ppm
- Keep best score per sequence across all scans

### Improvements over v1
- Proper charge state resolution (averagine pattern differs by charge)
- G-test shape penalty (rejects bad isotope distributions)
- Left-search (finds correct monoisotopic even when M+0 is weak)
- Peak dependency graph (resolves overlapping envelopes)
- Meaningful score (not raw intensity)

### Remaining limitations
- Still single-scan (best scan only, no chromatographic integration)
- No explicit M-1/M-2 check (left-search partially addresses this)
- No cosine similarity threshold on final output
- Sequential processing (slow: ~15 min for 3587 scans)

---

## v3: XIC Integration + Cosine + M-1 Penalty + Parallel

**Script:** `discover_peptides.py`

### Algorithm

#### Phase 1a: Read all MS1 scans
Read every MS1 scan from the mzML file, storing m/z array, intensity array, retention time, and scan number.

```
For each MS1 scan in mzML:
    Store: (scan_number, RT, mz_array, intensity_array)
```

#### Phase 1b: Parallel deconvolution
Distribute scans across 8 worker processes. Each worker runs ms_deisotope independently:

```python
with ProcessPoolExecutor(max_workers=8) as executor:
    results = executor.map(_decon_one_scan, work_items, chunksize=50)
```

Each worker returns serializable envelope dicts (not ms_deisotope objects, which can't be pickled across processes):

```python
{"score": 245.3, "charge": 2, "mz": 711.324, "nm": 1420.633,
 "n_peaks": 5, "env": [{"mz": 711.324, "int": 24506.0}, ...]}
```

The `chunksize=50` batches scans to reduce IPC overhead.

#### Phase 2: Match envelopes to candidates
For each deconvoluted envelope across all scans:

1. Look up its neutral mass in a sorted mass index (binary search)
2. If within ±0.02 Da AND ±15 ppm of a candidate peptide, record as a hit
3. Key: `(sequence, charge_state)`
4. Store: scan number, RT, score, m/z, ppm error, envelope peaks, raw spectrum

#### Phase 3: Build XIC features (multi-scan grouping)
For each (sequence, charge) group:

1. Sort hits by scan number
2. Split into consecutive features by RT gap:
   ```
   If RT_current - RT_previous > 0.5 min:
       Start new feature
   Else:
       Add to current feature
   ```
3. Discard features with fewer than `min_scans` (default 2) consecutive hits

This ensures each "discovered" peptide is backed by detection in multiple consecutive scans — not just a single lucky match.

#### Phase 4: Score each feature

For each XIC feature (a peptide detected across N consecutive scans):

**4a. Integrated score:**
```
score_integrated = Σ (ms_deisotope_score for each scan in the feature)
```
Sums across all scans in the elution window. A peptide detected in 10 scans with average score 200 gets integrated score 2000.

**4b. Cosine similarity against averagine:**
Using the best-scoring scan's envelope:

```python
# Theoretical distribution from averagine + brainpy
theo = averagine.isotopic_cluster(neutral_mass, charge=1, truncate_after=0.999)
theo_intensities = normalize_to_max(theo)

# Observed distribution
obs_intensities = normalize_to_max(best_envelope)

# Cosine similarity
cosine = dot(theo, obs) / (norm(theo) × norm(obs))
```

Values:
- 1.0 = perfect match to theoretical
- > 0.95 = excellent
- 0.85-0.95 = acceptable
- < 0.7 = likely false positive or interference

This catches Problem 4 (abnormal M+3 height, distorted pattern).

**4c. M-1/M-2 penalty:**
Check the raw spectrum at the best scan for peaks below the monoisotopic:

```python
For offset in [-1, -2]:
    mz_check = mz_mono + offset × 1.003355 / charge
    Find strongest peak within ±15 ppm of mz_check
    If found:
        ratio = intensity_M_minus / intensity_M0
        If ratio > 0.3:  # significant
            penalty *= (1.0 - ratio)
```

Values:
- 1.0 = clean (no M-1/M-2 peaks)
- 0.5 = moderate (M-1 at 50% of M+0 → penalty 0.5)
- 0.0 = severe (M-1 stronger than M+0 → penalty 0.0, peptide rejected)

This catches Problem 2 (false monoisotopic from overlapping heavier ion).

**4d. Final score:**
```
score_final = score_integrated × cosine_similarity × m_minus_penalty
```

This multiplicative combination means a peptide needs ALL THREE to score well:
- Detected in multiple scans (high integrated score)
- Isotope pattern matches theory (high cosine)
- No interfering M-1/M-2 peaks (penalty near 1.0)

#### Phase 5: Filter and deduplicate

1. Keep best charge state per sequence (highest score_final)
2. Filter: `score_final >= 50` AND `cosine_sim >= 0.7`
3. Sort by protein position

### Output columns

| Column | Description |
|---|---|
| `sequence` | Amino acid sequence |
| `start`, `end` | 1-based protein position |
| `charge` | Best-scoring charge state |
| `mz_mono` | Observed monoisotopic m/z |
| `neutral_mass` | Observed neutral mass |
| `theo_mass` | Theoretical neutral mass from sequence |
| `mass_error_ppm` | Average ppm error across all scans in the feature |
| `score_integrated` | Sum of per-scan ms_deisotope scores |
| `score_best_scan` | Highest single-scan score |
| `score_final` | integrated × cosine × M-1 penalty |
| `cosine_sim` | Cosine similarity vs averagine (0-1) |
| `m_minus_penalty` | M-1/M-2 penalty factor (0-1) |
| `n_scans` | Number of consecutive scans with detection |
| `n_peaks` | Isotope peaks in envelope at best scan |
| `rt_start` | Retention time of first detection (min) |
| `rt_apex` | RT of highest-scoring scan (min) |
| `rt_end` | RT of last detection (min) |
| `scan_apex` | Scan number of highest-scoring detection |

### Tunable parameters

| Parameter | Default | Effect |
|---|---|---|
| `--ppm` | 15.0 | Mass tolerance for matching. Lower = stricter, may miss poorly calibrated data |
| `--min-len` | 5 | Shortest peptide. Shorter = more candidates, more noise |
| `--max-len` | 40 | Longest peptide. Longer = more candidates, slower |
| `--charge-min` | 1 | Minimum charge state |
| `--charge-max` | 6 | Maximum charge state. Higher = slower deconvolution |
| `--min-score` | 50 | Minimum final score. Higher = fewer peptides, fewer false positives |
| `--min-cosine` | 0.7 | Minimum cosine similarity. 0.85 = stringent, 0.7 = permissive |
| `--min-scans` | 2 | Minimum consecutive scans. Higher = more confident but may miss low-abundance |
| `-j` / `--workers` | 4 | Parallel deconvolution workers (max 8) |

### Dependencies

| Library | Version | Purpose |
|---|---|---|
| ms_deisotope | 0.0.60 | Per-scan deconvolution (PenalizedMSDeconVFitter) |
| pyteomics | 4.7.5 | Neutral mass calculation, mzML reading |
| numpy | any | Array operations |
| pandas | any | DataFrame output |

### Physical constants

| Constant | Value | Source |
|---|---|---|
| PROTON | 1.007276466812 Da | CODATA 2018 |
| NEUTRON | 1.0033548378 Da | ¹³C - ¹²C mass difference |

---

## What is NOT done (future improvements)

1. **Full 3D feature detection** — currently we group by RT gap, but don't fit Gaussian elution profiles or use OpenMS FeatureFinder
2. **FDR estimation** — no decoy search. Would need scrambled sequences + target-decoy competition
3. **MS/MS confirmation** — no fragment ion (b/y) matching. Would need Sage/Comet when sagepy FASTA bug is fixed
4. **Back-exchange correction** — not relevant for the undeuterated mapping run, but needed for HDX timepoints
5. **Retention time prediction** — could filter by predicted RT (SSRCalc, Prosit) to reject peptides eluting at wrong times
