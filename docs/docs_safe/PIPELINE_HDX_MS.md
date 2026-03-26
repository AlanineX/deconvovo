# HDX-MS Pipeline

Hydrogen-deuterium exchange mass spectrometry analysis for Thermo Orbitrap data.

**Status:** Development (scripts in `scripts/hdx/` and `skills/`)

## Overview

```
Thermo .raw → mzML → Peptide Discovery → HDX Reconstruction → Kinetics Fitting
```

## Data: TTR (Transthyretin)

- **Protein:** Human transthyretin, 129 residues
- **Instrument:** Thermo Orbitrap
- **Digestion:** Online pepsin (non-specific cleavage)
- **Timepoints:** 10 min, 30 min, 3 h, 6 h, 12 h, 24 h + undeuterated mapping run
- **Known peptides:** 11 (from PEAKS/HDExaminer Result.xlsx)

---

## Step 1: RAW to mzML Conversion

Converts Thermo proprietary .RAW files to open mzML format.

**Script:** `skills/public/ms-install-transform/scripts/raw_to_mzml_parallel.py`

**Tool:** ThermoRawFileParser (Mono/.NET, runs on Linux)

**Input:** Directory of .raw files

**Output:** One .mzML file per .raw (77-94 MB each for TTR dataset)

**CLI:**
```bash
python skills/public/ms-install-transform/scripts/raw_to_mzml_parallel.py \
    --input-dir data_1_thermo/merged_20260217_TTR_TimeP/raw \
    --out-dir output/01_hdx_mzml \
    --jobs 8
```

**Parameters:**
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--thermo-parser` | `tools/ThermoRawFileParser-linux/ThermoRawFileParser` | Path to converter binary |
| `--format` | 1 | Output format (1 = mzML) |
| `--jobs` | 1 | Parallel workers (up to 8) |
| `--overwrite` | off | Reprocess existing outputs |

**Performance:** 81s sequential → 16s parallel (8 workers), 5x speedup

**How Thermo .RAW stores data:**
- Proprietary binary format, no public spec
- Scans are sequential: MS1 and MS2 interleaved (DDA mode)
- Each scan has a filter string encoding: analyzer type, polarity, MS level, mass range, fragmentation
- Retention time is stored per scan
- XICs extracted by querying m/z ± tolerance across all scans
- No native Linux reader — must convert to mzML first

---

## Step 2: Peptide Discovery

Find all detectable pepsin-digest peptides from the undeuterated mapping run.

**Script:** `scripts/hdx/discover_peptides.py`

**Approach:** In-silico pepsin digest + MS1 isotopic envelope matching

**Input:**
- Mapping run mzML
- Protein sequence (TTR, 129 residues)
- Optional: known peptide list for fidelity comparison

**Output:**
- `discovered_peptides.csv` — all matched peptides with scores
- `peptide_comparison.csv` — known vs discovered comparison
- `residue_coverage.csv` — per-residue coverage depth

**CLI:**
```bash
python scripts/hdx/discover_peptides.py \
    --mzml output/01_hdx_mzml/01_20260217_TTR_peptide_mapping.mzML \
    --known-peptides data_1_thermo/peptide_targets_sequence_start_end.csv \
    --out-dir output/03_hdx_peptide_discovery
```

**Parameters:**
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--ppm` | 15.0 | Mass tolerance for isotope matching |
| `--min-len` | 4 | Minimum peptide length |
| `--max-len` | 40 | Maximum peptide length |
| `--charge-min` | 1 | Minimum charge state to search |
| `--charge-max` | 6 | Maximum charge state to search |

**How it works:**
1. Generate ALL possible sub-peptides of length 4-40 from the protein sequence (non-specific cleavage simulating pepsin)
2. For each peptide × charge state, compute theoretical monoisotopic m/z: `(mass + z × 1.00728) / z`
3. Scan every MS1 spectrum for matching isotopic envelopes (5 peaks, weighted by position)
4. Score: `Σ(weight_i × matched_intensity_i)` where weights = [1.0, 0.8, 0.6, 0.4, 0.2]
5. Require ≥ 3 matched isotope peaks within ± ppm tolerance
6. Keep best score per unique sequence across all scans
7. Compare against known target list for fidelity assessment

**Results for TTR:**
- 2687 unique peptides discovered
- 11/11 known targets found (100% fidelity)
- 100% sequence coverage

---

## Step 3: HDX Reconstruction

Extract deuterium uptake from MS1 isotopic envelopes across all timepoints.

**Script:** `skills/public/hdx-process-extract/scripts/reconstruct_hdx_from_raw.py`

**Input:**
- mzML files (all timepoints)
- Peptide target list (Sequence, Start, End)

**Output:**
- `run_manifest.csv` — all runs with time labels
- `peptide_targets.csv` — input peptides with theoretical masses
- `reconstructed_peptide_timecourse.csv` — long format (one row per peptide per timepoint)
- `reconstructed_peptide_summary.csv` — summary per peptide

**CLI:**
```bash
python skills/public/hdx-process-extract/scripts/reconstruct_hdx_from_raw.py \
    --peptide-source data_1_thermo/peptide_targets_sequence_start_end.csv \
    --mzml-dir output/01_hdx_mzml \
    --out-dir output/02_hdx_reconstruction \
    --jobs 4
```

**Parameters:**
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--charge-min` | 1 | Minimum charge state |
| `--charge-max` | 8 | Maximum charge state |
| `--n-isotopes` | 6 | Isotope peaks to match (M+0 through M+5) |
| `--min-isotopes` | 2 | Minimum matched peaks to accept |
| `--ppm` | 15.0 | Mass tolerance |
| `--jobs` | 1 | Parallel workers (up to 8) |
| `--baseline-pattern` | "peptide_mapping" | Regex to identify undeuterated run |
| `--maxd-pattern` | "24h" | Regex to identify max-exchange run |

**How it works:**
1. **Time parsing:** Extract time labels from filenames (e.g., "10min", "3h", "24h")
2. **Mass calculation:** Neutral mass from sequence via pyteomics
3. **Envelope matching:** For each peptide × charge, scan all MS1 spectra. Score isotopic envelopes by weighted intensity: `weights = [1.0, 0.7, 0.5, 0.35, 0.25, 0.18]`
4. **Best hit selection:** Per peptide × charge × run, keep the scan with highest envelope score
5. **Neutral mass estimation:** Weighted average of `(m/z × z - z × PROTON - iso × NEUTRON)` across matched peaks
6. **Delta mass:** `Δm = mass(timepoint) - mass(baseline)`
7. **Exchange ratio:** `ratio = Δm(timepoint) / Δm(max_exchange)`

**Physical constants:**
| Constant | Value |
|----------|-------|
| PROTON | 1.007276466812 Da |
| NEUTRON | 1.0033548378 Da |

---

## Step 4: Kinetics Fitting

Fit exchange kinetics to deuterium uptake curves.

**Script:** `skills/public/hdx-kinetics-fit-plot/scripts/fit_hdx_kinetics.py`

**Input:** `reconstructed_peptide_timecourse.csv`

**Output:**
- `kinetics_fit.csv` — best model per peptide with parameters
- `kinetics_model_comparison.csv` — all models ranked by AICc
- `kinetics_fit_overlay.html` — interactive Plotly visualization
- `kinetics_predicted.csv` — fitted curves on fine time grid
- `kinetics_rate_summary.md` — summary statistics

**CLI:**
```bash
python skills/public/hdx-kinetics-fit-plot/scripts/fit_hdx_kinetics.py \
    --input output/02_hdx_reconstruction/reconstructed_peptide_timecourse.csv \
    --out-dir output/02_hdx_reconstruction/analysis_fit
```

**Parameters:**
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--min-points` | 4 | Minimum timepoints to attempt fit |
| `--max-grid-multiplier` | 20.0 | Extend prediction to 20× max observed time |

**Models tested (per peptide):**

| Model | Equation | Parameters |
|-------|----------|------------|
| Mono-exponential | `U = U_inf × (1 - e^(-k×t))` | U_inf, k |
| Stretched exponential | `U = U_inf × (1 - e^(-(k×t)^β))` | U_inf, k, β |
| Bi-exponential | `U = U_inf × (w × (1-e^(-k_fast×t)) + (1-w) × (1-e^(-k_slow×t)))` | U_inf, k_fast, k_slow, w |

**Model selection:** AICc (corrected Akaike Information Criterion). Lower is better. Delta AICc > 2 indicates meaningful difference.

**Fitting strategy:**
- Multiple initial guesses (grid search on k, U_inf)
- Bounds: k ∈ [1e-7, 10], U_inf ∈ [0, 2×max_uptake], β ∈ [0.2, 2.0]
- scipy.optimize.curve_fit with TRF method, maxfev=50000
- Collapses duplicate timepoints by mean before fitting

**Derived metrics:**
| Metric | Description |
|--------|-------------|
| t½ | Time to reach 50% of U_inf |
| t90 | Time to reach 90% of U_inf |
| k_app | Apparent rate = ln(2) / t½ |
| R² | Coefficient of determination |

**Results for TTR:**

| Peptide | Region | Model | R² | t½ (min) | Rate class |
|---------|--------|-------|----|----------|------------|
| GSGPTGTGESKCPLM | 1-15 | mono_exp | 0.005 | ~2 | Fast (N-term, flexible) |
| AVVTNPKE | 122-129 | mono_exp | 0.25 | ~3 | Fast (C-term, flexible) |
| KVEIDTKS | 72-79 | stretched | 0.99 | ~8 | Intermediate |
| FVEGIYKVEIDT | 66-77 | stretched | 0.97 | ~22 | Intermediate |
| ALLSPYSYSTT | 111-121 | mono_exp | 0.96 | ~75 | Slow (beta sheet) |
| VKVLDAVRGSPAIN | 16-29 | stretched | 0.97 | ~200 | Slow (buried core) |

---

## Step 5: Timecourse Analysis (optional)

Generate summary statistics, coverage plots, and heatmaps.

**Script:** `skills/public/hdx-kinetics-fit-plot/scripts/analyze_hdx_timecourse.py`

**Output:**
- `peptide_summary.csv` — early/late uptake, delta, AUC, t50
- `kinetics_lines.html` — line plot per peptide
- `uptake_heatmap.html` — heatmap (peptides × timepoints)

---

## Output Directory Structure

```
output/
├── 01_hdx_mzml/                        Step 1: mzML files
│   ├── 01_20260217_TTR_peptide_mapping.mzML
│   ├── 11_20260217_wt_TTR_10min.mzML
│   ├── 03_20260217_wt_TTR_30min.mzML
│   ├── 01_20260217_wt_TTR_3h.mzML
│   ├── 05_20260217_wt_TTR_6h.mzML
│   ├── 07_20260217_wt_TTR_12h.mzML
│   └── 09_20260217_wt_TTR_24h.mzML
├── 02_hdx_reconstruction/              Step 3: reconstruction + Step 4: kinetics
│   ├── run_manifest.csv
│   ├── peptide_targets.csv
│   ├── reconstructed_peptide_timecourse.csv
│   ├── reconstructed_peptide_summary.csv
│   └── analysis_fit/
│       ├── kinetics_fit.csv
│       ├── kinetics_model_comparison.csv
│       └── kinetics_fit_overlay.html
└── 03_hdx_peptide_discovery/           Step 2: discovery
    ├── discovered_peptides.csv
    ├── peptide_comparison.csv
    ├── residue_coverage.csv
    └── TTR.fasta
```
