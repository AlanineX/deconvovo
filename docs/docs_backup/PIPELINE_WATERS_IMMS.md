# Waters IM-MS Pipeline

Complete workflow for Waters SYNAPT travelling-wave ion mobility mass spectrometry data.

## Overview

```
.raw files → Convert → Text files → HTML Viewer → CCS Calibration → Analyte CCS
```

## Step 1: Convert .raw to Text

Extracts MS and IM data from Waters MassLynx .raw directories.

**Module:** `deconvovo.waters_convert`, `deconvovo.imms_convert`

**Tool:** CDCReader.exe (from UniDec, runs natively on Windows, via Wine on Linux)

**Input:** Directory containing `.raw` folders

**Output per run:**
- `{run}_ms.txt` — MS spectrum (2 columns: m/z, intensity)
- `{run}_im.txt` — IM data (3 columns: m/z, drift_bin, intensity)

**CLI:**
```bash
python -m deconvovo.waters_convert data/*.raw -o output/converted
python -m deconvovo.imms_convert -i data/ -o output/converted -j 4
```

**Parameters:**
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--fn` | 1 | MassLynx function number |
| `--ms-bin` | 0 | m/z bin size for MS (0 = raw resolution) |
| `--im-bin` | 0 | m/z bin size for IM (0 = raw resolution) |
| `-j` | 4 | Parallel workers |

**Notes:**
- CDCReader reads the binary `.DAT`/`.IDX`/`.STS` files inside each `.raw` directory
- The pusher period (drift time calibration) is read from `_FUNC001.STS` (stat code 76)
- EDC delay coefficient is read from `_extern.inf`
- Drift bins are integer indices; drift time = `(bin + 0.5) × pusher_μs / 1000` ms

---

## Step 2: Interactive 2D IM-MS Viewer

Generates self-contained HTML files with linked heatmap, drift profile, and m/z panels.

**Module:** `deconvovo.imms_plot`, `deconvovo.imms_html`

**Input:** Converted text directory + optional .raw directory (for pusher period)

**Output per run:**
- `{run}_2d_imms.html` — Interactive viewer (opens in browser)
- `{run}_2d_imms.csv` — Raw data export (m/z, drift_bin, drift_time_ms, intensity)

**CLI:**
```bash
python -m deconvovo.imms_plot -i output/converted -o output/html --raw-dir data/ -j 8
```

**Parameters:**
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--raw-dir` | None | Optional: .raw directory for pusher period → drift time axis |
| `--skip-existing` | off | Skip runs with existing HTML (checks file size > 0) |
| `-j` | 8 | Parallel workers |

**HTML viewer features:**
- 2D heatmap (center): drift time × m/z, log-scaled intensity
- Drift profile (left): linked to zoom, full-range TIC from heatmap
- m/z spectrum (bottom): raw MS data, independent of drift range
- Client-side smoothing: Gaussian, Savitzky-Golay, moving average
- Noise filtering: adjustable threshold (% of max)
- Colormap selection: Viridis, Plasma, Inferno, Hot, etc.
- Export: PNG screenshot, CSV for drift profile and m/z spectrum
- rAF-debounced zoom handler for smooth interaction

**Config:** `config/imms_plot_config.json` — presets for smoothing, colormaps, figure size

---

## Step 3: CCS Calibration

Builds a TW-IMS CCS calibration curve from known calibrant species.

**Module:** `deconvovo.imms_ccs_calibrate`

**Sub-modules:**
- `deconvovo.ccs_peak_pick` — drift profile extraction, peak detection, auto m/z window
- `deconvovo.ccs_fit` — power-law fitting (honest, no outlier rejection)
- `deconvovo.ccs_convert` — CCS conversion, resampling, smoothing
- `deconvovo.ccs_plot` — matplotlib plotting

**Input:**
- `--data-dir` — one or more directories with converted text files
- `--calibrant-csv` — species parameters (name, mw, z, mz, ccs, raw_path, ...)
- `--analyte-csv` — optional analyte species for CCS computation

**Output:**
- `calibration_curve.json` — fit parameters (X, A, R², EDC, gas_mw)
- `calibrant_summary.csv` — per-calibrant results
- `calibration_plots.png` — 2-panel fit plot (ln-ln + diagnostic)
- `cali_png/` — drift profile per calibrant
- `cali_csv/` — peak detection details per calibrant
- `analyte_summary.csv` — per-analyte CCS values
- `csv_raw/`, `csv_smoothed/` — CCS profiles
- `png_raw/`, `png_smoothed/`, `png_overlay/` — CCS profile plots

**CLI:**
```bash
python -m deconvovo.imms_ccs_calibrate \
    -o output/ccs \
    --data-dir output/converted_tunemix output/converted_adp \
    --calibrant-csv config/calibrants_tunemix.csv \
    --analyte-csv config/analytes_adp.csv \
    --conversion-method direct
```

**Calibrant CSV columns:**

| Column | Required | Default | Description |
|--------|----------|---------|-------------|
| `name` | Yes | — | Species name |
| `formula` | No | — | Molecular formula (for auto m/z window) |
| `mw` | Yes | — | Molecular weight (Da) |
| `z` | Yes | — | Charge state |
| `mz` | Yes | — | Observed m/z |
| `ccs` | Yes | — | Literature CCS (Å²) |
| `raw_path` | Yes | — | Path to .raw directory (for pusher + EDC) |
| `mz_window` | No | auto | Extraction window (Da) or "auto" |
| `peak_select` | No | highest | "highest" (tallest peak) or "last" (most extended conformer) |
| `min_ms_intensity` | No | 100 | MS detection threshold |
| `min_im_intensity` | No | 500 | IM peak intensity threshold |
| `peak_height_frac` | No | 0.10 | Minimum peak height as fraction of tallest |
| `gas_mw` | No | 28.014 | Drift gas molecular weight (N₂) |

**CCS methods:**
- **direct** (recommended): `CCS = A × t'_D^X × z × μ` — power-law, no intercept
- **twostep**: `CCS = slope × t''_D + intercept` — linear, has intercept artifact for small molecules

**Auto m/z window:** Computed from molecular formula isotopic envelope (3σ width) or averagine model. Fluorinated compounds (tune mix) get tighter windows than averagine predicts.

**Peak selection:** For denatured proteins with TW-IMS rollover, use `peak_select=last` to pick the most extended conformer at each charge state.

**Honest fitting:** No automated outlier rejection. User controls which points to include by editing the calibrant CSV. All provided points are fit.

---

## Step 4: Analyte CCS

Apply calibration to compute CCS for target analyte species.

Uses the same `imms_ccs_calibrate` module with `--analyte-csv`.

**Analyte CSV columns:**

| Column | Required | Default | Description |
|--------|----------|---------|-------------|
| `species` | Yes | — | Species name |
| `mw` | Yes | — | Molecular weight (Da) |
| `z` | Yes | — | Charge state |
| `mz` | Yes | — | Observed m/z |
| `raw_path` | Yes | — | Path to .raw directory |
| `mz_window` | No | auto | Extraction window |
| `min_ms_intensity` | No | 50 | MS detection threshold |

**CCS profile processing:**
1. Extract drift profile (sum intensity within ±mz_window per drift bin)
2. Convert each bin to CCS via calibration curve
3. Resample to uniform CCS grid (0.5 Å² step, anchored at multiples)
4. Gaussian smooth (σ = 3.0 Å²)
5. Output raw + smoothed profiles as CSV and PNG

---

## Physical Constants

| Constant | Value | Used for |
|----------|-------|----------|
| `H_MASS` | 1.00728 Da | m/z ↔ MW conversion |
| `GAS_MW_N2` | 28.014 Da | Reduced mass in CCS equation |
| `NEUTRON` | 1.00336 Da | Isotopic envelope spacing |

## Key Equations

**EDC correction:** `t'_D = t_D - C × √(m/z) / 1000`

**Reduced mass factor:** `μ = √(1/MW + 1/gas_MW)`

**Reduced CCS:** `Ω' = CCS / (z × μ)`

**Power-law fit:** `ln(Ω') = X × ln(t'_D) + ln(A)`

**Direct CCS:** `CCS = A × (t'_D)^X × z × μ`
