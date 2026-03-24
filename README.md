# DeconVoVo

**Portable IM-MS analysis suite for Waters SYNAPT travelling-wave ion mobility mass spectrometry data.**

DeconVoVo provides a complete pipeline from raw instrument files to collision cross-section (CCS) values, with an interactive desktop GUI — no Python installation required for Windows users.

---

## Features

### 1. Data Conversion
Convert Waters MassLynx `.raw` directories to analysis-ready text files using CDCReader (bundled).
- MS spectrum: 2-column (m/z, intensity)
- IM data: 3-column (m/z, drift_bin, intensity)

### 2. Interactive 2D IM-MS Viewer
Generate self-contained HTML files with:
- 2D heatmap (drift time × m/z) with linked marginal panels
- Client-side smoothing, noise filtering, colormap selection
- Zoom-linked drift profile and m/z spectrum panels
- Export to PNG and CSV

### 3. CCS Calibration
Build TW-IMS CCS calibration curves following the Ruotolo & Robinson (2008) protocol:
- **CSV-driven** — all parameters (m/z window, peak selection, thresholds) specified per-row
- **Auto m/z window** from molecular formula or averagine model
- **Honest fitting** — no automated outlier rejection; user controls point selection via CSV
- **Two methods**: direct (power-law, recommended) and two-step (linear diagnostic)
- Supports both protein calibrants and custom small-molecule calibrants (e.g., Agilent tune mix)

### 4. Analyte CCS Computation
Apply calibration to compute CCS for analyte species:
- Per-species drift profile extraction and CCS conversion
- Uniform CCS grid resampling with Gaussian smoothing
- Raw, smoothed, and overlay CCS profile plots
- Summary CSV with CCS values, drift times, and detection status

### 5. UniDec Deconvolution
Mass deconvolution via UniDec engine (bundled):
- Automated peak detection and mass assignment
- Configurable mass/charge ranges via CLI

### 6. Desktop GUI
Dark-themed scientific interface (PySide6/Qt):
- Sidebar navigation across all workflow steps
- Built-in CSV editor for calibrant and analyte tables
- Embedded matplotlib calibration plots
- Real-time log panel with color-coded messages
- Background execution with progress indicators
- Full pipeline mode: run everything end-to-end

---

## Installation

### Windows (portable executable)
Download `DeconVoVo-win64.zip` from [Releases](../../releases), extract, and run `DeconVoVo.exe`.
No Python or other software required.

### Python (development)
```bash
pip install -e ".[gui,deconv]"

# Launch GUI
python -m gui.app

# Or use CLI
python -m deconvovo.imms_ccs_calibrate \
    -o output/ccs \
    --calibrant-csv config/calibrants_tunemix.csv \
    --analyte-csv config/analytes_adp.csv
```

**Requirements:** Python ≥ 3.10, numpy, scipy, pandas, matplotlib, plotly. Optional: PySide6 (GUI), unidec (deconvolution).

---

## Quick Start

### CCS Calibration (CLI)
```bash
# 1. Convert .raw files
python -m deconvovo.waters_convert data/*.raw -o output/converted

# 2. Generate HTML viewers
python -m deconvovo.imms_plot -i output/converted -o output/html --raw-dir data/ -j 8

# 3. Run CCS calibration with analytes
python -m deconvovo.imms_ccs_calibrate \
    -o output/ccs_direct \
    --calibrant-csv config/calibrants_tunemix.csv \
    --analyte-csv config/analytes_adp.csv \
    --conversion-method direct
```

### Configuration CSVs

**Calibrant CSV** (required columns: `name`, `mw`, `z`, `mz`, `ccs`, `data_dir`, `raw_path`):
| Column | Description |
|--------|-------------|
| `name` | Species name |
| `formula` | Molecular formula (for auto m/z window) |
| `mw` | Molecular weight (Da) |
| `z` | Charge state |
| `mz` | Observed m/z |
| `ccs` | Literature CCS (Å²) |
| `data_dir` | Path to converted text files |
| `raw_path` | Path to .raw directory |
| `mz_window` | Extraction window (Da) or `auto` |
| `peak_select` | `highest` or `last` (most extended conformer) |
| `min_ms_intensity` | MS detection threshold |
| `min_im_intensity` | IM peak threshold |
| `peak_height_frac` | Minimum peak height as fraction of tallest |
| `gas_mw` | Drift gas molecular weight (default: 28.014 for N₂) |

See `config/` for example CSVs.

---

## Project Structure
```
deconvovo/           # Core library (vendor-agnostic)
  io.py              # Shared I/O, constants, m/z utilities
  vendors/waters.py  # Waters-specific binary readers
  ccs_peak_pick.py   # Drift profile extraction, peak detection
  ccs_fit.py         # Power-law calibration fitting
  ccs_convert.py     # CCS conversion, resampling, smoothing
  ccs_plot.py        # Matplotlib plotting
  imms_ccs_calibrate.py  # CCS orchestrator
  imms_html.py       # Interactive HTML viewer (800+ lines JS)
  imms_plot.py       # Batch HTML generation
  imms_deconv.py     # UniDec deconvolution
  waters_convert.py  # CDCReader wrapper
  smooth.py          # Shared smoothing utilities
  parallel.py        # Multiprocessing map
  taskqueue.py       # DAG task scheduler
  cli.py             # Full pipeline CLI

gui/                 # Desktop GUI (PySide6/Qt)
  app.py             # Main window
  theme.py           # Dark scientific theme
  widgets/           # Reusable widgets (CSV table, plot canvas, log, file picker)
  panels/            # Workflow panels (convert, HTML, CCS calibration, analyte, full)

config/              # Example CSV configurations
demo_data/           # Sample Waters .raw files for testing
docs/                # Mathematical derivations and method comparisons
build/               # PyInstaller spec for Windows packaging
```

---

## References

- Ruotolo, B.T. & Robinson, C.V. (2006). *Nat. Protoc.* **3**, 1139–1152. TW-IMS CCS calibration protocol.
- Bush, M.F. et al. (2010). *Anal. Chem.* **82**, 9557–9565. Unified CCS Compendium.

## License

MIT
