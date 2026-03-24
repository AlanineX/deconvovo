# DeconVoVo

An analysis suite for Waters SYNAPT travelling-wave ion mobility mass spectrometry (TW-IMS) data — from raw instrument files to collision cross-sections.

## Getting started

**Windows** — download `DeconVoVo-win64.zip` from [Releases](../../releases), extract, and double-click `DeconVoVo.exe`.

**Linux / macOS** — install with pip and launch the GUI or use the command line:
```bash
pip install -e ".[gui]"
python -m gui.app
```

## Workflow

1. **Convert** — extract MS and IM data from Waters `.raw` directories to text files
2. **Visualize** — generate interactive 2D IM-MS heatmaps (HTML) with linked drift and m/z panels, smoothing controls, and CSV export
3. **CCS Calibrate** — build TW-IMS CCS calibration curves from calibrant species. Supports both protein and small-molecule calibrants (e.g. Agilent tune mix). Direct (power-law) and two-step fitting methods
4. **Analyte CCS** — apply calibration to compute CCS profiles, with resampling and smoothing on a uniform CCS grid

All parameters are CSV-driven — m/z extraction windows, peak selection strategy, intensity thresholds, and drift gas mass are specified per-row in the config. No hardcoded assumptions. See `config/` for examples.

## Command line

```bash
# Convert .raw to text
python -m deconvovo.waters_convert data/*.raw -o converted/

# Generate interactive HTML viewers
python -m deconvovo.imms_plot -i converted/ -o html/ --raw-dir data/ -j 8

# CCS calibration + analyte CCS in one step
python -m deconvovo.imms_ccs_calibrate \
    -o output/ccs \
    --calibrant-csv config/calibrants_tunemix.csv \
    --analyte-csv config/analytes_adp.csv
```
