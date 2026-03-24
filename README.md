# DeconVoVo

Portable analysis suite for Waters SYNAPT travelling-wave ion mobility mass spectrometry (TW-IMS) data. Converts raw instrument files, visualizes 2D IM-MS data, and calibrates collision cross-sections (CCS).

## Download

**Windows:** grab `DeconVoVo-win64.zip` from [Releases](../../releases). Extract and run `DeconVoVo.exe` — nothing else to install.

## What it does

- **Convert** — Waters `.raw` → text (MS + IM) via bundled CDCReader
- **Visualize** — interactive 2D IM-MS heatmaps (HTML, opens in browser) with linked drift and m/z panels, client-side smoothing, and CSV export
- **CCS Calibrate** — build TW-IMS CCS calibration curves from a CSV of calibrants. Supports protein and small-molecule (tune mix) calibrants. Direct (power-law) and two-step methods
- **Analyte CCS** — apply calibration to compute CCS profiles for target species
- **Deconvolve** — intact mass deconvolution via bundled UniDec engine

All parameters are CSV-driven: m/z extraction windows, peak selection, intensity thresholds, and drift gas mass are specified per-row. No hardcoded assumptions.

## GUI

Dark-themed desktop interface with editable CSV tables, embedded calibration plots, real-time log, and background execution. Five workflow panels: Convert, Viewer, Calibration, Analyte CCS, Full Pipeline.

## CLI

```bash
pip install -e ".[gui,deconv]"

# CCS calibration
python -m deconvovo.imms_ccs_calibrate \
    -o output/ccs --calibrant-csv config/calibrants_tunemix.csv \
    --analyte-csv config/analytes_adp.csv

# HTML viewer
python -m deconvovo.imms_plot -i converted/ -o html/ --raw-dir data/ -j 8
```

## Project layout

`deconvovo/` — core library (IO, CCS fitting, plotting, vendor readers)
`gui/` — PySide6 desktop application
`config/` — example calibrant and analyte CSVs
`demo_data/` — sample Waters `.raw` files

## References

Ruotolo & Robinson (2008) *Nat. Protoc.* 3, 1139. Bush et al. (2010) *Anal. Chem.* 82, 9557.
