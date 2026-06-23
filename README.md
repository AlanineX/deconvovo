# DeconVoVo

An analysis suite for Waters SYNAPT travelling-wave ion mobility mass spectrometry (TW-IMS) data — from raw instrument files to collision cross-sections.

## Getting started

**Windows** — download `DeconVoVo-win64.zip` from [Releases](../../releases), extract, and double-click `DeconVoVo.exe`.

**Linux / macOS** — install with pip and launch the GUI or use the command line:
```bash
pip install -e ".[all]"
python -m gui.app
```

> Waters `.raw` conversion uses UniDec's `CDCReader.exe`. On Windows it runs
> natively; on Linux/macOS install `wine` first (`sudo apt install wine` /
> `brew install --cask wine-stable`). All other CLI steps are pure Python and
> work on any OS without Wine.

## Workflow

1. **Convert** — extract MS and IM data from Waters `.raw` directories to text files
2. **Visualize** — generate interactive 2D IM-MS heatmaps (HTML) with linked drift and m/z panels, smoothing controls, and CSV export
3. **CCS Calibrate** — build TW-IMS CCS calibration curves from calibrant species. Supports both protein and small-molecule calibrants (e.g. Agilent tune mix). Direct (power-law) and two-step fitting methods
4. **Analyte CCS** — apply calibration to compute CCS profiles, with resampling and smoothing on a uniform CCS grid

All parameters are CSV-driven — m/z extraction windows, peak selection strategy, intensity thresholds, and drift gas mass are specified per-row in the config. No hardcoded assumptions. See `config/` for examples.

## Command line

Every step runs as `python -m deconvovo.<module>`. Pass `--help` to any of them
for the full flag list.

```bash
# 1. Convert Waters .raw → text (_ms.txt + _im.txt)
python -m deconvovo.waters_convert data/*.raw -o converted/

# 2. Build interactive 2D IM-MS HTML viewers
python -m deconvovo.imms_plot -i converted/ -o html/ --raw-dir data/ -j 8

# 3. (optional) CCS calibration + analyte CCS in one step
python -m deconvovo.imms_ccs_calibrate \
    -o output/ccs \
    --data-dir data/ \
    --calibrant-csv config/calibrants_tunemix.csv \
    --analyte-csv config/analytes_adp.csv

# 4. Per-run peak summary CSV
python -m deconvovo.imms_summary -i html/
```

### One-shot pipeline

`deconvovo.cli` runs steps 1 + 2 end-to-end (and 3/4 if their inputs are
present):

```bash
python -m deconvovo.cli -i data/ -o output/ \
    --mass-range 5000 100000 \
    --charge-range 1 30 \
    --mass-bins 10 \
    -j 8
```


## Command line for me

Uses the custom config at `/home/alan/Downloads/carter_0623/imms_plot_config.json`
(starting m/z 350–1250, drift 0–15 ms, PuBu colormap, 0.3 σ 2D Gaussian smooth,
0.5 % 2D noise floor, 6400 m/z bins). Edit that file to change any default —
the dropdowns in the HTML still work for runtime tuning.

One-shot — convert + viewers in a single invocation:

```bash
python -m deconvovo.cli \
    -i /home/alan/Downloads/carter_0623/Myoglobin_Water \
    -o /home/alan/Downloads/carter_0623/Myoglobin_Water_output \
    --config /home/alan/Downloads/carter_0623/imms_plot_config.json \
    --skip-deconv \
    -j 8
```

Or step-by-step if you want to inspect the intermediate text files:

```bash
# 1. Convert .raw to text
python -m deconvovo.waters_convert \
    /home/alan/Downloads/carter_0623/Myoglobin_Water \
    -o /home/alan/Downloads/carter_0623/Myoglobin_Water_converted

# 2. Generate interactive HTML viewers (custom config)
python -m deconvovo.imms_plot \
    -i /home/alan/Downloads/carter_0623/Myoglobin_Water_converted \
    -o /home/alan/Downloads/carter_0623/Myoglobin_Water_html \
    --raw-dir /home/alan/Downloads/carter_0623/Myoglobin_Water \
    --config /home/alan/Downloads/carter_0623/imms_plot_config.json \
    -j 8

# 3. (optional) CCS calibration + analyte CCS
python -m deconvovo.imms_ccs_calibrate \
    -o /home/alan/Downloads/carter_0623/Myoglobin_Water_ccs \
    --data-dir /home/alan/Downloads/carter_0623/Myoglobin_Water \
    --calibrant-csv config/calibrants_tunemix.csv \
    --analyte-csv config/analytes_adp.csv
```