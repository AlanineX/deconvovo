# Waters plotting, mobiligrams and CCS

This folder includes the Waters plotting/calibration modules and complete real protein/ADP and tune-mix/Betaine examples. It can be copied independently of the main checkout. The default tutorial fits the protein standards and applies the curve to ADP. For the workflow map and historical recovery audit, see [WATERS_PIPELINE_MAP.md](WATERS_PIPELINE_MAP.md).

## Install

Python 3.10+ is required. RAW conversion also needs Wine on Linux/macOS; see the [conversion README](README.md).

From `portable_waters_convert/`, on Linux/macOS:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements-pipeline.txt
```

On Windows, create the environment with `py -3 -m venv .venv`, then use `.venv\Scripts\python` instead of `.venv/bin/python` in these commands. The five plotting dependencies are NumPy, pandas, SciPy, Plotly and matplotlib. UniDec installation is unnecessary. For the desktop Waters workspace, install `requirements-gui.txt` and run `run_gui.py`; see [GUI_README.md](GUI_README.md). Its demo selector loads [seven complete acquired production runs](examples/production/README.md), totaling 106 original files and 235.75 MiB, with protein/ADP as the default and tune-mix/Betaine as a separate example.

## Protein standards → ADP CCS

```bash
.venv/bin/python waters_convert.py examples/production/protein_adp/raw -o converted --fn 1 --ms-bin 0 --im-bin 0
.venv/bin/python -m deconvovo.imms_plot -i converted -o html --raw-dir examples/production/protein_adp/raw --config examples/production/protein_adp/viewer_config.json -j 1
.venv/bin/python -m deconvovo.imms_ccs_calibrate -o ccs --data-dir converted --calibrant-csv examples/production/protein_adp/calibrants.csv --analyte-csv examples/production/protein_adp/analytes.csv --conversion-method direct
```

Three unfolded protein acquisitions supply 25 reference rows; the curve is applied to 34 target rows from two ADP solution conditions. The inherited protein fit is poor and many ADP peaks lie below its reference drift interval. This is a processing and QC example, not a validated ADP CCS measurement.

The copied HTML viewer embeds the Plotly JavaScript matching the installed Python package, so the generated HTML opens **offline**. It provides linked 2D/drift/mass-spectrum panels, smoothing/noise/intensity controls, display binning, colormaps, and PNG/CSV export.

## Betaine acquisition: 1000–2000 visual detail

This supplementary view uses the complete Betaine acquisition and supplies no
CCS tables. Its high-mass features are unassigned; the Betaine ion near 118 lies
outside the focus. Convert this separate run, then generate its focused viewer:

```bash
.venv/bin/python waters_convert.py examples/production/small_molecule/raw/20260320_BETAINE_V2.raw -o betaine_converted --fn 1 --ms-bin 0 --im-bin 0
.venv/bin/python -m deconvovo.imms_plot -i betaine_converted -o html/betaine_1000_2000 --raw-dir examples/production/small_molecule/raw --config examples/production/small_molecule/betaine_1000_2000_viewer_config.json -j 1
```

| Display choice | Betaine detail |
| --- | --- |
| Starting mass window | 1000–2000 |
| Starting drift window | 4–12 ms |
| Intensity scaling | Cube root |
| Mass / drift bins | 6400 / native |
| 2D Gaussian smoothing | 0.1 Da mass; 0.07 ms drift |
| 2D noise cutoff / fade | 5% / 1% |
| Mass / drift profiles | Profile MS; unsmoothed drift |

The full-batch overview and Betaine focus are display configurations. The heatmap
domain is baked from `initial_view.mz_range`; regenerate with a wider range to
inspect outside it. Display settings do not change the converted TXT inputs.
Choose `smooth_2d.method: raw` and `render: off` for an unsmoothed display.

`-i` accepts a folder or an explicit list of TXT paths sharing one parent folder. Bash-expanded globs work; a quoted wildcard is not expanded by this module. `--raw-dir` is the parent of the matching `.raw` folders, allowing drift-time units to be recovered. Start with `-j 1`; increase workers only if memory permits, because each worker loads a complete IM file and builds its own heatmap.

## Extract one exact target-mass mobiligram

```bash
.venv/bin/python extract_mobiligram.py "/path/to/run_im.txt" --mz 700.1 --mz-window 1 --raw-dir "/path/to/raw_parent" -o mobiligrams/target
```

The script reuses `ccs_peak_pick.extract_drift_profile`; it exports the unsmoothed sum within the inclusive target window, rather than an HTML heatmap marginal. The required settings are `--mz`, `--mz-window` (half-width; default `2`), and output prefix `-o`. RAW metadata supplies the pusher period; without it, output stays in drift bins and has no milliseconds column.

For mass-to-charge values $x_i$, target $x_0$, half-width $w$, drift bin $d$ and intensity $I_i$, the shared extractor computes:

$$
P_d = \sum_{i:\,d_i=d,\,|x_i-x_0|\leq w} I_i
$$

With pusher period $p$ in microseconds, the bin-center time in milliseconds is:

$$
t_d = \frac{(d+0.5)p}{1000}
$$

The metadata helper reads `_FUNC001.STS`; verify that the period applies if extracting another acquisition function.

## CCS calibration and analyte profiles

The first example fits `examples/production/protein_adp/calibrants.csv` and
applies the curve to `examples/production/protein_adp/analytes.csv` using the
same converted batch. Its protein standards and ADP tables belong together.

The included rows retain historical species assignments and peak choices. For your own data, edit `raw_path`, standards, charge states and conditions; these references are not automatically suitable for another acquisition. Omit `--analyte-csv` for calibration alone. The standalone templates `calibrants_protein_all.csv` and `calibrants_protein_unfolded.csv` are also available.

| CSV/CLI setting | Meaning |
| --- | --- |
| `raw_path` | Complete RAW folder; required for pusher period and EDC metadata; relative paths resolve beside the CSV, with existing working-directory paths retained for compatibility |
| Calibrant `name,mw,z,mz,ccs` | Species label, neutral mass, charge, target mass-to-charge value, reference CCS |
| Analyte `species,mw,z,mz` | Species label, neutral mass, charge and extraction target |
| `mz_window` | Numeric extraction half-width, or `auto` for the existing isotope-envelope approximation |
| `peak_select` | `highest` or `last`; calibrant peak choice |
| `min_ms_intensity,min_im_intensity,peak_height_frac` | Signal/peak eligibility settings; review them in the CSV |
| `gas_mw` | Drift-gas mass; default nitrogen `28.014`; calibrant rows must agree |
| `--conversion-method` | `direct` power-law conversion, or `twostep` diagnostic linear conversion |

At least **three valid calibrant peaks** are required after EDC correction. Use calibrants acquired under applicable instrument/gas conditions; a historical curve is not automatically valid for another experiment.

Betaine alone supplies one reference point and cannot fit the free two-parameter
TW-IMS curve. The separate tune-mix/Betaine example fits six Hexakis references
and applies that curve to Betaine; it does not calibrate ADP from Betaine.

Outputs include `calibration_curve.json`, `calibrant_summary.csv`, `calibration_plots.png`, calibrant drift plots/peak CSVs, `analyte_summary.csv`, and analyte CCS CSV/PNG folders for raw, smoothed and overlay views. Uniform CCS resampling uses a fixed `0.5` square-angstrom step; Gaussian smoothing uses `sigma=3` grid points. Both are existing constants in `deconvovo/ccs_convert.py`.

The [production guide](examples/production/README.md) records the complete protein/ADP and auxiliary tune-mix/Betaine examples, their relative-path tables and scientific limits. [Retained validation](examples/production/validation.json) records their earlier Linux/Wine execution.

## Current limitations

- The older browser raw/processed CSV buttons are distinct from the desktop's **Export current view** menu. Current-view profile CSVs use the displayed traces and current ranges; current 2D CSV uses transformed/filtered heatmap values. Keep the generated point CSV and original TXT/RAW data for raw intensities.
- Use `extract_mobiligram.py` for an exact target-window sum from raw IM points, independent of display bins, filtering and smoothing.
- Initial and rebinned heatmaps, point CSVs and CCS profiles use drift-bin centres, including the half-bin timing offset. Matched profiles follow the axis ranges Plotly actually renders.
- The calibration orchestrator rejects inconsistent EDC delay coefficients across calibrant RAW folders. Confirm compatible acquisition metadata before fitting.
- Release builds execute the protein-to-ADP workflow from source and the frozen Windows application, including space-containing paths, offline viewing, alignment and current-view exports. Successful releases attach [native validation](https://github.com/AlanineX/deconvovo/releases/download/v0.3.2/native_smoke.json). The retained data contains seven acquisitions and 106 original files; macOS native runtime remains untested.
