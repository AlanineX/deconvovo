# DeconVoVo · Waters Instrument Console

Convert complete Waters `.raw` folders, inspect IM-MS data, and fit/apply CCS in
one local workspace. The desktop viewer and exported single-file HTML work offline.

Data flow: RAW folders → paired MS/IM tables → viewer → optional calibration → analyte CCS.

## Open the GUI

The v0.3.1 Windows package is awaiting final verification and publication on
[GitHub Releases](https://github.com/AlanineX/deconvovo/releases). For a matching
release, extract the whole `WatersConsole-win64.zip`, keep its executable and
`_internal` folder together, and open `WatersConsole/WatersConsole.exe`.
The compiled application needs no Python, Wine, or internet connection to run.

From source, use Python 3.10+ in your environment:

```bash
python -m pip install -e ".[gui]"
python -m gui.app
```

Windows runs the bundled reader natively. Linux/macOS RAW conversion needs
compatible x64 Wine; plotting and CSV processing use Python.

## Process and view

1. Choose a complete RAW folder, a parent with acquisitions, or a converted TXT folder.
2. Choose a writable **Output** folder outside the RAW acquisitions. For RAW sources, **Convert only** creates MS/IM tables without HTML, CCS, or reference CSVs.
3. Click **Process & view** for conversion and viewing, then select a run to inspect its measured plots.
4. Adjust the viewer controls and the visible **X · m/z** / **Y · drift time (ms)** minimum and maximum fields. **Advanced processing & defaults → Save defaults** keeps display choices across relaunches.
5. Use **Open in browser** for the full viewer, or **Export portable HTML** to save the selected run with current defaults. That one file includes Plotly and its data; copy it alone to another computer.
6. For CCS, open **CCS review**, load/edit the calibrant and analyte tables, then click **Calibrate & apply**. Alternatively, enable **Include CCS when processing** before starting the batch.

**Export current view** saves the visible full plot, heatmap, MS profile, or drift
profile as PNG; its CSV choices save the visible heatmap, MS, or drift values.
The 2D CSV contains displayed values after transforms/filters, rather than raw
point intensities. **Export portable HTML** regenerates an interactive document
with current defaults; it does not capture temporary zoom changes.

**Load real demo → Protein calibrants + ADP CCS** selects three protein standards
and two ADP acquisitions. The separate tune-mix/Betaine profile has two runs.
Together they contain **seven complete RAW acquisitions, 106 original files,
235.75 MiB**. [Demo details](demo_data/production/README.md) include the original-file provenance.

The protein fit is poor and many ADP peaks extrapolate below its reference
interval; review peak assignments and acquisition conditions before interpreting
CCS. Betaine also extrapolates. Betaine alone supplies one reference and cannot
fit the free two-parameter curve; this implementation needs at least three valid
reference peaks. Use the protein standards for the ADP example.

For the **Betaine acquisition's** high-mass detail, load **Betaine · 1000–2000 m/z
detail** or regenerate its viewer using the
[focused configuration](demo_data/production/small_molecule/betaine_1000_2000_viewer_config.json).
This visual-only view uses 1000–2000 and 4–12 ms; its high-mass features are
unassigned, and the Betaine ion near 118 lies outside the focus. Protein → ADP
remains the primary CCS example. Display settings do not alter calibration inputs.

## Write calibrant and analyte CSVs

Save UTF-8, comma-separated files with a header row. Each calibrant row describes
one reference ion/charge state with known CCS; each analyte row describes one
extraction target. The bundled [calibrant](demo_data/production/protein_adp/calibrants.csv)
and [analyte](demo_data/production/protein_adp/analytes.csv) tables are editable starting points.

| Required column | Calibrants | Analytes | Meaning |
| --- | --- | --- | --- |
| `name` | Yes | — | Reference label; not parsed to infer identity or charge. |
| `species` | — | Yes | Target label; not parsed to infer identity or charge. |
| `mw` | Yes | Yes | Neutral molecular mass in Da; used for mass correction and automatic extraction width. |
| `z` | Yes | Yes | Positive integer charge state, such as `8` or `1`. |
| `mz` | Yes | Yes | Target mass-to-charge value in Th; supplied explicitly. |
| `ccs` | Yes | — | Known reference CCS in Å²; analyte CCS is calculated from the fitted curve. |
| `raw_path` | Yes | Yes | Complete acquisition directory ending in `.raw`, including its metadata. |

Minimal row formats, taken from the demo:

```csv
name,mw,z,mz,ccs,raw_path
UBQ z=8,8565,8,1071.63,1622,raw/20260216_UBQ_DEN_NEW_V1.raw
```

```csv
species,mw,z,mz,raw_path
[ADP+H]+,427.2,1,428.21,raw/20260212_200UM_ADP_1M_EDDA_PH6_0.raw
```

The single calibrant row shows the format; it cannot fit a curve by itself.
Start with the supplied 25-reference table for the real ADP demo.

### RAW paths

Relative `raw_path` values resolve beside the CSV. If the CSV is in `my_batch/`,
`raw/sample.raw` points to `my_batch/raw/sample.raw`; `raw\sample.raw` also works.
The corresponding converted files must be named `sample_ms.txt` and `sample_im.txt`.

Windows absolute examples are `C:\Waters data\sample.raw` and
`\\server\share\sample.raw`. In a CSV, use literal backslashes, without Python-style
doubling. A field containing a comma must be enclosed in double quotes; spaces
may also be quoted: `"C:\Waters data\sample.raw"`. Quote complete paths with spaces
when passing CLI arguments. Keep the original RAW folders reachable for drift
period and EDC metadata; an MS/IM TXT pair alone does not supply that metadata.

### Optional columns

Omit an optional setting or leave it blank to use its listed default.
An omitted or blank `formula` uses the mass-based automatic window estimate;
when supplied, use an elemental formula for that estimate.

| Column | Applies to | Default | Meaning / accepted values |
| --- | --- | --- | --- |
| `mz_window` | Both | `auto` | Extraction half-width in Th, such as `2`, or `auto` for the isotope-envelope estimate. |
| `formula` | Both | Blank | Elemental formula used only when `mz_window` is `auto`; blank uses the mass-based estimate. |
| `min_ms_intensity` | Both | Calibrants `100`; analytes `50` | Minimum maximum MS intensity inside the target window, in reader intensity units. |
| `peak_select` | Calibrants | `highest` | `highest`: tallest accepted drift peak; `last`: accepted peak at the latest drift bin. |
| `min_im_intensity` | Calibrants | `500` | Minimum peak height of the summed IM drift profile, in reader intensity units. |
| `peak_height_frac` | Calibrants | `0.10` | Peak-detection height as a fraction of the profile maximum; normally 0–1. |
| `gas_mw` | Calibrants | `28.014` | Drift-gas molecular mass in Da; nitrogen default. All references must agree. |

The program uses the supplied `mw`, `z`, and `mz`; it does not calculate `mz`
from the label or mass. It measures drift peaks and intensities from the data,
reads drift timing/EDC from RAW metadata, fits the references, and calculates
analyte CCS. `min_im_intensity`, `peak_height_frac`, and `peak_select` do not
control analyte extraction. `direct` is the power-law CCS method; `twostep` uses
the diagnostic linear conversion.

In the GUI, **Load CSV**, direct cell editing, **+ Row**, and **- Row** prepare the
full-width **Calibrant sheet** and **Analyte sheet** tabs. Scroll horizontally,
drag column edges to resize, and hover over
headers/cells for field meanings and full values. **Save CSV** saves your
edits; processing otherwise uses a temporary snapshot without rewriting the
shipped reference CSV. Use **Save CSV** to a separate file for your own study.

## Brief CLI

Run from the repository; use `--help` for options. On Windows, `python` can be
replaced with `py -3`.

```bash
python -m deconvovo.waters_convert demo_data/production/protein_adp/raw -o output/converted --fn 1 --ms-bin 0 --im-bin 0
python -m deconvovo.imms_plot -i output/converted -o output/html --raw-dir demo_data/production/protein_adp/raw --config demo_data/production/protein_adp/viewer_config.json -j 1
python -m deconvovo.imms_ccs_calibrate -o output/ccs --data-dir output/converted --calibrant-csv demo_data/production/protein_adp/calibrants.csv --analyte-csv demo_data/production/protein_adp/analytes.csv --conversion-method direct
```

Function `1` selects `_FUNC001`; confirm its role in the acquisition method.
Zero mass-bin widths retain reader resolution. The TXT tables combine scans;
keep the original RAW data. The existing batch CLI is `python -m deconvovo.cli`;
optional UniDec deconvolution remains available through `.[deconv]`.

The [portable toolkit](portable_waters_convert/README.md) runs independently.
The complete [HTML viewer data flow and math](docs/html_viewer_data_flow_and_math.md)
explains binning, smoothing, intensity transforms, cutoff/fade, and linked profiles.
[Reader licenses](portable_waters_convert/bin/) retain the third-party terms.
