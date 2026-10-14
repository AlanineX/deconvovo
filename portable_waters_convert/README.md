# Portable Waters Instrument Console

Copy this complete folder to use Waters RAW conversion, the local viewer,
single-file offline HTML, and CCS calibration independently of DeconVoVo.

RAW folders → paired MS/IM tables → native plots → optional CCS / separate HTML export.

## Open the GUI

Download the Windows package from
[GitHub Releases](https://github.com/AlanineX/deconvovo/releases). Its compiled ZIP
bundles Python, Qt, Matplotlib, NumPy/SciPy, the HTML-export assets, and the reader.
Extract the whole folder and open
`WatersConsole/WatersConsole.exe`, keeping `_internal` beside it.

For this source toolkit, use Python 3.10+ in your environment:

```bash
python -m pip install -r requirements-gui.txt
python run_gui.py
```

Use `py -3` instead of `python` on Windows if needed. The included reader runs
natively on Windows x64; Linux/macOS RAW conversion requires compatible x64 Wine.
After installing dependencies, conversion, both viewers, and CCS run locally.

1. Choose a complete `.raw` directory, a parent with acquisitions, or a converted TXT folder.
2. Choose a writable output folder outside RAW acquisitions. For RAW sources, **Convert only** creates MS/IM tables without HTML, CCS, or reference CSVs; **Process & view** also opens native plots without generating HTML.
3. Select a run, adjust display controls and the **X · m/z** / **Y · drift time (ms)** minimum/maximum fields, and inspect locally. **Reset view** restores numeric ranges while keeping current filters.
4. **Export portable HTML** saves the selected run with current display defaults into one file containing Plotly and the data. Copy that file alone to another computer.
5. Open **CCS review**, load/edit the two tables, and **Calibrate & apply**; or enable CCS during processing. **Advanced processing & defaults → Save defaults** keeps display choices across relaunches.

**Export current view** offers PNG for the full plot, heatmap, MS, and drift,
plus CSV for the visible heatmap/MS/drift values. The heatmap CSV contains
displayed transformed/filtered values. Portable HTML keeps current defaults,
without saving temporary zoom changes.
Native viewing reads the converted TXT data and uses Matplotlib with NumPy/SciPy,
without a browser runtime. Native ranges can widen directly from those inputs.
The separate HTML export embeds Plotly and opens in an external browser offline.

**Load real demo → Protein calibrants + ADP CCS** uses three protein standards
and two ADP acquisitions. Tune mix/Betaine is a separate two-run example. The
[seven complete acquisitions](examples/production/README.md) contain **106 original
files, 235.75 MiB**. The protein fit is poor and ADP/Betaine extrapolate below their
reference drift intervals; review assignments and acquisition conditions before
interpreting CCS. Betaine alone provides one reference and cannot fit the free
two-parameter curve; at least three valid reference peaks are required.

For the **Betaine acquisition's** visual-only 1000–2000 detail, choose
**Betaine · 1000–2000 m/z detail** or use
[its focused configuration](examples/production/small_molecule/betaine_1000_2000_viewer_config.json).
The view spans 4–12 ms with cube-root intensity and 6400 mass bins. High-mass
features are unassigned; the Betaine ion near 118 lies outside this focus.
Protein → ADP remains the primary CCS example; its batch overview spans 100–2000.

## Write CSV inputs

Use UTF-8, comma-separated files with a header. One calibrant row identifies a
reference ion/charge state with known CCS; one analyte row identifies an extraction
target. Start with the bundled [protein calibrants](examples/production/protein_adp/calibrants.csv)
and [ADP analytes](examples/production/protein_adp/analytes.csv).

| Required column | Calibrants | Analytes | Meaning |
| --- | --- | --- | --- |
| `name` | Yes | — | Reference label; not parsed for identity or charge. |
| `species` | — | Yes | Target label; not parsed for identity or charge. |
| `mw` | Yes | Yes | Neutral molecular mass in Da, used for mass correction and automatic extraction width. |
| `z` | Yes | Yes | Positive integer charge state. |
| `mz` | Yes | Yes | Explicit target mass-to-charge value in Th. |
| `ccs` | Yes | — | Known reference CCS in Å²; analyte CCS is calculated. |
| `raw_path` | Yes | Yes | Complete `.raw` acquisition directory, including metadata. |

Minimal row formats from the demo:

```csv
name,mw,z,mz,ccs,raw_path
UBQ z=8,8565,8,1071.63,1622,raw/20260216_UBQ_DEN_NEW_V1.raw
```

```csv
species,mw,z,mz,raw_path
[ADP+H]+,427.2,1,428.21,raw/20260212_200UM_ADP_1M_EDDA_PH6_0.raw
```

One calibrant row cannot fit a curve. Use the supplied 25-reference protein table
for the real ADP workflow, then review the fit and reference peak selections.

### Paths and optional columns

Relative `raw_path` values resolve beside their CSV: `raw/sample.raw` means
`raw/sample.raw` under the CSV folder. `raw\sample.raw` also works. Converted data
must contain `sample_ms.txt` and `sample_im.txt` for that acquisition.
Windows absolute paths can be `C:\Waters data\sample.raw` or `\\server\share\sample.raw`.
Use literal backslashes in CSV; do not double them as Python escapes. Quote a
field containing a comma, and optionally a field with spaces:
`"C:\Waters data\sample.raw"`. Quote CLI path arguments containing spaces.
Keep RAW folders reachable for drift timing and EDC metadata.

| Optional column | Applies to | Default | Meaning / accepted values |
| --- | --- | --- | --- |
| `mz_window` | Both | `auto` | Target extraction half-width in Th, such as `2`, or isotope-envelope estimate `auto`. |
| `formula` | Both | Blank | Elemental formula used only by automatic window estimation; blank uses the mass-based estimate. |
| `min_ms_intensity` | Both | Calibrants `100`; analytes `50` | Minimum maximum MS intensity in the target window, in reader units. |
| `peak_select` | Calibrants | `highest` | `highest`: tallest accepted drift peak; `last`: latest accepted drift bin. |
| `min_im_intensity` | Calibrants | `500` | Minimum peak height of the summed IM drift profile, in reader units. |
| `peak_height_frac` | Calibrants | `0.10` | Detection threshold as a fraction of the profile maximum; normally 0–1. |
| `gas_mw` | Calibrants | `28.014` | Drift-gas molecular mass in Da; nitrogen default, identical across references. |

Omitted or blank settings use the listed defaults. The program uses the supplied
`mw`, `z`, and `mz` rather than
deriving `mz` from a name or mass. It reads
intensities/drift peaks and RAW metadata, fits known references, and calculates
analyte CCS. Calibrant-only thresholds and peak selection do not control analytes.
Use `direct` for power-law CCS conversion; `twostep` uses the diagnostic linear fit.

The GUI provides **Load CSV**, direct cell edits, **+ Row**, **- Row**, horizontal
scrolling, resizable columns, header/cell tooltips, and **Save CSV** in full-width
**Calibrant sheet** and **Analyte sheet** tabs.
A processing job snapshots edits without rewriting
the source table; use **Save CSV** to a separate file to retain your study settings.

## Brief CLI

Run these commands from this toolkit folder:

```bash
python waters_convert.py examples/production/protein_adp/raw -o output/converted --fn 1 --ms-bin 0 --im-bin 0
python -m deconvovo.imms_plot -i output/converted -o output/html --raw-dir examples/production/protein_adp/raw --config examples/production/protein_adp/viewer_config.json -j 1
python -m deconvovo.imms_ccs_calibrate -o output/ccs --data-dir output/converted --calibrant-csv examples/production/protein_adp/calibrants.csv --analyte-csv examples/production/protein_adp/analytes.csv --conversion-method direct
```

For the separate Betaine visual detail:

```bash
python waters_convert.py examples/production/small_molecule/raw/20260320_BETAINE_V2.raw -o output/betaine_converted --fn 1 --ms-bin 0 --im-bin 0
python -m deconvovo.imms_plot -i output/betaine_converted -o output/betaine_1000_2000 --raw-dir examples/production/small_molecule/raw --config examples/production/small_molecule/betaine_1000_2000_viewer_config.json -j 1
```

Function `1` selects `_FUNC001`; confirm its acquisition role. Zero mass-bin widths
retain reader resolution. The exported tables combine scans; preserve the RAW data.
The converter and `txt_to_csv.py` use only Python's standard library. Scientific
plotting/CCS uses `requirements-pipeline.txt`; use `--help` for individual options.

See [GUI details](GUI_README.md), [pipeline commands](PIPELINE_README.md), and the
complete [viewer data flow and math](docs/html_viewer_data_flow_and_math.md).
[Third-party reader licenses](bin/) remain included. Optional mzML/Rainbow routes
have separate requirements and are unnecessary for the Console.
