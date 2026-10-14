# Waters Instrument Console

The selected Instrument Console design keeps the acquisition browser, measured
plot, and display/CCS inspector visible together. Source and output paths are
shared across processing and analysis:

- **Process & view:** convert complete RAW acquisitions, inspect each run in native plots, and reload existing results. HTML is generated only on explicit export.
- **Convert only:** create MS/IM tables without generating HTML, running CCS, or requiring calibrant/analyte CSVs.
- **CCS review:** edit calibrant and analyte tables, fit and apply a calibration in one job, then inspect curves, reference mobility profiles, and analyte summaries.

RAW extraction → paired MS/IM tables → native plots → optional calibration and analyte CCS.
Portable offline HTML is a separate export from those converted tables.

## Install and launch

Run from this portable project directory. Linux/macOS:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements-gui.txt
.venv/bin/python run_gui.py
```

Windows PowerShell:

```powershell
py -3 -m venv .venv
.venv\Scripts\python -m pip install -r requirements-gui.txt
.venv\Scripts\python run_gui.py
```

Windows x64 runs the bundled reader natively. Linux RAW extraction uses Wine.
macOS also needs compatible x64 Wine/emulation. Python and the declared
dependencies are required for this source distribution.

## Use real acquired data

Choose a profile and click **Load real demo**, then **Process & view**. The
demo's source, display configuration, and available CCS tables load together.
Results are written to a user documents folder; change Output if needed.
The [production demo guide](examples/production/README.md) records the seven
complete acquisitions, 106 original files, 235.75 MiB, checks, and scientific limits.

| Profile | Complete RAW runs | Workflow |
| --- | ---: | --- |
| Protein calibrants + ADP CCS (default) | 5 | Three protein standards fit the curve for two ADP conditions; 25 reference rows and 34 inherited target assignments |
| Tune mix + Betaine | 2 | Batch viewing, six calibrant rows, one analyte |
| Betaine · 1000–2000 m/z detail | One existing Betaine run | Visual-only detail; no calibrant/analyte CSVs or CCS fit |

CCS starts off. Review the tables in **CCS review**, then choose
**Calibrate & apply**, or enable CCS during **Process & view**. Edited table
snapshots feed the job without rewriting the shipped reference tables.
Relative RAW paths resolve beside their original CSV.

The inherited protein calibration has $R^2 = 0.874117$ and needs peak-selection
and condition review. The tune-mix fit has $R^2 = 0.999915$, but Betaine lies
below its reference drift interval. ADP also extrapolates. These demonstrations
test processing and expose QC decisions; they do not establish quantitative CCS
accuracy. Betaine alone supplies one CCS reference and cannot fit the free
two-parameter TW-IMS curve; this implementation requires at least three valid
calibrant peaks. Use the protein standards for the ADP example. The tune-mix
calibration remains a separate small-molecule workflow.

For your data, select a complete `.raw` directory, a parent containing RAW
acquisitions, or an existing folder of paired `_ms.txt`/`_im.txt`. Nested
condition folders are supported. Double-click a completed viewer row to open
the run's native plots. **Open output folder** locates exported artifacts.

## Display and extraction settings

**Advanced processing & defaults** holds acquisition function and bin
widths plus display intensity, smoothing, cutoff/fade, rendering, native/display
bins, profile/stick mode, and initial mass/drift ranges. RAW extraction disables
reader-side MS smoothing; visual smoothing remains adjustable in the viewer.
Use **Save defaults** in this section to keep viewer settings across relaunches.
Saved configuration retains fields that the desktop does not edit; frozen builds
store it in the user's application configuration directory.

The separate **Betaine · 1000–2000 m/z detail** choice uses the same complete
Betaine acquisition with the included
[focused configuration](examples/production/small_molecule/betaine_1000_2000_viewer_config.json).
It spans 1000–2000 and 4–12 ms, with cube-root intensity, 6400 mass bins and
native drift bins. Its high-mass features are unassigned; the Betaine ion near
118 is outside the focus. The protein/ADP overview remains the default CCS example.
The [pipeline guide](PIPELINE_README.md) gives the separate visual-only command.

The visible **X · m/z** and **Y · drift time (ms)** controls provide minimum and
maximum fields. CCS review uses full-width **Calibrant sheet** and **Analyte sheet**
tabs, with separate calibration/profile/result tabs for inspection.

**Reuse unchanged RAW conversions** checks extraction settings and source-file
metadata. Changed functions, bin widths, or input files trigger fresh extraction.
Native plots update from the loaded tables with current display settings.
A failed reader run
stops downstream work and retains previous files without marking the new job
complete.

The central viewer renders the heatmap and linked profiles with Matplotlib and
NumPy/SciPy. It does not embed a browser or load HTML. Inspector controls update
the native plots; widening the X/Y ranges can use the source TXT directly.
**Reset view** restores the numeric ranges while retaining current filters.
Mouse zoom/pan updates both linked profiles; changing a numeric range replaces
that mouse viewport, while changing an unrelated display control retains it.

**Export current view** offers full-plot, heatmap, MS-profile and drift-profile
PNG, plus heatmap/MS/drift CSV. These exports use the current view and settings;
the 2D CSV contains displayed transformed/filtered values, not raw point intensities.

**Export portable HTML** regenerates only the selected converted run with the
current display defaults, then saves one complete HTML file. It embeds both
its data and matching Plotly JavaScript: copy that file alone and open it in a
browser on another computer. It needs no Python, server, companion files, or
internet connection. The full exported toolbar remains available.

Native viewing and separate HTML export work offline after installing the source dependencies, or
directly from the self-contained native application. Extraction and CCS also
run locally using the supplied reader and data.
Exact target-window extraction remains available through `extract_mobiligram.py`.
The desktop job uses one background thread and processes acquisitions in order,
keeping memory bounded. Controls lock while the job runs; progress identifies
the current run and stage. Use the existing CLIs for parallel batch work.

## Verification

The [retained production record](examples/production/validation.json) contains
earlier Linux/Wine conversions and offline viewers for all seven retained
acquisitions, plus both CCS demo profiles. Selecting those stored results did
not rerun the analysis. The two-run small-molecule workflow
also completed through the desktop's real background worker. Chrome checks with
network access disabled exercised bin controls, profile/stick mode, and linked
heatmap zoom.

Release builds run the protein-standards → ADP workflow from source and the
frozen native application using paths containing spaces. Gates check alignment,
all seven current-view exports, native fonts, source protection and separate
portable offline HTML. The builder rejects Chromium/WebEngine artifacts.
Successful releases attach [native validation](https://github.com/AlanineX/deconvovo/releases/download/v0.3.3/native_smoke.json),
[build provenance](https://github.com/AlanineX/deconvovo/releases/download/v0.3.3/build_provenance.json), and an
[actual GUI preview](https://github.com/AlanineX/deconvovo/releases/download/v0.3.3/native_console.png). These software checks do not
establish quantitative CCS accuracy. macOS native runtime remains untested.
