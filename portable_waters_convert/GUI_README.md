# Waters Instrument Console

The Instrument Console opens an acquisition, displays native plots, and keeps
its editable targets beside the spectrum. Calibration is optional and reusable:

- Set the two wide **Input/Output** rows, refresh the acquisition list, then **Process files**.
- **Targets:** enter an m/z interval; selecting its row focuses the 2D plot and exact drift profile.
- **Calibration** tab: process references, choose peaks, fit and review; saved curves can be loaded.

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

Your own data needs no named setup: choose Input/Output and **Process files**.
The viewer uses saved display styles and the acquisition's full range, with CCS off. Outputs go to
`Documents/DeconVoVo/Waters/<acquisition>`; change Output if needed.
Wide editable **Input** and **Output** rows span the panels; **Data** holds only **Examples**.
Choosing an example fills both paths, its display settings and available CCS
references/targets, and lists its RAW folders. It does not convert them yet.
Set Input/Output directly or use their **…** folder pickers. **Acquisitions**
lists filenames and **State**: **RAW** awaits processing; **Ready** has converted
files. Both columns resize. **Refresh** rescans Input for added/deleted RAW
folders and removes acquisitions whose source disappeared. **Process files**
processes the listed runs; select one acquisition for native viewing.
The integer **CPU** field limits simultaneous RAW readers; it defaults to one
and permits up to the available CPU cores.
The [production demo guide](examples/production/README.md) records the seven
complete acquisitions, 106 original files, 235.75 MiB, checks, and scientific limits.

| Profile | Complete RAW runs | Workflow |
| --- | ---: | --- |
| Protein standards + ADP (default) | 5 | Three protein standards fit the curve for two ADP conditions; 25 reference rows and 34 inherited target assignments |
| Tune mix + Betaine | 2 | Two acquisitions, six reference ions, one analyte |
| Betaine | One existing Betaine run | Visual-only detail; no calibrant/analyte CSVs or CCS fit |

Opening an example does not calculate CCS. In **Calibration**, process the
reference table, inspect its peaks, then fit/review before applying.
Edited table snapshots feed the job without rewriting the shipped CSV.
Relative reference RAW paths resolve beside their original CSV.

The inherited protein calibration has $R^2 = 0.874117$ and needs peak-selection
and condition review. The tune-mix fit has $R^2 = 0.999915$, but Betaine lies
below its reference drift interval. ADP also extrapolates. These demonstrations
test processing and expose QC decisions; they do not establish quantitative CCS
accuracy. Betaine alone supplies one CCS reference and cannot fit the free
two-parameter TW-IMS curve; this implementation requires at least three valid
calibrant peaks. Use the protein standards for the ADP example. The tune-mix
calibration remains a separate small-molecule workflow.

**Input** accepts complete `.raw` directories and parents containing runs.
Validated converted pairs are reused automatically; no TXT picker is needed.
Nested acquisition folders are supported. Select a run to view it.
**Open output folder** locates exported artifacts.

## Interest profiles and optional CCS

Under **Spectrum**, **Targets → Add** starts a row from the visible m/z interval.
**m/z min/max** are required, finite and increasing. Selecting a row focuses the
2D plot to its mass interval and drift bounds, then shows its exact profile
beside the heatmap and in **Mobiligram**. Auto bounds fill all valid rows before
selection, on data load/import or mass-bound edits. **Auto** encloses every positive
raw drift sample inside the interval, ignoring display filters, with equal 20%
padding before clipping to the acquisition. For a 5–11 ms signal this gives
2.8–13.2 ms before clipping. Uncheck **Auto** to enter finite increasing
drift bounds manually; their header tooltips give the units. **Full view** deselects the target and restores the
acquisition; **Home/Reset view** restores full acquisition axes, including zero
endpoints, while retaining display filters. ADP EDDA's full drift extent is
13.8 ms (200 bins, 69 µs pusher), rather than a preset's viewing limit.
Mass/charge fields stay hidden for profiling and appear for a CCS model. **Export → Drift PNG/CSV** saves
the selected profile inside the current drift view. **Import CSV** loads targets;
**Export CSV** saves target inputs, drift ranges, computed CCS results and flags.
Import restores results only for the same acquisition; otherwise they become
stale while inputs remain. These are separate from the exported profile samples. See the
[CSV reference](docs/csv_reference.md).

The central **Calibration** tab has linked plots above one editable reference
table; drag the vertical divider to adjust their heights. Use **Load CSV**, or
**Add standard RAW → + Row**. Enter explicit **m/z min/max**, mass, charge and
**CCS_reported**; legacy centre/window CSVs also load. **Extraction options**
exposes signal thresholds. **Process references** detects candidates without
fitting; watch processing status and inspect their native drift times.
Click an ion row to show its profile and circle its calibration point; clicking
the curve point selects that row. **Highest/Last** buttons immediately choose
the tallest/latest accepted peaks for **Selected ions** (Ctrl/Shift multiselect)
or **All ions**. For **Manual**, click a detected profile marker; no automatic
fallback is used. Chosen peak/time cells are read-only with full-precision tooltips.
Peak choices preview immediately; the previous curve becomes dashed and requires
refitting. **Fit calibration** is explicit. Inspect the new curve and warnings,
then click **Generate CCS**. The selected complete target gains an interactive
CCS profile below Axis ranges on the right. **Save calibration** stores the
current curve; **Load calibration** restores it without refitting. Loaded curves
also require **Generate CCS** before target processing. Complete targets then
refresh automatically when their inputs or acquisition change. Peak or fit edits
clear old CCS until refitting and generating again. Targets need mass and charge;
incomplete rows remain uncalculated. Saving preserves results.
Reprocess after extraction edits and reuse curves only under matching conditions.

## Display and viewing

The separate **Betaine** choice uses the same complete
Betaine acquisition with the included
[focused configuration](examples/production/small_molecule/detail_view.json).
It spans 1000–2000 and 4–12 ms, with cube-root intensity, 6400 mass bins and
native drift bins. Its high-mass features are unassigned; the Betaine ion near
118 is outside the focus. The protein/ADP overview remains the default CCS example.

The visible **X · m/z** and **Y · drift time (ms)** controls provide minimum and
maximum fields. They follow the current zoom; editing them sets the view range.
They are separate from a selected target's extraction interval.

Conversion reuse checks extraction settings and source-file metadata. Changed functions, bin widths, or input files trigger fresh extraction.
Native plots update from the loaded tables with current display settings.
A failed reader run
stops downstream work and retains previous files without marking the new job
complete.

The central viewer renders the heatmap and linked profiles with Matplotlib and
NumPy/SciPy. It does not embed a browser or load HTML. Inspector controls update
the native plots; widening the X/Y ranges can use the source TXT directly.
**Home/Reset view** restores the full acquired mass and physical drift ranges
while retaining current filters.
Drift-only zoom crops the native curve and rescales its intensity axis. All
acquired samples remain on the line at every zoom, with no automatic point markers.
Changing the extraction mass window recomputes the drift profile; selected targets
retain their explicit mass bounds. The heatmap's Gaussian smoothing and
cube-root/other intensity transforms affect only its 2D appearance.
For very narrow mass views, native Gaussian smoothing uses surrounding acquired
signal as a convolution margin. Axis bounds and target 1D extraction stay exact;
the backing heatmap grid can extend beyond the visible region.
Changing a numeric range replaces the mouse viewport; unrelated display
controls retain it.
Hover reports coordinates and the processed heatmap cell value in the toolbar.
Wheel zoom needs no mode selection; **Zoom** drags a rectangle and **Pan** moves
the view. **Back/Forward** restore previous bounds; **Home** restores the full acquisition.

**Export** offers full-plot, heatmap, MS-profile and drift-profile
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
A desktop worker coordinates separate reader processes, up to **CPU**, with
isolated staged output pairs. HTML and optional CCS stages run sequentially
after all conversions succeed. Controls lock during processing and progress
identifies each run/stage. Failed refreshes invalidate their conversion cache;
queued readers cancel, running readers finish, and downstream work stops.
Earlier files can remain available without marking a failed refresh successful.

## Verification

The [retained production record](examples/production/validation.json) contains
earlier Linux/Wine conversions and offline viewers for all seven retained
acquisitions, plus both CCS demo profiles. Selecting those stored results did
not rerun the analysis. The two-run small-molecule workflow
also completed through the desktop's real background worker. Chrome checks with
network access disabled exercised bin controls, profile/stick mode, and linked
heatmap zoom.

Earlier Windows x64 source and frozen application checks passed with two complete
tune-mix/Betaine RAW acquisitions, six calibrant peaks, one finite analyte CCS
result, and ten PNG outputs. The embedded viewer and a copied single-file HTML
rendered offline, and the console fitted a 1280 × 800 display with native Segoe UI.
Current retained-data verification covers seven RAW
acquisitions and 106 original files; new package checks are recorded separately.

Older `gui_validation.json` and `self_contained_validation.json` describe the
previous five-panel copy; [production validation](examples/production/validation.json)
retains the corresponding seven-run Linux/Wine execution results. macOS native runtime remains
untested.
