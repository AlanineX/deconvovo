# Waters workflow recovery and staleness map

The current examples use protein standards followed by ADP CCS, with tune mix and Betaine as a separate small-molecule check. [PIPELINE_README.md](PIPELINE_README.md) gives runnable commands. The repository table and Carter notes below retain the 2026-09-29 recovery audit; current acquisition provenance is in [demo_manifest.json](examples/production/demo_manifest.json).

## Historical repository snapshot: 2026-09-29

| Repository/branch | Live head | Commit date | Assessment |
| --- | --- | --- | --- |
| [Public DeconVoVo main](https://github.com/AlanineX/deconvovo) | `430475d` | 2026-06-24 | Contains the enhanced Waters viewer and modular CCS pipeline; the corresponding local Waters modules match this head |
| Public DeconVoVo dev | `5be539e` | 2026-03-25 | Older HDX/development snapshot; behind the enhanced viewer |
| [Private DeconVoVo main](https://github.com/AlanineX/deconvovo-dev) | `1ff67d0` | 2026-03-25 | Older viewer/configuration; its core CCS modules match the public/current copies |
| Private DeconVoVo dev | `81230d7` | 2026-03-25 | Older development/documentation snapshot |
| Local main | `6c45ceb` | 2026-07-20 | Seventeen commits ahead of public main, principally HDX work; additional modified/untracked work includes the newer Thermo suite |

At the audit date, live `git ls-remote` matched the cached remote heads. Repository listing identified these two DeconVoVo repositories; no additional Waters-named repository was found in that account/local script search. This table is a dated snapshot, not a current release-state check.

The audit recovered the enhanced Waters modules committed in June. The local README/package description also includes Thermo/native-MS work, while the portable project keeps the Waters workflow independent.

## Working routes, copied entry points and dependencies

```text
Waters RAW folder → waters_convert.py → paired MS/IM TXT
                                      ├→ imms_plot + imms_html → 2D HTML / point CSV
                                      ├→ target-window extraction → mobiligram CSV / PNG
                                      └→ calibration/analyte CSV + RAW metadata → CCS fit / profiles / plots
```

| Capability | Original source in main checkout | Portable entry point | Status |
| --- | --- | --- | --- |
| RAW → text | `deconvovo/waters_convert.py` | `waters_convert.py` | Already copied; standalone reader discovery, standard-library summary, explicit no MS smoothing, failure exit status |
| Enhanced 2D HTML | `deconvovo/imms_plot.py`, `imms_html.py` | `python -m deconvovo.imms_plot` | Current public code copied; one additional patch embeds matching Plotly JS for offline use |
| Specific mass-window mobiligram | `deconvovo/ccs_peak_pick.py::extract_drift_profile` | `extract_mobiligram.py` | Existing extraction function reused in a small new CSV/PNG command |
| Calibrant drift extraction/peak picking | `ccs_peak_pick.py`, `imms_ccs_calibrate.py` | `python -m deconvovo.imms_ccs_calibrate` | Original numerical routines copied unchanged |
| CCS fit/application/plots | `ccs_fit.py`, `ccs_convert.py`, `ccs_plot.py` | Same calibration entry point, optional `--analyte-csv` | Original numerical routines copied unchanged |
| Shared input/units/smoothing | `io.py`, `vendors/waters.py`, `smooth.py`, `parallel.py` | Imported by the copied modules | Included; package initializer selects the headless plotting backend |

Plotting/calibration needs NumPy, pandas, SciPy, Plotly and matplotlib. Wine is only needed for the Windows RAW reader on non-Windows hosts. The original GUI and UniDec mass-deconvolution engine are not dependencies of these plotting/CCS routes.

## Current protein standards → ADP workflow

| Item | Included location / role |
| --- | --- |
| RAW data | `examples/production/protein_adp/raw/`; three protein standards and two ADP solution conditions |
| Calibrant table | `examples/production/protein_adp/calibrants.csv`; 25 inherited protein-reference rows |
| Analyte table | `examples/production/protein_adp/analytes.csv`; 34 inherited ADP target rows |
| Batch overview | `examples/production/protein_adp/viewer_config.json` |
| Separate Betaine 1000–2000 detail | `examples/production/small_molecule/betaine_1000_2000_viewer_config.json`; visual-only view of the existing Betaine RAW, with no CCS tables |
| GUI route | **Load real demo → Protein calibrants + ADP CCS → Process & view**; review tables before enabling CCS |

The current corpus has seven complete acquisitions, 106 original files and 235.75 MiB: the five-run protein/ADP batch plus a two-run tune-mix/Betaine example. The inherited protein fit is poor and ADP extrapolates below its reference interval. Betaine alone supplies one reference and cannot fit the free two-parameter TW-IMS curve; at least three valid calibrant peaks are required. Use the protein standards for ADP. The separate small-molecule example fits six Hexakis references before applying the curve to Betaine.

The supplementary Betaine detail spans 1000–2000 and 4–12 ms. Its high-mass
features remain unassigned, and the Betaine ion near 118 is outside that focus.

## Historical Carter workflow recovered

Base directory: `/home/alan/Downloads/carter_0623/`. The shared custom configuration is `imms_plot_config.json`, copied to `config/carter_imms_plot_config.json`.

| Dataset prefix | RAW acquisitions | Paired converted MS files | Existing enhanced HTML |
| --- | ---: | ---: | ---: |
| `Myoglobin_Water` | 6 | 6 | 3 |
| `Myoglobin_D2O` | 6 | 6 | 3 |
| `Myoglobin_200mMAA_Water` | 6 | 6 | 3 |
| `Myoglobin_200mMAA_D2O` | 6 | 6 | 3 |
| `Myoglobin_200mMTEAA_Water` | 6 | 6 | 3 |
| `Myoglobin_200mMTEAA_D2O` | 6 | 6 | 3 |

For each prefix, converted files live in `<prefix>_converted/` and HTML in `<prefix>_html/`. The 18 existing viewers were modified 2026-06-24 UTC. Their source features match the enhanced generation path. The three-viewer subsets correspond to the 24-hour acquisitions; six acquisitions exist per condition. Some RAW directories also contain older conversion/UniDec sidecars; the dedicated converted folders give the clearest run pairing.

Carter's preset uses mass range `[350,1250]`, drift range `[0,15]`, 6400 mass bins, `W-Spectral`, square-root intensity, Gaussian smoothing, noise cutoff `15`, and fade `50`. These are display choices, not calibration or peak-identification settings.

## What is stale or limited

| Item | Evidence | Recommended use |
| --- | --- | --- |
| Private/public development branches | March heads lack the June enhanced viewer/config changes | Recover enhanced plotting from public main/current local modules |
| `imms_html_blur_diagnosis.md` | Written June 23; describes full-acquisition histogramming. The June 24 code clips the initial heatmap to the configured mass range and updates smoothing/controls | Keep as historical diagnosis; its claim about the current histogram domain is stale |
| `waters_unidec_pipeline.py` reference | Appears in a config comment; no matching current file or tracked history was found | Current modules are `imms_plot` / `imms_html`; the copied config comment is corrected |
| Former all-in-one conversion checks | Earlier cache checks used MS TXT alone and reader failures could be marked OK; current conversion requires the requested nonempty outputs and propagates failures | Use the current converter, batch CLI or merged desktop workflow; select explicit settings for a new acquisition |
| Historical HTML “self-contained” wording | Original HTML explicitly loads `https://cdn.plot.ly/plotly-2.35.2.min.js`; current installed Python Plotly is 6.6.0 with JS 3.4.0 | Copied module embeds the JS matching the installed Python package, fixing the offline dependency/version mismatch |
| Current-view and raw exports | Current-view CSV uses the displayed traces/ranges, including transforms and filters; older browser raw/processed buttons and the generated point CSV have separate meanings | Use the target-window extractor for exact raw mobiligrams; retain original TXT/RAW inputs |
| Drift coordinates and matched ranges | Initial/rebinned heatmaps and numerical exports now use bin centres; linked profiles read the ranges Plotly actually renders | The half-bin time convention and resolved matched-axis ranges keep the plotted axes aligned |
| Historical CCS outputs | Existing tune-mix and protein fits are stored under `output/legacy/`; these are experiment-specific results | Do not apply a curve to unrelated acquisition conditions |

The CCS routines date mainly to March, but they still import and execute. The CSV templates needed portable dataset paths: copied templates now expect matching RAW folders under `data/`; their species/settings are preserved. The calibration orchestrator now rejects differing EDC delay coefficients across calibrant RAW folders; compatible acquisition conditions still require review.

The selected design 01 Instrument Console replaces the five-panel Waters GUI with a three-pane workbench: acquisitions, measured plots, and display/CCS controls. **Process & view** combines conversion and viewing; **CCS review** combines editable calibrant/analyte tables, curve diagnostics, and application results using shared paths. **Advanced processing & defaults → Save defaults** preserves viewer choices across relaunches without dropping configuration fields the desktop does not edit. See [GUI_README.md](GUI_README.md). The original local Thermo workspace remains a separate application workflow.

## Verification and remaining limits

- The existing real RAW input was used to generate a new offline HTML/point CSV and a target-window mobiligram CSV/PNG in `examples/plots/`.
- Exact numerical extraction and half-bin time conversion are checked against the underlying IM points.
- The copied calibration route is exercised with six existing tune-mix calibrants and a Betaine analyte; [production notes/results](examples/production/README.md) include the two complete RAW acquisitions. The full example also passes a fresh install in a relocated project with a fresh Wine prefix. This checks execution, not general accuracy for unrelated acquisitions.
- Browser checks exercise the generated HTML with networking disabled; the validation record lists the controls checked.
- The copied plotting package is exercised from a separate folder to check that imports/config lookup do not depend on the original checkout.
- Release builds check the protein-to-ADP workflow from source and the frozen Windows application, including space-containing paths, offline viewing, alignment and current-view exports. Successful releases attach [native validation](https://github.com/AlanineX/deconvovo/releases/download/v0.3.1/native_smoke.json). Current retained data covers seven acquisitions and 106 original files; [production validation](examples/production/validation.json) retains earlier Linux/Wine execution. macOS native runtime remains untested.
