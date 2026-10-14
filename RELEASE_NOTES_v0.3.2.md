# Waters Instrument Console · v0.3.2

Download **WatersConsole-win64.zip**, extract the whole archive, and open **WatersConsole/WatersConsole.exe**. Conversion, the desktop viewer, and single-file HTML work offline.

The primary example uses **three protein standards to calibrate two ADP acquisitions**. The separate tune-mix/Betaine example remains an auxiliary check. The package contains **seven complete acquired RAW folders, 106 original files, 235.75 MiB**. The protein fit is poor and many ADP peaks extrapolate; review reference assignments and acquisition conditions before interpreting CCS.

The separate **Betaine · 1000–2000 m/z detail** view shows its existing acquisition over 4–12 ms. Those high-mass features are unassigned, and the Betaine ion near 118 is outside the focus. This visual-only view supplies no calibration tables; Betaine alone cannot fit the free two-parameter CCS curve.

The Windows alignment-test fixture now explicitly decodes UTF-8 HTML. The earlier unreleased candidate stopped at 117 passed tests and six fixture setup errors; no binary was published, and all alignment assertions remain enforced.

The native build omits unused Qt wrapper/debug resources while retaining the linked libraries, normal WebEngine resources, scientific tools and all complete demo acquisitions.

The update aligns drift and mass profiles with the heatmap, includes its rightmost mass edge, and uses the displayed axis ranges. Visible X/Y range fields control the local view. **Export current view** offers four PNG and three CSV choices; **Export portable HTML** keeps a complete offline interactive document.

**CCS review** provides full-width editable calibrant/analyte sheets, field help, resizable columns and Save CSV. **Convert only** extracts tables without requiring reference CSVs. Blank optional formula cells use the existing mass-based estimate. Shared path handling supports Windows drive/UNC paths and paths containing spaces.

The release build runs all eight focused test files, verifies the complete demo bytes, and executes the protein-to-ADP workflow from source and from the frozen application under paths containing spaces. Native gates also check offline viewing, aligned profiles, live exports, native fonts and the actual laptop layout. Permanent [native validation](https://github.com/AlanineX/deconvovo/releases/download/v0.3.2/native_smoke.json), [build provenance](https://github.com/AlanineX/deconvovo/releases/download/v0.3.2/build_provenance.json), and an [actual Windows GUI preview](https://github.com/AlanineX/deconvovo/releases/download/v0.3.2/native_console.png) accompany a successful release.

The full viewer math remains in the source documentation. Existing conversion, calibration and optional UniDec CLI interfaces remain available.
