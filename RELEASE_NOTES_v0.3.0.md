# Waters Instrument Console · v0.3.0

Download **WatersConsole-win64.zip**, extract the whole archive, and open **WatersConsole/WatersConsole.exe**. Python, Wine, UniDec and an internet connection are unnecessary for the Windows application.

The selected design 01 puts the acquisition list, large plots and settings together. **Process & view** converts a batch and opens its results; optional calibration and analyte application share paths and method settings.

Choose either the **local desktop viewer** or **Export portable HTML**. Both have the enhanced controls and work offline. Share one exported HTML file; it contains the plotting library and the selected run data.

Nine complete real acquisitions are included: five protein/ADP runs, two tune-mix/Betaine runs and two native HSP104 runs, totaling **300.29 MiB of original RAW data**. Viewer profiles and optional reference tables are supplied. The protein calibration fit is poor, ADP/Betaine results extrapolate beyond reference drift ranges, and HSP104 has no supplied compatible CCS calibration.

The new build runs focused source tests and verifies all 148 original demo files. Before and after freezing, it converts the tune-mix and Betaine RAW acquisitions, fits six CCS references, applies the curve to one analyte, and checks ten result PNGs. It also checks offline HTML, the embedded console, native fonts and laptop layout. The source and frozen Windows application passed these checks, including both offline viewing routes. Permanent [native validation](https://github.com/AlanineX/deconvovo/releases/download/v0.3.0/native_smoke.json), [build provenance](https://github.com/AlanineX/deconvovo/releases/download/v0.3.0/build_provenance.json), and an [actual Windows GUI preview](https://github.com/AlanineX/deconvovo/releases/download/v0.3.0/native_console.png) accompany the ZIP. All nine acquisitions have a separate Linux/Wine conversion record; the Windows smoke covers the two small-molecule runs. Existing command-line conversion, CCS and optional UniDec deconvolution interfaces remain available.

Keep the executable and `_internal` folder together. Outputs go to the writable folder you choose. Portable source installation instructions and scientific limitations are in the README and demo notes.
