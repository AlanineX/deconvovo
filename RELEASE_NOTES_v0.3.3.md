# Waters Instrument Console · v0.3.3

Desktop plots now use **native Qt/Matplotlib with NumPy/SciPy**. The application contains no Chromium/WebEngine runtime and does not generate or interpret HTML for ordinary viewing. Conversion, native viewing and CCS run offline.

**Export portable HTML** explicitly creates a separate, complete offline Plotly document with the selected run data and current defaults. Native PNG/CSV exports remain available directly from the GUI; a browser is needed only when opening the optional HTML document.

The selected design 01 Console retains visible X/Y limits, Reset view, seven current-view exports, full-width editable calibrant/analyte sheets and Convert only. Native plots can widen directly from the loaded TXT data. The full viewer math describes the native API, bin centres, smoothing, transforms, cutoff/fade, linked profiles and export precision.

The primary real example uses three protein standards to calibrate two ADP acquisitions. The package retains **seven complete acquired RAW folders, 106 original files, 235.75 MiB**. The separate **Betaine · 1000–2000 m/z detail** is visual-only; its high-mass features are unassigned. Poor protein fit and ADP extrapolation require review before quantitative interpretation.

Release gates run every focused test, real protein-to-ADP processing from source and the frozen application, quoted/space-containing paths, native profile alignment, all seven PNG/CSV exports, protected source files, and separate offline HTML generation. The builder rejects Chromium/WebEngine artifacts. Permanent [native validation](https://github.com/AlanineX/deconvovo/releases/download/v0.3.3/native_smoke.json), [build provenance](https://github.com/AlanineX/deconvovo/releases/download/v0.3.3/build_provenance.json), and an [actual GUI preview](https://github.com/AlanineX/deconvovo/releases/download/v0.3.3/native_console.png) accompany a successful release.

Download **WatersConsole-win64.zip**, extract the whole archive, and open **WatersConsole/WatersConsole.exe**. Python, Wine and internet access are unnecessary for the compiled Windows application.
