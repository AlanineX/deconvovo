# Changelog

## v0.3.0 · Waters Instrument Console

- Replaced the desktop's separate step pages with the selected design 01 workbench: run list, measured plots, shared settings and an optional activity log.
- Merged conversion and viewing into **Process & view**; calibration and analyte application share the same paths and method. Background work keeps the console responsive and locks affected inputs while a job runs.
- Added a local enhanced viewer and **Export portable HTML**. Both include local Plotly assets and work offline; export preserves the chosen settings.
- Preserved saved viewer defaults, exposed the enhanced controls and **Advanced → Save defaults**, and refreshed HTML when processing so older output cannot mask changed settings.
- Kept terminal progress compatible with standard Windows text encoding, including per-run pusher timing.
- Added nine complete real RAW acquisitions across protein/ADP, tune-mix/Betaine and native HSP104 workflows, with original bytes and portable configuration.
- Added a pinned Waters-only Windows build, byte verification for all 148 demo files, and source/frozen checks covering two RAW acquisitions, six CCS references, one analyte and ten result PNGs, plus offline rendering, native fonts and laptop layout.
- Added permanent native validation, build provenance and actual Windows GUI preview assets beside the application ZIP.
- Conversion errors fail the workflow instead of reporting successful empty results. Existing conversion, CCS and optional UniDec command-line interfaces remain available.

The demonstrations retain their scientific limitations: poor protein calibration fit, ADP/Betaine extrapolation, inherited candidate assignments, and no compatible HSP104 CCS curve. See [demo notes](demo_data/production/README.md).

## v0.2.3

The enhanced browser viewer was present, but the desktop omitted several settings, generation could discard saved defaults, and HTML loaded Plotly from the internet. Its Windows package did not include the new Console, offline desktop viewer or complete production demos.
