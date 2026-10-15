# Changelog

## v0.3.4 · Native acquisition and calibration workspace

- Refined both themes with native glass materials, readable selection/focus states and a visible dark gradient.
- Simplified RAW processing, conversion-cache reuse, resizable acquisition columns and complete real demos with simple names.
- Preserved raw 1D points through zoom and 2D bin changes; added target editing and linked native CCS profiles.
- Merged calibration review and reference editing, added bulk/manual peak choices, cached fitting and explicit CCS generation with stale-result handling.
- Validated source and frozen applications with complete real acquisitions on Windows and Linux; kept portable offline HTML separate.

## v0.3.3 · Native Waters Instrument Console

- Replaced the desktop browser engine with native Qt/Matplotlib plots backed by NumPy/SciPy; Chromium/WebEngine is excluded from the application.
- Separated portable offline HTML generation from ordinary native viewing.
- Preserved all display controls, seven live PNG/CSV exports, editable sheets, path handling and complete real data.
- Documented the native math API, serialization precision, SG normalization and physical drift smoothing widths.
- Verified intensity labels using internal mathtext, including after CCS plotting changes global font settings.

## v0.3.2 · Waters Instrument Console

- Fixed UTF-8 decoding in Windows alignment-test fixtures while retaining every assertion and native gate.
- Removed confirmed unused Qt wrapper/debug resources from the native build while retaining all seven complete acquisitions and required browser/scientific runtimes.
- Made protein standards → ADP CCS the primary complete example; seven real acquisitions retain 106 original files and 235.75 MiB.
- Added a separate visual-only Betaine 1000–2000 detail using the existing acquisition, without calibration inputs.
- Aligned heatmap, mass and drift profiles at bin centres and the rendered axis ranges, including the rightmost mass edge.
- Added visible X/Y limits and current-view PNG/CSV exports while retaining offline desktop and portable HTML viewing.
- Added full-width editable calibrant/analyte sheets with field help, column resizing and Save CSV; Convert only needs no reference tables.
- Preserved blank optional formula defaults and shared Windows/UNC/space-containing path handling.
- Extended source/frozen Windows gates to the protein-to-ADP example, path handling, alignment and all current-view exports; retained permanent validation, provenance and preview assets.
- Preserved the complete mathematical explanation in both source layouts.

The protein calibration fit is poor, many ADP results extrapolate, and the inherited assignments require review. The Betaine detail is a visual example with unassigned high-mass features. See [demo notes](demo_data/production/README.md).

## v0.3.1 · Unreleased candidate

Windows validation stopped at an implicit cp1252 test-fixture read of UTF-8 HTML: 117 tests passed and six setups failed. No release binary was published. The corrected fixture is included in v0.3.2; the failed tag remains in history.

## v0.3.0

Introduced the selected design 01 Console, merged processing/viewing, local enhanced viewer, portable offline HTML, saved defaults and native validation assets. Its generated distribution is superseded by v0.3.2.

## v0.2.3

Included the enhanced browser viewer. The later Console corrected omitted settings, saved defaults and the external Plotly dependency.
