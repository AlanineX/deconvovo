# Changelog

## v0.3.1 · Waters Instrument Console

- Made protein standards → ADP CCS the primary complete example; seven real acquisitions retain 106 original files and 235.75 MiB.
- Added a separate visual-only Betaine 1000–2000 detail using the existing acquisition, without calibration inputs.
- Aligned heatmap, mass and drift profiles at bin centres and the rendered axis ranges, including the rightmost mass edge.
- Added visible X/Y limits and current-view PNG/CSV exports while retaining offline desktop and portable HTML viewing.
- Added full-width editable calibrant/analyte sheets with field help, column resizing and Save CSV; Convert only needs no reference tables.
- Preserved blank optional formula defaults and shared Windows/UNC/space-containing path handling.
- Extended source/frozen Windows gates to the protein-to-ADP example, path handling, alignment and all current-view exports; retained permanent validation, provenance and preview assets.
- Preserved the complete mathematical explanation in both source layouts.

The protein calibration fit is poor, many ADP results extrapolate, and the inherited assignments require review. The Betaine detail is a visual example with unassigned high-mass features. See [demo notes](demo_data/production/README.md).

## v0.3.0

Introduced the selected design 01 Console, merged processing/viewing, local enhanced viewer, portable offline HTML, saved defaults and native validation assets. Its generated distribution is superseded by v0.3.1.

## v0.2.3

Included the enhanced browser viewer. The later Console corrected omitted settings, saved defaults and the external Plotly dependency.
