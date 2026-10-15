# Offline Waters CLI

Run these commands from this folder with Python 3.10+. Install
`requirements-pipeline.txt` first. Linux RAW conversion requires Wine64;
Windows uses the bundled reader. Processing and exported viewers work offline.

## Protein standards and ADP

```bash
python waters_convert.py examples/production/protein_adp/raw -o converted --fn 1 --ms-bin 0 --im-bin 0
python -m deconvovo.imms_ccs_calibrate -o ccs --data-dir converted --calibrant-csv examples/production/protein_adp/calibrants.csv --analyte-csv examples/production/protein_adp/analytes.csv
python -m deconvovo.imms_plot -i converted -o html --raw-dir examples/production/protein_adp/raw --config examples/production/protein_adp/viewer_config.json
```

Calibration is optional: omit the CCS command to convert or inspect data only.
Review reference peaks, acquisition conditions and extrapolation warnings before
interpreting CCS. The inherited protein assignments need review.

## Betaine detail

```bash
python waters_convert.py examples/production/small_molecule/raw/betaine.raw -o betaine_converted --fn 1 --ms-bin 0 --im-bin 0
python -m deconvovo.imms_plot -i betaine_converted -o betaine_html --raw-dir examples/production/small_molecule/raw/betaine.raw --config examples/production/small_molecule/detail_view.json
```

This displays unassigned high-mass features at 1000–2000 m/z. The Betaine ion near
118 m/z is outside that detail view.

Each HTML output embeds its data and Plotly JavaScript. Copy the HTML alone to
view it offline. The native GUI uses Qt plots; HTML is a separate export.
See the [GUI guide](GUI_README.md), [CSV columns](docs/csv_reference.md), and
[complete real demos](examples/production/README.md).
