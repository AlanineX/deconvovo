# Protein standards to ADP CCS

Run these commands from the project directory after installing the source dependencies. Windows runs the bundled reader natively; Linux RAW extraction uses Wine.

```sh
python -m deconvovo.waters_convert demo_data/production/protein_adp/raw -o output/adp/_converted --fn 1 --ms-bin 0 --im-bin 0
python -m deconvovo.imms_ccs_calibrate -o output/adp/ccs_direct --data-dir output/adp/_converted --calibrant-csv demo_data/production/protein_adp/calibrants.csv --analyte-csv demo_data/production/protein_adp/analytes.csv --conversion-method direct
```

In the GUI, select **Protein standards + ADP**, choose **Process files**, then open **Calibration**. Process references, review the numbered peaks, choose **Fit calibration**, then **Generate CCS**. Conversion and spectrum viewing work without CCS sheets.

The inherited protein reference fit and ADP extrapolation require review before quantitative use. See the [GUI and CSV guide](README.md) and [demo details](demo_data/production/README.md).
