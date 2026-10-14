# Protein standards to ADP CCS

Run these commands from the project directory after installing the source dependencies. Windows runs the bundled reader natively; Linux RAW extraction uses Wine.

```sh
python -m deconvovo.waters_convert demo_data/production/protein_adp/raw -o output/adp/_converted --fn 1 --ms-bin 0 --im-bin 0
python -m deconvovo.imms_ccs_calibrate -o output/adp/ccs_direct --data-dir output/adp/_converted --calibrant-csv demo_data/production/protein_adp/calibrants.csv --analyte-csv demo_data/production/protein_adp/analytes.csv --conversion-method direct
```

The GUI offers the same example under **Protein calibrants + ADP CCS**. Select **Include CCS when processing** to calculate CCS; leave it off for conversion and viewing. **Convert only** extracts the MS/IM tables without generating viewers or running calibration.

The inherited protein reference fit and ADP extrapolation require review before quantitative use. See the [GUI and CSV guide](README.md) and [demo details](demo_data/production/README.md).
