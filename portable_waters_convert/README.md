# Portable Waters toolkit

Copy this complete folder to use the Waters converter, enhanced viewer and Instrument Console independently of the rest of DeconVoVo. On Windows, the [compiled application](https://github.com/AlanineX/deconvovo/releases/tag/v0.3.0) bundles Python and its dependencies; extract it and open `WatersConsole.exe`.

## Desktop from source

From this folder, with Python 3.10+:

```bash
python -m venv .venv
# Activate .venv for your operating system.
python -m pip install -r requirements-gui.txt
python run_gui.py
```

Select a real demo or your RAW folder, choose a writable output folder, and click **Process & view**. Select a processed run to inspect it locally, open its offline HTML, or export one portable HTML file. The three demo profiles supply the appropriate viewer settings and optional CCS tables. [Nine complete real acquisitions](examples/production/README.md) contain 300.29 MiB of original RAW data.

The chosen Instrument Console keeps the run list, plots and controls together. CCS calibration and analyte application share paths and method settings. Both the embedded desktop viewer and exported HTML use bundled plotting assets and run without internet access. **Advanced → Save defaults** keeps viewer settings across relaunches.

## Converter alone

The converter and TXT-to-CSV scripts use only Python's standard library. The included `bin/CDCReader.exe`, `MassLynxRaw.dll` and `cdt.dll` run natively on Windows x64. Linux/macOS RAW conversion needs compatible x64 Wine; scientific plotting and CSV processing run in Python.

```bash
python waters_convert.py examples/production/small_molecule/raw -o output/converted --fn 1 --ms-bin 0 --im-bin 0
python txt_to_csv.py output/converted
```

The outputs are paired `*_ms.txt` and `*_im.txt` tables. The IM table stores mass-to-charge ratio, native drift-bin index and intensity; per-run pusher timing converts the drift axis to milliseconds. The tables combine scans; preserve every file in each original RAW folder.

## Scientific tools

```bash
python -m deconvovo.imms_plot -i output/converted -o output/html --raw-dir examples/production/small_molecule/raw --config examples/production/small_molecule/viewer_config.json
python -m deconvovo.imms_ccs_calibrate -o output/ccs --data-dir output/converted --calibrant-csv examples/production/small_molecule/calibrants.csv --analyte-csv examples/production/small_molecule/analytes.csv
```

`extract_mobiligram.py` extracts a chosen mass window through the same numerical modules. Use `--help` for its options. The optional mzML/Rainbow scripts have their own requirements files and are not needed for the Console.

The protein demo has a poor calibration fit and many ADP results extrapolate below the reference interval; Betaine also extrapolates despite a high tune-mix fit score. HSP104 has no compatible supplied CCS curve. Read [demo scientific limits](examples/production/README.md) before interpreting results.

## Native build

```bash
python -m pip install -r build/requirements.txt
python build/build_native.py --output-dir build/artifacts
```

Build on Windows x64 for a Windows application. The builder packages the reader, offline Plotly, Qt WebEngine and nine full demos, verifies their original bytes, and records provenance. [Third-party licenses](bin/) are included; UniDec requests citation of [Marty et al., Analytical Chemistry (2015)](https://doi.org/10.1021/acs.analchem.5b00140) for publications using its software.
