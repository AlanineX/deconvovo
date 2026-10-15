# Waters Instrument Console

Open Waters RAW, inspect native plots, and optionally calculate CCS. Everything runs locally.

[Downloads](https://github.com/AlanineX/deconvovo/releases): Windows—extract the whole ZIP and open `WatersConsole.exe`. Linux—extract the archive and run `WatersConsole` (Ubuntu 24.04 x64; install Wine for RAW conversion).

## GUI

1. Choose **Input** (RAW folder/parent) and **Output** in the two wide rows, then **Process files**. Valid unchanged conversions are reused.
2. **Targets → Add:** enter m/z min/max and select the row to focus the 2D plot and drift profile. Keep **Auto** or enter manual drift bounds.
3. Hover, wheel-zoom or use **Zoom/Pan**. Export the interest profile as PNG/CSV; **Export portable HTML** creates one offline interactive file.
4. **Calibration:** **Process references** → choose peaks → **Fit calibration** → inspect → **Generate CCS**. Complete targets then update automatically.

Optional **Examples**: **Protein standards + ADP**, **Tune mix + Betaine**, **Betaine**. Seven real acquisitions, 106 files, 235.75 MiB. [Demo limits](demo_data/production/README.md): poor protein fit; ADP/Betaine extrapolate.

## CSV inputs

| Field | Meaning |
| --- | --- |
| `name` / `species` | Calibrant label / analyte label. |
| `mw` | Neutral molecular mass, Da. |
| `z` | Positive integer charge. |
| `mz_min,mz_max` / `mz` | Inclusive bounds / legacy centre, Th. |
| `CCS_reported` | Known reference CCS, Å²; calibrants only. |
| `raw_path` | Complete `.raw` folder; relative to the CSV. |

```csv
name,mw,z,mz_min,mz_max,CCS_reported,raw_path
UBQ z=8,8565,8,1069.63,1073.63,1622,raw/ubiquitin.raw
```

```csv
species,mw,z,mz,raw_path
[ADP+H]+,427.2,1,428.21,raw/adp_edda.raw
```

A fit needs at least three valid reference peaks. [Full CSV reference](docs/csv_reference.md) covers optional fields, thresholds and Windows paths.

## CLI — replace paths as needed

```bash
python -m deconvovo.waters_convert raw -o converted
python -m deconvovo.imms_plot -i converted -o html --raw-dir raw
python -m deconvovo.imms_ccs_calibrate -o ccs --data-dir converted --calibrant-csv calibrants.csv --analyte-csv analytes.csv
```

[Source setup and GUI details](portable_waters_convert/GUI_README.md) · [Viewer math](docs/html_viewer_data_flow_and_math.md) · [Licenses](portable_waters_convert/bin/)
