# Waters CSV reference

Use UTF-8, comma-separated files with a header. Each calibrant row describes
one reference ion/charge state with known CCS; each analyte row describes one
extraction target. Start with the bundled
[protein calibrants](../examples/production/protein_adp/calibrants.csv) and
[ADP analytes](../examples/production/protein_adp/analytes.csv) for the ADP example.

These sheets are for CCS. The separate interest m/z → drift-profile workflow
requires only the m/z range, without molecular mass, charge or a CCS sheet.

## CCS required columns

| Column | Calibrants | Analytes | Meaning |
| --- | --- | --- | --- |
| `name` | Yes | — | Reference label; not parsed for identity or charge. |
| `species` | — | Yes | Target label; not parsed for identity or charge. |
| `mw` | Yes | Yes | Neutral molecular mass in Da; used for mass correction and automatic extraction width. |
| `z` | Yes | Yes | Positive integer charge state, such as `8` or `1`. |
| `mz_min`, `mz_max` | Preferred pair | Alternative pair | Inclusive extraction bounds in Th; calibrant bounds must be finite, positive and increasing. |
| `mz` | Legacy alternative | Yes, or derived from bounds | Legacy extraction centre in Th. |
| `CCS_reported` | Yes | — | Known reference CCS in Å²; analyte CCS is calculated from the fitted curve. |
| `raw_path` | Yes | Yes | Complete acquisition directory ending in `.raw`, including its metadata. |

```csv
name,mw,z,mz_min,mz_max,CCS_reported,raw_path
UBQ z=8,8565,8,1069.63,1073.63,1622,raw/ubiquitin.raw
```

```csv
species,mw,z,mz,raw_path
[ADP+H]+,427.2,1,428.21,raw/adp_edda.raw
```

The preferred reference header is `CCS_reported`; legacy `ccs` is also accepted
(case-sensitive). If both are supplied, their numeric values must agree. This
reference input is separate from calculated analyte CCS.

For calibrants, a supplied `mz_min`/`mz_max` pair takes precedence over legacy
`mz`/`mz_window`: the program derives their midpoint and half-width. Supply both
bounds or neither; a partial pair is rejected. With both omitted/blank, legacy
`mz` plus numeric or `auto` half-width remains accepted.

These rows illustrate the format. A fit requires at least three valid reference
peaks after signal thresholds and EDC correction; one row cannot fit a curve.
Use the supplied 25-reference table for the protein → ADP example. Betaine alone
supplies one reference and cannot fit the free two-parameter curve. Its separate
example fits six tune-mix references before applying the curve to Betaine.

## Optional columns

Omitted or blank settings use these defaults. A blank `formula` uses the
mass-based automatic window estimate.

| Column | Applies to | Default | Meaning / supported values |
| --- | --- | --- | --- |
| `mz_window` | Both | `auto` | Legacy extraction half-width in Th, such as `2`, or `auto`; explicit bounds take precedence. |
| `formula` | Both | Blank | Elemental formula used only when `mz_window` is `auto`. |
| `min_ms_intensity` | Both | Calibrants `100`; analytes `50` | Minimum maximum MS intensity inside the target window, in reader intensity units. |
| `peak_select` | Calibrants | `highest` | `highest`: greatest accepted height; `last`: latest time; `manual`: requires a detected native `peak_drift_bin`. All use signal thresholds; legacy `tallest` / `latest` aliases also work. |
| `peak_drift_bin` | Calibrants | Blank | Manual override: zero-based integer bin of a detected candidate passing the thresholds; takes priority over `peak_select`. |
| `min_im_intensity` | Calibrants | `500` | Minimum peak height of the summed IM drift profile, in reader intensity units. |
| `peak_height_frac` | Calibrants | `0.10` | Peak-detection height as a fraction of the profile maximum; greater than 0 and at most 1. |
| `gas_mw` | Calibrants | `28.014` | Drift-gas molecular mass in Da; nitrogen default. All references must agree. |

The program uses supplied `mw` and `z`, and the explicit-bounds midpoint or
legacy `mz`; it does not calculate mass-to-charge from a label or molecular mass. It measures intensities and drift peaks from data, reads
pusher timing and EDC from RAW metadata, fits references, and calculates analyte
CCS. Calibrant-only thresholds and peak selection do not control analyte extraction.
`direct` is the power-law CCS method; `twostep` uses the diagnostic linear conversion.

## RAW paths

Relative paths resolve beside the CSV. If the CSV is in `my_batch/`,
`raw/sample.raw` means `my_batch/raw/sample.raw`; `raw\sample.raw` also works.
The converted pair must be named `sample_ms.txt` and `sample_im.txt`.
Convert every referenced acquisition; `--data-dir` must contain the matching
MS/IM pairs for the calibrants and analytes.

Windows absolute paths include `C:\Waters data\sample.raw` and
`\\server\share\sample.raw`. Use literal backslashes in CSV, without Python-style
doubling. Quote a field containing a comma, and optionally one containing spaces:
`"C:\Waters data\sample.raw"`. Quote CLI path arguments containing spaces.
Keep original RAW folders reachable for drift timing and EDC metadata.

## Profile-target CSV

For the current acquisition, **Targets → Import CSV** accepts:

```csv
species,mz_min,mz_max,drift_mode,drift_min,drift_max
Interest 1,1000,1200,auto,,
Interest 2,1200,1400,manual,3,13
```

| Column | Meaning |
| --- | --- |
| `mz_min`, `mz_max` | Required finite increasing inclusive bounds in Th. |
| `species` | Optional label; `name` is a fallback. |
| `drift_mode` | `auto` (default) or `manual`; the GUI shows an Auto drift checkbox. |
| `drift_min`, `drift_max` | Finite increasing display bounds for manual mode, in ms or bins without timing. Auto fills all valid rows on data load/import or mass-bound edits, before selection. |
| `mw`, `z` | Optional molecular mass in Da and positive integer charge; required only for CCS. |
| `mz`, `mz_window` | Derived midpoint and half-width; retained for CCS compatibility. |
| `raw_path` | Bound to the currently open acquisition during import. |

Auto drift uses every positive original drift sample within the inclusive mass
interval, independently of display noise/smoothing. It applies equal padding of
20% of the larger absolute endpoint and clips to acquisition bounds. With no
positive signal it uses the full acquisition drift range. Manual bounds are used
unchanged. Selecting a row focuses both 2D axes and the drift profile;
**Full view** deselects the target and restores the full acquisition.

Existing analyte CSVs still import: `mz` plus a positive `mz_window` supplies the
bounds; `auto` needs `mw` and `z` and may use `formula`. The GUI recomputes midpoint
and half-width from explicit bounds. `drift_mode` values `false`, `0`, `no` or
`off` also select manual mode; use `auto`/`manual` for new sheets.

**Export CSV** saves inputs and any computed results, including `peak_ccs` (Å²),
`peak_drift_time_ms`, `status`, `extrapolated`, `calibration_review_required`, and
`conditions_status`; imported extra columns are retained. Import restores these
results as read-only values only when the original RAW path matches the current
acquisition. Mismatched or missing acquisition provenance clears those results
and marks them stale while preserving target inputs. Keep results with their
original calibration. This target-sheet export is separate from **Export → Drift CSV**,
which saves the currently plotted profile samples.

## GUI editing and reusable calibration

In **Spectrum → Targets**, use **Add**, edit the name and bounds, and select a row
for profiling. **Delete** removes that row; **Import/Export CSV** loads or saves
the target sheet, including computed results.
Mass and charge may remain blank until CCS is requested.

The **Calibration** tab has linked plots above one editable reference table.
Use **Load CSV**, or **Add standard RAW → + Row**; enter bounds and known
reference properties. **- Row** removes rows and **Save CSV** saves a copy.
**Extraction options** shows thresholds. Numeric cells show a compact value;
hover for full precision. Edits pass as a snapshot without rewriting the input CSV.

**Process references** detects candidates without fitting; follow processing
status and inspect each chosen $t_d$. Click a row to link its profile and circled
curve point; click the curve point to select the row. **Highest/Last** buttons
apply to **Selected ions** (Ctrl/Shift multiselect) or **All ions**, with immediate
peak previews. Manual choice clicks a detected profile marker and records its
zero-based native `peak_drift_bin`; an arbitrary time or unaccepted bin cannot
substitute. Manual never falls back to an automatic choice.

`selected_peak` is the 1-based candidate number in native drift-bin order;
`selected_td_ms` is its acquired bin-centre time in ms. Both are read-only display
fields with full-precision tooltips. Exported display fields are ignored when
loading a reference CSV and regenerated; they are not peak-selection inputs.

Peak/mode edits update previews immediately and invalidate the previous CCS results.
Click **Fit calibration**, inspect the selected peaks and fit warnings, then click
**Generate CCS** to accept the fit and show native CCS profiles below Axis ranges in the right column.
**Load calibration** restores a JSON curve; **Save calibration** stores the current
fit for reuse. Loaded curves also require **Generate CCS** before target processing.
After generation, complete targets update automatically when their inputs or the
acquisition change. Targets need mass and charge; incomplete rows stay uncalculated.
Extraction edits require **Process references** again. Saving preserves results.
The interactive CCS curves use the acquired samples; GUI updates skip PNG generation.

These demo assignments and reference values are inherited. The protein fit is
poor; ADP and Betaine extrapolate below their reference drift intervals. Review
peak assignments and acquisition conditions before treating CCS as quantitative.
