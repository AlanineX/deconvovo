# DeconVoVo · Waters Instrument Console

Process Waters SYNAPT travelling-wave ion-mobility mass spectrometry data in one workspace: browse acquisitions, convert them, inspect measured spectra, export interactive viewers, and optionally review CCS calibration.

The desktop follows the **Instrument Console (design 01)**: a dark workbench with a run list, large plots, and settings beside the data. Processing and viewing share one batch; calibration and analyte application share one set of paths and method settings.

## Windows: download and open

1. Download [WatersConsole-win64.zip](https://github.com/AlanineX/deconvovo/releases/download/v0.3.0/WatersConsole-win64.zip) from [v0.3.0](https://github.com/AlanineX/deconvovo/releases/tag/v0.3.0).
2. Extract the **whole ZIP** to a local folder.
3. Open `WatersConsole/WatersConsole.exe`.
4. Choose a demo or your Waters `.raw` folder, choose a writable output folder, then click **Process & view**.

Python, Wine, UniDec, and internet access are unnecessary for the Windows application. Keep the executable and its `_internal` folder together. Downloads and source dependency installation need a connection; conversion and both viewing routes run locally afterward.

## View locally or share a portable HTML

| Route | How to use it | What travels with it |
| --- | --- | --- |
| Desktop viewer | Select a processed run in the console | The application folder and your output data |
| Portable HTML | Choose **Export portable HTML** or open a generated `*_2d_imms.html` | One HTML file with Plotly and the run data embedded |

Both routes include the enhanced mass/drift bin controls, profile/stick display, smoothing, cutoff and fade, colormaps, intensity scaling, and linked profiles. The local mass profile and mobiligram also read the measured TXT data. HTML works in a modern browser without a server or internet access.

Existing converted tables can be reused; viewer settings are saved without discarding defaults that the desktop does not edit. Processing refreshes HTML so an old viewer cannot silently hide current settings. Use **Advanced → Save defaults** to keep display choices across relaunches. An error opens the activity log and prevents a failed job from being reported as complete.

Data flow: **Waters RAW folder → MS/IM tables → desktop viewer or portable HTML → optional calibration and analyte CCS**.

## Real production demos

The download includes **nine complete real acquisitions, 300.29 MiB of RAW data**. Load a demo from the console; it selects the matching viewer profile and optional CCS tables.

| Demo | Runs | RAW size | Main use |
| --- | ---: | ---: | --- |
| Protein standards + ADP | 5 | 195.07 MiB | Batch processing, 25 calibrant rows, and 34 target rows |
| Tune mix + Betaine | 2 | 40.68 MiB | Small-molecule calibration and an independent analyte |
| Native HSP104 pair | 2 | 64.55 MiB | High-mass viewing and differing drift timing |

[Demo details and scientific limits](demo_data/production/README.md) explain the inherited reference sheets and assignments. The protein calibration has a poor fit; 30 of 34 ADP peaks fall below its reference drift interval. Betaine is also an extrapolation below the tune-mix interval despite a high calibrant fit score. No compatible HSP104 calibration is supplied. Review assignments and acquisition conditions before treating exported CCS values as quantitative measurements.

The source retains both existing demo locations so the main and portable workflows keep their paths; Git stores identical RAW files once, and the compiled ZIP carries one corpus. File-level provenance is in each demo's [manifest](demo_data/production/demo_manifest.json).

## Run from source

Use Python 3.10+; Python 3.12 is used for the pinned native build. From this repository:

```bash
python -m venv .venv
# Activate .venv for your operating system.
python -m pip install -e ".[gui]"
python -m gui.app
```

On Windows, activate with `.venv\Scripts\activate`; on Linux/macOS, use `source .venv/bin/activate`. The bundled reader runs natively on Windows. RAW conversion on Linux/macOS needs an x64 Wine installation; plotting, table analysis and browser viewing do not.

The standalone [portable Waters folder](portable_waters_convert/README.md) can be copied independently. Its converter uses Python's standard library; its desktop and scientific tools use the listed dependencies. Optional UniDec deconvolution remains available through `python -m pip install -e ".[deconv]"` and the existing command-line modules.

## Build and validation

The [Windows build workflow](.github/workflows/build-windows.yml) uses [pinned dependencies](portable_waters_convert/build/requirements.txt), runs the focused Waters tests, verifies the original demo bytes, and checks real RAW conversion, offline HTML, and the embedded viewer before and after freezing. Its native smoke test also checks that the viewer fits a 1280 × 800 laptop display. Build and validation reports are retained as workflow artifacts.

```bash
python -m pip install -r portable_waters_convert/build/requirements.txt
python portable_waters_convert/build/build_native.py --output-dir waters-build
```

The resulting folder is `waters-build/dist/WatersConsole`. [CHANGELOG.md](CHANGELOG.md) describes the upgrade from v0.2.3. Third-party reader licenses remain in [portable_waters_convert/bin](portable_waters_convert/bin/).

## Command line

Supported conversion, plotting, calibration, and optional deconvolution commands retain their existing interfaces. Use `--help` for each module's full options.

```bash
python -m deconvovo.waters_convert demo_data/production/small_molecule/raw -o output/converted --fn 1 --ms-bin 0 --im-bin 0
python -m deconvovo.imms_plot -i output/converted -o output/html --raw-dir demo_data/production/small_molecule/raw --config demo_data/production/small_molecule/viewer_config.json
python -m deconvovo.imms_ccs_calibrate -o output/ccs --data-dir output/converted --calibrant-csv demo_data/production/small_molecule/calibrants.csv --analyte-csv demo_data/production/small_molecule/analytes.csv
```

For conversion plus viewers in the existing batch CLI:

```bash
python -m deconvovo.cli -i demo_data/production/hsp104/raw -o output/hsp104 --config demo_data/production/hsp104/viewer_config.json --skip-deconv -j 2
```

Omit `--skip-deconv` after installing the optional UniDec dependency when you need the existing deconvolution lane. Conversion and desktop controls use the selected acquisition function; function 1 means `_FUNC001`, so confirm its role in the instrument method. Zero mass bin widths retain the reader's original mass resolution; the extracted tables combine scans and are not a RAW-data archive. Keep original RAW acquisitions.

## HTML viewer data flow and math

The HTML viewer turns each run into an interactive image from a 3-column IM point
cloud. The three dimensions are m/z, drift, and intensity:

$$
P = \{(m_i, d_i, I_i)\}_{i=1}^{N}
$$

$m_i$ is m/z, $d_i$ is the native drift bin, and $I_i$ is intensity. When the
pusher period is known, drift bin is displayed as drift time:

$$
t_i = (d_i + 0.5)\frac{pusher\_us}{1000}
$$

### Step 1: raw files to text files

`waters_convert` reads each Waters `.raw` directory through UniDec's
`CDCReader.exe` and writes:

- `run_ms.txt`: profile MS data, columns `m/z, intensity`
- `run_im.txt`: IM point cloud, columns `m/z, drift_bin, intensity`

`imms_plot` pairs those two files by run name and calls `plot_im_data()`.

### Step 2: point cloud to heatmap bins

The 2D panel is a matrix:

$$
H[r, c]
$$

$r$ is a drift row and $c$ is an m/z column. If `initial_view.mz_range` is set,
the initial heatmap spends its m/z bins only on that starting range. For the
Carter config this is 350-1250 m/z.

For $N_m$ m/z bins over $[m_{min},m_{max}]$, the m/z bin width is:

$$
\Delta m = \frac{m_{max} - m_{min}}{N_m}
$$

Each IM point is assigned to one m/z column:

$$
c_i = \left\lfloor \frac{m_i - m_{min}}{\Delta m} \right\rfloor
$$

The heatmap pixel value before display processing is the sum of all intensities
that land in the same drift row and m/z column:

$$
H[r, c] = \sum_{i: d_i = r,\ c_i = c} I_i
$$

Concrete example: if `[350, 1250]` is shown with `6400` m/z bins, then:

$$
\Delta m = \frac{1250 - 350}{6400} = 0.140625\ \mathrm{Da/bin}
$$

A point at $\frac{m}{z}=700.10$ goes to:

$$
c = \left\lfloor \frac{700.10 - 350}{0.140625} \right\rfloor = 2489
$$

The `m/z Bins` tab rebuilds this matrix with a different $N_m$. The `Drift Bins`
tab groups native drift rows into fewer or more displayed rows.

### Step 3: 2D Smooth and Scale

The browser applies 2D smoothing to $H$ before intensity scaling and cutoff/fade
processing.
For Gaussian smoothing, the preset stores physical widths:

$$
\sigma_m\ \mathrm{in\ Da}, \quad \sigma_t\ \mathrm{in\ ms}
$$

Those physical widths are converted to heatmap-bin units:

$$
\sigma_{m,bins} = \frac{\sigma_m \cdot scale}{\Delta m}
$$

$$
\sigma_{t,bins} = \frac{\sigma_t}{\Delta t}
$$

`Scale` multiplies only the m/z-axis smoothing width. It changes peak shape in
the 2D plot without stretching the axes:

- `0.5`: thinner in m/z
- `1`: preset m/z width
- $2$: fatter in m/z
- `sqrt`: multiply m/z width by the square root of the heatmap aspect ratio
- `full`: multiply m/z width by the full heatmap aspect ratio

For the Carter default Gaussian preset, $\sigma_m=0.25\ \mathrm{Da}$ and
$\sigma_t=0.15\ \mathrm{ms}$. With the `6400`-bin example above:

$$
\sigma_{m,bins}(scale=1) = \frac{0.25}{0.140625} = 1.78\ \mathrm{bins}
$$

$$
\sigma_{m,bins}(scale=2) = \frac{0.25 \cdot 2}{0.140625} = 3.56\ \mathrm{bins}
$$

So `Scale = 2` makes features horizontally fatter than `Scale = 1`, while
`Scale = 0.5` makes them horizontally thinner.

### Step 4: Intensity transform

After smoothing, the `Intensity` tab transforms the smoothed heatmap:

$$
Z = f(H_{smooth})
$$

Available transforms include:

$$
f(x) = x
$$

$$
f(x) = \sqrt{x}
$$

$$
f(x) = \log_{10}(x + 1)
$$

The cutoff threshold uses $Z$, not the original raw intensity matrix. This makes
the noise controls match the visible color range.

### Step 5: Noise Cutoff and Cutoff Fade

This step takes each displayed heatmap value $Z[r,c]$ and returns a final
display value $Z_{out}[r,c]$ before color mapping. The two controls are connected
in the math, but they do different jobs:

- `Noise Cutoff`: where the cutoff threshold sits.
- `Cutoff Fade`: whether that threshold edge is abrupt or gradual.

`Noise Cutoff = 0%` disables the cutoff operation. In that case:

$$
Z_{out} = Z
$$

If `Noise Cutoff` is greater than zero, the viewer first computes the displayed
maximum:

$$
Z_{max} = \max(Z)
$$

Then it converts the `Noise Cutoff` percent into an absolute threshold:

$$
T = Z_{max}\frac{noise\_2d}{100}
$$

For example, if $Z_{max} = 100$ and `Noise Cutoff = 5%`:

$$
T = 100 \cdot \frac{5}{100} = 5
$$

Now every heatmap value is compared to $T$.

#### Hard cutoff

`Cutoff Fade = Off` means a hard cutoff at $T$. Values below the threshold are
removed. Values at or above the threshold are kept unchanged:

$$
Z_{out}(v) =
\begin{cases}
0, & v < T \\
v, & v \ge T
\end{cases}
$$

With $T=5$:

| Input value $v$ | Comparison | Output |
| ---: | --- | ---: |
| 2 | $2<5$ | 0 |
| 5 | $5\ge5$ | 5 |
| 8 | $8\ge5$ | 8 |

#### Fade cutoff

`Cutoff Fade` = $s$% creates a gradual ramp around $T$. First, the fade percent is
converted into a half-width:

$$
S = Z_{max}\frac{s}{100}
$$

The fade interval is:

$$
[T - S,\ T + S]
$$

Values below $T-S$ are fully removed. Values above $T+S$ are fully kept.
Values inside the interval are multiplied by a linear weight from 0 to 1:

$$
w(v) =
\begin{cases}
0, & v \le T - S \\
\frac{v - (T - S)}{2S}, & T - S < v < T + S \\
1, & v \ge T + S
\end{cases}
$$

The final value is:

$$
Z_{out}(v) = v \cdot w(v)
$$

Concrete example for $Z_{max} = 100$, `Noise Cutoff = 5%`, and
`Cutoff Fade = 20%`:

$$
T = 100 \cdot \frac{5}{100} = 5
$$

$$
S = 100 \cdot \frac{20}{100} = 20
$$

So the fade interval is:

$$
[T-S,\ T+S] = [-15,\ 25]
$$

Now compute weights and outputs:

| Input $v$ | Weight calculation | Weight $w$ | Output $v\cdot w$ |
| ---: | --- | ---: | ---: |
| 2 | $\frac{2-(-15)}{40}$ | 0.425 | 0.85 |
| 5 | $\frac{5-(-15)}{40}$ | 0.500 | 2.50 |
| 8 | $\frac{8-(-15)}{40}$ | 0.575 | 4.60 |

This is why larger `Cutoff Fade` values do not simply remove more signal. They
make the cutoff edge broader, so very low values can survive faintly while
values near the threshold are dimmed.

#### Side-by-side example

For $Z_{max} = 100$, `Noise Cutoff = 5%`, and three values $2$, $5$, and $8$:

| Setting | $Z=2$ | $Z=5$ | $Z=8$ |
| --- | ---: | ---: | ---: |
| `Cutoff Fade Off` | 0.00 | 5.00 | 8.00 |
| `Cutoff Fade 1%` | 0.00 | 2.50 | 8.00 |
| `Cutoff Fade 5%` | 0.40 | 2.50 | 6.40 |
| `Cutoff Fade 20%` | 0.85 | 2.50 | 4.60 |
| `Cutoff Fade 50%` | 0.94 | 2.50 | 4.24 |

Use `Noise Cutoff` when the threshold should move higher or lower. Use
`Cutoff Fade` when the same threshold should be sharper or more gradual.

### Step 6: value to screen pixel

After cutoff/fade processing, the colormap maps each final value to a color:

$$
color[r,c] = colormap(Z_{out}[r,c])
$$

The `Colormap` tab changes only the color lookup. The custom `W-` colormaps are
listed first, followed by the remaining colormaps in natural sorted order.

The `Render` tab controls Plotly's interpolation between heatmap cells:

- `Crisp`: no interpolation between cells
- `Smooth (fast)`: bilinear interpolation
- `Smooth (best)`: bicubic interpolation

Rendering changes how pixels are drawn on screen. It does not change $H$,
$Z$, exported CSV values, or the scientific binning.

### Step 7: marginal profiles

The left drift profile is recomputed from the current visible m/z window:

$$
D[r] = \sum_{c=c_{lo}}^{c_{hi}} H[r,c]
$$

Then `Drift Noise` thresholds that 1D profile and `Drift Smooth` smooths it.

The bottom m/z profile is rebuilt from the raw per-drift IM point lists inside
the current visible drift window. `m/z Noise` thresholds that profile. `MS Mode`
chooses whether to show the profile as a line or peak-pick it into sticks.
