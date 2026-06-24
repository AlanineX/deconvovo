# DeconVoVo

An analysis suite for Waters SYNAPT travelling-wave ion mobility mass spectrometry (TW-IMS) data — from raw instrument files to collision cross-sections.

## Getting started

**Windows** — download `DeconVoVo-win64.zip` from [Releases](../../releases), extract, and double-click `DeconVoVo.exe`.

**Linux / macOS** — install with pip and launch the GUI or use the command line:
```bash
pip install -e ".[all]"
python -m gui.app
```

> Waters `.raw` conversion uses UniDec's `CDCReader.exe`. On Windows it runs
> natively; on Linux/macOS install `wine` first (`sudo apt install wine` /
> `brew install --cask wine-stable`). All other CLI steps are pure Python and
> work on any OS without Wine.

## Workflow

1. **Convert** — extract MS and IM data from Waters `.raw` directories to text files
2. **Visualize** — generate interactive 2D IM-MS heatmaps (HTML) with linked drift and m/z panels, smoothing controls, and CSV export
3. **CCS Calibrate** — build TW-IMS CCS calibration curves from calibrant species. Supports both protein and small-molecule calibrants (e.g. Agilent tune mix). Direct (power-law) and two-step fitting methods
4. **Analyte CCS** — apply calibration to compute CCS profiles, with resampling and smoothing on a uniform CCS grid

All parameters are CSV-driven — m/z extraction windows, peak selection strategy, intensity thresholds, and drift gas mass are specified per-row in the config. No hardcoded assumptions. See `config/` for examples.

## HTML viewer data flow and math

The HTML viewer turns each run into an interactive image from a 3-column IM point
cloud. The three dimensions are m/z, drift, and intensity:

$$
P = \{(m_i, d_i, I_i)\}_{i=1}^{N}
$$

`m_i` is m/z, `d_i` is the native drift bin, and `I_i` is intensity. When the
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

`r` is a drift row and `c` is an m/z column. If `initial_view.mz_range` is set,
the initial heatmap spends its m/z bins only on that starting range. For the
Carter config this is 350-1250 m/z.

For `N_m` m/z bins over `[m_min, m_max]`, the m/z bin width is:

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

A point at `m/z = 700.10` goes to:

$$
c = \left\lfloor \frac{700.10 - 350}{0.140625} \right\rfloor = 2489
$$

The `m/z Bins` tab rebuilds this matrix with a different `N_m`. The `Drift Bins`
tab groups native drift rows into fewer or more displayed rows.

### Step 3: 2D Smooth and Scale

The browser applies 2D smoothing to `H` before intensity scaling and cutoff/fade
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
- `2`: fatter in m/z
- `sqrt`: multiply m/z width by the square root of the heatmap aspect ratio
- `full`: multiply m/z width by the full heatmap aspect ratio

For the Carter default Gaussian preset, `sigma = 0.25 Da` and
`sigma_dr = 0.15 ms`. With the `6400`-bin example above:

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

The cutoff threshold uses `Z`, not the original raw intensity matrix. This makes
the noise controls match the visible color range.

### Step 5: Noise Cutoff and Cutoff Fade

This step takes each displayed heatmap value `Z[r,c]` and returns a final
display value `Z_out[r,c]` before color mapping. The two controls are connected
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

For example, if `Z_max = 100` and `Noise Cutoff = 5%`:

$$
T = 100 \cdot \frac{5}{100} = 5
$$

Now every heatmap value is compared to `T`.

#### Hard cutoff

`Cutoff Fade = Off` means a hard cutoff at `T`. Values below the threshold are
removed. Values at or above the threshold are kept unchanged:

$$
Z_{out}(v) =
\begin{cases}
0, & v < T \\
v, & v \ge T
\end{cases}
$$

With `T = 5`:

| Input value `v` | Comparison | Output |
| ---: | --- | ---: |
| 2 | `2 < 5` | 0 |
| 5 | `5 >= 5` | 5 |
| 8 | `8 >= 5` | 8 |

#### Fade cutoff

`Cutoff Fade = s%` creates a gradual ramp around `T`. First, the fade percent is
converted into a half-width:

$$
S = Z_{max}\frac{s}{100}
$$

The fade interval is:

$$
[T - S,\ T + S]
$$

Values below `T - S` are fully removed. Values above `T + S` are fully kept.
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

Concrete example for `Z_max = 100`, `Noise Cutoff = 5%`, and
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

| Input `v` | Weight calculation | Weight `w` | Output `v*w` |
| ---: | --- | ---: | ---: |
| 2 | `(2 - (-15)) / 40` | 0.425 | 0.85 |
| 5 | `(5 - (-15)) / 40` | 0.500 | 2.50 |
| 8 | `(8 - (-15)) / 40` | 0.575 | 4.60 |

This is why larger `Cutoff Fade` values do not simply remove more signal. They
make the cutoff edge broader, so very low values can survive faintly while
values near the threshold are dimmed.

#### Side-by-side example

For `Z_max = 100`, `Noise Cutoff = 5%`, and three values `2`, `5`, and `8`:

| Setting | `Z = 2` | `Z = 5` | `Z = 8` |
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

Rendering changes how pixels are drawn on screen. It does not change `H`,
`Z`, exported CSV values, or the scientific binning.

### Step 7: marginal profiles

The left drift profile is recomputed from the current visible m/z window:

$$
D[r] = \sum_{c=c_{lo}}^{c_{hi}} H[r,c]
$$

Then `Drift Noise` thresholds that 1D profile and `Drift Smooth` smooths it.

The bottom m/z profile is rebuilt from the raw per-drift IM point lists inside
the current visible drift window. `m/z Noise` thresholds that profile. `MS Mode`
chooses whether to show the profile as a line or peak-pick it into sticks.

## Command line

Every step runs as `python -m deconvovo.<module>`. Pass `--help` to any of them
for the full flag list.

```bash
# 1. Convert Waters .raw → text (_ms.txt + _im.txt)
python -m deconvovo.waters_convert data/*.raw -o converted/

# 2. Build interactive 2D IM-MS HTML viewers
python -m deconvovo.imms_plot -i converted/ -o html/ --raw-dir data/ -j 8

# 3. (optional) CCS calibration + analyte CCS in one step
python -m deconvovo.imms_ccs_calibrate \
    -o output/ccs \
    --data-dir data/ \
    --calibrant-csv config/calibrants_tunemix.csv \
    --analyte-csv config/analytes_adp.csv

# 4. Per-run peak summary CSV
python -m deconvovo.imms_summary -i html/
```

### One-shot pipeline

`deconvovo.cli` runs steps 1 + 2 end-to-end (and 3/4 if their inputs are
present):

```bash
python -m deconvovo.cli -i data/ -o output/ \
    --mass-range 5000 100000 \
    --charge-range 1 30 \
    --mass-bins 10 \
    -j 8
```


## Command line for me

Uses the custom config at `/home/alan/Downloads/carter_0623/imms_plot_config.json`
(starting m/z 350–1250, drift 0–15 ms, W-Spectral colormap, 0.5 Da m/z /
0.18 ms drift 2D Gaussian smooth, sqrt m/z scale, 15% Noise Cutoff, 50% Cutoff
Fade, profile m/z mode, Smooth render, sqrt intensity, and 6400 m/z bins). Edit
that file to change any default — the dropdowns in the HTML still work for
runtime tuning.

One-shot — convert + viewers in a single invocation:

```bash
python -m deconvovo.cli \
    -i /home/alan/Downloads/carter_0623/Myoglobin_D2O \
    -o /home/alan/Downloads/carter_0623/Myoglobin_D2O_output \
    --config /home/alan/Downloads/carter_0623/imms_plot_config.json \
    --skip-deconv \
    -j 8
```

Or step-by-step if you want to inspect the intermediate text files:

```bash
# 1. Convert .raw to text
python -m deconvovo.waters_convert \
    /home/alan/Downloads/carter_0623/Myoglobin_D2O \
    -o /home/alan/Downloads/carter_0623/Myoglobin_D2O_converted

# 2. Generate interactive HTML viewers (custom config)
python -m deconvovo.imms_plot \
    -i /home/alan/Downloads/carter_0623/Myoglobin_D2O_converted/*_24H_* \
    -o /home/alan/Downloads/carter_0623/Myoglobin_D2O_html \
    --raw-dir /home/alan/Downloads/carter_0623/Myoglobin_D2O \
    --config /home/alan/Downloads/carter_0623/imms_plot_config.json \
    -j 8

# 3. (optional) CCS calibration + analyte CCS
python -m deconvovo.imms_ccs_calibrate \
    -o /home/alan/Downloads/carter_0623/Myoglobin_D2O_ccs \
    --data-dir /home/alan/Downloads/carter_0623/Myoglobin_D2O \
    --calibrant-csv config/calibrants_tunemix.csv \
    --analyte-csv config/analytes_adp.csv
```
