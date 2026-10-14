# HTML viewer data flow and math

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

This is the bin-centre convention in the initial heatmap, the point CSV,
the linked profiles, and CCS extraction. It also applies after display
rebinning. Native bin zero is centred at 0.055 ms when the pusher period is
110 microseconds; it is not drawn at zero milliseconds.

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

Each in-domain IM point is assigned to one m/z column. Cells include their
lower edge, and the last cell also includes the domain's upper edge:

$$
c_i =
\begin{cases}
N_m-1, & m_i=m_{max} \\
\left\lfloor \frac{m_i-m_{min}}{\Delta m} \right\rfloor,
  & m_{min}\leq m_i<m_{max} \\
\text{discard}, & \text{otherwise}
\end{cases}
$$

Python's initial histogram and JavaScript rebinning follow this same rule, so
a point exactly on the right boundary is counted once in the final cell.

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

For $N_D^{\mathrm{native}}$ native drift bins and $N_D$ displayed drift bins,
the displayed coordinate is the centre of the corresponding interval:

$$
\widetilde d_r = \left(r+\frac{1}{2}\right)
\frac{N_D^{\mathrm{native}}}{N_D},
\qquad
t_r = \widetilde d_r\frac{pusher\_us}{1000}
$$

Without pusher metadata, the viewer retains the bin-centre coordinate
$\widetilde d_r$ in drift-bin units. Rebinning retains the fixed baked m/z
domain; zoom changes the visible interval, rather than moving that domain.

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
f(x) = \sqrt[3]{x}
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

The viewer reads the resolved heatmap ranges from Plotly's `_fullLayout`.
Matched axes can leave stale ranges in the editable `layout`, so those saved
values are not used as the authority for the plotted viewport.

Visible cells are selected by their edges. A partially visible cell still
contributes to the projection. After intersecting the displayed m/z interval
$[x_{lo},x_{hi}]$ with the baked domain, the column bounds are:

$$
c_{lo} = \max\left(0,
\left\lfloor\frac{x_{lo}-m_{min}}{\Delta m}\right\rfloor\right),
\qquad
c_{hi} = \min\left(N_m-1,
\left\lceil\frac{x_{hi}-m_{min}}{\Delta m}\right\rceil-1\right)
$$

This applies only when the intersected interval has positive width. An
entirely out-of-domain viewport selects no cells. The same floor/ceiling rule
selects drift cells using their interval width. Both marginal plots therefore
follow the heatmap ranges actually drawn on screen.

The left drift profile is recomputed from the current visible m/z window:

$$
D[r] = \sum_{c=c_{lo}}^{c_{hi}} H[r,c]
$$

Then `Drift Noise` thresholds that 1D profile and `Drift Smooth` smooths it.

The bottom m/z profile is rebuilt from the raw per-drift IM point lists inside
the current visible drift window. `m/z Noise` thresholds that profile. `MS Mode`
chooses whether to show the profile as a line or peak-pick it into sticks.

### Step 8: current-view PNG and CSV

The desktop's **Export current view** menu uses the same rendered Plotly traces.
It waits for pending display changes, refreshes the linked profiles, and reads
the resolved viewport. It does not run a separate raw-data reconstruction.

| Export | Contents |
| --- | --- |
| Full plot PNG | Current heatmap and both linked profiles, with the current zoom and settings. |
| 2D PNG | Current heatmap trace, colorscale and resolved axis ranges. |
| m/z PNG | Current displayed mass profile or sticks and resolved axis ranges. |
| Drift PNG | Current displayed drift profile and resolved axis ranges. |
| 2D CSV | Visible overlapping heatmap cells, retaining their original centres and displayed values. |
| m/z CSV | Displayed mass-profile/stick vertices inside the current mass range. |
| Drift CSV | Displayed drift-profile vertices inside the current drift range. |

The heatmap CSV header is `mz,drift_time_ms,display_value,intensity_transform`
when pusher metadata is available; otherwise it uses `drift_bin` for the second
column. Its `display_value` is $Z_{out}[r,c]$: it includes current smoothing,
intensity transformation and cutoff/fade. It is not the original point
intensity $I_i$ or unprocessed heatmap sum $H[r,c]$.

Profile CSV headers are `mz,display_value` and
`drift_time_ms,display_value` (or `drift_bin,display_value` without timing).
They retain the profile's own noise/smoothing/display-mode settings. Blank rows
retain the breaks between sticks. A partially visible heatmap cell is exported
with its full displayed value and original centre; the export does not split or
apportion that cell's intensity.

The generated `run_2d_imms.csv` is a separate export of the original
positive-intensity points with
`mz,drift_bin,intensity` and optional `drift_time_ms`. Keep that file and the
original TXT/RAW inputs for numerical work on raw intensities. The older browser
CSV buttons provide their named raw/processed exports; they are distinct from
the desktop's current-view CSVs.

**Export portable HTML** regenerates one complete interactive document from the
converted tables and current defaults. It includes Plotly and the run data,
but it does not capture temporary zoom changes like current-view PNG/CSV does.

### Betaine 1000–2000 visual example

The focused example is the Betaine acquisition, separate from the protein → ADP
CCS workflow. Its visible high-mass features remain unassigned; the Betaine ion
near 118 is outside this focus. The display configuration uses 1000–2000 and
4–12 ms, cube-root intensity, 6400 mass bins, native drift bins, Gaussian 2D
widths of 0.1 Da and 0.07 ms, 5% cutoff and 1% fade. The mass profile stays in
profile mode and the drift profile is unsmoothed. These display choices do not
change the calibration inputs or establish a chemical assignment.
