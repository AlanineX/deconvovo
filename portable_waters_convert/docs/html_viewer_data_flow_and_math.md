# Waters viewer data flow and math

The desktop draws native Matplotlib plots from converted TXT data, using
NumPy/SciPy for the numerical transforms. It does not embed a browser or load
HTML. **Input/Output → Process files** prepares the tables and native plots,
reusing validated unchanged conversions. **Export portable HTML** produces the
separate Plotly document for an external browser.

The sections below retain the complete HTML data-flow explanation and describe
the corresponding native calculations. Both routes use the same physical
bin-centre and display conventions, with distinct rendering implementations.

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

The wide **Input/Output** rows and **Process files** discover the selected RAW input.
**CPU** limits independent reader processes, each writing an isolated staged
pair; display calculations still use one selected acquisition. HTML and optional
CCS stages run sequentially after all requested pairs succeed. A reader failure
cancels queued work and stops downstream stages; failed refresh caches remain
invalid even if earlier output files survive.

Native viewing reads the paired files directly. Its numerical interface is
[`deconvovo.imms_view`](../deconvovo/imms_view.py):
`load_view_data(im_path, ms_path=None, pusher_us=None, config=None)`, returning
`ViewData`; `ViewData.render(config=None, viewport=None)` produces a
`ViewSnapshot` for plotting and current-view export. Every native render derives
from the original TXT values at their available precision, rather than from a
previous display or an HTML payload. Range changes can widen the native view
directly from those source data. Without a configured mass window, the full
mass domain spans all original MS/IM coordinates, including zero-only outer
endpoints; an explicit starting window still clips it normally.
The desktop loads the run once in a background worker and reuses those immutable
arrays. A numeric range edit replaces a mouse-selected viewport; other display
edits retain it. **Home/Reset view** retains current filters and restores full original mass
coordinates (including zero endpoints) and the acquired physical drift extent.

Native hover reports m/z, drift coordinates and the processed heatmap cell
value $Z_{out}[r,c]$ in the toolbar. Wheel zoom needs no mode selection.
**Zoom** selects a dragged rectangle; **Pan** moves the view. Gestures on a
marginal change its physical coordinate rather than directly zooming intensity.
Linked projections recompute after drag release; wheel changes are debounced.
**Back/Forward** restore earlier physical bounds; **Home/Reset** restore the
full acquired ranges while retaining display filters.

For explicit HTML generation, `imms_plot` pairs the files by run name and calls
`plot_im_data()`. The embedded `rawIM` payload preserves the available floating-
point m/z coordinates and intensities; it does not round them to four decimal
places or integer counts. The embedded profile axis also retains acquired
zero-intensity sample coordinates from the MS and IM tables. Initial Python
histograms and browser rebinning use the same input precision. Native and HTML
sums can differ by floating-point accumulation order, without a reduced profile
coordinate axis or serialization rounding.

### Step 2: point cloud to heatmap bins

The 2D panel is a matrix:

$$
H[r, c]
$$

$r$ is a drift row and $c$ is an m/z column. Normally `initial_view.mz_range`
sets the heatmap's starting mass domain. For a very narrow Gaussian-smoothed
window, the internal processing domain expands to include real surrounding data
needed by the physical kernel. Its width also accounts for the finest original
mass spacing, subject to the available acquisition extent. The
visible limits remain exact; this changes only 2D convolution, without widening
1D selections, adding profile samples or disabling filters.

Here $[m_{min},m_{max}]$ denotes those internal processing bounds, which may
extend beyond the visible interval. For $N_m$ m/z bins, the cell width is:

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

Concrete example: if processing bounds are `[350, 1250]` with `6400` m/z bins:

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
$\widetilde d_r$ in drift-bin units. HTML rebinning derives its processing
context from the fixed requested mass bounds and current bin/Gaussian settings;
mouse zoom changes visible limits without changing that context. The native
renderer can rebuild its domain from original TXT when the requested window
widens. It does not require a new HTML file.

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

Here $\Delta t$ is the current displayed drift-cell width. It is the difference
between adjacent drift centres, or twice the first centre when there is one
display cell. For timed data:

$$
\Delta t = \frac{N_D^{\mathrm{native}}}{N_D}
\frac{pusher\_us}{1000}
$$

Using the native pusher step after drift rebinning would give the wrong
smoothing width. Both native and HTML routes now use the displayed cell width,
keeping the configured physical Gaussian width stable as drift-bin count changes.

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

#### Savitzky–Golay coefficient correction

For an odd window width $w$, signal length $n$, and polynomial order $k$, the smoothing coefficients
$c_j$ evaluate the local least-squares polynomial at its centre. They preserve
a constant signal before the nonnegative display clamp:

$$
\sum_{j=0}^{w-1} c_j = 1,
\qquad
\widehat y_i = \max\left(0,\sum_{j=0}^{w-1}c_j
y_{\operatorname{clip}(i+j-h,0,n-1)}\right),
\qquad h=\frac{w-1}{2}
$$

The former JavaScript implementation had an unused loop that repeated the
coefficient sum $k+1$ times, introducing a gain factor $k+1$ per pass. That loop
is removed, and the native implementation uses SciPy's correctly normalized
Savitzky–Golay coefficients with nearest-endpoint padding. This corrects a numerical scaling defect; it does not
change the acquired data or provide a calibration correction.

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

In exported HTML, the `Render` tab controls Plotly's interpolation between heatmap cells:

- `Crisp`: no interpolation between cells
- `Smooth (fast)`: bilinear interpolation
- `Smooth (best)`: bicubic interpolation

Rendering changes how pixels are drawn on screen. It does not change $H$,
$Z$, exported CSV values, or the scientific binning.
Native Matplotlib uses nearest-cell, bilinear or bicubic interpolation for the
same `off`, `fast` and `best` settings.

### Step 7: marginal profiles

The HTML viewer reads the resolved heatmap ranges from Plotly's `_fullLayout`.
Matched axes can leave stale ranges in the editable `layout`, so those saved
values are not used as the authority for the plotted viewport. Native plots use
shared Matplotlib axes and the native snapshot's viewport.

Both 1D profiles project original IM points, independently of the 2D histogram's
mass and drift bin counts. Shipped defaults use raw drift smoothing and zero
1D noise thresholds; explicit settings can enable processing. The drift axis
keeps every original drift bin,
including trailing zero-intensity bins. Its native coordinates are $t_d$ as
defined above, or $d+0.5$ without pusher timing. For the visible m/z interval:

$$
D[d] = \sum_{i:\,x_{lo}\leq m_i\leq x_{hi},\,d_i=d} I_i
$$

`Drift Noise` thresholds this full native profile and `Drift Smooth` smooths its
intensities while retaining its original coordinates. There is no additional
drift point limit or rebinning caused by either 2D display-bin control.

Drift-only zoom crops the native curve and rescales its intensity axis; sparse
views mark existing samples without resampling or changing exported sample values.
Changing the extraction mass window recomputes the drift profile; selected targets
retain their explicit mass bounds. The heatmap's Gaussian smoothing and
cube-root/other intensity transforms affect only its 2D appearance.

The bottom m/z profile uses the acquired coordinate axis. Let $x_j$ be the sorted
unique union of all m/z coordinates in the paired MS and IM tables, including
zero-intensity samples. The MS table contributes coordinates and gaps; its
intensities do not replace the drift-linked IM sums. Let $B_{vis}$ contain the
original drift bins whose cell edges overlap the current drift viewport. At each
original coordinate:

$$
M_j = \sum_{i:\,m_i=x_j,\,d_i\in B_{vis}} I_i
$$

Here $N_d$ is the original drift-bin count and $\Delta t$ is the native pusher
period in milliseconds, or one bin without timing. The selected bin bounds are:

$$
k_{lo}=\max\left(0,\left\lfloor\frac{t_{lo}}{\Delta t}\right\rfloor\right),
\qquad
k_{hi}=\min\left(N_d-1,\left\lceil\frac{t_{hi}}{\Delta t}\right\rceil-1\right)
$$

These bounds apply to a positive-width intersection with the original drift
domain; otherwise no bins contribute. Repeated IM points at the same coordinate
sum; coordinates with no selected signal remain zero. Native viewing caches the
full-drift projection; both viewers recompute selected-drift sums on the same
coordinate axis. The trace retains $x_{lo}\leq x_j\leq x_{hi}$ within the heatmap
domain, inclusively; equality after domain clipping retains an exact acquired
endpoint. An entirely out-of-domain mass interval produces an empty trace.
No interpolated points, adaptive grid, minimum spacing or hidden point budget is
used in the mass profile.

`m/z Noise` thresholds the visible profile as a percentage of its maximum.
There is no additional m/z-profile smoothing. Profile mode retains every
original sample, including zero gaps. Stick mode selects local maxima above 1%
of the visible maximum that are highest within $[x_j-0.08,\,x_j+0.08]$ Th;
stick positions remain original coordinates. Endpoints are excluded, and a flat
peak uses its middle sample, rounded down. Native Matplotlib and exported
Plotly lines disable path simplification, preserving all trace vertices. A
single original mass sample uses a point marker; no neighboring samples are
invented to draw a line.

### Target m/z → drift profile and focused view

A target's mass interval is an explicit extraction input, separate from CCS
reference sheets. Selecting a target now focuses the 2D axes to that interval
and its auto/manual drift bounds.
`ViewData.extract_drift_profile(mz_range, config=None)` selects original native
TXT points by an inclusive m/z interval, before mass-axis histogramming:

$$
R[d] = \sum_{i:\,m_{lo}\leq m_i\leq m_{hi},\,d_i=d} I_i
$$

It retains the same original native drift-bin coordinates regardless of 2D
binning. Only drift noise and drift smoothing apply to $R$. The 2D Gaussian
settings, intensity transform, cutoff/fade and mass-profile noise do not enter
this extraction. A selected **Targets** row drives both the drift panel beside
the heatmap and the **Mobiligram** tab. Its exact mass bounds rebuild the native
heatmap domain from the cached original TXT points; no RAW conversion is needed.
Both linked profiles then follow the focused physical viewport. **Full view**
deselects the target and restores the acquisition axes. **Home/Reset view** clears
preset bounds and mouse zoom to restore the complete original acquisition while
retaining display filters.

Auto drift bounds fill every valid target row on data load/import and mass-bound
edits, without selecting or focusing that row. Manual bounds remain unchanged.
**Auto drift** determines focus from every positive original raw drift sample
inside the target's inclusive mass bounds. Detection disables drift noise and
smoothing, so visual settings cannot trim the detected support. If the first
and last positive native sample centres are $a$ and $b$, equal padding is applied
at both ends before clipping to the acquired physical extent $[0,T_{max}]$:

$$
p=0.2\max(|a|,|b|),\qquad
 t_{min}=\max(0,a-p),\qquad t_{max}=\min(T_{max},b+p)
$$

Thus a 5–11 ms support gives 2.8–13.2 ms before acquisition clipping. With no
positive target signal, Auto uses the full acquired drift extent. Manual
`drift_min`/`drift_max` values are used unchanged. These bounds select the view;
the extracted profile still retains all original drift sample coordinates.
For $N_d$ original bins and native width $\Delta t$, the acquired extent is:

$$
T_{max}=N_d\Delta t
$$

For example, ADP EDDA has 200 bins at 69 µs, so its full drift extent is 13.8 ms.
Target-sheet **Export CSV** stores these inputs and any computed results;
current-view profile CSV exports the displayed sample trace within the view.

**Export → Drift PNG/CSV** captures the selected profile within the current drift
viewport; its filename contains the target name. The exported intensities and
coordinates are the extracted, drift-processed values. Molecular mass, charge,
reference CCS and analyte/calibrant sheets are unnecessary for profiling.

### Step 8: current-view PNG and CSV

The desktop's **Export** menu captures its native numerical
snapshot and Matplotlib plots. Exported HTML has its own current-view functions
that wait for pending display changes and read the rendered Plotly traces and
resolved viewport. Each route captures its displayed values, rather than
running a separate raw-data reconstruction during export.

| Export | Contents |
| --- | --- |
| Full plot PNG | Current heatmap and both linked profiles, with the current zoom and settings. |
| 2D PNG | Current heatmap, colorscale and axis ranges. |
| m/z PNG | Current displayed mass profile or sticks and resolved axis ranges. |
| Drift PNG | Current displayed drift profile and resolved axis ranges. |
| 2D CSV | Visible overlapping heatmap cells, retaining their original centres and displayed values. |
| m/z CSV | Displayed mass-profile/stick vertices inside the current mass range. |
| Drift CSV | Displayed drift-profile vertices inside the current drift range. |

For 2D CSV, visible cells are selected by their edges; a partly visible cell
retains its full value. After intersecting the visible m/z interval with the
heatmap domain, the column bounds are:

$$
c_{lo} = \max\left(0,
\left\lfloor\frac{x_{lo}-m_{min}}{\Delta m}\right\rfloor\right),
\qquad
c_{hi} = \min\left(N_m-1,
\left\lceil\frac{x_{hi}-m_{min}}{\Delta m}\right\rceil-1\right)
$$

This selects cells only for a positive-width intersection. The same rule
selects 2D drift cells using their displayed cell width. These export cell bounds
do not rebin either 1D profile.

The heatmap CSV header is `mz,drift_time_ms,display_value,intensity_transform`
when pusher metadata is available; otherwise it uses `drift_bin` for the second
column. Its `display_value` is $Z_{out}[r,c]$: it includes current smoothing,
intensity transformation and cutoff/fade. It is not the original point
intensity $I_i$ or unprocessed heatmap sum $H[r,c]$.

Profile CSV headers are `mz,display_value` and
`drift_time_ms,display_value` (or `drift_bin,display_value` without timing).
They retain the profile's own noise/smoothing/display-mode settings. The profile
CSVs contain the actual displayed trace coordinates and values, without export
resampling or a point limit. Blank rows retain the breaks between sticks.
A partially visible heatmap cell is exported with its full displayed value and original centre; the export does not split or
apportion that cell's intensity.

The generated `run_2d_imms.csv` is a separate export of the original
positive-intensity points with
`mz,drift_bin,intensity` and optional `drift_time_ms`. Keep that file and the
original TXT/RAW inputs for numerical work on raw intensities. The older browser
CSV buttons provide their named raw/processed exports. Their full-domain drift
exports also retain original bin centres and precision, without the former
0.02-ms resampling. The legacy raw MS CSV also preserves original precision
and includes only acquired points within its chosen mass interval. They are
distinct from
the desktop's native current-view CSVs. Native RAW opening does not create
this point CSV or HTML automatically; they belong to explicit HTML generation.

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
