# CCS Calibration Workflow for Waters SYNAPT G2 TW-IMS

Converts drift times to collision cross-sections (CCS) using the Ruotolo & Robinson
(Nature Protocols, 2008) travelling-wave calibration protocol.

---

## Big Picture: Why Two Fits?

In a **drift-tube** IMS, CCS relates to drift time by a simple, closed-form equation
derived from kinetic theory — you can invert it directly. Travelling-wave IMS (TWIMS)
has no such closed-form: ions are propelled by a periodic voltage wave, and the
relationship between drift time and CCS depends on wave velocity, wave height, gas
pressure, and ion properties in a complex, non-analytic way.

The Ruotolo & Robinson protocol solves this empirically with **two sequential fits**
against calibrant proteins whose CCS values are known from drift-tube experiments:

1. **Fit 1 (power law)** discovers the *shape* of the drift-time-to-CCS mapping —
   specifically, the exponent $X$ in $\Omega' \propto (t'_D)^X$. This exponent is
   instrument-specific and depends on the wave parameters.

2. **Fit 2 (linear)** uses $X$ to construct a linearized variable $t''_D$ that
   collapses all calibrant species (different proteins, different charge states) onto
   a single straight line against literature CCS. The slope and intercept of this line
   are the final calibration constants.

Neither fit alone is sufficient: Fit 1 operates on *reduced* quantities (charge and
mass factored out) to discover $X$, but the reduced cross-section $\Omega'$ is not
what we want to report. Fit 2 converts back to absolute CCS ($\text{\AA}^2$) using
the now-known $X$.

---

## Data Flow Overview

```
┌─────────────────────────────────────────────────────────────┐
│  INPUT FILES                                                │
│                                                             │
│  .raw/                     _extern.inf  → IMS settings, C   │
│  _FUNC001.STS              → pusher frequency f             │
│  _im.txt (converted)       → 2D data: (m/z, drift_bin, I)  │
│  ccs_calibrant_table.json  → literature Ω for each (protein, z) │
│  adp_analytes.csv          → species list: (name, m/z, z, MW)   │
└──────────────────────────────┬──────────────────────────────┘
                               │
          ┌────────────────────┼────────────────────┐
          ▼                    ▼                    ▼
   ┌─────────────┐    ┌──────────────┐    ┌──────────────────┐
   │ Pusher       │    │ IMS Settings │    │ Peak Picking     │
   │ Detection    │    │ Validation   │    │ (per calibrant   │
   │              │    │              │    │  protein × z)    │
   │ STS binary   │    │ _extern.inf  │    │                  │
   │ → f → Δt     │    │ → [HARD]/    │    │ _im.txt          │
   │              │    │   [WARN]     │    │ + expected m/z   │
   └──────┬───────┘    └──────┬───────┘    │ → drift profile  │
          │                   │            │ → peak t_D       │
          │                   │            └────────┬─────────┘
          │                   │                     │
          │    ┌──────────────┘                     │
          ▼    ▼                                    ▼
   ┌──────────────────────────────────────────────────────────┐
   │  CALIBRATION FITTING                                     │
   │                                                          │
   │  For each picked peak:                                   │
   │    t_D (from peak picking) + C (from _extern.inf)        │
   │    → t'_D  (EDC-corrected)                               │
   │    Ω_lit (from table) + z, MW, M_gas                     │
   │    → Ω'    (reduced cross-section)                       │
   │                                                          │
   │  Fit 1: ln(Ω') vs ln(t'_D) → X, A       [power law]    │
   │    with iterative outlier rejection + re-inclusion        │
   │                                                          │
   │  Fit 2: Ω_lit vs t''_D     → slope, intercept [linear]  │
   │    where t''_D = (t'_D)^X · z · μ                       │
   │                                                          │
   │  Output: calibration_curve.json                          │
   │          calibrant_summary.csv                           │
   │          calibration_plots.png                           │
   └──────────────────────────┬───────────────────────────────┘
                              │
                              ▼
   ┌──────────────────────────────────────────────────────────┐
   │  ANALYTE CCS CONVERSION                                  │
   │                                                          │
   │  For each analyte run × each species:                    │
   │    _im.txt + species m/z → drift profile I(bin)          │
   │    bin → t_D → t'_D → t''_D → Ω   (using fitted curve) │
   │    → csv_raw/   (native bins, non-uniform CCS spacing)   │
   │                                                          │
   │    Resample to uniform CCS grid (Δ = 0.5 Å²)            │
   │    Gaussian smooth (σ = 3.0) in CCS space                │
   │    → csv_smoothed/   (uniform grid)                      │
   │    → png_raw/ png_smoothed/ png_overlay/                 │
   └──────────────────────────────────────────────────────────┘
```

---

## Step-by-Step Theory with Worked Examples

### Step 1. Pusher period — from STS binary to time axis

**Where the data comes from:** the `_FUNC001.STS` file inside each `.raw` directory
stores per-scan instrument statistics. Stat code 76 is the pusher frequency $f$ in Hz.

**What we compute:**

$$\Delta t = \left\lfloor \frac{10^6}{f} \right\rfloor \;\mu\text{s}$$

**Why:** the IM separation is recorded as discrete bins. Each bin spans one pusher
cycle. To convert bin index $i$ to a physical drift time:

$$t_D = i \times \frac{\Delta t}{1000} \;\;\text{(ms)}$$

> **Example:** $f = 9090$ Hz $\;\Rightarrow\; \Delta t = \lfloor 10^6 / 9090 \rfloor = 110\;\mu\text{s/bin}$
>
> Bin $i = 31$ $\;\Rightarrow\; t_D = 31 \times 110 / 1000 = 3.410$ ms

**Data flow:** the pusher period is read once per run and cached. It feeds into both
the calibrant peak-picking step (converting peak bins to $t_D$) and the analyte CCS
conversion step (converting every bin to $t_D$).

### Step 2. IMS settings validation — gating which runs can use the calibration

**Where the data comes from:** `_extern.inf` in each `.raw` directory contains all
instrument parameters as key-value pairs.

**What we check:** the calibration curve is only valid for data acquired with the
same wave parameters. The pipeline compares each analyte run against the calibrant
reference:

| Setting | Mismatch action |
|---------|----------------|
| IMS Wave Velocity (m/s) | `[HARD]` — skip run |
| IMS Wave Height (V) | `[HARD]` — skip run |
| IMS Gas Flow (mL/min) | `[HARD]` — skip run |
| Trap Gas Flow (mL/min) | `[WARN]` — proceed with warning |
| Transfer Wave settings | `[WARN]` — proceed with warning |

**Why hard-skip on wave parameters:** the power-law exponent $X$ is determined by
the wave velocity and height. A calibration at 400 m/s / 25 V does not apply to data
at 650 m/s / 40 V — the exponent would be completely different.

**Why soft-warn on trap gas:** trap gas flow affects ion activation before the IM cell
but does not directly change the IM separation physics. Differences may indicate
slightly different conformer populations, but the calibration curve still applies.

**Data flow:** settings are read from `.raw` dirs via `--raw-dir`. Also extracts the
**EDC Delay Coefficient** $C$ used in Step 4.

### Step 3. Peak picking — matching calibrant charge states to their drift times

**Where the data comes from:**
- `ccs_calibrant_table.json` — Table 2 from the paper: for each protein, lists the
  literature CCS ($\Omega_{\text{lit}}$) at each charge state $z$
- `_im.txt` — the 2D IM data as $(m/z,\;\text{drift\_bin},\;\text{intensity})$ triples

**What we compute:** for each (protein, $z$) entry in the table:

1. Calculate the expected $m/z$:

$$m/z = \frac{M_W + z \cdot m_H}{z}, \quad m_H = 1.00728\;\text{Da}$$

> **Example (Myoglobin $z = 15$):**
> $m/z = \frac{16951.0 + 15 \times 1.00728}{15} = 1131.07$ Da

2. Extract a 1D drift profile by summing all 2D data within $\pm 2$ Da of target $m/z$:

$$I(i) = \sum_{\substack{j \;\text{where}\\ |m/z_j - m/z_{\text{target}}| \leq 2}} I(j, i)$$

3. Find the peak drift time as the intensity-weighted centroid ($\pm 5$ bins around
   the argmax). Reject flat/noisy profiles where peak/median $< 3$.

**Why centroid, not argmax:** with only ~5-10 drift bins spanning a peak, the argmax
has $\pm 0.5$ bin discretization error. The centroid gives sub-bin precision, which
matters because the CCS calibration amplifies small drift-time differences through the
power law.

**Data flow:** each successfully picked peak becomes a row in `calibrant_peaks.csv`
with fields `(protein, z, m/z, MW, drift_bin, drift_time_ms, peak_intensity,
literature_ccs)`. This table is the input to the fitting step.

### Step 4. EDC correction — removing the mass-dependent transfer delay

**Where the data comes from:** the EDC Delay Coefficient $C$ is read from
`_extern.inf` in Step 2.

**What we compute:**

$$t'_D = t_D - \frac{C \cdot \sqrt{m/z}}{1000}$$

**Why:** the SYNAPT's transfer optics introduce a mass-dependent delay before ions
reach the IM cell. Heavier ions (higher $m/z$) take longer to transfer. This delay
is *not* part of the IM separation and must be subtracted. The $\sqrt{m/z}$ dependence
comes from the time-of-flight relationship in the transfer region.

> **Example (MYO $z = 15$, $C = 1.58$):**
>
> $t'_D = 7.3184 - \frac{1.58 \times \sqrt{1131.07}}{1000} = 7.3184 - 0.0531 = 7.2653$ ms
>
> The correction is small (~0.05 ms) but systematic: without it, high-$m/z$ species
> appear to have inflated CCS.

**Data flow:** $t'_D$ is computed for every calibrant peak (fed to fitting) and later
for every analyte drift bin (fed to CCS conversion).

### Step 5. Reduced cross-section — factoring out charge and mass

**What we compute:**

$$\mu = \sqrt{\frac{1}{M_W} + \frac{1}{M_{\text{gas}}}}, \quad M_{\text{gas}} = 28.014\;\text{Da (N}_2\text{)}$$

$$\Omega' = \frac{\Omega_{\text{lit}}}{z \cdot \mu}$$

**Why:** different calibrant species have different molecular weights and charge states.
To discover the universal power law $\Omega' \propto (t'_D)^X$, we need to normalize
out the ion-specific factors. The reduced cross-section $\Omega'$ removes:
- **Charge $z$** — higher charge means higher drift velocity at the same CCS, so ions
  arrive earlier. Dividing by $z$ corrects for this.
- **Reduced mass $\mu$** — the collision dynamics depend on the ion-gas mass ratio.
  A 17 kDa myoglobin ion collides differently with N$_2$ than an 8.6 kDa ubiquitin ion.
  The $\mu$ factor normalizes this.

After this normalization, data from UBQ, CYTC, and MYO should all fall on the *same*
power-law curve — that's the whole point of the two-step approach.

> **Example (MYO $z = 15$, $\Omega_{\text{lit}} = 3080\;\text{\AA}^2$):**
>
> $\mu = \sqrt{\frac{1}{16951.0} + \frac{1}{28.014}} = 0.18909$
>
> $\Omega' = \frac{3080}{15 \times 0.18909} = 1085.90$

**Data flow:** $\Omega'$ is computed only for calibrant peaks (not analytes). It pairs
with $t'_D$ as input to Fit 1.

### Step 6. Fit 1 — Power-law regression: discovering the instrument exponent

**What we fit:**

$$\ln(\Omega') = X \cdot \ln(t'_D) + \ln(A)$$

which is the log-linearized form of:

$$\Omega' = A \cdot (t'_D)^X$$

**Why this is Fit 1 and not the final answer:** this fit operates on *reduced*
quantities ($\Omega'$, not $\Omega$). Its purpose is solely to discover the
instrument-specific exponent $X$, which encodes how the travelling wave's nonlinear
separation maps drift time to mobility. Once $X$ is known, Fit 2 uses it to build
the actual CCS conversion.

**Iterative outlier rejection:** under denaturing conditions, TWIMS can exhibit
"rollover" — high-charge-state ions have high mobility and can surf the wave, arriving
at drift times that don't follow the power law. These points would corrupt $X$.
The fitter:
1. Fits all points
2. Removes the point with the largest $|\text{residual}|$ in $\ln$-$\ln$ space
3. Refits
4. Repeats until $R^2 \geq 0.98$ or only 5 points remain
5. **Re-inclusion pass:** tries adding each excluded point back — if $R^2$ stays
   above threshold, the point is restored

> **Example (MYO $z = 15$) — one point on the regression:**
>
> $\ln(t'_D) = \ln(7.2653) = 1.9831$
>
> $\ln(\Omega') = \ln(1085.90) = 6.9902$

> **Result from this dataset:** 32 peaks picked $\to$ 23 excluded (TWIMS rollover)
> $\to$ **9 retained**, $X = 0.2980$, $A = 594.55$, $R^2 = 0.987$

**Data flow:** $X$ and $A$ are stored in `calibration_curve.json` and fed to Fit 2.

### Step 7. Fit 2 — Linear regression: converting to absolute CCS

**What we compute and fit:**

$$t''_D = (t'_D)^X \cdot z \cdot \mu$$

$$\Omega_{\text{lit}} = \text{slope} \cdot t''_D + \text{intercept}$$

**Why a second fit:** Fit 1 gave us $X$ but operates in reduced space. The variable
$t''_D$ re-introduces the charge $z$ and mass factor $\mu$ that were divided out in
Step 5, using the now-known power law to transform $t'_D$. The result is a quantity
that is linearly proportional to $\Omega_{\text{lit}}$ across all calibrant species —
different proteins, different charge states, all on one line.

The slope and intercept of this line are the **final calibration constants**. To convert
any measured drift time to CCS, you only need: $t_D$, $m/z$, $z$, $M_W$, and these
constants.

> **Example (MYO $z = 15$):**
>
> $t''_D = (7.2653)^{0.2980} \times 15 \times 0.18909 = 5.1222$
>
> $\Omega_{\text{cal}} = 573.02 \times 5.1222 + 122.83 = 3057.9\;\text{\AA}^2$
>
> $\text{Error} = \frac{3057.9 - 3080}{3080} \times 100 = -0.72\%$

> **Result from this dataset:** slope $= 573.02$, intercept $= 122.83$, $R^2 = 0.999$

**Data flow:** slope and intercept are stored in `calibration_curve.json`. Together
with $X$ and $C$, they define the complete calibration. The `calibrant_summary.csv`
applies this calibration back to *every* picked peak (including excluded ones) so you
can audit the error for each charge state.

### Step 8. Analyte CCS conversion — applying the calibration bin-by-bin

**Where the data comes from:**
- `adp_analytes.csv` — species list with $m/z$, $z$, $M_W$ for each analyte
- `_im.txt` — the analyte run's 2D IM data (same format as calibrants)
- `calibration_curve.json` — the fitted $X$, slope, intercept, $C$

**What we compute for each drift bin $i$:**

$$\boxed{\Omega(i) = \text{slope} \cdot \left[\left(i \cdot \frac{\Delta t}{1000} - \frac{C \cdot \sqrt{m/z}}{1000}\right)^{\!X} \cdot z \cdot \sqrt{\frac{1}{M_W} + \frac{1}{M_{\text{gas}}}}\;\right] + \text{intercept}}$$

This is the same chain as Steps 4 + 7, applied to every bin instead of just peak bins.

**Which runs are processed:** only genuine analyte runs. The pipeline automatically
excludes:
- Calibrant runs specified in `--calibrants`
- Protein standard runs (names containing UBQ/MYO/CYTC/CALI)
- Runs with `[HARD]` IMS setting mismatches

> **Worked example: $[\text{ADP}+\text{H}]^+$ at drift bin 31**
>
> | Parameter | Value | Source |
> |-----------|-------|--------|
> | $m/z$ | 428.21 | `adp_analytes.csv` |
> | $z$ | 1 | `adp_analytes.csv` |
> | $M_W$ | 427.20 Da | `adp_analytes.csv` |
> | $C$ | 1.58 | `_extern.inf` |
> | $\Delta t$ | 110 $\mu$s | `_FUNC001.STS` |
> | $X$ | 0.2980 | Fit 1 |
> | slope | 573.02 | Fit 2 |
> | intercept | 122.83 | Fit 2 |
>
> **1. Drift time:**
> $t_D = 31 \times 110 / 1000 = 3.4100$ ms
>
> **2. EDC correction:**
> $t'_D = 3.4100 - 1.58 \times \sqrt{428.21} / 1000 = 3.4100 - 0.0327 = 3.3773$ ms
>
> **3. Reduced-mass factor:**
> $\mu = \sqrt{1/427.20 + 1/28.014} = 0.19503$
>
> **4. Linearized drift time:**
> $t''_D = (3.3773)^{0.2980} \times 1 \times 0.19503 = 0.28031$
>
> **5. CCS:**
> $\Omega = 573.02 \times 0.28031 + 122.83 = \mathbf{283.4\;\text{\AA}^2}$

**Data flow:** the full bin-by-bin conversion is saved to `csv_raw/` with columns
`(species, mw, z, mz, drift_bin, drift_time_ms, ccs, intensity)`. This is the
native-resolution output with non-uniform CCS spacing.

### Step 9. Resampling and smoothing — uniform CCS grid for analysis

**Why resample:** the mapping $i \to \Omega$ is nonlinear (power law), so native bins
are not evenly spaced in CCS. Near the ADP peak, spacing is $\approx 1.5\;\text{\AA}^2$/bin,
but at low drift times it exceeds $15\;\text{\AA}^2$/bin. If you smooth in bin space,
a Gaussian with $\sigma = 3$ bins covers $\sim 4.5\;\text{\AA}^2$ near the peak but
$\sim 45\;\text{\AA}^2$ at the edges — the smoothing width is spatially non-uniform in CCS.

**What we do:**

1. **Resample** the raw $({\Omega_i, I_i})$ profile via linear interpolation onto a
   uniform grid:

   $$\Omega_k = k \times \Delta\Omega, \quad \Delta\Omega = 0.5\;\text{\AA}^2$$

   Grid points are anchored at multiples of $0.5$ so that **all species share aligned
   grid points** in overlapping CCS regions — this enables direct comparison, difference
   spectra, etc.

2. **Smooth** with a Gaussian kernel ($\sigma = 3.0$ grid points $= 1.5\;\text{\AA}^2$)
   on the uniform grid. Now the smoothing width is uniform in CCS space everywhere.

**Data flow:** the resampled + smoothed profile is saved to `csv_smoothed/`. Three
styles of PNG are generated: `png_raw/` (gray, native bins), `png_smoothed/` (blue,
uniform grid), `png_overlay/` (both together). All PNGs auto-zoom to the significant
intensity region with 15% padding.

> **Note on the HTML drift-time viewer:** the interactive viewer smooths in
> drift-bin-index space. This is correct there because bins are uniformly spaced in
> $t_D$ (each bin = one pusher period). The CCS pipeline uses the resample-then-smooth
> approach instead because the CCS mapping is nonlinear.

---

## Prerequisites

1. **Calibrant protein runs** acquired on the same instrument with matching IMS
   wave parameters as the analyte runs. Supported calibrants:
   - Ubiquitin (UBQ) — $z = 6$–$11$
   - Cytochrome c (CYTC) — $z = 7$–$18$
   - Myoglobin (MYO) — $z = 8$–$22$

2. **Converted text files** (`_im.txt`, `_ms.txt`) from the Waters `.raw` data.
   If you haven't converted yet, run the main pipeline first:
   ```bash
   deconvovo -i data_2_waters/20260216 -o output/waters_20260216
   ```

3. **Analyte species list** (CSV) with columns: `species,mz,z,mw`.
   See `scripts/adp_analytes.csv` for an example.

---

## Running the Calibration

### Standalone

```bash
python -m deconvovo.imms_ccs_calibrate \
  -i output/waters_20260216/_converted \
  -o output/ccs_calibration \
  --calibrants "ubiquitin:20260216_UBQ_DEN_NEW_V1,myoglobin:20260216_MYO_DEN_NEW_V1,cytochrome_c:20260216_CYTC_DEN_TUNE_V1" \
  --analytes scripts/adp_analytes.csv \
  --raw-dir data_2_waters/20260216
```

### From the main pipeline CLI

```bash
deconvovo -i data_2_waters/20260216 -o output/results \
  --ccs-calibrants "ubiquitin:20260216_UBQ_DEN_NEW_V1,myoglobin:20260216_MYO_DEN_NEW_V1,cytochrome_c:20260216_CYTC_DEN_TUNE_V1" \
  --analyte-csv scripts/adp_analytes.csv \
  --raw-dir data_2_waters/20260216
```

**Arguments:**

| Arg | Description |
|-----|-------------|
| `-i` | Directory containing `_im.txt` / `_ms.txt` files |
| `-o` | Output directory (created if needed) |
| `--calibrants` | Comma-separated `protein:run_name` pairs |
| `--analytes` | Path to analyte species CSV |
| `--raw-dir` | Path to `.raw` directories (reads IMS settings + pusher period) |
| `--min-r2` | Minimum $R^2$ for calibration fits (default: 0.98) |
| `--mz-tol` | $m/z$ tolerance for drift profile extraction (default: 2.0 Da) |

---

## Output Structure

```
output/ccs_calibration/
├── calibrant_peaks.csv          # All detected calibrant peaks
├── calibrant_summary.csv        # CCS_literature vs CCS_calibrated for every peak
├── calibration_curve.json       # Fit parameters (X, A, slope, intercept, R²)
├── calibration_plots.png        # 2-panel: ln-ln + linear fit (used ● + excluded ✕)
├── warnings.log                 # All warnings (settings mismatches, excluded peaks)
├── csv_raw/                     # Per-species raw CCS profiles (native drift bins)
│   └── {run}_{species}_ccs.csv
├── csv_smoothed/                # Per-species smoothed CCS profiles (uniform grid)
│   └── {run}_{species}_ccs.csv
├── png_raw/                     # Per-species raw profile plots (gray)
│   └── {run}_{species}_ccs.png
├── png_smoothed/                # Per-species smoothed profile plots (blue line)
│   └── {run}_{species}_ccs.png
└── png_overlay/                 # Per-species raw + smoothed overlay plots
    └── {run}_{species}_ccs.png
```

### `csv_raw` columns

| Column | Description |
|--------|-------------|
| `species` | Species label from analyte CSV |
| `mw` | Molecular weight $M_W$ (Da) |
| `z` | Charge state |
| `mz` | $m/z$ used for extraction |
| `drift_bin` | Native pusher bin index $i$ |
| `drift_time_ms` | $t_D = i \times \Delta t / 1000$ (ms) |
| `ccs` | $\Omega$ from calibration curve ($\text{\AA}^2$) — non-uniform spacing |
| `intensity` | Raw intensity at this bin |

### `csv_smoothed` columns

| Column | Description |
|--------|-------------|
| `species` | Species label from analyte CSV |
| `mw` | Molecular weight $M_W$ (Da) |
| `z` | Charge state |
| `mz` | $m/z$ used for extraction |
| `ccs` | $\Omega$ on uniform $0.5\;\text{\AA}^2$ grid — aligned across all species |
| `intensity` | Gaussian-smoothed intensity ($\sigma = 3.0$ in CCS space) |

### `calibrant_summary.csv` columns

| Column | Description |
|--------|-------------|
| `protein` | Calibrant protein name |
| `mw` | $M_W$ (Da) |
| `z` | Charge state |
| `mz` | Expected $m/z$ |
| `drift_time_ms` | Measured $t_D$ (ms) |
| `t_prime` | EDC-corrected $t'_D$ (ms) |
| `CCS_literature` | $\Omega_{\text{lit}}$ from Table 2 ($\text{\AA}^2$) |
| `CCS_calibrated` | $\Omega_{\text{cal}}$ from the calibration curve ($\text{\AA}^2$) |
| `CCS_error_pct` | $(\Omega_{\text{cal}} - \Omega_{\text{lit}}) / \Omega_{\text{lit}} \times 100$ (%) |
| `used_in_fit` | Whether this point was used in the final fit |

---

## Assessing Calibration Quality

1. **Check `calibration_plots.png`** — both panels should show a clear linear trend.
   Blue circles = used in fit. Red $\times$ = excluded (typically TWIMS rollover).

2. **Check `calibrant_summary.csv`** — for points with `used_in_fit=True`,
   `CCS_error_pct` should be $< 5\%$. Large errors on excluded points are expected.

   > **Note:** a point can have low `CCS_error_pct` but still be excluded if it
   > doesn't fit the $\ln$-$\ln$ power law. This happens because the two-step
   > calibration can compensate in the linear step: a point far from the power-law
   > line in $(\ln t'_D,\;\ln\Omega')$ space may still land near the correct CCS after
   > the linear transform. But including such a point would distort the power-law
   > exponent $X$, degrading accuracy for other points.

3. **Check $R^2$ values** in `calibration_curve.json`:
   - `r2_lnln` $\geq 0.98$ ($\ln$-$\ln$ fit)
   - `r2_linear` $\geq 0.98$ (linear fit)

4. **Check `warnings.log`** for:
   - `[HARD]` — critical mismatches (runs were skipped)
   - `[WARN]` — non-critical differences (runs processed with caveat)
   - `Excluded from fit` — charge states removed by outlier rejection

---

## Smoothing Library

CCS profiles are smoothed using `deconvovo.smooth`, a reusable Python library that
provides the same smoothing methods as the interactive HTML 2D drift-time viewer:

```python
from deconvovo.smooth import smooth1d

# Gaussian smoothing
y_smooth = smooth1d(y, method="gaussian", sigma=3.0)

# Savitzky-Golay (preserves peak shape)
y_smooth = smooth1d(y, method="sg", window=7, polyorder=2)

# Moving average
y_smooth = smooth1d(y, method="ma", window=5)

# With noise floor (zero out below 1% of max, then smooth)
y_smooth = smooth1d(y, method="gaussian", sigma=3.0, noise_pct=1.0)
```

---

## Troubleshooting

| Issue | Solution |
|-------|----------|
| $R^2$ too low after outlier rejection | Too few well-resolved charge states. Try different calibrant runs or lower `--min-r2` |
| All UBQ/CYTC points excluded | TWIMS rollover — check `calibrant_summary.csv` for `CCS_error_pct` |
| CCS values seem wrong | If `CCS_error_pct` $> 10\%$ for used points, the calibration is unreliable |
| "Cannot find pusher period" | Provide `--raw-dir` pointing to the `.raw` directories |
| No analyte profiles generated | Check that `_im.txt` files exist and aren't filtered by protein-name tags or `[HARD]` warnings |
| Analyte run skipped with `[HARD]` | Critical IMS mismatch — calibration does not apply to this run |

---

## References

- Ruotolo BT, Benesch JLP, Sandercock AM, Hyung SJ, Robinson CV. Ion mobility–mass
  spectrometry analysis of large protein complexes. *Nature Protocols* 3(7):1139–1152, 2008.
- Literature CCS values: `scripts/ccs_calibrant_table.json` (Table 2 from the paper,
  excluding metastable entries)
