# CCS Conversion Methods — Theory, Derivation, and Comparison

---

## 1. Physical Background: Why Calibration Is Needed

### Drift-tube IMS: direct CCS measurement

In a uniform-field **drift tube**, an ion drifts through a buffer gas under a constant
electric field $E$. Its drift velocity is:

$$v_d = K \cdot E$$

where $K$ is the ion mobility. The **Mason-Schamp equation** (from kinetic theory of
gas-phase ion-neutral collisions) relates $K$ directly to the collision cross-section
$\Omega$:

$$K = \frac{3\,z\,e}{16\,N_0} \cdot \frac{1}{\sqrt{\mu\,k_B\,T}} \cdot \frac{1}{\Omega}$$

where:
- $z$ = charge state, $e$ = elementary charge
- $N_0$ = buffer gas number density
- $k_B$ = Boltzmann constant, $T$ = temperature
- $\mu = \frac{m_{\text{ion}} \cdot M_{\text{gas}}}{m_{\text{ion}} + M_{\text{gas}}}$ = **reduced mass** of the ion–gas collision pair

Since $v_d = L / t_D$ (drift length / drift time), you can solve directly for $\Omega$
from measured $t_D$. No calibration needed — the physics gives a closed-form equation.

### Travelling-wave IMS (TWIMS): no closed-form equation

In TWIMS, ions are propelled by a **periodic voltage wave** that sweeps through the
mobility cell. The ion dynamics are fundamentally different — ions "surf" the wave or
fall behind it depending on their mobility relative to the wave speed.

Shvartsburg & Smith (2008) showed theoretically that in the low-mobility regime, the
ion transit velocity scales as:

$$v_{\text{transit}} \propto K^2 \cdot E^2 \quad \text{(TWIMS)}$$

compared to drift-tube:

$$v_{\text{drift}} \propto K \cdot E \quad \text{(DT-IMS)}$$

This means the transit time in TWIMS is:

$$t_{\text{transit}} \propto \frac{1}{K^2} \propto \Omega^2 \quad \text{(low-mobility limit)}$$

But this $K^2$ scaling is only exact for an idealized triangular waveform at low
mobility. For real waveforms and intermediate mobilities, the exponent deviates from 2
in a way that depends on wave velocity, wave height, gas pressure, and the waveform
profile. **There is no closed-form equation** relating $t_D$ to $\Omega$ in TWIMS.

This is why **empirical calibration** with known standards is required.

---

### Key variables used throughout this document

| Symbol | Name | Definition | Units |
|--------|------|------------|-------|
| $t_D$ | **Measured drift time** | Raw arrival time from the instrument: $t_D = \text{bin} \times \Delta t / 1000$. Includes both the IM separation and the transfer delay. | ms |
| $t'_D$ | **Corrected drift time** | $t_D$ with the mass-dependent transfer delay removed: $t'_D = t_D - C \sqrt{m/z} / 1000$. Represents only the time in the IM cell. This is the input to the power-law fit. | ms |
| $t''_D$ | **Corrected drift time** (Ruotolo's term) | $t''_D = (t'_D)^X \cdot z \cdot \mu$. Not a physical time — a dimensionless composite variable that combines the power-law-transformed drift time with charge and reduced mass. Constructed so that $\Omega \approx A \cdot t''_D$. **Used only in the twostep method** — the direct method does not need it. | (dimensionless) |
| $\Omega$ | **Collision cross-section (CCS)** | The physical property we report — the effective scattering area of the ion. | Å² |
| $\Omega'$ | **Reduced cross-section** | CCS normalized by charge and reduced mass: $\Omega' = \Omega / (z \cdot \mu)$. A temporary intermediate used only during calibration fitting to collapse all calibrant species onto one curve. | (dimensionless) |
| $\mu$ | **Reduced-mass factor** | $\mu = \sqrt{1/M_W + 1/M_{\text{gas}}}$, from the Mason-Schamp equation. | Da$^{-1/2}$ |
| $X$ | **Power-law exponent** | Slope of $\ln(\Omega')$ vs $\ln(t'_D)$. Instrument-specific, depends on wave parameters. | (dimensionless) |
| $A$ | **Power-law prefactor** | $A = e^{\text{intercept of ln-ln fit}}$. | (dimensionless) |
| $C$ | **EDC Delay Coefficient** | Instrument constant for the transfer delay correction. Typically 1.4–1.6. | ms·Da$^{-1/2}$ |

---

## 2. The EDC Correction: What It Is and Why $\sqrt{m/z}$

### Why can't the instrument just measure the IM drift time directly?

The SYNAPT records an **arrival time**: the interval from the moment the trap gate
releases a packet of ions into the IM cell, to the moment each ion hits the TOF
detector. Ions must physically travel through:

$$\text{Trap} \;\to\; \underbrace{\text{IM cell}}_{\text{separation happens here}} \;\to\; \text{Transfer cell} \;\to\; \text{TOF}$$

The instrument cannot detect the instant an ion exits the IM cell — there is no
detector there. It only knows when the ion arrives at the TOF. So the recorded
"drift time" is actually:

$$t_D^{\text{measured}} = t_{\text{IM}} + t_{\text{transfer}}$$

The transfer time $t_{\text{transfer}}$ is **not** part of the IM separation but is
baked into every measurement. This is a hardware constraint, not a design flaw — the
transfer cell is needed to guide ions from the IM cell into the TOF.

### Why the transfer time scales as $\sqrt{m/z}$

In the transfer region, ions are propelled by a travelling wave and eventually
accelerated into the TOF by a fixed voltage $V$. From energy conservation:

$$\frac{1}{2} m v^2 = z e V \quad \Rightarrow \quad v = \sqrt{\frac{2\,e\,V}{m/z}}$$

The transit time over a fixed path length $d$:

$$t_{\text{transfer}} = \frac{d}{v} \propto \sqrt{m/z}$$

Heavier ions (higher $m/z$) travel slower through the transfer optics. This is the
same physics as a time-of-flight mass analyzer.

### Why the correction doesn't cancel out in fitting

Calibrant charge states span a wide $m/z$ range (e.g., 680–2120 Da for our data).
Each charge state has a **different** transfer delay:

| Charge state | $m/z$ | $\sqrt{m/z}$ | Transfer delay |
|---|---|---|---|
| CYTC $z = 18$ | 687 | 26.2 | small |
| MYO $z = 13$ | 1305 | 36.1 | medium |
| MYO $z = 8$ | 2120 | 46.0 | large |

If you don't subtract the transfer delay, MYO $z = 8$ appears to have a ~0.03 ms
longer "drift time" than it actually does in the IM cell, while CYTC $z = 18$ has
almost no extra delay. This **differential** error distorts the power-law fit because
it's not constant — it varies with $m/z$ across your calibrant points. The correction
removes this $m/z$-dependent bias so the remaining $t'_D$ reflects only the IM
separation.

### The correction formula

$$t'_D = t_D - \frac{C \cdot \sqrt{m/z}}{1000}$$

where $C$ is the **Enhanced Duty Cycle (EDC) Delay Coefficient**, an
instrument-specific constant (typically 1.4–1.6) that absorbs the proportionality
constants ($d$, $V$, geometry, etc.). It is stored in `_extern.inf` and read
automatically by the pipeline. The name "Enhanced Duty Cycle" comes from a Waters
firmware feature that uses this same $\sqrt{m/z}$ relationship to synchronize the TOF
pusher with ion arrivals, improving sensitivity — the delay coefficient was already
characterized for this purpose, so it is reused for the IM correction.

> **Concrete example** (MYO $z = 15$, $m/z = 1131.14$, $C = 1.58$):
>
> $$t'_D = 7.418 - \frac{1.58 \times \sqrt{1131.14}}{1000} = 7.418 - 0.053 = 7.365\;\text{ms}$$
>
> The correction is small (~0.7%) but varies systematically across charge states.

---

## 3. The Reduced Cross-Section: Why Divide by $z \cdot \mu$

### Terminology: $\Omega$ vs $\Omega'$

- $\Omega$ = **collision cross-section (CCS)** in Å². This is the physical property
  we want to report. It represents the effective scattering area of the ion when
  colliding with buffer gas molecules.
- $\Omega'$ = a **charge-and-mass-normalized** intermediate variable used *only*
  during calibration. It has no standard name — the Ruotolo & Robinson paper calls
  it the "corrected cross-section." It is **not** the real CCS — $\Omega$ is. We
  compute $\Omega'$ temporarily to collapse all calibrant species onto a single
  curve, then convert back to $\Omega$ at the end.

### Why we need $\Omega'$: the calibrant diversity problem

Our calibrants are three different proteins at different charge states:
- UBQ at 8.6 kDa, $z = 8$–$11$
- CYTC at 12.4 kDa, $z = 10$–$18$
- MYO at 17.0 kDa, $z = 11$–$22$

A MYO ion at $z = 15$ and a UBQ ion at $z = 9$ have very different masses, charges,
and CCS values. To fit a single power law $\Omega' = A \cdot (t'_D)^X$ across **all**
of them, we must first normalize out the ion-specific factors (charge and mass) so
that only the shape/size dependence remains.

### Origin from the Mason-Schamp equation

The Mason-Schamp equation (kinetic theory of ion-gas collisions) relates CCS to
mobility:

$$\Omega = \frac{3\,z\,e}{16\,N_0\,K} \cdot \sqrt{\frac{2\pi}{\mu\,k_B\,T}}$$

where $\mu = \frac{m_{\text{ion}} \cdot M_{\text{gas}}}{m_{\text{ion}} + M_{\text{gas}}}$
is the **reduced mass** of the ion-gas collision pair. This equation shows that
$\Omega$ depends on three ion-specific quantities:

- $z$ (charge) — higher charge → higher drift velocity at the same CCS
- $\mu$ (reduced mass) — different ion masses → different collision dynamics
- $K$ (mobility) — what the IM cell actually separates by

To isolate the mobility-dependent part, we divide out $z$ and the mass factor:

$$\Omega' = \frac{\Omega}{z \cdot \sqrt{\frac{1}{M_W} + \frac{1}{M_{\text{gas}}}}}$$

### Why this specific form?

The reduced mass is $\mu = \frac{m \cdot M}{m + M}$, which gives:

$$\frac{1}{\sqrt{\mu}} = \sqrt{\frac{m + M}{m \cdot M}} = \sqrt{\frac{1}{m} + \frac{1}{M}}$$

So the factor $\sqrt{1/M_W + 1/M_{\text{gas}}}$ is simply $1/\sqrt{\mu_{\text{reduced}}}$.
It appears naturally when rearranging the Mason-Schamp equation to isolate $K$.

After this normalization, data from UBQ, CYTC, and MYO at all charge states should
fall on a **single curve** of $\Omega'$ vs $t'_D$ — that's the whole point.

> **Concrete example** (MYO $z = 15$, $\Omega_{\text{lit}} = 3230\;\text{\AA}^2$):
>
> $$\mu_{\text{factor}} = \sqrt{\frac{1}{16952} + \frac{1}{28.014}} = 0.18909$$
>
> $$\Omega' = \frac{3230}{15 \times 0.18909} = 1138.8$$
>
> This $\Omega'$ value can now be plotted alongside UBQ and CYTC charge states on the
> same calibration curve.

---

## 4. The Power Law: Why $\Omega' = A \cdot (t'_D)^X$

### Theoretical basis: why not just $n = 2$?

Shvartsburg & Smith (2008) derived from first principles that in TWIMS, the ion
transit velocity scales as:

$$v_{\text{transit}} \propto K^n \cdot E^n$$

For an **ideal triangular waveform** at the **low-mobility limit**, $n = 2$ exactly:

$$t_D \propto \frac{1}{K^2} \propto \Omega^2$$

So why not just use $\Omega = A \cdot (t'_D)^{0.5}$ (i.e., fix $X = 1/n = 0.5$)?

Because $n = 2$ is only exact under two conditions that don't hold in practice:

1. **Ideal triangular waveform** — real SYNAPT waveforms are closer to sinusoidal.
   Shvartsburg & Smith showed that for half-sinusoidal profiles, the exponent
   deviates from 2 and depends on the mobility range.
2. **Low-mobility limit** ($c \ll 1$, where $c = K E_{\max} / s$) — at intermediate
   and high $c$, the scaling becomes more gradual. In practice, protein and small-
   molecule ions span a range of $c$ values, so $n$ is not constant.

The empirical approach: treat $n$ (or equivalently $X = 1/n$) as a **free parameter**
determined by fitting against calibrant standards. This adapts the calibration to the
actual instrument conditions.

### How the fit works: what calibrant data goes in

For each calibrant charge state, we have a pair of known values:
- $t'_D$ — measured drift time (EDC-corrected, from the IM data)
- $\Omega'$ — literature CCS normalized by $z$ and $\mu$ (from Table 2)

We plot $\ln(\Omega')$ vs $\ln(t'_D)$ for all calibrant points. If the power law
$\Omega' = A \cdot (t'_D)^X$ holds, this should be a straight line:

$$\ln(\Omega') = X \cdot \ln(t'_D) + \ln(A)$$

Linear regression on these $(\ln t'_D,\; \ln\Omega')$ pairs gives:
- **slope** $= X$ (the power-law exponent)
- **intercept** $= \ln(A)$ (so $A = e^{\text{intercept}}$)
- **$R^2$** — how well the power law fits the data

### What $X$ means physically

$X = 1/n$ where $n$ is the effective mobility-to-transit-time exponent. It encodes
how the travelling wave's nonlinear separation maps drift time to mobility for a
specific set of wave parameters:

| $X$ value | $n = 1/X$ | Regime |
|-----------|-----------|--------|
| $X = 1.0$ | $n = 1$ | Drift-tube-like (linear; never exactly in TWIMS) |
| $X = 0.5$ | $n = 2$ | Ideal triangular wave, low mobility limit |
| $X \approx 0.6$–$0.9$ | $n \approx 1.1$–$1.7$ | Typical SYNAPT operating range |

$X$ is **instrument-specific** — it changes with wave velocity, wave height, and gas
pressure. It does **not** measure resolution or separation quality; it describes the
*shape* of the drift-time-to-CCS mapping. A different wave height gives a different
$X$ but may have the same resolution.

$X$ is identical for Direct method and Twostep method — both use the same Fit 1. It only depends
on the calibrant data and the regression, not on how the result is applied.

> **From our dataset:**
> - Auto mode (WH=25V): $X = 0.697$ (23 calibrant points after outlier rejection)
> - Unfolded mode (WH=25V): $X = 0.829$ (25 unfolded charge states, no rejection)
> - Excel reference (WH=18V): $X = 0.715$ (different instrument settings)

---

## 5. Direct Method — Power-Law (our pipeline)

### Formula

After Fit 1 gives $X$ and $A$, convert any drift time to CCS:

$$\Omega' = A \cdot (t'_D)^X$$

$$\boxed{\Omega = A \cdot (t'_D)^X \cdot z \cdot \sqrt{\frac{1}{M_W} + \frac{1}{M_{\text{gas}}}}}$$

**2 parameters** ($X$, $A$) from **1 fit**. No intercept.

### Worked example: $[\text{ADP}+\text{H}]^+$ at peak (bin 49, pusher 69 μs)

$$t_D = 49 \times \frac{69}{1000} = 3.381\;\text{ms}$$

$$t'_D = 3.381 - \frac{1.58 \times \sqrt{428.008}}{1000} = 3.381 - 0.0327 = 3.348\;\text{ms}$$

$$\Omega' = 227.31 \times (3.348)^{0.8287} = 227.31 \times 2.734 = 621.5$$

$$\mu = \sqrt{\frac{1}{427} + \frac{1}{28.014}} = 0.1951$$

$$\Omega = 621.5 \times 1 \times 0.1951 = \mathbf{121.3\;\text{\AA}^2}$$

### Properties

- **No intercept**: $\Omega \to 0$ as $t_D \to 0$ (physically correct)
- **Smooth extrapolation**: power law is well-behaved outside calibrant range
- **Same functional form** as the underlying physics ($\Omega' \propto (t'_D)^X$)

---

## 6. Twostep Method — Two-Step Linear (Ruotolo & Robinson paper)

### Rationale

Direct method converts from $\Omega'$ (normalized) back to $\Omega$ (absolute) by simply
multiplying: $\Omega = \Omega' \cdot z \cdot \mu$. This assumes the power law is
perfect. Twostep method instead runs a **second regression** to handle the conversion,
allowing a non-zero intercept to absorb systematic deviations.

### Where $t''_D$ comes from: derivation

$t''_D$ is not an independent physical quantity — it is constructed algebraically
from the Fit 1 result. The derivation:

From the direct method, we know:

$$\Omega = A \cdot (t'_D)^X \cdot z \cdot \mu$$

Ruotolo rearranges this by grouping the terms that depend on the measurement
($(t'_D)^X$) with the terms that depend on the ion ($z \cdot \mu$), defining:

$$t''_D \equiv (t'_D)^X \cdot z \cdot \mu$$

so that the direct method becomes simply $\Omega = A \cdot t''_D$. This is a linear
relationship between $\Omega$ and $t''_D$, forced through the origin with slope $A$.

The twostep method relaxes this by allowing a **free slope and intercept** instead of
assuming slope $= A$ and intercept $= 0$. The name $t''_D$ is sometimes called the
"corrected drift time" in the literature, though it is not a physical time — it is
a dimensionless composite variable that combines the power-law-transformed drift time
with the ion's charge and reduced mass.

### How the second fit works

For each calibrant point, compute $t''_D$ and pair it with the known $\Omega_{\text{lit}}$:

$$t''_D = (t'_D)^X \cdot z \cdot \sqrt{\frac{1}{M_W} + \frac{1}{M_{\text{gas}}}}$$

Plot $\Omega_{\text{lit}}$ vs $t''_D$ for all calibrant points. If the power law
were perfect, this would be a straight line through the origin with slope $= A$.
The twostep method fits:

$$\Omega_{\text{lit}} = \text{slope} \cdot t''_D + \text{intercept}$$

allowing the slope to differ from $A$ and the intercept to be nonzero. This absorbs
systematic deviations from the pure power law within the calibrant range.

### Conversion formula

$$\boxed{\Omega = \text{slope} \cdot \left[(t'_D)^X \cdot z \cdot \sqrt{\frac{1}{M_W} + \frac{1}{M_{\text{gas}}}}\right] + \text{intercept}}$$

**3 parameters** ($X$ from Fit 1, slope + intercept from Fit 2) from **2 fits**.

### Worked example: $[\text{ADP}+\text{H}]^+$ at peak (unfolded mode)

Using the same $t'_D = 3.348$ ms, with unfolded-mode Fit 2 parameters
(slope $= 217.76$, intercept $= 113.55$):

$$t''_D = (3.348)^{0.8287} \times 1 \times 0.1951 = 2.734 \times 0.1951 = 0.5334$$

$$\Omega = 217.76 \times 0.5334 + 113.55 = 116.2 + 113.6 = \mathbf{229.7\;\text{\AA}^2}$$

Compare Direct method on the same data: $\Omega = 121.3$ Å². The difference (108 Å²) comes
entirely from the intercept.

### When Twostep method works well

Within the calibrant range ($t''_D \approx 5$–$17$, CCS $\approx 1525$–$3815$ Å²),
the intercept is a minor correction. For example, at the median calibrant CCS of
~2700 Å² (MYO $z = 12$, the middle of our range):

$$\frac{113.55}{2700} \approx 4.2\%$$

The intercept shifts the predicted CCS by only ~4% — a small correction that can
improve accuracy if the power law systematically over- or under-predicts at certain
$t''_D$ values. Twostep method can outperform Direct method within this range because
the intercept absorbs real curvature that the power law misses (e.g., from waveform
non-ideality or velocity relaxation effects). This is why the Ruotolo & Robinson
paper recommends it — their use case is **protein complexes** calibrated against
proteins of similar size.

### When Twostep method fails

When the analyte $t''_D$ is far below the calibrant $t''_D$ range, the intercept
becomes a large fraction of the predicted CCS:

| Mode | Calibrant $t''_D$ range | ADP $t''_D$ | Intercept | Intercept / CCS |
|------|------------------------|-------------|-----------|-----------------|
| auto | 5.5 – 13.4 | 0.45 | +29.1 | 18.6% |
| unfolded | 7.1 – 16.7 | 0.53 | +113.5 | 49.5% |

The intercept was fitted using calibrant data where $t''_D$ is 10–30× larger than
the ADP $t''_D$. At the ADP range, the intercept has no physical basis — it is pure
extrapolation of a linear fit far outside its fitted domain. The effect is worse for
unfolded mode because the higher $X$ (0.83 vs 0.70) compresses the ADP $t''_D$ while
inflating the calibrant $t''_D$ range, increasing the intercept.

---

## 7. Comparison: All 4 Trials

We ran every combination of pick mode (auto / unfolded) and conversion method
(direct / twostep) on the same dataset (WH=25V, ADP+H$^+$ analyte):

| Trial | Pick mode | Conversion | $X$ | $R^2$ | Cal. pts | **ADP+H CCS** |
|-------|-----------|-----------|------|-------|----------|---------------|
| auto × direct | auto (outlier rejection) | direct | 0.697 | 0.991 | 23 | **128.8 Å²** |
| auto × twostep | auto (outlier rejection) | twostep | 0.697 | 0.999 | 23 | **156.3 Å²** |
| unfolded × direct | unfolded (no rejection) | direct | 0.829 | 0.919 | 25 | **120.7 Å²** |
| unfolded × twostep | unfolded (no rejection) | twostep | 0.829 | 0.971 | 25 | **229.2 Å²** |

### Observations

1. **Direct method is stable**: 121–129 Å² across pick modes (8 Å² spread, 6.3%)
2. **Twostep method is unstable**: 156–229 Å² across pick modes (73 Å² spread, 47%)
   — the intercept varies dramatically with $X$, amplifying the pick-mode dependence
3. **Within calibrant range**, both methods give ~equal accuracy (< 3% mean error)
4. **Auto mode has better $R^2$** (0.991 vs 0.919) because outlier rejection removes
   TWIMS-rollover points; but the retained set may include native conformers whose
   CCS doesn't match the extended literature values

### When to use which

| Use case | Recommended |
|----------|-------------|
| **Small-molecule analytes** (CCS far below calibrant range) | **Direct** — theoretically grounded, no extrapolation artifact |
| **Protein analytes** (CCS within calibrant range) | **Twostep** — intercept corrects systematic curvature |
| **Denatured calibrants** with known extended CCS | **Unfolded pick mode** — avoids conformational ambiguity |
| **Mixed or unknown conditions** | **Auto pick mode** — outlier rejection handles bad points |

### Error propagation

**Direct** (2 parameters, 1 fit):

$$\frac{\delta\Omega}{\Omega} \approx \sqrt{\left(\frac{\delta A}{A}\right)^2 + \left(\ln(t'_D) \cdot \delta X\right)^2}$$

**Twostep** (3 parameters, 2 fits):

$$\delta\Omega \approx \sqrt{(t''_D \cdot \delta m)^2 + (m \cdot t''_D \cdot \ln t'_D \cdot \delta X)^2 + (\delta b)^2}$$

where $m$ = slope, $b$ = intercept. The $\delta b$ term is a constant absolute error.
For large CCS ($\Omega \gg b$), $\delta b / \Omega$ is negligible. For small CCS
($\Omega \sim b$), it dominates.

---

## 8. Choice for This Pipeline

### Why the direct method is preferred (positive justification)

The direct method is not chosen "because twostep is bad" — it is preferred on its
own merits:

1. **Theoretical grounding.** The underlying TWIMS physics is a power law
   ($v \propto K^n$, Shvartsburg & Smith 2008). The direct method preserves this
   functional form exactly: $\Omega' = A \cdot (t'_D)^X$. The conversion
   $\Omega = \Omega' \cdot z \cdot \mu$ is not a fit — it is a direct algebraic
   reversal of the normalization applied in §3. No additional parameters are
   introduced.

2. **Correct boundary behavior.** As $t'_D \to 0$, $\Omega \to 0$. A species with
   zero drift time has zero cross-section. This is physically correct and holds by
   construction — the power law passes through the origin.

3. **Stability under different calibrant selections.** The direct method gives
   121–129 Å² for ADP regardless of whether we use auto or unfolded pick mode
   (6.3% spread). This stability indicates the result is robust to the choice of
   which calibrant points are included.

4. **Minimal parameter count.** Two parameters ($X$, $A$) from one fit. Each
   parameter has a clear physical meaning: $X$ is the wave-ion interaction exponent,
   $A$ is the proportionality constant.

### When twostep is appropriate

The twostep method is not wrong — it is the standard approach in the Ruotolo &
Robinson (2008) paper and has specific advantages:

- When analytes fall **within** the calibrant CCS range, the intercept corrects
  real systematic deviations from the power law (e.g., from waveform non-ideality)
- For **protein-vs-protein** calibration (similar sizes, similar charge states),
  the twostep method typically gives lower residuals within the fitted range

### Default for this pipeline

The **direct method** is the default because our primary use case is small-molecule
analytes (ADP, CCS ~ 120 Å²) calibrated against denatured proteins
(CCS ~ 1500–4000 Å²). The twostep method is available via `--conversion-method twostep`
and its results are always included in `calibrant_summary.csv` (columns `CCS_twostep`
and `err_twostep_pct`) for comparison.

Both methods share the same Fit 1 and therefore the same $X$ and $A$. They differ
only in the final conversion step.

**If your analytes are proteins with CCS in the 1000–5000 Å² range**, the twostep
method may be more appropriate — the intercept can correct real systematic deviations
within that
range.

---

## 9. Code Lineage

| Source | Method | Notes |
|--------|--------|-------|
| `g2_ccs_calibration_a.py` (old script) | **twostep** (active), direct (commented out) | Had both; used Fit 2 linear |
| Excel reference (`G2_CCS_Calib_ADP.xlsx`) | **direct** | Power-law in CCS Conv sheet |
| Ruotolo & Robinson (2008) paper | **twostep** | Describes two-step procedure |
| This pipeline (`imms_ccs_calibrate.py`) | **direct** (default) | Both available via `--conversion-method` |

---

## FAQ

### Q1. Why subtract the EDC term? If it's systematic, won't it cancel out in fitting?

No. The transfer delay is $\propto \sqrt{m/z}$, and each calibrant charge state has a
**different** $m/z$ (range 680–2120 Da in our data). MYO $z = 8$ at $m/z = 2120$ gets
a ~0.07 ms delay while CYTC $z = 18$ at $m/z = 687$ gets ~0.04 ms. This differential
bias distorts the power-law fit. The correction removes it so that only the true IM
separation time remains. *(See §2)*

### Q2. Why does the instrument include the transfer time in the measurement? Is it a design flaw?

It's a hardware constraint. The SYNAPT has no detector at the exit of the IM cell —
ions must physically travel through the transfer cell and into the TOF before they are
detected. The recorded "drift time" is the total gate-to-detector interval. There is
no way to measure the IM exit time directly. The EDC correction is the software fix
for this. *(See §2, "Why can't the instrument just measure the IM drift time
directly?")*

### Q3. What is $\Omega$ vs $\Omega'$? Is $\Omega'$ the real CCS?

No — $\Omega$ (CCS, in Å²) is the physical property we report. $\Omega'$ is a
temporary intermediate used only during calibration fitting. It normalizes out charge
($z$) and reduced mass ($\mu$) so that data from different proteins at different charge
states collapse onto a single power-law curve. After fitting, we convert back to
$\Omega$. Think of $\Omega'$ as a "standardized" CCS that removes ion-specific factors
to reveal the underlying mobility relationship. *(See §3, "Terminology")*

### Q4. Since theory says $n = 2$ for low mobility, why not just use $\Omega = A \cdot (t'_D)^{0.5}$?

Because $n = 2$ requires two idealizations that don't hold in practice: (1) a perfect
triangular waveform (real SYNAPT waves are sinusoidal), and (2) the low-mobility limit
where $c = KE_{\max}/s \ll 1$ (real proteins and small molecules span a range of $c$).
Shvartsburg & Smith (2008) showed that for realistic waveforms, the exponent deviates
from 2 in ways that depend on the specific wave parameters. Using $X$ as a free
parameter lets the calibration adapt to the actual instrument conditions. *(See §4,
"Why not just $n = 2$?")*

### Q5. What is $X$ from our fitting? Does it differ between the two methods?

$X$ comes from **Fit 1 only** — it is the slope of the $\ln(\Omega')$ vs $\ln(t'_D)$
regression. It is **identical** for Direct method and Twostep method because both use the same
Fit 1. The two methods differ only in how they apply $X$ to convert an unknown drift
time to CCS. Our values: $X = 0.697$ (auto mode), $X = 0.829$ (unfolded mode). These
differ because the two modes use different subsets of calibrant points. *(See §4,
"What $X$ means physically")*

### Q6. Does a higher $X$ mean better separation or resolution?

No. $X$ describes the *shape* of the drift-time-to-CCS mapping — the nonlinearity of
the travelling wave interaction. It does not measure how well two species are resolved.
Resolution depends on peak width vs peak separation, which is controlled by gas
pressure, path length, and wave parameters independently of $X$. A different wave
height gives a different $X$ but may have the same resolution. *(See §4)*

### Q7. How do you know the twostep method gives wrong results? Just because the number is larger?

Twostep method is not wrong in general — it is the standard approach in the Ruotolo &
Robinson paper and works well for its intended use case (protein complexes calibrated
against similar-sized proteins). The issue is specific to **extrapolation far below
the calibrant range**: the ADP $t''_D$ (~0.5) is 10–30× smaller than the calibrant
$t''_D$ range (~5–17), so the intercept (29–114 Å² depending on pick mode) becomes
18–50% of the predicted CCS. Direct method avoids this because the power law has no
intercept — it naturally goes to zero as $t'_D \to 0$. Both methods give equivalent
results within the calibrant range. *(See §6-§7)*

### Q8. Are there reported $X$ values for other instruments? How about cyclic IM?

Typical $X$ values for SYNAPT instruments range from ~0.5 to ~0.9 depending on wave
velocity, wave height, and gas pressure. Lower wave velocities and higher pressures
tend to give $X$ closer to 0.5 (more drift-tube-like). The cyclic IMS (Waters SELECT
SERIES) uses the same travelling-wave principle but with a longer effective path
length (multiple passes), so the same calibration approach applies — though the
longer path gives higher resolution without necessarily changing $X$. Drift-tube IMS
instruments (Agilent, Bruker timsTOF) do not need $X$ at all because $n = 1$ exactly
(linear $K$-to-$t_D$ relationship), and CCS is computed directly from the
Mason-Schamp equation. *(See §1 and §4)*

---

## References

- Shvartsburg AA, Smith RD. Fundamentals of Traveling Wave Ion Mobility Spectrometry.
  *Anal. Chem.* 80(24):9689-9699, 2008.
  [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC2761765/)

- Ruotolo BT, Benesch JLP, Sandercock AM, Hyung SJ, Robinson CV. Ion mobility–mass
  spectrometry analysis of large protein complexes. *Nature Protocols* 3(7):1139-1152, 2008.

- Bush MF, Hall Z, Giles K, Hoyes J, Robinson CV, Ruotolo BT. Collision Cross Sections
  of Proteins and Their Complexes: A Calibration Framework and Database for Gas-Phase
  Structural Biology. *Anal. Chem.* 82(22):9557-9565, 2010.

- Gabelica V et al. Recommendations for reporting ion mobility mass spectrometry
  measurements. *Mass Spectrom. Rev.* 38(3):291-320, 2019.

- Thalassinos K et al. Characterization of Phosphorylated Peptides Using Traveling
  Wave-Based and Drift Cell-Based Ion Mobility Mass Spectrometry. *Anal. Chem.*
  81(1):248-254, 2009.
  [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC3149990/)
