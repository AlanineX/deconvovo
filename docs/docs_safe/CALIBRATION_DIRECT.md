# Direct Method — Stepwise Derivation

CCS calibration for travelling-wave IMS using a single power-law fit.

---

## Starting Point: the Mason-Schamp equation

From kinetic theory of ion-gas collisions (Mason & McDaniel, 1988), the ion mobility
$K$ relates to the collision cross-section $\Omega$ by:

$$K = \frac{3\,z\,e}{16\,N_0} \cdot \sqrt{\frac{2\pi}{\mu\,k_B\,T}} \cdot \frac{1}{\Omega} \tag{1}$$

where $\mu = \frac{m_{\text{ion}} \cdot M_{\text{gas}}}{m_{\text{ion}} + M_{\text{gas}}}$
is the reduced mass of the ion-gas collision pair.

> **Source:** Mason EA, McDaniel EW. *Transport Properties of Ions in Gases.* Wiley, 1988.
> Equation (1) is the low-field limit of the kinetic theory result.

This equation applies directly in drift-tube IMS, where $K$ is measured from drift
velocity. In TWIMS, $K$ cannot be measured directly — we need a different route.

---

## Step 1: TWIMS transit time depends on mobility via a power law

Shvartsburg & Smith (2008) showed from first principles that in TWIMS, the ion transit
velocity scales as:

$$v_{\text{transit}} \propto K^n \tag{2}$$

where $n$ depends on the waveform shape and operating regime.

> **Source:** Shvartsburg AA, Smith RD. *Anal. Chem.* 80:9689, 2008.
> Equation (2) is their main theoretical result (Eq. 10 in the paper for
> triangular waves: $\bar{v} = (KE)^2/s$, giving $n = 2$).

Since transit time $t \propto 1/v$:

$$t_D \propto \frac{1}{K^n} \tag{3}$$

This is the fundamental difference from drift-tube IMS where $t_D \propto 1/K$ ($n = 1$).

---

## Step 2: Substitute the Mason-Schamp relation

From Eq. (1), at fixed temperature and pressure, $K \propto \frac{z}{\Omega \sqrt{\mu}}$.
Substituting into Eq. (3):

$$t_D \propto \left(\frac{\Omega \sqrt{\mu}}{z}\right)^n \tag{4}$$

Inverting to express $\Omega$ in terms of $t_D$:

$$\Omega \propto \frac{z}{\sqrt{\mu}} \cdot t_D^{1/n} \tag{5}$$

---

## Step 3: Correct for the transfer delay

The measured drift time includes a mass-dependent transfer delay (see §2 of
CCS_METHOD_COMPARISON.md):

$$t'_D = t_D - \frac{C \cdot \sqrt{m/z}}{1000} \tag{6}$$

> **Source:** Waters instrument documentation. $C$ is the EDC Delay Coefficient
> stored in `_extern.inf`. The $\sqrt{m/z}$ dependence comes from
> $t_{\text{transfer}} \propto \sqrt{m/z}$ in the transfer optics.

Replace $t_D$ with $t'_D$ in Eq. (5):

$$\Omega \propto \frac{z}{\sqrt{\mu}} \cdot (t'_D)^{1/n} \tag{7}$$

---

## Step 4: Define the reduced cross-section $\Omega'$

To fit a single curve across different proteins and charge states, we normalize
out the ion-specific factors $z$ and $\mu$. Define:

$$\Omega' \equiv \frac{\Omega}{z \cdot \sqrt{\frac{1}{M_W} + \frac{1}{M_{\text{gas}}}}} \tag{8}$$

Note that $\sqrt{1/M_W + 1/M_{\text{gas}}} = 1/\sqrt{\mu}$, so Eq. (8) is equivalent
to $\Omega' = \Omega \cdot \sqrt{\mu} / z$.

> **Source:** Ruotolo BT et al. *Nature Protocols* 3:1139, 2008, Step 12.
> This normalization is derived from rearranging Eq. (1) to isolate the
> mobility-dependent part.

Substituting into Eq. (7):

$$\Omega' \propto (t'_D)^{1/n} \tag{9}$$

All the ion-specific factors ($z$, $\mu$) cancel. This means $\Omega'$ depends
**only** on the drift time and the wave parameters — not on which protein or
charge state we're looking at.

---

## Step 5: Define $X$ and $A$, fit by linear regression

Set $X \equiv 1/n$ and include a proportionality constant $A$:

$$\Omega' = A \cdot (t'_D)^X \tag{10}$$

Take the natural logarithm:

$$\ln(\Omega') = X \cdot \ln(t'_D) + \ln(A) \tag{11}$$

This is a straight line in $(\ln t'_D, \ln \Omega')$ space.

> **Source:** Ruotolo et al. 2008, Steps 13-14.

**Fitting procedure:** For each calibrant charge state, we know:
- $t'_D$ from the measured drift time (Eq. 6)
- $\Omega'$ from the literature CCS and Eq. (8)

Plot all calibrant $(\ln t'_D, \ln \Omega')$ pairs and perform ordinary least-squares
linear regression. The slope gives $X$, the intercept gives $\ln(A)$.

---

## Step 6: Convert back to absolute CCS

For an unknown analyte with known $m/z$, $z$, and $M_W$:

1. Compute $t'_D$ from the measured drift time using Eq. (6)
2. Compute $\Omega'$ from the power law: $\Omega' = A \cdot (t'_D)^X$ (Eq. 10)
3. Reverse the normalization (Eq. 8) to get the absolute CCS:

$$\Omega = \Omega' \cdot z \cdot \sqrt{\frac{1}{M_W} + \frac{1}{M_{\text{gas}}}} \tag{12}$$

Combining Eqs. (10) and (12):

$$\boxed{\Omega = A \cdot (t'_D)^X \cdot z \cdot \sqrt{\frac{1}{M_W} + \frac{1}{M_{\text{gas}}}}} \tag{13}$$

This is the complete direct method. Two parameters ($X$, $A$) from one fit.

---

## Summary of data flow

```
For each calibrant (protein, z):
  Known: Ω_lit, MW, z, m/z
  Measured: t_D (drift time from IM data)
  ──────────────────────────────────────
  [Eq. 6]  t'_D = t_D − C·√(m/z)/1000          ← corrected drift time
  [Eq. 8]  Ω' = Ω_lit / (z · √(1/MW + 1/Mgas)) ← reduced CCS
  [Eq. 11] Plot ln(Ω') vs ln(t'_D)               ← linear regression → X, A

For each analyte species:
  Known: MW, z, m/z
  Measured: t_D
  ──────────────────────────────────────
  [Eq. 6]  t'_D = t_D − C·√(m/z)/1000
  [Eq. 13] Ω = A · (t'_D)^X · z · √(1/MW + 1/Mgas)  ← CCS result
```

---

## References

- Mason EA, McDaniel EW. *Transport Properties of Ions in Gases.* Wiley, 1988.
- Shvartsburg AA, Smith RD. Fundamentals of Traveling Wave Ion Mobility Spectrometry.
  *Anal. Chem.* 80:9689-9699, 2008. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC2761765/)
- Ruotolo BT, Benesch JLP, Sandercock AM, Hyung SJ, Robinson CV.
  *Nature Protocols* 3:1139-1152, 2008.
