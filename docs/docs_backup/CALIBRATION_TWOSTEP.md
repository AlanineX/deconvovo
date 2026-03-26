# Twostep Method — Stepwise Derivation

CCS calibration for travelling-wave IMS using two sequential fits,
as described in Ruotolo & Robinson (*Nature Protocols*, 2008), Steps 11-16.

---

## Steps 1-5: identical to the direct method

The twostep method shares the same first five steps as the direct method.
See `CALIBRATION_DIRECT.md` for the full derivation of:

- **Step 1:** TWIMS transit time $t_D \propto 1/K^n$ (Shvartsburg & Smith 2008)
- **Step 2:** Substituting Mason-Schamp: $\Omega \propto (z/\sqrt{\mu}) \cdot t_D^{1/n}$
- **Step 3:** EDC correction: $t'_D = t_D - C \sqrt{m/z}/1000$
- **Step 4:** Reduced cross-section: $\Omega' = \Omega / (z \cdot \mu_f)$ where
  $\mu_f = \sqrt{1/M_W + 1/M_{\text{gas}}}$
- **Step 5:** Power-law fit: $\ln(\Omega') = X \cdot \ln(t'_D) + \ln(A)$, giving $X$ and $A$

At this point, both methods have the same $X$ and $A$ from the same Fit 1.

---

## Step 6: Construct $t''_D$ (where the methods diverge)

The direct method proceeds to compute CCS as $\Omega = A \cdot (t'_D)^X \cdot z \cdot \mu_f$.
The twostep method takes a different path.

**Starting from the direct method equation:**

$$\Omega = A \cdot (t'_D)^X \cdot z \cdot \mu_f \tag{from Direct Eq. 13}$$

Ruotolo groups the measurement-dependent and ion-dependent terms:

$$\Omega = A \cdot \underbrace{(t'_D)^X \cdot z \cdot \mu_f}_{\equiv\; t''_D} \tag{14}$$

> **Why this grouping:** the factor $(t'_D)^X$ comes from the measurement, while
> $z \cdot \mu_f$ comes from the ion's identity. The product $t''_D$ combines them
> into a single variable that should be linearly proportional to $\Omega$.

**Definition:**

$$t''_D \equiv (t'_D)^X \cdot z \cdot \sqrt{\frac{1}{M_W} + \frac{1}{M_{\text{gas}}}} \tag{15}$$

> **Source:** Ruotolo et al. 2008, Step 15. The paper calls this the
> "corrected drift time" though it is not a physical time — it is a
> dimensionless composite variable.

If the power law (Eq. 10) were perfect, then $\Omega = A \cdot t''_D$ exactly —
a straight line through the origin with slope $A$.

---

## Step 7: Fit 2 — linear regression of $\Omega_{\text{lit}}$ vs $t''_D$

**The key insight of the twostep method:** instead of assuming the relationship
$\Omega = A \cdot t''_D$ holds exactly, allow a **free slope and intercept**:

$$\Omega_{\text{lit}} = \text{slope} \cdot t''_D + \text{intercept} \tag{16}$$

> **Source:** Ruotolo et al. 2008, Step 16.
> **Why a second fit:** the power law $\Omega' = A \cdot (t'_D)^X$ is an
> approximation. In practice, it may systematically over- or under-predict
> at certain drift times due to waveform non-ideality, velocity relaxation,
> or other effects not captured by a single power-law exponent. The intercept
> absorbs this systematic offset within the calibrant range.

**Fitting procedure:** For each calibrant charge state, compute $t''_D$ from
Eq. (15) using the $X$ obtained in Step 5. Pair it with the known
$\Omega_{\text{lit}}$ from the literature. Perform ordinary least-squares
linear regression on all $(t''_D, \Omega_{\text{lit}})$ pairs.

This gives:
- **slope** — the effective proportionality constant (replaces $A$)
- **intercept** — the systematic offset correction
- **$R^2$** — quality of the linear fit

If the power law were perfect: slope $= A$, intercept $= 0$.

---

## Step 8: Convert an unknown analyte to CCS

For an unknown analyte with known $m/z$, $z$, and $M_W$:

1. Compute $t'_D$ from the measured drift time:
   $$t'_D = t_D - \frac{C \cdot \sqrt{m/z}}{1000} \tag{Eq. 6, same as direct}$$

2. Compute $t''_D$ using $X$ from Fit 1:
   $$t''_D = (t'_D)^X \cdot z \cdot \sqrt{\frac{1}{M_W} + \frac{1}{M_{\text{gas}}}} \tag{Eq. 15}$$

3. Apply the Fit 2 linear equation:
   $$\boxed{\Omega = \text{slope} \cdot t''_D + \text{intercept}} \tag{17}$$

This is the complete twostep method. Three parameters ($X$ from Fit 1, slope
and intercept from Fit 2) from two fits.

---

## Comparison with the direct method at Step 6

The direct method's conversion is:

$$\Omega_{\text{direct}} = A \cdot t''_D \tag{from Eq. 13, rewritten}$$

The twostep method's conversion is:

$$\Omega_{\text{twostep}} = \text{slope} \cdot t''_D + \text{intercept} \tag{Eq. 17}$$

Both use the same $t''_D$. The difference is:

$$\Omega_{\text{twostep}} - \Omega_{\text{direct}} = (\text{slope} - A) \cdot t''_D + \text{intercept} \tag{18}$$

Within the calibrant range (large $t''_D$), this difference is small because
slope $\approx A$ and the intercept is a minor correction. Outside the calibrant
range (small $t''_D$ for small molecules), the intercept dominates:

$$\frac{\text{intercept}}{\Omega_{\text{twostep}}} \to 1 \quad \text{as} \quad t''_D \to 0 \tag{19}$$

This is why the twostep method is unsuitable for extrapolation to species with
CCS far below the calibrant range.

---

## Summary of data flow

```
FIT 1 (same as direct method):
  For each calibrant (protein, z):
    [Eq. 6]  t'_D = t_D − C·√(m/z)/1000
    [Eq. 8]  Ω' = Ω_lit / (z · μ_f)
    [Eq. 11] Plot ln(Ω') vs ln(t'_D) → linear regression → X, A

FIT 2 (twostep only):
  For each calibrant (protein, z):
    [Eq. 15] t''_D = (t'_D)^X · z · μ_f       ← using X from Fit 1
    [Eq. 16] Plot Ω_lit vs t''_D               ← linear regression → slope, intercept

CONVERSION:
  For each analyte species:
    [Eq. 6]  t'_D = t_D − C·√(m/z)/1000
    [Eq. 15] t''_D = (t'_D)^X · z · μ_f
    [Eq. 17] Ω = slope · t''_D + intercept     ← CCS result
```

---

## When to use the twostep method

The twostep method is appropriate when:

1. **Analyte CCS falls within the calibrant CCS range** — the intercept
   is a valid correction, not an extrapolation artifact
2. **The power law is imperfect** — e.g., systematic curvature in the
   $\ln(\Omega')$ vs $\ln(t'_D)$ plot that a single $X$ cannot capture
3. **Protein-vs-protein calibration** — calibrants and analytes have
   similar sizes, charges, and structural classes

The twostep method should NOT be used when analyte CCS is far below the
calibrant range (e.g., small molecules calibrated against proteins), because
the intercept has no physical basis at small $t''_D$ values.

---

## References

- Ruotolo BT, Benesch JLP, Sandercock AM, Hyung SJ, Robinson CV.
  *Nature Protocols* 3:1139-1152, 2008. (Steps 11-16 describe the two-step procedure)
- Shvartsburg AA, Smith RD. *Anal. Chem.* 80:9689-9699, 2008.
  [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC2761765/)
- Mason EA, McDaniel EW. *Transport Properties of Ions in Gases.* Wiley, 1988.
