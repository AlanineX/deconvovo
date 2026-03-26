# ADP CCS Literature Validation

Comparison of our TWIMS-calibrated ADP CCS values against experimentally measured
drift-tube IMS (DTIMS) reference values from the Unified CCS Compendium.

---

## Literature Reference Values

**Source:** Unified CCS Compendium (Picache et al., 2019; McLean Group, Vanderbilt)
([Zenodo](https://zenodo.org/records/6860818),
[Paper](https://pubs.rsc.org/en/content/articlelanding/2019/sc/c8sc04396e))

**Method:** Drift-tube IMS (DTIMS), direct measurement from Mason-Schamp equation.
No calibration needed — these are the ground-truth CCS values.

**Gas:** $\text{N}_2$

### ADP (Adenosine 5'-Diphosphate, MW = 427.20 Da)

| Adduct | $m/z$ | $z$ | CCS (Å²) | SD | Source |
|--------|-------|-----|-----------|-----|--------|
| $[\text{M+H}]^+$ | 428.03 | +1 | **183.2** | 0.15 | Adenosine-3',5'-DP |
| $[\text{M+H}]^+$ | 428.03 | +1 | **187.2** | 0.46 | Adenosine-5'-DP |
| $[\text{M-H}]^-$ | 426.02 | −1 | 184.3 | 0.10 | Adenosine-3',5'-DP |
| $[\text{M-H}]^-$ | 426.02 | −1 | 184.8 | 1.56 | Adenosine-5'-DP |
| $[\text{M+Na}]^+$ | 450.02 | +1 | 188.1 | 0.21 | Adenosine-3',5'-DP |
| $[\text{M+Na}]^+$ | 450.01 | +1 | 194.2 | 1.41 | Adenosine-5'-DP |

> Note: Two isomers exist (3',5'-diphosphate and 5'-diphosphate) with slightly
> different CCS. Our ADP is adenosine 5'-diphosphate. The reference value for
> $[\text{M+H}]^+$ is **187.2 ± 0.5 Å²**.

### Related Nucleotides (for context)

| Compound | Adduct | $m/z$ | CCS (Å²) |
|----------|--------|-------|-----------|
| Adenine | $[\text{M+H}]^+$ | 136.06 | 125.6 |
| Adenosine | $[\text{M+H}]^+$ | 268.10 | 155.5 |
| AMP | $[\text{M+H}]^+$ | 348.07 | 171.5 |
| **ADP** | $[\text{M+H}]^+$ | 428.03 | **187.2** |
| ATP | $[\text{M+H}]^+$ | 508.00 | 196.7 |

The CCS increases with molecular size: adenine (126) < adenosine (156) < AMP (172) <
ADP (187) < ATP (197). Each phosphate group adds ~15–25 Å² to the cross-section.

---

## Our Pipeline Results

All trials use the same data (WH=25V, EDC=1.58, ADP+H peak at bin 49, $t_D = 3.381$ ms).

| Trial | Pick mode | Conversion | ADP+H CCS | vs Lit (187.2) |
|-------|-----------|-----------|-----------|----------------|
| auto × direct | auto (23 pts) | direct | **128.8 Å²** | **−31.2%** |
| auto × twostep | auto (23 pts) | twostep | **156.3 Å²** | **−16.5%** |
| unfolded × direct | unfolded (25 pts) | direct | **120.7 Å²** | **−35.5%** |
| unfolded × twostep | unfolded (25 pts) | twostep | **229.2 Å²** | **+22.4%** |

**None of the 4 trials reproduce the literature value (187.2 Å²).** The closest is
auto × twostep at 156.3 Å² (−16.5%).

---

## Analysis: Why the Discrepancy?

### 1. Calibrant-to-analyte CCS gap (extrapolation)

Our calibrants (denatured proteins) have CCS = 1525–3815 Å². ADP has CCS ≈ 187 Å².
This is a **20× extrapolation** below the calibrant range. The power-law exponent $X$
was determined from data in the 1500–4000 range and may not be accurate at 187.

### 2. Different instrument settings

The literature DTIMS values are measured directly (no calibration needed). Our TWIMS
values depend on the calibration curve, which depends on wave parameters. The
calibration is only as good as its extrapolation.

### 3. The m/z we use may not match

Our analyte CSV uses $m/z = 428.21$ while the compendium reports $m/z = 428.03$
for $[\text{M+H}]^+$. The difference (0.18 Da) suggests our $m/z$ may be the
monoisotopic mass + 0.18 or a different adduct. This would affect which peak we
extract from the IM data, but not the CCS calculation itself.

### 4. Structural limitations of protein calibration for small molecules

The Ruotolo & Robinson protocol was designed for **protein** CCS calibration. Proteins
and small molecules have fundamentally different ion-gas interaction dynamics:
- Proteins are large, roughly spherical, with many charges
- ADP is small, planar, singly charged

The power-law exponent $X$ may be different for small molecules than for proteins in
the same TWIMS cell. This is a known limitation — see Bush et al. (2010): "calibrants
should be structurally similar to the analytes."

---

## What Would Be Needed for Accurate Small-Molecule CCS

1. **Small-molecule calibrants** — use polyalanine, tetraalkylammonium salts, or other
   small standards with known DTIMS CCS values in the 100–300 Å² range, rather than
   denatured proteins.

2. **Class-specific calibration** — fit the power law using calibrants from the same
   chemical class (nucleotides, peptides, lipids) as the analytes.

3. **Direct DTIMS measurement** — if available, measure ADP CCS directly on a
   drift-tube instrument, bypassing TWIMS calibration entirely.

---

## References

- Picache JA, Rose BS, Balber A, et al. Collision cross section compendium to annotate
  and predict multi-omic compound identities. *Chem. Sci.* 10:983-993, 2019.
  [Paper](https://pubs.rsc.org/en/content/articlelanding/2019/sc/c8sc04396e)
  [Data](https://zenodo.org/records/6860818)

- May JC, Goodwin CR, Lareau NM, et al. Conformational Ordering of Biomolecules in the
  Gas Phase: Nitrogen Collision Cross Sections Measured on a Prototype High Resolution
  Drift Tube Ion Mobility-Mass Spectrometer. *Anal. Chem.* 86(4):2107-2116, 2014.
  [Paper](https://pubs.acs.org/doi/10.1021/ac4038448)

- Zheng X, Aly NA, Zhou Y, et al. A structural examination and collision cross section
  database for over 500 metabolites and xenobiotics using drift tube ion mobility
  spectrometry. *Chem. Sci.* 8:7724-7736, 2017.
  [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC5853271/)

- Bush MF, et al. Collision Cross Sections of Proteins and Their Complexes: A
  Calibration Framework and Database for Gas-Phase Structural Biology. *Anal. Chem.*
  82:9557-9565, 2010.

- Stow SM, et al. An Interlaboratory Evaluation of Drift Tube Ion Mobility-Mass
  Spectrometry Collision Cross Section Measurements. *Anal. Chem.* 89:9048-9055, 2017.
  [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC5744684/)

- Dodds JN, May JC, McLean JA. Ion Mobility Derived Collision Cross Sections to Support
  Metabolomics Applications. *Anal. Chem.* 86(8):3842-3848, 2014.
  [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC4004193/)
