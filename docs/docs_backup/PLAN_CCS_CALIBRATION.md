# Plan: CCS Calibration Pipeline for Waters SYNAPT G2

## 1. What the Paper Says (Ruotolo & Robinson, Nature Protocols 2008)

### The Calibration Procedure (Steps 11–16)

**Step 11** — Correct observed drift time for mass-dependent TOF flight time:
```
t'_D = t_D − (C × √(m/z)) / 1000
```
- `t_D` = observed drift time (ms), from the IM data
- `C` = EDC delay coefficient (instrument-specific, from `_extern.inf`, typically 1.4–1.6)
- `m/z` = mass-to-charge ratio of the ion
- `t'_D` = corrected drift time (ms)

**Step 12** — Calculate charge- and mass-corrected CCS (Ω'):
```
Ω' = Ω / [z × √(1/μ)]
   = Ω / [z × √(1/MW + 1/m_gas)]
```
- `Ω` = literature CCS value (Å², from Table 2)
- `z` = charge state
- `μ` = reduced mass = (MW × m_gas) / (MW + m_gas)
- `m_gas` = neutral gas mass (28.014 Da for N₂)

**Step 13** — Plot ln(t'_D) vs ln(Ω') for all calibrant ions

**Step 14** — Fit the linear relationship:
```
ln(Ω') = X × ln(t'_D) + ln(A)
```
- `X` = slope (the "exponential factor")
- `A` = exp(intercept)
- R² should be > 0.98

**Step 15** — Calculate double-corrected drift time (t''_D):
```
t''_D = (t'_D)^X × z × √(1/MW + 1/m_gas)
```

**Step 16** — Fit Ω vs t''_D (linear):
```
Ω = slope₂ × t''_D + intercept₂
```
R² should be > 0.98. This is the final calibration curve.

### Applying to Unknowns
For an unknown ion at charge state z, observed m/z, observed drift time t_D:
1. `t'_D = t_D − C × √(m/z) / 1000`
2. `MW = m/z × z` (approximate, ignoring charge carrier mass)
3. `t''_D = (t'_D)^X × z × √(1/MW + 1/m_gas)`
4. `CCS = slope₂ × t''_D + intercept₂`

## 2. What Your Script Does (g2_ccs_calibration_a.py)

### Correct Steps
- t'_D calculation: ✓ matches paper
- Ω' calculation: ✓ (second version on line 71, which overwrites the wrong first version)
- First regression (ln vs ln): ✓ gets X and ln(A)
- t''_D calculation: ✓ (second version on line 100, overwrites wrong first version)
- Second regression (Ω vs t''_D): ✓ gets final calibration line
- User data application: ✓ applies calibration correctly

### Issues Found

| # | Issue | Severity |
|---|-------|----------|
| 1 | **Duplicate wrong formulas** — Lines 63–69 compute Ω' with √(μ) instead of √(1/μ), then lines 71–77 overwrite with the correct formula. Same for t''_D (lines 96–99 vs 100–103). The wrong versions use `1/(1/MW + 1/m_gas)` inside the sqrt — this is μ, not 1/μ. | Low (overwritten, but confusing) |
| 2 | **Hardcoded instrument parameters** — EDC coeff, wave velocity, pressure all hardcoded. Should read from `_extern.inf`. | Medium |
| 3 | **Manual Excel input** — User must manually extract drift times from MassLynx and enter them in Excel. Should automate from `_im.txt` data. | High |
| 4 | **No instrument setting validation** — Doesn't check that calibrant and analyte were acquired with matching IMS settings. | High |
| 5 | **mpmath/Decimal overkill** — Uses arbitrary-precision math for values that need ~6 significant figures. float64 has 15. | Low |
| 6 | **Hardcoded file paths** — Not portable. | Medium |
| 7 | **neutral_gas_mass = 28** — Should be 28.014 (molecular mass of N₂). Minor effect. | Low |
| 8 | **No error estimation** — Paper Step 20 describes error propagation. Script doesn't implement it. | Medium |

### Is the Paper Wrong?
The paper's procedure is sound and well-established. One subtlety: the EDC correction (Step 11) is described as an approximation for the m/z-dependent flight time from the Transfer T-wave to the pusher. The paper notes C "varies slightly from instrument to instrument" and values from the EDC setup are "only valid if the exit lens of the transfer T-wave guide and transfer ion optics are the same as used when setting up the EDC calibration." This is a practical concern but not a theoretical error.

The paper's two-step fitting (ln-ln first, then linear) is a valid approach to handle the nonlinear relationship between drift time and CCS in a travelling-wave instrument. The exponent X accounts for the fact that TWIMS doesn't have a simple linear drift-time-to-CCS relationship like drift-tube IMS does.

## 3. Instrument Setting Comparison

From comparing `_extern.inf` across runs:

| Parameter | UBQ/MYO/CYTC (analytes) | CALI_UBQ (calibrant) | Match? |
|---|---|---|---|
| IMS Wave Velocity | 400 m/s | 400 m/s | ✓ |
| IMS Wave Height | 25.0 V | 25.0 V | ✓ |
| IMS Gas Flow | 80.00 mL/min | 80.00 mL/min | ✓ |
| EDC Delay Coefficient | 1.5800 | 1.5800 | ✓ |
| **Trap Gas Flow** | **0.40 mL/min** | **0.20 mL/min** | **MISMATCH** |
| Trap Collision Energy | 4.0 | 4.0 | ✓ |
| Transfer Collision Energy | 0.0 | 0.0 | ✓ |
| HeliumCellGasFlow | 200.00 | 200.00 | ✓ |
| All DC voltages | identical | identical | ✓ |
| All RF offsets | identical | identical | ✓ |
| Mass range | 250–3000 | 50–2000 | differs (OK) |

**Critical concern**: Trap Gas Flow differs 2× (0.20 vs 0.40 mL/min). The paper (Step 10) explicitly states: "Use precisely the same instrument conditions (including pressures) for all elements downstream of the trapping ion guide." The trap is upstream of the IMS cell, but the trap gas flow affects the trap pressure (2.39e-2 vs 2.50e-2 mbar) which could affect ion activation before IM separation. This is a **warning**, not a hard stop — it may or may not significantly affect results depending on the ions.

## 4. Calibrant Data (Table 2 from Paper)

Three proteins used as calibrants, with literature CCS values:
- **Cytochrome c** (12,359 Da): charge states 7–18, CCS 1,247–2,766 Å²
- **Equine Myoglobin** (16,952 Da): charge states 8–22, CCS 1,673–3,815 Å²
- **Bovine Ubiquitin** (8,565 Da): charge states 6–11, CCS 1,041–1,802 Å²

Our data has: UBQ_DEN_NEW_V1, MYO_DEN_NEW_V1, CYTC_DEN_TUNE_V1 — all three calibrant proteins.

Note: Some charge states in Table 2 are marked with superscript 'a' (metastable distributions) — these should be excluded from calibration (the user's CSV already has `use` column for this).

## 5. Proposed Pipeline Architecture

### New Step: `imms_ccs_calibrate`

```
Inputs:
  - Calibrant .raw files (UBQ, MYO, CYTC) or their _im.txt/_ms.txt
  - Analyte .raw files or their _im.txt/_ms.txt
  - Literature CCS values (built-in Table 2, or user-provided CSV)

Outputs:
  - calibration_curve.json (X, A, slope, intercept, R², calibrant points)
  - calibration_plots.png (ln-ln plot + Ω vs t''_D plot)
  - Per-analyte CCS profiles (drift time → CCS converted)
  - CCS-calibrated 2D IM-MS HTML (optional: add CCS axis to existing viewer)
```

### Processing Steps

```
Step 0: Validate instrument settings
  ├── Read _extern.inf from all .raw files
  ├── Compare IMS-critical parameters
  ├── WARN if Trap Gas Flow differs
  └── HARD STOP if IMS Wave Velocity/Height/Gas Flow differ

Step 1: Peak picking from calibrant IM-MS data
  ├── Load _im.txt for each calibrant protein
  ├── For each known charge state (from Table 2):
  │   ├── Calculate expected m/z = (MW + z×1.00728) / z
  │   ├── Extract drift profile at that m/z (±tolerance)
  │   ├── Find peak drift time (Gaussian fit or centroid)
  │   └── Record: (protein, z, m/z, t_D_obs)
  └── Output: calibrant_peaks.csv

Step 2: Build calibration curve
  ├── For each calibrant ion:
  │   ├── t'_D = t_D - C × √(m/z) / 1000
  │   ├── Ω' = Ω_lit / [z × √(1/MW + 1/28.014)]
  │   ├── ln(t'_D), ln(Ω')
  │   └── t''_D = (t'_D)^X × z × √(1/MW + 1/28.014)
  ├── Fit 1: ln(Ω') = X × ln(t'_D) + ln(A)  → get X, A
  ├── Fit 2: Ω = slope × t''_D + intercept   → get calibration line
  ├── Validate: R² > 0.98 for both fits
  └── Output: calibration_curve.json, calibration_plots.png

Step 3: Apply to analyte
  ├── Load analyte IM-MS data
  ├── For each species of interest (user-specified m/z, z):
  │   ├── Extract drift profile
  │   ├── Convert each drift bin to CCS:
  │   │   t'_D = t_D - C × √(m/z) / 1000
  │   │   t''_D = (t'_D)^X × z × √(1/MW + 1/28.014)
  │   │   CCS = slope × t''_D + intercept
  │   └── Plot CCS profile (abundance vs CCS)
  └── Output: ccs_profiles.csv, ccs_profiles.png
```

### Integration with Current Pipeline

```
deconvovo/
  imms_ccs_calibrate.py    ← NEW step module

scripts/
  waters_unidec_pipeline.py  (existing — plot_im_data could gain CCS axis option)
  imms_ccs_config.json       ← NEW config with calibrant data + settings
  CCS_CALIBRANT_TABLE.csv    ← built-in Table 2 data

DAG integration in cli.py:
  convert → deconv ──→ summary
         → plot
         → ccs_calibrate → ccs_plot   ← NEW branch
```

The CCS step runs in parallel with plot/deconv since it only needs the converted text files.

### Config File (`imms_ccs_config.json`)

```json
{
    "neutral_gas": {"name": "N2", "mass": 28.014},
    "calibrant_proteins": ["ubiquitin", "myoglobin", "cytochrome_c"],
    "calibrant_table": "CCS_CALIBRANT_TABLE.csv",
    "peak_picking": {
        "mz_tolerance_da": 2.0,
        "min_drift_intensity": 0.05,
        "peak_method": "centroid"
    },
    "validation": {
        "min_r_squared": 0.98,
        "ims_settings_must_match": [
            "IMS Wave Velocity",
            "IMS Wave Height",
            "IMS Gas Flow",
            "HeliumCellGasFlow",
            "EDC Delay Coefficient"
        ],
        "ims_settings_warn_if_differ": [
            "Trap Gas Flow",
            "Trap Collision Energy"
        ]
    }
}
```

## 6. Key Design Decisions / Questions

1. **Automated peak picking** — How to identify charge state peaks from the 2D IM-MS data without manual intervention? Approach: for each calibrant protein, we know MW and expected charge states (from Table 2). We calculate expected m/z values, extract drift profiles at those m/z values, and find the dominant drift time peak. This is approximate but should work for well-separated calibrant ions.

2. **Which calibrant run to use?** — We have `CALI_UBQ_0` (UBQ only, different Trap Gas Flow) and `UBQ_DEN_NEW_V1`, `MYO_DEN_NEW_V1`, `CYTC_DEN_TUNE_V1` (all three proteins, matching settings). The protein runs should be used as calibrants since they have matching IMS settings. The CALI run has a Trap Gas Flow mismatch.

3. **Multiple wave heights** — The paper (Steps 9, 17) says to repeat calibration for each wave height used. Our data all uses 25.0 V fixed wave height, so one calibration curve suffices.

4. **Interface for the existing 2D viewer** — The CCS calibration could add a second y-axis (CCS in Å²) to the existing interactive HTML. Or generate separate CCS profile plots. The drift time → CCS conversion is monotonic (for fixed m/z and z), so it's a relabeling of the y-axis for a specific charge state.

5. **Error estimation** — Paper Step 20: total error = √(E_R² + E_Cal² + E_S²) where E_R = reproducibility (replicate stdev), E_Cal = calibration curve error, E_S ≈ 1% (N₂/He conversion). We should implement this.

6. **The proton mass correction** — The script uses MW = m/z × z, but technically MW = m/z × z − z × 1.00728 (subtracting proton masses). For large proteins this is negligible (<0.1%) but should be correct.

## 7. Implementation Order

1. **Phase 1**: Instrument validation + calibration curve from manually provided drift times
   - Read `_extern.inf`, validate settings
   - Build calibration from CSV (like current script, but automated)
   - Verify against existing reference output

2. **Phase 2**: Automated peak picking from calibrant IM-MS data
   - Load _im.txt, extract drift profiles at expected m/z
   - Find peak drift times automatically
   - No manual Excel step needed

3. **Phase 3**: Full integration
   - Add to DAG task queue
   - CCS-calibrated 2D viewer
   - Error estimation
   - Multi-wave-height support (if needed)

## 8. Questions for User

1. The CALI_UBQ_0 run has different Trap Gas Flow (0.20 vs 0.40 mL/min). Should we use the protein runs (UBQ/MYO/CYTC_DEN) as calibrants instead? They have matching IMS settings and contain the same calibrant proteins.

2. For the analyte species of interest — do you have specific m/z and charge states to target, or should the tool auto-detect major features from the 2D IM-MS data?

3. The current calibrant CSV has drift times from a different instrument/settings (IMS wave velocity 650 m/s, wave height 40 V — from the script's hardcoded values). Our data uses 400 m/s / 25 V. We need NEW drift times from OUR calibrant runs. The literature CCS values (Ω column) are universal, but the observed drift times must come from our instrument.

4. Should CCS calibration be a separate step (like `python -m deconvovo.imms_ccs_calibrate`) or integrated into the main pipeline?
