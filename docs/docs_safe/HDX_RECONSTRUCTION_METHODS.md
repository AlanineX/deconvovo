# HDX-MS Deuterium Uptake: Established Methods vs Our Implementation

Research date: 2026-03-25

This document reviews how established HDX-MS tools compute deuterium uptake
from MS1 data, then compares our `reconstruct_hdx_from_raw.py` script against
those methods. The goal is to determine whether our reconstruction pipeline is
scientifically valid.

---

## 1. How Deuterium Uptake Is Measured from an Isotopic Envelope

### 1.1 The Centroid Mass Shift Method

The standard method for measuring deuterium incorporation in HDX-MS is the
**centroid mass shift** (also called the intensity-weighted average mass
method). When a peptide incorporates deuterium, each D-for-H substitution adds
approximately 1.00628 Da. This shifts the entire isotopic envelope to higher
m/z. The centroid of that envelope tracks the average mass increase.

**Centroid on the m/z scale** (first moment of the isotopic distribution):

```
         sum_i( (m/z)_i * I_i )
C_mz = ---------------------------
            sum_i( I_i )
```

where the sum runs over **all peaks in the isotopic envelope** of the peptide
at a given charge state z, `(m/z)_i` is the m/z of the i-th isotopic peak,
and `I_i` is its intensity.

**Conversion to neutral mass** from centroid m/z:

```
M = C_mz * z  -  z * m_proton
```

where `m_proton = 1.00728 Da` (charge carrier mass).

This is the whole-envelope centroid -- not specific peaks. You use all resolved
isotopic peaks above a noise threshold, weighted by their intensities.

**References:**
- Zhang & Smith (1993) *Protein Science* 2:522-531 -- the foundational HDX-MS
  method paper that introduced peptide-level H/D exchange with LC-MS
- Engen & Smith (2001) *Anal Chem* 73:256A-265A -- established the centroid
  method as the standard for bottom-up HDX
- Masson et al. (2019) *Nature Methods* 16:595-602,
  doi:10.1038/s41592-019-0459-y -- community consensus guidelines

### 1.2 What the Centroid Represents

The centroid mass shift gives the **average number of deuterons incorporated**
across the population of molecules. Under EX2 kinetics (the common regime),
the isotopic envelope is a single broadened distribution that shifts gradually
to higher mass. The centroid tracks this shift directly:

```
Deuterium_uptake(t) = M(t) - M(0)    [in Daltons]
```

where `M(t)` is the centroid mass at labeling time t, and `M(0)` is the
centroid mass of the undeuterated control.

Under EX1 kinetics (cooperative unfolding), the envelope becomes bimodal, and
the centroid still gives the average but may mask the underlying dynamics.
Envelope deconvolution methods are needed for EX1.

### 1.3 How Established Tools Compute It

**DynamX (Waters):** Uses peptide lists from ProteinLynx Global SERVER (PLGS).
For each peptide/charge/timepoint, DynamX finds the isotopic envelope in the
MS1 scan, computes the centroid (intensity-weighted average m/z of the
envelope), converts to mass. The absolute uptake = centroid mass(t) - centroid
mass(undeuterated). DynamX does NOT perform back-exchange correction itself --
it exports raw uptake values.

**HDExaminer (Sierra Analytics / Trajan Scientific):** Fully automated.
Matches theoretical isotope patterns to observed envelopes, computes the
centroid of the matched envelope. Reports uptake as mass shift vs undeuterated.
Color-codes confidence based on how well observed vs theoretical patterns
match. Also does not perform back-exchange correction internally.

**PyHDX (Jhsmit/PyHDX, MIT license):** Does NOT compute centroids from raw
spectra. It imports pre-computed peptide uptake tables (typically DynamX "state
data" CSV exports). Its contribution is downstream: it corrects for
back-exchange, computes relative fractional uptake (RFU), and derives per-
residue Gibbs free energies of exchange using regularized fitting of
overlapping peptides. Input columns expected: `start`, `stop`, `sequence`,
`state`, `exposure`, `uptake`, `uptake_sd`.

**DECA (komiveslab/DECA, Python):** Also does NOT compute centroids from raw
spectra. It imports peptide uptake tables from DynamX or HDX Workbench. Its
role is post-processing: back-exchange correction, overlapping peptide
segmentation (OPS), statistical significance testing, and visualization.

**The Deuterium Calculator (OUWuLab):** This is the closest to what our script
tries to do -- it reads mzML files and computes deuterium uptake from raw
spectra. It uses:
1. Theoretical masses to locate peptides via Accurate Mass and Time (AMT)
2. Gaussian fitting to evaluate isotope distribution quality (R^2 cutoff)
3. Weighted average mass: `WeightedAveMass = sum(Em/z_n * I_n / I_T) * z - z * 1.00627`
4. Uptake = weighted average mass(deuterated) - weighted average mass(undeuterated)

**HaDeX (R package):** Reads DynamX cluster output. The `Center` column
contains the geometric centroid of the isotopic envelope. It converts to mass
via `exp_mass = Center * z - z * proton_mass`, then aggregates across charge
states using intensity-weighted mean:
`avg_exp_mass = weighted.mean(exp_mass, Inten)`.

---

## 2. Back-Exchange Correction

### 2.1 What It Is and Why It Is Needed

During the HDX-MS workflow, after deuterium labeling is quenched (pH 2.5,
0 degC), the protein is digested and peptides are separated by LC before MS
analysis. During these steps (typically 5-15 minutes at low pH), some
incorporated deuterium atoms exchange back to hydrogen from the aqueous LC
solvents. This is called **back-exchange**.

Typical back-exchange is 10-40% of incorporated deuterium. Without correction,
all uptake measurements are systematically underestimated.

### 2.2 The Correction Formula

The standard back-exchange correction (Zhang & Smith 1993, Equation 2 in
Konermann et al. 2011):

```
                    m(t) - m(0)
D_corrected(t) = ---------------  *  N
                  m(FD) - m(0)
```

where:
- `m(t)` = centroid mass of deuterated peptide at time t
- `m(0)` = centroid mass of the undeuterated control
- `m(FD)` = centroid mass of the fully deuterated control
- `N` = maximum number of exchangeable backbone amide hydrogens

The **relative fractional uptake (RFU)** is the normalized version without
multiplying by N:

```
                m(t) - m(0)
RFU(t)    = ---------------
              m(FD) - m(0)
```

This ranges from 0 (no exchange) to 1 (full exchange at all sites).

### 2.3 Calculating N (MaxD -- Maximum Exchangeable Amides)

N = number of residues in the peptide
    - number of proline residues (no backbone amide NH)
    - 1 (for the N-terminal residue, whose amine back-exchanges too fast)

Special cases:
- If the peptide has an N-terminal proline: N = length - prolines (not -1,
  because proline already lacks the amide)
- Some protocols subtract 2 instead of 1 (excluding the penultimate residue
  if it also back-exchanges rapidly)

Formula: `N = len(peptide) - count(Pro) - 1` (or -2 for fast back-exchangers)

### 2.4 The Fully Deuterated (maxD) Control

The fully deuterated control is a sample where the protein has been incubated
in high-concentration deuterium (typically >90% D2O) under denaturing
conditions for extended time (hours to days) to achieve complete exchange.
It is then quenched and analyzed under the SAME LC-MS conditions as the
experimental samples.

Purpose: Since back-exchange is sequence-dependent and varies by peptide,
the maxD control measures the actual back-exchange for each specific peptide
under the exact experimental conditions used. This peptide-specific correction
is far more accurate than a global estimate.

**Do all tools require it?**
- For comparative/differential HDX (state A vs state B): No. If both states
  experience the same back-exchange, the difference cancels. Masson et al.
  (2019) state that for relative comparisons, the maxD control is not
  strictly required.
- For absolute quantification of deuterium incorporation: Yes. Without maxD,
  you cannot know how many deuterons were truly incorporated.

### 2.5 D2O Percentage Correction

If the labeling buffer is not 100% D2O (typically 85-95% due to dilution),
an additional correction is needed:

```
                  m(t) - m(0)
D_corrected = -----------------  *  N
               m(FD) - m(0)

where m(FD) itself is affected by D2O%, or alternatively:

D_corrected = (m(t) - m(0)) / (D2O_fraction * N)
```

PyHDX implements this via its `d_percentage` parameter in `correct_d_uptake()`.

---

## 3. How Established Tools Handle the Peptide-Level Workflow

### 3.1 DynamX (Waters) Workflow

1. **Peptide identification:** PLGS performs data-independent acquisition
   (MSe) on undeuterated samples to build a peptide list with retention times.
2. **Isotope envelope detection:** For each peptide in each deuterated run,
   DynamX searches for the isotopic envelope at the expected m/z and RT
   (within a tolerance window).
3. **Centroid calculation:** Computes intensity-weighted centroid of the
   isotopic envelope for the best charge state.
4. **Uptake output:** Reports absolute mass shift (Da) per peptide per
   timepoint. No back-exchange correction.
5. **Export:** CSV "state data" files consumed by downstream tools (PyHDX,
   DECA, HaDeX, Deuteros, HDXBoxeR).

### 3.2 HDExaminer (Sierra Analytics / Trajan)

1. Imports peptide lists from Mascot, PLGS, Byonic, or manual entry.
2. Automatically finds isotopic envelopes across all charge states.
3. Scores envelope matches (high/medium/low confidence) based on fit to
   theoretical isotope distributions.
4. Computes centroid mass shift.
5. Allows manual inspection and correction of individual results.
6. Does NOT do back-exchange correction internally.

### 3.3 PyHDX (github.com/Jhsmit/PyHDX)

PyHDX does NOT process raw spectra. It reads pre-computed uptake data:

**Input:** DynamX "state data" CSV with columns: `start`, `stop`, `sequence`,
`state`, `exposure`, `uptake`, `uptake_sd`, `maxuptake`.

**Processing (from `process.py`):**
1. `correct_d_uptake()`: Adjusts exchangeable residue count:
   - Marks N-terminal residues (drop_first=1 or 2) as non-exchanging ('x')
   - Converts prolines to non-exchanging ('p')
   - Computes `ex_residues = len(seq) - count('x') - count('p')`
   - Adjusts for D2O percentage
2. `apply_control()`: Computes RFU with back-exchange correction:
   ```
   RFU = (uptake - nd_uptake) / (fd_uptake - nd_uptake)
   ```
   where fd = fully deuterated control, nd = non-deuterated control.
   Propagates uncertainty via standard error propagation.
3. Downstream: fits RFU curves to extract per-residue DeltaG using
   regularized least squares over overlapping peptides.

### 3.4 DECA (github.com/komiveslab/DECA)

**Input:** DynamX or HDX Workbench CSV exports.

**Processing:**
1. Back-exchange correction (global or time-dependent with LEAP robot
   adjustment)
2. Overlapping Peptide Segmentation (OPS): resolves overlapping peptides
   by subtraction: `new_uptake = peptide1_uptake - peptide2_uptake`
3. Double-exponential curve fitting: `y = a*(1 - e^(-bt)) + c*(1 - e^(-0.01t))`
4. Statistical tests (Welch's t-test, ANOVA, Tukey)
5. Visualization (heatmaps, uptake plots, PyMOL scripts)

### 3.5 Masson et al. (2019) Community Guidelines

Key recommendations from doi:10.1038/s41592-019-0459-y:

- **Undeuterated control:** Required for all experiments (defines m(0))
- **Fully deuterated control:** Recommended but not strictly required for
  differential comparisons
- **Replicates:** Minimum 3 technical replicates per condition
- **Timepoints:** Multiple timepoints spanning seconds to hours (log-spaced)
- **Reporting:** Must report peptide-level data including sequence, charge
  state, retention time, number of exchangeable amides, and uptake values
- **Back-exchange:** Should be reported; if maxD control used, specify
  preparation method
- **Statistical significance:** Use global or per-peptide significance
  thresholds

---

## 4. Correct Workflow for Computing Exchange from Raw mzML

### Step 1: Extract the Isotopic Envelope for a Known Peptide

For each target peptide at each charge state z:

a. Compute theoretical monoisotopic m/z:
   ```
   mono_mz = (M_theoretical + z * m_proton) / z
   ```

b. Generate expected isotopic peak positions:
   ```
   peak_i_mz = mono_mz + i * (neutron_mass / z)     for i = 0, 1, 2, ..., n
   ```
   where neutron_mass = 1.003355 Da (mass difference between 13C and 12C).

c. In each MS1 scan within the expected retention time window, search for
   peaks matching the expected m/z values within tolerance (typically 10-20
   ppm).

d. Validate: require a minimum number of isotopic peaks detected (typically
   3+), and optionally score the match quality vs theoretical isotope
   distribution (using averagine model or exact formula).

**Citations:**
- Senko et al. (1995) *JASMS* 6:229-233 -- averagine model for isotope
  distributions
- The Deuterium Calculator (2023) *J Proteome Res* -- AMT approach for
  mzML-based HDX extraction

### Step 2: Compute the Centroid Mass

For the identified isotopic envelope at charge state z:

```
C_mz = sum_i( mz_i * I_i ) / sum_i( I_i )

M_centroid = C_mz * z  -  z * m_proton
```

This is the **observed centroid mass** for this peptide at this timepoint.

Critical: the centroid should be computed from the OBSERVED m/z values
(not the theoretical ones). The whole point is that deuteration shifts
the observed envelope.

### Step 3: Compute Deuterium Uptake

**Absolute uptake (Da):**
```
Delta_m(t) = M_centroid(t) - M_centroid(undeuterated)
```

**Number of deuterons:**
```
D(t) = Delta_m(t) / 1.00628
```

(approximately, since D-H mass difference = 2.01410 - 1.00782 = 1.00628 Da)

In practice, `Delta_m(t)` in Da IS the number of deuterons to a good
approximation (since each D adds ~1.006 Da).

### Step 4: Normalize (Back-Exchange Correction)

**Fractional uptake (requires maxD control):**
```
                 M_centroid(t) - M_centroid(0)
Fractional(t) = ---------------------------------
                 M_centroid(FD) - M_centroid(0)
```

**Absolute corrected uptake:**
```
D_corrected(t) = Fractional(t) * N
```

where `N = len(peptide) - count(Pro) - 1` (or -2).

**Without maxD control (theoretical normalization, less accurate):**
```
                 M_centroid(t) - M_centroid(0)
Fractional(t) = ---------------------------------
                       N * 1.00628
```

This assumes zero back-exchange, which systematically underestimates the true
fractional uptake.

### Step 5: Handle Multiple Charge States

When the same peptide is observed at multiple charge states, the standard
approach (Houde et al. 2011, Weis et al. 2006) is:

a. Compute centroid mass independently for each charge state.
b. Combine using intensity-weighted average:
   ```
   M_combined = sum_z( I_z * M_z ) / sum_z( I_z )
   ```
   where `I_z` is the total envelope intensity at charge state z.

c. Alternatively, use only the most abundant charge state (simpler, used by
   some tools).

d. Check consistency: if the mass from different charge states disagrees by
   more than ~0.5 Da, flag as a potential misassignment or interference.

**References for charge state handling:**
- Weis (2006) *JASMS* 17:1700-1703
- Houde et al. (2011) *Anal Chem* 83:5789-5794
- Hageman & Bhatt (2024) *Chem Rev* -- 2024 comprehensive review,
  doi:10.1021/acs.chemrev.4c00438

---

## 5. What Our Current Script Does Wrong

Our script: `skills/public/hdx-process-extract/scripts/reconstruct_hdx_from_raw.py`

### Issue 1: Does NOT Compute the Centroid Mass (CRITICAL)

The script computes a `neutral_mass_est` but NOT via the centroid method:

```python
neutral_series = (matched_mz[used] * z - z * PROTON) - iso_idx[used] * NEUTRON
neutral_mass_est = float(np.average(neutral_series, weights=matched_int[used]))
```

This converts each matched isotopic peak BACK to the monoisotopic neutral
mass by subtracting `i * neutron`, then averages those estimates. This is
computing an **estimate of the monoisotopic mass**, not the centroid of the
deuterated envelope. It deliberately removes the isotope shift information.

**What the centroid method does differently:** It computes the
intensity-weighted average of the OBSERVED m/z values (or their
corresponding masses), WITHOUT subtracting the isotope index offsets. The
centroid INCLUDES the isotope shift -- that is the whole measurement.

Correct centroid:
```python
M_centroid = np.average(matched_mz[used] * z - z * PROTON, weights=matched_int[used])
```

Our script:
```python
# WRONG: subtracts isotope offsets, collapsing everything to monoisotopic mass
neutral_series = (matched_mz[used] * z - z * PROTON) - iso_idx[used] * NEUTRON
neutral_mass_est = np.average(neutral_series, weights=matched_int[used])
```

The `- iso_idx[used] * NEUTRON` term undoes the very signal we need to
measure. It strips out the information about how the envelope has shifted.
The resulting `mass_error_da` then reflects only mass accuracy noise, not
deuterium incorporation.

**This is the fundamental flaw: the script measures mass accuracy, not
deuterium exchange.**

### Issue 2: No Centroid Mass Recorded Per Timepoint

Because the script computes monoisotopic mass estimates (Issue 1), it has no
concept of centroid mass at different timepoints. The `delta_mass_vs_baseline`
metric compares monoisotopic mass estimates between runs, which should be
approximately zero for all timepoints if the mass measurement is accurate.
Any non-zero delta is noise from mass accuracy variation, not deuterium signal.

### Issue 3: No Undeuterated Control

The script uses `peptide_mapping` run as baseline, which may or may not be
an undeuterated control. In proper HDX-MS:
- The undeuterated control is protein that was NEVER exposed to D2O
- `peptide_mapping` in Waters workflows may refer to the identification run,
  which is indeed undeuterated
- But the script does not verify this or handle the case where it is absent

### Issue 4: No Fully Deuterated (maxD) Control

There is no maxD control handling. The script normalizes `exchange_ratio_vs_max`
by dividing by the longest timepoint's mass shift, which is NOT the same as
a fully deuterated control. The longest timepoint likely has incomplete exchange.
This makes the "exchange ratio" meaningless as a fractional uptake measure.

### Issue 5: No Back-Exchange Correction

Even if the centroid were computed correctly, there is no back-exchange
correction. All uptake values would be systematically underestimated by
10-40%.

### Issue 6: No Exchangeable Amide Count (N / MaxD)

The script does not compute N (number of exchangeable amide hydrogens) for
any peptide. It does not account for prolines or the N-terminal residue.
Without N, it is impossible to report fractional uptake or corrected absolute
uptake in the way the field expects.

### Issue 7: Picks Only the Best-Scoring Scan (Loses Chromatographic Information)

The script keeps only the single "best" scan for each peptide per run. In
proper HDX analysis:
- The peptide elutes over an LC peak spanning multiple scans
- The centroid should be computed from the summed/averaged envelope across
  the elution window (or from the apex scan with verification)
- Using only the highest-scoring scan ignores chromatographic consistency

### Issue 8: Arbitrary Scoring Weights

```python
weights = np.array([1.0, 0.7, 0.5, 0.35, 0.25, 0.18, 0.13, 0.1])
score = float((matched_int * w).sum() * (k / matched_int.size))
```

These hardcoded weights down-weight higher isotopic peaks, which biases toward
the monoisotopic peak. Established tools use either:
- Theoretical isotope distribution matching (chi-squared or cosine similarity)
- Gaussian fitting to evaluate envelope quality (The Deuterium Calculator)
- No weighting for the centroid calculation itself (all peaks weighted by
  their observed intensity)

### Issue 9: No Retention Time Constraint

The script searches ALL MS1 scans in the run for each peptide. In proper
HDX analysis, peptides are found at specific retention times (from the
identification run), and only scans within an RT window (+/- 30s typical)
are searched. Without RT constraints, the script may match noise features
or different peptides with similar m/z.

### Issue 10: No Support for Multiple Charge States Per Peptide

The script picks the single best charge state for each peptide. It does not
combine information across charge states. While using the best charge state
is a valid simplification used by some tools, the script should at least
verify consistency across charge states.

### Summary of What Needs to Change

| Aspect | Current Script | Correct Method |
|--------|---------------|----------------|
| Mass metric | Monoisotopic mass estimate | Centroid mass of envelope |
| Key calculation | `mass - i * neutron` (strips signal) | Weighted average of observed m/z |
| Undeuterated reference | `peptide_mapping` (assumed) | Explicit undeuterated control |
| maxD control | None (uses max timepoint) | Fully deuterated control |
| Back-exchange correction | None | `(m(t)-m(0))/(m(FD)-m(0)) * N` |
| Exchangeable amides N | Not computed | `len(seq) - prolines - 1` |
| Scan selection | Best single scan | RT-windowed, chromatographic peak |
| Charge state handling | Best only | Weighted average or consistency check |
| Envelope quality | Arbitrary weights | Theoretical pattern matching |

---

## 6. Open-Source Reference Implementations

### 6.1 PyHDX (github.com/Jhsmit/PyHDX)

- **What it does:** Derives per-residue DeltaG from peptide-level HDX data
- **Input:** Pre-computed uptake data (DynamX CSV export)
- **Does NOT read raw spectra** -- starts from peptide uptake tables
- **Key functions:** `correct_d_uptake()` (adjusts for back-exchange, prolines,
  D2O%), `apply_control()` (computes RFU from FD/ND controls)
- **RFU formula:** `RFU = (uptake - nd_uptake) / (fd_uptake - nd_uptake)`
- **License:** MIT
- **Citation:** Smit et al. (2021) *Anal Chem* 93:9497-9504

### 6.2 HDXBoxeR (CRAN, github.com/mkajano/HDXBoxeR)

- **What it does:** Statistical analysis and visualization of multiple HDX-MS
  protein states
- **Input:** HDExaminer exports
- **Computes:** Back-exchange statistics, peptide length distributions, Welch's
  T-test and Critical Interval framework
- **Does NOT compute centroids** from raw data
- **Citation:** Janowska et al. (2024) *Bioinformatics* 40(8):btae479

### 6.3 HaDeX (CRAN, github.com/hadexversum/HaDeX)

- **What it does:** Analysis and visualization of HDX-MS data
- **Input:** DynamX cluster data
- **Key formula:** `exp_mass = Center * z - z * proton_mass`, aggregated across
  charge states via `weighted.mean(exp_mass, Inten)`
- **Computes:** Deuteration levels, kinetic curves, Woods plots, differential
  analysis
- **Citation:** Puchala et al. (2020) *Bioinformatics* 36(16):4516-4518

### 6.4 DECA (github.com/komiveslab/DECA)

- **What it does:** Post-processing of HDX-MS data from DynamX/HDX Workbench
- **Input:** Pre-computed uptake tables (CSV)
- **Computes:** Back-exchange correction (global + LEAP robot), overlapping
  peptide segmentation, statistical significance, curve fitting
- **Does NOT compute centroids** from raw data
- **Citation:** Lumpkin & Komives (2019) *Mol Cell Proteomics* 18:2516-2523

### 6.5 The Deuterium Calculator (github.com/OUWuLab/TheDeuteriumCalculator)

- **What it does:** Computes deuterium uptake directly from mzML files
- **Input:** mzML raw data + peptide library CSV
- **Processing:** Accurate Mass and Time (AMT) matching, Gaussian fitting of
  isotopic distributions, weighted average mass calculation
- **Key formula:** `WeightedAveMass = sum(mz_n * I_n / I_total) * z - z * 1.00627`
- **Output:** IC-HDX standardized tables, DECA-compatible CSVs
- **This is the best reference implementation for what our script tries to do**
- **Citation:** Espada et al. (2023) *J Proteome Res* 22:585-592

### 6.6 Deuteros / Deuteros 2.0 (github.com/andymlau/Deuteros_2.0)

- **What it does:** Rapid analysis and visualization of differential HDX-MS
- **Input:** DynamX exports
- **Computes:** Woods plots, volcano plots, heatmaps with statistical cutoffs
- **Does NOT compute centroids** from raw data
- **Citation:** Lau et al. (2019) *Bioinformatics* 35(17):3171-3173

### 6.7 DeuteRater (github.com/PMSeitzer/DeuteRater)

- **What it does:** Kinetic curve calculation for proteins labeled with D2O
  (metabolic labeling, different from HDX-MS backbone exchange)
- **Not directly applicable** to backbone HDX-MS experiments

### 6.8 HDXer (github.com/Lucy-Forrest-Lab/HDXer)

- **What it does:** Computes theoretical HDX from MD simulations, compares to
  experiment, ensemble refinement
- **Not for experimental data processing** but useful for validation

---

## 7. Key Papers and References

| Paper | DOI | Relevance |
|-------|-----|-----------|
| Zhang & Smith (1993) *Protein Sci* | 10.1002/pro.5560020404 | Foundational HDX-MS method, back-exchange correction formula |
| Masson et al. (2019) *Nature Methods* | 10.1038/s41592-019-0459-y | Community consensus guidelines for HDX-MS |
| Engen (2009) *Anal Chem* | 10.1021/ac802405a | Review of HDX-MS analysis methods |
| Houde et al. (2011) *Anal Chem* | 10.1021/ac200258g | HDX-MS for biopharmaceutical characterization |
| Konermann et al. (2011) *Mass Spec Rev* | 10.1002/mas.20290 | Comprehensive review of HDX-MS theory |
| Weis et al. (2006) *JASMS* | 10.1016/j.jasms.2006.07.025 | Practical considerations for HDX data analysis |
| Hageman & Bhatt (2024) *Chem Rev* | 10.1021/acs.chemrev.4c00438 | 2024 comprehensive review of HDX-MS computational tools |
| Lumpkin & Komives (2019) *Mol Cell Proteomics* | 10.1074/mcp.TIR119.001731 | DECA software paper |
| Smit et al. (2021) *Anal Chem* | 10.1021/acs.analchem.1c02155 | PyHDX software paper |
| Espada et al. (2023) *J Proteome Res* | 10.1021/acs.jproteome.2c00558 | The Deuterium Calculator -- mzML-based HDX processing |
| Puchala et al. (2020) *Bioinformatics* | 10.1093/bioinformatics/btaa587 | HaDeX R package |
| Janowska et al. (2024) *Bioinformatics* | 10.1093/bioinformatics/btae479 | HDXBoxeR R package |

---

## 8. Conclusions and Recommendations

### The script's fundamental flaw is fixable

The core issue is that `best_feature_for_peptide()` subtracts isotope index
offsets when computing `neutral_mass_est`, collapsing the envelope back to
the monoisotopic mass. Removing the `- iso_idx[used] * NEUTRON` term and
instead computing the intensity-weighted average of the observed masses
would give a proper centroid.

### Minimum changes to make the script scientifically valid

1. **Compute centroid mass** (not monoisotopic mass): Remove the isotope
   offset subtraction. Compute `M_centroid = average(observed_mz * z - z * proton, weights=intensity)`.

2. **Require an undeuterated control run** explicitly labeled, and compute
   `uptake(t) = M_centroid(t) - M_centroid(undeuterated)` for each peptide.

3. **Compute N** (exchangeable amides) for each peptide:
   `N = len(seq) - count('P') - 1`.

4. **Report uptake in Da** as the primary metric (this is what the field
   uses for differential comparisons).

5. **Add optional back-exchange correction** if a maxD control is available:
   `fractional(t) = uptake(t) / uptake(FD)`, then `D_corrected = fractional * N`.

6. **Add retention time constraints** so the search only examines scans
   within the expected elution window.

### Recommended reference implementation to study

The Deuterium Calculator (github.com/OUWuLab/TheDeuteriumCalculator) is the
closest open-source tool to what we are trying to build. It reads mzML,
matches peptides by AMT, computes weighted average mass (centroid), and
outputs standard HDX tables. Its source code should be studied before
modifying our script.
