# HDX-MS Python Package Source Code Analysis

Research date: 2026-03-25

Analysis of open-source Python packages for HDX-MS (Hydrogen-Deuterium Exchange Mass
Spectrometry) data processing. Key question: does any established Python package go from
raw mzML spectra to deuterium uptake, or do all tools require pre-processed input from
vendor software (DynamX, HDExaminer)?

---

## Executive Summary

**No established Python package provides a complete, production-quality pipeline from raw
mzML to deuterium uptake.** The HDX-MS field is split:

1. **Vendor software** (Waters DynamX, Thermo HDExaminer/BioPharma Finder) handles raw
   data: peptide identification, chromatographic extraction, isotope envelope detection,
   and centroid mass calculation. These export CSV peptide tables.

2. **Open-source Python tools** (PyHDX, HaDeX, Deuteros, DECA, HDXBoxeR) consume those
   pre-processed CSV tables and perform downstream analysis: back-exchange correction,
   protection factor calculation, kinetic fitting, statistical comparison, visualization.

3. **Two Python tools attempt mzML reading** but with significant limitations:
   - **TheDeuteriumCalculator** — reads mzML, computes centroids, but requires a
     pre-identified peptide list (from database search) and is a single 1400-line script
   - **pyHXExpress** — has mzML reading utilities (via pyteomics) but primarily expects
     tabular spectral data exported from HDExaminer; mzML support is secondary

4. **ExMS2** reads mzML natively and is the closest to a complete pipeline, but it is
   written in **MATLAB** (compiled standalone, no MATLAB license required).

---

## 1. PyHDX — Derive deltaG from HDX-MS Data

**Repository:** https://github.com/Jhsmit/PyHDX
**Author:** Jochem Smit
**Language:** Python (pip install pyhdx)
**Citation:** Smit et al., Analytical Chemistry

### What It Does

PyHDX derives Gibbs free energy (deltaG) of exchange at single-residue resolution from
peptide-level HDX-MS data. It does NOT read raw spectra — it requires pre-processed
peptide uptake tables.

### Input Format

Accepts CSV files exported from:
- **Waters DynamX** (format `"DynamX_v3_state"`)
- **HDExaminer** (format `"HDExaminer_peptide_pool"`)

Standardized columns (from `hdxms-datasets` package):
```
start, end, sequence, state, exposure, centroid_mz, rt, rt_sd, uptake,
uptake_sd, maxuptake, fd_uptake, fd_uptake_sd
```

The `fileIO.py` module reads CSV with JSON metadata headers. The `batch_processing.py`
module's `StateParser` class handles multi-state experiments.

### Key Source Files

```
pyhdx/
├── process.py          # Core uptake correction, back-exchange, RFU calculation
├── models.py           # HDXTimepoint, HDXMeasurement, Coverage, PeptideUptakeModel
├── fileIO.py           # Read DynamX/HDExaminer CSV, write results
├── fitting.py          # scipy/symfit kinetic fitting
├── fitting_torch.py    # PyTorch-based deltaG fitting
├── fit_models.py       # Exponential association/dissociation models
├── batch_processing.py # Multi-state batch analysis
├── datasets.py         # Integration with hdxms-datasets package
├── support.py          # DataFrame intersection, time conversion utilities
├── alignment.py        # Sequence alignment for multi-protein comparison
├── config.py           # Configuration management
└── web/                # Panel-based web application
```

### How It Reads Peptide Data

From `pyhdx/fileIO.py`:
- `read_dynamx()` — legacy reader for DynamX CSV, handles time unit conversion
- `DataFile` dataclass — wraps format-specific readers from `hdxms_datasets.formats.FMT_REGISTRY`
- `csv_to_hdxm()` — converts CSV DataFrames into `HDXMeasurement` or `HDXMeasurementSet`

Column type casting:
```python
PEPTIDE_DTYPES = {"start": int, "end": int, "stop": int, "_start": int, "_stop": int}
```
The `"stop"` column = `end + 1` (half-open interval convention).

### How It Computes Deuterium Uptake / Back-Exchange Correction

From `pyhdx/process.py`:

**`apply_control()` — Relative Fractional Uptake (RFU):**
```python
def apply_control(experiment, fd_control, nd_control=None, deepcopy=False):
    # RFU = (uptake - nd_uptake) / (fd_uptake - nd_uptake)
    out_df["rfu"] = (u - n) / (f - n)

    # Error propagation:
    out_df["rfu_sd"] = np.sqrt(
        (1 / (f - n)) ** 2 * u_sd**2
        + ((u - f) / (f - n) ** 2) ** 2 * n_sd**2
        + ((n - u) / (f - n) ** 2) ** 2 * f_sd**2
    )
```

If no non-deuterated control is provided, nd_uptake defaults to 0.

**`correct_d_uptake()` — Exchangeable residue correction:**
```python
def correct_d_uptake(peptides, drop_first=1, d_percentage=100, deepcopy=False):
    # 1. Replace P with p (prolines don't exchange)
    # 2. Mark first `drop_first` residues as non-exchanging ("x")
    # 3. Count exchangeable residues = len - prolines - N-term drop
    # 4. Scale by d_percentage (D2O fraction in labeling buffer)
    ex_residues = (len(s) - s.count("x") - s.count("p")) * d_percentage / 100.0
    peptides["ex_residues"] = ex_residues

    # Corrected uptake = RFU * exchangeable residues
    if "rfu" in peptides:
        peptides["uptake_corrected"] = peptides["rfu"] * ex_residues
```

### How It Fits Exchange Kinetics

**Step 1: Intrinsic rate calculation** (via `hdxrate` package):
```python
from hdxrate import k_int_from_sequence
# Returns per-residue intrinsic exchange rates (s^-1)
```

**Step 2: deltaG fitting** (from `fitting_torch.py`):
```python
class DeltaGFit(nn.Module):
    def forward(self, temperature, X, k_int, timepoints):
        # Protection factor = exp(deltaG / RT)
        pfact = torch.exp(self.dG / (constants.R * temperature))
        # Predicted uptake = 1 - exp(-k_int / (1 + pfact) * t)
        uptake = 1 - torch.exp(-torch.matmul((k_int / (1 + pfact)), timepoints))
        # Map residue uptake to peptide uptake via coverage matrix X
        return torch.matmul(X, uptake)
```

The `X` matrix (Np x Nr, peptides x residues) maps residue-level predictions to
peptide-level observables. The `Z` matrix weights by 1/exchangeable_residues.

**Step 3: Kinetic models** (from `fit_models.py`):
- `TwoComponentAssociationModel`: `y = 1 - (r*exp(-k1*t) + (1-r)*exp(-k2*t))`
- `OneComponentAssociationModel`: `y = 1 - exp(-k1*t)`
- Two corresponding dissociation models

**Step 4: Error estimation** (from `fitting_torch.py`):
- Hessian-based: inverse Hessian of MSE loss
- Jacobian-based: SVD for robust covariance, chi-squared degrees of freedom correction

### Key Data Model Classes (from `models.py`)

- **`Coverage`** — peptide-residue correspondence matrix (X, Z matrices)
- **`HDXTimepoint`** — single state + exposure: computes `rfu_peptides`, `d_exp`,
  `rfu_residues`, `rfu_residues_sd` via weighted average
- **`HDXMeasurement`** — full kinetic series for one protein state; computes k_int from
  temperature/pH, provides `get_tensors()` for fitting
- **`HDXMeasurementSet`** — multiple states with padded arrays for simultaneous fitting
- **`PeptideUptakeModel`** — per-peptide kinetic model using analytical or numerical
  integration (scipy `solve_ivp`)

---

## 2. HDXrate — Intrinsic Exchange Rate Calculator

**Repository:** https://github.com/Jhsmit/HDXrate
**Author:** Jochem Smit
**Language:** Python (pip install hdxrate)

### What It Does

Calculates intrinsic amide hydrogen exchange rates from amino acid sequence, temperature,
and pH. Implements the Englander group methodology (Bai et al. 1993, Connelly et al. 1993).

### Source Structure

```
hdxrate/
├── __init__.py
├── hdxrate.py      # All calculation logic
└── constants.txt   # Side-chain inductive effect parameters
```

### Core Function: `k_int_from_sequence()`

```python
def k_int_from_sequence(sequence, temperature, pH_read, reference="poly",
                        exchange_type="HD", d_percentage=100.0,
                        ph_correction=True, wildcard="X", return_sum=True):
    """
    Returns numpy array of per-residue exchange rates (s^-1).
    Sequence: single-letter amino acid codes
    Temperature: Kelvin
    pH_read: measured pH (auto-corrected to pD for HD exchange)
    """
```

### Key Constants and Equations

```python
R = 1.987  # cal/(mol*K)

# Activation energies (cal/mol)
E_act = {"acid": 14000.0, "base": 17000.0, "water": 19000.0}

# Reference rates (poly-DL-alanine): k_acid, k_base, k_water
rates_cat = {
    ("HD", "poly"): tuple((10 ** np.array([1.62, 10.18, -1.5])) / 60),
    ...
}
```

**pH to pD correction:**
```python
pD = pH_read + 0.4 * d_percentage / 100  # Rubinson 2017
```

**Temperature-corrected rates (Arrhenius):**
```python
k_acid = k_acid_ref * np.exp(-E_act["acid"] * (1/T - 1/293) / R)
k_base = k_base_ref * np.exp(-E_act["base"] * (1/T - 1/293) / R)
k_water = k_water_ref * np.exp(-E_act["water"] * (1/T - 1/293) / R)
```

**Per-residue rate calculation:**
```python
# Side-chain inductive factors from neighbor residues
Fa = 10 ** (prev_rho_acid + curr_lambda_acid)  # acid pathway
Fb = 10 ** (curr_lambda_base + prev_rho_base)  # base pathway

# Total intrinsic rate = sum of three pathways
k_int = Fa * k_acid * [D+] + Fb * k_base * [OD-] + Fb * k_water
```

Prolines have zero exchange rate. The first residue (N-terminal amide) is set to infinity
(instantaneous exchange). Side-chain effects of D, E, H are pH/temperature-dependent.

---

## 3. TheDeuteriumCalculator — mzML to Deuterium Uptake

**Repository:** https://github.com/OUWuLab/TheDeuteriumCalculator
**Author:** Thomas Welborn, University of Oklahoma
**Language:** Python (single script, ~1400 lines)
**Citation:** Welborn et al., J. Proteome Res. 2023

### What It Does

**This is the only established Python tool that reads mzML and computes deuterium uptake
directly.** However, it still requires a pre-identified peptide list (from database search).

### Input Requirements

1. **Peptide identification CSV** — columns: ID, ScanNum, Precursor m/z, Charge, Peptide
   (from any search engine output)
2. **mzML files** — converted from vendor raw via ProteoWizard msConvert
3. **Parameters** — PPM tolerance (default 10), sliding window, retention time range
   (default ±30 sec), deuterium recovery rate, D2O fraction

### Key Algorithm: Centroid Mass from Isotope Envelope

**Weighted average mass (Equation 3 from paper):**
```
WeightedAve.mass = [SUM(Em/z_n * I_n / I_T) - 1.00627 Da] * z
```
Where:
- `Em/z_n` = experimental m/z for each isotopic peak
- `I_n` = intensity of individual peaks
- `I_T` = total intensity across envelope
- `z` = charge state

**Theoretical deuterated masses:**
```
Tmass_n = (m/z * z) + n * 1.00628 Da
```
Where n = number of deuterons incorporated (0 to max_D).

**Maximum deuterium:**
```python
def sequence_to_max_deuterium():
    # max_D = (sequence_length - proline_count) * recovery_rate * D2O_fraction
```

### Data Processing Pipeline

1. Read peptide ID list from CSV
2. For each mzML file, iterate MS1 scans via `read_mzml()`
3. Use sliding window to sum intensities across retention time range
4. Binary search (`compare()`) to match theoretical m/z values within PPM tolerance
5. Calculate intensity-weighted centroid m/z (`set_weighted_mass()`)
6. Convert to mass: `mass = centroid_mz * charge - hydrogen_mass`
7. Deuterium uptake = mass_deuterated - mass_undeuterated
8. Quality check: Gaussian fit to envelope, report R^2

### Back-Exchange Handling

The paper states back-exchange is minimized experimentally (subzero-temperature
chromatography) rather than corrected mathematically. The tool includes a
"Deuterium Recovery Rate" parameter that scales the maximum theoretical deuterium,
but there is no explicit fully-deuterated control correction.

### Dependencies
```
lxml, pyteomics, numpy, scipy, pandas, matplotlib
```

### Limitations
- Single monolithic script (not a pip-installable package)
- Requires pre-identified peptide list (does not do peptide identification)
- No built-in back-exchange correction from FD controls
- No multi-state differential analysis
- No kinetic fitting

---

## 4. pyHXExpress — Multimodal HDX-MS Spectral Fitting

**Repository:** https://github.com/tuttlelm/pyHXExpress
**Author:** Leland Tuttle
**Language:** Python
**Citation:** Tuttle et al., JASMS 2025

### What It Does

Python adaptation of HX-Express for automated multimodal (bimodal, trimodal) fitting of
HDX-MS isotope envelopes. Primarily designed for HDExaminer spectral exports but has
mzML reading capability.

### Source Structure

```
pyhxexpress/
├── __init__.py
├── config.py              # Fitting parameters
├── hxex.py                # Core: binomial fitting, peak picking, centroid
├── hxex_fragments.py      # Fragment ion handling
├── hxex_mzml_utils.py     # mzML reading via pyteomics
├── hxex_updating.py       # Updated analysis functions
└── hxex_utils.py          # Plotting, autoscale utilities
```

### mzML Reading (from `hxex_mzml_utils.py`)

```python
from pyteomics import mzml

def read_mzml(mzml_files, scan_idxs=[], ...):
    """Read spectra from mzML via pyteomics.mzml.MzML parser.
    Returns DataFrame with: file, id, idx, scan_start_time_sec, mz, Intensity, MS_Level
    """
    f = mzml.MzML(mzml_file)
    scan = f.get_by_id(id)
    rt = scan['scanList']['scan'][0]['scan start time']
    # Extracts m/z array and intensity array per scan
```

```python
def raw_to_mzml(raw_files, out_dir, ...):
    """Convert Thermo .raw to mzML via ThermoRawFileParser"""
    subprocess.run(['mono', trp_path, '-i=' + raw_file, '-b=' + out_file, ...])
```

```python
def spectra_peakpicker(metadf, spectradf, ...):
    """Extract isotope peaks from raw spectra at predicted m/z positions.
    Uses delta_m = [1.003355, 1.003355, 1.004816, ...] for isotope spacing.
    Calls hxex.peak_picker() for each peptide/charge combination."""
```

### Core Centroid Calculation (from `hxex.py`)

```python
centroid_j = sum(mz * y) / sum(y)  # intensity-weighted average m/z
```

### Binomial Fitting for Multimodal Detection

```python
def binom(bins, n, p):
    """Binomial distribution using real-valued nCk for non-integer deuteration"""

def binom_isotope(bins, n, p):
    """Convolve binomial with natural abundance isotope envelope"""
    # Accounts for 13C natural isotope pattern underlying deuterium distribution

def n_binomials(bins, *params):
    """Multi-population binomial: sum of N populations with different fD values
    params: [scaler, n_1..n_N, mu_1..mu_N, frac_1..frac_N]"""

def count_amides(peptide, count_sc=0.0):
    """Exchangeable amides = len - prolines - N-term drop + sidechain contribution"""
    proline = peptide[config.Nterm_subtract:].count('P')
    n_amides = max(0, len(peptide) - proline - config.Nterm_subtract + int(ex_sc * count_sc))
```

### Statistical Model Selection

F-test between N vs N-1 population fits:
```python
F = ((prev_rss - rss) / 3) / (rss / (n_bins + 1 - n_params))
p = 1.0 - stats.f.cdf(F, 3, n_bins + 1 - n_params)
# Requires p < 0.05 and population fraction > 3%
```

### Input Formats
1. **Primary:** Tabular m/z vs intensity data (HX-Express format or HDExaminer export)
2. **Secondary:** mzML files via `hxex_mzml_utils.read_mzml()` + metadata CSV specifying
   peptide, charge, modification for each target

---

## 5. DeuteRater — NOT HDX-MS (Metabolic Labeling)

**Repository:** https://github.com/JC-Price/DeuteRater
**Author:** JC Price lab, BYU

**Important clarification:** DeuteRater is for **metabolic deuterium labeling** (D2O
drinking studies) to measure **protein turnover rates**, NOT for HDX-MS conformational
dynamics. It reads mzML and tracks isotope envelope shifts over days/weeks to compute
protein synthesis/degradation kinetics. Different application entirely.

---

## 6. Other Notable Tools (Non-Python)

### ExMS2 (MATLAB)
- **The most complete raw-to-uptake pipeline**, but MATLAB-based
- Reads mzML, identifies peptides, extracts isotope envelopes, computes uptake
- Compiled standalone (no MATLAB license needed)
- Source: https://hx2.med.upenn.edu/

### HaDeX (R)
- R package for HDX-MS analysis and visualization
- Reads DynamX CSV tables (not raw mzML)
- CRAN: https://cran.r-project.org/web/packages/HaDeX/

### HDXer (Python, but for MD simulations)
- Computes predicted HDX from molecular dynamics ensembles
- Does NOT process experimental mass spectra
- Source: https://github.com/Lucy-Forrest-Lab/HDXer

---

## 7. The Critical Gap: mzML to Deuterium Uptake

### Standard HDX-MS Data Processing Pipeline

From the Chemical Reviews survey (Skinner et al. 2024):

1. **Peptide identification** — From undeuterated control, typically via PLGS (Waters) or
   Proteome Discoverer (Thermo). This step is universally done by vendor software.

2. **Peak assignment** — Match identified peptides in deuterated spectra by retention time
   and m/z. Vendor software or ExMS2.

3. **Centroid mass calculation** — Intensity-weighted average m/z of isotope envelope:
   ```
   centroid_mz = SUM(mz_i * I_i) / SUM(I_i)
   mass = centroid_mz * z - z * 1.00728  (charge-deconvolved)
   ```

4. **Deuterium uptake** — Mass shift from undeuterated:
   ```
   D_abs = mass(t) - mass(t=0)
   ```

5. **Back-exchange correction** — Normalize to fully-deuterated control:
   ```
   D_frac = [mass(t) - mass(undeuterated)] / [mass(fully_deuterated) - mass(undeuterated)]
   ```

6. **Residue-level resolution** — PyHDX's deltaG fitting, or overlapping peptide
   subtractive methods.

### What Exists in Python

| Step | Python Tool | Status |
|------|------------|--------|
| Raw conversion | pyteomics, pymzML | Mature |
| Peptide ID | None (PLGS/PD/MaxQuant) | Gap |
| Peak extraction | pyHXExpress (partial), TheDeutCalc | Basic |
| Centroid mass | pyHXExpress, TheDeutCalc | Basic |
| Back-exchange | PyHDX `apply_control()` | Mature |
| Kinetic fitting | PyHDX (torch), pyHXExpress (binomial) | Mature |
| deltaG derivation | PyHDX `DeltaGFit` | Mature |
| Intrinsic rates | HDXrate `k_int_from_sequence()` | Mature |

### The Bottleneck

The fundamental bottleneck is **peptide identification from MS/MS spectra**. No Python
HDX tool does this natively. All tools either:
- Require a pre-existing peptide list (TheDeuteriumCalculator, pyHXExpress)
- Require fully-processed uptake tables (PyHDX, HaDeX, DECA, Deuteros)

The closest to a complete open-source Python pipeline would be:
1. **pyteomics** or **pymzML** to parse mzML
2. **pyHXExpress** `read_mzml()` + `spectra_peakpicker()` to extract isotope envelopes
3. Centroid calculation (trivial: `sum(mz*I)/sum(I)`)
4. **PyHDX** `apply_control()` for back-exchange correction
5. **PyHDX** `DeltaGFit` for residue-level deltaG fitting
6. **HDXrate** for intrinsic exchange rates

But step 1 still requires an external peptide identification input.

---

## 8. Key Equations Reference

### Centroid Mass (intensity-weighted average)
```
m_centroid = [SUM_i(mz_i * I_i) / SUM_i(I_i)] * z - z * m_H
```
Where m_H = 1.00728 Da (proton mass).

### Deuterium Uptake (absolute)
```
D_abs(t) = m_centroid(t) - m_centroid(t=0)
```

### Relative Fractional Uptake (RFU)
```
RFU = (uptake - nd_uptake) / (fd_uptake - nd_uptake)
```
Where nd = non-deuterated control, fd = fully deuterated control.

### RFU Error Propagation (PyHDX implementation)
```
sigma_RFU = sqrt(
    (1/(f-n))^2 * sigma_u^2 +
    ((u-f)/(f-n)^2)^2 * sigma_n^2 +
    ((n-u)/(f-n)^2)^2 * sigma_f^2
)
```

### Exchangeable Residues
```
N_ex = (peptide_length - N_prolines - N_term_drop) * D2O_fraction / 100
```
Where N_term_drop = 1 (first two amides exchange too fast to measure).

### Intrinsic Exchange Rate (Englander model)
```
k_int = Fa * k_acid * [D+] + Fb * k_base * [OD-] + Fb * k_water
```
Where Fa, Fb are side-chain inductive factors from neighboring residues.

### Protection Factor
```
PF = k_int / k_obs = exp(deltaG / RT)
```

### Predicted Peptide Uptake (PyHDX forward model)
```
D_peptide(t) = X @ [1 - exp(-k_int / (1 + PF) * t)]
```
Where X is the peptide-residue coverage matrix (Np x Nr).

---

## 9. Relevance to deconvovo

For a Waters IM-MS pipeline that wants to add HDX-MS capability:

1. **Waters instruments use DynamX** for peptide identification and basic uptake
   calculation. There is no avoiding this step.

2. **Post-DynamX analysis** is where Python tools shine:
   - Read DynamX CSV export with PyHDX's `read_dynamx()` or `hdxms-datasets`
   - Apply back-exchange correction with `apply_control()`
   - Fit residue-level deltaG with PyTorch-based `DeltaGFit`
   - Calculate intrinsic rates with HDXrate

3. **If building from raw mzML** (e.g., for non-standard HDX workflows):
   - Use pyteomics/pymzML for mzML parsing
   - Implement centroid calculation (trivial math)
   - The hard part is peptide identification — would need to integrate with
     an existing search engine (Comet, MSFragger, etc.)
   - pyHXExpress shows how to do peak picking at known peptide m/z positions

4. **For Waters-specific raw data** (.raw → mzML):
   - Waters does not support ThermoRawFileParser
   - Would need msConvert (ProteoWizard) for Waters → mzML conversion
   - Or read Waters .raw directly with proprietary SDK

---

## Sources

- [Jhsmit/PyHDX](https://github.com/Jhsmit/PyHDX) — Main PyHDX repository
- [Jhsmit/HDXrate](https://github.com/Jhsmit/HDXrate) — Intrinsic exchange rate calculator
- [OUWuLab/TheDeuteriumCalculator](https://github.com/OUWuLab/TheDeuteriumCalculator) — mzML-reading HDX tool
- [tuttlelm/pyHXExpress](https://github.com/tuttlelm/pyHXExpress) — Multimodal HDX spectral fitting
- [JC-Price/DeuteRater](https://github.com/JC-Price/DeuteRater) — Metabolic labeling (NOT HDX-MS)
- [Lucy-Forrest-Lab/HDXer](https://github.com/Lucy-Forrest-Lab/HDXer) — MD simulation HDX prediction
- [hadexversum/HDX-MS-resources](https://github.com/hadexversum/HDX-MS-resources) — Curated HDX-MS tool list
- [Jhsmit/hdxms-datasets](https://github.com/Jhsmit/hdxms-datasets) — Standardized HDX-MS dataset format
- [Skinner et al. 2024, Chemical Reviews](https://pubs.acs.org/doi/10.1021/acs.chemrev.4c00438) — Comprehensive review of HDX-MS computational tools
- [Welborn et al. 2023, J. Proteome Res.](https://pubs.acs.org/doi/abs/10.1021/acs.jproteome.2c00558) — TheDeuteriumCalculator paper
- [Tuttle et al. 2025, JASMS](https://pmc.ncbi.nlm.nih.gov/articles/PMC11952558/) — pyHXExpress paper
