# HDX-MS Pipeline Planning Document

**Research survey compiled: 2026-03-24**
**Purpose: Design a complete open-source HDX-MS data processing pipeline**

---

## Table of Contents

1. [Standard HDX-MS Experimental Workflow](#1-standard-hdx-ms-experimental-workflow)
2. [Open-Source Tools for HDX-MS Data Processing](#2-open-source-tools-for-hdx-ms-data-processing)
3. [Thermo .RAW File Structure and Access](#3-thermo-raw-file-structure-and-access)
4. [Peptide Mapping for HDX-MS](#4-peptide-mapping-for-hdx-ms)
5. [Pipeline Architecture Recommendation](#5-pipeline-architecture-recommendation)

---

## 1. Standard HDX-MS Experimental Workflow

### 1.1 Overview

HDX-MS measures protein conformational dynamics by monitoring the rate at which backbone amide hydrogens exchange with deuterium from D2O solvent. Regions that are solvent-exposed or structurally dynamic exchange faster; buried or hydrogen-bonded regions exchange slower. The technique is widely used for epitope mapping, binding site identification, protein-ligand interactions, and conformational change studies.

**Key reference:** Masson et al., "Recommendations for performing, interpreting and reporting hydrogen deuterium exchange mass spectrometry (HDX-MS) experiments" (Nature Methods, 2019).

### 1.2 Step-by-Step Workflow

#### Step 1: Sample Preparation
- Protein of interest prepared at working concentration (typically 5-50 uM)
- Multiple conditions prepared: apo protein, protein + ligand, protein + binding partner, etc.
- All buffers and conditions must be matched between states being compared

#### Step 2: Deuterium Labeling
- Protein diluted into D2O-based buffer (typically 10-20x dilution, giving ~90% D2O final)
- Labeling carried out at multiple time points: typically 10s, 30s, 1min, 5min, 30min, 4h, 24h
- Undeuterated control (in H2O) and fully deuterated control (maxD) also prepared
- Temperature typically 20-25C during labeling

#### Step 3: Quenching
- Exchange reaction stopped by:
  - Lowering pH to ~2.5 (adding quench buffer with acid)
  - Lowering temperature to 0-4C
- At pH 2.5 and 0C, back-exchange half-life is ~30-120 min, giving a window for analysis
- Quenched samples may be flash-frozen for later analysis or analyzed immediately

#### Step 4: Proteolytic Digestion
- Online digestion using immobilized pepsin column at 0-20C
- Pepsin is active at low pH (unlike trypsin), making it ideal for HDX
- Pepsin is non-specific -- cleavage sites are not perfectly predictable from sequence alone, but are reproducible under consistent conditions
- Digestion typically 2-3 minutes
- Alternative proteases: nepenthesin, fungal protease XIII, protease type XVIII
- Dual protease strategies improve sequence coverage

#### Step 5: LC Separation
- Peptides trapped on C18 trap column at 0C
- Rapid reverse-phase UPLC gradient (typically 5-10 minutes)
- Short gradients minimize back-exchange during chromatography
- All LC components kept at 0C (HDX manager / refrigerated column compartment)

#### Step 6: Mass Spectrometry
- Electrospray ionization (ESI)
- High-resolution MS required to resolve isotopic envelopes
- MS1 data used for deuterium uptake measurement (centroid mass shift)
- MS/MS (MS2) data used for peptide identification (fragmentation spectra)
- Typical instruments:
  - **Waters:** Synapt G2/G2-Si/XS (QTOF with ion mobility), Xevo G2-S QTof
  - **Thermo:** Orbitrap Exploris 480, Orbitrap Eclipse Tribrid, Q Exactive series
- Automation: LEAP PAL robot / CTC PAL-based H/D Platform for reproducible sample handling

#### Step 7: Data Analysis
- Peptide identification from MS/MS data (undeuterated samples)
- Deuterium uptake extraction from MS1 isotopic envelopes (deuterated samples)
- Back-exchange correction using maxD control
- Statistical comparison between states
- Visualization: uptake plots, butterfly plots, Woods plots, heat maps, structural mapping

### 1.3 Software Used at Each Step (Traditional Workflow)

| Step | Waters Platform | Thermo Platform | Open-Source Alternative |
|------|----------------|-----------------|----------------------|
| Peptide ID | ProteinLynx Global Server (PLGS) | Proteome Discoverer + Sequest/Mascot | MSFragger, Comet, Sage |
| Uptake extraction | DynamX 3.0 | BioPharma Finder / HDExaminer PRO | PyHDX, DECA, Deuterium Calculator |
| Statistical analysis | DynamX + Deuteros | HDX Workbench | HDXBoxeR, MEMHDX, HaDeX |
| Visualization | DynamX | HDX Workbench / HDExaminer PRO | HDX-Viewer, Deuteros, HaDeX |
| Structural mapping | PyMOL export | PyMOL export | HDX-Viewer, PyMOL scripting |

### 1.4 Deuterium Uptake Calculation

**Basic centroid mass method:**
```
D_uptake = m(t) - m(0)
```
Where m(t) is the centroid mass of the isotopic envelope at time t, and m(0) is the undeuterated centroid mass.

**Relative fractional uptake (RFU):**
```
RFU = D_uptake / MaxUptake
MaxUptake = N - P - 1   (if N-terminal residue is not proline)
MaxUptake = N - P       (if N-terminal residue is proline)
```
Where N = number of amino acids in peptide, P = number of prolines (which lack amide hydrogens).

**Back-exchange correction:**
```
D_corrected = D_uptake * MaxUptake / (m(maxD) - m(0))
```
Where m(maxD) is the centroid mass of the fully deuterated control. Typical back-exchange: 25-35%.

### 1.5 Exchange Kinetics: EX1 vs EX2

- **EX2 kinetics** (most common): Each amide exchanges independently. Isotopic envelope shifts gradually to higher m/z. Centroid mass analysis is appropriate.
- **EX1 kinetics**: Cooperative unfolding exposes multiple amides simultaneously. Produces characteristic bimodal isotopic distribution. Requires binomial fitting or spectral deconvolution, not simple centroid analysis.
- **EXX regime**: Intermediate behavior; bimodal distribution where both populations shift over time.

---

## 2. Open-Source Tools for HDX-MS Data Processing

### 2.1 Reading Thermo .RAW Files

| Tool | Language | pip install | Linux | Vendor DLL | Active | Notes |
|------|----------|-------------|-------|------------|--------|-------|
| **ThermoRawFileParser** | C# (.NET) | N/A (CLI binary) | Yes (Mono/.NET) | Uses Thermo RawFileReader DLLs (bundled) | Yes (v1.4.4, May 2024) | Gold standard for RAW->mzML conversion. Available via bioconda, Docker. |
| **pyrawr** | Python | `pip install pyrawr` | Yes (requires ThermoRawFileParser) | Indirect (via TRFP) | Moderate | Python wrapper around ThermoRawFileParser CLI. Returns spectra, XICs, metadata as Python dicts. |
| **fisher-py** | Python | `pip install fisher-py` | Yes (requires .NET 8.0 SDK) | Yes (RawFileReader DLLs via pythonnet) | Moderate | Direct Python access to RAW files via pythonnet/.NET interop. Can read on Linux with .NET runtime. |
| **pymsfilereader** | Python | Not on PyPI | **Windows only** | Yes (MSFileReader COM DLLs) | No (legacy) | Requires Windows COM registration. Not viable for Linux. |
| **rawrr** | R | Bioconductor | Yes (with Mono) | Yes (RawFileReader) | Yes | R package; complete OS-independent access to RAW spectral data. Not Python. |

**Recommendation:** Use **ThermoRawFileParser** for RAW-to-mzML conversion (robust, well-maintained, Linux-native via Mono/.NET), then use Python mzML readers for downstream analysis. Alternatively, use **fisher-py** for direct Python access if .NET SDK is available.

### 2.2 Reading mzML / Open Format Files

| Tool | pip install | Version | Linux | Notes |
|------|-------------|---------|-------|-------|
| **pyopenms** | `pip install pyopenms` | 3.3.0 | Yes | Most comprehensive. Full OpenMS bindings. mzML, mzXML, TraML, mzIdentML. Signal processing, feature detection, XIC extraction. Large install (~500MB). |
| **pymzml** | `pip install pymzml` | 2.5.6 | Yes | Lightweight, fast mzML parser. MIT license. Only requires numpy+regex. Good for quick iteration. |
| **pyteomics** | `pip install pyteomics[XML]` | 5.0 | Yes | Versatile framework. Reads mzML, mzXML, MGF, pepXML, mzIdentML, FASTA. Also has in-silico digestion (parser.cleave). Pure Python + lxml. |
| **ms_deisotope** | `pip install ms-deisotope` | 0.0.60 | Yes | Built on pyteomics. Deisotoping and charge state deconvolution. MzMLLoader for indexed random access. C extensions for speed. |
| **matchms** | `pip install matchms` | Latest | Yes | Spectral similarity computation and processing. Good for MS/MS spectrum matching. |
| **spectrum_utils** | `pip install spectrum_utils` | Latest | Yes | MS data processing and visualization. Lightweight. Good for spectrum plotting. |

### 2.3 Peptide Identification (Search Engines)

| Tool | Language | Install | Linux | Speed | Notes |
|------|----------|---------|-------|-------|-------|
| **Sage** | Rust | Native binary; Python: `pip install sagepy` | Yes | Fastest (~5x MSFragger) | MIT-licensed. Cloud-ready. Includes RT prediction, LFQ, FDR control. sagepy provides Python API. |
| **MSFragger** | Java | FragPipe bundle or standalone JAR | Yes (Java) | Very fast | Academic license (free for academic use). Part of FragPipe ecosystem. Best for open search / modified peptides. |
| **Comet** | C++ | Binaries from comet-ms.sourceforge.net | Yes | Fast (fragment-ion indexing in 2024) | Open-source (Apache 2.0). Widely integrated into other tools. |
| **X!Tandem** | C++ | Binaries available | Yes | Moderate | Open-source. Mature but less actively developed. |
| **OMSSA** | C++ | Legacy | Yes | Moderate | NCBI tool. No longer actively maintained. |

**For HDX-MS specifically:** Peptide identification is typically done on undeuterated control samples only. The search database is the protein of interest (not a full proteome), and the protease is set to "non-specific" or "pepsin" (semi-specific cleavage). This makes the search space manageable even for non-specific digestion.

**Recommended approach for our pipeline:**
1. Convert RAW to mzML (ThermoRawFileParser)
2. Run Sage or Comet with non-specific cleavage against the target protein FASTA
3. Parse results with pyteomics (pepXML or mzIdentML reader)
4. Filter by FDR and build peptide list

### 2.4 HDX-MS Specific Analysis Software

#### 2.4.1 Open-Source / Free Tools

| Tool | Language | Install | Platform | Active | What It Does |
|------|----------|---------|----------|--------|-------------|
| **PyHDX** | Python | `pip install pyhdx` | Linux/Mac/Win + Web | Yes | Derives Gibbs free energy at residue level from overlapping peptides. Uses PyTorch for optimization. Web interface available. Reads DynamX/HDExaminer CSV exports. |
| **DECA** | Python | github.com/komiveslab/DECA | Linux/Mac/Win | Moderate | Comprehensive automated post-processing. Reads HDX Workbench or DynamX output. Python 2.7 and 3.7 compatible. |
| **The Deuterium Calculator** | Python | github.com/OUWuLab/TheDeuteriumCalculator | Linux/Mac/Win | Yes | Differential and non-differential analysis. Standardized output per IC-HDX recommendations. Generates Woods plots. |
| **HaDeX** | R | CRAN / github.com/hadexversum/HaDeX | Linux/Mac/Win + Web | Yes (v2) | Interactive analysis and visualization. Web server available (hadex.mslab-ibb.pl). |
| **HDXBoxeR** | R | github.com/mkajano/HDXBoxeR | Linux/Mac/Win | Yes (2024) | Statistical analysis across multiple states/time points. Volcano plots, CI-based significance. |
| **MEMHDX** | R / Web | memhdx.c3bi.pasteur.fr | Web | Yes | Mixed-effects model for statistical validation. Web application. |
| **Deuteros 2.0** | MATLAB | github.com/andymlau/Deuteros_2.0 | Win/Mac | Moderate | Visualization and significance testing. Designed for DynamX output. Requires MATLAB. |
| **HRaDeX** | R | github.com/hadexversum/HRaDeX | Linux/Mac/Win + Web | Yes (2024) | High-resolution deuterium uptake rates. Web server available (hradex.mslab-ibb.pl). |
| **HDX-Viewer** | Web (JS) | github.com/david-bouyssie/hdx-viewer | Web | Moderate | Interactive 3D visualization. Maps HDX data onto PDB structures via NGL viewer. |
| **HD-eXplosion** | Web | hd-explosion.utdallas.edu | Web | Moderate | Web-based analysis platform. |
| **PIGEON-FEATHER** | Python | github.com/glasgowlab/PIGEON-FEATHER | Linux/Mac/Win | Yes (2025) | Residue-level Gibbs free energy (DeltaG_op) from conventional HDX-MS data. State-of-the-art for residue-level resolution. |
| **ReX (RexMS)** | R | github.com/ococrook/RexMS | Linux/Mac/Win | Yes (2025) | Bayesian residue-level HDX inference. Leverages peptide overlap and temporal/sequence correlation. |
| **PyHDX (HDXrate)** | Python | github.com/Jhsmit/HDXrate | Linux/Mac/Win | Yes | Intrinsic exchange rate calculation. Companion to PyHDX. |
| **exPfact** | Python | github.com/skinnersp/exPfact | Linux/Mac/Win | Moderate | Protection factor extraction from HDX-MS data. |
| **HDXer** | Python | github.com/Lucy-Forrest-Lab/HDXer | Linux/Mac/Win | Moderate | Predicts HDX protection factors from MD simulations. |
| **QUDeX-MS** | Java | sourceforge.net/projects/qudex-ms | Linux/Mac/Win | Low | Resolved isotopic fine structure analysis. |

#### 2.4.2 Commercial Tools

| Tool | Vendor | Notes |
|------|--------|-------|
| **DynamX 3.0** | Waters | Standard for Waters HDX systems. Peptide assignment + uptake extraction. Exports CSV state/difference files. |
| **HDExaminer PRO** | Trajan (formerly Sierra Analytics) | Fully automated. Supports Thermo + Waters data. Heatmaps, volcano plots, PyMOL export. |
| **HDX Workbench** | (Academic/Commercial) | Feature-rich. Java-based. Supports Waters + Thermo. Integrates peptide search with uptake analysis. |
| **BioPharma Finder** | Thermo Fisher | Enterprise solution. Includes HDX analysis module (MassAnalyzer). |
| **Mass Spec Studio** | msstudio.ca | Most comprehensive integrated tool. AutoHX module for automated peptide search + deuteration. Trial version available (limited). |

### 2.5 Deuterium Uptake Calculation and Kinetics Fitting

#### Core Calculation Libraries

| Tool | Language | What It Computes |
|------|----------|-----------------|
| **PyHDX** | Python | RFU, DeltaG, DeltaDeltaG. SGD optimization via PyTorch. Residue-level resolution from peptide overlap. |
| **HDXrate** | Python | Intrinsic exchange rates (k_int) from sequence. Uses Bai/Englander parameters. |
| **DECA** | Python | Automated centroid fitting, uptake curves, differential plots. |
| **The Deuterium Calculator** | Python | Uptake values, back-exchange correction, statistical comparison. |
| **HRaDeX** | R | High-resolution uptake rates via weighted averaging of overlapping peptides. |
| **PIGEON-FEATHER** | Python | Residue-level DeltaG_op via Bayesian/optimization approach. |

#### Kinetics Fitting Approaches

- **Simple exponential:** Fit uptake vs time to `D(t) = N * (1 - exp(-k*t))`. Works for EX2 regime.
- **Stretched exponential:** `D(t) = N * (1 - exp(-(k*t)^beta))`. Accounts for heterogeneous exchange rates.
- **Multi-exponential:** Sum of exponentials for fast/medium/slow exchanging populations.
- **Bimodal fitting:** Required for EX1 kinetics. Fit isotopic envelope as sum of two binomial distributions.
- **Bayesian approaches:** ReX, PyHDX use Bayesian/optimization methods for residue-level deconvolution.

### 2.6 Visualization and Comparison Tools

| Tool | Type | Key Plots |
|------|------|-----------|
| **HaDeX** | R + Web | Uptake curves, differential plots, Woods plots, coverage maps, volcano plots |
| **HDX-Viewer** | Web | Interactive 3D protein structure colored by HDX data (NGL viewer) |
| **Deuteros 2.0** | MATLAB | Butterfly plots, significance maps, coverage plots |
| **HDXBoxeR** | R | Multi-state comparison, volcano plots, heatmaps, CI-based significance |
| **The Deuterium Calculator** | Python | Woods plots, uptake curves, statistical analysis |
| **spectrum_utils** | Python | MS/MS spectrum annotation and mirror plots |
| **matplotlib/plotly** | Python | Custom plotting (isotopic envelopes, uptake kinetics, heatmaps) |

---

## 3. Thermo .RAW File Structure and Access

### 3.1 File Format Overview

Thermo .RAW files are **proprietary binary files** with no publicly available specification. Access is provided through Thermo Fisher's **RawFileReader** .NET assemblies. Key facts:

- Binary format, not human-readable
- Contains all data from a single LC-MS/MS acquisition
- Includes: scan data (spectra), instrument method, tune data, status log, error log, UV/PDA data (if applicable)
- Two levels of abstraction: low-level binary structures and high-level API access

### 3.2 Scan Organization (MS1 vs MS2)

**Scan numbering:** Scans are numbered sequentially starting from 1. Each scan has:
- **Scan number** (integer, sequential)
- **Retention time** (minutes from injection)
- **MS level** (1 for MS1/survey, 2 for MS2/fragmentation, etc.)
- **Filter string** describing the scan parameters

**Filter string format** (e.g., `FTMS + p NSI Full ms [350.00-1800.00]`):
- Analyzer type: `FTMS` (Orbitrap), `ITMS` (ion trap)
- Polarity: `+` or `-`
- Scan mode: `Full ms` (MS1), `Full ms2` (MS2), `SIM`, `SRM`
- Mass range: `[low-high]`
- For MS2: includes precursor m/z and fragmentation type (CID, HCD, ETD, EThcD)

**Typical HDX-MS acquisition:**
- MS1 scans: profile mode, high resolution (e.g., 120,000 at m/z 200)
- MS2 scans: centroid mode (for peptide ID on undeuterated samples)
- Data-dependent acquisition (DDA) or data-independent acquisition (DIA)
- MS1 scans interleaved with MS2 scans during DDA

**Scan hierarchy in DDA:**
```
Scan 1: MS1 (survey scan, profile)
Scan 2: MS2 (fragmentation of top-N precursor #1, centroid)
Scan 3: MS2 (fragmentation of top-N precursor #2, centroid)
...
Scan N+1: MS2 (fragmentation of top-N precursor #N, centroid)
Scan N+2: MS1 (next survey scan)
...
```

### 3.3 Retention Time Structure

- Retention time is recorded for each scan in minutes (float)
- Not uniformly spaced -- depends on scan speed, number of MS2 scans between MS1 scans, lock mass calibrations
- Typical MS1 scan interval: 0.5-3 seconds (depends on resolution and cycle time)
- Total run time for HDX: typically 5-15 minutes (short gradients to minimize back-exchange)
- Retention time can be queried by scan number or used to find scans within a time window

### 3.4 Extracting XICs (Extracted Ion Chromatograms)

An XIC shows the intensity of a specific m/z value (or m/z range) across all scans over time:

**Parameters for XIC extraction:**
- Target m/z (center of extraction window)
- Mass tolerance (e.g., +/- 10 ppm or +/- 0.02 Da)
- MS level filter (typically MS1 only)
- Retention time range (optional, to limit extraction window)

**Methods to extract XICs:**

1. **ThermoRawFileParser** (direct from RAW):
   ```bash
   ThermoRawFileParser -i file.raw -x xic_input.json -o output/
   ```
   Where `xic_input.json` specifies m/z values and tolerances.

2. **pyrawr** (Python, via ThermoRawFileParser):
   ```python
   from pyrawr import RawFile
   raw = RawFile("file.raw")
   xic = raw.get_xic(mz=500.25, tolerance=10, tolerance_unit="ppm")
   ```

3. **fisher-py** (Python, direct .NET):
   ```python
   from fisher_py import RawFile
   raw = RawFile("file.raw")
   # Access chromatogram data programmatically
   ```

4. **pyopenms** (from mzML):
   ```python
   import pyopenms
   exp = pyopenms.MSExperiment()
   pyopenms.MzMLFile().load("file.mzML", exp)
   # Iterate spectra, extract intensities at target m/z
   ```

5. **pymzml** (from mzML):
   ```python
   import pymzml
   run = pymzml.run.Reader("file.mzML")
   for spectrum in run:
       if spectrum.ms_level == 1:
           mz_array = spectrum.mz
           intensity_array = spectrum.i
           # Extract intensity at target m/z
   ```

6. **pyteomics** (from mzML):
   ```python
   from pyteomics import mzml
   for spectrum in mzml.MzML("file.mzML"):
       if spectrum['ms level'] == 1:
           mz = spectrum['m/z array']
           intensity = spectrum['intensity array']
   ```

### 3.5 Python Libraries for Reading RAW/mzML on Linux

**Direct RAW file access on Linux:**

| Library | Mechanism | Requirements |
|---------|-----------|-------------|
| fisher-py | pythonnet -> .NET RawFileReader | .NET 8.0 SDK installed on Linux |
| pyrawr | subprocess -> ThermoRawFileParser | Mono or .NET runtime + ThermoRawFileParser binary |

**Recommended two-step approach (most reliable on Linux):**

1. **Convert:** `ThermoRawFileParser -i input.raw -o output/ -f 2` (format 2 = mzML)
2. **Read in Python:** Use pyteomics, pymzml, or pyopenms to read the mzML

This two-step approach avoids .NET runtime dependencies in the Python process and is the most robust path for Linux-based pipelines.

### 3.6 Aggregated Signal and Time Slicing

For HDX-MS, what matters is extracting the isotopic envelope of each peptide at each time point:

1. **Find the peptide's retention time** from the undeuterated peptide ID
2. **Extract MS1 scans** within a retention time window around the peptide's elution peak (e.g., +/- 0.2 min)
3. **Sum or average** the MS1 spectra in that window to improve signal-to-noise
4. **Extract the isotopic envelope** at the peptide's expected m/z (for each charge state)
5. **Calculate centroid mass** of the isotopic envelope
6. **Compare** to undeuterated centroid to get deuterium uptake

---

## 4. Peptide Mapping for HDX-MS

### 4.1 From Protein Sequence to Pepsin Digest Prediction

**Key challenge:** Pepsin is a non-specific protease. Unlike trypsin (cleaves after K/R), pepsin's cleavage sites cannot be precisely predicted from sequence alone.

**Pepsin cleavage preferences:**
- Preferentially cleaves at hydrophobic residues: F, L, W, Y (in P1 or P1' position)
- Also cleaves at M, A, E, D with lower frequency
- Avoids cleaving before or after R, K, H, G, P
- At pH > 2: broader specificity
- Cleavage is highly reproducible under consistent conditions (pH, temperature, time, enzyme:substrate ratio)

**In silico digestion tools:**

| Tool | Type | Pepsin Support | Notes |
|------|------|----------------|-------|
| **pyteomics.parser.cleave()** | Python function | Yes (via expasy_rules or custom regex) | Built-in expasy_rules dict includes pepsin rules. Can specify missed cleavages. Best for Python integration. |
| **PeptideCutter** (ExPASy) | Web tool | Yes | Predicts cleavage sites for 36+ proteases/chemicals. Web-only, no programmatic API. |
| **Rapid Peptides Generator (RPG)** | Web/standalone | Yes | Supports multiple simultaneous proteases and custom rules. Pepsin at pH > 2 with ~20% miscleavage. |
| **Protein Cleaver** | R Shiny app | Yes (via cleaver R package) | 36 proteolytic enzymes. Interactive web interface. |
| **DeepDigest** | Python (ML) | No (8 proteases, not pepsin) | Deep learning approach; not directly applicable to pepsin. |

**Practical approach for HDX-MS:**
In practice, pepsin digestion prediction is approximate. The real peptide map is determined experimentally from the undeuterated MS/MS data. In silico prediction is useful for:
- Estimating expected peptide sizes (typically 5-25 residues)
- Estimating sequence coverage
- Planning the search space for database searching

### 4.2 MS/MS Spectra to Peptide Identification

#### Step-by-step process:

1. **Acquire MS/MS data** on undeuterated protein digest (no deuterium, so peptide masses are unambiguous)

2. **Convert RAW to mzML:**
   ```bash
   ThermoRawFileParser -i undeuterated.raw -f 2 -o ./mzml/
   ```

3. **Create FASTA database** containing the target protein sequence (and optionally common contaminants)

4. **Run search engine** with non-specific or semi-specific cleavage:

   **Using Sage (fastest, open-source):**
   ```json
   {
     "database": {
       "fasta": "protein.fasta",
       "enzyme": {"cleave_at": "", "restrict": "", "missed_cleavages": 0, "min_len": 5, "max_len": 50}
     },
     "precursor_tol": {"ppm": [-10, 10]},
     "fragment_tol": {"ppm": [-20, 20]}
   }
   ```
   Note: Setting empty `cleave_at` enables non-specific search. For a single-protein database, this is computationally tractable.

   **Using Comet:**
   ```
   num_enzyme_termini = 0          # non-specific
   search_enzyme_number = 0        # no enzyme
   peptide_mass_tolerance = 10     # ppm
   ```

   **Using sagepy (Python API):**
   ```python
   import sagepy
   # Configure non-specific search against target protein
   # Returns PSMs with scores, q-values
   ```

5. **Parse search results** with pyteomics:
   ```python
   from pyteomics import pepxml, mzid
   # Read pepXML or mzIdentML output
   psms = list(pepxml.read('results.pep.xml'))
   # Filter by q-value < 0.01 (1% FDR)
   ```

6. **Build peptide list** with retention times, charge states, and m/z values

#### Key considerations for HDX-MS peptide ID:

- Search only the target protein (small database = fast, even with non-specific cleavage)
- Use tight precursor mass tolerance (5-10 ppm for Orbitrap)
- No enzyme specificity or semi-specific (pepsin rules)
- Variable modifications: typically only oxidation (M)
- Fixed modifications: none (no alkylation, since no reduction/alkylation in HDX workflow)
- Charge states: typically 2+ to 5+
- Filter for high confidence: q-value < 0.01, reasonable peptide length (5-50 residues)

### 4.3 Matching Experimental Peptides to Theoretical Digest

Once peptides are identified by MS/MS, they need to be validated and organized:

1. **Map peptides to protein sequence:**
   ```python
   from pyteomics import parser
   protein_seq = "MKTAYIAKQRQISFVKSHFSRQ..."

   # For each identified peptide, find its position
   for peptide in identified_peptides:
       start = protein_seq.find(peptide.sequence)
       end = start + len(peptide.sequence)
       # Record: (start, end, sequence, charge, rt, m/z)
   ```

2. **Calculate sequence coverage:**
   - Mark each residue as covered/uncovered
   - Typical HDX coverage: 85-98% with good pepsin digestion
   - Gaps are common around prolines and highly charged regions

3. **Assess redundancy (overlap):**
   - Multiple overlapping peptides covering the same region enable residue-level resolution
   - More overlap = higher spatial resolution in HDX measurements
   - Tools like PyHDX, PIGEON-FEATHER, and ReX exploit this overlap

4. **Build the peptide map:**
   - For each peptide: (start_residue, end_residue, sequence, charge_state, retention_time, theoretical_mz, max_deuterium_uptake)
   - This map is used to extract deuterium uptake from all deuterated time-point samples

5. **Calculate theoretical properties for each peptide:**
   ```python
   from pyteomics import mass

   # Theoretical m/z
   mz_theoretical = mass.calculate_mass(sequence=peptide_seq, charge=z)

   # Maximum exchangeable amides
   n_amides = len(peptide_seq) - 1  # minus N-terminal
   n_prolines = peptide_seq.count('P')
   max_uptake = n_amides - n_prolines

   # Intrinsic exchange rate (from HDXrate or Bai/Englander parameters)
   # Depends on sequence context, pH, temperature
   ```

### 4.4 Pyteomics In-Silico Digestion

Pyteomics includes built-in support for proteolytic cleavage:

```python
from pyteomics import parser

# Trypsin digestion (for reference)
peptides = parser.cleave("MKTAYIAKQRQISFVK", "trypsin", missed_cleavages=1)

# Pepsin digestion (using ExPASy rules)
# expasy_rules includes 'pepsin pH1.3' and 'pepsin pH>2'
pepsin_peptides = parser.cleave(sequence, parser.expasy_rules['pepsin pH>2'], missed_cleavages=2)

# Non-specific (generates all subsequences of given length range)
# More realistic for pepsin since rules are approximate
all_peptides = set()
for length in range(5, 30):
    for i in range(len(sequence) - length + 1):
        all_peptides.add(sequence[i:i+length])
```

---

## 5. Pipeline Architecture Recommendation

### 5.1 Proposed Open-Source HDX-MS Pipeline

```
Phase 1: RAW File Handling
    Thermo .RAW files
         |
         v
    ThermoRawFileParser (RAW -> mzML)
         |
         v
    mzML files (open format)

Phase 2: Peptide Identification (undeuterated samples only)
    mzML (undeuterated) + protein.fasta
         |
         v
    Sage (or Comet) -- non-specific cleavage search
         |
         v
    pepXML / mzIdentML / TSV results
         |
         v
    pyteomics -- parse, filter FDR < 1%
         |
         v
    Peptide Map: [(start, end, seq, charge, RT, m/z), ...]

Phase 3: Deuterium Uptake Extraction
    mzML (all timepoints) + Peptide Map
         |
         v
    For each peptide, each timepoint:
      - Extract MS1 scans in RT window
      - Sum/average spectra
      - Extract isotopic envelope at expected m/z
      - Calculate centroid mass
      - D_uptake = centroid(deuterated) - centroid(undeuterated)
         |
         v
    Uptake Table: peptide x timepoint -> D_uptake

Phase 4: Analysis and Visualization
    Uptake Table
         |
         +---> Back-exchange correction (using maxD)
         +---> RFU calculation
         +---> Kinetics fitting (exponential curves)
         +---> Differential analysis (state A vs state B)
         +---> Statistical significance (Woods plots, volcano plots)
         +---> Residue-level deconvolution (PyHDX / PIGEON-FEATHER)
         +---> Structural mapping (PyMOL / HDX-Viewer)
```

### 5.2 Recommended Python Package Stack

```
# Core dependencies
pip install pyteomics[XML]      # mzML reading, FASTA parsing, in-silico digestion, search result parsing
pip install pymzml              # Lightweight mzML reader (alternative)
pip install ms-deisotope        # Deisotoping and charge state deconvolution
pip install numpy scipy pandas  # Numerical computation
pip install matplotlib plotly   # Visualization

# RAW file handling
pip install pyrawr              # Python wrapper for ThermoRawFileParser
# Also install: ThermoRawFileParser (via bioconda or .NET binary)

# Peptide search (choose one)
pip install sagepy              # Python interface to Sage search engine
# Or: download Comet binary from comet-ms.sourceforge.net

# HDX-specific analysis
pip install pyhdx               # Residue-level analysis, DeltaG calculation
# pip install fisher-py         # Direct RAW access (optional, needs .NET SDK)

# Spectrum processing
pip install spectrum_utils      # MS/MS visualization
pip install matchms             # Spectral similarity
```

### 5.3 Key Design Decisions

1. **mzML as intermediate format:** Convert RAW to mzML early. All downstream tools read mzML. Avoids vendor lock-in and .NET dependencies in the main pipeline.

2. **Sage for peptide ID:** Fastest open-source search engine. MIT license. Python API via sagepy. Handles non-specific cleavage well against small databases.

3. **pyteomics as the backbone:** Reads mzML, parses search results (pepXML/mzIdentML), does in-silico digestion, calculates masses. Well-maintained, pure Python, pip-installable.

4. **Custom isotopic envelope extraction:** For deuterium uptake, we need to extract and analyze isotopic envelopes ourselves. ms_deisotope helps with deisotoping, but the centroid mass calculation for HDX is straightforward:
   - Find peaks in the isotopic envelope
   - Calculate weighted average m/z (centroid)
   - Compare to undeuterated centroid

5. **PyHDX for advanced analysis:** Once uptake values are calculated, PyHDX can derive residue-level resolution and protection factors.

### 5.4 Data Flow and File Formats

```
Input:
  - Thermo .RAW files (one per sample/timepoint)
  - protein.fasta (target protein sequence)
  - experiment.csv (sample metadata: condition, timepoint, replicate)

Intermediate:
  - .mzML files (converted from RAW)
  - search_results.pepxml or .tsv (peptide IDs)
  - peptide_map.csv (peptide list with positions, RT, m/z)

Output:
  - uptake_table.csv (peptide x timepoint deuterium uptake matrix)
  - kinetics_fits.csv (rate constants per peptide)
  - differential.csv (state comparison with p-values)
  - woods_plot.png/html (differential visualization)
  - coverage_map.png (sequence coverage)
  - structure_colored.pml (PyMOL script for structural mapping)
```

### 5.5 Challenges and Considerations

1. **Back-exchange during LC:** 25-35% of incorporated deuterium is lost during chromatography. Must be measured with maxD control and corrected.

2. **Isotopic envelope overlap:** At high charge states or for co-eluting peptides, isotopic envelopes can overlap. Deisotoping (ms_deisotope) helps, but careful XIC extraction at specific m/z is essential.

3. **Non-specific cleavage:** Pepsin's non-specificity means the peptide search space is large. Using a single-protein FASTA database keeps this manageable.

4. **Bimodal distributions (EX1):** Standard centroid analysis fails for EX1 kinetics. Detecting and handling bimodal isotopic envelopes requires spectral fitting, not just centroid calculation.

5. **Reproducibility:** All steps (labeling, quench, digestion, LC) must be highly reproducible. Pipeline should include QC checks (RT stability, intensity consistency, peptide recovery).

6. **Retention time alignment:** Across many samples, RT may drift. Alignment between deuterated and undeuterated runs is critical for correct peptide extraction.

---

## Appendix A: Key References

### Review Articles
- Stofella et al., "Computational Tools for Hydrogen-Deuterium Exchange Mass Spectrometry Data Analysis," Chemical Reviews, 2024, 124(21), 12242-12263. https://pubs.acs.org/doi/10.1021/acs.chemrev.4c00438
- Masson et al., "Recommendations for performing, interpreting and reporting HDX-MS experiments," Nature Methods, 2019. https://www.nature.com/articles/s41592-019-0459-y
- Engen et al., "Advances in Hydrogen/Deuterium Exchange Mass Spectrometry," Chemical Reviews, 2022. https://pubs.acs.org/doi/10.1021/acs.chemrev.1c00279

### Software Papers
- ThermoRawFileParser: Hulstaert et al., J Proteome Res, 2020. https://pmc.ncbi.nlm.nih.gov/articles/PMC7116465/
- PyHDX: Smit et al., Anal Chem, 2021. https://pubs.acs.org/doi/abs/10.1021/acs.analchem.1c02155
- Sage: Lazear, J Proteome Res, 2023. https://pubs.acs.org/doi/abs/10.1021/acs.jproteome.3c00486
- DECA: Lumpkin & Komives, MCP, 2020. https://www.sciencedirect.com/science/article/pii/S1535947620316534
- The Deuterium Calculator: Wu Lab, J Proteome Res, 2023. https://pubs.acs.org/doi/abs/10.1021/acs.jproteome.2c00558
- HDXBoxeR: Kajano et al., Bioinformatics, 2024. https://academic.oup.com/bioinformatics/article/40/8/btae479/7723990
- HRaDeX: hadexversum, J Proteome Res, 2024. https://pubs.acs.org/doi/10.1021/acs.jproteome.4c00700
- PIGEON-FEATHER: Glasgow Lab, Nature Chemical Biology, 2025. https://www.nature.com/articles/s41589-025-02049-1
- ReX: Crook et al., Communications Chemistry, 2025. https://www.nature.com/articles/s42004-025-01719-4

### Resource Lists
- HDX-MS Software Resources (curated list): https://github.com/hadexversum/HDX-MS-resources
- HDX-MS community resources: http://hdxms.net/resources/

## Appendix B: Quick-Start Installation (Linux)

```bash
# 1. Install ThermoRawFileParser (via bioconda)
conda install -c bioconda thermorawfileparser

# 2. Install Python packages
pip install pyteomics[XML] pymzml ms-deisotope numpy scipy pandas matplotlib plotly

# 3. Install search engine
pip install sagepy
# Or download Sage binary: https://github.com/lazear/sage/releases

# 4. Install HDX analysis
pip install pyhdx

# 5. Optional: direct RAW access
pip install fisher-py pyrawr

# 6. Verify
python -c "import pyteomics; import pymzml; import sagepy; print('All packages loaded')"
```
