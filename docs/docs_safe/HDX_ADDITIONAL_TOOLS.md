# Definitive Survey: HDX-MS Software Tools

**Research date:** 2026-03-25
**Purpose:** Comprehensive survey of all available tools for building an mzML-to-uptake HDX-MS pipeline

---

## Executive Summary

The HDX-MS software landscape is sharply divided:

1. **Raw data processing (mzML/raw to peptide uptake tables)** is dominated by commercial
   vendor software (Waters DynamX, Thermo HDExaminer/BioPharma Finder). Only three open-source
   tools attempt this: **TheDeuteriumCalculator** (Python), **ExMS2** (MATLAB standalone), and
   **OligoR** (R, oligonucleotide-focused). Mass Spec Studio (now CHRONECT Studio) was the
   most comprehensive open-source option but was acquired by Trajan in 2024 and commercialized.

2. **Downstream analysis (pre-processed CSV to statistics/visualization/protection factors)**
   has abundant open-source options in Python and R.

3. **No established Python package provides a complete, production-quality pipeline from raw
   mzML to deuterium uptake.** This is the gap we would fill.

Key review: Stofella et al., "Computational Tools for Hydrogen-Deuterium Exchange Mass
Spectrometry Data Analysis," *Chemical Reviews* 124(21), 2024. Covers 30+ tools across all
workflow stages.

Community resource list: https://github.com/hadexversum/HDX-MS-resources

---

## Table of Contents

1. [Tools That Process Raw/mzML Data](#1-tools-that-process-rawmzml-data)
2. [Downstream Analysis Tools (pre-processed input)](#2-downstream-analysis-tools-pre-processed-input)
3. [Protection Factor / Residue-Level Resolution Tools](#3-protection-factor--residue-level-resolution-tools)
4. [Multimodal / EX1 Kinetics Tools](#4-multimodal--ex1-kinetics-tools)
5. [Visualization-Only Tools](#5-visualization-only-tools)
6. [Simulation / Prediction Tools](#6-simulation--prediction-tools)
7. [General MS Libraries Useful for HDX](#7-general-ms-libraries-useful-for-hdx)
8. [Deep Dive: TheDeuteriumCalculator](#8-deep-dive-thedeuteriumcalculator)
9. [Deep Dive: ExMS2](#9-deep-dive-exms2)
10. [Deep Dive: ms_deisotope for HDX](#10-deep-dive-ms_deisotope-for-hdx)
11. [Commercial / Vendor Tools (for reference)](#11-commercial--vendor-tools-for-reference)
12. [Summary Table](#12-summary-table)

---

## 1. Tools That Process Raw/mzML Data

These are the rare tools that can start from mzML (or raw) files rather than requiring
pre-processed peptide tables from vendor software.

### TheDeuteriumCalculator

- **Language:** Python 3.8+
- **Input:** mzML files (converted from any vendor format via MSConvert)
- **Requires:** Pre-identified peptide list (from database search, not built-in)
- **Computes:** Weighted centroid mass, deuterium uptake, Woods plots, IC-HDX standardized output
- **GitHub:** https://github.com/OUWuLab/TheDeuteriumCalculator
- **Citation:** J. Proteome Res. 2023, doi:10.1021/acs.jproteome.2c00558
- **pip installable:** No (clone and run script)
- **Dependencies:** pyteomics, lxml, numpy, scipy, pandas, matplotlib
- **License:** Open source
- **See:** [Section 8 for deep dive](#8-deep-dive-thedeuteriumcalculator)

### ExMS2

- **Language:** MATLAB (compiled standalone available -- no MATLAB license needed)
- **Input:** mzML files
- **Computes:** Peptide identification + validation, centroid uptake, bimodal/multimodal analysis, ETD/ECD fragment analysis
- **Download:** Contact Englander lab, hx2.med.upenn.edu
- **Citation:** Kan et al., Anal. Chem. 2019, doi:10.1021/acs.analchem.9b01682
- **pip installable:** No (MATLAB compiled)
- **License:** Available on request
- **See:** [Section 9 for deep dive](#9-deep-dive-exms2)

### OligoR

- **Language:** R (Shiny web app)
- **Input:** mzML files (converted via MSConvert from Thermo .raw or Agilent .d)
- **Computes:** HDX uptake for **oligonucleotides** (DNA, not proteins), bimodal distribution deconvolution, native MS kinetics
- **GitHub:** https://github.com/EricLarG4/OligoR
- **Citation:** Anal. Chem. 2023, doi:10.1021/acs.analchem.3c01321
- **pip installable:** No (R/Shiny)
- **Note:** Protein HDX is out of scope for this tool, but its mzML-reading approach (via mzR package) and bimodal deconvolution algorithm are worth studying

### Mass Spec Studio / CHRONECT Studio (formerly open, now commercial)

- **Language:** C# (.NET)
- **Input:** Raw MS data from multiple vendors
- **Computes:** Full pipeline: automated peptide search (HXpipe), differential analysis, multimodal analysis (AutoHX)
- **URL:** https://www.msstudio.ca/hdx/
- **Citation:** Structure 2014; AutoHX in Nature Communications 2024
- **Status:** Acquired by Trajan in 2024, renamed CHRONECT Studio / HDExaminer PRO. Now commercial.
- **Note:** Was considered "the most comprehensive" tool in the 2024 Chemical Reviews survey

### Hexicon 2

- **Language:** Unknown (standalone application)
- **Input:** Raw MS data
- **Computes:** Automated peptide identification via NITPICK peak detection, chromatogram alignment, deuteration distribution estimation via iterative deconvolution
- **URL:** https://hx2.mr.mpg.de/ (available on request)
- **Citation:** J. Am. Soc. Mass Spectrom. 2014, doi:10.1007/s13361-014-0850-y
- **Status:** Available on request from Max Planck Institute

### HDX Workbench

- **Language:** Java
- **Input:** Raw MS data (Thermo)
- **Computes:** Full pipeline: data extraction, peptide matching, uptake calculation, visualization
- **URL:** https://sourceforge.net/projects/hdxwb/
- **Citation:** J. Am. Soc. Mass Spectrom. 2012, doi:10.1007/s13361-012-0419-6
- **Status:** Legacy tool, no longer actively maintained

### mhdx_pipeline (Rocklin Lab)

- **Language:** Python + Snakemake
- **Input:** mzML files (from LC-IMS-MS experiments)
- **Computes:** Isotopic cluster identification via factorization, path optimization across timepoints
- **GitHub:** https://github.com/Rocklin-Lab/mhdx_pipeline
- **Note:** Specific to **IMS-MS** (ion mobility) HDX experiments on protein libraries, not general-purpose bottom-up HDX. Uses IMTBX and Grppr for isotope peak processing.

---

## 2. Downstream Analysis Tools (pre-processed input)

These tools consume peptide uptake tables (CSV) exported from vendor software and perform
statistical analysis, back-exchange correction, visualization, and comparison.

### PyHDX

- **Language:** Python
- **Input:** Pre-processed CSV peptide tables (D-uptake/exposure format)
- **Computes:** Gibbs free energy (deltaG) at single-residue resolution via regularized fitting
- **GitHub:** https://github.com/Jhsmit/PyHDX
- **pip installable:** Yes (`pip install pyhdx`)
- **Citation:** Anal. Chem. 2021
- **Note:** Does NOT read mzML. Excellent for downstream deltaG derivation from peptide-level data.

### DECA (Deuterium Exchange Correction and Analysis)

- **Language:** Python
- **Input:** CSV from DynamX, HDX Workbench, or generic format
- **Computes:** Back-exchange correction (global and local), overlapping peptide segmentation, statistical analysis, Woods plots
- **GitHub:** https://github.com/komiveslab/DECA
- **Citation:** Mol. Cell. Proteomics 2019, doi:10.1074/mcp.TIR119.001731
- **pip installable:** No (clone and run)
- **Group:** Komives lab (UCSD)

### HaDeX

- **Language:** R
- **Input:** DynamX 2.0/3.0 cluster data CSV
- **Computes:** Differential analysis, uncertainty intervals, quality control, publication figures
- **GitHub:** https://github.com/hadexversum/HaDeX
- **CRAN:** Yes (`install.packages("HaDeX")`)
- **Web server:** Available
- **Citation:** Bioinformatics 2020, doi:10.1093/bioinformatics/btaa587

### HaDeX2

- **Language:** R
- **GitHub:** https://github.com/hadexversum/HaDeX2
- **Note:** Next-generation version, under development

### HDXBoxeR

- **Language:** R
- **Input:** Pre-processed peptide uptake CSV
- **Computes:** Welch T-test, critical interval framework, volcano plots, heat maps, back exchange calculation, PyMOL script generation
- **GitHub:** https://github.com/mkajano/HDXBoxeR
- **CRAN:** Yes (`install.packages("HDXBoxeR")`)
- **Citation:** Bioinformatics 2024, doi:10.1093/bioinformatics/btae479

### Deuteros 2.0

- **Language:** Python (standalone app, no MATLAB required)
- **Input:** Pre-processed data (designed to pair with Waters DynamX output)
- **Computes:** Peptide-level uptake kinetics, differential analysis, Woods plots, volcano plots, back-exchange correction, PyMOL export
- **GitHub:** https://github.com/andymlau/Deuteros_2.0
- **Citation:** Bioinformatics 2019 (v1), 2020 (v2)
- **License:** Apache 2.0
- **pip installable:** No (standalone binary)

### Kingfisher

- **Language:** R (Shiny web app)
- **Input:** Exports from common software packages
- **Computes:** Statistical analysis, high-resolution representations
- **Web server:** https://kingfisher.wustl.edu/
- **Citation:** Protein Sci. 2025, doi:10.1002/pro.70096
- **Note:** Very recent (March 2025), from Washington University in St. Louis

### HydroBot

- **Language:** Unknown (standalone)
- **Input:** Pre-processed HDX-MS data
- **Computes:** Multistate analysis, k-means/hierarchical clustering, three-state coordinate construction, volcano plots, trajectory visualization
- **Citation:** J. Am. Soc. Mass Spectrom. 2026, 37(1):341-345
- **Note:** Very new (published 2026). First tool for multistate/nonequilibrium HDX analysis.

### HD-eXplosion

- **Language:** Web-based
- **Input:** Pre-processed peptide data
- **Computes:** Differential analysis, t-test with hybrid thresholding
- **GitHub:** https://github.com/HD-Explosion
- **Note:** Web server, combines t-tests with manual thresholding for multiplicity correction

### MEMHDX

- **Language:** R (web server)
- **Input:** Pre-processed peptide data
- **Computes:** Mixed-effects model accounting for technical and biological replicates; only tool using mixed models for HDX
- **Note:** Web server only; publication available but limited maintenance info

### HDXanalyzer

- **Language:** R/Python hybrid
- **Input:** Pre-processed data
- **Computes:** Linear regression model, differential analysis
- **Note:** Older tool, limited current maintenance

### hdxmsqc

- **Language:** R (Bioconductor)
- **Input:** Pre-processed HDX-MS data
- **Computes:** Quality control metrics, experimental design encoding, QC visualization
- **GitHub:** https://github.com/ococrook/hdxmsqc
- **Bioconductor:** Yes
- **Note:** Focused purely on QC, not analysis

### hdxms-datasets

- **Language:** Python
- **Input/Output:** Curated benchmark HDX-MS datasets (CSV peptide tables)
- **GitHub:** https://github.com/Jhsmit/hdxms-datasets
- **pip installable:** Yes (`pip install hdxms-datasets`)
- **Note:** Not an analysis tool -- provides standardized test datasets for validating other tools. From the same author as PyHDX.

---

## 3. Protection Factor / Residue-Level Resolution Tools

These tools take peptide-level uptake data and compute single-residue resolution protection
factors or free energies.

### PIGEON-FEATHER

- **Language:** Python
- **Input:** HXMS format files (can be generated from BioPharma Finder, HDExaminer, DynamX, HDX Workbench via PFLink converter)
- **Computes:** Free energies of opening (deltaGop) at single/near-single amino acid resolution
- **GitHub:** https://github.com/glasgowlab/PIGEON-FEATHER
- **Dependencies:** Python 3.11, MDAnalysis, Numba, PyOpenMS, hdxrate, bayesian_hdx (modified fork)
- **pip installable:** No (clone + local install or Docker)
- **Citation:** Nat. Chem. Biol. 2025, doi:10.1038/s41589-025-02049-1
- **Group:** Glasgow lab
- **Note:** State of the art for residue-level energetics from HDX-MS

### exPfact

- **Language:** Python + R (mclust)
- **Input:** Pre-processed peptide uptake data
- **Computes:** Protection factors at single-amide resolution using regularized fitting with smoothness constraint
- **GitHub:** https://github.com/skinnersp/exPfact (also https://github.com/pacilab/exPfact)
- **Citation:** Skinner et al. 2019, doi:10.1016/j.bpj.2019.02.024
- **License:** GPL v2 (academic); commercial license available
- **Group:** Paci lab

### RexMS

- **Language:** R
- **Input:** Pre-processed peptide-level data
- **Computes:** Residue-level deuterium uptake via Bayesian non-parametric change-point model (reversible jump MCMC)
- **GitHub:** https://github.com/ococrook/RexMS
- **Install:** `remotes::install_github("ococrook/RexMS")`
- **Citation:** Nat. Commun. Chem. 2025, doi:10.1038/s42004-025-01719-4
- **Note:** Provides statistical confidence (probability of change) at residue level

### bayesian_hdx

- **Language:** Python
- **Input:** Pre-processed HDX data
- **Computes:** Residue-level protection factors via Bayesian sampling, includes intrinsic rate calculation (Bai/Englander 1993)
- **GitHub:** https://github.com/salilab/bayesian_hdx
- **Group:** Sali lab (UCSF)
- **Note:** V2.0 calculates explicit protection factors per residue rather than differences between states

### HDXmodeller

- **Language:** Web server
- **Input:** Low-resolution peptide-level HDX-MS data
- **Computes:** High-resolution exchange rates at residue level with auto-validation (99% accuracy via covariance matrix method)
- **URL:** Online webserver
- **Citation:** Commun. Biol. 2021, doi:10.1038/s42003-021-01709-x

### HRaDeX

- **Language:** R (package + web server)
- **Input:** Pre-processed peptide-level uptake curves
- **Computes:** High-resolution deuterium uptake rates considering overlapping peptide trajectories; achieves 7.15% average RMSE in reconstruction
- **GitHub:** https://github.com/hadexversum/HRaDeX
- **Citation:** J. Proteome Res. 2025, doi:10.1021/acs.jproteome.4c00700
- **Note:** Very recent (April 2025)

### HDfleX

- **Language:** MATLAB
- **Input:** Pre-processed data at any resolution level (peptide, ETD fragment, intact protein)
- **Computes:** Merges bottom-up and top-down HDX data, stretched exponential fitting of uptake curves
- **Citation:** Anal. Chem. 2022, doi:10.1021/acs.analchem.1c05339
- **Note:** Unique in combining multiple resolution levels (peptide + ETD fragments)

### PFNet

- **Language:** Python (neural network)
- **Input:** HDX-MS data
- **Computes:** Protection factor prediction via deep learning; automatic correction of experimental data
- **GitHub:** https://github.com/glasgowlab/PFNet
- **Hugging Face:** https://huggingface.co/spaces/glasgow-lab/PFNet
- **Citation:** J. Am. Soc. Mass Spectrom. 2023, doi:10.1021/jasms.3c00285
- **Group:** Glasgow lab

---

## 4. Multimodal / EX1 Kinetics Tools

Tools specifically for detecting and analyzing bimodal isotope distributions (EX1 kinetics,
conformational heterogeneity).

### pyHXExpress

- **Language:** Python
- **Input:** Spectral data (primarily from HDExaminer exports; has some mzML reading capability via pyteomics but this is secondary)
- **Computes:** Multimodal fits of HDX-MS spectra (unimodal, bimodal, trimodal); automated high-throughput analysis of entire peptide pools
- **GitHub:** https://github.com/tuttlelm/pyHXExpress
- **pip installable:** Yes (`pip install pyhxexpress`)
- **Citation:** BioRxiv 2025 (Tuttle, Klevit, Guttman)
- **License:** GPL-3.0
- **Note:** Python implementation of HX-Express v3 fitting subroutines. Under active development (last edit April 2025). Can work from Jupyter notebooks.

### deMix

- **Language:** Java (standalone with GUI)
- **Input:** Raw isotopic distributions
- **Computes:** Bimodal distribution detection and deconvolution using binomial distribution model; robust "Matched Peak Count" metric for noisy data
- **GitHub:** https://github.com/HanyangBISLab/deMix (also https://github.com/seungjinna/deMix)
- **Citation:** Sci. Rep. 2019, doi:10.1038/s41598-019-39512-8
- **Note:** Recent version includes AlphaFold structure mapping

### HX-Express (v3)

- **Language:** Microsoft Excel VBA
- **Input:** Mass spectral peak data
- **Computes:** Deuterium uptake, peak width plots, double binomial distribution fitting for bimodal detection
- **Status:** Original tool from 2006, v3 latest. Being superseded by pyHXExpress (Python port).

### ExMS2 (multimodal module)

- **See [ExMS2 deep dive](#9-deep-dive-exms2).** Has a dedicated bimodal/multimodal analysis module that quantifies population fractions and uptake of each conformational state.

---

## 5. Visualization-Only Tools

### HDX-Viewer

- **Language:** JavaScript/Web
- **Input:** Pre-processed HDX data
- **Computes:** Interactive 3D visualization of HDX data mapped onto protein structures
- **GitHub:** https://github.com/david-bouyssie/hdx-viewer
- **Citation:** Bioinformatics 2019, doi:10.1093/bioinformatics/btz550

---

## 6. Simulation / Prediction Tools

### HDXer

- **Language:** Python
- **Input:** Molecular dynamics trajectories + experimental HDX data
- **Computes:** Predicted H-D exchange rates from MD simulations; maximum-entropy ensemble reweighting to fit simulation to experiment
- **GitHub:** https://github.com/Lucy-Forrest-Lab/HDXer
- **Citation:** Biophys. J. 2020; LiveCoMS 2023
- **pip installable:** No (git clone)
- **Note:** Not for experimental data processing; for interpreting HDX results with structural models

### HDXrate

- **Language:** Python
- **Input:** Protein sequence
- **Computes:** Intrinsic (unfolded) amide hydrogen exchange rates based on Bai/Englander/Connelly/Molday parameters
- **GitHub:** https://github.com/Jhsmit/HDXrate
- **pip installable:** Yes
- **Note:** Pure calculation tool. Essential component for computing protection factors (PF = k_intrinsic / k_observed). Used as dependency by PIGEON-FEATHER and others.

---

## 7. General MS Libraries Useful for HDX

These are not HDX-specific but provide the building blocks for an mzML-to-uptake pipeline.

### pyteomics

- **Language:** Python
- **pip:** Yes (`pip install pyteomics`)
- **Relevant features:** mzML parsing, mass calculations, isotopic distribution computation
- **GitHub:** https://github.com/levitsky/pyteomics
- **Note:** Used by TheDeuteriumCalculator as its mzML reader. Mature, well-maintained.

### ms_deisotope

- **Language:** Python
- **pip:** Yes (`pip install ms_deisotope`)
- **GitHub:** https://github.com/mobiusklein/ms_deisotope
- **Relevant features:** Deisotoping, charge state deconvolution, isotopic pattern fitting (averagine model), mzML reading (via ms_peak_picker + brainpy)
- **HDX-specific functionality:** None explicitly. No "deuterium" or "HDX" code paths. However, its isotopic pattern fitting is the core algorithm needed for HDX envelope matching. The averagine model would need modification to account for deuterium incorporation shifting the isotope envelope.
- **See:** [Section 10 for details](#10-deep-dive-ms_deisotope-for-hdx)

### pyOpenMS

- **Language:** Python (C++ bindings)
- **pip:** Yes (`pip install pyopenms`)
- **Relevant features:** Peak picking, centroiding, mzML I/O, isotopic pattern generation, feature detection
- **Note:** Very comprehensive but heavy dependency. Used by PIGEON-FEATHER.

### pymzML

- **Language:** Python
- **pip:** Yes (`pip install pymzml`)
- **GitHub:** https://github.com/pymzml/pymzML
- **Relevant features:** Lightweight mzML reader. MIT license. Minimal dependencies (numpy only).

### brainpy

- **Language:** Python (Cython)
- **pip:** Bundled with ms_deisotope
- **Relevant features:** Fast theoretical isotopic distribution calculation using the Baffling Recursive Algorithm for Isotopic distributioN calculations (BRAIN)

---

## 8. Deep Dive: TheDeuteriumCalculator

**Repository:** https://github.com/OUWuLab/TheDeuteriumCalculator
**Paper:** PMC12510632

### Architecture

Single Python script (~1400 lines): `TheDeuteriumCalculator_V1.3.py` plus `PARAMETERS.py`
for configuration. No package structure. No tests.

### Algorithm: mzML to Deuterium Uptake

**Step 1: mzML reading.** Uses pyteomics to parse mzML. Processes all scans.

**Step 2: Peak detection via sliding window.** Sums individual m/z peak intensities across
scan ranges using a sliding window approach. The window size and overlap ("slide fraction")
are tunable. Recommended slide fraction < 4.

**Step 3: Retention time matching.** Uses "accurate mass and time" (AMT) strategy:
- Default RT window: +/- 30 seconds (adjustable)
- Default mass tolerance: 10 ppm (adjustable)
- If multiple peaks match, selects highest intensity

**Step 4: Theoretical mass generation.** For each peptide, generates theoretical masses for
all possible deuteration states:
```
Tmass,n = m/z * z + n * (1.00628 Da)
```

**Step 5: Centroid mass calculation.** Intensity-weighted average:
```
WeightedAve.mass = [SUM(Em/z_n * I_n / I_T) - 1.00627 Da] * z
```
Deuterium uptake = deuterated weighted mass - undeuterated weighted mass.

**Step 6: Quality assessment.** Fits isotope distribution to Gaussian function via
least-squares. R^2 > 0.9 cutoff rejects 91% of bad spectra while keeping good ones. Low-R^2
peptides flagged for manual validation.

### Limitations

1. **Requires pre-identified peptide list** -- cannot discover peptides from mzML alone
2. **Gaussian quality filter fails on bimodal distributions** (EX1 kinetics, conformational heterogeneity)
3. **No back-exchange correction** -- relies on experimental design (subzero UPLC) rather than computation
4. **No charge state deconvolution** -- uses charge from input peptide list only
5. **Single monolithic script** -- difficult to reuse components
6. **No EIC extraction** -- works on summed spectra, not extracted ion chromatograms
7. **No statistical framework** -- downstream tools (DECA, Deuteros) needed

### Validation Data

- Anthrax vaccine protective antigen (PA) + monoclonal antibodies (purified system)
- E. coli lysate (complex mixture)
- Triplicate HDX runs confirmed reproducibility

### What We Can Learn From It

- The AMT matching strategy (RT window + mass tolerance + highest intensity selection) is simple and effective
- Gaussian R^2 quality scoring is a useful first filter
- The centroid mass formula is standard and correct
- pyteomics is sufficient for mzML reading (no need for heavier libraries)

---

## 9. Deep Dive: ExMS2

**Lab:** S. Walter Englander, University of Pennsylvania (hx2.med.upenn.edu)
**Paper:** Kan et al., Anal. Chem. 2019

### Why It Is Considered Gold Standard

1. **Integrated pipeline:** Only tool that does peptide identification, validation, uptake
   calculation, multimodal analysis, AND ETD/ECD fragment analysis in one package.

2. **12 quality tests for peptide validation.** ExMS2 applies 12 sequential checks to each
   candidate peptide match, including: mass accuracy, retention time consistency, isotope
   pattern quality, charge state validation, signal-to-noise, and more. This rigorous
   validation is unmatched by other tools.

3. **Bimodal/multimodal detection.** ExMS2 can detect when a peptide's isotope distribution
   is bimodal (indicating EX1 kinetics or conformational heterogeneity). It quantifies:
   - Relative population fractions from MS amplitudes
   - Deuterium uptake of each population separately
   This is critical for studying protein folding/unfolding.

4. **ETD/ECD module.** Can analyze electron transfer dissociation / electron capture
   dissociation fragments for single-residue resolution. This is unique -- most tools cannot
   handle top-down or middle-down HDX data.

5. **No MATLAB license required.** Despite being written in MATLAB, it compiles to a
   standalone executable using the free MATLAB Runtime.

### Key Algorithms Not Found in Other Tools

- Multi-level validation cascade (12 tests) for peptide identification confidence
- Automatic bimodal distribution fitting with population quantification
- ETD fragment-level deuterium localization
- Integration of bottom-up and top-down HDX workflows

### Why No Python Port Exists

- The MATLAB code is complex and tightly integrated
- No publicly available source repository (distributed as compiled binary)
- The Englander lab has not released source code for community reimplementation
- **pyHXExpress partially reimplements the multimodal fitting** but not the full pipeline

### What a Python Reimplementation Would Need

1. mzML reader (pyteomics or pymzML)
2. Peptide identification from MS/MS or AMT matching
3. Isotope envelope extraction and scoring
4. Multi-criteria validation framework (12+ tests)
5. Centroid mass + uptake calculation with back-exchange correction
6. Bimodal distribution detection and deconvolution
7. Statistical framework for differential analysis

---

## 10. Deep Dive: ms_deisotope for HDX

**Repository:** https://github.com/mobiusklein/ms_deisotope

### HDX-Specific Functionality

**None.** There is no code in ms_deisotope that references "deuterium," "HDX," or "hydrogen
exchange." The library is designed for standard proteomics deisotoping.

### Why It Is Still Relevant

ms_deisotope's core algorithm -- fitting theoretical isotopic distributions to experimental
peaks to determine charge state and monoisotopic mass -- is exactly what HDX analysis needs,
with one critical modification: **the theoretical isotopic distribution must be shifted to
account for deuterium incorporation.**

In standard proteomics:
- Theoretical pattern = natural isotope distribution of peptide formula
- Goal: find monoisotopic peak and charge state

In HDX:
- Theoretical pattern = natural distribution + N deuterium atoms incorporated
- Goal: find the value of N (deuterium uptake) that best fits the observed envelope
- The envelope shifts right and broadens with increasing deuteration

### What Would Need to Change

1. **Isotopic model:** Replace averagine with peptide-specific formula + variable deuterium count
2. **Fitting objective:** Instead of finding monoisotopic mass, find best-fit deuterium count
3. **Envelope matching:** Score theoretical vs observed across a range of deuteration states
4. **Bimodal awareness:** Handle distributions that are sums of two differently-deuterated populations

The `brainpy` isotopic distribution calculator within ms_deisotope could be adapted for this
by modifying the elemental composition to include variable deuterium substitution at each
exchangeable position.

---

## 11. Commercial / Vendor Tools (for reference)

| Tool | Vendor | Input | Notes |
|------|--------|-------|-------|
| DynamX | Waters | Waters .raw | Industry standard for Waters instruments. Exports CSV cluster data. |
| HDExaminer PRO | Trajan (formerly Sierra Analytics) | Thermo .raw, Waters .raw | Formerly AutoHX/Mass Spec Studio HDX module. Now commercial. |
| BioPharma Finder | Thermo Fisher | Thermo .raw | Integrated with MassAnalyzer. |
| PLGS | Waters | Waters .raw | ProteinLynx Global SERVER. For peptide identification. |
| Protein Metrics (Byos) | Dotmatics | Multiple vendor formats | Commercial, comprehensive. |

---

## 12. Summary Table

### Tools That Can Read mzML/Raw Data Directly

| Tool | Language | mzML | Raw | Peptide ID | Centroid | Bimodal | Back-Exch | Open Source | pip |
|------|----------|------|-----|------------|----------|---------|-----------|-------------|-----|
| TheDeuteriumCalculator | Python | Yes | No | Needs input list | Yes | No | No | Yes | No |
| ExMS2 | MATLAB | Yes | No | Yes (12 tests) | Yes | Yes | Yes | On request | No |
| OligoR | R | Yes | No | N/A (oligo) | Yes | Yes | N/A | Yes | No |
| Mass Spec Studio | C# | Yes | Yes | Yes (HXpipe) | Yes | Yes | Yes | Was free, now commercial | No |
| Hexicon 2 | ? | Yes | Yes | Yes (NITPICK) | Yes | Yes | ? | On request | No |
| HDX Workbench | Java | No | Yes (Thermo) | Yes | Yes | No | Yes | Yes (SourceForge) | No |
| mhdx_pipeline | Python | Yes | No | Factorization | Yes | No | No | Yes | No |

### Downstream Analysis Tools (pre-processed CSV input)

| Tool | Language | deltaG/PF | Statistics | Bimodal | Visualization | Open Source | pip/CRAN |
|------|----------|-----------|------------|---------|---------------|-------------|----------|
| PyHDX | Python | Yes (deltaG) | Basic | No | Yes (web) | Yes | pip |
| PIGEON-FEATHER | Python | Yes (deltaGop) | Yes | No | Yes | Yes | No |
| exPfact | Python/R | Yes (PF) | Yes | No | Basic | Yes (GPL) | No |
| RexMS | R | Yes (residue) | Bayesian MCMC | No | Yes | Yes | GitHub |
| bayesian_hdx | Python | Yes (PF) | Bayesian | No | Basic | Yes | No |
| DECA | Python | No | Yes | No | Yes | Yes | No |
| HaDeX | R | No | Yes (t-test) | No | Yes | Yes | CRAN |
| HDXBoxeR | R | No | Yes (Welch T) | No | Yes | Yes | CRAN |
| Deuteros 2.0 | Python | No | Yes | No | Yes | Yes (Apache) | No |
| Kingfisher | R | No | Yes | No | Yes (web) | Yes | No |
| HydroBot | ? | No | Yes | No | Yes | ? | No |
| HD-eXplosion | Web | No | Yes (hybrid) | No | Yes | Yes | No |
| MEMHDX | R/Web | No | Mixed model | No | Yes | ? | No |
| pyHXExpress | Python | No | No | Yes | Yes | Yes (GPL) | pip |
| deMix | Java | No | No | Yes | Yes | Yes | No |
| HRaDeX | R | Yes (rates) | Yes | No | Yes | Yes | No |
| HDXmodeller | Web | Yes (rates) | Auto-valid | No | Yes | ? | No |
| PFNet | Python | Yes (PF, ML) | ML-based | No | Yes | Yes | No |
| HDXer | Python | Simulation | No | No | Yes | Yes | No |
| HDXrate | Python | Intrinsic only | No | No | No | Yes | pip |
| hdxms-datasets | Python | N/A (data) | N/A | N/A | N/A | Yes | pip |
| hdxmsqc | R | No | QC only | No | Yes | Yes | Bioconductor |
| HDX-Viewer | Web | No | No | No | 3D viz only | Yes | No |

---

## Key Takeaways for Building an mzML-to-Uptake Pipeline

1. **The gap is real.** No Python tool provides a complete, production-quality mzML-to-uptake
   pipeline. TheDeuteriumCalculator is closest but requires a peptide list and has no bimodal
   handling, no back-exchange correction, and is a single unstructured script.

2. **The algorithm is well-understood.** The centroid mass calculation, AMT matching, and
   isotope envelope scoring are standard across all tools. The hard parts are: peptide
   identification without vendor software, bimodal detection, and robust quality filtering.

3. **Key building blocks exist in Python:**
   - `pyteomics` or `pymzML` for mzML reading
   - `brainpy` (from ms_deisotope) for theoretical isotopic distributions
   - `HDXrate` for intrinsic exchange rates
   - `scipy.optimize` for envelope fitting

4. **ExMS2's 12-test validation cascade** is the benchmark for peptide identification quality.
   Reimplementing even a subset of these tests would significantly improve on
   TheDeuteriumCalculator's Gaussian-R^2-only approach.

5. **Bimodal detection is increasingly expected.** pyHXExpress, ExMS2, deMix, and OligoR all
   handle it. A modern pipeline should detect EX1 kinetics.

6. **For Waters IM-MS data specifically**, the mhdx_pipeline from Rocklin Lab shows that
   combining HDX with ion mobility in Python is feasible, though their approach is specialized
   for protein libraries rather than single-protein conformational analysis.
