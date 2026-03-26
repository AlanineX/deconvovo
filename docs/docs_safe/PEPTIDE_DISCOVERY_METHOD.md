# Peptide Discovery Method

Detailed documentation for `scripts/hdx/discover_peptides.py`.

## What it does

Finds all detectable pepsin-digest peptides from an HDX-MS undeuterated mapping run by searching MS1 spectra for isotopic envelopes matching theoretical peptide masses.

**This is NOT an MS/MS-based peptide identification.** There is no fragmentation, no b/y ion matching, no database search engine (Sage/Comet/MSFragger). It is purely an **MS1 precursor-level search** based on accurate mass and isotope pattern matching.

## Why MS1-only (not MS/MS)?

In HDX-MS, deuterium uptake is measured from MS1 isotopic envelope shifts — not from MS/MS fragmentation. The peptide mapping run is typically acquired with DDA (MS1 + MS2), but we only need to confirm that a peptide's isotopic envelope is detectable in MS1. If it's there, we can track its mass shift across timepoints.

For true peptide identification with sequence confirmation, you'd use MS/MS search engines (Sage, Comet, MSFragger) with fragment ion matching (b/y ions for CID/HCD). We attempted this with sagepy v0.4.4 but hit a FASTA parser bug. The MS1 approach is a pragmatic alternative that works well for single-protein experiments where the sequence is known.

## Algorithm

### Step 1: In-silico pepsin digest

Generate ALL possible sub-peptides of the protein sequence.

```
For start in [0, 1, ..., N-1]:
    For end in [start+min_len, ..., min(start+max_len, N)]:
        peptide = sequence[start:end]
```

**Why all sub-peptides?** Pepsin is non-specific — it cleaves at many positions (preferentially after F, L, W, Y, but not exclusively). Rather than predicting cleavage sites, we enumerate all possible peptides and let the MS1 data tell us which ones exist.

**Mass calculation:** `pyteomics.mass.fast_mass(sequence, ion_type="M", charge=0)` — computes monoisotopic neutral mass from amino acid composition. This is the standard approach using the monoisotopic masses of each residue (e.g., G=57.02146, A=71.03711, etc.).

**Complexity:** For a 129-residue protein with min_len=5, max_len=40: ~4,000 candidate peptides. With 6 charge states each: ~15,000 targets to search.

### Step 2: Compute theoretical m/z

For each peptide at each charge state z:

```
mz_monoisotopic = (neutral_mass + z × PROTON) / z
```

Where `PROTON = 1.007276466812 Da`.

Only targets with `300 < mz < 2000` are kept (typical MS1 acquisition range).

### Step 3: Search MS1 spectra

For each MS1 scan in the mzML file, for each target:

1. **Compute tolerance window:** `tol = mz_mono × ppm / 1e6`
2. **Check 5 isotope peaks** (M+0 through M+4):
   - `mz_M+i = mz_mono + i × NEUTRON / z` where `NEUTRON = 1.0033548378 Da`
   - Search the spectrum for a peak within ±tol of each theoretical isotope
   - Record the maximum intensity within the tolerance window
3. **Require minimum matches:** At least `min_isotopes` (default 3) of the 5 peaks must be present
4. **Score the match:**

```
score = Σ (weight_i × matched_intensity_i)  for i where intensity > 0

weights = [1.0, 0.8, 0.6, 0.4, 0.2]  (M+0 weighted highest)
```

5. **Keep best hit:** Per (sequence, charge), keep the scan with the highest score across all MS1 scans

### Step 4: Deduplicate

A peptide may be found at multiple charge states. Keep only the charge state with the highest score per unique sequence.

### Step 5: Compare against known targets

If a known peptide list is provided, check which known sequences were found and which new sequences were discovered.

## Output columns explained

| Column | Type | Description |
|--------|------|-------------|
| `sequence` | str | Amino acid sequence of the discovered peptide |
| `start` | int | 1-based start position in the protein |
| `end` | int | End position (inclusive) |
| `charge` | int | Best-scoring charge state |
| `mz_mono` | float | Theoretical monoisotopic m/z at that charge: `(mass + z×H) / z` |
| `neutral_mass` | float | Monoisotopic neutral mass from pyteomics |
| `score` | float | Weighted sum of matched isotope peak intensities (see formula above) |
| `matched_isotopes` | int | How many of the 5 theoretical isotope peaks were found (3-5) |
| `mono_intensity` | float | Intensity of the M+0 (monoisotopic) peak. Can be 0 if M+0 was not matched but M+1,M+2,M+3 were (common for larger peptides where M+0 is weak) |
| `scan` | int | MS1 scan number where the best match was found |

**What the score is NOT:**
- Not a statistical confidence (no p-value, no FDR)
- Not a fragmentation score (no b/y ions)
- Not a spectral library match
- Not normalized — absolute intensity units, depends on instrument sensitivity

**What `matched_isotopes` means:**
- 3 = minimum (marginal detection)
- 4 = good
- 5 = excellent (all isotope peaks found)
- Does NOT verify the isotope intensity pattern shape (e.g., doesn't check that M+1 > M+0 for large peptides as expected from the averagine model)

**Why `mono_intensity` can be 0:**
For peptides larger than ~1500 Da, the monoisotopic peak (M+0) is often the weakest in the envelope. The most abundant peak may be M+1 or M+2. The algorithm finds M+1, M+2, M+3 but not M+0, so `mono_intensity = 0` while the peptide is still confidently detected.

## What ion types are involved?

**None.** This is NOT a fragmentation-based search. There are no b ions, y ions, c ions, or z ions. All matching is done at the **intact precursor** level in MS1 — the isotopic envelope of the whole peptide `[M + zH]^z+`.

If you wanted b/y ion confirmation, you would need:
1. An MS/MS search engine (Sage, Comet, MSFragger)
2. Process the MS2 spectra (not MS1)
3. Score fragment ion matches against theoretical fragmentation patterns
4. Apply FDR control via target-decoy competition

## Libraries used

| Library | What it does here |
|---------|-------------------|
| `pyteomics.mass` | `fast_mass()` — computes monoisotopic neutral mass from amino acid sequence |
| `pyteomics.mzml` | Reads mzML files, iterates over spectra, extracts m/z and intensity arrays |
| `numpy` | Array operations for peak matching (vectorized tolerance check) |
| `pandas` | DataFrame for candidate management and CSV output |

**No other libraries.** No machine learning, no spectral libraries, no external databases.

## Physical constants

| Constant | Value | Source |
|----------|-------|--------|
| PROTON | 1.007276466812 Da | CODATA 2018 proton mass |
| NEUTRON | 1.0033548378 Da | ¹³C - ¹²C mass difference (not actual neutron mass) |

The "NEUTRON" constant is actually the mass difference between ¹³C and ¹²C isotopes, which determines the spacing between isotope peaks in a mass spectrum. It's called NEUTRON by convention in proteomics software.

## Hardcoded values

| Item | Value | Why | Should it be a parameter? |
|------|-------|-----|--------------------------|
| `TTR_SEQUENCE` | 129 aa | Default protein for this project | **Yes** — should use `--fasta` instead. Currently you can override with `--sequence` CLI arg |
| `n_isotopes` | 5 | Number of isotope peaks to check (M+0 to M+4) | Could be, but 5 is standard |
| `min_isotopes` | 3 | Minimum matched peaks to accept a hit | Could be. Lower = more false positives, higher = miss weak signals |
| `weights` | [1.0, 0.8, 0.6, 0.4, 0.2] | Scoring weights per isotope position | Arbitrary. Down-weights higher isotopes. No theoretical basis for these specific values |
| m/z range | 300-2000 | Typical Orbitrap MS1 range | Should match acquisition settings |

## Tunable parameters (CLI)

| Parameter | Default | Effect of changing |
|-----------|---------|-------------------|
| `--ppm` | 15.0 | **Lower (5-10):** fewer false matches, may miss if calibration drifts. **Higher (20-30):** more matches, more false positives. 15 ppm is typical for Orbitrap. |
| `--min-len` | 4 | **Higher (6-8):** fewer very short peptides (often noise). **Lower (3):** catches short fragments. |
| `--max-len` | 40 | **Higher (50):** catches longer peptides. **Lower (30):** faster, excludes large peptides. |
| `--charge-min` | 1 | Usually 1 for small peptides. Set to 2 if singly-charged noise is a problem. |
| `--charge-max` | 6 | For peptides up to ~5 kDa, z=1-6 covers most cases. Increase for larger peptides. |

## What's NOT reported (known gaps)

| Missing | Why it matters | How to add |
|---------|---------------|------------|
| **ppm error** | Can't assess mass accuracy per hit | Compute `(mz_observed - mz_theoretical) / mz_theoretical × 1e6` for each isotope peak |
| **Isotope pattern fit** | Score doesn't check if pattern shape matches averagine prediction | Compare observed vs theoretical isotope ratios |
| **Retention time** | Can't confirm elution behavior | Record RT from the best-matching scan |
| **FDR / decoy search** | No statistical confidence | Would need a decoy protein or scrambled peptides for target-decoy competition |
| **MS/MS confirmation** | No sequence-level evidence | Would need Sage/Comet search on MS2 data |
| **Back-exchange correction** | Irrelevant for mapping run (undeuterated) | Only needed for deuterated timepoints |

## Limitations

1. **False positives:** Any molecule in the sample with the right mass and isotope pattern will match. No sequence-level confirmation. For a single-protein digest this is usually fine; for complex mixtures it would be unreliable.

2. **No scoring normalization:** The score is raw intensity × weights. A peptide from a high-TIC scan scores higher than the same peptide from a low-TIC scan. This makes score comparison between peptides meaningless for abundance ranking.

3. **No decoy FDR:** Without a null model, we can't estimate the false discovery rate. The 2687 "discovered" peptides include an unknown number of false positives.

4. **Hardcoded sequence:** The TTR sequence is embedded in the script. Should accept a FASTA file instead.

5. **Brute-force search:** O(scans × targets × isotopes) — slow for large proteins. For TTR (129 aa, ~15K targets, ~3600 scans) it takes ~5 minutes. For a 500-residue protein it would be impractical without optimization (e.g., sorted m/z binary search).

## Recommended improvements

1. Add `--fasta` parameter, remove hardcoded sequence
2. Report ppm error per match: `ppm_error = (mz_best_match - mz_theoretical) / mz_theoretical × 1e6`
3. Report retention time of best-matching scan
4. Add isotope pattern scoring (cosine similarity against averagine model)
5. Add decoy search for FDR estimation (reverse sequence or scrambled peptides)
6. Use binary search for peak matching (currently linear scan)
7. Consider MS/MS search via Sage when the FASTA parser bug is fixed (sagepy > 0.4.4)
