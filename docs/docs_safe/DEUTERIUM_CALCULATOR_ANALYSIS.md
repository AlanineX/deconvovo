# TheDeuteriumCalculator: Complete Source Code Analysis

Research date: 2026-03-25

**Repository:** https://github.com/OUWuLab/TheDeuteriumCalculator
**Paper:** Cupp-Sutton et al. (2023) "The Deuterium Calculator: An Open-Source Tool
for Hydrogen-Deuterium Exchange Mass Spectrometry Analysis" *J. Proteome Res.*
DOI: 10.1021/acs.jproteome.2c00558
**Version analyzed:** 1.3.0 (TheDeuteriumCalculator_V1.3.py, ~1400 lines)
**Author:** Thomas Welborn, University of Oklahoma (Wu Lab)
**Dependencies:** numpy, scipy, pandas, matplotlib, pyteomics, lxml

---

## 1. Input Handling

### 1.1 Three Required Input Files

1. **Identification mzML file** -- The LC-MS/MS run used for peptide database
   searching (NOT the experimental HDX data). This provides scan numbers and
   retention times for identified peptides. Read via `pyteomics.mzml.read()`.

2. **Identification CSV file** -- Database search results (e.g., from MS-GF+)
   with columns:
   - `ID` -- unique ascending integer
   - `ScanNum` -- scan number where peptide was identified
   - `Peptide` -- amino acid sequence (may include flanking residues with dots,
     e.g., `K.PEPTIDER.A`)
   - `Precursor` -- monoisotopic m/z value
   - `Charge` -- charge state integer

3. **Protein sequence file** -- Plain text file containing the full protein
   sequence. Used to map peptide start/end positions via substring matching.

### 1.2 How mzML Files Are Read

The `ExperimentalRun.read_mzml()` method uses pyteomics:

```python
def read_mzml(self, file: str):
    total = 0
    with mzml.read(file) as f:
        for scan in f:
            if scan["ms level"] == 1:
                total += 1
    count = 0
    with mzml.read(file) as f:
        for scan in f:
            if scan["ms level"] == 1:
                retention_time = (scan["scanList"]["scan"][0]["scan start time"])
                retention_time *= CON.MINUTES_TO_SECONDS
                count += 1
                self.process_scan(scan)
    self.set_tuple_dictionary(self.sliding_window(retention_time))
```

Key observations:
- Only MS1 scans are processed (filters `ms level == 1`)
- Two passes: first counts total scans, second processes them
- `scan start time` is in minutes; converted to seconds
- Each scan is passed to `process_scan()` which extracts (m/z, intensity) tuples
  above a noise threshold

### 1.3 Scan Processing and Noise Filtering

```python
def process_scan(self, scan):
    retention_time = (scan["scanList"]["scan"][0]["scan start time"]
                      * CON.MINUTES_TO_SECONDS)
    tuple_list = []
    for j in range(len(scan["m/z array"])):
        if scan["intensity array"][j] > CON.NOISE_LIMIT:
            tuple_list.append((scan["m/z array"][j],
                             scan["intensity array"][j]))
    self.all_peaks.append({"retention time": retention_time,
                          "tuple list": tuple_list})
```

**Noise filtering is a simple hard threshold** (`NOISE_LIMIT = 10000` by
default). Any peak below this intensity is discarded. No local noise estimation,
no SNR ratio -- just a flat cutoff.

### 1.4 Parameters (from PARAMETERS.py)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `NOISE_LIMIT` | 10000 | Hard intensity threshold for peak inclusion |
| `PPM_MATCH_TOLERANCE` | 10 | PPM tolerance for matching peaks to theoretical masses |
| `SLIDING_WINDOW_PPM_TOLERANCE` | 1 | PPM tolerance for merging peaks within sliding window |
| `SLIDING_WINDOW_SIZE` | 30 s | Width of sliding window in seconds |
| `SLIDE_FRACTION` | 3 | Window moves by SIZE/FRACTION each step (overlap) |
| `RETENTION_TOLERANCE` | 30 s | +/- RT window around identified peptide |
| `DEUTERIUM_RECOVERY_RATE` | 1 | Back-exchange correction factor (0-1) |
| `DEUTERIUM_FRACTION` | 0.833 | D2O fraction in labeling buffer |
| `DEUTERIUM_MASS_DIFFERENCE` | 1.00628 | Mass difference D-H in Da |
| `MASS_OF_WATER` | 18.01528 | |
| `MASS_OF_HYDROGEN` | 1.007276 | Proton mass |
| `RETENTION_SHIFT_SLOPE` | 1 | RT calibration slope |
| `RETENTION_SHIFT_INTERCEPT` | 0 | RT calibration intercept |

**Note:** The `PEPTIDE_MASS_DICTIONARY` uses **average masses** (not
monoisotopic), e.g., `'A': 71.0779`. This is unusual -- most tools use
monoisotopic masses. This affects the `set_average_mass()` calculation in the
Peptide class but NOT the centroid calculation itself (which uses experimental
m/z values directly).

### 1.5 Interactive Runtime Inputs

The program prompts for:
- Number of time points and their values (in seconds)
- Whether this is a differential experiment (two conditions vs one)
- Number of replicates per condition
- Path to each experimental mzML file

---

## 2. Centroid Mass Calculation -- THE KEY FUNCTION

### 2.1 The Core Algorithm: `Peptide.set_weighted_mass()`

This is the function that computes the centroid mass from the isotopic envelope:

```python
def set_weighted_mass(self):
    total_intensity = 0
    mass = 0
    for det in self._deuterium_dictionary.keys():
        mz, intensity, ppm_error = self.get_deuterium(det)
        total_intensity += intensity
        mass += intensity * mz
    if total_intensity == 0:
        mass = 0
    else:
        mass /= total_intensity
    self._weighted_mass_to_charge = mass
    weighted_mass = self._weighted_mass_to_charge - CON.MASS_OF_HYDROGEN
    self._mass_shift = weighted_mass * self._charge - self._average_mass
    if total_intensity == 0 or self._fit == "Insufficient Match":
        self._weighted_mass_to_charge = -1
        self._mass_shift = -1
```

**Breaking this down:**

1. **Intensity-weighted average m/z**: `C_mz = sum(I_n * mz_n) / sum(I_n)`
   - Iterates over ALL deuteration slots (0 to maxD)
   - Uses whatever peaks were matched (unmatched slots have intensity=0, so
     they contribute nothing)

2. **Convert to mass shift**:
   ```
   mass_shift = (C_mz - m_proton) * z - average_mass
   ```
   where `average_mass` is computed from the amino acid sequence using the
   average mass dictionary.

3. **Sentinel values**: If no peaks matched (`total_intensity == 0`) or the
   Gaussian fit failed (`"Insufficient Match"`), both `weighted_mass_to_charge`
   and `mass_shift` are set to -1.

### 2.2 From Paper (Equation 3)

```
Weighted Ave. mass = [sum_{n=0}^{Dmax} (Em/z_n * I_n / I_T) - 1.00627 Da] * z
```

This matches the code: subtract proton mass, multiply by charge.

### 2.3 How the Envelope Boundaries Are Defined

The number of peaks in the isotope envelope is **NOT adaptive** -- it is
determined by `sequence_to_max_deuterium()`:

```python
def sequence_to_max_deuterium(sequence: str):
    max_deuterium = len(sequence) - 2
    for letter in sequence:
        if letter.lower() == 'p':
            max_deuterium -= 1
    return int(max_deuterium * CON.DEUTERIUM_RECOVERY_RATE
               * CON.DEUTERIUM_FRACTION) + 1
```

This computes:
- Start with `len(sequence) - 2` (removes N-terminal and amide backbone
  assumption: first two residues don't exchange at backbone amide)
- Subtract 1 for each proline (no backbone amide NH)
- Multiply by recovery rate and D2O fraction
- Add 1 for the M+0 peak
- Truncate to integer

**So the envelope width is fixed per peptide** based on sequence length, proline
count, and D2O fraction. The code then looks for peaks at:
```
m/z_theoretical(n) = (peptide_mass + n * 1.00628) / charge
```
for n = 0, 1, 2, ..., maxD.

This is a **fixed-width envelope** approach. It does NOT adapt based on what
peaks are actually observed. If the D2O fraction is 0.833, a 10-residue peptide
(no prolines) would look for (10-2) * 1 * 0.833 + 1 = 7.66 -> 7+1 = 8 peaks
(M+0 through M+7).

### 2.4 How Individual Peaks Are Matched: `compare()`

This is a **binary search** on the sorted peak list:

```python
def compare(target, charge, array, full_array):
    midpoint = int(len(array) / 2)
    try:
        if abs(get_ppm(target, array[midpoint][0] * charge)) <= CON.PPM_MATCH_TOLERANCE:
            return_list = [(array[midpoint][0], array[midpoint][1])]
            offset = 1
            # Check adjacent peaks in both directions
            while offset != 0 and (midpoint - offset) > 0:
                peak = array[midpoint - offset]
                ppm = abs(get_ppm(target, peak[0] * charge))
                if ppm <= CON.PPM_MATCH_TOLERANCE:
                    return_list.append((peak[0], peak[1]))
                    offset += 1
                else:
                    offset = 0
            offset = 1
            while offset != 0 and (midpoint + offset < len(array)):
                peak = array[midpoint + offset]
                ppm = abs(get_ppm(target, peak[0] * charge))
                if ppm <= CON.PPM_MATCH_TOLERANCE:
                    return_list.append((peak[0], peak[1]))
                    offset += 1
                else:
                    offset = 0
            # Pick highest-intensity match
            high_intensity = (0, 0)
            for key, value in return_list:
                if value > high_intensity[1]:
                    high_intensity = key, value
            ppm_error = abs(get_ppm(target, high_intensity[0] * charge))
            return ppm_error, high_intensity[0], high_intensity[1]
        elif len(array) == 1 or len(array) == 0:
            return 0, 0, 0
        elif array[midpoint][0] * charge <= target:
            return compare(target, charge, array[midpoint:], full_array)
        else:
            return compare(target, charge, array[0: midpoint], full_array)
    except IndexError:
        return 0, 0, 0
```

Key behavior:
- **Binary search** recursively narrows the sorted peak array
- When a match is found within PPM tolerance, expands outward to find all
  adjacent matches
- **Picks the highest-intensity peak** among all matches for a given deuteration
  slot
- Returns (ppm_error, m/z, intensity) or (0, 0, 0) if no match

### 2.5 Calling Sequence in `match_peptides()`

```python
def match_peptides(self, pep):
    start, end = pep.get_rt_start_end()
    charge = pep.get_charge()
    pep_mass = pep.get_mass_over_charge() * charge
    tuple_list = []
    for rt, tup_list in self.get_tuple_dictionary().items():
        if start <= rt[0] <= end or start <= rt[1] <= end:
            tuple_list.extend(tup_list)
    tuple_list.sort(key=lambda x: x[0])
    for det in range(pep.get_max_deuterium() + 1):
        ppm_error, mz, intensity = compare(
            pep_mass + det * CON.DEUTERIUM_MASS_DIFFERENCE,
            charge, tuple_list, tuple_list)
        if ppm_error != 0:
            pep.set_deuterium(det, mz, intensity, ppm_error)
    pep.set_fit()
    pep.set_weighted_mass()
```

The flow:
1. Collect all peaks from sliding windows that overlap the peptide's RT window
2. Sort by m/z
3. For each deuteration level (0 to maxD), compute the theoretical mass and
   binary-search for a matching peak
4. After all peaks matched, compute Gaussian fit and weighted mass

### 2.6 Gaussian Quality Assessment

```python
def fit_gaussian(y_data):
    x_data = list(range(len(y_data)))
    consecutive_non_zeroes = 0
    for y in y_data:
        if consecutive_non_zeroes == 3:
            break
        if y > 0:
            consecutive_non_zeroes += 1
        else:
            consecutive_non_zeroes = 0
    if consecutive_non_zeroes != 3:
        return "Insufficient Match"
    try:
        popt, pcov = curve_fit(gaussian_trend, x_data, y_data)
        r_squared = calculate_r_squared(x_data, y_data, *popt)
        if r_squared < 0:
            return 0
        return r_squared
    except RuntimeError:
        return "Couldn't Fit"
```

**Important details:**
- Requires at least 3 consecutive non-zero intensities or returns
  "Insufficient Match"
- Fits `a * exp(-(x - mean)^2 / (2*sigma^2))` to the intensity profile
- Reports R-squared; the paper recommends R^2 > 0.9 as the quality threshold
  (91% of such spectra were manually verified as acceptable)
- The Gaussian fit is a **quality metric only** -- it does NOT affect the
  centroid calculation itself

---

## 3. RT Window Handling

### 3.1 Retention Time Mapping

The `set_retention_times()` function reads the identification mzML and builds
a scan-number-to-RT dictionary from MS2 scans:

```python
def set_retention_times(file: str):
    retention_scan_dictionary = {}
    with mzml.read(file) as f:
        for scan in f:
            if scan["ms level"] == 2:
                scan_time = float(scan["scanList"]["scan"][0]["scan start time"])
                scan_time = (scan_time - CON.RETENTION_SHIFT_INTERCEPT) / CON.RETENTION_SHIFT_SLOPE
                scan_time *= CON.MINUTES_TO_SECONDS
                retention_scan_dictionary[scan["index"] + 1] = scan_time
    return retention_scan_dictionary
```

This maps each MS2 scan number to its retention time (in seconds), with an
optional linear calibration (slope/intercept) for RT alignment between the
identification run and experimental runs.

### 3.2 RT Window for Each Peptide

```python
def set_pep_retention_times(self, file: str):
    conversion_dictionary = set_retention_times(file)
    for pep in self.peptides:
        scan = pep.get_scan()
        rt = conversion_dictionary[scan]
        pep.set_retention_times(rt - CON.RETENTION_TOLERANCE,
                                rt + CON.RETENTION_TOLERANCE)
```

Each peptide gets a window of `[RT - 30s, RT + 30s]` (default). This is a
**fixed symmetric window** centered on the identification RT.

### 3.3 Sliding Window for Peak Aggregation

The sliding window merges peaks from multiple scans:

```python
def sliding_window(self, retention_time: float):
    window_count = int(int((retention_time // CON.SLIDING_WINDOW_SIZE) + 1) *
                      (CON.SLIDING_WINDOW_SIZE / SLIDE_AMOUNT)) - 1
    start = 0
    stop = CON.SLIDING_WINDOW_SIZE
    windows = []
    for i in range(window_count):
        windows.append((start, stop))
        start += SLIDE_AMOUNT
        stop += SLIDE_AMOUNT
    window_dictionary = {}
    for window in windows:
        key = window
        window_dictionary[key] = []
    for dictionary in self.all_peaks:
        if window[0] <= dictionary["retention time"] < window[1]:
            window_dictionary[key].extend(dictionary["tuple list"])
    # ... tuple_combine called on each window ...
```

**Note: There appears to be a bug here.** The loop `for dictionary in
self.all_peaks` only checks against the LAST `window` from the preceding loop
(the variable `window` is not re-bound inside this loop). This means only the
last time window gets peaks assigned to it. The correct logic would need to
iterate over all windows for each scan. **This is potentially a significant bug
in the V1.3 code.**

After collecting peaks per window, `tuple_combine()` merges peaks within 1 PPM
of each other (summing intensities, averaging m/z weighted by count):

```python
def tuple_combine(some_list):
    # Sort by m/z, then iteratively merge adjacent peaks within PPM tolerance
    # Weighted average for m/z, sum for intensity
    # Repeats until no more merges are possible
```

### 3.4 Integration vs Single Scan

The approach **integrates across all scans in the RT window**. It does NOT pick
a single scan. All peaks from all MS1 scans within the sliding windows that
overlap the peptide's RT window are pooled, merged by `tuple_combine()`, then
searched for matches. This is essentially a **summed spectrum** approach.

There is no XIC (extracted ion chromatogram) logic. The tool does not track the
elution profile of individual peptides. It simply pools everything within the
RT window.

---

## 4. Deuterium Uptake Calculation

### 4.1 Mass Shift Formula

The deuterium uptake for each peptide is the mass shift:

```
mass_shift = (C_mz - m_proton) * z - M_average
```

where:
- `C_mz` = intensity-weighted average m/z from matched isotope envelope
- `m_proton` = 1.007276 Da (MASS_OF_HYDROGEN)
- `z` = charge state
- `M_average` = sum of average amino acid residue masses + 18.01528 (water)

### 4.2 Non-Deuterated Control Subtraction

The program's `main()` function shows a critical step: **non-deuterated baseline
subtraction**. Before processing experimental data:

1. A non-deuterated mzML file is processed first (time = -10 s sentinel value)
2. The mass shift is computed for the undeuterated sample
3. For each experimental file, the undeuterated shift is subtracted:

```python
for j in range(Peptide_count):
    diff = df1.loc[j, "Shift"] - df2.loc[j, "Shift"]
    df1.loc[j, "Shift"] = diff
```

This means the final uptake = (deuterated centroid shift) - (undeuterated
centroid shift). This is the standard approach.

### 4.3 Back-Exchange Correction

**Minimal.** The code has two parameters:
- `DEUTERIUM_RECOVERY_RATE` (default 1.0 = no correction)
- `DEUTERIUM_FRACTION` (default 0.833)

These are used ONLY in `sequence_to_max_deuterium()` to limit the number of
envelope peaks searched. They do NOT correct the actual mass shift value.

The paper explicitly states they rely on experimental minimization of
back-exchange ("sub-zero temperature UPLC-HDX-MS platform") rather than
computational correction.

**There is no fully deuterated control processing.** No correction formula like:
```
D_corrected = D_measured * maxD_theoretical / (D_fully_deut - D_undeut)
```

### 4.4 Differential Uptake

For differential experiments (two conditions), the uptake difference is:
```
delta_D = uptake_condition2 - uptake_condition1
```

Statistical significance is assessed via confidence intervals using a two-tailed
Student's t-test with pooled standard deviation from replicates.

### 4.5 Fractional Uptake

```python
def max_deuterium_mlf(sequence: str):
    max_deuterium = len(sequence) - 2
    for letter in sequence:
        if letter.lower() == 'p':
            max_deuterium -= 1
    max_d_mlf = max_deuterium * CON.DEUTERIUM_FRACTION
    return max_d_mlf
```

Fractional uptake = mass_shift / max_deuterium_mlf, representing the fraction
of exchangeable amides that incorporated deuterium.

---

## 5. Output Format

### 5.1 Detailed Output (per experimental run)

CSV file with one row per (peptide, deuteration level) pair:

| Column | Description |
|--------|-------------|
| `Start` | Residue start position in protein |
| `End` | Residue end position |
| `Sequence` | Peptide sequence |
| `Charge` | Charge state |
| `SequenceMz` | Theoretical monoisotopic m/z |
| `Condition` | "Free" or "Complex" (condition label) |
| `Deuterium` | Deuteration level (0, 1, 2, ...) |
| `RT` | Retention time in minutes |
| `Mz` | Matched experimental m/z |
| `Intensity` | Matched peak intensity |
| `PpmError` | PPM error of match |
| `Average` | Weighted average m/z (centroid) |
| `Shift` | Mass shift (deuterium uptake in Da) |
| `Gaussian Fit` | R-squared of Gaussian fit |

### 5.2 Recommendation Table 1 (IC-HDX format)

Peptide-level summary with uptake at each timepoint:

| Column | Description |
|--------|-------------|
| `Start` | Residue start |
| `End` | Residue end |
| `Sequence` | Peptide sequence |
| `Peptide mass (Da)` | Average mass |
| `Retention time (min)` | Average RT |
| `Uptake [Condition] (D)` | Average uptake per timepoint (multiple columns) |
| `Uptake error (SD) [Condition] (D)` | Std dev per timepoint |

### 5.3 Recommendation Table 2 (long format, IC-HDX)

| Column | Description |
|--------|-------------|
| `Protein state` | Condition label |
| `Sequence` | Peptide sequence |
| `Start` | Residue start |
| `End` | Residue end |
| `Peptide mass (Da)` | Monoisotopic mass |
| `Retention time (min)` | Average RT |
| `HDX time (s)` | Timepoint |
| `Uptake (D)` | Average deuterium uptake |
| `Uptake SD (D)` | Standard deviation |

### 5.4 Compatibility with PyHDX

**Not directly compatible.** PyHDX expects Waters DynamX "state data" format
with different column names and conventions. The Recommendation Table 2 output
is close to what PyHDX needs but would require column renaming/reformatting:

- PyHDX wants: `start`, `end`, `sequence`, `state`, `exposure`, `uptake`, `uptake_sd`
- TheDeuteriumCalculator produces: `Start`, `End`, `Sequence`, `Protein state`,
  `HDX time (s)`, `Uptake (D)`, `Uptake SD (D)`

A simple column rename + case change would likely suffice. The tool also has a
DECA format conversion option (menu option 4) for downstream analysis.

---

## 6. Validation

### 6.1 Test Data

**No test data is included in the repository.** The repo contains only:
- `TheDeuteriumCalculator_V1.3.py`
- `PARAMETERS.py`
- `README.md`

The paper references supplementary files (antigen mzML, peptide library CSV,
antigen sequence file) deposited on Figshare, but these are not in the repo.

### 6.2 Published Validation Results

From the paper:
- **Gaussian quality**: R^2 > 0.9 threshold correctly identified 91% of
  manually validated acceptable spectra
- **Processing speed**: ~1 minute for a 300 MB mzML file on desktop hardware
- **Test system**: Anthrax vaccine antigen (rPA) and E. coli lysate digest
- **Reproducibility**: "high reproducibility observed for matched peptides
  among runs"

### 6.3 Comparison with Commercial Software

**No direct comparison with DynamX, HDExaminer, or other commercial tools is
presented.** The paper validates by manual spectral inspection, not by
cross-software comparison. Output can be exported to DECA for downstream
statistical analysis.

---

## 7. What We Can Learn and Adopt

### 7.1 Algorithms Worth Adopting

1. **Intensity-weighted centroid calculation** -- Their `set_weighted_mass()`
   function is the standard approach and matches our implementation in
   `HDX_RECONSTRUCTION_METHODS.md`. The formula is correct:
   ```
   C_mz = sum(I * mz) / sum(I)
   mass_shift = (C_mz - m_proton) * z - M_average
   ```

2. **Gaussian fit as quality metric** -- Fitting a Gaussian to the isotope
   envelope intensities and using R^2 as a quality score is a good, simple
   approach. R^2 > 0.9 is their recommended threshold. We should implement
   this as a data quality filter.

3. **Non-deuterated baseline subtraction** -- Critical step: always subtract
   the undeuterated centroid shift from the deuterated shift. This corrects
   for any systematic m/z bias and for the natural isotope distribution.

4. **Fixed envelope width from sequence** -- Computing maxD from sequence
   length minus 2 (N-term + first amide) minus prolines is the standard
   approach. Their formula:
   ```
   maxD = (len - 2 - n_proline) * recovery_rate * D2O_fraction + 1
   ```

5. **Binary search for peak matching** -- Efficient O(log n) search on sorted
   peak lists. Good for large datasets.

### 7.2 Pitfalls They Handle That We Should Check

1. **Minimum 3 consecutive peaks** -- They require at least 3 consecutive
   non-zero intensities before accepting an envelope. This prevents false
   matches from noise.

2. **Sentinel values for failed matches** -- Setting mass_shift = -1 when
   fitting fails or no peaks match. Downstream code must check for these.

3. **PPM error tracking per peak** -- They store and report the PPM error
   for each matched peak, enabling post-hoc quality assessment.

4. **Proline exclusion** -- Prolines have no backbone amide NH and cannot
   exchange; they correctly subtract these from maxD.

### 7.3 Pitfalls and Limitations in Their Code

1. **Possible sliding window bug** -- The loop that assigns peaks to windows
   appears to only populate the last window. If this is indeed a bug, it means
   the code may work only because `match_peptides()` collects peaks from ALL
   overlapping windows, and the window dictionary gets populated correctly
   through `tuple_combine()`. This needs investigation.

2. **No XIC / no elution profile** -- They pool all peaks in the RT window
   without tracking the peptide's elution profile. This means co-eluting
   interference is not detected.

3. **No back-exchange correction** -- Relying on experimental minimization
   only. For quantitative HDX, a fully deuterated control and correction
   formula would be better.

4. **Average masses in amino acid dictionary** -- Using average rather than
   monoisotopic masses for the reference mass. This is fine for computing
   uptake (since they subtract undeuterated control), but could cause issues
   if the reference mass is used for other purposes.

5. **No overlapping charge state handling** -- Each peptide is processed
   independently at its identified charge state. No merging of information
   across charge states.

6. **Hardcoded noise threshold** -- A flat intensity cutoff (10000) is too
   simplistic for instruments with varying signal levels across the m/z range.

7. **Interactive prompts** -- The script requires manual user input (file
   paths, parameters) at runtime, making it unsuitable for automated
   pipelines.

### 7.4 License Compatibility

**No license file exists in the repository.** The paper calls it "open-source"
and it is published in a public GitHub repository, but without an explicit
license, default copyright applies (all rights reserved). Strictly speaking,
we cannot copy code from this repository.

**However:** The algorithms are described in a published paper and represent
standard HDX-MS methodology (weighted centroid, Gaussian fitting). We can
implement the same algorithms independently based on the published equations
without any license concern. The specific code snippets in this document are
for reference and understanding only -- our implementation should be written
from scratch using the same mathematical formulas.

---

## 8. Summary: Key Takeaways for Our Reconstruction Script

### What Their Centroid Calculation Does (and We Should Match)

1. For each peptide, compute maxD from sequence (length - 2 - prolines)
2. For each deuteration level n = 0 to maxD:
   - Compute theoretical mass: `M_theo = M_mono + n * 1.00628`
   - Search for matching peak within PPM tolerance
   - Pick highest-intensity match if multiple candidates
3. Compute centroid: `C_mz = sum(I_n * mz_n) / sum(I_n)` over ALL matched peaks
4. Convert to uptake: `uptake = (C_mz - m_proton) * z - M_average`
5. Subtract non-deuterated baseline uptake
6. Validate with Gaussian fit R^2 > 0.9

### What They Do for RT/Scan Integration

- Pool ALL MS1 scans within +/- 30 seconds of identification RT
- Merge peaks from all scans using a sliding window approach
- Within each window, merge peaks within 1 PPM (weighted average m/z, summed
  intensity)
- This effectively creates a "summed spectrum" across the elution window

### What's Missing (and We May Need)

- Back-exchange correction using fully deuterated control
- XIC-based elution profile analysis
- Adaptive envelope width based on observed peaks
- Multi-charge-state merging
- Automated peak quality filtering beyond R^2 threshold
- Non-interactive batch processing

---

## References

- Cupp-Sutton KA, Welborn T, Fang M, Langford JB, Wang Z, Smith K, Wu S.
  "The Deuterium Calculator: An Open-Source Tool for Hydrogen-Deuterium Exchange
  Mass Spectrometry Analysis." *J Proteome Res.* 2023;22(3):817-825.
  DOI: 10.1021/acs.jproteome.2c00558. PMID: 36695755. PMC: PMC12510632.

- GitHub: https://github.com/OUWuLab/TheDeuteriumCalculator

- Zhang Z, Smith DL. "Determination of amide hydrogen exchange by mass
  spectrometry: A new tool for protein structure elucidation." *Protein Sci.*
  1993;2:522-531.

- Masson GR, et al. "Recommendations for performing, interpreting and reporting
  hydrogen deuterium exchange mass spectrometry (HDX-MS) experiments."
  *Nat Methods.* 2019;16:595-602.
