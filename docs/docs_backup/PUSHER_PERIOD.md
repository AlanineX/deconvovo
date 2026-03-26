# Pusher Period: Drift Time Calibration for Waters SYNAPT IM-MS

## What is the Pusher Period?

In Waters SYNAPT time-of-flight (TOF) mass spectrometers, the **pusher** is an
electrode that periodically pulses ions from the ion source into the TOF flight
tube. The time between consecutive pulses is the **pusher period** (also called
"pusher interval"). It is measured in microseconds (μs).

In ion mobility (IM) mode, the instrument records which pusher cycle each ion
arrives at — this is the **drift bin number** (an integer, typically 0–199).
Each drift bin corresponds to exactly one pusher cycle (when
`ADC Pushes Per IMS Increment = 1`, the typical setting).

The conversion from drift bins to physical drift time is:

```
drift_time_ms = bin_number × pusher_period_μs / 1000
```

For example, with a 69 μs pusher period:
- Bin 0 → 0.00 ms
- Bin 100 → 6.90 ms
- Bin 199 → 13.73 ms

## Why the Pusher Period Matters

The pusher period determines the **time resolution** of the IM separation axis.
It directly affects:
- Drift time values displayed on plots
- CCS (collisional cross-section) calculations derived from drift time
- Comparison of drift profiles between runs with different mass ranges

Runs with different m/z acquisition ranges will have **different pusher
periods** because heavier ions need longer flight times in the TOF tube. For
example:
- 50–2000 Da range → 69 μs pusher period
- 250–3000 Da range → 110 μs pusher period

Using the wrong pusher period gives incorrect drift times. A 32% error was
found when using a TOF physics approximation instead of the actual instrument
value (see "How We Found It" below).

## How the Pusher Period is Stored

The pusher period is stored in the Waters `.raw` directory as **stat code 76**
("Pusher Frequency") inside the `_FUNC001.STS` binary file. This is a per-scan
status/statistics file that records instrument parameters for each scan.

### STS File Format

```
Offset  Content
------  -------
0x00    8-byte header: [unknown(2), unknown(2), n_descriptors(2), n_scans(2)]
0x08    Padding (24 bytes)
0x20    Descriptor block: n_desc × 32 bytes each
        Each descriptor: [stat_code(2), type(2), data_offset(2), name(26)]
...     16-byte data section header
...     Per-scan records: n_desc bytes each (variable-width packed values)
```

Key descriptor for stat code 76:
```
code=76, type=1, data_offset=40, name="Pusher Frequency"
```

The value at `data_offset=40` in each scan record is the pusher frequency in
**Hz** (integer). The pusher period in μs is:

```
pusher_period_μs = floor(1,000,000 / pusher_frequency_Hz)
```

### Values Found in Our Data

| Dataset | Acq Range (Da) | Pusher Freq (Hz) | Pusher Period (μs) |
|---------|----------------|-------------------|--------------------|
| 20260212 ADP runs | 50–2000 | 14,440 | 69 |
| 20260212 UBQ calibrant | 50–2000 | 14,440 | 69 |
| 20260216 protein runs | 250–3000 | 9,070 | 110 |
| 20251105 MYO trypsin | 250–5000 | 9,070 | 110 |

## How Other Software Reads It

### UniDec (University of Arizona)

UniDec reads the pusher frequency via the MassLynx SDK (Windows DLL):

```python
# From unidec/modules/waters_import_wizard/data_importer.py
pusher = get_stat_name(file_path, "Transport RF")  # legacy name for stat code 76
if pusher is not None:
    pusher_interval_us = floor((1.0 / pusher) * 1000 * 1000)
```

The `get_stat_name()` function calls `MassLynxRawInfoReader.GetScanItems()` →
`MassLynxRawReader.massLynxDll.getScanItems()` — a Windows-only DLL. UniDec
then applies the conversion:

```python
# From unidec/modules/IM_functions.py
if config.pusher != 0:
    ygrid = np.array(ygrid) * config.pusher * 0.001  # bin → ms
```

UniDec's GUI tooltip documents: `AT = Pusher * Bin / 1000` where AT = arrival
time in ms, Pusher = interval in μs, Bin = drift bin number.

Note: UniDec labels this stat as **"Transport RF"** for historical reasons —
it is NOT the transport RF voltage. The name refers to the pusher repetition
frequency, which in early Waters firmware was derived from the transport RF
oscillator.

### ms_deisotope (J. Klein, Boston University)

ms_deisotope uses the MassLynx SDK's built-in conversion function:

```python
# From ms_deisotope/data_source/_vendor/masslynx/loader.py
drift_time = self.info_reader.GetDriftTime(function, scan)
```

`GetDriftTime()` calls `massLynxDll.getDriftTime()` — an opaque DLL function
that returns the drift time in ms directly. This is the most authoritative
source but requires Windows.

ms_deisotope also documents the full formula in a comment:

```
Cyclic arrival time = (bin number * pusher period * pushesperbin) + ADC start delay
```

### MassLynx (Waters, proprietary)

MassLynx displays drift time in ms on its IM-MS plots. It reads the pusher
period from the same internal source (stat code 76 via the SDK). The SDK
provides `GetDriftTime()` which handles the conversion including the ADC start
delay offset.

## How We Read It (Our Pipeline)

Since the MassLynx SDK DLL is Windows-only, we parse the `_FUNC001.STS` binary
directly on any platform:

```python
def _read_pusher_from_sts(sts_path):
    data = sts_path.read_bytes()
    n_desc = struct.unpack_from('<H', data, 4)[0]
    data_start = 0x20 + n_desc * 32
    ds = data[data_start:]

    # Find stat code 76 descriptor
    for i in range(n_desc):
        off = 0x20 + i * 32
        code = struct.unpack_from('<H', data, off)[0]
        if code == 76:
            data_offset = struct.unpack_from('<H', data, off + 4)[0]
            break

    # Read from first scan record (16-byte header in data section)
    header_size = len(ds) % n_desc
    freq_hz = struct.unpack_from('<i', ds, header_size + data_offset)[0]
    pusher_us = math.floor(1e6 / freq_hz)
    return pusher_us
```

This gives the same value UniDec would read via the Windows SDK, verified by
comparing against all runs in our datasets.

## What We Tried Before (and Why It Failed)

### TOF Physics Calculation (abandoned)

Initially we calculated the pusher period from the TOF flight time equation:

```
t_max = Lteff × sqrt(mz_max × u / (2 × Veff × e))
```

Where `Lteff`, `Veff` are from `_extern.inf` and `mz_max` is the acquisition
end mass. This gives the **maximum ion flight time** in the TOF tube, which
should approximate the pusher period.

**Results:**

| Dataset | TOF Calc (μs) | Actual (μs) | Error |
|---------|---------------|-------------|-------|
| 50–2000 Da | 68.17 | 69 | -1.2% |
| 250–3000 Da | 83.49 | 110 | **-31.8%** |
| 250–5000 Da | 107.78 | 110 | -2.0% |

The 250–3000 Da runs had a **32% error** because the actual pusher period
includes dead time, ADC overhead, and is often set conservatively by the
instrument firmware — it is NOT simply the maximum flight time. The 250-5000
runs happened to be close by coincidence.

### Fields NOT Found in _extern.inf

The following fields do NOT exist in the SYNAPT G2 `_extern.inf`:
- "Pusher Interval" — UniDec's code has a commented-out search for this
- "Pusher Frequency" — not in _extern.inf, only in STS binary
- "Pusher Period" — not stored as a named field anywhere

The `Pusher = 1900.0` field in `_extern.inf` is the pusher electrode **voltage**
in volts, not a timing parameter.

## Limitations

1. **ADC Start Delay**: ms_deisotope documents a formula that includes an
   "ADC start delay" offset: `arrival_time = (bin × pusher × pushesperbin) + delay`.
   Our conversion omits this offset (as does UniDec), meaning bin 0 maps to
   0.0 ms instead of the true start delay. For relative drift time comparisons
   this is fine; for absolute CCS calculations, the offset matters.

2. **Multiple Pushes Per Bin**: If `ADC Pushes Per IMS Increment > 1`, each
   drift bin spans multiple pusher cycles. The conversion becomes:
   `drift_time_ms = bin × pusher_us × pushes_per_bin / 1000`.
   All our data has pushes_per_bin = 1.

3. **Per-Scan Variation**: The STS file stores the pusher frequency per scan.
   We read only the first scan. In practice, the pusher frequency is constant
   across all scans within a run, but it could theoretically vary.

## Summary

| Approach | Source | Platform | Accuracy |
|----------|--------|----------|----------|
| MassLynx SDK `GetDriftTime()` | DLL | Windows only | Gold standard |
| UniDec: stat "Transport RF" via SDK | DLL | Windows only | Good (no ADC delay) |
| **Our pipeline: stat code 76 from STS binary** | Binary parse | **Cross-platform** | **Same as UniDec** |
| TOF physics from Lteff/Veff/EndMass | _extern.inf | Cross-platform | **Unreliable (up to 32% error)** |
