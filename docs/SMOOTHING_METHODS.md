# Smoothing Methods for 2D IM-MS Interactive Viewer

## Research Findings

### UniDec (Universal Deconvolution)
UniDec's primary smoothing method for IM-MS data is **Gaussian smoothing** applied
as a separable 2D filter. The key insight from UniDec's source code (`unidectools.py`)
is that Gaussian kernels have all-positive weights, making it mathematically impossible
to produce negative values from non-negative input data. This is critical for sparse
IM-MS heatmaps where many bins are zero.

### MassLynx (Waters)
MassLynx uses Savitzky-Golay (SG) for 1D chromatographic traces where data is dense
and continuous. For 2D IM-MS displays, it applies mean/median filters or Gaussian
smoothing to fill sparse regions.

### Literature
- Marty et al., Anal. Chem. 2015: UniDec uses Gaussian smoothing for charge state
  deconvolution of native MS data.
- Ruotolo et al., Science 2005: CCS measurements from IM-MS use Gaussian fitting
  for arrival time distributions.

## Method Comparison

| Method          | Negatives? | Peak shape    | Fill sparse gaps | Speed   |
|-----------------|------------|---------------|------------------|---------|
| Gaussian        | No         | Broadens      | Yes              | Fast    |
| Savitzky-Golay  | **Yes**    | Preserves     | No               | Fast    |
| Moving Average  | No         | Broadens      | Yes              | Fastest |

## Why Gaussian is Preferred for 2D Heatmaps

1. **No negative values**: All kernel weights are positive (exp(-x^2) > 0), so
   the convolution of non-negative data is always non-negative. SG polynomial
   fitting can overshoot into negatives, especially near sharp edges in sparse data.

2. **Fills sparse gaps**: IM-MS heatmaps are inherently sparse (many zero bins).
   Gaussian spreading naturally fills neighboring bins, producing a visually
   continuous distribution. SG preserves zeros between peaks.

3. **Matches physical model**: Ion arrival time distributions are approximately
   Gaussian, so Gaussian smoothing is physically motivated. The smoothing kernel
   shape matches the expected peak shape.

## SG Negative Value Explanation

Savitzky-Golay fits a polynomial (degree k) to each sliding window of width w.
When data has sharp transitions (e.g., a peak next to zeros in sparse IM-MS data),
the polynomial can overshoot below zero on the trailing edge. This is called
**Gibbs-like ringing** and is more severe with:
- Wider windows (large w)
- Higher polynomial order (large k)
- Sparser data (more zero bins)

The viewer clamps SG output to zero (`Math.max(0, s)`) to prevent negative display
values, but this creates flat-bottom artifacts. Gaussian smoothing avoids the
problem entirely.

## Preset Rationale

### 2D Smooth Presets
- **Raw**: No smoothing — see the actual binned data
- **Gaussian sigma=1**: Light smoothing, minimal broadening
- **Gaussian sigma=2**: Standard smoothing for typical IM-MS data (recommended)
- **Gaussian sigma=3**: Heavy smoothing for very sparse data
- **SG(7,2)**: Light polynomial smooth (preserved for comparison)
- **SG(15,3)**: Heavy polynomial smooth
- **Mov Avg 5**: Simple 5-point average
- **Mov Avg 11**: Wider average for noisy data

### Drift Smooth Presets
- **Raw**: No smoothing — see the actual drift profile
- **Gaussian sigma=1**: Light smoothing
- **Gaussian sigma=2**: Moderate smoothing
- **SG(7,2)**: Light polynomial (good for dense drift profiles)
- **SG(15,3)**: Heavy polynomial
- **Mov Avg 5**: Simple 5-point average

### Noise Threshold
Removes pixels below X% of the maximum intensity. Useful for cleaning up
background noise after smoothing.

## Implementation Details

- **Separable 2D**: All 2D smoothing methods apply 1D along columns (drift),
  then 1D along rows (m/z). This is mathematically equivalent to full 2D
  convolution for separable kernels (Gaussian, box/MA) and is O(n*w) vs O(n*w^2).
- **Reflected boundaries**: All methods use `reflect` mode at array edges,
  matching scipy's default (`mode='reflect'`).
- **JSON presets**: Each dropdown option stores `{"method":"...", ...params}`
  as its value. `JSON.parse()` extracts the preset, and dispatch functions
  route to the correct smoothing implementation.
