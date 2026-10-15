"""Native IM-MS display arrays from immutable converted Waters text data.

The desktop renderer consumes these arrays directly. Portable HTML generation
remains separate; it is never loaded or interpreted by this module.
"""
from __future__ import annotations

from dataclasses import dataclass
from copy import deepcopy
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.ndimage import correlate1d
from scipy.signal import find_peaks, savgol_coeffs

from deconvovo.paths import local_path

MAX_HEATMAP_CELLS = 8_000_000


@dataclass(frozen=True)
class ViewSnapshot:
    mz_centres: np.ndarray
    drift_centres: np.ndarray
    display2d: np.ndarray
    drift_intensity: np.ndarray
    ms_mz: np.ndarray
    ms_intensity: np.ndarray
    mz_range: tuple[float, float]
    drift_range: tuple[float, float]
    drift_unit: str
    intensity_label: str
    transform: str
    raw2d: np.ndarray
    visible_mz_slice: slice
    visible_drift_slice: slice
    baked_mz_range: tuple[float, float]
    ms_mode: str
    drift_profile_centres: np.ndarray | None = None


def _readonly(array):
    result = np.asarray(array, dtype=np.float64)
    result.setflags(write=False)
    return result


def _bounds(value, fallback, name):
    if value is None:
        return tuple(fallback)
    if len(value) != 2 or not all(math.isfinite(float(v)) for v in value) or value[0] >= value[1]:
        raise ValueError(f"The {name} minimum must be finite and lower than its maximum.")
    return tuple(map(float, value))


def _configuration(config):
    if config is not None:
        return deepcopy(config)
    base = Path(sys._MEIPASS) if getattr(sys, "frozen", False) else Path(__file__).parent.parent
    path = base / "config" / "imms_plot_config.json"
    if path.is_file():
        return json.loads(path.read_text(encoding="utf-8"))
    return {"defaults": {"mz_bins": 800, "smooth_2d": {"method": "raw"},
        "smooth_drift": {"method": "raw"}, "scale": "1", "z_scale": "log10"}}


def _native_bin_position(value, width):
    position = value / width
    nearest = round(position)
    # Decimal axis entries such as 0.3 / 0.1 can land just below bin 3.
    tolerance = 2 * sys.float_info.epsilon * max(1, abs(position))
    return nearest if abs(position - nearest) <= tolerance else position


def _filter_axis(values, kernel, axis, boundary="mirror"):
    """Use SciPy convolution, retaining the HTML's single-fold edge convention."""
    radius = len(kernel) // 2
    n = values.shape[axis]
    if boundary == "nearest" or radius < n:
        return correlate1d(values, kernel, axis=axis, mode=boundary)
    # The browser folds an out-of-range index once, then clamps it. SciPy's
    # repeated reflection differs only when a custom kernel exceeds the axis.
    indices = np.arange(-radius, n + radius)
    indices = np.where(indices < 0, -indices, indices)
    indices = np.where(indices >= n, 2 * n - 2 - indices, indices)
    padded = np.take(values, np.clip(indices, 0, n - 1), axis=axis)
    filtered = correlate1d(padded, kernel, axis=axis, mode="nearest")
    selection = [slice(None)] * values.ndim
    selection[axis] = slice(radius, radius + n)
    return filtered[tuple(selection)]


def _smooth_axis(values, preset, axis=0):
    method = preset.get("method", "raw")
    n = values.shape[axis]
    if method == "gaussian":
        sigma = float(preset.get("sigma", 0))
        if not math.isfinite(sigma):
            raise ValueError("The Gaussian smoothing width must be finite.")
        if sigma <= 0 or n < 3:
            return values.copy()
        radius = math.ceil(4 * sigma)
        if radius >= n and (n + 2 * radius) * values.size // n > MAX_HEATMAP_CELLS:
            raise ValueError("The Gaussian smoothing window is too large. Choose less smoothing or fewer bins.")
        offsets = np.arange(-radius, radius + 1, dtype=float)
        kernel = np.exp(-0.5 * (offsets / sigma) ** 2)
        kernel /= kernel.sum()
        return np.maximum(0, _filter_axis(values, kernel, axis))
    if method == "sg":
        window, order = int(preset.get("w", 0)), int(preset.get("k", 2))
        if window < 3 or window > n or window % 2 == 0:
            return values.copy()
        try:
            kernel = savgol_coeffs(window, order)
        except ValueError:
            return values.copy()  # singular polynomial fit leaves the input unchanged
        return np.maximum(0, _filter_axis(values, kernel, axis, "nearest"))
    if method == "ma":
        window = int(preset.get("w", 0))
        if window < 2 or n < window:
            return values.copy()
        size = 2 * (window // 2) + 1
        return np.maximum(0, _filter_axis(values, np.full(size, 1 / size), axis))
    return values.copy()


def _odd(value):
    rounded = math.floor(value + 0.5)  # JavaScript Math.round for nonnegative widths
    return 3 if rounded < 3 else rounded + (rounded % 2 == 0)


def _smooth_2d(raw, mz_domain, drift_centres, pusher_us, defaults):
    preset = defaults.get("smooth_2d", {"method": "raw"})
    method = preset.get("method", "raw")
    ratio = raw.shape[1] / raw.shape[0]
    scale = defaults.get("scale", "1")
    factor = ratio if scale == "full" else math.sqrt(ratio) if scale == "sqrt" else float(scale)
    if not math.isfinite(factor) or factor < 0:
        raise ValueError("The m/z smoothing scale must be finite and nonnegative.")
    if method == "gaussian":
        base_sigma = float(preset.get("sigma", 0)) / ((mz_domain[1] - mz_domain[0]) / raw.shape[1])
        drift_width = abs(drift_centres[1] - drift_centres[0]) if len(drift_centres) > 1 else 2 * abs(drift_centres[0])
        drift_sigma = float(preset["sigma_dr"]) / drift_width if preset.get("sigma_dr") is not None else base_sigma
        result = _smooth_axis(raw, {"method": "gaussian", "sigma": drift_sigma}, 0)
        return _smooth_axis(result, {"method": "gaussian", "sigma": base_sigma * factor}, 1)
    if method in {"sg", "ma"}:
        drift_preset = dict(preset)
        mz_preset = dict(preset, w=_odd(float(preset.get("w", 0)) * factor))
        result = _smooth_axis(raw, drift_preset, 0)
        return _smooth_axis(result, mz_preset, 1)
    return raw


def _display_values(raw, mz_domain, drift_centres, pusher_us, defaults):
    values = np.maximum(0, _smooth_2d(raw, mz_domain, drift_centres, pusher_us, defaults))
    transform = defaults.get("z_scale", "log10")
    labels = {"linear": "I", "sqrt": "√I", "cbrt": "∛I", "p4": "I^¼", "p5": "I^⅕",
        "log2": "log₂(I+1)", "log10": "log₁₀(I+1)"}
    if transform == "sqrt":
        values = np.sqrt(values)
    elif transform == "cbrt":
        values = np.cbrt(values)
    elif transform == "p4":
        values = np.power(values, 0.25)
    elif transform == "p5":
        values = np.power(values, 0.2)
    elif transform == "log2":
        values = np.log2(values + 1)
    elif transform != "linear":
        transform = "log10"
        values = np.log10(values + 1)
    cutoff = float(defaults.get("noise_2d", 0))
    if cutoff > 0:
        fade = defaults.get("denoise", "off")
        fade = 0 if fade == "off" else float(fade)
        peak = values.max(initial=0)
        threshold, spread = peak * cutoff / 100, peak * max(0, fade) / 100
        if spread <= 0:
            values = np.where(values < threshold, 0, values)
        else:
            values = values * np.clip((values - (threshold - spread)) / (2 * spread), 0, 1)
    return _readonly(values), transform, labels[transform]


class ViewData:
    """Immutable acquired values; each render is a pure function of the settings."""

    def __init__(self, mz, drift, intensity, ms_mz, ms_intensity, pusher_us, config,
                 profile_mz=None, native_drift_count=None):
        self.raw_mz = _readonly(mz)
        self.raw_drift = np.asarray(drift, dtype=np.int64)
        self.raw_drift.setflags(write=False)
        self.raw_intensity = _readonly(intensity)
        self.raw_ms_mz = _readonly(ms_mz)
        self.raw_ms_intensity = _readonly(ms_intensity)
        self.pusher_us = float(pusher_us or 0)
        if not math.isfinite(self.pusher_us) or self.pusher_us < 0:
            raise ValueError("The pusher period must be finite and nonnegative.")
        self.config = deepcopy(config)
        self.n_drift_native = int(drift.max()) + 1 if native_drift_count is None else int(native_drift_count)
        self.drift_unit = "ms" if self.pusher_us > 0 else "bins"
        self.full_drift_range = (0.0, self.n_drift_native * (self.pusher_us / 1000 if self.pusher_us > 0 else 1))
        self._native_drift_centres = _readonly((np.arange(self.n_drift_native) + 0.5)
            * (self.pusher_us / 1000 if self.pusher_us > 0 else 1))
        self._unique_mz = np.unique(mz)
        # Keep acquired m/z coordinates, including zero-intensity gap samples.
        # The paired MS axis supplies those gaps without imposing a new grid.
        self._profile_mz = _readonly(np.union1d(ms_mz, mz if profile_mz is None else profile_mz))
        self.full_mz_range = (float(self._profile_mz[0]), float(self._profile_mz[-1]))
        self._finest_mz_step = float(np.diff(self._profile_mz).min()) if len(self._profile_mz) > 1 else 0.0
        self._profile_indices = np.searchsorted(self._profile_mz, self.raw_mz)
        self._profile_indices.setflags(write=False)
        self._full_ms_profile = _readonly(np.bincount(self._profile_indices,
            weights=self.raw_intensity, minlength=len(self._profile_mz)))
        deltas = np.diff(self._unique_mz)
        deltas = deltas[(deltas > 0) & (deltas < 0.1)]
        self.native_mz_step = float(np.percentile(deltas, 95) * 1.3) if len(deltas) else 0.01
        self.baked_mz_range = self._domain(config)
        # These caches replace arrays atomically; they never mutate acquired
        # values or make a result depend on the order of earlier interactions.
        self._grid_cache = None
        self._display_cache = None

    def _domain(self, config):
        requested = _bounds((config.get("initial_view") or {}).get("mz_range"), self.full_mz_range, "m/z range")
        clipped = (max(self.full_mz_range[0], requested[0]), min(self.full_mz_range[1], requested[1]))
        # An out-of-domain numeric window is a correctly labelled empty plot.
        domain = clipped if clipped[1] > clipped[0] else requested
        if domain[1] <= domain[0]:
            # A one-m/z acquisition still needs a nonzero cell width.
            domain = (domain[0] - self.native_mz_step / 2, domain[1] + self.native_mz_step / 2)
        defaults = config.get("defaults", {})
        preset = defaults.get("smooth_2d", {})
        sigma = float(preset.get("sigma", 0)) if preset.get("method") == "gaussian" else 0
        if preset.get("method") == "gaussian" and math.isfinite(sigma) and sigma > 0 and self.full_mz_range[1] > self.full_mz_range[0]:
            roi, previous_bins = domain, None
            span = roi[1] - roi[0]
            drift_bins = defaults.get("drift_bins", "native")
            drift_bins = self.n_drift_native if drift_bins == "native" else int(drift_bins)
            while drift_bins > 0:
                mz_bins = defaults.get("mz_bins", 800)
                if mz_bins == "native":
                    mz_bins = min(25600, max(200, int(np.searchsorted(self._unique_mz, domain[1], "right")
                        - np.searchsorted(self._unique_mz, domain[0]))))
                else:
                    mz_bins = int(mz_bins)
                if mz_bins == previous_bins or mz_bins <= 0:
                    break
                previous_bins = mz_bins
                ratio = mz_bins / drift_bins
                scale = defaults.get("scale", "1")
                factor = ratio if scale == "full" else math.sqrt(ratio) if scale == "sqrt" else float(scale)
                if not math.isfinite(factor) or factor < 0:
                    break  # the existing smoothing validation reports the invalid setting
                support = sigma * max(factor, ratio if preset.get("sigma_dr") is None else 0)
                if span >= 8 * support:
                    break
                # Convolve with real surrounding data without widening the
                # requested view or any exact RAW 1D selection. Keep context
                # cells no finer than the finest acquired mass interval.
                padding = max(4 * support, (mz_bins * self._finest_mz_step - span) / 2)
                width = min(self.full_mz_range[1] - self.full_mz_range[0], span + 2 * padding)
                low = max(self.full_mz_range[0], min(roi[0] - padding, self.full_mz_range[1] - width))
                domain = (low, low + width)
        return domain

    def _drift_axis(self, bins):
        step = self.n_drift_native / bins
        edges = np.floor(np.arange(bins + 1) * step).astype(np.int64)
        edges[-1] = self.n_drift_native  # avoid losing the final bin to floating-point floor
        centres = (np.arange(bins) + 0.5) * step * (self.pusher_us / 1000 if self.pusher_us > 0 else 1)
        return edges, _readonly(centres)

    def extract_drift_profile(self, mz_range, config=None):
        """Inclusive RAW m/z ROI on every native acquired drift bin."""
        cfg = deepcopy(self.config if config is None else config)
        defaults = cfg.get("defaults", {})
        lo, hi = _bounds(mz_range, self.full_mz_range, "m/z region")
        return self._native_drift_centres, self._drift_profile(lo, hi, defaults), self.drift_unit

    def _drift_profile(self, lo, hi, defaults):
        selected = (self.raw_mz >= lo) & (self.raw_mz <= hi)
        profile = np.bincount(self.raw_drift[selected], weights=self.raw_intensity[selected],
            minlength=self.n_drift_native)
        noise = float(defaults.get("noise_drift", 0))
        if noise > 0:
            profile[profile < profile.max(initial=0) * noise / 100] = 0
        profile = _smooth_axis(profile, defaults.get("smooth_drift", {"method": "raw"}))
        return _readonly(profile)

    def _grid(self, config):
        defaults = config.get("defaults", {})
        domain = self._domain(config)
        native_count = min(25600, max(200, int(np.searchsorted(self._unique_mz, domain[1], "right") - np.searchsorted(self._unique_mz, domain[0]))))
        mz_bins = defaults.get("mz_bins", 800)
        drift_bins = defaults.get("drift_bins", "native")
        mz_bins = native_count if mz_bins == "native" else int(mz_bins)
        drift_bins = self.n_drift_native if drift_bins == "native" else int(drift_bins)
        if mz_bins < 1 or drift_bins < 1:
            raise ValueError("The m/z and drift bin counts must be positive or Native.")
        if mz_bins * drift_bins > MAX_HEATMAP_CELLS:
            raise ValueError(f"The requested heatmap has {mz_bins * drift_bins:,} cells. "
                f"Choose fewer m/z or drift bins (limit {MAX_HEATMAP_CELLS:,} cells).")
        key = (domain, mz_bins, drift_bins)
        cached = self._grid_cache
        if cached is not None and cached[0] == key:
            return key, *cached[1:]
        mz_step = (domain[1] - domain[0]) / mz_bins
        mz_indices = np.floor((self.raw_mz - domain[0]) / mz_step).astype(np.int64)
        mz_indices[self.raw_mz == domain[1]] = mz_bins - 1
        # JS groups native integer drift bins over [floor(d*s),floor((d+1)*s)).
        # searchsorted reproduces that rule when fractional/upsampled bins exist.
        drift_edges, drift_centres = self._drift_axis(drift_bins)
        drift_indices = np.searchsorted(drift_edges, self.raw_drift, side="right") - 1
        valid = (mz_indices >= 0) & (mz_indices < mz_bins) & (drift_indices >= 0) & (drift_indices < drift_bins)
        flat = drift_indices[valid] * mz_bins + mz_indices[valid]
        raw = np.bincount(flat, weights=self.raw_intensity[valid], minlength=drift_bins * mz_bins).reshape(drift_bins, mz_bins)
        mz_centres = domain[0] + (np.arange(mz_bins) + 0.5) * mz_step
        cached = (key, _readonly(raw), _readonly(mz_centres), _readonly(drift_centres))
        self._grid_cache = cached
        return cached

    def render(self, config=None, viewport=None):
        cfg = deepcopy(self.config if config is None else config)
        defaults = cfg.get("defaults", {})
        key, raw, mz_centres, drift_centres = self._grid(cfg)
        domain, mz_bins, drift_bins = key
        options = {name: defaults.get(name) for name in ("smooth_2d", "scale", "noise_2d", "denoise", "z_scale")}
        display_key = (key, json.dumps(options, sort_keys=True))
        displayed = self._display_cache
        if displayed is None or displayed[0] != display_key:
            values, transform, label = _display_values(raw, domain, drift_centres, self.pusher_us, defaults)
            displayed = (display_key, values, transform, label)
            self._display_cache = displayed
        _display_key, values, transform, label = displayed
        initial = cfg.get("initial_view") or {}
        viewport = viewport or {}
        mz_range = _bounds(viewport.get("mz_range", initial.get("mz_range")), domain, "m/z range")
        drift_range = _bounds(viewport.get("drift_range", initial.get("drift_range")), self.full_drift_range, "drift range")
        mz_lo, mz_hi = max(domain[0], mz_range[0]), min(domain[1], mz_range[1])
        dr_lo, dr_hi = max(0, drift_range[0]), min(self.full_drift_range[1], drift_range[1])
        mz_width = (domain[1] - domain[0]) / mz_bins
        drift_width = self.full_drift_range[1] / drift_bins
        m_lo = max(0, math.floor((mz_lo - domain[0]) / mz_width)) if mz_hi > mz_lo else 0
        m_hi = min(mz_bins - 1, math.ceil((mz_hi - domain[0]) / mz_width) - 1) if mz_hi > mz_lo else -1
        d_lo = max(0, math.floor(dr_lo / drift_width)) if dr_hi > dr_lo else 0
        d_hi = min(drift_bins - 1, math.ceil(dr_hi / drift_width) - 1) if dr_hi > dr_lo else -1
        mz_slice, drift_slice = slice(m_lo, m_hi + 1), slice(d_lo, d_hi + 1)
        drift_profile = self._drift_profile(mz_lo, mz_hi, defaults)
        native_width = self.pusher_us / 1000 if self.pusher_us > 0 else 1
        native_lo = max(0, math.floor(_native_bin_position(dr_lo, native_width))) if dr_hi > dr_lo else 0
        native_hi = min(self.n_drift_native - 1, math.ceil(_native_bin_position(dr_hi, native_width)) - 1) if dr_hi > dr_lo else -1
        if native_lo == 0 and native_hi == self.n_drift_native - 1:
            full_ms_profile = self._full_ms_profile
        else:
            selected = (self.raw_drift >= native_lo) & (self.raw_drift <= native_hi)
            full_ms_profile = np.bincount(self._profile_indices[selected],
                weights=self.raw_intensity[selected], minlength=len(self._profile_mz))
        start = np.searchsorted(self._profile_mz, mz_lo, side="left")
        stop = np.searchsorted(self._profile_mz, mz_hi, side="right") if mz_hi >= mz_lo else start
        ms_mz = self._profile_mz[start:stop]
        ms_intensity = full_ms_profile[start:stop].copy()
        noise = float(defaults.get("noise_mz", 0))
        if noise > 0:
            ms_intensity[ms_intensity < ms_intensity.max(initial=0) * noise / 100] = 0
        ms_mode = defaults.get("ms_mode", "profile")
        if ms_mode == "stick" and len(ms_mz) > 2:
            threshold = ms_intensity.max(initial=0) * 0.01
            peaks, _ = find_peaks(ms_intensity, height=threshold)
            # Native TOF spacing varies with m/z: select local peaks in a real
            # mass neighborhood rather than smoothing or merging the samples.
            lower = np.searchsorted(ms_mz, ms_mz[peaks] - 0.08, side="left")
            upper = np.searchsorted(ms_mz, ms_mz[peaks] + 0.08, side="right")
            peaks = np.asarray([peak for peak, lo, hi in zip(peaks, lower, upper)
                if ms_intensity[peak] > threshold
                and ms_intensity[peak] == ms_intensity[lo:hi].max()], dtype=np.int64)
            stick_mz, stick_intensity = np.full((len(peaks), 3), np.nan), np.full((len(peaks), 3), np.nan)
            stick_mz[:, :2] = ms_mz[peaks, None]
            stick_intensity[:, 0], stick_intensity[:, 1] = 0, ms_intensity[peaks]
            ms_mz, ms_intensity = stick_mz.ravel(), stick_intensity.ravel()
        return ViewSnapshot(mz_centres, drift_centres, values, _readonly(drift_profile),
            _readonly(ms_mz), _readonly(ms_intensity), mz_range, drift_range,
            self.drift_unit, label, transform, raw, mz_slice, drift_slice, domain, ms_mode,
            self._native_drift_centres)


def load_view_data(im_path, ms_path=None, pusher_us=None, config=None):
    """Read TXT once; retain full floating-point precision for all native views."""
    im = np.loadtxt(str(local_path(im_path)), ndmin=2)
    if im.shape[1] < 3 or not np.isfinite(im[:, :3]).all():
        raise ValueError("The IM table must contain finite m/z, drift-bin and intensity columns.")
    positive = im[:, 2] > 0
    if not positive.any():
        raise ValueError("The IM table contains no positive intensity.")
    mz, drift, intensity = im[positive, 0], im[positive, 1].astype(np.int64), im[positive, 2]
    if (drift < 0).any():
        raise ValueError("The IM table must use nonnegative drift bins.")
    ms_mz = ms_intensity = np.empty(0)
    if ms_path is not None and local_path(ms_path).is_file():
        ms = np.loadtxt(str(local_path(ms_path)), ndmin=2)
        if ms.shape[1] < 2 or not np.isfinite(ms[:, :2]).all():
            raise ValueError("The MS table must contain finite m/z and intensity columns.")
        ms_mz, ms_intensity = ms[:, 0], ms[:, 1]
    return ViewData(mz, drift, intensity, ms_mz, ms_intensity, pusher_us, _configuration(config),
        profile_mz=im[:, 0], native_drift_count=int(im[:, 1].max()) + 1)
