"""CCS calibration — apply calibration to drift profiles, resample, smooth."""
from __future__ import annotations

import math

import numpy as np

from deconvovo.smooth import smooth1d



def apply_ccs(drift_time_ms: float, mz: float, z: int, mw: float,
              cal: dict, method: str = "direct") -> float:
    """Convert a single drift time to CCS."""
    edc = cal["edc_coeff"]
    X = cal["X"]
    tp = drift_time_ms - edc * math.sqrt(mz) / 1000.0
    if tp <= 0:
        return float("nan")
    mu = math.sqrt(1.0 / mw + 1.0 / cal["gas_mw"])
    if method == "twostep":
        tdp = (tp ** X) * z * mu
        return cal["slope_diagnostic"] * tdp + cal["intercept_diagnostic"]
    else:
        return cal["A"] * (tp ** X) * z * mu


def convert_profile_to_ccs(profile: np.ndarray, pusher_us: float,
                            mz: float, z: int, mw: float, cal: dict,
                            method: str = "direct"):
    """Convert full drift profile to CCS space. Returns (bins, dt_ms, ccs, intensity)."""
    n = len(profile)
    bins = np.arange(n)
    dt_ms = (bins + 0.5) * pusher_us / 1000.0  # bin centers
    ccs_arr = np.full(n, np.nan)
    for i in range(n):
        ccs_arr[i] = apply_ccs(dt_ms[i], mz, z, mw, cal, method)
    valid = np.isfinite(ccs_arr)
    return bins[valid], dt_ms[valid], ccs_arr[valid], profile[valid]


def resample_to_uniform_ccs(ccs, intensity, step=0.5):
    """Resample to uniform CCS grid (anchored at multiples of step)."""
    if len(ccs) < 2:
        return ccs.copy(), intensity.copy()
    ccs_lo = math.floor(float(ccs.min()) / step) * step
    ccs_hi = math.ceil(float(ccs.max()) / step) * step
    ccs_uni = np.arange(ccs_lo, ccs_hi + step * 0.5, step)
    int_uni = np.interp(ccs_uni, ccs, intensity)
    return ccs_uni, int_uni


def smooth_ccs_profile(ccs_uni, int_uni, sigma=3.0):
    """Gaussian smooth on uniform CCS grid."""
    return smooth1d(int_uni, method="gaussian", sigma=sigma)


def significant_range(ccs, intensity, pad_pct=15.0):
    """Find display range covering significant intensity + padding."""
    if len(ccs) == 0 or intensity.max() == 0:
        return (0, 1)
    threshold = intensity.max() * 0.005
    sig = ccs[intensity > threshold]
    if len(sig) == 0:
        return (float(ccs.min()), float(ccs.max()))
    lo, hi = float(sig.min()), float(sig.max())
    span = hi - lo if hi > lo else 1.0
    pad = span * pad_pct / 100.0
    return (lo - pad, hi + pad)
