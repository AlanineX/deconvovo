"""Reusable 1D smoothing utilities for IM-MS profiles.

Mirrors the smoothing methods available in the interactive HTML viewer
(Gaussian, Savitzky-Golay, moving average) so that offline/batch processing
produces identical results.

Usage:
    from deconvovo.smooth import smooth1d, apply_noise_floor

    y_smooth = smooth1d(y, method="gaussian", sigma=1.5)
    y_clean  = apply_noise_floor(y, pct=1.0)
"""
from __future__ import annotations

import numpy as np
from scipy.ndimage import uniform_filter1d
from scipy.signal import savgol_filter


def gaussian1d(y: np.ndarray, sigma: float) -> np.ndarray:
    """Gaussian smoothing (matches JS gaussSmooth1D in the HTML viewer).

    Uses reflect-mode boundary handling so endpoints are not distorted.
    """
    if sigma <= 0 or len(y) < 3:
        return y.copy()
    from scipy.ndimage import gaussian_filter1d
    return np.maximum(0, gaussian_filter1d(y.astype(float), sigma=sigma, mode="reflect"))


def savgol1d(y: np.ndarray, window: int, polyorder: int = 2) -> np.ndarray:
    """Savitzky-Golay smoothing (matches JS sgSmooth1D).

    Preserves peak shape better than Gaussian at the cost of possible
    ringing on very noisy data.
    """
    if window < 3 or len(y) < window:
        return y.copy()
    if window % 2 == 0:
        window += 1
    polyorder = min(polyorder, window - 1)
    return np.maximum(0, savgol_filter(y.astype(float), window, polyorder))


def moving_average1d(y: np.ndarray, window: int) -> np.ndarray:
    """Moving average smoothing (matches JS maSmooth1D)."""
    if window < 2 or len(y) < window:
        return y.copy()
    return np.maximum(0, uniform_filter1d(y.astype(float), size=window, mode="reflect"))


def apply_noise_floor(y: np.ndarray, pct: float) -> np.ndarray:
    """Zero out values below pct% of the maximum (noise threshold).

    Applied BEFORE smoothing so boundaries get blurred naturally.
    """
    if pct <= 0:
        return y.copy()
    threshold = y.max() * pct / 100.0
    out = y.copy()
    out[out < threshold] = 0
    return out


def smooth1d(
    y: np.ndarray,
    method: str = "raw",
    sigma: float = 1.0,
    window: int = 5,
    polyorder: int = 2,
    noise_pct: float = 0.0,
) -> np.ndarray:
    """Dispatch smoothing by method name.

    This is the single entry point for all 1D smoothing. Method names
    match the HTML viewer's preset format.

    Args:
        y: Input 1D array.
        method: "raw", "gaussian", "sg", or "ma".
        sigma: Gaussian sigma (only for method="gaussian").
        window: Window size (for "sg" and "ma").
        polyorder: Polynomial order (for "sg").
        noise_pct: Noise floor as % of max, applied before smoothing.

    Returns:
        Smoothed 1D array (same length as y).
    """
    out = y.astype(float).copy()

    # Apply noise floor first (matches HTML viewer order)
    if noise_pct > 0:
        out = apply_noise_floor(out, noise_pct)

    if method == "gaussian":
        return gaussian1d(out, sigma)
    elif method == "sg":
        return savgol1d(out, window, polyorder)
    elif method == "ma":
        return moving_average1d(out, window)
    else:  # "raw"
        return out
