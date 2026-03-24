"""Shared I/O utilities for IM-MS data.

Generic functions that work with any vendor's converted text files.
Vendor-specific readers (STS binary, _extern.inf) are in deconvovo/vendors/.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


# =============================================================================
# Physical constants (shared across all modules)
# =============================================================================

H_MASS = 1.00728     # proton mass (Da) — for m/z ↔ MW conversion
GAS_MW_N2 = 28.014   # N2 default drift gas mass (Da)


def mz_from_mw(mw: float, z: int) -> float:
    """Compute m/z for protonated species: [M + zH]^z+."""
    return (mw + z * H_MASS) / z


def mw_from_mz(mz: float, z: int) -> float:
    """Compute neutral MW from m/z for protonated species: [M + zH]^z+."""
    return mz * z - z * H_MASS


# =============================================================================
# Generic text file I/O (vendor-agnostic — any IM-MS data in this format)
# =============================================================================

def read_im_txt(im_path: Path) -> pd.DataFrame:
    """Read IM data: 3-column text (m/z, drift_bin, intensity)."""
    data = np.loadtxt(str(im_path))
    if data.ndim != 2 or data.shape[1] < 3:
        return pd.DataFrame(columns=["mz", "drift_bin", "intensity"])
    return pd.DataFrame({
        "mz": data[:, 0],
        "drift_bin": data[:, 1].astype(int),
        "intensity": data[:, 2],
    })


def read_ms_txt(ms_path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Read MS data: 2-column text (m/z, intensity)."""
    data = np.loadtxt(str(ms_path))
    if data.ndim != 2 or data.shape[1] < 2:
        return np.empty(0), np.empty(0)
    return data[:, 0], data[:, 1]


def measure_ms_intensity(ms_file: Path, target_mz: float,
                         mz_window: float) -> float:
    """Max MS intensity within ±mz_window of target_mz."""
    if not ms_file.exists():
        return 0.0
    data = np.loadtxt(str(ms_file))
    if data.ndim != 2 or data.shape[1] < 2:
        return 0.0
    mz_arr, int_arr = data[:, 0], data[:, 1]
    mask = (mz_arr >= target_mz - mz_window) & (mz_arr <= target_mz + mz_window)
    return float(int_arr[mask].max()) if mask.any() else 0.0


# =============================================================================
# Vendor-dispatching helpers
# =============================================================================

def find_pusher_period(raw_path: Path) -> float | None:
    """Find pusher period from a .raw directory. Dispatches to vendor reader."""
    from deconvovo.vendors.waters import find_pusher_period as _waters_pp
    return _waters_pp(raw_path)


def read_extern_inf(raw_path: Path) -> dict[str, str]:
    """Read instrument parameters. Dispatches to vendor reader."""
    from deconvovo.vendors.waters import read_extern_inf as _waters_inf
    return _waters_inf(raw_path)


def find_pusher_for_run(data_dir: Path, run_name: str,
                        raw_dir: Path | None = None) -> float | None:
    """Try to find pusher period by searching common locations."""
    for search_dir in [data_dir, data_dir.parent]:
        rp = search_dir / f"{run_name}.raw"
        if rp.is_dir():
            pp = find_pusher_period(rp)
            if pp is not None:
                return pp
    if raw_dir:
        rp = raw_dir / f"{run_name}.raw"
        if rp.is_dir():
            pp = find_pusher_period(rp)
            if pp is not None:
                return pp
    return None
