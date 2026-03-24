"""CCS calibration — drift profile extraction and peak detection."""
from __future__ import annotations

import math
import re

import numpy as np
import pandas as pd



def parse_formula(formula: str) -> dict[str, int]:
    """Parse 'C10H15N5O10P2' → {'C':10, 'H':15, ...}."""
    elements = {}
    for m in re.finditer(r'([A-Z][a-z]?)(\d*)', formula):
        e, n = m.group(1), int(m.group(2)) if m.group(2) else 1
        if e:
            elements[e] = elements.get(e, 0) + n
    return elements


def compute_auto_mz_window(mw: float, z: int, formula: str = "") -> float:
    """Compute m/z extraction window from isotopic envelope (±3σ, min 2 Da).

    Uses exact atom counts if formula provided, else averagine model.
    """
    isotope_data = {
        "C": (0.0111, 1.00336), "H": (0.00012, 1.00628),
        "N": (0.00364, 0.99703), "O": (0.00038, 1.00422),
        "S": (0.0075, 0.99939), "Cl": (0.2423, 1.99705),
        "Br": (0.4931, 1.99796),
    }

    if formula and formula.strip():
        atoms = parse_formula(formula.strip())
        variance = 0.0
        for elem, count in atoms.items():
            if elem in isotope_data:
                p, dm = isotope_data[elem]
                if p > 0:
                    variance += count * p * (1 - p) * (dm / z) ** 2
        n_O = atoms.get("O", 0)
        if n_O > 0:
            variance += n_O * 0.00205 * 0.99795 * (2.00426 / z) ** 2
    else:
        n_atoms = mw / 14.0
        variance = n_atoms * 0.011 * 0.989 * (1.003 / z) ** 2

    return max(2.0, math.ceil(3.0 * math.sqrt(variance) * 2) / 2)


def extract_drift_profile(im_data: pd.DataFrame, target_mz: float,
                          mz_window: float) -> np.ndarray:
    """Sum intensity within ±mz_window of target_mz, grouped by drift bin."""
    mask = (im_data["mz"] >= target_mz - mz_window) & (im_data["mz"] <= target_mz + mz_window)
    subset = im_data[mask]
    max_bin = int(im_data["drift_bin"].max()) if not im_data.empty else 0
    if subset.empty:
        return np.zeros(max_bin + 1)
    profile = np.zeros(max_bin + 1)
    for _, row in subset.iterrows():
        b = int(row["drift_bin"])
        if 0 <= b <= max_bin:
            profile[b] += row["intensity"]
    return profile


def detect_peaks_in_profile(profile: np.ndarray, pusher_us: float,
                            min_height_frac: float = 0.10,
                            min_distance: int = 3) -> list[dict]:
    """Detect peaks in a 1D drift profile. Returns list sorted by drift bin."""
    from scipy.signal import find_peaks as _find_peaks
    if profile.max() == 0:
        return []
    threshold = profile.max() * min_height_frac
    peak_indices, _ = _find_peaks(profile, height=threshold, distance=min_distance)
    if len(peak_indices) == 0:
        return []
    results = []
    for idx in peak_indices:
        results.append({
            "bin_index": int(idx),
            "drift_bin": int(idx),
            "drift_time_ms": (float(idx) + 0.5) * pusher_us / 1000.0,
            "peak_intensity": float(profile[idx]),
        })
    results.sort(key=lambda r: r["drift_bin"])
    return results


def select_peak(all_peaks: list[dict], mode: str = "highest") -> dict:
    """Select one peak: 'highest' (tallest) or 'last' (most extended)."""
    if mode == "last":
        return all_peaks[-1]
    return max(all_peaks, key=lambda p: p["peak_intensity"])
