#!/usr/bin/env python3
"""
MassLynx .raw binary reader and analysis library for Waters SYNAPT G2 IM-MS data.

Reads MassLynx binary formats directly without vendor libraries:
  - IDX files: 22-byte scan index records (offset, npoints, flags)
  - DAT files: 8-byte (intensity_u32, channel_u32) spectral records
  - FUNC001: retention-time-resolved MS data (summed for overview spectrum)
  - FUNC002: drift-time-resolved IMS data (200 bins typical for SYNAPT)

Calibration: sqrt(m/z) = k * channel + d  (two-point linear in sqrt space)

Species databases are configurable: built-in presets for common analytes
or user-provided JSON files.
"""
from __future__ import annotations

import json
import re
import struct
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd

# =============================================================================
# Physical constants (monoisotopic)
# =============================================================================
H = 1.00728
NA = 22.98922
K = 38.96316
NH4 = 18.03437
H2O = 18.01056
EDA = 60.0531  # ethylenediamine C2H8N2


# =============================================================================
# Data classes
# =============================================================================

@dataclass
class ScanData:
    """One scan (or drift bin) worth of raw channel/intensity data."""
    channels: np.ndarray
    intensities: np.ndarray
    total_intensity: float


@dataclass
class Calibration:
    """Calibration object. With the rainbow-based reader, data is already
    calibrated by the vendor polynomial, so to_mz() is the identity.
    The k/d fields are kept for backward compatibility."""
    k: float = 1.0
    d: float = 0.0
    method: str = "vendor_polynomial"

    def to_mz(self, mz_values: np.ndarray) -> np.ndarray:
        """Identity — data from read_func_scans is already in calibrated m/z."""
        return mz_values


# =============================================================================
# Species database
# =============================================================================

def build_adp_species_db(
    neutral: float = 427.0294,
    eda_mass: float = EDA,
) -> dict[str, float]:
    """ADP species: free, salt adducts, EDA complexes, dimers."""
    db: dict[str, float] = {}
    # Monomers
    db["ADP+H"] = neutral + H
    db["ADP+Na"] = neutral + NA
    db["ADP+K"] = neutral + K
    db["ADP+NH4"] = neutral + NH4
    db["ADP+2Na-H"] = neutral + 2 * NA - H
    db["ADP+3Na-2H"] = neutral + 3 * NA - 2 * H
    db["ADP-H2O+H"] = neutral - H2O + H
    db["ADP-H2O+Na"] = neutral - H2O + NA
    # Doubly charged monomers
    db["ADP+2H [2+]"] = (neutral + 2 * H) / 2
    db["ADP+H+Na [2+]"] = (neutral + H + NA) / 2
    # EDA complexes
    db["ADP+EDA+H"] = neutral + eda_mass + H
    db["ADP+EDA+Na"] = neutral + eda_mass + NA
    db["ADP+EDA+K"] = neutral + eda_mass + K
    db["ADP+EDA+2H [2+]"] = (neutral + eda_mass + 2 * H) / 2
    db["ADP+EDA+H+Na [2+]"] = (neutral + eda_mass + H + NA) / 2
    db["ADP+2EDA+H"] = neutral + 2 * eda_mass + H
    db["ADP+2EDA+2H [2+]"] = (neutral + 2 * eda_mass + 2 * H) / 2
    db["ADP+2EDA+H+Na [2+]"] = (neutral + 2 * eda_mass + H + NA) / 2
    # Dimers
    db["2ADP+H"] = 2 * neutral + H
    db["2ADP+Na"] = 2 * neutral + NA
    db["2ADP+K"] = 2 * neutral + K
    db["2ADP+2Na-H"] = 2 * neutral + 2 * NA - H
    db["2ADP+2H [2+]"] = (2 * neutral + 2 * H) / 2
    db["2ADP+H+Na [2+]"] = (2 * neutral + H + NA) / 2
    db["2ADP+EDA+H"] = 2 * neutral + eda_mass + H
    db["2ADP+EDA+Na"] = 2 * neutral + eda_mass + NA
    db["2ADP+EDA+2H [2+]"] = (2 * neutral + eda_mass + 2 * H) / 2
    # Trimers
    db["3ADP+H"] = 3 * neutral + H
    db["3ADP+Na"] = 3 * neutral + NA
    # Free EDA
    db["EDA+H"] = eda_mass + H
    return db


def build_protein_species_db(
    name: str,
    neutral_mass: float,
    min_z: int = 5,
    max_z: int = 25,
) -> dict[str, float]:
    """Generate expected m/z for each charge state of a protein."""
    return {
        f"{name} [{z}+]": (neutral_mass + z * H) / z
        for z in range(min_z, max_z + 1)
    }


BUILTIN_PRESETS: dict[str, Callable[[], dict[str, float]]] = {
    "adp": build_adp_species_db,
    "ubiquitin": lambda: build_protein_species_db("UBQ", 8565.8, 5, 15),
    "myoglobin_apo": lambda: build_protein_species_db("MYO_apo", 16951.5, 8, 25),
    "myoglobin_holo": lambda: build_protein_species_db("MYO_holo", 17567.5, 8, 25),
    "cytochrome_c": lambda: build_protein_species_db("CYTC", 12384.0, 7, 20),
}


def load_species_db(source: str | Path) -> dict[str, float]:
    """Load species DB from a JSON file path or a named preset."""
    p = Path(source)
    if p.exists() and p.suffix == ".json":
        with open(p) as f:
            return json.load(f)
    name = str(source).lower().replace("-", "_")
    if name in BUILTIN_PRESETS:
        return BUILTIN_PRESETS[name]()
    raise ValueError(
        f"Unknown species DB '{source}'. "
        f"Use a JSON file or one of: {', '.join(BUILTIN_PRESETS)}"
    )


# =============================================================================
# Binary reading
# =============================================================================

def read_header(raw_dir: Path) -> dict[str, str]:
    """Parse _HEADER.TXT metadata."""
    path = raw_dir / "_HEADER.TXT"
    meta: dict[str, str] = {}
    if not path.exists():
        return meta
    with open(path, errors="replace") as f:
        for line in f:
            line = line.strip()
            if line.startswith("$$"):
                parts = line[2:].strip().split(":", 1)
                if len(parts) == 2:
                    meta[parts[0].strip()] = parts[1].strip()
    return meta


def read_extern_inf(raw_dir: Path) -> dict[str, str]:
    """Parse _extern.inf instrument parameters."""
    path = raw_dir / "_extern.inf"
    params: dict[str, str] = {}
    if not path.exists():
        return params
    with open(path, errors="replace") as f:
        for line in f:
            parts = [p.strip() for p in line.strip().split("\t") if p.strip()]
            if len(parts) >= 2:
                params[parts[0]] = parts[-1]
    return params


def parse_vendor_calibration(raw_dir: Path) -> list[float]:
    """Extract vendor calibration polynomial from _HEADER.TXT Cal Function 1."""
    header = read_header(raw_dir)
    cal_str = header.get("Cal Function 1", "")
    if not cal_str:
        return []
    coeffs: list[float] = []
    for part in cal_str.split(","):
        part = part.strip()
        try:
            coeffs.append(float(part))
        except ValueError:
            pass  # skip 'T1' suffix
    return coeffs


def read_ms_txt(ms_path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Read a CDCReader MS output file (2-column: m/z intensity)."""
    data = np.loadtxt(str(ms_path))
    if data.ndim != 2 or data.shape[1] < 2:
        return np.empty(0), np.empty(0)
    return data[:, 0], data[:, 1]


def read_im_txt(im_path: Path) -> pd.DataFrame:
    """Read a CDCReader IM output file (3-column: m/z drift_bin intensity)."""
    data = np.loadtxt(str(im_path))
    if data.ndim != 2 or data.shape[1] < 3:
        return pd.DataFrame(columns=["mz", "drift_bin", "intensity"])
    return pd.DataFrame({"mz": data[:, 0], "drift_bin": data[:, 1].astype(int), "intensity": data[:, 2]})


def read_im_as_scans(im_path: Path) -> list[ScanData]:
    """Read CDCReader IM output and group by drift bin into ScanData list."""
    df = read_im_txt(im_path)
    if df.empty:
        return []

    max_bin = int(df["drift_bin"].max())
    scans: list[ScanData] = []
    for i in range(max_bin + 1):
        subset = df[df["drift_bin"] == i]
        nonzero = subset[subset["intensity"] > 0]
        if len(nonzero) > 0:
            mz_arr = nonzero["mz"].values.astype(np.float64)
            int_arr = nonzero["intensity"].values.astype(np.float64)
        else:
            mz_arr = np.empty(0, dtype=np.float64)
            int_arr = np.empty(0, dtype=np.float64)
        scans.append(ScanData(mz_arr, int_arr, float(int_arr.sum())))
    return scans


def read_func_scans(raw_dir: Path, func: int = 1) -> list[ScanData]:
    """Read Waters data. Prefers CDCReader-converted text files if available.

    Looks for {raw_stem}_ms.txt and {raw_stem}_im.txt in the parent directory.
    Falls back to rainbow binary reader if text files are not found.
    """
    # Check for converted text files
    raw_stem = raw_dir.stem
    parent = raw_dir.parent

    # For func=2 (IMS), try IM text file
    if func == 2:
        im_path = parent / f"{raw_stem}_im.txt"
        if im_path.exists():
            return read_im_as_scans(im_path)

    # For func=1 (MS), the text file is a summed spectrum, not per-scan
    # Return a single ScanData with the full spectrum
    if func == 1:
        ms_path = parent / f"{raw_stem}_ms.txt"
        if ms_path.exists():
            mz, intensity = read_ms_txt(ms_path)
            if len(mz) > 0:
                return [ScanData(mz, intensity, float(intensity.sum()))]

    # Fallback: try rainbow binary reader
    try:
        from rainbow.waters.masslynx import parse_funcidx, parse_funcdat8

        idx_path = str(raw_dir / f"_FUNC{func:03d}.IDX")
        dat_path = str(raw_dir / f"_FUNC{func:03d}.DAT")
        if not Path(idx_path).exists() or not Path(dat_path).exists():
            return []

        calib = parse_vendor_calibration(raw_dir)
        times, pair_counts, bytes_per_pair = parse_funcidx(idx_path)
        if bytes_per_pair != 8:
            return []

        ylabels, data = parse_funcdat8(
            dat_path, pair_counts, prec=4,
            calib=calib if calib else None,
        )
        valid_cols = ylabels <= 5000.0
        ylabels = ylabels[valid_cols]
        data = data[:, valid_cols]

        scans: list[ScanData] = []
        for i in range(data.shape[0]):
            row = data[i]
            nz = row > 0
            mz_arr = ylabels[nz].astype(np.float64) if nz.any() else np.empty(0, dtype=np.float64)
            int_arr = row[nz].astype(np.float64) if nz.any() else np.empty(0, dtype=np.float64)
            scans.append(ScanData(mz_arr, int_arr, float(int_arr.sum())))
        return scans
    except Exception:
        return []


def sum_scans(scans: list[ScanData]) -> tuple[np.ndarray, np.ndarray]:
    """Sum all scans into a single (m/z, intensities) pair, averaged.

    NOTE: With the rainbow-based reader, channels are already calibrated m/z values.
    """
    mz_sum: dict[float, float] = {}
    for s in scans:
        for mz_val, intens in zip(s.channels, s.intensities):
            key = round(float(mz_val), 4)
            mz_sum[key] = mz_sum.get(key, 0.0) + intens

    if not mz_sum:
        return np.empty(0, dtype=np.float64), np.empty(0, dtype=np.float64)

    sorted_mz = sorted(mz_sum.keys())
    mz_arr = np.array(sorted_mz, dtype=np.float64)
    int_arr = np.array([mz_sum[m] for m in sorted_mz], dtype=np.float64)
    if len(scans) > 0:
        int_arr /= len(scans)
    return mz_arr, int_arr


# =============================================================================
# Calibration
# =============================================================================

def calibrate_two_point(ch1: float, mz1: float, ch2: float, mz2: float) -> Calibration:
    """Two-point sqrt(m/z) calibration."""
    k = (np.sqrt(mz2) - np.sqrt(mz1)) / (ch2 - ch1)
    d = np.sqrt(mz1) - k * ch1
    return Calibration(k, d, "two_point")


def _find_peak_clusters(
    channels: np.ndarray,
    intensities: np.ndarray,
    min_int: float = 100,
    gap_threshold: float = 500_000,
) -> list[dict]:
    """Group channels into intensity clusters separated by gaps."""
    nz = intensities > min_int
    if nz.sum() == 0:
        return []

    nz_idx = np.where(nz)[0]
    clusters: list[tuple[int, int]] = []
    start = nz_idx[0]
    prev = nz_idx[0]
    for idx in nz_idx[1:]:
        if channels[idx] - channels[prev] > gap_threshold:
            clusters.append((start, prev))
            start = idx
        prev = idx
    clusters.append((start, prev))

    result = []
    for s, e in clusters:
        sl = slice(s, e + 1)
        cl_ch = channels[sl]
        cl_int = intensities[sl]
        apex = int(np.argmax(cl_int))
        result.append({
            "center_ch": cl_ch[apex],
            "max_int": cl_int[apex],
            "total_int": float(np.sum(cl_int)),
            "n_points": len(cl_ch),
        })
    return sorted(result, key=lambda x: -x["max_int"])


def calibrate_adp(
    channels: np.ndarray,
    intensities: np.ndarray,
    anchor_mz: tuple[float, float] = (428.037, 450.019),
) -> Calibration:
    """Calibrate using two brightest clusters as ADP anchors."""
    clusters = _find_peak_clusters(channels, intensities, min_int=100)
    if len(clusters) < 2:
        return _calibrate_fallback(channels)

    top2 = sorted(clusters[:2], key=lambda c: c["center_ch"])
    return calibrate_two_point(
        top2[0]["center_ch"], anchor_mz[0],
        top2[1]["center_ch"], anchor_mz[1],
    )._replace_method("adp_two_point")


def calibrate_protein_selfconsistent(
    channels: np.ndarray,
    intensities: np.ndarray,
    expected_masses: list[float] | None = None,
    mass_range: tuple[float, float] = (5000, 50000),
) -> Calibration:
    """Self-consistent protein calibration using charge-state spacing."""
    clusters = _find_peak_clusters(channels, intensities, min_int=200, gap_threshold=300_000)
    if len(clusters) < 2:
        return _calibrate_fallback(channels)

    top_clusters = clusters[: min(15, len(clusters))]

    if expected_masses is None:
        expected_masses = [8565.8, 16951.5, 17567.5, 12384.0]

    trial_masses: list[float] = []
    for m in expected_masses:
        for delta in [-50, -20, 0, 20, 50]:
            trial_masses.append(m + delta)
    trial_masses.extend(np.arange(mass_range[0], mass_range[1], 500).tolist())

    best_result = None
    best_score = 0.0

    for M_trial in trial_masses:
        for i in range(min(5, len(top_clusters))):
            for j in range(min(5, len(top_clusters))):
                if i == j:
                    continue
                cl_a, cl_b = top_clusters[i], top_clusters[j]
                if cl_a["center_ch"] >= cl_b["center_ch"]:
                    continue

                for z_a in range(5, 25):
                    z_b = z_a - 1
                    if z_b < 3:
                        continue
                    mz_a = (M_trial + z_a * H) / z_a
                    mz_b = (M_trial + z_b * H) / z_b
                    if mz_a <= 0 or mz_b <= 0 or mz_a > 5000 or mz_b > 5000 or mz_a >= mz_b:
                        continue

                    dch = cl_b["center_ch"] - cl_a["center_ch"]
                    if dch <= 0:
                        continue
                    k = (np.sqrt(mz_b) - np.sqrt(mz_a)) / dch
                    d_ = np.sqrt(mz_a) - k * cl_a["center_ch"]
                    if k <= 0:
                        continue

                    n_consistent = 0
                    total_int = 0.0
                    for cl in top_clusters:
                        sq = k * cl["center_ch"] + d_
                        if sq <= 0:
                            continue
                        mz_cl = sq ** 2
                        for z_try in range(3, 30):
                            mz_exp = (M_trial + z_try * H) / z_try
                            if abs(mz_cl - mz_exp) < 2.0:
                                n_consistent += 1
                                total_int += cl["max_int"]
                                break

                    score = n_consistent * 1000 + total_int * 0.01
                    if n_consistent >= 3 and score > best_score:
                        best_score = score
                        best_result = (k, d_)

    if best_result:
        cal = Calibration(best_result[0], best_result[1], "protein_selfconsistent")
        # Sanity check: the calibrated m/z range must be reasonable.
        # A SYNAPT acquisition typically spans 50-2000+ m/z.
        # If the calibration compresses the full channel range into < 200 Da,
        # it found a spurious local minimum — reject and fall back.
        mz_test = cal.to_mz(channels)
        mz_valid = mz_test[(mz_test > 0) & (mz_test < 50000) & np.isfinite(mz_test)]
        if len(mz_valid) > 0:
            mz_span = float(mz_valid.max() - mz_valid.min())
            if mz_span < 200:
                return _calibrate_fallback(channels)
        return cal
    return _calibrate_fallback(channels)


def _calibrate_fallback(channels: np.ndarray) -> Calibration:
    """Linear calibration mapping full channel range to 50-2000 m/z (sqrt space)."""
    if len(channels) < 2:
        return Calibration(1e-9, 0.0, "fallback")
    k = (np.sqrt(2000) - np.sqrt(50)) / (channels.max() - channels.min())
    d = np.sqrt(50) - k * channels.min()
    return Calibration(k, d, "fallback")


def auto_calibrate(
    channels: np.ndarray,
    intensities: np.ndarray,
    run_class: str,
) -> Calibration:
    """Dispatch calibration by run class.

    For ADP runs: two-point calibration on ADP+H and ADP+Na peaks.
    For protein/calibrant runs: self-consistent charge-state calibration,
    with automatic fallback if the result is unreasonable.
    """
    if "protein" in run_class or "calibrant" in run_class:
        masses = {
            "protein_ubq": [8565.8],
            "calibrant_ubq": [8565.8],
            "protein_myo": [16951.5, 17567.5],
            "protein_cytc": [12384.0],
        }
        return calibrate_protein_selfconsistent(
            channels, intensities,
            expected_masses=masses.get(run_class),
        )
    return calibrate_adp(channels, intensities)


# Patch Calibration to support _replace_method
def _replace_method(self: Calibration, method: str) -> Calibration:
    return Calibration(self.k, self.d, method)
Calibration._replace_method = _replace_method  # type: ignore[attr-defined]


# =============================================================================
# Peak detection and assignment
# =============================================================================

def detect_peaks(
    mz: np.ndarray,
    intensity: np.ndarray,
    min_int_frac: float = 0.02,
    min_abs_int: float = 50,
) -> list[dict]:
    """Detect local maxima above threshold."""
    if len(mz) < 3:
        return []
    threshold = max(min_abs_int, min_int_frac * intensity.max())
    peaks = []
    for i in range(1, len(mz) - 1):
        if intensity[i] > threshold and intensity[i] >= intensity[i - 1] and intensity[i] > intensity[i + 1]:
            peaks.append({"mz": mz[i], "intensity": intensity[i], "idx": i})
    return peaks


def group_isotopes(
    peaks: list[dict],
    tol: float = 0.5,
    spacing: float = 1.003,
) -> list[dict]:
    """Group peaks into isotope envelopes."""
    if not peaks:
        return []
    peaks_s = sorted(peaks, key=lambda p: p["mz"])
    used: set[int] = set()
    groups = []

    for i, pk in enumerate(peaks_s):
        if i in used:
            continue
        group = [pk]
        used.add(i)
        last_mz = pk["mz"]
        for j in range(i + 1, len(peaks_s)):
            if j in used:
                continue
            delta = peaks_s[j]["mz"] - last_mz
            if delta < spacing - tol:
                continue
            if delta > spacing + tol:
                break
            group.append(peaks_s[j])
            used.add(j)
            last_mz = peaks_s[j]["mz"]

        groups.append({
            "mono_mz": group[0]["mz"],
            "max_int": max(p["intensity"] for p in group),
            "total_area": sum(p["intensity"] for p in group),
            "n_iso": len(group),
        })

    return sorted(groups, key=lambda g: -g["max_int"])


def assign_species(
    groups: list[dict],
    species_db: dict[str, float],
    tol: float = 0.8,
) -> list[dict]:
    """Assign known species to isotope groups. Works with any species_db."""
    assignments = []
    for g in groups:
        mono = g["mono_mz"]
        best_name: str | None = None
        best_delta = tol
        best_exp: float | None = None

        for name, exp_mz in species_db.items():
            delta = abs(mono - exp_mz)
            if delta < best_delta:
                best_delta = delta
                best_name = name
                best_exp = exp_mz

        charge = 2 if best_name and "[2+]" in best_name else 1
        assignments.append({
            "mono_mz": mono,
            "max_int": g["max_int"],
            "total_area": g["total_area"],
            "n_iso": g["n_iso"],
            "species": best_name or "unassigned",
            "expected_mz": best_exp,
            "mass_error": best_delta if best_name else None,
            "charge": charge,
        })
    return assignments


def deconvolve_protein(
    mz: np.ndarray,
    intensity: np.ndarray,
    expected_mass: float | None = None,
    min_z: int = 4,
    max_z: int = 25,
    mz_tol: float = 3.0,
) -> list[dict]:
    """Deconvolve protein charge-state envelope."""
    peaks = detect_peaks(mz, intensity, min_int_frac=0.05, min_abs_int=200)
    if len(peaks) < 3:
        return []

    peaks_by_int = sorted(peaks, key=lambda p: -p["intensity"])
    results = []
    used: set[int] = set()

    for seed in peaks_by_int[:20]:
        if id(seed) in used:
            continue

        best_series = None
        best_score = 0.0

        for z_seed in range(min_z, max_z + 1):
            M = seed["mz"] * z_seed - z_seed * H
            if expected_mass and abs(M - expected_mass) > 500:
                continue
            if M < 1000 or M > 100000:
                continue

            series = [(seed["mz"], z_seed, seed["intensity"])]
            for z in range(min_z, max_z + 1):
                if z == z_seed:
                    continue
                exp_mz = (M + z * H) / z
                if exp_mz < mz.min() or exp_mz > mz.max():
                    continue
                for pk in peaks:
                    if abs(pk["mz"] - exp_mz) < mz_tol:
                        series.append((pk["mz"], z, pk["intensity"]))
                        break

            if len(series) >= 3:
                score = len(series) * sum(s[2] for s in series)
                if score > best_score:
                    best_score = score
                    best_series = series

        if best_series:
            masses = [s[0] * s[1] - s[1] * H for s in best_series]
            results.append({
                "neutral_mass": float(np.mean(masses)),
                "mass_std": float(np.std(masses)),
                "n_charges": len(best_series),
                "charge_min": min(s[1] for s in best_series),
                "charge_max": max(s[1] for s in best_series),
                "max_intensity": max(s[2] for s in best_series),
                "total_intensity": sum(s[2] for s in best_series),
                "series": best_series,
            })
            for s in best_series:
                for pk in peaks:
                    if abs(pk["mz"] - s[0]) < 1.0:
                        used.add(id(pk))

    return results


# =============================================================================
# IMS: 2D data and per-species drift extraction
# =============================================================================

def build_2d_data(
    scans: list[ScanData],
    cal: Calibration,
    mz_range: tuple[float, float] = (50, 2000),
    n_bins: int = 800,
) -> tuple[np.ndarray, np.ndarray]:
    """Build 2D heatmap: drift_bins x mz_bins.

    Args:
        n_bins: Number of m/z bins. 800 gives ~1 Da bins for typical IM-MS
            ranges, which balances resolution with fill density. Use higher
            values (2000-4000) only for CSV export where visual fill doesn't
            matter.
    """
    n_drift = len(scans)
    mz_edges = np.linspace(mz_range[0], mz_range[1], n_bins + 1)
    n_mz = len(mz_edges) - 1
    heatmap = np.zeros((n_drift, n_mz), dtype=np.float64)

    for drift_idx, s in enumerate(scans):
        if len(s.channels) == 0:
            continue
        mz = cal.to_mz(s.channels)
        bin_idx = np.searchsorted(mz_edges, mz) - 1
        valid = (bin_idx >= 0) & (bin_idx < n_mz)
        for bi, intens in zip(bin_idx[valid], s.intensities[valid]):
            heatmap[drift_idx, bi] += intens

    mz_centers = (mz_edges[:-1] + mz_edges[1:]) / 2
    return heatmap, mz_centers


def build_2d_data_hires(
    scans: list[ScanData],
    cal: Calibration,
    mz_range: tuple[float, float] = (50, 2000),
    max_bins: int = 4000,
) -> tuple[np.ndarray, np.ndarray]:
    """Build high-resolution 2D data for CSV export (adaptive binning)."""
    n_drift = len(scans)
    all_mz: set[float] = set()
    for s in scans:
        if len(s.channels) > 0:
            vals = cal.to_mz(s.channels)
            in_range = vals[(vals >= mz_range[0]) & (vals <= mz_range[1])]
            all_mz.update(in_range.tolist())

    if len(all_mz) > 2:
        sorted_mz = np.sort(np.array(list(all_mz)))
        diffs = np.diff(sorted_mz)
        pos_diffs = diffs[diffs > 0]
        if len(pos_diffs) > 0:
            native_res = float(np.median(pos_diffs))
            target_width = max(2.0 * native_res, (mz_range[1] - mz_range[0]) / max_bins)
            actual_bins = max(500, min(int((mz_range[1] - mz_range[0]) / target_width), max_bins))
        else:
            actual_bins = 800
    else:
        actual_bins = 800

    return build_2d_data(scans, cal, mz_range, n_bins=actual_bins)


def extract_species_drift_profiles(
    scans: list[ScanData],
    cal: Calibration,
    species_db: dict[str, float],
    mz_tol: float = 1.0,
) -> pd.DataFrame:
    """Extract per-species intensity across drift bins.

    For each drift bin, sums intensity within mz_tol of each species' expected m/z.
    Returns DataFrame: drift_bin | total_tic | species_1 | species_2 | ...
    """
    n_drift = len(scans)
    tic = np.zeros(n_drift)
    profiles: dict[str, np.ndarray] = {name: np.zeros(n_drift) for name in species_db}

    for i, s in enumerate(scans):
        tic[i] = s.total_intensity
        if len(s.channels) == 0:
            continue
        mz_vals = cal.to_mz(s.channels)

        for name, target_mz in species_db.items():
            # Halve tolerance for doubly-charged (isotope spacing ~0.5 Da)
            tol = mz_tol / 2 if "[2+]" in name else mz_tol
            mask = np.abs(mz_vals - target_mz) <= tol
            profiles[name][i] = float(s.intensities[mask].sum())

    data: dict[str, np.ndarray | list] = {"drift_bin": np.arange(n_drift), "total_tic": tic}
    data.update(profiles)
    return pd.DataFrame(data)


# =============================================================================
# Plotting helpers
# =============================================================================

def plot_spectrum(
    mz: np.ndarray,
    intensity: np.ndarray,
    title: str,
    annotations: list[dict] | None = None,
    mz_range: tuple[float, float] | None = None,
) -> "go.Figure":
    import plotly.graph_objects as go

    if mz_range:
        mask = (mz >= mz_range[0]) & (mz <= mz_range[1])
    else:
        mask = (mz > 0) & (mz < 5000) & np.isfinite(mz) & np.isfinite(intensity)
    mz_p, int_p = mz[mask], intensity[mask]

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=mz_p, y=int_p, mode="lines",
                             line=dict(color="#1f77b4", width=0.8), name="Spectrum"))

    if annotations:
        for ann in annotations[:25]:
            if ann.get("species", "unassigned") != "unassigned":
                fig.add_annotation(
                    x=ann["mono_mz"], y=ann["max_int"],
                    text=f"{ann['species']}<br>{ann['mono_mz']:.1f}",
                    showarrow=True, arrowhead=2, arrowsize=0.5,
                    font=dict(size=7), bgcolor="rgba(255,255,255,0.7)",
                )

    fig.update_layout(title=title, xaxis_title="m/z", yaxis_title="Intensity",
                      template="plotly_white", width=1200, height=500)
    return fig


def plot_drift_profiles(
    drift_df: pd.DataFrame,
    run_name: str,
    top_n: int = 10,
) -> "go.Figure":
    """Plot per-species drift profiles as overlaid traces."""
    import plotly.graph_objects as go

    drift_bins = drift_df["drift_bin"].values
    fig = go.Figure()

    # TIC trace (background, gray)
    tic = drift_df["total_tic"].values
    fig.add_trace(go.Scatter(
        x=drift_bins, y=tic, mode="lines", name="Total TIC",
        line=dict(color="lightgray", width=2), fill="tozeroy",
        fillcolor="rgba(200,200,200,0.2)",
    ))

    # Species traces — pick top N by total intensity across all drift bins
    species_cols = [c for c in drift_df.columns if c not in ("drift_bin", "total_tic")]
    totals = {c: drift_df[c].sum() for c in species_cols}
    top_species = sorted(totals, key=lambda c: -totals[c])[:top_n]

    for name in top_species:
        vals = drift_df[name].values
        if vals.max() > 0:
            fig.add_trace(go.Scatter(
                x=drift_bins, y=vals, mode="lines", name=name,
                line=dict(width=1.5),
            ))

    fig.update_layout(
        title=f"{run_name} — Per-Species Drift Time Profiles",
        xaxis_title="Drift Time Bin",
        yaxis_title="Intensity",
        template="plotly_white",
        width=1200, height=600,
        legend=dict(font=dict(size=9)),
    )
    return fig


def plot_2d_imms_scatter(
    scans: list[ScanData],
    cal: Calibration,
    run_name: str,
    mz_range: tuple[float, float] | None = None,
    max_points: int = 300_000,
) -> "go.Figure":
    """2D IM-MS scatter plot of raw data — no binning, no smoothing.

    Each point = one (m/z, drift_bin) observation colored by log intensity.
    Shows the full data range by default. Use mz_range to override.
    """
    import plotly.graph_objects as go

    # Collect all raw data points
    mz_all: list[float] = []
    drift_all: list[int] = []
    int_all: list[float] = []

    for drift_idx, s in enumerate(scans):
        if len(s.channels) == 0:
            continue
        mz = cal.to_mz(s.channels)
        valid = np.isfinite(mz) & (mz > 0)
        if mz_range:
            valid &= (mz >= mz_range[0]) & (mz <= mz_range[1])
        mz_v = mz[valid]
        int_v = s.intensities[valid]
        mz_all.extend(mz_v.tolist())
        drift_all.extend([drift_idx] * len(mz_v))
        int_all.extend(int_v.tolist())

    n_points = len(mz_all)
    mz_arr = np.array(mz_all)
    drift_arr = np.array(drift_all)
    int_arr = np.array(int_all)

    # Downsample if too many points for browser performance
    if n_points > max_points:
        idx = np.random.default_rng(42).choice(n_points, max_points, replace=False)
        mz_arr, drift_arr, int_arr = mz_arr[idx], drift_arr[idx], int_arr[idx]
        n_shown = max_points
    else:
        n_shown = n_points

    log_int = np.log10(int_arr + 1)

    fig = go.Figure(data=go.Scattergl(
        x=mz_arr, y=drift_arr,
        mode="markers",
        marker=dict(
            size=2,
            color=log_int,
            colorscale="Viridis",
            colorbar=dict(title="log₁₀(I+1)"),
            opacity=0.6,
        ),
        hovertemplate="m/z: %{x:.2f}<br>Drift bin: %{y}<br>Intensity: %{customdata:.0f}<extra></extra>",
        customdata=int_arr,
    ))

    mz_lo = float(mz_arr.min()) if len(mz_arr) > 0 else 0
    mz_hi = float(mz_arr.max()) if len(mz_arr) > 0 else 1
    fig.update_layout(
        title=(f"{run_name} — 2D IM-MS (raw data)<br>"
               f"<sub>{n_shown:,} points shown"
               f"{f' (sampled from {n_points:,})' if n_shown < n_points else ''}"
               f", m/z {mz_lo:.1f}–{mz_hi:.1f}</sub>"),
        xaxis_title="m/z",
        yaxis_title="Drift Time Bin",
        template="plotly_white",
        width=1200, height=700,
        yaxis=dict(range=[-1, len(scans)]),
    )
    return fig


# =============================================================================
# Classification helpers
# =============================================================================

def classify_run(name: str) -> str:
    u = name.upper()
    if "CALI" in u and "UBQ" in u:
        return "calibrant_ubq"
    if "UBQ" in u:
        return "protein_ubq"
    if "MYO" in u:
        return "protein_myo"
    if "CYTC" in u:
        return "protein_cytc"
    if "ADP" in u:
        if "AMAC" in u:
            return "adp_amac"
        if "WATER" in u:
            return "adp_water"
        if "EDDA" in u:
            return "adp_edda"
        return "adp_other"
    return "unknown"


def extract_buffer_info(name: str) -> tuple[str, float | None, str | None]:
    """Extract buffer type, pH, and concentration from filename."""
    u = name.upper()
    buf = "unknown"
    ph: float | None = None
    conc: str | None = None

    if "1M_EDDA" in u:
        buf, conc = "EDDA", "1M"
    elif "200MM_EDDA" in u:
        buf, conc = "EDDA", "200mM"
    elif "AMAC" in u:
        buf = "AMAC"
    elif "WATER" in u:
        buf = "water"

    ph_match = re.search(r"PH[_]?(\d+[-.]?\d*)", u)
    if ph_match:
        try:
            ph = float(ph_match.group(1).replace("-", "."))
        except ValueError:
            pass

    return buf, ph, conc


def auto_species_db(run_class: str) -> dict[str, float]:
    """Auto-select species DB from run classification."""
    if "adp" in run_class:
        return build_adp_species_db()
    if "ubq" in run_class:
        return build_protein_species_db("UBQ", 8565.8, 5, 15)
    if "myo" in run_class:
        return build_protein_species_db("MYO", 16951.5, 8, 25)
    if "cytc" in run_class:
        return build_protein_species_db("CYTC", 12384.0, 7, 20)
    return {}
