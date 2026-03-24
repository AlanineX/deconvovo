"""CCS calibration — power-law fitting (honest, no outlier rejection)."""
from __future__ import annotations

import logging
import math

import numpy as np
from scipy import stats
from deconvovo.io import GAS_MW_N2

_log = logging.getLogger("ccs_calibrate")


def _fit_lnln(points: list[dict]) -> tuple[float, float, float]:
    """Fit ln(Ω') = X·ln(t'_D) + ln(A). Returns (X, A, R²)."""
    ln_t = np.array([math.log(p["t_prime"]) for p in points])
    ln_o = np.array([math.log(p["omega_prime"]) for p in points])
    sl, ic, r, _, _ = stats.linregress(ln_t, ln_o)
    return sl, math.exp(ic), r ** 2


def build_calibration_curve(peaks: list[dict], edc_coeff: float,
                            gas_mw: float = GAS_MW_N2) -> dict:
    """Build TW-IMS CCS calibration curve.

    Fits ALL provided points honestly — no outlier rejection. The user
    controls which points to include via the calibrant CSV.
    """
    if len(peaks) < 3:
        raise ValueError(f"Need >=3 calibrant peaks, got {len(peaks)}")

    # Prepare: compute t'_D and Omega'
    valid = []
    for p in peaks:
        tp = p["drift_time_ms"] - edc_coeff * math.sqrt(p["mz"]) / 1000.0
        if tp <= 0:
            _log.warning("t'_D <= 0 for %s, excluding", p["name"])
            continue
        mu = math.sqrt(1.0 / p["mw"] + 1.0 / gas_mw)
        op = p["ccs"] / (p["z"] * mu)
        valid.append({**p, "t_prime": tp, "omega_prime": op, "mu_factor": mu})

    if len(valid) < 3:
        raise ValueError(f"Need >=3 valid points after EDC correction, got {len(valid)}")

    # Honest fit — all points, no rejection
    X, A, r2 = _fit_lnln(valid)
    _log.info("Fit 1 (ln-ln): X=%.4f, A=%.4f, R²=%.6f (%d pts)", X, A, r2, len(valid))

    # Fit 2 (diagnostic)
    tdp_arr = np.array([(p["t_prime"] ** X) * p["z"] * p["mu_factor"] for p in valid])
    ccs_arr = np.array([p["ccs"] for p in valid])
    slope_d, intercept_d, r_lin, _, _ = stats.linregress(tdp_arr, ccs_arr)
    r2_lin = r_lin ** 2
    _log.info("Fit 2 (diagnostic): slope=%.4f, intercept=%.4f, R²=%.6f",
              slope_d, intercept_d, r2_lin)

    # Build summary
    summary = []
    for vp in valid:
        td = (vp["t_prime"] ** X) * vp["z"] * vp["mu_factor"]
        ccs_direct = A * (vp["t_prime"] ** X) * vp["z"] * vp["mu_factor"]
        ccs_twostep = slope_d * td + intercept_d
        summary.append({
            "name": vp["name"], "mw": vp["mw"], "z": vp["z"], "mz": vp["mz"],
            "ms_intensity": vp.get("ms_intensity", 0),
            "peak_intensity_im": vp.get("peak_intensity", 0),
            "drift_time_ms": vp["drift_time_ms"], "t_prime": vp["t_prime"],
            "t_double_prime": float(td),
            "CCS_literature": vp["ccs"],
            "CCS_direct": float(ccs_direct), "CCS_twostep": float(ccs_twostep),
            "err_direct_pct": (ccs_direct - vp["ccs"]) / vp["ccs"] * 100,
            "err_twostep_pct": (ccs_twostep - vp["ccs"]) / vp["ccs"] * 100,
        })

    return {
        "X": X, "A": A,
        "slope_diagnostic": slope_d, "intercept_diagnostic": intercept_d,
        "r2_lnln": r2, "r2_linear_diagnostic": r2_lin,
        "edc_coeff": edc_coeff, "gas_mw": gas_mw,
        "n_points": len(valid),
        "calibrant_summary": summary,
    }
