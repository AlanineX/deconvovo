"""CCS calibration — all plotting functions."""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from deconvovo.ccs_convert import significant_range


_FS = 16  # base font size — everything else scales from this

def _get_plt():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Liberation Sans", "Arial", "Helvetica"],
        "font.size": _FS,
        "axes.titlesize": _FS * 1.15,
        "axes.labelsize": _FS,
        "xtick.labelsize": _FS * 0.9,
        "ytick.labelsize": _FS * 0.9,
        "legend.fontsize": _FS * 0.85,
    })
    return plt


def plot_calibration(cal, out_path, summary):
    """2-panel calibration plot: ln-ln + linear diagnostic."""
    plt = _get_plt()
    pts = list(summary)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    X, A = cal["X"], cal["A"]

    lt = [math.log(p["t_prime"]) for p in pts]
    lo = [math.log(p["CCS_literature"] / (p["z"] * math.sqrt(1/p["mw"] + 1/cal["gas_mw"]))) for p in pts]
    ax1.scatter(lt, lo, c="steelblue", s=80, zorder=5, edgecolors="black",
                linewidth=0.5, label="Calibrant")
    for p, x, y in zip(pts, lt, lo):
        ax1.annotate(p['name'], (x, y), fontsize=_FS * 0.55,
                     textcoords="offset points", xytext=(5, 5), alpha=0.8)
    t_range = np.linspace(min(lt) - 0.15, max(lt) + 0.15, 100)
    ax1.plot(t_range, X * t_range + math.log(A), "r-", linewidth=1.5, alpha=0.8,
             label=f"y = {X:.3f}x + {math.log(A):.3f}")
    ax1.set_xlabel("ln(t'_D)"); ax1.set_ylabel("ln(Ω')")
    ax1.set_title(f"Fit 1: ln-ln (R² = {cal['r2_lnln']:.5f})")
    ax1.legend()

    td = [p["t_double_prime"] for p in pts]
    cc = [p["CCS_literature"] for p in pts]
    ax2.scatter(td, cc, c="steelblue", s=80, zorder=5, edgecolors="black",
                linewidth=0.5, label="Calibrant")
    for p, x, y in zip(pts, td, cc):
        ax2.annotate(p['name'], (x, y), fontsize=_FS * 0.55,
                     textcoords="offset points", xytext=(5, 5), alpha=0.8)
    td_range = np.linspace(min(td) * 0.9, max(td) * 1.1, 100)
    ax2.plot(td_range, cal["slope_diagnostic"] * td_range + cal["intercept_diagnostic"],
             "r-", linewidth=1.5, alpha=0.8,
             label=f"y = {cal['slope_diagnostic']:.1f}x + {cal['intercept_diagnostic']:.1f}")
    ax2.set_xlabel("t''_D"); ax2.set_ylabel("CCS (Å²)")
    ax2.set_title(f"Fit 2: diagnostic (R² = {cal['r2_linear_diagnostic']:.5f})")
    ax2.legend()

    fig.suptitle(f"TW-IMS CCS Calibration ({cal['n_points']} points)",
                 fontsize=_FS * 1.3, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_drift_profile(profile, pusher_us, all_peaks, selected, name, z, mz,
                       out_path, mw=0, mz_window=0):
    """Calibrant drift profile with all peaks labeled, selected in red."""
    plt = _get_plt()
    fig, ax = plt.subplots(figsize=(10, 4))
    dt = np.arange(len(profile)) * pusher_us / 1000.0
    ax.fill_between(dt, profile, alpha=0.25, color="gray")
    ax.plot(dt, profile, color="gray", linewidth=0.8)
    for pk in all_peaks:
        b = pk["bin_index"]
        is_sel = pk["bin_index"] == selected["bin_index"]
        c = "red" if is_sel else "steelblue"
        s = 80 if is_sel else 40
        ax.scatter(dt[b], profile[b], c=c, s=s, zorder=5,
                   edgecolors="black" if is_sel else c, linewidth=1.0)
        label = f"{pk['drift_time_ms']:.2f} ms"
        if is_sel:
            label += " [SEL]"
        ax.annotate(label, (dt[b], profile[b]),
                    fontsize=_FS * (0.65 if is_sel else 0.55),
                    fontweight="bold" if is_sel else "normal", color=c,
                    textcoords="offset points", xytext=(5, 8), alpha=0.9)
    ax.set_xlabel("Drift Time (ms)"); ax.set_ylabel("Intensity")
    ax.set_title(f"{name}  z={z}+  |  m/z {mz:.2f} ± {mz_window:.1f} Da",
                 fontsize=_FS * 0.9)
    nz = np.nonzero(profile)[0]
    if len(nz) > 0:
        ax.set_xlim(max(0, nz[0] - 5) * pusher_us / 1000,
                     min(len(profile) - 1, nz[-1] + 5) * pusher_us / 1000)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_ccs_single(ccs_raw, int_raw, ccs_sm, int_sm, species, run, style,
                    out_path, mw=0, z=0, mz=0, mz_window=0):
    """CCS profile: Raw, Smoothed, or Overlay."""
    plt = _get_plt()
    fig, ax = plt.subplots(figsize=(8, 4))
    has_r = ccs_raw is not None and len(ccs_raw) > 0
    has_s = ccs_sm is not None and len(ccs_sm) > 0
    if has_r:
        rm = int_raw.max() if int_raw.max() > 0 else 1.0
        if style == "Overlay":
            ax.fill_between(ccs_raw, int_raw / rm, alpha=0.2, color="gray", label="Raw")
            ax.plot(ccs_raw, int_raw / rm, color="gray", linewidth=0.6, alpha=0.6)
        else:
            ax.fill_between(ccs_raw, int_raw / rm, alpha=0.3, color="gray")
            ax.plot(ccs_raw, int_raw / rm, color="gray", linewidth=1.0)
    if has_s:
        sm = int_sm.max() if int_sm.max() > 0 else 1.0
        ax.plot(ccs_sm, int_sm / sm, color="steelblue", linewidth=1.5,
                label="Smoothed" if style == "Overlay" else None)
    ref_c = ccs_sm if has_s else ccs_raw
    ref_i = int_sm if has_s else int_raw
    lo, hi = significant_range(ref_c, ref_i)
    ax.set_xlim(lo, hi)
    ax.set_xlabel("CCS (Å²)"); ax.set_ylabel("Normalized Intensity")
    ax.set_title(f"{species}  z={z}  |  {run}  ({style})",
                 fontsize=_FS * 0.9)
    if style == "Overlay":
        ax.legend()
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
