"""CCS calibration orchestrator for TW-IMS data.

CSV-driven: calibrant and analyte configs specify species parameters.
Data directory (where converted _ms.txt/_im.txt live) is a single argument.

Usage:
    python -m deconvovo.imms_ccs_calibrate \
        -o output/ccs \
        --data-dir output/converted \
        --calibrant-csv config/calibrants_tunemix.csv \
        --analyte-csv config/analytes_adp.csv
"""
from __future__ import annotations

import argparse
import json
import logging
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from deconvovo.io import read_im_txt, read_extern_inf, find_pusher_period, measure_ms_intensity
from deconvovo.ccs_peak_pick import (
    compute_auto_mz_window, extract_drift_profile,
    detect_peaks_in_profile, select_peak,
)
from deconvovo.ccs_fit import build_calibration_curve
from deconvovo.ccs_convert import (
    apply_ccs, convert_profile_to_ccs,
    resample_to_uniform_ccs, smooth_ccs_profile,
)
from deconvovo.ccs_plot import plot_calibration, plot_drift_profile, plot_ccs_single

_log = logging.getLogger("ccs_calibrate")


def _setup_logging(out_dir: Path) -> None:
    _log.setLevel(logging.DEBUG)
    _log.handlers.clear()
    fmt = logging.Formatter("%(levelname)s: %(message)s")
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    ch.setFormatter(fmt)
    _log.addHandler(ch)
    fh = logging.FileHandler(out_dir / "warnings.log", mode="w")
    fh.setLevel(logging.WARNING)
    fh.setFormatter(logging.Formatter("%(asctime)s %(levelname)s: %(message)s",
                                       datefmt="%Y-%m-%d %H:%M:%S"))
    _log.addHandler(fh)


def _resolve_mz_window(row, mw, z):
    """Resolve mz_window from CSV row: 'auto', numeric, or NaN → auto."""
    raw = row["mz_window"]
    formula = str(row.get("formula", "")) if "formula" in row.index else ""
    if str(raw).strip().lower() == "auto" or pd.isna(raw):
        return compute_auto_mz_window(mw, z, formula)
    return float(raw)


# =============================================================================
# Orchestrator
# =============================================================================

def _find_text_file(data_dirs: list[Path], filename: str) -> Path | None:
    """Search multiple data directories for a text file."""
    for d in data_dirs:
        p = d / filename
        if p.exists():
            return p
    return None


def run(out_dir: Path, data_dirs: Path | list[Path], calibrant_csv: Path,
        analyte_csv: Path | None = None,
        conversion_method: str = "direct") -> dict:
    """Run CCS calibration.

    Args:
        out_dir: Output directory for results.
        data_dirs: One or more directories containing converted _ms.txt/_im.txt files.
        calibrant_csv: CSV with calibrant species (name, mw, z, mz, ccs, raw_path, ...).
        analyte_csv: Optional CSV with analyte species.
        conversion_method: "direct" (power-law) or "twostep" (linear).
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    _setup_logging(out_dir)

    # Normalize to list
    if isinstance(data_dirs, Path):
        data_dirs = [data_dirs]
    data_dirs = [d for d in data_dirs if d.is_dir()]
    if not data_dirs:
        raise FileNotFoundError("No valid data directories found.\n"
                                "Convert .raw files first, then point --data-dir at the output.")

    cal_df = pd.read_csv(calibrant_csv)
    _log.info("=== CCS Calibration ===")
    _log.info("  Calibrant CSV: %s (%d rows)", calibrant_csv, len(cal_df))
    _log.info("  Data dirs: %s", [str(d) for d in data_dirs])
    _log.info("  Conversion: %s", conversion_method)

    # Defaults for optional columns
    defaults = {"mz_window": "auto", "peak_select": "highest",
                "min_ms_intensity": 100.0, "min_im_intensity": 500.0,
                "peak_height_frac": 0.10,
                "gas_mw": 28.014}
    for col, default in defaults.items():
        if col not in cal_df.columns:
            cal_df[col] = default
        else:
            cal_df[col] = cal_df[col].fillna(default)

    # Get pusher periods + EDC from .raw dirs
    pusher_cache = {}
    edc_coeff = None
    for raw_path in cal_df["raw_path"].unique():
        rp = Path(raw_path)
        if rp.is_dir():
            pp = find_pusher_period(rp)
            if pp is not None:
                pusher_cache[raw_path] = pp
            params = read_extern_inf(rp)
            if "EDC Delay Coefficient" in params:
                edc_coeff = float(params["EDC Delay Coefficient"])
    if edc_coeff is None:
        raise ValueError("Cannot read EDC Delay Coefficient from any .raw directory. "
                         "Check raw_path column in calibrant CSV.")
    _log.info("  EDC: %.4f", edc_coeff)

    # --- Pick calibrant peaks ---
    _log.info("--- Picking calibrant peaks ---")
    cali_csv_dir = out_dir / "cali_csv"
    cali_png_dir = out_dir / "cali_png"
    cali_csv_dir.mkdir(exist_ok=True)
    cali_png_dir.mkdir(exist_ok=True)

    all_peaks = []
    for _, row in cal_df.iterrows():
        name = row["name"]
        mw = float(row["mw"]); z = int(row["z"]); mz = float(row["mz"])
        ccs = float(row["ccs"])
        raw_path = row["raw_path"]
        run = Path(raw_path).stem
        mz_window = _resolve_mz_window(row, mw, z)
        peak_sel = str(row["peak_select"])
        min_ms = float(row["min_ms_intensity"])
        min_im = float(row["min_im_intensity"])
        phf = float(row["peak_height_frac"])

        pp = pusher_cache.get(raw_path)
        if pp is None:
            _log.warning("No pusher for %s, skipping", name); continue

        ms_file = _find_text_file(data_dirs, f"{run}_ms.txt")
        im_file = _find_text_file(data_dirs, f"{run}_im.txt")

        if ms_file is None:
            _log.warning("  %s: %s_ms.txt not found in any data dir, skipping", name, run); continue

        ms_int = measure_ms_intensity(ms_file, mz, mz_window)
        if ms_int < min_ms:
            _log.info("  %s: MS int %.0f < %.0f, skip", name, ms_int, min_ms); continue

        im_data = read_im_txt(im_file) if im_file else pd.DataFrame()
        if im_data.empty:
            _log.warning("  %s: no IM data", name); continue

        profile = extract_drift_profile(im_data, mz, mz_window)
        detected = detect_peaks_in_profile(profile, pp, phf)
        if not detected:
            _log.warning("  %s: no peaks", name); continue

        detected = [p for p in detected if p["peak_intensity"] >= min_im]
        if not detected:
            _log.info("  %s: all peaks below IM threshold", name); continue

        selected = select_peak(detected, peak_sel)
        _log.info("  %s: t_D=%.3f ms, CCS_lit=%.1f, %d peaks (%s), ±%.1f Da",
                  name, selected["drift_time_ms"], ccs, len(detected), peak_sel, mz_window)

        all_peaks.append({
            "name": name, "mw": mw, "z": z, "mz": mz, "ccs": ccs,
            "drift_bin": selected["drift_bin"],
            "drift_time_ms": selected["drift_time_ms"],
            "peak_intensity": selected["peak_intensity"],
            "ms_intensity": ms_int,
        })

        safe = name.replace("/", "_").replace(" ", "_")
        pk_rows = [{"name": name, "z": z, "mz": mz, "ccs": ccs,
                     "peak_bin": pk["bin_index"], "peak_drift_ms": pk["drift_time_ms"],
                     "peak_intensity": pk["peak_intensity"],
                     "selected": pk["bin_index"] == selected["bin_index"]}
                    for pk in detected]
        pd.DataFrame(pk_rows).to_csv(cali_csv_dir / f"{run}_{safe}.csv",
                                      index=False, float_format="%.4f")
        plot_drift_profile(profile, pp, detected, selected, name, z, mz,
                           cali_png_dir / f"{run}_{safe}.png", mw, mz_window)

    if len(all_peaks) < 3:
        raise ValueError(f"Need >=3 calibrant peaks, got {len(all_peaks)}")

    # --- Fit calibration ---
    _log.info("--- Building calibration curve ---")
    gas_mw_vals = cal_df["gas_mw"].unique()
    if len(gas_mw_vals) > 1:
        raise ValueError(f"All calibrant rows must use the same gas_mw, got: {gas_mw_vals}")
    gas_mw = float(gas_mw_vals[0])
    cal = build_calibration_curve(all_peaks, edc_coeff, gas_mw)

    with open(out_dir / "calibration_curve.json", "w") as f:
        json.dump({k: v for k, v in cal.items() if k != "calibrant_summary"}, f, indent=2)
    summary = cal["calibrant_summary"]
    pd.DataFrame(summary).to_csv(out_dir / "calibrant_summary.csv",
                                  index=False, float_format="%.4f")
    plot_calibration(cal, out_dir / "calibration_plots.png", summary)

    # --- Apply to analytes ---
    if analyte_csv and analyte_csv.exists():
        _log.info("--- Applying CCS to analytes ---")
        ana_df = pd.read_csv(analyte_csv)
        defaults_ana = {"mz_window": "auto", "min_ms_intensity": 50.0}
        for col, default in defaults_ana.items():
            if col not in ana_df.columns:
                ana_df[col] = default
            else:
                ana_df[col] = ana_df[col].fillna(default)

        for d in ["csv_raw", "csv_smoothed", "png_raw", "png_smoothed", "png_overlay"]:
            (out_dir / d).mkdir(exist_ok=True)

        ana_summary = []

        # Group by raw_path to read each IM file only once
        for raw_path_str, group in ana_df.groupby("raw_path"):
            run = Path(raw_path_str).stem

            pp = pusher_cache.get(raw_path_str)
            if pp is None:
                pp = find_pusher_period(Path(raw_path_str)) if Path(raw_path_str).is_dir() else None
            if pp is None:
                _log.warning("No pusher for run %s, skipping", run)
                continue

            im_file = _find_text_file(data_dirs, f"{run}_im.txt")
            ms_file = _find_text_file(data_dirs, f"{run}_ms.txt")
            im_data = read_im_txt(im_file) if im_file else pd.DataFrame()
            _log.info("  %s: %d species", run, len(group))

            for _, row in group.iterrows():
                sp = row["species"]; mw = float(row["mw"]); z = int(row["z"])
                mz = float(row["mz"])
                mz_w = _resolve_mz_window(row, mw, z)
                min_ms = float(row["min_ms_intensity"])

                if ms_file is None:
                    continue
                ms_int = measure_ms_intensity(ms_file, mz, mz_w)
                if ms_int < min_ms:
                    ana_summary.append({"run": run, "species": sp, "mw": mw, "z": z,
                                        "mz": mz, "ms_intensity": ms_int,
                                        "detected": ms_int > 0, "above_threshold": False,
                                        "peak_ccs": None, "peak_drift_bin": None,
                                        "peak_drift_time_ms": None,
                                        "t_prime": None, "t_double_prime": None})
                    continue

                if im_data.empty:
                    continue

                profile = extract_drift_profile(im_data, mz, mz_w)
                if profile.max() == 0:
                    ana_summary.append({"run": run, "species": sp, "mw": mw, "z": z,
                                        "mz": mz, "ms_intensity": ms_int,
                                        "detected": True, "above_threshold": True,
                                        "peak_ccs": None, "peak_drift_bin": None,
                                        "peak_drift_time_ms": None,
                                        "t_prime": None, "t_double_prime": None})
                    continue

                bins, dt_ms, ccs_n, int_n = convert_profile_to_ccs(
                    profile, pp, mz, z, mw, cal, conversion_method)
                if len(ccs_n) == 0:
                    continue

                ccs_uni, int_uni = resample_to_uniform_ccs(ccs_n, int_n)
                int_sm = smooth_ccs_profile(ccs_uni, int_uni)

                peak_ccs = ccs_uni[np.argmax(int_sm)] if int_sm.max() > 0 else None
                pi = np.argmax(int_n)
                peak_dt = float(dt_ms[pi])
                tp = peak_dt - cal["edc_coeff"] * math.sqrt(mz) / 1000.0
                mu = math.sqrt(1.0 / mw + 1.0 / cal["gas_mw"])
                tdp = (tp ** cal["X"]) * z * mu if tp > 0 else None

                ana_summary.append({"run": run, "species": sp, "mw": mw, "z": z,
                                    "mz": mz, "ms_intensity": ms_int,
                                    "detected": True, "above_threshold": True,
                                    "peak_ccs": peak_ccs, "peak_drift_bin": int(bins[pi]),
                                    "peak_drift_time_ms": peak_dt,
                                    "t_prime": tp if tp > 0 else None,
                                    "t_double_prime": tdp})

                safe = sp.replace("/", "_").replace(" ", "_")
                stem = f"{run}_{safe}_ccs"

                pd.DataFrame({"species": sp, "mw": mw, "z": z, "mz": mz,
                               "drift_bin": bins.astype(int), "drift_time_ms": dt_ms,
                               "ccs": ccs_n, "intensity": int_n}
                             ).to_csv(out_dir / "csv_raw" / f"{stem}.csv",
                                      index=False, float_format="%.4f")
                pd.DataFrame({"species": sp, "mw": mw, "z": z, "mz": mz,
                               "ccs": ccs_uni, "intensity": int_sm}
                             ).to_csv(out_dir / "csv_smoothed" / f"{stem}.csv",
                                      index=False, float_format="%.4f")
                for style, cr, ir, cs, ism, d in [
                    ("Raw", ccs_n, int_n, None, None, "png_raw"),
                    ("Smoothed", None, None, ccs_uni, int_sm, "png_smoothed"),
                    ("Overlay", ccs_n, int_n, ccs_uni, int_sm, "png_overlay"),
                ]:
                    plot_ccs_single(cr, ir, cs, ism, sp, run, style,
                                    out_dir / d / f"{stem}.png", mw, z, mz, mz_w)

        if ana_summary:
            pd.DataFrame(ana_summary).to_csv(out_dir / "analyte_summary.csv",
                                              index=False, float_format="%.4f")
            n_above = sum(1 for r in ana_summary if r["above_threshold"])
            _log.info("Wrote analyte_summary.csv (%d rows, %d above threshold)",
                      len(ana_summary), n_above)

    _log.info("")
    _log.info("=== CCS Calibration complete ===")
    _log.info("  R²: %.5f, Points: %d",
              cal["r2_lnln"], cal["n_points"])
    return cal


def main() -> None:
    parser = argparse.ArgumentParser(description="CCS calibration for TW-IMS data.")
    parser.add_argument("-o", "--output", required=True, help="Output directory")
    parser.add_argument("--data-dir", required=True, nargs="+",
                        help="Directory(s) with converted _ms.txt/_im.txt files")
    parser.add_argument("--calibrant-csv", required=True, help="Calibrant CSV")
    parser.add_argument("--analyte-csv", default=None, help="Analyte CSV (optional)")
    parser.add_argument("--conversion-method", default="direct",
                        choices=["direct", "twostep"])
    args = parser.parse_args()
    run(Path(args.output).resolve(),
        [Path(d).resolve() for d in args.data_dir],
        Path(args.calibrant_csv).resolve(),
        Path(args.analyte_csv).resolve() if args.analyte_csv else None,
        args.conversion_method)


if __name__ == "__main__":
    main()
