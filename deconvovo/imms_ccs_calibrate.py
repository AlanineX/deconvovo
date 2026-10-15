"""CCS calibration orchestrator for TW-IMS data.

CSV-driven: calibrant and analyte configs specify species parameters.
Data directory (where converted _ms.txt/_im.txt live) is a single argument.

Usage:
    python -m deconvovo.imms_ccs_calibrate \
        -o output/ccs \
        --data-dir output/converted \
        --calibrant-csv demo_data/production/protein_adp/calibrants.csv \
        --analyte-csv demo_data/production/protein_adp/analytes.csv
"""
from __future__ import annotations

import argparse
import json
import logging
import math
import re
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from deconvovo.io import read_im_txt, read_extern_inf, find_pusher_period, measure_ms_intensity
from deconvovo.paths import local_path
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
    # Keep the GUI's log bridge while replacing our own output handlers.
    for handler in list(_log.handlers):
        if getattr(handler, "_waters_output_handler", False):
            _log.removeHandler(handler)
            handler.close()
    fmt = logging.Formatter("%(levelname)s: %(message)s")
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    ch.setFormatter(fmt)
    ch._waters_output_handler = True
    _log.addHandler(ch)
    fh = logging.FileHandler(out_dir / "warnings.log", mode="w", encoding="utf-8")
    fh.setLevel(logging.WARNING)
    fh.setFormatter(logging.Formatter("%(asctime)s %(levelname)s: %(message)s",
                                       datefmt="%Y-%m-%d %H:%M:%S"))
    fh._waters_output_handler = True
    _log.addHandler(fh)


def normalize_calibrant_columns(frame: pd.DataFrame) -> pd.DataFrame:
    """Accept the reported reference CCS header while retaining the fit's schema."""
    frame = frame.copy()
    if "CCS_reported" not in frame.columns:
        return frame
    if "ccs" in frame.columns:
        reported = pd.to_numeric(frame["CCS_reported"], errors="coerce")
        legacy = pd.to_numeric(frame["ccs"], errors="coerce")
        if not (reported.eq(legacy) | (reported.isna() & legacy.isna())).all():
            raise ValueError("Conflicting CCS_reported and ccs columns")
        return frame.drop(columns="CCS_reported")
    return frame.rename(columns={"CCS_reported": "ccs"})


def normalize_calibrant_mz(frame: pd.DataFrame) -> pd.DataFrame:
    """Derive the legacy centre/half-width from explicit extraction bounds."""
    frame = frame.copy()
    if not {"mz_min", "mz_max"} & set(frame.columns):
        return frame
    for index, row in frame.iterrows():
        lower, upper = row.get("mz_min"), row.get("mz_max")
        has_lower, has_upper = _has_value(lower), _has_value(upper)
        if not has_lower and not has_upper:
            continue
        if not has_lower or not has_upper:
            raise ValueError(f"Row {index + 1}: enter both mz_min and mz_max")
        lower, upper = _validate_mz_bounds(lower, upper, f"Row {index + 1}")
        frame.loc[index, "mz"] = (lower + upper) / 2
        frame.loc[index, "mz_window"] = (upper - lower) / 2
    return frame


def _has_value(value) -> bool:
    return value is not None and not pd.isna(value) and bool(str(value).strip())


def _validate_mz_bounds(lower, upper, label="m/z"):
    try:
        lower, upper = float(lower), float(upper)
    except (TypeError, ValueError):
        raise ValueError(f"{label}: m/z bounds must be finite, positive and increasing") from None
    if not math.isfinite(lower) or not math.isfinite(upper) or lower <= 0 or upper <= lower:
        raise ValueError(f"{label}: m/z bounds must be finite, positive and increasing")
    return lower, upper


def resolve_mz_bounds(row, mw, z):
    """Return explicit bounds, or resolve the legacy centre and half-width."""
    lower, upper = row.get("mz_min"), row.get("mz_max")
    has_lower, has_upper = _has_value(lower), _has_value(upper)
    if has_lower != has_upper:
        raise ValueError("Enter both mz_min and mz_max")
    if has_lower:
        return _validate_mz_bounds(lower, upper)
    width = _resolve_mz_window(row, mw, z)
    centre = float(row["mz"])
    return _validate_mz_bounds(centre - width, centre + width)


def resolve_calibrant_peak(detected, row, *, allow_unselected_manual=False):
    """Select a detected native peak; Manual never substitutes an automatic peak."""
    mode = str(row.get("peak_select", "highest")).strip().lower()
    mode = {"latest": "last", "tallest": "highest"}.get(mode, mode)
    if mode not in {"highest", "last", "manual"}:
        raise ValueError("peak_select must be highest, last or manual")
    override = row.get("peak_drift_bin")
    if _has_value(override):
        try:
            peak_bin = float(override)
        except (TypeError, ValueError):
            peak_bin = float("nan")
        candidate = next((peak for peak in detected if peak["bin_index"] == peak_bin), None)
        if not math.isfinite(peak_bin) or not peak_bin.is_integer() or candidate is None:
            raise ValueError(f"{row.get('name', 'Reference')}: peak_drift_bin must be one of the detected candidates "
                             f"{[peak['bin_index'] for peak in detected]}")
        return candidate
    if mode == "manual":
        if allow_unselected_manual:
            return None
        raise ValueError(f"{row.get('name', 'Reference')}: Manual selection needs peak_drift_bin; "
                         "process the reference and select a detected peak first")
    return select_peak(detected, mode) if detected else None


def _read_species_csv(csv_path: Path) -> pd.DataFrame:
    """Resolve portable acquisition paths beside the table, retaining old cwd paths."""
    csv_path = local_path(csv_path).resolve()
    frame = normalize_calibrant_mz(normalize_calibrant_columns(pd.read_csv(csv_path)))
    if "raw_path" not in frame.columns:
        raise ValueError(f"Missing raw_path column in {csv_path.name}")

    return _resolve_raw_paths(frame, csv_path.parent)


def _resolve_raw_paths(frame, base_dir=None, default_raw=None):
    frame = frame.copy()
    if "raw_path" not in frame.columns:
        frame["raw_path"] = default_raw

    def resolve(value):
        if pd.isna(value) or not str(value).strip():
            value = default_raw
        if value is None or not str(value).strip():
            raise ValueError("Empty raw_path; every row needs a RAW acquisition path")
        path = local_path(str(value))
        if path.is_absolute():
            return str(path)
        relative = (base_dir or Path.cwd()) / path
        if relative.is_dir() or not path.is_dir():
            return str(relative.resolve())
        return str(path.resolve())

    frame["raw_path"] = frame["raw_path"].map(resolve)
    return frame


def _validate_species(frame, *, calibrants=False):
    required = {"mw", "z", "mz"} | ({"name", "ccs"} if calibrants else {"species"})
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"Missing columns: {', '.join(sorted(missing))}")
    if frame.empty:
        raise ValueError("Add at least one target" if not calibrants else "Add calibrant references")
    for index, row in frame.iterrows():
        label = row["name" if calibrants else "species"]
        if pd.isna(label) or not str(label).strip():
            raise ValueError(f"Row {index + 1}: enter a name")
        for column in ("mw", "z", "mz", "ccs") if calibrants else ("mw", "z", "mz"):
            try:
                value = float(row[column])
            except (TypeError, ValueError):
                raise ValueError(f"Row {index + 1}: {column} must be a finite positive number") from None
            if not math.isfinite(value) or value <= 0:
                raise ValueError(f"Row {index + 1}: {column} must be a finite positive number")
            if column == "z" and not value.is_integer():
                raise ValueError(f"Row {index + 1}: z must be a positive integer")
        width = _resolve_mz_window(row, float(row["mw"]), int(float(row["z"])))
        if not math.isfinite(width) or width <= 0:
            raise ValueError(f"Row {index + 1}: mz_window must be finite and positive")
        for column in ("min_ms_intensity", "min_im_intensity", "peak_height_frac", "gas_mw"):
            if column not in frame.columns:
                continue
            value = float(row[column])
            if not math.isfinite(value) or value < 0 or (column == "gas_mw" and value == 0):
                raise ValueError(f"Row {index + 1}: invalid {column}")
            if column == "peak_height_frac" and not 0 < value <= 1:
                raise ValueError(f"Row {index + 1}: peak_height_frac must be greater than 0 and at most 1")
        if calibrants:
            resolve_mz_bounds(row, float(row["mw"]), int(float(row["z"])))
            mode = str(row["peak_select"]).strip().lower()
            if mode not in {"highest", "tallest", "last", "latest", "manual"}:
                raise ValueError(f"Row {index + 1}: peak_select must be highest, last or manual")
            if mode == "manual" and not _has_value(row.get("peak_drift_bin")):
                raise ValueError(f"Row {index + 1}: Manual selection needs peak_drift_bin; "
                                 "process the reference and select a detected peak first")
    for column in ("mw", "z", "mz", "ccs") if calibrants else ("mw", "z", "mz"):
        frame[column] = pd.to_numeric(frame[column], errors="raise")
    frame["z"] = frame["z"].astype(int)


_CONDITION_FIELDS = (
    "EDC Delay Coefficient", "IMS Wave Velocity (m/s)", "IMS Wave Height (V)",
    "IMS Gas Flow (mL/min)", "HeliumCellGasFlow", "IMS",
    "Manual Trap Collision Energy", "Trap Collision Energy",
    "Using Auto Trap Collision Energy (eV)", "Transfer Collision Energy",
    "Variable Wave Height Enabled", "Wave Height Start (V)", "Wave Height End (V)",
    "Variable Wave Velocity Enabled", "Wave Velocity Start (m/s)", "Wave Velocity End (m/s)",
)
_REQUIRED_CONDITIONS = (*_CONDITION_FIELDS[:6], "Acquisition function")


def _acquisition_metadata(raw_path, data_dirs=()):
    raw_path = local_path(raw_path).resolve()
    params = read_extern_inf(raw_path) if raw_path.is_dir() else {}
    conditions = {key: params[key] for key in _CONDITION_FIELDS if key in params}
    conversion = {}
    for directory in data_dirs:
        settings = directory / "conversion_settings.json"
        if not settings.is_file():
            continue
        try:
            recorded = json.loads(settings.read_text(encoding="utf-8")).get(raw_path.stem, {})
            conversion = {key: recorded[key] for key in ("function", "ms_bin", "im_bin", "ms_smoothing") if key in recorded}
        except (OSError, ValueError, TypeError, AttributeError):
            continue
        if conversion:
            if "function" in conversion:
                conditions["Acquisition function"] = conversion["function"]
            break
    return {"raw_path": str(raw_path), "pusher_period_us": find_pusher_period(raw_path),
            "conditions": conditions, "conversion_settings": conversion}


def _review_warnings(cal):
    warnings = []
    if not cal.get("accepted", False):
        warnings.append("Calibration has not been marked reviewed; inspect reference peaks and fit before use.")
    if cal.get("r2_lnln", 0) <= 0.95:
        warnings.append("Calibration fit requires review (R² at or below 0.95).")
    if not cal.get("calibrant_summary"):
        warnings.append("Saved calibration has no reference summary; its reference range cannot be verified.")
    return warnings


def load_saved_calibration(model_or_path) -> dict:
    """Read saved coefficients without extracting standards or fitting again."""
    if isinstance(model_or_path, dict):
        cal = dict(model_or_path)
    else:
        cal = json.loads(local_path(model_or_path).read_text(encoding="utf-8"))
    if not isinstance(cal, dict):
        raise ValueError("Calibration JSON must contain an object")
    if "accepted" in cal and not isinstance(cal["accepted"], bool):
        raise ValueError("Saved calibration review state must be true or false")
    method = cal.get("method", "direct")
    if method not in {"direct", "twostep"}:
        raise ValueError("CCS method must be direct or twostep")
    for key in ("X", "A", "edc_coeff", "gas_mw"):
        try:
            value = float(cal[key])
        except (KeyError, TypeError, ValueError):
            raise ValueError(f"Saved calibration needs a finite {key}") from None
        if not math.isfinite(value) or (key in {"X", "A", "gas_mw"} and value <= 0):
            raise ValueError(f"Saved calibration has invalid {key}")
        cal[key] = value
    for key in ("slope_diagnostic", "intercept_diagnostic", "r2_lnln", "r2_linear_diagnostic"):
        if key in cal:
            value = float(cal[key])
            if not math.isfinite(value):
                raise ValueError(f"Saved calibration has invalid {key}")
            cal[key] = value
    if method == "twostep" and any(key not in cal for key in ("slope_diagnostic", "intercept_diagnostic")):
        raise ValueError("Two-step calibration needs its diagnostic slope and intercept")
    if "n_points" in cal and (float(cal["n_points"]) < 3 or not float(cal["n_points"]).is_integer()):
        raise ValueError("Saved calibration needs at least three reference points")
    cal.setdefault("method", method)
    summary = cal.get("calibrant_summary", [])
    if summary:
        if not isinstance(summary, list):
            raise ValueError("Saved calibration reference summary must contain rows")
        for row in summary:
            if not isinstance(row, dict) or not str(row.get("name", "")).strip():
                raise ValueError("Saved calibration reference rows need a name")
            for key in ("mw", "z", "mz", "CCS_literature", "t_prime", "t_double_prime"):
                try:
                    value = float(row[key])
                except (KeyError, TypeError, ValueError):
                    raise ValueError(f"Saved calibration reference row needs a finite positive {key}") from None
                if not math.isfinite(value) or value <= 0 or (key == "z" and not value.is_integer()):
                    raise ValueError(f"Saved calibration reference row has invalid {key}")
        times = [float(row["t_prime"]) for row in summary]
        cal["reference_t_prime_range"] = [min(times), max(times)]
    elif "reference_t_prime_range" in cal:
        bounds = [float(value) for value in cal["reference_t_prime_range"]]
        if len(bounds) != 2 or not all(math.isfinite(value) and value > 0 for value in bounds) or bounds[0] > bounds[1]:
            raise ValueError("Saved calibration has invalid reference range")
        cal["reference_t_prime_range"] = bounds
    cal["review_warnings"] = _review_warnings(cal)
    return cal


def save_calibration(calibration, destination, *, accepted=True) -> Path:
    """Save coefficients and their review context in one portable JSON file."""
    cal = load_saved_calibration(calibration)
    cal["accepted"] = bool(accepted)
    cal["review_warnings"] = _review_warnings(cal)
    destination = local_path(destination).resolve()
    if destination.suffix.lower() != ".json":
        raise ValueError("Save calibration as a .json file")
    if any(parent.suffix.lower() == ".raw" and parent.is_dir() for parent in destination.parents):
        raise ValueError("Save calibration outside the RAW acquisition")
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(cal, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    return destination


def _same_condition(left, right):
    try:
        return math.isclose(float(left), float(right), rel_tol=1e-6, abs_tol=1e-6)
    except (TypeError, ValueError):
        return str(left).strip().casefold() == str(right).strip().casefold()


def _check_conditions(cal, current):
    references = cal.get("reference_acquisitions", [])
    fields = set(_REQUIRED_CONDITIONS)
    for reference in references:
        fields.update(reference.get("conditions", {}))
    differences, missing = [], []
    for field in sorted(fields):
        values = [reference.get("conditions", {}).get(field) for reference in references]
        value = current["conditions"].get(field)
        if not values or any(item is None for item in values) or value is None:
            missing.append(field)
        known = [item for item in values if item is not None]
        if value is not None and any(not _same_condition(item, value) for item in known):
            differences.append({"field": field, "reference_values": list(dict.fromkeys(known)),
                                "acquisition_value": value})
    status = "different" if differences else "unknown" if missing else "match"
    return {**current, "status": status, "differences": differences, "missing_fields": missing}


def _prepare_targets(analyte_targets, raw_path=None):
    base_dir = None
    if isinstance(analyte_targets, pd.DataFrame):
        frame = analyte_targets.copy().reset_index(drop=True)
    else:
        path = local_path(analyte_targets).resolve()
        base_dir = path.parent
        frame = pd.read_csv(path)
    if "species" not in frame:
        frame["species"] = frame["name"] if "name" in frame else [f"Target {i + 1}" for i in range(len(frame))]
    if {"mz_min", "mz_max"} <= set(frame.columns):
        lower = pd.to_numeric(frame["mz_min"], errors="raise")
        upper = pd.to_numeric(frame["mz_max"], errors="raise")
        if ((lower <= 0) | (upper <= lower) | ~np.isfinite(lower) | ~np.isfinite(upper)).any():
            raise ValueError("Target m/z bounds must be finite, positive and increasing")
        frame["mz"] = (lower + upper) / 2
        frame["mz_window"] = (upper - lower) / 2
    for column, default in {"mz_window": "auto", "min_ms_intensity": 50.0}.items():
        frame[column] = frame[column].fillna(default) if column in frame else default
    frame = _resolve_raw_paths(frame, base_dir, local_path(raw_path).resolve() if raw_path else None)
    _validate_species(frame)
    return frame


def _resolve_mz_window(row, mw, z):
    """Resolve mz_window from CSV row: 'auto', numeric, or NaN → auto."""
    raw = row.get("mz_window", "auto")
    formula = row.get("formula", "")
    formula = "" if pd.isna(formula) else str(formula).strip()
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


def _apply_analytes(out_dir, data_dirs, ana_df, cal, method, *, plot_profiles=True):
    """Apply existing coefficients to original-intensity drift profiles."""
    warnings = list(cal.get("review_warnings", _review_warnings(cal)))
    for directory in ("csv_raw", "csv_smoothed", "png_raw", "png_smoothed", "png_overlay"):
        (out_dir / directory).mkdir(exist_ok=True)
    ana_summary, condition_checks, used_stems, native_profiles = [], [], set(), []
    reference_range = cal.get("reference_t_prime_range")
    for raw_path_str, group in ana_df.groupby("raw_path", sort=False):
        run_name = Path(raw_path_str).stem
        current = _acquisition_metadata(raw_path_str, data_dirs)
        check = _check_conditions(cal, current)
        condition_checks.append(check)
        if check["differences"]:
            fields = ", ".join(item["field"] for item in check["differences"])
            warnings.append(f"{run_name}: acquisition conditions differ from references: {fields}; review applicability.")
        if check["missing_fields"]:
            fields = ", ".join(check["missing_fields"])
            warnings.append(f"{run_name}: matching acquisition conditions cannot be verified (missing: {fields}).")
        pp = current["pusher_period_us"]
        if pp is not None and (not math.isfinite(pp) or pp <= 0):
            raise ValueError(f"{run_name}: invalid acquisition pusher period")
        im_file = _find_text_file(data_dirs, f"{run_name}_im.txt")
        ms_file = _find_text_file(data_dirs, f"{run_name}_ms.txt")
        im_data = read_im_txt(im_file) if im_file else pd.DataFrame()
        if pp is None:
            warnings.append(f"{run_name}: no acquisition pusher period; CCS was not calculated.")
        if im_file is None or ms_file is None:
            warnings.append(f"{run_name}: converted MS/IM tables are missing; CCS was not calculated.")

        for _, row in group.iterrows():
            sp = str(row["species"])
            mw, z, mz = float(row["mw"]), int(row["z"]), float(row["mz"])
            mz_w = _resolve_mz_window(row, mw, z)
            ms_int = measure_ms_intensity(ms_file, mz, mz_w) if ms_file else 0.0
            above = ms_int >= float(row["min_ms_intensity"])
            summary = {"run": run_name, "species": sp, "mw": mw, "z": z, "mz": mz,
                       "mz_window": mz_w, "method": method, "ms_intensity": ms_int, "detected": ms_int > 0,
                       "above_threshold": above, "peak_ccs": None, "peak_drift_bin": None,
                       "peak_drift_time_ms": None, "t_prime": None, "t_double_prime": None,
                       "extrapolated": None, "calibration_review_required": bool(cal.get("review_warnings")) or check["status"] != "match",
                       "conditions_status": check["status"], "pusher_period_us": pp,
                       "status": "missing_data"}
            ana_summary.append(summary)
            if pp is None or ms_file is None or im_file is None:
                continue
            if not above:
                summary["status"] = "below_threshold"
                continue
            if im_data.empty:
                summary["status"] = "no_profile"
                continue
            profile = extract_drift_profile(im_data, mz, mz_w)
            if not len(profile) or profile.max() == 0:
                summary["status"] = "no_profile"
                continue

            bins, dt_ms, ccs_n, int_n = convert_profile_to_ccs(profile, pp, mz, z, mw, cal, method)
            if not len(ccs_n):
                summary["status"] = "no_valid_ccs"
                continue
            if np.any(np.diff(ccs_n) <= 0):
                raise ValueError("Saved calibration does not produce an increasing CCS axis")
            ccs_uni, int_uni = resample_to_uniform_ccs(ccs_n, int_n)
            int_sm = smooth_ccs_profile(ccs_uni, int_uni)
            peak_ccs = float(ccs_uni[np.argmax(int_sm)]) if int_sm.max() > 0 else None
            pi = np.argmax(int_n)
            peak_dt = float(dt_ms[pi])
            tp = peak_dt - cal["edc_coeff"] * math.sqrt(mz) / 1000.0
            mu = math.sqrt(1.0 / mw + 1.0 / cal["gas_mw"])
            extrapolated = not reference_range[0] <= tp <= reference_range[1] if reference_range else None
            summary.update({"detected": True, "peak_ccs": peak_ccs, "peak_drift_bin": int(bins[pi]),
                            "peak_drift_time_ms": peak_dt, "t_prime": tp,
                            "t_double_prime": (tp ** cal["X"]) * z * mu,
                            "extrapolated": extrapolated, "status": "complete"})
            if extrapolated:
                summary["calibration_review_required"] = True
                warnings.append(f"{run_name} / {sp}: target peak is outside the calibration reference drift range (extrapolated CCS).")
            safe = re.sub(r'[<>:"/\\|?*\x00-\x1f\s]', "_", sp).strip(".") or "target"
            stem = f"{run_name}_{safe}_ccs"
            original_stem, suffix = stem, 2
            while stem in used_stems:
                stem = f"{original_stem}_{suffix}"
                suffix += 1
            used_stems.add(stem)
            summary["profile_csv"] = f"csv_raw/{stem}.csv"
            summary["smoothed_profile_csv"] = f"csv_smoothed/{stem}.csv"
            if not plot_profiles:
                # Native GUI consumes the original transformed samples directly;
                # avoid re-reading rounded CSVs or substituting the smooth grid.
                native_profiles.append({"target": dict(summary), "ccs": ccs_n,
                                        "intensity": int_n, "drift_time_ms": dt_ms})
            pd.DataFrame({"species": sp, "mw": mw, "z": z, "mz": mz,
                          "drift_bin": bins.astype(int), "drift_time_ms": dt_ms,
                          "ccs": ccs_n, "intensity": int_n}).to_csv(
                out_dir / "csv_raw" / f"{stem}.csv", index=False, float_format="%.4f")
            pd.DataFrame({"species": sp, "mw": mw, "z": z, "mz": mz,
                          "ccs": ccs_uni, "intensity": int_sm}).to_csv(
                out_dir / "csv_smoothed" / f"{stem}.csv", index=False, float_format="%.4f")
            if plot_profiles:
                for style, cr, ir, cs, ism, directory in (
                    ("Raw", ccs_n, int_n, None, None, "png_raw"),
                    ("Smoothed", None, None, ccs_uni, int_sm, "png_smoothed"),
                    ("Overlay", ccs_n, int_n, ccs_uni, int_sm, "png_overlay"),
                ):
                    plot_ccs_single(cr, ir, cs, ism, sp, run_name, style,
                                    out_dir / directory / f"{stem}.png", mw, z, mz, mz_w)
    pd.DataFrame(ana_summary).to_csv(out_dir / "analyte_summary.csv", index=False, float_format="%.4f")
    warnings = list(dict.fromkeys(warnings))
    for warning in warnings:
        _log.warning(warning)
    (out_dir / "condition_checks.json").write_text(json.dumps(condition_checks, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return {"calibration": cal, "analyte_summary": ana_summary, "warnings": warnings,
            "condition_checks": condition_checks, "native_profiles": native_profiles}


def apply_saved_calibration(calibration, out_dir, data_dirs, analyte_targets, *,
                            raw_path=None, conversion_method=None, plot_profiles=True) -> dict:
    """Reuse a reviewed curve for targets in the current acquisition or a CSV."""
    cal = load_saved_calibration(calibration)
    method = conversion_method or cal["method"]
    if method not in {"direct", "twostep"}:
        raise ValueError("CCS method must be direct or twostep")
    if method == "twostep" and ("slope_diagnostic" not in cal or "intercept_diagnostic" not in cal):
        raise ValueError("Two-step calibration needs its diagnostic slope and intercept")
    cal["method"] = method
    directories = [data_dirs] if isinstance(data_dirs, (str, Path)) else data_dirs
    data_dirs = [local_path(directory).resolve() for directory in directories]
    if not data_dirs or any(not directory.is_dir() for directory in data_dirs):
        raise FileNotFoundError("Choose a directory containing the converted MS/IM tables")
    out_dir = local_path(out_dir).resolve()
    if any(out_dir == directory or out_dir in directory.parents for directory in data_dirs):
        raise ValueError("CCS output must not replace the converted input directory")
    targets = _prepare_targets(analyte_targets, raw_path)

    def build(staged):
        staged.mkdir(parents=True)
        _setup_logging(staged)
        save_calibration(cal, staged / "calibration_curve.json", accepted=cal.get("accepted", False))
        if cal.get("calibrant_summary"):
            pd.DataFrame(cal["calibrant_summary"]).to_csv(staged / "calibrant_summary.csv", index=False, float_format="%.4f")
        return _apply_analytes(staged, data_dirs, targets, cal, method, plot_profiles=plot_profiles)

    result = _publish_results(out_dir, build)
    result["output_dir"] = str(out_dir)
    return result


def _run_calibration(out_dir: Path, data_dirs: Path | list[Path], calibrant_csv: Path,
        analyte_csv: Path | None = None,
        conversion_method: str = "direct") -> dict:
    """Run CCS calibration.

    Args:
        out_dir: Output directory for results.
        data_dirs: One or more directories containing converted _ms.txt/_im.txt files.
        calibrant_csv: CSV with calibrant species (name, mw, z, mz, CCS_reported, raw_path, ...).
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

    if analyte_csv is not None and not Path(analyte_csv).is_file():
        raise FileNotFoundError(f"Analyte table does not exist: {analyte_csv}")
    cal_df = _read_species_csv(calibrant_csv)
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
    if conversion_method not in {"direct", "twostep"}:
        raise ValueError("CCS method must be direct or twostep")

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
                value = float(params["EDC Delay Coefficient"])
                if not math.isfinite(value):
                    raise ValueError("EDC Delay Coefficient must be finite")
                if edc_coeff is not None and not math.isclose(edc_coeff, value, rel_tol=1e-6, abs_tol=1e-6):
                    raise ValueError("Calibrant acquisitions have different EDC delay coefficients; "
                                     "use standards acquired under matching instrument settings")
                edc_coeff = value
    if edc_coeff is None:
        raise ValueError("Cannot read EDC Delay Coefficient from any .raw directory. "
                         "Check raw_path column in calibrant CSV.")
    _log.info("  EDC: %.4f", edc_coeff)
    _validate_species(cal_df, calibrants=True)

    # --- Pick calibrant peaks ---
    _log.info("--- Picking calibrant peaks ---")
    cali_csv_dir = out_dir / "cali_csv"
    cali_png_dir = out_dir / "cali_png"
    cali_csv_dir.mkdir(exist_ok=True)
    cali_png_dir.mkdir(exist_ok=True)

    all_peaks, extraction = [], []
    for row_index, row in cal_df.iterrows():
        name = row["name"]
        mw = float(row["mw"]); z = int(row["z"]); mz = float(row["mz"])
        ccs = float(row["ccs"])
        raw_path = row["raw_path"]
        run = Path(raw_path).stem
        mz_min, mz_max = resolve_mz_bounds(row, mw, z)
        mz_window = _resolve_mz_window(row, mw, z)
        peak_sel = str(row["peak_select"]).strip().lower()
        peak_sel = {"latest": "last", "tallest": "highest"}.get(peak_sel, peak_sel)
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

        selected = resolve_calibrant_peak(detected, row)
        _log.info("  %s: t_D=%.3f ms, CCS_lit=%.1f, %d peaks (%s), ±%.1f Da",
                  name, selected["drift_time_ms"], ccs, len(detected), peak_sel, mz_window)

        all_peaks.append({
            "name": name, "mw": mw, "z": z, "mz": mz, "ccs": ccs,
            "mz_min": mz_min, "mz_max": mz_max,
            "drift_bin": selected["drift_bin"],
            "drift_time_ms": selected["drift_time_ms"],
            "peak_intensity": selected["peak_intensity"],
            "ms_intensity": ms_int,
        })
        extraction.append({"row_index": int(row_index), "name": str(name), "raw_path": raw_path,
                           "mz": mz, "mz_min": mz_min, "mz_max": mz_max,
                           "z": z, "mw": mw, "ccs": ccs, "mz_window": mz_window,
                           "peak_select": peak_sel, "peak_drift_bin": selected["bin_index"],
                           "drift_time_ms": selected["drift_time_ms"], "pusher_period_us": pp,
                           "min_ms_intensity": min_ms, "min_im_intensity": min_im,
                           "peak_height_frac": phf, "candidates": detected})

        safe = re.sub(r'[<>:"/\\|?*\x00-\x1f\s]', "_", str(name)).strip(".") or "reference"
        if (cali_csv_dir / f"{run}_{safe}.csv").exists():
            safe += f"_{row_index + 1}"
        pk_rows = [{"name": name, "z": z, "mz": mz, "mz_min": mz_min, "mz_max": mz_max, "ccs": ccs,
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
    cal.update({"format_version": 1, "method": conversion_method, "accepted": False,
                "reference_acquisitions": [_acquisition_metadata(path, data_dirs) for path in cal_df["raw_path"].unique()],
                "extraction": extraction})
    cal = load_saved_calibration(cal)
    save_calibration(cal, out_dir / "calibration_curve.json", accepted=False)
    summary = cal["calibrant_summary"]
    pd.DataFrame(summary).to_csv(out_dir / "calibrant_summary.csv",
                                  index=False, float_format="%.4f")
    plot_calibration(cal, out_dir / "calibration_plots.png", summary)

    # --- Apply to analytes using the same saved-model conversion path ---
    if analyte_csv and analyte_csv.exists():
        result = _apply_analytes(out_dir, data_dirs, _prepare_targets(analyte_csv), cal, conversion_method)
        cal["application_warnings"] = result["warnings"]
        cal["condition_checks"] = result["condition_checks"]

    _log.info("")
    _log.info("=== CCS Calibration complete ===")
    _log.info("  R²: %.5f, Points: %d",
              cal["r2_lnln"], cal["n_points"])
    return cal


def _publish_results(out_dir, build):
    """Replace completed results together, retaining the previous set on failure."""
    out_dir = local_path(out_dir).resolve()
    if any(parent.suffix.lower() == ".raw" and parent.is_dir()
           for parent in (out_dir, *out_dir.parents)):
        raise ValueError("Choose an output directory outside the RAW acquisition")
    out_dir.parent.mkdir(parents=True, exist_ok=True)
    if out_dir.exists() and not out_dir.is_dir():
        raise ValueError(f"CCS output is not a directory: {out_dir}")
    with tempfile.TemporaryDirectory(prefix=".waters-ccs-", dir=out_dir.parent) as temporary:
        staged = Path(temporary) / "results"
        try:
            result = build(staged)
        finally:
            # Windows cannot rename/delete an open warnings.log file.
            for handler in list(_log.handlers):
                if getattr(handler, "_waters_output_handler", False):
                    _log.removeHandler(handler)
                    handler.close()
        previous = None
        if out_dir.exists():
            previous_root = Path(tempfile.mkdtemp(prefix=f".{out_dir.name}-previous-", dir=out_dir.parent))
            previous = previous_root / out_dir.name
            out_dir.rename(previous)
        try:
            staged.rename(out_dir)
        except OSError:
            if previous is not None:
                previous.rename(out_dir)
            raise
    return result


def run(out_dir: Path, data_dirs: Path | list[Path], calibrant_csv: Path,
        analyte_csv: Path | None = None,
        conversion_method: str = "direct") -> dict:
    """Fit and optionally apply a curve, preserving the established CSV/CLI API."""
    directories = [data_dirs] if isinstance(data_dirs, (str, Path)) else data_dirs
    data_dirs = [local_path(directory).resolve() for directory in directories]
    calibrant_csv = local_path(calibrant_csv).resolve()
    analyte_csv = local_path(analyte_csv).resolve() if analyte_csv is not None else None
    return _publish_results(out_dir, lambda staged: _run_calibration(
        staged, data_dirs, calibrant_csv, analyte_csv, conversion_method))


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
    run(local_path(args.output).resolve(),
        [local_path(d).resolve() for d in args.data_dir],
        local_path(args.calibrant_csv).resolve(),
        local_path(args.analyte_csv).resolve() if args.analyte_csv else None,
        args.conversion_method)


if __name__ == "__main__":
    main()
