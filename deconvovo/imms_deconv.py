"""UniDec deconvolution for IM-MS data.

Wraps UniDec to run mass deconvolution on MS text files.
Includes batch orchestration and single-file engine.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

# numpy 2.x compat for UniDec
if not hasattr(np, 'trapz'):
    np.trapz = np.trapezoid


# =============================================================================
# Engine — single-file deconvolution
# =============================================================================

def _build_config(mass_range: tuple | None, charge_range: tuple | None,
                  mass_bins: float | None) -> dict:
    """Return UniDec config. User provides ranges via CLI args."""
    return {
        "minmz": 50, "maxmz": 4000,
        "masslb": mass_range[0] if mass_range else 50,
        "massub": mass_range[1] if mass_range else 100000,
        "startz": charge_range[0] if charge_range else 1,
        "endz": charge_range[1] if charge_range else 50,
        "massbins": mass_bins or 1,
        "peakwindow": 50, "peakthresh": 0.01,
    }


def run_unidec_ms(ms_file: Path, run_name: str, out_dir: Path,
                   mass_range: tuple | None, charge_range: tuple | None,
                   mass_bins: float | None) -> dict:
    """Run UniDec deconvolution on a single MS text file."""
    from unidec import engine as eng

    u = eng.UniDec()
    u.open_file(str(ms_file))

    cfg = _build_config(mass_range, charge_range, mass_bins)
    for key, val in cfg.items():
        setattr(u.config, key, val)

    u.autorun(silent=True)

    result = {
        "run_name": run_name,
        "n_peaks": len(u.pks.peaks),
        "r_squared": getattr(u, 'rsquared', None),
    }

    peaks = []
    for pk in u.pks.peaks:
        peaks.append({
            "mass": pk.mass, "height": pk.height,
            "avgcharge": pk.avgcharge,
            "score": getattr(pk, 'dscore', None),
            "area": getattr(pk, 'area', None),
        })
    result["peaks"] = peaks

    if hasattr(u.data, 'massdat') and u.data.massdat is not None:
        md = u.data.massdat
        pd.DataFrame({"mass": md[:, 0], "intensity": md[:, 1]}).to_csv(
            out_dir / f"{run_name}_deconv_mass.csv", index=False, float_format="%.4f")

    if peaks:
        pd.DataFrame(peaks).to_csv(
            out_dir / f"{run_name}_peaks.csv", index=False, float_format="%.4f")

    if hasattr(u.data, 'data2') and u.data.data2 is not None:
        pd.DataFrame({"mz": u.data.data2[:, 0], "intensity": u.data.data2[:, 1]}).to_csv(
            out_dir / f"{run_name}_spectrum.csv", index=False, float_format="%.4f")

    try:
        u.gen_html_report()
        src_report = Path(u.config.dirname) / u.config.reportfile
        if src_report.exists():
            import shutil
            shutil.copy2(src_report, out_dir / f"{run_name}_unidec_report.html")
    except Exception:
        pass

    return result


# =============================================================================
# Batch orchestration
# =============================================================================

def run(data_dir: Path, out_dir: Path, mass_range: tuple | None = None,
        charge_range: tuple | None = None, mass_bins: float | None = None,
        skip_existing: bool = False) -> list[dict]:
    """Run UniDec deconvolution on all _ms.txt files."""
    out_dir.mkdir(parents=True, exist_ok=True)

    try:
        from unidec import engine as _eng  # noqa: F401
    except ImportError:
        print("  UniDec not installed — skipping")
        return []

    ms_files = sorted(data_dir.glob("*_ms.txt"))
    print(f"  {len(ms_files)} runs")

    results = []
    for ms_file in ms_files:
        run_name = ms_file.stem.replace("_ms", "")

        if skip_existing and (out_dir / f"{run_name}_peaks.csv").exists():
            continue

        print(f"    {run_name}", end="", flush=True)
        try:
            result = run_unidec_ms(ms_file, run_name, out_dir,
                                    mass_range, charge_range, mass_bins)
            print(f" — {result['n_peaks']} peaks")
            results.append(result)
        except Exception as e:
            print(f" — ERROR: {e}")
            results.append({"run_name": run_name, "error": str(e)})

    return results


def main() -> None:
    parser = argparse.ArgumentParser(description="UniDec deconvolution")
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--mass-range", type=float, nargs=2, default=None)
    parser.add_argument("--charge-range", type=int, nargs=2, default=None)
    parser.add_argument("--mass-bins", type=float, default=None)
    parser.add_argument("--skip-existing", action="store_true")
    args = parser.parse_args()
    run(Path(args.input).resolve(), Path(args.output).resolve(),
        args.mass_range, args.charge_range, args.mass_bins, args.skip_existing)


if __name__ == "__main__":
    main()
