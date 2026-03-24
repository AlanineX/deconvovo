#!/usr/bin/env python3
"""
Waters SYNAPT G2 IM-MS analysis pipeline.

End-to-end: read MassLynx .raw -> calibrate -> deconvolve -> IMS drift profiles
            -> per-species drift extraction -> 2D heatmaps -> report

Usage:
    # Full pipeline on a directory of .raw folders
    .venv-ms/bin/python scripts/waters_pipeline.py data_2_waters/20260216 -o output/20_waters_deconv

    # Single .raw directory with custom species DB
    .venv-ms/bin/python scripts/waters_pipeline.py data.raw -o out --species-db adp

    # Protein runs only, skip IMS
    .venv-ms/bin/python scripts/waters_pipeline.py data/ -o out --species-db ubiquitin --skip-ims

    # Custom species from JSON
    .venv-ms/bin/python scripts/waters_pipeline.py data/ -o out --species-db my_species.json
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd

# Import the core library (sibling in scripts/)
sys.path.insert(0, str(Path(__file__).resolve().parent))
import masslynx_reader as mlr


# =============================================================================
# Input collection
# =============================================================================

def collect_raw_dirs(inputs: list[str]) -> list[Path]:
    """Accept .raw directories or parent directories containing them."""
    results: list[Path] = []
    for item in inputs:
        p = Path(item).resolve()
        if p.is_dir() and p.suffix.lower() == ".raw":
            results.append(p)
        elif p.is_dir():
            for child in sorted(p.iterdir()):
                if child.is_dir() and child.suffix.lower() == ".raw":
                    results.append(child)
    # Deduplicate
    seen: set[Path] = set()
    deduped: list[Path] = []
    for r in results:
        if r not in seen:
            seen.add(r)
            deduped.append(r)
    return deduped


# =============================================================================
# Per-run analysis
# =============================================================================

def analyze_run(
    raw_dir: Path,
    out_dir: Path,
    species_db_source: str | None,
    drift_mz_tol: float,
    species_tol: float,
    skip_ims: bool,
    skip_deconv: bool,
    skip_plots: bool,
    skip_existing: bool,
) -> dict:
    """Full analysis of one .raw directory."""
    run_name = raw_dir.stem
    run_class = mlr.classify_run(run_name)
    buf, ph, conc = mlr.extract_buffer_info(run_name)
    is_protein = "protein" in run_class or "calibrant" in run_class

    # Skip check
    marker = out_dir / f"{run_name}_spectrum.csv"
    if skip_existing and marker.exists():
        print(f"  {run_name} — skipped (exists)")
        return {"run_name": run_name, "status": "skipped"}

    print(f"  {run_name} [{run_class}]", end="", flush=True)

    # Resolve species DB
    if species_db_source:
        species_db = mlr.load_species_db(species_db_source)
    else:
        species_db = mlr.auto_species_db(run_class)

    result: dict = {
        "run_name": run_name,
        "run_class": run_class,
        "buffer": buf,
        "pH": ph,
        "buffer_conc": conc,
        "status": "pending",
    }

    # --- Metadata ---
    header = mlr.read_header(raw_dir)
    params = mlr.read_extern_inf(raw_dir)
    result["description"] = header.get("Sample Description", "")
    result["instrument"] = header.get("Instrument", "")

    # --- FUNC001: overview spectrum ---
    # Data is already vendor-calibrated to m/z by the rainbow reader
    scans_f1 = mlr.read_func_scans(raw_dir, func=1)
    if not scans_f1:
        print(" — no FUNC001")
        result["status"] = "no_data"
        return result

    mz, intensities = mlr.sum_scans(scans_f1)
    result["n_scans_f1"] = len(scans_f1)

    # Calibration is handled by vendor polynomial in the reader
    cal = mlr.Calibration(method="vendor_polynomial")
    result["cal_method"] = cal.method

    valid = np.isfinite(mz) & np.isfinite(intensities) & (mz > 0) & (mz < 5000)
    mz, intensities = mz[valid], intensities[valid]

    if len(mz) == 0:
        print(" — no valid m/z data")
        result["status"] = "no_valid_data"
        return result

    result["mz_min"] = float(mz.min())
    result["mz_max"] = float(mz.max())
    result["max_intensity"] = float(intensities.max())

    # Save spectrum CSV
    pd.DataFrame({"mz": mz, "intensity": intensities}).to_csv(
        out_dir / f"{run_name}_spectrum.csv", index=False, float_format="%.4f")

    # --- Deconvolution ---
    assignments: list[dict] = []
    deconv_results: list[dict] = []

    if not skip_deconv:
        if is_protein:
            if "ubq" in run_class:
                deconv_results = mlr.deconvolve_protein(mz, intensities, expected_mass=8565.8)
            elif "myo" in run_class:
                deconv_results = mlr.deconvolve_protein(mz, intensities, expected_mass=16951.5)
            elif "cytc" in run_class:
                deconv_results = mlr.deconvolve_protein(mz, intensities, expected_mass=12384.0)
            else:
                deconv_results = mlr.deconvolve_protein(mz, intensities)
            result["n_protein_species"] = len(deconv_results)
        else:
            peaks = mlr.detect_peaks(mz, intensities)
            groups = mlr.group_isotopes(peaks)
            assignments = mlr.assign_species(groups, species_db, tol=species_tol)
            n_assigned = sum(1 for a in assignments if a["species"] != "unassigned")
            n_unassigned = sum(1 for a in assignments if a["species"] == "unassigned")
            result["n_assigned"] = n_assigned
            result["n_unassigned"] = n_unassigned

        # Save deconv CSV
        if assignments:
            pd.DataFrame([{
                "mono_mz": a["mono_mz"], "species": a["species"],
                "expected_mz": a["expected_mz"] or "", "charge": a["charge"],
                "mass_error": a["mass_error"] if a["mass_error"] is not None else "",
                "max_intensity": a["max_int"], "total_area": a["total_area"],
                "n_isotopes": a["n_iso"],
            } for a in assignments]).to_csv(
                out_dir / f"{run_name}_deconv.csv", index=False, float_format="%.4f")
        elif deconv_results:
            pd.DataFrame([{
                "neutral_mass": d["neutral_mass"], "mass_std": d["mass_std"],
                "n_charges": d["n_charges"],
                "charge_range": f"{d['charge_min']}-{d['charge_max']}",
                "max_intensity": d["max_intensity"],
            } for d in deconv_results]).to_csv(
                out_dir / f"{run_name}_deconv.csv", index=False, float_format="%.2f")

    # --- Plots (spectrum) ---
    if not skip_plots:
        fig = mlr.plot_spectrum(mz, intensities, f"{run_name} — Spectrum",
                                annotations=assignments if assignments else None)
        fig.write_html(str(out_dir / f"{run_name}_spectrum.html"))

    # --- FUNC002: IMS ---
    has_ims = False
    if not skip_ims:
        scans_f2 = mlr.read_func_scans(raw_dir, func=2)
        if scans_f2:
            total_tic = sum(s.total_intensity for s in scans_f2)
            if total_tic > 0:
                has_ims = True
                n_drift = len(scans_f2)
                result["n_drift_bins"] = n_drift
                result["total_ims_tic"] = total_tic

                # Drift TIC profile
                tic_arr = np.array([s.total_intensity for s in scans_f2])
                peak_bin = int(np.argmax(tic_arr))
                result["peak_drift_bin"] = peak_bin
                result["peak_drift_int"] = float(tic_arr[peak_bin])

                # FWHM
                hm = tic_arr[peak_bin] / 2
                above = tic_arr >= hm
                if above.any():
                    result["drift_fwhm"] = int(
                        len(above) - 1 - np.argmax(above[::-1]) - np.argmax(above))
                else:
                    result["drift_fwhm"] = 0

                # Per-species drift profiles
                drift_df = mlr.extract_species_drift_profiles(
                    scans_f2, cal, species_db, mz_tol=drift_mz_tol)
                drift_df.to_csv(out_dir / f"{run_name}_drift_profiles.csv",
                                index=False, float_format="%.2f")

                # Auto-detect m/z range from data
                all_mz: list[float] = []
                for s in scans_f2:
                    if len(s.channels) > 0:
                        v = cal.to_mz(s.channels)
                        good = v[(v > 0) & (v < 5000) & np.isfinite(v)]
                        if len(good) > 0:
                            all_mz.extend([float(good.min()), float(good.max())])

                if all_mz:
                    dmin, dmax = min(all_mz), max(all_mz)
                    pad = (dmax - dmin) * 0.05
                    mz_range_full = (max(0, dmin - pad), dmax + pad)
                else:
                    mz_range_full = (100, 1000) if not is_protein else (200, 3000)

                # 2D high-res CSV (for analysis)
                heatmap_hr, mz_hr = mlr.build_2d_data_hires(
                    scans_f2, cal, mz_range=mz_range_full)
                result["n_mz_bins"] = heatmap_hr.shape[1]

                rows_2d = []
                for di in range(heatmap_hr.shape[0]):
                    for mi in range(heatmap_hr.shape[1]):
                        if heatmap_hr[di, mi] > 0:
                            rows_2d.append((di, round(mz_hr[mi], 4),
                                            round(heatmap_hr[di, mi], 2)))
                pd.DataFrame(rows_2d, columns=["drift_bin", "mz", "intensity"]).to_csv(
                    out_dir / f"{run_name}_2d_imms.csv", index=False, float_format="%.4f")

                # Plots
                if not skip_plots:
                    # Drift profiles (per-species)
                    fig_drift = mlr.plot_drift_profiles(drift_df, run_name)
                    fig_drift.write_html(str(out_dir / f"{run_name}_drift_profiles.html"))

                    # 2D IM-MS scatter plot — raw data, full range, no smoothing
                    fig_2d = mlr.plot_2d_imms_scatter(
                        scans_f2, cal, run_name)
                    fig_2d.write_html(str(out_dir / f"{run_name}_2d_imms.html"))

                    # Zoomed scatter for ADP region
                    if not is_protein:
                        fig_zoom = mlr.plot_2d_imms_scatter(
                            scans_f2, cal, f"{run_name} (ADP region)",
                            mz_range=(400, 510))
                        fig_zoom.write_html(
                            str(out_dir / f"{run_name}_2d_imms_zoom.html"))

    result["has_ims"] = has_ims
    result["status"] = "complete"
    print(f" — done" + (f" (IMS: {result.get('n_drift_bins', 0)} bins)" if has_ims else ""))
    return result


# =============================================================================
# Report generation
# =============================================================================

def build_adp_inventory(all_results: list[dict]) -> pd.DataFrame:
    """Cross-run ADP adduct inventory."""
    # This is generated from the deconv CSVs after the pipeline runs
    rows = []
    for r in all_results:
        if "adp" not in r.get("run_class", ""):
            continue
        out_dir = r.get("_out_dir")
        if not out_dir:
            continue
        deconv_path = Path(out_dir) / f"{r['run_name']}_deconv.csv"
        if not deconv_path.exists():
            continue

        df = pd.read_csv(deconv_path)
        row: dict = {
            "run": r["run_name"],
            "buffer": r.get("buffer", ""),
            "pH": r.get("pH", ""),
            "buffer_conc": r.get("buffer_conc", ""),
        }

        for key in ["ADP+H", "ADP+Na", "ADP+K", "ADP+NH4", "ADP+2Na-H",
                     "ADP+EDA+H", "ADP+EDA+Na", "ADP+EDA+2H [2+]",
                     "ADP+EDA+H+Na [2+]", "ADP+2EDA+2H [2+]",
                     "2ADP+H", "2ADP+Na", "EDA+H", "unassigned"]:
            subset = df[df["species"] == key]
            row[f"{key}_count"] = len(subset)
            row[f"{key}_max_int"] = float(subset["max_intensity"].max()) if len(subset) > 0 else 0.0

        eda_subset = df[df["species"].str.contains("EDA", na=False)]
        row["total_EDA_complexes"] = len(eda_subset)
        row["max_EDA_intensity"] = float(eda_subset["max_intensity"].max()) if len(eda_subset) > 0 else 0.0

        assigned = df[df["species"] != "unassigned"]
        row["total_assigned"] = len(assigned)
        row["total_unassigned"] = len(df) - len(assigned)
        rows.append(row)

    return pd.DataFrame(rows)


def generate_report(all_results: list[dict], out_dir: Path) -> str:
    """Generate markdown report."""
    lines: list[str] = []
    lines.append("# Waters SYNAPT G2 IM-MS Pipeline Report\n")

    n_total = len(all_results)
    n_ok = sum(1 for r in all_results if r.get("status") == "complete")
    adp_runs = [r for r in all_results if "adp" in r.get("run_class", "")]
    protein_runs = [r for r in all_results if "protein" in r.get("run_class", "")]
    cal_runs = [r for r in all_results if "calibrant" in r.get("run_class", "")]

    lines.append(f"Processed **{n_ok}/{n_total}** runs: "
                 f"{len(adp_runs)} ADP, {len(protein_runs)} protein, {len(cal_runs)} calibrant\n")

    # Summary table
    lines.append("## Run Summary\n")
    lines.append("| Run | Class | Cal | IMS Bins | Peak Bin | FWHM | Status |")
    lines.append("|-----|-------|-----|----------|----------|------|--------|")
    for r in all_results:
        name = r["run_name"][:50]
        cls = r.get("run_class", "")
        cal_m = r.get("cal_method", "")
        n_d = r.get("n_drift_bins", "")
        pk = r.get("peak_drift_bin", "")
        fw = r.get("drift_fwhm", "")
        st = r.get("status", "")
        lines.append(f"| {name} | {cls} | {cal_m} | {n_d} | {pk} | {fw} | {st} |")
    lines.append("")

    # IMS drift profiles
    ims_runs = [r for r in all_results if r.get("has_ims")]
    if ims_runs:
        lines.append("## Ion Mobility Drift Profiles\n")
        lines.append("Per-species drift profiles are in `{run_name}_drift_profiles.csv` "
                      "and `{run_name}_drift_profiles.html`.\n")
        lines.append("| Run | Peak Bin | FWHM | IMS TIC |")
        lines.append("|-----|----------|------|---------|")
        for r in ims_runs:
            tic = r.get("total_ims_tic", 0)
            tic_str = f"{tic/1e6:.1f}M" if tic > 1e6 else f"{tic:.0f}"
            lines.append(f"| {r['run_name'][:50]} | {r.get('peak_drift_bin', '')} "
                         f"| {r.get('drift_fwhm', '')} | {tic_str} |")
        lines.append("")

    # Per-run output listing
    lines.append("## Output Files Per Run\n")
    lines.append("| File | Content |")
    lines.append("|------|---------|")
    lines.append("| `{name}_spectrum.csv/html` | Averaged m/z vs intensity (FUNC001) |")
    lines.append("| `{name}_deconv.csv` | Deconvolution: species assignments |")
    lines.append("| `{name}_drift_profiles.csv/html` | Per-species drift time profiles |")
    lines.append("| `{name}_2d_imms.csv/html` | 2D IM-MS heatmap (drift x m/z) |")
    lines.append("| `{name}_2d_imms_zoom.html` | Zoomed ADP region (m/z 400-510) |")
    lines.append("")

    return "\n".join(lines)


# =============================================================================
# Main
# =============================================================================

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Waters SYNAPT G2 IM-MS analysis pipeline.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("inputs", nargs="+",
                        help=".raw directories or parent dirs containing .raw")
    parser.add_argument("-o", "--out-dir", required=True,
                        help="Output directory")
    parser.add_argument("--species-db", default=None,
                        help="Species DB: preset name (adp, ubiquitin, ...) or JSON path")
    parser.add_argument("--species-tol", type=float, default=0.8,
                        help="m/z tolerance for species matching (Da, default: 0.8)")
    parser.add_argument("--drift-mz-tol", type=float, default=1.0,
                        help="m/z tolerance for per-species drift extraction (Da, default: 1.0)")
    parser.add_argument("--skip-ims", action="store_true",
                        help="Skip FUNC002 ion mobility analysis")
    parser.add_argument("--skip-deconv", action="store_true",
                        help="Skip deconvolution/assignment")
    parser.add_argument("--skip-plots", action="store_true",
                        help="Skip HTML plot generation (CSV only)")
    parser.add_argument("--skip-existing", action="store_true",
                        help="Skip runs whose outputs already exist")
    parser.add_argument("--jobs", type=int, default=1,
                        help="Parallel workers (default: 1, max 8)")

    args = parser.parse_args()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    raw_dirs = collect_raw_dirs(args.inputs)
    if not raw_dirs:
        print("No .raw directories found.")
        sys.exit(1)

    print("=" * 70)
    print("Waters SYNAPT G2 — IM-MS Pipeline")
    print("=" * 70)
    print(f"\n{len(raw_dirs)} .raw directories -> {out_dir}\n")

    all_results: list[dict] = []

    for raw_dir in raw_dirs:
        try:
            result = analyze_run(
                raw_dir, out_dir,
                species_db_source=args.species_db,
                drift_mz_tol=args.drift_mz_tol,
                species_tol=args.species_tol,
                skip_ims=args.skip_ims,
                skip_deconv=args.skip_deconv,
                skip_plots=args.skip_plots,
                skip_existing=args.skip_existing,
            )
            result["_out_dir"] = str(out_dir)
            all_results.append(result)
        except Exception as e:
            import traceback
            print(f"  ERROR: {e}")
            traceback.print_exc()
            all_results.append({"run_name": raw_dir.stem, "status": f"error: {e}"})

    # --- Aggregate outputs ---
    # Summary CSV
    summary_df = pd.DataFrame([{k: v for k, v in r.items() if k != "_out_dir"}
                                for r in all_results])
    summary_path = out_dir / "pipeline_summary.csv"
    summary_df.to_csv(summary_path, index=False)
    print(f"\nWrote: {summary_path}")

    # ADP adduct inventory
    adp_inv = build_adp_inventory(all_results)
    if not adp_inv.empty:
        inv_path = out_dir / "adp_adduct_inventory.csv"
        adp_inv.to_csv(inv_path, index=False)
        print(f"Wrote: {inv_path}")

    # Report
    report = generate_report(all_results, out_dir)
    report_path = out_dir / "PIPELINE_REPORT.md"
    report_path.write_text(report)
    print(f"Wrote: {report_path}")

    # Final stats
    n_ok = sum(1 for r in all_results if r.get("status") == "complete")
    n_ims = sum(1 for r in all_results if r.get("has_ims"))
    print(f"\n{'=' * 70}")
    print(f"Pipeline complete: {n_ok}/{len(all_results)} runs, {n_ims} with IMS data")
    print(f"Output: {out_dir}")


if __name__ == "__main__":
    main()
