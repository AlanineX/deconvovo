#!/usr/bin/env python3
"""Generate quicklook tables and plots from mzML files."""

from __future__ import annotations

import argparse
import csv
import gzip
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
from pyteomics import mzml


SCAN_RE = re.compile(r"scan=(\d+)")


def parse_scan_number(scan_id: str) -> int:
    match = SCAN_RE.search(scan_id or "")
    if not match:
        return -1
    return int(match.group(1))


def write_scan_table(scan_rows: list[dict], out_csv: Path) -> None:
    fieldnames = [
        "scan_index",
        "scan_number",
        "scan_id",
        "rt_min",
        "ms_level",
        "tic",
        "base_peak_mz",
        "base_peak_intensity",
        "n_peaks",
    ]
    with out_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(scan_rows)


def write_scan_points_csv(
    out_csv: Path,
    scan_number: int,
    rt_min: float,
    ms_level: int,
    mzs: np.ndarray,
    intensities: np.ndarray,
) -> None:
    with out_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["scan_number", "rt_min", "ms_level", "mz", "intensity"])
        for mz_val, int_val in zip(mzs, intensities):
            writer.writerow([scan_number, rt_min, ms_level, float(mz_val), float(int_val)])


def build_tic_plot(scan_rows: list[dict], out_html: Path, title: str) -> None:
    x = np.asarray([row["rt_min"] for row in scan_rows], dtype=float)
    y = np.asarray([row["tic"] for row in scan_rows], dtype=float)
    fig = go.Figure(
        data=[
            go.Scatter(
                x=x,
                y=y,
                mode="lines",
                line={"width": 2.2, "color": "#0f766e"},
                fill="tozeroy",
                fillcolor="rgba(15,118,110,0.16)",
                hovertemplate="RT: %{x:.3f} min<br>TIC: %{y:.3e}<extra></extra>",
                name="TIC",
            )
        ]
    )
    fig.update_layout(
        title=title,
        xaxis_title="Retention Time (min)",
        yaxis_title="Total Ion Current",
        template="plotly_white",
        font={"family": "Helvetica, Arial, sans-serif", "size": 13},
        paper_bgcolor="#f7fafc",
        plot_bgcolor="#fdfefe",
        margin={"l": 70, "r": 20, "t": 60, "b": 60},
    )
    fig.write_html(out_html, include_plotlyjs="cdn")


def build_scan_plot(
    scan_number: int,
    rt_min: float,
    mzs: np.ndarray,
    intensities: np.ndarray,
    out_html: Path,
    title_prefix: str,
) -> None:
    if len(mzs) == 0:
        return

    keep_n = 14000
    if len(mzs) > keep_n:
        selected = np.argsort(intensities)[-keep_n:]
        mzs = mzs[selected]
        intensities = intensities[selected]

    order = np.argsort(mzs)
    mzs = mzs[order]
    intensities = intensities[order]

    stems_x = np.column_stack((mzs, mzs, np.full_like(mzs, np.nan))).ravel()
    stems_y = np.column_stack((np.zeros_like(intensities), intensities, np.full_like(intensities, np.nan))).ravel()

    max_i = float(np.nanmax(intensities)) if len(intensities) else 1.0
    norm = intensities / max_i if max_i > 0 else intensities

    fig = go.Figure()
    fig.add_trace(
        go.Scattergl(
            x=stems_x,
            y=stems_y,
            mode="lines",
            line={"width": 1.1, "color": "rgba(37,99,235,0.35)"},
            hoverinfo="skip",
            name="Sticks",
        )
    )
    fig.add_trace(
        go.Scattergl(
            x=mzs,
            y=intensities,
            mode="markers",
            marker={
                "size": 4,
                "color": norm,
                "colorscale": "Turbo",
                "showscale": True,
                "colorbar": {"title": "Relative Intensity", "len": 0.7},
                "line": {"width": 0},
            },
            hovertemplate="m/z: %{x:.5f}<br>Intensity: %{y:.3e}<extra></extra>",
            name=f"Scan {scan_number}",
        )
    )
    fig.update_layout(
        title=f"{title_prefix}: Scan {scan_number} at {rt_min:.3f} min ({len(mzs)} peaks)",
        xaxis_title="m/z",
        yaxis_title="Intensity",
        template="plotly_white",
        font={"family": "Helvetica, Arial, sans-serif", "size": 13},
        paper_bgcolor="#f5f8ff",
        plot_bgcolor="#fcfdff",
        margin={"l": 70, "r": 30, "t": 65, "b": 60},
        legend={"orientation": "h", "y": 1.03, "x": 0},
    )
    fig.update_yaxes(type="log")
    fig.write_html(out_html, include_plotlyjs="cdn")


def process_mzml_file(
    path: Path,
    out_dir: Path,
    export_all_points: bool,
    all_points_ms_level: int | None,
) -> None:
    base = path.stem
    run_dir = out_dir / base
    run_dir.mkdir(parents=True, exist_ok=True)

    scan_rows: list[dict] = []
    first_spec = None
    apex_spec = None

    with mzml.read(str(path)) as reader:
        for idx, spectrum in enumerate(reader, start=1):
            scan_id = spectrum.get("id", "")
            scan_number = parse_scan_number(scan_id)
            rt_min = float(spectrum.get("scanList", {}).get("scan", [{}])[0].get("scan start time", 0.0))
            ms_level = int(spectrum.get("ms level", 0))
            tic = float(spectrum.get("total ion current", 0.0))
            base_peak_mz = float(spectrum.get("base peak m/z", 0.0))
            base_peak_intensity = float(spectrum.get("base peak intensity", 0.0))
            mzs = np.asarray(spectrum.get("m/z array", []), dtype=float)
            intensities = np.asarray(spectrum.get("intensity array", []), dtype=float)

            row = {
                "scan_index": idx,
                "scan_number": scan_number,
                "scan_id": scan_id,
                "rt_min": rt_min,
                "ms_level": ms_level,
                "tic": tic,
                "base_peak_mz": base_peak_mz,
                "base_peak_intensity": base_peak_intensity,
                "n_peaks": int(len(mzs)),
            }
            scan_rows.append(row)
            spec = {
                "scan_index": idx,
                "scan_number": scan_number,
                "rt_min": rt_min,
                "ms_level": ms_level,
                "mzs": mzs,
                "intensities": intensities,
                "tic": tic,
            }
            if first_spec is None:
                first_spec = spec
            if apex_spec is None or tic > apex_spec["tic"]:
                apex_spec = spec

    if not scan_rows:
        raise RuntimeError(f"No spectra found in {path}")
    if first_spec is None or apex_spec is None:
        raise RuntimeError(f"Could not identify first/apex spectra in {path}")

    scan_table_csv = run_dir / f"{base}.scan_table.csv"
    write_scan_table(scan_rows, scan_table_csv)

    middle_row = scan_rows[len(scan_rows) // 2]
    middle_spec = None

    need_second_pass = export_all_points or (
        middle_row["scan_index"] not in {first_spec["scan_index"], apex_spec["scan_index"]}
    )
    if need_second_pass:
        all_points_path = run_dir / f"{base}.all_mz_intensity.csv.gz"
        with mzml.read(str(path)) as reader:
            if export_all_points:
                handle = gzip.open(all_points_path, "wt", newline="", encoding="utf-8")
                writer = csv.writer(handle)
                writer.writerow(["scan_number", "rt_min", "ms_level", "mz", "intensity"])
            else:
                handle = None
                writer = None

            try:
                for idx, spectrum in enumerate(reader, start=1):
                    scan_number = parse_scan_number(spectrum.get("id", ""))
                    rt_min = float(
                        spectrum.get("scanList", {}).get("scan", [{}])[0].get("scan start time", 0.0)
                    )
                    ms_level = int(spectrum.get("ms level", 0))
                    mzs = np.asarray(spectrum.get("m/z array", []), dtype=float)
                    intensities = np.asarray(spectrum.get("intensity array", []), dtype=float)

                    if writer is not None and (
                        all_points_ms_level is None or ms_level == all_points_ms_level
                    ):
                        for mz_val, int_val in zip(mzs, intensities):
                            writer.writerow([scan_number, rt_min, ms_level, float(mz_val), float(int_val)])

                    if idx == middle_row["scan_index"] and middle_spec is None:
                        middle_spec = {
                            "scan_index": idx,
                            "scan_number": scan_number,
                            "rt_min": rt_min,
                            "ms_level": ms_level,
                            "mzs": mzs,
                            "intensities": intensities,
                        }

                    if writer is None and middle_spec is not None:
                        break
            finally:
                if handle is not None:
                    handle.close()
    else:
        if middle_row["scan_index"] == first_spec["scan_index"]:
            middle_spec = first_spec
        else:
            middle_spec = apex_spec

    if middle_spec is None:
        raise RuntimeError(f"Could not retrieve middle scan from {path}")

    write_scan_points_csv(
        run_dir / f"{base}.scan_{first_spec['scan_number']}_mz_intensity.csv",
        first_spec["scan_number"],
        first_spec["rt_min"],
        first_spec["ms_level"],
        first_spec["mzs"],
        first_spec["intensities"],
    )
    write_scan_points_csv(
        run_dir / f"{base}.scan_{middle_spec['scan_number']}_mz_intensity.csv",
        middle_spec["scan_number"],
        middle_spec["rt_min"],
        middle_spec["ms_level"],
        middle_spec["mzs"],
        middle_spec["intensities"],
    )
    write_scan_points_csv(
        run_dir / f"{base}.scan_{apex_spec['scan_number']}_mz_intensity.csv",
        apex_spec["scan_number"],
        apex_spec["rt_min"],
        apex_spec["ms_level"],
        apex_spec["mzs"],
        apex_spec["intensities"],
    )

    build_tic_plot(scan_rows, run_dir / f"{base}.tic.html", f"{base} TIC")
    build_scan_plot(
        apex_spec["scan_number"],
        apex_spec["rt_min"],
        apex_spec["mzs"],
        apex_spec["intensities"],
        run_dir / f"{base}.apex_scan_mz_intensity.html",
        title_prefix=base,
    )


def collect_input_files(inputs: list[str]) -> list[Path]:
    results: list[Path] = []
    for item in inputs:
        path = Path(item)
        if path.is_dir():
            for candidate in sorted(path.glob("*.mzML")):
                results.append(candidate.resolve())
        elif path.is_file() and path.suffix.lower() == ".mzml":
            results.append(path.resolve())
    dedup = []
    seen = set()
    for p in results:
        if p not in seen:
            dedup.append(p)
            seen.add(p)
    return dedup


def run_is_complete(run_dir: Path, base: str, export_all_points: bool) -> bool:
    if not run_dir.exists():
        return False
    required = [
        run_dir / f"{base}.scan_table.csv",
        run_dir / f"{base}.tic.html",
        run_dir / f"{base}.apex_scan_mz_intensity.html",
    ]
    if not all(p.exists() for p in required):
        return False
    scan_csvs = list(run_dir.glob(f"{base}.scan_*_mz_intensity.csv"))
    if len(scan_csvs) < 3:
        return False
    if export_all_points and not (run_dir / f"{base}.all_mz_intensity.csv.gz").exists():
        return False
    return True


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate mzML quicklook outputs.")
    parser.add_argument(
        "inputs",
        nargs="+",
        help="mzML file(s) or directory/directories containing mzML files",
    )
    parser.add_argument(
        "-o",
        "--out-dir",
        default="output/quicklook",
        help="Output directory (default: output/quicklook)",
    )
    parser.add_argument(
        "--export-all-points",
        action="store_true",
        help="Export full long-format m/z-intensity table as .csv.gz for each file",
    )
    parser.add_argument(
        "--all-points-ms-level",
        type=int,
        default=None,
        help="Optional MS level filter for --export-all-points (example: 1)",
    )
    parser.add_argument(
        "--skip-existing",
        action="store_true",
        help="Skip files whose quicklook outputs already exist",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=1,
        help="Parallel workers across mzML files (max 8)",
    )
    args = parser.parse_args()

    files = collect_input_files(args.inputs)
    if not files:
        raise SystemExit("No mzML files found in provided inputs.")

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    to_process: list[Path] = []
    for file_path in files:
        base = file_path.stem
        run_dir = out_dir / base
        if args.skip_existing and run_is_complete(
            run_dir=run_dir,
            base=base,
            export_all_points=args.export_all_points,
        ):
            print(f"Skipped {file_path}")
            continue
        to_process.append(file_path)

    if not to_process:
        print("Nothing to process.")
        return

    workers = min(max(1, args.jobs), 8, len(to_process))
    if workers == 1:
        for file_path in to_process:
            process_mzml_file(
                file_path,
                out_dir,
                export_all_points=args.export_all_points,
                all_points_ms_level=args.all_points_ms_level,
            )
            print(f"Processed {file_path}")
        return

    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(
                process_mzml_file,
                file_path,
                out_dir,
                args.export_all_points,
                args.all_points_ms_level,
            ): file_path
            for file_path in to_process
        }
        for fut in as_completed(futures):
            file_path = futures[fut]
            fut.result()
            print(f"Processed {file_path}")


if __name__ == "__main__":
    main()
