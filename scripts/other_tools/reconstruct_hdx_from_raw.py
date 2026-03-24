#!/usr/bin/env python3
"""Reconstruct peptide-level HDX metrics from raw-derived mzML only.

This script uses only peptide sequence/start/end as targets and searches
MS1 scans for isotopic envelopes to estimate:
- peptide abundance proxy (envelope intensity at best scan)
- neutral mass estimate
- delta mass vs baseline (default peptide_mapping run)
- exchange ratio vs max-timepoint reference
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from pyteomics import mass, mzml


PROTON = 1.007276466812
NEUTRON = 1.0033548378

TIME_RE = re.compile(
    r"(?P<value>\d+(?:\.\d+)?)\s*(?P<unit>min|h|hr|hrs|hour|hours|d|day|days)",
    re.IGNORECASE,
)


@dataclass(frozen=True)
class RunInfo:
    run_id: str
    mzml_path: Path
    time_label: str
    time_min: float | None
    is_mapping: bool


def parse_time_from_name(name: str) -> tuple[str, float | None, bool]:
    lower = name.lower()
    if "peptide_mapping" in lower:
        return "peptide_mapping", None, True
    m = TIME_RE.search(lower)
    if not m:
        return "unknown", None, False
    value = float(m.group("value"))
    unit = m.group("unit")
    if unit == "min":
        return f"{value:g}min", value, False
    if unit in {"h", "hr", "hrs", "hour", "hours"}:
        return f"{value:g}h", value * 60.0, False
    if unit in {"d", "day", "days"}:
        return f"{value:g}d", value * 1440.0, False
    return "unknown", None, False


def collect_runs(mzml_dir: Path) -> list[RunInfo]:
    runs: list[RunInfo] = []
    for p in sorted(mzml_dir.glob("*.mzML")):
        label, tmin, is_map = parse_time_from_name(p.name)
        runs.append(
            RunInfo(
                run_id=p.stem,
                mzml_path=p.resolve(),
                time_label=label,
                time_min=tmin,
                is_mapping=is_map,
            )
        )
    if not runs:
        raise ValueError(f"No .mzML files found in {mzml_dir}")
    return runs


def load_targets(peptide_source: Path) -> pd.DataFrame:
    ext = peptide_source.suffix.lower()
    if ext == ".xlsx":
        df = pd.read_excel(
            peptide_source,
            sheet_name=0,
            usecols=["Sequence", "Start", "End"],
            engine="openpyxl",
        )
    elif ext == ".csv":
        df = pd.read_csv(peptide_source, usecols=["Sequence", "Start", "End"])
    else:
        raise ValueError("peptide source must be .xlsx or .csv")

    df = df.rename(columns={"Sequence": "sequence", "Start": "start", "End": "end"})
    df = df.dropna(subset=["sequence", "start", "end"]).copy()
    df["sequence"] = df["sequence"].astype(str).str.strip()
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    df = df.dropna(subset=["start", "end"])
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    df = df[df["sequence"] != ""]
    df = df.drop_duplicates(subset=["sequence", "start", "end"]).reset_index(drop=True)
    if df.empty:
        raise ValueError("No peptide targets found after cleaning.")

    df["fragment_id"] = df["sequence"] + "|" + df["start"].astype(str) + "-" + df["end"].astype(str)
    df["neutral_mass_theoretical"] = df["sequence"].apply(lambda s: mass.calculate_mass(sequence=s))
    return df.sort_values(["start", "end", "sequence"]).reset_index(drop=True)


def match_targets(
    mz_arr: np.ndarray,
    int_arr: np.ndarray,
    targets: np.ndarray,
    ppm: float,
) -> tuple[np.ndarray, np.ndarray]:
    n = mz_arr.size
    idx = np.searchsorted(mz_arr, targets)
    out_mz = np.full(targets.shape[0], np.nan, dtype=float)
    out_int = np.zeros(targets.shape[0], dtype=float)

    for i, pos in enumerate(idx):
        tol = targets[i] * ppm * 1e-6
        best_j = -1
        best_i = -1.0
        for j in (pos - 1, pos, pos + 1):
            if j < 0 or j >= n:
                continue
            err = abs(mz_arr[j] - targets[i])
            if err <= tol and int_arr[j] > best_i:
                best_i = float(int_arr[j])
                best_j = j
        if best_j >= 0:
            out_mz[i] = float(mz_arr[best_j])
            out_int[i] = best_i
    return out_mz, out_int


def best_feature_for_peptide(
    mz_arr: np.ndarray,
    int_arr: np.ndarray,
    theo_by_charge: dict[int, np.ndarray],
    min_isotopes: int,
    ppm: float,
) -> dict | None:
    weights = np.array([1.0, 0.7, 0.5, 0.35, 0.25, 0.18, 0.13, 0.1], dtype=float)
    best: dict | None = None

    for z, theo in theo_by_charge.items():
        matched_mz, matched_int = match_targets(mz_arr, int_arr, theo, ppm=ppm)
        k = int((matched_int > 0).sum())
        if k < min_isotopes:
            continue
        w = weights[: matched_int.size]
        score = float((matched_int * w).sum() * (k / matched_int.size))
        if best is None or score > best["score"]:
            iso_idx = np.arange(matched_int.size, dtype=float)
            used = matched_int > 0
            neutral_series = (matched_mz[used] * z - z * PROTON) - iso_idx[used] * NEUTRON
            neutral_mass_est = float(np.average(neutral_series, weights=matched_int[used]))
            best = {
                "charge": int(z),
                "score": score,
                "matched_isotopes": k,
                "envelope_intensity": float(matched_int.sum()),
                "neutral_mass_est": neutral_mass_est,
                "mono_mz_matched": float(matched_mz[0]) if np.isfinite(matched_mz[0]) else np.nan,
            }
    return best


def search_run(
    run: RunInfo,
    targets: pd.DataFrame,
    charge_min: int,
    charge_max: int,
    n_isotopes: int,
    min_isotopes: int,
    ppm: float,
) -> list[dict]:
    theo: dict[str, dict[int, np.ndarray]] = {}
    for r in targets.itertuples():
        per_z: dict[int, np.ndarray] = {}
        for z in range(charge_min, charge_max + 1):
            mono = (r.neutral_mass_theoretical + z * PROTON) / z
            per_z[z] = mono + (np.arange(n_isotopes, dtype=float) * NEUTRON / z)
        theo[r.fragment_id] = per_z

    best_by_peptide: dict[str, dict] = {}
    with mzml.read(str(run.mzml_path)) as reader:
        for spec in reader:
            if int(spec.get("ms level", 0)) != 1:
                continue
            mz_arr = np.asarray(spec.get("m/z array", []), dtype=float)
            int_arr = np.asarray(spec.get("intensity array", []), dtype=float)
            if mz_arr.size == 0:
                continue
            scan_id = str(spec.get("id", ""))
            rt = float(spec.get("scanList", {}).get("scan", [{}])[0].get("scan start time", np.nan))
            for r in targets.itertuples():
                feat = best_feature_for_peptide(
                    mz_arr=mz_arr,
                    int_arr=int_arr,
                    theo_by_charge=theo[r.fragment_id],
                    min_isotopes=min_isotopes,
                    ppm=ppm,
                )
                if feat is None:
                    continue
                key = r.fragment_id
                prev = best_by_peptide.get(key)
                if prev is None or feat["score"] > prev["score"]:
                    best_by_peptide[key] = {
                        **feat,
                        "scan_id": scan_id,
                        "rt_min": rt,
                    }

    rows: list[dict] = []
    for r in targets.itertuples():
        b = best_by_peptide.get(r.fragment_id)
        if b is None:
            rows.append(
                {
                    "run_id": run.run_id,
                    "time_label": run.time_label,
                    "time_min": run.time_min,
                    "is_mapping": run.is_mapping,
                    "sequence": r.sequence,
                    "start": r.start,
                    "end": r.end,
                    "fragment_id": r.fragment_id,
                    "neutral_mass_theoretical": r.neutral_mass_theoretical,
                    "charge": np.nan,
                    "matched_isotopes": 0,
                    "score": 0.0,
                    "envelope_intensity": 0.0,
                    "neutral_mass_est": np.nan,
                    "mass_error_da": np.nan,
                    "mass_error_ppm": np.nan,
                    "scan_id": "",
                    "rt_min": np.nan,
                }
            )
            continue

        mass_error_da = b["neutral_mass_est"] - r.neutral_mass_theoretical
        mass_error_ppm = (mass_error_da / r.neutral_mass_theoretical) * 1e6
        rows.append(
            {
                "run_id": run.run_id,
                "time_label": run.time_label,
                "time_min": run.time_min,
                "is_mapping": run.is_mapping,
                "sequence": r.sequence,
                "start": r.start,
                "end": r.end,
                "fragment_id": r.fragment_id,
                "neutral_mass_theoretical": r.neutral_mass_theoretical,
                "charge": b["charge"],
                "matched_isotopes": b["matched_isotopes"],
                "score": b["score"],
                "envelope_intensity": b["envelope_intensity"],
                "neutral_mass_est": b["neutral_mass_est"],
                "mass_error_da": mass_error_da,
                "mass_error_ppm": mass_error_ppm,
                "scan_id": b["scan_id"],
                "rt_min": b["rt_min"],
            }
        )
    return rows


def finalize_metrics(df: pd.DataFrame, baseline_pattern: str, maxd_pattern: str) -> pd.DataFrame:
    out = df.copy()
    out["delta_mass_vs_baseline"] = np.nan
    out["exchange_ratio_vs_max"] = np.nan
    out["abundance_rel_to_peptide_max"] = np.nan

    for frag_id, g in out.groupby("fragment_id"):
        g = g.copy()

        baseline_rows = g[g["run_id"].str.contains(baseline_pattern, case=False, regex=True, na=False)]
        if baseline_rows.empty:
            # fallback to earliest numeric timepoint
            baseline_rows = g[g["time_min"].notna()].sort_values("time_min").head(1)
        baseline_mass = baseline_rows["neutral_mass_est"].dropna()
        baseline_mass = float(baseline_mass.iloc[0]) if not baseline_mass.empty else np.nan

        g["delta_mass_vs_baseline"] = g["neutral_mass_est"] - baseline_mass
        out.loc[g.index, "delta_mass_vs_baseline"] = g["delta_mass_vs_baseline"]

        max_ref_rows = g[g["run_id"].str.contains(maxd_pattern, case=False, regex=True, na=False)]
        if max_ref_rows.empty:
            max_ref_rows = g[g["time_min"].notna()].sort_values("time_min").tail(1)
        max_shift = max_ref_rows["delta_mass_vs_baseline"].dropna()
        max_shift = float(max_shift.iloc[0]) if not max_shift.empty else np.nan
        if np.isfinite(max_shift) and max_shift != 0:
            out.loc[g.index, "exchange_ratio_vs_max"] = g["delta_mass_vs_baseline"] / max_shift

        peak_max = g["envelope_intensity"].max()
        if np.isfinite(peak_max) and peak_max > 0:
            out.loc[g.index, "abundance_rel_to_peptide_max"] = g["envelope_intensity"] / peak_max

    return out


def build_summary(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for frag, g in df.groupby("fragment_id", sort=True):
        g = g.sort_values(["is_mapping", "time_min", "run_id"])
        non_map = g[g["is_mapping"] == False]
        if non_map.empty:
            continue
        first = non_map.sort_values("time_min").head(1).iloc[0]
        last = non_map.sort_values("time_min").tail(1).iloc[0]
        rows.append(
            {
                "fragment_id": frag,
                "sequence": first["sequence"],
                "start": int(first["start"]),
                "end": int(first["end"]),
                "early_time_min": float(first["time_min"]),
                "late_time_min": float(last["time_min"]),
                "early_abundance": float(first["envelope_intensity"]),
                "late_abundance": float(last["envelope_intensity"]),
                "early_exchange_ratio": float(first["exchange_ratio_vs_max"])
                if pd.notna(first["exchange_ratio_vs_max"])
                else np.nan,
                "late_exchange_ratio": float(last["exchange_ratio_vs_max"])
                if pd.notna(last["exchange_ratio_vs_max"])
                else np.nan,
                "max_delta_mass": float(non_map["delta_mass_vs_baseline"].max())
                if non_map["delta_mass_vs_baseline"].notna().any()
                else np.nan,
            }
        )
    return pd.DataFrame(rows).sort_values(["start", "end"]).reset_index(drop=True)


def search_run_worker(
    run: RunInfo,
    targets: pd.DataFrame,
    charge_min: int,
    charge_max: int,
    n_isotopes: int,
    min_isotopes: int,
    ppm: float,
) -> tuple[str, list[dict]]:
    rows = search_run(
        run=run,
        targets=targets,
        charge_min=charge_min,
        charge_max=charge_max,
        n_isotopes=n_isotopes,
        min_isotopes=min_isotopes,
        ppm=ppm,
    )
    return run.run_id, rows


def main() -> None:
    parser = argparse.ArgumentParser(description="Reconstruct HDX peptide metrics from mzML.")
    parser.add_argument(
        "--peptide-source",
        required=True,
        help="Peptide list source (.xlsx or .csv) with Sequence/Start/End",
    )
    parser.add_argument(
        "--mzml-dir",
        required=True,
        help="Directory containing raw-derived mzML runs",
    )
    parser.add_argument(
        "--out-dir",
        default="output/hdx_raw_reconstructed",
        help="Output directory",
    )
    parser.add_argument("--charge-min", type=int, default=1)
    parser.add_argument("--charge-max", type=int, default=8)
    parser.add_argument("--n-isotopes", type=int, default=6)
    parser.add_argument("--min-isotopes", type=int, default=2)
    parser.add_argument("--ppm", type=float, default=15.0)
    parser.add_argument("--jobs", type=int, default=1, help="Parallel workers across runs (max 8)")
    parser.add_argument("--baseline-pattern", default="peptide_mapping")
    parser.add_argument("--maxd-pattern", default="24h")
    args = parser.parse_args()

    peptide_source = Path(args.peptide_source).resolve()
    mzml_dir = Path(args.mzml_dir).resolve()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    targets = load_targets(peptide_source)
    runs = collect_runs(mzml_dir)
    workers = min(max(1, args.jobs), 8, len(runs))

    run_manifest = pd.DataFrame(
        [
            {
                "run_id": r.run_id,
                "mzml_path": str(r.mzml_path),
                "time_label": r.time_label,
                "time_min": r.time_min,
                "is_mapping": r.is_mapping,
            }
            for r in runs
        ]
    )
    run_manifest.to_csv(out_dir / "run_manifest.csv", index=False)
    targets.to_csv(out_dir / "peptide_targets.csv", index=False)

    all_rows: list[dict] = []
    if workers == 1:
        for r in runs:
            print(f"Searching run: {r.run_id}")
            rows = search_run(
                run=r,
                targets=targets,
                charge_min=args.charge_min,
                charge_max=args.charge_max,
                n_isotopes=args.n_isotopes,
                min_isotopes=args.min_isotopes,
                ppm=args.ppm,
            )
            all_rows.extend(rows)
    else:
        print(f"Searching {len(runs)} runs in parallel with {workers} workers")
        rows_by_run: dict[str, list[dict]] = {}
        with ProcessPoolExecutor(max_workers=workers) as executor:
            futures = {
                executor.submit(
                    search_run_worker,
                    r,
                    targets,
                    args.charge_min,
                    args.charge_max,
                    args.n_isotopes,
                    args.min_isotopes,
                    args.ppm,
                ): r.run_id
                for r in runs
            }
            for future in as_completed(futures):
                run_id = futures[future]
                done_run_id, rows = future.result()
                rows_by_run[done_run_id] = rows
                print(f"Finished run: {run_id}")
        for r in runs:
            all_rows.extend(rows_by_run[r.run_id])

    long_df = pd.DataFrame(all_rows)
    long_df = finalize_metrics(
        long_df,
        baseline_pattern=args.baseline_pattern,
        maxd_pattern=args.maxd_pattern,
    )
    long_df = long_df.sort_values(["start", "end", "is_mapping", "time_min", "run_id"]).reset_index(drop=True)
    long_out = out_dir / "reconstructed_peptide_timecourse.csv"
    long_df.to_csv(long_out, index=False)

    summary = build_summary(long_df)
    summary_out = out_dir / "reconstructed_peptide_summary.csv"
    summary.to_csv(summary_out, index=False)

    print(f"Wrote: {out_dir / 'run_manifest.csv'}")
    print(f"Wrote: {out_dir / 'peptide_targets.csv'}")
    print(f"Wrote: {long_out}")
    print(f"Wrote: {summary_out}")


if __name__ == "__main__":
    main()
