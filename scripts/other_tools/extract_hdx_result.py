#!/usr/bin/env python3
"""Extract reusable peptide-level HDX tables from Result.xlsx-style exports."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd


def norm_text(value: object) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip().lower()


def parse_time_to_min(value: object) -> float | None:
    if pd.isna(value):
        return None
    if isinstance(value, (int, float)):
        return float(value)

    s = str(value).strip().lower().replace(" ", "")
    if not s:
        return None
    if re.fullmatch(r"\d+(\.\d+)?", s):
        return float(s)

    m = re.fullmatch(r"(\d+(?:\.\d+)?)(s|sec|secs|second|seconds|min|mins|m|h|hr|hrs|hour|hours|d|day|days)", s)
    if not m:
        return None
    num = float(m.group(1))
    unit = m.group(2)
    if unit in {"s", "sec", "secs", "second", "seconds"}:
        return num / 60.0
    if unit in {"min", "mins", "m"}:
        return num
    if unit in {"h", "hr", "hrs", "hour", "hours"}:
        return num * 60.0
    if unit in {"d", "day", "days"}:
        return num * 1440.0
    return None


def find_col(top_headers: list[str], name: str) -> int:
    target = name.strip().lower()
    for i, col in enumerate(top_headers):
        if col == target:
            return i
    raise ValueError(f"Required column '{name}' not found in first header row.")


def to_numeric_safe(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def parse_sheet(df_raw: pd.DataFrame, sequence_contains: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    if df_raw.shape[0] < 3:
        raise ValueError("Sheet has fewer than 3 rows; cannot parse two-row header + data.")

    header_top = [norm_text(v) for v in df_raw.iloc[0].tolist()]
    header_sub_raw = df_raw.iloc[1].tolist()
    data = df_raw.iloc[2:].reset_index(drop=True).copy()

    seq_idx = find_col(header_top, "sequence")
    start_idx = find_col(header_top, "start")
    end_idx = find_col(header_top, "end")
    maxnd_idx = find_col(header_top, "maxnd")
    maxd_idx = find_col(header_top, "maxd")
    delta_mass_idx = find_col(header_top, "delta mass")

    deuterium_heading_indices = [i for i, h in enumerate(header_top) if "deuterium uptake" in h]
    if not deuterium_heading_indices:
        raise ValueError("Could not locate 'Deuterium uptake (%)' heading.")
    deuterium_start_idx = min(deuterium_heading_indices)

    # Mass timepoints are between Delta Mass and Deuterium uptake blocks.
    mass_time_cols: list[int] = []
    for i in range(delta_mass_idx + 1, deuterium_start_idx):
        t = parse_time_to_min(header_sub_raw[i])
        if t is not None:
            mass_time_cols.append(i)

    # Percent uptake columns are in/after the Deuterium uptake block.
    pct_time_cols: list[int] = []
    for i in range(deuterium_start_idx, len(header_sub_raw)):
        t = parse_time_to_min(header_sub_raw[i])
        if t is not None:
            pct_time_cols.append(i)

    if not mass_time_cols:
        raise ValueError("No mass uptake timepoint columns detected.")
    if not pct_time_cols:
        raise ValueError("No deuterium uptake (%) timepoint columns detected.")

    seq = data.iloc[:, seq_idx].astype(str)
    valid_rows = seq.str.strip().ne("") & data.iloc[:, seq_idx].notna()
    data = data.loc[valid_rows].reset_index(drop=True)

    if sequence_contains:
        pattern = "|".join(re.escape(x) for x in sequence_contains if x.strip())
        if pattern:
            keep = data.iloc[:, seq_idx].astype(str).str.contains(pattern, case=False, regex=True, na=False)
            data = data.loc[keep].reset_index(drop=True)

    wide_rows = []
    long_rows = []

    pct_col_by_min = {parse_time_to_min(header_sub_raw[i]): i for i in pct_time_cols}

    for _, row in data.iterrows():
        sequence = str(row.iloc[seq_idx]).strip()
        start = row.iloc[start_idx]
        end = row.iloc[end_idx]
        maxnd = row.iloc[maxnd_idx]
        maxd = row.iloc[maxd_idx]
        delta_mass = row.iloc[delta_mass_idx]

        base = {
            "sequence": sequence,
            "start": start,
            "end": end,
            "fragment_id": f"{sequence}|{start}-{end}",
            "maxND": maxnd,
            "maxD": maxd,
            "delta_mass": delta_mass,
        }

        wide_row = dict(base)

        for col in mass_time_cols:
            time_label = header_sub_raw[col]
            time_min = parse_time_to_min(time_label)
            mass_uptake = row.iloc[col]
            pct_col = pct_col_by_min.get(time_min)
            pct_uptake = row.iloc[pct_col] if pct_col is not None else pd.NA

            wide_row[f"mass_uptake_{time_min:g}min"] = mass_uptake
            wide_row[f"deuterium_uptake_pct_{time_min:g}min"] = pct_uptake

            maxnd_num = pd.to_numeric(pd.Series([maxnd]), errors="coerce").iloc[0]
            mass_num = pd.to_numeric(pd.Series([mass_uptake]), errors="coerce").iloc[0]
            pct_num = pd.to_numeric(pd.Series([pct_uptake]), errors="coerce").iloc[0]

            normalized_to_maxnd = pd.NA
            if pd.notna(maxnd_num) and maxnd_num != 0 and pd.notna(mass_num):
                normalized_to_maxnd = mass_num / maxnd_num

            exchange_ratio = pd.NA
            if pd.notna(pct_num):
                exchange_ratio = pct_num / 100.0

            long_rows.append(
                {
                    **base,
                    "time_label": time_label,
                    "time_min": time_min,
                    "mass_uptake": mass_uptake,
                    "deuterium_uptake_pct": pct_uptake,
                    "exchange_ratio": exchange_ratio,
                    "normalized_mass_to_maxND": normalized_to_maxnd,
                }
            )

        wide_rows.append(wide_row)

    wide_df = pd.DataFrame(wide_rows)
    long_df = pd.DataFrame(long_rows)

    for col in ("start", "end", "maxND", "maxD", "delta_mass"):
        if col in wide_df.columns:
            wide_df[col] = to_numeric_safe(wide_df[col])
        if col in long_df.columns:
            long_df[col] = to_numeric_safe(long_df[col])

    sort_cols = [c for c in ["start", "end", "sequence"] if c in wide_df.columns]
    if sort_cols:
        wide_df = wide_df.sort_values(sort_cols).reset_index(drop=True)
        long_df = long_df.sort_values(sort_cols + ["time_min"]).reset_index(drop=True)

    return wide_df, long_df


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract HDX peptide wide/long tables from Result.xlsx")
    parser.add_argument("--input", required=True, help="Input .xlsx path")
    parser.add_argument("--sheet", default=None, help="Sheet name (default: first sheet)")
    parser.add_argument(
        "--out-dir",
        default="output/hdx",
        help="Output directory (default: output/hdx)",
    )
    parser.add_argument(
        "--sequence-contains",
        action="append",
        default=[],
        help="Optional case-insensitive peptide sequence substring filter (can repeat)",
    )
    args = parser.parse_args()

    in_path = Path(args.input).resolve()
    if not in_path.exists():
        raise SystemExit(f"Input file not found: {in_path}")

    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    xls = pd.ExcelFile(in_path, engine="openpyxl")
    sheet = args.sheet or xls.sheet_names[0]
    if sheet not in xls.sheet_names:
        raise SystemExit(f"Sheet '{sheet}' not found. Available: {xls.sheet_names}")

    df_raw = pd.read_excel(in_path, sheet_name=sheet, header=None, engine="openpyxl")
    wide_df, long_df = parse_sheet(df_raw, args.sequence_contains)

    stem = in_path.stem
    filter_tag = ""
    if args.sequence_contains:
        safe = "_".join(re.sub(r"[^A-Za-z0-9]+", "-", s.strip()) for s in args.sequence_contains if s.strip())
        if safe:
            filter_tag = f".filter-{safe}"

    wide_out = out_dir / f"{stem}{filter_tag}.peptide_wide.csv"
    long_out = out_dir / f"{stem}{filter_tag}.peptide_timecourse_long.csv"

    wide_df.to_csv(wide_out, index=False)
    long_df.to_csv(long_out, index=False)

    print(f"Wrote: {wide_out}")
    print(f"Wrote: {long_out}")
    print(f"Peptides: {len(wide_df)}")
    print(f"Timepoint rows: {len(long_df)}")


if __name__ == "__main__":
    main()
