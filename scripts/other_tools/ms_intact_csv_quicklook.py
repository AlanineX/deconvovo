#!/usr/bin/env python3
"""Build quicklook outputs from intact-mass tabular exports."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go


MZ_CANDIDATES = ("mz", "m/z", "mass", "observed mass", "neutral mass")
INTENSITY_CANDIDATES = ("intensity", "abundance", "area", "height")


def load_table(path: Path) -> pd.DataFrame:
    head = path.read_text(encoding="utf-8", errors="ignore")[:1000].lower()
    if '#"+esi scan' in head and "#point,x(thomsons),y(counts)" in head:
        df = pd.read_csv(path, comment="#", header=None, names=["point", "mz", "intensity"])
        return df[["mz", "intensity"]]

    try:
        df = pd.read_csv(path)
    except Exception:
        df = pd.read_csv(path, header=None, engine="python")

    if df.shape[0] > 0:
        first_row = [str(value).strip().lower() for value in df.iloc[0].tolist()]
        if len(first_row) >= 2 and "x(thomsons)" in first_row[0] and "y(counts)" in first_row[1]:
            df = df.iloc[1:].reset_index(drop=True)
            df.columns = ["mz", "intensity"]
            return df

    if df.shape[1] == 2 and all(str(c).startswith("Unnamed:") or isinstance(c, int) for c in df.columns):
        df.columns = ["mz", "intensity"]
        return df

    if len(df.columns) == 2 and all(str(c).replace(".", "", 1).isdigit() for c in df.columns):
        first = str(df.columns[0])
        second = str(df.columns[1])
        data = pd.DataFrame([[float(first), float(second)]], columns=["mz", "intensity"])
        extra = df.copy()
        extra.columns = ["mz", "intensity"]
        return pd.concat([data, extra], ignore_index=True)

    rename: dict[str, str] = {}
    for col in df.columns:
        lower = str(col).strip().lower()
        if lower in MZ_CANDIDATES and "mz" not in rename.values():
            rename[col] = "mz"
        elif lower in INTENSITY_CANDIDATES and "intensity" not in rename.values():
            rename[col] = "intensity"
        elif "x(thomsons)" in lower and "mz" not in rename.values():
            rename[col] = "mz"
        elif "y(counts)" in lower and "intensity" not in rename.values():
            rename[col] = "intensity"

    df = df.rename(columns=rename)
    return df


def choose_columns(df: pd.DataFrame) -> tuple[str, str]:
    if {"mz", "intensity"}.issubset(df.columns):
        return "mz", "intensity"

    numeric_cols = [col for col in df.columns if pd.api.types.is_numeric_dtype(df[col])]
    if len(numeric_cols) >= 2:
        return numeric_cols[0], numeric_cols[1]

    raise ValueError("Could not identify m/z and intensity columns")


def build_plot(df: pd.DataFrame, mz_col: str, intensity_col: str, title: str, out_html: Path) -> None:
    plot_df = df[[mz_col, intensity_col]].dropna().copy()
    plot_df[mz_col] = pd.to_numeric(plot_df[mz_col], errors="coerce")
    plot_df[intensity_col] = pd.to_numeric(plot_df[intensity_col], errors="coerce")
    plot_df = plot_df.dropna().sort_values(mz_col)

    fig = go.Figure()
    fig.add_trace(
        go.Scattergl(
            x=plot_df[mz_col],
            y=plot_df[intensity_col],
            mode="lines",
            line={"width": 1.8, "color": "#0f766e"},
            fill="tozeroy",
            fillcolor="rgba(15,118,110,0.12)",
            hovertemplate="m/z=%{x:.5f}<br>intensity=%{y:.3e}<extra></extra>",
            name="spectrum",
        )
    )
    fig.update_layout(
        title=title,
        xaxis_title="m/z or Mass",
        yaxis_title="Intensity / Abundance",
        template="plotly_white",
        font={"family": "Helvetica, Arial, sans-serif", "size": 13},
        paper_bgcolor="#f8fafc",
        plot_bgcolor="#ffffff",
        margin={"l": 70, "r": 30, "t": 60, "b": 60},
    )
    fig.write_html(out_html, include_plotlyjs="cdn")


def main() -> None:
    parser = argparse.ArgumentParser(description="Quicklook for intact MS CSV exports")
    parser.add_argument("inputs", nargs="+", help="Input CSV file(s)")
    parser.add_argument("-o", "--out-dir", default="output/10_intact_csv_quicklook")
    parser.add_argument("--top-n", type=int, default=25)
    args = parser.parse_args()

    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    for item in args.inputs:
        path = Path(item).resolve()
        df = load_table(path)
        mz_col, intensity_col = choose_columns(df)

        work = df[[mz_col, intensity_col]].copy()
        work[mz_col] = pd.to_numeric(work[mz_col], errors="coerce")
        work[intensity_col] = pd.to_numeric(work[intensity_col], errors="coerce")
        work = work.dropna().sort_values(mz_col).reset_index(drop=True)

        stem = path.stem
        run_dir = out_dir / stem
        run_dir.mkdir(parents=True, exist_ok=True)

        cleaned_path = run_dir / f"{stem}.cleaned.csv"
        top_path = run_dir / f"{stem}.top_peaks.csv"
        html_path = run_dir / f"{stem}.quicklook.html"

        cleaned = work.rename(columns={mz_col: "mz", intensity_col: "intensity"})
        cleaned.to_csv(cleaned_path, index=False)
        cleaned.nlargest(args.top_n, "intensity").sort_values("intensity", ascending=False).to_csv(top_path, index=False)
        build_plot(cleaned, "mz", "intensity", stem, html_path)

        print(f"Processed {path}")
        print(f"Wrote: {cleaned_path}")
        print(f"Wrote: {top_path}")
        print(f"Wrote: {html_path}")


if __name__ == "__main__":
    main()
