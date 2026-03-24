#!/usr/bin/env python3
"""Summarize fragment-table exports into compact CSV and HTML quicklooks."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import plotly.express as px


def load_fragment_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    rename = {
        "Observed Mass": "observed_mass",
        "Theoretical Mass": "theoretical_mass",
        "Start AA": "start_aa",
        "End AA": "end_aa",
        "Frag Type": "frag_type",
        "Intensity": "intensity",
        "Error": "error_ppm",
        "Sequence": "sequence",
        "Frag_number": "fragment_number",
    }
    return df.rename(columns=rename)


def build_plot(df: pd.DataFrame, out_html: Path, title: str) -> None:
    plot_df = df.copy()
    for col in ("observed_mass", "intensity", "error_ppm", "start_aa", "end_aa"):
        if col in plot_df.columns:
            plot_df[col] = pd.to_numeric(plot_df[col], errors="coerce")
    plot_df = plot_df.dropna(subset=["observed_mass", "intensity"])

    hover_cols = [col for col in ("frag_type", "sequence", "start_aa", "end_aa", "error_ppm") if col in plot_df.columns]
    color_col = "frag_type" if "frag_type" in plot_df.columns else None

    fig = px.scatter(
        plot_df,
        x="observed_mass",
        y="intensity",
        color=color_col,
        hover_data=hover_cols,
        title=title,
        template="plotly_white",
    )
    fig.update_traces(marker={"size": 8, "line": {"width": 0.5, "color": "#0f172a"}})
    fig.update_layout(
        xaxis_title="Observed Mass",
        yaxis_title="Intensity",
        font={"family": "Helvetica, Arial, sans-serif", "size": 13},
        paper_bgcolor="#f8fafc",
        plot_bgcolor="#ffffff",
        margin={"l": 70, "r": 30, "t": 60, "b": 60},
    )
    fig.write_html(out_html, include_plotlyjs="cdn")


def main() -> None:
    parser = argparse.ArgumentParser(description="Quick summary for fragment CSV exports")
    parser.add_argument("inputs", nargs="+", help="Fragment CSV file(s)")
    parser.add_argument("-o", "--out-dir", default="output/11_fragment_csv_summary")
    parser.add_argument("--top-n", type=int, default=25)
    args = parser.parse_args()

    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    summary_rows: list[dict[str, object]] = []

    for item in args.inputs:
        path = Path(item).resolve()
        df = load_fragment_table(path)
        stem = path.stem
        run_dir = out_dir / stem
        run_dir.mkdir(parents=True, exist_ok=True)

        clean_path = run_dir / f"{stem}.cleaned.csv"
        top_path = run_dir / f"{stem}.top_fragments.csv"
        html_path = run_dir / f"{stem}.fragment_map.html"

        df.to_csv(clean_path, index=False)

        if "intensity" in df.columns:
            top_df = df.sort_values("intensity", ascending=False).head(args.top_n)
        else:
            top_df = df.head(args.top_n)
        top_df.to_csv(top_path, index=False)
        build_plot(df, html_path, stem)

        summary_rows.append(
            {
                "file": str(path),
                "fragments": int(len(df)),
                "top_intensity": float(pd.to_numeric(df.get("intensity"), errors="coerce").max())
                if "intensity" in df.columns
                else None,
                "min_observed_mass": float(pd.to_numeric(df.get("observed_mass"), errors="coerce").min())
                if "observed_mass" in df.columns
                else None,
                "max_observed_mass": float(pd.to_numeric(df.get("observed_mass"), errors="coerce").max())
                if "observed_mass" in df.columns
                else None,
            }
        )

        print(f"Processed {path}")
        print(f"Wrote: {clean_path}")
        print(f"Wrote: {top_path}")
        print(f"Wrote: {html_path}")

    summary_path = out_dir / "fragment_summary.csv"
    pd.DataFrame(summary_rows).to_csv(summary_path, index=False)
    print(f"Wrote: {summary_path}")


if __name__ == "__main__":
    main()
