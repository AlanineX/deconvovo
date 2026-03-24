#!/usr/bin/env python3
"""Analyze HDX peptide timecourse CSV and export reusable summary artifacts."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go


def estimate_t50(time_min: np.ndarray, uptake_pct: np.ndarray) -> float | None:
    if len(time_min) < 2:
        return None
    y_max = np.nanmax(uptake_pct)
    if not np.isfinite(y_max) or y_max <= 0:
        return None
    target = 0.5 * y_max
    for i in range(1, len(time_min)):
        y0, y1 = uptake_pct[i - 1], uptake_pct[i]
        t0, t1 = time_min[i - 1], time_min[i]
        if not np.isfinite(y0) or not np.isfinite(y1):
            continue
        if (y0 <= target <= y1) or (y1 <= target <= y0):
            if y1 == y0:
                return float(t0)
            frac = (target - y0) / (y1 - y0)
            return float(t0 + frac * (t1 - t0))
    return None


def trapezoid_auc(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 2:
        return float("nan")
    return float(np.trapezoid(y, x))


def build_line_plot(df: pd.DataFrame, out_html: Path) -> None:
    fig = px.line(
        df,
        x="time_min",
        y="deuterium_uptake_pct",
        color="fragment_id",
        markers=True,
        template="plotly_white",
        title="HDX Timecourse by Peptide",
    )
    fig.update_layout(
        xaxis_title="Time (min)",
        yaxis_title="Deuterium Uptake (%)",
        legend_title="Peptide",
        font={"family": "Helvetica, Arial, sans-serif", "size": 12},
        paper_bgcolor="#f7fafc",
        plot_bgcolor="#ffffff",
    )
    fig.write_html(out_html, include_plotlyjs="cdn")


def build_heatmap(df: pd.DataFrame, out_html: Path) -> None:
    pivot = (
        df.pivot_table(
            index="fragment_id",
            columns="time_min",
            values="deuterium_uptake_pct",
            aggfunc="mean",
        )
        .sort_index()
        .sort_index(axis=1)
    )
    x = [f"{c:g}" for c in pivot.columns]
    y = pivot.index.tolist()
    z = pivot.values

    fig = go.Figure(
        data=[
            go.Heatmap(
                x=x,
                y=y,
                z=z,
                colorscale="YlOrRd",
                colorbar={"title": "Uptake (%)"},
                hovertemplate="Peptide: %{y}<br>Time: %{x} min<br>Uptake: %{z:.2f}%<extra></extra>",
            )
        ]
    )
    fig.update_layout(
        title="HDX Uptake Heatmap",
        xaxis_title="Time (min)",
        yaxis_title="Peptide",
        template="plotly_white",
        font={"family": "Helvetica, Arial, sans-serif", "size": 12},
        paper_bgcolor="#f7fafc",
        plot_bgcolor="#ffffff",
    )
    fig.write_html(out_html, include_plotlyjs="cdn")


def main() -> None:
    parser = argparse.ArgumentParser(description="Analyze peptide-level HDX timecourse CSV.")
    parser.add_argument(
        "--input",
        required=True,
        help="Input peptide_timecourse_long.csv path",
    )
    parser.add_argument(
        "--out-dir",
        default="output/hdx_analysis",
        help="Output directory (default: output/hdx_analysis)",
    )
    args = parser.parse_args()

    in_path = Path(args.input).resolve()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(in_path)
    required_cols = {
        "fragment_id",
        "sequence",
        "start",
        "end",
        "time_min",
        "deuterium_uptake_pct",
        "exchange_ratio",
        "mass_uptake",
    }
    missing = required_cols - set(df.columns)
    if missing:
        raise SystemExit(f"Missing required columns: {sorted(missing)}")

    df = df.sort_values(["fragment_id", "time_min"]).reset_index(drop=True)
    time_min = float(df["time_min"].min())
    time_max = float(df["time_min"].max())

    per_rows = []
    for frag_id, g in df.groupby("fragment_id", sort=True):
        g = g.sort_values("time_min")
        t = g["time_min"].to_numpy(dtype=float)
        u = g["deuterium_uptake_pct"].to_numpy(dtype=float)

        early_row = g.iloc[0]
        late_row = g.iloc[-1]
        is_monotonic = bool(np.all(np.diff(u) >= -1e-9))
        per_rows.append(
            {
                "fragment_id": frag_id,
                "sequence": str(early_row["sequence"]),
                "start": int(early_row["start"]),
                "end": int(early_row["end"]),
                "uptake_pct_early": float(early_row["deuterium_uptake_pct"]),
                "uptake_pct_late": float(late_row["deuterium_uptake_pct"]),
                "exchange_ratio_early": float(early_row["exchange_ratio"]),
                "exchange_ratio_late": float(late_row["exchange_ratio"]),
                "mass_uptake_early": float(early_row["mass_uptake"]),
                "mass_uptake_late": float(late_row["mass_uptake"]),
                "uptake_delta_pct": float(np.nanmax(u) - np.nanmin(u)),
                "uptake_auc_pct_min": trapezoid_auc(t, u),
                "t50_min_est": estimate_t50(t, u),
                "is_monotonic": is_monotonic,
            }
        )

    per_df = pd.DataFrame(per_rows).sort_values(["start", "end"]).reset_index(drop=True)
    per_out = out_dir / "peptide_summary.csv"
    per_df.to_csv(per_out, index=False)

    top_fast = per_df.sort_values("uptake_pct_early", ascending=False).head(5)
    top_slow = per_df.sort_values("uptake_pct_early", ascending=True).head(5)
    top_dynamic = per_df.sort_values("uptake_delta_pct", ascending=False).head(5)

    max_pos = int(per_df["end"].max())
    coverage = np.zeros(max_pos + 1, dtype=int)
    for _, r in per_df.iterrows():
        coverage[int(r["start"]) : int(r["end"]) + 1] += 1
    covered_positions = int((coverage[1:] > 0).sum())

    report_lines = [
        "# HDX Analysis Summary",
        "",
        f"- Input: `{in_path}`",
        f"- Peptides: {per_df.shape[0]}",
        f"- Timepoints: {sorted(df['time_min'].unique().tolist())}",
        f"- Early timepoint: {time_min:g} min",
        f"- Late timepoint: {time_max:g} min",
        f"- Sequence coverage (from peptide intervals): {covered_positions}/{max_pos} ({covered_positions/max_pos:.2%})",
        f"- Monotonic peptides: {int(per_df['is_monotonic'].sum())}/{per_df.shape[0]}",
        "",
        "## Top Fast-Exchanging Peptides (Earliest Timepoint)",
        "",
    ]
    report_lines.extend(
        [f"- {r.fragment_id}: {r.uptake_pct_early:.2f}% at {time_min:g} min" for r in top_fast.itertuples()]
    )
    report_lines.extend(["", "## Slowest Early-Exchanging Peptides", ""])
    report_lines.extend(
        [f"- {r.fragment_id}: {r.uptake_pct_early:.2f}% at {time_min:g} min" for r in top_slow.itertuples()]
    )
    report_lines.extend(["", "## Most Dynamic Peptides (Max-Min Uptake)", ""])
    report_lines.extend(
        [f"- {r.fragment_id}: Δ{r.uptake_delta_pct:.2f}% uptake" for r in top_dynamic.itertuples()]
    )

    report_path = out_dir / "summary.md"
    report_path.write_text("\n".join(report_lines) + "\n", encoding="utf-8")

    build_line_plot(df, out_dir / "kinetics_lines.html")
    build_heatmap(df, out_dir / "uptake_heatmap.html")

    print(f"Wrote: {per_out}")
    print(f"Wrote: {report_path}")
    print(f"Wrote: {out_dir / 'kinetics_lines.html'}")
    print(f"Wrote: {out_dir / 'uptake_heatmap.html'}")


if __name__ == "__main__":
    main()
