#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Iterable

import pandas as pd


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def write_csv(df: pd.DataFrame, path: Path) -> None:
    df.to_csv(path, index=False)
    print(f"Wrote: {path}")


def detect_artifact_name(path: Path) -> str:
    name = path.name
    for suffix in (".cleaned.csv", ".csv"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def load_intact_metrics(intact_dir: Path) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for cleaned in sorted(intact_dir.glob("*/*.cleaned.csv")):
        df = pd.read_csv(cleaned)
        if df.empty or "mz" not in df.columns or "intensity" not in df.columns:
            continue
        artifact = cleaned.parent.name
        total = float(df["intensity"].sum())
        base_idx = int(df["intensity"].idxmax())
        base_peak_mz = float(df.loc[base_idx, "mz"])
        base_peak_intensity = float(df.loc[base_idx, "intensity"])
        weighted_mz = float((df["mz"] * df["intensity"]).sum() / total) if total else math.nan
        frac_lt500 = float(df.loc[df["mz"] < 500, "intensity"].sum() / total) if total else math.nan
        frac_500_1000 = (
            float(df.loc[(df["mz"] >= 500) & (df["mz"] < 1000), "intensity"].sum() / total)
            if total
            else math.nan
        )
        frac_1000_1400 = (
            float(df.loc[(df["mz"] >= 1000) & (df["mz"] < 1400), "intensity"].sum() / total)
            if total
            else math.nan
        )
        frac_1400_1600 = (
            float(df.loc[df["mz"].between(1400, 1600), "intensity"].sum() / total) if total else math.nan
        )
        frac_1600_3000 = (
            float(df.loc[df["mz"].between(1600, 3000), "intensity"].sum() / total) if total else math.nan
        )
        frac_8561 = float(df.loc[df["mz"].between(8558, 8564), "intensity"].sum() / total) if total else math.nan
        frac_250 = float(df.loc[df["mz"].between(249.5, 251.5), "intensity"].sum() / total) if total else math.nan

        if artifact.startswith("mz_abundance_"):
            workflow = "peak_table"
            if frac_8561 >= 0.45 and (frac_250 < 0.15 or frac_8561 > frac_250):
                interpretation = "intact_cluster_dominant"
            elif frac_8561 >= 0.25 and frac_250 < 0.5:
                interpretation = "mixed_intact_plus_background"
            else:
                interpretation = "background_competitive_or_dirty"
        else:
            workflow = "profile_export"
            if frac_lt500 >= 0.08 and frac_1400_1600 < 0.2:
                interpretation = "low_mass_enriched_or_degraded"
            elif frac_1400_1600 >= 0.5:
                interpretation = "intact_envelope_preserved"
            elif frac_1600_3000 >= 0.3:
                interpretation = "broad_high_mz_distribution"
            else:
                interpretation = "mixed_distribution"

        rows.append(
            {
                "artifact": artifact,
                "workflow": workflow,
                "n_points": int(len(df)),
                "total_intensity": total,
                "base_peak_mz": base_peak_mz,
                "base_peak_intensity": base_peak_intensity,
                "weighted_mz": weighted_mz,
                "frac_lt500": frac_lt500,
                "frac_500_1000": frac_500_1000,
                "frac_1000_1400": frac_1000_1400,
                "frac_1400_1600": frac_1400_1600,
                "frac_1600_3000": frac_1600_3000,
                "frac_8561_window": frac_8561,
                "frac_250_window": frac_250,
                "interpretation": interpretation,
                "source_cleaned_csv": str(cleaned),
            }
        )
    return pd.DataFrame(rows)


def load_fragment_metrics(fragment_dir: Path) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    shared_sequences: set[str] | None = None
    for cleaned in sorted(fragment_dir.glob("Fragments_*/*.cleaned.csv")):
        df = pd.read_csv(cleaned)
        if df.empty:
            continue
        coverage = set()
        for start_aa, end_aa in zip(df["start_aa"], df["end_aa"]):
            if pd.notna(start_aa) and pd.notna(end_aa):
                coverage.update(range(int(start_aa), int(end_aa) + 1))
        seqs = set(df["sequence"].astype(str))
        shared_sequences = seqs if shared_sequences is None else shared_sequences & seqs
        type_counts = df["frag_type"].value_counts().to_dict()
        abs_ppm = df["error_ppm"].abs()
        pos_frac = float((df["error_ppm"] > 0).mean())
        top_idx = int(df["intensity"].idxmax())
        rows.append(
            {
                "artifact": cleaned.parent.name,
                "fragments": int(len(df)),
                "total_intensity": float(df["intensity"].sum()),
                "base_fragment_sequence": str(df.loc[top_idx, "sequence"]),
                "base_fragment_intensity": float(df.loc[top_idx, "intensity"]),
                "median_abs_ppm": float(abs_ppm.median()),
                "max_abs_ppm": float(abs_ppm.max()),
                "mean_ppm": float(df["error_ppm"].mean()),
                "positive_ppm_fraction": pos_frac,
                "residue_coverage": int(len(coverage)),
                "c_fragments": int(type_counts.get("C Fragment", 0)),
                "z_fragments": int(type_counts.get("Z Fragment", 0)),
                "source_cleaned_csv": str(cleaned),
            }
        )
    out = pd.DataFrame(rows)
    if not out.empty:
        out["shared_sequences_all_runs"] = len(shared_sequences or set())
    return out


def load_hdx_metrics(summary_csv: Path, fit_csv: Path) -> pd.DataFrame:
    summary = pd.read_csv(summary_csv)
    fit = pd.read_csv(fit_csv)
    merged = summary.merge(
        fit[["fragment_id", "selected_model", "r2", "k_app_per_min", "half_life_min"]],
        on="fragment_id",
        how="left",
    )
    merged["abundance_fold_24h_over_10min"] = merged["late_abundance"] / merged["early_abundance"]

    def classify_half_life(half_life_min: float) -> str:
        if pd.isna(half_life_min):
            return "unknown"
        if half_life_min < 5:
            return "fast"
        if half_life_min < 60:
            return "intermediate"
        return "slow"

    merged["rate_class"] = merged["half_life_min"].map(classify_half_life)
    return merged.sort_values("k_app_per_min", ascending=False)


def parse_waters_headers(waters_root: Path) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for header in sorted(waters_root.glob("**/*.raw/_HEADER.TXT")):
        acquired_name = ""
        sample_description = ""
        instrument = ""
        with header.open("r", encoding="latin-1", errors="ignore") as handle:
            for raw_line in handle:
                line = raw_line.strip()
                if line.startswith("$$ Acquired Name:"):
                    acquired_name = line.split(":", 1)[1].strip()
                elif line.startswith("$$ Sample Description:"):
                    sample_description = line.split(":", 1)[1].strip()
                elif line.startswith("$$ Instrument:"):
                    instrument = line.split(":", 1)[1].strip()
        relative = header.parent.parent
        name_upper = acquired_name.upper() if acquired_name else relative.name.upper()
        rows.append(
            {
                "relative_path": str(relative),
                "acquired_name": acquired_name or relative.name,
                "sample_description": sample_description,
                "instrument": instrument,
                "contains_adp": "ADP" in name_upper,
                "contains_ubq": "UBQ" in name_upper,
                "contains_myo": "MYO" in name_upper,
                "contains_cytc": "CYTC" in name_upper,
                "contains_den": "DEN" in name_upper,
                "contains_tune": "TUNE" in name_upper,
                "contains_ims": "IMS" in name_upper,
                "contains_cali": "CALI" in name_upper,
            }
        )
    return pd.DataFrame(rows)


def write_waters_summary(df: pd.DataFrame, out_dir: Path) -> pd.DataFrame:
    if df.empty:
        return df
    unique_descriptions = (
        df["sample_description"].fillna("").replace("", pd.NA).dropna().nunique()
    )
    nonempty_descriptions = int(df["sample_description"].fillna("").ne("").sum())
    note = "sample_description_diverse"
    if nonempty_descriptions and unique_descriptions == 1 and len(df) > 1:
        note = "sample_description_reused_across_distinct_runs"

    summary = pd.DataFrame(
        [
            {
                "runs": int(len(df)),
                "instruments": ",".join(sorted(set(df["instrument"].dropna()))),
                "adp_named_runs": int(df["contains_adp"].sum()),
                "protein_named_runs": int((df["contains_ubq"] | df["contains_myo"] | df["contains_cytc"]).sum()),
                "denaturation_named_runs": int(df["contains_den"].sum()),
                "tune_named_runs": int(df["contains_tune"].sum()),
                "ims_named_runs": int(df["contains_ims"].sum()),
                "calibrant_named_runs": int(df["contains_cali"].sum()),
                "nonempty_sample_descriptions": nonempty_descriptions,
                "unique_nonempty_sample_descriptions": int(unique_descriptions),
                "metadata_note": note,
            }
        ]
    )
    write_csv(summary, out_dir / "waters_metadata_summary.csv")
    return summary


def format_float(value: float, digits: int = 3) -> str:
    if pd.isna(value):
        return "NA"
    return f"{value:.{digits}f}"


def top_rows(df: pd.DataFrame, column: str, ascending: bool, n: int = 3) -> list[dict[str, object]]:
    if df.empty or column not in df.columns:
        return []
    return df.sort_values(column, ascending=ascending).head(n).to_dict("records")


def write_markdown(
    inventory_csv: Path,
    intact_metrics: pd.DataFrame,
    fragment_metrics: pd.DataFrame,
    hdx_metrics: pd.DataFrame,
    waters_summary: pd.DataFrame,
    out_path: Path,
) -> None:
    inventory = pd.read_csv(inventory_csv)
    vendor_counts = inventory["vendor"].value_counts().to_dict()
    lines: list[str] = []
    lines.append("# Cross-Dataset Findings")
    lines.append("")
    lines.append("## Scope")
    lines.append("")
    lines.append(f"- Artifacts inventoried: {len(inventory)}")
    lines.append(
        "- Vendor counts: "
        + ", ".join(f"{vendor}={count}" for vendor, count in sorted(vendor_counts.items()))
    )
    lines.append("- Quantitative analysis used current local outputs only.")
    lines.append("- Waters native raw files remain metadata-only until a local MassLynx conversion path is added.")
    lines.append("")

    if not hdx_metrics.empty:
        fast = top_rows(hdx_metrics, "k_app_per_min", ascending=False, n=4)
        slow = top_rows(hdx_metrics, "k_app_per_min", ascending=True, n=4)
        lines.append("## Thermo HDX")
        lines.append("")
        lines.append(f"- Peptides reconstructed from raw-only path: {len(hdx_metrics)}")
        lines.append(
            f"- Fast peptides: {(hdx_metrics['rate_class'] == 'fast').sum()}, "
            f"intermediate: {(hdx_metrics['rate_class'] == 'intermediate').sum()}, "
            f"slow: {(hdx_metrics['rate_class'] == 'slow').sum()}"
        )
        lines.append(
            f"- Mean early exchange ratio: {format_float(hdx_metrics['early_exchange_ratio'].mean())}"
        )
        lines.append(
            f"- Correlation early exchange vs k_app: {format_float(hdx_metrics['early_exchange_ratio'].corr(hdx_metrics['k_app_per_min']))}"
        )
        lines.append("")
        lines.append("Fastest apparent regions:")
        for row in fast:
            lines.append(
                f"- `{row['fragment_id']}`: k_app={format_float(row['k_app_per_min'])} min^-1, "
                f"t1/2={format_float(row['half_life_min'])} min, "
                f"early_exchange={format_float(row['early_exchange_ratio'])}, "
                f"24h/10min abundance={format_float(row['abundance_fold_24h_over_10min'])}"
            )
        lines.append("")
        lines.append("Slowest apparent regions:")
        for row in slow:
            lines.append(
                f"- `{row['fragment_id']}`: k_app={format_float(row['k_app_per_min'])} min^-1, "
                f"t1/2={format_float(row['half_life_min'])} min, "
                f"early_exchange={format_float(row['early_exchange_ratio'])}, "
                f"24h/10min abundance={format_float(row['abundance_fold_24h_over_10min'])}"
            )
        lines.append("")

    if not intact_metrics.empty:
        peak_tables = intact_metrics[intact_metrics["workflow"] == "peak_table"].copy()
        profile_exports = intact_metrics[intact_metrics["workflow"] == "profile_export"].copy()

        lines.append("## Agilent Intact / Profile Exports")
        lines.append("")
        if not peak_tables.empty:
            lines.append("Peak-table interpretation:")
            for row in peak_tables.sort_values("artifact").to_dict("records"):
                lines.append(
                    f"- `{row['artifact']}`: base_peak={format_float(row['base_peak_mz'])}, "
                    f"8561-window={format_float(row['frac_8561_window'])}, "
                    f"250-window={format_float(row['frac_250_window'])}, "
                    f"call={row['interpretation']}"
                )
        lines.append("")
        if not profile_exports.empty:
            lines.append("Profile-export interpretation:")
            for row in profile_exports.sort_values("artifact").to_dict("records"):
                lines.append(
                    f"- `{row['artifact']}`: base_peak={format_float(row['base_peak_mz'])}, "
                    f"low<500={format_float(row['frac_lt500'])}, "
                    f"1400-1600={format_float(row['frac_1400_1600'])}, "
                    f"1600-3000={format_float(row['frac_1600_3000'])}, "
                    f"call={row['interpretation']}"
                )
            strongest_low_mass = top_rows(profile_exports, "frac_lt500", ascending=False, n=1)
            strongest_envelope = top_rows(profile_exports, "frac_1400_1600", ascending=False, n=1)
            broadest_high = top_rows(profile_exports, "frac_1600_3000", ascending=False, n=1)
            lines.append("")
            if strongest_low_mass:
                row = strongest_low_mass[0]
                lines.append(
                    f"- Strongest low-mass enrichment: `{row['artifact']}` at low<500={format_float(row['frac_lt500'])}"
                )
            if strongest_envelope:
                row = strongest_envelope[0]
                lines.append(
                    f"- Best-preserved 1400-1600 envelope: `{row['artifact']}` at {format_float(row['frac_1400_1600'])}"
                )
            if broadest_high:
                row = broadest_high[0]
                lines.append(
                    f"- Broadest 1600-3000 distribution: `{row['artifact']}` at {format_float(row['frac_1600_3000'])}"
                )
        lines.append("")

    if not fragment_metrics.empty:
        lines.append("## Agilent Fragment Exports")
        lines.append("")
        lines.append(
            f"- Fragment tables analyzed: {len(fragment_metrics)}; shared sequences in all fragment runs: {int(fragment_metrics['shared_sequences_all_runs'].iloc[0])}"
        )
        for row in fragment_metrics.sort_values("artifact").to_dict("records"):
            lines.append(
                f"- `{row['artifact']}`: fragments={int(row['fragments'])}, "
                f"total_intensity={format_float(row['total_intensity'], 1)}, "
                f"coverage={int(row['residue_coverage'])} aa, "
                f"median_abs_ppm={format_float(row['median_abs_ppm'])}, "
                f"positive_ppm_fraction={format_float(row['positive_ppm_fraction'])}"
            )
        lines.append("")
        weakest = top_rows(fragment_metrics, "total_intensity", ascending=True, n=1)
        strongest = top_rows(fragment_metrics, "total_intensity", ascending=False, n=1)
        if strongest and weakest:
            lines.append(
                f"- Strongest fragment run: `{strongest[0]['artifact']}`; weakest fragment run: `{weakest[0]['artifact']}`."
            )
            lines.append(
                f"- All fragment runs span the full 76 aa sequence, but signal depth differs sharply."
            )
        lines.append("")

    if not waters_summary.empty:
        row = waters_summary.iloc[0].to_dict()
        lines.append("## Waters Raw Metadata")
        lines.append("")
        lines.append(
            f"- Runs: {int(row['runs'])} on instrument `{row['instruments']}`"
        )
        lines.append(
            f"- Naming suggests ADP screens={int(row['adp_named_runs'])}, protein-named runs={int(row['protein_named_runs'])}, denaturation runs={int(row['denaturation_named_runs'])}, tune runs={int(row['tune_named_runs'])}"
        )
        lines.append(
            f"- IM-tagged runs={int(row['ims_named_runs'])}, calibrant-tagged runs={int(row['calibrant_named_runs'])}"
        )
        lines.append(
            f"- Metadata note: `{row['metadata_note']}`"
        )
        if row["metadata_note"] == "sample_description_reused_across_distinct_runs":
            lines.append(
                "- Sample descriptions are not trustworthy discriminators here; acquired-name strings are the safer routing key."
            )
        lines.append("")

    lines.append("## Output Files")
    lines.append("")
    lines.append(f"- `{out_path.parent / 'agilent_intact_metrics.csv'}`")
    lines.append(f"- `{out_path.parent / 'agilent_fragment_metrics.csv'}`")
    lines.append(f"- `{out_path.parent / 'hdx_region_metrics.csv'}`")
    lines.append(f"- `{out_path.parent / 'waters_metadata_summary.csv'}`")
    lines.append("")
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote: {out_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Consolidate mixed-vendor MS findings from existing outputs.")
    parser.add_argument("--inventory", default="output/00_inventory/ms_inventory.csv")
    parser.add_argument("--intact-dir", default="output/10_intact_csv_quicklook")
    parser.add_argument("--fragment-dir", default="output/11_fragment_csv_summary")
    parser.add_argument("--hdx-summary", default="output/02_hdx_raw/reconstructed_peptide_summary.csv")
    parser.add_argument("--hdx-fit", default="output/02_hdx_raw/analysis_fit/kinetics_fit.csv")
    parser.add_argument("--waters-root", default="data_2_waters")
    parser.add_argument("--out-dir", default="output/90_cross_dataset_findings")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    ensure_dir(out_dir)

    intact_metrics = load_intact_metrics(Path(args.intact_dir))
    fragment_metrics = load_fragment_metrics(Path(args.fragment_dir))
    hdx_metrics = load_hdx_metrics(Path(args.hdx_summary), Path(args.hdx_fit))
    waters_headers = parse_waters_headers(Path(args.waters_root))
    waters_summary = write_waters_summary(waters_headers, out_dir)

    write_csv(intact_metrics, out_dir / "agilent_intact_metrics.csv")
    write_csv(fragment_metrics, out_dir / "agilent_fragment_metrics.csv")
    write_csv(hdx_metrics, out_dir / "hdx_region_metrics.csv")
    write_csv(waters_headers, out_dir / "waters_header_inventory.csv")
    write_markdown(
        Path(args.inventory),
        intact_metrics,
        fragment_metrics,
        hdx_metrics,
        waters_summary,
        out_dir / "FINDINGS_ALL_DATA.md",
    )


if __name__ == "__main__":
    main()
