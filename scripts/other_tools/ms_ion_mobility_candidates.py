#!/usr/bin/env python3
"""Filter an MS inventory down to ion-mobility and calibrant candidates."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract ion-mobility candidates from ms_inventory.csv")
    parser.add_argument("--inventory", required=True, help="Path to ms_inventory.csv")
    parser.add_argument("-o", "--out-dir", default="output/12_ion_mobility_candidates")
    args = parser.parse_args()

    inventory = Path(args.inventory).resolve()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(inventory)
    subset = df[(df["ion_mobility"].fillna(False)) | (df["calibrant"].fillna(False))].copy()
    subset = subset.sort_values(["vendor", "workflow_hint", "relative_path"]).reset_index(drop=True)

    csv_path = out_dir / "ion_mobility_candidates.csv"
    md_path = out_dir / "ion_mobility_candidates.md"

    subset.to_csv(csv_path, index=False)

    lines = ["# Ion Mobility Candidates", ""]
    lines.append(f"- Candidate artifacts: {len(subset)}")
    lines.append("")
    if subset.empty:
        lines.append("No ion-mobility or calibrant candidates detected.")
    else:
        lines.append(
            subset[
                [
                    "relative_path",
                    "vendor",
                    "workflow_hint",
                    "analyte_hint",
                    "ion_mobility",
                    "calibrant",
                    "notes",
                ]
            ].to_csv(index=False)
        )
    lines.append("")

    md_path.write_text("\n".join(lines), encoding="utf-8")

    print(f"Wrote: {csv_path}")
    print(f"Wrote: {md_path}")


if __name__ == "__main__":
    main()
