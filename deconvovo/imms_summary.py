"""Step 4: Generate summary CSV from deconvolution results."""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def run(results: list[dict], out_dir: Path) -> None:
    """Write summary CSV from deconv result dicts."""
    if not results:
        print("  No results to summarize")
        return

    rows = []
    for r in results:
        row = {
            "run_name": r.get("run_name"),
            "n_peaks": r.get("n_peaks", 0),
        }
        peaks = r.get("peaks", [])
        if peaks:
            top = sorted(peaks, key=lambda p: -p["height"])[:3]
            row["top_masses"] = "; ".join(f"{p['mass']:.1f}" for p in top)
        rows.append(row)

    out_path = out_dir / "unidec_summary.csv"
    pd.DataFrame(rows).to_csv(out_path, index=False)
    print(f"  Wrote {out_path}")


def run_from_dir(out_dir: Path) -> None:
    """Rebuild summary by scanning _peaks.csv files in output dir."""
    peaks_files = sorted(out_dir.glob("*_peaks.csv"))
    if not peaks_files:
        print("  No _peaks.csv files found")
        return

    rows = []
    for pf in peaks_files:
        run_name = pf.stem.replace("_peaks", "")
        df = pd.read_csv(pf)
        top = df.nlargest(3, "height")
        rows.append({
            "run_name": run_name,
            "n_peaks": len(df),
            "top_masses": "; ".join(f"{m:.1f}" for m in top["mass"]),
        })

    out_path = out_dir / "unidec_summary.csv"
    pd.DataFrame(rows).to_csv(out_path, index=False)
    print(f"  Wrote {out_path} ({len(peaks_files)} runs)")


def main() -> None:
    parser = argparse.ArgumentParser(description="Step 4: Summary CSV")
    parser.add_argument("-i", "--input", required=True,
                        help="Output dir with _peaks.csv files")
    args = parser.parse_args()
    print("=== Step 4: Summary ===")
    run_from_dir(Path(args.input).resolve())


if __name__ == "__main__":
    main()
