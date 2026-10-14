#!/usr/bin/env python3
"""Extract an unsmoothed target-m/z mobiligram using the existing CCS extractor."""
from __future__ import annotations

import argparse
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="Converted *_im.txt file")
    parser.add_argument("--mz", required=True, type=float, help="Target mass-to-charge ratio")
    parser.add_argument("--mz-window", type=float, default=2.0, help="Extraction half-width (default: 2)")
    parser.add_argument("--raw-dir", type=Path, help="Parent folder containing the matching .raw acquisition")
    parser.add_argument("-o", "--output", required=True, type=Path, help="Output prefix for .csv and .png")
    args = parser.parse_args()
    if args.mz <= 0 or args.mz_window <= 0:
        parser.error("--mz and --mz-window must be positive")

    import numpy as np
    import pandas as pd
    from deconvovo.io import read_im_txt, find_pusher_for_run
    from deconvovo.ccs_peak_pick import extract_drift_profile
    import matplotlib.pyplot as plt

    data = read_im_txt(args.input)
    if data.empty:
        parser.error("No IM points found")
    profile = extract_drift_profile(data, args.mz, args.mz_window)
    if not np.any(profile > 0):
        parser.error("No positive intensity in the requested extraction window")
    bins = np.arange(len(profile))
    run_name = args.input.stem.removesuffix("_im")
    pusher = find_pusher_for_run(args.input.parent, run_name, args.raw_dir)
    columns = {"drift_bin": bins}
    if pusher is not None:
        axis = (bins + 0.5) * pusher / 1000.0
        columns["drift_time_ms"] = axis
        label = "Drift time (ms)"
    else:
        axis = bins
        label = "Drift bin"
    columns["intensity"] = profile
    args.output.parent.mkdir(parents=True, exist_ok=True)
    csv_path = Path(str(args.output) + ".csv")
    png_path = Path(str(args.output) + ".png")
    pd.DataFrame(columns).to_csv(csv_path, index=False)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(axis, profile, color="steelblue", linewidth=1.3)
    ax.set(xlabel=label, ylabel="Summed intensity",
           title=f"{run_name}\nm/z {args.mz:g} ± {args.mz_window:g}")
    fig.tight_layout()
    fig.savefig(png_path, dpi=150)
    plt.close(fig)
    print(csv_path)
    print(png_path)
    print(f"Pusher: {pusher} µs/bin" if pusher is not None else "Pusher unavailable: exported drift bins without time units")


if __name__ == "__main__":
    main()
