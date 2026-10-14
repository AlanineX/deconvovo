#!/usr/bin/env python3
"""Export supported low-resolution Waters functions to retention-time CSV grids."""
from __future__ import annotations

import argparse
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="Supported Waters .raw folder (not HRMS/IM-MS)")
    parser.add_argument("-o", "--out-dir", required=True, type=Path, help="Output folder")
    args = parser.parse_args()
    try:
        import rainbow
    except ImportError:
        parser.error("Install dependencies: python -m pip install -r requirements-rainbow.txt")

    try:
        data = rainbow.read(str(args.input))
    except Exception as error:
        parser.exit(1, f"Rainbow could not read this dataset: {error}\n"
                    "For Waters HRMS/IM-MS, use waters_convert.py or MSConvert.\n")
    if not data.datafiles:
        parser.error("Rainbow did not find any supported data files")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    for function in data.datafiles:
        target = args.out_dir / f"{args.input.stem}_{Path(function.name).stem}.csv"
        function.export_csv(str(target))
        print(target)


if __name__ == "__main__":
    main()
