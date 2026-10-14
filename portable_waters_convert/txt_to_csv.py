#!/usr/bin/env python3
"""Convert CDCReader MS/IM text tables to CSV without rounding the values."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def convert(source: Path) -> Path:
    target = source.with_suffix(".csv")
    with source.open(encoding="utf-8-sig") as infile:
        first = next((line.split() for line in infile if line.strip()), None)
        if first is None or len(first) not in (2, 3):
            raise ValueError(f"Expected a nonempty 2-column MS or 3-column IM table: {source}")
        header = ["mz", "intensity"] if len(first) == 2 else ["mz", "drift_bin", "intensity"]
        with target.open("w", newline="", encoding="utf-8") as outfile:
            writer = csv.writer(outfile)
            writer.writerow(header)
            writer.writerow(first)
            for line in infile:
                row = line.split()
                if not row:
                    continue
                if len(row) != len(header):
                    raise ValueError(f"Unexpected column count in {source}: {line.strip()}")
                writer.writerow(row)
    return target


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+", type=Path, help="MS/IM .txt files or output folders")
    args = parser.parse_args()
    files = []
    for path in args.inputs:
        files.extend(sorted(path.glob("*_ms.txt")) + sorted(path.glob("*_im.txt"))
                     if path.is_dir() else [path])
    if not files:
        parser.error("No *_ms.txt or *_im.txt files found")
    for source in dict.fromkeys(files):
        print(convert(source))


if __name__ == "__main__":
    main()
