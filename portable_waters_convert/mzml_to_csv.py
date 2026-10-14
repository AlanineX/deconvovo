#!/usr/bin/env python3
"""Export per-spectrum m/z and intensity from mzML, with RT and its units."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="Input .mzML file")
    parser.add_argument("-o", "--output", required=True, type=Path, help="Output CSV")
    args = parser.parse_args()
    try:
        from pyteomics import mzml
    except ImportError:
        parser.error("Install dependencies: python -m pip install -r requirements-mzml.txt")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with mzml.MzML(str(args.input)) as reader, args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["spectrum_id", "ms_level", "retention_time", "retention_time_unit", "mz", "intensity"])
        count = 0
        for spectrum in reader:
            scan = (spectrum.get("scanList", {}).get("scan") or [{}])[0]
            rt = scan.get("scan start time", "")
            mzs = spectrum["m/z array"]
            intensities = spectrum["intensity array"]
            if len(mzs) != len(intensities):
                raise ValueError(f"Array lengths differ in {spectrum['id']}")
            writer.writerows((spectrum["id"], spectrum.get("ms level", ""), rt,
                              getattr(rt, "unit_info", ""), mass, intensity)
                             for mass, intensity in zip(mzs, intensities))
            count += len(mzs)
    print(f"Wrote {count} peak rows to {args.output}")


if __name__ == "__main__":
    main()
