#!/usr/bin/env python3
"""Dispatch vendor-specific standardization actions to the best local route."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


def detect_format(path: Path) -> str:
    if path.is_file() and path.suffix.lower() == ".raw":
        return "thermo_raw_file"
    if path.is_dir() and path.suffix.lower() == ".raw":
        return "waters_masslynx_raw_dir"
    if path.is_dir() and path.suffix.lower() == ".d":
        return "agilent_d_dir"
    if path.is_file() and path.suffix.lower() == ".mzml":
        return "mzml"
    return "unknown"


def main() -> None:
    parser = argparse.ArgumentParser(description="Dispatch standardization for mixed-vendor MS inputs")
    parser.add_argument("inputs", nargs="+", help="Input file(s) or directories")
    parser.add_argument("-o", "--out-dir", default="output/01_standardized")
    parser.add_argument("--jobs", type=int, default=8)
    args = parser.parse_args()

    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    thermo_files: list[Path] = []
    passthrough_files: list[Path] = []
    unsupported: list[tuple[Path, str]] = []

    for item in args.inputs:
        path = Path(item).resolve()
        fmt = detect_format(path)
        if fmt == "thermo_raw_file":
            thermo_files.append(path)
        elif fmt == "mzml":
            passthrough_files.append(path)
        elif fmt in {"waters_masslynx_raw_dir", "agilent_d_dir"}:
            unsupported.append((path, fmt))
        else:
            unsupported.append((path, "unknown"))

    if thermo_files:
        thermo_dir = out_dir / "thermo_mzml"
        thermo_dir.mkdir(parents=True, exist_ok=True)
        raw_parent = thermo_files[0].parent
        if any(path.parent != raw_parent for path in thermo_files):
            raise SystemExit("Thermo RAW dispatch currently requires all RAW files to share one parent directory")

        cmd = [
            sys.executable,
            "scripts/raw_to_mzml_parallel.py",
            "--input-dir",
            str(raw_parent),
            "--out-dir",
            str(thermo_dir),
            "--jobs",
            str(args.jobs),
        ]
        subprocess.run(cmd, check=True)

    report_path = out_dir / "standardization_report.md"
    lines = ["# Standardization Dispatch Report", ""]

    if thermo_files:
        lines.append("## Thermo")
        lines.append("")
        lines.append(f"- Converted RAW files from `{thermo_files[0].parent}` into `{out_dir / 'thermo_mzml'}`")
        lines.append("")

    if passthrough_files:
        lines.append("## Already Standardized")
        lines.append("")
        for path in passthrough_files:
            lines.append(f"- `{path}`")
        lines.append("")

    if unsupported:
        lines.append("## Metadata-Only / Unsupported Local Standardization")
        lines.append("")
        for path, fmt in unsupported:
            if fmt == "waters_masslynx_raw_dir":
                lines.append(f"- `{path}`: Waters MassLynx RAW detected; local converter not installed here")
            elif fmt == "agilent_d_dir":
                lines.append(f"- `{path}`: Agilent `.d` detected; local converter not installed here")
            else:
                lines.append(f"- `{path}`: unsupported or unknown input")
        lines.append("")

    report_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote: {report_path}")


if __name__ == "__main__":
    main()
