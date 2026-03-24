#!/usr/bin/env python3
"""Convert Thermo RAW files to mzML with optional parallel workers."""

from __future__ import annotations

import argparse
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


def output_name_for(raw_file: Path, fmt: int) -> str:
    if fmt == 2:
        return f"{raw_file.stem}.mzML"
    if fmt == 1:
        return f"{raw_file.stem}.mzML"
    if fmt == 0:
        return f"{raw_file.stem}.mgf"
    if fmt == 3:
        return f"{raw_file.stem}.parquet"
    return raw_file.name


def convert_one(
    parser_bin: Path,
    raw_file: Path,
    out_dir: Path,
    fmt: int,
    metadata_fmt: int,
    log_level: int,
    overwrite: bool,
) -> tuple[str, str]:
    out_name = output_name_for(raw_file, fmt)
    out_path = out_dir / out_name
    if out_path.exists() and not overwrite:
        return str(raw_file), "skipped"

    cmd = [
        str(parser_bin),
        "-i",
        str(raw_file),
        "-o",
        str(out_dir),
        "-f",
        str(fmt),
        "-m",
        str(metadata_fmt),
        "-l",
        str(log_level),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        err = proc.stderr.strip() or proc.stdout.strip() or "unknown error"
        raise RuntimeError(f"Conversion failed for {raw_file.name}: {err}")
    return str(raw_file), "ok"


def main() -> None:
    parser = argparse.ArgumentParser(description="Parallel RAW -> mzML conversion helper")
    parser.add_argument("--input-dir", required=True, help="Directory containing .raw files")
    parser.add_argument("--out-dir", required=True, help="Output directory for mzML files")
    parser.add_argument(
        "--thermo-parser",
        default="tools/ThermoRawFileParser-linux/ThermoRawFileParser",
        help="Path to ThermoRawFileParser executable",
    )
    parser.add_argument("--format", type=int, default=1, help="Output format code (default: 1 mzML)")
    parser.add_argument("--metadata-format", type=int, default=2, help="Metadata format code (default: 2 none)")
    parser.add_argument("--log-level", type=int, default=2, help="Parser log level (default: 2)")
    parser.add_argument("--jobs", type=int, default=1, help="Parallel workers (max 8)")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing outputs")
    args = parser.parse_args()

    input_dir = Path(args.input_dir).resolve()
    out_dir = Path(args.out_dir).resolve()
    parser_bin = Path(args.thermo_parser).resolve()

    if not parser_bin.exists():
        raise SystemExit(f"ThermoRawFileParser not found: {parser_bin}")
    if not input_dir.exists():
        raise SystemExit(f"Input directory not found: {input_dir}")

    out_dir.mkdir(parents=True, exist_ok=True)
    raw_files = sorted(input_dir.glob("*.raw"))
    if not raw_files:
        raise SystemExit(f"No .raw files found in {input_dir}")

    workers = min(max(1, args.jobs), 8, len(raw_files))
    print(f"RAW files: {len(raw_files)} | workers: {workers}")

    ok = 0
    skipped = 0

    if workers == 1:
        for rf in raw_files:
            _, status = convert_one(
                parser_bin=parser_bin,
                raw_file=rf,
                out_dir=out_dir,
                fmt=args.format,
                metadata_fmt=args.metadata_format,
                log_level=args.log_level,
                overwrite=args.overwrite,
            )
            print(f"{rf.name}: {status}")
            if status == "ok":
                ok += 1
            else:
                skipped += 1
    else:
        with ThreadPoolExecutor(max_workers=workers) as executor:
            futures = {
                executor.submit(
                    convert_one,
                    parser_bin,
                    rf,
                    out_dir,
                    args.format,
                    args.metadata_format,
                    args.log_level,
                    args.overwrite,
                ): rf
                for rf in raw_files
            }
            for fut in as_completed(futures):
                rf = futures[fut]
                _, status = fut.result()
                print(f"{rf.name}: {status}")
                if status == "ok":
                    ok += 1
                else:
                    skipped += 1

    print(f"Done. converted={ok} skipped={skipped} out_dir={out_dir}")


if __name__ == "__main__":
    main()
