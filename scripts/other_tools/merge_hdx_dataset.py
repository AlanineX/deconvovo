#!/usr/bin/env python3
"""Merge split HDX data folders into one reusable dataset directory."""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import re
import shutil
from pathlib import Path


RAW_PATTERN = re.compile(r"(?P<time>\d+(?:min|h|hr|hrs|hour|hours|d|day|days))\.raw$", re.IGNORECASE)


def parse_time_to_min_from_name(name: str) -> tuple[str, float | None]:
    lower = name.lower()
    if "peptide_mapping" in lower:
        return "peptide_mapping", None

    m = RAW_PATTERN.search(name)
    if not m:
        return "unknown", None

    time_label = m.group("time").lower()
    m2 = re.match(r"(\d+)([a-z]+)", time_label)
    if not m2:
        return time_label, None

    value = float(m2.group(1))
    unit = m2.group(2)
    if unit == "min":
        return time_label, value
    if unit in {"h", "hr", "hrs", "hour", "hours"}:
        return time_label, value * 60.0
    if unit in {"d", "day", "days"}:
        return time_label, value * 1440.0
    return time_label, None


def sha256_of_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def unique_target_path(target_dir: Path, filename: str) -> Path:
    out = target_dir / filename
    if not out.exists():
        return out
    stem = Path(filename).stem
    suffix = Path(filename).suffix
    i = 2
    while True:
        candidate = target_dir / f"{stem}.dup{i}{suffix}"
        if not candidate.exists():
            return candidate
        i += 1


def link_or_copy(src: Path, dst: Path, mode: str) -> None:
    if mode == "copy":
        shutil.copy2(src, dst)
        return
    os.symlink(src.resolve(), dst)


def main() -> None:
    parser = argparse.ArgumentParser(description="Merge HDX split folders into one dataset directory.")
    parser.add_argument(
        "--input-root",
        default="data_1",
        help="Root folder to scan for .raw and Result.xlsx files (default: data_1)",
    )
    parser.add_argument(
        "--out-dir",
        default="data_1/merged_hdx_dataset",
        help="Output merged folder (default: data_1/merged_hdx_dataset)",
    )
    parser.add_argument(
        "--mode",
        choices=["symlink", "copy"],
        default="symlink",
        help="Materialize merged files as symlink or copy (default: symlink)",
    )
    args = parser.parse_args()

    input_root = Path(args.input_root).resolve()
    out_dir = Path(args.out_dir).resolve()
    raw_out = out_dir / "raw"
    result_out = out_dir / "result"
    raw_out.mkdir(parents=True, exist_ok=True)
    result_out.mkdir(parents=True, exist_ok=True)

    files = sorted(input_root.glob("**/*"))
    raw_files = [p for p in files if p.is_file() and p.suffix.lower() == ".raw"]
    xlsx_files = [p for p in files if p.is_file() and p.suffix.lower() == ".xlsx"]

    manifest_rows: list[dict] = []
    seen_sig_to_target: dict[tuple[int, str], Path] = {}

    for src in raw_files + xlsx_files:
        size = src.stat().st_size
        digest = sha256_of_file(src)
        sig = (size, digest)

        if src.suffix.lower() == ".raw":
            target_dir = raw_out
            time_label, time_min = parse_time_to_min_from_name(src.name)
            kind = "raw"
        else:
            target_dir = result_out
            time_label, time_min = "analysis_table", None
            kind = "xlsx"

        if sig in seen_sig_to_target:
            target = seen_sig_to_target[sig]
            status = "duplicate_skipped"
        else:
            target = unique_target_path(target_dir, src.name)
            link_or_copy(src, target, args.mode)
            seen_sig_to_target[sig] = target
            status = "linked" if args.mode == "symlink" else "copied"

        manifest_rows.append(
            {
                "source_path": str(src),
                "merged_path": str(target),
                "file_type": kind,
                "file_name": src.name,
                "time_label": time_label,
                "time_min": time_min,
                "size_bytes": size,
                "sha256": digest,
                "status": status,
            }
        )

    manifest = out_dir / "manifest.csv"
    with manifest.open("w", newline="", encoding="utf-8") as handle:
        fieldnames = [
            "source_path",
            "merged_path",
            "file_type",
            "file_name",
            "time_label",
            "time_min",
            "size_bytes",
            "sha256",
            "status",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(manifest_rows)

    print(f"Input root: {input_root}")
    print(f"Merged output: {out_dir}")
    print(f"RAW files seen: {len(raw_files)}")
    print(f"XLSX files seen: {len(xlsx_files)}")
    print(f"Manifest: {manifest}")


if __name__ == "__main__":
    main()
