"""Step 1: Convert Waters .raw directories to text files."""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

from deconvovo import waters_convert as wc


def run(input_dir: Path, output_dir: Path) -> Path:
    """Convert .raw dirs to _ms.txt/_im.txt. Returns directory with text files."""
    output_dir.mkdir(parents=True, exist_ok=True)
    raw_dirs = wc.collect_raw_dirs([str(input_dir)])
    if not raw_dirs:
        print("  No .raw directories found")
        return output_dir

    wine = wc.check_wine()
    cdcreader = wc.find_cdcreader()
    dlls = wc.find_support_dlls(cdcreader.parent)

    work_dir = output_dir / ".cdcreader"
    work_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(cdcreader, work_dir / "CDCReader.exe")
    for dll in dlls:
        shutil.copy2(dll, work_dir / dll.name)

    print(f"  {len(raw_dirs)} .raw dirs → {output_dir}")
    for raw_dir in raw_dirs:
        run_name = raw_dir.stem
        if (output_dir / f"{run_name}_ms.txt").exists():
            print(f"    {run_name} — cached")
            continue
        print(f"    {run_name} ...", end="", flush=True)
        try:
            wc.convert_one_raw(raw_dir, output_dir, work_dir)
            print(" OK")
        except Exception as e:
            print(f" ERROR: {e}")

    return output_dir


def main() -> None:
    parser = argparse.ArgumentParser(description="Step 1: Convert .raw → text")
    parser.add_argument("-i", "--input", required=True, help="Dir with .raw folders")
    parser.add_argument("-o", "--output", required=True, help="Output dir for text files")
    args = parser.parse_args()
    print("=== Step 1: Convert ===")
    run(Path(args.input).resolve(), Path(args.output).resolve())


if __name__ == "__main__":
    main()
