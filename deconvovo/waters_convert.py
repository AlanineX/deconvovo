#!/usr/bin/env python3
"""
Waters MassLynx .raw converter using UniDec's CDCReader.exe via Wine.

Converts Waters .raw directories to analysis-ready text files:
  - MS spectrum: 2-column (m/z, intensity)
  - IM data: 3-column (m/z, drift_bin, intensity) for ion mobility

Uses the vendor MassLynxRaw.dll + cdt.dll for correct calibration.

Usage:
    python -m deconvovo.waters_convert data_2_waters/20260216 -o output/22_waters_converted
    python -m deconvovo.waters_convert mydata.raw -o out --ms-bin 0 --im-bin 0
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd


def find_cdcreader() -> Path:
    """Locate CDCReader.exe from UniDec installation."""
    try:
        import unidec
        pkg_dir = Path(unidec.__file__).parent
        candidate = pkg_dir / "bin" / "CDCReader.exe"
        if candidate.exists():
            return candidate
    except ImportError:
        pass
    raise FileNotFoundError(
        "CDCReader.exe not found. Install UniDec: pip install unidec"
    )


def find_support_dlls(cdcreader_dir: Path) -> list[Path]:
    """Find MassLynxRaw.dll and cdt.dll needed by CDCReader."""
    dlls = []
    for name in ["MassLynxRaw.dll", "cdt.dll"]:
        # Check CDCReader directory first, then UniDec Waters importer
        for search_dir in [cdcreader_dir, cdcreader_dir.parent / "UniDecImporter" / "Waters"]:
            p = search_dir / name
            if p.exists():
                dlls.append(p)
                break
    return dlls


def check_wine() -> str:
    """Check that Wine is available. Returns wine binary path."""
    wine = shutil.which("wine")
    if not wine:
        raise FileNotFoundError(
            "Wine not found. Install wine to run CDCReader.exe on Linux."
        )
    return wine


def convert_one_raw(
    raw_dir: Path,
    out_dir: Path,
    work_dir: Path,
    func: int = 1,
    ms_bin: float = 0,
    im_bin: float = 0,
    skip_ms: bool = False,
    skip_im: bool = False,
) -> dict:
    """Convert one .raw directory using CDCReader.exe via Wine.

    Returns dict with paths to output files and metadata.
    """
    run_name = raw_dir.stem
    ms_out = out_dir / f"{run_name}_ms.txt"
    im_out = out_dir / f"{run_name}_im.txt"

    # Build Wine path for the raw directory
    raw_wine_path = f"Z:{raw_dir}"

    # Build CDCReader command
    cdcreader = work_dir / "CDCReader.exe"
    cmd = [
        "wine", str(cdcreader),
        "-r", raw_wine_path,
        "--fn", str(func),
    ]

    if not skip_ms:
        cmd += ["-m", f"Z:{ms_out}"]
    else:
        cmd += ["--skip_ms", "1"]

    if not skip_im:
        cmd += ["-i", f"Z:{im_out}"]
    else:
        cmd += ["--skip_im", "1"]

    cmd += ["--ms_bin", str(ms_bin)]
    cmd += ["--im_bin", str(im_bin)]

    # Run CDCReader
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        timeout=120,
        env={**os.environ, "WINEDEBUG": "-all"},
    )

    info: dict = {
        "run_name": run_name,
        "raw_dir": str(raw_dir),
        "func": func,
        "status": "complete" if result.returncode == 0 else f"error:{result.returncode}",
    }

    if result.returncode != 0:
        info["stderr"] = result.stderr[:500]
        return info

    # Parse output files
    if not skip_ms and ms_out.exists() and ms_out.stat().st_size > 0:
        info["ms_file"] = str(ms_out)
        info["ms_size"] = ms_out.stat().st_size
    if not skip_im and im_out.exists() and im_out.stat().st_size > 0:
        info["im_file"] = str(im_out)
        info["im_size"] = im_out.stat().st_size

    return info


def collect_raw_dirs(inputs: list[str]) -> list[Path]:
    """Accept .raw directories or parent directories containing them."""
    results: list[Path] = []
    for item in inputs:
        p = Path(item).resolve()
        if p.is_dir() and p.suffix.lower() == ".raw":
            results.append(p)
        elif p.is_dir():
            for child in sorted(p.iterdir()):
                if child.is_dir() and child.suffix.lower() == ".raw":
                    results.append(child)
    seen: set[Path] = set()
    return [r for r in results if not (r in seen or seen.add(r))]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Convert Waters MassLynx .raw to text (MS + IM) via CDCReader.",
    )
    parser.add_argument("inputs", nargs="+", help=".raw dirs or parent dirs")
    parser.add_argument("-o", "--out-dir", required=True, help="Output directory")
    parser.add_argument("--fn", type=int, default=1, help="Function number (default: 1)")
    parser.add_argument("--ms-bin", type=float, default=0,
                        help="m/z bin size for MS (0=raw, default: 0)")
    parser.add_argument("--im-bin", type=float, default=0,
                        help="m/z bin size for IM (0=raw, default: 0)")
    parser.add_argument("--skip-ms", action="store_true", help="Skip MS output")
    parser.add_argument("--skip-im", action="store_true", help="Skip IM output")
    parser.add_argument("--skip-existing", action="store_true",
                        help="Skip runs whose outputs already exist")

    args = parser.parse_args()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    # Setup
    wine = check_wine()
    cdcreader_path = find_cdcreader()
    dlls = find_support_dlls(cdcreader_path.parent)

    # Create work directory with CDCReader + DLLs
    work_dir = out_dir / ".cdcreader"
    work_dir.mkdir(exist_ok=True)
    shutil.copy2(cdcreader_path, work_dir / "CDCReader.exe")
    for dll in dlls:
        shutil.copy2(dll, work_dir / dll.name)

    raw_dirs = collect_raw_dirs(args.inputs)
    if not raw_dirs:
        print("No .raw directories found.")
        sys.exit(1)

    print(f"Waters Convert: {len(raw_dirs)} .raw dirs -> {out_dir}")
    print(f"CDCReader: {cdcreader_path}")
    print(f"Wine: {wine}")
    print()

    results = []
    for raw_dir in raw_dirs:
        run_name = raw_dir.stem
        ms_exists = (out_dir / f"{run_name}_ms.txt").exists()
        im_exists = (out_dir / f"{run_name}_im.txt").exists()

        if args.skip_existing and ms_exists and im_exists:
            print(f"  {run_name} — skipped (exists)")
            results.append({"run_name": run_name, "status": "skipped"})
            continue

        print(f"  {run_name}", end="", flush=True)
        try:
            info = convert_one_raw(
                raw_dir, out_dir, work_dir,
                func=args.fn,
                ms_bin=args.ms_bin,
                im_bin=args.im_bin,
                skip_ms=args.skip_ms,
                skip_im=args.skip_im,
            )
            ms_size = info.get("ms_size", 0)
            im_size = info.get("im_size", 0)
            print(f" — {info['status']}"
                  f" (MS: {ms_size/1024:.0f}K, IM: {im_size/1024:.0f}K)")
            results.append(info)
        except Exception as e:
            print(f" — ERROR: {e}")
            results.append({"run_name": run_name, "status": f"error: {e}"})

    # Summary
    pd.DataFrame(results).to_csv(out_dir / "convert_summary.csv", index=False)
    n_ok = sum(1 for r in results if r.get("status") == "complete")
    print(f"\nDone: {n_ok}/{len(results)} converted")
    print(f"Output: {out_dir}")


if __name__ == "__main__":
    main()
