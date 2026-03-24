#!/usr/bin/env python3
"""
Waters MassLynx .raw converter using UniDec's CDCReader.exe.

On Windows: runs CDCReader.exe natively.
On Linux/macOS: runs CDCReader.exe via Wine.

Converts Waters .raw directories to analysis-ready text files:
  - MS spectrum: 2-column (m/z, intensity)
  - IM data: 3-column (m/z, drift_bin, intensity) for ion mobility

Usage:
    python -m deconvovo.waters_convert data_2_waters/20260216 -o output/converted
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

IS_WINDOWS = sys.platform == "win32"


def _get_base_path() -> Path:
    """Get base path — handles both normal Python and PyInstaller frozen exe."""
    if getattr(sys, 'frozen', False):
        return Path(sys._MEIPASS)
    return Path(__file__).parent.parent


def find_cdcreader() -> Path:
    """Locate CDCReader.exe from UniDec installation or bundled location."""
    # Check PyInstaller bundle first
    if getattr(sys, 'frozen', False):
        candidate = Path(sys._MEIPASS) / "unidec" / "bin" / "CDCReader.exe"
        if candidate.exists():
            return candidate

    # Check UniDec package in current environment
    try:
        import unidec
        pkg_dir = Path(unidec.__file__).parent
        candidate = pkg_dir / "bin" / "CDCReader.exe"
        if candidate.exists():
            return candidate
    except ImportError:
        pass

    # Search venvs in project directory (e.g. .venv-ms/)
    project_root = Path(__file__).parent.parent
    for venv_dir in sorted(project_root.glob(".venv*")):
        for candidate in venv_dir.rglob("CDCReader.exe"):
            return candidate

    raise FileNotFoundError(
        "CDCReader.exe not found. Install UniDec (pip install unidec) "
        "or ensure a .venv* directory with UniDec exists in the project root."
    )


def find_support_dlls(cdcreader_dir: Path) -> list[Path]:
    """Find MassLynxRaw.dll and cdt.dll needed by CDCReader."""
    dlls = []
    for name in ["MassLynxRaw.dll", "cdt.dll"]:
        for search_dir in [cdcreader_dir, cdcreader_dir.parent / "UniDecImporter" / "Waters"]:
            p = search_dir / name
            if p.exists():
                dlls.append(p)
                break
    return dlls


def _to_native_path(p: Path) -> str:
    """Convert path for CDCReader: native on Windows, Z: prefix on Linux/Wine."""
    if IS_WINDOWS:
        return str(p)
    return f"Z:{p}"


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
    """Convert one .raw directory using CDCReader.exe.

    Returns dict with paths to output files and metadata.
    """
    run_name = raw_dir.stem
    ms_out = out_dir / f"{run_name}_ms.txt"
    im_out = out_dir / f"{run_name}_im.txt"

    cdcreader = work_dir / "CDCReader.exe"

    # Build command: native on Windows, via Wine on Linux
    if IS_WINDOWS:
        cmd = [str(cdcreader)]
    else:
        wine = shutil.which("wine")
        if not wine:
            raise FileNotFoundError(
                "Wine not found. Install wine to run CDCReader.exe on Linux.")
        cmd = ["wine", str(cdcreader)]

    cmd += ["-r", _to_native_path(raw_dir.resolve()), "--fn", str(func)]

    if not skip_ms:
        cmd += ["-m", _to_native_path(ms_out.resolve())]
    else:
        cmd += ["--skip_ms", "1"]

    if not skip_im:
        cmd += ["-i", _to_native_path(im_out.resolve())]
    else:
        cmd += ["--skip_im", "1"]

    cmd += ["--ms_bin", str(ms_bin)]
    cmd += ["--im_bin", str(im_bin)]

    # Environment: suppress Wine debug noise on Linux
    env = dict(os.environ)
    if not IS_WINDOWS:
        env["WINEDEBUG"] = "-all"

    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        timeout=120,
        env=env,
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


def setup_work_dir(out_dir: Path) -> Path:
    """Set up work directory with CDCReader + DLLs. Returns work_dir path."""
    cdcreader_path = find_cdcreader()
    dlls = find_support_dlls(cdcreader_path.parent)

    work_dir = out_dir / ".cdcreader"
    work_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(cdcreader_path, work_dir / "CDCReader.exe")
    for dll in dlls:
        shutil.copy2(dll, work_dir / dll.name)

    return work_dir


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

    work_dir = setup_work_dir(out_dir)
    cdcreader_path = find_cdcreader()

    raw_dirs = collect_raw_dirs(args.inputs)
    if not raw_dirs:
        print("No .raw directories found.")
        sys.exit(1)

    platform = "native" if IS_WINDOWS else "Wine"
    print(f"Waters Convert: {len(raw_dirs)} .raw dirs -> {out_dir}")
    print(f"CDCReader: {cdcreader_path} ({platform})")
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

    pd.DataFrame(results).to_csv(out_dir / "convert_summary.csv", index=False)
    n_ok = sum(1 for r in results if r.get("status") == "complete")
    print(f"\nDone: {n_ok}/{len(results)} converted")
    print(f"Output: {out_dir}")


if __name__ == "__main__":
    main()
