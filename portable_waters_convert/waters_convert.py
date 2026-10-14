#!/usr/bin/env python3
"""
Standalone Waters MassLynx .raw converter using the bundled CDCReader.exe.

On Windows: runs CDCReader.exe natively.
On Linux/macOS: runs CDCReader.exe via Wine.

Converts Waters .raw directories to analysis-ready text files:
  - MS spectrum: 2-column (m/z, intensity)
  - IM data: 3-column (m/z, drift_bin, intensity) for ion mobility

Usage:
    python waters_convert.py examples/input -o converted --fn 1
"""
from __future__ import annotations

import argparse
import csv
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from deconvovo.paths import local_path

IS_WINDOWS = sys.platform == "win32"


def find_cdcreader() -> Path:
    """Find the reader relative to this script, independent of the checkout."""
    candidate = Path(__file__).resolve().parent / "bin" / "CDCReader.exe"
    if not candidate.is_file():
        raise FileNotFoundError(f"Missing bundled reader: {candidate}")
    return candidate


def find_support_dlls(cdcreader_dir: Path) -> list[Path]:
    """Find MassLynxRaw.dll and cdt.dll needed by CDCReader."""
    dlls = [cdcreader_dir / name for name in ("MassLynxRaw.dll", "cdt.dll")]
    for dll in dlls:
        if not dll.is_file():
            raise FileNotFoundError(f"Missing bundled DLL: {dll}")
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
    if skip_ms and skip_im:
        raise ValueError("Choose at least one output; cannot skip both MS and IM")
    if func < 1 or ms_bin < 0 or im_bin < 0:
        raise ValueError("Function must be at least 1; bin widths must be nonnegative")
    raw_dir, out_dir, work_dir = (local_path(path).resolve() for path in (raw_dir, out_dir, work_dir))
    if not raw_dir.is_dir():
        raise FileNotFoundError(f"Waters .raw folder not found: {raw_dir}")
    if out_dir.is_relative_to(raw_dir) or any(
        parent.suffix.lower() == ".raw" and parent.is_dir()
        for parent in (out_dir, *out_dir.parents)
    ):
        raise ValueError("Choose an output directory outside the RAW acquisition")
    out_dir.mkdir(parents=True, exist_ok=True)
    run_name = raw_dir.stem

    cdcreader = work_dir / "CDCReader.exe"

    # Build command: native on Windows, via Wine on Linux
    if IS_WINDOWS:
        cmd = [str(cdcreader)]
    else:
        wine = shutil.which("wine")
        if not wine:
            raise FileNotFoundError(
                "Wine not found. Install wine to run CDCReader.exe on Linux.")
        cmd = [wine, str(cdcreader)]

    cmd += ["-r", _to_native_path(raw_dir.resolve()), "--fn", str(func)]

    cmd += ["--ms_bin", str(ms_bin)]
    cmd += ["--im_bin", str(im_bin)]
    # Keep the exported spectrum unaltered; smoothing belongs in the viewer.
    cmd += ["--ms_smooth_window", "0", "--ms_number_smooth", "0"]

    # Environment: suppress Wine debug noise on Linux
    env = dict(os.environ)
    if not IS_WINDOWS:
        env["WINEDEBUG"] = "-all"

    # A reader can exit successfully without writing data (e.g. a missing
    # function). Stage each run so old outputs cannot conceal that failure.
    with tempfile.TemporaryDirectory(prefix=".waters-", dir=out_dir) as staging:
        requested = {name: Path(staging) / f"{run_name}_{name}.txt"
                     for name, skipped in (("ms", skip_ms), ("im", skip_im))
                     if not skipped}
        for name, skipped in (("ms", skip_ms), ("im", skip_im)):
            cmd += ([f"--skip_{name}", "1"] if skipped else
                    [f"-{name[0]}", _to_native_path(requested[name].resolve())])
        result = subprocess.run(cmd, capture_output=True, text=True, errors="replace",
                                timeout=120, env=env)
        info = {"run_name": run_name, "raw_dir": str(raw_dir), "func": func,
                "status": "complete" if result.returncode == 0
                else f"error:{result.returncode}"}
        missing = [name for name, path in requested.items()
                   if not path.is_file() or path.stat().st_size == 0]
        if result.returncode != 0 or missing:
            if result.returncode == 0:
                info["status"] = "error: missing or empty " + ", ".join(missing) + " output"
            info["stderr"] = (result.stderr or result.stdout).strip()[:500]
            return info
        for name, path in requested.items():
            output = out_dir / path.name
            path.replace(output)
            info[f"{name}_file"] = str(output)
            info[f"{name}_size"] = output.stat().st_size
        return info


def collect_raw_dirs(inputs: list[str]) -> list[Path]:
    """Accept .raw directories or parent directories containing them."""
    results: list[Path] = []
    for item in inputs:
        p = local_path(item).resolve()
        if p.is_dir() and p.suffix.lower() == ".raw":
            results.append(p)
        elif p.is_dir():
            for child in sorted(p.iterdir()):
                if child.is_dir() and child.suffix.lower() == ".raw":
                    results.append(child)
    seen: set[Path] = set()
    return [r for r in results if not (r in seen or seen.add(r))]


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(errors="backslashreplace")
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
    if args.fn < 1 or args.ms_bin < 0 or args.im_bin < 0:
        parser.error("--fn must be at least 1; bin widths must be nonnegative")
    if args.skip_ms and args.skip_im:
        parser.error("Choose at least one output; cannot skip both MS and IM")

    raw_dirs = collect_raw_dirs(args.inputs)
    if not raw_dirs:
        parser.error("No .raw directories found; provide a complete .raw folder or its parent")

    cdcreader_path = find_cdcreader()
    find_support_dlls(cdcreader_path.parent)
    work_dir = cdcreader_path.parent
    out_dir = local_path(args.out_dir).resolve()
    if any(out_dir.is_relative_to(raw) for raw in raw_dirs) or any(
        parent.suffix.lower() == ".raw" and parent.is_dir()
        for parent in (out_dir, *out_dir.parents)
    ):
        parser.error("Choose an output directory outside the RAW acquisition")
    out_dir.mkdir(parents=True, exist_ok=True)

    platform = "native" if IS_WINDOWS else "Wine"
    print(f"Waters Convert: {len(raw_dirs)} .raw dirs -> {out_dir}")
    print(f"CDCReader: {cdcreader_path} ({platform})")
    print()

    results = []
    for raw_dir in raw_dirs:
        run_name = raw_dir.stem
        ms_out = out_dir / f"{run_name}_ms.txt"
        im_out = out_dir / f"{run_name}_im.txt"
        ms_exists = ms_out.is_file() and ms_out.stat().st_size > 0
        im_exists = im_out.is_file() and im_out.stat().st_size > 0

        if args.skip_existing and (args.skip_ms or ms_exists) and (args.skip_im or im_exists):
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
            if info.get("stderr"):
                print(f"    {info['stderr']}")
            results.append(info)
        except Exception as e:
            print(f" — ERROR: {e}")
            results.append({"run_name": run_name, "status": f"error: {e}"})

    fields = list(dict.fromkeys(key for result in results for key in result))
    with (out_dir / "convert_summary.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(results)
    n_ok = sum(1 for r in results if r.get("status") == "complete")
    print(f"\nDone: {n_ok}/{len(results)} converted")
    print(f"Output: {out_dir}")
    if any(str(result["status"]).startswith("error") for result in results):
        sys.exit(1)


if __name__ == "__main__":
    main()
