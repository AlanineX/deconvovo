"""Step 1: Convert Waters .raw directories to text files."""
from __future__ import annotations

import argparse
from pathlib import Path

from deconvovo import waters_convert as wc
from deconvovo.parallel import parallel_map


def _convert_one(args: dict) -> dict:
    """Worker: convert one .raw dir. Top-level for pickling."""
    raw_dir = Path(args["raw_dir"])
    out_dir = Path(args["out_dir"])
    work_dir = Path(args["work_dir"])
    run_name = raw_dir.stem

    if (out_dir / f"{run_name}_ms.txt").exists():
        return {"run_name": run_name, "status": "cached"}

    try:
        info = wc.convert_one_raw(raw_dir, out_dir, work_dir)
        return {"run_name": run_name, "status": info.get("status", "unknown")}
    except Exception as e:
        return {"run_name": run_name, "status": f"ERROR: {e}"}


def run(input_dir: Path, output_dir: Path, n_workers: int = 4) -> Path:
    """Convert .raw dirs to _ms.txt/_im.txt. Returns directory with text files."""
    output_dir.mkdir(parents=True, exist_ok=True)
    raw_dirs = wc.collect_raw_dirs([str(input_dir)])
    if not raw_dirs:
        print("  No .raw directories found")
        return output_dir

    work_dir = wc.setup_work_dir(output_dir)

    print(f"  {len(raw_dirs)} .raw dirs → {output_dir}")

    items = [{"raw_dir": str(rd), "out_dir": str(output_dir), "work_dir": str(work_dir)}
             for rd in raw_dirs]

    results = parallel_map(_convert_one, items, n_workers=n_workers)

    for r in results:
        print(f"    {r['run_name']} — {r['status']}")

    return output_dir


def main() -> None:
    parser = argparse.ArgumentParser(description="Step 1: Convert .raw → text")
    parser.add_argument("-i", "--input", required=True, help="Dir with .raw folders")
    parser.add_argument("-o", "--output", required=True, help="Output dir for text files")
    parser.add_argument("-j", "--workers", type=int, default=4, help="Parallel workers")
    args = parser.parse_args()
    print("=== Step 1: Convert ===")
    run(Path(args.input).resolve(), Path(args.output).resolve(), n_workers=args.workers)


if __name__ == "__main__":
    main()
