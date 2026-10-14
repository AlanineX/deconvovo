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

    outputs = [out_dir / f"{run_name}_{kind}.txt" for kind in ("ms", "im")]
    if all(path.is_file() and path.stat().st_size > 0 for path in outputs):
        return {"run_name": run_name, "status": "cached"}

    try:
        info = wc.convert_one_raw(raw_dir, out_dir, work_dir)
        return info
    except Exception as e:
        return {"run_name": run_name, "status": f"ERROR: {e}"}


def run(input_dir: Path, output_dir: Path, n_workers: int = 4) -> Path:
    """Convert .raw dirs to _ms.txt/_im.txt. Returns directory with text files."""
    output_dir.mkdir(parents=True, exist_ok=True)
    raw_dirs = wc.collect_raw_dirs([str(input_dir)])
    if not raw_dirs:
        raise FileNotFoundError(f"No .raw directories found in {input_dir}")

    work_dir = wc.setup_work_dir(output_dir)

    print(f"  {len(raw_dirs)} .raw dirs → {output_dir}")

    items = [{"raw_dir": str(rd), "out_dir": str(output_dir), "work_dir": str(work_dir)}
             for rd in raw_dirs]

    results = parallel_map(_convert_one, items, n_workers=n_workers)

    for r in results:
        print(f"    {r['run_name']} — {r['status']}")

    failed = [r for r in results if r.get("status") not in ("complete", "cached")]
    if failed:
        detail = "; ".join(f"{r['run_name']}: {r['status']}" for r in failed)
        raise RuntimeError(f"Waters conversion failed: {detail}")

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
