"""Step 3: Generate interactive 2D IM-MS HTML plots."""
from __future__ import annotations

import argparse
from pathlib import Path

from deconvovo.io import find_pusher_for_run
from deconvovo.imms_html import plot_im_data
from deconvovo.parallel import parallel_map


def _plot_one_run(args: dict) -> dict:
    """Worker: plot a single run. Top-level for pickling."""
    ms_file = Path(args["ms_file"])
    im_file = Path(args["im_file"])
    run_name = args["run_name"]
    out_dir = Path(args["out_dir"])

    result = {"run_name": run_name, "status": []}
    if not im_file.exists():
        result["status"].append("no _im.txt")
        return result
    try:
        plot_im_data(im_file, ms_file, run_name, out_dir,
                     pusher_us=args.get("pusher_us"))
        result["status"].append("IM")
    except Exception as e:
        result["status"].append(f"ERROR: {e}")
    return result


def run(data_dir: Path, out_dir: Path, skip_existing: bool = False,
        raw_dir: Path | None = None, n_workers: int = 8) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    ms_files = sorted(data_dir.glob("*_ms.txt"))
    print(f"  {len(ms_files)} runs")

    pusher_cache = {}
    for ms_file in ms_files:
        rn = ms_file.stem.replace("_ms", "")
        pp = find_pusher_for_run(data_dir, rn, raw_dir)
        if pp is not None:
            pusher_cache[rn] = pp

    if pusher_cache:
        from collections import defaultdict
        by_period = defaultdict(list)
        for rn, pp in pusher_cache.items():
            by_period[f"{pp:.2f}"].append(rn)
        for period, runs in sorted(by_period.items()):
            print(f"  Pusher: {period} μs — {runs[0]}" +
                  (f" (+{len(runs)-1} more)" if len(runs) > 1 else ""))

    run_args = []
    for ms_file in ms_files:
        run_name = ms_file.stem.replace("_ms", "")
        if skip_existing and (out_dir / f"{run_name}_2d_imms.html").exists():
            continue
        if not (data_dir / f"{run_name}_im.txt").exists():
            continue
        run_args.append({
            "ms_file": str(ms_file), "run_name": run_name,
            "im_file": str(data_dir / f"{run_name}_im.txt"),
            "out_dir": str(out_dir),
            "pusher_us": pusher_cache.get(run_name),
        })

    print(f"  Processing {len(run_args)} runs")
    results = parallel_map(_plot_one_run, run_args, n_workers=n_workers)
    for r in results:
        print(f"    {r['run_name']} {' — '.join(r.get('status', []))}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Step 3: Interactive IM-MS plots")
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--raw-dir", default=None)
    parser.add_argument("--skip-existing", action="store_true")
    parser.add_argument("-j", "--workers", type=int, default=8)
    args = parser.parse_args()
    run(Path(args.input).resolve(), Path(args.output).resolve(), args.skip_existing,
        raw_dir=Path(args.raw_dir).resolve() if args.raw_dir else None,
        n_workers=args.workers)


if __name__ == "__main__":
    main()
