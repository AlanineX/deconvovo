"""Step 3: Generate interactive 2D IM-MS HTML plots."""
from __future__ import annotations

import argparse
from pathlib import Path

from deconvovo.io import find_pusher_for_run
from deconvovo.imms_html import plot_im_data
from deconvovo.parallel import parallel_map


def _resolve_inputs(inputs: list[str]) -> tuple[Path, list[Path]]:
    """Accept a single directory OR a list of file paths (typically a shell
    glob like `_converted/*_24H_*_ms.txt`) and return (data_dir, ms_files).

    All MS files must live in the same parent directory.
    """
    paths = [Path(p).resolve() for p in inputs]
    # Single directory: enumerate _ms.txt inside
    if len(paths) == 1 and paths[0].is_dir():
        return paths[0], sorted(paths[0].glob("*_ms.txt"))
    # List of files (mixture of _ms.txt and _im.txt is OK; we keep only _ms.txt)
    ms_files = []
    for p in paths:
        if p.is_dir():
            ms_files.extend(p.glob("*_ms.txt"))
        elif p.name.endswith("_ms.txt"):
            ms_files.append(p)
        elif p.name.endswith("_im.txt"):
            ms = p.parent / p.name.replace("_im.txt", "_ms.txt")
            if ms.exists():
                ms_files.append(ms)
    ms_files = sorted(set(ms_files))
    if not ms_files:
        return paths[0], []
    parents = {f.parent for f in ms_files}
    if len(parents) > 1:
        raise SystemExit(
            f"all input files must share one parent dir; got: {sorted(parents)}"
        )
    return parents.pop(), ms_files


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
                     pusher_us=args.get("pusher_us"),
                     config_path=args.get("config_path"))
        result["status"].append("IM")
    except Exception as e:
        result["status"].append(f"ERROR: {e}")
    return result


def run(inputs, out_dir: Path, skip_existing: bool = False,
        raw_dir: Path | None = None, n_workers: int = 8,
        config_path: Path | None = None) -> None:
    """`inputs` is either a Path (legacy: a single data_dir) or a list of
    file/dir paths (typically from a shell glob)."""
    out_dir.mkdir(parents=True, exist_ok=True)
    if isinstance(inputs, (str, Path)):
        inputs = [str(inputs)]
    data_dir, ms_files = _resolve_inputs(list(inputs))
    print(f"  {len(ms_files)} runs from {data_dir}")

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
        html_out = out_dir / f"{run_name}_2d_imms.html"
        if skip_existing and html_out.exists() and html_out.stat().st_size > 0:
            continue
        if not (data_dir / f"{run_name}_im.txt").exists():
            continue
        run_args.append({
            "ms_file": str(ms_file), "run_name": run_name,
            "im_file": str(data_dir / f"{run_name}_im.txt"),
            "out_dir": str(out_dir),
            "pusher_us": pusher_cache.get(run_name),
            "config_path": str(config_path) if config_path else None,
        })

    print(f"  Processing {len(run_args)} runs")
    results = parallel_map(_plot_one_run, run_args, n_workers=n_workers)
    for r in results:
        print(f"    {r['run_name']} {' — '.join(r.get('status', []))}")
    return results


def main() -> None:
    parser = argparse.ArgumentParser(description="Step 3: Interactive IM-MS plots")
    parser.add_argument("-i", "--input", nargs="+", required=True,
                        help="A directory of *_ms.txt files, OR a list of files "
                             "(shell globs work, e.g. converted/*_24H_*_ms.txt)")
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--raw-dir", default=None)
    parser.add_argument("--config", default=None,
                        help="Path to custom imms_plot_config.json "
                             "(defaults / presets / initial mz+drift ranges)")
    parser.add_argument("--skip-existing", action="store_true")
    parser.add_argument("-j", "--workers", type=int, default=8)
    args = parser.parse_args()
    run(args.input, Path(args.output).resolve(), args.skip_existing,
        raw_dir=Path(args.raw_dir).resolve() if args.raw_dir else None,
        n_workers=args.workers,
        config_path=Path(args.config).resolve() if args.config else None)


if __name__ == "__main__":
    main()
