"""Wrapper CLI: runs all pipeline steps as a DAG with parallel workers."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from deconvovo.taskqueue import TaskQueue


def _get_wc():
    from deconvovo import waters_convert
    return waters_convert


# === Top-level worker functions (picklable) ===

def _convert_one(args: dict) -> dict:
    """Worker: convert one .raw dir to text files."""
    wc = _get_wc()
    run_name = args["run_name"]
    raw_dir = Path(args["raw_dir"])
    out_dir = Path(args["out_dir"])
    work_dir = Path(args["work_dir"])
    result = {"run_name": run_name, "step": "convert", "status": []}
    if (out_dir / f"{run_name}_ms.txt").exists():
        result["status"].append("cached")
        return result
    try:
        wc.convert_one_raw(raw_dir, out_dir, work_dir)
        result["status"].append("OK")
    except Exception as e:
        result["error"] = str(e)
        result["status"].append(f"ERROR: {e}")
    return result


def _deconv_one(args: dict) -> dict:
    """Worker: run UniDec deconvolution on one run."""
    from deconvovo.imms_deconv import run_unidec_ms
    ms_file = Path(args["ms_file"])
    run_name = args["run_name"]
    out_dir = Path(args["out_dir"])
    result = {"run_name": run_name, "step": "deconv", "status": []}
    try:
        r = run_unidec_ms(ms_file, run_name, out_dir,
                          args.get("mass_range"), args.get("charge_range"),
                          args.get("mass_bins"))
        result.update(r)
        result["status"].append(f"{r['n_peaks']} peaks")
    except Exception as e:
        result["error"] = str(e)
        result["status"].append(f"ERROR: {e}")
    return result


def _plot_one(args: dict) -> dict:
    """Worker: generate interactive IM-MS plot for one run."""
    from deconvovo.imms_html import plot_im_data
    ms_file = Path(args["ms_file"])
    im_file = Path(args["im_file"])
    run_name = args["run_name"]
    out_dir = Path(args["out_dir"])
    result = {"run_name": run_name, "step": "plot", "status": []}
    if not im_file.exists():
        result["status"].append("no _im.txt")
        return result
    try:
        plot_im_data(im_file, ms_file, run_name, out_dir,
                     pusher_us=args.get("pusher_us"))
        result["status"].append("IM")
    except Exception as e:
        result["error"] = str(e)
        result["status"].append(f"ERROR: {e}")
    return result



def _summary(args: dict) -> dict:
    """Worker: generate summary CSV."""
    from deconvovo import imms_summary
    results = args.get("results", [])
    out_dir = Path(args["out_dir"])
    imms_summary.run(results, out_dir)
    return {"step": "summary", "status": ["OK"]}


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Waters IM-MS analysis pipeline with interactive 2D viewer.",
    )
    parser.add_argument("-i", "--input", required=True,
                        help="Input directory (.raw folders or _ms.txt/_im.txt text files)")
    parser.add_argument("-o", "--output", required=True, help="Output directory")
    parser.add_argument("--mass-range", type=float, nargs=2, default=None,
                        help="Mass range for deconvolution (Da)")
    parser.add_argument("--charge-range", type=int, nargs=2, default=None,
                        help="Charge state range")
    parser.add_argument("--mass-bins", type=float, default=None,
                        help="Mass bin size (Da)")
    parser.add_argument("--skip-deconv", action="store_true",
                        help="Skip UniDec deconvolution")
    parser.add_argument("--skip-plots", action="store_true",
                        help="Skip HTML plot generation")
    parser.add_argument("--skip-existing", action="store_true",
                        help="Skip runs whose outputs already exist")
    parser.add_argument("--raw-dir", default=None,
                        help="Path to .raw directories (for drift time calibration)")
    parser.add_argument("-j", "--workers", type=int, default=8,
                        help="Number of parallel workers (default: 8)")

    args = parser.parse_args()
    input_dir = Path(args.input).resolve()
    out_dir = Path(args.output).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    # --- Detect input type ---
    ms_files = sorted(input_dir.glob("*_ms.txt"))
    has_raw = any(
        d.is_dir() and d.suffix.lower() == ".raw" for d in input_dir.iterdir()
    ) if input_dir.is_dir() else False

    raw_dir_path = Path(args.raw_dir).resolve() if args.raw_dir else None
    needs_convert = not ms_files and has_raw

    if needs_convert:
        converted_dir = out_dir / "_converted"
        converted_dir.mkdir(parents=True, exist_ok=True)
        if raw_dir_path is None:
            raw_dir_path = input_dir
        data_dir = converted_dir

        # Setup CDCReader work dir (shared, read-only after setup)
        wc = _get_wc()
        work_dir = wc.setup_work_dir(converted_dir)

        raw_dirs = sorted(
            d for d in input_dir.iterdir()
            if d.is_dir() and d.suffix.lower() == ".raw"
        )
        run_names = [d.stem for d in raw_dirs]
        print(f"=== Pipeline: {len(run_names)} runs, {args.workers} workers ===")
        print(f"Input: {input_dir}")
        print(f"Output: {out_dir}")
    elif ms_files:
        data_dir = input_dir
        run_names = [f.stem.replace("_ms", "") for f in ms_files]
        print(f"=== Pipeline: {len(run_names)} runs, {args.workers} workers ===")
        print(f"Input (text): {data_dir}")
        print(f"Output: {out_dir}")
    else:
        print(f"No _ms.txt files or .raw directories found in {input_dir}")
        sys.exit(1)

    # --- Find pusher periods ---
    from deconvovo.io import find_pusher_for_run
    pusher_cache = {}
    for rn in run_names:
        pp = find_pusher_for_run(data_dir, rn, raw_dir_path)
        if pp is not None:
            pusher_cache[rn] = pp

    if pusher_cache:
        from collections import defaultdict
        by_period = defaultdict(list)
        for rn, pp in pusher_cache.items():
            by_period[f"{pp:.2f}"].append(rn)
        for period, runs in sorted(by_period.items()):
            print(f"Pusher: {period} μs/bin — {runs[0]}" +
                  (f" (+{len(runs)-1} more)" if len(runs) > 1 else ""))
    print()

    # --- Build DAG ---
    q = TaskQueue(n_workers=args.workers)
    deconv_ids = []

    for rn in run_names:
        # Step 1: Convert (if needed)
        if needs_convert:
            raw_dir = input_dir / f"{rn}.raw"
            t_conv = q.add(_convert_one, {
                "run_name": rn, "raw_dir": str(raw_dir),
                "out_dir": str(data_dir), "work_dir": str(work_dir),
            })
            deps = [t_conv]
        else:
            deps = []

        ms_path = str(data_dir / f"{rn}_ms.txt")
        im_path = str(data_dir / f"{rn}_im.txt")

        # Step 2: Deconv (depends on convert)
        if not args.skip_deconv:
            if not (args.skip_existing and (out_dir / f"{rn}_peaks.csv").exists()):
                t_dec = q.add(_deconv_one, {
                    "ms_file": ms_path, "run_name": rn,
                    "out_dir": str(out_dir),
                    "mass_range": args.mass_range, "charge_range": args.charge_range,
                    "mass_bins": args.mass_bins,
                }, depends_on=deps)
                deconv_ids.append(t_dec)

        # Step 3: Plot (depends on convert, NOT on deconv)
        if not args.skip_plots:
            if not (args.skip_existing and (out_dir / f"{rn}_2d_imms.html").exists()):
                q.add(_plot_one, {
                    "ms_file": ms_path, "im_file": im_path,
                    "run_name": rn,
                    "out_dir": str(out_dir),
                    "pusher_us": pusher_cache.get(rn),
                }, depends_on=deps)

    # Run DAG
    print(f"Queued {len(q.tasks)} tasks")
    results = q.run()

    # Report
    by_run = {}
    for tid, r in sorted(results.items()):
        rn = r.get("run_name", "")
        step = r.get("step", "?")
        status = ", ".join(r.get("status", []))
        if rn:
            by_run.setdefault(rn, []).append(f"{step}:{status}")

    for rn, steps in by_run.items():
        print(f"  {rn}  {' | '.join(steps)}")

    # Summary
    if deconv_ids:
        deconv_results = [results[tid] for tid in deconv_ids if tid in results]
        if any(r.get("n_peaks") for r in deconv_results):
            from deconvovo import imms_summary
            imms_summary.run(deconv_results, out_dir)

    n_files = sum(1 for _ in out_dir.iterdir() if _.is_file())
    print(f"\nDone: {n_files} files in {out_dir}")


if __name__ == "__main__":
    main()
