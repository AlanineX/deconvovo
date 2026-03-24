"""General-purpose parallel execution for pipeline steps.

Usage:
    from deconvovo.parallel import parallel_map

    def process_item(args_dict):
        # do work, return result dict
        return {"name": args_dict["name"], "status": "ok"}

    results = parallel_map(process_item, items_list, n_workers=8)

Requirements for worker functions:
  - Must be a top-level function (not a closure/lambda) — pickle requirement
  - Must accept a single dict argument and return a dict
  - Must be self-contained (import dependencies inside if needed)
"""
from __future__ import annotations

from multiprocessing import Pool


def parallel_map(fn, items: list, n_workers: int = 8) -> list:
    """Run fn(item) for each item in parallel. Returns results in input order.

    Falls back to sequential execution when n_workers=1 or len(items)<=1.
    """
    n = len(items)
    if n == 0:
        return []
    n_workers = min(n_workers, n)
    if n_workers <= 1:
        return [fn(item) for item in items]
    with Pool(n_workers) as pool:
        return list(pool.map(fn, items))
