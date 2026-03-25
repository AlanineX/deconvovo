"""General-purpose DAG task queue with dependency tracking.

Workers pick up whatever task is ready next — no idle time between steps.
Tasks form a directed acyclic graph: each task can depend on others.

Usage:
    q = TaskQueue(n_workers=8)
    for run in runs:
        t1 = q.add(convert_fn, {"run": run})
        t2 = q.add(deconv_fn,  {"run": run}, depends_on=[t1])
        t3 = q.add(plot_fn,    {"run": run}, depends_on=[t1])  # plot needs convert, NOT deconv
        deconv_ids.append(t2)
    q.add(summary_fn, {"runs": runs}, depends_on=deconv_ids)
    results = q.run()

Requirements:
  - Worker functions must be top-level (picklable), accept a single dict, return a dict
  - Dependencies are by task ID (int returned by q.add)
"""
from __future__ import annotations

import os
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED


class TaskQueue:
    def __init__(self, n_workers: int = 8):
        # Force sequential in GUI to prevent spawning new windows
        if os.environ.get("DECONVOVO_GUI") == "1":
            n_workers = 1
        self.n_workers = n_workers
        self.tasks: list[tuple] = []  # (id, fn, args, deps)

    def add(self, fn, args: dict, depends_on: list[int] | None = None) -> int:
        """Add a task. Returns task ID for use in depends_on."""
        tid = len(self.tasks)
        self.tasks.append((tid, fn, args, depends_on or []))
        return tid

    def run(self) -> dict[int, dict]:
        """Execute all tasks respecting dependencies. Returns {task_id: result}."""
        n_tasks = len(self.tasks)
        if n_tasks == 0:
            return {}

        results = {}
        completed = set()
        submitted = set()
        n_workers = min(self.n_workers, n_tasks)

        # Sequential fallback
        if n_workers <= 1:
            for tid, fn, args, deps in self.tasks:
                try:
                    results[tid] = fn(args)
                except Exception as e:
                    results[tid] = {"error": str(e)}
                completed.add(tid)
            return results

        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            futures = {}  # future -> task_id

            while len(completed) < n_tasks:
                # Submit tasks whose dependencies are all completed
                for tid, fn, args, deps in self.tasks:
                    if tid not in submitted and all(d in completed for d in deps):
                        future = executor.submit(fn, args)
                        futures[future] = tid
                        submitted.add(tid)

                # Wait for at least one to finish
                if futures:
                    done, _ = wait(futures.keys(), return_when=FIRST_COMPLETED)
                    for f in done:
                        tid = futures.pop(f)
                        try:
                            results[tid] = f.result()
                        except Exception as e:
                            results[tid] = {"error": str(e)}
                        completed.add(tid)

        return results
