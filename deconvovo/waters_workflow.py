"""Run the Waters desktop workflow using the existing numerical modules."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Callable

from deconvovo.paths import local_path


def run(
    source: Path,
    out_dir: Path,
    *,
    convert: bool = True,
    html: bool = True,
    calibrant_csv: Path | None = None,
    analyte_csv: Path | None = None,
    method: str = "direct",
    config_path: Path | None = None,
    skip_existing: bool = False,
    func: int = 1,
    ms_bin: float = 0,
    im_bin: float = 0,
    progress: Callable[[int, int, str], None] | None = None,
) -> dict:
    """RAW folders → validated MS/IM pairs → offline viewers → optional CCS.

    Process one acquisition at a time, keeping batch memory bounded. Stop on
    reader or plotting failure before reporting a successful workflow.
    """
    try:
        from deconvovo import waters_convert as wc
    except ImportError:
        # The portable source distribution keeps its stdlib converter at root.
        import waters_convert as wc
    from deconvovo.imms_plot import run as plot_run

    source = local_path(source).resolve()
    out_dir = local_path(out_dir).resolve()
    calibrant_csv = local_path(calibrant_csv).resolve() if calibrant_csv is not None else None
    analyte_csv = local_path(analyte_csv).resolve() if analyte_csv is not None else None
    config_path = local_path(config_path).resolve() if config_path is not None else None
    if not source.is_dir():
        raise FileNotFoundError(f"Source directory does not exist: {source}")
    if method not in {"direct", "twostep"}:
        raise ValueError("CCS method must be direct or twostep")
    if analyte_csv and not calibrant_csv:
        raise ValueError("Applying CCS requires a calibrant table")
    for label, path in (("Viewer configuration", config_path),
                        ("Calibrant table", calibrant_csv),
                        ("Analyte table", analyte_csv)):
        if path is not None and not Path(path).is_file():
            raise FileNotFoundError(f"{label} does not exist: {path}")
    raw_dirs = []
    if convert:
        # Inspect all condition folders, even when some runs are direct children.
        raw_dirs = [source] if source.suffix.lower() == ".raw" else sorted(
            p for p in source.rglob("*") if p.is_dir() and p.suffix.lower() == ".raw")
    if any(out_dir.is_relative_to(raw) for raw in raw_dirs) or any(
        parent.suffix.lower() == ".raw" and parent.is_dir()
        for parent in (out_dir, *out_dir.parents)
    ):
        raise ValueError("Choose an output directory outside the RAW acquisition")
    if raw_dirs:
        stems = [p.stem for p in raw_dirs]
        if len(stems) != len({name.casefold() for name in stems}):
            raise ValueError("RAW acquisitions must have distinct run names")
        converted_dir = out_dir / "_converted"
    else:
        converted_dir = source
    names = [p.stem for p in raw_dirs] if raw_dirs else [
        p.name.removesuffix("_ms.txt") for p in sorted(source.glob("*_ms.txt"))
    ]
    if not names:
        raise FileNotFoundError("No RAW acquisitions or converted MS/IM tables found")
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "waters_workflow.json").unlink(missing_ok=True)
    converted_dir.mkdir(parents=True, exist_ok=True)
    total = len(names) * (int(bool(raw_dirs)) + int(html)) + int(bool(calibrant_csv))
    completed = 0

    def report(message: str) -> None:
        if progress:
            progress(completed, max(total, 1), message)

    report(f"Preparing {len(names)} acquisitions")
    runs = []
    settings = {"function": func, "ms_bin": ms_bin, "im_bin": im_bin,
                "ms_smoothing": False}
    settings_path = converted_dir / "conversion_settings.json"
    previous = {}
    if settings_path.is_file():
        try:
            previous = json.loads(settings_path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            pass
    if raw_dirs:
        work_dir = wc.setup_work_dir(converted_dir) if hasattr(wc, "setup_work_dir") else wc.find_cdcreader().parent
        wc.find_support_dlls(work_dir)
        for raw in raw_dirs:
            report(f"Converting {raw.stem}")
            outputs = [converted_dir / f"{raw.stem}_{kind}.txt" for kind in ("ms", "im")]
            run_settings = {**settings, "raw_source": str(raw.resolve()),
                            "source_files": [
                                [p.name, p.stat().st_size, p.stat().st_mtime_ns]
                                for p in sorted(raw.iterdir()) if p.is_file()
                            ]}
            cached = (skip_existing and previous.get(raw.stem) == run_settings
                      and all(p.is_file() and p.stat().st_size > 0 for p in outputs))
            if cached:
                info = {"run_name": raw.stem, "status": "cached"}
            else:
                # Invalidate this run's cache before a refresh; old outputs can
                # survive a failed reader run without becoming valid new data.
                previous.pop(raw.stem, None)
                settings_path.write_text(json.dumps(previous, indent=2) + "\n", encoding="utf-8")
                info = wc.convert_one_raw(raw, converted_dir, work_dir,
                                          func=func, ms_bin=ms_bin, im_bin=im_bin)
                if info.get("status") != "complete":
                    raise RuntimeError(f"{raw.stem}: {info.get('status', 'conversion failed')} {info.get('stderr', '')}")
                previous[raw.stem] = run_settings
                settings_path.write_text(json.dumps(previous, indent=2) + "\n", encoding="utf-8")
            runs.append(info)
            completed += 1
            report(f"Converted {raw.stem}")
    for name in names:
        for kind in ("ms", "im"):
            p = converted_dir / f"{name}_{kind}.txt"
            if not p.is_file() or p.stat().st_size == 0:
                raise FileNotFoundError(f"Missing or empty converted table: {p.name}")
    viewers = []
    if html:
        for name in names:
            report(f"Generating viewer: {name}")
            result = plot_run([str(converted_dir / f"{name}_ms.txt")], out_dir,
                              raw_dir=next((p.parent for p in raw_dirs if p.stem == name), None), n_workers=1,
                              config_path=config_path, skip_existing=False)
            failures = [r for r in (result or []) if any(str(s).startswith("ERROR") for s in r.get("status", []))]
            if failures:
                raise RuntimeError(f"Viewer generation failed: {failures}")
            viewer = out_dir / f"{name}_2d_imms.html"
            if not viewer.is_file() or viewer.stat().st_size == 0:
                raise RuntimeError(f"Viewer generation produced no output: {name}")
            viewers.append(str(viewer))
            completed += 1
            report(f"Viewer ready: {name}")
    calibration = None
    ccs_dir = None
    if calibrant_csv:
        report("Fitting calibration and applying analyte CCS")
        from deconvovo.imms_ccs_calibrate import run as ccs_run
        ccs_dir = out_dir / f"ccs_{method}"
        calibration = ccs_run(ccs_dir, converted_dir, Path(calibrant_csv),
                              analyte_csv=Path(analyte_csv) if analyte_csv else None,
                              conversion_method=method)
        completed += 1
        report("CCS results ready")
    result = {"status": "complete", "converted_dir": str(converted_dir),
              "viewers": viewers, "ccs_dir": str(ccs_dir) if ccs_dir else None,
              "runs": runs or [{"run_name": name, "status": "converted input"} for name in names],
              "calibration": calibration}
    (out_dir / "waters_workflow.json").write_text(
        json.dumps({k: v for k, v in result.items() if k != "calibration"}, indent=2) + "\n",
        encoding="utf-8")
    report("Workflow complete")
    return result
