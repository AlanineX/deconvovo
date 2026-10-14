"""Exercise the frozen reader, actual console, and single-file offline export."""
from __future__ import annotations

from copy import deepcopy
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time
import traceback


def _qt_diagnostics(application, window, panel, viewer) -> dict:
    """Record the actual runtime layout without changing fonts or widget sizes."""
    from PySide6.QtCore import QPoint, Qt, QLibraryInfo
    from PySide6.QtGui import QFontDatabase, QFontInfo, QFontMetricsF
    from PySide6.QtWidgets import QWidget

    def size(value):
        return {"width": value.width(), "height": value.height()}

    def rect(value):
        return {"x": value.x(), "y": value.y(), **size(value)}

    def font(value):
        actual = QFontInfo(value)
        metrics = QFontMetricsF(value)
        return {"requested": value.toString(), "actual_family": actual.family(),
                "point_size": actual.pointSizeF(), "pixel_size": actual.pixelSize(),
                "height": metrics.height(), "line_spacing": metrics.lineSpacing(),
                "average_char_width": metrics.averageCharWidth()}

    def widget(value):
        if value is None:
            return None
        position = value.mapToGlobal(QPoint())
        result = {"class": type(value).__name__, "object_name": value.objectName(),
                  "geometry": rect(value.geometry()), "frame_geometry": rect(value.frameGeometry()),
                  "size": size(value.size()),
                  "global_position": {"x": position.x(), "y": position.y()},
                  "size_hint": size(value.sizeHint()),
                  "minimum_size_hint": size(value.minimumSizeHint()),
                  "minimum_size": size(value.minimumSize()),
                  "maximum_size": size(value.maximumSize()),
                  "visible": value.isVisible(), "font": font(value.font()),
                  "device_pixel_ratio": value.devicePixelRatioF()}
        layout = value.layout()
        if layout is not None:
            margins = layout.contentsMargins()
            result["layout"] = {"geometry": rect(layout.geometry()), "spacing": layout.spacing(),
                                "margins": [margins.left(), margins.top(), margins.right(), margins.bottom()]}
            result["layout"]["items"] = []
            for index in range(layout.count()):
                item = layout.itemAt(index)
                child = item.widget()
                result["layout"]["items"].append({
                    "class": type(child).__name__ if child is not None else "layout_or_spacer",
                    "object_name": child.objectName() if child is not None else "",
                    "geometry": rect(item.geometry()), "size_hint": size(item.sizeHint()),
                    "minimum_size": size(item.minimumSize()),
                })
        if hasattr(value, "title"):
            result["title"] = value.title()
        if hasattr(value, "text"):
            result["text"] = value.text()
        return result

    targets = {"window": window, "central_widget": window.centralWidget(),
               "window_body": getattr(window, "body", None),
               "window_stack": getattr(window, "stack", None), "workspace": panel}
    if panel is not None:
        for name in ("console_splitter", "browser_panel", "sources_group", "plot_panel", "tabs",
                     "job_status", "job_progress", "output_summary", "display_group", "ccs_card"):
            targets[name] = getattr(panel, name, None)
        targets["plot_card"] = panel.tabs.parentWidget()
        targets["output_card"] = panel.output_summary.parentWidget()
        targets["plot_direct_children"] = [widget(child) for child in
            panel.plot_panel.findChildren(QWidget, options=Qt.FindDirectChildrenOnly)]
    if viewer is not None:
        targets.update(viewer=viewer, viewer_tabs=viewer.tabs, im_stack=viewer.im_stack,
                       native_canvas=viewer.canvas, viewer_status=viewer.status)
    direct_children = targets.pop("plot_direct_children", [])
    screen = window.screen()
    font_directory = Path(os.environ.get("QT_QPA_FONTDIR") or
                          (QLibraryInfo.path(QLibraryInfo.LibraryPath.LibrariesPath) + "/fonts"))
    families = QFontDatabase.families()
    result = {"platform": application.platformName(), "application_font": font(application.font()),
              "font_database": {"family_count": len(families), "families": families,
                                "offscreen_font_directory": str(font_directory),
                                "font_directory_exists": font_directory.is_dir()},
              "environment": {name: os.environ.get(name) for name in
                              ("QT_QPA_PLATFORM", "QT_QUICK_BACKEND", "QT_OPENGL", "QSG_RHI_BACKEND",
                               "QSG_RHI_PREFER_SOFTWARE_RENDERER", "QT_SCALE_FACTOR",
                               "QT_FONT_DPI", "QT_QPA_FONTDIR", "QTWEBENGINE_CHROMIUM_FLAGS")},
              "widgets": {name: widget(value) for name, value in targets.items()},
              "plot_direct_children": direct_children}
    if screen is not None:
        result["screen"] = {"name": screen.name(), "geometry": rect(screen.geometry()),
                            "available_geometry": rect(screen.availableGeometry()),
                            "logical_dpi": screen.logicalDotsPerInch(),
                            "physical_dpi": screen.physicalDotsPerInch(),
                            "device_pixel_ratio": screen.devicePixelRatio()}
    return result


def run(output: Path, root: Path) -> int:
    output.mkdir(parents=True, exist_ok=True)
    report = {"status": "running"}
    (output / "native_smoke.json").write_text(json.dumps(report) + "\n", encoding="utf-8")
    application = None
    window = None
    panel = None
    viewer = None
    try:
        # Explicit build-test opt-in; normal desktop startup keeps Qt's default
        # platform selection and Chromium sandbox behavior.
        if os.name == "nt":
            # Windows' offscreen plugin has no system font database. Validate
            # the same native platform and fonts used by the desktop app.
            os.environ["QT_QPA_PLATFORM"] = "windows"
        from gui.widgets.worker import set_gui_mode
        set_gui_mode()
        from PySide6.QtWidgets import QApplication
        from PySide6.QtGui import QFont, QFontDatabase, QFontInfo, QImage
        from gui.app import MainWindow
        from gui.theme import FONT_FAMILY, stylesheet
        application = QApplication([])
        application.setOrganizationName("DeconVoVo")
        application.setApplicationName("DeconVoVo Waters")
        application.setStyle("Fusion")
        application.setStyleSheet(stylesheet(dark=True))
        application.setFont(QFont(FONT_FAMILY, 10))
        from deconvovo.waters_workflow import run as process

        production = root / "examples" / "production"
        manifest = json.loads((production / "demo_manifest.json").read_text(encoding="utf-8"))
        profile = next(value for value in manifest["profiles"] if value["id"] == "protein_adp")
        adp_raw = production / profile["source"] / "20260212_200UM_ADP_200MM_AMAC.raw"
        config_path = production / profile["viewer_config"]
        config = json.loads(config_path.read_text(encoding="utf-8"))
        calibrant_path = production / profile["calibrant_csv"]
        analyte_path = production / profile["analyte_csv"]
        results = process(production / profile["source"], output / "real_demo", html=False, config_path=config_path,
                          calibrant_csv=calibrant_path, analyte_csv=analyte_path)
        if len(results["runs"]) != 5 or results["viewers"]:
            raise RuntimeError("Primary smoke requires three protein standards and two ADP acquisition outputs")
        ccs_dir = Path(results["ccs_dir"])
        with calibrant_path.open(encoding="utf-8", newline="") as handle:
            expected_references = list(csv.DictReader(handle))
        with (ccs_dir / "calibrant_summary.csv").open(encoding="utf-8", newline="") as handle:
            references = list(csv.DictReader(handle))
        with (ccs_dir / "analyte_summary.csv").open(encoding="utf-8", newline="") as handle:
            analytes = list(csv.DictReader(handle))
        if (len(expected_references) != 25 or len(references) != 25
                or results["calibration"]["n_points"] != 25 or len(analytes) != 34):
            raise RuntimeError("ADP CCS smoke requires all 25 valid protein references and 34 ADP targets")
        ccs_values = [float(row["peak_ccs"]) for row in analytes]
        if not all(math.isfinite(value) and value > 0 for value in ccs_values):
            raise RuntimeError("All 34 ADP software CCS results must be finite and positive")
        reference_times = [float(row["drift_time_ms"]) for row in references]
        peak_times = [float(row["peak_drift_time_ms"]) for row in analytes]
        safe = lambda name: name.replace("/", "_").replace(" ", "_")
        pngs = [ccs_dir / "calibration_plots.png"]
        pngs += [ccs_dir / "cali_png" / f"{Path(row['raw_path']).stem}_{safe(row['name'])}.png"
                 for row in expected_references]
        pngs += [ccs_dir / folder / f"{row['run']}_{safe(row['species'])}_ccs.png" for row in analytes
                 for folder in ("png_raw", "png_smoothed", "png_overlay")]
        if len(pngs) != 128 or any(not p.is_file() or not p.read_bytes().startswith(b"\x89PNG\r\n\x1a\n") for p in pngs):
            raise RuntimeError("ADP CCS smoke did not publish all 128 calibration/analyte PNGs")
        report["ccs_smoke"] = {"method": "direct", "profile": profile["id"], "raw_runs": 5,
                               "calibrants": 25, "analytes": 34,
                               "analyte_runs": sorted({row["run"] for row in analytes}),
                               "ccs_angstrom_squared_range": [min(ccs_values), max(ccs_values)],
                               "analyte_drift_range_ms": [min(peak_times), max(peak_times)],
                               "calibrant_drift_range_ms": [min(reference_times), max(reference_times)],
                               "extrapolated_analytes": sum(value < min(reference_times) or value > max(reference_times)
                                                            for value in peak_times),
                               "scientific_acceptance": "review_required_poor_fit_and_extrapolation_not_validated",
                               "log_fit_r2": results["calibration"]["r2_lnln"],
                               "pngs": [p.relative_to(ccs_dir).as_posix() for p in pngs]}
        # The Betaine detail view is a separate visual-only acquisition. Its
        # converted pair must also reopen without re-reading the RAW or CCS.
        profile = next(value for value in manifest["profiles"] if value["id"] == "betaine_1000_2000")
        raw = production / profile["source"]
        config_path = production / profile["viewer_config"]
        config = json.loads(config_path.read_text(encoding="utf-8"))
        existing_html = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                         for p in (output / "betaine_demo").glob("*.html")}
        converted = process(raw, output / "betaine_demo", html=False)
        if (len(converted["runs"]) != 1 or converted["viewers"] or converted["ccs_dir"]
                or existing_html != {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                     for p in (output / "betaine_demo").glob("*.html")}):
            raise RuntimeError("Betaine convert-only must publish one MS/IM pair without HTML or CCS")
        report["convert_only"] = {"run": raw.stem, "ms_im_pairs": 1, "viewers": 0, "ccs": False}
        report["preview_profile"] = profile["id"]
        report["preview_run"] = raw.stem
        report["raw_reader"] = converted["runs"][0]

        def wait_until(predicate, seconds=90):
            deadline = time.monotonic() + seconds
            while time.monotonic() < deadline:
                application.processEvents()
                if predicate():
                    return
                time.sleep(0.02)
            raise RuntimeError("Native GUI did not reach the expected state")

        window = MainWindow()
        window.resize(1280, 800)
        window.show()
        application.processEvents()
        window.resize(1280, 800)
        application.processEvents()
        if window.width() != 1280 or window.height() != 800:
            raise RuntimeError(f"Requested 1280x800 client, obtained {window.width()}x{window.height()}")
        actual_font = QFontInfo(application.font())
        report["native_gui_fonts"] = {"family_count": len(QFontDatabase.families()),
                                      "actual_family": actual_font.family(), "pixel_size": actual_font.pixelSize()}
        if not report["native_gui_fonts"]["family_count"] or not actual_font.family().strip() or actual_font.pixelSize() <= 0:
            raise RuntimeError("Qt could not load usable desktop fonts")
        panel = window.get_panel(0)
        panel.combo_demo.setCurrentIndex(next(index for index in range(panel.combo_demo.count())
            if (panel.combo_demo.itemData(index) or {}).get("id") == profile["id"]))
        panel._load_demo()
        for picker, expected_path in ((panel.calibration_panel.pick_csv, calibrant_path),
                                      (panel.analyte_panel.pick_ana, analyte_path)):
            picker.set_path(f'"{expected_path}"')
            if Path(picker.path()).resolve() != expected_path.resolve():
                raise RuntimeError("The GUI did not preserve an Explorer-quoted CSV path")
        report["quoted_csv_tables"] = {"calibrants": panel.calibration_panel.table.table.rowCount(),
                                       "analytes": panel.analyte_panel.table.table.rowCount()}
        if report["quoted_csv_tables"] != {"calibrants": 25, "analytes": 34}:
            raise RuntimeError("Quoted CSV paths did not load all 25 calibrants and 34 analytes")
        panel._load_demo()
        for picker, expected_path in ((panel.pick_source, raw), (panel.pick_out, output / "betaine_demo")):
            picker.set_path(f'"{expected_path}"')
            if Path(picker.path()).resolve() != expected_path.resolve():
                raise RuntimeError("The GUI did not preserve an Explorer-quoted source/output path")
        report["quoted_gui_paths"] = {"source_output_and_csv": "passed",
                                      "bundled_raw_path_contains_space": " " in str(raw),
                                      "output_path_contains_space": " " in str(output)}
        panel.reload_runs()
        panel.run_table.selectRow(next(index for index, row in enumerate(panel._rows) if row["run_name"] == raw.stem))
        viewer = panel.viewer_widget
        record = {"run_name": raw.stem,
                  "ms": Path(converted["converted_dir"]) / f"{raw.stem}_ms.txt",
                  "im": Path(converted["converted_dir"]) / f"{raw.stem}_im.txt"}
        render_events, export_events = [], []
        viewer.render_finished.connect(render_events.append)
        viewer.export_finished.connect(export_events.append)
        viewer.set_run(record, config, raw)
        wait_until(lambda: viewer.is_ready and viewer.snapshot is not None)
        import numpy as np
        from deconvovo.io import find_pusher_period
        acquired = np.loadtxt(record["im"])
        positive = acquired[acquired[:, 2] > 0]
        native_drift_count = int(positive[:, 1].max()) + 1
        raw_mz_domain = [float(positive[:, 0].min()), float(positive[:, 0].max())]
        pusher_us = find_pusher_period(raw)
        if not pusher_us or pusher_us <= 0:
            raise RuntimeError("Native alignment requires the actual RAW pusher period")

        def measure_alignment(settings, label, viewport=None):
            viewer.canvas.draw()
            snapshot = viewer.snapshot
            configured = settings["initial_view"]
            requested = configured.get("mz_range") or raw_mz_domain
            domain = [max(raw_mz_domain[0], requested[0]), min(raw_mz_domain[1], requested[1])]
            ranges = viewport or configured
            bins = settings["defaults"]["mz_bins"]
            if bins == "native":
                bins = min(25600, max(200, len(np.unique(positive[
                    (positive[:, 0] >= domain[0]) & (positive[:, 0] <= domain[1]), 0]))))
            drift_bins = settings["defaults"].get("drift_bins", "native")
            drift_bins = native_drift_count if drift_bins == "native" else int(drift_bins)
            upper_drift = native_drift_count * pusher_us / 1000
            expected_mz = domain[0] + (np.arange(int(bins)) + 0.5) * (domain[1] - domain[0]) / int(bins)
            expected_drift = (np.arange(drift_bins) + 0.5) * upper_drift / drift_bins
            def centres(actual, expected):
                error = float(np.max(np.abs(np.asarray(actual) - expected))) if len(actual) == len(expected) else math.inf
                return {"count": len(actual), "expected_count": len(expected), "max_error": error}
            measured = {"expected_from_raw_and_preset": {"mz_domain": domain, "mz_range": ranges["mz_range"],
                         "drift_range": ranges["drift_range"], "raw_pusher_us": pusher_us},
                        "heatmap_mz_centres": centres(snapshot.mz_centres, expected_mz),
                        "heatmap_drift_centres": centres(snapshot.drift_centres, expected_drift),
                        "profile_drift_centres": centres(viewer.ax_drift.lines[0].get_ydata(), expected_drift)}
            extent = list(map(float, viewer.ax_heatmap.images[0].get_extent()))
            measured["image_extent"] = extent
            if not np.allclose(extent, [*domain, 0, upper_drift], rtol=0, atol=1e-7):
                raise RuntimeError("Native heatmap pixels do not use the independent RAW-derived cell edges")
            for name, profile_axis, bounds, direction in (
                    ("mz", viewer.ax_ms, ranges["mz_range"], 0),
                    ("drift", viewer.ax_drift, ranges["drift_range"], 1)):
                heat = viewer.ax_heatmap
                samples = []
                for fraction in (0, .25, .5, .75, 1):
                    value = bounds[0] + fraction * (bounds[1] - bounds[0])
                    point = [value, 0] if direction == 0 else [0, value]
                    expected = (heat.bbox.x0 + fraction * heat.bbox.width if direction == 0
                                else heat.bbox.y0 + fraction * heat.bbox.height)
                    hp = float(heat.transData.transform(point)[direction])
                    pp = float(profile_axis.transData.transform(point)[direction])
                    samples.append({"coordinate": value, "expected_pixel": expected, "heatmap_pixel": hp,
                                    "profile_pixel": pp, "error_px": max(abs(hp - expected), abs(pp - expected))})
                measured[name] = {"samples": samples, "max_error_px": max(row["error_px"] for row in samples)}
                if measured[name]["max_error_px"] > 1:
                    raise RuntimeError(f"Native {name} profile differs from its physical pixel projection")
            for key in ("heatmap_mz_centres", "heatmap_drift_centres", "profile_drift_centres"):
                if measured[key]["max_error"] > 1e-7:
                    raise RuntimeError(f"{label} {key} differs from RAW-derived histogram centres")
            report.setdefault("profile_alignment", {})[label] = measured

        def useful_chart():
            width, height = viewer.canvas.width(), viewer.canvas.height()
            viewer.canvas.draw()
            axes_visible = all(ax.bbox.x0 >= 0 and ax.bbox.y0 >= 0 and ax.bbox.x1 <= width + 1 and ax.bbox.y1 <= height + 1
                               for ax in (viewer.ax_heatmap, viewer.ax_drift, viewer.ax_ms))
            report["console_chart_at_1280x800"] = {"renderer": "Matplotlib QtAgg", "width": width,
                                                  "height": height, "fully_visible": axes_visible}
            return width >= 300 and height >= 250 and axes_visible
        wait_until(useful_chart, 15)
        measure_alignment(config, "initial_detail")
        repaint_until = time.monotonic() + .75
        while time.monotonic() < repaint_until:
            application.processEvents(); time.sleep(.02)
        window.grab().save(str(output / "native_console.png"))
        changed = deepcopy(config)
        changed["defaults"].update({"mz_bins": 250, "noise_2d": 7.25, "ms_mode": "stick"})
        changed["initial_view"] = {"mz_range": [1100, 1800], "drift_range": [5, 10]}
        previous = len(render_events)
        viewer.apply_config(changed)
        wait_until(lambda: len(render_events) > previous and viewer.is_ready)
        measure_alignment(changed, "after_live_config_rebin")
        viewport = {"mz_range": [1200, 1600], "drift_range": [6, 9]}
        viewer.ax_heatmap.set_xlim(viewport["mz_range"])
        viewer.ax_heatmap.set_ylim(viewport["drift_range"])
        wait_until(lambda: viewer.is_ready and np.allclose(viewer.snapshot.mz_range, viewport["mz_range"])
                   and np.allclose(viewer.snapshot.drift_range, viewport["drift_range"]))
        measure_alignment(changed, "after_mouse_zoom", viewport)
        changed["defaults"]["noise_2d"] = 8
        captures = output / "current view exports with spaces"
        captures.mkdir(exist_ok=True)
        source_hashes = {kind: hashlib.sha256(record[kind].read_bytes()).hexdigest() for kind in ("ms", "im")}
        viewer.apply_config(changed)
        exported = {}
        for kind in ("full_png", "heatmap_png", "mz_png", "drift_png", "heatmap_csv", "mz_csv", "drift_csv"):
            destination = captures / (raw.stem + "_current_" + kind.replace("_", "."))
            previous = len(export_events)
            returned = viewer.export_current(kind, destination)
            if Path(returned).resolve() != destination.resolve():
                raise RuntimeError("Native export did not return the requested destination")
            wait_until(lambda: len(export_events) > previous and any(
                Path(event["path"]).resolve() == destination.resolve() for event in export_events[previous:]), 30)
            event = next(event for event in export_events[previous:] if Path(event["path"]).resolve() == destination.resolve())
            if event.get("error") or not destination.is_file() or destination.stat().st_size == 0:
                raise RuntimeError(f"Native current-view export failed: {event}")
            exported[kind] = destination
        measure_alignment(changed, "after_control_preserving_zoom", viewport)
        snapshot = viewer.snapshot
        with exported["heatmap_csv"].open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames != ["mz", "drift_time_ms", "display_value", "intensity_transform"]:
                raise RuntimeError("Native heatmap CSV must label coordinates and displayed transform")
            rows = list(reader)
        expected = [(float(mz), float(drift), float(value)) for drift, row in zip(
            snapshot.drift_centres[snapshot.visible_drift_slice],
            snapshot.display2d[snapshot.visible_drift_slice, snapshot.visible_mz_slice])
            for mz, value in zip(snapshot.mz_centres[snapshot.visible_mz_slice], row)]
        actual = [(float(row["mz"]), float(row["drift_time_ms"]), float(row["display_value"])) for row in rows]
        if actual != expected or not all(row["intensity_transform"] == snapshot.transform for row in rows):
            raise RuntimeError("Native heatmap CSV differs from the displayed coordinates and values")
        csv_counts = {"heatmap": len(rows)}
        for kind, column, coordinates, values, bounds in (
                ("mz_csv", "mz", snapshot.ms_mz, snapshot.ms_intensity, snapshot.mz_range),
                ("drift_csv", "drift_time_ms", snapshot.drift_centres, snapshot.drift_intensity, snapshot.drift_range)):
            with exported[kind].open(encoding="utf-8", newline="") as handle:
                reader = csv.DictReader(handle)
                if reader.fieldnames != [column, "display_value"]:
                    raise RuntimeError(f"Native {kind} has incorrect column labels")
                rows = list(reader)
            actual = [(float(row[column]), float(row["display_value"])) if row[column] and row["display_value"]
                      else (None, None) for row in rows]
            expected = [(float(x), float(y)) if np.isfinite(x) and np.isfinite(y) else (None, None)
                        for x, y in zip(coordinates, values) if not np.isfinite(x) or not np.isfinite(y) or bounds[0] <= x <= bounds[1]]
            if actual != expected:
                raise RuntimeError(f"Native {kind} differs from its displayed trace")
            csv_counts[kind.removesuffix("_csv")] = len(rows)
        images = {}
        for kind in ("full_png", "heatmap_png", "mz_png", "drift_png"):
            image = QImage(str(exported[kind]))
            if image.isNull() or image.width() < 1000 or image.height() < 700:
                raise RuntimeError(f"Native {kind} is not a valid rendered PNG")
            images[kind] = {"width": image.width(), "height": image.height()}
        for destination in (record["ms"], raw / "protected current view.csv"):
            try:
                viewer.export_current("heatmap_csv", destination)
            except ValueError:
                pass
            else:
                raise RuntimeError("Native current-view export failed to reject a protected source location")
        if any(hashlib.sha256(record[kind].read_bytes()).hexdigest() != digest for kind, digest in source_hashes.items()):
            raise RuntimeError("Native exports changed original MS/IM data")
        report["current_view_exports"] = {"formats": sorted(exported), "completed": 7, "csv_rows": csv_counts,
                                           "csv_matches_display": True, "png_dimensions": images,
                                           "source_guard": "passed", "completion": "export_finished"}
        report["native_viewer"] = {"renderer": type(viewer.canvas).__module__, "ready": viewer.is_ready,
                                    "needs_html_to_render": False, "drift_unit": snapshot.drift_unit,
                                    "mouse_zoom_preserved_after_controls": True}
        # The explicit workspace action generates HTML from TXT on its worker;
        # native plotting above has never loaded or interpreted an HTML file.
        standalone = output / "standalone_copy"
        standalone.mkdir(exist_ok=True)
        generated = standalone / "exported offline viewer.html"
        returned = panel.export_portable_html(generated)
        if returned is None:
            raise RuntimeError("Explicit portable HTML export could not start")
        wait_until(lambda: not panel._busy and generated.is_file(), 120)
        contents = generated.read_text(encoding="utf-8")
        if re.search(r'<(?:script|link|img)\b[^>]*(?:src|href)\s*=\s*[\"\']https?://', contents, re.I):
            raise RuntimeError("Portable HTML references an external resource")
        if "plotly.js" not in contents or "Plotly.newPlot" not in contents:
            raise RuntimeError("Portable HTML is missing embedded Plotly")
        if not all(f'id="{name}"' in contents for name in ("nbins", "dbins", "msmode", "noise2d", "denoise", "sm2d", "cmap", "zscale")):
            raise RuntimeError("Portable HTML is missing offline interactive controls")
        copied = standalone / "renamed offline viewer.html"
        shutil.copyfile(generated, copied)
        digest = hashlib.sha256(generated.read_bytes()).hexdigest()
        if hashlib.sha256(copied.read_bytes()).hexdigest() != digest:
            raise RuntimeError("Moved standalone HTML bytes changed")
        report["portable_html_bytes"] = generated.stat().st_size
        report["standalone_copy"] = {"sha256": digest, "bytes": copied.stat().st_size,
                                      "embedded_plotly": True, "offline_controls": 8, "native_browser_runtime": False}
        browser_modules = [name for name in sys.modules if "qtwebengine" in name.casefold()]
        if os.name == "nt":
            command = ["powershell.exe", "-NoProfile", "-Command",
                       "Get-CimInstance Win32_Process | Select-Object ProcessId,ParentProcessId,Name | ConvertTo-Json -Compress"]
            reply = subprocess.run(command, capture_output=True, text=True, check=True)
            processes = [(int(row["ProcessId"]), int(row["ParentProcessId"]), row["Name"])
                         for row in json.loads(reply.stdout)]
        else:
            reply = subprocess.run(["ps", "-eo", "pid=,ppid=,comm="], capture_output=True, text=True, check=True)
            processes = [(int(parts[0]), int(parts[1]), parts[2]) for line in reply.stdout.splitlines()
                         if len(parts := line.strip().split(None, 2)) == 3]
        family = {os.getpid()}
        while True:
            expanded = family | {pid for pid, parent, _name in processes if parent in family}
            if expanded == family:
                break
            family = expanded
        browser_processes = [name for pid, _parent, name in processes if pid in family and
                             any(word in name.casefold() for word in ("webengine", "chromium", "chrome", "msedge"))]
        report["browser_free_gui"] = {"loaded_webengine_modules": browser_modules,
                                      "browser_child_processes": browser_processes}
        if browser_modules or browser_processes:
            raise RuntimeError("Native plotting loaded a browser runtime or launched a browser process")
        report["status"] = "passed"
        return_code = 0
    except Exception:
        report["status"] = "failed"
        report["error"] = traceback.format_exc()
        return_code = 1
    finally:
        if window is not None:
            try:
                repaint_until = time.monotonic() + .75
                while time.monotonic() < repaint_until:
                    application.processEvents(); time.sleep(.02)
                if viewer is not None:
                    viewer.canvas.draw()
                final_image = output / "native_console_final_state.png"
                captured = window.grab()
                report["final_screenshot"] = {"path": final_image.name, "saved": captured.save(str(final_image), "PNG")}
                if report["status"] != "passed":
                    captured.save(str(output / "native_console.png"), "PNG")
                report["qt_layout_diagnostics"] = _qt_diagnostics(application, window, panel, viewer)
            except Exception:
                report["diagnostics_error"] = traceback.format_exc()
        if viewer is not None:
            viewer.wait_for_workers()
            viewer.close()
        if window is not None:
            window.close()
        if application is not None:
            application.processEvents()
        (output / "native_smoke.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return return_code
