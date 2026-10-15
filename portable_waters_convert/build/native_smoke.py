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
        # Validate the same desktop platform and fonts used by the native GUI.
        if os.name == "nt":
            # Windows' offscreen plugin has no system font database. Validate
            # the same native platform and fonts used by the desktop app.
            os.environ["QT_QPA_PLATFORM"] = "windows"
        from gui.widgets.worker import set_gui_mode
        set_gui_mode()
        from PySide6.QtWidgets import QApplication, QComboBox, QPushButton
        from PySide6.QtCore import QEvent, QPointF, QPoint, Qt
        from PySide6.QtGui import QFont, QFontDatabase, QFontInfo, QImage, QMouseEvent, QWheelEvent
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
        config_path = production / profile["viewer_config"]
        config = json.loads(config_path.read_text(encoding="utf-8"))
        calibrant_path = production / profile["calibrant_csv"]
        analyte_path = production / profile["analyte_csv"]
        results = process(production / profile["source"], output / "real_demo", html=False, jobs=2, config_path=config_path,
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
        if (not panel.btn_open_raw.text().casefold().startswith("open data")
                or any(button.isVisible() for button in (panel.btn_run, panel.btn_convert, panel.btn_text))):
            raise RuntimeError("The native console must expose one neutral Open data action")
        report["open_data_action"] = {"primary": panel.btn_open_raw.text(), "legacy_actions_prominent": False}
        example_counts = {}
        default_output = panel._default_output
        try:
            for example_id, count in (("protein_adp", 5), ("small_molecule", 2), ("betaine_1000_2000", 1)):
                queued_output = output / "queued examples" / example_id
                panel._default_output = lambda name, destination=queued_output: destination
                index = next(index for index in range(panel.combo_demo.count())
                    if (panel.combo_demo.itemData(index) or {}).get("id") == example_id)
                panel.combo_demo.setCurrentIndex(index)
                application.processEvents()
                if panel.is_busy or len(panel._rows) != count or (queued_output / "_converted").exists():
                    raise RuntimeError("Selecting an Examples family must list its RAWs without converting them")
                example_counts[example_id] = count
        finally:
            panel._default_output = default_output
        if not panel.btn_process.isVisible() or not panel.btn_reload.isVisible() or not panel.spin_jobs.isVisible():
            raise RuntimeError("Data must expose Process files, Refresh, and the CPU limit")
        panel.spin_jobs.setValue(2)
        report["data_queue"] = {"example_family_counts": example_counts, "conversion_requires_process": True,
                                "process_refresh_cpu_controls_visible": True, "cpu_limit": panel.spin_jobs.value()}
        refresh_source = output / "refresh source with spaces"
        refresh_source.mkdir(exist_ok=True)
        panel.pick_source.set_path(str(refresh_source))
        refreshed_raw = refresh_source / "added betaine.raw"
        shutil.copytree(raw, refreshed_raw)
        panel.btn_reload.click()
        if len(panel._rows) != 1 or panel._rows[0]["run_name"] != refreshed_raw.stem or panel.is_busy:
            raise RuntimeError("Refresh did not discover an added complete RAW folder without converting it")
        state_header = panel.run_table.horizontalHeader()
        if (panel.run_table.isColumnHidden(4) or state_header.sectionViewportPosition(4) < 0
                or state_header.sectionViewportPosition(4) + state_header.sectionSize(4) > panel.run_table.viewport().width()
                or state_header.sectionResizeMode(4) != state_header.ResizeMode.Interactive):
            raise RuntimeError("Acquisition State must be visible and user-resizable")
        shutil.rmtree(refreshed_raw)
        panel.btn_reload.click()
        if panel._rows or panel.is_busy:
            raise RuntimeError("Refresh did not remove a deleted RAW folder from the file list")
        report["data_queue"]["refresh_added_deleted_real_raw"] = "passed"
        report["data_queue"]["state_column_visible_and_resizable"] = True
        demo_index = next(index for index in range(panel.combo_demo.count())
            if (panel.combo_demo.itemData(index) or {}).get("id") == profile["id"])
        panel.combo_demo.setCurrentIndex(demo_index)
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
        custom_output = output / "custom RAW input with spaces"
        default_output = panel._default_output
        panel._default_output = lambda name: custom_output
        try:
            panel.pick_source.set_path("")
            panel.open_raw(f'"{raw}"')
        finally:
            panel._default_output = default_output
        wait_until(lambda: not panel._busy, 120)
        opened = panel._last_result
        if (opened.get("status") != "complete" or len(opened.get("runs", [])) != 1
                or opened.get("viewers") or opened.get("ccs_dir")
                or Path(opened["converted_dir"]).resolve() != (custom_output / "_converted").resolve()):
            raise RuntimeError("Open RAW did not automatically convert the custom input for native viewing")
        wait_until(lambda: panel.viewer_widget.is_ready and
                   panel.viewer_widget._record.get("run_name") == raw.stem and
                   Path(panel.viewer_widget._record.get("im", "")).resolve() ==
                   (custom_output / "_converted" / f"{raw.stem}_im.txt").resolve())
        custom_config = panel._viewer_config()
        if (custom_config.get("initial_view") != {"mz_range": None, "drift_range": None}
                or panel.chk_ccs.isChecked() or panel.calibration_panel.table.table.rowCount()
                or panel.analyte_panel.table.table.rowCount()):
            raise RuntimeError("Open RAW retained example ranges, CCS selection, or example CSV tables")
        report["custom_raw_gui"] = {"run": raw.stem, "automatically_processed": True,
                                     "native_ready": True, "html_outputs": 0, "ccs_enabled": False,
                                     "example_ranges_cleared": True, "example_tables_cleared": True}
        panel.combo_demo.setCurrentIndex(demo_index)
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
                  "raw": raw,
                  "ms": Path(converted["converted_dir"]) / f"{raw.stem}_ms.txt",
                  "im": Path(converted["converted_dir"]) / f"{raw.stem}_im.txt"}
        render_events, export_events = [], []
        viewer.render_finished.connect(render_events.append)
        viewer.export_finished.connect(export_events.append)
        viewer.set_run(record, config, raw)
        wait_until(lambda: viewer.is_ready and viewer.snapshot is not None)
        import numpy as np
        import pandas as pd
        from deconvovo.io import find_pusher_period
        acquired = np.loadtxt(record["im"])
        acquired_ms = np.loadtxt(record["ms"], ndmin=2)
        if (acquired[:, 1] < 0).any() or not np.array_equal(acquired[:, 1], np.floor(acquired[:, 1])):
            raise RuntimeError("Original acquired drift indices must be nonnegative integers")
        positive = acquired[acquired[:, 2] > 0]
        native_drift_count = int(acquired[:, 1].max()) + 1
        acquired_mz_axis = np.union1d(acquired_ms[:, 0], acquired[:, 0])
        raw_mz_domain = [float(acquired_mz_axis.min()), float(acquired_mz_axis.max())]
        pusher_us = find_pusher_period(raw)
        if not pusher_us or pusher_us <= 0:
            raise RuntimeError("Native alignment requires the actual RAW pusher period")

        def measure_alignment(settings, label, viewport=None, smoothing_context=None):
            viewer.canvas.draw()
            snapshot = viewer.snapshot
            configured = settings["initial_view"]
            requested = configured.get("mz_range") or raw_mz_domain
            domain = [max(raw_mz_domain[0], requested[0]), min(raw_mz_domain[1], requested[1])]
            if smoothing_context is not None:
                context = list(map(float, smoothing_context))
                if not (raw_mz_domain[0] <= context[0] <= domain[0] < domain[1] <= context[1] <= raw_mz_domain[1]):
                    raise RuntimeError("The smoothing context must enclose the exact target inside the acquired domain")
                domain = context
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
            expected_profile_drift = (np.arange(native_drift_count) + 0.5) * pusher_us / 1000
            def centres(actual, expected):
                error = float(np.max(np.abs(np.asarray(actual) - expected))) if len(actual) == len(expected) else math.inf
                return {"count": len(actual), "expected_count": len(expected), "max_error": error}
            measured = {"expected_from_raw_and_preset": {"mz_domain": domain, "mz_range": ranges["mz_range"],
                         "drift_range": ranges["drift_range"], "raw_pusher_us": pusher_us},
                        "heatmap_mz_centres": centres(snapshot.mz_centres, expected_mz),
                        "heatmap_drift_centres": centres(snapshot.drift_centres, expected_drift),
                        "profile_drift_centres": centres(snapshot.drift_profile_centres, expected_profile_drift),
                        "rendered_profile_drift_centres": centres(viewer.ax_drift.lines[0].get_ydata(), expected_profile_drift)}
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
            for key in ("heatmap_mz_centres", "heatmap_drift_centres", "profile_drift_centres", "rendered_profile_drift_centres"):
                if measured[key]["max_error"] > 1e-7:
                    raise RuntimeError(f"{label} {key} differs from RAW-derived histogram centres")
            report.setdefault("profile_alignment", {})[label] = measured

        def measure_native_profiles(label, original_im=acquired, original_ms=acquired_ms, period=pusher_us):
            snapshot = viewer.snapshot
            defaults = viewer._config["defaults"]
            points_original = original_im[original_im[:, 2] > 0]
            count = int(original_im[:, 1].max()) + 1
            mass_axis = np.union1d(original_ms[:, 0], original_im[:, 0])
            if (snapshot.ms_mode != "profile" or defaults.get("noise_mz", 0) != 0
                    or defaults.get("noise_drift", 0) != 0
                    or defaults.get("smooth_drift", {}).get("method", "raw") != "raw"):
                raise RuntimeError("Acquired-point acceptance requires explicit raw 1D profile settings")
            lo, hi = snapshot.mz_range
            mz_selected = (points_original[:, 0] >= lo) & (points_original[:, 0] <= hi)
            drift_values = np.bincount(points_original[mz_selected, 1].astype(int),
                weights=points_original[mz_selected, 2], minlength=count)
            drift_axis = (np.arange(count) + .5) * (period / 1000)
            low, high = snapshot.drift_range
            step = period / 1000
            edges = np.arange(count + 1) * step
            epsilon = 8 * np.finfo(float).eps * max(1, abs(low), abs(high), abs(edges[-1]))
            selected_bins = (edges[1:] > max(0, low) + epsilon) & (edges[:-1] < min(edges[-1], high) - epsilon)
            points = points_original[selected_bins[points_original[:, 1].astype(int)]]
            mass_values = np.bincount(np.searchsorted(mass_axis, points[:, 0]),
                weights=points[:, 2], minlength=len(mass_axis))
            mass_selected = (mass_axis >= lo) & (mass_axis <= hi)
            expected_axis, expected_values = mass_axis[mass_selected], mass_values[mass_selected]
            for actual, expected in (
                    (snapshot.ms_mz, expected_axis), (snapshot.ms_intensity, expected_values),
                    (viewer.ax_ms.lines[0].get_xdata(), expected_axis),
                    (viewer.ax_ms.lines[0].get_ydata(), expected_values),
                    (snapshot.drift_profile_centres, drift_axis), (snapshot.drift_intensity, drift_values),
                    (viewer.ax_drift.lines[0].get_ydata(), drift_axis),
                    (viewer.ax_drift.lines[0].get_xdata(), drift_values),
                    (viewer.mass_plot.ax().lines[0].get_xdata(), expected_axis),
                    (viewer.mobility_plot.ax().lines[0].get_xdata(), drift_axis)):
                if not np.array_equal(actual, expected):
                    raise RuntimeError(f"{label}: native 1D profile differs from original acquired coordinates or sums")
            report.setdefault("native_acquired_profiles", {})[label] = {
                "acquired_mz_coordinates": len(mass_axis), "displayed_mz_coordinates": len(expected_axis),
                "native_drift_bins": count, "heatmap_drift_bins": len(snapshot.drift_centres),
                "zero_only_native_drift_bins": int(np.count_nonzero(np.bincount(points_original[:, 1].astype(int),
                    minlength=count) == 0)), "exact_original_coordinates_and_values": True}

        def useful_chart():
            width, height = viewer.canvas.width(), viewer.canvas.height()
            viewer.canvas.draw()
            physical_width, physical_height = viewer.canvas.get_width_height(physical=True)
            axes_visible = all(ax.bbox.x0 >= 0 and ax.bbox.y0 >= 0 and ax.bbox.x1 <= physical_width + 1 and ax.bbox.y1 <= physical_height + 1
                               for ax in (viewer.ax_heatmap, viewer.ax_drift, viewer.ax_ms))
            report["console_chart_at_1280x800"] = {"renderer": "Matplotlib QtAgg", "width": width,
                                                  "height": height, "physical_width": physical_width,
                                                  "physical_height": physical_height,
                                                  "fully_visible": axes_visible}
            return width >= 300 and height >= 250 and axes_visible
        wait_until(useful_chart, 15)
        measure_alignment(config, "initial_detail")
        measure_native_profiles("initial_detail")
        if native_drift_count != 143 or len(viewer.snapshot.ms_mz) != 39881:
            raise RuntimeError("Real Betaine detail must retain all 143 drift bins and 39,881 acquired m/z samples")
        repaint_until = time.monotonic() + .75
        while time.monotonic() < repaint_until:
            application.processEvents(); time.sleep(.02)
        window.grab().save(str(output / "native_console.png"))
        original_profiles = viewer.snapshot
        coarse = deepcopy(config)
        coarse["defaults"].update({"mz_bins": 80, "drift_bins": 64})
        previous = len(render_events)
        panel.viewer_panel.spin_mzbins.setValue(80)
        panel.viewer_panel.spin_driftbins.setValue(64)
        wait_until(lambda: len(render_events) > previous and viewer.is_ready)
        measure_alignment(coarse, "coarse_2d_native_1d")
        measure_native_profiles("coarse_2d_native_1d")
        for name in ("ms_mz", "ms_intensity", "drift_profile_centres", "drift_intensity"):
            if not np.array_equal(getattr(viewer.snapshot, name), getattr(original_profiles, name)):
                raise RuntimeError("Coarse heatmap bins changed acquired 1D samples or values")
        previous = len(render_events)
        panel.viewer_panel.spin_mzbins.setValue(config["defaults"]["mz_bins"])
        panel.viewer_panel.spin_driftbins.setValue(0)
        wait_until(lambda: len(render_events) > previous and viewer.is_ready)
        # Measure actual Agg paths on a full native drift viewport, avoiding
        # legitimate clipping of sorted lines outside a mouse zoom.
        viewer.ax_heatmap.set_ylim(0, native_drift_count * pusher_us / 1000)
        wait_until(lambda: viewer.is_ready and np.allclose(viewer.snapshot.drift_range,
                   [0, native_drift_count * pusher_us / 1000]))
        measure_native_profiles("full_native_drift_render")
        from matplotlib.backends.backend_agg import RendererAgg
        original_draw_path = RendererAgg.draw_path
        rendered_profiles = []

        def capture_profile_path(renderer, gc, path, affine, rgbFace=None):
            if len(path.vertices) > 100:
                rendered_profiles.append((path.vertices.copy(), path.should_simplify, renderer.dpi))
            return original_draw_path(renderer, gc, path, affine, rgbFace)

        RendererAgg.draw_path = capture_profile_path
        try:
            full_snapshot = viewer.snapshot
            expected_paths = (
                np.column_stack((full_snapshot.ms_mz, full_snapshot.ms_intensity)),
                np.column_stack((full_snapshot.drift_intensity, full_snapshot.drift_profile_centres)),
                np.column_stack((full_snapshot.drift_profile_centres, full_snapshot.drift_intensity)))
            for canvas, relevant in ((viewer.canvas, expected_paths[:2]),
                                     (viewer.mass_plot.canvas, expected_paths[:1]),
                                     (viewer.mobility_plot.canvas, expected_paths[2:])):
                rendered_profiles.clear()
                canvas.draw()
                for expected in relevant:
                    matches = [simplify for vertices, simplify, _dpi in rendered_profiles
                               if vertices.shape == expected.shape and np.array_equal(vertices, expected)]
                    if not matches or any(matches):
                        raise RuntimeError("A native canvas omitted or simplified complete acquired 1D data")
            render_exports = output / "full resolution render proof"
            render_exports.mkdir(exist_ok=True)
            for kind in ("full_png", "mz_png", "drift_png"):
                rendered_profiles.clear()
                export_dpi = viewer.figure.dpi * 3
                viewer.export_current(kind, render_exports / (kind + ".png"))
                relevant = expected_paths[:2] if kind == "full_png" else (
                    expected_paths[0] if kind == "mz_png" else expected_paths[2],)
                for expected in relevant:
                    matches = [simplify for vertices, simplify, dpi in rendered_profiles
                               if dpi == export_dpi and vertices.shape == expected.shape
                               and np.array_equal(vertices, expected)]
                    if not matches or any(matches):
                        raise RuntimeError("The actual PNG renderer did not preserve an acquired profile path")
        finally:
            RendererAgg.draw_path = original_draw_path
        report["native_profile_rendering"] = {"complete_acquired_paths": True,
            "gui_and_png_simplification_disabled": True, "native_mz_samples": len(full_snapshot.ms_mz),
            "native_drift_samples": len(full_snapshot.drift_profile_centres)}
        panel._set_viewer_config(config)
        panel._apply_viewer_settings()
        viewer.ax_heatmap.set_ylim(config["initial_view"]["drift_range"])
        wait_until(lambda: viewer.is_ready and np.allclose(viewer.snapshot.drift_range, config["initial_view"]["drift_range"]))
        # Exercise real Qt input, rather than only setting axes through Python.
        # A drag lasts longer than the render debounce, exposing redraws that
        # incorrectly clear the axes while a user's mouse button is held.
        def pointer(fractions):
            viewer.canvas.draw()
            x, y = viewer.ax_heatmap.transAxes.transform(fractions)
            ratio = viewer.canvas.device_pixel_ratio
            local = QPointF(x / ratio, viewer.canvas.height() - y / ratio)
            global_point = viewer.canvas.mapToGlobal(local.toPoint())
            return local, QPointF(global_point)

        def mouse(kind, point, button, buttons):
            local, global_point = point
            QApplication.sendEvent(viewer.canvas, QMouseEvent(kind, local, global_point, button, buttons,
                                                               Qt.KeyboardModifier.NoModifier))

        def limits():
            return {"mz_range": list(viewer.snapshot.mz_range), "drift_range": list(viewer.snapshot.drift_range)}

        def range_fields_match(label):
            measured = {}
            for key, controls in (
                    ("mz_range", (panel.viewer_panel.spin_mz_min, panel.viewer_panel.spin_mz_max)),
                    ("drift_range", (panel.viewer_panel.spin_drift_min, panel.viewer_panel.spin_drift_max))):
                values = [control.value() for control in controls]
                expected = limits()[key]
                if any(abs(value - bound) > 0.5 * 10 ** -control.decimals() + 1e-9
                       for value, bound, control in zip(values, expected, controls)):
                    raise RuntimeError(f"{label}: numeric {key} fields do not match the native viewport")
                measured[key] = {"fields": values, "viewport": expected}
            report.setdefault("dynamic_range_fields", {})[label] = measured

        before = limits()
        range_fields_match("initial_detail")
        local, global_point = pointer((.5, .5))
        mouse(QEvent.Type.MouseMove, (local, global_point), Qt.MouseButton.NoButton, Qt.MouseButton.NoButton)
        hover = viewer.toolbar.locLabel.text()
        if "m/z" not in hover or "ms" not in hover or "cbrt(I)" not in hover:
            raise RuntimeError("Native hover did not show physical heatmap coordinates and its displayed value")
        previous = len(render_events)
        QApplication.sendEvent(viewer.canvas, QWheelEvent(local, global_point, QPoint(), QPoint(0, 120),
            Qt.MouseButton.NoButton, Qt.KeyboardModifier.NoModifier, Qt.ScrollPhase.NoScrollPhase, False))
        wait_until(lambda: len(render_events) > previous and viewer.is_ready)
        wheeled = limits()
        if any(not np.isclose(np.diff(wheeled[key])[0], .8 * np.diff(before[key])[0], rtol=1e-7)
               for key in before):
            raise RuntimeError("Wheel input did not zoom both physical heatmap dimensions")
        measure_alignment(config, "after_qt_wheel", wheeled)
        range_fields_match("after_qt_wheel")

        def drag(start_fraction, end_fraction):
            start, end = pointer(start_fraction), pointer(end_fraction)
            previous = len(render_events)
            mouse(QEvent.Type.MouseButtonPress, start, Qt.MouseButton.LeftButton, Qt.MouseButton.LeftButton)
            for fraction in (.25, .5, .75, 1):
                local = start[0] + (end[0] - start[0]) * fraction
                global_point = start[1] + (end[1] - start[1]) * fraction
                mouse(QEvent.Type.MouseMove, (local, global_point), Qt.MouseButton.NoButton, Qt.MouseButton.LeftButton)
                application.processEvents()
                time.sleep(.05)
            if len(render_events) != previous:
                raise RuntimeError("Native plot redrew and cleared axes during an active drag")
            mouse(QEvent.Type.MouseButtonRelease, end, Qt.MouseButton.LeftButton, Qt.MouseButton.NoButton)
            wait_until(lambda: len(render_events) > previous and viewer.is_ready)
            return limits()

        viewer.toolbar.pan()
        panned = drag((.5, .5), (.6, .6))
        viewer.toolbar.pan()
        if any(not np.isclose(np.diff(panned[key])[0], np.diff(wheeled[key])[0], rtol=1e-7)
               or np.allclose(panned[key], wheeled[key]) for key in before):
            raise RuntimeError("Native toolbar pan did not translate both axes while preserving their spans")
        measure_alignment(config, "after_qt_pan", panned)
        range_fields_match("after_qt_pan")
        viewer.toolbar.zoom()
        zoomed = drag((.25, .25), (.75, .75))
        viewer.toolbar.zoom()
        if any(np.diff(zoomed[key])[0] >= .6 * np.diff(panned[key])[0] for key in before):
            raise RuntimeError("Native toolbar rectangle zoom did not narrow both axes")
        measure_alignment(config, "after_qt_rectangle_zoom", zoomed)
        range_fields_match("after_qt_rectangle_zoom")
        application.processEvents()
        window.grab().save(str(output / "native_ranges_after_zoom.png"))
        viewer.toolbar.back()
        wait_until(lambda: viewer.is_ready and all(np.allclose(limits()[key], panned[key]) for key in before))
        range_fields_match("after_back")
        viewer.toolbar.forward()
        wait_until(lambda: viewer.is_ready and all(np.allclose(limits()[key], zoomed[key]) for key in before))
        range_fields_match("after_forward")
        viewer.toolbar.home()
        full_ranges = {"mz_range": raw_mz_domain, "drift_range": [0, native_drift_count * pusher_us / 1000]}
        wait_until(lambda: viewer.is_ready and all(np.allclose(limits()[key], full_ranges[key]) for key in before))
        range_fields_match("after_home")
        previous = len(render_events)
        local, global_point = pointer((.5, .5))
        QApplication.sendEvent(viewer.canvas, QWheelEvent(local, global_point, QPoint(), QPoint(0, 120),
            Qt.MouseButton.NoButton, Qt.KeyboardModifier.NoModifier, Qt.ScrollPhase.NoScrollPhase, False))
        wait_until(lambda: len(render_events) > previous and viewer.is_ready)
        if any(not np.isclose(np.diff(limits()[key])[0], .8 * np.diff(full_ranges[key])[0], rtol=1e-7)
               for key in before):
            raise RuntimeError("Reset acceptance did not start from a genuinely zoomed native viewport")
        panel.btn_viewer.click()
        wait_until(lambda: viewer.is_ready and all(np.allclose(limits()[key], full_ranges[key]) for key in before))
        range_fields_match("after_reset_action")
        report["native_qt_interaction"] = {"hover": hover, "wheel": "passed", "pan": "passed",
                                           "rectangle_zoom": "passed", "back_forward_home": "passed",
                                           "reset_action": "passed",
                                           "redraw_deferred_until_drag_release": True}
        changed = deepcopy(config)
        changed["defaults"].update({"mz_bins": 250, "noise_2d": 7.25, "ms_mode": "profile"})
        changed["initial_view"] = {"mz_range": [1100, 1800], "drift_range": [5, 10]}
        previous = len(render_events)
        panel.viewer_panel.spin_mzbins.setValue(250)
        panel.viewer_panel.spin_noise2d.setValue(7.25)
        panel.viewer_panel.combo_ms_mode.setCurrentIndex(panel.viewer_panel.combo_ms_mode.findData("profile"))
        for control, value in ((panel.viewer_panel.spin_mz_min, 1100), (panel.viewer_panel.spin_mz_max, 1800),
                               (panel.viewer_panel.spin_drift_min, 5), (panel.viewer_panel.spin_drift_max, 10)):
            control.setValue(value)
        wait_until(lambda: len(render_events) > previous and viewer.is_ready and
                   all(np.allclose(limits()[key], changed["initial_view"][key]) for key in before))
        if any(viewer._config["defaults"].get(key) != value
               for key, value in (("mz_bins", 250), ("noise_2d", 7.25), ("ms_mode", "profile"))):
            raise RuntimeError("Manual GUI settings did not reach the native renderer")
        measure_alignment(changed, "after_live_config_rebin")
        range_fields_match("after_manual_numeric_edit")
        viewport = {"mz_range": [1200, 1600], "drift_range": [6, 9]}
        viewer.ax_heatmap.set_xlim(viewport["mz_range"])
        viewer.ax_heatmap.set_ylim(viewport["drift_range"])
        wait_until(lambda: viewer.is_ready and np.allclose(viewer.snapshot.mz_range, viewport["mz_range"])
                   and np.allclose(viewer.snapshot.drift_range, viewport["drift_range"]))
        measure_alignment(changed, "after_mouse_zoom", viewport)
        range_fields_match("after_mouse_zoom")
        changed["defaults"]["noise_2d"] = 8
        captures = output / "current view exports with spaces"
        captures.mkdir(exist_ok=True)
        source_hashes = {kind: hashlib.sha256(record[kind].read_bytes()).hexdigest() for kind in ("ms", "im")}
        previous = len(render_events)
        panel.viewer_panel.spin_noise2d.setValue(8)
        wait_until(lambda: viewer._config["defaults"].get("noise_2d") == 8 and not viewer.is_ready)
        pending_render = previous
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
        if len(render_events) <= pending_render or not viewer.is_ready or viewer._config["defaults"].get("noise_2d") != 8:
            raise RuntimeError("The queued export did not wait for the applied live noise control")
        measure_alignment(changed, "after_control_preserving_zoom", viewport)
        measure_native_profiles("after_control_preserving_zoom")
        range_fields_match("after_control_preserving_zoom")
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
                ("drift_csv", "drift_time_ms", snapshot.drift_profile_centres, snapshot.drift_intensity, snapshot.drift_range)):
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
                                           "source_guard": "passed", "completion": "export_finished",
                                           "pending_render_export": True}
        report["native_viewer"] = {"renderer": type(viewer.canvas).__module__, "ready": viewer.is_ready,
                                    "needs_html_to_render": False, "drift_unit": snapshot.drift_unit,
                                    "mouse_zoom_preserved_after_controls": True}
        # Select an exact acquisition-local target, independently summing the
        # original positive TXT points instead of projecting coarse image bins.
        target_points = positive[(positive[:, 0] >= viewport["mz_range"][0]) &
                                 (positive[:, 0] <= viewport["mz_range"][1])]
        peak_mz = float(target_points[np.argmax(target_points[:, 2]), 0])
        roi_bounds = (peak_mz - 1e-6, peak_mz + 1e-6)
        editor = viewer.targets_editor
        buttons = {button.text(): button for button in editor.findChildren(QPushButton)}
        buttons["Add"].click()
        if editor.table.rowCount() != 1:
            raise RuntimeError("The target sheet Add action did not create its editable row")
        for column, value in ((0, "Exact peak"), (1, repr(roi_bounds[0])), (2, repr(roi_bounds[1]))):
            editor.table.item(0, column).setText(value)
        wait_until(lambda: viewer.is_ready and viewer.active_profile_range == roi_bounds)
        selected = positive[(positive[:, 0] >= roi_bounds[0]) & (positive[:, 0] <= roi_bounds[1])]
        expected_profile = np.bincount(selected[:, 1].astype(int), weights=selected[:, 2],
                                       minlength=native_drift_count)
        if (expected_profile.max() <= 0 or not np.array_equal(viewer.snapshot.drift_intensity, expected_profile)
                or not np.array_equal(viewer.ax_drift.lines[0].get_xdata(), expected_profile)
                or not np.array_equal(viewer.mobility_plot.ax().lines[0].get_ydata(), expected_profile)):
            raise RuntimeError("Selected target drift profiles differ from the independent exact TXT-point sum")
        centres_native = (np.arange(native_drift_count) + .5) * pusher_us / 1000
        support = centres_native[expected_profile > 0]
        pad = .2 * max(abs(support[0]), abs(support[-1]))
        focused_ranges = {"mz_range": list(roi_bounds), "drift_range": [max(0, support[0] - pad),
            min(native_drift_count * pusher_us / 1000, support[-1] + pad)]}
        if not all(np.allclose(limits()[key], focused_ranges[key], rtol=0, atol=1e-9) for key in focused_ranges):
            raise RuntimeError("Target selection did not focus exact mass bounds and padded native support")
        context = viewer.snapshot.baked_mz_range
        context_points = positive[(positive[:, 0] >= context[0]) & (positive[:, 0] <= context[1])]
        if not np.isclose(viewer.snapshot.raw2d.sum(), context_points[:, 2].sum(), rtol=0, atol=1e-9):
            raise RuntimeError("The target heatmap smoothing context changed the acquired RAW intensity")
        rows, columns = viewer.snapshot.raw2d.shape
        drift_edges = (np.arange(rows + 1) * native_drift_count) // rows
        mz_edges = np.linspace(context[0], context[1], columns + 1)
        independent_grid, _, _ = np.histogram2d(context_points[:, 1], context_points[:, 0],
            bins=(drift_edges, mz_edges), weights=context_points[:, 2])
        if not np.array_equal(independent_grid, viewer.snapshot.raw2d):
            raise RuntimeError("The target heatmap differs cellwise from the independent acquired-point histogram")
        range_fields_match("after_exact_target_selection")
        focused_config = deepcopy(changed)
        focused_config["initial_view"] = focused_ranges
        measure_alignment(focused_config, "after_exact_target_selection", smoothing_context=context)
        application.processEvents()
        window.grab().save(str(output / "native_target_profile.png"))
        for kind in ("drift_csv", "drift_png"):
            destination = captures / ("Exact peak target." + kind.rsplit("_", 1)[1])
            previous = len(export_events)
            viewer.export_current(kind, destination)
            wait_until(lambda: any(Path(event["path"]).resolve() == destination.resolve()
                       and not event.get("error") for event in export_events[previous:]))
        with (captures / "Exact peak target.csv").open(encoding="utf-8", newline="") as handle:
            roi_rows = list(csv.DictReader(handle))
        visible = (viewer.snapshot.drift_profile_centres >= focused_ranges["drift_range"][0]) & \
                  (viewer.snapshot.drift_profile_centres <= focused_ranges["drift_range"][1])
        if [(float(row["drift_time_ms"]), float(row["display_value"])) for row in roi_rows] != list(
                zip(viewer.snapshot.drift_profile_centres[visible], expected_profile[visible])):
            raise RuntimeError("Selected target mobiligram CSV differs from its exact displayed profile")
        if QImage(str(captures / "Exact peak target.png")).isNull():
            raise RuntimeError("Selected target mobiligram PNG did not render")
        editor.table.item(0, 3).setCheckState(Qt.Unchecked)
        editor.table.item(0, 4).setText("5.123456789")
        editor.table.item(0, 5).setText("10.987654321")
        wait_until(lambda: viewer.is_ready and np.allclose(viewer.snapshot.drift_range,
                   [5.123456789, 10.987654321], rtol=0, atol=1e-9))
        if not np.array_equal(viewer.snapshot.drift_intensity, expected_profile):
            raise RuntimeError("Manual target drift bounds changed the native acquired profile")
        range_fields_match("manual_target_drift")
        buttons["Delete"].click()
        wait_until(lambda: viewer.is_ready and viewer.active_profile_range is None)
        if editor.table.rowCount():
            raise RuntimeError("Deleting the selected target did not remove its row")
        target_csv = captures / "target sheet with spaces.csv"
        pd.DataFrame([{"name": "Imported exact peak", "mz_min": roi_bounds[0], "mz_max": roi_bounds[1]}]).to_csv(
            target_csv, index=False)
        editor.import_csv(target_csv)
        editor.table.selectRow(0)
        wait_until(lambda: viewer.is_ready and viewer.active_profile_range == roi_bounds)
        if not np.array_equal(viewer.snapshot.drift_intensity, expected_profile):
            raise RuntimeError("Imported target bounds did not preserve the exact selected drift profile")
        buttons["Full view"].click()
        wait_until(lambda: viewer.is_ready and viewer.active_profile_range is None and
                   all(np.allclose(limits()[key], full_ranges[key]) for key in full_ranges))
        measure_native_profiles("full_view_after_target")
        if panel._viewer_config().get("initial_view") != {"mz_range": None, "drift_range": None}:
            raise RuntimeError("Full view did not clear both configured ranges")
        report["target_profile"] = {"mz_bounds": list(roi_bounds), "raw_point_count": len(selected),
            "raw_intensity_sum": float(expected_profile.sum()), "independent_txt_sum": True,
            "exact_target_focus": focused_ranges, "automatic_and_manual_drift": "passed",
            "add_edit_delete_import_full_view": "passed", "mobiligram_png_csv": "passed",
            "mass_and_charge_required_for_visual_roi": False}
        # The explicit workspace action generates HTML from TXT on its worker;
        # native plotting above has never loaded or interpreted an HTML file.
        standalone = output / "standalone_copy"
        standalone.mkdir(exist_ok=True)
        generated = standalone / "exported offline viewer.html"
        returned = panel.export_portable_html(generated)
        if returned is None:
            raise RuntimeError("Explicit portable HTML export could not start")
        html_success, html_errors = [], []
        if panel._worker is None:
            raise RuntimeError("Portable HTML export did not create its worker")
        panel._worker.finished.connect(html_success.append)
        panel._worker.error.connect(html_errors.append)
        wait_until(lambda: not panel._busy, 120)
        if not html_success or html_errors or not generated.is_file():
            raise RuntimeError(f"Portable HTML worker did not complete successfully: {html_errors}")
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
        # The primary pipeline built the real 25-reference curve once. Exercise
        # the desktop review/save/reload/apply route using its existing ADP pair.
        protein_index = next(index for index in range(panel.combo_demo.count())
            if (panel.combo_demo.itemData(index) or {}).get("id") == "protein_adp")
        panel.combo_demo.setCurrentIndex(protein_index)
        if not panel._load_demo():
            raise RuntimeError("The single-acquisition ADP example could not load")
        panel.pick_out.set_path(str(output / "real_demo"))
        panel.reload_runs()
        if len(panel._rows) != 5:
            raise RuntimeError("The protein/ADP family must list all five complete RAW acquisitions")
        adp_index = next(index for index, row in enumerate(panel._rows) if row["run_name"] == "adp_edda")
        panel.run_table.selectRow(adp_index)
        panel._selected_run_changed()
        adp_record = dict(panel._rows[adp_index])
        wait_until(lambda: viewer.is_ready and viewer._record.get("run_name") == adp_record["run_name"]
                   and len(viewer.get_targets()) == 17)
        adp_im = np.loadtxt(adp_record["im"], ndmin=2)
        adp_ms = np.loadtxt(adp_record["ms"], ndmin=2)
        adp_period = find_pusher_period(adp_record["raw"])
        measure_native_profiles("adp_overview", adp_im, adp_ms, adp_period)
        if len(viewer.snapshot.drift_profile_centres) != 200:
            raise RuntimeError("The ADP example must retain all 200 acquired drift bins")
        viewer.targets_editor.table.selectRow(0)
        wait_until(lambda: viewer.is_ready and viewer.active_profile_range is not None)
        range_fields_match("saved_model_acquisition_loaded")
        workbench = panel.ccs_workbench
        def side_geometry():
            return {name: (widget.isVisible(), widget.mapToGlobal(QPoint()).x(),
                           widget.mapToGlobal(QPoint()).y(), widget.width(), widget.height())
                    for name, widget in (("data", panel.browser_panel), ("display", panel.inspector_scroll))}
        fixed_sides = side_geometry()
        panel.tabs.setCurrentIndex(1)
        application.processEvents()
        if side_geometry() != fixed_sides or not all(values[0] for values in fixed_sides.values()):
            raise RuntimeError("Calibration must keep the Data and Display panels visible and stationary")
        table = workbench.references.table
        columns = {table.horizontalHeaderItem(index).text(): index for index in range(table.columnCount())}
        for field in ("peak_select", "CCS_reported"):
            column = columns.get(field)
            if column is not None:
                table.scrollToItem(table.item(0, column))
                application.processEvents()
            if (column is None or table.isColumnHidden(column)
                    or table.horizontalHeader().sectionViewportPosition(column) < 0
                    or table.horizontalHeader().sectionViewportPosition(column) + table.columnWidth(column) > table.viewport().width()):
                raise RuntimeError(f"The calibration {field} control/header is not visible")
        peak_cell = table.item(0, columns["peak_select"])
        table.scrollToItem(peak_cell)
        application.processEvents()
        table.openPersistentEditor(peak_cell)
        application.processEvents()
        peak_editors = [editor for editor in table.findChildren(QComboBox) if editor.isVisible()
                        and {editor.itemData(index) for index in range(editor.count())} == {"highest", "last", "manual"}]
        if not peak_editors:
            raise RuntimeError("The calibration sheet has no visible Highest/Last peak selector")
        window.grab().save(str(output / "native_calibration_references.png"))
        table.closePersistentEditor(peak_cell)
        report["fixed_calibration_side_panels"] = {"geometry": fixed_sides,
            "highest_last_visible": True, "ccs_reported_visible": True}
        workbench.set_output_dir(output / "saved model reuse with spaces")
        workbench.load_calibration(ccs_dir / "calibration_curve.json")
        if workbench.get_calibration() is None or workbench.summary.table.rowCount() != 25:
            raise RuntimeError("The workbench could not load and review all 25 real calibration references")
        coefficients = {key: workbench.get_calibration()[key] for key in ("A", "X", "edc_coeff", "gas_mw", "method")}
        if not workbench.btn_save.isEnabled():
            raise RuntimeError("A fitted calibration must be saveable before generating target CCS")
        model_path = output / "saved calibration with spaces.json"
        workbench.save_calibration(model_path)
        if not model_path.is_file():
            raise RuntimeError("The fitted calibration was not saved")
        if json.loads(model_path.read_text(encoding="utf-8")).get("accepted"):
            raise RuntimeError("Saving before Generate CCS must not mark the fit accepted")
        workbench.clear_references(preserve_calibration=True)
        workbench.load_calibration(model_path)
        calibration = workbench.get_calibration()
        if (workbench._generated or workbench.references.table.rowCount()
                or any(calibration[key] != value for key, value in coefficients.items())
                or not calibration.get("review_warnings")):
            raise RuntimeError("Reloaded coefficients, review state, or poor-fit warnings were lost")
        panel.tabs.setCurrentIndex(1)
        application.processEvents()
        workbench.plot.canvas.draw()
        curve, placeholder = workbench.plot.ax(0, 0), workbench.plot.ax(0, 1)
        ratio = workbench.plot.canvas.device_pixel_ratio
        curve_height = float(curve.bbox.height / ratio)
        if curve_height < 120:
            raise RuntimeError("Calibration Review has less than 120 pixels of usable curve data height")
        from matplotlib.colors import to_rgba
        background = placeholder.get_facecolor()[:3]
        readable_text = any(text.get_text().strip() and max(abs(a - b) for a, b in
                            zip(to_rgba(text.get_color())[:3], background)) >= .3
                            for text in (placeholder.title, *placeholder.texts))
        if placeholder.axison or not readable_text:
            raise RuntimeError("Saved-reference placeholder must be readable and have no dummy axes")
        report["calibration_review_layout"] = {"curve_data_height": curve_height,
            "canvas_height": workbench.plot.canvas.height(), "saved_placeholder_axes_off": True,
            "saved_placeholder_readable": True}
        window.grab().save(str(output / "native_calibration_review.png"))
        before_apply = viewer.snapshot
        before_apply_roi = viewer.active_profile_range
        adp_source_hashes = {key: hashlib.sha256(Path(adp_record[key]).read_bytes()).hexdigest()
                            for key in ("ms", "im")}
        applied_events = []
        workbench.targets_applied.connect(applied_events.append)
        from deconvovo import imms_ccs_calibrate as calibration_backend
        build_curve = calibration_backend.build_calibration_curve

        def reject_refit(*args, **kwargs):
            raise RuntimeError("Applying the saved curve attempted to refit the standards")

        calibration_backend.build_calibration_curve = reject_refit
        try:
            workbench.btn_apply.click()
            if not workbench.is_busy:
                raise RuntimeError("Applying the saved curve to current ADP targets did not start")
            wait_until(lambda: not workbench.is_busy, 120)
        finally:
            calibration_backend.build_calibration_curve = build_curve
        if len(applied_events) != 1 or len(applied_events[0].get("analyte_summary", [])) != 17:
            raise RuntimeError("Saved-curve application did not publish all 17 current ADP target results")
        applied = applied_events[0]
        profiles = applied.get("native_profiles", [])
        if len(profiles) != 17 or viewer.ccs_plot.parentWidget() is not panel.ccs_profile_group:
            raise RuntimeError("Generate CCS did not deliver 17 native profiles to the right-column panel")
        if not panel.ccs_profile_group.isVisible() or panel.ccs_profile_group.y() <= panel.range_group.y():
            raise RuntimeError("The native CCS panel is not visible below Axis ranges")
        from deconvovo.ccs_convert import convert_profile_to_ccs
        from deconvovo.ccs_peak_pick import extract_drift_profile
        from deconvovo.io import read_im_txt
        acquired_im = read_im_txt(Path(adp_record["im"]))
        for profile in profiles:
            target = profile["target"]
            raw_profile = extract_drift_profile(acquired_im, target["mz"], target["mz_window"])
            _bins, drift, ccs, intensity = convert_profile_to_ccs(raw_profile,
                target["pusher_period_us"], target["mz"], target["z"], target["mw"],
                workbench.get_calibration(), target["method"])
            if not all(np.array_equal(profile[key], expected) for key, expected in
                       (("drift_time_ms", drift), ("ccs", ccs), ("intensity", intensity))):
                raise RuntimeError("A native target CCS plot changed the acquired profile or CCS transformation")
        workbench.save_calibration(model_path)
        if not json.loads(model_path.read_text(encoding="utf-8")).get("accepted"):
            raise RuntimeError("Generate CCS acceptance was not preserved when saving")
        workbench.load_calibration(model_path)
        if workbench._generated or viewer._ccs_profiles:
            raise RuntimeError("Loading a saved curve must require an explicit Generate CCS step")
        calibration_backend.build_calibration_curve = reject_refit
        try:
            workbench.btn_apply.click()
            wait_until(lambda: not workbench.is_busy, 120)
        finally:
            calibration_backend.build_calibration_curve = build_curve
        if len(applied_events) != 2 or len(viewer._ccs_profiles) != 17:
            raise RuntimeError("Reloaded saved calibration did not regenerate native CCS profiles")
        records = applied["analyte_summary"]
        baseline = {row["species"]: row for row in analytes if row["run"] == adp_record["run_name"]}
        for row in records:
            if (row["status"] != "complete" or not math.isfinite(float(row["peak_ccs"]))
                    or float(row["peak_ccs"]) <= 0 or not row["calibration_review_required"]
                    or not math.isclose(float(row["peak_ccs"]), float(baseline[row["species"]]["peak_ccs"]),
                                        rel_tol=0, abs_tol=.000050001)):
                raise RuntimeError("Saved-curve CCS results changed or lost the poor-fit review flag")
        conditions = applied.get("condition_checks", [])
        if not conditions:
            raise RuntimeError("Saved-curve acquisition-condition checks or warnings were lost")
        for check in conditions:
            required_warning = "conditions differ" if check["status"] == "different" else "cannot be verified"
            if check["status"] != "match" and not any(
                    required_warning in warning.casefold() for warning in applied.get("warnings", [])):
                raise RuntimeError("Saved-curve application did not retain its condition-specific warning")
        if any(row.get("extrapolated") for row in records) and not any(
                "extrapolat" in warning.casefold() for warning in applied.get("warnings", [])):
            raise RuntimeError("Saved-curve extrapolation was not reported")
        result_by_name = {row["species"]: row for row in records}
        for index, target in viewer.get_targets().iterrows():
            row = result_by_name[target["species"]]
            cells = [viewer.targets_editor.table.item(index, column) for column in (8, 9, 10)]
            if (not math.isclose(float(cells[0].text()), float(row["peak_ccs"]), rel_tol=0, abs_tol=.0050001)
                    or not math.isclose(float(cells[1].text()), float(row["peak_drift_time_ms"]), rel_tol=0, abs_tol=.0005001)
                    or cells[2].text() != "complete" or any(cell.flags() & Qt.ItemIsEditable for cell in cells)
                    or "Calibration review required: Yes" not in cells[0].toolTip()):
                raise RuntimeError("The current target sheet did not retain its read-only CCS results and review flags")
        if (viewer.active_profile_range != before_apply_roi or viewer.snapshot.mz_range != before_apply.mz_range
                or viewer.snapshot.drift_range != before_apply.drift_range):
            raise RuntimeError("Saved-curve application changed the selected ROI or viewing ranges")
        for name in ("display2d", "ms_mz", "ms_intensity", "drift_intensity"):
            if not np.array_equal(getattr(viewer.snapshot, name), getattr(before_apply, name), equal_nan=True):
                raise RuntimeError("Saved-curve application changed the native spectral data")
        if any(hashlib.sha256(Path(adp_record[key]).read_bytes()).hexdigest() != value
               for key, value in adp_source_hashes.items()):
            raise RuntimeError("Saved-curve application changed a converted source table")
        targets_saved = output / "saved target sheet with results.csv"
        viewer.targets_editor.export_csv(targets_saved)
        if not targets_saved.is_file():
            raise RuntimeError("The current target sheet and results were not saved")
        saved_frame = pd.read_csv(targets_saved)
        if len(saved_frame) != 17 or not {"drift_mode", "drift_min", "drift_max", "peak_ccs",
                "calibration_review_required", "conditions_status"} <= set(saved_frame.columns):
            raise RuntimeError("Target CSV omitted entered focus settings or calibration results/flags")
        viewer.targets_editor.import_csv(targets_saved)
        wait_until(lambda: viewer.is_ready)
        if any(viewer.targets_editor.table.item(index, 10).text() != "complete"
               for index in range(viewer.targets_editor.table.rowCount())):
            raise RuntimeError("Reimporting same-acquisition targets did not restore the read-only result status")
        viewer.targets_editor.table.selectRow(0)
        wait_until(lambda: viewer.is_ready and viewer.active_profile_range is not None)
        panel.tabs.setCurrentIndex(0)
        application.processEvents()
        if side_geometry() != fixed_sides:
            raise RuntimeError("Returning from Calibration moved the Data or Display panels")
        wait_until(useful_chart, 15)
        range_fields_match("after_saved_model_apply")
        header = viewer.targets_editor.table.horizontalHeader()
        if any(header.sectionViewportPosition(column) < 0 or
               header.sectionViewportPosition(column) + header.sectionSize(column) >
               viewer.targets_editor.table.viewport().width() for column in range(6)):
            raise RuntimeError("The six target focus input columns are not fully visible")
        footer_labels = [window.status_text]
        if panel.job_status.isVisible():
            footer_labels.append(panel.job_status)
        if any(label.height() > 50 or len(label.text()) > 200 for label in footer_labels):
            raise RuntimeError("Saved calibration warnings expanded the spectrum footer")
        window.grab().save(str(output / "native_saved_calibration_results.png"))
        calibration = workbench.get_calibration()
        report["saved_calibration_reuse"] = {"run": adp_record["run_name"], "references": 25, "targets": 17,
            "gui_load_save_generate_reload_generate": "passed", "native_ccs_profile_count": 17,
            "native_ccs_panel": "right below axis ranges", "refit_blocked": True,
            "coefficients_preserved": True, "reference_csv_not_required_for_reuse": True,
            "results_match_primary": True, "read_only_result_columns": True,
            "target_csv_inputs_results_flags_round_trip": True,
            "focus_input_columns_visible": 6, "concise_warning_footer": True,
            "roi_ranges_and_inputs_preserved": True, "review_warnings": calibration["review_warnings"],
            "extrapolated_targets": sum(bool(row.get("extrapolated")) for row in records),
            "condition_checks": conditions, "application_warnings": applied.get("warnings", []),
            "scientific_acceptance": "software_reuse_verified_poor_fit_not_validated_quantitation"}
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
