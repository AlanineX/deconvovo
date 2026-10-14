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
                       web_view=viewer.web_view, viewer_status=viewer.status)
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
        os.environ["DECONVOVO_ENABLE_HEADLESS_WEBENGINE"] = "1"
        os.environ.pop("DECONVOVO_DISABLE_WEBENGINE", None)
        if os.name == "nt":
            # Windows' offscreen plugin has no system font database. Validate
            # the same native platform and fonts used by the desktop app.
            os.environ["QT_QPA_PLATFORM"] = "windows"
        from gui.widgets.worker import set_gui_mode
        set_gui_mode()
        from PySide6.QtWidgets import QApplication
        from PySide6.QtCore import QUrl
        from PySide6.QtGui import QFont, QFontDatabase, QFontInfo
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
        small_molecule = production / "small_molecule"
        raw = small_molecule / "raw" / "20260320_BETAINE_V2.raw"
        config_path = small_molecule / "viewer_config.json"
        config = json.loads(config_path.read_text(encoding="utf-8"))
        results = process(small_molecule / "raw", output / "real_demo", config_path=config_path,
                          calibrant_csv=small_molecule / "calibrants.csv",
                          analyte_csv=small_molecule / "analytes.csv")
        html = output / "real_demo" / f"{raw.stem}_2d_imms.html"
        if len(results["runs"]) != 2 or str(html) not in results["viewers"]:
            raise RuntimeError("Small-molecule smoke requires tune mix and Betaine acquisition outputs")
        ccs_dir = Path(results["ccs_dir"])
        with (ccs_dir / "calibrant_summary.csv").open(encoding="utf-8", newline="") as handle:
            references = list(csv.DictReader(handle))
        with (ccs_dir / "analyte_summary.csv").open(encoding="utf-8", newline="") as handle:
            analytes = list(csv.DictReader(handle))
        if len(references) != 6 or results["calibration"]["n_points"] != 6 or len(analytes) != 1:
            raise RuntimeError("CCS smoke requires six valid references and one Betaine analyte")
        peak_ccs = float(analytes[0]["peak_ccs"])
        if not math.isfinite(peak_ccs) or peak_ccs <= 0:
            raise RuntimeError("Betaine CCS must be finite and positive")
        reference_times = [float(row["drift_time_ms"]) for row in references]
        peak_time = float(analytes[0]["peak_drift_time_ms"])
        pngs = [ccs_dir / "calibration_plots.png", *sorted((ccs_dir / "cali_png").glob("*.png"))]
        pngs += [ccs_dir / folder / f"{raw.stem}_Betaine_ccs.png" for folder in
                 ("png_raw", "png_smoothed", "png_overlay")]
        if len(pngs) != 10 or any(not p.is_file() or not p.read_bytes().startswith(b"\x89PNG\r\n\x1a\n") for p in pngs):
            raise RuntimeError("CCS smoke did not publish all ten calibration/analyte PNGs")
        report["ccs_smoke"] = {"method": "direct", "raw_runs": 2, "calibrants": 6, "analytes": 1,
                               "peak_ccs_angstrom_squared": peak_ccs, "peak_drift_time_ms": peak_time,
                               "calibrant_drift_range_ms": [min(reference_times), max(reference_times)],
                               "extrapolated": peak_time < min(reference_times) or peak_time > max(reference_times),
                               "scientific_acceptance": "software_demo_requires_review_not_independently_validated",
                               "log_fit_r2": results["calibration"]["r2_lnln"],
                               "pngs": [p.relative_to(ccs_dir).as_posix() for p in pngs]}
        contents = html.read_text(encoding="utf-8")
        if re.search(r'<(?:script|link|img)\b[^>]*(?:src|href)\s*=\s*[\"\']https?://', contents, re.I):
            raise RuntimeError("Exported HTML contains an external resource")
        if "plotly.js" not in contents or "Plotly.newPlot" not in contents:
            raise RuntimeError("Exported HTML is missing its embedded Plotly runtime")
        report["raw_reader"] = next(run for run in results["runs"] if run["run_name"] == raw.stem)
        report["portable_html_bytes"] = html.stat().st_size

        def wait_until(predicate, seconds=90):
            deadline = time.monotonic() + seconds
            while time.monotonic() < deadline:
                application.processEvents()
                if predicate():
                    return
                time.sleep(0.02)
            raise RuntimeError("Frozen GUI did not reach the expected state")

        def javascript(source):
            reply = []
            viewer.web_view.page().runJavaScript(source, lambda value: reply.append(value))
            wait_until(lambda: bool(reply), 15)
            return reply[0]

        window = MainWindow()
        window.resize(1280, 800)
        window.show()
        application.processEvents()
        # Initial Windows placement can shrink a window to the CI desktop.
        # Validate the claimed client size after that placement has finished.
        window.resize(1280, 800)
        application.processEvents()
        if window.width() != 1280 or window.height() != 800:
            raise RuntimeError(f"Requested 1280x800 client, obtained {window.width()}x{window.height()}")
        actual_font = QFontInfo(application.font())
        report["native_gui_fonts"] = {"family_count": len(QFontDatabase.families()),
                                      "actual_family": actual_font.family(),
                                      "pixel_size": actual_font.pixelSize()}
        if (not report["native_gui_fonts"]["family_count"] or not actual_font.family().strip()
                or actual_font.pixelSize() <= 0):
            raise RuntimeError("Qt could not load usable desktop fonts; native GUI labels cannot render.")
        panel = window.get_panel(0)
        panel.viewer_panel.set_config(config)
        panel.pick_source.set_path(str(raw))
        panel.pick_out.set_path(str(html.parent))
        panel.reload_runs()
        viewer = panel.viewer_widget
        if viewer.web_view is None:
            raise RuntimeError("Frozen console did not create its embedded browser")
        record = {"run_name": raw.stem, "viewer": html,
                  "ms": Path(results["converted_dir"]) / f"{raw.stem}_ms.txt",
                  "im": Path(results["converted_dir"]) / f"{raw.stem}_im.txt"}
        viewer.set_run(record, config, raw)
        wait_until(lambda: viewer._html_ready)
        state_script = """JSON.stringify((function(){var p=document.getElementById('plot');return {
          version:typeof Plotly==='undefined'?null:Plotly.version,
          traces:p&&p.data?p.data.length:0,
          ready:!!(p&&p._fullLayout),
          controls:['nbins','dbins','msmode','noise2d','denoise','sm2d','cmap','zscale']
            .every(function(id){return !!document.getElementById(id)}),
          mz_bins:document.getElementById('nbins').value,
          noise:document.getElementById('noise2d').value,
          mode:document.getElementById('msmode').value
        };})())"""
        state = {}

        def ready():
            nonlocal state
            value = javascript(state_script)
            state = json.loads(value) if value else {}
            return state.get("ready") and state.get("traces", 0) >= 3 and state.get("controls")

        wait_until(ready)
        geometry_script = """JSON.stringify((function(){
          var p=document.getElementById('plot'),r=p.getBoundingClientRect();return {
            styled:!!document.getElementById('waters-console-style')&&p.dataset.consoleStyled==='1',
            width:r.width,height:r.height,
            fully_visible:r.top>=-2&&r.left>=-2&&r.bottom<=innerHeight+2&&r.right<=innerWidth+2
          };})())"""
        geometry = {}

        def useful_chart():
            nonlocal geometry
            value = javascript(geometry_script)
            geometry = json.loads(value) if value else {}
            report["console_chart_at_1280x800"] = geometry
            report["embedded_viewer_before_geometry_check"] = state
            return (geometry.get("styled") and geometry.get("fully_visible")
                    and geometry.get("height", 0) >= 250 and geometry.get("width", 0) >= 300)

        wait_until(useful_chart, 15)
        report["console_chart_at_1280x800"] = geometry
        # JavaScript acknowledges layout before Chromium submits its next
        # frame. Let the compositor repaint before grabbing a validation image.
        repaint_until = time.monotonic() + 0.75
        while time.monotonic() < repaint_until:
            application.processEvents()
            time.sleep(0.02)
        window.grab().save(str(output / "native_console.png"))
        changed = deepcopy(config)
        changed["defaults"].update({"mz_bins": 250, "noise_2d": 7.25, "ms_mode": "stick"})
        viewer.apply_config(changed)

        def controls_applied():
            ready()
            return (state.get("mz_bins") == "250" and state.get("noise") == "7.25"
                    and state.get("mode") == "stick")

        wait_until(controls_applied)
        report["embedded_viewer"] = state
        report["embedded_viewer"]["url_scheme"] = viewer.web_view.url().scheme()
        report["embedded_viewer"]["remote_requests_during_render"] = viewer.blocked_requests
        if viewer.blocked_requests:
            raise RuntimeError("The local viewer attempted to load an external resource")
        standalone = output / "standalone_copy"
        standalone.mkdir(exist_ok=True)
        copied = viewer.export_html(standalone / "renamed offline viewer.html")
        digest = hashlib.sha256(html.read_bytes()).hexdigest()
        if hashlib.sha256(copied.read_bytes()).hexdigest() != digest:
            raise RuntimeError("Standalone export differs from the generated viewer")
        viewer.set_run({"run_name": "Moved standalone HTML", "viewer": copied}, config)
        wait_until(lambda: viewer._html_ready)
        wait_until(ready)
        report["standalone_copy"] = {"sha256": digest, "bytes": copied.stat().st_size,
                                      "rendered_from_file_url": viewer.web_view.url().isLocalFile(),
                                      "traces": state["traces"]}
        # Direct navigation reaches the interceptor before any DNS or HTTP.
        viewer.web_view.load(QUrl("https://offline-test.invalid/blocked"))
        wait_until(lambda: viewer.blocked_requests > 0, 15)
        report["remote_navigation_blocked"] = True
        viewer.set_run({"run_name": "Moved standalone HTML", "viewer": copied}, config)
        wait_until(lambda: viewer._html_ready)
        wait_until(ready)
        report["status"] = "passed"
        return_code = 0
    except Exception:
        report["status"] = "failed"
        report["error"] = traceback.format_exc()
        return_code = 1
    finally:
        if window is not None:
            try:
                repaint_until = time.monotonic() + 0.75
                while time.monotonic() < repaint_until:
                    application.processEvents()
                    time.sleep(0.02)
                final_image = output / "native_console_final_state.png"
                captured = window.grab()
                report["final_screenshot"] = {"path": final_image.name,
                                              "saved": captured.save(str(final_image), "PNG")}
                if report["status"] != "passed":
                    captured.save(str(output / "native_console.png"), "PNG")
                report["qt_layout_diagnostics"] = _qt_diagnostics(application, window, panel, viewer)
                if viewer is not None and viewer.web_view is not None:
                    value = javascript("""JSON.stringify((function(){var p=document.getElementById('plot');return {
                      inner_width:innerWidth,inner_height:innerHeight,device_pixel_ratio:devicePixelRatio,
                      plotly_ready:!!(p&&p._fullLayout),traces:p&&p.data?p.data.length:0,
                      plot_rect:p?{width:p.getBoundingClientRect().width,height:p.getBoundingClientRect().height}:null,
                      layout_width:p&&p._fullLayout?p._fullLayout.width:null,
                      layout_height:p&&p._fullLayout?p._fullLayout.height:null
                    };})())""")
                    report["browser_final_state"] = json.loads(value) if value else None
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
