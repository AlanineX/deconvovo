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
        from PySide6.QtWidgets import QApplication, QFileDialog
        from PySide6.QtCore import QUrl
        from PySide6.QtGui import QFont, QFontDatabase, QFontInfo, QImage
        from PySide6.QtWebEngineCore import QWebEngineDownloadRequest
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
        results = process(production / profile["source"], output / "real_demo", config_path=config_path,
                          calibrant_csv=calibrant_path, analyte_csv=analyte_path)
        adp_html = output / "real_demo" / f"{adp_raw.stem}_2d_imms.html"
        if len(results["runs"]) != 5 or str(adp_html) not in results["viewers"]:
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
        results = process(Path(converted["converted_dir"]), output / "betaine_demo",
                          convert=False, config_path=config_path)
        html = output / "betaine_demo" / f"{raw.stem}_2d_imms.html"
        if (len(results["viewers"]) != 1 or str(html) not in results["viewers"] or results["ccs_dir"]
                or len(results["runs"]) != 1 or results["runs"][0]["status"] != "converted input"):
            raise RuntimeError("Betaine viewer-only must render its converted pair without CCS")
        report["convert_only"] = {"run": raw.stem, "ms_im_pairs": 1, "viewers": 0, "ccs": False}
        report["viewer_only"] = {"run": raw.stem, "viewers": 1, "reconverted": False, "ccs": False}
        report["preview_profile"] = profile["id"]
        report["preview_run"] = raw.stem
        contents = html.read_text(encoding="utf-8")
        if re.search(r'<(?:script|link|img)\b[^>]*(?:src|href)\s*=\s*[\"\']https?://', contents, re.I):
            raise RuntimeError("Exported HTML contains an external resource")
        if "plotly.js" not in contents or "Plotly.newPlot" not in contents:
            raise RuntimeError("Exported HTML is missing its embedded Plotly runtime")
        report["raw_reader"] = converted["runs"][0]
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
        panel.combo_demo.setCurrentIndex(next(index for index in range(panel.combo_demo.count())
            if (panel.combo_demo.itemData(index) or {}).get("id") == profile["id"]))
        panel._load_demo()
        panel.viewer_panel.set_config(config)
        # Exercise quoted CSV paths on the ADP tables before returning the
        # viewer-only demo to its own empty calibration/analyte tables.
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
        pasted_paths = [(panel.pick_source, raw),
                        (panel.pick_out, html.parent)]
        for picker, expected_path in pasted_paths:
            picker.set_path(f'"{expected_path}"')
            if Path(picker.path()).resolve() != expected_path.resolve():
                raise RuntimeError("The GUI did not preserve an Explorer-quoted absolute path")
        report["quoted_gui_paths"] = {"source_output_and_csv": "passed",
                                      "bundled_raw_path_contains_space": " " in str(raw),
                                      "output_path_contains_space": " " in str(output)}
        panel.reload_runs()
        panel.run_table.selectRow(next(index for index, row in enumerate(panel._rows)
                                      if row["run_name"] == raw.stem))
        viewer = panel.viewer_widget
        if viewer.web_view is None:
            raise RuntimeError("Frozen console did not create its embedded browser")
        record = {"run_name": raw.stem, "viewer": html,
                  "ms": Path(results["converted_dir"]) / f"{raw.stem}_ms.txt",
                  "im": Path(results["converted_dir"]) / f"{raw.stem}_im.txt"}
        # Derive cell coordinates independently from the acquired tables and
        # RAW metadata, rather than trusting two axes with the same bad range.
        import numpy as np
        from deconvovo.io import find_pusher_period
        acquired = np.loadtxt(record["im"])
        positive = acquired[acquired[:, 2] > 0]
        native_drift_count = int(positive[:, 1].max()) + 1
        raw_mz_domain = [float(positive[:, 0].min()), float(positive[:, 0].max())]
        baked_mz_domain = [max(raw_mz_domain[0], config["initial_view"]["mz_range"][0]),
                           min(raw_mz_domain[1], config["initial_view"]["mz_range"][1])]
        pusher_us = find_pusher_period(raw)
        if not pusher_us or pusher_us <= 0:
            raise RuntimeError("Betaine alignment requires its actual RAW pusher period")

        def measure_alignment(settings, label):
            # A control can acknowledge its new value before Plotly finishes
            # repainting the traces; measure the committed render.
            javascript("window.watersSmokeSettled=false;Promise.resolve(window.watersPendingRender)"
                       ".then(function(){window.watersSmokeSettled=true;});true")
            wait_until(lambda: javascript("window.watersSmokeSettled===true"), 15)
            render_error = javascript("window.watersConsoleConfigError||''")
            if render_error:
                raise RuntimeError(f"Viewer configuration failed: {render_error}")
            view = settings["initial_view"]
            domain = baked_mz_domain
            bins = settings["defaults"]["mz_bins"]
            if bins == "native":
                bins = min(25600, max(200, len(np.unique(positive[
                    (positive[:, 0] >= domain[0]) & (positive[:, 0] <= domain[1]), 0]))))
            drift_bins = settings["defaults"].get("drift_bins", "native")
            expected = {"mz_range": view["mz_range"], "drift_range": view["drift_range"],
                        "mz_domain": domain, "mz_bins": int(bins),
                        "drift_bins": native_drift_count if drift_bins == "native" else int(drift_bins),
                        "drift_upper_ms": native_drift_count * pusher_us / 1000,
                        "raw_pusher_us": pusher_us}
            source = "var expected=" + json.dumps(expected) + ";" + """
              var p=document.getElementById('plot'),lo=p._fullLayout;
              function compare(heat,profile,range,vertical){
                var h=lo[heat],a=lo[profile],samples=[];
                for(var i=0;i<=4;i++){
                  var t=i/4,v=range[0]+t*(range[1]-range[0]);
                  var oracle=h._offset+h._length*(vertical?1-t:t);
                  var hp=h._offset+h.l2p(v),ap=a._offset+a.l2p(v);
                  samples.push({coordinate:v,expected_pixel:oracle,heatmap_pixel:hp,
                    profile_pixel:ap,error_px:Math.max(Math.abs(hp-oracle),Math.abs(ap-oracle))});
                }
                return {heatmap_axis:heat,profile_axis:profile,heatmap_type:h.type,profile_type:a.type,
                  heatmap_range:h.range,profile_range:a.range,samples:samples,
                  max_error_px:Math.max.apply(null,samples.map(function(s){return s.error_px}))};
              }
              function centres(values,count,lower,upper){
                if(values.length!==count)return {count:values.length,expected_count:count,max_error:null};
                var error=0;
                for(var i=0;i<count;i++)error=Math.max(error,Math.abs(values[i]-(lower+(i+0.5)*(upper-lower)/count)));
                return {count:values.length,expected_count:count,max_error:error};
              }
              return {expected_from_raw_and_preset:expected,
                mz:compare('xaxis2','xaxis4',expected.mz_range,false),
                drift:compare('yaxis2','yaxis',expected.drift_range,true),
                heatmap_mz_centres:centres(p.data[2].x,expected.mz_bins,expected.mz_domain[0],expected.mz_domain[1]),
                heatmap_drift_centres:centres(p.data[2].y,expected.drift_bins,0,expected.drift_upper_ms),
                profile_drift_centres:centres(p.data[0].y,expected.drift_bins,0,expected.drift_upper_ms)};
            """
            value = javascript("JSON.stringify((function(){" + source + "})())")
            measured = json.loads(value) if value else {}
            report.setdefault("profile_alignment", {})[label] = measured
            for name in ("mz", "drift"):
                axis = measured.get(name, {})
                if (axis.get("heatmap_type") != "linear" or axis.get("profile_type") != "linear"
                        or not math.isfinite(axis.get("max_error_px", math.inf))
                        or axis["max_error_px"] > 1):
                    raise RuntimeError(f"{label} {name} profile differs from the independent physical pixel projection")
            for name in ("heatmap_mz_centres", "heatmap_drift_centres", "profile_drift_centres"):
                cells = measured.get(name, {})
                if (cells.get("count") != cells.get("expected_count")
                        or cells.get("max_error") is None or cells["max_error"] > 1e-7):
                    raise RuntimeError(f"{label} {name} differ from RAW-derived histogram cell centres")

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
        measure_alignment(config, "initial_detail")
        # JavaScript acknowledges layout before Chromium submits its next
        # frame. Let the compositor repaint before grabbing a validation image.
        repaint_until = time.monotonic() + 0.75
        while time.monotonic() < repaint_until:
            application.processEvents()
            time.sleep(0.02)
        window.grab().save(str(output / "native_console.png"))
        changed = deepcopy(config)
        changed["defaults"].update({"mz_bins": 250, "noise_2d": 7.25, "ms_mode": "stick"})
        changed["initial_view"] = {"mz_range": [1100, 1800], "drift_range": [5, 10]}
        viewer.apply_config(changed)

        def controls_applied():
            ready()
            return (state.get("mz_bins") == "250" and state.get("noise") == "7.25"
                    and state.get("mode") == "stick")

        wait_until(controls_applied)
        measure_alignment(changed, "after_live_rebin_and_zoom")
        captures = output / "current view exports with spaces"
        captures.mkdir(exist_ok=True)
        downloads = []

        def track_download(download):
            event = {"name": download.downloadFileName(), "state": download.state().name}
            downloads.append(event)
            download.stateChanged.connect(lambda state: event.update(state=state.name))

        viewer._profile.downloadRequested.connect(track_download)
        source_hashes = {kind: hashlib.sha256(record[kind].read_bytes()).hexdigest() for kind in ("ms", "im")}
        original_dialog = QFileDialog.getSaveFileName
        try:
            QFileDialog.getSaveFileName = staticmethod(lambda parent, title, suggested, filters:
                                                       (str(captures / Path(suggested).name), filters))
            exported = {}
            for kind in ("full_png", "heatmap_png", "mz_png", "drift_png", "heatmap_csv", "mz_csv", "drift_csv"):
                previous = len(downloads)
                viewer.export_current(kind)
                wait_until(lambda: len(downloads) > previous and downloads[-1]["state"] in
                           {"DownloadCompleted", "DownloadCancelled", "DownloadInterrupted"}, 30)
                if downloads[-1]["state"] != "DownloadCompleted":
                    raise RuntimeError(f"Current-view {kind} download did not complete: {downloads[-1]}")
                destination = captures / downloads[-1]["name"]
                if not destination.is_file() or destination.stat().st_size == 0:
                    raise RuntimeError(f"Completed current-view {kind} download has no bytes")
                exported[kind] = destination
            displayed = json.loads(javascript("""JSON.stringify({heatmap:P.data[2],drift:P.data[0],
                mz:P.data[1],visible:getVisibleRange(),transform:document.getElementById('zscale').value})"""))
            visible = displayed["visible"]
            hm = displayed["heatmap"]
            with exported["heatmap_csv"].open(encoding="utf-8", newline="") as handle:
                reader = csv.DictReader(handle)
                if reader.fieldnames != ["mz", "drift_time_ms", "display_value", "intensity_transform"]:
                    raise RuntimeError("Current heatmap CSV must label its coordinates and displayed transform")
                rows = list(reader)
            actual = [(float(row["mz"]), float(row["drift_time_ms"]), float(row["display_value"])) for row in rows]
            expected = [(hm["x"][m], hm["y"][d], hm["z"][d][m])
                        for d in range(visible["dBl"], visible["dBh"] + 1)
                        for m in range(visible["mBl"], visible["mBh"] + 1)]
            if actual != expected or not all(row["intensity_transform"] == displayed["transform"] for row in rows):
                raise RuntimeError("Current heatmap CSV does not contain the displayed coordinates and values")
            csv_counts = {"heatmap": len(rows)}
            for kind, coordinate, trace, lower, upper in (
                    ("mz_csv", "mz", displayed["mz"], visible["mzLo"], visible["mzHi"]),
                    ("drift_csv", "drift_time_ms", displayed["drift"], visible["drLo"], visible["drHi"])):
                with exported[kind].open(encoding="utf-8", newline="") as handle:
                    reader = csv.DictReader(handle)
                    if reader.fieldnames != [coordinate, "display_value"]:
                        raise RuntimeError(f"Current {kind} has incorrect column labels")
                    rows = list(reader)
                actual = [(float(row[coordinate]), float(row["display_value"]))
                          if row[coordinate] and row["display_value"] else (None, None) for row in rows]
                coordinates, values = (trace["x"], trace["y"]) if kind == "mz_csv" else (trace["y"], trace["x"])
                expected = [(value, intensity) if value is not None and intensity is not None else (None, None)
                            for value, intensity in zip(coordinates, values)
                            if value is None or intensity is None or lower <= value <= upper]
                if actual != expected:
                    raise RuntimeError(f"Current {kind} differs from its displayed trace")
                csv_counts[kind.removesuffix("_csv")] = len(rows)
            images = {}
            for kind in ("full_png", "heatmap_png", "mz_png", "drift_png"):
                image = QImage(str(exported[kind]))
                if image.isNull() or image.width() < 1000 or image.height() < 700:
                    raise RuntimeError(f"Current {kind} is not a valid rendered PNG")
                images[kind] = {"width": image.width(), "height": image.height()}
            # Exercise the actual native download guard, then independently
            # confirm RAW containment and every original MS/IM byte.
            QFileDialog.getSaveFileName = staticmethod(lambda parent, title, suggested, filters:
                                                       (str(record["ms"]), filters))
            previous = len(downloads)
            viewer.export_current("heatmap_csv")
            wait_until(lambda: len(downloads) > previous and downloads[-1]["state"] == "DownloadCancelled", 15)
            try:
                viewer._output_path(raw / "protected current view.csv", {".csv"})
            except ValueError:
                pass
            else:
                raise RuntimeError("Current-view export did not protect the original RAW directory")
            if any(hashlib.sha256(record[kind].read_bytes()).hexdigest() != digest
                   for kind, digest in source_hashes.items()):
                raise RuntimeError("Current-view exports changed an original MS or IM table")
            report["current_view_exports"] = {"formats": sorted(exported), "completed": 7,
                                               "csv_rows": csv_counts, "csv_matches_display": True,
                                               "png_dimensions": images, "source_guard": "passed",
                                               "downloads": downloads}
        finally:
            QFileDialog.getSaveFileName = staticmethod(original_dialog)
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
        measure_alignment(config, "moved_standalone")
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
