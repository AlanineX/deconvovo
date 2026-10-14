"""Offline enhanced Waters HTML with real, unsmoothed TXT preview fallback."""
from __future__ import annotations

from copy import deepcopy
import json
import math
import os
from pathlib import Path
import platform
import shutil
import threading

import numpy as np
from PySide6.QtCore import QObject, QRunnable, QThreadPool, QTimer, Qt, QUrl, Signal, Slot
from PySide6.QtGui import QDesktopServices
from PySide6.QtWidgets import QFileDialog, QLabel, QPushButton, QStackedWidget, QTabWidget, QVBoxLayout, QWidget

from deconvovo.paths import local_path
from gui.widgets.plot_canvas import PlotCanvas


def _webengine_allowed() -> bool:
    if os.environ.get("DECONVOVO_DISABLE_WEBENGINE") == "1":
        return False
    if os.environ.get("DECONVOVO_ENABLE_HEADLESS_WEBENGINE") == "1":
        return True
    return not (
        os.environ.get("QT_QPA_PLATFORM", "").casefold() in {"offscreen", "minimal"}
        or (platform.system() == "Linux" and not os.environ.get("DISPLAY")
            and not os.environ.get("WAYLAND_DISPLAY"))
    )


try:
    if not _webengine_allowed():
        raise ImportError("WebEngine disabled for this process")
    from PySide6.QtWebEngineCore import QWebEnginePage, QWebEngineProfile, QWebEngineUrlRequestInterceptor
    from PySide6.QtWebEngineWidgets import QWebEngineView
    WEBENGINE_AVAILABLE = True
except ImportError:
    WEBENGINE_AVAILABLE = False


def is_local_resource(url: QUrl) -> bool:
    """The desktop viewer must never need or contact a remote resource."""
    return url.scheme().casefold() in {"file", "data", "blob", "qrc", "about"}


if WEBENGINE_AVAILABLE:
    class _OfflineRequests(QWebEngineUrlRequestInterceptor):
        def __init__(self, parent=None):
            super().__init__(parent)
            self.blocked_requests = 0

        def interceptRequest(self, info):
            if not is_local_resource(info.requestUrl()):
                self.blocked_requests += 1
                info.block(True)


# Presentation changes only. The original document keeps its own data,
# scientific transforms, linked axes, controls, and export handlers.
_CONSOLE_STYLE = """
(function() {
  var style = document.getElementById('waters-console-style');
  if (!style) { style=document.createElement('style'); style.id='waters-console-style';
    document.head.appendChild(style); }
  style.textContent = `
    html,body{height:100%;width:100%;overflow:hidden;background:#111a24;color:#e5edf5}
    body{display:flex;flex-direction:column;font-family:'Segoe UI',Arial,sans-serif}
    #title,#subtitle,#toolbar,#exportpanel{display:none!important}
    #main{flex:1;display:flex!important;align-items:stretch!important;flex-direction:column;min-height:0!important;
      height:auto!important;overflow:hidden}
    #plot{flex:1;height:auto!important;width:100%!important;min-width:0;min-height:0!important}
  `;
  var p=document.getElementById('plot');
  function resize() {
    if (!p || !p._fullLayout || typeof Plotly==='undefined') return;
    var width=p.clientWidth, height=p.clientHeight;
    if(width<200 || height<100) return; // hidden tabs are resized when visible
    if (p._fullLayout.width===width && p._fullLayout.height===height && p.dataset.consoleStyled) return;
    var layout={width:width,height:height,autosize:true,paper_bgcolor:'#111a24',
      plot_bgcolor:'#111a24','font.color':'#e5edf5','font.size':11,
      'margin.l':55,'margin.r':40,'margin.t':12,'margin.b':48};
    ['xaxis','xaxis2','xaxis3','xaxis4','yaxis','yaxis2','yaxis3','yaxis4'].forEach(function(a){
      layout[a+'.gridcolor']='#283442';layout[a+'.linecolor']='#536577';
      layout[a+'.minor.gridcolor']='#283442';
      layout[a+'.minor.showgrid']=false;
      layout[a+'.tickfont.size']=10;layout[a+'.tickfont.color']='#8b9bad';
      layout[a+'.title.font.size']=12;layout[a+'.title.font.color']='#e5edf5';
    });
    if(!p.dataset.consoleStyled) {
      var profiles=[];p.data.forEach(function(t,i){if(t.type==='scatter')profiles.push(i)});
      Plotly.restyle(p,{'line.color':'#8dcfc4','line.width':1.1},profiles);
    }
    p.dataset.consoleStyled='1'; Plotly.relayout(p,layout);
  }
  if (!window.watersConsoleResize) {
    window.watersConsoleResize=resize;
    new ResizeObserver(resize).observe(p);
    window.addEventListener('resize',resize);
  }
  resize();
  return !!(p && p._fullLayout);
})()
"""


def _config_script(config: dict) -> str:
    """Set real HTML controls, then use the existing viewer's apply functions."""
    values = json.dumps(config.get("defaults", {}), ensure_ascii=True)
    presets = json.dumps(config.get("presets", {}).get("colormap", []), ensure_ascii=True)
    initial = config.get("initial_view") or {}
    for key in ("mz_range", "drift_range"):
        bounds = initial.get(key)
        if bounds is not None and (not isinstance(bounds, (list, tuple)) or len(bounds) != 2
                or not all(isinstance(value, (int, float)) and math.isfinite(value) for value in bounds)
                or bounds[0] >= bounds[1]):
            raise ValueError(f"The {key.replace('_', ' ')} must contain a finite minimum lower than its maximum.")
    ranges = json.dumps(initial, ensure_ascii=True)
    return """
    (function(){
      var d=DEFAULTS, presets=PRESETS, ranges=RANGES, changed=false, rebinned=false;
      var map={smooth_2d:'sm2d',scale:'mzscale',noise_2d:'noise2d',denoise:'denoise',
        smooth_drift:'smdrift',noise_drift:'noisedrift',ms_mode:'msmode',noise_mz:'noisemz',
        mz_bins:'nbins',drift_bins:'dbins',colormap:'cmap',render:'render',z_scale:'zscale'};
      if (typeof applyAll!=='function') return false;
      var rawLo=typeof RAW_MZ_LO0==='undefined'?HM_MZ_LO0:RAW_MZ_LO0;
      var rawHi=typeof RAW_MZ_HI0==='undefined'?HM_MZ_HI0:RAW_MZ_HI0;
      if(ranges.mz_range && ((ranges.mz_range[0]<HM_MZ_LO0 && rawLo<HM_MZ_LO0) ||
                            (ranges.mz_range[1]>HM_MZ_HI0 && rawHi>HM_MZ_HI0))) {
        window.watersConsoleConfigError='The current HTML covers m/z '+HM_MZ_LO0+' to '+HM_MZ_HI0+
          '. Use Process & view to regenerate the requested range.';
        return window.watersConsoleConfigError;
      }
      window.watersConsoleConfigError='';
      Object.keys(map).forEach(function(key){
        if (!(key in d)) return;
        var el=document.getElementById(map[key]); if(!el)return;
        var requested=d[key], value=null;
        Array.from(el.options).some(function(option){
          var match;
          if(key==='colormap') match=option.textContent.trim()===String(requested);
          else if(typeof requested==='object') {
            try {var obj=JSON.parse(option.value);
              match=Object.keys(requested).every(function(k){return obj[k]===requested[k]}) &&
                Object.keys(obj).length===Object.keys(requested).length;
            } catch(e){match=false;}
          } else match=option.value===String(requested);
          if(match){value=option.value;return true;} return false;
        });
        if(value===null){
          if(key==='colormap') {
            var custom=presets.find(function(p){return typeof p==='object'&&p.name===requested});
            if(!custom)return; value=JSON.stringify(custom.scale);
          } else value=typeof requested==='object'?JSON.stringify(requested):String(requested);
          el.add(new Option(typeof requested==='object'?'Custom':String(requested),value));
        }
        if(el.value!==value){el.value=value;changed=true;
          if(key==='mz_bins'||key==='drift_bins')rebinned=true;}
      });
      var signature=JSON.stringify(ranges), patch={};
      var rangeChanged=Object.keys(ranges).length && signature!==window.watersConsoleRanges;
      if(rangeChanged) {
        [['mz_range',['xaxis2','xaxis4'],[HM_MZ_LO0,HM_MZ_HI0]],
         ['drift_range',['yaxis','yaxis2'],[0,nD_native*(pusherUs>0?pusherUs/1000:1)]]].forEach(function(item){
          if(!(item[0] in ranges))return;
          item[1].forEach(function(axis){
            patch[axis+'.range']=(ranges[item[0]]||item[2]).slice();
            patch[axis+'.autorange']=false;
          });
        });
        window.watersConsoleRanges=signature;
      }
      window.watersPendingRender=Promise.resolve(window.watersPendingRender).then(function(){
        if(changed)return rebinned?rebinFull():applyAll();
      }).then(function(){
        if(rangeChanged)return Plotly.relayout(P,patch);
      }).then(function(){
        if(rangeChanged)return updateMarginals();
      }).catch(function(error){window.watersConsoleConfigError=error.message;});
      return true;
    })()
    """.replace("DEFAULTS", values).replace("PRESETS", presets).replace("RANGES", ranges)


class _PreviewSignals(QObject):
    completed = Signal(object)


class _PreviewTask(QRunnable):
    def __init__(self, generation: int, record: dict, config: dict):
        super().__init__()
        self.generation, self.record, self.config = generation, record, config
        self.cancelled = threading.Event()
        self.signals = _PreviewSignals()

    def run(self):
        result = {"generation": self.generation, "errors": []}
        try:
            from deconvovo.io import read_im_txt, read_ms_txt
            from deconvovo.ccs_peak_pick import extract_drift_profile
            ms = self.record.get("ms")
            if ms and Path(ms).is_file():
                try:
                    result["ms"] = read_ms_txt(Path(ms))
                except Exception as exc:
                    result["errors"].append(f"Mass spectrum: {exc}")
            if self.cancelled.is_set():
                return
            im = self.record.get("im")
            if im and Path(im).is_file():
                try:
                    data = read_im_txt(Path(im))
                    if not data.empty:
                        view = (self.config.get("initial_view") or {}).get("mz_range")
                        lo, hi = view if view else (float(data.mz.min()), float(data.mz.max()))
                        visible = data[(data.mz >= lo) & (data.mz <= hi) & (data.intensity > 0)]
                        # This is a labelled raw overview, never a scientific rebin.
                        stride = max(1, (len(visible) + 19999) // 20000)
                        result["im"] = visible.iloc[::stride].to_numpy(copy=True)
                        result["im_count"] = len(visible)
                        result["stride"] = stride
                        result["mz_range"] = (lo, hi)
                        if not self.cancelled.is_set():
                            # Reuse the established extractor, retaining raw-bin
                            # units and unsmoothed intensity throughout the preview.
                            profile = extract_drift_profile(data, (lo + hi) / 2, (hi - lo) / 2)
                            result["mobiligram"] = (np.arange(len(profile)), profile)
                except Exception as exc:
                    result["errors"].append(f"Mobility data: {exc}")
        except Exception as exc:
            result["errors"].append(str(exc))
        if not self.cancelled.is_set():
            try:
                self.signals.completed.emit(result)
            except RuntimeError:
                pass  # A closing application may have destroyed the receiver.


class WatersViewerWidget(QWidget):
    """Selected acquisition inspection; large HTML is loaded from a file URL."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self._generation = 0
        self._closed = False
        self._record: dict = {}
        self._config: dict = {}
        self._viewer_path: Path | None = None
        self._tasks: list[_PreviewTask] = []
        self._pool = QThreadPool.globalInstance()
        self._html_ready = False
        self._html_generation = 0
        self._profile = None
        self.web_view = None
        self._config_timer = QTimer(self)
        self._config_timer.setSingleShot(True)
        self._config_timer.setInterval(120)
        self._config_timer.timeout.connect(self._apply_pending_config)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(7)
        self.tabs = QTabWidget()
        self.tabs.setDocumentMode(True)
        layout.addWidget(self.tabs, 1)
        self.im_stack = QStackedWidget()
        self.empty = QLabel()
        self.empty.setObjectName("watersViewerEmpty")
        self.empty.setWordWrap(True)
        self.empty.setAlignment(Qt.AlignCenter)
        self.empty.setMargin(24)
        self.im_stack.addWidget(self.empty)
        self.overview_plot = PlotCanvas(figsize=(7, 5), nrows=2)
        self.im_stack.addWidget(self.overview_plot)
        if WEBENGINE_AVAILABLE:
            self._profile = QWebEngineProfile(self)
            self._interceptor = _OfflineRequests(self._profile)
            self._profile.setUrlRequestInterceptor(self._interceptor)
            self._profile.downloadRequested.connect(self._download_file)
            self.web_view = QWebEngineView()
            self.web_view.setPage(QWebEnginePage(self._profile, self.web_view))
            self.web_view.loadFinished.connect(self._html_loaded)
            self.im_stack.addWidget(self.web_view)
        self.tabs.addTab(self.im_stack, "IM-MS")
        self.mass_plot = PlotCanvas(figsize=(7, 5))
        self.mobility_plot = PlotCanvas(figsize=(7, 5))
        self.tabs.addTab(self.mass_plot, "Mass profile")
        self.tabs.addTab(self.mobility_plot, "Mobiligram")
        self.status = QLabel()
        self.status.setWordWrap(True)
        self.status.setProperty("subtitle", True)
        layout.addWidget(self.status)
        self.btn_browser = QPushButton("Open enhanced viewer in browser")
        self.btn_browser.setProperty("secondary", True)
        self.btn_browser.clicked.connect(self.open_browser)
        layout.addWidget(self.btn_browser, 0, Qt.AlignRight)
        self.clear("Select an acquisition. Process & view creates its offline interactive viewer.")

    @property
    def blocked_requests(self) -> int:
        return self._interceptor.blocked_requests if self._profile is not None else 0

    def _download_file(self, download):
        suggested = download.downloadFileName()
        suffix = Path(suggested).suffix.lower() or ".csv"
        path, _ = QFileDialog.getSaveFileName(self, "Export viewer data or image", suggested,
                                             f"Viewer export (*{suffix})")
        if path:
            try:
                destination = self._output_path(path, {suffix})
            except ValueError as error:
                download.cancel()
                self.status.setText(str(error))
                self.status.show()
                return
            download.setDownloadDirectory(str(destination.parent))
            download.setDownloadFileName(destination.name)
            if self._html_ready:
                self.status.hide()
            download.accept()
        else:
            download.cancel()

    def _output_path(self, path: Path | str, suffixes: set[str]) -> Path:
        destination = local_path(path).resolve()
        if any(parent.suffix.casefold() == ".raw" for parent in (destination, *destination.parents)):
            raise ValueError("Choose an export location outside RAW acquisitions.")
        if any(self._record.get(key) and destination == Path(self._record[key]).resolve()
               for key in ("ms", "im")):
            raise ValueError("An export cannot replace a source MS or IM table.")
        if not destination.suffix:
            destination = destination.with_suffix(sorted(suffixes)[0])
        if destination.suffix.casefold() not in suffixes:
            raise ValueError("Choose an export filename ending in " + " or ".join(sorted(suffixes)) + ".")
        return destination

    def _cancel_previews(self):
        for task in self._tasks:
            task.cancelled.set()
        self._tasks.clear()

    def clear(self, message: str):
        self._generation += 1
        self._cancel_previews()
        self._record = {}
        self._viewer_path = None
        self._html_ready = False
        self._config_timer.stop()
        if self.web_view:
            self.web_view.stop()
        self.empty.setText(message)
        self.im_stack.setCurrentWidget(self.empty)
        self.status.setText(message)
        self.status.show()
        self.btn_browser.setEnabled(False)
        for plot in (self.overview_plot, self.mass_plot, self.mobility_plot):
            plot.clear()

    def set_run(self, record: dict, config: dict, source: Path | None = None):
        """Display a snapshot and reject results from an older selection."""
        if self._closed:
            return
        self.clear("Loading selected acquisition…")
        self._record = dict(record)
        self._config = deepcopy(config)
        self.setToolTip(str(record.get("run_name", "")))
        viewer = record.get("viewer")
        self._viewer_path = Path(viewer).resolve() if viewer and Path(viewer).is_file() else None
        self.btn_browser.setEnabled(self._viewer_path is not None)
        has_text = any(record.get(key) and Path(record[key]).is_file() for key in ("ms", "im"))
        if self._viewer_path and self.web_view:
            self._html_generation = self._generation
            self.im_stack.setCurrentWidget(self.web_view)
            # setHtml cannot reliably carry multi-megabyte real acquisitions.
            self.web_view.load(QUrl.fromLocalFile(str(self._viewer_path)))
            self.status.setText("Loading offline interactive viewer…")
        elif self._viewer_path:
            self.status.setText("Embedded browser unavailable. Open the offline enhanced HTML in your browser.")
            self.empty.setText("The enhanced viewer is ready. Open it in your browser, or inspect the raw TXT preview here.")
        elif not has_text:
            self.empty.setText("RAW acquisition selected. Process & view converts the data and creates an offline viewer.")
            self.status.setText("No converted spectrum available for this acquisition yet.")
        else:
            self.empty.setText("Reading converted spectra in the background…")
            self.status.setText("Raw preview loading. Process & view creates the full enhanced HTML.")
        if has_text:
            task = _PreviewTask(self._generation, dict(record), deepcopy(config))
            task.signals.completed.connect(self._preview_loaded)
            self._tasks.append(task)
            self._pool.start(task)

    @Slot(object)
    def _preview_loaded(self, result: dict):
        if self._closed or result["generation"] != self._generation:
            return
        self._tasks.clear()
        ms = result.get("ms")
        im = result.get("im")
        for plot in (self.overview_plot, self.mass_plot, self.mobility_plot):
            plot.clear()
        if ms is not None and len(ms[0]):
            for ax in (self.mass_plot.ax(), self.overview_plot.ax(1)):
                ax.plot(*ms, color="#8dcfc4", linewidth=0.85)
                ax.set(xlabel="m/z", ylabel="Raw intensity")
                if result.get("mz_range"):
                    ax.set_xlim(*result["mz_range"])
            self.mass_plot.ax().set_title("Raw mass spectrum · no smoothing", fontsize=10)
        else:
            self.mass_plot.ax().text(0.5, 0.5, "No readable mass spectrum", transform=self.mass_plot.ax().transAxes, ha="center")
        if im is not None and len(im):
            from matplotlib.colors import LogNorm
            self.overview_plot.ax().scatter(im[:, 0], im[:, 1], c=im[:, 2],
                s=2, cmap="viridis", norm=LogNorm(), linewidths=0, rasterized=True)
            self.overview_plot.ax().set(xlabel="m/z", ylabel="Raw drift bin")
            sampling = f" · every {result['stride']}th point" if result["stride"] > 1 else ""
            self.overview_plot.ax().set_title("Raw IM-MS scatter preview" + sampling, fontsize=10)
        else:
            self.overview_plot.ax().text(0.5, 0.5, "No readable mobility data", transform=self.overview_plot.ax().transAxes, ha="center")
        profile = result.get("mobiligram")
        if profile is not None:
            self.mobility_plot.ax().plot(*profile, color="#8dcfc4", linewidth=1.3)
            lo, hi = result["mz_range"]
            self.mobility_plot.ax().set(xlabel="Raw drift bin", ylabel="Summed raw intensity",
                title=f"Unsmoothed mobiligram · m/z {lo:g}–{hi:g}")
        else:
            self.mobility_plot.ax().text(0.5, 0.5, "No readable mobiligram", transform=self.mobility_plot.ax().transAxes, ha="center")
        for plot in (self.overview_plot, self.mass_plot, self.mobility_plot):
            plot.refresh()
        if self.im_stack.currentWidget() != self.web_view:
            if (ms is not None and len(ms[0])) or (im is not None and len(im)):
                self.im_stack.setCurrentWidget(self.overview_plot)
                self.status.setText("Raw preview · drift bins · no smoothing. " +
                    ("Full offline viewer available in browser." if self._viewer_path else "Process & view creates the enhanced viewer."))
            else:
                self.empty.setText("Converted files contain no readable spectrum. Process & view recreates the results.")
                self.status.setText("; ".join(result["errors"]) or "No readable spectrum in the selected files.")

    def apply_config(self, config: dict):
        """Debounce live controls and numeric ranges on the existing plot."""
        _config_script(config)  # validate before accepting a programmatic change
        self._config = deepcopy(config)
        if self._html_ready:
            self._config_timer.start()

    @Slot(bool)
    def _html_loaded(self, success: bool):
        if self._closed or self._html_generation != self._generation or not self._viewer_path:
            return
        loaded_url = self.web_view.url()
        if not loaded_url.isLocalFile() or Path(loaded_url.toLocalFile()).resolve() != self._viewer_path.resolve():
            return
        self._html_ready = success
        if success:
            self.web_view.page().runJavaScript(_CONSOLE_STYLE)
            self._apply_pending_config()
            caption = "Offline interactive viewer · controls and linked profiles use the exported data."
            self.status.setText(caption)
            self.tabs.setTabToolTip(0, caption)
            self.status.hide()
        else:
            self.status.setText("Embedded viewer could not load. Open the enhanced HTML in your browser.")
            self.status.show()
            self.im_stack.setCurrentWidget(self.overview_plot if self._record.get("ms") else self.empty)

    def _apply_pending_config(self):
        if self._closed or not self._html_ready or not self.web_view:
            return
        def applied(result):
            if isinstance(result, str) and result and not self._closed:
                self.status.setText(result)
                self.status.show()
            elif result is True and not self._closed and self.status.text().startswith("The current HTML covers m/z "):
                self.status.hide()

        self.web_view.page().runJavaScript(_config_script(self._config), applied)

    def open_browser(self):
        if self._viewer_path and self._viewer_path.is_file():
            QDesktopServices.openUrl(QUrl.fromLocalFile(str(self._viewer_path)))

    def export_current(self, kind: str):
        """Capture the live Plotly viewport through the offline HTML exporter."""
        if kind not in {"full_png", "heatmap_png", "mz_png", "drift_png",
                        "heatmap_csv", "mz_csv", "drift_csv"}:
            raise ValueError("Unknown current-view export.")
        if self._closed or not self._html_ready or not self.web_view:
            raise ValueError("Process & view before exporting the current view.")
        # Flush a pending control change before queueing the capture. Both
        # scripts execute on the same page; the exporter waits for Plotly.
        self._config_timer.stop()
        self._apply_pending_config()
        script = """
        (function(){
          if(typeof watersExportCurrent!=='function') return false;
          watersExportCurrent(KIND).catch(function(error){
            alert('Unable to export the current view: '+error.message);
          });
          return true;
        })()
        """.replace("KIND", json.dumps(kind))

        def started(available):
            if not available and not self._closed:
                self.status.setText("Regenerate this viewer to enable current-view PNG and CSV exports.")
                self.status.show()

        self.web_view.page().runJavaScript(script, started)

    def export_html(self, destination: Path) -> Path:
        """Copy the existing single-file offline artifact without regeneration."""
        if self._viewer_path is None or not self._viewer_path.is_file():
            raise ValueError("Process & view before exporting portable HTML.")
        destination = self._output_path(destination, {".html", ".htm"})
        if destination.resolve() != self._viewer_path.resolve():
            shutil.copyfile(self._viewer_path, destination)
        return destination

    def wait_for_workers(self, timeout_ms: int = 5000) -> bool:
        self._cancel_previews()
        return self._pool.waitForDone(timeout_ms)

    def closeEvent(self, event):
        self._closed = True
        self._generation += 1
        self._cancel_previews()
        self._config_timer.stop()
        if self.web_view:
            self.web_view.stop()
            # Page must be destroyed before its off-record profile.
            self.web_view.deleteLater()
            self.web_view = None
        super().closeEvent(event)
