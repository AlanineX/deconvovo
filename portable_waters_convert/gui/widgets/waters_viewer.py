"""Native Qt/Matplotlib Waters inspection; numerical work uses imms_view."""
from __future__ import annotations

from copy import deepcopy
import csv
import math
from pathlib import Path
import shutil
import tempfile
import threading

import numpy as np
from matplotlib import colormaps
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.figure import Figure
from matplotlib.ticker import EngFormatter
from PySide6.QtCore import QObject, QRunnable, QThreadPool, QTimer, Qt, QUrl, Signal, Slot
from PySide6.QtGui import QColor, QDesktopServices, QIcon, QPainter
from PySide6.QtWidgets import QFileDialog, QLabel, QStackedWidget, QTabWidget, QVBoxLayout, QWidget

from deconvovo.paths import local_path
from gui.theme import mpl_rc
from gui.widgets import plot_canvas
from gui.widgets.plot_canvas import PlotCanvas

_EXPORT_KINDS = {"full_png", "heatmap_png", "mz_png", "drift_png", "heatmap_csv", "mz_csv", "drift_csv"}
_INTENSITY_CAPTIONS = {
    "linear": "I", "sqrt": r"$\sqrt{I}$", "cbrt": r"$I^{\frac{1}{3}}$",
    "p4": r"$I^{\frac{1}{4}}$", "p5": r"$I^{\frac{1}{5}}$",
    "log2": r"$\log_{2}(I+1)$", "log10": r"$\log_{10}(I+1)$",
}


def _profile_color():
    return "#8dcfc4" if plot_canvas._dark_mode else "#257a72"


def _colormap(config):
    name = config.get("defaults", {}).get("colormap", "Viridis")
    scale = next((p["scale"] for p in config.get("presets", {}).get("colormap", [])
                  if isinstance(p, dict) and p.get("name") == name), None)
    if scale is None:
        for spelling in (name, str(name).lower()):
            if spelling in colormaps:
                return colormaps[spelling]
        from plotly.colors import get_colorscale
        scale = get_colorscale(name)
    from plotly.colors import unlabel_rgb
    stops = [(float(position), tuple(channel / 255 for channel in unlabel_rgb(color))
              if color.startswith("rgb") else color) for position, color in scale]
    return LinearSegmentedColormap.from_list(str(name), stops)


def _extent(centres, bounds):
    if len(centres) > 1:
        half = (centres[1] - centres[0]) / 2
        return float(centres[0] - half), float(centres[-1] + half)
    return bounds


def _heatmap(ax, snapshot, config):
    x = _extent(snapshot.mz_centres, tuple(snapshot.baked_mz_range))
    y = _extent(snapshot.drift_centres, (0., float(snapshot.drift_centres[0] * 2)))
    interpolation = {"off": "nearest", "fast": "bilinear", "best": "bicubic"}.get(
        config.get("defaults", {}).get("render", "off"), "nearest")
    return ax.imshow(snapshot.display2d, extent=(*x, *y), origin="lower", aspect="auto",
                     interpolation=interpolation, cmap=_colormap(config))


def _style_axes(figure, axes):
    rc = mpl_rc(dark=plot_canvas._dark_mode)
    figure.set_facecolor(rc["figure.facecolor"])
    for ax in axes:
        ax.set_facecolor(rc["axes.facecolor"])
        ax.tick_params(colors=rc["xtick.color"], labelsize=9)
        ax.xaxis.get_offset_text().set_color(rc["xtick.color"])
        ax.yaxis.get_offset_text().set_color(rc["ytick.color"])
        for label in (ax.xaxis.label, ax.yaxis.label):
            label.set_color(rc["axes.labelcolor"])
            label.set_size(10)
        for spine in ax.spines.values():
            spine.set_color(rc["axes.edgecolor"])


class _Signals(QObject):
    completed = Signal(object)


class _RenderTask(QRunnable):
    def __init__(self, generation, record, source, data, config, viewport):
        super().__init__()
        self.generation = generation
        self.record, self.source, self.data = record, source, data
        self.config, self.viewport = config, viewport
        self.cancelled = threading.Event()
        self.signals = _Signals()

    def run(self):
        result = {"task": self, "generation": self.generation}
        try:
            if not self.cancelled.is_set():
                from deconvovo.imms_view import load_view_data
                from deconvovo.io import find_pusher_for_run, find_pusher_period
                data = self.data
                if data is None:
                    raw = self.record.get("raw")
                    pusher = find_pusher_period(Path(raw)) if raw else None
                    if pusher is None:
                        pusher = find_pusher_for_run(Path(self.record["im"]).parent,
                            self.record.get("run_name", ""), self.source)
                    data = load_view_data(self.record["im"], self.record.get("ms"),
                                          pusher_us=pusher, config=self.config)
                if not self.cancelled.is_set():
                    result.update(data=data, snapshot=data.render(self.config, viewport=self.viewport))
        except Exception as exc:
            result["error"] = str(exc)
        result["cancelled"] = self.cancelled.is_set()
        try:
            self.signals.completed.emit(result)
        except RuntimeError:
            pass


def _write_csv(snapshot, kind, destination):
    drift_column = "drift_time_ms" if snapshot.drift_unit == "ms" else "drift_bin"
    with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", newline="", dir=destination.parent,
                                     suffix=".csv", delete=False) as stream:
        staged = Path(stream.name)
        try:
            writer = csv.writer(stream)
            if kind == "heatmap_csv":
                writer.writerow(["mz", drift_column, "display_value", "intensity_transform"])
                xs = snapshot.mz_centres[snapshot.visible_mz_slice]
                ys = snapshot.drift_centres[snapshot.visible_drift_slice]
                z = snapshot.display2d[snapshot.visible_drift_slice, snapshot.visible_mz_slice]
                for drift, row in zip(ys, z):
                    writer.writerows((format(mz, ".17g"), format(drift, ".17g"),
                                      format(value, ".17g"), snapshot.transform) for mz, value in zip(xs, row))
            else:
                drift = kind == "drift_csv"
                coordinate = snapshot.drift_centres if drift else snapshot.ms_mz
                values = snapshot.drift_intensity if drift else snapshot.ms_intensity
                low, high = snapshot.drift_range if drift else snapshot.mz_range
                writer.writerow([drift_column if drift else "mz", "display_value"])
                for x, y in zip(coordinate, values):
                    if np.isnan(x) or np.isnan(y):
                        writer.writerow(["", ""])
                    elif low <= x <= high:
                        writer.writerow([format(x, ".17g"), format(y, ".17g")])
            stream.flush()
        except BaseException:
            staged.unlink(missing_ok=True)
            raise
    try:
        staged.replace(destination)
    finally:
        staged.unlink(missing_ok=True)


class _CsvTask(QRunnable):
    def __init__(self, snapshot, kind, destination):
        super().__init__()
        self.snapshot, self.kind, self.destination = snapshot, kind, destination
        self.signals = _Signals()

    def run(self):
        result = {"task": self, "path": str(self.destination), "kind": self.kind}
        try:
            _write_csv(self.snapshot, self.kind, self.destination)
        except Exception as exc:
            result["error"] = str(exc)
        try:
            self.signals.completed.emit(result)
        except RuntimeError:
            pass


class _PlotToolbar(NavigationToolbar2QT):
    # Saving and layout editing must go through the application's guarded
    # exports and physically shared GridSpec, respectively.
    toolitems = tuple(item for item in NavigationToolbar2QT.toolitems
                      if item[0] in {"Home", "Back", "Forward", "Pan", "Zoom"})

    def __init__(self, canvas, parent, reset):
        self._reset = reset
        super().__init__(canvas, parent, coordinates=False)

    def home(self, *args):
        self._reset()

    def _icon(self, name):
        icon = super()._icon(name)
        pixmap = icon.pixmap(24, 24)
        painter = QPainter(pixmap)
        painter.setCompositionMode(QPainter.CompositionMode_SourceIn)
        painter.fillRect(pixmap.rect(), QColor("#e5edf5" if plot_canvas._dark_mode else "#273746"))
        painter.end()
        return QIcon(pixmap)

    def refresh_icons(self):
        for _text, _tooltip, image, callback in self.toolitems:
            self._actions[callback].setIcon(self._icon(image + ".png"))


class WatersViewerWidget(QWidget):
    """Native linked heatmap and profiles sharing one numerical snapshot."""

    ready_changed = Signal(bool)
    render_finished = Signal(object)
    export_finished = Signal(object)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._generation = 0
        self._closed = self._ready = False
        self._record, self._config = {}, {}
        self._source = self._data = self._snapshot = self._viewport = None
        self._updating_axes = False
        self._pending_exports = []
        self._tasks, self._exports = [], []
        self._pool = QThreadPool(self)
        self._pool.setMaxThreadCount(1)
        self._render_timer = QTimer(self)
        self._render_timer.setSingleShot(True)
        self._render_timer.setInterval(120)
        self._render_timer.timeout.connect(self._start_render)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(4)
        self.tabs = QTabWidget()
        self.tabs.setDocumentMode(True)
        layout.addWidget(self.tabs, 1)
        self.im_stack = QStackedWidget()
        self.empty = QLabel()
        self.empty.setWordWrap(True)
        self.empty.setAlignment(Qt.AlignCenter)
        self.empty.setMargin(20)
        self.im_stack.addWidget(self.empty)
        page = QWidget()
        plot_layout = QVBoxLayout(page)
        plot_layout.setContentsMargins(0, 0, 0, 0)
        plot_layout.setSpacing(0)
        self.plot = PlotCanvas(figsize=(7, 4))
        self.figure, self.canvas = self.plot.fig, self.plot.canvas
        self.figure.clear()
        self.figure.set_layout_engine(None)
        grid = self.figure.add_gridspec(2, 3, width_ratios=[1, 5, 0.15],
            height_ratios=[4, 1.2], left=0.095, right=0.93, bottom=0.14, top=0.96,
            hspace=0.045, wspace=0.06)
        self.ax_heatmap = self.figure.add_subplot(grid[0, 1])
        self.ax_drift = self.figure.add_subplot(grid[0, 0], sharey=self.ax_heatmap)
        self.ax_ms = self.figure.add_subplot(grid[1, 1], sharex=self.ax_heatmap)
        self.ax_colorbar = self.figure.add_subplot(grid[0, 2])
        self.plot.axes = np.array([[self.ax_drift, self.ax_heatmap, self.ax_ms, self.ax_colorbar]])
        plot_layout.addWidget(self.plot, 1)
        self.toolbar = _PlotToolbar(self.canvas, page, lambda: self.reset_view())
        self.toolbar.setObjectName("watersPlotTools")
        plot_layout.addWidget(self.toolbar)
        self.im_stack.addWidget(page)
        self._plot_page = page
        self.tabs.addTab(self.im_stack, "IM-MS")
        self.mass_plot = PlotCanvas(figsize=(7, 4))
        self.mobility_plot = PlotCanvas(figsize=(7, 4))
        self.tabs.addTab(self.mass_plot, "Mass profile")
        self.tabs.addTab(self.mobility_plot, "Mobiligram")
        self.status = QLabel()
        self.status.setWordWrap(True)
        self.status.setProperty("subtitle", True)
        layout.addWidget(self.status)
        for ax in (self.ax_heatmap, self.ax_drift, self.ax_ms):
            ax.callbacks.connect("xlim_changed", self._axes_changed)
            ax.callbacks.connect("ylim_changed", self._axes_changed)
        self.canvas.mpl_connect("scroll_event", self._scroll_zoom)
        self.clear("Select an acquisition. Process & view prepares the native IM-MS plots.")

    @property
    def is_ready(self):
        return self._ready

    @property
    def snapshot(self):
        return self._snapshot

    def _set_ready(self, ready):
        if ready != self._ready:
            self._ready = ready
            self.ready_changed.emit(ready)

    def _cancel_renders(self):
        for task in self._tasks:
            task.cancelled.set()

    def clear(self, message):
        self._generation += 1
        self._cancel_renders()
        self._render_timer.stop()
        self._data = self._snapshot = self._viewport = None
        self._record = {}
        for kind, destination in self._pending_exports:
            self.export_finished.emit({"kind": kind, "path": str(destination),
                                       "error": "Selection changed before export completed."})
        self._pending_exports.clear()
        self._set_ready(False)
        self.empty.setText(message)
        self.im_stack.setCurrentWidget(self.empty)
        self.status.setText(message)
        self.status.show()

    def set_run(self, record, config, source=None):
        if self._closed:
            return
        self.clear("Loading converted spectra…")
        self._record, self._config = dict(record), deepcopy(config)
        self._source = source
        self.setToolTip(str(record.get("run_name", "")))
        if not record.get("im") or not Path(record["im"]).is_file():
            self.empty.setText("Process & view converts this acquisition before native inspection.")
            self.status.setText("No converted IM table is available for the selected acquisition.")
            return
        self._start_render()

    def apply_config(self, config):
        if self._closed:
            return
        config = deepcopy(config)
        for name in ("mz_range", "drift_range"):
            bounds = (config.get("initial_view") or {}).get(name)
            if bounds is not None and (len(bounds) != 2 or
                    not all(math.isfinite(float(value)) for value in bounds) or bounds[0] >= bounds[1]):
                raise ValueError(f"The {name.replace('_', ' ')} needs a finite minimum lower than its maximum.")
        if config == self._config:
            return
        if config.get("initial_view", {}) != self._config.get("initial_view", {}):
            self._viewport = None
        self._config = config
        self._generation += 1
        self._cancel_renders()
        self._set_ready(False)
        self._render_timer.start()

    def reset_view(self):
        self._viewport = None
        self._render_timer.stop()
        if self._record.get("im"):
            self._start_render()

    def _start_render(self):
        if self._closed or not self._record.get("im"):
            return
        self._generation += 1
        self._cancel_renders()
        self._set_ready(False)
        task = _RenderTask(self._generation, dict(self._record), self._source,
                           self._data, deepcopy(self._config), deepcopy(self._viewport))
        task.signals.completed.connect(self._render_done)
        self._tasks.append(task)
        self._pool.start(task)

    @Slot(object)
    def _render_done(self, result):
        if result["task"] in self._tasks:
            self._tasks.remove(result["task"])
        if self._closed or result["generation"] != self._generation or result["cancelled"]:
            return
        if result.get("error"):
            self.status.setText(result["error"])
            self.status.show()
            for kind, destination in self._pending_exports:
                self.export_finished.emit({"kind": kind, "path": str(destination), "error": result["error"]})
            self._pending_exports.clear()
            return
        self._data, self._snapshot = result["data"], result["snapshot"]
        self._draw_snapshot()
        self.im_stack.setCurrentWidget(self._plot_page)
        if self._snapshot.drift_unit == "bins":
            self.status.setText("Native view · raw drift-bin centres; pusher timing metadata is unavailable.")
            self.status.show()
        else:
            self.status.hide()
        self._set_ready(True)
        self.render_finished.emit(self._snapshot)
        pending, self._pending_exports = self._pending_exports, []
        for kind, destination in pending:
            self._save_current(kind, destination)

    def _draw_snapshot(self):
        snapshot = self._snapshot
        self._updating_axes = True
        try:
            for ax in (self.ax_heatmap, self.ax_drift, self.ax_ms, self.ax_colorbar):
                ax.clear()
                if ax is not self.ax_colorbar:
                    ax.callbacks.connect("xlim_changed", self._axes_changed)
                    ax.callbacks.connect("ylim_changed", self._axes_changed)
            image = _heatmap(self.ax_heatmap, snapshot, self._config)
            self.ax_drift.plot(snapshot.drift_intensity, snapshot.drift_centres,
                               color=_profile_color(), linewidth=1.1)
            self.ax_ms.plot(snapshot.ms_mz, snapshot.ms_intensity, color=_profile_color(), linewidth=1)
            self.ax_drift.invert_xaxis()
            self.ax_heatmap.set_xlim(snapshot.mz_range)
            self.ax_heatmap.set_ylim(snapshot.drift_range)
            self.ax_drift.set_ylabel("Drift time (ms)" if snapshot.drift_unit == "ms" else "Drift bin")
            self.ax_ms.set_xlabel("m/z")
            self.ax_ms.set_ylabel("Intensity")
            self.ax_heatmap.tick_params(labelbottom=False, labelleft=False)
            self.ax_ms.yaxis.set_major_formatter(EngFormatter(sep=""))
            self.ax_drift.tick_params(axis="x", bottom=False, labelbottom=False)
            self.ax_drift.xaxis.get_offset_text().set_visible(False)
            colorbar = self.figure.colorbar(image, cax=self.ax_colorbar)
            colorbar.set_label(_INTENSITY_CAPTIONS[snapshot.transform], usetex=False)
            for ax in (self.ax_drift, self.ax_ms):
                ax.grid(True, color="#283442", linewidth=0.4, alpha=0.5)
            _style_axes(self.figure, (self.ax_heatmap, self.ax_drift, self.ax_ms, self.ax_colorbar))
            for plot, coordinate, intensity, label in (
                (self.mass_plot, snapshot.ms_mz, snapshot.ms_intensity, "m/z"),
                (self.mobility_plot, snapshot.drift_centres, snapshot.drift_intensity,
                 "Drift time (ms)" if snapshot.drift_unit == "ms" else "Drift bin")):
                plot.clear()
                plot.ax().plot(coordinate, intensity, color=_profile_color(), linewidth=1.1)
                plot.ax().set(xlabel=label, ylabel="Intensity")
                plot.ax().set_xlim(snapshot.mz_range if plot is self.mass_plot else snapshot.drift_range)
                plot.refresh()
            self.canvas.draw_idle()
        finally:
            self._updating_axes = False

    def refresh_theme(self):
        """Refresh presentation only; scientific arrays and viewport stay fixed."""
        for plot in (self.plot, self.mass_plot, self.mobility_plot):
            for ax in plot.fig.axes:
                for line in ax.lines:
                    line.set_color(_profile_color())
            _style_axes(plot.fig, plot.fig.axes)
            plot.refresh()
        self.toolbar.refresh_icons()

    def _axes_changed(self, _axis):
        if self._updating_axes or not self._ready or self._closed:
            return
        self._viewport = {"mz_range": list(self.ax_heatmap.get_xlim()),
                          "drift_range": list(self.ax_heatmap.get_ylim())}
        self._render_timer.start()

    def _scroll_zoom(self, event):
        if not self._ready or event.inaxes not in (self.ax_heatmap, self.ax_drift, self.ax_ms):
            return
        factor = 0.8 if event.button == "up" else 1.25
        self._updating_axes = True
        try:
            if event.inaxes in (self.ax_heatmap, self.ax_ms):
                low, high = self.ax_heatmap.get_xlim()
                anchor = event.xdata if event.xdata is not None else (low + high) / 2
                self.ax_heatmap.set_xlim(anchor + (low - anchor) * factor, anchor + (high - anchor) * factor)
            if event.inaxes in (self.ax_heatmap, self.ax_drift):
                low, high = self.ax_heatmap.get_ylim()
                anchor = event.ydata if event.ydata is not None else (low + high) / 2
                self.ax_heatmap.set_ylim(anchor + (low - anchor) * factor, anchor + (high - anchor) * factor)
        finally:
            self._updating_axes = False
        self._axes_changed(self.ax_heatmap)
        self.canvas.draw_idle()

    def _output_path(self, path, suffixes):
        destination = local_path(path).resolve()
        if any(parent.suffix.casefold() == ".raw" for parent in (destination, *destination.parents)):
            raise ValueError("Choose an export location outside RAW acquisitions.")
        if any(self._record.get(key) and destination == Path(self._record[key]).resolve() for key in ("ms", "im")):
            raise ValueError("An export cannot replace a source MS or IM table.")
        if not destination.suffix:
            destination = destination.with_suffix(sorted(suffixes)[0])
        if destination.suffix.casefold() not in suffixes:
            raise ValueError("Choose an export filename ending in " + " or ".join(sorted(suffixes)) + ".")
        return destination

    def export_current(self, kind, destination=None):
        if kind not in _EXPORT_KINDS:
            raise ValueError("Unknown current-view export.")
        if self._closed or not self._record.get("im"):
            raise ValueError("Process & view before exporting the current view.")
        suffix = ".csv" if kind.endswith("csv") else ".png"
        if destination is None:
            name = self._record.get("run_name", "waters") + "_current_" + kind.rsplit("_", 1)[0] + suffix
            selected, _ = QFileDialog.getSaveFileName(self, "Export current native view", name,
                                                       f"Viewer export (*{suffix})")
            if not selected:
                return None
            destination = Path(selected)
        destination = self._output_path(destination, {suffix})
        if self._render_timer.isActive() or not self._ready:
            self._pending_exports.append((kind, destination))
            if self._render_timer.isActive():
                self._render_timer.stop()
                self._start_render()
        else:
            self._save_current(kind, destination)
        return destination

    def _save_current(self, kind, destination):
        if kind.endswith("csv"):
            task = _CsvTask(self._snapshot, kind, destination)
            task.signals.completed.connect(self._export_done)
            self._exports.append(task)
            self._pool.start(task)
            return
        self.canvas.draw()
        if kind == "full_png":
            self.figure.savefig(destination, dpi=self.figure.dpi * 3)
        else:
            figure = Figure(figsize=self.figure.get_size_inches(), dpi=self.figure.dpi)
            FigureCanvasAgg(figure)
            ax = figure.add_subplot(111)
            snapshot = self._snapshot
            if kind == "heatmap_png":
                image = _heatmap(ax, snapshot, self._config)
                ax.set(xlim=snapshot.mz_range, ylim=snapshot.drift_range, xlabel="m/z",
                    ylabel="Drift time (ms)" if snapshot.drift_unit == "ms" else "Drift bin")
                colorbar = figure.colorbar(image, ax=ax)
                colorbar.set_label(_INTENSITY_CAPTIONS[snapshot.transform], usetex=False)
            elif kind == "mz_png":
                ax.plot(snapshot.ms_mz, snapshot.ms_intensity, color=_profile_color(), linewidth=1)
                ax.set(xlim=snapshot.mz_range, xlabel="m/z", ylabel="Intensity")
            else:
                ax.plot(snapshot.drift_centres, snapshot.drift_intensity, color=_profile_color(), linewidth=1)
                ax.set(xlim=snapshot.drift_range,
                    xlabel="Drift time (ms)" if snapshot.drift_unit == "ms" else "Drift bin", ylabel="Intensity")
            _style_axes(figure, figure.axes)
            figure.tight_layout()
            figure.savefig(destination, dpi=figure.dpi * 3)
            figure.clear()
        self.export_finished.emit({"path": str(destination), "kind": kind})

    @Slot(object)
    def _export_done(self, result):
        if result["task"] in self._exports:
            self._exports.remove(result["task"])
        if result.get("error") and not self._closed:
            self.status.setText(result["error"])
            self.status.show()
        self.export_finished.emit({key: value for key, value in result.items() if key != "task"})

    def open_browser(self):
        """Optional standalone HTML artifact; the GUI plot never uses it."""
        viewer = self._record.get("viewer")
        if viewer and Path(viewer).is_file():
            QDesktopServices.openUrl(QUrl.fromLocalFile(str(Path(viewer).resolve())))

    def export_html(self, destination):
        viewer = self._record.get("viewer")
        if not viewer or not Path(viewer).is_file():
            raise ValueError("Export portable HTML from the workspace's converted-data action.")
        destination = self._output_path(destination, {".html", ".htm"})
        if destination != Path(viewer).resolve():
            shutil.copyfile(viewer, destination)
        return destination

    def wait_for_workers(self, timeout_ms=5000):
        self._cancel_renders()
        return self._pool.waitForDone(timeout_ms)

    def closeEvent(self, event):
        self._closed = True
        self._generation += 1
        self._cancel_renders()
        self._render_timer.stop()
        for plot in (self.plot, self.mass_plot, self.mobility_plot):
            plot.canvas._draw_pending = False
        super().closeEvent(event)
