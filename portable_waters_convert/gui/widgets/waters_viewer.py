"""Native Qt/Matplotlib Waters inspection; numerical work uses imms_view."""
from __future__ import annotations

from copy import deepcopy
from dataclasses import replace
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
from matplotlib.patheffects import Normal
from matplotlib.ticker import EngFormatter
from PySide6.QtCore import QObject, QRunnable, QThreadPool, QTimer, Qt, QUrl, QSize, Signal, Slot
from PySide6.QtGui import QColor, QDesktopServices, QIcon, QPainter
from PySide6.QtWidgets import QFileDialog, QLabel, QSplitter, QStackedWidget, QTabWidget, QVBoxLayout, QWidget

from deconvovo.paths import local_path
from gui.theme import ACCENT, mpl_rc
from gui.widgets.waters_targets import WatersTargetsWidget, _matches
from gui.widgets import plot_canvas
from gui.widgets.plot_canvas import PlotCanvas

_EXPORT_KINDS = {"full_png", "heatmap_png", "mz_png", "drift_png", "heatmap_csv", "mz_csv", "drift_csv"}
_INTENSITY_CAPTIONS = {
    "linear": "I", "sqrt": r"$\sqrt{I}$", "cbrt": r"$I^{\frac{1}{3}}$",
    "p4": r"$I^{\frac{1}{4}}$", "p5": r"$I^{\frac{1}{5}}$",
    "log2": r"$\log_{2}(I+1)$", "log10": r"$\log_{10}(I+1)$",
}


def _profile_color():
    return ACCENT if plot_canvas._dark_mode else "#4b7454"


def _drift_profile_centres(snapshot):
    centres = snapshot.drift_profile_centres
    return snapshot.drift_centres if centres is None else centres


def _visible_profile_limits(coordinate, intensity, bounds):
    """Scale only the visible line; acquired arrays and export rows stay intact."""
    low, high = sorted(bounds)
    visible = intensity[(coordinate >= low) & (coordinate <= high)]
    if len(coordinate) and high >= coordinate[0] and low <= coordinate[-1]:
        edges = [max(low, coordinate[0]), min(high, coordinate[-1])]
        # A cropped line can cross the view boundary between acquired samples.
        # These graphical endpoint values are used only to set the axis limits.
        visible = np.concatenate((visible, np.interp(edges, coordinate, intensity)))
    minimum, maximum = visible.min(initial=0), visible.max(initial=0)
    padding = .05 * (maximum - minimum) if maximum > minimum else 1.
    return minimum - padding if minimum < 0 else 0., maximum + padding


def _automatic_drift_range(data, mz_range, config):
    """Enclose every positive acquired drift sample, with equal end padding."""
    raw_config = deepcopy(config)
    raw_config.setdefault("defaults", {}).update({"noise_drift": 0, "smooth_drift": {"method": "raw"}})
    centres, intensity, _unit = data.extract_drift_profile(mz_range, raw_config)
    observed = centres[intensity > 0]
    if not len(observed):
        return tuple(data.full_drift_range)
    low, high = float(observed[0]), float(observed[-1])
    padding = .2 * max(abs(low), abs(high))
    return max(data.full_drift_range[0], low - padding), min(data.full_drift_range[1], high + padding)


class _FullResolutionProfile(Normal):
    def draw_path(self, renderer, gc, path, affine, rgbFace=None):
        # Long sorted lines acquire fresh viewport paths during drawing.
        # Keep every acquired/profile sample on that actual render path too.
        path.should_simplify = False
        return super().draw_path(renderer, gc, path, affine, rgbFace)


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
        ax.title.set_color(rc["text.color"])
        for text in ax.texts:
            text.set_color(rc["text.color"])
        for spine in ax.spines.values():
            spine.set_color(rc["axes.edgecolor"])


class _Signals(QObject):
    completed = Signal(object)


class _RenderTask(QRunnable):
    def __init__(self, generation, record, source, data, config, viewport, profile_range=None,
                 automatic_regions=(), automatic_range_cache=None):
        super().__init__()
        self.generation = generation
        self.record, self.source, self.data = record, source, data
        self.config, self.viewport = config, viewport
        self.profile_range = profile_range
        self.automatic_regions = automatic_regions
        self.automatic_range_cache = dict(automatic_range_cache or {})
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
                    snapshot = data.render(self.config, viewport=self.viewport)
                    if self.profile_range is not None:
                        centres, intensity, unit = data.extract_drift_profile(self.profile_range, self.config)
                        snapshot = replace(snapshot, drift_profile_centres=centres, drift_intensity=intensity, drift_unit=unit)
                    ranges, cached_ranges = {}, self.automatic_range_cache
                    for row, mz_range in self.automatic_regions:
                        if self.cancelled.is_set():
                            break
                        if mz_range not in cached_ranges:
                            cached_ranges[mz_range] = _automatic_drift_range(data, mz_range, self.config)
                        ranges[row] = cached_ranges[mz_range]
                    result.update(data=data, snapshot=snapshot, automatic_drift_ranges=ranges,
                                  automatic_range_cache=cached_ranges)
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
                coordinate = _drift_profile_centres(snapshot) if drift else snapshot.ms_mz
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

    def __init__(self, canvas, parent, reset, physical_axes=None):
        self._reset = reset
        self._physical_axes = physical_axes or {}
        super().__init__(canvas, parent, coordinates=True)
        self.setObjectName("watersPlotTools")
        self.setMovable(False)
        self.setIconSize(QSize(18, 18))

    @staticmethod
    def _mouse_event_to_message(event):
        # Keep the native coordinate feedback on one line. The heatmap's
        # formatter already includes its displayed cell value.
        if event.inaxes and event.inaxes.get_navigate():
            return event.inaxes.format_coord(event.xdata, event.ydata)
        return ""

    def _constrain(self, event, axes):
        if axes and axes[0] in self._physical_axes:
            event.key = self._physical_axes[axes[0]]

    def drag_pan(self, event):
        self._constrain(event, self._pan_info.axes if self._pan_info else [])
        super().drag_pan(event)

    def drag_zoom(self, event):
        self._constrain(event, self._zoom_info.axes if self._zoom_info else [])
        super().drag_zoom(event)

    def release_zoom(self, event):
        self._constrain(event, self._zoom_info.axes if self._zoom_info else [])
        super().release_zoom(event)

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
    targets_changed = Signal(object)
    ranges_reset = Signal(object)
    target_focused = Signal(object)
    ccs_profile_visible = Signal(bool)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._generation = 0
        self._closed = self._ready = False
        self._record, self._config = {}, {}
        self._source = self._data = self._snapshot = self._viewport = None
        self._updating_axes = False
        self._drag_active = False
        self._navigation_history = None
        self._target_record = None
        self._target_focus_pending = False
        self._automatic_range_cache = {}
        self._ccs_profiles = []
        self._ccs_profile = None
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
        self.plot.setMinimumHeight(220)
        self.figure, self.canvas = self.plot.fig, self.plot.canvas
        self.figure.clear()
        self.figure.set_layout_engine(None)
        grid = self.figure.add_gridspec(2, 3, width_ratios=[1, 5, 0.15],
            height_ratios=[4, 1.2], left=0.095, right=0.93, bottom=0.14, top=0.96,
            hspace=0.045, wspace=0.06)
        self._heatmap_grid = grid
        self.ax_heatmap = self.figure.add_subplot(grid[0, 1])
        self.ax_drift = self.figure.add_subplot(grid[0, 0], sharey=self.ax_heatmap)
        self.ax_ms = self.figure.add_subplot(grid[1, 1], sharex=self.ax_heatmap)
        self.ax_colorbar = self.figure.add_subplot(grid[0, 2])
        self.canvas.mpl_connect("resize_event", self._layout_spectrum)
        self.plot.axes = np.array([[self.ax_drift, self.ax_heatmap, self.ax_ms, self.ax_colorbar]])
        plot_layout.addWidget(self.plot, 1)
        self.toolbar = _PlotToolbar(self.canvas, page, lambda: self.reset_view(),
                                     {self.ax_ms: "x", self.ax_drift: "y"})
        self.toolbar.setObjectName("watersPlotTools")
        plot_layout.addWidget(self.toolbar)
        self.im_stack.addWidget(page)
        self._plot_page = page
        self.tabs.addTab(self.im_stack, "IM-MS")
        self.mass_plot = PlotCanvas(figsize=(7, 4))
        self.mobility_plot = PlotCanvas(figsize=(7, 4))
        self.tabs.addTab(self.mass_plot, "Mass profile")
        self.tabs.addTab(self.mobility_plot, "Mobiligram")
        self.ccs_plot = PlotCanvas(self, figsize=(3, 2))
        self.ccs_plot.setMinimumWidth(200)
        self.ccs_plot.setMinimumHeight(170)
        self.ccs_toolbar = _PlotToolbar(self.ccs_plot.canvas, self.ccs_plot, self._draw_ccs_profile)
        self.ccs_plot.layout().addWidget(self.ccs_toolbar)
        self.ccs_plot.canvas.mpl_connect("scroll_event", self._scroll_ccs)
        self.ccs_plot.ax().callbacks.connect("xlim_changed", self._ccs_limits_changed)
        self.ccs_plot.hide()
        self.mass_toolbar = _PlotToolbar(self.mass_plot.canvas, self.mass_plot,
            lambda: self.reset_view(), {self.mass_plot.ax(): "x"})
        self.mobility_toolbar = _PlotToolbar(self.mobility_plot.canvas, self.mobility_plot,
            lambda: self.reset_view(), {self.mobility_plot.ax(): "x"})
        self.mass_plot.layout().addWidget(self.mass_toolbar)
        self.mobility_plot.layout().addWidget(self.mobility_toolbar)
        self._toolbars = (self.toolbar, self.mass_toolbar, self.mobility_toolbar)
        self.targets_editor = WatersTargetsWidget()
        self.targets_editor.changed.connect(self._targets_changed)
        self.targets_editor.selected.connect(self._profile_selected)
        self.targets_editor.full_view_requested.connect(self.reset_view)
        self.targets_editor.error.connect(self._profile_error)
        self.target_splitter = QSplitter(Qt.Vertical)
        self.target_splitter.setChildrenCollapsible(False)
        self.target_splitter.addWidget(self.tabs)
        self.target_splitter.addWidget(self.targets_editor)
        self.target_splitter.setStretchFactor(0, 1)
        self.target_splitter.setStretchFactor(1, 0)
        self.target_splitter.setSizes([420, self.targets_editor.preferred_height])
        layout.addWidget(self.target_splitter, 1)
        self.status = QLabel()
        self.status.setWordWrap(True)
        self.status.setProperty("subtitle", True)
        layout.addWidget(self.status)
        for ax in (self.ax_heatmap, self.ax_drift, self.ax_ms):
            ax.callbacks.connect("xlim_changed", self._axes_changed)
            ax.callbacks.connect("ylim_changed", self._axes_changed)
        for canvas in (self.canvas, self.mass_plot.canvas, self.mobility_plot.canvas):
            canvas.mpl_connect("scroll_event", self._scroll_zoom)
            canvas.mpl_connect("button_press_event", self._gesture_started)
            canvas.mpl_connect("button_release_event", self._gesture_finished)
        self.clear("Open RAW… to start, or open an example.")

    @property
    def is_ready(self):
        return self._ready

    def _layout_spectrum(self, _event=None):
        width, height = self.canvas.get_width_height()
        width, height = max(1, width), max(1, height)
        self._heatmap_grid.update(left=max(.095, min(.20, 58 / width)),
            right=min(.93, 1 - 48 / width), bottom=max(.14, min(.38, 56 / height)),
            hspace=.12 if height < 340 else .045)
        # GridSpec.update only repositions pyplot-managed figures in supported
        # Matplotlib versions; this native Figure is deliberately independent.
        for ax in (self.ax_heatmap, self.ax_drift, self.ax_ms, self.ax_colorbar):
            ax.set_position(ax.get_subplotspec().get_position(self.figure))
        self.ax_ms.yaxis.set_label_position("left")
        if width < 600:
            box = self.ax_ms.get_position()
            self.ax_ms.yaxis.set_label_coords(.035, box.y0 + box.height / 2, transform=self.figure.transFigure)
        else:
            self.ax_ms.yaxis.set_label_coords(-.10, .5)

    @property
    def snapshot(self):
        return self._snapshot

    @property
    def active_profile_range(self):
        return self.targets_editor.active_range

    def get_targets(self):
        return self.targets_editor.get_targets()

    def update_target_results(self, records, profiles=()):
        """Merge readonly CCS results without changing the plotted acquisition."""
        self.targets_editor.update_results(records)
        self._ccs_profiles = list(profiles)
        self._draw_ccs_profile()

    def invalidate_target_results(self):
        """Mark displayed CCS results stale after a calibration change."""
        self.targets_editor.invalidate_results()
        self._ccs_profiles = []
        self._ccs_profile = None
        self.ccs_plot.hide()
        self.ccs_profile_visible.emit(False)

    def _draw_ccs_profile(self):
        rows = self.targets_editor.table.selectionModel().selectedRows()
        targets = self.get_targets()
        target = targets.iloc[rows[0].row()] if rows and rows[0].row() < len(targets) else None
        profile = next((item for item in self._ccs_profiles if target is not None
                        and _matches(target, item["target"])), None)
        self._ccs_profile = profile
        if not self._ccs_profiles:
            self.ccs_plot.hide()
            self.ccs_profile_visible.emit(False)
            return
        self.ccs_plot.show()
        self.ccs_profile_visible.emit(True)
        ax = self.ccs_plot.ax()
        ax.clear()
        if profile is None:
            ax.set_axis_off()
            ax.text(.5, .5, "Select a target with CCS", transform=ax.transAxes,
                    ha="center", va="center", fontsize=10, wrap=True,
                    color=mpl_rc(dark=plot_canvas._dark_mode)["text.color"])
        else:
            ax.set_axis_on()
            x, y = profile["ccs"], profile["intensity"]
            ax.plot(x, y, color=_profile_color(), linewidth=1.15,
                    path_effects=[_FullResolutionProfile()])
            from deconvovo.ccs_convert import significant_range
            ax.set_xlim(significant_range(x, y))
            ax.set_ylim(_visible_profile_limits(x, y, ax.get_xlim()))
            ax.set_title(profile["target"]["species"], fontsize=10)
            ax.set_xlabel("CCS (Å²)")
            ax.set_ylabel("Intensity")
            ax.yaxis.set_major_formatter(EngFormatter(sep=""))
            ax.grid(alpha=.15)
        ax.callbacks.connect("xlim_changed", self._ccs_limits_changed)
        self.ccs_toolbar.update()
        self.ccs_plot.refresh()

    def _ccs_limits_changed(self, ax):
        if self._ccs_profile is not None:
            ax.set_ylim(_visible_profile_limits(self._ccs_profile["ccs"],
                                               self._ccs_profile["intensity"], ax.get_xlim()))

    def _scroll_ccs(self, event):
        if event.inaxes is not self.ccs_plot.ax() or self._ccs_profile is None:
            return
        ax = self.ccs_plot.ax()
        self.ccs_toolbar.push_current()
        low, high = ax.get_xlim()
        anchor = event.xdata if event.xdata is not None else (low + high) / 2
        factor = .8 if event.button == "up" else 1.25
        ax.set_xlim(anchor + (low - anchor) * factor, anchor + (high - anchor) * factor)
        self.ccs_toolbar.push_current()
        self.ccs_plot.canvas.draw_idle()

    def _refresh_profile_controls(self):
        self.targets_editor.set_context(self._record.get("raw"),
            tuple(self._snapshot.mz_range) if self._snapshot is not None else None)

    def _targets_changed(self, frame):
        self.targets_changed.emit(frame)
        self._draw_ccs_profile()
        self._request_profile_render()

    def _profile_selected(self, bounds):
        self._draw_ccs_profile()
        self._target_focus_pending = bounds is not None and self._data is None
        if bounds is not None:
            if self._data is not None:
                self._focus_target(bounds)
        self._request_profile_render()

    def _focus_target(self, bounds):
        if self.targets_editor.active_drift_mode == "auto":
            drift = self._automatic_range_cache.get(bounds)
            if drift is None:
                drift = _automatic_drift_range(self._data, bounds, self._config)
                self._automatic_range_cache[bounds] = drift
            self.targets_editor.set_detected_drift_range(drift, self._data.drift_unit)
        else:
            drift = self.targets_editor.active_drift_range
        if drift is None:
            return
        focused = {"mz_range": list(bounds), "drift_range": list(drift)}
        self._config.setdefault("initial_view", {}).update(focused)
        self._viewport = None
        # A target changes the heatmap domain just like an explicit numeric
        # range. Start navigation history on that domain, keeping all RAW
        # arrays available for Home and later targets.
        for toolbar in self._toolbars:
            toolbar.update()
        self._navigation_history = (None, False)
        self.target_focused.emit(deepcopy(focused))

    def _profile_error(self, message):
        self.status.setText(message)
        self.status.show()

    def _request_profile_render(self):
        if not self._record.get("im") or self._closed:
            return
        self._generation += 1
        self._cancel_renders()
        self._set_ready(False)
        self._render_timer.start()

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
        self._target_focus_pending = False
        self._automatic_range_cache = {}
        self.targets_editor.setEnabled(False)
        self._drag_active = False
        self._navigation_history = None
        for toolbar in self._toolbars:
            toolbar.update()
        self._record = {}
        self._ccs_profiles = []
        self._ccs_profile = None
        self.ccs_plot.hide()
        self.ccs_profile_visible.emit(False)
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
        key = str(Path(record.get("raw") or record.get("im") or record.get("run_name", ".")).resolve())
        if key != self._target_record:
            self.targets_editor.set_targets(self.get_targets().iloc[0:0])
            self._target_record = key
        self.clear("Loading converted spectra…")
        self._record, self._config = dict(record), deepcopy(config)
        self._target_focus_pending = self.active_profile_range is not None
        self.targets_editor.set_context(record.get("raw"), None)
        self._source = source
        self.setToolTip(str(record.get("run_name", "")))
        if not record.get("im") or not Path(record["im"]).is_file():
            self.empty.setText("Choose View to prepare and inspect this acquisition.")
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
        old_ranges, new_ranges = self._config.get("initial_view") or {}, config.get("initial_view") or {}
        changed_ranges = [name for name in ("mz_range", "drift_range") if old_ranges.get(name) != new_ranges.get(name)]
        if changed_ranges:
            if self._viewport is not None:
                for name in changed_ranges:
                    self._viewport.pop(name, None)
                self._viewport = self._viewport or None
            self._navigation_history = None
            for toolbar in self._toolbars:
                toolbar.update()
        self._config = config
        self._generation += 1
        self._cancel_renders()
        self._set_ready(False)
        self._render_timer.start()

    def reset_view(self, axes=None):
        self._record_navigation_start()
        self._navigation_history = (None, False)
        axes = ("mz_range", "drift_range") if axes is None else tuple(axes)
        if any(name not in {"mz_range", "drift_range"} for name in axes):
            raise ValueError("Reset axes must be mz_range or drift_range.")
        for name in axes:
            self._config.setdefault("initial_view", {})[name] = None
            if self._viewport is not None:
                self._viewport.pop(name, None)
        self._viewport = self._viewport or None
        self.ranges_reset.emit(axes)
        self._drag_active = False
        self._render_timer.stop()
        if self._record.get("im"):
            self._start_render()

    def _start_render(self):
        if self._closed or not self._record.get("im") or self._drag_active:
            return
        self._generation += 1
        self._cancel_renders()
        self._set_ready(False)
        task = _RenderTask(self._generation, dict(self._record), self._source,
                           self._data, deepcopy(self._config), deepcopy(self._viewport), self.active_profile_range,
                           self.targets_editor.automatic_regions(), self._automatic_range_cache)
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
        self._automatic_range_cache = result.get("automatic_range_cache", {})
        self.targets_editor.set_detected_drift_ranges(result.get("automatic_drift_ranges", {}), self._data.drift_unit)
        if self._target_focus_pending and self.active_profile_range is not None:
            self._target_focus_pending = False
            self._focus_target(self.active_profile_range)
            self._start_render()
            return
        self.targets_editor.setEnabled(True)
        self._draw_snapshot()
        self._refresh_profile_controls()
        if self._navigation_history is not None:
            active, already_pushed = self._navigation_history
            for toolbar in self._toolbars:
                if toolbar is not active or not already_pushed:
                    toolbar.push_current()
            self._navigation_history = None
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
        drift_centres = _drift_profile_centres(snapshot)
        step = float(np.median(np.diff(drift_centres))) if len(drift_centres) > 1 else 0
        self.plot.setToolTip(f"{len(drift_centres)} original drift samples · spacing {step:.6g} {snapshot.drift_unit}. "
                             "The side profile and Mobiligram use the same acquired points; heatmap binning is independent.")
        self._updating_axes = True
        try:
            for ax in (self.ax_heatmap, self.ax_drift, self.ax_ms, self.ax_colorbar):
                ax.clear()
                if ax is not self.ax_colorbar:
                    ax.callbacks.connect("xlim_changed", self._axes_changed)
                    ax.callbacks.connect("ylim_changed", self._axes_changed)
            image = _heatmap(self.ax_heatmap, snapshot, self._config)
            if self.active_profile_range is not None:
                self.ax_heatmap.axvspan(*self.active_profile_range, facecolor=_profile_color(),
                    edgecolor=_profile_color(), alpha=.18, linewidth=1.1)
            self.ax_drift.plot(snapshot.drift_intensity, drift_centres,
                               color=_profile_color(), linewidth=1.1,
                               path_effects=[_FullResolutionProfile()])
            self.ax_ms.plot(snapshot.ms_mz, snapshot.ms_intensity, color=_profile_color(), linewidth=1,
                            marker="o" if len(snapshot.ms_mz) == 1 else None, markersize=3.5,
                            path_effects=[_FullResolutionProfile()])
            profile_limits = _visible_profile_limits(drift_centres, snapshot.drift_intensity, snapshot.drift_range)
            self.ax_drift.set_xlim(profile_limits[::-1])
            self.ax_heatmap.set_xlim(snapshot.mz_range)
            self.ax_heatmap.set_ylim(snapshot.drift_range)
            self.ax_drift.set_ylabel("Drift time (ms)" if snapshot.drift_unit == "ms" else "Drift bin")
            self.ax_ms.set_xlabel("m/z")
            self.ax_ms.set_ylabel("Intensity")
            self.ax_heatmap.tick_params(labelbottom=False, labelleft=False)
            self.ax_heatmap.format_coord = self._heatmap_coordinates
            self.ax_ms.format_coord = lambda x, y: f"m/z {x:.5g} · I {y:.5g}"
            self.ax_drift.format_coord = lambda x, y: f"Drift {y:.5g} {snapshot.drift_unit} · I {x:.5g}"
            self.ax_colorbar.set_navigate(False)
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
                (self.mobility_plot, drift_centres, snapshot.drift_intensity,
                 "Drift time (ms)" if snapshot.drift_unit == "ms" else "Drift bin")):
                plot.clear()
                plot.ax().plot(coordinate, intensity, color=_profile_color(), linewidth=1.1,
                               marker="o" if plot is self.mass_plot and len(coordinate) == 1 else None,
                               markersize=3.5,
                               path_effects=[_FullResolutionProfile()])
                plot.ax().set(xlabel=label, ylabel="Intensity")
                if plot is self.mobility_plot and self.active_profile_range is not None:
                    low, high = self.active_profile_range
                    plot.ax().set_title(f"{self.targets_editor.active_name} · m/z {low:.12g}–{high:.12g}", fontsize=10)
                plot.ax().set_xlim(snapshot.mz_range if plot is self.mass_plot else snapshot.drift_range)
                if plot is self.mobility_plot:
                    plot.ax().set_ylim(profile_limits)
                plot.ax().callbacks.connect("xlim_changed", self._axes_changed)
                plot.ax().format_coord = lambda x, y, label=label: f"{label} {x:.5g} · I {y:.5g}"
                plot.refresh()
            self.canvas.draw_idle()
        finally:
            self._updating_axes = False
        self._layout_spectrum()

    def refresh_theme(self):
        """Refresh presentation only; scientific arrays and viewport stay fixed."""
        for plot in (self.plot, self.mass_plot, self.mobility_plot, self.ccs_plot):
            for ax in plot.fig.axes:
                for line in ax.lines:
                    line.set_color(_profile_color())
            _style_axes(plot.fig, plot.fig.axes)
            plot.refresh()
        for toolbar in self._toolbars:
            toolbar.refresh_icons()
        self.ccs_toolbar.refresh_icons()

    def _heatmap_coordinates(self, mz, drift):
        snapshot = self._snapshot
        if snapshot is None:
            return ""
        message = f"m/z {mz:.5g} · Drift {drift:.4g} {snapshot.drift_unit}"
        x0, x1 = _extent(snapshot.mz_centres, snapshot.baked_mz_range)
        y0, y1 = _extent(snapshot.drift_centres, (0., float(snapshot.drift_centres[0] * 2)))
        if x0 <= mz <= x1 and y0 <= drift <= y1:
            col = min(len(snapshot.mz_centres) - 1, int((mz - x0) / (x1 - x0) * len(snapshot.mz_centres)))
            row = min(len(snapshot.drift_centres) - 1, int((drift - y0) / (y1 - y0) * len(snapshot.drift_centres)))
            label = {"linear": "I", "sqrt": "sqrt(I)", "cbrt": "cbrt(I)",
                     "p4": "I^(1/4)", "p5": "I^(1/5)",
                     "log2": "log2(I+1)", "log10": "log10(I+1)"}[snapshot.transform]
            return message + f" · {label} {snapshot.display2d[row, col]:.5g}"
        return message + " · outside data"

    def _gesture_started(self, event):
        toolbar = event.canvas.toolbar
        if self._snapshot is not None and toolbar in self._toolbars and (toolbar._pan_info or toolbar._zoom_info):
            self._drag_active = True
            self._record_navigation_start()
            self._render_timer.stop()
            self._generation += 1
            self._cancel_renders()
            self._set_ready(False)

    def _record_navigation_start(self):
        if self._snapshot is not None:
            for toolbar in self._toolbars:
                if toolbar._nav_stack() is None:
                    toolbar.push_current()

    def _gesture_finished(self, event):
        if self._drag_active:
            self._drag_active = False
            self._navigation_history = (event.canvas.toolbar, True)
            self._axes_changed(self.ax_heatmap)

    def _axes_changed(self, _axis):
        if self._updating_axes or self._snapshot is None or self._closed:
            return
        self._updating_axes = True
        try:
            if _axis is self.mass_plot.ax():
                self.ax_heatmap.set_xlim(_axis.get_xlim())
            elif _axis is self.mobility_plot.ax():
                self.ax_heatmap.set_ylim(_axis.get_xlim())
        finally:
            self._updating_axes = False
        self._viewport = {"mz_range": list(self.ax_heatmap.get_xlim()),
                          "drift_range": list(self.ax_heatmap.get_ylim())}
        self._generation += 1
        self._cancel_renders()
        self._set_ready(False)
        if self._drag_active:
            self._render_timer.stop()
        else:
            self._render_timer.start()

    def _scroll_zoom(self, event):
        if self._snapshot is None or self._closed or event.inaxes not in (
                self.ax_heatmap, self.ax_drift, self.ax_ms, self.mass_plot.ax(), self.mobility_plot.ax()):
            return
        self._record_navigation_start()
        self._navigation_history = (event.canvas.toolbar, False)
        factor = 0.8 if event.button == "up" else 1.25
        self._updating_axes = True
        try:
            if event.inaxes in (self.ax_heatmap, self.ax_ms, self.mass_plot.ax()):
                low, high = self.ax_heatmap.get_xlim()
                anchor = event.xdata if event.xdata is not None else (low + high) / 2
                self.ax_heatmap.set_xlim(anchor + (low - anchor) * factor, anchor + (high - anchor) * factor)
            if event.inaxes in (self.ax_heatmap, self.ax_drift, self.mobility_plot.ax()):
                low, high = self.ax_heatmap.get_ylim()
                coordinate = event.xdata if event.inaxes is self.mobility_plot.ax() else event.ydata
                anchor = coordinate if coordinate is not None else (low + high) / 2
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
            raise ValueError("View the acquisition before exporting its current plot.")
        suffix = ".csv" if kind.endswith("csv") else ".png"
        if destination is None:
            name = self._record.get("run_name", "waters") + "_current_" + kind.rsplit("_", 1)[0] + suffix
            if kind.startswith("drift") and self.active_profile_range is not None:
                target = self.targets_editor.active_name
                safe_name = "".join(character if character.isalnum() or character in " ._-" else "_" for character in target).strip(" .") or "target"
                name = self._record.get("run_name", "waters") + "_" + safe_name + "_current_drift" + suffix
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
                if self.active_profile_range is not None:
                    ax.axvspan(*self.active_profile_range, facecolor=_profile_color(), edgecolor=_profile_color(), alpha=.18, linewidth=1.1)
                ax.set(xlim=snapshot.mz_range, ylim=snapshot.drift_range, xlabel="m/z",
                    ylabel="Drift time (ms)" if snapshot.drift_unit == "ms" else "Drift bin")
                colorbar = figure.colorbar(image, ax=ax)
                colorbar.set_label(_INTENSITY_CAPTIONS[snapshot.transform], usetex=False)
            elif kind == "mz_png":
                ax.plot(snapshot.ms_mz, snapshot.ms_intensity, color=_profile_color(), linewidth=1,
                        marker="o" if len(snapshot.ms_mz) == 1 else None, markersize=3.5,
                        path_effects=[_FullResolutionProfile()])
                ax.set(xlim=snapshot.mz_range, xlabel="m/z", ylabel="Intensity")
            else:
                ax.plot(_drift_profile_centres(snapshot), snapshot.drift_intensity, color=_profile_color(), linewidth=1,
                        path_effects=[_FullResolutionProfile()])
                ax.set(xlim=snapshot.drift_range,
                    xlabel="Drift time (ms)" if snapshot.drift_unit == "ms" else "Drift bin", ylabel="Intensity")
                ax.set_ylim(_visible_profile_limits(_drift_profile_centres(snapshot), snapshot.drift_intensity,
                                                  snapshot.drift_range))
                if self.active_profile_range is not None:
                    low, high = self.active_profile_range
                    ax.set_title(f"{self.targets_editor.active_name} · m/z {low:.12g}–{high:.12g}", fontsize=10)
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
