"""Optional CCS: build and review a curve, then reuse it on current targets."""
from __future__ import annotations

import math
import json
from pathlib import Path

import numpy as np
import pandas as pd
from PySide6.QtCore import Qt, Signal, QSignalBlocker, QTimer
from PySide6.QtGui import QColor, QFont, QIcon, QPainter, QPalette, QTextDocument
from PySide6.QtWidgets import (
    QAbstractItemView, QBoxLayout, QCheckBox, QComboBox, QFileDialog, QHeaderView, QHBoxLayout, QLabel, QPushButton,
    QProgressBar, QSizePolicy, QSplitter, QStyle, QStyleOptionHeader, QStyledItemDelegate, QVBoxLayout, QWidget,
)

from deconvovo.paths import local_path
from gui.theme import ACCENT, ERROR, SUCCESS, TEXT_DIM, WARNING, _DARK, _LIGHT, mpl_rc
from gui.widgets.csv_table import CsvTable
from gui.widgets.plot_canvas import PlotCanvas
from gui.widgets.worker import Worker


REFERENCE_COLUMNS = ["name", "mz_min", "mz_max", "peak_select", "selected_td_ms", "selected_peak",
                     "mw", "z", "CCS_reported", "raw_path"]
REFERENCE_REQUIRED = {"name", "mw", "z", "mz_min", "mz_max", "CCS_reported", "raw_path"}
RESULT_COLUMNS = {"selected_peak", "selected_td_ms"}


class _ReferenceHeader(QHeaderView):
    def paintSection(self, painter, rect, logical_index):
        key = self.model().headerData(logical_index, Qt.Horizontal)
        labels = {"name": "Name", "mz_min": "m/z min", "mz_max": "m/z max",
                  "peak_select": "Selection", "selected_peak": "Peak", "mw": "Mass (Da)",
                  "z": "Charge", "raw_path": "Standard RAW", "CCS_reported": "CCS_reported",
                  "selected_td_ms": "Selected t<sub>d</sub> (ms)"}
        labels["drift_time_ms"] = labels["selected_td_ms"]
        if key not in labels:
            return super().paintSection(painter, rect, logical_index)
        option = QStyleOptionHeader()
        self.initStyleOption(option)
        self.initStyleOptionForIndex(option, logical_index)
        option.rect = rect
        painter.save()
        self.style().drawControl(QStyle.CE_HeaderSection, option, painter, self)
        document = QTextDocument()
        document.setDocumentMargin(0)
        font = QFont(self.font())
        font.setPointSizeF(10)
        font.setWeight(QFont.DemiBold)
        document.setDefaultFont(font)
        color = option.palette.color(QPalette.ButtonText).name()
        document.setHtml(f'<span style="color:{color}">{labels[key]}</span>')
        document.adjustSize()
        size = document.size()
        painter.setClipRect(rect)
        painter.translate(rect.x() + (rect.width() - size.width()) / 2,
                          rect.y() + (rect.height() - size.height()) / 2)
        document.drawContents(painter)
        painter.restore()


class _DriftValueDelegate(QStyledItemDelegate):
    def displayText(self, value, _locale):
        try:
            return f"{float(value):.1f}" if str(value).strip() else ""
        except (TypeError, ValueError):
            return str(value)


class _PeakChoiceDelegate(QStyledItemDelegate):
    def createEditor(self, parent, _option, _index):
        editor = QComboBox(parent)
        editor.addItem("Highest", "highest")
        editor.addItem("Last", "last")
        editor.addItem("Manual", "manual")
        editor.setItemData(0, "Tallest detected drift peak that passes the thresholds.", Qt.ToolTipRole)
        editor.setItemData(1, "Latest detected drift peak that passes the thresholds.", Qt.ToolTipRole)
        editor.setItemData(2, "Process references, then click a detected drift peak to select it.", Qt.ToolTipRole)
        return editor

    def setEditorData(self, editor, index):
        mode = str(index.data() or "highest").lower()
        editor.setCurrentIndex(2 if mode == "manual" else 1 if mode in {"last", "latest"} else 0)

    def setModelData(self, editor, model, index):
        current = str(index.data()).lower()
        current = {"tallest": "highest", "latest": "last"}.get(current, current)
        if current != editor.currentData():
            model.setData(index, editor.currentData(), Qt.EditRole)

    def displayText(self, value, _locale):
        return {"highest": "Highest", "tallest": "Highest", "last": "Last", "latest": "Last", "manual": "Manual"}.get(
            str(value).lower(), str(value))


class _RawPathDelegate(QStyledItemDelegate):
    def initStyleOption(self, option, index):
        super().initStyleOption(option, index)
        option.textElideMode = Qt.ElideMiddle


class CcsWorkbench(QWidget):
    """Workspace-owned panel; acquisition and analyte paths come from the viewer."""

    targets_applied = Signal(object)
    calibration_changed = Signal(object)
    busy_changed = Signal(bool)

    def __init__(self, log_panel, parent=None):
        super().__init__(parent)
        self._log = log_panel
        self._worker = None
        self._busy = False
        self._dark = True
        self._calibration = None
        self._generated = False
        self._preview_calibration = None
        self._curve_points = []
        self._record = None
        self._data_dir = None
        self._out_dir = Path.home() / "DeconVoVo" / "CCS"
        self._targets = pd.DataFrame()
        self._last_target_request = None
        self._target_refresh = QTimer(self)
        self._target_refresh.setSingleShot(True)
        self._target_refresh.setInterval(250)
        self._target_refresh.timeout.connect(self._refresh_target_ccs)
        self._profiles = []
        self._processed_signature = None
        self._processed_directories = []
        self._processing_signature = None
        self._reference_raw = None
        self._toolbar_groups = []

        layout = QVBoxLayout(self)
        layout.setContentsMargins(8, 4, 8, 4)
        layout.setSpacing(4)
        help_text = QLabel("Build a curve from known standards, or load one measured under matching conditions.")
        help_text.setWordWrap(True)
        help_text.setStyleSheet(f"color: {TEXT_DIM};")

        action_bar = QWidget()
        action_bar.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Maximum)
        actions = QBoxLayout(QBoxLayout.LeftToRight, action_bar)
        actions.setContentsMargins(0, 0, 0, 0)
        self.btn_process = QPushButton("Process references")
        self.btn_process.setProperty("accentOutline", True)
        self._refresh_process_icon()
        self.btn_process.setProperty("compact", True)
        self.btn_process.setToolTip("Extract each ion's drift profile and detect its candidate peaks before fitting.")
        self.btn_process.clicked.connect(self.process_references)
        self.btn_build = QPushButton("Fit calibration")
        self.btn_build.setProperty("primary", True)
        self.btn_build.setProperty("compact", True)
        self.btn_build.clicked.connect(self.build_calibration)
        self.btn_load = QPushButton("Load calibration")
        self.btn_load.setProperty("secondary", True)
        self.btn_load.setProperty("compact", True)
        self.btn_load.clicked.connect(self.load_calibration)
        self.btn_save = QPushButton("Save calibration")
        self.btn_save.setProperty("secondary", True)
        self.btn_save.setProperty("compact", True)
        self.btn_save.clicked.connect(self.save_calibration)
        file_actions, fit_actions = QWidget(), QWidget()
        files, fits = QHBoxLayout(file_actions), QHBoxLayout(fit_actions)
        for row in (files, fits):
            row.setContentsMargins(0, 0, 0, 0)
        files.addWidget(self.btn_load)
        files.addWidget(self.btn_save)
        files.addStretch()
        fits.addStretch()
        fits.addWidget(self.btn_process)
        fits.addWidget(self.btn_build)
        actions.addWidget(file_actions, 1)
        actions.addWidget(fit_actions)
        self._toolbar_groups.append((actions, (file_actions, fit_actions)))
        layout.addWidget(action_bar)

        references_page = QWidget()
        references_layout = QVBoxLayout(references_page)
        references_layout.setContentsMargins(0, 4, 0, 0)
        references_layout.setSpacing(4)
        source_actions = QHBoxLayout()
        self.btn_add_raw = QPushButton("Add standard RAW")
        self.btn_add_raw.setProperty("secondary", True)
        self.btn_add_raw.setProperty("compact", True)
        self.btn_add_raw.clicked.connect(self.add_reference_raw)
        source_actions.addWidget(self.btn_add_raw)
        self.btn_add_ion = QPushButton("Add ion")
        self.btn_add_ion.setProperty("secondary", True)
        self.btn_add_ion.setProperty("compact", True)
        self.btn_add_ion.clicked.connect(self.add_reference_ion)
        self.btn_add_ion.hide()
        source_actions.addWidget(QLabel("Model:"))
        self.combo_method = QComboBox()
        self.combo_method.addItem("Power law", "direct")
        self.combo_method.addItem("Two-step", "twostep")
        self.combo_method.currentIndexChanged.connect(self._references_changed)
        source_actions.addWidget(self.combo_method)
        self.show_extraction = QCheckBox("Extraction options")
        self.show_extraction.setToolTip("Show legacy centre/width, thresholds, formula and the selected raw drift bin.")
        self.show_extraction.toggled.connect(self._show_reference_columns)
        source_actions.addWidget(self.show_extraction)
        source_actions.addStretch()
        self.reference_source = QLabel("Import a CSV, or add a standard acquisition and fill in its ions.")
        self.reference_source.setStyleSheet(f"color: {TEXT_DIM};")
        self.reference_source.setWordWrap(True)
        references_layout.addLayout(source_actions)
        self.reference_source.hide()
        self.references = CsvTable(editable=True)
        self.references._lbl.hide()
        self.references.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.references.table.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.references.table.itemSelectionChanged.connect(self._reference_selection_changed)
        self.references.table.currentCellChanged.connect(self._active_reference_changed)
        self.references.table.setHorizontalHeader(_ReferenceHeader(Qt.Horizontal, self.references.table))
        self.references.table.horizontalHeader().setStretchLastSection(False)
        self.references.table.horizontalHeader().setSectionResizeMode(QHeaderView.Interactive)
        self.references.table.horizontalHeader().setMinimumSectionSize(60)
        self.references.table.horizontalHeader().setMinimumHeight(34)
        self._raw_path_delegate = _RawPathDelegate(self.references.table)
        self._peak_choice_delegate = _PeakChoiceDelegate(self.references.table)
        self._drift_value_delegate = _DriftValueDelegate(self.references.table)
        self.references.load_dataframe(pd.DataFrame(columns=REFERENCE_COLUMNS))
        self._size_reference_columns()
        # New rows inherit the current standard, while allowing a blank sheet.
        for button in self.references.findChildren(QPushButton):
            if button.text() == "+ Row":
                button.clicked.disconnect()
                button.clicked.connect(self.add_reference_ion)
                self.btn_add_ion = button
            elif button.text() == "Load CSV":
                button.clicked.disconnect()
                button.clicked.connect(self.load_reference_csv)
        self.references.data_changed.connect(self._references_changed)
        self.references.table.cellChanged.connect(self._peak_mode_edited)
        self.references.table.cellClicked.connect(self._reference_clicked)
        selection_actions = self.references.layout().itemAt(0).layout()
        self.combo_selection_scope = QComboBox()
        self.combo_selection_scope.addItem("Selected ions", "selected")
        self.combo_selection_scope.addItem("All ions", "all")
        self.combo_selection_scope.setToolTip("Ctrl/Shift-click table rows to change several ions together.")
        self.combo_selection_scope.currentIndexChanged.connect(self._update_actions)
        self.btn_highest = QPushButton("Highest")
        self.btn_last = QPushButton("Last")
        for button, mode in ((self.btn_highest, "highest"), (self.btn_last, "last")):
            button.setProperty("secondary", True)
            button.setProperty("compact", True)
            button.clicked.connect(lambda _checked=False, choice=mode: self.set_peak_selection(choice))
        self.references.layout().removeItem(selection_actions)
        selection_actions.setParent(None)
        selection_bar = QWidget(self.references)
        selection_bar.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Maximum)
        selection_layout = QBoxLayout(QBoxLayout.LeftToRight, selection_bar)
        selection_layout.setContentsMargins(0, 0, 0, 0)
        peak_group, sheet_group = QWidget(), QWidget()
        peak_actions = QHBoxLayout(peak_group)
        peak_actions.setContentsMargins(0, 0, 0, 0)
        for control in (self.combo_selection_scope, self.btn_highest, self.btn_last):
            peak_actions.addWidget(control)
        peak_actions.addStretch()
        sheet_group.setLayout(selection_actions)
        selection_layout.addWidget(peak_group)
        selection_layout.addWidget(sheet_group, 1)
        self.references.layout().insertWidget(0, selection_bar)
        self._toolbar_groups.append((selection_layout, (peak_group, sheet_group)))
        self.references.table.setMinimumHeight(120)
        references_layout.addWidget(self.references, 1)

        review_page = QWidget()
        review_layout = QVBoxLayout(review_page)
        review_layout.setContentsMargins(0, 0, 0, 0)
        review_layout.setSpacing(4)
        peak_controls = QWidget(review_page)
        peak_actions = QHBoxLayout(peak_controls)
        peak_actions.setContentsMargins(0, 0, 0, 0)
        peak_actions.addWidget(QLabel("Reference:"))
        self.combo_reference = QComboBox()
        self.combo_reference.setMinimumWidth(100)
        self.combo_reference.currentIndexChanged.connect(self._reference_combo_changed)
        peak_actions.addWidget(self.combo_reference, 1)
        peak_actions.addWidget(QLabel("Drift peak:"))
        self.combo_peak = QComboBox()
        self.combo_peak.setMinimumWidth(100)
        peak_actions.addWidget(self.combo_peak)
        self.combo_peak.activated.connect(self.use_selected_peak)
        self.btn_use_peak = QPushButton("Use peak")
        self.btn_use_peak.setProperty("secondary", True)
        self.btn_use_peak.setProperty("compact", True)
        self.btn_use_peak.setToolTip("Choose this acquired peak in Manual mode, then fit the calibration.")
        self.btn_use_peak.clicked.connect(self.use_selected_peak)
        peak_actions.addWidget(self.btn_use_peak)
        peak_controls.hide()  # Table selection and direct peak clicks replace duplicate selectors.
        self.peak_status = QLabel("")
        self.peak_status.setWordWrap(True)
        self.peak_status.setTextFormat(Qt.RichText)
        self.peak_status.hide()
        self.plot = PlotCanvas(figsize=(9, 3.2), nrows=1, ncols=2)
        self.plot.canvas.mpl_connect("button_press_event", self._profile_clicked)
        self.plot.setMinimumHeight(200)
        review_layout.addWidget(self.plot, 2)
        self.summary = CsvTable(editable=False)
        self.summary._lbl.hide()
        self.summary.table.setHorizontalHeader(_ReferenceHeader(Qt.Horizontal, self.summary.table))
        self.summary.table.horizontalHeader().setSectionResizeMode(QHeaderView.Interactive)
        self.summary.table.horizontalHeader().setStretchLastSection(False)
        self.summary.table.horizontalHeader().setMinimumHeight(34)
        for button in self.summary.findChildren(QPushButton):
            if button.text() == "Load CSV":
                button.hide()
        self.summary.table.setMinimumHeight(90)
        self.summary.table.cellClicked.connect(self._summary_clicked)
        self.summary.hide()  # Kept as a result/export model; the editable table is the review table.
        self.splitter = QSplitter(Qt.Vertical)
        self.splitter.setChildrenCollapsible(False)
        self.splitter.addWidget(review_page)
        self.splitter.addWidget(references_page)
        self.splitter.setStretchFactor(0, 1)
        self.splitter.setStretchFactor(1, 1)
        self.splitter.setSizes([280, 360])
        layout.addWidget(self.splitter, 1)

        apply_actions = QHBoxLayout()
        self.target_status = QLabel("Open an acquisition and add targets below its 2D spectrum.")
        self.target_status.setWordWrap(True)
        self.target_status.setStyleSheet(f"color: {TEXT_DIM};")
        apply_actions.addWidget(self.target_status, 1)
        self.btn_apply = QPushButton("Generate CCS")
        self.btn_apply.setProperty("primary", True)
        self.btn_apply.setProperty("compact", True)
        self.btn_apply.setToolTip("Accept this fit and generate interactive CCS profiles for the current targets.")
        self.btn_apply.clicked.connect(self.generate_ccs)
        apply_actions.addWidget(self.btn_apply)
        layout.addLayout(apply_actions)
        self.status = QLabel("No calibration loaded.")
        self.status.setWordWrap(True)
        layout.addWidget(self.status)
        self.progress = QProgressBar()
        self.progress.setRange(0, 0)
        self.progress.setFixedHeight(5)
        self.progress.setTextVisible(False)
        self.progress.hide()
        layout.addWidget(self.progress)
        self._update_actions()

    def set_output_dir(self, directory):
        self._out_dir = local_path(directory).resolve()

    def set_acquisition(self, record, data_dir):
        self._record = dict(record) if record else None
        self._data_dir = local_path(data_dir).resolve() if data_dir else None
        self._update_actions()
        self._queue_target_ccs()

    def set_targets(self, frame):
        self._targets = frame.copy() if frame is not None else pd.DataFrame()
        self._update_actions()
        self._queue_target_ccs()

    def _queue_target_ccs(self, *_args):
        self._target_refresh.start()

    def _target_request_key(self, targets):
        raw = local_path(self._record["raw"])
        files = [self._data_dir / f"{raw.stem}_{kind}.txt" for kind in ("ms", "im")]
        sources = tuple((str(path), path.stat().st_size, path.stat().st_mtime_ns) for path in files)
        fields = [key for key in ("species", "mw", "z", "mz", "mz_window", "mz_min", "mz_max",
                                  "min_ms_intensity", "formula", "raw_path") if key in targets]
        return (str(raw), str(self._out_dir), sources,
                targets[fields].to_json(orient="split", index=False, double_precision=15))

    def _refresh_target_ccs(self):
        if self.is_busy or not self.isEnabled():
            self._target_refresh.start()
            return
        if (self._calibration is None or not self._generated or not self._record
                or not self._record.get("raw") or not self._data_dir or self._targets.empty):
            return
        required = ("mw", "z", "mz_min", "mz_max")
        if not set(required) <= set(self._targets.columns):
            return
        values = self._targets[list(required)].apply(pd.to_numeric, errors="coerce")
        valid = (np.isfinite(values).all(axis=1) & values.gt(0).all(axis=1)
                 & values["z"].mod(1).eq(0) & values["mz_max"].gt(values["mz_min"]))
        targets = self._targets.loc[valid].copy()
        if targets.empty:
            self.status.setText("Enter target mass and charge to calculate CCS automatically.")
            return
        try:
            if self._target_request_key(targets) != self._last_target_request:
                self.apply_targets(targets)
        except OSError:
            pass  # Converted tables are not ready yet; acquisition loading will queue another update.

    def set_references(self, frame, csv_path=None, *, preserve_calibration=False):
        self._profiles = []
        self._preview_calibration = None
        self._processed_directories = []
        self._processed_signature = None
        self._reference_raw = None
        self._populate_references()
        frame = frame.copy()
        for key in RESULT_COLUMNS:
            if key in frame:
                frame[key] = ""
        self.references.load_dataframe(self._display_references(frame), local_path(csv_path) if csv_path else None)
        if csv_path:
            self.reference_source.setText(local_path(csv_path).name)
            self.reference_source.setToolTip(str(local_path(csv_path)))
        self._size_reference_columns()
        self._show_reference_columns()
        if preserve_calibration:
            self._draw_curve()
            self._draw_profile()
            self._update_actions()
        else:
            self._references_changed()

    def clear_references(self, *, preserve_calibration=True):
        self._reference_raw = None
        self.reference_source.setText("Import a CSV, or add a standard acquisition and fill in its ions.")
        self.set_references(pd.DataFrame(columns=REFERENCE_COLUMNS),
                            preserve_calibration=preserve_calibration)

    def get_calibration(self):
        return self._calibration

    def _set_status_color(self, role):
        self._status_role = role
        color = (_DARK if self._dark else _LIGHT)[role]
        self.status.setStyleSheet(f"color: {color};")

    def refresh_theme(self, dark):
        self._dark = bool(dark)
        self._set_status_color(getattr(self, "_status_role", "text_dim"))
        self._refresh_process_icon()
        self._draw_curve()
        self._draw_profile()

    def _refresh_process_icon(self):
        pixmap = self.style().standardIcon(QStyle.SP_MediaPlay).pixmap(16, 16)
        painter = QPainter(pixmap)
        painter.setCompositionMode(QPainter.CompositionMode_SourceIn)
        painter.fillRect(pixmap.rect(), QColor(self._plot_colors()[1]))
        painter.end()
        self.btn_process.setIcon(QIcon(pixmap))

    def _plot_colors(self):
        return ("#a4b3d0", ACCENT, mpl_rc(True)["text.color"]) if self._dark else (
            "#657d69", "#4b7454", mpl_rc(False)["text.color"])

    @property
    def is_busy(self):
        return self._busy or bool(self._worker and self._worker.isRunning())

    def resizeEvent(self, event):
        super().resizeEvent(event)
        for layout, groups in self._toolbar_groups:
            required = sum(group.sizeHint().width() for group in groups) + layout.spacing()
            direction = QBoxLayout.TopToBottom if layout.parentWidget().width() < required else QBoxLayout.LeftToRight
            if layout.direction() != direction:
                layout.setDirection(direction)
            layout.parentWidget().setMinimumHeight(layout.sizeHint().height())

    def _update_actions(self, *_args):
        busy = self._busy
        ready = self._calibration is not None
        self.btn_process.setEnabled(not busy)
        self.btn_process.setText("Processing…" if busy and getattr(self, "_operation", "") == "_process" else "Process references")
        self.btn_build.setText("Fitting…" if busy and getattr(self, "_operation", "") == "_fit" else "Fit calibration")
        selected = [profile for profile in self._profiles if profile.get("selected_bin") is not None]
        pending = any(profile.get("selection_error") or
                      profile.get("mode") == "manual" and profile.get("selected_bin") is None
                      for profile in self._profiles)
        self.btn_build.setEnabled(not busy and len(selected) >= 3 and not pending)
        self.btn_build.setToolTip("Fit the selected drift peaks." if len(selected) >= 3 and not pending else
                                 "Process references and select at least three peaks; choose every Manual reference.")
        self.btn_load.setEnabled(not busy)
        self.btn_add_raw.setEnabled(not busy)
        self.btn_add_ion.setEnabled(not busy)
        self.combo_method.setEnabled(not busy)
        self.btn_save.setEnabled(not busy and ready)
        self.references.setEnabled(not busy)
        self.combo_selection_scope.setEnabled(not busy)
        scope_count = (self.references.table.rowCount() if self.combo_selection_scope.currentData() == "all"
                       else len(self._selected_reference_rows()))
        for button in (self.btn_highest, self.btn_last):
            button.setEnabled(not busy and scope_count > 0)
            button.setToolTip(f"Set {scope_count} {'ion' if scope_count == 1 else 'ions'} to {button.text()}; update peak previews immediately.")
        index = self.combo_reference.currentIndex()
        editable_peak = (0 <= index < len(self._profiles)
                         and self._profiles[index]["row_index"] < self.references.table.rowCount())
        self.btn_use_peak.setEnabled(not busy and self.combo_peak.currentIndex() >= 0 and editable_peak)
        self.btn_apply.setEnabled(not busy and ready and bool(self._record) and not self._targets.empty)
        if self.btn_process.property("accentOutline"):
            self.btn_process.setProperty("accentOutline", False)
            self.btn_process.style().unpolish(self.btn_process)
            self.btn_process.style().polish(self.btn_process)
        next_action = (self.btn_apply if not self._generated else None) if ready else (
            self.btn_build if len(selected) >= 3 and not pending else None if self._profiles else self.btn_process)
        for button in (self.btn_process, self.btn_build, self.btn_apply):
            primary = button is next_action and not busy
            if button.property("primary") != primary:
                button.setProperty("primary", primary)
                button.style().unpolish(button)
                button.style().polish(button)
        count = len(self._targets)
        acquisition = (self._record or {}).get("run_name", "")
        self.target_status.setText(f"{count} targets in {acquisition}" if acquisition else
                                   "Open an acquisition and add targets below its 2D spectrum.")

    def _references_changed(self):
        frame = self.references.to_dataframe()
        displayed = self._display_references(frame)
        if list(frame.columns) != list(displayed.columns):
            self.references.load_dataframe(displayed, self.references._csv_path)
            self._size_reference_columns()
        self._show_reference_columns()
        if self._processed_signature is not None and self._reference_signature(displayed) != self._processed_signature:
            self._profiles = []
            self._processed_directories = []
            self._processed_signature = None
            self._preview_calibration = None
            self.combo_reference.clear()
            self.summary.load_dataframe(pd.DataFrame())
            self.references._loading = True
            for key in RESULT_COLUMNS:
                column = displayed.columns.get_loc(key)
                for row in range(self.references.table.rowCount()):
                    self.references.table.item(row, column).setText("")
            self.references._loading = False
            self.status.setText("References changed. Process references again before fitting.")
        elif self._profiles:
            self._refresh_selections()
        if self._calibration is not None:
            self._generated = False
            self.status.setText("Peak preview updated · Fit calibration, then Generate CCS." if self._profiles else
                                "References changed · Process references before fitting.")
            self._set_status_color("warning")
            if self._profiles:
                self._preview_calibration = self._calibration
            self._calibration = None
            self.calibration_changed.emit(None)
        if self._profiles:
            missing_manual = sum(profile.get("mode") == "manual" and profile.get("selected_bin") is None
                                 for profile in self._profiles)
            self.status.setText(f"Choose a numbered drift peak for {missing_manual} Manual ion{'s' if missing_manual != 1 else ''} before fitting."
                                if missing_manual else "Peak preview updated · Fit calibration, then Generate CCS.")
            self._set_status_color("warning")
        self._draw_curve()
        self._draw_profile()
        self._update_actions()

    def _show_reference_columns(self, *_args):
        blocker = QSignalBlocker(self.references.table)
        for column in range(self.references.table.columnCount()):
            key = self.references.table.horizontalHeaderItem(column).text()
            self.references.table.setColumnHidden(column, key not in REFERENCE_COLUMNS
                                                  and not self.show_extraction.isChecked())
            self.references.table.setItemDelegateForColumn(
                column, self._raw_path_delegate if key == "raw_path" else
                self._peak_choice_delegate if key == "peak_select" else
                self._drift_value_delegate if key == "selected_td_ms" else None)
            if key in RESULT_COLUMNS:
                for row in range(self.references.table.rowCount()):
                    item = self.references.table.item(row, column)
                    item.setFlags(item.flags() & ~Qt.ItemIsEditable)
        del blocker

    def _size_reference_columns(self):
        widths = {"name": 115, "mz_min": 80, "mz_max": 80, "peak_select": 90,
                  "selected_td_ms": 130, "selected_peak": 60}
        for column in range(self.references.table.columnCount()):
            key = self.references.table.horizontalHeaderItem(column).text()
            if key in widths:
                self.references.table.setColumnWidth(column, widths[key])

    @staticmethod
    def _display_references(frame):
        frame = frame.copy()
        if "ccs" in frame and "CCS_reported" not in frame:
            frame = frame.rename(columns={"ccs": "CCS_reported"})
        if "peak_select" not in frame:
            frame["peak_select"] = "highest"
        else:
            frame["peak_select"] = frame["peak_select"].fillna("highest").replace("", "highest")
        # Resolve legacy centre/half-width imports once, making the actual band editable.
        from deconvovo.imms_ccs_calibrate import resolve_mz_bounds
        for key in REFERENCE_COLUMNS:
            if key not in frame:
                frame[key] = ""
        for index, row in frame.iterrows():
            missing = lambda value: pd.isna(value) or not str(value).strip()
            if missing(row.get("mz_min", "")) and missing(row.get("mz_max", "")):
                try:
                    lower, upper = resolve_mz_bounds(row, float(row["mw"]), int(float(row["z"])))
                    frame.loc[index, "mz_min"] = str(lower)
                    frame.loc[index, "mz_max"] = str(upper)
                except (TypeError, ValueError, KeyError):
                    pass  # Incomplete new rows remain editable, and are checked on Process.
        order = [key for key in REFERENCE_COLUMNS if key in frame]
        return frame[order + [key for key in frame if key not in order]]

    @staticmethod
    def _reference_signature(frame):
        keys = [key for key in frame if key not in RESULT_COLUMNS | {"peak_select", "peak_drift_bin"}]
        return tuple(tuple(str(value) for value in row) for row in frame[keys].fillna("").to_numpy())

    def _peak_mode_edited(self, row, column):
        if self.references._loading:
            return
        headers = [self.references.table.horizontalHeaderItem(i).text()
                   for i in range(self.references.table.columnCount())]
        if headers[column] == "peak_select" and "peak_drift_bin" in headers:
            override = self.references.table.item(row, headers.index("peak_drift_bin"))
            if override is not None and override.text():
                override.setText("")

    def _reference_clicked(self, row, _column):
        self._active_reference_changed(row, 0, -1, -1)

    def _active_reference_changed(self, row, _column, _previous_row, _previous_column):
        if self.references._loading:
            return
        index = next((i for i, profile in enumerate(self._profiles) if profile["row_index"] == row), None)
        if index is not None:
            blocker = QSignalBlocker(self.combo_reference)
            self.combo_reference.setCurrentIndex(index)
            del blocker
            self._draw_profile()
        self._draw_curve()

    def _reference_combo_changed(self, index):
        if 0 <= index < len(self._profiles):
            row = self._profiles[index]["row_index"]
            self.references.table.setCurrentCell(row, 0)
            self.references.table.selectRow(row)
        self._draw_profile()
        self._draw_curve()

    def _selected_reference_rows(self):
        rows = sorted({index.row() for index in self.references.table.selectedIndexes()})
        current = self.references.table.currentRow()
        return rows or ([current] if current >= 0 else [])

    def _reference_selection_changed(self):
        if self.references._loading:
            return
        self._draw_curve()
        self._update_actions()

    def set_peak_selection(self, mode):
        """Apply a preset to the explicitly selected rows or the whole sheet."""
        if self._busy:
            return
        rows = (list(range(self.references.table.rowCount()))
                if self.combo_selection_scope.currentData() == "all" else self._selected_reference_rows())
        if not rows:
            return
        headers = [self.references.table.horizontalHeaderItem(i).text()
                   for i in range(self.references.table.columnCount())]
        changed = False
        self.references._loading = True
        blocker = QSignalBlocker(self.references.table)
        for row in rows:
            if "peak_drift_bin" in headers:
                override = self.references.table.item(row, headers.index("peak_drift_bin"))
                changed |= bool(override.text().strip())
                override.setText("")
            item = self.references.table.item(row, headers.index("peak_select"))
            current = {"tallest": "highest", "latest": "last"}.get(item.text().lower(), item.text().lower())
            changed |= current != mode
            item.setText(mode)
        del blocker
        self.references._loading = False
        if changed:
            self._references_changed()

    def _summary_clicked(self, row, _column):
        if 0 <= row < len(self._profiles):
            self.combo_reference.setCurrentIndex(row)

    def _refresh_selections(self):
        from deconvovo.imms_ccs_calibrate import resolve_calibrant_peak
        frame = self.references.to_dataframe()
        for profile in self._profiles:
            row = frame.iloc[profile["row_index"]]
            profile["mode"] = str(row.get("peak_select", "highest")).lower()
            profile["selection_error"] = False
            try:
                selected = resolve_calibrant_peak(profile["candidates"], row, allow_unselected_manual=True) if profile["candidates"] else None
            except ValueError as error:
                selected = None
                profile["selection_error"] = True
                self.status.setText(str(error))
                self._set_status_color("error")
            profile["selected_bin"] = selected["bin_index"] if selected else None
        self._write_selected_values()
        self._refresh_summary()

    def _write_selected_values(self):
        self.references._loading = True
        headers = [self.references.table.horizontalHeaderItem(i).text()
                   for i in range(self.references.table.columnCount())]
        for profile in self._profiles:
            chosen = next(((i, peak) for i, peak in enumerate(profile["candidates"], 1)
                           if peak["bin_index"] == profile.get("selected_bin")), None)
            values = {"selected_peak": str(chosen[0]) if chosen else "",
                      "selected_td_ms": str(chosen[1]["drift_time_ms"]) if chosen else ""}
            for key, value in values.items():
                item = self.references.table.item(profile["row_index"], headers.index(key))
                item.setText(value)
                item.setToolTip(value + (" ms; acquired drift bin " + str(profile["selected_bin"]) if key == "selected_td_ms" and chosen else ""))
        self.references._loading = False

    def _refresh_summary(self):
        rows = []
        references = self.references.to_dataframe()
        calibration_rows = (self._calibration or {}).get("calibrant_summary", [])
        error_column = "err_twostep_pct" if (self._calibration or {}).get("method") == "twostep" else "err_direct_pct"
        for profile in self._profiles:
            chosen = next(((i, peak) for i, peak in enumerate(profile["candidates"], 1)
                           if peak["bin_index"] == profile.get("selected_bin")), None)
            item = {"name": profile["name"], "mz_min": profile["mz_min"], "mz_max": profile["mz_max"],
                    "peak_select": profile["mode"], "selected_peak": chosen[0] if chosen else "",
                    "selected_td_ms": chosen[1]["drift_time_ms"] if chosen else "",
                    "candidates": len(profile["candidates"]),
                    "CCS_reported": references.iloc[profile["row_index"]]["CCS_reported"]}
            fitted = next((row for row in calibration_rows if row["name"] == profile["name"]
                           and row["z"] == profile["z"] and math.isclose(row["mz"], profile["mz"])), {})
            if calibration_rows:
                item[error_column] = fitted.get(error_column, "")
            if profile.get("reason"):
                item["status"] = profile["reason"]
            rows.append(item)
        self.summary.load_dataframe(pd.DataFrame(rows))
        self._format_summary_columns()

    def _format_summary_columns(self):
        for column, key in enumerate(self.summary.to_dataframe().columns):
            self.summary.table.setItemDelegateForColumn(column, None)
            if key in {"selected_td_ms", "drift_time_ms"}:
                self.summary.table.setItemDelegateForColumn(column, self._drift_value_delegate)
            elif key == "peak_select":
                self.summary.table.setItemDelegateForColumn(column, self._peak_choice_delegate)

    def load_reference_csv(self, path=None):
        if not isinstance(path, (str, Path)):
            path, _ = QFileDialog.getOpenFileName(self, "Load reference CSV", "", "CSV files (*.csv)")
        if not path:
            return
        try:
            self.set_references(pd.read_csv(path), path)
        except (OSError, ValueError) as error:
            self._error(str(error))

    def add_reference_raw(self, path=None):
        if not isinstance(path, (str, Path)):
            path = QFileDialog.getExistingDirectory(self, "Choose a standard RAW acquisition")
        if not path:
            return
        raw = local_path(path).resolve()
        if not raw.is_dir() or raw.suffix.lower() != ".raw":
            self._error("Choose a complete .raw acquisition folder.")
            return
        self._reference_raw = raw
        self.reference_source.setText(raw.name)
        self.reference_source.setToolTip(str(raw))
        self.add_reference_ion()

    def add_reference_ion(self):
        frame = self.references.to_dataframe()
        widths = {str(frame.columns[column]): self.references.table.columnWidth(column)
                  for column in range(len(frame.columns))}
        source = self._reference_raw
        if source is None and not frame.empty and "raw_path" in frame:
            selected = self.references.table.currentRow()
            text = str(frame.iloc[selected if selected >= 0 else -1]["raw_path"])
            if text.strip():
                source = local_path(text)
        row = {key: "" for key in frame.columns}
        row.update(name=f"Ion {len(frame) + 1}", raw_path=str(source) if source else "", peak_select="highest")
        self.references.load_dataframe(pd.concat([frame, pd.DataFrame([row])], ignore_index=True),
                                       self.references._csv_path)
        self._size_reference_columns()
        for column, key in enumerate(self.references.to_dataframe().columns):
            if str(key) in widths:
                self.references.table.setColumnWidth(column, widths[str(key)])
        self.references.table.setCurrentCell(len(frame), 1)
        self._references_changed()

    def _start(self, function, done, *args):
        self._set_status_color("text_dim")
        self._busy = True
        self._operation = function.__name__
        self.busy_changed.emit(True)
        self._worker = Worker(function, *args)
        self._worker.finished.connect(done)
        self._worker.error.connect(self._error)
        self.progress.show()
        self._worker.start()
        self._update_actions()

    def _check_output(self):
        if any(parent.suffix.lower() == ".raw" and parent.is_dir()
               for parent in (self._out_dir, *self._out_dir.parents)):
            raise ValueError("Choose an output directory outside the RAW acquisition.")

    def _reference_snapshot(self):
        self._check_output()
        self.references.validated_dataframe(REFERENCE_REQUIRED, "Reference")
        snapshots = self._out_dir / "reference_inputs"
        snapshots.mkdir(parents=True, exist_ok=True)
        csv_path = self.references.write_snapshot(snapshots / "references.csv", REFERENCE_REQUIRED, "Reference")
        from deconvovo.imms_ccs_calibrate import normalize_calibrant_columns, normalize_calibrant_mz
        frame = normalize_calibrant_mz(normalize_calibrant_columns(pd.read_csv(csv_path)))
        frame.to_csv(csv_path, index=False)
        return frame, csv_path

    def process_references(self):
        try:
            frame, _csv_path = self._reference_snapshot()
            self._processing_signature = self._reference_signature(self.references.to_dataframe())
            self._calibration = None
            self._preview_calibration = None
            self._generated = False
            self.calibration_changed.emit(None)
            self.status.setText("Processing references… extracting drift profiles and detecting peaks.")
            self._start(self._process, self._processed, frame, self._out_dir, self._data_dir)
        except (OSError, ValueError) as error:
            self._error(str(error))

    def build_calibration(self):
        try:
            if not self._profiles or not self._processed_directories:
                raise ValueError("Process references first, then review or select their drift peaks.")
            if any(profile.get("mode") == "manual" and profile.get("selected_bin") is None for profile in self._profiles):
                raise ValueError("Choose a detected peak for every Manual reference before fitting.")
            frame, csv_path = self._reference_snapshot()
            if sum(profile.get("selected_bin") is not None for profile in self._profiles) < 3:
                raise ValueError("A calibration needs at least three selected reference peaks.")
            self.status.setText("Fitting calibration to the selected drift peaks…")
            self._start(self._fit, self._built, frame, csv_path, self._out_dir,
                        self._processed_directories, self.combo_method.currentData(), self._profiles)
        except (OSError, ValueError) as error:
            self._error(str(error))

    @staticmethod
    def _process(frame, out_dir, current_data):
        from deconvovo.waters_workflow import run as convert
        directories = []
        sources = [local_path(path).resolve() for path in frame["raw_path"].unique()]
        if len({source.stem.casefold() for source in sources}) != len(sources):
            raise ValueError("Standard acquisitions must have distinct file names.")
        for index, raw in enumerate(sources):
            if not raw.is_dir() or raw.suffix.lower() != ".raw":
                raise ValueError(f"Standard acquisition does not exist: {raw}")
            # The same file name can belong to a different acquisition. Reuse
            # the current pair only when its conversion provenance matches.
            reusable = False
            if current_data:
                settings = current_data / "conversion_settings.json"
                try:
                    previous = json.loads(settings.read_text(encoding="utf-8")).get(raw.stem, {})
                    signature = [[p.name, p.stat().st_size, p.stat().st_mtime_ns]
                                 for p in sorted(raw.iterdir()) if p.is_file()]
                    reusable = (previous.get("raw_source") == str(raw)
                                and previous.get("source_files") == signature
                                and all((current_data / f"{raw.stem}_{kind}.txt").is_file()
                                        for kind in ("ms", "im")))
                except (OSError, ValueError, TypeError):
                    pass
            if reusable:
                if current_data not in directories:
                    directories.append(current_data)
            else:
                result = convert(raw, out_dir / "standards" / str(index + 1),
                                 html=False, skip_existing=True)
                directories.append(Path(result["converted_dir"]))
        profiles = CcsWorkbench._load_profiles(frame, directories)
        return {"profiles": profiles, "directories": directories}

    @staticmethod
    def _fit(frame, csv_path, out_dir, directories, method="direct", profiles=None):
        from deconvovo.ccs_fit import build_calibration_curve
        from deconvovo.io import read_extern_inf
        from deconvovo.imms_ccs_calibrate import _acquisition_metadata, load_saved_calibration, save_calibration

        profiles = profiles if profiles is not None else CcsWorkbench._load_profiles(frame, directories)
        edc = None
        for raw in frame["raw_path"].unique():
            value = read_extern_inf(local_path(raw)).get("EDC Delay Coefficient")
            if value is None:
                continue
            value = float(value)
            if not math.isfinite(value):
                raise ValueError("EDC Delay Coefficient must be finite")
            if edc is not None and not math.isclose(edc, value, rel_tol=1e-6, abs_tol=1e-6):
                raise ValueError("Calibrant acquisitions have different EDC delay coefficients.")
            edc = value
        if edc is None:
            raise ValueError("Cannot read EDC Delay Coefficient from the standard RAW acquisitions.")
        gases = pd.to_numeric(frame.get("gas_mw", pd.Series(28.014, index=frame.index)), errors="raise").fillna(28.014).unique()
        if len(gases) != 1:
            raise ValueError("All calibrant rows must use the same gas_mw.")
        peaks, extraction = [], []
        for profile in profiles:
            selected = next((peak for peak in profile["candidates"]
                             if peak["bin_index"] == profile.get("selected_bin")), None)
            if selected is None:
                continue
            index = profile["row_index"]
            row = frame.loc[index]
            common = {"name": str(row["name"]), "mw": float(row["mw"]), "z": int(row["z"]),
                      "mz": profile["mz"], "ccs": float(row["ccs"]),
                      "mz_min": profile["mz_min"], "mz_max": profile["mz_max"],
                      "drift_time_ms": selected["drift_time_ms"]}
            peaks.append({**common, "drift_bin": selected["drift_bin"],
                          "peak_intensity": selected["peak_intensity"], "ms_intensity": profile["ms_intensity"]})
            extraction.append({**common, "row_index": int(index), "raw_path": str(row["raw_path"]),
                               "mz_window": (profile["mz_max"] - profile["mz_min"]) / 2,
                               "peak_select": profile["mode"], "peak_drift_bin": selected["bin_index"],
                               "pusher_period_us": profile["pusher_period_us"],
                               **{key: profile[key] for key in ("min_ms_intensity", "min_im_intensity", "peak_height_frac")},
                               "candidates": profile["candidates"]})
        cal = build_calibration_curve(peaks, edc, float(gases[0]))
        cal.update(format_version=1, method=method, accepted=False, extraction=extraction,
                   reference_acquisitions=[_acquisition_metadata(raw, directories) for raw in frame["raw_path"].unique()])
        cal = load_saved_calibration(cal)
        destination = out_dir / "calibration"
        destination.mkdir(parents=True, exist_ok=True)
        save_calibration(cal, destination / "calibration_curve.json", accepted=False)
        pd.DataFrame(cal["calibrant_summary"]).to_csv(destination / "calibrant_summary.csv", index=False)
        return {"calibration": cal, "profiles": profiles}

    @staticmethod
    def _build(frame, csv_path, out_dir, current_data, method="direct"):
        processed = CcsWorkbench._process(frame, out_dir, current_data)
        return CcsWorkbench._fit(frame, csv_path, out_dir, processed["directories"], method)

    @staticmethod
    def _load_profiles(frame, directories, calibration=None):
        from deconvovo.io import find_pusher_period, read_im_txt
        from deconvovo.imms_ccs_calibrate import resolve_mz_bounds, resolve_calibrant_peak
        from deconvovo.ccs_peak_pick import detect_peaks_in_profile, extract_drift_profile
        from deconvovo.io import read_ms_txt
        profiles = []
        for raw_text, group in frame.groupby("raw_path", sort=False):
            raw = local_path(raw_text)
            pusher = find_pusher_period(raw)
            im_file = next((directory / f"{raw.stem}_im.txt" for directory in directories
                            if (directory / f"{raw.stem}_im.txt").is_file()), None)
            if im_file is None or pusher is None:
                raise ValueError(f"Cannot read drift data or pusher period for {raw.name}.")
            ms_file = next((directory / f"{raw.stem}_ms.txt" for directory in directories
                            if (directory / f"{raw.stem}_ms.txt").is_file()), None)
            if ms_file is None:
                raise ValueError(f"Cannot find the paired mass spectrum for {raw.name}.")
            data = read_im_txt(im_file)
            ms_mz, ms_values = read_ms_txt(ms_file)
            for index, row in group.iterrows():
                row = row.copy()
                lower, upper = resolve_mz_bounds(row, float(row["mw"]), int(row["z"]))
                mz, width = (lower + upper) / 2, (upper - lower) / 2
                intensity = extract_drift_profile(data, mz, width)
                fraction = row.get("peak_height_frac", 0.1)
                fraction = 0.1 if pd.isna(fraction) else float(fraction)
                candidates = detect_peaks_in_profile(intensity, pusher, fraction)
                minimum = row.get("min_im_intensity", 500)
                minimum = 500 if pd.isna(minimum) else float(minimum)
                candidates = [peak for peak in candidates if peak["peak_intensity"] >= minimum]
                min_ms = row.get("min_ms_intensity", 100)
                min_ms = 100 if pd.isna(min_ms) else float(min_ms)
                ms_intensity = float(ms_values[(ms_mz >= lower) & (ms_mz <= upper)].max(initial=0))
                reason = ""
                if ms_intensity < min_ms:
                    candidates, reason = [], "MS intensity below threshold"
                elif not candidates:
                    reason = "No drift peaks pass the thresholds"
                selected = resolve_calibrant_peak(candidates, row, allow_unselected_manual=True) if candidates else None
                extraction = next((item for item in (calibration or {}).get("extraction", [])
                                   if item.get("row_index") == int(index)), {})
                profiles.append({"row_index": int(index), "name": str(row["name"]),
                                 "mz": mz, "mz_min": lower, "mz_max": upper, "z": int(row["z"]),
                                 "mode": str(row.get("peak_select", "highest")).lower(), "reason": reason,
                                 "ms_intensity": ms_intensity, "pusher_period_us": pusher,
                                 "min_ms_intensity": min_ms, "min_im_intensity": minimum, "peak_height_frac": fraction,
                                 "drift": (np.arange(len(intensity)) + 0.5) * pusher / 1000,
                                 "intensity": intensity, "candidates": candidates,
                                 "selected_bin": extraction.get("peak_drift_bin", selected["bin_index"] if selected else None)})
        return profiles

    def _processed(self, result):
        self._busy = False
        self.busy_changed.emit(False)
        self.progress.hide()
        self._profiles = result["profiles"]
        self._processed_directories = result["directories"]
        self._processed_signature = self._processing_signature
        self._populate_references()
        self._write_selected_values()
        self._refresh_summary()
        if self._profiles and self.references.table.currentRow() < 0:
            self.references.table.selectRow(self._profiles[0]["row_index"])
            self.references.table.setCurrentCell(self._profiles[0]["row_index"], 0)
        self._draw_curve()
        self._draw_profile()
        pending = sum(profile["mode"] == "manual" and profile["selected_bin"] is None for profile in self._profiles)
        self.status.setText(f"{len(self._profiles)} references processed" +
                            (f" · Choose {pending} Manual peak{'s' if pending != 1 else ''}." if pending else " · Review peaks, then fit calibration."))
        self._set_status_color("text_dim")
        self._update_actions()
        self._pulse_status()

    def _populate_references(self):
        previous = self.combo_reference.currentIndex()
        self.combo_reference.blockSignals(True)
        self.combo_reference.clear()
        for item in self._profiles:
            self.combo_reference.addItem(f"{item['name']} · {item['mz_min']:.1f}–{item['mz_max']:.1f} m/z")
        if self._profiles:
            self.combo_reference.setCurrentIndex(max(0, min(previous, len(self._profiles) - 1)))
        self.combo_reference.blockSignals(False)

    def _built(self, result):
        self._busy = False
        self.busy_changed.emit(False)
        self._profiles = result["profiles"]
        self._set_calibration(result["calibration"])
        self._generated = False
        self.progress.hide()
        self._update_actions()
        self._pulse_status()

    def _pulse_status(self):
        effects = getattr(self.window(), "_console_effects", None)
        if effects is not None:
            effects.pulse(self.status, warning=bool((self._calibration or {}).get("review_warnings")))

    def _set_calibration(self, calibration):
        self._calibration = calibration
        self._generated = False
        self._last_target_request = None
        self._preview_calibration = None
        self.combo_method.blockSignals(True)
        self.combo_method.setCurrentIndex(max(0, self.combo_method.findData(calibration.get("method", "direct"))))
        self.combo_method.blockSignals(False)
        r2 = float(calibration.get("r2_lnln", 0))
        self.status.setText(f"{calibration.get('n_points', 0)} references · R² {r2:.4f}" +
                            (" · Fit needs review." if r2 <= 0.95 else " · Review peaks, then Generate CCS."))
        self._set_status_color("warning" if r2 <= 0.95 else "text_dim")
        summary = pd.DataFrame(calibration.get("calibrant_summary", []))
        error_column = "err_twostep_pct" if calibration.get("method") == "twostep" else "err_direct_pct"
        fields = [key for key in ("name", "mz", "z", "drift_time_ms", "CCS_literature", error_column) if key in summary]
        self.summary.load_dataframe((summary[fields] if fields else summary).rename(columns={"CCS_literature": "CCS_reported"}))
        self._format_summary_columns()
        self._draw_curve()
        self._populate_references()
        if self._profiles:
            self._write_selected_values()
            self._refresh_summary()
        self._draw_profile()
        self.calibration_changed.emit(calibration)
        self._queue_target_ccs()

    def _draw_curve(self):
        primary, accent, _text = self._plot_colors()
        ax = self.plot.ax(0, 0)
        ax.clear()
        ax.tick_params(labelsize=9)
        ax.xaxis.get_offset_text().set_fontsize(9)
        ax.yaxis.get_offset_text().set_fontsize(9)
        cal = self._calibration or self._preview_calibration
        self._curve_points = []
        pending = self._calibration is None and cal is not None
        rows = cal.get("calibrant_summary", []) if cal else []
        if not rows:
            ax.set_axis_off()
        if rows:
            two_step = cal.get("method") == "twostep"
            if self._profiles:
                # Profiles retain the sheet's row_index even when RAW grouping
                # changes their order or two references share a display name.
                references = self.references.to_dataframe()
                for profile in self._profiles:
                    peak = next((peak for peak in profile["candidates"]
                                 if peak["bin_index"] == profile.get("selected_bin")), None)
                    if peak is None:
                        continue
                    row_index = profile["row_index"]
                    row = references.iloc[row_index]
                    try:
                        mass, charge, ccs = float(row["mw"]), int(float(row["z"])), float(row["CCS_reported"])
                        corrected = peak["drift_time_ms"] - cal["edc_coeff"] * math.sqrt(profile["mz"]) / 1000
                        if corrected <= 0:
                            continue
                        reduced_mass = math.sqrt(1 / mass + 1 / cal["gas_mw"])
                        x = corrected ** cal["X"] * charge * reduced_mass if two_step else corrected
                        y = ccs if two_step else ccs / (charge * reduced_mass)
                        if math.isfinite(x) and math.isfinite(y):
                            self._curve_points.append((row_index, x, y))
                    except (TypeError, ValueError, KeyError, ZeroDivisionError):
                        continue
            else:
                # A loaded curve can be used without its original acquisitions.
                # Do not associate it with unrelated rows already in the editor.
                for row in rows:
                    x = float(row["t_double_prime" if two_step else "t_prime"])
                    y = float(row["CCS_literature"]) if two_step else float(row["CCS_literature"]) / (
                        float(row["z"]) * math.sqrt(1 / float(row["mw"]) + 1 / cal["gas_mw"]))
                    self._curve_points.append((None, x, y))
            if self._curve_points:
                xs = np.array([point[1] for point in self._curve_points])
                ys = np.array([point[2] for point in self._curve_points])
                ax.scatter(xs, ys, color=primary, s=24, zorder=3)
                selected_rows = set(self._selected_reference_rows())
                marked = [point for point in self._curve_points if point[0] in selected_rows]
                if marked:
                    ax.scatter([point[1] for point in marked], [point[2] for point in marked],
                               facecolors="none", edgecolors=accent, linewidths=1.7, s=125, zorder=5)
                line = np.linspace(xs.min(), xs.max(), 150)
                fit = (cal["slope_diagnostic"] * line + cal["intercept_diagnostic"] if two_step else
                       cal["A"] * line ** cal["X"])
                ax.plot(line, fit, color=accent, linestyle="--" if pending else "-")
            ax.set_xlabel("Transformed drift time" if two_step else "Corrected drift time (ms)", fontsize=10)
            ax.set_ylabel("CCS (Å²)" if two_step else "Reduced CCS", fontsize=10)
            ax.grid(alpha=0.2)
        title = ("Calibration curve · fit required" if pending else
                 f"Calibration · R² {cal['r2_lnln']:.4f}" if cal and "r2_lnln" in cal else "Calibration curve")
        ax.set_title(title, fontsize=11)
        self.plot.refresh()

    def _draw_profile(self, *_args):
        primary, accent, text = self._plot_colors()
        ax = self.plot.ax(0, 1)
        ax.clear()
        ax.tick_params(labelsize=9)
        ax.xaxis.get_offset_text().set_fontsize(9)
        ax.yaxis.get_offset_text().set_fontsize(9)
        self.combo_peak.clear()
        self.peak_status.setText("")
        index = self.combo_reference.currentIndex()
        if 0 <= index < len(self._profiles):
            profile = self._profiles[index]
            line, = ax.plot(profile["drift"], profile["intensity"], color=primary)
            line.get_path().should_simplify = False
            selected_index = -1
            for number, peak in enumerate(profile["candidates"], 1):
                bin_index = peak["bin_index"]
                self.combo_peak.addItem(f"Peak {number} · {peak['drift_time_ms']:.1f} ms", bin_index)
                self.combo_peak.setItemData(number - 1, f"t_d {peak['drift_time_ms']:.6g} ms · bin {bin_index} · intensity {peak['peak_intensity']:.6g}", Qt.ToolTipRole)
                selected = bin_index == profile.get("selected_bin")
                ax.scatter([peak["drift_time_ms"]], [peak["peak_intensity"]],
                           color=accent if selected else primary, s=40 if selected else 20,
                           marker="*" if selected else "o", zorder=4)
                ax.annotate(str(number), (peak["drift_time_ms"], peak["peak_intensity"]),
                            xytext=(0, 7), textcoords="offset points", ha="center", fontsize=8, color=text)
                if bin_index == profile.get("selected_bin"):
                    selected_index = number - 1
                    ax.axvline(peak["drift_time_ms"], color=accent, linestyle="--")
            ax.margins(y=0.16)
            self.combo_peak.setCurrentIndex(selected_index)
            self.combo_peak.setPlaceholderText("Choose peak")
            mode = {"tallest": "Highest", "latest": "Last"}.get(profile["mode"], profile["mode"].capitalize())
            chosen = profile["candidates"][selected_index] if selected_index >= 0 else None
            self.peak_status.setText(f"{mode} · m/z {profile['mz_min']:.1f}–{profile['mz_max']:.1f} · " +
                                     (f"Peak {selected_index + 1} · t<sub>d</sub> {chosen['drift_time_ms']:.1f} ms"
                                      if chosen else profile.get("reason") or "Click a detected peak to select its t<sub>d</sub>."))
            self.peak_status.setToolTip(f"Selected acquired bin {chosen['bin_index']}; t_d = {chosen['drift_time_ms']:.12g} ms. Click another detected peak to select it manually." if chosen else "Click a detected drift peak to select it manually.")
            ax.set_xlabel("t$_d$ (ms)", fontsize=10)
            ax.set_ylabel("Intensity", fontsize=10)
            ax.set_title(profile["name"], fontsize=11)
            ax.grid(alpha=0.2)
        else:
            ax.set_axis_off()
            message = ("Saved curve loaded.\nProcess references to\ninspect their peaks."
                       if self._calibration else "Process references\nto inspect their drift peaks.")
            # Titles already follow PlotCanvas's dark/light theme on refresh.
            ax.set_title(message, y=0.5, va="center", fontsize=11)
        self.plot.refresh()
        self._update_actions()

    def _profile_clicked(self, event):
        if self._busy or event.button != 1:
            return
        if event.inaxes is self.plot.ax(0, 0):
            points = [point for point in self._curve_points if point[0] is not None]
            if points:
                ax = self.plot.ax(0, 0)
                distances = [np.hypot(*(ax.transData.transform((x, y)) - np.array([event.x, event.y])))
                             for _row, x, y in points]
                closest = int(np.argmin(distances))
                if distances[closest] <= 24:
                    row = points[closest][0]
                    self.references.table.setCurrentCell(row, 0)
                    self.references.table.selectRow(row)
            return
        if event.inaxes is not self.plot.ax(0, 1):
            return
        index = self.combo_reference.currentIndex()
        if not 0 <= index < len(self._profiles):
            return
        candidates = self._profiles[index]["candidates"]
        if not candidates:
            return
        ax = self.plot.ax(0, 1)
        distances = [np.hypot(*(ax.transData.transform((peak["drift_time_ms"], peak["peak_intensity"])) -
                                 np.array([event.x, event.y]))) for peak in candidates]
        closest = int(np.argmin(distances))
        if distances[closest] <= 24:
            self.combo_peak.setCurrentIndex(closest)
            self.use_selected_peak()

    def use_selected_peak(self, *_args):
        index = self.combo_reference.currentIndex()
        if index < 0 or self.combo_peak.currentIndex() < 0:
            return
        frame = self.references.to_dataframe()
        row = self._profiles[index]["row_index"]
        if row not in frame.index:
            return
        if "peak_drift_bin" not in frame:
            frame["peak_drift_bin"] = ""
        frame.loc[row, "peak_select"] = "manual"
        frame.loc[row, "peak_drift_bin"] = self.combo_peak.currentData()
        selected_rows = self._selected_reference_rows()
        widths = [self.references.table.columnWidth(i) for i in range(self.references.table.columnCount())]
        self.references.load_dataframe(frame, self.references._csv_path)
        self._size_reference_columns()
        self._show_reference_columns()
        for column, width in enumerate(widths):
            self.references.table.setColumnWidth(column, width)
        self.references.table.setCurrentCell(row, 0)
        selection = self.references.table.selectionModel()
        for selected_row in selected_rows:
            selection.select(self.references.table.model().index(selected_row, 0),
                             selection.SelectionFlag.Select | selection.SelectionFlag.Rows)
        self._references_changed()
        self.status.setText("Manual peak preview updated · Fit calibration, then Generate CCS.")

    def load_calibration(self, path=None):
        if not isinstance(path, (str, Path)):
            path, _ = QFileDialog.getOpenFileName(self, "Load calibration", "", "Calibration (*.json)")
        if not path:
            return
        try:
            from deconvovo.imms_ccs_calibrate import load_saved_calibration
            calibration = load_saved_calibration(path)
            self._profiles = []
            self._processed_signature = None
            self._processed_directories = []
            self._set_calibration(calibration)
            self._generated = False
            warnings = calibration.get("review_warnings", [])
            self.status.setText(f"Loaded curve · {calibration.get('n_points', 0)} references · "
                                f"R² {calibration.get('r2_lnln', 0):.4f}" +
                                (" · Fit needs review." if warnings else "."))
            self.status.setToolTip(f"{local_path(path).name}\nUse matching acquisition conditions." +
                                  ("\n" + "\n".join(warnings) if warnings else ""))
            self._set_status_color("warning" if warnings else "text_dim")
            self._update_actions()
        except (OSError, ValueError, KeyError) as error:
            self._error(str(error))

    def save_calibration(self, path=None):
        if self._calibration is None:
            self._error("Fit or load a calibration before saving it.")
            return
        if not isinstance(path, (str, Path)):
            path, _ = QFileDialog.getSaveFileName(self, "Save calibration", "calibration.json", "Calibration (*.json)")
        if not path:
            return
        try:
            from deconvovo.imms_ccs_calibrate import load_saved_calibration, save_calibration
            destination = save_calibration(self._calibration, path, accepted=self._generated)
            self._calibration = load_saved_calibration(destination)
            self.status.setText(f"Saved {destination.name}. Reuse it for acquisitions measured under matching conditions.")
            self._set_status_color("success")
        except (OSError, ValueError) as error:
            self._error(str(error))

    def generate_ccs(self):
        if self._calibration is None:
            return
        self._generated = True
        self._last_target_request = None
        self._refresh_target_ccs()
        self._update_actions()

    def apply_targets(self, targets=None):
        if self._calibration is None or not self._generated:
            self._error("Load an accepted calibration, or review the current peaks and fit first.")
            return
        if not self._record or not self._record.get("raw") or not self._data_dir or self._targets.empty:
            self._error("Open an acquisition and add analyte targets below its 2D spectrum.")
            return
        try:
            self._check_output()
            automatic = isinstance(targets, pd.DataFrame)
            targets = targets.copy() if automatic else self._targets.copy()
            from deconvovo.imms_ccs_calibrate import load_saved_calibration, save_calibration
            saved = save_calibration(self._calibration, self._out_dir / "accepted_calibration.json", accepted=True)
            self._calibration = load_saved_calibration(saved)
            self._last_target_request = self._target_request_key(targets)
            self.status.setText("Applying the saved calibration to current targets…")
            self._start(self._apply, self._applied, self._calibration,
                        self._out_dir / ("target_values" if automatic else "analytes"),
                        self._data_dir, targets, self._record.get("raw"), False)
        except (OSError, ValueError) as error:
            self._error(str(error))

    @staticmethod
    def _apply(calibration, out_dir, data_dir, targets, raw_path, plot_profiles=True):
        from deconvovo.imms_ccs_calibrate import apply_saved_calibration
        return apply_saved_calibration(calibration, out_dir, data_dir, targets, raw_path=raw_path,
                                       plot_profiles=plot_profiles)

    def _applied(self, result):
        self._busy = False
        self.busy_changed.emit(False)
        self.progress.hide()
        warnings = result.get("warnings", [])
        targets = result.get("analyte_summary", [])
        notes = []
        if (self._calibration or {}).get("review_warnings"):
            notes.append("fit needs review")
        extrapolated = sum(bool(target.get("extrapolated")) for target in targets)
        if extrapolated:
            notes.append(f"{extrapolated} extrapolated")
        if any(check.get("status") != "match" for check in result.get("condition_checks", [])):
            notes.append("conditions need review")
        if warnings and not notes:
            notes.append("review required")
        self.status.setText(f"CCS ready for {len(targets)} targets" +
                            (" · " + " · ".join(notes) if notes else "."))
        self.status.setToolTip("\n".join(warnings))
        for warning in warnings:
            self._log.log(warning, "warning")
        self._set_status_color("warning" if warnings else "success")
        self._log.log(self.status.text(), "warning" if warnings else "success")
        self.targets_applied.emit(result)
        self._update_actions()

    def _error(self, message):
        self._last_target_request = None
        self._busy = False
        self.busy_changed.emit(False)
        self.progress.hide()
        self.status.setText(str(message).splitlines()[0])
        self._set_status_color("error")
        self._log.log(str(message), "error")
        self._update_actions()

    def wait_for_worker(self):
        if self._worker and self._worker.isRunning():
            self._worker.wait()

    def closeEvent(self, event):
        self._target_refresh.stop()
        self.wait_for_worker()
        super().closeEvent(event)
