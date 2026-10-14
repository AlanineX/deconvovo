"""Native Waters instrument console backed by the existing scientific workflow."""
from __future__ import annotations

import json
import logging
import re
import shutil
import sys
import tempfile
from pathlib import Path

import pandas as pd
from PySide6.QtCore import QStandardPaths, Qt, QTimer, QUrl
from PySide6.QtGui import QDesktopServices, QFontMetrics
from PySide6.QtWidgets import (
    QAbstractItemView, QCheckBox, QComboBox, QDoubleSpinBox, QFileDialog,
    QFormLayout, QGridLayout, QGroupBox, QHBoxLayout, QHeaderView, QLabel,
    QProgressBar, QPushButton, QScrollArea, QSizePolicy, QSpinBox, QSplitter,
    QTabWidget, QTableWidget, QTableWidgetItem, QToolButton, QVBoxLayout, QWidget,
)

from gui.panels.ccs_analyte import CcsAnalytePanel
from gui.panels.ccs_calibrate import CcsCalibrationPanel
from gui.panels.html_viewer import HtmlViewerPanel
from gui.theme import ERROR, SUCCESS, TEXT_DIM, WARNING
from gui.widgets.file_picker import DirPicker
from gui.widgets.log_panel import LogPanel
from gui.widgets.waters_viewer import WatersViewerWidget
from gui.widgets.worker import Worker


def _bundle_root() -> Path:
    if getattr(sys, "frozen", False):
        return Path(sys._MEIPASS)
    return Path(__file__).resolve().parents[2]


def _demo_manifest() -> Path | None:
    root = _bundle_root()
    for relative in ("demo_data/production/demo_manifest.json", "examples/production/demo_manifest.json"):
        path = root / relative
        if path.is_file():
            return path
    return None


def discover_runs(source: Path | None, out_dir: Path | None = None) -> list[dict]:
    """Inspect run names and file presence without loading spectral arrays."""
    rows: dict[str, dict] = {}
    roots = []
    if source and source.exists():
        root = source.parent if source.is_file() else source
        roots.append(root)
        raw_dirs = [root] if root.suffix.casefold() == ".raw" else sorted(
            p for p in root.rglob("*") if p.is_dir() and p.suffix.casefold() == ".raw"
        )
        for raw in raw_dirs:
            name = raw.stem
            rows.setdefault(name, {"run_name": name})["raw"] = raw
    if out_dir and out_dir.is_dir():
        roots.append(out_dir / "_converted")
        roots.append(out_dir)
    for root in roots:
        if not root.is_dir() or root.suffix.casefold() == ".raw":
            continue
        for suffix, key in (("_ms.txt", "ms"), ("_im.txt", "im")):
            for file in sorted(root.glob(f"*{suffix}")):
                name = file.name[:-len(suffix)]
                if root != roots[0] and source and source.exists() and name not in rows:
                    continue
                record = rows.setdefault(name, {"run_name": name})
                if file.stat().st_size:
                    record[key] = file
    if out_dir and out_dir.is_dir():
        for file in sorted(out_dir.glob("*_2d_imms.html")):
            if file.stat().st_size:
                name = file.name[:-len("_2d_imms.html")]
                if source and source.exists() and name not in rows:
                    continue
                rows.setdefault(name, {"run_name": name})["viewer"] = file
    return [rows[name] for name in sorted(rows, key=str.casefold)]


class WatersWorkspacePanel(QWidget):
    """Acquisition browser, measured viewer, and settings in one console."""

    def __init__(self, log_panel: LogPanel, parent=None):
        super().__init__(parent)
        self.setObjectName("watersConsole")
        self._log = log_panel
        self._worker: Worker | None = None
        self._busy = False
        self._job_files: tempfile.TemporaryDirectory | None = None
        self._rows: list[dict] = []
        self._last_result: dict = {}
        self._auto_output = ""
        self._manifest_path = _demo_manifest()
        self._selected_run_title = "No acquisition selected"
        self._refresh_view_timer = QTimer(self)
        self._refresh_view_timer.setSingleShot(True)
        self._refresh_view_timer.setInterval(180)
        self._refresh_view_timer.timeout.connect(self._apply_viewer_settings)
        self.viewer_panel = HtmlViewerPanel(log_panel=self._log, parent=self)
        self.calibration_panel = CcsCalibrationPanel(log_panel=self._log, parent=self)
        self.analyte_panel = CcsAnalytePanel(log_panel=self._log, parent=self)
        self.viewer_panel.hide()

        layout = QVBoxLayout(self)
        layout.setContentsMargins(20, 16, 20, 12)
        layout.setSpacing(12)
        header = QHBoxLayout()
        heading = QVBoxLayout()
        eyebrow = QLabel("INSTRUMENT CONSOLE  /  WATERS")
        eyebrow.setProperty("eyebrow", True)
        heading.addWidget(eyebrow)
        title = QLabel("Waters IM-MS workbench")
        title.setProperty("title", True)
        heading.addWidget(title)
        subtitle = QLabel("One batch. Conversion, visual inspection, and CCS together.")
        subtitle.setProperty("subtitle", True)
        heading.addWidget(subtitle)
        header.addLayout(heading, 1)
        self.queue_summary = QLabel("No inputs loaded")
        self.queue_summary.setProperty("chip", True)
        header.addWidget(self.queue_summary, alignment=Qt.AlignVCenter)
        self.btn_run = QPushButton("Process && view")
        self.btn_run.setMinimumHeight(40)
        self.btn_run.setMinimumWidth(150)
        self.btn_run.clicked.connect(self.process_and_view)
        header.addWidget(self.btn_run)
        layout.addLayout(header)

        self.console_splitter = QSplitter(Qt.Horizontal)
        self.console_splitter.setChildrenCollapsible(False)
        self.console_splitter.setHandleWidth(8)
        self._build_acquisition_browser()
        self._build_plot_workspace()
        self._build_inspector()
        self.console_splitter.setStretchFactor(0, 0)
        self.console_splitter.setStretchFactor(1, 1)
        self.console_splitter.setStretchFactor(2, 0)
        self.console_splitter.setSizes([250, 700, 300])
        layout.addWidget(self.console_splitter, 1)
        self.job_status = QLabel("Choose an acquisition folder or load a real demo to begin.")
        self.job_status.setWordWrap(True)
        self.job_status.setProperty("subtitle", True)
        layout.addWidget(self.job_status)
        self.job_progress = QProgressBar()
        self.job_progress.setVisible(False)
        self.job_progress.setMaximumHeight(8)
        self.job_progress.setTextVisible(False)
        layout.addWidget(self.job_progress)
        self._load_demo_choices()
        self.pick_source.path_changed.connect(self._source_changed)
        self.pick_out.path_changed.connect(self._paths_changed)
        self.calibration_panel.pick_csv.path_changed.connect(self._sync_ccs_paths)
        self.combo_method.currentTextChanged.connect(self._paths_changed)
        self.calibration_panel.table.data_changed.connect(self._update_ccs_cards)
        self.analyte_panel.table.data_changed.connect(self._update_ccs_cards)
        self.calibration_panel.pick_csv.path_changed.connect(self._update_ccs_cards)
        self.analyte_panel.pick_ana.path_changed.connect(self._update_ccs_cards)
        controls = [control for kind in (QSpinBox, QDoubleSpinBox, QComboBox, QCheckBox)
                    for control in self._settings_group.findChildren(kind)] + self.display_controls
        for control in controls:
            signal = control.currentIndexChanged if isinstance(control, QComboBox) else control.toggled if isinstance(control, QCheckBox) else control.valueChanged
            signal.connect(lambda *_: self._refresh_view_timer.start())
        self.viewer_panel.presets_editor.textChanged.connect(self._refresh_view_timer.start)
        self.reload_runs()

    @staticmethod
    def _button(text, action):
        button = QPushButton(text)
        button.setProperty("secondary", True)
        button.clicked.connect(action)
        return button

    @staticmethod
    def _scroll(content):
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QScrollArea.NoFrame)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        scroll.setWidget(content)
        return scroll

    def _compact_picker(self, picker, caption, layout):
        label = QLabel(caption)
        label.setProperty("subtitle", True)
        label.setProperty("controlLabel", True)
        layout.addWidget(label)
        picker.layout().itemAt(0).layout().itemAt(0).widget().hide()
        button = picker.findChild(QPushButton)
        button.setText("…")
        button.setFixedWidth(32)
        button.setToolTip(picker._title)
        picker.line.setMinimumWidth(0)
        picker.line.setToolTip(caption)
        layout.addWidget(picker)

    def _build_acquisition_browser(self):
        self.browser_panel = QWidget()
        self.browser_panel.setMinimumWidth(220)
        self.browser_panel.setMaximumWidth(330)
        column = QVBoxLayout(self.browser_panel)
        column.setContentsMargins(0, 0, 0, 0)
        column.setSpacing(12)
        self.sources_group = QGroupBox("BATCH SOURCE")
        paths = QVBoxLayout(self.sources_group)
        paths.setContentsMargins(0, 0, 0, 0)
        paths.setSpacing(6)
        self.combo_demo = QComboBox()
        self.combo_demo.setMinimumWidth(0)
        self.combo_demo.setSizeAdjustPolicy(QComboBox.AdjustToMinimumContentsLengthWithIcon)
        self.combo_demo.setMinimumContentsLength(12)
        demo = QHBoxLayout()
        demo.addWidget(self.combo_demo, 1)
        self.btn_demo = self._button("Load", self._load_demo)
        self.btn_demo.setToolTip("Load complete real acquisitions and their matching settings")
        self.btn_demo.setProperty("compact", True)
        self.btn_demo.setStyleSheet("padding: 6px;")
        self.btn_demo.setFixedWidth(72)
        demo.addWidget(self.btn_demo)
        paths.addLayout(demo)
        self.pick_source = DirPicker("Source:", dialog_title="Select a complete Waters .raw folder or parent")
        self._compact_picker(self.pick_source, "Source folder", paths)
        self.pick_source.line.setToolTip("Complete .raw folder, a parent containing runs, or a folder of paired converted TXT files")
        self.btn_text = self._button("Use converted TXT…", self._pick_text)
        paths.addWidget(self.btn_text)
        self.pick_out = DirPicker("Output:", dialog_title="Select a writable output folder")
        self._compact_picker(self.pick_out, "Output folder", paths)
        column.addWidget(self.sources_group)

        runs = QGroupBox("Acquisitions")
        runs_layout = QVBoxLayout(runs)
        runs_layout.setContentsMargins(0, 0, 0, 0)
        self.overview = QLabel()
        self.overview.setWordWrap(True)
        self.overview.setProperty("subtitle", True)
        runs_layout.addWidget(self.overview)
        self.run_table = QTableWidget(0, 5)
        self.run_table.setHorizontalHeaderLabels(["Acquisition", "Input", "Data", "Viewer", "State"])
        for col in (1, 2, 3):
            self.run_table.hideColumn(col)
        self.run_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.run_table.setSelectionMode(QAbstractItemView.SingleSelection)
        self.run_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.run_table.setWordWrap(False)
        self.run_table.setShowGrid(False)
        self.run_table.verticalHeader().setVisible(False)
        self.run_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
        self.run_table.horizontalHeader().setSectionResizeMode(4, QHeaderView.Fixed)
        self.run_table.setColumnWidth(4, 80)
        self.run_table.verticalHeader().setDefaultSectionSize(55)
        self.run_table.setMinimumHeight(110)
        self.run_table.itemSelectionChanged.connect(self._selected_run_changed)
        self.run_table.itemDoubleClicked.connect(lambda _item: self.open_viewer())
        runs_layout.addWidget(self.run_table, 1)
        self.btn_reload = self._button("Reload existing results", self.reload_runs)
        runs_layout.addWidget(self.btn_reload)
        column.addWidget(runs, 1)
        self.console_splitter.addWidget(self.browser_panel)

    def _build_plot_workspace(self):
        self.plot_panel = QWidget()
        self.plot_panel.setMinimumWidth(360)
        column = QVBoxLayout(self.plot_panel)
        column.setContentsMargins(0, 0, 0, 0)
        column.setSpacing(12)
        plot_card = QGroupBox()
        plot_card.setObjectName("watersPlotCard")
        plot_card.setStyleSheet("QGroupBox#watersPlotCard { margin-top: 0; padding: 12px; }")
        plot_layout = QVBoxLayout(plot_card)
        plot_layout.setContentsMargins(0, 0, 0, 0)
        selected = QHBoxLayout()
        labels = QVBoxLayout()
        heading = QLabel("SELECTED ACQUISITION")
        heading.setProperty("eyebrow", True)
        labels.addWidget(heading)
        self.selected_run_label = QLabel(self._selected_run_title)
        self.selected_run_label.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Preferred)
        self.selected_run_label.setProperty("section", True)
        labels.addWidget(self.selected_run_label)
        self.selected_state_label = QLabel("Select a run in the acquisition browser")
        self.selected_state_label.setProperty("subtitle", True)
        labels.addWidget(self.selected_state_label)
        selected.addLayout(labels, 1)
        self.btn_viewer = self._button("Open in browser", self.open_viewer)
        selected.addWidget(self.btn_viewer)
        plot_layout.addLayout(selected)
        self.tabs = QTabWidget()
        self.viewer_widget = WatersViewerWidget(self)
        self.viewer_widget.btn_browser.hide()  # the selected-run header owns this action
        self.tabs.addTab(self.viewer_widget, "Local viewer · offline")
        self._build_ccs_page()
        plot_layout.addWidget(self.tabs, 1)
        column.addWidget(plot_card, 1)

        output = QGroupBox("Batch outputs")
        output_layout = QVBoxLayout(output)
        output_layout.setContentsMargins(0, 0, 0, 0)
        row = QHBoxLayout()
        self.output_summary = QLabel("No outputs yet")
        self.output_summary.setWordWrap(False)
        self.output_summary.setToolTip("MS / IM tables and offline viewers appear after processing. CCS summaries appear after optional calibration.")
        self.output_summary.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Preferred)
        self.output_summary.setProperty("subtitle", True)
        row.addWidget(self.output_summary, 1)
        self.btn_folder = self._button("Open output folder", self.open_output_folder)
        row.addWidget(self.btn_folder)
        self.btn_export = self._button("Export portable HTML", self.export_portable_html)
        self.btn_export.setToolTip("Export the selected converted run with the current display defaults as one offline HTML file, including all interactive controls. No RAW reconversion is needed; temporary zoom changes are not saved.")
        for button in (self.btn_folder, self.btn_export):
            button.setProperty("compact", True)
        row.addWidget(self.btn_export)
        output_layout.addLayout(row)
        column.addWidget(output)
        self.console_splitter.addWidget(self.plot_panel)

    def _build_inspector(self):
        content = QWidget()
        column = QVBoxLayout(content)
        column.setContentsMargins(0, 0, 0, 0)
        column.setSpacing(12)
        self.display_group = QGroupBox("Viewer defaults")
        display = QVBoxLayout(self.display_group)
        display.setContentsMargins(0, 0, 0, 0)
        note = QLabel("Controls update the displayed viewer. Process & view saves them in fresh offline HTML.")
        note.setWordWrap(True)
        note.setProperty("subtitle", True)
        display.addWidget(note)
        self._settings_group = next(group for group in self.viewer_panel.findChildren(QGroupBox)
                                    if group.title() == "Viewer Configuration")
        grid = QGridLayout()
        grid.setHorizontalSpacing(10)
        grid.setVerticalSpacing(7)
        controls = [
            ("m/z bins", self.viewer_panel.spin_mzbins),
            ("Drift bins", self.viewer_panel.spin_driftbins),
            ("MS display", self.viewer_panel.combo_ms_mode),
            ("Intensity", self.viewer_panel.combo_zscale),
            ("Colour map", self.viewer_panel.combo_cmap),
            ("2D smoothing", self.viewer_panel.combo_smooth2d),
            ("Noise cutoff (%)", self.viewer_panel.spin_noise2d),
            ("Cutoff fade (%)", self.viewer_panel.spin_denoise),
        ]
        self.display_controls = [control for _, control in controls]
        self.viewer_panel.spin_denoise.setSpecialValueText("Off")
        for i, (caption, control) in enumerate(controls):
            for form in self._settings_group.findChildren(QFormLayout):
                if form.getWidgetPosition(control)[0] >= 0:
                    removed = form.takeRow(control)
                    if removed.labelItem and removed.labelItem.widget():
                        removed.labelItem.widget().hide()
                    break
            label = QLabel(caption)
            label.setProperty("subtitle", True)
            label.setProperty("controlLabel", True)
            grid.addWidget(label, (i // 2) * 2, i % 2)
            control.setMinimumWidth(0)
            if isinstance(control, QComboBox):
                control.setSizeAdjustPolicy(QComboBox.AdjustToMinimumContentsLengthWithIcon)
                control.setMinimumContentsLength(6)
            grid.addWidget(control, (i // 2) * 2 + 1, i % 2)
        display.addLayout(grid)
        column.addWidget(self.display_group)

        self.ccs_card = QGroupBox("Calibration + analytes")
        ccs = QVBoxLayout(self.ccs_card)
        ccs.setContentsMargins(0, 0, 0, 0)
        cards = QHBoxLayout()
        self.btn_review_cal = self._button("Calibrants\nReview table →", lambda: self.show_ccs_review(0))
        self.btn_review_ana = self._button("Analytes\nReview table →", lambda: self.show_ccs_review(1))
        for button in (self.btn_review_cal, self.btn_review_ana):
            button.setMinimumHeight(74)
            cards.addWidget(button, 1)
        ccs.addLayout(cards)
        self.chk_ccs = QCheckBox("Include CCS when processing")
        ccs.addWidget(self.chk_ccs)
        self.demo_description = QLabel("Load calibrant and analyte tables in CCS review.")
        self.demo_description.setWordWrap(True)
        self.demo_description.setProperty("subtitle", True)
        ccs.addWidget(self.demo_description)
        column.addWidget(self.ccs_card)

        self.btn_advanced = QToolButton()
        self.btn_advanced.setText("Advanced processing & defaults")
        self.btn_advanced.setToolButtonStyle(Qt.ToolButtonTextBesideIcon)
        self.btn_advanced.setArrowType(Qt.RightArrow)
        self.btn_advanced.setCheckable(True)
        column.addWidget(self.btn_advanced)
        self.advanced = QWidget()
        advanced_layout = QVBoxLayout(self.advanced)
        advanced_layout.setContentsMargins(0, 0, 0, 0)
        extraction = QGroupBox("RAW extraction")
        form = QFormLayout(extraction)
        form.setRowWrapPolicy(QFormLayout.WrapAllRows)
        self.spin_func = QSpinBox()
        self.spin_func.setRange(1, 100)
        self.spin_func.setValue(1)
        self.spin_ms_bin = QDoubleSpinBox()
        self.spin_im_bin = QDoubleSpinBox()
        for control in (self.spin_ms_bin, self.spin_im_bin):
            control.setRange(0, 1000)
            control.setDecimals(4)
            control.setSpecialValueText("Native resolution")
        self.chk_reuse = QCheckBox("Reuse unchanged RAW conversions")
        self.chk_reuse.setChecked(True)
        self.chk_reuse.setToolTip("Only unchanged acquisitions and extraction settings are reused. Display settings always regenerate HTML.")
        form.addRow("Acquisition function:", self.spin_func)
        form.addRow("MS bin width:", self.spin_ms_bin)
        form.addRow("IM bin width:", self.spin_im_bin)
        form.addRow(self.chk_reuse)
        advanced_layout.addWidget(extraction)
        self._settings_group.setTitle("Additional viewer defaults")
        self.viewer_panel.layout().removeWidget(self._settings_group)
        for form in self._settings_group.findChildren(QFormLayout):
            form.setRowWrapPolicy(QFormLayout.WrapAllRows)
            form.setLabelAlignment(Qt.AlignLeft)
        controls = [control for kind in (QSpinBox, QDoubleSpinBox, QComboBox)
                    for control in self._settings_group.findChildren(kind)]
        for control in controls:
            control.setMinimumWidth(0)
            if isinstance(control, QComboBox):
                control.setSizeAdjustPolicy(QComboBox.AdjustToMinimumContentsLengthWithIcon)
                control.setMinimumContentsLength(7)
        self._settings_group.setMinimumHeight(350)
        advanced_layout.addWidget(self._settings_group)
        self.viewer_panel.chk_skip.hide()
        self.viewer_panel.spin_workers.hide()
        for label in self._settings_group.findChildren(QLabel):
            if label.text() == "Workers:":
                label.hide()
        for button in self._settings_group.findChildren(QPushButton):
            if button.text() == "Save Config":
                button.setText("Save defaults")
                button.clicked.disconnect()
                button.clicked.connect(self._save_defaults)
                self.btn_save_defaults = button
        self.advanced.hide()
        self.btn_advanced.toggled.connect(self._toggle_advanced)
        column.addWidget(self.advanced)
        column.addStretch()
        self.inspector_scroll = self._scroll(content)
        self.inspector_scroll.setMinimumWidth(285)
        self.inspector_scroll.setMaximumWidth(370)
        self.console_splitter.addWidget(self.inspector_scroll)

    def _build_ccs_page(self):
        content = QWidget()
        page = QVBoxLayout(content)
        page.setContentsMargins(8, 8, 8, 8)
        help_text = QLabel("Review tables and calibration together. Shared batch paths feed the analysis automatically; calibrants must match instrument conditions.")
        help_text.setWordWrap(True)
        page.addWidget(help_text)
        self.ccs_demo_description = QLabel("")
        self.ccs_demo_description.setWordWrap(True)
        self.ccs_demo_description.setProperty("subtitle", True)
        page.addWidget(self.ccs_demo_description)
        row = QHBoxLayout()
        row.addWidget(QLabel("Method:"))
        self.combo_method = QComboBox()
        self.combo_method.addItems(["direct", "twostep"])
        row.addWidget(self.combo_method)
        row.addStretch()
        self.btn_ccs = QPushButton("Calibrate && apply")
        self.btn_ccs.clicked.connect(self.calibrate_and_apply)
        row.addWidget(self.btn_ccs)
        page.addLayout(row)
        self.ccs_editors = QTabWidget()
        for panel in (self.calibration_panel, self.analyte_panel):
            for picker in (panel.pick_data, panel.pick_out):
                picker.hide()
            panel.combo_method.hide()
            panel.btn_run.hide()
            panel.progress.hide()
            for label in panel.findChildren(QLabel, options=Qt.FindDirectChildrenOnly):
                if label.property("title") or label.property("subtitle") or label.text() == "Method:":
                    label.hide()
            panel.setMinimumHeight(400)
            splitter = panel.findChild(QSplitter)
            if splitter:
                splitter.setOrientation(Qt.Vertical)
                splitter.setSizes([150, 250])
            for picker in panel.findChildren(DirPicker):
                picker.line.setMinimumWidth(0)
        self.analyte_panel.pick_cal.hide()
        self.ccs_editors.addTab(self.calibration_panel, "Calibrants && curve")
        self.ccs_editors.addTab(self.analyte_panel, "Analytes && results")
        drift_tabs = next(tab for tab in self.calibration_panel.findChildren(QTabWidget)
                          if tab.indexOf(self.calibration_panel.plot_drift) >= 0)
        drift_index = drift_tabs.indexOf(self.calibration_panel.plot_drift)
        drift_tabs.removeTab(drift_index)
        self.drift_scroll = self._scroll(self.calibration_panel.plot_drift)
        drift_tabs.insertTab(drift_index, self.drift_scroll, "Drift Profiles")
        page.addWidget(self.ccs_editors, 1)
        self.tabs.addTab(self._scroll(content), "CCS review")

    def show_ccs_review(self, index=0):
        self.tabs.setCurrentIndex(1)
        self.ccs_editors.setCurrentIndex(index)

    def _save_defaults(self):
        if self.viewer_panel._save_config():
            self._log.log("Defaults saved.", "success")

    def _update_ccs_cards(self, *_):
        self.btn_review_cal.setText(f"Calibrants\n{self.calibration_panel.table.table.rowCount()} references\nReview →")
        self.btn_review_ana.setText(f"Analytes\n{self.analyte_panel.table.table.rowCount()} targets\nReview →")

    def _selected_run_changed(self):
        self._update_open_buttons()
        index = self.run_table.currentRow()
        if not 0 <= index < len(self._rows):
            self._selected_run_title = "No acquisition selected"
            self.selected_state_label.setText("Select a run in the acquisition browser")
            self.viewer_widget.clear("Select a batch or acquisition folder, then Process & view to inspect measured data.")
        else:
            record = self._rows[index]
            self._selected_run_title = record["run_name"]
            self.selected_state_label.setText(self.run_table.item(index, 4).text())
            source = Path(self.pick_source.path()) if self.pick_source.path() else None
            try:
                config = self.viewer_panel.current_config()
            except ValueError:
                config = self.viewer_panel._cfg
            self.viewer_widget.set_run(record, config, source)
        self._update_selected_title()

    def export_portable_html(self, destination=None):
        if self._busy:
            return None
        index = self.run_table.currentRow()
        if not 0 <= index < len(self._rows) or not all(self._rows[index].get(key) for key in ("ms", "im")):
            self._problem("Process & view first to create paired MS / IM data for portable export.")
            return None
        record = dict(self._rows[index])
        output = Path(self.pick_out.path()) if self.pick_out.path() else Path(record["ms"]).parent
        proposed = output / f"{record['run_name']}_2d_imms.html"
        if not destination or isinstance(destination, bool):
            selected, _ = QFileDialog.getSaveFileName(self, "Export portable offline viewer", str(proposed), "HTML viewer (*.html *.htm)")
            if not selected:
                return None
            destination = Path(selected)
        destination = Path(destination).expanduser().resolve()
        if destination in {Path(record[key]).resolve() for key in ("ms", "im")}:
            self._problem("Choose a separate HTML file. Converted MS / IM source tables cannot be replaced by export.")
            return None
        raw = Path(record["raw"]).resolve() if record.get("raw") else None
        if any(parent.suffix.casefold() == ".raw" for parent in destination.parents) or (raw and destination.is_relative_to(raw)):
            self._problem("Choose an export location outside RAW acquisition folders.")
            return None
        if not destination.suffix:
            destination = destination.with_suffix(".html")
        elif destination.suffix.casefold() not in {".html", ".htm"}:
            self._problem("Portable export requires an .html or .htm filename. Choose a separate HTML file.")
            return None
        self._job_files = tempfile.TemporaryDirectory(prefix="deconvovo-waters-export-")
        try:
            config_path = Path(self._job_files.name) / "viewer_config.json"
            config_path.write_text(json.dumps(self.viewer_panel.current_config(), indent=2), encoding="utf-8")
        except (OSError, ValueError) as error:
            self._cleanup_job_files()
            self._problem(str(error))
            return None
        self._set_busy(True)
        self.job_status.setText(f"Exporting portable HTML: {record['run_name']}…")
        self.job_status.setStyleSheet(f"color: {TEXT_DIM};")
        self.job_progress.setVisible(True)
        self.job_progress.setRange(0, 0)
        source = Path(self.pick_source.path()) if self.pick_source.path() else None
        self._worker = Worker(self._do_export, record, config_path, destination, source)
        self._worker.finished.connect(self._on_export_done)
        self._worker.error.connect(self._on_export_error)
        self._worker.start()
        return destination

    @staticmethod
    def _do_export(record, config_path, destination, source):
        from deconvovo.imms_html import plot_im_data
        from deconvovo.io import find_pusher_for_run, find_pusher_period
        pusher = find_pusher_period(Path(record["raw"])) if record.get("raw") else None
        if pusher is None:
            pusher = find_pusher_for_run(Path(record["ms"]).parent, record["run_name"], source)
        generated_dir = config_path.parent / "html"
        generated_dir.mkdir()
        plot_im_data(Path(record["im"]), Path(record["ms"]), record["run_name"], generated_dir,
                     pusher_us=pusher, config_path=config_path)
        generated = generated_dir / f"{record['run_name']}_2d_imms.html"
        if not generated.is_file() or not generated.stat().st_size:
            raise RuntimeError("Portable HTML generation produced no viewer")
        # Replace only after the complete single-file document exists.
        with tempfile.NamedTemporaryFile(dir=destination.parent, suffix=".html", delete=False) as handle:
            staged = Path(handle.name)
        try:
            shutil.copyfile(generated, staged)
            staged.replace(destination)
        finally:
            staged.unlink(missing_ok=True)
        return {"path": str(destination), "run_name": record["run_name"]}

    def _on_export_done(self, result):
        self.wait_for_worker()
        self._cleanup_job_files()
        self._set_busy(False)
        self.job_progress.setVisible(False)
        # Keep pipeline results intact and display the newly exported artifact.
        for index, record in enumerate(self._rows):
            if record["run_name"] == result["run_name"]:
                record["viewer"] = Path(result["path"])
                record["display_status"] = "Viewer ready"
                self.run_table.item(index, 3).setText("Offline HTML")
                self.run_table.item(index, 4).setText("Viewer ready")
                break
        self._selected_run_changed()
        self.job_status.setText(f"Exported one-file offline viewer: {Path(result['path']).name}")
        self.job_status.setStyleSheet(f"color: {SUCCESS};")
        self._log.log(f"Exported portable HTML: {result['path']}", "success")

    def _on_export_error(self, message):
        self.wait_for_worker()
        self._cleanup_job_files()
        self._set_busy(False)
        self.job_progress.setVisible(False)
        self.job_status.setText("Portable export failed. Review the log; existing results are retained.")
        self.job_status.setStyleSheet(f"color: {ERROR};")
        self._log.log(f"Portable export failed: {message}", "error")

    def _update_selected_title(self):
        self.selected_run_label.setToolTip(self._selected_run_title)
        self.selected_run_label.setText(QFontMetrics(self.selected_run_label.font()).elidedText(self._selected_run_title, Qt.ElideRight, max(100, self.selected_run_label.width())))

    def _apply_viewer_settings(self):
        try:
            config = self.viewer_panel.current_config()
        except ValueError:
            return
        if hasattr(self.viewer_widget, "apply_config"):
            self.viewer_widget.apply_config(config)

    def resizeEvent(self, event):
        super().resizeEvent(event)
        self._update_selected_title()

    def closeEvent(self, event):
        if self._busy:
            event.ignore()
            return
        self._refresh_view_timer.stop()
        self.viewer_widget.close()
        self.viewer_widget.wait_for_workers()
        super().closeEvent(event)

    def _toggle_advanced(self, checked):
        self.advanced.setVisible(checked)
        self.btn_advanced.setArrowType(Qt.DownArrow if checked else Qt.RightArrow)

    def _load_demo_choices(self):
        self.combo_demo.addItem("Choose a real acquisition demo…", None)
        if self._manifest_path:
            try:
                manifest = json.loads(self._manifest_path.read_text(encoding="utf-8"))
                for profile in manifest.get("profiles", []):
                    if profile.get("bundled", True):
                        self.combo_demo.addItem(profile.get("name", profile.get("label", profile["id"])), profile)
                index = next((i for i in range(1, self.combo_demo.count()) if self.combo_demo.itemData(i).get("id") == manifest.get("default_profile")), 0)
                self.combo_demo.setCurrentIndex(index)
            except (OSError, ValueError, KeyError) as error:
                self._log.log(f"Could not load demo catalog: {error}", "warning")
        self.btn_demo.setEnabled(self.combo_demo.count() > 1)
        if self.combo_demo.count() == 1:
            self.combo_demo.setItemText(0, "No bundled demo catalog found")

    def _load_demo(self):
        profile = self.combo_demo.currentData()
        if not profile or not self._manifest_path:
            return
        base = self._manifest_path.parent
        source = base / profile.get("source", profile.get("input_dir", ""))
        if not source.is_dir():
            self._problem(f"Demo acquisition folder is missing: {source}")
            return
        config = profile.get("viewer_config", profile.get("plot_config"))
        try:
            if config:
                self.viewer_panel.set_config(json.loads((base / config).read_text(encoding="utf-8")))
        except (OSError, ValueError) as error:
            self._problem(f"Could not load demo viewer settings: {error}")
            return
        self.pick_source.set_path(str(source.resolve()))
        self.pick_out.set_path(str(self._default_output(profile.get("id", source.name))))
        for picker, key in ((self.calibration_panel.pick_csv, "calibrant_csv"), (self.analyte_panel.pick_ana, "analyte_csv")):
            picker.set_path(str((base / profile[key]).resolve()) if profile.get(key) else "")
        if not profile.get("calibrant_csv"):
            self.calibration_panel.table.load_dataframe(pd.DataFrame())
        if not profile.get("analyte_csv"):
            self.analyte_panel.table.load_dataframe(pd.DataFrame())
        self.combo_method.setCurrentText(profile.get("method", "direct"))
        self.chk_ccs.setChecked(False)
        description = profile.get("description", "")
        if profile.get("id") == "protein_adp":
            fit = re.search(r"\((0\.\d+)\)", description)
            scope = "Review needed: protein calibration fit" + (f" {fit.group(1)}" if fit else " needs review") + "; ADP CCS extrapolates."
        elif profile.get("id") == "small_molecule":
            scope = "CCS extrapolates: Betaine lies below the calibration drift range."
        elif profile.get("ccs_scope") == "viewer_only":
            scope = "Viewer only: no compatible CCS calibration is bundled."
        else:
            scope = description
        count = profile.get("raw_count")
        self.demo_description.setText(scope + (f"\n{count} complete real acquisitions." if count else ""))
        self.demo_description.setToolTip(description)
        self.ccs_demo_description.setText(profile.get("description", ""))
        self._update_ccs_cards()
        self.job_status.setText("Real acquisition demo loaded. Process & view creates fresh offline HTML in your output folder.")
        self.tabs.setCurrentIndex(0)
        self._sync_ccs_paths()
        self.reload_runs()

    @staticmethod
    def _default_output(name):
        documents = QStandardPaths.writableLocation(QStandardPaths.DocumentsLocation)
        return Path(documents or str(Path.home())) / "DeconVoVo" / "Waters" / name

    def _pick_text(self):
        file, _ = QFileDialog.getOpenFileName(self, "Select an existing converted Waters MS or IM file", self.pick_source.path(), "Waters converted TXT (*_ms.txt *_im.txt);;Text files (*.txt)")
        if file:
            self.pick_source.set_path(str(Path(file).parent))

    def _source_changed(self, _text=""):
        text = self.pick_source.path()
        if text and (not self.pick_out.path() or self.pick_out.path() == self._auto_output):
            self._auto_output = str(self._default_output(Path(text).stem or "results"))
            self.pick_out.set_path(self._auto_output)
        self._last_result = {}
        self.demo_description.setText("")
        self.ccs_demo_description.setText("")
        self._paths_changed()

    def _paths_changed(self, _text=""):
        self._sync_ccs_paths()
        self.reload_runs()

    def _sync_ccs_paths(self, _text=""):
        source = Path(self.pick_source.path()) if self.pick_source.path() else None
        out = Path(self.pick_out.path()) if self.pick_out.path() else None
        converted = Path(self._last_result["converted_dir"]) if self._last_result.get("converted_dir") else (out / "_converted" if out else None)
        if source and source.is_dir() and any(source.glob("*_im.txt")):
            converted = source
        ccs = out / f"ccs_{self.combo_method.currentText()}" if out else None
        for panel in (self.calibration_panel, self.analyte_panel):
            panel.pick_data.set_path(str(converted) if converted else "")
            panel.pick_out.set_path(str(ccs) if ccs else "")
            panel.combo_method.setCurrentText(self.combo_method.currentText())
        self.analyte_panel.pick_cal.set_path(self.calibration_panel.pick_csv.path())

    def reload_runs(self):
        current = self.run_table.currentRow()
        selected_name = self._rows[current]["run_name"] if 0 <= current < len(self._rows) else None
        source = Path(self.pick_source.path()) if self.pick_source.path() else None
        out = Path(self.pick_out.path()) if self.pick_out.path() else None
        try:
            self._rows = discover_runs(source, out)
        except OSError as error:
            self._rows = []
            self._problem(f"Could not inspect input: {error}")
        self.run_table.blockSignals(True)
        self.run_table.setRowCount(len(self._rows))
        returned = {run.get("run_name", run.get("name", "")): run for run in self._last_result.get("runs", []) if isinstance(run, dict)}
        for i, record in enumerate(self._rows):
            paired = bool(record.get("ms") and record.get("im"))
            data_status = "MS + IM" if paired else "Missing IM" if record.get("ms") else "Missing MS" if record.get("im") else "Not converted"
            status = "Viewer ready" if record.get("viewer") else "Ready to view" if paired else "Needs conversion" if record.get("raw") else "Incomplete TXT pair"
            info = returned.get(record["run_name"], {})
            detail = info.get("status")
            if detail:
                detail = " · ".join(map(str, detail)) if isinstance(detail, list) else str(detail)
                if any(token in detail.casefold() for token in ("error", "fail", "missing")):
                    status = detail
            values = [record["run_name"], "RAW acquisition" if record.get("raw") else "Converted TXT", data_status, "Offline HTML" if record.get("viewer") else "—", status]
            record["display_status"] = status
            for column, value in enumerate(values):
                item = QTableWidgetItem(value)
                item.setToolTip("\n".join(str(record[key]) for key in ("raw", "ms", "im", "viewer") if record.get(key)))
                self.run_table.setItem(i, column, item)
            name = record["run_name"]
            short_name = name[9:] if len(name) > 9 and name[:8].isdigit() and name[8] == "_" else name
            self.run_table.item(i, 0).setText(short_name + "\n" + ("RAW acquisition" if record.get("raw") else "Converted TXT"))
        raw_count = sum(bool(row.get("raw")) for row in self._rows)
        pair_count = sum(bool(row.get("ms") and row.get("im")) for row in self._rows)
        viewer_count = sum(bool(row.get("viewer")) for row in self._rows)
        self.overview.setText(f"{len(self._rows)} runs  ·  {raw_count} RAW acquisitions  ·  {pair_count} MS/IM pairs  ·  {viewer_count} viewers" if self._rows else "No runs loaded. Choose a source folder or load a real demo above.")
        self.queue_summary.setText(f"{len(self._rows)} inputs ready" if self._rows else "No inputs loaded")
        self.output_summary.setText(f"{len(self._rows)} runs · {viewer_count} viewers" if self._rows else "No outputs yet")
        self.output_summary.setToolTip(f"{pair_count} converted MS / IM pairs; {viewer_count} offline HTML viewers. CCS summaries appear after optional calibration.")
        if self._rows:
            selected = next((i for i, record in enumerate(self._rows) if record["run_name"] == selected_name), 0)
            self.run_table.selectRow(selected)
        self.run_table.blockSignals(False)
        self._selected_run_changed()
        self._update_ccs_cards()
        self._reload_ccs_results()

    def _update_open_buttons(self):
        index = self.run_table.currentRow()
        self.btn_viewer.setEnabled(0 <= index < len(self._rows) and bool(self._rows[index].get("viewer")))
        self.btn_export.setEnabled(not self._busy and 0 <= index < len(self._rows) and all(self._rows[index].get(key) for key in ("ms", "im")))
        self.btn_folder.setEnabled(bool(self.pick_out.path()) and Path(self.pick_out.path()).is_dir())

    def open_viewer(self):
        index = self.run_table.currentRow()
        if 0 <= index < len(self._rows) and self._rows[index].get("viewer"):
            QDesktopServices.openUrl(QUrl.fromLocalFile(str(self._rows[index]["viewer"].resolve())))

    def open_output_folder(self):
        path = Path(self.pick_out.path()) if self.pick_out.path() else None
        if path and path.is_dir():
            QDesktopServices.openUrl(QUrl.fromLocalFile(str(path.resolve())))

    @property
    def is_busy(self):
        return self._busy

    def wait_for_worker(self):
        if self._worker and self._worker.isRunning():
            self._worker.wait()

    def process_and_view(self):
        self._start(html=True, ccs=self.chk_ccs.isChecked())

    def calibrate_and_apply(self):
        self._start(html=False, ccs=True)

    def _problem(self, message):
        self.job_status.setText(message)
        self.job_status.setStyleSheet(f"color: {ERROR};")
        self._log.log(message, "error")

    def _snapshot_table(self, panel, picker, filename, required):
        frame = panel.table.to_dataframe()
        missing = required - set(frame.columns)
        if frame.empty or missing:
            label = "Calibrant" if filename.startswith("calibrant") else "Analyte"
            raise ValueError(f"{label} table needs rows and columns: {', '.join(sorted(required))}.")
        if filename.startswith("calibrant") and len(frame) < 3:
            raise ValueError("CCS calibration needs at least three calibrant rows.")
        original = panel.table._csv_path or (Path(picker.path()) if picker.path() else None)
        if "raw_path" in frame.columns:
            def absolute_reference(value):
                text = str(value).strip()
                if not text:
                    return text
                path = Path(text)
                if original and not path.is_absolute():
                    relative = original.parent / path
                    path = relative if relative.exists() or not path.exists() else path
                    return str(path.resolve())
                return text
            frame["raw_path"] = frame["raw_path"].map(absolute_reference)
        path = Path(self._job_files.name) / filename
        frame.to_csv(path, index=False)
        return path

    def _start(self, *, html, ccs):
        if self._busy:
            return
        source_text, out_text = self.pick_source.path(), self.pick_out.path()
        source = Path(source_text) if source_text else None
        if not source or not source.exists():
            self._problem("Choose an existing complete Waters RAW folder or converted TXT source.")
            return
        if not out_text:
            self._problem("Choose an output folder.")
            return
        out = Path(out_text)
        if out.exists() and not out.is_dir():
            self._problem("The output path is a file. Choose a folder.")
            return
        try:
            rows = discover_runs(source)
        except OSError as error:
            self._problem(f"Could not inspect the source folder: {error}")
            return
        if not rows:
            self._problem("No Waters RAW acquisitions or MS/IM text pairs were found in the source.")
            return
        incomplete = [row["run_name"] for row in rows if not row.get("raw") and not (row.get("ms") and row.get("im"))]
        if incomplete:
            self._problem(f"Incomplete converted data: {', '.join(incomplete)}. Each run needs both _ms.txt and _im.txt.")
            return
        self._job_files = tempfile.TemporaryDirectory(prefix="deconvovo-waters-")
        try:
            config = self.viewer_panel.current_config() if html else None
            config_path = Path(self._job_files.name) / "viewer_config.json" if html else None
            if config_path:
                config_path.write_text(json.dumps(config, indent=2), encoding="utf-8")
            cal_csv = self._snapshot_table(self.calibration_panel, self.calibration_panel.pick_csv, "calibrants.csv", {"name", "mw", "z", "mz", "ccs", "raw_path"}) if ccs else None
            analyte_frame = self.analyte_panel.table.to_dataframe()
            ana_csv = self._snapshot_table(self.analyte_panel, self.analyte_panel.pick_ana, "analytes.csv", {"species", "mw", "z", "mz", "raw_path"}) if ccs and not analyte_frame.empty else None
            kwargs = dict(convert=any(row.get("raw") for row in rows), html=html, calibrant_csv=cal_csv, analyte_csv=ana_csv, method=self.combo_method.currentText(), config_path=config_path, skip_existing=self.chk_reuse.isChecked(), func=self.spin_func.value(), ms_bin=self.spin_ms_bin.value(), im_bin=self.spin_im_bin.value())
        except (OSError, ValueError) as error:
            self._cleanup_job_files()
            self._problem(str(error))
            return
        self._last_result = {}
        self._set_busy(True)
        self.job_status.setText("Processing Waters acquisitions…")
        self.job_status.setStyleSheet(f"color: {TEXT_DIM};")
        self.job_progress.setVisible(True)
        self.job_progress.setRange(0, 0)
        self._worker = Worker(self._do_run, source, out, kwargs)
        self._worker.progress.connect(self._on_progress)
        self._worker.log.connect(self._on_job_log)
        self._worker.finished.connect(self._on_done)
        self._worker.error.connect(self._on_error)
        self._worker.start()

    def _set_busy(self, busy):
        self._busy = busy
        for widget in (self.sources_group, self.btn_run, self.btn_ccs, self.chk_ccs, self.display_group, self.advanced, self.btn_advanced, self.btn_reload, self.ccs_editors, self.combo_method):
            widget.setEnabled(not busy)
        self._update_open_buttons()

    def _do_run(self, source, out, kwargs):
        from deconvovo.waters_workflow import run
        logger = logging.getLogger("ccs_calibrate")
        handler = self._log.get_handler()
        logger.addHandler(handler)
        def progress(current, total, message=""):
            self._worker.progress.emit(int(current), int(total))
            if message:
                self._worker.log.emit(str(message), "info")
        try:
            return run(source, out, progress=progress, **kwargs)
        finally:
            logger.removeHandler(handler)

    def _on_progress(self, current, total):
        self.job_progress.setRange(0, max(1, total))
        self.job_progress.setValue(current)

    def _on_job_log(self, message, level):
        self.job_status.setText(message)
        states = (("Converting ", "Converting"), ("Converted ", "Converted"),
                  ("Generating viewer: ", "Creating HTML"), ("Viewer ready: ", "Viewer ready"))
        for prefix, state in states:
            if message.startswith(prefix):
                name = message.removeprefix(prefix)
                for index, record in enumerate(self._rows):
                    if record["run_name"] == name:
                        record["display_status"] = state
                        self.run_table.item(index, 4).setText(state)
                        if index == self.run_table.currentRow():
                            self.selected_state_label.setText(state)
                        break
                break
        self._log.log(message, level)

    def _cleanup_job_files(self):
        if self._job_files:
            self._job_files.cleanup()
            self._job_files = None

    def _on_done(self, result):
        self.wait_for_worker()
        self._last_result = result or {}
        self._cleanup_job_files()
        self._set_busy(False)
        self.job_progress.setVisible(False)
        self._sync_ccs_paths()
        self.reload_runs()
        calibration = self._last_result.get("calibration")
        if calibration:
            self.calibration_panel._cal = calibration
            self.calibration_panel._draw_calibration(calibration)
            self._reload_drift_profiles()
            self._show_fit_status(calibration)
        count = sum(bool(row.get("viewer")) for row in self._rows)
        message = f"Complete: {len(self._rows)} runs" + (f", {count} offline viewers." if count else ".")
        if self._last_result.get("ccs_dir"):
            message += " CCS results are ready in CCS analysis."
            if calibration and calibration.get("r2_lnln", 0) <= 0.95:
                message += " Calibration fit needs review before interpreting analyte CCS."
        self.job_status.setText(message)
        color = WARNING if calibration and calibration.get("r2_lnln", 0) <= 0.95 else SUCCESS
        self.job_status.setStyleSheet(f"color: {color};")
        self._log.log(message, "success")

    def _on_error(self, message):
        self.wait_for_worker()
        self._cleanup_job_files()
        self._set_busy(False)
        self.job_progress.setVisible(False)
        self.reload_runs()
        if self.calibration_panel._cal:
            self.calibration_panel.lbl_status.setText(
                "Previous CCS results; the latest workflow failed. "
                + self.calibration_panel.lbl_status.text())
        self.job_status.setText("Processing failed. Review the log below; completed files remain available.")
        self.job_status.setStyleSheet(f"color: {ERROR};")
        self._log.log(f"Waters processing failed: {message}", "error")

    def _reload_ccs_results(self):
        if not self.pick_out.path():
            return
        ccs = Path(self.pick_out.path()) / f"ccs_{self.combo_method.currentText()}"
        for panel, filename in ((self.calibration_panel, "calibrant_summary.csv"), (self.analyte_panel, "analyte_summary.csv")):
            path = ccs / filename
            if path.is_file():
                try:
                    frame = pd.read_csv(path)
                    panel.summary_table.load_dataframe(frame, path)
                except (OSError, ValueError) as error:
                    self._log.log(f"Could not reload {filename}: {error}", "warning")
            else:
                panel.summary_table.load_dataframe(pd.DataFrame())
        curve = ccs / "calibration_curve.json"
        summary = ccs / "calibrant_summary.csv"
        if curve.is_file() and summary.is_file():
            try:
                calibration = json.loads(curve.read_text(encoding="utf-8"))
                calibration["calibrant_summary"] = pd.read_csv(summary).to_dict("records")
                self.calibration_panel._cal = calibration
                self.calibration_panel._draw_calibration(calibration)
                self._show_fit_status(calibration)
                self._reload_drift_profiles()
            except (OSError, ValueError, KeyError) as error:
                self._log.log(f"Could not reload calibration diagnostics: {error}", "warning")
        else:
            self.calibration_panel._cal = None
            self.calibration_panel.lbl_status.setText("")
            self.calibration_panel.plot_cal.clear()
            self.calibration_panel.plot_drift.clear()

    def _show_fit_status(self, calibration):
        r2 = calibration.get("r2_lnln", 0)
        color = SUCCESS if r2 > 0.99 else WARNING if r2 > 0.95 else ERROR
        text = f"R² = {r2:.5f}  ·  {calibration.get('n_points', 0)} calibrants"
        if r2 <= 0.95:
            text += "  ·  Review peak selection and calibration before using CCS"
        self.calibration_panel.lbl_status.setText(text)
        self.calibration_panel.lbl_status.setWordWrap(True)
        self.calibration_panel.lbl_status.setStyleSheet(f"color: {color};")

    def _reload_drift_profiles(self):
        # Keep every calibrant diagnostic available without enlarging the page.
        import matplotlib.image as mpimg
        panel = self.calibration_panel
        directory = Path(panel.pick_out.path()) / "cali_png"
        files = sorted(directory.glob("*.png"))
        if not files:
            panel.plot_drift.clear()
            return
        columns = min(3, len(files))
        rows = (len(files) + columns - 1) // columns
        panel.plot_drift.fig.clear()
        panel.plot_drift.axes = panel.plot_drift.fig.subplots(rows, columns, squeeze=False)
        panel.plot_drift.setMinimumHeight(max(360, rows * 210))
        for i, file in enumerate(files):
            axis = panel.plot_drift.ax(*divmod(i, columns))
            axis.imshow(mpimg.imread(str(file)))
            axis.set_axis_off()
        for i in range(len(files), rows * columns):
            panel.plot_drift.ax(*divmod(i, columns)).set_axis_off()
        panel.plot_drift.refresh()
