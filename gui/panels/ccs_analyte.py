"""Analyte CCS panel — apply calibration to analytes."""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QGroupBox,
    QPushButton, QProgressBar, QTabWidget, QComboBox, QSplitter,
)
from PySide6.QtCore import Qt

from gui.widgets.file_picker import FilePicker, DirPicker
from gui.widgets.csv_table import CsvTable
from gui.widgets.log_panel import LogPanel
from gui.widgets.worker import Worker
from gui.theme import TEXT_DIM, SUCCESS, ERROR


class CcsAnalytePanel(QWidget):
    def __init__(self, log_panel: LogPanel, parent=None):
        super().__init__(parent)
        self._log = log_panel
        self._worker = None

        layout = QVBoxLayout(self)
        layout.setContentsMargins(12, 12, 12, 12)
        layout.setSpacing(8)

        title = QLabel("Analyte CCS")
        title.setProperty("title", True)
        layout.addWidget(title)
        sub = QLabel("Apply a CCS calibration to compute analyte collision cross-sections. "
                     "Requires a completed calibration (from CCS Calibration panel or JSON file).")
        sub.setProperty("subtitle", True)
        sub.setWordWrap(True)
        layout.addWidget(sub)

        # Controls
        ctrl = QVBoxLayout()
        row1 = QHBoxLayout()
        self.pick_cal_csv = FilePicker("Calibrant CSV:", dialog_title="Select calibrant CSV")
        row1.addWidget(self.pick_cal_csv, stretch=1)
        ctrl.addLayout(row1)

        row2 = QHBoxLayout()
        self.pick_ana_csv = FilePicker("Analyte CSV:", dialog_title="Select analyte CSV")
        row2.addWidget(self.pick_ana_csv, stretch=1)
        ctrl.addLayout(row2)

        row3 = QHBoxLayout()
        self.pick_out = DirPicker("Output:", dialog_title="Select output directory")
        row3.addWidget(self.pick_out, stretch=1)
        row3.addWidget(QLabel("Method:"))
        self.combo_method = QComboBox()
        self.combo_method.addItems(["direct", "twostep"])
        self.combo_method.setFixedWidth(100)
        row3.addWidget(self.combo_method)
        ctrl.addLayout(row3)

        layout.addLayout(ctrl)

        # Splitter: analyte table | results
        splitter = QSplitter(Qt.Horizontal)

        self.table = CsvTable(editable=True)
        self.pick_ana_csv.path_changed.connect(self._on_ana_csv_changed)
        splitter.addWidget(self.table)

        # Results tabs
        tabs = QTabWidget()
        self.summary_table = CsvTable(editable=False)
        tabs.addTab(self.summary_table, "Analyte Summary")
        splitter.addWidget(tabs)

        splitter.setStretchFactor(0, 2)
        splitter.setStretchFactor(1, 3)
        layout.addWidget(splitter, stretch=1)

        # Run
        bot = QHBoxLayout()
        self.lbl_status = QLabel("")
        self.lbl_status.setStyleSheet(f"color: {TEXT_DIM};")
        bot.addWidget(self.lbl_status)
        bot.addStretch()
        self.btn_run = QPushButton("Compute Analyte CCS")
        self.btn_run.setFixedWidth(200)
        self.btn_run.setFixedHeight(36)
        self.btn_run.clicked.connect(self._run)
        bot.addWidget(self.btn_run)
        layout.addLayout(bot)

        self.progress = QProgressBar()
        self.progress.setVisible(False)
        layout.addWidget(self.progress)

    def _on_ana_csv_changed(self, path: str):
        if path and Path(path).exists():
            try:
                df = pd.read_csv(path)
                self.table.load_dataframe(df, Path(path))
                self._log.log(f"Loaded {len(df)} analyte rows from {Path(path).name}")
            except Exception as e:
                self._log.log(f"Failed to load analyte CSV: {e}", "error")

    def _run(self):
        cal_csv = self.pick_cal_csv.path()
        ana_csv = self.pick_ana_csv.path()
        out = self.pick_out.path()

        if not cal_csv or not ana_csv or not out:
            self._log.log("Set calibrant CSV, analyte CSV, and output directory.", "warning")
            return

        method = self.combo_method.currentText()
        self.btn_run.setEnabled(False)
        self.progress.setVisible(True)
        self.progress.setRange(0, 0)
        self._log.log(f"Computing analyte CCS ({method})...", "accent")

        self._worker = Worker(self._do_run, cal_csv, ana_csv, out, method)
        self._worker.finished.connect(self._on_done)
        self._worker.error.connect(self._on_error)
        self._worker.start()

    def _do_run(self, cal_csv, ana_csv, out, method):
        import logging
        logger = logging.getLogger("ccs_calibrate")
        handler = self._log.get_handler()
        logger.addHandler(handler)
        try:
            from deconvovo.imms_ccs_calibrate import run as ccs_run
            return ccs_run(
                Path(out), Path(cal_csv),
                analyte_csv=Path(ana_csv),
                conversion_method=method,
            )
        finally:
            logger.removeHandler(handler)

    def _on_done(self, cal):
        self.btn_run.setEnabled(True)
        self.progress.setVisible(False)

        out = Path(self.pick_out.path())
        summary_csv = out / "analyte_summary.csv"
        if summary_csv.exists():
            df = pd.read_csv(summary_csv)
            self.summary_table.load_dataframe(df, summary_csv)
            n = len(df)
            above = df["above_threshold"].sum() if "above_threshold" in df.columns else n
            self.lbl_status.setText(f"{n} analytes, {above} above threshold")
            self.lbl_status.setStyleSheet(f"color: {SUCCESS}; font-weight: bold;")
            self._log.log(f"Analyte CCS complete: {n} rows, {above} above threshold", "success")
        else:
            self._log.log("Analyte CCS complete (no summary file found).", "success")

    def _on_error(self, msg):
        self.btn_run.setEnabled(True)
        self.progress.setVisible(False)
        self.lbl_status.setText("Failed")
        self.lbl_status.setStyleSheet(f"color: {ERROR};")
        self._log.log(f"Analyte CCS failed: {msg}", "error")
