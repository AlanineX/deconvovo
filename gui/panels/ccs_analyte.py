"""Analyte CCS panel — apply calibration to analytes."""
from __future__ import annotations

from pathlib import Path
import pandas as pd

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel,
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
        sub = QLabel("Apply a CCS calibration to compute analyte collision cross-sections.")
        sub.setProperty("subtitle", True)
        layout.addWidget(sub)

        # Controls
        self.pick_cal = FilePicker("Calibrant CSV:", dialog_title="Select calibrant CSV")
        layout.addWidget(self.pick_cal)
        self.pick_ana = FilePicker("Analyte CSV:", dialog_title="Select analyte CSV")
        self.pick_ana.path_changed.connect(self._on_csv_changed)
        layout.addWidget(self.pick_ana)
        self.pick_data = DirPicker("Converted dir:", dialog_title="Select directory with _ms.txt/_im.txt")
        layout.addWidget(self.pick_data)

        out_row = QHBoxLayout()
        self.pick_out = DirPicker("Output:", dialog_title="Select output directory")
        out_row.addWidget(self.pick_out, stretch=1)
        out_row.addWidget(QLabel("Method:"))
        self.combo_method = QComboBox()
        self.combo_method.addItems(["direct", "twostep"])
        self.combo_method.setFixedWidth(100)
        out_row.addWidget(self.combo_method)
        layout.addLayout(out_row)

        # Splitter: table | summary
        splitter = QSplitter(Qt.Horizontal)
        self.table = CsvTable(editable=True)
        splitter.addWidget(self.table)
        tabs = QTabWidget()
        self.summary_table = CsvTable(editable=False)
        tabs.addTab(self.summary_table, "Analyte Summary")
        splitter.addWidget(tabs)
        splitter.setStretchFactor(0, 2)
        splitter.setStretchFactor(1, 3)
        layout.addWidget(splitter, stretch=1)

        # Bottom
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

    def _on_csv_changed(self, path):
        if path and Path(path).exists():
            try:
                self.table.load_dataframe(pd.read_csv(path), Path(path))
                self._log.log(f"Loaded analyte CSV: {Path(path).name}")
            except Exception as e:
                self._log.log(f"Failed to load: {e}", "error")

    def _run(self):
        cal, ana = self.pick_cal.path(), self.pick_ana.path()
        data, out = self.pick_data.path(), self.pick_out.path()
        if not all([cal, ana, data, out]):
            self._log.log("Set all fields: calibrant CSV, analyte CSV, converted dir, output.", "warning")
            return
        method = self.combo_method.currentText()
        self.btn_run.setEnabled(False)
        self.progress.setVisible(True)
        self.progress.setRange(0, 0)
        self._log.log(f"Computing analyte CCS ({method})...", "accent")
        self._worker = Worker(self._do_run, cal, ana, data, out, method)
        self._worker.finished.connect(self._on_done)
        self._worker.error.connect(self._on_error)
        self._worker.start()

    def _do_run(self, cal_csv, ana_csv, data_dir, out, method):
        import logging
        logger = logging.getLogger("ccs_calibrate")
        handler = self._log.get_handler()
        logger.addHandler(handler)
        try:
            from deconvovo.imms_ccs_calibrate import run as ccs_run
            return ccs_run(Path(out), Path(data_dir), Path(cal_csv),
                           analyte_csv=Path(ana_csv), conversion_method=method)
        finally:
            logger.removeHandler(handler)

    def _on_done(self, cal):
        self.btn_run.setEnabled(True)
        self.progress.setVisible(False)
        csv = Path(self.pick_out.path()) / "analyte_summary.csv"
        if csv.exists():
            df = pd.read_csv(csv)
            self.summary_table.load_dataframe(df, csv)
            above = df["above_threshold"].sum() if "above_threshold" in df.columns else len(df)
            self.lbl_status.setText(f"{len(df)} analytes, {above} above threshold")
            self.lbl_status.setStyleSheet(f"color: {SUCCESS}; font-weight: bold;")
            self._log.log(f"Analyte CCS complete: {len(df)} rows", "success")

    def _on_error(self, msg):
        self.btn_run.setEnabled(True)
        self.progress.setVisible(False)
        self.lbl_status.setText("Failed")
        self.lbl_status.setStyleSheet(f"color: {ERROR};")
        self._log.log(f"Analyte CCS failed: {msg}", "error")
