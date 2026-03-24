"""Full pipeline panel — run everything end-to-end."""
from __future__ import annotations

from pathlib import Path

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QGroupBox,
    QPushButton, QProgressBar, QCheckBox, QSpinBox, QComboBox,
)

from gui.widgets.file_picker import DirPicker, FilePicker
from gui.widgets.log_panel import LogPanel
from gui.widgets.worker import Worker
from gui.theme import TEXT_DIM, SUCCESS, ERROR


class FullPipelinePanel(QWidget):
    def __init__(self, log_panel: LogPanel, parent=None):
        super().__init__(parent)
        self._log = log_panel
        self._worker = None

        layout = QVBoxLayout(self)
        layout.setContentsMargins(16, 16, 16, 16)
        layout.setSpacing(12)

        title = QLabel("Full Pipeline")
        title.setProperty("title", True)
        layout.addWidget(title)
        sub = QLabel("Run the complete workflow: Convert .raw → HTML Viewer → CCS Calibration → Analyte CCS.")
        sub.setProperty("subtitle", True)
        sub.setWordWrap(True)
        layout.addWidget(sub)

        # Paths
        grp = QGroupBox("Input")
        gl = QVBoxLayout(grp)
        self.pick_raw = DirPicker("Raw .raw dir:", dialog_title="Select directory containing .raw folders")
        self.pick_out = DirPicker("Output:", dialog_title="Select output directory")
        gl.addWidget(self.pick_raw)
        gl.addWidget(self.pick_out)
        layout.addWidget(grp)

        # CCS config
        grp2 = QGroupBox("CCS Configuration")
        gl2 = QVBoxLayout(grp2)
        self.pick_cal_csv = FilePicker("Calibrant CSV:")
        gl2.addWidget(self.pick_cal_csv)
        self.pick_ana_csv = FilePicker("Analyte CSV:")
        gl2.addWidget(self.pick_ana_csv)
        row_opts = QHBoxLayout()
        row_opts.addWidget(QLabel("Method:"))
        self.combo_method = QComboBox()
        self.combo_method.addItems(["direct", "twostep"])
        self.combo_method.setFixedWidth(100)
        row_opts.addWidget(self.combo_method)
        row_opts.addStretch()
        gl2.addLayout(row_opts)
        layout.addWidget(grp2)

        # Steps
        grp3 = QGroupBox("Steps")
        sl = QHBoxLayout(grp3)
        self.chk_convert = QCheckBox("Convert .raw")
        self.chk_convert.setChecked(True)
        sl.addWidget(self.chk_convert)
        self.chk_html = QCheckBox("HTML Viewer")
        self.chk_html.setChecked(True)
        sl.addWidget(self.chk_html)
        self.chk_ccs = QCheckBox("CCS Calibration")
        self.chk_ccs.setChecked(True)
        sl.addWidget(self.chk_ccs)
        self.chk_analyte = QCheckBox("Analyte CCS")
        self.chk_analyte.setChecked(True)
        sl.addWidget(self.chk_analyte)
        sl.addStretch()
        sl.addWidget(QLabel("Workers:"))
        self.spin_workers = QSpinBox()
        self.spin_workers.setRange(1, 16)
        self.spin_workers.setValue(8)
        sl.addWidget(self.spin_workers)
        layout.addWidget(grp3)

        # Run
        run_row = QHBoxLayout()
        self.lbl_status = QLabel("")
        self.lbl_status.setStyleSheet(f"color: {TEXT_DIM};")
        run_row.addWidget(self.lbl_status)
        run_row.addStretch()
        self.btn_run = QPushButton("Run Full Pipeline")
        self.btn_run.setFixedWidth(200)
        self.btn_run.setFixedHeight(40)
        self.btn_run.clicked.connect(self._run)
        run_row.addWidget(self.btn_run)
        layout.addLayout(run_row)

        self.progress = QProgressBar()
        self.progress.setVisible(False)
        layout.addWidget(self.progress)

        layout.addStretch()

    def _run(self):
        raw_dir = self.pick_raw.path()
        out = self.pick_out.path()
        if not raw_dir or not out:
            self._log.log("Set raw directory and output directory.", "warning")
            return

        cal_csv = self.pick_cal_csv.path() or None
        ana_csv = self.pick_ana_csv.path() or None
        method = self.combo_method.currentText()
        workers = self.spin_workers.value()

        steps = {
            "convert": self.chk_convert.isChecked(),
            "html": self.chk_html.isChecked(),
            "ccs": self.chk_ccs.isChecked(),
            "analyte": self.chk_analyte.isChecked(),
        }

        if steps["ccs"] and not cal_csv:
            self._log.log("CCS calibration requires a calibrant CSV.", "warning")
            return
        if steps["analyte"] and not ana_csv:
            self._log.log("Analyte CCS requires an analyte CSV.", "warning")
            return

        self.btn_run.setEnabled(False)
        self.progress.setVisible(True)
        self.progress.setRange(0, 0)
        self.lbl_status.setText("Running pipeline...")
        self._log.log("Starting full pipeline...", "accent")

        self._worker = Worker(
            self._do_run, raw_dir, out, cal_csv, ana_csv,
            method, steps, workers)
        self._worker.finished.connect(self._on_done)
        self._worker.error.connect(self._on_error)
        self._worker.start()

    def _do_run(self, raw_dir, out, cal_csv, ana_csv,
                method, steps, workers):
        import logging
        logger = logging.getLogger("ccs_calibrate")
        handler = self._log.get_handler()
        logger.addHandler(handler)

        out = Path(out)
        out.mkdir(parents=True, exist_ok=True)
        converted_dir = out / "_converted"

        try:
            # Step 1: Convert .raw → text
            if steps["convert"]:
                self._log._sig.message.emit("Step 1: Converting .raw files...", "accent")
                from deconvovo.imms_convert import run as convert_run
                convert_run(Path(raw_dir), converted_dir)

            data_dir = converted_dir if converted_dir.exists() else Path(raw_dir)

            # Step 2: HTML viewers
            if steps["html"]:
                self._log._sig.message.emit("Step 2: Generating HTML viewers...", "accent")
                from deconvovo.imms_plot import run as plot_run
                plot_run(
                    data_dir, out,
                    skip_existing=True,
                    raw_dir=Path(raw_dir),
                    n_workers=workers,
                )

            # Step 3: CCS calibration (+ analytes)
            if steps["ccs"] and cal_csv:
                self._log._sig.message.emit("Step 3: CCS calibration...", "accent")
                from deconvovo.imms_ccs_calibrate import run as ccs_run
                ccs_out = out / f"ccs_{method}"
                ccs_run(
                    ccs_out, Path(cal_csv),
                    analyte_csv=Path(ana_csv) if ana_csv and steps["analyte"] else None,
                    conversion_method=method,
                )

            return {"status": "complete"}
        finally:
            logger.removeHandler(handler)

    def _on_done(self, result):
        self.btn_run.setEnabled(True)
        self.progress.setVisible(False)
        self.lbl_status.setText("Pipeline complete")
        self.lbl_status.setStyleSheet(f"color: {SUCCESS}; font-weight: bold;")
        self._log.log("Full pipeline complete.", "success")

    def _on_error(self, msg):
        self.btn_run.setEnabled(True)
        self.progress.setVisible(False)
        self.lbl_status.setText("Pipeline failed")
        self.lbl_status.setStyleSheet(f"color: {ERROR};")
        self._log.log(f"Pipeline failed: {msg}", "error")
