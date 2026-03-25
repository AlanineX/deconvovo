"""CCS Calibration panel — the core scientific workflow."""
from __future__ import annotations

import logging
import math
from pathlib import Path

import numpy as np
import pandas as pd

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QGroupBox,
    QPushButton, QProgressBar, QTabWidget, QComboBox, QSplitter,
)
from PySide6.QtCore import Qt

from gui.widgets.file_picker import FilePicker, DirPicker
from gui.widgets.csv_table import CsvTable
from gui.widgets.plot_canvas import PlotCanvas
from gui.widgets.log_panel import LogPanel
from gui.widgets.worker import Worker
from gui.theme import TEXT_BRIGHT, TEXT_DIM, ACCENT, SUCCESS, WARNING, ERROR


class CcsCalibrationPanel(QWidget):
    def __init__(self, log_panel: LogPanel, parent=None):
        super().__init__(parent)
        self._log = log_panel
        self._worker = None
        self._cal = None  # calibration result dict

        layout = QVBoxLayout(self)
        layout.setContentsMargins(12, 12, 12, 12)
        layout.setSpacing(8)

        # Title
        title = QLabel("CCS Calibration")
        title.setProperty("title", True)
        layout.addWidget(title)
        sub = QLabel("Build a TW-IMS CCS calibration curve from calibrant CSV. "
                     "Edit the table directly or load a CSV file.")
        sub.setProperty("subtitle", True)
        sub.setWordWrap(True)
        layout.addWidget(sub)

        # Top: controls
        ctrl = QHBoxLayout()
        self.pick_csv = FilePicker("Calibrant CSV:", dialog_title="Select calibrant CSV")
        self.pick_csv.path_changed.connect(self._on_csv_changed)
        ctrl.addWidget(self.pick_csv, stretch=1)

        ctrl.addWidget(QLabel("Method:"))
        self.combo_method = QComboBox()
        self.combo_method.addItems(["direct", "twostep"])
        self.combo_method.setFixedWidth(100)
        ctrl.addWidget(self.combo_method)

        layout.addLayout(ctrl)

        row2 = QHBoxLayout()
        self.pick_data = DirPicker("Converted dir:", dialog_title="Select directory with converted _ms.txt/_im.txt files")
        row2.addWidget(self.pick_data, stretch=1)
        layout.addLayout(row2)

        out_row = QHBoxLayout()
        self.pick_out = DirPicker("Output:", dialog_title="Select output directory")
        out_row.addWidget(self.pick_out, stretch=1)
        layout.addLayout(out_row)

        # Main splitter: table left, plots right
        splitter = QSplitter(Qt.Horizontal)

        # Left: calibrant table
        self.table = CsvTable(editable=True)
        splitter.addWidget(self.table)

        # Right: plots in tabs
        plot_tabs = QTabWidget()

        # Calibration plot (2-panel: ln-ln + linear)
        self.plot_cal = PlotCanvas(figsize=(8, 3.5), nrows=1, ncols=2)
        plot_tabs.addTab(self.plot_cal, "Calibration Curve")

        # Drift profiles — grid of all calibrants
        self.plot_drift = PlotCanvas(figsize=(8, 6), nrows=3, ncols=3)
        plot_tabs.addTab(self.plot_drift, "Drift Profiles")

        # Summary table
        self.summary_table = CsvTable(editable=False)
        plot_tabs.addTab(self.summary_table, "Calibrant Summary")

        splitter.addWidget(plot_tabs)
        splitter.setStretchFactor(0, 2)
        splitter.setStretchFactor(1, 3)

        layout.addWidget(splitter, stretch=1)

        # Bottom: run + progress
        bot = QHBoxLayout()
        self.lbl_status = QLabel("")
        self.lbl_status.setStyleSheet(f"color: {TEXT_DIM};")
        bot.addWidget(self.lbl_status)
        bot.addStretch()
        self.btn_run = QPushButton("Run Calibration")
        self.btn_run.setFixedWidth(180)
        self.btn_run.setFixedHeight(36)
        self.btn_run.clicked.connect(self._run)
        bot.addWidget(self.btn_run)
        layout.addLayout(bot)

        self.progress = QProgressBar()
        self.progress.setVisible(False)
        layout.addWidget(self.progress)

    def _on_csv_changed(self, path: str):
        if path and Path(path).exists():
            try:
                df = pd.read_csv(path)
                self.table.load_dataframe(df, Path(path))
                self._log.log(f"Loaded {len(df)} calibrants from {Path(path).name}", "info")
            except Exception as e:
                self._log.log(f"Failed to load CSV: {e}", "error")

    def _run(self):
        data = self.pick_data.path()
        out = self.pick_out.path()
        if not data or not out:
            self._log.log("Set data directory and output directory.", "warning")
            return

        # Get current table as DataFrame
        df = self.table.to_dataframe()
        if len(df) < 3:
            self._log.log("Need at least 3 calibrant rows.", "warning")
            return

        # Check required columns
        required = {"name", "mw", "z", "mz", "ccs", "data_dir", "raw_path"}
        missing = required - set(df.columns)
        if missing:
            self._log.log(f"Missing required columns: {missing}", "error")
            return

        method = self.combo_method.currentText()

        self.btn_run.setEnabled(False)
        self.progress.setVisible(True)
        self.progress.setRange(0, 0)
        self.lbl_status.setText("Running calibration...")
        self._log.log(f"Starting CCS calibration ({method})...", "accent")

        # Save table to temp CSV and run
        import tempfile
        tmp = Path(tempfile.mktemp(suffix=".csv"))
        df.to_csv(tmp, index=False)

        self._worker = Worker(self._do_calibrate, str(tmp), data, out, method)
        self._worker.finished.connect(lambda r: self._on_done(r, tmp))
        self._worker.error.connect(lambda e: self._on_error(e, tmp))
        self._worker.start()

    def _do_calibrate(self, csv_path: str, data_dir: str, out_dir: str, method: str):
        import logging
        logger = logging.getLogger("ccs_calibrate")
        handler = self._log.get_handler()
        logger.addHandler(handler)
        try:
            from deconvovo.imms_ccs_calibrate import run as ccs_run
            return ccs_run(Path(out_dir), Path(data_dir), Path(csv_path),
                           conversion_method=method)
        finally:
            logger.removeHandler(handler)

    def _on_done(self, cal, tmp_csv: Path):
        tmp_csv.unlink(missing_ok=True)
        self._cal = cal
        self.btn_run.setEnabled(True)
        self.progress.setVisible(False)

        r2 = cal["r2_lnln"]
        n = cal["n_points"]
        X = cal["X"]
        A = cal["A"]

        color = SUCCESS if r2 > 0.99 else WARNING if r2 > 0.95 else ERROR
        self.lbl_status.setText(
            f"R² = {r2:.5f}  •  {n} points  •  "
            f"Ω' = {A:.2f} · t'_D^{X:.4f}")
        self.lbl_status.setStyleSheet(f"color: {color}; font-weight: bold;")
        self._log.log(f"Calibration complete: R²={r2:.5f}, {n} points", "success")

        # Update calibration plot
        self._draw_calibration(cal)

        # Draw drift profiles from saved PNGs
        self._draw_drift_profiles()

        # Load summary table
        out_dir = self.pick_out.path()
        summary_csv = Path(out_dir) / "calibrant_summary.csv"
        if summary_csv.exists():
            self.summary_table.load_dataframe(pd.read_csv(summary_csv), summary_csv)

    def _on_error(self, msg, tmp_csv: Path):
        tmp_csv.unlink(missing_ok=True)
        self.btn_run.setEnabled(True)
        self.progress.setVisible(False)
        self.lbl_status.setText("Calibration failed")
        self.lbl_status.setStyleSheet(f"color: {ERROR};")
        self._log.log(f"Calibration failed: {msg}", "error")

    def _draw_calibration(self, cal: dict):
        """Draw 2-panel calibration plot in the embedded canvas."""
        self.plot_cal.clear()
        ax1 = self.plot_cal.ax(0, 0)
        ax2 = self.plot_cal.ax(0, 1)

        summary = cal.get("calibrant_summary", [])
        if not summary:
            self.plot_cal.refresh()
            return

        X, A = cal["X"], cal["A"]
        gas_mw = cal["gas_mw"]

        # Panel 1: ln-ln
        lt = [math.log(p["t_prime"]) for p in summary]
        lo = [math.log(p["CCS_literature"] / (p["z"] * math.sqrt(1/p["mw"] + 1/gas_mw)))
              for p in summary]
        ax1.scatter(lt, lo, c="#5b8af5", s=50, zorder=5, edgecolors="white", linewidth=0.5)
        for p, x, y in zip(summary, lt, lo):
            ax1.annotate(p["name"], (x, y), fontsize=7,
                         textcoords="offset points", xytext=(4, 4), alpha=0.8)
        t_range = np.linspace(min(lt) - 0.15, max(lt) + 0.15, 100)
        ax1.plot(t_range, X * t_range + math.log(A), color="#e8584a",
                 linewidth=1.5, alpha=0.8)
        ax1.set_xlabel("ln(t'_D)")
        ax1.set_ylabel("ln(Ω')")
        ax1.set_title(f"ln-ln  R² = {cal['r2_lnln']:.5f}", fontsize=10)
        ax1.grid(True, alpha=0.3)

        # Panel 2: linear diagnostic
        td = [p["t_double_prime"] for p in summary]
        cc = [p["CCS_literature"] for p in summary]
        ax2.scatter(td, cc, c="#5b8af5", s=50, zorder=5, edgecolors="white", linewidth=0.5)
        for p, x, y in zip(summary, td, cc):
            ax2.annotate(p["name"], (x, y), fontsize=7,
                         textcoords="offset points", xytext=(4, 4), alpha=0.8)
        td_range = np.linspace(min(td) * 0.9, max(td) * 1.1, 100)
        sl = cal["slope_diagnostic"]
        ic = cal["intercept_diagnostic"]
        ax2.plot(td_range, sl * td_range + ic, color="#e8584a",
                 linewidth=1.5, alpha=0.8)
        ax2.set_xlabel("t''_D")
        ax2.set_ylabel("CCS (Å²)")
        ax2.set_title(f"Diagnostic  R² = {cal['r2_linear_diagnostic']:.5f}", fontsize=10)
        ax2.grid(True, alpha=0.3)

        self.plot_cal.refresh()

    def _draw_drift_profiles(self):
        """Load calibrant drift profile PNGs and display in grid."""
        out_dir = Path(self.pick_out.path())
        png_dir = out_dir / "cali_png"
        if not png_dir.exists():
            return

        pngs = sorted(png_dir.glob("*.png"))
        if not pngs:
            return

        # Resize grid to fit the number of profiles
        n = len(pngs)
        ncols = min(3, n)
        nrows = math.ceil(n / ncols)

        # Recreate the canvas with correct grid
        from gui.widgets.plot_canvas import PlotCanvas
        old = self.plot_drift
        parent_tabs = old.parent()
        idx = parent_tabs.indexOf(old)

        self.plot_drift = PlotCanvas(figsize=(8, max(3, nrows * 2.5)), nrows=nrows, ncols=ncols)
        parent_tabs.removeTab(idx)
        parent_tabs.insertTab(idx, self.plot_drift, "Drift Profiles")

        import matplotlib.image as mpimg
        for i, png_path in enumerate(pngs):
            r, c = divmod(i, ncols)
            ax = self.plot_drift.ax(r, c)
            img = mpimg.imread(str(png_path))
            ax.imshow(img)
            ax.set_axis_off()

        # Hide unused axes
        for i in range(n, nrows * ncols):
            r, c = divmod(i, ncols)
            self.plot_drift.ax(r, c).set_visible(False)

        self.plot_drift.refresh()

    def get_calibration(self) -> dict | None:
        """Return current calibration dict (used by Analyte panel)."""
        return self._cal
