"""CCS Calibration panel."""
from __future__ import annotations

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
from gui.theme import TEXT_DIM, ACCENT, SUCCESS, WARNING, ERROR


def _scatter_annotate(ax, xs, ys, names, color="#5b8af5", fit_x=None, fit_y=None):
    """Scatter points with name labels and optional fit line."""
    ax.scatter(xs, ys, c=color, s=50, zorder=5, edgecolors="white", linewidth=0.5)
    for name, x, y in zip(names, xs, ys):
        ax.annotate(name, (x, y), fontsize=7, textcoords="offset points",
                    xytext=(4, 4), alpha=0.8)
    if fit_x is not None:
        ax.plot(fit_x, fit_y, color="#e8584a", linewidth=1.5, alpha=0.8)
    ax.grid(True, alpha=0.3)


class CcsCalibrationPanel(QWidget):
    def __init__(self, log_panel: LogPanel, parent=None):
        super().__init__(parent)
        self._log = log_panel
        self._worker = None
        self._cal = None

        layout = QVBoxLayout(self)
        layout.setContentsMargins(12, 12, 12, 12)
        layout.setSpacing(8)

        title = QLabel("CCS Calibration")
        title.setProperty("title", True)
        layout.addWidget(title)
        sub = QLabel("Build a TW-IMS CCS calibration curve from calibrant CSV.")
        sub.setProperty("subtitle", True)
        layout.addWidget(sub)

        # Controls
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

        self.pick_data = DirPicker("Converted dir:", dialog_title="Select directory with _ms.txt/_im.txt")
        layout.addWidget(self.pick_data)
        self.pick_out = DirPicker("Output:", dialog_title="Select output directory")
        layout.addWidget(self.pick_out)

        # Splitter: table | plots
        splitter = QSplitter(Qt.Horizontal)
        self.table = CsvTable(editable=True)
        splitter.addWidget(self.table)

        tabs = QTabWidget()
        self.plot_cal = PlotCanvas(figsize=(8, 3.5), nrows=1, ncols=2)
        tabs.addTab(self.plot_cal, "Calibration Curve")
        self.plot_drift = PlotCanvas(figsize=(8, 6), nrows=3, ncols=3)
        tabs.addTab(self.plot_drift, "Drift Profiles")
        self.summary_table = CsvTable(editable=False)
        tabs.addTab(self.summary_table, "Calibrant Summary")
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
        self.btn_run = QPushButton("Run Calibration")
        self.btn_run.setFixedWidth(180)
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
                df = pd.read_csv(path)
                self.table.load_dataframe(df, Path(path))
                self._log.log(f"Loaded {len(df)} calibrants from {Path(path).name}")
            except Exception as e:
                self._log.log(f"Failed to load CSV: {e}", "error")

    def _run(self):
        data, out = self.pick_data.path(), self.pick_out.path()
        if not data or not out:
            self._log.log("Set converted directory and output directory.", "warning")
            return
        df = self.table.to_dataframe()
        if len(df) < 3:
            self._log.log("Need at least 3 calibrant rows.", "warning")
            return
        required = {"name", "mw", "z", "mz", "ccs", "raw_path"}
        missing = required - set(df.columns)
        if missing:
            self._log.log(f"Missing columns: {missing}", "error")
            return

        method = self.combo_method.currentText()
        self.btn_run.setEnabled(False)
        self.progress.setVisible(True)
        self.progress.setRange(0, 0)
        self._log.log(f"Starting CCS calibration ({method})...", "accent")

        import tempfile
        tmp = Path(tempfile.mktemp(suffix=".csv"))
        df.to_csv(tmp, index=False)
        self._worker = Worker(self._do_calibrate, str(tmp), data, out, method)
        self._worker.finished.connect(lambda r: self._on_done(r, tmp))
        self._worker.error.connect(lambda e: self._on_error(e, tmp))
        self._worker.start()

    def _do_calibrate(self, csv_path, data_dir, out_dir, method):
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

    def _on_done(self, cal, tmp_csv):
        tmp_csv.unlink(missing_ok=True)
        self._cal = cal
        self.btn_run.setEnabled(True)
        self.progress.setVisible(False)

        r2, n, X, A = cal["r2_lnln"], cal["n_points"], cal["X"], cal["A"]
        color = SUCCESS if r2 > 0.99 else WARNING if r2 > 0.95 else ERROR
        self.lbl_status.setText(f"R² = {r2:.5f}  |  {n} points  |  "
                                f"Ω' = {A:.2f} · t'_D^{X:.4f}")
        self.lbl_status.setStyleSheet(f"color: {color}; font-weight: bold;")
        self._log.log(f"Calibration complete: R²={r2:.5f}, {n} points", "success")

        self._draw_calibration(cal)
        self._draw_drift_profiles()
        out_dir = Path(self.pick_out.path())
        csv = out_dir / "calibrant_summary.csv"
        if csv.exists():
            self.summary_table.load_dataframe(pd.read_csv(csv), csv)

    def _on_error(self, msg, tmp_csv):
        tmp_csv.unlink(missing_ok=True)
        self.btn_run.setEnabled(True)
        self.progress.setVisible(False)
        self.lbl_status.setText("Calibration failed")
        self.lbl_status.setStyleSheet(f"color: {ERROR};")
        self._log.log(f"Calibration failed: {msg}", "error")

    def _draw_calibration(self, cal):
        self.plot_cal.clear()
        summary = cal.get("calibrant_summary", [])
        if not summary:
            self.plot_cal.refresh(); return

        X, A, gas_mw = cal["X"], cal["A"], cal["gas_mw"]
        names = [p["name"] for p in summary]

        # ln-ln
        lt = [math.log(p["t_prime"]) for p in summary]
        lo = [math.log(p["CCS_literature"] / (p["z"] * math.sqrt(1/p["mw"] + 1/gas_mw)))
              for p in summary]
        t_fit = np.linspace(min(lt) - 0.15, max(lt) + 0.15, 100)
        ax1 = self.plot_cal.ax(0, 0)
        _scatter_annotate(ax1, lt, lo, names, fit_x=t_fit, fit_y=X * t_fit + math.log(A))
        ax1.set_xlabel("ln(t'_D)"); ax1.set_ylabel("ln(Ω')")
        ax1.set_title(f"ln-ln  R² = {cal['r2_lnln']:.5f}", fontsize=10)

        # Diagnostic
        td = [p["t_double_prime"] for p in summary]
        cc = [p["CCS_literature"] for p in summary]
        sl, ic = cal["slope_diagnostic"], cal["intercept_diagnostic"]
        td_fit = np.linspace(min(td) * 0.9, max(td) * 1.1, 100)
        ax2 = self.plot_cal.ax(0, 1)
        _scatter_annotate(ax2, td, cc, names, fit_x=td_fit, fit_y=sl * td_fit + ic)
        ax2.set_xlabel("t''_D"); ax2.set_ylabel("CCS (Å²)")
        ax2.set_title(f"Diagnostic  R² = {cal['r2_linear_diagnostic']:.5f}", fontsize=10)

        self.plot_cal.refresh()

    def _draw_drift_profiles(self):
        out_dir = Path(self.pick_out.path())
        pngs = sorted((out_dir / "cali_png").glob("*.png")) if (out_dir / "cali_png").exists() else []
        if not pngs: return

        n = len(pngs)
        ncols = min(3, n)
        nrows = math.ceil(n / ncols)

        from gui.widgets.plot_canvas import PlotCanvas
        old = self.plot_drift
        tabs = old.parent()
        idx = tabs.indexOf(old)
        self.plot_drift = PlotCanvas(figsize=(8, max(3, nrows * 2.5)), nrows=nrows, ncols=ncols)
        tabs.removeTab(idx)
        tabs.insertTab(idx, self.plot_drift, "Drift Profiles")

        import matplotlib.image as mpimg
        for i, png in enumerate(pngs):
            ax = self.plot_drift.ax(*divmod(i, ncols))
            ax.imshow(mpimg.imread(str(png)))
            ax.set_axis_off()
        for i in range(n, nrows * ncols):
            self.plot_drift.ax(*divmod(i, ncols)).set_visible(False)
        self.plot_drift.refresh()

    def get_calibration(self):
        return self._cal
