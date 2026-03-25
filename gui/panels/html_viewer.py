"""2D IM-MS HTML viewer generation panel."""
from __future__ import annotations

import json
import sys
from pathlib import Path

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QGroupBox,
    QPushButton, QProgressBar, QCheckBox, QSpinBox, QComboBox,
)

from gui.widgets.file_picker import DirPicker
from gui.widgets.log_panel import LogPanel
from gui.widgets.worker import Worker


def _load_config():
    """Load imms_plot_config.json for populating GUI defaults."""
    if getattr(sys, 'frozen', False):
        base = Path(sys._MEIPASS)
    else:
        base = Path(__file__).parent.parent.parent
    cfg_path = base / "config" / "imms_plot_config.json"
    if cfg_path.exists():
        with open(cfg_path) as f:
            return json.load(f)
    return {}


class HtmlViewerPanel(QWidget):
    def __init__(self, log_panel: LogPanel, parent=None):
        super().__init__(parent)
        self._log = log_panel
        self._worker = None

        layout = QVBoxLayout(self)
        layout.setContentsMargins(16, 16, 16, 16)
        layout.setSpacing(12)

        title = QLabel("2D IM-MS Interactive Viewer")
        title.setProperty("title", True)
        layout.addWidget(title)
        sub = QLabel("Generate interactive HTML heatmaps with linked drift and m/z panels.")
        sub.setProperty("subtitle", True)
        layout.addWidget(sub)

        # Input
        grp_in = QGroupBox("Input")
        gl_in = QVBoxLayout(grp_in)
        self.pick_data = DirPicker("Text dir:",
            dialog_title="Select directory with _ms.txt / _im.txt files",
            hint="Directory with converted text files (_ms.txt and _im.txt)")
        self.pick_raw = DirPicker("Raw dir:",
            dialog_title="Select directory with .raw folders",
            hint="Provides pusher period for drift time axis (leave empty to use bin index)")
        gl_in.addWidget(self.pick_data)
        gl_in.addWidget(self.pick_raw)
        layout.addWidget(grp_in)

        # Output
        grp_out = QGroupBox("Output")
        gl_out = QVBoxLayout(grp_out)
        self.pick_out = DirPicker("Output dir:",
            dialog_title="Select output directory for HTML files",
            hint="HTML viewer files and companion CSVs will be written here")
        gl_out.addWidget(self.pick_out)
        layout.addWidget(grp_out)

        # Options — load defaults from config
        cfg = _load_config()
        defaults = cfg.get("defaults", {})
        presets = cfg.get("presets", {})

        opts = QGroupBox("Options")
        ol = QVBoxLayout(opts)

        # Row 1: skip, workers
        row1 = QHBoxLayout()
        self.chk_skip = QCheckBox("Skip existing")
        self.chk_skip.setChecked(True)
        row1.addWidget(self.chk_skip)
        row1.addStretch()
        row1.addWidget(QLabel("Workers:"))
        self.spin_workers = QSpinBox()
        self.spin_workers.setRange(1, 16)
        self.spin_workers.setValue(8)
        row1.addWidget(self.spin_workers)
        ol.addLayout(row1)

        # Row 2: m/z bins, colormap
        row2 = QHBoxLayout()
        row2.addWidget(QLabel("m/z bins:"))
        self.spin_mzbins = QSpinBox()
        self.spin_mzbins.setRange(100, 25600)
        self.spin_mzbins.setSingleStep(200)
        self.spin_mzbins.setValue(defaults.get("mz_bins", 3200))
        row2.addWidget(self.spin_mzbins)
        row2.addSpacing(16)
        row2.addWidget(QLabel("Colormap:"))
        self.combo_cmap = QComboBox()
        cmap_names = []
        for c in presets.get("colormap", ["Viridis"]):
            if isinstance(c, str):
                cmap_names.append(c)
            elif isinstance(c, dict) and "name" in c:
                cmap_names.append(c["name"])
        self.combo_cmap.addItems(cmap_names)
        default_cmap = defaults.get("colormap", "Viridis")
        if default_cmap in cmap_names:
            self.combo_cmap.setCurrentText(default_cmap)
        row2.addWidget(self.combo_cmap)
        row2.addStretch()
        ol.addLayout(row2)

        # Row 3: figure size
        row3 = QHBoxLayout()
        fig_cfg = cfg.get("figure", {})
        row3.addWidget(QLabel("Width:"))
        self.spin_width = QSpinBox()
        self.spin_width.setRange(600, 2400)
        self.spin_width.setSingleStep(100)
        self.spin_width.setValue(fig_cfg.get("width", 1190))
        row3.addWidget(self.spin_width)
        row3.addSpacing(16)
        row3.addWidget(QLabel("Height:"))
        self.spin_height = QSpinBox()
        self.spin_height.setRange(400, 1600)
        self.spin_height.setSingleStep(100)
        self.spin_height.setValue(fig_cfg.get("height", 680))
        row3.addWidget(self.spin_height)
        row3.addStretch()
        ol.addLayout(row3)

        layout.addWidget(opts)

        run_row = QHBoxLayout()
        run_row.addStretch()
        self.btn_run = QPushButton("Generate HTML Plots")
        self.btn_run.setFixedWidth(200)
        self.btn_run.setFixedHeight(36)
        self.btn_run.clicked.connect(self._run)
        run_row.addWidget(self.btn_run)
        layout.addLayout(run_row)

        self.progress = QProgressBar()
        self.progress.setVisible(False)
        layout.addWidget(self.progress)

        layout.addStretch()

    def _run(self):
        data = self.pick_data.path()
        out = self.pick_out.path()
        if not data or not out:
            self._log.log("Set text directory and output directory.", "warning")
            return
        raw = self.pick_raw.path() or None
        skip = self.chk_skip.isChecked()
        workers = self.spin_workers.value()
        mz_bins = self.spin_mzbins.value()

        self.btn_run.setEnabled(False)
        self.progress.setVisible(True)
        self.progress.setRange(0, 0)
        self._log.log(f"Generating HTML plots ({mz_bins} m/z bins)...", "accent")

        self._worker = Worker(self._do_plot, data, out, raw, skip, workers, mz_bins)
        self._worker.finished.connect(self._on_done)
        self._worker.error.connect(self._on_error)
        self._worker.start()

    def _do_plot(self, data, out, raw, skip, workers, mz_bins):
        data_p = Path(data)
        im_files = list(data_p.glob("*_im.txt"))
        self._log._sig.message.emit(
            f"Found {len(im_files)} IM files in {data_p.name}", "info")
        if not im_files:
            raise FileNotFoundError(f"No _im.txt files found in {data}")

        out_p = Path(out)
        for stale in out_p.glob("*_2d_imms.html"):
            if stale.stat().st_size == 0:
                stale.unlink()
                self._log._sig.message.emit(f"  Removed stale 0-byte: {stale.name}", "warning")

        from deconvovo.imms_plot import run as plot_run
        results = plot_run(
            data_p, out_p,
            skip_existing=skip,
            raw_dir=Path(raw) if raw else None,
            n_workers=workers,
        )
        errors = []
        if results:
            for r in results:
                status = " | ".join(r.get("status", []))
                level = "error" if "ERROR" in status else "info"
                self._log._sig.message.emit(f"  {r['run_name']}: {status}", level)
                if "ERROR" in status:
                    errors.append(f"{r['run_name']}: {status}")
        if errors:
            raise RuntimeError("HTML generation failed:\n" + "\n".join(errors))
        return results

    def _on_done(self, results):
        self.btn_run.setEnabled(True)
        self.progress.setVisible(False)
        out = Path(self.pick_out.path())
        htmls = list(out.glob("*_2d_imms.html"))
        n_ok = sum(1 for h in htmls if h.stat().st_size > 0)
        n_empty = len(htmls) - n_ok
        msg = f"Done: {n_ok} HTML files"
        if n_empty:
            msg += f" ({n_empty} empty)"
        self._log.log(msg, "success" if not n_empty else "warning")

    def _on_error(self, msg):
        self.btn_run.setEnabled(True)
        self.progress.setVisible(False)
        self._log.log(f"HTML generation failed: {msg}", "error")
