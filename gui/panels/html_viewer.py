"""2D IM-MS HTML viewer generation panel."""
from __future__ import annotations

import json
import sys
from pathlib import Path

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QGroupBox,
    QPushButton, QProgressBar, QCheckBox, QSpinBox, QComboBox,
    QDoubleSpinBox, QFormLayout,
)

from gui.widgets.file_picker import DirPicker
from gui.widgets.log_panel import LogPanel
from gui.widgets.worker import Worker


def _config_path():
    if getattr(sys, 'frozen', False):
        return Path(sys._MEIPASS) / "config" / "imms_plot_config.json"
    return Path(__file__).parent.parent.parent / "config" / "imms_plot_config.json"


def _load_config():
    p = _config_path()
    if p.exists():
        with open(p) as f:
            return json.load(f)
    return {}


def _save_config(cfg):
    p = _config_path()
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w") as f:
        json.dump(cfg, f, indent=2, ensure_ascii=True)


class HtmlViewerPanel(QWidget):
    def __init__(self, log_panel: LogPanel, parent=None):
        super().__init__(parent)
        self._log = log_panel
        self._worker = None
        self._cfg = _load_config()

        layout = QVBoxLayout(self)
        layout.setContentsMargins(16, 16, 16, 16)
        layout.setSpacing(10)

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

        # Config options — full imms_plot_config.json editor
        defaults = self._cfg.get("defaults", {})
        presets = self._cfg.get("presets", {})
        fig_cfg = self._cfg.get("figure", {})

        opts = QGroupBox("Viewer Configuration")
        ol = QVBoxLayout(opts)

        # Row 1: processing
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

        # Row 2: m/z bins, drift bins, colormap
        row2 = QHBoxLayout()
        row2.addWidget(QLabel("m/z bins:"))
        self.spin_mzbins = QSpinBox()
        self.spin_mzbins.setRange(100, 25600)
        self.spin_mzbins.setSingleStep(200)
        self.spin_mzbins.setValue(defaults.get("mz_bins", 3200))
        row2.addWidget(self.spin_mzbins)

        row2.addWidget(QLabel("Colormap:"))
        self.combo_cmap = QComboBox()
        cmap_names = []
        for c in presets.get("colormap", ["Viridis"]):
            name = c if isinstance(c, str) else c.get("name", "?")
            cmap_names.append(name)
        self.combo_cmap.addItems(cmap_names)
        default_cmap = defaults.get("colormap", "Viridis")
        if default_cmap in cmap_names:
            self.combo_cmap.setCurrentText(default_cmap)
        row2.addWidget(self.combo_cmap)

        row2.addWidget(QLabel("Scale:"))
        self.combo_scale = QComboBox()
        self.combo_scale.addItems([str(s) for s in presets.get("scale", ["1"])])
        self.combo_scale.setCurrentText(str(defaults.get("scale", "1")))
        row2.addWidget(self.combo_scale)
        row2.addStretch()
        ol.addLayout(row2)

        # Row 3: smoothing defaults
        row3 = QHBoxLayout()
        row3.addWidget(QLabel("2D smooth:"))
        self.combo_smooth2d = QComboBox()
        for p in presets.get("smooth_2d", []):
            self.combo_smooth2d.addItem(p.get("label", "Raw"))
        s2d = defaults.get("smooth_2d", {})
        for i, p in enumerate(presets.get("smooth_2d", [])):
            if {k: v for k, v in p.items() if k != "label"} == s2d:
                self.combo_smooth2d.setCurrentIndex(i)
        row3.addWidget(self.combo_smooth2d)

        row3.addWidget(QLabel("Drift smooth:"))
        self.combo_smooth_drift = QComboBox()
        for p in presets.get("smooth_drift", []):
            self.combo_smooth_drift.addItem(p.get("label", "Raw"))
        sd = defaults.get("smooth_drift", {})
        for i, p in enumerate(presets.get("smooth_drift", [])):
            if {k: v for k, v in p.items() if k != "label"} == sd:
                self.combo_smooth_drift.setCurrentIndex(i)
        row3.addWidget(self.combo_smooth_drift)
        row3.addStretch()
        ol.addLayout(row3)

        # Row 4: noise thresholds
        row4 = QHBoxLayout()
        row4.addWidget(QLabel("Noise 2D (%):"))
        self.spin_noise2d = QDoubleSpinBox()
        self.spin_noise2d.setRange(0, 10)
        self.spin_noise2d.setSingleStep(0.1)
        self.spin_noise2d.setValue(defaults.get("noise_2d", 0.5))
        row4.addWidget(self.spin_noise2d)

        row4.addWidget(QLabel("Drift:"))
        self.spin_noise_drift = QDoubleSpinBox()
        self.spin_noise_drift.setRange(0, 10)
        self.spin_noise_drift.setSingleStep(0.5)
        self.spin_noise_drift.setValue(defaults.get("noise_drift", 0))
        row4.addWidget(self.spin_noise_drift)

        row4.addWidget(QLabel("m/z:"))
        self.spin_noise_mz = QDoubleSpinBox()
        self.spin_noise_mz.setRange(0, 10)
        self.spin_noise_mz.setSingleStep(0.5)
        self.spin_noise_mz.setValue(defaults.get("noise_mz", 0))
        row4.addWidget(self.spin_noise_mz)
        row4.addStretch()
        ol.addLayout(row4)

        # Row 5: figure size + save
        row5 = QHBoxLayout()
        row5.addWidget(QLabel("Fig width:"))
        self.spin_width = QSpinBox()
        self.spin_width.setRange(600, 2400)
        self.spin_width.setSingleStep(100)
        self.spin_width.setValue(fig_cfg.get("width", 1190))
        row5.addWidget(self.spin_width)

        row5.addWidget(QLabel("Height:"))
        self.spin_height = QSpinBox()
        self.spin_height.setRange(400, 1600)
        self.spin_height.setSingleStep(100)
        self.spin_height.setValue(fig_cfg.get("height", 680))
        row5.addWidget(self.spin_height)

        row5.addWidget(QLabel("Font:"))
        self.spin_font = QSpinBox()
        self.spin_font.setRange(10, 32)
        self.spin_font.setValue(fig_cfg.get("font_size", 20))
        row5.addWidget(self.spin_font)

        row5.addStretch()
        btn_save = QPushButton("Save Config")
        btn_save.setProperty("secondary", True)
        btn_save.clicked.connect(self._save_config)
        row5.addWidget(btn_save)
        ol.addLayout(row5)

        layout.addWidget(opts)

        # Run
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

    def _save_config(self):
        """Write current UI values back to imms_plot_config.json."""
        presets = self._cfg.get("presets", {})

        # Read smoothing preset from combo index
        s2d_idx = self.combo_smooth2d.currentIndex()
        s2d_presets = presets.get("smooth_2d", [])
        s2d = {k: v for k, v in s2d_presets[s2d_idx].items() if k != "label"} if s2d_idx < len(s2d_presets) else {"method": "raw"}

        sd_idx = self.combo_smooth_drift.currentIndex()
        sd_presets = presets.get("smooth_drift", [])
        sd = {k: v for k, v in sd_presets[sd_idx].items() if k != "label"} if sd_idx < len(sd_presets) else {"method": "raw"}

        self._cfg["defaults"] = {
            "smooth_2d": s2d,
            "scale": self.combo_scale.currentText(),
            "smooth_drift": sd,
            "noise_2d": self.spin_noise2d.value(),
            "noise_drift": self.spin_noise_drift.value(),
            "noise_mz": self.spin_noise_mz.value(),
            "mz_bins": self.spin_mzbins.value(),
            "colormap": self.combo_cmap.currentText(),
        }
        self._cfg["figure"] = {
            "font_size": self.spin_font.value(),
            "width": self.spin_width.value(),
            "height": self.spin_height.value(),
        }
        _save_config(self._cfg)
        self._log.log("Config saved to imms_plot_config.json", "success")

    def _run(self):
        data = self.pick_data.path()
        out = self.pick_out.path()
        if not data or not out:
            self._log.log("Set text directory and output directory.", "warning")
            return
        # Save config before running so the viewer uses current settings
        self._save_config()

        raw = self.pick_raw.path() or None
        skip = self.chk_skip.isChecked()
        workers = self.spin_workers.value()

        self.btn_run.setEnabled(False)
        self.progress.setVisible(True)
        self.progress.setRange(0, 0)
        self._log.log(f"Generating HTML plots...", "accent")

        self._worker = Worker(self._do_plot, data, out, raw, skip, workers)
        self._worker.finished.connect(self._on_done)
        self._worker.error.connect(self._on_error)
        self._worker.start()

    def _do_plot(self, data, out, raw, skip, workers):
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
        self._log.log(f"Done: {n_ok} HTML files", "success")

    def _on_error(self, msg):
        self.btn_run.setEnabled(True)
        self.progress.setVisible(False)
        self._log.log(f"HTML generation failed: {msg}", "error")
