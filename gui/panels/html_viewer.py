"""2D IM-MS HTML viewer generation panel."""
from __future__ import annotations

import json
import sys
from pathlib import Path

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QGroupBox,
    QPushButton, QProgressBar, QCheckBox, QSpinBox, QComboBox,
    QDoubleSpinBox, QTextEdit, QTabWidget, QFormLayout,
)
from PySide6.QtCore import Qt
from PySide6.QtGui import QFont

from gui.widgets.file_picker import DirPicker
from gui.widgets.log_panel import LogPanel
from gui.widgets.worker import Worker

_MIN_COMBO_W = 200


def _config_path():
    if getattr(sys, 'frozen', False):
        return Path(sys._MEIPASS) / "config" / "imms_plot_config.json"
    return Path(__file__).parent.parent.parent / "config" / "imms_plot_config.json"


def _load_config():
    p = _config_path()
    if p.exists():
        with open(p, encoding="utf-8") as f:
            return json.load(f)
    return {}


def _save_config(cfg):
    p = _config_path()
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w", encoding="utf-8") as f:
        json.dump(cfg, f, indent=2, ensure_ascii=False)


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

        # ---- Input ----
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

        # ---- Output ----
        grp_out = QGroupBox("Output")
        gl_out = QVBoxLayout(grp_out)
        self.pick_out = DirPicker("Output dir:",
            dialog_title="Select output directory for HTML files",
            hint="HTML viewer files and companion CSVs will be written here")
        gl_out.addWidget(self.pick_out)
        layout.addWidget(grp_out)

        # ---- Configuration ----
        defaults = self._cfg.get("defaults", {})
        presets = self._cfg.get("presets", {})
        fig_cfg = self._cfg.get("figure", {})

        cfg_grp = QGroupBox("Viewer Configuration")
        cfg_layout = QVBoxLayout(cfg_grp)
        cfg_tabs = QTabWidget()

        # ==== Tab 1: Defaults ====
        defaults_w = QWidget()
        dl = QVBoxLayout(defaults_w)
        dl.setContentsMargins(8, 8, 8, 8)
        dl_sub = QLabel("Initial viewer settings. Users can change them interactively in the HTML.")
        dl_sub.setProperty("subtitle", True)
        dl.addWidget(dl_sub)

        # Processing row
        proc_row = QHBoxLayout()
        self.chk_skip = QCheckBox("Skip existing")
        self.chk_skip.setChecked(True)
        proc_row.addWidget(self.chk_skip)
        proc_row.addStretch()
        proc_row.addWidget(QLabel("Workers:"))
        self.spin_workers = QSpinBox()
        self.spin_workers.setRange(1, 16)
        self.spin_workers.setValue(8)
        proc_row.addWidget(self.spin_workers)
        dl.addLayout(proc_row)

        # Form layout for defaults
        form = QFormLayout()
        form.setLabelAlignment(Qt.AlignRight)
        form.setHorizontalSpacing(12)
        form.setVerticalSpacing(8)

        self.combo_smooth2d = QComboBox()
        self.combo_smooth2d.setMinimumWidth(_MIN_COMBO_W)
        s2d_presets = presets.get("smooth_2d", [])
        for p in s2d_presets:
            self.combo_smooth2d.addItem(p.get("label", "Raw"))
        s2d_def = defaults.get("smooth_2d", {})
        for i, p in enumerate(s2d_presets):
            if {k: v for k, v in p.items() if k != "label"} == s2d_def:
                self.combo_smooth2d.setCurrentIndex(i)
        form.addRow("2D Smoothing:", self.combo_smooth2d)

        self.combo_smooth_drift = QComboBox()
        self.combo_smooth_drift.setMinimumWidth(_MIN_COMBO_W)
        sd_presets = presets.get("smooth_drift", [])
        for p in sd_presets:
            self.combo_smooth_drift.addItem(p.get("label", "Raw"))
        sd_def = defaults.get("smooth_drift", {})
        for i, p in enumerate(sd_presets):
            if {k: v for k, v in p.items() if k != "label"} == sd_def:
                self.combo_smooth_drift.setCurrentIndex(i)
        form.addRow("Drift Smoothing:", self.combo_smooth_drift)

        self.combo_cmap = QComboBox()
        self.combo_cmap.setMinimumWidth(_MIN_COMBO_W)
        cmap_names = []
        for c in presets.get("colormap", ["Viridis"]):
            cmap_names.append(c if isinstance(c, str) else c.get("name", "?"))
        self.combo_cmap.addItems(cmap_names)
        self.combo_cmap.setCurrentText(defaults.get("colormap", "Viridis"))
        form.addRow("Colormap:", self.combo_cmap)

        self.combo_scale = QComboBox()
        self.combo_scale.setMinimumWidth(_MIN_COMBO_W)
        self.combo_scale.addItems([str(s) for s in presets.get("scale", ["1"])])
        self.combo_scale.setCurrentText(str(defaults.get("scale", "1")))
        form.addRow("Intensity Scale:", self.combo_scale)

        self.spin_mzbins = QSpinBox()
        self.spin_mzbins.setRange(100, 25600)
        self.spin_mzbins.setSingleStep(200)
        self.spin_mzbins.setValue(defaults.get("mz_bins", 3200))
        self.spin_mzbins.setMinimumWidth(_MIN_COMBO_W)
        form.addRow("m/z Bins:", self.spin_mzbins)

        # Noise thresholds
        self.spin_noise2d = QDoubleSpinBox()
        self.spin_noise2d.setRange(0, 10); self.spin_noise2d.setSingleStep(0.1)
        self.spin_noise2d.setDecimals(2); self.spin_noise2d.setSuffix(" %")
        self.spin_noise2d.setValue(defaults.get("noise_2d", 0.5))
        self.spin_noise2d.setMinimumWidth(_MIN_COMBO_W)
        form.addRow("Noise 2D:", self.spin_noise2d)

        self.spin_noise_drift = QDoubleSpinBox()
        self.spin_noise_drift.setRange(0, 10); self.spin_noise_drift.setSingleStep(0.5)
        self.spin_noise_drift.setDecimals(2); self.spin_noise_drift.setSuffix(" %")
        self.spin_noise_drift.setValue(defaults.get("noise_drift", 0))
        self.spin_noise_drift.setMinimumWidth(_MIN_COMBO_W)
        form.addRow("Noise Drift:", self.spin_noise_drift)

        self.spin_noise_mz = QDoubleSpinBox()
        self.spin_noise_mz.setRange(0, 10); self.spin_noise_mz.setSingleStep(0.5)
        self.spin_noise_mz.setDecimals(2); self.spin_noise_mz.setSuffix(" %")
        self.spin_noise_mz.setValue(defaults.get("noise_mz", 0))
        self.spin_noise_mz.setMinimumWidth(_MIN_COMBO_W)
        form.addRow("Noise m/z:", self.spin_noise_mz)

        # Figure dimensions
        fig_label = QLabel("Figure dimensions:")
        fig_label.setProperty("subtitle", True)
        form.addRow(fig_label)

        self.spin_width = QSpinBox()
        self.spin_width.setRange(600, 2400); self.spin_width.setSingleStep(100)
        self.spin_width.setValue(fig_cfg.get("width", 1190)); self.spin_width.setSuffix(" px")
        self.spin_width.setMinimumWidth(_MIN_COMBO_W)
        form.addRow("Width:", self.spin_width)

        self.spin_height = QSpinBox()
        self.spin_height.setRange(400, 1600); self.spin_height.setSingleStep(100)
        self.spin_height.setValue(fig_cfg.get("height", 680)); self.spin_height.setSuffix(" px")
        self.spin_height.setMinimumWidth(_MIN_COMBO_W)
        form.addRow("Height:", self.spin_height)

        self.spin_font = QSpinBox()
        self.spin_font.setRange(10, 32)
        self.spin_font.setValue(fig_cfg.get("font_size", 20)); self.spin_font.setSuffix(" pt")
        self.spin_font.setMinimumWidth(_MIN_COMBO_W)
        form.addRow("Font Size:", self.spin_font)

        dl.addLayout(form)
        cfg_tabs.addTab(defaults_w, "Defaults")

        # ==== Tab 2: Presets ====
        presets_w = QWidget()
        pl = QVBoxLayout(presets_w)
        pl.setContentsMargins(8, 8, 8, 8)
        pl_sub = QLabel("Dropdown options for each HTML control. Edit JSON to add or remove presets.")
        pl_sub.setProperty("subtitle", True)
        pl.addWidget(pl_sub)

        self.presets_editor = QTextEdit()
        mono = QFont("Consolas", 10)
        mono.setStyleHint(QFont.Monospace)
        self.presets_editor.setFont(mono)
        self.presets_editor.setTabStopDistance(24)
        self.presets_editor.setLineWrapMode(QTextEdit.NoWrap)
        self.presets_editor.setMinimumHeight(200)
        self.presets_editor.setPlainText(json.dumps(presets, indent=2, ensure_ascii=False))
        pl.addWidget(self.presets_editor, stretch=1)

        btn_row = QHBoxLayout()
        btn_format = QPushButton("Format JSON")
        btn_format.setProperty("secondary", True)
        btn_format.clicked.connect(self._format_json)
        btn_row.addWidget(btn_format)
        btn_reset = QPushButton("Reset to Defaults")
        btn_reset.setProperty("secondary", True)
        btn_reset.clicked.connect(self._reset_presets)
        btn_row.addWidget(btn_reset)
        btn_row.addStretch()
        pl.addLayout(btn_row)

        cfg_tabs.addTab(presets_w, "Presets")

        cfg_layout.addWidget(cfg_tabs)

        # Save button
        save_row = QHBoxLayout()
        save_row.addStretch()
        btn_save = QPushButton("Save Config")
        btn_save.setProperty("secondary", True)
        btn_save.clicked.connect(self._save_config)
        save_row.addWidget(btn_save)
        cfg_layout.addLayout(save_row)

        layout.addWidget(cfg_grp, stretch=1)

        # ---- Run ----
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

    def _format_json(self):
        try:
            obj = json.loads(self.presets_editor.toPlainText())
            self.presets_editor.setPlainText(json.dumps(obj, indent=2, ensure_ascii=False))
            self._log.log("JSON formatted.", "info")
        except json.JSONDecodeError as e:
            self._log.log(f"Invalid JSON: {e}", "error")

    def _reset_presets(self):
        default_cfg = _load_config()
        presets = default_cfg.get("presets", {})
        self.presets_editor.setPlainText(json.dumps(presets, indent=2, ensure_ascii=False))
        self._log.log("Presets reset to saved config.", "info")

    def _save_config(self):
        try:
            presets = json.loads(self.presets_editor.toPlainText())
        except json.JSONDecodeError as e:
            self._log.log(f"Invalid presets JSON: {e}", "error")
            return

        s2d_presets = presets.get("smooth_2d", [])
        sd_presets = presets.get("smooth_drift", [])
        s2d_idx = self.combo_smooth2d.currentIndex()
        sd_idx = self.combo_smooth_drift.currentIndex()
        s2d = {k: v for k, v in s2d_presets[s2d_idx].items() if k != "label"} if s2d_idx < len(s2d_presets) else {"method": "raw"}
        sd = {k: v for k, v in sd_presets[sd_idx].items() if k != "label"} if sd_idx < len(sd_presets) else {"method": "raw"}

        self._cfg["defaults"] = {
            "smooth_2d": s2d, "scale": self.combo_scale.currentText(),
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
        self._cfg["presets"] = presets
        _save_config(self._cfg)
        self._log.log("Config saved.", "success")

    def _run(self):
        data = self.pick_data.path()
        out = self.pick_out.path()
        if not data or not out:
            self._log.log("Set text directory and output directory.", "warning")
            return
        self._save_config()
        raw = self.pick_raw.path() or None
        skip = self.chk_skip.isChecked()
        workers = self.spin_workers.value()

        self.btn_run.setEnabled(False)
        self.progress.setVisible(True)
        self.progress.setRange(0, 0)
        self._log.log("Generating HTML plots...", "accent")

        self._worker = Worker(self._do_plot, data, out, raw, skip, workers)
        self._worker.finished.connect(self._on_done)
        self._worker.error.connect(self._on_error)
        self._worker.start()

    def _do_plot(self, data, out, raw, skip, workers):
        data_p = Path(data)
        im_files = list(data_p.glob("*_im.txt"))
        self._log._sig.message.emit(f"Found {len(im_files)} IM files", "info")
        if not im_files:
            raise FileNotFoundError(f"No _im.txt files found in {data}")
        out_p = Path(out)
        for stale in out_p.glob("*_2d_imms.html"):
            if stale.stat().st_size == 0:
                stale.unlink()
        from deconvovo.imms_plot import run as plot_run
        results = plot_run(data_p, out_p, skip_existing=skip,
                           raw_dir=Path(raw) if raw else None, n_workers=workers)
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
        n = sum(1 for h in out.glob("*_2d_imms.html") if h.stat().st_size > 0)
        self._log.log(f"Done: {n} HTML files", "success")

    def _on_error(self, msg):
        self.btn_run.setEnabled(True)
        self.progress.setVisible(False)
        self._log.log(f"HTML generation failed: {msg}", "error")
