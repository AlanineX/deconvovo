"""2D IM-MS HTML viewer generation panel."""
from __future__ import annotations

import json
import sys
from copy import deepcopy
from pathlib import Path

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QGroupBox,
    QPushButton, QProgressBar, QCheckBox, QSpinBox, QComboBox,
    QDoubleSpinBox, QTextEdit, QTabWidget, QFormLayout, QScrollArea, QFileDialog,
)
from PySide6.QtCore import Qt, QUrl, QStandardPaths
from PySide6.QtGui import QFont, QDesktopServices

from gui.widgets.file_picker import DirPicker
from gui.widgets.log_panel import LogPanel
from gui.widgets.worker import Worker

_MIN_COMBO_W = 200


def _config_path():
    if getattr(sys, 'frozen', False):
        return Path(QStandardPaths.writableLocation(QStandardPaths.AppConfigLocation)) / "imms_plot_config.json"
    return Path(__file__).parent.parent.parent / "config" / "imms_plot_config.json"


def _load_config():
    p = _config_path()
    if not p.exists() and getattr(sys, 'frozen', False):
        p = Path(sys._MEIPASS) / "config" / "imms_plot_config.json"
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
        sub = QLabel("Create offline HTML viewers with linked drift and m/z panels. Open them in your browser.")
        sub.setProperty("subtitle", True)
        sub.setWordWrap(True)
        layout.addWidget(sub)

        # ---- Input ----
        grp_in = QGroupBox("Input")
        gl_in = QVBoxLayout(grp_in)
        self.pick_data = DirPicker("Text dir:",
            dialog_title="Select directory with _ms.txt / _im.txt files",
            hint="Directory with converted text files (_ms.txt and _im.txt)")
        self.pick_raw = DirPicker("Raw dir:",
            dialog_title="Select directory with .raw folders",
            hint="Reads pusher period from .raw to convert drift bins to time (ms). Without this, y-axis shows bin numbers.")
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
        self.config_tabs = cfg_tabs

        # ==== Tab 1: Defaults ====
        defaults_w = QWidget()
        dl = QVBoxLayout(defaults_w)
        dl.setContentsMargins(8, 8, 8, 8)
        dl_sub = QLabel("Initial viewer settings. Users can change them interactively in the HTML.")
        dl_sub.setProperty("subtitle", True)
        dl.addWidget(dl_sub)

        # Processing row
        proc_row = QHBoxLayout()
        self.chk_skip = QCheckBox("Reuse existing HTML")
        self.chk_skip.setChecked(False)
        self.chk_skip.setToolTip("Reuses files only when they were generated with these settings. Leave off to refresh older viewers.")
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
        if not any({k: v for k, v in p.items() if k != "label"} == s2d_def for p in s2d_presets):
            self.combo_smooth2d.addItem("Custom (saved)")
            self.combo_smooth2d.setCurrentIndex(len(s2d_presets))
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
        if not any({k: v for k, v in p.items() if k != "label"} == sd_def for p in sd_presets):
            self.combo_smooth_drift.addItem("Custom (saved)")
            self.combo_smooth_drift.setCurrentIndex(len(sd_presets))
        form.addRow("Drift Smoothing:", self.combo_smooth_drift)

        self.combo_cmap = QComboBox()
        self.combo_cmap.setMinimumWidth(_MIN_COMBO_W)
        cmap_names = []
        for c in presets.get("colormap", ["Viridis"]):
            cmap_names.append(c if isinstance(c, str) else c.get("name", "?"))
        self.combo_cmap.addItems(cmap_names)
        if defaults.get("colormap", "Viridis") not in cmap_names:
            self.combo_cmap.addItem(defaults.get("colormap", "Viridis"))
        self.combo_cmap.setCurrentText(defaults.get("colormap", "Viridis"))
        form.addRow("Colormap:", self.combo_cmap)

        self.combo_scale = QComboBox()
        self.combo_scale.setMinimumWidth(_MIN_COMBO_W)
        self.combo_scale.addItems([str(s) for s in presets.get("scale", ["1"])])
        self.combo_scale.setCurrentText(str(defaults.get("scale", "1")))
        form.addRow("m/z smoothing scale:", self.combo_scale)
        self.combo_scale.setToolTip("Multiplier for m/z-axis 2D smoothing. Lower values make peaks thinner; higher values make them wider.")

        self.combo_zscale = QComboBox()
        self.combo_zscale.setMinimumWidth(_MIN_COMBO_W)
        self.combo_zscale.addItems([str(s) for s in presets.get("z_scale", ["linear", "sqrt", "cbrt", "p4", "p5", "log2", "log10"])])
        self.combo_zscale.setCurrentText(str(defaults.get("z_scale", "sqrt")))
        form.addRow("Heatmap intensity:", self.combo_zscale)

        self.combo_ms_mode = QComboBox()
        self.combo_ms_mode.addItem("Profile", "profile")
        self.combo_ms_mode.addItem("Sticks (peak picked)", "stick")
        self.combo_ms_mode.setCurrentIndex(max(0, self.combo_ms_mode.findData(defaults.get("ms_mode", "profile"))))
        form.addRow("m/z display:", self.combo_ms_mode)

        dl.addLayout(form)
        dl.addStretch()
        cfg_tabs.addTab(defaults_w, "Display")

        # Advanced settings stay available without crowding the main display choices.
        advanced_w = QWidget()
        advanced_layout = QVBoxLayout(advanced_w)
        form = QFormLayout()
        form.setLabelAlignment(Qt.AlignRight)
        form.setHorizontalSpacing(12)
        form.setVerticalSpacing(8)

        self.spin_mzbins = QSpinBox()
        self.spin_mzbins.setRange(0, 25600)
        self.spin_mzbins.setSpecialValueText("Native")
        self.spin_mzbins.setSingleStep(200)
        mz_bins = defaults.get("mz_bins", 3200)
        self.spin_mzbins.setValue(0 if str(mz_bins).lower() == "native" else int(mz_bins))
        self.spin_mzbins.setMinimumWidth(_MIN_COMBO_W)
        self.spin_mzbins.setToolTip("Native keeps the available m/z sampling (capped at 25,600 bins). Lower counts reduce HTML size.")
        form.addRow("m/z Bins:", self.spin_mzbins)

        self.spin_driftbins = QSpinBox()
        self.spin_driftbins.setRange(0, 25600)
        self.spin_driftbins.setSpecialValueText("Native")
        drift_bins = defaults.get("drift_bins", "native")
        self.spin_driftbins.setValue(0 if str(drift_bins).lower() == "native" else int(drift_bins))
        form.addRow("Drift bins:", self.spin_driftbins)

        self.combo_render = QComboBox()
        for label, value in [("Crisp", "off"), ("Smooth (fast)", "fast"), ("Smooth (best)", "best")]:
            self.combo_render.addItem(label, value)
        self.combo_render.setCurrentIndex(max(0, self.combo_render.findData(defaults.get("render", "off"))))
        form.addRow("Heatmap rendering:", self.combo_render)

        # Noise thresholds
        self.spin_noise2d = QDoubleSpinBox()
        self.spin_noise2d.setRange(0, 100); self.spin_noise2d.setSingleStep(0.1)
        self.spin_noise2d.setDecimals(2); self.spin_noise2d.setSuffix(" %")
        self.spin_noise2d.setValue(defaults.get("noise_2d", 0.5))
        self.spin_noise2d.setMinimumWidth(_MIN_COMBO_W)
        form.addRow("Noise 2D:", self.spin_noise2d)

        self.spin_denoise = QDoubleSpinBox()
        self.spin_denoise.setRange(0, 100)
        self.spin_denoise.setDecimals(2)
        self.spin_denoise.setSpecialValueText("Off (hard cutoff)")
        self.spin_denoise.setSuffix(" %")
        denoise = defaults.get("denoise", 5)
        self.spin_denoise.setValue(0 if str(denoise) == "off" else float(denoise))
        self.spin_denoise.setToolTip("Softens the 2D noise cutoff. Requires a positive Noise 2D threshold.")
        form.addRow("Cutoff fade:", self.spin_denoise)

        self.spin_noise_drift = QDoubleSpinBox()
        self.spin_noise_drift.setRange(0, 100); self.spin_noise_drift.setSingleStep(0.5)
        self.spin_noise_drift.setDecimals(2); self.spin_noise_drift.setSuffix(" %")
        self.spin_noise_drift.setValue(defaults.get("noise_drift", 0))
        self.spin_noise_drift.setMinimumWidth(_MIN_COMBO_W)
        form.addRow("Noise Drift:", self.spin_noise_drift)

        self.spin_noise_mz = QDoubleSpinBox()
        self.spin_noise_mz.setRange(0, 100); self.spin_noise_mz.setSingleStep(0.5)
        self.spin_noise_mz.setDecimals(2); self.spin_noise_mz.setSuffix(" %")
        self.spin_noise_mz.setValue(defaults.get("noise_mz", 0))
        self.spin_noise_mz.setMinimumWidth(_MIN_COMBO_W)
        form.addRow("Noise m/z:", self.spin_noise_mz)

        initial = self._cfg.get("initial_view", {}) or {}
        self.chk_mz_range = QCheckBox("Set initial m/z window")
        self.chk_drift_range = QCheckBox("Set initial drift window")
        self.spin_mz_min = QDoubleSpinBox()
        self.spin_mz_max = QDoubleSpinBox()
        self.spin_drift_min = QDoubleSpinBox()
        self.spin_drift_max = QDoubleSpinBox()
        for check, lo, hi, key, fallback in [
            (self.chk_mz_range, self.spin_mz_min, self.spin_mz_max, "mz_range", (0, 10000)),
            (self.chk_drift_range, self.spin_drift_min, self.spin_drift_max, "drift_range", (0, 100)),
        ]:
            enabled = bool(initial.get(key))
            bounds = initial.get(key) or fallback
            check.setChecked(enabled)
            for spin, value in zip((lo, hi), bounds):
                spin.setRange(0, 1000000)
                spin.setDecimals(4)
                spin.setValue(float(value))
                spin.setEnabled(enabled)
                check.toggled.connect(spin.setEnabled)
            row = QHBoxLayout()
            row.addWidget(lo)
            row.addWidget(QLabel("to"))
            row.addWidget(hi)
            form.addRow(check, row)
        self.chk_drift_range.setToolTip("Milliseconds when a RAW pusher period is available; otherwise drift bin numbers.")

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

        advanced_layout.addLayout(form)
        advanced_layout.addStretch()
        advanced_scroll = QScrollArea()
        advanced_scroll.setWidgetResizable(True)
        advanced_scroll.setFrameShape(QScrollArea.NoFrame)
        advanced_scroll.setWidget(advanced_w)
        cfg_tabs.addTab(advanced_scroll, "Advanced")

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

        cfg_tabs.addTab(presets_w, "Preset JSON")

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
        self.btn_open_html = QPushButton("Open HTML…")
        self.btn_open_html.setProperty("secondary", True)
        self.btn_open_html.clicked.connect(self._open_html)
        run_row.addWidget(self.btn_open_html)
        self.btn_open_out = QPushButton("Open output folder")
        self.btn_open_out.setProperty("secondary", True)
        self.btn_open_out.clicked.connect(self._open_output)
        run_row.addWidget(self.btn_open_out)
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
        self.set_config(self._cfg)

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

    def current_config(self):
        """Collect the displayed settings without writing the application defaults."""
        try:
            presets = json.loads(self.presets_editor.toPlainText())
        except json.JSONDecodeError as e:
            raise ValueError(f"Invalid presets JSON: {e}") from e
        if not isinstance(presets, dict):
            raise ValueError("Preset JSON must contain an object of dropdown options.")
        cfg = deepcopy(self._cfg)
        cfg.setdefault("defaults", {}).update({
            "smooth_2d": deepcopy(self.combo_smooth2d.currentData()),
            "smooth_drift": deepcopy(self.combo_smooth_drift.currentData()),
            "scale": self.combo_scale.currentText(),
            "noise_2d": self.spin_noise2d.value(),
            "noise_drift": self.spin_noise_drift.value(),
            "noise_mz": self.spin_noise_mz.value(),
            "mz_bins": self.spin_mzbins.value() or "native",
            "drift_bins": self.spin_driftbins.value() or "native",
            "colormap": self.combo_cmap.currentText(),
            "z_scale": self.combo_zscale.currentText(),
            "render": self.combo_render.currentData(),
            "ms_mode": self.combo_ms_mode.currentData(),
            "denoise": self.spin_denoise.value() or "off",
        })
        cfg.setdefault("figure", {}).update({
            "font_size": self.spin_font.value(),
            "width": self.spin_width.value(),
            "height": self.spin_height.value(),
        })
        cfg.setdefault("presets", {}).update(presets)
        if not isinstance(cfg.get("initial_view"), dict):
            cfg["initial_view"] = {}
        for check, lo, hi, key in [
            (self.chk_mz_range, self.spin_mz_min, self.spin_mz_max, "mz_range"),
            (self.chk_drift_range, self.spin_drift_min, self.spin_drift_max, "drift_range"),
        ]:
            if check.isChecked() and lo.value() >= hi.value():
                raise ValueError(f"The {key.replace('_', ' ')} minimum must be lower than its maximum.")
            cfg.setdefault("initial_view", {})[key] = [lo.value(), hi.value()] if check.isChecked() else None
        return cfg

    def set_config(self, cfg):
        """Load workspace/preset settings into the existing editor widgets."""
        self._cfg = deepcopy(cfg)
        defaults = cfg.get("defaults", {})
        presets = cfg.get("presets", {})
        for combo, key in [(self.combo_smooth2d, "smooth_2d"), (self.combo_smooth_drift, "smooth_drift")]:
            combo.clear()
            default = defaults.get(key, {"method": "raw"})
            for preset in presets.get(key, []):
                value = {k: v for k, v in preset.items() if k != "label"}
                combo.addItem(preset.get("label", str(value)), value)
            index = combo.findData(default)
            if index < 0:
                combo.addItem("Custom (saved)", default)
                index = combo.count() - 1
            combo.setCurrentIndex(index)
        for combo, key, fallback in [
            (self.combo_scale, "scale", ["1"]),
            (self.combo_zscale, "z_scale", ["linear", "sqrt", "cbrt", "p4", "p5", "log2", "log10"]),
            (self.combo_cmap, "colormap", ["Viridis"]),
        ]:
            combo.clear()
            options = presets.get(key, fallback)
            names = [str(value.get("name", "?")) if isinstance(value, dict) else str(value) for value in options]
            default = str(defaults.get(key, {"scale": "1", "z_scale": "sqrt", "colormap": "Viridis"}[key]))
            combo.addItems(names)
            if default not in names:
                combo.addItem(default)
            combo.setCurrentText(default)
        for combo, key, fallback in [(self.combo_render, "render", "off"), (self.combo_ms_mode, "ms_mode", "profile")]:
            combo.setCurrentIndex(max(0, combo.findData(defaults.get(key, fallback))))
        for spin, key, fallback in [
            (self.spin_mzbins, "mz_bins", 3200), (self.spin_driftbins, "drift_bins", "native"),
        ]:
            value = defaults.get(key, fallback)
            spin.setValue(0 if str(value).lower() == "native" else int(value))
        for spin, key, fallback in [
            (self.spin_noise2d, "noise_2d", 0.5), (self.spin_noise_drift, "noise_drift", 0),
            (self.spin_noise_mz, "noise_mz", 0), (self.spin_denoise, "denoise", 5),
        ]:
            value = defaults.get(key, fallback)
            spin.setValue(0 if str(value) == "off" else float(value))
        figure = cfg.get("figure", {})
        for spin, key, fallback in [(self.spin_font, "font_size", 20), (self.spin_width, "width", 1190), (self.spin_height, "height", 680)]:
            spin.setValue(figure.get(key, fallback))
        initial = cfg.get("initial_view", {}) or {}
        for check, lo, hi, key, fallback in [
            (self.chk_mz_range, self.spin_mz_min, self.spin_mz_max, "mz_range", (0, 10000)),
            (self.chk_drift_range, self.spin_drift_min, self.spin_drift_max, "drift_range", (0, 100)),
        ]:
            bounds = initial.get(key) or fallback
            lo.setValue(float(bounds[0]))
            hi.setValue(float(bounds[1]))
            check.setChecked(bool(initial.get(key)))
        self.presets_editor.setPlainText(json.dumps(presets, indent=2, ensure_ascii=False))

    def _save_config(self):
        try:
            cfg = self.current_config()
            _save_config(cfg)
        except (ValueError, OSError) as e:
            self._log.log(str(e), "error")
            return False
        self._cfg = cfg
        self._log.log("Config saved.", "success")
        return True

    def _run(self):
        data = self.pick_data.path()
        out = self.pick_out.path()
        if not data or not out:
            self._log.log("Set text directory and output directory.", "warning")
            return
        if not self._save_config():
            return
        raw = self.pick_raw.path() or None
        skip = self.chk_skip.isChecked()
        workers = self.spin_workers.value()

        self.btn_run.setEnabled(False)
        self.progress.setVisible(True)
        self.progress.setRange(0, 0)
        self._log.log("Generating HTML plots...", "accent")

        self._worker = Worker(self._do_plot, data, out, raw, skip, workers, deepcopy(self._cfg))
        self._worker.finished.connect(self._on_done)
        self._worker.error.connect(self._on_error)
        self._worker.start()

    def _do_plot(self, data, out, raw, skip, workers, config=None):
        data_p = Path(data)
        im_files = list(data_p.glob("*_im.txt"))
        self._log._sig.message.emit(f"Found {len(im_files)} IM files", "info")
        if not im_files:
            raise FileNotFoundError(f"No _im.txt files found in {data}")
        out_p = Path(out)
        out_p.mkdir(parents=True, exist_ok=True)
        config = config if config is not None else deepcopy(self._cfg)
        used_config_path = out_p / "viewer_config.used.json"
        if skip:
            try:
                previous = json.loads(used_config_path.read_text(encoding="utf-8"))
            except (OSError, ValueError):
                previous = None
            if previous != config:
                skip = False
                self._log._sig.message.emit("Viewer settings changed or previous settings are unknown; regenerating existing HTML.", "info")
        config_path = out_p / "viewer_config.json"
        config_path.write_text(json.dumps(config, indent=2, ensure_ascii=False), encoding="utf-8")
        for stale in out_p.glob("*_2d_imms.html"):
            if stale.stat().st_size == 0:
                stale.unlink()
        from deconvovo.imms_plot import run as plot_run
        results = plot_run(data_p, out_p, skip_existing=skip,
                           raw_dir=Path(raw) if raw else None, n_workers=workers,
                           config_path=config_path)
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
        used_config_path.write_text(json.dumps(config, indent=2, ensure_ascii=False), encoding="utf-8")
        return results

    def _open_html(self):
        output = self.pick_out.path()
        if not output:
            self._log.log("Choose an output directory first.", "warning")
            return
        viewers = sorted(Path(output).glob("*_2d_imms.html"))
        if not viewers:
            self._log.log("Generate an HTML viewer first.", "warning")
            return
        if len(viewers) == 1:
            path = str(viewers[0])
        else:
            path, _ = QFileDialog.getOpenFileName(self, "Open generated HTML viewer", output, "HTML viewers (*.html)")
        if path:
            QDesktopServices.openUrl(QUrl.fromLocalFile(str(Path(path).resolve())))

    def _open_output(self):
        output = self.pick_out.path()
        if output and Path(output).is_dir():
            QDesktopServices.openUrl(QUrl.fromLocalFile(str(Path(output).resolve())))
        else:
            self._log.log("Choose an existing output directory first.", "warning")

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
