"""2D IM-MS HTML viewer generation panel."""
from __future__ import annotations

from pathlib import Path

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QGroupBox,
    QPushButton, QProgressBar, QCheckBox, QSpinBox,
)

from gui.widgets.file_picker import DirPicker
from gui.widgets.log_panel import LogPanel
from gui.widgets.worker import Worker


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
        sub = QLabel("Generate interactive HTML heatmaps with linked drift and m/z panels. "
                     "Opens in your web browser.")
        sub.setProperty("subtitle", True)
        sub.setWordWrap(True)
        layout.addWidget(sub)

        grp = QGroupBox("Paths")
        gl = QVBoxLayout(grp)
        self.pick_data = DirPicker("Text data:", dialog_title="Select directory with _ms.txt / _im.txt")
        self.pick_raw = DirPicker("Raw data:", dialog_title="Select directory with .raw folders (for pusher period)")
        self.pick_out = DirPicker("Output:", dialog_title="Select output directory for HTML files")
        gl.addWidget(self.pick_data)
        gl.addWidget(self.pick_raw)
        gl.addWidget(self.pick_out)
        layout.addWidget(grp)

        opts = QGroupBox("Options")
        ol = QHBoxLayout(opts)
        self.chk_skip = QCheckBox("Skip existing")
        self.chk_skip.setChecked(True)
        ol.addWidget(self.chk_skip)
        ol.addStretch()
        ol.addWidget(QLabel("Workers:"))
        self.spin_workers = QSpinBox()
        self.spin_workers.setRange(1, 16)
        self.spin_workers.setValue(8)
        ol.addWidget(self.spin_workers)
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
            self._log.log("Set data and output directories.", "warning")
            return

        raw = self.pick_raw.path() or None
        skip = self.chk_skip.isChecked()
        workers = self.spin_workers.value()

        self.btn_run.setEnabled(False)
        self.progress.setVisible(True)
        self.progress.setRange(0, 0)
        self._log.log(f"Generating HTML plots from {data}...", "accent")

        self._worker = Worker(self._do_plot, data, out, raw, skip, workers)
        self._worker.finished.connect(self._on_done)
        self._worker.error.connect(self._on_error)
        self._worker.start()

    def _do_plot(self, data, out, raw, skip, workers):
        from deconvovo.imms_plot import run as plot_run
        plot_run(
            Path(data), Path(out),
            skip_existing=skip,
            raw_dir=Path(raw) if raw else None,
            n_workers=workers,
        )

    def _on_done(self, _):
        self.btn_run.setEnabled(True)
        self.progress.setVisible(False)
        self._log.log("HTML generation complete.", "success")

    def _on_error(self, msg):
        self.btn_run.setEnabled(True)
        self.progress.setVisible(False)
        self._log.log(f"HTML generation failed: {msg}", "error")
