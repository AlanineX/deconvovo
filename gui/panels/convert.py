"""Convert .raw → text panel."""
from __future__ import annotations

from pathlib import Path

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QGroupBox,
    QPushButton, QProgressBar, QCheckBox, QSpinBox,
)
from PySide6.QtCore import Qt

from gui.widgets.file_picker import DirPicker
from gui.widgets.log_panel import LogPanel
from gui.widgets.worker import Worker
from gui.theme import TEXT_BRIGHT


class ConvertPanel(QWidget):
    def __init__(self, log_panel: LogPanel, parent=None):
        super().__init__(parent)
        self._log = log_panel
        self._worker = None

        layout = QVBoxLayout(self)
        layout.setContentsMargins(16, 16, 16, 16)
        layout.setSpacing(12)

        # Title
        title = QLabel("Convert Waters .raw → Text")
        title.setProperty("title", True)
        layout.addWidget(title)
        sub = QLabel("Extract MS and IM data from Waters MassLynx .raw directories using CDCReader.")
        sub.setProperty("subtitle", True)
        sub.setWordWrap(True)
        layout.addWidget(sub)

        # Input / Output
        grp = QGroupBox("Paths")
        gl = QVBoxLayout(grp)
        self.pick_input = DirPicker("Input (.raw):", dialog_title="Select directory containing .raw folders")
        self.pick_output = DirPicker("Output:", dialog_title="Select output directory")
        gl.addWidget(self.pick_input)
        gl.addWidget(self.pick_output)
        layout.addWidget(grp)

        # Options
        grp2 = QGroupBox("Options")
        ol = QHBoxLayout(grp2)
        self.chk_skip = QCheckBox("Skip existing")
        self.chk_skip.setChecked(True)
        ol.addWidget(self.chk_skip)
        self.chk_skip_ms = QCheckBox("Skip MS")
        ol.addWidget(self.chk_skip_ms)
        self.chk_skip_im = QCheckBox("Skip IM")
        ol.addWidget(self.chk_skip_im)
        ol.addStretch()
        ol.addWidget(QLabel("Function:"))
        self.spin_func = QSpinBox()
        self.spin_func.setRange(1, 10)
        self.spin_func.setValue(1)
        ol.addWidget(self.spin_func)
        layout.addWidget(grp2)

        # Run
        run_row = QHBoxLayout()
        run_row.addStretch()
        self.btn_run = QPushButton("Convert")
        self.btn_run.setFixedWidth(160)
        self.btn_run.setFixedHeight(36)
        self.btn_run.clicked.connect(self._run)
        run_row.addWidget(self.btn_run)
        layout.addLayout(run_row)

        # Progress
        self.progress = QProgressBar()
        self.progress.setVisible(False)
        layout.addWidget(self.progress)

        layout.addStretch()

    def _run(self):
        inp = self.pick_input.path()
        out = self.pick_output.path()
        if not inp or not out:
            self._log.log("Set input and output directories first.", "warning")
            return

        self.btn_run.setEnabled(False)
        self.progress.setVisible(True)
        self.progress.setRange(0, 0)  # indeterminate
        self._log.log(f"Converting {inp} → {out}", "accent")

        self._worker = Worker(self._do_convert, inp, out)
        self._worker.finished.connect(self._on_done)
        self._worker.error.connect(self._on_error)
        self._worker.start()

    def _do_convert(self, inp: str, out: str):
        from deconvovo.imms_convert import run as convert_run
        return convert_run(Path(inp), Path(out))

    def _on_done(self, result):
        self.btn_run.setEnabled(True)
        self.progress.setVisible(False)
        self._log.log("Conversion complete.", "success")

    def _on_error(self, msg):
        self.btn_run.setEnabled(True)
        self.progress.setVisible(False)
        self._log.log(f"Conversion failed: {msg}", "error")
