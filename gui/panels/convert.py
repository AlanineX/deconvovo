"""Convert .raw to text panel."""
from __future__ import annotations

from pathlib import Path

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QGroupBox,
    QPushButton, QProgressBar, QCheckBox, QSpinBox,
)

from gui.widgets.file_picker import DirPicker
from gui.widgets.log_panel import LogPanel
from gui.widgets.worker import Worker


class ConvertPanel(QWidget):
    def __init__(self, log_panel: LogPanel, parent=None):
        super().__init__(parent)
        self._log = log_panel
        self._worker = None

        layout = QVBoxLayout(self)
        layout.setContentsMargins(16, 16, 16, 16)
        layout.setSpacing(12)

        title = QLabel("Convert Waters .raw to Text")
        title.setProperty("title", True)
        layout.addWidget(title)
        sub = QLabel("Extract MS and IM data from Waters MassLynx .raw directories using CDCReader.")
        sub.setProperty("subtitle", True)
        sub.setWordWrap(True)
        layout.addWidget(sub)

        grp = QGroupBox("Paths")
        gl = QVBoxLayout(grp)
        self.pick_input = DirPicker("Input dir:",
            dialog_title="Select directory containing .raw folders",
            hint="Directory containing one or more .raw folders from MassLynx acquisition")
        self.pick_output = DirPicker("Output dir:",
            dialog_title="Select output directory",
            hint="Directory where converted _ms.txt and _im.txt files will be written")
        gl.addWidget(self.pick_input)
        gl.addWidget(self.pick_output)
        layout.addWidget(grp)

        grp2 = QGroupBox("Options")
        ol = QHBoxLayout(grp2)
        self.chk_skip = QCheckBox("Skip existing")
        self.chk_skip.setChecked(True)
        ol.addWidget(self.chk_skip)
        ol.addStretch()
        ol.addWidget(QLabel("Workers:"))
        self.spin_workers = QSpinBox()
        self.spin_workers.setRange(1, 16)
        self.spin_workers.setValue(4)
        ol.addWidget(self.spin_workers)
        layout.addWidget(grp2)

        run_row = QHBoxLayout()
        run_row.addStretch()
        self.btn_run = QPushButton("Convert")
        self.btn_run.setFixedWidth(160)
        self.btn_run.setFixedHeight(36)
        self.btn_run.clicked.connect(self._run)
        run_row.addWidget(self.btn_run)
        layout.addLayout(run_row)

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
        workers = self.spin_workers.value()
        self.btn_run.setEnabled(False)
        self.progress.setVisible(True)
        self.progress.setRange(0, 0)
        self._log.log(f"Converting {inp} (workers={workers})", "accent")
        self._worker = Worker(self._do_convert, inp, out, workers)
        self._worker.finished.connect(self._on_done)
        self._worker.error.connect(self._on_error)
        self._worker.start()

    def _do_convert(self, inp, out, workers):
        from deconvovo.imms_convert import run as convert_run
        return convert_run(Path(inp), Path(out), n_workers=workers)

    def _on_done(self, result):
        self.btn_run.setEnabled(True)
        self.progress.setVisible(False)
        self._log.log("Conversion complete.", "success")

    def _on_error(self, msg):
        self.btn_run.setEnabled(True)
        self.progress.setVisible(False)
        self._log.log(f"Conversion failed: {msg}", "error")
