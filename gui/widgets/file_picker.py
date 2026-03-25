"""Directory and file picker widgets."""
from __future__ import annotations

from PySide6.QtWidgets import (
    QWidget, QHBoxLayout, QLineEdit, QPushButton, QFileDialog, QLabel,
)
from PySide6.QtCore import Signal


class _BasePicker(QWidget):
    path_changed = Signal(str)

    def __init__(self, label, placeholder, parent=None):
        super().__init__(parent)
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(6)
        lbl = QLabel(label)
        lbl.setFixedWidth(100)
        layout.addWidget(lbl)
        self.line = QLineEdit()
        self.line.setPlaceholderText(placeholder)
        self.line.textChanged.connect(self.path_changed.emit)
        layout.addWidget(self.line, stretch=1)
        btn = QPushButton("  Browse  ")
        btn.setProperty("secondary", True)
        btn.setFixedHeight(28)
        btn.clicked.connect(self._browse)
        layout.addWidget(btn)

    def _browse(self): ...
    def path(self) -> str: return self.line.text().strip()
    def set_path(self, p: str): self.line.setText(p)


class DirPicker(_BasePicker):
    def __init__(self, label="Directory:", parent=None,
                 dialog_title="Select Directory"):
        super().__init__(label, "Select a directory...", parent)
        self._title = dialog_title

    def _browse(self):
        d = QFileDialog.getExistingDirectory(self, self._title, self.line.text())
        if d: self.line.setText(d)


class FilePicker(_BasePicker):
    def __init__(self, label="File:", parent=None,
                 dialog_title="Select File",
                 file_filter="CSV Files (*.csv);;All Files (*)"):
        super().__init__(label, "Select a file...", parent)
        self._title = dialog_title
        self._filter = file_filter

    def _browse(self):
        f, _ = QFileDialog.getOpenFileName(self, self._title, self.line.text(), self._filter)
        if f: self.line.setText(f)
