"""Directory and file picker widgets."""
from __future__ import annotations

from pathlib import Path

from PySide6.QtWidgets import (
    QWidget, QHBoxLayout, QLineEdit, QPushButton, QFileDialog, QLabel,
)
from PySide6.QtCore import Signal

from gui.theme import TEXT_DIM


class DirPicker(QWidget):
    """Directory browser with label and text field."""

    path_changed = Signal(str)

    def __init__(self, label: str = "Directory:", parent=None,
                 dialog_title: str = "Select Directory"):
        super().__init__(parent)
        self._dialog_title = dialog_title
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(6)

        lbl = QLabel(label)
        lbl.setFixedWidth(100)
        layout.addWidget(lbl)

        self.line = QLineEdit()
        self.line.setPlaceholderText("Select a directory...")
        self.line.textChanged.connect(lambda t: self.path_changed.emit(t))
        layout.addWidget(self.line, stretch=1)

        btn = QPushButton("Browse")
        btn.setProperty("secondary", True)
        btn.setFixedHeight(28)
        btn.setFixedWidth(70)
        btn.clicked.connect(self._browse)
        layout.addWidget(btn)

    def _browse(self):
        d = QFileDialog.getExistingDirectory(self, self._dialog_title, self.line.text())
        if d:
            self.line.setText(d)

    def path(self) -> str:
        return self.line.text().strip()

    def set_path(self, p: str):
        self.line.setText(p)


class FilePicker(QWidget):
    """File browser with label and text field."""

    path_changed = Signal(str)

    def __init__(self, label: str = "File:", parent=None,
                 dialog_title: str = "Select File",
                 file_filter: str = "CSV Files (*.csv);;All Files (*)"):
        super().__init__(parent)
        self._dialog_title = dialog_title
        self._filter = file_filter
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(6)

        lbl = QLabel(label)
        lbl.setFixedWidth(100)
        layout.addWidget(lbl)

        self.line = QLineEdit()
        self.line.setPlaceholderText("Select a file...")
        self.line.textChanged.connect(lambda t: self.path_changed.emit(t))
        layout.addWidget(self.line, stretch=1)

        btn = QPushButton("Browse")
        btn.setProperty("secondary", True)
        btn.setFixedHeight(28)
        btn.setFixedWidth(70)
        btn.clicked.connect(self._browse)
        layout.addWidget(btn)

    def _browse(self):
        f, _ = QFileDialog.getOpenFileName(
            self, self._dialog_title, self.line.text(), self._filter)
        if f:
            self.line.setText(f)

    def path(self) -> str:
        return self.line.text().strip()

    def set_path(self, p: str):
        self.line.setText(p)
