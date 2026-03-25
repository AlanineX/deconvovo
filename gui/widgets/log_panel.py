"""Real-time log panel with color-coded messages."""
from __future__ import annotations

import logging
from datetime import datetime

from PySide6.QtWidgets import QTextEdit, QVBoxLayout, QWidget, QHBoxLayout, QPushButton, QLabel
from PySide6.QtCore import Qt, Signal, QObject

from gui.theme import SUCCESS, WARNING, ERROR, TEXT_DIM, TEXT, ACCENT


class _LogSignal(QObject):
    """Thread-safe signal for log messages from worker threads."""
    message = Signal(str, str)  # (html, level)


class LogPanel(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.setSpacing(2)

        # Header
        header = QHBoxLayout()
        lbl = QLabel("Log")
        lbl.setStyleSheet(f"font-weight: bold; color: {TEXT}; font-size: 9pt;")
        header.addWidget(lbl)
        header.addStretch()
        btn_clear = QPushButton("Clear Log")
        btn_clear.setProperty("secondary", True)
        btn_clear.setFixedHeight(24)
        btn_clear.clicked.connect(self.clear)
        header.addWidget(btn_clear)
        layout.addLayout(header)

        self.text = QTextEdit()
        self.text.setReadOnly(True)
        self.text.setLineWrapMode(QTextEdit.NoWrap)
        layout.addWidget(self.text)

        # Thread-safe signal
        self._sig = _LogSignal()
        self._sig.message.connect(self._append)

        # Python logging handler
        self._handler = _QtLogHandler(self)
        self._handler.setFormatter(logging.Formatter("%(message)s"))

    def get_handler(self) -> logging.Handler:
        return self._handler

    def log(self, msg: str, level: str = "info"):
        """Thread-safe log. Can be called from any thread."""
        self._sig.message.emit(msg, level)

    def _append(self, msg: str, level: str):
        color = {
            "info": TEXT_DIM,
            "success": SUCCESS,
            "warning": WARNING,
            "error": ERROR,
            "accent": ACCENT,
        }.get(level, TEXT_DIM)
        ts = datetime.now().strftime("%H:%M:%S")
        self.text.append(f'<span style="color:{TEXT_DIM}">{ts}</span> '
                         f'<span style="color:{color}">{msg}</span>')
        self.text.verticalScrollBar().setValue(self.text.verticalScrollBar().maximum())

    def clear(self):
        self.text.clear()


class _QtLogHandler(logging.Handler):
    """Bridge Python logging → LogPanel."""
    def __init__(self, panel: LogPanel):
        super().__init__()
        self._panel = panel

    def emit(self, record):
        level = {
            logging.DEBUG: "info",
            logging.INFO: "info",
            logging.WARNING: "warning",
            logging.ERROR: "error",
            logging.CRITICAL: "error",
        }.get(record.levelno, "info")
        self._panel.log(self.format(record), level)
