"""Real-time log panel with color-coded messages."""
from __future__ import annotations

import logging
import html
from datetime import datetime

from PySide6.QtWidgets import QTextEdit, QVBoxLayout, QWidget, QHBoxLayout, QPushButton, QLabel
from PySide6.QtCore import Signal, QObject
from PySide6.QtGui import QColor, QTextCharFormat, QTextCursor

from gui.theme import _DARK, _LIGHT


class _Sig(QObject):
    message = Signal(str, str)


class LogPanel(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self._colors = _DARK
        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.setSpacing(2)

        header = QHBoxLayout()
        lbl = QLabel("Log")
        lbl.setStyleSheet("font-weight: bold; font-size: 9pt;")
        header.addWidget(lbl)
        header.addStretch()
        btn = QPushButton("Clear Log")
        btn.setProperty("secondary", True)
        btn.setFixedHeight(24)
        btn.clicked.connect(lambda: self.text.clear())
        header.addWidget(btn)
        layout.addLayout(header)

        self.text = QTextEdit()
        self.text.setReadOnly(True)
        self.text.setLineWrapMode(QTextEdit.NoWrap)
        layout.addWidget(self.text)

        self._sig = _Sig()
        self._sig.message.connect(self._append)
        self._handler = _QtLogHandler(self)
        self._handler.setFormatter(logging.Formatter("%(message)s"))

    def get_handler(self): return self._handler

    def log(self, msg, level="info"):
        self._sig.message.emit(msg, level)

    def refresh_theme(self, dark):
        self._colors = _DARK if dark else _LIGHT
        roles = ("text_dim", "success", "warning", "error", "accent")
        replacements = {palette[role]: self._colors[role] for palette in (_DARK, _LIGHT) for role in roles}
        vertical, horizontal = self.text.verticalScrollBar(), self.text.horizontalScrollBar()
        position, side = vertical.value(), horizontal.value()
        at_end = position == vertical.maximum()
        # Change native text formats only; paths and message text stay exact.
        changes = []
        block = self.text.document().begin()
        while block.isValid():
            fragments = block.begin()
            while not fragments.atEnd():
                fragment = fragments.fragment()
                if fragment.isValid():
                    color = fragment.charFormat().foreground().color().name()
                    if color in replacements:
                        changes.append((fragment.position(), fragment.length(), replacements[color]))
                fragments += 1
            block = block.next()
        cursor = QTextCursor(self.text.document())
        cursor.beginEditBlock()
        for start, length, color in changes:
            cursor.setPosition(start)
            cursor.setPosition(start + length, QTextCursor.KeepAnchor)
            style = QTextCharFormat()
            style.setForeground(QColor(color))
            cursor.mergeCharFormat(style)
        cursor.endEditBlock()
        vertical.setValue(vertical.maximum() if at_end else position)
        horizontal.setValue(side)

    def _append(self, msg, level):
        role = {"info": "text_dim", "success": "success", "warning": "warning",
                "error": "error", "accent": "accent"}.get(level, "text_dim")
        c = self._colors[role]
        ts = datetime.now().strftime("%H:%M:%S")
        self.text.append(f'<span style="color:{self._colors["text_dim"]}">{ts}</span> '
                         f'<span style="color:{c}">{html.escape(str(msg))}</span>')
        self.text.verticalScrollBar().setValue(self.text.verticalScrollBar().maximum())


class _QtLogHandler(logging.Handler):
    def __init__(self, panel):
        super().__init__()
        self._panel = panel

    def emit(self, record):
        levels = {logging.WARNING: "warning", logging.ERROR: "error", logging.CRITICAL: "error"}
        self._panel.log(self.format(record), levels.get(record.levelno, "info"))
