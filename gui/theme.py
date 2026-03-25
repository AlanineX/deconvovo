"""Theme for DeconVoVo GUI — modern dark/light with depth layering.

Design principles from PyDracula, QFluentWidgets, qt-material:
- Solid backgrounds with 4-layer depth hierarchy (no OS transparency)
- 8px grid spacing system
- Single accent color used sparingly
- Card-based grouping with rounded corners
- Typography hierarchy: title > section > body > caption
"""
from __future__ import annotations

import atexit
import os
import sys
import tempfile

# -- Accent + status colors (same in dark and light) --
ACCENT = "#4a9eff"
ACCENT_HOVER = "#6ab4ff"
ACCENT_PRESSED = "#3580d8"
ACCENT_DIM = "rgba(74,158,255,40)"
SUCCESS = "#3ddc84"
WARNING = "#ffb74d"
ERROR = "#ef5350"

# -- Fonts --
FONT_FAMILY = "Segoe UI"
FONT_SIZE = 10
FONT_SIZE_SMALL = 9
FONT_SIZE_TITLE = 14
FONT_SIZE_SECTION = 11
FONT_SIZE_CAPTION = 8

# -- 4-layer dark palette --
_DARK = dict(
    L0="#16161e", L1="#1e1e2e", L2="#272740", L3="#323252",
    border="#3a3a5a", border_hi="#4a4a6a",
    text="#e0e0f0", text_dim="#7a7a9a", text_hi="#ffffff",
    plot_bg="#16161e", plot_face="#1e1e2e", plot_hex="#1e1e2e",
    grid="#2e2e4e", log_bg="#12121a",
)
_LIGHT = dict(
    L0="#f0f0f5", L1="#ffffff", L2="#e8e8f0", L3="#d8d8e4",
    border="#c4c4d4", border_hi="#a0a0b8",
    text="#1a1a2e", text_dim="#6a6a80", text_hi="#000000",
    plot_bg="#f8f8fc", plot_face="#ffffff", plot_hex="#ffffff",
    grid="#d0d0e0", log_bg="#f4f4f8",
)

TEXT_DIM = _DARK["text_dim"]
TEXT = _DARK["text"]

# -- Checkmark PNG --
_checkmark_path: str | None = None

def _ensure_checkmark_png() -> str:
    """Create a 16x16 white checkmark PNG. Cached for session lifetime."""
    global _checkmark_path
    if _checkmark_path and os.path.isfile(_checkmark_path):
        return _checkmark_path
    from PySide6.QtCore import Qt
    from PySide6.QtGui import QColor, QPainter, QPen, QPixmap

    pm = QPixmap(16, 16)
    pm.fill(QColor(0, 0, 0, 0))
    painter = QPainter(pm)
    painter.setRenderHint(QPainter.RenderHint.Antialiasing)
    pen = QPen(QColor(255, 255, 255), 2.5)
    pen.setCapStyle(Qt.PenCapStyle.RoundCap)
    pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
    painter.setPen(pen)
    painter.drawLine(3, 8, 6, 11)
    painter.drawLine(6, 11, 13, 4)
    painter.end()

    # Write next to the exe (frozen) or in temp (dev)
    if getattr(sys, 'frozen', False):
        d = os.path.dirname(sys.executable)
    else:
        d = tempfile.gettempdir()
    path = os.path.join(d, "dvv_checkmark.png")
    pm.save(path, "PNG")
    _checkmark_path = path
    atexit.register(lambda p=path: os.remove(p) if os.path.isfile(p) else None)
    return path


def mpl_rc(dark=True):
    p = _DARK if dark else _LIGHT
    return {
        "figure.facecolor": p["plot_bg"], "axes.facecolor": p["plot_face"],
        "axes.edgecolor": p["border"], "axes.labelcolor": p["text"],
        "axes.titlecolor": p["text_hi"], "xtick.color": p["text_dim"],
        "ytick.color": p["text_dim"], "text.color": p["text"],
        "grid.color": p["grid"], "grid.alpha": 0.3,
        "font.family": "sans-serif",
        "font.sans-serif": ["Segoe UI", "Liberation Sans", "Arial"],
        "font.size": 10, "legend.facecolor": p["plot_hex"],
        "legend.edgecolor": p["border"], "figure.dpi": 100,
    }


def stylesheet(dark=True):
    p = _DARK if dark else _LIGHT
    A, AH, AP, AD = ACCENT, ACCENT_HOVER, ACCENT_PRESSED, ACCENT_DIM
    chk = _ensure_checkmark_png().replace("\\", "/")
    return f"""

    /* ===== BASE ===== */
    QWidget {{
        background-color: {p['L0']};
        color: {p['text']};
        font-family: "{FONT_FAMILY}", "Liberation Sans", "Arial", sans-serif;
        font-size: {FONT_SIZE}pt;
        border: none;
    }}
    QMainWindow {{ background-color: {p['L0']}; }}

    /* ===== SIDEBAR ===== */
    QListWidget#nav {{
        background-color: {p['L1']};
        border: none;
        border-right: 1px solid {p['border']};
        outline: none;
        padding: 4px 0;
        font-size: {FONT_SIZE_SECTION}pt;
    }}
    QListWidget#nav::item {{
        padding: 14px 20px;
        border: none;
        border-left: 3px solid transparent;
        color: {p['text_dim']};
        margin: 2px 4px;
        border-radius: 6px;
    }}
    QListWidget#nav::item:selected {{
        background-color: {AD};
        color: {p['text_hi']};
        border-left: 3px solid {A};
    }}
    QListWidget#nav::item:hover:!selected {{
        background-color: {p['L2']};
        color: {p['text']};
    }}

    /* ===== BUTTONS ===== */
    QPushButton {{
        background-color: {A};
        color: {p['text_hi']};
        border: none;
        border-radius: 6px;
        padding: 8px 24px;
        font-weight: 600;
        font-size: {FONT_SIZE}pt;
        min-height: 20px;
    }}
    QPushButton:hover {{ background-color: {AH}; }}
    QPushButton:pressed {{ background-color: {AP}; }}
    QPushButton:disabled {{ background-color: {p['L2']}; color: {p['text_dim']}; }}
    QPushButton[secondary="true"] {{
        background-color: {p['L2']};
        color: {p['text']};
        border: 1px solid {p['border']};
        font-weight: normal;
    }}
    QPushButton[secondary="true"]:hover {{
        background-color: {p['L3']};
        border-color: {p['border_hi']};
    }}

    /* ===== INPUTS ===== */
    QLineEdit, QSpinBox, QDoubleSpinBox, QComboBox {{
        background-color: {p['L2']};
        color: {p['text']};
        border: 1px solid {p['border']};
        border-radius: 6px;
        padding: 6px 10px;
        font-size: {FONT_SIZE}pt;
        selection-background-color: {AD};
    }}
    QLineEdit:focus, QSpinBox:focus, QDoubleSpinBox:focus, QComboBox:focus {{
        border: 1px solid {A};
    }}
    QComboBox::drop-down {{ border: none; padding-right: 10px; }}
    QComboBox QAbstractItemView {{
        background-color: {p['L2']};
        color: {p['text']};
        selection-background-color: {A};
        border: 1px solid {p['border']};
        border-radius: 4px;
    }}

    /* ===== LABELS ===== */
    QLabel {{ background: transparent; color: {p['text']}; }}
    QLabel[title="true"] {{
        font-size: {FONT_SIZE_TITLE}pt;
        font-weight: 700;
        color: {p['text_hi']};
        padding-bottom: 2px;
    }}
    QLabel[subtitle="true"] {{
        font-size: {FONT_SIZE_SMALL}pt;
        color: {p['text_dim']};
        padding-bottom: 8px;
    }}

    /* ===== CARDS (GroupBox) ===== */
    QGroupBox {{
        background-color: {p['L1']};
        border: 1px solid {p['border']};
        border-radius: 8px;
        margin-top: 16px;
        padding: 20px 16px 12px 16px;
        font-weight: 600;
        font-size: {FONT_SIZE_SMALL}pt;
        color: {p['text_dim']};
    }}
    QGroupBox::title {{
        subcontrol-origin: margin;
        left: 16px;
        padding: 0 8px;
        color: {p['text_dim']};
    }}

    /* ===== TABLE ===== */
    QTableWidget {{
        background-color: {p['L1']};
        alternate-background-color: {p['L2']};
        gridline-color: {p['border']};
        border: 1px solid {p['border']};
        border-radius: 6px;
        font-size: {FONT_SIZE_SMALL}pt;
    }}
    QTableWidget::item {{ padding: 4px 8px; }}
    QTableWidget::item:selected {{ background-color: {A}; color: {p['text_hi']}; }}
    QHeaderView::section {{
        background-color: {p['L3']};
        color: {p['text_hi']};
        border: none;
        border-right: 1px solid {p['border']};
        border-bottom: 1px solid {p['border']};
        padding: 6px 10px;
        font-weight: 600;
        font-size: {FONT_SIZE_SMALL}pt;
    }}

    /* ===== SCROLLBAR ===== */
    QScrollBar:vertical {{
        background: transparent;
        width: 8px;
        border: none;
        margin: 4px 2px;
    }}
    QScrollBar::handle:vertical {{
        background-color: {p['border']};
        border-radius: 4px;
        min-height: 32px;
    }}
    QScrollBar::handle:vertical:hover {{ background-color: {p['text_dim']}; }}
    QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {{ height: 0; }}
    QScrollBar:horizontal {{
        background: transparent;
        height: 8px;
        border: none;
        margin: 2px 4px;
    }}
    QScrollBar::handle:horizontal {{
        background-color: {p['border']};
        border-radius: 4px;
        min-width: 32px;
    }}

    /* ===== PROGRESS BAR ===== */
    QProgressBar {{
        background-color: {p['L2']};
        border: none;
        border-radius: 4px;
        text-align: center;
        color: {p['text_hi']};
        font-size: {FONT_SIZE_SMALL}pt;
        min-height: 6px;
        max-height: 6px;
    }}
    QProgressBar::chunk {{
        background-color: {A};
        border-radius: 3px;
    }}

    /* ===== TABS ===== */
    QTabWidget::pane {{
        border: 1px solid {p['border']};
        background-color: {p['L1']};
        border-radius: 6px;
        top: -1px;
    }}
    QTabBar::tab {{
        background-color: transparent;
        color: {p['text_dim']};
        border: none;
        border-bottom: 2px solid transparent;
        padding: 8px 20px;
        margin-right: 4px;
        font-size: {FONT_SIZE_SMALL}pt;
    }}
    QTabBar::tab:selected {{
        color: {p['text_hi']};
        border-bottom: 2px solid {A};
    }}
    QTabBar::tab:hover:!selected {{
        color: {p['text']};
        border-bottom: 2px solid {p['border']};
    }}

    /* ===== SPLITTER ===== */
    QSplitter::handle {{
        background-color: {p['border']};
    }}
    QSplitter::handle:hover {{
        background-color: {A};
    }}
    QSplitter::handle:horizontal {{ width: 3px; margin: 0 2px; }}
    QSplitter::handle:vertical {{ height: 3px; margin: 2px 0; }}

    /* ===== LOG PANEL ===== */
    QTextEdit[readOnly="true"] {{
        background-color: {p['log_bg']};
        color: {p['text_dim']};
        border: 1px solid {p['border']};
        border-radius: 6px;
        font-family: "Cascadia Code", "Consolas", "Liberation Mono", monospace;
        font-size: {FONT_SIZE_SMALL}pt;
        padding: 8px;
    }}

    /* ===== TOOLTIP ===== */
    QToolTip {{
        background-color: {p['L3']};
        color: {p['text']};
        border: 1px solid {p['border']};
        border-radius: 4px;
        padding: 6px 10px;
        font-size: {FONT_SIZE_SMALL}pt;
    }}

    /* ===== CHECKBOX ===== */
    QCheckBox {{ background: transparent; spacing: 8px; }}
    QCheckBox::indicator {{
        width: 18px; height: 18px;
        border: 2px solid {p['border']};
        border-radius: 4px;
        background-color: {p['L2']};
    }}
    QCheckBox::indicator:hover {{ border-color: {A}; }}
    QCheckBox::indicator:checked {{
        border-color: {A};
        background-color: {A};
        image: url("{chk}");
    }}
    """
