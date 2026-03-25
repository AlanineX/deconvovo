"""Theme for DeconVoVo GUI — dark and light modes.

Color palette inspired by scientific instrument software.
"""
from __future__ import annotations

# -- Shared colors (work on both backgrounds) --
ACCENT = "#5b8af5"
ACCENT_HOVER = "#7aa4ff"
ACCENT_PRESSED = "#4070d8"
SUCCESS = "#2ea870"
WARNING = "#c89a20"
ERROR = "#d04840"

# -- Fonts --
FONT_FAMILY = "Segoe UI"
FONT_SIZE = 10
FONT_SIZE_SMALL = 9
FONT_SIZE_TITLE = 13
FONT_SIZE_HEADER = 11

# -- Dark palette --
# Qt uses rgba for translucency, matplotlib needs hex
_DARK = {
    "bg":       "rgba(30, 30, 46, 230)",
    "bg_mid":   "rgba(37, 37, 56, 240)",
    "bg_mid_hex": "#252538",   # hex fallback for mpl
    "bg_input": "#2e2e44",
    "bg_head":  "#363654",
    "border":   "#3e3e5e",
    "text":     "#d4d4e8",
    "text_dim": "#8888aa",
    "text_hi":  "#f0f0ff",
    "plot_bg":  "#1a1a2e",
    "plot_face":"#22223a",
    "grid":     "#3a3a5a",
    "log_bg":   "#1a1a2e",
}

# -- Light palette --
_LIGHT = {
    "bg":       "rgba(244, 244, 248, 230)",
    "bg_mid":   "rgba(255, 255, 255, 240)",
    "bg_mid_hex": "#ffffff",
    "bg_input": "#e8e8f0",
    "bg_head":  "#dcdce8",
    "border":   "#c0c0d0",
    "text":     "#2a2a3a",
    "text_dim": "#6a6a80",
    "text_hi":  "#1a1a2a",
    "plot_bg":  "#f8f8fc",
    "plot_face":"#ffffff",
    "grid":     "#d0d0e0",
    "log_bg":   "#f0f0f6",
}

# -- Current mode (mutable) --
_current_mode = "dark"

def _pal():
    return _DARK if _current_mode == "dark" else _LIGHT


# Module-level accessors for panels that import these
def _get(key):
    return _pal()[key]

# Convenience — panels import these directly
BG_DARK = property(lambda self: _pal()["bg"])
BG_MID = property(lambda self: _pal()["bg_mid"])
TEXT = "#d4d4e8"       # default; panels should use stylesheet, not inline
TEXT_DIM = "#8888aa"
TEXT_BRIGHT = "#f0f0ff"
PLOT_BG = "#1a1a2e"
PLOT_FACE = "#22223a"


def mpl_rc(dark: bool = True) -> dict:
    p = _DARK if dark else _LIGHT
    return {
        "figure.facecolor": p["plot_bg"],
        "axes.facecolor": p["plot_face"],
        "axes.edgecolor": p["border"],
        "axes.labelcolor": p["text"],
        "axes.titlecolor": p["text_hi"],
        "xtick.color": p["text_dim"],
        "ytick.color": p["text_dim"],
        "text.color": p["text"],
        "grid.color": p["grid"],
        "grid.alpha": 0.4,
        "font.family": "sans-serif",
        "font.sans-serif": ["Segoe UI", "Liberation Sans", "Arial"],
        "font.size": 10,
        "legend.facecolor": p["bg_mid_hex"],
        "legend.edgecolor": p["border"],
        "figure.dpi": 100,
    }

# Keep MPL_RC for backward compat
MPL_RC = mpl_rc(dark=True)


def stylesheet(dark: bool = True) -> str:
    """Return the full Qt stylesheet for the given mode."""
    p = _DARK if dark else _LIGHT
    return f"""
    QWidget {{
        background-color: {p['bg']};
        color: {p['text']};
        font-family: "{FONT_FAMILY}", "Liberation Sans", "Arial", sans-serif;
        font-size: {FONT_SIZE}pt;
    }}
    QMainWindow {{ background-color: {p['bg']}; }}

    QListWidget {{
        background-color: {p['bg_mid']};
        border: none; border-right: 1px solid {p['border']};
        outline: none; padding: 4px 0px; font-size: {FONT_SIZE_HEADER}pt;
    }}
    QListWidget::item {{ padding: 12px 16px; border: none; color: {p['text_dim']}; }}
    QListWidget::item:selected {{
        background-color: {ACCENT}; color: {p['text_hi']};
        border-left: 3px solid {ACCENT_HOVER};
    }}
    QListWidget::item:hover:!selected {{ background-color: {p['bg_head']}; color: {p['text']}; }}

    QPushButton {{
        background-color: {ACCENT}; color: {p['text_hi']};
        border: none; border-radius: 4px; padding: 8px 20px;
        font-weight: bold; font-size: {FONT_SIZE}pt; min-height: 18px;
    }}
    QPushButton:hover {{ background-color: {ACCENT_HOVER}; }}
    QPushButton:pressed {{ background-color: {ACCENT_PRESSED}; }}
    QPushButton:disabled {{ background-color: {p['bg_input']}; color: {p['text_dim']}; }}
    QPushButton[secondary="true"] {{
        background-color: {p['bg_input']}; color: {p['text']};
        border: 1px solid {p['border']};
    }}
    QPushButton[secondary="true"]:hover {{ background-color: {p['bg_head']}; }}

    QLineEdit, QSpinBox, QDoubleSpinBox, QComboBox {{
        background-color: {p['bg_input']}; color: {p['text']};
        border: 1px solid {p['border']}; border-radius: 3px;
        padding: 5px 8px; font-size: {FONT_SIZE}pt;
    }}
    QLineEdit:focus, QSpinBox:focus, QDoubleSpinBox:focus, QComboBox:focus {{
        border: 1px solid {ACCENT};
    }}
    QComboBox::drop-down {{ border: none; padding-right: 8px; }}
    QComboBox QAbstractItemView {{
        background-color: {p['bg_input']}; color: {p['text']};
        selection-background-color: {ACCENT}; border: 1px solid {p['border']};
    }}

    QLabel {{ background-color: transparent; color: {p['text']}; }}
    QLabel[title="true"] {{ font-size: {FONT_SIZE_TITLE}pt; font-weight: bold; color: {p['text_hi']}; }}
    QLabel[subtitle="true"] {{ font-size: {FONT_SIZE_SMALL}pt; color: {p['text_dim']}; }}

    QGroupBox {{
        background-color: {p['bg_mid']}; border: 1px solid {p['border']};
        border-radius: 5px; margin-top: 12px; padding-top: 18px;
        font-weight: bold; color: {p['text_hi']};
    }}
    QGroupBox::title {{ subcontrol-origin: margin; left: 12px; padding: 0 6px; }}

    QTableWidget {{
        background-color: {p['bg_mid']}; alternate-background-color: {p['bg_input']};
        gridline-color: {p['border']}; border: 1px solid {p['border']};
        border-radius: 3px; font-size: {FONT_SIZE_SMALL}pt;
    }}
    QTableWidget::item {{ padding: 3px 6px; }}
    QTableWidget::item:selected {{ background-color: {ACCENT}; color: {p['text_hi']}; }}
    QHeaderView::section {{
        background-color: {p['bg_head']}; color: {p['text_hi']};
        border: none; border-right: 1px solid {p['border']};
        border-bottom: 1px solid {p['border']};
        padding: 5px 8px; font-weight: bold; font-size: {FONT_SIZE_SMALL}pt;
    }}

    QScrollBar:vertical {{
        background-color: {p['bg']}; width: 10px; border: none;
    }}
    QScrollBar::handle:vertical {{
        background-color: {p['border']}; border-radius: 5px; min-height: 30px;
    }}
    QScrollBar::handle:vertical:hover {{ background-color: {p['text_dim']}; }}
    QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {{ height: 0px; }}
    QScrollBar:horizontal {{
        background-color: {p['bg']}; height: 10px; border: none;
    }}
    QScrollBar::handle:horizontal {{
        background-color: {p['border']}; border-radius: 5px; min-width: 30px;
    }}

    QProgressBar {{
        background-color: {p['bg_input']}; border: 1px solid {p['border']};
        border-radius: 4px; text-align: center; color: {p['text_hi']};
        font-size: {FONT_SIZE_SMALL}pt; min-height: 20px;
    }}
    QProgressBar::chunk {{ background-color: {ACCENT}; border-radius: 3px; }}

    QTabWidget::pane {{
        border: 1px solid {p['border']}; background-color: {p['bg_mid']};
        border-radius: 3px;
    }}
    QTabBar::tab {{
        background-color: {p['bg_input']}; color: {p['text_dim']};
        border: 1px solid {p['border']}; border-bottom: none;
        padding: 7px 16px; margin-right: 2px;
        border-top-left-radius: 4px; border-top-right-radius: 4px;
    }}
    QTabBar::tab:selected {{
        background-color: {p['bg_mid']}; color: {p['text_hi']};
        border-bottom: 2px solid {ACCENT};
    }}
    QTabBar::tab:hover:!selected {{ background-color: {p['bg_head']}; color: {p['text']}; }}

    QSplitter::handle {{ background-color: {p['border']}; }}
    QSplitter::handle:horizontal {{ width: 1px; }}
    QSplitter::handle:vertical {{ height: 1px; }}

    QTextEdit[readOnly="true"] {{
        background-color: {p['log_bg']}; color: {p['text_dim']};
        border: 1px solid {p['border']}; border-radius: 3px;
        font-family: "Consolas", "Liberation Mono", monospace;
        font-size: {FONT_SIZE_SMALL}pt;
    }}

    QToolTip {{
        background-color: {p['bg_head']}; color: {p['text']};
        border: 1px solid {p['border']}; padding: 4px;
    }}

    QCheckBox {{ background-color: transparent; spacing: 6px; }}
    QCheckBox::indicator {{
        width: 18px; height: 18px; border: 2px solid {p['border']};
        border-radius: 4px; background-color: {p['bg_input']};
    }}
    QCheckBox::indicator:hover {{
        border-color: {ACCENT};
    }}
    QCheckBox::indicator:checked {{
        border-color: {ACCENT}; background-color: {ACCENT};
        image: url("data:image/svg+xml;utf8,<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 16 16'><path d='M3 8l3 3 7-7' stroke='white' stroke-width='2.5' fill='none' stroke-linecap='round' stroke-linejoin='round'/></svg>");
    }}
    """
