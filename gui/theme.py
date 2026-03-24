"""Dark scientific theme for DeconVoVo GUI.

Color palette inspired by scientific instrument software (MassLynx, Xcalibur).
"""
from __future__ import annotations

# -- Color palette --
BG_DARK = "#1e1e2e"       # main background
BG_MID = "#252538"         # panel background
BG_LIGHT = "#2e2e44"       # input fields, table cells
BG_HEADER = "#363654"      # table headers, sidebar hover
BORDER = "#3e3e5e"         # borders, separators
TEXT = "#d4d4e8"           # primary text
TEXT_DIM = "#8888aa"       # secondary / placeholder
TEXT_BRIGHT = "#f0f0ff"    # titles, emphasis
ACCENT = "#5b8af5"         # primary accent (buttons, selected)
ACCENT_HOVER = "#7aa4ff"   # hover state
ACCENT_PRESSED = "#4070d8" # pressed state
SUCCESS = "#4ec990"        # green — completed, good R²
WARNING = "#e8b84a"        # yellow — warnings
ERROR = "#e8584a"          # red — errors
PLOT_BG = "#1a1a2e"        # matplotlib figure background
PLOT_FACE = "#22223a"      # matplotlib axes background
GRID = "#3a3a5a"           # matplotlib grid

# -- Fonts --
FONT_FAMILY = "Segoe UI"   # Windows default; fallback in stylesheet
FONT_SIZE = 10
FONT_SIZE_SMALL = 9
FONT_SIZE_TITLE = 13
FONT_SIZE_HEADER = 11

# -- Matplotlib rc overrides --
MPL_RC = {
    "figure.facecolor": PLOT_BG,
    "axes.facecolor": PLOT_FACE,
    "axes.edgecolor": BORDER,
    "axes.labelcolor": TEXT,
    "axes.titlecolor": TEXT_BRIGHT,
    "xtick.color": TEXT_DIM,
    "ytick.color": TEXT_DIM,
    "text.color": TEXT,
    "grid.color": GRID,
    "grid.alpha": 0.4,
    "font.family": "sans-serif",
    "font.sans-serif": ["Segoe UI", "Liberation Sans", "Arial"],
    "font.size": 10,
    "legend.facecolor": BG_MID,
    "legend.edgecolor": BORDER,
    "figure.dpi": 100,
}


def stylesheet() -> str:
    """Return the full Qt stylesheet."""
    return f"""
    /* ---- Global ---- */
    QWidget {{
        background-color: {BG_DARK};
        color: {TEXT};
        font-family: "{FONT_FAMILY}", "Liberation Sans", "Arial", sans-serif;
        font-size: {FONT_SIZE}pt;
    }}

    /* ---- Main window ---- */
    QMainWindow {{
        background-color: {BG_DARK};
    }}

    /* ---- Sidebar ---- */
    QListWidget {{
        background-color: {BG_MID};
        border: none;
        border-right: 1px solid {BORDER};
        outline: none;
        padding: 4px 0px;
        font-size: {FONT_SIZE_HEADER}pt;
    }}
    QListWidget::item {{
        padding: 12px 16px;
        border: none;
        color: {TEXT_DIM};
    }}
    QListWidget::item:selected {{
        background-color: {ACCENT};
        color: {TEXT_BRIGHT};
        border-left: 3px solid {ACCENT_HOVER};
    }}
    QListWidget::item:hover:!selected {{
        background-color: {BG_HEADER};
        color: {TEXT};
    }}

    /* ---- Buttons ---- */
    QPushButton {{
        background-color: {ACCENT};
        color: {TEXT_BRIGHT};
        border: none;
        border-radius: 4px;
        padding: 8px 20px;
        font-weight: bold;
        font-size: {FONT_SIZE}pt;
        min-height: 18px;
    }}
    QPushButton:hover {{
        background-color: {ACCENT_HOVER};
    }}
    QPushButton:pressed {{
        background-color: {ACCENT_PRESSED};
    }}
    QPushButton:disabled {{
        background-color: {BG_LIGHT};
        color: {TEXT_DIM};
    }}
    QPushButton[secondary="true"] {{
        background-color: {BG_LIGHT};
        color: {TEXT};
        border: 1px solid {BORDER};
    }}
    QPushButton[secondary="true"]:hover {{
        background-color: {BG_HEADER};
    }}

    /* ---- Input fields ---- */
    QLineEdit, QSpinBox, QDoubleSpinBox, QComboBox {{
        background-color: {BG_LIGHT};
        color: {TEXT};
        border: 1px solid {BORDER};
        border-radius: 3px;
        padding: 5px 8px;
        font-size: {FONT_SIZE}pt;
    }}
    QLineEdit:focus, QSpinBox:focus, QDoubleSpinBox:focus, QComboBox:focus {{
        border: 1px solid {ACCENT};
    }}
    QComboBox::drop-down {{
        border: none;
        padding-right: 8px;
    }}
    QComboBox QAbstractItemView {{
        background-color: {BG_LIGHT};
        color: {TEXT};
        selection-background-color: {ACCENT};
        border: 1px solid {BORDER};
    }}

    /* ---- Labels ---- */
    QLabel {{
        background-color: transparent;
        color: {TEXT};
    }}
    QLabel[title="true"] {{
        font-size: {FONT_SIZE_TITLE}pt;
        font-weight: bold;
        color: {TEXT_BRIGHT};
    }}
    QLabel[subtitle="true"] {{
        font-size: {FONT_SIZE_SMALL}pt;
        color: {TEXT_DIM};
    }}

    /* ---- Group box ---- */
    QGroupBox {{
        background-color: {BG_MID};
        border: 1px solid {BORDER};
        border-radius: 5px;
        margin-top: 12px;
        padding-top: 18px;
        font-weight: bold;
        color: {TEXT_BRIGHT};
    }}
    QGroupBox::title {{
        subcontrol-origin: margin;
        left: 12px;
        padding: 0 6px;
    }}

    /* ---- Table ---- */
    QTableWidget {{
        background-color: {BG_MID};
        alternate-background-color: {BG_LIGHT};
        gridline-color: {BORDER};
        border: 1px solid {BORDER};
        border-radius: 3px;
        font-size: {FONT_SIZE_SMALL}pt;
    }}
    QTableWidget::item {{
        padding: 3px 6px;
    }}
    QTableWidget::item:selected {{
        background-color: {ACCENT};
        color: {TEXT_BRIGHT};
    }}
    QHeaderView::section {{
        background-color: {BG_HEADER};
        color: {TEXT_BRIGHT};
        border: none;
        border-right: 1px solid {BORDER};
        border-bottom: 1px solid {BORDER};
        padding: 5px 8px;
        font-weight: bold;
        font-size: {FONT_SIZE_SMALL}pt;
    }}

    /* ---- Scrollbar ---- */
    QScrollBar:vertical {{
        background-color: {BG_DARK};
        width: 10px;
        border: none;
    }}
    QScrollBar::handle:vertical {{
        background-color: {BORDER};
        border-radius: 5px;
        min-height: 30px;
    }}
    QScrollBar::handle:vertical:hover {{
        background-color: {TEXT_DIM};
    }}
    QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {{
        height: 0px;
    }}
    QScrollBar:horizontal {{
        background-color: {BG_DARK};
        height: 10px;
        border: none;
    }}
    QScrollBar::handle:horizontal {{
        background-color: {BORDER};
        border-radius: 5px;
        min-width: 30px;
    }}

    /* ---- Progress bar ---- */
    QProgressBar {{
        background-color: {BG_LIGHT};
        border: 1px solid {BORDER};
        border-radius: 4px;
        text-align: center;
        color: {TEXT_BRIGHT};
        font-size: {FONT_SIZE_SMALL}pt;
        min-height: 20px;
    }}
    QProgressBar::chunk {{
        background-color: {ACCENT};
        border-radius: 3px;
    }}

    /* ---- Tab widget ---- */
    QTabWidget::pane {{
        border: 1px solid {BORDER};
        background-color: {BG_MID};
        border-radius: 3px;
    }}
    QTabBar::tab {{
        background-color: {BG_LIGHT};
        color: {TEXT_DIM};
        border: 1px solid {BORDER};
        border-bottom: none;
        padding: 7px 16px;
        margin-right: 2px;
        border-top-left-radius: 4px;
        border-top-right-radius: 4px;
    }}
    QTabBar::tab:selected {{
        background-color: {BG_MID};
        color: {TEXT_BRIGHT};
        border-bottom: 2px solid {ACCENT};
    }}
    QTabBar::tab:hover:!selected {{
        background-color: {BG_HEADER};
        color: {TEXT};
    }}

    /* ---- Splitter ---- */
    QSplitter::handle {{
        background-color: {BORDER};
    }}
    QSplitter::handle:horizontal {{
        width: 1px;
    }}
    QSplitter::handle:vertical {{
        height: 1px;
    }}

    /* ---- Log panel ---- */
    QTextEdit[readOnly="true"] {{
        background-color: {PLOT_BG};
        color: {TEXT_DIM};
        border: 1px solid {BORDER};
        border-radius: 3px;
        font-family: "Consolas", "Liberation Mono", monospace;
        font-size: {FONT_SIZE_SMALL}pt;
    }}

    /* ---- Tooltips ---- */
    QToolTip {{
        background-color: {BG_HEADER};
        color: {TEXT};
        border: 1px solid {BORDER};
        padding: 4px;
    }}
    """
