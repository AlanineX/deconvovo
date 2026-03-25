"""Theme for DeconVoVo GUI — dark and light modes."""
from __future__ import annotations

# -- Shared accent colors --
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

# -- Palettes (Qt uses rgba for translucency, mpl needs hex) --
_DARK = dict(
    bg="rgba(30,30,46,230)", bg_mid="rgba(37,37,56,240)", bg_mid_hex="#252538",
    bg_input="#2e2e44", bg_head="#363654", border="#3e3e5e",
    text="#d4d4e8", text_dim="#8888aa", text_hi="#f0f0ff",
    plot_bg="#1a1a2e", plot_face="#22223a", grid="#3a3a5a", log_bg="#1a1a2e",
)
_LIGHT = dict(
    bg="rgba(244,244,248,230)", bg_mid="rgba(255,255,255,240)", bg_mid_hex="#ffffff",
    bg_input="#e8e8f0", bg_head="#dcdce8", border="#c0c0d0",
    text="#2a2a3a", text_dim="#6a6a80", text_hi="#1a1a2a",
    plot_bg="#f8f8fc", plot_face="#ffffff", grid="#d0d0e0", log_bg="#f0f0f6",
)

# Stable module-level refs for panels that import color constants
TEXT_DIM = _DARK["text_dim"]
TEXT = _DARK["text"]


def mpl_rc(dark=True):
    p = _DARK if dark else _LIGHT
    return {
        "figure.facecolor": p["plot_bg"], "axes.facecolor": p["plot_face"],
        "axes.edgecolor": p["border"], "axes.labelcolor": p["text"],
        "axes.titlecolor": p["text_hi"], "xtick.color": p["text_dim"],
        "ytick.color": p["text_dim"], "text.color": p["text"],
        "grid.color": p["grid"], "grid.alpha": 0.4,
        "font.family": "sans-serif",
        "font.sans-serif": ["Segoe UI", "Liberation Sans", "Arial"],
        "font.size": 10, "legend.facecolor": p["bg_mid_hex"],
        "legend.edgecolor": p["border"], "figure.dpi": 100,
    }


def stylesheet(dark=True):
    p = _DARK if dark else _LIGHT
    A, AH, AP = ACCENT, ACCENT_HOVER, ACCENT_PRESSED
    FS, FSS, FST, FSH = FONT_SIZE, FONT_SIZE_SMALL, FONT_SIZE_TITLE, FONT_SIZE_HEADER
    FF = FONT_FAMILY
    return f"""
    QWidget {{ background-color:{p['bg']}; color:{p['text']};
        font-family:"{FF}","Liberation Sans","Arial",sans-serif; font-size:{FS}pt; }}
    QMainWindow {{ background-color:{p['bg']}; }}

    QListWidget {{ background-color:{p['bg_mid']}; border:none;
        border-right:1px solid {p['border']}; outline:none; padding:4px 0; font-size:{FSH}pt; }}
    QListWidget::item {{ padding:12px 16px; border:none; color:{p['text_dim']}; }}
    QListWidget::item:selected {{ background-color:{A}; color:{p['text_hi']};
        border-left:3px solid {AH}; }}
    QListWidget::item:hover:!selected {{ background-color:{p['bg_head']}; color:{p['text']}; }}

    QPushButton {{ background-color:{A}; color:{p['text_hi']}; border:none;
        border-radius:4px; padding:8px 20px; font-weight:bold; font-size:{FS}pt; min-height:18px; }}
    QPushButton:hover {{ background-color:{AH}; }}
    QPushButton:pressed {{ background-color:{AP}; }}
    QPushButton:disabled {{ background-color:{p['bg_input']}; color:{p['text_dim']}; }}
    QPushButton[secondary="true"] {{ background-color:{p['bg_input']}; color:{p['text']};
        border:1px solid {p['border']}; }}
    QPushButton[secondary="true"]:hover {{ background-color:{p['bg_head']}; }}

    QLineEdit, QSpinBox, QDoubleSpinBox, QComboBox {{ background-color:{p['bg_input']};
        color:{p['text']}; border:1px solid {p['border']}; border-radius:3px;
        padding:5px 8px; font-size:{FS}pt; }}
    QLineEdit:focus, QSpinBox:focus, QDoubleSpinBox:focus, QComboBox:focus {{
        border:1px solid {A}; }}
    QComboBox::drop-down {{ border:none; padding-right:8px; }}
    QComboBox QAbstractItemView {{ background-color:{p['bg_input']}; color:{p['text']};
        selection-background-color:{A}; border:1px solid {p['border']}; }}

    QLabel {{ background-color:transparent; color:{p['text']}; }}
    QLabel[title="true"] {{ font-size:{FST}pt; font-weight:bold; color:{p['text_hi']}; }}
    QLabel[subtitle="true"] {{ font-size:{FSS}pt; color:{p['text_dim']}; }}

    QGroupBox {{ background-color:{p['bg_mid']}; border:1px solid {p['border']};
        border-radius:5px; margin-top:12px; padding-top:18px;
        font-weight:bold; color:{p['text_hi']}; }}
    QGroupBox::title {{ subcontrol-origin:margin; left:12px; padding:0 6px; }}

    QTableWidget {{ background-color:{p['bg_mid']}; alternate-background-color:{p['bg_input']};
        gridline-color:{p['border']}; border:1px solid {p['border']};
        border-radius:3px; font-size:{FSS}pt; }}
    QTableWidget::item {{ padding:3px 6px; }}
    QTableWidget::item:selected {{ background-color:{A}; color:{p['text_hi']}; }}
    QHeaderView::section {{ background-color:{p['bg_head']}; color:{p['text_hi']};
        border:none; border-right:1px solid {p['border']};
        border-bottom:1px solid {p['border']}; padding:5px 8px;
        font-weight:bold; font-size:{FSS}pt; }}

    QScrollBar:vertical {{ background-color:{p['bg']}; width:10px; border:none; }}
    QScrollBar::handle:vertical {{ background-color:{p['border']}; border-radius:5px; min-height:30px; }}
    QScrollBar::handle:vertical:hover {{ background-color:{p['text_dim']}; }}
    QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {{ height:0; }}
    QScrollBar:horizontal {{ background-color:{p['bg']}; height:10px; border:none; }}
    QScrollBar::handle:horizontal {{ background-color:{p['border']}; border-radius:5px; min-width:30px; }}

    QProgressBar {{ background-color:{p['bg_input']}; border:1px solid {p['border']};
        border-radius:4px; text-align:center; color:{p['text_hi']};
        font-size:{FSS}pt; min-height:20px; }}
    QProgressBar::chunk {{ background-color:{A}; border-radius:3px; }}

    QTabWidget::pane {{ border:1px solid {p['border']}; background-color:{p['bg_mid']};
        border-radius:3px; }}
    QTabBar::tab {{ background-color:{p['bg_input']}; color:{p['text_dim']};
        border:1px solid {p['border']}; border-bottom:none; padding:7px 16px;
        margin-right:2px; border-top-left-radius:4px; border-top-right-radius:4px; }}
    QTabBar::tab:selected {{ background-color:{p['bg_mid']}; color:{p['text_hi']};
        border-bottom:2px solid {A}; }}
    QTabBar::tab:hover:!selected {{ background-color:{p['bg_head']}; color:{p['text']}; }}

    QSplitter::handle {{ background-color:{p['border']}; }}
    QSplitter::handle:horizontal {{ width:1px; }}
    QSplitter::handle:vertical {{ height:1px; }}

    QTextEdit[readOnly="true"] {{ background-color:{p['log_bg']}; color:{p['text_dim']};
        border:1px solid {p['border']}; border-radius:3px;
        font-family:"Consolas","Liberation Mono",monospace; font-size:{FSS}pt; }}

    QToolTip {{ background-color:{p['bg_head']}; color:{p['text']};
        border:1px solid {p['border']}; padding:4px; }}

    QCheckBox {{ background-color:transparent; spacing:6px; }}
    QCheckBox::indicator {{ width:18px; height:18px; border:2px solid {p['border']};
        border-radius:4px; background-color:{p['bg_input']}; }}
    QCheckBox::indicator:hover {{ border-color:{A}; }}
    QCheckBox::indicator:checked {{ border-color:{A}; background-color:{A}; }}
    """
