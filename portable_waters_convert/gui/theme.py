"""Native obsidian/sapphire and parchment/sage glass themes for the desktop console."""
from __future__ import annotations

import atexit
from html import escape
import os
from random import Random
import sys
import tempfile

from PySide6.QtCore import QPointF, QRectF, QRect, Qt
from PySide6.QtGui import QColor, QImage, QPainter, QPixmap, QRadialGradient, QLinearGradient, QPen
from PySide6.QtWidgets import QGraphicsBlurEffect, QGraphicsPixmapItem, QGraphicsScene, QWidget


class ConsoleBackdrop(QWidget):
    """One cached scene and one blurred copy, shared by all glass panels."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.dark = True
        self._scene = self._frosted = None
        self.cache_builds = 0

    def set_dark(self, dark):
        if self.dark != bool(dark):
            self.dark = bool(dark)
            self._scene = self._frosted = None
            self.update()

    def resizeEvent(self, event):
        self._scene = self._frosted = None
        super().resizeEvent(event)

    def background(self, frosted=False):
        width, height = max(1, self.width()), max(1, self.height())
        if self._scene is None:
            ratio = self.devicePixelRatioF()
            image = QImage(round(width * ratio), round(height * ratio), QImage.Format_ARGB32_Premultiplied)
            image.setDevicePixelRatio(ratio)
            image.fill(QColor(_DARK["L0"] if self.dark else _LIGHT["L0"]))
            painter = QPainter(image)
            painter.setRenderHint(QPainter.Antialiasing)
            if self.dark:
                gradient = QLinearGradient(0, 0, width, height)
                for position, color in ((0, "#303f5b"), (.38, "#1b2b43"), (.74, "#101c2e"), (1, "#080c15")):
                    gradient.setColorAt(position, QColor(color))
                painter.fillRect(QRectF(0, 0, width, height), gradient)
            lights = ((.06, .02, .85, "#344967"), (.96, .98, .9, "#080b13")) if self.dark else (
                (.10, .01, .8, "#fff9e8"), (.88, .88, .8, "#cbd6c2"))
            for x, y, radius, color in lights:
                glow = QRadialGradient(QPointF(width * x, height * y), width * radius)
                glow.setColorAt(0, QColor(color))
                glow.setColorAt(1, QColor(0, 0, 0, 0))
                painter.fillRect(QRectF(0, 0, width, height), glow)
            # Large, quiet contours belong to the material beneath the glass.
            painter.setPen(QPen(QColor(171, 187, 224, 12) if self.dark else QColor(103, 128, 91, 14), 1))
            painter.drawEllipse(QRectF(-width * .55, -height * .7, width * 1.5, height * 1.6))
            painter.drawEllipse(QRectF(width * .42, height * .38, width * .9, height * 1.1))
            grain = QPixmap(64, 64)
            grain.fill(Qt.transparent)
            texture = QPainter(grain)
            random = Random(17)
            texture.setPen(QColor(215, 225, 255, 3) if self.dark else QColor(107, 99, 72, 2))
            for _ in range(384):
                texture.drawPoint(random.randrange(64), random.randrange(64))
            if not self.dark:
                # Linen threads are confined to the background; panel frosting
                # softens them and scientific surfaces remain opaque.
                texture.setPen(QPen(QColor(98, 91, 66, 2), .6))
                for position in (3, 12, 24, 35, 47, 58):
                    texture.drawLine(position, 0, position, 64)
                    texture.drawLine(0, position + 2, 64, position + 2)
            texture.end()
            painter.drawTiledPixmap(QRect(0, 0, width, height), grain)
            if self.dark:
                painter.setPen(QColor(212, 224, 255, 18))
                for _ in range(max(1, width * height // 60000)):
                    painter.drawPoint(random.randrange(width), random.randrange(height))
            painter.end()
            self._scene = QPixmap.fromImage(image)
            scene = QGraphicsScene()
            item = QGraphicsPixmapItem(self._scene)
            blur = QGraphicsBlurEffect()
            blur.setBlurRadius(24)
            item.setGraphicsEffect(blur)
            scene.addItem(item)
            self._frosted = QPixmap(image.size())
            self._frosted.setDevicePixelRatio(ratio)
            self._frosted.fill(QColor(_DARK["L0"] if self.dark else _LIGHT["L0"]))
            painter = QPainter(self._frosted)
            scene.render(painter, QRectF(0, 0, width, height), QRectF(0, 0, width, height))
            painter.end()
            self.cache_builds += 1
        return self._frosted if frosted else self._scene

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.drawPixmap(0, 0, self.background())
        painter.end()

# -- Accent + status colors (same in dark and light) --
ACCENT = "#8aa8f2"
ACCENT_HOVER = "#adc3ff"
ACCENT_PRESSED = "#718ed2"
ACCENT_DIM = "rgba(138,168,242,25)"
SUCCESS = "#77bbaa"
WARNING = "#e9bc74"
ERROR = "#e88883"

# -- Fonts --
FONT_FAMILY = "Segoe UI"
FONT_SIZE = 10
FONT_SIZE_SMALL = 9
FONT_SIZE_TITLE = 14
FONT_SIZE_SECTION = 11
FONT_SIZE_CAPTION = 8

# -- 4-layer dark palette --
_DARK = dict(
    L0="#0b0e14", L1="#171c27", L2="#111722", L3="#26334e",
    border="#303b51", border_hi="#657aa3",
    text="#e5edf3", text_dim="#9aabba", text_hi="#f5f8fb",
    plot_bg="#111722", plot_face="#111722", plot_hex="#171c27",
    grid="#313b4e", log_bg="#101923",
    card="rgba(45,59,82,208)", card_bottom="rgba(13,20,34,224)",
    card_highlight="#526382", chip="rgba(38,49,72,230)", chip_border="#4b5d82",
    header="#202a3b", table="#131b28", alternate="#182231",
    table_grid="#263248", selected="#2a3958", selected_border="#7f9edf",
    button="#273247", button_bottom="#1d2738", button_hover="#34425e",
    splitter="#2b3549", focus="#253550", hover="#1e2b41", accent=ACCENT,
    accent_hover=ACCENT_HOVER, accent_pressed=ACCENT_PRESSED,
    accent_dim=ACCENT_DIM, on_accent="#101a30", warning=WARNING, error=ERROR, success=SUCCESS,
)
_LIGHT = dict(
    L0="#eeeade", L1="#fffdf6", L2="#f5f2e8", L3="#e5eddf",
    border="#c9cebd", border_hi="#91a38b",
    text="#303d2e", text_dim="#5e6958", text_hi="#253e29",
    plot_bg="#fffdf6", plot_face="#fffdf6", plot_hex="#fffdf6",
    grid="#ced9d0", log_bg="#f7f5ef",
    card="rgba(255,253,246,216)", card_bottom="rgba(241,239,225,226)",
    card_highlight="#ffffff", chip="rgba(229,237,218,232)", chip_border="#bccdb3",
    header="#eaefe8", table="#fffefa", alternate="#f8f8f2",
    table_grid="#e7eae3", selected="#dcede3", selected_border="#76a88e",
    button="#fffefa", button_bottom="#f0f2eb", button_hover="#edf4ec",
    splitter="#d9dfd6", focus="#e7f4eb", hover="#f4f8f1", accent="#4b7454",
    accent_hover="#5d8663", accent_pressed="#3a6143",
    accent_dim="rgba(75,116,84,28)", on_accent="#fffefa",
    warning="#93601d", error="#b34039", success="#287357",
)

TEXT_DIM = _DARK["text_dim"]
TEXT = _DARK["text"]

# -- Checkmark PNG --
_checkmark_path: str | None = None
_icon_paths: dict[tuple[str, str], str] = {}


def _icon_png(kind: str, color: str) -> str:
    """Small native icons use the same writable, session-owned cache as ticks."""
    key = (kind, color)
    if key in _icon_paths and os.path.isfile(_icon_paths[key]):
        return _icon_paths[key]
    from PySide6.QtCore import QPointF, Qt
    from PySide6.QtGui import QColor, QPainter, QPen, QPixmap, QPolygonF

    pixmap = QPixmap(16, 16)
    pixmap.fill(Qt.transparent)
    painter = QPainter(pixmap)
    painter.setRenderHint(QPainter.Antialiasing)
    pen = QPen(QColor(color), 1.6)
    pen.setCapStyle(Qt.RoundCap)
    pen.setJoinStyle(Qt.RoundJoin)
    painter.setPen(pen)
    if kind in ("up", "down", "back", "forward"):
        points = {"up": [(4, 10), (8, 6), (12, 10)], "down": [(4, 6), (8, 10), (12, 6)],
                  "back": [(10, 3), (5, 8), (10, 13)], "forward": [(6, 3), (11, 8), (6, 13)]}[kind]
        painter.drawPolyline(QPolygonF([QPointF(x, y) for x, y in points]))
    elif kind == "parent":
        painter.drawPolyline(QPolygonF([QPointF(3, 7), QPointF(8, 2), QPointF(13, 7)]))
        painter.drawLine(8, 3, 8, 14)
    elif kind in ("list", "detail"):
        for y in (4, 8, 12):
            painter.drawPoint(3, y)
            painter.drawLine(6, y, 13, y)
            if kind == "detail":
                painter.drawLine(10, y - 1, 10, y + 1)
    elif kind == "ranges":
        for x, y, dx, dy in ((3, 3, 4, 4), (13, 3, -4, 4), (3, 13, 4, -4), (13, 13, -4, -4)):
            painter.drawLine(x, y, x + dx, y)
            painter.drawLine(x, y, x, y + dy)
    elif kind == "display":
        for y, x in ((4, 6), (8, 11), (12, 5)):
            painter.drawLine(3, y, 13, y)
            painter.drawLine(x, y - 2, x, y + 2)
    elif kind == "ccs":
        painter.drawLine(3, 13, 3, 3)
        painter.drawLine(3, 13, 13, 13)
        painter.drawPolyline(QPolygonF([QPointF(4, 11), QPointF(7, 8), QPointF(10, 6), QPointF(13, 4)]))
    elif kind == "acquisitions":
        for y in (4, 8, 12):
            painter.drawPoint(3, y)
            painter.drawLine(6, y, 13, y)
    else:
        painter.drawRoundedRect(2, 4, 12, 9, 1, 1)
        painter.drawLine(3, 4, 3, 2)
        painter.drawLine(3, 2, 7, 2)
        painter.drawLine(7, 2, 9, 4)
        if kind == "folder-add":
            painter.drawLine(6, 8, 10, 8)
            painter.drawLine(8, 6, 8, 10)
    painter.end()
    with tempfile.NamedTemporaryFile(prefix="dvv_icon_", suffix=".png", delete=False) as icon:
        path = icon.name
    pixmap.save(path, "PNG")
    _icon_paths[key] = path
    atexit.register(lambda p=path: os.remove(p) if os.path.isfile(p) else None)
    return path


def decorate_section(group, kind: str) -> None:
    """Use the existing GroupBox title space; do not add a layout row."""
    from PySide6.QtCore import Qt
    from PySide6.QtWidgets import QLabel

    if hasattr(group, "_section_chip"):
        return
    group.setProperty("glassSection", True)
    group.setProperty("sectionKind", kind)
    chip = QLabel(group)
    chip.setProperty("sectionChip", True)
    chip.setTextFormat(Qt.RichText)
    chip.setAttribute(Qt.WA_TransparentForMouseEvents)
    chip.setAccessibleName(group.title())
    group._section_chip = chip
    chip.move(12, 0)
    chip.show()


def refresh_console_style(root, dark: bool) -> None:
    """Retint native title icons and the shallow shadows on static cards."""
    from PySide6.QtCore import QUrl
    from PySide6.QtGui import QColor, QPalette
    from PySide6.QtWidgets import QApplication, QAbstractSpinBox, QGraphicsDropShadowEffect, QGroupBox, QLabel

    p = _DARK if dark else _LIGHT
    palette = QPalette()
    for role, key in ((QPalette.Window, "L1"), (QPalette.WindowText, "text"),
            (QPalette.Base, "table"), (QPalette.AlternateBase, "alternate"),
            (QPalette.Text, "text"), (QPalette.Button, "button"), (QPalette.ButtonText, "text"),
            (QPalette.Highlight, "selected"), (QPalette.HighlightedText, "text_hi"),
            (QPalette.ToolTipBase, "L1"), (QPalette.ToolTipText, "text"),
            (QPalette.PlaceholderText, "text_dim"), (QPalette.Link, "accent")):
        palette.setColor(role, QColor(p[key]))
    for role in (QPalette.Text, QPalette.WindowText, QPalette.ButtonText):
        palette.setColor(QPalette.Disabled, role, QColor(p["text_dim"]))
    QApplication.instance().setPalette(palette)
    effects = getattr(root, "_console_effects", None)
    if effects is not None:
        effects.set_dark(dark)
    for backdrop in root.findChildren(ConsoleBackdrop):
        backdrop.set_dark(dark)
    for spin in root.findChildren(QAbstractSpinBox):
        spin.setButtonSymbols(QAbstractSpinBox.UpDownArrows)
    for group in root.findChildren(QGroupBox):
        if not hasattr(group, "_section_chip"):
            continue
        path = QUrl.fromLocalFile(_icon_png(group.property("sectionKind"), p["accent"])).toString()
        chip = group._section_chip
        chip.setText(f'<img src="{path}" width="12" height="12">&nbsp; {escape(group.title())}')
        chip.adjustSize()
        chip.raise_()
        effect = group.graphicsEffect()
        if effect is None:
            effect = QGraphicsDropShadowEffect(group)
            effect.setBlurRadius(16)
            effect.setOffset(0, 3)
            group.setGraphicsEffect(effect)
        effect.setColor(QColor(0, 0, 0, 58 if dark else 26))
    # Existing panels use these named colours in a few inline status labels.
    # Retint them in place so warnings and captions stay readable in light mode.
    pairs = ((TEXT_DIM, _LIGHT["text_dim"], p["text_dim"]),
             (WARNING, _LIGHT["warning"], p["warning"]),
             (ERROR, _LIGHT["error"], p["error"]),
             (SUCCESS, _LIGHT["success"], p["success"]))
    for label in root.findChildren(QLabel):
        css = label.styleSheet()
        for dark_color, light_color, current in pairs:
            css = css.replace(dark_color, current).replace(light_color, current)
        if css != label.styleSheet():
            label.setStyleSheet(css)

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

    # Installed application folders may be read-only. Give every process its
    # own writable icon file so simultaneous app sessions cannot delete it.
    with tempfile.NamedTemporaryFile(prefix="dvv_checkmark_", suffix=".png", delete=False) as icon:
        path = icon.name
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
    A, AH, AP, AD = p["accent"], p["accent_hover"], p["accent_pressed"], p["accent_dim"]
    heading_font = f'"{FONT_FAMILY}", "Liberation Sans", sans-serif' if dark else '"Georgia", "Liberation Serif", serif'
    chk = _ensure_checkmark_png().replace("\\", "/")
    up = _icon_png("up", p["text"]).replace("\\", "/")
    down = _icon_png("down", p["text"]).replace("\\", "/")
    up_dim = _icon_png("up", p["text_dim"]).replace("\\", "/")
    down_dim = _icon_png("down", p["text_dim"]).replace("\\", "/")
    return f"""

    /* ===== BASE ===== */
    QWidget {{
        background-color: transparent;
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
        background: qlineargradient(x1:0, y1:0, x2:0, y2:1,
                                   stop:0 {p['button']}, stop:1 {p['button_bottom']});
        color: {p['text']};
        border: 1px solid {p['border']};
        border-radius: 6px;
        padding: 7px 23px;
        font-weight: 600;
        font-size: {FONT_SIZE}pt;
        min-height: 20px;
    }}
    QPushButton:hover {{
        background: qlineargradient(x1:0, y1:0, x2:0, y2:1, stop:0 {p['button']}, stop:1 {p['button_bottom']});
        border-color: {p['border_hi']};
    }}
    QPushButton:pressed {{ border-color: {p['border_hi']}; }}
    QPushButton:focus {{ border-color: {A}; }}
    QPushButton:disabled {{
        background: {p['L2']}; color: {p['text_dim']};
        border-color: {p['table_grid']};
    }}
    QPushButton[secondary="true"] {{ font-weight: normal; }}
    QPushButton[accentOutline="true"] {{
        background: {AD}; color: {A}; border-color: {p['selected_border']};
    }}
    QPushButton[accentOutline="true"]:hover {{
        background: {AD}; border-color: {A};
    }}
    QPushButton[primary="true"] {{
        background: {A}; color: {p['on_accent']}; border-color: {A};
    }}
    QPushButton[primary="true"]:hover {{ background: {A}; border-color: {AH}; }}
    QPushButton[primary="true"]:pressed {{ background: {A}; border-color: {AP}; }}
    QPushButton[primary="true"]:disabled, QPushButton[accentOutline="true"]:disabled {{
        background: {p['L2']}; color: {p['text_dim']}; border-color: {p['table_grid']};
    }}

    /* ===== INPUTS ===== */
    QLineEdit, QSpinBox, QDoubleSpinBox, QComboBox {{
        background-color: {p['L2']};
        color: {p['text']};
        border: 1px solid {p['border']};
        border-radius: 6px;
        padding: 6px 10px;
        font-size: {FONT_SIZE}pt;
        selection-background-color: {p['selected']};
    }}
    QLineEdit:focus, QSpinBox:focus, QDoubleSpinBox:focus, QComboBox:focus {{
        border: 1px solid {A};
        background-color: {p['focus']};
        selection-background-color: {A};
        selection-color: {p['on_accent']};
    }}
    QLineEdit:hover:!focus, QSpinBox:hover:!focus,
    QDoubleSpinBox:hover:!focus, QComboBox:hover:!focus {{
        border-color: {p['border_hi']};
    }}
    QSpinBox, QDoubleSpinBox {{ padding-right: 22px; }}
    QSpinBox::up-button, QDoubleSpinBox::up-button {{
        subcontrol-origin: border;
        subcontrol-position: top right;
        width: 18px;
        border: 1px solid {p['border']};
        border-left: none;
        border-bottom: none;
        border-top-right-radius: 5px;
        background-color: {p['L1']};
    }}
    QSpinBox::down-button, QDoubleSpinBox::down-button {{
        subcontrol-origin: border;
        subcontrol-position: bottom right;
        width: 18px;
        border: 1px solid {p['border']};
        border-left: none;
        border-bottom-right-radius: 5px;
        background-color: {p['L1']};
    }}
    QSpinBox::up-button:hover, QDoubleSpinBox::up-button:hover,
    QSpinBox::down-button:hover, QDoubleSpinBox::down-button:hover {{
        background-color: {p['L3']};
    }}
    QSpinBox::up-arrow, QDoubleSpinBox::up-arrow {{
        image: url("{up}"); width: 12px; height: 12px;
    }}
    QSpinBox::down-arrow, QDoubleSpinBox::down-arrow {{
        image: url("{down}"); width: 12px; height: 12px;
    }}
    QSpinBox::up-arrow:disabled, QDoubleSpinBox::up-arrow:disabled {{ image: url("{up_dim}"); }}
    QSpinBox::down-arrow:disabled, QDoubleSpinBox::down-arrow:disabled {{ image: url("{down_dim}"); }}
    QComboBox::down-arrow {{ image: url("{down}"); width: 12px; height: 12px; }}
    QComboBox::drop-down {{
        subcontrol-origin: border;
        subcontrol-position: top right;
        width: 22px;
        border: none;
    }}
    QGroupBox#displayControls QComboBox {{ padding: 6px; }}
    QGroupBox#displayControls QComboBox::drop-down {{ width: 18px; }}
    QComboBox QAbstractItemView {{
        background-color: {p['L2']};
        color: {p['text']};
        selection-background-color: {A};
        border: 1px solid {p['border']};
        border-radius: 4px;
    }}

    QDialog {{
        background: qlineargradient(x1:0, y1:0, x2:1, y2:1, stop:0 {p['L1']}, stop:1 {p['L2']});
        color: {p['text']}; border: 1px solid {p['border_hi']};
    }}
    QMenu {{
        background: qlineargradient(x1:0, y1:0, x2:0, y2:1, stop:0 {p['L1']}, stop:1 {p['L2']});
        color: {p['text']}; margin: 14px;
        border: 1px solid {p['border_hi']}; border-radius: 8px; padding: 6px;
    }}
    QMenu::item {{ padding: 7px 24px 7px 12px; border-radius: 4px; }}
    QMenu::item:selected {{ background: {p['selected']}; color: {p['text_hi']}; }}
    QMenu::item:disabled {{ color: {p['text_dim']}; }}
    QMenu::separator {{ height: 1px; background: {p['border']}; margin: 5px 6px; }}

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
        background: qlineargradient(x1:0, y1:0, x2:0.8, y2:1,
                                   stop:0 {p['card']}, stop:1 {p['card_bottom']});
        border: 1px solid {p['border']};
        border-top-color: {p['card_highlight']};
        border-radius: 10px;
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
    QGroupBox[glassSection="true"]::title {{ color: transparent; }}
    QLabel[sectionChip="true"] {{
        background-color: {p['chip']};
        color: {p['text_hi']};
        border: 1px solid {p['chip_border']};
        border-radius: 6px;
        padding: 3px 7px;
        font-weight: 600;
        font-size: 10pt;
        font-family: {heading_font};
    }}

    /* ===== TABLE ===== */
    QTableWidget {{
        background-color: {p['table']};
        alternate-background-color: {p['alternate']};
        gridline-color: transparent;
        border: 1px solid {p['border']};
        border-radius: 6px;
        font-size: {FONT_SIZE_SMALL}pt;
    }}
    QTableWidget::item {{ padding: 3px 8px; border-bottom: 1px solid {p['table_grid']}; }}
    QTableWidget::item:hover {{ background-color: {p['hover']}; }}
    QTableWidget::item:selected {{ background-color: {p['selected']}; color: {p['text_hi']}; }}
    QTableWidget::item:selected:!active {{ background-color: {p['selected']}; color: {p['text_hi']}; }}
    QTableCornerButton::section {{
        background-color: {p['header']};
        border: none;
        border-bottom: 1px solid {p['border']};
    }}
    QTableWidget QLineEdit {{
        padding: 1px 5px;
        border-radius: 4px;
        font-size: {FONT_SIZE}pt;
    }}
    QHeaderView::section {{
        background-color: {p['header']};
        color: {p['text']};
        border: none;
        border-right: 1px solid {p['table_grid']};
        border-bottom: 1px solid {p['border']};
        padding: 6px 10px;
        font-weight: 600;
        font-size: {FONT_SIZE_SMALL}pt;
    }}

    /* ===== SCROLLBAR ===== */
    QScrollBar:vertical {{
        background: transparent;
        width: 11px;
        border: none;
        margin: 4px 2px;
    }}
    QScrollBar::handle:vertical {{
        background-color: {p['border_hi']};
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
        background-color: {p['card']};
        border-radius: 8px;
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
        background: {p['button_bottom']};
        color: {p['text_hi']};
        font-weight: 600;
        border-bottom: 2px solid {A};
    }}
    QTabBar::tab:hover:!selected {{
        color: {p['text']};
        border-bottom: 2px solid {p['border']};
    }}

    /* ===== SPLITTER ===== */
    QSplitter::handle {{
        background-color: {p['splitter']};
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
        width: 16px; height: 16px;
        border: 1px solid {p['border']};
        border-radius: 4px;
        background-color: {p['L2']};
    }}
    QCheckBox::indicator:hover {{ border-color: {A}; }}
    QCheckBox::indicator:checked {{
        border-color: {A};
        background-color: {A};
        image: url("{chk}");
    }}
    QGroupBox#axisRanges QCheckBox::indicator {{
        width: 14px; height: 14px;
        border-width: 1px;
        border-radius: 3px;
    }}

    QWidget#appHeader {{
        background-color: {p['card']};
        border-bottom: 1px solid {p['border']};
    }}
    QLabel#appBrand {{
        font-family: {heading_font}; color: {A}; font-size: 14pt; font-weight: 700; }}
    QLabel#acquisitionTitle {{ color: {p['text_hi']}; font-size: 12pt; font-weight: 600; font-family: {heading_font}; }}
    QLabel#panelCaption {{ color: {p['text_hi']}; font-size: 10pt; font-weight: 600; }}
    QToolBar#watersPlotTools {{
        background: transparent;
        border: none;
        border-top: 1px solid {p['table_grid']};
        spacing: 3px;
    }}
    QToolBar#watersPlotTools QToolButton {{
        background: transparent;
        border: 1px solid transparent;
        border-radius: 5px;
        padding: 4px;
    }}
    QToolBar#watersPlotTools QToolButton:hover {{
        background: {p['hover']}; border-color: {p['border']};
    }}
    QToolBar#watersPlotTools QToolButton:checked {{
        background: {AD}; border-color: {p['selected_border']};
    }}
    QLabel#appStatus {{ color: {p['text_dim']}; font-size: 9pt; }}
    QLabel#appStatus[statusLevel="warning"] {{ color: {p['warning']}; }}
    QLabel#appStatus[statusLevel="error"] {{ color: {p['error']}; }}
    QListWidget#modeNav {{
        background: transparent;
        border: none;
        outline: none;
        font-size: 9pt;
    }}
    QListWidget#modeNav::item {{
        background: transparent;
        color: {p['text_dim']};
        padding: 7px 14px;
        margin: 4px;
        border-radius: 5px;
    }}
    QListWidget#modeNav::item:selected {{
        background: {p['L3']};
        color: {A};
    }}
    QWidget#appFooter {{
        border-top: 1px solid {p['border']};
        background: {p['L0']};
    }}
    QWidget#consoleCard, QFrame#consoleCard {{
        background: qlineargradient(x1:0, y1:0, x2:0, y2:1,
                                   stop:0 {p['card']}, stop:1 {p['card_bottom']});
        border: 1px solid {p['border']};
        border-radius: 8px;
    }}
    QLabel[eyebrow="true"] {{
        color: {p['text_dim']};
        font-size: 8pt;
        font-weight: 600;
    }}
    QLabel[chip="true"] {{
        background: {p['L3']};
        color: {A};
        border-radius: 10px;
        padding: 5px 10px;
        font-size: 9pt;
    }}
    QPushButton[compact="true"] {{
        padding: 4px 10px;
        min-height: 14px;
        font-size: 10pt;
    }}
    """
