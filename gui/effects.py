"""Native glass materials and short, local interaction feedback."""
from __future__ import annotations

from PySide6.QtCore import QEasingCurve, QEvent, QPoint, QPointF, QRectF, QSize, Qt, QVariantAnimation, QObject
from PySide6.QtGui import QColor, QIcon, QLinearGradient, QPainter, QPainterPath, QPen
from PySide6.QtWidgets import (
    QAbstractButton, QAbstractSpinBox, QApplication, QComboBox, QGraphicsOpacityEffect,
    QDialog, QGraphicsDropShadowEffect, QGroupBox, QLineEdit, QMenu, QStyle, QToolButton,
    QStyledItemDelegate, QTableWidget, QWidget,
)
from shiboken6 import isValid

from gui.theme import ConsoleBackdrop, _DARK, _LIGHT, _icon_png


class GlassSurface(QWidget):
    """Sample the shared frosted scene without capturing or blurring data widgets."""

    def __init__(self, group, backdrop):
        super().__init__(group)
        self.group, self.backdrop = group, backdrop
        self.setAttribute(Qt.WA_TransparentForMouseEvents)
        self.setAttribute(Qt.WA_NoSystemBackground)
        ancestor = group
        while ancestor is not None and ancestor is not backdrop:
            ancestor.installEventFilter(self)
            ancestor = ancestor.parentWidget()
        self.setGeometry(group.rect())
        self.lower()
        self.show()

    def eventFilter(self, obj, event):
        if not isValid(self.group) or not isValid(self.backdrop) or QApplication.closingDown():
            return False
        if event.type() in (QEvent.Resize, QEvent.Move, QEvent.Show):
            self.setGeometry(self.group.rect())
            self.update()
        return False

    def paintEvent(self, event):
        p = _DARK if self.backdrop.dark else _LIGHT
        framed = isinstance(self.group, QGroupBox)
        top = 16 if framed and self.group.title() else 0
        radius = 9 if framed else 0
        rect = QRectF(self.rect()).adjusted(1, top + 1, -1, -1)
        if rect.isEmpty():
            return
        painter = QPainter(self)
        painter.setRenderHint(QPainter.Antialiasing)
        path = QPainterPath()
        path.addRoundedRect(rect, radius, radius)
        painter.setClipPath(path)
        origin = self.mapTo(self.backdrop, QPoint())
        painter.drawPixmap(-origin.x(), -origin.y(), self.backdrop.background(frosted=True))
        tint = QLinearGradient(rect.topLeft(), rect.bottomRight())
        tint.setColorAt(0, QColor(45, 59, 82, 208) if self.backdrop.dark else QColor(255, 253, 246, 216))
        tint.setColorAt(1, QColor(13, 20, 34, 224) if self.backdrop.dark else QColor(241, 239, 225, 226))
        painter.fillPath(path, tint)
        # Re-sample only the narrow glass edge. Data surfaces are never warped.
        inner = QPainterPath()
        inner.addRoundedRect(rect.adjusted(2, 2, -2, -2), max(0, radius - 2), max(0, radius - 2))
        if framed:
            painter.save()
            painter.setClipPath(path.subtracted(inner))
            painter.setOpacity(.18 if self.backdrop.dark else .12)
            painter.drawPixmap(-origin.x() + 1, -origin.y() + 2, self.backdrop.background(frosted=True))
            painter.restore()
        sheen = QLinearGradient(rect.topLeft(), rect.bottomRight())
        sheen.setColorAt(0, QColor(187, 205, 255, 86) if self.backdrop.dark else QColor(255, 255, 255, 205))
        sheen.setColorAt(.35, QColor(187, 205, 255, 15) if self.backdrop.dark else QColor(255, 255, 255, 38))
        sheen.setColorAt(1, QColor(0, 0, 0, 25 if self.backdrop.dark else 13))
        if framed:
            painter.setPen(QPen(sheen, 1))
            painter.drawRoundedRect(rect.adjusted(1, 1, -1, -1), max(0, radius - 1), max(0, radius - 1))
        painter.setPen(QPen(QColor(p["border"]), .6))
        painter.drawLine(rect.bottomLeft() + QPointF(9, 0), rect.bottomRight() - QPointF(9, 0))
        painter.end()


class SelectionIndicator(QStyledItemDelegate):
    def __init__(self, controller, parent):
        super().__init__(parent)
        self.controller = controller

    def paint(self, painter, option, index):
        super().paint(painter, option, index)
        if index.column() == 0 and option.state & QStyle.State_Selected:
            painter.save()
            painter.fillRect(option.rect.adjusted(1, 3, -option.rect.width() + 4, -3),
                             QColor(self.controller.palette["accent"]))
            painter.restore()


class FloatingEffect(QGraphicsDropShadowEffect):
    """Keep the native shadow and content fade in one Qt graphics effect."""

    def __init__(self, parent, dark):
        super().__init__(parent)
        self.opacity = 1.
        self.setBlurRadius(20)
        self.setOffset(0, 5)
        self.set_dark(dark)

    def set_dark(self, dark):
        self.setColor(QColor(0, 0, 0, 130 if dark else 65))

    def setOpacity(self, value):
        self.opacity = float(value)
        self.update()

    def draw(self, painter):
        painter.save()
        painter.setOpacity(painter.opacity() * self.opacity)
        super().draw(painter)
        painter.restore()


class InteractionLayer(QWidget):
    """One transparent overlay; animations repaint only the affected controls."""

    def __init__(self, root, dark):
        super().__init__(root)
        self.dark = dark
        self.targets, self.values, self.animations = {}, {}, {}
        self.pulse_warning = False
        self.setAttribute(Qt.WA_TransparentForMouseEvents)
        self.setAttribute(Qt.WA_NoSystemBackground)
        self.setGeometry(root.rect())
        self.show()
        self.raise_()

    def bounds(self, widget):
        if widget is None or not isValid(widget) or not widget.isVisible():
            return QRectF()
        origin = self.mapFromGlobal(widget.mapToGlobal(QPoint()))
        return QRectF(origin.x(), origin.y(), widget.width(), widget.height())

    def animate(self, channel, target, end, duration):
        old = self.targets.get(channel)
        if old is target and self.animations.get(channel) is not None:
            animation = self.animations[channel]
            if animation.endValue() == end and animation.state() == QVariantAnimation.Running:
                return
        dirty = self.bounds(old).adjusted(-7, -7, 7, 7)
        start = self.values.get(channel, 0.) if old is target else 0.
        self.targets[channel] = target
        animation = self.animations.get(channel)
        if animation is None:
            animation = QVariantAnimation(self)
            animation.valueChanged.connect(lambda value, key=channel: self._advance(key, value))
            self.animations[channel] = animation
        animation.stop()
        animation.setDuration(duration)
        animation.setEasingCurve(QEasingCurve.OutCubic)
        animation.setStartValue(start)
        animation.setEndValue(end)
        animation.start()
        self.update(dirty.toAlignedRect())

    def _advance(self, channel, value):
        self.values[channel] = float(value)
        self.update(self.bounds(self.targets.get(channel)).adjusted(-7, -7, 7, 7).toAlignedRect())

    def pulse(self, widget, warning=False):
        self.pulse_warning = warning
        self.animate("pulse", widget, 0., 420)
        animation = self.animations["pulse"]
        animation.stop()
        animation.setKeyValues([(0., 0.), (.18, 1.), (1., 0.)])
        animation.start()

    def paintEvent(self, event):
        p = _DARK if self.dark else _LIGHT
        painter = QPainter(self)
        painter.setRenderHint(QPainter.Antialiasing)
        for channel, widget in self.targets.items():
            value = self.values.get(channel, 0.)
            rect = self.bounds(widget).adjusted(.5, .5, -.5, -.5)
            if value <= .001 or rect.isEmpty() or not rect.intersects(QRectF(event.rect()).adjusted(-7, -7, 7, 7)):
                continue
            color = QColor(p["warning"] if channel == "pulse" and self.pulse_warning else p["accent"])
            if channel in ("hover", "press"):
                surface = QPainterPath()
                surface.addRoundedRect(rect.adjusted(1, 1, -1, -1), 5, 5)
                wash = QLinearGradient(rect.topLeft(), rect.bottomLeft())
                if channel == "press":
                    wash.setColorAt(0, QColor(0, 0, 0, round(38 * value)))
                    wash.setColorAt(1, QColor(0, 0, 0, round(10 * value)))
                else:
                    wash_color = QColor("#ffffff" if self.dark else p["accent"])
                    wash_color.setAlpha(round(18 * value))
                    wash.setColorAt(0, wash_color)
                    wash_color.setAlpha(round(3 * value))
                    wash.setColorAt(1, wash_color)
                painter.fillPath(surface, wash)
            if channel == "press":
                color = QColor(0, 0, 0, round(100 * value))
                painter.setPen(QPen(color, 1.4))
                painter.drawRoundedRect(rect.adjusted(1, 1, -1, -1), 5, 5)
                continue
            strength = 1. if channel == "focus" else .38 if channel == "hover" else .55
            for width, alpha in ((10, 5), (7, 10), (4, 24), (1.2, 155)):
                color.setAlpha(round(alpha * strength * value))
                painter.setPen(QPen(color, width))
                painter.drawRoundedRect(rect, 6, 6)
        painter.end()


class ConsoleEffects(QObject):
    def __init__(self, root, dark):
        super().__init__(root)
        self.root, self.dark = root, dark
        self.palette = _DARK if dark else _LIGHT
        self.surfaces = []
        self.popups = []
        self.dialog_layers = []
        backdrop = root.centralWidget()
        for group in root.findChildren(QGroupBox):
            if group.property("glassSection") or group.objectName() in ("watersPlotCard", "consoleCard"):
                self.surfaces.append(GlassSurface(group, backdrop))
        for strip in (root.findChild(QWidget, "appHeader"), root.findChild(QWidget, "appFooter")):
            if strip is not None:
                self.surfaces.append(GlassSurface(strip, backdrop))
        for menu in root.findChildren(QMenu):
            self._prepare_menu(menu)
        self.layer = InteractionLayer(backdrop, dark)
        for table in root.findChildren(QTableWidget):
            font = table.font()
            if hasattr(font, "setFeature"):
                font.setFeature(font.Tag.fromString("tnum"), 1)
                table.setFont(font)
            if type(table.itemDelegate()) is QStyledItemDelegate:
                table.setItemDelegate(SelectionIndicator(self, table))
        QApplication.instance().installEventFilter(self)

    def set_dark(self, dark):
        self.dark = bool(dark)
        self.palette = _DARK if dark else _LIGHT
        self.layer.dark = self.dark
        self.layer.update()
        for surface in self.surfaces:
            surface.update()
        for popup in self.popups:
            if isValid(popup) and isinstance(popup.graphicsEffect(), FloatingEffect):
                popup.graphicsEffect().set_dark(self.dark)
                if isinstance(popup, QDialog):
                    self._retint_dialog_icons(popup)
        for layer in self.dialog_layers:
            if isValid(layer):
                layer.dark = self.dark
                layer.update()

    def pulse(self, widget, warning=False):
        self.layer.pulse(widget, warning)

    def _control(self, obj):
        if not isinstance(obj, QWidget) or not isValid(obj) or not obj.isEnabled():
            return None
        window = obj.window()
        if window is not self.root and not hasattr(window, "_material_layer"):
            return None
        if isinstance(obj, QLineEdit) and isinstance(obj.parentWidget(), (QAbstractSpinBox, QComboBox)):
            obj = obj.parentWidget()
        return obj if isinstance(obj, (QAbstractButton, QAbstractSpinBox, QComboBox, QLineEdit)) else None

    def shutdown(self):
        QApplication.instance().removeEventFilter(self)
        for layer in [self.layer, *self.dialog_layers]:
            if isValid(layer):
                for animation in layer.animations.values():
                    animation.stop()

    def _owns_popup(self, popup):
        owner = popup.parent()
        while owner is not None and owner is not self.root:
            owner = owner.parent()
        return owner is self.root

    @staticmethod
    def _prepare_menu(menu):
        menu.setAttribute(Qt.WA_TranslucentBackground)
        menu.setWindowFlag(Qt.NoDropShadowWindowHint)

    def _retint_dialog_icons(self, dialog):
        kinds = {"backButton": "back", "forwardButton": "forward", "toParentButton": "parent",
                 "newFolderButton": "folder-add", "listModeButton": "list", "detailModeButton": "detail"}
        for button in dialog.findChildren(QToolButton):
            kind = kinds.get(button.objectName())
            if kind is None:
                continue
            icon = QIcon()
            for mode, color in ((QIcon.Normal, self.palette["text"]), (QIcon.Active, self.palette["accent"]),
                                (QIcon.Disabled, self.palette["border_hi"])):
                icon.addFile(_icon_png(kind, color), QSize(16, 16), mode)
            button.setIcon(icon)
            button.setIconSize(QSize(16, 16))

    def eventFilter(self, obj, event):
        if not isValid(self.root) or QApplication.closingDown():
            return False
        kind = event.type()
        if kind == QEvent.Polish and isinstance(obj, QMenu) and self._owns_popup(obj):
            self._prepare_menu(obj)
        if obj is self.root.centralWidget() and kind == QEvent.Resize:
            self.layer.setGeometry(obj.rect())
            self.layer.raise_()
            for surface in self.surfaces:
                surface.update()
        if kind == QEvent.Resize and hasattr(obj, "_material_layer"):
            obj._material_layer.setGeometry(obj.rect())
            obj._material_layer.raise_()
        if kind in (QEvent.Hide, QEvent.EnabledChange) and isinstance(obj, QWidget):
            if kind == QEvent.Hide or not obj.isEnabled():
                for layer in [self.layer, *self.dialog_layers]:
                    if isValid(layer):
                        for channel in ("hover", "focus", "press"):
                            if layer.targets.get(channel) is obj:
                                layer.animate(channel, obj, 0., 120)
        if (kind == QEvent.Show and isinstance(obj, QWidget)
                and (obj.windowType() == Qt.Popup or isinstance(obj, QDialog)) and self._owns_popup(obj)):
            if obj.graphicsEffect() is None:
                effect = FloatingEffect(obj, self.dark) if isinstance(obj, (QMenu, QDialog)) else QGraphicsOpacityEffect(obj)
                obj.setGraphicsEffect(effect)
                self.popups.append(obj)
                animation = QVariantAnimation(obj)
                animation.setDuration(140)
                animation.setEasingCurve(QEasingCurve.OutCubic)
                animation.valueChanged.connect(effect.setOpacity)
                obj._material_fade = animation
            if hasattr(obj, "_material_fade"):
                obj._material_fade.stop()
                obj._material_fade.setStartValue(.15)
                obj._material_fade.setEndValue(1.)
                obj._material_fade.start()
            if isinstance(obj, QDialog):
                if not hasattr(obj, "_material_layer"):
                    obj._material_layer = InteractionLayer(obj, self.dark)
                    self.dialog_layers.append(obj._material_layer)
                self._retint_dialog_icons(obj)
        if kind not in (QEvent.Enter, QEvent.Leave, QEvent.FocusIn, QEvent.FocusOut,
                        QEvent.MouseButtonPress, QEvent.MouseButtonRelease, QEvent.KeyPress, QEvent.KeyRelease):
            return False
        control = self._control(obj)
        if control is None:
            return False
        layer = getattr(control.window(), "_material_layer", self.layer)
        if kind == QEvent.Enter and control is obj:
            layer.animate("hover", control, 1., 120)
        elif kind == QEvent.Leave and control is obj and layer.targets.get("hover") is control:
            layer.animate("hover", control, 0., 120)
        elif kind == QEvent.FocusIn:
            layer.animate("focus", control, 1., 160)
        elif kind == QEvent.FocusOut and layer.targets.get("focus") is control:
            layer.animate("focus", control, 0., 160)
        elif isinstance(control, QAbstractButton):
            if kind == QEvent.MouseButtonPress or kind == QEvent.KeyPress and event.key() in (Qt.Key_Space, Qt.Key_Return):
                layer.animate("press", control, 1., 90)
            elif kind == QEvent.MouseButtonRelease or kind == QEvent.KeyRelease and event.key() in (Qt.Key_Space, Qt.Key_Return):
                layer.animate("press", control, 0., 110)
        return False


def install_console_effects(root, dark=True):
    root._console_effects = ConsoleEffects(root, dark)
