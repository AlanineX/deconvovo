"""Embedded matplotlib canvas for Qt panels."""
from __future__ import annotations

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure

from PySide6.QtWidgets import QVBoxLayout, QWidget

from gui.theme import mpl_rc

# Module-level flag — toggled by app._toggle_theme()
_dark_mode = True


def set_plot_dark(dark: bool):
    global _dark_mode
    _dark_mode = dark


class PlotCanvas(QWidget):
    """Embeddable matplotlib figure. Follows app dark/light mode."""

    def __init__(self, parent=None, figsize=(6, 3.5), nrows=1, ncols=1):
        super().__init__(parent)
        self.fig = Figure(figsize=figsize, constrained_layout=True)
        self.axes = self.fig.subplots(nrows, ncols, squeeze=False)
        self._apply_theme()
        self.canvas = FigureCanvasQTAgg(self.fig)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.canvas)

    def _apply_theme(self):
        rc = mpl_rc(dark=_dark_mode)
        self.fig.set_facecolor(rc["figure.facecolor"])
        for row in self.axes:
            for ax in row:
                ax.set_facecolor(rc["axes.facecolor"])
                ax.tick_params(colors=rc["xtick.color"])
                ax.xaxis.label.set_color(rc["axes.labelcolor"])
                ax.yaxis.label.set_color(rc["axes.labelcolor"])
                ax.title.set_color(rc["axes.titlecolor"])
                for spine in ax.spines.values():
                    spine.set_color(rc["axes.edgecolor"])

    def ax(self, row=0, col=0):
        return self.axes[row][col]

    def clear(self):
        for row in self.axes:
            for a in row:
                a.clear()
        self._apply_theme()
        self.canvas.draw_idle()

    def refresh(self):
        self._apply_theme()
        self.canvas.draw_idle()
