"""Embedded matplotlib canvas for Qt panels."""
from __future__ import annotations

import matplotlib
matplotlib.use("QtAgg")

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure

from PySide6.QtWidgets import QVBoxLayout, QWidget

from gui.theme import mpl_rc


class PlotCanvas(QWidget):
    """Embeddable matplotlib figure with dark theme.

    Uses per-figure styling — does NOT modify global rcParams,
    so PNG output from ccs_plot.py keeps its own (light) theme.
    """
    def __init__(self, parent=None, figsize=(6, 3.5), nrows=1, ncols=1):
        super().__init__(parent)
        rc = mpl_rc(dark=True)
        self.fig = Figure(figsize=figsize, constrained_layout=True,
                          facecolor=rc["figure.facecolor"])
        self.axes = self.fig.subplots(nrows, ncols, squeeze=False)
        # Apply dark theme per-axes, not globally
        for row in self.axes:
            for ax in row:
                ax.set_facecolor(rc["axes.facecolor"])
                ax.tick_params(colors=rc["xtick.color"])
                ax.xaxis.label.set_color(rc["axes.labelcolor"])
                ax.yaxis.label.set_color(rc["axes.labelcolor"])
                ax.title.set_color(rc["axes.titlecolor"])
                for spine in ax.spines.values():
                    spine.set_color(rc["axes.edgecolor"])

        self.canvas = FigureCanvasQTAgg(self.fig)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.canvas)

    def ax(self, row=0, col=0):
        return self.axes[row][col]

    def clear(self):
        rc = mpl_rc(dark=True)
        for row in self.axes:
            for a in row:
                a.clear()
                a.set_facecolor(rc["axes.facecolor"])
        self.canvas.draw_idle()

    def refresh(self):
        self.canvas.draw_idle()
