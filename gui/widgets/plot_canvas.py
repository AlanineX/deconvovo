"""Embedded matplotlib canvas for Qt panels."""
from __future__ import annotations

import matplotlib
matplotlib.use("QtAgg")

import matplotlib.pyplot as plt
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure

from PySide6.QtWidgets import QVBoxLayout, QWidget

from gui.theme import mpl_rc


class PlotCanvas(QWidget):
    """Embeddable matplotlib figure with dark theme."""

    def __init__(self, parent=None, figsize=(6, 3.5), nrows=1, ncols=1):
        super().__init__(parent)
        plt.rcParams.update(mpl_rc(dark=True))

        self.fig = Figure(figsize=figsize, constrained_layout=True)
        self.axes = self.fig.subplots(nrows, ncols, squeeze=False)
        self.canvas = FigureCanvasQTAgg(self.fig)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.canvas)

    def ax(self, row=0, col=0):
        return self.axes[row][col]

    def clear(self):
        for row in self.axes:
            for a in row:
                a.clear()
        self.canvas.draw_idle()

    def refresh(self):
        self.canvas.draw_idle()
