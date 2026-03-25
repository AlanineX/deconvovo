"""DeconVoVo main application window."""
from __future__ import annotations

import os
import sys

from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QHBoxLayout, QVBoxLayout,
    QListWidget, QStackedWidget, QSplitter, QLabel, QPushButton, QFrame,
)
from PySide6.QtCore import Qt
from PySide6.QtGui import QFont

from gui.theme import stylesheet, ACCENT, ACCENT_DIM, FONT_FAMILY
from gui.widgets.log_panel import LogPanel
from gui.panels.convert import ConvertPanel
from gui.panels.html_viewer import HtmlViewerPanel
from gui.panels.ccs_calibrate import CcsCalibrationPanel
from gui.panels.ccs_analyte import CcsAnalytePanel
from gui.panels.full_pipeline import FullPipelinePanel


PANELS = [
    ("Convert .raw",    "Convert Waters .raw directories to text files",   ConvertPanel),
    ("2D IM-MS Viewer", "Generate interactive HTML IM-MS heatmaps",        HtmlViewerPanel),
    ("CCS Calibration", "Build TW-IMS CCS calibration curve",             CcsCalibrationPanel),
    ("Analyte CCS",     "Apply calibration to compute analyte CCS",        CcsAnalytePanel),
    ("Full Pipeline",   "Run the complete workflow end-to-end",            FullPipelinePanel),
]


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self._dark = True
        self.setWindowTitle("DeconVoVo")
        self.setMinimumSize(960, 640)
        from PySide6.QtGui import QGuiApplication
        screen = QGuiApplication.primaryScreen().availableGeometry()
        w = min(int(screen.width() * 0.85), 1440)
        h = min(int(screen.height() * 0.85), 900)
        self.resize(w, h)

        central = QWidget()
        self.setCentralWidget(central)
        root = QHBoxLayout(central)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)

        # ---- Sidebar ----
        sidebar = QWidget()
        sidebar.setFixedWidth(220)
        sb = QVBoxLayout(sidebar)
        sb.setContentsMargins(0, 0, 0, 0)
        sb.setSpacing(0)

        # Brand
        brand = QWidget()
        bl = QVBoxLayout(brand)
        bl.setContentsMargins(20, 24, 20, 16)
        bl.setSpacing(2)
        lbl_title = QLabel("DeconVoVo")
        lbl_title.setStyleSheet(f"color: {ACCENT}; font-size: 18pt; font-weight: 700; background: transparent;")
        lbl_sub = QLabel("IM-MS Analysis Suite")
        lbl_sub.setProperty("subtitle", True)
        lbl_sub.setStyleSheet("padding-bottom: 0;")
        bl.addWidget(lbl_title)
        bl.addWidget(lbl_sub)
        sb.addWidget(brand)

        # Separator
        sep = QFrame()
        sep.setFrameShape(QFrame.HLine)
        sep.setStyleSheet(f"color: {ACCENT_DIM}; margin: 0 16px;")
        sb.addWidget(sep)

        # Navigation
        self.nav = QListWidget()
        self.nav.setObjectName("nav")
        for name, _, _ in PANELS:
            self.nav.addItem(name)
        self.nav.setCurrentRow(0)
        self.nav.currentRowChanged.connect(self._switch_panel)
        sb.addWidget(self.nav, stretch=1)

        # Bottom: theme + version
        bot = QWidget()
        bot_l = QHBoxLayout(bot)
        bot_l.setContentsMargins(16, 8, 16, 12)
        self.btn_theme = QPushButton("Light Mode")
        self.btn_theme.setProperty("secondary", True)
        self.btn_theme.setFixedHeight(28)
        self.btn_theme.clicked.connect(self._toggle_theme)
        bot_l.addWidget(self.btn_theme)
        bot_l.addStretch()
        from gui import __version__
        ver = QLabel(f"v{__version__}")
        ver.setProperty("subtitle", True)
        ver.setStyleSheet("padding-bottom: 0;")
        bot_l.addWidget(ver)
        sb.addWidget(bot)

        root.addWidget(sidebar)

        # ---- Content ----
        right = QSplitter(Qt.Vertical)

        self.stack = QStackedWidget()
        self.log_panel = LogPanel()

        self._panels = []
        for name, desc, cls in PANELS:
            panel = cls(log_panel=self.log_panel)
            self._panels.append(panel)
            self.stack.addWidget(panel)

        right.addWidget(self.stack)
        right.addWidget(self.log_panel)
        right.setStretchFactor(0, 5)
        right.setStretchFactor(1, 1)
        right.setHandleWidth(1)

        root.addWidget(right, stretch=1)

    def _switch_panel(self, idx: int):
        if 0 <= idx < self.stack.count():
            self.stack.setCurrentIndex(idx)
            name, desc, _ = PANELS[idx]
            self.log_panel.log(f"--- {name}: {desc} ---")

    def _toggle_theme(self):
        self._dark = not self._dark
        QApplication.instance().setStyleSheet(stylesheet(dark=self._dark))
        self.btn_theme.setText("Dark Mode" if not self._dark else "Light Mode")
        # Update all embedded plot canvases
        from gui.widgets.plot_canvas import PlotCanvas, set_plot_dark
        set_plot_dark(self._dark)
        for canvas in self.findChildren(PlotCanvas):
            canvas.refresh()

    def get_panel(self, idx: int):
        return self._panels[idx]


def main():
    # Signal pipeline code to avoid multiprocessing
    from gui.widgets.worker import set_gui_mode
    set_gui_mode()

    app = QApplication(sys.argv)
    app.setStyle("Fusion")
    app.setStyleSheet(stylesheet(dark=True))
    app.setFont(QFont(FONT_FAMILY, 10))

    win = MainWindow()
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
