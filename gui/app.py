"""DeconVoVo main application window."""
from __future__ import annotations

import sys

from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QHBoxLayout, QVBoxLayout,
    QListWidget, QStackedWidget, QSplitter, QLabel,
)
from PySide6.QtCore import Qt, QSize
from PySide6.QtGui import QFont, QIcon

from gui.theme import stylesheet, BG_MID, TEXT_BRIGHT, TEXT_DIM, ACCENT, FONT_FAMILY
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
        self.setWindowTitle("DeconVoVo — IM-MS Analysis Suite")
        self.setMinimumSize(1200, 780)
        self.resize(1400, 860)

        central = QWidget()
        self.setCentralWidget(central)
        root = QHBoxLayout(central)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)

        # ---- Sidebar ----
        sidebar_w = QWidget()
        sidebar_l = QVBoxLayout(sidebar_w)
        sidebar_l.setContentsMargins(0, 0, 0, 0)
        sidebar_l.setSpacing(0)

        # Logo / title
        title_w = QWidget()
        title_w.setStyleSheet(f"background-color: {BG_MID};")
        title_l = QVBoxLayout(title_w)
        title_l.setContentsMargins(16, 16, 16, 12)
        lbl_title = QLabel("DeconVoVo")
        lbl_title.setProperty("title", True)
        lbl_title.setStyleSheet(f"color: {ACCENT}; font-size: 16pt; font-weight: bold;")
        lbl_sub = QLabel("IM-MS Analysis Suite")
        lbl_sub.setProperty("subtitle", True)
        title_l.addWidget(lbl_title)
        title_l.addWidget(lbl_sub)
        sidebar_l.addWidget(title_w)

        self.nav = QListWidget()
        self.nav.setFixedWidth(200)
        for name, desc, _ in PANELS:
            self.nav.addItem(name)
        self.nav.setCurrentRow(0)
        self.nav.currentRowChanged.connect(self._switch_panel)
        sidebar_l.addWidget(self.nav, stretch=1)

        # Version label
        from gui import __version__
        ver = QLabel(f"  v{__version__}")
        ver.setProperty("subtitle", True)
        ver.setStyleSheet(f"background-color: {BG_MID}; padding: 8px;")
        sidebar_l.addWidget(ver)

        root.addWidget(sidebar_w)

        # ---- Right side: panels + log ----
        right = QSplitter(Qt.Vertical)

        # Panel stack
        self.stack = QStackedWidget()
        self.log_panel = LogPanel()

        self._panels = []
        for name, desc, cls in PANELS:
            panel = cls(log_panel=self.log_panel)
            self._panels.append(panel)
            self.stack.addWidget(panel)

        right.addWidget(self.stack)
        right.addWidget(self.log_panel)
        right.setStretchFactor(0, 4)
        right.setStretchFactor(1, 1)

        root.addWidget(right, stretch=1)

    def _switch_panel(self, idx: int):
        if 0 <= idx < self.stack.count():
            self.stack.setCurrentIndex(idx)
            name, desc, _ = PANELS[idx]
            self.log_panel.log(f"--- {name}: {desc} ---")

    def get_panel(self, idx: int):
        return self._panels[idx]


def main():
    app = QApplication(sys.argv)
    app.setStyle("Fusion")
    app.setStyleSheet(stylesheet())
    app.setFont(QFont(FONT_FAMILY, 10))

    win = MainWindow()
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
