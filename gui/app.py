"""DeconVoVo main application window."""
from __future__ import annotations

import os
import sys

# Force non-interactive backend BEFORE any matplotlib import
# This MUST be before any deconvovo import that touches matplotlib
os.environ["MPLBACKEND"] = "Agg"
import matplotlib
matplotlib.use("Agg")

from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QHBoxLayout, QVBoxLayout,
    QListWidget, QStackedWidget, QSplitter, QLabel, QPushButton,
)
from PySide6.QtCore import Qt
from PySide6.QtGui import QFont

from gui.theme import stylesheet, ACCENT, FONT_FAMILY
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
        self.setMinimumSize(900, 600)
        from PySide6.QtGui import QGuiApplication
        screen = QGuiApplication.primaryScreen().availableGeometry()
        w = min(int(screen.width() * 0.85), 1400)
        h = min(int(screen.height() * 0.85), 860)
        self.resize(w, h)

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

        # Title
        title_w = QWidget()
        title_l = QVBoxLayout(title_w)
        title_l.setContentsMargins(16, 16, 16, 12)
        lbl_title = QLabel("DeconVoVo")
        lbl_title.setStyleSheet(f"color: {ACCENT}; font-size: 16pt; font-weight: bold; background: transparent;")
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

        # Bottom: theme toggle + version
        bot_w = QWidget()
        bot_l = QHBoxLayout(bot_w)
        bot_l.setContentsMargins(8, 4, 8, 8)
        self.btn_theme = QPushButton("Light Mode")
        self.btn_theme.setProperty("secondary", True)
        self.btn_theme.setFixedHeight(26)
        self.btn_theme.clicked.connect(self._toggle_theme)
        bot_l.addWidget(self.btn_theme)
        bot_l.addStretch()
        from gui import __version__
        ver = QLabel(f"v{__version__}")
        ver.setProperty("subtitle", True)
        bot_l.addWidget(ver)
        sidebar_l.addWidget(bot_w)

        root.addWidget(sidebar_w)

        # ---- Right side: panels + log ----
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
        right.setStretchFactor(0, 4)
        right.setStretchFactor(1, 1)

        root.addWidget(right, stretch=1)

        # Apply Windows transparency after window is shown
        if sys.platform == "win32":
            from PySide6.QtCore import QTimer
            QTimer.singleShot(0, self._apply_win11_backdrop)

    def _apply_win11_backdrop(self):
        """Apply Windows 11 Mica/Acrylic backdrop. No-op if unsupported."""
        try:
            import ctypes
            hwnd = int(self.winId())
            dwm = ctypes.windll.dwmapi

            # DWMWA_USE_IMMERSIVE_DARK_MODE = 20
            val = ctypes.c_int(1 if self._dark else 0)
            dwm.DwmSetWindowAttribute(hwnd, 20, ctypes.byref(val), 4)

            # DWMWA_SYSTEMBACKDROP_TYPE = 38 (Win11 22H2+)
            # 2 = Mica, 3 = Acrylic, 4 = Mica Alt
            backdrop = ctypes.c_int(3)  # Acrylic
            dwm.DwmSetWindowAttribute(hwnd, 38, ctypes.byref(backdrop), 4)
        except Exception:
            pass

    def _switch_panel(self, idx: int):
        if 0 <= idx < self.stack.count():
            self.stack.setCurrentIndex(idx)
            name, desc, _ = PANELS[idx]
            self.log_panel.log(f"--- {name}: {desc} ---")

    def _toggle_theme(self):
        self._dark = not self._dark
        QApplication.instance().setStyleSheet(stylesheet(dark=self._dark))
        self.btn_theme.setText("Dark Mode" if not self._dark else "Light Mode")
        if sys.platform == "win32":
            self._apply_win11_backdrop()

    def get_panel(self, idx: int):
        return self._panels[idx]


def main():
    # Signal pipeline code to avoid multiprocessing (prevents extra windows)
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
