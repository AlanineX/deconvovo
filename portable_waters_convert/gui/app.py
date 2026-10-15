"""DeconVoVo desktop shell for the Instrument Console workspaces."""
from __future__ import annotations

import sys

from PySide6.QtCore import Qt
from PySide6.QtGui import QFont, QGuiApplication
from PySide6.QtWidgets import (
    QApplication, QGroupBox, QHBoxLayout, QLabel, QListView, QListWidget,
    QMainWindow, QPushButton, QSizePolicy, QSplitter, QStackedWidget,
    QVBoxLayout, QWidget,
)

from gui.theme import ConsoleBackdrop, FONT_FAMILY, decorate_section, refresh_console_style, stylesheet
from gui.widgets.log_panel import LogPanel
from gui.panels.waters_workspace import WatersWorkspacePanel

PANELS = [
    ("Waters IM-MS", "Process runs, inspect viewers, and review CCS", WatersWorkspacePanel),
]


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self._dark = True
        QApplication.setAttribute(Qt.AA_DontUseNativeDialogs, True)
        self.setWindowTitle("DeconVoVo · Instrument Console")
        self.setMinimumSize(1100, 700)
        screen = QGuiApplication.primaryScreen()
        if screen:
            bounds = screen.availableGeometry()
            self.resize(min(1440, int(bounds.width() * 0.92)), min(960, int(bounds.height() * 0.92)))
        else:
            self.resize(1440, 900)

        central = ConsoleBackdrop()
        self.setCentralWidget(central)
        layout = QVBoxLayout(central)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)

        header = QWidget()
        header.setObjectName("appHeader")
        header.setFixedHeight(58)
        row = QHBoxLayout(header)
        row.setContentsMargins(20, 8, 20, 8)
        row.setSpacing(14)
        brand = QLabel("DeconVoVo")
        brand.setObjectName("appBrand")
        row.addWidget(brand)
        self.nav = QListWidget()
        self.nav.setObjectName("modeNav")
        self.nav.setFlow(QListView.LeftToRight)
        self.nav.setWrapping(False)
        self.nav.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.nav.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.nav.setFixedHeight(40)
        self.nav.setFixedWidth(350 if len(PANELS) > 1 else 160)
        for name, description, _ in PANELS:
            self.nav.addItem(name)
            self.nav.item(self.nav.count() - 1).setToolTip(description)
        row.addWidget(self.nav)
        self.nav.setVisible(len(PANELS) > 1)
        row.addStretch()
        self.btn_theme = QPushButton("Light mode")
        self.btn_theme.setProperty("secondary", True)
        self.btn_theme.setProperty("compact", True)
        self.btn_theme.clicked.connect(self._toggle_theme)
        row.addWidget(self.btn_theme)
        layout.addWidget(header)

        self.stack = QStackedWidget()
        self.log_panel = LogPanel()
        self.log_panel.setMinimumHeight(110)
        self.log_panel.hide()
        self._panels = []
        for _name, _description, cls in PANELS:
            panel = cls(log_panel=self.log_panel)
            if isinstance(panel, WatersWorkspacePanel):
                panel.job_status.hide()  # Status is already shown in the footer and activity log.
            self._panels.append(panel)
            self.stack.addWidget(panel)
        self.body = QSplitter(Qt.Vertical)
        self.body.setChildrenCollapsible(False)
        self.body.addWidget(self.stack)
        self.body.addWidget(self.log_panel)
        self.body.setStretchFactor(0, 1)
        self.body.setStretchFactor(1, 0)
        layout.addWidget(self.body, 1)

        footer = QWidget()
        footer.setObjectName("appFooter")
        footer.setFixedHeight(38)
        bottom = QHBoxLayout(footer)
        bottom.setContentsMargins(20, 4, 20, 4)
        self.status_text = QLabel("Ready · choose an input folder or an example.")
        self.status_text.setObjectName("appStatus")
        self.status_text.setTextFormat(Qt.PlainText)
        self.status_text.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Preferred)
        bottom.addWidget(self.status_text, 1)
        self.btn_log = QPushButton("Activity log")
        self.btn_log.setProperty("secondary", True)
        self.btn_log.setProperty("compact", True)
        self.btn_log.setCheckable(True)
        self.btn_log.toggled.connect(self._show_log)
        bottom.addWidget(self.btn_log)
        from gui import __version__
        version = QLabel(f"v{__version__}")
        version.setProperty("eyebrow", True)
        bottom.addWidget(version)
        layout.addWidget(footer)
        self.log_panel._sig.message.connect(self._status_message)
        self.nav.currentRowChanged.connect(self._switch_panel)
        self.nav.setCurrentRow(next(i for i, entry in enumerate(PANELS) if entry[2] is WatersWorkspacePanel))
        sections = {"Data": "data", "Acquisitions": "acquisitions", "Display": "display",
                    "Axis ranges": "ranges", "CCS profile": "ccs", "Optional CCS": "ccs"}
        for group in self.findChildren(QGroupBox):
            if group.title() in sections:
                decorate_section(group, sections[group.title()])
        refresh_console_style(self, self._dark)
        from gui.effects import install_console_effects
        install_console_effects(self, self._dark)

    def _show_log(self, visible):
        self.log_panel.setVisible(visible)
        self.btn_log.setText("Hide activity log" if visible else "Activity log")
        if visible:
            self.body.setSizes([max(400, self.body.height() - 160), 160])

    def _status_message(self, message, level):
        if not message:
            return
        self.status_text.setText(message.splitlines()[0])
        self.status_text.setToolTip(message)
        self.status_text.setProperty("statusLevel", level)
        self.status_text.style().unpolish(self.status_text)
        self.status_text.style().polish(self.status_text)
        if not self._dark:
            refresh_console_style(self, False)
        if level == "error":
            self.btn_log.setChecked(True)
        elif hasattr(self, "_console_effects") and (level == "success" or message.startswith("CCS ready")):
            self._console_effects.pulse(self.status_text, warning=level == "warning")

    def _switch_panel(self, index):
        if 0 <= index < self.stack.count():
            self.stack.setCurrentIndex(index)
            self.log_panel.log(f"{PANELS[index][0]} workspace", "info")

    def _toggle_theme(self):
        self._dark = not self._dark
        QApplication.instance().setStyleSheet(stylesheet(dark=self._dark))
        self.btn_theme.setText("Light mode" if self._dark else "Dark mode")
        refresh_console_style(self, self._dark)
        self.log_panel.refresh_theme(self._dark)
        from gui.widgets.plot_canvas import PlotCanvas, set_plot_dark
        set_plot_dark(self._dark)
        for canvas in self.findChildren(PlotCanvas):
            canvas.refresh()
        for panel in self._panels:
            if isinstance(panel, WatersWorkspacePanel):
                panel.viewer_widget.refresh_theme()
                panel.ccs_workbench.refresh_theme(self._dark)

    def get_panel(self, index):
        return self._panels[index]

    def closeEvent(self, event):
        if any(isinstance(panel, WatersWorkspacePanel) and panel.is_busy for panel in self._panels):
            self.log_panel.log("A Waters job is running. Wait for it to finish before closing.", "warning")
            event.ignore()
            return
        for panel in self._panels:
            panel.close()
        self._console_effects.shutdown()
        super().closeEvent(event)


def main():
    from gui.widgets.worker import set_gui_mode
    set_gui_mode()
    app = QApplication(sys.argv)
    app.setOrganizationName("DeconVoVo")
    app.setApplicationName("DeconVoVo Waters")
    app.setStyle("Fusion")
    app.setStyleSheet(stylesheet(dark=True))
    app.setFont(QFont(FONT_FAMILY, 10))
    window = MainWindow()
    window.show()
    raise SystemExit(app.exec())


if __name__ == "__main__":
    main()
