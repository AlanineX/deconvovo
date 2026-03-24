"""Editable CSV table widget backed by pandas DataFrame."""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from PySide6.QtWidgets import (
    QTableWidget, QTableWidgetItem, QVBoxLayout, QHBoxLayout,
    QWidget, QPushButton, QFileDialog, QHeaderView, QLabel,
)
from PySide6.QtCore import Qt, Signal

from gui.theme import TEXT_DIM


class CsvTable(QWidget):
    """Editable table loaded from / saved to CSV."""

    data_changed = Signal()  # emitted when user edits or loads

    def __init__(self, parent=None, editable: bool = True):
        super().__init__(parent)
        self._editable = editable
        self._df = pd.DataFrame()

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(4)

        # Toolbar
        toolbar = QHBoxLayout()
        self._lbl = QLabel("")
        self._lbl.setStyleSheet(f"color: {TEXT_DIM}; font-size: 9pt;")
        toolbar.addWidget(self._lbl)
        toolbar.addStretch()

        btn_load = QPushButton("Load CSV")
        btn_load.setProperty("secondary", True)
        btn_load.setFixedHeight(26)
        btn_load.clicked.connect(self._load_csv)
        toolbar.addWidget(btn_load)

        btn_save = QPushButton("Save CSV")
        btn_save.setProperty("secondary", True)
        btn_save.setFixedHeight(26)
        btn_save.clicked.connect(self._save_csv)
        toolbar.addWidget(btn_save)

        btn_add = QPushButton("+ Row")
        btn_add.setProperty("secondary", True)
        btn_add.setFixedHeight(26)
        btn_add.clicked.connect(self._add_row)
        toolbar.addWidget(btn_add)

        btn_del = QPushButton("- Row")
        btn_del.setProperty("secondary", True)
        btn_del.setFixedHeight(26)
        btn_del.clicked.connect(self._del_row)
        toolbar.addWidget(btn_del)

        layout.addLayout(toolbar)

        self.table = QTableWidget()
        self.table.setAlternatingRowColors(True)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.verticalHeader().setDefaultSectionSize(26)
        if editable:
            self.table.cellChanged.connect(self._on_cell_changed)
        layout.addWidget(self.table)

        self._csv_path: Path | None = None
        self._loading = False

    def load_dataframe(self, df: pd.DataFrame, path: Path | None = None):
        self._loading = True
        self._df = df.copy()
        self._csv_path = path

        self.table.setRowCount(len(df))
        self.table.setColumnCount(len(df.columns))
        self.table.setHorizontalHeaderLabels(list(df.columns))

        for r in range(len(df)):
            for c in range(len(df.columns)):
                val = df.iloc[r, c]
                txt = "" if pd.isna(val) else str(val)
                item = QTableWidgetItem(txt)
                if not self._editable:
                    item.setFlags(item.flags() & ~Qt.ItemIsEditable)
                self.table.setItem(r, c, item)

        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self._lbl.setText(f"{len(df)} rows" + (f"  •  {path.name}" if path else ""))
        self._loading = False

    def to_dataframe(self) -> pd.DataFrame:
        rows = []
        cols = [self.table.horizontalHeaderItem(c).text()
                for c in range(self.table.columnCount())]
        for r in range(self.table.rowCount()):
            row = {}
            for c in range(self.table.columnCount()):
                item = self.table.item(r, c)
                row[cols[c]] = item.text() if item else ""
            rows.append(row)
        return pd.DataFrame(rows, columns=cols)

    def _on_cell_changed(self, row, col):
        if not self._loading:
            self.data_changed.emit()

    def _load_csv(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Load CSV", "", "CSV Files (*.csv);;All Files (*)")
        if path:
            df = pd.read_csv(path)
            self.load_dataframe(df, Path(path))
            self.data_changed.emit()

    def _save_csv(self):
        default = str(self._csv_path) if self._csv_path else ""
        path, _ = QFileDialog.getSaveFileName(
            self, "Save CSV", default, "CSV Files (*.csv)")
        if path:
            self.to_dataframe().to_csv(path, index=False)
            self._csv_path = Path(path)
            self._lbl.setText(f"{self.table.rowCount()} rows  •  {Path(path).name}")

    def _add_row(self):
        r = self.table.rowCount()
        self.table.insertRow(r)
        for c in range(self.table.columnCount()):
            self.table.setItem(r, c, QTableWidgetItem(""))
        self.data_changed.emit()

    def _del_row(self):
        rows = sorted(set(idx.row() for idx in self.table.selectedIndexes()), reverse=True)
        for r in rows:
            self.table.removeRow(r)
        if not rows and self.table.rowCount() > 0:
            self.table.removeRow(self.table.rowCount() - 1)
        self.data_changed.emit()
