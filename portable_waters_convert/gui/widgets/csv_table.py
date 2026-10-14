"""Editable CSV table widget backed by pandas DataFrame."""
from __future__ import annotations

from pathlib import Path
import math

import pandas as pd

from PySide6.QtWidgets import (
    QTableWidget, QTableWidgetItem, QVBoxLayout, QHBoxLayout,
    QWidget, QPushButton, QFileDialog, QHeaderView, QLabel,
)
from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QFontMetrics

from gui.theme import TEXT_DIM
from deconvovo.paths import local_path


FIELD_HELP = {
    "name": "Calibrant display name.",
    "species": "Analyte display name.",
    "mw": "Molecular mass in Da; used for mass correction and automatic extraction width.",
    "z": "Positive integer charge state.",
    "mz": "Target mass-to-charge ratio; centre of the extraction window.",
    "ccs": "Calibrant reference collision cross-section in Å².",
    "mz_window": "Extraction half-width: 3 means target ±3. 'auto' estimates the isotope width.",
    "peak_select": "Drift peak selection: 'highest' is the tallest peak; 'last' is the latest peak.",
    "min_ms_intensity": "Minimum maximum MS intensity inside the extraction window.",
    "min_im_intensity": "Minimum drift-peak intensity after summing the m/z extraction window.",
    "peak_height_frac": "Peak-detection fraction of the drift-profile maximum: 0.1 means 10%.",
    "raw_path": "Complete Waters RAW acquisition folder. Relative paths use this CSV's folder.",
    "formula": "Optional molecular formula for automatic isotope-width estimation.",
    "gas_mw": "Gas molecular mass in Da; the default is 28.014.",
}

FIELD_WIDTHS = {
    "name": 180, "species": 240, "mw": 110, "z": 65, "mz": 115, "ccs": 115,
    "mz_window": 130, "peak_select": 130, "min_ms_intensity": 160,
    "min_im_intensity": 160, "peak_height_frac": 155, "raw_path": 420,
    "formula": 150, "gas_mw": 115,
}


class CsvTable(QWidget):
    """Editable table loaded from / saved to CSV."""

    data_changed = Signal()  # emitted when user edits or loads
    csv_path_changed = Signal(str)

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
        btn_load.setProperty("compact", True)
        btn_load.setMinimumHeight(30)
        btn_load.clicked.connect(self._load_csv)
        toolbar.addWidget(btn_load)

        btn_save = QPushButton("Save CSV")
        btn_save.setProperty("secondary", True)
        btn_save.setProperty("compact", True)
        btn_save.setMinimumHeight(30)
        btn_save.clicked.connect(self._save_csv)
        toolbar.addWidget(btn_save)

        btn_add = QPushButton("+ Row")
        btn_add.setProperty("secondary", True)
        btn_add.setProperty("compact", True)
        btn_add.setMinimumHeight(30)
        btn_add.setVisible(editable)
        btn_add.clicked.connect(self._add_row)
        toolbar.addWidget(btn_add)

        btn_del = QPushButton("- Row")
        btn_del.setProperty("secondary", True)
        btn_del.setProperty("compact", True)
        btn_del.setMinimumHeight(30)
        btn_del.setVisible(editable)
        btn_del.clicked.connect(self._del_row)
        toolbar.addWidget(btn_del)

        layout.addLayout(toolbar)

        self.table = QTableWidget()
        self.table.setAlternatingRowColors(True)
        self.table.setStyleSheet("QTableWidget, QHeaderView::section { font-size: 10pt; }")
        self.table.setWordWrap(False)
        self.table.setTextElideMode(Qt.ElideRight)
        self.table.setHorizontalScrollMode(QTableWidget.ScrollPerPixel)
        self.table.setVerticalScrollMode(QTableWidget.ScrollPerPixel)
        self.table.horizontalHeader().setStretchLastSection(False)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.Interactive)
        self.table.horizontalHeader().setMinimumSectionSize(60)
        self.table.horizontalHeader().setMinimumHeight(34)
        self.table.verticalHeader().setDefaultSectionSize(34)
        self.table.setMinimumHeight(200)
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
                item.setToolTip(txt)
                if not self._editable:
                    item.setFlags(item.flags() & ~Qt.ItemIsEditable)
                self.table.setItem(r, c, item)

        metrics = QFontMetrics(self.table.font())
        for column, key in enumerate(df.columns):
            header = self.table.horizontalHeaderItem(column)
            header.setToolTip(FIELD_HELP.get(str(key), f"CSV field: {key}") + "\nDrag the column edge to resize; hover over a cell for its full value.")
            width = FIELD_WIDTHS.get(str(key), min(260, max(130, metrics.horizontalAdvance(str(key)) + 30)))
            self.table.setColumnWidth(column, width)
        self._lbl.setText(f"{len(df)} rows · {len(df.columns)} columns" + (f"  •  {path.name}" if path else ""))
        self._lbl.setToolTip("Double-click a cell to edit. Scroll horizontally for additional fields; hover over headers for field definitions.")
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
            item = self.table.item(row, col)
            if item and item.toolTip() != item.text():
                item.setToolTip(item.text())
            self.data_changed.emit()

    def validated_dataframe(self, required: set[str], label: str = "CSV") -> pd.DataFrame:
        """Reject invalid required numbers before submitting edited CCS rows."""
        frame = self.to_dataframe()
        missing = required - set(frame.columns)
        if frame.empty or missing:
            raise ValueError(f"{label} table needs rows and columns: {', '.join(sorted(required))}.")
        for key in ("mw", "z", "mz", "ccs"):
            if key not in required:
                continue
            values = pd.to_numeric(frame[key], errors="coerce")
            valid = values.map(lambda value: math.isfinite(value) and value > 0)
            if key == "z":
                valid &= values.mod(1).eq(0)
            if not valid.all():
                row = int(valid[~valid].index[0])
                self.table.setCurrentCell(row, frame.columns.get_loc(key))
                kind = "positive integer" if key == "z" else "positive finite number"
                raise ValueError(f"{label} row {row + 1}: {key} must be a {kind}.")
        return frame

    def write_snapshot(self, destination: Path, required: set[str], label: str) -> Path:
        """Submit current cells without moving relative RAW references beside a temp file."""
        frame = self.validated_dataframe(required, label)
        if "raw_path" in frame:
            def absolute_reference(value):
                path = local_path(str(value))
                if not path.is_absolute() and self._csv_path:
                    relative = self._csv_path.parent / path
                    path = relative if relative.exists() or not path.exists() else path
                return str(path.resolve())
            frame["raw_path"] = frame["raw_path"].map(absolute_reference)
        frame.to_csv(destination, index=False)
        return destination

    def _load_csv(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Load CSV", "", "CSV Files (*.csv);;All Files (*)")
        if path:
            df = pd.read_csv(path)
            self.load_dataframe(df, Path(path))
            self.csv_path_changed.emit(str(path))
            self.data_changed.emit()

    def _save_csv(self):
        default = str(self._csv_path) if self._csv_path else ""
        path, _ = QFileDialog.getSaveFileName(
            self, "Save CSV", default, "CSV Files (*.csv)")
        if path:
            self.to_dataframe().to_csv(path, index=False)
            self._csv_path = Path(path)
            self._lbl.setText(f"{self.table.rowCount()} rows  •  {Path(path).name}")
            self.csv_path_changed.emit(str(path))

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
