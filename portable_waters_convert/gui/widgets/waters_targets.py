"""Compact editable m/z targets for the current Waters acquisition."""
from __future__ import annotations

import math
from pathlib import Path
import tempfile

import pandas as pd
from PySide6.QtCore import Qt, Signal, QSize, QSignalBlocker
from PySide6.QtGui import QFont, QPalette, QTextDocument
from PySide6.QtWidgets import (
    QAbstractItemView, QFileDialog, QHeaderView, QHBoxLayout, QLabel,
    QPushButton, QStyle, QStyleOptionHeader, QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget,
)

from deconvovo.ccs_peak_pick import compute_auto_mz_window
from deconvovo.paths import local_path


FIELDS = ("species", "mw", "z", "mz", "mz_window", "mz_min", "mz_max", "raw_path",
          "drift_mode", "drift_min", "drift_max")
VISIBLE_FIELDS = ("species", "mz_min", "mz_max", "drift_mode", "drift_min", "drift_max", "mw", "z")
RESULT_FIELDS = ("peak_ccs", "peak_drift_time_ms", "status")
INPUT_COLUMNS = {field: column for column, field in enumerate(VISIBLE_FIELDS)}
RESULT_COLUMNS = {field: len(VISIBLE_FIELDS) + column for column, field in enumerate(RESULT_FIELDS)}
MATCH_FIELDS = ("species", "mz", "mw", "z", "mz_window")
HEADERS = ("Name", "m/z min", "m/z max", "Auto", "Min td", "Max td",
           "Mass (Da)", "Charge", "CCS (Å²)", "Drift peak (ms)", "Status")


class _TargetsHeader(QHeaderView):
    """Use Qt's native header chrome and rich text for the two drift labels."""

    def paintSection(self, painter, rect, logical_index):
        labels = {INPUT_COLUMNS["drift_min"]: "Min", INPUT_COLUMNS["drift_max"]: "Max"}
        if logical_index not in labels:
            return super().paintSection(painter, rect, logical_index)
        option = QStyleOptionHeader()
        self.initStyleOption(option)
        self.initStyleOptionForIndex(option, logical_index)
        option.rect = rect
        painter.save()
        self.style().drawControl(QStyle.CE_HeaderSection, option, painter, self)
        document = QTextDocument()
        document.setDocumentMargin(0)
        font = QFont(self.font())
        font.setPointSizeF(10)
        font.setWeight(QFont.DemiBold)
        document.setDefaultFont(font)
        color = option.palette.color(QPalette.ButtonText).name()
        document.setHtml(f'<span style="color: {color}">{labels[logical_index]} t<sub>d</sub></span>')
        document.adjustSize()
        size = document.size()
        painter.setClipRect(rect)
        painter.translate(rect.x() + (rect.width() - size.width()) / 2,
                          rect.y() + (rect.height() - size.height()) / 2)
        document.drawContents(painter)
        painter.restore()


def _number(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def _bounds(lower, upper):
    lower, upper = _number(lower), _number(upper)
    return (lower, upper) if math.isfinite(lower) and math.isfinite(upper) and lower < upper else None


def _text(value):
    return "" if pd.isna(value) else str(value)


def _missing(value):
    return pd.isna(value) or isinstance(value, str) and not value.strip()


def _matches(row, record):
    if not set(MATCH_FIELDS) <= record.keys() or not set(RESULT_FIELDS) & record.keys():
        return False
    for field in MATCH_FIELDS:
        left, right = row[field], record[field]
        if _missing(left) or _missing(right):
            if not (_missing(left) and _missing(right)):
                return False
        elif field == "species":
            if str(left) != str(right):
                return False
        else:
            left, right = _number(left), _number(right)
            if not (math.isfinite(left) and math.isfinite(right) and math.isclose(left, right, rel_tol=1e-9, abs_tol=1e-9)):
                return False
    return True


def _flag(value):
    if _missing(value):
        return "Unknown"
    return "Yes" if str(value).strip().casefold() in {"true", "1", "yes"} else "No"


class WatersTargetsWidget(QWidget):
    """Editable target regions with automatic or explicit drift bounds."""

    changed = Signal(object)
    selected = Signal(object)
    error = Signal(str)
    full_view_requested = Signal()
    preferred_height = 120

    def __init__(self, parent=None):
        super().__init__(parent)
        self._frame = pd.DataFrame(columns=FIELDS)
        self._results: list[dict] = []
        self._raw_path: Path | None = None
        self._default_range: tuple[float, float] | None = None
        self._loading = False
        self.setMinimumHeight(60)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(3)
        toolbar = QHBoxLayout()
        toolbar.setSpacing(4)
        self.label = QLabel("Targets")
        self.label.setObjectName("panelCaption")
        self.label.setStyleSheet("font-size: 10pt;")
        toolbar.addWidget(self.label)
        toolbar.addStretch()
        for name, callback in (("Full view", self._full_view), ("Add", self._add),
                               ("Delete", self._delete), ("Import CSV", self.import_csv),
                               ("Export CSV", self.export_csv)):
            button = QPushButton(name)
            button.setProperty("secondary", True)
            button.setProperty("compact", True)
            button.setMinimumHeight(24)
            button.setToolTip({
                "Full view": "Deselect the target and restore the full acquisition range.",
                "Add": "Add a target from the visible m/z range.",
                "Delete": "Delete the selected target.",
                "Import CSV": "Import m/z bounds or legacy CCS targets for this acquisition.",
                "Export CSV": "Save target inputs, drift ranges and any computed CCS results.",
            }[name])
            button.clicked.connect(lambda checked=False, callback=callback: callback())
            toolbar.addWidget(button)
        layout.addLayout(toolbar)

        self.table = QTableWidget(0, len(HEADERS))
        self.table.setHorizontalHeader(_TargetsHeader(Qt.Horizontal, self.table))
        self.table.setHorizontalHeaderLabels(HEADERS)
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SingleSelection)
        self.table.setAlternatingRowColors(True)
        self.table.setWordWrap(False)
        self.table.setStyleSheet("QTableWidget, QHeaderView::section { font-size: 10pt; }")
        self.table.setHorizontalScrollMode(QAbstractItemView.ScrollPerPixel)
        self.table.setVerticalScrollMode(QAbstractItemView.ScrollPerPixel)
        self.table.setMinimumHeight(30)
        self.table.verticalHeader().setVisible(False)
        self.table.verticalHeader().setDefaultSectionSize(26)
        header = self.table.horizontalHeader()
        header.setMinimumHeight(25)
        header.setMinimumSectionSize(55)
        header.setSectionResizeMode(QHeaderView.Interactive)
        for column, width in enumerate((120, 90, 90, 55, 80, 80, 100, 60, 90, 124, 140)):
            self.table.setColumnWidth(column, width)
        for column, help_text in enumerate((
            "Target name. Select a row to focus the spectrum on its m/z and drift ranges.",
            "Required exact lower m/z bound; enter any decimal precision.",
            "Required exact upper m/z bound; must be greater than the lower bound.",
            "Automatic drift range. Checked: include every positive raw drift sample, with equal 20% padding. Uncheck to enter a manual range.",
            "Lower drift bound in milliseconds. Automatically detected when Auto drift is checked.",
            "Upper drift bound in milliseconds. Automatically detected when Auto drift is checked.",
            "Optional molecular mass in Da; required only for CCS calculation.",
            "Optional charge state; required only for CCS calculation.",
            "Read-only CCS peak in Å². Hover for calibration and extrapolation flags.",
            "Read-only drift-time peak in milliseconds. Hover for the full result.",
            "Read-only calculation status. Editing inputs makes an existing result stale.",
        )):
            self.table.horizontalHeaderItem(column).setToolTip(help_text)
        for field, description in (("drift_min", "Minimum drift time"), ("drift_max", "Maximum drift time")):
            self.table.horizontalHeaderItem(INPUT_COLUMNS[field]).setData(Qt.AccessibleTextRole, description)
        self.set_ccs_fields_visible(False)
        self.table.itemChanged.connect(self._edited)
        self.table.itemSelectionChanged.connect(self._selection_changed)
        layout.addWidget(self.table)

    def sizeHint(self):
        return QSize(520, self.preferred_height)

    def selected_range(self):
        rows = self.table.selectionModel().selectedRows()
        if not rows:
            return None
        row = self._frame.iloc[rows[0].row()]
        return _bounds(row["mz_min"], row["mz_max"])

    def set_ccs_fields_visible(self, visible):
        """Show optional CCS inputs without changing target data or results."""
        order = ((0, 1, 2, 8, 3, 4, 5, 9, 10, 6, 7) if visible else tuple(range(len(HEADERS))))
        header = self.table.horizontalHeader()
        for position, column in enumerate(order):
            header.moveSection(header.visualIndex(column), position)
        for field in ("mw", "z"):
            self.table.setColumnHidden(INPUT_COLUMNS[field], not bool(visible))

    @property
    def active_range(self):
        return self.selected_range()

    @property
    def active_name(self):
        rows = self.table.selectionModel().selectedRows()
        return _text(self._frame.iloc[rows[0].row()]["species"]) if rows else "Full view"

    @property
    def active_drift_mode(self):
        rows = self.table.selectionModel().selectedRows()
        return self._frame.iloc[rows[0].row()]["drift_mode"] if rows else None

    @property
    def active_drift_range(self):
        rows = self.table.selectionModel().selectedRows()
        if not rows:
            return None
        row = self._frame.iloc[rows[0].row()]
        return _bounds(row["drift_min"], row["drift_max"])

    def set_detected_drift_range(self, bounds, unit="ms"):
        """Display detected limits without triggering an input-edit feedback loop."""
        rows = self.table.selectionModel().selectedRows()
        self.set_detected_drift_ranges({rows[0].row(): bounds} if rows else {}, unit)

    def automatic_regions(self):
        """Return detached, valid m/z regions needing automatic drift limits."""
        return tuple((index, bounds) for index, row in self._frame.iterrows()
                     if row["drift_mode"] == "auto" and (bounds := _bounds(row["mz_min"], row["mz_max"])) is not None)

    def set_detected_drift_ranges(self, ranges, unit="ms"):
        """Update all computed limits silently; manual rows are never overwritten."""
        for field, name in (("drift_min", "Min"), ("drift_max", "Max")):
            header = self.table.horizontalHeaderItem(INPUT_COLUMNS[field])
            header.setText(f"{name} td")
            header.setToolTip(f"{'Lower' if field == 'drift_min' else 'Upper'} drift bound in "
                             f"{'milliseconds' if unit == 'ms' else 'raw drift bins'}. "
                             "Automatically detected when Auto is checked; editable in manual mode.")
        for row, bounds in ranges.items():
            if row < len(self._frame) and self._frame.iloc[row]["drift_mode"] == "auto":
                self._frame.at[row, "drift_min"], self._frame.at[row, "drift_max"] = bounds
                self._refresh_drift_cells(row)

    def get_targets(self) -> pd.DataFrame:
        """Return a detached snapshot, including unedited imported calculation fields."""
        return self._frame.copy(deep=True)

    def set_context(self, raw_path: Path | None, default_range: tuple | None):
        self._default_range = _bounds(*default_range) if default_range is not None else None
        path = Path(raw_path).resolve() if raw_path is not None else None
        if path != self._raw_path:
            self._raw_path = path
            self._frame["raw_path"] = str(path) if path is not None else ""
            for row in range(len(self._frame)):
                self._stale_result(row)
            if not self._frame.empty:
                self.changed.emit(self.get_targets())

    def set_targets(self, frame: pd.DataFrame):
        if not frame.columns.is_unique:
            raise ValueError("Target CSV contains duplicate columns.")
        frame = frame.copy(deep=True).reset_index(drop=True)
        if "species" not in frame:
            frame["species"] = frame["name"] if "name" in frame else [f"Target {i + 1}" for i in range(len(frame))]
        explicit_bounds = {"mz_min", "mz_max"} <= set(frame.columns)
        for field in FIELDS:
            if field not in frame:
                frame[field] = "auto" if field == "drift_mode" else math.nan if field not in {"species", "raw_path"} else ""
        frame["drift_mode"] = frame["drift_mode"].map(
            lambda value: "manual" if str(value).strip().casefold() in {"manual", "false", "0", "no", "off"} else "auto")
        for field in ("mw", "z", "mz", "mz_min", "mz_max", "drift_min", "drift_max"):
            frame[field] = frame[field].map(_number)
        for field in FIELDS:
            frame[field] = frame[field].astype(object)
        if not explicit_bounds:
            for index, row in frame.iterrows():
                centre, width = _number(row["mz"]), _number(row["mz_window"])
                if pd.isna(row["mz_window"]) or str(row["mz_window"]).strip().casefold() in {"", "auto"}:
                    mass, charge = _number(row["mw"]), _number(row["z"])
                    if math.isfinite(mass) and mass > 0 and math.isfinite(charge) and charge > 0:
                        width = compute_auto_mz_window(mass, charge, _text(row.get("formula", "")))
                if math.isfinite(centre) and math.isfinite(width) and width > 0:
                    frame.at[index, "mz_min"] = centre - width
                    frame.at[index, "mz_max"] = centre + width
        result_fields = (*RESULT_FIELDS, "extrapolated", "calibration_review_required", "conditions_status")
        self._results = []
        for index, row in frame.iterrows():
            result = {field: row[field] for field in result_fields if field in row and not _missing(row[field])}
            if result:
                source = row["raw_path"]
                try:
                    same_acquisition = (self._raw_path is not None and not _missing(source)
                                        and local_path(source).resolve() == self._raw_path)
                except (OSError, TypeError, ValueError):
                    same_acquisition = False
                if not same_acquisition:
                    result = {"status": "Stale"}
                    for field in result_fields:
                        if field in frame:
                            frame.at[index, field] = "Stale" if field == "status" else math.nan
            self._results.append(result)
        frame["raw_path"] = str(self._raw_path) if self._raw_path is not None else ""
        self._frame = frame[list(FIELDS) + [field for field in frame if field not in FIELDS]]
        for index in range(len(self._frame)):
            self._update_window(index)
        self._populate()
        self.changed.emit(self.get_targets())
        self.selected.emit(None)

    def update_results(self, records):
        """Attach CCS summaries without editing inputs or triggering processing signals."""
        records = records.to_dict("records") if isinstance(records, pd.DataFrame) else list(records)
        remaining = list(records)
        for row, target in self._frame.iterrows():
            for index, record in enumerate(remaining):
                if _matches(target, record):
                    self._results[row] = dict(record)
                    self._refresh_result_cells(row)
                    remaining.pop(index)
                    break

    def invalidate_results(self):
        """Mark prior CCS results stale after calibration changes, without feedback."""
        for row in range(len(self._results)):
            self._stale_result(row)

    def import_csv(self, path: Path | None = None):
        if path is None:
            filename, _ = QFileDialog.getOpenFileName(self, "Import targets for current acquisition", "", "CSV files (*.csv);;All files (*)")
            if not filename:
                return
            path = Path(filename)
        try:
            # Keep names, identifiers, formulae and optional calculation fields verbatim.
            # set_targets converts only the numeric target fields using Python's float parser.
            frame = pd.read_csv(local_path(path), dtype=str, keep_default_na=False)
            for field in ("ccs", "min_ms_intensity", "min_im_intensity", "peak_height_frac", "gas_mw"):
                if field in frame:
                    frame[field] = frame[field].mask(frame[field].eq(""), math.nan)
            if not {"mz_min", "mz_max"} <= set(frame.columns) and "mz" not in frame:
                raise ValueError("Target CSV needs mz_min and mz_max, or legacy mz and mz_window fields.")
            self.set_targets(frame)
        except (OSError, ValueError, pd.errors.ParserError) as exc:
            self._auto()
            self.error.emit(f"Cannot import targets: {exc}")

    def export_csv(self, path: Path | None = None):
        if path is None:
            filename, _ = QFileDialog.getSaveFileName(self, "Export targets", "targets.csv", "CSV files (*.csv)")
            if not filename:
                return
            path = filename
        staged = None
        try:
            destination = local_path(path).resolve()
            if not destination.suffix:
                destination = destination.with_suffix(".csv")
            if any(parent.suffix.casefold() == ".raw" for parent in (destination, *destination.parents)):
                raise ValueError("Choose an export location outside RAW acquisitions.")
            frame = self.get_targets()
            for row, result in enumerate(self._results):
                if not result:
                    continue
                for field in RESULT_FIELDS:
                    frame.at[row, field] = result.get(field, math.nan)
                for field in ("extrapolated", "calibration_review_required", "conditions_status"):
                    if field in result:
                        frame.at[row, field] = result[field]
            with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", newline="", dir=destination.parent,
                                             suffix=".csv", delete=False) as stream:
                staged = Path(stream.name)
                frame.to_csv(stream, index=False, float_format="%.17g")
            staged.replace(destination)
            return destination
        except (OSError, ValueError) as exc:
            self.error.emit(f"Cannot export targets: {exc}")
        finally:
            if staged is not None:
                staged.unlink(missing_ok=True)

    def _update_window(self, index):
        row = self._frame.iloc[index]
        bounds = _bounds(row["mz_min"], row["mz_max"])
        self._frame.at[index, "mz"] = (bounds[0] + bounds[1]) / 2 if bounds else math.nan
        self._frame.at[index, "mz_window"] = (bounds[1] - bounds[0]) / 2 if bounds else math.nan

    def _populate(self, selected_row=None):
        self._loading = True
        self.table.clearSelection()
        self.table.setCurrentItem(None)
        self.table.setRowCount(len(self._frame))
        for row in range(len(self._frame)):
            for column, field in enumerate(VISIBLE_FIELDS):
                text = "" if field == "drift_mode" else _text(self._frame.iloc[row][field])
                item = QTableWidgetItem(text)
                item.setToolTip(text)
                if field == "drift_mode":
                    item.setFlags((item.flags() | Qt.ItemIsUserCheckable) & ~Qt.ItemIsEditable)
                    item.setCheckState(Qt.Checked if self._frame.iloc[row][field] == "auto" else Qt.Unchecked)
                if column:
                    item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
                self.table.setItem(row, column, item)
            self._refresh_drift_cells(row)
            self._refresh_result_cells(row)
        self._loading = False
        if selected_row is not None:
            self.table.selectRow(selected_row)
        self.label.setText("Targets")

    def _edited(self, item):
        if self._loading or item.column() >= len(VISIBLE_FIELDS):
            return
        field, text = VISIBLE_FIELDS[item.column()], item.text().strip()
        value = ("auto" if item.checkState() == Qt.Checked else "manual") if field == "drift_mode" else text if field == "species" else _number(text)
        self._frame.at[item.row(), field] = value
        self._update_window(item.row())
        if field not in {"drift_mode", "drift_min", "drift_max"}:
            self._stale_result(item.row())
        if field == "drift_mode":
            self._refresh_drift_cells(item.row())
        self._loading = True
        item.setToolTip(item.text())
        self._loading = False
        self.changed.emit(self.get_targets())
        if any(index.row() == item.row() for index in self.table.selectionModel().selectedRows()):
            self._selection_changed()

    def _refresh_drift_cells(self, row):
        auto = self._frame.iloc[row]["drift_mode"] == "auto"
        with QSignalBlocker(self.table):
            for field in ("drift_min", "drift_max"):
                item = self.table.item(row, INPUT_COLUMNS[field])
                value = self._frame.iloc[row][field]
                numeric = _number(value)
                item.setText(f"{numeric:.1f}" if auto and math.isfinite(numeric) else _text(value))
                item.setToolTip(_text(value))
                item.setFlags(item.flags() & ~Qt.ItemIsEditable if auto else item.flags() | Qt.ItemIsEditable)

    def _refresh_result_cells(self, row):
        result = self._results[row]
        if result.get("status") == "Stale":
            tooltip = "Inputs changed; recalculate CCS for this target."
        elif result:
            tooltip = "\n".join((
                f"CCS peak: {_text(result.get('peak_ccs'))} Å²",
                f"Drift peak: {_text(result.get('peak_drift_time_ms'))} ms",
                f"Status: {_text(result.get('status'))}",
                f"Extrapolated: {_flag(result.get('extrapolated'))}",
                f"Calibration review required: {_flag(result.get('calibration_review_required'))}",
                f"Conditions: {_text(result.get('conditions_status')) or 'Unknown'}",
            ))
        else:
            tooltip = "Run CCS calculation for this target."
        with QSignalBlocker(self.table):
            for column, field in enumerate(RESULT_FIELDS, start=len(VISIBLE_FIELDS)):
                value = result.get(field)
                if field == "status":
                    text = _text(value)
                else:
                    numeric = _number(value)
                    text = (f"{numeric:.2f}" if field == "peak_ccs" else f"{numeric:.3f}") if math.isfinite(numeric) else ""
                item = self.table.item(row, column)
                if item is None:
                    item = QTableWidgetItem()
                    item.setFlags(item.flags() & ~Qt.ItemIsEditable)
                    self.table.setItem(row, column, item)
                item.setText(text)
                item.setToolTip(tooltip)
                if field != "status":
                    item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)

    def _stale_result(self, row):
        if self._results[row]:
            self._results[row] = {"status": "Stale"}
            self._refresh_result_cells(row)

    def _selection_changed(self):
        if self._loading:
            return
        self.label.setText("Targets")
        bounds = self.active_range
        self.selected.emit(bounds)
        if self.table.selectionModel().selectedRows() and bounds is None:
            self.error.emit("Target m/z bounds must be finite and increasing.")
        elif bounds is not None and self.active_drift_mode == "manual" and self.active_drift_range is None:
            self.error.emit("Manual drift bounds must be finite and increasing.")

    def _full_view(self):
        self._auto()
        self.full_view_requested.emit()

    def _auto(self):
        self._loading = True
        self.table.clearSelection()
        self.table.setCurrentItem(None)
        self._loading = False
        self._selection_changed()

    def _add(self):
        row = {field: math.nan for field in self._frame.columns}
        row["species"] = f"Target {len(self._frame) + 1}"
        row["drift_mode"] = "auto"
        row["raw_path"] = str(self._raw_path) if self._raw_path is not None else ""
        if self._default_range is not None:
            row["mz_min"], row["mz_max"] = self._default_range
        self._frame.loc[len(self._frame)] = row
        self._results.append({})
        self._update_window(len(self._frame) - 1)
        self._populate(len(self._frame) - 1)
        self.changed.emit(self.get_targets())

    def _delete(self):
        rows = self.table.selectionModel().selectedRows()
        if not rows:
            return
        self._frame = self._frame.drop(rows[0].row()).reset_index(drop=True)
        self._results.pop(rows[0].row())
        self._populate()
        self.changed.emit(self.get_targets())
        self.selected.emit(None)
