"""Blank optional formula cells should behave like an omitted formula column."""
from io import StringIO

import pandas as pd
import pytest

from deconvovo.imms_ccs_calibrate import _resolve_mz_window


@pytest.mark.parametrize("cell", ["", " ", "NA"])
def test_csv_blank_formula_uses_mass_fallback(cell):
    row = pd.read_csv(StringIO(f'formula,mz_window\n"{cell}",auto\n')).iloc[0]
    omitted = pd.Series({"mz_window": "auto"})
    assert _resolve_mz_window(row, 12327.0, 10) == _resolve_mz_window(omitted, 12327.0, 10)
