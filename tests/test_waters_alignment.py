"""Coordinate and viewport contracts shared by the Waters heatmap and profiles."""
from __future__ import annotations

import json
import os
from pathlib import Path
import re

import numpy as np
import pandas as pd
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtQml import QJSEngine
from PySide6.QtWidgets import QApplication

from deconvovo.imms_html import plot_im_data


@pytest.fixture(params=[None, 100.0])
def viewer(tmp_path: Path, request):
    app = QApplication.instance() or QApplication([])
    cfg = json.loads((Path(__file__).parents[1] / "config/imms_plot_config.json").read_text(encoding="utf-8"))
    cfg["defaults"].update({"mz_bins": 2, "smooth_2d": {"method": "raw"}, "smooth_drift": {"method": "raw"}, "noise_2d": 0, "denoise": "off"})
    cfg["initial_view"] = {"mz_range": [1000, 2000], "drift_range": [0, 0.3] if request.param else [0, 3]}
    config = tmp_path / "config.json"
    config.write_text(json.dumps(cfg))
    im = tmp_path / "sample_im.txt"
    ms = tmp_path / "sample_ms.txt"
    np.savetxt(im, [[1000, 0, 10], [1500, 1, 20], [2000, 2, 30]])
    np.savetxt(ms, [[1000, 10], [1500, 20], [2000, 30]])
    plot_im_data(im, ms, "sample", tmp_path, config_path=config, pusher_us=request.param)
    content = (tmp_path / "sample_2d_imms.html").read_text(encoding="utf-8")
    bootstrap = "var H=" + re.search(r"\nvar H=(.*?)\nvar figData=", content, re.S).group(1)
    engine = QJSEngine()
    value = engine.evaluate("var document={getElementById:function(id){return {value:id==='nbins'?'2':'native'}}};function applyAll(){};" + bootstrap)
    assert not value.isError(), value.toString()
    yield engine, content, tmp_path, request.param
    app.processEvents()


def read_json(engine, expression):
    value = engine.evaluate(f"JSON.stringify({expression})")
    assert not value.isError(), value.toString()
    return json.loads(value.toString())


def test_native_centres_and_histogram_boundaries_survive_rebin(viewer):
    engine, _content, directory, pusher = viewer
    initial = read_json(engine, "({H:H,DR:DR,MZ:MZ})")
    unit = pusher / 1000 if pusher else 1
    assert initial["DR"] == pytest.approx([0.5 * unit, 1.5 * unit, 2.5 * unit])
    assert initial["H"] == [[10, 0], [0, 20], [0, 30]]
    if pusher:
        exported = pd.read_csv(directory / "sample_2d_imms.csv")
        assert exported["drift_time_ms"].tolist() == pytest.approx(initial["DR"])
    value = engine.evaluate("rebinFull()")
    assert not value.isError(), value.toString()
    after = read_json(engine, "({H:H,DR:DR,MZ:MZ})")
    assert after == initial
    assert sum(map(sum, after["H"])) == 60


def test_partially_visible_cells_are_included_and_outside_drift_is_empty(viewer):
    engine, content, _directory, pusher = viewer
    function = re.search(r"function getVisibleRange\(\) \{.*?\n\}\n\n// === Apply 2D", content, re.S).group(0).split("\n\n// === Apply 2D")[0]
    unit = pusher / 1000 if pusher else 1
    # Neither m/z centre lies inside this viewport, but portions of both cells do.
    code = f"var P={{layout:{{xaxis2:{{range:[1300,1700]}},yaxis2:{{range:[{1.2 * unit},{1.8 * unit}]}}}}}};" + function
    value = engine.evaluate(code)
    assert not value.isError(), value.toString()
    selected = read_json(engine, "getVisibleRange()")
    assert (selected["mBl"], selected["mBh"]) == (0, 1)
    assert (selected["dBl"], selected["dBh"]) == (1, 1)
    assert (selected["ndBl"], selected["ndBh"]) == (1, 1)
    engine.evaluate(f"P.layout.yaxis2.range=[{4 * unit},{5 * unit}]")
    empty = read_json(engine, "getVisibleRange()")
    assert empty["dBh"] == -1
    assert empty["ndBh"] == -1
    engine.evaluate("P.layout.xaxis2.range=[2100,2200]")
    empty = read_json(engine, "getVisibleRange()")
    assert empty["mBh"] == -1
    # The editable axes can differ from Plotly's resolved matched-axis ranges.
    engine.evaluate(f"P._fullLayout={{xaxis2:{{range:[1300,1700]}},yaxis2:{{range:[{1.2 * unit},{1.8 * unit}]}}}}")
    rendered = read_json(engine, "getVisibleRange()")
    assert (rendered["mBl"], rendered["mBh"]) == (0, 1)
    assert (rendered["ndBl"], rendered["ndBh"]) == (1, 1)


def test_linked_axes_share_domains_and_initial_ranges(viewer):
    _engine, content, _directory, _pusher = viewer
    figure = json.loads(re.search(r"\nvar figData=(.*?);\nvar P=", content, re.S).group(1))
    layout = figure["layout"]
    assert layout["xaxis2"]["domain"] == layout["xaxis4"]["domain"]
    assert layout["xaxis2"]["range"] == layout["xaxis4"]["range"]
    assert layout["xaxis2"]["matches"] == "x4"
    assert layout["yaxis"]["domain"] == layout["yaxis2"]["domain"]
    assert layout["yaxis"]["range"] == layout["yaxis2"]["range"]
    assert layout["yaxis2"]["matches"] == "y"
