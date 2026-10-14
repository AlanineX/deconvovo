"""Current-view capture and numeric range controls on the real offline viewer."""
from __future__ import annotations

import json
import os
from pathlib import Path
import re
import subprocess
import sys

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("DECONVOVO_DISABLE_WEBENGINE", "1")

from PySide6.QtWidgets import QApplication
from PySide6.QtQml import QJSEngine

from gui.widgets.waters_viewer import WatersViewerWidget


@pytest.fixture
def widget():
    app = QApplication.instance() or QApplication([])
    viewer = WatersViewerWidget()
    yield viewer
    viewer.wait_for_workers()
    viewer.close()
    viewer.deleteLater()
    app.processEvents()


@pytest.mark.parametrize("bounds", [[1, 1], [2, 1], [0, float("inf")], [float("nan"), 2], [1]])
def test_numeric_range_requires_finite_ordered_bounds(widget, bounds):
    with pytest.raises(ValueError):
        widget.apply_config({"initial_view": {"mz_range": bounds}})


def test_current_export_requires_a_live_viewer_and_known_format(widget):
    with pytest.raises(ValueError, match="Unknown"):
        widget.export_current("raw_csv")
    with pytest.raises(ValueError):
        widget.export_current("mz_csv")


def test_export_normalizes_pasted_paths_and_keeps_source_guards(widget, tmp_path):
    folder = tmp_path / "Folder with spaces"
    folder.mkdir()
    target = folder / "current plot"
    assert widget._output_path('"' + str(target) + '"', {".png"}) == target.with_suffix(".png")
    source = folder / "sample_im.txt"
    source.write_text("1000 0 10\n")
    widget._record = {"im": source}
    alias = folder / "alias.csv"
    try:
        alias.symlink_to(source)
    except OSError:
        alias = source
    with pytest.raises(ValueError, match="source MS or IM"):
        widget._output_path(alias, {".csv"})
    with pytest.raises(ValueError, match="outside RAW"):
        widget._output_path(folder / "original.RAW" / "current.png", {".png"})
    assert source.read_text() == "1000 0 10\n"


def test_csv_uses_current_display_values_bin_units_and_stick_breaks(tmp_path):
    from deconvovo.imms_html import plot_im_data

    app = QApplication.instance() or QApplication([])
    im = tmp_path / "sample_im.txt"
    np.savetxt(im, [[1000, 0, 10], [1500, 1, 20], [2000, 2, 30]])
    plot_im_data(im, None, "sample", tmp_path)
    html = (tmp_path / "sample_2d_imms.html").read_text(encoding="utf-8")
    function = re.search(r"function watersCurrentCSV\(kind\) \{.*?\n\}\nasync function", html, re.S).group(0)
    function = function.rsplit("\nasync function", 1)[0]
    engine = QJSEngine()
    value = engine.evaluate("""
      var pusherUs=0;
      var document={getElementById:function(){return {value:'sqrt'}}};
      function getVisibleRange(){return {mzLo:1450,mzHi:1700,drLo:1.2,drHi:1.8,
        mBl:1,mBh:2,dBl:1,dBh:1};}
      var P={data:[{x:[10,0,8],y:[0.5,1.5,2.5]},
        {x:[1400,1400,null,1600,1600,null],y:[0,9,null,0,5,null]},
        {x:[1167,1500,1833],y:[0.5,1.5,2.5],z:[[1,2,3],[4,5,6],[7,8,9]]}]};
    """ + function)
    assert not value.isError(), value.toString()
    assert engine.evaluate("watersCurrentCSV('drift_csv')").toString() == "drift_bin,display_value\n1.5,0"
    assert engine.evaluate("watersCurrentCSV('mz_csv')").toString() == "mz,display_value\n,\n1600,0\n1600,5\n,"
    assert engine.evaluate("watersCurrentCSV('heatmap_csv').join('')").toString() == (
        "mz,drift_bin,display_value,intensity_transform\n1500,1.5,5,sqrt\n1833,1.5,6,sqrt\n")
    app.processEvents()



@pytest.mark.parametrize("module_path", ["gui/widgets/waters_viewer.py", "portable_waters_convert/gui/widgets/waters_viewer.py"])
def test_native_current_exports_capture_settled_controls_and_numeric_ranges(tmp_path, module_path, monkeypatch):
    """The real Qt canvas captures native arrays without creating an HTML viewer."""
    import csv
    import importlib.util
    import time
    from PIL import Image
    import deconvovo.io as io

    app = QApplication.instance() or QApplication([])
    spec = importlib.util.spec_from_file_location("_native_capture_test", Path(__file__).parents[1] / module_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    monkeypatch.setattr(io, "find_pusher_for_run", lambda *a: 100.0)
    config = json.loads((Path(__file__).parents[1] / "config/imms_plot_config.json").read_text(encoding="utf-8"))
    config["defaults"].update({"mz_bins": 4, "drift_bins": "native", "smooth_2d": {"method": "raw"},
        "smooth_drift": {"method": "raw"}, "noise_2d": 0, "denoise": "off", "noise_mz": 0,
        "noise_drift": 0, "ms_mode": "profile", "z_scale": "linear"})
    config["initial_view"] = {"mz_range": [1000, 2000], "drift_range": [0, 0.4]}
    im, ms = tmp_path / "sample_im.txt", tmp_path / "sample_ms.txt"
    points = np.array([[950, 0, 50], [1000, 0, 4], [1250, 0, 9], [1500, 1, 25],
        [1750, 2, 100], [2000, 3, 16], [2100, 3, 50]])
    np.savetxt(im, points)
    np.savetxt(ms, points[:, [0, 2]])
    original = im.read_bytes(), ms.read_bytes()
    out = tmp_path / "captures"
    out.mkdir()
    viewer = module.WatersViewerWidget()
    viewer.resize(900, 650)
    viewer.show()
    def wait(predicate):
        deadline = time.monotonic() + 15
        while time.monotonic() < deadline:
            app.processEvents()
            if predicate():
                return
            time.sleep(0.01)
        raise AssertionError("Native rendering or capture timed out: " + viewer.status.text())
    try:
        viewer.set_run({"run_name": "sample", "im": im, "ms": ms}, config)
        wait(lambda: viewer.is_ready)
        first = out / "before.png"
        viewer.export_current("full_png", first)
        wait(first.is_file)
        before = first.read_bytes()
        config["defaults"].update({"mz_bins": 8, "drift_bins": 2, "z_scale": "sqrt", "noise_2d": 30,
            "denoise": "off", "noise_drift": 20, "noise_mz": 10})
        config["initial_view"] = {"mz_range": [1200, 1800], "drift_range": [0.1, 0.3]}
        viewer.apply_config(config)
        # Export immediately: the queued control change must be included.
        heatmap = out / "heatmap.csv"
        viewer.export_current("heatmap_csv", heatmap)
        wait(lambda: heatmap.is_file() and viewer.is_ready)
        snapshot = viewer.snapshot
        assert len(snapshot.mz_centres) == 8 and len(snapshot.drift_centres) == 2
        assert viewer.ax_heatmap.get_xlim() == pytest.approx((1200, 1800))
        assert viewer.ax_ms.get_xlim() == pytest.approx((1200, 1800))
        assert viewer.ax_heatmap.get_ylim() == pytest.approx((0.1, 0.3))
        assert viewer.ax_drift.get_ylim() == pytest.approx((0.1, 0.3))
        with heatmap.open(encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        expected = [(m, d, value) for d, values in zip(snapshot.drift_centres[snapshot.visible_drift_slice],
            snapshot.display2d[snapshot.visible_drift_slice, snapshot.visible_mz_slice])
            for m, value in zip(snapshot.mz_centres[snapshot.visible_mz_slice], values)]
        actual = [tuple(float(row[key]) for key in ("mz", "drift_time_ms", "display_value")) for row in rows]
        np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-12)
        assert all(row["intensity_transform"] == "sqrt" for row in rows)
        for kind in ("mz_csv", "drift_csv", "full_png", "heatmap_png", "mz_png", "drift_png"):
            destination = out / kind.replace("_", ".")
            viewer.export_current(kind, destination)
            wait(destination.is_file)
            if kind.endswith("png"):
                with Image.open(destination) as image:
                    image.verify()
                    assert image.width >= 500 and image.height >= 300
        assert (out / "full.png").read_bytes() != before
        with (out / "mz.csv").open(encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        expected = [(x, y) for x, y in zip(snapshot.ms_mz, snapshot.ms_intensity) if 1200 <= x <= 1800]
        np.testing.assert_allclose([(float(r["mz"]), float(r["display_value"])) for r in rows], expected, rtol=1e-12)
        # Native bounds can expand directly from TXT, including beyond initial bounds.
        config["initial_view"] = {"mz_range": [900, 2200], "drift_range": [0, 0.4]}
        viewer.apply_config(config)
        wait(lambda: viewer.is_ready and viewer.ax_heatmap.get_xlim() == (900, 2200))
        config["initial_view"] = {"mz_range": None, "drift_range": None}
        viewer.apply_config(config)
        wait(lambda: viewer.is_ready and viewer.ax_heatmap.get_xlim() == (950, 2100))
        assert not hasattr(viewer, "web_view")
        assert not list(tmp_path.glob("*.html"))
        assert (im.read_bytes(), ms.read_bytes()) == original
    finally:
        viewer.close()
        viewer.wait_for_workers()
        viewer.deleteLater()
        app.processEvents()
