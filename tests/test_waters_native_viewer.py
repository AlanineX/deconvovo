"""Native Waters canvas behavior, real shared axes, and asynchronous selection."""
from __future__ import annotations

from copy import deepcopy
import json
import os
from pathlib import Path
import threading
import time
import warnings
from types import SimpleNamespace

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtCore import QCoreApplication, QEvent
from PySide6.QtWidgets import QApplication
import deconvovo.imms_view as engine
import deconvovo.io as io
from gui.widgets.waters_viewer import WatersViewerWidget


@pytest.fixture
def config():
    config = json.loads((Path(__file__).parents[1] / "config/imms_plot_config.json").read_text(encoding="utf-8"))
    config["defaults"].update({"mz_bins": 8, "drift_bins": "native", "smooth_2d": {"method": "raw"},
        "smooth_drift": {"method": "raw"}, "noise_2d": 0, "denoise": "off", "noise_mz": 0,
        "noise_drift": 0, "ms_mode": "profile", "z_scale": "linear", "colormap": "Viridis"})
    config["initial_view"] = {"mz_range": [100, 104], "drift_range": [0, .4]}
    return config


@pytest.fixture
def widget(monkeypatch):
    app = QApplication.instance() or QApplication([])
    monkeypatch.setattr(io, "find_pusher_for_run", lambda *args: 100.0)
    panel = WatersViewerWidget()
    panel.resize(640, 480)
    panel.show()
    yield panel
    panel.close()
    assert panel.wait_for_workers(10000)
    panel.deleteLater()
    QCoreApplication.sendPostedEvents(None, QEvent.DeferredDelete)
    app.processEvents()


def record(tmp_path, name="sample"):
    im, ms = tmp_path / f"{name}_im.txt", tmp_path / f"{name}_ms.txt"
    points = np.array([[99, 0, 1], [100, 0, 2], [101, 0, 4], [102, 1, 16],
                       [103, 2, 8], [104, 3, 4], [105, 3, 2]])
    np.savetxt(im, points)
    np.savetxt(ms, points[:, [0, 2]])
    return {"run_name": name, "im": im, "ms": ms}


def wait(predicate, timeout=6):
    deadline = time.monotonic() + timeout
    while not predicate() and time.monotonic() < deadline:
        QApplication.instance().processEvents()
        time.sleep(.01)
    QApplication.instance().processEvents()
    assert predicate()


def test_native_viewer_reads_txt_without_opening_or_parsing_html(widget, config, tmp_path, monkeypatch):
    item = record(tmp_path)
    html = tmp_path / "invalid.html"
    html.write_text("this is deliberately not a viewer", encoding="utf-8")
    item["viewer"] = html
    from gui.widgets.waters_viewer import QDesktopServices
    monkeypatch.setattr(QDesktopServices, "openUrl", lambda *args: pytest.fail("Native plot opened a browser"))
    read_text = Path.read_text
    def guarded(path, *args, **kwargs):
        if path == html:
            pytest.fail("Native plot parsed HTML instead of reading its TXT source")
        return read_text(path, *args, **kwargs)
    monkeypatch.setattr(Path, "read_text", guarded)
    widget.set_run(item, config)
    wait(lambda: widget.is_ready)
    assert not hasattr(widget, "web_view")
    assert widget.snapshot.display2d.shape == (4, 8)
    assert widget.im_stack.currentWidget() is widget._plot_page


def test_shared_axis_projection_matches_independent_physical_oracle(widget, config, tmp_path):
    widget.set_run(record(tmp_path), config)
    wait(lambda: widget.is_ready)
    widget.canvas.draw()
    for fraction in np.linspace(0, 1, 5):
        mz = 100 + fraction * 4
        drift = fraction * .4
        x_expected = widget.ax_heatmap.bbox.x0 + fraction * widget.ax_heatmap.bbox.width
        y_expected = widget.ax_heatmap.bbox.y0 + fraction * widget.ax_heatmap.bbox.height
        x_heat = widget.ax_heatmap.transData.transform((mz, drift))[0]
        x_mass = widget.ax_ms.transData.transform((mz, 0))[0]
        y_heat = widget.ax_heatmap.transData.transform((mz, drift))[1]
        y_drift = widget.ax_drift.transData.transform((0, drift))[1]
        np.testing.assert_allclose([x_heat, x_mass], x_expected, atol=1e-7)
        np.testing.assert_allclose([y_heat, y_drift], y_expected, atol=1e-7)
    assert widget.ax_drift.xaxis_inverted()
    assert widget.canvas.height() >= 250
    np.testing.assert_array_equal(widget.ax_heatmap.images[0].get_array(), widget.snapshot.display2d)
    np.testing.assert_array_equal(widget.ax_drift.lines[0].get_ydata(), widget.snapshot.drift_centres)
    assert not widget.ax_drift.xaxis.get_offset_text().get_visible()
    assert not widget.ax_drift.get_xticklabels()
    assert widget.ax_ms.yaxis.get_major_formatter()(10000) == "10k"
    assert widget.mobility_plot.ax().get_xticklabels()


def test_numeric_ranges_override_zoom_but_filters_preserve_it_and_reset_keeps_controls(widget, config, tmp_path):
    widget.set_run(record(tmp_path), config)
    wait(lambda: widget.is_ready)
    widget.ax_heatmap.set_xlim(101, 103)
    widget.ax_heatmap.set_ylim(.1, .3)
    wait(lambda: widget.is_ready and widget.snapshot.mz_range == (101, 103))
    changed = deepcopy(config)
    changed["defaults"]["z_scale"] = "sqrt"
    widget.apply_config(changed)
    wait(lambda: widget.is_ready and widget.snapshot.transform == "sqrt")
    assert widget.snapshot.mz_range == (101, 103)
    assert widget.snapshot.drift_range == pytest.approx((.1, .3))
    changed["initial_view"]["mz_range"] = [99, 105]
    widget.apply_config(changed)
    wait(lambda: widget.is_ready and widget.snapshot.mz_range == (99, 105))
    widget.ax_heatmap.set_xlim(100, 102)
    wait(lambda: widget.is_ready and widget.snapshot.mz_range == (100, 102))
    widget.reset_view()
    wait(lambda: widget.is_ready and widget.snapshot.mz_range == (99, 105))
    assert widget.snapshot.transform == "sqrt"


def test_scroll_zoom_recomputes_both_linked_profiles(widget, config, tmp_path):
    widget.set_run(record(tmp_path), config)
    wait(lambda: widget.is_ready)
    widget._scroll_zoom(SimpleNamespace(inaxes=widget.ax_heatmap, button="up", xdata=102., ydata=.2))
    wait(lambda: widget.is_ready and widget.snapshot.mz_range != (100, 104))
    assert widget.snapshot.mz_range == pytest.approx((100.4, 103.6))
    assert widget.snapshot.drift_range == pytest.approx((.04, .36))
    assert widget.ax_ms.get_xlim() == pytest.approx(widget.snapshot.mz_range)
    assert widget.ax_drift.get_ylim() == pytest.approx(widget.snapshot.drift_range)


def test_one_cached_data_load_handles_debounced_control_changes(widget, config, tmp_path, monkeypatch):
    original = engine.load_view_data
    calls = []
    def load(*args, **kwargs):
        calls.append(args[0])
        return original(*args, **kwargs)
    monkeypatch.setattr(engine, "load_view_data", load)
    widget.set_run(record(tmp_path), config)
    wait(lambda: widget.is_ready)
    for bins in (16, 32, 64):
        config["defaults"]["mz_bins"] = bins
        widget.apply_config(config)
    wait(lambda: widget.is_ready and len(widget.snapshot.mz_centres) == 64)
    assert len(calls) == 1


def test_stale_slow_selection_cannot_replace_new_run(widget, config, tmp_path, monkeypatch):
    slow, latest = record(tmp_path, "old"), record(tmp_path, "new")
    entered, release = threading.Event(), threading.Event()
    original = engine.load_view_data
    def load(im, *args, **kwargs):
        if Path(im) == slow["im"]:
            entered.set()
            release.wait(3)
        return original(im, *args, **kwargs)
    monkeypatch.setattr(engine, "load_view_data", load)
    widget.set_run(slow, config)
    wait(entered.is_set)
    widget.set_run(latest, config)
    release.set()
    wait(lambda: widget.is_ready)
    assert widget._record["run_name"] == "new"
    assert widget.snapshot is not None
    assert not widget._tasks


def test_raw_only_and_invalid_txt_have_clear_states(widget, config, tmp_path):
    raw = tmp_path / "original.raw"
    raw.mkdir()
    widget.set_run({"run_name": "raw", "raw": raw}, config)
    assert not widget.is_ready
    assert "Process & view" in widget.empty.text()
    item = record(tmp_path)
    item["im"].write_text("invalid", encoding="utf-8")
    widget.set_run(item, config)
    wait(lambda: bool(widget.status.text()) and "invalid" in widget.status.text().lower())
    assert not widget.is_ready


def test_ready_signal_and_tabs_preserve_run_inspection_choice(widget, config, tmp_path):
    states = []
    widget.ready_changed.connect(states.append)
    widget.tabs.setCurrentIndex(2)
    widget.set_run(record(tmp_path), config)
    wait(lambda: widget.is_ready)
    assert states[-1] is True and widget.tabs.currentIndex() == 2
    widget.clear("Select another run")
    assert states[-1] is False
    assert widget.tabs.currentIndex() == 2


def test_one_drift_cell_retains_physical_extent_when_viewport_crops_it(widget, config, tmp_path):
    config["defaults"]["drift_bins"] = 1
    config["initial_view"]["drift_range"] = [.1, .3]
    widget.set_run(record(tmp_path), config)
    wait(lambda: widget.is_ready)
    assert widget.snapshot.drift_centres == pytest.approx([.2])
    assert widget.ax_heatmap.images[0].get_extent()[2:] == pytest.approx((0, .4))
    assert widget.ax_heatmap.get_ylim() == pytest.approx((.1, .3))


def test_toolbar_keeps_navigation_and_routes_home_to_current_config(widget, config, tmp_path):
    names = [action.text() for action in widget.toolbar.actions()]
    assert "Save" not in names and "Subplots" not in names and "Customize" not in names
    assert {"Home", "Back", "Forward", "Pan", "Zoom"} <= set(names)
    widget.set_run(record(tmp_path), config)
    wait(lambda: widget.is_ready)
    widget.ax_heatmap.set_xlim(101, 103)
    wait(lambda: widget.is_ready and widget.snapshot.mz_range == (101, 103))
    widget.toolbar.home()
    wait(lambda: widget.is_ready and widget.snapshot.mz_range == (100, 104))


def test_config_change_immediately_invalidates_old_render_and_keeps_multiple_exports(widget, config, tmp_path, monkeypatch):
    widget.set_run(record(tmp_path), config)
    wait(lambda: widget.is_ready)
    original_render = engine.ViewData.render
    entered, release = threading.Event(), threading.Event()
    def slow(data, cfg=None, viewport=None):
        if (cfg or {}).get("defaults", {}).get("mz_bins") == 16:
            entered.set()
            release.wait(3)
        return original_render(data, cfg, viewport)
    monkeypatch.setattr(engine.ViewData, "render", slow)
    config["defaults"]["mz_bins"] = 16
    widget.apply_config(config)
    assert not widget.is_ready
    wait(entered.is_set)
    config["defaults"]["mz_bins"] = 32
    widget.apply_config(config)
    assert not widget.is_ready
    a, b = tmp_path / "heatmap.csv", tmp_path / "drift.csv"
    widget.export_current("heatmap_csv", a)
    widget.export_current("drift_csv", b)
    release.set()
    wait(lambda: widget.is_ready and a.is_file() and b.is_file())
    assert len(widget.snapshot.mz_centres) == 32
    assert not widget._pending_exports


def test_light_theme_refresh_changes_presentation_without_rendering_data_again(widget, config, tmp_path):
    from gui.widgets.plot_canvas import set_plot_dark
    widget.set_run(record(tmp_path), config)
    wait(lambda: widget.is_ready)
    snapshot = widget.snapshot
    limits = widget.ax_heatmap.get_xlim(), widget.ax_heatmap.get_ylim()
    try:
        set_plot_dark(False)
        widget.refresh_theme()
        assert widget.ax_ms.lines[0].get_color() == "#257a72"
        assert widget.figure.get_facecolor()[0] > .9
        assert widget.snapshot is snapshot
        assert (widget.ax_heatmap.get_xlim(), widget.ax_heatmap.get_ylim()) == limits
    finally:
        set_plot_dark(True)
        widget.refresh_theme()


@pytest.mark.parametrize("transform", ["linear", "sqrt", "cbrt", "p4", "p5", "log2", "log10"])
def test_native_colorbar_and_png_labels_render_without_missing_sans_glyphs(widget, config, tmp_path, transform):
    import matplotlib
    config["defaults"]["z_scale"] = transform
    with matplotlib.rc_context({"font.family": "Liberation Sans"}), warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always")
        widget.set_run(record(tmp_path), config)
        wait(lambda: widget.is_ready)
        snapshot = widget.snapshot
        widget.canvas.draw()
        widget.refresh_theme()
        widget.canvas.draw()
        target = tmp_path / f"{transform}.png"
        widget.export_current("heatmap_png", target)
        assert target.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
        assert widget.snapshot is snapshot and snapshot.transform == transform
        assert not widget.ax_colorbar.yaxis.label.get_usetex()
        if transform == "cbrt":
            assert snapshot.intensity_label == "∛I"
        assert not [str(item.message) for item in captured if "Glyph" in str(item.message)]
