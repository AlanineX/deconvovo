"""Waters viewer configuration preservation and offline rendering regressions."""
from __future__ import annotations

import json
import os
from copy import deepcopy
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtWidgets import QApplication

import gui.panels.html_viewer as viewer_panel
from deconvovo.imms_html import _load_plot_config, plot_im_data
from gui.widgets.log_panel import LogPanel


@pytest.fixture
def config():
    cfg = json.loads((Path(__file__).parents[1] / "config" / "imms_plot_config.json").read_text(encoding="utf-8"))
    cfg["defaults"].update({"mz_bins": "native", "drift_bins": "native", "ms_mode": "profile", "retained_default": {"value": "α"}})
    cfg["figure"]["retained_figure"] = "keep"
    cfg["presets"]["retained_preset"] = ["keep"]
    cfg["initial_view"] = {"mz_range": [100, 103], "drift_range": [0, 0.3], "retained_view": "keep"}
    cfg["retained_root"] = {"value": "keep"}
    return cfg


@pytest.fixture
def panel(tmp_path, monkeypatch, config):
    app = QApplication.instance() or QApplication([])
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(config, ensure_ascii=False), encoding="utf-8")
    monkeypatch.setattr(viewer_panel, "_config_path", lambda: config_path)
    widget = viewer_panel.HtmlViewerPanel(LogPanel())
    yield widget
    widget.close()
    widget.deleteLater()
    app.processEvents()


def test_editor_preserves_unknown_nested_settings(panel, config, tmp_path):
    panel.spin_noise2d.setValue(7.25)
    # An incomplete preset edit must retain keys the editor does not expose.
    presets = deepcopy(config["presets"])
    del presets["retained_preset"]
    panel.presets_editor.setPlainText(json.dumps(presets))
    assert panel._save_config()
    saved = json.loads((tmp_path / "config.json").read_text(encoding="utf-8"))
    assert saved["defaults"]["retained_default"] == config["defaults"]["retained_default"]
    assert saved["figure"]["retained_figure"] == "keep"
    assert saved["presets"]["retained_preset"] == ["keep"]
    assert saved["initial_view"]["retained_view"] == "keep"
    assert saved["retained_root"] == config["retained_root"]
    assert saved["defaults"]["mz_bins"] == "native"
    assert saved["defaults"]["noise_2d"] == 7.25


def test_custom_smoothing_and_workspace_config_round_trip(panel, config):
    config["defaults"]["smooth_2d"] = {"method": "gaussian", "sigma": 0.1234, "sigma_dr": 0.22}
    config["defaults"].update({"z_scale": "linear", "render": "best", "ms_mode": "stick", "drift_bins": 75, "denoise": "off"})
    panel.set_config(config)
    collected = panel.current_config()
    assert collected["defaults"]["smooth_2d"] == config["defaults"]["smooth_2d"]
    for key in ("mz_bins", "z_scale", "render", "ms_mode", "drift_bins", "denoise"):
        assert collected["defaults"][key] == config["defaults"][key]
    assert collected["initial_view"] == config["initial_view"]


def test_invalid_editor_settings_do_not_save_or_launch(panel, tmp_path, monkeypatch):
    original = (tmp_path / "config.json").read_bytes()
    panel.presets_editor.setPlainText("{ invalid")
    assert panel._save_config() is False
    assert (tmp_path / "config.json").read_bytes() == original
    panel.pick_data.set_path(str(tmp_path))
    panel.pick_out.set_path(str(tmp_path / "output"))
    monkeypatch.setattr(viewer_panel, "Worker", lambda *args, **kwargs: pytest.fail("Invalid config launched a worker"))
    panel._run()


def test_invalid_initial_range_is_rejected(panel):
    panel.spin_mz_min.setValue(200)
    panel.spin_mz_max.setValue(100)
    with pytest.raises(ValueError, match="minimum must be lower"):
        panel.current_config()


def test_frozen_config_uses_user_storage_and_bundled_fallback(tmp_path, monkeypatch, config):
    import sys

    bundle = tmp_path / "bundle"
    (bundle / "config").mkdir(parents=True)
    bundled = bundle / "config" / "imms_plot_config.json"
    bundled.write_text(json.dumps(config), encoding="utf-8")
    user_dir = tmp_path / "user_config"
    monkeypatch.setattr(sys, "frozen", True, raising=False)
    monkeypatch.setattr(sys, "_MEIPASS", str(bundle), raising=False)
    monkeypatch.setattr(viewer_panel.QStandardPaths, "writableLocation", lambda kind: str(user_dir))
    assert viewer_panel._load_config() == config
    updated = deepcopy(config)
    updated["defaults"]["z_scale"] = "linear"
    viewer_panel._save_config(updated)
    assert viewer_panel._config_path() == user_dir / "imms_plot_config.json"
    assert viewer_panel._load_config() == updated
    assert json.loads(bundled.read_text(encoding="utf-8")) == config


def test_changed_settings_regenerate_old_html(panel, tmp_path, monkeypatch):
    import deconvovo.imms_plot

    data = tmp_path / "input"
    data.mkdir()
    (data / "sample_im.txt").write_text("100 0 1\n101 0 2\n", encoding="utf-8")
    out = tmp_path / "output"
    out.mkdir()
    (out / "sample_2d_imms.html").write_text("old viewer", encoding="utf-8")
    calls = []

    def fake_run(inputs, output, **kwargs):
        calls.append(kwargs)
        return [{"run_name": "sample", "status": ["IM"]}]

    monkeypatch.setattr(deconvovo.imms_plot, "run", fake_run)
    panel._do_plot(str(data), str(out), None, True, 1)
    assert calls[-1]["skip_existing"] is False
    assert json.loads(calls[-1]["config_path"].read_text(encoding="utf-8")) == panel.current_config()
    panel._do_plot(str(data), str(out), None, True, 1)
    assert calls[-1]["skip_existing"] is True
    panel.spin_noise2d.setValue(8)
    panel._do_plot(str(data), str(out), None, True, 1, panel.current_config())
    assert calls[-1]["skip_existing"] is False


def test_generated_viewer_embeds_plotly_and_custom_defaults(tmp_path, config):
    from plotly.offline import get_plotlyjs

    config["defaults"].update({"mz_bins": 250, "z_scale": "linear", "render": "best", "ms_mode": "stick", "drift_bins": 3, "noise_2d": 7.25})
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(config, ensure_ascii=False), encoding="utf-8")
    assert _load_plot_config(config_path)["defaults"]["retained_default"]["value"] == "α"
    im_file = tmp_path / "sample_im.txt"
    ms_file = tmp_path / "sample_ms.txt"
    im_file.write_text("100 0 1\n101 0 2\n102 1 3\n103 1 4\n100 2 2\n103 2 1\n", encoding="utf-8")
    ms_file.write_text("100 3\n101 2\n102 3\n103 5\n", encoding="utf-8")
    plot_im_data(im_file, ms_file, "sample", tmp_path, pusher_us=100, config_path=config_path)
    html = (tmp_path / "sample_2d_imms.html").read_text(encoding="utf-8")
    assert "<script src=" not in html
    assert get_plotlyjs() in html
    for value in ("250", "linear", "best", "stick", "3", "7.25"):
        assert f'<option value="{value}" selected>' in html
    for control in ("nbins", "dbins", "msmode", "zscale", "render", "denoise"):
        assert f'id="{control}"' in html
