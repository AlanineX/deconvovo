"""Waters GUI workflow boundaries, table edits, demos, and background jobs."""
from __future__ import annotations

import json
import os
from pathlib import Path
import threading
import time

import pandas as pd
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtWidgets import QApplication

from gui.panels.waters_workspace import WatersWorkspacePanel, discover_runs
from gui.widgets.log_panel import LogPanel


@pytest.fixture
def panel():
    app = QApplication.instance() or QApplication([])
    widget = WatersWorkspacePanel(LogPanel())
    yield widget
    widget.wait_for_worker()
    widget.close()
    app.processEvents()
    widget.deleteLater()
    app.processEvents()


def _pairs(folder: Path, name="sample"):
    folder.mkdir(parents=True, exist_ok=True)
    for suffix in ("ms", "im"):
        (folder / f"{name}_{suffix}.txt").write_text("100 1\n", encoding="utf-8")
    return folder


def _wait_until(predicate, timeout=3.0):
    app = QApplication.instance()
    deadline = time.monotonic() + timeout
    while not predicate() and time.monotonic() < deadline:
        app.processEvents()
        time.sleep(0.005)
    app.processEvents()
    assert predicate(), "Qt worker did not finish in time"


def test_discovery_handles_nested_raw_single_raw_and_existing_results(tmp_path):
    raw = tmp_path / "acquisitions" / "condition" / "sample.raw"
    raw.mkdir(parents=True)
    output = tmp_path / "results"
    _pairs(output / "_converted")
    viewer = output / "sample_2d_imms.html"
    viewer.write_text("<html>offline viewer</html>", encoding="utf-8")
    assert discover_runs(raw)[0]["raw"] == raw
    rows = discover_runs(tmp_path / "acquisitions", output)
    assert len(rows) == 1
    assert set(rows[0]) == {"run_name", "raw", "ms", "im", "viewer"}
    assert rows[0]["viewer"] == viewer


def test_real_demo_loads_full_batch_settings_and_ccs_tables(panel):
    if panel.combo_demo.count() < 2:
        pytest.skip("Real production corpus is not present in this checkout")
    profile = panel.combo_demo.currentData()
    assert profile["id"] == "protein_adp"
    panel._load_demo()
    assert panel.run_table.rowCount() == 5
    assert panel.calibration_panel.table.table.rowCount() == 25
    assert panel.analyte_panel.table.table.rowCount() == 34
    assert not panel.chk_ccs.isChecked()
    assert "DeconVoVo" in panel.pick_out.path()
    assert Path(panel.pick_out.path()).is_relative_to(Path.home())
    demo_config = json.loads((panel._manifest_path.parent / profile["viewer_config"]).read_text(encoding="utf-8"))
    assert panel.viewer_panel.current_config()["initial_view"]["mz_range"] == demo_config["initial_view"]["mz_range"]
    assert not panel.calibration_panel.pick_data.isVisible()
    assert not panel.analyte_panel.pick_cal.isVisible()
    panel.combo_demo.setCurrentIndex(next(i for i in range(panel.combo_demo.count()) if (panel.combo_demo.itemData(i) or {}).get("id") == "hsp104"))
    panel._load_demo()
    assert panel.run_table.rowCount() == 2
    assert panel.calibration_panel.table.to_dataframe().empty
    assert panel.analyte_panel.table.to_dataframe().empty
    assert not panel.chk_ccs.isChecked()


def test_shared_paths_follow_source_output_and_method(panel, tmp_path):
    source = _pairs(tmp_path / "converted")
    output = tmp_path / "results"
    panel.pick_source.set_path(str(source))
    panel.pick_out.set_path(str(output))
    panel.combo_method.setCurrentText("twostep")
    for editor in (panel.calibration_panel, panel.analyte_panel):
        assert editor.pick_data.path() == str(source)
        assert editor.pick_out.path() == str(output / "ccs_twostep")
        assert editor.combo_method.currentText() == "twostep"
    assert panel.tabs.count() == 2
    assert panel.chk_reuse.isChecked()


def test_converted_processing_is_one_background_job_and_locks_inputs(panel, tmp_path, monkeypatch):
    import deconvovo.waters_workflow as workflow
    source = _pairs(tmp_path / "converted")
    out = tmp_path / "results"
    panel.pick_source.set_path(str(source))
    panel.pick_out.set_path(str(out))
    started = threading.Event()
    release = threading.Event()
    captured = []
    def fake_run(source_path, output_path, **kwargs):
        captured.append((source_path, output_path, kwargs))
        assert kwargs["convert"] is False
        assert kwargs["html"] is True
        assert kwargs["skip_existing"] is True
        assert json.loads(kwargs["config_path"].read_text(encoding="utf-8"))["defaults"]
        started.set()
        release.wait(2.0)
        output_path.mkdir()
        viewer = output_path / "sample_2d_imms.html"
        viewer.write_text("<html>offline</html>", encoding="utf-8")
        kwargs["progress"](1, 1, "Viewer ready: sample")
        return {"converted_dir": str(source_path), "viewers": [str(viewer)], "runs": [], "status": "complete"}
    monkeypatch.setattr(workflow, "run", fake_run)
    panel.process_and_view()
    try:
        _wait_until(started.is_set)
        assert panel.is_busy
        assert not panel.btn_run.isEnabled()
        assert not panel.btn_ccs.isEnabled()
        assert not panel.sources_group.isEnabled()
        assert not panel.ccs_editors.isEnabled()
        panel.process_and_view()
        assert len(captured) == 1
    finally:
        release.set()
    _wait_until(lambda: not panel.is_busy)
    assert panel.btn_run.isEnabled()
    assert panel.btn_ccs.isEnabled()
    assert panel.sources_group.isEnabled()
    assert panel.job_progress.value() == 1
    assert "1 offline viewers" in panel.job_status.text()
    assert panel.btn_viewer.isEnabled()
    assert not captured[0][2]["config_path"].exists()


def test_ccs_snapshots_table_edits_and_resolves_raw_paths_without_rewriting_source(panel, tmp_path, monkeypatch):
    import deconvovo.waters_workflow as workflow
    source = _pairs(tmp_path / "converted")
    original = tmp_path / "tables"
    original.mkdir()
    cal_path = original / "calibrants.csv"
    ana_path = original / "analytes.csv"
    cal_rows = [{"name": str(i), "mw": 100, "z": 1, "mz": 101, "ccs": 120, "raw_path": "raw/sample.raw"} for i in range(3)]
    pd.DataFrame(cal_rows).to_csv(cal_path, index=False)
    pd.DataFrame([{"species": "target", "mw": 100, "z": 1, "mz": 101, "raw_path": "raw/sample.raw"}]).to_csv(ana_path, index=False)
    saved_cal = cal_path.read_bytes()
    panel.pick_source.set_path(str(source))
    panel.pick_out.set_path(str(tmp_path / "output"))
    panel.calibration_panel.pick_csv.set_path(str(cal_path))
    panel.analyte_panel.pick_ana.set_path(str(ana_path))
    table = panel.analyte_panel.table.table
    mz_col = next(i for i in range(table.columnCount()) if table.horizontalHeaderItem(i).text() == "mz")
    table.item(0, mz_col).setText("202")
    captured = []
    def fake_run(_source, _output, **kwargs):
        assert kwargs["html"] is False
        assert kwargs["convert"] is False
        cal = pd.read_csv(kwargs["calibrant_csv"])
        ana = pd.read_csv(kwargs["analyte_csv"])
        assert ana.iloc[0]["mz"] == 202
        assert set(cal["raw_path"]) == {str((original / "raw/sample.raw").resolve())}
        assert ana.iloc[0]["raw_path"] == str((original / "raw/sample.raw").resolve())
        captured.append(kwargs)
        return {"converted_dir": str(source), "runs": [], "status": "complete"}
    monkeypatch.setattr(workflow, "run", fake_run)
    panel.calibrate_and_apply()
    _wait_until(lambda: not panel.is_busy)
    assert len(captured) == 1
    assert cal_path.read_bytes() == saved_cal
    assert not captured[0]["calibrant_csv"].exists()


def test_incomplete_pairs_and_missing_calibrants_fail_before_start(panel, tmp_path):
    source = tmp_path / "converted"
    source.mkdir()
    (source / "sample_ms.txt").write_text("100 1\n", encoding="utf-8")
    panel.pick_source.set_path(str(source))
    panel.pick_out.set_path(str(tmp_path / "output"))
    panel.process_and_view()
    assert not panel.is_busy
    assert panel._worker is None
    assert "Incomplete converted data" in panel.job_status.text()
    _pairs(source)
    panel.calibrate_and_apply()
    assert not panel.is_busy
    assert panel._worker is None
    assert "Calibrant table needs rows" in panel.job_status.text()
    assert panel._job_files is None


def test_error_restores_controls_and_exposes_completed_results(panel, tmp_path, monkeypatch):
    import deconvovo.waters_workflow as workflow
    source = _pairs(tmp_path / "converted")
    output = tmp_path / "output"
    output.mkdir()
    (output / "sample_2d_imms.html").write_text("<html>partial result</html>", encoding="utf-8")
    panel.pick_source.set_path(str(source))
    panel.pick_out.set_path(str(output))
    def fail(*_args, **_kwargs):
        raise RuntimeError("Reader failed: sample")
    monkeypatch.setattr(workflow, "run", fail)
    panel.process_and_view()
    _wait_until(lambda: not panel.is_busy)
    assert panel.btn_run.isEnabled()
    assert panel.btn_ccs.isEnabled()
    assert panel.sources_group.isEnabled()
    assert panel.btn_viewer.isEnabled()
    assert "Processing failed" in panel.job_status.text()
    assert "Reader failed: sample" in panel._log.text.toPlainText()
    assert panel._job_files is None


def test_viewer_opens_local_file_url(panel, tmp_path, monkeypatch):
    from gui.panels import waters_workspace
    source = _pairs(tmp_path / "converted")
    output = tmp_path / "output"
    output.mkdir()
    viewer = output / "sample_2d_imms.html"
    viewer.write_text("<html>offline</html>", encoding="utf-8")
    panel.pick_source.set_path(str(source))
    panel.pick_out.set_path(str(output))
    opened = []
    monkeypatch.setattr(waters_workspace.QDesktopServices, "openUrl", lambda url: opened.append(url.toLocalFile()))
    panel.open_viewer()
    assert [Path(path).resolve() for path in opened] == [viewer.resolve()]


def test_low_quality_fit_is_explicitly_marked_for_review(panel):
    panel._show_fit_status({"r2_lnln": 0.874117, "n_points": 25})
    assert "Review peak selection" in panel.calibration_panel.lbl_status.text()
    assert "0.87412" in panel.calibration_panel.lbl_status.text()


def test_discovery_excludes_other_output_runs_and_marks_empty_tables(tmp_path):
    source = _pairs(tmp_path / "source", "wanted")
    (source / "broken_ms.txt").write_text("", encoding="utf-8")
    (source / "broken_im.txt").write_text("100 1\n", encoding="utf-8")
    output = tmp_path / "output"
    _pairs(output / "_converted", "unrelated")
    (output / "unrelated_2d_imms.html").write_text("<html>old</html>", encoding="utf-8")
    rows = discover_runs(source, output)
    assert [row["run_name"] for row in rows] == ["broken", "wanted"]
    assert "ms" not in rows[0]
    assert "im" in rows[0]


def test_reload_existing_ccs_curve_and_summaries(panel, tmp_path):
    source = _pairs(tmp_path / "source")
    output = tmp_path / "output"
    ccs = output / "ccs_direct"
    ccs.mkdir(parents=True)
    calibration = {"r2_lnln": 0.999, "r2_linear_diagnostic": 0.999, "n_points": 3,
                   "X": 0.5, "A": 100, "gas_mw": 28.014,
                   "slope_diagnostic": 100, "intercept_diagnostic": 1}
    (ccs / "calibration_curve.json").write_text(json.dumps(calibration), encoding="utf-8")
    pd.DataFrame([{"name": str(i), "mw": 100 + i, "z": 1, "t_prime": 4 + i,
                   "t_double_prime": 2 + i, "CCS_literature": 200 + i}
                  for i in range(3)]).to_csv(ccs / "calibrant_summary.csv", index=False)
    pd.DataFrame([{"species": "target", "peak_ccs": 120, "above_threshold": True}]).to_csv(ccs / "analyte_summary.csv", index=False)
    panel.pick_source.set_path(str(source))
    panel.pick_out.set_path(str(output))
    panel.reload_runs()
    assert panel.calibration_panel.get_calibration()["n_points"] == 3
    assert len(panel.calibration_panel.summary_table.to_dataframe()) == 3
    assert len(panel.analyte_panel.summary_table.to_dataframe()) == 1
    assert "0.99900" in panel.calibration_panel.lbl_status.text()
    panel.combo_method.setCurrentText("twostep")
    assert panel.calibration_panel.get_calibration() is None
    assert panel.calibration_panel.summary_table.to_dataframe().empty
    assert panel.analyte_panel.summary_table.to_dataframe().empty


def test_console_preserves_selected_acquisition_and_reuses_real_renderer(panel, tmp_path, monkeypatch):
    source = _pairs(tmp_path / "source", "alpha")
    _pairs(source, "beta")
    calls = []
    monkeypatch.setattr(panel.viewer_widget, "set_run", lambda record, config, source: calls.append(record["run_name"]))
    panel.pick_source.set_path(str(source))
    panel.run_table.selectRow(1)
    assert panel.selected_run_label.toolTip() == "beta"
    panel.reload_runs()
    assert panel.run_table.currentRow() == 1
    assert panel.selected_run_label.toolTip() == "beta"
    assert calls[-1] == "beta"
    assert panel.console_splitter.count() == 3
    assert panel.run_table.isColumnHidden(1)
    panel.show_ccs_review(1)
    assert panel.tabs.currentIndex() == 1
    assert panel.ccs_editors.currentIndex() == 1


def test_selected_export_regenerates_current_defaults_in_background_and_preserves_workflow(panel, tmp_path, monkeypatch):
    import deconvovo.imms_html as html
    source = _pairs(tmp_path / "source", "chosen")
    output = tmp_path / "output"
    output.mkdir()
    destination = tmp_path / "portable.html"
    panel.pick_source.set_path(str(source))
    panel.pick_out.set_path(str(output))
    panel.viewer_panel.spin_noise2d.setValue(1.25)
    panel._last_result = {"status": "complete", "marker": "keep workflow"}
    entered = threading.Event()
    release = threading.Event()
    calls = []
    def plot(im_file, ms_file, name, directory, **kwargs):
        config = json.loads(Path(kwargs["config_path"]).read_text(encoding="utf-8"))
        calls.append((im_file, ms_file, name, config))
        entered.set()
        release.wait(2)
        (directory / f"{name}_2d_imms.html").write_text("<html><script>/* embedded Plotly */</script></html>", encoding="utf-8")
    monkeypatch.setattr(html, "plot_im_data", plot)
    assert panel.export_portable_html(destination) == destination
    try:
        _wait_until(entered.is_set)
        assert panel.is_busy
        assert not panel.btn_export.isEnabled()
        assert not panel.btn_run.isEnabled()
    finally:
        release.set()
    _wait_until(lambda: not panel.is_busy)
    assert len(calls) == 1
    assert calls[0][2] == "chosen"
    assert calls[0][3]["defaults"]["noise_2d"] == 1.25
    assert destination.read_text(encoding="utf-8").startswith("<html><script>")
    assert panel._last_result["marker"] == "keep workflow"
    assert panel.btn_export.isEnabled()
    assert panel._job_files is None
    assert panel._rows[0]["viewer"] == destination


def test_export_failure_preserves_existing_html_and_restores_controls(panel, tmp_path, monkeypatch):
    import deconvovo.imms_html as html
    source = _pairs(tmp_path / "source")
    destination = tmp_path / "portable.html"
    destination.write_text("previous valid portable viewer", encoding="utf-8")
    panel.pick_source.set_path(str(source))
    panel.pick_out.set_path(str(tmp_path / "output"))
    monkeypatch.setattr(html, "plot_im_data", lambda *_args, **_kwargs: (_ for _ in ()).throw(ValueError("bad selected data")))
    panel.export_portable_html(destination)
    _wait_until(lambda: not panel.is_busy)
    assert destination.read_text(encoding="utf-8") == "previous valid portable viewer"
    assert panel.btn_run.isEnabled()
    assert panel.btn_export.isEnabled()
    assert "Portable export failed" in panel.job_status.text()
    assert panel._job_files is None


def test_display_controls_debounce_live_renderer_config(panel, monkeypatch):
    calls = []
    monkeypatch.setattr(panel.viewer_widget, "apply_config", lambda config: calls.append(config))
    panel.viewer_panel.spin_noise2d.setValue(1.1)
    panel.viewer_panel.spin_noise2d.setValue(1.2)
    panel.viewer_panel.spin_noise2d.setValue(1.3)
    _wait_until(lambda: bool(calls))
    assert len(calls) == 1
    assert calls[0]["defaults"]["noise_2d"] == 1.3


@pytest.mark.parametrize('kind', ['ms', 'im'])
def test_export_cannot_replace_converted_source_tables(panel, tmp_path, kind):
    source = _pairs(tmp_path / 'converted')
    destination = source / f'sample_{kind}.txt'
    before = destination.read_bytes()
    panel.pick_source.set_path(str(source))
    panel.pick_out.set_path(str(tmp_path / 'output'))
    assert panel.export_portable_html(destination) is None
    assert destination.read_bytes() == before
    assert panel._worker is None
    assert panel._job_files is None
    assert not panel.is_busy
    assert 'source tables cannot be replaced' in panel.job_status.text()


@pytest.mark.parametrize('relative', ['_FUNC001.DAT', 'nested/viewer.html'])
def test_export_cannot_write_anywhere_inside_raw_acquisitions(panel, tmp_path, relative):
    source = _pairs(tmp_path / 'converted')
    raw = tmp_path / 'independent.RAW'
    destination = raw / relative
    destination.parent.mkdir(parents=True)
    destination.write_bytes(b'original acquisition bytes')
    panel.pick_source.set_path(str(source))
    panel.pick_out.set_path(str(tmp_path / 'output'))
    assert panel.export_portable_html(destination) is None
    assert destination.read_bytes() == b'original acquisition bytes'
    assert panel._worker is None
    assert panel._job_files is None
    assert not panel.is_busy
    assert 'outside RAW acquisition folders' in panel.job_status.text()


def test_export_rejects_other_file_formats_before_starting_worker(panel, tmp_path):
    source = _pairs(tmp_path / 'converted')
    destination = tmp_path / 'scientific_results.csv'
    destination.write_text('original,results\n', encoding='utf-8')
    panel.pick_source.set_path(str(source))
    panel.pick_out.set_path(str(tmp_path / 'output'))
    assert panel.export_portable_html(destination) is None
    assert destination.read_text(encoding="utf-8") == 'original,results\n'
    assert panel._worker is None
    assert panel._job_files is None
    assert not panel.is_busy
    assert '.html or .htm filename' in panel.job_status.text()


def test_export_resolves_aliases_before_checking_source_protection(panel, tmp_path):
    source = _pairs(tmp_path / 'converted')
    alias = tmp_path / 'dangerous.html'
    try:
        alias.symlink_to(source / 'sample_ms.txt')
    except OSError:
        pytest.skip('File symlinks are unavailable in this environment')
    panel.pick_source.set_path(str(source))
    panel.pick_out.set_path(str(tmp_path / 'output'))
    assert panel.export_portable_html(alias) is None
    assert alias.is_symlink()
    assert (source / 'sample_ms.txt').read_text(encoding="utf-8") == '100 1\n'
    assert panel._worker is None
    assert not panel.is_busy
