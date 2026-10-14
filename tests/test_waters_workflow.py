"""Failure, cache, and path contracts for the merged Waters workflow."""
from pathlib import Path
import json
import os

import pandas as pd
import pytest

from deconvovo import waters_convert, imms_plot, imms_ccs_calibrate, waters_workflow
from deconvovo.io import find_pusher_for_run
from deconvovo.paths import clean_path_text, local_path


def _raw(tmp_path):
    raw = tmp_path / "source" / "sample.raw"
    raw.mkdir(parents=True)
    (raw / "_FUNC001.DAT").write_bytes(b"unchanged source")
    return raw


def _mock_reader(monkeypatch, calls, status="complete"):
    monkeypatch.setattr(waters_convert, "setup_work_dir", lambda p: p)
    monkeypatch.setattr(waters_convert, "find_support_dlls", lambda p: [])

    def convert(raw, out, work, **settings):
        calls.append(settings)
        if status == "complete":
            (out / f"{raw.stem}_ms.txt").write_text("500 10\n")
            (out / f"{raw.stem}_im.txt").write_text("500 1 10\n")
        return {"run_name": raw.stem, "status": status}

    monkeypatch.setattr(waters_convert, "convert_one_raw", convert)


def test_reader_failure_stops_before_any_downstream_step(tmp_path, monkeypatch):
    raw = _raw(tmp_path)
    _mock_reader(monkeypatch, [], "error: missing im")
    downstream = []
    monkeypatch.setattr(imms_plot, "run", lambda *a, **k: downstream.append("viewer"))
    monkeypatch.setattr(imms_ccs_calibrate, "run", lambda *a, **k: downstream.append("ccs"))
    out = tmp_path / "out"
    with pytest.raises(RuntimeError, match="missing im"):
        waters_workflow.run(raw, out)
    assert not downstream
    assert not (out / "waters_workflow.json").exists()


def test_plotting_progress_works_with_windows_legacy_console_encoding(tmp_path, monkeypatch):
    import io
    import sys

    source = tmp_path / "converted"
    source.mkdir()
    (source / "sample_ms.txt").write_text("500 10\n")
    (source / "sample_im.txt").write_text("500 1 10\n")
    monkeypatch.setattr(imms_plot, "find_pusher_for_run", lambda *args: 110.0)
    monkeypatch.setattr(imms_plot, "_plot_one_run", lambda args: {"run_name": "sample", "status": ["IM"]})
    buffer = io.BytesIO()
    console = io.TextIOWrapper(buffer, encoding="cp1252")
    monkeypatch.setattr(sys, "stdout", console)
    result = imms_plot.run(source, tmp_path / "html", n_workers=1)
    console.flush()
    assert result[0]["status"] == ["IM"]
    assert "Pusher: 110.00 us" in buffer.getvalue().decode("cp1252")


def test_parent_batch_processes_direct_and_nested_acquisitions(tmp_path, monkeypatch):
    source = tmp_path / "batch"
    (source / "direct.raw").mkdir(parents=True)
    (source / "condition/nested.raw").mkdir(parents=True)
    _mock_reader(monkeypatch, [])
    result = waters_workflow.run(source, tmp_path / "out", html=False)
    assert {r["run_name"] for r in result["runs"]} == {"direct", "nested"}


def test_output_cannot_be_inside_any_original_acquisition(tmp_path, monkeypatch):
    raw = _raw(tmp_path)
    _mock_reader(monkeypatch, [])
    with pytest.raises(ValueError, match="outside the RAW"):
        waters_workflow.run(raw.parent, raw / "results", html=False)
    assert not (raw / "results").exists()


def test_ccs_only_output_cannot_create_anything_inside_raw(tmp_path):
    raw = _raw(tmp_path)
    output = raw / "new folder/ccs"
    with pytest.raises(ValueError, match="outside the RAW"):
        imms_ccs_calibrate.run(output, tmp_path, tmp_path / "unused.csv")
    assert not (raw / "new folder").exists()


def test_batch_names_cannot_collide_on_windows(tmp_path, monkeypatch):
    source = tmp_path / "batch"
    (source / "one/sample.raw").mkdir(parents=True)
    (source / "two/Sample.raw").mkdir(parents=True)
    with pytest.raises(ValueError, match="distinct run names"):
        waters_workflow.run(source, tmp_path / "out", html=False)
    assert not (tmp_path / "out").exists()


def test_conversion_reuse_requires_matching_parameters_and_source(tmp_path, monkeypatch):
    raw = _raw(tmp_path)
    calls = []
    _mock_reader(monkeypatch, calls)
    out = tmp_path / "out"
    waters_workflow.run(raw, out, html=False)
    waters_workflow.run(raw, out, html=False, skip_existing=True)
    assert len(calls) == 1
    waters_workflow.run(raw, out, html=False, skip_existing=True, func=2)
    assert len(calls) == 2
    (raw / "_FUNC001.DAT").write_bytes(b"changed acquisition")
    waters_workflow.run(raw, out, html=False, skip_existing=True, func=2)
    assert len(calls) == 3


def test_failed_refresh_invalidates_previous_conversion_cache(tmp_path, monkeypatch):
    raw = _raw(tmp_path)
    calls = []
    _mock_reader(monkeypatch, calls)
    out = tmp_path / "out"
    waters_workflow.run(raw, out, html=False)
    _mock_reader(monkeypatch, calls, "error: missing im")
    with pytest.raises(RuntimeError):
        waters_workflow.run(raw, out, html=False)
    metadata = json.loads((out / "_converted/conversion_settings.json").read_text())
    assert raw.stem not in metadata
    assert not (out / "waters_workflow.json").exists()


def test_existing_converted_pairs_refresh_viewer_and_report_progress(tmp_path, monkeypatch):
    source = tmp_path / "converted"
    source.mkdir()
    for name in ["first", "second"]:
        (source / f"{name}_ms.txt").write_text("500 10\n")
        (source / f"{name}_im.txt").write_text("500 1 10\n")
    calls = []

    def plot(paths, output, **kwargs):
        calls.append(kwargs)
        name = Path(paths[0]).name.removesuffix("_ms.txt")
        (output / f"{name}_2d_imms.html").write_text("<html>fresh</html>")
        return [{"run_name": name, "status": ["IM"]}]

    monkeypatch.setattr(imms_plot, "run", plot)
    events = []
    result = waters_workflow.run(source, tmp_path / "out", skip_existing=True,
                                progress=lambda *event: events.append(event))
    assert len(result["viewers"]) == 2
    assert all(c["skip_existing"] is False for c in calls)
    assert events[-1] == (2, 2, "Workflow complete")


def test_viewer_error_is_not_reported_as_complete(tmp_path, monkeypatch):
    source = tmp_path / "converted"
    source.mkdir()
    for kind in ["ms", "im"]:
        (source / f"sample_{kind}.txt").write_text("500 10\n")
    monkeypatch.setattr(imms_plot, "run", lambda *a, **k: [{"status": ["ERROR: bad table"]}])
    with pytest.raises(RuntimeError, match="Viewer generation failed"):
        waters_workflow.run(source, tmp_path / "out")


def test_csv_paths_remain_portable_and_old_cwd_paths_work(tmp_path, monkeypatch):
    folder = tmp_path / "demo"
    raw = folder / "raw/sample.raw"
    raw.mkdir(parents=True)
    csv = folder / "calibrants.csv"
    pd.DataFrame({"raw_path": ["raw/sample.raw"]}).to_csv(csv, index=False)
    monkeypatch.chdir(tmp_path)
    assert imms_ccs_calibrate._read_species_csv(csv)["raw_path"][0] == str(raw)
    pd.DataFrame({"raw_path": ["demo/raw/sample.raw"]}).to_csv(csv, index=False)
    assert imms_ccs_calibrate._read_species_csv(csv)["raw_path"][0] == str(raw)


def test_csv_windows_relative_paths_quotes_spaces_and_utf8_bom(tmp_path, monkeypatch):
    folder = tmp_path / "portable demo"
    raw = folder / "raw data/sample.raw"
    raw.mkdir(parents=True)
    csv = folder / "calibrants.csv"
    pd.DataFrame({"raw_path": ['"raw data\\sample.raw"']}).to_csv(csv, index=False, encoding="utf-8-sig")
    monkeypatch.chdir(tmp_path)
    frame = imms_ccs_calibrate._read_species_csv(csv)
    assert Path(frame["raw_path"][0]) == raw
    assert not (tmp_path / "raw data\\sample.raw").exists()


def test_pasted_paths_expand_quotes_home_and_environment(tmp_path, monkeypatch):
    monkeypatch.setenv("WATERS_PATH_TEST", str(tmp_path / "real data"))
    variable = "%WATERS_PATH_TEST%" if os.name == "nt" else "$WATERS_PATH_TEST"
    assert local_path(f' "{variable}" ') == tmp_path / "real data"
    assert local_path('"~"') == Path.home()
    assert clean_path_text("  ") == ""


@pytest.mark.parametrize("text", [r'C:\Data With Spaces\sample.raw', r'\\server\share\Data With Spaces\sample.raw'])
def test_windows_absolute_paths_never_become_fake_posix_paths(text):
    assert clean_path_text(f'"{text}"') == text
    if os.name == "nt":
        assert local_path(f'"{text}"').is_absolute()
        assert local_path(text).drive
    else:
        with pytest.raises(ValueError, match="needs Windows"):
            local_path(text)


def test_workflow_quoted_paths_preserve_cache_and_original_raw(tmp_path, monkeypatch):
    raw = tmp_path / "Raw Data/sample.raw"
    raw.mkdir(parents=True)
    original = raw / "_FUNC001.DAT"
    original.write_bytes(b"untouched RAW")
    calls = []
    _mock_reader(monkeypatch, calls)
    out = tmp_path / "Output Data"
    config = tmp_path / "Viewer Config.json"
    config.write_text("{}", encoding="utf-8")
    for _ in range(2):
        waters_workflow.run(f'"{raw}"', f'"{out}"', html=False,
                            config_path=f'"{config}"', skip_existing=True)
    assert len(calls) == 1
    metadata = json.loads((out / "_converted/conversion_settings.json").read_text(encoding="utf-8"))
    assert Path(metadata[raw.stem]["raw_source"]) == raw
    assert original.read_bytes() == b"untouched RAW"
    with pytest.raises(ValueError, match="outside the RAW"):
        waters_workflow.run(f'"{raw}"', f'"{raw / "results"}"', html=False)
    assert not (raw / "results").exists()


def test_csv_empty_raw_path_is_reported_before_calibration(tmp_path):
    csv = tmp_path / "calibrants.csv"
    pd.DataFrame({"raw_path": [None]}).to_csv(csv, index=False)
    with pytest.raises(ValueError, match="Empty raw_path"):
        imms_ccs_calibrate._read_species_csv(csv)


def test_mixed_edc_coefficients_reject_calibration(tmp_path, monkeypatch):
    raws = [tmp_path / "one.raw", tmp_path / "two.raw"]
    for raw in raws:
        raw.mkdir()
    csv = tmp_path / "calibrants.csv"
    pd.DataFrame({"raw_path": [str(p) for p in raws]}).to_csv(csv, index=False)
    monkeypatch.setattr(imms_ccs_calibrate, "find_pusher_period", lambda p: 100)
    monkeypatch.setattr(imms_ccs_calibrate, "read_extern_inf", lambda p: {
        "EDC Delay Coefficient": "1.58" if p.name == "one.raw" else "1.72"
    })
    with pytest.raises(ValueError, match="different EDC"):
        imms_ccs_calibrate.run(tmp_path / "out", tmp_path, csv)


def test_reopened_converted_output_recovers_drift_time_metadata(tmp_path, monkeypatch):
    raw = _raw(tmp_path)
    converted = tmp_path / "out/_converted"
    converted.mkdir(parents=True)
    (converted / "conversion_settings.json").write_text(json.dumps({
        "sample": {"raw_source": str(raw)}
    }))
    monkeypatch.setattr("deconvovo.io.find_pusher_period", lambda p: 110 if p == raw else None)
    assert find_pusher_for_run(converted, "sample") == 110


def test_ccs_refresh_removes_obsolete_outputs_without_destroying_prior_results(tmp_path, monkeypatch):
    out = tmp_path / "ccs"
    out.mkdir()
    (out / "analyte_summary.csv").write_text("old analytes")
    (out / "cali_png").mkdir()
    (out / "cali_png/obsolete.png").write_bytes(b"old plot")

    def fresh(staged, *args):
        staged.mkdir()
        (staged / "calibration_curve.json").write_text("{}")
        return {"r2_lnln": 0.9}

    monkeypatch.setattr(imms_ccs_calibrate, "_run_calibration", fresh)
    imms_ccs_calibrate.run(out, tmp_path, tmp_path / "unused.csv")
    assert (out / "calibration_curve.json").exists()
    assert not (out / "analyte_summary.csv").exists()
    assert not (out / "cali_png/obsolete.png").exists()
    assert any(tmp_path.glob(".ccs-previous-*/ccs/analyte_summary.csv"))


def test_failed_ccs_refresh_keeps_prior_complete_outputs(tmp_path, monkeypatch):
    out = tmp_path / "ccs"
    out.mkdir()
    (out / "calibration_curve.json").write_text("old result")

    def fail(staged, *args):
        staged.mkdir()
        (staged / "calibration_curve.json").write_text("partial")
        raise ValueError("not enough peaks")

    monkeypatch.setattr(imms_ccs_calibrate, "_run_calibration", fail)
    with pytest.raises(ValueError, match="not enough peaks"):
        imms_ccs_calibrate.run(out, tmp_path, tmp_path / "unused.csv")
    assert (out / "calibration_curve.json").read_text() == "old result"
    assert not list(tmp_path.glob(".waters-ccs-*"))
