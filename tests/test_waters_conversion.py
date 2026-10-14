from __future__ import annotations

import csv
from pathlib import Path
import subprocess
import sys

import pytest

from deconvovo import cli, imms_convert, waters_convert
from deconvovo.taskqueue import TaskQueue
from portable_waters_convert import waters_convert as portable_converter


@pytest.fixture(params=[waters_convert, portable_converter], ids=["installed", "portable"])
def converter(request, monkeypatch):
    module = request.param
    monkeypatch.setattr(module, "IS_WINDOWS", True)
    return module


def _paths(tmp_path):
    raw = tmp_path / "sample.raw"
    raw.mkdir()
    output = tmp_path / "converted"
    output.mkdir()
    return raw, output, tmp_path / "reader"


def _reader_result(command, *, outputs=("ms", "im"), empty=(), returncode=0):
    for name in outputs:
        flag = "-m" if name == "ms" else "-i"
        path = Path(command[command.index(flag) + 1])
        path.write_text("" if name in empty else f"{name} data\n")
    return subprocess.CompletedProcess(command, returncode, "reader detail", "")


@pytest.mark.parametrize("empty", [(), ("im",)])
def test_zero_exit_missing_or_empty_output_is_failure(converter, monkeypatch, tmp_path, empty):
    raw, output, work = _paths(tmp_path)
    # Existing data must not disguise a reader that writes no valid IM output.
    for name in ("ms", "im"):
        (output / f"sample_{name}.txt").write_text("old data")
    writes = ("ms", "im") if empty else ("ms",)
    monkeypatch.setattr(converter.subprocess, "run", lambda command, **kwargs:
                        _reader_result(command, outputs=writes, empty=empty))
    result = converter.convert_one_raw(raw, output, work)
    assert result["status"] == "error: missing or empty im output"
    assert result["stderr"] == "reader detail"
    assert (output / "sample_ms.txt").read_text() == "old data"
    assert (output / "sample_im.txt").read_text() == "old data"
    assert not list(output.glob(".waters-*"))


def test_reader_failure_preserves_existing_outputs(converter, monkeypatch, tmp_path):
    raw, output, work = _paths(tmp_path)
    (output / "sample_ms.txt").write_text("old data")
    monkeypatch.setattr(converter.subprocess, "run", lambda command, **kwargs:
                        _reader_result(command, returncode=5))
    result = converter.convert_one_raw(raw, output, work)
    assert result["status"] == "error:5"
    assert (output / "sample_ms.txt").read_text() == "old data"
    assert not (output / "sample_im.txt").exists()


def test_valid_output_promoted_without_ms_smoothing(converter, monkeypatch, tmp_path):
    raw, output, work = _paths(tmp_path)
    commands = []

    def reader(command, **kwargs):
        commands.append(command)
        return _reader_result(command)

    monkeypatch.setattr(converter.subprocess, "run", reader)
    result = converter.convert_one_raw(raw, output, work, func=2)
    assert result["status"] == "complete"
    assert Path(result["ms_file"]).read_text() == "ms data\n"
    assert Path(result["im_file"]).read_text() == "im data\n"
    command = commands[0]
    assert command[command.index("--fn") + 1] == "2"
    assert command[command.index("--ms_smooth_window") + 1] == "0"
    assert command[command.index("--ms_number_smooth") + 1] == "0"


def test_only_requested_output_is_required(converter, monkeypatch, tmp_path):
    raw, output, work = _paths(tmp_path)

    def reader(command, **kwargs):
        assert "-i" not in command
        assert command[command.index("--skip_im") + 1] == "1"
        return _reader_result(command, outputs=("ms",))

    monkeypatch.setattr(converter.subprocess, "run", reader)
    result = converter.convert_one_raw(raw, output, work, skip_im=True)
    assert result["status"] == "complete"
    assert "im_file" not in result


def test_timeout_removes_staging_and_preserves_outputs(converter, monkeypatch, tmp_path):
    raw, output, work = _paths(tmp_path)
    (output / "sample_ms.txt").write_text("old data")

    def reader(command, **kwargs):
        _reader_result(command, outputs=("ms",))
        raise subprocess.TimeoutExpired(command, 120)

    monkeypatch.setattr(converter.subprocess, "run", reader)
    with pytest.raises(subprocess.TimeoutExpired):
        converter.convert_one_raw(raw, output, work)
    assert (output / "sample_ms.txt").read_text() == "old data"
    assert not list(output.glob(".waters-*"))


def _configure_main(converter, monkeypatch, raw, output, flags):
    monkeypatch.setattr(sys, "argv", ["waters_convert", str(raw), "-o", str(output), *flags])
    monkeypatch.setattr(converter, "find_cdcreader", lambda: output / "reader" / "CDCReader.exe")
    if converter is waters_convert:
        monkeypatch.setattr(converter, "setup_work_dir", lambda path: output / "reader")
    else:
        monkeypatch.setattr(converter, "find_support_dlls", lambda path: [])


def test_cli_skip_existing_honors_output_selection(converter, monkeypatch, tmp_path):
    raw, output, work = _paths(tmp_path)
    (output / "sample_ms.txt").write_text("data")
    _configure_main(converter, monkeypatch, raw, output, ["--skip-existing", "--skip-im"])

    def unexpected(*args, **kwargs):
        pytest.fail("Complete requested output must be reused")

    monkeypatch.setattr(converter, "convert_one_raw", unexpected)
    converter.main()
    with (output / "convert_summary.csv").open() as handle:
        assert list(csv.DictReader(handle))[0]["status"] == "skipped"


def test_cli_retries_empty_cache_and_exits_nonzero_on_failure(converter, monkeypatch, tmp_path):
    raw, output, work = _paths(tmp_path)
    (output / "sample_ms.txt").touch()
    (output / "sample_im.txt").write_text("data")
    _configure_main(converter, monkeypatch, raw, output, ["--skip-existing"])
    monkeypatch.setattr(converter, "convert_one_raw", lambda *args, **kwargs:
                        {"run_name": "sample", "status": "error: missing ms output"})
    with pytest.raises(SystemExit) as error:
        converter.main()
    assert error.value.code == 1
    with (output / "convert_summary.csv").open() as handle:
        assert list(csv.DictReader(handle))[0]["status"].startswith("error")


def test_legacy_conversion_requires_both_cached_outputs(monkeypatch, tmp_path):
    raw, output, work = _paths(tmp_path)
    (output / "sample_ms.txt").write_text("old data")
    monkeypatch.setattr(waters_convert, "convert_one_raw", lambda *args:
                        {"run_name": "sample", "status": "error: missing im output"})
    result = imms_convert._convert_one({"raw_dir": str(raw), "out_dir": str(output),
                                      "work_dir": str(work)})
    assert result["status"].startswith("error")
    monkeypatch.setattr(waters_convert, "setup_work_dir", lambda path: work)
    with pytest.raises(RuntimeError, match="sample: error: missing im output"):
        imms_convert.run(raw, output, n_workers=1)


def test_legacy_cli_records_converter_failure(monkeypatch, tmp_path):
    raw, output, work = _paths(tmp_path)
    monkeypatch.setattr(waters_convert, "convert_one_raw", lambda *args:
                        {"status": "error: missing im output", "stderr": "reader detail"})
    result = cli._convert_one({"run_name": "sample", "raw_dir": str(raw),
                               "out_dir": str(output), "work_dir": str(work)})
    assert result["error"] == "error: missing im output: reader detail"
    assert result["status"][0].startswith("ERROR")


def _task(args):
    if args.get("fail"):
        return {"status": "error: reader wrote no data"}
    Path(args["marker"]).touch()
    return {"status": ["OK"]}


@pytest.mark.parametrize("workers", [1, 2])
def test_failed_dependencies_skip_descendants_and_keep_independent_tasks(tmp_path, workers):
    queue = TaskQueue(n_workers=workers)
    first = queue.add(_task, {"fail": True})
    second = queue.add(_task, {"marker": str(tmp_path / "dependent")}, depends_on=[first])
    third = queue.add(_task, {"marker": str(tmp_path / "descendant")}, depends_on=[second])
    independent = queue.add(_task, {"marker": str(tmp_path / "independent")})
    results = queue.run()
    assert "dependencies failed" in results[second]["error"]
    assert "dependencies failed" in results[third]["error"]
    assert results[independent]["status"] == ["OK"]
    assert not (tmp_path / "dependent").exists()
    assert not (tmp_path / "descendant").exists()
    assert (tmp_path / "independent").exists()
