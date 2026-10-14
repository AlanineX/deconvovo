"""Standalone offline HTML remains an optional artifact, outside the native GUI."""
from pathlib import Path
import os
import re

import numpy as np
from deconvovo.imms_html import plot_im_data

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")


def test_standalone_html_keeps_embedded_runtime_and_complete_controls(tmp_path):
    im = tmp_path / "sample_im.txt"
    np.savetxt(im, [[100, 0, 1], [101, 0, 2], [102, 1, 3], [103, 1, 4]])
    plot_im_data(im, None, "sample", tmp_path, pusher_us=100)
    html = (tmp_path / "sample_2d_imms.html").read_text(encoding="utf-8")
    assert "<script src=" not in html
    assert "Plotly.newPlot" in html
    for control in ("sm2d", "mzscale", "noise2d", "denoise", "smdrift", "noisedrift", "msmode",
                    "noisemz", "nbins", "dbins", "cmap", "render", "zscale"):
        assert f'id="{control}"' in html
    assert "watersExportCurrent" in html
    assert "waters-console-style" not in html
    assert (tmp_path / "sample_2d_imms.csv").is_file()


def test_standalone_sg_smoothing_preserves_constant_signal_without_extra_gain(tmp_path):
    from PySide6.QtWidgets import QApplication
    from PySide6.QtQml import QJSEngine
    app = QApplication.instance() or QApplication([])
    im = tmp_path / "sample_im.txt"
    np.savetxt(im, [[100, 0, 1], [101, 0, 2], [102, 1, 3]])
    plot_im_data(im, None, "sample", tmp_path)
    html = (tmp_path / "sample_2d_imms.html").read_text(encoding="utf-8")
    function = re.search(r"function sgSmooth1D\(y,w,ord\)\{.*?\n\}\nfunction sgSmooth2D", html, re.S).group(0)
    function = function.rsplit("\nfunction sgSmooth2D", 1)[0]
    engine = QJSEngine()
    assert not engine.evaluate(function).isError()
    actual = engine.evaluate("JSON.stringify(sgSmooth1D([7,7,7,7,7,7,7],5,2))").toString()
    import json
    np.testing.assert_allclose(json.loads(actual), np.full(7, 7), atol=1e-12)
    app.processEvents()


def test_native_module_contains_no_browser_or_javascript_runtime():
    root = Path(__file__).parents[1]
    for path in (root / "gui/widgets/waters_viewer.py", root / "portable_waters_convert/gui/widgets/waters_viewer.py"):
        source = path.read_text(encoding="utf-8")
        assert "QtWebEngine" not in source and "runJavaScript" not in source and "QJSEngine" not in source
