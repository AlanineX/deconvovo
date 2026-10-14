"""Real Waters data previews and offline enhanced-viewer integration."""
from __future__ import annotations

import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("DECONVOVO_DISABLE_WEBENGINE", "1")

import numpy as np
import pytest
from PySide6.QtCore import QUrl
from PySide6.QtWidgets import QApplication, QWidget

import gui.widgets.waters_viewer as viewer


@pytest.fixture
def widget():
    app = QApplication.instance() or QApplication([])
    panel = viewer.WatersViewerWidget()
    yield panel
    panel.wait_for_workers()
    panel.close()
    panel.deleteLater()
    app.processEvents()


def wait_until(predicate, timeout=4):
    deadline = time.monotonic() + timeout
    while not predicate() and time.monotonic() < deadline:
        QApplication.instance().processEvents()
        time.sleep(0.01)
    assert predicate()


def text_record(tmp_path, name="sample"):
    ms, im = tmp_path / f"{name}_ms.txt", tmp_path / f"{name}_im.txt"
    ms.write_text("100 1\n101 5\n102 7\n103 3\n")
    im.write_text("100 0 1\n101 0 2\n102 1 3\n103 1 4\n100 2 2\n103 2 1\n")
    return {"run_name": name, "ms": ms, "im": im}


def test_raw_input_has_honest_empty_state_and_no_browser(widget, tmp_path):
    raw = tmp_path / "sample.raw"
    raw.mkdir()
    widget.set_run({"run_name": "sample", "raw": raw}, {})
    assert widget.web_view is None
    assert widget.im_stack.currentWidget() is widget.empty
    assert "Process & view" in widget.empty.text()
    assert not widget.btn_browser.isEnabled()
    assert [widget.tabs.tabText(i) for i in range(3)] == ["IM-MS", "Mass profile", "Mobiligram"]


def test_pending_plot_draw_is_cancelled_when_qt_destroys_canvas(widget):
    import shiboken6
    from gui.widgets.plot_canvas import PlotCanvas

    plot = PlotCanvas()
    canvas = plot.canvas
    canvas.draw_idle()
    assert canvas._draw_pending
    shiboken6.delete(canvas)
    assert not canvas._draw_pending
    # This is the callback already queued by draw_idle. It must return without
    # calling width/height on the deleted C++ widget.
    canvas._draw_idle()
    plot.close()


def test_preview_uses_raw_ms_and_existing_drift_extractor(widget, tmp_path):
    record = text_record(tmp_path)
    widget.set_run(record, {})
    wait_until(lambda: len(widget.mass_plot.ax().lines) == 1)
    np.testing.assert_array_equal(widget.mass_plot.ax().lines[0].get_ydata(), [1, 5, 7, 3])
    np.testing.assert_array_equal(widget.mobility_plot.ax().lines[0].get_ydata(), [3, 7, 3])
    assert widget.im_stack.currentWidget() is widget.overview_plot
    assert "drift bins" in widget.status.text()
    assert "Raw drift bin" == widget.mobility_plot.ax().get_xlabel()


def test_window_preview_reuses_extractor_and_preserves_selected_tab(widget, tmp_path):
    widget.tabs.setCurrentIndex(2)
    widget.set_run(text_record(tmp_path), {"initial_view": {"mz_range": [101, 102]}})
    wait_until(lambda: len(widget.mobility_plot.ax().lines) == 1)
    np.testing.assert_array_equal(widget.mobility_plot.ax().lines[0].get_ydata(), [2, 3, 0])
    assert widget.tabs.currentIndex() == 2


def test_stale_read_cannot_replace_current_selection(widget, tmp_path):
    widget.set_run({"run_name": "new", "raw": tmp_path}, {})
    before = widget.status.text()
    widget._preview_loaded({"generation": widget._generation - 1, "ms": (np.arange(2), np.ones(2)), "errors": []})
    assert widget.status.text() == before
    assert not widget.mass_plot.ax().lines
    assert widget._record["run_name"] == "new"


def test_empty_or_invalid_txt_is_reported_without_invented_plots(widget, tmp_path):
    record = text_record(tmp_path)
    record["ms"].write_text("invalid")
    record["im"].write_text("100 0 1\n")
    widget.set_run(record, {})
    wait_until(lambda: "no readable spectrum" in widget.empty.text().casefold())
    assert not widget.mass_plot.ax().lines
    assert "Mass spectrum:" in widget.status.text()


def test_browser_and_export_use_exact_existing_local_html(widget, tmp_path, monkeypatch):
    html = tmp_path / "real_2d_imms.html"
    html.write_text("<html><script>offline=true</script></html>")
    urls = []
    monkeypatch.setattr(viewer.QDesktopServices, "openUrl", lambda url: urls.append(url))
    widget.set_run({"run_name": "real", "viewer": html}, {})
    widget.open_browser()
    assert Path(urls[0].toLocalFile()).resolve() == html.resolve()
    target = tmp_path / "standalone.html"
    assert widget.export_html(target) == target
    assert target.read_bytes() == html.read_bytes()
    assert widget.btn_browser.isEnabled()


def test_html_copy_cannot_replace_source_data_or_enter_raw_folder(widget, tmp_path):
    record = text_record(tmp_path)
    html = tmp_path / "view.html"
    html.write_text("<html>offline</html>")
    record["viewer"] = html
    widget.set_run(record, {})
    source = record["ms"].read_bytes()
    with pytest.raises(ValueError, match="source MS"):
        widget.export_html(record["ms"])
    assert record["ms"].read_bytes() == source
    raw = tmp_path / "original.RAW"
    raw.mkdir()
    with pytest.raises(ValueError, match="outside RAW"):
        widget.export_html(raw / "view.html")
    assert not (raw / "view.html").exists()


def test_native_download_rejects_source_path_before_accept(widget, tmp_path, monkeypatch):
    record = text_record(tmp_path)
    widget.set_run(record, {})
    before = record["im"].read_bytes()
    monkeypatch.setattr(viewer.QFileDialog, "getSaveFileName", lambda *a: (str(record["im"]), ""))
    class Download:
        accepted = False
        cancelled = False
        def downloadFileName(self): return "mobiligram.csv"
        def cancel(self): self.cancelled = True
        def accept(self): self.accepted = True
        def setDownloadDirectory(self, value): pytest.fail("Protected download reached filesystem setup")
        def setDownloadFileName(self, value): pytest.fail("Protected download reached filesystem setup")
    download = Download()
    widget._download_file(download)
    assert download.cancelled and not download.accepted
    assert record["im"].read_bytes() == before


def test_native_download_keeps_correct_suffix_and_writable_output(widget, tmp_path, monkeypatch):
    destination = tmp_path / "plot"
    monkeypatch.setattr(viewer.QFileDialog, "getSaveFileName", lambda *a: (str(destination), ""))
    class Download:
        accepted = False
        def downloadFileName(self): return "plot.png"
        def setDownloadDirectory(self, value): self.directory = value
        def setDownloadFileName(self, value): self.name = value
        def cancel(self): pytest.fail("Valid PNG export was cancelled")
        def accept(self): self.accepted = True
    download = Download()
    widget._download_file(download)
    assert download.accepted
    assert download.directory == str(tmp_path)
    assert download.name == "plot.png"


def test_large_html_is_loaded_by_url_without_reading_or_set_html(widget, tmp_path):
    class FakePage:
        def runJavaScript(self, script):
            pass
    class FakeView(QWidget):
        def __init__(self):
            super().__init__()
            self.urls = []
        def stop(self):
            pass
        def load(self, url):
            self.urls.append(url)
        def page(self):
            return FakePage()
        def setHtml(self, _html):
            pytest.fail("Multi-megabyte Waters HTML must be loaded from a file URL")
    widget.web_view = FakeView()
    widget.im_stack.addWidget(widget.web_view)
    html = tmp_path / "large_2d_imms.html"
    with html.open("wb") as stream:
        stream.seek(37_000_000)
        stream.write(b" ")
    widget.set_run({"run_name": "large", "viewer": html}, {})
    assert widget.web_view.urls[0].isLocalFile()
    assert Path(widget.web_view.urls[0].toLocalFile()).resolve() == html.resolve()


def test_loaded_file_url_is_normalized_and_stale_or_remote_completions_are_ignored(widget, tmp_path):
    class FakePage:
        def __init__(self):
            self.calls = []
        def runJavaScript(self, script, callback=None):
            self.calls.append(script)
            if callback:
                callback(True)
    class FakeView(QWidget):
        def __init__(self):
            super().__init__()
            self.current_url = QUrl()
            self.document = FakePage()
        def stop(self):
            pass
        def load(self, url):
            self.current_url = url
        def url(self):
            return self.current_url
        def page(self):
            return self.document
    widget.web_view = FakeView()
    widget.im_stack.addWidget(widget.web_view)
    html = tmp_path / "selected viewer.html"
    html.write_text("<html>offline</html>", encoding="utf-8")
    widget.set_run({"run_name": "selected", "viewer": html}, {})
    # Qt uses portable separators on Windows. Normalize paths rather than
    # comparing Qt's path string to the native pathlib spelling.
    widget.web_view.current_url = QUrl.fromLocalFile(str(tmp_path / "sub" / ".." / html.name))
    widget._html_loaded(True)
    assert widget._html_ready
    assert len(widget.web_view.document.calls) == 2
    widget.clear("A newer selection is pending")
    widget._viewer_path = (tmp_path / "new.html").resolve()
    widget._html_generation = widget._generation
    widget._html_loaded(True)
    assert not widget._html_ready
    assert len(widget.web_view.document.calls) == 2
    widget.web_view.current_url = QUrl("https://example.invalid/new.html")
    widget._html_loaded(True)
    assert not widget._html_ready


@pytest.mark.parametrize("url,allowed", [("file:///tmp/view.html", True), ("data:text/plain,ok", True),
    ("blob:file:///tmp/id", True), ("qrc:/qtwebchannel/qwebchannel.js", True),
    ("about:blank", True), ("https://cdn.plot.ly/plotly.js", False),
    ("http://localhost:8000/data", False), ("wss://example.org", False)])
def test_offline_policy_blocks_remote_requests(url, allowed):
    assert viewer.is_local_resource(QUrl(url)) is allowed


def test_live_config_delegates_to_original_viewer_algorithms():
    config = {"defaults": {"smooth_2d": {"method": "gaussian", "sigma": 0.1234},
        "mz_bins": 400, "noise_2d": 3, "ms_mode": "profile"}}
    script = viewer._config_script(config)
    assert json.dumps(config["defaults"]) in script
    assert "rebinFull()" in script
    assert "applyAll()" in script
    assert "gaussSmooth" not in script
    assert not re.search(r"\bHM_MZ_LO0\s*=(?!=)", script)


def test_disable_webengine_overrides_explicit_headless_opt_in(monkeypatch):
    monkeypatch.setenv("DECONVOVO_DISABLE_WEBENGINE", "1")
    monkeypatch.setenv("DECONVOVO_ENABLE_HEADLESS_WEBENGINE", "1")
    assert viewer._webengine_allowed() is False


def test_real_webengine_loads_offline_and_live_controls_reuse_original_html(tmp_path):
    from deconvovo.imms_html import plot_im_data

    record = text_record(tmp_path)
    config = json.loads((Path(__file__).parents[1] / "config/imms_plot_config.json").read_text(encoding="utf-8"))
    config["defaults"].update({"mz_bins": 200, "drift_bins": "native", "ms_mode": "profile", "noise_2d": 0})
    config["initial_view"] = {"mz_range": [100, 103]}
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(config))
    plot_im_data(record["im"], record["ms"], "sample", tmp_path, config_path=config_path)
    html = tmp_path / "sample_2d_imms.html"
    assert "<script src=" not in html.read_text(encoding="utf-8")
    script = r'''
import json,sys,time
from pathlib import Path
from PySide6.QtCore import QUrl,QCoreApplication,QEvent
from PySide6.QtWidgets import QApplication
from gui.widgets.waters_viewer import WatersViewerWidget,WEBENGINE_AVAILABLE
app=QApplication([])
assert WEBENGINE_AVAILABLE
w=WatersViewerWidget();w.resize(550,400);w.show()
cfg=json.loads(Path(sys.argv[2]).read_text(encoding="utf-8"))
load_events=[]
w.web_view.loadFinished.connect(lambda ok:load_events.append({'ok':ok,'url':w.web_view.url().toString()}))
w.set_run({'run_name':'sample','viewer':Path(sys.argv[1])},cfg)
def js(expression):
    result=[];w.web_view.page().runJavaScript(expression,result.append)
    deadline=time.monotonic()+5
    while not result and time.monotonic()<deadline:
        app.processEvents();time.sleep(.01)
    return result[0] if result else None
deadline=time.monotonic()+25
plot_ready=False
while time.monotonic()<deadline:
    app.processEvents();time.sleep(.02)
    plot_ready=js("typeof P!=='undefined' && !!P._fullLayout")
    if w._html_ready and plot_ready:break
else:
    diagnostics={'load_events':load_events,'html_ready':w._html_ready,'plot_ready':plot_ready,
        'loaded_file':w.web_view.url().toLocalFile(),'expected_file':str(w._viewer_path),
        'blocked_requests':w.blocked_requests,'platform':app.platformName(),
        'document':js("JSON.stringify({state:document.readyState,plotly:typeof Plotly,plot:typeof P})")}
    raise AssertionError('Offline viewer readiness failed: '+json.dumps(diagnostics))
assert w.web_view.url().isLocalFile()
assert w.blocked_requests==0
assert not w.status.isVisible()
assert 'Offline' in w.tabs.tabToolTip(0)
cfg['defaults'].update({'mz_bins':400,'noise_2d':3,'ms_mode':'stick','colormap':'Plasma'})
w.apply_config(cfg)
deadline=time.monotonic()+10
while time.monotonic()<deadline:
    app.processEvents();time.sleep(.02)
    if js("document.getElementById('nbins').value==='400' && H[0].length===400"):break
else:raise AssertionError('Qt controls did not reach existing rebinFull')
assert js("document.getElementById('noise2d').value==='3'")
assert js("document.getElementById('msmode').value==='stick'")
assert js("document.getElementById('cmap').selectedOptions[0].textContent==='Plasma'")
width=js('P._fullLayout.width');assert 300<=width<=550
geometry=json.loads(js("JSON.stringify({toolbar:getComputedStyle(document.getElementById('toolbar')).display,exports:getComputedStyle(document.getElementById('exportpanel')).display,plot:P.getBoundingClientRect().height,layout:P._fullLayout.height,viewport:innerHeight,traces:P.data.length})"))
assert geometry['toolbar']=='none' and geometry['exports']=='none',geometry
assert geometry['traces']==3 and geometry['layout']==geometry['plot'],geometry
assert geometry['plot']<=geometry['viewport'] and geometry['plot']>=250,geometry
assert js("P._fullLayout.yaxis4._offset+P._fullLayout.yaxis4._length <= innerHeight")
w.web_view.load(QUrl('https://example.invalid/no-network'))
deadline=time.monotonic()+5
while not w.blocked_requests and time.monotonic()<deadline:
    app.processEvents();time.sleep(.02)
assert w.blocked_requests>=1
w.close();w.deleteLater();QCoreApplication.sendPostedEvents(None,QEvent.DeferredDelete)
app.processEvents();print('offline-controls-and-clean-close-ok')
'''
    env = os.environ.copy()
    env.pop("DECONVOVO_DISABLE_WEBENGINE", None)
    env.update(QT_QPA_PLATFORM="offscreen", DECONVOVO_ENABLE_HEADLESS_WEBENGINE="1",
        QTWEBENGINE_DISABLE_SANDBOX="1", QTWEBENGINE_CHROMIUM_FLAGS="--disable-gpu")
    result = subprocess.run([sys.executable, "-c", script, str(html), str(config_path)],
        cwd=Path(__file__).parents[1], env=env, capture_output=True, text=True, timeout=40)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "offline-controls-and-clean-close-ok" in result.stdout
    assert "Release of profile requested but WebEnginePage still not deleted" not in result.stderr
