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

from gui.widgets.waters_viewer import WatersViewerWidget, _config_script


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
def test_numeric_range_requires_finite_ordered_bounds(bounds):
    with pytest.raises(ValueError, match="finite minimum"):
        _config_script({"initial_view": {"mz_range": bounds}})


def test_current_export_requires_a_live_viewer_and_known_format(widget):
    with pytest.raises(ValueError, match="Unknown"):
        widget.export_current("raw_csv")
    with pytest.raises(ValueError, match="Process & view"):
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


@pytest.mark.parametrize("module_path", ["deconvovo/imms_html.py", "portable_waters_convert/deconvovo/imms_html.py"])
def test_native_current_exports_capture_settled_controls_and_numeric_ranges(tmp_path, module_path):
    """Exercise real Plotly, Qt downloads, CSV bytes and rendered PNGs together."""
    import importlib.util

    spec = importlib.util.spec_from_file_location("_waters_html_export_test", Path(__file__).parents[1] / module_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    config = json.loads((Path(__file__).parents[1] / "config/imms_plot_config.json").read_text(encoding="utf-8"))
    config["defaults"].update({"mz_bins": 4, "drift_bins": "native", "smooth_2d": {"method": "raw"},
        "smooth_drift": {"method": "raw"}, "noise_2d": 0, "denoise": "off", "noise_mz": 0,
        "noise_drift": 0, "ms_mode": "profile", "z_scale": "linear"})
    config["initial_view"] = {"mz_range": [1000, 2000], "drift_range": [0, 0.4]}
    cfg = tmp_path / "config.json"
    cfg.write_text(json.dumps(config), encoding="utf-8")
    im, ms = tmp_path / "sample_im.txt", tmp_path / "sample_ms.txt"
    points = [[950, 0, 50], [1000, 0, 4], [1250, 0, 9], [1500, 1, 25],
        [1750, 2, 100], [2000, 3, 16], [2100, 3, 50]]
    np.savetxt(im, points)
    np.savetxt(ms, [[p[0], p[2]] for p in points])
    original = (im.read_bytes(), ms.read_bytes())
    module.plot_im_data(im, ms, "sample", tmp_path, config_path=cfg, pusher_us=100)
    html = tmp_path / "sample_2d_imms.html"
    script = r'''
import csv,json,sys,time
from pathlib import Path
from PIL import Image
from PySide6.QtCore import QCoreApplication,QEvent
from PySide6.QtWidgets import QApplication,QFileDialog
from PySide6.QtWebEngineCore import QWebEnginePage,QWebEngineDownloadRequest
import gui.widgets.waters_viewer as viewer
errors=[]
class CheckedPage(QWebEnginePage):
    def javaScriptConsoleMessage(self,level,message,line,source):
        if level==QWebEnginePage.JavaScriptConsoleMessageLevel.ErrorMessageLevel:
            errors.append({'message':message,'line':line,'source':source})
viewer.QWebEnginePage=CheckedPage
WatersViewerWidget=viewer.WatersViewerWidget;WEBENGINE_AVAILABLE=viewer.WEBENGINE_AVAILABLE
folder=Path(sys.argv[1]);out=folder/'captures';out.mkdir()
app=QApplication([]);assert WEBENGINE_AVAILABLE
w=WatersViewerWidget();w.resize(900,650);w.show()
completed={}
def track(download):
    def changed(state):
        if state==QWebEngineDownloadRequest.DownloadState.DownloadCompleted:
            name=download.downloadFileName();completed[name]=completed.get(name,0)+1
    download.stateChanged.connect(changed)
w._profile.downloadRequested.connect(track)
cfg=json.loads((folder/'config.json').read_text(encoding='utf-8'))
w.set_run({'run_name':'sample','viewer':folder/'sample_2d_imms.html',
    'ms':folder/'sample_ms.txt','im':folder/'sample_im.txt'},cfg)
QFileDialog.getSaveFileName=staticmethod(lambda parent,title,suggested,filter:(str(out/Path(suggested).name),filter))
def js(code):
    result=[];w.web_view.page().runJavaScript(code,result.append)
    deadline=time.monotonic()+5
    while not result and time.monotonic()<deadline:app.processEvents();time.sleep(.01)
    assert result,code
    return result[0]
def wait(predicate,seconds=15):
    deadline=time.monotonic()+seconds
    while time.monotonic()<deadline:
        app.processEvents();time.sleep(.02)
        if predicate():return
    raise AssertionError('Timed out waiting for native current-view export: '+json.dumps(errors))
def captured(name):return (out/name).exists() and completed.get(name,0)>0
wait(lambda:w._html_ready and js("typeof P!=='undefined' && !!P._fullLayout"),25)
wait(lambda:js("P.data[2].z.length===4"))
js("window.pngAxes=[];var originalToImage=Plotly.toImage;Plotly.toImage=function(gd,opts){if(gd!==P)pngAxes.push({xTicks:gd.querySelectorAll('.xtick text').length,yTicks:gd.querySelectorAll('.ytick text').length});return originalToImage.apply(this,arguments)};true")
w.export_current('full_png');full=out/'sample_current_full.png';wait(lambda:captured(full.name))
first=full.read_bytes()
cfg['defaults'].update({'mz_bins':8,'drift_bins':2,'z_scale':'sqrt','noise_2d':30,'denoise':'off',
    'noise_drift':20,'noise_mz':10})
cfg['initial_view']={'mz_range':[1200,1800],'drift_range':[0.1,0.3]}
w.apply_config(cfg)
# Capture immediately, before the debounce interval: export must flush and wait.
w.export_current('heatmap_csv');heatmap=out/'sample_current_heatmap.csv';wait(lambda:captured(heatmap.name))
state=json.loads(js("JSON.stringify({hm:P.data[2],drift:P.data[0],mz:P.data[1],vr:getVisibleRange(),layout:{xr:P._fullLayout.xaxis2.range,xr2:P._fullLayout.xaxis4.range,yr:P._fullLayout.yaxis.range,yr2:P._fullLayout.yaxis2.range}})"))
assert state['layout']['xr']==state['layout']['xr2']==[1200,1800],state['layout']
assert state['layout']['yr']==state['layout']['yr2']==[0.1,0.3],state['layout']
assert len(state['hm']['x'])==8 and len(state['hm']['y'])==2
rows=list(csv.DictReader(heatmap.open(encoding='utf-8')))
assert set(rows[0])=={'mz','drift_time_ms','display_value','intensity_transform'}
vr=state['vr'];expected=[]
for d in range(vr['dBl'],vr['dBh']+1):
    for m in range(vr['mBl'],vr['mBh']+1):
        expected.append((state['hm']['x'][m],state['hm']['y'][d],state['hm']['z'][d][m]))
actual=[tuple(float(row[key]) for key in ('mz','drift_time_ms','display_value')) for row in rows]
assert actual==expected,(actual,expected)
assert all(row['intensity_transform']=='sqrt' for row in rows)
assert any(row[2]==0 for row in actual) and any(row[2]>0 for row in actual)
for kind in ['mz_csv','drift_csv','heatmap_png','mz_png','drift_png']:
    w.export_current(kind);wait(lambda k=kind:captured('sample_current_'+k.replace('_','.')))
mzrows=list(csv.DictReader((out/'sample_current_mz.csv').open(encoding='utf-8')))
assert [(float(r['mz']),float(r['display_value'])) for r in mzrows]==[
    (x,y) for x,y in zip(state['mz']['x'],state['mz']['y']) if x is not None and 1200<=x<=1800]
driftrows=list(csv.DictReader((out/'sample_current_drift.csv').open(encoding='utf-8')))
assert [(float(r['drift_time_ms']),float(r['display_value'])) for r in driftrows]==[
    (y,x) for x,y in zip(state['drift']['x'],state['drift']['y']) if y is not None and .1<=y<=.3]
w.export_current('full_png');wait(lambda:completed.get(full.name,0)==2)
assert full.read_bytes()!=first
for name in ['full','heatmap','mz','drift']:
    image=Image.open(out/('sample_current_'+name+'.png'));image.verify()
    assert image.width>=1000 and image.height>=700,(name,image.size)
pngAxes=json.loads(js('JSON.stringify(pngAxes)'))
assert len(pngAxes)==3 and all(axes['xTicks']>0 and axes['yTicks']>0 for axes in pngAxes),pngAxes
assert w.blocked_requests==0
assert not errors,errors
# Changing a display control must preserve a zoom made with the mouse.
js("Plotly.relayout(P,{'xaxis2.range':[1300,1600],'xaxis4.range':[1300,1600]});true")
wait(lambda:js("P._fullLayout.xaxis2.range[0]===1300"))
cfg['defaults']['noise_2d']=50;w.apply_config(cfg)
wait(lambda:js("document.getElementById('noise2d').value==='50'"))
assert json.loads(js('JSON.stringify(P._fullLayout.xaxis2.range)'))==[1300,1600]
# Clearing the numeric range explicitly resets both linked dimensions.
cfg['initial_view']={'mz_range':None,'drift_range':None};w.apply_config(cfg)
wait(lambda:js("P._fullLayout.xaxis2.range[0]===1000 && P._fullLayout.yaxis2.range[0]===0"))
assert json.loads(js('JSON.stringify(P._fullLayout.xaxis4.range)'))==[1000,2000]
assert json.loads(js('JSON.stringify(P._fullLayout.yaxis.range)'))==[0,.4]
cfg['initial_view']['mz_range']=[900,2100];w.apply_config(cfg)
wait(lambda:w.status.isVisible() and 'regenerate' in w.status.text())
assert json.loads(js('JSON.stringify(P._fullLayout.xaxis2.range)'))==[1000,2000]
cfg['initial_view']['mz_range']=[1200,1800];w.apply_config(cfg)
wait(lambda:js("P._fullLayout.xaxis2.range[0]===1200") and not w.status.isVisible())
assert not errors,errors
w.close();w.deleteLater();QCoreApplication.sendPostedEvents(None,QEvent.DeferredDelete)
app.processEvents();print('settled-current-exports-and-numeric-ranges-ok')
'''
    env = os.environ.copy()
    env.pop("DECONVOVO_DISABLE_WEBENGINE", None)
    env.update(QT_QPA_PLATFORM="offscreen", DECONVOVO_ENABLE_HEADLESS_WEBENGINE="1",
        QTWEBENGINE_DISABLE_SANDBOX="1", QTWEBENGINE_CHROMIUM_FLAGS="--disable-gpu")
    result = subprocess.run([sys.executable, "-c", script, str(tmp_path)],
        cwd=Path(__file__).parents[1], env=env, capture_output=True, text=True, timeout=120)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "settled-current-exports-and-numeric-ranges-ok" in result.stdout
    assert (im.read_bytes(), ms.read_bytes()) == original
    assert "Release of profile requested but WebEnginePage still not deleted" not in result.stderr
