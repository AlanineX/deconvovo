"""Native scientific arrays checked against independent vectors and HTML JS."""
from __future__ import annotations

from copy import deepcopy
import json
import hashlib
import os
from pathlib import Path
import re

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtQml import QJSEngine
from PySide6.QtWidgets import QApplication

from deconvovo.imms_view import load_view_data, _smooth_axis
from deconvovo.imms_html import plot_im_data


def config(**defaults):
    values = {"mz_bins": 8, "drift_bins": "native", "smooth_2d": {"method": "raw"},
        "smooth_drift": {"method": "raw"}, "noise_2d": 0, "noise_drift": 0,
        "noise_mz": 0, "denoise": "off", "scale": "1", "z_scale": "linear", "ms_mode": "profile"}
    values.update(defaults)
    cfg = json.loads((Path(__file__).parents[1] / "config/imms_plot_config.json").read_text(encoding="utf-8"))
    cfg["defaults"].update(values)
    cfg["initial_view"] = {"mz_range": [1000, 1003], "drift_range": [0, 0.9]}
    return cfg


@pytest.fixture(scope="module")
def application():
    app = QApplication.instance() or QApplication([])
    yield app
    app.processEvents()


def generated_html(directory, im, cfg, pusher=100):
    directory.mkdir(exist_ok=True)
    path = directory / "config.json"
    path.write_text(json.dumps(cfg), encoding="utf-8")
    plot_im_data(im, None, "sample", directory, pusher_us=pusher, config_path=path)
    return (directory / "sample_2d_imms.html").read_text(encoding="utf-8")


def js_engine(html, cfg, viewport=None):
    """Execute the existing generated numerical functions with a small Plotly stub."""
    bootstrap = "var H=" + re.search(r"\nvar H=(.*?)\nvar figData=", html, re.S).group(1)
    smoothing = "function sgSmooth1D" + html.split("function sgSmooth1D", 1)[1].split("// === Get visible zoom range ===", 1)[0]
    calculation = "function getVisibleRange()" + html.split("function getVisibleRange()", 1)[1].split("var _mzGuardState", 1)[0]
    dfl = cfg["defaults"]
    controls = {"nbins": str(dfl["mz_bins"]), "dbins": str(dfl.get("drift_bins", "native")),
        "sm2d": json.dumps(dfl["smooth_2d"]), "smdrift": json.dumps(dfl["smooth_drift"]),
        "mzscale": str(dfl["scale"]), "noise2d": str(dfl["noise_2d"]),
        "noisedrift": str(dfl["noise_drift"]), "noisemz": str(dfl["noise_mz"]),
        "denoise": str(dfl["denoise"]), "msmode": dfl["ms_mode"], "zscale": dfl["z_scale"],
        "render": "off", "cmap": json.dumps([[0, "#000000"], [1, "#ffffff"]])}
    view = viewport if viewport is not None else cfg.get("initial_view", {})
    setup = """
    var window={},controls=CONTROLS;
    var document={getElementById:function(id){return {value:controls[id]}}};
    var P={layout:{xaxis2:{range:MZRANGE},yaxis2:{range:DRRANGE}},data:[{},{},{}]};
    var Plotly={restyle:function(plot,changes,indices){
      if(!Array.isArray(indices))indices=[indices];
      Object.keys(changes).forEach(function(key){
        if(key.indexOf('.')>=0)return;
        indices.forEach(function(index,j){plot.data[index][key]=changes[key][j]});
      });return Promise.resolve();}};
    function setMzPanelGuard(){return Promise.resolve();}
    function applyAll(){apply2D();updateMarginals();return Promise.resolve();}
    """.replace("CONTROLS", json.dumps(controls)).replace("MZRANGE", json.dumps(view.get("mz_range"))).replace("DRRANGE", json.dumps(view.get("drift_range")))
    engine = QJSEngine()
    evaluated = engine.evaluate(setup + bootstrap + smoothing + calculation)
    assert not evaluated.isError(), evaluated.toString()
    return engine


def js_snapshot(engine, rebin=False):
    evaluated = engine.evaluate("rebinFull()" if rebin else "applyAll()")
    assert not evaluated.isError(), evaluated.toString()
    result = engine.evaluate("JSON.stringify({raw:H,mz:MZ,drift:DR,z:P.data[2].z,dp:P.data[0].x,ms_mz:P.data[1].x,ms_i:P.data[1].y})")
    assert not result.isError(), result.toString()
    return json.loads(result.toString())


@pytest.fixture
def sample(tmp_path):
    im = tmp_path / "sample_im.txt"
    points = [[1000, 0, 1], [1003, 8, 2]]
    for drift in range(9):
        for i in range(61):
            mz = 1000.2 + i * 0.04
            intensity = 2 + (30 if i in {12, 13, 14} else 0) + drift * (i % 7)
            points.append([mz, drift, intensity])
    np.savetxt(im, points, fmt="%.4f")
    return im


@pytest.mark.parametrize("defaults", [
    {},
    {"smooth_2d": {"method": "gaussian", "sigma": 0.18, "sigma_dr": 0.15},
     "smooth_drift": {"method": "gaussian", "sigma": 2}, "noise_2d": 30, "denoise": 15, "z_scale": "sqrt"},
    {"smooth_2d": {"method": "sg", "w": 7, "k": 2},
     "smooth_drift": {"method": "sg", "w": 5, "k": 2}, "scale": "sqrt", "noise_drift": 5, "z_scale": "cbrt"},
    {"smooth_2d": {"method": "ma", "w": 5}, "smooth_drift": {"method": "ma", "w": 4},
     "scale": "full", "noise_mz": 15, "z_scale": "log2", "ms_mode": "stick"},
    {"smooth_2d": {"method": "gaussian", "sigma": 2, "sigma_dr": 2},
     "smooth_drift": {"method": "gaussian", "sigma": 5}, "z_scale": "p4"},
    {"mz_bins": 17, "drift_bins": 5, "z_scale": "p5", "noise_2d": 5},
    {"drift_bins": 15, "noise_2d": 20, "denoise": "off", "z_scale": "log10"},
    {"drift_bins": 5, "smooth_2d": {"method": "gaussian", "sigma": 0.2, "sigma_dr": 0.17},
     "scale": "0.5"},
])
def test_native_arrays_match_independent_existing_html_js(application, sample, tmp_path, defaults):
    cfg = config(**defaults)
    html = generated_html(tmp_path / "html", sample, cfg)
    data = load_view_data(sample, pusher_us=100, config=cfg)
    native = data.render(cfg)
    reference = js_snapshot(js_engine(html, cfg), rebin=cfg["defaults"].get("drift_bins") != "native")
    for actual, key in [(native.raw2d, "raw"), (native.mz_centres, "mz"),
                        (native.drift_centres, "drift"), (native.display2d, "z"),
                        (native.drift_intensity, "dp"), (native.ms_mz, "ms_mz"),
                        (native.ms_intensity, "ms_i")]:
        np.testing.assert_allclose(actual, np.asarray(reference[key], dtype=float), rtol=1e-10, atol=1e-10, equal_nan=True)


@pytest.mark.parametrize("viewport", [
    {"mz_range": [1000.7, 1000.8], "drift_range": [0.12, 0.18]},
    {"mz_range": [1000, 1003], "drift_range": [2, 3]},
    {"mz_range": [999, 999.5], "drift_range": [0, 0.9]},
])
def test_partial_and_empty_viewports_match_browser_semantics(application, sample, tmp_path, viewport):
    cfg = config(mz_bins=4)
    html = generated_html(tmp_path / "html", sample, cfg)
    native = load_view_data(sample, pusher_us=100, config=cfg).render(cfg, viewport)
    reference = js_snapshot(js_engine(html, cfg, viewport))
    np.testing.assert_allclose(native.drift_intensity, reference["dp"])
    np.testing.assert_allclose(native.ms_mz, reference["ms_mz"])
    np.testing.assert_allclose(native.ms_intensity, reference["ms_i"])
    assert native.mz_range == tuple(viewport["mz_range"])
    assert native.drift_range == tuple(viewport["drift_range"])


def test_half_bin_centres_final_edge_full_precision_and_immutable_arrays(tmp_path):
    im = tmp_path / "sample_im.txt"
    np.savetxt(im, [[1000, 0, 0.125], [1500, 1, 20.75], [2000, 2, 30.5]], fmt="%.9f")
    cfg = config(mz_bins=2)
    cfg["initial_view"] = {"mz_range": [1000, 2000]}
    data = load_view_data(im, pusher_us=100, config=cfg)
    snapshot = data.render(cfg)
    np.testing.assert_allclose(snapshot.drift_centres, [0.05, 0.15, 0.25])
    np.testing.assert_allclose(snapshot.raw2d, [[0.125, 0], [0, 20.75], [0, 30.5]])
    assert snapshot.raw2d.sum() == 51.375
    assert snapshot.drift_unit == "ms"
    for array in (data.raw_mz, data.raw_drift, data.raw_intensity, snapshot.raw2d, snapshot.display2d):
        assert not array.flags.writeable
    bins = load_view_data(im, config=cfg).render(cfg)
    assert bins.drift_unit == "bins"
    np.testing.assert_array_equal(bins.drift_centres, [0.5, 1.5, 2.5])


def test_numeric_broadening_and_reset_rebin_original_txt_without_history(sample):
    cfg = config(mz_bins=8)
    cfg["initial_view"]["mz_range"] = [1000, 1001]
    data = load_view_data(sample, pusher_us=100, config=cfg)
    narrow = data.render(cfg)
    wide_cfg = deepcopy(cfg)
    wide_cfg["initial_view"]["mz_range"] = [1000, 1003]
    wide = data.render(wide_cfg)
    assert wide.raw2d.sum() > narrow.raw2d.sum()
    wide_cfg["defaults"]["mz_bins"] = 25
    data.render(wide_cfg)
    repeated = data.render(cfg)
    np.testing.assert_array_equal(repeated.raw2d, narrow.raw2d)
    np.testing.assert_array_equal(repeated.display2d, narrow.display2d)
    reset = deepcopy(cfg)
    reset["initial_view"]["mz_range"] = None
    full = data.render(reset)
    assert full.baked_mz_range == data.full_mz_range
    assert full.raw2d.sum() == pytest.approx(data.raw_intensity.sum())


def test_sg_constant_gain_and_known_impulse_are_corrected_in_both_runtimes(application, sample, tmp_path):
    cfg = config()
    engine = js_engine(generated_html(tmp_path / "html", sample, cfg), cfg)
    constant = np.full(31, 7.0)
    impulse = np.zeros(31)
    impulse[15] = 1
    expected = np.zeros(31)
    expected[12:19] = np.maximum(0, np.array([-2, 3, 6, 7, 6, 3, -2]) / 21)
    for signal, truth in [(constant, constant), (impulse, expected)]:
        np.testing.assert_allclose(_smooth_axis(signal, {"method": "sg", "w": 7, "k": 2}), truth, atol=1e-12)
        result = engine.evaluate("JSON.stringify(sgSmooth1D(" + json.dumps(signal.tolist()) + ",7,2))")
        assert not result.isError(), result.toString()
        np.testing.assert_allclose(json.loads(result.toString()), truth, atol=1e-12)
    # A constant-one signal directly checks that the implemented coefficients
    # sum to one, including replicated boundary samples.
    result = engine.evaluate("JSON.stringify(sgSmooth1D(" + json.dumps(np.ones(31).tolist()) + ",7,2))")
    np.testing.assert_allclose(json.loads(result.toString()), 1, atol=1e-12)


def test_physical_drift_gaussian_width_stays_stable_when_rebinned(application, tmp_path):
    im = tmp_path / "impulse_im.txt"
    np.savetxt(im, [[1000, 0, 1], [1500, 400, 1000], [2000, 799, 1]])
    cfg = config(mz_bins=8, smooth_2d={"method": "gaussian", "sigma": 0, "sigma_dr": 1.2})
    cfg["initial_view"] = {"mz_range": [1000, 2000]}
    data = load_view_data(im, pusher_us=100, config=cfg)
    html = generated_html(tmp_path / "html", im, cfg)
    widths = []
    for bins in ["native", 400, 200]:
        cfg["defaults"]["drift_bins"] = bins
        snapshot = data.render(cfg)
        profile = snapshot.display2d[:, 4]
        reference = js_snapshot(js_engine(html, cfg), rebin=bins != "native")
        np.testing.assert_allclose(snapshot.display2d, reference["z"], atol=1e-9)
        mean = np.sum(snapshot.drift_centres * profile) / profile.sum()
        width = np.sqrt(np.sum((snapshot.drift_centres - mean) ** 2 * profile) / profile.sum())
        widths.append(width)
    assert widths == pytest.approx([1.2] * 3, rel=0.002)
    assert max(widths) - min(widths) < 0.002


def test_large_heatmap_is_rejected_before_allocation(sample, monkeypatch):
    cfg = config(mz_bins=1_000_000, drift_bins=100)
    data = load_view_data(sample, config=cfg)
    monkeypatch.setattr(np, "bincount", lambda *a, **k: pytest.fail("Oversize heatmap was allocated"))
    with pytest.raises(ValueError, match="Choose fewer"):
        data.render(cfg)


def test_pathological_gaussian_padding_is_rejected_without_allocation():
    with pytest.raises(ValueError, match="window is too large"):
        _smooth_axis(np.ones((3, 3000)), {"method": "gaussian", "sigma": 1000}, 0)


def test_real_betaine_native_precision_and_html_parity(application, tmp_path):
    root = Path(__file__).parents[1]
    raw_folder = root / "portable_waters_convert/examples/production/small_molecule/raw/20260320_BETAINE_V2.raw"
    original_hashes = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in raw_folder.iterdir() if p.is_file()}
    im = root / "build/waters_console_validation/editor_readability/betaine_detail/_converted/20260320_BETAINE_V2_im.txt"
    if not im.is_file():
        from deconvovo.waters_workflow import run

        converted = tmp_path / "real_betaine"
        run(raw_folder, converted, html=False)
        im = converted / "_converted/20260320_BETAINE_V2_im.txt"
        assert not list(converted.glob("*.html"))
    cfg = config(mz_bins=200, smooth_2d={"method": "gaussian", "sigma": 0.1, "sigma_dr": 0.07},
        noise_2d=5, denoise=1, z_scale="cbrt")
    cfg["initial_view"] = {"mz_range": [1000, 2000], "drift_range": [4, 12]}
    data = load_view_data(im, pusher_us=110, config=cfg)
    native = data.render(cfg)
    reference = js_snapshot(js_engine(generated_html(tmp_path / "real_html", im, cfg, 110), cfg))
    np.testing.assert_allclose(native.raw2d, reference["raw"], atol=1e-10)
    np.testing.assert_allclose(native.display2d, reference["z"], atol=1e-10)
    np.testing.assert_allclose(native.drift_intensity, reference["dp"], atol=1e-10)
    np.testing.assert_allclose(native.ms_mz, reference["ms_mz"], atol=1e-10)
    assert native.ms_intensity.sum() == pytest.approx(sum(reference["ms_i"]), rel=0.001)
    difference = np.abs(native.ms_intensity - reference["ms_i"]).sum() / max(1, native.ms_intensity.sum())
    assert difference < 0.005  # HTML's 4-dp serialization can move a boundary sample
    raw = np.loadtxt(im)
    np.testing.assert_array_equal(data.raw_mz, raw[raw[:, 2] > 0, 0])
    np.testing.assert_array_equal(data.raw_intensity, raw[raw[:, 2] > 0, 2])
    assert original_hashes == {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in raw_folder.iterdir() if p.is_file()}
