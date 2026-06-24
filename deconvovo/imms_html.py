"""Interactive 2D IM-MS HTML viewer with linked marginal panels.

Generates a single self-contained HTML file with:
- 2D heatmap (center), drift profile (left), m/z spectrum (bottom)
- Client-side smoothing, noise filtering, colormap selection
- Zoom-linked marginals with rAF debouncing
- Export buttons for PNG + CSV (raw + processed)
"""
from __future__ import annotations

import json as _json
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.colors as _pc
from plotly.subplots import make_subplots
from scipy.signal import savgol_filter
from scipy.ndimage import gaussian_filter

def _load_plot_config(config_path: Path | str | None = None) -> dict:
    """Load imms_plot_config.json. Falls back to defaults if missing.

    If `config_path` is given (absolute or relative), that file is used.
    Otherwise looks for `config/imms_plot_config.json` in the repo / frozen bundle.
    """
    import sys
    if config_path is not None:
        cfg_path = Path(config_path).expanduser().resolve()
    else:
        if getattr(sys, 'frozen', False):
            base = Path(sys._MEIPASS)
        else:
            base = Path(__file__).parent.parent
        cfg_path = base / "config" / "imms_plot_config.json"
    default = {
        "defaults": {
            "smooth_2d": {"method": "raw"}, "scale": "1",
            "smooth_drift": {"method": "raw"}, "noise_pct": 0,
            "mz_bins": 800, "colormap": "Viridis",
        },
        "presets": {
            "smooth_2d": [{"method": "raw", "label": "Raw"}],
            "scale": ["1", "2", "sqrt", "full"],
            "smooth_drift": [{"method": "raw", "label": "Raw"}],
            "noise_pct": [0, 0.5, 1, 3, 5],
            "mz_bins": [200, 400, 800, 1600, 3200],
            "colormap": ["Viridis", "Plasma", "Inferno", "Cividis",
                         "Hot", "Turbo", "Electric", "Bluered"],
        },
        "figure": {"font_size": 20, "width": 1190, "height": 680},
    }
    if cfg_path.exists():
        import json
        try:
            with open(cfg_path) as f:
                return json.load(f)
        except Exception:
            pass
    return default


def plot_im_data(im_file: Path, ms_file: Path | None, run_name: str, out_dir: Path,
                  n_mz_bins: int | None = None, pusher_us: float | None = None,
                  config_path: Path | str | None = None):
    """Interactive 2D IM-MS viewer with linked marginals.

    Layout:
        [drift profile (left)]  [2D heatmap (center)]
                                [m/z spectrum (bottom)]

    Key design:
      - m/z panel uses RAW _ms.txt data (full isotopic resolution), NOT the binned heatmap
      - Drift and 2D smoothing are INDEPENDENT dropdowns
      - Savitzky-Golay smoothing for 1D traces (preserves peak shape)
      - Gaussian smoothing for 2D heatmap (fills sparse gaps)
      - Colormap chooser
      - All axes linked for zoom/pan
    """
    from scipy.signal import savgol_filter
    from plotly.subplots import make_subplots
    import json as _json

    # Load config (custom path overrides the repo default)
    cfg = _load_plot_config(config_path)
    dfl = cfg["defaults"]
    prs = cfg["presets"]

    def _as_float(value, fallback=0.0):
        try:
            return float(value)
        except (TypeError, ValueError):
            return fallback

    if _as_float(dfl.get("noise_2d")) <= 0 < _as_float(dfl.get("denoise")):
        first_positive_noise = next(
            (_as_float(v) for v in prs.get("noise_2d", []) if _as_float(v) > 0),
            0.5,
        )
        dfl["noise_2d"] = first_positive_noise

    fig_cfg = cfg.get("figure", {})
    FS = fig_cfg.get("font_size", 20)
    FIG_W = fig_cfg.get("width", 1190)
    FIG_H = fig_cfg.get("height", 680)
    if n_mz_bins is None:
        n_mz_bins = dfl.get("mz_bins", 800)
    # "native" = use as many m/z bins as the data actually supports — counted
    # below from the unique m/z values in the IM data, capped to keep heatmap
    # serialization tractable.

    # --- Load IM data ---
    data = np.loadtxt(str(im_file))
    if data.ndim != 2 or data.shape[1] < 3:
        return

    mz_raw, drift_raw, int_raw = data[:, 0], data[:, 1].astype(int), data[:, 2]
    nz = int_raw > 0
    mz_nz, drift_nz, int_nz = mz_raw[nz], drift_raw[nz], int_raw[nz]
    n_points = len(mz_nz)
    if n_points == 0:
        return

    # Compute drift time if pusher period is known
    if pusher_us:
        drift_time_nz = (drift_nz.astype(float) + 0.5) * pusher_us / 1000.0
        pd.DataFrame({"mz": mz_nz, "drift_bin": drift_nz,
                       "drift_time_ms": drift_time_nz, "intensity": int_nz}).to_csv(
            out_dir / f"{run_name}_2d_imms.csv", index=False, float_format="%.4f")
    else:
        pd.DataFrame({"mz": mz_nz, "drift_bin": drift_nz, "intensity": int_nz}).to_csv(
            out_dir / f"{run_name}_2d_imms.csv", index=False, float_format="%.4f")

    # --- Build heatmap grid ---
    # The heatmap m/z range is CLIPPED to initial_view.mz_range when set, so
    # the same N bins serve only the visible window — 5x sharper than spreading
    # them over the full data range. Out-of-window samples are dropped from the
    # heatmap but stay in `rawIM` for client-side re-binning if the user zooms
    # to a wider window.
    mz_lo, mz_hi = float(mz_nz.min()), float(mz_nz.max())
    _iv_pre = cfg.get("initial_view") or {}
    if _iv_pre.get("mz_range"):
        _hm_lo, _hm_hi = _iv_pre["mz_range"]
        hm_mz_range = (max(mz_lo, float(_hm_lo)), min(mz_hi, float(_hm_hi)))
        _hm_mask = (mz_nz >= hm_mz_range[0]) & (mz_nz <= hm_mz_range[1])
    else:
        hm_mz_range = (mz_lo, mz_hi)
        _hm_mask = slice(None)
    n_drift = int(drift_nz.max()) + 1

    # Native bin count = unique m/z values inside the heatmap m/z range, capped
    # at 25 600 so the embedded heatmap stays under ~5 M cells × 8 bytes ≈ 40 MB
    # JSON. Clipping first means "Native" is native at the visible scale, not
    # diluted by samples outside the view.
    n_mz_native = int(min(25600, max(200, len(np.unique(mz_nz[_hm_mask])))))
    if isinstance(n_mz_bins, str) and n_mz_bins.lower() == "native":
        n_mz_bins = n_mz_native
    mz_edges = np.linspace(hm_mz_range[0], hm_mz_range[1], n_mz_bins + 1)
    mz_centers = (mz_edges[:-1] + mz_edges[1:]) / 2
    heatmap = np.zeros((n_drift, n_mz_bins), dtype=np.float64)
    bin_idx = np.searchsorted(mz_edges, mz_nz) - 1
    valid = (bin_idx >= 0) & (bin_idx < n_mz_bins)
    for d, bi, intensity in zip(drift_nz[valid], bin_idx[valid], int_nz[valid]):
        heatmap[int(d), int(bi)] += intensity

    # No auto-trim — show full data range, user can zoom interactively
    sub_raw = heatmap
    sub_mz = mz_centers
    bw = (sub_mz[-1] - sub_mz[0]) / len(sub_mz) if len(sub_mz) > 1 else 0
    drift_bins_int = np.arange(n_drift)
    if pusher_us:
        drift_bins = drift_bins_int * pusher_us / 1000  # convert to ms
        drift_label = "Drift Time (ms)"
    else:
        drift_bins = drift_bins_int.astype(float)
        drift_label = "Drift Time (bins)"

    # --- Load raw MS spectrum for m/z panel (full isotopic resolution) ---
    # This is the CDCReader _ms.txt with ~0.01 Da spacing — NOT the binned heatmap.
    # Profile MS data naturally has peak widths from the TOF instrument resolution.
    # Zoom in to ~1-2 Da range to see individual isotope peaks.
    if ms_file and Path(ms_file).exists():
        ms_data = np.loadtxt(str(ms_file))
        ms_mz, ms_int = ms_data[:, 0], ms_data[:, 1]
    else:
        ms_mz, ms_int = sub_mz, sub_raw.sum(axis=0)

    # --- Raw IM data grouped by drift bin (for independent drift profile) ---
    # Analogous to _ms.txt for the m/z panel: raw data, not rebinned heatmap.
    # rawIM[d] = sorted array of [mz, intensity] pairs at drift bin d.
    from collections import defaultdict
    _drift_groups = defaultdict(list)
    for _mz, _dr, _it in zip(mz_nz, drift_nz, int_nz):
        _drift_groups[int(_dr)].append([round(float(_mz), 4), round(float(_it))])
    _raw_im = []
    for _d in range(n_drift):
        pts = _drift_groups.get(_d, [])
        pts.sort(key=lambda p: p[0])  # sort by m/z for fast range filtering
        _raw_im.append(pts)

    # Native m/z sampling interval — the floor for the m/z panel's profile grid.
    # Profile mode must NOT bin finer than this or the connecting line dips to
    # zero between native samples and renders as a comb of spikes (looks like
    # sticks). Use the 95th percentile of the within-peak sample spacing
    # (diffs < 0.1 Da exclude the large inter-peak gaps); ×1.3 safety for the
    # m/z-dependent growth of TOF sampling. Falls back to 0.01 Da.
    _umz = np.unique(mz_nz)
    _dmz = np.diff(_umz)
    _dmz = _dmz[(_dmz > 0) & (_dmz < 0.1)]
    native_mz_step = float(np.percentile(_dmz, 95) * 1.3) if len(_dmz) else 0.01

    # --- Raw drift profile ---
    drift_tic_raw = sub_raw.sum(axis=1)

    # --- Build figure with MINIMAL traces ---
    # All smoothing/noise is done CLIENT-SIDE in JS from the embedded raw data.
    # This eliminates the smoothing×noise cross-product problem and enables
    # truly custom parameter values.
    #
    # Traces:  [0] drift profile  [1] m/z spectrum  [2] 2D heatmap
    # All three get updated by JS when controls change.

    fig = make_subplots(
        rows=2, cols=2,
        row_heights=[0.76, 0.24],
        column_widths=[0.14, 0.86],
        shared_xaxes="columns",
        shared_yaxes="rows",
        horizontal_spacing=0.012,
        vertical_spacing=0.012,
    )

    DRIFT_IDX = 0
    MZ_IDX = 1
    HM_IDX = 2

    # Drift profile (left)
    fig.add_trace(go.Scatter(
        x=drift_tic_raw, y=drift_bins, mode="lines",
        line=dict(color="black", width=0.7), showlegend=False,
        hovertemplate="Drift: %{y:.3f}<br>TIC: %{x:.0f}<extra></extra>",
    ), row=1, col=1)

    # m/z spectrum (bottom) — initial: profile trace (drift-summed, native resolution)
    # so first paint shows data even before any user interaction. updateMarginals()
    # replaces this whenever the user pans/zooms or toggles MS Mode.
    #
    # CLIP to the initial m/z view if set. Otherwise Plotly's y-axis autoscale
    # uses the FULL trace's max (potentially a peak way outside the visible
    # window — Myoglobin native at m/z 1953 falsely collapses visible 350–1250
    # peaks to ~20% of the axis, looking like nothing is shown).
    _iv = cfg.get("initial_view") or {}
    _ms_mz_init = ms_mz
    _ms_int_init = ms_int
    if _iv.get("mz_range"):
        _lo, _hi = _iv["mz_range"]
        _mask = (ms_mz >= _lo) & (ms_mz <= _hi)
        if _mask.any():
            _ms_mz_init = ms_mz[_mask]
            _ms_int_init = ms_int[_mask]
    fig.add_trace(go.Scatter(
        x=_ms_mz_init.tolist(), y=_ms_int_init.tolist(), mode="lines",
        line=dict(color="#333", width=0.6), showlegend=False,
        hovertemplate="m/z: %{x:.4f}<br>I: %{y:.0f}<extra></extra>",
    ), row=2, col=2)

    # 2D heatmap (center) — initial: RAW (no smoothing)
    # zsmooth controls Plotly's inter-cell interpolation:
    #   "off"  -> False, each cell is a discrete pixel (crisp; best at high mz_bins)
    #   "fast" -> bilinear interpolation
    #   "best" -> bicubic interpolation (smooth but blurry for high-res data)
    _render = dfl.get("render", "off")
    _zsmooth = False if _render == "off" else _render
    _zscale = dfl.get("z_scale", "log10")
    _zpos = np.maximum(0, sub_raw)
    if _zscale == "linear":
        _z_init = _zpos; _z_title = "I"; _hover_z = "I: %{z:.0f}"
    elif _zscale == "sqrt":
        _z_init = np.sqrt(_zpos); _z_title = "√I"; _hover_z = "√I: %{z:.2f}"
    elif _zscale == "cbrt":
        _z_init = np.cbrt(_zpos); _z_title = "∛I"; _hover_z = "∛I: %{z:.2f}"
    elif _zscale == "p4":
        _z_init = np.power(_zpos, 0.25); _z_title = "I^¼"; _hover_z = "I^¼: %{z:.2f}"
    elif _zscale == "p5":
        _z_init = np.power(_zpos, 0.2); _z_title = "I^⅕"; _hover_z = "I^⅕: %{z:.2f}"
    elif _zscale == "log2":
        _z_init = np.log2(_zpos + 1); _z_title = "log₂(I+1)"; _hover_z = "log₂(I+1): %{z:.2f}"
    else:  # log10
        _z_init = np.log10(_zpos + 1); _z_title = "log₁₀(I+1)"; _hover_z = "log(I+1): %{z:.2f}"
    # Resolve the default colormap to a Plotly scale. It may be a CUSTOM name
    # (e.g. "W-Rainbow") defined in presets["colormap"] as a {name, scale} dict —
    # Plotly wouldn't recognise the bare name. init applyAll() re-applies the
    # dropdown's resolved scale on first paint regardless; this keeps the
    # pre-JS frame correct too.
    _default_cmap = dfl.get("colormap", "Viridis")
    _default_scale = _default_cmap
    for _cm in prs.get("colormap", []):
        if isinstance(_cm, dict) and _cm.get("name") == _default_cmap:
            _default_scale = _cm["scale"]
            break
    fig.add_trace(go.Heatmap(
        z=_z_init, x=sub_mz, y=drift_bins,
        colorscale=_default_scale,
        colorbar=dict(title=dict(text=_z_title, side="top"),
                      len=0.76, y=0.62, thickness=22),
        zsmooth=_zsmooth,
        hovertemplate=f"m/z: %{{x:.2f}}<br>Drift: %{{y:.3f}}<br>{_hover_z}<extra></extra>",
    ), row=1, col=2)

    # No Plotly dropdowns — all controls are HTML elements wired to JS
    cmaps = prs.get("colormap", ["Viridis", "Plasma", "Inferno", "Cividis",
                                    "Hot", "Turbo", "Electric", "Bluered"])

    def _cmap_label(cm):
        return str(cm.get("name", "")) if isinstance(cm, dict) else str(cm)

    def _natural_sort_key(text: str):
        import re
        return [int(p) if p.isdigit() else p.lower()
                for p in re.split(r"(\d+)", text)]

    cmaps = sorted(
        cmaps,
        key=lambda cm: (
            0 if _cmap_label(cm).startswith("W-") else 1,
            _natural_sort_key(_cmap_label(cm)),
        ),
    )

    # =================================================================
    # Publication-quality styling
    # =================================================================
    FF = "Liberation Sans, Arial, Helvetica, sans-serif"
    # FS, FIG_W, FIG_H loaded from config above
    LC = "#333"; AC = "#444"; GC = "#e0e0e0"

    fig.update_layout(
        font=dict(family=FF, size=FS, color=LC),
        template="plotly_white", width=FIG_W, height=FIG_H,
        margin=dict(l=80, r=100, t=30, b=70),
        plot_bgcolor="white", paper_bgcolor="white",
    )
    ax = dict(showline=True, linewidth=1, linecolor=AC, mirror=True)
    gm = dict(showgrid=True, gridwidth=0.5, gridcolor=GC, griddash="dot")
    go_ = dict(showgrid=False)
    tf = dict(family=FF, size=int(FS * 1.2))  # axis title font
    tkf = dict(size=FS)                        # axis tick font
    tkfs = dict(size=int(FS * 0.9))            # smaller tick font

    fig.update_xaxes(row=2, col=2, title_text="<i>m/z</i>",
                     title_font=tf, tickfont=tkf, **ax, **gm,
                     minor=dict(showgrid=True, gridwidth=0.3, gridcolor="#eee", griddash="dot"))
    fig.update_yaxes(row=2, col=2, rangemode="tozero",
                     matches=None,   # un-share from yaxis3 (the empty corner,
                                     # which would force this to [0, 1] and
                                     # hide the m/z intensity spectrum)
                     fixedrange=True,  # lock intensity to autoscale; users
                                       # zoom the m/z panel along x only
                                       # (matches the original behaviour)
                     title_text="Intensity", title_font=tf,
                     tickfont=tkfs, **ax, **go_)
    fig.update_yaxes(row=1, col=1, title_text=drift_label,
                     title_font=tf, tickfont=tkf, **ax, **gm,
                     minor=dict(showgrid=True, gridwidth=0.3, gridcolor="#eee", griddash="dot"))
    fig.update_xaxes(row=1, col=1, autorange="reversed",
                     autorangeoptions=dict(include=[0]), tickfont=tkfs, **ax, **go_)
    fig.update_xaxes(row=1, col=2, showticklabels=False, **ax)
    fig.update_yaxes(row=1, col=2, showticklabels=False, **ax)

    # Optional initial-view ranges from config (interactive zoom still works
    # — these just set the starting m/z and drift-time windows).
    iv = cfg.get("initial_view") or {}
    if iv.get("mz_range"):
        lo, hi = iv["mz_range"]
        fig.update_xaxes(row=1, col=2, range=[lo, hi])
        fig.update_xaxes(row=2, col=2, range=[lo, hi])
    if iv.get("drift_range"):
        lo, hi = iv["drift_range"]
        fig.update_yaxes(row=1, col=1, range=[lo, hi])
        fig.update_yaxes(row=1, col=2, range=[lo, hi])
    fig.update_xaxes(row=2, col=1, showticklabels=False, fixedrange=True, showline=False, showgrid=False)
    fig.update_yaxes(row=2, col=1, showticklabels=False, fixedrange=True, showline=False, showgrid=False)

    # =================================================================
    # Write custom HTML: Plotly figure + HTML toolbar + all-JS smoothing
    # =================================================================

    # --- Build dropdown HTML from config ---
    def _smooth_options(presets, default):
        """Generate <option> tags for smoothing presets."""
        opts = []
        for p in presets:
            d = {k: v for k, v in p.items() if k != "label"}
            val = _json.dumps(d)
            label = p.get("label", str(d))
            # Replace ALL non-ASCII chars with HTML entities
            label = "".join(f"&#{ord(c)};" if ord(c) > 127 else c for c in label)
            sel = " selected" if d == default else ""
            opts.append(f"<option value='{val}'{sel}>{label}</option>")
        return "\n      ".join(opts)

    def _simple_options(values, default, fmt=None):
        # When no custom fmt is given, fall back to the legacy m/z-scale label
        # map (used by the Scale dropdown). When a fmt is given, trust it
        # entirely — otherwise "sqrt" gets mis-labeled as "√ratio" in
        # unrelated dropdowns (e.g. Intensity).
        legacy = {"0": "Off", "sqrt": "&radic;ratio", "full": "ratio"}
        opts = []
        for v in values:
            sel = " selected" if str(v) == str(default) else ""
            if fmt is not None:
                label = fmt(v)
            else:
                label = legacy.get(str(v), str(v))
            opts.append(f'<option value="{v}"{sel}>{label}</option>')
        return "\n      ".join(opts)

    sm2d_opts = _smooth_options(prs["smooth_2d"], dfl["smooth_2d"])
    smdrift_opts = _smooth_options(prs["smooth_drift"], dfl["smooth_drift"])
    scale_opts = _simple_options(prs["scale"], dfl["scale"])
    nfmt = lambda v: "Off" if v == 0 else f"{v}%"
    noise2d_opts = _simple_options(prs.get("noise_2d", [0, 0.5, 1, 3, 5]),
                                    dfl.get("noise_2d", 0), nfmt)
    # Cutoff Fade = off/hard threshold, or ramp half-width (% of max).
    denoise_opts = _simple_options(prs.get("denoise", ["off", 1, 5, 10, 20, 50]),
                                    dfl.get("denoise", 5),
                                    lambda v: "Off" if str(v) == "off" else f"{v}%")
    noise_dr_opts = _simple_options(prs.get("noise_drift", [0, 0.5, 1, 3, 5]),
                                     dfl.get("noise_drift", 0), nfmt)
    noise_mz_opts = _simple_options(prs.get("noise_mz", [0, 0.5, 1, 3, 5]),
                                     dfl.get("noise_mz", 0), nfmt)
    bins_opts = _simple_options(
        prs["mz_bins"], dfl["mz_bins"],
        lambda v: f"Native ({n_mz_native:,})" if str(v) == "native" else str(v),
    )
    # Intensity-transform dropdown (linear / sqrt / log10) — addresses the
    # diagnosis's "log10 compresses strong peaks" point.
    _zscale_labels = {"linear": "Linear (I)", "sqrt": "Sqrt (√I)", "cbrt": "Cbrt (∛I)",
                      "p4": "4th root (I^¼)", "p5": "5th root (I^⅕)",
                      "log2": "log₂(I+1)", "log10": "log₁₀(I+1)"}
    zscale_opts = _simple_options(
        prs.get("z_scale", ["linear", "sqrt", "cbrt", "p4", "p5", "log2", "log10"]),
        dfl.get("z_scale", "sqrt"),
        lambda v: _zscale_labels.get(str(v), str(v)),
    )
    # Render-mode dropdown: maps to Plotly Heatmap.zsmooth
    _render_labels = {"off": "Crisp", "fast": "Smooth (fast)", "best": "Smooth (best)"}
    render_opts = _simple_options(
        prs.get("render", ["off", "fast", "best"]),
        dfl.get("render", "off"),
        lambda v: _render_labels.get(str(v), str(v)),
    )
    _ms_mode_labels = {"stick": "Sticks", "profile": "Profile"}
    ms_mode_opts = _simple_options(
        prs.get("ms_mode", ["stick", "profile"]),
        dfl.get("ms_mode", "profile"),
        lambda v: _ms_mode_labels.get(str(v), str(v)),
    )
    dbins_opts = _simple_options(
        prs.get("drift_bins", ["native", 50, 100, 200]),
        dfl.get("drift_bins", "native"),
        lambda v: f"{v}" if v != "native" else "Native"
    )
    # Colormap options: resolve ALL names to arrays so Plotly.js doesn't need name lookup
    import plotly.colors as _pc

    def _resolve_colorscale(name: str, n: int = 32):
        """Plotly built-in first, then matplotlib by name; fall back to the
        string so Plotly.js can try at runtime."""
        try:
            return _pc.get_colorscale(name)
        except Exception:
            pass
        try:
            try:
                import matplotlib as _mpl
                cmap = _mpl.colormaps[name]
            except (KeyError, AttributeError):
                import matplotlib.cm as _mcm
                cmap = _mcm.get_cmap(name)
            from matplotlib.colors import to_hex as _to_hex
            return [[i / (n - 1), _to_hex(cmap(i / (n - 1)))] for i in range(n)]
        except Exception:
            return name

    cmap_opts_list = []
    for cm in cmaps:
        if isinstance(cm, dict):
            scale = cm["scale"]
            label = cm["name"]
            sel = " selected" if cm["name"] == dfl["colormap"] else ""
        else:
            scale = _resolve_colorscale(cm)
            label = cm
            sel = " selected" if cm == dfl["colormap"] else ""
        val = _json.dumps(scale)
        cmap_opts_list.append(f"<option value='{val}'{sel}>{label}</option>")
    cmap_opts = "\n      ".join(cmap_opts_list)

    hm_json = _json.dumps(sub_raw.tolist())
    mz_json = _json.dumps(sub_mz.tolist())
    dr_json = _json.dumps(drift_bins.tolist())
    ms_mz_j = _json.dumps(ms_mz.tolist())
    ms_int_j = _json.dumps(ms_int.tolist())
    raw_im_j = _json.dumps(_raw_im)
    fig_json = fig.to_json()
    cmaps_json = _json.dumps(cmaps)

    html = f"""<!DOCTYPE html><html><head><meta charset="utf-8">
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<style>
*{{box-sizing:border-box;margin:0;padding:0}}
body{{font-family:Liberation Sans,Arial,sans-serif;background:#fff}}
#toolbar{{display:flex;align-items:center;gap:18px;padding:8px 16px;
  background:#f8f8f8;border-bottom:1px solid #ddd;flex-wrap:wrap}}
#toolbar .ctl{{display:flex;align-items:center;gap:5px}}
#toolbar label{{font-size:{FS}px;color:#444;font-weight:600;white-space:nowrap}}
#toolbar select,#toolbar input{{font-size:{FS}px;padding:4px 8px;border:1px solid #bbb;
  border-radius:3px;background:#fff}}
#toolbar input[type=number]{{width:80px}}
#toolbar button{{font-size:{FS}px;padding:6px 18px;border:1px solid #999;border-radius:3px;
  background:#eee;cursor:pointer}}
#toolbar button:hover{{background:#ddd}}
#title{{font-size:{int(FS*1.3)}px;font-weight:700;padding:6px 16px 2px;color:#333}}
#subtitle{{font-size:{FS}px;color:#888;padding:0 16px 4px}}
#main{{height:calc(100vh - 140px);min-height:700px}}
#plot{{width:100%;height:100%}}
#exportpanel button:hover{{background:#ddd}}
</style></head><body>
<div id="title">{run_name}</div>
<div id="subtitle">{n_points:,} IM points &middot; {len(ms_mz):,} MS points &middot; {f'{pusher_us:.2f} &mu;s/bin' if pusher_us else 'uncalibrated'} &middot;
  Double-click to reset zoom &middot; &#x1F3E0; resets all</div>
<div id="toolbar">
  <div class="ctl">
    <label>2D Smooth:</label>
    <select id="sm2d">
      {sm2d_opts}
    </select>
  </div>
  <div class="ctl">
    <label>Scale:</label>
    <select id="mzscale" title="Multiplier for m/z-axis 2D smoothing. Lower values make peaks thinner in m/z; higher values make them fatter in m/z.">
      {scale_opts}
    </select>
  </div>
  <span style="color:#ccc">|</span>
  <div class="ctl">
    <label>Noise Cutoff:</label>
    <select id="noise2d" title="Threshold level as % of displayed max. 0% disables the cutoff.">
      {noise2d_opts}
    </select>
  </div>
  <div class="ctl">
    <label>Cutoff Fade:</label>
    <select id="denoise" title="Off = hard cutoff at the Noise Cutoff threshold. Percent values set transition half-width as % of max. Needs Noise Cutoff &gt; 0 to take effect.">
      {denoise_opts}
    </select>
  </div>
  <span style="color:#ccc">|</span>
  <div class="ctl">
    <label>Drift Smooth:</label>
    <select id="smdrift">
      {smdrift_opts}
    </select>
  </div>
  <div class="ctl">
    <label>Drift Noise:</label>
    <select id="noisedrift">
      {noise_dr_opts}
    </select>
  </div>
  <span style="color:#ccc">|</span>
  <div class="ctl">
    <label>m/z Mode:</label>
    <select id="msmode" title="Sticks: peak picks shown as vertical lines (cleaner). Profile: raw continuous trace.">
      {ms_mode_opts}
    </select>
  </div>
  <div class="ctl">
    <label>m/z Noise:</label>
    <select id="noisemz">
      {noise_mz_opts}
    </select>
  </div>
  <span style="color:#ccc">|</span>
  <div class="ctl">
    <label>m/z Bins:</label>
    <select id="nbins">
      {bins_opts}
    </select>
  </div>
  <div class="ctl">
    <label>Drift Bins:</label>
    <select id="dbins">
      {dbins_opts}
    </select>
  </div>
  <span style="color:#ccc">|</span>
  <div class="ctl">
    <label>Colormap:</label>
    <select id="cmap">
      {cmap_opts}
    </select>
  </div>
  <div class="ctl">
    <label>Render:</label>
    <select id="render" title="Crisp = pixel-sharp; Smooth = bilinear/bicubic interpolation">
      {render_opts}
    </select>
  </div>
  <div class="ctl">
    <label>Intensity:</label>
    <select id="zscale" title="Heatmap intensity transform. log₁₀ compresses dynamic range but accentuates low-level halos; Linear / Sqrt preserve relative peak heights better.">
      {zscale_opts}
    </select>
  </div>
  <span style="color:#ccc">|</span>
  <button id="apply">Apply</button>
</div>
<div id="main" style="display:flex;align-items:flex-start">
  <div id="plot" style="flex:1;min-width:0"></div>
  <div id="exportpanel" style="display:flex;flex-direction:column;gap:8px;padding:12px 14px;
    background:#f8f8f8;border-left:1px solid #ddd;min-width:120px">
    <label style="font-size:{int(FS*0.9)}px;font-weight:700;color:#444">Export</label>
    <button id="expfull" title="Export full plot as high-res PNG"
      style="font-size:{int(FS*0.85)}px;padding:5px 10px;border:1px solid #999;border-radius:3px;
      background:#eee;cursor:pointer">Full PNG</button>
    <button id="exp2d" title="Export 2D heatmap only as high-res PNG"
      style="font-size:{int(FS*0.85)}px;padding:5px 10px;border:1px solid #999;border-radius:3px;
      background:#eee;cursor:pointer">2D PNG</button>
    <button id="expdrift" title="Export drift profile as CSV"
      style="font-size:{int(FS*0.85)}px;padding:5px 10px;border:1px solid #999;border-radius:3px;
      background:#eee;cursor:pointer">Drift CSV</button>
    <button id="expmz" title="Export m/z spectrum as CSV"
      style="font-size:{int(FS*0.85)}px;padding:5px 10px;border:1px solid #999;border-radius:3px;
      background:#eee;cursor:pointer">m/z CSV</button>
  </div>
</div>
<script>
var H={hm_json},MZ={mz_json},DR={dr_json},msMz={ms_mz_j},msI={ms_int_j};
var rawIM={raw_im_j};
// hmMzLo/hmMzHi = the heatmap m/z domain. HM_MZ_LO0/HI0 is the IMMUTABLE baked
// domain; hmMzLo/hmMzHi are kept equal to it (never narrowed by zoom).
var HM_MZ_LO0={hm_mz_range[0]},HM_MZ_HI0={hm_mz_range[1]};
var hmMzLo=HM_MZ_LO0,hmMzHi=HM_MZ_HI0;
var nD=H.length,nM=MZ.length;
var nD_native=nD;
var nM_native={n_mz_native};
var nativeMzStep={native_mz_step};
var pusherUs={pusher_us if pusher_us else 0};
var axRatio=nM/nD;
function ensureOdd(v){{v=Math.round(v);return v<3?3:(v%2===0?v+1:v);}}

// --- Rebin heatmap from raw IM data (both axes) ---
// Reads the CURRENT visible m/z range from the heatmap's x-axis so the same
// nbMz bins serve only what the user is looking at — the dominant fix for
// the "blurry overview" problem identified in the blur diagnosis doc. Zooming
// then picking nbMz spends those bins on the current window, not the full
// data range.
// === Rebuild the heatmap H from immutable rawIM ===
// ORDER-INDEPENDENCE: H is rebuilt over the FIXED baked m/z domain
// [HM_MZ_LO0, HM_MZ_HI0] — NEVER the live zoom. So H is a pure function of
// (nbMz, nbDr): it has no zoom history and cannot ratchet. Zooming changes only
// the visible axis range; the heatmap is already at nbMz resolution across the
// whole domain. hmMzLo/hmMzHi stay constant = the baked domain.
function rebinFull() {{
  var mbVal = document.getElementById('nbins').value;
  var nbMz = (mbVal === 'native') ? nM_native : parseInt(mbVal);
  var dbVal = document.getElementById('dbins').value;
  var nbDr = (dbVal === 'native') ? nD_native : parseInt(dbVal);
  var mzLo = HM_MZ_LO0, mzHi = HM_MZ_HI0;   // immutable baked domain
  hmMzLo = mzLo; hmMzHi = mzHi;
  var mzStep = (mzHi - mzLo) / nbMz;
  var drStep = nD_native / nbDr;
  var newMZ = new Array(nbMz);
  for (var i = 0; i < nbMz; i++) newMZ[i] = mzLo + (i + 0.5) * mzStep;
  var newDR = new Array(nbDr);
  for (var i = 0; i < nbDr; i++) {{
    var binCenter = (i + 0.5) * drStep;
    newDR[i] = pusherUs > 0 ? binCenter * pusherUs / 1000 : binCenter;
  }}
  var newH = [];
  for (var d = 0; d < nbDr; d++) {{
    var row = new Array(nbMz).fill(0);
    var dLo = Math.floor(d * drStep), dHi = Math.floor((d + 1) * drStep);
    if (dHi > nD_native) dHi = nD_native;
    for (var dd = dLo; dd < dHi; dd++) {{
      var pts = rawIM[dd];
      for (var i = 0; i < pts.length; i++) {{
        var bi = Math.floor((pts[i][0] - mzLo) / mzStep);
        if (bi >= 0 && bi < nbMz) row[bi] += pts[i][1];
      }}
    }}
    newH.push(row);
  }}
  H = newH; MZ = newMZ; DR = newDR; nD = nbDr; nM = nbMz;
  axRatio = nM / nD;
  applyAll();
}}
var figData={fig_json};
var P=document.getElementById('plot');
Plotly.newPlot(P,figData.data,figData.layout,{{responsive:true}}).then(function(){{
  // One-shot render so the FIRST paint already reflects the default smoothing /
  // noise / intensity controls. Without this the baked traces are raw and the
  // heatmap visibly "jumps" the first time any control is touched.
  applyAll();
}});

// --- Savitzky-Golay in JS (fits polynomial to sliding window) ---
function sgSmooth1D(y,w,ord){{
  if(w<3||w>y.length||w%2===0)return y.slice();
  var half=Math.floor(w/2),n=y.length,out=new Array(n);
  // Build Vandermonde and pseudo-inverse for the window
  var V=[];for(var i=0;i<w;i++){{var row=[];var x=i-half;var v=1;
    for(var j=0;j<=ord;j++){{row.push(v);v*=x;}}V.push(row);}}
  // Solve V^T V c = V^T y for each window position
  // Precompute (V^T V)^-1 V^T = C (the convolution coefficients)
  var p=ord+1;
  // V^T V
  var VtV=[];for(var i=0;i<p;i++){{VtV.push(new Array(p).fill(0));
    for(var j=0;j<p;j++)for(var k=0;k<w;k++)VtV[i][j]+=V[k][i]*V[k][j];}}
  // Invert VtV (small matrix, use Gauss-Jordan)
  var aug=[];for(var i=0;i<p;i++){{aug.push([]);for(var j=0;j<p;j++)aug[i].push(VtV[i][j]);
    for(var j=0;j<p;j++)aug[i].push(i===j?1:0);}}
  for(var i=0;i<p;i++){{var mx=i;for(var j=i+1;j<p;j++)if(Math.abs(aug[j][i])>Math.abs(aug[mx][i]))mx=j;
    [aug[i],aug[mx]]=[aug[mx],aug[i]];var d=aug[i][i];if(Math.abs(d)<1e-12)return y.slice();
    for(var j=0;j<2*p;j++)aug[i][j]/=d;
    for(var j=0;j<p;j++)if(j!==i){{var f=aug[j][i];for(var k=0;k<2*p;k++)aug[j][k]-=f*aug[i][k];}}}}
  var inv=[];for(var i=0;i<p;i++){{inv.push([]);for(var j=0;j<p;j++)inv[i].push(aug[i][p+j]);}}
  // C = inv * V^T => coeffs[j] = sum_i inv[0][i] * V[j][i]  (we only need c[0] = smoothed value)
  var coeffs=new Array(w);
  for(var j=0;j<w;j++){{var s=0;for(var i=0;i<p;i++)for(var k=0;k<p;k++)s+=inv[0][k]*V[j][k];coeffs[j]=s;}}
  // Apply
  for(var i=0;i<n;i++){{var s=0;for(var j=0;j<w;j++){{var idx=i-half+j;
    if(idx<0)idx=0;if(idx>=n)idx=n-1;s+=coeffs[j]*y[idx];}}out[i]=Math.max(0,s);}}
  return out;
}}
function sgSmooth2D(h,wd,wm,od,om){{
  var nd=h.length,nm=h[0].length;
  // Smooth along drift (axis 0)
  var t1=[];for(var m=0;m<nm;m++){{var col=[];for(var d=0;d<nd;d++)col.push(h[d][m]);
    var sc=(wd>=3&&wd<=nd)?sgSmooth1D(col,wd,od):col;
    for(var d=0;d<nd;d++){{if(!t1[d])t1[d]=[];t1[d][m]=Math.max(0,sc[d]);}}}}
  // Smooth along mz (axis 1)
  var t2=[];for(var d=0;d<nd;d++){{
    var row=t1[d];var sr=(wm>=3&&wm<=nm)?sgSmooth1D(row,wm,om):row;
    t2.push(sr.map(function(v){{return Math.max(0,v)}}));}}
  return t2;
}}
// --- Gaussian smoothing (all-positive kernel, cannot produce negatives) ---
function gaussSmooth1D(y, sigma) {{
  var n = y.length;
  if (sigma <= 0 || n < 3) return y.slice();
  var r = Math.ceil(4 * sigma);
  var kern = [], ksum = 0;
  for (var i = -r; i <= r; i++) {{
    var v = Math.exp(-0.5 * (i / sigma) * (i / sigma));
    kern.push(v); ksum += v;
  }}
  for (var i = 0; i < kern.length; i++) kern[i] /= ksum;
  var out = new Array(n);
  for (var i = 0; i < n; i++) {{
    var s = 0;
    for (var j = 0; j < kern.length; j++) {{
      var idx = i - r + j;
      if (idx < 0) idx = -idx;
      if (idx >= n) idx = 2 * n - 2 - idx;
      if (idx < 0) idx = 0;
      if (idx >= n) idx = n - 1;
      s += kern[j] * y[idx];
    }}
    out[i] = Math.max(0, s);
  }}
  return out;
}}
function gaussSmooth2D(h, sigmaD, sigmaM) {{
  var nd = h.length, nm = h[0].length;
  var t1 = [];
  for (var m = 0; m < nm; m++) {{
    var col = [];
    for (var d = 0; d < nd; d++) col.push(h[d][m]);
    var sc = gaussSmooth1D(col, sigmaD);
    for (var d = 0; d < nd; d++) {{ if (!t1[d]) t1[d] = []; t1[d][m] = sc[d]; }}
  }}
  var t2 = [];
  for (var d = 0; d < nd; d++) t2.push(gaussSmooth1D(t1[d], sigmaM));
  return t2;
}}

// --- Moving average smoothing ---
function maSmooth1D(y, w) {{
  var n = y.length;
  if (w < 2 || n < w) return y.slice();
  var half = Math.floor(w / 2);
  var out = new Array(n);
  for (var i = 0; i < n; i++) {{
    var s = 0, c = 0;
    for (var j = -half; j <= half; j++) {{
      var idx = i + j;
      if (idx < 0) idx = -idx;
      if (idx >= n) idx = 2 * n - 2 - idx;
      if (idx < 0) idx = 0;
      if (idx >= n) idx = n - 1;
      s += y[idx]; c++;
    }}
    out[i] = Math.max(0, s / c);
  }}
  return out;
}}
function maSmooth2D(h, wD, wM) {{
  var nd = h.length, nm = h[0].length;
  var t1 = [];
  for (var m = 0; m < nm; m++) {{
    var col = [];
    for (var d = 0; d < nd; d++) col.push(h[d][m]);
    var sc = maSmooth1D(col, wD);
    for (var d = 0; d < nd; d++) {{ if (!t1[d]) t1[d] = []; t1[d][m] = sc[d]; }}
  }}
  var t2 = [];
  for (var d = 0; d < nd; d++) t2.push(maSmooth1D(t1[d], wM));
  return t2;
}}

// --- Dispatch: route preset object to smoothing function ---
function smooth1D(y, preset) {{
  if (preset.method === 'gaussian') return gaussSmooth1D(y, preset.sigma);
  if (preset.method === 'sg') return sgSmooth1D(y, preset.w, preset.k);
  if (preset.method === 'ma') return maSmooth1D(y, preset.w);
  return y.slice();
}}
function smooth2D(h, preset) {{
  // Compute m/z scale factor from dropdown
  var sm = document.getElementById('mzscale').value;
  var sf = (sm === 'full') ? axRatio
         : (sm === 'sqrt') ? Math.sqrt(axRatio)
         : parseFloat(sm);  // '1' → 1, '2' → 2
  if (preset.method === 'gaussian') {{
    // PHYSICAL-UNIT sigma: preset.sigma is the base m/z smoothing width in Da,
    // and preset.sigma_dr is the drift width in ms. Scale multiplies the m/z
    // width only, so the dropdown again changes the feature anisotropy:
    // <1 = thinner in m/z, >1 = fatter in m/z. We convert to BIN sigma using
    // the CURRENT heatmap geometry so the physical width remains stable when
    // m/z bin count changes.
    var nmHere = h[0] ? h[0].length : nM;
    var binWidthDa = (hmMzHi - hmMzLo) / nmHere;
    var baseMzSig = preset.sigma / binWidthDa;             // Da -> bins
    var mzSig = baseMzSig * sf;                            // Scale affects m/z width
    var drMsPerBin = pusherUs > 0 ? pusherUs / 1000
                   : (DR.length > 1 ? Math.abs(DR[1] - DR[0]) : 1);
    var drSig = (preset.sigma_dr != null)
              ? preset.sigma_dr / drMsPerBin               // ms -> bins
              : baseMzSig;                                 // fallback: unscaled drift width
    return gaussSmooth2D(h, drSig, mzSig);
  }}
  if (preset.method === 'sg') {{
    var wM = ensureOdd(preset.w * sf);
    return sgSmooth2D(h, preset.w, wM, preset.k, preset.k);
  }}
  if (preset.method === 'ma') {{
    var wM = ensureOdd(preset.w * sf);
    return maSmooth2D(h, preset.w, wM);
  }}
  return h;
}}

// === Get visible zoom range ===
function getVisibleRange() {{
  var layout = P.layout;
  var xr = null, yr = null;
  // xaxis2 = heatmap x (shared with xaxis4 = m/z panel x)
  if (layout.xaxis2 && layout.xaxis2.range) xr = layout.xaxis2.range.slice();
  if (layout.yaxis2 && layout.yaxis2.range) yr = layout.yaxis2.range.slice();
  // Fallback to full range. CLAMP the m/z window to the heatmap domain so the
  // m/z marginal can never project rawIM peaks that lie outside what the heatmap
  // shows (e.g. panning past the [350,1250] clip). This keeps the bottom panel
  // and the 2D heatmap consistent.
  var mzLo = xr ? Math.max(hmMzLo, xr[0]) : MZ[0];
  var mzHi = xr ? Math.min(hmMzHi, xr[1]) : MZ[nM-1];
  if (mzHi <= mzLo) {{ mzLo = MZ[0]; mzHi = MZ[nM-1]; }}
  var drLo = yr ? yr[0] : DR[0],  drHi = yr ? yr[1] : DR[nD-1];
  // current heatmap-bin indices in DR (drift)
  var dBl = 0, dBh = nD - 1;
  if (yr) {{
    for (var i = 0; i < nD; i++) if (DR[i] >= drLo) {{ dBl = i; break; }}
    for (var i = nD - 1; i >= 0; i--) if (DR[i] <= drHi) {{ dBh = i; break; }}
  }}
  // m/z bin indices
  var mBl = 0, mBh = nM-1;
  for (var i = 0; i < nM; i++) if (MZ[i] >= mzLo) {{ mBl = i; break; }}
  for (var i = nM-1; i >= 0; i--) if (MZ[i] <= mzHi) {{ mBh = i; break; }}
  // Native drift bin range for rawIM indexing (rawIM has nD_native entries)
  var drStep = nD_native / nD;
  var ndBl = Math.max(0, Math.floor(dBl * drStep));
  var ndBh = Math.min(nD_native - 1, Math.floor((dBh + 1) * drStep) - 1);
  return {{ mzLo:mzLo, mzHi:mzHi, drLo:drLo, drHi:drHi,
           dBl:dBl, dBh:dBh, ndBl:ndBl, ndBh:ndBh, mBl:mBl, mBh:mBh }};
}}

// === Apply 2D heatmap (reads ONLY 2D controls) ===
// Order: smooth → noise threshold → intensity transform.
// SMOOTH FIRST (was: noise-then-smooth) so that thresholding doesn't leave
// soft halos around clipped regions — the blur diagnosis identified the old
// order as a "softened edges around thresholded regions" contributor.
function apply2D() {{
  var preset = JSON.parse(document.getElementById('sm2d').value);
  var ncut = parseFloat(document.getElementById('noise2d').value);
  var cmap = JSON.parse(document.getElementById('cmap').value);
  var renderVal = document.getElementById('render').value;
  var zsmooth = (renderVal === 'off') ? false : renderVal;
  var zmode = document.getElementById('zscale').value;  // "linear"/"sqrt"/"log10"
  // 1. Smooth the raw heatmap first
  var sm = smooth2D(H, preset);
  // 2. Intensity transform → DISPLAY space (this is what the colormap sees).
  var transform, title;
  if (zmode === 'linear')      {{ transform = function(v){{ return Math.max(0, v); }}; title = 'I'; }}
  else if (zmode === 'sqrt')   {{ transform = function(v){{ return Math.sqrt(Math.max(0, v)); }}; title = '√I'; }}
  else if (zmode === 'cbrt')   {{ transform = function(v){{ return Math.cbrt(Math.max(0, v)); }}; title = '∛I'; }}
  else if (zmode === 'p4')     {{ transform = function(v){{ return Math.pow(Math.max(0, v), 0.25); }}; title = 'I^¼'; }}
  else if (zmode === 'p5')     {{ transform = function(v){{ return Math.pow(Math.max(0, v), 0.2); }}; title = 'I^⅕'; }}
  else if (zmode === 'log2')   {{ transform = function(v){{ return Math.log2(Math.max(0, v) + 1); }}; title = 'log₂(I+1)'; }}
  else /* log10 */             {{ transform = function(v){{ return Math.log10(Math.max(0, v) + 1); }}; title = 'log₁₀(I+1)'; }}
  var z = sm.map(function(r) {{ return r.map(transform); }});
  // 3. Cutoff in DISPLAY space (on the transformed values the colormap shows).
  // Two connected but independent controls:
  //    "Noise Cutoff" = threshold level T (% of displayed max — where the cut sits)
  //    "Cutoff Fade"  = off/hard step, or ramp half-width S (% of displayed max)
  // Transfer: weight ramps 0→1 across [T−S, T+S], times the display value.
  //   Fade Off → S = 0 → hard step at T (sharp)
  //   S large → wider gradual fade
  if (ncut > 0) {{
    var denoiseMode = document.getElementById('denoise').value;
    var sPct = (denoiseMode === 'off') ? 0 : parseFloat(denoiseMode);
    var zmax = 0;
    for (var d = 0; d < z.length; d++) for (var m = 0; m < z[d].length; m++)
      if (z[d][m] > zmax) zmax = z[d][m];
    var T = zmax * ncut / 100, S = zmax * (sPct > 0 ? sPct : 0) / 100;
    var dfn;
    if (S <= 0) {{ dfn = function(v){{ return v < T ? 0 : v; }}; }}        // hard step
    else {{ var lo = T - S, span = 2 * S;
      dfn = function(v){{ var w = (v - lo) / span; w = w < 0 ? 0 : (w > 1 ? 1 : w); return v * w; }}; }}
    z = z.map(function(r) {{ return r.map(dfn); }});
  }}
  // Colorbar title rides on the heatmap trace (per-trace colorbar, NOT a
  // layout coloraxis) so it must be set via restyle, not relayout.
  Plotly.restyle(P, {{z:[z], x:[MZ], y:[DR], colorscale:[cmap], zsmooth:[zsmooth],
    'colorbar.title.text':[title]}}, {HM_IDX});
}}

// === Binary search on sorted array ===
function bsearch(arr, val) {{
  var lo = 0, hi = arr.length - 1;
  while (lo <= hi) {{ var mid = (lo + hi) >> 1; if (arr[mid] < val) lo = mid + 1; else hi = mid - 1; }}
  return lo;
}}

// === Update both marginals — ONE Plotly.restyle, no chains ===
function updateMarginals() {{
  var vr = getVisibleRange();

  // Drift: sum heatmap in visible m/z range (linked to m/z zoom)
  var dp = new Array(nD).fill(0);
  for (var d = 0; d < nD; d++) for (var m = vr.mBl; m <= vr.mBh; m++) dp[d] += H[d][m];
  var ncut_d = parseFloat(document.getElementById('noisedrift').value);
  if (ncut_d > 0) {{
    var pk = 0; for (var i = 0; i < nD; i++) if (dp[i] > pk) pk = dp[i];
    var thr = pk * ncut_d / 100;
    for (var i = 0; i < nD; i++) if (dp[i] < thr) dp[i] = 0;
  }}
  var preset = JSON.parse(document.getElementById('smdrift').value);
  if (preset.method !== 'raw') dp = smooth1D(dp, preset);

  // m/z: drift-zoom-aware profile from rawIM (stored at 4dp = 0.0001 Da, finer than the
  // native ~8 mDa Agilent / ~12 mDa Waters TOF sample spacing).
  // ADAPTIVE m/z grid step: chosen so the visible window has at most ~12k points
  // (keeps profile rendering fast at wide zooms while preserving near-native
  // resolution at narrow zooms). Native step floor = 0.0001 Da.
  var _range = vr.mzHi - vr.mzLo;
  var _niceSteps = [0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0];
  var MZ_STEP = _niceSteps[_niceSteps.length - 1];
  for (var _ni = 0; _ni < _niceSteps.length; _ni++) {{
    if (_range / _niceSteps[_ni] <= 12000) {{ MZ_STEP = _niceSteps[_ni]; break; }}
  }}
  // Floor the grid step at the native sampling interval. A FINER grid injects
  // zeros between native samples, so the connected profile line collapses to a
  // comb of spikes (which looked like sticks). With the floor, profile mode
  // draws true smooth peak humps; stick mode peak-picks from the same profile.
  if (MZ_STEP < nativeMzStep) MZ_STEP = nativeMzStep;
  var nMz = Math.max(1, Math.ceil(_range / MZ_STEP) + 1);
  var mx = new Array(nMz), my = new Array(nMz).fill(0);
  for (var i = 0; i < nMz; i++) mx[i] = vr.mzLo + i * MZ_STEP;
  for (var d = vr.ndBl; d <= vr.ndBh; d++) {{
    var pts = rawIM[d];
    if (!pts || pts.length === 0) continue;
    var lo = 0, hi = pts.length - 1, i0 = pts.length;
    while (lo <= hi) {{ var mid = (lo+hi) >> 1;
      if (pts[mid][0] < vr.mzLo) lo = mid + 1; else {{ i0 = mid; hi = mid - 1; }} }}
    for (var i = i0; i < pts.length && pts[i][0] <= vr.mzHi; i++) {{
      var bi = Math.round((pts[i][0] - vr.mzLo) / MZ_STEP);
      if (bi >= 0 && bi < nMz) my[bi] += pts[i][1];
    }}
  }}
  var ncut_m = parseFloat(document.getElementById('noisemz').value);
  if (ncut_m > 0 && my.length > 0) {{
    var pk = 0; for (var i = 0; i < my.length; i++) if (my[i] > pk) pk = my[i];
    var thr = pk * ncut_m / 100;
    my = my.map(function(v) {{ return v < thr ? 0 : v; }});
  }}
  // Peak of the m/z window BEFORE stick conversion — drives the empty-window
  // guard below (a zoom into an m/z gap with no signal otherwise collapses
  // yaxis4's tozero autorange to ~[0,0] and the panel renders blank).
  var _mzPeak = 0; for (var i = 0; i < my.length; i++) if (my[i] > _mzPeak) _mzPeak = my[i];
  // Apply MS display mode: stick (peak-picked) vs profile (raw)
  var msmode = document.getElementById('msmode').value;
  if (msmode === 'stick' && mx.length > 2) {{
    var pk = 0; for (var i = 0; i < my.length; i++) if (my[i] > pk) pk = my[i];
    var thr = pk * 0.01;  // 1% of in-window max — keeps real peaks, drops shoulders
    // Local-max window ≈ one peak WIDTH (~0.08 Da), NOT the inter-peak spacing.
    // (The grid is now floored at the native sampling so the profile is smooth —
    // no comb — so a narrow window finds one stick per real peak without the old
    // over-wide ±0.5 Da window that merged adjacent isotopes into a single stick.)
    var win = Math.max(2, Math.round(0.08 / MZ_STEP));
    var sx = [], sy = [];
    for (var i = win; i < my.length - win; i++) {{
      if (my[i] <= thr) continue;
      var isMax = true;
      for (var j = i - win; j <= i + win; j++) {{
        if (j !== i && my[j] > my[i]) {{ isMax = false; break; }}
      }}
      if (isMax) {{ sx.push(mx[i], mx[i], null); sy.push(0, my[i], null); }}
    }}
    mx = sx; my = sy;
  }}

  // Single restyle — Plotly auto-scales axes (fixedrange removed)
  Plotly.restyle(P, {{x: [dp, mx], y: [DR, my]}}, [{DRIFT_IDX}, {MZ_IDX}]);
  // Empty-window guard: keep the m/z panel visible (with a sane fallback range)
  // instead of collapsing to a degenerate axis. Only relayout on a state CHANGE,
  // and set a flag so our own relayout doesn't re-trigger the handler (no loop).
  setMzPanelGuard(_mzPeak > 0);
}}
var _mzGuardState = null, _guardingY = false;
function setMzPanelGuard(hasSignal) {{
  if (_mzGuardState === hasSignal) return;  // no transition → nothing to do
  _mzGuardState = hasSignal;
  _guardingY = true;
  var p = hasSignal ? Plotly.relayout(P, {{'yaxis4.autorange': true}})
                    : Plotly.relayout(P, {{'yaxis4.range': [0, 1]}});
  p.then(function() {{ _guardingY = false; }});
}}

// === Master apply (called by Apply button) ===
function applyAll() {{
  apply2D();
  updateMarginals();
}}
document.getElementById('apply').onclick = applyAll;

// Wire dropdown changes
document.getElementById('sm2d').onchange = apply2D;
document.getElementById('mzscale').onchange = apply2D;
document.getElementById('noise2d').onchange = apply2D;
document.getElementById('denoise').onchange = apply2D;
document.getElementById('noisedrift').onchange = updateMarginals;
document.getElementById('noisemz').onchange = updateMarginals;
document.getElementById('msmode').onchange = updateMarginals;
document.getElementById('cmap').onchange = apply2D;
document.getElementById('render').onchange = apply2D;
document.getElementById('zscale').onchange = apply2D;
document.getElementById('smdrift').onchange = updateMarginals;
document.getElementById('nbins').onchange = rebinFull;
document.getElementById('dbins').onchange = rebinFull;

// === Export helpers ===
function downloadCSV(filename, csvText) {{
  var blob = new Blob([csvText], {{type:'text/csv'}});
  var a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = filename;
  a.click();
  URL.revokeObjectURL(a.href);
}}
function downloadPNG(filename, dataUrl) {{
  var a = document.createElement('a');
  a.href = dataUrl;
  a.download = filename;
  a.click();
}}
document.getElementById('expfull').onclick = function() {{
  // Export full plot at layout dimensions × 3 DPI
  var lo = P.layout;
  Plotly.toImage(P, {{format:'png', width:lo.width, height:lo.height, scale:3}}).then(function(url) {{
    downloadPNG('{run_name}_full_plot.png', url);
  }});
}};
document.getElementById('exp2d').onclick = function() {{
  // Export 2D heatmap only — clone the heatmap trace into a standalone figure
  // that matches the HTML appearance exactly
  var hm = P.data[{HM_IDX}];
  var lo = P.layout;
  var tmpDiv = document.createElement('div');
  tmpDiv.style.cssText = 'position:absolute;left:-9999px';
  document.body.appendChild(tmpDiv);
  // Clone axis style from the main figure's heatmap axes (xaxis2/yaxis2)
  var xSrc = lo.xaxis2 || {{}};
  var ySrc = lo.yaxis2 || {{}};
  var ax = {{showline:true, linewidth:1, linecolor:'#444', mirror:true,
    showgrid:true, gridwidth:0.5, gridcolor:'#e0e0e0', griddash:'dot',
    minor:{{showgrid:true, gridwidth:0.3, gridcolor:'#eee', griddash:'dot'}}}};
  var lay = {{
    font: lo.font, width:lo.width, height:lo.height,
    margin:{{l:80, r:100, t:20, b:70}},
    template:'plotly_white',
    plot_bgcolor:'white', paper_bgcolor:'white',
    xaxis:Object.assign({{title:'<i>m/z</i>',
      range: xSrc.range ? xSrc.range.slice() : null}}, ax),
    yaxis:Object.assign({{title:'{drift_label}',
      range: ySrc.range ? ySrc.range.slice() : null}}, ax)
  }};
  // Deep-clone the heatmap trace with all current visual state (incl. zsmooth)
  // and the LIVE colorbar title (tracks the current intensity transform, not a
  // hardcoded log10).
  var _cbTitle = (hm.colorbar && hm.colorbar.title && hm.colorbar.title.text)
               ? hm.colorbar.title.text : 'I';
  var trace = {{type:'heatmap', z:hm.z, x:hm.x, y:hm.y,
    colorscale:hm.colorscale, zsmooth:hm.zsmooth,
    colorbar:{{title:{{text:_cbTitle, side:'top'}}, thickness:22}},
    hovertemplate:hm.hovertemplate}};
  Plotly.newPlot(tmpDiv, [trace], lay).then(function() {{
    return Plotly.toImage(tmpDiv, {{format:'png', width:lo.width, height:lo.height, scale:3}});
  }}).then(function(url) {{
    downloadPNG('{run_name}_2d_heatmap.png', url);
    Plotly.purge(tmpDiv);
    document.body.removeChild(tmpDiv);
  }});
}};
document.getElementById('expdrift').onclick = function() {{
  // Export TWO drift CSVs: raw (native bins) + smoothed (0.02ms grid)
  // 1. Raw — native drift bins, full range TIC
  var dp = new Array(nD).fill(0);
  for (var d = 0; d < nD; d++) for (var m = 0; m < nM; m++) dp[d] += H[d][m];
  var lines = ['drift_time_ms,intensity'];
  for (var d = 0; d < nD; d++) lines.push(DR[d].toFixed(4) + ',' + dp[d].toFixed(1));
  downloadCSV('{run_name}_drift_raw.csv', lines.join('\\n'));

  // 2. Smoothed — apply current settings, resample to 0.02ms grid
  var preset = JSON.parse(document.getElementById('smdrift').value);
  var ncut = parseFloat(document.getElementById('noisedrift').value);
  var dpS = dp.slice();
  if (ncut > 0) {{
    var pk = 0; for (var i = 0; i < dpS.length; i++) if (dpS[i] > pk) pk = dpS[i];
    var thr = pk * ncut / 100;
    dpS = dpS.map(function(v) {{ return v < thr ? 0 : v; }});
  }}
  if (preset.method !== 'raw') dpS = smooth1D(dpS, preset);
  // Resample to 0.02ms uniform grid via linear interpolation
  var dtLo = DR[0], dtHi = DR[nD-1], step = 0.02;
  var lines2 = ['drift_time_ms,intensity'];
  for (var t = dtLo; t <= dtHi; t += step) {{
    // Find bracketing bins
    var idx = 0;
    while (idx < nD-1 && DR[idx+1] < t) idx++;
    var v;
    if (idx >= nD-1) v = dpS[nD-1];
    else {{
      var frac = (t - DR[idx]) / (DR[idx+1] - DR[idx]);
      v = dpS[idx] + frac * (dpS[idx+1] - dpS[idx]);
    }}
    lines2.push(t.toFixed(4) + ',' + Math.max(0, v).toFixed(1));
  }}
  downloadCSV('{run_name}_drift_smoothed.csv', lines2.join('\\n'));
}};
document.getElementById('expmz').onclick = function() {{
  // Export TWO m/z CSVs: raw + denoised (noise floor applied)
  var vr = getVisibleRange();
  var iLo = bsearch(msMz, vr.mzLo);
  var iHi = bsearch(msMz, vr.mzHi);
  if (iHi >= msMz.length) iHi = msMz.length - 1;

  // 1. Raw — visible range, no processing
  var lines = ['mz,intensity'];
  for (var i = iLo; i <= iHi; i++) lines.push(msMz[i].toFixed(6) + ',' + msI[i].toFixed(1));
  downloadCSV('{run_name}_mz_raw.csv', lines.join('\\n'));

  // 2. Denoised — noise floor applied
  var mx = msMz.slice(iLo, iHi + 1);
  var my = msI.slice(iLo, iHi + 1);
  var ncut = parseFloat(document.getElementById('noisemz').value);
  if (ncut > 0) {{
    var peak = 0; for (var i = 0; i < my.length; i++) if (my[i] > peak) peak = my[i];
    var thr = peak * ncut / 100;
    my = my.map(function(v) {{ return v < thr ? 0 : v; }});
  }}
  var lines2 = ['mz,intensity'];
  for (var i = 0; i < mx.length; i++) lines2.push(mx[i].toFixed(6) + ',' + my[i].toFixed(1));
  downloadCSV('{run_name}_mz_denoised.csv', lines2.join('\\n'));
}};

// === Recompute marginals on any zoom/pan (rAF-debounced) ===
var _rafId = 0;
P.on('plotly_relayout', function(ed) {{
  if (_guardingY) return;  // ignore the relayout our own empty-window guard fired
  cancelAnimationFrame(_rafId);
  _rafId = requestAnimationFrame(updateMarginals);
}});
</script></body></html>"""

    (out_dir / f"{run_name}_2d_imms.html").write_text(html, encoding="utf-8")

    # Standalone drift profile removed — the 2D viewer's drift panel
    # serves this purpose with zoom-linked m/z filtering.
