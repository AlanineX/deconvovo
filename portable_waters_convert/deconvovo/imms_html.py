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
            with open(cfg_path, encoding="utf-8") as f:
                return json.load(f)
        except Exception:
            pass
    return default


def plot_im_data(im_file: Path, ms_file: Path | None, run_name: str, out_dir: Path,
                  n_mz_bins: int | str | None = None, pusher_us: float | None = None,
                  config_path: Path | str | None = None):
    """Interactive 2D IM-MS viewer with linked marginals.

    Layout:
        [drift profile (left)]  [2D heatmap (center)]
                                [m/z spectrum (bottom)]

    Key design:
      - m/z panel sums RAW IM on the original MS/IM coordinates, NOT the binned heatmap
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

    # Read the original MS axis before defining the full viewing domain.
    # Zero-only outer coordinates still belong to the acquisition.
    if ms_file and Path(ms_file).exists():
        ms_data = np.atleast_2d(np.loadtxt(str(ms_file)))
        ms_mz, ms_int = ms_data[:, 0], ms_data[:, 1]
    else:
        ms_mz = np.unique(mz_raw)
        ms_int = np.bincount(np.searchsorted(ms_mz, mz_nz), weights=int_nz,
                             minlength=len(ms_mz))

    # The linked spectrum sums IM intensities on the complete original axis.
    # A mismatched MS table adds its coordinates without replacing IM counts.
    profile_mz = np.union1d(ms_mz, mz_raw)
    profile_int = np.bincount(np.searchsorted(profile_mz, mz_nz), weights=int_nz,
                             minlength=len(profile_mz))

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
    # The heatmap normally uses the requested mass window. A window narrower
    # than a Gaussian's support needs surrounding acquired data for filtering;
    # this processing context does not change the visible limits or 1D ROI.
    mz_lo, mz_hi = float(profile_mz[0]), float(profile_mz[-1])
    _iv_pre = cfg.get("initial_view") or {}
    if _iv_pre.get("mz_range"):
        _hm_lo, _hm_hi = _iv_pre["mz_range"]
        hm_mz_range = (max(mz_lo, float(_hm_lo)), min(mz_hi, float(_hm_hi)))
        if hm_mz_range[1] <= hm_mz_range[0]:
            hm_mz_range = (float(_hm_lo), float(_hm_hi))
        _hm_mask = (mz_nz >= hm_mz_range[0]) & (mz_nz <= hm_mz_range[1])
    else:
        hm_mz_range = (mz_lo, mz_hi)
        _hm_mask = slice(None)
    requested_mz_domain = hm_mz_range
    n_drift = int(drift_raw.max()) + 1

    # Native bin count = unique m/z values inside the heatmap m/z range, capped
    # at 25 600 so the embedded heatmap stays under ~5 M cells × 8 bytes ≈ 40 MB
    # JSON. Clipping first means "Native" is native at the visible scale, not
    # diluted by samples outside the view.
    _unique_im_mz = np.unique(mz_nz)
    def _native_count(domain):
        return min(25600, max(200, int(np.searchsorted(_unique_im_mz, domain[1], "right")
            - np.searchsorted(_unique_im_mz, domain[0]))))
    _native_mode = isinstance(n_mz_bins, str) and n_mz_bins.lower() == "native"
    n_mz_native = _native_count(hm_mz_range)
    if _native_mode:
        n_mz_bins = n_mz_native
    _spacing = np.diff(profile_mz)
    _spacing = _spacing[_spacing > 0]
    finest_mz_step = float(_spacing.min()) if len(_spacing) else 0.01
    _preset = dfl.get("smooth_2d", {})
    _sigma = float(_preset.get("sigma", 0)) if _preset.get("method") == "gaussian" else 0
    if _preset.get("method") == "gaussian" and np.isfinite(_sigma) and _sigma > 0 and mz_hi > mz_lo:
        _dr_count = dfl.get("drift_bins", "native")
        _dr_count = n_drift if _dr_count == "native" else int(_dr_count)
        _previous_bins = None
        _span = requested_mz_domain[1] - requested_mz_domain[0]
        while _dr_count > 0:
            if _native_mode:
                n_mz_bins = _native_count(hm_mz_range)
            if n_mz_bins == _previous_bins or n_mz_bins <= 0:
                break
            _previous_bins = n_mz_bins
            _ratio = n_mz_bins / _dr_count
            _scale = dfl.get("scale", "1")
            _factor = _ratio if _scale == "full" else np.sqrt(_ratio) if _scale == "sqrt" else float(_scale)
            if not np.isfinite(_factor) or _factor < 0:
                break
            _support = _sigma * max(_factor, _ratio if _preset.get("sigma_dr") is None else 0)
            if _span >= 8 * _support:
                break
            _padding = max(4 * _support, (n_mz_bins * finest_mz_step - _span) / 2)
            _width = min(mz_hi - mz_lo, _span + 2 * _padding)
            _low = max(mz_lo, min(requested_mz_domain[0] - _padding, mz_hi - _width))
            hm_mz_range = (_low, _low + _width)
    n_mz_native = _native_count(hm_mz_range)
    if n_mz_bins * n_drift > 8_000_000:
        raise ValueError("The requested heatmap is too large. Choose fewer m/z or drift bins.")
    mz_edges = np.linspace(hm_mz_range[0], hm_mz_range[1], n_mz_bins + 1)
    mz_centers = (mz_edges[:-1] + mz_edges[1:]) / 2
    heatmap = np.zeros((n_drift, n_mz_bins), dtype=np.float64)
    # Histogram cells include their lower edge; the final cell also includes
    # the domain's upper edge, matching the browser's rebinning convention.
    bin_idx = np.searchsorted(mz_edges, mz_nz, side="right") - 1
    bin_idx[mz_nz == mz_edges[-1]] = n_mz_bins - 1
    valid = (bin_idx >= 0) & (bin_idx < n_mz_bins)
    for d, bi, intensity in zip(drift_nz[valid], bin_idx[valid], int_nz[valid]):
        heatmap[int(d), int(bi)] += intensity

    # No auto-trim — show full data range, user can zoom interactively
    sub_raw = heatmap
    sub_mz = mz_centers
    bw = (sub_mz[-1] - sub_mz[0]) / len(sub_mz) if len(sub_mz) > 1 else 0
    drift_bins_int = np.arange(n_drift)
    drift_centers = drift_bins_int.astype(float) + 0.5
    if pusher_us:
        drift_bins = drift_centers * pusher_us / 1000  # bin centres in ms
        drift_label = "Drift Time (ms)"
    else:
        drift_bins = drift_centers
        drift_label = "Drift Time (bins)"

    # --- Raw IM data grouped by drift bin (for independent drift profile) ---
    # Analogous to _ms.txt for the m/z panel: raw data, not rebinned heatmap.
    # rawIM[d] = sorted array of [mz, intensity] pairs at drift bin d.
    from collections import defaultdict
    _drift_groups = defaultdict(list)
    for _mz, _dr, _it in zip(mz_nz, drift_nz, int_nz):
        _drift_groups[int(_dr)].append([float(_mz), float(_it)])
    _raw_im = []
    for _d in range(n_drift):
        pts = _drift_groups.get(_d, [])
        pts.sort(key=lambda p: p[0])  # sort by m/z for fast range filtering
        _raw_im.append(pts)

    # --- Raw drift profile ---
    drift_tic_raw = np.bincount(drift_nz[valid], weights=int_nz[valid],
                                minlength=n_drift)

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
        x=drift_tic_raw, y=drift_bins,
        mode="lines+markers" if len(drift_bins) == 1 else "lines",
        marker=dict(color="black", size=5),
        line=dict(color="black", width=0.7, simplify=False), showlegend=False,
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
    _ms_mz_init = profile_mz
    _ms_int_init = profile_int
    if _iv.get("mz_range"):
        _lo, _hi = _iv["mz_range"]
        _mask = (profile_mz >= _lo) & (profile_mz <= _hi)
        _ms_mz_init = profile_mz[_mask]
        _ms_int_init = profile_int[_mask]
    fig.add_trace(go.Scatter(
        x=_ms_mz_init.tolist(), y=_ms_int_init.tolist(),
        mode="lines+markers" if len(_ms_mz_init) == 1 else "lines",
        marker=dict(color="#333", size=5),
        line=dict(color="#333", width=0.6, simplify=False), showlegend=False,
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
        presets = list(presets)
        if not any({k: v for k, v in p.items() if k != "label"} == default
                   for p in presets):
            presets.append({**default, "label": "Custom"})
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
        values = list(values)
        if not any(str(v) == str(default) for v in values):
            values.append(default)
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
    profile_mz_j = _json.dumps(profile_mz.tolist())
    native_im_mz_j = _json.dumps(_unique_im_mz.tolist())
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
    <button id="expdrift" title="Export full-domain drift intensity and a separately smoothed profile at every original drift bin; these are not a current-view capture"
      style="font-size:{int(FS*0.85)}px;padding:5px 10px;border:1px solid #999;border-radius:3px;
      background:#eee;cursor:pointer">Drift raw + smoothed CSV</button>
    <button id="expmz" title="Export the original MS table in the visible m/z window and a separately thresholded copy; these are not the drift-linked displayed trace"
      style="font-size:{int(FS*0.85)}px;padding:5px 10px;border:1px solid #999;border-radius:3px;
      background:#eee;cursor:pointer">m/z raw + denoised CSV</button>
  </div>
</div>
<script>
var H={hm_json},MZ={mz_json},DR={dr_json},msMz={ms_mz_j},msI={ms_int_j};
var DR_PROFILE=DR.slice();
var rawIM={raw_im_j};
var profileMz={profile_mz_j};
var nativeImMz={native_im_mz_j};
var profileIndex=new Map();
for(var i=0;i<profileMz.length;i++)profileIndex.set(profileMz[i],i);
// The initial processing domain can include Gaussian context. Rebinning derives
// context from the fixed requested bounds and controls, never from mouse zoom.
var HM_MZ_LO0={hm_mz_range[0]},HM_MZ_HI0={hm_mz_range[1]};
var REQUEST_MZ_LO0={requested_mz_domain[0]},REQUEST_MZ_HI0={requested_mz_domain[1]};
var RAW_MZ_LO0={mz_lo},RAW_MZ_HI0={mz_hi};
var finestMzStep={finest_mz_step};
var hmMzLo=HM_MZ_LO0,hmMzHi=HM_MZ_HI0;
var nD=H.length,nM=MZ.length;
var nD_native=nD;
var nM_native={n_mz_native};
var pusherUs={pusher_us if pusher_us else 0};
var axRatio=nM/nD;
function ensureOdd(v){{v=Math.round(v);return v<3?3:(v%2===0?v+1:v);}}

// Keep numeric viewing limits separate from Gaussian processing context.
function nativeMzCount(lo,hi) {{
  var right=bsearch(nativeImMz,hi);
  if(right<nativeImMz.length&&nativeImMz[right]===hi)right++;
  return Math.min(25600,Math.max(200,right-bsearch(nativeImMz,lo)));
}}
function processingMzRange(nbMz,nbDr,nativeMode) {{
  var lo=REQUEST_MZ_LO0,hi=REQUEST_MZ_HI0,span=hi-lo;
  var preset=JSON.parse(document.getElementById('sm2d').value),sigma=Number(preset.sigma)||0;
  if(preset.method==='gaussian'&&Number.isFinite(sigma)&&sigma>0&&RAW_MZ_HI0>RAW_MZ_LO0) {{
    var previous=null;
    while(nbDr>0) {{
      if(nativeMode)nbMz=nativeMzCount(lo,hi);
      if(nbMz===previous||nbMz<=0)break;
      previous=nbMz;
      var ratio=nbMz/nbDr,scale=document.getElementById('mzscale').value;
      var factor=scale==='full'?ratio:scale==='sqrt'?Math.sqrt(ratio):parseFloat(scale);
      if(!Number.isFinite(factor)||factor<0)break;
      var support=sigma*Math.max(factor,preset.sigma_dr==null?ratio:0);
      if(span>=8*support)break;
      var padding=Math.max(4*support,(nbMz*finestMzStep-span)/2);
      var width=Math.min(RAW_MZ_HI0-RAW_MZ_LO0,span+2*padding);
      lo=Math.max(RAW_MZ_LO0,Math.min(REQUEST_MZ_LO0-padding,RAW_MZ_HI0-width));hi=lo+width;
    }}
  }}
  if(nativeMode)nbMz=nativeMzCount(lo,hi);
  return [lo,hi,nbMz];
}}

// --- Rebin heatmap from raw IM data (both axes) ---
// === Rebuild the heatmap H from immutable rawIM ===
// H is a pure function of the original requested bounds, bin counts and
// Gaussian settings. Mouse zoom changes visible limits and linked profiles;
// it never ratchets the processing context or replaces acquired coordinates.
function rebinFull() {{
  var mbVal = document.getElementById('nbins').value;
  var nbMz = (mbVal === 'native') ? nM_native : parseInt(mbVal);
  var dbVal = document.getElementById('dbins').value;
  var nbDr = (dbVal === 'native') ? nD_native : parseInt(dbVal);
  if(!(nbMz>0&&nbDr>0)||!Number.isFinite(nbMz*nbDr)||nbMz*nbDr>8000000)
    throw new Error('The requested heatmap is too large. Choose fewer m/z or drift bins.');
  var domain=processingMzRange(nbMz,nbDr,mbVal==='native'),mzLo=domain[0],mzHi=domain[1];
  nbMz=domain[2];
  if(nbMz*nbDr>8000000)throw new Error('The requested heatmap is too large. Choose fewer m/z or drift bins.');
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
    if (d === nbDr - 1 || dHi > nD_native) dHi = nD_native;
    for (var dd = dLo; dd < dHi; dd++) {{
      var pts = rawIM[dd];
      for (var i = 0; i < pts.length; i++) {{
        var bi = Math.floor((pts[i][0] - mzLo) / mzStep);
        if (pts[i][0] === mzHi) bi = nbMz - 1;
        if (bi >= 0 && bi < nbMz) row[bi] += pts[i][1];
      }}
    }}
    newH.push(row);
  }}
  H = newH; MZ = newMZ; DR = newDR; nD = nbDr; nM = nbMz;
  axRatio = nM / nD;
  return applyAll();
}}
var figData={fig_json};
var P=document.getElementById('plot');
window.watersPendingRender=Plotly.newPlot(P,figData.data,figData.layout,{{responsive:true}}).then(function(){{
  // First paint applies both saved bin counts and the filtering context.
  return rebinFull();
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
  for(var j=0;j<w;j++){{var s=0;for(var k=0;k<p;k++)s+=inv[0][k]*V[j][k];coeffs[j]=s;}}
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
  if(!Number.isFinite(r)||n+2*r>8000000)
    throw new Error('The Gaussian smoothing window is too large. Choose less smoothing or fewer bins.');
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
  var rd=Math.ceil(4*sigmaD),rm=Math.ceil(4*sigmaM);
  if(!Number.isFinite(rd)||!Number.isFinite(rm)||
      (rd>=nd&&(nd+2*rd)*nm>8000000)||(rm>=nm&&(nm+2*rm)*nd>8000000))
    throw new Error('The Gaussian smoothing window is too large. Choose less smoothing or fewer bins.');
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
    var drMsPerBin = DR.length > 1 ? Math.abs(DR[1] - DR[0]) : 2 * Math.abs(DR[0]);
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
  // Matched axes can retain stale editable ranges after a relayout; project
  // the data over the ranges Plotly actually rendered for the shared axes.
  var layout = P._fullLayout || P.layout;
  var xr = null, yr = null;
  // xaxis2 = heatmap x (shared with xaxis4 = m/z panel x)
  if (layout.xaxis2 && layout.xaxis2.range) xr = layout.xaxis2.range.slice();
  if (layout.yaxis2 && layout.yaxis2.range) yr = layout.yaxis2.range.slice();
  // Select visible cells by their edges, rather than only their centres:
  // a partially visible coloured cell must still contribute to the projection.
  // An entirely out-of-domain viewport selects no cells, never the full run.
  var mzLo = xr ? Math.max(hmMzLo, Math.min(xr[0], xr[1])) : hmMzLo;
  var mzHi = xr ? Math.min(hmMzHi, Math.max(xr[0], xr[1])) : hmMzHi;
  var driftHi = nD_native * (pusherUs > 0 ? pusherUs / 1000 : 1);
  var drLo = yr ? Math.max(0, Math.min(yr[0], yr[1])) : 0;
  var drHi = yr ? Math.min(driftHi, Math.max(yr[0], yr[1])) : driftHi;
  var mzStep = (hmMzHi - hmMzLo) / nM, drWidth = driftHi / nD;
  var mBl = mzHi > mzLo ? Math.max(0, Math.floor((mzLo - hmMzLo) / mzStep)) : 0;
  var mBh = mzHi > mzLo ? Math.min(nM - 1, Math.ceil((mzHi - hmMzLo) / mzStep) - 1) : -1;
  var dBl = drHi > drLo ? Math.max(0, Math.floor(drLo / drWidth)) : 0;
  var dBh = drHi > drLo ? Math.min(nD - 1, Math.ceil(drHi / drWidth) - 1) : -1;
  // Raw IM drift selection uses physical native-bin edges, independently of
  // the coarse heatmap. Partially visible native bins retain their samples.
  var nativeWidth=pusherUs>0?pusherUs/1000:1;
  function nativePosition(value) {{
    var position=value/nativeWidth,nearest=Math.round(position);
    return Math.abs(position-nearest)<=2*Number.EPSILON*Math.max(1,Math.abs(position))
      ?nearest:position;
  }}
  var ndBl=drHi>drLo?Math.max(0,Math.floor(nativePosition(drLo))):0;
  var ndBh=drHi>drLo?Math.min(nD_native-1,Math.ceil(nativePosition(drHi))-1):-1;
  return {{ mzLo:mzLo, mzHi:mzHi, drLo:drLo, drHi:drHi,
           dBl:dBl, dBh:dBh, ndBl:ndBl, ndBh:ndBh, mBl:mBl, mBh:mBh }};
}}

// === Apply 2D heatmap (reads ONLY 2D controls) ===
// Order: smooth → intensity transform → display-space noise cutoff/fade.
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
  return Plotly.restyle(P, {{z:[z], x:[MZ], y:[DR], colorscale:[cmap], zsmooth:[zsmooth],
    'colorbar.title.text':[title]}}, {HM_IDX});
}}

// === Binary search on sorted array ===
function bsearch(arr, val) {{
  var lo = 0, hi = arr.length - 1;
  while (lo <= hi) {{ var mid = (lo + hi) >> 1; if (arr[mid] < val) lo = mid + 1; else hi = mid - 1; }}
  return lo;
}}

// Sum exact RAW m/z points on every original drift bin, without 2D rebinning.
function rawDriftProfile(mzLo,mzHi) {{
  var dp=new Array(nD_native).fill(0);
  if(mzHi<mzLo)return dp;
  for(var d=0;d<nD_native;d++) {{
    var pts=rawIM[d],lo=0,hi=pts.length-1,i0=pts.length;
    while(lo<=hi) {{ var mid=(lo+hi)>>1;
      if(pts[mid][0]<mzLo)lo=mid+1;else {{i0=mid;hi=mid-1;}} }}
    for(var i=i0;i<pts.length&&pts[i][0]<=mzHi;i++)dp[d]+=pts[i][1];
  }}
  return dp;
}}

// === Update both marginals — ONE Plotly.restyle, no chains ===
function updateMarginals() {{
  var vr = getVisibleRange();

  // Drift: exact visible m/z selection at all original acquisition drift bins.
  var dp=rawDriftProfile(vr.mzLo,vr.mzHi);
  var ncut_d = parseFloat(document.getElementById('noisedrift').value);
  if (ncut_d > 0) {{
    var pk = 0; for (var i = 0; i < dp.length; i++) if (dp[i] > pk) pk = dp[i];
    var thr = pk * ncut_d / 100;
    for (var i = 0; i < dp.length; i++) if (dp[i] < thr) dp[i] = 0;
  }}
  var preset = JSON.parse(document.getElementById('smdrift').value);
  if (preset.method !== 'raw') dp = smooth1D(dp, preset);

  // Exact acquired m/z coordinates, including the original zero-gap samples.
  // Zooming selects coordinates and drift bins; it never resamples the profile.
  var iLo=bsearch(profileMz,vr.mzLo),iHi=bsearch(profileMz,vr.mzHi);
  if(iHi<profileMz.length&&profileMz[iHi]===vr.mzHi)iHi++;
  if(vr.mzHi<vr.mzLo)iHi=iLo;
  var mx=profileMz.slice(iLo,iHi),my=new Array(mx.length).fill(0);
  for (var d = vr.ndBl; d <= vr.ndBh; d++) {{
    var pts = rawIM[d];
    if (!pts || pts.length === 0) continue;
    var lo = 0, hi = pts.length - 1, i0 = pts.length;
    while (lo <= hi) {{ var mid = (lo+hi) >> 1;
      if (pts[mid][0] < vr.mzLo) lo = mid + 1; else {{ i0 = mid; hi = mid - 1; }} }}
    for (var i = i0; i < pts.length && pts[i][0] <= vr.mzHi; i++) {{
      var bi=profileIndex.get(pts[i][0])-iLo;
      if(bi>=0&&bi<my.length)my[bi]+=pts[i][1];
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
  if (msmode === 'stick') {{
    var thr=_mzPeak*0.01,sx=[],sy=[];
    // Find local maxima on the original coordinates; the midpoint of a flat
    // plateau matches scipy.find_peaks. Endpoints are not local maxima.
    for(var i=1;i<my.length-1;i++) {{
      if(my[i]<=my[i-1])continue;
      var right=i;
      while(right+1<my.length&&my[right+1]===my[i])right++;
      if(right<my.length-1&&my[right]>my[right+1]) {{
        var peak=Math.floor((i+right)/2),isMax=my[peak]>thr;
        var left=bsearch(mx,mx[peak]-0.08);
        for(var j=left;isMax&&j<mx.length&&mx[j]<=mx[peak]+0.08;j++)
          if(my[j]>my[peak])isMax=false;
        if(isMax) {{ sx.push(mx[peak],mx[peak],null);sy.push(0,my[peak],null); }}
      }}
      i=right;
    }}
    mx=sx;my=sy;
  }}

  // Single restyle — Plotly auto-scales axes (fixedrange removed)
  var render=Plotly.restyle(P, {{x: [dp, mx], y: [DR_PROFILE, my],
    mode:[dp.length===1?'lines+markers':'lines',mx.length===1?'lines+markers':'lines']}}, [{DRIFT_IDX}, {MZ_IDX}]);
  // Empty-window guard: keep the m/z panel visible (with a sane fallback range)
  // instead of collapsing to a degenerate axis. Only relayout on a state CHANGE,
  // and set a flag so our own relayout doesn't re-trigger the handler (no loop).
  return Promise.all([render,setMzPanelGuard(_mzPeak > 0)]);
}}
var _mzGuardState = null, _guardingY = false;
function setMzPanelGuard(hasSignal) {{
  if (_mzGuardState === hasSignal) return;  // no transition → nothing to do
  _mzGuardState = hasSignal;
  _guardingY = true;
  var p = hasSignal ? Plotly.relayout(P, {{'yaxis4.autorange': true}})
                    : Plotly.relayout(P, {{'yaxis4.range': [0, 1]}});
  return p.then(function() {{ _guardingY = false; }});
}}

// === Master apply (called by Apply button) ===
function applyAll() {{
  return Promise.all([apply2D(),updateMarginals()]);
}}
document.getElementById('apply').onclick = applyAll;

// Wire dropdown changes
document.getElementById('sm2d').onchange = rebinFull;
document.getElementById('mzscale').onchange = rebinFull;
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
  var blob = new Blob(Array.isArray(csvText)?csvText:[csvText], {{type:'text/csv'}});
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
// Current-view CSVs use the rendered traces, never a second raw-data pipeline.
function watersCurrentCSV(kind) {{
  var vr=getVisibleRange(), lines=[], trace, coordinate, value;
  var driftColumn=pusherUs>0?'drift_time_ms':'drift_bin';
  if(kind==='heatmap_csv') {{
    trace=P.data[{HM_IDX}];
    var chunks=['mz,'+driftColumn+',display_value,intensity_transform\\n'];
    var transform=document.getElementById('zscale').value;
    // Keep a partially visible cell with its original centre and value.
    for(var d=vr.dBl;d<=vr.dBh;d++){{
      lines=[];
      for(var m=vr.mBl;m<=vr.mBh;m++)
        lines.push([trace.x[m],trace.y[d],trace.z[d][m],transform].join(','));
      if(lines.length)chunks.push(lines.join('\\n')+'\\n');
    }}
    return chunks; // bound the temporary line list to one drift row
  }} else if(kind==='mz_csv'||kind==='drift_csv') {{
    var drift=kind==='drift_csv';trace=P.data[drift?{DRIFT_IDX}:{MZ_IDX}];
    coordinate=drift?trace.y:trace.x;value=drift?trace.x:trace.y;
    var lo=drift?vr.drLo:vr.mzLo,hi=drift?vr.drHi:vr.mzHi;
    lines.push((drift?driftColumn:'mz')+',display_value');
    for(var i=0;i<coordinate.length;i++) {{
      // Blank rows preserve the breaks in a stick-mode line trace.
      if(coordinate[i]===null||value[i]===null){{lines.push(',');continue;}}
      if(coordinate[i]>=lo&&coordinate[i]<=hi)
        lines.push(coordinate[i]+','+value[i]);
    }}
  }} else throw new Error('Unknown current-view CSV');
  return lines.join('\\n');
}}
async function watersExportCurrent(kind) {{
  await Promise.resolve(window.watersPendingRender);
  if(window.watersConsoleConfigError)throw new Error(window.watersConsoleConfigError);
  // Finish a zoom's pending linked profiles before reading the live traces.
  cancelAnimationFrame(_rafId);
  await updateMarginals();
  await Plotly.redraw(P);
  await new Promise(function(resolve){{requestAnimationFrame(resolve);}});
  var prefix={_json.dumps(run_name)}+'_current_';
  if(kind.endsWith('_csv')) {{
    downloadCSV(prefix+kind.replace('_csv','.csv'),watersCurrentCSV(kind));
    return;
  }}
  var lo=P.layout, full=P._fullLayout;
  var image={{format:'png',width:full.width,height:full.height,scale:3}};
  if(kind==='full_png'){{downloadPNG(prefix+'full.png',await Plotly.toImage(P,image));return;}}
  var axes, index;
  if(kind==='heatmap_png'){{index={HM_IDX};axes=['xaxis2','yaxis2'];}}
  else if(kind==='mz_png'){{index={MZ_IDX};axes=['xaxis4','yaxis4'];}}
  else if(kind==='drift_png'){{index={DRIFT_IDX};axes=['xaxis','yaxis'];}}
  else throw new Error('Unknown current-view PNG');
  function axis(name,anchor){{
    return Object.assign({{}},lo[name],{{domain:[0,1],anchor:anchor,matches:null,
      range:full[name].range.slice(),autorange:false,showticklabels:true}});
  }}
  var trace=Object.assign({{}},P.data[index],{{xaxis:'x',yaxis:'y'}});
  if(index!=={HM_IDX})trace.line=Object.assign({{}},trace.line,{{simplify:false}});
  if(index==={HM_IDX})trace.colorbar=Object.assign({{}},trace.colorbar,{{len:1,y:0.5}});
  var layout={{font:lo.font,width:full.width,height:full.height,
    margin:{{l:80,r:index==={HM_IDX}?100:40,t:20,b:70}},
    plot_bgcolor:lo.plot_bgcolor,paper_bgcolor:lo.paper_bgcolor,
    xaxis:axis(axes[0],'y'),yaxis:axis(axes[1],'x')}};
  if(index==={HM_IDX}){{
    layout.xaxis.title=Object.assign({{}},layout.xaxis.title,{{text:'<i>m/z</i>'}});
    layout.yaxis.title=Object.assign({{}},layout.yaxis.title,{{text:'{drift_label}'}});
  }}
  if(index==={DRIFT_IDX})layout.xaxis.title=Object.assign({{}},layout.xaxis.title,{{text:'Intensity'}});
  var temporary=document.createElement('div');
  temporary.style.cssText='position:absolute;left:-9999px';document.body.appendChild(temporary);
  try{{
    await Plotly.newPlot(temporary,[trace],layout,{{staticPlot:true}});
    downloadPNG(prefix+kind.replace('_png','.png'),await Plotly.toImage(temporary,image));
  }}finally{{Plotly.purge(temporary);temporary.remove();}}
}}
document.getElementById('expfull').onclick = function() {{watersExportCurrent('full_png');}};
document.getElementById('exp2d').onclick = function() {{watersExportCurrent('heatmap_png');}};
document.getElementById('expdrift').onclick = function() {{
  // Both legacy CSVs retain every original coordinate over the baked domain.
  var dp=rawDriftProfile(REQUEST_MZ_LO0,REQUEST_MZ_HI0);
  var driftColumn=pusherUs>0?'drift_time_ms':'drift_bin';
  var lines=[driftColumn+',intensity'];
  for(var d=0;d<nD_native;d++)lines.push(DR_PROFILE[d]+','+dp[d]);
  downloadCSV('{run_name}_drift_raw.csv',lines.join('\\n'));

  var preset=JSON.parse(document.getElementById('smdrift').value);
  var ncut=parseFloat(document.getElementById('noisedrift').value),dpS=dp.slice();
  if(ncut>0) {{
    var pk=0;for(var i=0;i<dpS.length;i++)if(dpS[i]>pk)pk=dpS[i];
    var thr=pk*ncut/100;
    dpS=dpS.map(function(v) {{return v<thr?0:v;}});
  }}
  if(preset.method!=='raw')dpS=smooth1D(dpS,preset);
  var lines2=[driftColumn+',intensity'];
  for(var d=0;d<nD_native;d++)lines2.push(DR_PROFILE[d]+','+Math.max(0,dpS[d]));
  downloadCSV('{run_name}_drift_smoothed.csv',lines2.join('\\n'));
}};
document.getElementById('expmz').onclick = function() {{
  // Export TWO m/z CSVs: raw + denoised (noise floor applied)
  var vr = getVisibleRange();
  var iLo = bsearch(msMz, vr.mzLo);
  var iHi = bsearch(msMz, vr.mzHi);
  if (iHi < msMz.length && msMz[iHi] === vr.mzHi) iHi++;
  if (vr.mzHi < vr.mzLo) iHi = iLo;

  // 1. Raw — visible range, no processing
  var lines = ['mz,intensity'];
  for (var i = iLo; i < iHi; i++) lines.push(msMz[i] + ',' + msI[i]);
  downloadCSV('{run_name}_mz_raw.csv', lines.join('\\n'));

  // 2. Denoised — noise floor applied
  var mx = msMz.slice(iLo, iHi);
  var my = msI.slice(iLo, iHi);
  var ncut = parseFloat(document.getElementById('noisemz').value);
  if (ncut > 0) {{
    var peak = 0; for (var i = 0; i < my.length; i++) if (my[i] > peak) peak = my[i];
    var thr = peak * ncut / 100;
    my = my.map(function(v) {{ return v < thr ? 0 : v; }});
  }}
  var lines2 = ['mz,intensity'];
  for (var i = 0; i < mx.length; i++) lines2.push(mx[i] + ',' + my[i]);
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

    # Embed the JavaScript matching the installed Plotly Python package so the
    # viewer renders typed arrays correctly and opens without internet access.
    from plotly.offline import get_plotlyjs
    html = html.replace('<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>',
                        '<script>' + get_plotlyjs() + '</script>')
    (out_dir / f"{run_name}_2d_imms.html").write_text(html, encoding="utf-8")

    # Standalone drift profile removed — the 2D viewer's drift panel
    # serves this purpose with zoom-linked m/z filtering.
