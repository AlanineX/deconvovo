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

def _load_plot_config() -> dict:
    """Load imms_plot_config.json from config/ dir. Returns defaults if missing."""
    import sys
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
                  n_mz_bins: int | None = None, pusher_us: float | None = None):
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

    # Load config
    cfg = _load_plot_config()
    dfl = cfg["defaults"]
    prs = cfg["presets"]
    fig_cfg = cfg.get("figure", {})
    FS = fig_cfg.get("font_size", 20)
    FIG_W = fig_cfg.get("width", 1190)
    FIG_H = fig_cfg.get("height", 680)
    if n_mz_bins is None:
        n_mz_bins = dfl.get("mz_bins", 800)

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

    pd.DataFrame({"mz": mz_nz, "drift_bin": drift_nz, "intensity": int_nz}).to_csv(
        out_dir / f"{run_name}_2d_imms.csv", index=False, float_format="%.4f")

    # --- Build heatmap grid ---
    mz_lo, mz_hi = float(mz_nz.min()), float(mz_nz.max())
    hm_mz_range = (mz_lo, mz_hi)
    n_drift = int(drift_nz.max()) + 1

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
        _drift_groups[int(_dr)].append([round(float(_mz), 2), round(float(_it))])
    _raw_im = []
    for _d in range(n_drift):
        pts = _drift_groups.get(_d, [])
        pts.sort(key=lambda p: p[0])  # sort by m/z for fast range filtering
        _raw_im.append(pts)

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

    # m/z spectrum (bottom) — raw MS data (no drift-dependent scaling)
    fig.add_trace(go.Scatter(
        x=ms_mz, y=ms_int, mode="lines",
        line=dict(color="black", width=0.5), showlegend=False,
        hovertemplate="m/z: %{x:.4f}<br>Int: %{y:.0f}<extra></extra>",
    ), row=2, col=2)

    # 2D heatmap (center) — initial: RAW (no smoothing)
    fig.add_trace(go.Heatmap(
        z=np.log10(sub_raw + 1), x=sub_mz, y=drift_bins,
        colorscale=dfl.get("colormap", "Viridis"),
        colorbar=dict(title=dict(text="log₁₀(I+1)", side="top"),
                      len=0.76, y=0.62, thickness=22),
        zsmooth="best",
        hovertemplate="m/z: %{x:.2f}<br>Drift: %{y:.3f}<br>log(I+1): %{z:.2f}<extra></extra>",
    ), row=1, col=2)

    # No Plotly dropdowns — all controls are HTML elements wired to JS
    cmaps = prs.get("colormap", ["Viridis", "Plasma", "Inferno", "Cividis",
                                    "Hot", "Turbo", "Electric", "Bluered"])

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
                     title_text="Intensity", title_font=tf,
                     tickfont=tkfs, **ax, **go_)
    fig.update_yaxes(row=1, col=1, title_text=drift_label,
                     title_font=tf, tickfont=tkf, **ax, **gm,
                     minor=dict(showgrid=True, gridwidth=0.3, gridcolor="#eee", griddash="dot"))
    fig.update_xaxes(row=1, col=1, autorange="reversed",
                     autorangeoptions=dict(include=[0]), tickfont=tkfs, **ax, **go_)
    fig.update_xaxes(row=1, col=2, showticklabels=False, **ax)
    fig.update_yaxes(row=1, col=2, showticklabels=False, **ax)
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
            # Use HTML entity for sigma
            label = label.replace("σ", "&sigma;")
            sel = " selected" if d == default else ""
            opts.append(f"<option value='{val}'{sel}>{label}</option>")
        return "\n      ".join(opts)

    def _simple_options(values, default, fmt=str):
        opts = []
        for v in values:
            sel = " selected" if str(v) == str(default) else ""
            label = {"0": "Off", "sqrt": "&radic;ratio", "full": "ratio"}.get(str(v), fmt(v))
            opts.append(f'<option value="{v}"{sel}>{label}</option>')
        return "\n      ".join(opts)

    sm2d_opts = _smooth_options(prs["smooth_2d"], dfl["smooth_2d"])
    smdrift_opts = _smooth_options(prs["smooth_drift"], dfl["smooth_drift"])
    scale_opts = _simple_options(prs["scale"], dfl["scale"])
    nfmt = lambda v: "Off" if v == 0 else f"{v}%"
    noise2d_opts = _simple_options(prs.get("noise_2d", [0, 0.5, 1, 3, 5]),
                                    dfl.get("noise_2d", 0), nfmt)
    noise_dr_opts = _simple_options(prs.get("noise_drift", [0, 0.5, 1, 3, 5]),
                                     dfl.get("noise_drift", 0), nfmt)
    noise_mz_opts = _simple_options(prs.get("noise_mz", [0, 0.5, 1, 3, 5]),
                                     dfl.get("noise_mz", 0), nfmt)
    bins_opts = _simple_options(prs["mz_bins"], dfl["mz_bins"])
    dbins_opts = _simple_options(
        prs.get("drift_bins", ["native", 50, 100, 200]),
        dfl.get("drift_bins", "native"),
        lambda v: f"{v}" if v != "native" else "Native"
    )
    # Colormap options: resolve ALL names to arrays so Plotly.js doesn't need name lookup
    import plotly.colors as _pc
    cmap_opts_list = []
    for cm in cmaps:
        if isinstance(cm, dict):
            scale = cm["scale"]
            label = cm["name"]
            sel = " selected" if cm["name"] == dfl["colormap"] else ""
        else:
            try:
                scale = _pc.get_colorscale(cm)
            except Exception:
                scale = cm  # fallback to string name
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
<div id="subtitle">{n_points:,} IM points · {len(ms_mz):,} MS points · {f'{pusher_us:.2f} μs/bin' if pusher_us else 'uncalibrated'} ·
  Double-click to reset zoom · 🏠 resets all</div>
<div id="toolbar">
  <div class="ctl">
    <label>2D Smooth:</label>
    <select id="sm2d">
      {sm2d_opts}
    </select>
  </div>
  <div class="ctl">
    <label>Scale:</label>
    <select id="mzscale" title="How m/z axis smoothing scales relative to drift axis">
      {scale_opts}
    </select>
  </div>
  <span style="color:#ccc">|</span>
  <div class="ctl">
    <label>2D Noise:</label>
    <select id="noise2d">
      {noise2d_opts}
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
var hmMzLo={hm_mz_range[0]},hmMzHi={hm_mz_range[1]};
var nD=H.length,nM=MZ.length;
var nD_native=nD;
var pusherUs={pusher_us if pusher_us else 0};
var axRatio=nM/nD;
function ensureOdd(v){{v=Math.round(v);return v<3?3:(v%2===0?v+1:v);}}

// --- Rebin heatmap from raw IM data (both axes) ---
function rebinFull() {{
  var nbMz = parseInt(document.getElementById('nbins').value);
  var dbVal = document.getElementById('dbins').value;
  var nbDr = (dbVal === 'native') ? nD_native : parseInt(dbVal);
  var mzStep = (hmMzHi - hmMzLo) / nbMz;
  var drStep = nD_native / nbDr;
  var newMZ = new Array(nbMz);
  for (var i = 0; i < nbMz; i++) newMZ[i] = hmMzLo + (i + 0.5) * mzStep;
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
        var bi = Math.floor((pts[i][0] - hmMzLo) / mzStep);
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
Plotly.newPlot(P,figData.data,figData.layout,{{responsive:true}});

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
    return gaussSmooth2D(h, preset.sigma, preset.sigma * sf);
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
  // Fallback to full range
  var mzLo = xr ? xr[0] : MZ[0], mzHi = xr ? xr[1] : MZ[nM-1];
  var dLo  = yr ? Math.max(0, Math.floor(yr[0])) : 0;
  var dHi  = yr ? Math.min(nD-1, Math.ceil(yr[1])) : nD-1;
  // m/z bin indices
  var mBl = 0, mBh = nM-1;
  for (var i = 0; i < nM; i++) if (MZ[i] >= mzLo) {{ mBl = i; break; }}
  for (var i = nM-1; i >= 0; i--) if (MZ[i] <= mzHi) {{ mBh = i; break; }}
  return {{ mzLo:mzLo, mzHi:mzHi, dLo:dLo, dHi:dHi, mBl:mBl, mBh:mBh }};
}}

// === Apply 2D heatmap (reads ONLY 2D controls) ===
// Order: noise threshold (on raw) → smooth → log transform
// Noise BEFORE smoothing so the boundary gets blurred (no rectangular artifacts)
function apply2D() {{
  var preset = JSON.parse(document.getElementById('sm2d').value);
  var ncut = parseFloat(document.getElementById('noise2d').value);
  var cmap = JSON.parse(document.getElementById('cmap').value);
  // 1. Noise threshold on raw data
  var src = H;
  if (ncut > 0) {{
    var mx = 0;
    for (var d = 0; d < nD; d++) for (var m = 0; m < nM; m++) if (H[d][m] > mx) mx = H[d][m];
    var thr = mx * ncut / 100;
    src = H.map(function(r) {{ return r.map(function(v) {{ return v < thr ? 0 : v; }}); }});
  }}
  // 2. Smooth (blurs the noise boundary edges)
  var sm = smooth2D(src, preset);
  // 3. Log transform
  var z = sm.map(function(r) {{ return r.map(function(v) {{
    return Math.log10(Math.max(0, v) + 1);
  }}); }});
  Plotly.restyle(P, {{z:[z], x:[MZ], y:[DR], colorscale:[cmap]}}, {HM_IDX});
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

  // m/z: raw _ms.txt in visible range, with noise
  var iLo = bsearch(msMz, vr.mzLo), iHi = bsearch(msMz, vr.mzHi);
  if (iHi >= msMz.length) iHi = msMz.length - 1;
  var mx = msMz.slice(iLo, iHi + 1);
  var my = msI.slice(iLo, iHi + 1);
  var ncut_m = parseFloat(document.getElementById('noisemz').value);
  if (ncut_m > 0 && my.length > 0) {{
    var pk = 0; for (var i = 0; i < my.length; i++) if (my[i] > pk) pk = my[i];
    var thr = pk * ncut_m / 100;
    my = my.map(function(v) {{ return v < thr ? 0 : v; }});
  }}

  // Single restyle — Plotly auto-scales axes (fixedrange removed)
  Plotly.restyle(P, {{x: [dp, mx], y: [DR, my]}}, [{DRIFT_IDX}, {MZ_IDX}]);
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
document.getElementById('noisedrift').onchange = updateMarginals;
document.getElementById('noisemz').onchange = updateMarginals;
document.getElementById('cmap').onchange = apply2D;
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
  // Deep-clone the heatmap trace with all current visual state
  var trace = {{type:'heatmap', z:hm.z, x:hm.x, y:hm.y,
    colorscale:hm.colorscale, zsmooth:'best',
    colorbar:{{title:{{text:'log\u2081\u2080(I+1)', side:'top'}}, thickness:22}},
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
  cancelAnimationFrame(_rafId);
  _rafId = requestAnimationFrame(updateMarginals);
}});
</script></body></html>"""

    (out_dir / f"{run_name}_2d_imms.html").write_text(html, encoding="utf-8")

    # Standalone drift profile removed — the 2D viewer's drift panel
    # serves this purpose with zoom-linked m/z filtering.

