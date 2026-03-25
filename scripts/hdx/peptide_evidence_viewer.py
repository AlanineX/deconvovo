#!/usr/bin/env python3
"""Interactive peptide evidence viewer for HDX-MS discovery results.

Features:
- Clickable TIC chromatogram — click to load that scan's spectrum
- Full MS1 spectrum — zoomable, shows peptide markers on selection
- Sortable peptide table with full sequences
- Isotopic envelope detail with ppm errors

Usage:
    python scripts/hdx/peptide_evidence_viewer.py \
        --mzml output/01_hdx_mzml/01_20260217_TTR_peptide_mapping.mzML \
        --peptides output/03_hdx_peptide_discovery/discovered_peptides.csv \
        --out output/03_hdx_peptide_discovery/evidence_viewer.html
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from pyteomics import mzml

PROTON = 1.007276466812
NEUTRON = 1.0033548378


def _downsample(mz, intensity, n=5000):
    if len(mz) <= n:
        return mz, intensity
    top = np.argsort(intensity)[-n:]
    top = np.sort(top)
    return mz[top], intensity[top]


def extract_evidence(mzml_path: Path, peptides: pd.DataFrame,
                     ppm: float = 15.0, n_isotopes: int = 5,
                     n_explore_scans: int = 50) -> tuple[list[dict], dict]:
    """Extract per-peptide evidence + stored scans for TIC navigation."""
    # Pass 1: collect all scan metadata
    scan_meta = []  # [{scan, rt, tic, mz, intensity}, ...]
    print("  Reading MS1 scans...")
    for spec in mzml.MzML(str(mzml_path)):
        if spec.get("ms level", 1) != 1:
            continue
        mz_arr = spec.get("m/z array", np.array([]))
        int_arr = spec.get("intensity array", np.array([]))
        rt = spec.get("scanList", {}).get("scan", [{}])[0].get("scan start time", 0)
        scan_meta.append({
            "scan": len(scan_meta) + 1,
            "rt": float(rt),
            "tic": float(int_arr.sum()) if len(int_arr) > 0 else 0,
            "mz": mz_arr, "intensity": int_arr,
        })
    print(f"    {len(scan_meta)} MS1 scans")

    # Determine which scans to store
    peptide_scans = set(int(s) for s in peptides["scan"].unique() if s <= len(scan_meta))
    apex_scan = max(scan_meta, key=lambda s: s["tic"])["scan"]
    step = max(1, len(scan_meta) // n_explore_scans)
    explore_scans = set(range(1, len(scan_meta) + 1, step))
    store_scans = peptide_scans | explore_scans | {apex_scan}
    print(f"    Storing {len(store_scans)} scans ({len(peptide_scans)} peptide + {len(explore_scans)} explore)")

    # Build stored scan spectra (downsampled)
    scan_spectra = {}
    for sm in scan_meta:
        if sm["scan"] in store_scans:
            ds_mz, ds_int = _downsample(sm["mz"], sm["intensity"])
            scan_spectra[sm["scan"]] = {
                "mz": [round(float(m), 4) for m in ds_mz],
                "int": [round(float(i), 1) for i in ds_int],
            }

    # TIC arrays
    tic_scans = [s["scan"] for s in scan_meta]
    tic_rt = [round(s["rt"], 3) for s in scan_meta]
    tic_int = [round(s["tic"], 0) for s in scan_meta]
    stored_list = sorted(store_scans)

    global_data = {
        "tic_scans": tic_scans, "tic_rt": tic_rt, "tic_int": tic_int,
        "stored_scans": stored_list,
        "scan_spectra": scan_spectra,
        "apex_scan": apex_scan,
    }

    # Pass 2: per-peptide evidence
    evidence = []
    scan_lookup = {s["scan"]: s for s in scan_meta}

    for _, row in peptides.iterrows():
        scan_id = int(row["scan"])
        sm = scan_lookup.get(scan_id)
        if sm is None:
            continue
        mz_arr, int_arr = sm["mz"], sm["intensity"]
        z = int(row["charge"])
        mz_mono = float(row["mz_mono"])

        # Region for envelope plot
        mz_lo = mz_mono - 1.5
        mz_hi = mz_mono + (n_isotopes + 1) * NEUTRON / z + 1.5
        mask = (mz_arr >= mz_lo) & (mz_arr <= mz_hi)

        matched = []
        for iso in range(n_isotopes):
            mz_theo = mz_mono + iso * NEUTRON / z
            tol = mz_theo * ppm / 1e6
            pmask = np.abs(mz_arr - mz_theo) <= tol
            if pmask.any():
                bi = np.argmax(int_arr * pmask)
                matched.append({"iso": iso, "mz_t": round(mz_theo, 4),
                    "mz_o": round(float(mz_arr[bi]), 4),
                    "i": round(float(int_arr[bi]), 1),
                    "ppm": round((float(mz_arr[bi]) - mz_theo) / mz_theo * 1e6, 2),
                    "ok": True})
            else:
                matched.append({"iso": iso, "mz_t": round(mz_theo, 4),
                    "mz_o": None, "i": 0, "ppm": None, "ok": False})

        n_ok = sum(1 for p in matched if p["ok"])
        avg_ppm = np.mean([p["ppm"] for p in matched if p["ok"]])

        evidence.append({
            "seq": row["sequence"],
            "s": int(row["start"]), "e": int(row["end"]),
            "z": z, "mz": round(mz_mono, 4),
            "sc": round(float(row["score"]), 1),
            "ni": n_ok, "scan": scan_id,
            "ppm": round(float(avg_ppm), 2) if not np.isnan(avg_ppm) else None,
            "rmz": [round(float(m), 4) for m in mz_arr[mask]],
            "ri": [round(float(i), 1) for i in int_arr[mask]],
            "m": matched,
        })

    print(f"    {len(evidence)} peptides with evidence")
    return evidence, global_data


def build_html(evidence, gd, title="Peptide Evidence Viewer"):
    ev_j = json.dumps(evidence)
    gd_j = json.dumps(gd)
    return f"""<!DOCTYPE html>
<html><head><meta charset="utf-8"><title>{title}</title>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<style>
*{{box-sizing:border-box}}
body{{font-family:'Segoe UI',Arial,sans-serif;margin:0;padding:10px;background:#f8f8fc}}
h2{{font-size:17px;margin:0 0 8px;color:#333}}
.row{{display:flex;gap:10px;margin-bottom:10px}}
.plot{{background:#fff;border:1px solid #ddd;border-radius:6px}}
#tic{{flex:1;height:200px}} #fullspec{{flex:1;height:200px}}
#main{{display:flex;gap:10px}}
#tbl-wrap{{width:460px;flex-shrink:0}}
#tbl-scroll{{max-height:50vh;overflow-y:auto;border:1px solid #ddd;border-radius:6px;background:#fff}}
table{{width:100%;border-collapse:collapse;font-size:11px;cursor:pointer}}
th{{background:#e0e0ee;padding:5px 6px;text-align:left;position:sticky;top:0;cursor:pointer;user-select:none;white-space:nowrap}}
th:hover{{background:#c8c8dd}} th.asc::after{{content:' \\25B2'}} th.desc::after{{content:' \\25BC'}}
td{{padding:3px 6px;border-bottom:1px solid #eee;white-space:nowrap}}
tr:hover{{background:#e8f0ff}} tr.sel{{background:#b8d0ff}} tr.scan-hl{{background:#fff8e0}}
#rhs{{flex:1;display:flex;flex-direction:column;gap:10px}}
#env{{height:320px}}
#pinfo{{font-size:11px}}
#pinfo td{{padding:2px 5px}}
.ok{{color:#2e7d32}} .no{{color:#c62828}}
#status{{font-size:11px;color:#666;padding:4px 0}}
.info{{font-size:11px;color:#888;padding:2px 0}}
</style></head><body>
<h2>{title}</h2>
<div id="status"></div>

<div class="row">
  <div id="tic" class="plot"></div>
  <div id="fullspec" class="plot"></div>
</div>

<div id="main">
<div id="tbl-wrap">
  <div class="info">Click headers to sort. Click row to view envelope. Click TIC to change scan.</div>
  <div id="tbl-scroll">
    <table><thead><tr>
      <th data-c="seq">Sequence</th><th data-c="s">Start</th><th data-c="e">End</th>
      <th data-c="z">z</th><th data-c="mz">m/z</th><th data-c="sc">Score</th>
      <th data-c="ni">Iso</th><th data-c="ppm">ppm</th><th data-c="scan">Scan</th>
    </tr></thead><tbody id="tb"></tbody></table>
  </div>
</div>
<div id="rhs">
  <div id="env" class="plot"></div>
  <div id="pinfo"></div>
</div>
</div>

<script>
var D={ev_j};
var G={gd_j};
var N={NEUTRON};
var curScan=G.apex_scan, sortC=null, sortA=true;
var ord=D.map(function(_,i){{return i}});

// ---- TIC ----
var ticTrace={{x:G.tic_rt, y:G.tic_int, type:'scatter', mode:'lines',
  line:{{color:'#4a9eff',width:1}}, hovertemplate:'RT:%{{x:.2f}} min<br>TIC:%{{y:.0f}}<extra></extra>'}};
var ticLine={{x:[0,0],y:[0,1],type:'scatter',mode:'lines',
  line:{{color:'red',width:1.5,dash:'dot'}},showlegend:false,yaxis:'y'}};
Plotly.newPlot('tic',[ticTrace,ticLine],{{
  margin:{{t:25,b:30,l:50,r:10}},title:'TIC Chromatogram (click to select scan)',
  xaxis:{{title:'RT (min)'}},yaxis:{{title:'TIC'}},
}},{{responsive:true}});

document.getElementById('tic').on('plotly_click',function(ed){{
  if(!ed||!ed.points||!ed.points[0])return;
  var pi=ed.points[0].pointIndex;
  var sn=G.tic_scans[pi];
  // Snap to nearest stored scan
  var best=G.stored_scans[0],bd=Math.abs(sn-best);
  G.stored_scans.forEach(function(s){{var d=Math.abs(sn-s);if(d<bd){{best=s;bd=d}}}});
  selectScan(best);
}});

function selectScan(sn){{
  curScan=sn;
  // Move TIC line
  var idx=G.tic_scans.indexOf(sn);
  if(idx<0)idx=0;
  var rt=G.tic_rt[idx];
  var ymax=Math.max.apply(null,G.tic_int);
  Plotly.restyle('tic',{{x:[[rt,rt]],y:[[0,ymax]]}},1);

  // Load spectrum
  var sp=G.scan_spectra[sn];
  if(sp){{
    Plotly.react('fullspec',[{{
      x:sp.mz,y:sp['int'],type:'scatter',mode:'lines',
      line:{{color:'#555',width:0.5}},
      hovertemplate:'m/z:%{{x:.4f}}<br>Int:%{{y:.0f}}<extra></extra>',
    }}],{{
      margin:{{t:25,b:30,l:50,r:10}},
      title:'MS1 Scan '+sn+' (RT:'+G.tic_rt[G.tic_scans.indexOf(sn)].toFixed(2)+' min)',
      xaxis:{{title:'m/z'}},yaxis:{{title:'Intensity'}},
    }},{{responsive:true}});
  }}

  // Highlight table rows matching this scan
  document.querySelectorAll('tr.scan-hl').forEach(function(r){{r.classList.remove('scan-hl')}});
  document.querySelectorAll('tr[data-scan="'+sn+'"]').forEach(function(r){{r.classList.add('scan-hl')}});

  var nAtScan=D.filter(function(d){{return d.scan===sn}}).length;
  document.getElementById('status').textContent=
    D.length+' peptides | Scan '+sn+' | RT '+G.tic_rt[G.tic_scans.indexOf(sn)].toFixed(2)+
    ' min | '+nAtScan+' peptides at this scan';
}}

// ---- Table ----
var tb=document.getElementById('tb');
function renderTable(){{
  tb.innerHTML='';
  ord.forEach(function(i){{
    var d=D[i],tr=document.createElement('tr');
    tr.dataset.idx=i; tr.dataset.scan=d.scan;
    tr.innerHTML='<td>'+d.seq+'</td><td>'+d.s+'</td><td>'+d.e+'</td><td>'+d.z+
      '</td><td>'+d.mz.toFixed(2)+'</td><td>'+d.sc.toFixed(0)+
      '</td><td>'+d.ni+'</td><td>'+(d.ppm!==null?d.ppm.toFixed(1):'-')+
      '</td><td>'+d.scan+'</td>';
    tr.onclick=function(){{showPeptide(i)}};
    if(d.scan===curScan)tr.classList.add('scan-hl');
    tb.appendChild(tr);
  }});
}}
document.querySelectorAll('th[data-c]').forEach(function(th){{
  th.onclick=function(){{
    var c=this.dataset.c;
    if(sortC===c)sortA=!sortA; else{{sortC=c;sortA=true}}
    document.querySelectorAll('th').forEach(function(h){{h.className=''}});
    this.className=sortA?'asc':'desc';
    ord.sort(function(a,b){{
      var va=D[a][c],vb=D[b][c];
      if(va===null)return 1;if(vb===null)return -1;
      if(typeof va==='string')return sortA?va.localeCompare(vb):vb.localeCompare(va);
      return sortA?va-vb:vb-va;
    }});
    renderTable();
  }};
}});
renderTable();

// ---- Envelope ----
function showPeptide(idx){{
  document.querySelectorAll('tr.sel').forEach(function(r){{r.classList.remove('sel')}});
  var tr=document.querySelector('tr[data-idx="'+idx+'"]');
  if(tr)tr.classList.add('sel');
  var d=D[idx];

  // Mark on full spectrum
  var sp=G.scan_spectra[d.scan];
  if(sp){{
    var markers=d.m.map(function(p){{return p.mz_t}});
    var mi=d.m.map(function(p){{return p.ok?p.i:0}});
    Plotly.react('fullspec',[
      {{x:sp.mz,y:sp['int'],type:'scatter',mode:'lines',line:{{color:'#555',width:0.5}},
        hovertemplate:'m/z:%{{x:.4f}}<br>Int:%{{y:.0f}}<extra></extra>',name:'Spectrum'}},
      {{x:markers,y:mi,type:'scatter',mode:'markers',
        marker:{{color:'red',size:8,symbol:'triangle-down'}},
        name:d.seq+' z='+d.z+'+'}},
    ],{{
      margin:{{t:25,b:30,l:50,r:10}},
      title:'MS1 Scan '+d.scan,
      xaxis:{{title:'m/z',range:[d.mz-3,d.mz+8*N/d.z+3]}},
      yaxis:{{title:'Intensity'}},
      showlegend:true,legend:{{x:0.65,y:1}},
    }},{{responsive:true}});
  }}

  // Envelope detail
  Plotly.newPlot('env',[
    {{x:d.rmz,y:d.ri,type:'bar',width:0.003,marker:{{color:'#aab'}},name:'MS1',
      hovertemplate:'m/z:%{{x:.4f}}<br>Int:%{{y:.0f}}<extra></extra>'}},
    {{x:d.m.map(function(p){{return p.mz_t}}),
      y:d.m.map(function(p){{return p.ok?p.i:0}}),
      type:'scatter',mode:'markers+text',
      marker:{{color:d.m.map(function(p){{return p.ok?'#1565c0':'#c62828'}}),size:9,symbol:'diamond'}},
      text:d.m.map(function(p){{return 'M+'+p.iso+(p.ok?' '+p.ppm.toFixed(1)+'ppm':' miss')}}),
      textposition:'top center',textfont:{{size:8}},name:'Theoretical'}},
  ],{{
    title:d.seq+'  z='+d.z+'+  scan:'+d.scan,
    xaxis:{{title:'m/z'}},yaxis:{{title:'Intensity'}},
    margin:{{t:30,b:35,l:50,r:10}},bargap:0,showlegend:true,legend:{{x:0.6,y:1}},
  }},{{responsive:true}});

  // Peak table
  var h='<table><tr><th>Iso</th><th>m/z theo</th><th>m/z obs</th><th>Int</th><th>ppm</th></tr>';
  d.m.forEach(function(p){{
    h+='<tr class="'+(p.ok?'ok':'no')+'"><td>M+'+p.iso+'</td><td>'+p.mz_t.toFixed(4)+
      '</td><td>'+(p.mz_o?p.mz_o.toFixed(4):'-')+'</td><td>'+(p.i>0?p.i.toFixed(0):'-')+
      '</td><td>'+(p.ppm!==null?p.ppm.toFixed(2):'-')+'</td></tr>';
  }});
  h+='</table><div class="info">Avg ppm: '+(d.ppm!==null?d.ppm.toFixed(2):'-')+
    ' | Score: '+d.sc.toFixed(1)+' | Matched: '+d.ni+'/5</div>';
  document.getElementById('pinfo').innerHTML=h;

  // Also select this scan on TIC
  selectScan(d.scan);
}}

// Init
selectScan(G.apex_scan);
if(D.length>0)showPeptide(ord[0]);
</script></body></html>"""


def main():
    parser = argparse.ArgumentParser(description="Peptide evidence viewer")
    parser.add_argument("--mzml", required=True)
    parser.add_argument("--peptides", required=True)
    parser.add_argument("--out", default=None)
    parser.add_argument("--ppm", type=float, default=15.0)
    parser.add_argument("--max-peptides", type=int, default=500)
    args = parser.parse_args()

    peptides = pd.read_csv(args.peptides)
    peptides = peptides.nlargest(args.max_peptides, "score").sort_values(["start", "end"])
    print(f"=== Peptide Evidence Viewer ({len(peptides)} peptides) ===")

    evidence, gd = extract_evidence(Path(args.mzml), peptides, ppm=args.ppm)

    out = Path(args.out) if args.out else Path(args.peptides).parent / "evidence_viewer.html"
    out.write_text(build_html(evidence, gd), encoding="utf-8")
    size_mb = out.stat().st_size / 1024 / 1024
    print(f"  Wrote {out} ({size_mb:.1f} MB)")


if __name__ == "__main__":
    main()
