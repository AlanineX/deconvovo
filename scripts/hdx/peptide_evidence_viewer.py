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
                     n_explore_scans: int = 50,
                     scan_col: str = "scan", score_col: str = "score") -> tuple[list[dict], dict]:
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
    peptide_scans = set(int(s) for s in peptides[scan_col].unique() if s <= len(scan_meta))
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
        scan_id = int(row[scan_col])
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

        # XIC: extract monoisotopic intensity across all scans
        xic_rt, xic_int = [], []
        tol_xic = mz_mono * ppm / 1e6
        for sm2 in scan_meta:
            xic_rt.append(round(sm2["rt"], 3))
            mz2, int2 = sm2["mz"], sm2["intensity"]
            if len(mz2) > 0:
                xmask = np.abs(mz2 - mz_mono) <= tol_xic
                xic_int.append(round(float(int2[xmask].max()), 1) if xmask.any() else 0)
            else:
                xic_int.append(0)

        # Use score_col-safe access
        sc_val = float(row.get("score_final", row.get("score", 0)))

        evidence.append({
            "seq": row["sequence"],
            "s": int(row["start"]), "e": int(row["end"]),
            "z": z, "mz": round(mz_mono, 4),
            "sc": round(sc_val, 1),
            "ni": n_ok, "scan": scan_id,
            "ppm": round(float(avg_ppm), 2) if not np.isnan(avg_ppm) else None,
            "rmz": [round(float(m), 4) for m in mz_arr[mask]],
            "ri": [round(float(i), 1) for i in int_arr[mask]],
            "m": matched,
            "xrt": xic_rt, "xi": xic_int,
        })

    print(f"    {len(evidence)} peptides with evidence")
    return evidence, global_data


def _stick_trace(mz, intensity):
    """Convert arrays to zero-interleaved stick spectrum (vertical lines)."""
    x, y = [], []
    for m, i in zip(mz, intensity):
        x.extend([m, m, None])
        y.extend([0, i, None])
    return x, y


def build_html(evidence, gd, title="Peptide Evidence Viewer"):
    ev_j = json.dumps(evidence)
    gd_j = json.dumps(gd)
    N = NEUTRON
    return f"""<!DOCTYPE html>
<html><head><meta charset="utf-8"><title>{title}</title>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<style>
*{{box-sizing:border-box}}
body{{font-family:'Liberation Sans',Arial,sans-serif;margin:0;padding:8px;background:#fff}}
#wrap{{display:flex;gap:8px;height:calc(100vh - 56px)}}
#plots{{flex:1;display:flex;flex-direction:column;gap:4px;min-width:0}}
#plots .p{{background:#fff;border:1px solid #e0e0e0;border-radius:4px;flex:1;min-height:0}}
#right{{width:440px;flex-shrink:0;display:flex;flex-direction:column;gap:6px}}
#tbl-scroll{{flex:1;overflow-y:auto;border:1px solid #e0e0e0;border-radius:4px;background:#fff}}
table{{width:100%;border-collapse:collapse;font-size:11px;cursor:pointer}}
th{{background:#f4f4f8;padding:5px 6px;text-align:left;position:sticky;top:0;
    cursor:pointer;user-select:none;white-space:nowrap;border-bottom:2px solid #ddd}}
th:hover{{background:#e0e0ee}} th.asc::after{{content:' \\25B2'}} th.desc::after{{content:' \\25BC'}}
td{{padding:3px 6px;border-bottom:1px solid #f0f0f0;white-space:nowrap;font-size:10px}}
tr:hover{{background:#e8f0ff}} tr.sel{{background:#c0d8ff}} tr.mz-hl{{background:#fff8e0}}
#pinfo{{font-size:10px;border:1px solid #e0e0e0;border-radius:4px;padding:6px;background:#fafafe}}
#pinfo td{{padding:1px 5px}}
.ok{{color:#2e7d32;font-weight:600}} .no{{color:#c62828}}
#status{{font-size:11px;color:#555;padding:4px 8px;background:#f8f8fc;border-radius:4px;
         border:1px solid #e0e0e0}}
h2{{font-size:16px;margin:0 0 6px;color:#333;font-weight:600}}
.hint{{font-size:10px;color:#999;padding:2px 0}}
</style></head><body>
<h2>{title}</h2>
<div id="status"></div>

<div id="wrap">
<div id="plots">
  <div id="tic" class="p"></div>
  <div id="xic" class="p"></div>
  <div id="ms1" class="p"></div>
  <div id="env" class="p"></div>
</div>
<div id="right">
  <div class="hint">Click headers to sort. Click row to view. TIC: click to change scan.</div>
  <div id="tbl-scroll">
    <table><thead><tr>
      <th data-c="seq">Sequence</th><th data-c="s">Start</th><th data-c="e">End</th>
      <th data-c="z">z</th><th data-c="mz">m/z</th><th data-c="sc">Score</th>
      <th data-c="ni">Iso</th><th data-c="ppm">ppm</th>
    </tr></thead><tbody id="tb"></tbody></table>
  </div>
  <div id="pinfo"></div>
</div>
</div>

<script>
var D={ev_j};
var G={gd_j};
var N={N};
var curScan=G.apex_scan, curPep=-1, sortC=null, sortA=true;
var ord=D.map(function(_,i){{return i}});

// ---- Font config (IMMS pattern) ----
var FS=13, FF="Liberation Sans, Arial, sans-serif", LC="#333";
var tf={{family:FF,size:Math.round(FS*1.1),color:LC}};
var tkf={{family:FF,size:FS,color:'#666'}};
var gridStyle={{gridwidth:0.5,gridcolor:'#e8e8e8',griddash:'dot'}};
var ax={{mirror:true,linewidth:1,linecolor:'#ccc',zeroline:false}};

// ---- Stick spectrum helper ----
function stickData(mz,inten){{
  var x=[],y=[];
  for(var i=0;i<mz.length;i++){{x.push(mz[i],mz[i],null);y.push(0,inten[i],null);}}
  return {{x:x,y:y}};
}}

// ---- Binary search (from IMMS viewer) ----
function bsearch(arr,val){{
  var lo=0,hi=arr.length-1;
  while(lo<=hi){{var mid=(lo+hi)>>1;if(arr[mid]<val)lo=mid+1;else hi=mid-1;}}
  return lo;
}}

// ---- 1. TIC ----
var rtMin=G.tic_rt[0],rtMax=G.tic_rt[G.tic_rt.length-1];
Plotly.newPlot('tic',[
  {{x:G.tic_rt,y:G.tic_int,type:'scatter',mode:'lines',
    line:{{color:'#333',width:0.8}},hovertemplate:'RT: %{{x:.2f}} min<br>TIC: %{{y:.0f}}<extra></extra>'}},
  {{x:[0,0],y:[0,1],type:'scatter',mode:'lines',
    line:{{color:'#e53935',width:1.5,dash:'dot'}},showlegend:false}},
],{{
  margin:{{t:24,b:24,l:55,r:8}},
  title:{{text:'TIC Chromatogram',font:tf}},
  xaxis:Object.assign({{title:'RT (min)',titlefont:tkf,tickfont:tkf,range:[rtMin,rtMax]}},ax,gridStyle),
  yaxis:Object.assign({{title:'TIC',titlefont:tkf,tickfont:tkf,fixedrange:true}},ax,gridStyle),
  showlegend:false,
}},{{responsive:true,displayModeBar:false}});

document.getElementById('tic').on('plotly_click',function(ed){{
  if(!ed||!ed.points||!ed.points[0])return;
  var sn=G.tic_scans[ed.points[0].pointIndex];
  var best=G.stored_scans[0],bd=Math.abs(sn-best);
  G.stored_scans.forEach(function(s){{var d=Math.abs(sn-s);if(d<bd){{best=s;bd=d}}}});
  selectScan(best);
}});

// ---- 2. XIC (same RT scale as TIC) ----
Plotly.newPlot('xic',[
  {{x:[],y:[],type:'scatter',mode:'lines',line:{{color:'#1565c0',width:1}},
    fill:'tozeroy',fillcolor:'rgba(21,101,192,0.08)',
    hovertemplate:'RT: %{{x:.2f}} min<br>Int: %{{y:.0f}}<extra></extra>'}},
],{{
  margin:{{t:24,b:24,l:55,r:8}},
  title:{{text:'XIC (select a peptide)',font:tf}},
  xaxis:Object.assign({{title:'RT (min)',titlefont:tkf,tickfont:tkf,range:[rtMin,rtMax]}},ax,gridStyle),
  yaxis:Object.assign({{title:'Intensity',titlefont:tkf,tickfont:tkf,fixedrange:true}},ax,gridStyle),
  showlegend:false,
}},{{responsive:true,displayModeBar:false}});

// ---- 3. MS1 scan (stick spectrum, y fixed) ----
Plotly.newPlot('ms1',[
  {{x:[],y:[],type:'scatter',mode:'lines',line:{{color:'#333',width:0.6}},
    hovertemplate:'m/z: %{{x:.4f}}<br>Int: %{{y:.0f}}<extra></extra>',name:'MS1'}},
  {{x:[],y:[],type:'scatter',mode:'markers',marker:{{color:'#e53935',size:7,symbol:'triangle-down'}},
    name:'Peptide',hovertemplate:'%{{text}}<extra></extra>',text:[]}},
],{{
  margin:{{t:24,b:24,l:55,r:8}},
  title:{{text:'MS1 Scan',font:tf}},
  xaxis:Object.assign({{title:'m/z',titlefont:tkf,tickfont:tkf}},ax,gridStyle),
  yaxis:Object.assign({{title:'Intensity',titlefont:tkf,tickfont:tkf,fixedrange:true}},ax,gridStyle),
  showlegend:true,legend:{{x:0.75,y:0.95,font:{{size:10}}}},
}},{{responsive:true}});

// rAF debounce on MS1 zoom -> highlight table
var _raf=0;
document.getElementById('ms1').on('plotly_relayout',function(ed){{
  cancelAnimationFrame(_raf);
  _raf=requestAnimationFrame(highlightMzRange);
}});
function highlightMzRange(){{
  var la=document.getElementById('ms1')._fullLayout;
  if(!la||!la.xaxis)return;
  var r=la.xaxis.range;
  if(!r)return;
  document.querySelectorAll('tr.mz-hl').forEach(function(t){{t.classList.remove('mz-hl')}});
  D.forEach(function(d,i){{
    if(d.mz>=r[0]&&d.mz<=r[1]){{
      var t=document.querySelector('tr[data-idx="'+i+'"]');
      if(t&&!t.classList.contains('sel'))t.classList.add('mz-hl');
    }}
  }});
}}

// ---- 4. Envelope (stick spectrum) ----
Plotly.newPlot('env',[
  {{x:[],y:[],type:'scatter',mode:'lines',line:{{color:'#888',width:0.8}},name:'MS1',
    hovertemplate:'m/z: %{{x:.4f}}<br>Int: %{{y:.0f}}<extra></extra>'}},
  {{x:[],y:[],type:'scatter',mode:'markers+text',
    marker:{{size:8,symbol:'diamond'}},textposition:'top center',
    textfont:{{size:9,family:FF}},name:'Theoretical'}},
],{{
  margin:{{t:24,b:28,l:55,r:8}},
  title:{{text:'Isotopic Envelope',font:tf}},
  xaxis:Object.assign({{title:'m/z',titlefont:tkf,tickfont:tkf}},ax,gridStyle),
  yaxis:Object.assign({{title:'Intensity',titlefont:tkf,tickfont:tkf,fixedrange:true}},ax,gridStyle),
  showlegend:true,legend:{{x:0.65,y:0.95,font:{{size:10}}}},bargap:0,
}},{{responsive:true,displayModeBar:false}});

// ==== selectScan ====
function selectScan(sn){{
  curScan=sn;
  var idx=G.tic_scans.indexOf(sn); if(idx<0)idx=0;
  var rt=G.tic_rt[idx];
  var ymax=Math.max.apply(null,G.tic_int);
  Plotly.restyle('tic',{{x:[[rt,rt]],y:[[0,ymax]]}},1);

  var sp=G.scan_spectra[sn];
  if(sp){{
    var st=stickData(sp.mz,sp['int']);
    Plotly.restyle('ms1',{{x:[st.x],y:[st.y]}},0);
    Plotly.relayout('ms1',{{'title.text':'MS1 Scan '+sn+' (RT: '+rt.toFixed(2)+' min)'}});
  }}

  var nAt=D.filter(function(d){{return d.scan===sn}}).length;
  document.getElementById('status').textContent=
    D.length+' peptides | Scan '+sn+' | RT '+rt.toFixed(2)+' min | '+nAt+' at this scan';
}}

// ==== Table ====
var tb=document.getElementById('tb');
function renderTable(){{
  tb.innerHTML='';
  ord.forEach(function(i){{
    var d=D[i],tr=document.createElement('tr');
    tr.dataset.idx=i; tr.dataset.scan=d.scan;
    tr.innerHTML='<td>'+d.seq+'</td><td>'+d.s+'</td><td>'+d.e+'</td><td>'+d.z+
      '</td><td>'+d.mz.toFixed(2)+'</td><td>'+d.sc.toFixed(0)+
      '</td><td>'+d.ni+'</td><td>'+(d.ppm!==null?d.ppm.toFixed(1):'-')+'</td>';
    tr.onclick=function(){{showPeptide(i)}};
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

// ==== showPeptide ====
function showPeptide(idx){{
  curPep=idx;
  document.querySelectorAll('tr.sel').forEach(function(r){{r.classList.remove('sel')}});
  var tr=document.querySelector('tr[data-idx="'+idx+'"]');
  if(tr){{tr.classList.add('sel');tr.scrollIntoView({{block:'nearest'}});}}
  var d=D[idx];

  // XIC (same RT axis as TIC)
  if(d.xrt&&d.xi){{
    Plotly.restyle('xic',{{x:[d.xrt],y:[d.xi]}},0);
    Plotly.relayout('xic',{{'title.text':'XIC: '+d.seq+' (m/z '+d.mz.toFixed(2)+', z='+d.z+'+)'}});
  }}

  // MS1 scan — stick spectrum + peptide markers, auto-zoom to peptide m/z
  var sp=G.scan_spectra[d.scan];
  if(sp){{
    var st=stickData(sp.mz,sp['int']);
    var mk=d.m.map(function(p){{return p.mz_t}});
    var mi=d.m.map(function(p){{return p.ok?p.i:0}});
    var mt=d.m.map(function(p){{return p.ok?'M+'+p.iso:''}});
    Plotly.restyle('ms1',{{x:[st.x,mk],y:[st.y,mi],text:[null,mt]}},[0,1]);
    Plotly.relayout('ms1',{{
      'title.text':'MS1 Scan '+d.scan+' (RT: '+G.tic_rt[G.tic_scans.indexOf(d.scan)].toFixed(2)+' min)',
      'xaxis.range':[d.mz-2, d.mz+7*N/d.z+2],
    }});
  }}

  // Envelope — stick spectrum + theoretical diamonds
  var est=stickData(d.rmz,d.ri);
  var emc=d.m.map(function(p){{return p.ok?'#1565c0':'#c62828'}});
  var emt=d.m.map(function(p){{return 'M+'+p.iso+(p.ok?' '+p.ppm.toFixed(1)+'ppm':' miss')}});
  Plotly.restyle('env',{{
    x:[est.x, d.m.map(function(p){{return p.mz_t}})],
    y:[est.y, d.m.map(function(p){{return p.ok?p.i:0}})],
    'marker.color':[null, emc],
    text:[null, emt],
  }},[0,1]);
  Plotly.relayout('env',{{'title.text':d.seq+'  z='+d.z+'+  scan:'+d.scan}});

  // Peak info table
  var h='<table style="width:100%"><tr><th>Iso</th><th>m/z theo</th><th>m/z obs</th><th>Int</th><th>ppm</th></tr>';
  d.m.forEach(function(p){{
    h+='<tr class="'+(p.ok?'ok':'no')+'"><td>M+'+p.iso+'</td><td>'+p.mz_t.toFixed(4)+
      '</td><td>'+(p.mz_o?p.mz_o.toFixed(4):'-')+'</td><td>'+(p.i>0?p.i.toFixed(0):'-')+
      '</td><td>'+(p.ppm!==null?p.ppm.toFixed(2):'-')+'</td></tr>';
  }});
  h+='</table><div style="margin-top:4px;color:#666;font-size:10px">Avg ppm: '+(d.ppm!==null?d.ppm.toFixed(2):'-')+
    ' | Score: '+d.sc.toFixed(1)+' | Matched: '+d.ni+'/5</div>';
  document.getElementById('pinfo').innerHTML=h;

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
    score_col = "score_final" if "score_final" in peptides.columns else "score"
    scan_col = "scan_apex" if "scan_apex" in peptides.columns else "scan"
    peptides = peptides.nlargest(args.max_peptides, score_col).sort_values(["start", "end"])
    print(f"=== Peptide Evidence Viewer ({len(peptides)} peptides) ===")

    evidence, gd = extract_evidence(Path(args.mzml), peptides, ppm=args.ppm,
                                     scan_col=scan_col, score_col=score_col)

    out = Path(args.out) if args.out else Path(args.peptides).parent / "evidence_viewer.html"
    out.write_text(build_html(evidence, gd), encoding="utf-8")
    size_mb = out.stat().st_size / 1024 / 1024
    print(f"  Wrote {out} ({size_mb:.1f} MB)")


if __name__ == "__main__":
    main()
