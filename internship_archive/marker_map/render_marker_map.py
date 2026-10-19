#!/usr/bin/env python3
"""
render_marker_map.py (v4) — Render the Retrotransposon Marker Map. CAN-NOT-CRASH build.

v4 hardening (after 'trouble in init_rank' + 'triangulation failed' crash on full view):
  * EDGE LABELS -> xlabel: research-confirmed fix for the init_rank bug class
    (labeled edges + clusters + rank=same). Applied from rung 1, spine included.
  * FALLBACK LADDER per view, one risk removed per rung, stderr-classified retries:
      r1 newrank=true  + rows + polyline + xlabel
      r2 newrank=false + rows + polyline + xlabel
      r3 rows dropped  + polyline + xlabel
      r4 splines=line  + xlabel
      r5 labels hidden + line
    A rung "succeeds" only if: exit 0, output exists, non-empty, parses as SVG,
    AND stderr contains no fatal markers (init_rank / triangulation / lost edge...).
  * COMPOSITE FALLBACK for 'full': if every rung fails, the full map is assembled
    from the spine + per-section SVGs (which render independently and reliably)
    as nested <svg> panels — XML prologs stripped, IDs prefixed, no reruns needed.
  * JOB ORDER: sections render BEFORE full, so the composite always has material.
  * Failed variants' DOT + stderr are preserved as {name}.failed.{rung}.dot/.log.

Usage:
  python3 render_marker_map.py                       # everything + HTML viewer
  python3 render_marker_map.py --view full           # one view
  python3 render_marker_map.py --view section_A
  python3 render_marker_map.py --dpi 200 --outdir ~/maps
"""
from __future__ import annotations

import argparse, csv, html, json, re, subprocess, sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
import marker_map_data as D  # noqa: E402

# ---------------------------------------------------------------- constants

SPINE_ORDER = ["meth", "rna", "orf2p", "en", "tprt",
               "position", "junction", "lineage", "reseal", "clinical"]

SECTION_ORDER = ["SPINE", "A", "B", "C", "D", "E", "G", "X"]
LEGEND_VIEWS = {"full", "spine"}

SECTION_STYLE = {   # fill, border, title  (Okabe-Ito derived, colorblind-safe)
    "SPINE": ("#ECECEC", "#616161", "CAUSAL SPINE"),
    "A":     ("#EBF8FD", "#0072B2", "A · TARGETING — where it lands"),
    "B":     ("#E6EFEC", "#009E73", "B · CREDIBILITY — is the call real"),
    "C":     ("#FEFCE1", "#B7950B", "C · ELEMENT STATE — is it armed"),
    "D":     ("#FAE5E6", "#D55E00", "D · MEDICAL & DISEASE"),
    "E":     ("#F9EBF6", "#CC79A7", "E · RESOURCES — known sites"),
    "G":     ("#FCF5E6", "#E69F00", "G · OUR DATA"),
    "X":     ("#FDEDEC", "#C0392B", "X · VERIFIED GAPS"),
}
SPINE_NODE_FILL = "#30343B"

EDGE_STYLE = {      # tier -> (style, color, penwidth)
    "demonstrated": ("solid",  "#1A1A1A", 1.6),
    "associated":   ("dashed", "#0072B2", 1.2),
    "proposed":     ("dotted", "#B7950B", 1.2),
    "unknown":      ("dotted", "#8A8A8A", 1.0),
    "gap":          ("dashed", "#C0392B", 1.2),
}
SPINE_EDGE_COLOR, SPINE_EDGE_PW = "#30343B", 2.6

ROW_SIZE = 6
GRID_PREFIXES = {"A": ("a_",), "B": ("b_",), "C": ("c_",), "D": ("d_",),
                 "E": ("e_",), "G": ("g_",), "X": ("x_",)}

ANCHOR_CANDIDATES = {
    "A": ["a_133", "a_motif7", "a_sense"], "B": ["b_triad", "b_callers"],
    "C": ["c_clock", "c_bulk", "c_tubio"], "D": ["d_germline", "d_cancerdrv", "d_tissue"],
    "E": ["e_counts", "e_l1hs", "e_ages"], "G": ["g_hg", "g_lr", "g_plan"],
    "X": ["x_tsdcomp", "x_tss", "x_fda"],
}
ATTACH = {
    "A": ["en", "tprt", "position"],          "B": ["junction", "tprt", "position"],
    "C": ["reseal", "meth", "junction"],      "D": ["clinical", "reseal"],
    "E": ["position", "junction", "lineage"], "G": ["position", "junction", "meth"],
    "X": ["position", "junction", "meth", "reseal"],
}

# ---------------------------------------------------------------- fallback ladder

VARIANTS = [   # rung -> (newrank, rows, splines, labels)
    ("r1", dict(newrank=True,  rows=True,  spl="polyline", labels=True)),
    ("r2", dict(newrank=False, rows=True,  spl="polyline", labels=True)),
    ("r3", dict(newrank=False, rows=False, spl="polyline", labels=True)),
    ("r4", dict(newrank=False, rows=False, spl="line",     labels=True)),
    ("r5", dict(newrank=False, rows=False, spl="line",     labels=False)),
]
FATAL_PATTERNS = [
    "trouble in init_rank", "triangulation failed", "not in any triangle",
    "routesplines", "Pshortestpath", "Error: lost",
]

# ---------------------------------------------------------------- helpers

def esc(s: str) -> str:
    return html.escape(s, quote=True)

def esc_attr(s: str) -> str:
    return s.replace("\\", "\\\\").replace('"', '\\"').replace("\n", " ")

def html_label(n) -> str:
    nid, label, section, kind, evidence, source, views = n
    lines = label.split("\n")
    title, body = lines[0], lines[1:]
    if section == "SPINE":
        tcol, bcol, ccol = "#FFFFFF", "#E8E8E8", "#C9C9C9"
    elif kind == "gap":
        tcol, bcol, ccol = "#C0392B", "#333333", "#777777"
    else:
        tcol, bcol, ccol = "#111111", "#333333", "#777777"
    p = ['<<TABLE BORDER="0" CELLBORDER="0" CELLSPACING="0" CELLPADDING="2">']
    p.append(f'<TR><TD ALIGN="LEFT"><FONT POINT-SIZE="12" COLOR="{tcol}"><B>{esc(title)}</B></FONT></TD></TR>')
    if body:
        p.append(f'<TR><TD ALIGN="LEFT"><FONT POINT-SIZE="9" COLOR="{bcol}">'
                 f'{"<BR ALIGN=\"LEFT\"/>".join(esc(b) for b in body)}</FONT></TD></TR>')
    if source:
        p.append(f'<TR><TD ALIGN="LEFT"><FONT POINT-SIZE="8" COLOR="{ccol}"><I>{esc(source)}</I></FONT></TD></TR>')
    p.append('</TABLE>>')
    return "".join(p)

def node_stmt(n) -> str:
    nid, label, section, kind, evidence, source, views = n
    fill, border, _ = SECTION_STYLE[section]
    if section == "SPINE":
        fill = SPINE_NODE_FILL
    tooltip = f"{source} [{evidence}]" if source else evidence
    return (f'    "{nid}" [label={html_label(n)}, shape={D.KIND_SHAPE[kind]}, '
            f'style="filled", fillcolor="{fill}", color="{border}", '
            f'fontname="Helvetica", tooltip="{esc_attr(tooltip)}"];')

def edge_stmt(e, var, heavy=False) -> str:
    src, dst, tier, label, views = e
    style, color, pw = EDGE_STYLE[tier]
    if heavy:
        color, pw = SPINE_EDGE_COLOR, SPINE_EDGE_PW
    a = [f"style={style}", f'color="{color}"', f"penwidth={pw}",
         "arrowsize=0.7", 'fontname="Helvetica"']
    if heavy:
        a.append("weight=8")
    if label and var["labels"]:
        # xlabel, never label: avoids the init_rank bug class and layout distortion
        a += [f'xlabel="{esc_attr(label)}"', "fontsize=8", 'fontcolor="#555555"']
    return f'  "{src}" -> "{dst}" [{", ".join(a)}];'

# ---------------------------------------------------------------- validation

def validate():
    nodes, edges = list(D.NODES), list(D.EDGES)
    ids, errors, warnings = set(), [], []
    for n in nodes:
        if n[0] in ids: errors.append(f"duplicate node id {n[0]!r}")
        ids.add(n[0])
        if n[2] not in SECTION_STYLE: errors.append(f"{n[0]}: unknown section {n[2]!r}")
        if n[3] not in D.KIND_SHAPE:  errors.append(f"{n[0]}: unknown kind {n[3]!r}")
        if n[4] not in EDGE_STYLE:    errors.append(f"{n[0]}: unknown evidence {n[4]!r}")
        if not set(n[6]) <= {"full", "spine"}: errors.append(f"{n[0]}: bad views {n[6]}")
    for e in edges:
        if e[0] not in ids: errors.append(f"edge src {e[0]!r} unknown")
        if e[1] not in ids: errors.append(f"edge dst {e[1]!r} unknown")
        if e[2] not in EDGE_STYLE: errors.append(f"edge {e[0]}->{e[1]}: bad tier {e[2]!r}")
        if not set(e[4]) <= {"full", "spine"}: errors.append(f"edge {e[0]}->{e[1]}: bad views {e[4]}")
    by_id = {n[0]: n for n in nodes}
    for s, cands in ANCHOR_CANDIDATES.items():
        if not any(by_id[c][2] == s for c in cands if c in by_id):
            warnings.append(f"no preferred anchor for {s}; using first node of section")
    return nodes, edges, errors, warnings

# ---------------------------------------------------------------- views / grid

def view_nodes(nodes, edges, view):
    if view == "full":
        return [n[0] for n in nodes if "full" in n[6]]
    if view == "spine":
        return [n[0] for n in nodes if "spine" in n[6]]
    sec = view.split("_", 1)[1]
    inc = {n[0] for n in nodes if n[2] == sec}
    spine_ids = set(SPINE_ORDER)
    for e in edges:
        if e[0] in inc and e[1] in spine_ids: inc.add(e[1])
        if e[1] in inc and e[0] in spine_ids: inc.add(e[0])
    return [n[0] for n in nodes if n[0] in inc]

def view_edges(edges, view, node_set):
    ns = set(node_set)
    if view == "full":
        return [e for e in edges if "full" in e[4]]
    if view == "spine":
        return [e for e in edges if "spine" in e[4]]
    return [e for e in edges if e[0] in ns and e[1] in ns]

def build_rows(nodes, view_ids, sec):
    pref = GRID_PREFIXES.get(sec)
    if not pref:
        return []
    members = [n[0] for n in nodes if n[2] == sec and n[0] in view_ids
               and n[0].startswith(pref)]
    return [members[i:i + ROW_SIZE] for i in range(0, len(members), ROW_SIZE)] if members else []

def resolve_anchor(nodes, sec, ids):
    by_id = {n[0]: n for n in nodes}
    for c in ANCHOR_CANDIDATES.get(sec, []):
        if c in by_id and by_id[c][2] == sec and c in ids:
            return c
    for n in nodes:
        if n[2] == sec and n[0] in ids:
            return n[0]
    return None

def resolve_attach(nodes, sec, ids):
    for c in ATTACH.get(sec, []):
        if c in ids:
            return c
    for s in SPINE_ORDER:
        if s in ids:
            return s
    return None

# ---------------------------------------------------------------- dot emission

def build_dot(nodes, edges, view, var, with_legend=True):
    ids = view_nodes(nodes, edges, view)
    idset = set(ids)
    eset = view_edges(edges, view, ids)
    by_id = {n[0]: n for n in nodes}
    section_of = {n[0]: n[2] for n in nodes}
    chain_pairs = set(zip(SPINE_ORDER, SPINE_ORDER[1:]))

    gattrs = (f'rankdir=TB, compound=true, splines={var["spl"]}, '
              'nodesep=0.45, ranksep=0.75, pad="0.3", fontname="Helvetica", bgcolor="white"')
    if var["newrank"]:
        gattrs += ", newrank=true"
    if var["labels"]:
        gattrs += ", forcelabels=true"

    L = ['digraph marker_map {', f'  graph [{gattrs}];',
         '  node [fontname="Helvetica"];',
         f'  label="Retrotransposon Marker Map — {view}"; labelloc="t"; fontsize=18;']

    for sec in SECTION_ORDER:
        members = [i for i in ids if section_of[i] == sec]
        if not members:
            continue
        fill, border, title = SECTION_STYLE[sec]
        L.append(f'  subgraph cluster_{sec} {{')
        L.append(f'    label=<<FONT POINT-SIZE="13"><B>{esc(title)}</B></FONT>>; labeljust="l";')
        L.append(f'    style="rounded,filled"; fillcolor="{fill}"; color="{border}"; '
                 f'penwidth=1.4; fontname="Helvetica";')
        for i in members:
            L.append(node_stmt(by_id[i]))
        if sec != "SPINE" and var["rows"]:
            rows = build_rows(nodes, ids, sec)
            for row in rows:
                L.append("    { rank=same; " + " ".join(f'"{i}";' for i in row) + " }")
            for r0, r1 in zip(rows, rows[1:]):
                L.append(f'    "{r0[0]}" -> "{r1[0]}" [style=invis, constraint=true, '
                         f'minlen=1, weight=4];')
        L.append("  }")

    if view != "spine":
        secs = sorted({section_of[i] for i in ids} - {"SPINE"})
        for sec in secs:
            a, t = resolve_anchor(nodes, sec, idset), resolve_attach(nodes, sec, idset)
            if a and t:
                L.append(f'  "{a}" -> "{t}" [style=invis, constraint=true, weight=6, minlen=1];')

    for e in eset:
        if (e[0], e[1]) in chain_pairs:
            L.append(edge_stmt(e, var, heavy=True))
    for e in eset:
        if (e[0], e[1]) not in chain_pairs:
            L.append(edge_stmt(e, var, heavy=False))

    if with_legend:
        L += legend_stmts()
    L.append("}")
    return "\n".join(L) + "\n"

def legend_stmts():
    rows = ['<TR><TD COLSPAN="2"><B>Edge evidence tiers</B></TD></TR>']
    for tier, (st, col, _p) in EDGE_STYLE.items():
        rows.append(f'<TR><TD BGCOLOR="{col}" WIDTH="26"> </TD>'
                    f'<TD ALIGN="LEFT"><B>{tier}</B> — {st}</TD></TR>')
    rows.append('<TR><TD COLSPAN="2"><B>Node kinds</B>: box chain · ellipse marker · '
                'folder resource · hexagon gap · box3d clinical · diamond structure</TD></TR>')
    rows.append('<TR><TD COLSPAN="2"><B>Flags</B>: [V]=verified · [S]=snippet · '
                '[C]=computed · [W]=our data — in each node\'s citation line</TD></TR>')
    table = ('<<TABLE BORDER="1" CELLBORDER="0" CELLSPACING="0" CELLPADDING="4">'
             + "".join(rows) + "</TABLE>>")
    return ['  subgraph cluster_legend {',
            '    label=<<FONT POINT-SIZE="13"><B>LEGEND</B></FONT>>; style="rounded"; color="#999999";',
            f'    "legend" [shape=plaintext, label={table}];', '  }']

# ---------------------------------------------------------------- hardened rendering

def classify_failure(stderr: str) -> str | None:
    for pat in FATAL_PATTERNS:
        if pat in stderr:
            return pat
    return None

def _dot_once(dot_path: Path, out_path: Path, fmt: str, dpi: int) -> tuple[bool, str]:
    cmd = ["dot", f"-T{fmt}"]
    if fmt == "png":
        cmd.append(f"-Gdpi={dpi}")
    cmd += [str(dot_path), "-o", str(out_path)]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    except subprocess.TimeoutExpired:
        return False, "timeout after 600s"
    except FileNotFoundError:
        sys.exit("FATAL: `dot` not found on PATH — install graphviz (brew install graphviz).")
    err = proc.stderr or ""
    if proc.returncode != 0:
        return False, err.strip()[:300] or f"exit {proc.returncode}"
    fatal = classify_failure(err)
    if fatal:
        return False, f"stderr flagged: {fatal}"
    if not out_path.exists() or out_path.stat().st_size == 0:
        return False, "no output written"
    if fmt == "svg":
        head = out_path.read_text(errors="replace").lstrip()[:5]
        if not head.startswith("<"):
            return False, "output is not SVG"
    return True, ""

def render_view(nodes, edges, view, outdir: Path, name: str, formats, dpi: int,
                with_legend: bool) -> str:
    """Run the fallback ladder. Returns a status string. Never raises except on
    environment-level failure (no dot). Final fallback for 'full' = composite."""
    failures = []
    chosen = None
    for fmt in formats:
        got = chosen
        rungs = ([chosen] if chosen else [v for v, _ in VARIANTS])
        for rung, var in [(r, v) for r, v in VARIANTS if not got or r == got]:
            dot = build_dot(nodes, edges, view, var, with_legend=with_legend)
            dp = outdir / f"{name}.dot" if rung == VARIANTS[0][0] else outdir / f"{name}.{rung}.dot"
            dp.write_text(dot)
            ok, why = _dot_once(dp, outdir / f"{name}.{fmt}", fmt, dpi)
            if ok:
                if got is None:
                    got = rung
                    (outdir / f"{name}.dot").write_text(dot)   # canonical = winning rung
                break
            failures.append(f"{fmt}/{rung}: {why.splitlines()[0] if why else 'unknown'}")
            (outdir / f"{name}.failed.{rung}.dot").write_text(dot)
        if got is None:
            continue
        chosen = got
        if fmt == formats[-1]:
            return f"rung {chosen}"
    # ladder exhausted for this view
    if view == "full":
        comp = compose_full_svg(outdir)
        if comp:
            (outdir / "marker_map_full.svg").write_text(comp)
            return "composite fallback"
    sys.exit(f"FATAL: '{name}' failed every fallback rung. Evidence:\n  " +
             "\n  ".join(failures[:10]) +
             f"\nPreserved DOTs: {outdir}/{name}.failed.*.dot")

# ---------------------------------------------------------------- composite fallback

def _read_svg(path: Path) -> str | None:
    if not path.exists():
        return None
    t = path.read_text(errors="replace")
    t = re.sub(r"<\?xml[^>]*\?>", "", t)
    t = re.sub(r"<!DOCTYPE[^>]*>", "", t, flags=re.S)
    return t if "<svg" in t else None

def _prefix_ids(text: str, pfx: str) -> str:
    text = re.sub(r'\bid="([^"]+)"', lambda m: f'id="{pfx}_{m.group(1)}"', text)
    text = re.sub(r"url\(#([^)]+)\)", rf"url(#{pfx}_\1)", text)
    text = re.sub(r'(xlink:href)="#([^"]+)"', rf'\1="#{pfx}_\2"', text)
    text = re.sub(r'(?<![:\w])href="#([^"]+)"', rf'href="#{pfx}_\1"', text)
    return text

def _place(text: str, x, y, w, h) -> str:
    m = re.search(r"<svg\b([^>]*)>", text)
    attrs = m.group(1)
    attrs = re.sub(r'\s(width|height)="[^"]*"', "", attrs)
    new = (f'<svg{attrs} x="{x}" y="{y}" width="{w}" height="{h}" '
           f'preserveAspectRatio="xMidYMid meet">')
    return text[:m.start()] + new + text[m.end():]

def compose_full_svg(outdir: Path) -> str | None:
    CELL_W, CELL_H, M = 5200, 1600, 40
    panels = [("SPINE", "Causal spine", "marker_map_spine.svg")] + \
             [(s, SECTION_STYLE[s][2], f"marker_map_section_{s}.svg")
              for s in SECTION_ORDER if s != "SPINE"]
    body, used = [], []
    for idx, (key, title, fname) in enumerate(panels):
        col, row = idx % 2, idx // 2
        x, y = M + col * (CELL_W + M), 70 + M + row * (CELL_H + M)
        svg = _read_svg(outdir / fname)
        if svg:
            svg = _prefix_ids(svg, key)
            svg = _place(svg, x, y, CELL_W, CELL_H)
            body.append(svg)
            used.append(title)
        else:
            body.append(f'<rect x="{x}" y="{y}" width="{CELL_W}" height="{CELL_H}" '
                        f'fill="#FDEDEC" stroke="#C0392B"/>'
                        f'<text x="{x+40}" y="{y+80}" font-size="40" fill="#C0392B">'
                        f'{esc(title)} — render failed (see failed DOT)</text>')
    W = 2 * CELL_W + 3 * M
    H = 70 + 4 * (CELL_H + M) + M
    head = (f'<svg xmlns="http://www.w3.org/2000/svg" '
            f'xmlns:xlink="http://www.w3.org/1999/xlink" viewBox="0 0 {W} {H}" '
            f'width="{W}pt" height="{H}pt">'
            f'<rect x="0" y="0" width="{W}" height="{H}" fill="white"/>'
            f'<text x="{W/2}" y="50" text-anchor="middle" font-size="34" '
            f'font-family="Helvetica">Retrotransposon Marker Map — composite '
            f'(spine + sections)</text>')
    return head + "".join(body) + "</svg>"

# ---------------------------------------------------------------- HTML viewer

HTML = r"""<!DOCTYPE html>
<html><head><meta charset="utf-8"><title>Retrotransposon Marker Map</title>
<style>
body{margin:0;display:flex;height:100vh;font-family:-apple-system,Helvetica,Arial,sans-serif}
nav{width:210px;padding:14px;background:#1e2430;color:#eee;box-sizing:border-box;overflow:auto}
nav h1{font-size:14px;margin:0 0 12px;color:#9fc3ff}
nav button{display:block;width:100%;text-align:left;margin:4px 0;padding:8px 10px;border:0;
 border-radius:6px;background:#2c3547;color:#dde;cursor:pointer;font-size:13px}
nav button:hover{background:#3b475e} nav p{font-size:11px;color:#aab;line-height:1.5}
main{flex:1;display:flex;flex-direction:column;overflow:hidden}
.viewport{flex:1;position:relative;overflow:hidden;background:#fafafa;touch-action:none;cursor:grab}
.viewport.dragging{cursor:grabbing}
.stage{position:absolute;left:0;top:0;transform-origin:0 0}.stage svg{display:block}
.hud{position:absolute;top:10px;right:10px;display:flex;gap:6px;z-index:5}
.hud button,.viewtitle{border:0;border-radius:6px;background:#ffffffd9;padding:6px 11px;
 font-size:13px;cursor:pointer}.viewtitle{position:absolute;top:10px;left:12px;
 font-weight:600;color:#333;cursor:default;z-index:5}
</style></head><body>
<nav><h1>MARKER MAP</h1><div id="nav"></div>
<p>Wheel = zoom at cursor<br>Drag = pan<br>Double-click = fit<br>Hover nodes = citations</p></nav>
<main><div id="viewer" class="viewport"><div class="stage" id="stage"></div>
<div class="viewtitle" id="vtitle"></div>
<div class="hud"><button onclick="zoomBy(1.3)">+</button><button onclick="zoomBy(0.77)">&minus;</button>
<button onclick="fit()">fit</button></div></div></main>
<script>
const VIEWS=__VIEWS__;
let cur=0,st={x:0,y:0,s:1},stageW=1,stageH=1;
const vp=document.getElementById('viewer'),stage=document.getElementById('stage'),
titleEl=document.getElementById('vtitle'),nav=document.getElementById('nav');
function clampS(s){return Math.min(12,Math.max(0.1,s));}
function apply(){stage.style.transform=`translate(${st.x}px,${st.y}px) scale(${st.s})`;}
function fit(){st.s=Math.min(vp.clientWidth/stageW,vp.clientHeight/stageH,1);
 st.x=(vp.clientWidth-stageW*st.s)/2;st.y=(vp.clientHeight-stageH*st.s)/2;apply();}
function zoomAt(px,py,f){const s2=clampS(st.s*f),k=s2/st.s;
 st.x=px-(px-st.x)*k;st.y=py-(py-st.y)*k;st.s=s2;apply();}
function zoomBy(f){zoomAt(vp.clientWidth/2,vp.clientHeight/2,f);}
function load(i){cur=i;stage.innerHTML=VIEWS[i].svg;
 const svg=stage.querySelector('svg'),vb=svg.viewBox.baseVal;
 stageW=vp.clientWidth;stageH=stageW*vb.height/vb.width;
 svg.style.width=stageW+'px';svg.style.height=stageH+'px';
 nav.querySelectorAll('button').forEach((b,j)=>b.style.background=j===i?'#4b6cb7':'');
 titleEl.textContent=VIEWS[i].title;fit();}
vp.addEventListener('wheel',e=>{e.preventDefault();const r=vp.getBoundingClientRect();
 zoomAt(e.clientX-r.left,e.clientY-r.top,Math.exp(-e.deltaY*0.0015));},{passive:false});
let drag=null;
vp.addEventListener('pointerdown',e=>{drag={x:e.clientX,y:e.clientY,sx:st.x,sy:st.y};
 vp.setPointerCapture(e.pointerId);vp.classList.add('dragging');});
vp.addEventListener('pointermove',e=>{if(!drag)return;
 st.x=drag.sx+(e.clientX-drag.x);st.y=drag.sy+(e.clientY-drag.y);apply();});
vp.addEventListener('pointerup',()=>{drag=null;vp.classList.remove('dragging');});
vp.addEventListener('dblclick',fit);
window.addEventListener('resize',()=>load(cur));
VIEWS.forEach((v,i)=>{const b=document.createElement('button');b.textContent=v.title;
 b.onclick=()=>load(i);nav.appendChild(b);});
load(0);
</script></body></html>"""

def write_html(outdir: Path, rendered: list):
    views = []
    for name, title in rendered:
        p = outdir / f"{name}.svg"
        if p.exists():
            views.append({"key": name, "title": title, "svg": p.read_text()})
    (outdir / "index.html").write_text(HTML.replace("__VIEWS__", json.dumps(views)))

def write_csvs(outdir: Path, nodes, edges):
    with (outdir / "marker_map_nodes.csv").open("w", newline="") as f:
        w = csv.writer(f); w.writerow(["id", "label", "section", "kind", "evidence", "source", "views"])
        for n in nodes:
            w.writerow([n[0], n[1].replace("\n", " | ")] + list(n[2:]))
    with (outdir / "marker_map_edges.csv").open("w", newline="") as f:
        w = csv.writer(f); w.writerow(["src", "dst", "tier", "label", "views"])
        for e in edges:
            w.writerow(list(e[:4]) + [",".join(e[4])])

# ---------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--view", default="all",
                    help="all | spine | full | section_A ... section_X")
    ap.add_argument("--outdir", default=str(SCRIPT_DIR / "marker_map_out"))
    ap.add_argument("--dpi", type=int, default=150)
    ap.add_argument("--no-html", dest="no_html", action="store_true")
    ap.add_argument("--legend", action="store_true", help="force legend in section views")
    a = ap.parse_args()

    nodes, edges, errors, warnings = validate()
    print("Validating data...")
    for e in errors:   print(f"  ERROR: {e}")
    for w in warnings: print(f"  warn:  {w}")
    if errors:
        sys.exit(f"validation failed with {len(errors)} error(s)")
    for s in SECTION_ORDER:
        cnt = sum(1 for n in nodes if n[2] == s)
        if cnt: print(f"  {SECTION_STYLE[s][2]:<42} {cnt}")
    print(f"OK: {len(nodes)} nodes, {len(edges)} edges validated")

    outdir = Path(a.outdir).expanduser(); outdir.mkdir(parents=True, exist_ok=True)

    # Order: spine -> sections -> full LAST (composite fallback needs section SVGs)
    if a.view == "all":
        jobs = [("spine", "svg,png,pdf")] + \
               [(s, "svg,png") for s in SECTION_ORDER if s != "SPINE"] + \
               [("full", "svg,png,pdf")]
    else:
        key = a.view.replace("section_", "")
        fmts = "svg,png,pdf" if key in ("spine", "full") else "svg,png"
        jobs = [(key, fmts)]

    results = {}   # name -> (title, status)
    for sec, fmts in jobs:
        if sec == "spine":
            view, name, title = "spine", "marker_map_spine", "Causal spine"
        elif sec == "full":
            view, name, title = "full", "marker_map_full", "Full map"
        else:
            view, name, title = f"section_{sec}", f"marker_map_section_{sec}", SECTION_STYLE[sec][2]
        legend = view in LEGEND_VIEWS or a.legend
        status = render_view(nodes, edges, view, outdir, name,
                             tuple(fmts.split(",")), a.dpi, with_legend=legend)
        results[name] = (title, status)
        print(f"  rendered {name}  ({fmts}{', legend' if legend else ''}, {status})")

    write_csvs(outdir, nodes, edges)
    print("  marker_map_nodes.csv / marker_map_edges.csv")
    if not a.no_html:
        order = ["marker_map_spine"] + \
                ([k for k in ("marker_map_full",) if k in results] if
                 results.get("marker_map_full", ("", ""))[1] != "composite fallback" else []) + \
                [f"marker_map_section_{s}" for s in SECTION_ORDER if s != "SPINE"]
        # if full fell back to composite, show it too (it exists as SVG)
        disp = [(n, results[n][0]) for n in order if n in results]
        if "marker_map_full" in results and \
           results["marker_map_full"][1] == "composite fallback":
            disp.insert(1, ("marker_map_full", "Full map (composite)"))
        write_html(outdir, disp)
        print("  index.html  (offline pan/zoom viewer)")
    print(f"\nOpen: open {outdir}/index.html")

if __name__ == "__main__":
    main()
