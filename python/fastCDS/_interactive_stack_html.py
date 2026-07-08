"""Multi-isoform variant of the TFRegDB2 viewer.

Ported from `tfregdb2/src/components/tf/IsoformViewer.tsx` + `GenomicStructureView.tsx`.

Key idea: all isoforms render on a **single shared axis** built from the
union of every isoform's features. Compact mode collapses the gaps
between union-exon-blocks to 80 virtual bp; gaps where a given isoform
is missing an exon appear as empty space aligned with the same
coordinate column on every other row, so you can compare exon usage at
a glance.

The single-isoform viewer in `_tfregdb_html.py` is untouched — this
module emits a separate self-contained HTML page for the stack case.
"""

from __future__ import annotations

import html
import json

from ._interactive_html import (
    COMPACT_INTRON_BP, DOMAIN_COLOR, _segments_to_payload,
)


def _payloads(segs_by_id: dict[str, list]) -> list[dict]:
    """Convert {input_id: [Segment, ...]} → list of payload dicts.

    Empty isoforms are dropped silently (no features = nothing to plot).
    """
    out = []
    for k, segs in segs_by_id.items():
        if not segs:
            continue
        p = _segments_to_payload(segs)
        if not p:
            continue
        p["inputId"] = k
        out.append(p)
    return out


HTML_TEMPLATE = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>{title_safe}</title>
<style>
  body {{
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
    background: #f8fafc; color: #0f172a; margin: 0; padding: 24px;
  }}
  .panel {{
    background: white; border: 1px solid #e2e8f0; border-radius: 8px;
    padding: 24px; max-width: 1280px; margin: 0 auto;
    box-shadow: 0 1px 2px rgba(15, 23, 42, 0.04);
  }}
  h1 {{ font-size: 16px; margin: 0 0 6px 0; }}
  .meta {{ font-size: 12px; color: #475569; margin-bottom: 14px; }}
  .meta b {{ color: #0f172a; }}
  .controls {{ display: flex; flex-wrap: wrap; gap: 12px; align-items: center;
              margin-bottom: 14px; font-size: 12px; }}
  .controls label {{ display: inline-flex; align-items: center; cursor: pointer; user-select: none; }}
  .controls input[type=radio], .controls input[type=checkbox] {{ margin-right: 4px; }}
  .controls button {{
    font-size: 11px; padding: 3px 10px; border: 1px solid #cbd5e1;
    background: #f8fafc; border-radius: 4px; cursor: pointer;
  }}
  .controls button:hover {{ background: #e2e8f0; }}
  .legend {{ display: flex; flex-wrap: wrap; gap: 14px; font-size: 11px;
            color: #475569; margin-bottom: 12px; }}
  .legend .swatch {{ display: inline-block; vertical-align: middle; margin-right: 4px;
                    border-radius: 2px; border: 1px solid rgba(0,0,0,0.1); }}
  .legend .sw-utr     {{ width: 14px; height: 8px;  background: #f0c078; }}
  .legend .sw-cds     {{ width: 14px; height: 16px; background: #1f77b4; }}
  .legend .sw-cdsd    {{ width: 14px; height: 16px; background: #d62728; }}
  .legend .sw-line    {{ width: 16px; height: 2px; background: #64748b; vertical-align: middle;
                       display: inline-block; margin-right: 4px; }}
  .legend .sw-redline {{ width: 16px; height: 2px; background: #d62728; vertical-align: middle;
                       display: inline-block; margin-right: 4px; }}
  .strand-arrow {{ width: 100%; height: 22px; position: relative; margin-bottom: 6px;
                  color: #475569; }}
  .strand-arrow svg {{ width: 100%; height: 100%; display: block; }}
  .strand-arrow .end-label {{
    position: absolute; top: 50%; transform: translateY(-50%);
    font-size: 11px; font-weight: 600; color: #475569;
    background: #f1f5f9; padding: 2px 6px; border-radius: 3px;
    border: 1px solid #e2e8f0; line-height: 1;
    font-family: ui-monospace, 'SF Mono', Menlo, monospace;
  }}
  .strand-arrow .end-label.left  {{ left: 0; }}
  .strand-arrow .end-label.right {{ right: 0; }}
  /* Stack of isoform rows — labels on the left, tracks on the right
     all sharing a single axis. */
  .iso-stack {{ display: flex; flex-direction: column; gap: 8px; }}
  .iso-row {{ display: grid; grid-template-columns: 180px 1fr; gap: 12px; align-items: center; }}
  .iso-label {{
    font-size: 11px; color: #475569; text-align: right;
    font-family: ui-monospace, 'SF Mono', Menlo, monospace;
    line-height: 1.2; overflow: hidden;
  }}
  .iso-label .tx-id  {{ color: #0f172a; font-weight: 600; display: block; }}
  .iso-label .pr-id  {{ display: block; color: #64748b; }}
  .iso-label .badge {{
    display: inline-block; padding: 1px 4px; border-radius: 999px;
    background: #eff6ff; color: #2563eb; font-size: 9px; margin-top: 2px;
    font-family: -apple-system, BlinkMacSystemFont, sans-serif;
  }}
  .iso-label .badge.mane {{ background: #ecfdf5; color: #059669; }}
  .iso-row.canonical .iso-label .tx-id {{ color: #2563eb; }}
  .iso-plot {{
    position: relative; width: 100%; height: {plot_height}px;
    background: linear-gradient(180deg, #fbfdff 0%, #f3f7fb 100%);
    border: 1px solid #e2e8f0; border-radius: 4px;
    overflow: hidden; user-select: none; cursor: crosshair;
  }}
  .iso-plot.dragging {{ cursor: grabbing; }}
  .iso-plot .baseline {{
    position: absolute; left: 0; right: 0; top: 50%; height: 1px;
    background: #cbd5e1;
  }}
  .feat {{ position: absolute; border-radius: 2px; box-shadow: 0 1px 1px rgba(0,0,0,0.06); cursor: help; }}
  .feat.cds        {{ background: #1f77b4; height: 45%; top: 27.5%; }}
  .feat.cds-domain {{ background: #d62728; height: 45%; top: 27.5%; }}
  .feat.utr        {{ background: #f0c078; height: 20%; top: 40%; }}
  .feat.intron {{
    position: absolute; height: 22.5%; top: 38.75%;
    display: flex; align-items: center; justify-content: space-around;
    pointer-events: none;
  }}
  .feat.intron::before {{
    content: ''; position: absolute; left: 0; right: 0; top: 50%; height: 1px;
    background: #94a3b8;
  }}
  .feat.intron.in-domain::before {{ background: #d62728; height: 2px; }}
  .chevron {{ width: 12px; height: 12px; color: #64748b; position: relative; z-index: 1; }}
  .feat.intron.in-domain .chevron {{ color: #d62728; }}
  .box-shade {{ position: absolute; top: 0; bottom: 0;
               background: rgba(15, 23, 42, 0.40);
               pointer-events: none; z-index: 50; }}
  .box-rect {{ position: absolute; top: 0; bottom: 0;
              border-left: 2px solid #4f46e5;
              border-right: 2px solid #4f46e5;
              background: rgba(79, 70, 229, 0.10);
              box-shadow: inset 0 0 0 1px rgba(255,255,255,0.4);
              pointer-events: none; z-index: 51; }}
  .box-readout {{
    position: absolute; top: 4px;
    background: #4f46e5; color: white;
    font-size: 10.5px; font-weight: 600;
    padding: 2px 6px; border-radius: 3px;
    pointer-events: none; z-index: 52; white-space: nowrap;
    box-shadow: 0 1px 4px rgba(15, 23, 42, 0.25);
  }}
  /* Ticks row aligned to the track column (180px left gutter). */
  .ticks-wrap {{ display: grid; grid-template-columns: 180px 1fr; gap: 12px; margin-top: 6px; }}
  .ticks {{ position: relative; height: 18px; width: 100%;
           font-size: 10px; color: #64748b; }}
  .tick {{
    position: absolute; top: 0; transform: translateX(-50%);
    border-left: 1px solid #cbd5e1; padding-left: 2px; padding-top: 4px;
    white-space: nowrap;
  }}
  .tick.major {{ color: #0f172a; }}
  /* Minimap — shared, sits in the track column. */
  .minimap-wrap {{ display: grid; grid-template-columns: 180px 1fr; gap: 12px;
                  margin-top: 14px; padding-top: 8px; border-top: 1px solid #e2e8f0; }}
  .minimap-label {{ font-size: 11px; color: #64748b; margin-bottom: 4px;
                   display: flex; justify-content: space-between; }}
  .minimap {{ position: relative; width: 100%; height: 36px;
             background: #f8fafc; border: 1px solid #e2e8f0; border-radius: 3px;
             overflow: hidden; cursor: pointer; }}
  .mini-iso-row {{ position: relative; }}
  .mini-feat {{ position: absolute; border-radius: 1px; }}
  .mini-feat.cds        {{ background: #1f77b4; }}
  .mini-feat.cds-domain {{ background: #d62728; }}
  .viewport {{
    position: absolute; top: 0; height: 100%;
    box-sizing: border-box;   /* keep the 2px side borders inside width:100% so
                                 the right edge isn't pushed past the clipped
                                 minimap and rendered thinner than the left */
    background: rgba(79, 70, 229, 0.16);
    border-left: 2px solid rgba(79, 70, 229, 0.7);
    border-right: 2px solid rgba(79, 70, 229, 0.7);
    cursor: grab;
  }}
  .viewport.dragging {{ cursor: grabbing; }}
  .viewport-handle {{
    position: absolute; top: 0; height: 100%; width: 6px;
    background: rgba(79, 70, 229, 0.7); cursor: ew-resize;
  }}
  /* flush to the inner edges so neither handle is clipped at the extremes */
  .viewport-handle.left  {{ left: 0; }}
  .viewport-handle.right {{ right: 0; }}
  .help-text {{ font-size: 10.5px; color: #94a3b8; margin-top: 6px; }}
  .tooltip {{
    position: fixed; background: white; border: 1px solid #cbd5e1; border-radius: 6px;
    padding: 8px 10px; font-size: 12px; box-shadow: 0 4px 12px rgba(15,23,42,0.15);
    pointer-events: none; z-index: 100; max-width: 260px;
  }}
  .tooltip-title {{ font-weight: 600; margin-bottom: 2px; color: #0f172a; }}
  .tooltip-subtitle {{ color: #475569; }}
  h1 a {{ color: #2563eb; text-decoration: underline; font-size: 12px;
         margin-left: 8px; font-weight: normal; }}
</style>
</head>
<body>
<div class="panel">
  <h1 id="header"></h1>
  <div class="meta" id="meta"></div>
  <div class="legend">
    <span><span class="swatch sw-utr"></span>UTR</span>
    <span><span class="swatch sw-cds"></span>CDS</span>
    <span><span class="swatch sw-cdsd"></span>CDS (in domain)</span>
    <span><span class="sw-line"></span>intron</span>
    <span><span class="sw-redline"></span>intron (inside domain span)</span>
  </div>
  <div class="controls">
    <strong>Layout:</strong>
    <label><input type="radio" name="mode" value="compact" checked> Compact (introns = {compact_intron_bp} bp)</label>
    <label><input type="radio" name="mode" value="genomic"> True genomic scale</label>
    <label style="margin-left: 12px;"><input type="checkbox" id="show-utr" checked> Show UTRs</label>
    <button id="btn-reset">Reset zoom</button>
    <button id="btn-fit-domain">Fit to domain</button>
  </div>
  <div style="display: grid; grid-template-columns: 180px 1fr; gap: 12px;">
    <div></div>
    <div class="strand-arrow" id="strand-arrow"></div>
  </div>
  <div class="iso-stack" id="iso-stack"></div>
  <div class="ticks-wrap">
    <div></div>
    <div class="ticks" id="ticks"></div>
  </div>
  <div class="minimap-wrap">
    <div></div>
    <div>
      <div class="minimap-label">
        <span>Overview (union of all isoforms)</span><span id="viewport-readout"></span>
      </div>
      <div class="minimap" id="minimap"></div>
      <div class="help-text">Drag on any isoform row to box-zoom · scroll to zoom around cursor · double-click any row to reset · drag the indigo minimap rectangle to pan/resize · click empty minimap to recenter</div>
    </div>
  </div>
</div>
<div class="tooltip" id="tooltip" style="display:none"></div>
<script>
const STRUCTURES = {payloads_json};
const COMPACT_INTRON_BP = {compact_intron_bp};

/* ============================================================
 * Build a shared axis from the union of every isoform's features.
 * ============================================================ */
function buildSharedAxis(structs, mode, includeUtr) {{
  let lo = Infinity, hi = -Infinity;
  const ranges = [];
  for (const s of structs) {{
    for (const f of s.features) {{
      if (f.type === 'intron') continue;
      if (!includeUtr && (f.type === '5utr' || f.type === '3utr')) continue;
      if (f.start < lo) lo = f.start;
      if (f.end   > hi) hi = f.end;
      ranges.push({{ start: f.start, end: f.end }});
    }}
  }}
  if (!isFinite(lo) || !isFinite(hi)) {{
    return {{ mode, lo: 0, hi: 1, totalVirtual: 1, posToV: () => 0, ticks: [] }};
  }}
  const span = Math.max(1, hi - lo + 1);
  if (mode === 'genomic') {{
    const posToV = (pos) => Math.max(0, Math.min(span, pos - lo));
    return {{
      mode, lo, hi, totalVirtual: span, posToV,
      ticks: [{{ pos: lo, v: 0, major: true }}, {{ pos: hi, v: span, major: true }}],
    }};
  }}
  ranges.sort((a, b) => a.start - b.start || a.end - b.end);
  const merged = [];
  for (const r of ranges) {{
    const last = merged[merged.length - 1];
    if (!last || r.start > last.end + 1) merged.push({{ ...r }});
    else last.end = Math.max(last.end, r.end);
  }}
  const segs = [];
  let cursor = 0;
  for (let i = 0; i < merged.length; i++) {{
    if (i > 0) {{
      const gs = merged[i-1].end + 1, ge = merged[i].start - 1;
      segs.push({{ kind: 'gap', gStart: gs, gEnd: ge,
                   vStart: cursor, vEnd: cursor + COMPACT_INTRON_BP }});
      cursor += COMPACT_INTRON_BP;
    }}
    const r = merged[i];
    const realLen = Math.max(1, r.end - r.start + 1);
    segs.push({{ kind: 'exon', gStart: r.start, gEnd: r.end,
                 vStart: cursor, vEnd: cursor + realLen }});
    cursor += realLen;
  }}
  const totalVirtual = Math.max(1, cursor);
  const posToV = (pos) => {{
    if (pos <= lo) return 0;
    if (pos >= hi) return totalVirtual;
    for (const sg of segs) {{
      if (pos >= sg.gStart && pos <= sg.gEnd) {{
        if (sg.kind === 'gap')
          return sg.vStart + (sg.vEnd - sg.vStart) / 2;
        const rl = sg.gEnd - sg.gStart + 1;
        const rel = rl > 0 ? (pos - sg.gStart) / rl : 0;
        return sg.vStart + rel * (sg.vEnd - sg.vStart);
      }}
    }}
    return 0;
  }};
  const exonSegs = segs.filter(s => s.kind === 'exon');
  const MIN_GAP_V = totalVirtual * 0.055;
  const ticks = [];
  for (let i = 0; i < exonSegs.length; i++) {{
    const e = exonSegs[i];
    if (i > 0 && ticks.length && e.vStart - ticks[ticks.length-1].v < MIN_GAP_V) continue;
    ticks.push({{ pos: e.gStart, v: e.vStart, major: i === 0 }});
  }}
  const last = exonSegs[exonSegs.length - 1];
  if (last) {{
    if (ticks.length && last.vEnd - ticks[ticks.length-1].v < MIN_GAP_V) ticks.pop();
    ticks.push({{ pos: last.gEnd, v: last.vEnd, major: true }});
  }}
  return {{ mode, lo, hi, totalVirtual, posToV, ticks }};
}}

/* ============================================================
 * State + DOM refs
 * ============================================================ */
let axis = buildSharedAxis(STRUCTURES, 'compact', true);
let viewport = [0, axis.totalVirtual];
const stackEl   = document.getElementById('iso-stack');
const ticksEl   = document.getElementById('ticks');
const minimapEl = document.getElementById('minimap');
const readoutEl = document.getElementById('viewport-readout');
const tooltipEl = document.getElementById('tooltip');

const fmt = (n) => n.toLocaleString();

function attachTooltip(el, title, subtitle) {{
  el.addEventListener('mouseenter', () => {{
    tooltipEl.innerHTML =
      '<div class="tooltip-title">' + title + '</div>' +
      (subtitle ? '<div class="tooltip-subtitle">' + subtitle + '</div>' : '');
    tooltipEl.style.display = 'block';
  }});
  el.addEventListener('mousemove', (e) => {{
    tooltipEl.style.left = (e.clientX + 12) + 'px';
    tooltipEl.style.top  = (e.clientY + 12) + 'px';
  }});
  el.addEventListener('mouseleave', () => {{ tooltipEl.style.display = 'none'; }});
}}

/* ============================================================
 * Header (from "primary" structure — MANE Select if any, else first)
 * ============================================================ */
function pickPrimary() {{
  for (const s of STRUCTURES) if (s.isManeSelect) return s;
  for (const s of STRUCTURES) if (s.isCanonical) return s;
  return STRUCTURES[0];
}}
function renderHeader() {{
  const s = pickPrimary();
  if (!s) return;
  const bits = [];
  if (s.geneName) bits.push('<b>' + s.geneName + '</b>');
  if (s.domainId) bits.push('domain=' + s.domainId);
  bits.push(STRUCTURES.length + ' isoforms');
  let header = bits.join(' · ');
  if (s.externalLink) {{
    header += ' <a href="' + s.externalLink + '" target="_blank" rel="noopener">↗ link</a>';
  }}
  document.getElementById('header').innerHTML = header;
  const metaBits = [];
  metaBits.push('chrom <b>' + (s.chrom || '?') + '</b> strand <b>' + s.strand + '</b>');
  document.getElementById('meta').innerHTML = metaBits.join(' · ');
}}

function renderStrandArrow(strand) {{
  const wrap = document.getElementById('strand-arrow');
  wrap.innerHTML = '';
  const NS = 'http://www.w3.org/2000/svg';
  const svg = document.createElementNS(NS, 'svg');
  svg.setAttribute('viewBox', '0 0 1000 22');
  svg.setAttribute('preserveAspectRatio', 'none');
  const labelPad = 38;
  const yMid = 11;
  const stroke = '#94a3b8';
  const shaft = document.createElementNS(NS, 'line');
  shaft.setAttribute('x1', labelPad);
  shaft.setAttribute('x2', 1000 - labelPad);
  shaft.setAttribute('y1', yMid); shaft.setAttribute('y2', yMid);
  shaft.setAttribute('stroke', stroke); shaft.setAttribute('stroke-width', '1.5');
  svg.appendChild(shaft);
  const headSize = 8;
  const head = document.createElementNS(NS, 'path');
  if (strand === '+') {{
    const hx = 1000 - labelPad;
    head.setAttribute('d',
      'M' + (hx - headSize) + ' ' + (yMid - headSize) +
      ' L' + hx + ' ' + yMid +
      ' L' + (hx - headSize) + ' ' + (yMid + headSize));
  }} else {{
    const hx = labelPad;
    head.setAttribute('d',
      'M' + (hx + headSize) + ' ' + (yMid - headSize) +
      ' L' + hx + ' ' + yMid +
      ' L' + (hx + headSize) + ' ' + (yMid + headSize));
  }}
  head.setAttribute('fill', 'none'); head.setAttribute('stroke', stroke);
  head.setAttribute('stroke-width', '2'); head.setAttribute('stroke-linejoin', 'round');
  head.setAttribute('stroke-linecap', 'round');
  svg.appendChild(head);
  wrap.appendChild(svg);
  const mkLabel = (txt, side) => {{
    const s = document.createElement('span');
    s.className = 'end-label ' + side;
    s.textContent = txt;
    wrap.appendChild(s);
  }};
  if (strand === '+') {{ mkLabel("5'", 'left'); mkLabel("3'", 'right'); }}
  else                 {{ mkLabel("3'", 'left'); mkLabel("5'", 'right'); }}
}}

/* ============================================================
 * Mapping helpers
 * ============================================================ */
function vToVisPct(v) {{
  const [a, b] = viewport;
  const span = Math.max(1, b - a);
  return ((Math.max(a, Math.min(b, v)) - a) / span) * 100;
}}
function posToVisPct(pos) {{ return vToVisPct(axis.posToV(pos)); }}
function vToFullPct(v) {{ return (v / Math.max(1, axis.totalVirtual)) * 100; }}
function posToFullPct(pos) {{ return vToFullPct(axis.posToV(pos)); }}

function vToGenomic(v) {{
  const ts = axis.ticks;
  if (!ts.length) return axis.lo;
  if (v <= ts[0].v) return ts[0].pos;
  if (v >= ts[ts.length-1].v) return ts[ts.length-1].pos;
  for (let i = 0; i < ts.length - 1; i++) {{
    if (v >= ts[i].v && v <= ts[i+1].v) {{
      const span = ts[i+1].v - ts[i].v;
      if (span <= 0) return ts[i].pos;
      const rel = (v - ts[i].v) / span;
      return ts[i].pos + rel * (ts[i+1].pos - ts[i].pos);
    }}
  }}
  return ts[ts.length-1].pos;
}}

function updateReadout() {{
  const [vLo, vHi] = viewport;
  const a = Math.round(vToGenomic(vLo)), b = Math.round(vToGenomic(vHi));
  readoutEl.textContent = fmt(a) + '–' + fmt(b)
    + '  (' + fmt(b - a + 1) + ' bp visible)';
}}

/* ============================================================
 * Build the isoform stack (one row per structure).
 * ============================================================ */
function buildStack() {{
  stackEl.innerHTML = '';
  for (const s of STRUCTURES) {{
    const row = document.createElement('div');
    row.className = 'iso-row' + (s.isCanonical || s.isManeSelect ? ' canonical' : '');
    const label = document.createElement('div');
    label.className = 'iso-label';
    const inputId = s.inputId || s.transcriptId || '';
    let html = '';
    if (s.transcriptId)
      html += '<span class="tx-id">' + s.transcriptId + '</span>';
    if (s.proteinId)
      html += '<span class="pr-id">' + s.proteinId + '</span>';
    if (s.isManeSelect)
      html += '<span class="badge mane">MANE</span>';
    else if (s.isCanonical)
      html += '<span class="badge">Canonical</span>';
    label.innerHTML = html || ('<span class="tx-id">' + inputId + '</span>');
    const plot = document.createElement('div');
    plot.className = 'iso-plot';
    plot.dataset.idx = String(STRUCTURES.indexOf(s));
    const baseline = document.createElement('div');
    baseline.className = 'baseline';
    plot.appendChild(baseline);
    row.appendChild(label);
    row.appendChild(plot);
    stackEl.appendChild(row);
    attachGestures(plot);
  }}
}}

/* ============================================================
 * Per-isoform render
 * ============================================================ */
function renderIsoform(plot, s) {{
  while (plot.children.length > 1) plot.removeChild(plot.lastChild);
  const showUtr = document.getElementById('show-utr').checked;
  let cdsLo = Infinity, cdsHi = -Infinity;
  for (const f of s.features) if (f.type === 'cds') {{
    if (f.start < cdsLo) cdsLo = f.start;
    if (f.end   > cdsHi) cdsHi = f.end;
  }}
  const visible = s.features.filter(f => {{
    if (!showUtr && (f.type === '5utr' || f.type === '3utr')) return false;
    if (!showUtr && f.type === 'intron')
      return f.end >= cdsLo && f.start <= cdsHi;
    return true;
  }});
  for (const f of visible) {{
    const left  = posToVisPct(f.start);
    const right = posToVisPct(f.end);
    const w = Math.max(0.10, right - left);
    if (left >= 100 || right <= 0) continue;
    if (f.type === 'intron') {{
      const inDomain = f.plotGroup === 'intron_domain_span';
      const intron = document.createElement('div');
      intron.className = 'feat intron' + (inDomain ? ' in-domain' : '');
      intron.style.left = left + '%'; intron.style.width = w + '%';
      const nMarks = w > 20 ? 3 : w > 8 ? 2 : 1;
      for (let k = 0; k < nMarks; k++) {{
        const cx = ((k + 1) / (nMarks + 1)) * 100;
        const svg = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
        svg.setAttribute('viewBox', '0 0 24 24');
        svg.setAttribute('width', '12'); svg.setAttribute('height', '12');
        svg.setAttribute('class', 'chevron');
        svg.style.position = 'absolute';
        svg.style.left = cx + '%'; svg.style.top = '50%';
        svg.style.transform = 'translate(-50%, -50%)';
        const path = document.createElementNS('http://www.w3.org/2000/svg', 'path');
        path.setAttribute('d', s.strand === '+' ? 'M8 5l8 7-8 7' : 'M16 5l-8 7 8 7');
        path.setAttribute('fill', 'none');
        path.setAttribute('stroke', 'currentColor');
        path.setAttribute('stroke-width', '3');
        path.setAttribute('stroke-linecap', 'round');
        path.setAttribute('stroke-linejoin', 'round');
        svg.appendChild(path);
        intron.appendChild(svg);
      }}
      plot.appendChild(intron);
      continue;
    }}
    const isUtr = f.type === '5utr' || f.type === '3utr';
    const isCdsDomain = f.type === 'cds' && f.plotGroup === 'CDS_domain';
    const div = document.createElement('div');
    div.className = 'feat ' + (isUtr ? 'utr' : (isCdsDomain ? 'cds-domain' : 'cds'));
    div.style.left = left + '%'; div.style.width = w + '%';
    const featLabel = f.type === '5utr' ? "5' UTR"
      : f.type === '3utr' ? "3' UTR"
      : (isCdsDomain ? 'CDS (in domain)' : 'CDS');
    const aaTxt = (f.type === 'cds' && f.aaStart && f.aaEnd)
      ? ' · aa ' + f.aaStart + '–' + f.aaEnd : '';
    attachTooltip(div,
      featLabel + ' · ' + fmt(f.end - f.start + 1) + ' bp',
      (s.transcriptId || '') + ' · ' + fmt(f.start) + '–' + fmt(f.end) + aaTxt);
    plot.appendChild(div);
  }}
}}

function renderAll() {{
  const mode = document.querySelector('input[name=mode]:checked').value;
  axis = buildSharedAxis(STRUCTURES, mode, document.getElementById('show-utr').checked);
  viewport = [0, axis.totalVirtual];
  // (Re-)build the stack so per-row plots get freshly-bound gesture handlers.
  buildStack();
  renderAllIsoforms();
  renderTicks();
  renderMinimap();
  updateReadout();
}}

function renderAllIsoforms() {{
  const rows = stackEl.querySelectorAll('.iso-plot');
  rows.forEach((plot, i) => renderIsoform(plot, STRUCTURES[i]));
}}

function renderTicks() {{
  ticksEl.innerHTML = '';
  const [vLo, vHi] = viewport;
  const tList = [];
  tList.push({{ v: vLo, pos: vToGenomic(vLo), major: true }});
  for (const t of axis.ticks) {{
    if (t.v > vLo + (vHi - vLo) * 0.04 && t.v < vHi - (vHi - vLo) * 0.04) {{
      tList.push(t);
    }}
  }}
  tList.push({{ v: vHi, pos: vToGenomic(vHi), major: true }});
  const seen = new Set();
  tList.forEach(t => {{
    const xp = vToVisPct(t.v);
    const key = Math.round(xp * 10);
    if (seen.has(key)) return;
    seen.add(key);
    const el = document.createElement('div');
    el.className = 'tick' + (t.major ? ' major' : '');
    el.style.left = xp + '%';
    // Keep edge ticks inside the track: anchor the first by its left edge and
    // the last by its right edge; interior ticks stay centred on the mark.
    el.style.transform = xp <= 0.5 ? 'translateX(0)'
                       : xp >= 99.5 ? 'translateX(-100%)' : 'translateX(-50%)';
    el.textContent = fmt(Math.round(t.pos));
    ticksEl.appendChild(el);
  }});
}}

/* ============================================================
 * Minimap — overlays a thin row per isoform at the full extent.
 * ============================================================ */
function renderMinimap() {{
  minimapEl.innerHTML = '';
  const n = STRUCTURES.length;
  const rowH = Math.max(2, Math.floor((36 - 4) / n));
  STRUCTURES.forEach((s, i) => {{
    const rowEl = document.createElement('div');
    rowEl.style.position = 'absolute';
    rowEl.style.top = (2 + i * rowH) + 'px';
    rowEl.style.height = rowH + 'px';
    rowEl.style.left = '0'; rowEl.style.right = '0';
    for (const f of s.features) {{
      if (f.type === 'intron' || f.type === '5utr' || f.type === '3utr') continue;
      const isCdsDomain = f.plotGroup === 'CDS_domain';
      const div = document.createElement('div');
      div.className = 'mini-feat ' + (isCdsDomain ? 'cds-domain' : 'cds');
      div.style.position = 'absolute';
      const l = posToFullPct(f.start), r = posToFullPct(f.end);
      div.style.left = l + '%';
      div.style.width = Math.max(0.05, r - l) + '%';
      div.style.top = '0'; div.style.height = '100%';
      rowEl.appendChild(div);
    }}
    minimapEl.appendChild(rowEl);
  }});
  // Viewport rectangle.
  const vp = document.createElement('div');
  vp.className = 'viewport'; vp.id = 'viewport-rect';
  const [vLo, vHi] = viewport;
  vp.style.left  = vToFullPct(vLo) + '%';
  vp.style.width = (vToFullPct(vHi) - vToFullPct(vLo)) + '%';
  const handleL = document.createElement('div');
  handleL.className = 'viewport-handle left';
  const handleR = document.createElement('div');
  handleR.className = 'viewport-handle right';
  vp.appendChild(handleL); vp.appendChild(handleR);
  minimapEl.appendChild(vp);
  attachMinimapDrag(vp,      'pan');
  attachMinimapDrag(handleL, 'left');
  attachMinimapDrag(handleR, 'right');
}}

/* ============================================================
 * Gestures: per-row box-zoom, wheel-zoom, dblclick reset, shift-drag pan.
 * ============================================================ */
function attachGestures(plot) {{
  let boxState = null, panState = null, didDrag = false;
  let shadeL, shadeR, boxRect, boxReadout;

  function ensureBoxNodes() {{
    if (shadeL && shadeL.isConnected && shadeR && shadeR.isConnected
        && boxRect && boxRect.isConnected && boxReadout && boxReadout.isConnected) return;
    if (shadeL && shadeL.parentNode)    shadeL.parentNode.removeChild(shadeL);
    if (shadeR && shadeR.parentNode)    shadeR.parentNode.removeChild(shadeR);
    if (boxRect && boxRect.parentNode)  boxRect.parentNode.removeChild(boxRect);
    if (boxReadout && boxReadout.parentNode) boxReadout.parentNode.removeChild(boxReadout);
    shadeL = document.createElement('div'); shadeL.className = 'box-shade';
    shadeR = document.createElement('div'); shadeR.className = 'box-shade';
    boxRect = document.createElement('div'); boxRect.className = 'box-rect';
    boxReadout = document.createElement('div'); boxReadout.className = 'box-readout';
    shadeL.style.display = shadeR.style.display = 'none';
    boxRect.style.display = boxReadout.style.display = 'none';
    plot.appendChild(shadeL); plot.appendChild(shadeR);
    plot.appendChild(boxRect); plot.appendChild(boxReadout);
  }}
  function drawBox(px0, px1, mpW) {{
    ensureBoxNodes();
    const lo = Math.max(0, Math.min(px0, px1));
    const hi = Math.min(mpW, Math.max(px0, px1));
    shadeL.style.left  = '0px';      shadeL.style.width = lo + 'px';     shadeL.style.display = 'block';
    shadeR.style.left  = hi + 'px';  shadeR.style.width = Math.max(0, mpW - hi) + 'px'; shadeR.style.display = 'block';
    boxRect.style.left = lo + 'px';  boxRect.style.width = Math.max(0, hi - lo) + 'px'; boxRect.style.display = 'block';
    const [a, b] = viewport; const w = b - a;
    const bpLo = Math.round(vToGenomic(a + (lo / mpW) * w));
    const bpHi = Math.round(vToGenomic(a + (hi / mpW) * w));
    boxReadout.textContent = fmt(bpLo) + '–' + fmt(bpHi)
      + '  (' + fmt(Math.max(0, bpHi - bpLo + 1)) + ' bp)';
    const cx = Math.max(0, Math.min(mpW, (lo + hi) / 2));
    boxReadout.style.left = cx + 'px'; boxReadout.style.transform = 'translateX(-50%)';
    boxReadout.style.display = 'block';
  }}
  function clearBox() {{
    if (shadeL)     shadeL.style.display     = 'none';
    if (shadeR)     shadeR.style.display     = 'none';
    if (boxRect)    boxRect.style.display    = 'none';
    if (boxReadout) boxReadout.style.display = 'none';
  }}

  plot.addEventListener('mousedown', (e) => {{
    if (e.button !== 0) return;
    const rect = plot.getBoundingClientRect();
    const x = e.clientX - rect.left;
    if (e.shiftKey) {{
      panState = {{ startX: e.clientX, mpW: rect.width,
                    startVLo: viewport[0], startVHi: viewport[1] }};
      plot.classList.add('dragging');
    }} else {{
      boxState = {{ startPx: x, endPx: x, mpW: rect.width, mpLeft: rect.left }};
      didDrag = false;
      drawBox(x, x, rect.width);
    }}
    e.preventDefault();
  }});

  window.addEventListener('mousemove', (e) => {{
    if (boxState) {{
      const x = Math.max(0, Math.min(boxState.mpW, e.clientX - boxState.mpLeft));
      boxState.endPx = x;
      if (Math.abs(x - boxState.startPx) >= 1) didDrag = true;
      drawBox(boxState.startPx, boxState.endPx, boxState.mpW);
      return;
    }}
    if (!panState) return;
    const dx = e.clientX - panState.startX;
    const w = panState.startVHi - panState.startVLo;
    const dV = -(dx / panState.mpW) * w;
    let a = Math.max(0, Math.min(axis.totalVirtual - w, panState.startVLo + dV));
    viewport = [a, a + w];
    renderAllIsoforms(); renderTicks(); renderMinimap(); updateReadout();
  }});

  window.addEventListener('mouseup', () => {{
    if (boxState) {{
      const {{ startPx, endPx, mpW }} = boxState;
      boxState = null;
      clearBox();
      if (didDrag) {{
        const lo = Math.max(0, Math.min(startPx, endPx));
        const hi = Math.min(mpW, Math.max(startPx, endPx));
        if (hi - lo >= 6) {{
          const [a, b] = viewport;
          const w = b - a;
          const newA = a + (lo / mpW) * w;
          const newB = a + (hi / mpW) * w;
          const minSpan = Math.max(20, axis.totalVirtual * 0.0005);
          if (newB - newA >= minSpan) {{
            viewport = [newA, newB];
            renderAllIsoforms(); renderTicks(); renderMinimap(); updateReadout();
          }}
        }}
      }}
      didDrag = false;
    }}
    if (panState) {{ plot.classList.remove('dragging'); panState = null; }}
  }});

  plot.addEventListener('wheel', (e) => {{
    e.preventDefault();
    const rect = plot.getBoundingClientRect();
    const relX = (e.clientX - rect.left) / rect.width;
    const [a, b] = viewport;
    const w = b - a;
    const center = a + relX * w;
    const factor = e.deltaY < 0 ? 0.8 : 1.25;
    let newW = Math.max(20, Math.min(axis.totalVirtual, w * factor));
    let newA = center - relX * newW;
    let newB = newA + newW;
    if (newA < 0)                  {{ newA = 0; newB = newW; }}
    if (newB > axis.totalVirtual)  {{ newB = axis.totalVirtual; newA = newB - newW; }}
    viewport = [newA, newB];
    renderAllIsoforms(); renderTicks(); renderMinimap(); updateReadout();
  }}, {{ passive: false }});

  plot.addEventListener('dblclick', () => {{
    viewport = [0, axis.totalVirtual];
    renderAllIsoforms(); renderTicks(); renderMinimap(); updateReadout();
  }});
}}

/* ============================================================
 * Minimap drag (pan / resize).
 * ============================================================ */
let dragState = null;
function attachMinimapDrag(el, kind) {{
  el.addEventListener('mousedown', (e) => {{
    const rect = minimapEl.getBoundingClientRect();
    dragState = {{ kind, startX: e.clientX, mmW: rect.width,
                   startVLo: viewport[0], startVHi: viewport[1] }};
    el.classList.add('dragging');
    e.preventDefault(); e.stopPropagation();
  }});
}}
minimapEl.addEventListener('click', (e) => {{
  if (e.target.closest('.viewport')) return;
  const rect = minimapEl.getBoundingClientRect();
  const rel = (e.clientX - rect.left) / rect.width;
  const v = rel * axis.totalVirtual;
  const w = viewport[1] - viewport[0];
  let a = v - w / 2, b = v + w / 2;
  if (a < 0)                   {{ a = 0; b = w; }}
  if (b > axis.totalVirtual)   {{ b = axis.totalVirtual; a = b - w; }}
  viewport = [a, b];
  renderAllIsoforms(); renderTicks(); renderMinimap(); updateReadout();
}});
window.addEventListener('mousemove', (e) => {{
  if (!dragState) return;
  const dx = e.clientX - dragState.startX;
  const dV = (dx / dragState.mmW) * axis.totalVirtual;
  let [a, b] = [dragState.startVLo, dragState.startVHi];
  const minSpan = Math.max(20, axis.totalVirtual * 0.005);
  if (dragState.kind === 'pan') {{
    const w = b - a;
    a = Math.max(0, Math.min(axis.totalVirtual - w, a + dV));
    b = a + w;
  }} else if (dragState.kind === 'left') {{
    a = Math.max(0, Math.min(b - minSpan, a + dV));
  }} else if (dragState.kind === 'right') {{
    b = Math.min(axis.totalVirtual, Math.max(a + minSpan, b + dV));
  }}
  viewport = [a, b];
  renderAllIsoforms(); renderTicks(); renderMinimap(); updateReadout();
}});
window.addEventListener('mouseup', () => {{
  if (dragState) {{
    document.querySelectorAll('.dragging').forEach(el => el.classList.remove('dragging'));
    dragState = null;
  }}
}});

document.getElementById('btn-reset').addEventListener('click', () => {{
  viewport = [0, axis.totalVirtual];
  renderAllIsoforms(); renderTicks(); renderMinimap(); updateReadout();
}});
document.getElementById('btn-fit-domain').addEventListener('click', () => {{
  // Use the domain envelope from any structure that has one (they should
  // largely agree across isoforms of the same gene).
  for (const s of STRUCTURES) {{
    if (!s.domains || !s.domains.length) continue;
    const d = s.domains[0];
    if (d.genomicStart == null || d.genomicEnd == null) continue;
    const a = axis.posToV(d.genomicStart);
    const b = axis.posToV(d.genomicEnd);
    const pad = (b - a) * 0.10;
    let vLo = Math.max(0, a - pad);
    let vHi = Math.min(axis.totalVirtual, b + pad);
    if (vHi - vLo < 20) {{ vLo = Math.max(0, vLo - 10); vHi = vLo + 40; }}
    viewport = [vLo, vHi];
    renderAllIsoforms(); renderTicks(); renderMinimap(); updateReadout();
    return;
  }}
}});

/* ============================================================
 * Boot + auto-resize handshake (for embedding in an iframe).
 * ============================================================ */
renderHeader();
renderStrandArrow(pickPrimary().strand);
renderAll();
document.querySelectorAll('input[name=mode]').forEach(r =>
  r.addEventListener('change', renderAll));
document.getElementById('show-utr').addEventListener('change', renderAll);

(function () {{
  if (window.parent === window) return;
  const broadcast = () => {{
    const h = Math.max(
      document.documentElement.scrollHeight,
      document.body.scrollHeight,
      document.documentElement.offsetHeight,
      document.body.offsetHeight
    );
    window.parent.postMessage({{ type: 'fastCDS-resize', height: h }}, '*');
  }};
  window.addEventListener('load', () => setTimeout(broadcast, 30));
  document.querySelectorAll('input').forEach(el =>
    el.addEventListener('change', () => setTimeout(broadcast, 30)));
  if (typeof ResizeObserver !== 'undefined') {{
    new ResizeObserver(broadcast).observe(document.body);
  }}
}})();
</script>
</body>
</html>
"""


def _render_stack_to_string(segs_by_id: dict[str, list], *,
                            link_template: str | None = None,
                            plot_height: int = 40) -> str:
    """Build the multi-isoform HTML string."""
    payloads = _payloads(segs_by_id)
    if not payloads:
        raise ValueError("no isoforms with mappable segments")

    if link_template:
        for p in payloads:
            try:
                fields = {
                    "protein_id":    p.get("proteinId", ""),
                    "gene_name":     p.get("geneName", ""),
                    "transcript_id": p.get("transcriptId", ""),
                    "chrom":         p.get("chrom", ""),
                }
                p["externalLink"] = link_template.format(**fields)
            except (KeyError, IndexError):
                pass

    # Header title from the primary (MANE / canonical / first) structure.
    primary = next((p for p in payloads if p.get("isManeSelect")), None) \
        or next((p for p in payloads if p.get("isCanonical")), None) \
        or payloads[0]
    title_bits = [primary.get("geneName", ""), primary.get("domainId", ""),
                  f"{len(payloads)} isoforms"]
    title = " · ".join(b for b in title_bits if b) or "fastCDS"

    return HTML_TEMPLATE.format(
        title_safe=html.escape(title),
        payloads_json=json.dumps(payloads),
        compact_intron_bp=COMPACT_INTRON_BP,
        plot_height=int(plot_height),
    )


def render_interactive_html_stack(segs_by_id: dict[str, list],
                               out: str,
                               *,
                               link_template: str | None = None,
                               plot_height: int = 40) -> None:
    """Write a multi-isoform stack HTML viewer to ``out``.

    All isoforms share a single axis built from the union of their
    features, so exon presence/absence aligns across rows.

    Parameters
    ----------
    segs_by_id : dict[str, list[Segment]]
        One entry per isoform — use ``_segments_from_dataframe(df)``
        or ``load_isoform_tsv(path)`` to build this.
    out : str
        Output HTML path.
    link_template : str, optional
        Same as the single-isoform viewer (``{protein_id}`` etc.).
    plot_height : int
        Per-isoform track height in pixels (default 40, suited for
        4–8 isoforms — bump for taller rows).
    """
    html_out = _render_stack_to_string(
        segs_by_id, link_template=link_template, plot_height=plot_height,
    )
    with open(out, "w", encoding="utf-8") as f:
        f.write(html_out)


_JUPYTER_STACK_SEQ = 0


def render_interactive_jupyter_stack(segs_by_id: dict[str, list],
                                  *,
                                  height: int | None = None,
                                  plot_height: int = 40,
                                  link_template: str | None = None):
    """Embed the multi-isoform stack viewer in a Jupyter notebook."""
    from IPython.display import HTML  # noqa: F401 — needed only in notebooks
    html_str = _render_stack_to_string(
        segs_by_id, link_template=link_template, plot_height=plot_height,
    )
    safe = (html_str
            .replace("&", "&amp;")
            .replace('"', "&quot;"))

    global _JUPYTER_STACK_SEQ
    _JUPYTER_STACK_SEQ += 1
    iframe_id = f"p2e-tfregdb2-stack-{_JUPYTER_STACK_SEQ}"
    # Generous fallback; postMessage will tighten to the actual height.
    initial_h = int(height) if height is not None else max(
        320, 200 + plot_height * len(segs_by_id) + 60)

    iframe = (
        f'<iframe id="{iframe_id}" srcdoc="{safe}" '
        f'style="width: 100%; height: {initial_h}px; border: none; '
        f'border-radius: 6px; display: block;" '
        f'sandbox="allow-scripts"></iframe>'
    )
    pin = height is not None
    listener = "" if pin else f"""
<script>
(function () {{
  var ifr = document.getElementById({iframe_id!r});
  if (!ifr) return;
  window.addEventListener('message', function (e) {{
    if (!e.data || e.data.type !== 'fastCDS-resize') return;
    if (e.source !== ifr.contentWindow) return;
    var h = Math.max(200, Math.min(4000, e.data.height + 4));
    ifr.style.height = h + 'px';
  }});
}})();
</script>"""
    return HTML(
        f'<div class="fastCDS-viewer" style="margin: 0;">{iframe}{listener}</div>'
    )
