"""Standalone interactive HTML renderer for a single isoform.

Draws a slate gene track with chevron introns, lane-packed coloured domain
overlays, and a compact 80 bp intron mode, plus a minimap with box-zoom,
drag-pan and wheel-zoom.

Implemented as a single self-contained HTML file (vanilla JS, inline CSS, no
React, no CDN) so the output works offline and stays embed friendly.
"""

from __future__ import annotations

import html
import json
from typing import Any


COMPACT_INTRON_BP = 80  # mirrors GenomicStructureView.tsx

# Mirrors the React palette — keep in sync if upstream changes.
DOMAIN_COLOR = {
    "AD":           "#059669",
    "RD":           "#e11d48",
    "Bifunctional": "#7c3aed",
    "DNA":          "#2563eb",
    "Oligo":        "#d97706",
    "Other":        "#64748b",
}


def _segments_to_payload(segs) -> dict:
    """Translate fastCDS Segment objects into the JS payload."""
    if not segs:
        return {}
    s0 = segs[0]
    features = []
    domain_g_start = None
    domain_g_end = None
    aa_min = None
    aa_max = None
    for s in segs:
        ft_map = {
            "CDS": "cds",
            "intron": "intron",
            "five_prime_UTR": "5utr",
            "three_prime_UTR": "3utr",
        }
        ft = ft_map.get(s.feature_type)
        if ft is None:
            continue
        aa_s = None
        aa_e = None
        if hasattr(s, "aa_start_encoded"):
            try:
                aa_s = int(s.aa_start_encoded) if s.aa_start_encoded not in (None, "", "NA") else None
                aa_e = int(s.aa_end_encoded) if s.aa_end_encoded not in (None, "", "NA") else None
            except (TypeError, ValueError):
                aa_s = aa_e = None
        features.append({
            "type": ft,
            "start": int(s.feature_genomic_start),
            "end": int(s.feature_genomic_end),
            "order": int(getattr(s, "feature_order_transcript", 0) or 0),
            "aaStart": aa_s,
            "aaEnd": aa_e,
            "plotGroup": getattr(s, "plot_group", "") or "",
        })
        if s.plot_group in ("CDS_domain", "intron_domain_span"):
            if domain_g_start is None or s.feature_genomic_start < domain_g_start:
                domain_g_start = int(s.feature_genomic_start)
            if domain_g_end is None or s.feature_genomic_end > domain_g_end:
                domain_g_end = int(s.feature_genomic_end)
            if s.plot_group == "CDS_domain" and aa_s is not None and aa_e is not None:
                aa_min = aa_s if aa_min is None else min(aa_min, aa_s)
                aa_max = aa_e if aa_max is None else max(aa_max, aa_e)

    domains = []
    if domain_g_start is not None and domain_g_end is not None:
        domains.append({
            "type": "Other",
            "source": getattr(s0, "domain_id", None) or None,
            "aaStart": aa_min,
            "aaEnd": aa_max,
            "genomicStart": domain_g_start,
            "genomicEnd": domain_g_end,
        })

    return {
        "transcriptId":  getattr(s0, "transcript_id", "") or "",
        "proteinId":     getattr(s0, "protein_id", "") or "",
        "geneName":      getattr(s0, "gene_name", "") or "",
        "domainId":      getattr(s0, "domain_id", "") or "",
        "chrom":         getattr(s0, "chrom", "") or "",
        "strand":        s0.strand or "+",
        "isManeSelect":  getattr(s0, "is_mane_select", "NA") == "true",
        "isCanonical":   getattr(s0, "is_ensembl_canonical", "NA") == "true",
        "features":      features,
        "domains":       domains,
    }


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
  .meta {{ font-size: 12px; color: #475569; margin-bottom: 14px; }}
  .meta b {{ color: #0f172a; }}
  .controls {{ display: flex; flex-wrap: wrap; gap: 12px; align-items: center; margin-bottom: 14px; font-size: 12px; }}
  .controls label {{ display: inline-flex; align-items: center; cursor: pointer; user-select: none; }}
  .controls input[type=radio], .controls input[type=checkbox] {{ margin-right: 4px; }}
  .controls button {{
    font-size: 11px; padding: 3px 10px; border: 1px solid #cbd5e1;
    background: #f8fafc; border-radius: 4px; cursor: pointer;
  }}
  .controls button:hover {{ background: #e2e8f0; }}
  .legend {{
    display: flex; flex-wrap: wrap; gap: 14px; font-size: 11px;
    color: #475569; margin-bottom: 12px;
  }}
  .legend .swatch {{
    display: inline-block; vertical-align: middle; margin-right: 4px;
    border-radius: 2px; border: 1px solid rgba(0,0,0,0.1);
  }}
  .legend .sw-utr      {{ width: 14px; height: 8px;  background: #f0c078; }}
  .legend .sw-cds      {{ width: 14px; height: 16px; background: #1f77b4; }}
  .legend .sw-cdsd     {{ width: 14px; height: 16px; background: #d62728; }}
  .legend .sw-line     {{ width: 16px; height: 2px; background: #64748b; vertical-align: middle; display: inline-block; margin-right: 4px; }}
  .legend .sw-redline  {{ width: 16px; height: 2px; background: #d62728; vertical-align: middle; display: inline-block; margin-right: 4px; }}
  .strand-arrow {{
    width: 100%; height: 22px; position: relative; margin-bottom: 6px;
    color: #475569;
  }}
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
  /* Main track — height is configurable via the renderer's plot_height kwarg. */
  .main-plot {{
    position: relative; width: 100%; height: {plot_height}px;
    background: linear-gradient(180deg, #fbfdff 0%, #f3f7fb 100%);
    border: 1px solid #e2e8f0; border-radius: 4px;
    overflow: hidden; user-select: none; cursor: crosshair;
  }}
  .main-plot.dragging {{ cursor: grabbing; }}
  .main-plot .baseline {{
    position: absolute; left: 0; right: 0; top: 50%; height: 1px;
    background: #cbd5e1;
  }}
  /* Feature heights/positions are percentages of .main-plot so the track
     scales cleanly when plot_height is increased for Jupyter embedding. */
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
  /* Box-zoom overlay (vCRE-style: shade outside the drag selection). */
  .box-shade {{
    position: absolute; top: 0; bottom: 0;
    background: rgba(15, 23, 42, 0.40);
    pointer-events: none; z-index: 50;
  }}
  .box-rect {{
    position: absolute; top: 0; bottom: 0;
    border-left:  2px solid #4f46e5;
    border-right: 2px solid #4f46e5;
    background: rgba(79, 70, 229, 0.10);
    box-shadow: inset 0 0 0 1px rgba(255,255,255,0.4);
    pointer-events: none; z-index: 51;
  }}
  .box-readout {{
    position: absolute; top: 4px;
    background: #4f46e5; color: white;
    font-size: 10.5px; font-weight: 600;
    padding: 2px 6px; border-radius: 3px;
    pointer-events: none; z-index: 52;
    white-space: nowrap;
    box-shadow: 0 1px 4px rgba(15, 23, 42, 0.25);
  }}
  .ticks {{
    position: relative; height: 18px; width: 100%; margin-top: 6px;
    font-size: 10px; color: #64748b;
  }}
  .tick {{
    position: absolute; top: 0; transform: translateX(-50%);
    border-left: 1px solid #cbd5e1; padding-left: 2px; padding-top: 4px;
    white-space: nowrap;
  }}
  .tick.major {{ color: #0f172a; }}
  /* Minimap — overview at full extent, with a draggable viewport rect. */
  .minimap-wrap {{
    margin-top: 14px; padding-top: 8px; border-top: 1px solid #e2e8f0;
  }}
  .minimap-label {{
    font-size: 11px; color: #64748b; margin-bottom: 4px;
    display: flex; justify-content: space-between;
  }}
  .minimap {{
    position: relative; width: 100%; height: 32px;
    background: #f8fafc; border: 1px solid #e2e8f0; border-radius: 3px;
    overflow: hidden; cursor: pointer;
  }}
  .mini-baseline {{
    position: absolute; left: 0; right: 0; top: 50%; height: 1px;
    background: #cbd5e1;
  }}
  .mini-feat {{ position: absolute; border-radius: 1px; }}
  .mini-feat.cds        {{ background: #1f77b4; height: 12px; top: 10px; }}
  .mini-feat.cds-domain {{ background: #d62728; height: 12px; top: 10px; }}
  .mini-feat.utr        {{ background: #f0c078; height: 5px; top: 14px; }}
  .mini-feat.intron-d   {{ background: #d62728; height: 1px; top: 16px; }}
  .viewport {{
    position: absolute; top: 0; height: 100%;
    box-sizing: border-box;   /* keep the 2px side borders inside width:100% so
                                 the right edge isn't clipped thinner than left */
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
  .tooltip {{
    position: fixed; background: white; border: 1px solid #cbd5e1; border-radius: 6px;
    padding: 8px 10px; font-size: 12px; box-shadow: 0 4px 12px rgba(15,23,42,0.15);
    pointer-events: none; z-index: 10; max-width: 260px;
  }}
  .tooltip-title {{ font-weight: 600; margin-bottom: 2px; color: #0f172a; }}
  .tooltip-subtitle {{ color: #475569; }}
  h1 {{ font-size: 16px; margin: 0 0 6px 0; }}
  h1 a {{ color: #2563eb; text-decoration: underline; font-size: 12px; margin-left: 8px; font-weight: normal; }}
  .badge {{
    display: inline-block; padding: 2px 6px; border-radius: 999px;
    background: #f1f5f9; color: #475569; font-size: 11px; margin-right: 4px;
  }}
  .badge.green {{ background: #ecfdf5; color: #059669; }}
  .badge.blue  {{ background: #eff6ff; color: #2563eb; }}
  .help-text {{ font-size: 10.5px; color: #94a3b8; margin-top: 6px; }}
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
  <div class="strand-arrow" id="strand-arrow"></div>
  <div class="main-plot" id="main-plot">
    <div class="baseline"></div>
  </div>
  <div class="ticks" id="ticks"></div>
  <div class="minimap-wrap">
    <div class="minimap-label">
      <span>Overview</span><span id="viewport-readout"></span>
    </div>
    <div class="minimap" id="minimap">
      <div class="mini-baseline"></div>
    </div>
    <div class="help-text">Drag on the main plot to box-zoom (shift-drag to pan) · scroll to zoom around the cursor · double-click to reset · drag the indigo minimap rectangle (or its edges) to pan/resize · click empty minimap to recenter</div>
  </div>
</div>
<div class="tooltip" id="tooltip" style="display:none"></div>
<script>
const STRUCTURE = {payload_json};
const DOMAIN_COLOR = {domain_colors_json};
const COMPACT_INTRON_BP = {compact_intron_bp};

/* ============================================================
 * Shared-axis: maps genomic position -> virtual position in [0, totalVirtual]
 * Compact mode collapses introns to COMPACT_INTRON_BP virtual bp.
 * ============================================================ */
function buildSharedAxis(s, mode, includeUtr) {{
  const ranges = [];
  let lo = Infinity, hi = -Infinity;
  for (const f of s.features) {{
    if (f.type === 'intron') continue;
    if (!includeUtr && (f.type === '5utr' || f.type === '3utr')) continue;
    if (f.start < lo) lo = f.start;
    if (f.end > hi) hi = f.end;
    ranges.push({{ start: f.start, end: f.end }});
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
  // Ticks at exon boundaries.
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
let axis = buildSharedAxis(STRUCTURE, 'compact', true);
// viewport: [vLo, vHi] subset of [0, axis.totalVirtual] that's currently
// shown in the main plot. Minimap always shows the full extent.
let viewport = [0, axis.totalVirtual];
const mainPlot   = document.getElementById('main-plot');
const minimapEl  = document.getElementById('minimap');
const ticksEl    = document.getElementById('ticks');
const readoutEl  = document.getElementById('viewport-readout');
const tooltipEl  = document.getElementById('tooltip');

/* ============================================================
 * Tooltip
 * ============================================================ */
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

const fmt = (n) => n.toLocaleString();

/* ============================================================
 * Header + strand arrow
 * ============================================================ */
function renderHeader(s) {{
  const bits = [];
  if (s.geneName)     bits.push('<b>' + s.geneName + '</b>');
  if (s.proteinId)    bits.push(s.proteinId);
  if (s.transcriptId) bits.push(s.transcriptId);
  if (s.domainId)     bits.push('domain=' + s.domainId);
  let header = bits.join(' · ');
  if (s.externalLink) {{
    header += ' <a href="' + s.externalLink + '" target="_blank" rel="noopener">↗ link</a>';
  }}
  document.getElementById('header').innerHTML = header;
  const metaBits = [];
  metaBits.push('chrom <b>' + (s.chrom || '?') + '</b> strand <b>' + s.strand + '</b>');
  if (s.isManeSelect) metaBits.push('<span class="badge green">MANE Select</span>');
  if (s.isCanonical)  metaBits.push('<span class="badge blue">Ensembl canonical</span>');
  document.getElementById('meta').innerHTML = metaBits.join(' · ');
}}

function renderStrandArrow(strand) {{
  const wrap = document.getElementById('strand-arrow');
  wrap.innerHTML = '';
  // SVG arrow that spans the full plot width, with chevrons evenly placed
  // along the line so the transcription direction is obvious at a glance.
  const NS = 'http://www.w3.org/2000/svg';
  const svg = document.createElementNS(NS, 'svg');
  svg.setAttribute('viewBox', '0 0 1000 22');
  svg.setAttribute('preserveAspectRatio', 'none');
  const labelPad = 38;            // px reserved on each side for 5'/3' badges
  const yMid = 11;
  const stroke = '#94a3b8';
  // Main shaft.
  const shaft = document.createElementNS(NS, 'line');
  shaft.setAttribute('x1', labelPad);
  shaft.setAttribute('x2', 1000 - labelPad);
  shaft.setAttribute('y1', yMid);
  shaft.setAttribute('y2', yMid);
  shaft.setAttribute('stroke', stroke);
  shaft.setAttribute('stroke-width', '1.5');
  svg.appendChild(shaft);
  // Big arrowhead at the 3' end.
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
  head.setAttribute('fill', 'none');
  head.setAttribute('stroke', stroke);
  head.setAttribute('stroke-width', '2');
  head.setAttribute('stroke-linejoin', 'round');
  head.setAttribute('stroke-linecap', 'round');
  svg.appendChild(head);
  // Repeated faint chevrons along the shaft pointing in the strand direction.
  const nChev = 6;
  const x0 = labelPad + 30, x1 = (1000 - labelPad) - 30;
  for (let i = 0; i < nChev; i++) {{
    const cx = x0 + (i + 0.5) * (x1 - x0) / nChev;
    const ch = document.createElementNS(NS, 'path');
    if (strand === '+') {{
      ch.setAttribute('d',
        'M' + (cx - 4) + ' ' + (yMid - 4) +
        ' L' + (cx + 2) + ' ' + yMid +
        ' L' + (cx - 4) + ' ' + (yMid + 4));
    }} else {{
      ch.setAttribute('d',
        'M' + (cx + 4) + ' ' + (yMid - 4) +
        ' L' + (cx - 2) + ' ' + yMid +
        ' L' + (cx + 4) + ' ' + (yMid + 4));
    }}
    ch.setAttribute('fill', 'none');
    ch.setAttribute('stroke', '#cbd5e1');
    ch.setAttribute('stroke-width', '1.5');
    ch.setAttribute('stroke-linejoin', 'round');
    ch.setAttribute('stroke-linecap', 'round');
    svg.appendChild(ch);
  }}
  wrap.appendChild(svg);
  // 5' / 3' labels as small pills floating above the shaft.
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
  // virtual-coord -> percent of the main-plot viewport. Outside viewport is clamped.
  const [a, b] = viewport;
  const span = Math.max(1, b - a);
  return ((Math.max(a, Math.min(b, v)) - a) / span) * 100;
}}
function posToVisPct(pos) {{ return vToVisPct(axis.posToV(pos)); }}
function vToFullPct(v) {{
  return (v / Math.max(1, axis.totalVirtual)) * 100;
}}
function posToFullPct(pos) {{ return vToFullPct(axis.posToV(pos)); }}

/* ============================================================
 * Main plot render
 * ============================================================ */
function render() {{
  const mode = document.querySelector('input[name=mode]:checked').value;
  const showUtr = document.getElementById('show-utr').checked;
  axis = buildSharedAxis(STRUCTURE, mode, showUtr);
  // Reset viewport on layout/utr change to "full".
  viewport = [0, axis.totalVirtual];
  renderMain();
  renderMinimap();
  renderTicks();
  updateReadout();
}}

function renderMain() {{
  // Clear children except the baseline.
  while (mainPlot.children.length > 1) mainPlot.removeChild(mainPlot.lastChild);
  const features = STRUCTURE.features.slice().sort((a, b) => a.start - b.start);
  const showUtr = document.getElementById('show-utr').checked;
  let cdsLo = Infinity, cdsHi = -Infinity;
  for (const f of features) if (f.type === 'cds') {{
    if (f.start < cdsLo) cdsLo = f.start;
    if (f.end > cdsHi) cdsHi = f.end;
  }}
  const visible = features.filter(f => {{
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
        svg.setAttribute('width', '12');
        svg.setAttribute('height', '12');
        svg.setAttribute('class', 'chevron');
        svg.style.position = 'absolute';
        svg.style.left = cx + '%';
        svg.style.top = '50%';
        svg.style.transform = 'translate(-50%, -50%)';
        const path = document.createElementNS('http://www.w3.org/2000/svg', 'path');
        path.setAttribute('d', STRUCTURE.strand === '+' ? 'M8 5l8 7-8 7' : 'M16 5l-8 7 8 7');
        path.setAttribute('fill', 'none');
        path.setAttribute('stroke', 'currentColor');
        path.setAttribute('stroke-width', '3');
        path.setAttribute('stroke-linecap', 'round');
        path.setAttribute('stroke-linejoin', 'round');
        svg.appendChild(path);
        intron.appendChild(svg);
      }}
      mainPlot.appendChild(intron);
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
    attachTooltip(
      div,
      featLabel + ' · ' + fmt(f.end - f.start + 1) + ' bp',
      fmt(f.start) + '–' + fmt(f.end) + aaTxt,
    );
    mainPlot.appendChild(div);
  }}
}}

function renderTicks() {{
  ticksEl.innerHTML = '';
  // Ticks at the viewport's left edge, right edge, plus any axis tick inside.
  const [vLo, vHi] = viewport;
  const ticks = [];
  // Always anchor the visible bounds.
  ticks.push({{ v: vLo, pos: vToGenomic(vLo), major: true }});
  for (const t of axis.ticks) {{
    if (t.v > vLo + (vHi - vLo) * 0.04 && t.v < vHi - (vHi - vLo) * 0.04) {{
      ticks.push(t);
    }}
  }}
  ticks.push({{ v: vHi, pos: vToGenomic(vHi), major: true }});
  // De-dupe by xPct.
  const seen = new Set();
  ticks.forEach(t => {{
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

// Inverse of axis.posToV — given a virtual coord, return the genomic pos at
// that point. Approximation good enough for ticks.
function vToGenomic(v) {{
  // We don't carry the segs from buildSharedAxis out; quick re-derive by
  // linear search through axis.ticks.
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

/* ============================================================
 * Minimap render (always full extent; viewport rect on top)
 * ============================================================ */
function renderMinimap() {{
  // Clear except the baseline.
  while (minimapEl.children.length > 1) minimapEl.removeChild(minimapEl.lastChild);
  for (const f of STRUCTURE.features) {{
    if (f.type === 'intron') {{
      if (f.plotGroup !== 'intron_domain_span') continue;
      const div = document.createElement('div');
      div.className = 'mini-feat intron-d';
      const l = posToFullPct(f.start), r = posToFullPct(f.end);
      div.style.left = l + '%'; div.style.width = Math.max(0.05, r - l) + '%';
      minimapEl.appendChild(div);
      continue;
    }}
    const isUtr = f.type === '5utr' || f.type === '3utr';
    const isCdsDomain = f.type === 'cds' && f.plotGroup === 'CDS_domain';
    const div = document.createElement('div');
    div.className = 'mini-feat ' + (isUtr ? 'utr' : (isCdsDomain ? 'cds-domain' : 'cds'));
    const l = posToFullPct(f.start), r = posToFullPct(f.end);
    div.style.left = l + '%'; div.style.width = Math.max(0.05, r - l) + '%';
    minimapEl.appendChild(div);
  }}
  // Viewport rectangle.
  const viewportDiv = document.createElement('div');
  viewportDiv.className = 'viewport';
  viewportDiv.id = 'viewport-rect';
  const [vLo, vHi] = viewport;
  viewportDiv.style.left  = vToFullPct(vLo) + '%';
  viewportDiv.style.width = (vToFullPct(vHi) - vToFullPct(vLo)) + '%';
  const handleL = document.createElement('div');
  handleL.className = 'viewport-handle left';
  const handleR = document.createElement('div');
  handleR.className = 'viewport-handle right';
  viewportDiv.appendChild(handleL); viewportDiv.appendChild(handleR);
  minimapEl.appendChild(viewportDiv);

  attachMinimapDrag(viewportDiv, 'pan');
  attachMinimapDrag(handleL,     'left');
  attachMinimapDrag(handleR,     'right');
}}

function updateReadout() {{
  const [vLo, vHi] = viewport;
  const a = Math.round(vToGenomic(vLo)), b = Math.round(vToGenomic(vHi));
  readoutEl.textContent = fmt(a) + '–' + fmt(b)
    + '  (' + fmt(b - a + 1) + ' bp visible)';
}}

/* ============================================================
 * Interaction: minimap drag (pan / resize) + main plot wheel/dbl-click
 * ============================================================ */
let dragState = null;

function attachMinimapDrag(el, kind) {{
  el.addEventListener('mousedown', (e) => {{
    const rect = minimapEl.getBoundingClientRect();
    dragState = {{
      kind, startX: e.clientX, mmW: rect.width, mmLeft: rect.left,
      startVLo: viewport[0], startVHi: viewport[1],
    }};
    el.classList.add('dragging');
    e.preventDefault(); e.stopPropagation();
  }});
}}

minimapEl.addEventListener('click', (e) => {{
  if (e.target.closest('.viewport')) return;  // ignore clicks on rect / handles
  const rect = minimapEl.getBoundingClientRect();
  const rel = (e.clientX - rect.left) / rect.width;
  const v = rel * axis.totalVirtual;
  const w = viewport[1] - viewport[0];
  let a = v - w / 2, b = v + w / 2;
  if (a < 0)                    {{ a = 0; b = w; }}
  if (b > axis.totalVirtual)    {{ b = axis.totalVirtual; a = b - w; }}
  viewport = [a, b];
  renderMain(); renderMinimap(); renderTicks(); updateReadout();
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
  renderMain(); renderMinimap(); renderTicks(); updateReadout();
}});
window.addEventListener('mouseup', () => {{
  if (dragState) {{
    document.querySelectorAll('.dragging').forEach(el => el.classList.remove('dragging'));
    dragState = null;
  }}
}});

// Mouse-wheel zoom on main plot, centered on cursor.
mainPlot.addEventListener('wheel', (e) => {{
  e.preventDefault();
  const rect = mainPlot.getBoundingClientRect();
  const relX = (e.clientX - rect.left) / rect.width;
  const [a, b] = viewport;
  const w = b - a;
  const center = a + relX * w;
  // Zoom in on positive wheel-up, out on wheel-down.
  const factor = e.deltaY < 0 ? 0.8 : 1.25;
  let newW = Math.max(20, Math.min(axis.totalVirtual, w * factor));
  let newA = center - relX * newW;
  let newB = newA + newW;
  if (newA < 0)                  {{ newA = 0; newB = newW; }}
  if (newB > axis.totalVirtual)  {{ newB = axis.totalVirtual; newA = newB - newW; }}
  viewport = [newA, newB];
  renderMain(); renderMinimap(); renderTicks(); updateReadout();
}}, {{ passive: false }});

// Box-zoom (default): drag a rectangle on the main plot. The unselected
// area is shaded; releasing zooms the viewport to the selection. Mirrors
// the vCRE-vis-js Dashboard.svelte behaviour (onPlotMouseDown/Up).
// Shift-drag on the main plot pans instead.
let boxState = null;  // {{ startPx, endPx, mpW, mpLeft }}
let panState = null;
let didDrag  = false;
let shadeL, shadeR, boxRect, boxReadout;

function ensureBoxNodes() {{
  // Re-create if missing or detached. renderMain() purges all children of
  // mainPlot except the baseline, which would otherwise leave our refs
  // dangling after a wheel-zoom or layout change.
  if (shadeL && shadeL.isConnected && shadeR && shadeR.isConnected
      && boxRect && boxRect.isConnected && boxReadout && boxReadout.isConnected) return;
  if (shadeL && shadeL.parentNode)    shadeL.parentNode.removeChild(shadeL);
  if (shadeR && shadeR.parentNode)    shadeR.parentNode.removeChild(shadeR);
  if (boxRect && boxRect.parentNode)  boxRect.parentNode.removeChild(boxRect);
  if (boxReadout && boxReadout.parentNode) boxReadout.parentNode.removeChild(boxReadout);
  shadeL     = document.createElement('div'); shadeL.className     = 'box-shade';
  shadeR     = document.createElement('div'); shadeR.className     = 'box-shade';
  boxRect    = document.createElement('div'); boxRect.className    = 'box-rect';
  boxReadout = document.createElement('div'); boxReadout.className = 'box-readout';
  shadeL.style.display = shadeR.style.display = 'none';
  boxRect.style.display = boxReadout.style.display = 'none';
  mainPlot.appendChild(shadeL);
  mainPlot.appendChild(shadeR);
  mainPlot.appendChild(boxRect);
  mainPlot.appendChild(boxReadout);
}}

function drawBoxOverlay(px0, px1, mpW) {{
  ensureBoxNodes();
  const lo = Math.max(0, Math.min(px0, px1));
  const hi = Math.min(mpW, Math.max(px0, px1));
  shadeL.style.left  = '0px';
  shadeL.style.width = lo + 'px';
  shadeL.style.display = 'block';
  shadeR.style.left  = hi + 'px';
  shadeR.style.width = Math.max(0, mpW - hi) + 'px';
  shadeR.style.display = 'block';
  boxRect.style.left  = lo + 'px';
  boxRect.style.width = Math.max(0, hi - lo) + 'px';
  boxRect.style.display = 'block';
  // Show the genomic range of the selection as a small badge above the rect.
  const [a, b] = viewport;
  const w = b - a;
  const bpLo = Math.round(vToGenomic(a + (lo / mpW) * w));
  const bpHi = Math.round(vToGenomic(a + (hi / mpW) * w));
  boxReadout.textContent = fmt(bpLo) + '–' + fmt(bpHi)
    + '  (' + fmt(Math.max(0, bpHi - bpLo + 1)) + ' bp)';
  // Anchor the badge to the center of the rect, clamped to the plot edges.
  // Width is unknown until rendered, so position via transform: translateX(-50%).
  const cx = Math.max(0, Math.min(mpW, (lo + hi) / 2));
  boxReadout.style.left = cx + 'px';
  boxReadout.style.transform = 'translateX(-50%)';
  boxReadout.style.display = 'block';
}}

function clearBoxOverlay() {{
  if (shadeL)     shadeL.style.display     = 'none';
  if (shadeR)     shadeR.style.display     = 'none';
  if (boxRect)    boxRect.style.display    = 'none';
  if (boxReadout) boxReadout.style.display = 'none';
}}

mainPlot.addEventListener('mousedown', (e) => {{
  if (e.button !== 0) return;
  // Accept any descendant of mainPlot as a drag origin — feature boxes
  // cover most of the plot, so requiring empty space made the gesture
  // hard to discover.
  const rect = mainPlot.getBoundingClientRect();
  const x = e.clientX - rect.left;
  if (e.shiftKey) {{
    panState = {{
      startX: e.clientX, mpW: rect.width,
      startVLo: viewport[0], startVHi: viewport[1],
    }};
    mainPlot.classList.add('dragging');
  }} else {{
    boxState = {{ startPx: x, endPx: x, mpW: rect.width, mpLeft: rect.left }};
    didDrag = false;
    // Show a zero-width preview at the click position so the user sees
    // *something* even before they've moved the mouse.
    drawBoxOverlay(x, x, rect.width);
  }}
  e.preventDefault();
}});

window.addEventListener('mousemove', (e) => {{
  if (boxState) {{
    const x = Math.max(0, Math.min(boxState.mpW, e.clientX - boxState.mpLeft));
    boxState.endPx = x;
    if (Math.abs(x - boxState.startPx) >= 1) didDrag = true;
    // Re-render every move so the rect tracks the cursor smoothly.
    drawBoxOverlay(boxState.startPx, boxState.endPx, boxState.mpW);
    return;
  }}
  if (!panState) return;
  const dx = e.clientX - panState.startX;
  const w = panState.startVHi - panState.startVLo;
  const dV = -(dx / panState.mpW) * w;
  let a = Math.max(0, Math.min(axis.totalVirtual - w, panState.startVLo + dV));
  viewport = [a, a + w];
  renderMain(); renderMinimap(); renderTicks(); updateReadout();
}});

window.addEventListener('mouseup', () => {{
  if (boxState) {{
    const {{ startPx, endPx, mpW }} = boxState;
    boxState = null;
    clearBoxOverlay();
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
          renderMain(); renderMinimap(); renderTicks(); updateReadout();
        }}
      }}
    }}
    didDrag = false;
  }}
  if (panState) {{ mainPlot.classList.remove('dragging'); panState = null; }}
}});

mainPlot.addEventListener('dblclick', () => {{
  viewport = [0, axis.totalVirtual];
  renderMain(); renderMinimap(); renderTicks(); updateReadout();
}});

document.getElementById('btn-reset').addEventListener('click', () => {{
  viewport = [0, axis.totalVirtual];
  renderMain(); renderMinimap(); renderTicks(); updateReadout();
}});

document.getElementById('btn-fit-domain').addEventListener('click', () => {{
  if (!STRUCTURE.domains || !STRUCTURE.domains.length) return;
  const d = STRUCTURE.domains[0];
  if (d.genomicStart == null || d.genomicEnd == null) return;
  const a = axis.posToV(d.genomicStart);
  const b = axis.posToV(d.genomicEnd);
  // Add 10 % padding on each side.
  const pad = (b - a) * 0.10;
  let vLo = Math.max(0, a - pad);
  let vHi = Math.min(axis.totalVirtual, b + pad);
  if (vHi - vLo < 20) {{ vLo = Math.max(0, vLo - 10); vHi = vLo + 40; }}
  viewport = [vLo, vHi];
  renderMain(); renderMinimap(); renderTicks(); updateReadout();
}});

/* ============================================================
 * Boot
 * ============================================================ */
renderHeader(STRUCTURE);
renderStrandArrow(STRUCTURE.strand);
render();
document.querySelectorAll('input[name=mode]').forEach(r =>
  r.addEventListener('change', render));
document.getElementById('show-utr').addEventListener('change', render);

/* ============================================================
 * Auto-resize handshake (only fires when embedded in an iframe).
 * The parent listens for postMessage and resizes the iframe to fit.
 * ============================================================ */
(function () {{
  if (window.parent === window) return;  // standalone page, nothing to do.
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
  // Re-broadcast when controls change (e.g., layout toggle adds/removes ticks).
  document.querySelectorAll('input').forEach(el =>
    el.addEventListener('change', () => setTimeout(broadcast, 30)));
  // Also catch any reflow we missed (e.g., font load, viewer redraw).
  if (typeof ResizeObserver !== 'undefined') {{
    new ResizeObserver(broadcast).observe(document.body);
  }}
}})();
</script>
</body>
</html>
"""


def _render_to_string(segs, *,
                      link_template: str | None = None,
                      plot_height: int = 80) -> str:
    """Render the interactive viewer to an HTML string (no file I/O)."""
    payload = _segments_to_payload(segs)
    if not payload:
        raise ValueError("no segments to plot")
    s0 = segs[0]
    title_bits = [getattr(s0, "gene_name", "") or "",
                  getattr(s0, "protein_id", "") or "",
                  getattr(s0, "transcript_id", "") or ""]
    title = " · ".join(b for b in title_bits if b) or "fastCDS"

    # A user template wins; otherwise auto-link from the ID (Ensembl / RefSeq /
    # UniProt). Custom-GTF IDs resolve to None, so no link is shown.
    from .plot import _resolve_link   # lazy: avoids an import cycle
    url = _resolve_link(link_template, s0)
    if url:
        payload["externalLink"] = url

    return HTML_TEMPLATE.format(
        title_safe=html.escape(title),
        payload_json=json.dumps(payload),
        domain_colors_json=json.dumps(DOMAIN_COLOR),
        compact_intron_bp=COMPACT_INTRON_BP,
        plot_height=int(plot_height),
    )


def _coerce_segments(source, input_id: str | None = None) -> list:
    """Accept either a ready ``Segment`` list (used as-is) or a high-level
    source — a :class:`MappingResult`, an isoform DataFrame, or a path to an
    ``isoform_structure.tsv`` — and return one isoform's ``Segment`` list.

    With several isoforms, pass ``input_id`` to pick one; a lone isoform is
    chosen automatically. Keeps callers off the private
    ``_segments_from_dataframe`` helper.
    """
    if isinstance(source, list):
        return source
    from .plot import _segments_from_source   # lazy: avoids an import cycle
    by_id = _segments_from_source(source)
    if input_id is not None:
        try:
            return by_id[input_id]
        except KeyError:
            raise KeyError(
                f"input_id {input_id!r} not found; available: {sorted(by_id)}"
            ) from None
    if len(by_id) == 1:
        return next(iter(by_id.values()))
    raise ValueError(
        f"source has {len(by_id)} isoforms {sorted(by_id)}; "
        "pass input_id= to choose one."
    )


def render_interactive_html(source, out: str, *,
                         input_id: str | None = None,
                         link_template: str | None = None,
                         plot_height: int = 80) -> None:
    """Write a standalone HTML viewer to `out`.

    Vanilla JS + inline CSS — no CDN, no external deps. Includes a vCRE-style
    minimap with a draggable viewport rectangle, drag-to-box-zoom on the main
    plot, and wheel-zoom around the cursor.

    `source` is either a ``Segment`` list or a high-level source (a
    :class:`MappingResult`, an isoform DataFrame, or a path to
    ``isoform_structure.tsv``); pass ``input_id`` to pick one isoform when the
    source holds several.
    """
    segs = _coerce_segments(source, input_id)
    html_out = _render_to_string(
        segs, link_template=link_template, plot_height=plot_height,
    )
    with open(out, "w", encoding="utf-8") as f:
        f.write(html_out)


_JUPYTER_IFRAME_SEQ = 0


def render_interactive_jupyter(source, *,
                            input_id: str | None = None,
                            height: int | None = None,
                            plot_height: int = 140,
                            link_template: str | None = None):
    """Embed the viewer in a Jupyter notebook.

    Returns an ``IPython.display.HTML`` value that the notebook frontend
    renders inline. The viewer runs inside a sandboxed iframe (via
    ``srcdoc``) so its CSS and global JS don't leak into the host
    notebook. By default the iframe auto-resizes to the viewer's natural
    content height — no inner scrollbar — via a postMessage handshake
    between the iframe and a small listener injected next to it.

    Parameters
    ----------
    source : list[Segment] or MappingResult or DataFrame or path
        Either a ready ``Segment`` list, or a high-level source (a
        :class:`MappingResult`, an isoform DataFrame, or a path to an
        ``isoform_structure.tsv``). No need to touch internal helpers.
    input_id : str, optional
        Which isoform to draw when ``source`` holds several; a lone
        isoform is picked automatically.
    height : int, optional
        Pin the iframe height in pixels. By default the iframe sizes
        itself to fit the viewer exactly. Pass an int to override (e.g.,
        for nbconvert/static exports where postMessage may not fire).
    plot_height : int
        Main-track height in pixels (default 140). Higher values make the
        CDS/UTR boxes more prominent — useful in a notebook where the
        viewer competes for attention with other cells.
    link_template : str, optional
        Override the header linkout, which is otherwise auto-derived from the
        ID (Ensembl / RefSeq / UniProt; none for custom IDs). A URL template
        with placeholders ``{protein_id}``, ``{gene_name}``,
        ``{transcript_id}``, ``{chrom}``, ``{start}``, ``{end}``.
    """
    from IPython.display import HTML  # local import — only needed in notebooks
    segs = _coerce_segments(source, input_id)
    html_str = _render_to_string(
        segs, link_template=link_template, plot_height=plot_height,
    )
    # srcdoc requires HTML-attribute escaping (quotes + ampersands).
    safe = (html_str
            .replace("&", "&amp;")
            .replace('"', "&quot;"))

    global _JUPYTER_IFRAME_SEQ
    _JUPYTER_IFRAME_SEQ += 1
    iframe_id = f"fastcds-viewer-{_JUPYTER_IFRAME_SEQ}"
    # Generous fallback height so the panel never shows half-rendered if the
    # postMessage handshake is blocked (nbconvert static exports, some
    # JupyterLab CSP configs). The JS will shrink/grow as needed.
    initial_h = int(height) if height is not None else 560

    iframe = (
        f'<iframe id="{iframe_id}" srcdoc="{safe}" '
        f'style="width: 100%; height: {initial_h}px; border: none; '
        f'border-radius: 6px; display: block;" '
        f'sandbox="allow-scripts"></iframe>'
    )
    # The wrapper div is here for two reasons: (1) it dodges IPython's
    # "Consider using IFrame instead" warning (which only fires when the
    # payload starts with <iframe …>); (2) the inline <script> next to the
    # iframe listens for the postMessage and resizes the iframe to its
    # content height, so the viewer fits without an inner scrollbar.
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
