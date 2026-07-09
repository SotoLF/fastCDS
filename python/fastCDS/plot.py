#!/usr/bin/env python3
"""Render `isoform_structure.tsv` into a transcript-architecture plot.

Usage examples
--------------
# Single domain to a static figure (format follows the extension)
fastCDS plot --isoform isoform_structure.tsv --input-id RD1 --out RD1.pdf

# Every query in the TSV to a multipage PDF
fastCDS plot --isoform isoform_structure.tsv --all --out queries.pdf

# Interactive HTML viewer (self-contained vanilla JS by default)
fastCDS plot --isoform isoform_structure.tsv --input-id RD1 --out RD1.html

# ...or the plotly engine for the same .html output
fastCDS plot --isoform isoform_structure.tsv --input-id RD1 \
                    --out RD1.html --engine plotly

One output target, one rule: the ``--out`` extension picks the format.
``.pdf`` / ``.png`` / ``.svg`` render a static matplotlib figure; ``.html``
renders the interactive viewer, whose engine (``js`` or ``plotly``) is chosen
with ``--engine``. This mirrors the paper's three output types -- the BED12
track (written by ``fastCDS map``), the static figure, and the interactive
viewer -- where ``js`` and ``plotly`` are two engines of the same viewer, not
separate output types.

Design notes
------------
- The C++ binary writes a tidy "one row per structural segment" table; this
  script just groups by `input_id` and draws boxes. We never re-derive
  coordinates from the genome here.
- Strand is honored: genomic coordinates are always on the X axis (so the plot
  matches IGV), but a small arrow at the top indicates the 5'→3' direction of
  translation.
- Color scheme mirrors `plot_group`:
    five_prime_UTR        -> light grey
    three_prime_UTR       -> light grey
    CDS / CDS_no_domain   -> blue
    CDS_domain            -> red
    intron / intron_domain_span
                          -> thin line, dashed if outside domain span
- `--no-introns` collapses the X axis to a concatenated CDS-only view (the
  "spliced transcript" style) using `feature_order_transcript` ordering.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from collections import defaultdict
from dataclasses import dataclass
from typing import Iterable

# matplotlib is required; plotly is optional and only imported when the
# interactive viewer is requested with --engine plotly.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages


# --------------------------------------------------------------------------- #
# Row parsing
# --------------------------------------------------------------------------- #

# A single row of isoform_structure.tsv. We only keep the columns the plotter
# actually uses. NA-able numerics become `None`.
@dataclass
class Segment:
    input_id: str
    gene_name: str
    transcript_id: str
    protein_id: str
    domain_id: str
    chrom: str
    strand: str
    feature_type: str
    feature_id: str
    feature_part: int
    feature_genomic_start: int
    feature_genomic_end: int
    feature_order_transcript: int
    overlaps_domain: str
    plot_group: str
    is_mane_select: str = "NA"
    is_ensembl_canonical: str = "NA"


def _i_or_none(x: str) -> int | None:
    if x == "NA" or x == "":
        return None
    try:
        return int(x)
    except ValueError:
        return None


def load_isoform_tsv(path: str) -> dict[str, list[Segment]]:
    """Group rows by input_id, preserving file order within each group."""
    by_id: dict[str, list[Segment]] = defaultdict(list)
    with open(path, "r", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {
            "input_id", "feature_type", "feature_genomic_start",
            "feature_genomic_end", "plot_group", "strand", "chrom",
        }
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise SystemExit(
                f"isoform tsv {path!r} is missing columns: {sorted(missing)}"
            )
        for row in reader:
            seg = Segment(
                input_id=row.get("input_id", ""),
                gene_name=row.get("gene_name", "") or "",
                transcript_id=row.get("transcript_id", "") or "",
                protein_id=row.get("protein_id", "") or "",
                domain_id=row.get("domain_id", "") or "",
                chrom=row.get("chrom", "") or "",
                strand=row.get("strand", "+") or "+",
                feature_type=row.get("feature_type", "") or "",
                feature_id=row.get("feature_id", "") or "",
                feature_part=_i_or_none(row.get("feature_part", "1")) or 1,
                feature_genomic_start=int(row["feature_genomic_start"]),
                feature_genomic_end=int(row["feature_genomic_end"]),
                feature_order_transcript=(
                    _i_or_none(row.get("feature_order_transcript", "0")) or 0
                ),
                overlaps_domain=row.get("overlaps_domain", "NA") or "NA",
                plot_group=row.get("plot_group", "") or "",
                is_mane_select=row.get("is_mane_select", "NA") or "NA",
                is_ensembl_canonical=row.get("is_ensembl_canonical", "NA") or "NA",
            )
            by_id[seg.input_id].append(seg)

    # `isoform_structure.tsv` is the plotter's input, but the sibling
    # `domain_cds_segments.tsv` / `domain_introns.tsv` share its exact header,
    # so they load without a column error and would draw a partial view. A real
    # transcript structure always has CDS rows; if there are none, the file is
    # almost certainly an introns-only / non-coding output, not an isoform TSV.
    types = {s.feature_type for segs in by_id.values() for s in segs}
    if by_id and "CDS" not in types:
        print(
            f"note: {path!r} has no CDS rows (feature types {sorted(types)}); "
            "this does not look like an isoform_structure.tsv. `fastCDS plot` "
            "renders isoform_structure.tsv (from `--output isoform` or the "
            "default `--output all`); the coding/introns/span outputs are "
            "genome-browser tracks, not plot input.",
            file=sys.stderr,
        )
    return by_id


# --------------------------------------------------------------------------- #
# Colors / heights
# --------------------------------------------------------------------------- #

# Palette. UTRs are deliberately tan/warm so they read clearly against the
# grey intron line — keeping everything in greyscale (the old default) made
# UTR and intron visually merge.
PLOT_GROUP_COLORS = {
    "five_prime_UTR":     "#F0C078",   # warm tan
    "three_prime_UTR":    "#F0C078",
    "CDS":                "#1F77B4",
    "CDS_no_domain":      "#1F77B4",
    "CDS_domain":         "#D62728",
    "intron":             "#666666",
    "intron_domain_span": "#D62728",
}

# Height (in y units) for each feature kind. UTRs are slightly thinner than
# CDS, introns are a thin line.
def feature_height(feature_type: str) -> float:
    if feature_type in ("CDS",):
        return 0.6
    if feature_type in ("five_prime_UTR", "three_prime_UTR"):
        return 0.4
    return 0.02  # intron line


# --------------------------------------------------------------------------- #
# Drawing
# --------------------------------------------------------------------------- #

def _title_for(segs: list[Segment]) -> str:
    s = segs[0]
    bits = []
    if s.gene_name: bits.append(s.gene_name)
    if s.protein_id: bits.append(s.protein_id)
    if s.transcript_id: bits.append(s.transcript_id)
    if s.domain_id: bits.append(f"domain={s.domain_id}")
    flags = []
    if s.is_mane_select == "true": flags.append("MANE_Select")
    if s.is_ensembl_canonical == "true": flags.append("Ensembl_canonical")
    if flags: bits.append("[" + ",".join(flags) + "]")
    bits.append(f"({s.chrom} {s.strand})")
    return "  ".join(bits)


def _draw_genomic(ax, segs: list[Segment], *, show_introns: bool,
                  show_utr: bool, highlight_domain: bool) -> None:
    """Draw on genomic-coordinate X axis (matches IGV)."""
    y_center = 0.5
    for s in segs:
        if s.feature_type == "intron" and not show_introns:
            continue
        if s.feature_type in ("five_prime_UTR", "three_prime_UTR") and not show_utr:
            continue
        color = PLOT_GROUP_COLORS.get(s.plot_group, "#777777")
        if not highlight_domain:
            # Collapse the domain coloring to the plain feature color.
            if s.plot_group == "CDS_domain":
                color = PLOT_GROUP_COLORS["CDS"]
            elif s.plot_group == "intron_domain_span":
                color = PLOT_GROUP_COLORS["intron"]
        h = feature_height(s.feature_type)
        x = s.feature_genomic_start - 0.5  # center on integer base
        w = s.feature_genomic_end - s.feature_genomic_start + 1
        if s.feature_type == "intron":
            # Single thin horizontal line spanning the intron.
            ax.plot([x, x + w], [y_center, y_center], color=color, lw=1.2,
                    solid_capstyle="butt", zorder=1)
            continue
        rect = mpatches.Rectangle(
            (x, y_center - h / 2.0), w, h,
            facecolor=color, edgecolor="black", linewidth=0.4, zorder=2,
        )
        ax.add_patch(rect)

    # Strand arrow.
    xmin = min(s.feature_genomic_start for s in segs)
    xmax = max(s.feature_genomic_end for s in segs)
    arrow_y = y_center + 0.55
    pad = (xmax - xmin) * 0.04 if xmax > xmin else 1
    if (segs[0].strand or "+") == "+":
        ax.annotate("", xy=(xmax - pad, arrow_y), xytext=(xmin + pad, arrow_y),
                    arrowprops=dict(arrowstyle="->", color="#444444", lw=1))
        ax.text(xmin, arrow_y + 0.05, "5'", color="#444444", va="bottom")
        ax.text(xmax, arrow_y + 0.05, "3'", color="#444444", va="bottom",
                ha="right")
    else:
        ax.annotate("", xy=(xmin + pad, arrow_y), xytext=(xmax - pad, arrow_y),
                    arrowprops=dict(arrowstyle="->", color="#444444", lw=1))
        ax.text(xmax, arrow_y + 0.05, "5'", color="#444444", va="bottom",
                ha="right")
        ax.text(xmin, arrow_y + 0.05, "3'", color="#444444", va="bottom")

    ax.set_xlim(xmin - (xmax - xmin) * 0.02 if xmax > xmin else xmin - 1,
                xmax + (xmax - xmin) * 0.02 if xmax > xmin else xmax + 1)
    ax.set_ylim(0, 1.6)
    ax.set_yticks([])
    ax.set_xlabel(f"genomic position ({segs[0].chrom})")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)


def _draw_compact_genomic(ax, segs: list[Segment], *, show_utr: bool,
                          highlight_domain: bool) -> None:
    """Genomic-style layout with introns clamped to a fixed visual width.

    Use this when one long intron dwarfs the rest of the transcript — full
    genomic coords would render every exon as a hair-thin sliver. We keep
    CDS / UTR at their true bp scale (so domain proportions stay honest) and
    replace each intron's width with `intron_disp_w`. A small "//" marker on
    the intron line flags compression, and the x-axis ticks carry the true
    genomic coordinate at each segment boundary.
    """
    # Genomic-order draw; sentinel introns get the clamped width.
    ordered = sorted(segs, key=lambda s: s.feature_genomic_start)
    if not show_utr:
        ordered = [s for s in ordered if s.feature_type not in
                   ("five_prime_UTR", "three_prime_UTR")]

    # Display width of an intron — fixed fraction of the largest non-intron
    # feature so the figure stays balanced.
    nonintron_widths = [
        s.feature_genomic_end - s.feature_genomic_start + 1
        for s in ordered if s.feature_type != "intron"
    ]
    if not nonintron_widths:
        return
    intron_disp_w = max(max(nonintron_widths) * 0.15, 50.0)

    y_center = 0.5
    cursor = 0.0
    # display_x_of_boundary[i] = display cursor at the *start* of ordered[i].
    boundaries = []
    for s in ordered:
        boundaries.append((cursor, s))
        true_w = s.feature_genomic_end - s.feature_genomic_start + 1
        disp_w = intron_disp_w if s.feature_type == "intron" else true_w
        color = PLOT_GROUP_COLORS.get(s.plot_group, "#777777")
        if not highlight_domain:
            if s.plot_group == "CDS_domain":
                color = PLOT_GROUP_COLORS["CDS"]
            elif s.plot_group == "intron_domain_span":
                color = PLOT_GROUP_COLORS["intron"]
        if s.feature_type == "intron":
            ax.plot([cursor, cursor + disp_w], [y_center, y_center],
                    color=color, lw=1.2, solid_capstyle="butt", zorder=1)
            # "//" compression marker — only when we actually compressed.
            if disp_w < true_w:
                ax.text(cursor + disp_w / 2.0, y_center + 0.04, "//",
                        ha="center", va="bottom", fontsize=8, color="#444444",
                        zorder=3)
        else:
            h = feature_height(s.feature_type)
            rect = mpatches.Rectangle(
                (cursor, y_center - h / 2.0), disp_w, h,
                facecolor=color, edgecolor="black", linewidth=0.4, zorder=2,
            )
            ax.add_patch(rect)
        cursor += disp_w

    # Strand arrow at the top.
    xmin, xmax = 0.0, cursor
    arrow_y = y_center + 0.55
    pad = (xmax - xmin) * 0.04 if xmax > xmin else 1
    if (ordered[0].strand or "+") == "+":
        ax.annotate("", xy=(xmax - pad, arrow_y), xytext=(xmin + pad, arrow_y),
                    arrowprops=dict(arrowstyle="->", color="#444444", lw=1))
        ax.text(xmin, arrow_y + 0.05, "5'", color="#444444", va="bottom")
        ax.text(xmax, arrow_y + 0.05, "3'", color="#444444", va="bottom",
                ha="right")
    else:
        ax.annotate("", xy=(xmin + pad, arrow_y), xytext=(xmax - pad, arrow_y),
                    arrowprops=dict(arrowstyle="->", color="#444444", lw=1))
        ax.text(xmax, arrow_y + 0.05, "5'", color="#444444", va="bottom",
                ha="right")
        ax.text(xmin, arrow_y + 0.05, "3'", color="#444444", va="bottom")

    # X-axis ticks at segment boundaries, labelled with the true genomic
    # coord (5' end of each feature in display order; plus the 3' end of the
    # last feature so the right edge is anchored).
    ticks, labels = [], []
    for (disp_x, seg) in boundaries:
        ticks.append(disp_x)
        labels.append(str(seg.feature_genomic_start))
    ticks.append(cursor)
    labels.append(str(ordered[-1].feature_genomic_end))
    # Thin out if too many ticks (cap at ~10).
    if len(ticks) > 10:
        step = max(1, len(ticks) // 10)
        ticks = ticks[::step]
        labels = labels[::step]
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels, fontsize=7, rotation=30, ha="right")
    ax.set_xlim(xmin - (xmax - xmin) * 0.02, xmax + (xmax - xmin) * 0.02)
    ax.set_ylim(0, 1.6)
    ax.set_yticks([])
    ax.set_xlabel(f"genomic position, introns compressed ({ordered[0].chrom})")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)


def _draw_spliced(ax, segs: list[Segment], *, show_utr: bool,
                  highlight_domain: bool) -> None:
    """Concatenate non-intron features in translation order. Useful when
    introns dwarf the rest of the figure."""
    kept = [s for s in segs if s.feature_type != "intron"]
    if not show_utr:
        kept = [s for s in kept if s.feature_type not in
                ("five_prime_UTR", "three_prime_UTR")]
    kept = sorted(kept, key=lambda s: s.feature_order_transcript)

    y_center = 0.5
    cursor = 0.0
    for s in kept:
        color = PLOT_GROUP_COLORS.get(s.plot_group, "#777777")
        if not highlight_domain and s.plot_group == "CDS_domain":
            color = PLOT_GROUP_COLORS["CDS"]
        w = s.feature_genomic_end - s.feature_genomic_start + 1
        h = feature_height(s.feature_type)
        rect = mpatches.Rectangle(
            (cursor, y_center - h / 2.0), w, h,
            facecolor=color, edgecolor="black", linewidth=0.4,
        )
        ax.add_patch(rect)
        cursor += w
    ax.set_xlim(0, cursor)
    ax.set_ylim(0, 1.6)
    ax.set_yticks([])
    ax.set_xlabel("spliced position (nt, translation order)")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)


def _legend(fig, ax, *, show_introns: bool, show_utr: bool,
            highlight_domain: bool) -> None:
    """Legend at the bottom of the *figure* (not the axes). With
    constrained_layout enabled, matplotlib reserves space for it so it never
    collides with the xlabel."""
    from matplotlib.lines import Line2D
    handles = []
    if show_utr:
        handles.append(mpatches.Patch(facecolor=PLOT_GROUP_COLORS["five_prime_UTR"],
                                      edgecolor="black", label="UTR"))
    handles.append(mpatches.Patch(facecolor=PLOT_GROUP_COLORS["CDS"],
                                  edgecolor="black", label="CDS"))
    if highlight_domain:
        handles.append(mpatches.Patch(facecolor=PLOT_GROUP_COLORS["CDS_domain"],
                                      edgecolor="black", label="CDS (domain)"))
    if show_introns:
        # Render intron as a thin line in the legend to match how it's drawn.
        handles.append(Line2D([0], [0], color=PLOT_GROUP_COLORS["intron"],
                              lw=1.5, label="intron"))
        if highlight_domain:
            handles.append(Line2D([0], [0],
                                  color=PLOT_GROUP_COLORS["intron_domain_span"],
                                  lw=1.5, label="intron (domain span)"))
    # `loc="outside lower center"` is a constrained_layout-aware placement —
    # matplotlib reserves a row at the bottom of the figure for the legend so
    # it never collides with the xlabel.
    try:
        fig.legend(handles=handles, loc="outside lower center",
                   ncol=len(handles), frameon=False)
    except (ValueError, TypeError):
        # Older matplotlib (<3.6): fall back to a manually anchored legend.
        fig.legend(handles=handles, loc="lower center",
                   ncol=len(handles), frameon=False,
                   bbox_to_anchor=(0.5, 0.0))


def render_one(segs: list[Segment], *, out: str | None = None,
               width: float = 12.0, height: float = 2.6,
               title: str | None = None,
               show_introns: bool = True,
               show_utr: bool = True,
               highlight_domain: bool = True,
               spliced: bool = False,
               compact_genomic: bool = False,
               fig=None):
    """Render one query. If `fig` is provided we draw into it (PdfPages).

    Uses constrained_layout so the bottom legend, the xlabel and the title
    don't compete for space. constrained_layout supersedes tight_layout and
    handles fig.legend() positioning correctly out of the box.

    Layout modes (mutually exclusive — first true one wins):
        spliced=True          -> introns dropped, exons concatenated
        compact_genomic=True  -> genomic order, introns clamped to a fixed
                                 display width (best for long-intron genes)
        otherwise             -> full genomic layout (1 bp = 1 unit)
    """
    own_fig = fig is None
    if own_fig:
        fig, ax = plt.subplots(figsize=(width, height), constrained_layout=True)
    else:
        # The caller (e.g. PdfPages) already created the figure. Make sure it
        # uses constrained_layout so the same logic applies.
        try:
            fig.set_layout_engine("constrained")
        except Exception:
            pass
        ax = fig.add_subplot(1, 1, 1)
    if spliced:
        _draw_spliced(ax, segs, show_utr=show_utr,
                      highlight_domain=highlight_domain)
    elif compact_genomic:
        _draw_compact_genomic(ax, segs, show_utr=show_utr,
                              highlight_domain=highlight_domain)
    else:
        _draw_genomic(ax, segs, show_introns=show_introns, show_utr=show_utr,
                      highlight_domain=highlight_domain)
    _legend(fig, ax,
            show_introns=show_introns and not spliced,
            show_utr=show_utr,
            highlight_domain=highlight_domain)
    ax.set_title(title or _title_for(segs))
    if own_fig and out:
        fig.savefig(out)
        plt.close(fig)
        print(f"Wrote {out}", file=sys.stderr)
    return fig


# --------------------------------------------------------------------------- #
# Optional plotly HTML export
# --------------------------------------------------------------------------- #

from ._interactive_html import render_interactive_html


def _expand_link_template(template: str, s0: Segment) -> str | None:
    """Fill `{protein_id}`, `{gene_name}`, `{transcript_id}`, `{chrom}`,
    `{start}`, `{end}` placeholders from the first plottable segment. Returns
    None if any referenced placeholder is empty (so we don't generate a dead
    link to TFRegDB2 etc.)."""
    fields = {
        "protein_id":   s0.protein_id,
        "gene_name":    s0.gene_name,
        "transcript_id": s0.transcript_id,
        "chrom":        s0.chrom,
        "start":        s0.feature_genomic_start,
        "end":          s0.feature_genomic_end,
    }
    try:
        url = template.format(**fields)
    except (KeyError, IndexError):
        return None
    # If any used field was empty, the template will contain `''` substrings
    # that defeat the purpose of a linkout. Heuristic: bail if the URL still
    # looks broken (`//`, trailing `/`, etc. after host).
    if "/{" in url or "}" in url:
        return None
    return url


# UniProt accession grammar (the official 6/10-char pattern).
_UNIPROT_RE = re.compile(
    r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$"
)


def _auto_link(s0: Segment) -> str | None:
    """Build a sensible external linkout from the ID form alone, so a plain
    `fastCDS plot` gets a working link with no `--link-template`.

    Recognises Ensembl stable IDs (ENSP/ENST/ENSG and their species-prefixed
    forms), NCBI RefSeq (NP_/XP_/YP_ proteins, NM_/XM_/NR_/XR_ transcripts),
    and UniProt accessions. Custom GTF IDs match nothing, so no link is
    emitted (returns None) rather than a dead one.
    """
    pid = (s0.protein_id or "").split(".")[0]
    tid = (s0.transcript_id or "").split(".")[0]
    # Ensembl: the species-agnostic /id/ resolver redirects a protein ID to its
    # protein page and a transcript ID to its transcript page, so we never need
    # to know the species. Prefer the protein ID, fall back to the transcript.
    for x in (pid, tid):
        if re.fullmatch(r"ENS[A-Z]*[GTP][0-9]+", x):
            return f"https://www.ensembl.org/id/{x}"
    if re.fullmatch(r"[NXY]P_[0-9]+", pid):          # RefSeq protein
        return f"https://www.ncbi.nlm.nih.gov/protein/{pid}"
    if re.fullmatch(r"[NX][MR]_[0-9]+", tid):        # RefSeq transcript / mRNA
        return f"https://www.ncbi.nlm.nih.gov/nuccore/{tid}"
    if _UNIPROT_RE.match(pid):                        # UniProt accession
        return f"https://www.uniprot.org/uniprotkb/{pid}/entry"
    return None


def _resolve_link(link_template: str | None, s0: Segment) -> str | None:
    """A user `--link-template` wins; otherwise fall back to the automatic
    ID-based linkout. Returns a ready URL, or None for no link."""
    if link_template:
        return _expand_link_template(link_template, s0)
    return _auto_link(s0)


def render_html(segs: list[Segment], out: str, *, highlight_domain: bool = True,
                show_introns: bool = True, show_utr: bool = True,
                link_template: str | None = None) -> None:
    """Render an interactive HTML view with a bottom rangeslider.

    The output uses ``go.Bar`` (horizontal) for UTR/CDS segments and
    ``go.Scattergl`` for introns, with a rangeslider attached to the x-axis.
    Drag the rangeslider handles at the bottom to zoom / pan the top view in
    real time. Hover on any segment to see its feature_id, type, coordinates,
    strand, and (for CDS in a domain) the aa range it encodes.

    Why this differs from the matplotlib path: shapes (``add_shape``) don't
    respond to hover and won't move with the rangeslider, so we draw the
    same gene model with ``go.Bar`` traces — one per visual category — so
    each gets its own legend entry and its own hover tooltip.
    """
    try:
        import plotly.graph_objects as go
    except ImportError:
        raise SystemExit(
            "plotly is required for --engine plotly. Install it with "
            "'pip install plotly', or use the default --engine js viewer "
            "(no extra dependencies)."
        )

    # ---- filter & group ------------------------------------------------ #
    def keep(s: Segment) -> bool:
        if s.feature_type == "intron" and not show_introns:
            return False
        if s.feature_type in ("five_prime_UTR", "three_prime_UTR") and not show_utr:
            return False
        return True

    segs_keep = [s for s in segs if keep(s)]
    if not segs_keep:
        raise ValueError("no segments left to plot after filtering")

    # Bucket segments by the legend category they belong to. We collapse
    # five_prime_UTR / three_prime_UTR into a single "UTR" legend entry (they
    # share a color); CDS_no_domain and any bare "CDS" go into "CDS". When
    # --no-highlight is in effect we also pull CDS_domain → CDS and
    # intron_domain_span → intron so the highlighting disappears entirely.
    CATEGORIES = [
        ("UTR",                  "#F0C078", 0.4),
        ("CDS",                  PLOT_GROUP_COLORS["CDS"],         0.6),
        ("CDS (domain)",         PLOT_GROUP_COLORS["CDS_domain"],  0.6),
        ("intron",               PLOT_GROUP_COLORS["intron"],      None),  # line
        ("intron (domain span)", PLOT_GROUP_COLORS["intron_domain_span"], None),
    ]
    buckets: dict[str, list[Segment]] = {cat[0]: [] for cat in CATEGORIES}

    for s in segs_keep:
        pg = s.plot_group
        ft = s.feature_type
        if not highlight_domain:
            if pg == "CDS_domain": pg = "CDS_no_domain"
            if pg == "intron_domain_span": pg = "intron"
        if ft in ("five_prime_UTR", "three_prime_UTR"):
            buckets["UTR"].append(s)
        elif ft == "intron":
            buckets["intron (domain span)" if pg == "intron_domain_span"
                    else "intron"].append(s)
        elif ft == "CDS":
            buckets["CDS (domain)" if pg == "CDS_domain" else "CDS"].append(s)

    # ---- build the traces --------------------------------------------- #
    fig = go.Figure()

    for name, color, height in CATEGORIES:
        group = buckets.get(name) or []
        if not group:
            continue

        if height is None:
            # Intron — one Scattergl trace with `None` breaks between segments.
            xs: list[float | None] = []
            ys: list[float | None] = []
            for s in group:
                xs += [s.feature_genomic_start - 0.5,
                       s.feature_genomic_end + 0.5, None]
                ys += [0, 0, None]
            fig.add_trace(go.Scattergl(
                x=xs, y=ys,
                mode="lines",
                line=dict(color=color, width=1.5),
                name=name,
                hoverinfo="skip",
            ))
            continue

        # UTR / CDS — one Bar trace with a row per segment.
        widths_x = [s.feature_genomic_end - s.feature_genomic_start + 1
                    for s in group]
        bases    = [s.feature_genomic_start - 0.5 for s in group]
        ys       = [0] * len(group)

        # Hover columns. aa coords are NA on UTR rows; we encode "NA" so the
        # tooltip stays readable.
        def _aa(v):
            try:
                vi = int(v) if v not in (None, "", "NA") else None
            except (TypeError, ValueError):
                vi = None
            return "" if vi is None else str(vi)

        custom = [
            [s.feature_id, s.feature_type,
             s.chrom, s.feature_genomic_start, s.feature_genomic_end,
             s.strand or "+",
             s.plot_group]
            for s in group
        ]
        fig.add_trace(go.Bar(
            x=widths_x, y=ys, base=bases,
            orientation="h",
            width=height,
            marker=dict(color=color, line=dict(color="black", width=0.5)),
            customdata=custom,
            hovertemplate=(
                "<b>%{customdata[0]}</b> · %{customdata[1]}<br>"
                "%{customdata[2]}:%{customdata[3]:,}–%{customdata[4]:,}  "
                "(%{customdata[5]} strand)<br>"
                "plot_group: %{customdata[6]}"
                "<extra></extra>"
            ),
            name=name,
        ))

    # ---- layout, rangeslider, strand arrow ---------------------------- #
    s0 = segs_keep[0]
    strand = s0.strand or "+"

    xmin = min(s.feature_genomic_start for s in segs_keep)
    xmax = max(s.feature_genomic_end for s in segs_keep)
    pad  = max(1, int((xmax - xmin) * 0.02))

    # Strand arrow annotation. xref/ayref="x"/"x" so the arrow scales with
    # the rangeslider zoom — it always points the full visible span.
    if strand == "+":
        arrow_from, arrow_to = "x domain", "x domain"
        ax, x = 0.02, 0.98
        label_xref, label_x = "paper", 0.0
    else:
        arrow_from, arrow_to = "x domain", "x domain"
        ax, x = 0.98, 0.02
        label_xref, label_x = "paper", 1.0

    # Title — same fields the matplotlib path shows.
    title_bits = []
    if s0.gene_name:     title_bits.append(s0.gene_name)
    if s0.protein_id:    title_bits.append(s0.protein_id)
    if s0.transcript_id: title_bits.append(s0.transcript_id)
    if s0.domain_id:     title_bits.append(f"domain={s0.domain_id}")
    flags = []
    if s0.is_mane_select == "true":       flags.append("MANE_Select")
    if s0.is_ensembl_canonical == "true": flags.append("Ensembl_canonical")
    if flags: title_bits.append("[" + ",".join(flags) + "]")
    title_bits.append(f"({s0.chrom} {strand})")
    # Optional clickable linkout (Ensembl / RefSeq / UniProt, auto by ID).
    # Rendered as an `<a href>` inside the plotly title — plotly title text
    # accepts a small set of HTML tags including `<a>`.
    url = _resolve_link(link_template, s0)
    if url:
        title_bits.append(
            f'<a href="{url}" target="_blank" '
            f'style="color:#1f77b4;text-decoration:underline">↗ link</a>'
        )
    title = "   ".join(title_bits)

    fig.update_layout(
        title=dict(text=title, x=0.01, xanchor="left",
                   font=dict(size=13)),
        barmode="overlay",
        plot_bgcolor="white",
        paper_bgcolor="white",
        height=420,
        margin=dict(l=60, r=40, t=70, b=70),
        legend=dict(orientation="h", yanchor="bottom", y=1.02,
                    xanchor="right", x=1.0, font=dict(size=11)),
        annotations=[
            dict(
                xref=arrow_from, yref="y",
                axref=arrow_to,  ayref="y",
                x=x, ax=ax, y=0.55, ay=0.55,
                showarrow=True,
                arrowhead=2, arrowsize=1, arrowwidth=1.2, arrowcolor="#555",
            ),
            dict(
                xref="paper", yref="paper",
                x=0.0 if strand == "+" else 1.0, y=1.04,
                text=("5'" if strand == "+" else "3'"),
                showarrow=False, font=dict(size=11, color="#444"),
            ),
            dict(
                xref="paper", yref="paper",
                x=1.0 if strand == "+" else 0.0, y=1.04,
                text=("3'" if strand == "+" else "5'"),
                showarrow=False, font=dict(size=11, color="#444"),
            ),
        ],
        xaxis=dict(
            title=f"genomic position ({s0.chrom})",
            range=[xmin - pad, xmax + pad],
            # The rangeslider is the "slider in the bottom where you can zoom,
            # extend the boxes, etc." — drag its handles to update the main
            # view in real time. Plotly renders a mini-thumbnail of the bars
            # inside it, which doubles as an overview map.
            #
            # Pin its range to the same window as the main axis (autorange off),
            # otherwise plotly re-derives a slightly wider extent for the slider
            # and the minimap thumbnail no longer lines up under the main plot.
            rangeslider=dict(
                visible=True,
                autorange=False,
                range=[xmin - pad, xmax + pad],
                thickness=0.10,
                bgcolor="#f5f5f5",
                bordercolor="#ccc",
                borderwidth=1,
            ),
            showline=True, linecolor="black", linewidth=1, mirror=False,
            tickformat=",",   # genomic coords get thousands separators
        ),
        yaxis=dict(
            visible=False,
            range=[-0.5, 0.8],
            fixedrange=True,
        ),
    )

    fig.write_html(out, include_plotlyjs="cdn")
    print(f"Wrote {out}", file=sys.stderr)


# --------------------------------------------------------------------------- #
# High-level Python API
# --------------------------------------------------------------------------- #

def _segments_from_dataframe(df) -> dict[str, list[Segment]]:
    """Adapt a pandas isoform DataFrame into the internal Segment dicts.

    The plotter's drawing routines were written against a list of small
    dataclass Segments because that was the original CSV path. Reusing them
    means the matplotlib code stays untouched no matter how the rows arrived.
    """
    by_id: dict[str, list[Segment]] = defaultdict(list)
    # We tolerate missing optional columns (older outputs without MANE flags).
    def _get(row, key, default=""):
        v = row.get(key, default) if hasattr(row, "get") else default
        if v is None: return default
        # Avoid pandas NaN sneaking into string slots.
        try:
            import math
            if isinstance(v, float) and math.isnan(v): return default
        except Exception:
            pass
        return v

    for r in df.to_dict(orient="records"):
        try:
            seg = Segment(
                input_id=str(_get(r, "input_id")),
                gene_name=str(_get(r, "gene_name")),
                transcript_id=str(_get(r, "transcript_id")),
                protein_id=str(_get(r, "protein_id")),
                domain_id=str(_get(r, "domain_id")),
                chrom=str(_get(r, "chrom")),
                strand=str(_get(r, "strand") or "+"),
                feature_type=str(_get(r, "feature_type")),
                feature_id=str(_get(r, "feature_id")),
                feature_part=int(_get(r, "feature_part", 1) or 1),
                feature_genomic_start=int(r["feature_genomic_start"]),
                feature_genomic_end=int(r["feature_genomic_end"]),
                feature_order_transcript=int(_get(r, "feature_order_transcript", 0) or 0),
                overlaps_domain=str(_get(r, "overlaps_domain", "NA")),
                plot_group=str(_get(r, "plot_group")),
                is_mane_select=str(_get(r, "is_mane_select", "NA")),
                is_ensembl_canonical=str(_get(r, "is_ensembl_canonical", "NA")),
            )
            by_id[seg.input_id].append(seg)
        except (KeyError, ValueError, TypeError):
            # Skip rows with NaN coordinates etc. — they're produced by the
            # mapper for unmapped queries and aren't drawable anyway.
            continue
    return by_id


def _segments_from_source(source) -> dict[str, list[Segment]]:
    """Accept any of: MappingResult, isoform DataFrame, path-to-tsv,
    or a pre-built dict[input_id -> list[Segment]]."""
    # MappingResult — duck-typed to avoid an import cycle.
    if hasattr(source, "isoform") and not isinstance(source, dict):
        return _segments_from_dataframe(source.isoform)
    if isinstance(source, dict):
        return source
    try:
        import pandas as pd
        if isinstance(source, pd.DataFrame):
            return _segments_from_dataframe(source)
    except ImportError:
        pass
    if isinstance(source, (str, os.PathLike)):
        return load_isoform_tsv(os.fspath(source))
    raise TypeError(
        f"Unsupported source type for plot(): {type(source).__name__}. "
        f"Expected MappingResult, pandas DataFrame, or path to "
        f"isoform_structure.tsv."
    )


# Static image extensions handled by matplotlib. Everything else is routed by
# `_out_kind`; only ``.html`` (the interactive viewer) is also supported.
_STATIC_EXTS = (".pdf", ".png", ".svg")
_ENGINES = ("js", "plotly")


def _out_kind(out: str) -> str:
    """Classify an output path by extension: ``"static"`` or ``"html"``.

    Raises ValueError (with the offending value and the valid set) for any
    other extension, so callers fail fast with an actionable message.
    """
    ext = os.path.splitext(out)[1].lower()
    if ext in _STATIC_EXTS:
        return "static"
    if ext == ".html":
        return "html"
    raise ValueError(
        f"unsupported output extension {ext or '(none)'!r} for out={out!r}. "
        f"Use a static image ({', '.join(_STATIC_EXTS)}) or an interactive "
        f"viewer (.html)."
    )


def _render_one_target(segs, out, *, engine, width, height, title,
                       show_introns, show_utr, highlight_domain,
                       spliced, compact_genomic, link_template, fig=None):
    """Render one query's segments to a single ``out`` file, dispatching on the
    extension. Returns the matplotlib Figure for static output, else None.

    ``engine``/``spliced``/``compact_genomic``/``width``/``height`` only affect
    the path they belong to (engine -> .html; the rest -> static); they are
    ignored for the other format.
    """
    if _out_kind(out) == "static":
        return render_one(segs, out=out, width=width, height=height,
                          title=title, show_introns=show_introns,
                          show_utr=show_utr, highlight_domain=highlight_domain,
                          spliced=spliced, compact_genomic=compact_genomic,
                          fig=fig)
    # .html interactive viewer.
    if engine not in _ENGINES:
        raise ValueError(
            f"unknown engine {engine!r}; choose one of {_ENGINES}.")
    if engine == "plotly":
        render_html(segs, out, highlight_domain=highlight_domain,
                    show_introns=show_introns, show_utr=show_utr,
                    link_template=link_template)
    else:  # "js" — self-contained, no CDN, no plotly dependency.
        render_interactive_html(segs, out, link_template=link_template)
    return None


def plot(source, *, input_id: str | None = None,
         out: str | None = None,
         engine: str = "js",
         title: str | None = None,
         width: float = 12.0,
         height: float = 2.6,
         show_introns: bool = True,
         show_utr: bool = True,
         highlight_domain: bool = True,
         spliced: bool = False,
         compact_genomic: bool = False,
         link_template: str | None = None):
    """Render one query (or every query) from a mapping result.

    The ``out`` extension picks the format: ``.pdf`` / ``.png`` / ``.svg``
    render a static matplotlib figure, and ``.html`` renders the interactive
    viewer. For ``.html``, ``engine`` selects the renderer -- ``"js"`` (the
    default, a self-contained vanilla-JS viewer with no dependencies) or
    ``"plotly"`` (needs ``pip install plotly``).

    Parameters
    ----------
    source : MappingResult | pandas.DataFrame | str
        Either a :class:`MappingResult` (we use its ``.isoform`` table), a
        DataFrame in the same schema, or the path to an
        ``isoform_structure.tsv``.
    input_id : str | None
        Which query to draw. If ``None`` and ``source`` carries a single
        query, that one is used; otherwise pass ``input_id`` explicitly.
        Use :func:`plot_all` for batches.
    out : str | None
        Output file. Extension drives the format (see above). If ``None``,
        a static figure is built and returned without being saved.
    engine : {"js", "plotly"}
        Interactive renderer for ``.html`` output. Ignored for static output.
    title : str | None
        Override the auto-generated title.
    width, height : float
        Static figure size in inches (matplotlib path only).
    show_introns, show_utr, highlight_domain, spliced, compact_genomic : bool
        Same toggles as the CLI flags.

    Returns
    -------
    matplotlib.figure.Figure | None
        The Figure for static output (or when ``out`` is ``None``); ``None``
        for ``.html`` output.

    Examples
    --------
    >>> import fastCDS as fc
    >>> r = fc.map_query("ENSP00000269305", 10, 50, "AD1", index="human.idx")
    >>> fc.plot(r, out="AD1.pdf")                  # static
    >>> fc.plot(r, out="AD1.html")                 # interactive (vanilla JS)
    >>> fc.plot(r, out="AD1.html", engine="plotly")
    >>> # From a DataFrame:
    >>> fc.plot(r.isoform, input_id="AD1", out="AD1.pdf")
    """
    by_id = _segments_from_source(source)
    if not by_id:
        raise ValueError("source has no plottable rows.")
    if input_id is None:
        ids = list(by_id.keys())
        if len(ids) != 1:
            raise ValueError(
                f"source has {len(ids)} input_ids - pass input_id=... "
                f"or call plot_all(). Available: {ids[:5]}{'...' if len(ids)>5 else ''}")
        input_id = ids[0]
    if input_id not in by_id:
        keys = ", ".join(sorted(by_id)[:10])
        raise KeyError(f"input_id {input_id!r} not in source. Sample: {keys}")

    if out is None:
        # No file requested — hand back the static figure for programmatic use.
        return render_one(by_id[input_id], out=None, width=width,
                          height=height, title=title,
                          show_introns=show_introns, show_utr=show_utr,
                          highlight_domain=highlight_domain, spliced=spliced,
                          compact_genomic=compact_genomic)
    return _render_one_target(
        by_id[input_id], out, engine=engine, width=width, height=height,
        title=title, show_introns=show_introns, show_utr=show_utr,
        highlight_domain=highlight_domain, spliced=spliced,
        compact_genomic=compact_genomic, link_template=link_template)


def plot_all(source, *, out: str,
             engine: str = "js",
             width: float = 12.0,
             height: float = 2.6,
             show_introns: bool = True,
             show_utr: bool = True,
             highlight_domain: bool = True,
             spliced: bool = False,
             compact_genomic: bool = False,
             link_template: str | None = None) -> int:
    """Render every input_id in ``source``. The ``out`` extension drives the
    format exactly like :func:`plot`. A ``.pdf`` target becomes a single
    multipage PDF (one query per page); any other target (``.png`` / ``.svg``
    / ``.html``) writes one file per query named ``base.<input_id>.ext``.

    Returns the number of pages / files written.
    """
    by_id = _segments_from_source(source)
    ids = list(by_id.keys())
    if not ids:
        raise ValueError("source has no plottable rows.")

    kind = _out_kind(out)  # validate up front

    if kind == "static" and out.lower().endswith(".pdf"):
        with PdfPages(out) as pdf:
            for i in ids:
                fig = plt.figure(figsize=(width, height),
                                 constrained_layout=True)
                render_one(by_id[i], fig=fig, width=width, height=height,
                           show_introns=show_introns, show_utr=show_utr,
                           highlight_domain=highlight_domain, spliced=spliced,
                           compact_genomic=compact_genomic)
                pdf.savefig(fig)
                plt.close(fig)
        print(f"Wrote {out} ({len(ids)} pages)", file=sys.stderr)
    else:
        base, ext = os.path.splitext(out)
        for i in ids:
            _render_one_target(
                by_id[i], f"{base}.{i}{ext}", engine=engine,
                width=width, height=height, title=None,
                show_introns=show_introns, show_utr=show_utr,
                highlight_domain=highlight_domain, spliced=spliced,
                compact_genomic=compact_genomic, link_template=link_template)
    return len(ids)


# Compatibility alias matching the user-facing schematic name.
plot_isoform = plot


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #

def _argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="fastCDS plot",
        description="Render isoform_structure.tsv. The --out extension picks "
                    "the format: .pdf/.png/.svg (static matplotlib) or .html "
                    "(interactive viewer; engine via --engine).",
    )
    p.add_argument("--isoform", required=True,
                   help="Path to isoform_structure.tsv (the plot-ready table).")
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--input-id", help="A single input_id to render.")
    g.add_argument("--all", action="store_true",
                   help="Render every input_id (multipage PDF if --out ends in .pdf).")
    p.add_argument("--out", required=True,
                   help="Output file. Format follows the extension: "
                        ".pdf/.png/.svg -> static figure, .html -> interactive "
                        "viewer.")
    p.add_argument("--engine", choices=_ENGINES, default="js",
                   help="Interactive engine for .html output: 'js' (default, "
                        "self-contained) or 'plotly'. Both ship with fastCDS. "
                        "Ignored for static output.")
    p.add_argument("--title", help="Override the figure title.")
    p.add_argument("--width", type=float, default=12.0, help="Figure width in inches.")
    p.add_argument("--height", type=float, default=2.6, help="Figure height in inches.")
    p.add_argument("--no-highlight", dest="highlight", action="store_false",
                   help="Do not color CDS_domain segments differently.")
    p.add_argument("--no-introns", dest="introns", action="store_false",
                   help="Hide intron lines.")
    p.add_argument("--no-utr", dest="utr", action="store_false",
                   help="Hide UTR boxes.")
    p.add_argument("--spliced", action="store_true",
                   help="Concatenate non-intron features in translation order.")
    p.add_argument("--compact-genomic", action="store_true",
                   help="Keep genomic order, but clamp every intron to a "
                        "fixed display width (best for long-intron genes).")
    p.add_argument("--link-template", default=None,
                   help="Override the external linkout shown next to the title "
                        "in .html output. By default the link is chosen "
                        "automatically from the ID (Ensembl / NCBI RefSeq / "
                        "UniProt; custom-GTF IDs get no link). Pass a URL "
                        "template to override, with placeholders {protein_id}, "
                        "{gene_name}, {transcript_id}, {chrom}, {start}, {end}. "
                        "Example: 'https://www.ensembl.org/Homo_sapiens/"
                        "Transcript/ProteinSummary?p={protein_id}'")
    p.set_defaults(highlight=True, introns=True, utr=True)
    return p


def main(argv: list[str] | None = None) -> int:
    args = _argparser().parse_args(argv)
    if args.spliced and args.compact_genomic:
        print("error: --spliced and --compact-genomic are mutually exclusive",
              file=sys.stderr)
        return 2
    by_id = load_isoform_tsv(args.isoform)
    if not by_id:
        print(f"No rows in {args.isoform}", file=sys.stderr)
        return 1

    if args.input_id:
        if args.input_id not in by_id:
            keys = ", ".join(sorted(by_id)[:10])
            print(f"input_id {args.input_id!r} not found. Sample: {keys}",
                  file=sys.stderr)
            return 2
        ids = [args.input_id]
    else:
        ids = list(by_id.keys())

    # Validate the output extension once, up front, with a clear message.
    try:
        kind = _out_kind(args.out)
    except ValueError as e:
        print(f"error: {e}", file=sys.stderr)
        return 2

    if kind == "static" and args.out.lower().endswith(".pdf") and len(ids) > 1:
        # Multipage PDF — one query per page.
        with PdfPages(args.out) as pdf:
            for i in ids:
                fig = plt.figure(figsize=(args.width, args.height),
                                 constrained_layout=True)
                render_one(by_id[i], fig=fig,
                           width=args.width, height=args.height,
                           show_introns=args.introns, show_utr=args.utr,
                           highlight_domain=args.highlight,
                           spliced=args.spliced,
                           compact_genomic=args.compact_genomic)
                pdf.savefig(fig)
                plt.close(fig)
        print(f"Wrote {args.out} ({len(ids)} pages)", file=sys.stderr)
    elif len(ids) > 1:
        # Any other target holds one query per file: base.<input_id>.ext.
        base, ext = os.path.splitext(args.out)
        for i in ids:
            _render_one_target(
                by_id[i], f"{base}.{i}{ext}", engine=args.engine,
                width=args.width, height=args.height, title=None,
                show_introns=args.introns, show_utr=args.utr,
                highlight_domain=args.highlight, spliced=args.spliced,
                compact_genomic=args.compact_genomic,
                link_template=args.link_template)
    else:
        _render_one_target(
            by_id[ids[0]], args.out, engine=args.engine,
            width=args.width, height=args.height, title=args.title,
            show_introns=args.introns, show_utr=args.utr,
            highlight_domain=args.highlight, spliced=args.spliced,
            compact_genomic=args.compact_genomic,
            link_template=args.link_template)
    return 0


if __name__ == "__main__":
    sys.exit(main())
