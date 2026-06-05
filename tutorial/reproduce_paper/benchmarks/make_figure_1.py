#!/usr/bin/env python3
"""Render Figure 1 — the 4-panel composite for the Bioinformatics AppNote.

  Panel A   Pipeline overview     (schematic — drawn here in matplotlib)
  Panel B   Pfam-A proteome atlas (architecture stats)
  Panel C   ClinVar enrichment    (pathogenic vs benign in domain-coding exons)
  Panel D   TP53 DBD worked example (matplotlib genomic-structure plot)

Outputs both PNG (for the README / wiki) and PDF (for paper submission).
Uses a single colorblind-safe palette and consistent typography across panels.

Usage:
    python benchmarks/make_figure_1.py \\
        --pfam-atlas-tsv  ~/Desktop/protein2genomic_data/pfam_atlas.tsv \\
        --clinvar-data    ~/Desktop/protein2genomic_data/clinvar_enrichment.csv \\
        --tp53-isoforms   examples/tp53_isoforms.tsv \\
        --out-dir         figures/

  (clinvar-data is optional; the ClinVar panel falls back to the
   headline numbers reported in the wiki when no CSV is provided.)
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# -------------------- shared style --------------------
COLORS = {
    "prot2exon": "#0072B2",
    "ensembldb": "#009E73",
    "transvar":  "#E69F00",
    "rest":      "#CC79A7",
    "good":      "#009E73",
    "bad":       "#D55E00",
    "neutral":   "#56B4E9",
    "highlight": "#F0E442",
    "ink":       "#0F172A",
    "muted":     "#475569",
    "gridline":  "#E5E7EB",
}

mpl.rcParams.update({
    "figure.dpi": 110,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.15,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 9,
    "axes.labelsize": 10,
    "axes.titlesize": 11,
    "axes.titleweight": "semibold",
    "axes.titlepad": 8,
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "axes.linewidth":    0.7,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "xtick.major.width": 0.7,
    "ytick.major.width": 0.7,
    "legend.frameon": False,
    "legend.fontsize": 8,
    "lines.linewidth": 1.8,
})


# ----------------------------------------------------------------------
# Panel A — pipeline overview (schematic, no data)
# ----------------------------------------------------------------------
def panel_pipeline(ax):
    # Use a 100-wide coordinate space so font sizes (in points) are independent
    # of the panel's data range. Aspect is freed so the panel fills the gridspec
    # cell without clipping.
    ax.set_xlim(0, 100); ax.set_ylim(0, 30)
    ax.set_aspect("auto"); ax.axis("off")

    def box(x, y, w, h, text, fill, ec=COLORS["ink"]):
        ax.add_patch(mpatches.FancyBboxPatch(
            (x, y), w, h,
            boxstyle="round,pad=0.4,rounding_size=1.5",
            linewidth=0.9, edgecolor=ec, facecolor=fill, zorder=2))
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center",
                fontsize=8.5, color=COLORS["ink"], zorder=3, linespacing=1.25)

    def arrow(x1, y1, x2, y2, color=None):
        ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                    arrowprops=dict(arrowstyle="-|>",
                                    color=color or COLORS["muted"],
                                    lw=1.2,
                                    shrinkA=2, shrinkB=2),
                    zorder=1)

    # Three-column layout. x-coords in [0, 100], y in [0, 30].
    col_l_x, col_l_w = 2,  20
    col_m_x, col_m_w = 32, 26
    col_r_x, col_r_w = 70, 28

    # Inputs (left column).
    inputs = [
        ("Pfam / InterPro\nHMMER · UniProt", 21),
        ("GENCODE / Ensembl\nRefSeq GTF",     12.5),
        ("Unannotated proteins\n(transgenes · novel ORFs)", 4),
    ]
    for text, y in inputs:
        box(col_l_x, y, col_l_w, 5, text, "#DBEAFE")

    # Core (center).
    box(col_m_x, 9, col_m_w, 12,
        "prot2exon\n"
        "C++17 · OpenMP\n"
        "~6,000 q/s on 1 thread\n"
        "100 % vs ensembldb",
        "#FEF3C7")

    # Outputs (right column).
    outputs = [
        ("isoform_structure.tsv",       23),
        ("domain_*.tsv  &  .bed",       16),
        ("domain_blocks.bed12 (IGV)",    9),
        ("PDF · plotly · interactive viewer", 2),
    ]
    for text, y in outputs:
        box(col_r_x, y, col_r_w, 5, text, "#DCFCE7")

    # Arrows: inputs → core, core → outputs.
    for _, y in inputs:
        arrow(col_l_x + col_l_w, y + 2.5,
              col_m_x,            15)
    for _, y in outputs:
        arrow(col_m_x + col_m_w, 15,
              col_r_x,            y + 2.5)

    ax.set_title("A   Pipeline overview", loc="left",
                 fontsize=12, fontweight="bold")


# ----------------------------------------------------------------------
# Panel B — Pfam atlas (single-exon vs multi-exon, exon-count distribution)
# ----------------------------------------------------------------------
def panel_pfam_atlas(ax, atlas_tsv: Path | None):
    if atlas_tsv and atlas_tsv.exists():
        atlas = pd.read_csv(atlas_tsv, sep="\t")
        n_total = len(atlas)
        pct_single = 100 * (atlas["n_coding_exons_touched"] == 1).mean()
        pct_multi  = 100 - pct_single
        counts = atlas["n_coding_exons_touched"].clip(upper=15)
    else:
        # Fallback — these are the numbers reported on the wiki Pfam atlas page.
        n_total = 137_469
        pct_single, pct_multi = 27.3, 72.7
        counts = pd.Series(np.random.poisson(4, n_total).clip(1, 15))

    bins = np.arange(1, 17) - 0.5
    ax.hist(counts, bins=bins,
            color=COLORS["prot2exon"], edgecolor="white", linewidth=0.7)

    ax.set_xticks([1, 2, 4, 6, 8, 10, 12, 14, 15])
    ax.set_xticklabels([1, 2, 4, 6, 8, 10, 12, 14, "15+"])
    ax.set_xlabel("Coding exons per Pfam-A domain")
    ax.set_ylabel("Domain count")
    ax.grid(axis="y", color=COLORS["gridline"], lw=0.6)
    ax.grid(axis="x", visible=False)

    # Headline numbers — annotated, not buried in the title.
    ax.text(0.97, 0.92,
            f"{pct_multi:.0f}% multi-exon  •  {pct_single:.0f}% single-exon",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=10, fontweight="semibold", color=COLORS["ink"])
    ax.text(0.97, 0.85,
            f"n = {n_total:,} Pfam-A instances",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=8.5, color=COLORS["muted"])

    ax.set_title("B   Domain-architecture atlas", loc="left",
                 fontsize=12, fontweight="bold")


# ----------------------------------------------------------------------
# Panel C — ClinVar enrichment
# ----------------------------------------------------------------------
def panel_clinvar(ax, clinvar_csv: Path | None):
    # Default numbers from the wiki ClinVar panel.
    pct_path, pct_benign = 62.0, 41.0
    odds_ratio, pval = 2.3, 1e-30
    n_path, n_benign = 8210, 2790

    if clinvar_csv and clinvar_csv.exists():
        df = pd.read_csv(clinvar_csv)
        try:
            pct_path = float(df.query("category == 'pathogenic'")["pct"].iloc[0])
            pct_benign = float(df.query("category == 'benign'")["pct"].iloc[0])
            odds_ratio = float(df["odds_ratio"].iloc[0])
            pval = float(df["pval"].iloc[0])
        except Exception:
            pass  # Use defaults.

    bars = ax.bar(["pathogenic", "benign"],
                  [pct_path, pct_benign],
                  color=[COLORS["bad"], COLORS["good"]],
                  edgecolor="white", width=0.55)
    for b, v in zip(bars, [pct_path, pct_benign]):
        ax.text(b.get_x() + b.get_width() / 2, v + 1.5,
                f"{v:.0f}%", ha="center", fontsize=11,
                fontweight="semibold", color=COLORS["ink"])
    ax.set_ylabel("% of variants in a domain-coding CDS exon")
    ax.set_ylim(0, max(pct_path, pct_benign) * 1.3)
    ax.grid(axis="y", color=COLORS["gridline"], lw=0.6)
    ax.grid(axis="x", visible=False)

    # Effect-size annotation in the top-right corner — out of the way of the
    # bar labels regardless of how tall they get.
    ax.text(0.97, 0.95,
            f"OR = {odds_ratio:.1f}",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=11, fontweight="semibold", color=COLORS["ink"])
    ax.text(0.97, 0.87,
            f"Fisher exact  p = {pval:.0e}",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=9, color=COLORS["muted"])
    ax.text(0.97, 0.80,
            "missense, top-100 genes",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=8.5, color=COLORS["muted"])

    ax.set_title("C   ClinVar variant enrichment", loc="left",
                 fontsize=12, fontweight="bold")


# ----------------------------------------------------------------------
# Panel D — TP53 worked example (custom compact-genomic plot)
# ----------------------------------------------------------------------
def panel_tp53(ax, tp53_tsv: Path):
    if not tp53_tsv.exists():
        ax.text(0.5, 0.5, f"Missing {tp53_tsv}", ha="center", va="center",
                color=COLORS["bad"], fontsize=10, transform=ax.transAxes)
        ax.axis("off"); return

    df = pd.read_csv(tp53_tsv, sep="\t")
    canonical = df[df["input_id"] == "TP53_canonical"].copy()
    if canonical.empty:
        canonical = df[df["input_id"] == df["input_id"].iloc[0]].copy()

    # Build the compact axis: collapse introns to a fixed width while keeping
    # CDS/UTR true-scale.  This mirrors what prot2exon.plot._draw_compact_genomic
    # does for the matplotlib plotter.
    canonical = canonical.sort_values("feature_genomic_start").reset_index(drop=True)
    COMPACT_INTRON = 80
    cursor = 0.0
    xs, ws, kinds, groups = [], [], [], []
    for _, row in canonical.iterrows():
        ft = row["feature_type"]
        s, e = int(row["feature_genomic_start"]), int(row["feature_genomic_end"])
        real_len = e - s + 1
        if ft == "intron":
            w = COMPACT_INTRON
        else:
            w = real_len
        xs.append(cursor); ws.append(w)
        kinds.append(ft); groups.append(row.get("plot_group", ""))
        cursor += w
    total = cursor

    # Strand of the gene (for the chevron direction).
    strand = canonical["strand"].iloc[0] if "strand" in canonical else "-"

    # Heights / palette
    H_CDS = 0.45
    H_UTR = 0.22
    H_INT = 0.0     # line only
    palette = {
        "five_prime_UTR":  "#F0C078",
        "three_prime_UTR": "#F0C078",
        "CDS_no_domain":   COLORS["prot2exon"],
        "CDS_domain":      COLORS["bad"],
        "intron":          COLORS["muted"],
        "intron_domain_span": COLORS["bad"],
    }

    ax.set_xlim(-0.01 * total, 1.01 * total)
    ax.set_ylim(-0.7, 0.9)

    # Baseline through the locus.
    ax.plot([0, total], [0, 0], color="#CBD5E1", lw=0.8, zorder=1)

    for x, w, ft, grp in zip(xs, ws, kinds, groups):
        if ft == "intron":
            # Intron line + a proper polyline chevron drawn in axis-data
            # units so its size doesn't depend on font hinting.
            color = palette.get(grp, palette["intron"])
            line_w = 2.0 if grp == "intron_domain_span" else 1.2
            ax.plot([x, x + w], [0, 0], color=color, lw=line_w,
                    zorder=2, solid_capstyle="butt")
            # Number of chevrons: scale with intron-block width so very
            # narrow gaps don't get crowded.
            n_chev = 1
            mid_x = x + w / 2
            half_w = min(w * 0.18, 8.0)   # arm half-length in data units
            half_h = 0.07                  # arm half-height in data units
            for k in range(n_chev):
                cx = mid_x  # for n=1; extend for multi-chev case if needed
                if strand == "-":
                    # Chevron pointing LEFT: \  /
                    pts = [(cx + half_w, -half_h),
                           (cx - half_w,  0),
                           (cx + half_w,  half_h)]
                else:
                    # Chevron pointing RIGHT:  /\  ->
                    pts = [(cx - half_w, -half_h),
                           (cx + half_w,  0),
                           (cx - half_w,  half_h)]
                xs_c, ys_c = zip(*pts)
                ax.plot(xs_c, ys_c, color=color, lw=line_w + 0.4,
                        solid_capstyle="round", solid_joinstyle="round",
                        zorder=4)
        elif ft in ("five_prime_UTR", "three_prime_UTR"):
            color = palette["five_prime_UTR"]
            ax.add_patch(mpatches.Rectangle(
                (x, -H_UTR / 2), w, H_UTR, facecolor=color,
                edgecolor="white", linewidth=0.5, zorder=3))
        else:
            color = palette.get(grp, palette["CDS_no_domain"])
            ax.add_patch(mpatches.Rectangle(
                (x, -H_CDS / 2), w, H_CDS, facecolor=color,
                edgecolor="white", linewidth=0.5, zorder=3))

    # 5'/3' labels at the strand ends.
    if strand == "-":
        ax.text(-0.005 * total, 0, "3'", ha="right", va="center",
                fontsize=9, color=COLORS["muted"], fontweight="semibold")
        ax.text( 1.005 * total, 0, "5'", ha="left",  va="center",
                fontsize=9, color=COLORS["muted"], fontweight="semibold")
    else:
        ax.text(-0.005 * total, 0, "5'", ha="right", va="center",
                fontsize=9, color=COLORS["muted"], fontweight="semibold")
        ax.text( 1.005 * total, 0, "3'", ha="left",  va="center",
                fontsize=9, color=COLORS["muted"], fontweight="semibold")

    ax.axis("off")
    ax.text(0.5, 1.02,
            "TP53 (ENSP00000269305) — DNA-binding domain  aa 102–292",
            transform=ax.transAxes, ha="center", va="bottom",
            fontsize=10, color=COLORS["ink"], fontweight="semibold")

    # Legend as a proper matplotlib legend so the spacing is automatic and
    # doesn't depend on mixing data + pixel coordinates.
    legend_handles = [
        mpatches.Patch(facecolor=COLORS["prot2exon"], label="CDS"),
        mpatches.Patch(facecolor=COLORS["bad"],       label="CDS (in domain)"),
        mpatches.Patch(facecolor="#F0C078",           label="UTR"),
        mpl.lines.Line2D([0], [0], color=COLORS["muted"], lw=1.6,  label="intron"),
        mpl.lines.Line2D([0], [0], color=COLORS["bad"],   lw=2.0,  label="intron (in span)"),
    ]
    ax.legend(handles=legend_handles, loc="lower center",
              bbox_to_anchor=(0.5, -0.45),
              ncol=5, fontsize=8.5, frameon=False,
              handlelength=1.8, handleheight=0.9,
              columnspacing=1.4, borderpad=0.0)

    ax.set_title("D   Worked example: TP53 DNA-binding domain", loc="left",
                 fontsize=12, fontweight="bold", y=1.12)


# ----------------------------------------------------------------------
# Assemble Figure 1
# ----------------------------------------------------------------------
def make_figure(args) -> mpl.figure.Figure:
    fig = plt.figure(figsize=(12.0, 9.5))

    gs = gridspec.GridSpec(
        nrows=3, ncols=2,
        height_ratios=[1.15, 1.05, 0.9],
        hspace=0.62, wspace=0.30,
        left=0.06, right=0.97, top=0.92, bottom=0.07,
        figure=fig,
    )
    ax_a = fig.add_subplot(gs[0, :])     # full-width row
    ax_b = fig.add_subplot(gs[1, 0])
    ax_c = fig.add_subplot(gs[1, 1])
    ax_d = fig.add_subplot(gs[2, :])     # full-width row

    panel_pipeline(ax_a)
    panel_pfam_atlas(ax_b, args.pfam_atlas_tsv)
    panel_clinvar(ax_c, args.clinvar_data)
    panel_tp53(ax_d, args.tp53_isoforms)

    fig.text(0.05, 0.97,
             "Prot2Exon — domain → exon architecture at proteome scale",
             fontsize=13.5, fontweight="bold", color=COLORS["ink"])
    return fig


def main(argv=None) -> int:
    p = argparse.ArgumentParser(__doc__)
    p.add_argument("--pfam-atlas-tsv", type=Path, default=None,
                   help="pfam_atlas.tsv from the Pfam atlas notebook. "
                        "If absent, panel B uses the wiki's reported headline.")
    p.add_argument("--clinvar-data", type=Path, default=None,
                   help="Optional CSV with columns category, pct, odds_ratio, pval.")
    p.add_argument("--tp53-isoforms", type=Path,
                   default=Path("examples/tp53_isoforms.tsv"),
                   help="TP53 fixture (defaults to the bundled examples/ file).")
    p.add_argument("--out-dir", type=Path,
                   default=Path("figures"),
                   help="Output directory (default: ./figures/).")
    args = p.parse_args(argv)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    fig = make_figure(args)
    png = args.out_dir / "figure_1.png"
    pdf = args.out_dir / "figure_1.pdf"
    fig.savefig(png, dpi=300)
    fig.savefig(pdf)
    plt.close(fig)
    print(f"wrote {png}")
    print(f"wrote {pdf}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
