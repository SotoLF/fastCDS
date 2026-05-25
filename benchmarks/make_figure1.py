"""Compose the 4-panel Figure 1 of the prot2exon paper (PLAN.txt lines 153-161).

  Panel A — Pipeline overview (InterPro/Pfam -> prot2exon -> atlas) [drawn here]
  Panel B — Histogram of n_coding_exons per domain        [from Notebook 1 output]
  Panel C — ClinVar pathogenic enrichment plot            [from Notebook 2 output]
  Panel D — TP53 DBD with prot2exon.plot() overlay        [make_panel_d.py output]

Inputs are all PNGs already produced upstream; this script just lays them out
and renders the pipeline diagram for Panel A.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.image import imread
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch


def render_pipeline_panel(ax) -> None:
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 4)
    ax.set_axis_off()

    boxes = [
        (0.5, 1.5, 2.2, 1.0, "InterPro / Pfam\ndomain calls\n(protein coords)", "#cfe2ff"),
        (3.5, 1.5, 2.2, 1.0, "prot2exon\n(C++ + Python API)", "#ffd9a6"),
        (6.7, 2.6, 2.7, 0.9, "Proteome atlas\n(137K domains)", "#d4edda"),
        (6.7, 0.4, 2.7, 0.9, "Per-domain exon /\nintron architecture", "#d4edda"),
    ]
    for x, y, w, h, label, color in boxes:
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.05",
                                    fc=color, ec="black", lw=1.2))
        ax.text(x + w / 2, y + h / 2, label, ha="center", va="center", fontsize=9)

    arrows = [
        ((2.7, 2.0), (3.5, 2.0)),
        ((5.7, 2.3), (6.7, 3.05)),
        ((5.7, 1.7), (6.7, 0.85)),
    ]
    for start, end in arrows:
        ax.annotate("", xy=end, xytext=start,
                    arrowprops=dict(arrowstyle="->", lw=1.4, color="black"))

    ax.text(5.0, 3.6, "Pipeline overview", ha="center", fontsize=11, fontweight="bold")


def render_image_panel(ax, image_path: Path, title: str) -> None:
    ax.imshow(imread(image_path))
    ax.set_axis_off()
    ax.set_title(title, fontsize=11, fontweight="bold", pad=6)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel-b", required=True, type=Path,
                    help="Notebook 1 atlas figure (multi-panel histograms)")
    ap.add_argument("--panel-c", required=True, type=Path,
                    help="Notebook 2 ClinVar enrichment figure")
    ap.add_argument("--panel-d", required=True, type=Path,
                    help="TP53 DBD plot from make_panel_d.py")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--out-pdf", type=Path, default=None)
    args = ap.parse_args()

    fig = plt.figure(figsize=(16, 10))
    # 2x2 grid: A (top-left), B (top-right), C (bottom-left), D (bottom-right)
    gs = fig.add_gridspec(2, 2, hspace=0.25, wspace=0.15)

    ax_a = fig.add_subplot(gs[0, 0])
    render_pipeline_panel(ax_a)
    ax_a.text(-0.04, 1.02, "A", transform=ax_a.transAxes,
              fontsize=18, fontweight="bold", va="bottom")

    ax_b = fig.add_subplot(gs[0, 1])
    render_image_panel(ax_b, args.panel_b,
                       "Domain encoding architecture across 137K Pfam-A instances")
    ax_b.text(-0.04, 1.02, "B", transform=ax_b.transAxes,
              fontsize=18, fontweight="bold", va="bottom")

    ax_c = fig.add_subplot(gs[1, 0])
    render_image_panel(ax_c, args.panel_c,
                       "ClinVar pathogenic enrichment in Pfam-coding CDS exons")
    ax_c.text(-0.04, 1.02, "C", transform=ax_c.transAxes,
              fontsize=18, fontweight="bold", va="bottom")

    ax_d = fig.add_subplot(gs[1, 1])
    render_image_panel(ax_d, args.panel_d,
                       "TP53 DBD: domain projected onto the gene structure (prot2exon.plot)")
    ax_d.text(-0.04, 1.02, "D", transform=ax_d.transAxes,
              fontsize=18, fontweight="bold", va="bottom")

    fig.suptitle("Figure 1 — prot2exon: pipeline, atlas, clinical enrichment, worked example",
                 fontsize=13, fontweight="bold", y=0.99)
    fig.savefig(args.out, dpi=200, bbox_inches="tight")
    print(f"wrote {args.out}")
    if args.out_pdf:
        fig.savefig(args.out_pdf, bbox_inches="tight")
        print(f"wrote {args.out_pdf}")


if __name__ == "__main__":
    main()
