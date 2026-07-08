#!/usr/bin/env python3
"""Render the same isoform three ways (matplotlib, plotly, vanilla JS) plus a
static multi-isoform stack.

Uses the in-repo TP53 fixture, so it runs with no index and no network:

    python tutorial/examples/make_plot_gallery.py            # -> ./plot_gallery/
    python tutorial/examples/make_plot_gallery.py --out-dir /tmp/gallery

Outputs:
  * plot_matplotlib.png            static single-isoform figure
  * plot_js.html                   self-contained vanilla-JS viewer
  * plot_plotly.html               plotly viewer
  * plot_isoform_stack_static.png  every isoform stacked on one genomic axis
"""
from __future__ import annotations

import argparse
from pathlib import Path

import fastCDS as fc


def render_static_stack(src: str, out_png: Path) -> None:
    """Draw every isoform in `src` stacked on one shared genomic axis (the
    static matplotlib counterpart of the interactive stack viewer)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from fastCDS.plot import load_isoform_tsv, _draw_genomic, _legend

    by_id = load_isoform_tsv(src)
    ids = list(by_id)                      # fixture order: canonical first
    all_segs = [s for i in ids for s in by_id[i]]
    xmin = min(s.feature_genomic_start for s in all_segs)
    xmax = max(s.feature_genomic_end for s in all_segs)
    pad = (xmax - xmin) * 0.02 or 1
    n = len(ids)
    fig, axes = plt.subplots(n, 1, figsize=(12, 1.15 * n + 0.9),
                             sharex=True, constrained_layout=True)
    axes = [axes] if n == 1 else list(axes)
    for ax, i in zip(axes, ids):
        segs = by_id[i]
        _draw_genomic(ax, segs, show_introns=True, show_utr=True,
                      highlight_domain=True)
        ax.set_xlim(xmin - pad, xmax + pad)   # shared axis (override per-isoform)
        ax.set_xlabel("")
        ax.set_ylabel(segs[0].transcript_id, rotation=0, ha="right",
                      va="center", fontsize=8)
    axes[-1].set_xlabel(f"genomic position ({all_segs[0].chrom})")
    _legend(fig, axes[0], show_introns=True, show_utr=True, highlight_domain=True)
    fig.suptitle(f"{all_segs[0].gene_name} - DBD across {n} isoforms",
                 fontsize=11)
    fig.savefig(out_png, dpi=150)
    plt.close(fig)

FIXTURE = Path(__file__).with_name("tp53_isoforms.tsv")
INPUT_ID = "TP53_canonical"
LINK = ("https://www.ensembl.org/Homo_sapiens/Transcript/"
        "ProteinSummary?p={protein_id}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", default="plot_gallery",
                    help="where to write the three files (default: ./plot_gallery)")
    ap.add_argument("--isoform", default=str(FIXTURE),
                    help="isoform_structure.tsv to render (default: bundled TP53)")
    ap.add_argument("--input-id", default=INPUT_ID,
                    help=f"which query to draw (default: {INPUT_ID})")
    args = ap.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    src, iid = args.isoform, args.input_id

    # 1. static matplotlib — PNG, embed as-is
    fc.plot(src, input_id=iid, out=str(out / "plot_matplotlib.png"))

    # 2. interactive, self-contained vanilla JS (default engine)
    fc.plot(src, input_id=iid, out=str(out / "plot_js.html"), link_template=LINK)

    # 3. interactive, plotly engine
    try:
        fc.plot(src, input_id=iid, out=str(out / "plot_plotly.html"),
                engine="plotly", link_template=LINK)
    except SystemExit as e:   # plotly not installed
        print(f"skipped plotly panel: {e}")

    # 4. static multi-isoform stack (all isoforms in the TSV on one axis)
    render_static_stack(src, out / "plot_isoform_stack_static.png")

    print(f"wrote gallery to {out.resolve()}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
