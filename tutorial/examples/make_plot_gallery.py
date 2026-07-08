#!/usr/bin/env python3
"""Render the same isoform three ways — matplotlib, plotly, vanilla JS.

Produces the three panels shown on the Plotting wiki page so they can be
screenshotted (static PNG) or screen-recorded (the two interactive HTMLs).
Uses the in-repo TP53 fixture, so it runs with no index and no network:

    python tutorial/examples/make_plot_gallery.py            # -> ./plot_gallery/
    python tutorial/examples/make_plot_gallery.py --out-dir /tmp/gallery

Then:
  * plot_matplotlib.png  — ready to embed as-is
  * plot_js.html         — open in a browser, screenshot or record a GIF
  * plot_plotly.html     - same (plotly ships with fastCDS)
"""
from __future__ import annotations

import argparse
from pathlib import Path

import fastCDS as fc

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

    print(f"wrote gallery to {out.resolve()}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
