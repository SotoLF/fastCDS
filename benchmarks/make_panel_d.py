"""Build Figure 1 Panel D (PLAN.txt line 161): TP53 with Pfam domains
overlaid, generated through the prot2exon.plot() Python API.

This is the canonical "show the plotter works" deliverable.

Inputs:
  --index <human.idx>     prot2exon binary index (Phase 1 deliverable)

Outputs:
  --out <panel_d.png>     single-panel figure of TP53 isoform with domains
  --out-multi <panel_d_multi.png>  same data, multi-track (all 3 TP53 domains stacked)
"""

from __future__ import annotations

import argparse
import tempfile
from pathlib import Path

import prot2exon as p2g


# TP53 canonical isoform in Ensembl v86: ENSP00000269305 (TP53-201).
# Three Pfam-A domains (EnsDb.Hsapiens.v86 protein_domain table):
TP53_PROTEIN = "ENSP00000269305"
TP53_DOMAINS = [
    ("TAD",       "PF08563", 5, 29),    # transactivation domain
    ("DBD",       "PF00870", 95, 288),  # DNA-binding domain (the famous one)
    ("TETRAMER",  "PF07710", 319, 357), # tetramerization domain
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--index", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--out-multi", type=Path, default=None)
    ap.add_argument("--out-tsv", type=Path, default=None,
                    help="Optional: dump the underlying mapping result")
    args = ap.parse_args()

    # 1) Run prot2exon on all three TP53 domains via the Python wrapper.
    mapper = p2g.Mapper(args.index)
    queries = [{"protein_id": TP53_PROTEIN, "aa_start": s, "aa_end": e, "domain_id": name}
               for name, _pfam, s, e in TP53_DOMAINS]
    print(f"running prot2exon on {len(queries)} TP53 domains ...")
    result = mapper.map_batch(queries)

    print(f"got result with input_ids: {sorted(result.isoform['input_id'].unique())}")

    if args.out_tsv:
        result.isoform.to_csv(args.out_tsv, sep="\t", index=False)
        print(f"wrote {args.out_tsv}")

    # 2) Render the DBD (most paper-worthy) as a single panel.
    fig = p2g.plot(result, input_id="DBD", out=str(args.out),
                   title="TP53 — DNA-binding domain (PF00870, aa 95-288)",
                   width=10.0, height=2.6,
                   show_introns=True, show_utr=True, highlight_domain=True)
    print(f"wrote {args.out}")

    # 3) Optional: render all three domains as separate PNGs then stitch.
    if args.out_multi:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.image import imread
        out_multi_dir = args.out_multi.parent / "_panel_d_parts"
        out_multi_dir.mkdir(exist_ok=True)
        part_paths = []
        for name, pfam, s, e in TP53_DOMAINS:
            part = out_multi_dir / f"{name}.png"
            p2g.plot(result, input_id=name, out=str(part),
                     title=f"TP53 — {name} ({pfam}, aa {s}-{e})",
                     width=10.0, height=2.6,
                     show_introns=True, show_utr=True, highlight_domain=True)
            part_paths.append(part)
        # Stack vertically.
        imgs = [imread(p) for p in part_paths]
        fig, axes = plt.subplots(len(imgs), 1, figsize=(10, 2.6 * len(imgs)))
        for ax, img in zip(axes, imgs):
            ax.imshow(img)
            ax.set_axis_off()
        fig.tight_layout()
        fig.savefig(args.out_multi, dpi=150, bbox_inches="tight")
        print(f"wrote {args.out_multi}")


if __name__ == "__main__":
    main()
