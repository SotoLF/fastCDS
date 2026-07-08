#!/usr/bin/env python3
"""Build a full-transcript BED12 track from an `isoform_structure.tsv`.

`fastCDS map` writes a BED12 for the *domain* (`domain_blocks.bed`); this
companion turns the full isoform structure into a **transcript** BED12 — the
whole gene model with UTRs drawn thin and CDS drawn thick — so you can load it
in IGV/UCSC next to the domain track.

    python tutorial/examples/isoform_to_transcript_bed12.py \
        results/isoform_structure.tsv --out transcript.bed

One BED12 row per transcript (grouped by input_id). Genomic coordinates in the
TSV are 1-based inclusive; BED is 0-based half-open, handled here.
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict

EXONIC = {"five_prime_UTR", "three_prime_UTR", "CDS"}


def _merge_exons(feats):
    """Merge genomically contiguous exonic features (UTR+CDS in the same exon
    are adjacent; introns leave a gap) into (start, end) exon blocks, 1-based
    inclusive."""
    feats = sorted(feats, key=lambda f: f[0])
    blocks = []
    for start, end in feats:
        if blocks and start <= blocks[-1][1] + 1:
            blocks[-1][1] = max(blocks[-1][1], end)
        else:
            blocks.append([start, end])
    return blocks


def transcript_rows(path, *, rgb="30,120,180"):
    by_id = defaultdict(list)
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            by_id[row["input_id"]].append(row)

    for input_id, rows in by_id.items():
        exonic = [(int(r["feature_genomic_start"]), int(r["feature_genomic_end"]))
                  for r in rows if r["feature_type"] in EXONIC]
        if not exonic:
            continue
        cds = [(int(r["feature_genomic_start"]), int(r["feature_genomic_end"]))
               for r in rows if r["feature_type"] == "CDS"]
        r0 = rows[0]
        chrom, strand = r0["chrom"], r0["strand"] or "+"

        blocks = _merge_exons(exonic)                 # 1-based inclusive
        chrom_start = min(b[0] for b in blocks) - 1   # -> 0-based
        chrom_end = max(b[1] for b in blocks)         # -> 0-based half-open
        if cds:
            thick_start = min(c[0] for c in cds) - 1
            thick_end = max(c[1] for c in cds)
        else:                                         # non-coding: no thick part
            thick_start = thick_end = chrom_start

        sizes = [b[1] - b[0] + 1 for b in blocks]
        starts = [(b[0] - 1) - chrom_start for b in blocks]
        name = f"{r0.get('gene_name') or ''}_{r0['transcript_id']}".strip("_") \
            or input_id

        yield [
            chrom, chrom_start, chrom_end, name, 0, strand,
            thick_start, thick_end, rgb, len(blocks),
            ",".join(map(str, sizes)) + ",",
            ",".join(map(str, starts)) + ",",
        ]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("isoform", help="path to isoform_structure.tsv")
    ap.add_argument("--out", default="-", help="output .bed (default: stdout)")
    ap.add_argument("--rgb", default="30,120,180",
                    help="itemRgb for the transcript blocks (default: blue)")
    ap.add_argument("--track-name", default="transcript (fastCDS)",
                    help="IGV/UCSC track line name")
    args = ap.parse_args()

    out = sys.stdout if args.out == "-" else open(args.out, "w")
    out.write(f'track name="{args.track_name}" itemRgb="On"\n')
    for row in transcript_rows(args.isoform, rgb=args.rgb):
        out.write("\t".join(map(str, row)) + "\n")
    if out is not sys.stdout:
        out.close()
        print(f"wrote {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
