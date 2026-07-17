"""geneplot route: per-gene transcript-to-genome walk.

geneplot has no index. It builds a gffutils SQLite database of the genome on
first use, then re-reads the domain file for every gene, so the per-query cost
stays flat and high. Timed end-to-end like every other runner.

Usage:
    python run_geneplot.py <N> <ids.txt> --gff h86.gff3 --ipr h86.ipr \
        --ensp-enst ensp_enst.tsv [--geneplot-src DIR]

Build the --gff / --ipr / --ensp-enst inputs with build_human_tool_inputs.py.
"""

from __future__ import annotations

import argparse
import contextlib
import io
import logging
import sys
from pathlib import Path


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("n", type=int, help="Number of queries to map")
    ap.add_argument("ids", type=Path, help="ENSP ids, one per line")
    ap.add_argument("--gff", required=True, type=Path, help="GFF3 annotation")
    ap.add_argument("--ipr", required=True, type=Path, help="InterProScan .ipr domains")
    ap.add_argument("--ensp-enst", required=True, type=Path,
                    help="TSV: ENSP <tab> ENST")
    ap.add_argument("--geneplot-src", type=Path, default=None,
                    help="geneplot checkout to import from (default: installed package)")
    args = ap.parse_args()

    for p, what in [(args.ids, "ids"), (args.gff, "--gff"),
                    (args.ipr, "--ipr"), (args.ensp_enst, "--ensp-enst")]:
        if not p.exists():
            raise SystemExit(f"run_geneplot.py: {what} not found: {p}")

    if args.geneplot_src:
        sys.path.insert(0, str(args.geneplot_src))
    logging.getLogger("geneplot").setLevel(logging.ERROR)
    import geneplot as gp

    ids = [l.strip() for l in open(args.ids) if l.strip()][: args.n]

    ensp_to_enst = {}
    for line in open(args.ensp_enst):
        parts = line.split()
        if len(parts) >= 2:
            ensp_to_enst[parts[0]] = parts[1]

    # geneplot prints progress to stdout per call; silence it so the log stays
    # readable without changing what is executed.
    with contextlib.redirect_stdout(io.StringIO()):
        genome = gp.genome(str(args.gff), iprfile=str(args.ipr), vcffiles="./")

    mapped = 0
    for pid in ids:
        enst = ensp_to_enst.get(pid)
        if not enst:
            continue
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                gene = genome.gene(mRNAid="transcript:" + enst, proteinid=pid)
                gene._proteindoms(str(args.ipr), pid)
                gene._transcriptpos_to_genomepos()
            mapped += 1
        except Exception:
            pass

    print(f"geneplot mapped {mapped} / {len(ids)}")


if __name__ == "__main__":
    main()
