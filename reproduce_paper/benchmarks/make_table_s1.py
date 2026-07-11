"""Build Supplementary Table S1 — per-tool, per-category agreement with fastCDS.

Assembles one combined table across the four comparators from the per-category
agreement tables written by classify_external.py (one per tool), applying two
methodology rules:

  * TransVar returns a single enclosing genomic span, not the individual CDS
    intervals, so its coordinate agreement can only be scored at fastCDS's
    resolution for categories in which EVERY query maps to a single CDS block.
    Multi-block categories are reported as NA. Answerability is derived from the
    fastCDS output (number of coding_overlap blocks per query).

  * Ensembl REST is broken out in full: exact, off-by-one (<= 2 nt), and
    no-mapping (REST returned an error / no result). ensembldb and
    GenomicFeatures agree exactly everywhere, so only their percentage is shown.

The categories are the nine EXCLUSIVE strata from sample_validation_queries.py
(cds_incomplete is exclusive; the other eight are complete-CDS only).

Usage:
    python make_table_s1.py \
        --fastcds-cds       exclusive_v115/p2e/domain_cds_segments.tsv \
        --ensembldb-table   exclusive_v86/ensembldb_agreement.tsv \
        --genomicfeatures-table exclusive_v86/genomicfeatures_agreement.tsv \
        --transvar-table    exclusive_v95/transvar_agreement.tsv \
        --rest-table        exclusive_v115/rest_agreement.tsv \
        --out               table_S1_agreement.tsv
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

CATEGORY_ORDER = [
    "single_exon_domain", "multi_exon_domain", "codon_split_boundary",
    "plus_strand_gene", "minus_strand_gene", "cds_incomplete",
    "selenoprotein", "single_exon_gene", "many_exon_gene", "OVERALL",
]


def read_agreement(path: Path) -> dict:
    """category -> row dict from a classify_external.py agreement table
    (skips the leading '# Tool: ...' comment line)."""
    out = {}
    if path is None:
        return out
    with open(path) as f:
        lines = [ln for ln in f if not ln.startswith("#")]
    for row in csv.DictReader(lines, delimiter="\t"):
        out[row["category"]] = row
    return out


def transvar_answerable_categories(fastcds_cds: Path) -> set:
    """A category is TransVar-answerable iff every one of its queries maps to a
    single CDS block (so the single span equals the per-exon structure)."""
    blocks = defaultdict(set)
    cat_of = {}
    with open(fastcds_cds) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r.get("overlaps_domain") != "coding_overlap":
                continue
            s, e = r["domain_overlap_genomic_start"], r["domain_overlap_genomic_end"]
            if s == "NA" or e == "NA":
                continue
            qid = r["input_id"]
            blocks[qid].add((r["chrom"], s, e))
            cat = qid
            # category is the suffix after Q000000_
            if "_" in qid:
                cat = qid.split("_", 1)[1]
            cat_of[qid] = cat
    multiblock_cats = set()
    all_cats = set()
    for qid, blk in blocks.items():
        c = cat_of[qid]
        all_cats.add(c)
        if len(blk) > 1:
            multiblock_cats.add(c)
    return all_cats - multiblock_cats


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fastcds-cds", required=True, type=Path,
                    help="fastCDS domain_cds_segments.tsv (for TransVar answerability)")
    ap.add_argument("--ensembldb-table", required=True, type=Path)
    ap.add_argument("--genomicfeatures-table", required=True, type=Path)
    ap.add_argument("--transvar-table", required=True, type=Path)
    ap.add_argument("--rest-table", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    ens = read_agreement(args.ensembldb_table)
    gf = read_agreement(args.genomicfeatures_table)
    tv = read_agreement(args.transvar_table)
    rest = read_agreement(args.rest_table)
    tv_ok = transvar_answerable_categories(args.fastcds_cds)

    def pct(tbl, cat):
        r = tbl.get(cat)
        return "NA" if r is None else str(round(float(r["exact_pct"]), 1))

    header = ["category", "n", "ensembldb_pct", "GenomicFeatures_pct",
              "TransVar_pct", "REST_exact", "REST_off", "REST_no_map", "REST_pct"]
    rows = [header]
    for cat in CATEGORY_ORDER:
        r = rest.get(cat)
        if r is None:
            continue
        n = r["n"]
        # TransVar: only categories where every query is single-block (plus OVERALL is NA).
        tv_pct = pct(tv, cat) if (cat in tv_ok and cat != "OVERALL") else "NA"
        rows.append([
            cat, n, pct(ens, cat), pct(gf, cat), tv_pct,
            r["exact_match"], r["off_by_one"], r["only_fastCDS"],
            str(round(float(r["exact_pct"]), 2)),
        ])

    with open(args.out, "w") as f:
        for row in rows:
            f.write("\t".join(map(str, row)) + "\n")
    print(f"wrote {args.out}\n")
    widths = [max(len(str(row[i])) for row in rows) for i in range(len(header))]
    for row in rows:
        print("  ".join(str(c).ljust(w) for c, w in zip(row, widths)))


if __name__ == "__main__":
    main()
