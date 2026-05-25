"""Classify an external tool's per-query intervals against prot2exon's intervals
(taking prot2exon as the reference). Outputs the same Table-1-style agreement
breakdown as validate_vs_ensembldb.py, but for any tool whose output is in the
common 6-column TSV (query_id, chrom, start, end, strand, status).

Useful for comparing TransVar and Ensembl REST against the same reference
that ensembldb was compared to in Phase 2.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path


def normalize_chrom(c: str) -> str:
    return c[3:] if c.startswith("chr") else c


def interval_set(intervals):
    return {(normalize_chrom(c), s, e) for c, s, e in intervals}


def total_length(intervals):
    return sum(e - s + 1 for _, s, e in intervals)


def classify(ours, theirs) -> str:
    o_set = interval_set(ours)
    t_set = interval_set(theirs)
    if not o_set and not t_set:
        return "neither_mapped"
    if not t_set:
        return "only_prot2exon"
    if not o_set:
        return "only_external"
    if o_set == t_set:
        return "exact_match"
    if abs(total_length(ours) - total_length(theirs)) <= 2:
        return "off_by_one"
    return "structural_mismatch"


def load_p2e_intervals(cds_tsv: Path) -> dict:
    by_qid: dict = defaultdict(list)
    with open(cds_tsv) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["overlaps_domain"] != "coding_overlap":
                continue
            s, e = row["domain_overlap_genomic_start"], row["domain_overlap_genomic_end"]
            if s == "NA" or e == "NA":
                continue
            by_qid[row["input_id"]].append((row["chrom"], int(s), int(e)))
    return by_qid


def envelope(intervals):
    """Collapse a list of intervals to one (chrom, min_start, max_end). For
    multi-chrom (impossible here), keeps the first chrom."""
    if not intervals:
        return []
    chrom = intervals[0][0]
    return [(chrom, min(s for _, s, _ in intervals), max(e for _, _, e in intervals))]


def load_external(path: Path) -> dict:
    by_qid: dict = defaultdict(list)
    with open(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            qid = row["query_id"]
            if row["status"] != "ok":
                by_qid.setdefault(qid, [])
                continue
            by_qid[qid].append((row["chrom"], int(row["start"]), int(row["end"])))
    return by_qid


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--p2e-cds", required=True, type=Path,
                    help="prot2exon domain_cds_segments.tsv")
    ap.add_argument("--external", required=True, type=Path,
                    help="External tool TSV (query_id, chrom, start, end, strand, status)")
    ap.add_argument("--queries-meta", type=Path, default=None,
                    help="Optional queries metadata for category stratification")
    ap.add_argument("--tool-name", default="external")
    ap.add_argument("--envelope-only", action="store_true",
                    help="Collapse both sides to (chrom, min_start, max_end). Use for tools "
                         "that only emit a single genomic span per query (e.g. TransVar).")
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    ours = load_p2e_intervals(args.p2e_cds)
    theirs = load_external(args.external)
    if args.envelope_only:
        ours = {k: envelope(v) for k, v in ours.items()}
        theirs = {k: envelope(v) for k, v in theirs.items()}
    all_qids = list(theirs.keys())  # external query set defines the comparison

    cats: dict = {}
    if args.queries_meta:
        with open(args.queries_meta) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                cats[row["query_id"]] = row["category"]

    by_cat: dict = defaultdict(Counter)
    for qid in all_qids:
        bucket = classify(ours.get(qid, []), theirs.get(qid, []))
        cat = cats.get(qid, "UNKNOWN")
        by_cat["OVERALL"][bucket] += 1
        by_cat[cat][bucket] += 1

    BUCKETS = ["exact_match", "off_by_one", "structural_mismatch",
               "only_prot2exon", "only_external", "neither_mapped"]
    with open(args.out, "w") as f:
        f.write(f"# Tool: {args.tool_name}  |  ref: prot2exon  |  total queries: {len(all_qids):,}\n")
        f.write("category\tn\t" + "\t".join(BUCKETS) + "\texact_pct\n")
        for cat in ["OVERALL"] + sorted(k for k in by_cat if k != "OVERALL"):
            counts = by_cat[cat]
            n = sum(counts.values())
            considered = n - counts["neither_mapped"]
            exact_pct = (100.0 * counts["exact_match"] / considered) if considered else 0.0
            f.write(f"{cat}\t{n}\t" + "\t".join(str(counts[b]) for b in BUCKETS)
                    + f"\t{exact_pct:.2f}\n")
    print(f"wrote {args.out}")
    print(open(args.out).read())


if __name__ == "__main__":
    main()
