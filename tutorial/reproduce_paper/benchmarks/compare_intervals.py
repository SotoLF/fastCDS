#!/usr/bin/env python3
"""Compare per-query genomic intervals across protein-to-genome mappers.

Each tool reduces to: query_id -> set of (chrom, start, end, strand) coding
segments. Two tools "agree" on a query when those sets are identical. Reports
pairwise exact-segment agreement over the queries both tools mapped.

Usage:
  compare_intervals.py ensembldb=<tsv> genomicfeatures=<tsv> fastCDS=<coding_tsv>

The ensembldb/genomicfeatures TSVs have columns
  query_id chrom start end strand status
The fastCDS arg is a domain_cds_segments.tsv; its coding_overlap rows supply
(input_id, chrom, domain_overlap_genomic_start, domain_overlap_genomic_end, strand).
"""
import sys, csv
from collections import defaultdict


def load_rgranges(path):
    d = defaultdict(set)
    with open(path) as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            if row["status"] != "ok":
                continue
            d[row["query_id"]].add(
                (row["chrom"], int(row["start"]), int(row["end"]), row["strand"]))
    return d


def load_fastCDS(path):
    d = defaultdict(set)
    with open(path) as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            if row.get("overlaps_domain") != "coding_overlap":
                continue
            d[row["input_id"]].add(
                (row["chrom"], int(row["domain_overlap_genomic_start"]),
                 int(row["domain_overlap_genomic_end"]), row["strand"]))
    return d


def main():
    tools = {}
    for arg in sys.argv[1:]:
        name, path = arg.split("=", 1)
        tools[name] = load_fastCDS(path) if name == "fastCDS" else load_rgranges(path)

    names = list(tools)
    print(f"loaded: " + ", ".join(f"{n}={len(tools[n])} queries" for n in names))
    print()
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            a, b = names[i], names[j]
            da, db = tools[a], tools[b]
            common = set(da) & set(db)
            if not common:
                print(f"{a} vs {b}: no common queries"); continue
            agree = sum(1 for q in common if da[q] == db[q])
            pct = 100.0 * agree / len(common)
            print(f"{a:16s} vs {b:16s}: {agree}/{len(common)} exact-segment match ({pct:.2f}%)")
            # show up to 3 disagreements
            diffs = [q for q in common if da[q] != db[q]][:3]
            for q in diffs:
                print(f"    DIFF {q}: {a}={sorted(da[q])}  {b}={sorted(db[q])}")


if __name__ == "__main__":
    main()
