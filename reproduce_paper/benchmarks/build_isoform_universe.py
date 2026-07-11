#!/usr/bin/env python3
"""Build the protein-coding isoform universe for the conservation analysis.

Reads an Ensembl/GENCODE GTF and, for every protein-coding transcript that has a
CDS, records one row: gene, transcript, protein, the MANE-Select / Ensembl-
canonical flags, CDS length and CDS-exon count. The conservation analysis then
keeps only genes with >= 2 such isoforms (one "source" + >=1 "target").

Why parse the GTF directly (instead of the prebuilt index)? We need the gene ->
isoforms grouping plus the MANE / canonical tags to choose each domain's source
isoform; those tags live in the GTF attributes.

Usage:
  build_isoform_universe.py --gtf F.gtf --out-isoforms iso.tsv --out-genes genes.tsv
"""
from __future__ import annotations

import argparse
import re
from collections import defaultdict
from pathlib import Path

# attribute helpers -- GTF column 9 is `key "value"; key "value"; ...`
_ATTR = re.compile(r'(\w+) "([^"]*)"')


def _attrs(field: str) -> dict:
    return dict(_ATTR.findall(field))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gtf", required=True, type=Path)
    ap.add_argument("--out-isoforms", required=True, type=Path,
                    help="One row per protein-coding isoform")
    ap.add_argument("--out-genes", required=True, type=Path,
                    help="One row per multi-isoform gene (summary)")
    ap.add_argument("--min-isoforms", type=int, default=2,
                    help="Keep genes with at least this many coding isoforms")
    ap.add_argument("--appris", type=Path, default=None,
                    help="Optional, not used in the paper (the default source pick is "
                         "MANE Select / Ensembl canonical). APPRIS principal-isoforms TSV "
                         "(appris_data.principal.txt). "
                         "When given, the source isoform is the APPRIS PRINCIPAL "
                         "transcript (best-ranked), falling back to MANE/canonical/"
                         "longest only for genes APPRIS does not cover.")
    args = ap.parse_args()

    appris = _load_appris(args.appris) if args.appris else {}

    # transcript_id -> record; CDS lengths/exon counts accumulate as we stream
    tx = {}
    cds_len = defaultdict(int)
    cds_exons = defaultdict(int)
    protein = {}

    with open(args.gtf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9:
                continue
            feature = c[2]
            if feature == "transcript":
                a = _attrs(c[8])
                if a.get("transcript_biotype") != "protein_coding":
                    continue
                tags = c[8]  # tags repeat; substring test is enough
                tx[a["transcript_id"]] = dict(
                    gene_id=a.get("gene_id", ""),
                    gene_name=a.get("gene_name", ""),
                    transcript_id=a["transcript_id"],
                    chrom=c[0], strand=c[6],
                    is_mane_select='tag "MANE_Select"' in tags,
                    is_ensembl_canonical='tag "Ensembl_canonical"' in tags,
                )
            elif feature == "CDS":
                a = _attrs(c[8])
                tid = a.get("transcript_id")
                if tid is None:
                    continue
                cds_len[tid] += int(c[4]) - int(c[3]) + 1
                cds_exons[tid] += 1
                if tid not in protein and "protein_id" in a:
                    protein[tid] = a["protein_id"]

    # keep only translatable transcripts (have CDS + protein_id)
    iso_by_gene = defaultdict(list)
    rows = []
    for tid, rec in tx.items():
        if tid not in cds_len or tid not in protein:
            continue
        rec = dict(rec, protein_id=protein[tid],
                   cds_len=cds_len[tid], n_cds_exons=cds_exons[tid])
        rows.append(rec)
        iso_by_gene[rec["gene_id"]].append(rec)

    multi = {g: v for g, v in iso_by_gene.items() if len(v) >= args.min_isoforms}

    cols = ["gene_id", "gene_name", "transcript_id", "protein_id", "chrom",
            "strand", "is_mane_select", "is_ensembl_canonical", "cds_len",
            "n_cds_exons"]
    with open(args.out_isoforms, "w") as out:
        out.write("\t".join(cols) + "\n")
        for r in rows:
            if r["gene_id"] not in multi:
                continue
            out.write("\t".join(str(r[k]) for k in cols) + "\n")

    # gene summary: how many isoforms, and which one is the chosen source
    with open(args.out_genes, "w") as out:
        out.write("gene_id\tgene_name\tn_isoforms\tn_mane\tn_canonical\t"
                  "source_protein_id\tsource_basis\n")
        for g, v in multi.items():
            source, basis = _pick_source(v, appris)
            out.write(f"{g}\t{v[0]['gene_name']}\t{len(v)}\t"
                      f"{sum(x['is_mane_select'] for x in v)}\t"
                      f"{sum(x['is_ensembl_canonical'] for x in v)}\t"
                      f"{source}\t{basis}\n")

    n_iso = sum(len(v) for v in multi.values())
    print(f"protein-coding isoforms (multi-isoform genes only): {n_iso:,}")
    print(f"multi-isoform genes (>= {args.min_isoforms}): {len(multi):,}")
    print(f"  isoforms TSV: {args.out_isoforms}")
    print(f"  genes  TSV:   {args.out_genes}")


def _load_appris(path):
    """transcript_id (unversioned ENST) -> APPRIS label (PRINCIPAL:n / ALTERNATIVE:n)."""
    out = {}
    with open(path) as f:
        next(f, None)  # header row
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) < 5:
                continue
            tid = c[2].split(".")[0]
            out[tid] = c[4]
    return out


def _appris_rank(label):
    """Lower is more-principal: PRINCIPAL:1 < ... < PRINCIPAL:M < ALTERNATIVE < none."""
    if not label:
        return 99
    kind, _, n = label.partition(":")
    base = 0 if kind == "PRINCIPAL" else 10
    return base + (6 if n == "M" else int(n) if n.isdigit() else 9)


def _pick_source(isoforms, appris=None):
    """Choose the reference isoform.

    With APPRIS: the best-ranked PRINCIPAL transcript (tie-break longest CDS).
    Genes APPRIS does not cover fall back to MANE Select > Ensembl canonical >
    longest CDS. Without APPRIS: MANE > canonical > longest.
    """
    if appris:
        ranked = [(_appris_rank(appris.get(r["transcript_id"], "")), -r["cds_len"], r)
                  for r in isoforms]
        best_rank = min(ranked, key=lambda t: (t[0], t[1]))
        if best_rank[0] < 99:  # at least one isoform has an APPRIS label
            r = best_rank[2]
            return r["protein_id"], "appris_" + appris[r["transcript_id"]].replace(":", "").lower()
    for r in isoforms:
        if r["is_mane_select"]:
            return r["protein_id"], "mane_select"
    for r in isoforms:
        if r["is_ensembl_canonical"]:
            return r["protein_id"], "ensembl_canonical"
    longest = max(isoforms, key=lambda r: r["cds_len"])
    return longest["protein_id"], "longest_cds"


if __name__ == "__main__":
    main()
