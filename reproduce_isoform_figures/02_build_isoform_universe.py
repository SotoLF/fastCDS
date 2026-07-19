#!/usr/bin/env python3
"""Step 2 - list protein-coding isoforms from the Ensembl GTF.

Keeps genes that have 2+ coding isoforms and picks one reference isoform per
gene: MANE Select, else Ensembl canonical, else the longest CDS.

Outputs
  isoforms.tsv  one row per coding isoform
  genes.tsv     one row per gene, with its reference isoform

Run
  python 02_build_isoform_universe.py \
      --gtf Homo_sapiens.GRCh38.115.gtf.gz \
      --out-isoforms isoforms.tsv --out-genes genes.tsv
"""
import argparse
import gzip
import re
from collections import defaultdict
from pathlib import Path

# GTF column 9 is `key "value"; key "value"; ...`
ATTR = re.compile(r'(\w+) "([^"]*)"')


def attrs(field):
    return dict(ATTR.findall(field))


def opener(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def pick_reference(isoforms):
    """MANE Select > Ensembl canonical > longest CDS."""
    for r in isoforms:
        if r["is_mane_select"]:
            return r["protein_id"], "mane_select"
    for r in isoforms:
        if r["is_ensembl_canonical"]:
            return r["protein_id"], "ensembl_canonical"
    longest = max(isoforms, key=lambda r: r["cds_len"])
    return longest["protein_id"], "longest_cds"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gtf", required=True, type=Path)
    ap.add_argument("--out-isoforms", required=True, type=Path)
    ap.add_argument("--out-genes", required=True, type=Path)
    ap.add_argument("--min-isoforms", type=int, default=2)
    args = ap.parse_args()

    tx = {}                        # transcript_id -> record
    cds_len = defaultdict(int)     # transcript_id -> total CDS length
    cds_exons = defaultdict(int)   # transcript_id -> number of CDS exons
    protein = {}                   # transcript_id -> protein_id

    # read the GTF once, adding up CDS length and exon count per transcript
    with opener(args.gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9:
                continue
            if c[2] == "transcript":
                a = attrs(c[8])
                if a.get("transcript_biotype") != "protein_coding":
                    continue
                tx[a["transcript_id"]] = dict(
                    gene_id=a.get("gene_id", ""), gene_name=a.get("gene_name", ""),
                    transcript_id=a["transcript_id"], chrom=c[0], strand=c[6],
                    is_mane_select='tag "MANE_Select"' in c[8],
                    is_ensembl_canonical='tag "Ensembl_canonical"' in c[8])
            elif c[2] == "CDS":
                a = attrs(c[8])
                tid = a.get("transcript_id")
                if tid is None:
                    continue
                cds_len[tid] += int(c[4]) - int(c[3]) + 1
                cds_exons[tid] += 1
                protein.setdefault(tid, a.get("protein_id"))

    # keep only transcripts that translate (have CDS + protein_id), group by gene
    by_gene = defaultdict(list)
    for tid, rec in tx.items():
        if tid not in cds_len or not protein.get(tid):
            continue
        rec = dict(rec, protein_id=protein[tid],
                   cds_len=cds_len[tid], n_cds_exons=cds_exons[tid])
        by_gene[rec["gene_id"]].append(rec)
    multi = {g: v for g, v in by_gene.items() if len(v) >= args.min_isoforms}

    cols = ["gene_id", "gene_name", "transcript_id", "protein_id", "chrom", "strand",
            "is_mane_select", "is_ensembl_canonical", "cds_len", "n_cds_exons"]
    with open(args.out_isoforms, "w") as out:
        out.write("\t".join(cols) + "\n")
        for v in multi.values():
            for r in v:
                out.write("\t".join(str(r[k]) for k in cols) + "\n")

    with open(args.out_genes, "w") as out:
        out.write("gene_id\tgene_name\tn_isoforms\tsource_protein_id\tsource_basis\n")
        for g, v in multi.items():
            src, basis = pick_reference(v)
            out.write(f"{g}\t{v[0]['gene_name']}\t{len(v)}\t{src}\t{basis}\n")

    n_iso = sum(len(v) for v in multi.values())
    print(f"coding isoforms: {n_iso:,}   genes (2+ isoforms): {len(multi):,}")


if __name__ == "__main__":
    main()
