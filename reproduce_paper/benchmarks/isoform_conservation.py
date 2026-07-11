#!/usr/bin/env python3
"""Project Pfam domains across isoforms and score conservation (Stages 2-3-4).

For every Pfam domain that sits on a gene's *source* isoform (MANE Select /
Ensembl canonical), project it to genomic coordinates with fastCDS, then ask of
every *other* protein-coding isoform of that gene: are the domain's coding bases
retained as CDS, and in the same reading frame?

  coverage          = retained coding bp / domain bp
  inframe_coverage  = retained AND same-frame coding bp / domain bp
  conservation_class= conserved (>=0.80) / partial (0.50-0.80) / lost (<0.50)
  mechanism         = why it's not conserved (see classify_mechanism.py)

The frame test is cheap: over any contiguous source-domain x target-CDS overlap,
delta = (target cds_nt) - (source cds_nt) is constant, so one `delta % 3 == 0`
check classifies the whole sub-interval.

Usage:
  isoform_conservation.py --isoforms iso.tsv --genes genes.tsv \
      --pfam-meta pfam_human_v115_meta.tsv --index ensembl_v115_human.idx \
      --out conservation_long.tsv [--genes-limit N] [--threads 4]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd
import fastCDS as fc

sys.path.insert(0, str(Path(__file__).resolve().parent))
import classify_mechanism as cm


def _cds_at(g, G0, G1, Clo, strand):
    """CDS-nt coordinate of genomic base g inside a CDS exon [G0,G1] (Clo = its
    smallest cds_nt). Genomic ascending == cds ascending on '+', descending on '-'."""
    return Clo + (g - G0) if strand == "+" else Clo + (G1 - g)


def _segments(seg_df):
    """query_id -> dict(strand, intervals=[(g0,g1)], segs=[(g0,g1,clo)])."""
    out = {}
    sub = seg_df[["input_id", "strand", "domain_overlap_genomic_start",
                  "domain_overlap_genomic_end", "domain_overlap_cds_nt_start"]]
    for qid, strand, g0, g1, clo in sub.itertuples(index=False):
        d = out.get(qid)
        if d is None:
            d = out[qid] = dict(strand=strand, intervals=[], segs=[])
        g0, g1, clo = int(g0), int(g1), int(clo)
        d["intervals"].append((g0, g1))
        d["segs"].append((g0, g1, clo))
    return out


def _structures(struct_df):
    """target ensp -> dict(strand, cds_frame=[(g0,g1,clo)], cds/five_utr/three_utr/
    intron=[(g0,g1)], tx_start, tx_end)."""
    out = {}
    keep = struct_df[struct_df.feature_type.isin(
        ["CDS", "five_prime_UTR", "three_prime_UTR", "intron"])]
    sub = keep[["input_id", "strand", "feature_type", "feature_genomic_start",
                "feature_genomic_end", "cds_nt_start"]]
    for ensp, strand, ft, g0, g1, clo in sub.itertuples(index=False):
        d = out.get(ensp)
        if d is None:
            d = out[ensp] = dict(strand=strand, cds_frame=[], cds=[], five_utr=[],
                                 three_utr=[], intron=[], tx_start=10**18, tx_end=0)
        g0, g1 = int(g0), int(g1)
        d["tx_start"] = min(d["tx_start"], g0)
        d["tx_end"] = max(d["tx_end"], g1)
        if ft == "CDS":
            d["cds"].append((g0, g1))
            d["cds_frame"].append((g0, g1, int(clo)))
        elif ft == "five_prime_UTR":
            d["five_utr"].append((g0, g1))
        elif ft == "three_prime_UTR":
            d["three_utr"].append((g0, g1))
        elif ft == "intron":
            d["intron"].append((g0, g1))
    return out


def _score(dom, target):
    """Return (domain_bp, covered_bp, inframe_bp) for one domain vs one isoform."""
    strand = dom["strand"]
    domain_bp = sum(g1 - g0 + 1 for g0, g1 in dom["intervals"])
    covered = inframe = 0
    for g0, g1, clo in dom["segs"]:
        for G0, G1, Clo in target["cds_frame"]:
            lo, hi = max(g0, G0), min(g1, G1)
            if hi < lo:
                continue
            n = hi - lo + 1
            covered += n
            delta = _cds_at(lo, G0, G1, Clo, strand) - _cds_at(lo, g0, g1, clo, strand)
            if delta % 3 == 0:
                inframe += n
    return domain_bp, covered, inframe


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--isoforms", required=True, type=Path)
    ap.add_argument("--genes", required=True, type=Path)
    ap.add_argument("--pfam-meta", required=True, type=Path)
    ap.add_argument("--index", required=True, type=str)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--genes-limit", type=int, default=0,
                    help="For testing: restrict to the first N genes")
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args()

    genes = pd.read_csv(args.genes, sep="\t")
    isos = pd.read_csv(args.isoforms, sep="\t", dtype={"chrom": str})
    meta = pd.read_csv(args.pfam_meta, sep="\t")

    if args.genes_limit:
        keep = set(genes.gene_id.head(args.genes_limit))
        genes = genes[genes.gene_id.isin(keep)]
        isos = isos[isos.gene_id.isin(keep)]

    source_of = dict(zip(genes.gene_id, genes.source_protein_id))
    name_of = dict(zip(genes.gene_id, genes.gene_name))
    ensp2gene = dict(zip(isos.protein_id, isos.gene_id))
    src_set = set(genes.source_protein_id)

    # ---- source domains (Pfam on source isoforms) ----
    sd = meta[meta.protein_id.isin(src_set)].copy()
    qid2pfam = dict(zip(sd.query_id, sd.pfam_id))
    qid2ensp = dict(zip(sd.query_id, sd.protein_id))
    src_q = [dict(id=r.protein_id, aa_start=int(r.aa_start), aa_end=int(r.aa_end),
                  domain_id=r.query_id) for r in sd.itertuples()]
    print(f"source domains to project: {len(src_q):,}", file=sys.stderr)

    mp = fc.Mapper(args.index, threads=args.threads)
    src_res = mp.map_batch(src_q, output="coding")
    seg = src_res.cds_segments
    seg = seg[seg.overlaps_domain == "coding_overlap"]
    dom_by_gene = {}
    for qid, d in _segments(seg).items():
        g = ensp2gene.get(qid2ensp.get(qid))
        if g is None:
            continue
        dom_by_gene.setdefault(g, []).append((qid, d))

    # ---- target isoforms (every non-source isoform), structure-only ----
    tgt = isos[isos.protein_id != isos.gene_id.map(source_of)]
    tgt_q = [dict(id=e, domain_id=e) for e in tgt.protein_id.unique()]
    print(f"target isoforms to map: {len(tgt_q):,}", file=sys.stderr)
    tgt_res = mp.map_batch(tgt_q, output="isoform")
    struct = _structures(tgt_res.isoform)
    tgt_by_gene = {}
    for gene_id, ensp in tgt[["gene_id", "protein_id"]].itertuples(index=False):
        if ensp in struct:
            tgt_by_gene.setdefault(gene_id, []).append(ensp)

    # ---- intersect domain x target within each gene ----
    rows = []
    for g, doms in dom_by_gene.items():
        targets = tgt_by_gene.get(g, [])
        if not targets:
            continue
        for qid, dom in doms:
            for ensp in targets:
                tgt_struct = struct[ensp]
                domain_bp, covered, inframe = _score(dom, tgt_struct)
                if domain_bp == 0:
                    continue
                cov = covered / domain_bp
                ifc = inframe / domain_bp
                mech = cm.classify_mechanism(dom["intervals"], tgt_struct, cov, ifc)
                rows.append((qid, qid2pfam.get(qid, ""), g, name_of.get(g, ""),
                             source_of.get(g, ""), ensp, domain_bp, covered, inframe,
                             round(cov, 4), round(ifc, 4),
                             cm.conservation_class(ifc), mech))

    cols = ["query_id", "pfam_id", "gene_id", "gene_name", "source_ensp",
            "target_ensp", "domain_bp", "covered_bp", "inframe_bp",
            "coverage", "inframe_coverage", "conservation_class", "mechanism"]
    out = pd.DataFrame(rows, columns=cols)
    out.to_csv(args.out, sep="\t", index=False)
    print(f"wrote {len(out):,} (domain x target) rows -> {args.out}")
    if len(out):
        print(out["conservation_class"].value_counts().to_string())


if __name__ == "__main__":
    main()
