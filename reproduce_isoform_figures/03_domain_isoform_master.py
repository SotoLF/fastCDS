#!/usr/bin/env python3
"""Step 3 - map each Pfam domain across a gene's isoforms and score it.

Every Pfam domain on a gene's reference isoform is placed on the genome with
fastCDS, then checked in each other isoform of that gene:

  coverage           = domain coding bases still present in the isoform
  inframe_coverage   = those bases that also keep the reading frame
  conservation_class = conserved (inframe >= 0.80) / partial (0.50-0.80) / lost

Outputs
  domain_isoform_master.tsv     one row per (domain, isoform)
  domain_genomic_intervals.tsv  one row per domain coding-exon segment

Needs the fastCDS package + a prebuilt index, the Pfam hits table, and
Pfam-A.clans (for clan names).

Run
  python 03_domain_isoform_master.py \
      --isoforms isoforms.tsv --genes genes.tsv \
      --pfam-meta pfam_human_v115_meta.tsv --clans Pfam-A.clans.tsv \
      --index ensembl_v115_human.idx --threads 4 \
      --out-master domain_isoform_master.tsv \
      --out-intervals domain_genomic_intervals.tsv
"""
import argparse
from pathlib import Path

import pandas as pd
import fastCDS as fc

CONSERVED_T, PARTIAL_T = 0.80, 0.50   # thresholds on in-frame coverage


def conservation_class(inframe_coverage):
    if inframe_coverage >= CONSERVED_T:
        return "conserved"
    if inframe_coverage >= PARTIAL_T:
        return "partial"
    return "lost"


def cds_at(g, G0, G1, Clo, strand):
    """CDS-nt coordinate of genomic base g inside CDS exon [G0, G1]."""
    return Clo + (g - G0) if strand == "+" else Clo + (G1 - g)


def score(segs, intervals, strand, target_cds):
    """Return (domain_bp, covered_bp, inframe_bp) for one domain vs one isoform.

    Over any overlap the CDS-nt shift (target - source) is constant, so one
    `shift % 3 == 0` test tells whether the whole overlap stays in frame.
    """
    domain_bp = sum(g1 - g0 + 1 for g0, g1 in intervals)
    covered = inframe = 0
    for g0, g1, clo in segs:
        for G0, G1, Clo in target_cds:
            lo, hi = max(g0, G0), min(g1, G1)
            if hi < lo:
                continue
            n = hi - lo + 1
            covered += n
            if (cds_at(lo, G0, G1, Clo, strand) - cds_at(lo, g0, g1, clo, strand)) % 3 == 0:
                inframe += n
    return domain_bp, covered, inframe


def load_clans(path):
    """Pfam-A.clans.tsv: pfam_acc, clan_acc, clan_name, pfam_name, description."""
    c = pd.read_csv(path, sep="\t", header=None,
                    names=["pfam_id", "clan_acc", "clan_name", "pfam_name", "desc"])
    return (dict(zip(c.pfam_id, c.clan_acc.fillna(""))),
            dict(zip(c.pfam_id, c.clan_name.fillna(""))))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--isoforms", required=True, type=Path)
    ap.add_argument("--genes", required=True, type=Path)
    ap.add_argument("--pfam-meta", required=True, type=Path)
    ap.add_argument("--clans", required=True, type=Path)
    ap.add_argument("--index", required=True)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--out-master", required=True, type=Path)
    ap.add_argument("--out-intervals", required=True, type=Path)
    args = ap.parse_args()

    genes = pd.read_csv(args.genes, sep="\t")
    isos = pd.read_csv(args.isoforms, sep="\t", dtype={"chrom": str})
    meta = pd.read_csv(args.pfam_meta, sep="\t")
    clan_acc, clan_name = load_clans(args.clans)

    source_of = dict(zip(genes.gene_id, genes.source_protein_id))
    name_of = dict(zip(genes.gene_id, genes.gene_name))
    ensp2gene = dict(zip(isos.protein_id, isos.gene_id))
    ref_set = set(genes.source_protein_id)

    # ---- 1. Pfam domains on the reference isoforms ----
    sd = meta[meta.protein_id.isin(ref_set)].copy()
    qid2pfam = dict(zip(sd.query_id, sd.pfam_id))
    qid2ensp = dict(zip(sd.query_id, sd.protein_id))
    src_q = [{"protein_id": r.protein_id, "aa_start": int(r.aa_start),
              "aa_end": int(r.aa_end), "domain_id": r.query_id}
             for r in sd.itertuples()]
    print(f"domains to project: {len(src_q):,}")

    mp = fc.Mapper(args.index, threads=args.threads)

    # ---- 2. place each domain on the genome (its CDS segments) ----
    seg = mp.map_batch(src_q, output="coding").cds_segments
    seg = seg[seg.overlaps_domain == "coding_overlap"].copy()

    gi = pd.DataFrame({
        "domain_id": seg.domain_id,
        "source_protein_id": seg.protein_id,
        "chrom": seg.chrom, "strand": seg.strand,
        "exon_number": seg.exon_number,
        "genomic_start": seg.domain_overlap_genomic_start.astype(int),
        "genomic_end": seg.domain_overlap_genomic_end.astype(int),
    })
    gi["seg_bp"] = gi.genomic_end - gi.genomic_start + 1
    gi.to_csv(args.out_intervals, sep="\t", index=False)
    print(f"wrote {len(gi):,} coding segments -> {args.out_intervals}")

    # per-domain genome geometry, grouped by gene, ready for scoring
    dom_by_gene = {}
    for did, g in seg.groupby("domain_id"):
        segs = list(zip(g.domain_overlap_genomic_start.astype(int),
                        g.domain_overlap_genomic_end.astype(int),
                        g.domain_overlap_cds_nt_start.astype(int)))
        d = dict(strand=g.strand.iloc[0], segs=segs,
                 intervals=[(s[0], s[1]) for s in segs])
        gene = ensp2gene.get(qid2ensp.get(did))
        if gene:
            dom_by_gene.setdefault(gene, []).append((did, d))

    # ---- 3. read the CDS structure of every non-reference isoform ----
    tgt = isos[isos.protein_id != isos.gene_id.map(source_of)]
    tgt_q = [{"protein_id": e, "domain_id": e} for e in tgt.protein_id.unique()]
    print(f"isoforms to read: {len(tgt_q):,}")
    iso = mp.map_batch(tgt_q, output="isoform").isoform
    cds = iso[iso.feature_type == "CDS"]
    tgt_cds = {}
    for ensp, g in cds.groupby("protein_id"):
        tgt_cds[ensp] = list(zip(g.feature_genomic_start.astype(int),
                                 g.feature_genomic_end.astype(int),
                                 g.cds_nt_start.astype(int)))
    tgt_by_gene = {}
    for gene, ensp in tgt[["gene_id", "protein_id"]].itertuples(index=False):
        if ensp in tgt_cds:
            tgt_by_gene.setdefault(gene, []).append(ensp)

    # ---- 4. score each domain against each isoform of its gene ----
    rows = []
    for gene, doms in dom_by_gene.items():
        for ensp in tgt_by_gene.get(gene, []):
            target = tgt_cds[ensp]
            for did, d in doms:
                domain_bp, covered, inframe = score(d["segs"], d["intervals"],
                                                    d["strand"], target)
                if domain_bp == 0:
                    continue
                cov, ifc = covered / domain_bp, inframe / domain_bp
                pf = qid2pfam.get(did, "")
                rows.append((did, pf, clan_acc.get(pf, ""), clan_name.get(pf, ""),
                             gene, name_of.get(gene, ""), source_of.get(gene, ""), ensp,
                             domain_bp, covered, inframe, round(cov, 4), round(ifc, 4),
                             conservation_class(ifc)))

    cols = ["domain_id", "pfam_id", "clan", "clan_name", "gene_id", "gene_name",
            "source_ensp", "isoform_id", "domain_bp", "covered_bp", "inframe_bp",
            "coverage", "inframe_coverage", "conservation_class"]
    master = pd.DataFrame(rows, columns=cols)
    master.to_csv(args.out_master, sep="\t", index=False)
    print(f"wrote {len(master):,} (domain, isoform) rows -> {args.out_master}")
    print(master.conservation_class.value_counts().to_string())


if __name__ == "__main__":
    main()
