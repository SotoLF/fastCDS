"""Sample a stratified validation set of (protein_id, aa_start, aa_end) queries
from a GENCODE/Ensembl GTF, per PLAN.txt Phase 2.

The 9 strata (totaling 5,000 by default):
    1,000 single-exon DOMAINS      — aa range fits inside one CDS exon
    1,000 multi-exon DOMAINS       — aa range spans >= 2 CDS exons
      500 codon-split-boundary DOMAINS — aa range covers a codon that straddles an exon boundary
    1,000 plus-strand GENES        — random aa range on a + strand transcript
    1,000 minus-strand GENES       — random aa range on a - strand transcript
      200 CDS-incomplete proteins  — cds_start_NF / cds_end_NF tagged
      100 selenoproteins           — curated gene-name list
      100 single-exon GENES        — transcripts with 1 CDS exon
      100 many-exon GENES          — transcripts with > 20 CDS exons

Outputs:
  queries.bed       — 5,000 lines: protein_id\\taa_start\\taa_end\\tQUERY_<i>\\t<category>
  queries_meta.tsv  — one row per query: query_id, category, protein_id, transcript_id,
                      strand, n_cds_exons, gene_name, aa_start, aa_end

The metadata TSV is what the validation script uses for stratified reporting
in Table 1.
"""

from __future__ import annotations

import argparse
import gzip
import random
import re
import sys
from pathlib import Path

# Canonical human selenoprotein gene-symbol list (25 known).
SELENOPROTEINS = {
    "SELENOP", "SELENOO", "SELENOI", "SELENOH", "SELENOK", "SELENOM",
    "SELENON", "SELENOS", "SELENOT", "SELENOV", "SELENOW", "SELENOF",
    "GPX1", "GPX2", "GPX3", "GPX4",
    "TXNRD1", "TXNRD2", "TXNRD3",
    "DIO1", "DIO2", "DIO3",
    "MSRB1", "SEPHS2",
    # Two former names that sometimes still appear:
    "SEPP1",  # legacy for SELENOP
}

ATTR_RE = re.compile(r'(\w+)\s+"([^"]*)"')


def parse_attrs(field: str) -> dict[str, str]:
    return {k: v for k, v in ATTR_RE.findall(field)}


def parse_gtf(gtf_path: Path):
    """Yield transcript records with CDS layouts.

    Returns dicts:
        transcript_id, protein_id, gene_name, strand, chrom,
        cds_lengths (in translation order — ascending genomic on '+', reversed on '-'),
        n_cds_exons, has_cds_NF (True if cds_start_NF or cds_end_NF tagged),
        protein_length_aa (sum(cds_lengths) // 3),
        codon_split_aas (set of aa positions where a codon straddles an exon boundary)
    """
    # Aggregate per transcript: collect CDS rows; emit when we move on.
    current_tx: dict[str, list] = {}
    transcript_meta: dict[str, dict] = {}

    opener = gzip.open if str(gtf_path).endswith(".gz") else open

    with opener(gtf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            feat = parts[2]
            if feat not in ("CDS", "transcript"):
                continue
            attrs = parse_attrs(parts[8])
            tid = attrs.get("transcript_id")
            if not tid:
                continue
            # Strip versions for cross-tool compatibility (prot2exon strips, ensembldb does not).
            tid_unv = tid.split(".")[0]

            if feat == "transcript":
                # Tag list — GENCODE writes multiple `tag "X"` entries; ATTR_RE keeps the last.
                tags = re.findall(r'tag "([^"]+)"', parts[8])
                transcript_meta[tid_unv] = {
                    "chrom": parts[0],
                    "strand": parts[6],
                    "gene_name": attrs.get("gene_name", ""),
                    "protein_id": attrs.get("protein_id", "").split(".")[0],  # GENCODE puts it here; Ensembl doesn't
                    "has_cds_NF": ("cds_start_NF" in tags) or ("cds_end_NF" in tags),
                    # GENCODE uses transcript_type; Ensembl uses transcript_biotype.
                    "transcript_type": attrs.get("transcript_type") or attrs.get("transcript_biotype", ""),
                }
            else:  # CDS
                current_tx.setdefault(tid_unv, []).append(
                    (int(parts[3]), int(parts[4]), parts[6])
                )
                # Ensembl puts protein_id only on CDS rows — backfill if the transcript
                # record didn't already carry it.
                meta = transcript_meta.get(tid_unv)
                if meta is not None and not meta["protein_id"]:
                    pid = attrs.get("protein_id", "").split(".")[0]
                    if pid:
                        meta["protein_id"] = pid

    # Build the per-transcript records.
    for tid, cds_list in current_tx.items():
        meta = transcript_meta.get(tid)
        if not meta:
            continue
        if meta["transcript_type"] != "protein_coding":
            continue
        if not meta["protein_id"]:
            continue

        strand = meta["strand"]
        # CDS rows in translation order: ascending on '+', descending on '-'.
        sorted_genomic = sorted(cds_list, key=lambda r: r[0])
        ordered = sorted_genomic if strand == "+" else list(reversed(sorted_genomic))
        cds_lengths = [(e - s + 1) for (s, e, _) in ordered]

        total_nt = sum(cds_lengths)
        if total_nt < 9:  # need at least 3 codons to be useful
            continue
        protein_length_aa = total_nt // 3

        # Codon-split aas: at each internal cumulative-nt boundary not divisible
        # by 3, the codon straddling that boundary is split.
        codon_split_aas: set[int] = set()
        cum = 0
        for L in cds_lengths[:-1]:
            cum += L
            if cum % 3 != 0:
                # The split codon's aa number is ceil(cum / 3).
                codon_split_aas.add((cum + 2) // 3)

        yield {
            "transcript_id": tid,
            "protein_id": meta["protein_id"],
            "gene_name": meta["gene_name"],
            "strand": strand,
            "chrom": meta["chrom"],
            "n_cds_exons": len(cds_lengths),
            "cds_lengths": cds_lengths,
            "has_cds_NF": meta["has_cds_NF"],
            "protein_length_aa": protein_length_aa,
            "codon_split_aas": codon_split_aas,
        }


def pick_single_exon_domain_range(tx, rng: random.Random) -> tuple[int, int] | None:
    """Pick an aa range entirely contained within one CDS exon (in translation order)."""
    eligible = [(i, L) for i, L in enumerate(tx["cds_lengths"]) if L >= 6]  # need >= 2 aa to pick range
    if not eligible:
        return None
    i, L = rng.choice(eligible)
    # aa range encoded by exon i: depends on cumulative nt before it.
    cum_before = sum(tx["cds_lengths"][:i])
    # aa fully inside exon i: aa with both first and last nt inside the exon's nt span.
    # First aa fully inside: ceil((cum_before + 1) / 3) if (cum_before % 3) == 0, else higher.
    # Simpler: enumerate candidate aas and keep ones whose nt range [3*(aa-1)+1, 3*aa] is inside the exon's nt range.
    exon_nt_start = cum_before + 1
    exon_nt_end = cum_before + L
    first_aa = (exon_nt_start + 2) // 3  # may start mid-codon — bump to next aa fully inside
    while first_aa <= tx["protein_length_aa"] and (first_aa - 1) * 3 + 1 < exon_nt_start:
        first_aa += 1
    last_aa = exon_nt_end // 3
    while last_aa > 0 and last_aa * 3 > exon_nt_end:
        last_aa -= 1
    if last_aa - first_aa < 1:
        return None
    aa_start = rng.randint(first_aa, last_aa - 1)
    aa_end = rng.randint(aa_start, last_aa)
    return aa_start, aa_end


def pick_multi_exon_domain_range(tx, rng: random.Random) -> tuple[int, int] | None:
    """Pick an aa range spanning >= 2 CDS exons."""
    if tx["n_cds_exons"] < 2:
        return None
    # Pick two distinct exons in translation order; build a range that spans them.
    i, j = sorted(rng.sample(range(tx["n_cds_exons"]), 2))
    cum_before_i = sum(tx["cds_lengths"][:i])
    cum_after_j = sum(tx["cds_lengths"][:j + 1])
    # Take an aa near the start of exon i and an aa near the end of exon j.
    aa_start = max(1, (cum_before_i + 1 + 2) // 3)
    aa_end = min(tx["protein_length_aa"], cum_after_j // 3)
    if aa_end <= aa_start:
        return None
    return aa_start, aa_end


def pick_codon_split_range(tx, rng: random.Random) -> tuple[int, int] | None:
    """Pick an aa range that includes at least one codon-split aa."""
    if not tx["codon_split_aas"]:
        return None
    split_aa = rng.choice(list(tx["codon_split_aas"]))
    # Embed the split codon in a small range around it.
    half = rng.randint(2, 10)
    aa_start = max(1, split_aa - half)
    aa_end = min(tx["protein_length_aa"], split_aa + half)
    return aa_start, aa_end


def pick_random_range(tx, rng: random.Random) -> tuple[int, int] | None:
    """Pick any aa range within the protein."""
    if tx["protein_length_aa"] < 5:
        return None
    aa_start = rng.randint(1, tx["protein_length_aa"] - 1)
    aa_end = rng.randint(aa_start, tx["protein_length_aa"])
    return aa_start, aa_end


STRATA_SPEC = [
    # (category_name, target_count, transcript_filter, range_picker)
    ("single_exon_domain", 1000, lambda tx: True, pick_single_exon_domain_range),
    ("multi_exon_domain", 1000, lambda tx: tx["n_cds_exons"] >= 2, pick_multi_exon_domain_range),
    ("codon_split_boundary", 500, lambda tx: len(tx["codon_split_aas"]) > 0, pick_codon_split_range),
    ("plus_strand_gene", 1000, lambda tx: tx["strand"] == "+", pick_random_range),
    ("minus_strand_gene", 1000, lambda tx: tx["strand"] == "-", pick_random_range),
    ("cds_incomplete", 200, lambda tx: tx["has_cds_NF"], pick_random_range),
    ("selenoprotein", 100, lambda tx: tx["gene_name"] in SELENOPROTEINS, pick_random_range),
    ("single_exon_gene", 100, lambda tx: tx["n_cds_exons"] == 1, pick_random_range),
    ("many_exon_gene", 100, lambda tx: tx["n_cds_exons"] > 20, pick_random_range),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gtf", required=True, type=Path)
    ap.add_argument("--out-bed", required=True, type=Path)
    ap.add_argument("--out-meta", required=True, type=Path)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    print(f"parsing {args.gtf} ...", file=sys.stderr)
    transcripts = list(parse_gtf(args.gtf))
    print(f"  {len(transcripts):,} protein-coding transcripts", file=sys.stderr)

    rng = random.Random(args.seed)

    bed_rows = []
    meta_rows = []
    qid = 0
    for cat, target, tx_filter, range_picker in STRATA_SPEC:
        pool = [tx for tx in transcripts if tx_filter(tx)]
        if not pool:
            print(f"  WARN: no transcripts match category {cat}", file=sys.stderr)
            continue
        attempts = 0
        accepted = 0
        # Sample with replacement, picking aa range freshly each time.
        # Reject if the range picker returns None for this transcript.
        while accepted < target:
            attempts += 1
            if attempts > target * 50:
                print(f"  WARN: only sampled {accepted}/{target} for {cat} "
                      f"after {attempts} attempts", file=sys.stderr)
                break
            tx = rng.choice(pool)
            rng_pick = range_picker(tx, rng)
            if rng_pick is None:
                continue
            aa_start, aa_end = rng_pick
            qid += 1
            query_id = f"Q{qid:06d}_{cat}"
            bed_rows.append(f"{tx['protein_id']}\t{aa_start}\t{aa_end}\t{query_id}\n")
            meta_rows.append((query_id, cat, tx["protein_id"], tx["transcript_id"],
                              tx["strand"], tx["n_cds_exons"], tx["gene_name"],
                              aa_start, aa_end))
            accepted += 1
        print(f"  {cat}: {accepted}/{target} (pool of {len(pool)} transcripts)", file=sys.stderr)

    args.out_bed.write_text("".join(bed_rows))
    with open(args.out_meta, "w") as f:
        f.write("query_id\tcategory\tprotein_id\ttranscript_id\tstrand\t"
                "n_cds_exons\tgene_name\taa_start\taa_end\n")
        for row in meta_rows:
            f.write("\t".join(str(v) for v in row) + "\n")
    print(f"wrote {len(bed_rows):,} queries -> {args.out_bed}", file=sys.stderr)
    print(f"wrote metadata -> {args.out_meta}", file=sys.stderr)


if __name__ == "__main__":
    main()
