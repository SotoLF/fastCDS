#!/usr/bin/env python3
"""Generate GTF rows for custom protein/transcript entries.

Use case: proteins you care about aren't in the reference annotation —
transgenes, non-reference ORFs, manually curated alternative isoforms,
novel CDSs from a focused study. To map domains onto them with prot2exon
you need them in the GTF before you build the index (`prot2exon index`).

Input is a TSV with one row per transcript:

    protein_id      transcript_id   gene_id    gene_name   chrom    strand   blocks                            source
    NP_NOVEL_1      NM_NOVEL_1      G_NOV1     NOVEL_1     chr_X    +        2309-3221                          custom
    NP_NOVEL_2      NM_NOVEL_2      G_NOV2     NOVEL_2     chr_X    +        102449-102544;103176-103508        custom
    custom_orf_1    custom_orf_1_tx custom_orf_1_g custom_orf_1 chrSYN +    1-300                                manual

`blocks` is a semicolon-separated list of `start-end` pairs (1-based, inclusive,
in genomic order — the script handles strand-aware exon numbering itself).

The script emits one `transcript`, N `exon`, and N `CDS` lines per row, all
with the attribute keys prot2exon's parser reads (gene_id, transcript_id,
protein_id, gene_name, exon_number, plus a `source` tag for provenance).

Append the output to your existing GTF and rebuild the index:

    python3 scripts/append_custom_proteins.py --in custom.tsv \
        >> /path/to/working_copy.gtf
    ./build/prot2exon index --gtf /path/to/working_copy.gtf \
        --out human_plus_custom.idx
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


REQUIRED_COLS = ["protein_id", "transcript_id", "gene_id", "gene_name",
                 "chrom", "strand", "blocks"]


def parse_blocks(spec: str) -> list[tuple[int, int]]:
    """`100-149;200-250` → [(100,149), (200,250)] sorted ascending by start."""
    out = []
    for chunk in spec.split(";"):
        chunk = chunk.strip()
        if not chunk:
            continue
        if "-" not in chunk:
            raise ValueError(f"bad block (need start-end): {chunk!r}")
        s_str, e_str = chunk.split("-", 1)
        s, e = int(s_str), int(e_str)
        if s > e:
            raise ValueError(f"start > end in block {chunk!r}")
        out.append((s, e))
    if not out:
        raise ValueError("no blocks parsed")
    return sorted(out)


def fmt_attrs(pairs: list[tuple[str, str]]) -> str:
    parts = [f'{k} "{v}"' for k, v in pairs if v is not None and v != ""]
    return "; ".join(parts) + ";"


def emit_rows(row: dict, source_tag: str) -> list[str]:
    pid = row["protein_id"]
    tid = row["transcript_id"]
    gid = row["gene_id"]
    gname = row["gene_name"]
    chrom = row["chrom"]
    strand = row["strand"]
    if strand not in ("+", "-"):
        raise ValueError(f"{tid}: strand must be + or -, got {strand!r}")
    blocks = parse_blocks(row["blocks"])
    src = row.get("source") or source_tag

    tx_start = min(s for s, _ in blocks)
    tx_end = max(e for _, e in blocks)

    # Exon numbering follows GENCODE: 1..N in translation order.
    if strand == "-":
        numbered = list(zip(reversed(blocks),
                            range(1, len(blocks) + 1)))
    else:
        numbered = list(zip(blocks, range(1, len(blocks) + 1)))
    exon_num_by_block = {b: n for b, n in numbered}

    lines = []
    base = [("gene_id", gid), ("transcript_id", tid),
            ("gene_name", gname), ("protein_id", pid)]

    # transcript row
    lines.append("\t".join([
        chrom, src, "transcript", str(tx_start), str(tx_end),
        ".", strand, ".",
        fmt_attrs(base + [("transcript_type", "protein_coding")])
    ]))

    # exon rows
    for (s, e) in blocks:
        n = exon_num_by_block[(s, e)]
        lines.append("\t".join([
            chrom, src, "exon", str(s), str(e),
            ".", strand, ".",
            fmt_attrs(base + [("exon_number", str(n)),
                              ("exon_id", f"{tid}_E{n}")])
        ]))

    # CDS rows (same blocks; we assume the whole transcript codes — typical
    # for non-canonical / synthetic ORFs. If you need 5'/3' UTR padding, just expand
    # the exon block past the CDS block in the input).
    for (s, e) in blocks:
        n = exon_num_by_block[(s, e)]
        lines.append("\t".join([
            chrom, src, "CDS", str(s), str(e),
            ".", strand, "0",
            fmt_attrs(base + [("exon_number", str(n))])
        ]))

    return lines


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--in", dest="infile", required=True,
                    help="TSV with columns: " + ", ".join(REQUIRED_COLS) +
                         " [, source]")
    ap.add_argument("--source-tag", default="custom",
                    help="Default value for the GTF source column when the "
                         "input row has no `source` column (default: custom).")
    args = ap.parse_args()

    src_path = Path(args.infile)
    if not src_path.exists():
        sys.exit(f"error: input not found: {src_path}")

    # The GTF rows go to stdout (pipe/redirect them onto your reference GTF);
    # status and per-row errors go to stderr so they don't pollute the output.
    with open(src_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing = [c for c in REQUIRED_COLS if c not in (reader.fieldnames or [])]
        if missing:
            sys.exit(f"error: input missing columns: {missing}\n"
                     f"got: {reader.fieldnames}")
        n_rows = 0
        n_errors = 0
        for i, row in enumerate(reader, start=2):  # data rows start at line 2
            try:
                for line in emit_rows(row, args.source_tag):
                    sys.stdout.write(line + "\n")
                n_rows += 1
            except Exception as e:
                n_errors += 1
                sys.stderr.write(f"line {i}: {e}\n")
    sys.stderr.write(f"wrote {n_rows} transcripts ({n_errors} errors)\n")
    return 1 if n_errors else 0


if __name__ == "__main__":
    sys.exit(main())
