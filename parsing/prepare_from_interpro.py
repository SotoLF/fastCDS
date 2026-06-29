#!/usr/bin/env python3
"""Convert an InterProScan TSV into the BED-like format fastCDS eats.

InterProScan `-f TSV` produces one row per signature hit, with columns:

    1  protein_accession  (whatever was in the input FASTA: UniProt acc,
                            ENSP, your own ID, ...)
    2  md5
    3  length
    4  analysis             (Pfam, SMART, PANTHER, SUPERFAMILY, ...)
    5  signature_accession  (e.g. PF00069, SM00220, ...)
    6  signature_description
    7  start                (1-based aa)
    8  stop                 (1-based aa, inclusive)
    9  score / e-value
    10 status               (T)
    11 date
    12 interpro_accession   (e.g. IPR000719)  — may be "-"
    13 interpro_description                    — may be "-"
    14 GO annotations (optional)
    15 pathways      (optional)

We pass columns 1, 7, 8 straight through (protein → aa_start → aa_end). The
domain_id column we synthesize from a combination of the InterPro accession
(or signature accession as fallback) and the analysis source, so it stays
unique even when the same protein has overlapping Pfam + SMART hits for the
same domain. The original description goes into the 5th BED column for
traceability.

Mapping UniProt → ENSP
----------------------

If the input protein IDs in column 1 are already ENSPs, nothing to do — pass
`--id-type ensp` (or rely on auto-detection).

Otherwise the parser needs a mapping. See `scripts/_mapping.py` for the
Ensembl xref TSV format. Usage:

    python3 scripts/prepare_from_interpro.py \\
        --in proteome.interpro.tsv \\
        --mapping Homo_sapiens.GRCh38.110.uniprot.tsv \\
        --out domains.bed

A protein that doesn't resolve is dropped and reported on stderr. Use
`--keep-unmapped` to write a side TSV of the rejected rows for later review.
"""

from __future__ import annotations

import argparse
import os
import re
import sys

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(THIS_DIR, "..", "python"))
from fastCDS.prepare._mapping import (
    UniProtToEnsp, looks_like_uniprot, looks_like_ensp,
    open_text, dedup_and_sort_rows, write_bed_row,
)
from fastCDS.prepare._interpro import parse, detect_id_type


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description="InterProScan TSV → fastCDS BED-like.")
    p.add_argument("--in", dest="in_path", required=True,
                   help="InterProScan TSV (column 1 must be the protein id).")
    p.add_argument("--out", default="-",
                   help="Output BED-like (default stdout). Existing files are overwritten.")
    p.add_argument("--mapping",
                   help="Ensembl UniProt xref TSV (only required when input IDs are UniProt).")
    p.add_argument("--simple-mapping",
                   help="Two-column TSV mapping (uniprot\\tensp). Alternative to --mapping.")
    p.add_argument("--id-type", choices=("auto", "ensp", "uniprot"), default="auto",
                   help="Type of identifier in column 1.")
    p.add_argument("--analyses",
                   help="Comma-separated list of analyses to keep (default: keep all). "
                        "Common: Pfam,SMART,PROSITE_PROFILES,SUPERFAMILY,Gene3D.")
    p.add_argument("--min-length", type=int, default=5,
                   help="Drop hits shorter than this many aa.")
    p.add_argument("--keep-unmapped",
                   help="Write rejected rows to this file (with reason in column 2).")
    args = p.parse_args(argv)

    mapping = None
    if args.mapping:
        mapping = UniProtToEnsp.from_ensembl_xref_tsv(args.mapping)
    elif args.simple_mapping:
        mapping = UniProtToEnsp.from_simple_tsv(args.simple_mapping)

    analyses = set(args.analyses.split(",")) if args.analyses else None

    rows, rejected = parse(args.in_path,
                           analyses=analyses,
                           min_length=args.min_length,
                           id_type=args.id_type,
                           mapping=mapping,
                           source_filter=None)
    rows = dedup_and_sort_rows(rows)

    out_fh = sys.stdout if args.out == "-" else open(args.out, "w")
    try:
        out_fh.write("# fastCDS BED-like, generated from InterProScan\n")
        out_fh.write("# columns: ENSP\\taa_start\\taa_end\\tdomain_id\\tsignature_description\n")
        for r in rows:
            write_bed_row(out_fh, ensp=r[0], aa_start=r[1], aa_end=r[2],
                          domain_id=r[3], source=r[4])
    finally:
        if out_fh is not sys.stdout:
            out_fh.close()
    print(f"[interpro] wrote {len(rows)} rows "
          f"({len(rejected)} rejected) to {args.out}", file=sys.stderr)
    if args.keep_unmapped and rejected:
        with open(args.keep_unmapped, "w") as f:
            f.write("raw_row\treason\n")
            for raw, why in rejected:
                f.write(f"{raw}\t{why}\n")
        print(f"[interpro] rejected rows in {args.keep_unmapped}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
