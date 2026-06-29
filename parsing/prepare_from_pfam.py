#!/usr/bin/env python3
"""Convert HMMER domtblout (Pfam scan) to the BED-like format fastCDS eats.

HMMER's `--domtblout` is the right file: one row per *domain hit* (not per
sequence hit), so we get per-domain coordinates. Format reference:
http://eddylab.org/software/hmmer/Userguide.pdf §10.

Important columns (space-separated, ragged):

    1   target name           (depends on direction, see below)
    2   target accession
    3   tlen                  (target length)
    4   query name
    5   query accession
    6   qlen                  (query length)
    7   E-value (full)        full-sequence E-value
    8   score (full)          full-sequence bit score
    9   bias (full)
    10  #                     this-domain index
    11  of                    total domains for this protein
    12  c-Evalue              conditional E-value of this domain
    13  i-Evalue              independent E-value of this domain
    14  score (this)          this-domain bit score
    15  bias (this)
    16  hmm from
    17  hmm to
    18  ali from              ← we use these
    19  ali to                ← we use these
    20  env from
    21  env to
    22  acc
    23+ description (rest of the line)

Direction
---------

If you ran `hmmscan` with the Pfam HMM database as the database (typical):
    hmmscan --domtblout out.dom Pfam-A.hmm proteins.fa
then **target = HMM (Pfam family), query = your protein**. Pass `--mode scan`
(this is the default).

If you ran `hmmsearch` with your protein FASTA as the database:
    hmmsearch --domtblout out.dom one_hmm.hmm proteins.fa
then **target = your protein, query = HMM**. Pass `--mode search`.

We compute the protein-coordinate range from `ali from / ali to` of the
**protein** column (not the HMM). hmm_from / hmm_to are HMM-coordinate
positions and would be wrong here.

ID mapping
----------

If your protein FASTA used ENSPs already (the GENCODE-derived FASTA is the
typical case), no mapping is needed — `--id-type ensp`. If you used UniProt
accessions, supply `--mapping` (Ensembl xref TSV) or `--simple-mapping`.
"""

from __future__ import annotations

import argparse
import os
import re
import sys

# Make the package importable when running directly from a repo checkout.
THIS_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(THIS_DIR, "..", "python"))
from fastCDS.prepare._mapping import (
    UniProtToEnsp, looks_like_ensp, looks_like_uniprot,
    open_text, dedup_and_sort_rows, write_bed_row, strip_version,
)
from fastCDS.prepare._pfam import parse


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description="HMMER domtblout (Pfam scan) → fastCDS BED-like.")
    p.add_argument("--in", dest="in_path", required=True,
                   help="HMMER --domtblout file.")
    p.add_argument("--out", default="-")
    p.add_argument("--mode", choices=("scan", "search"), default="scan",
                   help="hmmscan: target=HMM, query=protein (default). "
                        "hmmsearch: target=protein, query=HMM.")
    p.add_argument("--min-score", type=float, default=0.0,
                   help="Drop hits below this per-domain bit score.")
    p.add_argument("--min-length", type=int, default=5)
    p.add_argument("--id-type", choices=("auto", "ensp", "uniprot"), default="auto")
    p.add_argument("--mapping",
                   help="Ensembl UniProt xref TSV (only required when protein IDs are UniProt).")
    p.add_argument("--simple-mapping",
                   help="Two-column TSV mapping (uniprot\\tensp). Alternative.")
    p.add_argument("--keep-unmapped",
                   help="Write rejected rows to this file.")
    args = p.parse_args(argv)

    mapping = None
    if args.mapping:
        mapping = UniProtToEnsp.from_ensembl_xref_tsv(args.mapping)
    elif args.simple_mapping:
        mapping = UniProtToEnsp.from_simple_tsv(args.simple_mapping)

    rows, rejected = parse(args.in_path,
                           mode=args.mode,
                           min_score=args.min_score,
                           min_length=args.min_length,
                           id_type=args.id_type,
                           mapping=mapping)
    rows = dedup_and_sort_rows(rows)

    out_fh = sys.stdout if args.out == "-" else open(args.out, "w")
    try:
        out_fh.write("# fastCDS BED-like, generated from HMMER domtblout (Pfam)\n")
        out_fh.write("# columns: ENSP\\taa_start\\taa_end\\tdomain_id\\tdescription\n")
        for r in rows:
            write_bed_row(out_fh, ensp=r[0], aa_start=r[1], aa_end=r[2],
                          domain_id=r[3], source=r[4])
    finally:
        if out_fh is not sys.stdout:
            out_fh.close()
    print(f"[pfam] wrote {len(rows)} rows "
          f"({len(rejected)} rejected) to {args.out}", file=sys.stderr)
    if args.keep_unmapped and rejected:
        with open(args.keep_unmapped, "w") as f:
            f.write("raw_row\treason\n")
            for raw, why in rejected:
                f.write(f"{raw}\t{why}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
