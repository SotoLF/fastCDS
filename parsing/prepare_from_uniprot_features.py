#!/usr/bin/env python3
"""Convert UniProt feature annotations into the BED-like format fastCDS eats.

Two input formats are supported, autodetected by extension or via `--format`:

* `.dat` / `.txt` / `.gz` — the classic UniProtKB **flat-file** with `FT`
  lines:

      FT   DOMAIN          93..312
      FT                   /note="Protein kinase"
      FT                   /evidence="ECO:0000255|PROSITE-ProRule:PRU00159"
      FT   ZN_FING         50..73
      FT                   /note="C4-type"

* `.json` — UniProtKB REST output (e.g.
  `https://rest.uniprot.org/uniprotkb/P04637.json`, or a JSONLines dump). We
  read the `features` array. Each feature has `type`, `location`,
  `description` and optionally `featureId`.

The `FT` flat-file format is the one most papers and most pipelines that
predate the REST API still use, so the parser handles it natively.

By default the parser keeps **DOMAIN / REPEAT / ZN_FING / DNA_BIND /
COILED / TRANSMEM / REGION** features (configurable via `--feature-types`).

Mapping UniProt → ENSP is handled by `_mapping.UniProtToEnsp`. UniProt entries
typically carry their own Ensembl cross-references in `DR Ensembl;` lines —
the parser will use those *first* (no external mapping needed) and only fall
back to the Ensembl xref TSV / `--simple-mapping` for entries that don't
carry an Ensembl xref.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(THIS_DIR, "..", "python"))
from fastCDS.prepare._mapping import (
    UniProtToEnsp, looks_like_ensp, looks_like_uniprot,
    open_text, dedup_and_sort_rows, write_bed_row, strip_version,
)
from fastCDS.prepare._uniprot import (
    parse_dat, parse_json, DEFAULT_FEATURE_TYPES,
)


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #

def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description="UniProt feature table (.dat or REST .json) → "
                    "fastCDS BED-like.")
    p.add_argument("--in", dest="in_path", required=True,
                   help="UniProt .dat / .dat.gz / .txt / .json file.")
    p.add_argument("--format", choices=("auto", "dat", "json"), default="auto")
    p.add_argument("--out", default="-")
    p.add_argument("--feature-types",
                   help=f"Comma-separated. Default: {','.join(sorted(DEFAULT_FEATURE_TYPES))}")
    p.add_argument("--mapping",
                   help="Ensembl UniProt xref TSV. Used as a fallback when an "
                        "entry has no DR Ensembl cross-reference.")
    p.add_argument("--simple-mapping",
                   help="Two-column TSV mapping (uniprot\\tensp). Alternative.")
    p.add_argument("--min-length", type=int, default=5)
    p.add_argument("--keep-unmapped",
                   help="Write rejected/unmapped rows here.")
    args = p.parse_args(argv)

    feat_types = (set(s.strip().upper() for s in args.feature_types.split(","))
                  if args.feature_types else DEFAULT_FEATURE_TYPES)

    fmt = args.format
    if fmt == "auto":
        low = args.in_path.lower().rstrip(".gz")
        fmt = "json" if low.endswith(".json") else "dat"
        print(f"[uniprot] auto-detected format = {fmt}", file=sys.stderr)

    if fmt == "dat":
        raw, dr_by_acc, rejected = parse_dat(args.in_path,
                                             feature_types=feat_types)
    else:
        raw, dr_by_acc, rejected = parse_json(args.in_path,
                                              feature_types=feat_types)

    # Optional fallback mapping for entries with no DR Ensembl.
    fallback = None
    if args.mapping:
        fallback = UniProtToEnsp.from_ensembl_xref_tsv(args.mapping)
    elif args.simple_mapping:
        fallback = UniProtToEnsp.from_simple_tsv(args.simple_mapping)

    out_rows: list[tuple[str, int, int, str, str]] = []
    unresolved = 0
    for acc, s, e, did, desc in raw:
        if (e - s + 1) < args.min_length:
            continue
        ensps = dr_by_acc.get(acc, [])
        if not ensps and fallback is not None:
            ensps = fallback.lookup(acc)
        if not ensps:
            unresolved += 1
            rejected.append((f"{acc} {s}-{e} {did}", "no ENSP"))
            continue
        for ensp in ensps:
            out_rows.append((ensp, s, e, did, desc))

    out_rows = dedup_and_sort_rows(out_rows)

    out_fh = sys.stdout if args.out == "-" else open(args.out, "w")
    try:
        out_fh.write("# fastCDS BED-like, generated from UniProt features\n")
        out_fh.write("# columns: ENSP\\taa_start\\taa_end\\tdomain_id\\tdescription\n")
        for r in out_rows:
            write_bed_row(out_fh, ensp=r[0], aa_start=r[1], aa_end=r[2],
                          domain_id=r[3], source=r[4])
    finally:
        if out_fh is not sys.stdout:
            out_fh.close()
    print(f"[uniprot] wrote {len(out_rows)} rows ({unresolved} unresolved, "
          f"{len(rejected)} total rejected) to {args.out}", file=sys.stderr)
    if args.keep_unmapped and rejected:
        with open(args.keep_unmapped, "w") as f:
            f.write("raw\treason\n")
            for raw, why in rejected:
                f.write(f"{raw}\t{why}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
