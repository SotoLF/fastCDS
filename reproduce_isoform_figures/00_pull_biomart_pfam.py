#!/usr/bin/env python3
"""Pull Pfam-A domain instances on the human proteome from Ensembl BioMart.

Queries BioMart per chromosome for (ensembl_peptide_id, pfam, pfam_start,
pfam_end) and writes the per-domain table the conservation pipeline consumes as
`--pfam-meta`:

  query_id  protein_id  pfam_id  aa_start  aa_end     (meta TSV)
  protein_id  aa_start  aa_end  query_id             (BED, for fastCDS map)

query_id is a stable `PFAM#######` running id assigned after sorting by
(protein_id, aa_start).

BioMart at www.ensembl.org serves the CURRENT release, so a plain run reflects
whatever release is live now. To pin to a specific release (e.g. the release-115
set used in the paper), point --host at that release's archive BioMart, e.g.
  --host https://<month><year>.archive.ensembl.org/biomart/martservice
(find the archive host for the release on the Ensembl "Archives" page).

Usage:
  00_pull_biomart_pfam.py --out-meta pfam_human_v115_meta.tsv \\
      --out-bed pfam_human_v115.bed
"""
from __future__ import annotations

import argparse
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

DEFAULT_HOST = "https://www.ensembl.org/biomart/martservice"
CHROMS = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]


def query_xml(chrom: str) -> str:
    return ('<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query>'
            '<Query virtualSchemaName="default" formatter="TSV" header="0" '
            'uniqueRows="1" count="" datasetConfigVersion="0.6">'
            '<Dataset name="hsapiens_gene_ensembl" interface="default">'
            f'<Filter name="chromosome_name" value="{chrom}"/>'
            '<Attribute name="ensembl_peptide_id"/>'
            '<Attribute name="pfam"/>'
            '<Attribute name="pfam_start"/>'
            '<Attribute name="pfam_end"/>'
            '</Dataset></Query>')


def fetch(host: str, chrom: str, retries: int = 4) -> str:
    url = host + "?" + urllib.parse.urlencode({"query": query_xml(chrom)})
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=180) as r:
                return r.read().decode()
        except Exception as e:                       # noqa: BLE001 - report + retry
            if attempt == retries - 1:
                print(f"chr{chrom} FAILED: {e}", file=sys.stderr)
                return ""
            time.sleep(5)
    return ""


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-meta", required=True, type=Path,
                    help="query_id/protein_id/pfam_id/aa_start/aa_end TSV")
    ap.add_argument("--out-bed", required=True, type=Path,
                    help="protein_id/aa_start/aa_end/query_id BED (fastCDS map input)")
    ap.add_argument("--host", default=DEFAULT_HOST,
                    help="BioMart martservice URL; swap for a release archive to pin")
    ap.add_argument("--chroms", nargs="+", default=CHROMS,
                    help="chromosomes to query (default: 1-22, X, Y, MT)")
    args = ap.parse_args()

    rows: list[tuple[str, str, int, int]] = []
    seen: set[tuple[str, str, str, str]] = set()
    for c in args.chroms:
        txt = fetch(args.host, c)
        n = 0
        for ln in txt.splitlines():
            f = ln.split("\t")
            if len(f) != 4:
                continue
            pid, pf, s, e = f
            if not pid or not pf.startswith("PF") or not s or not e:
                continue
            key = (pid, pf, s, e)
            if key in seen:
                continue
            seen.add(key)
            rows.append((pid, pf, int(s), int(e)))
            n += 1
        print(f"chr{c}: {n} pfam rows (total {len(rows)})", flush=True)

    if not rows:
        raise SystemExit(
            "00_pull_biomart_pfam.py: BioMart returned no Pfam rows. Check the "
            f"--host URL ({args.host}) is a reachable martservice endpoint.")

    rows.sort(key=lambda r: (r[0], r[2]))
    args.out_meta.parent.mkdir(parents=True, exist_ok=True)
    args.out_bed.parent.mkdir(parents=True, exist_ok=True)
    with args.out_bed.open("w") as bed, args.out_meta.open("w") as meta:
        meta.write("query_id\tprotein_id\tpfam_id\taa_start\taa_end\n")
        for i, (pid, pf, s, e) in enumerate(rows):
            if s < 1 or e < s:
                continue
            qid = f"PFAM{i:07d}"
            bed.write(f"{pid}\t{s}\t{e}\t{qid}\n")
            meta.write(f"{qid}\t{pid}\t{pf}\t{s}\t{e}\n")

    n_prot = len({r[0] for r in rows})
    print(f"\nwrote {len(rows)} Pfam-A instances across {n_prot} proteins "
          f"-> {args.out_meta}, {args.out_bed}")


if __name__ == "__main__":
    main()

