#!/usr/bin/env python3
"""Step 1 - list the catalytic Pfam families.

A Pfam family counts as catalytic if InterPro pfam2go maps it to GO:0003824
(catalytic activity) or to any descendant of that term. The descendant set is
read from EBI QuickGO; the Pfam-to-GO map is the InterPro pfam2go file.

Output
  pfam_catalytic.tsv   one row per catalytic Pfam family (column: pfam_id)

Input
  pfam2go   InterPro Pfam-to-GO map
            (current.geneontology.org/ontology/external2go/pfam2go)

Run
  python 01_pfam_catalytic.py --pfam2go pfam2go --out pfam_catalytic.tsv
"""
import argparse
import json
import re
import urllib.request
from pathlib import Path

QUICKGO = ("https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/"
           "GO:0003824/descendants?relations=is_a")


def catalytic_go_terms():
    """GO:0003824 (catalytic activity) and all of its is_a descendants."""
    req = urllib.request.Request(QUICKGO, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req) as r:
        data = json.load(r)
    terms = {"GO:0003824"}
    for res in data.get("results", []):
        terms.update(res.get("descendants") or [])
    return terms


def pfam_to_go(path):
    """InterPro pfam2go: each line maps one Pfam family to one GO id."""
    m = {}
    for line in open(path):
        if line.startswith("!"):
            continue
        pf = re.match(r"Pfam:(PF\d+)", line)
        go = re.search(r"(GO:\d+)\s*$", line.rstrip())
        if pf and go:
            m.setdefault(pf.group(1), set()).add(go.group(1))
    return m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pfam2go", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    catalytic = catalytic_go_terms()
    pfam2go = pfam_to_go(args.pfam2go)
    families = sorted(pf for pf, gos in pfam2go.items() if gos & catalytic)

    with open(args.out, "w") as f:
        f.write("pfam_id\n")
        f.write("\n".join(families) + "\n")

    print(f"catalytic GO terms: {len(catalytic):,}")
    print(f"catalytic Pfam families: {len(families):,} -> {args.out}")


if __name__ == "__main__":
    main()
