"""Extract Pfam-A domain instances from an EnsDb sqlite to BED for prot2exon.

Used to feed Phase 4 notebook 1 (Pfam proteome atlas). EnsDb.Hsapiens.v86
ships ~150K Pfam-A instances across ~19K proteins — no separate InterPro/UniProt
mapping required.
"""

from __future__ import annotations

import argparse
import sqlite3
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ensdb", required=True, type=Path)
    ap.add_argument("--out-bed", required=True, type=Path)
    ap.add_argument("--out-meta", required=True, type=Path,
                    help="One row per query with category/interpro info")
    args = ap.parse_args()

    con = sqlite3.connect(str(args.ensdb))
    cur = con.cursor()
    rows = cur.execute("""
        SELECT protein_id, protein_domain_id, interpro_accession,
               prot_dom_start, prot_dom_end
        FROM protein_domain
        WHERE protein_domain_source = 'pfam'
        ORDER BY protein_id, prot_dom_start
    """).fetchall()
    con.close()

    bed_lines = []
    meta_lines = ["query_id\tprotein_id\tpfam_id\tinterpro_id\taa_start\taa_end\n"]
    skipped = 0
    for i, (pid, pfam_id, ipr_id, s, e) in enumerate(rows):
        if s is None or e is None or s < 1 or e < s:
            skipped += 1
            continue
        qid = f"PFAM{i:07d}"
        bed_lines.append(f"{pid}\t{s}\t{e}\t{qid}\n")
        meta_lines.append(f"{qid}\t{pid}\t{pfam_id}\t{ipr_id or ''}\t{s}\t{e}\n")

    args.out_bed.write_text("".join(bed_lines))
    args.out_meta.write_text("".join(meta_lines))
    print(f"wrote {len(bed_lines):,} Pfam-A queries (skipped {skipped} with bad coords)")
    print(f"  BED:  {args.out_bed}")
    print(f"  meta: {args.out_meta}")


if __name__ == "__main__":
    main()
