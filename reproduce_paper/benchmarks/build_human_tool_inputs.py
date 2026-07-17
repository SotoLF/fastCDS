#!/usr/bin/env python3
"""Build geneplot's inputs for the human Ensembl-86 set.

geneplot ships only fruit-fly example data, but it is general (any GFF plus a
domain file), so it can be run on the same human workload as fastCDS,
ensembldb and GenomicFeatures. It wants an InterProScan `.ipr` table and, for
the benchmark harness, an ENSP -> ENST map. Both are derived here from inputs
we already use elsewhere, so no new data source is introduced.

Outputs (in --out-dir):
  h86.ipr        InterProScan-format domain table, one row per Pfam hit
  ensp_enst.tsv  ENSP <tab> ENST, one row per coding protein

The `.ipr` layout is InterProScan's 13-column TSV. geneplot reads the protein,
the signature accession, the aa start/end and the InterPro accession; the
remaining columns are filled with the literals InterProScan emits for a Pfam
hit carrying no description or GO terms.

geneplot also needs the Ensembl-86 GFF3 itself, which it turns into a gffutils
SQLite database on first use. That file is downloaded from Ensembl release-86
as-is and needs no preparation here.

Provenance note: run against `Homo_sapiens.GRCh38.86.chr.gtf` this reproduces
h86.ipr byte for byte, and 94,347 of the 94,384 ENSP -> ENST pairs used in the
published run. The 37 it does not emit sit on non-chromosomal scaffolds, which
the `.chr` GTF excludes by definition; pass a GTF that includes scaffolds to
recover them. They are 0.04% of the set and carry no timing weight.

Usage:
    python build_human_tool_inputs.py \\
        --gtf   Homo_sapiens.GRCh38.86.chr.gtf \\
        --pfam  pfam_human_v86_meta.tsv \\
        --out-dir "$FASTCDS_DATA/human_tool_bench"

Then benchmark with ../matched/run_geneplot.py (see matched/README.md).
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
from pathlib import Path

_PID = re.compile(r'protein_id "([^"]+)"')
_TID = re.compile(r'transcript_id "([^"]+)"')


def write_ensp_enst(gtf: Path, out: Path) -> int:
    """ENSP -> ENST from the GTF's CDS rows (Ensembl puts protein_id there)."""
    opener = gzip.open if gtf.suffix == ".gz" else open
    pairs: dict[str, str] = {}
    with opener(gtf, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.split("\t", 9)
            if len(cols) < 9 or cols[2] != "CDS":
                continue
            pid, tid = _PID.search(cols[8]), _TID.search(cols[8])
            if pid and tid:
                pairs.setdefault(pid.group(1), tid.group(1))
    if not pairs:
        raise SystemExit(
            f"build_human_tool_inputs.py: no CDS rows with protein_id in {gtf}. "
            f"Expected an Ensembl GTF (GENCODE puts protein_id on transcript rows).")
    with open(out, "w") as f:
        for pid, tid in sorted(pairs.items()):
            f.write(f"{pid}\t{tid}\n")
    return len(pairs)


def write_ipr(pfam: Path, out: Path) -> int:
    """Pfam-on-v86 table -> InterProScan 13-column .ipr."""
    need = {"protein_id", "pfam_id", "interpro_id", "aa_start", "aa_end"}
    rows = 0
    with open(pfam, newline="") as fh, open(out, "w") as f:
        rd = csv.DictReader(fh, delimiter="\t")
        missing = need - set(rd.fieldnames or [])
        if missing:
            raise SystemExit(
                f"build_human_tool_inputs.py: --pfam {pfam} is missing "
                f"column(s): {', '.join(sorted(missing))}. "
                f"Expected a TSV with: {', '.join(sorted(need))}.")
        for r in rd:
            f.write("\t".join([
                r["protein_id"],
                "-",                                  # md5 (unused by geneplot)
                "0",                                  # sequence length (unused)
                "Pfam",                               # analysis
                r["pfam_id"],                         # signature accession
                f"{r['pfam_id']} domain",             # signature description
                r["aa_start"], r["aa_end"],
                "-",                                  # e-value (unused)
                "T",                                  # match status
                "-",                                  # date (unused)
                r["interpro_id"],
                "-",                                  # InterPro description (unused)
            ]) + "\n")
            rows += 1
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gtf", required=True, type=Path,
                    help="Ensembl release-86 GTF (.gtf or .gtf.gz)")
    ap.add_argument("--pfam", required=True, type=Path,
                    help="pfam_human_v86_meta.tsv (query_id, protein_id, pfam_id, "
                         "interpro_id, aa_start, aa_end)")
    ap.add_argument("--out-dir", required=True, type=Path)
    args = ap.parse_args()

    for p, flag in [(args.gtf, "--gtf"), (args.pfam, "--pfam")]:
        if not p.exists():
            raise SystemExit(f"build_human_tool_inputs.py: {flag} not found: {p}")

    out = args.out_dir.expanduser()
    out.mkdir(parents=True, exist_ok=True)

    n = write_ensp_enst(args.gtf, out / "ensp_enst.tsv")
    print(f"wrote {out / 'ensp_enst.tsv'}  ({n:,} proteins)")
    n = write_ipr(args.pfam, out / "h86.ipr")
    print(f"wrote {out / 'h86.ipr'}  ({n:,} domain rows)")


if __name__ == "__main__":
    main()
