"""Build the shared query ladder every tool in the matched benchmark runs on.

Pulls every (protein_id, transcript_id) pair with a CDS out of the GTF, sorts
them for determinism, and cycles that list to fill each N. Cycling matters above
~94K: human v86 has fewer coding proteins than the top of the ladder, and the
alternative (sampling with replacement) would give each tool a different query
mix. Each ladder step is a prefix of the next, so N = 100 is the first 100 rows
of N = 1,000, and no tool is handed an easier subset than another.

Every query asks for the same aa 1..50 window, so the benchmark measures the
mapping itself rather than how much sequence each tool happened to receive.

Emits, per N:
    q_<N>.bed      protein_id, aa_start, aa_end, domain_id   (fastCDS)
    ids_<N>.txt    ENSP per line                             (ensembldb, geneplot, REST)
    enst_<N>.txt   ENST per line                             (GenomicFeatures route)

Usage:
    python prepare_matched_inputs.py --gtf Homo_sapiens.GRCh38.86.chr.gtf \
        --out-dir "$FASTCDS_DATA/matched"
"""

from __future__ import annotations

import argparse
import gzip
import re
from itertools import islice, cycle
from pathlib import Path

LADDER = [100, 1000, 10000, 100000, 1000000]

_PID = re.compile(r'protein_id "([^"]+)"')
_TID = re.compile(r'transcript_id "([^"]+)"')


def read_pairs(gtf: Path) -> list[tuple[str, str]]:
    """Unique (ENSP, ENST) pairs from the GTF's CDS rows, sorted by ENSP.

    Ensembl puts protein_id only on CDS rows, so those are the rows to read.
    """
    opener = gzip.open if gtf.suffix == ".gz" else open
    pairs: dict[str, str] = {}
    with opener(gtf, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.split("\t", 9)
            if len(cols) < 9 or cols[2] != "CDS":
                continue
            pid = _PID.search(cols[8])
            tid = _TID.search(cols[8])
            if pid and tid:
                pairs.setdefault(pid.group(1), tid.group(1))
    return sorted(pairs.items())


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gtf", required=True, type=Path, help="Ensembl GTF (.gtf or .gtf.gz)")
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--sizes", nargs="+", type=int, default=LADDER)
    ap.add_argument("--aa-start", type=int, default=1)
    ap.add_argument("--aa-end", type=int, default=50)
    args = ap.parse_args()

    if not args.gtf.exists():
        raise SystemExit(f"prepare_matched_inputs.py: --gtf not found: {args.gtf}")

    pairs = read_pairs(args.gtf)
    if not pairs:
        raise SystemExit(
            f"prepare_matched_inputs.py: no CDS rows with protein_id in {args.gtf}. "
            f"Expected an Ensembl GTF (GENCODE puts protein_id on transcript rows).")
    print(f"{len(pairs):,} unique coding proteins in {args.gtf.name}")

    out = args.out_dir.expanduser()
    out.mkdir(parents=True, exist_ok=True)
    width = len(str(max(args.sizes) - 1))

    for n in sorted(args.sizes):
        rows = list(islice(cycle(pairs), n))
        with open(out / f"q_{n}.bed", "w") as bed, \
             open(out / f"ids_{n}.txt", "w") as ids, \
             open(out / f"enst_{n}.txt", "w") as enst:
            for i, (pid, tid) in enumerate(rows):
                bed.write(f"{pid}\t{args.aa_start}\t{args.aa_end}\tQ{i:0{width}d}\n")
                ids.write(pid + "\n")
                enst.write(tid + "\n")
        print(f"  N={n:<9,} -> q_{n}.bed, ids_{n}.txt, enst_{n}.txt")

    print("wrote to", out)


if __name__ == "__main__":
    main()
