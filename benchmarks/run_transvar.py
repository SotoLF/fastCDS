"""Query TransVar's `panno` for a batch of (protein_id, aa_start, aa_end, query_id)
rows and emit a TSV in the same shape as ensembldb_query.R so the comparison
logic in validate_vs_ensembldb.py can reuse the classifier.

TransVar accepts HGVS protein notation: `<ENSP>:p.<AA><start>_<AA><end>`. We
look up the actual AA letters at the two endpoints from EnsDb's `protein` table
(no need for full sequences elsewhere).

Usage:
    python run_transvar.py <queries.bed> <out.tsv> <ensdb.sqlite> \
        [--transvar /path/to/transvar] [--db ensembl] [--batch 100]
"""

from __future__ import annotations

import argparse
import shutil
import sqlite3
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path


def fetch_protein_seqs(ensdb: Path, ids: list[str], by: str = "protein_id") -> dict[str, str]:
    """Look up protein_sequence keyed by `by` (either "protein_id" or "tx_id")."""
    assert by in ("protein_id", "tx_id")
    con = sqlite3.connect(str(ensdb))
    cur = con.cursor()
    seqs: dict[str, str] = {}
    chunk = 500
    for i in range(0, len(ids), chunk):
        slug = ids[i:i + chunk]
        placeholders = ",".join("?" * len(slug))
        rows = cur.execute(
            f"SELECT {by}, protein_sequence FROM protein WHERE {by} IN ({placeholders})",
            slug,
        ).fetchall()
        for k, seq in rows:
            seqs[k] = seq
    con.close()
    return seqs


def aa_at(seq: str, pos: int) -> str | None:
    """1-based residue lookup, returning None if out of range or unavailable."""
    if not seq or pos < 1 or pos > len(seq):
        return None
    return seq[pos - 1]


def parse_transvar_coordinates(s: str) -> list[tuple[str, int, int]]:
    """Pull (chrom, start, end) intervals from TransVar's coordinates field.

    Format examples we see:
        chr3:g.179218303G>A/c.1633G>A/p.E545K
        chr17:g.7676086_7676272/c.1_187/p.M1_S62
        chr17:g.(7676086_7676215)_(7676219_7676272)/...
    """
    # We only need the genomic span. Take the part before the first '/'.
    g = s.split("/", 1)[0]
    if not g.startswith("chr") and ":" not in g:
        return []
    chrom, _, rest = g.partition(":")
    # Strip leading 'g.'
    if rest.startswith("g."):
        rest = rest[2:]
    # Drop any base-change suffix like '534G>A' -> keep coordinate range only.
    # Forms: 'NNN_MMM' (range) or 'NNN' (point) or '(a_b)_(c_d)' (ambiguous endpoints).
    out: list[tuple[str, int, int]] = []
    # Handle ambiguous parentheses by stripping them and taking the min start, max end.
    cleaned = rest.replace("(", "").replace(")", "")
    if "_" in cleaned:
        # 'a_b' or 'a_b_c_d' etc.
        parts = cleaned.split("_")
        try:
            nums = [int(p.split(">")[0].rstrip("ATCGatcg")) for p in parts if p and p[0].isdigit()]
        except ValueError:
            nums = []
        if len(nums) >= 2:
            out.append((chrom, min(nums), max(nums)))
    else:
        try:
            pos = int(cleaned.split(">")[0].rstrip("ATCGatcg"))
            out.append((chrom, pos, pos))
        except ValueError:
            pass
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("bed", type=Path)
    ap.add_argument("out", type=Path)
    ap.add_argument("--ensdb", type=Path, required=True,
                    help="EnsDb sqlite to look up AA letters from (must have protein table)")
    ap.add_argument("--queries-meta", type=Path, default=None,
                    help="Optional meta TSV mapping query_id -> transcript_id. "
                         "Required because TransVar's panno keys on ENST, not ENSP.")
    ap.add_argument("--transvar", default=shutil.which("transvar") or "transvar")
    ap.add_argument("--db", default="ensembl", choices=["ensembl", "refseq", "ucsc"])
    ap.add_argument("--refversion", default="hg38")
    ap.add_argument("--limit", type=int, default=10000)
    args = ap.parse_args()

    # Optional ENSP->ENST mapping from queries metadata.
    qid_to_enst: dict[str, str] = {}
    if args.queries_meta:
        import csv as _csv
        with open(args.queries_meta) as f:
            for row in _csv.DictReader(f, delimiter="\t"):
                qid_to_enst[row["query_id"]] = row["transcript_id"]

    # Load queries (cap at --limit).
    rows = []
    with open(args.bed) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            rows.append((parts[0], int(parts[1]), int(parts[2]), parts[3]))
            if len(rows) >= args.limit:
                break
    print(f"{len(rows):,} queries loaded", file=sys.stderr)

    # Look up AA letters via the tx_id (ENST) since that's what TransVar uses.
    if qid_to_enst:
        unique_ensts = sorted({qid_to_enst.get(r[3], r[0]) for r in rows})
        key = "tx_id"
        print(f"  using ENST from queries_meta ({len(unique_ensts):,} unique)", file=sys.stderr)
    else:
        unique_ensts = sorted({r[0] for r in rows})
        key = "protein_id"
    print(f"  fetching {len(unique_ensts):,} protein sequences from EnsDb ...", file=sys.stderr)
    seqs = fetch_protein_seqs(args.ensdb, unique_ensts, by=key)
    print(f"  got {len(seqs):,} sequences ({len(unique_ensts) - len(seqs):,} missing)", file=sys.stderr)

    # Build HGVS strings, recording skipped queries with reason.
    hgvs_pairs = []   # (query_id, hgvs_str)
    skipped = {}      # query_id -> reason
    for pid, s, e, qid in rows:
        # Use ENST if we have it; TransVar requires ENST for ensembl mode.
        ident = qid_to_enst.get(qid, pid)
        seq = seqs.get(ident)
        if seq is None:
            skipped[qid] = "protein_not_in_ensdb"
            continue
        aa_s = aa_at(seq, s)
        aa_e = aa_at(seq, e)
        if aa_s is None or aa_e is None:
            skipped[qid] = "aa_pos_out_of_range"
            continue
        if s == e:
            hgvs = f"{ident}:p.{aa_s}{s}"
        else:
            hgvs = f"{ident}:p.{aa_s}{s}_{aa_e}{e}"
        hgvs_pairs.append((qid, hgvs))
    print(f"  built {len(hgvs_pairs):,} HGVS strings ({len(skipped):,} skipped)", file=sys.stderr)

    # Run TransVar via stdin (one HGVS per line). Use --noheader for clean parsing.
    cmd = [args.transvar, "panno", "-l", "/dev/stdin",
           f"--{args.db}", "--refversion", args.refversion, "--noheader"]
    print(f"running: {' '.join(cmd)}", file=sys.stderr)
    stdin_str = "\n".join(h for _, h in hgvs_pairs) + "\n"
    t0 = time.perf_counter()
    proc = subprocess.run(cmd, input=stdin_str, capture_output=True, text=True)
    wall_s = time.perf_counter() - t0
    print(f"transvar finished in {wall_s:.1f}s (rc={proc.returncode})", file=sys.stderr)
    if proc.returncode != 0:
        sys.stderr.write(proc.stderr[:1000] + "\n")
        # Continue — TransVar often returns non-zero on partial failures but still writes results.

    # Map HGVS back to query_id.
    hgvs_to_qid = {h: q for q, h in hgvs_pairs}

    # Parse output. With --noheader, columns are tab-separated: input, transcript, gene, strand, coordinates, ...
    with open(args.out, "w") as f:
        f.write("query_id\tchrom\tstart\tend\tstrand\tstatus\n")
        for qid, reason in skipped.items():
            f.write(f"{qid}\tNA\tNA\tNA\tNA\t{reason}\n")
        # Group results by input HGVS string. TransVar can emit multiple lines per input
        # (one per matching transcript). We keep all and pick the first non-error one.
        by_input: dict = defaultdict(list)
        for line in proc.stdout.splitlines():
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 5:
                continue
            hgvs_in = parts[0]
            coords = parts[4]
            by_input[hgvs_in].append((parts, coords))

        for qid, hgvs in hgvs_pairs:
            results = by_input.get(hgvs, [])
            if not results:
                f.write(f"{qid}\tNA\tNA\tNA\tNA\tno_result\n")
                continue
            # Prefer rows where coordinates parse non-trivially.
            picked = None
            for parts, coords in results:
                intervals = parse_transvar_coordinates(coords)
                if intervals:
                    picked = (parts, intervals)
                    break
            if picked is None:
                f.write(f"{qid}\tNA\tNA\tNA\tNA\terror\n")
                continue
            parts, intervals = picked
            strand = parts[3] if len(parts) > 3 else "."
            for c, gs, ge in intervals:
                f.write(f"{qid}\t{c}\t{gs}\t{ge}\t{strand}\tok\n")

    print(f"wrote {args.out}", file=sys.stderr)
    print(f"WALL_S {wall_s:.3f}", file=sys.stderr)


if __name__ == "__main__":
    main()
