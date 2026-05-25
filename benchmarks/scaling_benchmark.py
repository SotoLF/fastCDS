"""Phase 3 scaling benchmark: prot2exon vs ensembldb at increasing query scales.

Builds query subsets at N = 100, 1k, 10k, 100k (and 1M for prot2exon only),
runs each tool, captures wall time and peak RSS, and writes timings.tsv.

Source queries: existing 5K stratified set under <data>/queries_v86.bed.
Larger sizes (>5K) are produced by sampling with replacement.

ensembldb is much slower than prot2exon; the script caps ensembldb at
--ensembldb-max-n (default 100,000) and runs only --ensembldb-reps reps
for the largest sizes to keep total wall time tractable.
"""

from __future__ import annotations

import argparse
import os
import random
import re
import subprocess
import sys
import time
from pathlib import Path


PEAK_RSS_RE = re.compile(r"^BENCH_PEAK_RSS_MB\s+(\d+)", re.M)


def run_with_time(cmd: list[str], *, env=None) -> dict:
    """Run cmd; return wall_s + peak_rss_mb (via os.wait4 — per-child max RSS)."""
    t0 = time.perf_counter()
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            text=True, env=env)
    pid, status, rusage = os.wait4(proc.pid, 0)
    wall_s = time.perf_counter() - t0
    stdout = proc.stdout.read() if proc.stdout else ""
    stderr = proc.stderr.read() if proc.stderr else ""
    rc = os.waitstatus_to_exitcode(status)
    if rc != 0:
        raise SystemExit(
            f"command failed (rc={rc}): {cmd}\n"
            f"stdout: {stdout}\nstderr: {stderr}"
        )
    # Prefer prot2exon's self-reported peak if present (it's the binary's RSS,
    # not the wrapping shell's). Otherwise use ru_maxrss (KB on Linux).
    m = PEAK_RSS_RE.search(stderr)
    if m:
        return {"wall_s": wall_s, "peak_rss_mb": int(m.group(1))}
    return {"wall_s": wall_s, "peak_rss_mb": rusage.ru_maxrss // 1024}


def sample_bed(src_bed: Path, n: int, dst_bed: Path, seed: int):
    """Sample n queries from src_bed (with replacement if n > source size)."""
    lines = [l for l in src_bed.read_text().splitlines() if l and not l.startswith("#")]
    rng = random.Random(seed)
    if n <= len(lines):
        picked = rng.sample(lines, n)
    else:
        picked = [rng.choice(lines) for _ in range(n)]
    # Rewrite query_id so each row is unique (collisions in p2e summary would lose rows).
    out = []
    for i, line in enumerate(picked):
        parts = line.split("\t")
        # parts = [protein_id, aa_start, aa_end, query_id, ...]
        if len(parts) >= 4:
            parts[3] = f"BENCH{i:07d}"
        out.append("\t".join(parts) + "\n")
    dst_bed.write_text("".join(out))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True, type=Path, help="prot2exon binary")
    ap.add_argument("--p2e-index", required=True, type=Path)
    ap.add_argument("--ensdb", required=True, type=Path, help="EnsDb sqlite")
    ap.add_argument("--rscript", required=True, type=Path)
    ap.add_argument("--r-helper", required=True, type=Path, help="ensembldb_query.R")
    ap.add_argument("--source-bed", required=True, type=Path, help="5K stratified queries")
    ap.add_argument("--work-dir", required=True, type=Path)
    ap.add_argument("--sizes", type=int, nargs="+",
                    default=[100, 1000, 10000, 100000, 1000000])
    ap.add_argument("--ensembldb-max-n", type=int, default=100000,
                    help="Skip ensembldb for sizes above this (it's too slow)")
    ap.add_argument("--p2e-reps", type=int, default=3)
    ap.add_argument("--ensembldb-reps", type=int, default=2)
    ap.add_argument("--ensembldb-reps-large", type=int, default=1,
                    help="Reps for ensembldb at N >= 100,000")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    args.work_dir.mkdir(parents=True, exist_ok=True)
    rows = []

    for n in args.sizes:
        bed = args.work_dir / f"queries_n{n}.bed"
        sample_bed(args.source_bed, n, bed, args.seed + n)
        print(f"[N={n:>7,}] BED ready ({bed.stat().st_size:,} bytes)", file=sys.stderr)

        # --- prot2exon ---
        p2e_out = args.work_dir / f"p2e_n{n}"
        for rep in range(1, args.p2e_reps + 1):
            # Fresh outdir per rep so we time mapping, not "overwrite warnings".
            if p2e_out.exists():
                subprocess.run(["rm", "-rf", str(p2e_out)], check=True)
            r = run_with_time([
                str(args.bin),
                "--index", str(args.p2e_index),
                "--bed", str(bed),
                "--out-dir", str(p2e_out),
                "--output", "coding",
                "--threads", "1",
            ])
            rows.append({
                "tool": "prot2exon", "n": n, "rep": rep, "threads": 1,
                "wall_s": round(r["wall_s"], 3), "peak_rss_mb": r["peak_rss_mb"],
            })
            print(f"  prot2exon rep {rep}: {r['wall_s']:.2f}s, {r['peak_rss_mb']} MB",
                  file=sys.stderr)

        # --- ensembldb (capped at --ensembldb-max-n) ---
        if n > args.ensembldb_max_n:
            print(f"  ensembldb: skipped (N > {args.ensembldb_max_n})", file=sys.stderr)
        else:
            ens_out = args.work_dir / f"ens_n{n}.tsv"
            reps = args.ensembldb_reps_large if n >= 100000 else args.ensembldb_reps
            for rep in range(1, reps + 1):
                r = run_with_time([
                    str(args.rscript), str(args.r_helper),
                    str(args.ensdb), str(bed), str(ens_out),
                ])
                rows.append({
                    "tool": "ensembldb", "n": n, "rep": rep, "threads": 1,
                    "wall_s": round(r["wall_s"], 3), "peak_rss_mb": r["peak_rss_mb"],
                })
                print(f"  ensembldb rep {rep}: {r['wall_s']:.2f}s, {r['peak_rss_mb']} MB",
                      file=sys.stderr)

    # Write the per-rep TSV.
    with open(args.out, "w") as f:
        f.write("tool\tn\trep\tthreads\twall_s\tpeak_rss_mb\n")
        for row in rows:
            f.write(f"{row['tool']}\t{row['n']}\t{row['rep']}\t{row['threads']}\t"
                    f"{row['wall_s']}\t{row['peak_rss_mb']}\n")
    print(f"wrote {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
