"""Combined threads × batch-size grid for `prot2exon map`.

Supersedes the separate parallel (threads) and batch-size sweeps: it varies both
axes together so you can read wall-time *and* peak RAM off one grid. `--threads`
trades wall time; `--batch-size` caps peak RSS by streaming results to disk in
chunks (0 / omitted = one-shot, hold everything in memory).

    python benchmarks/threads_batch_grid.py \
        --bin build/prot2exon --index human_v86.idx \
        --bed bench/queries_n200000.bed --work-dir bench/grid \
        --threads 1 2 4 8 --batch-sizes 0 10000 50000 \
        --out bench/threads_batch_grid.tsv
"""
from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
import time
from pathlib import Path

PEAK_RSS_RE = re.compile(r"^BENCH_PEAK_RSS_MB\s+(\d+)", re.M)


def run_once(cmd: list[str]) -> dict:
    t0 = time.perf_counter()
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            text=True)
    _pid, status, rusage = os.wait4(proc.pid, 0)
    wall_s = time.perf_counter() - t0
    stderr = proc.stderr.read() if proc.stderr else ""
    if os.waitstatus_to_exitcode(status) != 0:
        raise SystemExit(f"failed: {' '.join(cmd)}\n{stderr}")
    m = PEAK_RSS_RE.search(stderr)
    rss = int(m.group(1)) if m else rusage.ru_maxrss // 1024
    return {"wall_s": wall_s, "peak_rss_mb": rss}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True, type=Path)
    ap.add_argument("--index", required=True, type=Path)
    ap.add_argument("--bed", required=True, type=Path)
    ap.add_argument("--work-dir", required=True, type=Path)
    ap.add_argument("--threads", type=int, nargs="+", default=[1, 2, 4, 8])
    ap.add_argument("--batch-sizes", type=int, nargs="+", default=[0, 10000, 50000],
                    help="0 = one-shot (no --batch-size)")
    ap.add_argument("--output", default="coding")
    ap.add_argument("--reps", type=int, default=1)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    n_queries = sum(1 for _ in open(args.bed))
    args.work_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for t in args.threads:
        for b in args.batch_sizes:
            best = None
            for rep in range(1, args.reps + 1):
                out_dir = args.work_dir / f"t{t}_b{b}_r{rep}"
                if out_dir.exists():
                    subprocess.run(["rm", "-rf", str(out_dir)], check=True)
                cmd = [str(args.bin), "map", "--index", str(args.index),
                       "--bed", str(args.bed), "--out-dir", str(out_dir),
                       "--output", args.output, "--threads", str(t)]
                if b > 0:
                    cmd += ["--batch-size", str(b)]
                r = run_once(cmd)
                subprocess.run(["rm", "-rf", str(out_dir)], check=True)
                # keep the fastest rep (least noisy)
                if best is None or r["wall_s"] < best["wall_s"]:
                    best = r
            label = "one-shot" if b == 0 else f"{b:,}"
            qps = n_queries / best["wall_s"]
            rows.append({"threads": t, "batch_size": b, "batch_label": label,
                         **best, "qps": qps})
            sys.stderr.write(
                f"threads={t} batch={label:>8}  "
                f"wall={best['wall_s']:6.1f}s  "
                f"rss={best['peak_rss_mb']:6d}MB  "
                f"{qps:8.0f} q/s\n")

    # write TSV
    cols = ["threads", "batch_size", "batch_label", "wall_s", "peak_rss_mb", "qps"]
    with open(args.out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(f"{r[c]:.3f}" if c in ("wall_s", "qps")
                               else str(r[c]) for c in cols) + "\n")
    sys.stderr.write(f"\nwrote {args.out}  (N = {n_queries:,} queries, "
                     f"output={args.output})\n")


if __name__ == "__main__":
    main()
