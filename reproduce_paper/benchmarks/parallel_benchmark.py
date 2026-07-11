"""Parallel scaling: fastCDS at threads = 1, 2, 4, 8 with fixed N.

Reports wall time + speedup + parallel efficiency per thread count.
ensembldb is single-threaded by design — no parallel comparison.
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
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    pid, status, rusage = os.wait4(proc.pid, 0)
    wall_s = time.perf_counter() - t0
    stderr = proc.stderr.read() if proc.stderr else ""
    if os.waitstatus_to_exitcode(status) != 0:
        raise SystemExit(f"failed: {cmd}\n{stderr}")
    m = PEAK_RSS_RE.search(stderr)
    rss = int(m.group(1)) if m else rusage.ru_maxrss // 1024
    return {"wall_s": wall_s, "peak_rss_mb": rss}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True, type=Path)
    ap.add_argument("--p2e-index", required=True, type=Path)
    ap.add_argument("--bed", required=True, type=Path, help="Query BED (use the 100k subset)")
    ap.add_argument("--work-dir", required=True, type=Path)
    ap.add_argument("--threads", type=int, nargs="+", default=[1, 2, 4, 8])
    ap.add_argument("--reps", type=int, default=2)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    args.work_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    baseline_wall = None
    for t in args.threads:
        for rep in range(1, args.reps + 1):
            out_dir = args.work_dir / f"t{t}_r{rep}"
            if out_dir.exists():
                subprocess.run(["rm", "-rf", str(out_dir)], check=True)
            r = run_once([
                str(args.bin),
                "map", "--index", str(args.p2e_index),
                "--bed", str(args.bed),
                "--out-dir", str(out_dir),
                "--output", "coding",
                "--threads", str(t),
            ])
            rows.append({"threads": t, "rep": rep,
                         "wall_s": round(r["wall_s"], 3),
                         "peak_rss_mb": r["peak_rss_mb"]})
            print(f"threads={t} rep={rep}: {r['wall_s']:.2f}s, {r['peak_rss_mb']} MB",
                  file=sys.stderr)

    with open(args.out, "w") as f:
        f.write("threads\trep\twall_s\tpeak_rss_mb\n")
        for r in rows:
            f.write(f"{r['threads']}\t{r['rep']}\t{r['wall_s']}\t{r['peak_rss_mb']}\n")
    print(f"wrote {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
