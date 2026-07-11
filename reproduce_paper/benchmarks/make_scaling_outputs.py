"""Generate paper benchmark outputs from the raw timing TSVs.

Inputs:
  --scaling-tsv     output of scaling_benchmark.py (per-rep timings)
  --parallel-tsv    output of parallel_benchmark.py (per-rep timings at varying threads)
  --p2e-index-size  bytes — `ls -l human_v86.idx` (or any fastCDS index)
  --ensdb-size      bytes — EnsDb sqlite file size
  --agreement       OVERALL exact-match percentage from matched-annotation validation (e.g. 100.00)

Outputs in --out-dir:
  paper speed/memory summary (now Supplementary Table S2)
  scaling.png        log-log runtime + parallel efficiency, two panels
  parallel.tsv       per-thread median wall + speedup + efficiency
"""

from __future__ import annotations

import argparse
import statistics
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_tsv(path: Path) -> list[dict]:
    rows = []
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        for line in f:
            parts = line.rstrip("\n").split("\t")
            row = {}
            for k, v in zip(header, parts):
                try:
                    row[k] = float(v) if "." in v else int(v)
                except ValueError:
                    row[k] = v
            rows.append(row)
    return rows


def median_by(rows, *keys, value="wall_s"):
    """Group rows by `keys`, return {key_tuple: median_value, ...}."""
    buckets: dict = {}
    for r in rows:
        k = tuple(r[x] for x in keys)
        buckets.setdefault(k, []).append(r[value])
    return {k: statistics.median(v) for k, v in buckets.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scaling-tsv", required=True, type=Path)
    ap.add_argument("--parallel-tsv", required=True, type=Path)
    ap.add_argument("--p2e-index-size", type=int, required=True, help="bytes")
    ap.add_argument("--ensdb-size", type=int, required=True, help="bytes")
    ap.add_argument("--agreement", type=float, default=100.0,
                    help="OVERALL exact-match percentage from matched-annotation validation")
    ap.add_argument("--out-dir", required=True, type=Path)
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    scaling = read_tsv(args.scaling_tsv)
    parallel = read_tsv(args.parallel_tsv)

    # Median wall and RSS per (tool, n).
    med_wall = median_by(scaling, "tool", "n", value="wall_s")
    med_rss = median_by(scaling, "tool", "n", value="peak_rss_mb")

    # ----- speed/memory summary (Supplementary Table S2) -----
    # Pick reference N = 10,000 for "Runtime 10K" since both tools have it.
    REF_N = 10000

    def fmt_runtime(tool, n):
        return f"{med_wall.get((tool, n), float('nan')):.1f} s"

    def fmt_rss(tool, n):
        v = med_rss.get((tool, n))
        return f"{v} MB" if v else "—"

    def fmt_throughput(tool, n):
        w = med_wall.get((tool, n))
        return f"{int(n / w):,} q/s" if w else "—"

    # Largest N each tool covered.
    p2e_max_n = max((n for (t, n) in med_wall if t == "fastCDS"), default=0)
    ens_max_n = max((n for (t, n) in med_wall if t == "ensembldb"), default=0)

    table_path = args.out_dir / "table1.tsv"
    with open(table_path, "w") as f:
        f.write("metric\tfastCDS\tensembldb\n")
        f.write(f"Exact agreement vs ensembldb\t100.00% (ref)\t{args.agreement:.2f}%\n")
        f.write(f"Runtime N={REF_N:,} (1 thread, median)\t{fmt_runtime('fastCDS', REF_N)}\t{fmt_runtime('ensembldb', REF_N)}\n")
        f.write(f"Peak RSS N={REF_N:,}\t{fmt_rss('fastCDS', REF_N)}\t{fmt_rss('ensembldb', REF_N)}\n")
        f.write(f"Throughput N={REF_N:,}\t{fmt_throughput('fastCDS', REF_N)}\t{fmt_throughput('ensembldb', REF_N)}\n")
        f.write(f"Largest N benchmarked\t{p2e_max_n:,}\t{ens_max_n:,}\n")
        f.write(f"Index size on disk\t{args.p2e_index_size/(1024*1024):.1f} MB\t{args.ensdb_size/(1024*1024):.1f} MB\n")
        f.write(f"Parallelism (OpenMP)\tYes\tNo\n")
        f.write(f"Plot-ready output\tYes\tNo\n")
        f.write(f"Multi-species support\tYes (any GTF)\tYes (any Ensembl release)\n")
    print(f"wrote {table_path}")
    print(open(table_path).read())

    # ----- Parallel scaling -----
    # median_by returns {tuple_key: value}; for single-key grouping the keys are 1-tuples.
    par_med_raw = median_by(parallel, "threads", value="wall_s")
    par_med = {k[0]: v for k, v in par_med_raw.items()}
    threads_sorted = sorted(par_med)
    baseline = par_med[threads_sorted[0]]  # T=1 wall
    par_path = args.out_dir / "parallel.tsv"
    with open(par_path, "w") as f:
        f.write("threads\twall_s\tspeedup\tefficiency\n")
        for t in threads_sorted:
            sp = baseline / par_med[t] if par_med[t] else 0
            eff = sp / t
            f.write(f"{t}\t{par_med[t]:.3f}\t{sp:.2f}\t{eff:.2f}\n")
    print(f"wrote {par_path}")
    print(open(par_path).read())

    # ----- Scaling figure -----
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))

    # Panel A: runtime vs N (log-log), one line per tool.
    for tool, marker in [("fastCDS", "o"), ("ensembldb", "s")]:
        xs = sorted(n for (t, n) in med_wall if t == tool)
        ys = [med_wall[(tool, n)] for n in xs]
        if xs:
            ax1.loglog(xs, ys, marker=marker, label=tool, linewidth=2, markersize=8)
    ax1.set_xlabel("Queries (N)")
    ax1.set_ylabel("Wall time (s)")
    ax1.set_title("Single-thread scaling (log-log)")
    ax1.grid(True, which="both", alpha=0.3)
    ax1.legend()

    # Panel B: parallel efficiency vs thread count.
    speedups = [baseline / par_med[t] for t in threads_sorted]
    efficiencies = [s / t for s, t in zip(speedups, threads_sorted)]
    ax2.plot(threads_sorted, speedups, "o-", label="speedup", linewidth=2, markersize=8)
    ax2.plot(threads_sorted, threads_sorted, "--", color="gray", alpha=0.5, label="ideal")
    ax2.set_xlabel("Threads")
    ax2.set_ylabel("Speedup vs single-thread")
    ax2_b = ax2.twinx()
    ax2_b.plot(threads_sorted, efficiencies, "v--", color="C2", label="efficiency", markersize=6)
    ax2_b.set_ylabel("Parallel efficiency (speedup / threads)")
    ax2_b.set_ylim(0, 1.1)
    ax2.set_title(f"fastCDS parallel scaling (N=100k)")
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc="upper left")
    ax2_b.legend(loc="lower right")

    fig.tight_layout()
    png_path = args.out_dir / "scaling.png"
    fig.savefig(png_path, dpi=150)
    print(f"wrote {png_path}")


if __name__ == "__main__":
    main()
