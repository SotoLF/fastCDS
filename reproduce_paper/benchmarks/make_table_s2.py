"""Build Supplementary Table S2 — mapping speed and peak memory.

Single-core performance for each tool: throughput (queries per second) and peak
resident set size, each measured at the tool's benchmark query count. fastCDS,
GenomicFeatures, and ensembldb are benchmarked at N=10,000; VisProDom, geneplot,
and the Ensembl REST API at N=1,000 (VisProDom/geneplot are slower and REST is
rate-limited). Wall time is measured from process launch until all results are
written, including index/database/annotation loading. REST's peak RSS reflects
only the local client (mapping runs on the Ensembl server).

Reads the scaling TSV written by scaling_benchmark.py (columns: tool, n,
wall_s, peak_rss_mb) and emits the table sorted by throughput.

(The per-category coordinate-agreement table is Supplementary Table S1; build it
with make_table_s1.py.)

Usage:
    python make_table_s2.py --scaling-tsv scaling.tsv --out table_S2_speed_memory.tsv
"""

from __future__ import annotations

import argparse
import statistics
from pathlib import Path

# Query count at which each tool is reported (matches the manuscript caption).
BENCH_N = {
    "fastCDS": 10000, "GenomicFeatures": 10000, "ensembldb": 10000,
    "VisProDom": 1000, "geneplot": 1000, "Ensembl_REST": 1000,
}
DISPLAY = {"Ensembl_REST": "Ensembl REST", "GenomicFeatures": "GenomicFeatures (GRanges)"}


def read_tsv(path: Path) -> list[dict]:
    rows = []
    with open(path) as f:
        header = None
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            header = line.rstrip("\n").split("\t")
            break
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            r = {}
            for k, v in zip(header, parts):
                try:
                    r[k] = float(v) if "." in v else int(v)
                except ValueError:
                    r[k] = v
            rows.append(r)
    return rows


def median_by(rows, *keys, value):
    buckets: dict = {}
    for r in rows:
        buckets.setdefault(tuple(r[k] for k in keys), []).append(r[value])
    return {k: statistics.median(v) for k, v in buckets.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scaling-tsv", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    rows = read_tsv(args.scaling_tsv)
    med_wall = median_by(rows, "tool", "n", value="wall_s")
    med_rss = median_by(rows, "tool", "n", value="peak_rss_mb")

    out = [("tool", "queries_per_s", "peak_rss_mb", "n")]
    recs = []
    for tool, n in BENCH_N.items():
        wall = med_wall.get((tool, n))
        if wall is None:
            continue
        qps = n / wall
        rss = med_rss.get((tool, n), 0)
        recs.append((DISPLAY.get(tool, tool), qps, int(round(rss)), n))
    recs.sort(key=lambda r: r[1], reverse=True)
    for name, qps, rss, n in recs:
        out.append((name, f"{qps:.1f}" if qps < 100 else f"{qps:.0f}", str(rss), str(n)))

    with open(args.out, "w") as f:
        for row in out:
            f.write("\t".join(map(str, row)) + "\n")
    print(f"wrote {args.out}\n")
    widths = [max(len(str(r[i])) for r in out) for i in range(len(out[0]))]
    for r in out:
        print("  ".join(str(c).ljust(w) for c, w in zip(r, widths)))


if __name__ == "__main__":
    main()
