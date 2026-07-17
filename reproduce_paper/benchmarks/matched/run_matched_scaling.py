"""Matched scaling benchmark: every tool, same machine, same queries, one thread.

Runs each tool over the query ladder (N = 100 ... 1,000,000) under a wait4()
harness that records wall time and peak RSS, appends every measurement to
results.tsv, then assembles the tidy table the scaling_and_ram notebook reads.

Each tool is pushed to the largest N it can finish in a practical time/memory
budget and stopped there, so the per-tool ladders differ (see --sizes).

Usage:
    # full sweep, then assemble
    python run_matched_scaling.py --work-dir "$FASTCDS_DATA/matched" \
        --bin build/fastCDS --index "$FASTCDS_DATA/human_v86.idx" \
        --ensdb "$EnsDb_v86_path"

    # one tool only (re-measure), then re-assemble
    python run_matched_scaling.py --work-dir ... --tools ensembldb --sizes 1000
    python run_matched_scaling.py --work-dir ... --assemble-only

Inputs expected in --work-dir (built by prepare_matched_inputs.py):
    q_<N>.bed        fastCDS query BEDs        (protein_id, aa_start, aa_end, domain_id)
    ids_<N>.txt      ENSP ids, one per line    (ensembldb / geneplot / REST)
    enst_<N>.txt     ENST ids, one per line    (GenomicFeatures GRanges route)
    cds_by_tx.rds    CDS-by-transcript GRangesList (GenomicFeatures route)
"""

from __future__ import annotations

import argparse
import os
import resource
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent

# The N ladder each tool can actually finish. fastCDS reaches 1e6; the
# Bioconductor routes and geneplot stall at 1e4; REST is network-bound at 1e3.
LADDER = {
    "fastCDS":         [100, 1000, 10000, 100000, 1000000],
    "ensembldb":       [100, 1000, 10000],
    "genomicfeatures": [100, 1000, 10000],
    "geneplot":        [100, 1000, 10000],
    "rest":            [100, 1000],
}

# Tools whose curve is reported in the paper. Anything else measured into
# results.tsv is ignored at assemble time.
PAPER_TOOLS = list(LADDER)


def measure(tool: str, n: int, cmd: list[str], work: Path) -> dict:
    """Run cmd to completion; return wall seconds and peak RSS in MB.

    Uses fork + wait4 so the child's own peak RSS is charged to the tool and
    not to this driver (a subprocess.run + /proc poll would miss short peaks).
    """
    log = work / f"log_{tool}_{n}.txt"
    t0 = time.time()
    with open(log, "w") as lf:
        pid = os.fork()
        if pid == 0:
            os.dup2(lf.fileno(), 1)
            os.dup2(lf.fileno(), 2)
            os.chdir(work)
            try:
                os.execvp(cmd[0], cmd)
            except Exception:
                os._exit(127)
    _, status, ru = os.wait4(pid, 0)
    wall = time.time() - t0
    ok = os.WIFEXITED(status) and os.WEXITSTATUS(status) == 0
    return {
        "tool": tool,
        "n": n,
        "wall_s": round(wall, 3),
        "rss_mb": round(ru.ru_maxrss / 1024.0),   # Linux ru_maxrss is KB
        "ok": int(ok),
    }


def build_cmd(tool: str, n: int, args, work: Path) -> list[str]:
    """The command line that maps N queries with `tool`, single-threaded."""
    py = sys.executable
    if tool == "fastCDS":
        return [args.bin, "map", "--index", str(args.index),
                "--bed", str(work / f"q_{n}.bed"),
                "--out-dir", str(work / f"out_fastcds_{n}"),
                "--output", "coding", "--threads", "1"]
    if tool == "ensembldb":
        return [args.rscript, str(HERE / "run_ensembldb.R"),
                str(work / f"ids_{n}.txt"), str(args.ensdb)]
    if tool == "genomicfeatures":
        return [args.rscript, str(HERE / "run_gf_granges.R"),
                str(work / f"enst_{n}.txt"), str(work / "cds_by_tx.rds")]
    if tool == "geneplot":
        return [py, str(HERE / "run_geneplot.py"), str(n),
                str(work / f"ids_{n}.txt"), "--gff", str(args.gff),
                "--ipr", str(args.ipr), "--ensp-enst", str(args.ensp_enst)]
    if tool == "rest":
        return [py, str(HERE / "run_rest.py"), str(n), str(work / f"ids_{n}.txt")]
    raise SystemExit(
        f"run_matched_scaling.py: unknown tool {tool!r}. "
        f"Valid tools: {', '.join(LADDER)}")


def assemble(work: Path, out: Path) -> None:
    """Collapse results.tsv into the tidy per-tool curve the notebook plots."""
    res = work / "results.tsv"
    if not res.exists():
        raise SystemExit(
            f"run_matched_scaling.py: no measurements at {res}. "
            f"Run a sweep first (drop --assemble-only).")
    r = pd.read_csv(res, sep="\t")
    r = r[r.ok == 1]
    r = r[r.tool.isin(PAPER_TOOLS)]
    # Re-measured (tool, n) pairs: keep the most recent.
    r = r.drop_duplicates(["tool", "n"], keep="last").sort_values(["tool", "n"])
    out.parent.mkdir(parents=True, exist_ok=True)
    r[["tool", "n", "wall_s", "rss_mb"]].to_csv(out, sep="\t", index=False)
    print("wrote", out)
    print(r.pivot_table(index="n", columns="tool", values="wall_s").round(1).to_string())


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--work-dir", required=True, type=Path,
                    help="Directory holding the query ladder; results land here")
    ap.add_argument("--out", type=Path, default=None,
                    help="Assembled table (default: <work-dir>/../bench/scaling_matched.tsv)")
    ap.add_argument("--tools", nargs="+", default=PAPER_TOOLS,
                    help=f"Subset to run. Valid: {', '.join(LADDER)}")
    ap.add_argument("--sizes", nargs="+", type=int, default=None,
                    help="Override the per-tool N ladder")
    ap.add_argument("--assemble-only", action="store_true",
                    help="Skip measuring; just rebuild the tidy table")
    ap.add_argument("--bin", default="build/fastCDS", help="fastCDS binary")
    ap.add_argument("--index", type=Path, help="fastCDS index (.idx)")
    ap.add_argument("--ensdb", type=Path, help="EnsDb SQLite path (ensembldb route)")
    ap.add_argument("--rscript", default="Rscript")
    ap.add_argument("--gff", type=Path, help="GFF3 for geneplot")
    ap.add_argument("--ipr", type=Path, help="InterProScan .ipr for geneplot")
    ap.add_argument("--ensp-enst", type=Path, help="ENSP -> ENST map for geneplot")
    args = ap.parse_args()

    work = args.work_dir.expanduser().resolve()
    out = args.out or work.parent / "bench" / "scaling_matched.tsv"

    if not args.assemble_only:
        work.mkdir(parents=True, exist_ok=True)
        unknown = [t for t in args.tools if t not in LADDER]
        if unknown:
            raise SystemExit(
                f"run_matched_scaling.py: unknown tool(s) {', '.join(unknown)}. "
                f"Valid tools: {', '.join(LADDER)}")
        results = work / "results.tsv"
        if not results.exists():
            results.write_text("tool\tn\twall_s\trss_mb\tok\n")
        for tool in args.tools:
            for n in (args.sizes or LADDER[tool]):
                row = measure(tool, n, build_cmd(tool, n, args, work), work)
                with open(results, "a") as f:
                    f.write(f"{row['tool']}\t{row['n']}\t{row['wall_s']}\t"
                            f"{row['rss_mb']}\t{row['ok']}\n")
                print(f"  {tool:<16} n={n:<8} wall={row['wall_s']:>9.2f}s  "
                      f"rss={row['rss_mb']:>6}MB  ok={bool(row['ok'])}", flush=True)

    assemble(work, out)


if __name__ == "__main__":
    main()
