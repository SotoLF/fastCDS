"""Build the 4-tool paper Table 1 (prot2exon, ensembldb, TransVar, Ensembl REST).

Reads the measurements collected by the various Phase 3 scripts plus the
agreement tables produced by classify_external.py, and emits a single combined
table in the shape PLAN.txt specifies.

Hardcodes a handful of "measured-on-this-machine" numbers from the previous
runs so the script is reproducible from cached artifacts alone.
"""

from __future__ import annotations

import argparse
import csv
import statistics
from pathlib import Path


def median_by(rows, *keys, value="wall_s"):
    buckets: dict = {}
    for r in rows:
        k = tuple(r[x] for x in keys)
        buckets.setdefault(k, []).append(r[value])
    return {k: statistics.median(v) for k, v in buckets.items()}


def read_tsv(path: Path) -> list[dict]:
    rows = []
    with open(path) as f:
        # Skip blank and comment lines before reading the header.
        header = None
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            header = line.rstrip("\n").split("\t")
            break
        if header is None:
            return rows
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            row = {}
            for k, v in zip(header, parts):
                try:
                    row[k] = float(v) if "." in v else int(v)
                except ValueError:
                    row[k] = v
            rows.append(row)
    return rows


def parse_agreement(table_path: Path) -> dict:
    """Return {n_considered, exact, off_by_one, structural, only_p2e, only_ext}."""
    rows = read_tsv(table_path)
    overall = next((r for r in rows if r["category"] == "OVERALL"), None)
    if not overall:
        return {}
    return overall


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scaling-tsv", required=True, type=Path)
    ap.add_argument("--p2e-index-size", type=int, required=True)
    ap.add_argument("--ensdb-size", type=int, required=True)
    ap.add_argument("--ensembldb-agreement-table", required=True, type=Path,
                    help="Phase 2 v86 validation table1.tsv")
    ap.add_argument("--transvar-agreement-table", required=True, type=Path)
    ap.add_argument("--rest-agreement-table", required=True, type=Path)
    ap.add_argument("--transvar-wall-s", type=float, required=True, help="Wall sec at N=10k")
    ap.add_argument("--transvar-rss-mb", type=int, required=True)
    ap.add_argument("--rest-wall-s", type=float, required=True,
                    help="Wall sec for the REST run (1k queries, rate-limited)")
    ap.add_argument("--rest-n", type=int, default=1000)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    scaling = read_tsv(args.scaling_tsv)
    med_wall = median_by(scaling, "tool", "n", value="wall_s")
    med_rss = median_by(scaling, "tool", "n", value="peak_rss_mb")
    p2e_wall = med_wall.get(("prot2exon", 10000), float("nan"))
    p2e_rss = med_rss.get(("prot2exon", 10000), 0)
    ens_wall = med_wall.get(("ensembldb", 10000), float("nan"))
    ens_rss = med_rss.get(("ensembldb", 10000), 0)

    ensembldb_ag = parse_agreement(args.ensembldb_agreement_table)
    transvar_ag = parse_agreement(args.transvar_agreement_table)
    rest_ag = parse_agreement(args.rest_agreement_table)

    def agreement_pct(ag: dict, denominator_strategy: str = "non_neither") -> str:
        """Format as 'XX.XX% (n/m)' where m is the relevant denominator."""
        if not ag:
            return "—"
        exact = ag["exact_match"]
        # Consider only queries where the external tool returned data (omit only_prot2exon).
        # That's the relevant question for the paper: when both tools answer, do they agree?
        considered = ag["n"] - ag["only_prot2exon"] - ag["neither_mapped"]
        if considered == 0:
            return "0 / 0"
        pct = 100.0 * exact / considered
        return f"{pct:.2f}% ({exact:,}/{considered:,})"

    rows = [
        ("Tool",                          "prot2exon",                            "ensembldb",                                  "TransVar",                                   "Ensembl REST"),
        ("Exact agreement vs prot2exon",  "ref",                                  agreement_pct(ensembldb_ag),                  agreement_pct(transvar_ag),                   agreement_pct(rest_ag)),
        ("Runtime N=10,000 (1 thread)",   f"{p2e_wall:.2f} s",                    f"{ens_wall:.0f} s",                           f"{args.transvar_wall_s:.2f} s",              f"rate-limited (~{10000/15:.0f}s @ 15 q/s cap)"),
        ("Peak RSS @ N=10,000",           f"{p2e_rss} MB",                         f"{ens_rss} MB",                               f"{args.transvar_rss_mb} MB",                  "N/A (HTTP client)"),
        ("Throughput @ N=10,000 (q/s)",   f"{int(10000/p2e_wall):,}",              f"{10000/ens_wall:.0f}" if ens_wall else "—",  f"{int(10000/args.transvar_wall_s):,}",        f"{args.rest_n/args.rest_wall_s:.2f} (network-bound)"),
        ("Parallelism (OpenMP / threads)","Yes",                                   "No",                                          "No",                                          "N/A"),
        ("Plot-ready output schema",      "Yes",                                   "No",                                          "No",                                          "No"),
        ("Multi-species support",         "Yes (any GTF)",                         "Yes (any Ensembl release)",                   "Yes (hg19/hg38/mm9/mm10/etc.)",               "Yes (Ensembl-supported species)"),
        ("Largest N tested here",         f"{max(n for (t,n) in med_wall if t=='prot2exon'):,}",
                                          f"{max(n for (t,n) in med_wall if t=='ensembldb'):,}",
                                          "10,000",
                                          f"{args.rest_n:,} (rate cap)"),
        ("Index/DB size on disk",         f"{args.p2e_index_size/(1024*1024):.0f} MB binary",
                                          f"{args.ensdb_size/(1024*1024):.0f} MB sqlite",
                                          "236 MB (transvardb + 3 GB fasta)",
                                          "N/A (remote)"),
    ]

    with open(args.out, "w") as f:
        for row in rows:
            f.write("\t".join(row) + "\n")
    print(f"wrote {args.out}\n")
    # Pretty-print for stdout.
    widths = [max(len(row[i]) for row in rows) for i in range(len(rows[0]))]
    for row in rows:
        print(" | ".join(c.ljust(w) for c, w in zip(row, widths)))


if __name__ == "__main__":
    main()
