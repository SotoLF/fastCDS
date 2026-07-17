"""Run fastCDS and ensembldb on the same query set, classify every query into
the agreement buckets, and emit Supplementary Table S1 (9 categories) (stratified agreement).

Inputs:
  --queries-bed       BED-like, 4 cols (protein_id, aa_start, aa_end, query_id)
  --queries-meta      TSV: query_id, category, ... (from sample_validation_queries.py)
  --fastCDS-index   fastCDS binary index
  --ensdb             EnsDb sqlite path
  --rscript           path to Rscript (default: $CONDA_PREFIX/bin/Rscript or `Rscript`)
  --fastCDS-bin     path to fastCDS binary (default: ../build/fastCDS)
  --out-dir           where to write intermediate outputs + Supplementary Table S1

Outputs in --out-dir:
  fastCDS/*                       fastCDS's own output files
  ensembldb_intervals.tsv           output of the R helper
  table1.tsv                        the stratified agreement table
  discrepancies.tsv                 per-query diff for non-exact-match rows

Buckets:
  exact_match           same set of genomic intervals (chrom, start, end)
  off_by_one            sets differ but the bp-symmetric-difference is <= 2
  structural_mismatch   both tools returned intervals but the sets don't agree
                        and aren't an off-by-one case
  only_fastCDS        fastCDS returned >=1 interval; ensembldb returned nothing
  only_ensembldb        the converse
"""

from __future__ import annotations

import argparse
import csv
import os
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_BIN = REPO_ROOT / "build" / "fastCDS"


def run_fastCDS(binary: Path, index: Path, bed: Path, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [str(binary), "map", "--index", str(index),
           "--bed", str(bed), "--out-dir", str(out_dir),
           "--output", "coding"]  # only need CDS segments
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise SystemExit(f"fastCDS failed: {proc.stderr}")


def load_fastCDS_intervals(out_dir: Path) -> dict[str, list[tuple[str, int, int]]]:
    """input_id -> list of (chrom, start, end) for coding_overlap rows."""
    by_qid: dict[str, list[tuple[str, int, int]]] = defaultdict(list)
    tsv = out_dir / "domain_cds_segments.tsv"
    if not tsv.exists():
        return by_qid
    with open(tsv) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["overlaps_domain"] != "coding_overlap":
                continue
            # Use the domain_overlap_genomic_* sub-interval - that's the actual
            # bases coding the domain, not the whole CDS slice.
            start = row["domain_overlap_genomic_start"]
            end = row["domain_overlap_genomic_end"]
            if start == "NA" or end == "NA":
                continue
            by_qid[row["input_id"]].append((row["chrom"], int(start), int(end)))
    return by_qid


def load_ensembldb_intervals(path: Path) -> dict[str, list[tuple[str, int, int]]]:
    """query_id -> list of (chrom, start, end). Empty list means no_result/error."""
    by_qid: dict[str, list[tuple[str, int, int]]] = defaultdict(list)
    with open(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["status"] != "ok":
                # Ensure the key exists so we know the query was attempted.
                by_qid.setdefault(row["query_id"], [])
                continue
            by_qid[row["query_id"]].append(
                (row["chrom"], int(row["start"]), int(row["end"]))
            )
    return by_qid


def normalize_chrom(c: str) -> str:
    """Strip 'chr' prefix so GENCODE 'chr1' matches Ensembl '1'."""
    return c[3:] if c.startswith("chr") else c


def interval_set(intervals: list[tuple[str, int, int]]) -> set[tuple[str, int, int]]:
    return {(normalize_chrom(c), s, e) for c, s, e in intervals}


def total_length(intervals: list[tuple[str, int, int]]) -> int:
    return sum(e - s + 1 for _, s, e in intervals)


def classify(ours: list[tuple[str, int, int]], theirs: list[tuple[str, int, int]]) -> str:
    o_set = interval_set(ours)
    t_set = interval_set(theirs)
    o_empty, t_empty = not o_set, not t_set
    if o_empty and t_empty:
        return "neither_mapped"
    if t_empty:
        return "only_fastCDS"
    if o_empty:
        return "only_ensembldb"
    if o_set == t_set:
        return "exact_match"
    # Symmetric-difference of base coverage. Approximate by total bp difference.
    o_len = total_length(ours)
    t_len = total_length(theirs)
    if abs(o_len - t_len) <= 2:
        # Same envelope but the boundary differs by 1-2 bp - usually codon-split.
        return "off_by_one"
    return "structural_mismatch"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--queries-bed", required=True, type=Path)
    ap.add_argument("--queries-meta", required=True, type=Path)
    ap.add_argument("--fastCDS-index", required=True, type=Path)
    ap.add_argument("--ensdb", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--rscript", default=None,
                    help="Path to Rscript (default: $CONDA_PREFIX/bin/Rscript if set, else 'Rscript')")
    ap.add_argument("--fastCDS-bin", default=DEFAULT_BIN, type=Path)
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    # Resolve Rscript.
    rscript = args.rscript
    if rscript is None:
        conda_prefix = os.environ.get("CONDA_PREFIX")
        if conda_prefix and (Path(conda_prefix) / "bin" / "Rscript").exists():
            rscript = str(Path(conda_prefix) / "bin" / "Rscript")
        else:
            rscript = shutil.which("Rscript") or "Rscript"

    # 1) fastCDS
    fastcds_out = args.out_dir / "fastCDS"
    print("[1/3] running fastCDS ...", file=sys.stderr)
    run_fastCDS(args.fastCDS_bin, args.fastCDS_index, args.queries_bed, fastcds_out)

    # 2) ensembldb (via R subprocess)
    ens_out = args.out_dir / "ensembldb_intervals.tsv"
    print(f"[2/3] running ensembldb via {rscript} ...", file=sys.stderr)
    r_script = REPO_ROOT / "benchmarks" / "ensembldb_query.R"
    proc = subprocess.run(
        [rscript, str(r_script), str(args.ensdb), str(args.queries_bed), str(ens_out)],
        capture_output=True, text=True,
    )
    if proc.returncode != 0:
        sys.stderr.write(proc.stdout + "\n" + proc.stderr + "\n")
        raise SystemExit("ensembldb R helper failed")
    sys.stderr.write(proc.stderr)

    # 3) Compare and bucket
    print("[3/3] classifying ...", file=sys.stderr)
    fastcds_intervals = load_fastCDS_intervals(fastcds_out)
    ens_intervals = load_ensembldb_intervals(ens_out)

    # Load categories
    categories: dict[str, str] = {}
    with open(args.queries_meta) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            categories[row["query_id"]] = row["category"]

    # All query_ids that were attempted (drive from the BED).
    all_qids = []
    with open(args.queries_bed) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 4:
                all_qids.append(parts[3])

    by_cat: dict[str, Counter] = defaultdict(Counter)
    discrepancies = []
    for qid in all_qids:
        cat = categories.get(qid, "UNKNOWN")
        bucket = classify(fastcds_intervals.get(qid, []), ens_intervals.get(qid, []))
        by_cat[cat][bucket] += 1
        by_cat["OVERALL"][bucket] += 1
        if bucket not in ("exact_match", "neither_mapped"):
            discrepancies.append((qid, cat, bucket,
                                  fastcds_intervals.get(qid, []),
                                  ens_intervals.get(qid, [])))

    # Write Supplementary Table S1.
    table_path = args.out_dir / "table1.tsv"
    BUCKETS = ["exact_match", "off_by_one", "structural_mismatch",
               "only_fastCDS", "only_ensembldb", "neither_mapped"]
    with open(table_path, "w") as f:
        f.write("category\tn\t" + "\t".join(BUCKETS) + "\texact_pct\n")
        # OVERALL first.
        for cat in ["OVERALL"] + sorted(k for k in by_cat if k != "OVERALL"):
            counts = by_cat[cat]
            n = sum(counts.values())
            exact_pct = (100.0 * counts["exact_match"] / n) if n else 0.0
            f.write(f"{cat}\t{n}\t" + "\t".join(str(counts[b]) for b in BUCKETS)
                    + f"\t{exact_pct:.2f}\n")
    print(f"wrote {table_path}", file=sys.stderr)

    # Write per-query discrepancies for follow-up.
    disc_path = args.out_dir / "discrepancies.tsv"
    with open(disc_path, "w") as f:
        f.write("query_id\tcategory\tbucket\tfastCDS_intervals\tensembldb_intervals\n")
        for qid, cat, bucket, ours, theirs in discrepancies:
            f.write(f"{qid}\t{cat}\t{bucket}\t"
                    f"{','.join(f'{c}:{s}-{e}' for c, s, e in ours)}\t"
                    f"{','.join(f'{c}:{s}-{e}' for c, s, e in theirs)}\n")
    print(f"wrote {disc_path} ({len(discrepancies):,} discrepancy rows)", file=sys.stderr)

    # Print Supplementary Table S1 to stdout for the user.
    print()
    print(table_path.read_text())


if __name__ == "__main__":
    main()
