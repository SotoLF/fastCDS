# Benchmarks

> **Scope.** This page covers the head-to-head comparison against other protein-to-genome tools. For prot2exon's own perf knobs (`--threads`, `--batch-size`, RAM cap, single-tool scaling), see [[Performance and RAM]].

Prot2Exon, ensembldb, TransVar, Ensembl REST measured on the same 5,000-query stratified set. See [[Validation]] for the correctness side of the comparison.

## Table 1 — the 4-tool comparison

| Metric | prot2exon | ensembldb | TransVar | Ensembl REST |
|---|---|---|---|---|
| Exact agreement vs prot2exon | ref | 100.00 % (5,000 / 5,000) | 100.00 % (1,761 / 1,761) | 98.30 % (983 / 1,000) |
| Runtime @ N = 10,000 (1 thread) | **1.71 s** | 1,558 s | 7.54 s | rate-limited (~667 s @ 15 q/s cap; observed 9,180 s) |
| Peak RSS @ N = 10,000 | 788 MB | 1,252 MB | **284 MB** | N/A (HTTP client) |
| Throughput @ N = 10,000 (q/s) | **5,847** | 6 | 1,326 | 1.09 (network-bound) |
| Parallelism (OpenMP / threads) | Yes | No | No | N/A |
| Plot-ready output schema | Yes | No | No | No |
| Multi-species support | Yes (any GTF) | Yes (any Ensembl release) | Yes (hg19/hg38/mm9/mm10/…) | Yes (Ensembl-supported) |
| Largest N tested | **1,000,000** | 10,000 | 10,000 | 1,000 (rate cap) |
| Index / DB size on disk | 87 MB binary | 333 MB sqlite | 236 MB transvardb + 3 GB fasta | N/A (remote) |

**~900× speedup over ensembldb at N = 10,000** with a smaller index, **~4.4× faster than TransVar with no FASTA required**, identical genomic intervals against every tool that returned data.

The agreement denominators differ because each tool answers a slightly different question — we report them honestly:

- **ensembldb**: 100 % on a stratified set built from the same Ensembl 86 annotation it indexes. No annotation drift; strongest claim.
- **TransVar**: 100 % on the 1,761 queries whose ENST is present in TransVar's bundled annotation. The other 3,239 fall into `only_prot2exon` (annotation-drift story, same as the v113 EnsDb).
- **Ensembl REST**: 98.30 % on 1,000 queries (rate limit caps practical N). The 17 disagreements are all **off-by-one** (codon-split convention) — zero structural mismatches. The other tools agree with prot2exon, not REST, on those rows.

## Why these three comparators

Each represents a different use pattern users are choosing between:

| Tool | What it represents |
|---|---|
| **ensembldb** | The R/Bioconductor canonical; widely cited |
| **TransVar** | Variant-annotation perspective, HGVS-based; popular with clinical/variant teams |
| **Ensembl REST** | The "no install" zero-overhead path |

Others we considered and rejected: GeneMANIA (no per-domain output), peptidomics tools (different problem), VEP (genome→protein, the opposite direction).

## Raw scaling (single thread)

| N | prot2exon wall (median, s) | prot2exon RSS (MB) | ensembldb wall (s) | ensembldb RSS (MB) |
|---:|---:|---:|---:|---:|
| 100 | 1.33 | 659 | 23.3 | 979 |
| 1,000 | 1.37 | 671 | 168.8 | 988 |
| 10,000 | 1.71 | 788 | 1,558.2 | 1,252 |
| 100,000 | 4.66 | 1,959 | (skipped) | — |
| 1,000,000 | 129.4 | 11,045 | (skipped) | — |

ensembldb was capped at N = 10K — linear extrapolation puts its N = 100K at ~4.3 h, N = 1M at ~43 h. Continuing past 10K would have contributed nothing beyond time burned.

The prot2exon side at small N is dominated by **~1.3 s of one-time index load**. The actual mapping work for the first 1,000 queries is essentially free; wall time at N = 100 ≈ N = 1,000. That's why "throughput per query" is misleading at small N — it counts index load against per-query work.

## Parallel scaling (N = 1,000,000)

| Threads | Wall (median, s) | Speedup | Efficiency |
|---:|---:|---:|---:|
| 1 | 107.5 | 1.00 | 1.00 |
| 2 | 65.0 | 1.65 | 0.83 |
| 4 | 71.4 | 1.51 | 0.38 |
| 8 | 73.9 | 1.46 | 0.18 |

Mapping itself is OpenMP-parallel; the seven output files are written concurrently via `parallel sections`; **stdio write buffers are bumped to 1 MiB per file** (`pubsetbuf`) — the latter cut single-thread wall time by ~17 % vs the default 8 KB buffer. Past 4 threads the bottleneck is OS page-cache flushing + raw disk bandwidth, not formatting.

For RAM-bounded large-N runs, see `--batch-size` on [[Performance and RAM]] — the same 1M benchmark shows 11× peak-RSS reduction with `--batch-size 10000`.

## What's in `benchmarks/`

```
benchmarks/sample_validation_queries.py   # GTF parser + 9-stratum sampler
benchmarks/ensembldb_query.R              # Batch proteinToGenome via R subprocess
benchmarks/validate_vs_ensembldb.py       # Validator — runs both, classifies, emits Table 1
benchmarks/scaling_benchmark.py           # prot2exon vs ensembldb at N=100..1M
benchmarks/parallel_benchmark.py          # prot2exon at threads 1, 2, 4, 8
benchmarks/run_ensembl_rest.py            # Rate-limited REST client
benchmarks/run_transvar.py                # Builds HGVS from EnsDb seqs, drives TransVar
benchmarks/classify_external.py           # Bucket-classifier with --envelope-only
benchmarks/make_scaling_outputs.py        # 2-tool Table 1 + scaling.png
benchmarks/make_table_1.py                # 4-tool Table 1
benchmarks/make_figure_1.py               # 4-panel composite (paper Figure 1)
```

## Reproducing the scaling table

```bash
conda activate prot2exon-val

python benchmarks/scaling_benchmark.py \
    --bin build/prot2exon \
    --p2e-index ~/Desktop/protein2genomic_data/human.idx \
    --ensdb $CONDA_PREFIX/lib/R/library/EnsDb.Hsapiens.v86/extdata/EnsDb.Hsapiens.v86.sqlite \
    --rscript $CONDA_PREFIX/bin/Rscript \
    --r-helper benchmarks/ensembldb_query.R \
    --source-bed queries.bed \
    --work-dir bench \
    --sizes 100 1000 10000 100000 1000000 \
    --p2e-reps 2 --ensembldb-reps 1 --ensembldb-max-n 10000 \
    --out bench/timings.tsv

python benchmarks/parallel_benchmark.py \
    --bin build/prot2exon \
    --p2e-index ~/Desktop/protein2genomic_data/human.idx \
    --bed bench/queries_n100000.bed \
    --work-dir bench/parallel \
    --out bench/parallel.tsv
```

Total wall: scaling ~50 min (ensembldb N = 10K alone is ~26 min), parallel ~1 min.

## Notes on the external comparators

- **Ensembl REST is network-bound, not rate-limited.** The 15 q/s cap isn't the bottleneck — per-request HTTP RTT is (~900 ms). At 1.09 q/s effective, a concurrent client (HTTP keep-alive + async pool) could push closer to the 15 cap, but at the cost of moving away from how anyone actually writes a REST script. We measured the canonical pattern, not the speed-of-light pattern.
- **TransVar measures envelope only.** Its `coordinates` field is one `chrN:g.start_end` per query — the genomic envelope, introns included, not per-CDS intervals. We compare via `classify_external.py --envelope-only` (collapse both sides to `(chrom, min_start, max_end)` before set comparison). Apples-to-apples. prot2exon retains the per-exon decomposition, the harder computation; TransVar is being asked an easier question.
- **TransVar throws away `protein_id` on ENSP input — keys on ENST.** Feeding it ENSPs returns silently empty results. We feed it ENSTs via the `queries_meta.tsv` `transcript_id` column.
- **`transvar config` is interactive on first run.** Pipe `echo "" | transvar config …` to satisfy the FASTA prompt.
