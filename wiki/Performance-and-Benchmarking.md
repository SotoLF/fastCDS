# Performance and benchmarking

This page covers prot2exon's own performance characteristics (threading, memory, scaling), its coordinate accuracy validated against ensembldb, and head-to-head speed comparisons against the other common protein-to-genome tools. For how to run a mapping job see [[Mapping]]; for how the binary index is built see [[Index]].

## Speed and memory usage

All single-tool numbers below are measured against the pre-built human GENCODE index (3.17 GB GTF → ~298 MB binary index). The index encodes per-protein exon and CDS intervals plus `gene_id`, `gene_name`, and `exon_number`; build it once with `prot2exon index` and reuse it until you upgrade the annotation (see [[Index]]). Reproduce the measurements with the scripts in `tutorial/reproduce_paper/benchmarks/` ([`run_benchmark.py`](https://github.com/SotoLF/Prot2Exon/blob/main/tutorial/reproduce_paper/benchmarks/run_benchmark.py), [`scaling_benchmark.py`](https://github.com/SotoLF/Prot2Exon/blob/main/tutorial/reproduce_paper/benchmarks/scaling_benchmark.py), [`parallel_benchmark.py`](https://github.com/SotoLF/Prot2Exon/blob/main/tutorial/reproduce_paper/benchmarks/parallel_benchmark.py)) or the [`scaling_and_ram.ipynb`](https://github.com/SotoLF/Prot2Exon/blob/main/tutorial/reproduce_paper/end_to_end/scaling_and_ram.ipynb) notebook.

### Single-threaded query performance

Median of 3 reps, single core:

| Queries | Wall (all output) | Peak RSS | Wall (isoform output) | Peak RSS |
|---|---|---|---|---|
| 100 | 1.7 s | 847 MB | 1.6 s | 847 MB |
| 1,000 | 2.1 s | 858 MB | 1.7 s | 859 MB |
| 10,000 | 5.5 s | 974 MB | 3.3 s | 974 MB |
| 100,000 | 37 s | 2,126 MB | 17.6 s | 2,125 MB |

Effective throughput at 100k is roughly 2,800 q/s producing all outputs and 5,700 q/s for a single isoform TSV; the difference is that the full-output path writes four TSVs and three BEDs versus one TSV. The constant ~1.5 s floor at small N is index loading, per-query cost is roughly linear and dominated by output formatting, and memory grows because output strings accumulate before being written.

### Parallelism and memory: `--threads` × `--batch-size`

Two knobs, two independent axes — one trades wall time, the other trades peak
RAM, so we sweep them together rather than in separate tables:

- **`--threads`** runs per-query processing in an OpenMP loop (and the per-file
  writes via OpenMP `sections`). It cuts **wall time**.
- **`--batch-size N`** streams results to disk in chunks of `N` and frees each
  chunk before the next (`0` / omitted = one-shot, hold everything in memory).
  It caps **peak RSS** at `O(N × per-query result size)`, byte-identical output.

**Wall time (s)** at each (threads, batch) cell — fastest of 2 reps,
**N = 1,000,000** queries, `--output coding`, on a quiet 32-core / 125 GB
workstation (enough RAM that nothing swaps, so this is the *pure* tradeoff):

| Threads | one-shot | `--batch-size 100000` | `--batch-size 50000` | `--batch-size 10000` |
|---:|---:|---:|---:|---:|
| 1  | 23.7 | 18.9 | 18.6 | 18.8 |
| 4  | 12.7 | 13.1 | 13.2 | 13.5 |
| 8  | 11.3 | 12.2 | 12.3 | 12.4 |
| 16 | **10.8** | 12.3 | 12.2 | 12.3 |
| 32 | **10.8** | 12.2 | 12.0 | 12.2 |

**Peak RSS depends only on the batch size, not the thread count** (it's flat down
every column): one-shot **13.3 GB**, `--batch-size 100000 / 50000 / 10000` =
**2.4 / 1.6 / 0.9 GB**. So the two knobs really are orthogonal — `--threads` for
wall time, `--batch-size` for memory.

Reading the grid down a column: **`--threads` cuts wall time ~2.2×** by 16
threads (23.7 → 10.8 s one-shot), then plateaus as the single-threaded TSV writer
and memory bandwidth take over — set `--threads` to your physical core count.
Across a row: once you have **≥ 4 threads, bigger batches are slightly faster**
(closer to one-shot, since fewer flush cycles) but use proportionally more RAM —
a genuine speed/memory dial. The one exception is **single-threaded**, where
one-shot is the *slowest* cell (23.7 s): holding all 1 M results (13.3 GB) in
memory thrashes the allocator and cache, and any batch fixes it.

**Rule of thumb:** one-shot is fastest when you have spare cores *and* spare RAM;
otherwise `--batch-size 10000` gives ~14× less memory (0.9 vs 13.3 GB) for
~10–15 % more wall, with byte-identical output. Batch when N is large or RAM is
tight.

```bash
prot2exon map --index human.idx --bed q.bed --out-dir out --threads 8 --batch-size 10000
```

Reproduce the grid with [`tutorial/reproduce_paper/benchmarks/threads_batch_grid.py`](https://github.com/SotoLF/Prot2Exon/blob/main/tutorial/reproduce_paper/benchmarks/threads_batch_grid.py).

## Accuracy vs other tools

Coordinate correctness is validated against ensembldb, the Bioconductor canonical for protein-to-genome mapping (~800 paper citations). Because ensembldb is an independent R/SQL implementation on top of EnsDb, an agreement is genuine cross-validation rather than testing the same code twice. The validator drives the prot2exon binary and shells out to `Rscript` for ensembldb; reproduce it with [`validate_vs_ensembldb.py`](https://github.com/SotoLF/Prot2Exon/blob/main/tutorial/reproduce_paper/benchmarks/validate_vs_ensembldb.py) and the [`validation.ipynb`](https://github.com/SotoLF/Prot2Exon/blob/main/tutorial/reproduce_paper/end_to_end/validation.ipynb) notebook.

The headline result is 100.00% exact match against ensembldb on a 5,000-query stratified set — zero off-by-ones and zero structural mismatches. Random sampling would underweight the corner cases that matter, so a 9-stratum sampler ([`sample_validation_queries.py`](https://github.com/SotoLF/Prot2Exon/blob/main/tutorial/reproduce_paper/benchmarks/sample_validation_queries.py)) ensures every condition that historically breaks these tools is represented:

| Stratum | n | What it stresses |
|---|---:|---|
| `single_exon_domain` | 1,000 | The common shape (high test coverage) |
| `multi_exon_domain` | 1,000 | The hard case — find all CDS pieces |
| `codon_split_boundary` | 500 | Codons straddling exon boundaries, where off-by-ones hide |
| `plus_strand_gene` | 1,000 | Strand-handling A/B |
| `minus_strand_gene` | 1,000 | Strand-handling A/B (most bugs live here) |
| `cds_incomplete` | 200 | `cds_start_NF` / `cds_end_NF` transcripts |
| `selenoprotein` | 100 | 25-gene curated list (UGA → Sec recoding) |
| `single_exon_gene` | 100 | Boundary case for the multi-exon path |
| `many_exon_gene` | 100 | > 20 CDS exons — exon-walker stress |

Every stratum returns 100% exact agreement on the matched-annotation (v86) path:

```
category               n     exact  off_by_one  structural  only_p2e  only_ens  exact_pct
OVERALL              5000   5000          0           0         0         0    100.00
cds_incomplete        200    200          0           0         0         0    100.00
codon_split_boundary  500    500          0           0         0         0    100.00
many_exon_gene        100    100          0           0         0         0    100.00
minus_strand_gene    1000   1000          0           0         0         0    100.00
multi_exon_domain    1000   1000          0           0         0         0    100.00
plus_strand_gene     1000   1000          0           0         0         0    100.00
selenoprotein         100    100          0           0         0         0    100.00
single_exon_domain   1000   1000          0           0         0         0    100.00
single_exon_gene      100    100          0           0         0         0    100.00
```

The agreement is bucketed by whether both tools return the same `(chrom, start, end)` set (`exact_match`), differ by ≤ 2 bp under the codon-split convention (`off_by_one`), disagree in a way that is neither (`structural_mismatch`), or return data from only one side (`only_prot2exon` / `only_ensembldb`); rows where both tools return nothing are excluded from the denominator.

The cross-validation also extends to **`GenomicFeatures::proteinToGenome`**, a second independent Bioconductor implementation (its GRangesList method maps against a CDS-by-transcript `GRangesList` held in memory, rather than querying the EnsDb SQLite). On a 1,000-query v86 subset all three tools agree on every query — prot2exon ≡ ensembldb ≡ GenomicFeatures at **1,000 / 1,000 exact-segment match**, pairwise. Reproduce with [`proteintogenome_bench.R`](https://github.com/SotoLF/Prot2Exon/blob/main/tutorial/reproduce_paper/benchmarks/proteintogenome_bench.R) and [`compare_intervals.py`](https://github.com/SotoLF/Prot2Exon/blob/main/tutorial/reproduce_paper/benchmarks/compare_intervals.py).

A second run tests annotation drift, comparing prot2exon on a current GENCODE GTF against an Ensembl 113 EnsDb pulled via AnnotationHub. There it is 100% exact match on the 2,264 / 5,000 queries where both tools return data. The remaining 2,736 fall into `only_prot2exon` because those ENSPs exist in the v113 EnsDb's protein table but the EnsDb has no CDS linkage for them, so `proteinToGenome()` reports "No CDS found." That is a gap in the EnsDb rather than a prot2exon disagreement, confirmed by running prot2exon on the same GTF the ENSPs came from and getting complete answers.

A few practical points matter when reproducing this. Installing ensembldb is best done through BiocManager from CRAN/Bioconductor rather than bioconda, whose `r-rjson` recipe pins an ancient `r=3.3.1` that conflicts with modern r-base; the R compile also needs `gcc ≤ 14` (gcc 15+ rejects `rtracklayer`'s bundled UCSC C code) and the development `libxml2` headers (`libxml2=2.13.9`, not the runtime-only conda-forge 2.15). The sampler also handles both GTF dialects — GENCODE uses `transcript_type` and puts `protein_id` on transcript rows, while Ensembl uses `transcript_biotype` and puts `protein_id` only on CDS rows — and skipping that distinction silently halves the number of usable queries.

## Speed vs other tools

The same 5,000-query stratified set drives a head-to-head against the three tools users most often choose between: ensembldb (the R/Bioconductor canonical), TransVar (the HGVS-based variant-annotation perspective, popular with clinical teams), and Ensembl REST (the no-install zero-overhead path). Tools considered and rejected include GeneMANIA (no per-domain output), peptidomics tools (a different problem), and VEP (genome→protein, the opposite direction). Reproduce with [`scaling_benchmark.py`](https://github.com/SotoLF/Prot2Exon/blob/main/tutorial/reproduce_paper/benchmarks/scaling_benchmark.py), [`run_transvar.py`](https://github.com/SotoLF/Prot2Exon/blob/main/tutorial/reproduce_paper/benchmarks/run_transvar.py), [`run_ensembl_rest.py`](https://github.com/SotoLF/Prot2Exon/blob/main/tutorial/reproduce_paper/benchmarks/run_ensembl_rest.py), and the [`software_comparison.ipynb`](https://github.com/SotoLF/Prot2Exon/blob/main/tutorial/reproduce_paper/end_to_end/software_comparison.ipynb) notebook.

> **The one number to quote: ~970× faster than ensembldb.** It is measured
> **end-to-end at N = 10,000, single thread** — total wall time from process
> start until all results are written, *including the one-time index load*,
> divided by N. That's the number a user actually waits for, and a fair one (it
> charges prot2exon for its own load). The only caveat is that you must **state
> N**, because the ratio grows with N: prot2exon's ~1.2 s index load is amortized
> over more queries, so prot2exon-vs-ensembldb is ~130× at N = 1,000 and ~970× at
> N = 10,000. We use N = 10,000 everywhere as the headline.

| Metric | prot2exon | ensembldb | TransVar | Ensembl REST |
|---|---|---|---|---|
| Exact agreement vs prot2exon | ref | 100.00% (5,000 / 5,000) | 100.00% (1,761 / 1,761) | 98.30% (983 / 1,000) |
| Runtime @ N = 10,000 (1 thread) | 1.71 s | 1,558 s | 7.54 s | rate-limited (~667 s @ 15 q/s cap; observed 9,180 s) |
| Peak RSS @ N = 10,000 | 788 MB | 1,252 MB | 284 MB | N/A (HTTP client) |
| Throughput @ N = 10,000 (q/s) | 5,847 | 6 | 1,326 | 1.09 (network-bound) |
| Parallelism (OpenMP / threads) | Yes | No | No | N/A |
| Plot-ready output schema | Yes | No | No | No |
| Multi-species support | Yes (any GTF) | Yes (any Ensembl release) | Yes (hg19/hg38/mm9/mm10/…) | Yes (Ensembl-supported) |
| Largest N tested | 1,000,000 | 10,000 | 10,000 | 1,000 (rate cap) |
| Index / DB size on disk | 87 MB binary | 333 MB sqlite | 236 MB transvardb + 3 GB fasta | N/A (remote) |

That is a ~970× end-to-end speedup over ensembldb at N = 10,000 (5,847 vs 6 q/s) with a smaller index, about 4.4× faster than TransVar with no FASTA required, and identical genomic intervals against every tool that returned data. The agreement denominators differ because each tool answers a slightly different question: ensembldb is 100% on a stratified set built from the same annotation it indexes (no drift, the strongest claim); TransVar is 100% on the 1,761 queries whose ENST is in its bundled annotation, with the other 3,239 falling into `only_prot2exon` from annotation drift; and Ensembl REST is 98.30% on the 1,000 queries the rate limit allows, where all 17 disagreements are off-by-one under the codon-split convention with zero structural mismatches, and the other tools agree with prot2exon rather than REST on those rows.

Raw single-thread scaling:

| N | prot2exon wall (median, s) | prot2exon RSS (MB) | ensembldb wall (s) | ensembldb RSS (MB) |
|---:|---:|---:|---:|---:|
| 100 | 1.33 | 659 | 23.3 | 979 |
| 1,000 | 1.37 | 671 | 168.8 | 988 |
| 10,000 | 1.71 | 788 | 1,558.2 | 1,252 |
| 100,000 | 4.66 | 1,959 | (skipped) | — |
| 1,000,000 | 129.4 | 11,045 | (skipped) | — |

ensembldb was capped at N = 10K — linear extrapolation puts its N = 100K at ~4.3 h and N = 1M at ~43 h, so continuing past 10K would have contributed nothing beyond time burned. The prot2exon side at small N is dominated by the ~1.3 s one-time index load, so the actual mapping work for the first 1,000 queries is essentially free and wall time at N = 100 ≈ N = 1,000; this is why per-query throughput is misleading at small N, since it counts index load against per-query work.

The parallel scaling matches the numbers in the **Parallelism and memory** section above (~1.4–1.65× by 2–4 threads, plateauing past 4 as page-cache flushing and disk bandwidth dominate). For RAM-bounded large-N runs, the same 1M benchmark shows the large peak-RSS reduction from `--batch-size 10000` described there.

A few notes on the external comparators. Ensembl REST is network-bound rather than rate-limited — the 15 q/s cap is not the bottleneck, the ~900 ms per-request HTTP round-trip is, and a concurrent keep-alive client could push closer to the cap but at the cost of measuring something other than how anyone actually writes a REST script. TransVar reports only the genomic envelope (one `chrN:g.start_end` per query, introns included) rather than per-CDS intervals, so the comparison collapses both sides to `(chrom, min_start, max_end)` via [`classify_external.py`](https://github.com/SotoLF/Prot2Exon/blob/main/tutorial/reproduce_paper/benchmarks/classify_external.py) `--envelope-only` for an apples-to-apples result, meaning prot2exon is doing the harder per-exon decomposition while TransVar answers an easier question. TransVar also keys on ENST and silently returns empty for ENSP input, so it is fed ENSTs from the query metadata, and `transvar config` is interactive on first run (pipe an empty string to satisfy the FASTA prompt).

### GenomicFeatures::proteinToGenome — the GRanges path

`proteinToGenome` exists in **two** Bioconductor packages: ensembldb (EnsDb/SQLite-backed) and GenomicFeatures (a GRangesList method that maps in memory against a CDS-by-transcript `GRangesList`). The GenomicFeatures path trades a one-time setup (`cdsBy(edb, "tx")` + a protein→transcript map) for SQLite-free mapping, so it is meaningfully faster than ensembldb, but it is still an R-level per-query loop, far short of the indexed C++ path. On 1,000 v86 queries, one machine, one thread (setup + map = end-to-end total):

| Tool | Setup (load) | Map (1,000) | Total (end-to-end) | Peak RSS |
|---|---:|---:|---:|---:|
| **prot2exon** | 1.19 s (index load) | **0.014 s** | **1.20 s** | 674 MB |
| GenomicFeatures::proteinToGenome | 5.99 s (`cdsBy`) | 37.71 s | 43.70 s | 1,364 MB |
| ensembldb::proteinToGenome | n/a | 160.12 s | 160.12 s | 1,163 MB |

At this N = 1,000 the end-to-end speedup is ~36× over GenomicFeatures (1.20 s vs 43.70 s) — smaller than the N = 10,000 headline only because prot2exon's 1.19 s load isn't amortized yet. GenomicFeatures is ~3.7× faster than ensembldb end-to-end, confirming the GRanges representation is the lighter of the two, and all three return identical coordinates on the 1,000-query set (see Accuracy above). One subtlety when reproducing: loading an EnsDb pulls in ensembldb's own `proteinToGenome` GRangesList method (which requires protein-sequence metadata on the CDS), so the runner fetches GenomicFeatures' method explicitly via `getMethod(..., where = asNamespace("GenomicFeatures"))`.

### Other domain-on-gene tools (geneplot, VisProDom)

Two other tools that "show domains on gene structure" are sometimes mentioned alongside these. Both are plot-oriented rather than proteome-scale mappers, so they sit outside the tables above, but tested here both do map coordinates internally:

- **geneplot** (Python; `gffutils` + BioPython) draws one gene at a time from a GFF3 plus an InterProScan domain file, and it *does* map internally: its `_transcriptpos_to_genomepos()` walks the gene's CDS to turn protein coordinates into genomic ones. Run on the human Ensembl-86 set, the end-to-end cost is ~14 genes/s: it first builds a `gffutils` SQLite database of the whole genome (minutes), then for each gene queries that database and re-reads the InterProScan domain file. It is built for single-gene figures, not proteome-scale batches.
- **VisProDom** (R/Shiny) is a domain viewer, but its `CreDat(gff, annofile)` function *is* a real batch mapper (pure R/dplyr, cumulative-CDS arithmetic — the same idea prot2exon implements in C++). It has no prebuilt index, so it is **O(genome) per call**: it rebuilds every transcript's CDS layout on each invocation. On its bundled maize example it finishes in seconds, but on the human Ensembl-86 set the genome rebuild dominates: mapping the 1,000-query set takes **42.7 s at 2.3 GB**, since `CreDat` rebuilds every transcript's CDS layout on every call regardless of how many domains you ask for. That batch-recompute model is exactly what prot2exon's persistent binary index removes. Reproduce with [`visprodom_bench.R`](https://github.com/SotoLF/Prot2Exon/blob/main/tutorial/reproduce_paper/benchmarks/visprodom_bench.R).
