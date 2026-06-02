# Performance and benchmarking

This page covers prot2exon's own performance characteristics (threading, memory, scaling), its coordinate accuracy validated against ensembldb, and head-to-head speed comparisons against the other common protein-to-genome tools. For how to run a mapping job see [[Mapping]]; for how the binary index is built see [[Index]].

## Speed and memory usage

All single-tool numbers below are measured against the pre-built human GENCODE index (3.17 GB GTF → ~298 MB binary index). The index encodes per-protein exon and CDS intervals plus `gene_id`, `gene_name`, and `exon_number`; build it once with `prot2exon index` and reuse it until you upgrade the annotation (see [[Index]]). Reproduce the measurements with the scripts in `benchmarks/` ([`run_benchmark.py`](https://github.com/SotoLF/Prot2Exon/blob/main/benchmarks/run_benchmark.py), [`scaling_benchmark.py`](https://github.com/SotoLF/Prot2Exon/blob/main/benchmarks/scaling_benchmark.py), [`parallel_benchmark.py`](https://github.com/SotoLF/Prot2Exon/blob/main/benchmarks/parallel_benchmark.py)) or the [`scaling_and_ram.ipynb`](https://github.com/SotoLF/Prot2Exon/blob/main/notebooks/scaling_and_ram.ipynb) notebook.

### Single-threaded query performance

Median of 3 reps, single core:

| Queries | Wall (all output) | Peak RSS | Wall (isoform output) | Peak RSS |
|---|---|---|---|---|
| 100 | 1.7 s | 847 MB | 1.6 s | 847 MB |
| 1,000 | 2.1 s | 858 MB | 1.7 s | 859 MB |
| 10,000 | 5.5 s | 974 MB | 3.3 s | 974 MB |
| 100,000 | 37 s | 2,126 MB | 17.6 s | 2,125 MB |

Effective throughput at 100k is roughly 2,800 q/s producing all outputs and 5,700 q/s for a single isoform TSV; the difference is that the full-output path writes four TSVs and three BEDs versus one TSV. The constant ~1.5 s floor at small N is index loading, per-query cost is roughly linear and dominated by output formatting, and memory grows because output strings accumulate before being written.

### `--threads` (OpenMP)

Per-query processing runs in an OpenMP parallel loop, and the output writes for the seven files run concurrently via OpenMP `sections`. Peak RSS is largely unaffected by the thread count.

```bash
prot2exon map --index human.idx --bed q.bed --out-dir out --threads 8
```

Parallel scaling at N = 1,000,000:

| Threads | Wall (median, s) | Speedup | Efficiency |
|---:|---:|---:|---:|
| 1 | 107.5 | 1.00 | 1.00 |
| 2 | 65.0 | 1.65 | 0.83 |
| 4 | 71.4 | 1.51 | 0.38 |
| 8 | 73.9 | 1.46 | 0.18 |

Mapping scales well to 2 threads (0.83 efficiency), then speedup plateaus past 4 as the OS page cache and disk bandwidth take over from per-row formatting as the bottleneck. The writer uses 1 MiB stdio buffers per file (`pubsetbuf`, about 100× fewer syscalls than the default 8 KB), which cut single-threaded wall time by roughly 17% over the legacy buffer size. Pushing past the 4-thread plateau would require a different layer entirely (direct `io_uring`, on-the-fly compression, or smaller TSV output). Set `--threads` to your physical core count rather than the hyperthread count.

### `--batch-size` (RAM cap)

By default the binary holds every query's full result set in memory before writing, which is fast when it fits but at roughly 10 KB per query a 1M-query run needs about 10 GB. Pass `--batch-size N` to stream results to disk in chunks of N and free each chunk before the next, making peak RAM `O(N × per-query result size)`. Outputs are byte-identical to the one-shot path; batching is purely a memory optimisation.

1M-query benchmark (`--threads 4`, full output, NVMe SSD, 16 GB RAM machine):

| Mode | Wall | Peak RSS | Notes |
|---|---:|---:|---|
| `--batch-size 0` (one-shot) | 122.5 s | 10.4 GB | swap-bound on this box; faster with 32 GB+ |
| `--batch-size 10000` | 58.5 s | 961 MB | byte-identical outputs to one-shot |

On the 16 GB machine, batching cuts peak RAM about 11× and wall time about 2× because the one-shot run is swap-bound. On a machine with enough RAM to avoid swap, one-shot is slightly faster, so choose `--batch-size 0` for maximum throughput when RAM allows and `--batch-size N` to cap the working set. `--batch-size 10000` is a sensible default for 100k+ queries, and on HDDs the write phase becomes I/O-bound but `--threads` still helps because the writes are parallelised.

```bash
prot2exon map --index human.idx --bed q.bed --out-dir out --threads 4 --batch-size 10000
```

## Accuracy vs other tools

Coordinate correctness is validated against ensembldb, the Bioconductor canonical for protein-to-genome mapping (~800 paper citations). Because ensembldb is an independent R/SQL implementation on top of EnsDb, an agreement is genuine cross-validation rather than testing the same code twice. The validator drives the prot2exon binary and shells out to `Rscript` for ensembldb; reproduce it with [`validate_vs_ensembldb.py`](https://github.com/SotoLF/Prot2Exon/blob/main/benchmarks/validate_vs_ensembldb.py) and the [`validation.ipynb`](https://github.com/SotoLF/Prot2Exon/blob/main/notebooks/validation.ipynb) notebook.

The headline result is 100.00% exact match against ensembldb on a 5,000-query stratified set — zero off-by-ones and zero structural mismatches. Random sampling would underweight the corner cases that matter, so a 9-stratum sampler ([`sample_validation_queries.py`](https://github.com/SotoLF/Prot2Exon/blob/main/benchmarks/sample_validation_queries.py)) ensures every condition that historically breaks these tools is represented:

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

The cross-validation also extends to **`GenomicFeatures::proteinToGenome`**, a second independent Bioconductor implementation (its GRangesList method maps against a CDS-by-transcript `GRangesList` held in memory, rather than querying the EnsDb SQLite). On a 1,000-query v86 subset all three tools agree on every query — prot2exon ≡ ensembldb ≡ GenomicFeatures at **1,000 / 1,000 exact-segment match**, pairwise. Reproduce with [`proteintogenome_bench.R`](https://github.com/SotoLF/Prot2Exon/blob/main/benchmarks/proteintogenome_bench.R) and [`compare_intervals.py`](https://github.com/SotoLF/Prot2Exon/blob/main/benchmarks/compare_intervals.py).

A second run tests annotation drift, comparing prot2exon on a current GENCODE GTF against an Ensembl 113 EnsDb pulled via AnnotationHub. There it is 100% exact match on the 2,264 / 5,000 queries where both tools return data. The remaining 2,736 fall into `only_prot2exon` because those ENSPs exist in the v113 EnsDb's protein table but the EnsDb has no CDS linkage for them, so `proteinToGenome()` reports "No CDS found." That is a gap in the EnsDb rather than a prot2exon disagreement, confirmed by running prot2exon on the same GTF the ENSPs came from and getting complete answers.

A few practical points matter when reproducing this. Installing ensembldb is best done through BiocManager from CRAN/Bioconductor rather than bioconda, whose `r-rjson` recipe pins an ancient `r=3.3.1` that conflicts with modern r-base; the R compile also needs `gcc ≤ 14` (gcc 15+ rejects `rtracklayer`'s bundled UCSC C code) and the development `libxml2` headers (`libxml2=2.13.9`, not the runtime-only conda-forge 2.15). The sampler also handles both GTF dialects — GENCODE uses `transcript_type` and puts `protein_id` on transcript rows, while Ensembl uses `transcript_biotype` and puts `protein_id` only on CDS rows — and skipping that distinction silently halves the number of usable queries.

## Speed vs other tools

The same 5,000-query stratified set drives a head-to-head against the three tools users most often choose between: ensembldb (the R/Bioconductor canonical), TransVar (the HGVS-based variant-annotation perspective, popular with clinical teams), and Ensembl REST (the no-install zero-overhead path). Tools considered and rejected include GeneMANIA (no per-domain output), peptidomics tools (a different problem), and VEP (genome→protein, the opposite direction). Reproduce with [`scaling_benchmark.py`](https://github.com/SotoLF/Prot2Exon/blob/main/benchmarks/scaling_benchmark.py), [`run_transvar.py`](https://github.com/SotoLF/Prot2Exon/blob/main/benchmarks/run_transvar.py), [`run_ensembl_rest.py`](https://github.com/SotoLF/Prot2Exon/blob/main/benchmarks/run_ensembl_rest.py), and the [`software_comparison.ipynb`](https://github.com/SotoLF/Prot2Exon/blob/main/notebooks/software_comparison.ipynb) notebook.

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

That is roughly a 900× speedup over ensembldb at N = 10,000 with a smaller index, about 4.4× faster than TransVar with no FASTA required, and identical genomic intervals against every tool that returned data. The agreement denominators differ because each tool answers a slightly different question: ensembldb is 100% on a stratified set built from the same annotation it indexes (no drift, the strongest claim); TransVar is 100% on the 1,761 queries whose ENST is in its bundled annotation, with the other 3,239 falling into `only_prot2exon` from annotation drift; and Ensembl REST is 98.30% on the 1,000 queries the rate limit allows, where all 17 disagreements are off-by-one under the codon-split convention with zero structural mismatches, and the other tools agree with prot2exon rather than REST on those rows.

Raw single-thread scaling:

| N | prot2exon wall (median, s) | prot2exon RSS (MB) | ensembldb wall (s) | ensembldb RSS (MB) |
|---:|---:|---:|---:|---:|
| 100 | 1.33 | 659 | 23.3 | 979 |
| 1,000 | 1.37 | 671 | 168.8 | 988 |
| 10,000 | 1.71 | 788 | 1,558.2 | 1,252 |
| 100,000 | 4.66 | 1,959 | (skipped) | — |
| 1,000,000 | 129.4 | 11,045 | (skipped) | — |

ensembldb was capped at N = 10K — linear extrapolation puts its N = 100K at ~4.3 h and N = 1M at ~43 h, so continuing past 10K would have contributed nothing beyond time burned. The prot2exon side at small N is dominated by the ~1.3 s one-time index load, so the actual mapping work for the first 1,000 queries is essentially free and wall time at N = 100 ≈ N = 1,000; this is why per-query throughput is misleading at small N, since it counts index load against per-query work.

The parallel scaling matches the single-tool numbers in the [Speed and memory usage](#--threads-openmp) section above (1.65× at 2 threads, plateauing past 4 as page-cache flushing and disk bandwidth dominate). For RAM-bounded large-N runs, the same 1M benchmark shows the ~11× peak-RSS reduction from `--batch-size 10000` described there.

A few notes on the external comparators. Ensembl REST is network-bound rather than rate-limited — the 15 q/s cap is not the bottleneck, the ~900 ms per-request HTTP round-trip is, and a concurrent keep-alive client could push closer to the cap but at the cost of measuring something other than how anyone actually writes a REST script. TransVar reports only the genomic envelope (one `chrN:g.start_end` per query, introns included) rather than per-CDS intervals, so the comparison collapses both sides to `(chrom, min_start, max_end)` via [`classify_external.py`](https://github.com/SotoLF/Prot2Exon/blob/main/benchmarks/classify_external.py) `--envelope-only` for an apples-to-apples result, meaning prot2exon is doing the harder per-exon decomposition while TransVar answers an easier question. TransVar also keys on ENST and silently returns empty for ENSP input, so it is fed ENSTs from the query metadata, and `transvar config` is interactive on first run (pipe an empty string to satisfy the FASTA prompt).

### GenomicFeatures::proteinToGenome — the GRanges path

`proteinToGenome` exists in **two** Bioconductor packages: ensembldb (EnsDb/SQLite-backed) and GenomicFeatures (a GRangesList method that maps in memory against a CDS-by-transcript `GRangesList`). The GenomicFeatures path trades a one-time setup (`cdsBy(edb, "tx")` + a protein→transcript map) for SQLite-free mapping, so it is meaningfully faster than ensembldb — but it is still an R-level per-query loop, far short of the indexed C++ path. On 1,000 v86 queries, one machine, one thread:

| Tool | Setup | Map (1,000) | Total | Map throughput | Peak RSS |
|---|---:|---:|---:|---:|---:|
| **prot2exon** | 1.5 s (index load) | **0.02 s** | **1.6 s** | ~50,000 q/s | 734 MB |
| GenomicFeatures::proteinToGenome | 9.2 s (`cdsBy`) | 48.5 s | 57.7 s | ~21 q/s | 1,365 MB |
| ensembldb::proteinToGenome | — | 201.9 s | 201.9 s | ~5 q/s | 1,163 MB |

GenomicFeatures is ~3.5× faster than ensembldb end-to-end, confirming the GRanges representation is the lighter of the two — but prot2exon's mapping is still ~2,400× faster than GenomicFeatures' and ~10,000× faster than ensembldb's per query, at roughly half the RAM, and all three return identical coordinates (see Accuracy above). One subtlety when reproducing: loading an EnsDb pulls in ensembldb's own `proteinToGenome` GRangesList method (which requires protein-sequence metadata on the CDS), so the runner fetches GenomicFeatures' method explicitly via `getMethod(..., where = asNamespace("GenomicFeatures"))`.

### Visualization tools (geneplot, VisProDom)

Two other tools that "show domains on gene structure" are sometimes mentioned alongside these, but they solve the *visualization* problem, not the coordinate-mapping one, so they are not in the tables above:

- **geneplot** (Python; `gffutils` + BioPython `GenomeDiagram`) consumes coordinates that are *already* mapped — a GFF3 plus an InterProScan domain file — and renders them. It performs no de-novo protein→genome mapping, so it is a peer of prot2exon's [[Plotting]], not its mapper.
- **VisProDom** (R/Shiny) is a domain viewer, but its `CreDat(gff, annofile)` function *is* a real batch mapper (pure R/dplyr, cumulative-CDS arithmetic — the same idea prot2exon implements in C++). It has no prebuilt index, so it is **O(genome) per call**: it rebuilds every transcript's CDS layout on each invocation. On its bundled maize example (619K-row GFF, ~76K domain hits) it maps 262K domain-segments in **17.9 s at 780 MB**, and because the cost is dominated by the genome rebuild rather than the query count, wall time is essentially flat — ~6–8 s — whether you map 100 domains or 10,000. That batch-recompute model is exactly what prot2exon's persistent binary index removes. Reproduce with [`visprodom_bench.R`](https://github.com/SotoLF/Prot2Exon/blob/main/benchmarks/visprodom_bench.R).
