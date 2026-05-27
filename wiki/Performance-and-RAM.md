# Performance and RAM

> **Scope.** This page covers prot2exon's own performance knobs — `--threads`, `--batch-size`, peak RSS, scaling with N. For the head-to-head against other tools (ensembldb / TransVar / Ensembl REST), see [[Benchmarks]].

Measured on a single core against the pre-built human GENCODE v49 index (3.17 GB GTF → 128 MB binary index). Reproduce with `python3 benchmarks/run_benchmark.py && python3 benchmarks/plot_benchmark.py`.

## Index build (one-time per annotation)

| Metric | Value |
|---|---|
| GTF size (uncompressed) | 3,173 MB |
| Lines parsed | ~5.8 M |
| Wall time | **18.1 s** |
| Peak RSS | 2.1 GB |
| Output index size | 128 MB |
| Compression ratio | ~25× |

The index encodes per-protein exon and CDS intervals plus `gene_id`, `gene_name`, `exon_number`. Build it once and reuse it until you upgrade the annotation.

## Query performance (single-threaded, median of 3 reps)

| Queries | `--output all` (wall / peak RSS) | `--output isoform` (wall / peak RSS) |
|---|---|---|
| 100 | 1.7 s / 847 MB | 1.6 s / 847 MB |
| 1,000 | 2.1 s / 858 MB | 1.7 s / 859 MB |
| 10,000 | 5.5 s / 974 MB | 3.3 s / 974 MB |
| 100,000 | 37 s / 2,126 MB | 17.6 s / 2,125 MB |

Effective throughput at 100k: ~2,800 q/s in `all`, ~5,700 q/s in `isoform`. The `all` mode is slower because it writes four TSVs and three BEDs vs. one TSV in `isoform`.

The constant ~1.5 s floor at small N is index loading. Per-query cost is roughly linear and dominated by output formatting. Memory grows because output strings accumulate before being written.

## `--threads` (OpenMP)

Per-query processing runs in an OpenMP parallel loop; the output writes for the seven files also run concurrently via OpenMP `sections`.

```bash
./build/prot2exon --index human.idx --bed q.bed --out-dir out --output all --threads 8
```

Peak RSS is largely unaffected by `--threads`.

### Parallel scaling at N = 1,000,000

| Threads | Wall (median, s) | Speedup | Efficiency |
|---:|---:|---:|---:|
| 1 | 107.5 | 1.00 | 1.00 |
| 2 | 65.0 | 1.65 | 0.83 |
| 4 | 71.4 | 1.51 | 0.38 |
| 8 | 73.9 | 1.46 | 0.18 |

Mapping scales well — 0.83 efficiency at 2 threads. Speedup peaks at 2 threads and plateaus past 4 as **the OS page cache + disk bandwidth** take over from per-row formatting as the bottleneck. Mapping is OpenMP-parallel and the seven output files are written concurrently via `parallel sections`; what's left is the raw I/O cost.

The writer uses **1 MiB stdio buffers per file** (`pubsetbuf`, ~100× fewer syscalls than the default 8 KB) which cut single-threaded wall time by ~17 % over the legacy buffer size. Pushing past the 4-thread plateau would need a different layer entirely — direct `io_uring`, on-the-fly compression, or smaller TSV output — none of which are clearly worth the complexity until users complain.

## `--batch-size` (RAM cap)

By default the binary holds every query's full result set in memory before writing. Fast when it fits, but at ~10 KB per query a 1 M-query run needs ~10 GB.

Pass `--batch-size N` (Python: `Mapper(..., batch_size=N)`) to stream results to disk in chunks of N and free each chunk before processing the next. Peak RAM becomes `O(N × per-query result size)`.

### 1 M-query benchmark

Human GENCODE v86 index, `--output all --threads 4`, NVMe SSD, 16 GB RAM machine:

| Mode | Wall | Peak RSS | Notes |
|---|---:|---:|---|
| `--batch-size 0` (one-shot) | 122.5 s | 10.4 GB | swap-bound on this box; would be faster with 32 GB+ |
| `--batch-size 10000` | **58.5 s** | **961 MB** | byte-identical outputs to one-shot |

At 1 M queries on a 16 GB machine, batching cuts peak RAM ~11× and wall time ~2×. On a machine with enough RAM to avoid swap, one-shot is slightly faster — choose `--batch-size 0` for max throughput when RAM allows, `--batch-size N` to cap the working set.

`--batch-size 10000` is a sensible default for large runs. Outputs are byte-identical to the one-shot path; batching is purely a memory optimisation.

## Tuning checklist

- Use a Release build (`cmake -DCMAKE_BUILD_TYPE=Release`). Debug builds are ~5× slower.
- Persist the index — never re-parse the GTF per run.
- Set `--threads` to your physical core count, not hyperthread count.
- For 100k+ queries, prefer `--batch-size 10000` over the default `0` unless you've benchmarked and have plenty of headroom.
- On HDDs, the BED+TSV write phase becomes I/O-bound; `--threads` still helps because the writes are parallelised.

## Reproducing the benchmark

```bash
python3 benchmarks/run_benchmark.py \
    --sizes 100 1000 10000 100000 \
    --modes all isoform \
    --reps 3 --threads 1

python3 benchmarks/plot_benchmark.py
# results.tsv, index_build.json, scaling.png in benchmarks/
```

Add `--sizes 1000000` for a million-query stress test (allow ~5–10 min).
