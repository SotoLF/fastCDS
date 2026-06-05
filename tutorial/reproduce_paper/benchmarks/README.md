# tutorial/reproduce_paper/benchmarks/

Reproduction harness for the correctness + speed comparisons reported in [[Validation]] and [[Benchmarks]] on the wiki.

## At a glance

| Question | Headline |
|---|---|
| Are the coordinates correct? | **100.00 % exact match vs ensembldb _and_ GenomicFeatures::proteinToGenome** — three-way agreement, zero off-by-ones, zero structural mismatches. |
| Is it fast enough? | **5,847 q/s** on one thread (N = 10,000, end-to-end) — **~970× ensembldb, ~4.4× TransVar, ~5,400× Ensembl REST**. |

### How "speedup" is measured

All speedup ratios below are **end-to-end at a fixed N**: total wall time from
process start until every result is written, **including the one-time index /
database load**, divided by N. We report N = 10,000 on one thread. This is the
number a user actually waits for, and it avoids cherry-picking — so quote *this*
ratio, not a warm-cache one.

Two consequences worth stating plainly so the numbers don't look contradictory:

- **The ratio grows with N.** prot2exon's ~1.2 s index load is amortized over
  more queries as N rises, while ensembldb's per-query cost dominates throughout.
  So prot2exon-vs-ensembldb is ~130× at N = 1,000 but ~970× at N = 10,000. Always
  state N with a speedup.
- **Warm mapping-only is faster still but biased**, so we don't headline it: with
  the index already loaded, prot2exon maps at ~71,000 q/s, which is ~11,000×
  ensembldb's per-query rate — but that figure drops prot2exon's load cost while
  keeping the comparators' steady-state, so it flatters us. Use it only to explain
  *where* the speed comes from (the mapping is near-free; I/O is the floor).

See [`wiki/Performance-and-Benchmarking.md`](../wiki/Performance-and-Benchmarking.md) for the full result tables, design rationale (9-stratum sampler, why each comparator), the GenomicFeatures/ensembldb/VisProDom comparison, and the practical notes worth knowing if you reproduce.

## Scripts

| File | Purpose |
|---|---|
| `sample_validation_queries.py` | GTF parser + 9-stratum sampler (1 K plus / 1 K minus / 1 K single-exon / 1 K multi-exon / 500 codon-split / 200 incomplete CDS / 100 selenoprotein / 100 single-exon-gene / 100 many-exon-gene) |
| `ensembldb_query.R` | Batch `proteinToGenome` via R subprocess |
| `proteintogenome_bench.R` | Speed + RAM for `ensembldb::proteinToGenome` (SQLite-backed) vs `GenomicFeatures::proteinToGenome` (GRangesList, in-memory) on the same queries |
| `visprodom_bench.R` | Speed + RAM for VisProDom's `CreDat()` batch mapper (pure R/dplyr; O(genome) per call) |
| `compare_intervals.py` | Per-query exact-segment agreement across mappers (ensembldb / GenomicFeatures / prot2exon) |
| `validate_vs_ensembldb.py` | Runs prot2exon + ensembldb, classifies into 6 buckets, emits Table 1 |
| `scaling_benchmark.py` | prot2exon vs ensembldb at N = 100 … 1 M |
| `threads_batch_grid.py` | combined `--threads` × `--batch-size` grid (wall + peak RSS in one sweep) |
| `parallel_benchmark.py` | threads-only sweep (superseded by `threads_batch_grid.py`; kept for back-compat) |
| `run_ensembl_rest.py` | Rate-limited REST client |
| `run_transvar.py` | Builds HGVS from EnsDb sequences, drives TransVar |
| `classify_external.py` | Bucket-classifier (use `--envelope-only` for TransVar) |
| `make_scaling_outputs.py` | 2-tool Table 1 + scaling.png |
| `make_table_1.py` | 4-tool Table 1 |
| `make_figure_1.py` | 4-panel composite (paper Figure 1) |

## Reproducing in one block

```bash
# 0) Environment (R + Bioconductor — conda is the path of least friction)
conda env create -f tutorial/reproduce_paper/benchmarks/environment.yml
conda activate prot2exon-val
Rscript -e 'install.packages("BiocManager", repos="https://cran.r-project.org"); \
            BiocManager::install(c("ensembldb", "GenomicFeatures", \
                                    "EnsDb.Hsapiens.v86", "AnnotationHub"), \
                                  ask=FALSE, update=FALSE)'

# 1) Build the prot2exon index
./build/prot2exon index --gtf Homo_sapiens.GRCh38.86.chr.gtf --out human_v86.idx

# 2) 5,000 stratified queries
python tutorial/reproduce_paper/benchmarks/sample_validation_queries.py \
    --gtf Homo_sapiens.GRCh38.86.chr.gtf \
    --out-bed queries_v86.bed --out-meta queries_v86_meta.tsv

# 3) Correctness validation
EnsDb_v86_path=$(Rscript -e 'cat(system.file("extdata/EnsDb.Hsapiens.v86.sqlite",
                                              package="EnsDb.Hsapiens.v86"))')
python tutorial/reproduce_paper/benchmarks/validate_vs_ensembldb.py \
    --queries-bed queries_v86.bed --queries-meta queries_v86_meta.tsv \
    --prot2exon-index human_v86.idx --ensdb "$EnsDb_v86_path" \
    --out-dir validation_v86

# 4) Scaling + parallel
python tutorial/reproduce_paper/benchmarks/scaling_benchmark.py \
    --bin build/prot2exon --p2e-index human_v86.idx \
    --ensdb "$EnsDb_v86_path" \
    --rscript $CONDA_PREFIX/bin/Rscript --r-helper tutorial/reproduce_paper/benchmarks/ensembldb_query.R \
    --source-bed queries_v86.bed --work-dir bench \
    --sizes 100 1000 10000 100000 1000000 \
    --p2e-reps 2 --ensembldb-reps 1 --ensembldb-max-n 10000 \
    --out bench/timings.tsv

# combined threads × batch-size grid (wall + peak RSS in one sweep).
# one-shot at N=1M holds ~13 GB in RAM — drop the `0` batch on a small-RAM box.
python tutorial/reproduce_paper/benchmarks/threads_batch_grid.py \
    --bin build/prot2exon --index human_v86.idx \
    --bed bench/queries_n1000000.bed --work-dir bench/grid \
    --threads 1 4 8 16 32 --batch-sizes 0 10000 50000 100000 --reps 2 \
    --out bench/threads_batch_grid.tsv

# 5) ensembldb vs GenomicFeatures::proteinToGenome (speed + RAM, same queries)
head -1000 queries_v86.bed > q1k.bed
for tool in ensembldb genomicfeatures; do
  Rscript tutorial/reproduce_paper/benchmarks/proteintogenome_bench.R $tool "$EnsDb_v86_path" \
      q1k.bed ${tool}_intervals.tsv ${tool}_timing.tsv
done
build/prot2exon map --index human_v86.idx --bed q1k.bed --out-dir p2e_q1k --output coding
python tutorial/reproduce_paper/benchmarks/compare_intervals.py \
    ensembldb=ensembldb_intervals.tsv \
    genomicfeatures=genomicfeatures_intervals.tsv \
    prot2exon=p2e_q1k/domain_cds_segments.tsv

# 6) VisProDom CreDat() batch-mapper characterization (its own maize example data)
git clone --depth 1 https://github.com/whweve/VisProDom /tmp/VisProDom
Rscript tutorial/reproduce_paper/benchmarks/visprodom_bench.R /tmp/VisProDom
```

Total wall: validation ~5 min, scaling ~50 min (ensembldb N = 10K is ~26 min of that), parallel ~1 min, proteinToGenome head-to-head ~5 min (ensembldb N = 1K is ~3.5 min of that), VisProDom ~1 min.

Headline results (one machine, one thread): on 1,000 v86 queries prot2exon ≡ ensembldb ≡ GenomicFeatures at 1,000 / 1,000 exact-segment match; mapping time prot2exon 0.014 s vs GenomicFeatures 37.71 s vs ensembldb 160.12 s; see [`proteintogenome_results.tsv`](proteintogenome_results.tsv).

## Other tools on human data (VisProDom, geneplot)

VisProDom and geneplot ship only non-human example data (maize, fruit-fly), but
both are general (any GFF + domain file), so we run them on the **same human
Ensembl-86 set** for an apples-to-apples comparison. [`build_human_tool_inputs.py`](build_human_tool_inputs.py)
turns the Ensembl-86 GFF3 + the Pfam-on-v86 domains into each tool's format
(geneplot: InterProScan `.ipr`; VisProDom: a Phytozome-style GFF + RPS-BLAST
`annofile`), then [`geneplot_human.py`](geneplot_human.py) and
[`visprodom_human.R`](visprodom_human.R) run them.

Neither has an index, and on the human genome it shows: VisProDom rebuilds the
whole genome on every `CreDat` call (**42.7 s, 2.3 GB** for the 1,000-query set,
~19 proteins/s), and geneplot builds a `gffutils` SQLite database of the genome
first (minutes) then re-reads the domain file per gene (**~14 genes/s**
end-to-end). On their bundled examples both finish in seconds; on human they
land in the slow-tool range (vs prot2exon's ~830 q/s end-to-end on the same
1,000-query set, index load included).

## Reproduction notes

1. **`conda install -c bioconda bioconductor-ensembldb` never solves.** Bioconda's `r-rjson` recipe declares `r=3.3.1` — metadata bug, conflicts with any modern r-base. Use BiocManager instead.
2. **R compile in conda needs `gcc ≤ 14`.** gcc 15+ rejects `rtracklayer`'s bundled UCSC C code.
3. **conda-forge `libxml2 2.15` is runtime-only.** R's `XML` package needs headers — use `libxml2=2.13.9`.
4. **`ensembldb::ensDbFromGtf` does not import protein annotations.** Use packaged EnsDbs (e.g. `EnsDb.Hsapiens.v86`) or `AnnotationHub`.
5. **AnnotationHub's v113 EnsDb has many ENSPs without CDS linkage.** That's why the v113 path shows 2,736 / 5,000 `only_prot2exon` — an EnsDb gap, not a mapping bug.
6. **GTF dialect mismatch.** GENCODE: `transcript_type` + `protein_id` on transcript rows. Ensembl: `transcript_biotype` + `protein_id` only on CDS rows. The sampler handles both — skip it and you silently get half the usable queries.
7. **TransVar throws away `protein_id` on ENSP input — keys on ENST.** Feeding ENSPs returns empty silently. We feed ENSTs via the `queries_meta.tsv` `transcript_id` column.
8. **`transvar config` is interactive on first run.** Pipe `echo "" | transvar config …` to satisfy the FASTA prompt.
