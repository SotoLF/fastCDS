# reproduce_paper/benchmarks/

Reproduction harness for the correctness + speed comparisons reported in [[Performance-and-Benchmarking]] on the wiki.

## Results

| Question | Result |
|---|---|
| Are the coordinates correct? | 100.00% exact match against both Bioconductor `proteinToGenome` methods (ensembldb and GenomicFeatures) on the 5,000-query set: zero off-by-ones, zero structural mismatches. |
| Is it fast enough? | 5,886 q/s on one thread (N = 10,000, end-to-end): ~980x ensembldb, >200x GenomicFeatures, ~4,560x Ensembl REST (REST measured at N = 1,000, the largest it finishes). |

### How "speedup" is measured

Speedup is measured end-to-end at a fixed N: total wall time from process start
until every result is written, including the one-time index or database load,
divided by N. The reported figure is **N = 10,000 on one thread**, which charges
fastCDS for its own index load.

State N alongside the ratio, because it grows with N: fastCDS's ~1.5 s index load
is amortized over more queries as N rises, while ensembldb's per-query cost stays
constant, so fastCDS-vs-ensembldb is ~130x at N = 1,000 and ~980x at N = 10,000.
N = 10,000 is used as the headline figure throughout.

See [`wiki/Performance-and-Benchmarking.md`](../../wiki/Performance-and-Benchmarking.md) for the full result tables, the design rationale (9-category sampler, why each comparator), the GenomicFeatures / ensembldb comparison, and notes for reproducing.

## Scripts

| File | Purpose |
|---|---|
| `sample_validation_queries.py` | GTF parser + 9-category sampler (1 K plus / 1 K minus / 1 K single-exon / 1 K multi-exon / 500 codon-split / 200 incomplete CDS / 100 selenoprotein / 100 single-exon-gene / 100 many-exon-gene). `cds_incomplete` is **exclusive**: the other eight categories draw only from complete-CDS transcripts, so every `cds_start_NF`/`cds_end_NF` query lands in `cds_incomplete`. |
| `ensembldb_query.R` | Batch `ensembldb::proteinToGenome` (EnsDb method) via R subprocess |
| `genomicfeatures_query.R` | Batch `GenomicFeatures::proteinToGenome` (GRangesList method, in-memory) via R subprocess - the independent second reference implementation |
| `proteintogenome_bench.R` | Speed + RAM for `ensembldb::proteinToGenome` (SQLite-backed) vs `GenomicFeatures::proteinToGenome` (GRangesList, in-memory) on the same queries |
| `compare_intervals.py` | Per-query exact-segment agreement across mappers (ensembldb / GenomicFeatures / fastCDS) |
| `validate_vs_ensembldb.py` | Runs fastCDS + ensembldb, classifies into the agreement buckets, emits Supplementary Table S1 |
| `scaling_benchmark.py` | fastCDS vs ensembldb at N = 100 to 1 M |
| `threads_batch_grid.py` | combined `--threads` x `--batch-size` grid (wall + peak RSS in one sweep) |
| `parallel_benchmark.py` | threads-only sweep (superseded by `threads_batch_grid.py`; kept for back-compat) |
| `run_ensembl_rest.py` | Rate-limited REST client |
| `run_transvar.py` | Builds HGVS from EnsDb sequences, drives TransVar |
| `classify_external.py` | Bucket-classifier (use `--envelope-only` for TransVar) |
| `make_scaling_outputs.py` | scaling.png + wall-time summary |
| `make_table_s1.py` | **Table S1** - per-tool, per-category agreement (ensembldb, GenomicFeatures, TransVar, REST across the 9 categories; TransVar reported NA where a category has multi-CDS-block queries; REST split into exact / off-by-one / no-mapping) |
| `make_table_s2.py` | **Table S2** - mapping speed and peak memory across the five tools (fastCDS, GenomicFeatures, geneplot, ensembldb, REST) |
| `build_human_tool_inputs.py` | Builds geneplot's human inputs (`h86.ipr`, `ensp_enst.tsv`) from the Ensembl-86 GTF + the Pfam-on-v86 table |
| [`matched/`](matched/) | The scaling harness: every tool, same machine, same queries, one thread, N = 10² to 10⁶. See [`matched/README.md`](matched/README.md). |

## Reproducing in one block

```bash
# 0) Environment (R + Bioconductor - conda is the path of least friction)
conda env create -f reproduce_paper/benchmarks/environment.yml
conda activate fastCDS
Rscript -e 'install.packages("BiocManager", repos="https://cran.r-project.org"); \
            BiocManager::install(c("ensembldb", "GenomicFeatures", \
                                    "EnsDb.Hsapiens.v86", "AnnotationHub"), \
                                  ask=FALSE, update=FALSE)'

# 1) Build the fastCDS index
./build/fastCDS index --gtf Homo_sapiens.GRCh38.86.chr.gtf --out human_v86.idx

# 2) 5,000 stratified queries
python reproduce_paper/benchmarks/sample_validation_queries.py \
    --gtf Homo_sapiens.GRCh38.86.chr.gtf \
    --out-bed queries_v86.bed --out-meta queries_v86_meta.tsv

# 3) Correctness validation
EnsDb_v86_path=$(Rscript -e 'cat(system.file("extdata/EnsDb.Hsapiens.v86.sqlite",
                                              package="EnsDb.Hsapiens.v86"))')
python reproduce_paper/benchmarks/validate_vs_ensembldb.py \
    --queries-bed queries_v86.bed --queries-meta queries_v86_meta.tsv \
    --fastCDS-index human_v86.idx --ensdb "$EnsDb_v86_path" \
    --out-dir validation_v86

# 4) Scaling + parallel
python reproduce_paper/benchmarks/scaling_benchmark.py \
    --bin build/fastCDS --fastcds-index human_v86.idx \
    --ensdb "$EnsDb_v86_path" \
    --rscript $CONDA_PREFIX/bin/Rscript --r-helper reproduce_paper/benchmarks/ensembldb_query.R \
    --source-bed queries_v86.bed --work-dir bench \
    --sizes 100 1000 10000 100000 1000000 \
    --fastcds-reps 2 --ensembldb-reps 1 --ensembldb-max-n 10000 \
    --out bench/timings.tsv

# combined threads x batch-size grid (wall + peak RSS in one sweep).
# one-shot at N=1M holds ~13 GB in RAM - drop the `0` batch on a small-RAM box.
python reproduce_paper/benchmarks/threads_batch_grid.py \
    --bin build/fastCDS --index human_v86.idx \
    --bed bench/queries_n1000000.bed --work-dir bench/grid \
    --threads 1 4 8 16 32 --batch-sizes 0 10000 50000 100000 --reps 2 \
    --out bench/threads_batch_grid.tsv

# 5) ensembldb vs GenomicFeatures::proteinToGenome (speed + RAM, same queries)
head -1000 queries_v86.bed > q1k.bed
for tool in ensembldb genomicfeatures; do
  Rscript reproduce_paper/benchmarks/proteintogenome_bench.R $tool "$EnsDb_v86_path" \
      q1k.bed ${tool}_intervals.tsv ${tool}_timing.tsv
done
build/fastCDS map --index human_v86.idx --bed q1k.bed --out-dir fastcds_q1k --output coding
python reproduce_paper/benchmarks/compare_intervals.py \
    ensembldb=ensembldb_intervals.tsv \
    genomicfeatures=genomicfeatures_intervals.tsv \
    fastCDS=fastcds_q1k/domain_cds_segments.tsv
```

Total wall: validation ~5 min, scaling ~50 min (ensembldb N = 10K is ~26 min of that), parallel ~1 min, proteinToGenome head-to-head ~5 min (ensembldb N = 1K is ~3.5 min of that).

Headline results (one machine, one thread): on 1,000 v86 queries fastCDS = ensembldb = GenomicFeatures at 1,000 / 1,000 exact-segment match; mapping time fastCDS 0.014 s vs GenomicFeatures 37.71 s vs ensembldb 160.12 s; see [`proteintogenome_results.tsv`](proteintogenome_results.tsv).

## Matched-pairwise agreement

Table S1 compares each external tool against fastCDS on the **same Ensembl
release** (fastCDS re-indexed to match) and the same 9-stratum sample, so a
disagreement means a real coordinate difference and not an annotation-version
gap. Per tool: build the release's index, sample the strata with
[`sample_validation_queries.py`](sample_validation_queries.py), map with both
sides, then classify.

```bash
# TransVar pins to Ensembl 95, so that release is its match. It also keys on
# ENST and silently returns empty for ENSP, which is why run_transvar.py needs
# --queries-meta (it reads the transcript_id column from there).
./build/fastCDS index --gtf Homo_sapiens.GRCh38.95.gtf --out human_v95.idx
python sample_validation_queries.py --gtf Homo_sapiens.GRCh38.95.gtf \
    --out-bed strat_v95.bed --out-meta strat_v95_meta.tsv
./build/fastCDS map --index human_v95.idx --bed strat_v95.bed \
    --out-dir strat_v95_fastcds --output coding
python run_transvar.py strat_v95.bed tv_v95_raw.tsv \
    --ensdb "$EnsDb_v95_path" --queries-meta strat_v95_meta.tsv
# TransVar reports one genomic envelope per query, not per-CDS intervals, so
# both sides collapse to (chrom, min start, max end) before comparing.
python classify_external.py --external tv_v95_raw.tsv \
    --fastcds-cds strat_v95_fastcds/domain_cds_segments.tsv \
    --queries-meta strat_v95_meta.tsv --tool-name transvar \
    --envelope-only --out transvar_v95.tsv

# Ensembl REST tracks the current release. --workers keeps the 5,000-query run
# to minutes while --qps holds the aggregate under the public ceiling.
python run_ensembl_rest.py strat_v115.bed rest_v115_raw.tsv --limit 5000 --workers 12
python classify_external.py --external rest_v115_raw.tsv \
    --fastcds-cds strat_v115_fastcds/domain_cds_segments.tsv \
    --queries-meta strat_v115_meta.tsv --tool-name rest --out rest_v115.tsv

python make_table_s1.py --help   # for the full per-tool arg list
```

`software_comparison.ipynb` reads the resulting per-tool tables from
`$FASTCDS_DATA/validation_matched/`.

## geneplot on human data

geneplot ships only non-human example data (fruit-fly), but it is general (any
GFF + domain file), so we run it on the **same human Ensembl-86 set** every
other tool sees. [`build_human_tool_inputs.py`](build_human_tool_inputs.py)
turns the Ensembl-86 GTF + the Pfam-on-v86 domains into geneplot's InterProScan
`.ipr` format, then [`matched/run_geneplot.py`](matched/run_geneplot.py) runs it
under the same timing harness as everything else.

geneplot has no index, and on the human genome it shows: it builds a `gffutils`
SQLite database of the genome first (minutes) then re-reads the domain file per
gene (**~14 genes/s** end-to-end). On its bundled example it finishes in seconds;
on human it lands in the slow-tool range (vs fastCDS's ~830 q/s end-to-end on the
same 1,000-query set, index load included).

## Reproduction notes

1. **`conda install -c bioconda bioconductor-ensembldb` never solves.** Bioconda's `r-rjson` recipe declares `r=3.3.1` - metadata bug, conflicts with any modern r-base. Use BiocManager instead.
2. **R compile in conda needs `gcc ≤ 14`.** gcc 15+ rejects `rtracklayer`'s bundled UCSC C code.
3. **conda-forge `libxml2 2.15` is runtime-only.** R's `XML` package needs headers - use `libxml2=2.13.9`.
4. **`ensembldb::ensDbFromGtf` does not import protein annotations.** Use packaged EnsDbs (e.g. `EnsDb.Hsapiens.v86`) or `AnnotationHub`.
5. **AnnotationHub's v113 EnsDb has many ENSPs without CDS linkage.** That's why the v113 path shows 2,736 / 5,000 `only_fastCDS` - an EnsDb gap, not a mapping bug.
6. **GTF dialect mismatch.** GENCODE: `transcript_type` + `protein_id` on transcript rows. Ensembl: `transcript_biotype` + `protein_id` only on CDS rows. The sampler handles both - skip it and you silently get half the usable queries.
7. **TransVar throws away `protein_id` on ENSP input - keys on ENST.** Feeding ENSPs returns empty silently. We feed ENSTs via the `queries_meta.tsv` `transcript_id` column.
8. **`transvar config` is interactive on first run.** Pipe `echo "" | transvar config ...` to satisfy the FASTA prompt.
