# benchmarks/

Reproduction harness for the correctness + speed comparisons reported in [[Validation]] and [[Benchmarks]] on the wiki.

## At a glance

| Question | Headline |
|---|---|
| Are the coordinates correct? | **100.00 % exact match vs ensembldb** on 5,000 stratified queries — zero off-by-ones, zero structural mismatches. |
| Is it fast enough? | **5,847 q/s** on one thread — **~900× ensembldb, ~4.4× TransVar, ~5,400× Ensembl REST**. |

See [`wiki/Validation.md`](../wiki/Validation.md) and [`wiki/Benchmarks.md`](../wiki/Benchmarks.md) for the full result tables, design rationale (9-stratum sampler, why each comparator), and the gotchas worth knowing if you reproduce.

## Scripts

| File | Purpose |
|---|---|
| `sample_validation_queries.py` | GTF parser + 9-stratum sampler (1 K plus / 1 K minus / 1 K single-exon / 1 K multi-exon / 500 codon-split / 200 incomplete CDS / 100 selenoprotein / 100 single-exon-gene / 100 many-exon-gene) |
| `ensembldb_query.R` | Batch `proteinToGenome` via R subprocess |
| `validate_vs_ensembldb.py` | Runs prot2exon + ensembldb, classifies into 6 buckets, emits Table 1 |
| `scaling_benchmark.py` | prot2exon vs ensembldb at N = 100 … 1 M |
| `parallel_benchmark.py` | prot2exon at threads 1, 2, 4, 8 |
| `run_ensembl_rest.py` | Rate-limited REST client |
| `run_transvar.py` | Builds HGVS from EnsDb sequences, drives TransVar |
| `classify_external.py` | Bucket-classifier (use `--envelope-only` for TransVar) |
| `make_scaling_outputs.py` | 2-tool Table 1 + scaling.png |
| `make_table_1.py` | 4-tool Table 1 |
| `make_figure_1.py` | 4-panel composite (paper Figure 1) |

## Reproducing in one block

```bash
# 0) Environment (R + Bioconductor — conda is the path of least friction)
conda env create -f benchmarks/environment.yml
conda activate prot2exon-val
Rscript -e 'install.packages("BiocManager", repos="https://cran.r-project.org"); \
            BiocManager::install(c("ensembldb", "EnsDb.Hsapiens.v86", "AnnotationHub"), \
                                  ask=FALSE, update=FALSE)'

# 1) Build the prot2exon index
./build/prot2exon --gtf Homo_sapiens.GRCh38.86.chr.gtf \
    --build-index --index human_v86.idx

# 2) 5,000 stratified queries
python benchmarks/sample_validation_queries.py \
    --gtf Homo_sapiens.GRCh38.86.chr.gtf \
    --out-bed queries_v86.bed --out-meta queries_v86_meta.tsv

# 3) Correctness validation
EnsDb_v86_path=$(Rscript -e 'cat(system.file("extdata/EnsDb.Hsapiens.v86.sqlite",
                                              package="EnsDb.Hsapiens.v86"))')
python benchmarks/validate_vs_ensembldb.py \
    --queries-bed queries_v86.bed --queries-meta queries_v86_meta.tsv \
    --prot2exon-index human_v86.idx --ensdb "$EnsDb_v86_path" \
    --out-dir validation_v86

# 4) Scaling + parallel
python benchmarks/scaling_benchmark.py \
    --bin build/prot2exon --p2e-index human_v86.idx \
    --ensdb "$EnsDb_v86_path" \
    --rscript $CONDA_PREFIX/bin/Rscript --r-helper benchmarks/ensembldb_query.R \
    --source-bed queries_v86.bed --work-dir bench \
    --sizes 100 1000 10000 100000 1000000 \
    --p2e-reps 2 --ensembldb-reps 1 --ensembldb-max-n 10000 \
    --out bench/timings.tsv

python benchmarks/parallel_benchmark.py \
    --bin build/prot2exon --p2e-index human_v86.idx \
    --bed bench/queries_n100000.bed --work-dir bench/parallel \
    --out bench/parallel.tsv
```

Total wall: validation ~5 min, scaling ~50 min (ensembldb N = 10K is ~26 min of that), parallel ~1 min.

## Gotchas (each one cost real time)

1. **`conda install -c bioconda bioconductor-ensembldb` never solves.** Bioconda's `r-rjson` recipe declares `r=3.3.1` — metadata bug, conflicts with any modern r-base. Use BiocManager instead.
2. **R compile in conda needs `gcc ≤ 14`.** gcc 15+ rejects `rtracklayer`'s bundled UCSC C code.
3. **conda-forge `libxml2 2.15` is runtime-only.** R's `XML` package needs headers — use `libxml2=2.13.9`.
4. **`ensembldb::ensDbFromGtf` does not import protein annotations.** Use packaged EnsDbs (e.g. `EnsDb.Hsapiens.v86`) or `AnnotationHub`.
5. **AnnotationHub's v113 EnsDb has many ENSPs without CDS linkage.** That's why the v113 path shows 2,736 / 5,000 `only_prot2exon` — an EnsDb gap, not a mapping bug.
6. **GTF dialect mismatch.** GENCODE: `transcript_type` + `protein_id` on transcript rows. Ensembl: `transcript_biotype` + `protein_id` only on CDS rows. The sampler handles both — skip it and you silently get half the usable queries.
7. **TransVar throws away `protein_id` on ENSP input — keys on ENST.** Feeding ENSPs returns empty silently. We feed ENSTs via the `queries_meta.tsv` `transcript_id` column.
8. **`transvar config` is interactive on first run.** Pipe `echo "" | transvar config …` to satisfy the FASTA prompt.
