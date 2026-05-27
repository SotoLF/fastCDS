# Validation

## Headline

| Question | Answer |
|---|---|
| Are the coordinates correct? | **100.00 % exact match vs ensembldb** on 5,000 stratified queries — *zero* off-by-ones, *zero* structural mismatches. |
| Is it fast enough at proteome scale? | **5,847 queries/s** on a single thread — ~900× ensembldb, ~4.4× TransVar, ~5,400× Ensembl REST. |
| Same on a current annotation? | Yes — 100 % on the 2,264/5,000 GENCODE-v49 queries where ensembldb's v113 EnsDb has the protein indexed (the remaining 2,736 are an EnsDb gap, not a prot2exon disagreement). |

## How the validation set is built

Random sampling would underweight the corner cases that actually matter. We use a **9-stratum sampler** so every condition that historically breaks protein-to-genome tools is represented:

| Stratum | n | What it stresses |
|---|---:|---|
| `single_exon_domain` | 1,000 | The common shape (high test coverage) |
| `multi_exon_domain` | 1,000 | The hard case — find all CDS pieces |
| `codon_split_boundary` | 500 | Codons straddling exon boundaries — where off-by-ones hide |
| `plus_strand_gene` | 1,000 | Strand-handling A/B |
| `minus_strand_gene` | 1,000 | Strand-handling A/B (most bugs live here) |
| `cds_incomplete` | 200 | `cds_start_NF` / `cds_end_NF` transcripts |
| `selenoprotein` | 100 | 25-gene curated list (UGA → Sec recoding) |
| `single_exon_gene` | 100 | Boundary case for the multi-exon path |
| `many_exon_gene` | 100 | > 20 CDS exons — exon-walker stress |

Selenoprotein gene list is hardcoded in `benchmarks/sample_validation_queries.py` (canonical UniProt selenoprotein curation, 25 genes).

## Result detail (matched annotation, v86 path)

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

## Annotation-drift run (current-release path)

Same comparison logic, **GENCODE v49 GTF for prot2exon, Ensembl 113 EnsDb via AnnotationHub for ensembldb**:

- **100 % exact match** on the 2,264 / 5,000 queries where both tools return data.
- The other 2,736 fall into `only_prot2exon`: those ENSPs are present in the v113 EnsDb's `protein` table but the EnsDb has no CDS linkage for them, so `proteinToGenome()` reports "No CDS found". That's an EnsDb gap — not a prot2exon bug. (Verified by running prot2exon on the same v49 GTF the ENSPs came from and getting a complete answer.)

## Classifier buckets

| Bucket | Definition |
|---|---|
| `exact_match` | Same `(chrom, start, end)` set from both tools |
| `off_by_one` | Sets differ but total bp differs by ≤ 2 (codon-split convention) |
| `structural_mismatch` | Both returned intervals; sets disagree and aren't off-by-one |
| `only_prot2exon` | prot2exon returned ≥ 1 interval; ensembldb returned nothing |
| `only_ensembldb` | Converse |
| `neither_mapped` | Both empty (excluded from the percentage denominator) |

## Why ensembldb is the reference comparator

- **Bioconductor canonical** for protein-to-genome mapping.
- ~800 paper citations.
- **Independent implementation** (R/SQL on top of EnsDb), so an agreement is genuine cross-validation, not testing the same code twice.

## Reproducing

The validator drives the C++ binary and shells out to `Rscript` for ensembldb. Both must be on `$PATH`.

```bash
# 0) Create the conda env (yml ships under benchmarks/)
conda env create -f benchmarks/environment.yml
conda activate prot2exon-val
Rscript -e 'install.packages("BiocManager", repos="https://cran.r-project.org"); \
            BiocManager::install(c("ensembldb", "EnsDb.Hsapiens.v86", "AnnotationHub"), \
                                  ask=FALSE, update=FALSE)'

# 1) Annotation + index
wget http://ftp.ensembl.org/pub/release-86/gtf/homo_sapiens/Homo_sapiens.GRCh38.86.chr.gtf.gz
gunzip Homo_sapiens.GRCh38.86.chr.gtf.gz
./build/prot2exon --gtf Homo_sapiens.GRCh38.86.chr.gtf --build-index --index human_v86.idx

# 2) 5,000 stratified queries
python benchmarks/sample_validation_queries.py \
    --gtf Homo_sapiens.GRCh38.86.chr.gtf \
    --out-bed queries_v86.bed \
    --out-meta queries_v86_meta.tsv

# 3) Run the validation (emits table1.tsv + discrepancies.tsv)
EnsDb_v86_path=$(Rscript -e 'cat(system.file("extdata/EnsDb.Hsapiens.v86.sqlite",
                                              package="EnsDb.Hsapiens.v86"))')
python benchmarks/validate_vs_ensembldb.py \
    --queries-bed queries_v86.bed \
    --queries-meta queries_v86_meta.tsv \
    --prot2exon-index human_v86.idx \
    --ensdb "$EnsDb_v86_path" \
    --out-dir validation_v86
```

Outputs:

- `validation_v86/prot2exon/` — raw prot2exon output (4 TSVs + 3 BEDs + BED12)
- `validation_v86/ensembldb_intervals.tsv` — output of the R helper
- `validation_v86/table1.tsv` — the agreement table
- `validation_v86/discrepancies.tsv` — per-query diff for non-exact-match rows

See [[Benchmarks]] for the speed comparisons against ensembldb / TransVar / Ensembl REST.

## Gotchas worth knowing if you try to reproduce

1. **`conda install -c bioconda bioconductor-ensembldb` never solves.** Bioconda's `r-rjson` recipe declares `r=3.3.1` as a dependency — a metadata bug that conflicts with any modern r-base. Use BiocManager via CRAN/Bioconductor (compiles from source but works).
2. **R compile inside conda needs the env activated AND `gcc ≤ 14`.** gcc 15+ rejects `rtracklayer`'s bundled UCSC C code. Pin `gcc_linux-64=13`.
3. **conda-forge `libxml2 2.15` is runtime-only** (46 KB). The R `XML` package's configure script needs headers — use `libxml2=2.13.9` (681 KB).
4. **`ensembldb::ensDbFromGtf` does not import protein annotations.** Only transcripts. Use packaged EnsDbs (`EnsDb.Hsapiens.v86`) or `AnnotationHub` for newer releases.
5. **The AnnotationHub v113 EnsDb has many ENSPs without CDS linkage** — that's why the v113 run shows `only_prot2exon` for 2,736 / 5,000 queries.
6. **GTF dialect mismatch.** GENCODE uses `transcript_type` and puts `protein_id` on transcript rows; Ensembl uses `transcript_biotype` and puts `protein_id` only on CDS rows. The sampler handles both — skip this and you silently get half as many usable queries.
