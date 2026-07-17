# matched/ - scaling benchmark, same machine and same queries for every tool

Backs the scaling curves in `notebooks/scaling_and_ram.ipynb` (**Figure 1B**,
**Figure S2**) and the speed / memory rows of **Table S2**.

Every tool here is run on the **same machine**, over the **same query sets**,
**single-threaded**, and timed **end-to-end**: wall clock from process start
until the last result is written, with peak RSS taken from `wait4()`. Nothing is
excluded from a tool's time, including its own index or database load, which is
the only way the ratios mean anything.

Each tool climbs the ladder N = 10² to 10⁶ and stops at the last N it can finish
in a practical time and memory budget:

| Tool | Route | Stops at | Why |
|---|---|---|---|
| fastCDS | `fastCDS map --threads 1` | 1,000,000 | - |
| ensembldb | `proteinToGenome(IRanges, EnsDb)` | 10,000 | SQLite hit per query, ~28 min at 10⁴ |
| GenomicFeatures | `proteinToGenome,GRangesList` on a `cdsBy()` list | 10,000 | in-memory, ~5x faster per query than ensembldb, still far off fastCDS |
| geneplot | per-gene transcript-to-genome walk | 10,000 | no index, ~0.12 s/query, 10⁵ would take hours |
| Ensembl REST | one `/map/translation` call per query | 1,000 | network-bound, ~900 ms per round-trip |

## Files

| File | Purpose |
|---|---|
| `prepare_matched_inputs.py` | Builds the shared query ladder (`q_<N>.bed`, `ids_<N>.txt`, `enst_<N>.txt`) from the GTF. Each step is a prefix of the next. |
| `run_matched_scaling.py` | The driver: runs every tool at every N under a `wait4()` wall + peak-RSS harness, appends to `results.tsv`, and assembles `scaling_matched.tsv` for the notebook. |
| `run_ensembldb.R` | ensembldb route (EnsDb method, SQLite-backed). |
| `run_gf_granges.R` | GenomicFeatures route (GRangesList method, in-memory). |
| `run_geneplot.py` | geneplot route. |
| `run_rest.py` | Ensembl REST route (timing only; the agreement runs use `../run_ensembl_rest.py`). |

## Reproducing

```bash
# 0) One-time: the GenomicFeatures route needs a CDS-by-transcript GRangesList.
Rscript -e 'library(ensembldb); library(EnsDb.Hsapiens.v86);
            saveRDS(cdsBy(EnsDb.Hsapiens.v86, by="tx"), "'"$FASTCDS_DATA"'/matched/cds_by_tx.rds")'

# 1) Shared query ladder
python prepare_matched_inputs.py \
    --gtf Homo_sapiens.GRCh38.86.chr.gtf \
    --out-dir "$FASTCDS_DATA/matched"

# 2) Full sweep + assemble
EnsDb_v86_path=$(Rscript -e 'cat(system.file("extdata/EnsDb.Hsapiens.v86.sqlite",
                                              package="EnsDb.Hsapiens.v86"))')
python run_matched_scaling.py \
    --work-dir "$FASTCDS_DATA/matched" \
    --bin ../../../build/fastCDS --index "$FASTCDS_DATA/human_v86.idx" \
    --ensdb "$EnsDb_v86_path" \
    --gff "$FASTCDS_DATA/human_tool_bench/h86.gff3" \
    --ipr "$FASTCDS_DATA/human_tool_bench/h86.ipr" \
    --ensp-enst "$FASTCDS_DATA/human_tool_bench/ensp_enst.tsv"

# re-measure one tool, then rebuild the table
python run_matched_scaling.py --work-dir "$FASTCDS_DATA/matched" --tools ensembldb --sizes 1000
python run_matched_scaling.py --work-dir "$FASTCDS_DATA/matched" --assemble-only
```

The full sweep is dominated by the slow tools: ensembldb and GenomicFeatures at
N = 10⁴ are ~28 min and ~6 min, geneplot at 10⁴ is ~21 min, and REST at 10³ is
~13 min of round-trips. fastCDS finishes the whole ladder, 10⁶ included, in
under 4 min.

`results.tsv` is append-only and keyed on `(tool, n)`; a re-measured pair
supersedes the earlier row at assemble time, so re-running a single tool is
safe.

The geneplot inputs (`--gff`, `--ipr`, `--ensp-enst`) come from
[`../build_human_tool_inputs.py`](../build_human_tool_inputs.py), which converts
the Ensembl-86 GFF3 plus Pfam-on-v86 domains into geneplot's InterProScan
format. geneplot ships only fruit-fly example data, but it is general (any GFF
plus domain file), so running it on the same human set is what makes it
comparable here.
