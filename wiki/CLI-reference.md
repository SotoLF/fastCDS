# CLI reference

```
USAGE
  prot2exon --gtf FILE --build-index --index FILE
  prot2exon (--gtf FILE | --index FILE) --bed FILE --out-dir DIR [--output KIND]
  prot2exon plot --isoform FILE (--input-id ID | --all) [--out F] [--html F] [--html-interactive F]
  prot2exon fetch <subcommand> [options]
```

## Mapping

| Option | Effect |
|---|---|
| `--bed FILE` | Query file (TSV / BED-like). Rows can be `id  aa_start  aa_end  [domain_id]`, or just `id` for whole-transcript-structure mode. `id` can be an ENSP or ENST. |
| `--out-dir DIR` | Output directory (created if missing). |
| `--gtf FILE` | GTF annotation, parsed on the fly. |
| `--index FILE` | Pre-built binary index (faster, recommended). |
| `--output KIND` | One of `coding`, `introns`, `span`, `isoform`, `bed12`, `all`. Default: `all`. |
| `--build-index` | Build a binary index from `--gtf` into `--index`. |
| `--threads NUM` | Process queries in parallel via OpenMP. Default: 1. |
| `--batch-size N` | Stream results to disk in chunks of N to cap peak RAM (default 0 = unbounded, fastest when working set fits). See [Performance and RAM](Performance-and-RAM). |
| `--verbose` | Log progress to stderr. |
| `--version` | Print version and exit. |
| `--help` | Full help including all output schemas. |

Run `prot2exon --help` for the full schema reference at the terminal.

## Plotting

`prot2exon plot` renders `isoform_structure.tsv` to a figure. Every flag is documented on [Plotting and viewers](Plotting-and-viewers).

```bash
# One domain to PDF (matplotlib)
prot2exon plot --isoform results/isoform_structure.tsv \
    --input-id TP53_DBD --out tp53_dbd.pdf

# All input_ids in the TSV to a multipage PDF
prot2exon plot --isoform results/isoform_structure.tsv \
    --all --out queries.pdf

# Interactive plotly HTML
prot2exon plot --isoform results/isoform_structure.tsv \
    --input-id TP53_DBD --html tp53_dbd.html

# Self-contained vanilla JS HTML (no CDN)
prot2exon plot --isoform results/isoform_structure.tsv \
    --input-id TP53_DBD --html-interactive tp53_dbd.html
```

## Fetching indices

### Pre-built indexes (instant, from Zenodo)

```bash
# Discover what's available — both pre-built indexes AND GTF-build presets.
prot2exon fetch list

# Download a pre-built index by name — URL + sha256 are looked up internally.
prot2exon fetch index --preset human-v49

# Or pass a custom URL with --url (any Zenodo / HTTPS / file:// source).
prot2exon fetch index \
    --url https://zenodo.org/record/<RECORD>/files/gencode_v49_human.idx \
    --out ~/.cache/prot2exon/human.idx \
    --sha256 <hex-from-MANIFEST.tsv>
```

Available presets (also from `prot2exon fetch list`):

- `human-v49`  · GENCODE v49 (~298 MB, current human)
- `mouse-vM34` · GENCODE vM34 (~73 MB)
- `human-v86`  · Ensembl 86 (~87 MB, for validation reproducibility)
- `yeast`      · NCBI RefSeq S. cerevisiae R64 (~1.4 MB)

With `--preset`, `--sha256` is baked in and `--out` defaults to `~/.cache/prot2exon/<preset>.idx`. With `--url`, both are user-supplied (`--sha256` is optional but strongly recommended).

### Build from a GTF

```bash
# One-command download + build for stock species
prot2exon fetch human --release 49      # GENCODE human v49
prot2exon fetch mouse --release M34     # GENCODE mouse vM34
prot2exon fetch refseq --preset yeast   # NCBI RefSeq S. cerevisiae R64
prot2exon fetch refseq --preset ecoli   # NCBI RefSeq E. coli K-12 MG1655
prot2exon fetch ensembl --species danio_rerio --assembly GRCz11 --release 115

# List built-in presets
prot2exon fetch list

# Use a custom URL (any host)
prot2exon fetch human --release 49 \
    --gtf-url https://your.host/path/custom.gtf.gz
```

Indices land in `~/.cache/prot2exon/` by default (override with `--cache-dir` or `--out`). The downloaded `.gtf.gz` is cached so a rebuild is cheap; the uncompressed `.gtf` is removed after indexing unless you pass `--keep-gtf`. The command prints the final index path on stdout so it pipes cleanly into the mapper:

```bash
IDX=$(prot2exon fetch refseq --preset yeast)
prot2exon --index "$IDX" --bed queries.bed --out-dir results --output all
```

Re-running with the index already cached prints `using cached index: …  (--force to rebuild)`. Pass `--force` to download + rebuild from upstream.
