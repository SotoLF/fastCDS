# Building an index

Every `prot2exon map` run needs a binary index of a genome's GTF annotation. You either build one from a GTF — obtain the GTF, then index it with `prot2exon index` — or skip all of that and retrieve a pre-built index from Zenodo with `prot2exon fetch`. Once you have an index, see [[Mapping]]. For installation, see [[Installation]].

## 1. Obtain a GTF

Download a GTF for your species and annotation source. prot2exon reads GENCODE, Ensembl, and NCBI RefSeq GTFs interchangeably.

- **GENCODE** (human and mouse): <https://www.gencodegenes.org/human/> — e.g. `gencode.v49.primary_assembly.annotation.gtf.gz`.
- **Ensembl** (any species): <https://www.ensembl.org/info/data/ftp/index.html> — e.g. `Homo_sapiens.GRCh38.110.gtf.gz`.
- **NCBI RefSeq**: via the NCBI genomes FTP — e.g. `GCF_000146045.2_R64_genomic.gtf.gz` (yeast).

```bash
# Example: GENCODE human v49
curl -O https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.primary_assembly.annotation.gtf.gz
gunzip gencode.v49.primary_assembly.annotation.gtf.gz
```

The three dialects differ only slightly, and prot2exon handles each: GENCODE and Ensembl share the format (`gene_name`, `protein_id`, `transcript_id`, plus `tag` for MANE Select / Ensembl_canonical); RefSeq uses `gene` instead of `gene_name` and carries no MANE tags (those columns report `NA`). IDs are matched with the `.version` suffix stripped on both the GTF and BED sides, so versioned and unversioned IDs interoperate.

### Custom GTF

If the proteins you care about aren't in the reference — transgenes, non-reference ORFs, manually curated isoforms — add them to the GTF before indexing. The simplest case is just concatenating GTFs you already have:

```bash
cat reference.gtf custom_proteins.gtf > combined.gtf
```

Keep every `transcript_id` unique across the combined file. To inject custom proteins from a small table of genomic blocks instead of hand-writing GTF lines, `scripts/append_custom_proteins.py` does that (strand-aware exon numbering, one TSV row per protein) — see the script's header for the column format. Either way you end up with one GTF to index in the next step, and your custom protein IDs then behave exactly like reference ones.

## 2. Index a local GTF

Turn the GTF into a binary index with `prot2exon index`:

```bash
prot2exon index --gtf combined.gtf --out human.idx
```

| | |
|---|---|
| **Input** | a GTF file — `--gtf your.gtf` |
| **Output** | a binary index — `--out your.idx` (`--index` is an accepted alias) |

The `.idx` is a binary serialisation of the parsed GTF (chromosome names, transcript records, CDS / exon vectors, attribute lookups). The format is versioned (`INDEX_FORMAT_VERSION = 3`); loading an index built by an older prot2exon returns an explicit error asking you to rebuild, so rebuild after upgrading.

From Python, `build_index` is the mirror of `prot2exon index` — it indexes a local GTF and returns the `Path` to the `.idx`:

```python
import prot2exon as p2e
idx = p2e.build_index("combined.gtf", out="human.idx")
```

See [[Python API]] for using the resulting index programmatically.

## 3. Retrieve a pre-built index from Zenodo

To skip the GTF download + build entirely, `prot2exon fetch <target>` pulls a ready-to-use binary index from the Zenodo deposit — a single sha256-verified HTTPS download, cached in `~/.cache/prot2exon/`:

```bash
prot2exon fetch list      # see every target
prot2exon fetch human     # GENCODE v49 index   -> ~/.cache/prot2exon/human.idx
prot2exon fetch mouse     # GENCODE vM34 index
prot2exon fetch yeast     # RefSeq R64 index
```

Available pre-built indexes:

| Target | Index binary | Source annotation |
|---|---|---|
| `human` | `gencode_v49_human.idx` (~298 MB) | GENCODE v49 basic, GRCh38 — current human |
| `mouse` | `gencode_vM34_mouse.idx` (~73 MB) | GENCODE vM34 basic, GRCm39 |
| `yeast` | `refseq_R64_yeast.idx` (~1.4 MB) | NCBI RefSeq *S. cerevisiae* R64 |
| `human-v86` | `ensembl_v86_human.idx` (~87 MB) | Ensembl 86, matches `EnsDb.Hsapiens.v86` (validation) |
| `human-v115` | `ensembl_v15_human.idx` (~87 MB) | Ensembl 115 


When the pinned Zenodo release isn't what you need, point `fetch` at a different source — it then downloads that GTF and builds the index:

```bash
prot2exon fetch human --release 50                                   # a specific GENCODE release
prot2exon fetch ensembl --species danio_rerio --assembly GRCz11 --release 115
prot2exon fetch human --gtf-url https://your.host/custom.gtf.gz      # any GTF URL
```

`prot2exon fetch` at a glance:

| | |
|---|---|
| **Input** | a target name (`human`, `mouse`, …); optionally `--release` / `--gtf-url` to override the source |
| **Parameters** | `--out` (default `~/.cache/prot2exon/<target>.idx`), `--cache-dir`, `--force`, `--keep-gtf` |
| **Output** | a ready `.idx`; the path is printed on stdout so it pipes into the mapper |

From Python, the same retrieval returns the `Path`:

```python
import prot2exon as p2e
idx = p2e.fetch_index("human")                # pre-built, from Zenodo
idx = p2e.fetch_index("human", release="50")  # or build a specific release
```
