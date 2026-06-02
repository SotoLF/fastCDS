# Building an index (`prot2exon index`, `prot2exon fetch`)

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

#### GTF Concatenation

If the proteins you care about aren't in the reference — transgenes, non-reference ORFs, manually curated isoforms — add them to the GTF before indexing. The simplest case is just concatenating GTFs you already have:

```bash
cat reference.gtf custom_proteins.gtf > combined.gtf
```

Make sure to keep every `transcript_id` unique across the combined file. 

#### Input TSV

To inject custom proteins from a small table of genomic blocks instead of hand-writing GTF lines, `scripts/append_custom_proteins.py` does strand-aware exon numbering: one TSV row per transcript.


A tab-separated table with one row per transcript:

```
protein_id   transcript_id   gene_id  gene_name  chrom    strand  blocks
NP_NOVEL_1   NM_NOVEL_1      G_NOV1   NOVEL_1    chr_X    +       100-150;200-280;350-410
NP_NOVEL_2   NM_NOVEL_2      G_NOV2   NOVEL_2    chr_X    -       5000-5100;4800-4900;4600-4700
```

| Column | Meaning |
|---|---|
| `protein_id`, `transcript_id`, `gene_id`, `gene_name` | IDs the rest of the pipeline will see. Pick whatever scheme you like; just keep them unique. |
| `chrom`, `strand` | Genomic placement. Strand drives exon numbering. |
| `blocks` | Semicolon-separated `start-end` genomic ranges (1-based inclusive, GTF style), written in genomic order (ascending `start`). The script assigns exon numbers strand-aware: on `-` strand the highest-coordinate block becomes exon 1. |

#### Run

The script emits the custom `transcript` / `exon` / `CDS` rows to stdout — append them to a copy of your reference GTF, then index it:

```bash
cp gencode.v49.primary_assembly.annotation.gtf combined.gtf
python3 scripts/append_custom_proteins.py --in my_custom_proteins.tsv >> combined.gtf
```

Use `--out combined_rows.gtf` to write to a file instead of stdout, and `--source-tag <text>` to set the GTF `source` column for these rows (default `custom`).

Either way you end up with one GTF to index in the next step, and your custom protein IDs then behave exactly like reference ones.

## 2. Index a local GTF

Turn the GTF into a binary index with `prot2exon index`:

```bash
prot2exon index --gtf gencode.v49.primary_assembly.annotation.gtf --out human.idx
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

`prot2exon fetch <target>` downloads a ready-to-use binary index straight from the Zenodo deposit — one sha256-verified HTTPS download — so you skip both the GTF download **and** the build. Point `--out` at wherever you want the `.idx`:

```bash
prot2exon fetch list                     # see every target
prot2exon fetch human --out human.idx    # pre-built GENCODE v49 index -> ./human.idx
prot2exon fetch mouse --out mouse.idx    # GENCODE vM34 index
prot2exon fetch yeast --out yeast.idx    # RefSeq R64 index
```

Without `--out` the index lands in the cache (`~/.cache/prot2exon/<target>.idx`) and that path is printed on stdout, so it still pipes straight into the mapper. Re-running reuses the cached file — pass `--force` to re-download.

Available pre-built indexes:

| Target | Index binary | Source annotation |
|---|---|---|
| `human` | `gencode_v49_human.idx` (~298 MB) | GENCODE v49 basic, GRCh38 — current human |
| `mouse` | `gencode_vM34_mouse.idx` (~73 MB) | GENCODE vM34 basic, GRCm39 |
| `yeast` | `refseq_R64_yeast.idx` (~1.4 MB) | NCBI RefSeq *S. cerevisiae* R64 |
| `human-v86` | `ensembl_v86_human.idx` (~87 MB) | Ensembl 86, matches `EnsDb.Hsapiens.v86` (validation) |


`prot2exon fetch` at a glance:

| | |
|---|---|
| **Input** | a target name (`human`, `mouse`, …) |
| **Parameters** | `--out` (where to write the `.idx`; default `~/.cache/prot2exon/<target>.idx`), `--cache-dir`, `--force` |
| **Output** | a ready `.idx`; its path is printed on stdout so it pipes into the mapper |

**Need a release that isn't in the table above?** There's no pre-built index to fetch, so make one the normal way: [obtain that GTF](#1-obtain-a-gtf) and then [index it](#2-index-a-local-gtf) — that's exactly what `fetch` saves you when a pre-built index *does* exist. (For convenience `fetch` can also download-and-build a specific release in one step via `--release` / `--gtf-url`; see `prot2exon fetch --help`. It's the same obtain-GTF-then-index route, automated.)

From Python, the same retrieval returns the `Path`:

```python
import prot2exon as p2e
idx = p2e.fetch_index("human")    # pre-built index, from Zenodo
```
