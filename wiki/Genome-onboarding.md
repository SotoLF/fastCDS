# Onboarding a new genome

There are two paths: **download a pre-built index** (instant, from our Zenodo deposit) or **build from a GTF** (handles any species).

## Fastest path: pre-built index from Zenodo

Skip the GTF parse + build (which can take ~20 s for human, less for smaller organisms) by downloading a pre-built `.idx` straight from the Zenodo data deposit.

```bash
# See what's available
prot2exon fetch list

# Download by preset name — URL + sha256 are baked in, --out defaults to ~/.cache/prot2exon/<preset>.idx
prot2exon fetch index --preset human-v49

# Or pass a custom URL (any HTTPS / Zenodo / file:// source)
prot2exon fetch index \
    --url https://zenodo.org/record/<RECORD>/files/gencode_v49_human.idx \
    --out ~/.cache/prot2exon/human.idx \
    --sha256 ed848d78125dc795fa86a0af5402cb08ad679626fb153dda7a8ff2d6b47844f7
```

Re-running with `--out` already populated returns the cached path (pass `--force` to redownload). With `--url`, the `--sha256` flag is optional but strongly recommended — it verifies the downloaded file against the manifest.

Available pre-built indexes in the Zenodo bundle (full list + checksums in the deposit's `MANIFEST.tsv`):

| File | Source annotation |
|---|---|
| `gencode_v49_human.idx` | GENCODE v49 basic, GRCh38 (current human) |
| `gencode_vM34_mouse.idx` | GENCODE vM34 basic, GRCm39 |
| `ensembl_v86_human.idx` | Ensembl 86 — matches `EnsDb.Hsapiens.v86` (for validation reproducibility) |
| `refseq_R64_yeast.idx` | NCBI RefSeq *S. cerevisiae* R64 |

## One-command rebuild from upstream: `prot2exon fetch`

```bash
prot2exon fetch list                 # show what's built-in
prot2exon fetch human --release 49   # GENCODE human 49
prot2exon fetch mouse --release M34  # GENCODE mouse M34
prot2exon fetch ensembl --species homo_sapiens --release 110
prot2exon fetch refseq --preset yeast
```

The script:

1. Downloads the right `.gtf.gz` to `~/.cache/prot2exon/` (override with `--cache-dir`).
2. Gunzips it.
3. Runs the binary's `--build-index`.
4. Prints the resulting `.idx` path on stdout (so you can pipe it into shell scripts).

Re-running is a no-op — already-cached downloads and indexes are reused. Pass `--force` to rebuild.

Override the URL if you need a release prot2exon doesn't know about:

```bash
prot2exon fetch human --gtf-url https://ftp.ensembl.org/.../Homo_sapiens.GRCh38.110.gtf.gz
```

## Manual recipes

If you'd rather drive it yourself or are on a corp network that blocks the upstream URLs:

```bash
# GENCODE human
curl -O https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.primary_assembly.annotation.gtf.gz
gunzip gencode.v49.primary_assembly.annotation.gtf.gz
./build/prot2exon --gtf gencode.v49.primary_assembly.annotation.gtf --build-index --index gencode_v49.idx

# Ensembl human
curl -O https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz
gunzip Homo_sapiens.GRCh38.110.gtf.gz
./build/prot2exon --gtf Homo_sapiens.GRCh38.110.gtf --build-index --index ensembl_110.idx

# NCBI RefSeq (Yeast)
curl -O https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.gtf.gz
gunzip GCF_000146045.2_R64_genomic.gtf.gz
./build/prot2exon --gtf GCF_000146045.2_R64_genomic.gtf --build-index --index yeast.idx
```

## GTF compatibility

prot2exon handles all three major GTF dialects:

| Source | What we use |
|---|---|
| **GENCODE** | `gene_name`, `protein_id`, `transcript_id`, `tag` (MANE Select, Ensembl_canonical, CCDS, etc.) |
| **Ensembl**  | Same as GENCODE — they share the format. |
| **NCBI RefSeq** | `gene` (instead of `gene_name`), `protein_id`, `transcript_id`. MANE-Select tags are absent — those columns are `NA`. |

Two specific parser bugs to be aware of (both fixed):

- **RefSeq's `gene "PAU8"` vs GENCODE's `gene_name "PAU8"`** — the parser now falls back to `gene` when `gene_name` is empty, so `gene_name` is populated for both dialects.
- **Substring-vs-attribute matching on `tag`** — the parser uses word-boundary attribute matching, so `locus_tag "..."` no longer gets misread as a `tag` annotation.

If you find an annotation source where prot2exon doesn't extract the expected columns, please open an issue with a short reproducer (a few GTF lines is usually enough).

## Parser requirements

- `exon` and `CDS` feature lines (other features like `gene`, `transcript`, `start_codon` are read past).
- Attributes in standard GTF form: `key "value";` (GFF3's `key=value;...` does **not** work).
- `gene_id` and `transcript_id` on every relevant line. `protein_id` is optional but enables ENSP / NP lookups.

The tool normalises IDs by stripping the `.<version>` suffix on both sides (GTF and BED), so a BED with versioned IDs works against an unversioned index, and vice versa.

## Index format

The `.idx` is a binary serialisation of the parsed GTF — chromosome names, transcript records, CDS / exon vectors, attribute lookups. The format is versioned (`INDEX_FORMAT_VERSION = 3`); loading an older index returns an explicit error asking you to rebuild it.

Version 3 adds:

- the transcript → protein reverse map (so ENST queries resolve in the same lookup pass);
- per-transcript MANE Select / Ensembl Canonical flags.

On GTFs that don't carry `tag "..."` attributes at all (e.g. base Ensembl, RefSeq), those two flag columns report `NA` instead of `false`.
