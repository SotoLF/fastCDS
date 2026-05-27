# Adding custom proteins to an existing GTF

Sometimes the proteins you care about aren't in the reference annotation — transgenes, non-reference ORFs, manually curated alternative isoforms, novel CDS from a focused RNA-seq run. The `scripts/append_custom_proteins.py` helper injects custom CDS rows into an existing GTF before you build the index, so the rest of the pipeline (`map_batch`, plotting, BED12 output, …) keeps working unchanged.

## Input TSV

```
protein_id   transcript_id   gene_id  gene_name  chrom    strand  blocks                              source
NP_NOVEL_1   NM_NOVEL_1      G_NOV1   NOVEL_1    chr_X    +       100-150,200-280,350-410             custom
NP_NOVEL_2   NM_NOVEL_2      G_NOV2   NOVEL_2    chr_X    -       5000-5100,4800-4900,4600-4700       custom
```

| Column | Meaning |
|---|---|
| `protein_id`, `transcript_id`, `gene_id`, `gene_name` | IDs the rest of the pipeline will see. Pick whatever scheme you like; just keep them unique. |
| `chrom`, `strand` | Genomic placement. Strand drives exon numbering. |
| `blocks` | Comma-separated `start-end` genomic ranges (1-based inclusive, GTF style). |
| `source` | Free text — surfaces in the GTF as the `source` column for these rows. |

## Run

```bash
python3 scripts/append_custom_proteins.py \
    --gtf gencode.v49.primary_assembly.annotation.gtf \
    --tsv my_custom_proteins.tsv \
    --out gencode.v49.with_custom.gtf
```

Then build an index against the augmented GTF:

```bash
./build/prot2exon --gtf gencode.v49.with_custom.gtf --build-index --index custom.idx
```

…and map BED queries that reference your custom protein IDs:

```
NP_NOVEL_1   1   50   NOVEL_1_signal_peptide
NP_NOVEL_2   100 300  NOVEL_2_envelope_domain
```

## Strand-aware exon numbering

`blocks` are written in genomic order — *increasing* `start` first, regardless of strand. The script handles strand-aware numbering automatically:

- **`+` strand** — block 0 (lowest coord) becomes exon 1.
- **`−` strand** — block N (highest coord) becomes exon 1.

So `5000-5100,4800-4900,4600-4700` on the `−` strand still parses as a valid 3-exon CDS with exon 1 = `5000-5100`.

## Validations the script enforces

- `strand` must be `+` or `−`.
- Every block must contain a `-`; both ends must parse to int.
- `end >= start` for every block.
- `protein_id` and `transcript_id` must be unique within the TSV.

Failures are reported with the offending TSV row number.

## Gotchas worth knowing

- The `chrom` value in your custom rows must match how it will be referenced downstream (IGV / the fasta you'd visualise against). Use a non-colliding name if your custom locus isn't on a real reference chromosome.
- Synthetic / non-canonical CDSs may not be a multiple of 3 (programmed frameshifts, premature stops, etc.). The mapper handles this — the row is still mapped, but `cds_length_mismatch = true` and `status` becomes `ok_cds_mismatch`. That's a feature, not a bug.
- Each `transcript_id` must be unique across the *combined* GTF (canonical + your additions). If you re-run with a tweaked custom file but reuse old IDs, your existing index is stale silently — rebuild it.
- The same `--build-index` command works regardless of whether the GTF has custom rows — the parser auto-handles versioned / unversioned IDs and `gene_name` / `gene` symbol attributes.

## Workflow recap

1. Identify each protein's genomic coordinates (manual curation, alignment output, or an alternative GFF).
2. Write a TSV with one row per protein.
3. Run `append_custom_proteins.py` to fold the rows into the GTF as new transcripts.
4. Build an index on the augmented GTF.
5. Map BED queries against your custom protein IDs exactly like reference proteins — same output schemas, same plotters.
