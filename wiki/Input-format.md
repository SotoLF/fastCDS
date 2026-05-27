# Input format

prot2exon eats whitespace-separated, BED-like text. Lines starting with `#` are ignored, as are blank lines.

```
# rows with a domain (ENSP or ENST, with or without version)
ENSP00000269305      10    50   AD1     TF1
ENST00000269305.9    10    50   AD1_ENST  TF1   # same answer as the line above

# row without a domain (whole-transcript structure only)
ENSP00000418960
```

## Columns

| # | Required | Meaning |
|---|---|---|
| 1 | yes | `id` — `ENSP*` (protein) or `ENST*` (transcript). Versioned or unversioned; the suffix is stripped on both sides. RefSeq (`NP_*`, `NM_*`) and custom IDs also work as long as the GTF has them. |
| 2 | no  | `aa_start` — 1-based inclusive. Omit (or set to 0) for no-domain mode. |
| 3 | no  | `aa_end` — 1-based inclusive. |
| 4 | no  | `domain_id` (used as `input_id` for tracking through outputs). |
| 5+ | no | Ignored — free space for human-readable metadata. |

If column 4 is missing and a domain is present, `input_id` falls back to `id:aa_start-aa_end`. In no-domain mode it falls back to `id`. Every row stays uniquely identifiable in the outputs.

## ENST vs ENSP

- ENST resolves to the same intervals as the matching ENSP — the two give identical mapping output.
- A *non-coding* ENST (no CDS records in the GTF) is reported in `unmapped_domains.tsv` with `reason = no_CDS_for_protein` and `protein_id = NA`.
- The summary's `input_id_type` column records which form was supplied (`ENSP` / `ENST`).

## No-domain mode

A row that's just `ENSP00000418960` (no aa range) is processed as **whole-transcript structure with no domain**:

- The overlap columns are `NA`.
- The companion BEDs are empty for that row.
- The `isoform_structure.tsv` table is still populated, so the transcript can be plotted as-is.

This is useful for plotting transcript architecture before any domain annotation exists.

## Where do BEDs of aa ranges come from?

Many real-world workflows start from a domain database. Helpers in `scripts/` convert the common formats:

| Source | Helper |
|---|---|
| InterProScan TSV | `scripts/prepare_from_interpro.py` |
| UniProt feature table (REST JSON or `.dat`) | `scripts/prepare_from_uniprot_features.py` |
| HMMER `--domtblout` (e.g. Pfam scan) | `scripts/prepare_from_pfam.py` |

All three emit the BED-like format above, including the UniProt-accession → ENSP cross-reference step where needed.

## Unannotated proteins

If your protein doesn't exist in any standard GTF (transgenes, non-reference ORFs, manually curated isoforms…), see [[Custom proteins]] — there's an `append_custom_proteins.py` helper that injects custom CDS rows into an existing GTF before you build the index.
