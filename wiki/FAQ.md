# FAQ

## "index already exists" — what do I do?

This is informational, not an error. `fastCDS fetch` caches indexes under `~/.cache/fastCDS/`; on the second invocation it just prints `using cached index: /path/to/yeast.idx  (--force to rebuild)` and prints the path. Pass `--force` if you actually want to rebuild from the upstream GTF.

## My TSV row has `cds_length_mismatch = true` — bug?

Usually not. It flags cases where the sum of CDS lengths isn't a clean multiple of 3, which happens legitimately for:

- **Selenoproteins (Sec)** — UGA stop codon read through as selenocysteine.
- **Programmed readthrough / stop-codon recoding.**
- **Incomplete CDS** — transcripts annotated as `cds_start_NF`, `cds_end_NF`, or `mRNA_end_NF`.

`cds_nt_remainder` (1 or 2) tells you which side. The mapping is still emitted; you just lose codon-precision at the affected end.

## A query has `reason = no_CDS_for_protein` — what's going on?

You probably passed a transcript ID for a *non-coding* transcript (lncRNA, processed transcript, …). There's no CDS to map onto, so the row goes to `unmapped_domains.tsv`. The summary's `input_id_type` column distinguishes ENSP / ENST queries.

## How do I build a query BED from a domain database (Pfam, InterProScan, UniProt)?

Use the `fastCDS.prepare` helpers — `from_pfam()`, `from_interproscan()`, `from_uniprot_features()` — which turn those tools' standard outputs into a DataFrame that `Mapper.map_batch()` consumes directly (the `parsing/` CLI wrappers do the same from the command line). See the [[Python API]] page for the exact signatures.

## How does fastCDS handle codons split across an exon boundary?

Codons that straddle two CDS exons (1+2 or 2+1 patterns) are handled correctly — the mapper tracks `cds_nt_start` / `cds_nt_end` across the exon boundary so aa coordinates round-trip cleanly. CDS exons are split only when *the domain* partially covers them, never on the codon boundary.

## MANE Select and Ensembl canonical — where are they?

`is_mane_select` and `is_ensembl_canonical` columns are present on every TSV. They're booleans (`true` / `false`) parsed from GENCODE/Ensembl `tag` annotations. RefSeq GTFs don't carry MANE Select tags, so both columns are `NA` there.

## The plot looks empty / I see no UTRs

Two common causes:

- You loaded `domain_cds_segments.tsv` instead of `isoform_structure.tsv`. The CDS table doesn't include UTR rows by design — only the full isoform table does.
- You passed `--no-utr` to the plotter, or unchecked **Show UTRs** in the interactive viewer.

## The viewer shows long stretches of empty intron — fix?

Use compact mode. CLI: `--compact-genomic` (matplotlib) or the **Compact (introns = 80 bp)** radio in the interactive viewer. Long human genes (TP53 spans ~19 kb) become legible at exon resolution.

## My `--engine plotly` HTML file is huge

That's expected — plotly bundles its full JS engine in every standalone HTML (~3 MB). If size matters, drop `--engine plotly` and use the default `js` engine (`--out x.html`) — that viewer is ~35 KB and has no JS dependency.

## Can I run fastCDS on Windows?

The Python wrapper works on Windows. The C++ binary should build cleanly under MSVC 2019+ or MinGW, but is mainly tested on Linux. WSL is the path of least friction.

## How do I cite fastCDS?

A formal citation will land when the accompanying manuscript is posted. Until then, please cite the repository URL.

## Where do I report a bug?

Open an issue on the GitHub repo. A minimal reproducer (the GTF lines or BED row that triggers it, plus the exact CLI you ran) makes triage 10× faster.

## I want to contribute — where do I start?

- **Bug fixes** — open a PR with a regression test under `software_tests/` (pytest; pick the module that fits — `test_correctness.py` for mapping output, `test_compat.py` for GTF dialects, etc.).
- **New plot styles** — both matplotlib (`plot.py`) and the JS viewer (`_interactive_html.py`) are templated; copy an existing draw function and add a CLI flag.
- **New input adapters** — `parsing/prepare_from_*.py` is the conventional spot for "turn external format X into fastCDS BED".
- **New genome onboarding presets** — extend `_resolve_url` in `python/fastCDS/fetch.py`.
