# Mapping (`fastCDS map`)

`fastCDS map` projects protein-domain amino-acid coordinates onto a transcript's genomic CDS / UTR / intron structure. You give it a query BED and an index (or GTF), pick how you want the answer shaped with `--output`, and it writes a set of TSVs and BEDs describing where the domain lives in the genome.

## Input: the query BED

fastCDS reads whitespace-separated, BED-like text. Lines starting with `#` and blank lines are ignored.

```
# rows with a domain (ENSP or ENST, with or without version)
ENSP00000269305      10    50   AD1     TF1
ENST00000269305.9    10    50   AD1_ENST  TF1   # same answer as the line above

# row without a domain (whole-transcript structure only)
ENSP00000418960
```

| # | Required | Meaning |
|---|---|---|
| 1 | yes | `id` — `ENSP*` (protein) or `ENST*` (transcript). Versioned or unversioned; the version suffix is stripped on both sides before lookup. RefSeq (`NP_*`, `NM_*`) and custom IDs also work as long as the GTF/index has them. |
| 2 | no  | `aa_start` — 1-based inclusive. Omit (or set to 0) for no-domain mode. |
| 3 | no  | `aa_end` — 1-based inclusive. |
| 4 | no  | `domain_id` (used as `input_id` for tracking through the outputs). |
| 5+ | no | Ignored — free space for human-readable metadata. |

Without column 4, `input_id` falls back to `id:aa_start-aa_end` (or just `id` in no-domain mode).

An ENST resolves to the same intervals as its matching ENSP, so the two produce identical mapping output; the summary's `input_id_type` column records which form you supplied. A *non-coding* ENST (no CDS in the GTF) is reported in `unmapped_domains.tsv` with `reason = no_CDS_for_protein`. A row that is just an id with no aa range is processed in **no-domain (structure-only) mode**: the overlap columns are `NA` and the companion BEDs are empty, but `isoform_structure.tsv` is still fully populated so the transcript can be plotted as-is.

The index that `map` reads is built with `fastCDS index` — see [[Index]]. (Building a query BED from a domain database such as Pfam or InterProScan? See the [[FAQ]].)

## Mapping a domain

Every `--output` mode describes the *same* mapping; they are complementary views of where a domain sits in the genome. Pick the one that answers your question, or use the default `all` to get everything (see [Output description](#output-description)).

### Coding exons (`--output coding`)

Answers *which CDS exons code the domain?* — every CDS exon of the transcript, classified by whether it overlaps the domain.

| Output File | Contents |
|---|---|
| `domain_cds_segments.tsv` | one row per CDS exon of the transcript, each with an `overlaps_domain` column |
| `domain_cds_segments.bed` | subset: only the CDS exons that code the domain (`coding_overlap`) |

```bash
fastCDS map --index human.idx --bed q.bed --out-dir results --output coding
```

The same in Python:

```python
import fastCDS as fc
mapper = fc.Mapper(index="human.idx")
result = mapper.map_batch(
    [{"protein_id": "ENSP00000269305", "aa_start": 102, "aa_end": 292, "domain_id": "TP53_DBD"}],
    output="coding",
)
result.cds_segments   # DataFrame of CDS rows
```

`overlaps_domain` is `coding_overlap` for a CDS exon that codes the domain and `no` otherwise (`NA` in no-domain mode); partial-overlap rows also report the coding sub-interval in genomic / CDS-nt / aa space. Full column layout in [Output description](#output-description).

### Introns (`--output introns`)

Answers *which introns fall inside the domain's genomic span?* — every intron of the transcript, classified by whether it lies within the domain envelope.

| Output File | Contents |
|---|---|
| `domain_introns.tsv` | one row per intron of the transcript, each with an `overlaps_domain` column |
| `domain_introns.bed` | subset: only the introns inside the domain span (`inside_domain_genomic_span`) |

```bash
fastCDS map --index human.idx --bed q.bed --out-dir results --output introns
```

The same in Python:

```python
result = mapper.map_batch(queries, output="introns")
result.introns        # DataFrame of intron rows
```

`overlaps_domain` is `inside_domain_genomic_span` for an intron between two domain-coding CDS exons, `no` for an intron elsewhere, and `NA` in no-domain mode. Full schema in [Output description](#output-description).

### Genomic span (`--output span`)

Answers *what is the single genomic envelope of the domain, introns included?* — one interval per domain, from its first coding base to its last.

| Output File | Contents |
|---|---|
| `domain_span_with_introns.bed` | one row per domain: first → last coding base, spanning any introns in between |

```bash
fastCDS map --index human.idx --bed q.bed --out-dir results --output span
```

The same in Python:

```python
result = mapper.map_batch(queries, output="span")
result.span_bed       # DataFrame, one row per domain
```

See [Output description](#output-description) for the BED column meanings.

### Also emitting BED12 (`--bed12`)

`--bed12` is an add-on flag, not an output mode: it writes one extra file, `domain_blocks.bed12`, *in addition to* whatever `--output` you chose. Use it to get an IGV-ready BED12 alongside a narrower view such as `coding`. It is a no-op under `--output all` or `--output bed12`, which already write the file.

| Output File | Contents |
|---|---|
| `domain_blocks.bed12` | one IGV/UCSC-ready BED12 row per domain — the domain envelope drawn thick, with the coding CDS slices as blocks and the in-domain introns as gaps |

```bash
fastCDS map --index human.idx --bed q.bed --out-dir results --output coding --bed12
```

The same in Python:

```python
# The Python results carry the bed12 rows whenever they're produced:
result = mapper.map_batch(queries, output="all")
result.bed12          # DataFrame, one row per domain
```

Field-by-field meaning in [Output description](#output-description).

## Get the transcript structure (`--output isoform`)

`--output isoform` writes the plot-ready `isoform_structure.tsv`: one row per structural feature of the transcript — `five_prime_UTR`, `CDS`, `three_prime_UTR`, and inferred `intron` rows. CDS exons that the domain only partially covers are *split* into separate rows (the overlapping and non-overlapping portions), and every row carries a `plot_group` string for direct colour mapping.

```bash
fastCDS map --index human.idx --bed q.bed --out-dir results --output isoform
```

The same in Python:

```python
result = mapper.map_batch(queries, output="isoform")
result.isoform        # DataFrame, one row per structural feature
```

This is the table you feed to a renderer; see [[Plotting]] for turning it into a transcript diagram. The column layout (including `feature_id` / `feature_part` and the `plot_group` values) is in [Output description](#output-description).

## Output description

`domain_mapping_summary.tsv` is written for **every** command (one row per input query, mapped or not), and `unmapped_domains.tsv` appears whenever at least one query fails. The per-feature TSVs and BEDs depend on `--output`. This section documents every file's schema.

### The `--output all` parameter and what each mode writes

`--output` selects which feature files are produced; the default is `all`.

| `--output` | Question | TSV(s) | BED(s) |
|---|---|---|---|
| `coding`  | Where are the CDS exons, and which code the domain? | `domain_cds_segments.tsv` | `domain_cds_segments.bed` |
| `introns` | Where are the introns, and which lie in the domain span? | `domain_introns.tsv` | `domain_introns.bed` |
| `span`    | What is the domain's genomic envelope (introns included)? | — | `domain_span_with_introns.bed` |
| `isoform` | How is the whole transcript organised? | `isoform_structure.tsv` | — |
| `bed12`   | One IGV-ready BED12 row per domain. | — | `domain_blocks.bed12` |
| `all` (default) | Everything above. | all 4 TSVs | all 4 BEDs |

`all` additionally writes `run_metadata.json`. The `--bed12` flag adds `domain_blocks.bed12` to any of the first four modes. Regardless of mode, `domain_mapping_summary.tsv` is always written and `unmapped_domains.tsv` is written on any failure.

### Coordinate conventions

| Output | System |
|---|---|
| `*.bed` | 0-based half-open (BED standard). |
| `*.tsv` | 1-based inclusive (matches GTF). |

`NA` means *not applicable to this row* — e.g. CDS-nt fields on a UTR row, or overlap fields in no-domain mode. The BED and the TSV describe the same intervals in different conventions, which is why a coordinate can differ by 1 between them: a CDS at GTF `7676219..7676272` (1-based inclusive, length 54) is BED `7676218..7676272` (0-based half-open, still length 54) — `end − start` matches, only `start` shifts.

### `domain_mapping_summary.tsv`

One row per input query, written for every `--output` mode.

| Column | Type | Meaning |
|---|---|---|
| `input_id` | string | User-supplied identifier |
| `protein_id` | string | Normalised (version suffix stripped) |
| `transcript_id` | string | Transcript that this protein belongs to |
| `gene_id` | string | Ensembl gene id |
| `gene_name` | string | HGNC-style gene symbol |
| `domain_id` | string | BED column 4 (if any) |
| `chrom`, `strand` | | Chromosome and strand of the transcript |
| `aa_start`, `aa_end` | int / NA | Input domain bounds; `NA` in no-domain mode |
| `domain_length_aa` | int / NA | `aa_end − aa_start + 1` |
| `domain_length_nt` | int / NA | `domain_length_aa × 3` |
| `protein_length_aa` | int | Total CDS length / 3 for this protein |
| `domain_genomic_start` / `_end` | int / NA | Genomic envelope of the domain |
| `n_coding_segments` | int / NA | Number of CDS exon slices the domain spans |
| `fully_mapped` | bool | `true` if the entire aa range fits inside the CDS |
| `no_domain_mode` | bool | `true` if the BED row had no aa range |
| `input_id_type` | `ENSP` / `ENST` / `NA` | Which kind of id the user supplied |
| `is_mane_select` | `true` / `false` / `NA` | MANE Select transcript? `NA` if the GTF lacks `tag` attributes. |
| `is_ensembl_canonical` | `true` / `false` / `NA` | Ensembl_canonical transcript? `NA` if the GTF lacks `tag` attributes. |
| `cds_length_mismatch` | bool | `true` if `sum(CDS_nt) % 3 != 0` (Sec, readthrough, incomplete) |
| `cds_nt_remainder` | int (0, 1, 2) | `sum(CDS_nt) % 3` |
| `n_coding_exons_touched` | int / NA | Distinct CDS exons the domain overlaps (split CDS rows count once) |
| `n_introns_spanned` | int / NA | Introns between two domain-coding CDS rows |
| `is_single_exon_domain` | bool / NA | `n_coding_exons_touched == 1` |
| `fraction_domain_in_largest_exon` | float / NA | Max over CDS exons of summed `domain_overlap_fraction_of_domain`; in [0, 1] |
| `intron_burden_nt` | int / NA | Sum of lengths of the introns counted in `n_introns_spanned` |
| `status` | string | `ok` / `ok_cds_mismatch` / `partial` / `partial_cds_mismatch` / `structure_only` / unmapped reason |

### Feature TSVs (`domain_cds_segments.tsv`, `domain_introns.tsv`, `isoform_structure.tsv`)

All three share one column layout; they differ only in which feature rows they hold — `domain_cds_segments.tsv` the CDS rows, `domain_introns.tsv` the introns, `isoform_structure.tsv` every 5′UTR / CDS / 3′UTR / intron row.

| Column(s) | Meaning |
|---|---|
| `input_id`, `gene_id`, `gene_name`, `transcript_id`, `protein_id`, `domain_id` | the query and its IDs |
| `chrom`, `strand` | chromosome and `+` / `−` |
| `feature_genomic_start` / `_end`, `feature_length_nt` | genomic span (1-based inclusive) and its length (`end − start + 1`) |
| `feature_type` | `five_prime_UTR`, `CDS`, `three_prime_UTR`, or `intron` |
| `feature_id` | stable id in translation order (`CDS_1` = most 5′ CDS; UTRs/introns numbered separately). **Unchanged when a CDS is split** |
| `feature_part` | `1..K` for the pieces of a CDS split by partial domain overlap (same `feature_id`); always `1` otherwise |
| `exon_number` | source GTF `exon_number` (UTR / CDS rows); `NA` for introns |
| `feature_order_genomic` / `feature_order_transcript` | position `1..N` along the chromosome / in translation order (equal on `+` strand, reversed on `−`) |
| `cds_nt_start` / `_end`, `aa_start_encoded` / `aa_end_encoded` | this slice's CDS-relative nt offsets and the aa it encodes (`aa = ⌈cds_nt / 3⌉`); `NA` on UTR / intron rows |
| `overlaps_domain` | `coding_overlap` (CDS interval codes the domain), `inside_domain_genomic_span` (intron between two coding rows), `no` (outside the domain; UTRs are always `no`), or `NA` (no-domain query) |
| `domain_overlap_genomic_start` / `_end`, `domain_overlap_cds_nt_start` / `_end`, `domain_overlap_aa_start` / `_end` | the sub-interval of this row that codes the domain, in genomic / CDS-nt / aa coordinates (filled only on `coding_overlap` rows) |
| `domain_overlap_fraction_of_feature` / `_of_domain` | overlap ÷ feature length, and overlap ÷ domain length (the latter sums to `1.0` across a domain's coding rows) |
| `plot_group` | single string for colour mapping: `CDS_domain`, `CDS_no_domain`, `intron_domain_span`, `five_prime_UTR`, `three_prime_UTR`, `intron` (just the feature type in no-domain mode) |

A CDS exon that a domain only partially covers is split into rows that share a `feature_id` and differ by `feature_part` (e.g. a domain ending mid-`CDS_2` gives `CDS_2` part 1 = `coding_overlap`, part 2 = `no`). Group by `(input_id, feature_id)` to re-aggregate the original exon. For `−` strand genes, plot on `feature_genomic_start/end` but read 5′→3′ order from `feature_order_transcript`.

### Companion BEDs

All three are 6-column standard BED (`chrom`, `start_0based`, `end`, `name`, `score=0`, `strand`).

| File | Rows |
|---|---|
| `domain_cds_segments.bed` | CDS rows with `overlaps_domain == coding_overlap` |
| `domain_introns.bed` | intron rows with `overlaps_domain == inside_domain_genomic_span` |
| `domain_span_with_introns.bed` | one row per domain (genomic envelope, introns included) |

`name` is `protein_id[_domain_id]_<aa_start>-<aa_end>`. No-domain queries contribute zero rows to any companion BED.

### `domain_blocks.bed12`

One BED12 row per domain. The whole feature is drawn thick in IGV; blocks are the CDS slices that code the domain. Empty in no-domain mode.

| BED12 field | Meaning here |
|---|---|
| `chrom` | chromosome |
| `chromStart` (0-based) | first base of the domain genomic envelope |
| `chromEnd` | last base + 1 of the envelope (so `end − start` = envelope length, introns included) |
| `name` | `protein_id[_domain_id]_<aa_start>-<aa_end>` |
| `score` / `strand` | `0` / transcript strand |
| `thickStart` / `thickEnd` | equal to `chromStart` / `chromEnd` — IGV draws the whole feature thick |
| `itemRgb` | `255,0,0` (red) |
| `blockCount` | number of CDS slices that code the domain (`coding_overlap` rows) |
| `blockSizes` / `blockStarts` | comma-separated, in genomic order (starts are offsets from `chromStart`) |

### `unmapped_domains.tsv`

Written only when at least one query failed.

| Column | Meaning |
|---|---|
| `input_id`, `protein_id`, `aa_start`, `aa_end`, `domain_id` | identity (`protein_id = NA` for non-coding ENST queries) |
| `reason` | one of: `protein_not_in_index`, `no_CDS_for_protein`, `domain_beyond_protein_length`, `no_overlap` |

### `run_metadata.json`

Written only with `--output all`. Records `tool` / `version` / `timestamp_utc`, `output_kind`, `annotation_source` (the GTF or index used), `index_format_version`, `coordinate_conventions`, `query_counts` (`{ total, mapped, unmapped, no_domain_mode }`), and `cli` (the full invocation).

### MANE Select / Ensembl Canonical

GENCODE GTFs (since v34) carry `tag "MANE_Select"` and `tag "Ensembl_canonical"`. The tool surfaces these in `is_mane_select` / `is_ensembl_canonical` (`true` / `false` / `NA`) on `domain_mapping_summary.tsv` and every feature TSV. `NA` is reserved for GTFs that carry no `tag` attribute anywhere (typical of base Ensembl) — distinct from `false` ("tags exist, but not on this transcript"). Filter to MANE Select queries with:

```bash
awk -F'\t' 'NR==1 || $20=="true"' domain_mapping_summary.tsv   # col 20 = is_mane_select
```

### CDS-length mismatch (Sec, readthrough, incomplete CDS)

The tool maps every query but flags the rare cases where the CDS isn't a multiple of 3: `cds_length_mismatch` (`true` if `sum(CDS_nt) % 3 != 0`), `cds_nt_remainder` (`0`/`1`/`2`), and a `_cds_mismatch` suffix on `status`. This catches selenoproteins (internal `TGA` re-coded as Sec, plus a C-terminal extension), stop-codon readthrough, and incomplete `cds_start_NF` / `cds_end_NF` transcripts. The aa↔nt math uses ceiling division (`aa = ⌈cds_nt / 3⌉`), so a domain at the very C-terminus of an incomplete CDS may be clipped by 1–2 aa; the flag is your hint to check.

### Codons split across exons

A codon that straddles an exon boundary is mapped correctly in both phases (1+2 and 2+1) — the math is cumulative-nt-based, not exon-by-codon. For a domain on the split aa, both the upstream and downstream CDS slices report a `coding_overlap` row, the intron between them is marked `inside_domain_genomic_span`, and the two `domain_overlap_fraction_of_domain` values (≈1/3 and ≈2/3) sum to `1.0`. Two synthetic cases in `software_tests/` exercise this.

### Worked example

Input BED:

```
ENSP00000306245    5    100    RD1    TF2
```

`ENSP00000306245` is a `+` strand transcript; domain `RD1` (aa 5..100) spans the end of one CDS slice and most of the next. The relevant rows of `isoform_structure.tsv` (selected columns):

| feature_type | feature_id | feature_part | feature_genomic_start | feature_genomic_end | aa_start_encoded | aa_end_encoded | overlaps_domain | plot_group |
|---|---|---|---|---|---|---|---|---|
| CDS | CDS_1 | 1 | 75278983 | 75278994 | 1 | 4 | no | CDS_no_domain |
| CDS | CDS_2 | 1 | 75278995 | 75279123 | 5 | 47 | coding_overlap | CDS_domain |
| intron | intron_1 | 1 | 75279124 | 75279876 | NA | NA | inside_domain_genomic_span | intron_domain_span |
| CDS | CDS_3 | 1 | 75279877 | 75280035 | 48 | 100 | coding_overlap | CDS_domain |
| CDS | CDS_4 | 1 | 75280036 | 75280128 | 101 | 131 | no | CDS_no_domain |
| three_prime_UTR | three_prime_UTR_1 | 1 | 75281422 | 75282230 | NA | NA | no | three_prime_UTR |

To plot, fill by `plot_group`; to highlight the domain, take `plot_group ∈ {CDS_domain, intron_domain_span}`; to get just the domain-coding bases, take `domain_cds_segments.bed`. The same id on its own line (no aa columns) emits the whole-transcript rows with `overlaps_domain = NA`.
