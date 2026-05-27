# Output modes

Pass `--output KIND` together with `--out-dir DIR`. Every mode also writes `domain_mapping_summary.tsv` (one row per input, ok or not) and `unmapped_domains.tsv` (only when at least one row failed).

| `--output` | Question | TSV(s) | BED(s) |
|---|---|---|---|
| `coding`  | Where are the CDS exons of this transcript, and which overlap the domain? | `domain_cds_segments.tsv` | `domain_cds_segments.bed` (subset: `coding_overlap` rows) |
| `introns` | Where are the introns, and which lie within the domain span? | `domain_introns.tsv` | `domain_introns.bed` (subset: `inside_domain_genomic_span` rows) |
| `span`    | What is the single genomic envelope of the domain (introns included)? | — | `domain_span_with_introns.bed` (1 row / domain) |
| `isoform` | How is the whole transcript organised, and where does the domain fall on it? | `isoform_structure.tsv` | — |
| `bed12`   | One IGV-ready BED12 row per domain (blocks = CDS slices coding the domain). | — | `domain_blocks.bed12` |
| `all` (default) | Everything above. | all 4 TSVs | all 4 BEDs + `run_metadata.json` |

## What's in each TSV

The full schemas live in the [README](https://github.com/SotoLF/Prot2Exon/blob/main/README.md#file-by-file-schemas); the quick mental model:

### `domain_cds_segments.tsv`

One row per CDS exon of the transcript. Every row has an `overlaps_domain` column:

- `coding_overlap` — this CDS exon partially or fully codes the domain.
- `no_overlap` — this CDS exon is in the same transcript but not in the domain.
- `NA` — no-domain mode.

For a CDS exon that the domain only partially covers, the row also reports the overlap as `domain_overlap_{genomic,cds_nt,aa}_{start,end}` plus `domain_overlap_fraction_of_feature` and `…_of_domain`.

### `domain_introns.tsv`

Same shape, but rows are introns. `overlaps_domain` here means:

- `inside_domain_genomic_span` — this intron sits between two CDS exons that code the domain.
- `outside` — this intron is in the transcript but outside the domain envelope.
- `NA` — no-domain mode.

### `isoform_structure.tsv`

The plot-ready table. One row per structural feature:

| `feature_type` | What it means |
|---|---|
| `five_prime_UTR` | 5′ UTR exon |
| `CDS` | a CDS exon (or a *part* of a CDS exon — see below) |
| `three_prime_UTR` | 3′ UTR exon |
| `intron` | inferred between exons |

**CDS rows that the domain partially covers are split** so the domain-overlapping portion and the non-overlapping portion are separate rows. The split rows keep the original CDS number — only `feature_part` distinguishes them. `plot_group` (`CDS_no_domain`, `CDS_domain`, `intron_domain_span`, …) is what the plotter uses for colouring.

### `domain_blocks.bed12`

One IGV-ready row per domain. The 12-column BED's blocks are the CDS slices that code the domain — so dropping the file into IGV shows the domain as a multi-block alignment block on the gene.

## Coordinate conventions

| Output | System |
|---|---|
| `*.bed` | 0-based half-open (BED standard). |
| `*.tsv` | 1-based inclusive (matches GTF). |

`NA` means *not applicable to this row* (e.g. CDS-nt fields on a UTR row, or overlap fields in no-domain mode).

> **Why does the BED differ by 1 from the TSV?** Same interval, different convention. A CDS at GTF positions `7676219..7676272` (1-based inclusive, length 54) is BED `7676218..7676272` (0-based half-open, length still 54). `end − start` matches in both; only `start` shifts.

---

## File-by-file schemas

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
| `n_introns_spanned` | int / NA | Introns between two domain-coding CDS rows (i.e. with `inside_domain_genomic_span`) |
| `is_single_exon_domain` | bool / NA | `n_coding_exons_touched == 1` |
| `fraction_domain_in_largest_exon` | float / NA | Max over CDS exons of summed `domain_overlap_fraction_of_domain`; in [0, 1] |
| `intron_burden_nt` | int / NA | Sum of lengths of the introns counted in `n_introns_spanned` |
| `status` | string | `ok` / `ok_cds_mismatch` / `partial` / `partial_cds_mismatch` / `structure_only` / unmapped reason |

### Feature TSVs (`domain_cds_segments.tsv`, `domain_introns.tsv`, `isoform_structure.tsv`)

All three share the same column layout — they differ only in which feature types they include:

| File | Rows |
|---|---|
| `domain_cds_segments.tsv` | every CDS row |
| `domain_introns.tsv`      | every intron row |
| `isoform_structure.tsv`   | every 5′UTR / CDS / 3′UTR / intron row |

You can therefore `cat` or `join` them on `input_id` / `feature_id` interchangeably.

**Identity columns:** `input_id`, `gene_id`, `gene_name`, `transcript_id`, `protein_id`, `domain_id`.

**Location columns:**

| Column | Meaning |
|---|---|
| `chrom` | Chromosome |
| `strand` | `+` or `−` |
| `feature_genomic_start` / `_end` | 1-based inclusive |
| `feature_length_nt` | `end − start + 1` |

**Feature-type columns:**

| Column | Meaning |
|---|---|
| `feature_type` | One of `five_prime_UTR`, `CDS`, `three_prime_UTR`, `intron` |
| `feature_id` | **Stable across splits.** Numbered in translation order: `CDS_1` is the most 5′ CDS, then `CDS_2`, etc. UTRs and introns are numbered separately (`five_prime_UTR_1`, `intron_1`, …) |
| `feature_part` | 1..K when a CDS row was split by partial domain overlap; pieces of the same original CDS share the same `feature_id`, differ by `feature_part`. Always `1` for UTR / intron rows |
| `exon_number` | Source GTF `exon_number` for UTR / CDS rows; `NA` for introns |

**CDS splitting — `feature_id` and `feature_part`.** A single GTF CDS exon can be split into multiple isoform-structure rows when the domain only covers part of it. The pieces always share the same `feature_id` — only `feature_part` differs.

Example: domain ends in the middle of `CDS_2` (translation order), so `CDS_2` produces two rows:

| feature_id | feature_part | start | end | overlaps_domain |
|---|---|---|---|---|
| CDS_2 | 1 | 75278995 | 75279120 | coding_overlap |
| CDS_2 | 2 | 75279121 | 75279123 | no |

To re-aggregate the full original CDS, group by `(input_id, feature_id)`. To plot the overlap shape, fill by `plot_group` and ignore `feature_part`.

**Ordering columns:**

| Column | Meaning |
|---|---|
| `feature_order_genomic` | 1..N along the chromosome (low → high coord). Always equals row position. |
| `feature_order_transcript` | 1..N in translation direction. Equals `feature_order_genomic` on `+` strand; reversed on `−` strand. |

For `−` strand genes, plot using `feature_genomic_start/end` on the X axis (genomic coords), but use `feature_order_transcript` to interpret biological order (5′ → 3′ of the protein).

**CDS-coordinate columns** (NA on UTR / intron rows):

| Column | Meaning |
|---|---|
| `cds_nt_start` | CDS-relative nt offset (1-based) of this slice's first base |
| `cds_nt_end` | CDS-relative nt offset of the last base |
| `aa_start_encoded` | First aa that this slice encodes (1-based) |
| `aa_end_encoded` | Last aa that this slice encodes |

The mapping is `aa = ⌈cds_nt / 3⌉`. A 1-nt CDS slice containing only the third base of an aa still reports that aa in `aa_start_encoded` / `aa_end_encoded`.

**Domain-overlap columns.** `overlaps_domain` is **not** a yes/no flag — it discriminates *coding* overlap from *intronic* overlap inside the domain envelope:

| Value | Meaning |
|---|---|
| `no` | The row is outside the domain entirely. UTR rows always carry `no` |
| `coding_overlap` | CDS row whose genomic interval overlaps a domain-coding range |
| `inside_domain_genomic_span` | Intron located between two `coding_overlap` CDS rows |
| `NA` | No-domain query (no domain to compare against) |

The companion columns are filled only for `coding_overlap` rows:

| Column | Meaning |
|---|---|
| `domain_overlap_genomic_start` / `_end` | Sub-interval of this row that codes the domain (1-based inclusive) |
| `domain_overlap_cds_nt_start` / `_end` | Same overlap projected to CDS-nt (1-based) |
| `domain_overlap_aa_start` / `_end` | Same overlap projected to aa |
| `domain_overlap_fraction_of_feature` | `overlap_length / feature_length_nt` |
| `domain_overlap_fraction_of_domain` | `overlap_length / domain_length_nt`. Sums to `1.0` across all `coding_overlap` rows of the same domain |

**Plotting column.** `plot_group` is a single string suitable for direct colour mapping.

With a domain:

| Value | Maps to |
|---|---|
| `five_prime_UTR` / `three_prime_UTR` | UTR exon segment |
| `CDS_no_domain` | CDS outside the domain |
| `CDS_domain` | CDS that encodes the domain |
| `intron` | Intron outside the domain genomic span |
| `intron_domain_span` | Intron between two `CDS_domain` rows |

In no-domain mode, `plot_group` is just the feature type: `five_prime_UTR`, `CDS`, `three_prime_UTR`, `intron`.

ggplot / ggtranscript-style code:

```r
ggplot(rows, aes(xmin = feature_genomic_start, xmax = feature_genomic_end,
                 y = transcript_id, fill = plot_group)) +
  geom_rect()
```

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
| `score` | `0` |
| `strand` | transcript strand |
| `thickStart` / `thickEnd` | equal to `chromStart` / `chromEnd` — IGV draws the whole feature thick |
| `itemRgb` | `255,0,0` (red) |
| `blockCount` | number of CDS slices that code the domain (`coding_overlap` rows) |
| `blockSizes` | comma-separated, in genomic order |
| `blockStarts` | comma-separated offsets from `chromStart`, in genomic order |

Drop the file directly into IGV / UCSC and the domain's exonic blocks render with in-domain introns as BED12 gaps.

### `unmapped_domains.tsv`

Written only when at least one row failed.

| Column | Meaning |
|---|---|
| `input_id`, `protein_id`, `aa_start`, `aa_end`, `domain_id` | identity (`protein_id = NA` for non-coding ENST queries) |
| `reason` | one of: `protein_not_in_index`, `no_CDS_for_protein`, `domain_beyond_protein_length`, `no_overlap` |

### `run_metadata.json`

Written only with `--output all`. Records:

- `tool`, `version`, `timestamp_utc`
- `output_kind` — the value of `--output`
- `annotation_source` — path to GTF or index used
- `index_format_version`
- `coordinate_conventions` — restated explicitly
- `query_counts` — `{ total, mapped, unmapped, no_domain_mode }`
- `cli` — full argv of the invocation

---

## Worked example

Input BED:

```
ENSP00000306245    5    100    RD1    TF2
```

`ENSP00000306245` is a `+` strand transcript with this layout:

```
[CDS slice 1]   intron 1   [CDS slice 2]   intron 2   [CDS slice 3]   ...
75278983..94    75279124..  75279877..035   ...
```

Domain `RD1` (aa 5..100) spans the end of slice 1 and most of slice 2. The relevant rows of `isoform_structure.tsv` look like (selected columns):

| feature_type | feature_id | feature_part | feature_genomic_start | feature_genomic_end | aa_start_encoded | aa_end_encoded | overlaps_domain | plot_group |
|---|---|---|---|---|---|---|---|---|
| CDS | CDS_1 | 1 | 75278983 | 75278994 | 1 | 4 | no | CDS_no_domain |
| CDS | CDS_2 | 1 | 75278995 | 75279123 | 5 | 47 | coding_overlap | CDS_domain |
| intron | intron_1 | 1 | 75279124 | 75279876 | NA | NA | inside_domain_genomic_span | intron_domain_span |
| CDS | CDS_3 | 1 | 75279877 | 75280035 | 48 | 100 | coding_overlap | CDS_domain |
| CDS | CDS_4 | 1 | 75280036 | 75280128 | 101 | 131 | no | CDS_no_domain |
| intron | intron_2 | 1 | 75280129 | 75280559 | NA | NA | no | intron |
| ... | ... | ... | ... | ... | ... | ... | ... | ... |
| three_prime_UTR | three_prime_UTR_1 | 1 | 75281422 | 75282230 | NA | NA | no | three_prime_UTR |

To plot, fill by `plot_group`. To highlight the domain, the rows you want are `plot_group ∈ {CDS_domain, intron_domain_span}`. To get just the genomic bases that code the domain, take `domain_cds_segments.bed` (it already contains only the `coding_overlap` rows).

If the same input file had `ENSP00000306245` on its own line (no aa columns), an additional set of rows for the whole transcript would be emitted, all with `overlaps_domain = NA` and `plot_group = CDS / intron / five_prime_UTR / three_prime_UTR`.

---

## MANE Select / Ensembl Canonical

GENCODE GTFs (since v34) carry `tag "MANE_Select"` and `tag "Ensembl_canonical"` on every feature of the relevant transcripts. The tool parses these and surfaces them in two columns on `domain_mapping_summary.tsv` and every feature TSV:

| Column | Values |
|---|---|
| `is_mane_select` | `true` / `false` / `NA` |
| `is_ensembl_canonical` | `true` / `false` / `NA` |

`NA` is reserved for the case where the GTF carries **no** `tag` attribute anywhere — typical of base Ensembl GTFs. We treat absence-of-information differently from "tag attributes exist but this transcript doesn't carry MANE_Select" (which is `false`). Filter to MANE Select queries with:

```bash
awk -F'\t' 'NR==1 || $20=="true"' domain_mapping_summary.tsv   # col 20 = is_mane_select
```

## CDS-length mismatch (Sec, readthrough, incomplete CDS)

The tool maps every query, but flags the rare cases where the CDS isn't a multiple of 3.

| Column | Meaning |
|---|---|
| `cds_length_mismatch` | `true` if `sum(CDS_nt) % 3 != 0` |
| `cds_nt_remainder` | `0`, `1`, or `2` (the remainder) |
| `status` | gains a `_cds_mismatch` suffix (`ok_cds_mismatch`, `partial_cds_mismatch`) |

This catches:

- **Selenoproteins (GPX4, SELENO\*)** — internal `TGA` re-coded as Sec, plus an SECIS-driven C-terminal extension that makes the annotated CDS longer than `protein_length × 3`.
- **Stop-codon readthrough** — annotated CDS extends past a `TGA` / `TAG` that's bypassed by the ribosome.
- **Incomplete 5′ or 3′** GENCODE transcripts with `cds_start_NF` / `cds_end_NF`.

The aa↔nt math still uses ceiling division (`aa = ⌈cds_nt / 3⌉`), so a domain at the very C-terminus of an incomplete CDS may be clipped by 1–2 aa; the `cds_length_mismatch` flag is your hint to check.

## Codons split across exons (1+2 vs 2+1)

A codon that straddles an exon boundary is mapped correctly in both phases — the underlying math is purely cumulative-nt-based, not exon-by-codon. Two synthetic test cases included in `tests/` exercise this:

```
#                            aa 1  aa 2  aa 3
# 1+2 split: CDS_1 = 1 nt of codon 2, CDS_2 = remaining 2 nt
CDS_1   pos 100..103   (4 nt → aa 1 + first base of aa 2)
intron  pos 104..199
CDS_2   pos 200..204   (5 nt → last 2 bases of aa 2 + aa 3)
domain aa 2..2 ⇒ CDS_1 row with aa overlap 2..2 (overlap_fraction_of_domain ≈ 1/3)
                ⇒ CDS_2 row with aa overlap 2..2 (overlap_fraction_of_domain ≈ 2/3)
                ⇒ intron in between marked inside_domain_genomic_span
```

The 2+1 case is symmetric (`CDS_1 = 5 nt`, `CDS_2 = 4 nt`). In both phases, summing `domain_overlap_fraction_of_domain` over `coding_overlap` rows equals `1.0`.
