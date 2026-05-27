# Quickstart

End-to-end in five commands. Assumes you've built the binary per [[Installation]].

## 1. Fetch a stock index

```bash
prot2exon fetch human --release 49 --cache-dir ~/.cache/prot2exon
# prints the resulting .idx path on stdout
```

Same thing from Python (smoother for notebook workflows):

```python
import prot2exon as p2e
idx = p2e.fetch_index("human", release="49")     # Path to .idx
```

Or build manually from any GTF:

```bash
./build/prot2exon \
    --gtf gencode.v49.primary_assembly.annotation.gtf \
    --build-index --index human.idx
```

See [[Genome onboarding]] for non-Ensembl annotations.

## 2. Prepare a BED of domain queries

Plain text, whitespace-separated. Lines starting with `#` are ignored.

```
# rows with a domain (ENSP or ENST, with or without version)
ENSP00000269305   10    50   AD1     TP53_AD1
ENSP00000269305  102   292   DBD     TP53_DBD

# row without a domain (whole-transcript structure only)
ENSP00000418960
```

Columns:

| # | Required | Meaning |
|---|---|---|
| 1 | yes | `ENSP*` or `ENST*` (versioned or not — the suffix is stripped on both sides) |
| 2 | no  | 1-based inclusive aa start (omit for no-domain mode) |
| 3 | no  | 1-based inclusive aa end |
| 4 | no  | `domain_id` (used as `input_id` for tracking) |

See [[Input format]] for prep helpers (UniProt features, InterProScan, HMMER/Pfam → BED).

Or build the queries directly in Python from any of those sources — the same parsing logic the CLI scripts use is also exposed as a clean function-returning-DataFrame API:

```python
queries = p2e.prepare.from_pfam("pfam_hits.dom", mode="scan", id_type="ensp")
# or
queries = p2e.prepare.from_interproscan("interpro.tsv", analyses={"Pfam"})
# or
queries = p2e.prepare.from_uniprot_features("uniprot.dat")
```

Each returns a DataFrame with `protein_id, aa_start, aa_end, domain_id, description` columns that `Mapper.map_batch()` consumes directly.

## 3. Map the BED

```bash
./build/prot2exon \
    --index human.idx \
    --bed queries.bed \
    --out-dir results \
    --output all \
    --threads $(nproc)
```

Outputs (under `results/`):

- `domain_mapping_summary.tsv` — one row per input, ok or not
- `domain_cds_segments.tsv` + `.bed` — CDS exons; subset BED = exons that touch the domain
- `domain_introns.tsv` + `.bed` — intron rows; subset BED = introns inside the domain span
- `domain_span_with_introns.bed` — one row per domain (genomic envelope)
- `isoform_structure.tsv` — full UTR/CDS/intron layout (plot-ready)
- `domain_blocks.bed12` — one IGV-ready row per domain
- `run_metadata.json` — counts and timings
- `unmapped_domains.tsv` — only if at least one row failed

See [[Output modes]] for the per-mode question each table answers.

## 4. Plot a single domain

```bash
prot2exon plot \
    --isoform results/isoform_structure.tsv \
    --input-id TP53_DBD \
    --out tp53_dbd.pdf
```

Or the interactive HTML viewer:

```bash
prot2exon plot \
    --isoform results/isoform_structure.tsv \
    --input-id TP53_DBD \
    --html-interactive tp53_dbd.html \
    --link-template 'https://www.ensembl.org/Homo_sapiens/Transcript/ProteinSummary?p={protein_id}'
```

See [[Plotting and viewers]] for all the rendering options and gestures.

## 5. Or use the Python API

```python
import prot2exon as p2e

mapper = p2e.Mapper(index="human.idx")
result = mapper.map_batch([
    {"protein_id": "ENSP00000269305", "aa_start": 102, "aa_end": 292, "domain_id": "TP53_DBD"},
])
print(result.summary)

# Plot a PDF (matplotlib)
p2e.plot(result, input_id="TP53_DBD", out="tp53_dbd.pdf")

# Or embed the interactive viewer in a Jupyter cell
p2e.render_interactive_jupyter(
    p2e.plot._segments_from_dataframe(result.isoform)["TP53_DBD"],
    plot_height=160,
)
```

Full Python API on [[Python API]].
