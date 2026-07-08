# tutorial/reproduce_paper/ — reproduce the paper's example analyses

These notebooks reproduce the **example application** in the manuscript
(Section 4: the human Pfam atlas, functional-architecture analyses, ClinVar
enrichment, and the AlphaFold pLDDT-junction analysis) **end to end, from data
download to figure**. Every heavy input is fetched by the notebook itself
(nothing needs to be staged by hand), and each step is **cache-aware**, so a
re-run reuses what is already on disk.

For the **speed / accuracy benchmarks** (vs `ensembldb`, `GenomicFeatures`,
TransVar, Ensembl REST, VisProDom, geneplot) see
[`benchmarks/README.md`](benchmarks/README.md), which has its own one-block
reproduction harness and covers the `software_comparison.ipynb` and
`scaling_and_ram.ipynb` notebooks.

## What lands where

All notebooks read and write a single data directory,
`~/Desktop/protein2genomic_data/` (set once at the top of each notebook as
`DATA` — edit that line to relocate). The first run downloads the inputs there;
later runs and the other notebooks reuse them.

## Environment

```bash
pip install fastCDS pandas numpy scipy matplotlib jupyter
# the fastCDS C++ binary ships with the wheel; no separate build needed
```

(The benchmark notebooks additionally need an R/Bioconductor + TransVar
environment — see [`benchmarks/README.md`](benchmarks/README.md).)

## Data sources (all downloaded automatically)

| Data | Source | Fetched by |
|---|---|---|
| Ensembl release-115 GTF -> binary index | `ftp.ensembl.org/pub/release-115` | `domain_functional_atlas`, `alphafold_plddt_junctions` |
| Pfam-A domains per protein (ENSP, aa range) | Ensembl **BioMart** (release 115) | `domain_functional_atlas` |
| Pfam -> GO functional classes | `current.geneontology.org` `external2go/pfam2go` | `domain_functional_atlas` |
| ClinVar missense variants | `ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38` | `domain_functional_atlas` |
| ENSP -> UniProt accession | Ensembl **BioMart** (release 115) | `alphafold_plddt_junctions` |
| Per-residue pLDDT (AlphaFold DB v6) | `alphafold.ebi.ac.uk` | `alphafold_plddt_junctions` |

## Run order

Run these in order — the first builds the shared release-115 index and the
mapped Pfam atlas that the second reuses:

| # | Notebook | Produces (manuscript) | Approx. wall time* |
|---|---|---|---|
| 1 | **`end_to_end/domain_functional_atlas.ipynb`** | Builds the r115 index + maps the Pfam atlas, then domain fragmentation and completeness across isoforms, the source charts for the Figure 1 domain-retention panels. | ~5 min (first run: + GTF download/index ~3 min, BioMart pull ~6 min) |
| 2 | `end_to_end/alphafold_plddt_junctions.ipynb` | AlphaFold per-residue pLDDT vs distance to the nearest exon-exon junction, across the whole canonical human proteome (Fig 1G). | ~8 min (first run: + ~18 k AlphaFold JSON fetches ~6 min) |

\* on a desktop with a warm network; first-run download times added in italics.

`../walkthrough_end_to_end.ipynb` is a standalone zero-to-figure tour of the API
and does not depend on the atlas.

## Where the figures land

Each notebook writes its figure as a vector PDF into
[`figures_pdf/`](figures_pdf/), named for the manuscript panel it backs (for
example `Figure_1G_plddt_junctions.pdf`, `Figure_S2A_scaling.pdf`,
`Figure_S1C_isoform_viewer.pdf`).

## Counting convention

The functional-architecture analyses count **one representative (Ensembl
canonical) isoform per gene**, so a domain on a multi-isoform gene is not counted
once per splice variant. The ClinVar interval set instead uses the **union of
domain-coding exons across all isoforms** (a variant counts if it falls in any
isoform's domain-coding exon). Both choices are stated in the relevant notebook.
