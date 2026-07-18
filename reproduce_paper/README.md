# reproduce_paper/ - reproduce the paper's example analyses

These notebooks reproduce the manuscript's example application end to end. Every
heavy input is fetched by the notebook itself (nothing needs to be staged by
hand), and each step is **cache-aware**, so a re-run reuses what is already on
disk.

For the **speed / accuracy benchmarks** (speed vs `ensembldb`,
`GenomicFeatures`, geneplot and Ensembl REST; coordinate agreement additionally
vs TransVar) see
[`benchmarks/README.md`](benchmarks/README.md), which has its own one-block
reproduction harness and covers `software_comparison.ipynb` and
`scaling_and_ram.ipynb`.

## Notebooks

| Notebook | Backs (manuscript) |
|---|---|
| `notebooks/isoform_domain_conservation_analysis.ipynb` | Projects every source-isoform Pfam domain onto the gene's alternative protein-coding isoforms and scores retention (intact / partially trimmed / skipped). The analysis behind **Figure 1C-F** (the final panels are composed from this output). |
| `notebooks/scaling_and_ram.ipynb` | Runtime scaling, throughput, and the thread x batch-size sweep (**Figure 1B**, **Figure S2**). |
| `notebooks/software_comparison.ipynb` | Per-tool, per-category coordinate agreement (**Table S1**) and mapping speed / peak memory (**Table S2**). |

`../tutorial/walkthrough_end_to_end.ipynb` is a standalone zero-to-figure tour of
the API and does not depend on these.

## Data directory and outputs

All notebooks read and write a single data directory,
`~/Desktop/fastCDS_data/` (set once at the top of each notebook as
`DATA` - edit that line, or set `FASTCDS_DATA`, to relocate). The first run
downloads the inputs there;
later runs and the other notebooks reuse them. Figures are written as vector PDFs
into [`figures/`](figures/), named for the manuscript panel they back (for
example `Figure_S2A_scaling.pdf`, `Figure_S2B_throughput.pdf`).

## Environment

```bash
pip install fastCDS pandas numpy scipy matplotlib jupyter
# the fastCDS C++ binary ships with the wheel; no separate build needed
```

(The benchmark notebooks additionally need an R/Bioconductor environment, plus
TransVar for the agreement table - see
[`benchmarks/README.md`](benchmarks/README.md).)

## Data sources (all downloaded automatically)

| Data | Source | Fetched by |
|---|---|---|
| Ensembl release-115 GTF -> binary index | `ftp.ensembl.org/pub/release-115` | conservation |
| Pfam-A domains per protein (ENSP, aa range) | Ensembl **BioMart** (release 115) | conservation |

## Counting convention

The conservation analysis compares each Pfam domain on a gene's **reference
isoform** (MANE Select / Ensembl canonical) against the gene's other
protein-coding isoforms, and reports retention and reading-frame preservation
(`intact` >= 0.80, `partial` 0.50-0.80, `lost` < 0.50 of the domain's coding
bases). The convention is stated at the top of the notebook.
