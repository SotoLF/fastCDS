# prot2exon

Map protein-domain amino-acid coordinates to their underlying genomic CDS/UTR/intron structure, using any Ensembl, GENCODE, or NCBI RefSeq GTF.

For each input query (a `protein_id` **or** a `transcript_id`, optionally with an aa range), prot2exon answers two related but distinct questions:

1. **Mapping** — *which exact genomic bases code this domain?*
2. **Structure / visualization** — *how is the whole transcript organised in 5′UTR / CDS / 3′UTR / intron, and where does the domain fall on it?*

A C++17 binary does the heavy lifting (≤ 1 µs per query on a warm index), a Python wrapper gives you DataFrames, and a `plot` subcommand renders matplotlib / plotly / a vanilla-JS standalone HTML viewer.

## Wiki navigation

| Page | What's on it |
|---|---|
| [[Installation]] | Build the C++ binary, install the Python wrapper, Docker. |
| [[Quickstart]] | End-to-end in five commands. |
| [[Input format]] | BED-like input, ENSP vs ENST, no-domain mode, prep scripts. |
| [[Output modes]] | `coding` / `introns` / `span` / `isoform` / `bed12` / `all`. |
| [[Plotting and viewers]] | matplotlib, plotly, the interactive HTML viewer, Jupyter embed. |
| [[Genome onboarding]] | `prot2exon fetch`, manual recipes, GTF compatibility. |
| [[Custom proteins]] | Append unannotated proteins (transgenes, non-reference ORFs, …) to an existing GTF. |
| [[Python API]] | `Mapper`, `MappingResult`, `plot`, `render_interactive_jupyter`. |
| [[CLI reference]] | Every CLI flag at a glance: map, plot, fetch. |
| [[Notebooks]] | The three worked-example notebooks shipped under `notebooks/` — walkthrough, Pfam atlas, ClinVar enrichment. |
| [[Validation]] | 100 % exact match vs ensembldb on 5,000 stratified queries — design, results, gotchas. |
| [[Benchmarks]] | 4-tool head-to-head (prot2exon, ensembldb, TransVar, Ensembl REST), scaling, parallel scaling. |
| [[Performance and RAM]] | `--threads`, `--batch-size`, 1 M-query benchmark. |
| [[Architecture]] | What lives where: parser, mapper, writer, index format. |
| [[FAQ]] | Common gotchas: CDS-length mismatch, MANE Select, "index already exists", ... |

## Project links

- Repository: <https://github.com/SotoLF/Prot2Exon>
- Owner: Luis F. Soto Ugaldi ([@SotoLF](https://github.com/SotoLF))
- Collaborator: George D. Muñoz Esquivel ([@george123ya](https://github.com/george123ya))
- Issue tracker: open a ticket on the same repo.
- Citation: see the bottom of [[FAQ]].
