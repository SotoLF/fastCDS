# Prot2Exon

<p align="center">
  <a href="https://pypi.org/project/prot2exon/"><img alt="PyPI"
      src="https://img.shields.io/pypi/v/prot2exon?logo=pypi&logoColor=white&color=3775A9"></a>
  <a href="https://bioconda.github.io/recipes/prot2exon/README.html"><img alt="Bioconda"
      src="https://img.shields.io/conda/vn/bioconda/prot2exon?label=bioconda&logo=anaconda&logoColor=white&color=3C7E22"></a>
  <a href="https://pixi.sh/"><img alt="pixi"
      src="https://img.shields.io/badge/pixi-add%20prot2exon-yellow?logo=python&logoColor=white"></a>
  <a href="https://github.com/SotoLF/Prot2Exon/wiki"><img alt="Docs"
      src="https://img.shields.io/badge/docs-wiki-blueviolet?logo=gitbook&logoColor=white"></a>
  <a href="LICENSE"><img alt="License"
      src="https://img.shields.io/github/license/SotoLF/Prot2Exon?color=success"></a>
</p>

Map protein-domain amino-acid coordinates to their underlying genomic CDS / UTR / intron structure, using any Ensembl, GENCODE, or NCBI RefSeq GTF.

For each input query (a `protein_id` **or** a `transcript_id`, optionally with an aa range), Prot2Exon answers two related but distinct questions:

1. **Mapping** — *which exact genomic bases code this domain?*
2. **Structure** — *how is the whole transcript organised in 5′UTR / CDS / 3′UTR / intron, and where does the domain fall on it?*

A C++17 binary does the heavy lifting: once the index is loaded, the mapping itself is ~1 µs per query, so end-to-end throughput is dominated by output formatting (~5,800 queries/s for a single-isoform TSV, ~2,800 q/s writing the full `--output all` set). A Python wrapper gives you pandas DataFrames and three plot styles (matplotlib, plotly, and a vanilla-JS standalone HTML viewer).

📖 **[Full documentation lives in the wiki →](https://github.com/SotoLF/Prot2Exon/wiki)**

## Install

```bash
pip install prot2exon                                 # pre-built wheel: binary + wrapper (Linux/macOS)
mamba install -c bioconda -c conda-forge prot2exon    # or via conda/mamba
pixi add -c bioconda -c conda-forge prot2exon         # or via pixi
```

The pip wheel bundles the compiled binary, so all four commands work immediately. On platforms without a wheel (e.g. Windows) pip builds from source — see [Installation](https://github.com/SotoLF/Prot2Exon/wiki/Installation) for the toolchain and a from-source build.

## Quickstart

The workflow is four commands — `index` (or `fetch`) to get an index, `map` to project queries, `plot` to render:

```bash
# 1. Get an index (one-time per annotation)
prot2exon index --gtf my.gtf --out human.idx   # build from a GTF you already have
# or skip the GTF — grab a pre-built index from Zenodo:
prot2exon fetch human --out human.idx

# 2. Map a BED of domain queries
prot2exon map \
    --index human.idx \
    --bed queries.bed --out-dir results --output all --threads 8

# 3. Plot a single domain
prot2exon plot \
    --isoform results/isoform_structure.tsv \
    --input-id TP53_DBD \
    --html-interactive tp53_dbd.html
```

Run `prot2exon fetch list` to see every target and source.

The same workflow from Python:

```python
import prot2exon as p2e

# 1. Get an index (one-time)
idx = p2e.build_index("my.gtf", out="human.idx")    # build from a GTF; Path-returning
# idx = p2e.fetch_index("human")                    # or grab a pre-built one from Zenodo

# 2. Map queries
mapper = p2e.Mapper(index=str(idx))
result = mapper.map_batch([
    {"protein_id": "ENSP00000269305", "aa_start": 102, "aa_end": 292, "domain_id": "TP53_DBD"},
])
result.summary       # one-row DataFrame
result.isoform       # plot-ready DataFrame

# 3. Plot — three rendering paths, same data
p2e.plot(result, input_id="TP53_DBD", out="tp53_dbd.pdf")                # matplotlib
p2e.plot(result, input_id="TP53_DBD", html="tp53_dbd.html")              # plotly
p2e.plot(result, input_id="TP53_DBD", html_interactive="tp53_dbd.html")  # vanilla JS
```

Full reference: [Building an index](https://github.com/SotoLF/Prot2Exon/wiki/Index), [Mapping](https://github.com/SotoLF/Prot2Exon/wiki/Mapping), [Plotting](https://github.com/SotoLF/Prot2Exon/wiki/Plotting), [Python API](https://github.com/SotoLF/Prot2Exon/wiki/Python-API).

## Validation + benchmarks

- **100.00 % exact match vs ensembldb** on 5,000 stratified queries (9 strata covering single/multi-exon, both strands, codon-split, selenoproteins, incomplete CDS) — zero off-by-ones, zero structural mismatches.
- **~970× faster than ensembldb** (also ~4.4× TransVar, ~5,400× Ensembl REST) — measured *end-to-end*: total wall time from process start to all results written, including the one-time index load, over N = 10,000 queries on one thread (prot2exon 5,847 q/s vs ensembldb 6 q/s). The gap grows with N because prot2exon amortizes its ~1.2 s index load; the per-query mapping itself, once loaded, is faster still.
- Full design, numbers, and scaling curves: [Performance and benchmarking](https://github.com/SotoLF/Prot2Exon/wiki/Performance-and-Benchmarking). Reproduce via [`tutorial/reproduce_paper/benchmarks/README.md`](tutorial/reproduce_paper/benchmarks/README.md).

## Notebooks

Worked examples under [`tutorial/`](tutorial/) — each opens in [Colab](https://colab.research.google.com/github/SotoLF/Prot2Exon) or [nbviewer](https://nbviewer.org/github/SotoLF/Prot2Exon). They reproduce the manuscript's example analyses **end to end, from data download to figure**; see [`tutorial/reproduce_paper/README.md`](tutorial/reproduce_paper/README.md) for the data sources, run order, and runtimes.

| Notebook | What it covers |
|---|---|
| [`walkthrough_end_to_end.ipynb`](tutorial/walkthrough_end_to_end.ipynb) ([view on nbviewer](https://nbviewer.org/github/SotoLF/Prot2Exon/blob/main/tutorial/walkthrough_end_to_end.ipynb)) | Zero-to-figure tour: `fetch_index` → BED prep → `map_batch` → all plot styles. The interactive HTML viewers only render on **nbviewer**, not GitHub. |
| [`domain_functional_atlas.ipynb`](tutorial/reproduce_paper/end_to_end/domain_functional_atlas.ipynb) | Builds the Ensembl-r115 Pfam atlas and the three functional-architecture analyses: single-exon fraction by domain function (Fig 1D), domain position in the transcript, and completeness. |
| [`pfam_proteome_atlas.ipynb`](tutorial/reproduce_paper/end_to_end/pfam_proteome_atlas.ipynb) | Map every Pfam-A domain on the human proteome; single- vs multi-exon architecture, intron burden, dominant-exon fraction. |
| [`clinvar_pathogenic.ipynb`](tutorial/reproduce_paper/end_to_end/clinvar_pathogenic.ipynb) | Pathogenic-variant enrichment in domain-coding exons (ClinVar missense). |
| [`alphafold_plddt_junctions.ipynb`](tutorial/reproduce_paper/end_to_end/alphafold_plddt_junctions.ipynb) | AlphaFold per-residue pLDDT as a function of distance to exon-exon junctions, across the canonical human proteome. |
| [`validation.ipynb`](tutorial/reproduce_paper/end_to_end/validation.ipynb) | **Accuracy** vs `ensembldb` (and `GenomicFeatures::proteinToGenome`): 100% exact match across the 9-stratum, 5,000-query set, with the per-stratum figure. |
| [`software_comparison.ipynb`](tutorial/reproduce_paper/end_to_end/software_comparison.ipynb) | **Speed** vs `ensembldb` / TransVar / Ensembl REST: agreement + throughput + RSS head-to-head at N = 10,000. |
| [`scaling_and_ram.ipynb`](tutorial/reproduce_paper/end_to_end/scaling_and_ram.ipynb) | prot2exon measured against itself (not other tools): wall-clock + peak-RSS scaling curves, OpenMP speedup, and the `--batch-size` RAM cap at N = 1 M. |

## Citation

A formal citation will land when the accompanying manuscript is posted. Until then, please cite the repository URL: <https://github.com/SotoLF/Prot2Exon>. See also [`CITATION.cff`](CITATION.cff).

## License

MIT — see [LICENSE](LICENSE).

## Maintainers

- **Owner**: Luis F. Soto Ugaldi ([@SotoLF](https://github.com/SotoLF))
- **Collaborator**: George D. Muñoz Esquivel ([@george123ya](https://github.com/george123ya))

Issues and pull requests welcome.
