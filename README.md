<p align="center">
  <img alt="fastCDS" width="300"
       src="https://raw.githubusercontent.com/SotoLF/fastCDS/main/wiki/images/fastcds_icon.png">
</p>

# fastCDS

<p align="center">
  <a href="https://pypi.org/project/fastCDS/"><img alt="PyPI"
      src="https://img.shields.io/pypi/v/fastCDS?logo=pypi&logoColor=white&color=3775A9"></a>
  <a href="https://bioconda.github.io/recipes/fastCDS/README.html"><img alt="Bioconda"
      src="https://img.shields.io/conda/vn/bioconda/fastCDS?label=bioconda&logo=anaconda&logoColor=white&color=3C7E22"></a>
  <a href="https://pixi.sh/"><img alt="pixi"
      src="https://img.shields.io/badge/pixi-add%20fastCDS-yellow?logo=python&logoColor=white"></a>
  <a href="https://github.com/SotoLF/fastCDS/wiki"><img alt="Docs"
      src="https://img.shields.io/badge/docs-wiki-blueviolet?logo=gitbook&logoColor=white"></a>
  <a href="LICENSE"><img alt="License"
      src="https://img.shields.io/github/license/SotoLF/fastCDS?color=success"></a>
</p>

<p align="center">
  <img alt="fastCDS interactive isoform viewer"
       src="https://raw.githubusercontent.com/SotoLF/fastCDS/main/wiki/images/plot_isoform_stack.gif"
       width="840">
  <br>
  <em>The TP53 DNA-binding domain mapped across four isoforms on a shared genomic axis (self-contained interactive viewer).</em>
</p>

Map protein-domain amino-acid coordinates to their underlying genomic CDS / UTR / intron structure, using any Ensembl, GENCODE, or NCBI RefSeq GTF.

For each input query (a `protein_id` **or** a `transcript_id`, optionally with an aa range), fastCDS answers two related but distinct questions:

1. **Mapping** - *which exact genomic bases code this domain?*
2. **Structure** - *how is the whole transcript organised in 5'UTR / CDS / 3'UTR / intron, and where does the domain fall on it?*

A C++17 binary does the heavy lifting: once the index is loaded, the mapping itself is ~1 us per query, so end-to-end throughput is dominated by output formatting (~5,800 queries/s for a single-isoform TSV, ~2,800 q/s writing the full `--output all` set). fastCDS produces three kinds of output: a **BED12 track** for genome browsers, a **static figure** (matplotlib PDF/PNG/SVG), and an **interactive HTML viewer**. A Python wrapper adds pandas DataFrames; for the viewer you pick the engine (self-contained vanilla JS by default, or plotly).

📖 **[Full documentation lives in the wiki](https://github.com/SotoLF/fastCDS/wiki)**

## Install

```bash
pip install fastCDS                                 # pre-built wheel: binary + wrapper (Linux/macOS)
mamba install -c bioconda -c conda-forge fastCDS    # or via conda/mamba
pixi add -c bioconda -c conda-forge fastCDS         # or via pixi
```

The pip wheel bundles the compiled binary, so all four commands work immediately. On platforms without a wheel (e.g. Windows) pip builds from source - see [Installation](https://github.com/SotoLF/fastCDS/wiki/Installation) for the toolchain and a from-source build.

## Quickstart

The workflow is three steps - get an index (`index` from a GTF, or `fetch` a pre-built one), `map` your domain queries onto it, then `plot`:

```mermaid
flowchart LR
    GTF[GTF annotation] --> INDEX([fastCDS index])
    ZEN[Zenodo] --> FETCH([fastCDS fetch])
    INDEX --> IDX[index<br/>human.idx]
    FETCH --> IDX
    IDX --> MAP([fastCDS map])
    BED[query BED<br/>protein + aa range] --> MAP
    MAP --> TSV[isoform_structure.tsv]
    MAP --> B12[BED tracks<br/>BED12 + coding / introns / span]
    TSV --> PLOT([fastCDS plot])
    PLOT --> STATIC[static figure<br/>.pdf / .png / .svg]
    PLOT --> INTER[interactive viewer<br/>.html: js or plotly]
    B12 --> IGV[genome browser<br/>IGV / UCSC]

    classDef cmd fill:#2f6db0,color:#ffffff,stroke:#1c4a7d,stroke-width:1px;
    classDef file fill:#eef2ff,color:#111111,stroke:#9aa7d0,stroke-width:1px;
    class INDEX,FETCH,MAP,PLOT cmd;
    class GTF,ZEN,IDX,BED,TSV,B12,STATIC,INTER,IGV file;
```

```bash
# 1. Get an index (one-time per annotation)
fastCDS index --gtf my.gtf --out human.idx   # build from a GTF you already have
# or skip the GTF - grab a pre-built index from Zenodo:
fastCDS fetch human --out human.idx

# 2. Map a BED of domain queries
fastCDS map \
    --index human.idx \
    --bed queries.bed --out-dir results --output all --threads 8

# 3. Plot a single domain - the --out extension picks the format
fastCDS plot \
    --isoform results/isoform_structure.tsv \
    --input-id TP53_DBD \
    --out tp53_dbd.pdf              # static figure (.pdf/.png/.svg)
#   --out tp53_dbd.html            # interactive viewer (vanilla JS by default)
#   --out tp53_dbd.html --engine plotly
```

Run `fastCDS fetch list` to see every target and source.

The same workflow from Python:

```python
import fastCDS as fc

# 1. Get an index (one-time)
idx = fc.build_index("my.gtf", out="human.idx")    # build from a GTF; Path-returning
# idx = fc.fetch_index("human")                    # or grab a pre-built one from Zenodo

# 2. Map queries
mapper = fc.Mapper(index=str(idx))
result = mapper.map_batch([
    {"protein_id": "ENSP00000269305", "aa_start": 102, "aa_end": 292, "domain_id": "TP53_DBD"},
])
result.summary       # one-row DataFrame
result.isoform       # plot-ready DataFrame

# 3. Plot - the out extension picks the format; engine picks the HTML renderer
fc.plot(result, input_id="TP53_DBD", out="tp53_dbd.pdf")                    # static figure
fc.plot(result, input_id="TP53_DBD", out="tp53_dbd.html")                   # interactive (vanilla JS)
fc.plot(result, input_id="TP53_DBD", out="tp53_dbd.html", engine="plotly")  # interactive (plotly)
```

Full reference: [Building an index](https://github.com/SotoLF/fastCDS/wiki/Index), [Mapping](https://github.com/SotoLF/fastCDS/wiki/Mapping), [Plotting](https://github.com/SotoLF/fastCDS/wiki/Plotting), [Python API](https://github.com/SotoLF/fastCDS/wiki/Python-API).

## Validation + benchmarks

- **100.00 % exact match vs ensembldb** on 5,000 stratified queries (9 strata covering single/multi-exon, both strands, codon-split, selenoproteins, incomplete CDS) - zero off-by-ones, zero structural mismatches.
- **~970x faster than ensembldb** (also ~4.4x TransVar, ~5,400x Ensembl REST) - measured *end-to-end*: total wall time from process start to all results written, including the one-time index load, over N = 10,000 queries on one thread (fastCDS 5,847 q/s vs ensembldb 6 q/s). The gap grows with N because fastCDS amortizes its ~1.2 s index load; the per-query mapping itself, once loaded, is faster still.
- Full design, numbers, and scaling curves: [Performance and benchmarking](https://github.com/SotoLF/fastCDS/wiki/Performance-and-Benchmarking). Reproduce via [`tutorial/reproduce_paper/benchmarks/README.md`](tutorial/reproduce_paper/benchmarks/README.md).

## Notebooks

Worked examples under [`tutorial/`](tutorial/) - each opens in [Colab](https://colab.research.google.com/github/SotoLF/fastCDS) or [nbviewer](https://nbviewer.org/github/SotoLF/fastCDS). They reproduce the manuscript's example analyses **end to end, from data download to figure**; see [`tutorial/reproduce_paper/README.md`](tutorial/reproduce_paper/README.md) for the data sources, run order, and runtimes.

| Notebook | What it covers |
|---|---|
| [`walkthrough_end_to_end.ipynb`](tutorial/walkthrough_end_to_end.ipynb) ([view on nbviewer](https://nbviewer.org/github/SotoLF/fastCDS/blob/main/tutorial/walkthrough_end_to_end.ipynb)) | Zero-to-figure tour: `fetch_index` -> BED prep -> `map_batch` -> all plot styles. The interactive HTML viewers only render on **nbviewer**, not GitHub. |
| [`domain_functional_atlas.ipynb`](tutorial/reproduce_paper/end_to_end/domain_functional_atlas.ipynb) | The Ensembl-r115 Pfam atlas: single-exon fraction + domain size by function (Fig 1D), domain position along the protein (Q1/Q50/Q100), completeness, all per-family, **and the ClinVar pathogenic-variant enrichment** in domain-coding exons. |
| [`alphafold_plddt_junctions.ipynb`](tutorial/reproduce_paper/end_to_end/alphafold_plddt_junctions.ipynb) | AlphaFold per-residue pLDDT vs distance to the nearest exon-exon junction (2D density), across the canonical human proteome. |
| [`software_comparison.ipynb`](tutorial/reproduce_paper/end_to_end/software_comparison.ipynb) | **Speed + accuracy** vs `ensembldb` / GenomicFeatures / TransVar / Ensembl REST / VisProDom / geneplot: agreement + end-to-end throughput + RSS, all on the same human set. |
| [`scaling_and_ram.ipynb`](tutorial/reproduce_paper/end_to_end/scaling_and_ram.ipynb) | fastCDS measured against itself (not other tools): wall-clock + peak-RSS scaling curves, OpenMP speedup, and the `--batch-size` RAM cap at N = 1 M. |

## Citation

A formal citation will land when the accompanying manuscript is posted. Until then, please cite the repository URL: <https://github.com/SotoLF/fastCDS>. See also [`CITATION.cff`](CITATION.cff).

## License

MIT - see [LICENSE](LICENSE).

## Maintainers

- **Owner**: Luis F. Soto Ugaldi ([@SotoLF](https://github.com/SotoLF))
- **Collaborator**: George D. Muñoz Esquivel ([@george123ya](https://github.com/george123ya))

Issues and pull requests welcome.
