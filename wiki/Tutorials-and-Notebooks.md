# Tutorials and notebooks

This page walks you from a clean checkout to a rendered figure, then points you at the worked example notebooks that live under [`tutorial/`](https://github.com/SotoLF/Prot2Exon/tree/main/notebooks) in the repo. Start with the quickstart, then dive into whichever notebook matches your use case.

## Quickstart: a real run you can copy-paste

From zero to a figure in about five lines: fetch a human index, write a tiny BED with one TP53 domain, map it onto genomic exon structure, and plot. The CLI exposes exactly four commands — `fetch`, `index`, `map`, `plot` — see [[Index]], [[Mapping]], and [[Plotting]] for the full reference.

```bash
# 1. Get a human index (pre-built, from Zenodo -> ~/.cache/prot2exon/human.idx)
prot2exon fetch human

# 2. Write a tiny query: the TP53 DNA-binding domain
printf 'ENSP00000269305\t102\t292\tTP53_DBD\n' > queries.bed

# 3. Map the domain onto its genomic exon structure
prot2exon map --index ~/.cache/prot2exon/human.idx --bed queries.bed --out-dir results --output all

# 4. Plot it
prot2exon plot --isoform results/isoform_structure.tsv --input-id TP53_DBD --out tp53.pdf
```

The equivalent in Python:

```python
import prot2exon as p2e

idx = p2e.fetch_index("human")
mapper = p2e.Mapper(index=idx)
result = mapper.map_batch([
    {"protein_id": "ENSP00000269305", "aa_start": 102, "aa_end": 292, "domain_id": "TP53_DBD"},
])
p2e.plot(result, input_id="TP53_DBD", out="tp53.pdf")
```

Both paths produce the same TP53 DBD figure. The index is downloaded and cached on first use; later runs reuse it.

## End-to-end walkthrough notebook

The headline tour takes you from an empty working directory to a domain mapped onto its genomic exon structure and rendered as a figure — the whole `fetch → map → plot` surface in one place, every output shown inline:

1. **Get an index** — `p2e.fetch_index("yeast")` pulls a small index, so the notebook runs in seconds with nothing else installed.
2. **Build queries** — compose a query BED by hand, or generate one from a domain database with `p2e.prepare.from_pfam` / `from_interproscan` / `from_uniprot_features`.
3. **Map** — `mapper.map_batch(...)` returns the mapping as pandas DataFrames (`summary`, `isoform`, …).
4. **Plot** — render a bundled TP53 fixture every way prot2exon offers: static PDF, compact-genomic layout, the standalone interactive HTML viewer, an inline Jupyter embed, and a Plotly comparison.

It's the guided, executed version of the [Quickstart](#quickstart-a-real-run-you-can-copy-paste) above — a good place to start if you want to see each step's output before adapting it to your own data.

Notebook: [`walkthrough_end_to_end.ipynb`](https://github.com/SotoLF/Prot2Exon/blob/main/tutorial/walkthrough_end_to_end.ipynb) · [**view on nbviewer**](https://nbviewer.org/github/SotoLF/Prot2Exon/blob/main/tutorial/walkthrough_end_to_end.ipynb) (the inline interactive HTML viewers render on nbviewer, not on GitHub's notebook view)

## Pfam proteome atlas

This notebook maps every Pfam-A domain across the human proteome (~150 K domains over ~19 K proteins) onto genomic exon structure and computes architecture statistics: the share of domains encoded by exactly one CDS exon vs two or more, the distribution of coding exons touched per domain, the median and maximum intronic span within each domain envelope, and the distribution of `fraction_domain_in_largest_exon`. The result is the proteome-wide atlas figure (Panel B of the accompanying paper's Figure 1).

Notebook: [`pfam_proteome_atlas.ipynb`](https://github.com/SotoLF/Prot2Exon/blob/main/tutorial/reproduce_paper/end_to_end/pfam_proteome_atlas.ipynb)

## ClinVar enrichment

This notebook tests whether pathogenic missense variants are enriched in CDS exons that code for Pfam-A domains, relative to benign controls. It downloads the latest ClinVar VCF (GRCh38), filters to missense SNVs with a clear clinical significance, picks the top genes by pathogenic count, asks for each variant whether it falls in a CDS exon coding at least one Pfam-A domain, and runs a Fisher's exact test on the enrichment. The result is Panel C of the figure.

Notebook: [`clinvar_pathogenic.ipynb`](https://github.com/SotoLF/Prot2Exon/blob/main/tutorial/reproduce_paper/end_to_end/clinvar_pathogenic.ipynb)
