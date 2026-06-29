# Tutorials and notebooks

This page walks you from a clean checkout to a rendered figure, then points you at the worked example notebooks that live under [`tutorial/`](https://github.com/SotoLF/fastCDS/tree/main/notebooks) in the repo. Start with the quickstart, then dive into whichever notebook matches your use case.

## Quickstart: a real run you can copy-paste

From zero to a figure in about five lines: fetch a human index, write a tiny BED with one TP53 domain, map it onto genomic exon structure, and plot. The CLI exposes exactly four commands — `fetch`, `index`, `map`, `plot` — see [[Index]], [[Mapping]], and [[Plotting]] for the full reference.

```bash
# 1. Get a human index (pre-built, from Zenodo -> ~/.cache/fastCDS/human.idx)
fastCDS fetch human

# 2. Write a tiny query: the TP53 DNA-binding domain
printf 'ENSP00000269305\t102\t292\tTP53_DBD\n' > queries.bed

# 3. Map the domain onto its genomic exon structure
fastCDS map --index ~/.cache/fastCDS/human.idx --bed queries.bed --out-dir results --output all

# 4. Plot it
fastCDS plot --isoform results/isoform_structure.tsv --input-id TP53_DBD --out tp53.pdf
```

The equivalent in Python:

```python
import fastCDS as fc

idx = fc.fetch_index("human")
mapper = fc.Mapper(index=idx)
result = mapper.map_batch([
    {"protein_id": "ENSP00000269305", "aa_start": 102, "aa_end": 292, "domain_id": "TP53_DBD"},
])
fc.plot(result, input_id="TP53_DBD", out="tp53.pdf")
```

Both paths produce the same TP53 DBD figure. The index is downloaded and cached on first use; later runs reuse it.

## End-to-end walkthrough notebook

The headline tour takes you from an empty working directory to a domain mapped onto its genomic exon structure and rendered as a figure — the whole `fetch → map → plot` surface in one place, every output shown inline:

1. **Get an index** — `fc.fetch_index("yeast")` pulls a small index, so the notebook runs in seconds with nothing else installed.
2. **Build queries** — compose a query BED by hand, or generate one from a domain database with `fc.prepare.from_pfam` / `from_interproscan` / `from_uniprot_features`.
3. **Map** — `mapper.map_batch(...)` returns the mapping as pandas DataFrames (`summary`, `isoform`, …).
4. **Plot** — render a bundled TP53 fixture every way fastCDS offers: static PDF, compact-genomic layout, the standalone interactive HTML viewer, an inline Jupyter embed, and a Plotly comparison.

It's the guided, executed version of the [Quickstart](#quickstart-a-real-run-you-can-copy-paste) above — a good place to start if you want to see each step's output before adapting it to your own data.

Notebook: [`walkthrough_end_to_end.ipynb`](https://github.com/SotoLF/fastCDS/blob/main/tutorial/walkthrough_end_to_end.ipynb) · [**view on nbviewer**](https://nbviewer.org/github/SotoLF/fastCDS/blob/main/tutorial/walkthrough_end_to_end.ipynb) (the inline interactive HTML viewers render on nbviewer, not on GitHub's notebook view)

## Domain functional atlas (and ClinVar)

This notebook maps every Pfam-A domain on the human proteome (Ensembl release 115) onto genomic exon structure and, counting one canonical isoform per gene, reports: the single-exon fraction **and** domain size by function (Fig 1D), where each domain sits along its protein (start / middle / end landmarks, in amino-acid space so UTRs play no part), and how fully a domain fills the exons that encode it (completeness) — each also broken down per Pfam family. It then runs the **ClinVar** test on the same atlas: whether pathogenic missense variants are enriched in CDS exons that code a Pfam-A domain (Fisher's exact, Panel C).

Notebook: [`domain_functional_atlas.ipynb`](https://github.com/SotoLF/fastCDS/blob/main/tutorial/reproduce_paper/end_to_end/domain_functional_atlas.ipynb)

## AlphaFold pLDDT at splice junctions

Per-residue AlphaFold pLDDT against distance (in nucleotides) to the nearest exon-exon junction, across the canonical human proteome — a 2D density showing model confidence dipping right at the splice site.

Notebook: [`alphafold_plddt_junctions.ipynb`](https://github.com/SotoLF/fastCDS/blob/main/tutorial/reproduce_paper/end_to_end/alphafold_plddt_junctions.ipynb)
