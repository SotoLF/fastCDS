# Notebooks

Worked examples live under [`notebooks/`](https://github.com/SotoLF/Prot2Exon/tree/main/notebooks) in the repo.

## `walkthrough_end_to_end.ipynb`

The headline tour — five steps from "I have nothing on disk" to "I have a domain mapped onto its genomic exon structure and rendered as a figure". Uses **yeast** for the live mapping demo (fast, network-light) and the bundled **TP53 fixture** for the plot showcases (richer human exon structure).

Sections:

1. **Setup** — `import prot2exon as p2e` (assumes pip/bioconda/pixi install; one comment for dev checkouts).
2. **`p2e.fetch_index("yeast")`** — one-shot index download + build.
3. **Build queries** — either compose a BED directly, or use `p2e.prepare.from_pfam` / `from_interproscan` / `from_uniprot_features` to convert real-world database outputs into a DataFrame.
4. **`mapper.map_batch(queries)`** — live yeast mapping.
5. **TP53 fixture** — load `examples/tp53_isoforms.tsv` (canonical TP53 isoform pre-mapped against the DBD).
6. **matplotlib PDF** — `p2e.plot(tp53, input_id="TP53_canonical", out=...)`.
7. **Compact-genomic mode** — `p2e.plot(..., compact_genomic=True)`.
8. **Interactive HTML** — `p2e.plot(..., html_interactive=...)`.
9. **Inline embed** — `p2e.render_interactive_jupyter(...)`.
10. **Plotly comparison** — `p2e.plot(..., html=...)`.

Run with `jupyter nbconvert --to notebook --execute --inplace` to refresh outputs.

## `validation.ipynb`

Loads `benchmarks/validate_vs_ensembldb.py`'s output (`table1.tsv` under `validation_v86/` and `validation_v113/`) and renders the per-stratum exact-match figure. These are the numbers behind the [[Validation]] wiki page — keep them in lockstep.

Sections:

1. Matched-annotation path (v86): 100 % on every stratum, headline barplot.
2. Annotation-drift path (v113): same logic, GENCODE v49 vs Ensembl 113 EnsDb. The lower headline number is an EnsDb gap (ENSPs without CDS linkage), not a prot2exon bug.
3. Takeaways.

## `software_comparison.ipynb`

The 4-tool head-to-head from [[Benchmarks]]. Loads `rest_table.tsv`, `transvar_table.tsv`, and the v86 ensembldb table; renders agreement + throughput + peak-RSS bars side-by-side.

Sections:

1. Per-tool agreement vs prot2exon (honest denominators).
2. Throughput @ N = 10,000 (log-scale q/s) + peak-RSS bars.
3. Takeaways.

## `benchmarking.ipynb`

Visualises `benchmarks/scaling_benchmark.py` + `parallel_benchmark.py` outputs.

Sections:

1. Scaling curves — wall-clock and peak-RSS vs N, log-log, prot2exon and ensembldb on the same axes.
2. Parallel speedup at N = 100,000 (with the "ideal linear" reference) — explains why it plateaus at ~1.34× (single-threaded TSV writer).
3. `--batch-size` cap at N = 1M (peak RAM 10.4 GB → 0.96 GB).

## `pfam_proteome_atlas.ipynb`

Map every Pfam-A domain on the human proteome (`EnsDb.Hsapiens.v86`, ~150 K domains across ~19 K proteins) and compute architecture statistics:

- % of domains encoded by exactly 1 CDS exon vs ≥ 2 CDS exons
- Distribution of `n_coding_exons_touched` per domain
- Median / max **intronic span** within the domain envelope
- Distribution of `fraction_domain_in_largest_exon`

This produces **Panel B** of the accompanying paper's Figure 1.

## `clinvar_pathogenic.ipynb`

Test the hypothesis "pathogenic missense variants are enriched in CDS exons that code for Pfam-A domains, relative to benign controls" on ClinVar data:

1. Download the latest ClinVar VCF (GRCh38), filter to missense SNVs with `CLNSIG ∈ {Pathogenic, Likely_pathogenic, Benign, Likely_benign}`.
2. Pick the top N genes by pathogenic count.
3. For each variant, ask: does it fall in a CDS exon that codes at least one Pfam-A domain?
4. Fisher's exact test on pathogenic enrichment in domain-coding exons.

This produces **Panel C** of the figure.

## Regenerating + executing

`notebooks/generate_notebooks.py` is the single source of truth — it overwrites each `.ipynb` file every run, which intentionally wipes any embedded outputs (the .ipynb on disk is treated as a derived artefact, not a manual document).

To keep outputs embedded, pass `--run` so the script also `jupyter nbconvert`-executes each notebook in place:

```bash
# Regenerate AND execute every notebook
python3 notebooks/generate_notebooks.py --run

# Subset
python3 notebooks/generate_notebooks.py --run validation benchmarking

# Bump per-cell timeout (default 600 s) for the slow data notebooks
python3 notebooks/generate_notebooks.py --run --timeout 1800

# Continue past cell errors (useful when iterating on plot code)
python3 notebooks/generate_notebooks.py --run --allow-errors
```

The available names are: `walkthrough_end_to_end`, `validation`, `software_comparison`, `benchmarking`, `pfam_proteome_atlas`, `clinvar_pathogenic`.

## Running the notebooks

The walkthrough is self-contained (yeast index is ~2 MB, the TP53 fixture lives in the repo). The Pfam atlas + ClinVar notebooks need:

- A built human index (`p2e.fetch_index("human", release="49")` or similar).
- An EnsDb-derived Pfam BED (`benchmarks/extract_pfam_domains.py`).
- The ClinVar VCF (downloaded by the notebook on first run).

All three are tested by the repo's CI via `python3 notebooks/generate_notebooks.py` (regenerates the .ipynb files from a single `generate_notebooks.py` source), so refactors don't silently break them.
