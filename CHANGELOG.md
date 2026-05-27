# Changelog

All notable changes to Prot2Exon. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses [semver](https://semver.org/).

## [Unreleased]

### Added
- `prot2exon fetch index --url URL --out PATH [--sha256 HEX]` subcommand — download a pre-built `.idx` from any URL (including the Zenodo deposit) and skip the GTF parse + build step. Supports cache hits, sha256 verification, and `--force` re-download.
- Per-thread-aware 1 MiB write buffer on every TSV/BED writer (`pubsetbuf`) — cuts single-thread wall time on the 1 M-query benchmark by ~17 % (107 s vs 130 s).
- Multi-isoform stack viewer: `render_interactive_html_stack()` and `render_interactive_jupyter_stack()`. All isoforms render on a single shared axis built from the union of their features, so skipped exons appear as empty space lined up across rows.
- `prot2exon.fetch_index()` Python API mirroring the CLI; smoother for notebook workflows.
- `prot2exon.prepare` submodule with `from_pfam()`, `from_interproscan()`, `from_uniprot_features()` returning DataFrames that `Mapper.map_batch()` consumes directly. The standalone scripts under `scripts/` are now thin CLI wrappers around this submodule, so the Python API and CLI never drift.
- New worked-example notebooks: `validation.ipynb`, `software_comparison.ipynb`, `benchmarking.ipynb` (visualise the data behind the Validation + Benchmarks wiki pages).
- New **interactive HTML viewer** (`prot2exon.render_interactive_html`) — self-contained vanilla JS, no CDN, with vCRE-style minimap + box-zoom + drag-pan + wheel-zoom; also embeddable in Jupyter via `render_interactive_jupyter()`. (Formerly named "TFRegDB2 viewer" after the project the design was ported from; the public API kept old names as deprecated aliases for one release.)
- `--compact-genomic` plot mode (matplotlib) that clamps long introns to a fixed display width while keeping CDS/UTR true-scale.
- `--link-template URL` plot flag for clickable external links in HTML output.
- `examples/tp53_isoforms.tsv` fixture (four pre-mapped TP53 isoforms) used by the walkthrough notebook.
- Zenodo release bundle (`zenodo_release/` — 510 MB, 4 pre-built indexes + benchmark BEDs + result tables + sha256 manifest + README).

### Changed
- Renamed from `protein2genomic` to **Prot2Exon** (GitHub repo + Docker label + bioconda recipe + Python wrapper README; local clone directory can stay as-is).
- README slimmed from 733 lines to ~110 with badges + install table + quickstart + wiki link. The full reference now lives across the 13 wiki pages.
- Notebook generator (`notebooks/generate_notebooks.py`) gained `--run` (executes after generating, so outputs persist) and `--out-dir` (test-suite uses a tempdir so the regression check can't wipe embedded outputs).

### Fixed
- Destruction-order bug in the buffered `ofstream` (the 1 MiB buffer member was being freed before the base `~ofstream` flushed through it, writing garbage to disk).
- GTF parser: `gene_name` falls back to `gene` for NCBI RefSeq GTFs (which use the latter).
- GTF parser: word-boundary attribute matching — `locus_tag "..."` no longer trips a `tag` substring check.
- `prot2exon fetch` no longer prints "index already exists / pass --force to rebuild" (looked alarming); says "using cached index: …" instead.
- interactive viewer: pfam atlas figure now displays inline in `nbconvert --execute` runs (forced `%matplotlib inline` backend; `display(fig)` before `savefig(bbox_inches="tight")`).

## [2.2.0] - prior to changelog

Initial public-facing release. Baseline functionality:

- C++17 binary (`prot2exon`) with GTF parsing, binary index, OpenMP parallel mapping, 7+ output files (4 TSVs, 3 BEDs, BED12, metadata JSON).
- Python wrapper (`pip install -e python/`) with `Mapper` / `MappingResult` / `plot()`.
- Phase 1: 10-test golden-file suite (now 118 tests).
- Phase 2: 100.00 % exact match against ensembldb's v86 on 5,000 stratified queries.
- Phase 3: 4-tool head-to-head (ensembldb / TransVar / Ensembl REST) — 5,847 q/s single-threaded.
- Phase 4: Pfam-A proteome atlas + ClinVar enrichment notebooks.
- Phase 6: `prepare_from_interpro.py`, `prepare_from_pfam.py`, `prepare_from_uniprot_features.py`.

[Unreleased]: https://github.com/SotoLF/Prot2Exon/compare/v2.2.0...HEAD
[2.2.0]: https://github.com/SotoLF/Prot2Exon/releases/tag/v2.2.0
