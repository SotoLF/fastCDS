# Changelog

All notable changes to fastCDS. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses [semver](https://semver.org/).

## [Unreleased]

### Added
- **`pip install fastCDS` now ships the compiled binary.** The build moved to scikit-build-core, which runs CMake during the wheel build and bundles the binary as `fastCDS/_bin/fastCDS-core`. cibuildwheel publishes pre-built `py3-none` wheels for Linux (manylinux x86_64) and macOS (Intel + Apple Silicon) on each version tag, so all four commands work straight from pip with no compiler. A new cross-platform `fastCDS` console-script (Python dispatcher, mirrors `bin/fastCDS`) is the entry point; binary discovery is centralised in `fastCDS._binary.find_binary` and checks the bundled `_bin/` first. The CMake default no longer hard-codes `-march=native` (it would break portable wheels) - opt in with `-DFASTCDS_NATIVE=ON` for a tuned local build.
- **Subcommand CLI.** The binary now takes one of four commands: `fastCDS index --gtf FILE --out FILE` (build an index; `--index` accepted as an alias for `--out`) and `fastCDS map (--index|--gtf) --bed FILE --out-dir DIR [...]` (map queries), joining the existing `fetch` and `plot` commands. The old bare-flag forms (`--build-index`, `--index ... --bed ...`) still work as an undocumented fallback so existing scripts don't break.
- `fastCDS map --bed12` - emit `domain_blocks.bed` alongside any `--output` mode (no-op under `--output all`/`bed12`), so you can get the IGV-ready BED12 without switching to `--output bed12`.
- **`fastCDS fetch <target>` downloads a pre-built index from Zenodo** (`fetch human` / `mouse` / `yeast` / `human-v86`) with a single sha256-verified download - no GTF parse + build. `fetch` is now Zenodo-only: for any other release, species, or a custom GTF, build the index yourself with `fastCDS index` (a one-time ~15 s step). The earlier `fetch index --preset/--url` sub-subcommand and the `--release` / `--gtf-url` / `ensembl` build paths have been removed to keep `fetch` to a single, predictable job.
- Per-thread-aware 1 MiB write buffer on every TSV/BED writer (`pubsetbuf`) - cuts single-thread wall time on the 1 M-query benchmark by ~17 % (107 s vs 130 s).
- Multi-isoform stack viewer: `render_interactive_html_stack()` and `render_interactive_jupyter_stack()`. All isoforms render on a single shared axis built from the union of their features, so skipped exons appear as empty space lined up across rows.
- `fastCDS.fetch_index()` Python API mirroring the CLI; smoother for notebook workflows.
- `fastCDS.prepare` submodule with `from_pfam()`, `from_interproscan()`, `from_uniprot_features()` returning DataFrames that `Mapper.map_batch()` consumes directly. The standalone scripts under `parsing/` are now thin CLI wrappers around this submodule, so the Python API and CLI never drift.
- New worked-example notebooks: `validation.ipynb`, `software_comparison.ipynb`, `scaling_and_ram.ipynb` (visualise the data behind the Performance and Benchmarking wiki page).
- New **interactive HTML viewer** (`fastCDS.render_interactive_html`) - self-contained vanilla JS, no CDN, with a minimap plus box-zoom, drag-pan and wheel-zoom; also embeddable in Jupyter via `render_interactive_jupyter()`.
- `--compact-genomic` plot mode (matplotlib) that clamps long introns to a fixed display width while keeping CDS/UTR true-scale.
- `--link-template URL` plot flag for clickable external links in HTML output.
- `tutorial/examples/tp53_isoforms.tsv` fixture (four pre-mapped TP53 isoforms) used by the walkthrough notebook.
- Zenodo release bundle (`zenodo_release/` - 510 MB, 4 pre-built indexes + benchmark BEDs + result tables + sha256 manifest + README).
- **Scaling-benchmark harness under [`reproduce_paper/benchmarks/matched/`](reproduce_paper/benchmarks/matched/)**: `run_matched_scaling.py` sweeps every tool under a `wait4()` wall-clock + peak-RSS harness and assembles the table the notebook plots, `prepare_matched_inputs.py` rebuilds the shared query ladder straight from the GTF, and one small runner per tool covers the ensembldb, GenomicFeatures, geneplot and REST routes. Every path is an argument, so the whole benchmark reproduces from a clean checkout.
- `build_human_tool_inputs.py` builds geneplot's human inputs (`h86.ipr`, `ensp_enst.tsv`) from the Ensembl-86 GTF and the Pfam-on-v86 table, so geneplot can be timed on the same human set as everything else. On `Homo_sapiens.GRCh38.86.chr.gtf` it reproduces the published `.ipr` byte for byte.
- `run_ensembl_rest.py` gains `--workers` (paced thread pool, aggregate rate still capped by `--qps`) and `--envelope` (one row per query collapsed to the genomic envelope, for lining REST up against envelope-only tools).

### Changed
- **`fastCDS plot` output is now driven by a single `--out` extension.** `.pdf`/`.png`/`.svg` render a static matplotlib figure and `.html` renders the interactive viewer; the viewer's engine is chosen with `--engine {js,plotly}` (default `js`, the self-contained vanilla-JS renderer). This replaces the separate `--out` / `--html` / `--html-interactive` flags (and the Python `out=` / `html=` / `html_interactive=` kwargs), which have been **removed** - pass `--out x.html [--engine plotly]` / `out="x.html", engine="plotly"` instead. This mirrors the paper's three output types (BED12 track, static figure, interactive viewer), with `js` and `plotly` as two engines of the one viewer rather than separate outputs. An unsupported `--out` extension now fails fast with an actionable error. plotly is now a core dependency, so both viewer engines work from a single `pip install fastCDS` with no extras.
- **Renamed the project to `fastCDS`** (previously `Prot2Exon`, originally `protein2genomic`) - GitHub repo, Docker label, bioconda recipe, PyPI distribution, Python package/module (`import fastCDS`), CLI command, and local clone directory all updated. The Python docs' import alias is now `fc` (was `p2e`).
- README slimmed to ~95 lines (install one-liners + four-command quickstart in CLI and Python + validation/benchmarks + summarized notebooks); dropped the logo/figure images, the wiki-topics table, and the status section. The full reference lives in the wiki.
- **Speed and memory are reported from one run throughout.** README, the wiki and Supplementary Table S2 now quote the same matched benchmark: 5,886 q/s at N = 10,000 single-threaded, ~980x ensembldb, with Ensembl REST quoted at N = 1,000 (the largest N it finishes) and its coordinate agreement at 99.06% over the full 5,000-query set.
- Benchmark flags follow the project name: `--p2e-index` -> `--fastcds-index`, `--p2e-reps` -> `--fastcds-reps`, `--p2e-cds` -> `--fastcds-cds`, `--p2e-index-size` -> `--fastcds-index-size`. The conda env in `environment.yml` is now `fastCDS`.
- The notebooks' default data directory is now `~/Desktop/fastCDS_data`, overridable with `FASTCDS_DATA`.
- geneplot is timed through `benchmarks/matched/run_geneplot.py` under the shared harness, replacing the standalone `geneplot_human.py`.
- **Wiki reorganized into a sequential workflow**: Installation -> Building an index -> Mapping -> Plotting, then Tutorials and Notebooks, Performance and Benchmarking, and the Python API / Architecture / FAQ reference. Merged the old Genome-onboarding + Custom-proteins into *Building an index*, Output-modes + Input-format into *Mapping*, Performance-and-RAM + Validation + Benchmarks into *Performance and Benchmarking*; removed the standalone Quickstart and CLI-reference pages.

### Fixed
- The BED12 output file is now named `domain_blocks.bed` (was `domain_blocks.bed12`). BED12 is a 12-column BED, so IGV/UCSC and other genome browsers expect the standard `.bed` extension; the `.bed12` suffix stopped some tools from loading it. The `--bed12` flag and `--output bed12` mode names are unchanged (they name the format, not the file).
- Destruction-order bug in the buffered `ofstream` (the 1 MiB buffer member was being freed before the base `~ofstream` flushed through it, writing garbage to disk).
- GTF parser: `gene_name` falls back to `gene` for NCBI RefSeq GTFs (which use the latter).
- GTF parser: word-boundary attribute matching - `locus_tag "..."` no longer trips a `tag` substring check.
- `fastCDS fetch` no longer prints "index already exists / pass --force to rebuild" (looked alarming); says "using cached index: ..." instead.
- interactive viewer: pfam atlas figure now displays inline in `nbconvert --execute` runs (forced `%matplotlib inline` backend; `display(fig)` before `savefig(bbox_inches="tight")`).

## [2.2.0] - prior to changelog

Initial public-facing release. Baseline functionality:

- C++17 binary (`fastCDS`) with GTF parsing, binary index, OpenMP parallel mapping, 7+ output files (4 TSVs, 3 BEDs, BED12, metadata JSON).
- Python wrapper (`pip install -e python/`) with `Mapper` / `MappingResult` / `plot()`.
- Phase 1: 10-test golden-file suite (now 118 tests).
- Phase 2: 100.00 % exact match against ensembldb's v86 on 5,000 stratified queries.
- Phase 3: head-to-head vs ensembldb and Ensembl REST - 5,886 q/s single-threaded.
- Phase 4: Pfam-A proteome atlas + ClinVar enrichment notebooks.
- Phase 6: `prepare_from_interpro.py`, `prepare_from_pfam.py`, `prepare_from_uniprot_features.py`.

[Unreleased]: https://github.com/SotoLF/fastCDS/compare/v2.2.0...HEAD
[2.2.0]: https://github.com/SotoLF/fastCDS/releases/tag/v2.2.0
