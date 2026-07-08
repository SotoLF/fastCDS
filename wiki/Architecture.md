# Architecture

## Repo layout

```
fastCDS/
├── src/                       C++ binary
│   ├── main.cpp              CLI entry point
│   ├── gtf_parser.{cpp,hpp}  GTF -> in-memory index (+ index serialise/deserialise)
│   ├── domain_mapper.{cpp,hpp}  per-query mapping; emits per-feature rows
│   ├── output_writer.{cpp,hpp}  TSV/BED writers + StreamingWriter
│   └── utils.{cpp,hpp}       string utilities, attribute helpers
├── include/common.hpp        shared types (Result, Segment, ...)
├── python/fastCDS/         Python wrapper
│   ├── _client.py            Mapper (subprocess wrapper)
│   ├── _result.py            MappingResult dataclass
│   ├── plot.py               matplotlib + plotly renderer
│   ├── _interactive_html.py      standalone interactive HTML viewer (+ Jupyter wrapper)
│   └── fetch.py              `fastCDS fetch` subcommand
├── parsing/                  input adapters: prepare_from_{interpro,uniprot,pfam}, append_custom_proteins
├── tutorial/                 worked examples
│   ├── walkthrough_end_to_end.ipynb   zero-to-figure API tour
│   ├── examples/             small fixtures (tp53_isoforms.tsv, ...)
│   └── reproduce_paper/      reproduce the paper end to end
│       ├── benchmarks/       1M-query + scaling + validation harness
│       └── end_to_end/       atlas / ClinVar / AlphaFold / validation notebooks
├── software_tests/           pytest suite (golden-file + integration)
└── wiki/                     these pages
```

## Pipeline

```
GTF ──▶ gtf_parser ──▶ index (.idx)
                              │
        BED queries ──▶ domain_mapper ──▶ vector<Result>
                                                  │
                                       output_writer ──▶ TSVs + BEDs + run_metadata.json
```

### `gtf_parser`

- Streams the GTF line-by-line; no full DOM in memory.
- Resolves `gene_name` (GENCODE/Ensembl) with a fallback to `gene` (NCBI RefSeq).
- Word-boundary attribute matching to avoid false-positive tag detection (e.g. `locus_tag "..."` no longer trips a `tag` substring check).
- Extracts MANE Select, Ensembl_canonical, and CCDS tags when present.
- Serialises the parsed index to a single binary file (~25x smaller than the source GTF).

### `domain_mapper`

- Loads the index from disk (mmap-friendly, ~1.5 s for human).
- Per query, looks up the transcript, walks its CDS+intron+UTR features, and emits one `Result` per feature plus the global summary row.
- Annotates each row with whether it overlaps the domain (`coding_overlap`, `inside_domain_genomic_span`, `outside`, `NA`).
- **CDS splitting**: when a domain partially covers a CDS exon, the exon is emitted as two rows that share `feature_id`/`exon_number` but differ in `feature_part`. This keeps the original CDS number stable while letting the plotter colour the two pieces differently.
- Supports streaming via `process_domains_streaming(batch_size, callback)` so the writer can drain results in chunks (see `--batch-size`).

### `output_writer`

- Per-result append helpers (`append_summary_row`, `append_isoform_rows`, ...) feed both the one-shot path and the streaming path.
- `StreamingWriter` flushes each batch immediately and frees its memory before the next batch starts.
- `write_all()` wraps the per-file write loops in `#pragma omp parallel sections` so the seven outputs are written concurrently when `--threads > 1`.

## Python wrapper

- `Mapper` shells out to the binary per call. It writes the BED to a temp dir, invokes the binary with the requested `--output` flags, and reads the TSVs back as DataFrames.
- `MappingResult` owns the temp dir while alive - when the result goes out of scope, the dir is cleaned up (unless you call `result.write(...)` to persist).
- The plotter (`plot.py`) takes either a `MappingResult`, a DataFrame, or a path; it derives `list[Segment]` per `input_id` and dispatches to matplotlib, plotly, or `_interactive_html.render_interactive_html`.

## Three renderers, one input

All three renderers consume the same `Segment` list (from `_segments_from_dataframe` or `load_isoform_tsv`). `_out_kind(out)` classifies the `--out` extension and `_render_one_target` dispatches:

- `.pdf` / `.png` / `.svg` -> matplotlib (`_draw_genomic`, `_draw_compact_genomic`, or `_draw_spliced`)
- `.html` + `engine="js"` (default) -> `render_interactive_html` (vanilla JS template in `_interactive_html.py`)
- `.html` + `engine="plotly"` -> plotly (`render_html` with `Bar` + `Scattergl` + rangeslider)

The JS template renders a shared-axis view (compact-mode collapses introns to 80 virtual bp) with a vCRE-style minimap and box-zoom on the main plot. Both file output (`render_interactive_html`) and Jupyter embedding (`render_interactive_jupyter`) share the same `_render_to_string` helper.

## Testing

The suite is pytest, split by concern under `software_tests/`:

- **`test_correctness.py`** - golden-file diffs of the mapper's TSV/BED outputs (`software_tests/golden/`), regenerated with `pytest --update-goldens`.
- **`test_errors.py` / `test_schema.py`** - `status`/`reason` values and the coordinate conventions (BED 0-based vs TSV 1-based).
- **`test_compat.py`** - no-tag GTFs, RefSeq dialect, custom-protein injection (`parsing/append_custom_proteins.py`).
- **`test_bed12.py`** - BED12 block geometry, the `--bed12` add-on, `--batch-size` equivalence.
- **`test_plotting.py`** - plot flags and the plotly/interactive/Jupyter renderers.
- **`test_fetch.py`** - `fastCDS fetch` (offline Zenodo-index download, sha256 verify, `fetch list`).
- **`test_python_api.py`** - the Python wrapper API (`Mapper`, `map_batch`, DataFrame helpers).

`conftest.py` regenerates the synthetic GTFs (`software_tests/make_synthetic_gtf.py`), builds the shared indexes, and maps a fixed query set once per session (the `out_all` fixture). Run with:

```bash
pip install -r software_tests/requirements-dev.txt && pip install -e python/
pytest -q
```
