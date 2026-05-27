# Architecture

## Repo layout

```
Prot2Exon/
├── src/                       C++ binary
│   ├── main.cpp              CLI entry point
│   ├── gtf_parser.{cpp,hpp}  GTF → in-memory index (+ index serialise/deserialise)
│   ├── domain_mapper.{cpp,hpp}  per-query mapping; emits per-feature rows
│   ├── output_writer.{cpp,hpp}  TSV/BED writers + StreamingWriter
│   └── utils.{cpp,hpp}       string utilities, attribute helpers
├── include/common.hpp        shared types (Result, Segment, ...)
├── python/prot2exon/         Python wrapper
│   ├── _client.py            Mapper (subprocess wrapper)
│   ├── _result.py            MappingResult dataclass
│   ├── plot.py               matplotlib + plotly renderer
│   ├── _interactive_html.py      standalone interactive HTML viewer (+ Jupyter wrapper)
│   └── fetch.py              `prot2exon fetch` subcommand
├── scripts/                  helpers: prepare_from_interpro/uniprot/pfam, append_custom_proteins
├── notebooks/                walkthrough + paper figure notebooks
├── tests/run_tests.py        end-to-end test suite (~109 assertions)
├── benchmarks/               1M-query + scaling harness
├── examples/                 small fixtures (tp53_isoforms.tsv, ...)
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
- Serialises the parsed index to a single binary file (~25× smaller than the source GTF).

### `domain_mapper`

- Loads the index from disk (mmap-friendly, ~1.5 s for human).
- Per query, looks up the transcript, walks its CDS+intron+UTR features, and emits one `Result` per feature plus the global summary row.
- Annotates each row with whether it overlaps the domain (`coding_overlap`, `inside_domain_genomic_span`, `outside`, `NA`).
- **CDS splitting**: when a domain partially covers a CDS exon, the exon is emitted as two rows that share `feature_id`/`exon_number` but differ in `feature_part`. This keeps the original CDS number stable while letting the plotter colour the two pieces differently.
- Supports streaming via `process_domains_streaming(batch_size, callback)` so the writer can drain results in chunks (see `--batch-size`).

### `output_writer`

- Per-result append helpers (`append_summary_row`, `append_isoform_rows`, …) feed both the one-shot path and the streaming path.
- `StreamingWriter` flushes each batch immediately and frees its memory before the next batch starts.
- `write_all()` wraps the per-file write loops in `#pragma omp parallel sections` so the seven outputs are written concurrently when `--threads > 1`.

## Python wrapper

- `Mapper` shells out to the binary per call. It writes the BED to a temp dir, invokes the binary with the requested `--output` flags, and reads the TSVs back as DataFrames.
- `MappingResult` owns the temp dir while alive — when the result goes out of scope, the dir is cleaned up (unless you call `result.write(...)` to persist).
- The plotter (`plot.py`) takes either a `MappingResult`, a DataFrame, or a path; it derives `list[Segment]` per `input_id` and dispatches to matplotlib, plotly, or `_interactive_html.render_interactive_html`.

## Three renderers, one input

All three plot paths consume the same `Segment` list (from `_segments_from_dataframe` or `load_isoform_tsv`). The renderer is chosen by which output flag you pass:

- `--out` → matplotlib (`_draw_genomic`, `_draw_compact_genomic`, or `_draw_spliced`)
- `--html` → plotly (`render_html` with `Bar` + `Scattergl` + rangeslider)
- `--html-interactive` → `render_interactive_html` (vanilla JS template in `_interactive_html.py`)

The JS template renders a shared-axis view (compact-mode collapses introns to 80 virtual bp) with a vCRE-style minimap and box-zoom on the main plot. Both file output (`render_interactive_html`) and Jupyter embedding (`render_interactive_jupyter`) share the same `_render_to_string` helper.

## Testing

`tests/run_tests.py` is a single-file end-to-end suite that:

1. Generates synthetic GTFs (`scripts/make_synthetic_gtf.py`).
2. Builds an index.
3. Runs the binary against handcrafted BED queries.
4. Asserts on every output (~109 assertions covering coordinate conventions, CDS splitting, MANE Select, tag extraction, RefSeq dialect, custom-protein injection, plotly/interactive/jupyter renderers, `--batch-size` equivalence, …).

Run with:

```bash
python3 tests/run_tests.py
```

Expect `109 passed, 2 failed` — the two failures require matplotlib + pandas in the system python and aren't regressions.
