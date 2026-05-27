# Python API

The wrapper is a thin layer over the C++ binary. It writes inputs to a temp dir, shells out, and reads the TSVs back as pandas DataFrames. The C++ binary remains the source of truth for every mapping decision — Python only assembles BED rows, runs the binary, and parses the outputs.

## `Mapper`

```python
import prot2exon as p2e

mapper = p2e.Mapper(
    index="human.idx",
    binary=None,           # auto-discovered via PROT2EXON_BIN / repo / $PATH
    threads=None,          # default: all cores
    batch_size=0,          # 0 = unbounded (one-shot); >0 enables streaming
    extra_args=(),
)
```

## `Mapper.map(...)` — single query

```python
result = mapper.map(
    "ENSP00000269305",
    aa_start=102, aa_end=292, domain_id="TP53_DBD",
)
result.summary       # one-row DataFrame
result.isoform       # plot-ready DataFrame
result.bed12         # IGV-ready DataFrame
```

ENST also works — it resolves to the same intervals as the matching ENSP, so the two give identical mapping output:

```python
mapper.map("ENST00000269305", aa_start=102, aa_end=292, domain_id="TP53_DBD")
```

## `Mapper.map_batch(...)` — many queries in one call

```python
result = mapper.map_batch([
    {"protein_id": "ENSP00000269305", "aa_start": 102, "aa_end": 292, "domain_id": "TP53_DBD"},
    {"protein_id": "ENSP00000269305", "aa_start": 10,  "aa_end": 50,  "domain_id": "TP53_AD1"},
])
```

For very large query sets (millions), bound RAM by setting `batch_size`:

```python
mapper = p2e.Mapper(index="human.idx", batch_size=10000, threads=8)
result = mapper.map_batch(million_queries)
```

See [[Performance and RAM]] for the 1 M-query benchmark.

## `map_query(...)` — one-off (creates a `Mapper` internally)

```python
result = p2e.map_query(
    "ENSP00000269305", 102, 292, "TP53_DBD",
    index="human.idx",
)
```

Use `Mapper(...).map_batch(...)` if you have many queries — each `map_query` call reloads the index from disk.

## Persisting outputs (`keep_outputs=`)

By default `Mapper.map` / `map_batch` use a temp dir that's deleted when the `MappingResult` goes out of scope. Pass `keep_outputs="path/"` to keep the TSVs and BEDs on disk:

```python
result = mapper.map("ENSP00000269305", 102, 292, "TP53_DBD",
                    keep_outputs="my_results/")
# my_results/{domain_mapping_summary.tsv, isoform_structure.tsv, ...}
```

Same files the CLI would write under `--out-dir`.

## Reading an existing run

```python
result = p2e.read_results_dir("my_results/")
result.summary       # pandas DataFrame
result.isoform       # ...
```

`read_results_dir` parses the output directory of a previous run into a fresh `MappingResult` — useful for picking up where a long batch left off, or analysing CLI runs from Python.

## Binary discovery

`Mapper(index=...)` finds the C++ binary by checking, in order:

1. The `binary=` constructor argument.
2. `$PROT2EXON_BIN`.
3. `<repo>/build/prot2exon` (development checkouts).
4. `<repo>/bin/prot2exon` (the shell wrapper).
5. `prot2exon-core` on `$PATH`, then `prot2exon` on `$PATH`.

For installed users, ship the compiled binary on `$PATH` (pip / conda installs do this automatically) or set `PROT2EXON_BIN` explicitly.

## `MappingResult`

Dataclass with these attributes (all DataFrames):

| Attribute | Source TSV |
|---|---|
| `summary` | `domain_mapping_summary.tsv` |
| `isoform` | `isoform_structure.tsv` |
| `cds_segments` | `domain_cds_segments.tsv` |
| `introns` | `domain_introns.tsv` |
| `cds_bed`, `introns_bed`, `span_bed`, `bed12` | the four `*.bed` files |
| `unmapped` | `unmapped_domains.tsv` (empty if all rows mapped) |
| `metadata` | parsed `run_metadata.json` (only when `--output all`) |
| `out_dir` | path to the temp dir while alive |

Properties: `n_total`, `n_mapped`, `n_unmapped`.

Persist outputs to a permanent location:

```python
result.write("results/my_run/")
```

Filter to one input_id:

```python
sub = result.by_input_id("TP53_DBD")
```

## `plot(...)`

```python
p2e.plot(
    source,                # a MappingResult, an isoform DataFrame, or a path to isoform_structure.tsv
    input_id="TP53_DBD",
    out="tp53_dbd.pdf",    # matplotlib
    html=None,             # plotly HTML
    html_interactive=None,     # interactive HTML
    show_introns=True,
    show_utr=True,
    compact_genomic=False, # clamp introns to 80 bp
    spliced=False,         # drop introns entirely (mutex with compact_genomic)
    link_template=None,    # URL template with {protein_id}, {gene_name}, {transcript_id}, {chrom}
    title=None,
)
```

To render every `input_id` in the source, pass `all=True` instead of `input_id`, or use the explicit helper:

```python
p2e.plot_all(result, out="all_queries.pdf")    # multipage PDF (matplotlib)
p2e.plot_all(result, html="all_queries.html")  # one plotly file per input_id
```

`source` can be a `MappingResult`, an isoform DataFrame, or a path to `isoform_structure.tsv` — all three are accepted:

```python
p2e.plot(result, input_id="TP53_DBD", out="tp53.pdf")
p2e.plot(result.isoform, input_id="TP53_DBD", out="tp53.pdf")
p2e.plot("results/isoform_structure.tsv", input_id="TP53_DBD", out="tp53.pdf")
```

## Interactive viewer helpers

```python
# Standalone HTML file (any browser, offline)
p2e.render_interactive_html(segs, "out.html",
                         link_template="https://...",
                         plot_height=80)

# Inline in a Jupyter notebook
p2e.render_interactive_jupyter(segs,
                            height=None,     # default: auto-resize via postMessage
                            plot_height=140, # main-track height in px
                            link_template="https://...")
```

`segs` is a `list[Segment]`. Build it from a `MappingResult`:

```python
from prot2exon.plot import _segments_from_dataframe
segs_by_id = _segments_from_dataframe(result.isoform)
segs = segs_by_id["TP53_DBD"]
```

Or load straight from disk:

```python
from prot2exon.plot import load_isoform_tsv
segs_by_id = load_isoform_tsv("results/isoform_structure.tsv")
```

See [[Plotting and viewers]] for renderer-specific gestures and flags.
