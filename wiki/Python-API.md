# Python API

The wrapper is a thin layer over the C++ binary. It writes inputs to a temp dir, shells out, and reads the TSVs back as pandas DataFrames. The C++ binary remains the source of truth for every mapping decision — Python only assembles BED rows, runs the binary, and parses the outputs.

## Getting an index

Two functions mirror the `index` / `fetch` commands and return the `Path` to a `.idx` (see [[Index]]):

```python
import fastCDS as fc

idx = fc.build_index("combined.gtf", out="human.idx")  # build from a local GTF (`fastCDS index`)
idx = fc.fetch_index("human")                           # pre-built, from Zenodo (`fastCDS fetch`)
```

`fetch_index` only serves the pre-built Zenodo targets (`"human"`, `"mouse"`, `"mouse-vm25"`, `"human-v86"`, `"human-v95"`, `"human-v115"`, `"yeast"` — run `fastCDS fetch list` to see them). For any other release, species, or a custom GTF, build the index yourself with `build_index`.

`build_index(gtf, out=None, *, binary=None, force=False)` caches: if `out` exists it's returned untouched unless `force=True`.

## `Mapper`

```python
import fastCDS as fc

mapper = fc.Mapper(
    index="human.idx",
    binary=None,           # auto-discovered via FASTCDS_BIN / repo / $PATH
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
mapper = fc.Mapper(index="human.idx", batch_size=10000, threads=8)
result = mapper.map_batch(million_queries)
```

See [[Performance and Benchmarking]] for the 1 M-query benchmark.

## `map_query(...)` — one-off (creates a `Mapper` internally)

```python
result = fc.map_query(
    "ENSP00000269305",
    aa_start=102, aa_end=292, domain_id="TP53_DBD",
    index="human.idx",
)
```

Use `Mapper(...).map_batch(...)` if you have many queries — each `map_query` call reloads the index from disk.

## Persisting outputs (`keep_outputs=`)

By default `Mapper.map` / `map_batch` use a temp dir that's deleted when the `MappingResult` goes out of scope. Pass `keep_outputs="path/"` to keep the TSVs and BEDs on disk:

```python
result = mapper.map("ENSP00000269305",
                    aa_start=102, aa_end=292, domain_id="TP53_DBD",
                    keep_outputs="my_results/")
# my_results/{domain_mapping_summary.tsv, isoform_structure.tsv, ...}
```

Same files the CLI would write under `--out-dir`.

## Reading an existing run

```python
result = fc.read_results_dir("my_results/")
result.summary       # pandas DataFrame
result.isoform       # ...
```

`read_results_dir` parses the output directory of a previous run into a fresh `MappingResult` — useful for picking up where a long batch left off, or analysing CLI runs from Python.

## Binary discovery

`Mapper(index=...)` finds the C++ binary by checking, in order:

1. The `binary=` constructor argument.
2. `$FASTCDS_BIN`.
3. `<repo>/build/fastCDS` (development checkouts).
4. `<repo>/bin/fastCDS` (the shell wrapper).
5. `fastCDS-core` on `$PATH`, then `fastCDS` on `$PATH`.

For installed users, ship the compiled binary on `$PATH` (pip / conda installs do this automatically) or set `FASTCDS_BIN` explicitly.

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

One output argument, `out`, whose extension picks the format: `.pdf`/`.png`/`.svg` render a static matplotlib figure, `.html` renders the interactive viewer (engine chosen by `engine=`).

```python
# Static figure — format follows the extension
fc.plot(result, input_id="TP53_DBD", out="tp53_dbd.pdf")   # or .png / .svg

# Interactive viewer — .html; the default engine is the self-contained vanilla JS
fc.plot(result, input_id="TP53_DBD", out="tp53_dbd.html")
fc.plot(result, input_id="TP53_DBD", out="tp53_dbd.html", engine="plotly")

# No out= → nothing is written; you get the matplotlib Figure back to tweak
fig = fc.plot(result, input_id="TP53_DBD")
```

`plot()` returns the matplotlib `Figure` for static output (and when `out` is `None`), or `None` for `.html` output.

**`source` (first argument)** accepts three things interchangeably — a `MappingResult`, its isoform DataFrame, or a path to an `isoform_structure.tsv` on disk:

```python
fc.plot(result,            input_id="TP53_DBD", out="tp53.pdf")   # a MappingResult
fc.plot(result.isoform,    input_id="TP53_DBD", out="tp53.pdf")   # its isoform DataFrame
fc.plot("run/isoform_structure.tsv", input_id="TP53_DBD", out="tp53.pdf")  # a file on disk
```

**Full keyword list** (all keyword-only after `source`):

```python
fc.plot(
    result,                # MappingResult | isoform DataFrame | path to isoform_structure.tsv
    input_id="TP53_DBD",   # which query to draw (optional if the source has exactly one)
    out="tp53_dbd.pdf",    # extension picks the format: .pdf/.png/.svg static, .html interactive
    engine="js",           # interactive renderer for .html: "js" (default) or "plotly"
    show_introns=True,     # draw intron lines
    show_utr=True,         # draw UTR boxes
    highlight_domain=True, # colour the CDS-in-domain segments red
    compact_genomic=False, # keep genomic order but clamp long introns to a fixed display width
    spliced=False,         # drop introns entirely, concatenate exons (mutex with compact_genomic)
    link_template=None,    # linkout URL; {protein_id},{gene_name},{transcript_id},{chrom},{start},{end}
    title=None,            # override the auto-generated title
    width=12.0, height=2.6,  # static figure size in inches
)
```

To render **every** `input_id` in the source, use `plot_all` — a `.pdf` target becomes one multipage PDF, any other extension writes one file per query (`base.<input_id>.ext`):

```python
fc.plot_all(result, out="all_queries.pdf")                    # multipage PDF (static)
fc.plot_all(result, out="all_queries.png")                    # all_queries.TP53_DBD.png, ...
fc.plot_all(result, out="all_queries.html")                   # one vanilla-JS file per input_id
fc.plot_all(result, out="all_queries.html", engine="plotly")  # one plotly file per input_id
```

## Interactive viewer helpers

`fc.plot(..., out="x.html")` is all most people need. These lower-level helpers exist for two cases the top-level `plot()` doesn't cover: embedding the viewer **inline in a Jupyter notebook**, and controlling the viewer's pixel `plot_height`.

They take a `segs` argument — a `list[Segment]` for one query. Build it first, from a `MappingResult` or straight from a TSV on disk:

```python
from fastCDS.plot import load_isoform_tsv, _segments_from_dataframe

# from a MappingResult already in memory
segs = _segments_from_dataframe(result.isoform)["TP53_DBD"]

# ...or load from disk
segs = load_isoform_tsv("results/isoform_structure.tsv")["TP53_DBD"]
```

Then render it:

```python
# Standalone HTML file (any browser, offline)
fc.render_interactive_html(segs, "out.html",
                         link_template="https://...",
                         plot_height=80)

# Inline in a Jupyter notebook (returns an IPython.display.HTML)
fc.render_interactive_jupyter(segs,
                            height=None,     # default: auto-resize via postMessage
                            plot_height=140, # main-track height in px
                            link_template="https://...")
```

See [[Plotting]] for renderer-specific gestures and flags.
