# Plotting (`fastCDS plot`)

Once you have an `isoform_structure.tsv` from [[Mapping]], `fastCDS plot` renders it as a static figure (matplotlib) or as an interactive viewer (plotly or a self-contained vanilla-JS viewer). The plotter reads the TSV directly and never re-derives coordinates from the genome, so anything you can express by editing the table is plottable.

## Static figures

### CLI

Render a single query to a vector PDF:

```bash
fastCDS plot \
    --isoform results/isoform_structure.tsv \
    --input-id TP53_DBD \
    --out tp53.pdf
```

The output format follows the file extension (`.pdf`, `.png`, or `.svg`). To render every `input_id` in the TSV at once, use `--all`; with a `.pdf` target this becomes a single multipage PDF (one query per page), while a `.png`/`.svg` target writes one file per query named `base.<input_id>.ext`:

```bash
fastCDS plot --isoform results/isoform_structure.tsv --all --out queries.pdf   # multipage PDF
fastCDS plot --isoform results/isoform_structure.tsv --all --out queries.png   # queries.TP53_DBD.png, ...
```

### Python

From a `MappingResult` (or an isoform DataFrame, or a path to the TSV):

```python
import fastCDS as fc

result = fc.map_query("ENSP00000269305", 10, 50, "TP53_DBD", index="human.idx")
fc.plot(result, input_id="TP53_DBD", out="tp53.pdf")
```

`fc.plot()` returns the matplotlib `Figure` (useful when you pass neither `out` nor an HTML target). For batches, `fc.plot_all(source, out="queries.pdf")` mirrors `--all`. See [[Python API]] for the full client.

### Arguments

| Flag / kwarg | Default | Effect |
|---|---|---|
| `--isoform FILE` | required | Path to the `isoform_structure.tsv` to plot. |
| `--input-id ID` | — | Render a single query (mutually exclusive with `--all`). |
| `--all` | — | Render every `input_id`; multipage PDF if `--out` ends in `.pdf`, else one file per query. |
| `--out FILE` | — | Static matplotlib output; format from the `.pdf`/`.png`/`.svg` extension. |
| `--title STR` | derived | Override the auto-generated title. |
| `--width`, `--height` | 12, 2.6 | Figure size in inches. |
| `--no-highlight` | off | Don't color CDS-in-domain segments red. |
| `--no-introns` | off | Hide intron lines. |
| `--no-utr` | off | Hide UTR boxes. |
| `--spliced` | off | Concatenate non-intron features in translation order (no introns). |
| `--compact-genomic` | off | Genomic order, but clamp each intron to a fixed display width (best for long-intron genes); CDS/UTR stay at true bp scale with a `//` compression mark. Mutually exclusive with `--spliced`. |

The Python `plot()` / `plot_all()` functions accept the same toggles as keyword arguments (`input_id`, `out`, `title`, `width`, `height`, `show_introns`, `show_utr`, `highlight_domain`, `spliced`, `compact_genomic`).

## Interactive viewers

### CLI

Two interactive paths share the same command; the renderer is chosen by the output flag:

```bash
# plotly HTML (CDN-backed; hover tooltips and a bottom rangeslider)
fastCDS plot --isoform results/isoform_structure.tsv --input-id TP53_DBD \
    --html tp53.html

# self-contained vanilla-JS viewer (no CDN, single offline file)
fastCDS plot --isoform results/isoform_structure.tsv --input-id TP53_DBD \
    --html-interactive tp53.html \
    --link-template 'https://www.ensembl.org/Homo_sapiens/Transcript/ProteinSummary?p={protein_id}'
```

`--html` needs plotly (`pip install plotly`); `--html-interactive` has no extra dependencies. You can pass both alongside `--out` in one run. The standalone viewer supports box-zoom (drag to zoom into a genomic range), shift-drag to pan, mouse-wheel zoom, double-click to reset, a draggable minimap, UTR rendering with strand arrows, and a Compact / True-genomic layout toggle. `--link-template` adds a clickable linkout next to the title using the placeholders `{protein_id}`, `{gene_name}`, `{transcript_id}`, `{chrom}`, `{start}`, `{end}`.

With `--all`, each interactive flag writes one file per query (`base.<input_id>.ext`).

### Python

For a single isoform, obtain the segments from a result and embed the viewer inline in a notebook:

```python
import fastCDS as fc
from fastCDS.plot import _segments_from_dataframe

result = fc.map_query("ENSP00000269305", 10, 50, "TP53_DBD", index="human.idx")
segments = _segments_from_dataframe(result.isoform)["TP53_DBD"]

fc.render_interactive_jupyter(segments, plot_height=160)
```

For several isoforms stacked in one viewer, pass the whole `dict[input_id -> segments]` to the stack variant:

```python
segs_by_id = _segments_from_dataframe(result.isoform)
fc.render_interactive_jupyter_stack(segs_by_id, plot_height=40)
```

To write standalone HTML files instead of embedding inline, use the file builders `fc.render_interactive_html(segments, "tp53.html")` and `fc.render_interactive_html_stack(segs_by_id, "stack.html")`. See [[Tutorials and Notebooks]] for end-to-end notebook examples.

### Arguments

| Flag / kwarg | Default | Effect |
|---|---|---|
| `--html FILE` | — | plotly interactive HTML (CDN-backed); requires plotly. |
| `--html-interactive FILE` | — | Self-contained vanilla-JS viewer (no CDN, offline-capable). |
| `--link-template URL` / `link_template=` | — | External linkout next to the title (HTML viewers only); placeholders `{protein_id}`, `{gene_name}`, `{transcript_id}`, `{chrom}`, `{start}`, `{end}`. |
| `--height N` (CLI) | 2.6 | matplotlib figure height in inches (static path). |
| `plot_height=` (Python) | 140 single / 40 stack | Main-track height in pixels for the Jupyter / standalone viewers. |
| `height=` (Python) | auto | Pin the Jupyter iframe height in px (for static exports where the auto-resize handshake can't fire). |

The hidden `--html-tfregdb2` CLI flag and the `render_tfregdb2_*` Python functions are deprecated aliases for the `--html-interactive` / `render_interactive_*` names, kept only for backwards compatibility and slated for removal.
