# Plotting and viewers

prot2exon ships three rendering paths for `isoform_structure.tsv`:

| Renderer | When to use it | File output |
|---|---|---|
| **matplotlib** | Paper figures, batch PDFs, no JS | PDF / PNG / SVG |
| **plotly** | Interactive in-browser, hover tooltips, rangeslider | one self-contained HTML |
| **interactive viewer (vanilla JS)** | Same interactivity as plotly, no CDN, smaller payload, polished domain-aware look | one self-contained HTML |

All three are exposed via the same `prot2exon plot` CLI and the same Python `plot()` call. The renderer is chosen by which output flag you pass.

## CLI

```bash
prot2exon plot \
    --isoform results/isoform_structure.tsv \
    --input-id TP53_DBD \
    --out tp53_dbd.pdf \            # matplotlib
    --html tp53_dbd_plotly.html \   # plotly
    --html-interactive tp53_dbd.html   # the interactive viewer
```

You can pass any subset. `--all` renders every `input_id` in the TSV (matplotlib only: multi-panel PDF; plotly / interactive: one HTML per input_id with `.{id}.html` filename pattern).

## Full flag reference

| Flag | Default | Effect |
|---|---|---|
| `--input-id ID` / `--all` | — | Pick one query, or render every `input_id` in the TSV. |
| `--out FILE` | — | matplotlib PDF / PNG / SVG. With `--all`, PDF becomes multipage. |
| `--html FILE` | — | Interactive plotly HTML. Requires `plotly` (install with `pip install "prot2exon[html]"`). |
| `--html-interactive FILE` | — | Self-contained vanilla-JS HTML viewer (no CDN, no JS deps). |
| `--no-highlight` | off | Don't colour CDS-in-domain segments red. |
| `--no-introns` | off | Hide intron lines (X axis stays genomic — introns become gaps). |
| `--no-utr` | off | Hide UTR boxes. |
| `--spliced` | off | matplotlib only: concatenate non-intron features in translation order (no introns at all). |
| `--compact-genomic` | off | Genomic order, but clamp each intron to a fixed display width. Use for human genes where one long intron dwarfs the exons; CDS / UTR stay at true bp scale and a `//` mark flags compression. Mutually exclusive with `--spliced`. |
| `--link-template URL` | — | (HTML only) Clickable linkout next to the title. Placeholders: `{protein_id}`, `{gene_name}`, `{transcript_id}`, `{chrom}`, `{start}`, `{end}`. |
| `--width`, `--height` | 12, 2.2 | matplotlib figure size in inches. |
| `--title STR` | derived | Override the auto-generated title. |

The plot reads the TSV directly — it doesn't re-derive coordinates from the genome — so any feature you can express by editing `isoform_structure.tsv` (filter to a specific transcript, join with extra annotation, etc.) is plottable.

`--link-template` examples:

- **Ensembl** (recommended for ENSPs): `https://www.ensembl.org/Homo_sapiens/Transcript/ProteinSummary?p={protein_id}`
- **UCSC**: `https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position={chrom}:{start}-{end}`
- **Custom internal LIMS**: `https://lims.example.com/protein/{protein_id}`

> **Don't use the UniProt entry endpoint with ENSP IDs.** `uniprot.org/uniprotkb/{id}/entry` only resolves UniProt accessions; an ENSP returns 404. Either use the UniProt **search** URL (`https://www.uniprot.org/uniprotkb?query={protein_id}`) or convert the column to UniProt accessions before plotting.

## The standalone interactive HTML viewer

Self-contained vanilla JS — no CDN, no JS deps. Single ~35 KB file that works offline and embeds in any host page.

Gestures:

| Gesture | Effect |
|---|---|
| Drag on the main plot | Box-zoom (the area outside the rectangle shades dark; the badge shows the genomic range live) |
| Shift + drag on the main plot | Pan |
| Mouse wheel | Zoom around the cursor |
| Double-click | Reset zoom |
| Drag the indigo minimap rectangle | Pan |
| Drag minimap edges | Zoom (resize the viewport) |
| Click empty minimap | Recenter viewport |
| Layout radio (Compact / True genomic) | Switch between intron-clamped and true-scale |
| Show UTRs checkbox | Toggle UTRs |
| Reset zoom / Fit to domain buttons | Quick presets |

## Embedding in a Jupyter notebook

```python
import prot2exon as p2e
from prot2exon.plot import _segments_from_dataframe

segs = _segments_from_dataframe(result.isoform)["TP53_DBD"]
p2e.render_interactive_jupyter(
    segs,
    plot_height=160,                            # main track height in px
    link_template="https://www.ensembl.org/Homo_sapiens/Transcript/ProteinSummary?p={protein_id}",
)
```

By default the iframe auto-resizes to its content (no inner scrollbar) via a postMessage handshake between the iframe and a small inline listener. For static exports where postMessage can't fire (e.g. nbconvert PDFs), pass `height=...` to pin a fixed iframe height instead.

## Picking a renderer

- **Static figure for a paper:** matplotlib (`--out file.pdf`). Vector output is print-ready.
- **Interactive embed in a docs site / dashboard:** the interactive viewer (`--html-interactive file.html`). Smallest payload, no CDN, no plotly bundle.
- **Quick exploration with hover tooltips and the plotly UI:** plotly (`--html file.html`).
- **Notebook-friendly:** `render_interactive_jupyter()` (inline) or `--out file.png` then `IPython.display.Image`.
