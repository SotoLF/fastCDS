"""Plot helpers and the standalone HTML viewers.

Drives the `fastCDS.plot` entry point and the interactive-HTML renderers
against the shared `out_all` isoform table. Visual fidelity isn't asserted —
these check that each path runs, honours its flags, and (for the offline
viewer) produces a self-contained file.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest

from conftest import REPO_ROOT, WRAPPER

# The wrapper shells out to `python3 -m fastCDS.plot`; make sure that resolves
# to the same interpreter running the tests (which has pandas/matplotlib) rather
# than a bare system python3 on PATH.
_WRAPPER_ENV = {**os.environ,
                "PATH": os.pathsep.join([str(Path(sys.executable).parent),
                                         os.environ.get("PATH", "")])}


def _isoform_tsv(out_all):
    return str(out_all / "isoform_structure.tsv")


def test_compact_genomic_pdf(out_all, tmp_path):
    """`--compact-genomic` clamps introns to a fixed display width; here we just
    confirm it runs and writes a non-empty PDF."""
    from fastCDS.plot import main
    pdf = tmp_path / "compact_genomic.pdf"
    rc = main(["--isoform", _isoform_tsv(out_all), "--input-id", "Q1_ENSP",
               "--out", str(pdf), "--compact-genomic"])
    assert rc == 0
    assert pdf.exists() and pdf.stat().st_size > 0


def test_compact_genomic_spliced_mutually_exclusive(out_all, tmp_path):
    """`--spliced` and `--compact-genomic` are mutually exclusive — the CLI
    rejects the combination with exit code 2."""
    from fastCDS.plot import main
    rc = main(["--isoform", _isoform_tsv(out_all), "--input-id", "Q1_ENSP",
               "--out", str(tmp_path / "mx.pdf"), "--spliced", "--compact-genomic"])
    assert rc == 2


def test_link_template_html(out_all, tmp_path):
    """`--out x.html --engine plotly --link-template` renders a plotly figure
    with the `{protein_id}` placeholder expanded into the linkout URL."""
    pytest.importorskip("plotly")
    from fastCDS.plot import main
    html = tmp_path / "link.html"
    rc = main(["--isoform", _isoform_tsv(out_all), "--input-id", "Q1_ENSP",
               "--out", str(html), "--engine", "plotly",
               "--link-template", "https://example.com/{protein_id}/entry"])
    assert rc == 0
    assert html.exists() and html.stat().st_size > 0
    text = html.read_text()
    assert "example.com" in text          # host present
    assert "ENSP1" in text                 # placeholder expanded
    assert "{protein_id}" not in text      # raw placeholder gone


def test_html_engine_js_is_self_contained(out_all, tmp_path):
    """`--out x.html` (default `--engine js`) writes the vanilla-JS viewer: no
    plotly CDN, the transcript id embedded in the JS payload, and a populated
    domains array."""
    from fastCDS.plot import main
    html = tmp_path / "interactive.html"
    rc = main(["--isoform", _isoform_tsv(out_all), "--input-id", "Q1_ENSP",
               "--out", str(html)])
    assert rc == 0
    assert html.exists() and html.stat().st_size > 0
    text = html.read_text()
    assert "cdn.plot.ly" not in text and "esm.sh" not in text
    assert "ENST1" in text
    assert '"domains":' in text and '"type": "Other"' in text


def test_unsupported_out_extension_is_rejected(out_all, tmp_path):
    """An out path with an unknown extension fails fast (rc 2), not a crash."""
    from fastCDS.plot import main
    rc = main(["--isoform", _isoform_tsv(out_all), "--input-id", "Q1_ENSP",
               "--out", str(tmp_path / "nope.jpeg")])
    assert rc == 2


def test_render_to_string_plot_height(out_all):
    """`_render_to_string` backs both the file writer and the Jupyter wrapper;
    the `plot_height` kwarg must reach the emitted CSS."""
    from fastCDS._interactive_html import _render_to_string
    from fastCDS.plot import load_isoform_tsv
    segs = load_isoform_tsv(_isoform_tsv(out_all))["Q1_ENSP"]
    out = _render_to_string(segs, plot_height=200)
    assert "height: 200px" in out
    assert "top: 27.5%" in out


def test_render_interactive_html_stack(out_all, tmp_path):
    """The stack viewer renders a dict of {input_id: segments} on one shared
    axis and embeds a STRUCTURES array."""
    from fastCDS import render_interactive_html_stack
    from fastCDS.plot import load_isoform_tsv
    by_id = load_isoform_tsv(_isoform_tsv(out_all))
    stack = tmp_path / "stack.html"
    render_interactive_html_stack(by_id, str(stack))
    assert stack.exists() and stack.stat().st_size > 0
    text = stack.read_text()
    assert "const STRUCTURES" in text and "iso-stack" in text


def test_render_interactive_jupyter_is_exported():
    """The Jupyter wrapper is part of the public API."""
    import fastCDS
    assert hasattr(fastCDS, "render_interactive_jupyter")


def test_cli_plotter_smoke(out_all, tmp_path):
    """The `fastCDS plot` CLI (wrapper binary) writes a non-empty PDF."""
    pdf = tmp_path / "Q1.pdf"
    proc = subprocess.run(
        [str(WRAPPER), "plot", "--isoform", _isoform_tsv(out_all),
         "--input-id", "Q1_ENSP", "--out", str(pdf)],
        capture_output=True, text=True, env=_WRAPPER_ENV,
    )
    assert proc.returncode == 0, proc.stderr
    assert pdf.exists() and pdf.stat().st_size > 0


def test_make_figure_1_renders(tmp_path):
    """make_figure_1.py renders the headline figure (PNG + PDF)
    without crashing, using the in-repo TP53 fixture."""
    fig_dir = tmp_path / "figures"
    bench = REPO_ROOT / "tutorial" / "reproduce_paper" / "benchmarks"
    proc = subprocess.run(
        [sys.executable, str(bench / "make_figure_1.py"),
         "--tp53-isoforms", str(REPO_ROOT / "tutorial" / "examples" / "tp53_isoforms.tsv"),
         "--out-dir", str(fig_dir)],
        capture_output=True, text=True,
    )
    assert proc.returncode == 0, proc.stderr
    png, pdf = fig_dir / "figure_1.png", fig_dir / "figure_1.pdf"
    assert png.exists() and png.stat().st_size > 10_000
    assert pdf.exists() and pdf.stat().st_size > 5_000
