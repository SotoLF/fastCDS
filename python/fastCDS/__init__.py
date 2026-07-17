"""Python client for the fastCDS C++ mapper.

Quick start
-----------

>>> import fastCDS as fc
>>> mapper = fc.Mapper(index="human.idx")
>>> result = mapper.map("ENSP00000269305", aa_start=10, aa_end=50,
...                     domain_id="AD1")
>>> result.summary
>>> fc.plot(result, input_id="AD1", out="AD1.pdf")

For one-off calls (creates a Mapper internally):

>>> result = fc.map_query("ENSP00000269305", aa_start=10, aa_end=50,
...                       domain_id="AD1", index="human.idx")

Batch:

>>> queries = [
...     {"protein_id": "ENSP00000269305", "aa_start": 10, "aa_end": 50,
...      "domain_id": "AD1"},
...     {"protein_id": "ENST00000303562", "aa_start": 5, "aa_end": 100,
...      "domain_id": "RD1"},
... ]
>>> result = mapper.map_batch(queries)

The wrapper shells out to the compiled C++ binary for every call. For many
queries, prefer `map_batch` over many `map` calls: each call reloads the
index from disk.
"""

from ._client import Mapper, map_query, build_index
from ._result import MappingResult, read_results_dir
from .plot import plot, plot_all, plot_isoform
from ._interactive_html import render_interactive_html, render_interactive_jupyter
from ._interactive_stack_html import (
    render_interactive_html_stack, render_interactive_jupyter_stack,
)
from .fetch import fetch_index
from . import prepare


__version__ = "2.3.0"
__all__ = [
    "Mapper",
    "MappingResult",
    "map_query",
    "build_index",
    "plot",
    "plot_all",
    "plot_isoform",
    "read_results_dir",
    "render_interactive_html",
    "render_interactive_jupyter",
    "render_interactive_html_stack",
    "render_interactive_jupyter_stack",
    "fetch_index",
    "prepare",
]
