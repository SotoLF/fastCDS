"""Python client for the prot2exon C++ mapper.

Quick start
-----------

>>> import prot2exon as p2g
>>> mapper = p2g.Mapper(index="human.idx")
>>> result = mapper.map("ENSP00000269305", aa_start=10, aa_end=50,
...                     domain_id="AD1")
>>> result.summary
>>> p2g.plot(result, input_id="AD1", out="AD1.pdf")

For one-off calls (creates a Mapper internally):

>>> result = p2g.map_query("ENSP00000269305", 10, 50, "AD1",
...                        index="human.idx")

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

from ._client import Mapper, map_query
from ._result import MappingResult, read_results_dir
from .plot import plot, plot_isoform

__version__ = "2.2.0"
__all__ = [
    "Mapper",
    "MappingResult",
    "map_query",
    "plot",
    "plot_isoform",
    "read_results_dir",
]
