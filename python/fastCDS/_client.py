"""Thin client around the fastCDS C++ binary.

The binary remains the source of truth for all mapping logic. This module
just:

  1. assembles a BED file from Python-side queries,
  2. invokes the binary with the right flags,
  3. reads the produced TSV / BED / BED12 files back into pandas DataFrames.

By design, nothing is reimplemented in Python — if you change the C++ behavior,
the Python wrapper picks it up automatically because it only ever reads the
files the binary writes.

Binary discovery order (so `Mapper(index=...)` "just works"):

  1. The `binary=` constructor argument, if given.
  2. The `FASTCDS_BIN` env var.
  3. `<repo>/build/fastCDS`, relative to this file (development).
  4. The `bin/fastCDS` shell wrapper, relative to this file.
  5. PATH lookup of `fastCDS-core`, then `fastCDS`.

If none of the above resolves to an executable file, a clear FileNotFoundError
is raised with the discovery chain that was tried.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any, Iterable

import pandas as pd

from ._result import MappingResult, read_results_dir


_VALID_OUTPUTS = {"coding", "introns", "span", "isoform", "bed12", "all"}


def _discover_binary() -> str:
    """Find the C++ binary (shared logic in ``fastCDS._binary``)."""
    from ._binary import find_binary
    return find_binary()


def build_index(gtf, out=None, *, binary: str | None = None,
                force: bool = False) -> Path:
    """Build a binary index from a local GTF — Python mirror of ``fastCDS index``.

    Parameters
    ----------
    gtf : str | os.PathLike
        Path to the (uncompressed) GTF file to index.
    out : str | os.PathLike | None
        Output ``.idx`` path. Defaults to ``<gtf>`` with a ``.idx`` suffix.
    binary : str | None
        Path to the fastCDS binary. Auto-discovered if omitted.
    force : bool
        Rebuild even if ``out`` already exists. Default False (cached).

    Returns
    -------
    Path : path to the built (or cached) index.

    Examples
    --------
    >>> import fastCDS as fc
    >>> idx = fc.build_index("combined.gtf", out="human.idx")
    """
    gtf = Path(gtf)
    if not gtf.exists():
        raise FileNotFoundError(f"GTF not found: {gtf}")
    out = Path(out) if out else gtf.with_suffix(".idx")
    if out.exists() and not force:
        return out
    binary = binary or _discover_binary()
    proc = subprocess.run(
        [binary, "index", "--gtf", str(gtf), "--out", str(out)],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
    )
    if proc.returncode != 0:
        raise RuntimeError(
            f"`fastCDS index` exited {proc.returncode}\n{proc.stderr}")
    return out


def _nan_to_none(x: Any) -> Any:
    """Treat pandas NaN as missing. When a DataFrame is converted via
    ``to_dict(orient="records")``, absent values surface as float NaN — which
    is truthy in Python and does not equal None, so the rest of this module
    would otherwise mis-handle it (silently stringifying "nan" as an id, or
    blowing up in ``int(NaN)``)."""
    return None if isinstance(x, float) and x != x else x


def _query_to_bed_row(q: dict[str, Any]) -> str:
    """Turn one query dict into a BED-like line.

    Accepts both `protein_id` (the canonical name in our schema) and `id`
    (because users naturally type that). Either ENSP or ENST.
    """
    pid = (_nan_to_none(q.get("protein_id"))
           or _nan_to_none(q.get("id"))
           or _nan_to_none(q.get("transcript_id")))
    if not pid:
        raise ValueError(
            f"query is missing protein_id / id / transcript_id: {q!r}")
    aa_start = _nan_to_none(q.get("aa_start"))
    aa_end = _nan_to_none(q.get("aa_end"))
    domain_id = _nan_to_none(q.get("domain_id"))

    # No-domain (structure-only): aa_start == aa_end == 0 OR both omitted.
    if aa_start is None and aa_end is None:
        aa_start, aa_end = 0, 0
    if (aa_start in (0, None)) != (aa_end in (0, None)):
        raise ValueError(
            f"either both aa_start and aa_end must be set, or neither: {q!r}")

    parts = [str(pid), str(int(aa_start or 0)), str(int(aa_end or 0))]
    # The C++ parser is whitespace-tokenizing, so we must always have 3 numeric
    # columns when a domain_id is present (otherwise it'd shift left).
    if domain_id:
        parts.append(str(domain_id))
    return "\t".join(parts)


class Mapper:
    """Repeatedly map domain queries against a fastCDS index.

    The index file is reloaded by the C++ binary on every invocation, so
    repeated `.map()` calls each pay the index-load cost (~1.5 s on the
    human GENCODE index). Prefer `.map_batch()` when you have multiple
    queries.

    Parameters
    ----------
    index : str | os.PathLike
        Path to the binary index built with `fastCDS index`.
    binary : str | None
        Path to the C++ binary. Auto-discovered if omitted (see module
        docstring for the search order).
    threads : int
        Forwarded to `--threads N` for parallel batch mapping.
    batch_size : int
        Forwarded to `--batch-size N`. 0 (default) processes all queries
        at once with the lowest overhead and the highest peak RAM. A
        positive N streams results to disk in chunks of N queries, bounding
        peak memory at roughly O(N * per-query result size). Useful for
        million-query runs on memory-constrained machines.
    verbose : bool
        If True, the binary's stderr is streamed to this process's stderr;
        otherwise it's captured and only surfaced on failure.
    """

    def __init__(self, index, binary=None, *, threads: int = 1,
                 batch_size: int = 0,
                 verbose: bool = False) -> None:
        self.index = os.fspath(index)
        if not os.path.exists(self.index):
            raise FileNotFoundError(f"Index file not found: {self.index}")
        self.binary = binary or _discover_binary()
        if not os.access(self.binary, os.X_OK):
            raise FileNotFoundError(
                f"Binary not executable: {self.binary}")
        self.threads = int(threads)
        if int(batch_size) < 0:
            raise ValueError(f"batch_size must be >= 0, got {batch_size!r}")
        self.batch_size = int(batch_size)
        self.verbose = bool(verbose)

    # ------------------------------------------------------------------ #
    # Public API
    # ------------------------------------------------------------------ #
    def map(self, protein_id: str,
            aa_start: int | None = None,
            aa_end: int | None = None,
            domain_id: str | None = None,
            *, output: str = "all",
            keep_outputs: str | bool = False) -> MappingResult:
        """Single-query convenience wrapper around `.map_batch`."""
        return self.map_batch([{
            "protein_id": protein_id,
            "aa_start": aa_start,
            "aa_end": aa_end,
            "domain_id": domain_id,
        }], output=output, keep_outputs=keep_outputs)

    def map_batch(self,
                  queries: Iterable[dict[str, Any]] | pd.DataFrame,
                  *, output: str = "all",
                  keep_outputs: str | bool = False) -> MappingResult:
        """Map a batch of queries in a single binary invocation.

        Parameters
        ----------
        queries : iterable of dict | DataFrame
            Each query supports keys ``protein_id`` (or ``id`` / ``transcript_id``),
            ``aa_start``, ``aa_end``, ``domain_id``. A query without aa coordinates
            (or with both set to 0) is processed in no-domain / structure-only mode.
            A DataFrame is unpacked row-wise.
        output : {"coding","introns","span","isoform","bed12","all"}
            Forwarded to ``--output``. Default: ``all``.
        keep_outputs : bool | str
            False (default) — write the binary's outputs to a tempdir,
                read them back, then clean up.
            True            — write to a tempdir and leave the files in place;
                the returned MappingResult's ``out_dir`` field points at it.
            str             — use the given path as the output directory
                (created if missing) and leave the files in place.
        """
        if output not in _VALID_OUTPUTS:
            raise ValueError(
                f"output must be one of {sorted(_VALID_OUTPUTS)}, got {output!r}")

        if isinstance(queries, pd.DataFrame):
            queries = queries.to_dict(orient="records")
        queries = list(queries)
        if not queries:
            raise ValueError("queries is empty")

        # Build the BED.
        bed_lines = [_query_to_bed_row(q) for q in queries]

        # Resolve the output dir.
        if keep_outputs is False:
            out_dir_ctx: Any = tempfile.TemporaryDirectory(prefix="p2g_")
            cleanup = True
        else:
            cleanup = False
            if keep_outputs is True:
                out_dir = tempfile.mkdtemp(prefix="p2g_")
            else:
                out_dir = os.fspath(keep_outputs)
                os.makedirs(out_dir, exist_ok=True)
            out_dir_ctx = _NullCtx(out_dir)

        with out_dir_ctx as out_dir, \
             tempfile.NamedTemporaryFile(
                 "w", suffix=".bed", delete=False, prefix="p2g_query_") as bed_fh:
            bed_path = bed_fh.name
            bed_fh.write("# fastCDS batch query (auto-generated)\n")
            bed_fh.write("\n".join(bed_lines) + "\n")

        try:
            self._run(bed_path=bed_path, out_dir=out_dir, output=output)
            result = read_results_dir(out_dir)
            # When we don't keep_outputs, blank out the out_dir field because
            # the directory is about to disappear and a stale path would
            # mislead the user.
            if cleanup:
                result.out_dir = None
            return result
        finally:
            try:
                os.unlink(bed_path)
            except OSError:
                pass
            if cleanup:
                out_dir_ctx.cleanup()

    # ------------------------------------------------------------------ #
    # Internals
    # ------------------------------------------------------------------ #
    def _run(self, *, bed_path: str, out_dir: str, output: str) -> None:
        cmd = [
            self.binary,
            "map",
            "--index", self.index,
            "--bed", bed_path,
            "--out-dir", out_dir,
            "--output", output,
        ]
        if self.threads > 1:
            cmd += ["--threads", str(self.threads)]
        if self.batch_size > 0:
            cmd += ["--batch-size", str(self.batch_size)]
        if self.verbose:
            cmd += ["--verbose"]
        proc = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=None if self.verbose else subprocess.PIPE,
            text=True,
        )
        if proc.returncode != 0:
            err = proc.stderr or "(stderr captured nowhere)"
            raise RuntimeError(
                f"fastCDS exited {proc.returncode}\n"
                f"command: {' '.join(cmd)}\n"
                f"stderr:\n{err}"
            )


class _NullCtx:
    """Pretends to be a context manager for a path we don't own. Allows the
    same `with ... as out_dir:` block to handle both the tempdir-managed and
    user-supplied cases without branching."""

    def __init__(self, path: str) -> None:
        self.path = path

    def __enter__(self) -> str:
        return self.path

    def __exit__(self, *exc) -> None:
        pass

    def cleanup(self) -> None:  # parity with TemporaryDirectory
        pass


def map_query(protein_id: str,
              aa_start: int | None = None,
              aa_end: int | None = None,
              domain_id: str | None = None,
              *, index, binary=None,
              output: str = "all",
              keep_outputs: str | bool = False) -> MappingResult:
    """One-shot mapping without creating a Mapper instance yourself.

    Equivalent to ``Mapper(index=index, binary=binary).map(...)``.
    """
    return Mapper(index=index, binary=binary).map(
        protein_id, aa_start=aa_start, aa_end=aa_end, domain_id=domain_id,
        output=output, keep_outputs=keep_outputs)
