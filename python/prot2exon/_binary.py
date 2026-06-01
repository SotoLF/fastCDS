"""Locate the prot2exon C++ binary across the install layouts we support.

Search order:
  1. ``$PROT2EXON_BIN`` (explicit override)
  2. the wheel-bundled binary at ``prot2exon/_bin/prot2exon-core`` (pip install)
  3. ``<repo>/build/prot2exon`` then ``<repo>/bin/prot2exon`` (source checkout)
  4. ``prot2exon-core`` / ``prot2exon`` on ``$PATH`` (conda, manual install)
"""

from __future__ import annotations

import os
import shutil
import sys
from pathlib import Path


def find_binary(*, required: bool = True) -> str | None:
    chain: list[str] = []

    env = os.environ.get("PROT2EXON_BIN")
    if env:
        chain.append(f"PROT2EXON_BIN={env}")
        if os.access(env, os.X_OK):
            return env

    here = Path(__file__).resolve()
    pkg = here.parent                       # .../prot2exon
    candidates = [
        pkg / "_bin" / "prot2exon-core",        # wheel-bundled (POSIX)
        pkg / "_bin" / "prot2exon-core.exe",    # wheel-bundled (Windows)
        here.parents[2] / "build" / "prot2exon",  # dev build
        here.parents[2] / "bin" / "prot2exon",     # repo wrapper
    ]
    for cand in candidates:
        chain.append(str(cand))
        if cand.exists() and os.access(cand, os.X_OK):
            return str(cand)

    for name in ("prot2exon-core", "prot2exon"):
        path = shutil.which(name)
        chain.append(f"PATH:{name}")
        if path:
            return path

    if required:
        raise FileNotFoundError(
            "Could not locate the prot2exon C++ binary. Searched (in order): "
            + ", ".join(chain)
            + ". Install a wheel (`pip install prot2exon`), set PROT2EXON_BIN, "
            "or build the binary from source."
        )
    return None
