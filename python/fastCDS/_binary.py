"""Locate the fastCDS C++ binary across the install layouts we support.

Search order:
  1. ``$FASTCDS_BIN`` (explicit override)
  2. the wheel-bundled binary at ``fastCDS/_bin/fastCDS-core`` (pip install)
  3. ``<repo>/build/fastCDS`` then ``<repo>/bin/fastCDS`` (source checkout)
  4. ``fastCDS-core`` / ``fastCDS`` on ``$PATH`` (conda, manual install)
"""

from __future__ import annotations

import os
import shutil
import sys
from pathlib import Path


def find_binary(*, required: bool = True) -> str | None:
    chain: list[str] = []

    env = os.environ.get("FASTCDS_BIN")
    if env:
        chain.append(f"FASTCDS_BIN={env}")
        if os.access(env, os.X_OK):
            return env

    here = Path(__file__).resolve()
    pkg = here.parent                       # .../fastCDS
    candidates = [
        pkg / "_bin" / "fastCDS-core",        # wheel-bundled (POSIX)
        pkg / "_bin" / "fastCDS-core.exe",    # wheel-bundled (Windows)
        here.parents[2] / "build" / "fastCDS",  # dev build
        here.parents[2] / "bin" / "fastCDS",     # repo wrapper
    ]
    for cand in candidates:
        chain.append(str(cand))
        if cand.exists() and os.access(cand, os.X_OK):
            return str(cand)

    for name in ("fastCDS-core", "fastCDS"):
        path = shutil.which(name)
        chain.append(f"PATH:{name}")
        if path:
            return path

    if required:
        raise FileNotFoundError(
            "Could not locate the fastCDS C++ binary. Searched (in order): "
            + ", ".join(chain)
            + ". Install a wheel (`pip install fastCDS`), set FASTCDS_BIN, "
            "or build the binary from source."
        )
    return None
