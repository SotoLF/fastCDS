"""Unified ``prot2exon`` command — the cross-platform dispatcher.

This is the console-script entry point installed by pip. It mirrors the repo's
``bin/prot2exon`` shell wrapper but in pure Python, so it works the same on
Linux / macOS / Windows:

    prot2exon index ...   -> C++ binary
    prot2exon map   ...   -> C++ binary
    prot2exon fetch ...   -> Python (prot2exon.fetch)
    prot2exon plot  ...   -> Python (prot2exon.plot)
    prot2exon --version | --help | <anything else>  -> C++ binary
"""

from __future__ import annotations

import subprocess
import sys

from ._binary import find_binary


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)

    if argv:
        if argv[0] == "plot":
            from .plot import main as plot_main
            return plot_main(argv[1:])
        if argv[0] == "fetch":
            from .fetch import main as fetch_main
            return fetch_main(argv[1:])

    # index / map / --version / --help / no args -> the C++ binary.
    binary = find_binary()
    return subprocess.call([binary, *argv])


if __name__ == "__main__":
    sys.exit(main())
