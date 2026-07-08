"""Download a ready-to-query fastCDS index from Zenodo.

`fastCDS fetch <target>` pulls a pre-built binary index from our Zenodo
deposit — one sha256-verified HTTPS download, no GTF parse + build:

    fastCDS fetch human                 # pre-built GENCODE v49 index
    fastCDS fetch mouse                 # pre-built GENCODE vM34 index
    fastCDS fetch yeast                 # pre-built RefSeq R64 index
    fastCDS fetch list                  # what's available

`fetch` is *only* the "pull a pre-built index off Zenodo" path. To use any
other annotation — a different release, a non-model species, a custom GTF —
download that GTF and build the index yourself with `fastCDS index` (see the
"Building an index" wiki page); it's a one-time ~15 s step.

Indexes land in `~/.cache/fastCDS/` by default (override with `--out` or
`--cache-dir`). The final index path is printed on stdout so it pipes straight
into the mapper.
"""

from __future__ import annotations

import argparse
import os
import shutil
import sys
import urllib.request
from dataclasses import dataclass
from pathlib import Path


def _default_cache_dir() -> Path:
    xdg = os.environ.get("XDG_CACHE_HOME")
    base = Path(xdg) if xdg else Path.home() / ".cache"
    return base / "fastCDS"


# ------------------------------------------------------------------------ #
# Pre-built indexes hosted on Zenodo — the only source `fetch` uses.
# A single HTTPS download of the binary index (typically ~10x smaller than the
# source GTF); ready to use immediately, no parse + build.
#
# Published Zenodo deposit: record 21266614 (version DOI 10.5281/zenodo.21266614,
# concept DOI 10.5281/zenodo.21266613). The per-file sha256 below are verified
# after download; a mismatch aborts with a pointer to `fastCDS index`.
# ------------------------------------------------------------------------ #

_ZENODO_BASE = "https://zenodo.org/records/21266614/files"


@dataclass(frozen=True)
class ZenodoIndex:
    filename: str
    sha256: str
    short_desc: str

    @property
    def url(self) -> str:
        return f"{_ZENODO_BASE}/{self.filename}"

    @property
    def published(self) -> bool:
        return "<RECORD>" not in self.url


ZENODO_IDX: dict[str, ZenodoIndex] = {
    "human": ZenodoIndex(
        "gencode_v49_human.idx",
        "ed848d78125dc795fa86a0af5402cb08ad679626fb153dda7a8ff2d6b47844f7",
        "GENCODE human v49 (~298 MB) — current human, most users"),
    "mouse": ZenodoIndex(
        "gencode_vM34_mouse.idx",
        "a8b22d9e229643903fc2ae9e7b867c7b9a72a07187dd26874f8a331f4213d8e9",
        "GENCODE mouse vM34 (~73 MB) — current mouse (GRCm39)"),
    "mouse-vm25": ZenodoIndex(
        "gencode_vM25_mouse.idx",
        "14610af24f0fe24f1d0282f3903f7cba7f657d8afc3834ae65b09ae7b7197ce9",
        "GENCODE mouse vM25 (~72 MB) — last GRCm38/mm10 release"),
    "human-v86": ZenodoIndex(
        "ensembl_v86_human.idx",
        "5999c3c4fdfb16517b0a687d3cb2ecff424ee4a3fa4019ab7825321c4bb6f25a",
        "Ensembl 86 human (~87 MB) — matches EnsDb.Hsapiens.v86, validation"),
    "human-v95": ZenodoIndex(
        "ensembl_v95_human.idx",
        "1d8e531cb23de01538c7390f51091966890dff945be32b6abc39cdd2c2274cac",
        "Ensembl 95 human (~92 MB) — matches TransVar's annotation, validation"),
    "human-v115": ZenodoIndex(
        "ensembl_v115_human.idx",
        "0b25b5ce07ac8fcf6116644c01c2b61f23dad94f8b4421aebf57cdb225d5f3a8",
        "Ensembl 115 human (~265 MB) — current Ensembl, REST validation + Pfam atlas"),
    "yeast": ZenodoIndex(
        "refseq_R64_yeast.idx",
        "201aeff2539d7b54ad82d09da69bc4ed0c2cf2e97454f1a3577b8df99d3b490b",
        "NCBI RefSeq S. cerevisiae R64 (~1.4 MB) — walkthrough notebook"),
}

# Targets the CLI exposes as subcommands (besides `list`).
TARGETS = tuple(ZENODO_IDX)


# ------------------------------------------------------------------------ #
# Download + verify
# ------------------------------------------------------------------------ #

def _download(url: str, dest: Path) -> None:
    sys.stderr.write(f"downloading {url}\n  -> {dest}\n")
    dest.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(url) as resp, open(dest, "wb") as out:
        shutil.copyfileobj(resp, out, length=1024 * 1024)
    sys.stderr.write(f"  done ({dest.stat().st_size / 1024 / 1024:.1f} MB)\n")


def _verify_sha256(path: Path, expected: str) -> None:
    import hashlib
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    got = h.hexdigest()
    if got != expected:
        sys.exit(
            f"error: sha256 mismatch for {path}\n"
            f"  expected: {expected}\n"
            f"  got:      {got}\n"
            f"Re-run with --force to redownload."
        )
    sys.stderr.write(f"sha256 OK ({got[:16]}...)\n")


_UNPUBLISHED_HINT = (
    "The Zenodo deposit isn't published yet (the record id is still a "
    "placeholder).\nBuild the index locally instead — download the GTF and run "
    "`fastCDS index`:\n"
    "  see https://github.com/SotoLF/fastCDS/wiki/Index"
)


def _cmd_get(args: argparse.Namespace) -> Path:
    """Download the pre-built index for `args.target` and return its .idx path."""
    target = args.target
    cache = Path(args.cache_dir or _default_cache_dir())
    force = args.force

    zen = ZENODO_IDX.get(target)
    if zen is None:
        sys.exit(f"error: no pre-built index for {target!r}. "
                 f"Available: {', '.join(sorted(ZENODO_IDX))}. "
                 f"For any other annotation, build it with `fastCDS index`.")
    if not zen.published:
        sys.exit(f"error: {target!r} — {_UNPUBLISHED_HINT}")

    out = Path(args.out) if args.out else (cache / f"{target}.idx")
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.exists() and not force:
        sys.stderr.write(f"using cached index: {out}  (--force to re-download)\n")
        _verify_sha256(out, zen.sha256)
        print(out)
        return out
    _download(zen.url, out)
    _verify_sha256(out, zen.sha256)
    sys.stderr.write(f"ready: {out}\n")
    print(out)
    return out


def _cmd_list() -> int:
    print("Pre-built indexes available from Zenodo (`fastCDS fetch <target>`):")
    print("=" * 72)
    for name in sorted(ZENODO_IDX):
        mark = "" if ZENODO_IDX[name].published else "   [not published yet]"
        print(f"  {name:<12s}  {ZENODO_IDX[name].short_desc}{mark}")
    print()
    print("Need a different release, species, or a custom GTF? There's no")
    print("pre-built index to fetch — download that GTF and build one with")
    print("`fastCDS index` (see the 'Building an index' wiki page).")
    return 0


# ------------------------------------------------------------------------ #
# Python API
# ------------------------------------------------------------------------ #

def fetch_index(
    target: str,
    *,
    cache_dir: str | Path | None = None,
    out: str | Path | None = None,
    force: bool = False,
    quiet: bool = False,
) -> Path:
    """Download a pre-built fastCDS index from Zenodo and return its path.

    Python mirror of ``fastCDS fetch``. Only pre-built Zenodo targets are
    available (``fastCDS fetch list``); for any other annotation, build an
    index from a GTF with :func:`build_index` / ``fastCDS index``.

    Parameters
    ----------
    target : {"human", "mouse", "mouse-vm25", "human-v86", "human-v95", "human-v115", "yeast"}
        A built-in pre-built target (see ``fastCDS fetch list``).
    cache_dir, out, force, quiet
        As per the CLI flags.

    Examples
    --------
    >>> idx = fetch_index("human")     # pre-built GENCODE v49 index, from Zenodo
    >>> idx = fetch_index("yeast")     # pre-built RefSeq R64 index
    """
    ns = argparse.Namespace(
        target=target,
        cache_dir=str(cache_dir) if cache_dir else None,
        out=str(out) if out else None,
        force=force,
    )
    if quiet:
        import io
        old_err, old_out = sys.stderr, sys.stdout
        sys.stderr, sys.stdout = io.StringIO(), io.StringIO()
        try:
            return _cmd_get(ns)
        finally:
            sys.stderr, sys.stdout = old_err, old_out
    return _cmd_get(ns)


# ------------------------------------------------------------------------ #
# CLI
# ------------------------------------------------------------------------ #

def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="fastCDS fetch",
        description="Download a pre-built index from Zenodo. For any other "
                    "annotation, build one from a GTF with `fastCDS index`.")
    sub = p.add_subparsers(dest="target", required=True)

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--out", help="Output .idx path (default: cache dir)")
    common.add_argument("--cache-dir",
                        help="Override cache dir (default: ~/.cache/fastCDS)")
    common.add_argument("--force", action="store_true",
                        help="Re-download even if the index already exists")

    sub.add_parser("list", help="Show available pre-built indexes")
    for name in TARGETS:
        sub.add_parser(name, parents=[common], help=ZENODO_IDX[name].short_desc)
    return p


def main(argv: list[str] | None = None) -> int:
    parser = build_argparser()
    args = parser.parse_args(argv)
    if args.target == "list":
        return _cmd_list()
    _cmd_get(args)   # prints the path; raises SystemExit on error
    return 0


if __name__ == "__main__":
    sys.exit(main())
