"""One-command download + build for prot2exon indices.

Wraps the canonical `curl` + `gunzip` + `--build-index` sequence so users can go
from "nothing on disk" to "ready-to-query index" with one invocation:

    prot2exon fetch human --release 49
    prot2exon fetch mouse --release M34
    prot2exon fetch ensembl --species danio_rerio --release 115
    prot2exon fetch refseq --preset yeast
    prot2exon fetch list

By default the built index lands in `~/.cache/prot2exon/`. Override with
`--out`. The downloaded `.gtf.gz` is removed after the index is built unless
`--keep-gtf` is passed.

This is intentionally a thin shell wrapper around the same `curl` and
`prot2exon --build-index` you'd run by hand — no proprietary URL list, no
remote service to depend on. If a preset URL goes stale, you can pass
`--gtf-url URL` to override.
"""

from __future__ import annotations

import argparse
import gzip
import os
import shutil
import subprocess
import sys
import urllib.request
from dataclasses import dataclass
from pathlib import Path


def _default_cache_dir() -> Path:
    xdg = os.environ.get("XDG_CACHE_HOME")
    base = Path(xdg) if xdg else Path.home() / ".cache"
    return base / "prot2exon"


def _default_binary() -> str:
    # Mirror the search order used by bin/prot2exon: PROT2EXON_BIN env var,
    # then the dev build path, then PATH.
    env = os.environ.get("PROT2EXON_BIN")
    if env and os.access(env, os.X_OK):
        return env
    here = Path(__file__).resolve()
    dev = here.parents[2] / "build" / "prot2exon"
    if dev.is_file() and os.access(dev, os.X_OK):
        return str(dev)
    on_path = shutil.which("prot2exon-core") or shutil.which("prot2exon")
    if on_path:
        return on_path
    sys.exit("error: could not locate the prot2exon binary "
             "(set PROT2EXON_BIN or build it under ./build/)")


@dataclass(frozen=True)
class Preset:
    name: str
    url: str
    short_desc: str


# Canonical URLs we've verified end-to-end. Generic species lookups use the
# Ensembl URL pattern below and don't need a preset entry.
PRESETS: dict[str, Preset] = {
    # GENCODE human — most users
    "human": Preset(
        name="human",
        url="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"
            "release_{release}/gencode.v{release}.basic.annotation.gtf.gz",
        short_desc="GENCODE human (basic annotation). Release placeholder: e.g. 49"),
    # GENCODE mouse
    "mouse": Preset(
        name="mouse",
        url="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/"
            "release_{release}/gencode.v{release}.basic.annotation.gtf.gz",
        short_desc="GENCODE mouse (basic annotation). Release placeholder: e.g. M34"),
    # NCBI RefSeq presets — accession is fixed per assembly, so no release knob.
    "yeast": Preset(
        name="yeast",
        url="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/"
            "GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.gtf.gz",
        short_desc="NCBI RefSeq Saccharomyces cerevisiae R64 (~2 MB compressed)"),
    "ecoli": Preset(
        name="ecoli",
        url="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/"
            "GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gtf.gz",
        short_desc="NCBI RefSeq Escherichia coli K-12 MG1655 (~1 MB compressed)"),
}

# Generic Ensembl URL pattern. `species` is the Ensembl lowercase name with
# underscores (e.g. danio_rerio). `cap_species` is the capitalised assembly
# prefix (e.g. Danio_rerio). Release is numeric.
ENSEMBL_PATTERN = (
    "https://ftp.ensembl.org/pub/release-{release}/gtf/{species}/"
    "{cap_species}.{assembly}.{release}.gtf.gz"
)


def _resolve_url(args: argparse.Namespace) -> tuple[str, str]:
    """Returns (url, suggested_filename_stem) based on CLI args."""
    if args.gtf_url:
        # When the user overrides the URL, still respect the subcommand for
        # the cache filename stem so a `fetch human --release 49 --gtf-url ...`
        # uses `human_v49.idx` rather than a generic name.
        if args.subcmd in ("human", "mouse") and getattr(args, "release", None):
            stem = f"{args.subcmd}_v{args.release}"
        elif args.subcmd == "ensembl" and getattr(args, "species", None):
            stem = f"{args.species}_v{getattr(args, 'release', 'custom')}"
        elif args.subcmd == "refseq" and getattr(args, "preset", None):
            stem = args.preset
        else:
            stem = "custom"
        return args.gtf_url, stem
    if args.subcmd in ("human", "mouse"):
        preset = PRESETS[args.subcmd]
        if not args.release:
            sys.exit(f"error: {args.subcmd} requires --release "
                     f"(e.g. 49 for human, M34 for mouse)")
        url = preset.url.format(release=args.release)
        stem = f"{args.subcmd}_v{args.release}"
        return url, stem
    if args.subcmd == "ensembl":
        if not (args.species and args.release and args.assembly):
            sys.exit("error: ensembl requires --species, --assembly, --release "
                     "(e.g. --species danio_rerio --assembly GRCz11 --release 115)")
        cap_species = args.species[:1].upper() + args.species[1:]
        url = ENSEMBL_PATTERN.format(
            release=args.release, species=args.species.lower(),
            cap_species=cap_species, assembly=args.assembly)
        stem = f"{args.species}_v{args.release}"
        return url, stem
    if args.subcmd == "refseq":
        if not args.preset or args.preset not in PRESETS:
            sys.exit(f"error: refseq --preset must be one of "
                     f"{sorted(k for k, p in PRESETS.items() if 'NCBI' in p.short_desc)}, "
                     f"or pass --gtf-url to download a custom RefSeq GTF")
        preset = PRESETS[args.preset]
        return preset.url, args.preset
    sys.exit(f"error: unknown subcommand: {args.subcmd}")


def _download(url: str, dest: Path) -> None:
    sys.stderr.write(f"downloading {url}\n  -> {dest}\n")
    dest.parent.mkdir(parents=True, exist_ok=True)
    # Use urllib so this works with zero external deps. Stream to disk to
    # avoid holding the entire GTF in memory.
    with urllib.request.urlopen(url) as resp, open(dest, "wb") as out:
        shutil.copyfileobj(resp, out, length=1024 * 1024)
    sys.stderr.write(f"  done ({dest.stat().st_size / 1024 / 1024:.1f} MB)\n")


def _gunzip(src: Path, dst: Path) -> None:
    sys.stderr.write(f"gunzipping {src} -> {dst}\n")
    with gzip.open(src, "rb") as fin, open(dst, "wb") as fout:
        shutil.copyfileobj(fin, fout, length=8 * 1024 * 1024)


def _build_index(binary: str, gtf: Path, index: Path) -> None:
    sys.stderr.write(f"building index -> {index}\n")
    proc = subprocess.run(
        [binary, "--gtf", str(gtf), "--build-index", "--index", str(index)],
        stderr=None,  # surface the binary's progress directly
    )
    if proc.returncode != 0:
        sys.exit(f"error: --build-index exited {proc.returncode}")


def _cmd_list() -> int:
    print("Built-in presets (also see `prot2exon fetch <preset> --help`):")
    print()
    for name in sorted(PRESETS):
        p = PRESETS[name]
        print(f"  {name:<8s}  {p.short_desc}")
    print()
    print("Generic Ensembl species (any organism with a GTF on ftp.ensembl.org):")
    print("  prot2exon fetch ensembl --species danio_rerio "
          "--assembly GRCz11 --release 115")
    print()
    print("Any GTF URL (override the preset):")
    print("  prot2exon fetch human --release 49 "
          "--gtf-url https://your.host/path/custom.gtf.gz")
    return 0


def fetch_index(
    species: str,
    *,
    release: str | None = None,
    assembly: str | None = None,
    gtf_url: str | None = None,
    cache_dir: str | Path | None = None,
    out: str | Path | None = None,
    binary: str | None = None,
    force: bool = False,
    keep_gtf: bool = False,
    quiet: bool = False,
) -> Path:
    """Download a GTF and build an index. Returns the path to the .idx.

    Smooth Python equivalent of ``prot2exon fetch``. Re-uses the cached
    download / gunzip / index whenever it can (pass ``force=True`` to
    re-fetch from upstream).

    Parameters
    ----------
    species : {"human", "mouse", "yeast", "ecoli", ...}
        Built-in presets are listed by ``prot2exon.fetch.PRESETS``. Any other
        value is treated as an Ensembl species name (e.g. ``"danio_rerio"``)
        and resolved against the generic Ensembl URL pattern — pass
        ``assembly`` and ``release`` in that case.
    release : str
        Required for ``human`` / ``mouse`` (GENCODE release: e.g. "49", "M34")
        and Ensembl species (release number: e.g. "115").
    assembly : str
        Required for Ensembl species (e.g. "GRCz11" for zebrafish).
    gtf_url : str
        Override the upstream URL entirely. Useful for mirrors or custom
        annotations.
    cache_dir : str | Path
        Where to keep the downloaded GTF / built index. Defaults to
        ``~/.cache/prot2exon/``.
    out : str | Path
        Override the index path (default: ``{cache_dir}/{species}_v{release}.idx``).
    binary : str
        Path to the ``prot2exon`` binary. Auto-discovered if omitted.
    force : bool
        Re-download and re-build even if cached. Default False.
    keep_gtf : bool
        Don't delete the uncompressed .gtf after building. Default False —
        the .gtf.gz stays cached for cheap rebuilds, the .gtf is removed.
    quiet : bool
        Suppress the progress lines that the CLI prints to stderr.

    Returns
    -------
    Path : path to the built (or cached) index.

    Examples
    --------
    >>> idx = fetch_index("human", release="49")
    >>> idx = fetch_index("mouse", release="M34")
    >>> idx = fetch_index("yeast")                 # RefSeq preset
    >>> idx = fetch_index("danio_rerio", release="115", assembly="GRCz11")
    """
    # Build a minimal argparse-like Namespace so _resolve_url can reuse its
    # routing logic. Note _resolve_url's subcmd values: human, mouse, ensembl, refseq.
    if species in ("human", "mouse"):
        subcmd = species
        preset = None
    elif species in PRESETS and species not in ("human", "mouse"):
        # RefSeq presets (yeast, ecoli) — pass through `refseq --preset <species>`.
        subcmd = "refseq"
        preset = species
    else:
        subcmd = "ensembl"
        preset = None

    ns = argparse.Namespace(
        subcmd=subcmd,
        preset=preset,
        species=species if subcmd == "ensembl" else None,
        release=release,
        assembly=assembly,
        gtf_url=gtf_url,
        cache_dir=str(cache_dir) if cache_dir else None,
        out=str(out) if out else None,
        binary=binary,
        force=force,
        keep_gtf=keep_gtf,
    )

    if quiet:
        # Swap stderr + stdout for sinks so the internal _download / _gunzip
        # / _build_index helpers stay quiet without us having to thread a
        # `quiet=` kwarg through every one of them.
        import io
        old_err, old_out = sys.stderr, sys.stdout
        sys.stderr, sys.stdout = io.StringIO(), io.StringIO()
        try:
            _cmd_fetch(ns)
        finally:
            sys.stderr, sys.stdout = old_err, old_out
    else:
        _cmd_fetch(ns)

    url, stem = _resolve_url(ns)
    cache = Path(cache_dir or _default_cache_dir())
    return Path(out) if out else (cache / f"{stem}.idx")


def _cmd_fetch(args: argparse.Namespace) -> int:
    binary = args.binary or _default_binary()
    url, stem = _resolve_url(args)
    cache = Path(args.cache_dir or _default_cache_dir())
    cache.mkdir(parents=True, exist_ok=True)
    gtf_gz = cache / f"{stem}.gtf.gz"
    gtf = cache / f"{stem}.gtf"
    out_idx = Path(args.out) if args.out else (cache / f"{stem}.idx")

    if out_idx.exists() and not args.force:
        sys.stderr.write(f"using cached index: {out_idx}  (--force to rebuild)\n")
        print(out_idx)
        return 0

    if not gtf_gz.exists() or args.force:
        _download(url, gtf_gz)
    else:
        sys.stderr.write(f"reusing cached download: {gtf_gz}\n")

    if not gtf.exists() or args.force:
        _gunzip(gtf_gz, gtf)
    else:
        sys.stderr.write(f"reusing cached gunzip:   {gtf}\n")

    _build_index(binary, gtf, out_idx)

    if not args.keep_gtf:
        # The .gtf (uncompressed) can be huge (3 GB for human); the .gtf.gz
        # stays in cache so a rebuild is cheap.
        try:
            gtf.unlink()
            sys.stderr.write(f"removed {gtf} (pass --keep-gtf to retain)\n")
        except FileNotFoundError:
            pass

    sys.stderr.write(f"ready: {out_idx}\n")
    print(out_idx)
    return 0


def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="prot2exon fetch",
        description="Download and build a prot2exon index in one command.")
    sub = p.add_subparsers(dest="subcmd", required=True)

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--out", help="Output .idx path (default: cache dir)")
    common.add_argument("--cache-dir",
                        help="Override cache dir (default: ~/.cache/prot2exon)")
    common.add_argument("--binary",
                        help="Path to the prot2exon binary (default: auto-discover)")
    common.add_argument("--gtf-url",
                        help="Override the preset URL with a custom GTF location")
    common.add_argument("--keep-gtf", action="store_true",
                        help="Keep the uncompressed GTF after indexing")
    common.add_argument("--force", action="store_true",
                        help="Re-download / rebuild even if outputs exist")

    sub.add_parser("list", help="Show built-in presets")

    sp_idx = sub.add_parser("index",
        help="Download a pre-built index directly (skip GTF parse + build)")
    sp_idx.add_argument("--url", required=True,
        help="URL of the pre-built .idx (http://, https://, or file://)")
    sp_idx.add_argument("--out", required=True,
        help="Output path for the downloaded index")
    sp_idx.add_argument("--sha256",
        help="Optional sha256 to verify the download against")
    sp_idx.add_argument("--force", action="store_true",
        help="Re-download even if --out already exists")

    for sp_name in ("human", "mouse"):
        sp = sub.add_parser(sp_name, parents=[common],
                            help=PRESETS[sp_name].short_desc)
        sp.add_argument("--release", required=True,
                        help="Release tag (e.g. 49 for human, M34 for mouse)")

    sp = sub.add_parser("ensembl", parents=[common],
                        help="Generic Ensembl species lookup")
    sp.add_argument("--species", required=True,
                    help="Ensembl species name (e.g. danio_rerio)")
    sp.add_argument("--assembly", required=True,
                    help="Assembly tag (e.g. GRCz11)")
    sp.add_argument("--release", required=True, help="Ensembl release number")

    sp = sub.add_parser("refseq", parents=[common],
                        help="NCBI RefSeq presets (yeast, ecoli) or custom URL")
    sp.add_argument("--preset", help="Preset name (yeast, ecoli, ...)")

    return p


def _cmd_index(args: argparse.Namespace) -> int:
    """Download a pre-built .idx directly. No GTF parse step."""
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    if out.exists() and not args.force:
        sys.stderr.write(f"using cached index: {out}  (--force to re-download)\n")
        if args.sha256:
            _verify_sha256(out, args.sha256)
        print(out)
        return 0

    _download(args.url, out)

    if args.sha256:
        _verify_sha256(out, args.sha256)

    sys.stderr.write(f"ready: {out}\n")
    print(out)
    return 0


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
    sys.stderr.write(f"sha256 OK ({got[:16]}…)\n")


def main(argv: list[str] | None = None) -> int:
    parser = build_argparser()
    args = parser.parse_args(argv)
    if args.subcmd == "list":
        return _cmd_list()
    if args.subcmd == "index":
        return _cmd_index(args)
    return _cmd_fetch(args)


if __name__ == "__main__":
    sys.exit(main())
