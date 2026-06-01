"""Get a ready-to-query prot2exon index in one command.

`prot2exon fetch <target>` downloads a pre-built index from our Zenodo deposit
by default — no GTF parse + build step. You only reach for a GTF when you want
something the deposit doesn't carry:

    prot2exon fetch human                 # pre-built GENCODE v49 index from Zenodo
    prot2exon fetch mouse                 # pre-built GENCODE vM34 index
    prot2exon fetch yeast                 # pre-built RefSeq R64 index
    prot2exon fetch list                  # what's available

    # Override the source when you need a different annotation:
    prot2exon fetch human --release 50              # build from upstream GENCODE v50
    prot2exon fetch human --gtf-url https://.../x.gtf.gz   # build from any GTF
    prot2exon fetch ensembl --species danio_rerio --assembly GRCz11 --release 115

Building from a local GTF you already have is the separate `prot2exon index`
command. `fetch` is specifically the "pull it off the network" path.

Indexes land in `~/.cache/prot2exon/` by default (override with `--out` or
`--cache-dir`). On a build, the `.gtf.gz` stays cached for cheap rebuilds and
the uncompressed `.gtf` is removed unless `--keep-gtf` is given. The final
index path is printed on stdout so it pipes straight into the mapper.
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
    # Shared discovery: env var, wheel-bundled binary, dev build, then PATH.
    from ._binary import find_binary
    b = find_binary(required=False)
    if not b:
        sys.exit("error: could not locate the prot2exon binary "
                 "(install a wheel, set PROT2EXON_BIN, or build under ./build/)")
    return b


# ------------------------------------------------------------------------ #
# Pre-built indexes hosted on Zenodo — the default source for `fetch`.
# A single HTTPS download of the binary index (typically ~10x smaller than
# the source GTF); ready to use immediately, no parse + build.
#
# The `<RECORD>` placeholder is filled in on first publish (one sed over the
# URLs + the per-file sha256 from the deposit's MANIFEST.tsv). Until then,
# `fetch <target>` transparently falls back to building from the upstream GTF
# for any target that has a GTF source below.
# ------------------------------------------------------------------------ #

_ZENODO_BASE = "https://zenodo.org/record/<RECORD>/files"


@dataclass(frozen=True)
class ZenodoIndex:
    filename: str
    sha256: str
    short_desc: str

    @property
    def url(self) -> str:
        return f"{_ZENODO_BASE}/{self.filename}"


ZENODO_IDX: dict[str, ZenodoIndex] = {
    "human": ZenodoIndex(
        "gencode_v49_human.idx",
        "ed848d78125dc795fa86a0af5402cb08ad679626fb153dda7a8ff2d6b47844f7",
        "GENCODE human v49 (~298 MB) — current human, most users"),
    "mouse": ZenodoIndex(
        "gencode_vM34_mouse.idx",
        "a8b22d9e229643903fc2ae9e7b867c7b9a72a07187dd26874f8a331f4213d8e9",
        "GENCODE mouse vM34 (~73 MB)"),
    "human-v86": ZenodoIndex(
        "ensembl_v86_human.idx",
        "5999c3c4fdfb1651"  # truncated — full hash filled in after upload
        "",
        "Ensembl 86 human (~87 MB) — matches EnsDb.Hsapiens.v86, validation"),
    "yeast": ZenodoIndex(
        "refseq_R64_yeast.idx",
        "201aeff2539d7b54ad82d09da69bc4ed0c2cf2e97454f1a3577b8df99d3b490b",
        "NCBI RefSeq S. cerevisiae R64 (~1.4 MB) — walkthrough notebook"),
}


# ------------------------------------------------------------------------ #
# Upstream GTF sources, used when the user overrides with --release / falls
# back from an unpublished Zenodo record. `{release}` is filled per call.
# ------------------------------------------------------------------------ #

@dataclass(frozen=True)
class GtfSource:
    url: str
    needs_release: bool
    short_desc: str


GTF_BUILD: dict[str, GtfSource] = {
    "human": GtfSource(
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"
        "release_{release}/gencode.v{release}.basic.annotation.gtf.gz",
        True, "GENCODE human — needs --release (e.g. 49)"),
    "mouse": GtfSource(
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/"
        "release_{release}/gencode.v{release}.basic.annotation.gtf.gz",
        True, "GENCODE mouse — needs --release (e.g. M34)"),
    "yeast": GtfSource(
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/"
        "GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.gtf.gz",
        False, "NCBI RefSeq Saccharomyces cerevisiae R64"),
    "ecoli": GtfSource(
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/"
        "GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gtf.gz",
        False, "NCBI RefSeq Escherichia coli K-12 MG1655"),
}

# Generic Ensembl URL pattern for `fetch ensembl --species ...` (always builds).
ENSEMBL_PATTERN = (
    "https://ftp.ensembl.org/pub/release-{release}/gtf/{species}/"
    "{cap_species}.{assembly}.{release}.gtf.gz"
)

# Targets the CLI exposes as subcommands (besides `list` and `ensembl`).
TARGETS = ("human", "mouse", "yeast", "ecoli", "human-v86")


# ------------------------------------------------------------------------ #
# Network + build helpers
# ------------------------------------------------------------------------ #

def _download(url: str, dest: Path) -> None:
    sys.stderr.write(f"downloading {url}\n  -> {dest}\n")
    dest.parent.mkdir(parents=True, exist_ok=True)
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
        [binary, "index", "--gtf", str(gtf), "--out", str(index)],
        stderr=None,  # surface the binary's progress directly
    )
    if proc.returncode != 0:
        sys.exit(f"error: `prot2exon index` exited {proc.returncode}")


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


# ------------------------------------------------------------------------ #
# Source resolution
# ------------------------------------------------------------------------ #

def _gtf_build_url(target: str, args: argparse.Namespace) -> tuple[str, str]:
    """Return (gtf_url, cache_stem) for a GTF-build of `target`."""
    if target == "ensembl":
        if not (args.species and args.assembly and args.release):
            sys.exit("error: ensembl requires --species, --assembly, --release "
                     "(e.g. --species danio_rerio --assembly GRCz11 --release 115)")
        cap = args.species[:1].upper() + args.species[1:]
        url = ENSEMBL_PATTERN.format(release=args.release, species=args.species.lower(),
                                     cap_species=cap, assembly=args.assembly)
        return url, f"{args.species}_v{args.release}"
    src = GTF_BUILD.get(target)
    if src is None:
        sys.exit(f"error: no GTF source known for {target!r}; "
                 f"pass --gtf-url to build from a custom GTF")
    if src.needs_release:
        if not args.release:
            sys.exit(f"error: building {target!r} from a GTF needs --release "
                     f"(e.g. {'49' if target == 'human' else 'M34'}), "
                     f"or pass --gtf-url with a GTF location")
        return src.url.format(release=args.release), f"{target}_v{args.release}"
    return src.url, target


def _do_idx_download(url: str, out: Path, sha256: str | None, force: bool) -> Path:
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.exists() and not force:
        sys.stderr.write(f"using cached index: {out}  (--force to re-download)\n")
        if sha256:
            _verify_sha256(out, sha256)
        print(out)
        return out
    _download(url, out)
    if sha256:
        _verify_sha256(out, sha256)
    sys.stderr.write(f"ready: {out}\n")
    print(out)
    return out


def _do_gtf_build(binary: str, gtf_url: str, stem: str, cache: Path,
                  out_idx: Path, force: bool, keep_gtf: bool) -> Path:
    cache.mkdir(parents=True, exist_ok=True)
    gtf_gz = cache / f"{stem}.gtf.gz"
    gtf = cache / f"{stem}.gtf"

    if out_idx.exists() and not force:
        sys.stderr.write(f"using cached index: {out_idx}  (--force to rebuild)\n")
        print(out_idx)
        return out_idx

    if not gtf_gz.exists() or force:
        _download(gtf_url, gtf_gz)
    else:
        sys.stderr.write(f"reusing cached download: {gtf_gz}\n")

    if not gtf.exists() or force:
        _gunzip(gtf_gz, gtf)
    else:
        sys.stderr.write(f"reusing cached gunzip:   {gtf}\n")

    _build_index(binary, gtf, out_idx)

    if not keep_gtf:
        try:
            gtf.unlink()
            sys.stderr.write(f"removed {gtf} (pass --keep-gtf to retain)\n")
        except FileNotFoundError:
            pass

    sys.stderr.write(f"ready: {out_idx}\n")
    print(out_idx)
    return out_idx


def _cmd_get(args: argparse.Namespace) -> Path:
    """Resolve a source for `args.target`, download or build, return the .idx path."""
    target = args.target
    binary = args.binary or _default_binary()
    cache = Path(args.cache_dir or _default_cache_dir())
    force = args.force

    # 1. Explicit GTF URL → build.
    if getattr(args, "gtf_url", None):
        if target == "ensembl" and getattr(args, "species", None):
            stem = f"{args.species}_v{args.release}"
        elif getattr(args, "release", None):
            stem = f"{target}_v{args.release}"
        else:
            stem = target
        out = Path(args.out) if args.out else (cache / f"{stem}.idx")
        return _do_gtf_build(binary, args.gtf_url, stem, cache, out,
                             force, args.keep_gtf)

    # 2. Explicit release, or the ensembl target → build from upstream.
    if getattr(args, "release", None) or target == "ensembl":
        gtf_url, stem = _gtf_build_url(target, args)
        out = Path(args.out) if args.out else (cache / f"{stem}.idx")
        return _do_gtf_build(binary, gtf_url, stem, cache, out,
                             force, args.keep_gtf)

    # 3. Default: pull the pre-built index from Zenodo (sha256 verified).
    zen = ZENODO_IDX.get(target)
    if zen is not None and "<RECORD>" not in zen.url:
        out = Path(args.out) if args.out else (cache / f"{target}.idx")
        return _do_idx_download(zen.url, out, zen.sha256, force)

    # 3b. Zenodo not published yet → fall back to a GTF build if we can.
    if target in GTF_BUILD or target == "ensembl":
        if zen is not None:
            sys.stderr.write(
                f"note: Zenodo deposit not published yet — building {target!r} "
                f"from the upstream GTF instead.\n")
        gtf_url, stem = _gtf_build_url(target, args)
        out = Path(args.out) if args.out else (cache / f"{stem}.idx")
        return _do_gtf_build(binary, gtf_url, stem, cache, out,
                             force, args.keep_gtf)

    # 3c. Zenodo-only target (e.g. human-v86) with an unpublished record.
    sys.exit(
        f"error: {target!r} is only distributed as a pre-built Zenodo index, "
        f"but the deposit isn't published yet (URL still has `<RECORD>`).\n"
        f"Build an equivalent from the upstream GTF with `--gtf-url`, or wait "
        f"for the published release."
    )


def _cmd_list() -> int:
    print("Targets default to a pre-built index from Zenodo. Override the source")
    print("with --release / --gtf-url to build from a GTF instead.")
    print()
    zenodo_ready = not any("<RECORD>" in z.url for z in ZENODO_IDX.values())
    if not zenodo_ready:
        print("  ⚠  The Zenodo record id isn't published yet — `<RECORD>` is still")
        print("     a placeholder. Until then, `fetch <target>` falls back to")
        print("     building from the upstream GTF where one is available.")
        print()
    print("=" * 72)
    print("Pre-built indexes  (the default `fetch <target>` source)")
    print("=" * 72)
    for name in sorted(ZENODO_IDX):
        print(f"  {name:<12s}  {ZENODO_IDX[name].short_desc}")
    print()
    print("=" * 72)
    print("Build from GTF  (--release / --gtf-url, or the default for targets")
    print("                 with no pre-built index such as ecoli)")
    print("=" * 72)
    for name in sorted(GTF_BUILD):
        print(f"  {name:<12s}  {GTF_BUILD[name].short_desc}")
    print()
    print("Generic Ensembl species (always built from the GTF):")
    print("    prot2exon fetch ensembl --species danio_rerio "
          "--assembly GRCz11 --release 115")
    print()
    print("Custom GTF:")
    print("    prot2exon fetch human --gtf-url https://your.host/x.gtf.gz")
    return 0


# ------------------------------------------------------------------------ #
# Python API
# ------------------------------------------------------------------------ #

def fetch_index(
    target: str,
    *,
    release: str | None = None,
    assembly: str | None = None,
    species: str | None = None,
    gtf_url: str | None = None,
    cache_dir: str | Path | None = None,
    out: str | Path | None = None,
    binary: str | None = None,
    force: bool = False,
    keep_gtf: bool = False,
    quiet: bool = False,
) -> Path:
    """Get a prot2exon index and return its path. Python mirror of ``fetch``.

    With no source override, downloads the pre-built index for ``target`` from
    Zenodo (falling back to a GTF build while the deposit is unpublished).
    Override the source with ``release`` or ``gtf_url`` to build from a GTF.

    Parameters
    ----------
    target : {"human", "mouse", "yeast", "ecoli", "human-v86", "ensembl", ...}
        A built-in target (see ``prot2exon fetch list``), or an Ensembl species
        name when paired with ``assembly`` + ``release``.
    release, assembly, species
        Select an upstream GTF build. ``release`` triggers a build for
        human/mouse; Ensembl species also need ``assembly``.
    gtf_url
        Build from this GTF URL instead of the default source.
    cache_dir, out, binary, force, keep_gtf, quiet
        As per the CLI flags.

    Examples
    --------
    >>> idx = fetch_index("human")                       # pre-built, from Zenodo
    >>> idx = fetch_index("human", release="49")         # build GENCODE v49
    >>> idx = fetch_index("yeast")                        # pre-built RefSeq R64
    >>> idx = fetch_index("danio_rerio", release="115", assembly="GRCz11")
    """
    # An unknown target paired with assembly/release is an Ensembl species.
    is_known = target in ZENODO_IDX or target in GTF_BUILD or target == "ensembl"
    if not is_known and (assembly or species):
        species = species or target
        eff_target = "ensembl"
    elif target == "ensembl":
        eff_target = "ensembl"
        species = species or None
    else:
        eff_target = target

    ns = argparse.Namespace(
        target=eff_target,
        species=species,
        assembly=assembly,
        release=release,
        gtf_url=gtf_url,
        cache_dir=str(cache_dir) if cache_dir else None,
        out=str(out) if out else None,
        binary=binary,
        force=force,
        keep_gtf=keep_gtf,
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
        prog="prot2exon fetch",
        description="Get a ready-to-query index: pre-built from Zenodo by "
                    "default, or built from a GTF you point it at.")
    sub = p.add_subparsers(dest="target", required=True)

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--out", help="Output .idx path (default: cache dir)")
    common.add_argument("--cache-dir",
                        help="Override cache dir (default: ~/.cache/prot2exon)")
    common.add_argument("--binary",
                        help="Path to the prot2exon binary (default: auto-discover)")
    common.add_argument("--gtf-url",
                        help="Build from this GTF instead of the default source")
    common.add_argument("--keep-gtf", action="store_true",
                        help="Keep the uncompressed GTF after a build")
    common.add_argument("--force", action="store_true",
                        help="Re-download / rebuild even if outputs exist")

    sub.add_parser("list", help="Show available targets and sources")

    for name in ("human", "mouse"):
        sp = sub.add_parser(name, parents=[common],
                            help=ZENODO_IDX[name].short_desc)
        sp.add_argument("--release",
                        help="Build from this GENCODE release instead of "
                             "the pre-built index (e.g. 49 / M34)")

    for name in ("yeast", "ecoli", "human-v86"):
        desc = (ZENODO_IDX[name].short_desc if name in ZENODO_IDX
                else GTF_BUILD[name].short_desc)
        sub.add_parser(name, parents=[common], help=desc)

    sp = sub.add_parser("ensembl", parents=[common],
                        help="Build any Ensembl species from its GTF")
    sp.add_argument("--species", required=True,
                    help="Ensembl species name (e.g. danio_rerio)")
    sp.add_argument("--assembly", required=True, help="Assembly tag (e.g. GRCz11)")
    sp.add_argument("--release", required=True, help="Ensembl release number")

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
