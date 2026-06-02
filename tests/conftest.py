"""Shared pytest fixtures and helpers for the prot2exon test suite.

Most tests drive the C++ binary directly and golden-diff its outputs; others
cover the Python wrapper, the plot helpers, `prot2exon fetch`, and the
packaging scripts. The shared `out_all` fixture maps a fixed multi-query BED
once per session so the integration tests can read its outputs without each
re-running the mapper.

Pass `--update-goldens` to regenerate the expected files under tests/golden/.
"""

from __future__ import annotations

import filecmp
import shutil
import subprocess
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
BIN = REPO_ROOT / "build" / "prot2exon"
WRAPPER = REPO_ROOT / "bin" / "prot2exon"
TESTS_DIR = Path(__file__).resolve().parent
GOLDEN_DIR = TESTS_DIR / "golden"
WITH_TAGS_GTF = TESTS_DIR / "with_tags.gtf"
NO_TAGS_GTF = TESTS_DIR / "no_tags.gtf"

# A fixed multi-query BED exercising the mapper's main paths: ENSP vs ENST,
# MANE flags, structure-only, beyond-length, negative strand, codon splits,
# a selenoprotein-like CDS mismatch, a non-coding transcript, and a missing
# protein. The `out_all` fixture maps this once per session; several tests
# assert on specific input_ids (Q1_ENSP, Q7_SEC, Q8_SPLIT12, …) from it.
BED1 = """\
# 1) ENSP query, MANE+canonical, plus strand, aa 1..10 = CDS_1 fully
ENSP1\t1\t10\tQ1_ENSP
# 2) ENST query for the same transcript — should produce identical mapping
ENST1\t1\t10\tQ2_ENST
# 3) Versioned ENSP — version is stripped
ENSP1.5\t1\t10\tQ3_VER
# 4) Structure-only (no aa)
ENSP1\t0\t0\tQ4_STRUCT
ENSP2\t0\t0\tQ4b_STRUCT_M
# 5) Domain beyond protein length (ENSP1 has 44 aa)
ENSP1\t100\t110\tQ5_BEYOND
# 6) Negative strand, full coverage of CDS_1 in translation order (aa 1..17)
ENSP2\t1\t17\tQ6_NEG
# 7) Selenoprotein-like, domain at the end (aa 8..8)
ENSP3\t8\t8\tQ7_SEC
# 8) Codon split 1+2 — domain on the split aa (aa 2..2)
ENSP4\t2\t2\tQ8_SPLIT12
# 9) Codon split 2+1 — domain on the split aa (aa 2..2)
ENSP5\t2\t2\tQ9_SPLIT21
# 10) Non-coding ENST — no CDS
ENST6\t1\t10\tQ10_NONCODING
# 11) Protein not in index
ENSP99\t1\t10\tQ11_NOTFOUND
"""


def pytest_addoption(parser):
    parser.addoption(
        "--update-goldens",
        action="store_true",
        default=False,
        help="Overwrite golden files with the current binary output.",
    )


@pytest.fixture(scope="session")
def update_goldens(pytestconfig) -> bool:
    return pytestconfig.getoption("--update-goldens")


@pytest.fixture(scope="session")
def binary() -> Path:
    if not BIN.exists():
        pytest.fail(
            f"binary not found: {BIN}. Build it first:\n"
            f"  cd {REPO_ROOT} && cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j"
        )
    return BIN


@pytest.fixture(scope="session")
def synthetic_gtfs(tmp_path_factory):
    """Regenerate the synthetic GTFs from make_synthetic_gtf.py once per session."""
    subprocess.run(
        ["python3", str(TESTS_DIR / "make_synthetic_gtf.py")],
        check=True, capture_output=True,
    )
    return {"with_tags": WITH_TAGS_GTF, "no_tags": NO_TAGS_GTF}


@pytest.fixture(scope="session")
def with_tags_index(binary, synthetic_gtfs, tmp_path_factory) -> Path:
    idx_dir = tmp_path_factory.mktemp("idx")
    idx = idx_dir / "with_tags.idx"
    subprocess.run(
        [str(binary), "index", "--gtf", str(synthetic_gtfs["with_tags"]),
         "--out", str(idx)],
        check=True, capture_output=True, text=True,
    )
    return idx


@pytest.fixture(scope="session")
def no_tags_index(binary, synthetic_gtfs, tmp_path_factory) -> Path:
    """Index built from the GTF that carries no tag attributes at all — used to
    verify the MANE/canonical flags come back as `NA` (vs `false` when the GTF
    has tags but this transcript lacks them)."""
    idx_dir = tmp_path_factory.mktemp("idx_no_tags")
    idx = idx_dir / "no_tags.idx"
    subprocess.run(
        [str(binary), "index", "--gtf", str(synthetic_gtfs["no_tags"]),
         "--out", str(idx)],
        check=True, capture_output=True, text=True,
    )
    return idx


@pytest.fixture(scope="session")
def out_all(binary, with_tags_index, tmp_path_factory) -> Path:
    """Map the fixed BED1 query set once with `--output all`; return the output
    directory. Shared by the BED12, batch-size, and plotting tests."""
    out = tmp_path_factory.mktemp("out_all")
    bed = out / "queries.bed"
    bed.write_text(BED1)
    subprocess.run(
        [str(binary), "map", "--index", str(with_tags_index),
         "--bed", str(bed), "--out-dir", str(out), "--output", "all"],
        check=True, capture_output=True, text=True,
    )
    return out


def run_mapping(binary: Path, index: Path, bed_text: str, out_dir: Path,
                *, output_kind: str = "all") -> None:
    """Write `bed_text` to <out_dir>/queries.bed, then map into out_dir."""
    out_dir.mkdir(parents=True, exist_ok=True)
    bed = out_dir / "queries.bed"
    bed.write_text(bed_text)
    proc = subprocess.run(
        [str(binary), "map", "--index", str(index),
         "--bed", str(bed), "--out-dir", str(out_dir),
         "--output", output_kind],
        capture_output=True, text=True,
    )
    if proc.returncode != 0:
        raise AssertionError(
            f"prot2exon failed (rc={proc.returncode})\n"
            f"stdout: {proc.stdout}\nstderr: {proc.stderr}"
        )


def read_tsv(path: Path) -> list[dict]:
    """Read a TSV into a list of dict rows (header-keyed)."""
    import csv
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def summary_by_id(out_dir: Path) -> dict[str, dict]:
    """Read `domain_mapping_summary.tsv` from `out_dir`, keyed by input_id."""
    return {r["input_id"]: r for r in read_tsv(out_dir / "domain_mapping_summary.tsv")}


def assert_matches_golden(produced: Path, test_name: str,
                          filenames: list[str], *, update: bool) -> None:
    """Compare each file in `filenames` under `produced` against the golden
    copy under tests/golden/<test_name>/. With update=True, copy current
    outputs over the goldens and skip the diff."""
    golden = GOLDEN_DIR / test_name
    if update:
        golden.mkdir(parents=True, exist_ok=True)
        for name in filenames:
            src = produced / name
            assert src.exists(), f"binary did not produce {name}"
            shutil.copy2(src, golden / name)
        pytest.skip(f"goldens updated under {golden.relative_to(REPO_ROOT)}")

    if not golden.exists():
        pytest.fail(
            f"golden dir missing: {golden.relative_to(REPO_ROOT)}.\n"
            f"Run with --update-goldens to create it."
        )
    diffs = []
    for name in filenames:
        prod_path = produced / name
        gold_path = golden / name
        if not prod_path.exists():
            diffs.append(f"  {name}: not produced by binary")
            continue
        if not gold_path.exists():
            diffs.append(f"  {name}: missing golden (run --update-goldens)")
            continue
        if not filecmp.cmp(prod_path, gold_path, shallow=False):
            diffs.append(_format_diff(prod_path, gold_path, name))
    if diffs:
        raise AssertionError("golden mismatch:\n" + "\n".join(diffs))


def _format_diff(produced: Path, golden: Path, name: str) -> str:
    """Inline a short unified-diff snippet so failures are debuggable in CI."""
    import difflib
    prod = produced.read_text().splitlines(keepends=True)
    gold = golden.read_text().splitlines(keepends=True)
    diff = "".join(difflib.unified_diff(
        gold, prod,
        fromfile=f"golden/{name}", tofile=f"produced/{name}",
        n=2,
    ))
    snippet = "\n".join(diff.splitlines()[:30])
    return f"  {name}:\n{snippet}"
