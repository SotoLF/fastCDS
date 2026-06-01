"""Shared pytest fixtures and helpers for the prot2exon golden-file test suite.

The suite drives the C++ binary directly. The Python wrapper is out of scope
here — these tests verify the mapper, not the wrapper.

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
TESTS_DIR = Path(__file__).resolve().parent
GOLDEN_DIR = TESTS_DIR / "golden"
WITH_TAGS_GTF = TESTS_DIR / "with_tags.gtf"
NO_TAGS_GTF = TESTS_DIR / "no_tags.gtf"


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
