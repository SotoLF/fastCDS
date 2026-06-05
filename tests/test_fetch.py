"""`prot2exon fetch` — the pre-built-index downloader.

`fetch` downloads a sha256-verified binary index from Zenodo; it does not build
from a GTF (that is `prot2exon index`). Exercised offline against `file://`
URLs so nothing depends on a live host:
  * a target downloads to `--out` and verifies its baked-in sha256,
  * a sha256 mismatch aborts rather than keeping a corrupt index,
  * an unpublished (`<RECORD>`) target errors with a build-locally pointer,
  * `fetch list` advertises the pre-built targets.
"""

from __future__ import annotations

import hashlib
import os
import subprocess
import sys

import pytest

from conftest import REPO_ROOT


def _fake_target(monkeypatch, tmp_path, *, name="testtarget", sha=None):
    """Stage a `file://` Zenodo deposit with one fake index and register it."""
    from prot2exon import fetch as _fetch
    payload = b"synthetic index payload" * 100
    depot = tmp_path / "zenodo"
    depot.mkdir(exist_ok=True)
    (depot / f"{name}.idx").write_bytes(payload)
    real_sha = hashlib.sha256(payload).hexdigest()
    monkeypatch.setattr(_fetch, "_ZENODO_BASE", depot.resolve().as_uri())
    monkeypatch.setitem(_fetch.ZENODO_IDX, name,
                        _fetch.ZenodoIndex(f"{name}.idx", sha or real_sha, "test"))
    return _fetch, payload


def test_fetch_downloads_and_verifies(tmp_path, monkeypatch):
    """A target downloads to `--out` and the file matches the expected sha256."""
    _fetch, payload = _fake_target(monkeypatch, tmp_path)
    out = tmp_path / "got.idx"
    got = _fetch.fetch_index("testtarget", out=out, quiet=True)
    assert got == out and out.exists()
    assert out.read_bytes() == payload


def test_fetch_second_run_is_cached(tmp_path, monkeypatch):
    """Re-fetching reuses the cached file instead of re-downloading."""
    _fetch, _ = _fake_target(monkeypatch, tmp_path)
    out = tmp_path / "got.idx"
    _fetch.fetch_index("testtarget", out=out, quiet=True)
    import io
    buf = io.StringIO()
    old = sys.stderr
    sys.stderr = buf
    try:
        _fetch.fetch_index("testtarget", out=out)
    finally:
        sys.stderr = old
    assert "cached index" in buf.getvalue()


def test_fetch_bad_sha_aborts(tmp_path, monkeypatch):
    """A sha256 mismatch aborts with a clear message."""
    _fetch, _ = _fake_target(monkeypatch, tmp_path, sha="0" * 64)
    with pytest.raises(SystemExit) as exc:
        _fetch.fetch_index("testtarget", out=tmp_path / "x.idx", quiet=True)
    assert "sha256 mismatch" in str(exc.value)


def test_fetch_unpublished_target_points_to_index(tmp_path):
    """An unpublished (`<RECORD>`) target errors with a build-locally pointer."""
    proc = subprocess.run(
        [sys.executable, "-m", "prot2exon.fetch", "human-v86",
         "--out", str(tmp_path / "v86.idx")],
        env={**os.environ, "PYTHONPATH": str(REPO_ROOT / "python")},
        capture_output=True, text=True,
    )
    assert proc.returncode != 0
    assert "prot2exon index" in proc.stderr


def test_fetch_unknown_target_rejected(tmp_path):
    """An unknown target is rejected by argparse (only known subcommands exist)."""
    proc = subprocess.run(
        [sys.executable, "-m", "prot2exon.fetch", "nope"],
        env={**os.environ, "PYTHONPATH": str(REPO_ROOT / "python")},
        capture_output=True, text=True,
    )
    assert proc.returncode != 0


def test_fetch_list(tmp_path):
    """`fetch list` shows the pre-built indexes including a concrete target."""
    proc = subprocess.run(
        [sys.executable, "-m", "prot2exon.fetch", "list"],
        env={**os.environ, "PYTHONPATH": str(REPO_ROOT / "python")},
        capture_output=True, text=True,
    )
    assert proc.returncode == 0, proc.stderr
    assert "Pre-built indexes" in proc.stdout
    assert "yeast" in proc.stdout
    assert "prot2exon index" in proc.stdout   # points users to the build path
