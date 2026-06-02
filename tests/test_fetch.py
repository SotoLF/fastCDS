"""`prot2exon fetch` — the index retrieval/build helper.

Exercised offline against `file://` URLs so nothing depends on a live host:
  * `--gtf-url` downloads + gunzips + builds an index (and re-runs as a no-op),
  * the pre-built-index path verifies a baked-in sha256 (and aborts on mismatch),
  * `fetch list` advertises both pre-built targets and GTF-build sources.
"""

from __future__ import annotations

import gzip
import hashlib
import os
import subprocess
import sys

import pytest

from conftest import BIN, REPO_ROOT, WITH_TAGS_GTF, run_mapping, summary_by_id

FETCH_PY = REPO_ROOT / "python" / "prot2exon" / "fetch.py"


def _gzip_gtf(dest):
    with gzip.open(dest, "wb") as g:
        g.write(WITH_TAGS_GTF.read_bytes())
    return dest


def test_fetch_gtf_url_builds_working_index(binary, tmp_path):
    """`fetch <target> --gtf-url file://...` runs the download + gunzip + build
    pipeline, prints the index path, removes the uncompressed GTF (no
    --keep-gtf), and the index it builds actually maps queries."""
    fetch_url = _gzip_gtf(tmp_path / "fetch_src.gtf.gz").resolve().as_uri()
    cache = tmp_path / "cache"
    proc = subprocess.run(
        [sys.executable, str(FETCH_PY), "human", "--release", "TEST",
         "--gtf-url", fetch_url, "--cache-dir", str(cache), "--binary", str(binary)],
        capture_output=True, text=True,
    )
    assert proc.returncode == 0, proc.stderr
    idx_path = proc.stdout.strip().splitlines()[-1]
    assert idx_path and os.path.exists(idx_path)
    assert not (cache / "human_vTEST.gtf").exists()  # uncompressed GTF cleaned up

    # The freshly built index produces a working mapping.
    run_mapping(binary, idx_path, "ENSP1\t1\t10\tFETCH_TEST\n", tmp_path)
    assert summary_by_id(tmp_path)["FETCH_TEST"]["status"] == "ok"


def test_fetch_gtf_url_second_run_is_cached(binary, tmp_path):
    """Re-running the same fetch hits the cache instead of rebuilding."""
    fetch_url = _gzip_gtf(tmp_path / "fetch_src.gtf.gz").resolve().as_uri()
    cache = tmp_path / "cache"
    args = [sys.executable, str(FETCH_PY), "human", "--release", "TEST",
            "--gtf-url", fetch_url, "--cache-dir", str(cache), "--binary", str(binary)]
    first = subprocess.run(args, capture_output=True, text=True)
    assert first.returncode == 0, first.stderr
    second = subprocess.run(args, capture_output=True, text=True)
    assert second.returncode == 0, second.stderr
    assert "cached index" in second.stderr


def test_idx_download_verifies_sha256(tmp_path):
    """The pre-built-index path downloads to `--out` and verifies the file's
    sha256 against the expected hash."""
    from prot2exon import fetch as _fetch
    src = tmp_path / "src.idx"
    src.write_bytes(b"synthetic index payload" * 100)
    sha = hashlib.sha256(src.read_bytes()).hexdigest()
    url = src.resolve().as_uri()

    dest = tmp_path / "fetched.idx"
    got = _fetch._do_idx_download(url, dest, sha, force=True)
    assert os.path.exists(got)
    assert dest.stat().st_size == src.stat().st_size


def test_idx_download_bad_sha_aborts(tmp_path):
    """A sha256 mismatch aborts with a clear message rather than keeping a
    corrupt index."""
    from prot2exon import fetch as _fetch
    src = tmp_path / "src.idx"
    src.write_bytes(b"synthetic index payload" * 100)
    url = src.resolve().as_uri()

    with pytest.raises(SystemExit) as exc:
        _fetch._do_idx_download(url, tmp_path / "fetched.idx", "0" * 64, force=True)
    assert "sha256 mismatch" in str(exc.value)


def test_fetch_list(tmp_path):
    """`fetch list` shows the pre-built indexes, the GTF-build sources, and at
    least one concrete target (yeast)."""
    proc = subprocess.run(
        [sys.executable, "-m", "prot2exon.fetch", "list"],
        env={**os.environ, "PYTHONPATH": str(REPO_ROOT / "python")},
        capture_output=True, text=True,
    )
    assert proc.returncode == 0, proc.stderr
    assert "Pre-built indexes" in proc.stdout
    assert "Build from GTF" in proc.stdout
    assert "yeast" in proc.stdout


def test_fetch_zenodo_only_target_message(tmp_path):
    """A Zenodo-only target (human-v86, no GTF fallback) gives a meaningful
    message: the `<RECORD>` placeholder pre-publish, or a clean success once the
    deposit is live. Both are acceptable."""
    proc = subprocess.run(
        [sys.executable, "-m", "prot2exon.fetch", "human-v86",
         "--out", str(tmp_path / "v86.idx")],
        env={**os.environ, "PYTHONPATH": str(REPO_ROOT / "python")},
        capture_output=True, text=True,
    )
    assert "<RECORD>" in proc.stderr or proc.returncode == 0
