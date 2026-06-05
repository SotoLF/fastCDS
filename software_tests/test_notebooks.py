"""The notebook generator still runs.

Cheap regression: `tutorial/generate_notebooks.py` should import and emit the
expected .ipynb files. We pass `--out-dir` so the run writes to a tempdir
rather than overwriting any executed notebooks the user keeps under tutorial/.
"""

from __future__ import annotations

import subprocess
import sys

import pytest

from conftest import REPO_ROOT

NB_GEN = REPO_ROOT / "tutorial" / "generate_notebooks.py"


@pytest.mark.skipif(not NB_GEN.exists(), reason="notebook generator not present")
def test_generate_notebooks(tmp_path):
    out = tmp_path / "nb_out"
    proc = subprocess.run(
        [sys.executable, str(NB_GEN), "--out-dir", str(out)],
        capture_output=True, text=True,
    )
    assert proc.returncode == 0, proc.stderr
    assert "walkthrough_end_to_end.ipynb" in proc.stdout
    assert (out / "walkthrough_end_to_end.ipynb").exists()
    assert (out / "validation.ipynb").exists()
