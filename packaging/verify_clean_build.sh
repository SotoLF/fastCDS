#!/bin/bash
# Phase 5 clean-build verification.
#
# Performs the canonical user install from scratch in a throwaway directory
# (mirrors PLAN.txt Phase 5 line 187: "Verify git clone && cmake && make
# workflow on clean Ubuntu and macOS").
#
# Steps:
#   1. git clone into /tmp/<random>
#   2. cmake + make
#   3. ./build/prot2exon --help
#   4. pip install python/ into a throwaway venv
#   5. import prot2exon
#   6. pytest tests/
#
# Records build time + binary size + test count.

set -euo pipefail

# Source selection:
#   verify_clean_build.sh                 -> use this repo's working tree (copy via tar)
#   verify_clean_build.sh <url> [ref]     -> git clone the URL @ ref
SRC_MODE="${1:-local}"

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORK="$(mktemp -d -t prot2exon_verify_XXXXXX)"
trap 'rm -rf "${WORK}"' EXIT
echo "workdir: ${WORK}"

t_clone_0=$(date +%s)
if [[ "${SRC_MODE}" == "local" ]]; then
    echo "source: local working tree (${REPO_ROOT})"
    mkdir -p "${WORK}/src"
    # Use tar to honor .dockerignore-style excludes without depending on rsync.
    # NOTE: do not exclude tests/golden — they ARE the pytest fixtures.
    tar -C "${REPO_ROOT}" \
        --exclude='./build' --exclude='./.venv' --exclude='./.git' \
        --exclude='*.idx' --exclude='*.pdf' --exclude='*.gtf' --exclude='*.gtf.gz' \
        --exclude='*.fa' --exclude='*.fa.*' \
        --exclude='__pycache__' --exclude='*.egg-info' \
        --exclude='./python/dist' \
        -cf - . | tar -C "${WORK}/src" -xf -
else
    REMOTE="${SRC_MODE}"
    REF="${2:-main}"
    echo "source: git clone ${REMOTE} @ ${REF}"
    git clone --quiet --depth 1 --branch "${REF}" "${REMOTE}" "${WORK}/src"
fi
t_clone=$(( $(date +%s) - t_clone_0 ))
echo "[1/6] source staged (${t_clone}s)"

cd "${WORK}/src"

t_cmake_0=$(date +%s)
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release > "${WORK}/cmake.log" 2>&1
t_cmake=$(( $(date +%s) - t_cmake_0 ))
echo "[2/6] cmake OK (${t_cmake}s)"

t_make_0=$(date +%s)
cmake --build build --parallel >> "${WORK}/cmake.log" 2>&1
t_make=$(( $(date +%s) - t_make_0 ))
binary_size=$(stat -c%s build/prot2exon)
echo "[3/6] make OK (${t_make}s; binary=${binary_size} bytes)"

./build/prot2exon --help > /dev/null
./build/prot2exon --version
echo "[4/6] binary smoke test OK"

# Python wrapper in throwaway venv. Full install (with pandas + matplotlib)
# so the import smoke test exercises the real wrapper. If you want to skip
# the dep download, point PROT2EXON_VERIFY_VENV at an existing venv that
# already has pandas/matplotlib and we'll install --no-deps into a copy.
t_pip_0=$(date +%s)
if [[ -n "${PROT2EXON_VERIFY_VENV:-}" && -x "${PROT2EXON_VERIFY_VENV}/bin/python" ]]; then
    python3 -m venv --system-site-packages --copies "${WORK}/venv"
    # Bridge to the user-provided venv by appending its site-packages to PYTHONPATH.
    SP=$("${PROT2EXON_VERIFY_VENV}/bin/python" -c "import site; print(site.getsitepackages()[0])")
    export PYTHONPATH="${SP}${PYTHONPATH:+:${PYTHONPATH}}"
    echo "  using deps from ${PROT2EXON_VERIFY_VENV}"
    PIP_FLAGS="--no-deps"
else
    python3 -m venv "${WORK}/venv"
    PIP_FLAGS=""
fi
# shellcheck source=/dev/null
source "${WORK}/venv/bin/activate"
if ! pip install ${PIP_FLAGS} --quiet python/ > "${WORK}/pip.log" 2>&1; then
    echo "pip install FAILED — log:" >&2
    cat "${WORK}/pip.log" >&2
    exit 1
fi
t_pip=$(( $(date +%s) - t_pip_0 ))
python -c "import prot2exon; print('  wrapper version:', prot2exon.__version__)"
echo "[5/6] python wrapper install + import OK (${t_pip}s)"

# Run the project's pytest if available via `python -m pytest`.
if [[ -d tests ]] && python -c "import pytest" 2>/dev/null; then
    t_test_0=$(date +%s)
    # The tests need the prot2exon binary on PATH (the wrapper auto-discovers).
    export PROT2EXON_BIN="${WORK}/src/build/prot2exon"
    if python -m pytest -q tests/ > "${WORK}/pytest.log" 2>&1; then
        n_tests=$(grep -oE '[0-9]+ passed' "${WORK}/pytest.log" | head -1 || echo "?")
        t_test=$(( $(date +%s) - t_test_0 ))
        echo "[6/6] pytest: ${n_tests} (${t_test}s)"
    else
        echo "[6/6] pytest: FAILED — see ${WORK}/pytest.log" >&2
        tail -20 "${WORK}/pytest.log" >&2
        exit 1
    fi
else
    echo "[6/6] pytest skipped (pytest not available; would have run tests/)"
fi

echo
echo "=== VERIFICATION SUMMARY ==="
echo "  source:       ${SRC_MODE}"
echo "  stage:        ${t_clone}s"
echo "  cmake:        ${t_cmake}s"
echo "  make:         ${t_make}s"
echo "  pip install:  ${t_pip}s"
[[ -n "${t_test:-}" ]] && echo "  pytest:       ${t_test}s"
echo "  binary size:  ${binary_size} bytes"
echo "  total wall:   $((t_clone + t_cmake + t_make + t_pip + ${t_test:-0}))s"
echo "  workdir (will be cleaned): ${WORK}"
