# Installation

## Build from source

Requirements:

- C++17 toolchain (g++ ≥ 9, clang ≥ 10, or MSVC ≥ 2019)
- CMake ≥ 3.16
- OpenMP (optional — parallelises per-query processing; the binary still works without it, just single-threaded)

```bash
git clone https://github.com/SotoLF/Prot2Exon.git
cd Prot2Exon
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$(nproc)
```

You'll find the binary at `build/prot2exon`.

## Python wrapper

The Python package shells out to the C++ binary. Install with the wrapper requirements:

```bash
pip install -e python/         # editable install from the repo
# or, if you build a wheel: pip install dist/prot2exon-*.whl
```

The wrapper auto-discovers the binary via, in order:

1. `PROT2EXON_BIN` environment variable
2. `./build/prot2exon` relative to the repo root
3. `prot2exon` on `$PATH`

To pin a specific binary:

```bash
export PROT2EXON_BIN=/path/to/build/prot2exon
```

## Docker

A `Dockerfile` is shipped at the repo root. Build:

```bash
docker build -t prot2exon .
docker run --rm -v $(pwd):/work prot2exon \
    --index /work/human.idx --bed /work/queries.bed --out-dir /work/out --output all
```

The image bundles both the C++ binary and the Python wrapper.

## Smoke test

```bash
./build/prot2exon --version
python3 -m prot2exon.plot --help     # if the wrapper is installed
```

Then run the end-to-end test suite:

```bash
python3 tests/run_tests.py
```

You should see `109 passed, 2 failed` (the two failures need matplotlib + pandas in system python and are not regressions).

See [[Performance and RAM]] for tuning flags once you're up and running.
