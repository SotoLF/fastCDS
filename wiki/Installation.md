# How to install

prot2exon ships as a single package containing both a **C++ binary** (the `index` and `map` commands) and a **Python wrapper** (the `fetch` / `plot` commands, the `Mapper` API, and the DataFrame helpers). pip and bioconda both give you the whole thing.

| Method | Ships the binary | Best for |
|---|---|---|
| **pip** | ✅ pre-built wheels (Linux + macOS) | Most users. |
| **bioconda / pixi** | ✅ built by conda | conda/mamba users; pinned scientific stacks. |
| **build from source** | ✅ you compile it | Development, custom builds, Windows (via WSL). |

## pip

```bash
pip install prot2exon
# optional extras:
pip install "prot2exon[html]"        # plotly interactive HTML
pip install "prot2exon[all]"         # html + benchmarks + notebooks
```

The wheel bundles the compiled binary (as `prot2exon/_bin/prot2exon-core`), so all four commands work immediately — nothing else to install:

```bash
prot2exon --version
prot2exon fetch list
```

Pre-built `py3-none` wheels are published for **Linux (manylinux x86_64)** and **macOS (Intel + Apple Silicon)**, and work on any Python 3.9+. On platforms without a wheel — notably **Windows** — pip falls back to the source distribution, which compiles the C++ on install and therefore needs a C++17 toolchain + CMake (the smoothest Windows path is WSL or conda). To point the wrapper at a binary other than the bundled one, set `PROT2EXON_BIN`:

```bash
export PROT2EXON_BIN=/path/to/prot2exon-core
```

Discovery order: `$PROT2EXON_BIN`, then the wheel-bundled `_bin/prot2exon-core`, then `./build/prot2exon` in a source checkout, then `prot2exon-core` on `$PATH`.

## bioconda

```bash
conda install -c bioconda -c conda-forge prot2exon
# or, faster:
mamba install -c bioconda -c conda-forge prot2exon
```

conda compiles the binary as part of the recipe, so all four commands land on your `PATH` just like the pip install.

## pixi

[pixi](https://pixi.sh) resolves from the same conda channels:

```bash
pixi add -c bioconda -c conda-forge prot2exon      # in a project
pixi global install -c bioconda -c conda-forge prot2exon   # as a global tool
```

## Build from source

Requirements:

- C++17 toolchain (g++ ≥ 9, clang ≥ 10, or MSVC ≥ 2019)
- CMake ≥ 3.16
- OpenMP (optional — enables `--threads`; the binary still runs single-threaded without it)

```bash
git clone https://github.com/SotoLF/Prot2Exon.git
cd Prot2Exon
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$(nproc)
```

The binary lands at `build/prot2exon`. Then install the Python wrapper editable so `fetch`/`plot` and the API are available too:

```bash
cd ..
pip install -e .
```

The repo's `bin/prot2exon` wrapper finds `build/prot2exon` automatically, so you can run the four commands straight from the checkout.

## Docker

A `Dockerfile` at the repo root bundles both the binary and the wrapper:

```bash
docker build -t prot2exon .
docker run --rm -v "$(pwd):/work" prot2exon \
    map --index /work/human.idx --bed /work/queries.bed \
        --out-dir /work/out --output all
```

## Smoke test

```bash
prot2exon --version             # C++ binary, prints the index format version
prot2exon plot --help           # Python plotter
prot2exon fetch list            # pre-built indexes + GTF-build presets
```

Then run the end-to-end test suite from a source checkout:

```bash
python3 tests/run_tests.py      # expects "N passed, 0 failed"
```

Once installed, head to [[Building an index|Index]]. Tuning flags for large runs live on [[Performance and Benchmarking]].
