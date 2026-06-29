# How to install

fastCDS ships as a single package containing both a **C++ binary** (the `index` and `map` commands) and a **Python wrapper** (the `fetch` / `plot` commands, the `Mapper` API, and the DataFrame helpers). pip and bioconda both give you the whole thing.

## pip

```bash
pip install fastCDS
```

The wheel bundles the compiled binary (as `fastCDS/_bin/fastCDS-core`), so all four commands work immediately — nothing else to install:

```bash
fastCDS --version
fastCDS fetch list
```

Ready-to-use wheels cover **Linux** and **macOS** (Intel + Apple Silicon), Python 3.9+. **Windows is not supported directly** — there's no Windows build, so run fastCDS inside [WSL](https://learn.microsoft.com/windows/wsl/install) (Windows Subsystem for Linux), where it installs and behaves exactly like on Linux.

## bioconda

```bash
conda install -c bioconda -c conda-forge fastCDS
# or, faster:
mamba install -c bioconda -c conda-forge fastCDS
```

conda compiles the binary as part of the recipe, so all four commands land on your `PATH` just like the pip install.

## pixi

[pixi](https://pixi.sh) resolves from the same conda channels:

```bash
pixi add -c bioconda -c conda-forge fastCDS      # in a project
pixi global install -c bioconda -c conda-forge fastCDS   # as a global tool
```

## Build from source

Requirements:

- C++17 toolchain (g++ ≥ 9, clang ≥ 10, or MSVC ≥ 2019)
- CMake ≥ 3.16
- OpenMP (optional)

```bash
git clone https://github.com/SotoLF/fastCDS.git
cd fastCDS
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$(nproc)


pip install -e .
```

The repo's `bin/fastCDS` wrapper finds `build/fastCDS` automatically, so you can run the four commands straight from the checkout.

## Docker

A `Dockerfile` at the repo root builds an image with the binary and wrapper ready to go:

```bash
docker build -t fastCDS .
docker run --rm -v "$(pwd):/work" fastCDS \
    map --index /work/human.idx --bed /work/queries.bed --out-dir /work/out --output all
```
