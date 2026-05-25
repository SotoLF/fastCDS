# prot2exon reproducible build (Phase 5 deliverable).
#
# Build:
#   docker build -t prot2exon:dev .
# Run (mounted local data):
#   docker run --rm -v "$(pwd):/data" prot2exon:dev \
#       prot2exon --gtf /data/annotation.gtf --bed /data/queries.bed --out-dir /data/out
#
# CPU portability note: the upstream CMakeLists uses -O3 -march=native which
# pins the binary to the build machine's instruction set. For a redistributable
# image we override to -march=x86-64-v3 (Haswell+) so the resulting binary
# runs on any modern x86_64 host without rebuilds.

FROM ubuntu:24.04 AS builder

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        cmake \
        git \
        ca-certificates \
        libomp-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /src
COPY CMakeLists.txt ./
COPY src/ ./src/
COPY include/ ./include/

# Build with portable AVX2 baseline (override the repo's -march=native).
RUN cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_CXX_FLAGS_RELEASE="-O3 -DNDEBUG -march=x86-64-v3 -mtune=generic" \
    && cmake --build build --parallel \
    && ./build/prot2exon --help > /dev/null


FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
        bash \
        ca-certificates \
        libgomp1 \
        python3 \
        python3-pip \
        python3-venv \
    && rm -rf /var/lib/apt/lists/*

# Drop the C++ binary in /usr/local/bin as prot2exon-core; the bash wrapper
# (installed below) is what users invoke as `prot2exon`.
COPY --from=builder /src/build/prot2exon /usr/local/bin/prot2exon-core
COPY bin/prot2exon /usr/local/bin/prot2exon
RUN chmod +x /usr/local/bin/prot2exon /usr/local/bin/prot2exon-core

# Install the Python wrapper (provides `prot2exon plot` + `prot2exon` Python module).
RUN python3 -m venv /opt/prot2exon-venv
ENV PATH="/opt/prot2exon-venv/bin:${PATH}"
COPY README.md /opt/prot2exon/README.md
COPY python/ /opt/prot2exon/python/
RUN ln -sf /opt/prot2exon/README.md /opt/prot2exon/python/README.md \
    && sed -i 's|readme = "../README.md"|readme = "README.md"|' /opt/prot2exon/python/pyproject.toml \
    && pip install --no-cache-dir /opt/prot2exon/python

WORKDIR /data

# Smoke-test the image. Fails the build if either path is broken.
RUN prot2exon-core --help > /dev/null \
    && python3 -c "import prot2exon; print('prot2exon python OK', prot2exon.__version__)"

LABEL org.opencontainers.image.title="prot2exon" \
      org.opencontainers.image.description="Map protein domain coords to genomic/transcript structure." \
      org.opencontainers.image.licenses="MIT" \
      org.opencontainers.image.source="https://github.com/SotoLF/protein2genomic"

ENTRYPOINT ["prot2exon"]
CMD ["--help"]
