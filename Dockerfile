# fastCDS reproducible build (Phase 5 deliverable).
#
# Build:
#   docker build -t fastCDS:dev .
# Run (mounted local data):
#   docker run --rm -v "$(pwd):/data" fastCDS:dev \
#       fastCDS map --gtf /data/annotation.gtf --bed /data/queries.bed --out-dir /data/out
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
    && ./build/fastCDS --help > /dev/null


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

# Drop the C++ binary in /usr/local/bin as fastCDS-core; the bash wrapper
# (installed below) is what users invoke as `fastCDS`.
COPY --from=builder /src/build/fastCDS /usr/local/bin/fastCDS-core
COPY bin/fastCDS /usr/local/bin/fastCDS
RUN chmod +x /usr/local/bin/fastCDS /usr/local/bin/fastCDS-core

# Install the Python wrapper (provides `fastCDS plot` + `fastCDS` Python module).
RUN python3 -m venv /opt/fastCDS-venv
ENV PATH="/opt/fastCDS-venv/bin:${PATH}"
COPY README.md /opt/fastCDS/README.md
COPY python/ /opt/fastCDS/python/
RUN ln -sf /opt/fastCDS/README.md /opt/fastCDS/python/README.md \
    && sed -i 's|readme = "../README.md"|readme = "README.md"|' /opt/fastCDS/python/pyproject.toml \
    && pip install --no-cache-dir /opt/fastCDS/python

WORKDIR /data

# Smoke-test the image. Fails the build if either path is broken.
RUN fastCDS-core --help > /dev/null \
    && python3 -c "import fastCDS; print('fastCDS python OK', fastCDS.__version__)"

LABEL org.opencontainers.image.title="fastCDS" \
      org.opencontainers.image.description="Map protein domain coords to genomic/transcript structure." \
      org.opencontainers.image.licenses="MIT" \
      org.opencontainers.image.source="https://github.com/SotoLF/fastCDS"

ENTRYPOINT ["fastCDS"]
CMD ["--help"]
