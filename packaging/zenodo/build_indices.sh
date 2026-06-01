#!/bin/bash
# Build the reproducible pre-built indices that should be uploaded to Zenodo
# for end-user `wget -O human.idx <zenodo-url>` consumption.
#
# Outputs:
#   <out-dir>/human.idx       (GENCODE v49 primary assembly)
#   <out-dir>/mouse.idx       (GENCODE M34 primary assembly)
#   <out-dir>/MANIFEST.tsv    (annotation source, sha256, idx size, build time)
#
# Usage:
#   packaging/zenodo/build_indices.sh [out-dir]

set -euo pipefail

OUT_DIR="${1:-${HOME}/Desktop/protein2genomic_data/zenodo_indices}"
mkdir -p "${OUT_DIR}"

# Locate the binary (prefer build/, then PATH).
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
BIN="${REPO_ROOT}/build/prot2exon"
if [[ ! -x "${BIN}" ]]; then
    BIN="$(command -v prot2exon-core || command -v prot2exon)"
fi
echo "binary: ${BIN}"

HUMAN_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.primary_assembly.annotation.gtf.gz"
MOUSE_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M34/gencode.vM34.primary_assembly.annotation.gtf.gz"

build_one() {
    local label="$1"; shift
    local url="$1";   shift
    local out_idx="$1"; shift
    local gz="${OUT_DIR}/$(basename "${url}")"
    local gtf="${gz%.gz}"

    if [[ ! -f "${gtf}" ]]; then
        echo "[${label}] downloading ${url}"
        curl -sSL -o "${gz}" "${url}"
        gunzip -kf "${gz}"
    fi
    local t0=$(date +%s)
    echo "[${label}] building index -> ${out_idx}"
    "${BIN}" index --gtf "${gtf}" --out "${out_idx}"
    local t1=$(date +%s)
    local size=$(stat -c%s "${out_idx}")
    local sha=$(sha256sum "${out_idx}" | cut -d' ' -f1)
    printf "%s\t%s\t%d\t%s\t%d\n" "${label}" "${url}" "${size}" "${sha}" "$((t1 - t0))" \
        >> "${OUT_DIR}/MANIFEST.tsv"
}

printf "label\tsource_url\tidx_bytes\tsha256\tbuild_seconds\n" > "${OUT_DIR}/MANIFEST.tsv"

build_one "human_gencode_v49" "${HUMAN_URL}" "${OUT_DIR}/human.idx"
build_one "mouse_gencode_M34" "${MOUSE_URL}" "${OUT_DIR}/mouse.idx"

echo
echo "Done. Manifest:"
cat "${OUT_DIR}/MANIFEST.tsv"
echo
echo "Next: upload ${OUT_DIR}/human.idx and ${OUT_DIR}/mouse.idx to Zenodo,"
echo "      then paste the resulting DOIs into README.md's quickstart."
