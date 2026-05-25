# Zenodo — pre-built indices

PLAN.txt Phase 5 calls for downloadable, version-pinned indices on Zenodo so
end users can `wget` instead of rebuilding from a 3 GB GTF.

## Build the indices reproducibly

```bash
packaging/zenodo/build_indices.sh ~/zenodo_indices
```

This downloads GENCODE v49 (human) + M34 (mouse) primary-assembly GTFs from
EBI, builds two `.idx` files, and writes a `MANIFEST.tsv` with the source
URLs, file sizes, sha256 sums, and build time. Indices are reproducible bit
for bit given the same GTF (the binary uses a deterministic layout).

## Upload

1. Create a Zenodo deposit (or use an existing one — Zenodo supports versioned
   updates that keep a parent DOI).
2. Drag `human.idx` and `mouse.idx` in.
3. In the metadata: title `prot2exon pre-built indices — GENCODE v49 / M34`,
   keywords `prot2exon, GENCODE, Pfam, protein domain, genomic coordinates`,
   linked to the GitHub repo.
4. Publish; copy the two file URLs (each file gets its own download link).

## Wire the DOIs into README.md

Replace the placeholder URLs in the project's `README.md` quickstart with the
real Zenodo links:

```bash
# Human GENCODE v49 primary assembly (≈300 MB)
wget -O human.idx https://zenodo.org/records/<RECORD_ID>/files/human.idx

# Mouse GENCODE M34 primary assembly (≈250 MB)
wget -O mouse.idx https://zenodo.org/records/<RECORD_ID>/files/mouse.idx
```

## Why pre-built?

Building the human v49 index locally takes ~70 s and 6 GB of peak memory
(the raw 3.1 GB GTF expands during parse). A pre-built `.idx` removes that
cost for users who only need to run queries.
