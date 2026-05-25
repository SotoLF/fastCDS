#!/usr/bin/env Rscript
# Build a custom EnsDb sqlite from a GENCODE/Ensembl GTF. Run once per annotation.
#
# Usage:
#   Rscript build_ensdb.R <gtf_path> <output_dir>
#
# The output_dir will contain a single .sqlite file that loadEnsDb() can open.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("usage: build_ensdb.R <gtf_path> <output_dir>")
}
gtf <- args[1]
outdir <- args[2]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(ensembldb)
})

cat("building EnsDb from", gtf, "...\n")
# GENCODE filenames don't match Ensembl's <organism>.<genome>.<version> scheme,
# so we pass the metadata explicitly. v49 corresponds to Ensembl 113 and
# GRCh38.p14 genome assembly.
db_path <- ensDbFromGtf(
  gtf = gtf,
  outfile = file.path(outdir, "ensdb.sqlite"),
  organism = "Homo_sapiens",
  genomeVersion = "GRCh38",
  version = 113
)
cat("wrote", db_path, "\n")

# Smoke-test: open and report metadata.
edb <- EnsDb(db_path)
print(edb)
