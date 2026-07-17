# GenomicFeatures route: the proteinToGenome,GRangesList method run on a
# CDS-by-transcript GRangesList precomputed once with cdsBy(). Mapping is then
# pure in-memory GRanges work, ~5x faster per query than the ensembldb route,
# with identical coordinates.
#
# Note: calling GenomicFeatures::proteinToGenome directly on an EnsDb dispatches
# back to the ensembldb method (same generic), so it measures the same thing as
# run_ensembldb.R. The speedup comes from the GRangesList input, not the
# namespace, which is why this runner takes the .rds and not an EnsDb.
#
# Build the .rds once:
#   library(ensembldb); library(EnsDb.Hsapiens.v86)
#   saveRDS(cdsBy(EnsDb.Hsapiens.v86, by = "tx"), "cds_by_tx.rds")
#
# Usage: Rscript run_gf_granges.R <enst_ids.txt> <cds_by_tx.rds>
suppressMessages({library(GenomicFeatures); library(IRanges)})

a <- commandArgs(TRUE)
if (length(a) < 2) stop("run_gf_granges.R: need <enst_ids.txt> <cds_by_tx.rds>")
ids <- readLines(a[1])
# The .rds load is inside the timed run on purpose: it is this route's index
# build, the counterpart of fastCDS loading its .idx.
cgl <- readRDS(a[2])

prng <- IRanges(start = rep(1, length(ids)), end = rep(50, length(ids)), names = ids)
gf <- getMethod("proteinToGenome", "GRangesList")
res <- suppressWarnings(gf(prng, cgl))
cat(sprintf("genomicfeatures mapped=%d/%d\n", sum(lengths(res) > 0), length(ids)))
