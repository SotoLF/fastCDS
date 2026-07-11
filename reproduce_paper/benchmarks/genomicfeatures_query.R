#!/usr/bin/env Rscript
# GenomicFeatures-NATIVE proteinToGenome via the GRangesList method (independent
# of ensembldb's EnsDb method). Maps aa ranges through an in-memory CDS-by-tx
# GRangesList. Emits the common classifier TSV: query_id chrom start end strand status
#
# Usage: Rscript gf_granges_query.R <cds_by_tx.rds> <queries_tsv> <out_tsv>
# queries_tsv (no header): transcript_id  aa_start  aa_end  query_id
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3)
suppressPackageStartupMessages({ library(GenomicFeatures); library(IRanges); library(GenomicRanges) })

cgl <- readRDS(args[1])
q <- read.table(args[2], header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                col.names = c("tx", "aa_start", "aa_end", "query_id"))
cat("loaded", nrow(q), "queries;", sum(q$tx %in% names(cgl)), "transcripts in index\n")

ir <- IRanges(start = q$aa_start, end = q$aa_end)
names(ir) <- q$tx

gf <- getMethod("proteinToGenome", "GRangesList")   # GenomicFeatures' own mapper
t0 <- Sys.time()
res <- suppressWarnings(gf(ir, cgl))
cat("gf GRangesList proteinToGenome done in", format(Sys.time() - t0), "\n")

con <- file(args[3], open = "w")
writeLines("query_id\tchrom\tstart\tend\tstrand\tstatus", con)
for (i in seq_along(res)) {
  qid <- q$query_id[i]
  gr <- res[[i]]
  if (is.null(gr) || length(gr) == 0) {
    writeLines(paste(qid, "NA", "NA", "NA", "NA", "no_result", sep = "\t"), con); next
  }
  chroms <- as.character(seqnames(gr)); starts <- start(gr); ends <- end(gr)
  strands <- as.character(strand(gr))
  for (j in seq_along(starts))
    writeLines(paste(qid, chroms[j], starts[j], ends[j], strands[j], "ok", sep = "\t"), con)
}
close(con)
cat("wrote", args[3], "\n")
