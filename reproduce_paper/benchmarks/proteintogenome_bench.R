#!/usr/bin/env Rscript
# Benchmark a protein-to-genome mapper on a batch of (protein_id, aa_start,
# aa_end, query_id) queries, emitting per-segment genomic intervals plus a
# timing/RAM row. Run one tool per process so peak RSS is attributable.
#
# Usage:
#   Rscript proteintogenome_bench.R <tool> <ensdb_sqlite> <queries_tsv> \
#                                   <out_intervals_tsv> <timing_tsv>
#   tool in {ensembldb, genomicfeatures}
#
# ensembldb      : ensembldb::proteinToGenome(IRanges-by-protein, EnsDb) - queries
#                  the EnsDb SQLite per call.
# genomicfeatures: GenomicFeatures::proteinToGenome(IRanges-by-tx, GRangesList) -
#                  builds a CDS-by-transcript GRangesList once (setup), then maps
#                  entirely in memory (no SQLite during mapping).

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5)
  stop("usage: proteintogenome_bench.R <tool> <ensdb> <queries> <out> <timing>")
tool <- args[1]; ensdb_path <- args[2]; queries_path <- args[3]
out_path <- args[4]; timing_path <- args[5]

suppressPackageStartupMessages({
  library(ensembldb); library(GenomicFeatures)
  library(IRanges);   library(GenomicRanges)
})

peak_rss_mb <- function()
  as.numeric(sub("VmHWM:\\s*([0-9]+).*", "\\1",
                 grep("VmHWM", readLines("/proc/self/status"), value = TRUE))) / 1024

edb <- EnsDb(ensdb_path)
q <- read.table(queries_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                col.names = c("protein_id", "aa_start", "aa_end", "query_id"))
cat("tool:", tool, " loaded", nrow(q), "queries\n")

setup_s <- 0
write_intervals <- function(res, qids, drop_qids = character(0), drop_status = "no_result") {
  con <- file(out_path, "w")
  writeLines("query_id\tchrom\tstart\tend\tstrand\tstatus", con)
  for (i in seq_along(res)) {
    gr <- res[[i]]; qid <- qids[i]
    if (length(gr) == 0) { writeLines(paste(qid,"NA","NA","NA","NA","no_result",sep="\t"), con); next }
    ch <- as.character(seqnames(gr)); st <- start(gr); en <- end(gr); sd <- as.character(strand(gr))
    for (j in seq_along(st)) writeLines(paste(qid, ch[j], st[j], en[j], sd[j], "ok", sep="\t"), con)
  }
  for (qid in drop_qids) writeLines(paste(qid,"NA","NA","NA","NA",drop_status,sep="\t"), con)
  close(con)
}

if (tool == "ensembldb") {
  ir <- IRanges(start = q$aa_start, end = q$aa_end); names(ir) <- q$protein_id
  t0 <- Sys.time()
  res <- proteinToGenome(ir, edb, idType = "protein_id")
  map_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  write_intervals(res, q$query_id)

} else if (tool == "genomicfeatures") {
  # one-time setup: CDS-by-transcript GRangesList + protein_id -> tx_id map
  t0 <- Sys.time()
  cds_grl <- cdsBy(edb, by = "tx")
  pmap <- ensembldb::select(edb, keys = unique(q$protein_id),
                            keytype = "PROTEINID", columns = c("PROTEINID", "TXID"))
  setup_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  p2tx <- setNames(pmap$TXID, pmap$PROTEINID)
  q$tx_id <- p2tx[q$protein_id]
  keep <- !is.na(q$tx_id) & q$tx_id %in% names(cds_grl)
  ir <- IRanges(start = q$aa_start[keep], end = q$aa_end[keep]); names(ir) <- q$tx_id[keep]
  # ensembldb also defines a proteinToGenome method for GRangesList (which
  # requires protein-sequence metadata on the CDS); loading EnsDb pulls it in
  # and it can shadow GenomicFeatures'. Fetch GenomicFeatures' coordinate-only
  # method explicitly so we benchmark *its* implementation, as intended.
  gf_p2g <- getMethod("proteinToGenome", signature(db = "GRangesList"),
                      where = asNamespace("GenomicFeatures"))
  t1 <- Sys.time()
  res <- suppressWarnings(gf_p2g(ir, cds_grl))
  map_s <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
  write_intervals(res, q$query_id[keep], q$query_id[!keep], "no_protein_map")

} else stop("unknown tool: ", tool)

rss <- peak_rss_mb()
cat(sprintf("RESULT tool=%s setup_s=%.3f map_s=%.3f total_s=%.3f peak_rss_mb=%.0f\n",
            tool, setup_s, map_s, setup_s + map_s, rss))
writeLines(c("tool\tsetup_s\tmap_s\ttotal_s\tn_queries\tpeak_rss_mb",
             sprintf("%s\t%.3f\t%.3f\t%.3f\t%d\t%.0f",
                     tool, setup_s, map_s, setup_s + map_s, nrow(q), rss)),
           timing_path)
cat("wrote", out_path, "and", timing_path, "\n")
