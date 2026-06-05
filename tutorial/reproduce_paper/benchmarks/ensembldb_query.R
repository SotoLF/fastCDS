#!/usr/bin/env Rscript
# Batch-query ensembldb::proteinToGenome for a set of (protein_id, aa_start, aa_end, query_id)
# rows and emit a TSV of the resulting genomic intervals.
#
# Usage:
#   Rscript ensembldb_query.R <ensdb_sqlite> <queries_tsv> <out_tsv>
#
# queries_tsv: tab-separated, no header. Columns: protein_id  aa_start  aa_end  query_id
# out_tsv:     one row per (query_id, genomic interval). Columns:
#                  query_id  chrom  start  end  strand  status
#              status is one of: ok, no_result, error
#              For queries with no result, a single row with NA intervals + status is emitted.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("usage: ensembldb_query.R <ensdb_sqlite> <queries_tsv> <out_tsv>")
}
ensdb_path <- args[1]
queries_path <- args[2]
out_path <- args[3]

suppressPackageStartupMessages({
  library(ensembldb)
  library(IRanges)
})

edb <- EnsDb(ensdb_path)

q <- read.table(queries_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                col.names = c("protein_id", "aa_start", "aa_end", "query_id"))
cat("loaded", nrow(q), "queries\n")

# Build one IRanges per query, named with the protein_id (ensembldb uses names for the lookup).
ir <- IRanges(start = q$aa_start, end = q$aa_end)
names(ir) <- q$protein_id

cat("calling proteinToGenome ...\n")
t0 <- Sys.time()
res <- tryCatch(
  proteinToGenome(ir, edb, idType = "protein_id"),
  error = function(e) {
    cat("proteinToGenome error:", conditionMessage(e), "\n")
    NULL
  }
)
cat("proteinToGenome done in", format(Sys.time() - t0), "\n")

con <- file(out_path, open = "w")
writeLines("query_id\tchrom\tstart\tend\tstrand\tstatus", con)

if (is.null(res)) {
  # Whole batch failed — emit error rows so downstream classifies them.
  for (i in seq_len(nrow(q))) {
    writeLines(paste(q$query_id[i], "NA", "NA", "NA", "NA", "error", sep = "\t"), con)
  }
  close(con)
  quit(save = "no", status = 0)
}

for (i in seq_along(res)) {
  qid <- q$query_id[i]
  gr <- res[[i]]
  if (length(gr) == 0 || (is.list(gr) && all(sapply(gr, length) == 0))) {
    writeLines(paste(qid, "NA", "NA", "NA", "NA", "no_result", sep = "\t"), con)
    next
  }
  # gr may itself be a GRanges (rare) or have a structure with chrom/start/end accessors.
  chroms  <- as.character(seqnames(gr))
  starts  <- start(gr)
  ends    <- end(gr)
  strands <- as.character(strand(gr))
  for (j in seq_along(starts)) {
    writeLines(paste(qid, chroms[j], starts[j], ends[j], strands[j], "ok", sep = "\t"), con)
  }
}
close(con)
cat("wrote", out_path, "\n")
