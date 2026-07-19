# ensembldb route: proteinToGenome(IRanges, EnsDb) -> the EnsDb method, which
# hits the SQLite database once per query. This is the naive one-call route a
# user writes first, and it is what the "ensembldb" curve measures.
#
# Usage: Rscript run_ensembldb.R <ids.txt> [ensdb.sqlite]
#   ids.txt      one ENSP per line
#   ensdb.sqlite optional EnsDb SQLite; defaults to the EnsDb.Hsapiens.v86 package
suppressMessages({library(ensembldb); library(IRanges)})

a <- commandArgs(TRUE)
if (length(a) < 1) stop("run_ensembldb.R: need <ids.txt> [ensdb.sqlite]")
ids <- readLines(a[1])

if (length(a) >= 2 && nzchar(a[2]) && file.exists(a[2])) {
  edb <- EnsDb(a[2])
} else {
  suppressMessages(library(EnsDb.Hsapiens.v86))
  edb <- EnsDb.Hsapiens.v86
}

# Same 50-aa window for every query, so the comparison measures the mapping and
# not how much sequence each tool happens to be handed.
prng <- IRanges(start = rep(1, length(ids)), end = rep(50, length(ids)), names = ids)
res <- proteinToGenome(prng, edb)
cat(sprintf("ensembldb mapped=%d/%d\n", sum(lengths(res) > 0), length(ids)))
