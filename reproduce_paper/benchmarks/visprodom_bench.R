#!/usr/bin/env Rscript
# Benchmark VisProDom's CreDat() — the package's batch domain->genome mapper
# (pure R/dplyr; maps InterProScan/RPS-BLAST domain aa-coords onto a GFF gene
# model). Reports wall time, mapped rows, and peak RSS on the bundled example
# (maize proteome), and a 100/1k/10k scaling slice.
#
# Usage: Rscript visprodom_bench.R <visprodom_repo_dir>
args <- commandArgs(trailingOnly = TRUE)
repo <- if (length(args) >= 1) args[1] else "/tmp/VisProDom"
suppressMessages({library(dplyr); library(data.table)})
source(file.path(repo, "R/CreDat.r"))
load(file.path(repo, "data/gff.rda")); load(file.path(repo, "data/annofile.rda"))
peak <- function() as.numeric(sub("VmHWM:\\s*([0-9]+).*","\\1",
                 grep("VmHWM", readLines("/proc/self/status"), value=TRUE)))/1024
qidx <- grep("^QUERY", annofile); hdr <- annofile[seq_len(qidx[1]-1)]
run <- function(af, label) {
  t0 <- Sys.time(); res <- suppressWarnings(suppressMessages(CreDat(gff, af))); t1 <- Sys.time()
  w <- as.numeric(difftime(t1,t0,units="secs"))
  cat(sprintf("VPD %-8s wall_s=%.2f mapped_rows=%d peak_rss_mb=%.0f\n",
              label, w, sum(!is.na(res$VVV4)), peak()))
}
for (K in c(100, 1000, 10000)) {
  cut <- if (K < length(qidx)) qidx[K+1]-1 else length(annofile)
  run(c(hdr, annofile[qidx[1]:cut]), paste0("N=", K))
}
run(annofile, "FULL")
