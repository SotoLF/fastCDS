# VisProDom (CreDat) on the human set (offline benchmark harness).
# Measures VisProDom speed/RAM mapping the human Ensembl-86 query set.
# Assumes the offline harness layout under /tmp/p2gbench/ (VisProDom checkout
# + visprodom_human.gff input); see benchmarks/README.md.
suppressMessages({library(dplyr); library(data.table)})
source("/tmp/p2gbench/VisProDom/R/CreDat.r")
gff <- as.data.frame(fread("visprodom_human.gff", header=FALSE, sep="\t", quote=""))
annofile <- readLines("visprodom_human.annofile")
peak <- function() as.numeric(sub("VmHWM:\\s*([0-9]+).*","\\1", grep("VmHWM", readLines("/proc/self/status"), value=TRUE)))/1024
cat("gff rows:", nrow(gff), " annofile lines:", length(annofile), "\n")
t0 <- Sys.time()
res <- suppressWarnings(suppressMessages(CreDat(gff, annofile)))
w <- as.numeric(difftime(Sys.time(), t0, units="secs"))
cat(sprintf("RESULT VPD_HUMAN wall_s=%.2f mapped_rows=%d peak_rss_mb=%.0f\n",
            w, sum(!is.na(res$VVV4)), peak()))
