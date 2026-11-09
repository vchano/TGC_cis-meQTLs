#!/usr/bin/env Rscript
############################################################
# TreeGeneClimate (TGC) — TMS
# Step 5b1: methylKit import of one sample (all three contexts)
#
# Reads a sorted, deduplicated Bismark BAM with processBismarkAln and saves
# one methylRaw object per cytosine context (CpG, CHG, CHH). Import settings
# are the methylKit defaults: per-sample minimum coverage 10 reads and
# minimum base quality 20. No other filtering is applied here; the SNP mask,
# coverage limits and uniting are applied in Steps 5b2 and 5b3.
#
# USAGE
#   Rscript 5b1.tgc.tms.methylkit.import.R <sample_bam>
#   (output folder can be changed with the environment variable TGC_RAW_DIR)
#   (run as a SLURM array by 5b1.tgc.tms.methylkit.import.sh)
#
# INPUT
#   DATA/TMS/MAPPED.FILES.TMS/<SAMPLE>_R1_p_bismark_bt2_pe.dedup.bam  (Step 3b2)
#
# OUTPUT
#   RESULTS/TMS/METHYLKIT_RAW/<SAMPLE>_<CONTEXT>.rds
#   RESULTS/TMS/METHYLKIT_RAW/COUNTS/<SAMPLE>.tsv   (sites per context)
############################################################

suppressPackageStartupMessages(library(methylKit))

# === USER CONFIGURATION ===
PROJECT_ROOT <- Sys.getenv("TGC_PROJECT_ROOT", unset = "/path/to/your/project")
MIN_COV  <- 10L   # methylKit default
MIN_QUAL <- 20L   # methylKit default
ASSEMBLY <- "Picab02"
# ===========================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1 || !file.exists(args[1])) stop("Usage: 5b1...R <sample_bam>")
bam <- args[1]
sample_id <- sub("_R1_p_bismark_bt2_pe.*\\.bam$", "", basename(bam))

OUT_DIR <- Sys.getenv("TGC_RAW_DIR", unset = file.path(PROJECT_ROOT, "RESULTS/TMS/METHYLKIT_RAW"))
dir.create(file.path(OUT_DIR, "COUNTS"), recursive = TRUE, showWarnings = FALSE)

counts <- data.frame(sample = sample_id, context = c("CpG", "CHG", "CHH"), n_sites = NA_integer_)
for (i in seq_len(nrow(counts))) {
  ctx <- counts$context[i]
  mr <- processBismarkAln(location = bam, sample.id = sample_id, assembly = ASSEMBLY,
                          read.context = ctx, mincov = MIN_COV, minqual = MIN_QUAL,
                          save.folder = NULL, verbose = 0)
  saveRDS(mr, file.path(OUT_DIR, sprintf("%s_%s.rds", sample_id, ctx)))
  counts$n_sites[i] <- nrow(mr)
  rm(mr); gc()
}
write.table(counts, file.path(OUT_DIR, "COUNTS", paste0(sample_id, ".tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sample_id, ":", paste(counts$context, counts$n_sites, collapse = "; "), "\n")
cat("methylKit", as.character(packageVersion("methylKit")), "\n")
