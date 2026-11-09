#!/usr/bin/env Rscript
############################################################
# TreeGeneClimate (TGC) — TMS
# Step 5b2: Population-level methylation datasets (one cohort x context)
#
# Used for methylation levels (6b), epigenomic PCA (7b), SVMP heatmaps (8b)
# and genomic-epigenomic concordance (11ab).
#
#   1) Remove sites with an ECS SNP (Step 4a3; alternative allele called at
#      least twice) at the cytosine or at the bases defining its context on the
#      same strand (CpG: C and G; CHG/CHH: C and the two following bases)
#   2) Per-sample coverage filter (methylKit filterByCoverage): at least 10
#      reads (as at import, Step 5b1) and at most the sample's 99.9th coverage
#      percentile, which removes pile-ups without penalising samples sequenced
#      more deeply (sequencing depth differs up to five-fold among trees)
#   3) unite(): keep sites covered in at least 80% (CpG), 70% (CHG) or
#      50% (CHH) of the cohort's samples (strand-specific, destrand = FALSE)
#   4) Methylated-fraction filter: keep sites whose pooled methylated fraction
#      across samples (sum numCs / sum coverage) is > 0.05
#
# USAGE
#   Rscript 5b2.tgc.tms.methylkit.population.datasets.R <BREEDING|NATURAL> <CpG|CHG|CHH>
#
# INPUTS
#   RESULTS/TMS/METHYLKIT_RAW/<SAMPLE>_<CONTEXT>.rds              (Step 5b1)
#   RESULTS/ECS/VARIANT.CALLING/MASK/tgc.ecs.cytosine.mask.snps.tsv.gz (Step 4a3)
#
# OUTPUTS (RESULTS/TMS/RANALYSIS/METHYLKIT_OBJECTS/)
#   methylBase_<cohort>_<ctx>_cov10_p99.9_mpg<k>.rds          (after unite)
#   methylBase_<cohort>_<ctx>_cov10_p99.9_mpg<k>_mef0.05.rds  (after step 4)
#   summary_<cohort>_<ctx>_population.tsv
############################################################

# === USER CONFIGURATION ===
PROJECT_ROOT <- Sys.getenv("TGC_PROJECT_ROOT", unset = "/path/to/your/project")
COV_LO <- 10L
COV_HI_PERC <- 99.9
MPG_FRAC <- c(CpG = 0.80, CHG = 0.70, CHH = 0.50)
METH_FRAC_MIN <- 0.05
# ===========================

source(file.path(PROJECT_ROOT, "SCRIPTS/TMS/tgc.tms.functions.R"))

args <- commandArgs(trailingOnly = TRUE)
COHORT <- toupper(args[1]); CTX <- args[2]
stopifnot(COHORT %in% c("BREEDING", "NATURAL"), CTX %in% names(MPG_FRAC))

RAW_DIR   <- file.path(PROJECT_ROOT, "RESULTS/TMS/METHYLKIT_RAW")
MASK_FILE <- file.path(PROJECT_ROOT, "RESULTS/ECS/VARIANT.CALLING/MASK/tgc.ecs.cytosine.mask.snps.tsv.gz")
OUT_DIR   <- file.path(PROJECT_ROOT, "RESULTS/TMS/RANALYSIS/METHYLKIT_OBJECTS")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
msg <- function(...) cat("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n", sep = "")

files <- raw_files(RAW_DIR, COHORT, CTX)
msg(COHORT, " / ", CTX, ": ", length(files), " samples")
mask <- read_snp_mask(MASK_FILE)
msg("Mask: ", nrow(mask), " SNP positions")

# 1) import + mask
n_raw <- n_masked <- integer(length(files))
objs <- vector("list", length(files))
for (i in seq_along(files)) {
  mr <- readRDS(files[i])
  n_raw[i] <- nrow(mr)
  mk <- mask_methylRaw(mr, mask, CTX)
  objs[[i]] <- mk$object
  n_masked[i] <- mk$n_removed
}
mrl <- methylRawList(objs, treatment = rep(0L, length(objs)))
rm(objs); gc()

# 2) coverage
mrl_cov <- filterByCoverage(mrl, lo.count = COV_LO, hi.perc = COV_HI_PERC)
n_cov <- vapply(as.list(mrl_cov), nrow, integer(1))
rm(mrl); gc()

# 3) unite
mpg <- as.integer(ceiling(MPG_FRAC[[CTX]] * length(files)))
mb <- unite(mrl_cov, destrand = FALSE, min.per.group = mpg)
rm(mrl_cov); gc()
msg("Sites after unite (min.per.group = ", mpg, "): ", nrow(mb))
stem <- sprintf("methylBase_%s_%s_cov%d_p%s_mpg%d", tolower(COHORT), tolower(CTX), COV_LO, COV_HI_PERC, mpg)
saveRDS(mb, file.path(OUT_DIR, paste0(stem, ".rds")))

# 4) pooled methylated fraction
pmf <- pooled_methylated_fraction(mb)
mb_mef <- mb[!is.na(pmf) & pmf > METH_FRAC_MIN, ]
msg("Sites after methylated-fraction filter (> ", METH_FRAC_MIN, "): ", nrow(mb_mef))
saveRDS(mb_mef, file.path(OUT_DIR, sprintf("%s_mef%.2f.rds", stem, METH_FRAC_MIN)))

summ <- data.frame(
  cohort = COHORT, context = CTX, n_samples = length(files),
  raw_sites_median = median(n_raw), masked_sites_median = median(n_masked),
  cov_sites_median = median(n_cov), min_per_group = mpg,
  n_sites_after_unite = nrow(mb), n_sites_after_mef = nrow(mb_mef)
)
write.table(summ, file.path(OUT_DIR, sprintf("summary_%s_%s_population.tsv", tolower(COHORT), tolower(CTX))),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(summ)
msg("Step 5b2 done")
