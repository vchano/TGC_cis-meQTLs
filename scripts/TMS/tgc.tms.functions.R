############################################################
# TreeGeneClimate (TGC) — TMS
# Shared functions for Steps 5b1-5b3 (sourced, not run directly)
############################################################

suppressPackageStartupMessages({
  library(methylKit)
  library(data.table)
})

# Read the SNP table written by Step 4a3 (chr, pos, ref, alt, qual, ac, an) and
# keep SNPs whose alternative allele(s) were called at least min_ac times
# (summed over alleles for multi-allelic sites). No call-rate filter: with about
# 7 ECS reads per tree many real SNPs are genotyped in only part of the trees,
# and such SNPs at a cytosine or its context bases still turn methylation
# differences into genotype differences.
read_snp_mask <- function(file, min_ac = 2L) {
  m <- fread(file, select = c("chr", "pos", "ac"), colClasses = list(character = "ac"))
  m[, ac_sum := vapply(strsplit(ac, ",", fixed = TRUE), function(x) sum(as.numeric(x)), numeric(1))]
  m <- m[ac_sum >= min_ac, .(chr, pos)]
  m <- unique(m)
  setkey(m, chr, pos)
  m
}

# Bases that define each cytosine context, counted from the cytosine on its own
# strand: CpG = C and the following G; CHG and CHH = C and the two following bases
MASK_WINDOW <- c(CpG = 1L, CHG = 2L, CHH = 2L)

# TRUE for sites with a SNP at the cytosine or at the bases defining its context
# on the same strand: [pos, pos + window] for "+" and [pos - window, pos] for "-".
# Sites with unknown strand ("*") are checked on both sides.
sites_hit_by_mask <- function(chr, pos, strand, mask, window = 2L) {
  n <- length(pos)
  hit <- logical(n)
  sgn <- ifelse(strand == "-", -1L, 1L)
  offsets <- 0:window
  for (o in offsets) {
    q <- data.table(chr = chr, pos = pos + sgn * o)
    hit <- hit | !is.na(mask[q, on = .(chr, pos), which = TRUE, mult = "first"])
    if (any(strand == "*")) {
      u <- which(strand == "*")
      q2 <- data.table(chr = chr[u], pos = pos[u] - o)
      hit[u] <- hit[u] | !is.na(mask[q2, on = .(chr, pos), which = TRUE, mult = "first"])
    }
  }
  hit
}

# Remove masked sites from a methylRaw object; returns list(object, n_removed)
mask_methylRaw <- function(mr, mask, context) {
  d <- getData(mr)
  hit <- sites_hit_by_mask(as.character(d$chr), as.integer(d$start), as.character(d$strand),
                           mask, window = MASK_WINDOW[[context]])
  list(object = mr[!hit, ], n_removed = sum(hit))
}

# Per-sample methylRaw files written by Step 5b1
raw_files <- function(raw_dir, cohort, context) {
  pattern <- if (cohort == "BREEDING") "^P00[1-3]_" else "^P00[4-8]_"
  f <- list.files(raw_dir, pattern = paste0(pattern, ".*_", context, "\\.rds$"), full.names = TRUE)
  f <- sort(f)
  names(f) <- sub(paste0("_", context, "\\.rds$"), "", basename(f))
  f
}

# Pooled methylated fraction per site across samples (sum numCs / sum coverage);
# Step 5b2 keeps sites whose pooled methylated fraction exceeds a threshold
# (methylated-fraction filter; output files carry the suffix _mef<threshold>)
pooled_methylated_fraction <- function(mb) {
  d <- getData(mb)
  cs <- rowSums(as.matrix(d[, grep("^numCs[0-9]+$", names(d))]), na.rm = TRUE)
  ts <- rowSums(as.matrix(d[, grep("^numTs[0-9]+$", names(d))]), na.rm = TRUE)
  ifelse(cs + ts > 0, cs / (cs + ts), NA_real_)
}
