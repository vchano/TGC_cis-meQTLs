#!/usr/bin/env Rscript
############################################################
# TreeGeneClimate (TGC) — JOINT ECS + TBS
# Step 12ab1: MatrixEQTL5 cis-meQTL mapping
#
# KEY DIFFERENCES vs 15ab2 (MatrixEQTL4):
#   - Inputs from 12ab0 (non-imputed GDS, GCTA GRM PCs, unfiltered methylation)
#   - PC1..PC10 for BOTH cohorts (vs 10/3 in meQTL4)
#   - No quantile-normalised M-values
#   - Cis window: 100 kb (1e5) — same as original ECS pipeline
#   - HWE group imputation of missing genotypes
#   - Lambda (genomic inflation) reported per panel in summary
#
# MODEL
#   M-value ~ SNP + PC1..PC10   (linear, no GRM random effect)
#
# USAGE
# SLURM single panel:
#   TGC_COHORT=BREEDING TGC_CONTEXT=CpG Rscript --vanilla 12ab1.R
# Interactive (all 6 panels):
#   source this script in RStudio (no env vars set)
#
# INPUTS (from 12ab0 — RESULTS/JOINT/MQTL5/INPUTS/<COHORT>/<CONTEXT>/)
#   methylation_mvalues_matrix.rds
#   methylation_site_annot.rds
#   shared_sample_ids.rds
#   pcs_shared.rds
#   sample_groups.rds
#   snp_variant_annot.rds
#   snp_gds_sample_ids.rds
#
# OUTPUTS (RESULTS/JOINT/MATRIXEQTL5/<COHORT>/<CONTEXT>/)
#   cis_meqtl_all_results.rds / .tsv.gz
#   cis_meqtl_significant.tsv        (FDR < 0.05)
#   cis_meqtl_top_per_site.tsv
#   cis_meqtl_panel_summary.tsv      (includes lambda)
#   qq_<cohort>_<context>.tiff
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(MatrixEQTL)
  library(SNPRelate)
  library(gdsfmt)
})

options(stringsAsFactors = FALSE)

############################################################
# 1) SETTINGS
############################################################

CIS_WINDOW   <- 1e5       # 100 kb — same as original ECS pipeline
MIN_SAMPLES  <- 10L
FDR_THRESH   <- 0.05
# Keep all cis pairs (p <= 1) so that genome-wide BH correction is applied
# post-hoc over the complete set rather than a pre-filtered subset
PVAL_CIS_OUT <- 1

############################################################
# 2) PATHS
############################################################

# === USER CONFIGURATION ===
# Set TGC_PROJECT_ROOT as an environment variable, or edit the fallback path below
PROJECT_ROOT <- Sys.getenv("TGC_PROJECT_ROOT",
  unset = "/path/to/your/project")
# ===========================

MQTL5_ROOT  <- file.path(PROJECT_ROOT, "RESULTS", "JOINT", "MQTL5")
INPUTDIR    <- file.path(MQTL5_ROOT, "INPUTS")

MEQTL_ROOT  <- file.path(PROJECT_ROOT, "RESULTS", "JOINT", "MATRIXEQTL5")
OUTROOT     <- MEQTL_ROOT
LOGDIR      <- file.path(MEQTL_ROOT, "LOGS")
SUMDIR      <- file.path(MEQTL_ROOT, "SUMMARIES")

for (d in c(OUTROOT, LOGDIR, SUMDIR))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Non-imputed GDS files (same source as used in 12ab0)
RDATA_DIR <- "/mnt/vast-standard/home/chano/u15584/treegeneclimate/2025/ECS/RDATA"
GDS_FILES <- list(
  BREEDING = file.path(RDATA_DIR, "breeding.snp.gds"),
  NATURAL  = file.path(RDATA_DIR, "natural.snp.gds")
)

############################################################
# 3) PHASE DETECTION
############################################################

# When run via SLURM with env vars set, process a single panel.
# When run interactively (no env vars), iterate over all 6 panels.
env_cohort  <- Sys.getenv("TGC_COHORT",  unset = "")
env_context <- Sys.getenv("TGC_CONTEXT", unset = "")
RUN_MAPPING <- nzchar(env_cohort) && nzchar(env_context)
SLURM_MODE  <- RUN_MAPPING

if (RUN_MAPPING) {
  stopifnot(env_cohort  %in% c("BREEDING", "NATURAL"))
  stopifnot(env_context %in% c("CpG", "CHG", "CHH"))
  PANELS <- data.frame(cohort = env_cohort, context = env_context,
                       stringsAsFactors = FALSE)
} else {
  PANELS <- expand.grid(cohort  = c("BREEDING", "NATURAL"),
                        context = c("CpG", "CHG", "CHH"),
                        stringsAsFactors = FALSE)
}

############################################################
# 4) HELPERS
############################################################

CURRENT_LOGFILE <- NULL

log_msg <- function(...) {
  txt <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ",
                paste0(..., collapse = ""))
  cat(txt, "\n")
  if (!is.null(CURRENT_LOGFILE))
    cat(txt, "\n", file = CURRENT_LOGFILE, append = TRUE)
}
sep_line <- function() log_msg(paste(rep("=", 70), collapse = ""))

normalize_id <- function(x) {
  x <- as.character(x); x <- trimws(x)
  x <- gsub("\\.0$", "", x); x <- gsub("^X", "", x); x <- gsub("-", "_", x); x
}

# HWE group imputation: fill missing genotypes with 2*p within each group
# (family for BREEDING, population for NATURAL). Falls back to cohort-wide
# allele frequency if a group has no coverage for a given SNP.
impute_hwe <- function(geno_snp_x_samp, groups) {
  mode(geno_snp_x_samp) <- "numeric"
  geno_snp_x_samp[geno_snp_x_samp > 2 | geno_snp_x_samp < 0] <- NA_real_
  # Cohort-wide allele frequency (fallback for groups with no data)
  p_all <- rowSums(geno_snp_x_samp, na.rm = TRUE) /
           (2 * rowSums(!is.na(geno_snp_x_samp)))
  # Remove SNPs with no valid genotype calls
  keep  <- is.finite(p_all)
  g     <- geno_snp_x_samp[keep, , drop = FALSE]
  p_all <- p_all[keep]
  for (grp in unique(groups)) {
    idx <- which(groups == grp)
    # Within-group allele frequency
    pg  <- rowSums(g[, idx, drop = FALSE], na.rm = TRUE) /
           (2 * rowSums(!is.na(g[, idx, drop = FALSE])))
    pg[!is.finite(pg)] <- p_all[!is.finite(pg)]  # fallback
    miss <- which(is.na(g[, idx, drop = FALSE]), arr.ind = TRUE)
    if (nrow(miss)) g[, idx][miss] <- 2 * pg[miss[, 1]]
  }
  # Final sweep: any remaining NAs filled with cohort-wide expectation
  if (anyNA(g)) {
    fill2 <- matrix(2 * p_all, nrow = nrow(g), ncol = ncol(g))
    g[is.na(g)] <- fill2[is.na(g)]
  }
  list(geno = g, keep = keep)
}

# Genomic inflation factor: ratio of observed to expected median chi-squared statistic
compute_lambda <- function(pvals) {
  pvals <- pvals[is.finite(pvals) & pvals > 0 & pvals <= 1]
  if (length(pvals) < 2) return(NA_real_)
  median(qchisq(1 - pvals, df = 1), na.rm = TRUE) / qchisq(0.5, df = 1)
}

############################################################
# 5) MAPPING LOOP
############################################################

panel_summaries <- list()

for (p in seq_len(nrow(PANELS))) {
  cohort  <- PANELS$cohort[p]
  context <- PANELS$context[p]

  panel_tag    <- paste0(tolower(cohort), "_", tolower(context))
  panel_input  <- file.path(INPUTDIR, cohort, context)
  panel_output <- file.path(OUTROOT, cohort, context)
  dir.create(panel_output, recursive = TRUE, showWarnings = FALSE)

  # Skip completed panels to allow safe re-runs / partial SLURM restarts
  if (file.exists(file.path(panel_output, "cis_meqtl_all_results.rds"))) {
    cat("[", format(Sys.time(), "%H:%M:%S"), "] Panel ", cohort, "/", context,
        " already done — skipping\n", sep = "")
    next
  }

  CURRENT_LOGFILE <- file.path(LOGDIR, paste0("step12ab1_", panel_tag, ".log"))
  SUMMARY_FILE    <- file.path(SUMDIR, paste0("step12ab1_", panel_tag, "_summary.tsv"))
  SESSION_FILE    <- file.path(SUMDIR, paste0("step12ab1_", panel_tag, "_sessionInfo.txt"))
  if (file.exists(CURRENT_LOGFILE)) file.remove(CURRENT_LOGFILE)

  sep_line()
  log_msg("Step 12ab1 — MatrixEQTL5 cis-meQTL mapping")
  log_msg("Panel: ", cohort, " / ", context)
  log_msg("Cis window: +/-", CIS_WINDOW / 1e3, " kb | Model: M ~ SNP + PC1..10")

  # --- Check inputs ---
  req <- c("methylation_mvalues_matrix.rds", "methylation_site_annot.rds",
           "shared_sample_ids.rds", "pcs_shared.rds", "sample_groups.rds",
           "snp_variant_annot.rds")
  missing_f <- req[!file.exists(file.path(panel_input, req))]
  if (length(missing_f)) {
    log_msg("Missing inputs: ", paste(missing_f, collapse = ", "))
    fwrite(data.table(cohort = cohort, context = context,
                      status = "missing_inputs", lambda = NA_real_),
           SUMMARY_FILE, sep = "\t")
    if (SLURM_MODE) quit(save = "no", status = 1L); next
  }
  gds_path <- GDS_FILES[[cohort]]
  if (!file.exists(gds_path)) {
    log_msg("Missing GDS: ", gds_path)
    fwrite(data.table(cohort = cohort, context = context,
                      status = "missing_gds", lambda = NA_real_),
           SUMMARY_FILE, sep = "\t")
    if (SLURM_MODE) quit(save = "no", status = 1L); next
  }

  # --- Load panel inputs ---
  log_msg("Loading panel inputs...")
  m_mat      <- readRDS(file.path(panel_input, "methylation_mvalues_matrix.rds"))
  site_annot <- as.data.table(readRDS(file.path(panel_input, "methylation_site_annot.rds")))
  shared_dt  <- readRDS(file.path(panel_input, "shared_sample_ids.rds"))
  pcs_dt     <- as.data.table(readRDS(file.path(panel_input, "pcs_shared.rds")))
  groups_dt  <- as.data.table(readRDS(file.path(panel_input, "sample_groups.rds")))
  var_annot  <- as.data.table(readRDS(file.path(panel_input, "snp_variant_annot.rds")))

  shared_ids <- shared_dt$sample_id
  n_shared   <- length(shared_ids)
  n_sites    <- ncol(m_mat)
  n_snps     <- nrow(var_annot)
  log_msg("Shared samples: ", n_shared, " | Sites: ", n_sites, " | SNPs: ", n_snps)

  if (n_shared < MIN_SAMPLES) {
    fwrite(data.table(cohort = cohort, context = context,
                      status = "too_few_samples", lambda = NA_real_,
                      n_samples = n_shared),
           SUMMARY_FILE, sep = "\t")
    if (SLURM_MODE) quit(save = "no", status = 1L); next
  }

  # --- Read genotypes from GDS ---
  log_msg("Reading genotypes from GDS...")
  gds           <- snpgdsOpen(gds_path, readonly = TRUE)
  gds_samp_ids  <- read.gdsn(index.gdsn(gds, "sample.id"))
  gds_samp_norm <- normalize_id(gds_samp_ids)
  # Extract only the raw IDs that map to shared normalised IDs
  raw_shared    <- gds_samp_ids[gds_samp_norm %in% shared_ids]
  geno_list     <- snpgdsGetGeno(gds, sample.id = raw_shared,
                                 snpfirstdim = TRUE, with.id = TRUE)
  snpgdsClose(gds)

  # Re-order columns to match shared_ids ordering for alignment with M-value matrix
  samp_reorder <- match(shared_ids, normalize_id(geno_list$sample.id))
  geno_mat     <- geno_list$genotype[, samp_reorder, drop = FALSE]  # SNPs x samples
  rm(geno_list, samp_reorder); gc()
  log_msg("Genotype matrix: ", nrow(geno_mat), " SNPs x ", ncol(geno_mat), " samples")

  # --- HWE group imputation ---
  log_msg("Imputing missing genotypes by group...")
  groups_ord <- groups_dt$group[match(shared_ids, groups_dt$sample_id)]
  groups_ord[is.na(groups_ord)] <- "Unknown"
  imp      <- impute_hwe(geno_mat, groups_ord)
  geno_imp <- imp$geno
  var_filt <- var_annot[imp$keep]  # variant annotation subset to polymorphic SNPs
  rm(geno_mat, imp); gc()

  # Drop zero-variance SNPs (monomorphic after imputation)
  snp_sd   <- apply(geno_imp, 1, sd, na.rm = TRUE)
  keep_var <- is.finite(snp_sd) & snp_sd > 0
  geno_imp <- geno_imp[keep_var, , drop = FALSE]
  var_filt <- var_filt[keep_var]
  log_msg("SNPs after imputation + variance filter: ", nrow(geno_imp))

  # --- Build MatrixEQTL SlicedData objects ---
  log_msg("Building MatrixEQTL SlicedData objects...")
  m_reorder  <- match(shared_ids, normalize_id(rownames(m_mat)))
  m_ord      <- m_mat[m_reorder, , drop = FALSE]

  rownames(geno_imp) <- as.character(var_filt$snp_id)
  colnames(geno_imp) <- shared_ids
  rownames(m_ord)    <- shared_ids

  # MatrixEQTL expects features x samples orientation
  snpsSD <- SlicedData$new(); snpsSD$CreateFromMatrix(geno_imp)
  geneSD <- SlicedData$new(); geneSD$CreateFromMatrix(t(m_ord))   # sites x samples

  pc_cols <- grep("^PC[0-9]+$", names(pcs_dt), value = TRUE)
  pcs_ord <- pcs_dt[match(shared_ids, pcs_dt$sample_id), ..pc_cols]
  cvrtSD  <- SlicedData$new()
  cvrtSD$CreateFromMatrix(t(as.matrix(pcs_ord)))
  log_msg("Covariates: ", paste(pc_cols, collapse = ", "))

  # Position tables required by MatrixEQTL to define the cis window
  snpspos <- data.frame(
    snp = as.character(var_filt$snp_id),
    chr = as.character(var_filt$chr),
    pos = as.integer(var_filt$pos),
    stringsAsFactors = FALSE
  )
  genepos <- data.frame(
    gene  = colnames(m_ord),
    chr   = as.character(site_annot$chr),
    start = as.integer(site_annot$start),
    end   = as.integer(site_annot$end),
    stringsAsFactors = FALSE
  )

  # --- Run MatrixEQTL ---
  # pvOutputThreshold = 0 suppresses trans output entirely (cis-only analysis)
  cis_out_file <- file.path(panel_output,
                             paste0(panel_tag, "_MatrixEQTL_cis.txt"))
  log_msg("Running Matrix_eQTL_main (cis, window=", CIS_WINDOW, " bp)...")
  me <- tryCatch(
    Matrix_eQTL_main(
      snps                  = snpsSD,
      gene                  = geneSD,
      cvrt                  = cvrtSD,
      output_file_name      = "",
      pvOutputThreshold     = 0,
      useModel              = modelLINEAR,
      errorCovariance       = numeric(),
      verbose               = TRUE,
      output_file_name.cis  = cis_out_file,
      pvOutputThreshold.cis = PVAL_CIS_OUT,
      snpspos               = snpspos,
      genepos               = genepos,
      cisDist               = CIS_WINDOW
    ),
    error = function(e) { log_msg("MatrixEQTL error: ", conditionMessage(e)); NULL }
  )

  if (is.null(me)) {
    fwrite(data.table(cohort = cohort, context = context,
                      status = "matrixeqtl_error", lambda = NA_real_),
           SUMMARY_FILE, sep = "\t")
    if (SLURM_MODE) quit(save = "no", status = 1L); next
  }

  cis_res <- me$cis$eqtls
  if (is.null(cis_res) || nrow(cis_res) == 0L) {
    log_msg("No cis results for ", cohort, " / ", context)
    fwrite(data.table(cohort = cohort, context = context,
                      status = "no_cis_results", lambda = NA_real_,
                      n_samples = n_shared, n_sites = n_sites,
                      n_snps = nrow(geno_imp), n_pairs = 0L),
           SUMMARY_FILE, sep = "\t")
    if (SLURM_MODE) quit(save = "no", status = 1L); next
  }

  # --- Post-process ---
  cis_dt <- as.data.table(cis_res)
  setnames(cis_dt, c("snps", "gene"), c("snp", "site"))
  # BH correction applied across all cis pairs for this panel
  cis_dt[, p_FDR := p.adjust(pvalue, method = "BH")]

  # Annotate each pair with genomic positions and physical distance
  snp_pos_map  <- setNames(snpspos$pos,   snpspos$snp)
  site_pos_map <- setNames(genepos$start, genepos$gene)
  cis_dt[, snp_pos  := snp_pos_map[snp]]
  cis_dt[, site_pos := site_pos_map[site]]
  cis_dt[, distance := abs(snp_pos - site_pos)]

  lambda <- compute_lambda(cis_dt$pvalue)
  log_msg("Lambda (genomic inflation): ", round(lambda, 4))

  sig_dt <- cis_dt[p_FDR < FDR_THRESH][order(p_FDR)]
  # Best SNP per methylation site (minimum p-value) for downstream summarisation
  top_dt <- cis_dt[, .SD[which.min(pvalue)], by = site]

  log_msg("Total cis pairs: ", nrow(cis_dt),
          " | FDR<0.05: ", nrow(sig_dt),
          " | Top per site: ", nrow(top_dt))

  # --- Save outputs ---
  saveRDS(cis_dt, file.path(panel_output, "cis_meqtl_all_results.rds"))
  fwrite(cis_dt[, .(snp, site, statistic, pvalue, FDR = p_FDR,
                    snp_pos, site_pos, distance)],
         file.path(panel_output, "cis_meqtl_all_results.tsv.gz"),
         sep = "\t", compress = "gzip")
  fwrite(sig_dt, file.path(panel_output, "cis_meqtl_significant.tsv"),  sep = "\t")
  fwrite(top_dt, file.path(panel_output, "cis_meqtl_top_per_site.tsv"), sep = "\t")

  # --- QQ plot (TIFF) ---
  pv_qq <- cis_dt$pvalue
  pv_qq <- pv_qq[is.finite(pv_qq) & pv_qq > 0]
  if (length(pv_qq) >= 2) {
    qq_file <- file.path(panel_output, paste0("qq_", panel_tag, ".tiff"))
    n_qq <- length(pv_qq)
    qq_title <- sprintf("MatrixEQTL5 QQ — %s %s\nn=%d, lambda=%.3f",
                        cohort, context, n_qq, lambda)
    draw_qq <- function() {
      # Compare observed -log10(p) distribution against uniform expectation
      plot(-log10(ppoints(n_qq)), -log10(sort(pv_qq)),
           pch = 20, cex = 0.4, col = "#2b8cbe",
           xlab = "Expected -log10(p)", ylab = "Observed -log10(p)",
           main = qq_title)
      abline(0, 1, col = "red", lty = 2)
    }
    tiff(qq_file, width = 2400, height = 2400, res = 300, compression = "lzw")
    draw_qq(); dev.off()
    pdf(sub("\\.tiff$", ".pdf", qq_file), width = 8, height = 8)
    draw_qq(); dev.off()
    cairo_ps(sub("\\.tiff$", ".eps", qq_file), width = 8, height = 8)
    draw_qq(); dev.off()
    png(sub("\\.tiff$", ".png", qq_file), width = 800, height = 800)
    draw_qq(); dev.off()
    log_msg("QQ plot: ", qq_file)
  }

  # --- Summary with lambda ---
  sum_dt <- data.table(
    cohort        = cohort,
    context       = context,
    status        = "ok",
    n_samples     = n_shared,
    n_sites       = n_sites,
    n_snps_tested = nrow(geno_imp),
    n_cis_pairs   = nrow(cis_dt),
    n_top_sites   = nrow(top_dt),
    n_sig_fdr05   = nrow(sig_dt),
    lambda        = round(lambda, 4)
  )
  fwrite(sum_dt, SUMMARY_FILE, sep = "\t")
  panel_summaries[[panel_tag]] <- sum_dt
  writeLines(capture.output(sessionInfo()), SESSION_FILE)

  rm(geno_imp, snpsSD, geneSD, cvrtSD, cis_dt, cis_res, me); gc()
  sep_line()
}

# --- Master summary ---
if (length(panel_summaries) > 0) {
  master_dt <- rbindlist(panel_summaries, fill = TRUE)
  fwrite(master_dt,
         file.path(SUMDIR, "step12ab1_all_panels_summary.tsv"), sep = "\t")
  cat("\nMatrixEQTL5 — all panels summary:\n")
  print(master_dt)
}
log_msg("Step 12ab1 finished. Outputs: ", MEQTL_ROOT)
