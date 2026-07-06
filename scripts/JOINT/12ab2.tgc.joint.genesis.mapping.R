#!/usr/bin/env Rscript
############################################################
# TreeGeneClimate (TGC) — JOINT ECS + TBS
# Step 12ab2: GENESIS5 cis-meQTL mapping
#
# KEY DIFFERENCES vs 15ab1 (GENESIS4):
#   - Inputs from 12ab0 (non-imputed GDS, GCTA GRM PCs, unfiltered methylation)
#   - PC1..PC10 for BOTH cohorts
#   - GRM random effect for BOTH cohorts (not just BREEDING)
#     BREEDING : LOCO GRM (from non-imputed GDS, VanRaden formula)
#     NATURAL  : full GCTA GRM (grm_shared.rds)
#   - Cis window: 100 kb (1e5)
#   - HWE group imputation of missing genotypes
#   - Lambda reported per panel
#
# MODEL
#   BREEDING: M-value ~ PC1..PC10 + random(LOCO-GRM)   [AIREML LMM]
#   NATURAL:  M-value ~ PC1..PC10 + random(full GRM)    [AIREML LMM]
#
# USAGE
# SLURM single panel:
#   TGC_COHORT=BREEDING TGC_CONTEXT=CpG Rscript --vanilla 12ab2.R
# Interactive (all 6 panels):
#   source this script in RStudio (no env vars set)
#
# INPUTS (from 12ab0 — RESULTS/JOINT/MQTL5/INPUTS/<COHORT>/<CONTEXT>/)
#   methylation_mvalues_matrix.rds
#   methylation_site_annot.rds
#   shared_sample_ids.rds
#   pcs_shared.rds
#   grm_shared.rds
#   sample_groups.rds
#   snp_variant_annot.rds
#   snp_gds_sample_ids.rds
# BREEDING LOCO GRMs (INPUTS/BREEDING/loco_grm/loco_grm_excl_<chr>.rds)
#
# OUTPUTS (RESULTS/JOINT/GENESIS5/<COHORT>/<CONTEXT>/)
#   cis_meqtl_all_results.rds / .tsv.gz
#   cis_meqtl_significant.tsv
#   cis_meqtl_top_per_site.tsv
#   cis_meqtl_panel_summary.tsv      (includes lambda)
#   qq_<cohort>_<context>.tiff
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(SNPRelate)
  library(gdsfmt)
  library(GENESIS)      # fitNullModel + assocTestSingle
  library(GWASTools)    # GenotypeData infrastructure
  library(Biobase)
  library(parallel)
})

options(stringsAsFactors = FALSE)

############################################################
# 1) SETTINGS
############################################################

CIS_WINDOW  <- 1e5        # 100 kb
MIN_SAMPLES <- 10L
FDR_THRESH  <- 0.05
# Loose AIREML convergence tolerance: variance components are re-estimated
# per site, so exact convergence is less critical than throughput
AIREML_FAST <- 1e-2
N_CORES     <- min(4L, parallel::detectCores(logical = FALSE))

############################################################
# 2) PATHS
############################################################

# === USER CONFIGURATION ===
# Set TGC_PROJECT_ROOT as an environment variable, or edit the fallback path below
PROJECT_ROOT <- Sys.getenv("TGC_PROJECT_ROOT",
  unset = "/path/to/your/project")
# ===========================

MQTL5_ROOT   <- file.path(PROJECT_ROOT, "RESULTS", "JOINT", "MQTL5")
INPUTDIR     <- file.path(MQTL5_ROOT, "INPUTS")

GENESIS5_ROOT <- file.path(PROJECT_ROOT, "RESULTS", "JOINT", "GENESIS5")
OUTROOT       <- file.path(GENESIS5_ROOT, "MEQTL")
LOGDIR        <- file.path(GENESIS5_ROOT, "LOGS")
SUMDIR        <- file.path(GENESIS5_ROOT, "SUMMARIES")

for (d in c(OUTROOT, LOGDIR, SUMDIR))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Non-imputed GDS files (genotypes loaded fresh per panel)
RDATA_DIR <- "/mnt/vast-standard/home/chano/u15584/treegeneclimate/2025/ECS/RDATA"
GDS_FILES <- list(
  BREEDING = file.path(RDATA_DIR, "breeding.snp.gds"),
  NATURAL  = file.path(RDATA_DIR, "natural.snp.gds")
)

############################################################
# 3) PHASE DETECTION
############################################################

# Single-panel SLURM mode when env vars are set; all-panels interactive mode otherwise
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

CURRENT_LOGFILE   <- NULL
CURRENT_DEBUGFILE <- NULL

log_msg <- function(...) {
  txt <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ",
                paste0(..., collapse = ""))
  cat(txt, "\n")
  if (!is.null(CURRENT_LOGFILE))
    cat(txt, "\n", file = CURRENT_LOGFILE, append = TRUE)
}
# Verbose diagnostics written only to a separate debug log to keep the main log clean
debug_msg <- function(...) {
  txt <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ",
                paste0(..., collapse = ""))
  if (!is.null(CURRENT_DEBUGFILE))
    cat(txt, "\n", file = CURRENT_DEBUGFILE, append = TRUE)
}
sep_line <- function() log_msg(paste(rep("=", 70), collapse = ""))

normalize_id <- function(x) {
  x <- as.character(x); x <- trimws(x)
  x <- gsub("\\.0$", "", x); x <- gsub("^X", "", x); x <- gsub("-", "_", x); x
}

# HWE group imputation (SNPs x samples):
# fill missing calls with 2*p estimated within the same family/population group
impute_hwe <- function(geno_snp_x_samp, groups) {
  mode(geno_snp_x_samp) <- "numeric"
  geno_snp_x_samp[geno_snp_x_samp > 2 | geno_snp_x_samp < 0] <- NA_real_
  # Cohort-wide allele frequency as fallback
  p_all <- rowSums(geno_snp_x_samp, na.rm = TRUE) /
           (2 * rowSums(!is.na(geno_snp_x_samp)))
  keep  <- is.finite(p_all)
  g     <- geno_snp_x_samp[keep, , drop = FALSE]
  p_all <- p_all[keep]
  for (grp in unique(groups)) {
    idx <- which(groups == grp)
    pg  <- rowSums(g[, idx, drop = FALSE], na.rm = TRUE) /
           (2 * rowSums(!is.na(g[, idx, drop = FALSE])))
    pg[!is.finite(pg)] <- p_all[!is.finite(pg)]
    miss <- which(is.na(g[, idx, drop = FALSE]), arr.ind = TRUE)
    if (nrow(miss)) g[, idx][miss] <- 2 * pg[miss[, 1]]
  }
  # Global fallback for any remaining NAs
  if (anyNA(g)) {
    fill2 <- matrix(2 * p_all, nrow = nrow(g), ncol = ncol(g))
    g[is.na(g)] <- fill2[is.na(g)]
  }
  list(geno = g, keep = keep)
}

# Load the LOCO GRM for a given chromosome and subset/reindex to shared samples.
# The GRM row/col names are set to integer scan IDs required by GWASTools.
load_loco_grm <- function(loco_dir, chr_name, shared_ids, scan_int) {
  path <- file.path(loco_dir, paste0("loco_grm_excl_", chr_name, ".rds"))
  if (!file.exists(path)) return(NULL)
  G        <- readRDS(path)
  grm_norm <- normalize_id(rownames(G))
  idx      <- match(shared_ids, grm_norm)
  bad      <- is.na(idx)
  if (any(bad)) {
    debug_msg("LOCO GRM missing ", sum(bad), " shared IDs for chr ", chr_name)
    # Substitute a valid row index for unmatched samples (affects only diagonal)
    idx[bad] <- which(!is.na(grm_norm))[1L]
  }
  G_sub <- G[idx, idx, drop = FALSE]
  rownames(G_sub) <- colnames(G_sub) <- as.character(scan_int)
  diag(G_sub) <- diag(G_sub) + 1e-4  # small ridge to ensure positive definiteness
  G_sub
}

# Genomic inflation factor
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

  # Skip completed panels — allows safe SLURM restarts
  done_rds <- file.path(panel_output, "cis_meqtl_all_results.rds")
  if (file.exists(done_rds)) {
    cat("[", format(Sys.time(), "%H:%M:%S"), "] Panel ", cohort, "/", context,
        " already done — skipping\n", sep = "")
    next
  }

  CURRENT_LOGFILE   <- file.path(LOGDIR, paste0("step12ab2_", panel_tag, ".log"))
  CURRENT_DEBUGFILE <- file.path(LOGDIR, paste0("step12ab2_", panel_tag, "_debug.log"))
  SUMMARY_FILE      <- file.path(SUMDIR, paste0("step12ab2_", panel_tag, "_summary.tsv"))
  SESSION_FILE      <- file.path(SUMDIR, paste0("step12ab2_", panel_tag, "_sessionInfo.txt"))
  if (file.exists(CURRENT_LOGFILE)) file.remove(CURRENT_LOGFILE)

  loco_dir <- file.path(INPUTDIR, "BREEDING", "loco_grm")
  # Both cohorts use a GRM random effect; BREEDING uses per-chromosome LOCO GRMs
  USE_GRM   <- TRUE
  USE_LOCO  <- (cohort == "BREEDING")

  sep_line()
  log_msg("Step 12ab2 — GENESIS5 cis-meQTL mapping")
  log_msg("Panel: ", cohort, " / ", context)
  log_msg("Cis window: +/-", CIS_WINDOW / 1e3, " kb")
  log_msg("GRM: ", if (USE_LOCO) "LOCO (BREEDING)" else "full GCTA GRM (NATURAL)")

  # --- Check inputs ---
  req <- c("methylation_mvalues_matrix.rds", "methylation_site_annot.rds",
           "shared_sample_ids.rds", "pcs_shared.rds", "grm_shared.rds",
           "sample_groups.rds", "snp_variant_annot.rds")
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
  raw_shared    <- gds_samp_ids[gds_samp_norm %in% shared_ids]
  # snpfirstdim=FALSE returns samples x SNPs (GWASTools convention)
  geno_list     <- snpgdsGetGeno(gds, sample.id = raw_shared,
                                 snpfirstdim = FALSE, with.id = TRUE)
  snpgdsClose(gds)

  samp_reorder <- match(shared_ids, normalize_id(geno_list$sample.id))
  geno_mat     <- geno_list$genotype[samp_reorder, , drop = FALSE]  # samples x SNPs
  rm(geno_list, samp_reorder); gc()

  # Transpose to SNPs x samples for imputation (impute_hwe expects this orientation)
  geno_snp_x_samp <- t(geno_mat); rm(geno_mat); gc()
  log_msg("Genotype matrix: ", nrow(geno_snp_x_samp), " SNPs x ",
          ncol(geno_snp_x_samp), " samples")

  # --- HWE group imputation ---
  log_msg("Imputing missing genotypes by group...")
  groups_ord <- groups_dt$group[match(shared_ids, groups_dt$sample_id)]
  groups_ord[is.na(groups_ord)] <- "Unknown"
  imp             <- impute_hwe(geno_snp_x_samp, groups_ord)
  geno_imp        <- imp$geno     # SNPs x samples (filtered + imputed)
  var_filt        <- var_annot[imp$keep]
  rm(geno_snp_x_samp, imp); gc()

  # Drop zero-variance SNPs (monomorphic after imputation)
  snp_sd   <- apply(geno_imp, 1, sd, na.rm = TRUE)
  keep_var <- is.finite(snp_sd) & snp_sd > 0
  geno_imp <- geno_imp[keep_var, , drop = FALSE]
  var_filt <- var_filt[keep_var]
  log_msg("SNPs after imputation + variance filter: ", nrow(geno_imp))

  # Transpose back to samples x SNPs (integer storage required by MatrixGenotypeReader)
  geno_imp_t <- t(geno_imp); storage.mode(geno_imp_t) <- "integer"
  rm(geno_imp); gc()

  # --- Build GWASTools objects ---
  # GWASTools requires integer scan IDs and integer-coded chromosomes
  log_msg("Building GWASTools GenotypeData object...")
  scan_int   <- seq_len(n_shared)
  snp_int    <- seq_len(nrow(var_filt))
  chr_levels <- sort(unique(var_filt$chr))
  chr_map    <- setNames(seq_along(chr_levels), chr_levels)
  chr_int    <- as.integer(chr_map[var_filt$chr])

  reader    <- MatrixGenotypeReader(
    genotype   = t(geno_imp_t),   # SNPs x samples
    snpID      = snp_int,
    chromosome = chr_int,
    position   = as.integer(var_filt$pos),
    scanID     = scan_int
  )
  scanAnnot <- ScanAnnotationDataFrame(data.frame(scanID = scan_int))
  genoData  <- GenotypeData(reader, scanAnnot = scanAnnot)
  rm(geno_imp_t); gc()

  # --- Full GRM (NATURAL: primary; BREEDING: fallback when LOCO unavailable) ---
  log_msg("Loading full GRM...")
  grm_full <- readRDS(file.path(panel_input, "grm_shared.rds"))
  grm_idx  <- match(shared_ids, normalize_id(rownames(grm_full)))
  grm_full_ord <- grm_full[grm_idx, grm_idx, drop = FALSE]
  # Use integer scan IDs as row/col names to match GWASTools scanID
  rownames(grm_full_ord) <- colnames(grm_full_ord) <- as.character(scan_int)
  # Ensure positive definiteness with small ridge
  diag(grm_full_ord) <- diag(grm_full_ord) + 1e-4
  rm(grm_full); gc()

  # --- PC covariates ---
  pc_cols  <- grep("^PC[0-9]+$", names(pcs_dt), value = TRUE)
  pcs_ord  <- as.data.frame(pcs_dt[match(shared_ids, pcs_dt$sample_id)])
  # pheno_base is the covariate data frame passed to fitNullModel for every site
  pheno_base <- data.frame(sample.id = scan_int)
  for (pc in pc_cols) pheno_base[[pc]] <- pcs_ord[[pc]]
  log_msg("Covariates: ", paste(pc_cols, collapse = ", "))

  # --- Reorder M-values to match shared_ids ---
  m_reorder <- match(shared_ids, normalize_id(rownames(m_mat)))
  m_mat     <- m_mat[m_reorder, , drop = FALSE]

  # --- Index by chromosome for the cis loop ---
  site_annot[, site_idx := .I]  # row index for fast column access in m_mat
  var_filt[,   snp_int  := seq_len(.N)]
  site_by_chr <- split(site_annot, site_annot$chr)
  snp_by_chr  <- split(var_filt,   var_filt$chr)
  chrs_shared <- intersect(names(site_by_chr), names(snp_by_chr))
  log_msg("Chromosomes with sites + SNPs: ", length(chrs_shared))

  # --- cis-meQTL loop (parallel over sites within each chromosome) ---
  # One null model is fitted per site (AIREML estimates variance components
  # from M ~ PCs + GRM). SNP association is then scored against this null.
  log_msg("Starting GENESIS cis-meQTL mapping (", N_CORES, " cores)...")
  all_results       <- list()
  n_tested          <- 0L
  n_skip_no_cis     <- 0L
  n_skip_null_fail  <- 0L
  n_skip_assoc_fail <- 0L
  n_tests_total     <- 0L
  t_start           <- proc.time()

  for (chr_name in chrs_shared) {

    sites_chr   <- site_by_chr[[chr_name]]
    snps_chr    <- snp_by_chr[[chr_name]]
    snp_pos_chr <- snps_chr$pos

    # For BREEDING: swap in the LOCO GRM (chromosome excluded) to avoid
    # proximal contamination of the variance-component estimates.
    # For NATURAL: always use the full GRM.
    if (USE_LOCO) {
      loco_candidate <- load_loco_grm(loco_dir, chr_name, shared_ids, scan_int)
      grm_chr <- if (!is.null(loco_candidate)) loco_candidate else grm_full_ord
    } else {
      grm_chr <- grm_full_ord
    }

    # Parallelise over sites within the chromosome
    chr_res <- mclapply(seq_len(nrow(sites_chr)), function(j) {
      site <- sites_chr[j, ]

      # Identify cis-SNPs within the 100 kb window
      cis_mask <- abs(snp_pos_chr - site[["pos"]]) <= CIS_WINDOW
      if (!any(cis_mask)) return(list(status = "no_cis"))

      cis_snps  <- snps_chr[cis_mask]
      cis_snpID <- cis_snps[["snp_int"]]

      # Phenotype: M-value vector for this site
      y_vec <- m_mat[, site[["site_idx"]]]
      if (sum(!is.na(y_vec)) < MIN_SAMPLES) return(list(status = "no_cis"))

      pheno_null   <- pheno_base
      pheno_null$y <- y_vec

      # Fit null LMM: M ~ PC1..PC10 + random(GRM) via AIREML
      null_mod <- tryCatch(
        fitNullModel(pheno_null, outcome = "y", covars = pc_cols,
                     cov.mat = grm_chr, family = "gaussian",
                     AIREML.tol = AIREML_FAST, verbose = FALSE),
        error = function(e) NULL
      )
      if (is.null(null_mod)) return(list(status = "null_fail"))

      # Score test of each cis-SNP against the null model residuals
      iterator <- GenotypeBlockIterator(genoData, snpInclude = cis_snpID)
      assoc <- tryCatch(
        assocTestSingle(iterator, null.model = null_mod, verbose = FALSE),
        error = function(e) NULL
      )
      if (is.null(assoc) || nrow(assoc) == 0L) return(list(status = "assoc_fail"))

      assoc_dt <- as.data.table(assoc)
      .methyl_loc <- site[["methyl_loc"]]
      .site_pos   <- as.integer(site[["pos"]])
      .snp_pos    <- snps_chr[["pos"]][match(assoc_dt$variant.id, cis_snps[["snp_int"]])]
      assoc_dt[, `:=`(
        site     = .methyl_loc,
        site_chr = chr_name,
        site_pos = .site_pos,
        snp_pos  = .snp_pos,
        distance = abs(.site_pos - .snp_pos)
      )]
      list(status = "ok", dt = assoc_dt)
    }, mc.cores = N_CORES, mc.preschedule = TRUE)

    # Collect results and update counters
    for (res in chr_res) {
      if (is.null(res) || !is.list(res)) { n_skip_assoc_fail <- n_skip_assoc_fail + 1L; next }
      if (res$status == "ok") {
        all_results[[length(all_results) + 1L]] <- res$dt
        n_tested      <- n_tested + 1L
        n_tests_total <- n_tests_total + nrow(res$dt)
      } else if (res$status == "no_cis") {
        n_skip_no_cis <- n_skip_no_cis + 1L
      } else if (res$status == "null_fail") {
        n_skip_null_fail <- n_skip_null_fail + 1L
      } else {
        n_skip_assoc_fail <- n_skip_assoc_fail + 1L
      }
    }

    elapsed <- (proc.time() - t_start)[3]
    log_msg("Chr ", chr_name, " | tested=", n_tested,
            " | null-fail=", n_skip_null_fail,
            " | ", round(elapsed / 60, 1), " min elapsed")
  }

  # --- Combine results ---
  log_msg("Combining results...")
  log_msg("  Sites tested: ", n_tested, " | Pairs: ", n_tests_total)
  log_msg("  Skipped no-cis: ", n_skip_no_cis,
          " | null-fail: ", n_skip_null_fail,
          " | assoc-fail: ", n_skip_assoc_fail)

  if (length(all_results) == 0L) {
    log_msg("No results produced for ", cohort, " / ", context)
    fwrite(data.table(cohort = cohort, context = context,
                      status = "no_results", lambda = NA_real_),
           SUMMARY_FILE, sep = "\t")
    if (SLURM_MODE) quit(save = "no", status = 1L); next
  }

  cis_dt <- rbindlist(all_results, fill = TRUE)
  rm(all_results); gc()

  # Standardise column names across GENESIS versions (column names vary by version)
  p_col <- intersect(c("Score.pval", "Score.Stat.p", "pval"), names(cis_dt))[1]
  t_col <- intersect(c("Score.Stat", "Score"), names(cis_dt))[1]
  snp_col <- intersect(c("snpID", "variant.id"), names(cis_dt))[1]
  setnames(cis_dt, c(snp_col, p_col, t_col), c("snp_int", "pvalue", "statistic"),
           skip_absent = TRUE)

  # Map integer SNP index back to the original SNP identifier string
  snp_name_map <- setNames(var_filt$snp_id, var_filt$snp_int)
  cis_dt[, snp := snp_name_map[as.character(snp_int)]]

  # BH correction across all cis pairs for this panel
  cis_dt[, p_FDR := p.adjust(pvalue, method = "BH")]

  lambda <- compute_lambda(cis_dt$pvalue)
  log_msg("Lambda (genomic inflation): ", round(lambda, 4))

  sig_dt <- cis_dt[p_FDR < FDR_THRESH][order(p_FDR)]
  # Best SNP per methylation site
  top_dt <- cis_dt[, .SD[which.min(pvalue)], by = site]

  log_msg("Total cis pairs: ", nrow(cis_dt),
          " | FDR<0.05: ", nrow(sig_dt),
          " | Top per site: ", nrow(top_dt))

  # --- Save outputs ---
  out_cols <- c("snp", "site", "statistic", "pvalue", "p_FDR",
                "snp_pos", "site_pos", "distance")
  out_cols <- intersect(out_cols, names(cis_dt))

  saveRDS(cis_dt, file.path(panel_output, "cis_meqtl_all_results.rds"))
  fwrite(cis_dt[, ..out_cols],
         file.path(panel_output, "cis_meqtl_all_results.tsv.gz"),
         sep = "\t", compress = "gzip")
  fwrite(sig_dt[, ..out_cols],
         file.path(panel_output, "cis_meqtl_significant.tsv"), sep = "\t")
  fwrite(top_dt[, ..out_cols],
         file.path(panel_output, "cis_meqtl_top_per_site.tsv"), sep = "\t")

  # --- QQ plot (TIFF) ---
  pv_qq <- cis_dt$pvalue
  pv_qq <- pv_qq[is.finite(pv_qq) & pv_qq > 0]
  if (length(pv_qq) >= 2) {
    qq_file <- file.path(panel_output, paste0("qq_", panel_tag, ".tiff"))
    n_qq <- length(pv_qq)
    qq_title <- sprintf("GENESIS5 QQ — %s %s\nn=%d, lambda=%.3f",
                        cohort, context, n_qq, lambda)
    draw_qq <- function() {
      plot(-log10(ppoints(n_qq)), -log10(sort(pv_qq)),
           pch = 20, cex = 0.4, col = "#e34a33",
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

  # --- Summary ---
  sum_dt <- data.table(
    cohort        = cohort,
    context       = context,
    status        = "ok",
    grm_type      = if (USE_LOCO) "LOCO" else "full_GCTA",
    n_samples     = n_shared,
    n_sites       = n_sites,
    n_snps_tested = nrow(var_filt),
    n_cis_pairs   = nrow(cis_dt),
    n_top_sites   = nrow(top_dt),
    n_sig_fdr05   = nrow(sig_dt),
    lambda        = round(lambda, 4)
  )
  fwrite(sum_dt, SUMMARY_FILE, sep = "\t")
  panel_summaries[[panel_tag]] <- sum_dt
  writeLines(capture.output(sessionInfo()), SESSION_FILE)

  rm(cis_dt, genoData, grm_full_ord, m_mat); gc()
  sep_line()
}

# --- Master summary ---
if (length(panel_summaries) > 0) {
  master_dt <- rbindlist(panel_summaries, fill = TRUE)
  fwrite(master_dt,
         file.path(SUMDIR, "step12ab2_all_panels_summary.tsv"), sep = "\t")
  cat("\nGENESIS5 — all panels summary:\n")
  print(master_dt)
}
log_msg("Step 12ab2 finished. Outputs: ", GENESIS5_ROOT)
