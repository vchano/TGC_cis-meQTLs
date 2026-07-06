#!/usr/bin/env Rscript
############################################################
# TreeGeneClimate (TGC) — JOINT ECS + TBS
# Step 12ab0: Prepare meQTL5 inputs (OLD-style inputs)
#
# KEY DIFFERENCES vs 15ab0 (meQTL4):
#   - GDS     : non-imputed  breeding.snp.gds / natural.snp.gds
#   - GRM     : GCTA GRM     breeding_grm_gcta.rds / natural_grm_gcta.rds
#   - Methylation: unfiltered {cohort}_{ctx}_methylkit.rds (no MEF filter)
#   - PCs     : 10 for BOTH cohorts, derived from GCTA GRM eigendecomposition
#   - NO quantile normalisation of M-values
#   - Group assignment (Family / Population) saved per panel for HWE imputation
#   - LOCO GRMs for BREEDING computed from non-imputed GDS (VanRaden formula)
#   - Full GCTA GRM saved for NATURAL (used by GENESIS5)
#
# OUTPUTS (under RESULTS/JOINT/MQTL5/INPUTS/<COHORT>/<CONTEXT>/)
#   methylation_mvalues_matrix.rds / .tsv
#   methylation_site_annot.rds / .tsv
#   shared_sample_ids.rds / .tsv
#   pcs_shared.rds / .tsv            (PC1..PC10, from GRM eigen)
#   grm_shared.rds                   (GCTA GRM subset to shared samples)
#   sample_groups.rds / .tsv         (group assignments for HWE imputation)
#   snp_gds_sample_ids.rds / .tsv
#   snp_variant_annot.rds / .tsv
# Per-cohort LOCO GRMs (BREEDING only):
#   INPUTS/BREEDING/loco_grm/loco_grm_excl_<chr>.rds
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(SNPRelate)
  library(gdsfmt)
  library(methylKit)
})

options(stringsAsFactors = FALSE)
options(datatable.fread.datatable = TRUE)

############################################################
# 1) SETTINGS
############################################################

# 10 PCs used as fixed-effect covariates in both MatrixEQTL and GENESIS models
N_PCS    <- 10L

# Pseudo-count for logit transformation: M = log((pct + eps) / (100 - pct + eps))
# eps = 0.5 avoids log(0) at boundary values without over-shrinking
EPSILON_M <- 0.5

COHORTS  <- c("BREEDING", "NATURAL")
CONTEXTS <- c("CpG", "CHG", "CHH")

############################################################
# 2) PATHS
############################################################

# === USER CONFIGURATION ===
PROJECT_ROOT <- "/path/to/your/project"  # <-- set this
# ===========================
OUTROOT      <- file.path(PROJECT_ROOT, "RESULTS", "JOINT", "MQTL5")

LOGDIR   <- file.path(OUTROOT, "LOGS")
SUMDIR   <- file.path(OUTROOT, "SUMMARIES")
DBGDIR   <- file.path(OUTROOT, "DEBUG")
INPUTDIR <- file.path(OUTROOT, "INPUTS")

for (d in c(OUTROOT, LOGDIR, SUMDIR, DBGDIR, INPUTDIR))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

LOGFILE      <- file.path(LOGDIR, "step12ab0.log")
SUMMARY_FILE <- file.path(SUMDIR, "step12ab0_summary.tsv")
SESSION_FILE <- file.path(SUMDIR, "step12ab0_sessionInfo.txt")
OVERLAP_FILE <- file.path(DBGDIR, "step12ab0_overlap_summary.tsv")

# Clear any logs from a previous run before starting
for (f in c(LOGFILE, SUMMARY_FILE, OVERLAP_FILE))
  if (file.exists(f)) file.remove(f)

# --- OLD-STYLE input paths ---
RDATA_DIR   <- "/mnt/vast-standard/home/chano/u15584/treegeneclimate/2025/ECS/RDATA"
TBS_RDS_DIR <- "/scratch-scc/users/u15584/TGC/TBS/2025_NEW.ANALYSIS/METHYLKIT.FILES.TBS"

# Non-imputed GDS files — imputed GDS is used only in step 11ab
GDS_FILES <- list(
  BREEDING = file.path(RDATA_DIR, "breeding.snp.gds"),
  NATURAL  = file.path(RDATA_DIR, "natural.snp.gds")
)
# GCTA GRMs pre-computed from the full SNP set
GRM_FILES <- list(
  BREEDING = file.path(RDATA_DIR, "breeding_grm_gcta.rds"),
  NATURAL  = file.path(RDATA_DIR, "natural_grm_gcta.rds")
)
# Sample metadata (IID, Family/Population) for group-aware HWE imputation
ANNOT_FILES <- list(
  BREEDING = file.path(RDATA_DIR, "breeding_sample_annotation.rds"),
  NATURAL  = file.path(RDATA_DIR, "natural_sample_annotation.rds")
)

############################################################
# 3) HELPERS
############################################################

log_msg <- function(...) {
  txt <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ",
                paste0(..., collapse = ""))
  cat(txt, "\n")
  cat(txt, "\n", file = LOGFILE, append = TRUE)
}
sep_line <- function() log_msg(paste(rep("=", 70), collapse = ""))

# Canonical sample ID normalisation: remove trailing ".0" (Excel artefact),
# leading "X" (R's make.names prefix), replace "-" with "_"
normalize_id <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("\\.0$", "", x)
  x <- gsub("^X",   "", x)
  x <- gsub("-",    "_", x)
  x
}

# Convert a methylKit methylBase object to a logit M-value matrix (samples x sites).
# Uses non-quantile-normalised values; sites with zero variance are dropped.
methylkit_to_mvalues <- function(meth_path, eps = 0.5) {
  if (!file.exists(meth_path)) stop("Missing methylKit RDS: ", meth_path)
  mb   <- readRDS(meth_path)
  d    <- methylKit::getData(mb)
  cs_cols <- grep("^numCs[0-9]+$", colnames(d), value = TRUE)
  ts_cols <- grep("^numTs[0-9]+$", colnames(d), value = TRUE)
  stopifnot(length(cs_cols) > 0, length(cs_cols) == length(ts_cols))
  sample_ids <- mb@sample.ids
  stopifnot(length(sample_ids) == length(cs_cols))
  site_annot <- data.table(
    chr   = as.character(d$chr),
    start = as.integer(d$start),
    end   = as.integer(d$end)
  )
  site_annot[, pos       := start]
  site_annot[, methyl_loc := paste0(chr, ":", start, "-", end)]
  # Compute percentage methylation: 100 * C / (C + T)
  perc_mat <- vapply(seq_along(cs_cols), function(i) {
    cs  <- d[[cs_cols[i]]]; ts <- d[[ts_cols[i]]]
    pct <- 100 * cs / (cs + ts)
    pct[is.nan(pct)] <- NA_real_
    pct
  }, numeric(nrow(d)))
  colnames(perc_mat) <- sample_ids
  # Remove sites with zero variance or all-NA
  keep_var <- apply(perc_mat, 1, function(x) sum(!is.na(x)) > 1 && sd(x, na.rm = TRUE) > 0)
  perc_mat  <- perc_mat[keep_var, , drop = FALSE]
  site_annot <- site_annot[keep_var]
  # Logit transform: M-values are approximately normally distributed and
  # better suited for linear models than beta/percentage values
  m_mat <- log((perc_mat + eps) / (100 - perc_mat + eps))
  rownames(m_mat) <- site_annot$methyl_loc
  m_mat <- t(m_mat)   # samples x sites
  rownames(m_mat) <- sample_ids
  list(m_matrix = m_mat, site_annot = site_annot)
}

# Derive population structure PCs by eigendecomposition of the GRM.
# This avoids re-running GCTA --pca and keeps PCs on the same scale as the GRM.
grm_pcs <- function(G, sample_norm, n_pcs) {
  # G: matrix with rownames (raw or normalised IDs)
  # Returns data.table: sample_id_raw, sample_id, PC1..PCn_pcs
  ev  <- eigen(G, symmetric = TRUE)
  pcs <- ev$vectors[, seq_len(min(n_pcs, ncol(ev$vectors))), drop = FALSE]
  colnames(pcs) <- paste0("PC", seq_len(ncol(pcs)))
  rownames(pcs) <- sample_norm
  dt  <- as.data.table(pcs, keep.rownames = "sample_id")
  dt
}

############################################################
# 4) LOCO GRM — BREEDING (VanRaden formula, non-imputed GDS)
############################################################

# LOCO (leave-one-chromosome-out) GRMs are required by the GENESIS LMM to avoid
# proximal contamination: the chromosome being tested is excluded from the kinship
# matrix, preventing the GRM from inadvertently capturing the tested signal.
# Computed only for BREEDING because relatedness structure is stronger in family material.
compute_loco_grms <- function(gds_path, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  log_msg("Computing LOCO GRMs from: ", gds_path)
  gds <- snpgdsOpen(gds_path, readonly = TRUE)
  on.exit(snpgdsClose(gds), add = TRUE)

  sample_ids <- read.gdsn(index.gdsn(gds, "sample.id"))
  all_chr    <- as.character(read.gdsn(index.gdsn(gds, "snp.chromosome")))
  chrs       <- sort(unique(all_chr))
  log_msg("  Samples: ", length(sample_ids), " | Chromosomes: ", length(chrs))

  log_msg("  Reading all genotypes...")
  gl   <- snpgdsGetGeno(gds, snpfirstdim = TRUE, with.id = TRUE)
  geno <- gl$genotype   # SNPs x samples
  sample_norm <- normalize_id(gl$sample.id)
  log_msg("  Genotype matrix: ", nrow(geno), " SNPs x ", ncol(geno), " samples")

  # VanRaden (2008) GRM: G = Z Z' / (2 sum p_i (1-p_i))
  # where Z = (g - 2p) is the mean-centered genotype matrix
  vanraden_grm <- function(geno_sub) {
    p    <- rowSums(geno_sub, na.rm = TRUE) / (2 * rowSums(!is.na(geno_sub)))
    poly <- p > 0 & p < 1 & !is.na(p)  # exclude fixed/monomorphic SNPs
    if (sum(poly) < 2L) return(NULL)
    g2 <- geno_sub[poly, , drop = FALSE]
    p2 <- p[poly]
    Xc <- g2 - 2 * p2   # center by 2*allele frequency
    Xc[is.na(Xc)] <- 0  # set missing genotypes to zero after centering
    denom <- 2 * sum(p2 * (1 - p2))
    if (!is.finite(denom) || denom <= 0) return(NULL)
    G <- (t(Xc) %*% Xc) / denom
    rownames(G) <- colnames(G) <- sample_norm
    G
  }

  n_done <- 0L
  for (excl_chr in chrs) {
    out_path <- file.path(out_dir, paste0("loco_grm_excl_", excl_chr, ".rds"))
    if (file.exists(out_path)) {
      log_msg("  LOCO chr ", excl_chr, ": already exists — skipping")
      n_done <- n_done + 1L; next
    }
    # Build GRM from all SNPs except those on the chromosome being tested
    mask        <- all_chr != excl_chr
    n_snps_used <- sum(mask)
    if (n_snps_used < 10L) {
      log_msg("  LOCO chr ", excl_chr, ": only ", n_snps_used, " SNPs — skipping"); next
    }
    G <- vanraden_grm(geno[mask, , drop = FALSE])
    if (is.null(G)) {
      log_msg("  LOCO chr ", excl_chr, ": GRM failed — skipping"); next
    }
    saveRDS(G, out_path)
    log_msg("  LOCO chr ", excl_chr, ": ", nrow(G), "x", ncol(G),
            " (", n_snps_used, " SNPs) — saved")
    n_done <- n_done + 1L
  }
  log_msg("LOCO GRMs: ", n_done, " / ", length(chrs), " chromosomes saved")
  invisible(out_dir)
}

############################################################
# 5) LOAD COHORT-LEVEL OBJECTS
############################################################

sep_line()
log_msg("Step 12ab0 — meQTL5 input preparation (OLD-style inputs)")
log_msg("Output root: ", OUTROOT)

# Pre-load per-cohort objects once (GDS metadata, GRM, PCs, annotations)
# to avoid re-reading for each of the three contexts
gds_meta <- list()
grm_pcs_meta <- list()
grm_meta <- list()
annot_meta <- list()

for (cohort in COHORTS) {
  sep_line()
  log_msg("Loading cohort-level objects: ", cohort)

  # Read SNP and sample metadata from GDS (genotypes are not loaded here)
  gds_path <- GDS_FILES[[cohort]]
  stopifnot(file.exists(gds_path))
  gds      <- snpgdsOpen(gds_path, readonly = TRUE)
  samp_ids <- read.gdsn(index.gdsn(gds, "sample.id"))
  snp_ids  <- read.gdsn(index.gdsn(gds, "snp.id"))
  chr_ids  <- as.character(read.gdsn(index.gdsn(gds, "snp.chromosome")))
  pos_ids  <- as.integer(read.gdsn(index.gdsn(gds, "snp.position")))
  snpgdsClose(gds)

  gds_meta[[cohort]] <- list(
    sample_ids    = data.table(sample_id_raw = as.character(samp_ids),
                               sample_id     = normalize_id(samp_ids)),
    variant_annot = data.table(snp_id = as.character(snp_ids),
                               chr    = chr_ids,
                               pos    = pos_ids)
  )
  log_msg("  GDS: ", nrow(gds_meta[[cohort]]$sample_ids), " samples, ",
          nrow(gds_meta[[cohort]]$variant_annot), " SNPs")

  # GCTA GRM → derive 10 PCs via eigendecomposition
  grm_path <- GRM_FILES[[cohort]]
  stopifnot(file.exists(grm_path))
  G <- readRDS(grm_path)
  stopifnot(is.matrix(G) || inherits(G, "Matrix"))
  stopifnot(!is.null(rownames(G)))
  grm_meta[[cohort]] <- G

  grm_ids_norm <- normalize_id(rownames(G))
  pcs_dt <- grm_pcs(G, grm_ids_norm, N_PCS)
  grm_pcs_meta[[cohort]] <- pcs_dt
  log_msg("  GCTA GRM: ", nrow(G), "x", ncol(G),
          " | PCs derived: PC1..PC", N_PCS)

  # Group column differs by cohort: Family (BREEDING) or Population (NATURAL)
  # used later for within-group HWE imputation of missing genotypes
  annot_path <- ANNOT_FILES[[cohort]]
  stopifnot(file.exists(annot_path))
  annot <- readRDS(annot_path)
  group_col <- if (cohort == "BREEDING") "Family" else "Population"
  stopifnot("IID" %in% names(annot), group_col %in% names(annot))
  annot_dt <- data.table(
    sample_id = normalize_id(as.character(annot$IID)),
    group     = as.character(annot[[group_col]])
  )
  annot_meta[[cohort]] <- annot_dt
  log_msg("  Annotation: ", nrow(annot_dt), " samples, group col = ", group_col)
}

############################################################
# 6) LOCO GRM — BREEDING ONLY
############################################################

sep_line()
log_msg("Computing LOCO GRMs for BREEDING cohort (non-imputed GDS)...")
loco_dir <- file.path(INPUTDIR, "BREEDING", "loco_grm")
compute_loco_grms(GDS_FILES[["BREEDING"]], loco_dir)

############################################################
# 7) PANEL-SPECIFIC INPUTS
############################################################

summary_rows <- list()
overlap_rows <- list()

for (cohort in COHORTS) {
  for (context in CONTEXTS) {
    sep_line()
    log_msg("Panel: ", cohort, " / ", context)
    panel_dir <- file.path(INPUTDIR, cohort, context)
    dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)

    # --- Load unfiltered methylKit RDS → M-values (NO QN) ---
    ctx_lower <- tolower(context)
    coh_lower <- tolower(cohort)
    meth_path <- file.path(TBS_RDS_DIR,
                           paste0(coh_lower, "_", context, "_methylkit.rds"))
    if (!file.exists(meth_path)) {
      log_msg("  Missing methylKit RDS — skipping: ", meth_path)
      summary_rows[[paste(cohort, context)]] <- data.table(
        cohort = cohort, context = context, status = "missing_methylbase")
      next
    }
    log_msg("  methylKit file: ", meth_path)
    conv       <- methylkit_to_mvalues(meth_path, eps = EPSILON_M)
    m_mat      <- conv$m_matrix
    site_annot <- conv$site_annot
    log_msg("  M-value matrix: ", nrow(m_mat), " samples x ", ncol(m_mat), " sites")

    # --- Sample harmonisation ---
    # Find the intersection of samples present in all four data sources
    meth_ids_norm <- normalize_id(rownames(m_mat))
    gds_ids  <- gds_meta[[cohort]]$sample_ids$sample_id
    pcs_ids  <- grm_pcs_meta[[cohort]]$sample_id
    grm_ids  <- normalize_id(rownames(grm_meta[[cohort]]))
    shared_ids <- sort(Reduce(intersect, list(gds_ids, pcs_ids, grm_ids, meth_ids_norm)))

    overlap_rows[[paste(cohort, context)]] <- data.table(
      cohort   = cohort, context = context,
      gds_n    = length(gds_ids),  pcs_n   = length(pcs_ids),
      grm_n    = length(grm_ids),  methyl_n = length(meth_ids_norm),
      shared_n = length(shared_ids)
    )
    log_msg("  Shared samples: ", length(shared_ids),
            " (gds=", length(gds_ids), " pcs=", length(pcs_ids),
            " grm=", length(grm_ids), " meth=", length(meth_ids_norm), ")")

    if (length(shared_ids) < 10L) {
      log_msg("  Too few shared samples — skipping panel")
      summary_rows[[paste(cohort, context)]] <- data.table(
        cohort = cohort, context = context, status = "too_few_samples")
      next
    }

    # --- Restrict matrices ---
    meth_raw  <- rownames(m_mat)
    met_idx   <- match(shared_ids, meth_ids_norm)
    m_shared  <- m_mat[met_idx, , drop = FALSE]
    rownames(m_shared) <- meth_raw[met_idx]

    # Subset PCs to shared samples in the same order
    pcs_shared <- copy(grm_pcs_meta[[cohort]])[
      match(shared_ids, grm_pcs_meta[[cohort]]$sample_id)]

    # Subset GRM to shared samples (symmetric submatrix)
    G      <- grm_meta[[cohort]]
    grm_idx <- match(shared_ids, normalize_id(rownames(G)))
    grm_shared <- G[grm_idx, grm_idx, drop = FALSE]
    rownames(grm_shared) <- colnames(grm_shared) <- shared_ids

    # Group assignments for HWE imputation (used in steps 12ab1/12ab2)
    groups_dt <- annot_meta[[cohort]][match(shared_ids, sample_id)]
    groups_dt[is.na(group), group := "Unknown"]
    groups_dt[, sample_id := shared_ids]

    shared_ids_dt <- data.table(sample_id = shared_ids)

    # --- Write RDS outputs ---
    saveRDS(m_shared,      file.path(panel_dir, "methylation_mvalues_matrix.rds"))
    saveRDS(site_annot,    file.path(panel_dir, "methylation_site_annot.rds"))
    saveRDS(shared_ids_dt, file.path(panel_dir, "shared_sample_ids.rds"))
    saveRDS(pcs_shared,    file.path(panel_dir, "pcs_shared.rds"))
    saveRDS(grm_shared,    file.path(panel_dir, "grm_shared.rds"))
    saveRDS(groups_dt,     file.path(panel_dir, "sample_groups.rds"))
    saveRDS(gds_meta[[cohort]]$sample_ids,    file.path(panel_dir, "snp_gds_sample_ids.rds"))
    saveRDS(gds_meta[[cohort]]$variant_annot, file.path(panel_dir, "snp_variant_annot.rds"))

    # --- Write TSV outputs (human-readable, for inspection) ---
    fwrite(as.data.table(m_shared, keep.rownames = "sample_id_raw"),
           file.path(panel_dir, "methylation_mvalues_matrix.tsv"), sep = "\t")
    fwrite(site_annot,    file.path(panel_dir, "methylation_site_annot.tsv"),  sep = "\t")
    fwrite(shared_ids_dt, file.path(panel_dir, "shared_sample_ids.tsv"),       sep = "\t")
    fwrite(pcs_shared,    file.path(panel_dir, "pcs_shared.tsv"),              sep = "\t")
    fwrite(groups_dt,     file.path(panel_dir, "sample_groups.tsv"),           sep = "\t")
    fwrite(gds_meta[[cohort]]$sample_ids,
           file.path(panel_dir, "snp_gds_sample_ids.tsv"), sep = "\t")
    fwrite(gds_meta[[cohort]]$variant_annot,
           file.path(panel_dir, "snp_variant_annot.tsv"),  sep = "\t")

    summary_rows[[paste(cohort, context)]] <- data.table(
      cohort         = cohort, context = context, status = "ok",
      shared_samples = length(shared_ids),
      methyl_sites   = ncol(m_shared),
      n_pcs          = N_PCS
    )
    log_msg("  Panel inputs written | sites=", ncol(m_shared), " | pcs=", N_PCS)
    rm(m_mat, m_shared, grm_shared, pcs_shared, site_annot, conv, groups_dt); gc()
  }
}

############################################################
# 8) WRITE SUMMARIES
############################################################

sep_line()
log_msg("Writing summaries...")
summary_dt <- rbindlist(summary_rows, fill = TRUE)
overlap_dt <- rbindlist(overlap_rows, fill = TRUE)
fwrite(summary_dt, SUMMARY_FILE, sep = "\t")
fwrite(overlap_dt, OVERLAP_FILE, sep = "\t")
writeLines(capture.output(sessionInfo()), SESSION_FILE)

sep_line()
log_msg("Step 12ab0 finished")
log_msg("Inputs: ",    INPUTDIR)
log_msg("LOCO GRMs: ", loco_dir)
log_msg("Summary: ",   SUMMARY_FILE)
print(summary_dt)
sep_line()
