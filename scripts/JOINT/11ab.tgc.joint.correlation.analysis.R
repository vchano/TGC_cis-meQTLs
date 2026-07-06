#!/usr/bin/env Rscript
############################################################
# TreeGeneClimate (TGC) — JOINT (ECS + TBS)
# Step 11ab: Distance-based congruence between
#            genomic (ECS SNPs) and epigenomic (TBS methylation)
#
# PURPOSE
# - For each cohort (BREEDING, NATURAL) and each methylation context (CpG/CHG/CHH):
#   1) Build a genomic distance matrix from ECS imputed GDS genotypes (IBS distance = 1 - IBS)
#   2) Build an epigenomic distance matrix from TBS methylBase (Euclidean on filtered+imputed % methylation)
#   3) Align samples robustly (trimws + tolower), save diagnostics for mismatches
#   4) Compare distance matrices:
#       - Mantel test
#       - Procrustes / protest (PCoA ordinations)
#       - RV coefficient
#   5) Save per-comparison tables + combined 3x2 Procrustes panel + aggregated summary
#
# INPUTS
# - ECS:
#   RESULTS/ECS/RANALYSIS/RDATA/breeding.imputed.snp.gds
#   RESULTS/ECS/RANALYSIS/RDATA/natural.imputed.snp.gds
# - TBS:
#   RESULTS/TBS/RANALYSIS/METHYLKIT_OBJECTS/
#     methylBase_<cohort>_<context>_cov5_50_mpg*_mef0.05.rds
# - Metadata:
#   DATA/METADATA/breeding_sample2family.txt
#   DATA/METADATA/natural_sample2pop.txt
#
# OUTPUTS
# - /user/chano/u15584/TGC_project/RESULTS/JOINT/
#     FIGURES/11ab/
#     TABLES/11ab/
#
# NOTES
# - PERMUTATIONS fixed at 9999.
# - Run in a fresh R session if GDS files were previously opened.
############################################################

suppressPackageStartupMessages({
  library(SNPRelate)
  library(gdsfmt)
  library(methylKit)
  library(dplyr)
  library(tibble)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(vegan)   # mantel, procrustes, protest
  library(ade4)    # RV coefficient
  library(patchwork)
  library(grid)
})

options(stringsAsFactors = FALSE)
set.seed(1)

############################################################
# 1) FIXED PATHS
############################################################
# === USER CONFIGURATION ===
PROJECT_ROOT <- "/path/to/your/project"  # <-- set this
# ===========================

# ECS (GDS produced by 8a)
ECS_RDATA_DIR <- file.path(PROJECT_ROOT, "RESULTS/ECS/RANALYSIS/RDATA")
GDS_BREED <- file.path(ECS_RDATA_DIR, "breeding.imputed.snp.gds")
GDS_NATUR <- file.path(ECS_RDATA_DIR, "natural.imputed.snp.gds")

# TBS methylKit objects
TBS_RDS_DIR <- file.path(PROJECT_ROOT, "RESULTS/TBS/RANALYSIS/METHYLKIT_OBJECTS")

# Metadata maps
META_DIR <- file.path(PROJECT_ROOT, "DATA/METADATA")
MAP_BREED <- file.path(META_DIR, "breeding_sample2family.txt")
MAP_NATUR <- file.path(META_DIR, "natural_sample2pop.txt")

# JOINT outputs
JOINT_ROOT <- "/user/chano/u15584/TGC_project/RESULTS/JOINT"
OUT_FIG <- file.path(JOINT_ROOT, "FIGURES", "11ab")
OUT_TAB <- file.path(JOINT_ROOT, "TABLES",  "11ab")
OUT_DIAG <- file.path(OUT_TAB, "DIAGNOSE")

dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIAG, recursive = TRUE, showWarnings = FALSE)

############################################################
# 2) PARAMETERS
############################################################
PERMUTATIONS <- 9999
# Standard SNP QC thresholds used throughout the ECS pipeline
ECS_MAF_MIN  <- 0.05
ECS_MISS_MAX <- 0.10

# Context-specific NA tolerance: CHH sites have higher missing rates because
# non-CpG methylation is sparse and more variable across individuals
max_na_frac_by_ctx <- function(ctx) {
  switch(ctx,
         "CpG" = 0.20,
         "CHG" = 0.30,
         "CHH" = 0.50,
         0.30)
}

############################################################
# 3) HELPERS
############################################################
ensure_file <- function(p) if (!file.exists(p)) stop("Missing file: ", p, call. = FALSE)

# Normalise sample IDs to a canonical form for cross-dataset matching
norm_id <- function(x) tolower(trimws(as.character(x)))

# Select the TBS RDS with the highest mpg threshold (most stringent coverage filter)
pick_tbs_mef_rds <- function(cohort, ctx) {
  pat <- sprintf("^methylBase_%s_%s_cov5_50_mpg[0-9]+_mef0\\.05\\.rds$",
                 tolower(cohort), tolower(ctx))
  files <- list.files(TBS_RDS_DIR, pattern = pat, full.names = TRUE)
  if (!length(files)) stop("No TBS MEF RDS found for ", cohort, " / ", ctx, " under ", TBS_RDS_DIR)
  mpg_num <- as.integer(sub(".*_mpg([0-9]+)_mef0\\.05\\.rds$", "\\1", basename(files)))
  files[which.max(mpg_num)]
}

# ---- ECS genomic distance: IBS -> distance (1-IBS) ----
ecs_distance <- function(gds) {

  g <- snpgdsOpen(gds)
  on.exit(try(snpgdsClose(g), silent = TRUE), add = TRUE)

  stat <- snpgdsSNPRateFreq(g, with.id = TRUE)

  maf  <- stat$MinorFreq
  miss <- stat$MissingRate

  # Apply MAF and missingness QC before computing IBS
  keep <- is.finite(maf) & is.finite(miss) & maf >= 0.05 & miss <= 0.10

  if (sum(keep) == 0) stop("No ECS SNPs passed MAF/missingness filters.")

  ibs <- snpgdsIBS(
    g,
    snp.id = stat$snp.id[keep],
    num.thread = 8,
    autosome.only = FALSE  # spruce has no canonical autosomes in the reference
  )

  M <- ibs$ibs
  rownames(M) <- colnames(M) <- ibs$sample.id

  # Convert IBS similarity to a distance matrix
  D <- as.dist(1 - M)

  list(
    dist = D,
    ids = ibs$sample.id,
    n_snps = sum(keep)
  )
}

# ---- TBS epigenomic distance: Euclidean on %methylation (sites filtered+imputed) ----
methylbase_to_percent_matrix <- function(mb) {
  d <- getData(mb)
  numCs_cols <- grep("^numCs[0-9]+$", colnames(d), value = TRUE)
  numTs_cols <- grep("^numTs[0-9]+$", colnames(d), value = TRUE)
  if (length(numCs_cols) == 0 || length(numTs_cols) == 0 || length(numCs_cols) != length(numTs_cols)) {
    stop("Could not find matching numCs/numTs columns in methylBase.")
  }
  sample_ids <- mb@sample.ids
  if (is.null(sample_ids) || length(sample_ids) != length(numCs_cols)) {
    sample_ids <- paste0("S", seq_along(numCs_cols))
  }

  # Compute per-site percentage methylation for each sample
  perc_mat <- vapply(seq_along(numCs_cols), function(i) {
    numCs <- d[[numCs_cols[i]]]
    numTs <- d[[numTs_cols[i]]]
    p <- 100 * (numCs / (numCs + numTs))
    p[is.nan(p)] <- NA_real_  # sites with zero coverage yield NaN
    p
  }, numeric(nrow(d)))

  colnames(perc_mat) <- sample_ids
  perc_mat
}

filter_and_impute_tbs <- function(X_samples_x_sites, max_na_frac) {
  # Remove sites that are missing in too many samples
  na_frac <- colMeans(is.na(X_samples_x_sites))
  keep1 <- na_frac <= max_na_frac
  X2 <- X_samples_x_sites[, keep1, drop = FALSE]
  if (ncol(X2) < 10) stop("Too few TBS sites after NA filter (n=", ncol(X2), ").")

  # Mean-impute remaining NAs within each site (column-wise)
  for (j in seq_len(ncol(X2))) {
    idx <- is.na(X2[, j])
    if (any(idx)) {
      mu <- mean(X2[, j], na.rm = TRUE)
      if (!is.finite(mu)) mu <- 0
      X2[idx, j] <- mu
    }
  }

  # Drop invariant sites — they contribute nothing to Euclidean distance
  sds <- apply(X2, 2, sd)
  keep2 <- is.finite(sds) & sds > 0
  X3 <- X2[, keep2, drop = FALSE]
  if (ncol(X3) < 10) stop("Too few variable TBS sites after SD filter (n=", ncol(X3), ").")

  X3
}

compute_tbs_euclid_dist <- function(tbs_mef_rds, ctx) {
  ensure_file(tbs_mef_rds)
  mb <- readRDS(tbs_mef_rds)

  perc_sites_x_samples <- methylbase_to_percent_matrix(mb)
  # Transpose so rows = samples, columns = sites
  X <- t(perc_sites_x_samples)
  rownames(X) <- colnames(perc_sites_x_samples)

  max_na <- max_na_frac_by_ctx(ctx)
  X2 <- filter_and_impute_tbs(X, max_na_frac = max_na)
  # Euclidean distance across all filtered CpG/CHG/CHH sites (high-dimensional space)
  D <- dist(X2, method = "euclidean")
  list(dist = D, ids = rownames(X2), n_sites = ncol(X2))
}

# ---- Align by normalized IDs and save diagnostics ----
align_dist_by_ids <- function(D1, ids1, D2, ids2, out_diag_prefix) {
  n1 <- norm_id(ids1)
  n2 <- norm_id(ids2)

  # Detect duplicate IDs within each dataset before intersection
  dup1 <- duplicated(n1)
  dup2 <- duplicated(n2)
  if (any(dup1) || any(dup2)) {
    write_tsv(
      tibble(
        side = c(rep("ECS", sum(dup1)), rep("TBS", sum(dup2))),
        id   = c(ids1[dup1], ids2[dup2]),
        norm = c(n1[dup1], n2[dup2])
      ),
      file.path(OUT_DIAG, paste0(out_diag_prefix, "_duplicate_norm_ids.tsv"))
    )
    stop("Duplicate normalized IDs detected. See diagnostics: ", out_diag_prefix, "_duplicate_norm_ids.tsv")
  }

  map1 <- setNames(ids1, n1)
  map2 <- setNames(ids2, n2)
  common_norm <- intersect(names(map1), names(map2))

  # Log which samples are present in only one dataset for QC review
  only_1 <- setdiff(names(map1), common_norm)
  only_2 <- setdiff(names(map2), common_norm)

  write_tsv(tibble(side = "ECS_only", norm_id = only_1, example_id = map1[only_1]),
            file.path(OUT_DIAG, paste0(out_diag_prefix, "_missing_ECSonly.tsv")))
  write_tsv(tibble(side = "TBS_only", norm_id = only_2, example_id = map2[only_2]),
            file.path(OUT_DIAG, paste0(out_diag_prefix, "_missing_TBSonly.tsv")))

  if (length(common_norm) < 10) stop("Too few overlapping samples after ID normalization: ", length(common_norm))

  ids1_keep <- unname(map1[common_norm])
  ids2_keep <- unname(map2[common_norm])

  # Subset and re-order both distance matrices to the common sample set
  M1 <- as.matrix(D1); rownames(M1) <- colnames(M1) <- ids1
  M2 <- as.matrix(D2); rownames(M2) <- colnames(M2) <- ids2

  M1s <- M1[ids1_keep, ids1_keep, drop = FALSE]
  M2s <- M2[ids2_keep, ids2_keep, drop = FALSE]

  # Use normalised IDs as the canonical row/col names for both matrices
  ord <- common_norm
  rownames(M1s) <- colnames(M1s) <- ord
  rownames(M2s) <- colnames(M2s) <- ord

  list(D1 = as.dist(M1s), D2 = as.dist(M2s), common_norm = ord, n = length(ord))
}

# ---- Ordination (2 axes) ----
# PCoA (= classical MDS) on a distance matrix; add=TRUE corrects negative eigenvalues
pcoa2 <- function(D) {
  m <- cmdscale(D, k = 2, eig = TRUE, add = TRUE)
  X <- as.data.frame(m$points)
  colnames(X) <- c("Axis1", "Axis2")
  X
}

# ---- Procrustes plot object for panel ----
make_procrustes_plot <- function(X, Y, title_label = "a)") {
  # Symmetric Procrustes: both configurations are scaled and rotated optimally
  pro <- vegan::procrustes(X, Y, symmetric = TRUE)
  Yrot <- pro$Yrot

  df <- tibble(
    id = rownames(X),
    X1 = X[,1], X2 = X[,2],
    Y1 = Yrot[,1], Y2 = Yrot[,2]
  )

  # Segments connect each individual's genomic (SNP) position to its epigenomic (SMP) position
  p <- ggplot(df) +
    geom_segment(aes(x = Y1, y = Y2, xend = X1, yend = X2), alpha = 0.45, color = "grey50") +
    geom_point(aes(x = X1, y = X2, color = "SNPs"), size = 1.7) +
    geom_point(aes(x = Y1, y = Y2, color = "SMPs"), size = 1.7, alpha = 0.9) +
    scale_color_manual(values = c("SNPs" = "#4e79a7", "SMPs" = "#f28e2b"), name = NULL) + guides(color = guide_legend(override.aes = list(size = 6))) +
    labs(x = "PC1", y = "PC2") +
    annotate("text", x = -Inf, y = Inf, label = title_label,
             hjust = -0.2, vjust = 1.3, size = 7) +
    theme_minimal(base_size = 18) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.direction = "horizontal",
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 12),
      legend.text = element_text(size = 16),
      legend.key.size = unit(4, "lines")
    )

  list(plot = p, pro = pro)
}

############################################################
# 4) SANITY CHECKS
############################################################
ensure_file(GDS_BREED); ensure_file(GDS_NATUR)
ensure_file(MAP_BREED); ensure_file(MAP_NATUR)

############################################################
# 5) RUN: 2 cohorts × 3 contexts
############################################################
cohorts <- list(
  BREEDING = list(gds = GDS_BREED, map = MAP_BREED),
  NATURAL  = list(gds = GDS_NATUR, map = MAP_NATUR)
)
contexts <- c("CpG", "CHG", "CHH")

plot_list <- list()
summary_rows <- list()

# Panel labels follow publication convention (a–f, top-left to bottom-right)
panel_labels <- c("a)", "b)", "c)", "d)", "e)", "f)")
panel_i <- 1

for (coh in names(cohorts)) {
  message("\n==============================")
  message("COHORT: ", coh)
  message("==============================")

  # Compute IBS distance once per cohort (shared across all three contexts)
  ecs <- ecs_distance(cohorts[[coh]]$gds)
  message("ECS: IBS distance computed (SNPs passing QC = ", ecs$n_snps, ")")

  for (ctx in contexts) {
    message("\n--- ", coh, " / ", ctx, " ---")

    tbs_rds <- pick_tbs_mef_rds(coh, ctx)
    tbs <- compute_tbs_euclid_dist(tbs_rds, ctx)
    message("TBS: Euclidean distance computed (sites used = ", tbs$n_sites, ")")

    prefix <- paste0("11ab_", tolower(coh), "_", tolower(ctx))
    al <- align_dist_by_ids(ecs$dist, ecs$ids, tbs$dist, tbs$ids, out_diag_prefix = prefix)
    message("Aligned samples: n = ", al$n)

    # Mantel test: correlation between the two distance matrices (permutation-based)
    man <- vegan::mantel(al$D1, al$D2, method = "pearson", permutations = PERMUTATIONS)

    # PCoA of each distance matrix; Procrustes / protest then compares the ordinations
    X <- pcoa2(al$D1)
    Y <- pcoa2(al$D2)
    rownames(X) <- rownames(Y) <- al$common_norm

    pro_test <- vegan::protest(X, Y, permutations = PERMUTATIONS)
    t0 <- unname(pro_test$t0)   # Procrustes statistic (lower = better alignment)
    p_pro <- unname(pro_test$signif)

    # RV coefficient: multivariate analogue of the squared correlation
    rv_res <- ade4::RV.rtest(as.data.frame(X), as.data.frame(Y), nrepet = PERMUTATIONS)
    rv <- unname(rv_res$obs)
    p_rv <- unname(rv_res$pvalue)

    row <- tibble(
      cohort = coh,
      context = ctx,
      n_samples = al$n,
      ecs_snps_qc = ecs$n_snps,
      tbs_sites_used = tbs$n_sites,
      mantel_r = unname(man$statistic),
      mantel_p = unname(man$signif),
      procrustes_t0 = t0,
      procrustes_p = p_pro,
      rv = rv,
      rv_p = p_rv,
      tbs_rds = tbs_rds
    )

    out_sum <- file.path(OUT_TAB, paste0(prefix, "_summary.tsv"))
    write_tsv(row, out_sum)

    pro_obj <- make_procrustes_plot(X, Y, title_label = panel_labels[panel_i])
    plot_list[[panel_i]] <- pro_obj$plot

    # Per-sample Procrustes residuals indicate how well each individual is matched
    Yrot <- pro_obj$pro$Yrot
    residuals <- sqrt(rowSums((as.matrix(X) - as.matrix(Yrot))^2))
    out_res <- file.path(OUT_TAB, paste0(prefix, "_procrustes_residuals.tsv"))
    write_tsv(tibble(sample_norm = al$common_norm, residual = residuals), out_res)

    summary_rows[[prefix]] <- row
    message("Saved: ", basename(out_sum), " | ", basename(out_res))

    panel_i <- panel_i + 1
  }
}

############################################################
# 6) COMBINED 3x2 PROCRUSTES PANEL WITH COMMON LEGEND
############################################################
# Constrain y-axis for NATURAL panels (d–f) which have tighter ordination spread
plot_list[[4]] <- plot_list[[4]] + coord_cartesian(ylim = c(-0.05, 0.05))
plot_list[[5]] <- plot_list[[5]] + coord_cartesian(ylim = c(-0.05, 0.05))
plot_list[[6]] <- plot_list[[6]] + coord_cartesian(ylim = c(-0.05, 0.05))

# Row 1: BREEDING (CpG / CHG / CHH); Row 2: NATURAL (CpG / CHG / CHH)
combined_panel <-
  ((plot_list[[1]] | plot_list[[2]] | plot_list[[3]]) /
     (plot_list[[4]] | plot_list[[5]] | plot_list[[6]])) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal"
  )

out_panel <- file.path(OUT_FIG, "Figure11_ECS_TBS_Procrustes_panel.tiff")

# Save in all required publication formats
tiff(out_panel, width = 34, height = 24, units = "cm", res = 600, compression = "lzw")
print(combined_panel)
dev.off()
ggsave(sub("\\.tiff$", ".pdf", out_panel), plot = combined_panel, width = 34, height = 24, units = "cm")
ggsave(sub("\\.tiff$", ".eps", out_panel), plot = combined_panel, width = 34, height = 24, units = "cm", device = cairo_ps)
ggsave(sub("\\.tiff$", ".png", out_panel), plot = combined_panel, width = 34, height = 24, units = "cm", dpi = 150, device = "png")

############################################################
# 7) AGGREGATED SUMMARY + BH CORRECTION
############################################################
# Apply BH correction across all 6 panels jointly for each test type
all_df <- bind_rows(summary_rows) %>%
  mutate(
    mantel_p_adj = p.adjust(mantel_p, method = "BH"),
    procrustes_p_adj = p.adjust(procrustes_p, method = "BH"),
    rv_p_adj = p.adjust(rv_p, method = "BH")
  ) %>%
  arrange(cohort, context)

out_all <- file.path(OUT_TAB, "11ab_ecs_tbs_distance_summary_all.tsv")
write_tsv(all_df, out_all)

message("\nDONE 11ab.")
message("Combined figure: ", out_panel)
message("Tables: ", OUT_TAB)
message("Master summary: ", out_all)
sessionInfo()
