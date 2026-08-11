#!/usr/bin/env Rscript
############################################################
# TreeGeneClimate (TGC) — JOINT ECS + TBS
# Step 16ab7: Per-site methylation heritability
#
# Estimates SNP-based heritability (h²_SNP) at every
# MEF-retained cytosine site using per-site GREML via
# eigendecomposition of the GRM (EMMA; Kang et al. 2008,
# Genetics 178:1709-1723). The GRM eigendecomposition is
# performed once per cohort; per-site REML is a fast
# one-dimensional optimisation over the variance ratio.
#
# For the breeding cohort, pedigree heritability (h²_ped)
# is additionally estimated as 2 × ICC from a full-sib
# family random-effects model (lme4 REML).
#
# INPUTS
#   <RDATA_DIR>/<cohort>_grm_gcta.rds         — dense GRM (n × n)
#   RESULTS/JOINT/MQTL5/INPUTS/<C>/<X>/
#     methylation_mvalues_matrix.rds            — M-value matrix
#   ECS/breeding_sample2family.txt             — sample → family map
#   RESULTS/JOINT/COMBINED5/overlap/tables/
#     robust_markers_{breeding,natural}.tsv    — meQTL site lists
#
# OUTPUTS (RESULTS/JOINT/HERITABILITY/)
#   h2_all_sites.tsv.gz      — per-site h²_SNP [and h²_ped] results
#   h2_summary.tsv           — aggregate stats per cohort × context
#   h2_meqtl_vs_all.*        — figure: meQTL sites vs genome background
#   h2_grm_vs_ped.*          — figure: h²_SNP vs h²_ped (breeding)
#   h2_combined_panel.*      — combined A/B panel figure
#
# NOTE: This script runs sequentially (all sites in memory per
# cohort × context). For very large datasets (>100k sites) use
# the companion SLURM script to distribute chunks across nodes;
# combine chunk TSVs, then skip to section 2 (aggregation).
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

options(stringsAsFactors = FALSE)

# === USER CONFIGURATION ===
PROJECT_ROOT <- "/path/to/your/project"   # <-- set this
RDATA_DIR    <- "/path/to/your/rdata"     # directory with <cohort>_grm_gcta.rds
# ===========================

COHORTS  <- c("BREEDING", "NATURAL")
CONTEXTS <- c("CpG", "CHG", "CHH")

H2_ROOT     <- file.path(PROJECT_ROOT, "RESULTS/JOINT/HERITABILITY")
ROBUST_FILE <- file.path(PROJECT_ROOT,
  "RESULTS/JOINT/COMBINED5/overlap/tables/robust_markers_breeding.tsv")
ROBUST_NAT  <- file.path(PROJECT_ROOT,
  "RESULTS/JOINT/COMBINED5/overlap/tables/robust_markers_natural.tsv")
FIG_DIR     <- file.path(PROJECT_ROOT, "RESULTS/DRAFT")

dir.create(H2_ROOT, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

msg <- function(...) {
  cat(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n"))
  flush(stdout())
}

############################################################
# CORE FUNCTIONS
############################################################

# Per-site REML via EMMA eigendecomposition (Kang et al. 2008)
# U, D_pos: pre-computed eigenvectors / positive eigenvalues of GRM
# Returns: h²_SNP = Vg / (Vg + Ve)
emma_reml <- function(y, U, D_pos) {
  na_idx <- is.na(y)
  if (any(na_idx)) y[na_idx] <- mean(y, na.rm = TRUE)
  y  <- y - mean(y)
  Uy <- drop(crossprod(U, y))
  n  <- length(D_pos)
  obj <- function(ld) {
    d  <- D_pos + exp(ld)
    vg <- sum(Uy^2 / d) / n
    if (vg <= 0) return(1e15)
    n * log(vg) + sum(log(d))
  }
  opt <- tryCatch(optimize(obj, c(-20, 20), tol = 1e-8), error = function(e) NULL)
  if (is.null(opt)) return(c(h2 = NA_real_, Vg = NA_real_, Ve = NA_real_))
  delta <- exp(opt$minimum)
  d     <- D_pos + delta
  Vg    <- sum(Uy^2 / d) / n
  Ve    <- Vg * delta
  Vp    <- Vg + Ve
  c(h2 = if (Vp > 1e-15) Vg / Vp else NA_real_, Vg = Vg, Ve = Ve)
}

# Pedigree h² via full-sib family random effect
# h²_ped ≈ 2 × ICC (assumes dominance and shared environment ≈ 0)
pedigree_h2 <- function(y, family_vec) {
  na_idx <- is.na(y)
  if (any(na_idx)) y[na_idx] <- mean(y, na.rm = TRUE)
  df  <- data.frame(y = y, fam = factor(family_vec))
  fit <- tryCatch(
    suppressWarnings(suppressMessages(
      lme4::lmer(y ~ 1 + (1 | fam), data = df, REML = TRUE,
                 control = lme4::lmerControl(
                   optimizer = "nloptwrap",
                   optCtrl   = list(maxfun = 2e5),
                   check.conv.singular = "ignore")))),
    error = function(e) NULL)
  if (is.null(fit))
    return(c(h2_ped = NA_real_, ICC = NA_real_, Vf = NA_real_, Ve_ped = NA_real_))
  vc  <- as.data.frame(lme4::VarCorr(fit))
  Vf  <- vc$vcov[vc$grp == "fam"]
  Ve2 <- vc$vcov[vc$grp == "Residual"]
  if (length(Vf) == 0) Vf <- 0
  Vp  <- Vf + Ve2
  ICC <- if (Vp > 1e-15) Vf / Vp else NA_real_
  c(h2_ped = if (!is.na(ICC)) 2 * ICC else NA_real_,
    ICC = ICC, Vf = Vf, Ve_ped = Ve2)
}

############################################################
# 1) PER-SITE HERITABILITY
############################################################

all_list <- list()

for (cohort in COHORTS) {
  DO_PEDIGREE <- (cohort == "BREEDING")

  msg("Loading GRM: ", cohort)
  K <- as.matrix(readRDS(
    file.path(RDATA_DIR, sprintf("%s_grm_gcta.rds", tolower(cohort)))))

  fam_dt <- NULL
  if (DO_PEDIGREE) {
    fam_file <- file.path(PROJECT_ROOT, "ECS/breeding_sample2family.txt")
    fam_dt   <- fread(fam_file, header = FALSE, col.names = c("sample", "family"))
  }

  for (ctx in CONTEXTS) {
    mval_file <- file.path(PROJECT_ROOT, "RESULTS/JOINT/MQTL5/INPUTS",
                           cohort, ctx, "methylation_mvalues_matrix.rds")
    if (!file.exists(mval_file)) {
      msg("  MISSING: ", mval_file, " — skipping"); next
    }

    Mmat   <- as.matrix(readRDS(mval_file))
    common <- intersect(rownames(K), rownames(Mmat))
    K_sub  <- K[common, common, drop = FALSE]
    M_sub  <- Mmat[common, , drop = FALSE]
    n      <- nrow(K_sub)
    msg("  ", cohort, "/", ctx, ": ", n, " samples | ", ncol(M_sub), " sites")

    eig   <- eigen(K_sub, symmetric = TRUE)
    pos   <- eig$values > 1e-10
    U     <- eig$vectors[, pos, drop = FALSE]
    D_pos <- eig$values[pos]
    rm(eig); gc()

    fam_vec <- if (DO_PEDIGREE)
      factor(fam_dt$family[match(common, fam_dt$sample)]) else NULL

    n_sites <- ncol(M_sub)
    res_lst <- vector("list", n_sites)
    for (i in seq_len(n_sites)) {
      y          <- M_sub[, i]
      gr_res     <- emma_reml(y, U, D_pos)
      res_lst[[i]] <- if (DO_PEDIGREE) c(gr_res, pedigree_h2(y, fam_vec)) else gr_res
    }

    dt <- as.data.table(do.call(rbind, res_lst))
    dt[, `:=`(site = colnames(M_sub), cohort = cohort, context = ctx)]
    all_list[[paste(cohort, ctx)]] <- dt
    rm(M_sub, K_sub, U, D_pos, res_lst); gc()
  }
  rm(K); gc()
}

h2_all <- rbindlist(all_list, use.names = TRUE, fill = TRUE)
setorder(h2_all, cohort, context, site)
out_gz <- file.path(H2_ROOT, "h2_all_sites.tsv.gz")
fwrite(h2_all, out_gz, sep = "\t", compress = "gzip")
msg("Saved: h2_all_sites.tsv.gz (", nrow(h2_all), " rows)")

############################################################
# 2) SUMMARY TABLE
############################################################

summary_dt <- h2_all[, .(
  n_sites      = .N,
  n_valid      = sum(!is.na(h2)),
  mean_h2      = mean(h2, na.rm = TRUE),
  median_h2    = median(h2, na.rm = TRUE),
  sd_h2        = sd(h2, na.rm = TRUE),
  pct_h2_gt0.3 = mean(h2 > 0.30, na.rm = TRUE) * 100,
  pct_h2_gt0.5 = mean(h2 > 0.50, na.rm = TRUE) * 100,
  mean_h2_ped  = if ("h2_ped" %in% names(.SD))
                   mean(h2_ped, na.rm = TRUE) else NA_real_
), by = .(cohort, context)]
setorder(summary_dt, cohort, context)
fwrite(summary_dt, file.path(H2_ROOT, "h2_summary.tsv"), sep = "\t")
print(summary_dt[, .(cohort, context, n_valid, mean_h2, pct_h2_gt0.3, mean_h2_ped)])

############################################################
# 3) FIGURES
############################################################

meqtl_sites <- character(0)
for (rf in c(ROBUST_FILE, ROBUST_NAT)) {
  if (file.exists(rf))
    meqtl_sites <- union(meqtl_sites, unique(fread(rf, select = "site")$site))
}
msg("Robust meQTL sites: ", length(meqtl_sites))
h2_all[, is_meqtl := site %in% meqtl_sites]

save_fig <- function(base, plot, width, height) {
  ggsave(paste0(base, ".tiff"), plot, width = width, height = height,
         dpi = 300, compression = "lzw")
  ggsave(paste0(base, ".pdf"),  plot, width = width, height = height)
  ggsave(paste0(base, ".png"),  plot, width = width, height = height, dpi = 300)
  tryCatch(ggsave(paste0(base, ".eps"), plot, width = width, height = height),
           error = function(e) msg("EPS skipped: ", conditionMessage(e)))
}

ctx_cols   <- c(CpG = "#E63946", CHG = "#457B9D", CHH = "#2A9D8F")
coh_labels <- c(BREEDING = "Breeding (n=209)", NATURAL = "Natural (n=393)")

h2_plot <- h2_all[!is.na(h2)]
h2_plot[, context_f := factor(context, levels = c("CpG", "CHG", "CHH"))]
h2_plot[, cohort_f  := factor(cohort, levels = c("BREEDING", "NATURAL"),
                               labels = coh_labels)]

# Panel A: h²_SNP — meQTL sites vs genome background
h2_meqtl <- h2_plot[, .(
  group = ifelse(is_meqtl, "Robust meQTL sites", "All sites"),
  h2, context_f, cohort_f)]

p_meqtl <- ggplot(h2_meqtl, aes(x = context_f, y = h2, fill = group, colour = group)) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.5,
              position = position_dodge(0.8), linewidth = 0.3) +
  geom_boxplot(width = 0.06, outlier.shape = NA, fill = NA,
               position = position_dodge(0.8), linewidth = 0.4, coef = 0) +
  facet_wrap(~cohort_f, nrow = 1) +
  scale_fill_manual(values   = c("All sites" = "grey60",
                                  "Robust meQTL sites" = "#E63946"), name = NULL) +
  scale_colour_manual(values = c("All sites" = "grey40",
                                  "Robust meQTL sites" = "#9B1A22"), name = NULL) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(x = "Methylation context", y = expression(h[SNP]^2)) +
  theme_bw(base_size = 11) +
  theme(strip.background = element_rect(fill = "grey90"),
        legend.position  = "bottom",
        panel.grid.minor = element_blank())
save_fig(file.path(H2_ROOT, "h2_meqtl_vs_all"), p_meqtl, width = 8, height = 4.5)

# Panel B: h²_SNP vs h²_ped (breeding cohort only)
breed_h2 <- h2_all[cohort == "BREEDING" & !is.na(h2) & !is.na(h2_ped)]
if (nrow(breed_h2) > 0) {
  breed_h2[, context_f := factor(context, levels = c("CpG", "CHG", "CHH"))]
  set.seed(42)
  breed_plot <- breed_h2[, .SD[sample(.N, min(.N, 50000L))], by = context_f]

  p_ped <- ggplot(breed_plot, aes(x = h2, y = h2_ped, colour = context_f)) +
    geom_point(alpha = 0.08, size = 0.4, stroke = 0) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                colour = "black", linewidth = 0.5) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
    facet_wrap(~context_f, nrow = 1) +
    scale_colour_manual(values = ctx_cols, guide = "none") +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(x = expression(h[SNP]^2 * " (GREML / GRM)"),
         y = expression(h[Ped]^2 * " (2 × ICC, full-sib families)")) +
    theme_bw(base_size = 11) +
    theme(strip.background = element_rect(fill = "grey90"),
          panel.grid.minor = element_blank())
  save_fig(file.path(H2_ROOT, "h2_grm_vs_ped"), p_ped, width = 8, height = 3.5)

  # Combined two-panel figure
  combined <- (p_meqtl + labs(tag = "A") +
                 theme(plot.tag = element_text(face = "bold", size = 14))) /
              (p_ped   + labs(tag = "B") +
                 theme(plot.tag = element_text(face = "bold", size = 14))) +
    plot_layout(heights = c(4.5, 3.5))
  save_fig(file.path(FIG_DIR, "h2_combined_panel"), combined, width = 8, height = 8)
  msg("Saved: h2_combined_panel.*")
} else {
  msg("Skipping GRM-vs-pedigree figure (no breeding pedigree data)")
}

msg("Step 16ab7 complete. Outputs in: ", H2_ROOT)
