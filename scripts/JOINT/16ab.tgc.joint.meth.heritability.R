#!/usr/bin/env Rscript
############################################################
# TreeGeneClimate (TGC) — JOINT ECS + TMS
# Step 16ab: Per-site heritability of DNA methylation and whole-methylome
#            Mantel test (GRM kinship vs methylation similarity)
#
# For every methylation site entering cis-meQTL mapping (Step 12ab0):
#   - h2_SNP: REML via eigendecomposition of the GRM (EMMA)
#             on M-values; h2 = Vg / (Vg + Ve)
#   - h2_Ped (breeding only): 2 x ICC from lmer(M ~ 1 + (1 | family))
# Whole-methylome Mantel test per cohort x context: Pearson correlation between
#   pairwise GRM kinship and pairwise methylation similarity 1 - d_ij / max(d)
#   (Euclidean distance over all sites entering cis-meQTL mapping, computed
#   for each pair of trees over the sites observed in both, without imputing
#   missing values; trees observed at < 20% of the sites left out), 9,999
#   symmetric permutations of the kinship matrix.
# GRM: GCTA-estimator GRM from called genotypes (Step 8a), the same used for
#   mapping.
#
# MODES
#   chunk   : TGC_MODE=chunk TGC_COHORT=BREEDING TGC_CONTEXT=CHH TGC_CHUNK=7
#             (sites (TGC_CHUNK-1)*CHUNK_SIZE+1 ... TGC_CHUNK*CHUNK_SIZE)
#   collect : TGC_MODE=collect  (merges chunks, summary, heritability panel, Mantel test)
#
# INPUTS
#   RESULTS/JOINT/MEQTL/INPUTS/<COHORT>/<CTX>/methylation_mvalues_matrix.rds
#   RESULTS/JOINT/MEQTL/COMBINED/robust_pairs_<cohort>.tsv     (Step 13ab)
#   RESULTS/ECS/RANALYSIS/RDATA/<cohort>_grm_gcta.rds           (Step 8a)
#   DATA/METADATA/breeding_sample2family.txt
#
# OUTPUTS (RESULTS/JOINT/HERITABILITY/)
#   CHUNKS/<COHORT>_<CTX>_<k>.tsv.gz
#   h2_all_sites.tsv.gz, h2_summary.tsv, h2_meqtl_vs_background.tsv
#   mantel_grm_methylation.tsv
#   FIGURES/heritability_panel.{tiff,pdf,png,eps}
#   FIGURES/heritability_panel_h2_values.tsv      h2_SNP of every site, robust cis-meQTL or other (panel a)
#   FIGURES/heritability_panel_h2_ped_values.tsv  h2_SNP and h2_Ped of the breeding sites shown (panel b)
#
# R packages: data.table, lme4, ggplot2, patchwork
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})
options(stringsAsFactors = FALSE)

# === USER CONFIGURATION ===
PROJECT_ROOT <- Sys.getenv("TGC_PROJECT_ROOT", unset = "/path/to/your/project")
CHUNK_SIZE <- 2000L
N_PERM     <- 9999L
MIN_OBS_MANTEL <- 0.20   # Mantel test: trees observed at < 20% of sites are left out
# ===========================

MODE <- Sys.getenv("TGC_MODE", unset = "collect")
COHORTS  <- c("BREEDING", "NATURAL")
CONTEXTS <- c("CpG", "CHG", "CHH")
IN_ROOT  <- file.path(PROJECT_ROOT, "RESULTS/JOINT/MEQTL/INPUTS")
ROB_DIR  <- file.path(PROJECT_ROOT, "RESULTS/JOINT/MEQTL/COMBINED")
RDATA    <- file.path(PROJECT_ROOT, "RESULTS/ECS/RANALYSIS/RDATA")
H2_ROOT  <- file.path(PROJECT_ROOT, "RESULTS/JOINT/HERITABILITY")
FIG_DIR  <- file.path(H2_ROOT, "FIGURES")
for (d in c(file.path(H2_ROOT, "CHUNKS"), FIG_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
msg <- function(...) { cat(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n")); flush(stdout()) }

load_inputs <- function(cohort, ctx) {
  M <- as.matrix(readRDS(file.path(IN_ROOT, cohort, ctx, "methylation_mvalues_matrix.rds")))
  K <- as.matrix(readRDS(file.path(RDATA, sprintf("%s_grm_gcta.rds", tolower(cohort)))))
  common <- intersect(rownames(M), rownames(K))
  list(M = M[common, , drop = FALSE], K = K[common, common, drop = FALSE])
}

# Per-site REML via eigendecomposition (EMMA)
emma_reml <- function(y, U, D_pos) {
  y[is.na(y)] <- mean(y, na.rm = TRUE)
  y  <- y - mean(y)
  Uy <- drop(crossprod(U, y))
  n  <- length(D_pos)
  obj <- function(ld) { d <- D_pos + exp(ld); vg <- sum(Uy^2 / d) / n
                        if (vg <= 0) return(1e15); n * log(vg) + sum(log(d)) }
  opt <- tryCatch(optimize(obj, c(-20, 20), tol = 1e-8), error = function(e) NULL)
  if (is.null(opt)) return(c(h2 = NA_real_, Vg = NA_real_, Ve = NA_real_))
  delta <- exp(opt$minimum); d <- D_pos + delta
  Vg <- sum(Uy^2 / d) / n; Ve <- Vg * delta
  c(h2 = if (Vg + Ve > 1e-15) Vg / (Vg + Ve) else NA_real_, Vg = Vg, Ve = Ve)
}

# Pedigree h2 = 2 x ICC from a full-sib family random effect
pedigree_h2 <- function(y, fam) {
  y[is.na(y)] <- mean(y, na.rm = TRUE)
  fit <- tryCatch(suppressWarnings(suppressMessages(
    lme4::lmer(y ~ 1 + (1 | fam), data = data.frame(y = y, fam = fam), REML = TRUE,
               control = lme4::lmerControl(optimizer = "nloptwrap", optCtrl = list(maxfun = 2e5),
                                           check.conv.singular = "ignore")))),
    error = function(e) NULL)
  if (is.null(fit)) return(c(h2_ped = NA_real_, ICC = NA_real_))
  vc <- as.data.frame(lme4::VarCorr(fit))
  Vf <- vc$vcov[vc$grp == "fam"]; if (!length(Vf)) Vf <- 0
  Ve <- vc$vcov[vc$grp == "Residual"]
  ICC <- if (Vf + Ve > 1e-15) Vf / (Vf + Ve) else NA_real_
  c(h2_ped = 2 * ICC, ICC = ICC)
}

############################################################
# CHUNK MODE
############################################################
if (MODE == "chunk") {
  cohort <- Sys.getenv("TGC_COHORT"); ctx <- Sys.getenv("TGC_CONTEXT")
  k <- as.integer(Sys.getenv("TGC_CHUNK"))
  inp <- load_inputs(cohort, ctx)
  idx <- seq.int((k - 1L) * CHUNK_SIZE + 1L, min(k * CHUNK_SIZE, ncol(inp$M)))
  if (idx[1] > ncol(inp$M)) { msg("chunk ", k, " beyond ", ncol(inp$M), " sites"); quit(save = "no") }
  eig <- eigen(inp$K, symmetric = TRUE); pos <- eig$values > 1e-10
  U <- eig$vectors[, pos, drop = FALSE]; D <- eig$values[pos]
  fam <- NULL
  if (cohort == "BREEDING") {
    fm <- fread(file.path(PROJECT_ROOT, "DATA/METADATA/breeding_sample2family.txt"),
                header = FALSE, col.names = c("sample", "family"))
    fam <- factor(fm$family[match(rownames(inp$M), fm$sample)])
  }
  res <- rbindlist(lapply(idx, function(i) {
    y <- inp$M[, i]
    r <- as.list(emma_reml(y, U, D))
    if (!is.null(fam)) r <- c(r, as.list(pedigree_h2(y, fam)))
    as.data.table(r)
  }), fill = TRUE)
  res[, `:=`(site = colnames(inp$M)[idx], cohort = cohort, context = ctx)]
  fwrite(res, file.path(H2_ROOT, "CHUNKS", sprintf("%s_%s_%04d.tsv.gz", cohort, ctx, k)), sep = "\t")
  msg(cohort, "/", ctx, " chunk ", k, ": ", nrow(res), " sites")
  quit(save = "no")
}

############################################################
# COLLECT MODE
############################################################
h2_all <- rbindlist(lapply(list.files(file.path(H2_ROOT, "CHUNKS"), "\\.tsv\\.gz$", full.names = TRUE), fread), fill = TRUE)
for (co in COHORTS) for (cx in CONTEXTS) {
  n_exp <- ncol(readRDS(file.path(IN_ROOT, co, cx, "methylation_mvalues_matrix.rds")))
  n_got <- h2_all[cohort == co & context == cx, .N]
  if (n_got != n_exp) stop(co, "/", cx, ": ", n_got, " of ", n_exp, " sites in chunks")
}
setorder(h2_all, cohort, context, site)

# robust meQTL sites, cohort-specific
rob <- rbindlist(lapply(COHORTS, function(co)
  unique(fread(file.path(ROB_DIR, sprintf("robust_pairs_%s.tsv", tolower(co))), select = c("site", "context"))[, cohort := co])))
h2_all[, is_meqtl := FALSE]
h2_all[rob, on = .(cohort, context, site), is_meqtl := TRUE]
fwrite(h2_all, file.path(H2_ROOT, "h2_all_sites.tsv.gz"), sep = "\t")

summ <- h2_all[, .(n_sites = .N, n_valid = sum(!is.na(h2)), mean_h2 = mean(h2, na.rm = TRUE),
                   median_h2 = median(h2, na.rm = TRUE), pct_h2_gt0.3 = 100 * mean(h2 > 0.3, na.rm = TRUE),
                   mean_h2_ped = if ("h2_ped" %in% names(.SD)) mean(h2_ped, na.rm = TRUE) else NA_real_),
               by = .(cohort, context)]
setorder(summ, cohort, context)
fwrite(summ, file.path(H2_ROOT, "h2_summary.tsv"), sep = "\t")
mvb <- h2_all[, .(n = .N, mean_h2 = mean(h2, na.rm = TRUE), median_h2 = median(h2, na.rm = TRUE)),
              by = .(cohort, context, group = ifelse(is_meqtl, "robust_meqtl", "background"))]
p_w <- h2_all[, {
  w <- tryCatch(wilcox.test(h2[is_meqtl], h2[!is_meqtl]), error = function(e) NULL)
  .(wilcox_W = if (is.null(w)) NA_real_ else unname(w$statistic),
    wilcox_p = if (is.null(w)) NA_real_ else w$p.value)
}, by = .(cohort, context)]
mvb <- merge(mvb, p_w, by = c("cohort", "context"))
fwrite(mvb, file.path(H2_ROOT, "h2_meqtl_vs_background.tsv"), sep = "\t")
print(summ); print(mvb)

# ---- Mantel test: GRM kinship vs methylation similarity ----
set.seed(42)
mantel_rows <- list()
for (cohort in COHORTS) for (ctx in CONTEXTS) {
  inp <- load_inputs(cohort, ctx)
  M <- inp$M
  obs <- rowMeans(!is.na(M))
  keep <- obs >= MIN_OBS_MANTEL           # leave out failed / very low-complexity libraries
  M <- M[keep, , drop = FALSE]; inp$K <- inp$K[keep, keep]
  d <- as.matrix(dist(M))   # missing values excluded pair by pair (dist scales to all sites)
  sim <- 1 - d / max(d)
  lt <- lower.tri(sim)
  sim_vec <- sim[lt]; K <- inp$K
  r_obs <- cor(K[lt], sim_vec)
  n <- nrow(K)
  perm_r <- vapply(seq_len(N_PERM), function(i) { o <- sample.int(n); cor(K[o, o][lt], sim_vec) }, numeric(1))
  mantel_rows[[paste(cohort, ctx)]] <- data.table(cohort = cohort, context = ctx, n_samples = n,
    n_sites = ncol(M), mantel_r = r_obs, p_value = (sum(perm_r >= r_obs) + 1) / (N_PERM + 1), n_perm = N_PERM)
  msg("Mantel ", cohort, "/", ctx, ": r = ", round(r_obs, 3))
}
mantel <- rbindlist(mantel_rows)
fwrite(mantel, file.path(H2_ROOT, "mantel_grm_methylation.tsv"), sep = "\t")
print(mantel)

# ---- Heritability panel ----
save_fig <- function(base, plot, width, height) {
  ggsave(paste0(base, ".tiff"), plot, width = width, height = height, dpi = 300, compression = "lzw")
  ggsave(paste0(base, ".pdf"),  plot, width = width, height = height)
  ggsave(paste0(base, ".png"),  plot, width = width, height = height, dpi = 300)
  ggsave(paste0(base, ".eps"),  plot, width = width, height = height, device = cairo_ps)
}
# Methylation-context colours, shared by all figures (magenta, blue, turquoise; no red-green pair)
ctx_cols   <- c(CpG = "#D81B60", CHG = "#1F4E9C", CHH = "#17BECF")
coh_labels <- c(BREEDING = "Breeding (n = 209)", NATURAL = "Natural (n = 393)")
hp <- h2_all[!is.na(h2)]
hp[, context_f := factor(context, levels = CONTEXTS)]
hp[, cohort_f  := factor(cohort, levels = COHORTS, labels = coh_labels)]
hp[, group := factor(ifelse(is_meqtl, "Robust meQTL sites", "Other sites"), levels = c("Other sites", "Robust meQTL sites"))]

p_a <- ggplot(hp, aes(x = context_f, y = h2, fill = group, colour = group)) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.5, position = position_dodge(0.8), linewidth = 0.3) +
  geom_boxplot(width = 0.06, outlier.shape = NA, fill = NA, position = position_dodge(0.8), linewidth = 0.4, coef = 0) +
  facet_wrap(~cohort_f, nrow = 1) +
  scale_fill_manual(values = c("Other sites" = "grey60", "Robust meQTL sites" = "#E63946"), name = NULL) +
  scale_colour_manual(values = c("Other sites" = "grey40", "Robust meQTL sites" = "#9B1A22"), name = NULL) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(x = "Methylation context", y = expression(italic(h)[SNP]^2)) +
  theme_bw(base_size = 11) +
  theme(strip.background = element_rect(fill = "grey90"), legend.position = "bottom", panel.grid.minor = element_blank())

bp <- h2_all[cohort == "BREEDING" & !is.na(h2) & !is.na(h2_ped)]
bp[, context_f := factor(context, levels = CONTEXTS)]
set.seed(42)
bp <- bp[, .SD[sample(.N, min(.N, 50000L))], by = context_f]
p_b <- ggplot(bp, aes(x = h2, y = h2_ped, colour = context_f)) +
  geom_point(alpha = 0.08, size = 0.4, stroke = 0) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "black", linewidth = 0.5) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.8) +
  facet_wrap(~context_f, nrow = 1) +
  scale_colour_manual(values = ctx_cols, guide = "none") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(x = expression(italic(h)[SNP]^2), y = expression(italic(h)[Ped]^2)) +
  theme_bw(base_size = 11) +
  theme(strip.background = element_rect(fill = "grey90"), panel.grid.minor = element_blank())

panel_plot <- (p_a + labs(tag = "a") + theme(plot.tag = element_text(face = "bold", size = 14))) /
        (p_b + labs(tag = "b") + theme(plot.tag = element_text(face = "bold", size = 14))) +
        plot_layout(heights = c(4.5, 3.5))
save_fig(file.path(FIG_DIR, "heritability_panel"), panel_plot, width = 8, height = 8)

# Values shown in panels a and b
fwrite(hp[, .(panel = "a", cohort, context, site, robust_cis_meqtl = is_meqtl, h2_snp = h2)],
       file.path(FIG_DIR, "heritability_panel_h2_values.tsv"), sep = "\t")
fwrite(bp[, .(panel = "b", cohort, context, site, h2_snp = h2, h2_ped)],
       file.path(FIG_DIR, "heritability_panel_h2_ped_values.tsv"), sep = "\t")
msg("Step 16ab done")
