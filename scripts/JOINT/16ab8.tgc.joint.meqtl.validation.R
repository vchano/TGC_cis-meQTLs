#!/usr/bin/env Rscript
############################################################
# TreeGeneClimate (TGC) — JOINT ECS + TMS
# Step 16ab8: cis-meQTL validation analyses
#
# Three complementary analyses to validate and contextualise
# the cis-meQTL results:
#
# SECTION 1 — LD deconvolution of hotspot SNP counts
#   For each gene with ≥5 robust cis-meQTL SNPs, counts the
#   number of LD-independent signals using PLINK
#   --indep-pairwise (r² < 0.1, 1 Mb sliding window).
#   Requires PLINK v1.9 in PATH.
#   OUTPUT: RESULTS/JOINT/LD_CLUMP/ld_clump_results*.csv
#
# SECTION 2 — Whole-methylome Mantel test
#   Correlates pairwise GRM kinship with pairwise methylation
#   similarity across all MEF-retained sites (Mantel r,
#   9,999 permutations). Provides a continuous-kinship
#   complement to the categorical family/stand analyses.
#   OUTPUT: RESULTS/JOINT/MANTEL/mantel_results.csv + figures
#
# SECTION 3 — Permutation FDR validation
#   Runs MatrixEQTL on 100 label-swap permutations (CpG
#   context only) to empirically validate the BH-FDR
#   calibration of the real cis-meQTL results.
#   OUTPUT: RESULTS/JOINT/PERMUTATIONS/perm_*.csv + figures
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(MatrixEQTL)
  library(SNPRelate)
  library(gdsfmt)
})

options(stringsAsFactors = FALSE)

# === USER CONFIGURATION ===
PROJECT_ROOT <- "/path/to/your/project"   # <-- set this
RDATA_DIR    <- "/path/to/your/rdata"     # directory with GRM .rds and GDS files
# ===========================

COHORTS  <- c("BREEDING", "NATURAL")
CONTEXTS <- c("CpG", "CHG", "CHH")

msg <- function(...) {
  cat(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n"))
  flush(stdout())
}

############################################################
# SECTION 1 — LD deconvolution of hotspot SNP counts
############################################################

msg("=== Section 1: LD deconvolution ===")

LD_OUTDIR <- file.path(PROJECT_ROOT, "RESULTS/JOINT/LD_CLUMP")
LD_TMPDIR <- file.path(LD_OUTDIR, "tmp_plink")
dir.create(LD_OUTDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LD_TMPDIR, recursive = TRUE, showWarnings = FALSE)

BFILES <- list(
  BREEDING = file.path(PROJECT_ROOT, "ECS/STRUCTURE.ECS/BREEDING/tgc.ecs.breeding.gwas"),
  NATURAL  = file.path(PROJECT_ROOT, "ECS/STRUCTURE.ECS/NATURAL/tgc.ecs.natural.gwas")
)
ANN_FILES <- list(
  BREEDING = file.path(PROJECT_ROOT,
    "RESULTS/JOINT/ANNOTATION17/robust_breeding_snps_annotated.tsv"),
  NATURAL  = file.path(PROJECT_ROOT,
    "RESULTS/JOINT/ANNOTATION17/robust_natural_snps_annotated.tsv")
)
GENE_SUMMARY <- file.path(PROJECT_ROOT,
  "RESULTS/JOINT/ANNOTATION17/gene_summary_annotated.tsv")

MIN_SNPS  <- 5L
R2_THRESH <- 0.1
WIN_KB    <- 1000

ld_results <- list()

for (cohort in COHORTS) {
  msg("  Cohort: ", cohort)
  snp_dt <- fread(ANN_FILES[[cohort]], select = c("gene_id", "snp", "chr", "pos"))
  snp_dt[, plink_id := paste0(chr, ":", pos)]
  gene_snp    <- snp_dt[, .(plink_ids = list(unique(plink_id)),
                             n_snps    = uniqueN(plink_id)), by = gene_id]
  genes_to_run <- gene_snp[n_snps >= MIN_SNPS, gene_id]
  msg("  Genes with >= ", MIN_SNPS, " SNPs: ", length(genes_to_run))

  for (gene in genes_to_run) {
    snps   <- gene_snp[gene_id == gene, plink_ids][[1]]
    prefix <- file.path(LD_TMPDIR, paste0(cohort, "_", gene))
    writeLines(snps, paste0(prefix, ".extract"))
    ret <- system(paste("plink --bfile", BFILES[[cohort]],
                        "--extract", paste0(prefix, ".extract"),
                        "--allow-extra-chr",
                        "--indep-pairwise", WIN_KB, "1", R2_THRESH,
                        "--out", prefix, "--silent"),
                  ignore.stdout = TRUE, ignore.stderr = TRUE)
    prune_in  <- paste0(prefix, ".prune.in")
    n_indep   <- if (file.exists(prune_in)) length(readLines(prune_in)) else NA_integer_
    ld_results[[length(ld_results) + 1]] <- data.table(
      gene_id = gene, cohort = cohort,
      n_snps_raw = length(snps), n_independent = n_indep,
      r2_threshold = R2_THRESH, window_kb = WIN_KB)
    invisible(file.remove(list.files(LD_TMPDIR,
      pattern = paste0("^", cohort, "_", gene, "\\."), full.names = TRUE)))
  }
}

ld_dt <- rbindlist(ld_results)
setorder(ld_dt, cohort, -n_snps_raw)
fwrite(ld_dt, file.path(LD_OUTDIR, "ld_clump_results.csv"))
msg("Saved: ld_clump_results.csv (", nrow(ld_dt), " rows)")

if (file.exists(GENE_SUMMARY)) {
  ann <- fread(GENE_SUMMARY, select = c("gene_id", "gene_chr", "gene_start",
                                         "gene_end", "best_annotation_class",
                                         "eggnog_description"))
  ld_ann <- merge(ld_dt, ann, by = "gene_id", all.x = TRUE)
  setorder(ld_ann, cohort, -n_snps_raw)
  fwrite(ld_ann, file.path(LD_OUTDIR, "ld_clump_results_annotated.csv"))
  msg("Saved: ld_clump_results_annotated.csv")
}
unlink(LD_TMPDIR, recursive = TRUE)

############################################################
# SECTION 2 — Whole-methylome Mantel test
############################################################

msg("=== Section 2: Whole-methylome Mantel test ===")

MANTEL_OUTDIR <- file.path(PROJECT_ROOT, "RESULTS/JOINT/MANTEL")
dir.create(MANTEL_OUTDIR, recursive = TRUE, showWarnings = FALSE)

set.seed(42)
N_PERM <- 9999

mantel_test <- function(kin_vec, sim_vec, n_perm = 9999) {
  r_obs <- cor(kin_vec, sim_vec, method = "pearson")
  if (is.na(r_obs)) return(list(r = NA_real_, p = NA_real_))
  n     <- round((1 + sqrt(1 + 8 * length(kin_vec))) / 2)
  perm_r <- numeric(n_perm)
  for (i in seq_len(n_perm)) {
    idx <- sample(n)
    k_perm <- matrix(0, n, n)
    k_perm[lower.tri(k_perm)] <- kin_vec
    k_perm <- k_perm + t(k_perm)
    k_perm <- k_perm[idx, idx]
    perm_r[i] <- cor(k_perm[lower.tri(k_perm)], sim_vec, method = "pearson")
  }
  list(r = r_obs, p = (sum(perm_r >= r_obs) + 1) / (n_perm + 1))
}

mantel_results <- list()

for (cohort in COHORTS) {
  msg("  Loading GRM: ", cohort)
  K <- as.matrix(readRDS(
    file.path(RDATA_DIR, sprintf("%s_grm_gcta.rds", tolower(cohort)))))

  for (ctx in CONTEXTS) {
    mval_path <- file.path(PROJECT_ROOT, "RESULTS/JOINT/MQTL5/INPUTS",
                           cohort, ctx, "methylation_mvalues_matrix.rds")
    if (!file.exists(mval_path)) { msg("  MISSING: ", mval_path); next }
    Mmat   <- as.matrix(readRDS(mval_path))
    common <- intersect(rownames(K), rownames(Mmat))
    K_sub  <- K[common, common, drop = FALSE]
    M_sub  <- Mmat[common, , drop = FALSE]
    n      <- nrow(K_sub)
    msg("  ", cohort, "/", ctx, ": ", n, " samples | ", ncol(M_sub), " sites")

    meth_dist <- as.matrix(dist(M_sub, method = "euclidean"))
    max_d     <- max(meth_dist)
    meth_sim  <- if (max_d > 0) 1 - meth_dist / max_d else matrix(1, n, n)
    lt        <- lower.tri(K_sub)
    res       <- mantel_test(K_sub[lt], meth_sim[lt], n_perm = N_PERM)
    msg("  r = ", round(res$r, 4), " | p = ", signif(res$p, 3))

    mantel_results[[length(mantel_results) + 1]] <- data.table(
      cohort = cohort, context = ctx, n_samples = n,
      n_pairs = sum(lt), mantel_r = res$r, mantel_p = res$p, n_perm = N_PERM)
  }
}

mantel_dt <- rbindlist(mantel_results)
setorder(mantel_dt, cohort, context)
fwrite(mantel_dt, file.path(MANTEL_OUTDIR, "mantel_results.csv"))
msg("Saved: mantel_results.csv")
print(mantel_dt[, .(cohort, context, n_samples, mantel_r, mantel_p)])

############################################################
# SECTION 3 — Permutation FDR validation (CpG, 100 permutations)
############################################################

msg("=== Section 3: Permutation FDR validation ===")

PERM_OUTDIR <- file.path(PROJECT_ROOT, "RESULTS/JOINT/PERMUTATIONS")
dir.create(PERM_OUTDIR, recursive = TRUE, showWarnings = FALSE)

N_PERM_FDR  <- 100L
CIS_WINDOW  <- 1e5

GDS_FILES <- list(
  BREEDING = file.path(RDATA_DIR, "breeding.snp.gds"),
  NATURAL  = file.path(RDATA_DIR, "natural.snp.gds")
)
OBS_FILES <- list(
  BREEDING = file.path(PROJECT_ROOT,
    "RESULTS/JOINT/MATRIXEQTL5/BREEDING/CpG/cis_meqtl_significant.tsv"),
  NATURAL  = file.path(PROJECT_ROOT,
    "RESULTS/JOINT/MATRIXEQTL5/NATURAL/CpG/cis_meqtl_significant.tsv")
)

normalize_id <- function(x) {
  x <- as.character(x); x <- trimws(x)
  x <- gsub("\\.0$", "", x); x <- gsub("^X", "", x); x <- gsub("-", "_", x); x
}

impute_hwe <- function(geno_snp_x_samp, groups) {
  mode(geno_snp_x_samp) <- "numeric"
  geno_snp_x_samp[geno_snp_x_samp > 2 | geno_snp_x_samp < 0] <- NA_real_
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
  if (anyNA(g)) {
    fill2 <- matrix(2 * p_all, nrow = nrow(g), ncol = ncol(g))
    g[is.na(g)] <- fill2[is.na(g)]
  }
  list(geno = g, keep = keep)
}

perm_results <- list()

for (cohort in COHORTS) {
  ctx       <- "CpG"
  panel_dir <- file.path(PROJECT_ROOT, "RESULTS/JOINT/MQTL5/INPUTS", cohort, ctx)
  perm_dir  <- file.path(PERM_OUTDIR, cohort, ctx)
  dir.create(perm_dir, recursive = TRUE, showWarnings = FALSE)

  msg("  Loading panel inputs: ", cohort, "/", ctx)
  m_mat      <- readRDS(file.path(panel_dir, "methylation_mvalues_matrix.rds"))
  site_annot <- as.data.table(readRDS(file.path(panel_dir, "methylation_site_annot.rds")))
  shared_dt  <- readRDS(file.path(panel_dir, "shared_sample_ids.rds"))
  pcs_dt     <- as.data.table(readRDS(file.path(panel_dir, "pcs_shared.rds")))
  groups_dt  <- as.data.table(readRDS(file.path(panel_dir, "sample_groups.rds")))
  var_annot  <- as.data.table(readRDS(file.path(panel_dir, "snp_variant_annot.rds")))
  shared_ids <- shared_dt$sample_id
  n_shared   <- length(shared_ids)

  msg("  Loading genotypes from GDS...")
  gds          <- snpgdsOpen(GDS_FILES[[cohort]], readonly = TRUE)
  gds_samp_ids <- read.gdsn(index.gdsn(gds, "sample.id"))
  gds_norm     <- normalize_id(gds_samp_ids)
  raw_shared   <- gds_samp_ids[gds_norm %in% shared_ids]
  geno_list    <- snpgdsGetGeno(gds, sample.id = raw_shared,
                                 snpfirstdim = TRUE, with.id = TRUE)
  snpgdsClose(gds)
  samp_reorder <- match(shared_ids, normalize_id(geno_list$sample.id))
  geno_mat     <- geno_list$genotype[, samp_reorder, drop = FALSE]
  rm(geno_list, samp_reorder); gc()

  groups_ord <- groups_dt$group[match(shared_ids, groups_dt$sample_id)]
  groups_ord[is.na(groups_ord)] <- "Unknown"
  imp      <- impute_hwe(geno_mat, groups_ord)
  geno_imp <- imp$geno
  var_filt <- var_annot[imp$keep]
  snp_sd   <- apply(geno_imp, 1, sd, na.rm = TRUE)
  keep_var <- is.finite(snp_sd) & snp_sd > 0
  geno_imp <- geno_imp[keep_var, , drop = FALSE]
  var_filt <- var_filt[keep_var]
  rownames(geno_imp) <- as.character(var_filt$snp_id)
  colnames(geno_imp) <- shared_ids
  rm(geno_mat, imp); gc()

  pc_cols  <- grep("^PC[0-9]+$", names(pcs_dt), value = TRUE)
  pcs_ord  <- pcs_dt[match(shared_ids, pcs_dt$sample_id), ..pc_cols]
  snpspos  <- data.frame(snp = as.character(var_filt$snp_id),
                         chr = as.character(var_filt$chr),
                         pos = as.integer(var_filt$pos),
                         stringsAsFactors = FALSE)
  genepos  <- data.frame(gene  = colnames(m_mat),
                         chr   = as.character(site_annot$chr),
                         start = as.integer(site_annot$start),
                         end   = as.integer(site_annot$end),
                         stringsAsFactors = FALSE)
  m_reorder <- match(shared_ids, normalize_id(rownames(m_mat)))
  m_ord     <- m_mat[m_reorder, , drop = FALSE]

  msg("  Running ", N_PERM_FDR, " permutations...")
  for (perm_id in seq_len(N_PERM_FDR)) {
    out_file <- file.path(perm_dir, sprintf("perm_%03d.tsv", perm_id))
    if (file.exists(out_file)) next

    set.seed(perm_id)
    m_perm <- m_ord[sample(n_shared), , drop = FALSE]
    rownames(m_perm) <- shared_ids

    snpsSD <- SlicedData$new(); snpsSD$CreateFromMatrix(geno_imp)
    geneSD <- SlicedData$new(); geneSD$CreateFromMatrix(t(m_perm))
    cvrtSD <- SlicedData$new(); cvrtSD$CreateFromMatrix(t(as.matrix(pcs_ord)))

    tmp_out <- tempfile(fileext = ".txt")
    me <- tryCatch(
      Matrix_eQTL_main(
        snps = snpsSD, gene = geneSD, cvrt = cvrtSD,
        output_file_name = "", pvOutputThreshold = 0,
        useModel = modelLINEAR, errorCovariance = numeric(),
        verbose = FALSE,
        output_file_name.cis = tmp_out, pvOutputThreshold.cis = 1,
        snpspos = snpspos, genepos = genepos, cisDist = CIS_WINDOW),
      error = function(e) NULL)
    if (file.exists(tmp_out)) file.remove(tmp_out)

    if (is.null(me)) next
    cis_dt    <- as.data.table(me$cis$eqtls)
    cis_dt[, p_FDR := p.adjust(pvalue, method = "BH")]
    n_sig     <- sum(cis_dt$p_FDR < 0.05, na.rm = TRUE)
    lambda    <- median(qchisq(1 - cis_dt$pvalue[is.finite(cis_dt$pvalue) &
                               cis_dt$pvalue > 0], df = 1)) / qchisq(0.5, df = 1)
    fwrite(data.table(cohort = cohort, context = ctx, perm_id = perm_id,
                      n_pairs = nrow(cis_dt), n_sig_BH05 = n_sig,
                      lambda = lambda, status = "ok"),
           out_file, sep = "\t")
    if (perm_id %% 10 == 0) msg("  Perm ", perm_id, "/", N_PERM_FDR,
                                 " — n_sig=", n_sig)
  }

  perm_files <- list.files(perm_dir, pattern = "^perm_\\d+\\.tsv$",
                            full.names = TRUE)
  perm_dt    <- rbindlist(lapply(perm_files, fread), use.names = TRUE, fill = TRUE)
  perm_dt    <- perm_dt[status == "ok"]
  perm_results[[cohort]] <- perm_dt
  msg("  Loaded ", nrow(perm_dt), " completed permutations for ", cohort)
}

# Aggregate and compute empirical FDR
all_perm <- rbindlist(perm_results, use.names = TRUE, fill = TRUE)
fwrite(all_perm, file.path(PERM_OUTDIR, "perm_summary.csv"))

obs_dt <- rbindlist(lapply(COHORTS, function(coh) {
  f <- OBS_FILES[[coh]]
  n <- if (file.exists(f)) nrow(fread(f, select = "pvalue")) else NA_integer_
  data.table(cohort = coh, context = "CpG", n_sig_obs = n)
}))

efdr_dt <- all_perm[, .(
  n_perm          = .N,
  mean_n_sig_perm = mean(n_sig_BH05, na.rm = TRUE),
  sd_n_sig_perm   = sd(n_sig_BH05, na.rm = TRUE)
), by = .(cohort, context)]
efdr_dt <- merge(efdr_dt, obs_dt, by = c("cohort", "context"), all.x = TRUE)
efdr_dt[, empirical_FDR := mean_n_sig_perm / n_sig_obs]
fwrite(efdr_dt, file.path(PERM_OUTDIR, "perm_empirical_fdr.csv"))
msg("Empirical FDR:")
print(efdr_dt[, .(cohort, context, n_sig_obs, mean_n_sig_perm, empirical_FDR)])

# Figure: null distribution
coh_labels <- c(BREEDING = "Breeding (n=209)", NATURAL = "Natural (n=393)")
plot_dt  <- all_perm[context == "CpG"]
obs_line <- obs_dt[context == "CpG"]
plot_dt[,  cohort_f := factor(cohort, levels = c("BREEDING","NATURAL"), labels = coh_labels)]
obs_line[, cohort_f := factor(cohort, levels = c("BREEDING","NATURAL"), labels = coh_labels)]

p_perm <- ggplot(plot_dt, aes(x = n_sig_BH05)) +
  geom_histogram(binwidth = 1, fill = "grey70", colour = "grey50", linewidth = 0.3) +
  geom_vline(data = obs_line, aes(xintercept = n_sig_obs),
             colour = "#E63946", linewidth = 1, linetype = "dashed") +
  geom_text(data = obs_line,
            aes(x = n_sig_obs, y = Inf,
                label = paste0("Observed:\n", format(n_sig_obs, big.mark = ","))),
            hjust = -0.1, vjust = 1.3, size = 3, colour = "#E63946") +
  facet_wrap(~cohort_f, scales = "free_x") +
  labs(x = "Significant pairs (BH FDR < 0.05) under permuted null",
       y = "Count (permutations)",
       subtitle = paste0(N_PERM_FDR, " label-swap permutations; dashed line = observed")) +
  theme_bw(base_size = 11) +
  theme(strip.background = element_rect(fill = "grey90"),
        panel.grid.minor  = element_blank())

ggsave(file.path(PERM_OUTDIR, "perm_null_distribution.tiff"), p_perm,
       width = 8, height = 4, dpi = 300, compression = "lzw")
ggsave(file.path(PERM_OUTDIR, "perm_null_distribution.pdf"),  p_perm, width = 8, height = 4)
ggsave(file.path(PERM_OUTDIR, "perm_null_distribution.png"),  p_perm, width = 8, height = 4, dpi = 300)

msg("Step 16ab8 complete.")
