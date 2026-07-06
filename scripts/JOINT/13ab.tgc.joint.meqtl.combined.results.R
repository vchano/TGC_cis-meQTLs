#!/usr/bin/env Rscript
############################################################
# TreeGeneClimate (TGC) — JOINT ECS + TBS
# Step 13ab: Combined results — comprehensive summary, QQ plots, sig-site tables
#
# PURPOSE
# Reads cis-meQTL results from GENESIS5 and MatrixEQTL5, produces:
#   - Comprehensive summary (lambda, n_sites, n_sig, cis window — all panels)
#   - QQ plots (one per panel x tool, 12 total, TIFF)
#   - Long-format site FDR table (site, site_chr, context, min_FDR)
#   - Significant site lists (p_FDR < 5e-8 and p_FDR < 1e-10)
#
# SIGNIFICANCE THRESHOLDS (BH-adjusted p-values):
#   FDR_LOOSE  = 5e-8
#   FDR_STRICT = 1e-10
#
# USAGE
#   Rscript --vanilla 13ab.R
#
# INPUTS
#   RESULTS/JOINT/GENESIS5/MEQTL/<COHORT>/<CONTEXT>/cis_meqtl_all_results.rds
#   RESULTS/JOINT/MATRIXEQTL5/<COHORT>/<CONTEXT>/cis_meqtl_all_results.rds
#   RESULTS/JOINT/MQTL5/INPUTS/<COHORT>/<CONTEXT>/snp_variant_annot.rds
#
# OUTPUTS (RESULTS/JOINT/COMBINED5/)
#   comprehensive_summary.tsv
#   qq/qqplot_<tool>_<cohort>_<context>.tiff
#   tables/all_sites_<tool>_<cohort>.tsv   — long format: site, site_chr, context, min_FDR
#   sig_sites/sig_p5e8_<tool>_<cohort>.tsv
#   sig_sites/sig_p1e10_<tool>_<cohort>.tsv
############################################################

suppressPackageStartupMessages({
  library(data.table)
})

options(stringsAsFactors = FALSE)

############################################################
# 1) SETTINGS
############################################################

COHORTS       <- c("BREEDING", "NATURAL")
CONTEXTS      <- c("CpG", "CHG", "CHH")
TOOLS         <- c("GENESIS5", "MATRIXEQTL5")
FDR_LOOSE     <- 5e-8   # genome-wide FDR threshold (loose; used for sig-site tables)
FDR_STRICT    <- 1e-10  # strict threshold used in all downstream analyses
CIS_WINDOW_KB <- 100L   # cis window radius applied during meQTL mapping

############################################################
# 2) PATHS
############################################################

# === USER CONFIGURATION ===
PROJECT_ROOT <- "/path/to/your/project"  # <-- set this
# ===========================

# GENESIS5 results are stored under MEQTL/; MatrixEQTL5 directly under MATRIXEQTL5/
RESULT_ROOTS <- list(
  GENESIS5    = file.path(PROJECT_ROOT, "RESULTS", "JOINT", "GENESIS5",    "MEQTL"),
  MATRIXEQTL5 = file.path(PROJECT_ROOT, "RESULTS", "JOINT", "MATRIXEQTL5")
)

# Variant annotation RDS files produced by the MQTL5 input-prep step (12ab0.R)
ANNOT_ROOT <- file.path(PROJECT_ROOT, "RESULTS", "JOINT", "MQTL5", "INPUTS")

OUT_ROOT <- file.path(PROJECT_ROOT, "RESULTS", "JOINT", "COMBINED5")
QQ_DIR   <- file.path(OUT_ROOT, "qq")
TAB_DIR  <- file.path(OUT_ROOT, "tables")
SIG_DIR  <- file.path(OUT_ROOT, "sig_sites")
for (d in c(OUT_ROOT, QQ_DIR, TAB_DIR, SIG_DIR))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

############################################################
# 3) HELPERS
############################################################

msg <- function(...) cat(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n"))

# Genomic inflation factor lambda — ratio of observed to expected chi-squared median.
# Values > 1 indicate p-value inflation relative to null; used as QC metric.
compute_lambda <- function(pvals) {
  pvals <- pvals[is.finite(pvals) & pvals > 0 & pvals <= 1]
  if (length(pvals) < 2) return(NA_real_)
  median(qchisq(1 - pvals, df = 1), na.rm = TRUE) / qchisq(0.5, df = 1)
}

# Standardise columns from either tool's RDS output.
# Returns a data.table with columns:
#   snp, snp_chr, snp_pos, site, site_chr, site_pos,
#   context, statistic, beta, pvalue, p_FDR, distance
clean_result <- function(dt, tool, context, va) {
  dt <- as.data.table(dt)

  if (tool == "GENESIS5") {
    # snp_int is the raw GENESIS variant.id; snp (added by 12ab2.R) is the
    # canonical integer matching snp_variant_annot$snp_id — keep snp, drop snp_int
    drop_cols <- intersect(names(dt),
                           c("snp_int", "chr", "pos", "n.obs", "freq", "MAC",
                             "Score", "Score.SE", "Est.SE", "PVE"))
    if (length(drop_cols)) dt[, (drop_cols) := NULL]
    # Est -> beta
    if ("Est" %in% names(dt) && !"beta" %in% names(dt))
      setnames(dt, "Est", "beta")
    # Ensure p_FDR is present; compute from raw p-values if missing
    if (!"p_FDR" %in% names(dt) && "FDR" %in% names(dt))
      setnames(dt, "FDR", "p_FDR")
    if (!"p_FDR" %in% names(dt))
      dt[, p_FDR := p.adjust(pvalue, method = "BH")]
  }

  if (tool == "MATRIXEQTL5") {
    # FDR (MatrixEQTL native) and p_FDR (added by 12ab1.R) are identical — keep p_FDR
    if ("FDR" %in% names(dt) && "p_FDR" %in% names(dt))
      dt[, FDR := NULL]
    else if ("FDR" %in% names(dt) && !"p_FDR" %in% names(dt))
      setnames(dt, "FDR", "p_FDR")
    if (!"p_FDR" %in% names(dt))
      dt[, p_FDR := p.adjust(pvalue, method = "BH")]
    if (!"beta" %in% names(dt))
      dt[, beta := NA_real_]
  }

  # site_chr: parse from site name if not already present
  # e.g. "PA_chr03:790910176-790910176" -> "PA_chr03"
  if (!"site_chr" %in% names(dt))
    dt[, site_chr := sub(":.*", "", site)]

  # snp_chr: join from variant annotation on snp = snp_id
  if (!is.null(va) && "snp" %in% names(dt)) {
    if ("chr" %in% names(dt)) dt[, chr := NULL]   # drop any stale chr column
    va_sub <- va[, .(snp_id, snp_chr = chr)]
    dt <- merge(dt, va_sub, by.x = "snp", by.y = "snp_id", all.x = TRUE)
  } else {
    dt[, snp_chr := NA_character_]
  }

  # Add context column
  dt[, context := context]

  # Ensure beta exists
  if (!"beta" %in% names(dt)) dt[, beta := NA_real_]

  # Select and reorder final columns
  keep <- c("snp", "snp_chr", "snp_pos", "site", "site_chr", "site_pos",
            "context", "statistic", "beta", "pvalue", "p_FDR", "distance")
  keep <- intersect(keep, names(dt))
  dt[, .SD, .SDcols = keep]
}

############################################################
# 4) COLLECT RESULTS
############################################################

summary_rows    <- list()
sig_loose_list  <- list()   # p_FDR < 5e-8
sig_strict_list <- list()   # p_FDR < 1e-10

# Outer loops: tool → cohort → context (6 panels per tool)
for (tool in TOOLS) {
  rroot <- RESULT_ROOTS[[tool]]

  for (cohort in COHORTS) {
    long_list <- list()   # per-context long-format entries for all_sites table

    for (context in CONTEXTS) {
      panel_tag <- paste0(tolower(cohort), "_", tolower(context))
      rds_path  <- file.path(rroot, cohort, context, "cis_meqtl_all_results.rds")

      # Record missing panels in the summary so the output table is always complete
      if (!file.exists(rds_path)) {
        msg("Missing: ", rds_path, " — skipping")
        summary_rows[[paste(tool, panel_tag)]] <- data.table(
          tool = tool, cohort = cohort, context = context,
          cis_window_kb = CIS_WINDOW_KB,
          n_pairs = NA_integer_, n_sites_tested = NA_integer_,
          n_snps_tested = NA_integer_,
          n_sig_pairs_p1e10 = NA_integer_, n_sig_sites_p1e10 = NA_integer_,
          n_sig_snps_p1e10 = NA_integer_,
          lambda = NA_real_, status = "missing")
        next
      }

      msg("Reading: ", tool, " ", cohort, " ", context)
      raw_dt <- readRDS(rds_path)

      # Load SNP variant annotation (snp_id, chr, pos)
      va_path <- file.path(ANNOT_ROOT, cohort, context, "snp_variant_annot.rds")
      va <- if (file.exists(va_path)) as.data.table(readRDS(va_path)) else NULL

      # Harmonise column names and types across both tools
      dt <- clean_result(raw_dt, tool, context, va)
      rm(raw_dt); gc()

      # Per-panel summary statistics
      lambda         <- compute_lambda(dt$pvalue)
      n_sites        <- dt[, uniqueN(site)]
      n_snps         <- dt[, uniqueN(snp)]
      n_sig_l        <- dt[p_FDR < FDR_LOOSE,  uniqueN(site)]
      n_sig_s_sites  <- dt[p_FDR < FDR_STRICT, uniqueN(site)]
      n_sig_s_snps   <- dt[p_FDR < FDR_STRICT, uniqueN(snp)]
      n_sig_s_pairs  <- dt[p_FDR < FDR_STRICT, .N]

      msg("  ", panel_tag,
          " | pairs=", nrow(dt),
          " | sites=", n_sites,
          " | snps=", n_snps,
          " | lambda=", round(lambda, 4),
          " | sig_1e-10: pairs=", n_sig_s_pairs,
          "  sites=", n_sig_s_sites,
          "  snps=", n_sig_s_snps)

      # One row per panel in the comprehensive summary
      summary_rows[[paste(tool, panel_tag)]] <- data.table(
        tool              = tool,
        cohort            = cohort,
        context           = context,
        cis_window_kb     = CIS_WINDOW_KB,
        n_pairs           = nrow(dt),
        n_sites_tested    = n_sites,
        n_snps_tested     = n_snps,
        n_sig_pairs_p1e10 = n_sig_s_pairs,
        n_sig_sites_p1e10 = n_sig_s_sites,
        n_sig_snps_p1e10  = n_sig_s_snps,
        lambda            = round(lambda, 4),
        status            = "ok"
      )

      # QQ plots: TIFF only for render performance; EPS/PNG skipped (see note below).
      # Each plot uses all raw p-values (not just significant ones) to show the
      # full null distribution and quantify inflation via lambda.
      pv_qq <- dt$pvalue
      pv_qq <- pv_qq[is.finite(pv_qq) & pv_qq > 0]
      if (length(pv_qq) >= 2) {
        qq_file <- file.path(QQ_DIR,
          sprintf("qqplot_%s_%s_%s.tiff",
                  tolower(tool), tolower(cohort), tolower(context)))
        n_qq    <- length(pv_qq)
        # Colour distinguishes tools: red = GENESIS5 (LMM), blue = MatrixEQTL5 (LM)
        col_pt  <- if (grepl("GENESIS", tool)) "#e34a33" else "#2b8cbe"
        panel_labels <- c(
          BREEDING_CpG = "a)", BREEDING_CHG = "b)", BREEDING_CHH = "c)",
          NATURAL_CpG  = "d)", NATURAL_CHG  = "e)", NATURAL_CHH  = "f)"
        )
        qq_label <- panel_labels[paste0(cohort, "_", context)]
        draw_qq <- function() {
          # ppoints() generates uniform quantiles for expected distribution;
          # sort(pv_qq) orders observed values for the quantile-quantile comparison
          plot(-log10(ppoints(n_qq)), -log10(sort(pv_qq)),
               pch = 20, cex = 0.4, col = col_pt,
               xlab = "Expected -log10(p)", ylab = "Observed -log10(p)",
               main = "",
               cex.lab = 2.0, cex.axis = 2.0)
          mtext(qq_label, side = 3, adj = 0, line = 0.5, cex = 3.0, font = 1)
          abline(0, 1, col = "red", lty = 2)  # null expectation: observed = expected
        }
        # 2400×2400 px at 300 dpi gives an 8×8 inch print-ready panel
        tiff(qq_file, width = 2400, height = 2400, res = 300, compression = "lzw")
        draw_qq(); dev.off()
        pdf(sub("\\.tiff$", ".pdf", qq_file), width = 8, height = 8)
        draw_qq(); dev.off()
        # EPS and PNG skipped: EPS OOMs on large natural-context datasets (>1 GB vector
        # file); PNG render of millions of scatter points is prohibitively slow.
        # 16ab.R only needs the TIFF; 19ab exports PDF from native source.
      }

      # Significant-site tables (all SNP-site pairs passing threshold)
      sig_l <- dt[p_FDR < FDR_LOOSE][order(p_FDR)]
      sig_s <- dt[p_FDR < FDR_STRICT][order(p_FDR)]
      sig_loose_list[[paste(tool, panel_tag)]]  <- sig_l
      sig_strict_list[[paste(tool, panel_tag)]] <- sig_s

      # Long-format all-sites entry: minimum FDR per methylation site
      # (collapsing all SNP associations per site to its best p-value)
      min_fdr_dt <- dt[, .(
        site_chr = site_chr[1L],
        min_FDR  = min(p_FDR, na.rm = TRUE)
      ), by = site]
      min_fdr_dt[, context := context]
      long_list[[context]] <- min_fdr_dt

      rm(dt); gc()
    }  # context loop

    # Write long-format all-sites table once all three contexts are processed
    if (length(long_list) > 0) {
      long_dt <- rbindlist(long_list, fill = TRUE)
      setcolorder(long_dt, c("site", "site_chr", "context", "min_FDR"))
      setorder(long_dt, context, min_FDR)
      fwrite(long_dt,
             file.path(TAB_DIR,
               sprintf("all_sites_%s_%s.tsv", tolower(tool), tolower(cohort))),
             sep = "\t")
    }

    # Write significant-pair tables for each FDR threshold, combining all contexts
    cohort_key <- paste0(tool, " ", tolower(cohort))
    for (thresh_tag in c("p5e8", "p1e10")) {
      src  <- if (thresh_tag == "p5e8") sig_loose_list else sig_strict_list
      keys <- names(src)[grepl(cohort_key, names(src), ignore.case = TRUE)]
      if (!length(keys)) next
      combined <- rbindlist(src[keys], fill = TRUE)
      if (nrow(combined) > 0)
        fwrite(combined,
               file.path(SIG_DIR,
                 sprintf("sig_%s_%s_%s.tsv",
                         thresh_tag, tolower(tool), tolower(cohort))),
               sep = "\t")
    }
  }  # cohort loop
}  # tool loop

############################################################
# 5) WRITE COMPREHENSIVE SUMMARY
############################################################

# Bind all 12 panel summary rows (2 tools × 2 cohorts × 3 contexts)
summary_dt <- rbindlist(summary_rows, fill = TRUE)

fwrite(summary_dt, file.path(OUT_ROOT, "comprehensive_summary.tsv"), sep = "\t")

msg("Comprehensive summary saved: ", file.path(OUT_ROOT, "comprehensive_summary.tsv"))
cat("\nComprehensive summary (ok panels):\n")
print(summary_dt[status == "ok",
                 .(tool, cohort, context, cis_window_kb,
                   n_pairs, n_sites_tested, n_snps_tested,
                   n_sig_pairs_p1e10, n_sig_sites_p1e10, n_sig_snps_p1e10, lambda)])

msg("Step 13ab finished. Outputs: ", OUT_ROOT)

sessionInfo()
