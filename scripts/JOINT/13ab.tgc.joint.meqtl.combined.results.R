#!/usr/bin/env Rscript
############################################################
# TreeGeneClimate (TGC) — JOINT ECS + TMS
# Step 13ab: Combined cis-meQTL results
#
#   - Panel summary for both tools: tested pairs/sites/SNPs,
#     significant pairs, unique sites and unique SNPs at FDR < 1e-10, lambda
#   - Significant pairs per tool and cohort (FDR < 1e-10, all contexts)
#   - Robust pairs: the same SNP-site pair at FDR < 1e-10 in BOTH GENESIS and
#     MatrixEQTL (pair-level replication), with genomic coordinates
#   - Robust counts per cohort and context
#   - QQ panel: a, b MatrixEQTL and c, d GENESIS (breeding, natural), with the three
#     methylation contexts in each plot; all pairs with p < 1e-3 and a random sample
#     of 20,000 of the remaining pairs per tool, cohort and context
#
# INPUTS
#   RESULTS/JOINT/MEQTL/{MATRIXEQTL,GENESIS}/<COHORT>/<CTX>/cis_meqtl_all_results.rds
#   RESULTS/JOINT/MEQTL/{MATRIXEQTL,GENESIS}/<COHORT>/<CTX>/cis_meqtl_panel_summary.tsv
#
# OUTPUTS (RESULTS/JOINT/MEQTL/COMBINED/)
#   meqtl_panel_summary.tsv
#   sig_fdr1e10_<tool>_<cohort>.tsv
#   robust_pairs_<cohort>.tsv
#   robust_context_summary.tsv
#   qq_panel.{tiff,pdf}
#   qq_panel_values.tsv  (values shown in the QQ panel)
#
# R packages: data.table
############################################################

suppressPackageStartupMessages(library(data.table))
options(stringsAsFactors = FALSE)

# === USER CONFIGURATION ===
PROJECT_ROOT <- Sys.getenv("TGC_PROJECT_ROOT", unset = "/path/to/your/project")
FDR_STRICT <- 1e-10
# ===========================

MEQTL <- file.path(PROJECT_ROOT, "RESULTS/JOINT/MEQTL")
OUT   <- file.path(MEQTL, "COMBINED")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
COHORTS  <- c("BREEDING", "NATURAL")
CONTEXTS <- c("CpG", "CHG", "CHH")
TOOLS    <- c(MatrixEQTL = "MATRIXEQTL", GENESIS = "GENESIS")
msg <- function(...) cat("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n", sep = "")
set.seed(1)   # thinning of points in the QQ plots

panel_rows <- list(); sig <- list(); qq_data <- list()
for (tool in names(TOOLS)) for (cohort in COHORTS) for (ctx in CONTEXTS) {
  d <- file.path(MEQTL, TOOLS[[tool]], cohort, ctx)
  res <- readRDS(file.path(d, "cis_meqtl_all_results.rds"))
  ps  <- fread(file.path(d, "cis_meqtl_panel_summary.tsv"))
  s <- res[p_FDR < FDR_STRICT]
  s[, `:=`(tool = tool, cohort = cohort, context = ctx)]
  sig[[paste(tool, cohort, ctx)]] <- s
  panel_rows[[paste(tool, cohort, ctx)]] <- data.table(
    tool = tool, cohort = cohort, context = ctx, n_samples = ps$n_samples,
    sites_tested = ps$n_sites_with_cis, snps_tested = ps$n_snps_with_cis, pairs_tested = nrow(res),
    sig_pairs = nrow(s), sig_sites = uniqueN(s$site), sig_snps = uniqueN(s$snp), lambda = ps$lambda)
  # QQ data: all points with p < 1e-3, a 20,000-point sample of the rest
  p <- res$pvalue[res$pvalue > 0]
  n <- length(p); p <- sort(p); e <- ppoints(n)
  keep <- which(p < 1e-3)
  rest <- setdiff(seq_len(n), keep)
  if (length(rest) > 2e4) rest <- sort(sample(rest, 2e4))
  idx <- sort(c(keep, rest))
  qq_data[[paste(tool, cohort, ctx)]] <- list(o = -log10(p[idx]), e = -log10(e[idx]))
  msg(tool, " ", cohort, " ", ctx, ": ", nrow(res), " pairs; ", nrow(s), " at FDR < ", FDR_STRICT)
  rm(res); gc()
}
panel <- rbindlist(panel_rows)
fwrite(panel, file.path(OUT, "meqtl_panel_summary.tsv"), sep = "\t")
sig <- rbindlist(sig, fill = TRUE)
for (tl in names(TOOLS)) for (co in COHORTS)
  fwrite(sig[tool == tl & cohort == co], file.path(OUT, sprintf("sig_fdr1e10_%s_%s.tsv", tolower(tl), tolower(co))), sep = "\t")

# ---- robust pairs ----
keep_cols <- c("snp", "site", "context", "snp_chr", "snp_pos", "site_pos", "distance", "beta", "statistic", "pvalue", "p_FDR")
rob_summary <- list()
for (co in COHORTS) {
  m <- sig[tool == "MatrixEQTL" & cohort == co, ..keep_cols]
  g <- sig[tool == "GENESIS"    & cohort == co, ..keep_cols]
  rob <- merge(m, g, by = c("snp", "site", "context", "snp_chr", "snp_pos", "site_pos", "distance"),
               suffixes = c("_MatrixEQTL", "_GENESIS"))
  rob[, site_chr := sub(":.*$", "", site)]
  setcolorder(rob, c("context", "snp", "snp_chr", "snp_pos", "site", "site_chr", "site_pos", "distance"))
  setorder(rob, context, snp_chr, snp_pos)
  rob[, cohort := co]
  fwrite(rob, file.path(OUT, sprintf("robust_pairs_%s.tsv", tolower(co))), sep = "\t")
  rob_summary[[co]] <- rbind(
    rob[, .(pairs = .N, sites = uniqueN(site), snps = uniqueN(snp)), by = context],
    rob[, .(context = "All", pairs = .N, sites = uniqueN(site), snps = uniqueN(snp))])[, cohort := co]
  msg(co, ": ", nrow(rob), " robust pairs, ", uniqueN(rob$site), " sites, ", uniqueN(rob$snp), " SNPs")
}
rob_summary <- rbindlist(rob_summary)
setcolorder(rob_summary, c("cohort", "context"))
fwrite(rob_summary, file.path(OUT, "robust_context_summary.tsv"), sep = "\t")
print(panel); print(rob_summary)

# ---- QQ panel: a, b MatrixEQTL and c, d GENESIS (breeding, natural), the three contexts in each plot ----
CTX_COLORS <- c(CpG = "#D81B60", CHG = "#1F4E9C", CHH = "#17BECF")
draw_qq <- function() {
  layout(matrix(c(1, 2, 3, 4, 5, 5), nrow = 3, byrow = TRUE), heights = c(1, 1, 0.10))
  par(cex = 0.83, mar = c(4.2, 4.5, 2, 1))
  i <- 0
  for (tool in names(TOOLS)) for (cohort in COHORTS) {
    i <- i + 1
    q <- lapply(CONTEXTS, function(ctx) qq_data[[paste(tool, cohort, ctx)]])
    plot(range(unlist(lapply(q, `[[`, "e"))), range(unlist(lapply(q, `[[`, "o"))), type = "n",
         xlab = expression(Expected ~ -log[10](italic(p))), ylab = expression(Observed ~ -log[10](italic(p))))
    abline(0, 1, col = "grey30", lty = 2)
    for (k in seq_along(CONTEXTS)) points(q[[k]]$e, q[[k]]$o, pch = 20, cex = 0.35, col = CTX_COLORS[[CONTEXTS[k]]])
    mtext(letters[i], side = 3, adj = -0.18, line = 0.6, font = 2, cex = 0.83)
  }
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend("center", legend = CONTEXTS, pch = 19, col = CTX_COLORS[CONTEXTS], horiz = TRUE, bty = "n", cex = 1)
}
tiff(file.path(OUT, "qq_panel.tiff"), width = 20, height = 21, units = "cm", res = 300, compression = "lzw"); draw_qq(); dev.off()
pdf(file.path(OUT, "qq_panel.pdf"), width = 20 / 2.54, height = 21 / 2.54); draw_qq(); dev.off()

# Values shown in the QQ panel: expected and observed -log10(p) of each point
qq_rows <- list()
i <- 0
for (tool in names(TOOLS)) for (cohort in COHORTS) {
  i <- i + 1
  for (ctx in CONTEXTS) {
    q <- qq_data[[paste(tool, cohort, ctx)]]
    qq_rows[[paste(tool, cohort, ctx)]] <- data.table(panel = letters[i], tool = tool, cohort = cohort, context = ctx,
                                                      expected_log10p = q$e, observed_log10p = q$o)
  }
}
fwrite(rbindlist(qq_rows), file.path(OUT, "qq_panel_values.tsv"), sep = "\t")
msg("Step 13ab done")
