#!/usr/bin/env Rscript
############################################################
# TreeGeneClimate (TGC) — JOINT ECS + TMS
# Step 12ab1: cis-meQTL mapping with MatrixEQTL (one cohort x context panel)
#
# MODEL   M-value ~ SNP + PC1..PC10   (ordinary least squares)
# CIS     SNP within 100 kb of the methylation site
# OUTPUT  all cis pairs; Benjamini-Hochberg FDR over all cis pairs of the panel
#
# USAGE
#   TGC_COHORT=BREEDING TGC_CONTEXT=CpG Rscript 12ab1.tgc.joint.matrixeqtl.mapping.R
#
# INPUTS (Step 12ab0; RESULTS/JOINT/MEQTL/INPUTS/<COHORT>/)
#   genotypes.rds, snp_variant_annot.tsv
#   <CTX>/methylation_mvalues_matrix.rds, methylation_site_annot.tsv,
#   <CTX>/shared_sample_ids.tsv, pcs_shared.tsv
#
# OUTPUTS (RESULTS/JOINT/MEQTL/MATRIXEQTL/<COHORT>/<CTX>/)
#   cis_meqtl_all_results.rds / .tsv.gz
#     snp site snp_chr snp_pos site_pos distance beta statistic pvalue p_FDR
#   cis_meqtl_significant.tsv   (FDR < 0.05)
#   cis_meqtl_top_per_site.tsv
#   cis_meqtl_panel_summary.tsv (includes genomic inflation lambda)
#   qq_<cohort>_<ctx>.{tiff,pdf}
#
# R packages: MatrixEQTL, data.table
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(MatrixEQTL)
})
options(stringsAsFactors = FALSE)

# === USER CONFIGURATION ===
PROJECT_ROOT <- Sys.getenv("TGC_PROJECT_ROOT", unset = "/path/to/your/project")
CIS_WINDOW <- 1e5
FDR_THRESH <- 0.05
# ===========================

COHORT  <- Sys.getenv("TGC_COHORT");  CONTEXT <- Sys.getenv("TGC_CONTEXT")
stopifnot(COHORT %in% c("BREEDING", "NATURAL"), CONTEXT %in% c("CpG", "CHG", "CHH"))
IN_C  <- file.path(PROJECT_ROOT, "RESULTS/JOINT/MEQTL/INPUTS", COHORT)
IN_P  <- file.path(IN_C, CONTEXT)
OUT_P <- file.path(PROJECT_ROOT, "RESULTS/JOINT/MEQTL/MATRIXEQTL", COHORT, CONTEXT)
dir.create(OUT_P, recursive = TRUE, showWarnings = FALSE)
tag <- paste0(tolower(COHORT), "_", tolower(CONTEXT))
msg <- function(...) cat("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n", sep = "")

lambda_gc <- function(p) {
  p <- p[is.finite(p) & p > 0 & p <= 1]
  median(qchisq(p, df = 1, lower.tail = FALSE)) / qchisq(0.5, df = 1)
}

# ---- inputs ----
shared <- fread(file.path(IN_P, "shared_sample_ids.tsv"))$sample_id
M      <- readRDS(file.path(IN_P, "methylation_mvalues_matrix.rds"))[shared, , drop = FALSE]
sites  <- fread(file.path(IN_P, "methylation_site_annot.tsv"))
pcs    <- fread(file.path(IN_P, "pcs_shared.tsv"))
va     <- fread(file.path(IN_C, "snp_variant_annot.tsv"))
geno   <- readRDS(file.path(IN_C, "genotypes.rds"))[, shared, drop = FALSE]
stopifnot(identical(colnames(M), sites$methyl_loc), identical(pcs$sample_id, shared),
          identical(rownames(geno), as.character(va$snp_id)))
msg(COHORT, "/", CONTEXT, ": ", length(shared), " samples, ", ncol(M), " sites, ", nrow(geno), " SNPs")

snps <- SlicedData$new(); snps$CreateFromMatrix(geno * 1.0)
gene <- SlicedData$new(); gene$CreateFromMatrix(t(M))
cvrt <- SlicedData$new()
cvrt$CreateFromMatrix(t(as.matrix(pcs[, grep("^PC[0-9]+$", names(pcs)), with = FALSE])))
rm(geno); gc()

snpspos <- data.frame(snp = as.character(va$snp_id), chr = va$chr, pos = va$pos)
genepos <- data.frame(gene = sites$methyl_loc, chr = sites$chr, start = sites$start, end = sites$end)

# ---- mapping ----
me <- Matrix_eQTL_main(
  snps = snps, gene = gene, cvrt = cvrt,
  output_file_name = "", pvOutputThreshold = 0,
  useModel = modelLINEAR, errorCovariance = numeric(), verbose = FALSE,
  output_file_name.cis = "", pvOutputThreshold.cis = 1,
  snpspos = snpspos, genepos = genepos, cisDist = CIS_WINDOW,
  pvalue.hist = FALSE, min.pv.by.genesnp = FALSE, noFDRsaveMemory = FALSE)

res <- as.data.table(me$cis$eqtls)
setnames(res, c("snps", "gene"), c("snp", "site"))
res[, snp := as.character(snp)]; res[, site := as.character(site)]
res[, p_FDR := p.adjust(pvalue, method = "BH")]
res[, FDR := NULL]
res[, `:=`(snp_chr = va$chr[match(snp, as.character(va$snp_id))],
           snp_pos = va$pos[match(snp, as.character(va$snp_id))],
           site_pos = sites$pos[match(site, sites$methyl_loc)])]
res[, distance := abs(site_pos - snp_pos)]
stopifnot(max(res$distance) <= CIS_WINDOW)
setcolorder(res, c("snp", "site", "snp_chr", "snp_pos", "site_pos", "distance",
                   "beta", "statistic", "pvalue", "p_FDR"))

lam <- lambda_gc(res$pvalue)
sig <- res[p_FDR < FDR_THRESH][order(p_FDR)]
top <- res[res[, .I[which.min(pvalue)], by = site]$V1]
msg("cis pairs: ", nrow(res), " | FDR < ", FDR_THRESH, ": ", nrow(sig), " | lambda = ", round(lam, 4))

saveRDS(res, file.path(OUT_P, "cis_meqtl_all_results.rds"))
fwrite(res, file.path(OUT_P, "cis_meqtl_all_results.tsv.gz"), sep = "\t")
fwrite(sig, file.path(OUT_P, "cis_meqtl_significant.tsv"), sep = "\t")
fwrite(top, file.path(OUT_P, "cis_meqtl_top_per_site.tsv"), sep = "\t")
fwrite(data.table(tool = "MatrixEQTL", cohort = COHORT, context = CONTEXT,
                  n_samples = length(shared), n_sites = ncol(M), n_snps = nrow(va),
                  n_sites_with_cis = uniqueN(res$site), n_snps_with_cis = uniqueN(res$snp),
                  n_pairs = nrow(res), n_fdr05 = nrow(sig),
                  n_fdr1e10 = sum(res$p_FDR < 1e-10), lambda = round(lam, 4)),
       file.path(OUT_P, "cis_meqtl_panel_summary.tsv"), sep = "\t")

# ---- QQ plot (TIFF + PDF) ----
pv <- sort(res$pvalue[res$pvalue > 0])
qq <- function() {
  plot(-log10(ppoints(length(pv))), -log10(pv), pch = 20, cex = 0.4, col = "#2b8cbe",
       xlab = expression(Expected ~ -log[10](italic(p))), ylab = expression(Observed ~ -log[10](italic(p))))
  abline(0, 1, col = "red", lty = 2)
  mtext(sprintf("lambda = %.3f", lam), side = 3, adj = 0, line = 0.3)
}
tiff(file.path(OUT_P, sprintf("qq_%s.tiff", tag)), width = 2400, height = 2400, res = 300, compression = "lzw"); qq(); dev.off()
pdf(file.path(OUT_P, sprintf("qq_%s.pdf", tag)), width = 8, height = 8); qq(); dev.off()
msg("Step 12ab1 done")
