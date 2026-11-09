#!/usr/bin/env Rscript
############################################################
# TreeGeneClimate (TGC) — JOINT ECS + TMS
# Step 12ab2: cis-meQTL mapping with GENESIS (one cohort x context panel)
#
# MODEL   per methylation site, null linear mixed model
#           M-value ~ PC1..PC10 + g,  g ~ N(0, sigma2_g * GRM)
#         (AI-REML; BREEDING: LOCO GRM excluding the site's chromosome;
#          NATURAL: genome-wide GRM), then a score test for every SNP
#         within 100 kb of the site
# OUTPUT  all cis pairs; Benjamini-Hochberg FDR over all cis pairs of the panel
#
# The panel can be split into chunks of chromosomes/scaffolds that run in
# parallel, followed by a combine run:
#   TGC_COHORT=BREEDING TGC_CONTEXT=CHH TGC_NCHUNKS=8 TGC_CHUNK=3 Rscript 12ab2...R
#   TGC_COHORT=BREEDING TGC_CONTEXT=CHH TGC_NCHUNKS=8 TGC_CHUNK=combine Rscript 12ab2...R
# Without TGC_NCHUNKS the whole panel runs in one job.
#
# INPUTS (Step 12ab0; RESULTS/JOINT/MEQTL/INPUTS/<COHORT>/)
#   genotypes.rds, snp_variant_annot.tsv, loco_grm/ (BREEDING)
#   <CTX>/methylation_mvalues_matrix.rds, methylation_site_annot.tsv,
#   <CTX>/shared_sample_ids.tsv, pcs_shared.tsv, grm_shared.rds
#
# OUTPUTS (RESULTS/JOINT/MEQTL/GENESIS/<COHORT>/<CTX>/)
#   cis_meqtl_all_results.rds / .tsv.gz
#     snp site snp_chr snp_pos site_pos distance beta beta_se statistic
#     pvalue p_FDR freq MAC n_obs
#   cis_meqtl_significant.tsv, cis_meqtl_top_per_site.tsv
#   cis_meqtl_panel_summary.tsv, qq_<cohort>_<ctx>.{tiff,pdf}
#
# R packages: GENESIS, GWASTools, data.table, parallel
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(GENESIS)
  library(GWASTools)
  library(parallel)
})
options(stringsAsFactors = FALSE)

# === USER CONFIGURATION ===
PROJECT_ROOT <- Sys.getenv("TGC_PROJECT_ROOT", unset = "/path/to/your/project")
CIS_WINDOW  <- 1e5
FDR_THRESH  <- 0.05
AIREML_TOL  <- 1e-2
MIN_SAMPLES <- 10L
N_CORES     <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "4"))
# ===========================

COHORT  <- Sys.getenv("TGC_COHORT");  CONTEXT <- Sys.getenv("TGC_CONTEXT")
NCHUNKS <- as.integer(Sys.getenv("TGC_NCHUNKS", unset = "1"))
CHUNK   <- Sys.getenv("TGC_CHUNK", unset = "1")
stopifnot(COHORT %in% c("BREEDING", "NATURAL"), CONTEXT %in% c("CpG", "CHG", "CHH"))
IN_C  <- file.path(PROJECT_ROOT, "RESULTS/JOINT/MEQTL/INPUTS", COHORT)
IN_P  <- file.path(IN_C, CONTEXT)
OUT_P <- file.path(PROJECT_ROOT, "RESULTS/JOINT/MEQTL/GENESIS", COHORT, CONTEXT)
dir.create(file.path(OUT_P, "chunks"), recursive = TRUE, showWarnings = FALSE)
tag <- paste0(tolower(COHORT), "_", tolower(CONTEXT))
msg <- function(...) cat("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n", sep = "")

lambda_gc <- function(p) {
  p <- p[is.finite(p) & p > 0 & p <= 1]
  median(qchisq(p, df = 1, lower.tail = FALSE)) / qchisq(0.5, df = 1)
}

shared <- fread(file.path(IN_P, "shared_sample_ids.tsv"))$sample_id
sites  <- fread(file.path(IN_P, "methylation_site_annot.tsv"))
va     <- fread(file.path(IN_C, "snp_variant_annot.tsv"))

# Assign chromosomes to chunks, balancing the number of sites
chr_sites <- sites[, .N, by = chr][order(-N)]
chr_sites[, chunk := 0L]
load <- integer(NCHUNKS)
for (i in seq_len(nrow(chr_sites))) { k <- which.min(load); chr_sites$chunk[i] <- k; load[k] <- load[k] + chr_sites$N[i] }

if (CHUNK != "combine") {
  k <- as.integer(CHUNK)
  my_chr <- chr_sites[chunk == k, chr]
  msg(COHORT, "/", CONTEXT, " chunk ", k, "/", NCHUNKS, ": ", length(my_chr), " chromosomes, ",
      sum(sites$chr %in% my_chr), " sites, ", N_CORES, " cores")

  M    <- readRDS(file.path(IN_P, "methylation_mvalues_matrix.rds"))[shared, , drop = FALSE]
  pcs  <- fread(file.path(IN_P, "pcs_shared.tsv"))
  geno <- readRDS(file.path(IN_C, "genotypes.rds"))[, shared, drop = FALSE]
  stopifnot(identical(colnames(M), sites$methyl_loc), identical(pcs$sample_id, shared))

  # GWASTools genotype container (SNPs x samples, integer 0/1/2)
  va[, snp_int := .I]
  chr_levels <- sort(unique(va$chr))
  scan_int <- seq_along(shared)
  reader <- MatrixGenotypeReader(genotype = geno, snpID = va$snp_int,
                                 chromosome = match(va$chr, chr_levels),
                                 position = va$pos, scanID = scan_int)
  genoData <- GenotypeData(reader, scanAnnot = ScanAnnotationDataFrame(data.frame(scanID = scan_int)))
  rm(geno); gc()

  pc_cols <- grep("^PC[0-9]+$", names(pcs), value = TRUE)
  pheno_base <- data.frame(sample.id = scan_int, pcs[, ..pc_cols])
  grm_full <- readRDS(file.path(IN_P, "grm_shared.rds"))[shared, shared]
  prep_grm <- function(G) { G <- G[shared, shared]; dimnames(G) <- list(scan_int, scan_int); diag(G) <- diag(G) + 1e-4; G }
  grm_full <- prep_grm(grm_full)

  out <- list()
  for (ch in my_chr) {
    grm <- grm_full
    if (COHORT == "BREEDING") {
      f <- file.path(IN_C, "loco_grm", sprintf("loco_grm_excl_%s.rds", ch))
      if (file.exists(f)) grm <- prep_grm(readRDS(f))
    }
    s_chr <- sites[chr == ch]
    v_chr <- va[chr == ch]
    if (!nrow(v_chr)) next
    res_chr <- mclapply(seq_len(nrow(s_chr)), function(j) {
      cis <- v_chr[abs(pos - s_chr$pos[j]) <= CIS_WINDOW]
      if (!nrow(cis)) return(NULL)
      y <- M[, s_chr$methyl_loc[j]]
      if (sum(!is.na(y)) < MIN_SAMPLES) return(NULL)
      ph <- pheno_base; ph$y <- y
      nm <- tryCatch(fitNullModel(ph, outcome = "y", covars = pc_cols, cov.mat = grm,
                                  family = "gaussian", AIREML.tol = AIREML_TOL, verbose = FALSE),
                     error = function(e) NULL)
      if (is.null(nm)) return(data.table(site = s_chr$methyl_loc[j], null_failed = TRUE))
      it <- GenotypeBlockIterator(genoData, snpInclude = cis$snp_int)
      a <- tryCatch(assocTestSingle(it, null.model = nm, verbose = FALSE), error = function(e) NULL)
      if (is.null(a) || !nrow(a)) return(NULL)
      a <- as.data.table(a)
      data.table(snp = as.character(va$snp_id[a$variant.id]), site = s_chr$methyl_loc[j],
                 beta = a$Est, beta_se = a$Est.SE, statistic = a$Score.Stat, pvalue = a$Score.pval,
                 freq = a$freq, MAC = a$MAC, n_obs = a$n.obs, null_failed = FALSE)
    }, mc.cores = N_CORES, mc.preschedule = TRUE)
    out[[ch]] <- rbindlist(res_chr, fill = TRUE)
  }
  res <- rbindlist(out, fill = TRUE)
  saveRDS(res, file.path(OUT_P, "chunks", sprintf("chunk_%02d_of_%02d.rds", k, NCHUNKS)))
  msg("chunk ", k, " done: ", nrow(res[null_failed == FALSE]), " pairs; ",
      sum(res$null_failed, na.rm = TRUE), " sites with failed null model")
  if (NCHUNKS > 1) quit(save = "no")
}

# ---- combine (or single-chunk run) ----
files <- list.files(file.path(OUT_P, "chunks"), pattern = sprintf("_of_%02d\\.rds$", NCHUNKS), full.names = TRUE)
stopifnot(length(files) == NCHUNKS)
res <- rbindlist(lapply(files, readRDS), fill = TRUE)
n_null_failed <- sum(res$null_failed, na.rm = TRUE)
res <- res[null_failed == FALSE][, null_failed := NULL]
res[, p_FDR := p.adjust(pvalue, method = "BH")]
res[, `:=`(snp_chr = va$chr[match(snp, as.character(va$snp_id))],
           snp_pos = va$pos[match(snp, as.character(va$snp_id))],
           site_pos = sites$pos[match(site, sites$methyl_loc)])]
res[, distance := abs(site_pos - snp_pos)]
stopifnot(max(res$distance) <= CIS_WINDOW)
setcolorder(res, c("snp", "site", "snp_chr", "snp_pos", "site_pos", "distance",
                   "beta", "beta_se", "statistic", "pvalue", "p_FDR", "freq", "MAC", "n_obs"))

lam <- lambda_gc(res$pvalue)
sig <- res[p_FDR < FDR_THRESH][order(p_FDR)]
top <- res[res[, .I[which.min(pvalue)], by = site]$V1]
msg("cis pairs: ", nrow(res), " | FDR < ", FDR_THRESH, ": ", nrow(sig), " | lambda = ", round(lam, 4))

saveRDS(res, file.path(OUT_P, "cis_meqtl_all_results.rds"))
fwrite(res, file.path(OUT_P, "cis_meqtl_all_results.tsv.gz"), sep = "\t")
fwrite(sig, file.path(OUT_P, "cis_meqtl_significant.tsv"), sep = "\t")
fwrite(top, file.path(OUT_P, "cis_meqtl_top_per_site.tsv"), sep = "\t")
fwrite(data.table(tool = "GENESIS", cohort = COHORT, context = CONTEXT,
                  n_samples = length(shared), n_sites = nrow(sites), n_snps = nrow(va),
                  n_sites_with_cis = uniqueN(res$site), n_snps_with_cis = uniqueN(res$snp),
                  n_sites_null_failed = n_null_failed,
                  n_pairs = nrow(res), n_fdr05 = nrow(sig),
                  n_fdr1e10 = sum(res$p_FDR < 1e-10), lambda = round(lam, 4)),
       file.path(OUT_P, "cis_meqtl_panel_summary.tsv"), sep = "\t")

pv <- sort(res$pvalue[res$pvalue > 0])
qq <- function() {
  plot(-log10(ppoints(length(pv))), -log10(pv), pch = 20, cex = 0.4, col = "#e34a33",
       xlab = expression(Expected ~ -log[10](italic(p))), ylab = expression(Observed ~ -log[10](italic(p))))
  abline(0, 1, col = "red", lty = 2)
  mtext(sprintf("lambda = %.3f", lam), side = 3, adj = 0, line = 0.3)
}
tiff(file.path(OUT_P, sprintf("qq_%s.tiff", tag)), width = 2400, height = 2400, res = 300, compression = "lzw"); qq(); dev.off()
pdf(file.path(OUT_P, sprintf("qq_%s.pdf", tag)), width = 8, height = 8); qq(); dev.off()
msg("Step 12ab2 done")
