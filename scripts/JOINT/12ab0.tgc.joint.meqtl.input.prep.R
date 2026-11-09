#!/usr/bin/env Rscript
############################################################
# TreeGeneClimate (TGC) — JOINT ECS + TMS
# Step 12ab0: Inputs for cis-meQTL mapping (both cohorts, three contexts)
#
# GENOTYPES (per cohort, shared by Steps 12ab1 and 12ab2)
#   - BEAGLE-imputed genotypes (Step 8a GDS; best-guess calls 0/1/2)
#   - SNPs with within-cohort MAF >= 0.05 and missing rate <= 0.20
#   - Residual missing calls (variants on single-variant contigs, which
#     BEAGLE does not impute) set to the rounded expected genotype from the
#     allele frequency of the tree's family (breeding) or stand (natural).
#     Both mapping tools read the same integer genotype matrix.
# COVARIATES
#   - PC1-PC10 of the PCA on LD-pruned SNPs (Step 5a/8a)
# RELATEDNESS (GENESIS)
#   - GCTA-estimator GRMs from called (non-imputed) genotypes of SNPs called
#     in >= 99% of the trees, as in Step 8a
#   - BREEDING: leave-one-chromosome-out (LOCO) GRMs, all SNPs except those on
#     the tested chromosome/scaffold (exact, by subtracting its components)
#   - NATURAL: genome-wide GRM (Step 8a)
# METHYLATION (per cohort x context)
#   - Step 5b3 datasets; sites with zero variance removed
#   - M = log2[(beta x 100 + 0.5) / (100 - beta x 100 + 0.5)]
#
# INPUTS
#   RESULTS/ECS/RANALYSIS/RDATA/<cohort>.imputed.snp.gds, <cohort>.called.snp.gds
#   RESULTS/ECS/RANALYSIS/RDATA/<cohort>_grm_gcta.rds
#   RESULTS/ECS/RANALYSIS/RDATA/<cohort>_pca_pruned.rds
#   RESULTS/ECS/RANALYSIS/RDATA/<cohort>_sample_annotation.rds
#   RESULTS/TMS/RANALYSIS/METHYLKIT_OBJECTS_MEQTL/methylBase_<cohort>_<ctx>_meqtl.rds
#
# OUTPUTS (RESULTS/JOINT/MEQTL/INPUTS/<COHORT>/)
#   genotypes.rds                 integer matrix, SNPs x samples
#   snp_variant_annot.tsv         snp_id chr pos ref alt maf missing_rate
#   loco_grm/loco_grm_excl_<chr>.rds   (BREEDING)
#   <CTX>/methylation_mvalues_matrix.rds   samples x sites
#   <CTX>/methylation_site_annot.tsv       methyl_loc chr start end strand pos
#   <CTX>/shared_sample_ids.tsv, pcs_shared.tsv, grm_shared.rds
#   ../step12ab0_summary.tsv
#
# R packages: SNPRelate, gdsfmt, methylKit, data.table
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(SNPRelate)
  library(gdsfmt)
  library(methylKit)
})
options(stringsAsFactors = FALSE)

# === USER CONFIGURATION ===
PROJECT_ROOT <- Sys.getenv("TGC_PROJECT_ROOT", unset = "/path/to/your/project")
N_PCS    <- 10L
MIN_MAF  <- 0.05
MAX_MISS <- 0.20
EPS_M    <- 0.5
# ===========================
source(file.path(PROJECT_ROOT, "SCRIPTS/ECS/tgc.ecs.functions.R"))

ECS_RDATA <- file.path(PROJECT_ROOT, "RESULTS/ECS/RANALYSIS/RDATA")
METH_DIR  <- file.path(PROJECT_ROOT, "RESULTS/TMS/RANALYSIS/METHYLKIT_OBJECTS_MEQTL")
OUT_ROOT  <- file.path(PROJECT_ROOT, "RESULTS/JOINT/MEQTL/INPUTS")
dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)
COHORTS  <- c("BREEDING", "NATURAL")
CONTEXTS <- c("CpG", "CHG", "CHH")
GROUP_COL <- c(BREEDING = "Family", NATURAL = "Population")

msg <- function(...) cat("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n", sep = "")

# Missing genotypes (SNPs x samples) -> 2 x group allele frequency
fill_by_group <- function(geno, grp, round_values = FALSE) {
  p_all <- rowMeans(geno, na.rm = TRUE) / 2
  for (g in unique(grp)) {
    idx  <- which(grp == g)
    sub  <- geno[, idx, drop = FALSE]
    miss <- is.na(sub)
    if (!any(miss)) next
    pg <- rowMeans(sub, na.rm = TRUE) / 2
    pg[!is.finite(pg)] <- p_all[!is.finite(pg)]
    fill <- 2 * pg
    if (round_values) fill <- round(fill)
    sub[miss] <- fill[row(sub)[miss]]
    geno[, idx] <- sub
  }
  geno
}

# Methylation % (sites x samples) from a methylBase object
perc_methylation <- function(mb) {
  d  <- getData(mb)
  cs <- as.matrix(d[, grep("^numCs[0-9]+$", names(d))])
  ts <- as.matrix(d[, grep("^numTs[0-9]+$", names(d))])
  p  <- 100 * cs / (cs + ts)
  p[!is.finite(p)] <- NA_real_
  colnames(p) <- mb@sample.ids
  p
}

summary_rows <- list()

for (cohort in COHORTS) {
  coh <- tolower(cohort)
  out_c <- file.path(OUT_ROOT, cohort)
  dir.create(out_c, recursive = TRUE, showWarnings = FALSE)
  msg("=== ", cohort, " ===")

  # ---------------- genotypes ----------------
  annot <- as.data.table(readRDS(file.path(ECS_RDATA, sprintf("%s_sample_annotation.rds", coh))))
  gf <- snpgdsOpen(file.path(ECS_RDATA, sprintf("%s.imputed.snp.gds", coh)))
  samples <- read.gdsn(index.gdsn(gf, "sample.id"))
  va <- data.table(snp_id = read.gdsn(index.gdsn(gf, "snp.id")),
                   chr    = as.character(read.gdsn(index.gdsn(gf, "snp.chromosome"))),
                   pos    = as.integer(read.gdsn(index.gdsn(gf, "snp.position"))),
                   allele = read.gdsn(index.gdsn(gf, "snp.allele")))
  geno <- snpgdsGetGeno(gf, snpfirstdim = TRUE, with.id = FALSE, verbose = FALSE)
  snpgdsClose(gf)
  mode(geno) <- "numeric"
  geno[geno < 0 | geno > 2] <- NA
  grp <- annot[[GROUP_COL[[cohort]]]][match(samples, annot$IID)]
  grp[is.na(grp)] <- "Unknown"

  va[, missing_rate := rowMeans(is.na(geno))]
  af <- rowMeans(geno, na.rm = TRUE) / 2
  va[, maf := pmin(af, 1 - af)]
  va[, c("ref", "alt") := tstrsplit(allele, "/", fixed = TRUE)]
  keep <- is.finite(va$maf) & va$maf >= MIN_MAF & va$missing_rate <= MAX_MISS
  msg(nrow(va), " imputed SNPs; kept ", sum(keep), " (MAF >= ", MIN_MAF, ", missing <= ", MAX_MISS, ")")

  # LOCO GRMs (breeding): GCTA estimator on called genotypes, as the full GRM of Step 8a
  if (cohort == "BREEDING") {
    msg("LOCO GRMs ...")
    gc_ <- snpgdsOpen(file.path(ECS_RDATA, sprintf("%s.called.snp.gds", coh)))
    stopifnot(identical(read.gdsn(index.gdsn(gc_, "sample.id")), samples))
    chr_c <- as.character(read.gdsn(index.gdsn(gc_, "snp.chromosome")))
    Xc <- snpgdsGetGeno(gc_, snpfirstdim = TRUE, with.id = FALSE, verbose = FALSE)
    snpgdsClose(gc_)
    mode(Xc) <- "numeric"
    Xc[Xc < 0 | Xc > 2] <- NA
    keep_g <- grm_snps(Xc)                       # SNPs called in >= 99% of trees, as Step 8a
    Xc <- Xc[keep_g, , drop = FALSE]; chr_c <- chr_c[keep_g]
    total <- gcta_components(Xc)
    loco_dir <- file.path(out_c, "loco_grm")
    dir.create(loco_dir, showWarnings = FALSE)
    for (ch in unique(chr_c)) {
      part <- gcta_components(Xc[chr_c == ch, , drop = FALSE])
      saveRDS(gcta_grm(gcta_combine(total, part), samples),
              file.path(loco_dir, sprintf("loco_grm_excl_%s.rds", ch)))
    }
    msg("LOCO GRMs written for ", length(unique(chr_c)), " chromosomes/scaffolds")
    rm(Xc, total, part); gc()
  }

  geno <- geno[keep, , drop = FALSE]
  va   <- va[keep]
  geno <- fill_by_group(geno, grp, round_values = TRUE)
  storage.mode(geno) <- "integer"
  dimnames(geno) <- list(as.character(va$snp_id), samples)
  saveRDS(geno, file.path(out_c, "genotypes.rds"))
  fwrite(va[, .(snp_id, chr, pos, ref, alt, maf = round(maf, 5), missing_rate = round(missing_rate, 5))],
         file.path(out_c, "snp_variant_annot.tsv"), sep = "\t")
  rm(geno); gc()

  # ---------------- covariates and GRM ----------------
  pca <- readRDS(file.path(ECS_RDATA, sprintf("%s_pca_pruned.rds", coh)))
  pcs <- as.data.table(pca$scores)[, c("IID", paste0("PC", seq_len(N_PCS))), with = FALSE]
  setnames(pcs, "IID", "sample_id")
  G_full <- readRDS(file.path(ECS_RDATA, sprintf("%s_grm_gcta.rds", coh)))

  # ---------------- methylation panels ----------------
  for (ctx in CONTEXTS) {
    out_p <- file.path(out_c, ctx)
    dir.create(out_p, recursive = TRUE, showWarnings = FALSE)
    mb <- readRDS(file.path(METH_DIR, sprintf("methylBase_%s_%s_meqtl.rds", coh, tolower(ctx))))
    d  <- getData(mb)
    pm <- perc_methylation(mb)
    n_united <- nrow(pm)
    keep_var <- rowSums(!is.na(pm)) > 1 & apply(pm, 1, sd, na.rm = TRUE) > 0
    pm <- pm[keep_var, , drop = FALSE]
    site_annot <- data.table(chr = as.character(d$chr), start = as.integer(d$start),
                             end = as.integer(d$end), strand = as.character(d$strand))[keep_var]
    site_annot[, pos := start]
    site_annot[, methyl_loc := paste0(chr, ":", start, "-", end)]
    M <- t(log2((pm + EPS_M) / (100 - pm + EPS_M)))   # samples x sites
    colnames(M) <- site_annot$methyl_loc

    shared <- sort(Reduce(intersect, list(samples, pcs$sample_id, rownames(G_full), rownames(M))))
    M <- M[shared, , drop = FALSE]
    saveRDS(M, file.path(out_p, "methylation_mvalues_matrix.rds"))
    fwrite(site_annot[, .(methyl_loc, chr, start, end, strand, pos)],
           file.path(out_p, "methylation_site_annot.tsv"), sep = "\t")
    fwrite(data.table(sample_id = shared), file.path(out_p, "shared_sample_ids.tsv"), sep = "\t")
    fwrite(pcs[match(shared, sample_id)], file.path(out_p, "pcs_shared.tsv"), sep = "\t")
    saveRDS(G_full[shared, shared], file.path(out_p, "grm_shared.rds"))

    summary_rows[[paste(cohort, ctx)]] <- data.table(
      cohort = cohort, context = ctx, n_samples = length(shared),
      n_sites_united = n_united, n_sites_tested = ncol(M), n_snps = nrow(va))
    msg(cohort, "/", ctx, ": ", length(shared), " samples, ", ncol(M), " sites (", n_united, " united)")
    rm(mb, d, pm, M); gc()
  }
}

summ <- rbindlist(summary_rows)
fwrite(summ, file.path(OUT_ROOT, "step12ab0_summary.tsv"), sep = "\t")
print(summ)
msg("Step 12ab0 done")
sessionInfo()
