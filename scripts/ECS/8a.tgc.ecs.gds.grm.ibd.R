#!/usr/bin/env Rscript
############################################################
# TreeGeneClimate (TGC) — ECS
# Step 8a: Genotype GDS files, genomic relationships, IBD and PCA import
#
#   1) Convert the cohort VCFs of Step 6a to SNPRelate GDS: called genotypes
#      (before imputation) and BEAGLE-imputed genotypes (best-guess calls 0/1/2)
#   2) Genomic relationship matrix (GRM) per cohort with the GCTA estimator
#      from the called genotypes of the within-cohort polymorphic SNPs called
#      in >= 99% of the trees, each pair of trees over the SNPs genotyped in
#      both (tgc.ecs.functions.R); kinship K = GRM / 2
#   3) Import PLINK IBD (Step 7a) and PLINK PCA on LD-pruned SNPs (Step 5a)
#   4) Save R objects and tables used by Steps 9a, 10a, 11ab, 12ab0 and 16ab
#
# The GRM is computed from called rather than imputed genotypes, and from SNPs
# with few missing calls: with a mean ECS depth of about 7 reads per tree, the
# share of missing calls varies among trees (0.8-36%), and imputed genotypes or
# SNPs with many missing calls make poorly covered trees appear related to all
# others.
#
# INPUTS
#   RESULTS/ECS/VCF_SPLIT/tgc.ecs.<cohort>.call.filt.maf05.snvs.poly.vcf.gz
#   RESULTS/ECS/VCF_SPLIT/tgc.ecs.<cohort>.call.filt.maf05.snvs.poly.imputed.vcf.gz
#   RESULTS/ECS/POPGEN/RELATEDNESS/IBD/tgc.ecs.<cohort>.pruned.ibd.genome
#   RESULTS/ECS/POPGEN/STRUCTURE/<COHORT>/PCA/tgc.ecs.<cohort>.pca.pruned.{eigenvec,eigenval}
#   DATA/METADATA/breeding_sample2family.txt, natural_sample2pop.txt
#
# OUTPUTS (RESULTS/ECS/RANALYSIS/)
#   RDATA/<cohort>.called.snp.gds, <cohort>.imputed.snp.gds
#   RDATA/<cohort>_grm_gcta.rds, <cohort>_kinship_gcta_half.rds
#   RDATA/<cohort>_sample_annotation.rds
#   RDATA/<cohort>_ibd_plink_long.rds, <cohort>_ibd_plink_mat.rds
#   RDATA/<cohort>_pca_pruned.rds
#   TABLES/{grm_kinship,ibd,pca}/*.tsv (grm_kinship/grm_snp_counts.tsv: SNPs used for the GRM)
#
# R packages: SNPRelate, gdsfmt, data.table
############################################################

suppressPackageStartupMessages({
  library(SNPRelate)
  library(gdsfmt)
  library(data.table)
})
options(stringsAsFactors = FALSE)

# === USER CONFIGURATION ===
PROJECT_ROOT <- Sys.getenv("TGC_PROJECT_ROOT", unset = "/path/to/your/project")
# ===========================
source(file.path(PROJECT_ROOT, "SCRIPTS/ECS/tgc.ecs.functions.R"))

VCF_SPLIT_DIR <- file.path(PROJECT_ROOT, "RESULTS/ECS/VCF_SPLIT")
STRUCT_DIR    <- file.path(PROJECT_ROOT, "RESULTS/ECS/POPGEN/STRUCTURE")
IBD_DIR       <- file.path(PROJECT_ROOT, "RESULTS/ECS/POPGEN/RELATEDNESS/IBD")
META_DIR      <- file.path(PROJECT_ROOT, "DATA/METADATA")
RANA_DIR      <- file.path(PROJECT_ROOT, "RESULTS/ECS/RANALYSIS")
RDATA_DIR     <- file.path(RANA_DIR, "RDATA")
TAB_DIR       <- file.path(RANA_DIR, "TABLES")
for (d in c(RDATA_DIR, file.path(TAB_DIR, c("grm_kinship", "ibd", "pca"))))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

COHORTS <- list(
  breeding = list(map = "breeding_sample2family.txt", group = "Family"),
  natural  = list(map = "natural_sample2pop.txt",     group = "Population")
)

msg <- function(...) cat("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n", sep = "")
ensure_file <- function(x) if (!file.exists(x)) stop("Missing file: ", x)

vcf2gds <- function(vcf, gds) {
  ensure_file(vcf)
  if (file.exists(gds)) file.remove(gds)
  snpgdsVCF2GDS(vcf.fn = vcf, out.fn = gds, method = "biallelic.only", snpfirstdim = TRUE, verbose = FALSE)
}

read_plink_genome <- function(file, ids) {
  ensure_file(file)
  g <- fread(file, select = c("IID1", "IID2", "PI_HAT"))
  mat <- matrix(NA_real_, length(ids), length(ids), dimnames = list(ids, ids))
  diag(mat) <- 1
  ok <- g$IID1 %in% ids & g$IID2 %in% ids
  mat[cbind(g$IID1[ok], g$IID2[ok])] <- g$PI_HAT[ok]
  mat[cbind(g$IID2[ok], g$IID1[ok])] <- g$PI_HAT[ok]
  list(long = g, mat = mat)
}

grm_counts <- list()
for (cohort in names(COHORTS)) {
  cfg <- COHORTS[[cohort]]
  COH <- toupper(cohort)
  msg("=== ", COH, " ===")

  # 1) VCF -> GDS (called and imputed genotypes)
  stem <- file.path(VCF_SPLIT_DIR, sprintf("tgc.ecs.%s.call.filt.maf05.snvs.poly", cohort))
  gds_called <- file.path(RDATA_DIR, sprintf("%s.called.snp.gds", cohort))
  vcf2gds(paste0(stem, ".vcf.gz"), gds_called)
  vcf2gds(paste0(stem, ".imputed.vcf.gz"), file.path(RDATA_DIR, sprintf("%s.imputed.snp.gds", cohort)))

  # 2) GRM (GCTA estimator, called genotypes)
  map <- fread(file.path(META_DIR, cfg$map), header = FALSE, col.names = c("IID", "group"))
  gf  <- snpgdsOpen(gds_called)
  ids <- read.gdsn(index.gdsn(gf, "sample.id"))
  geno <- snpgdsGetGeno(gf, snpfirstdim = TRUE, with.id = FALSE, verbose = FALSE)
  snpgdsClose(gf)
  mode(geno) <- "numeric"
  geno[geno < 0 | geno > 2] <- NA
  grp <- map$group[match(ids, map$IID)]
  grp[is.na(grp)] <- "Unknown"
  msg(nrow(geno), " SNPs x ", ncol(geno), " samples; missing calls: ",
      round(100 * mean(is.na(geno)), 2), "%")
  keep <- grm_snps(geno)
  comp <- gcta_components(geno[keep, , drop = FALSE])
  G <- gcta_grm(comp, ids)
  msg("GRM from ", comp$n_snps, " SNPs; mean diagonal ", round(mean(diag(G)), 3))
  grm_counts[[cohort]] <- data.table(cohort = cohort, n_snps_polymorphic = nrow(geno), n_snps = comp$n_snps)
  rm(geno, comp); gc()

  annot <- data.table(IID = ids, group = grp)
  setnames(annot, "group", cfg$group)
  saveRDS(G,     file.path(RDATA_DIR, sprintf("%s_grm_gcta.rds", cohort)))
  saveRDS(G / 2, file.path(RDATA_DIR, sprintf("%s_kinship_gcta_half.rds", cohort)))
  saveRDS(annot, file.path(RDATA_DIR, sprintf("%s_sample_annotation.rds", cohort)))
  fwrite(data.table(IID = ids, G), file.path(TAB_DIR, "grm_kinship", sprintf("%s_grm_gcta.tsv", cohort)), sep = "\t")
  fwrite(annot, file.path(TAB_DIR, "grm_kinship", sprintf("%s_sample_annotation.tsv", cohort)), sep = "\t")

  # 3a) IBD
  ibd <- read_plink_genome(file.path(IBD_DIR, sprintf("tgc.ecs.%s.pruned.ibd.genome", cohort)), ids)
  saveRDS(ibd$long, file.path(RDATA_DIR, sprintf("%s_ibd_plink_long.rds", cohort)))
  saveRDS(ibd$mat,  file.path(RDATA_DIR, sprintf("%s_ibd_plink_mat.rds", cohort)))
  fwrite(ibd$long, file.path(TAB_DIR, "ibd", sprintf("%s_ibd_plink_long.tsv", cohort)), sep = "\t")
  msg("IBD pairs: ", nrow(ibd$long))

  # 3b) PCA (PLINK, LD-pruned SNPs)
  pca_prefix <- file.path(STRUCT_DIR, COH, "PCA", sprintf("tgc.ecs.%s.pca.pruned", cohort))
  ensure_file(paste0(pca_prefix, ".eigenvec"))
  ev <- fread(paste0(pca_prefix, ".eigenvec"), header = FALSE)
  setnames(ev, c("FID", "IID", paste0("PC", seq_len(ncol(ev) - 2))))
  eval <- scan(paste0(pca_prefix, ".eigenval"), quiet = TRUE)
  saveRDS(list(scores = ev, eigenval = eval), file.path(RDATA_DIR, sprintf("%s_pca_pruned.rds", cohort)))
  fwrite(ev, file.path(TAB_DIR, "pca", sprintf("%s_pca_pruned_scores.tsv", cohort)), sep = "\t")
  fwrite(data.table(PC = seq_along(eval), eigenval = eval),
         file.path(TAB_DIR, "pca", sprintf("%s_pca_pruned_eigenval.tsv", cohort)), sep = "\t")
}

fwrite(rbindlist(grm_counts), file.path(TAB_DIR, "grm_kinship", "grm_snp_counts.tsv"), sep = "\t")
msg("Step 8a done")
sessionInfo()
