############################################################
# TreeGeneClimate (TGC) — ECS
# Shared functions for Steps 8a and 12ab0 (sourced, not run directly)
#
# Genomic relationship matrix (GRM) with the GCTA estimator from called
# (non-imputed) genotypes. As in GCTA,
# each pair of trees is evaluated over the SNPs genotyped in both, so trees
# with more missing genotypes are not pulled towards the population mean (as
# happens when missing genotypes are imputed or set to their expected value).
#
#   off-diagonal  A_jk = sum_i z_ij z_ik / N_jk,  z_ij = (x_ij - 2p_i) / sqrt(2p_i(1 - p_i))
#   diagonal      A_jj = 1 + sum_i [x_ij^2 - (1 + 2p_i) x_ij + 2p_i^2] / [2p_i(1 - p_i)] / N_jj
#
# Only SNPs called in at least 99% of the cohort's trees are used (GRM_MAX_MISS):
# with about 7 reads per tree, genotype calls of poorly covered trees are biased
# (missing calls, heterozygotes called as homozygotes), and SNPs with many
# missing calls make the GRM reflect sequencing depth rather than relatedness.
#
# The sums are additive over SNPs, so leave-one-chromosome-out (LOCO) GRMs are
# obtained exactly by subtracting the components of one chromosome.
############################################################

GRM_MAX_MISS <- 0.01

# Rows (SNPs) of X kept for the GRM
grm_snps <- function(X, max_miss = GRM_MAX_MISS) rowMeans(is.na(X)) <= max_miss

# Components of the GCTA GRM from a genotype matrix X (SNPs x samples, values
# 0/1/2 or NA). p: allele frequencies (default: from X itself)
gcta_components <- function(X, p = rowMeans(X, na.rm = TRUE) / 2) {
  ok <- is.finite(p) & p > 0 & p < 1
  X <- X[ok, , drop = FALSE]; p <- p[ok]
  het <- 2 * p * (1 - p)
  called <- !is.na(X)
  Z <- (X - 2 * p) / sqrt(het)
  W <- (X^2 - (1 + 2 * p) * X + 2 * p^2) / het
  Z[!called] <- 0; W[!called] <- 0
  storage.mode(called) <- "double"
  list(num = crossprod(Z), cnt = crossprod(called), dnum = colSums(W), n_snps = nrow(X))
}

# Subtract (sign = -1) or add components
gcta_combine <- function(a, b, sign = -1) {
  list(num = a$num + sign * b$num, cnt = a$cnt + sign * b$cnt,
       dnum = a$dnum + sign * b$dnum, n_snps = a$n_snps + sign * b$n_snps)
}

# GRM from components; ids: sample names
gcta_grm <- function(comp, ids) {
  A <- comp$num / comp$cnt
  diag(A) <- 1 + comp$dnum / diag(comp$cnt)
  dimnames(A) <- list(ids, ids)
  A
}
