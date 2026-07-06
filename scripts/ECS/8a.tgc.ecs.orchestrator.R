############################################################
# TreeGeneClimate (TGC) — ECS
# Step 8a: R Master Orchestrator (CEPH-HDD layout)
#
# PURPOSE
# - One entry-point to:
#   1) Read cohort genotypes (Step 6a imputed VCFs) and metadata
#   2) Convert VCF -> GDS (once) for fast downstream access
#   3) Build GRM + Kinship (VanRaden; K = GRM/2) from imputed genotypes
#   4) Load PLINK IBD (.genome) produced in Step 7a
#   5) Import PCA + ADMIXTURE results from Step 5a (pruned sets)
#   6) Save reusable R objects + tidy tables under RESULTS/ECS/RANALYSIS
#
# NOTES
# - This orchestrator uses GT from the Beagle-imputed VCFs (hard calls).
#   You *can* use DS (dosage) later for GWAS/GRM, but SNPRelate's VCF->GDS
#   path reads GT. With near-zero missingness after imputation,
#   GT-based GRM is acceptable.
# - No plotting here. Plotting lives in separate scripts.
############################################################

suppressPackageStartupMessages({
  library(SNPRelate)  # GDS I/O and SNP utilities
  library(gdsfmt)     # low-level GDS file handle management
  library(dplyr)      # data wrangling
  library(tibble)     # tidy data frames
  library(readr)      # fast TSV output
})

options(stringsAsFactors = FALSE)
set.seed(1)

############################################################
# 0) TOGGLES
# Each block can be run independently once its dependencies exist.
# Set to FALSE to skip a step when rerunning a partial analysis.
############################################################
RUN_GDS_AND_GRM      <- TRUE   # VCF->GDS + GRM/Kinship
RUN_LOAD_IBD         <- TRUE   # load PLINK IBD from step 7a
RUN_IMPORT_PCA       <- TRUE   # read PLINK PCA outputs from step 5a
RUN_IMPORT_ADMIXTURE <- TRUE   # read ADMIXTURE Q/P + CV summaries from step 5a
RUN_SAVE_DAPC_INPUT  <- TRUE   # save a QC-filtered imputed dosage matrix for DAPC (no DAPC run)

# MAF and missingness thresholds applied when preparing the DAPC dosage matrix.
# These are intentionally lenient: DAPC is a descriptive ordination, not a test.
DAPC_MAF_MIN  <- 0.05
DAPC_MISS_MAX <- 0.10

############################################################
# 1) PATHS (CEPH-HDD)
############################################################
# === USER CONFIGURATION ===
PROJECT_ROOT <- "/path/to/your/project"  # <-- set this
# ===========================

# Genotypes and population-genetics outputs from the pipeline
VCF_SPLIT_DIR <- file.path(PROJECT_ROOT, "RESULTS/ECS/VCF_SPLIT")
STRUCT_DIR    <- file.path(PROJECT_ROOT, "RESULTS/ECS/POPGEN/STRUCTURE")
IBD_DIR       <- file.path(PROJECT_ROOT, "RESULTS/ECS/POPGEN/RELATEDNESS/IBD")

# R analysis outputs
RANA_DIR   <- file.path(PROJECT_ROOT, "RESULTS/ECS/RANALYSIS")
RDATA_DIR  <- file.path(RANA_DIR, "RDATA")   # RDS objects consumed by downstream scripts
TABLES_DIR <- file.path(RANA_DIR, "TABLES")  # human-readable TSV exports

dir.create(RDATA_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(TABLES_DIR, recursive = TRUE, showWarnings = FALSE)

# Subdirectories under TABLES/ mirror the analysis blocks below
TAB_GRM   <- file.path(TABLES_DIR, "grm_kinship")
TAB_IBD   <- file.path(TABLES_DIR, "ibd")
TAB_PCA   <- file.path(TABLES_DIR, "pca")
TAB_ADMIX <- file.path(TABLES_DIR, "admixture")
TAB_DAPC  <- file.path(TABLES_DIR, "dapc_inputs")

dir.create(TAB_GRM,   recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_IBD,   recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_PCA,   recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_ADMIX, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DAPC,  recursive = TRUE, showWarnings = FALSE)

# Step 6a outputs (imputed)
VCF_BREED_IMP <- file.path(VCF_SPLIT_DIR, "tgc.ecs.breeding.call.filt.maf05.snvs.poly.imputed.vcf.gz")
VCF_NATUR_IMP <- file.path(VCF_SPLIT_DIR, "tgc.ecs.natural.call.filt.maf05.snvs.poly.imputed.vcf.gz")

# Step 7a outputs (IBD on pruned sets)
IBD_BREED <- file.path(IBD_DIR, "tgc.ecs.breeding.pruned.ibd.genome")
IBD_NATUR <- file.path(IBD_DIR, "tgc.ecs.natural.pruned.ibd.genome")

# Step 5a PCA outputs (pruned)
PCA_BREED_EIGENVEC <- file.path(STRUCT_DIR, "BREEDING", "PCA", "tgc.ecs.breeding.pca.pruned.eigenvec")
PCA_BREED_EIGENVAL <- file.path(STRUCT_DIR, "BREEDING", "PCA", "tgc.ecs.breeding.pca.pruned.eigenval")
PCA_NATUR_EIGENVEC <- file.path(STRUCT_DIR, "NATURAL",  "PCA", "tgc.ecs.natural.pca.pruned.eigenvec")
PCA_NATUR_EIGENVAL <- file.path(STRUCT_DIR, "NATURAL",  "PCA", "tgc.ecs.natural.pca.pruned.eigenval")

# Step 5a ADMIXTURE outputs (pruned bed/bim/fam live in STRUCTURE/*)
ADMIX_BREED_DIR <- file.path(STRUCT_DIR, "BREEDING", "ADMIXTURE")
ADMIX_NATUR_DIR <- file.path(STRUCT_DIR, "NATURAL",  "ADMIXTURE")

# Metadata (expected under PROJECT_ROOT/DATA/METADATA)
META_DIR <- file.path(PROJECT_ROOT, "DATA/METADATA")

# Two-column flat files mapping sample IDs to family/population group labels
BREED_MAP_FILE <- file.path(META_DIR, "breeding_sample2family.txt")
NATUR_MAP_FILE <- file.path(META_DIR, "natural_sample2pop.txt")

PHENO_BREED_FILE <- file.path(META_DIR, "tgc.breeding.phenotypes.txt")
PHENO_NATUR_FILE <- file.path(META_DIR, "tgc.natural.phenotypes.txt")

############################################################
# 2) HELPERS
############################################################
# Halt with a clear message if an expected input is absent
ensure_file <- function(x) if (!file.exists(x)) stop("Missing file: ", x)

# SNPRelate keeps an internal registry of open GDS handles; close any
# dangling ones before opening new files to avoid handle conflicts.
close_open_gds <- function() {
  try({
    lst <- gdsfmt::showfile.gds()
    if (!is.null(lst) && length(lst)) {
      for (i in seq_along(lst)) try(gdsfmt::closefn.gds(lst[[i]]), silent = TRUE)
    }
  }, silent = TRUE)
}

# Convert VCF to GDS only when the GDS does not yet exist (idempotent).
# GDS stores genotypes in a compressed columnar format enabling fast random
# access by sample or SNP ID — essential for large spruce datasets.
vcf_to_gds_if_needed <- function(vcf, gds) {
  if (!file.exists(gds)) {
    message("VCF->GDS: ", basename(vcf), " -> ", basename(gds))
    snpgdsVCF2GDS(vcf.fn = vcf, out.fn = gds, method = "biallelic.only", snpfirstdim = TRUE)
  }
  gds
}

# Read a two-column whitespace-delimited map file (no header) into a
# data frame with standardised column names IID and <col2_name>.
read_map2 <- function(file, col2_name) {
  ensure_file(file)
  df <- read.table(file, header = FALSE, stringsAsFactors = FALSE)
  if (ncol(df) < 2) stop("Map file must have at least 2 columns: ", file)
  df <- df[, 1:2]
  colnames(df) <- c("IID", col2_name)
  df$IID <- as.character(df$IID)
  df[[col2_name]] <- as.character(df[[col2_name]])
  df
}

# Compute the VanRaden (2008) genomic relationship matrix.
# X is an n-samples x p-SNPs matrix of allele counts (0/1/2).
# Each SNP is mean-centred by 2*p_j (expected count under HWE);
# the result is scaled by the sum of heterozygosities so that
# diagonal elements approximate 1 for outbred individuals.
grm_vanraden <- function(X) {
  p <- colMeans(X, na.rm = TRUE) / 2          # allele frequency per SNP
  denom <- 2 * sum(p * (1 - p))               # total expected heterozygosity (scaling factor)
  if (!is.finite(denom) || denom <= 0) stop("Non-positive VanRaden denominator.")
  Xc <- sweep(X, 2, 2 * p, "-")              # centre each SNP column
  tcrossprod(Xc) / denom                      # GRM = Z Z' / denom
}

# Parse a PLINK .genome file (all-pairs IBD estimates) and return
# both a long-format tibble (for plotting) and a square matrix (for mixed models).
read_plink_genome <- function(genome_file, ids = NULL) {
  ensure_file(genome_file)
  df <- read.table(genome_file, header = TRUE, stringsAsFactors = FALSE)
  stopifnot(all(c("IID1", "IID2", "PI_HAT") %in% names(df)))
  long <- tibble(IID1 = df$IID1, IID2 = df$IID2, PI_HAT = df$PI_HAT)
  all_ids <- unique(c(long$IID1, long$IID2))
  if (!is.null(ids)) all_ids <- ids
  # Initialise an n x n matrix; PLINK outputs only the upper triangle,
  # so fill both [a,b] and [b,a] to produce a symmetric matrix.
  mat <- matrix(NA_real_, nrow = length(all_ids), ncol = length(all_ids),
                dimnames = list(all_ids, all_ids))
  diag(mat) <- 1
  for (i in seq_len(nrow(long))) {
    a <- long$IID1[i]; b <- long$IID2[i]; v <- long$PI_HAT[i]
    if (a %in% all_ids && b %in% all_ids) { mat[a, b] <- v; mat[b, a] <- v }
  }
  list(long = long, mat = mat)
}

# Read PLINK PCA output, assigning PC1..PCk column names regardless of
# whether the eigenvec file was written with or without a header row.
read_plink_pca <- function(eigenvec_file, eigenval_file) {
  ensure_file(eigenvec_file); ensure_file(eigenval_file)
  ev <- read.table(eigenvec_file, header = FALSE, stringsAsFactors = FALSE)
  colnames(ev)[1:2] <- c("FID", "IID")
  pcs <- paste0("PC", seq_len(ncol(ev) - 2))
  colnames(ev)[3:ncol(ev)] <- pcs
  eval <- scan(eigenval_file, quiet = TRUE)
  list(scores = as_tibble(ev), eigenval = eval)
}

# Parse the ADMIXTURE cross-validation error log produced by step 5a.
# Each line has the form "CV error (K=<k>): <value>"; regex extracts both fields.
read_cv_errors <- function(file) {
  ensure_file(file)
  x <- readLines(file, warn = FALSE)
  tibble(line = x) |>
    mutate(
      K = as.integer(sub(".*\\(K=([0-9]+)\\).*", "\\1", line)),
      CV = as.numeric(sub(".*:\\s*", "", line))
    ) |>
    select(K, CV) |>
    arrange(K)
}

# Collect the paths and K values of all *.Q files produced by ADMIXTURE,
# capping at max_k to exclude exploratory high-K runs if present.
collect_q_files <- function(structure_dir, admix_dir, prefix, max_k = 30) {
  q <- list.files(structure_dir, pattern = paste0("^", prefix, "\\.[0-9]+\\.Q$"), full.names = TRUE)
  if (!length(q)) return(tibble())
  tibble(Q_file = q) |>
    mutate(K = as.integer(sub(".*\\.(\\d+)\\.Q$", "\\1", Q_file))) |>
    filter(K <= max_k) |>
    arrange(K)
}

# Read a single ADMIXTURE Q file (one row per individual, one column per cluster)
# and attach standard column names Q1..QK.
read_q_matrix <- function(q_file) {
  Q <- as.matrix(read.table(q_file, header = FALSE))
  colnames(Q) <- paste0("Q", seq_len(ncol(Q)))
  Q
}

############################################################
# 3) SANITY CHECKS
# Verify critical inputs before any computation starts so failures
# are caught immediately with informative messages.
############################################################
ensure_file(VCF_BREED_IMP)
ensure_file(VCF_NATUR_IMP)
ensure_file(file.path(STRUCT_DIR, "BREEDING", "tgc.ecs.breeding.admix.pruned.fam"))
ensure_file(file.path(STRUCT_DIR, "NATURAL",  "tgc.ecs.natural.admix.pruned.fam"))

# Group-label maps are optional; the GRM is computed regardless.
# Missing labels appear as "Unknown" in downstream figures.
if (!file.exists(BREED_MAP_FILE)) message("NOTE: breeding map not found at ", BREED_MAP_FILE, " (GRM still runs; group labels will be 'Unknown').")
if (!file.exists(NATUR_MAP_FILE)) message("NOTE: natural map not found at ", NATUR_MAP_FILE, " (GRM still runs; group labels will be 'Unknown').")

############################################################
# 4) GDS + GRM/KINSHIP (from imputed VCFs)
############################################################
if (RUN_GDS_AND_GRM) {
  close_open_gds()

  message("=== BREEDING: GDS + GRM/Kinship (from imputed VCF) ===")
  gds_breed <- file.path(RDATA_DIR, "breeding.imputed.snp.gds")
  vcf_to_gds_if_needed(VCF_BREED_IMP, gds_breed)

  breed_map <- if (file.exists(BREED_MAP_FILE)) read_map2(BREED_MAP_FILE, "Family") else tibble(IID = character(), Family = character())

  # local() confines the GDS file handle and large intermediate matrices to a
  # temporary environment, releasing memory automatically when the block exits.
  local({
    gf <- snpgdsOpen(gds_breed, allow.duplicate = FALSE)
    on.exit(try(snpgdsClose(gf), silent = TRUE), add = TRUE)

    sample_ids <- read.gdsn(index.gdsn(gf, "sample.id"))
    snp_ids    <- read.gdsn(index.gdsn(gf, "snp.id"))

    # Extract hard-call genotypes as a SNP x sample integer matrix (0/1/2).
    # snpfirstdim=TRUE keeps SNPs in rows for efficient column-wise allele
    # frequency calculation; mode() coercion converts to numeric for arithmetic.
    geno <- snpgdsGetGeno(gf, sample.id = sample_ids, snp.id = snp_ids, snpfirstdim = TRUE, with.id = FALSE)
    mode(geno) <- "numeric"
    geno[geno > 2 | geno < 0] <- NA_real_   # guard against unexpected encoding values

    # Assign each sample to its family group; samples absent from the map
    # receive "Unknown" so they remain in the GRM without inflating any group.
    grp <- rep("Unknown", length(sample_ids))
    if (nrow(breed_map) > 0) {
      grp2 <- breed_map$Family[match(sample_ids, breed_map$IID)]
      grp2[is.na(grp2)] <- "Unknown"
      grp <- grp2
    }
    grp <- factor(grp)

    # Drop SNPs with undefined allele frequency (e.g. all-missing loci) before
    # imputing residual missing values.
    p_all <- rowMeans(geno, na.rm = TRUE) / 2
    keep_snp <- is.finite(p_all)
    geno2 <- geno[keep_snp, , drop = FALSE]
    p_all <- p_all[keep_snp]

    # Within-group mean imputation: replace residual missing genotypes with 2*p_g
    # where p_g is the group-specific allele frequency.  This is equivalent to
    # setting the missing call to the within-group expectation and avoids
    # artificially inflating between-family relatedness.
    lev <- levels(grp)
    for (g in lev) {
      idx <- which(grp == g)
      if (!length(idx)) next
      pg <- rowMeans(geno2[, idx, drop = FALSE], na.rm = TRUE) / 2
      pg[!is.finite(pg)] <- p_all[!is.finite(pg)]  # fall back to global freq if group is monomorphic
      miss <- is.na(geno2[, idx, drop = FALSE])
      if (any(miss)) {
        fill <- matrix(2 * pg, nrow = nrow(geno2), ncol = length(idx))
        geno2[, idx][miss] <- fill[miss]
      }
    }
    # Global fallback for any remaining NAs (e.g. samples with no group assignment)
    if (anyNA(geno2)) {
      fill2 <- matrix(2 * p_all, nrow = nrow(geno2), ncol = ncol(geno2))
      geno2[is.na(geno2)] <- fill2[is.na(geno2)]
    }
    geno2[geno2 < 0] <- 0; geno2[geno2 > 2] <- 2  # clip numerical noise

    # Transpose to samples x SNPs for VanRaden function, then drop invariant
    # SNPs (sd == 0) which contribute nothing to the GRM and could cause
    # numerical instability in downstream matrix inversions.
    X <- t(geno2)
    rownames(X) <- sample_ids
    sdv <- apply(X, 2, sd)
    keep_var <- is.finite(sdv) & sdv > 0
    if (!all(keep_var)) X <- X[, keep_var, drop = FALSE]

    G <- grm_vanraden(X); rownames(G) <- colnames(G) <- sample_ids
    K <- G / 2  # kinship matrix: GRM/2 converts to the probability-of-IBD scale

    saveRDS(G, file.path(RDATA_DIR, "breeding_grm_vanraden.rds"))
    saveRDS(K, file.path(RDATA_DIR, "breeding_kinship_vanraden_half.rds"))
    write_tsv(as.data.frame(G) |> rownames_to_column("IID"), file.path(TAB_GRM, "breeding_grm_vanraden.tsv"))
    write_tsv(as.data.frame(K) |> rownames_to_column("IID"), file.path(TAB_GRM, "breeding_kinship_vanraden_half.tsv"))

    annot <- tibble(IID = sample_ids, Family = as.character(grp))
    saveRDS(annot, file.path(RDATA_DIR, "breeding_sample_annotation.rds"))
    write_tsv(annot, file.path(TAB_GRM, "breeding_sample_annotation.tsv"))

    message("BREEDING GRM dims: ", nrow(G), " x ", ncol(G))
  })

  message("=== NATURAL: GDS + GRM/Kinship (from imputed VCF) ===")
  gds_natur <- file.path(RDATA_DIR, "natural.imputed.snp.gds")
  vcf_to_gds_if_needed(VCF_NATUR_IMP, gds_natur)

  natur_map <- if (file.exists(NATUR_MAP_FILE)) read_map2(NATUR_MAP_FILE, "Population") else tibble(IID = character(), Population = character())

  # Identical workflow to BREEDING above; local() again scopes memory to this block.
  local({
    gf <- snpgdsOpen(gds_natur, allow.duplicate = FALSE)
    on.exit(try(snpgdsClose(gf), silent = TRUE), add = TRUE)

    sample_ids <- read.gdsn(index.gdsn(gf, "sample.id"))
    snp_ids    <- read.gdsn(index.gdsn(gf, "snp.id"))

    geno <- snpgdsGetGeno(gf, sample.id = sample_ids, snp.id = snp_ids, snpfirstdim = TRUE, with.id = FALSE)
    mode(geno) <- "numeric"
    geno[geno > 2 | geno < 0] <- NA_real_

    # Map sample IDs to Finnish stand (population) labels
    grp <- rep("Unknown", length(sample_ids))
    if (nrow(natur_map) > 0) {
      grp2 <- natur_map$Population[match(sample_ids, natur_map$IID)]
      grp2[is.na(grp2)] <- "Unknown"
      grp <- grp2
    }
    grp <- factor(grp)

    p_all <- rowMeans(geno, na.rm = TRUE) / 2
    keep_snp <- is.finite(p_all)
    geno2 <- geno[keep_snp, , drop = FALSE]
    p_all <- p_all[keep_snp]

    # Within-population mean imputation (same rationale as BREEDING block above)
    lev <- levels(grp)
    for (g in lev) {
      idx <- which(grp == g)
      if (!length(idx)) next
      pg <- rowMeans(geno2[, idx, drop = FALSE], na.rm = TRUE) / 2
      pg[!is.finite(pg)] <- p_all[!is.finite(pg)]
      miss <- is.na(geno2[, idx, drop = FALSE])
      if (any(miss)) {
        fill <- matrix(2 * pg, nrow = nrow(geno2), ncol = length(idx))
        geno2[, idx][miss] <- fill[miss]
      }
    }
    if (anyNA(geno2)) {
      fill2 <- matrix(2 * p_all, nrow = nrow(geno2), ncol = ncol(geno2))
      geno2[is.na(geno2)] <- fill2[is.na(geno2)]
    }
    geno2[geno2 < 0] <- 0; geno2[geno2 > 2] <- 2

    X <- t(geno2)
    rownames(X) <- sample_ids
    sdv <- apply(X, 2, sd)
    keep_var <- is.finite(sdv) & sdv > 0
    if (!all(keep_var)) X <- X[, keep_var, drop = FALSE]

    G <- grm_vanraden(X); rownames(G) <- colnames(G) <- sample_ids
    K <- G / 2

    saveRDS(G, file.path(RDATA_DIR, "natural_grm_vanraden.rds"))
    saveRDS(K, file.path(RDATA_DIR, "natural_kinship_vanraden_half.rds"))
    write_tsv(as.data.frame(G) |> rownames_to_column("IID"), file.path(TAB_GRM, "natural_grm_vanraden.tsv"))
    write_tsv(as.data.frame(K) |> rownames_to_column("IID"), file.path(TAB_GRM, "natural_kinship_vanraden_half.tsv"))

    annot <- tibble(IID = sample_ids, Population = as.character(grp))
    saveRDS(annot, file.path(RDATA_DIR, "natural_sample_annotation.rds"))
    write_tsv(annot, file.path(TAB_GRM, "natural_sample_annotation.tsv"))

    message("NATURAL GRM dims: ", nrow(G), " x ", ncol(G))
  })
}

############################################################
# 5) IBD (PLINK .genome from step 7a)
# The PLINK MoM estimator (PI_HAT) was computed on LD-pruned SNPs
# (step 7a) and is loaded here for relatedness QC and plotting.
# Using the GRM sample order as the canonical ID set ensures the
# square matrix has consistent dimensions with downstream objects.
############################################################
if (RUN_LOAD_IBD) {
  message("=== IBD: load PLINK .genome (step 7a) ===")
  ensure_file(IBD_BREED)
  ensure_file(IBD_NATUR)

  # Derive the canonical sample order from the already-saved GRM
  breed_ids <- readRDS(file.path(RDATA_DIR, "breeding_grm_vanraden.rds")) |> rownames()
  ibd_b <- read_plink_genome(IBD_BREED, ids = breed_ids)
  saveRDS(ibd_b$long, file.path(RDATA_DIR, "breeding_ibd_plink_long.rds"))
  saveRDS(ibd_b$mat,  file.path(RDATA_DIR, "breeding_ibd_plink_mat.rds"))
  write_tsv(ibd_b$long, file.path(TAB_IBD, "breeding_ibd_plink_long.tsv"))

  natur_ids <- readRDS(file.path(RDATA_DIR, "natural_grm_vanraden.rds")) |> rownames()
  ibd_n <- read_plink_genome(IBD_NATUR, ids = natur_ids)
  saveRDS(ibd_n$long, file.path(RDATA_DIR, "natural_ibd_plink_long.rds"))
  saveRDS(ibd_n$mat,  file.path(RDATA_DIR, "natural_ibd_plink_mat.rds"))
  write_tsv(ibd_n$long, file.path(TAB_IBD, "natural_ibd_plink_long.tsv"))

  message("IBD imported: breeding pairs=", nrow(ibd_b$long), " natural pairs=", nrow(ibd_n$long))
}

############################################################
# 6) PCA (PLINK pruned PCA from step 5a)
# Eigenvectors computed by PLINK on LD-pruned SNPs are imported
# here rather than recomputed, to keep PCA consistent with the
# ADMIXTURE analysis which used the same pruned SNP set.
############################################################
if (RUN_IMPORT_PCA) {
  message("=== PCA: import PLINK eigenvec/eigenval (step 5a) ===")

  p_b <- read_plink_pca(PCA_BREED_EIGENVEC, PCA_BREED_EIGENVAL)
  saveRDS(p_b, file.path(RDATA_DIR, "breeding_pca_pruned.rds"))
  write_tsv(p_b$scores, file.path(TAB_PCA, "breeding_pca_pruned_scores.tsv"))
  # Store eigenvalues with their PC index so variance-explained plots can be
  # computed without reloading the full GDS.
  write_tsv(tibble(PC = seq_along(p_b$eigenval), eigenval = p_b$eigenval),
            file.path(TAB_PCA, "breeding_pca_pruned_eigenval.tsv"))

  p_n <- read_plink_pca(PCA_NATUR_EIGENVEC, PCA_NATUR_EIGENVAL)
  saveRDS(p_n, file.path(RDATA_DIR, "natural_pca_pruned.rds"))
  write_tsv(p_n$scores, file.path(TAB_PCA, "natural_pca_pruned_scores.tsv"))
  write_tsv(tibble(PC = seq_along(p_n$eigenval), eigenval = p_n$eigenval),
            file.path(TAB_PCA, "natural_pca_pruned_eigenval.tsv"))

  message("PCA imported.")
}

############################################################
# 7) ADMIXTURE (step 5a outputs)
# Only file paths and CV error summaries are stored here;
# the Q matrices themselves are large and read on demand by
# plotting scripts via the paths in the saved tibbles.
############################################################
if (RUN_IMPORT_ADMIXTURE) {
  message("=== ADMIXTURE: import Q matrices + CV summaries (step 5a) ===")

  cv_b_file <- file.path(ADMIX_BREED_DIR, "breeding_cv_errors.txt")
  cv_n_file <- file.path(ADMIX_NATUR_DIR, "natural_cv_errors.txt")
  # CV errors are generated by ADMIXTURE's --cv flag and are used to
  # select the optimal number of clusters K (minimum CV = best K).
  if (file.exists(cv_b_file)) {
    cv_b <- read_cv_errors(cv_b_file)
    saveRDS(cv_b, file.path(RDATA_DIR, "breeding_admixture_cv.rds"))
    write_tsv(cv_b, file.path(TAB_ADMIX, "breeding_admixture_cv.tsv"))
  } else {
    message("NOTE: missing breeding CV file: ", cv_b_file)
  }
  if (file.exists(cv_n_file)) {
    cv_n <- read_cv_errors(cv_n_file)
    saveRDS(cv_n, file.path(RDATA_DIR, "natural_admixture_cv.rds"))
    write_tsv(cv_n, file.path(TAB_ADMIX, "natural_admixture_cv.tsv"))
  } else {
    message("NOTE: missing natural CV file: ", cv_n_file)
  }

  # Prefix used when naming ADMIXTURE output files (must match step 5a convention)
  breed_prefix <- "tgc.ecs.breeding.admix.pruned"
  natur_prefix <- "tgc.ecs.natural.admix.pruned"

  # Collect paths and K values for all Q files; downstream plotting scripts
  # iterate over this tibble to build structure bar charts for each K.
  q_b <- collect_q_files(file.path(STRUCT_DIR, "BREEDING"), ADMIX_BREED_DIR, breed_prefix, max_k = 30)
  q_n <- collect_q_files(file.path(STRUCT_DIR, "NATURAL"),  ADMIX_NATUR_DIR, natur_prefix, max_k = 30)

  saveRDS(q_b, file.path(RDATA_DIR, "breeding_admixture_q_files.rds"))
  saveRDS(q_n, file.path(RDATA_DIR, "natural_admixture_q_files.rds"))
  write_tsv(q_b, file.path(TAB_ADMIX, "breeding_admixture_q_files.tsv"))
  write_tsv(q_n, file.path(TAB_ADMIX, "natural_admixture_q_files.tsv"))

  message("ADMIXTURE imported (file indices + CV if available).")
}

############################################################
# 8) DAPC INPUT (QC-filtered imputed dosages; no DAPC run)
#    UPDATED: now also saves SNP map (loc/chr/pos) aligned to X
#
# The DAPC dosage matrix (X) uses a relaxed MAF/missingness filter
# relative to the GRM because DAPC is a descriptive clustering method
# rather than a linear mixed model requiring a positive-definite matrix.
# Invariant SNPs (sd == 0) are still removed because they carry no
# discriminant information and can destabilise the DA step.
############################################################
if (RUN_SAVE_DAPC_INPUT) {
  message("=== DAPC INPUT: save QC-filtered imputed dosage matrices ===")

  # -------------------- BREEDING --------------------
  gds_breed <- file.path(RDATA_DIR, "breeding.imputed.snp.gds")
  ensure_file(gds_breed)
  breed_map <- if (file.exists(BREED_MAP_FILE)) read_map2(BREED_MAP_FILE, "Family") else tibble(IID = character(), Family = character())

  local({
    gf <- snpgdsOpen(gds_breed, allow.duplicate = FALSE)
    on.exit(try(snpgdsClose(gf), silent = TRUE), add = TRUE)

    samp <- read.gdsn(index.gdsn(gf, "sample.id"))

    # Read the full SNP annotation vectors once; used later to build the
    # aligned SNP map after filtering, avoiding multiple GDS traversals.
    # cache full SNP vectors once (minimal overhead; avoids ambiguous mapping)
    snp_id_all <- read.gdsn(index.gdsn(gf, "snp.id"))
    chr_all    <- read.gdsn(index.gdsn(gf, "snp.chromosome"))
    pos_all    <- read.gdsn(index.gdsn(gf, "snp.position"))

    # Compute per-SNP MAF and missing-call rate across all samples in one pass
    stat <- snpgdsSNPRateFreq(gf, with.id = TRUE, sample.id = samp)
    maf  <- stat$MinorFreq
    miss <- stat$MissingRate
    snp_id <- stat$snp.id

    keep <- is.finite(maf) & is.finite(miss) & maf >= DAPC_MAF_MIN & miss <= DAPC_MISS_MAX
    if (sum(keep) == 0) stop("BREEDING DAPC: 0 SNPs pass QC (maf/miss thresholds).")

    snp_id_keep <- snp_id[keep]

    geno <- snpgdsGetGeno(gf, sample.id = samp, snp.id = snp_id_keep, snpfirstdim = TRUE, with.id = FALSE)
    mode(geno) <- "numeric"
    geno[geno > 2 | geno < 0] <- NA_real_

    grp <- rep("Unknown", length(samp))
    if (nrow(breed_map) > 0) {
      g2 <- breed_map$Family[match(samp, breed_map$IID)]
      g2[is.na(g2)] <- "Unknown"
      grp <- g2
    }
    grp <- factor(grp)

    # Remove SNPs that remain with undefined allele frequency after MAF filter
    p_all <- rowMeans(geno, na.rm = TRUE) / 2
    keep2 <- is.finite(p_all)
    geno2 <- geno[keep2, , drop = FALSE]
    p_all <- p_all[keep2]

    snp_id_keep2 <- snp_id_keep[keep2]  # track SNP IDs through each filtering step

    # Within-family mean imputation (same rationale as GRM block in Section 4)
    for (g in levels(grp)) {
      idx <- which(grp == g)
      if (!length(idx)) next
      pg <- rowMeans(geno2[, idx, drop = FALSE], na.rm = TRUE) / 2
      pg[!is.finite(pg)] <- p_all[!is.finite(pg)]
      missm <- is.na(geno2[, idx, drop = FALSE])
      if (any(missm)) {
        fill <- matrix(2 * pg, nrow = nrow(geno2), ncol = length(idx))
        geno2[, idx][missm] <- fill[missm]
      }
    }
    if (anyNA(geno2)) {
      fill2 <- matrix(2 * p_all, nrow = nrow(geno2), ncol = ncol(geno2))
      geno2[is.na(geno2)] <- fill2[is.na(geno2)]
    }
    geno2[geno2 < 0] <- 0; geno2[geno2 > 2] <- 2

    X <- t(geno2); rownames(X) <- samp
    sdv <- apply(X, 2, sd)
    keep_var <- is.finite(sdv) & sdv > 0
    if (!all(keep_var)) X <- X[, keep_var, drop = FALSE]

    snp_id_final <- snp_id_keep2[keep_var]  # final post-QC SNP set

    # Build a SNP annotation table (locus ID, chromosome, position) aligned
    # column-for-column with X, so downstream DAPC loadings can be mapped
    # back to genomic coordinates.
    # build snp map aligned to X columns
    idx_map <- match(snp_id_final, snp_id_all)
    snp_df <- tibble(
      loc = as.character(snp_id_final),
      chr = as.character(chr_all[idx_map]),
      pos = as.integer(pos_all[idx_map])
    )

    # ensure loadings rownames == loc
    colnames(X) <- snp_df$loc

    # Bundle all components needed by downstream DAPC scripts into one RDS,
    # including QC parameters for reproducibility documentation.
    out <- list(
      X = X,
      snp = snp_df,
      group = tibble(IID = samp, Group = as.character(grp)),
      qc = list(maf_min = DAPC_MAF_MIN, miss_max = DAPC_MISS_MAX),
      note = "Imputed dosages derived from GT in Beagle-imputed VCF via SNPRelate GDS."
    )
    out_file <- file.path(RDATA_DIR, sprintf("breeding_dapc_input_maf%.2f_miss%.2f.rds", DAPC_MAF_MIN, DAPC_MISS_MAX))
    saveRDS(out, out_file)
    write_tsv(out$group, file.path(TAB_DAPC, "breeding_groups.tsv"))
    write_tsv(out$snp,   file.path(TAB_DAPC, "breeding_snps_loc_chr_pos.tsv"))
    message("Saved: ", out_file, " [", nrow(X), " x ", ncol(X), "]")
  })

  # -------------------- NATURAL --------------------
  gds_natur <- file.path(RDATA_DIR, "natural.imputed.snp.gds")
  ensure_file(gds_natur)
  natur_map <- if (file.exists(NATUR_MAP_FILE)) read_map2(NATUR_MAP_FILE, "Population") else tibble(IID = character(), Population = character())

  # Identical workflow to BREEDING above
  local({
    gf <- snpgdsOpen(gds_natur, allow.duplicate = FALSE)
    on.exit(try(snpgdsClose(gf), silent = TRUE), add = TRUE)

    samp <- read.gdsn(index.gdsn(gf, "sample.id"))

    snp_id_all <- read.gdsn(index.gdsn(gf, "snp.id"))
    chr_all    <- read.gdsn(index.gdsn(gf, "snp.chromosome"))
    pos_all    <- read.gdsn(index.gdsn(gf, "snp.position"))

    stat <- snpgdsSNPRateFreq(gf, with.id = TRUE, sample.id = samp)
    maf  <- stat$MinorFreq
    miss <- stat$MissingRate
    snp_id <- stat$snp.id

    keep <- is.finite(maf) & is.finite(miss) & maf >= DAPC_MAF_MIN & miss <= DAPC_MISS_MAX
    if (sum(keep) == 0) stop("NATURAL DAPC: 0 SNPs pass QC (maf/miss thresholds).")

    snp_id_keep <- snp_id[keep]

    geno <- snpgdsGetGeno(gf, sample.id = samp, snp.id = snp_id_keep, snpfirstdim = TRUE, with.id = FALSE)
    mode(geno) <- "numeric"
    geno[geno > 2 | geno < 0] <- NA_real_

    grp <- rep("Unknown", length(samp))
    if (nrow(natur_map) > 0) {
      g2 <- natur_map$Population[match(samp, natur_map$IID)]
      g2[is.na(g2)] <- "Unknown"
      grp <- g2
    }
    grp <- factor(grp)

    p_all <- rowMeans(geno, na.rm = TRUE) / 2
    keep2 <- is.finite(p_all)
    geno2 <- geno[keep2, , drop = FALSE]
    p_all <- p_all[keep2]

    snp_id_keep2 <- snp_id_keep[keep2]

    for (g in levels(grp)) {
      idx <- which(grp == g)
      if (!length(idx)) next
      pg <- rowMeans(geno2[, idx, drop = FALSE], na.rm = TRUE) / 2
      pg[!is.finite(pg)] <- p_all[!is.finite(pg)]
      missm <- is.na(geno2[, idx, drop = FALSE])
      if (any(missm)) {
        fill <- matrix(2 * pg, nrow = nrow(geno2), ncol = length(idx))
        geno2[, idx][missm] <- fill[missm]
      }
    }
    if (anyNA(geno2)) {
      fill2 <- matrix(2 * p_all, nrow = nrow(geno2), ncol = ncol(geno2))
      geno2[is.na(geno2)] <- fill2[is.na(geno2)]
    }
    geno2[geno2 < 0] <- 0; geno2[geno2 > 2] <- 2

    X <- t(geno2); rownames(X) <- samp
    sdv <- apply(X, 2, sd)
    keep_var <- is.finite(sdv) & sdv > 0
    if (!all(keep_var)) X <- X[, keep_var, drop = FALSE]

    snp_id_final <- snp_id_keep2[keep_var]

    idx_map <- match(snp_id_final, snp_id_all)
    snp_df <- tibble(
      loc = as.character(snp_id_final),
      chr = as.character(chr_all[idx_map]),
      pos = as.integer(pos_all[idx_map])
    )

    colnames(X) <- snp_df$loc

    out <- list(
      X = X,
      snp = snp_df,
      group = tibble(IID = samp, Group = as.character(grp)),
      qc = list(maf_min = DAPC_MAF_MIN, miss_max = DAPC_MISS_MAX),
      note = "Imputed dosages derived from GT in Beagle-imputed VCF via SNPRelate GDS."
    )
    out_file <- file.path(RDATA_DIR, sprintf("natural_dapc_input_maf%.2f_miss%.2f.rds", DAPC_MAF_MIN, DAPC_MISS_MAX))
    saveRDS(out, out_file)
    write_tsv(out$group, file.path(TAB_DAPC, "natural_groups.tsv"))
    write_tsv(out$snp,   file.path(TAB_DAPC, "natural_snps_loc_chr_pos.tsv"))
    message("Saved: ", out_file, " [", nrow(X), " x ", ncol(X), "]")
  })
}

message("=== 8a DONE. Outputs under: ", RANA_DIR)
sessionInfo()
