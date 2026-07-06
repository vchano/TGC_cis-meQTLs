#!/bin/bash
#-------------------------------------------------------------------------------
# TreeGeneClimate (TGC) — TBS
# Step 5b: methylKit import from sorted Bismark BAMs (by cohort + context)
#   - Ensure BAM indexes exist (CSI/BAI)
#   - Import methylation calls from BAMs (processBismarkAln)
#   - Build methylRawList (proper constructor + treatment)
#   - Filter by coverage (5 <= cov <= 50)
#   - Unite into methylBase (context-aware min.per.group; integer!)
#   - Save TWO RDS per condition:
#       1) methylBase after unite (cov-filtered only)
#       2) methylBase after MEF filter (cov + MEF)
#   - Write summary CSV with site counts at each stage
#
# INPUT:
#   DATA/TBS/MAPPED.FILES.TBS/*_R1_p_bismark_bt2_pe.sorted.bam
#
# OUTPUT:
#   RESULTS/TBS/RANALYSIS/METHYLKIT_OBJECTS/
#     methylBase_<cohort>_<context>_cov5_50_mpg<mpg>.rds
#     methylBase_<cohort>_<context>_cov5_50_mpg<mpg>_mef0.05.rds
#     summary_<cohort>_<context>_cov5_50_mpg<mpg>_mef0.05.csv
#
# JOB ARRAY MAP (6 tasks):
#   0: BREEDING / CpG
#   1: BREEDING / CHG
#   2: BREEDING / CHH
#   3: NATURAL  / CpG
#   4: NATURAL  / CHG
#   5: NATURAL  / CHH
#-------------------------------------------------------------------------------

#SBATCH -p YOUR_PARTITION
#SBATCH -t 02-00:00:00
#SBATCH -N 1
#SBATCH -c 24
#SBATCH --mem=120G
#SBATCH --job-name=TBS.METHYLKIT
#SBATCH --output=/path/to/your/project/LOGS/%x_%A_%a.out
#SBATCH --error=/path/to/your/project/LOGS/%x_%A_%a.err
#SBATCH --array=0-5%2
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=YOUR_EMAIL

set -euo pipefail

# -----------------------------
# Modules / environment
# -----------------------------
module purge
module load gcc/14.2.0
module load samtools/1.21
module load r/4.5.2

export LC_ALL=C.UTF-8
export LANG=C.UTF-8

# IMPORTANT: point to the library you want SLURM jobs to use
export R_LIBS_USER="/mnt/vast-standard/home/chano/u15584/Rlibs/4.5.2"

# -----------------------------
# Config
# -----------------------------
# === USER CONFIGURATION ===
PROJECT_ROOT="/path/to/your/project"  # <-- set this
# ===========================
BAM_DIR="${PROJECT_ROOT}/DATA/TBS/MAPPED.FILES.TBS"
OUT_BASE="${PROJECT_ROOT}/RESULTS/TBS/RANALYSIS/METHYLKIT_OBJECTS"
TMPDIR="${OUT_BASE}/_tmp_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"

mkdir -p "${OUT_BASE}" "${TMPDIR}" "${R_LIBS_USER}"

# -----------------------------
# Map array task -> cohort/context
# -----------------------------
contexts=(CpG CHG CHH)
cohorts=(BREEDING NATURAL)

task="${SLURM_ARRAY_TASK_ID}"
cohort="${cohorts[$(( task / 3 ))]}"
ctx="${contexts[$(( task % 3 ))]}"

echo "============================================================"
echo "TGC — TBS — Step 5b"
echo "Task:   ${task}"
echo "Cohort: ${cohort}"
echo "Ctx:    ${ctx}"
echo "Node:   $(hostname)"
echo "Start:  $(date)"
echo "R:      $(which R)"
echo "R_LIBS_USER: ${R_LIBS_USER}"
echo "============================================================"

# -----------------------------
# Select BAMs by cohort
# -----------------------------
shopt -s nullglob
if [[ "${cohort}" == "BREEDING" ]]; then
  bams=( "${BAM_DIR}"/P00[1-3]_*_R1_p_bismark_bt2_pe.sorted.bam )
else
  bams=( "${BAM_DIR}"/P00[4-8]_*_R1_p_bismark_bt2_pe.sorted.bam )
fi

if [[ ${#bams[@]} -eq 0 ]]; then
  echo "ERROR: No BAMs found for ${cohort} under: ${BAM_DIR}"
  exit 1
fi

echo "Found BAMs: ${#bams[@]}"

# -----------------------------
# Ensure BAM indexes exist (CSI or BAI)
# -----------------------------
echo "Checking/creating BAM indexes (.csi OR .bai) ..."
missing_idx=()
for b in "${bams[@]}"; do
  if [[ ! -f "${b}.csi" && ! -f "${b%.bam}.csi" && ! -f "${b}.bai" && ! -f "${b%.bam}.bai" ]]; then
    missing_idx+=( "$b" )
  fi
done

if [[ ${#missing_idx[@]} -gt 0 ]]; then
  echo "Missing indexes: ${#missing_idx[@]} (creating CSI with samtools index -c)"
  printf "%s\n" "${missing_idx[@]}" \
    | xargs -n 1 -P 6 -I {} samtools index -c -@ 2 "{}"
else
  echo "All indexes present."
fi

# -----------------------------
# Install/check methylKit once (avoid array races)
# -----------------------------
PKG_LOCK="${OUT_BASE}/.pkg_install_lock"
(
  flock -n 9 || exit 0
  Rscript --vanilla - <<'RS'
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
if (!requireNamespace("methylKit", quietly = TRUE)) {
  message("methylKit not found -> installing via BiocManager...")
  install.packages("BiocManager", repos="https://cloud.r-project.org")
  BiocManager::install("methylKit", ask=FALSE, update=FALSE)
}
suppressPackageStartupMessages(library(methylKit))
message("methylKit OK: ", as.character(packageVersion("methylKit")))
RS
) 9>"${PKG_LOCK}"

# -----------------------------
# Write and run R script
# -----------------------------
R_SCRIPT="${TMPDIR}/step5b_${cohort}_${ctx}.R"

cat > "${R_SCRIPT}" << 'RSCRIPT'
#!/usr/bin/env Rscript
############################################################
# TGC — TBS — Step 5b
# Save TWO methylBase RDS per cohort/context:
#   1) after cov + unite
#   2) after cov + unite + MEF
############################################################

.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressPackageStartupMessages(library(methylKit))

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT")
BAM_DIR      <- Sys.getenv("BAM_DIR")
OUT_BASE     <- Sys.getenv("OUT_BASE")
COHORT       <- Sys.getenv("COHORT")
CTX          <- Sys.getenv("CTX")

if (any(c(PROJECT_ROOT, BAM_DIR, OUT_BASE, COHORT, CTX) == "")) {
  stop("Missing env vars: PROJECT_ROOT/BAM_DIR/OUT_BASE/COHORT/CTX")
}

assembly_label <- "Picab02" # label only
mincov <- 5L
maxcov <- 50L
mef_threshold <- 0.05

# Context-aware unite threshold (fraction of samples required at a site)
mpg_frac <- switch(
  CTX,
  "CpG" = 0.80,
  "CHG" = 0.70,
  "CHH" = 0.50,
  0.70
)

# --- MEF filter on methylBase
mef_filter_methylBase <- function(mb, thr = 0.05) {
  df <- getData(mb)
  numCs_cols <- grep("^numCs[0-9]+$", colnames(df), value = TRUE)
  numTs_cols <- grep("^numTs[0-9]+$", colnames(df), value = TRUE)
  if (length(numCs_cols) == 0 || length(numTs_cols) == 0 || length(numCs_cols) != length(numTs_cols)) {
    stop("Could not find matching numCs/numTs columns in methylBase.")
  }
  totalCs <- rowSums(df[, numCs_cols, drop = FALSE], na.rm = TRUE)
  totalTs <- rowSums(df[, numTs_cols, drop = FALSE], na.rm = TRUE)
  total   <- totalCs + totalTs
  mef <- ifelse(total > 0, totalCs / total, NA_real_)
  keep <- !is.na(mef) & (mef > thr)
  list(filtered = mb[keep, ], n_before = nrow(mb), n_after = sum(keep))
}

# --- BAM discovery for cohort
pattern <- if (COHORT == "BREEDING") "^P00[1-3]_" else "^P00[4-8]_"
bam_files <- list.files(
  BAM_DIR,
  pattern = paste0(pattern, ".*_R1_p_bismark_bt2_pe\\.sorted\\.bam$"),
  full.names = TRUE
)
if (length(bam_files) == 0) stop("No BAMs found for cohort: ", COHORT)

bam_files <- sort(bam_files)
sample_ids <- sub("_R1_p_bismark_bt2_pe\\.sorted\\.bam$", "", basename(bam_files))

cat("\n============================================================\n")
cat("Step 5b — ", COHORT, " / ", CTX, "\n", sep="")
cat("Samples: ", length(bam_files), "\n", sep="")
cat("============================================================\n")

# --- Import each BAM
mrl_list <- lapply(seq_along(bam_files), function(i) {
  processBismarkAln(
    location     = bam_files[i],
    sample.id    = sample_ids[i],
    assembly     = assembly_label,
    read.context = CTX
  )
})

# methylRawList (dummy treatment; grouping is used later in 6b)
treat <- rep(0L, length(mrl_list))
mrl <- methylRawList(mrl_list, treatment = treat)

# Counts: raw sites per sample
n_sites_raw <- vapply(mrl_list, nrow, integer(1))
raw_min    <- min(n_sites_raw)
raw_median <- as.integer(stats::median(n_sites_raw))
raw_max    <- max(n_sites_raw)
cat("Raw sites/sample (min/median/max): ", raw_min, "/", raw_median, "/", raw_max, "\n", sep="")

# --- Coverage filter per sample
mrl_cov <- filterByCoverage(mrl, lo.count = mincov, hi.count = maxcov)

# Counts: cov sites per sample
mrl_cov_list <- as.list(mrl_cov)
n_sites_cov <- vapply(mrl_cov_list, nrow, integer(1))
cov_min    <- min(n_sites_cov)
cov_median <- as.integer(stats::median(n_sites_cov))
cov_max    <- max(n_sites_cov)
cat("Cov-filter sites/sample (min/median/max): ", cov_min, "/", cov_median, "/", cov_max, "\n", sep="")

# --- Unite (must be integer)
n <- length(mrl_cov)
min_per_group <- as.integer(max(1L, ceiling(mpg_frac * n)))
cat("Using min.per.group = ", min_per_group, " (", mpg_frac, " of n=", n, ")\n", sep="")

mb_unite <- unite(mrl_cov, destrand = FALSE, min.per.group = min_per_group)
n_after_unite <- nrow(mb_unite)
cat("Sites after unite(): ", format(n_after_unite, big.mark=","), "\n", sep="")

if (n_after_unite == 0) {
  stop(
    "unite() returned 0 sites.\n",
    "Lower mpg_frac (especially for CHH) or relax coverage.\n",
    "Current: mincov=", mincov, " maxcov=", maxcov,
    " mpg_frac=", mpg_frac, " min.per.group=", min_per_group, " n=", n
  )
}

# --- Save RDS #1: after unite (cov only)
out_unite <- file.path(
  OUT_BASE,
  sprintf("methylBase_%s_%s_cov%d_%d_mpg%d.rds",
          tolower(COHORT), tolower(CTX), mincov, maxcov, min_per_group)
)
saveRDS(mb_unite, out_unite)
cat("Saved methylBase (cov+unite): ", out_unite, "\n", sep="")

# --- MEF filter
mf <- mef_filter_methylBase(mb_unite, mef_threshold)
mb_mef <- mf$filtered
n_after_mef <- nrow(mb_mef)
cat("Sites after MEF >", mef_threshold, ": ", format(n_after_mef, big.mark=","), "\n", sep="")

# --- Save RDS #2: after MEF
out_mef <- file.path(
  OUT_BASE,
  sprintf("methylBase_%s_%s_cov%d_%d_mpg%d_mef%.2f.rds",
          tolower(COHORT), tolower(CTX), mincov, maxcov, min_per_group, mef_threshold)
)
saveRDS(mb_mef, out_mef)
cat("Saved methylBase (cov+unite+MEF): ", out_mef, "\n", sep="")

# --- Summary CSV
summary_csv <- file.path(
  OUT_BASE,
  sprintf("summary_%s_%s_cov%d_%d_mpg%d_mef%.2f.csv",
          tolower(COHORT), tolower(CTX), mincov, maxcov, min_per_group, mef_threshold)
)

summary_tab <- data.frame(
  cohort = COHORT,
  context = CTX,
  n_samples = length(bam_files),

  raw_sites_min = raw_min,
  raw_sites_median = raw_median,
  raw_sites_max = raw_max,

  cov_sites_min = cov_min,
  cov_sites_median = cov_median,
  cov_sites_max = cov_max,

  mpg_fraction = mpg_frac,
  min_per_group_used = min_per_group,

  n_sites_after_unite = n_after_unite,
  n_sites_after_mef = n_after_mef,

  rds_cov_unite = out_unite,
  rds_cov_unite_mef = out_mef,
  stringsAsFactors = FALSE
)

write.csv(summary_tab, summary_csv, row.names = FALSE)
cat("Saved summary CSV: ", summary_csv, "\n", sep="")
RSCRIPT

chmod +x "${R_SCRIPT}"

export PROJECT_ROOT BAM_DIR OUT_BASE
export COHORT="${cohort}"
export CTX="${ctx}"

Rscript --vanilla "${R_SCRIPT}"

echo "============================================================"
echo "Done:  $(date)"
echo "Cohort/context finished: ${cohort} / ${ctx}"
echo "Outputs in: ${OUT_BASE}"
echo "============================================================"