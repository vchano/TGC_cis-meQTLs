#!/bin/bash
#-------------------------------------------------------------------------------
# TreeGeneClimate (TGC) — JOINT ECS + TBS
# Step 16ab: Orchestrator — meQTL5 full pipeline (SLURM job array)
#
# Pipeline order (all within this single script):
#   1. Step 12ab0.R  — input prep (GRM-based PCs, LOCO GRM, unfiltered M-values)
#                      run once (task 0); other tasks wait for a flag file
#   2. Step 12ab1.R  — MatrixEQTL5 cis-meQTL mapping (one panel per task)
#   3. Step 12ab2.R  — GENESIS5 cis-meQTL mapping (same panel)
#   4. Step 13ab.R  — combined results, lambda table, QQ plots
#                      run once after all 6 tasks finish
#   5. Step 14ab.R  — circular Manhattan plots (both tools)
#                      run once after step 13ab completes
#   6. Step 16ab.R  — QQ + circular panel assembly (Figure 5, SuppFig 2-4)
#                      run once after step 14ab completes
#
# JOB ARRAY MAP (6 tasks: 2 cohorts x 3 contexts):
#   0: BREEDING / CpG    3: NATURAL / CpG
#   1: BREEDING / CHG    4: NATURAL / CHG
#   2: BREEDING / CHH    5: NATURAL / CHH
#
# KEY DIFFERENCES vs 15ab6 (meQTL4):
#   - Input prep uses non-imputed GDS, GCTA GRM, unfiltered methylation
#   - 10 PCs for BOTH cohorts
#   - MatrixEQTL5 uses 100 kb cis window
#   - GENESIS5 uses GRM for BOTH cohorts
#   - Outputs to MATRIXEQTL5, GENESIS5, COMBINED5
#
# USAGE
#   sbatch 16ab.sh
#
# OUTPUT roots
#   RESULTS/JOINT/MQTL5/          (prep inputs)
#   RESULTS/JOINT/MATRIXEQTL5/   (MatrixEQTL5 results)
#   RESULTS/JOINT/GENESIS5/       (GENESIS5 results)
#   RESULTS/JOINT/COMBINED5/      (combined tables, QQ plots, circular plots, panels)
#-------------------------------------------------------------------------------

#SBATCH -p YOUR_PARTITION
#SBATCH -t 04:00:00
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH --job-name=TGC.meQTL5
#SBATCH --output=/path/to/your/project/LOGS/%x_%A_%a.out
#SBATCH --error=/path/to/your/project/LOGS/%x_%A_%a.err
#SBATCH --array=0
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=YOUR_EMAIL

set -euo pipefail

# ---------------------------------------------------------------------------
# Environment
# ---------------------------------------------------------------------------
module purge
module load gcc/14.2.0
module load r/4.5.2
module load imagemagick/7.1.1-39

export LC_ALL=C.UTF-8
export LANG=C.UTF-8
export R_LIBS_USER="/mnt/vast-standard/home/chano/u15584/Rlibs/4.5.2"
mkdir -p "${R_LIBS_USER}"

# ---------------------------------------------------------------------------
# Install missing R packages (idempotent)
# ---------------------------------------------------------------------------
Rscript --vanilla -e '
  lib <- Sys.getenv("R_LIBS_USER")
  dir.create(lib, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(lib, .libPaths()))
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", lib = lib, repos = "https://cloud.r-project.org")
  bioc_needed <- c("SNPRelate", "gdsfmt", "GENESIS", "GWASTools", "Biobase", "methylKit")
  for (pkg in bioc_needed) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message("Installing: ", pkg)
      BiocManager::install(pkg, lib = lib, update = FALSE, ask = FALSE)
    }
  }
  cran_needed <- c("data.table", "MatrixEQTL", "circlize", "RColorBrewer", "magick")
  for (pkg in cran_needed) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message("Installing: ", pkg)
      install.packages(pkg, lib = lib, repos = "https://cloud.r-project.org")
    }
  }
'

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
# === USER CONFIGURATION ===
PROJECT_ROOT="/path/to/your/project"  # <-- set this
# ===========================
SCRIPTS="${PROJECT_ROOT}/SCRIPTS/JOINT"

PREP_SUMDIR="${PROJECT_ROOT}/RESULTS/JOINT/MQTL5/SUMMARIES"
ME5_SUMDIR="${PROJECT_ROOT}/RESULTS/JOINT/MATRIXEQTL5/SUMMARIES"
G5_SUMDIR="${PROJECT_ROOT}/RESULTS/JOINT/GENESIS5/SUMMARIES"

mkdir -p "${PREP_SUMDIR}" "${ME5_SUMDIR}" "${G5_SUMDIR}"
mkdir -p "${PROJECT_ROOT}/LOGS"

# ---------------------------------------------------------------------------
# Array task -> cohort / context
# ---------------------------------------------------------------------------
cohorts=(BREEDING NATURAL)
contexts=(CpG CHG CHH)

task="${SLURM_ARRAY_TASK_ID}"
cohort="${cohorts[$(( task / 3 ))]}"
ctx="${contexts[$(( task % 3 ))]}"
panel_tag="${cohort,,}_${ctx,,}"

echo "============================================================"
echo "TGC — JOINT — Step 16ab — meQTL5 pipeline"
echo "Task:    ${task}"
echo "Cohort:  ${cohort}"
echo "Context: ${ctx}"
echo "Node:    $(hostname)"
echo "CPUs:    ${SLURM_CPUS_ON_NODE:-4}"
echo "Memory:  ${SLURM_MEM_PER_NODE:-?} MB"
echo "Start:   $(date)"
echo "============================================================"

export TGC_PROJECT_ROOT="${PROJECT_ROOT}"

# ---------------------------------------------------------------------------
# STEP 1 — Input preparation (12ab0.R)
# Only the first task to acquire the lock runs prep; others wait.
# ---------------------------------------------------------------------------
PREP_FLAG="${PREP_SUMDIR}/step12ab0_done.flag"
PREP_LOCK="${PREP_SUMDIR}/.prep_lock"
mkdir -p "${PREP_SUMDIR}"

(
  flock 9
  if [[ ! -f "${PREP_FLAG}" ]]; then
    echo "Task ${task}: Running 12ab0.R (input prep)..."
    Rscript --vanilla "${SCRIPTS}/12ab0.R"
    touch "${PREP_FLAG}"
    echo "Task ${task}: 12ab0.R finished."
  else
    echo "Task ${task}: Prep already done (flag found)."
  fi
) 9>"${PREP_LOCK}"

# All tasks wait until prep flag exists
until [[ -f "${PREP_FLAG}" ]]; do
  echo "Task ${task}: Waiting for prep to complete..."
  sleep 15
done
echo "Task ${task}: Prep confirmed — proceeding to mapping."

# ---------------------------------------------------------------------------
# STEP 2 — MatrixEQTL5 mapping (12ab1.R)
# ---------------------------------------------------------------------------
echo "------------------------------------------------------------"
echo "Task ${task}: Running 12ab1.R (MatrixEQTL5) for ${cohort}/${ctx}"
echo "------------------------------------------------------------"
export TGC_COHORT="${cohort}"
export TGC_CONTEXT="${ctx}"

Rscript --vanilla "${SCRIPTS}/12ab1.R"
echo "Task ${task}: 12ab1.R done — $(date)"

# ---------------------------------------------------------------------------
# STEP 3 — GENESIS5 mapping (12ab2.R)
# ---------------------------------------------------------------------------
echo "------------------------------------------------------------"
echo "Task ${task}: Running 12ab2.R (GENESIS5) for ${cohort}/${ctx}"
echo "------------------------------------------------------------"

Rscript --vanilla "${SCRIPTS}/12ab2.R"
echo "Task ${task}: 12ab2.R done — $(date)"

unset TGC_COHORT
unset TGC_CONTEXT

# Write per-panel completion flags
touch "${ME5_SUMDIR}/.done_${panel_tag}"
touch "${G5_SUMDIR}/.done_${panel_tag}"

echo "Task ${task}: Mapping complete for ${cohort}/${ctx} — $(date)"

# ---------------------------------------------------------------------------
# STEP 4 & 5 — Post-processing (13ab.R + 14ab.R)
# The last task to finish all 6 panels runs post-processing.
# ---------------------------------------------------------------------------
MERGE_LOCK="${G5_SUMDIR}/.merge_lock"

(
  flock 9

  n_done=$(ls "${G5_SUMDIR}"/.done_* 2>/dev/null | wc -l)

  if [[ "${n_done}" -ge 6 ]]; then
    echo "============================================================"
    echo "All 6 panels done — running post-processing (task ${task})"
    echo "============================================================"

    echo "Running 13ab.R (combined results, lambda, QQ plots)..."
    Rscript --vanilla "${SCRIPTS}/13ab.R"
    echo "13ab.R done — $(date)"

    echo "Running 14ab.R (circular Manhattan plots)..."
    Rscript --vanilla "${SCRIPTS}/14ab.R"
    echo "14ab.R done — $(date)"

    echo "Running 16ab.R (QQ + circular panel assembly)..."
    Rscript --vanilla "${SCRIPTS}/16ab.R"
    echo "16ab.R done — $(date)"

    echo "Post-processing complete: $(date)"
  else
    echo "Task ${task}: ${n_done}/6 panels done — post-processing deferred."
  fi
) 9>"${MERGE_LOCK}"

echo "============================================================"
echo "Task ${task} finished: $(date)"
echo "============================================================"
