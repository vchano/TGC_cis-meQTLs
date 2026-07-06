#!/bin/bash
#-------------------------------------------------------------------------------
# TreeGeneClimate (TGC) — JOINT ECS + TBS
# Step 15ab: meQTL overlap analysis — Venn diagrams + panel assembly
#
# Reads significant meQTL results (p_FDR < 1e-10) from 13ab.R and produces:
#   - 6 Venn diagrams (tool comparison per cohort×context; cohort vs context comparisons)
#   - SuppFig5 panel: 3×2 tool-comparison Venns
#   - Figure 6 panel: 1×5 contexts + cohorts Venns
#   - overlap_summary.tsv
#
# REQUIRES
#   R packages: data.table, ggvenn, ggplot2, patchwork
#   (installed automatically if absent)
#
# USAGE
#   sbatch 15ab.sh
#
# INPUTS
#   RESULTS/JOINT/COMBINED5/sig_sites/sig_p1e10_*.tsv  (from 13ab.R)
#
# OUTPUTS
#   RESULTS/JOINT/COMBINED5/overlap/
#-------------------------------------------------------------------------------

#SBATCH -p YOUR_PARTITION
#SBATCH -t 00:30:00
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH --job-name=TGC.overlap16
#SBATCH --output=/path/to/your/project/LOGS/%x_%j.out
#SBATCH --error=/path/to/your/project/LOGS/%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=YOUR_EMAIL

set -euo pipefail

# ---------------------------------------------------------------------------
# Environment
# ---------------------------------------------------------------------------
module purge
module load gcc/14.2.0
module load r/4.5.2

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
  cran_needed <- c("data.table", "ggvenn", "ggplot2", "patchwork")
  for (pkg in cran_needed) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message("Installing: ", pkg)
      install.packages(pkg, lib = lib, repos = "https://cloud.r-project.org")
    }
  }
'

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
# === USER CONFIGURATION ===
PROJECT_ROOT="/path/to/your/project"  # <-- set this
# ===========================
SCRIPTS="${PROJECT_ROOT}/SCRIPTS/JOINT"

mkdir -p "${PROJECT_ROOT}/LOGS"

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
echo "============================================================"
echo "TGC — JOINT — Step 15ab — meQTL overlap analysis"
echo "Node:   $(hostname)"
echo "Start:  $(date)"
echo "============================================================"

Rscript --vanilla "${SCRIPTS}/15ab.R"

echo "============================================================"
echo "Step 15ab finished: $(date)"
echo "============================================================"
