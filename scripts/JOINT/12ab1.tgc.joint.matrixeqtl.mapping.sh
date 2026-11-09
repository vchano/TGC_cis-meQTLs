#!/bin/bash
#SBATCH --job-name=JOINT.12ab1.MEQTL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --array=0-5
#SBATCH --output=LOGS/%x_%A_%a.out
#SBATCH --error=LOGS/%x_%A_%a.err
#-------------------------------------------------------------------------------
# TreeGeneClimate (TGC) — JOINT
# Step 12ab1: MatrixEQTL cis-meQTL mapping, one array task per cohort x context
#-------------------------------------------------------------------------------
set -euo pipefail

# === USER CONFIGURATION ===
export TGC_PROJECT_ROOT="${TGC_PROJECT_ROOT:-/path/to/your/project}"
module purge 2>/dev/null || true
module load gcc/14.2.0 r/4.5.2 2>/dev/null || true
export LC_ALL=C.UTF-8 LANG=C.UTF-8
# ===========================

cohorts=(BREEDING BREEDING BREEDING NATURAL NATURAL NATURAL)
contexts=(CpG CHG CHH CpG CHG CHH)
i="${SLURM_ARRAY_TASK_ID:-0}"
export TGC_COHORT="${cohorts[$i]}" TGC_CONTEXT="${contexts[$i]}"
Rscript --vanilla "${TGC_PROJECT_ROOT}/SCRIPTS/JOINT/12ab1.tgc.joint.matrixeqtl.mapping.R"
