#!/bin/bash
#SBATCH --job-name=JOINT.12ab2.GENESIS
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --array=0-47
#SBATCH --output=LOGS/12ab2/%x_%A_%a.out
#SBATCH --error=LOGS/12ab2/%x_%A_%a.err
#-------------------------------------------------------------------------------
# TreeGeneClimate (TGC) — JOINT
# Step 12ab2: GENESIS cis-meQTL mapping
#   Mapping:  48 array tasks = 6 panels (cohort x context) x 8 chunks of
#             chromosomes/scaffolds:   sbatch 12ab2...sh
#   Combine:  one task per panel (FDR over the whole panel), after mapping:
#             sbatch --array=0-5 --export=ALL,TGC_CHUNK=combine 12ab2...sh
# Create LOGS/12ab2 before submitting.
#-------------------------------------------------------------------------------
set -euo pipefail

# === USER CONFIGURATION ===
export TGC_PROJECT_ROOT="${TGC_PROJECT_ROOT:-/path/to/your/project}"
NCHUNKS=8
module purge 2>/dev/null || true
module load gcc/14.2.0 r/4.5.2 2>/dev/null || true
export LC_ALL=C.UTF-8 LANG=C.UTF-8
# ===========================

cohorts=(BREEDING BREEDING BREEDING NATURAL NATURAL NATURAL)
contexts=(CpG CHG CHH CpG CHG CHH)
i="${SLURM_ARRAY_TASK_ID:-0}"
export TGC_NCHUNKS="${NCHUNKS}"
if [[ "${TGC_CHUNK:-}" == "combine" ]]; then
  p="${i}"
else
  p=$(( i / NCHUNKS )); export TGC_CHUNK=$(( i % NCHUNKS + 1 ))
fi
export TGC_COHORT="${cohorts[$p]}" TGC_CONTEXT="${contexts[$p]}"
echo "[$(date)] ${TGC_COHORT} ${TGC_CONTEXT} chunk ${TGC_CHUNK}/${NCHUNKS}"
Rscript --vanilla "${TGC_PROJECT_ROOT}/SCRIPTS/JOINT/12ab2.tgc.joint.genesis.mapping.R"
