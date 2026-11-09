#!/bin/bash
#SBATCH --job-name=TMS.5b2_5b3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --array=0-11
#SBATCH --output=LOGS/%x_%A_%a.out
#SBATCH --error=LOGS/%x_%A_%a.err
#-------------------------------------------------------------------------------
# TreeGeneClimate (TGC) — TMS
# Step 5b2_5b3: Steps 5b2 and 5b3 as one SLURM array (12 tasks):
#   tasks 0-5  : 5b2 population-level datasets  (cohort x context)
#   tasks 6-11 : 5b3 meQTL datasets             (cohort x context)
# Order within each block: BREEDING CpG, CHG, CHH, NATURAL CpG, CHG, CHH
#-------------------------------------------------------------------------------
set -euo pipefail

# === USER CONFIGURATION ===
export TGC_PROJECT_ROOT="${TGC_PROJECT_ROOT:-/path/to/your/project}"
module purge 2>/dev/null || true
module load gcc/14.2.0 r/4.5.2 2>/dev/null || true
export LC_ALL=C.UTF-8 LANG=C.UTF-8
# ===========================

i="${SLURM_ARRAY_TASK_ID:-0}"
cohorts=(BREEDING BREEDING BREEDING NATURAL NATURAL NATURAL)
contexts=(CpG CHG CHH CpG CHG CHH)
if [[ "${i}" -lt 6 ]]; then
  script="5b2.tgc.tms.methylkit.population.datasets.R"; k="${i}"
else
  script="5b3.tgc.tms.methylkit.meqtl.datasets.R"; k=$(( i - 6 ))
fi
Rscript --vanilla "${TGC_PROJECT_ROOT}/SCRIPTS/TMS/${script}" "${cohorts[$k]}" "${contexts[$k]}"
