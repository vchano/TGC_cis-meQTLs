#!/bin/bash
#SBATCH --job-name=ECS.7a.IBD
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=LOGS/%x_%j.out
#SBATCH --error=LOGS/%x_%j.err
#-------------------------------------------------------------------------------
# TreeGeneClimate (TGC) — ECS
# Step 7a: Pairwise identity by descent (PLINK --genome full) per cohort,
#          on the LD-pruned SNP sets
#
# Input (Step 5a):
#   RESULTS/ECS/POPGEN/STRUCTURE/<COHORT>/tgc.ecs.<cohort>.pruned.{bed,bim,fam}
#
# Output:
#   RESULTS/ECS/POPGEN/RELATEDNESS/IBD/tgc.ecs.<cohort>.pruned.ibd.genome
#
# Software: PLINK 1.9
#-------------------------------------------------------------------------------
set -euo pipefail

# === USER CONFIGURATION ===
PROJECT_ROOT="${TGC_PROJECT_ROOT:-/path/to/your/project}"
THREADS="${SLURM_CPUS_PER_TASK:-4}"
# Software environment (edit for your system)
module purge 2>/dev/null || true
module load gcc/14.2.0 plink/1.9 2>/dev/null || true
# ===========================

STRUCT_DIR="${PROJECT_ROOT}/RESULTS/ECS/POPGEN/STRUCTURE"
OUT_DIR="${PROJECT_ROOT}/RESULTS/ECS/POPGEN/RELATEDNESS/IBD"
mkdir -p "${OUT_DIR}"
echo "[$(date)] Step 7a started; $(plink --version 2>&1 | head -n 1)"

for cohort in breeding natural; do
  COH=$(echo "${cohort}" | tr '[:lower:]' '[:upper:]')
  bfile="${STRUCT_DIR}/${COH}/tgc.ecs.${cohort}.pruned"
  for ext in bed bim fam; do [[ -s "${bfile}.${ext}" ]] || { echo "ERROR: missing ${bfile}.${ext}"; exit 1; }; done
  out="${OUT_DIR}/tgc.ecs.${cohort}.pruned.ibd"
  plink --bfile "${bfile}" --allow-extra-chr --threads "${THREADS}" --genome full --out "${out}"
  echo "[$(date)] ${COH}: $(($(wc -l < "${out}.genome") - 1)) pairs"
done

echo "[$(date)] Step 7a done"
