#!/bin/bash
#-------------------------------------------------------------------------------
# TreeGeneClimate (TGC) — ECS
# Step 7a: Pairwise relatedness (IBD/PI_HAT) per cohort using PLINK --genome
#          Inputs: LD-pruned BED sets created in Step 5a (admix.pruned)
#          Outputs: .genome files for downstream R (orchestrator becomes Step 8a)
#
# Project root:
#   /path/to/your/project
#
# Inputs (from Step 5a):
#   RESULTS/ECS/POPGEN/STRUCTURE/BREEDING/tgc.ecs.breeding.admix.pruned.{bed,bim,fam}
#   RESULTS/ECS/POPGEN/STRUCTURE/NATURAL/tgc.ecs.natural.admix.pruned.{bed,bim,fam}
#
# Outputs:
#   RESULTS/ECS/POPGEN/RELATEDNESS/IBD/
#     - tgc.ecs.breeding.pruned.ibd.genome
#     - tgc.ecs.natural.pruned.ibd.genome
#     - (plus .log/.nosex)
#-------------------------------------------------------------------------------

#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --job-name=ECS.IBD.PLINK
#SBATCH --output=/path/to/your/project/LOGS/%x_%j.out
#SBATCH --error=/path/to/your/project/LOGS/%x_%j.err
#SBATCH --time=24:00:00
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=YOUR_EMAIL

set -euo pipefail

ts() { date +"%a %b %d %H:%M:%S %Z %Y"; }
die() { echo "[$(ts)] ERROR: $*" >&2; exit 1; }
need_cmd() { command -v "$1" >/dev/null 2>&1 || die "Command not found in PATH: $1"; }
need_file() { [[ -s "$1" ]] || die "Missing/empty file: $1"; }

echo "================================================================================"
echo "JobID = ${SLURM_JOB_ID:-NA}"
echo "User = ${SLURM_JOB_USER:-$USER}, Account = ${SLURM_JOB_ACCOUNT:-NA}"
echo "Partition = ${SLURM_JOB_PARTITION:-NA}, Nodelist = ${SLURM_JOB_NODELIST:-NA}"
echo "================================================================================"
echo "[$(ts)] SLURM job started: ${SLURM_JOB_NAME:-no_slurm}"
echo "[$(ts)] Host: $(hostname)"
echo "[$(ts)] CPUs: ${SLURM_CPUS_PER_TASK:-1}"

# --- modules ---
module purge
module load gcc/14.2.0
module load plink/1.9

need_cmd plink

echo "[$(ts)] plink: $(plink --version 2>&1 | head -n 1)"

# --- paths ---
# === USER CONFIGURATION ===
PROJECT_ROOT="/path/to/your/project"  # <-- set this
# ===========================
LOGS="${PROJECT_ROOT}/LOGS"

STRUCT_DIR="${PROJECT_ROOT}/RESULTS/ECS/POPGEN/STRUCTURE"
BREED_BFILE="${STRUCT_DIR}/BREEDING/tgc.ecs.breeding.admix.pruned"
NATUR_BFILE="${STRUCT_DIR}/NATURAL/tgc.ecs.natural.admix.pruned"

OUT_DIR="${PROJECT_ROOT}/RESULTS/ECS/POPGEN/RELATEDNESS/IBD"
mkdir -p "${LOGS}" "${OUT_DIR}"

# Sanity checks (BED trio)
need_file "${BREED_BFILE}.bed"; need_file "${BREED_BFILE}.bim"; need_file "${BREED_BFILE}.fam"
need_file "${NATUR_BFILE}.bed"; need_file "${NATUR_BFILE}.bim"; need_file "${NATUR_BFILE}.fam"

THREADS="${SLURM_CPUS_PER_TASK:-1}"

run_ibd() {
  local cohort="$1"
  local bfile="$2"
  local outprefix="$3"

  echo "--------------------------------------------------------------------------------"
  echo "[$(ts)] COHORT: ${cohort}"
  echo "[$(ts)] Input bfile: ${bfile}"
  echo "[$(ts)] Output: ${outprefix}.genome"
  echo "--------------------------------------------------------------------------------"

  # --genome full computes PI_HAT + Z0/Z1/Z2 etc.
  # --allow-extra-chr because contig/chrom names are non-standard (PA_cUP..., etc.)
  plink \
    --bfile "${bfile}" \
    --allow-extra-chr \
    --threads "${THREADS}" \
    --genome full \
    --out "${outprefix}"

  need_file "${outprefix}.genome"

  # quick summary
  echo "[$(ts)] ${cohort} .genome rows: $(($(wc -l < "${outprefix}.genome") - 1))"
  echo "[$(ts)] ${cohort} PI_HAT quick peek:"
  awk 'NR==1{print;next} NR<=6{print}' "${outprefix}.genome" | column -t || true
}

run_ibd "BREEDING" "${BREED_BFILE}" "${OUT_DIR}/tgc.ecs.breeding.pruned.ibd"
run_ibd "NATURAL"  "${NATUR_BFILE}" "${OUT_DIR}/tgc.ecs.natural.pruned.ibd"

echo "[$(ts)] ALL DONE."
echo "**************************************************"
exit 0
