#!/bin/bash
#-------------------------------------------------------------------------------
# TreeGeneClimate (TGC) — ECS
# Step 1a: FASTQC + MULTIQC on RAWDATA (ECS)
#
# Project root:
#   /path/to/your/project
#
# Input (data, frozen):
#   DATA/ECS/RAWDATA.ECS/*fastq.gz
#
# Output (results):
#   RESULTS/ECS/QC/RAWDATA/FASTQC/
#   RESULTS/ECS/QC/RAWDATA/MULTIQC/
#-------------------------------------------------------------------------------

#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH -n 48
#SBATCH -N 1
#SBATCH --job-name=ECS.FQC1
#SBATCH --output=/path/to/your/project/LOGS/%x_%j.out
#SBATCH --error=/path/to/your/project/LOGS/%x_%j.err
#SBATCH --ntasks-per-socket 24
#SBATCH --time=48:00:00
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=YOUR_EMAIL

set -euo pipefail

echo "[$(date)] SLURM job started: ${SLURM_JOB_NAME:-no_slurm}"

module purge
module load fastqc/0.11.4
module load multiqc/1.27.1

# === USER CONFIGURATION ===
PROJECT_ROOT="/path/to/your/project"  # <-- set this
# ===========================

INPUT="${PROJECT_ROOT}/DATA/ECS/RAWDATA.ECS"
QC_BASE="${PROJECT_ROOT}/RESULTS/ECS/QC/RAWDATA"
QC_FASTQC="${QC_BASE}/FASTQC"
QC_MULTIQC="${QC_BASE}/MULTIQC"
LOGS="${PROJECT_ROOT}/LOGS"

mkdir -p "${QC_FASTQC}" "${QC_MULTIQC}" "${LOGS}"

shopt -s nullglob
raw_fastq=( "${INPUT}"/*.fastq.gz )
if (( ${#raw_fastq[@]} == 0 )); then
  echo "ERROR: no FASTQ found in: ${INPUT}"
  exit 1
fi

fastqc "${raw_fastq[@]}" --outdir "${QC_FASTQC}" --threads 48

source activate multiqc
multiqc "${QC_FASTQC}" -o "${QC_MULTIQC}"
conda deactivate

echo "[$(date)] Done."
exit 0
