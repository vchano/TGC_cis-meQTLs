#!/bin/bash
#-------------------------------------------------------------------------------
# TreeGeneClimate (TGC) — ECS
# Step 2a: TRIMMOMATIC (paired-end) + FASTQC + MULTIQC on TRIMMED FASTQ (ECS)
#
# Project root:
#   /path/to/your/project
#
# Input (data, frozen):
#   DATA/ECS/RAWDATA.ECS/*_R1.fastq.gz
#   DATA/ECS/RAWDATA.ECS/*_R2.fastq.gz
#
# Output (data, frozen):
#   DATA/ECS/TRIMMED.FASTQ.ECS/           (paired reads)
#   DATA/ECS/TRIMMED.FASTQ.ECS/UNPAIRED/  (unpaired reads)
#
# Output (results):
#   RESULTS/ECS/QC/TRIMMED/FASTQC/
#   RESULTS/ECS/QC/TRIMMED/MULTIQC/
#
# Adapters:
#   DATA/METADATA/adapters.fa
#-------------------------------------------------------------------------------

#SBATCH -n 48
#SBATCH -N 1
#SBATCH --job-name=ECS.TRIM
#SBATCH --output=LOGS/%x_%j.out
#SBATCH --error=LOGS/%x_%j.err
#SBATCH --ntasks-per-socket 24
#SBATCH --time=48:00:00

set -euo pipefail

echo "[$(date)] SLURM job started: ${SLURM_JOB_NAME:-no_slurm}"

module purge
module load fastqc/0.11.4
module load anaconda3/2020.11
module load trimmomatic/0.36

# === USER CONFIGURATION ===
PROJECT_ROOT="${TGC_PROJECT_ROOT:-/path/to/your/project}"
# ===========================

INPUT="${PROJECT_ROOT}/DATA/ECS/RAWDATA.ECS"
TRIMMED="${PROJECT_ROOT}/DATA/ECS/TRIMMED.FASTQ.ECS"
UNPAIRED="${TRIMMED}/UNPAIRED"

QC_BASE="${PROJECT_ROOT}/RESULTS/ECS/QC/TRIMMED"
QC_FASTQC="${QC_BASE}/FASTQC"
QC_MULTIQC="${QC_BASE}/MULTIQC"

ADAPTERS="${PROJECT_ROOT}/DATA/METADATA/adapters.fa"
LOGS="${PROJECT_ROOT}/LOGS"

mkdir -p "${TRIMMED}" "${UNPAIRED}" "${QC_FASTQC}" "${QC_MULTIQC}" "${LOGS}"

if [[ ! -s "${ADAPTERS}" ]]; then
  echo "ERROR: adapters file not found or empty: ${ADAPTERS}"
  exit 1
fi

# ------------------------------------------------------------------------------
# Sample list (DATA/METADATA)
list=$(cat "${PROJECT_ROOT}/DATA/METADATA/samples_final_602.txt")

for sample in ${list}; do
  r1="${INPUT}/${sample}_R1.fastq.gz"
  r2="${INPUT}/${sample}_R2.fastq.gz"

  if [[ ! -s "${r1}" || ! -s "${r2}" ]]; then
    echo "WARNING: missing pair for ${sample}; skipping"
    continue
  fi

  java -jar "${TRIMMOMATIC_JAR:-/path/to/trimmomatic-0.36.jar}" PE \
    -threads 48 -phred33 \
    "${r1}" "${r2}" \
    "${TRIMMED}/${sample}_R1_p.fastq.gz" "${UNPAIRED}/${sample}_R1_u.fastq.gz" \
    "${TRIMMED}/${sample}_R2_p.fastq.gz" "${UNPAIRED}/${sample}_R2_u.fastq.gz" \
    ILLUMINACLIP:"${ADAPTERS}":2:30:10 \
    CROP:138 HEADCROP:12 SLIDINGWINDOW:5:20 MINLEN:30
done

shopt -s nullglob
trim_fastq=( "${TRIMMED}"/*_p.fastq.gz )
if (( ${#trim_fastq[@]} == 0 )); then
  echo "ERROR: no trimmed paired FASTQ produced in: ${TRIMMED}"
  exit 1
fi

fastqc "${trim_fastq[@]}" --outdir "${QC_FASTQC}" --threads 48

source activate multiqc
multiqc "${QC_FASTQC}" -o "${QC_MULTIQC}"
conda deactivate

echo "[$(date)] Done."
exit 0

