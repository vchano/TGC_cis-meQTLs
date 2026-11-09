#!/bin/bash
#-------------------------------------------------------------------------------
# TreeGeneClimate (TGC) — TMS
# Step 2b: TRIMMOMATIC (paired-end) + FASTQC + MULTIQC on TRIMMED FASTQ (TMS)
#
# Project root:
#   /path/to/your/project
#
# Input (data, frozen):
#   DATA/TMS/RAWDATA.TMS/HEL_101202_*_R1.fastq.gz
#   DATA/TMS/RAWDATA.TMS/HEL_101202_*_R2.fastq.gz
#
# Output (data, frozen):
#   DATA/TMS/TRIMMED.FASTQ.TMS/           (paired reads)
#   DATA/TMS/TRIMMED.FASTQ.TMS/UNPAIRED/  (unpaired reads)
#
# Output (results):
#   RESULTS/TMS/QC/TRIMMED/FASTQC/
#   RESULTS/TMS/QC/TRIMMED/MULTIQC/
#
# Adapters:
#   DATA/METADATA/adapters.fa
#-------------------------------------------------------------------------------

#SBATCH -n 48
#SBATCH -N 1
#SBATCH --job-name=TMS.TRIM
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

INPUT="${PROJECT_ROOT}/DATA/TMS/RAWDATA.TMS"
TRIMMED="${PROJECT_ROOT}/DATA/TMS/TRIMMED.FASTQ.TMS"
UNPAIRED="${TRIMMED}/UNPAIRED"

QC_BASE="${PROJECT_ROOT}/RESULTS/TMS/QC/TRIMMED"
QC_FASTQC="${QC_BASE}/FASTQC"
QC_MULTIQC="${QC_BASE}/MULTIQC"

ADAPTERS="${PROJECT_ROOT}/DATA/METADATA/adapters.fa"
LOGS="${PROJECT_ROOT}/LOGS"

mkdir -p "${TRIMMED}" "${UNPAIRED}" "${QC_FASTQC}" "${QC_MULTIQC}" "${LOGS}"

if [[ ! -s "${ADAPTERS}" ]]; then
  echo "ERROR: adapters file not found or empty: ${ADAPTERS}"
  exit 1
fi

list=$(cat "${PROJECT_ROOT}/DATA/METADATA/samples_tms_sequenced_620.txt")

for sample in ${list}; do
  r1="${INPUT}/HEL_101202_${sample}_R1.fastq.gz"
  r2="${INPUT}/HEL_101202_${sample}_R2.fastq.gz"

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

