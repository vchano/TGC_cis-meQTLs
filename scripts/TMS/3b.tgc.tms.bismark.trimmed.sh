#!/bin/bash
#-------------------------------------------------------------------------------
# TreeGeneClimate (TGC) — TMS
# Step 3b: BISMARK mapping on TRIMMED FASTQ (TMS), plates P001–P008
#
# Project root:
#   /path/to/your/project
#
# Input (data, frozen):
#   DATA/TMS/TRIMMED.FASTQ.TMS/*_R1_p.fastq.gz
#   DATA/TMS/TRIMMED.FASTQ.TMS/*_R2_p.fastq.gz
#
# Reference (frozen):
#   REFERENCE/Picab02_genome.assembly/bismark/
#
# Output (data, frozen):
#   DATA/TMS/MAPPED.FILES.TMS/<PLATE>/
#-------------------------------------------------------------------------------

#SBATCH -n 96
#SBATCH -N 1
#SBATCH --job-name=TMS.BISMARK
#SBATCH --output=LOGS/%x_%j.out
#SBATCH --error=LOGS/%x_%j.err
#SBATCH --exclusive
#SBATCH --time=48:00:00

set -euo pipefail

echo "[$(date)] SLURM job started: ${SLURM_JOB_NAME:-no_slurm}"

module purge
module load miniforge3/24.3.0-0
module load gcc/14.2.0
module load bowtie2/2.5.4
module load samtools/1.21
source activate bismark

# === USER CONFIGURATION ===
PROJECT_ROOT="${TGC_PROJECT_ROOT:-/path/to/your/project}"
# ===========================

INPUT="${PROJECT_ROOT}/DATA/TMS/TRIMMED.FASTQ.TMS"
OUTBASE="${PROJECT_ROOT}/DATA/TMS/MAPPED.FILES.TMS"
REF="${PROJECT_ROOT}/REFERENCE/Picab02_genome.assembly/bismark"

LOGS="${PROJECT_ROOT}/LOGS"
mkdir -p "${OUTBASE}" "${LOGS}"

list=$(cat "${PROJECT_ROOT}/DATA/METADATA/samples_tms_sequenced_620.txt")

for sample in ${list}; do
  plate="${sample%%_*}"
  outdir="${OUTBASE}/${plate}"

  r1="${INPUT}/${sample}_R1_p.fastq.gz"
  r2="${INPUT}/${sample}_R2_p.fastq.gz"

  if [[ ! -s "${r1}" || ! -s "${r2}" ]]; then
    echo "WARNING: missing pair for ${sample}; skipping"
    continue
  fi

  mkdir -p "${outdir}"

  bismark -q -L 15 -N 1 --score_min L,-1,-1 \
    --parallel 4 --un --ambiguous \
    -o "${outdir}" \
    --genome "${REF}" \
    -1 "${r1}" -2 "${r2}"
done

conda deactivate
echo "[$(date)] Done."
exit 0

