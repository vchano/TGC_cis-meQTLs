#!/bin/bash
#-------------------------------------------------------------------------------
# TreeGeneClimate (TGC) — ECS
# Step 3a: BOWTIE2 mapping + SAMTOOLS processing on TRIMMED FASTQ (ECS)
#
# Project root:
#   /path/to/your/project
#
# Input (data, frozen):
#   DATA/ECS/TRIMMED.FASTQ.ECS/*_R1_p.fastq.gz
#   DATA/ECS/TRIMMED.FASTQ.ECS/*_R2_p.fastq.gz
#
# Reference (frozen):
#   REFERENCE/Pabies2.0/   (bowtie2 index prefix: Pabies2.0)
#
# Output (data, frozen):
#   DATA/ECS/MAPPED.FILES.ECS/
#-------------------------------------------------------------------------------

#SBATCH -n 96
#SBATCH -N 1
#SBATCH --job-name=ECS.BT2
#SBATCH --output=LOGS/%x_%j.out
#SBATCH --error=LOGS/%x_%j.err
#SBATCH --ntasks-per-socket 24
#SBATCH --time=48:00:00

set -euo pipefail

echo "[$(date)] SLURM job started: ${SLURM_JOB_NAME:-no_slurm}"

module purge
module load gcc/14.2.0
module load bowtie2/2.5.4
module load samtools/1.21

# === USER CONFIGURATION ===
PROJECT_ROOT="${TGC_PROJECT_ROOT:-/path/to/your/project}"
# ===========================

INPUT="${PROJECT_ROOT}/DATA/ECS/TRIMMED.FASTQ.ECS"
OUTPUT="${PROJECT_ROOT}/DATA/ECS/MAPPED.FILES.ECS"

REF_BASE="${PROJECT_ROOT}/REFERENCE/Pabies2.0"
REF_PREFIX="${REF_BASE}/Pabies2.0"

TMP="${OUTPUT}/tmp"
OTHER="${OUTPUT}/OTHER_FILES"
LOGS="${PROJECT_ROOT}/LOGS"

mkdir -p "${OUTPUT}" "${TMP}" "${OTHER}" "${LOGS}"

# ------------------------------------------------------------------------------
# Sample list (DATA/METADATA)
list=$(cat "${PROJECT_ROOT}/DATA/METADATA/samples_final_602.txt")

for sample in ${list}; do
  r1="${INPUT}/${sample}_R1_p.fastq.gz"
  r2="${INPUT}/${sample}_R2_p.fastq.gz"

  if [[ ! -s "${r1}" || ! -s "${r2}" ]]; then
    echo "WARNING: missing pair for ${sample}; skipping"
    continue
  fi

  bowtie2 -p 96 -N 1 -L 25 -i S,1,2.00 --local \
    --un-gz "${OTHER}/${sample}.un.fastq.gz" \
    --al-gz "${OTHER}/${sample}.al.fastq.gz" \
    --un-conc-gz "${OTHER}/${sample}.un.conc.fastq.gz" \
    --al-conc-gz "${OTHER}/${sample}.al.conc.fastq.gz" \
    -x "${REF_PREFIX}" \
    -1 "${r1}" -2 "${r2}" \
    -S "${OUTPUT}/${sample}_bt2.sam"

  samtools view -S -b "${OUTPUT}/${sample}_bt2.sam" > "${OUTPUT}/${sample}_bt2.bam"

  samtools sort -l 1 -@ 96 -n \
    -o "${OUTPUT}/${sample}_bt2.nsrt.bam" \
    -T "${TMP}/${sample}_tmp" \
    "${OUTPUT}/${sample}_bt2.bam"

  samtools fixmate -@ 96 -O bam,level=1 -m \
    "${OUTPUT}/${sample}_bt2.nsrt.bam" \
    "${OUTPUT}/${sample}_bt2.nsrt.fix.bam"

  samtools sort -l 1 -@ 96 \
    -o "${OUTPUT}/${sample}_bt2.fix.psrt.bam" \
    -T "${TMP}/${sample}_tmp" \
    "${OUTPUT}/${sample}_bt2.nsrt.fix.bam"

  samtools markdup -@ 96 -r -s -O bam,level=1 \
    "${OUTPUT}/${sample}_bt2.fix.psrt.bam" \
    "${OUTPUT}/${sample}_bt2.fix.psrt.dedup.bam"

  rm -f "${OUTPUT}/${sample}_bt2.sam"
done

echo "[$(date)] Done."
exit 0

