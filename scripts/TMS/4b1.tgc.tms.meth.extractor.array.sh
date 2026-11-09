#!/bin/bash
#-------------------------------------------------------------------------------
# TreeGeneClimate (TGC) — TMS
# Step 4b1: Bismark methylation extractor (CpG, CHG, CHH)
# Array-based, 1 sample per task
#
# Input:
#   DATA/TMS/MAPPED.FILES.TMS/*_sorted.bam
#
# Output:
#   RESULTS/TMS/METHYLATION_CALLS/<SAMPLE>/
#-------------------------------------------------------------------------------

#SBATCH -t 00-48:00:00
#SBATCH -N 1
#SBATCH -c 24
#SBATCH --mem=100G
#SBATC -C ssd
#SBATCH --job-name=TMS.MEX
#SBATCH --output=LOGS/%x_%A_%a.out
#SBATCH --error=LOGS/%x_%A_%a.err
#SBATCH --array=20-601%4

set -euo pipefail

module purge
module load gcc/14.2.0
module load bismark/0.25.1
module load bowtie2/2.5.4
module load samtools/1.21

# === USER CONFIGURATION ===
PROJECT_ROOT="${TGC_PROJECT_ROOT:-/path/to/your/project}"
# ===========================

INPUT="${PROJECT_ROOT}/DATA/TMS/MAPPED.FILES.TMS"
OUTPUT_BASE="${PROJECT_ROOT}/RESULTS/TMS/METHYLATION_CALLS"
REF="${PROJECT_ROOT}/REFERENCE/Picab02_genome.assembly/bismark"

mkdir -p "${OUTPUT_BASE}"

# -----------------------------
# Full 602-sample list
# -----------------------------
mapfile -t samples < "${PROJECT_ROOT}/DATA/METADATA/samples_final_602.txt"

SAMPLE="${samples[$SLURM_ARRAY_TASK_ID]}"

if [[ -z "${SAMPLE:-}" ]]; then
  echo "ERROR: No sample found for index ${SLURM_ARRAY_TASK_ID}"
  exit 1
fi

BAM="${INPUT}/${SAMPLE}_R1_p_bismark_bt2_pe.bam"
OUTDIR="${OUTPUT_BASE}/${SAMPLE}"

echo "============================================================"
echo "Sample: ${SAMPLE}"
echo "BAM: ${BAM}"
echo "Node: $(hostname)"
echo "Start: $(date)"
echo "============================================================"

if [[ ! -f "$BAM" ]]; then
  echo "ERROR: BAM not found: $BAM"
  exit 1
fi

mkdir -p "$OUTDIR"

bismark_methylation_extractor \
  --multicore 24 \
  --buffer_size 8G \
  --gzip \
  --bedGraph \
  --no_overlap \
  --CX_context \
  --genome_folder "$REF" \
  -o "$OUTDIR" \
  "$BAM"

echo "Finished: $(date)"
exit 0

