#!/bin/bash
#SBATCH --job-name=TMS.5b1.IMPORT
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --array=0-601%60
#SBATCH --output=LOGS/5b1/%x_%A_%a.out
#SBATCH --error=LOGS/5b1/%x_%A_%a.err
#-------------------------------------------------------------------------------
# TreeGeneClimate (TGC) — TMS
# Step 5b1: methylKit import, one array task per sample (602 samples)
#
# Input:  DATA/TMS/MAPPED.FILES.TMS/<SAMPLE>_R1_p_bismark_bt2_pe.dedup.bam  (Step 3b2)
# Output: RESULTS/TMS/METHYLKIT_RAW/<SAMPLE>_<CONTEXT>.rds
#
# Create LOGS/5b1 before submitting.
# Software: R 4.5.2, methylKit 1.36.0
#-------------------------------------------------------------------------------
set -euo pipefail

# === USER CONFIGURATION ===
PROJECT_ROOT="${TGC_PROJECT_ROOT:-/path/to/your/project}"
# Software environment (edit for your system)
module purge 2>/dev/null || true
module load gcc/14.2.0 r/4.5.2 2>/dev/null || true
# export R_LIBS_USER=/path/to/your/Rlibs
export LC_ALL=C.UTF-8 LANG=C.UTF-8
# ===========================

BAM_DIR="${PROJECT_ROOT}/DATA/TMS/MAPPED.FILES.TMS"
BAM_SUFFIX="${TGC_BAM_SUFFIX:-_R1_p_bismark_bt2_pe.dedup.bam}"
i="${SLURM_ARRAY_TASK_ID:-0}"
sample=$(sed -n "$((i + 1))p" "${PROJECT_ROOT}/DATA/METADATA/samples_final_602.txt")
[[ -n "${sample}" ]] || { echo "Task ${i}: no sample"; exit 0; }

export TGC_PROJECT_ROOT="${PROJECT_ROOT}"
Rscript --vanilla "${PROJECT_ROOT}/SCRIPTS/TMS/5b1.tgc.tms.methylkit.import.R" "${BAM_DIR}/${sample}${BAM_SUFFIX}"
