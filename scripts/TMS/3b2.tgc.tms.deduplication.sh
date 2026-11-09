#!/bin/bash
#SBATCH --job-name=TMS.3b2.DEDUP
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --array=0-601%60
#SBATCH --output=LOGS/3b2/%x_%A_%a.out
#SBATCH --error=LOGS/3b2/%x_%A_%a.err
#-------------------------------------------------------------------------------
# TreeGeneClimate (TGC) — TMS
# Step 3b2: Removal of PCR duplicates from Bismark alignments (one sample per task)
#
#   samtools sort -n  ->  samtools fixmate -m  ->  samtools sort  ->  samtools markdup -r
#   (pairs with identical outer coordinates and orientation are collapsed to one)
#
# Input (Step 3b; Bismark paired-end alignments, any sort order):
#   DATA/TMS/MAPPED.FILES.TMS/<SAMPLE>${IN_SUFFIX}
# Output:
#   DATA/TMS/MAPPED.FILES.TMS/<SAMPLE>_R1_p_bismark_bt2_pe.dedup.bam (+ .csi)
#   RESULTS/TMS/QC/DEDUP/<SAMPLE>.markdup.txt     (samtools markdup statistics)
# Summary (TGC_MODE=summary, run once after the array):
#   RESULTS/TMS/QC/DEDUP/dedup_summary.tsv            duplicate reads per sample
#   RESULTS/TMS/QC/DEDUP/dedup_summary_by_cohort.tsv  median percentage per cohort
#
# Samples: DATA/METADATA/samples_final_602.txt (array task i -> line i+1)
# Create LOGS/3b2 before submitting.
# Software: SAMtools 1.21
#-------------------------------------------------------------------------------
set -euo pipefail

# === USER CONFIGURATION ===
PROJECT_ROOT="${TGC_PROJECT_ROOT:-/path/to/your/project}"
IN_SUFFIX="${TGC_BAM_IN_SUFFIX:-_R1_p_bismark_bt2_pe.sorted.bam}"
THREADS="${SLURM_CPUS_PER_TASK:-4}"
module purge 2>/dev/null || true
module load gcc/14.2.0 samtools/1.21 2>/dev/null || true
# ===========================

BAM_DIR="${PROJECT_ROOT}/DATA/TMS/MAPPED.FILES.TMS"
QC_DIR="${PROJECT_ROOT}/RESULTS/TMS/QC/DEDUP"
mkdir -p "${QC_DIR}"

# Summary mode: percentage of reads removed as duplicates per sample, and the median
# per cohort (breeding: the samples of breeding_sample2family.txt; natural: the others)
if [[ "${TGC_MODE:-}" == "summary" ]]; then
  map="${PROJECT_ROOT}/DATA/METADATA/breeding_sample2family.txt"
  out="${QC_DIR}/dedup_summary.tsv"
  {
    printf "sample\tcohort\treads\tduplicates\tduplicate_pct\n"
    for f in "${QC_DIR}"/*.markdup.txt; do
      awk -v s="$(basename "${f}" .markdup.txt)" -v map="${map}" '
        BEGIN { FS = ": *"; while ((getline l < map) > 0) { split(l, a, "\t"); b[a[1]] = 1 } }
        /^READ:/ { r = $2 }
        /^DUPLICATE TOTAL:/ { d = $2 }
        END { printf "%s\t%s\t%d\t%d\t%.6f\n", s, (s in b) ? "BREEDING" : "NATURAL", r, d, 100 * d / r }' "${f}"
    done
  } > "${out}"
  {
    printf "cohort\tsamples\tmedian_duplicate_pct\n"
    for c in BREEDING NATURAL; do
      awk -F'\t' -v c="${c}" '$2 == c { print $5 }' "${out}" | sort -g |
        awk -v c="${c}" '{ v[NR] = $1 } END { m = (NR % 2) ? v[(NR + 1) / 2] : (v[NR / 2] + v[NR / 2 + 1]) / 2
                                              printf "%s\t%d\t%.6f\n", c, NR, m }'
    done
  } > "${QC_DIR}/dedup_summary_by_cohort.tsv"
  echo "[$(date)] summary: ${out}"
  exit 0
fi
i="${SLURM_ARRAY_TASK_ID:-0}"
sample=$(sed -n "$((i + 1))p" "${PROJECT_ROOT}/DATA/METADATA/samples_final_602.txt")
[[ -n "${sample}" ]] || { echo "Task ${i}: no sample"; exit 0; }
in_bam="${BAM_DIR}/${sample}${IN_SUFFIX}"
out_bam="${BAM_DIR}/${sample}_R1_p_bismark_bt2_pe.dedup.bam"
[[ -s "${in_bam}" ]] || { echo "ERROR: missing ${in_bam}"; exit 1; }
tmp=$(mktemp -d "${TMPDIR:-/tmp}/dedup_${sample}_XXXX")
trap 'rm -rf "${tmp}"' EXIT

samtools sort -n -@ "${THREADS}" -T "${tmp}/n" -o "${tmp}/n.bam" "${in_bam}"
samtools fixmate -m -@ "${THREADS}" "${tmp}/n.bam" "${tmp}/fix.bam"
samtools sort -@ "${THREADS}" -T "${tmp}/c" -o "${tmp}/c.bam" "${tmp}/fix.bam"
samtools markdup -r -s -@ "${THREADS}" "${tmp}/c.bam" "${out_bam}" 2> "${QC_DIR}/${sample}.markdup.txt"
samtools index -c "${out_bam}"
echo "[$(date)] ${sample}: $(grep -E '^READ:|^DUPLICATE TOTAL:' "${QC_DIR}/${sample}.markdup.txt" | tr '\n' ' ')"
