#!/bin/bash
#-------------------------------------------------------------------------------
# TreeGeneClimate (TGC) — TBS
# Step 4b2: Coverage filtering of Bismark coverage files
#   Keep sites with 5 <= coverage <= 50
#   Separate outputs by cohort (BREEDING vs NATURAL)
#
# Input:
#   RESULTS/TBS/METHYLATION_CALLS/<SAMPLE>/*.bismark.cov.gz
#
# Output:
#   RESULTS/TBS/METHYLATION_FILTERED/<COHORT>/<SAMPLE>.bismark.min5.max50.cov.gz
#-------------------------------------------------------------------------------

#SBATCH -p YOUR_PARTITION
#SBATCH -t 00-12:00:00
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH --job-name=TBS.COVFILTER
#SBATCH --output=/path/to/your/project/LOGS/%x_%j.out
#SBATCH --error=/path/to/your/project/LOGS/%x_%j.err
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=YOUR_EMAIL

set -euo pipefail

# === USER CONFIGURATION ===
PROJECT_ROOT="/path/to/your/project"  # <-- set this
# ===========================
INPUT="${PROJECT_ROOT}/RESULTS/TBS/METHYLATION_CALLS"
OUTBASE="${PROJECT_ROOT}/RESULTS/TBS/METHYLATION_FILTERED"

MINCOV=5
MAXCOV=50

mkdir -p "${OUTBASE}/BREEDING" "${OUTBASE}/NATURAL"

echo "Starting coverage filtering: ${MINCOV} <= cov <= ${MAXCOV}"
echo "Input:  ${INPUT}"
echo "Output: ${OUTBASE}"
echo "----------------------------------"

shopt -s nullglob
files=( ${INPUT}/*/*.bismark.cov.gz )

if [[ ${#files[@]} -eq 0 ]]; then
  echo "ERROR: No .bismark.cov.gz files found under ${INPUT}/*/"
  exit 1
fi

for file in "${files[@]}"; do
  base=$(basename "$file" .bismark.cov.gz)   # e.g. P001_WA02_R1_p_bismark_bt2_pe
  sample="${base%%_R1_p_bismark_bt2_pe}"     # e.g. P001_WA02  (adjust if needed)

  if [[ "$sample" =~ ^P00[1-3]_ ]]; then
    cohort="BREEDING"
  elif [[ "$sample" =~ ^P00[4-8]_ ]]; then
    cohort="NATURAL"
  else
    echo "Skipping unknown cohort: ${sample} (from ${base})"
    continue
  fi

  outfile="${OUTBASE}/${cohort}/${base}.min${MINCOV}.max${MAXCOV}.cov.gz"

  echo "Processing ${base}  ->  ${cohort}"

  zcat "$file" \
    | awk -v lo="${MINCOV}" -v hi="${MAXCOV}" 'BEGIN{OFS="\t"} {cov=$5+$6; if(cov>=lo && cov<=hi) print $0}' \
    | gzip > "$outfile"
done

echo "Filtering completed."
exit 0
