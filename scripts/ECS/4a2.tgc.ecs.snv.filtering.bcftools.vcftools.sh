#!/bin/bash
#-------------------------------------------------------------------------------
# TreeGeneClimate (TGC) — ECS
# Step 4a2: SNV filtering + basic stats (bcftools) from merged unfiltered VCF
#
# Project root:
#   /path/to/your/project
#
# Input (results):
#   RESULTS/ECS/VARIANT.CALLING/tgc.ecs.allsamples.call.unfilt.snvs.vcf.gz
#
# Output (results):
#   RESULTS/ECS/VARIANT.CALLING/tgc.ecs.allsamples.call.filt.maf10.snvs.vcf.gz
#   RESULTS/ECS/VARIANT.CALLING/tgc.ecs.allsamples.call.filt.maf01.snvs.vcf.gz
#   RESULTS/ECS/VARIANT.CALLING/STATS_ALL/
#-------------------------------------------------------------------------------

#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH -n 96
#SBATCH -N 1
#SBATCH --job-name=ECS.FILT
#SBATCH --output=/path/to/your/project/LOGS/%x_%j.out
#SBATCH --error=/path/to/your/project/LOGS/%x_%j.err
#SBATCH --time=48:00:00
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=YOUR_EMAIL

set -euo pipefail

echo "[$(date)] SLURM job started: ${SLURM_JOB_NAME:-no_slurm}"

module purge
module load gcc/14.2.0
module load bcftools/1.19
module load samtools/1.21

# === USER CONFIGURATION ===
PROJECT_ROOT="/path/to/your/project"  # <-- set this
# ===========================

VCF_BASE="${PROJECT_ROOT}/RESULTS/ECS/VARIANT.CALLING"
STATS_DIR="${VCF_BASE}/STATS_ALL"
LOGS_DIR="${PROJECT_ROOT}/LOGS"

UNFILT_VCF="${VCF_BASE}/tgc.ecs.allsamples.call.unfilt.snvs.renamed.vcf.gz"
FILT_MAF10_VCF="${VCF_BASE}/tgc.ecs.allsamples.call.filt.maf10.snvs.vcf.gz"
FILT_MAF01_VCF="${VCF_BASE}/tgc.ecs.allsamples.call.filt.maf01.snvs.vcf.gz"

mkdir -p "${STATS_DIR}" "${LOGS_DIR}"

if [[ ! -s "${UNFILT_VCF}" ]]; then
  echo "ERROR: Input VCF not found: ${UNFILT_VCF}"
  exit 1
fi

echo "[$(date)] Counting unfiltered SNVs..."
bcftools view -H "${UNFILT_VCF}" | wc -l

# ------------------------------------------------------------------------------
echo "[$(date)] Filtering (MAF >= 0.10)..."
bcftools +fill-tags "${UNFILT_VCF}" -Ou -- -t MAF,F_MISSING,AF \
| bcftools view \
    --types snps \
    --min-alleles 2 \
    --max-alleles 2 \
    --include 'F_MISSING<=0.2 && MAF>=0.10 && INFO/DP>=10' \
    --threads 48 \
    -Oz -o "${FILT_MAF10_VCF}"

# CSI index (safe for large contigs/positions)
bcftools index -c "${FILT_MAF10_VCF}"

echo "[$(date)] Filtering (MAF >= 0.01)..."
bcftools +fill-tags "${UNFILT_VCF}" -Ou -- -t MAF,F_MISSING,AF \
| bcftools view \
    --types snps \
    --min-alleles 2 \
    --max-alleles 2 \
    --include 'F_MISSING<=0.2 && MAF>=0.01 && INFO/DP>=10' \
    --threads 48 \
    -Oz -o "${FILT_MAF01_VCF}"

# CSI index (safe for large contigs/positions)
bcftools index -c "${FILT_MAF01_VCF}"

echo "[$(date)] Counting filtered SNVs (MAF >= 0.10)..."
bcftools view -H "${FILT_MAF10_VCF}" | wc -l

echo "[$(date)] Counting filtered SNVs (MAF >= 0.01)..."
bcftools view -H "${FILT_MAF01_VCF}" | wc -l

# ------------------------------------------------------------------------------
echo "[$(date)] Writing stats tables (unfiltered vs MAF01 filtered)..."

bcftools +fill-tags "${UNFILT_VCF}" -- -t AF,MAF,F_MISSING \
| bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\t%MAF\t%F_MISSING\n' \
> "${STATS_DIR}/unfiltered.af_maf_missing.txt"

bcftools +fill-tags "${FILT_MAF01_VCF}" -- -t AF,MAF,F_MISSING \
| bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\t%MAF\t%F_MISSING\n' \
> "${STATS_DIR}/filtered_maf01.af_maf_missing.txt"

bcftools query -f '%CHROM\t%POS\t%QUAL\n' "${UNFILT_VCF}" > "${STATS_DIR}/unfiltered.qual.txt"
bcftools query -f '%CHROM\t%POS\t%QUAL\n' "${FILT_MAF01_VCF}" > "${STATS_DIR}/filtered_maf01.qual.txt"

bcftools query -f '%CHROM\t%POS\t%INFO/DP\n' "${UNFILT_VCF}" > "${STATS_DIR}/unfiltered.site_dp.txt"
bcftools query -f '%CHROM\t%POS\t%INFO/DP\n' "${FILT_MAF1_VCF}" > "${STATS_DIR}/filtered_maf01.site_dp.txt"

bcftools stats -s - "${UNFILT_VCF}" > "${STATS_DIR}/unfiltered.bcftools.stats.txt"
bcftools stats -s - "${FILT_MAF01_VCF}" > "${STATS_DIR}/filtered_maf01.bcftools.stats.txt"

echo "[$(date)] Outputs:"
echo "  Unfiltered: ${UNFILT_VCF}"
echo "  Filtered MAF10: ${FILT_MAF10_VCF}"
echo "  Filtered MAF01: ${FILT_MAF01_VCF}"
echo "  Stats: ${STATS_DIR}"

echo "[$(date)] Done."
exit 0
