#!/bin/bash
#SBATCH --job-name=ECS.4a2.FILT
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=LOGS/%x_%j.out
#SBATCH --error=LOGS/%x_%j.err
#-------------------------------------------------------------------------------
# TreeGeneClimate (TGC) — ECS
# Step 4a2: Sample renaming and SNP filtering of the joint call set
#
#   (i)   Rename samples in the merged VCF from BAM paths to sample IDs
#         (e.g. .../P001_WA02_bt2.fix.psrt.dedup.bam -> P001_WA02)
#   (ii)  Keep biallelic SNPs with
#           F_MISSING <= 0.20, MAF >= 0.05, INFO/DP >= 10 (summed site depth)
#         computed over all 602 samples (both cohorts jointly)
#   (iii) Site counts at each stage
#
# Submit from the project root (sbatch SCRIPTS/ECS/4a2...), with
# TGC_PROJECT_ROOT set or PROJECT_ROOT edited below.
#
# Input (Step 4a1):
#   RESULTS/ECS/VARIANT.CALLING/tgc.ecs.allsamples.call.unfilt.snvs.vcf.gz
#
# Output:
#   RESULTS/ECS/VARIANT.CALLING/tgc.ecs.allsamples.call.unfilt.snvs.renamed.vcf.gz
#   RESULTS/ECS/VARIANT.CALLING/tgc.ecs.allsamples.call.filt.maf05.snvs.vcf.gz
#   RESULTS/ECS/VARIANT.CALLING/STATS/
#     unfiltered.bcftools.stats.txt, filtered_maf05.bcftools.stats.txt
#     ecs_filtering_counts.tsv
#
# Software: BCFtools 1.19
#-------------------------------------------------------------------------------
set -euo pipefail

# === USER CONFIGURATION ===
PROJECT_ROOT="${TGC_PROJECT_ROOT:-/path/to/your/project}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"
# Software environment (edit for your system)
module purge 2>/dev/null || true
module load gcc/14.2.0 bcftools/1.19 2>/dev/null || true
# ===========================

VCF_DIR="${PROJECT_ROOT}/RESULTS/ECS/VARIANT.CALLING"
STATS_DIR="${VCF_DIR}/STATS"
RAW_VCF="${VCF_DIR}/tgc.ecs.allsamples.call.unfilt.snvs.vcf.gz"
REN_VCF="${VCF_DIR}/tgc.ecs.allsamples.call.unfilt.snvs.renamed.vcf.gz"
FILT_VCF="${VCF_DIR}/tgc.ecs.allsamples.call.filt.maf05.snvs.vcf.gz"
COUNTS="${STATS_DIR}/ecs_filtering_counts.tsv"

mkdir -p "${STATS_DIR}"
[[ -s "${RAW_VCF}" ]] || { echo "ERROR: missing input ${RAW_VCF}"; exit 1; }
echo "[$(date)] Step 4a2 started (threads=${THREADS})"

# ------------------------------------------------------------------------------
# (i) Rename samples: strip directory and BAM suffix
# ------------------------------------------------------------------------------
NAMES="${STATS_DIR}/sample_rename.tsv"
bcftools query -l "${RAW_VCF}" \
  | awk 'BEGIN{OFS="\t"} {n=$0; sub(/^.*\//, "", n); sub(/_bt2\.fix\.psrt\.dedup\.bam$/, "", n); print $0, n}' \
  > "${NAMES}"
bcftools reheader --samples "${NAMES}" --threads "${THREADS}" -o "${REN_VCF}" "${RAW_VCF}"
bcftools index -c --threads "${THREADS}" "${REN_VCF}"
echo "[$(date)] Renamed $(wc -l < "${NAMES}") samples"

# ------------------------------------------------------------------------------
# (ii) Joint filtering
# ------------------------------------------------------------------------------
bcftools +fill-tags "${REN_VCF}" -Ou -- -t MAF,F_MISSING,AF \
| bcftools view \
    --types snps --min-alleles 2 --max-alleles 2 \
    --include 'F_MISSING<=0.2 && MAF>=0.05 && INFO/DP>=10' \
    --threads "${THREADS}" -Oz -o "${FILT_VCF}"
bcftools index -c --threads "${THREADS}" "${FILT_VCF}"
echo "[$(date)] Filtered VCF written"

# ------------------------------------------------------------------------------
# (iii) Counts
# ------------------------------------------------------------------------------
bcftools stats --threads "${THREADS}" "${REN_VCF}"  > "${STATS_DIR}/unfiltered.bcftools.stats.txt"
bcftools stats --threads "${THREADS}" "${FILT_VCF}" > "${STATS_DIR}/filtered_maf05.bcftools.stats.txt"

sn () { awk -F'\t' -v k="$2" '$1=="SN" && $3==k {print $4}' "$1"; }
{
  echo -e "stage\trecords\tsnps\tmultiallelic_snp_sites"
  for s in unfiltered filtered_maf05; do
    f="${STATS_DIR}/${s}.bcftools.stats.txt"
    echo -e "${s}\t$(sn "$f" 'number of records:')\t$(sn "$f" 'number of SNPs:')\t$(sn "$f" 'number of multiallelic SNP sites:')"
  done
} > "${COUNTS}"
cat "${COUNTS}"

echo "[$(date)] Step 4a2 done"
