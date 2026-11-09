#!/bin/bash
#SBATCH --job-name=ECS.5a.SPLIT
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=LOGS/%x_%j.out
#SBATCH --error=LOGS/%x_%j.err
#-------------------------------------------------------------------------------
# TreeGeneClimate (TGC) — ECS
# Step 5a: Split the filtered joint VCF into BREEDING and NATURAL cohorts, then
#          (i)   unpruned PLINK BED set per cohort
#          (ii)  LD pruning (PLINK --indep-pairwise 50 10 0.2)
#          (iii) PCA on pruned SNPs (PLINK)
#          (iv)  pruned PLINK BED set per cohort (input to Step 7a)
#
# Input (Step 4a2):
#   RESULTS/ECS/VARIANT.CALLING/tgc.ecs.allsamples.call.filt.maf05.snvs.vcf.gz
# Cohort sample lists:
#   DATA/METADATA/all_breeding_samples.txt
#   DATA/METADATA/all_natural_samples.txt
#
# Output:
#   RESULTS/ECS/VCF_SPLIT/tgc.ecs.<cohort>.call.filt.maf05.snvs.vcf.gz
#   RESULTS/ECS/POPGEN/STRUCTURE/<COHORT>/
#     tgc.ecs.<cohort>.unpruned.{bed,bim,fam}
#     tgc.ecs.<cohort>.prune.prune.{in,out}
#     tgc.ecs.<cohort>.pruned.{bed,bim,fam}
#     PCA/tgc.ecs.<cohort>.pca.pruned.{eigenvec,eigenval}
#
# Software: BCFtools 1.19, PLINK 1.9
#-------------------------------------------------------------------------------
set -euo pipefail

# === USER CONFIGURATION ===
PROJECT_ROOT="${TGC_PROJECT_ROOT:-/path/to/your/project}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"
PRUNE_WINDOW=50   # SNPs
PRUNE_STEP=10     # SNPs
PRUNE_R2=0.2
# Software environment (edit for your system)
module purge 2>/dev/null || true
module load gcc/14.2.0 bcftools/1.19 plink/1.9 2>/dev/null || true
# ===========================

VCF_ALL="${PROJECT_ROOT}/RESULTS/ECS/VARIANT.CALLING/tgc.ecs.allsamples.call.filt.maf05.snvs.vcf.gz"
META_DIR="${PROJECT_ROOT}/DATA/METADATA"
SPLIT_DIR="${PROJECT_ROOT}/RESULTS/ECS/VCF_SPLIT"
STRUCT_DIR="${PROJECT_ROOT}/RESULTS/ECS/POPGEN/STRUCTURE"

[[ -s "${VCF_ALL}" ]] || { echo "ERROR: missing input ${VCF_ALL}"; exit 1; }
mkdir -p "${SPLIT_DIR}"
echo "[$(date)] Step 5a started (threads=${THREADS})"

plink_vcf () {  # common PLINK options for VCF input
  plink --vcf "$1" --double-id --allow-extra-chr --set-missing-var-ids @:# --threads "${THREADS}" "${@:2}"
}

for cohort in breeding natural; do
  COH=$(echo "${cohort}" | tr '[:lower:]' '[:upper:]')
  OUT="${STRUCT_DIR}/${COH}"
  mkdir -p "${OUT}/PCA"
  samples="${META_DIR}/all_${cohort}_samples.txt"
  [[ -s "${samples}" ]] || { echo "ERROR: missing ${samples}"; exit 1; }
  vcf="${SPLIT_DIR}/tgc.ecs.${cohort}.call.filt.maf05.snvs.vcf.gz"

  echo "[$(date)] ${COH}: split"
  bcftools view --threads "${THREADS}" -S "${samples}" -Oz -o "${vcf}" "${VCF_ALL}"
  bcftools index -c --threads "${THREADS}" "${vcf}"

  echo "[$(date)] ${COH}: unpruned BED set"
  plink_vcf "${vcf}" --make-bed --out "${OUT}/tgc.ecs.${cohort}.unpruned"

  echo "[$(date)] ${COH}: LD pruning (${PRUNE_WINDOW} ${PRUNE_STEP} ${PRUNE_R2})"
  plink_vcf "${vcf}" --indep-pairwise "${PRUNE_WINDOW}" "${PRUNE_STEP}" "${PRUNE_R2}" \
    --out "${OUT}/tgc.ecs.${cohort}.prune"

  echo "[$(date)] ${COH}: pruned BED set and PCA"
  plink_vcf "${vcf}" --extract "${OUT}/tgc.ecs.${cohort}.prune.prune.in" \
    --make-bed --out "${OUT}/tgc.ecs.${cohort}.pruned"
  plink_vcf "${vcf}" --extract "${OUT}/tgc.ecs.${cohort}.prune.prune.in" \
    --pca --out "${OUT}/PCA/tgc.ecs.${cohort}.pca.pruned"

  echo "[$(date)] ${COH}: $(wc -l < "${OUT}/tgc.ecs.${cohort}.prune.prune.in") SNPs after pruning"
done

echo "[$(date)] Step 5a done"
