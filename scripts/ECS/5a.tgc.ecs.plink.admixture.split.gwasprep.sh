#!/bin/bash
#-------------------------------------------------------------------------------
# TreeGeneClimate (TGC) — ECS
# Step 5a: Split final filtered VCF into BREEDING vs NATURAL cohorts, then:
#          (i) GWAS input BED sets (unpruned)
#          (ii) LD pruning (for structure analyses)
#          (iii) PCA on pruned SNPs
#          (iv) ADMIXTURE on pruned SNPs
#
# Project root:
#   /path/to/your/project
#
# Input (results):
#   RESULTS/ECS/VARIANT.CALLING/tgc.ecs.allsamples.call.filt.maf05.snvs.vcf.gz
#
# Cohort sample lists (frozen metadata):
#   DATA/METADATA/ECS/all_breeding_samples.txt
#   DATA/METADATA/ECS/all_natural_samples.txt
#
# Output (results):
#   RESULTS/ECS/VCF_SPLIT/                     (cohort VCFs)
#   RESULTS/ECS/GWAS/PLINK_INPUT/             (unpruned BED sets)
#   RESULTS/ECS/POPGEN/STRUCTURE/             (prune + PCA + ADMIXTURE inputs/outputs)
#-------------------------------------------------------------------------------

#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --job-name=ECS.PLINK.ADMIX
#SBATCH --output=/path/to/your/project/LOGS/%x_%j.out
#SBATCH --error=/path/to/your/project/LOGS/%x_%j.err
#SBATCH --time=48:00:00
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=YOUR_EMAIL

set -euo pipefail

echo "[$(date)] SLURM job started: ${SLURM_JOB_NAME:-no_slurm}"

# ------------------------------------------------------------------------------
# Load modules
module purge
module load gcc/14.2.0
module load bcftools/1.19
module load plink/1.9
module load miniforge3/24.3.0-0
source activate admixture

# ------------------------------------------------------------------------------
# Paths
# === USER CONFIGURATION ===
PROJECT_ROOT="/path/to/your/project"  # <-- set this
# ===========================
LOGS="${PROJECT_ROOT}/LOGS"

VCF_CALL_DIR="${PROJECT_ROOT}/RESULTS/ECS/VARIANT.CALLING"
VCF_ALL="${VCF_CALL_DIR}/tgc.ecs.allsamples.call.filt.maf05.snvs.vcf.gz"

META_DIR="${PROJECT_ROOT}/DATA/METADATA"
BREEDING_SAMPLES="${META_DIR}/all_breeding_samples.txt"
NATURAL_SAMPLES="${META_DIR}/all_natural_samples.txt"

SPLIT_DIR="${PROJECT_ROOT}/RESULTS/ECS/VCF_SPLIT"
VCF_BREEDING="${SPLIT_DIR}/tgc.ecs.breeding.call.filt.maf05.snvs.vcf.gz"
VCF_NATURAL="${SPLIT_DIR}/tgc.ecs.natural.call.filt.maf05.snvs.vcf.gz"

GWAS_DIR="${PROJECT_ROOT}/RESULTS/ECS/GWAS/PLINK_INPUT"
GWAS_BREEDING_DIR="${GWAS_DIR}/BREEDING"
GWAS_NATURAL_DIR="${GWAS_DIR}/NATURAL"

STRUCT_DIR="${PROJECT_ROOT}/RESULTS/ECS/POPGEN/STRUCTURE"
BREEDING_DIR="${STRUCT_DIR}/BREEDING"
NATURAL_DIR="${STRUCT_DIR}/NATURAL"

BREEDING_PCA_DIR="${BREEDING_DIR}/PCA"
NATURAL_PCA_DIR="${NATURAL_DIR}/PCA"

BREEDING_ADMIX_DIR="${BREEDING_DIR}/ADMIXTURE"
NATURAL_ADMIX_DIR="${NATURAL_DIR}/ADMIXTURE"

mkdir -p \
  "${LOGS}" \
  "${SPLIT_DIR}" \
  "${GWAS_BREEDING_DIR}" "${GWAS_NATURAL_DIR}" \
  "${BREEDING_PCA_DIR}" "${NATURAL_PCA_DIR}" \
  "${BREEDING_ADMIX_DIR}" "${NATURAL_ADMIX_DIR}"

# ------------------------------------------------------------------------------
# Sanity checks (fail fast)
if [[ ! -s "${VCF_ALL}" ]]; then
  echo "ERROR: missing input VCF: ${VCF_ALL}"
  exit 1
fi
if [[ ! -s "${BREEDING_SAMPLES}" ]]; then
  echo "ERROR: missing breeding sample list: ${BREEDING_SAMPLES}"
  exit 1
fi
if [[ ! -s "${NATURAL_SAMPLES}" ]]; then
  echo "ERROR: missing natural sample list: ${NATURAL_SAMPLES}"
  exit 1
fi

# ------------------------------------------------------------------------------
# Index VCF (CSI; robust for large coordinates)
if [[ ! -s "${VCF_ALL}.csi" ]]; then
  echo "[$(date)] Indexing input VCF (CSI)"
  bcftools index -c --threads 48 "${VCF_ALL}"
fi

# ------------------------------------------------------------------------------
# Split cohorts
echo "[$(date)] Splitting cohorts"
bcftools view --threads 48 -S "${BREEDING_SAMPLES}" -Oz -o "${VCF_BREEDING}" "${VCF_ALL}"
bcftools view --threads 48 -S "${NATURAL_SAMPLES}"  -Oz -o "${VCF_NATURAL}"  "${VCF_ALL}"

bcftools index -c --threads 48 "${VCF_BREEDING}"
bcftools index -c --threads 48 "${VCF_NATURAL}"

# ------------------------------------------------------------------------------
# GWAS input BED sets (UNPRUNED)
echo "[$(date)] Creating GWAS BED sets (unpruned)"
plink --vcf "${VCF_BREEDING}" --double-id --allow-extra-chr \
  --set-missing-var-ids @:# \
  --threads 48 \
  --make-bed --out "${GWAS_BREEDING_DIR}/tgc.ecs.breeding.gwas.unpruned"

plink --vcf "${VCF_NATURAL}" --double-id --allow-extra-chr \
  --set-missing-var-ids @:# \
  --threads 48 \
  --make-bed --out "${GWAS_NATURAL_DIR}/tgc.ecs.natural.gwas.unpruned"

# ------------------------------------------------------------------------------
# LD pruning (STRUCTURE analyses)
echo "[$(date)] LD pruning (STRUCTURE analyses)"
plink --vcf "${VCF_BREEDING}" --double-id --allow-extra-chr \
  --set-missing-var-ids @:# \
  --threads 48 \
  --indep-pairwise 50 10 0.2 \
  --out "${BREEDING_DIR}/tgc.ecs.breeding.prune"

plink --vcf "${VCF_NATURAL}" --double-id --allow-extra-chr \
  --set-missing-var-ids @:# \
  --threads 48 \
  --indep-pairwise 50 10 0.2 \
  --out "${NATURAL_DIR}/tgc.ecs.natural.prune"

# ------------------------------------------------------------------------------
# PCA on PRUNED SNP sets
echo "[$(date)] PCA on pruned SNP sets"
plink --vcf "${VCF_BREEDING}" --double-id --allow-extra-chr \
  --set-missing-var-ids @:# \
  --threads 48 \
  --extract "${BREEDING_DIR}/tgc.ecs.breeding.prune.prune.in" \
  --pca --out "${BREEDING_PCA_DIR}/tgc.ecs.breeding.pca.pruned"

plink --vcf "${VCF_NATURAL}" --double-id --allow-extra-chr \
  --set-missing-var-ids @:# \
  --threads 48 \
  --extract "${NATURAL_DIR}/tgc.ecs.natural.prune.prune.in" \
  --pca --out "${NATURAL_PCA_DIR}/tgc.ecs.natural.pca.pruned"

# ------------------------------------------------------------------------------
# ADMIXTURE inputs (PRUNED SNP sets) + runs
echo "[$(date)] Preparing ADMIXTURE inputs (pruned SNP sets)"
plink --vcf "${VCF_BREEDING}" --double-id --allow-extra-chr \
  --set-missing-var-ids @:# \
  --threads 48 \
  --extract "${BREEDING_DIR}/tgc.ecs.breeding.prune.prune.in" \
  --make-bed --out "${BREEDING_DIR}/tgc.ecs.breeding.admix.pruned"

plink --vcf "${VCF_NATURAL}" --double-id --allow-extra-chr \
  --set-missing-var-ids @:# \
  --threads 48 \
  --extract "${NATURAL_DIR}/tgc.ecs.natural.prune.prune.in" \
  --make-bed --out "${NATURAL_DIR}/tgc.ecs.natural.admix.pruned"

# ADMIXTURE expects numeric family IDs; set FID=0 in .fam
awk 'BEGIN{OFS="\t"}{$1="0"; print}' "${BREEDING_DIR}/tgc.ecs.breeding.admix.pruned.fam" > "${BREEDING_DIR}/tgc.ecs.breeding.admix.pruned.fam.tmp" \
  && mv "${BREEDING_DIR}/tgc.ecs.breeding.admix.pruned.fam.tmp" "${BREEDING_DIR}/tgc.ecs.breeding.admix.pruned.fam"

awk 'BEGIN{OFS="\t"}{$1="0"; print}' "${NATURAL_DIR}/tgc.ecs.natural.admix.pruned.fam" > "${NATURAL_DIR}/tgc.ecs.natural.admix.pruned.fam.tmp" \
  && mv "${NATURAL_DIR}/tgc.ecs.natural.admix.pruned.fam.tmp" "${NATURAL_DIR}/tgc.ecs.natural.admix.pruned.fam"

# Breeding: K=1..20
echo "[$(date)] ADMIXTURE (breeding) K=1..20"
for K in $(seq 1 20); do
  admixture -j96 --seed=123 "${BREEDING_DIR}/tgc.ecs.breeding.admix.pruned.bed" "${K}" \
    > "${BREEDING_ADMIX_DIR}/K${K}.out" 2> "${BREEDING_ADMIX_DIR}/K${K}.err"
done
grep -h "CV error" "${BREEDING_ADMIX_DIR}"/K*.out > "${BREEDING_ADMIX_DIR}/breeding_cv_errors.txt" || true

# Natural: K=1..10
echo "[$(date)] ADMIXTURE (natural) K=1..10"
for K in $(seq 1 10); do
  admixture -j96 --seed=123 "${NATURAL_DIR}/tgc.ecs.natural.admix.pruned.bed" "${K}" \
    > "${NATURAL_ADMIX_DIR}/K${K}.out" 2> "${NATURAL_ADMIX_DIR}/K${K}.err"
done
grep -h "CV error" "${NATURAL_ADMIX_DIR}"/K*.out > "${NATURAL_ADMIX_DIR}/natural_cv_errors.txt" || true

# ------------------------------------------------------------------------------
conda deactivate || true

echo "[$(date)] DONE."
echo '**************************************************'
exit 0
