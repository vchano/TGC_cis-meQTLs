#!/bin/bash
#-------------------------------------------------------------------------------
# TreeGeneClimate (TGC) — ECS
# Step 4a1: SNV calling (bcftools mpileup + call) from deduplicated BAMs
#
# Project root:
#   /path/to/your/project
#
# Input (data, frozen):
#   DATA/ECS/MAPPED.FILES.ECS/*_bt2.fix.psrt.dedup.bam
#
# Reference (frozen):
#   REFERENCE/Pabies2.0/Picab02_chromosomes_and_unplaced.fa
#
# Output (results):
#   RESULTS/ECS/VARIANT.CALLING/regions/chr_*.vcf.gz
#   RESULTS/ECS/VARIANT.CALLING/tgc.ecs.allsamples.call.unfilt.snvs.vcf.gz
#   RESULTS/ECS/VARIANT.CALLING/bam.list
#   RESULTS/ECS/VARIANT.CALLING/chrom.list
#-------------------------------------------------------------------------------

#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --job-name=ECS.CALL
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

BAM_DIR="${PROJECT_ROOT}/DATA/ECS/MAPPED.FILES.ECS"
REF_FASTA="${PROJECT_ROOT}/REFERENCE/Pabies2.0/Picab02_chromosomes_and_unplaced.fa"

VCF_BASE="${PROJECT_ROOT}/RESULTS/ECS/VARIANT.CALLING"
REGIONS_DIR="${VCF_BASE}/regions"
LOGS_DIR="${PROJECT_ROOT}/LOGS"

MERGED_VCF="${VCF_BASE}/tgc.ecs.allsamples.call.unfilt.snvs.vcf.gz"
BAM_LIST="${VCF_BASE}/bam.list"
CHROM_LIST="${VCF_BASE}/chrom.list"
FAIL_LOG="${VCF_BASE}/calling.fail.log"

mkdir -p "${VCF_BASE}" "${REGIONS_DIR}" "${LOGS_DIR}"
: > "${FAIL_LOG}"

# ------------------------------------------------------------------------------
# CPU layout (full-node, no oversubscription)
CPUS="${SLURM_CPUS_PER_TASK:-192}"
MAX_JOBS=24
THREADS_PER_JOB=$(( CPUS / MAX_JOBS ))   # 192/24 = 8
if [[ "${THREADS_PER_JOB}" -lt 1 ]]; then THREADS_PER_JOB=1; fi

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

TMPDIR_BASE="${SLURM_TMPDIR:-/tmp}"
TMP="${TMPDIR_BASE}/${SLURM_JOB_ID}"
mkdir -p "${TMP}"

echo "[$(date)] CPUs per task: ${CPUS}"
echo "[$(date)] Parallel jobs: ${MAX_JOBS}"
echo "[$(date)] Threads per job (bgzip): ${THREADS_PER_JOB}"
echo "[$(date)] TMP: ${TMP}"

# ------------------------------------------------------------------------------
# Reference index
if [[ ! -s "${REF_FASTA}.fai" ]]; then
  samtools faidx "${REF_FASTA}"
fi
cut -f1 "${REF_FASTA}.fai" > "${CHROM_LIST}"

# ------------------------------------------------------------------------------
# BAM list
find "${BAM_DIR}" -maxdepth 1 -type f -name '*_bt2.fix.psrt.dedup.bam' | sort > "${BAM_LIST}"

if [[ ! -s "${BAM_LIST}" ]]; then
  echo "ERROR: No BAMs found in: ${BAM_DIR}"
  exit 1
fi

BAM_COUNT="$(wc -l < "${BAM_LIST}")"
echo "[$(date)] BAMs found: ${BAM_COUNT}"

# ------------------------------------------------------------------------------
# Parallel per-chromosome calling
JOB_COUNT=0

while read -r chr; do
  (
    echo "[$(date)] Calling: ${chr}"

    out_vcf="${REGIONS_DIR}/chr_${chr}.vcf.gz"

    bcftools mpileup -Ou \
      -f "${REF_FASTA}" \
      -r "${chr}" \
      -b "${BAM_LIST}" \
      -a AD,DP,SP \
      2> "${REGIONS_DIR}/${chr}_mpileup.err" \
    | bcftools call -mv -f GQ,GP \
      --threads "${THREADS_PER_JOB}" \
      -Oz -o "${out_vcf}" \
      2> "${REGIONS_DIR}/${chr}_call.err"

    if [[ ! -s "${out_vcf}" ]]; then
      echo "FAILED: ${chr}" >> "${FAIL_LOG}"
      exit 1
    fi

    # CSI index (safe for large contigs/positions)
    bcftools index -c --threads "${THREADS_PER_JOB}" "${out_vcf}" \
      2> "${REGIONS_DIR}/${chr}_index.err" || true
  ) &

  JOB_COUNT=$((JOB_COUNT + 1))
  if [[ "${JOB_COUNT}" -ge "${MAX_JOBS}" ]]; then
    wait
    JOB_COUNT=0
  fi
done < "${CHROM_LIST}"

wait
echo "[$(date)] Per-chromosome calling complete."

# ------------------------------------------------------------------------------
# Merge per-chromosome VCFs (ordering from chrom.list)
VCF_FILES=()
while read -r chr; do
  f="${REGIONS_DIR}/chr_${chr}.vcf.gz"
  if [[ -s "${f}" ]]; then
    VCF_FILES+=("${f}")
  else
    echo "MISSING_VCF: ${chr}" >> "${FAIL_LOG}"
  fi
done < "${CHROM_LIST}"

if [[ "${#VCF_FILES[@]}" -eq 0 ]]; then
  echo "ERROR: No per-chrom VCFs were produced."
  exit 1
fi

bcftools concat --threads "${CPUS}" -Oz -o "${MERGED_VCF}" "${VCF_FILES[@]}"

# CSI index for merged VCF
bcftools index -c --threads "${CPUS}" "${MERGED_VCF}" || true

echo "[$(date)] Unfiltered merged VCF: ${MERGED_VCF}"
echo "[$(date)] Done."
exit 0
