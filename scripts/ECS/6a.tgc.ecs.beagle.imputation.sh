#!/bin/bash
#-------------------------------------------------------------------------------
# TreeGeneClimate (TGC) — ECS
# Step 6a: (i) Remove monomorphic sites within each cohort (poly-within-cohort)
#          (ii) BEAGLE imputation per cohort
#               - robust to contigs with only 1 variant (skip those for Beagle,
#                 then append them back unchanged)
#
# Project root:
#   /path/to/your/project
#
# Input (from Step 5a split):
#   RESULTS/ECS/VCF_SPLIT/tgc.ecs.breeding.call.filt.maf05.snvs.vcf.gz
#   RESULTS/ECS/VCF_SPLIT/tgc.ecs.natural.call.filt.maf05.snvs.vcf.gz
#
# Output:
#   RESULTS/ECS/VCF_SPLIT/*.poly.vcf.gz
#   RESULTS/ECS/VCF_SPLIT/*.poly.imputed.vcf.gz         (final, merged)
#   RESULTS/ECS/VCF_SPLIT/BEAGLE_TMP/<COHORT>/...        (intermediates)
#-------------------------------------------------------------------------------

#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --job-name=ECS.POLY.BEAGLE
#SBATCH --output=/path/to/your/project/LOGS/%x_%j.out
#SBATCH --error=/path/to/your/project/LOGS/%x_%j.err
#SBATCH --time=48:00:00
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=YOUR_EMAIL

set -euo pipefail

# ----------------------------- helpers ----------------------------------------
ts() { date +"%a %b %d %H:%M:%S %Z %Y"; }

die() { echo "[$(ts)] ERROR: $*" >&2; exit 1; }

need_file() { [[ -s "$1" ]] || die "Missing/empty file: $1"; }

need_cmd() { command -v "$1" >/dev/null 2>&1 || die "Command not found in PATH: $1"; }

count_sites() {
  # fast if indexed; otherwise falls back to counting lines (slower)
  local vcf="$1"
  if [[ -s "${vcf}.csi" || -s "${vcf}.tbi" ]]; then
    bcftools index -n "$vcf"
  else
    bcftools view -H "$vcf" | wc -l
  fi
}

missing_genotypes_sn() {
  # Extract "number of missing genotypes" from bcftools stats (SN line).
  # Returns integer, or "NA" if not found.
  local vcf="$1"
  local val
  val=$(bcftools stats -s - "$vcf" 2>/dev/null | awk -F'\t' '$1=="SN" && $3=="number of missing genotypes:" {print $4; exit}')
  [[ -n "${val:-}" ]] && echo "$val" || echo "NA"
}

make_contig_lists_by_variant_count() {
  # Writes two files:
  #   contigs_ge2.txt : contigs with >=2 variants
  #   contigs_eq1.txt : contigs with exactly 1 variant
  local vcf="$1"
  local out_ge2="$2"
  local out_eq1="$3"

  bcftools query -f '%CHROM\n' "$vcf" \
    | sort \
    | uniq -c \
    | awk '
        $1==1 {print $2 > eq1}
        $1>=2 {print $2 > ge2}
      ' ge2="$out_ge2" eq1="$out_eq1"
}

# ----------------------------- environment ------------------------------------
echo "================================================================================"
echo "JobID = ${SLURM_JOB_ID:-NA}"
echo "User = ${SLURM_JOB_USER:-$USER}, Account = ${SLURM_JOB_ACCOUNT:-NA}"
echo "Partition = ${SLURM_JOB_PARTITION:-NA}, Nodelist = ${SLURM_JOB_NODELIST:-NA}"
echo "================================================================================"

echo "[$(ts)] SLURM job started: ${SLURM_JOB_NAME:-no_slurm}"
echo "[$(ts)] Host: $(hostname)"
echo "[$(ts)] CPUs: ${SLURM_CPUS_PER_TASK:-1}"

module purge
module load gcc/14.2.0
module load bcftools/1.19

module load miniforge3/24.3.0-0
# robust conda init in batch
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate beagle

need_cmd bcftools
need_cmd beagle

echo "[$(ts)] bcftools: $(bcftools --version | head -n 1)"
echo "[$(ts)] beagle:   $(which beagle)"
echo "[$(ts)] beagle bin: $(file -b "$(which beagle)" || true)"

# ----------------------------- config -----------------------------------------
# === USER CONFIGURATION ===
PROJECT_ROOT="/path/to/your/project"  # <-- set this
# ===========================
LOGS="${PROJECT_ROOT}/LOGS"
SPLIT_DIR="${PROJECT_ROOT}/RESULTS/ECS/VCF_SPLIT"
TMP_BASE="${SPLIT_DIR}/BEAGLE_TMP"

mkdir -p "${LOGS}" "${TMP_BASE}"

BREEDING_IN="${SPLIT_DIR}/tgc.ecs.breeding.call.filt.maf05.snvs.vcf.gz"
NATURAL_IN="${SPLIT_DIR}/tgc.ecs.natural.call.filt.maf05.snvs.vcf.gz"

need_file "${BREEDING_IN}"
need_file "${NATURAL_IN}"

THREADS="${SLURM_CPUS_PER_TASK:-1}"

# Java memory: Beagle wrapper uses Java; _JAVA_OPTIONS is honored by Java.
JAVA_MEM="700g"
export _JAVA_OPTIONS="-Xmx${JAVA_MEM}"

echo "[$(ts)] Java opts: ${_JAVA_OPTIONS}"

# ----------------------------- main -------------------------------------------
run_cohort() {
  local cohort="$1"          # BREEDING / NATURAL
  local vcf_in="$2"          # input VCF (split, cohort)
  local out_prefix="$3"      # output prefix (full path, without extensions)

  local tmp_dir="${TMP_BASE}/${cohort}"
  mkdir -p "${tmp_dir}"

  echo "--------------------------------------------------------------------------------"
  echo "[$(ts)] COHORT: ${cohort}"
  echo "[$(ts)] Input:  ${vcf_in}"
  echo "--------------------------------------------------------------------------------"

  echo "[$(ts)] Indexing input VCF (CSI if missing)"
  if [[ ! -s "${vcf_in}.csi" ]]; then
    bcftools index -c --threads "${THREADS}" "${vcf_in}"
  fi

  echo "[$(ts)] Counting sites BEFORE within-cohort polymorphic filtering"
  local n_before
  n_before=$(count_sites "${vcf_in}")
  echo "  ${cohort} (all sites in split VCF): ${n_before}"

  # 1) poly-within-cohort filter
  local vcf_poly="${out_prefix}.poly.vcf.gz"
  echo "[$(ts)] Filtering to polymorphic-within-cohort sites (AC>0 && AC<AN)"
  bcftools view --threads "${THREADS}" -Oz \
    --include 'AC>0 && AC<AN' \
    -o "${vcf_poly}" \
    "${vcf_in}"

  bcftools index -c --threads "${THREADS}" "${vcf_poly}"

  echo "[$(ts)] Counting sites AFTER within-cohort polymorphic filtering"
  local n_poly
  n_poly=$(count_sites "${vcf_poly}")
  echo "  ${cohort} (poly): ${n_poly}"

  # Missingness before imputation (on poly set)
  echo "[$(ts)] Missing genotypes BEFORE imputation (poly set)"
  local miss_pre
  miss_pre=$(missing_genotypes_sn "${vcf_poly}")
  echo "  ${cohort} missing GT count (poly, pre-impute): ${miss_pre}"

  # 2) prevent Beagle crash on contigs with only 1 variant
  local contigs_ge2="${tmp_dir}/${cohort}.contigs_ge2.txt"
  local contigs_eq1="${tmp_dir}/${cohort}.contigs_eq1.txt"

  echo "[$(ts)] Identifying contigs with >=2 variants vs exactly 1 variant (poly set)"
  : > "${contigs_ge2}"
  : > "${contigs_eq1}"
  make_contig_lists_by_variant_count "${vcf_poly}" "${contigs_ge2}" "${contigs_eq1}"

  local n_contig_ge2 n_contig_eq1
  n_contig_ge2=$(wc -l < "${contigs_ge2}" || echo 0)
  n_contig_eq1=$(wc -l < "${contigs_eq1}" || echo 0)

  echo "  ${cohort} contigs with >=2 variants: ${n_contig_ge2}"
  echo "  ${cohort} contigs with  1 variant : ${n_contig_eq1}"

  # Build subset VCFs
  local vcf_ge2="${tmp_dir}/${cohort}.poly.ge2.vcf.gz"
  local vcf_eq1="${tmp_dir}/${cohort}.poly.eq1.vcf.gz"

  if [[ "${n_contig_ge2}" -gt 0 ]]; then
    echo "[$(ts)] Subsetting to contigs with >=2 variants (for Beagle)"
    # IMPORTANT: contigs_ge2 is a 1-column contig list, so use -r with comma-separated contigs
    local regions_ge2
    regions_ge2=$(paste -sd, "${contigs_ge2}")
    bcftools view --threads "${THREADS}" -Oz \
      -r "${regions_ge2}" \
      -o "${vcf_ge2}" \
      "${vcf_poly}"
    bcftools index -c --threads "${THREADS}" "${vcf_ge2}"
  else
    die "${cohort}: No contigs with >=2 variants found; nothing to impute."
  fi

  if [[ "${n_contig_eq1}" -gt 0 ]]; then
    echo "[$(ts)] Subsetting to contigs with exactly 1 variant (will NOT be imputed; appended back later)"
    # IMPORTANT: contigs_eq1 is a 1-column contig list, so use -r with comma-separated contigs
    local regions_eq1
    regions_eq1=$(paste -sd, "${contigs_eq1}")
    bcftools view --threads "${THREADS}" -Oz \
      -r "${regions_eq1}" \
      -o "${vcf_eq1}" \
      "${vcf_poly}"
    bcftools index -c --threads "${THREADS}" "${vcf_eq1}"
  else
    echo "[$(ts)] No 1-variant contigs for ${cohort} (good)."
  fi

  # missingness specifically in the set that will be imputed
  echo "[$(ts)] Missing genotypes BEFORE imputation (subset sent to Beagle, >=2 variants/contig)"
  local miss_pre_ge2
  miss_pre_ge2=$(missing_genotypes_sn "${vcf_ge2}")
  echo "  ${cohort} missing GT count (ge2 subset, pre-impute): ${miss_pre_ge2}"

  # 3) Beagle on ge2 subset
  local beagle_out_prefix="${tmp_dir}/${cohort}.poly.ge2.imputed"
  local vcf_ge2_imputed="${beagle_out_prefix}.vcf.gz"

  echo "[$(ts)] Running BEAGLE on >=2-variant contigs"
  echo "  Threads: ${THREADS}"
  echo "  Input:   ${vcf_ge2}"
  echo "  Output:  ${vcf_ge2_imputed}"

  beagle \
    gt="${vcf_ge2}" \
    out="${beagle_out_prefix}" \
    nthreads="${THREADS}"

  [[ -s "${vcf_ge2_imputed}" ]] || die "${cohort}: Beagle did not produce output VCF: ${vcf_ge2_imputed}"

  bcftools index -c --threads "${THREADS}" "${vcf_ge2_imputed}"

  echo "[$(ts)] Missing genotypes AFTER imputation (ge2 subset)"
  local miss_post_ge2
  miss_post_ge2=$(missing_genotypes_sn "${vcf_ge2_imputed}")
  echo "  ${cohort} missing GT count (ge2 subset, post-impute): ${miss_post_ge2}"

  # 4) Merge imputed ge2 subset + untouched eq1 subset back into a final poly.imputed VCF
  local vcf_final="${out_prefix}.poly.imputed.vcf.gz"

  if [[ "${n_contig_eq1}" -gt 0 ]]; then
    echo "[$(ts)] Merging imputed (ge2) + untouched (eq1) and sorting"
    bcftools concat -a -Oz \
      "${vcf_ge2_imputed}" \
      "${vcf_eq1}" \
      | bcftools sort -Oz -o "${vcf_final}" -
  else
    echo "[$(ts)] No eq1 subset; final = imputed ge2 (sorted anyway)"
    bcftools sort -Oz -o "${vcf_final}" "${vcf_ge2_imputed}"
  fi

  bcftools index -c --threads "${THREADS}" "${vcf_final}"

  echo "[$(ts)] Final counts and missingness (poly.imputed)"
  local n_final miss_post_all
  n_final=$(count_sites "${vcf_final}")
  miss_post_all=$(missing_genotypes_sn "${vcf_final}")
  echo "  ${cohort} sites (final poly.imputed): ${n_final}"
  echo "  ${cohort} missing GT count (final poly.imputed): ${miss_post_all}"

  if [[ "${n_contig_eq1}" -gt 0 ]]; then
    local n_eq1_sites miss_eq1
    n_eq1_sites=$(count_sites "${vcf_eq1}")
    miss_eq1=$(missing_genotypes_sn "${vcf_eq1}")
    echo "  ${cohort} sites on 1-variant contigs (not imputed): ${n_eq1_sites}"
    echo "  ${cohort} missing GT count on 1-variant contigs:      ${miss_eq1}"
  fi

  echo "[$(ts)] COHORT ${cohort} DONE: ${vcf_final}"
}

echo "[$(ts)] Starting Step 6a: within-cohort poly filtering + Beagle imputation"

BREEDING_PREFIX="${SPLIT_DIR}/tgc.ecs.breeding.call.filt.maf05.snvs"
NATURAL_PREFIX="${SPLIT_DIR}/tgc.ecs.natural.call.filt.maf05.snvs"

run_cohort "BREEDING" "${BREEDING_IN}" "${BREEDING_PREFIX}"
run_cohort "NATURAL"  "${NATURAL_IN}"  "${NATURAL_PREFIX}"

conda deactivate || true

echo "[$(ts)] ALL DONE."
echo "**************************************************"
exit 0
