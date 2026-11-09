#!/bin/bash
#SBATCH --job-name=ECS.6a.BEAGLE
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=LOGS/%x_%j.out
#SBATCH --error=LOGS/%x_%j.err
#-------------------------------------------------------------------------------
# TreeGeneClimate (TGC) — ECS
# Step 6a: Within-cohort polymorphic sites and BEAGLE imputation, per cohort
#   (i)   keep sites polymorphic within the cohort (AC > 0 && AC < AN)
#   (ii)  impute and phase with BEAGLE (contigs with >= 2 variants)
#   (iii) append variants on single-variant contigs unchanged (BEAGLE cannot
#         process contigs with one variant) and sort
#
# BEAGLE results depend on the number of threads (default seed); nthreads is fixed
# at 48 so that the imputed genotypes are reproducible.
#
# Input (Step 5a):
#   RESULTS/ECS/VCF_SPLIT/tgc.ecs.<cohort>.call.filt.maf05.snvs.vcf.gz
#
# Output:
#   RESULTS/ECS/VCF_SPLIT/tgc.ecs.<cohort>.call.filt.maf05.snvs.poly.vcf.gz
#   RESULTS/ECS/VCF_SPLIT/tgc.ecs.<cohort>.call.filt.maf05.snvs.poly.imputed.vcf.gz
#   RESULTS/ECS/VCF_SPLIT/imputation_counts.tsv
#   RESULTS/ECS/VCF_SPLIT/BEAGLE_TMP/<COHORT>/   (intermediates)
#
# Software: BCFtools 1.19, BEAGLE 5 (conda environment "beagle")
#-------------------------------------------------------------------------------
set -euo pipefail

# === USER CONFIGURATION ===
PROJECT_ROOT="${TGC_PROJECT_ROOT:-/path/to/your/project}"
THREADS="${SLURM_CPUS_PER_TASK:-16}"   # bcftools
BEAGLE_THREADS=48                      # fixed: results depend on the number of threads
export _JAVA_OPTIONS="-Xmx100g"
# Software environment (edit for your system)
module purge 2>/dev/null || true
module load gcc/14.2.0 bcftools/1.19 miniforge3/24.3.0-0 2>/dev/null || true
source "$(conda info --base)/etc/profile.d/conda.sh" 2>/dev/null || true
conda activate beagle 2>/dev/null || true
# ===========================

SPLIT_DIR="${PROJECT_ROOT}/RESULTS/ECS/VCF_SPLIT"
TMP_BASE="${SPLIT_DIR}/BEAGLE_TMP"
COUNTS="${SPLIT_DIR}/imputation_counts.tsv"
mkdir -p "${TMP_BASE}"

die () { echo "[$(date)] ERROR: $*" >&2; exit 1; }
command -v beagle >/dev/null || die "beagle not found in PATH"
n_sites () { bcftools index -n "$1"; }

echo "[$(date)] Step 6a started (BEAGLE threads=${BEAGLE_THREADS}); $(beagle 2>&1 | head -n 1)"
echo -e "cohort\tsites_split\tsites_poly\tcontigs_imputed\tcontigs_single_variant\tsites_single_variant\tsites_final" > "${COUNTS}"

for cohort in breeding natural; do
  COH=$(echo "${cohort}" | tr '[:lower:]' '[:upper:]')
  prefix="${SPLIT_DIR}/tgc.ecs.${cohort}.call.filt.maf05.snvs"
  vcf_in="${prefix}.vcf.gz"; vcf_poly="${prefix}.poly.vcf.gz"; vcf_final="${prefix}.poly.imputed.vcf.gz"
  tmp="${TMP_BASE}/${COH}"; mkdir -p "${tmp}"
  [[ -s "${vcf_in}" ]] || die "missing ${vcf_in}"

  # (i) polymorphic within cohort (AC/AN were recomputed for the cohort in Step 5a)
  bcftools view --threads "${THREADS}" --include 'AC>0 && AC<AN' -Oz -o "${vcf_poly}" "${vcf_in}"
  bcftools index -c --threads "${THREADS}" "${vcf_poly}"

  # (ii) contigs with >= 2 variants go to BEAGLE; single-variant contigs are kept as they are
  bcftools query -f '%CHROM\n' "${vcf_poly}" | sort | uniq -c \
    | awk -v ge2="${tmp}/contigs_ge2.txt" -v eq1="${tmp}/contigs_eq1.txt" \
        '$1>=2 {print $2 > ge2} $1==1 {print $2 > eq1}'
  touch "${tmp}/contigs_ge2.txt" "${tmp}/contigs_eq1.txt"
  bcftools view --threads "${THREADS}" -r "$(paste -sd, "${tmp}/contigs_ge2.txt")" \
    -Oz -o "${tmp}/poly.ge2.vcf.gz" "${vcf_poly}"
  bcftools index -c "${tmp}/poly.ge2.vcf.gz"

  beagle gt="${tmp}/poly.ge2.vcf.gz" out="${tmp}/poly.ge2.imputed" nthreads="${BEAGLE_THREADS}"
  [[ -s "${tmp}/poly.ge2.imputed.vcf.gz" ]] || die "${COH}: BEAGLE produced no output"
  bcftools index -c "${tmp}/poly.ge2.imputed.vcf.gz"

  # (iii) append single-variant contigs and sort
  n_eq1=$(wc -l < "${tmp}/contigs_eq1.txt"); s_eq1=0
  if [[ "${n_eq1}" -gt 0 ]]; then
    bcftools view --threads "${THREADS}" -r "$(paste -sd, "${tmp}/contigs_eq1.txt")" \
      -Oz -o "${tmp}/poly.eq1.vcf.gz" "${vcf_poly}"
    bcftools index -c "${tmp}/poly.eq1.vcf.gz"
    s_eq1=$(n_sites "${tmp}/poly.eq1.vcf.gz")
    bcftools concat -a "${tmp}/poly.ge2.imputed.vcf.gz" "${tmp}/poly.eq1.vcf.gz" -Ou \
      | bcftools sort -T "${tmp}" -Oz -o "${vcf_final}"
  else
    bcftools sort -T "${tmp}" -Oz -o "${vcf_final}" "${tmp}/poly.ge2.imputed.vcf.gz"
  fi
  bcftools index -c --threads "${THREADS}" "${vcf_final}"

  echo -e "${cohort}\t$(n_sites "${vcf_in}")\t$(n_sites "${vcf_poly}")\t$(wc -l < "${tmp}/contigs_ge2.txt")\t${n_eq1}\t${s_eq1}\t$(n_sites "${vcf_final}")" >> "${COUNTS}"
  echo "[$(date)] ${COH} done"
done

cat "${COUNTS}"
echo "[$(date)] Step 6a done"
