#!/bin/bash
#SBATCH --job-name=ECS.4a3.MASK
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=LOGS/%x_%j.out
#SBATCH --error=LOGS/%x_%j.err
#-------------------------------------------------------------------------------
# TreeGeneClimate (TGC) — ECS
# Step 4a3: SNP table used to mask methylation calls
#
# In enzymatic methyl-seq, an unmethylated cytosine is read as T. A C>T SNP on
# either strand (REF/ALT C/T on the forward strand, G/A for cytosines on the
# reverse strand) therefore looks like loss of methylation, and SNPs next to a
# cytosine change its sequence context (CpG/CHG/CHH). Every SNP called with
# reasonable confidence in the joint ECS call set is listed here, with its
# allele counts. Steps 5b2 and 5b3 keep the SNPs whose alternative allele was
# called at least twice (any call rate) and remove methylation sites with such
# a SNP at the cytosine or at the bases defining its context on the same strand
# (CpG: C and G; CHG/CHH: C and the two following bases).
#
# The table is built from the joint call set before MAF filtering, so that rare
# variants are also masked (filters for genotype analyses and for masking have
# different purposes).
#
# Selection: SNP records (bi- or multi-allelic), QUAL >= 20, INFO/DP >= 10,
#            at least one alternative allele called (INFO/AC > 0).
#
# Input (Step 4a2):
#   RESULTS/ECS/VARIANT.CALLING/tgc.ecs.allsamples.call.unfilt.snvs.renamed.vcf.gz
#
# Output:
#   RESULTS/ECS/VARIANT.CALLING/MASK/tgc.ecs.cytosine.mask.snps.tsv.gz
#     columns: chr pos ref alt qual ac an
#   RESULTS/ECS/VARIANT.CALLING/MASK/mask_counts.tsv
#
# Software: BCFtools 1.19
#-------------------------------------------------------------------------------
set -euo pipefail

# === USER CONFIGURATION ===
PROJECT_ROOT="${TGC_PROJECT_ROOT:-/path/to/your/project}"
THREADS="${SLURM_CPUS_PER_TASK:-4}"
MIN_QUAL=20
MIN_SITE_DP=10
# Software environment (edit for your system)
module purge 2>/dev/null || true
module load gcc/14.2.0 bcftools/1.19 2>/dev/null || true
# ===========================

VCF_DIR="${PROJECT_ROOT}/RESULTS/ECS/VARIANT.CALLING"
IN_VCF="${VCF_DIR}/tgc.ecs.allsamples.call.unfilt.snvs.renamed.vcf.gz"
OUT_DIR="${VCF_DIR}/MASK"
OUT="${OUT_DIR}/tgc.ecs.cytosine.mask.snps.tsv.gz"
mkdir -p "${OUT_DIR}"
[[ -s "${IN_VCF}" ]] || { echo "ERROR: missing input ${IN_VCF}"; exit 1; }
echo "[$(date)] Step 4a3 started"

bcftools view --threads "${THREADS}" --types snps \
    --include "QUAL>=${MIN_QUAL} && INFO/DP>=${MIN_SITE_DP} && INFO/AC>0" -Ou "${IN_VCF}" \
| bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%INFO/AC\t%INFO/AN\n' \
| awk 'BEGIN{OFS="\t"; print "chr","pos","ref","alt","qual","ac","an"} {print}' \
| gzip -c > "${OUT}"

{
  echo -e "class\tn_snps"
  zcat "${OUT}" | awk -F'\t' 'NR>1 {
      n++
      a = $3"/"$4
      if (a=="C/T" || a=="G/A") ct++
      else if (a=="T/C" || a=="A/G") tc++
      else if ($4 ~ /,/) multi++
      else other++
    } END {
      print "all\t" n; print "C/T_or_G/A\t" ct; print "T/C_or_A/G\t" tc
      print "multiallelic\t" multi; print "other_biallelic\t" other
    }'
} > "${OUT_DIR}/mask_counts.tsv"
cat "${OUT_DIR}/mask_counts.tsv"

echo "[$(date)] Step 4a3 done"
