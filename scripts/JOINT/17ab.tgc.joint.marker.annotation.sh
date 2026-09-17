#!/bin/bash
#-------------------------------------------------------------------------------
# TreeGeneClimate (TGC) — JOINT ECS + TMS
# Step 17ab: Marker annotation against reference genome GFF3
#            (Table 2, Table 3, Supplementary Table S3, S6, S7)
#
# Annotates two sets of markers with gene features from the GFF3:
#   1. TMS SVMPs — breeding top-150 balanced + natural 6 formal SVMPs (step 8b)
#   2. meQTL     — significant SNP + methylation site positions
#                  (p_FDR < 1e-10, both tools, both cohorts, all contexts)
#
# For each marker the script reports:
#   - Overlapping gene (distance_bp = 0) or nearest gene with distance
#   - annotation_class: genic | proximal_intergenic | distal_intergenic
#   - REF, ALT, AF, DR2 from the imputed VCF (SNP markers only)
#
# Section 9 additionally back-fills missing REF/ALT/AF in the master
# supplementary-tables workbook (Supplementary Table S6) from raw/unfiltered
# VCFs, for positions absent from the primary filtered VCF.
#
# REQUIRES (adjust module names for your HPC environment):
#   bedtools >= 2.27   (loaded below — CHECK MODULE NAME)
#   bcftools >= 1.19   (confirmed: module bcftools/1.19)
#   R >= 4.5           (confirmed: module r/4.5.2)
#
# USAGE
#   sbatch 17ab.sh
#
# INPUTS
#   RESULTS/TMS/RANALYSIS/TABLES/heatmap_markers_8B/  (needs 8b.R first)
#   RESULTS/JOINT/COMBINED5/overlap/tables/robust_markers_*.tsv  (needs 15ab.R first)
#   RESULTS/ECS/VCF_SPLIT/tgc.ecs.*.imputed.vcf.gz
#   REFERENCE/Pabies2.0/Picab02_230926_at01_all_sorted.gff3
#
# OUTPUTS (RESULTS/JOINT/ANNOTATION17/)
#   tms_svmp_{breeding,natural}_annotated.tsv   -> Supplementary Table S3
#   robust_{breeding,natural}_{snps,sites}_annotated.tsv  -> Supp. Table S6, S7
#   all_markers_annotated.tsv
#   gene_summary_annotated.tsv                  -> Table 2 (source genes), Table 3 (target genes)
#-------------------------------------------------------------------------------

#SBATCH -p YOUR_PARTITION
#SBATCH -t 02:00:00
#SBATCH -N 1
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH --job-name=TGC.annot17
#SBATCH --output=/path/to/your/project/LOGS/%x_%j.out
#SBATCH --error=/path/to/your/project/LOGS/%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=YOUR_EMAIL

set -euo pipefail

# ---------------------------------------------------------------------------
# Environment
# ---------------------------------------------------------------------------
module purge
module load gcc/14.2.0
module load r/4.5.2
module load bcftools/1.19

# bedtools2/2.31.1 confirmed available (gcc/14.2.0 must be loaded first)
module load bedtools2/2.31.1

export LC_ALL=C.UTF-8
export LANG=C.UTF-8
# export R_LIBS_USER="/path/to/your/Rlibs"  # uncomment and set if needed

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
# === USER CONFIGURATION ===
PROJECT_ROOT="/path/to/your/project"  # <-- set this
# ===========================
SCRIPTS="${PROJECT_ROOT}/SCRIPTS/JOINT"

mkdir -p "${PROJECT_ROOT}/LOGS"

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
echo "============================================================"
echo "TGC — JOINT — Step 17ab — Marker GFF3 annotation"
echo "Node:   $(hostname)"
echo "CPUs:   ${SLURM_CPUS_ON_NODE:-2}"
echo "Memory: ${SLURM_MEM_PER_NODE:-?} MB"
echo "Start:  $(date)"
echo "============================================================"

# Verify bedtools and bcftools are reachable
if ! command -v bedtools &>/dev/null; then
  echo "ERROR: bedtools not found in PATH. Adjust the module name above." >&2
  exit 1
fi
if ! command -v bcftools &>/dev/null; then
  echo "ERROR: bcftools not found in PATH." >&2
  exit 1
fi

echo "bedtools: $(bedtools --version | head -1)"
echo "bcftools: $(bcftools --version | head -1)"

Rscript --vanilla "${SCRIPTS}/17ab.R"

echo "============================================================"
echo "Step 17ab finished: $(date)"
echo "============================================================"
