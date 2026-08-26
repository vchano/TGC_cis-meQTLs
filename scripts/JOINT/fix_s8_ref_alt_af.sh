#!/bin/bash
#-------------------------------------------------------------------------------
# Fix S8: populate missing ref/alt/af in Supplementary Table 6
# Indexes VCF files, queries missing values, writes corrected xlsx
#
# USAGE
#   sbatch fix_s8_ref_alt_af.sh
#-------------------------------------------------------------------------------

#SBATCH -p scc-cpu
#SBATCH -t 01:30:00
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH --job-name=TGC.fix.s8
#SBATCH --output=/path/to/your/project/LOGS/%x_%j.out
#SBATCH --error=/path/to/your/project/LOGS/%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=YOUR_EMAIL

set -euo pipefail

module purge
module load gcc/14.2.0
module load r/4.5.2
module load bcftools/1.19

export LC_ALL=C.UTF-8
export LANG=C.UTF-8
export R_LIBS_USER="/path/to/your/Rlibs"
mkdir -p "${R_LIBS_USER}"

PROJECT_ROOT="/path/to/your/project"
SCRIPTS="${PROJECT_ROOT}/SCRIPTS/JOINT"
mkdir -p "${PROJECT_ROOT}/LOGS"

echo "============================================================"
echo "TGC — Fix S8: populate missing ref/alt/af"
echo "Node:  $(hostname)"
echo "Start: $(date)"
echo "============================================================"

Rscript --vanilla "${SCRIPTS}/fix_s8_ref_alt_af.R"

echo "============================================================"
echo "Done: $(date)"
echo "Output: ${PROJECT_ROOT}/RESULTS/DRAFT/NATURE.GENETICS/260819_Chano.etal.2026_tgc_supp.tables_corrected_s8.xlsx"
echo "============================================================"
