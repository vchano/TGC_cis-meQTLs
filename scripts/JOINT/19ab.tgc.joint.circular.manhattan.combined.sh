#!/bin/bash
#SBATCH -p scc-cpu
#SBATCH -t 02:00:00
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=48G
#SBATCH --job-name=TGC.manhattan.v4
#SBATCH --output=/path/to/your/project/LOGS/%x_%j.out
#SBATCH --error=/path/to/your/project/LOGS/%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=YOUR_EMAIL

set -euo pipefail

module purge
module load gcc/14.2.0
module load r/4.5.2
module load imagemagick/7.1.1-39

export LC_ALL=C.UTF-8
export LANG=C.UTF-8
export R_LIBS_USER="/path/to/your/Rlibs"
mkdir -p "${R_LIBS_USER}"

PROJECT_ROOT="/path/to/your/project"
SCRIPTS="${PROJECT_ROOT}/SCRIPTS/JOINT"
mkdir -p "${PROJECT_ROOT}/LOGS"

echo "============================================================"
echo "TGC — Circular Manhattan v4 (24cm, wider rings, tighter margins)"
echo "Node:  $(hostname)"
echo "Start: $(date)"
echo "============================================================"

echo "--- Step 16ab4.R: individual _4 circular plots ---"
Rscript --vanilla "${SCRIPTS}/16ab4.R"
echo "16ab4.R done — $(date)"

echo "--- Step 19ab_manhattan_combined_v4.R: combined _4 panels ---"
Rscript --vanilla "${SCRIPTS}/19ab_manhattan_combined_v4.R"
echo "19ab_manhattan_combined_v4.R done — $(date)"

echo "============================================================"
echo "All done: $(date)"
echo "Outputs (v4): ${PROJECT_ROOT}/RESULTS/CORRECTED/FIGURES/NEW/"
echo "============================================================"
