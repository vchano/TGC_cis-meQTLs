#!/bin/bash
#SBATCH -p scc-cpu
#SBATCH -t 06:00:00
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=96G
#SBATCH --job-name=TGC.14ab
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=vmchano.gaug@gmail.com

# 14ab.tgc.joint.manhattan.plots.sh
#
# Runs step 14ab in a single pass:
#   Part 1 — 8 individual circular Manhattan plots
#             (GENESIS5 + MatrixEQTL5) x (BREEDING + NATURAL) x (raw-p + BH-FDR)
#             Formats: TIFF, PDF, SVG, EPS + standalone legend TIFFs
#   Part 2 — 4 combined two-panel figures (BREEDING | NATURAL + legend strip)
#             Figure3_circular_panel_genesis5_5
#             EDF4_circular_panel_matrixeqtl5_5
#             Figure3_circular_panel_genesis5_5_FDR
#             EDF4_circular_panel_matrixeqtl5_5_FDR
#             Formats: TIFF (ImageMagick), PDF, SVG, PNG, EPS
#             Copies all to RESULTS/CORRECTED/FIGURES/NEW/
#
# Memory note: FDR_AXIS=TRUE on NATURAL/CHH (~500 M tests) requires ~96 GB RAM.

set -euo pipefail

module purge
module load gcc/14.2.0
module load r/4.5.2
module load imagemagick/7.1.1-39

export LC_ALL=C.UTF-8
export LANG=C.UTF-8
# export R_LIBS_USER="/path/to/your/Rlibs"  # uncomment and set if needed

PROJECT_ROOT="/path/to/your/project"  # <-- set this
SCRIPT="${PROJECT_ROOT}/SCRIPTS/JOINT/14ab.tgc.joint.manhattan.plots.R"

mkdir -p "${PROJECT_ROOT}/LOGS"
export TGC_PROJECT_ROOT="${PROJECT_ROOT}"

echo "============================================================"
echo "TGC — Circular Manhattan plots + Combined panels (step 14ab)"
echo "Node:  $(hostname)"
echo "Start: $(date)"
echo "RAM:   $(free -h | awk '/^Mem:/{print $2}')"
echo "============================================================"

Rscript --vanilla "${SCRIPT}"

echo "============================================================"
echo "All done: $(date)"

PANEL_DIR="${PROJECT_ROOT}/RESULTS/JOINT/COMBINED5/panels"
echo ""
echo "Combined panels:"
for stem in \
  "Figure3_circular_panel_genesis5_5" \
  "EDF4_circular_panel_matrixeqtl5_5" \
  "Figure3_circular_panel_genesis5_5_FDR" \
  "EDF4_circular_panel_matrixeqtl5_5_FDR"; do
  for fmt in tiff pdf svg png eps; do
    f="${PANEL_DIR}/${stem}.${fmt}"
    if [[ -f "${f}" ]]; then
      echo "  OK: ${stem}.${fmt}  ($(stat -c%s "${f}") bytes)"
    else
      echo "  MISSING: ${stem}.${fmt}"
    fi
  done
done

echo ""
echo "Copies in: ${PROJECT_ROOT}/RESULTS/CORRECTED/FIGURES/NEW/"
echo "============================================================"
