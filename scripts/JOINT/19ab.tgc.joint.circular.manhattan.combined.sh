#!/bin/bash
#SBATCH -p YOUR_PARTITION
#SBATCH -t 06:00:00
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=96G
#SBATCH --job-name=TGC.manhattan.v5
#SBATCH --output=/path/to/your/project/LOGS/%x_%j.out
#SBATCH --error=/path/to/your/project/LOGS/%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=YOUR_EMAIL

# 19ab.tgc.joint.circular.manhattan.combined.sh
#
# Runs in two steps:
#   Step 1: 16ab.tgc.joint.panel.assembly.R
#           Renders 8 individual circular plots (2 tools x 2 cohorts x 2 axis modes)
#           as TIFF + PDF + SVG + EPS. Also writes legend_circular_5.tiff (36 cm).
#
#   Step 2: 19ab.tgc.joint.circular.manhattan.combined.R
#           Assembles 4 combined panels:
#             Figure3_circular_panel_genesis5_5        -- GENESIS    p-value  A(Breeding)+B(Natural)
#             EDF4_circular_panel_matrixeqtl5_5        -- MatrixEQTL p-value  A(Breeding)+B(Natural)
#             Figure3_circular_panel_genesis5_5_FDR    -- GENESIS    FDR      A(Breeding)+B(Natural)
#             EDF4_circular_panel_matrixeqtl5_5_FDR    -- MatrixEQTL FDR      A(Breeding)+B(Natural)
#           Each panel in: tiff, pdf, svg, png, eps
#           Copies all formats to RESULTS/CORRECTED/FIGURES/NEW/
#
# v5 rendering specs vs v4:
#   Canvas: 18 cm (was 24)        CHR_LABEL_CEX: 1.20 (was 0.975)
#   Dot cex: 0.22 (was 0.15)     Margins: c(0.5,2,1.0,2) (was c(2,2,2.5,2))
#   Legend height: 2.0 cm (was 3.5)
#   Y-axis tick/title cex: 0.55/0.65 (unchanged -- larger caused overlap)
#
# NOTE: FDR_AXIS=TRUE on NATURAL/CHH (~500M tests) requires ~96 GB RAM.
#
# USAGE
#   sbatch 19ab.tgc.joint.circular.manhattan.combined.sh

set -euo pipefail

module purge
module load gcc/14.2.0
module load r/4.5.2
module load imagemagick/7.1.1-39

export LC_ALL=C.UTF-8
export LANG=C.UTF-8
export R_LIBS_USER="/path/to/your/Rlibs"   # <-- set this
mkdir -p "${R_LIBS_USER}"

# === USER CONFIGURATION ===
PROJECT_ROOT="/path/to/your/project"        # <-- set this
# ===========================
SCRIPTS="${PROJECT_ROOT}/SCRIPTS/JOINT"

mkdir -p "${PROJECT_ROOT}/LOGS"

echo "============================================================"
echo "TGC -- Circular Manhattan v5 (18 cm, FDR variant, 4 panels)"
echo "Node:  $(hostname)"
echo "Start: $(date)"
echo "RAM:   $(free -h | awk '/^Mem:/{print $2}')"
echo "============================================================"

# ── Step 1: individual circular plots (8 TIFFs + other formats + legends) ────
echo "--- 16ab.tgc.joint.panel.assembly.R: 8 individual plots ---"
Rscript --vanilla "${SCRIPTS}/16ab.tgc.joint.panel.assembly.R"
echo "16ab done -- $(date)"

# Verify key individual TIFFs exist before assembly
echo "--- Verifying individual TIFFs ---"
NEEDED=(
  "manhattan_circular_genesis5_breeding_5.tiff"
  "manhattan_circular_genesis5_natural_5.tiff"
  "manhattan_circular_matrixeqtl5_breeding_5.tiff"
  "manhattan_circular_matrixeqtl5_natural_5.tiff"
  "manhattan_circular_genesis5_breeding_5fdr.tiff"
  "manhattan_circular_genesis5_natural_5fdr.tiff"
  "manhattan_circular_matrixeqtl5_breeding_5fdr.tiff"
  "manhattan_circular_matrixeqtl5_natural_5fdr.tiff"
  "legend_circular_5.tiff"
)
MAN_DIR="${PROJECT_ROOT}/RESULTS/JOINT/COMBINED5/manhattan"
ALL_OK=1
for f in "${NEEDED[@]}"; do
  if [[ ! -f "${MAN_DIR}/${f}" ]]; then
    echo "  MISSING: ${f}" >&2
    ALL_OK=0
  else
    echo "  OK: ${f}  ($(stat -c%s "${MAN_DIR}/${f}") bytes)"
  fi
done
if [[ "${ALL_OK}" -eq 0 ]]; then
  echo "ERROR: one or more individual TIFFs missing -- aborting assembly" >&2
  exit 1
fi

# ── Step 2: assemble 4 combined panels (all formats) ─────────────────────────
echo "--- 19ab.tgc.joint.circular.manhattan.combined.R: 4 combined panels ---"
Rscript --vanilla "${SCRIPTS}/19ab.tgc.joint.circular.manhattan.combined.R"
echo "19ab done -- $(date)"

# ── Summary ──────────────────────────────────────────────────────────────────
echo "============================================================"
echo "All done: $(date)"
echo ""
echo "Panels:"
PANEL_DIR="${PROJECT_ROOT}/RESULTS/JOINT/COMBINED5/panels"
for stem in \
  "Figure3_circular_panel_genesis5_5" \
  "EDF4_circular_panel_matrixeqtl5_5" \
  "Figure3_circular_panel_genesis5_5_FDR" \
  "EDF4_circular_panel_matrixeqtl5_5_FDR"; do
  for fmt in tiff pdf svg png eps; do
    f="${PANEL_DIR}/${stem}.${fmt}"
    if [[ -f "${f}" ]]; then
      echo "  ${stem}.${fmt}  ($(stat -c%s "${f}") bytes)"
    else
      echo "  MISSING: ${stem}.${fmt}"
    fi
  done
done
echo ""
echo "Copies in: ${PROJECT_ROOT}/RESULTS/CORRECTED/FIGURES/NEW/"
echo "============================================================"
