#!/bin/bash
############################################################
# TreeGeneClimate (TGC) — JOINT ECS + TMS
# Step 16ab8: cis-meQTL validation — SLURM wrappers
#
# Three stages, each submitted as a separate SLURM job.
# Run each stage after the previous one completes.
#
# STAGE 1 — LD deconvolution (requires PLINK in PATH):
#   sbatch --job-name=TGC.ld.clump 16ab8.tgc.joint.meqtl.validation.sh ld
#
# STAGE 2 — Whole-methylome Mantel test:
#   sbatch --job-name=TGC.mantel 16ab8.tgc.joint.meqtl.validation.sh mantel
#
# STAGE 3a — Permutation FDR (sequential, one cohort at a time):
#   sbatch --export=ALL,TGC_COHORT=BREEDING,TGC_CONTEXT=CpG,PERM_START=1,PERM_END=100 \
#          --job-name=TGC.perm.B 16ab8.tgc.joint.meqtl.validation.sh perm
#   sbatch --export=ALL,TGC_COHORT=NATURAL,TGC_CONTEXT=CpG,PERM_START=1,PERM_END=100 \
#          --job-name=TGC.perm.N 16ab8.tgc.joint.meqtl.validation.sh perm
#
# STAGE 3b — Aggregate permutation results and generate figures:
#   sbatch --job-name=TGC.perm.collect 16ab8.tgc.joint.meqtl.validation.sh perm_collect
############################################################

set -euo pipefail

# === USER CONFIGURATION ===
PROJECT_ROOT="/path/to/your/project"   # <-- set this
R_LIBS_USER="/path/to/your/Rlibs"     # <-- set this
PARTITION="YOUR_PARTITION"
ACCOUNT="YOUR_ACCOUNT"
EMAIL="YOUR_EMAIL"
# ===========================

SCRIPTS="${PROJECT_ROOT}/SCRIPTS/JOINT"
LOGS="${PROJECT_ROOT}/LOGS"
mkdir -p "${LOGS}"

MODE="${1:-help}"

module_setup() {
  module purge
  module load gcc/14.2.0
  module load r/4.5.2
  export R_LIBS_USER="${R_LIBS_USER}"
  export LC_ALL=C.UTF-8
}

# ---- Stage 1: LD deconvolution ----
if [ "${MODE}" = "ld" ]; then
  #SBATCH -p YOUR_PARTITION
  #SBATCH -t 02:00:00
  #SBATCH -N 1 -c 1
  #SBATCH --mem=8G
  #SBATCH --mail-type=END,FAIL
  #SBATCH --mail-user=YOUR_EMAIL
  #SBATCH --output=YOUR_LOGS_DIR/ld_clump_%j.out
  #SBATCH --error=YOUR_LOGS_DIR/ld_clump_%j.err
  module_setup
  module load plink/1.9   # adjust to your HPC module name
  echo "TGC Step 16ab8 — LD clumping"
  echo "Start: $(date)"
  Rscript --vanilla "${SCRIPTS}/16ab8.tgc.joint.meqtl.validation.R" --section ld
  echo "Done: $(date)"
  exit 0
fi

# ---- Stage 2: Mantel test ----
if [ "${MODE}" = "mantel" ]; then
  #SBATCH -p YOUR_PARTITION
  #SBATCH -t 04:00:00
  #SBATCH -N 1 -c 1
  #SBATCH --mem=24G
  #SBATCH --mail-type=END,FAIL
  #SBATCH --mail-user=YOUR_EMAIL
  #SBATCH --output=YOUR_LOGS_DIR/mantel_%j.out
  #SBATCH --error=YOUR_LOGS_DIR/mantel_%j.err
  module_setup
  echo "TGC Step 16ab8 — Mantel test"
  echo "Start: $(date)"
  Rscript --vanilla "${SCRIPTS}/16ab8.tgc.joint.meqtl.validation.R" --section mantel
  echo "Done: $(date)"
  exit 0
fi

# ---- Stage 3a: Permutation FDR (sequential loop) ----
# Usage: sbatch --export=ALL,TGC_COHORT=BREEDING,TGC_CONTEXT=CpG,PERM_START=1,PERM_END=100 ...
if [ "${MODE}" = "perm" ]; then
  #SBATCH -p YOUR_PARTITION
  #SBATCH -t 08:00:00
  #SBATCH -N 1 -c 4
  #SBATCH --mem=32G
  #SBATCH --mail-type=END,FAIL
  #SBATCH --mail-user=YOUR_EMAIL
  #SBATCH --output=YOUR_LOGS_DIR/perm_seq_%j.out
  #SBATCH --error=YOUR_LOGS_DIR/perm_seq_%j.err
  : "${TGC_COHORT:?TGC_COHORT not set}"
  : "${TGC_CONTEXT:?TGC_CONTEXT not set}"
  : "${PERM_START:=1}"
  : "${PERM_END:=100}"
  module_setup
  echo "TGC Step 16ab8 — Permutation FDR"
  echo "Cohort=${TGC_COHORT}  Context=${TGC_CONTEXT}  Range=${PERM_START}-${PERM_END}"
  echo "Start: $(date)"
  for i in $(seq "${PERM_START}" "${PERM_END}"); do
    export TGC_PERM_ID="${i}"
    echo "--- Perm ${i}: $(date) ---"
    Rscript --vanilla "${SCRIPTS}/16ab8.tgc.joint.meqtl.validation.R" --section perm
  done
  echo "All perms done: $(date)"
  exit 0
fi

# ---- Stage 3b: Aggregate permutation results ----
if [ "${MODE}" = "perm_collect" ]; then
  #SBATCH -p YOUR_PARTITION
  #SBATCH -t 00:30:00
  #SBATCH -N 1 -c 1
  #SBATCH --mem=8G
  #SBATCH --mail-type=END,FAIL
  #SBATCH --mail-user=YOUR_EMAIL
  #SBATCH --output=YOUR_LOGS_DIR/perm_collect_%j.out
  #SBATCH --error=YOUR_LOGS_DIR/perm_collect_%j.err
  module_setup
  echo "TGC Step 16ab8 — Permutation collect"
  echo "Start: $(date)"
  Rscript --vanilla "${SCRIPTS}/16ab8.tgc.joint.meqtl.validation.R" --section perm_collect
  echo "Done: $(date)"
  exit 0
fi

echo "Usage: $0 {ld|mantel|perm|perm_collect}"
echo "  ld           — LD deconvolution of hotspot SNP counts"
echo "  mantel       — Whole-methylome Mantel test"
echo "  perm         — Permutation FDR worker (set TGC_COHORT, TGC_CONTEXT, PERM_START, PERM_END)"
echo "  perm_collect — Aggregate permutation results and generate figures"
exit 1
