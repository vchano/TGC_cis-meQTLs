#!/bin/bash
############################################################
# TreeGeneClimate (TGC) — JOINT ECS + TMS
# Step 16ab7: Methylation heritability — SLURM submission
#
# For small datasets: run the R script directly.
#   Rscript 16ab7.tgc.joint.meth.heritability.R
#
# For large datasets (>100k sites per context): distribute
# across cluster nodes using the array below, then collect.
#
# STAGE 1 — submit array (one task per chunk × cohort × context):
#   bash 16ab7.tgc.joint.meth.heritability.sh
#
# STAGE 2 — collect results after all array tasks finish:
#   sbatch --job-name=TGC.h2.collect \
#          --wrap="Rscript 16ab7.tgc.joint.meth.heritability.R" \
#          16ab7.tgc.joint.meth.heritability.sh collect
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
OUTDIR="${PROJECT_ROOT}/RESULTS/JOINT/HERITABILITY"
LOGS="${PROJECT_ROOT}/LOGS"
mkdir -p "${OUTDIR}" "${LOGS}"

CHUNK_SIZE=1000
MODE="${1:-submit}"

# ---- Collect mode: aggregate chunks and generate figures ----
if [ "${MODE}" = "collect" ]; then
  #SBATCH -p YOUR_PARTITION
  #SBATCH -t 01:00:00
  #SBATCH -N 1 -c 2
  #SBATCH --mem=32G
  #SBATCH --job-name=TGC.h2.collect
  #SBATCH --output=YOUR_LOGS_DIR/h2_collect_%j.out
  #SBATCH --error=YOUR_LOGS_DIR/h2_collect_%j.err
  #SBATCH --mail-type=END,FAIL
  #SBATCH --mail-user=YOUR_EMAIL
  module purge
  module load r/4.5.2
  export R_LIBS_USER="${R_LIBS_USER}"
  export LC_ALL=C.UTF-8
  Rscript --vanilla "${SCRIPTS}/16ab7.tgc.joint.meth.heritability.R"
  exit 0
fi

# ---- Submit mode: launch one SLURM array per cohort × context ----
# Site counts (adjust to match your dataset)
declare -A N_SITES
N_SITES["BREEDING_CpG"]=36662
N_SITES["BREEDING_CHG"]=89636
N_SITES["BREEDING_CHH"]=273738
N_SITES["NATURAL_CpG"]=66872
N_SITES["NATURAL_CHG"]=156957
N_SITES["NATURAL_CHH"]=499222

echo "TGC Step 16ab7 — heritability array submission"
echo "Chunk size: ${CHUNK_SIZE} sites"

for COHORT in BREEDING NATURAL; do
  for CONTEXT in CpG CHG CHH; do
    KEY="${COHORT}_${CONTEXT}"
    N="${N_SITES[${KEY}]}"
    N_CHUNKS=$(( (N + CHUNK_SIZE - 1) / CHUNK_SIZE ))
    OUTDIR_CTX="${OUTDIR}/${COHORT}/${CONTEXT}"
    mkdir -p "${OUTDIR_CTX}"
    echo "${COHORT}/${CONTEXT}: ${N} sites → ${N_CHUNKS} chunks"

    JOBID=$(sbatch --parsable \
      -p "${PARTITION}" \
      -t 00:30:00 -N 1 -c 1 --mem=8G \
      --job-name="TGC.h2.${COHORT:0:1}.${CONTEXT}" \
      --output="${LOGS}/h2_${COHORT}_${CONTEXT}_%A_%a.out" \
      --error="${LOGS}/h2_${COHORT}_${CONTEXT}_%A_%a.err" \
      --mail-type=END,FAIL --mail-user="${EMAIL}" \
      --array="1-${N_CHUNKS}%50" \
      --wrap="
        set -euo pipefail
        module purge && module load r/4.5.2
        export R_LIBS_USER='${R_LIBS_USER}'
        export LC_ALL=C.UTF-8
        CHUNK_START=\$(( (SLURM_ARRAY_TASK_ID - 1) * ${CHUNK_SIZE} + 1 ))
        CHUNK_END=\$(( SLURM_ARRAY_TASK_ID * ${CHUNK_SIZE} ))
        OUT_FILE='${OUTDIR_CTX}/chunk_'\${SLURM_ARRAY_TASK_ID}'.tsv'
        Rscript --vanilla '${SCRIPTS}/16ab7.tgc.joint.meth.heritability.R' \
          ${COHORT} ${CONTEXT} \${CHUNK_START} \${CHUNK_END} \${OUT_FILE}
      ")
    echo "  Submitted: job ${JOBID} (${N_CHUNKS} tasks)"
  done
done

echo ""
echo "When all tasks complete, run:"
echo "  bash $0 collect"
