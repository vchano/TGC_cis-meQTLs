#!/bin/bash
#SBATCH --job-name=JOINT.16ab.H2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --output=LOGS/16ab/%x_%A_%a.out
#SBATCH --error=LOGS/16ab/%x_%A_%a.err
#-------------------------------------------------------------------------------
# TreeGeneClimate (TGC) — JOINT
# Step 16ab: per-site heritability in chunks of 2,000 sites, then collection
#
#   bash SCRIPTS/JOINT/16ab.tgc.joint.meth.heritability.sh submit
#       computes the number of chunks per cohort x context from the Step 12ab0
#       summary and submits (i) the chunk array and (ii) the collect job
#       (summary, Mantel test, heritability panel) after the array has finished.
# Create LOGS/16ab before submitting.
#-------------------------------------------------------------------------------
set -euo pipefail

# === USER CONFIGURATION ===
export TGC_PROJECT_ROOT="${TGC_PROJECT_ROOT:-/path/to/your/project}"
CHUNK_SIZE=2000          # must match CHUNK_SIZE in the R script
EXTRA_SBATCH="${EXTRA_SBATCH:-}"   # e.g. "--partition=xyz"
COLLECT_AFTER="${COLLECT_AFTER:-}" # optional job id(s) the collect step also waits for (Step 13ab)
module purge 2>/dev/null || true
module load gcc/14.2.0 r/4.5.2 2>/dev/null || true
export LC_ALL=C.UTF-8 LANG=C.UTF-8
# ===========================

SUMMARY="${TGC_PROJECT_ROOT}/RESULTS/JOINT/MEQTL/INPUTS/step12ab0_summary.tsv"
# panels and number of chunks: "COHORT CONTEXT NCHUNKS"
panels () {
  awk -F'\t' -v cs="${CHUNK_SIZE}" 'NR==1{for(i=1;i<=NF;i++) h[$i]=i; next}
    {n=$h["n_sites_tested"]; print $h["cohort"], $h["context"], int((n + cs - 1) / cs)}' "${SUMMARY}"
}

if [[ "${1:-}" == "submit" ]]; then
  total=$(panels | awk '{s+=$3} END{print s}')
  jid=$(sbatch --parsable ${EXTRA_SBATCH} --array=0-$((total - 1))%100 \
        --export=ALL,TGC_MODE=chunk "$0")
  sbatch ${EXTRA_SBATCH} --dependency=afterok:"${jid}${COLLECT_AFTER:+:${COLLECT_AFTER}}" --job-name=JOINT.16ab.COLLECT \
        --mem=32G --cpus-per-task=2 --time=06:00:00 --array=0 \
        --export=ALL,TGC_MODE=collect "$0"
  echo "Submitted ${total} chunk tasks (job ${jid}) and the collect job"
  exit 0
fi

if [[ "${TGC_MODE:-chunk}" == "chunk" ]]; then
  i="${SLURM_ARRAY_TASK_ID:-0}"
  read -r TGC_COHORT TGC_CONTEXT TGC_CHUNK < <(panels | awk -v i="${i}" \
    '{ if (i < $3) { print $1, $2, i + 1; exit } i -= $3 }')
  export TGC_COHORT TGC_CONTEXT TGC_CHUNK
fi
Rscript --vanilla "${TGC_PROJECT_ROOT}/SCRIPTS/JOINT/16ab.tgc.joint.meth.heritability.R"
