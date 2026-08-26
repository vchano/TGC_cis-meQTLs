#!/bin/bash
#-------------------------------------------------------------------------------
# TGC — Compute per-sample ECS read counts from trimmed FASTQ files
# Outputs: RESULTS/JOINT/COMBINED5/tables/ecs_read_counts.tsv
#          (sample, trimmed_reads, mapped_reads, dup_rate)
#-------------------------------------------------------------------------------

#SBATCH -p scc-cpu
#SBATCH -t 02:00:00
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=8G
#SBATCH --job-name=TGC.ecs.readcounts
#SBATCH --output=/path/to/your/project/LOGS/%x_%j.out
#SBATCH --error=/path/to/your/project/LOGS/%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=YOUR_EMAIL

set -euo pipefail

module purge
module load gcc/14.2.0
module load samtools/1.21

export LC_ALL=C.UTF-8

PROJECT_ROOT="/path/to/your/project"
TRIM_DIR="${PROJECT_ROOT}/DATA/ECS/TRIMMED.FASTQ.ECS"
BAM_DIR="${PROJECT_ROOT}/DATA/ECS/MAPPED.FILES.ECS"
OUT_DIR="${PROJECT_ROOT}/RESULTS/JOINT/COMBINED5/tables"
mkdir -p "${OUT_DIR}"
OUT_TSV="${OUT_DIR}/ecs_read_counts.tsv"

echo "sample	trimmed_reads	mapped_reads	dup_reads	dedup_mapped_reads" > "${OUT_TSV}"

echo "[$(date)] Starting per-sample read count computation..."

for r1 in "${TRIM_DIR}"/*_R1_p.fastq.gz; do
  sample=$(basename "${r1}" _R1_p.fastq.gz)

  # Trimmed read count: lines/4 from R1 file (paired reads → 1 count per pair)
  trim_reads=$(zcat "${r1}" | wc -l | awk '{print $1/4}')

  # Mapped reads from dedup BAM flagstat
  dedup_bam="${BAM_DIR}/${sample}_bt2.fix.psrt.dedup.bam"
  if [[ -f "${dedup_bam}" ]]; then
    flagstat=$(samtools flagstat -@ 2 "${dedup_bam}")
    total_mapped=$(echo "${flagstat}" | grep "^[0-9]* + [0-9]* mapped" | awk '{print $1}')
    total_reads=$(echo "${flagstat}" | grep "^[0-9]* + [0-9]* in total" | awk '{print $1}')
    dup_reads=$(echo "${flagstat}" | grep "^[0-9]* + [0-9]* duplicates" | awk '{print $1}')
    dedup_mapped=$((total_mapped - dup_reads))
  else
    total_reads="NA"
    dup_reads="NA"
    dedup_mapped="NA"
  fi

  echo "${sample}	${trim_reads}	${total_reads}	${dup_reads}	${dedup_mapped}" >> "${OUT_TSV}"
  echo "[$(date)]   ${sample}: trimmed=${trim_reads}"
done

echo "[$(date)] Done. Output: ${OUT_TSV}"
echo "Samples counted: $(wc -l < "${OUT_TSV}") (incl. header)"
