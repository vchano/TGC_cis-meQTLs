# Code Archive — Cis-meQTL mapping reveals genetic architecture of DNA methylation in a giga-genome species

**Victor Chano¹·², Konstantin V. Krutovsky¹·², Kai Wang³, Matti Haapanen⁴, Fred O. Asiegbu³, Oliver Gailing¹·²**

¹ Department of Forest Genetics and Forest Tree Breeding, University of Göttingen. Büsgenweg 2, 37077 Göttingen, Germany.
² Center for Integrated Breeding Research (CiBreed), University of Göttingen. Albrecht-Thaer-Weg 3, 37075 Göttingen, Germany.
³ Department of Forest Science, University of Helsinki. FIN-00014 Helsinki, Finland.
⁴ Natural Resources Institute Finland (LUKE). Latokartanonkaari 9, Helsinki 00790, Finland.

\* Corresponding authors: Victor Chano, Oliver Gailing

**ORCIDs**
Victor Chano: 0000-0002-3423-9090
Konstantin V. Krutovsky: 0000-0002-8819-7084
Kai Wang: 0000-0002-5205-7346
Matti Haapanen: 0000-0003-3294-501X
Fred O. Asiegbu: 0000-0003-0223-7194
Oliver Gailing: 0000-0002-4572-2408

---

**GitHub repository:** https://github.com/vchano/TGC_cis-meQTLs

---

**Use of AI coding assistance**
Claude Code (Anthropic; claude.ai/code) was used to assist with the development of pipeline scripts and R code. All code suggestions were reviewed, tested, and validated by the first author prior to inclusion.

---

## Description

This archive contains the data-processing and analysis scripts of the study, from sequencing reads to cis-meQTL maps, heritability estimates and marker annotation. Scripts are organised into three sequential modules:

- **ECS** (Exome Capture Sequencing): raw data QC, adapter trimming, alignment, SNP calling and filtering, the SNP mask for methylation data, cohort splitting and LD pruning, imputation, IBD estimation, genotype files, GRM and PCA (Steps 1a–10a).
- **TMS** (Targeted Enzymatic Methylation Sequencing): raw data QC, trimming, alignment, PCR-duplicate removal, methylation extraction, methylation import and filtering (population-level and cis-meQTL datasets), methylation levels, epigenomic PCA and SVMP heatmaps (Steps 1b–8b).
- **JOINT**: genomic–epigenomic concordance, cis-meQTL inputs and mapping (MatrixEQTL and GENESIS), combined results and robust pairs, circular Manhattan plots, Venn overlaps, methylation heritability and whole-methylome Mantel test, and marker annotation (Steps 11ab–17ab).

R scripts are run with `Rscript --vanilla <script>` (R ≥ 4.5.2) from the project root, with `TGC_PROJECT_ROOT` set to it (or `PROJECT_ROOT` at the top of each script); the `.sh` scripts are SLURM job scripts. Results are written as tab-separated tables; each figure is written with a table of the values it shows (`*_values.tsv`).

**R packages:** data.table, methylKit, ggplot2, patchwork, ComplexHeatmap, circlize, SNPRelate, gdsfmt, GENESIS, GWASTools, MatrixEQTL, lme4, vegan, ade4, dunn.test, multcompView, ggvenn, dplyr, tibble
**Bioinformatics tools:** FastQC, MultiQC, Trimmomatic, Bowtie2, SAMtools, BCFtools, PLINK 1.9, BEAGLE 5, Bismark, BEDTools, ImageMagick

---

## ECS — Exome Capture Sequencing

---

### `1a.tgc.ecs.fastqc.rawdata.sh`

```bash
#!/bin/bash
#-------------------------------------------------------------------------------
# TreeGeneClimate (TGC) — ECS
# Step 1a: FASTQC + MULTIQC on RAWDATA (ECS)
#
# Project root:
#   /path/to/your/project
#
# Input (data, frozen):
#   DATA/ECS/RAWDATA.ECS/*fastq.gz
#
# Output (results):
#   RESULTS/ECS/QC/RAWDATA/FASTQC/
#   RESULTS/ECS/QC/RAWDATA/MULTIQC/
#-------------------------------------------------------------------------------

#SBATCH -n 48
#SBATCH -N 1
#SBATCH --job-name=ECS.FQC1
#SBATCH --output=LOGS/%x_%j.out
#SBATCH --error=LOGS/%x_%j.err
#SBATCH --ntasks-per-socket 24
#SBATCH --time=48:00:00

set -euo pipefail

echo "[$(date)] SLURM job started: ${SLURM_JOB_NAME:-no_slurm}"

module purge
module load fastqc/0.11.4
module load multiqc/1.27.1

# === USER CONFIGURATION ===
PROJECT_ROOT="${TGC_PROJECT_ROOT:-/path/to/your/project}"
# ===========================

INPUT="${PROJECT_ROOT}/DATA/ECS/RAWDATA.ECS"
QC_BASE="${PROJECT_ROOT}/RESULTS/ECS/QC/RAWDATA"
QC_FASTQC="${QC_BASE}/FASTQC"
QC_MULTIQC="${QC_BASE}/MULTIQC"
LOGS="${PROJECT_ROOT}/LOGS"

mkdir -p "${QC_FASTQC}" "${QC_MULTIQC}" "${LOGS}"

shopt -s nullglob
raw_fastq=( "${INPUT}"/*.fastq.gz )
if (( ${#raw_fastq[@]} == 0 )); then
  echo "ERROR: no FASTQ found in: ${INPUT}"
  exit 1
fi

fastqc "${raw_fastq[@]}" --outdir "${QC_FASTQC}" --threads 48

source activate multiqc
multiqc "${QC_FASTQC}" -o "${QC_MULTIQC}"
conda deactivate

echo "[$(date)] Done."
exit 0
```

---

### `2a.tgc.ecs.trimmomatic.and.fastqc.trimmed.sh`

```bash
#!/bin/bash
#-------------------------------------------------------------------------------
# TreeGeneClimate (TGC) — ECS
# Step 2a: TRIMMOMATIC (paired-end) + FASTQC + MULTIQC on TRIMMED FASTQ (ECS)
#
# Project root:
#   /path/to/your/project
#
# Input (data, frozen):
#   DATA/ECS/RAWDATA.ECS/*_R1.fastq.gz
#   DATA/ECS/RAWDATA.ECS/*_R2.fastq.gz
#
# Output (data, frozen):
#   DATA/ECS/TRIMMED.FASTQ.ECS/           (paired reads)
#   DATA/ECS/TRIMMED.FASTQ.ECS/UNPAIRED/  (unpaired reads)
#
# Output (results):
#   RESULTS/ECS/QC/TRIMMED/FASTQC/
#   RESULTS/ECS/QC/TRIMMED/MULTIQC/
#
# Adapters:
#   DATA/METADATA/adapters.fa
#-------------------------------------------------------------------------------

#SBATCH -n 48
#SBATCH -N 1
#SBATCH --job-name=ECS.TRIM
#SBATCH --output=LOGS/%x_%j.out
#SBATCH --error=LOGS/%x_%j.err
#SBATCH --ntasks-per-socket 24
#SBATCH --time=48:00:00

set -euo pipefail

echo "[$(date)] SLURM job started: ${SLURM_JOB_NAME:-no_slurm}"

module purge
module load fastqc/0.11.4
module load anaconda3/2020.11
module load trimmomatic/0.36

# === USER CONFIGURATION ===
PROJECT_ROOT="${TGC_PROJECT_ROOT:-/path/to/your/project}"
# ===========================

INPUT="${PROJECT_ROOT}/DATA/ECS/RAWDATA.ECS"
TRIMMED="${PROJECT_ROOT}/DATA/ECS/TRIMMED.FASTQ.ECS"
UNPAIRED="${TRIMMED}/UNPAIRED"

QC_BASE="${PROJECT_ROOT}/RESULTS/ECS/QC/TRIMMED"
QC_FASTQC="${QC_BASE}/FASTQC"
QC_MULTIQC="${QC_BASE}/MULTIQC"

ADAPTERS="${PROJECT_ROOT}/DATA/METADATA/adapters.fa"
LOGS="${PROJECT_ROOT}/LOGS"

mkdir -p "${TRIMMED}" "${UNPAIRED}" "${QC_FASTQC}" "${QC_MULTIQC}" "${LOGS}"

if [[ ! -s "${ADAPTERS}" ]]; then
  echo "ERROR: adapters file not found or empty: ${ADAPTERS}"
  exit 1
fi

# ------------------------------------------------------------------------------
# Sample list (DATA/METADATA)
list=$(cat "${PROJECT_ROOT}/DATA/METADATA/samples_final_602.txt")

for sample in ${list}; do
  r1="${INPUT}/${sample}_R1.fastq.gz"
  r2="${INPUT}/${sample}_R2.fastq.gz"

  if [[ ! -s "${r1}" || ! -s "${r2}" ]]; then
    echo "WARNING: missing pair for ${sample}; skipping"
    continue
  fi

  java -jar "${TRIMMOMATIC_JAR:-/path/to/trimmomatic-0.36.jar}" PE \
    -threads 48 -phred33 \
    "${r1}" "${r2}" \
    "${TRIMMED}/${sample}_R1_p.fastq.gz" "${UNPAIRED}/${sample}_R1_u.fastq.gz" \
    "${TRIMMED}/${sample}_R2_p.fastq.gz" "${UNPAIRED}/${sample}_R2_u.fastq.gz" \
    ILLUMINACLIP:"${ADAPTERS}":2:30:10 \
    CROP:138 HEADCROP:12 SLIDINGWINDOW:5:20 MINLEN:30
done

shopt -s nullglob
trim_fastq=( "${TRIMMED}"/*_p.fastq.gz )
if (( ${#trim_fastq[@]} == 0 )); then
  echo "ERROR: no trimmed paired FASTQ produced in: ${TRIMMED}"
  exit 1
fi

fastqc "${trim_fastq[@]}" --outdir "${QC_FASTQC}" --threads 48

source activate multiqc
multiqc "${QC_FASTQC}" -o "${QC_MULTIQC}"
conda deactivate

echo "[$(date)] Done."
exit 0
```

---

### `3a.tgc.ecs.bowtie.trimmed.sh`

```bash
#!/bin/bash
#-------------------------------------------------------------------------------
# TreeGeneClimate (TGC) — ECS
# Step 3a: BOWTIE2 mapping + SAMTOOLS processing on TRIMMED FASTQ (ECS)
#
# Project root:
#   /path/to/your/project
#
# Input (data, frozen):
#   DATA/ECS/TRIMMED.FASTQ.ECS/*_R1_p.fastq.gz
#   DATA/ECS/TRIMMED.FASTQ.ECS/*_R2_p.fastq.gz
#
# Reference (frozen):
#   REFERENCE/Pabies2.0/   (bowtie2 index prefix: Pabies2.0)
#
# Output (data, frozen):
#   DATA/ECS/MAPPED.FILES.ECS/
#-------------------------------------------------------------------------------

#SBATCH -n 96
#SBATCH -N 1
#SBATCH --job-name=ECS.BT2
#SBATCH --output=LOGS/%x_%j.out
#SBATCH --error=LOGS/%x_%j.err
#SBATCH --ntasks-per-socket 24
#SBATCH --time=48:00:00

set -euo pipefail

echo "[$(date)] SLURM job started: ${SLURM_JOB_NAME:-no_slurm}"

module purge
module load gcc/14.2.0
module load bowtie2/2.5.4
module load samtools/1.21

# === USER CONFIGURATION ===
PROJECT_ROOT="${TGC_PROJECT_ROOT:-/path/to/your/project}"
# ===========================

INPUT="${PROJECT_ROOT}/DATA/ECS/TRIMMED.FASTQ.ECS"
OUTPUT="${PROJECT_ROOT}/DATA/ECS/MAPPED.FILES.ECS"

REF_BASE="${PROJECT_ROOT}/REFERENCE/Pabies2.0"
REF_PREFIX="${REF_BASE}/Pabies2.0"

TMP="${OUTPUT}/tmp"
OTHER="${OUTPUT}/OTHER_FILES"
LOGS="${PROJECT_ROOT}/LOGS"

mkdir -p "${OUTPUT}" "${TMP}" "${OTHER}" "${LOGS}"

# ------------------------------------------------------------------------------
# Sample list (DATA/METADATA)
list=$(cat "${PROJECT_ROOT}/DATA/METADATA/samples_final_602.txt")

for sample in ${list}; do
  r1="${INPUT}/${sample}_R1_p.fastq.gz"
  r2="${INPUT}/${sample}_R2_p.fastq.gz"

  if [[ ! -s "${r1}" || ! -s "${r2}" ]]; then
    echo "WARNING: missing pair for ${sample}; skipping"
    continue
  fi

  bowtie2 -p 96 -N 1 -L 25 -i S,1,2.00 --local \
    --un-gz "${OTHER}/${sample}.un.fastq.gz" \
    --al-gz "${OTHER}/${sample}.al.fastq.gz" \
    --un-conc-gz "${OTHER}/${sample}.un.conc.fastq.gz" \
    --al-conc-gz "${OTHER}/${sample}.al.conc.fastq.gz" \
    -x "${REF_PREFIX}" \
    -1 "${r1}" -2 "${r2}" \
    -S "${OUTPUT}/${sample}_bt2.sam"

  samtools view -S -b "${OUTPUT}/${sample}_bt2.sam" > "${OUTPUT}/${sample}_bt2.bam"

  samtools sort -l 1 -@ 96 -n \
    -o "${OUTPUT}/${sample}_bt2.nsrt.bam" \
    -T "${TMP}/${sample}_tmp" \
    "${OUTPUT}/${sample}_bt2.bam"

  samtools fixmate -@ 96 -O bam,level=1 -m \
    "${OUTPUT}/${sample}_bt2.nsrt.bam" \
    "${OUTPUT}/${sample}_bt2.nsrt.fix.bam"

  samtools sort -l 1 -@ 96 \
    -o "${OUTPUT}/${sample}_bt2.fix.psrt.bam" \
    -T "${TMP}/${sample}_tmp" \
    "${OUTPUT}/${sample}_bt2.nsrt.fix.bam"

  samtools markdup -@ 96 -r -s -O bam,level=1 \
    "${OUTPUT}/${sample}_bt2.fix.psrt.bam" \
    "${OUTPUT}/${sample}_bt2.fix.psrt.dedup.bam"

  rm -f "${OUTPUT}/${sample}_bt2.sam"
done

echo "[$(date)] Done."
exit 0
```

---

### `4a1.tgc.ecs.snv.calling.bcftools.sh`

```bash
#!/bin/bash
#-------------------------------------------------------------------------------
# TreeGeneClimate (TGC) — ECS
# Step 4a1: SNV calling (bcftools mpileup + call) from deduplicated BAMs
#
# Project root:
#   /path/to/your/project
#
# Input (data, frozen):
#   DATA/ECS/MAPPED.FILES.ECS/*_bt2.fix.psrt.dedup.bam
#
# Reference (frozen):
#   REFERENCE/Pabies2.0/Picab02_chromosomes_and_unplaced.fa
#
# Output (results):
#   RESULTS/ECS/VARIANT.CALLING/regions/chr_*.vcf.gz
#   RESULTS/ECS/VARIANT.CALLING/tgc.ecs.allsamples.call.unfilt.snvs.vcf.gz
#   RESULTS/ECS/VARIANT.CALLING/bam.list
#   RESULTS/ECS/VARIANT.CALLING/chrom.list
#-------------------------------------------------------------------------------

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --job-name=ECS.CALL
#SBATCH --output=LOGS/%x_%j.out
#SBATCH --error=LOGS/%x_%j.err
#SBATCH --time=48:00:00

set -euo pipefail

echo "[$(date)] SLURM job started: ${SLURM_JOB_NAME:-no_slurm}"

module purge
module load gcc/14.2.0
module load bcftools/1.19
module load samtools/1.21

# === USER CONFIGURATION ===
PROJECT_ROOT="${TGC_PROJECT_ROOT:-/path/to/your/project}"
# ===========================

BAM_DIR="${PROJECT_ROOT}/DATA/ECS/MAPPED.FILES.ECS"
REF_FASTA="${PROJECT_ROOT}/REFERENCE/Pabies2.0/Picab02_chromosomes_and_unplaced.fa"

VCF_BASE="${PROJECT_ROOT}/RESULTS/ECS/VARIANT.CALLING"
REGIONS_DIR="${VCF_BASE}/regions"
LOGS_DIR="${PROJECT_ROOT}/LOGS"

MERGED_VCF="${VCF_BASE}/tgc.ecs.allsamples.call.unfilt.snvs.vcf.gz"
BAM_LIST="${VCF_BASE}/bam.list"
CHROM_LIST="${VCF_BASE}/chrom.list"
FAIL_LOG="${VCF_BASE}/calling.fail.log"

mkdir -p "${VCF_BASE}" "${REGIONS_DIR}" "${LOGS_DIR}"
: > "${FAIL_LOG}"

# ------------------------------------------------------------------------------
# CPU layout (full-node, no oversubscription)
CPUS="${SLURM_CPUS_PER_TASK:-192}"
MAX_JOBS=24
THREADS_PER_JOB=$(( CPUS / MAX_JOBS ))   # 192/24 = 8
if [[ "${THREADS_PER_JOB}" -lt 1 ]]; then THREADS_PER_JOB=1; fi

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

TMPDIR_BASE="${SLURM_TMPDIR:-/tmp}"
TMP="${TMPDIR_BASE}/${SLURM_JOB_ID}"
mkdir -p "${TMP}"

echo "[$(date)] CPUs per task: ${CPUS}"
echo "[$(date)] Parallel jobs: ${MAX_JOBS}"
echo "[$(date)] Threads per job (bgzip): ${THREADS_PER_JOB}"
echo "[$(date)] TMP: ${TMP}"

# ------------------------------------------------------------------------------
# Reference index
if [[ ! -s "${REF_FASTA}.fai" ]]; then
  samtools faidx "${REF_FASTA}"
fi
cut -f1 "${REF_FASTA}.fai" > "${CHROM_LIST}"

# ------------------------------------------------------------------------------
# BAM list
find "${BAM_DIR}" -maxdepth 1 -type f -name '*_bt2.fix.psrt.dedup.bam' | sort > "${BAM_LIST}"

if [[ ! -s "${BAM_LIST}" ]]; then
  echo "ERROR: No BAMs found in: ${BAM_DIR}"
  exit 1
fi

BAM_COUNT="$(wc -l < "${BAM_LIST}")"
echo "[$(date)] BAMs found: ${BAM_COUNT}"

# ------------------------------------------------------------------------------
# Parallel per-chromosome calling
JOB_COUNT=0

while read -r chr; do
  (
    echo "[$(date)] Calling: ${chr}"

    out_vcf="${REGIONS_DIR}/chr_${chr}.vcf.gz"

    bcftools mpileup -Ou \
      -f "${REF_FASTA}" \
      -r "${chr}" \
      -b "${BAM_LIST}" \
      -a AD,DP,SP \
      2> "${REGIONS_DIR}/${chr}_mpileup.err" \
    | bcftools call -mv -f GQ,GP \
      --threads "${THREADS_PER_JOB}" \
      -Oz -o "${out_vcf}" \
      2> "${REGIONS_DIR}/${chr}_call.err"

    if [[ ! -s "${out_vcf}" ]]; then
      echo "FAILED: ${chr}" >> "${FAIL_LOG}"
      exit 1
    fi

    # CSI index (safe for large contigs/positions)
    bcftools index -c --threads "${THREADS_PER_JOB}" "${out_vcf}" \
      2> "${REGIONS_DIR}/${chr}_index.err" || true
  ) &

  JOB_COUNT=$((JOB_COUNT + 1))
  if [[ "${JOB_COUNT}" -ge "${MAX_JOBS}" ]]; then
    wait
    JOB_COUNT=0
  fi
done < "${CHROM_LIST}"

wait
echo "[$(date)] Per-chromosome calling complete."

# ------------------------------------------------------------------------------
# Merge per-chromosome VCFs (ordering from chrom.list)
VCF_FILES=()
while read -r chr; do
  f="${REGIONS_DIR}/chr_${chr}.vcf.gz"
  if [[ -s "${f}" ]]; then
    VCF_FILES+=("${f}")
  else
    echo "MISSING_VCF: ${chr}" >> "${FAIL_LOG}"
  fi
done < "${CHROM_LIST}"

if [[ "${#VCF_FILES[@]}" -eq 0 ]]; then
  echo "ERROR: No per-chrom VCFs were produced."
  exit 1
fi

bcftools concat --threads "${CPUS}" -Oz -o "${MERGED_VCF}" "${VCF_FILES[@]}"

# CSI index for merged VCF
bcftools index -c --threads "${CPUS}" "${MERGED_VCF}" || true

echo "[$(date)] Unfiltered merged VCF: ${MERGED_VCF}"
echo "[$(date)] Done."
exit 0
```

---

### `4a2.tgc.ecs.snv.filtering.bcftools.sh`

```bash
#!/bin/bash
#SBATCH --job-name=ECS.4a2.FILT
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=LOGS/%x_%j.out
#SBATCH --error=LOGS/%x_%j.err
#-------------------------------------------------------------------------------
# TreeGeneClimate (TGC) — ECS
# Step 4a2: Sample renaming and SNP filtering of the joint call set
#
#   (i)   Rename samples in the merged VCF from BAM paths to sample IDs
#         (e.g. .../P001_WA02_bt2.fix.psrt.dedup.bam -> P001_WA02)
#   (ii)  Keep biallelic SNPs with
#           F_MISSING <= 0.20, MAF >= 0.05, INFO/DP >= 10 (summed site depth)
#         computed over all 602 samples (both cohorts jointly)
#   (iii) Site counts at each stage
#
# Submit from the project root (sbatch SCRIPTS/ECS/4a2...), with
# TGC_PROJECT_ROOT set or PROJECT_ROOT edited below.
#
# Input (Step 4a1):
#   RESULTS/ECS/VARIANT.CALLING/tgc.ecs.allsamples.call.unfilt.snvs.vcf.gz
#
# Output:
#   RESULTS/ECS/VARIANT.CALLING/tgc.ecs.allsamples.call.unfilt.snvs.renamed.vcf.gz
#   RESULTS/ECS/VARIANT.CALLING/tgc.ecs.allsamples.call.filt.maf05.snvs.vcf.gz
#   RESULTS/ECS/VARIANT.CALLING/STATS/
#     unfiltered.bcftools.stats.txt, filtered_maf05.bcftools.stats.txt
#     ecs_filtering_counts.tsv
#
# Software: BCFtools 1.19
#-------------------------------------------------------------------------------
set -euo pipefail

# === USER CONFIGURATION ===
PROJECT_ROOT="${TGC_PROJECT_ROOT:-/path/to/your/project}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"
# Software environment (edit for your system)
module purge 2>/dev/null || true
module load gcc/14.2.0 bcftools/1.19 2>/dev/null || true
# ===========================

VCF_DIR="${PROJECT_ROOT}/RESULTS/ECS/VARIANT.CALLING"
STATS_DIR="${VCF_DIR}/STATS"
RAW_VCF="${VCF_DIR}/tgc.ecs.allsamples.call.unfilt.snvs.vcf.gz"
REN_VCF="${VCF_DIR}/tgc.ecs.allsamples.call.unfilt.snvs.renamed.vcf.gz"
FILT_VCF="${VCF_DIR}/tgc.ecs.allsamples.call.filt.maf05.snvs.vcf.gz"
COUNTS="${STATS_DIR}/ecs_filtering_counts.tsv"

mkdir -p "${STATS_DIR}"
[[ -s "${RAW_VCF}" ]] || { echo "ERROR: missing input ${RAW_VCF}"; exit 1; }
echo "[$(date)] Step 4a2 started (threads=${THREADS})"

# ------------------------------------------------------------------------------
# (i) Rename samples: strip directory and BAM suffix
# ------------------------------------------------------------------------------
NAMES="${STATS_DIR}/sample_rename.tsv"
bcftools query -l "${RAW_VCF}" \
  | awk 'BEGIN{OFS="\t"} {n=$0; sub(/^.*\//, "", n); sub(/_bt2\.fix\.psrt\.dedup\.bam$/, "", n); print $0, n}' \
  > "${NAMES}"
bcftools reheader --samples "${NAMES}" --threads "${THREADS}" -o "${REN_VCF}" "${RAW_VCF}"
bcftools index -c --threads "${THREADS}" "${REN_VCF}"
echo "[$(date)] Renamed $(wc -l < "${NAMES}") samples"

# ------------------------------------------------------------------------------
# (ii) Joint filtering
# ------------------------------------------------------------------------------
bcftools +fill-tags "${REN_VCF}" -Ou -- -t MAF,F_MISSING,AF \
| bcftools view \
    --types snps --min-alleles 2 --max-alleles 2 \
    --include 'F_MISSING<=0.2 && MAF>=0.05 && INFO/DP>=10' \
    --threads "${THREADS}" -Oz -o "${FILT_VCF}"
bcftools index -c --threads "${THREADS}" "${FILT_VCF}"
echo "[$(date)] Filtered VCF written"

# ------------------------------------------------------------------------------
# (iii) Counts
# ------------------------------------------------------------------------------
bcftools stats --threads "${THREADS}" "${REN_VCF}"  > "${STATS_DIR}/unfiltered.bcftools.stats.txt"
bcftools stats --threads "${THREADS}" "${FILT_VCF}" > "${STATS_DIR}/filtered_maf05.bcftools.stats.txt"

sn () { awk -F'\t' -v k="$2" '$1=="SN" && $3==k {print $4}' "$1"; }
{
  echo -e "stage\trecords\tsnps\tmultiallelic_snp_sites"
  for s in unfiltered filtered_maf05; do
    f="${STATS_DIR}/${s}.bcftools.stats.txt"
    echo -e "${s}\t$(sn "$f" 'number of records:')\t$(sn "$f" 'number of SNPs:')\t$(sn "$f" 'number of multiallelic SNP sites:')"
  done
} > "${COUNTS}"
cat "${COUNTS}"

echo "[$(date)] Step 4a2 done"
```

---

### `4a3.tgc.ecs.cytosine.snp.mask.sh`

```bash
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
```

---

### `5a.tgc.ecs.plink.split.pruning.sh`

```bash
#!/bin/bash
#SBATCH --job-name=ECS.5a.SPLIT
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=LOGS/%x_%j.out
#SBATCH --error=LOGS/%x_%j.err
#-------------------------------------------------------------------------------
# TreeGeneClimate (TGC) — ECS
# Step 5a: Split the filtered joint VCF into BREEDING and NATURAL cohorts, then
#          (i)   unpruned PLINK BED set per cohort
#          (ii)  LD pruning (PLINK --indep-pairwise 50 10 0.2)
#          (iii) PCA on pruned SNPs (PLINK)
#          (iv)  pruned PLINK BED set per cohort (input to Step 7a)
#
# Input (Step 4a2):
#   RESULTS/ECS/VARIANT.CALLING/tgc.ecs.allsamples.call.filt.maf05.snvs.vcf.gz
# Cohort sample lists:
#   DATA/METADATA/all_breeding_samples.txt
#   DATA/METADATA/all_natural_samples.txt
#
# Output:
#   RESULTS/ECS/VCF_SPLIT/tgc.ecs.<cohort>.call.filt.maf05.snvs.vcf.gz
#   RESULTS/ECS/POPGEN/STRUCTURE/<COHORT>/
#     tgc.ecs.<cohort>.unpruned.{bed,bim,fam}
#     tgc.ecs.<cohort>.prune.prune.{in,out}
#     tgc.ecs.<cohort>.pruned.{bed,bim,fam}
#     PCA/tgc.ecs.<cohort>.pca.pruned.{eigenvec,eigenval}
#
# Software: BCFtools 1.19, PLINK 1.9
#-------------------------------------------------------------------------------
set -euo pipefail

# === USER CONFIGURATION ===
PROJECT_ROOT="${TGC_PROJECT_ROOT:-/path/to/your/project}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"
PRUNE_WINDOW=50   # SNPs
PRUNE_STEP=10     # SNPs
PRUNE_R2=0.2
# Software environment (edit for your system)
module purge 2>/dev/null || true
module load gcc/14.2.0 bcftools/1.19 plink/1.9 2>/dev/null || true
# ===========================

VCF_ALL="${PROJECT_ROOT}/RESULTS/ECS/VARIANT.CALLING/tgc.ecs.allsamples.call.filt.maf05.snvs.vcf.gz"
META_DIR="${PROJECT_ROOT}/DATA/METADATA"
SPLIT_DIR="${PROJECT_ROOT}/RESULTS/ECS/VCF_SPLIT"
STRUCT_DIR="${PROJECT_ROOT}/RESULTS/ECS/POPGEN/STRUCTURE"

[[ -s "${VCF_ALL}" ]] || { echo "ERROR: missing input ${VCF_ALL}"; exit 1; }
mkdir -p "${SPLIT_DIR}"
echo "[$(date)] Step 5a started (threads=${THREADS})"

plink_vcf () {  # common PLINK options for VCF input
  plink --vcf "$1" --double-id --allow-extra-chr --set-missing-var-ids @:# --threads "${THREADS}" "${@:2}"
}

for cohort in breeding natural; do
  COH=$(echo "${cohort}" | tr '[:lower:]' '[:upper:]')
  OUT="${STRUCT_DIR}/${COH}"
  mkdir -p "${OUT}/PCA"
  samples="${META_DIR}/all_${cohort}_samples.txt"
  [[ -s "${samples}" ]] || { echo "ERROR: missing ${samples}"; exit 1; }
  vcf="${SPLIT_DIR}/tgc.ecs.${cohort}.call.filt.maf05.snvs.vcf.gz"

  echo "[$(date)] ${COH}: split"
  bcftools view --threads "${THREADS}" -S "${samples}" -Oz -o "${vcf}" "${VCF_ALL}"
  bcftools index -c --threads "${THREADS}" "${vcf}"

  echo "[$(date)] ${COH}: unpruned BED set"
  plink_vcf "${vcf}" --make-bed --out "${OUT}/tgc.ecs.${cohort}.unpruned"

  echo "[$(date)] ${COH}: LD pruning (${PRUNE_WINDOW} ${PRUNE_STEP} ${PRUNE_R2})"
  plink_vcf "${vcf}" --indep-pairwise "${PRUNE_WINDOW}" "${PRUNE_STEP}" "${PRUNE_R2}" \
    --out "${OUT}/tgc.ecs.${cohort}.prune"

  echo "[$(date)] ${COH}: pruned BED set and PCA"
  plink_vcf "${vcf}" --extract "${OUT}/tgc.ecs.${cohort}.prune.prune.in" \
    --make-bed --out "${OUT}/tgc.ecs.${cohort}.pruned"
  plink_vcf "${vcf}" --extract "${OUT}/tgc.ecs.${cohort}.prune.prune.in" \
    --pca --out "${OUT}/PCA/tgc.ecs.${cohort}.pca.pruned"

  echo "[$(date)] ${COH}: $(wc -l < "${OUT}/tgc.ecs.${cohort}.prune.prune.in") SNPs after pruning"
done

echo "[$(date)] Step 5a done"
```

---

### `6a.tgc.ecs.beagle.imputation.sh`

```bash
#!/bin/bash
#SBATCH --job-name=ECS.6a.BEAGLE
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=LOGS/%x_%j.out
#SBATCH --error=LOGS/%x_%j.err
#-------------------------------------------------------------------------------
# TreeGeneClimate (TGC) — ECS
# Step 6a: Within-cohort polymorphic sites and BEAGLE imputation, per cohort
#   (i)   keep sites polymorphic within the cohort (AC > 0 && AC < AN)
#   (ii)  impute and phase with BEAGLE (contigs with >= 2 variants)
#   (iii) append variants on single-variant contigs unchanged (BEAGLE cannot
#         process contigs with one variant) and sort
#
# BEAGLE results depend on the number of threads (default seed); nthreads is fixed
# at 48 so that the imputed genotypes are reproducible.
#
# Input (Step 5a):
#   RESULTS/ECS/VCF_SPLIT/tgc.ecs.<cohort>.call.filt.maf05.snvs.vcf.gz
#
# Output:
#   RESULTS/ECS/VCF_SPLIT/tgc.ecs.<cohort>.call.filt.maf05.snvs.poly.vcf.gz
#   RESULTS/ECS/VCF_SPLIT/tgc.ecs.<cohort>.call.filt.maf05.snvs.poly.imputed.vcf.gz
#   RESULTS/ECS/VCF_SPLIT/imputation_counts.tsv
#   RESULTS/ECS/VCF_SPLIT/BEAGLE_TMP/<COHORT>/   (intermediates)
#
# Software: BCFtools 1.19, BEAGLE 5 (conda environment "beagle")
#-------------------------------------------------------------------------------
set -euo pipefail

# === USER CONFIGURATION ===
PROJECT_ROOT="${TGC_PROJECT_ROOT:-/path/to/your/project}"
THREADS="${SLURM_CPUS_PER_TASK:-16}"   # bcftools
BEAGLE_THREADS=48                      # fixed: results depend on the number of threads
export _JAVA_OPTIONS="-Xmx100g"
# Software environment (edit for your system)
module purge 2>/dev/null || true
module load gcc/14.2.0 bcftools/1.19 miniforge3/24.3.0-0 2>/dev/null || true
source "$(conda info --base)/etc/profile.d/conda.sh" 2>/dev/null || true
conda activate beagle 2>/dev/null || true
# ===========================

SPLIT_DIR="${PROJECT_ROOT}/RESULTS/ECS/VCF_SPLIT"
TMP_BASE="${SPLIT_DIR}/BEAGLE_TMP"
COUNTS="${SPLIT_DIR}/imputation_counts.tsv"
mkdir -p "${TMP_BASE}"

die () { echo "[$(date)] ERROR: $*" >&2; exit 1; }
command -v beagle >/dev/null || die "beagle not found in PATH"
n_sites () { bcftools index -n "$1"; }

echo "[$(date)] Step 6a started (BEAGLE threads=${BEAGLE_THREADS}); $(beagle 2>&1 | head -n 1)"
echo -e "cohort\tsites_split\tsites_poly\tcontigs_imputed\tcontigs_single_variant\tsites_single_variant\tsites_final" > "${COUNTS}"

for cohort in breeding natural; do
  COH=$(echo "${cohort}" | tr '[:lower:]' '[:upper:]')
  prefix="${SPLIT_DIR}/tgc.ecs.${cohort}.call.filt.maf05.snvs"
  vcf_in="${prefix}.vcf.gz"; vcf_poly="${prefix}.poly.vcf.gz"; vcf_final="${prefix}.poly.imputed.vcf.gz"
  tmp="${TMP_BASE}/${COH}"; mkdir -p "${tmp}"
  [[ -s "${vcf_in}" ]] || die "missing ${vcf_in}"

  # (i) polymorphic within cohort (AC/AN were recomputed for the cohort in Step 5a)
  bcftools view --threads "${THREADS}" --include 'AC>0 && AC<AN' -Oz -o "${vcf_poly}" "${vcf_in}"
  bcftools index -c --threads "${THREADS}" "${vcf_poly}"

  # (ii) contigs with >= 2 variants go to BEAGLE; single-variant contigs are kept as they are
  bcftools query -f '%CHROM\n' "${vcf_poly}" | sort | uniq -c \
    | awk -v ge2="${tmp}/contigs_ge2.txt" -v eq1="${tmp}/contigs_eq1.txt" \
        '$1>=2 {print $2 > ge2} $1==1 {print $2 > eq1}'
  touch "${tmp}/contigs_ge2.txt" "${tmp}/contigs_eq1.txt"
  bcftools view --threads "${THREADS}" -r "$(paste -sd, "${tmp}/contigs_ge2.txt")" \
    -Oz -o "${tmp}/poly.ge2.vcf.gz" "${vcf_poly}"
  bcftools index -c "${tmp}/poly.ge2.vcf.gz"

  beagle gt="${tmp}/poly.ge2.vcf.gz" out="${tmp}/poly.ge2.imputed" nthreads="${BEAGLE_THREADS}"
  [[ -s "${tmp}/poly.ge2.imputed.vcf.gz" ]] || die "${COH}: BEAGLE produced no output"
  bcftools index -c "${tmp}/poly.ge2.imputed.vcf.gz"

  # (iii) append single-variant contigs and sort
  n_eq1=$(wc -l < "${tmp}/contigs_eq1.txt"); s_eq1=0
  if [[ "${n_eq1}" -gt 0 ]]; then
    bcftools view --threads "${THREADS}" -r "$(paste -sd, "${tmp}/contigs_eq1.txt")" \
      -Oz -o "${tmp}/poly.eq1.vcf.gz" "${vcf_poly}"
    bcftools index -c "${tmp}/poly.eq1.vcf.gz"
    s_eq1=$(n_sites "${tmp}/poly.eq1.vcf.gz")
    bcftools concat -a "${tmp}/poly.ge2.imputed.vcf.gz" "${tmp}/poly.eq1.vcf.gz" -Ou \
      | bcftools sort -T "${tmp}" -Oz -o "${vcf_final}"
  else
    bcftools sort -T "${tmp}" -Oz -o "${vcf_final}" "${tmp}/poly.ge2.imputed.vcf.gz"
  fi
  bcftools index -c --threads "${THREADS}" "${vcf_final}"

  echo -e "${cohort}\t$(n_sites "${vcf_in}")\t$(n_sites "${vcf_poly}")\t$(wc -l < "${tmp}/contigs_ge2.txt")\t${n_eq1}\t${s_eq1}\t$(n_sites "${vcf_final}")" >> "${COUNTS}"
  echo "[$(date)] ${COH} done"
done

cat "${COUNTS}"
echo "[$(date)] Step 6a done"
```

---

### `7a.tgc.ecs.plink.ibd.sh`

```bash
#!/bin/bash
#SBATCH --job-name=ECS.7a.IBD
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=LOGS/%x_%j.out
#SBATCH --error=LOGS/%x_%j.err
#-------------------------------------------------------------------------------
# TreeGeneClimate (TGC) — ECS
# Step 7a: Pairwise identity by descent (PLINK --genome full) per cohort,
#          on the LD-pruned SNP sets
#
# Input (Step 5a):
#   RESULTS/ECS/POPGEN/STRUCTURE/<COHORT>/tgc.ecs.<cohort>.pruned.{bed,bim,fam}
#
# Output:
#   RESULTS/ECS/POPGEN/RELATEDNESS/IBD/tgc.ecs.<cohort>.pruned.ibd.genome
#
# Software: PLINK 1.9
#-------------------------------------------------------------------------------
set -euo pipefail

# === USER CONFIGURATION ===
PROJECT_ROOT="${TGC_PROJECT_ROOT:-/path/to/your/project}"
THREADS="${SLURM_CPUS_PER_TASK:-4}"
# Software environment (edit for your system)
module purge 2>/dev/null || true
module load gcc/14.2.0 plink/1.9 2>/dev/null || true
# ===========================

STRUCT_DIR="${PROJECT_ROOT}/RESULTS/ECS/POPGEN/STRUCTURE"
OUT_DIR="${PROJECT_ROOT}/RESULTS/ECS/POPGEN/RELATEDNESS/IBD"
mkdir -p "${OUT_DIR}"
echo "[$(date)] Step 7a started; $(plink --version 2>&1 | head -n 1)"

for cohort in breeding natural; do
  COH=$(echo "${cohort}" | tr '[:lower:]' '[:upper:]')
  bfile="${STRUCT_DIR}/${COH}/tgc.ecs.${cohort}.pruned"
  for ext in bed bim fam; do [[ -s "${bfile}.${ext}" ]] || { echo "ERROR: missing ${bfile}.${ext}"; exit 1; }; done
  out="${OUT_DIR}/tgc.ecs.${cohort}.pruned.ibd"
  plink --bfile "${bfile}" --allow-extra-chr --threads "${THREADS}" --genome full --out "${out}"
  echo "[$(date)] ${COH}: $(($(wc -l < "${out}.genome") - 1)) pairs"
done

echo "[$(date)] Step 7a done"
```

---

### `8a.tgc.ecs.gds.grm.ibd.R`

```r
#!/usr/bin/env Rscript
############################################################
# TreeGeneClimate (TGC) — ECS
# Step 8a: Genotype GDS files, genomic relationships, IBD and PCA import
#
#   1) Convert the cohort VCFs of Step 6a to SNPRelate GDS: called genotypes
#      (before imputation) and BEAGLE-imputed genotypes (best-guess calls 0/1/2)
#   2) Genomic relationship matrix (GRM) per cohort with the GCTA estimator
#      from the called genotypes of the within-cohort polymorphic SNPs called
#      in >= 99% of the trees, each pair of trees over the SNPs genotyped in
#      both (tgc.ecs.functions.R); kinship K = GRM / 2
#   3) Import PLINK IBD (Step 7a) and PLINK PCA on LD-pruned SNPs (Step 5a)
#   4) Save R objects and tables used by Steps 9a, 10a, 11ab, 12ab0 and 16ab
#
# The GRM is computed from called rather than imputed genotypes, and from SNPs
# with few missing calls: with a mean ECS depth of about 7 reads per tree, the
# share of missing calls varies among trees (0.8-36%), and imputed genotypes or
# SNPs with many missing calls make poorly covered trees appear related to all
# others.
#
# INPUTS
#   RESULTS/ECS/VCF_SPLIT/tgc.ecs.<cohort>.call.filt.maf05.snvs.poly.vcf.gz
#   RESULTS/ECS/VCF_SPLIT/tgc.ecs.<cohort>.call.filt.maf05.snvs.poly.imputed.vcf.gz
#   RESULTS/ECS/POPGEN/RELATEDNESS/IBD/tgc.ecs.<cohort>.pruned.ibd.genome
#   RESULTS/ECS/POPGEN/STRUCTURE/<COHORT>/PCA/tgc.ecs.<cohort>.pca.pruned.{eigenvec,eigenval}
#   DATA/METADATA/breeding_sample2family.txt, natural_sample2pop.txt
#
# OUTPUTS (RESULTS/ECS/RANALYSIS/)
#   RDATA/<cohort>.called.snp.gds, <cohort>.imputed.snp.gds
#   RDATA/<cohort>_grm_gcta.rds, <cohort>_kinship_gcta_half.rds
#   RDATA/<cohort>_sample_annotation.rds
#   RDATA/<cohort>_ibd_plink_long.rds, <cohort>_ibd_plink_mat.rds
#   RDATA/<cohort>_pca_pruned.rds
#   TABLES/{grm_kinship,ibd,pca}/*.tsv (grm_kinship/grm_snp_counts.tsv: SNPs used for the GRM)
#
# R packages: SNPRelate, gdsfmt, data.table
############################################################

suppressPackageStartupMessages({
  library(SNPRelate)
  library(gdsfmt)
  library(data.table)
})
options(stringsAsFactors = FALSE)

# === USER CONFIGURATION ===
PROJECT_ROOT <- Sys.getenv("TGC_PROJECT_ROOT", unset = "/path/to/your/project")
# ===========================
source(file.path(PROJECT_ROOT, "SCRIPTS/ECS/tgc.ecs.functions.R"))

VCF_SPLIT_DIR <- file.path(PROJECT_ROOT, "RESULTS/ECS/VCF_SPLIT")
STRUCT_DIR    <- file.path(PROJECT_ROOT, "RESULTS/ECS/POPGEN/STRUCTURE")
IBD_DIR       <- file.path(PROJECT_ROOT, "RESULTS/ECS/POPGEN/RELATEDNESS/IBD")
META_DIR      <- file.path(PROJECT_ROOT, "DATA/METADATA")
RANA_DIR      <- file.path(PROJECT_ROOT, "RESULTS/ECS/RANALYSIS")
RDATA_DIR     <- file.path(RANA_DIR, "RDATA")
TAB_DIR       <- file.path(RANA_DIR, "TABLES")
for (d in c(RDATA_DIR, file.path(TAB_DIR, c("grm_kinship", "ibd", "pca"))))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

COHORTS <- list(
  breeding = list(map = "breeding_sample2family.txt", group = "Family"),
  natural  = list(map = "natural_sample2pop.txt",     group = "Population")
)

msg <- function(...) cat("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n", sep = "")
ensure_file <- function(x) if (!file.exists(x)) stop("Missing file: ", x)

vcf2gds <- function(vcf, gds) {
  ensure_file(vcf)
  if (file.exists(gds)) file.remove(gds)
  snpgdsVCF2GDS(vcf.fn = vcf, out.fn = gds, method = "biallelic.only", snpfirstdim = TRUE, verbose = FALSE)
}

read_plink_genome <- function(file, ids) {
  ensure_file(file)
  g <- fread(file, select = c("IID1", "IID2", "PI_HAT"))
  mat <- matrix(NA_real_, length(ids), length(ids), dimnames = list(ids, ids))
  diag(mat) <- 1
  ok <- g$IID1 %in% ids & g$IID2 %in% ids
  mat[cbind(g$IID1[ok], g$IID2[ok])] <- g$PI_HAT[ok]
  mat[cbind(g$IID2[ok], g$IID1[ok])] <- g$PI_HAT[ok]
  list(long = g, mat = mat)
}

grm_counts <- list()
for (cohort in names(COHORTS)) {
  cfg <- COHORTS[[cohort]]
  COH <- toupper(cohort)
  msg("=== ", COH, " ===")

  # 1) VCF -> GDS (called and imputed genotypes)
  stem <- file.path(VCF_SPLIT_DIR, sprintf("tgc.ecs.%s.call.filt.maf05.snvs.poly", cohort))
  gds_called <- file.path(RDATA_DIR, sprintf("%s.called.snp.gds", cohort))
  vcf2gds(paste0(stem, ".vcf.gz"), gds_called)
  vcf2gds(paste0(stem, ".imputed.vcf.gz"), file.path(RDATA_DIR, sprintf("%s.imputed.snp.gds", cohort)))

  # 2) GRM (GCTA estimator, called genotypes)
  map <- fread(file.path(META_DIR, cfg$map), header = FALSE, col.names = c("IID", "group"))
  gf  <- snpgdsOpen(gds_called)
  ids <- read.gdsn(index.gdsn(gf, "sample.id"))
  geno <- snpgdsGetGeno(gf, snpfirstdim = TRUE, with.id = FALSE, verbose = FALSE)
  snpgdsClose(gf)
  mode(geno) <- "numeric"
  geno[geno < 0 | geno > 2] <- NA
  grp <- map$group[match(ids, map$IID)]
  grp[is.na(grp)] <- "Unknown"
  msg(nrow(geno), " SNPs x ", ncol(geno), " samples; missing calls: ",
      round(100 * mean(is.na(geno)), 2), "%")
  keep <- grm_snps(geno)
  comp <- gcta_components(geno[keep, , drop = FALSE])
  G <- gcta_grm(comp, ids)
  msg("GRM from ", comp$n_snps, " SNPs; mean diagonal ", round(mean(diag(G)), 3))
  grm_counts[[cohort]] <- data.table(cohort = cohort, n_snps_polymorphic = nrow(geno), n_snps = comp$n_snps)
  rm(geno, comp); gc()

  annot <- data.table(IID = ids, group = grp)
  setnames(annot, "group", cfg$group)
  saveRDS(G,     file.path(RDATA_DIR, sprintf("%s_grm_gcta.rds", cohort)))
  saveRDS(G / 2, file.path(RDATA_DIR, sprintf("%s_kinship_gcta_half.rds", cohort)))
  saveRDS(annot, file.path(RDATA_DIR, sprintf("%s_sample_annotation.rds", cohort)))
  fwrite(data.table(IID = ids, G), file.path(TAB_DIR, "grm_kinship", sprintf("%s_grm_gcta.tsv", cohort)), sep = "\t")
  fwrite(annot, file.path(TAB_DIR, "grm_kinship", sprintf("%s_sample_annotation.tsv", cohort)), sep = "\t")

  # 3a) IBD
  ibd <- read_plink_genome(file.path(IBD_DIR, sprintf("tgc.ecs.%s.pruned.ibd.genome", cohort)), ids)
  saveRDS(ibd$long, file.path(RDATA_DIR, sprintf("%s_ibd_plink_long.rds", cohort)))
  saveRDS(ibd$mat,  file.path(RDATA_DIR, sprintf("%s_ibd_plink_mat.rds", cohort)))
  fwrite(ibd$long, file.path(TAB_DIR, "ibd", sprintf("%s_ibd_plink_long.tsv", cohort)), sep = "\t")
  msg("IBD pairs: ", nrow(ibd$long))

  # 3b) PCA (PLINK, LD-pruned SNPs)
  pca_prefix <- file.path(STRUCT_DIR, COH, "PCA", sprintf("tgc.ecs.%s.pca.pruned", cohort))
  ensure_file(paste0(pca_prefix, ".eigenvec"))
  ev <- fread(paste0(pca_prefix, ".eigenvec"), header = FALSE)
  setnames(ev, c("FID", "IID", paste0("PC", seq_len(ncol(ev) - 2))))
  eval <- scan(paste0(pca_prefix, ".eigenval"), quiet = TRUE)
  saveRDS(list(scores = ev, eigenval = eval), file.path(RDATA_DIR, sprintf("%s_pca_pruned.rds", cohort)))
  fwrite(ev, file.path(TAB_DIR, "pca", sprintf("%s_pca_pruned_scores.tsv", cohort)), sep = "\t")
  fwrite(data.table(PC = seq_along(eval), eigenval = eval),
         file.path(TAB_DIR, "pca", sprintf("%s_pca_pruned_eigenval.tsv", cohort)), sep = "\t")
}

fwrite(rbindlist(grm_counts), file.path(TAB_DIR, "grm_kinship", "grm_snp_counts.tsv"), sep = "\t")
msg("Step 8a done")
sessionInfo()
```

---

### `9a.tgc.ecs.pca.R`

```r
#!/usr/bin/env Rscript
############################################################
# TreeGeneClimate (TGC) — ECS
# Step 9a: Genetic PCA panel (breeding + natural cohorts)
#
# OUTPUT:
#   RESULTS/ECS/RANALYSIS/PCA/pca_panel.tiff (+ pdf, eps, png)
#   RESULTS/ECS/RANALYSIS/PCA/pca_panel_values.tsv  (values shown in the panel)
############################################################

suppressPackageStartupMessages({
  library(dplyr)      # data wrangling
  library(tibble)     # tidy data frames
  library(ggplot2)    # plotting
  library(patchwork)  # multi-panel layout
  library(grid)       # low-level graphics (used implicitly by patchwork)
  library(scales)     # axis formatting (percent_format, pretty_breaks)
})

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------
# === USER CONFIGURATION ===
PROJECT_ROOT <- Sys.getenv("TGC_PROJECT_ROOT", unset = "/path/to/your/project")
# ===========================
rdata_dir <- file.path(PROJECT_ROOT, "RESULTS/ECS/RANALYSIS/RDATA")
fig_dir   <- file.path(PROJECT_ROOT, "RESULTS/ECS/RANALYSIS/PCA")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# PCA outputs from Step 5a (PLINK PCA on pruned SNP sets)
pca_dir <- file.path(PROJECT_ROOT, "RESULTS/ECS/POPGEN/STRUCTURE")

# PLINK writes one eigenvec and one eigenval file per run; construct both paths
# from a shared prefix so the naming convention is enforced consistently.
breed_prefix <- file.path(pca_dir, "BREEDING", "PCA", "tgc.ecs.breeding.pca.pruned")
natur_prefix <- file.path(pca_dir, "NATURAL",  "PCA", "tgc.ecs.natural.pca.pruned")

eigvec_b_file <- paste0(breed_prefix, ".eigenvec")
eigval_b_file <- paste0(breed_prefix, ".eigenval")
eigvec_n_file <- paste0(natur_prefix, ".eigenvec")
eigval_n_file <- paste0(natur_prefix, ".eigenval")

# Sample annotations (from orchestrator; must exist in RDATA/)
ann_breed <- readRDS(file.path(rdata_dir, "breeding_sample_annotation.rds"))  # IID, Family
ann_nat   <- readRDS(file.path(rdata_dir, "natural_sample_annotation.rds"))   # IID, Population

# ------------------------------------------------------------------------------
# Palettes (centralized colors)
# Fixed, named colour vectors ensure that each family/population always maps
# to the same colour across all figure panels, regardless of plotting order.
# ------------------------------------------------------------------------------
# 17 full-sib families (BREEDING cohort)
colors.17 <- c(
  "Family_16"="dodgerblue2","Family_27"="#E31A1C","Family_32"="green4",
  "Family_33"="#6A3D9A","Family_38"="#FF7F00","Family_39"="black",
  "Family_40"="gold1","Family_41"="skyblue2","Family_42"="#FB9A99",
  "Family_43"="palegreen2","Family_44"="gray70","Family_47"="khaki2",
  "Family_48"="orchid1","Family_50"="deeppink1","Family_51"="blue1",
  "Family_52"="steelblue4","Family_53"="darkturquoise"
)

# 25 Finnish stands (NATURAL cohort)
colors.25 <- c(
  "Asikkala"="dodgerblue2","Jämsä"="#E31A1C","Kauhajoki"="green4",
  "Koski"="#6A3D9A","Kuopio"="#FF7F00","Laihia"="black",
  "Lammi"="gold1","Leppävirta"="skyblue2","Loppi"="#FB9A99",
  "Luopioinen"="palegreen2","Mäntyharju"="#CAB2D6",
  "Marttila"="#FDBF6F","Miehikkälä"="gray80","Mikkeli"="khaki2",
  "Multia"="maroon","Muurame"="orchid1","Orivesi"="deeppink1",
  "Pälkäne"="blue1","Petäjävesi"="steelblue4","Punkaharju"="green1",
  "Punkalaidun"="yellow4","Puumala"="yellow3",
  "Rautalampi"="darkorange4","Savonlinna"="brown","Somero"="grey40"
)

# Restrict a named palette to the groups actually present in the data and
# assign grey70 to any group not covered by the predefined palette.
subset_palette <- function(pal_named, groups) {
  present <- unique(as.character(groups))
  pal <- pal_named[names(pal_named) %in% present]
  missing <- setdiff(present, names(pal))
  if (length(missing)) pal <- c(pal, setNames(rep("grey70", length(missing)), missing))
  pal
}

ensure_file <- function(p) if (!file.exists(p)) stop("Missing file: ", p, call. = FALSE)

# ------------------------------------------------------------------------------
# Readers (compatible across readr/base R versions)
# ------------------------------------------------------------------------------
read_plink_eigen <- function(eigvec_file, eigval_file) {
  ensure_file(eigval_file); ensure_file(eigvec_file)

  # Eigenvalues: robust single-column reader
  eval <- tryCatch(scan(eigval_file, what = numeric(), quiet = TRUE),
                   error = function(e) stop("Failed to read eigenvalues: ", eigval_file))

  # PLINK 1.9 writes eigenvec without a header; PLINK 2 writes one.
  # Try header=TRUE first; if it lacks FID/IID columns, fall back to header=FALSE.
  # Eigenvectors: try header=TRUE first, else header=FALSE
  ev_try <- try(read.table(eigvec_file, header = TRUE, stringsAsFactors = FALSE), silent = TRUE)
  if (inherits(ev_try, "try-error") || !all(c("FID","IID") %in% names(ev_try))) {
    ev <- read.table(eigvec_file, header = FALSE, stringsAsFactors = FALSE)
    nPC <- ncol(ev) - 2
    if (nPC < 1) stop("No PC columns detected in eigenvec: ", eigvec_file)
    colnames(ev) <- c("FID","IID", paste0("PC", seq_len(nPC)))
  } else {
    ev <- ev_try
    pc_cols <- setdiff(names(ev), c("FID","IID"))
    # Normalise column names to PC1, PC2, ... regardless of PLINK's labelling
    if (!all(grepl("^PC\\d+$", pc_cols))) {
      colnames(ev) <- c("FID","IID", paste0("PC", seq_len(length(pc_cols))))
    }
  }

  ev$IID <- as.character(ev$IID)
  list(eigvec = ev, eigval = eval)
}

# ------------------------------------------------------------------------------
# Plot builder
# ------------------------------------------------------------------------------
# PCA scatter for a user-specified pair of PCs (x, y).
# Axis labels include the percentage of variance explained by each PC,
# computed from the full eigenvalue vector so the proportion is global.
make_pca_scatter <- function(df_scores, eigval, x = 1, y = 2, group_col = "Group",
                             palette_named, title_label = "a") {
  stopifnot(paste0("PC", x) %in% names(df_scores), paste0("PC", y) %in% names(df_scores))
  pve <- eigval / sum(eigval)

  lx <- paste0("PC", x, " (", percent(pve[x], accuracy = 0.1), ")")
  ly <- paste0("PC", y, " (", percent(pve[y], accuracy = 0.1), ")")

  # Subset the palette to groups present in this dataset (avoids legend clutter)
  pal <- subset_palette(palette_named, df_scores[[group_col]])

  # .data[[]] tidy-evaluation allows passing column names as strings
  ggplot(df_scores, aes(x = .data[[paste0("PC", x)]],
                        y = .data[[paste0("PC", y)]],
                        color = .data[[group_col]])) +
    geom_point(size = 1.8, alpha = 0.9) +
    scale_color_manual(values = pal, name = group_col) +
    labs(title = title_label, x = lx, y = ly) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0, face = "bold"),
      panel.grid.minor = element_blank(),
      legend.title = element_blank()
    )
}

# ------------------------------------------------------------------------------
# Load PCA results and annotations
# ------------------------------------------------------------------------------
# Breeding
b <- read_plink_eigen(eigvec_b_file, eigval_b_file)
eigvec_b <- b$eigvec
eigval_b <- b$eigval

# Join eigenvectors with family labels; left_join preserves all samples and
# fills unmatched IIDs with "Unknown" rather than silently dropping them.
scores_b <- eigvec_b %>%
  left_join(ann_breed %>% dplyr::select(IID, Group = Family), by = "IID") %>%
  mutate(Group = ifelse(is.na(Group), "Unknown", as.character(Group)))

# Natural
n <- read_plink_eigen(eigvec_n_file, eigval_n_file)
eigvec_n <- n$eigvec
eigval_n <- n$eigval

scores_n <- eigvec_n %>%
  left_join(ann_nat %>% dplyr::select(IID, Group = Population), by = "IID") %>%
  mutate(Group = ifelse(is.na(Group), "Unknown", as.character(Group)))

# ------------------------------------------------------------------------------
# Single 6-plot panel (a-f), TIFF 34 cm x 26 cm @ 600 dpi
#   Row 1: breeding  a) PC1v2  b) PC1v3  c) PC2v3
#   Row 2: natural   d) PC1v2  e) PC1v3  f) PC2v3
#
# PC1 vs PC2, PC1 vs PC3 and PC2 vs PC3 show the three leading axes.
# ------------------------------------------------------------------------------
p_b_12 <- make_pca_scatter(scores_b, eigval_b, x = 1, y = 2, group_col = "Group",
                           palette_named = colors.17, title_label = "a")
p_b_13 <- make_pca_scatter(scores_b, eigval_b, x = 1, y = 3, group_col = "Group",
                           palette_named = colors.17, title_label = "b")
p_b_23 <- make_pca_scatter(scores_b, eigval_b, x = 2, y = 3, group_col = "Group",
                           palette_named = colors.17, title_label = "c")

p_n_12 <- make_pca_scatter(scores_n, eigval_n, x = 1, y = 2, group_col = "Group",
                           palette_named = colors.25, title_label = "d")
p_n_13 <- make_pca_scatter(scores_n, eigval_n, x = 1, y = 3, group_col = "Group",
                           palette_named = colors.25, title_label = "e")
p_n_23 <- make_pca_scatter(scores_n, eigval_n, x = 2, y = 3, group_col = "Group",
                           palette_named = colors.25, title_label = "f")

# patchwork: '|' composes panels horizontally, '/' stacks rows, '&' applies
# theme modifications to all panels in the assembled layout simultaneously.
# Guides are collected per row so breeding and natural legends remain separate
# and align to the top of their respective row's legend area.
pca_row_b <- (p_b_12 | p_b_13 | p_b_23) + plot_layout(guides = "collect") &
  theme(legend.position = "right", legend.justification = "top")
pca_row_n <- (p_n_12 | p_n_13 | p_n_23) + plot_layout(guides = "collect") &
  theme(legend.position = "right", legend.justification = "top")
pca_panel <- pca_row_b / pca_row_n

pca_tiff <- file.path(fig_dir, "pca_panel.tiff")
tiff(pca_tiff, width = 34, height = 26, units = "cm", res = 600, compression = "lzw")
print(pca_panel)
dev.off()
ggsave(sub("\\.tiff$", ".pdf", pca_tiff), plot = pca_panel, width = 34, height = 26, units = "cm")
ggsave(sub("\\.tiff$", ".eps", pca_tiff), plot = pca_panel, device = cairo_ps, width = 34, height = 26, units = "cm")
ggsave(sub("\\.tiff$", ".png", pca_tiff), plot = pca_panel, device = "png",      width = 34, height = 26, units = "cm", dpi = 150)

# Values shown in the panel: PC1-PC3 scores of every tree and the variance explained by each PC
sd_scores <- function(scores, eigval, cohort) {
  pve <- 100 * eigval / sum(eigval)
  data.frame(cohort = cohort, tree = scores$IID, group = scores$Group,
             PC1 = scores$PC1, PC2 = scores$PC2, PC3 = scores$PC3,
             PC1_variance_pct = pve[1], PC2_variance_pct = pve[2], PC3_variance_pct = pve[3])
}
write.table(rbind(sd_scores(scores_b, eigval_b, "BREEDING"), sd_scores(scores_n, eigval_n, "NATURAL")),
            file.path(fig_dir, "pca_panel_values.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat("Saved: ", pca_tiff, "  [PCA panel]\n", sep = "")

sessionInfo()
```

---

### `10a.tgc.ecs.grm.ibd.R`

```r
#!/usr/bin/env Rscript
############################################################
# TreeGeneClimate (TGC) — ECS
# Step 10a: GRM heatmaps + IBD violins
#           (breeding and natural cohorts)
#
# OUTPUT (FIXED; DO NOT CHANGE):
#   /path/to/your/project/RESULTS/ECS/RANALYSIS/FIGURES/GRM_IBD
#
# FIGURE FILES (saved individually):
#   grm_breeding.tiff   (GRM heatmap, breeding)
#   ibd_breeding.tiff   (IBD violin, PI_HAT, breeding)
#   grm_natural.tiff    (GRM heatmap, natural)
#   ibd_natural.tiff    (IBD violin, PI_HAT, natural)
#   grm_ibd_panel.tiff  (combined 2x2 panel, a-d)
#
# VALUES SHOWN IN THE PANEL (same folder):
#   grm_ibd_panel_grm_values.tsv    GRM value of every pair of trees (panels a, c)
#   grm_ibd_panel_pihat_values.tsv  PI_HAT of every pair of trees (panels b, d)
#
# NOTES
# - GRM legend: keep real data range, but print "nice" key labels (few decimals).
# - IBD violins: use PLINK PI_HAT (from 8a outputs), show thresholds:
#     Full-sib ≈ 0.50, Half-sib ≈ 0.25
# - π-hat symbol is used in y-axis label via expression(hat(pi)).
############################################################

suppressPackageStartupMessages({
  library(dplyr)        # data manipulation (joins, mutate, filter)
  library(tibble)       # modern data frames
  library(ggplot2)      # base plotting
  library(scales)       # axis formatting (comma(), pretty_breaks())
  library(grid)         # low-level grid graphics for panel labels and viewports

  library(ComplexHeatmap) # hierarchically clustered heatmaps with annotation tracks
  library(circlize)       # colorRamp2() for continuous colour scales
  library(viridisLite)    # perceptually uniform colour palettes (magma used for GRM)
})

options(stringsAsFactors = FALSE)
set.seed(1)

# ==============================================================================
# 1) FIXED OUTPUT LOCATIONS (DO NOT CHANGE)
# ==============================================================================
# === USER CONFIGURATION ===
# Set TGC_PROJECT_ROOT as an environment variable, or edit the fallback path below
PROJECT_ROOT <- Sys.getenv("TGC_PROJECT_ROOT", unset = "/path/to/your/project")
# ===========================
RANA_DIR     <- file.path(PROJECT_ROOT, "RESULTS/ECS/RANALYSIS")
RDATA_DIR    <- file.path(RANA_DIR, "RDATA")   # RDS objects written by step 8a

FIG_DIR      <- file.path(RANA_DIR, "FIGURES", "GRM_IBD")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# Helper: stop with an informative message if a required file is absent
ensure_file <- function(p) if (!file.exists(p)) stop("Missing file: ", p, call. = FALSE)

# ==============================================================================
# 2) INPUTS (from 8a)
# ==============================================================================
# GCTA GRM matrices from called genotypes (Step 8a; sample × sample, symmetric)
grm_b_file   <- file.path(RDATA_DIR, "breeding_grm_gcta.rds")
grm_n_file   <- file.path(RDATA_DIR, "natural_grm_gcta.rds")

# Sample annotation tables: IID + Family (breeding) or Population (natural)
ann_b_file   <- file.path(RDATA_DIR, "breeding_sample_annotation.rds")
ann_n_file   <- file.path(RDATA_DIR, "natural_sample_annotation.rds")

# PLINK pairwise IBD estimates (long format: IID1, IID2, PI_HAT)
ibd_b_long   <- file.path(RDATA_DIR, "breeding_ibd_plink_long.rds")
ibd_n_long   <- file.path(RDATA_DIR, "natural_ibd_plink_long.rds")

# Verify all required files exist before running any analysis
ensure_file(grm_b_file); ensure_file(grm_n_file)
ensure_file(ann_b_file); ensure_file(ann_n_file)
ensure_file(ibd_b_long); ensure_file(ibd_n_long)

# ==============================================================================
# 3) PALETTES
# ==============================================================================
# Fixed colour assignment for each of the 17 breeding families — named vector
# so colours remain consistent across all figures regardless of plotting order
colors.17 <- c(
  "Family_16"="dodgerblue2","Family_27"="#E31A1C","Family_32"="green4",
  "Family_33"="#6A3D9A","Family_38"="#FF7F00","Family_39"="black",
  "Family_40"="gold1","Family_41"="skyblue2","Family_42"="#FB9A99",
  "Family_43"="palegreen2","Family_44"="gray70","Family_47"="khaki2",
  "Family_48"="orchid1","Family_50"="deeppink1","Family_51"="blue1",
  "Family_52"="steelblue4","Family_53"="darkturquoise"
)

# Fixed colour assignment for each of the 25 Finnish natural stands
colors.25 <- c(
  "Asikkala"="dodgerblue2","Jämsä"="#E31A1C","Kauhajoki"="green4",
  "Koski"="#6A3D9A","Kuopio"="#FF7F00","Laihia"="black",
  "Lammi"="gold1","Leppävirta"="skyblue2","Loppi"="#FB9A99",
  "Luopioinen"="palegreen2","Mäntyharju"="#CAB2D6",
  "Marttila"="#FDBF6F","Miehikkälä"="gray80","Mikkeli"="khaki2",
  "Multia"="maroon","Muurame"="orchid1","Orivesi"="deeppink1",
  "Pälkäne"="blue1","Petäjävesi"="steelblue4","Punkaharju"="green1",
  "Punkalaidun"="yellow4","Puumala"="yellow3",
  "Rautalampi"="darkorange4","Savonlinna"="brown","Somero"="grey40"
)

# Return only the palette entries that actually appear in the data;
# assign grey70 to any group not present in the named palette
subset_palette <- function(pal_named, groups) {
  present <- unique(as.character(groups))
  pal <- pal_named[names(pal_named) %in% present]
  missing <- setdiff(present, names(pal))
  if (length(missing)) pal <- c(pal, setNames(rep("grey70", length(missing)), missing))
  pal
}

# ==============================================================================
# 4) FIGURE STYLE
# ==============================================================================
FS_BASE  <- 12  # base font size for axis text and labels
FS_LABEL <- 16  # panel label font size (A), B), C) …)

# ==============================================================================
# 5) HELPERS
# ==============================================================================

# ---- GRM heatmap with "nice" key labels (few decimals) ----
# Produces a ComplexHeatmap with hierarchical clustering and cohort-colour
# annotation bars on both axes. The legend shows evenly spaced "pretty" breaks
# rather than exact quantiles so tick labels are readable.
plot_grm_heatmap <- function(G, groups, palette_named, out_file, title_label = "a",
                             dpi = 600, w_cm = 20, h_cm = 16,
                             n_breaks = 5, key_digits = 2) {
  G <- as.matrix(G)
  ids <- rownames(G)

  # Match each sample to its group label and build the colour palette
  grp <- as.character(groups[match(ids, names(groups))])
  grp[is.na(grp)] <- "Unknown"
  pal <- subset_palette(palette_named, grp)

  # Annotation tracks shown along the top (columns) and left (rows) of the heatmap.
  # The track name (= legend title) depends on cohort size: 17 → breeding families,
  # 25 → natural stands.
  annot_title <- if (length(pal) == 17L) "Family" else "Natural stand"
  ha_top <- do.call(HeatmapAnnotation, setNames(
    list(grp, setNames(list(pal), annot_title), FALSE),
    c(annot_title, "col", "show_annotation_name")
  ))
  ra_left <- do.call(rowAnnotation, setNames(
    list(grp, setNames(list(pal), annot_title), FALSE),
    c(annot_title, "col", "show_annotation_name")
  ))

  # Colour scale: use 1st–99th percentile of off-diagonal values to avoid
  # distortion from the inflated diagonal (self-relatedness ≈ 1)
  G_offdiag <- G[row(G) != col(G)]
  rng <- quantile(G_offdiag, c(0.01, 0.99), na.rm = TRUE)
  col_fun <- circlize::colorRamp2(seq(rng[1], rng[2], length.out = 5),
                                  viridisLite::viridis(5, option = "magma"))

  # Generate evenly spaced legend breaks within the data range and round to
  # key_digits decimals so labels are compact on the printed figure
  breaks <- pretty(rng, n = n_breaks)
  breaks <- breaks[breaks >= rng[1] & breaks <= rng[2]]
  breaks <- round(breaks, key_digits)

  ht <- Heatmap(
    G,
    name = "GRM",
    col = col_fun,
    cluster_rows = TRUE,       # hierarchical clustering of trees
    cluster_columns = TRUE,
    show_row_names = FALSE,    # sample names omitted (too many to display)
    show_column_names = FALSE,
    top_annotation = ha_top,
    left_annotation = ra_left,
    show_heatmap_legend = TRUE,
    heatmap_legend_param = list(
      title = NULL,
      at = breaks,
      labels = format(breaks, nsmall = key_digits, trim = TRUE),
      labels_gp = gpar(fontsize = 7),
      grid_height = unit(3.0, "cm"),
      grid_width  = unit(0.35, "cm")
    )
  )

  # Nested draw function so the same panel label can be added consistently
  # across TIFF, PDF, EPS and PNG outputs
  draw_heatmap <- function() {
    draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right",
         merge_legend = TRUE, show_annotation_legend = TRUE)
    # Panel label (e.g. "A") placed in the top-left corner using grid coordinates
    grid.text(title_label, x = unit(0.35, "cm"),
              y = unit(1, "npc") - unit(0.35, "cm"),
              just = c("left","top"), gp = gpar(fontsize = FS_LABEL, fontface = "bold"))
  }
  # Write TIFF at 600 dpi with lossless LZW compression
  tiff(out_file, width = w_cm, height = h_cm, units = "cm", res = dpi, compression = "lzw")
  draw_heatmap()
  dev.off()
  cairo_pdf(sub("\\.tiff$", ".pdf", out_file), width = w_cm / 2.54, height = h_cm / 2.54)
  draw_heatmap()
  dev.off()
  cairo_ps(sub("\\.tiff$", ".eps", out_file), width = w_cm / 2.54, height = h_cm / 2.54)
  draw_heatmap()
  dev.off()
  png(sub("\\.tiff$", ".png", out_file), width = w_cm, height = h_cm, units = "cm", res = dpi)
  draw_heatmap()
  dev.off()
  invisible(list(file = out_file, ht = ht))  # return heatmap object for panel assembly
}

# ---- Prepare Within/Between using annotation (for PI_HAT pairs) ----
# Joins the pairwise IBD table with sample annotations to classify each pair
# as Within-group (same family or population) or Between-group
prep_pairs_pihat <- function(df_long, annot_df, group_col) {
  stopifnot(all(c("IID1","IID2","PI_HAT") %in% names(df_long)))
  stopifnot(all(c("IID", group_col) %in% names(annot_df)))

  # Extract group labels for each member of the pair separately, then rejoin
  a1 <- annot_df %>% dplyr::select(IID, g1 = !!rlang::sym(group_col))
  a2 <- annot_df %>% dplyr::select(IID, g2 = !!rlang::sym(group_col))

  df_long %>%
    transmute(s1 = as.character(IID1), s2 = as.character(IID2), PI_HAT = as.numeric(PI_HAT)) %>%
    left_join(a1, by = c("s1" = "IID")) %>%
    left_join(a2, by = c("s2" = "IID")) %>%
    mutate(
      g1 = ifelse(is.na(g1), "Unknown", as.character(g1)),
      g2 = ifelse(is.na(g2), "Unknown", as.character(g2)),
      # Within = both individuals belong to the same known group
      PairType = ifelse(g1 == g2 & g1 != "Unknown", "Within", "Between")
    )
}

# ---- PI_HAT violin (16x16), π-hat label, FS/Half-sib refs ----
# Violin + embedded boxplot comparing IBD (PI_HAT) distributions between
# within-group and between-group pairs. Reference lines at 0.25 (half-sib)
# and 0.50 (full-sib) mark the expected values.
plot_pihat_violin <- function(df, out_file, title_label = "b",
                              dpi = 600, w_cm = 16, h_cm = 16,
                              y_limits = c(-0.05, 0.70)) {

  df <- df %>% filter(is.finite(PI_HAT))  # remove rare NA/Inf from PLINK output
  df$PairType <- factor(df$PairType, levels = c("Between", "Within"))

  # Build sample-size labels (e.g. "Between (n=1,234)") for the x-axis
  counts <- df %>% count(PairType) %>%
    mutate(lab = paste0(PairType, " (n=", scales::comma(n), ")"))
  lbl <- setNames(counts$lab, counts$PairType)

  p <- ggplot(df, aes(x = PairType, y = PI_HAT, fill = PairType)) +
    geom_violin(trim = FALSE, alpha = 0.7, color = "grey25", width = 0.9) +
    # Narrow boxplot embedded inside each violin; outliers hidden (shown by violin)
    geom_boxplot(width = 0.18, outlier.alpha = 0, fill = "white", color = "grey25") +
    # Horizontal reference lines for expected relatedness in full-sib vs half-sib pairs
    geom_hline(yintercept = c(0.25, 0.50),
               linetype = c("dotted","dashed"),
               color = c("grey60","grey50")) +
    # Annotations use parse = TRUE to render the π-hat (hat(pi)) mathematical symbol
    annotate("text", x = 1.5, y = 0.50,
             label = "paste('≈ Full-sib (', hat(pi), ' = 0.50)')",
             parse = TRUE, vjust = -0.6, size = 3.2) +
    annotate("text", x = 1.5, y = 0.25,
             label = "paste('≈ Half-sib (', hat(pi), ' = 0.25)')",
             parse = TRUE, vjust = -0.6, size = 3.2) +

    scale_fill_manual(values = c(Between = "#9ecae1", Within = "#3182bd")) +
    coord_cartesian(ylim = y_limits, expand = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.03))) +
    scale_x_discrete(labels = function(x) lbl[x]) +  # replace axis labels with n-annotated versions
    theme_minimal(base_size = FS_BASE) +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.title.y = element_text(margin = margin(r = 8)),
      axis.text.x  = element_text(margin = margin(t = 6)),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 8, r = 12, b = 8, l = 12)
    ) +
    # y-axis label uses R expression() to render the π-hat symbol correctly
    labs(y = expression("Proportion of genome shared IBD ("*hat(pi)*")"))

  # Nested draw function adds the panel label after the plot is printed
  draw_labeled <- function() {
    print(p)
    grid.text(title_label,
              x = unit(0.35, "cm"),
              y = unit(1, "npc") - unit(0.35, "cm"),
              just = c("left","top"),
              gp = gpar(fontsize = FS_LABEL, fontface = "bold"))
  }
  tiff(out_file, width = w_cm, height = h_cm, units = "cm", res = dpi, compression = "lzw")
  draw_labeled()
  dev.off()
  cairo_pdf(sub("\\.tiff$", ".pdf", out_file), width = w_cm / 2.54, height = h_cm / 2.54)
  draw_labeled()
  dev.off()
  cairo_ps(sub("\\.tiff$", ".eps", out_file), width = w_cm / 2.54, height = h_cm / 2.54)
  draw_labeled()
  dev.off()
  png(sub("\\.tiff$", ".png", out_file), width = w_cm, height = h_cm, units = "cm", res = dpi)
  draw_labeled()
  dev.off()
  invisible(list(file = out_file, plot = p))  # return plot object for panel assembly
}

# ---- Combined 2×2 panel: A) GRM breeding, B) IBD breeding,
#                          C) GRM natural,  D) IBD natural ----
# Uses grid viewports to tile four sub-figures into a single 32×32 cm canvas.
# 2x2 panel layout.
save_grm_ibd_panel <- function(ht_b, ht_n, p_ibd_b, p_ibd_n, out_file,
                             dpi = 600, w_cm = 38, h_cm = 32) {
  draw_panel <- function() {
    grid.newpage()

    # A) GRM breeding — top-left quadrant (x: 0–0.5, y: 0.5–1.0)
    pushViewport(viewport(x = 0, y = 0.5, width = 0.5, height = 0.5,
                          just = c("left", "bottom")))
    draw(ht_b, newpage = FALSE,
         heatmap_legend_side = "right", annotation_legend_side = "right",
         merge_legend = TRUE, show_annotation_legend = TRUE)
    grid.text("a", x = unit(0.35, "cm"), y = unit(1, "npc") - unit(0.35, "cm"),
              just = c("left", "top"), gp = gpar(fontsize = FS_LABEL, fontface = "bold"))
    popViewport()

    # B) IBD breeding — top-right quadrant
    pushViewport(viewport(x = 0.5, y = 0.5, width = 0.5, height = 0.5,
                          just = c("left", "bottom")))
    grid.draw(ggplotGrob(p_ibd_b))  # convert ggplot to grob for viewport drawing
    grid.text("b", x = unit(0.35, "cm"), y = unit(1, "npc") - unit(0.35, "cm"),
              just = c("left", "top"), gp = gpar(fontsize = FS_LABEL, fontface = "bold"))
    popViewport()

    # C) GRM natural — bottom-left quadrant
    pushViewport(viewport(x = 0, y = 0, width = 0.5, height = 0.5,
                          just = c("left", "bottom")))
    draw(ht_n, newpage = FALSE,
         heatmap_legend_side = "right", annotation_legend_side = "right",
         merge_legend = TRUE, show_annotation_legend = TRUE)
    grid.text("c", x = unit(0.35, "cm"), y = unit(1, "npc") - unit(0.35, "cm"),
              just = c("left", "top"), gp = gpar(fontsize = FS_LABEL, fontface = "bold"))
    popViewport()

    # D) IBD natural — bottom-right quadrant
    pushViewport(viewport(x = 0.5, y = 0, width = 0.5, height = 0.5,
                          just = c("left", "bottom")))
    grid.draw(ggplotGrob(p_ibd_n))
    grid.text("d", x = unit(0.35, "cm"), y = unit(1, "npc") - unit(0.35, "cm"),
              just = c("left", "top"), gp = gpar(fontsize = FS_LABEL, fontface = "bold"))
    popViewport()
  }

  tiff(out_file, width = w_cm, height = h_cm, units = "cm", res = dpi, compression = "lzw")
  draw_panel(); dev.off()
  cairo_pdf(sub("\\.tiff$", ".pdf", out_file), width = w_cm / 2.54, height = h_cm / 2.54)
  draw_panel(); dev.off()
  cairo_ps(sub("\\.tiff$", ".eps", out_file), width = w_cm / 2.54, height = h_cm / 2.54)
  draw_panel(); dev.off()
  png(sub("\\.tiff$", ".png", out_file), width = w_cm, height = h_cm, units = "cm", res = dpi)
  draw_panel(); dev.off()
  invisible(out_file)
}

# ==============================================================================
# 6) LOAD DATA
# ==============================================================================
G_b   <- readRDS(grm_b_file)   # GRM matrix, breeding cohort
G_n   <- readRDS(grm_n_file)   # GRM matrix, natural cohort
ann_b <- readRDS(ann_b_file)   # sample annotation: IID, Family
ann_n <- readRDS(ann_n_file)   # sample annotation: IID, Population

ibd_b <- readRDS(ibd_b_long)   # pairwise IBD, breeding (IID1, IID2, PI_HAT)
ibd_n <- readRDS(ibd_n_long)   # pairwise IBD, natural

# Named vectors IID → group label, used for annotation and palette subsetting
grp_b_vec <- setNames(as.character(ann_b$Family), ann_b$IID)
grp_n_vec <- setNames(as.character(ann_n$Population), ann_n$IID)

# ==============================================================================
# 7) a / c : GRM heatmaps with family/stand annotation legends
# ==============================================================================
out_a <- file.path(FIG_DIR, "grm_breeding.tiff")
out_c <- file.path(FIG_DIR, "grm_natural.tiff")
res_a <- plot_grm_heatmap(G_b, grp_b_vec, colors.17, out_a, title_label = "a", key_digits = 2)
res_c <- plot_grm_heatmap(G_n, grp_n_vec, colors.25, out_c, title_label = "c", key_digits = 2,
                            w_cm = 24)

# ==============================================================================
# 8) b / d : IBD violins using PLINK PI_HAT (16x16 cm)
# ==============================================================================
# Classify each pairwise comparison as Within or Between family/population
pihat_b <- prep_pairs_pihat(ibd_b, ann_b, "Family")
pihat_n <- prep_pairs_pihat(ibd_n, ann_n, "Population")

out_b <- file.path(FIG_DIR, "ibd_breeding.tiff")
out_d <- file.path(FIG_DIR, "ibd_natural.tiff")

# y_limits include slight negative slack to prevent clipping of "Between" violin tails
res_b <- plot_pihat_violin(pihat_b, out_b, title_label = "b", y_limits = c(-0.05, 0.70))
res_d <- plot_pihat_violin(pihat_n, out_d, title_label = "d", y_limits = c(-0.05, 0.70))

# ==============================================================================
# 9) Combined 2×2 panel (GRM + IBD, breeding + natural)
# ==============================================================================
out_panel <- file.path(FIG_DIR, "grm_ibd_panel.tiff")
save_grm_ibd_panel(res_a$ht, res_c$ht, res_b$plot, res_d$plot, out_panel,
                dpi = 600, w_cm = 32, h_cm = 32)

cat(
  "Saved figures in:\n  ", FIG_DIR, "\n\n",
  "  - ", basename(out_a), "\n",
  "  - ", basename(out_b), "\n",
  "  - ", basename(out_c), "\n",
  "  - ", basename(out_d), "\n",
  "  - ", basename(out_panel), "  [combined panel]\n",
  sep = ""
)

# ==============================================================================
# 10) VALUES SHOWN IN THE PANEL
# ==============================================================================
# GRM: upper triangle with the diagonal (the matrices are symmetric)
grm_long <- function(G, groups, cohort, panel) {
  G <- as.matrix(G)
  ut <- which(upper.tri(G, diag = TRUE), arr.ind = TRUE)
  ids <- rownames(G)
  data.frame(panel = panel, cohort = cohort,
             tree_1 = ids[ut[, 1]], group_1 = unname(groups[ids[ut[, 1]]]),
             tree_2 = ids[ut[, 2]], group_2 = unname(groups[ids[ut[, 2]]]), grm = G[ut])
}
pihat_long <- function(d, cohort, panel) {
  d <- d[is.finite(d$PI_HAT), ]
  data.frame(panel = panel, cohort = cohort, tree_1 = d$s1, group_1 = d$g1, tree_2 = d$s2,
             group_2 = d$g2, pair_type = d$PairType, pi_hat = d$PI_HAT)
}
write.table(rbind(grm_long(G_b, grp_b_vec, "BREEDING", "a"), grm_long(G_n, grp_n_vec, "NATURAL", "c")),
            file.path(FIG_DIR, "grm_ibd_panel_grm_values.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(rbind(pihat_long(pihat_b, "BREEDING", "b"), pihat_long(pihat_n, "NATURAL", "d")),
            file.path(FIG_DIR, "grm_ibd_panel_pihat_values.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

sessionInfo()
```

---

### `tgc.ecs.functions.R`

```r
############################################################
# TreeGeneClimate (TGC) — ECS
# Shared functions for Steps 8a and 12ab0 (sourced, not run directly)
#
# Genomic relationship matrix (GRM) with the GCTA estimator from called
# (non-imputed) genotypes. As in GCTA,
# each pair of trees is evaluated over the SNPs genotyped in both, so trees
# with more missing genotypes are not pulled towards the population mean (as
# happens when missing genotypes are imputed or set to their expected value).
#
#   off-diagonal  A_jk = sum_i z_ij z_ik / N_jk,  z_ij = (x_ij - 2p_i) / sqrt(2p_i(1 - p_i))
#   diagonal      A_jj = 1 + sum_i [x_ij^2 - (1 + 2p_i) x_ij + 2p_i^2] / [2p_i(1 - p_i)] / N_jj
#
# Only SNPs called in at least 99% of the cohort's trees are used (GRM_MAX_MISS):
# with about 7 reads per tree, genotype calls of poorly covered trees are biased
# (missing calls, heterozygotes called as homozygotes), and SNPs with many
# missing calls make the GRM reflect sequencing depth rather than relatedness.
#
# The sums are additive over SNPs, so leave-one-chromosome-out (LOCO) GRMs are
# obtained exactly by subtracting the components of one chromosome.
############################################################

GRM_MAX_MISS <- 0.01

# Rows (SNPs) of X kept for the GRM
grm_snps <- function(X, max_miss = GRM_MAX_MISS) rowMeans(is.na(X)) <= max_miss

# Components of the GCTA GRM from a genotype matrix X (SNPs x samples, values
# 0/1/2 or NA). p: allele frequencies (default: from X itself)
gcta_components <- function(X, p = rowMeans(X, na.rm = TRUE) / 2) {
  ok <- is.finite(p) & p > 0 & p < 1
  X <- X[ok, , drop = FALSE]; p <- p[ok]
  het <- 2 * p * (1 - p)
  called <- !is.na(X)
  Z <- (X - 2 * p) / sqrt(het)
  W <- (X^2 - (1 + 2 * p) * X + 2 * p^2) / het
  Z[!called] <- 0; W[!called] <- 0
  storage.mode(called) <- "double"
  list(num = crossprod(Z), cnt = crossprod(called), dnum = colSums(W), n_snps = nrow(X))
}

# Subtract (sign = -1) or add components
gcta_combine <- function(a, b, sign = -1) {
  list(num = a$num + sign * b$num, cnt = a$cnt + sign * b$cnt,
       dnum = a$dnum + sign * b$dnum, n_snps = a$n_snps + sign * b$n_snps)
}

# GRM from components; ids: sample names
gcta_grm <- function(comp, ids) {
  A <- comp$num / comp$cnt
  diag(A) <- 1 + comp$dnum / diag(comp$cnt)
  dimnames(A) <- list(ids, ids)
  A
}
```

---

## TMS — Targeted Enzymatic Methylation Sequencing

---

### `1b.tgc.tms.fastqc.rawdata.sh`

```bash
#!/bin/bash
#-------------------------------------------------------------------------------
# TreeGeneClimate (TGC) — TMS
# Step 1b: FASTQC + MULTIQC on RAW FASTQ (TMS)
#
# Project root:
#   /path/to/your/project
#
# Input (data, frozen):
#   DATA/TMS/RAWDATA.TMS/*.fastq.gz
#
# Output (results):
#   RESULTS/TMS/QC/RAWDATA/FASTQC/
#   RESULTS/TMS/QC/RAWDATA/MULTIQC/
#-------------------------------------------------------------------------------

#SBATCH -n 48
#SBATCH -N 1
#SBATCH --job-name=TMS.FQC1
#SBATCH --output=LOGS/%x_%j.out
#SBATCH --error=LOGS/%x_%j.err
#SBATCH --ntasks-per-socket 24
#SBATCH --time=48:00:00

set -euo pipefail

echo "[$(date)] SLURM job started: ${SLURM_JOB_NAME:-no_slurm}"

module load fastqc/0.11.4
module load anaconda3/2020.11

# === USER CONFIGURATION ===
PROJECT_ROOT="${TGC_PROJECT_ROOT:-/path/to/your/project}"
# ===========================

INPUT="${PROJECT_ROOT}/DATA/TMS/RAWDATA.TMS"

QC_BASE="${PROJECT_ROOT}/RESULTS/TMS/QC/RAWDATA"
QC_FASTQC="${QC_BASE}/FASTQC"
QC_MULTIQC="${QC_BASE}/MULTIQC"

LOGS="${PROJECT_ROOT}/LOGS"
mkdir -p "${QC_FASTQC}" "${QC_MULTIQC}" "${LOGS}"

shopt -s nullglob
raw_fastq=( "${INPUT}"/*.fastq.gz )
if (( ${#raw_fastq[@]} == 0 )); then
  echo "ERROR: no FASTQ found in: ${INPUT}"
  exit 1
fi

fastqc "${raw_fastq[@]}" --outdir "${QC_FASTQC}" --threads 48

source activate multiqc
multiqc "${QC_FASTQC}" -o "${QC_MULTIQC}"
conda deactivate

echo "[$(date)] Done."
exit 0
```

---

### `2b.tgc.tms.trimmomatic.and.fastqc.trimmed.sh`

```bash
#!/bin/bash
#-------------------------------------------------------------------------------
# TreeGeneClimate (TGC) — TMS
# Step 2b: TRIMMOMATIC (paired-end) + FASTQC + MULTIQC on TRIMMED FASTQ (TMS)
#
# Project root:
#   /path/to/your/project
#
# Input (data, frozen):
#   DATA/TMS/RAWDATA.TMS/HEL_101202_*_R1.fastq.gz
#   DATA/TMS/RAWDATA.TMS/HEL_101202_*_R2.fastq.gz
#
# Output (data, frozen):
#   DATA/TMS/TRIMMED.FASTQ.TMS/           (paired reads)
#   DATA/TMS/TRIMMED.FASTQ.TMS/UNPAIRED/  (unpaired reads)
#
# Output (results):
#   RESULTS/TMS/QC/TRIMMED/FASTQC/
#   RESULTS/TMS/QC/TRIMMED/MULTIQC/
#
# Adapters:
#   DATA/METADATA/adapters.fa
#-------------------------------------------------------------------------------

#SBATCH -n 48
#SBATCH -N 1
#SBATCH --job-name=TMS.TRIM
#SBATCH --output=LOGS/%x_%j.out
#SBATCH --error=LOGS/%x_%j.err
#SBATCH --ntasks-per-socket 24
#SBATCH --time=48:00:00

set -euo pipefail

echo "[$(date)] SLURM job started: ${SLURM_JOB_NAME:-no_slurm}"

module purge
module load fastqc/0.11.4
module load anaconda3/2020.11
module load trimmomatic/0.36

# === USER CONFIGURATION ===
PROJECT_ROOT="${TGC_PROJECT_ROOT:-/path/to/your/project}"
# ===========================

INPUT="${PROJECT_ROOT}/DATA/TMS/RAWDATA.TMS"
TRIMMED="${PROJECT_ROOT}/DATA/TMS/TRIMMED.FASTQ.TMS"
UNPAIRED="${TRIMMED}/UNPAIRED"

QC_BASE="${PROJECT_ROOT}/RESULTS/TMS/QC/TRIMMED"
QC_FASTQC="${QC_BASE}/FASTQC"
QC_MULTIQC="${QC_BASE}/MULTIQC"

ADAPTERS="${PROJECT_ROOT}/DATA/METADATA/adapters.fa"
LOGS="${PROJECT_ROOT}/LOGS"

mkdir -p "${TRIMMED}" "${UNPAIRED}" "${QC_FASTQC}" "${QC_MULTIQC}" "${LOGS}"

if [[ ! -s "${ADAPTERS}" ]]; then
  echo "ERROR: adapters file not found or empty: ${ADAPTERS}"
  exit 1
fi

list=$(cat "${PROJECT_ROOT}/DATA/METADATA/samples_tms_sequenced_620.txt")

for sample in ${list}; do
  r1="${INPUT}/HEL_101202_${sample}_R1.fastq.gz"
  r2="${INPUT}/HEL_101202_${sample}_R2.fastq.gz"

  if [[ ! -s "${r1}" || ! -s "${r2}" ]]; then
    echo "WARNING: missing pair for ${sample}; skipping"
    continue
  fi

  java -jar "${TRIMMOMATIC_JAR:-/path/to/trimmomatic-0.36.jar}" PE \
    -threads 48 -phred33 \
    "${r1}" "${r2}" \
    "${TRIMMED}/${sample}_R1_p.fastq.gz" "${UNPAIRED}/${sample}_R1_u.fastq.gz" \
    "${TRIMMED}/${sample}_R2_p.fastq.gz" "${UNPAIRED}/${sample}_R2_u.fastq.gz" \
    ILLUMINACLIP:"${ADAPTERS}":2:30:10 \
    CROP:138 HEADCROP:12 SLIDINGWINDOW:5:20 MINLEN:30
done

shopt -s nullglob
trim_fastq=( "${TRIMMED}"/*_p.fastq.gz )
if (( ${#trim_fastq[@]} == 0 )); then
  echo "ERROR: no trimmed paired FASTQ produced in: ${TRIMMED}"
  exit 1
fi

fastqc "${trim_fastq[@]}" --outdir "${QC_FASTQC}" --threads 48

source activate multiqc
multiqc "${QC_FASTQC}" -o "${QC_MULTIQC}"
conda deactivate

echo "[$(date)] Done."
exit 0
```

---

### `3b.tgc.tms.bismark.trimmed.sh`

```bash
#!/bin/bash
#-------------------------------------------------------------------------------
# TreeGeneClimate (TGC) — TMS
# Step 3b: BISMARK mapping on TRIMMED FASTQ (TMS), plates P001–P008
#
# Project root:
#   /path/to/your/project
#
# Input (data, frozen):
#   DATA/TMS/TRIMMED.FASTQ.TMS/*_R1_p.fastq.gz
#   DATA/TMS/TRIMMED.FASTQ.TMS/*_R2_p.fastq.gz
#
# Reference (frozen):
#   REFERENCE/Picab02_genome.assembly/bismark/
#
# Output (data, frozen):
#   DATA/TMS/MAPPED.FILES.TMS/<PLATE>/
#-------------------------------------------------------------------------------

#SBATCH -n 96
#SBATCH -N 1
#SBATCH --job-name=TMS.BISMARK
#SBATCH --output=LOGS/%x_%j.out
#SBATCH --error=LOGS/%x_%j.err
#SBATCH --exclusive
#SBATCH --time=48:00:00

set -euo pipefail

echo "[$(date)] SLURM job started: ${SLURM_JOB_NAME:-no_slurm}"

module purge
module load miniforge3/24.3.0-0
module load gcc/14.2.0
module load bowtie2/2.5.4
module load samtools/1.21
source activate bismark

# === USER CONFIGURATION ===
PROJECT_ROOT="${TGC_PROJECT_ROOT:-/path/to/your/project}"
# ===========================

INPUT="${PROJECT_ROOT}/DATA/TMS/TRIMMED.FASTQ.TMS"
OUTBASE="${PROJECT_ROOT}/DATA/TMS/MAPPED.FILES.TMS"
REF="${PROJECT_ROOT}/REFERENCE/Picab02_genome.assembly/bismark"

LOGS="${PROJECT_ROOT}/LOGS"
mkdir -p "${OUTBASE}" "${LOGS}"

list=$(cat "${PROJECT_ROOT}/DATA/METADATA/samples_tms_sequenced_620.txt")

for sample in ${list}; do
  plate="${sample%%_*}"
  outdir="${OUTBASE}/${plate}"

  r1="${INPUT}/${sample}_R1_p.fastq.gz"
  r2="${INPUT}/${sample}_R2_p.fastq.gz"

  if [[ ! -s "${r1}" || ! -s "${r2}" ]]; then
    echo "WARNING: missing pair for ${sample}; skipping"
    continue
  fi

  mkdir -p "${outdir}"

  bismark -q -L 15 -N 1 --score_min L,-1,-1 \
    --parallel 4 --un --ambiguous \
    -o "${outdir}" \
    --genome "${REF}" \
    -1 "${r1}" -2 "${r2}"
done

conda deactivate
echo "[$(date)] Done."
exit 0
```

---

### `3b2.tgc.tms.deduplication.sh`

```bash
#!/bin/bash
#SBATCH --job-name=TMS.3b2.DEDUP
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --array=0-601%60
#SBATCH --output=LOGS/3b2/%x_%A_%a.out
#SBATCH --error=LOGS/3b2/%x_%A_%a.err
#-------------------------------------------------------------------------------
# TreeGeneClimate (TGC) — TMS
# Step 3b2: Removal of PCR duplicates from Bismark alignments (one sample per task)
#
#   samtools sort -n  ->  samtools fixmate -m  ->  samtools sort  ->  samtools markdup -r
#   (pairs with identical outer coordinates and orientation are collapsed to one)
#
# Input (Step 3b; Bismark paired-end alignments, any sort order):
#   DATA/TMS/MAPPED.FILES.TMS/<SAMPLE>${IN_SUFFIX}
# Output:
#   DATA/TMS/MAPPED.FILES.TMS/<SAMPLE>_R1_p_bismark_bt2_pe.dedup.bam (+ .csi)
#   RESULTS/TMS/QC/DEDUP/<SAMPLE>.markdup.txt     (samtools markdup statistics)
# Summary (TGC_MODE=summary, run once after the array):
#   RESULTS/TMS/QC/DEDUP/dedup_summary.tsv            duplicate reads per sample
#   RESULTS/TMS/QC/DEDUP/dedup_summary_by_cohort.tsv  median percentage per cohort
#
# Samples: DATA/METADATA/samples_final_602.txt (array task i -> line i+1)
# Create LOGS/3b2 before submitting.
# Software: SAMtools 1.21
#-------------------------------------------------------------------------------
set -euo pipefail

# === USER CONFIGURATION ===
PROJECT_ROOT="${TGC_PROJECT_ROOT:-/path/to/your/project}"
IN_SUFFIX="${TGC_BAM_IN_SUFFIX:-_R1_p_bismark_bt2_pe.sorted.bam}"
THREADS="${SLURM_CPUS_PER_TASK:-4}"
module purge 2>/dev/null || true
module load gcc/14.2.0 samtools/1.21 2>/dev/null || true
# ===========================

BAM_DIR="${PROJECT_ROOT}/DATA/TMS/MAPPED.FILES.TMS"
QC_DIR="${PROJECT_ROOT}/RESULTS/TMS/QC/DEDUP"
mkdir -p "${QC_DIR}"

# Summary mode: percentage of reads removed as duplicates per sample, and the median
# per cohort (breeding: the samples of breeding_sample2family.txt; natural: the others)
if [[ "${TGC_MODE:-}" == "summary" ]]; then
  map="${PROJECT_ROOT}/DATA/METADATA/breeding_sample2family.txt"
  out="${QC_DIR}/dedup_summary.tsv"
  {
    printf "sample\tcohort\treads\tduplicates\tduplicate_pct\n"
    for f in "${QC_DIR}"/*.markdup.txt; do
      awk -v s="$(basename "${f}" .markdup.txt)" -v map="${map}" '
        BEGIN { FS = ": *"; while ((getline l < map) > 0) { split(l, a, "\t"); b[a[1]] = 1 } }
        /^READ:/ { r = $2 }
        /^DUPLICATE TOTAL:/ { d = $2 }
        END { printf "%s\t%s\t%d\t%d\t%.6f\n", s, (s in b) ? "BREEDING" : "NATURAL", r, d, 100 * d / r }' "${f}"
    done
  } > "${out}"
  {
    printf "cohort\tsamples\tmedian_duplicate_pct\n"
    for c in BREEDING NATURAL; do
      awk -F'\t' -v c="${c}" '$2 == c { print $5 }' "${out}" | sort -g |
        awk -v c="${c}" '{ v[NR] = $1 } END { m = (NR % 2) ? v[(NR + 1) / 2] : (v[NR / 2] + v[NR / 2 + 1]) / 2
                                              printf "%s\t%d\t%.6f\n", c, NR, m }'
    done
  } > "${QC_DIR}/dedup_summary_by_cohort.tsv"
  echo "[$(date)] summary: ${out}"
  exit 0
fi
i="${SLURM_ARRAY_TASK_ID:-0}"
sample=$(sed -n "$((i + 1))p" "${PROJECT_ROOT}/DATA/METADATA/samples_final_602.txt")
[[ -n "${sample}" ]] || { echo "Task ${i}: no sample"; exit 0; }
in_bam="${BAM_DIR}/${sample}${IN_SUFFIX}"
out_bam="${BAM_DIR}/${sample}_R1_p_bismark_bt2_pe.dedup.bam"
[[ -s "${in_bam}" ]] || { echo "ERROR: missing ${in_bam}"; exit 1; }
tmp=$(mktemp -d "${TMPDIR:-/tmp}/dedup_${sample}_XXXX")
trap 'rm -rf "${tmp}"' EXIT

samtools sort -n -@ "${THREADS}" -T "${tmp}/n" -o "${tmp}/n.bam" "${in_bam}"
samtools fixmate -m -@ "${THREADS}" "${tmp}/n.bam" "${tmp}/fix.bam"
samtools sort -@ "${THREADS}" -T "${tmp}/c" -o "${tmp}/c.bam" "${tmp}/fix.bam"
samtools markdup -r -s -@ "${THREADS}" "${tmp}/c.bam" "${out_bam}" 2> "${QC_DIR}/${sample}.markdup.txt"
samtools index -c "${out_bam}"
echo "[$(date)] ${sample}: $(grep -E '^READ:|^DUPLICATE TOTAL:' "${QC_DIR}/${sample}.markdup.txt" | tr '\n' ' ')"
```

---

### `4b1.tgc.tms.meth.extractor.array.sh`

```bash
#!/bin/bash
#-------------------------------------------------------------------------------
# TreeGeneClimate (TGC) — TMS
# Step 4b1: Bismark methylation extractor (CpG, CHG, CHH)
# Array-based, 1 sample per task
#
# Input:
#   DATA/TMS/MAPPED.FILES.TMS/*_sorted.bam
#
# Output:
#   RESULTS/TMS/METHYLATION_CALLS/<SAMPLE>/
#-------------------------------------------------------------------------------

#SBATCH -t 00-48:00:00
#SBATCH -N 1
#SBATCH -c 24
#SBATCH --mem=100G
#SBATC -C ssd
#SBATCH --job-name=TMS.MEX
#SBATCH --output=LOGS/%x_%A_%a.out
#SBATCH --error=LOGS/%x_%A_%a.err
#SBATCH --array=20-601%4

set -euo pipefail

module purge
module load gcc/14.2.0
module load bismark/0.25.1
module load bowtie2/2.5.4
module load samtools/1.21

# === USER CONFIGURATION ===
PROJECT_ROOT="${TGC_PROJECT_ROOT:-/path/to/your/project}"
# ===========================

INPUT="${PROJECT_ROOT}/DATA/TMS/MAPPED.FILES.TMS"
OUTPUT_BASE="${PROJECT_ROOT}/RESULTS/TMS/METHYLATION_CALLS"
REF="${PROJECT_ROOT}/REFERENCE/Picab02_genome.assembly/bismark"

mkdir -p "${OUTPUT_BASE}"

# -----------------------------
# Full 602-sample list
# -----------------------------
mapfile -t samples < "${PROJECT_ROOT}/DATA/METADATA/samples_final_602.txt"

SAMPLE="${samples[$SLURM_ARRAY_TASK_ID]}"

if [[ -z "${SAMPLE:-}" ]]; then
  echo "ERROR: No sample found for index ${SLURM_ARRAY_TASK_ID}"
  exit 1
fi

BAM="${INPUT}/${SAMPLE}_R1_p_bismark_bt2_pe.bam"
OUTDIR="${OUTPUT_BASE}/${SAMPLE}"

echo "============================================================"
echo "Sample: ${SAMPLE}"
echo "BAM: ${BAM}"
echo "Node: $(hostname)"
echo "Start: $(date)"
echo "============================================================"

if [[ ! -f "$BAM" ]]; then
  echo "ERROR: BAM not found: $BAM"
  exit 1
fi

mkdir -p "$OUTDIR"

bismark_methylation_extractor \
  --multicore 24 \
  --buffer_size 8G \
  --gzip \
  --bedGraph \
  --no_overlap \
  --CX_context \
  --genome_folder "$REF" \
  -o "$OUTDIR" \
  "$BAM"

echo "Finished: $(date)"
exit 0
```

---

### `4b2.tgc.tms.meth.filtering.array.sh`

```bash
#!/bin/bash
#-------------------------------------------------------------------------------
# TreeGeneClimate (TGC) — TMS
# Step 4b2: Coverage filtering of Bismark coverage files
#   Keep sites with 5 <= coverage <= 50
#   Separate outputs by cohort (BREEDING vs NATURAL)
#
# Input:
#   RESULTS/TMS/METHYLATION_CALLS/<SAMPLE>/*.bismark.cov.gz
#
# Output:
#   RESULTS/TMS/METHYLATION_FILTERED/<COHORT>/<SAMPLE>.bismark.min5.max50.cov.gz
#-------------------------------------------------------------------------------

#SBATCH -t 00-12:00:00
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH --job-name=TMS.COVFILTER
#SBATCH --output=LOGS/%x_%j.out
#SBATCH --error=LOGS/%x_%j.err

set -euo pipefail

# === USER CONFIGURATION ===
PROJECT_ROOT="${TGC_PROJECT_ROOT:-/path/to/your/project}"
# ===========================
INPUT="${PROJECT_ROOT}/RESULTS/TMS/METHYLATION_CALLS"
OUTBASE="${PROJECT_ROOT}/RESULTS/TMS/METHYLATION_FILTERED"

MINCOV=5
MAXCOV=50

mkdir -p "${OUTBASE}/BREEDING" "${OUTBASE}/NATURAL"

echo "Starting coverage filtering: ${MINCOV} <= cov <= ${MAXCOV}"
echo "Input:  ${INPUT}"
echo "Output: ${OUTBASE}"
echo "----------------------------------"

shopt -s nullglob
files=( ${INPUT}/*/*.bismark.cov.gz )

if [[ ${#files[@]} -eq 0 ]]; then
  echo "ERROR: No .bismark.cov.gz files found under ${INPUT}/*/"
  exit 1
fi

for file in "${files[@]}"; do
  base=$(basename "$file" .bismark.cov.gz)   # e.g. P001_WA02_R1_p_bismark_bt2_pe
  sample="${base%%_R1_p_bismark_bt2_pe}"     # e.g. P001_WA02  (adjust if needed)

  if [[ "$sample" =~ ^P00[1-3]_ ]]; then
    cohort="BREEDING"
  elif [[ "$sample" =~ ^P00[4-8]_ ]]; then
    cohort="NATURAL"
  else
    echo "Skipping unknown cohort: ${sample} (from ${base})"
    continue
  fi

  outfile="${OUTBASE}/${cohort}/${base}.min${MINCOV}.max${MAXCOV}.cov.gz"

  echo "Processing ${base}  ->  ${cohort}"

  zcat "$file" \
    | awk -v lo="${MINCOV}" -v hi="${MAXCOV}" 'BEGIN{OFS="\t"} {cov=$5+$6; if(cov>=lo && cov<=hi) print $0}' \
    | gzip > "$outfile"
done

echo "Filtering completed."
exit 0
```

---

### `5b1.tgc.tms.methylkit.import.R`

```r
#!/usr/bin/env Rscript
############################################################
# TreeGeneClimate (TGC) — TMS
# Step 5b1: methylKit import of one sample (all three contexts)
#
# Reads a sorted, deduplicated Bismark BAM with processBismarkAln and saves
# one methylRaw object per cytosine context (CpG, CHG, CHH). Import settings
# are the methylKit defaults: per-sample minimum coverage 10 reads and
# minimum base quality 20. No other filtering is applied here; the SNP mask,
# coverage limits and uniting are applied in Steps 5b2 and 5b3.
#
# USAGE
#   Rscript 5b1.tgc.tms.methylkit.import.R <sample_bam>
#   (output folder can be changed with the environment variable TGC_RAW_DIR)
#   (run as a SLURM array by 5b1.tgc.tms.methylkit.import.sh)
#
# INPUT
#   DATA/TMS/MAPPED.FILES.TMS/<SAMPLE>_R1_p_bismark_bt2_pe.dedup.bam  (Step 3b2)
#
# OUTPUT
#   RESULTS/TMS/METHYLKIT_RAW/<SAMPLE>_<CONTEXT>.rds
#   RESULTS/TMS/METHYLKIT_RAW/COUNTS/<SAMPLE>.tsv   (sites per context)
############################################################

suppressPackageStartupMessages(library(methylKit))

# === USER CONFIGURATION ===
PROJECT_ROOT <- Sys.getenv("TGC_PROJECT_ROOT", unset = "/path/to/your/project")
MIN_COV  <- 10L   # methylKit default
MIN_QUAL <- 20L   # methylKit default
ASSEMBLY <- "Picab02"
# ===========================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1 || !file.exists(args[1])) stop("Usage: 5b1...R <sample_bam>")
bam <- args[1]
sample_id <- sub("_R1_p_bismark_bt2_pe.*\\.bam$", "", basename(bam))

OUT_DIR <- Sys.getenv("TGC_RAW_DIR", unset = file.path(PROJECT_ROOT, "RESULTS/TMS/METHYLKIT_RAW"))
dir.create(file.path(OUT_DIR, "COUNTS"), recursive = TRUE, showWarnings = FALSE)

counts <- data.frame(sample = sample_id, context = c("CpG", "CHG", "CHH"), n_sites = NA_integer_)
for (i in seq_len(nrow(counts))) {
  ctx <- counts$context[i]
  mr <- processBismarkAln(location = bam, sample.id = sample_id, assembly = ASSEMBLY,
                          read.context = ctx, mincov = MIN_COV, minqual = MIN_QUAL,
                          save.folder = NULL, verbose = 0)
  saveRDS(mr, file.path(OUT_DIR, sprintf("%s_%s.rds", sample_id, ctx)))
  counts$n_sites[i] <- nrow(mr)
  rm(mr); gc()
}
write.table(counts, file.path(OUT_DIR, "COUNTS", paste0(sample_id, ".tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sample_id, ":", paste(counts$context, counts$n_sites, collapse = "; "), "\n")
cat("methylKit", as.character(packageVersion("methylKit")), "\n")
```

---

### `5b1.tgc.tms.methylkit.import.sh`

```bash
#!/bin/bash
#SBATCH --job-name=TMS.5b1.IMPORT
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --array=0-601%60
#SBATCH --output=LOGS/5b1/%x_%A_%a.out
#SBATCH --error=LOGS/5b1/%x_%A_%a.err
#-------------------------------------------------------------------------------
# TreeGeneClimate (TGC) — TMS
# Step 5b1: methylKit import, one array task per sample (602 samples)
#
# Input:  DATA/TMS/MAPPED.FILES.TMS/<SAMPLE>_R1_p_bismark_bt2_pe.dedup.bam  (Step 3b2)
# Output: RESULTS/TMS/METHYLKIT_RAW/<SAMPLE>_<CONTEXT>.rds
#
# Create LOGS/5b1 before submitting.
# Software: R 4.5.2, methylKit 1.36.0
#-------------------------------------------------------------------------------
set -euo pipefail

# === USER CONFIGURATION ===
PROJECT_ROOT="${TGC_PROJECT_ROOT:-/path/to/your/project}"
# Software environment (edit for your system)
module purge 2>/dev/null || true
module load gcc/14.2.0 r/4.5.2 2>/dev/null || true
# export R_LIBS_USER=/path/to/your/Rlibs
export LC_ALL=C.UTF-8 LANG=C.UTF-8
# ===========================

BAM_DIR="${PROJECT_ROOT}/DATA/TMS/MAPPED.FILES.TMS"
BAM_SUFFIX="${TGC_BAM_SUFFIX:-_R1_p_bismark_bt2_pe.dedup.bam}"
i="${SLURM_ARRAY_TASK_ID:-0}"
sample=$(sed -n "$((i + 1))p" "${PROJECT_ROOT}/DATA/METADATA/samples_final_602.txt")
[[ -n "${sample}" ]] || { echo "Task ${i}: no sample"; exit 0; }

export TGC_PROJECT_ROOT="${PROJECT_ROOT}"
Rscript --vanilla "${PROJECT_ROOT}/SCRIPTS/TMS/5b1.tgc.tms.methylkit.import.R" "${BAM_DIR}/${sample}${BAM_SUFFIX}"
```

---

### `5b2.tgc.tms.methylkit.population.datasets.R`

```r
#!/usr/bin/env Rscript
############################################################
# TreeGeneClimate (TGC) — TMS
# Step 5b2: Population-level methylation datasets (one cohort x context)
#
# Used for methylation levels (6b), epigenomic PCA (7b), SVMP heatmaps (8b)
# and genomic-epigenomic concordance (11ab).
#
#   1) Remove sites with an ECS SNP (Step 4a3; alternative allele called at
#      least twice) at the cytosine or at the bases defining its context on the
#      same strand (CpG: C and G; CHG/CHH: C and the two following bases)
#   2) Per-sample coverage filter (methylKit filterByCoverage): at least 10
#      reads (as at import, Step 5b1) and at most the sample's 99.9th coverage
#      percentile, which removes pile-ups without penalising samples sequenced
#      more deeply (sequencing depth differs up to five-fold among trees)
#   3) unite(): keep sites covered in at least 80% (CpG), 70% (CHG) or
#      50% (CHH) of the cohort's samples (strand-specific, destrand = FALSE)
#   4) Methylated-fraction filter: keep sites whose pooled methylated fraction
#      across samples (sum numCs / sum coverage) is > 0.05
#
# USAGE
#   Rscript 5b2.tgc.tms.methylkit.population.datasets.R <BREEDING|NATURAL> <CpG|CHG|CHH>
#
# INPUTS
#   RESULTS/TMS/METHYLKIT_RAW/<SAMPLE>_<CONTEXT>.rds              (Step 5b1)
#   RESULTS/ECS/VARIANT.CALLING/MASK/tgc.ecs.cytosine.mask.snps.tsv.gz (Step 4a3)
#
# OUTPUTS (RESULTS/TMS/RANALYSIS/METHYLKIT_OBJECTS/)
#   methylBase_<cohort>_<ctx>_cov10_p99.9_mpg<k>.rds          (after unite)
#   methylBase_<cohort>_<ctx>_cov10_p99.9_mpg<k>_mef0.05.rds  (after step 4)
#   summary_<cohort>_<ctx>_population.tsv
############################################################

# === USER CONFIGURATION ===
PROJECT_ROOT <- Sys.getenv("TGC_PROJECT_ROOT", unset = "/path/to/your/project")
COV_LO <- 10L
COV_HI_PERC <- 99.9
MPG_FRAC <- c(CpG = 0.80, CHG = 0.70, CHH = 0.50)
METH_FRAC_MIN <- 0.05
# ===========================

source(file.path(PROJECT_ROOT, "SCRIPTS/TMS/tgc.tms.functions.R"))

args <- commandArgs(trailingOnly = TRUE)
COHORT <- toupper(args[1]); CTX <- args[2]
stopifnot(COHORT %in% c("BREEDING", "NATURAL"), CTX %in% names(MPG_FRAC))

RAW_DIR   <- file.path(PROJECT_ROOT, "RESULTS/TMS/METHYLKIT_RAW")
MASK_FILE <- file.path(PROJECT_ROOT, "RESULTS/ECS/VARIANT.CALLING/MASK/tgc.ecs.cytosine.mask.snps.tsv.gz")
OUT_DIR   <- file.path(PROJECT_ROOT, "RESULTS/TMS/RANALYSIS/METHYLKIT_OBJECTS")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
msg <- function(...) cat("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n", sep = "")

files <- raw_files(RAW_DIR, COHORT, CTX)
msg(COHORT, " / ", CTX, ": ", length(files), " samples")
mask <- read_snp_mask(MASK_FILE)
msg("Mask: ", nrow(mask), " SNP positions")

# 1) import + mask
n_raw <- n_masked <- integer(length(files))
objs <- vector("list", length(files))
for (i in seq_along(files)) {
  mr <- readRDS(files[i])
  n_raw[i] <- nrow(mr)
  mk <- mask_methylRaw(mr, mask, CTX)
  objs[[i]] <- mk$object
  n_masked[i] <- mk$n_removed
}
mrl <- methylRawList(objs, treatment = rep(0L, length(objs)))
rm(objs); gc()

# 2) coverage
mrl_cov <- filterByCoverage(mrl, lo.count = COV_LO, hi.perc = COV_HI_PERC)
n_cov <- vapply(as.list(mrl_cov), nrow, integer(1))
rm(mrl); gc()

# 3) unite
mpg <- as.integer(ceiling(MPG_FRAC[[CTX]] * length(files)))
mb <- unite(mrl_cov, destrand = FALSE, min.per.group = mpg)
rm(mrl_cov); gc()
msg("Sites after unite (min.per.group = ", mpg, "): ", nrow(mb))
stem <- sprintf("methylBase_%s_%s_cov%d_p%s_mpg%d", tolower(COHORT), tolower(CTX), COV_LO, COV_HI_PERC, mpg)
saveRDS(mb, file.path(OUT_DIR, paste0(stem, ".rds")))

# 4) pooled methylated fraction
pmf <- pooled_methylated_fraction(mb)
mb_mef <- mb[!is.na(pmf) & pmf > METH_FRAC_MIN, ]
msg("Sites after methylated-fraction filter (> ", METH_FRAC_MIN, "): ", nrow(mb_mef))
saveRDS(mb_mef, file.path(OUT_DIR, sprintf("%s_mef%.2f.rds", stem, METH_FRAC_MIN)))

summ <- data.frame(
  cohort = COHORT, context = CTX, n_samples = length(files),
  raw_sites_median = median(n_raw), masked_sites_median = median(n_masked),
  cov_sites_median = median(n_cov), min_per_group = mpg,
  n_sites_after_unite = nrow(mb), n_sites_after_mef = nrow(mb_mef)
)
write.table(summ, file.path(OUT_DIR, sprintf("summary_%s_%s_population.tsv", tolower(COHORT), tolower(CTX))),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(summ)
msg("Step 5b2 done")
```

---

### `5b2_5b3.tgc.tms.methylkit.datasets.sh`

```bash
#!/bin/bash
#SBATCH --job-name=TMS.5b2_5b3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --array=0-11
#SBATCH --output=LOGS/%x_%A_%a.out
#SBATCH --error=LOGS/%x_%A_%a.err
#-------------------------------------------------------------------------------
# TreeGeneClimate (TGC) — TMS
# Step 5b2_5b3: Steps 5b2 and 5b3 as one SLURM array (12 tasks):
#   tasks 0-5  : 5b2 population-level datasets  (cohort x context)
#   tasks 6-11 : 5b3 meQTL datasets             (cohort x context)
# Order within each block: BREEDING CpG, CHG, CHH, NATURAL CpG, CHG, CHH
#-------------------------------------------------------------------------------
set -euo pipefail

# === USER CONFIGURATION ===
export TGC_PROJECT_ROOT="${TGC_PROJECT_ROOT:-/path/to/your/project}"
module purge 2>/dev/null || true
module load gcc/14.2.0 r/4.5.2 2>/dev/null || true
export LC_ALL=C.UTF-8 LANG=C.UTF-8
# ===========================

i="${SLURM_ARRAY_TASK_ID:-0}"
cohorts=(BREEDING BREEDING BREEDING NATURAL NATURAL NATURAL)
contexts=(CpG CHG CHH CpG CHG CHH)
if [[ "${i}" -lt 6 ]]; then
  script="5b2.tgc.tms.methylkit.population.datasets.R"; k="${i}"
else
  script="5b3.tgc.tms.methylkit.meqtl.datasets.R"; k=$(( i - 6 ))
fi
Rscript --vanilla "${TGC_PROJECT_ROOT}/SCRIPTS/TMS/${script}" "${cohorts[$k]}" "${contexts[$k]}"
```

---

### `5b3.tgc.tms.methylkit.meqtl.datasets.R`

```r
#!/usr/bin/env Rscript
############################################################
# TreeGeneClimate (TGC) — TMS
# Step 5b3: Methylation datasets for cis-meQTL mapping (one cohort x context)
#
# Larger site sets than Step 5b2: no upper coverage limit and no
# methylated-fraction filter, and sites are united family by family
# (breeding) or stand by stand (natural), so that every group is represented.
#
#   1) Remove sites with an ECS SNP at the cytosine or at the bases defining
#      its context (as in Step 5b2)
#   2) Per-sample coverage >= 10 reads (applied at import, Step 5b1)
#   3) unite() with families (breeding) or stands (natural) as groups:
#      a site is kept if it is covered in at least m samples of EVERY group,
#      where m is the size of the smallest group (6 trees for breeding
#      families, 2 for natural stands); strand-specific (destrand = FALSE)
#
# USAGE
#   Rscript 5b3.tgc.tms.methylkit.meqtl.datasets.R <BREEDING|NATURAL> <CpG|CHG|CHH>
#
# INPUTS
#   RESULTS/TMS/METHYLKIT_RAW/<SAMPLE>_<CONTEXT>.rds              (Step 5b1)
#   RESULTS/ECS/VARIANT.CALLING/MASK/tgc.ecs.cytosine.mask.snps.tsv.gz (Step 4a3)
#   DATA/METADATA/breeding_sample2family.txt, natural_sample2pop.txt
#
# OUTPUTS (RESULTS/TMS/RANALYSIS/METHYLKIT_OBJECTS_MEQTL/)
#   methylBase_<cohort>_<ctx>_meqtl.rds   (attribute "groups": sample, group)
#   summary_<cohort>_<ctx>_meqtl.tsv
############################################################

# === USER CONFIGURATION ===
PROJECT_ROOT <- Sys.getenv("TGC_PROJECT_ROOT", unset = "/path/to/your/project")
# ===========================

source(file.path(PROJECT_ROOT, "SCRIPTS/TMS/tgc.tms.functions.R"))

args <- commandArgs(trailingOnly = TRUE)
COHORT <- toupper(args[1]); CTX <- args[2]
stopifnot(COHORT %in% c("BREEDING", "NATURAL"), CTX %in% c("CpG", "CHG", "CHH"))

RAW_DIR   <- file.path(PROJECT_ROOT, "RESULTS/TMS/METHYLKIT_RAW")
MASK_FILE <- file.path(PROJECT_ROOT, "RESULTS/ECS/VARIANT.CALLING/MASK/tgc.ecs.cytosine.mask.snps.tsv.gz")
MAP_FILE  <- file.path(PROJECT_ROOT, "DATA/METADATA",
                       if (COHORT == "BREEDING") "breeding_sample2family.txt" else "natural_sample2pop.txt")
OUT_DIR   <- file.path(PROJECT_ROOT, "RESULTS/TMS/RANALYSIS/METHYLKIT_OBJECTS_MEQTL")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
msg <- function(...) cat("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n", sep = "")

files <- raw_files(RAW_DIR, COHORT, CTX)
groups <- fread(MAP_FILE, header = FALSE, col.names = c("sample", "group"))
groups <- groups[match(names(files), sample)]
if (anyNA(groups$group)) stop("Samples without group: ", paste(names(files)[is.na(groups$group)], collapse = ", "))
group_index <- match(groups$group, sort(unique(groups$group)))
min_group <- min(table(group_index))
msg(COHORT, " / ", CTX, ": ", length(files), " samples in ", length(unique(group_index)),
    " groups; smallest group = ", min_group)

mask <- read_snp_mask(MASK_FILE)
n_raw <- n_masked <- integer(length(files))
objs <- vector("list", length(files))
for (i in seq_along(files)) {
  mr <- readRDS(files[i])
  n_raw[i] <- nrow(mr)
  mk <- mask_methylRaw(mr, mask, CTX)
  objs[[i]] <- mk$object
  n_masked[i] <- mk$n_removed
}
mrl <- methylRawList(objs, treatment = group_index)
rm(objs); gc()

mb <- unite(mrl, destrand = FALSE, min.per.group = as.integer(min_group))
rm(mrl); gc()
attr(mb, "groups") <- data.frame(sample = names(files), group = groups$group)
saveRDS(mb, file.path(OUT_DIR, sprintf("methylBase_%s_%s_meqtl.rds", tolower(COHORT), tolower(CTX))))
msg("Sites after unite: ", nrow(mb))

summ <- data.frame(
  cohort = COHORT, context = CTX, n_samples = length(files),
  n_groups = length(unique(group_index)), min_per_group = min_group,
  raw_sites_median = median(n_raw), masked_sites_median = median(n_masked),
  n_sites_after_unite = nrow(mb)
)
write.table(summ, file.path(OUT_DIR, sprintf("summary_%s_%s_meqtl.tsv", tolower(COHORT), tolower(CTX))),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(summ)
msg("Step 5b3 done")
```

---

### `6b.tgc.tms.methylation.levels.R`

```r
#!/usr/bin/env Rscript
############################################################
# TreeGeneClimate (TGC) — TMS
# Step 6b: Plot per-sample mean methylation (%) by group
#   - Breeding cohort: per family (17 families) -> panels a, c, e
#   - Natural cohort:  per stand  (25 stands)   -> panels b, d, f
#   - One plot per context (CpG/CHG/CHH)
#   - Kruskal-Wallis test per cohort x context and Dunn's post-hoc pairwise
#     comparisons (two-sided, Benjamini-Hochberg adjusted); letters from the
#     adjusted p-values (groups sharing a letter do not differ)
#
# STYLE:
#   - Panel titles: single bold lowercase letter (a–f), no extra text
#   - No legend
#   - Kruskal-Wallis p-value as subtitle; post-hoc letters above the boxes
#   - Posthoc letters angled; spaced above boxes to reduce overlap
#
# PANEL LAYOUT:
#   Top row:    a (breeding CpG)  |  b (natural CpG)
#   Middle row: c (breeding CHG)  |  d (natural CHG)
#   Bottom row: e (breeding CHH)  |  f (natural CHH)
#
# INPUT:
#   RDS from Step 5b2 (after unite, before the methylated-fraction filter):
#     RESULTS/TMS/RANALYSIS/METHYLKIT_OBJECTS/methylBase_<cohort>_<ctx>_cov10_p99.9_mpg<k>.rds
#   DATA/METADATA/breeding_sample2family.txt, natural_sample2pop.txt
#
# OUTPUT:
#   Figures:
#     RESULTS/TMS/RANALYSIS/FIGURES/METHYLATION_LEVELS/
#   Logs, posthoc tables and summary of the tests:
#     RESULTS/TMS/RANALYSIS/METHYLATION.LEVELS/
#       methylation_level_tests_summary.tsv, <cohort>_<ctx>_posthoc.tsv
#   Values shown in the panel:
#     RESULTS/TMS/RANALYSIS/FIGURES/METHYLATION_LEVELS/methylation_levels_panel_values.tsv
#
# R packages: methylKit, ggplot2, dunn.test, multcompView, patchwork
############################################################

## ---------------------------
## Packages
## ---------------------------
suppressPackageStartupMessages({
  library(methylKit)
  library(ggplot2)
  library(multcompView)
  library(patchwork)
  library(dunn.test)    # Dunn's post-hoc test
})

## ---------------------------
## Paths
## ---------------------------
# === USER CONFIGURATION ===
# Set TGC_PROJECT_ROOT as an environment variable, or edit the fallback path below
PROJECT_ROOT <- Sys.getenv("TGC_PROJECT_ROOT", unset = "/path/to/your/project")
# ===========================

rds_dir  <- file.path(PROJECT_ROOT, "RESULTS/TMS/RANALYSIS/METHYLKIT_OBJECTS")
fig_dir  <- file.path(PROJECT_ROOT, "RESULTS/TMS/RANALYSIS/FIGURES/METHYLATION_LEVELS")
log_dir  <- file.path(PROJECT_ROOT, "RESULTS/TMS/RANALYSIS/METHYLATION.LEVELS")

map_file_breeding <- file.path(PROJECT_ROOT, "DATA/METADATA/breeding_sample2family.txt")
map_file_natural  <- file.path(PROJECT_ROOT, "DATA/METADATA/natural_sample2pop.txt")

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)

log_file <- file.path(log_dir, "Step6b_methylation_level_stats_posthoc.log")

## ---------------------------
## Colour codes (families / stands)
## ---------------------------
colors.17 <- c(
  "Family_16" = "dodgerblue2", "Family_27" = "#E31A1C", "Family_32" = "green4",
  "Family_33" = "#6A3D9A", "Family_38" = "#FF7F00", "Family_39" = "black",
  "Family_40" = "gold1", "Family_41" = "skyblue2", "Family_42" = "#FB9A99",
  "Family_43" = "palegreen2", "Family_44" = "gray70", "Family_47" = "khaki2",
  "Family_48" = "orchid1", "Family_50" = "deeppink1", "Family_51" = "blue1",
  "Family_52" = "steelblue4", "Family_53" = "darkturquoise"
)

colors.25 <- c(
  "Asikkala"="dodgerblue2","Jämsä"="#E31A1C","Kauhajoki"="green4",
  "Koski"="#6A3D9A","Kuopio" ="#FF7F00","Laihia" ="black",
  "Lammi"="gold1","Leppävirta"="skyblue2","Loppi"="#FB9A99",
  "Luopioinen"="palegreen2","Mäntyharju"="#CAB2D6",
  "Marttila"="#FDBF6F","Miehikkälä"="gray80","Mikkeli"="khaki2",
  "Multia"="maroon","Muurame"="orchid1","Orivesi"="deeppink1",
  "Pälkäne"="blue1","Petäjävesi"="steelblue4","Punkaharju"="green1",
  "Punkalaidun"="yellow4","Puumala"="yellow3",
  "Rautalampi"="darkorange4","Savonlinna"="brown","Somero"="grey40"
)

## ---------------------------
## Helpers
## ---------------------------

# Read mapping (no header): Sample \t Group
read_map_noheader <- function(path) {
  x <- read.table(path, header = FALSE, sep = "", stringsAsFactors = FALSE)
  if (ncol(x) < 2) stop("Mapping file must have >=2 columns: sample_id <tab/space> group. File: ", path)
  colnames(x)[1:2] <- c("Sample", "Group")
  x[, c("Sample", "Group")]
}

# Find Step5b RDS for a cohort/context.
# If multiple match (e.g. reruns), it picks the most recently modified.
find_step5b_rds <- function(cohort, ctx, prefer_mef = TRUE) {
  cohort <- tolower(cohort)
  ctx <- tolower(ctx)
  
  # with the methylated-fraction filter
  pat_mef <- sprintf("^methylBase_%s_%s_cov10_p99\\.9_mpg[0-9]+_mef0\\.05\\.rds$", cohort, ctx)
  # cov+unite only
  pat_cov <- sprintf("^methylBase_%s_%s_cov10_p99\\.9_mpg[0-9]+\\.rds$", cohort, ctx)
  
  if (prefer_mef) {
    hits <- list.files(rds_dir, pattern = pat_mef, full.names = TRUE)
    if (length(hits) == 0) {
      # fall back to cov-only
      hits <- list.files(rds_dir, pattern = pat_cov, full.names = TRUE)
    }
  } else {
    hits <- list.files(rds_dir, pattern = pat_cov, full.names = TRUE)
    if (length(hits) == 0) {
      # fall back to the methylated-fraction-filtered file
      hits <- list.files(rds_dir, pattern = pat_mef, full.names = TRUE)
    }
  }
  
  if (length(hits) == 0) {
    stop("No Step5b RDS found for cohort=", cohort, " ctx=", ctx,
         " in: ", rds_dir, "\nExpected patterns: ", pat_cov, " OR ", pat_mef)
  }
  
  if (length(hits) > 1) {
    # pick most recent file
    mt <- file.info(hits)$mtime
    hits <- hits[order(mt, decreasing = TRUE)]
    warning("Multiple RDS matches for ", cohort, "/", ctx, ". Using newest:\n  ", hits[1])
  }
  
  hits[1]
}

# Compute per-sample mean methylation (%) over sites
methylbase_to_sample_means <- function(mb) {
  d <- getData(mb)
  
  numCs_cols <- grep("^numCs[0-9]+$", colnames(d), value = TRUE)
  numTs_cols <- grep("^numTs[0-9]+$", colnames(d), value = TRUE)
  
  if (length(numCs_cols) == 0 || length(numTs_cols) == 0 || length(numCs_cols) != length(numTs_cols)) {
    stop("Could not find matching numCs#/numTs# columns in methylBase.")
  }
  
  n <- length(numCs_cols)
  sample_ids <- mb@sample.ids
  if (length(sample_ids) != n) {
    warning("Length mismatch: sample.ids vs numCs columns. Using numCs column count.")
    sample_ids <- paste0("S", seq_len(n))
  }
  
  perc_mat <- vapply(seq_len(n), function(i) {
    numCs <- d[[numCs_cols[i]]]
    numTs <- d[[numTs_cols[i]]]
    perc <- 100 * (numCs / (numCs + numTs))
    perc[is.nan(perc)] <- NA_real_
    perc
  }, numeric(nrow(d)))
  
  colnames(perc_mat) <- sample_ids
  colMeans(perc_mat, na.rm = TRUE)
}

# Groups are compared with rank-based tests because per-sample mean methylation is
# not normally distributed within groups (Shapiro-Wilk test on the residuals from
# the group means, reported): Kruskal-Wallis test, then Dunn's test for all pairs
# of groups (two-sided p-values, Benjamini-Hochberg adjustment).
run_group_tests <- function(df) {
  df <- df[!is.na(df$Group) & !is.na(df$Methylation), ]
  df$Group <- droplevels(df$Group)

  resid <- df$Methylation - ave(df$Methylation, df$Group)
  sh_p  <- tryCatch(shapiro.test(resid)$p.value, error = function(e) NA_real_)
  kw    <- kruskal.test(Methylation ~ Group, data = df)
  kw_p  <- kw$p.value

  invisible(capture.output(
    dn <- dunn.test::dunn.test(df$Methylation, df$Group, method = "bh", altp = TRUE,
                               kw = FALSE, table = FALSE, list = FALSE)
  ))
  pairs <- data.frame(comparison = dn$comparisons, Z = dn$Z,
                      p = dn$altP, p_adj = dn$altP.adjusted)
  p_adj <- setNames(pairs$p_adj, gsub(" - ", "-", pairs$comparison, fixed = TRUE))
  letters <- multcompView::multcompLetters(p_adj, threshold = 0.05)$Letters

  list(
    method    = "Kruskal-Wallis + Dunn (BH)",
    shapiro_p = sh_p,
    global_p  = kw_p,
    kw_H      = unname(kw$statistic),
    kw_df     = unname(kw$parameter),
    posthoc   = pairs,
    letters   = letters
  )
}

# p-value formatting (only p-value in plot)
fmt_p <- function(p) {
  if (!is.finite(p)) return("p = NA")
  if (p < 1e-4) paste0("p = ", format(p, scientific = TRUE, digits = 2))
  else paste0("p = ", format(p, digits = 3))
}

# Save in all four formats independently (no conversions; each format rendered natively)
save_all_formats <- function(p, path_tiff, w_cm, h_cm) {
  ggsave(path_tiff,
         plot = p, width = w_cm, height = h_cm, units = "cm",
         dpi = 300, device = "tiff", compression = "lzw")
  ggsave(sub("\\.tiff$", ".pdf", path_tiff),
         plot = p, width = w_cm, height = h_cm, units = "cm",
         device = "pdf")
  ggsave(sub("\\.tiff$", ".eps", path_tiff),
         plot = p, width = w_cm, height = h_cm, units = "cm",
         device = cairo_ps)
  ggsave(sub("\\.tiff$", ".png", path_tiff),
         plot = p, width = w_cm, height = h_cm, units = "cm",
         dpi = 150, device = "png")
}

# Box plot of per-sample means by group
make_boxplot_clean <- function(df, color_vec, title_letter,
                               ylab = "Percent methylation",
                               letter_angle = 45,
                               letter_yfactor = 0.12) {
  
  present <- levels(df$Group)
  color_vec <- color_vec[names(color_vec) %in% present]
  
  ggplot(df, aes(x = Group, y = Methylation, fill = Group)) +
    geom_boxplot(outlier.shape = NA, width = 0.7) +
    geom_jitter(position = position_jitter(width = 0.2, seed = 1), size = 1, color = "black", alpha = 0.7) +
    scale_fill_manual(values = color_vec, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.03, 0.18))) +   # room for the letters
    coord_cartesian(clip = "off") +
    labs(title = title_letter, x = NULL, y = ylab) +
    theme_minimal(base_size = 18) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 13),
      axis.title.y = element_text(size = 16),
      plot.title  = element_text(size = 22, hjust = 0, vjust = 1, face = "bold"),
      legend.position = "none"
    )
}

# Analyze one cohort+context (build df, tests, plot, save fig + posthoc)
analyze_one_clean <- function(cohort, ctx, map_path, color_vec,
                              fig_prefix, title_letter,
                              prefer_mef = TRUE,
                              letter_angle = 45,
                              letter_yfactor = 0.12,
                              letter_size = 5) {
  
  rds_path <- find_step5b_rds(cohort, ctx, prefer_mef = prefer_mef)
  mb <- readRDS(rds_path)
  
  map <- read_map_noheader(map_path)
  
  sample_ids <- mb@sample.ids
  grp <- map$Group[match(sample_ids, map$Sample)]
  
  if (any(is.na(grp))) {
    missing <- sample_ids[is.na(grp)]
    warning(cohort, " / ", ctx, ": ", length(missing),
            " samples not found in mapping file. They will be dropped.\nExamples: ",
            paste(head(missing, 8), collapse = ", "))
  }
  
  means <- methylbase_to_sample_means(mb)
  
  df <- data.frame(
    Sample = names(means),
    Group = grp[match(names(means), sample_ids)],
    Methylation = as.numeric(means),
    stringsAsFactors = FALSE
  )
  df <- df[!is.na(df$Group), ]
  df$Group <- factor(df$Group, levels = unique(df$Group))
  
  tests <- run_group_tests(df)
  
  # Plot base
  p <- make_boxplot_clean(df, color_vec, title_letter,
                          ylab = "Percent methylation",
                          letter_angle = letter_angle,
                          letter_yfactor = letter_yfactor)
  
  # Only p-value annotation
  y_span <- diff(range(df$Methylation, na.rm = TRUE))
  if (!is.finite(y_span) || y_span == 0) y_span <- 1
  
  p <- p + labs(subtitle = fmt_p(tests$global_p)) +
    theme(plot.subtitle = element_text(size = 13, face = "italic", hjust = 0))
  
  # Posthoc letters angled + spaced
  letters_vec <- tests$letters
  letters_vec <- letters_vec[levels(df$Group)]
  
  box_max <- aggregate(Methylation ~ Group, data = df, max, na.rm = TRUE)
  box_max$Letter <- unname(letters_vec[as.character(box_max$Group)])
  box_max$y <- box_max$Methylation + letter_yfactor * y_span
  
  p <- p + geom_text(
    data = box_max,
    aes(x = Group, y = y, label = Letter),
    size = letter_size, fontface = 1, vjust = 0,
    angle = letter_angle
  )
  
  # Save figure (all four formats independently)
  out_fig <- file.path(fig_dir, sprintf("%s_%s_boxplot.tiff", fig_prefix, ctx))
  save_all_formats(p, out_fig, w_cm = 26, h_cm = 14)
  
  # Save the post-hoc table (Dunn's test, BH-adjusted)
  posthoc_tsv <- file.path(log_dir, sprintf("%s_%s_posthoc.tsv", fig_prefix, ctx))
  write.table(tests$posthoc, posthoc_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
  
  list(df = df, plot = p, tests = tests, fig = out_fig, posthoc = posthoc_tsv, rds = rds_path)
}

## ---------------------------
## RUN ALL 6
## ---------------------------
sink(log_file)
cat("TGC — TMS — Step 6b\n")
cat("Timestamp: ", format(Sys.time()), "\n\n", sep = "")

cat("RDS dir:  ", rds_dir, "\n", sep = "")
cat("Fig dir:  ", fig_dir, "\n", sep = "")
cat("Log dir:  ", log_dir, "\n\n", sep = "")

# Global methylation levels must be estimated from all covered loci (cov+unite),
# NOT from the methylated-fraction-filtered subset, which drops sites unmethylated in all trees
# and would inflate the means; that subset is used for PCA (7b), SVMPs (8b) and 11ab.
prefer_mef <- FALSE

# BREEDING -> panels a/c/e (left column, per context)
res_b_CpG <- analyze_one_clean("BREEDING", "CpG", map_file_breeding, colors.17, "BREEDING", "a",
                               prefer_mef = prefer_mef, letter_angle = 45, letter_yfactor = 0.12)
res_b_CHG <- analyze_one_clean("BREEDING", "CHG", map_file_breeding, colors.17, "BREEDING", "c",
                               prefer_mef = prefer_mef, letter_angle = 45, letter_yfactor = 0.12)
res_b_CHH <- analyze_one_clean("BREEDING", "CHH", map_file_breeding, colors.17, "BREEDING", "e",
                               prefer_mef = prefer_mef, letter_angle = 45, letter_yfactor = 0.12)

# NATURAL -> panels b/d/f (right column, per context)
res_n_CpG <- analyze_one_clean("NATURAL", "CpG", map_file_natural, colors.25, "NATURAL", "b",
                               prefer_mef = prefer_mef, letter_angle = 45, letter_yfactor = 0.12, letter_size = 4)
res_n_CHG <- analyze_one_clean("NATURAL", "CHG", map_file_natural, colors.25, "NATURAL", "d",
                               prefer_mef = prefer_mef, letter_angle = 45, letter_yfactor = 0.12, letter_size = 4)
res_n_CHH <- analyze_one_clean("NATURAL", "CHH", map_file_natural, colors.25, "NATURAL", "f",
                               prefer_mef = prefer_mef, letter_angle = 45, letter_yfactor = 0.12, letter_size = 4)

cat("\n\n====================\nSUMMARY\n====================\n")
summ_line <- function(res, cohort, ctx) {
  cat(cohort, " / ", ctx, "\n", sep = "")
  cat("  Method: ", res$tests$method, "\n", sep = "")
  cat("  Shapiro-Wilk p (residuals): ", format(res$tests$shapiro_p, digits = 4), "\n", sep = "")
  cat("  Kruskal-Wallis p: ", format(res$tests$global_p, digits = 8), "\n", sep = "")
  cat("  RDS: ", res$rds, "\n", sep = "")
  cat("  Figure: ", res$fig, "\n", sep = "")
  cat("  Posthoc TSV: ", res$posthoc, "\n\n", sep = "")
}
summ_line(res_b_CpG, "BREEDING", "CpG")
summ_line(res_b_CHG, "BREEDING", "CHG")
summ_line(res_b_CHH, "BREEDING", "CHH")
summ_line(res_n_CpG, "NATURAL",  "CpG")
summ_line(res_n_CHG, "NATURAL",  "CHG")
summ_line(res_n_CHH, "NATURAL",  "CHH")

sink()

# Combined 2x3 panel: left = breeding, right = natural; rows = CpG / CHG / CHH
panel_combined <-
  (res_b_CpG$plot | res_n_CpG$plot) /
  (res_b_CHG$plot | res_n_CHG$plot) /
  (res_b_CHH$plot | res_n_CHH$plot) +
  plot_layout(guides = "collect") &
  theme(legend.position = "none")

save_all_formats(
  panel_combined,
  file.path(fig_dir, "methylation_levels_panel.tiff"),
  w_cm = 48, h_cm = 42
)

# Summary table of the tests per cohort and context
summ_row <- function(res, cohort, ctx) {
  ph <- res$tests$posthoc
  data.frame(cohort = cohort, context = ctx, n_samples = nrow(res$df),
             n_groups = nlevels(res$df$Group),
             shapiro_p = res$tests$shapiro_p, kruskal_H = res$tests$kw_H,
             kruskal_df = res$tests$kw_df, kruskal_p = res$tests$global_p,
             dunn_sig_pairs = sum(ph$p_adj < 0.05), dunn_total_pairs = nrow(ph),
             mean_methylation = mean(res$df$Methylation))
}
summ <- rbind(summ_row(res_b_CpG, "BREEDING", "CpG"), summ_row(res_b_CHG, "BREEDING", "CHG"),
              summ_row(res_b_CHH, "BREEDING", "CHH"), summ_row(res_n_CpG, "NATURAL", "CpG"),
              summ_row(res_n_CHG, "NATURAL", "CHG"), summ_row(res_n_CHH, "NATURAL", "CHH"))
write.table(summ, file.path(log_dir, "methylation_level_tests_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(summ)

# Values shown in the panel: per-tree mean methylation, with the Dunn's test letters
sd_rows <- function(res, cohort, ctx, panel) {
  d <- res$df
  data.frame(panel = panel, cohort = cohort, context = ctx, tree = d$Sample,
             group = as.character(d$Group), percent_methylation = d$Methylation,
             dunn_letter = unname(res$tests$letters[as.character(d$Group)]),
             kruskal_p = res$tests$global_p)
}
write.table(rbind(sd_rows(res_b_CpG, "BREEDING", "CpG", "a"), sd_rows(res_n_CpG, "NATURAL", "CpG", "b"),
                  sd_rows(res_b_CHG, "BREEDING", "CHG", "c"), sd_rows(res_n_CHG, "NATURAL", "CHG", "d"),
                  sd_rows(res_b_CHH, "BREEDING", "CHH", "e"), sd_rows(res_n_CHH, "NATURAL", "CHH", "f")),
            file.path(fig_dir, "methylation_levels_panel_values.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

message("DONE.\nLog: ", log_file, "\nFigures: ", fig_dir, "\nPosthoc tables: ", log_dir)
```

---

### `7b.tgc.tms.pca.R`

```r
#!/usr/bin/env Rscript
############################################################
# TreeGeneClimate (TGC) — TMS
# Step 7b: Epigenomic PCA panel (breeding + natural cohorts)
#
# Save PCA as ONE panel (a-f) with shared legend per row:
#   top row = breeding (CpG, CHG, CHH) + legend at right of c)
#   bottom  = natural  (CpG, CHG, CHH) + legend at right of f)
# Panel letters are drawn inside plots (ggplot title).
# Epimarker IDs use genomic labels chr:start (chr:pos) instead of V123 / V1 etc.
#
# INPUT:
#   RESULTS/TMS/RANALYSIS/METHYLKIT_OBJECTS/
#     methylBase_<cohort>_<context>_cov10_p99.9_mpg*_mef0.05.rds   (Step 5b2)
#
# OUTPUT:
#   RESULTS/TMS/RANALYSIS/FIGURES/PCA/
#     pca_panel.tiff (+ pdf, eps, png)
#     pca_panel_values.tsv  (values shown in the panel)
############################################################

suppressPackageStartupMessages({
  library(methylKit)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(scales)
  library(grid)
  library(patchwork)
})

options(stringsAsFactors = FALSE)

# ==============================================================================
# 1) PATHS
# ==============================================================================
# === USER CONFIGURATION ===
# Set TGC_PROJECT_ROOT as an environment variable, or edit the fallback path below
PROJECT_ROOT <- Sys.getenv("TGC_PROJECT_ROOT",
  unset = "/path/to/your/project")
# ===========================
rds_dir <- file.path(PROJECT_ROOT, "RESULTS/TMS/RANALYSIS/METHYLKIT_OBJECTS")
fig_dir <- file.path(PROJECT_ROOT, "RESULTS/TMS/RANALYSIS/FIGURES/PCA")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

map_file_breeding <- file.path(PROJECT_ROOT, "DATA/METADATA/breeding_sample2family.txt")
map_file_natural  <- file.path(PROJECT_ROOT, "DATA/METADATA/natural_sample2pop.txt")

ensure_file <- function(p) if (!file.exists(p)) stop("Missing file: ", p, call. = FALSE)
ensure_file(map_file_breeding)
ensure_file(map_file_natural)

# ==============================================================================
# 2) PALETTES
# ==============================================================================
# Named vectors: consistent group-to-color mapping used across all PCA panels
colors.17 <- c(
  "Family_16"="dodgerblue2","Family_27"="#E31A1C","Family_32"="green4",
  "Family_33"="#6A3D9A","Family_38"="#FF7F00","Family_39"="black",
  "Family_40"="gold1","Family_41"="skyblue2","Family_42"="#FB9A99",
  "Family_43"="palegreen2","Family_44"="gray70","Family_47"="khaki2",
  "Family_48"="orchid1","Family_50"="deeppink1","Family_51"="blue1",
  "Family_52"="steelblue4","Family_53"="darkturquoise"
)

colors.25 <- c(
  "Asikkala"="dodgerblue2","Jämsä"="#E31A1C","Kauhajoki"="green4",
  "Koski"="#6A3D9A","Kuopio"="#FF7F00","Laihia"="black",
  "Lammi"="gold1","Leppävirta"="skyblue2","Loppi"="#FB9A99",
  "Luopioinen"="palegreen2","Mäntyharju"="#CAB2D6",
  "Marttila"="#FDBF6F","Miehikkälä"="gray80","Mikkeli"="khaki2",
  "Multia"="maroon","Muurame"="orchid1","Orivesi"="deeppink1",
  "Pälkäne"="blue1","Petäjävesi"="steelblue4","Punkaharju"="green1",
  "Punkalaidun"="yellow4","Puumala"="yellow3",
  "Rautalampi"="darkorange4","Savonlinna"="brown","Somero"="grey40"
)

# Restrict palette to groups actually present (avoids grey fallback entries in legend)
subset_palette <- function(pal_named, groups) {
  present <- unique(as.character(groups))
  pal <- pal_named[names(pal_named) %in% present]
  missing <- setdiff(present, names(pal))
  if (length(missing)) pal <- c(pal, setNames(rep("grey70", length(missing)), missing))
  pal
}

# ==============================================================================
# 3) HELPERS
# ==============================================================================

read_map_noheader <- function(path) {
  x <- read.table(path, header = FALSE, sep = "", stringsAsFactors = FALSE)
  if (ncol(x) < 2) stop("Mapping file must have >=2 columns: sample_id <tab/space> group. File: ", path)
  colnames(x)[1:2] <- c("Sample", "Group")
  x[, c("Sample", "Group")]
}

# Select the methylated-fraction-filtered RDS with the highest mpg (most sites)
pick_rds_mef <- function(cohort, ctx) {
  pat <- sprintf("^methylBase_%s_%s_cov10_p99\\.9_mpg[0-9]+_mef0\\.05\\.rds$",
                 tolower(cohort), tolower(ctx))
  files <- list.files(rds_dir, pattern = pat, full.names = TRUE)
  if (length(files) == 0) stop("No methylated-fraction-filtered RDS found for ", cohort, " / ", ctx, " under ", rds_dir)

  # extract the mpg number and pick the highest (largest site set)
  mpg_num <- as.integer(sub(".*_mpg([0-9]+)_mef0\\.05\\.rds$", "\\1", basename(files)))
  files[which.max(mpg_num)]
}

# Convert methylBase to methylation % matrix (sites x samples) + site info table
# Row names of the returned matrix are set to genomic "loc" IDs (chr:start-end)
# so that PCA loadings are interpretable rather than V1, V2, ...
methylbase_to_matrix <- function(mb) {
  d <- getData(mb)

  numCs_cols <- grep("^numCs[0-9]+$", colnames(d), value = TRUE)
  numTs_cols <- grep("^numTs[0-9]+$", colnames(d), value = TRUE)
  if (length(numCs_cols) == 0 || length(numTs_cols) == 0 || length(numCs_cols) != length(numTs_cols)) {
    stop("Could not find matching numCs#/numTs# columns in methylBase.")
  }

  sample_ids <- mb@sample.ids
  if (length(sample_ids) != length(numCs_cols)) {
    warning("sample.ids length mismatch; falling back to numCs columns count.")
    sample_ids <- paste0("S", seq_along(numCs_cols))
  }

  site_df <- d %>%
    as.data.frame() %>%
    dplyr::select(chr, start, end) %>%
    mutate(
      # "loc" is unique marker ID, same idea as ECS
      loc = paste0(chr, ":", start, "-", end),
      pos = start,
      chrpos = paste0(chr, ":", start)
    )

  perc_mat <- vapply(seq_along(numCs_cols), function(i) {
    numCs <- d[[numCs_cols[i]]]
    numTs <- d[[numTs_cols[i]]]
    p <- 100 * (numCs / (numCs + numTs))
    p[is.nan(p)] <- NA_real_
    p
  }, numeric(nrow(d)))

  rownames(perc_mat) <- site_df$loc
  colnames(perc_mat) <- sample_ids

  list(mat_sites_x_samples = perc_mat, site_df = site_df)
}

# panel label inside plot — uses ggplot title to match ECS PCA format
add_panel_label <- function(p, lab) {
  p + labs(title = lab) +
    theme(plot.title = element_text(hjust = 0, face = "bold"))
}

# context-aware NA thresholds:
# CHH is the most sparsely methylated context in conifers, so a more
# permissive missingness threshold is needed to retain enough sites.
max_na_frac_by_ctx <- function(ctx) {
  switch(ctx,
         "CpG" = 0.20,
         "CHG" = 0.30,
         "CHH" = 0.50,
         0.30
  )
}

# Remove high-missingness sites, then impute remaining NAs by site mean.
# Sites with near-zero variance after imputation are also removed (would
# inflate artificial PCs).
filter_and_impute_sites <- function(X_samples_x_sites,
                                    max_na_frac = 0.20,
                                    min_sd = 1e-8) {

  na_frac <- colMeans(is.na(X_samples_x_sites))
  keep1 <- na_frac <= max_na_frac
  X2 <- X_samples_x_sites[, keep1, drop = FALSE]
  if (ncol(X2) < 3) {
    stop("Too few sites after missingness filter: ", ncol(X2),
         " (max_na_frac=", max_na_frac, ")")
  }

  # impute remaining NAs by site mean
  for (j in seq_len(ncol(X2))) {
    idx <- is.na(X2[, j])
    if (any(idx)) {
      mu <- mean(X2[, j], na.rm = TRUE)
      if (!is.finite(mu)) mu <- 0
      X2[idx, j] <- mu
    }
  }

  # drop near-zero variance sites (uninformative for ordination)
  sds <- apply(X2, 2, sd)
  keep2 <- is.finite(sds) & (sds > min_sd)
  X3 <- X2[, keep2, drop = FALSE]
  if (ncol(X3) < 3) stop("Too few variable sites after SD filter: ", ncol(X3))

  list(X = X3)
}

run_pca <- function(X_samples_x_sites, max_na_frac = 0.20) {
  X2 <- filter_and_impute_sites(X_samples_x_sites, max_na_frac = max_na_frac)$X
  # Center and scale so that all sites contribute equally regardless of mean methylation level
  pr <- prcomp(X2, center = TRUE, scale. = TRUE)
  pve <- (pr$sdev^2) / sum(pr$sdev^2)   # proportion of variance explained per PC
  list(pr = pr, pve = pve, X_used = X2)
}

make_pca_scatter <- function(scores_df, pve, palette_named) {
  pal <- subset_palette(palette_named, scores_df$Group)

  ggplot(scores_df, aes(x = PC1, y = PC2, color = Group)) +
    geom_point(size = 1.8, alpha = 0.9) +
    scale_color_manual(values = pal) +
    labs(
      x = paste0("PC1 (", scales::percent(pve[1], accuracy = 0.1), ")"),
      y = paste0("PC2 (", scales::percent(pve[2], accuracy = 0.1), ")")
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(),
      legend.title = element_blank()
    )
}

# Save panel as TIFF (LZW compression), PDF, EPS and PNG
save_panel_tiff <- function(p, out_file, w_cm = 34, h_cm = 26, dpi = 600) {
  tiff(out_file, width = w_cm, height = h_cm, units = "cm", res = dpi, compression = "lzw")
  print(p); dev.off()
  ggsave(sub("\\.tiff$", ".pdf", out_file), plot = p, width = w_cm, height = h_cm, units = "cm")
  ggsave(sub("\\.tiff$", ".eps", out_file), plot = p, width = w_cm, height = h_cm, units = "cm", device = cairo_ps)
  ggsave(sub("\\.tiff$", ".png", out_file), plot = p, width = w_cm, height = h_cm, units = "cm", dpi = 150, device = "png")
  invisible(out_file)
}

# ==============================================================================
# 4) CORE RUNNER
# ==============================================================================
analyze_one <- function(cohort, ctx, map_path, palette_named) {

  rds_path <- pick_rds_mef(cohort, ctx)
  message("Using RDS: ", rds_path)
  mb <- readRDS(rds_path)

  map <- read_map_noheader(map_path)
  samp <- mb@sample.ids
  grp <- map$Group[match(samp, map$Sample)]

  if (any(is.na(grp))) {
    warning(cohort, " / ", ctx, ": ", sum(is.na(grp)), " samples missing in mapping; dropping them.")
  }
  keep_samp <- !is.na(grp)
  samp2 <- samp[keep_samp]
  grp2  <- grp[keep_samp]
  group_factor <- factor(as.character(grp2), levels = unique(as.character(grp2)))

  # Convert methylBase to percent matrix; row names = loc IDs for downstream traceability
  mm <- methylbase_to_matrix(mb)
  perc_sites_x_samples <- mm$mat_sites_x_samples
  site_df <- mm$site_df

  # subset samples (columns) to mapped samples
  perc_sites_x_samples <- perc_sites_x_samples[, samp2, drop = FALSE]

  # transpose to samples x sites; because rownames(perc_mat)=loc, colnames(X)=loc
  X <- t(perc_sites_x_samples)

  max_na_frac <- max_na_frac_by_ctx(ctx)

  # PCA on the full filtered site set
  pca <- run_pca(X, max_na_frac = max_na_frac)
  pve <- pca$pve
  scores <- as.data.frame(pca$pr$x[, 1:2, drop = FALSE])
  colnames(scores) <- c("PC1","PC2")
  scores$Group <- as.character(group_factor)
  p_pca <- make_pca_scatter(scores, pve, palette_named)

  list(
    cohort = cohort,
    ctx = ctx,
    pca_plot = p_pca,
    scores = data.frame(tree = rownames(X), group = scores$Group, PC1 = scores$PC1, PC2 = scores$PC2,
                        PC1_variance_pct = 100 * pve[1], PC2_variance_pct = 100 * pve[2])
  )
}

# ==============================================================================
# 5) RUN ALL 6
# ==============================================================================
res_b_cpg <- analyze_one("BREEDING", "CpG", map_file_breeding, colors.17)
res_b_chg <- analyze_one("BREEDING", "CHG", map_file_breeding, colors.17)
res_b_chh <- analyze_one("BREEDING", "CHH", map_file_breeding, colors.17)

res_n_cpg <- analyze_one("NATURAL", "CpG", map_file_natural, colors.25)
res_n_chg <- analyze_one("NATURAL", "CHG", map_file_natural, colors.25)
res_n_chh <- analyze_one("NATURAL", "CHH", map_file_natural, colors.25)

# ==============================================================================
# 6) PCA PANEL (a-f) with shared legend per row
# ==============================================================================
# Top row: breeding (A=CpG, B=CHG, C=CHH); bottom row: natural (D=CpG, E=CHG, F=CHH)
pca_a <- add_panel_label(res_b_cpg$pca_plot, "a")
pca_b <- add_panel_label(res_b_chg$pca_plot, "b")
pca_c <- add_panel_label(res_b_chh$pca_plot, "c")

pca_d <- add_panel_label(res_n_cpg$pca_plot, "d")
pca_e <- add_panel_label(res_n_chg$pca_plot, "e")
pca_f <- add_panel_label(res_n_chh$pca_plot, "f")

# Collect legend per row so breeding and natural families/populations have separate legends
pca_row1 <- (pca_a | pca_b | pca_c) + plot_layout(guides = "collect") & theme(legend.position = "right", legend.justification = "top")
pca_row2 <- (pca_d | pca_e | pca_f) + plot_layout(guides = "collect") & theme(legend.position = "right", legend.justification = "top")

pca_panel <- pca_row1 / pca_row2

out_pca_panel <- file.path(fig_dir, "pca_panel.tiff")
save_panel_tiff(pca_panel, out_pca_panel, w_cm = 34, h_cm = 26, dpi = 600)

# Values shown in the panel: PC1 and PC2 scores of every tree
res_all <- list(a = res_b_cpg, b = res_b_chg, c = res_b_chh, d = res_n_cpg, e = res_n_chg, f = res_n_chh)
write.table(do.call(rbind, lapply(names(res_all), function(k)
              cbind(panel = k, cohort = res_all[[k]]$cohort, context = res_all[[k]]$ctx, res_all[[k]]$scores))),
            file.path(fig_dir, "pca_panel_values.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ==============================================================================
# 7) FINISH
# ==============================================================================
cat("\nDONE Step 7b.\n\n")
cat("Saved: ", out_pca_panel, "  [PCA panel]\n", sep = "")
sessionInfo()
```

---

### `8b.tgc.tms.heatmaps.R`

```r
#!/usr/bin/env Rscript
############################################################
# TreeGeneClimate (TGC) — TMS
# Step 8b: Significantly variable methylated positions (SVMPs) and heatmaps per cohort,
#          on the population-level datasets (coverage and methylated-fraction filtered)
#
# INPUT (same as 7b):
#   RESULTS/TMS/RANALYSIS/METHYLKIT_OBJECTS/
#     methylBase_<cohort>_<context>_cov10_p99.9_mpg*_mef0.05.rds
#
# METHOD:
#   - For each cohort and each context (CpG/CHG/CHH):
#       * Convert methylBase to methylation % matrix (samples x loci)
#       * Filter loci by missingness (context-aware)
#       * Locus-by-locus Kruskal–Wallis test: %meth ~ Group (Family / Stand)
#       * BH adjust within context
#   - Select loci for heatmap by BH-adjusted p-value:
#       * Breeding: the 50 loci with the smallest adjusted p-value in each context (150)
#       * Natural: all loci with BH-adjusted p < 0.05 (SVMPs)
#   - Heatmap values are methylation % (0–100)
#   - Clustering requires a numeric matrix without too many NAs:
#       * Remaining NAs are imputed by locus mean (ONLY for heatmap/clustering)
#
# OUTPUT:
#   RESULTS/TMS/RANALYSIS/FIGURES/SVMP_HEATMAPS/
#     heatmap_breeding_top150_KW_BH.tiff
#     heatmap_natural_formal_SVMPs_KW_BH.tiff
#     svmp_heatmaps_panel.tiff  (combined A+B panel)
#
#   RESULTS/TMS/RANALYSIS/TABLES/heatmap_markers_8B/
#     TMS_8B_locus_tests_breeding_<ctx>.tsv
#     TMS_8B_locus_tests_natural_<ctx>.tsv
#     TMS_8B_selected_markers_breeding_top150.tsv
#     TMS_8B_selected_markers_natural_formal_SVMPs.tsv
#
#   RESULTS/TMS/RANALYSIS/FIGURES/SVMP_HEATMAPS/svmp_heatmaps_panel_values.tsv
#     values shown in the combined panel (tree x locus; values filled with the
#     locus mean are flagged)
############################################################

suppressPackageStartupMessages({
  library(methylKit)
  library(dplyr)
  library(tidyr)
  library(tibble)

  library(ComplexHeatmap)   # heatmaps with annotations
  library(circlize)         # colorRamp2 for continuous color scales
  library(viridisLite)      # colorblind-safe palettes
  library(grid)
})

options(stringsAsFactors = FALSE)
set.seed(1)   # reproducible hierarchical clustering (if ties in distance matrix)

# ==============================================================================
# 1) PATHS
# ==============================================================================
# === USER CONFIGURATION ===
# Set TGC_PROJECT_ROOT as an environment variable, or edit the fallback path below
PROJECT_ROOT <- Sys.getenv("TGC_PROJECT_ROOT",
  unset = "/path/to/your/project")
# ===========================

rds_dir <- file.path(PROJECT_ROOT, "RESULTS/TMS/RANALYSIS/METHYLKIT_OBJECTS")

fig_dir <- file.path(PROJECT_ROOT, "RESULTS/TMS/RANALYSIS/FIGURES/SVMP_HEATMAPS")
tab_dir <- file.path(PROJECT_ROOT, "RESULTS/TMS/RANALYSIS/TABLES/heatmap_markers_8B")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

map_file_breeding <- file.path(PROJECT_ROOT, "DATA/METADATA/breeding_sample2family.txt")
map_file_natural  <- file.path(PROJECT_ROOT, "DATA/METADATA/natural_sample2pop.txt")

ensure_file <- function(p) if (!file.exists(p)) stop("Missing file: ", p, call. = FALSE)
ensure_file(map_file_breeding)
ensure_file(map_file_natural)

# ==============================================================================
# 2) PALETTES
# ==============================================================================
# Named color vectors: consistent group-to-color mapping across heatmap annotations
colors.17 <- c(
  "Family_16"="dodgerblue2","Family_27"="#E31A1C","Family_32"="green4",
  "Family_33"="#6A3D9A","Family_38"="#FF7F00","Family_39"="black",
  "Family_40"="gold1","Family_41"="skyblue2","Family_42"="#FB9A99",
  "Family_43"="palegreen2","Family_44"="gray70","Family_47"="khaki2",
  "Family_48"="orchid1","Family_50"="deeppink1","Family_51"="blue1",
  "Family_52"="steelblue4","Family_53"="darkturquoise"
)

colors.25 <- c(
  "Asikkala"="dodgerblue2","Jämsä"="#E31A1C","Kauhajoki"="green4",
  "Koski"="#6A3D9A","Kuopio"="#FF7F00","Laihia"="black",
  "Lammi"="gold1","Leppävirta"="skyblue2","Loppi"="#FB9A99",
  "Luopioinen"="palegreen2","Mäntyharju"="#CAB2D6",
  "Marttila"="#FDBF6F","Miehikkälä"="gray80","Mikkeli"="khaki2",
  "Multia"="maroon","Muurame"="orchid1","Orivesi"="deeppink1",
  "Pälkäne"="blue1","Petäjävesi"="steelblue4","Punkaharju"="green1",
  "Punkalaidun"="yellow4","Puumala"="yellow3",
  "Rautalampi"="darkorange4","Savonlinna"="brown","Somero"="grey40"
)

subset_palette <- function(pal_named, groups) {
  present <- unique(as.character(groups))
  pal <- pal_named[names(pal_named) %in% present]
  missing <- setdiff(present, names(pal))
  if (length(missing)) pal <- c(pal, setNames(rep("grey70", length(missing)), missing))
  pal
}

# Methylation-context colours, shared by all figures (magenta, blue, turquoise; no red-green pair)
context_cols <- c(
  CpG = "#D81B60",
  CHG = "#1F4E9C",
  CHH = "#17BECF"
)

# ==============================================================================
# 3) HELPERS
# ==============================================================================

read_map_noheader <- function(path) {
  x <- read.table(path, header = FALSE, sep = "", stringsAsFactors = FALSE)
  if (ncol(x) < 2) stop("Mapping file must have >=2 columns: sample_id <tab/space> group. File: ", path)
  colnames(x)[1:2] <- c("Sample", "Group")
  x[, c("Sample", "Group")]
}

# Select the methylated-fraction-filtered RDS with the highest mpg (most sites)
pick_rds_mef <- function(cohort, ctx) {
  pat <- sprintf("^methylBase_%s_%s_cov10_p99\\.9_mpg[0-9]+_mef0\\.05\\.rds$",
                 tolower(cohort), tolower(ctx))
  files <- list.files(rds_dir, pattern = pat, full.names = TRUE)
  if (length(files) == 0) stop("No methylated-fraction-filtered RDS found for ", cohort, " / ", ctx, " under ", rds_dir)

  mpg_num <- as.integer(sub(".*_mpg([0-9]+)_mef0\\.05\\.rds$", "\\1", basename(files)))
  files[which.max(mpg_num)]
}

# Convert methylBase to methylation % matrix (samples x loci) + locus table
methylbase_to_percent_matrix <- function(mb) {
  d <- getData(mb)

  numCs_cols <- grep("^numCs[0-9]+$", colnames(d), value = TRUE)
  numTs_cols <- grep("^numTs[0-9]+$", colnames(d), value = TRUE)
  if (length(numCs_cols) == 0 || length(numTs_cols) == 0 || length(numCs_cols) != length(numTs_cols)) {
    stop("Could not find matching numCs#/numTs# columns in methylBase.")
  }

  sample_ids <- mb@sample.ids
  if (length(sample_ids) != length(numCs_cols)) {
    warning("sample.ids length mismatch; falling back to numCs columns count.")
    sample_ids <- paste0("S", seq_along(numCs_cols))
  }

  # percent methylation per locus per sample (loci x samples)
  perc_mat_loci_x_samples <- vapply(seq_along(numCs_cols), function(i) {
    numCs <- d[[numCs_cols[i]]]
    numTs <- d[[numTs_cols[i]]]
    p <- 100 * (numCs / (numCs + numTs))
    p[is.nan(p)] <- NA_real_
    p
  }, numeric(nrow(d)))

  colnames(perc_mat_loci_x_samples) <- sample_ids

  loci <- d %>%
    as.data.frame() %>%
    dplyr::select(chr, start, end) %>%
    mutate(loc = paste0(chr, ":", start, "-", end))

  # return samples x loci (transpose for downstream convenience)
  list(
    X = t(perc_mat_loci_x_samples),
    loci = loci
  )
}

# Context-aware missingness thresholds aligned with mpg (same as 7b logic)
# CHH has inherently lower methylation and more zero-coverage sites in conifers
max_na_frac_by_ctx <- function(ctx) {
  switch(ctx,
         "CpG" = 0.20,
         "CHG" = 0.30,
         "CHH" = 0.50,
         0.30
  )
}

# Filter loci by missingness; return filtered X and loci
filter_loci_by_missingness <- function(X, loci, max_na_frac) {
  na_frac <- colMeans(is.na(X))
  keep <- na_frac <= max_na_frac
  X2 <- X[, keep, drop = FALSE]
  loci2 <- loci[keep, , drop = FALSE]
  list(X = X2, loci = loci2, na_frac = na_frac[keep], kept = keep)
}

# Impute remaining NAs by locus mean (ONLY for heatmap/clustering)
# Imputation is intentionally withheld from the statistical tests (KW below)
# to avoid inflating significance; it is applied only to enable hierarchical clustering.
impute_locus_mean <- function(X) {
  X2 <- X
  for (j in seq_len(ncol(X2))) {
    idx <- is.na(X2[, j])
    if (any(idx)) {
      mu <- mean(X2[, j], na.rm = TRUE)
      if (!is.finite(mu)) mu <- 0
      X2[idx, j] <- mu
    }
  }
  X2
}

# Kruskal–Wallis per locus + epsilon-squared effect size
# KW is used because per-locus methylation % is bounded (0–100)
# and often bimodal; distributional assumptions of ANOVA are violated at individual loci.
kw_per_locus <- function(x, g) {
  ok <- is.finite(x) & !is.na(g)
  x <- x[ok]
  g <- g[ok]
  if (length(x) < 3 || length(unique(g)) < 2) {
    return(c(p = NA_real_, H = NA_real_, eps2 = NA_real_, n = length(x), k = length(unique(g))))
  }
  kt <- kruskal.test(x ~ g)
  H <- unname(kt$statistic)
  n <- length(x)
  k <- length(unique(g))
  # epsilon-squared for KW: (H - k + 1) / (n - k)
  eps2 <- (H - k + 1) / (n - k)
  c(p = unname(kt$p.value), H = H, eps2 = eps2, n = n, k = k)
}

# Balanced selection: EXACT equal per context when possible (e.g., top_n=150 => 50/50/50)
# This ensures all three methylation contexts are represented in the heatmap.
# NOTE: base indexing only (robust across dplyr versions).
select_top_loci_balanced <- function(df_all, top_n = 150L) {
  df_all <- df_all %>%
    dplyr::filter(is.finite(padj)) %>%
    dplyr::arrange(padj, dplyr::desc(eps2))
  if (nrow(df_all) == 0) return(df_all)

  ctxs <- intersect(c("CpG","CHG","CHH"), unique(df_all$Context))
  if (length(ctxs) == 0) return(df_all[0, ])

  if (top_n %% length(ctxs) != 0) {
    stop("top_n must be divisible by number of contexts. top_n=", top_n, " contexts=", length(ctxs))
  }
  per_ctx <- as.integer(top_n / length(ctxs))

  picked_list <- lapply(ctxs, function(ctx) {
    tmp <- df_all[df_all$Context == ctx, , drop = FALSE]
    if (nrow(tmp) == 0) return(tmp)
    tmp[seq_len(min(per_ctx, nrow(tmp))), , drop = FALSE]
  })
  picked <- dplyr::bind_rows(picked_list)

  # If any context has too few loci, fill remaining with best overall remaining
  target_n <- min(top_n, nrow(df_all))
  if (nrow(picked) < target_n) {
    already <- picked$marker_id
    fill_n <- target_n - nrow(picked)
    extra <- df_all[!df_all$marker_id %in% already, , drop = FALSE]
    if (nrow(extra) > 0) {
      extra <- extra[seq_len(min(fill_n, nrow(extra))), , drop = FALSE]
      picked <- dplyr::bind_rows(picked, extra)
    }
  }

  picked %>% dplyr::arrange(Context, padj, dplyr::desc(eps2))
}

# Overall selection: top_n by padj (then eps2), regardless of context
# NOTE: base indexing only (robust across dplyr versions).
select_top_loci_overall <- function(df_all, top_n = 150L) {
  df2 <- df_all %>%
    dplyr::filter(is.finite(padj)) %>%
    dplyr::arrange(padj, dplyr::desc(eps2))
  if (nrow(df2) == 0) return(df2)
  df2[seq_len(min(top_n, nrow(df2))), , drop = FALSE]
}

summarize_sig_counts <- function(df) {
  tibble(
    n_tested = sum(is.finite(df$p)),
    p_lt_0.05  = sum(df$p < 0.05,  na.rm = TRUE),
    p_lt_0.001 = sum(df$p < 0.001, na.rm = TRUE),
    padj_lt_0.05  = sum(df$padj < 0.05,  na.rm = TRUE),
    padj_lt_0.001 = sum(df$padj < 0.001, na.rm = TRUE)
  )
}

# ---- Heatmap builder (returns a Heatmap object; no device writing here)
# Rows = samples, columns = loci; both clustered with ward.D2 (minimises total
# within-cluster variance, suitable for methylation data which can form sharp groups).
build_heatmap_object <- function(mat_samples_x_markers,
                                 group_vec, group_palette,
                                 context_vec,
                                 group_legend_title  = "Group",
                                 show_legends        = TRUE,
                                 cluster_cols        = TRUE,
                                 show_annot_names    = TRUE,
                                 heatmap_name        = "%meth") {

  # align group vector to matrix rows (named factor indexed by sample ID)
  grp <- as.character(group_vec[rownames(mat_samples_x_markers)])
  grp[is.na(grp)] <- "Unknown"
  pal_grp <- subset_palette(group_palette, grp)

  # Create row annotation with the desired legend title (robust; no slot hacking)
  ha_row <- do.call(
    rowAnnotation,
    c(
      setNames(list(grp), group_legend_title),
      list(
        col                  = setNames(list(pal_grp), group_legend_title),
        show_annotation_name = show_annot_names,
        annotation_name_gp   = gpar(fontsize = 10)
      )
    )
  )

  # Subset context palette to only contexts present — avoids phantom legend entries
  present_ctxs <- intersect(c("CpG","CHG","CHH"), as.character(context_vec))
  ctx_col_sub  <- context_cols[present_ctxs]

  # Top annotation bar: cytosine context of each locus column
  ha_top <- HeatmapAnnotation(
    `Methylation context` = context_vec,
    col                   = list(`Methylation context` = ctx_col_sub),
    annotation_name_gp    = gpar(fontsize = 10),
    show_annotation_name  = show_annot_names
  )

  # 0–100 methylation scale (colorblind-friendly: cividis)
  col_fun <- circlize::colorRamp2(c(0, 50, 100), viridisLite::viridis(3, option = "cividis"))

  Heatmap(
    as.matrix(mat_samples_x_markers),
    name = heatmap_name,
    col = col_fun,
    cluster_rows             = TRUE,
    cluster_columns          = cluster_cols,
    clustering_method_rows   = "ward.D2",
    clustering_method_columns = "ward.D2",
    show_row_names  = FALSE,    # sample labels omitted for visual clarity
    show_column_names = FALSE,  # locus IDs too dense to display per column
    top_annotation = ha_top,
    left_annotation = ha_row,
    show_heatmap_legend = show_legends,
    heatmap_legend_param = list(
      title = NULL,
      at = c(0, 25, 50, 75, 100),
      labels = c("0", "25", "50", "75", "100"),
      labels_gp = gpar(fontsize = 8),
      grid_height = unit(3.2, "cm"),
      grid_width  = unit(0.38, "cm")
    )
  )
}

# ---- Save a single heatmap TIFF (20 x 20 cm by default)
# ComplexHeatmap requires a draw() call inside the graphics device;
# ggsave is not compatible, so each format is opened and closed manually.
save_heatmap_tiff <- function(mat_samples_x_markers,
                              group_vec, group_palette,
                              context_vec,
                              group_legend_title = "Group",
                              out_file, title_label = "a",
                              cluster_cols      = TRUE,
                              show_annot_names  = TRUE,
                              w_cm = 20, h_cm = 20, dpi = 600) {

  ht <- build_heatmap_object(
    mat_samples_x_markers = mat_samples_x_markers,
    group_vec             = group_vec,
    group_palette         = group_palette,
    context_vec           = context_vec,
    group_legend_title    = group_legend_title,
    show_legends          = TRUE,
    cluster_cols          = cluster_cols,
    show_annot_names      = show_annot_names,
    heatmap_name          = "%meth"
  )

  draw_heatmap <- function() {
    grid.newpage()
    draw(
      ht,
      heatmap_legend_side = "right",
      annotation_legend_side = "right",
      merge_legend = TRUE,        # combine row and top annotation legends
      show_annotation_legend = TRUE
    )
    # Panel letter positioned at top-left inside the graphics area
    grid.text(title_label,
              x = unit(0.35, "cm"),
              y = unit(1, "npc") - unit(0.35, "cm"),
              just = c("left","top"),
              gp = gpar(fontsize = 16, fontface = "bold"))
  }
  tiff(out_file, width = w_cm, height = h_cm, units = "cm", res = dpi, compression = "lzw")
  draw_heatmap(); dev.off()
  pdf(sub("\\.tiff$", ".pdf", out_file), width = w_cm / 2.54, height = h_cm / 2.54)
  draw_heatmap(); dev.off()
  cairo_ps(sub("\\.tiff$", ".eps", out_file), width = w_cm / 2.54, height = h_cm / 2.54)
  draw_heatmap(); dev.off()
  png(sub("\\.tiff$", ".png", out_file), width = w_cm, height = h_cm, units = "cm", res = dpi)
  draw_heatmap(); dev.off()

  invisible(out_file)
}

# ==============================================================================
# 4) CORE: process one cohort (combines CpG/CHG/CHH into one heatmap)
# ==============================================================================

process_cohort_for_heatmap <- function(cohort, map_path, palette_named,
                                       top_n = 150L) {

  map <- read_map_noheader(map_path)

  log_file <- file.path(tab_dir, sprintf("TMS_8B_LOG_%s.txt", tolower(cohort)))
  log_con <- file(log_file, open = "wt")
  on.exit(close(log_con), add = TRUE)

  writeLines(paste0("TMS Step 8b log — cohort: ", cohort), con = log_con)
  writeLines(paste0("Selection: Kruskal–Wallis per locus, BH correction within context."), con = log_con)
  writeLines(paste0("TopN for heatmap: ", top_n), con = log_con)
  writeLines("", con = log_con)

  contexts <- c("CpG","CHG","CHH")
  res_list <- list()
  test_tables <- list()

  for (ctx in contexts) {
    rds_path <- pick_rds_mef(cohort, ctx)
    message("Using RDS: ", rds_path)
    mb <- readRDS(rds_path)

    # group vector aligned to methylBase samples
    samp <- mb@sample.ids
    grp <- map$Group[match(samp, map$Sample)]

    # drop samples not in mapping
    keep_samp <- !is.na(grp)
    samp2 <- samp[keep_samp]
    grp2  <- grp[keep_samp]
    group_factor <- factor(as.character(grp2), levels = unique(as.character(grp2)))
    names(group_factor) <- samp2   # named so it can be indexed by sample ID later

    # methylation percent matrix
    mm <- methylbase_to_percent_matrix(mb)
    X <- mm$X
    loci <- mm$loci

    # subset to mapped samples
    X <- X[samp2, , drop = FALSE]

    # missingness filter (no imputation yet — KW runs on observed values only)
    max_na <- max_na_frac_by_ctx(ctx)
    flt <- filter_loci_by_missingness(X, loci, max_na_frac = max_na)
    Xf <- flt$X
    locif <- flt$loci

    if (ncol(Xf) < 5) {
      warning(cohort, "/", ctx, ": too few loci after missingness filter (", ncol(Xf), "). Skipping.")
      next
    }

    # locus-by-locus KW tests on filtered matrix (NO imputation for tests)
    kw_mat <- t(vapply(seq_len(ncol(Xf)), function(j) kw_per_locus(Xf[, j], group_factor),
                       numeric(5)))
    kw_df <- as.data.frame(kw_mat)
    kw_df$Context <- ctx
    kw_df$chr   <- locif$chr
    kw_df$start <- locif$start
    kw_df$end   <- locif$end
    kw_df$loc   <- locif$loc

    # Prefix context to marker_id so IDs remain unique when combined across contexts
    kw_df$marker_id <- paste0(ctx, "|", kw_df$loc)

    # BH correction applied within context (not across all three combined)
    # to control FDR independently per context given their different power profiles
    kw_df$padj <- p.adjust(kw_df$p, method = "BH")

    cnt <- summarize_sig_counts(kw_df)
    writeLines(paste0("Context: ", ctx), con = log_con)
    writeLines(paste0("  Loci tested: ", cnt$n_tested), con = log_con)
    writeLines(paste0("  Before adj:  p<0.05=", cnt$p_lt_0.05, " | p<0.001=", cnt$p_lt_0.001), con = log_con)
    writeLines(paste0("  After  BH:   q<0.05=", cnt$padj_lt_0.05, " | q<0.001=", cnt$padj_lt_0.001), con = log_con)
    writeLines("", con = log_con)

    # save per-context table
    out_ctx_tsv <- file.path(tab_dir, sprintf("TMS_8B_locus_tests_%s_%s.tsv",
                                              tolower(cohort), tolower(ctx)))
    kw_df_out <- kw_df %>%
      dplyr::select(marker_id, Context, chr, start, end, loc, p, padj, H, eps2, n, k)
    write.table(kw_df_out, out_ctx_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

    test_tables[[ctx]] <- kw_df_out

    # Impute NAs by locus mean for clustering only (not used in KW tests above)
    Xh <- impute_locus_mean(Xf)
    # Prefix column names with context to avoid ambiguity when matrices are combined
    colnames(Xh) <- paste0(ctx, "|", locif$loc)

    na_mask <- is.na(Xf)
    colnames(na_mask) <- colnames(Xh)

    res_list[[ctx]] <- list(
      X_heat = Xh,
      na_mask = na_mask,
      group = group_factor,
      loci = locif,
      rds = rds_path,
      max_na = max_na
    )

    message(sprintf("%s/%s: samples=%d | loci_before=%d | loci_after_missing=%d (max_na=%.2f)",
                    cohort, ctx, nrow(Xf), ncol(X), ncol(Xf), max_na))
  }

  if (length(res_list) == 0) stop("No contexts produced usable data for cohort: ", cohort)

  # Intersect sample sets across contexts: only samples present in all three matrices enter the heatmap
  samp_common <- Reduce(intersect, lapply(res_list, function(z) rownames(z$X_heat)))
  if (length(samp_common) < 10) stop("Too few common samples across contexts for ", cohort)

  # Concatenate per-context matrices column-wise to form a single combined matrix
  X_all <- NULL
  NA_all <- NULL
  df_all <- NULL

  for (ctx in names(res_list)) {
    Xc <- res_list[[ctx]]$X_heat[samp_common, , drop = FALSE]
    X_all <- if (is.null(X_all)) Xc else cbind(X_all, Xc)
    Nc <- res_list[[ctx]]$na_mask[samp_common, , drop = FALSE]
    NA_all <- if (is.null(NA_all)) Nc else cbind(NA_all, Nc)

    dfc <- test_tables[[ctx]]
    df_all <- if (is.null(df_all)) dfc else bind_rows(df_all, dfc)
  }

  # ---- Selection A: BALANCED (e.g., top150 => 50/50/50 when possible)
  selected_bal <- select_top_loci_balanced(df_all, top_n = top_n)
  out_sel_bal_tsv <- file.path(tab_dir, sprintf("TMS_8B_selected_markers_%s_top%d_BALANCED.tsv",
                                                tolower(cohort), top_n))
  write.table(selected_bal, out_sel_bal_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

  # ---- Selection B: OVERALL top_n by padj regardless of context
  selected_all <- select_top_loci_overall(df_all, top_n = top_n)
  out_sel_all_tsv <- file.path(tab_dir, sprintf("TMS_8B_selected_markers_%s_top%d_OVERALL.tsv",
                                                tolower(cohort), top_n))
  write.table(selected_all, out_sel_all_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

  # ---- Log selection summaries
  writeLines("Selected markers — BALANCED:", con = log_con)
  writeLines(paste(capture.output(print(selected_bal %>% dplyr::count(Context, name = "n_selected"))), collapse = "\n"), con = log_con)
  writeLines("", con = log_con)
  writeLines("Top 10 BALANCED:", con = log_con)
  writeLines(paste(capture.output(print(utils::head(selected_bal, 10))), collapse = "\n"), con = log_con)
  writeLines("", con = log_con)

  writeLines("Selected markers — OVERALL:", con = log_con)
  writeLines(paste(capture.output(print(selected_all %>% dplyr::count(Context, name = "n_selected"))), collapse = "\n"), con = log_con)
  writeLines("", con = log_con)
  writeLines("Top 10 OVERALL:", con = log_con)
  writeLines(paste(capture.output(print(utils::head(selected_all, 10))), collapse = "\n"), con = log_con)
  writeLines("", con = log_con)

  # ---- Build matrices for both selections by indexing combined matrix columns
  make_X_sel <- function(selected_df) {
    sel_ids <- selected_df$marker_id
    sel_ids <- sel_ids[sel_ids %in% colnames(X_all)]
    if (length(sel_ids) < 10) stop("Too few selected markers present in matrix for ", cohort, ": ", length(sel_ids))
    X_all[, sel_ids, drop = FALSE]
  }

  X_sel_bal <- make_X_sel(selected_bal)
  X_sel_all <- make_X_sel(selected_all)

  # Extract context label from the "ctx|loc" column name prefix for annotation bar
  ctx_vec_bal <- factor(sub("\\|.*$", "", colnames(X_sel_bal)), levels = c("CpG","CHG","CHH"))
  ctx_vec_all <- factor(sub("\\|.*$", "", colnames(X_sel_all)), levels = c("CpG","CHG","CHH"))

  # group vector from any context (named by sample)
  group_vec <- res_list[[1]]$group
  group_vec <- group_vec[samp_common]

  list(
    cohort = cohort,

    X_sel_bal = X_sel_bal,
    ctx_vec_bal = ctx_vec_bal,
    selected_bal_tsv = out_sel_bal_tsv,

    X_sel_all = X_sel_all,
    ctx_vec_all = ctx_vec_all,
    selected_all_tsv = out_sel_all_tsv,

    group_vec = group_vec,
    log_file  = log_file,
    df_all    = df_all,
    X_all     = X_all,
    NA_all    = NA_all
  )
}

# ==============================================================================
# 5) RUN HEATMAPS (per cohort)
# ==============================================================================
topN <- 150L   # number of loci shown per heatmap (50 per context when balanced)

# ---- BREEDING (panel a: top-150 balanced)
res_b <- process_cohort_for_heatmap("BREEDING", map_file_breeding, colors.17, top_n = topN)

out_b_bal <- file.path(fig_dir, sprintf("heatmap_breeding_top%d_KW_BH.tiff", topN))

save_heatmap_tiff(res_b$X_sel_bal, res_b$group_vec, colors.17, res_b$ctx_vec_bal,
                  group_legend_title = "Family",
                  out_file = out_b_bal, title_label = "a",
                  w_cm = 24, h_cm = 20, dpi = 600)

# ---- NATURAL (locus tests only; panel b uses the formal-SVMP subset below)
res_n <- process_cohort_for_heatmap("NATURAL", map_file_natural, colors.25, top_n = topN)

# ==============================================================================
# 6) NATURAL: FORMAL SVMPs ONLY (padj < 0.05) + COMBINED PANEL (a + b)
# ==============================================================================
# SVMPs (significantly variable methylation positions) = loci with significant
# inter-population variation in methylation level after BH correction

sig_df <- res_n$df_all[!is.na(res_n$df_all$padj) & res_n$df_all$padj < 0.05, ]
cat(sprintf("\nNatural formal SVMPs (padj < 0.05): %d loci\n", nrow(sig_df)))
if (nrow(sig_df) > 0) {
  ctx_counts <- table(sig_df$Context)
  cat(paste(names(ctx_counts), ctx_counts, sep = ": ", collapse = " | "), "\n")
}

if (nrow(sig_df) >= 2) {

  # Keep only IDs present in the full matrix; order by context then padj
  sig_ids <- sig_df$marker_id[sig_df$marker_id %in% colnames(res_n$X_all)]
  sig_df2 <- sig_df[match(sig_ids, sig_df$marker_id), ]
  ord     <- order(match(sig_df2$Context, c("CpG","CHG","CHH")), sig_df2$padj)
  sig_ids <- sig_ids[ord]
  sig_df2 <- sig_df2[ord, ]

  X_sig       <- res_n$X_all[, sig_ids, drop = FALSE]
  ctx_vec_sig <- factor(sub("\\|.*$", "", colnames(X_sig)), levels = c("CpG","CHG","CHH"))

  # Save the selected-marker table
  write.table(sig_df2,
              file.path(tab_dir, "TMS_8B_selected_markers_natural_formal_SVMPs.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  # ---- Panel b: natural formal SVMPs standalone
  out_n_sig <- file.path(fig_dir, "heatmap_natural_formal_SVMPs_KW_BH.tiff")
  save_heatmap_tiff(X_sig, res_n$group_vec, colors.25, ctx_vec_sig,
                    group_legend_title = "Natural stand",
                    out_file         = out_n_sig,
                    title_label      = "b",
                    cluster_cols     = TRUE,
                    show_annot_names = FALSE,   # annotation names omitted for the narrow panel
                    w_cm = 16, h_cm = 20, dpi = 600)

  # ---- Combined panel (a + b)
  # Breeding balanced (150 markers) gets 70% of width;
  # Natural formal SVMPs gets 30% — extra width ensures "Methylation context" legend fits.
  panel_w_cm <- 32
  panel_h_cm <- 20
  frac_a     <- 0.70
  label_gp   <- gpar(fontsize = 16, fontface = "bold")

  # Build heatmap objects for the combined panel (legends merged per sub-panel)
  ht_a_panel <- build_heatmap_object(
    mat_samples_x_markers = res_b$X_sel_bal,
    group_vec             = res_b$group_vec,
    group_palette         = colors.17,
    context_vec           = res_b$ctx_vec_bal,
    group_legend_title    = "Family",
    show_legends          = TRUE,
    cluster_cols          = TRUE,
    show_annot_names      = FALSE   # no sidebar / topbar labels in the plot
  )
  ht_b_panel <- build_heatmap_object(
    mat_samples_x_markers = X_sig,
    group_vec             = res_n$group_vec,
    group_palette         = colors.25,
    context_vec           = ctx_vec_sig,
    group_legend_title    = "Natural stand",
    show_legends          = TRUE,   # show color key in b) as well
    cluster_cols          = TRUE,
    show_annot_names      = FALSE   # no sidebar / topbar labels in the plot
    # ctx_vec_sig has only CpG + CHH — CHG is automatically absent from legend
  )

  # Use grid viewports to place panel a (70%) and panel b (30%) side by side
  draw_panel <- function() {
    grid.newpage()

    # Panel a — breeding
    pushViewport(viewport(x = 0, y = 0, width = frac_a, height = 1,
                          just = c("left", "bottom")))
    draw(ht_a_panel,
         heatmap_legend_side    = "right",
         annotation_legend_side = "right",
         merge_legend           = TRUE,
         newpage                = FALSE)
    grid.text("a",
              x    = unit(0.35, "cm"),
              y    = unit(1, "npc") - unit(0.35, "cm"),
              just = c("left", "top"),
              gp   = label_gp)
    popViewport()

    # Panel B — natural formal SVMPs
    pushViewport(viewport(x = frac_a, y = 0, width = 1 - frac_a, height = 1,
                          just = c("left", "bottom")))
    draw(ht_b_panel,
         heatmap_legend_side    = "right",
         annotation_legend_side = "right",
         merge_legend           = TRUE,
         newpage                = FALSE)
    grid.text("b",
              x    = unit(0.35, "cm"),
              y    = unit(1, "npc") - unit(0.35, "cm"),
              just = c("left", "top"),
              gp   = label_gp)
    popViewport()
  }

  panel_tiff <- file.path(fig_dir, "svmp_heatmaps_panel.tiff")
  tiff(panel_tiff, width = panel_w_cm, height = panel_h_cm,
       units = "cm", res = 600, compression = "lzw")
  draw_panel(); dev.off()
  pdf(sub("\\.tiff$", ".pdf", panel_tiff),
      width = panel_w_cm / 2.54, height = panel_h_cm / 2.54)
  draw_panel(); dev.off()
  cairo_ps(sub("\\.tiff$", ".eps", panel_tiff),
           width = panel_w_cm / 2.54, height = panel_h_cm / 2.54)
  draw_panel(); dev.off()
  png(sub("\\.tiff$", ".png", panel_tiff),
      width = panel_w_cm, height = panel_h_cm, units = "cm", res = 600)
  draw_panel(); dev.off()

  # ---- Values shown in the combined panel, one row per tree and locus
  heat_long <- function(X, na_mask, groups, cohort, panel) {
    data.frame(panel = panel, cohort = cohort,
               tree = rep(rownames(X), times = ncol(X)),
               group = rep(as.character(groups[rownames(X)]), times = ncol(X)),
               context = rep(sub("\\|.*$", "", colnames(X)), each = nrow(X)),
               locus = rep(sub("^[^|]*\\|", "", colnames(X)), each = nrow(X)),
               percent_methylation = as.vector(X),
               filled_with_locus_mean = as.vector(na_mask[rownames(X), colnames(X), drop = FALSE]))
  }
  write.table(rbind(heat_long(res_b$X_sel_bal, res_b$NA_all, res_b$group_vec, "BREEDING", "a"),
                    heat_long(X_sig, res_n$NA_all, res_n$group_vec, "NATURAL", "b")),
              file.path(fig_dir, "svmp_heatmaps_panel_values.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  cat("Panel b standalone (formal SVMPs): ", out_n_sig, "\n")
  cat("Combined panel:                        ", panel_tiff, "\n")

} else {
  cat("Fewer than 2 formal SVMPs in natural cohort — skipping 6-SVMP heatmap and panel.\n")
}

cat("\nDONE Step 8b.\n\n")
cat("Heatmaps saved in:\n  ", fig_dir, "\n\n", sep = "")
cat("Marker test tables + selected marker lists + logs saved in:\n  ", tab_dir, "\n\n", sep = "")
cat("Selected markers (BALANCED):\n  ", res_b$selected_bal_tsv, "\n  ", res_n$selected_bal_tsv, "\n", sep = "")
cat("Selected markers (OVERALL):\n  ", res_b$selected_all_tsv, "\n  ", res_n$selected_all_tsv, "\n", sep = "")
cat("Logs:\n  ", res_b$log_file, "\n  ", res_n$log_file, "\n", sep = "")
sessionInfo()
```

---

### `tgc.tms.functions.R`

```r
############################################################
# TreeGeneClimate (TGC) — TMS
# Shared functions for Steps 5b1-5b3 (sourced, not run directly)
############################################################

suppressPackageStartupMessages({
  library(methylKit)
  library(data.table)
})

# Read the SNP table written by Step 4a3 (chr, pos, ref, alt, qual, ac, an) and
# keep SNPs whose alternative allele(s) were called at least min_ac times
# (summed over alleles for multi-allelic sites). No call-rate filter: with about
# 7 ECS reads per tree many real SNPs are genotyped in only part of the trees,
# and such SNPs at a cytosine or its context bases still turn methylation
# differences into genotype differences.
read_snp_mask <- function(file, min_ac = 2L) {
  m <- fread(file, select = c("chr", "pos", "ac"), colClasses = list(character = "ac"))
  m[, ac_sum := vapply(strsplit(ac, ",", fixed = TRUE), function(x) sum(as.numeric(x)), numeric(1))]
  m <- m[ac_sum >= min_ac, .(chr, pos)]
  m <- unique(m)
  setkey(m, chr, pos)
  m
}

# Bases that define each cytosine context, counted from the cytosine on its own
# strand: CpG = C and the following G; CHG and CHH = C and the two following bases
MASK_WINDOW <- c(CpG = 1L, CHG = 2L, CHH = 2L)

# TRUE for sites with a SNP at the cytosine or at the bases defining its context
# on the same strand: [pos, pos + window] for "+" and [pos - window, pos] for "-".
# Sites with unknown strand ("*") are checked on both sides.
sites_hit_by_mask <- function(chr, pos, strand, mask, window = 2L) {
  n <- length(pos)
  hit <- logical(n)
  sgn <- ifelse(strand == "-", -1L, 1L)
  offsets <- 0:window
  for (o in offsets) {
    q <- data.table(chr = chr, pos = pos + sgn * o)
    hit <- hit | !is.na(mask[q, on = .(chr, pos), which = TRUE, mult = "first"])
    if (any(strand == "*")) {
      u <- which(strand == "*")
      q2 <- data.table(chr = chr[u], pos = pos[u] - o)
      hit[u] <- hit[u] | !is.na(mask[q2, on = .(chr, pos), which = TRUE, mult = "first"])
    }
  }
  hit
}

# Remove masked sites from a methylRaw object; returns list(object, n_removed)
mask_methylRaw <- function(mr, mask, context) {
  d <- getData(mr)
  hit <- sites_hit_by_mask(as.character(d$chr), as.integer(d$start), as.character(d$strand),
                           mask, window = MASK_WINDOW[[context]])
  list(object = mr[!hit, ], n_removed = sum(hit))
}

# Per-sample methylRaw files written by Step 5b1
raw_files <- function(raw_dir, cohort, context) {
  pattern <- if (cohort == "BREEDING") "^P00[1-3]_" else "^P00[4-8]_"
  f <- list.files(raw_dir, pattern = paste0(pattern, ".*_", context, "\\.rds$"), full.names = TRUE)
  f <- sort(f)
  names(f) <- sub(paste0("_", context, "\\.rds$"), "", basename(f))
  f
}

# Pooled methylated fraction per site across samples (sum numCs / sum coverage);
# Step 5b2 keeps sites whose pooled methylated fraction exceeds a threshold
# (methylated-fraction filter; output files carry the suffix _mef<threshold>)
pooled_methylated_fraction <- function(mb) {
  d <- getData(mb)
  cs <- rowSums(as.matrix(d[, grep("^numCs[0-9]+$", names(d))]), na.rm = TRUE)
  ts <- rowSums(as.matrix(d[, grep("^numTs[0-9]+$", names(d))]), na.rm = TRUE)
  ifelse(cs + ts > 0, cs / (cs + ts), NA_real_)
}
```

---

## JOINT — Integrated ECS + TMS Analysis

---

### `11ab.tgc.joint.correlation.analysis.R`

```r
#!/usr/bin/env Rscript
############################################################
# TreeGeneClimate (TGC) — JOINT (ECS + TMS)
# Step 11ab: Distance-based congruence between
#            genomic (ECS SNPs) and epigenomic (TMS methylation)
#
# PURPOSE
# - For each cohort (BREEDING, NATURAL) and each methylation context (CpG/CHG/CHH):
#   1) Build a genomic distance matrix from the called (non-imputed) ECS genotypes
#      (IBS distance = 1 - IBS; missing genotypes excluded pair by pair)
#   2) Build an epigenomic distance matrix from TMS methylBase (Euclidean distance on
#      % methylation over the sites observed in both trees, scaled to all sites)
#   3) Align samples robustly (trimws + tolower), save diagnostics for mismatches
#   4) Compare distance matrices:
#       - Mantel test
#       - Procrustes / protest (PCoA ordinations)
#       - RV coefficient
#   5) Save per-comparison tables + combined 3x2 Procrustes panel + aggregated summary
#
# INPUTS
# - ECS:
#   RESULTS/ECS/RANALYSIS/RDATA/breeding.called.snp.gds   (Step 8a)
#   RESULTS/ECS/RANALYSIS/RDATA/natural.called.snp.gds
# - TMS:
#   RESULTS/TMS/RANALYSIS/METHYLKIT_OBJECTS/
#     methylBase_<cohort>_<context>_cov10_p99.9_mpg*_mef0.05.rds
# - Metadata:
#   DATA/METADATA/breeding_sample2family.txt
#   DATA/METADATA/natural_sample2pop.txt
#
# OUTPUTS
# - PROJECT_ROOT/RESULTS/JOINT/
#     FIGURES/11ab/
#     TABLES/11ab/
#     FIGURES/11ab/procrustes_panel_values.tsv (values shown in the Procrustes panel)
#
# NOTES
# - PERMUTATIONS fixed at 9999.
# - Neither distance uses imputed values. Filling missing genotypes (BEAGLE) or
#   missing methylation values (site mean) pulls trees with more missing data
#   towards the centre of both distance matrices; because ECS and TMS libraries
#   were made from the same DNA, missingness is correlated between the two data
#   sets and would by itself create genomic-epigenomic concordance.
# - Run in a fresh R session if GDS files were previously opened.
############################################################

suppressPackageStartupMessages({
  library(SNPRelate)
  library(gdsfmt)
  library(methylKit)
  library(dplyr)
  library(tibble)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(vegan)   # mantel, procrustes, protest
  library(ade4)    # RV coefficient
  library(patchwork)
  library(grid)
})

options(stringsAsFactors = FALSE)
set.seed(1)

############################################################
# 1) FIXED PATHS
############################################################
# === USER CONFIGURATION ===
# Set TGC_PROJECT_ROOT as an environment variable, or edit the fallback path below
PROJECT_ROOT <- Sys.getenv("TGC_PROJECT_ROOT", unset = "/path/to/your/project")
# ===========================

# ECS (GDS produced by 8a)
ECS_RDATA_DIR <- file.path(PROJECT_ROOT, "RESULTS/ECS/RANALYSIS/RDATA")
GDS_BREED <- file.path(ECS_RDATA_DIR, "breeding.called.snp.gds")
GDS_NATUR <- file.path(ECS_RDATA_DIR, "natural.called.snp.gds")

# TMS methylKit objects
TMS_RDS_DIR <- file.path(PROJECT_ROOT, "RESULTS/TMS/RANALYSIS/METHYLKIT_OBJECTS")

# Metadata maps
META_DIR <- file.path(PROJECT_ROOT, "DATA/METADATA")
MAP_BREED <- file.path(META_DIR, "breeding_sample2family.txt")
MAP_NATUR <- file.path(META_DIR, "natural_sample2pop.txt")

# JOINT outputs
JOINT_ROOT <- file.path(PROJECT_ROOT, "RESULTS/JOINT")
OUT_FIG <- file.path(JOINT_ROOT, "FIGURES", "11ab")
OUT_TAB <- file.path(JOINT_ROOT, "TABLES",  "11ab")
OUT_DIAG <- file.path(OUT_TAB, "DIAGNOSE")

dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIAG, recursive = TRUE, showWarnings = FALSE)

############################################################
# 2) PARAMETERS
############################################################
PERMUTATIONS <- 9999
# SNPs for the IBS distance: within-cohort MAF >= 0.05 and called in >= 99% of
# the trees (as for the GRM, Step 8a), so that differences in call rate between
# trees do not shape the distances
ECS_MAF_MIN  <- 0.05
ECS_MISS_MAX <- 0.01
# Trees observed at fewer than this fraction of the retained methylation sites
# (failed or very low-complexity TMS libraries) are left out of the comparison
TMS_MIN_OBS  <- 0.20

# Context-specific NA tolerance: CHH sites have higher missing rates because
# non-CpG methylation is sparse and more variable across individuals
max_na_frac_by_ctx <- function(ctx) {
  switch(ctx,
         "CpG" = 0.20,
         "CHG" = 0.30,
         "CHH" = 0.50,
         0.30)
}

############################################################
# 3) HELPERS
############################################################
ensure_file <- function(p) if (!file.exists(p)) stop("Missing file: ", p, call. = FALSE)

# Normalise sample IDs to a canonical form for cross-dataset matching
norm_id <- function(x) tolower(trimws(as.character(x)))

# Select the TMS RDS with the highest mpg threshold (most stringent coverage filter)
pick_tms_mef_rds <- function(cohort, ctx) {
  pat <- sprintf("^methylBase_%s_%s_cov10_p99\\.9_mpg[0-9]+_mef0\\.05\\.rds$",
                 tolower(cohort), tolower(ctx))
  files <- list.files(TMS_RDS_DIR, pattern = pat, full.names = TRUE)
  if (!length(files)) stop("No methylated-fraction-filtered TMS RDS found for ", cohort, " / ", ctx, " under ", TMS_RDS_DIR)
  mpg_num <- as.integer(sub(".*_mpg([0-9]+)_mef0\\.05\\.rds$", "\\1", basename(files)))
  files[which.max(mpg_num)]
}

# ---- ECS genomic distance: IBS -> distance (1-IBS) ----
ecs_distance <- function(gds) {

  g <- snpgdsOpen(gds)
  on.exit(try(snpgdsClose(g), silent = TRUE), add = TRUE)

  stat <- snpgdsSNPRateFreq(g, with.id = TRUE)

  maf  <- stat$MinorFreq
  miss <- stat$MissingRate

  # Apply MAF and missingness QC before computing IBS
  keep <- is.finite(maf) & is.finite(miss) & maf >= ECS_MAF_MIN & miss <= ECS_MISS_MAX

  if (sum(keep) == 0) stop("No ECS SNPs passed MAF/missingness filters.")

  ibs <- snpgdsIBS(
    g,
    snp.id = stat$snp.id[keep],
    maf = NaN, missing.rate = NaN,  # filters applied above, no further defaults
    num.thread = 8,
    autosome.only = FALSE  # spruce has no canonical autosomes in the reference
  )

  M <- ibs$ibs
  rownames(M) <- colnames(M) <- ibs$sample.id

  # Convert IBS similarity to a distance matrix
  D <- as.dist(1 - M)

  list(
    dist = D,
    ids = ibs$sample.id,
    n_snps = sum(keep)
  )
}

# ---- TMS epigenomic distance: Euclidean on %methylation (sites filtered, no imputation) ----
methylbase_to_percent_matrix <- function(mb) {
  d <- getData(mb)
  numCs_cols <- grep("^numCs[0-9]+$", colnames(d), value = TRUE)
  numTs_cols <- grep("^numTs[0-9]+$", colnames(d), value = TRUE)
  if (length(numCs_cols) == 0 || length(numTs_cols) == 0 || length(numCs_cols) != length(numTs_cols)) {
    stop("Could not find matching numCs/numTs columns in methylBase.")
  }
  sample_ids <- mb@sample.ids
  if (is.null(sample_ids) || length(sample_ids) != length(numCs_cols)) {
    sample_ids <- paste0("S", seq_along(numCs_cols))
  }

  # Compute per-site percentage methylation for each sample
  perc_mat <- vapply(seq_along(numCs_cols), function(i) {
    numCs <- d[[numCs_cols[i]]]
    numTs <- d[[numTs_cols[i]]]
    p <- 100 * (numCs / (numCs + numTs))
    p[is.nan(p)] <- NA_real_  # sites with zero coverage yield NaN
    p
  }, numeric(nrow(d)))

  colnames(perc_mat) <- sample_ids
  perc_mat
}

# Keep sites missing in at most max_na_frac of the trees and variable across trees;
# remaining missing values are left as NA (no imputation)
filter_tms_sites <- function(X_samples_x_sites, max_na_frac) {
  na_frac <- colMeans(is.na(X_samples_x_sites))
  X2 <- X_samples_x_sites[, na_frac <= max_na_frac, drop = FALSE]
  sds <- apply(X2, 2, sd, na.rm = TRUE)
  X3 <- X2[, is.finite(sds) & sds > 0, drop = FALSE]
  if (ncol(X3) < 10) stop("Too few variable TMS sites after NA/SD filter (n=", ncol(X3), ").")
  X3
}

compute_tms_euclid_dist <- function(tms_mef_rds, ctx) {
  ensure_file(tms_mef_rds)
  mb <- readRDS(tms_mef_rds)

  perc_sites_x_samples <- methylbase_to_percent_matrix(mb)
  # Transpose so rows = samples, columns = sites
  X <- t(perc_sites_x_samples)
  rownames(X) <- colnames(perc_sites_x_samples)

  max_na <- max_na_frac_by_ctx(ctx)
  X2 <- filter_tms_sites(X, max_na_frac = max_na)
  obs <- rowMeans(!is.na(X2))
  if (any(obs < TMS_MIN_OBS))
    message("TMS: trees observed at < ", TMS_MIN_OBS * 100, "% of sites left out: ",
            paste(rownames(X2)[obs < TMS_MIN_OBS], collapse = ", "))
  X2 <- X2[obs >= TMS_MIN_OBS, , drop = FALSE]
  # Euclidean distance across all filtered sites; dist() excludes missing values
  # pair by pair and scales the sum of squares up to the full number of sites
  D <- dist(X2, method = "euclidean")
  if (anyNA(D)) stop("Some tree pairs share no TMS site")
  list(dist = D, ids = rownames(X2), n_sites = ncol(X2), excluded = sum(obs < TMS_MIN_OBS))
}

# ---- Align by normalized IDs and save diagnostics ----
align_dist_by_ids <- function(D1, ids1, D2, ids2, out_diag_prefix) {
  n1 <- norm_id(ids1)
  n2 <- norm_id(ids2)

  # Detect duplicate IDs within each dataset before intersection
  dup1 <- duplicated(n1)
  dup2 <- duplicated(n2)
  if (any(dup1) || any(dup2)) {
    write_tsv(
      tibble(
        side = c(rep("ECS", sum(dup1)), rep("TMS", sum(dup2))),
        id   = c(ids1[dup1], ids2[dup2]),
        norm = c(n1[dup1], n2[dup2])
      ),
      file.path(OUT_DIAG, paste0(out_diag_prefix, "_duplicate_norm_ids.tsv"))
    )
    stop("Duplicate normalized IDs detected. See diagnostics: ", out_diag_prefix, "_duplicate_norm_ids.tsv")
  }

  map1 <- setNames(ids1, n1)
  map2 <- setNames(ids2, n2)
  common_norm <- intersect(names(map1), names(map2))

  # Log which samples are present in only one dataset for QC review
  only_1 <- setdiff(names(map1), common_norm)
  only_2 <- setdiff(names(map2), common_norm)

  write_tsv(tibble(side = "ECS_only", norm_id = only_1, example_id = map1[only_1]),
            file.path(OUT_DIAG, paste0(out_diag_prefix, "_missing_ECSonly.tsv")))
  write_tsv(tibble(side = "TMS_only", norm_id = only_2, example_id = map2[only_2]),
            file.path(OUT_DIAG, paste0(out_diag_prefix, "_missing_TMSonly.tsv")))

  if (length(common_norm) < 10) stop("Too few overlapping samples after ID normalization: ", length(common_norm))

  ids1_keep <- unname(map1[common_norm])
  ids2_keep <- unname(map2[common_norm])

  # Subset and re-order both distance matrices to the common sample set
  M1 <- as.matrix(D1); rownames(M1) <- colnames(M1) <- ids1
  M2 <- as.matrix(D2); rownames(M2) <- colnames(M2) <- ids2

  M1s <- M1[ids1_keep, ids1_keep, drop = FALSE]
  M2s <- M2[ids2_keep, ids2_keep, drop = FALSE]

  # Use normalised IDs as the canonical row/col names for both matrices
  ord <- common_norm
  rownames(M1s) <- colnames(M1s) <- ord
  rownames(M2s) <- colnames(M2s) <- ord

  list(D1 = as.dist(M1s), D2 = as.dist(M2s), common_norm = ord, n = length(ord))
}

# ---- Ordination (2 axes) ----
# PCoA (= classical MDS) on a distance matrix; add=TRUE corrects negative eigenvalues
pcoa2 <- function(D) {
  m <- cmdscale(D, k = 2, eig = TRUE, add = TRUE)
  X <- as.data.frame(m$points)
  colnames(X) <- c("Axis1", "Axis2")
  X
}

# ---- Procrustes plot object for panel ----
make_procrustes_plot <- function(X, Y, title_label = "a") {
  # Symmetric Procrustes: both configurations are scaled and rotated optimally
  pro <- vegan::procrustes(X, Y, symmetric = TRUE)
  Yrot <- pro$Yrot

  df <- tibble(
    id = rownames(X),
    X1 = X[,1], X2 = X[,2],
    Y1 = Yrot[,1], Y2 = Yrot[,2]
  )

  # Segments connect each individual's genomic (SNP) position to its epigenomic (SMP) position
  p <- ggplot(df) +
    geom_segment(aes(x = Y1, y = Y2, xend = X1, yend = X2), alpha = 0.45, color = "grey50") +
    geom_point(aes(x = X1, y = X2, color = "SNPs"), size = 1.7) +
    geom_point(aes(x = Y1, y = Y2, color = "SMPs"), size = 1.7, alpha = 0.9) +
    scale_color_manual(values = c("SNPs" = "#4e79a7", "SMPs" = "#f28e2b"), name = NULL) + guides(color = guide_legend(override.aes = list(size = 6))) +
    labs(x = "PC1", y = "PC2", tag = title_label) +   # panel letter outside the plotting area
    theme_minimal(base_size = 18) +
    theme(
      panel.grid.minor = element_blank(),
      plot.tag = element_text(face = "bold", size = 20),
      legend.position = "bottom",
      legend.direction = "horizontal",
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 12),
      legend.text = element_text(size = 16),
      legend.key.size = unit(4, "lines")
    )

  list(plot = p, pro = pro, df = df)
}

############################################################
# 4) SANITY CHECKS
############################################################
ensure_file(GDS_BREED); ensure_file(GDS_NATUR)
ensure_file(MAP_BREED); ensure_file(MAP_NATUR)

############################################################
# 5) RUN: 2 cohorts × 3 contexts
############################################################
cohorts <- list(
  BREEDING = list(gds = GDS_BREED, map = MAP_BREED),
  NATURAL  = list(gds = GDS_NATUR, map = MAP_NATUR)
)
contexts <- c("CpG", "CHG", "CHH")

plot_list <- list()
summary_rows <- list()
coord_rows <- list()

# Panel labels a–f, top-left to bottom-right
panel_labels <- c("a", "b", "c", "d", "e", "f")
panel_i <- 1

for (coh in names(cohorts)) {
  message("\n==============================")
  message("COHORT: ", coh)
  message("==============================")

  # Compute IBS distance once per cohort (shared across all three contexts)
  ecs <- ecs_distance(cohorts[[coh]]$gds)
  message("ECS: IBS distance computed (SNPs passing QC = ", ecs$n_snps, ")")

  for (ctx in contexts) {
    message("\n--- ", coh, " / ", ctx, " ---")

    tms_rds <- pick_tms_mef_rds(coh, ctx)
    tem <- compute_tms_euclid_dist(tms_rds, ctx)
    message("TMS: Euclidean distance computed (sites used = ", tem$n_sites, ")")

    prefix <- paste0("11ab_", tolower(coh), "_", tolower(ctx))
    al <- align_dist_by_ids(ecs$dist, ecs$ids, tem$dist, tem$ids, out_diag_prefix = prefix)
    message("Aligned samples: n = ", al$n)

    # Mantel test: correlation between the two distance matrices (permutation-based)
    man <- vegan::mantel(al$D1, al$D2, method = "pearson", permutations = PERMUTATIONS)

    # PCoA of each distance matrix; Procrustes / protest then compares the ordinations
    X <- pcoa2(al$D1)
    Y <- pcoa2(al$D2)
    rownames(X) <- rownames(Y) <- al$common_norm

    pro_test <- vegan::protest(X, Y, permutations = PERMUTATIONS)
    t0 <- unname(pro_test$t0)   # Procrustes statistic (lower = better alignment)
    p_pro <- unname(pro_test$signif)

    # RV coefficient: multivariate analogue of the squared correlation
    rv_res <- ade4::RV.rtest(as.data.frame(X), as.data.frame(Y), nrepet = PERMUTATIONS)
    rv <- unname(rv_res$obs)
    p_rv <- unname(rv_res$pvalue)

    row <- tibble(
      cohort = coh,
      context = ctx,
      n_samples = al$n,
      ecs_snps_qc = ecs$n_snps,
      tms_sites_used = tem$n_sites,
      mantel_r = unname(man$statistic),
      mantel_p = unname(man$signif),
      procrustes_t0 = t0,
      procrustes_p = p_pro,
      rv = rv,
      rv_p = p_rv,
      tms_rds = tms_rds
    )

    out_sum <- file.path(OUT_TAB, paste0(prefix, "_summary.tsv"))
    write_tsv(row, out_sum)

    pro_obj <- make_procrustes_plot(X, Y, title_label = panel_labels[panel_i])
    plot_list[[panel_i]] <- pro_obj$plot
    coord_rows[[prefix]] <- data.frame(panel = panel_labels[panel_i], cohort = coh, context = ctx,
                                       tree = pro_obj$df$id,
                                       snp_PC1 = pro_obj$df$X1, snp_PC2 = pro_obj$df$X2,
                                       smp_PC1_rotated = pro_obj$df$Y1, smp_PC2_rotated = pro_obj$df$Y2)

    # Per-sample Procrustes residuals indicate how well each individual is matched
    Yrot <- pro_obj$pro$Yrot
    residuals <- sqrt(rowSums((as.matrix(X) - as.matrix(Yrot))^2))
    out_res <- file.path(OUT_TAB, paste0(prefix, "_procrustes_residuals.tsv"))
    write_tsv(tibble(sample_norm = al$common_norm, residual = residuals), out_res)

    summary_rows[[prefix]] <- row
    message("Saved: ", basename(out_sum), " | ", basename(out_res))

    panel_i <- panel_i + 1
  }
}

############################################################
# 6) COMBINED 3x2 PROCRUSTES PANEL WITH COMMON LEGEND
############################################################
# Constrain y-axis for NATURAL panels (D–F) which have tighter ordination spread
plot_list[[4]] <- plot_list[[4]] + coord_cartesian(ylim = c(-0.05, 0.05))
plot_list[[5]] <- plot_list[[5]] + coord_cartesian(ylim = c(-0.05, 0.05))
plot_list[[6]] <- plot_list[[6]] + coord_cartesian(ylim = c(-0.05, 0.05))

# Row 1: BREEDING (CpG / CHG / CHH); Row 2: NATURAL (CpG / CHG / CHH)
combined_panel <-
  ((plot_list[[1]] | plot_list[[2]] | plot_list[[3]]) /
     (plot_list[[4]] | plot_list[[5]] | plot_list[[6]])) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal"
  )

out_panel <- file.path(OUT_FIG, "procrustes_panel.tiff")

# Save as TIFF, PDF, EPS and PNG
tiff(out_panel, width = 34, height = 24, units = "cm", res = 600, compression = "lzw")
print(combined_panel)
dev.off()
ggsave(sub("\\.tiff$", ".pdf", out_panel), plot = combined_panel, width = 34, height = 24, units = "cm")
ggsave(sub("\\.tiff$", ".eps", out_panel), plot = combined_panel, width = 34, height = 24, units = "cm", device = cairo_ps)
ggsave(sub("\\.tiff$", ".png", out_panel), plot = combined_panel, width = 34, height = 24, units = "cm", dpi = 150, device = "png")

# Values shown in the panel: ordination coordinates of every tree
write_tsv(bind_rows(coord_rows), file.path(OUT_FIG, "procrustes_panel_values.tsv"))

############################################################
# 7) AGGREGATED SUMMARY + BH CORRECTION
############################################################
# Apply BH correction across all 6 panels jointly for each test type
all_df <- bind_rows(summary_rows) %>%
  mutate(
    mantel_p_adj = p.adjust(mantel_p, method = "BH"),
    procrustes_p_adj = p.adjust(procrustes_p, method = "BH"),
    rv_p_adj = p.adjust(rv_p, method = "BH")
  ) %>%
  arrange(cohort, context)

out_all <- file.path(OUT_TAB, "11ab_ecs_tms_distance_summary_all.tsv")
write_tsv(all_df, out_all)

message("\nDONE 11ab.")
message("Combined figure: ", out_panel)
message("Tables: ", OUT_TAB)
message("Master summary: ", out_all)
sessionInfo()
```

---

### `12ab0.tgc.joint.meqtl.input.prep.R`

```r
#!/usr/bin/env Rscript
############################################################
# TreeGeneClimate (TGC) — JOINT ECS + TMS
# Step 12ab0: Inputs for cis-meQTL mapping (both cohorts, three contexts)
#
# GENOTYPES (per cohort, shared by Steps 12ab1 and 12ab2)
#   - BEAGLE-imputed genotypes (Step 8a GDS; best-guess calls 0/1/2)
#   - SNPs with within-cohort MAF >= 0.05 and missing rate <= 0.20
#   - Residual missing calls (variants on single-variant contigs, which
#     BEAGLE does not impute) set to the rounded expected genotype from the
#     allele frequency of the tree's family (breeding) or stand (natural).
#     Both mapping tools read the same integer genotype matrix.
# COVARIATES
#   - PC1-PC10 of the PCA on LD-pruned SNPs (Step 5a/8a)
# RELATEDNESS (GENESIS)
#   - GCTA-estimator GRMs from called (non-imputed) genotypes of SNPs called
#     in >= 99% of the trees, as in Step 8a
#   - BREEDING: leave-one-chromosome-out (LOCO) GRMs, all SNPs except those on
#     the tested chromosome/scaffold (exact, by subtracting its components)
#   - NATURAL: genome-wide GRM (Step 8a)
# METHYLATION (per cohort x context)
#   - Step 5b3 datasets; sites with zero variance removed
#   - M = log2[(beta x 100 + 0.5) / (100 - beta x 100 + 0.5)]
#
# INPUTS
#   RESULTS/ECS/RANALYSIS/RDATA/<cohort>.imputed.snp.gds, <cohort>.called.snp.gds
#   RESULTS/ECS/RANALYSIS/RDATA/<cohort>_grm_gcta.rds
#   RESULTS/ECS/RANALYSIS/RDATA/<cohort>_pca_pruned.rds
#   RESULTS/ECS/RANALYSIS/RDATA/<cohort>_sample_annotation.rds
#   RESULTS/TMS/RANALYSIS/METHYLKIT_OBJECTS_MEQTL/methylBase_<cohort>_<ctx>_meqtl.rds
#
# OUTPUTS (RESULTS/JOINT/MEQTL/INPUTS/<COHORT>/)
#   genotypes.rds                 integer matrix, SNPs x samples
#   snp_variant_annot.tsv         snp_id chr pos ref alt maf missing_rate
#   loco_grm/loco_grm_excl_<chr>.rds   (BREEDING)
#   <CTX>/methylation_mvalues_matrix.rds   samples x sites
#   <CTX>/methylation_site_annot.tsv       methyl_loc chr start end strand pos
#   <CTX>/shared_sample_ids.tsv, pcs_shared.tsv, grm_shared.rds
#   ../step12ab0_summary.tsv
#
# R packages: SNPRelate, gdsfmt, methylKit, data.table
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(SNPRelate)
  library(gdsfmt)
  library(methylKit)
})
options(stringsAsFactors = FALSE)

# === USER CONFIGURATION ===
PROJECT_ROOT <- Sys.getenv("TGC_PROJECT_ROOT", unset = "/path/to/your/project")
N_PCS    <- 10L
MIN_MAF  <- 0.05
MAX_MISS <- 0.20
EPS_M    <- 0.5
# ===========================
source(file.path(PROJECT_ROOT, "SCRIPTS/ECS/tgc.ecs.functions.R"))

ECS_RDATA <- file.path(PROJECT_ROOT, "RESULTS/ECS/RANALYSIS/RDATA")
METH_DIR  <- file.path(PROJECT_ROOT, "RESULTS/TMS/RANALYSIS/METHYLKIT_OBJECTS_MEQTL")
OUT_ROOT  <- file.path(PROJECT_ROOT, "RESULTS/JOINT/MEQTL/INPUTS")
dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)
COHORTS  <- c("BREEDING", "NATURAL")
CONTEXTS <- c("CpG", "CHG", "CHH")
GROUP_COL <- c(BREEDING = "Family", NATURAL = "Population")

msg <- function(...) cat("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n", sep = "")

# Missing genotypes (SNPs x samples) -> 2 x group allele frequency
fill_by_group <- function(geno, grp, round_values = FALSE) {
  p_all <- rowMeans(geno, na.rm = TRUE) / 2
  for (g in unique(grp)) {
    idx  <- which(grp == g)
    sub  <- geno[, idx, drop = FALSE]
    miss <- is.na(sub)
    if (!any(miss)) next
    pg <- rowMeans(sub, na.rm = TRUE) / 2
    pg[!is.finite(pg)] <- p_all[!is.finite(pg)]
    fill <- 2 * pg
    if (round_values) fill <- round(fill)
    sub[miss] <- fill[row(sub)[miss]]
    geno[, idx] <- sub
  }
  geno
}

# Methylation % (sites x samples) from a methylBase object
perc_methylation <- function(mb) {
  d  <- getData(mb)
  cs <- as.matrix(d[, grep("^numCs[0-9]+$", names(d))])
  ts <- as.matrix(d[, grep("^numTs[0-9]+$", names(d))])
  p  <- 100 * cs / (cs + ts)
  p[!is.finite(p)] <- NA_real_
  colnames(p) <- mb@sample.ids
  p
}

summary_rows <- list()

for (cohort in COHORTS) {
  coh <- tolower(cohort)
  out_c <- file.path(OUT_ROOT, cohort)
  dir.create(out_c, recursive = TRUE, showWarnings = FALSE)
  msg("=== ", cohort, " ===")

  # ---------------- genotypes ----------------
  annot <- as.data.table(readRDS(file.path(ECS_RDATA, sprintf("%s_sample_annotation.rds", coh))))
  gf <- snpgdsOpen(file.path(ECS_RDATA, sprintf("%s.imputed.snp.gds", coh)))
  samples <- read.gdsn(index.gdsn(gf, "sample.id"))
  va <- data.table(snp_id = read.gdsn(index.gdsn(gf, "snp.id")),
                   chr    = as.character(read.gdsn(index.gdsn(gf, "snp.chromosome"))),
                   pos    = as.integer(read.gdsn(index.gdsn(gf, "snp.position"))),
                   allele = read.gdsn(index.gdsn(gf, "snp.allele")))
  geno <- snpgdsGetGeno(gf, snpfirstdim = TRUE, with.id = FALSE, verbose = FALSE)
  snpgdsClose(gf)
  mode(geno) <- "numeric"
  geno[geno < 0 | geno > 2] <- NA
  grp <- annot[[GROUP_COL[[cohort]]]][match(samples, annot$IID)]
  grp[is.na(grp)] <- "Unknown"

  va[, missing_rate := rowMeans(is.na(geno))]
  af <- rowMeans(geno, na.rm = TRUE) / 2
  va[, maf := pmin(af, 1 - af)]
  va[, c("ref", "alt") := tstrsplit(allele, "/", fixed = TRUE)]
  keep <- is.finite(va$maf) & va$maf >= MIN_MAF & va$missing_rate <= MAX_MISS
  msg(nrow(va), " imputed SNPs; kept ", sum(keep), " (MAF >= ", MIN_MAF, ", missing <= ", MAX_MISS, ")")

  # LOCO GRMs (breeding): GCTA estimator on called genotypes, as the full GRM of Step 8a
  if (cohort == "BREEDING") {
    msg("LOCO GRMs ...")
    gc_ <- snpgdsOpen(file.path(ECS_RDATA, sprintf("%s.called.snp.gds", coh)))
    stopifnot(identical(read.gdsn(index.gdsn(gc_, "sample.id")), samples))
    chr_c <- as.character(read.gdsn(index.gdsn(gc_, "snp.chromosome")))
    Xc <- snpgdsGetGeno(gc_, snpfirstdim = TRUE, with.id = FALSE, verbose = FALSE)
    snpgdsClose(gc_)
    mode(Xc) <- "numeric"
    Xc[Xc < 0 | Xc > 2] <- NA
    keep_g <- grm_snps(Xc)                       # SNPs called in >= 99% of trees, as Step 8a
    Xc <- Xc[keep_g, , drop = FALSE]; chr_c <- chr_c[keep_g]
    total <- gcta_components(Xc)
    loco_dir <- file.path(out_c, "loco_grm")
    dir.create(loco_dir, showWarnings = FALSE)
    for (ch in unique(chr_c)) {
      part <- gcta_components(Xc[chr_c == ch, , drop = FALSE])
      saveRDS(gcta_grm(gcta_combine(total, part), samples),
              file.path(loco_dir, sprintf("loco_grm_excl_%s.rds", ch)))
    }
    msg("LOCO GRMs written for ", length(unique(chr_c)), " chromosomes/scaffolds")
    rm(Xc, total, part); gc()
  }

  geno <- geno[keep, , drop = FALSE]
  va   <- va[keep]
  geno <- fill_by_group(geno, grp, round_values = TRUE)
  storage.mode(geno) <- "integer"
  dimnames(geno) <- list(as.character(va$snp_id), samples)
  saveRDS(geno, file.path(out_c, "genotypes.rds"))
  fwrite(va[, .(snp_id, chr, pos, ref, alt, maf = round(maf, 5), missing_rate = round(missing_rate, 5))],
         file.path(out_c, "snp_variant_annot.tsv"), sep = "\t")
  rm(geno); gc()

  # ---------------- covariates and GRM ----------------
  pca <- readRDS(file.path(ECS_RDATA, sprintf("%s_pca_pruned.rds", coh)))
  pcs <- as.data.table(pca$scores)[, c("IID", paste0("PC", seq_len(N_PCS))), with = FALSE]
  setnames(pcs, "IID", "sample_id")
  G_full <- readRDS(file.path(ECS_RDATA, sprintf("%s_grm_gcta.rds", coh)))

  # ---------------- methylation panels ----------------
  for (ctx in CONTEXTS) {
    out_p <- file.path(out_c, ctx)
    dir.create(out_p, recursive = TRUE, showWarnings = FALSE)
    mb <- readRDS(file.path(METH_DIR, sprintf("methylBase_%s_%s_meqtl.rds", coh, tolower(ctx))))
    d  <- getData(mb)
    pm <- perc_methylation(mb)
    n_united <- nrow(pm)
    keep_var <- rowSums(!is.na(pm)) > 1 & apply(pm, 1, sd, na.rm = TRUE) > 0
    pm <- pm[keep_var, , drop = FALSE]
    site_annot <- data.table(chr = as.character(d$chr), start = as.integer(d$start),
                             end = as.integer(d$end), strand = as.character(d$strand))[keep_var]
    site_annot[, pos := start]
    site_annot[, methyl_loc := paste0(chr, ":", start, "-", end)]
    M <- t(log2((pm + EPS_M) / (100 - pm + EPS_M)))   # samples x sites
    colnames(M) <- site_annot$methyl_loc

    shared <- sort(Reduce(intersect, list(samples, pcs$sample_id, rownames(G_full), rownames(M))))
    M <- M[shared, , drop = FALSE]
    saveRDS(M, file.path(out_p, "methylation_mvalues_matrix.rds"))
    fwrite(site_annot[, .(methyl_loc, chr, start, end, strand, pos)],
           file.path(out_p, "methylation_site_annot.tsv"), sep = "\t")
    fwrite(data.table(sample_id = shared), file.path(out_p, "shared_sample_ids.tsv"), sep = "\t")
    fwrite(pcs[match(shared, sample_id)], file.path(out_p, "pcs_shared.tsv"), sep = "\t")
    saveRDS(G_full[shared, shared], file.path(out_p, "grm_shared.rds"))

    summary_rows[[paste(cohort, ctx)]] <- data.table(
      cohort = cohort, context = ctx, n_samples = length(shared),
      n_sites_united = n_united, n_sites_tested = ncol(M), n_snps = nrow(va))
    msg(cohort, "/", ctx, ": ", length(shared), " samples, ", ncol(M), " sites (", n_united, " united)")
    rm(mb, d, pm, M); gc()
  }
}

summ <- rbindlist(summary_rows)
fwrite(summ, file.path(OUT_ROOT, "step12ab0_summary.tsv"), sep = "\t")
print(summ)
msg("Step 12ab0 done")
sessionInfo()
```

---

### `12ab1.tgc.joint.matrixeqtl.mapping.R`

```r
#!/usr/bin/env Rscript
############################################################
# TreeGeneClimate (TGC) — JOINT ECS + TMS
# Step 12ab1: cis-meQTL mapping with MatrixEQTL (one cohort x context panel)
#
# MODEL   M-value ~ SNP + PC1..PC10   (ordinary least squares)
# CIS     SNP within 100 kb of the methylation site
# OUTPUT  all cis pairs; Benjamini-Hochberg FDR over all cis pairs of the panel
#
# USAGE
#   TGC_COHORT=BREEDING TGC_CONTEXT=CpG Rscript 12ab1.tgc.joint.matrixeqtl.mapping.R
#
# INPUTS (Step 12ab0; RESULTS/JOINT/MEQTL/INPUTS/<COHORT>/)
#   genotypes.rds, snp_variant_annot.tsv
#   <CTX>/methylation_mvalues_matrix.rds, methylation_site_annot.tsv,
#   <CTX>/shared_sample_ids.tsv, pcs_shared.tsv
#
# OUTPUTS (RESULTS/JOINT/MEQTL/MATRIXEQTL/<COHORT>/<CTX>/)
#   cis_meqtl_all_results.rds / .tsv.gz
#     snp site snp_chr snp_pos site_pos distance beta statistic pvalue p_FDR
#   cis_meqtl_significant.tsv   (FDR < 0.05)
#   cis_meqtl_top_per_site.tsv
#   cis_meqtl_panel_summary.tsv (includes genomic inflation lambda)
#   qq_<cohort>_<ctx>.{tiff,pdf}
#
# R packages: MatrixEQTL, data.table
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(MatrixEQTL)
})
options(stringsAsFactors = FALSE)

# === USER CONFIGURATION ===
PROJECT_ROOT <- Sys.getenv("TGC_PROJECT_ROOT", unset = "/path/to/your/project")
CIS_WINDOW <- 1e5
FDR_THRESH <- 0.05
# ===========================

COHORT  <- Sys.getenv("TGC_COHORT");  CONTEXT <- Sys.getenv("TGC_CONTEXT")
stopifnot(COHORT %in% c("BREEDING", "NATURAL"), CONTEXT %in% c("CpG", "CHG", "CHH"))
IN_C  <- file.path(PROJECT_ROOT, "RESULTS/JOINT/MEQTL/INPUTS", COHORT)
IN_P  <- file.path(IN_C, CONTEXT)
OUT_P <- file.path(PROJECT_ROOT, "RESULTS/JOINT/MEQTL/MATRIXEQTL", COHORT, CONTEXT)
dir.create(OUT_P, recursive = TRUE, showWarnings = FALSE)
tag <- paste0(tolower(COHORT), "_", tolower(CONTEXT))
msg <- function(...) cat("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n", sep = "")

lambda_gc <- function(p) {
  p <- p[is.finite(p) & p > 0 & p <= 1]
  median(qchisq(p, df = 1, lower.tail = FALSE)) / qchisq(0.5, df = 1)
}

# ---- inputs ----
shared <- fread(file.path(IN_P, "shared_sample_ids.tsv"))$sample_id
M      <- readRDS(file.path(IN_P, "methylation_mvalues_matrix.rds"))[shared, , drop = FALSE]
sites  <- fread(file.path(IN_P, "methylation_site_annot.tsv"))
pcs    <- fread(file.path(IN_P, "pcs_shared.tsv"))
va     <- fread(file.path(IN_C, "snp_variant_annot.tsv"))
geno   <- readRDS(file.path(IN_C, "genotypes.rds"))[, shared, drop = FALSE]
stopifnot(identical(colnames(M), sites$methyl_loc), identical(pcs$sample_id, shared),
          identical(rownames(geno), as.character(va$snp_id)))
msg(COHORT, "/", CONTEXT, ": ", length(shared), " samples, ", ncol(M), " sites, ", nrow(geno), " SNPs")

snps <- SlicedData$new(); snps$CreateFromMatrix(geno * 1.0)
gene <- SlicedData$new(); gene$CreateFromMatrix(t(M))
cvrt <- SlicedData$new()
cvrt$CreateFromMatrix(t(as.matrix(pcs[, grep("^PC[0-9]+$", names(pcs)), with = FALSE])))
rm(geno); gc()

snpspos <- data.frame(snp = as.character(va$snp_id), chr = va$chr, pos = va$pos)
genepos <- data.frame(gene = sites$methyl_loc, chr = sites$chr, start = sites$start, end = sites$end)

# ---- mapping ----
me <- Matrix_eQTL_main(
  snps = snps, gene = gene, cvrt = cvrt,
  output_file_name = "", pvOutputThreshold = 0,
  useModel = modelLINEAR, errorCovariance = numeric(), verbose = FALSE,
  output_file_name.cis = "", pvOutputThreshold.cis = 1,
  snpspos = snpspos, genepos = genepos, cisDist = CIS_WINDOW,
  pvalue.hist = FALSE, min.pv.by.genesnp = FALSE, noFDRsaveMemory = FALSE)

res <- as.data.table(me$cis$eqtls)
setnames(res, c("snps", "gene"), c("snp", "site"))
res[, snp := as.character(snp)]; res[, site := as.character(site)]
res[, p_FDR := p.adjust(pvalue, method = "BH")]
res[, FDR := NULL]
res[, `:=`(snp_chr = va$chr[match(snp, as.character(va$snp_id))],
           snp_pos = va$pos[match(snp, as.character(va$snp_id))],
           site_pos = sites$pos[match(site, sites$methyl_loc)])]
res[, distance := abs(site_pos - snp_pos)]
stopifnot(max(res$distance) <= CIS_WINDOW)
setcolorder(res, c("snp", "site", "snp_chr", "snp_pos", "site_pos", "distance",
                   "beta", "statistic", "pvalue", "p_FDR"))

lam <- lambda_gc(res$pvalue)
sig <- res[p_FDR < FDR_THRESH][order(p_FDR)]
top <- res[res[, .I[which.min(pvalue)], by = site]$V1]
msg("cis pairs: ", nrow(res), " | FDR < ", FDR_THRESH, ": ", nrow(sig), " | lambda = ", round(lam, 4))

saveRDS(res, file.path(OUT_P, "cis_meqtl_all_results.rds"))
fwrite(res, file.path(OUT_P, "cis_meqtl_all_results.tsv.gz"), sep = "\t")
fwrite(sig, file.path(OUT_P, "cis_meqtl_significant.tsv"), sep = "\t")
fwrite(top, file.path(OUT_P, "cis_meqtl_top_per_site.tsv"), sep = "\t")
fwrite(data.table(tool = "MatrixEQTL", cohort = COHORT, context = CONTEXT,
                  n_samples = length(shared), n_sites = ncol(M), n_snps = nrow(va),
                  n_sites_with_cis = uniqueN(res$site), n_snps_with_cis = uniqueN(res$snp),
                  n_pairs = nrow(res), n_fdr05 = nrow(sig),
                  n_fdr1e10 = sum(res$p_FDR < 1e-10), lambda = round(lam, 4)),
       file.path(OUT_P, "cis_meqtl_panel_summary.tsv"), sep = "\t")

# ---- QQ plot (TIFF + PDF) ----
pv <- sort(res$pvalue[res$pvalue > 0])
qq <- function() {
  plot(-log10(ppoints(length(pv))), -log10(pv), pch = 20, cex = 0.4, col = "#2b8cbe",
       xlab = expression(Expected ~ -log[10](italic(p))), ylab = expression(Observed ~ -log[10](italic(p))))
  abline(0, 1, col = "red", lty = 2)
  mtext(sprintf("lambda = %.3f", lam), side = 3, adj = 0, line = 0.3)
}
tiff(file.path(OUT_P, sprintf("qq_%s.tiff", tag)), width = 2400, height = 2400, res = 300, compression = "lzw"); qq(); dev.off()
pdf(file.path(OUT_P, sprintf("qq_%s.pdf", tag)), width = 8, height = 8); qq(); dev.off()
msg("Step 12ab1 done")
```

---

### `12ab1.tgc.joint.matrixeqtl.mapping.sh`

```bash
#!/bin/bash
#SBATCH --job-name=JOINT.12ab1.MEQTL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --array=0-5
#SBATCH --output=LOGS/%x_%A_%a.out
#SBATCH --error=LOGS/%x_%A_%a.err
#-------------------------------------------------------------------------------
# TreeGeneClimate (TGC) — JOINT
# Step 12ab1: MatrixEQTL cis-meQTL mapping, one array task per cohort x context
#-------------------------------------------------------------------------------
set -euo pipefail

# === USER CONFIGURATION ===
export TGC_PROJECT_ROOT="${TGC_PROJECT_ROOT:-/path/to/your/project}"
module purge 2>/dev/null || true
module load gcc/14.2.0 r/4.5.2 2>/dev/null || true
export LC_ALL=C.UTF-8 LANG=C.UTF-8
# ===========================

cohorts=(BREEDING BREEDING BREEDING NATURAL NATURAL NATURAL)
contexts=(CpG CHG CHH CpG CHG CHH)
i="${SLURM_ARRAY_TASK_ID:-0}"
export TGC_COHORT="${cohorts[$i]}" TGC_CONTEXT="${contexts[$i]}"
Rscript --vanilla "${TGC_PROJECT_ROOT}/SCRIPTS/JOINT/12ab1.tgc.joint.matrixeqtl.mapping.R"
```

---

### `12ab2.tgc.joint.genesis.mapping.R`

```r
#!/usr/bin/env Rscript
############################################################
# TreeGeneClimate (TGC) — JOINT ECS + TMS
# Step 12ab2: cis-meQTL mapping with GENESIS (one cohort x context panel)
#
# MODEL   per methylation site, null linear mixed model
#           M-value ~ PC1..PC10 + g,  g ~ N(0, sigma2_g * GRM)
#         (AI-REML; BREEDING: LOCO GRM excluding the site's chromosome;
#          NATURAL: genome-wide GRM), then a score test for every SNP
#         within 100 kb of the site
# OUTPUT  all cis pairs; Benjamini-Hochberg FDR over all cis pairs of the panel
#
# The panel can be split into chunks of chromosomes/scaffolds that run in
# parallel, followed by a combine run:
#   TGC_COHORT=BREEDING TGC_CONTEXT=CHH TGC_NCHUNKS=8 TGC_CHUNK=3 Rscript 12ab2...R
#   TGC_COHORT=BREEDING TGC_CONTEXT=CHH TGC_NCHUNKS=8 TGC_CHUNK=combine Rscript 12ab2...R
# Without TGC_NCHUNKS the whole panel runs in one job.
#
# INPUTS (Step 12ab0; RESULTS/JOINT/MEQTL/INPUTS/<COHORT>/)
#   genotypes.rds, snp_variant_annot.tsv, loco_grm/ (BREEDING)
#   <CTX>/methylation_mvalues_matrix.rds, methylation_site_annot.tsv,
#   <CTX>/shared_sample_ids.tsv, pcs_shared.tsv, grm_shared.rds
#
# OUTPUTS (RESULTS/JOINT/MEQTL/GENESIS/<COHORT>/<CTX>/)
#   cis_meqtl_all_results.rds / .tsv.gz
#     snp site snp_chr snp_pos site_pos distance beta beta_se statistic
#     pvalue p_FDR freq MAC n_obs
#   cis_meqtl_significant.tsv, cis_meqtl_top_per_site.tsv
#   cis_meqtl_panel_summary.tsv, qq_<cohort>_<ctx>.{tiff,pdf}
#
# R packages: GENESIS, GWASTools, data.table, parallel
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(GENESIS)
  library(GWASTools)
  library(parallel)
})
options(stringsAsFactors = FALSE)

# === USER CONFIGURATION ===
PROJECT_ROOT <- Sys.getenv("TGC_PROJECT_ROOT", unset = "/path/to/your/project")
CIS_WINDOW  <- 1e5
FDR_THRESH  <- 0.05
AIREML_TOL  <- 1e-2
MIN_SAMPLES <- 10L
N_CORES     <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "4"))
# ===========================

COHORT  <- Sys.getenv("TGC_COHORT");  CONTEXT <- Sys.getenv("TGC_CONTEXT")
NCHUNKS <- as.integer(Sys.getenv("TGC_NCHUNKS", unset = "1"))
CHUNK   <- Sys.getenv("TGC_CHUNK", unset = "1")
stopifnot(COHORT %in% c("BREEDING", "NATURAL"), CONTEXT %in% c("CpG", "CHG", "CHH"))
IN_C  <- file.path(PROJECT_ROOT, "RESULTS/JOINT/MEQTL/INPUTS", COHORT)
IN_P  <- file.path(IN_C, CONTEXT)
OUT_P <- file.path(PROJECT_ROOT, "RESULTS/JOINT/MEQTL/GENESIS", COHORT, CONTEXT)
dir.create(file.path(OUT_P, "chunks"), recursive = TRUE, showWarnings = FALSE)
tag <- paste0(tolower(COHORT), "_", tolower(CONTEXT))
msg <- function(...) cat("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n", sep = "")

lambda_gc <- function(p) {
  p <- p[is.finite(p) & p > 0 & p <= 1]
  median(qchisq(p, df = 1, lower.tail = FALSE)) / qchisq(0.5, df = 1)
}

shared <- fread(file.path(IN_P, "shared_sample_ids.tsv"))$sample_id
sites  <- fread(file.path(IN_P, "methylation_site_annot.tsv"))
va     <- fread(file.path(IN_C, "snp_variant_annot.tsv"))

# Assign chromosomes to chunks, balancing the number of sites
chr_sites <- sites[, .N, by = chr][order(-N)]
chr_sites[, chunk := 0L]
load <- integer(NCHUNKS)
for (i in seq_len(nrow(chr_sites))) { k <- which.min(load); chr_sites$chunk[i] <- k; load[k] <- load[k] + chr_sites$N[i] }

if (CHUNK != "combine") {
  k <- as.integer(CHUNK)
  my_chr <- chr_sites[chunk == k, chr]
  msg(COHORT, "/", CONTEXT, " chunk ", k, "/", NCHUNKS, ": ", length(my_chr), " chromosomes, ",
      sum(sites$chr %in% my_chr), " sites, ", N_CORES, " cores")

  M    <- readRDS(file.path(IN_P, "methylation_mvalues_matrix.rds"))[shared, , drop = FALSE]
  pcs  <- fread(file.path(IN_P, "pcs_shared.tsv"))
  geno <- readRDS(file.path(IN_C, "genotypes.rds"))[, shared, drop = FALSE]
  stopifnot(identical(colnames(M), sites$methyl_loc), identical(pcs$sample_id, shared))

  # GWASTools genotype container (SNPs x samples, integer 0/1/2)
  va[, snp_int := .I]
  chr_levels <- sort(unique(va$chr))
  scan_int <- seq_along(shared)
  reader <- MatrixGenotypeReader(genotype = geno, snpID = va$snp_int,
                                 chromosome = match(va$chr, chr_levels),
                                 position = va$pos, scanID = scan_int)
  genoData <- GenotypeData(reader, scanAnnot = ScanAnnotationDataFrame(data.frame(scanID = scan_int)))
  rm(geno); gc()

  pc_cols <- grep("^PC[0-9]+$", names(pcs), value = TRUE)
  pheno_base <- data.frame(sample.id = scan_int, pcs[, ..pc_cols])
  grm_full <- readRDS(file.path(IN_P, "grm_shared.rds"))[shared, shared]
  prep_grm <- function(G) { G <- G[shared, shared]; dimnames(G) <- list(scan_int, scan_int); diag(G) <- diag(G) + 1e-4; G }
  grm_full <- prep_grm(grm_full)

  out <- list()
  for (ch in my_chr) {
    grm <- grm_full
    if (COHORT == "BREEDING") {
      f <- file.path(IN_C, "loco_grm", sprintf("loco_grm_excl_%s.rds", ch))
      if (file.exists(f)) grm <- prep_grm(readRDS(f))
    }
    s_chr <- sites[chr == ch]
    v_chr <- va[chr == ch]
    if (!nrow(v_chr)) next
    res_chr <- mclapply(seq_len(nrow(s_chr)), function(j) {
      cis <- v_chr[abs(pos - s_chr$pos[j]) <= CIS_WINDOW]
      if (!nrow(cis)) return(NULL)
      y <- M[, s_chr$methyl_loc[j]]
      if (sum(!is.na(y)) < MIN_SAMPLES) return(NULL)
      ph <- pheno_base; ph$y <- y
      nm <- tryCatch(fitNullModel(ph, outcome = "y", covars = pc_cols, cov.mat = grm,
                                  family = "gaussian", AIREML.tol = AIREML_TOL, verbose = FALSE),
                     error = function(e) NULL)
      if (is.null(nm)) return(data.table(site = s_chr$methyl_loc[j], null_failed = TRUE))
      it <- GenotypeBlockIterator(genoData, snpInclude = cis$snp_int)
      a <- tryCatch(assocTestSingle(it, null.model = nm, verbose = FALSE), error = function(e) NULL)
      if (is.null(a) || !nrow(a)) return(NULL)
      a <- as.data.table(a)
      data.table(snp = as.character(va$snp_id[a$variant.id]), site = s_chr$methyl_loc[j],
                 beta = a$Est, beta_se = a$Est.SE, statistic = a$Score.Stat, pvalue = a$Score.pval,
                 freq = a$freq, MAC = a$MAC, n_obs = a$n.obs, null_failed = FALSE)
    }, mc.cores = N_CORES, mc.preschedule = TRUE)
    out[[ch]] <- rbindlist(res_chr, fill = TRUE)
  }
  res <- rbindlist(out, fill = TRUE)
  saveRDS(res, file.path(OUT_P, "chunks", sprintf("chunk_%02d_of_%02d.rds", k, NCHUNKS)))
  msg("chunk ", k, " done: ", nrow(res[null_failed == FALSE]), " pairs; ",
      sum(res$null_failed, na.rm = TRUE), " sites with failed null model")
  if (NCHUNKS > 1) quit(save = "no")
}

# ---- combine (or single-chunk run) ----
files <- list.files(file.path(OUT_P, "chunks"), pattern = sprintf("_of_%02d\\.rds$", NCHUNKS), full.names = TRUE)
stopifnot(length(files) == NCHUNKS)
res <- rbindlist(lapply(files, readRDS), fill = TRUE)
n_null_failed <- sum(res$null_failed, na.rm = TRUE)
res <- res[null_failed == FALSE][, null_failed := NULL]
res[, p_FDR := p.adjust(pvalue, method = "BH")]
res[, `:=`(snp_chr = va$chr[match(snp, as.character(va$snp_id))],
           snp_pos = va$pos[match(snp, as.character(va$snp_id))],
           site_pos = sites$pos[match(site, sites$methyl_loc)])]
res[, distance := abs(site_pos - snp_pos)]
stopifnot(max(res$distance) <= CIS_WINDOW)
setcolorder(res, c("snp", "site", "snp_chr", "snp_pos", "site_pos", "distance",
                   "beta", "beta_se", "statistic", "pvalue", "p_FDR", "freq", "MAC", "n_obs"))

lam <- lambda_gc(res$pvalue)
sig <- res[p_FDR < FDR_THRESH][order(p_FDR)]
top <- res[res[, .I[which.min(pvalue)], by = site]$V1]
msg("cis pairs: ", nrow(res), " | FDR < ", FDR_THRESH, ": ", nrow(sig), " | lambda = ", round(lam, 4))

saveRDS(res, file.path(OUT_P, "cis_meqtl_all_results.rds"))
fwrite(res, file.path(OUT_P, "cis_meqtl_all_results.tsv.gz"), sep = "\t")
fwrite(sig, file.path(OUT_P, "cis_meqtl_significant.tsv"), sep = "\t")
fwrite(top, file.path(OUT_P, "cis_meqtl_top_per_site.tsv"), sep = "\t")
fwrite(data.table(tool = "GENESIS", cohort = COHORT, context = CONTEXT,
                  n_samples = length(shared), n_sites = nrow(sites), n_snps = nrow(va),
                  n_sites_with_cis = uniqueN(res$site), n_snps_with_cis = uniqueN(res$snp),
                  n_sites_null_failed = n_null_failed,
                  n_pairs = nrow(res), n_fdr05 = nrow(sig),
                  n_fdr1e10 = sum(res$p_FDR < 1e-10), lambda = round(lam, 4)),
       file.path(OUT_P, "cis_meqtl_panel_summary.tsv"), sep = "\t")

pv <- sort(res$pvalue[res$pvalue > 0])
qq <- function() {
  plot(-log10(ppoints(length(pv))), -log10(pv), pch = 20, cex = 0.4, col = "#e34a33",
       xlab = expression(Expected ~ -log[10](italic(p))), ylab = expression(Observed ~ -log[10](italic(p))))
  abline(0, 1, col = "red", lty = 2)
  mtext(sprintf("lambda = %.3f", lam), side = 3, adj = 0, line = 0.3)
}
tiff(file.path(OUT_P, sprintf("qq_%s.tiff", tag)), width = 2400, height = 2400, res = 300, compression = "lzw"); qq(); dev.off()
pdf(file.path(OUT_P, sprintf("qq_%s.pdf", tag)), width = 8, height = 8); qq(); dev.off()
msg("Step 12ab2 done")
```

---

### `12ab2.tgc.joint.genesis.mapping.sh`

```bash
#!/bin/bash
#SBATCH --job-name=JOINT.12ab2.GENESIS
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --array=0-47
#SBATCH --output=LOGS/12ab2/%x_%A_%a.out
#SBATCH --error=LOGS/12ab2/%x_%A_%a.err
#-------------------------------------------------------------------------------
# TreeGeneClimate (TGC) — JOINT
# Step 12ab2: GENESIS cis-meQTL mapping
#   Mapping:  48 array tasks = 6 panels (cohort x context) x 8 chunks of
#             chromosomes/scaffolds:   sbatch 12ab2...sh
#   Combine:  one task per panel (FDR over the whole panel), after mapping:
#             sbatch --array=0-5 --export=ALL,TGC_CHUNK=combine 12ab2...sh
# Create LOGS/12ab2 before submitting.
#-------------------------------------------------------------------------------
set -euo pipefail

# === USER CONFIGURATION ===
export TGC_PROJECT_ROOT="${TGC_PROJECT_ROOT:-/path/to/your/project}"
NCHUNKS=8
module purge 2>/dev/null || true
module load gcc/14.2.0 r/4.5.2 2>/dev/null || true
export LC_ALL=C.UTF-8 LANG=C.UTF-8
# ===========================

cohorts=(BREEDING BREEDING BREEDING NATURAL NATURAL NATURAL)
contexts=(CpG CHG CHH CpG CHG CHH)
i="${SLURM_ARRAY_TASK_ID:-0}"
export TGC_NCHUNKS="${NCHUNKS}"
if [[ "${TGC_CHUNK:-}" == "combine" ]]; then
  p="${i}"
else
  p=$(( i / NCHUNKS )); export TGC_CHUNK=$(( i % NCHUNKS + 1 ))
fi
export TGC_COHORT="${cohorts[$p]}" TGC_CONTEXT="${contexts[$p]}"
echo "[$(date)] ${TGC_COHORT} ${TGC_CONTEXT} chunk ${TGC_CHUNK}/${NCHUNKS}"
Rscript --vanilla "${TGC_PROJECT_ROOT}/SCRIPTS/JOINT/12ab2.tgc.joint.genesis.mapping.R"
```

---

### `13ab.tgc.joint.meqtl.combined.results.R`

```r
#!/usr/bin/env Rscript
############################################################
# TreeGeneClimate (TGC) — JOINT ECS + TMS
# Step 13ab: Combined cis-meQTL results
#
#   - Panel summary for both tools: tested pairs/sites/SNPs,
#     significant pairs, unique sites and unique SNPs at FDR < 1e-10, lambda
#   - Significant pairs per tool and cohort (FDR < 1e-10, all contexts)
#   - Robust pairs: the same SNP-site pair at FDR < 1e-10 in BOTH GENESIS and
#     MatrixEQTL (pair-level replication), with genomic coordinates
#   - Robust counts per cohort and context
#   - QQ panel: a, b MatrixEQTL and c, d GENESIS (breeding, natural), with the three
#     methylation contexts in each plot; all pairs with p < 1e-3 and a random sample
#     of 20,000 of the remaining pairs per tool, cohort and context
#
# INPUTS
#   RESULTS/JOINT/MEQTL/{MATRIXEQTL,GENESIS}/<COHORT>/<CTX>/cis_meqtl_all_results.rds
#   RESULTS/JOINT/MEQTL/{MATRIXEQTL,GENESIS}/<COHORT>/<CTX>/cis_meqtl_panel_summary.tsv
#
# OUTPUTS (RESULTS/JOINT/MEQTL/COMBINED/)
#   meqtl_panel_summary.tsv
#   sig_fdr1e10_<tool>_<cohort>.tsv
#   robust_pairs_<cohort>.tsv
#   robust_context_summary.tsv
#   qq_panel.{tiff,pdf}
#   qq_panel_values.tsv  (values shown in the QQ panel)
#
# R packages: data.table
############################################################

suppressPackageStartupMessages(library(data.table))
options(stringsAsFactors = FALSE)

# === USER CONFIGURATION ===
PROJECT_ROOT <- Sys.getenv("TGC_PROJECT_ROOT", unset = "/path/to/your/project")
FDR_STRICT <- 1e-10
# ===========================

MEQTL <- file.path(PROJECT_ROOT, "RESULTS/JOINT/MEQTL")
OUT   <- file.path(MEQTL, "COMBINED")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
COHORTS  <- c("BREEDING", "NATURAL")
CONTEXTS <- c("CpG", "CHG", "CHH")
TOOLS    <- c(MatrixEQTL = "MATRIXEQTL", GENESIS = "GENESIS")
msg <- function(...) cat("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n", sep = "")
set.seed(1)   # thinning of points in the QQ plots

panel_rows <- list(); sig <- list(); qq_data <- list()
for (tool in names(TOOLS)) for (cohort in COHORTS) for (ctx in CONTEXTS) {
  d <- file.path(MEQTL, TOOLS[[tool]], cohort, ctx)
  res <- readRDS(file.path(d, "cis_meqtl_all_results.rds"))
  ps  <- fread(file.path(d, "cis_meqtl_panel_summary.tsv"))
  s <- res[p_FDR < FDR_STRICT]
  s[, `:=`(tool = tool, cohort = cohort, context = ctx)]
  sig[[paste(tool, cohort, ctx)]] <- s
  panel_rows[[paste(tool, cohort, ctx)]] <- data.table(
    tool = tool, cohort = cohort, context = ctx, n_samples = ps$n_samples,
    sites_tested = ps$n_sites_with_cis, snps_tested = ps$n_snps_with_cis, pairs_tested = nrow(res),
    sig_pairs = nrow(s), sig_sites = uniqueN(s$site), sig_snps = uniqueN(s$snp), lambda = ps$lambda)
  # QQ data: all points with p < 1e-3, a 20,000-point sample of the rest
  p <- res$pvalue[res$pvalue > 0]
  n <- length(p); p <- sort(p); e <- ppoints(n)
  keep <- which(p < 1e-3)
  rest <- setdiff(seq_len(n), keep)
  if (length(rest) > 2e4) rest <- sort(sample(rest, 2e4))
  idx <- sort(c(keep, rest))
  qq_data[[paste(tool, cohort, ctx)]] <- list(o = -log10(p[idx]), e = -log10(e[idx]))
  msg(tool, " ", cohort, " ", ctx, ": ", nrow(res), " pairs; ", nrow(s), " at FDR < ", FDR_STRICT)
  rm(res); gc()
}
panel <- rbindlist(panel_rows)
fwrite(panel, file.path(OUT, "meqtl_panel_summary.tsv"), sep = "\t")
sig <- rbindlist(sig, fill = TRUE)
for (tl in names(TOOLS)) for (co in COHORTS)
  fwrite(sig[tool == tl & cohort == co], file.path(OUT, sprintf("sig_fdr1e10_%s_%s.tsv", tolower(tl), tolower(co))), sep = "\t")

# ---- robust pairs ----
keep_cols <- c("snp", "site", "context", "snp_chr", "snp_pos", "site_pos", "distance", "beta", "statistic", "pvalue", "p_FDR")
rob_summary <- list()
for (co in COHORTS) {
  m <- sig[tool == "MatrixEQTL" & cohort == co, ..keep_cols]
  g <- sig[tool == "GENESIS"    & cohort == co, ..keep_cols]
  rob <- merge(m, g, by = c("snp", "site", "context", "snp_chr", "snp_pos", "site_pos", "distance"),
               suffixes = c("_MatrixEQTL", "_GENESIS"))
  rob[, site_chr := sub(":.*$", "", site)]
  setcolorder(rob, c("context", "snp", "snp_chr", "snp_pos", "site", "site_chr", "site_pos", "distance"))
  setorder(rob, context, snp_chr, snp_pos)
  rob[, cohort := co]
  fwrite(rob, file.path(OUT, sprintf("robust_pairs_%s.tsv", tolower(co))), sep = "\t")
  rob_summary[[co]] <- rbind(
    rob[, .(pairs = .N, sites = uniqueN(site), snps = uniqueN(snp)), by = context],
    rob[, .(context = "All", pairs = .N, sites = uniqueN(site), snps = uniqueN(snp))])[, cohort := co]
  msg(co, ": ", nrow(rob), " robust pairs, ", uniqueN(rob$site), " sites, ", uniqueN(rob$snp), " SNPs")
}
rob_summary <- rbindlist(rob_summary)
setcolorder(rob_summary, c("cohort", "context"))
fwrite(rob_summary, file.path(OUT, "robust_context_summary.tsv"), sep = "\t")
print(panel); print(rob_summary)

# ---- QQ panel: a, b MatrixEQTL and c, d GENESIS (breeding, natural), the three contexts in each plot ----
CTX_COLORS <- c(CpG = "#D81B60", CHG = "#1F4E9C", CHH = "#17BECF")
draw_qq <- function() {
  layout(matrix(c(1, 2, 3, 4, 5, 5), nrow = 3, byrow = TRUE), heights = c(1, 1, 0.10))
  par(cex = 0.83, mar = c(4.2, 4.5, 2, 1))
  i <- 0
  for (tool in names(TOOLS)) for (cohort in COHORTS) {
    i <- i + 1
    q <- lapply(CONTEXTS, function(ctx) qq_data[[paste(tool, cohort, ctx)]])
    plot(range(unlist(lapply(q, `[[`, "e"))), range(unlist(lapply(q, `[[`, "o"))), type = "n",
         xlab = expression(Expected ~ -log[10](italic(p))), ylab = expression(Observed ~ -log[10](italic(p))))
    abline(0, 1, col = "grey30", lty = 2)
    for (k in seq_along(CONTEXTS)) points(q[[k]]$e, q[[k]]$o, pch = 20, cex = 0.35, col = CTX_COLORS[[CONTEXTS[k]]])
    mtext(letters[i], side = 3, adj = -0.18, line = 0.6, font = 2, cex = 0.83)
  }
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend("center", legend = CONTEXTS, pch = 19, col = CTX_COLORS[CONTEXTS], horiz = TRUE, bty = "n", cex = 1)
}
tiff(file.path(OUT, "qq_panel.tiff"), width = 20, height = 21, units = "cm", res = 300, compression = "lzw"); draw_qq(); dev.off()
pdf(file.path(OUT, "qq_panel.pdf"), width = 20 / 2.54, height = 21 / 2.54); draw_qq(); dev.off()

# Values shown in the QQ panel: expected and observed -log10(p) of each point
qq_rows <- list()
i <- 0
for (tool in names(TOOLS)) for (cohort in COHORTS) {
  i <- i + 1
  for (ctx in CONTEXTS) {
    q <- qq_data[[paste(tool, cohort, ctx)]]
    qq_rows[[paste(tool, cohort, ctx)]] <- data.table(panel = letters[i], tool = tool, cohort = cohort, context = ctx,
                                                      expected_log10p = q$e, observed_log10p = q$o)
  }
}
fwrite(rbindlist(qq_rows), file.path(OUT, "qq_panel_values.tsv"), sep = "\t")
msg("Step 13ab done")
```

---

### `14ab.tgc.joint.manhattan.plots.R`

```r
#!/usr/bin/env Rscript
############################################################
# TreeGeneClimate (TGC) — JOINT ECS + TBS
# Step 14ab: Circular Manhattan plots (GENESIS and MatrixEQTL)
#            18 cm canvas per cohort; raw p-value or BH-FDR axis
#            + Combined 2-panel figures (BREEDING | NATURAL + legend)
#
# PURPOSE
# Produces circular Manhattan plots from the cis-meQTL results.
# One plot per tool x cohort x axis mode (8 total):
#   GENESIS      BREEDING/NATURAL  x  raw-p / BH-FDR
#   MATRIXEQTL   BREEDING/NATURAL  x  raw-p / BH-FDR
# Three rings: CpG (outer), CHG (middle), CHH (inner).
# Dashed threshold line at BH-FDR < 1e-10.
#
# Rendering: 18 cm x 18 cm canvas per cohort, 36 cm for the combined panel. Each ring
# has its own radial axis; an axis with the range of the widest ring and the axis
# title are drawn in the centre of the circle. Context colours are shared with the
# other figures (magenta, blue, turquoise; no red-green pair).
#
# USAGE
#   Rscript --vanilla 14ab.tgc.joint.manhattan.plots.R
#   # or source from another script: SOURCED_14AB <- TRUE; source("14ab.R")
#
# INPUTS
#   RESULTS/JOINT/MEQTL/GENESIS/<COHORT>/<CONTEXT>/cis_meqtl_all_results.rds
#   RESULTS/JOINT/MEQTL/MATRIXEQTL/<COHORT>/<CONTEXT>/cis_meqtl_all_results.rds
#   RESULTS/JOINT/MEQTL/INPUTS/<COHORT>/snp_variant_annot.tsv
#   RESULTS/ECS/RANALYSIS/RDATA/<cohort>.imputed.snp.gds (chromosome layout)
#   DATA/METADATA/chrom.list  (optional; auto-built if absent)
#
# OUTPUTS (RESULTS/JOINT/MEQTL/FIGURES/manhattan/)
#   manhattan_circular_<tool>_<cohort>.{tiff,pdf,svg,eps}      -- raw p-value
#   manhattan_circular_<tool>_<cohort>_fdr.{tiff,pdf,svg,eps}   -- BH-FDR axis
#   legend_circular.tiff      (36 cm, for 2-panel combined figure)
#   legend_circular_wide.tiff (72 cm, for 4-panel comparison)
# OUTPUTS (RESULTS/JOINT/MEQTL/FIGURES/)
#   RESULTS/JOINT/MEQTL/FIGURES/circular_panel_genesis.{tiff,pdf,svg,png,eps}
#   RESULTS/JOINT/MEQTL/FIGURES/circular_panel_matrixeqtl.{tiff,pdf,svg,png,eps}
#   circular_panel_genesis_FDR.{tiff,pdf,svg,png,eps}
#   circular_panel_matrixeqtl_FDR.{tiff,pdf,svg,png,eps}
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(SNPRelate)
  library(gdsfmt)
  library(circlize)
  library(grid)
})

options(stringsAsFactors = FALSE)

############################################################
# 1) SETTINGS
############################################################

COHORTS        <- c("BREEDING", "NATURAL")
CTX_ORDER      <- c("CpG", "CHG", "CHH")   # outer -> inner rings
TOOLS          <- c("GENESIS", "MATRIXEQTL")
FDR_STRICT     <- 1e-10

# FDR_AXIS: when TRUE, y-axis shows -log10(BH-FDR); threshold constant at 10.
# Set before sourcing, or it defaults to FALSE.
if (!exists("FDR_AXIS")) FDR_AXIS <- FALSE

# Rendering settings
OUT_W_CM          <- 18; OUT_H_CM <- 18; OUT_RES <- 600
AXIS_GAP_FRAC     <- 0.80    # gap for the radial axes
AXIS_SECTOR_SCALE <- 1.10    # width multiplier of the axis sector (axis_panel2)
TRACK_HEIGHT      <- 0.25
CTX_COLORS        <- c(CpG = "#D81B60", CHG = "#1F4E9C", CHH = "#17BECF")
AXIS_TICK_CEX     <- 0.55    # radial axis tick labels
AXIS_TITLE_CEX    <- 0.85    # axis title, drawn in the centre of the circle
PANEL_LETTER_CEX  <- 1.49    # panel letters
LEGEND_CEX        <- 1.30    # context legend
CHR_LABEL_CEX     <- 1.02    # chromosome labels ("Un" fits its sector)

# Combined panel dimensions (BREEDING | NATURAL side-by-side + legend strip)
LEG_H_CM  <- 2.0                          # legend strip height
COMB_W_CM <- OUT_W_CM * 2                 # 36 cm total width
COMB_H_CM <- OUT_W_CM + LEG_H_CM          # 20 cm total height

############################################################
# 2) PATHS  -- set PROJECT_ROOT and RDATA_DIR to your paths
############################################################

# === USER CONFIGURATION ===
# Set TGC_PROJECT_ROOT as an environment variable, or edit the fallback path below
PROJECT_ROOT <- Sys.getenv("TGC_PROJECT_ROOT", unset = "/path/to/your/project")
# ===========================

RESULT_ROOTS <- list(
  GENESIS    = file.path(PROJECT_ROOT, "RESULTS", "JOINT", "MEQTL", "GENESIS"),
  MATRIXEQTL = file.path(PROJECT_ROOT, "RESULTS", "JOINT", "MEQTL", "MATRIXEQTL")
)

ANNOT_ROOT <- file.path(PROJECT_ROOT, "RESULTS", "JOINT", "MEQTL", "INPUTS")

OUT_MAN       <- file.path(PROJECT_ROOT, "RESULTS", "JOINT", "MEQTL", "FIGURES", "manhattan")
PANEL_DIR     <- file.path(PROJECT_ROOT, "RESULTS", "JOINT", "MEQTL", "FIGURES")
LOG_DIR       <- file.path(PROJECT_ROOT, "RESULTS", "JOINT", "MEQTL", "LOGS")
for (d in c(OUT_MAN, PANEL_DIR, LOG_DIR))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

LOGFILE <- file.path(LOG_DIR, "step14ab.log")
if (file.exists(LOGFILE)) file.remove(LOGFILE)

# Imputed genotype GDS files (Step 8a), used only for the contig layout
RDATA_DIR <- file.path(PROJECT_ROOT, "RESULTS", "ECS", "RANALYSIS", "RDATA")
GDS_FILES <- list(
  BREEDING = file.path(RDATA_DIR, "breeding.imputed.snp.gds"),
  NATURAL  = file.path(RDATA_DIR, "natural.imputed.snp.gds")
)
CHROM_MAP_FILE <- file.path(PROJECT_ROOT, "DATA", "METADATA", "chrom.list")

############################################################
# 3) HELPERS
############################################################

log_msg <- function(...) {
  txt <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ",
                paste0(..., collapse = ""))
  cat(txt, "\n")
  cat(txt, "\n", file = LOGFILE, append = TRUE)
}
sep_line <- function() log_msg(paste(rep("=", 70), collapse = ""))

MAGICK_BIN <- Sys.which("magick")
if (!nzchar(MAGICK_BIN)) MAGICK_BIN <- Sys.which("convert")
HAS_MAGICK <- nzchar(MAGICK_BIN)

run_magick <- function(args_str) {
  cmd <- paste(shQuote(MAGICK_BIN), args_str)
  ret <- system(cmd, ignore.stdout = TRUE, ignore.stderr = FALSE)
  if (ret != 0L) stop("magick failed:\n  ", cmd)
}

check_nonempty <- function(path, min_bytes = 5000L) {
  file.exists(path) && file.info(path)$size >= min_bytes
}

pval_column <- function(res_dt) {
  candidates <- c("pvalue", "pval", "Score.pval",
                  grep("pval$", names(res_dt), value = TRUE))
  candidates <- candidates[candidates %in% names(res_dt)]
  if (!length(candidates)) stop("No p-value column found")
  candidates[1]
}

find_pval_at_fdr <- function(pvals, fdr_thresh) {
  pvals <- sort(pvals[is.finite(pvals) & !is.na(pvals) & pvals > 0])
  if (!length(pvals)) return(NA_real_)
  padj <- p.adjust(pvals, method = "BH")
  idx  <- which(padj < fdr_thresh)
  if (!length(idx)) return(NA_real_)
  pvals[max(idx)]
}

read_gds_snppos <- function(gds_path) {
  if (!file.exists(gds_path)) return(NULL)
  gf <- snpgdsOpen(gds_path); on.exit(snpgdsClose(gf))
  data.table(
    snp = as.character(read.gdsn(index.gdsn(gf, "snp.id"))),
    chr = as.character(read.gdsn(index.gdsn(gf, "snp.chromosome"))),
    pos = as.integer (read.gdsn(index.gdsn(gf, "snp.position")))
  )
}

build_contig_lengths <- function(gds_paths) {
  all_dt <- rbindlist(lapply(gds_paths, read_gds_snppos), fill = TRUE)
  if (!nrow(all_dt)) stop("No contig positions from GDS files")
  cl <- all_dt[, .(len = max(pos, na.rm = TRUE)), by = chr]
  cl[!is.na(chr)]
}

read_chrom_map <- function(path) {
  if (!file.exists(path)) return(NULL)
  cm <- fread(path, data.table = FALSE)
  if (!all(c("Real_chrom", "ChromN") %in% colnames(cm))) return(NULL)
  cm$Real_chrom <- as.character(cm$Real_chrom)
  cm$ChromN     <- as.character(cm$ChromN)
  cm$ChromN[cm$ChromN == "Chr13"] <- "ChrUn"
  cm
}

build_default_chrom_map <- function(contig_len_dt) {
  dt <- copy(contig_len_dt); setorder(dt, -len)
  n  <- min(12L, nrow(dt))
  dt[, ChromN := ifelse(seq_len(.N) <= n, sprintf("Chr%02d", seq_len(.N)), "ChrUn")]
  data.frame(Real_chrom = dt$chr, ChromN = dt$ChromN, stringsAsFactors = FALSE)
}

prepare_layout <- function(contig_len_dt, chrom_map_df, gap_frac = AXIS_GAP_FRAC) {
  cm_dt <- as.data.table(chrom_map_df)
  dt    <- merge(contig_len_dt, cm_dt, by.x = "chr", by.y = "Real_chrom",
                 all.x = TRUE, sort = FALSE)
  dt[is.na(ChromN), ChromN := "ChrUn"]
  desired <- c("ChrUn", sprintf("Chr%02d", 1:12))
  lvls    <- c(intersect(desired, unique(dt$ChromN)),
               sort(setdiff(unique(dt$ChromN), desired)))
  dt[, ChromN := factor(ChromN, levels = lvls)]; setorder(dt, ChromN, chr)
  dt[, offset_within_chromN := {
    tmp <- cumsum(len); if (.N == 0) numeric(0) else c(0, tmp[-.N])
  }, by = ChromN]
  clen <- dt[, .(chromN_len = sum(len)), by = ChromN]
  clen[, ChromN := factor(ChromN, levels = lvls)]; setorder(clen, ChromN)
  gap_len <- max(clen$chromN_len) * gap_frac * 0.5 * AXIS_SECTOR_SCALE
  idx_un  <- which(lvls == "ChrUn")
  if (length(idx_un)) {
    lvls_g <- append(lvls, "axis_panel2", after = idx_un)
    sl_lst <- as.list(setNames(clen$chromN_len, as.character(clen$ChromN)))
    sl_lst <- append(sl_lst, list(axis_panel2 = gap_len), after = idx_un)
  } else {
    lvls_g <- c("axis_panel2", lvls)
    sl_lst <- c(list(axis_panel2 = gap_len),
                as.list(setNames(clen$chromN_len, as.character(clen$ChromN))))
  }
  sl <- unlist(sl_lst[lvls_g])
  if ("Chr01" %in% names(sl)) {
    idx <- which(names(sl) == "Chr01")
    nms <- names(sl); sl <- sl[c(nms[idx:length(nms)], nms[seq_len(idx - 1)])]
  }
  list(contig_info = dt, chrom_levels = lvls, sector_lens = sl)
}

map_to_chromN <- function(res_dt, contig_info, chrom_map_df, pval_col) {
  dt <- copy(res_dt)
  if (pval_col != "pval" && pval_col %in% names(dt)) dt[, pval := get(pval_col)]
  if (!"snp_pos" %in% names(dt)) return(data.table())
  if ("snp_chr" %in% names(dt)) {
    dt[, chr := snp_chr]
  } else {
    dt[, chr := NA_character_]
  }
  dt[, pos := as.numeric(snp_pos)]
  dt <- dt[!is.na(pos) & is.finite(pval) & pval > 0 & !is.na(chr)]
  if (!nrow(dt)) return(data.table())
  cm_dt <- as.data.table(chrom_map_df)
  dt[, Real_chrom := chr]
  dt <- merge(dt, cm_dt, by.x = "Real_chrom", by.y = "Real_chrom",
              all.x = TRUE, sort = FALSE)
  dt[is.na(ChromN), ChromN := "ChrUn"]
  setkey(contig_info, chr)
  out <- contig_info[dt, on = .(chr)]
  out[is.na(offset_within_chromN), offset_within_chromN := 0]
  out[, pos_within_chromN := offset_within_chromN + pos]
  # FDR_AXIS: use pre-computed BH-FDR column if available, else raw p
  if (FDR_AXIS && "fdr__" %in% names(out)) {
    out[, logp := -log10(fdr__)]
  } else {
    out[, logp := -log10(pval)]
  }
  out
}

############################################################
# 4) DATA LOADING (separated from plotting to allow multiple device types)
############################################################

load_panel_data <- function(tool, cohort) {
  ctx_colors <- CTX_COLORS[CTX_ORDER]

  panel_data       <- list()
  threshold_logp   <- list()
  max_logp_per_ctx <- list()

  for (ctx in CTX_ORDER) {
    rpath <- file.path(RESULT_ROOTS[[tool]], cohort, ctx, "cis_meqtl_all_results.rds")
    if (!file.exists(rpath)) {
      log_msg("  Missing RDS: ", rpath)
      panel_data[[ctx]] <- data.table()
      threshold_logp[[ctx]] <- NA_real_
      max_logp_per_ctx[[ctx]] <- 10
      next
    }

    log_msg("  Loading: ", tool, " / ", cohort, " / ", ctx,
            " [FDR_AXIS=", FDR_AXIS, "]")
    res_dt   <- as.data.table(readRDS(rpath))
    pval_col <- tryCatch(pval_column(res_dt), error = function(e) NA_character_)
    if (is.na(pval_col)) {
      log_msg("  No p-value column in ", rpath)
      panel_data[[ctx]] <- data.table()
      threshold_logp[[ctx]] <- NA_real_
      max_logp_per_ctx[[ctx]] <- 10
      next
    }
    log_msg("  p-value column: '", pval_col, "' | nrow=", nrow(res_dt))

    pvals <- res_dt[[pval_col]]

    # Compute BH-FDR on the full context test set before position filtering.
    # NOTE: FDR_AXIS=TRUE on NATURAL/CHH (~500M tests) requires ~96 GB RAM.
    if (FDR_AXIS) {
      log_msg("  Computing BH-FDR on ", length(pvals), " tests...")
      fdrs <- p.adjust(pvals, method = "BH")
      res_dt[, fdr__ := fdrs]
      rm(fdrs); gc()
      log_msg("  BH-FDR done")
    }

    # Rebuild snp_chr and snp_pos from snp_variant_annot.rds (authoritative coords).
    for (col_drop in intersect(c("chr", "snp_chr", "pos", "snp_pos"), names(res_dt)))
      res_dt[, (col_drop) := NULL]

    var_path <- file.path(ANNOT_ROOT, cohort, "snp_variant_annot.tsv")
    if (file.exists(var_path) && "snp" %in% names(res_dt)) {
      va <- fread(var_path)
      log_msg("  snp_variant_annot cols: ", paste(names(va), collapse = ", "))
      res_dt[, snp   := as.character(snp)]
      va[,   snp_id := as.character(snp_id)]
      if ("pos" %in% names(va)) {
        res_dt <- merge(res_dt,
                        va[, .(snp_id, snp_chr = chr, snp_pos = as.numeric(pos))],
                        by.x = "snp", by.y = "snp_id", all.x = TRUE)
        log_msg("  snp_pos matched: ", sum(!is.na(res_dt$snp_pos)), " / ", nrow(res_dt))
      } else {
        log_msg("  WARNING: snp_variant_annot has no 'pos' column")
        res_dt <- merge(res_dt, va[, .(snp_id, snp_chr = chr)],
                        by.x = "snp", by.y = "snp_id", all.x = TRUE)
      }
      n_matched <- sum(!is.na(res_dt$snp_chr))
      log_msg("  snp_chr matched: ", n_matched, " / ", nrow(res_dt),
              " (", round(100 * n_matched / max(1, nrow(res_dt)), 1), "%)")
    } else {
      if (!file.exists(var_path))
        log_msg("  WARNING: snp_variant_annot.tsv not found: ", var_path)
      if (!"snp" %in% names(res_dt))
        log_msg("  WARNING: no 'snp' column -- cannot add snp_chr")
      res_dt[, snp_chr := NA_character_]
    }

    # Threshold logp
    if (FDR_AXIS) {
      threshold_logp[[ctx]] <- -log10(FDR_STRICT)   # constant = 10 for all contexts
    } else {
      p_thr <- find_pval_at_fdr(pvals, FDR_STRICT)
      threshold_logp[[ctx]] <- if (is.na(p_thr)) NA_real_ else -log10(p_thr)
    }
    log_msg("  Threshold logp: ",
            if (is.na(threshold_logp[[ctx]])) "none" else round(threshold_logp[[ctx]], 2))

    mapped <- tryCatch(
      map_to_chromN(res_dt, contig_info_dt, chrom_map_df, pval_col),
      error = function(e) {
        log_msg("  map_to_chromN ERROR: ", conditionMessage(e))
        data.table()
      }
    )
    log_msg("  Mapped points: ", nrow(mapped))
    if (nrow(mapped) == 0 && "snp_chr" %in% names(res_dt)) {
      log_msg("  DIAG snp_chr sample  : ",
              paste(head(unique(na.omit(res_dt$snp_chr)), 5), collapse = ", "))
      log_msg("  DIAG contig_info chr : ",
              paste(head(unique(contig_info_dt$chr), 5), collapse = ", "))
    }

    max_logp_per_ctx[[ctx]] <- if (nrow(mapped) > 0 && any(is.finite(mapped$logp))) {
      max(mapped$logp[is.finite(mapped$logp)])
    } else { 10 }

    panel_data[[ctx]] <- mapped
    rm(res_dt, pvals); gc()
  }

  list(panel_data       = panel_data,
       threshold_logp   = threshold_logp,
       max_logp_per_ctx = max_logp_per_ctx,
       ctx_colors       = ctx_colors)
}

############################################################
# 5) CIRCULAR PLOT FUNCTION
############################################################

draw_circular_manhattan <- function(dat, tool, cohort, out_path = NULL,
                                    open_device = TRUE, label_cex = PANEL_LETTER_CEX) {

  panel_data       <- dat$panel_data
  threshold_logp   <- dat$threshold_logp
  max_logp_per_ctx <- dat$max_logp_per_ctx
  ctx_colors       <- dat$ctx_colors

  if (open_device) {
    if (grepl("\\.pdf$", out_path, ignore.case = TRUE)) {
      pdf(out_path, width = OUT_W_CM / 2.54, height = OUT_H_CM / 2.54)
    } else if (grepl("\\.svg$", out_path, ignore.case = TRUE)) {
      svg(out_path, width = OUT_W_CM / 2.54, height = OUT_H_CM / 2.54)
    } else if (grepl("\\.eps$", out_path, ignore.case = TRUE)) {
      cairo_ps(out_path, width = OUT_W_CM / 2.54, height = OUT_H_CM / 2.54)
    } else if (grepl("\\.png$", out_path, ignore.case = TRUE)) {
      png(out_path, width = OUT_W_CM, height = OUT_H_CM, units = "cm", res = 150)
    } else {
      tiff(out_path, width = OUT_W_CM, height = OUT_H_CM,
           units = "cm", res = OUT_RES, compression = "lzw")
    }
    on.exit({ circos.clear(); dev.off() }, add = TRUE)
  } else {
    on.exit(circos.clear(), add = TRUE)
  }

  par(mar = c(0.5, 2, 1.0, 2))   # reduced: less space above A/B label, below legend

  circos.par(
    cell.padding            = c(0, 0, 0, 0),
    track.margin            = c(0.005, 0.005),
    gap.after               = rep(0.5, length(sector_lens)),
    start.degree            = 90,
    clock.wise              = TRUE,
    points.overflow.warning = FALSE
  )
  circos.initialize(factors = names(sector_lens),
                    xlim    = cbind(rep(0, length(sector_lens)),
                                   unname(sector_lens)))

  # Chromosome label track
  circos.track(ylim = c(0, 1), track.height = 0.05,
               bg.border = NA, bg.col = "white",
    panel.fun = function(x, y) {
      sn <- CELL_META$sector.index
      if (sn == "axis_panel2") return(invisible())
      lbl <- if (sn == "ChrUn") "Un" else sn
      circos.text(CELL_META$xcenter, 0.5, lbl,
                  facing = "bending.inside", niceFacing = TRUE,
                  cex = CHR_LABEL_CEX, col = "black")
    })

  yaxis_label <- if (FDR_AXIS) expression(-log[10](FDR)) else expression(-log[10](p))

  # Explicit per-track calls (not a loop) to avoid R closure-in-loop issues.
  add_data_track <- function(dat_ctx, col_ctx, thr_logp, ylim_top) {
    circos.track(ylim = c(0, ylim_top), track.height = TRACK_HEIGHT,
                 bg.border = ifelse(get.all.sector.index() == "axis_panel2", NA, "grey55"), bg.col = "white",
      panel.fun = function(x, y) {
        sn <- CELL_META$sector.index

        if (sn == "axis_panel2") {
          ax_x     <- CELL_META$xlim[2]
          at_vals  <- pretty(c(0, ylim_top), n = 3)
          at_vals  <- at_vals[at_vals >= 0 & at_vals <= ylim_top]
          tick_len <- CELL_META$xrange * 0.16

          circos.lines(c(ax_x, ax_x), c(0, ylim_top), lwd = 0.7, straight = TRUE)
          for (v in at_vals) {
            circos.lines(c(ax_x - tick_len, ax_x), c(v, v), lwd = 0.6, straight = TRUE)
            if (v > 0) circos.text(ax_x - tick_len * 2.0, v,
                        labels = as.character(v),
                        adj = c(1, 0.5), cex = AXIS_TICK_CEX, font = 1, col = "black",
                        facing = "bending.inside", niceFacing = TRUE)
          }
          return(invisible())
        }

        if (!nrow(dat_ctx)) return(invisible())
        sub <- dat_ctx[ChromN == sn & !is.na(pos_within_chromN) & !is.na(logp)]
        if (!nrow(sub)) return(invisible())

        circos.points(sub$pos_within_chromN, sub$logp,
                      pch = 16, cex = 0.22,
                      col = adjustcolor(col_ctx, alpha.f = 0.4))

        if (is.finite(thr_logp) && thr_logp > 0 && thr_logp <= ylim_top)
          circos.lines(c(0, CELL_META$xlim[2]),
                       c(thr_logp, thr_logp),
                       lty = 2, lwd = 0.8, col = "black")
      })
  }

  ymax_cpg <- max_logp_per_ctx[["CpG"]]; if (!is.finite(ymax_cpg) || ymax_cpg <= 0) ymax_cpg <- 10
  ymax_chg <- max_logp_per_ctx[["CHG"]]; if (!is.finite(ymax_chg) || ymax_chg <= 0) ymax_chg <- 10
  ymax_chh <- max_logp_per_ctx[["CHH"]]; if (!is.finite(ymax_chh) || ymax_chh <= 0) ymax_chh <- 10

  add_data_track(panel_data[["CpG"]], ctx_colors[["CpG"]], threshold_logp[["CpG"]], ymax_cpg * 1.05)
  add_data_track(panel_data[["CHG"]], ctx_colors[["CHG"]], threshold_logp[["CHG"]], ymax_chg * 1.05)
  add_data_track(panel_data[["CHH"]], ctx_colors[["CHH"]], threshold_logp[["CHH"]], ymax_chh * 1.05)

  # Axis in the centre of the circle: range of the widest ring, with the axis title
  ymax_all <- max(ymax_cpg, ymax_chg, ymax_chh) * 1.05
  at_c <- pretty(c(0, ymax_all), n = 3)
  at_c <- at_c[at_c >= 0 & at_c <= ymax_all]
  r_hole <- get.cell.meta.data("cell.bottom.radius", sector.index = names(sector_lens)[1], track.index = 4)
  y_c <- function(v) -0.70 * r_hole + 1.40 * r_hole * v / ymax_all
  x_c <- 0.50 * r_hole
  segments(x_c, y_c(0), x_c, y_c(ymax_all), lwd = 0.7)
  segments(x_c - 0.08 * r_hole, y_c(at_c), x_c, y_c(at_c), lwd = 0.6)
  text(x_c - 0.13 * r_hole, y_c(at_c), labels = at_c, adj = c(1, 0.5), cex = AXIS_TICK_CEX)
  text(-0.25 * r_hole, 0, yaxis_label, srt = 90, adj = c(0.5, 0.5), cex = AXIS_TITLE_CEX)

  # Panel label: A = BREEDING, B = NATURAL
  fig_region  <- par("fig")
  fig_x1 <- fig_region[1]; fig_x2 <- fig_region[2]
  fig_y1 <- fig_region[3]; fig_y2 <- fig_region[4]
  label_x_npc <- fig_x1 + 0.02 * (fig_x2 - fig_x1)
  label_y_npc <- fig_y1 + 0.98 * (fig_y2 - fig_y1)
  upViewport(0)

  panel_label <- if (cohort == "BREEDING") "a" else "b"
  grid.text(panel_label,
            x = unit(label_x_npc, "npc"), y = unit(label_y_npc, "npc"),
            just = c("left", "top"),
            gp = gpar(cex = label_cex, fontface = "bold"))

  invisible(out_path)
}

############################################################
# 6) BUILD CHROMOSOME LAYOUT
############################################################

# Guard: skip main loop when sourced by another script (functions defined above remain available)
if (!exists("SOURCED_14AB") || !isTRUE(SOURCED_14AB)) {

sep_line()
log_msg("Step 14ab -- Circular Manhattan plots")
log_msg("Canvas: ", OUT_W_CM, " cm x ", OUT_H_CM, " cm @ ", OUT_RES, " dpi")
log_msg("Building chromosome layout from GDS files...")

contig_len_dt <- build_contig_lengths(unlist(GDS_FILES))
log_msg("Contigs found: ", nrow(contig_len_dt))

chrom_map_df <- read_chrom_map(CHROM_MAP_FILE)
if (is.null(chrom_map_df)) {
  log_msg("No chrom.list -- building default map from GDS (sorted by length)")
  chrom_map_df <- build_default_chrom_map(contig_len_dt)
  fwrite(as.data.table(chrom_map_df),
         file.path(OUT_MAN, "chrom_map_auto.tsv"), sep = "\t")
}

layout_info    <- prepare_layout(contig_len_dt, chrom_map_df)
contig_info_dt <- layout_info$contig_info
sector_lens    <- layout_info$sector_lens
log_msg("Sector order: ", paste(names(sector_lens), collapse = ", "))

############################################################
# 7) DRAW ALL PLOTS -- both p-value and FDR axis modes
############################################################

for (fdr_mode in c(FALSE, TRUE)) {
  FDR_AXIS <- fdr_mode
  suffix   <- if (fdr_mode) "_fdr" else ""
  sep_line()
  log_msg("=== Axis mode: ", if (fdr_mode) "BH-FDR" else "raw p-value",
          " | suffix: ", suffix, " ===")

  for (tool in TOOLS) {
    for (cohort in COHORTS) {
      sep_line()
      log_msg("Processing: ", tool, " | ", cohort)

      dat <- load_panel_data(tool, cohort)

      out_base <- file.path(OUT_MAN,
        paste0("manhattan_circular_", tolower(tool), "_", tolower(cohort), suffix))

      for (fmt in c("tiff", "pdf", "svg", "eps")) {
        out_path <- paste0(out_base, ".", fmt)
        log_msg("  Drawing (", fmt, "): ", basename(out_path))
        tryCatch(
          draw_circular_manhattan(dat, tool, cohort, out_path),
          error = function(e) log_msg("  ERROR (", fmt, "): ", conditionMessage(e))
        )
        log_msg("  Saved: ", out_path)
      }
    }
  }
}

############################################################
# 8) LEGENDS
############################################################

ctx_colors_global <- CTX_COLORS[CTX_ORDER]

draw_legend <- function(path, width_cm, height_cm = 2.0,
                        pt_cex = 2.0, text_cex = LEGEND_CEX) {
  tiff(path, width = width_cm, height = height_cm,
       units = "cm", res = OUT_RES, compression = "lzw")
  par(mar = c(0, 0, 0, 0))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  legend("center", legend = CTX_ORDER, pch = 16, col = ctx_colors_global,
         pt.cex = pt_cex, cex = text_cex, horiz = TRUE, bty = "n",
         x.intersp = 1.2)
  dev.off()
}

leg_36 <- file.path(OUT_MAN, "legend_circular.tiff")
draw_legend(leg_36, width_cm = OUT_W_CM * 2)
log_msg("Legend (36 cm) saved: ", leg_36)

leg_72 <- file.path(OUT_MAN, "legend_circular_wide.tiff")
draw_legend(leg_72, width_cm = OUT_W_CM * 4)
log_msg("Legend (72 cm) saved: ", leg_72)

writeLines(capture.output(sessionInfo()),
           file.path(LOG_DIR, "step14ab_sessionInfo.txt"))

sep_line()
log_msg("Individual plots done. Building combined panels...")

############################################################
# 9) DRAW COMBINED PANELS (BREEDING | NATURAL + legend)
############################################################

# Logfile already opened above (step14ab.log)
if (HAS_MAGICK) log_msg("ImageMagick: ", MAGICK_BIN) else
  log_msg("WARNING: ImageMagick not found — TIFF/PNG/EPS combined panels skipped")

ctx_colors_combined <- ctx_colors_global   # alias; defined in section 8

draw_combined_panel <- function(tool, out_base, fdr_mode) {

  FDR_AXIS <<- fdr_mode
  axis_tag  <- if (fdr_mode) "_fdr" else ""

  log_msg("  Loading data: ", tool, " BREEDING + NATURAL [FDR_AXIS=", fdr_mode, "]")
  dat_b <- load_panel_data(tool, "BREEDING")
  dat_n <- load_panel_data(tool, "NATURAL")

  draw_both <- function(dev_open_fn) {
    dev_open_fn()
    on.exit(dev.off(), add = TRUE)
    layout(matrix(c(1L, 2L, 3L, 3L), nrow = 2L, byrow = TRUE),
           widths  = c(1, 1),
           heights = c(OUT_W_CM, LEG_H_CM))
    par(cex = 1)   # layout() shrinks text by 0.83; keep the sizes of the single-cohort plots
    draw_circular_manhattan(dat_b, tool, "BREEDING",
                            open_device = FALSE, label_cex = PANEL_LETTER_CEX)
    draw_circular_manhattan(dat_n, tool, "NATURAL",
                            open_device = FALSE, label_cex = PANEL_LETTER_CEX)
    par(mar = c(0, 0, 0, 0))
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1))
    legend("center", legend = CTX_ORDER, pch = 16, col = ctx_colors_combined,
           pt.cex = 2.0, cex = LEGEND_CEX, horiz = TRUE, bty = "n", x.intersp = 1.2)
  }

  # PDF (cairo_pdf — vectorized)
  pdf_out <- paste0(out_base, ".pdf")
  log_msg("  Drawing PDF: ", basename(pdf_out))
  tryCatch(
    draw_both(function()
      cairo_pdf(pdf_out, width = COMB_W_CM / 2.54, height = COMB_H_CM / 2.54)),
    error = function(e) log_msg("  ERROR (pdf): ", conditionMessage(e))
  )
  if (check_nonempty(pdf_out)) log_msg("  PDF OK (", file.info(pdf_out)$size, " bytes)") else
    log_msg("  WARNING: PDF may be empty")

  # SVG (vectorized)
  svg_out <- paste0(out_base, ".svg")
  log_msg("  Drawing SVG: ", basename(svg_out))
  tryCatch(
    draw_both(function()
      svg(svg_out, width = COMB_W_CM / 2.54, height = COMB_H_CM / 2.54)),
    error = function(e) log_msg("  ERROR (svg): ", conditionMessage(e))
  )
  if (check_nonempty(svg_out)) log_msg("  SVG OK (", file.info(svg_out)$size, " bytes)") else
    log_msg("  WARNING: SVG may be empty")

  # TIFF (ImageMagick — stitch from individual per-cohort TIFFs)
  if (HAS_MAGICK) {
    tool_l  <- tolower(tool)
    img_b   <- file.path(OUT_MAN,
                 sprintf("manhattan_circular_%s_breeding%s.tiff", tool_l, axis_tag))
    img_n   <- file.path(OUT_MAN,
                 sprintf("manhattan_circular_%s_natural%s.tiff",  tool_l, axis_tag))
    leg_tif <- file.path(OUT_MAN, "legend_circular.tiff")

    if (all(file.exists(img_b, img_n, leg_tif))) {
      TMP      <- file.path(PANEL_DIR, ".tmp14ab")
      dir.create(TMP, recursive = TRUE, showWarnings = FALSE)
      pair_tmp <- file.path(TMP, sprintf("%s%s_circs.tiff", tool_l, axis_tag))
      run_magick(paste(shQuote(img_b), shQuote(img_n),
                       "+append -compress lzw", shQuote(pair_tmp)))
      tiff_out <- paste0(out_base, ".tiff")
      run_magick(paste(shQuote(pair_tmp), shQuote(leg_tif),
                       "-append -compress lzw", shQuote(tiff_out)))
      log_msg("  TIFF assembled: ", file.info(tiff_out)$size, " bytes")
      png_out <- paste0(out_base, ".png")
      run_magick(paste("-density 150", shQuote(tiff_out), shQuote(png_out)))
      log_msg("  PNG: ", basename(png_out))
      eps_out <- paste0(out_base, ".eps")
      run_magick(paste("-density 300", shQuote(tiff_out), "-compress zip", shQuote(paste0("eps3:", eps_out))))
      log_msg("  EPS: ", basename(eps_out))
      unlink(TMP, recursive = TRUE)
    } else {
      log_msg("  WARNING: individual TIFFs missing — TIFF/PNG/EPS skipped")
    }
  }

  invisible(out_base)
}

figure_specs <- list(
  list(tool = "GENESIS",    fdr = FALSE,
       stem = "circular_panel_genesis"),
  list(tool = "MATRIXEQTL", fdr = FALSE,
       stem = "circular_panel_matrixeqtl"),
  list(tool = "GENESIS",    fdr = TRUE,
       stem = "circular_panel_genesis_FDR"),
  list(tool = "MATRIXEQTL", fdr = TRUE,
       stem = "circular_panel_matrixeqtl_FDR")
)

for (spec in figure_specs) {
  sep_line()
  log_msg("Panel: ", spec$stem)
  draw_combined_panel(spec$tool,
                      file.path(PANEL_DIR, spec$stem),
                      spec$fdr)
}

sep_line()
log_msg("Step 14ab finished")
log_msg("  Individual plots : ", OUT_MAN)
log_msg("  Combined panels  : ", PANEL_DIR)
sep_line()

} # end SOURCED_14AB guard
```

---

### `15ab.tgc.joint.venn.overlap.R`

```r
#!/usr/bin/env Rscript
############################################################
# TreeGeneClimate (TGC) — JOINT ECS + TMS
# Step 15ab: Venn-diagram overlap analysis of significant meQTL markers
#
# "Robust" = SNP-site pair significant at p_FDR < 1e-10 in BOTH tools.
#
# OUTPUTS
# -------
# Tool comparison (6 Venns, labelled a)–f)):
#   tools/venn_tools_<cohort>_<ctx>.{tiff,pdf}
#   Order: BREEDING/CpG (a), BREEDING/CHG (b), BREEDING/CHH (c),
#          NATURAL/CpG  (d), NATURAL/CHG  (e), NATURAL/CHH  (f)
#   Sets: GENESIS vs MatrixEQTL — unit: SNPs
#
# Cohort comparison by context (3 Venns, c)–e)):
#   comparisons/venn_cohorts_<ctx>.{tiff,pdf}
#   Sets: BREEDING vs NATURAL — unit: SNP positions (robust)
#
# Context comparison by cohort (2 Venns, a)–b)):
#   comparisons/venn_contexts_<cohort>.{tiff,pdf}
#   Sets: CpG / CHG / CHH — unit: SNPs (robust)
#
# Tables (one per cohort):
#   tables/robust_markers_<cohort>.tsv
#
# Counts shown in the Venn diagrams:
#   panels/venn_panels_values.tsv
#
# USAGE
#   Rscript --vanilla 15ab.R   (or source() in RStudio)
#
# INPUTS
#   RESULTS/JOINT/MEQTL/COMBINED/sig_fdr1e10_<tool>_<cohort>.tsv  (Step 13ab)
#   RESULTS/JOINT/MEQTL/INPUTS/<COHORT>/snp_variant_annot.tsv     (Step 12ab0)
#
# R packages: data.table, ggvenn, ggplot2, patchwork
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggvenn)
  library(ggplot2)
  library(patchwork)
})

options(stringsAsFactors = FALSE)

############################################################
# 1) SETTINGS
############################################################

COHORTS  <- c("BREEDING", "NATURAL")
CONTEXTS <- c("CpG", "CHG", "CHH")
TOOLS    <- c("GENESIS", "MATRIXEQTL")

# Per-plot fill colours — each plot within a group gets a distinct colour pair/triple.
# Names match original set identifiers (before display-renaming).

# All three Venn types use SNPs as the unit:
#   Tool comparison   : GENESIS SNPs vs MatrixEQTL SNPs  → overlap = robust SNPs
#   Cohort comparison : BREEDING robust SNPs vs NATURAL robust SNPs (per context)
#   Context comparison: CpG vs CHG vs CHH robust SNPs (per cohort)
# "Robust SNPs" = SNPs significant in both tools for the same cohort × context
# (SNP-level intersection; independent of which methylation site is associated).
# The separate robust MARKER TABLE (section 9) uses the stricter pair-level criterion
# (same SNP–site pair in both tools) as required for downstream annotation.

# Methylation-context colours, shared by all figures (magenta, blue, turquoise;
# no red-green pair), and a light tint of each for the first set of a two-set Venn
CTX_COLORS <- c(CpG = "#D81B60", CHG = "#1F4E9C", CHH = "#17BECF")
CTX_TINTS  <- c(CpG = "#EC8DB0", CHG = "#8FA7CE", CHH = "#8BDFE7")

# Tool comparison (6 plots, GENESIS vs MatrixEQTL): tint / full colour of the context
TOOL_FILL_COLORS <- setNames(
  lapply(rep(CONTEXTS, times = 2), function(ctx) c(GENESIS = CTX_TINTS[[ctx]], MATRIXEQTL = CTX_COLORS[[ctx]])),
  paste0(rep(COHORTS, each = 3), "_", CONTEXTS))

# Cohort comparison (3 plots, BREEDING vs NATURAL): tint / full colour of the context
COHORT_FILL_COLORS <- setNames(
  lapply(CONTEXTS, function(ctx) c(BREEDING = CTX_TINTS[[ctx]], NATURAL = CTX_COLORS[[ctx]])), CONTEXTS)

# Context comparison (2 plots, CpG/CHG/CHH): the context colours
CTX_FILL_COLORS <- list(BREEDING = CTX_COLORS, NATURAL = CTX_COLORS)

# Tool comparison: BREEDING first (a–c), NATURAL second (d–f)
TOOL_LABELS <- setNames(
  paste0(letters[1:6], ")"),
  c("BREEDING_CpG", "BREEDING_CHG", "BREEDING_CHH",
    "NATURAL_CpG",  "NATURAL_CHG",  "NATURAL_CHH")
)
CTX_PANEL_LABELS    <- setNames(paste0(letters[3:5], ")"), CONTEXTS)
COHORT_PANEL_LABELS <- setNames(paste0(letters[1:2], ")"), COHORTS)

OUT_W <- 14; OUT_H <- 14; OUT_RES <- 300
OUT_W_2SET <- 20; OUT_H_2SET <- 12   # wider canvas for 2-set horizontal Venns

############################################################
# 2) PATHS
############################################################

# === USER CONFIGURATION ===
# Set TGC_PROJECT_ROOT as an environment variable, or edit the fallback path below
PROJECT_ROOT <- Sys.getenv("TGC_PROJECT_ROOT", unset = "/path/to/your/project")
# ===========================

SIG_DIR  <- file.path(PROJECT_ROOT, "RESULTS", "JOINT", "MEQTL", "COMBINED")
OUT_ROOT <- file.path(PROJECT_ROOT, "RESULTS", "JOINT", "MEQTL", "OVERLAP")
TOOL_DIR      <- file.path(OUT_ROOT, "tools")
CMP_DIR       <- file.path(OUT_ROOT, "comparisons")
TAB_DIR       <- file.path(OUT_ROOT, "tables")
PANEL_DIR     <- file.path(OUT_ROOT, "panels")
LOG_DIR       <- file.path(PROJECT_ROOT, "RESULTS", "JOINT", "MEQTL", "LOGS")

for (d in c(TOOL_DIR, CMP_DIR, TAB_DIR, PANEL_DIR, LOG_DIR))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

LOGFILE <- file.path(LOG_DIR, "step15ab.log")
if (file.exists(LOGFILE)) file.remove(LOGFILE)

############################################################
# 3) HELPERS
############################################################

log_msg <- function(...) {
  txt <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", paste0(..., collapse = ""))
  cat(txt, "\n")
  cat(txt, "\n", file = LOGFILE, append = TRUE)
}

save_plot <- function(gg, base_path,
                      width_cm = OUT_W, height_cm = OUT_H) {
  if (is.null(gg)) return(invisible(NULL))
  w_in <- width_cm  / 2.54
  h_in <- height_cm / 2.54
  for (fmt in c("tiff", "pdf", "eps", "png")) {
    out <- paste0(base_path, ".", fmt)
    tryCatch({
      if (fmt == "eps") {
        grDevices::postscript(out, width = w_in, height = h_in,
                              horizontal = FALSE, paper = "special", onefile = FALSE)
        print(gg)
        grDevices::dev.off()
      } else {
        dpi <- if (fmt %in% c("tiff", "png")) OUT_RES else 150
        ggsave(out, plot = gg, width = width_cm, height = height_cm, units = "cm",
               device = fmt, dpi = dpi)
      }
    }, error = function(e) log_msg("  ERROR (", fmt, "): ", conditionMessage(e)))
    log_msg("  Saved: ", out)
  }
}

# Counts of the Venn regions; a region is the combination of sets
# that contain an element and no other set
venn_counts <- list()
venn_regions <- function(set_list, figure, panel, cohort, context) {
  nms <- names(set_list)
  els <- unique(unlist(set_list))
  memb <- vapply(set_list, function(s) els %in% s, logical(length(els)))
  if (!is.matrix(memb)) memb <- matrix(memb, ncol = length(nms), dimnames = list(NULL, nms))
  region <- apply(memb, 1, function(r) paste(nms[r], collapse = " & "))
  tab <- table(region)
  data.table(figure = figure, panel = sub("[)]$", "", panel), cohort = cohort, context = context,
             sets = paste(nms, collapse = " | "), region = names(tab), n = as.integer(tab))
}

# Human-readable labels for set names
SET_DISPLAY <- c(
  BREEDING    = "Breeding cohort",
  NATURAL     = "Natural cohort",
  GENESIS    = "GENESIS",
  MATRIXEQTL = "MatrixEQTL"
)

# Draw a Venn diagram using ggvenn (ggplot2-native, true circles).
# Returns a ggplot object.
# set_list   : named list of character vectors (original names used for logic)
# fill_colors: named character vector of fill colours, keyed by original names
# panel_label: panel letter string, e.g. "a)"
draw_venn <- function(set_list, fill_colors, panel_label = NULL, sname_size = 7, sd = NULL) {

  set_list <- set_list[lengths(set_list) > 0]
  if (!is.null(sd) && length(set_list) >= 2)
    venn_counts[[length(venn_counts) + 1L]] <<- venn_regions(set_list, sd$figure, panel_label, sd$cohort, sd$context)
  log_msg("    Sets: ",
          paste(names(set_list), lengths(set_list), sep = "=", collapse = " | "))

  n <- length(set_list)

  white_bg <- theme(plot.background  = element_rect(fill = "white", color = NA),
                    panel.background = element_rect(fill = "white", color = NA))

  if (n < 2) {
    log_msg("    WARNING: fewer than 2 non-empty sets — skipping Venn")
    p <- ggplot() +
      annotate("text", x = 0.5, y = 0.5,
               label = paste0("< 2 non-empty sets\n",
                              paste(names(set_list), lengths(set_list),
                                    sep = "=", collapse = ", ")),
               size = 5, color = "grey40") +
      theme_void() + white_bg
    if (!is.null(panel_label))
      p <- p + labs(tag = panel_label) +
        theme(plot.tag = element_text(size = 22, face = "plain"))
    return(p)
  }

  # For 3-set context (CpG/CHG/CHH) Venns: reorder BEFORE renaming so
  # ggvenn places CpG at top-centre, CHG bottom-left, CHH bottom-right.
  if (n == 3 && all(c("CpG", "CHG", "CHH") %in% names(set_list)))
    set_list <- set_list[c("CpG", "CHG", "CHH")]

  orig_names <- names(set_list)

  # Resolve fill colours in set order (fallback grey if not in fill_colors)
  fcolors <- unname(fill_colors[orig_names])
  fcolors[is.na(fcolors)] <- "grey80"

  # Apply human-readable display names
  names(set_list) <- ifelse(orig_names %in% names(SET_DISPLAY),
                            SET_DISPLAY[orig_names], orig_names)

  p <- ggvenn(set_list,
              fill_color      = fcolors,
              fill_alpha      = 0.50,
              stroke_color    = "black",
              stroke_size     = 1.2,
              set_name_size   = sname_size,
              text_size       = 6,    # ~50% larger than default 4
              show_percentage = FALSE) +
    theme(plot.margin = margin(15, 15, 15, 15)) +
    white_bg

  if (!is.null(panel_label)) {
    lbl_up <- sub("[)]$", "", panel_label)
    p <- p + labs(tag = lbl_up) +
      theme(plot.tag = element_text(size = 22, face = "bold", hjust = 0))
  }

  p
}

############################################################
# 3b) OUTPUT DIMENSIONS (updated)
############################################################

# Tool comparison (6 Venns, a)–f)):
#   tools/venn_tools_<cohort>_<ctx>.{tiff,pdf}
#   Sets: GENESIS vs MatrixEQTL — unit: SNPs
#
# Cohort comparison by context (3 Venns, c)–e)):
#   comparisons/venn_cohorts_<ctx>.{tiff,pdf}
#   Sets: BREEDING vs NATURAL — unit: SNPs (robust)
#
# Context comparison by cohort (2 Venns, a)–b)):
#   comparisons/venn_contexts_<cohort>.{tiff,pdf}
#   Sets: CpG / CHG / CHH — unit: SNPs (robust)

############################################################
# 4) LOAD SIGNIFICANT DATA
############################################################

log_msg("Step 15ab — Venn overlap analysis")
log_msg("ggvenn version: ", as.character(packageVersion("ggvenn")))
log_msg("Reading FDR < 1e-10 tables from: ", SIG_DIR)

sig_raw <- list()
for (tool in TOOLS) {
  sig_raw[[tool]] <- list()
  for (cohort in COHORTS) {
    fname <- file.path(SIG_DIR,
      sprintf("sig_fdr1e10_%s_%s.tsv", tolower(tool), tolower(cohort)))
    if (!file.exists(fname)) {
      log_msg("  MISSING: ", fname)
      sig_raw[[tool]][[cohort]] <- data.table()
      next
    }
    dt <- fread(fname, showProgress = FALSE)
    if ("snp"  %in% names(dt)) dt[, snp  := as.character(snp)]
    if ("site" %in% names(dt)) dt[, site := as.character(site)]
    sig_raw[[tool]][[cohort]] <- dt
    log_msg("  ", tool, "/", cohort, ": ", nrow(dt), " pairs | cols: ",
            paste(names(dt), collapse = ", "))
    if ("context" %in% names(dt))
      log_msg("    context values: ",
              paste(sort(unique(dt$context)), collapse = ", "))
  }
}

############################################################
# 4b) LOAD SNP POSITION ANNOTATION
#     Positions of SNP identifiers from the cohort variant annotation (Step 12ab0);
#     cross-cohort comparisons use chr:pos, never integer SNP identifiers.
############################################################

ANNOT_ROOT <- file.path(PROJECT_ROOT, "RESULTS", "JOINT", "MEQTL", "INPUTS")

va_map <- list()   # va_map[[cohort]][[ctx]]: named char vector  pos_string[snp_id]
for (cohort in COHORTS) {
  va_map[[cohort]] <- list()
  for (ctx in CONTEXTS) {
    va_path <- file.path(ANNOT_ROOT, cohort, "snp_variant_annot.tsv")
    if (file.exists(va_path)) {
      va <- fread(va_path)
      va_map[[cohort]][[ctx]] <- setNames(paste0(va$chr, ":", va$pos),
                                          as.character(va$snp_id))
      log_msg("  va_map loaded: ", cohort, "/", ctx,
              " (", length(va_map[[cohort]][[ctx]]), " SNPs)")
    } else {
      log_msg("  va_map MISSING: ", va_path)
      va_map[[cohort]][[ctx]] <- character(0)
    }
  }
}

# Helper: convert integer SNP IDs to unique genomic position strings via va_map
ids_to_pos <- function(ids, va_lookup) {
  pos <- va_lookup[as.character(ids)]
  unique(pos[!is.na(pos) & nchar(pos) > 0])
}

############################################################
# 5) BUILD ROBUST SETS  (intersection of both tools, per cohort × context)
############################################################

log_msg("Building robust sets (intersection of both tools)...")

robust         <- list()   # pair-level robust (same snp+site in both tools) — used for tables
robust_snps    <- list()   # SNP-level robust  (same SNP in both tools)     — used for Venns
robust_snp_pos <- list()   # position-key robust ("chr:pos") for cross-cohort comparison

for (cohort in COHORTS) {
  robust[[cohort]]         <- list()
  robust_snps[[cohort]]    <- list()
  robust_snp_pos[[cohort]] <- list()
  for (ctx in CONTEXTS) {
    g5  <- sig_raw[["GENESIS"   ]][[cohort]]
    me5 <- sig_raw[["MATRIXEQTL"]][[cohort]]

    empty <- data.table()
    if (!nrow(g5) || !nrow(me5) ||
        !"context" %in% names(g5) || !"context" %in% names(me5) ||
        !"p_FDR"   %in% names(g5) || !"p_FDR"   %in% names(me5)) {
      robust[[cohort]][[ctx]]         <- empty
      robust_snps[[cohort]][[ctx]]    <- character(0)
      robust_snp_pos[[cohort]][[ctx]] <- character(0)
      next
    }

    g5_sub  <- g5 [context == ctx]
    me5_sub <- me5[context == ctx]

    if (!nrow(g5_sub) || !nrow(me5_sub)) {
      log_msg("  ", cohort, "/", ctx, ": one or both tools have 0 pairs")
      robust[[cohort]][[ctx]]         <- empty
      robust_snps[[cohort]][[ctx]]    <- character(0)
      robust_snp_pos[[cohort]][[ctx]] <- character(0)
      next
    }

    g5_keys  <- g5_sub [, .(snp, site, FDR_GENESIS    = p_FDR)]
    me5_keys <- me5_sub[, .(snp, site, FDR_MATRIXEQTL = p_FDR)]

    both <- merge(g5_keys, me5_keys, by = c("snp", "site"))
    both[, context := ctx]

    # Attach position columns from GENESIS if available
    pos_cols <- intersect(c("snp_chr", "snp_pos", "site_chr", "site_pos"),
                          names(g5_sub))
    if (length(pos_cols)) {
      g5_pos <- g5_sub[, c("snp", "site", pos_cols), with = FALSE]
      both   <- merge(both, g5_pos, by = c("snp", "site"), all.x = TRUE)
    }

    robust[[cohort]][[ctx]] <- both
    log_msg("  ", cohort, "/", ctx, ": ", nrow(both), " robust pairs | ",
            uniqueN(both$site), " sites | ", uniqueN(both$snp), " SNPs (pair-level)")

    # Pair-level positions: true genomic positions of SNPs in robust pairs only.
    # These are used for the Venn diagrams and cross-cohort/context comparisons.
    robust_snp_pos[[cohort]][[ctx]] <- ids_to_pos(unique(both$snp), va_map[[cohort]][[ctx]])

    # Integer-ID robust (kept for pair-level table joins within a single cohort)
    robust_snps[[cohort]][[ctx]] <- unique(both$snp)
    log_msg("  ", cohort, "/", ctx, ": ",
            length(robust_snps[[cohort]][[ctx]]), " unique SNPs (pair-level) | ",
            length(robust_snp_pos[[cohort]][[ctx]]), " unique positions")
  }
}

############################################################
# 6) TOOL COMPARISON — 6 Venns: GENESIS vs MatrixEQTL per cohort × context
#    Panel labels a)–f): BREEDING (a–c), NATURAL (d–f)
#    Unit: unique genomic positions (SNP chr:pos)
############################################################

# Storage for panel assembly (populated in sections 6-8)
tool_plots   <- list()   # 6 tool-comparison Venns  (key = "COHORT_ctx")
cohort_plots <- list()   # 3 cohort-comparison Venns (key = context)
ctx_plots    <- list()   # 2 context-comparison Venns (key = cohort)

log_msg("--- Tool comparison (6 Venns) ---")

get_sites <- function(dt, ctx) {
  if (!nrow(dt) || !"context" %in% names(dt)) return(character(0))
  unique(dt[context == ctx, site])
}

for (cohort in COHORTS) {
  for (ctx in CONTEXTS) {
    key <- paste0(cohort, "_", ctx)
    lbl <- TOOL_LABELS[[key]]
    log_msg("  Panel ", lbl, "  [", cohort, " / ", ctx, "]")

    # Tool-comparison Venns: pair-level SNP counts.
    # The intersection shown = SNPs targeting the SAME methylation site in both tools
    # (robust pair-level), NOT just the SNP-level set intersection.
    # Synthetic integer sets are built so ggvenn computes the correct counts:
    #   set_A (GENESIS)    : |A| = n_g5,  |A ∩ B| = n_rob
    #   set_B (MatrixEQTL) : |B| = n_me5, |A ∩ B| = n_rob
    # robust IDs = 1..n_rob; G5-only = (n_rob+1)..(n_rob+n_g5_only)
    # ME5-only starts after G5-only to guarantee G5-only ∩ ME5-only = ∅
    g5_sub_ctx  <- sig_raw[["GENESIS"   ]][[cohort]]
    me5_sub_ctx <- sig_raw[["MATRIXEQTL"]][[cohort]]
    if ("context" %in% names(g5_sub_ctx))  g5_sub_ctx  <- g5_sub_ctx [context == ctx]
    if ("context" %in% names(me5_sub_ctx)) me5_sub_ctx <- me5_sub_ctx[context == ctx]

    n_g5  <- uniqueN(g5_sub_ctx$snp)
    n_me5 <- uniqueN(me5_sub_ctx$snp)
    n_rob <- uniqueN(robust[[cohort]][[ctx]]$snp)
    n_g5_only  <- max(0L, n_g5  - n_rob)
    n_me5_only <- max(0L, n_me5 - n_rob)

    set_a <- seq_len(n_rob + n_g5_only)
    set_b <- c(seq_len(n_rob), seq_len(n_me5_only) + n_rob + n_g5_only)

    set_list <- list()
    if (n_g5  > 0L) set_list[["GENESIS"]]    <- set_a
    if (n_me5 > 0L) set_list[["MATRIXEQTL"]] <- set_b

    log_msg("  Pair-level counts — G5: ", n_g5, " | ME5: ", n_me5,
            " | Robust: ", n_rob)

    base      <- file.path(TOOL_DIR,
      paste0("venn_tools_", tolower(cohort), "_", tolower(ctx)))

    gg <- draw_venn(set_list, TOOL_FILL_COLORS[[key]], lbl, sname_size = 5.5,
                    sd = list(figure = "tools", cohort = cohort, context = ctx))
    tool_plots[[key]] <- gg
    save_plot(gg, base,      width_cm = 18, height_cm = 12)
  }
}

############################################################
# 7) COHORT COMPARISON — 3 Venns: BREEDING vs NATURAL per context (c)–e))
#    Unit: unique genomic positions (SNP chr:pos)
############################################################

log_msg("--- Cohort comparison by context (3 Venns) ---")

for (ctx in CONTEXTS) {
  lbl <- CTX_PANEL_LABELS[[ctx]]
  log_msg("  Panel ", lbl, "  [", ctx, "]")

  breed_snps   <- robust_snp_pos[["BREEDING"]][[ctx]]
  natural_snps <- robust_snp_pos[["NATURAL" ]][[ctx]]

  set_list <- list()
  if (length(breed_snps))   set_list[["BREEDING"]] <- breed_snps
  if (length(natural_snps)) set_list[["NATURAL"]]  <- natural_snps

  base      <- file.path(CMP_DIR,      paste0("venn_cohorts_", tolower(ctx)))

  gg <- draw_venn(set_list, COHORT_FILL_COLORS[[ctx]], lbl, sname_size = 5.5,
                  sd = list(figure = "cohorts_by_context", cohort = "BREEDING | NATURAL", context = ctx))
  cohort_plots[[ctx]] <- gg
  save_plot(gg, base,      width_cm = 14, height_cm = 10)
}

############################################################
# 8) CONTEXT COMPARISON — 2 Venns: CpG vs CHG vs CHH per cohort (a)–b))
#    Unit: unique genomic positions (SNP chr:pos)
############################################################

log_msg("--- Context comparison by cohort (2 Venns) ---")

for (cohort in COHORTS) {
  lbl <- COHORT_PANEL_LABELS[[cohort]]
  log_msg("  Panel ", lbl, "  [", cohort, "]")

  snp_sets <- list()
  for (ctx in CONTEXTS) {
    snps <- robust_snp_pos[[cohort]][[ctx]]
    if (length(snps)) snp_sets[[ctx]] <- snps
  }

  base      <- file.path(CMP_DIR,      paste0("venn_contexts_", tolower(cohort)))

  gg <- draw_venn(snp_sets, CTX_FILL_COLORS[[cohort]], lbl,
                  sd = list(figure = "contexts_by_cohort", cohort = cohort, context = "CpG | CHG | CHH"))
  ctx_plots[[cohort]] <- gg
  save_plot(gg, base,      width_cm = 21, height_cm = 14)
}

############################################################
# 8b) PANEL ASSEMBLY
#   Comparison panel: row1 = ctx_plots (a, b — context overlap per cohort)
#                     row2 = cohort_plots (c-e — cross-cohort position overlap per context)
#   Tool panel      : 2×3 grid, GENESIS vs MatrixEQTL per cohort×context
############################################################

log_msg("--- Panel assembly ---")

row1_valid <- Filter(Negate(is.null), ctx_plots)
row2_valid <- Filter(Negate(is.null), cohort_plots)

if (length(row1_valid) >= 1 && length(row2_valid) >= 1) {
  row1 <- Reduce(`+`, row1_valid) + plot_layout(ncol = length(row1_valid))
  row2 <- Reduce(`+`, row2_valid) + plot_layout(ncol = length(row2_valid))
  panel_plot <- row1 / row2 + plot_layout(heights = c(1, 1))
  save_plot(panel_plot,
            file.path(PANEL_DIR, "venn_contexts_cohorts_panel"),
            width_cm = 42, height_cm = 28)
  log_msg("  comparison panel saved → ", PANEL_DIR)
} else {
  log_msg("  WARNING: insufficient Venn plots for comparison panel")
}

row1s <- Filter(Negate(is.null), tool_plots[paste0("BREEDING_", CONTEXTS)])
row2s <- Filter(Negate(is.null), tool_plots[paste0("NATURAL_",  CONTEXTS)])

if (length(row1s) >= 1 && length(row2s) >= 1) {
  row1 <- Reduce(`+`, row1s) + plot_layout(ncol = 3)
  row2 <- Reduce(`+`, row2s) + plot_layout(ncol = 3)
  fs4  <- row1 / row2 + plot_layout(heights = c(1, 1))
  save_plot(fs4,
            file.path(PANEL_DIR, "venn_tools_panel"),
            width_cm = 54, height_cm = 24)
  log_msg("  tool panel saved → ", PANEL_DIR)
} else {
  log_msg("  WARNING: insufficient Venn plots for tool panel")
}

############################################################
# 9) TABLES — one per cohort, robust pairs, all contexts combined
############################################################

log_msg("--- Saving robust marker tables ---")

for (cohort in COHORTS) {
  rows <- lapply(CONTEXTS, function(ctx) robust[[cohort]][[ctx]])
  rows <- rows[sapply(rows, nrow) > 0]

  if (!length(rows)) {
    log_msg("  No robust pairs for ", cohort, " — no table written")
    next
  }

  tab <- rbindlist(rows, fill = TRUE)
  col_order <- intersect(
    c("snp", "snp_chr", "snp_pos", "site", "site_chr", "site_pos",
      "context", "FDR_GENESIS", "FDR_MATRIXEQTL"),
    names(tab))
  setcolorder(tab, col_order)
  setorder(tab, context, snp_chr, snp_pos)

  out <- file.path(TAB_DIR, paste0("robust_markers_", tolower(cohort), ".tsv"))
  fwrite(tab, out, sep = "\t")
  log_msg("  ", cohort, ": ", nrow(tab), " pairs | ",
          uniqueN(tab$site), " sites | ", uniqueN(tab$snp),
          " SNPs → ", basename(out))
}

############################################################
# 9b) ROBUST CONTEXT SUMMARY — per cohort × context counts
############################################################

log_msg("--- Saving robust context summary ---")

rob_ctx_rows <- list()
for (cohort in COHORTS) {
  for (ctx in CONTEXTS) {
    rob <- robust[[cohort]][[ctx]]
    n_pairs      <- if (is.data.table(rob) && nrow(rob) > 0) nrow(rob)          else 0L
    n_sites      <- if (is.data.table(rob) && nrow(rob) > 0) uniqueN(rob$site)  else 0L
    n_snps_pair  <- if (is.data.table(rob) && nrow(rob) > 0) uniqueN(rob$snp)   else 0L
    n_snps_snp   <- length(robust_snps[[cohort]][[ctx]])
    n_pos        <- length(robust_snp_pos[[cohort]][[ctx]])
    rob_ctx_rows[[paste(cohort, ctx)]] <- data.table(
      cohort                       = cohort,
      context                      = ctx,
      robust_pairs                 = n_pairs,
      robust_unique_sites          = n_sites,
      robust_unique_snps_pairlevel = n_snps_pair,
      robust_unique_snps_snplevel  = n_snps_snp,
      robust_unique_positions      = n_pos
    )
  }
}
robust_ctx_summary <- rbindlist(rob_ctx_rows)
fwrite(robust_ctx_summary,
       file.path(TAB_DIR, "robust_context_summary.tsv"), sep = "\t")
log_msg("Robust context summary saved: ", file.path(TAB_DIR, "robust_context_summary.tsv"))
print(robust_ctx_summary)

# Counts of every Venn region
fwrite(rbindlist(venn_counts), file.path(PANEL_DIR, "venn_panels_values.tsv"), sep = "\t")

############################################################
# 10) SESSION INFO
############################################################

writeLines(capture.output(sessionInfo()),
           file.path(LOG_DIR, "step15ab_sessionInfo.txt"))

log_msg("Step 15ab finished — outputs in: ", OUT_ROOT)
```

---

### `16ab.tgc.joint.meth.heritability.R`

```r
#!/usr/bin/env Rscript
############################################################
# TreeGeneClimate (TGC) — JOINT ECS + TMS
# Step 16ab: Per-site heritability of DNA methylation and whole-methylome
#            Mantel test (GRM kinship vs methylation similarity)
#
# For every methylation site entering cis-meQTL mapping (Step 12ab0):
#   - h2_SNP: REML via eigendecomposition of the GRM (EMMA)
#             on M-values; h2 = Vg / (Vg + Ve)
#   - h2_Ped (breeding only): 2 x ICC from lmer(M ~ 1 + (1 | family))
# Whole-methylome Mantel test per cohort x context: Pearson correlation between
#   pairwise GRM kinship and pairwise methylation similarity 1 - d_ij / max(d)
#   (Euclidean distance over all sites entering cis-meQTL mapping, computed
#   for each pair of trees over the sites observed in both, without imputing
#   missing values; trees observed at < 20% of the sites left out), 9,999
#   symmetric permutations of the kinship matrix.
# GRM: GCTA-estimator GRM from called genotypes (Step 8a), the same used for
#   mapping.
#
# MODES
#   chunk   : TGC_MODE=chunk TGC_COHORT=BREEDING TGC_CONTEXT=CHH TGC_CHUNK=7
#             (sites (TGC_CHUNK-1)*CHUNK_SIZE+1 ... TGC_CHUNK*CHUNK_SIZE)
#   collect : TGC_MODE=collect  (merges chunks, summary, heritability panel, Mantel test)
#
# INPUTS
#   RESULTS/JOINT/MEQTL/INPUTS/<COHORT>/<CTX>/methylation_mvalues_matrix.rds
#   RESULTS/JOINT/MEQTL/COMBINED/robust_pairs_<cohort>.tsv     (Step 13ab)
#   RESULTS/ECS/RANALYSIS/RDATA/<cohort>_grm_gcta.rds           (Step 8a)
#   DATA/METADATA/breeding_sample2family.txt
#
# OUTPUTS (RESULTS/JOINT/HERITABILITY/)
#   CHUNKS/<COHORT>_<CTX>_<k>.tsv.gz
#   h2_all_sites.tsv.gz, h2_summary.tsv, h2_meqtl_vs_background.tsv
#   mantel_grm_methylation.tsv
#   FIGURES/heritability_panel.{tiff,pdf,png,eps}
#   FIGURES/heritability_panel_h2_values.tsv      h2_SNP of every site, robust cis-meQTL or other (panel a)
#   FIGURES/heritability_panel_h2_ped_values.tsv  h2_SNP and h2_Ped of the breeding sites shown (panel b)
#
# R packages: data.table, lme4, ggplot2, patchwork
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})
options(stringsAsFactors = FALSE)

# === USER CONFIGURATION ===
PROJECT_ROOT <- Sys.getenv("TGC_PROJECT_ROOT", unset = "/path/to/your/project")
CHUNK_SIZE <- 2000L
N_PERM     <- 9999L
MIN_OBS_MANTEL <- 0.20   # Mantel test: trees observed at < 20% of sites are left out
# ===========================

MODE <- Sys.getenv("TGC_MODE", unset = "collect")
COHORTS  <- c("BREEDING", "NATURAL")
CONTEXTS <- c("CpG", "CHG", "CHH")
IN_ROOT  <- file.path(PROJECT_ROOT, "RESULTS/JOINT/MEQTL/INPUTS")
ROB_DIR  <- file.path(PROJECT_ROOT, "RESULTS/JOINT/MEQTL/COMBINED")
RDATA    <- file.path(PROJECT_ROOT, "RESULTS/ECS/RANALYSIS/RDATA")
H2_ROOT  <- file.path(PROJECT_ROOT, "RESULTS/JOINT/HERITABILITY")
FIG_DIR  <- file.path(H2_ROOT, "FIGURES")
for (d in c(file.path(H2_ROOT, "CHUNKS"), FIG_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
msg <- function(...) { cat(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n")); flush(stdout()) }

load_inputs <- function(cohort, ctx) {
  M <- as.matrix(readRDS(file.path(IN_ROOT, cohort, ctx, "methylation_mvalues_matrix.rds")))
  K <- as.matrix(readRDS(file.path(RDATA, sprintf("%s_grm_gcta.rds", tolower(cohort)))))
  common <- intersect(rownames(M), rownames(K))
  list(M = M[common, , drop = FALSE], K = K[common, common, drop = FALSE])
}

# Per-site REML via eigendecomposition (EMMA)
emma_reml <- function(y, U, D_pos) {
  y[is.na(y)] <- mean(y, na.rm = TRUE)
  y  <- y - mean(y)
  Uy <- drop(crossprod(U, y))
  n  <- length(D_pos)
  obj <- function(ld) { d <- D_pos + exp(ld); vg <- sum(Uy^2 / d) / n
                        if (vg <= 0) return(1e15); n * log(vg) + sum(log(d)) }
  opt <- tryCatch(optimize(obj, c(-20, 20), tol = 1e-8), error = function(e) NULL)
  if (is.null(opt)) return(c(h2 = NA_real_, Vg = NA_real_, Ve = NA_real_))
  delta <- exp(opt$minimum); d <- D_pos + delta
  Vg <- sum(Uy^2 / d) / n; Ve <- Vg * delta
  c(h2 = if (Vg + Ve > 1e-15) Vg / (Vg + Ve) else NA_real_, Vg = Vg, Ve = Ve)
}

# Pedigree h2 = 2 x ICC from a full-sib family random effect
pedigree_h2 <- function(y, fam) {
  y[is.na(y)] <- mean(y, na.rm = TRUE)
  fit <- tryCatch(suppressWarnings(suppressMessages(
    lme4::lmer(y ~ 1 + (1 | fam), data = data.frame(y = y, fam = fam), REML = TRUE,
               control = lme4::lmerControl(optimizer = "nloptwrap", optCtrl = list(maxfun = 2e5),
                                           check.conv.singular = "ignore")))),
    error = function(e) NULL)
  if (is.null(fit)) return(c(h2_ped = NA_real_, ICC = NA_real_))
  vc <- as.data.frame(lme4::VarCorr(fit))
  Vf <- vc$vcov[vc$grp == "fam"]; if (!length(Vf)) Vf <- 0
  Ve <- vc$vcov[vc$grp == "Residual"]
  ICC <- if (Vf + Ve > 1e-15) Vf / (Vf + Ve) else NA_real_
  c(h2_ped = 2 * ICC, ICC = ICC)
}

############################################################
# CHUNK MODE
############################################################
if (MODE == "chunk") {
  cohort <- Sys.getenv("TGC_COHORT"); ctx <- Sys.getenv("TGC_CONTEXT")
  k <- as.integer(Sys.getenv("TGC_CHUNK"))
  inp <- load_inputs(cohort, ctx)
  idx <- seq.int((k - 1L) * CHUNK_SIZE + 1L, min(k * CHUNK_SIZE, ncol(inp$M)))
  if (idx[1] > ncol(inp$M)) { msg("chunk ", k, " beyond ", ncol(inp$M), " sites"); quit(save = "no") }
  eig <- eigen(inp$K, symmetric = TRUE); pos <- eig$values > 1e-10
  U <- eig$vectors[, pos, drop = FALSE]; D <- eig$values[pos]
  fam <- NULL
  if (cohort == "BREEDING") {
    fm <- fread(file.path(PROJECT_ROOT, "DATA/METADATA/breeding_sample2family.txt"),
                header = FALSE, col.names = c("sample", "family"))
    fam <- factor(fm$family[match(rownames(inp$M), fm$sample)])
  }
  res <- rbindlist(lapply(idx, function(i) {
    y <- inp$M[, i]
    r <- as.list(emma_reml(y, U, D))
    if (!is.null(fam)) r <- c(r, as.list(pedigree_h2(y, fam)))
    as.data.table(r)
  }), fill = TRUE)
  res[, `:=`(site = colnames(inp$M)[idx], cohort = cohort, context = ctx)]
  fwrite(res, file.path(H2_ROOT, "CHUNKS", sprintf("%s_%s_%04d.tsv.gz", cohort, ctx, k)), sep = "\t")
  msg(cohort, "/", ctx, " chunk ", k, ": ", nrow(res), " sites")
  quit(save = "no")
}

############################################################
# COLLECT MODE
############################################################
h2_all <- rbindlist(lapply(list.files(file.path(H2_ROOT, "CHUNKS"), "\\.tsv\\.gz$", full.names = TRUE), fread), fill = TRUE)
for (co in COHORTS) for (cx in CONTEXTS) {
  n_exp <- ncol(readRDS(file.path(IN_ROOT, co, cx, "methylation_mvalues_matrix.rds")))
  n_got <- h2_all[cohort == co & context == cx, .N]
  if (n_got != n_exp) stop(co, "/", cx, ": ", n_got, " of ", n_exp, " sites in chunks")
}
setorder(h2_all, cohort, context, site)

# robust meQTL sites, cohort-specific
rob <- rbindlist(lapply(COHORTS, function(co)
  unique(fread(file.path(ROB_DIR, sprintf("robust_pairs_%s.tsv", tolower(co))), select = c("site", "context"))[, cohort := co])))
h2_all[, is_meqtl := FALSE]
h2_all[rob, on = .(cohort, context, site), is_meqtl := TRUE]
fwrite(h2_all, file.path(H2_ROOT, "h2_all_sites.tsv.gz"), sep = "\t")

summ <- h2_all[, .(n_sites = .N, n_valid = sum(!is.na(h2)), mean_h2 = mean(h2, na.rm = TRUE),
                   median_h2 = median(h2, na.rm = TRUE), pct_h2_gt0.3 = 100 * mean(h2 > 0.3, na.rm = TRUE),
                   mean_h2_ped = if ("h2_ped" %in% names(.SD)) mean(h2_ped, na.rm = TRUE) else NA_real_),
               by = .(cohort, context)]
setorder(summ, cohort, context)
fwrite(summ, file.path(H2_ROOT, "h2_summary.tsv"), sep = "\t")
mvb <- h2_all[, .(n = .N, mean_h2 = mean(h2, na.rm = TRUE), median_h2 = median(h2, na.rm = TRUE)),
              by = .(cohort, context, group = ifelse(is_meqtl, "robust_meqtl", "background"))]
p_w <- h2_all[, {
  w <- tryCatch(wilcox.test(h2[is_meqtl], h2[!is_meqtl]), error = function(e) NULL)
  .(wilcox_W = if (is.null(w)) NA_real_ else unname(w$statistic),
    wilcox_p = if (is.null(w)) NA_real_ else w$p.value)
}, by = .(cohort, context)]
mvb <- merge(mvb, p_w, by = c("cohort", "context"))
fwrite(mvb, file.path(H2_ROOT, "h2_meqtl_vs_background.tsv"), sep = "\t")
print(summ); print(mvb)

# ---- Mantel test: GRM kinship vs methylation similarity ----
set.seed(42)
mantel_rows <- list()
for (cohort in COHORTS) for (ctx in CONTEXTS) {
  inp <- load_inputs(cohort, ctx)
  M <- inp$M
  obs <- rowMeans(!is.na(M))
  keep <- obs >= MIN_OBS_MANTEL           # leave out failed / very low-complexity libraries
  M <- M[keep, , drop = FALSE]; inp$K <- inp$K[keep, keep]
  d <- as.matrix(dist(M))   # missing values excluded pair by pair (dist scales to all sites)
  sim <- 1 - d / max(d)
  lt <- lower.tri(sim)
  sim_vec <- sim[lt]; K <- inp$K
  r_obs <- cor(K[lt], sim_vec)
  n <- nrow(K)
  perm_r <- vapply(seq_len(N_PERM), function(i) { o <- sample.int(n); cor(K[o, o][lt], sim_vec) }, numeric(1))
  mantel_rows[[paste(cohort, ctx)]] <- data.table(cohort = cohort, context = ctx, n_samples = n,
    n_sites = ncol(M), mantel_r = r_obs, p_value = (sum(perm_r >= r_obs) + 1) / (N_PERM + 1), n_perm = N_PERM)
  msg("Mantel ", cohort, "/", ctx, ": r = ", round(r_obs, 3))
}
mantel <- rbindlist(mantel_rows)
fwrite(mantel, file.path(H2_ROOT, "mantel_grm_methylation.tsv"), sep = "\t")
print(mantel)

# ---- Heritability panel ----
save_fig <- function(base, plot, width, height) {
  ggsave(paste0(base, ".tiff"), plot, width = width, height = height, dpi = 300, compression = "lzw")
  ggsave(paste0(base, ".pdf"),  plot, width = width, height = height)
  ggsave(paste0(base, ".png"),  plot, width = width, height = height, dpi = 300)
  ggsave(paste0(base, ".eps"),  plot, width = width, height = height, device = cairo_ps)
}
# Methylation-context colours, shared by all figures (magenta, blue, turquoise; no red-green pair)
ctx_cols   <- c(CpG = "#D81B60", CHG = "#1F4E9C", CHH = "#17BECF")
coh_labels <- c(BREEDING = "Breeding (n = 209)", NATURAL = "Natural (n = 393)")
hp <- h2_all[!is.na(h2)]
hp[, context_f := factor(context, levels = CONTEXTS)]
hp[, cohort_f  := factor(cohort, levels = COHORTS, labels = coh_labels)]
hp[, group := factor(ifelse(is_meqtl, "Robust meQTL sites", "Other sites"), levels = c("Other sites", "Robust meQTL sites"))]

p_a <- ggplot(hp, aes(x = context_f, y = h2, fill = group, colour = group)) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.5, position = position_dodge(0.8), linewidth = 0.3) +
  geom_boxplot(width = 0.06, outlier.shape = NA, fill = NA, position = position_dodge(0.8), linewidth = 0.4, coef = 0) +
  facet_wrap(~cohort_f, nrow = 1) +
  scale_fill_manual(values = c("Other sites" = "grey60", "Robust meQTL sites" = "#E63946"), name = NULL) +
  scale_colour_manual(values = c("Other sites" = "grey40", "Robust meQTL sites" = "#9B1A22"), name = NULL) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(x = "Methylation context", y = expression(italic(h)[SNP]^2)) +
  theme_bw(base_size = 11) +
  theme(strip.background = element_rect(fill = "grey90"), legend.position = "bottom", panel.grid.minor = element_blank())

bp <- h2_all[cohort == "BREEDING" & !is.na(h2) & !is.na(h2_ped)]
bp[, context_f := factor(context, levels = CONTEXTS)]
set.seed(42)
bp <- bp[, .SD[sample(.N, min(.N, 50000L))], by = context_f]
p_b <- ggplot(bp, aes(x = h2, y = h2_ped, colour = context_f)) +
  geom_point(alpha = 0.08, size = 0.4, stroke = 0) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "black", linewidth = 0.5) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.8) +
  facet_wrap(~context_f, nrow = 1) +
  scale_colour_manual(values = ctx_cols, guide = "none") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(x = expression(italic(h)[SNP]^2), y = expression(italic(h)[Ped]^2)) +
  theme_bw(base_size = 11) +
  theme(strip.background = element_rect(fill = "grey90"), panel.grid.minor = element_blank())

panel_plot <- (p_a + labs(tag = "a") + theme(plot.tag = element_text(face = "bold", size = 14))) /
        (p_b + labs(tag = "b") + theme(plot.tag = element_text(face = "bold", size = 14))) +
        plot_layout(heights = c(4.5, 3.5))
save_fig(file.path(FIG_DIR, "heritability_panel"), panel_plot, width = 8, height = 8)

# Values shown in panels a and b
fwrite(hp[, .(panel = "a", cohort, context, site, robust_cis_meqtl = is_meqtl, h2_snp = h2)],
       file.path(FIG_DIR, "heritability_panel_h2_values.tsv"), sep = "\t")
fwrite(bp[, .(panel = "b", cohort, context, site, h2_snp = h2, h2_ped)],
       file.path(FIG_DIR, "heritability_panel_h2_ped_values.tsv"), sep = "\t")
msg("Step 16ab done")
```

---

### `16ab.tgc.joint.meth.heritability.sh`

```bash
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
```

---

### `17ab.tgc.joint.marker.annotation.R`

```r
#!/usr/bin/env Rscript
############################################################
# TreeGeneClimate (TGC) — JOINT ECS + TMS
# Step 17ab: Marker annotation against the reference genome GFF3
#
# Annotates:
#   - TMS KW SVMPs selected in Step 8b (breeding: 50 per context; natural: all)
#     -> tms_svmp_<cohort>_annotated.tsv
#   - Robust cis-meQTL markers: SNP and site positions of the pairs significant at
#     p_FDR < 1e-10 in both GENESIS and MatrixEQTL (robust_pairs_<cohort>.tsv, Step 13ab)
#     -> robust_<cohort>_snps_annotated.tsv, robust_<cohort>_sites_annotated.tsv
#   - Gene-centric summary -> gene_summary_annotated.tsv
#   - Genes by methylation-context combination of their robust SNPs (source genes)
#     or robust sites (target genes), with counts per cohort
#     -> source_genes_by_context.tsv, target_genes_by_context.tsv
#
# Each marker is annotated with:
#   - gene_id, gene_start, gene_end, gene_strand (GFF3 gene feature)
#   - distance_bp (0 = overlaps gene; positive = nearest gene distance)
#   - annotation_class (genic / proximal_intergenic / distal_intergenic, or
#     unassigned for contigs without gene models)
#   - ref, alt, af, dr2 (from imputed VCF; SNP markers only)
# The marker tables have one row per marker, context and gene model (a marker
# overlapping two gene models has one row per gene); columns empty in every row
# and the constant marker_type, source and tool columns are not written.
#
# Requires BEDTools (bedtools closest) and BCFtools in the PATH.
# R packages: data.table
############################################################

suppressPackageStartupMessages({
  library(data.table)
})

options(stringsAsFactors = FALSE)

############################################################
# 1) SETTINGS
############################################################

COHORTS  <- c("BREEDING", "NATURAL")
CONTEXTS <- c("CpG", "CHG", "CHH")
TOOLS    <- c("GENESIS", "MATRIXEQTL")

# Methylation classification thresholds (context-specific).
# CHG/CHH methylation is generally much lower than CpG in conifers.
BETA_HYPO  <- c(CpG = 0.30, CHG = 0.20, CHH = 0.10)
BETA_HYPER <- c(CpG = 0.70, CHG = 0.50, CHH = 0.30)

############################################################
# 2) PATHS
############################################################

# === USER CONFIGURATION ===
# Set TGC_PROJECT_ROOT as an environment variable, or edit the fallback path below
PROJECT_ROOT <- Sys.getenv("TGC_PROJECT_ROOT", unset = "/path/to/your/project")
# ===========================

# Picea abies v2.0 reference annotation
GFF3_FILE <- file.path(PROJECT_ROOT,
  "REFERENCE/Pabies2.0/Picab02_230926_at01_all_sorted.gff3")

# Beagle-imputed VCFs (AF + DR2 available); used for allele-frequency annotation
VCF_FILES <- list(
  BREEDING = file.path(PROJECT_ROOT,
    "RESULTS/ECS/VCF_SPLIT/tgc.ecs.breeding.call.filt.maf05.snvs.poly.imputed.vcf.gz"),
  NATURAL  = file.path(PROJECT_ROOT,
    "RESULTS/ECS/VCF_SPLIT/tgc.ecs.natural.call.filt.maf05.snvs.poly.imputed.vcf.gz")
)

# Robust pair tables produced by Step 13ab
ROBUST_ROOT   <- file.path(PROJECT_ROOT, "RESULTS/JOINT/MEQTL/COMBINED")
SVMP_ROOT     <- file.path(PROJECT_ROOT, "RESULTS/TMS/RANALYSIS/TABLES/heatmap_markers_8B")

# Directory with mean methylation beta files (one per cohort × context).
# Expected filename: mean_beta_<cohort>_<context>.tsv  (columns: site, mean_beta)
# Computed automatically from the M-value matrices of Step 12ab0 if not present.
BETA_ROOT <- file.path(PROJECT_ROOT, "RESULTS/TMS/METHYLATION/mean_betas")

# M-value matrix inputs from Step 12ab0 — used to compute mean betas.
# Path: MEQTL_INPUTDIR/<COHORT>/<CONTEXT>/methylation_mvalues_matrix.rds
# Rows = samples, columns = site IDs (chr:start-end).
MEQTL_INPUTDIR <- file.path(PROJECT_ROOT, "RESULTS/JOINT/MEQTL/INPUTS")

# eggNOG-mapper functional annotation (isoform-level; best mRNA selected per gene).
# Columns used: id (1), eggnog_description (5), eggnog_go (8), eggnog_KEGG_ko (9),
#               interpro_description/pfam (15), interpro_panther_description (20).
FUNC_ANNOT_FILE <- file.path(PROJECT_ROOT,
  "REFERENCE/Pabies2.0/Picab02_230926_at01_all_isoform_annotations_merged_sorted_non_redundant_panthers.tsv")

# List of mRNA IDs classified as TE-derived gene models (to build is_te_gene flag).
TE_IDS_FILE <- file.path(PROJECT_ROOT,
  "REFERENCE/Pabies2.0/TE_IDs_REMOVED.txt")

# Cohort VCFs before imputation — fallback for allele information if a SNP is
# absent from the imputed VCF (not expected: mapping uses imputed SNPs only).
RAW_VCF_FILES <- list(
  BREEDING = file.path(PROJECT_ROOT,
    "RESULTS/ECS/VCF_SPLIT/tgc.ecs.breeding.call.filt.maf05.snvs.poly.vcf.gz"),
  NATURAL  = file.path(PROJECT_ROOT,
    "RESULTS/ECS/VCF_SPLIT/tgc.ecs.natural.call.filt.maf05.snvs.poly.vcf.gz")
)

OUT_ROOT <- file.path(PROJECT_ROOT, "RESULTS/JOINT/ANNOTATION")
TMP_DIR  <- file.path(OUT_ROOT, "tmp")
dir.create(TMP_DIR, recursive = TRUE, showWarnings = FALSE)

############################################################
# 3) HELPERS
############################################################

msg <- function(...) cat(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n"))

# Write a marker table: unique rows, without the columns empty in every row and
# without the constant marker_type, source and tool columns.
write_marker_table <- function(dt, path) {
  out <- unique(dt)
  empty <- names(out)[vapply(out, function(x) all(is.na(x) | as.character(x) == ""), logical(1))]
  out <- out[, setdiff(names(out), c(empty, "marker_type", "source", "tool")), with = FALSE]
  fwrite(out, path, sep = "\t")
  invisible(nrow(out))
}

# Write a sorted 0-based BED file for a set of markers.
# chr_col, pos_col (1-based genomic), id_col are column names in dt.
write_marker_bed <- function(dt, chr_col, pos_col, id_col, path) {
  bed <- data.table(
    chr   = dt[[chr_col]],
    start = as.numeric(dt[[pos_col]]) - 1,  # BED is 0-based half-open
    end   = as.numeric(dt[[pos_col]]),
    id    = dt[[id_col]]
  )
  bed <- bed[!is.na(chr) & !is.na(start)]
  setorder(bed, chr, start)
  fwrite(bed, path, sep = "\t", col.names = FALSE)
  invisible(path)
}

# Run bedtools closest -d against genes BED and return annotated data.table.
# Marker BED has 4 columns: chr, start, end, id.
# Gene BED has 6 columns: chr, start, end, gene_id, ., strand.
run_bedtools_closest <- function(marker_bed, genes_bed) {
  if (!file.exists(marker_bed) || file.info(marker_bed)$size == 0)
    return(data.table())

  out_file <- tempfile(tmpdir = TMP_DIR, fileext = ".closest.tsv")
  cmd <- sprintf("bedtools closest -a %s -b %s -d > %s",
                 marker_bed, genes_bed, out_file)
  ret <- system(cmd)
  if (ret != 0 || !file.exists(out_file) || file.info(out_file)$size == 0) {
    msg("  WARNING: bedtools closest returned no results for ", marker_bed)
    return(data.table())
  }

  # Cols 1-4: marker; cols 5-10: gene (chr,start,end,id,.,strand); col 11: distance
  dt <- fread(out_file, header = FALSE, sep = "\t", fill = TRUE,
              col.names = c("chr","pos_start","pos_end","marker_id",
                            "gene_chr","gene_start","gene_end","gene_id",
                            "gene_score","gene_strand","distance_bp"))
  file.remove(out_file)

  dt[, gene_score := NULL]
  # distance = -1 means no gene on that chromosome at all
  dt[distance_bp < 0, `:=`(gene_id     = NA_character_,
                             gene_chr    = NA_character_,
                             gene_start  = NA_integer_,
                             gene_end    = NA_integer_,
                             gene_strand = NA_character_,
                             distance_bp = NA_integer_)]
  # bedtools returns "." when there is no overlapping feature
  dt[gene_chr == ".", gene_chr := NA_character_]
  dt
}

# Query a VCF for REF, ALT, AF and optionally DR2 at given marker positions.
# has_dr2 = FALSE for raw (non-imputed) VCFs that lack the Beagle DR2 INFO field.
query_vcf <- function(markers_dt, vcf_file, tmp_prefix, has_dr2 = TRUE) {
  empty <- data.table(chr = character(), pos = integer(),
                      ref = character(), alt = character(),
                      af  = numeric(),  dr2 = numeric())
  if (!file.exists(vcf_file)) {
    msg("  WARNING: VCF not found: ", vcf_file)
    return(empty)
  }
  markers_dt <- unique(markers_dt[!is.na(chr) & !is.na(pos)])
  if (nrow(markers_dt) == 0) return(empty)

  # Write 0-based BED regions for bcftools -R
  reg_file <- paste0(tmp_prefix, ".regions.bed")
  reg <- markers_dt[, .(chr, start = pos - 1L, end = pos)]
  setorder(reg, chr, start)
  fwrite(reg, reg_file, sep = "\t", col.names = FALSE)

  out_file <- paste0(tmp_prefix, ".vcf_query.tsv")
  # vcf_fmt is substituted as a value (not a sprintf format), so use single %
  vcf_fmt  <- if (has_dr2)
    "'%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/AF\\t%INFO/DR2\\n'"
  else
    "'%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/AF\\n'"
  vcf_cols <- if (has_dr2)
    c("chr","pos","ref","alt","af","dr2")
  else
    c("chr","pos","ref","alt","af")
  cmd <- sprintf("bcftools query -R %s -f %s %s > %s 2>/dev/null",
    reg_file, vcf_fmt, vcf_file, out_file)
  system(cmd)
  file.remove(reg_file)

  if (!file.exists(out_file) || file.info(out_file)$size == 0) {
    return(empty)
  }
  # na.strings="." converts bcftools missing-value dots to NA.
  # Explicit numeric cast prevents fread from typing all-NA columns as logical,
  # which would coerce numeric AF values to TRUE/FALSE when assigned later.
  res <- fread(out_file, header = FALSE, sep = "\t", fill = TRUE,
               col.names = vcf_cols, na.strings = ".")
  res[, pos := as.integer(pos)]
  if ("af"  %in% names(res)) res[, af  := as.numeric(af)]
  if ("dr2" %in% names(res)) res[, dr2 := as.numeric(dr2)]
  if (!"dr2" %in% names(res)) res[, dr2 := NA_real_]
  res
}

# For SNPs with NA ref/alt or NA af, query the raw (non-imputed) VCF.
# The imputed VCF strips INFO for directly genotyped markers, so AF is only
# available in the raw VCF for those sites. DR2 is never in the raw VCF.
# Returns dt with ref/alt/af filled where previously NA.
fill_from_raw_vcf <- function(dt, cohort, tmp_prefix) {
  raw_vcf <- RAW_VCF_FILES[[cohort]]
  if (is.null(raw_vcf) || !file.exists(raw_vcf)) return(dt)
  if (!all(c("ref","alt","af","chr","pos") %in% names(dt))) return(dt)
  missing <- unique(dt[is.na(ref) | is.na(af), .(chr, pos = as.integer(pos))])
  if (nrow(missing) == 0) return(dt)
  msg("  Raw-VCF fill for ", nrow(missing),
      " SNPs with NA REF/ALT or NA AF (imputed VCF strips INFO for typed markers)...")
  vi <- query_vcf(missing, raw_vcf, tmp_prefix, has_dr2 = FALSE)
  if (nrow(vi) == 0) return(dt)
  dt <- merge(dt,
    vi[, .(chr, pos, ref_raw = ref, alt_raw = alt, af_raw = af)],
    by = c("chr", "pos"), all.x = TRUE)
  dt[is.na(ref), `:=`(ref = ref_raw, alt = alt_raw)]
  dt[is.na(af) & !is.na(af_raw), af := af_raw]
  dt[, c("ref_raw","alt_raw","af_raw") := NULL]
  msg("    After fill — REF known: ", dt[!is.na(ref), .N],
      " | AF known: ", dt[!is.na(af), .N])
  dt
}

# Load eggNOG-mapper functional annotation.
# Returns a data.table keyed by gene_id with columns:
#   eggnog_description, pfam_domain, panther_description, go_terms, kegg_ko.
# One row per gene (best/first mRNA isoform).
load_functional_annotation <- function() {
  empty <- data.table(gene_id = character(), eggnog_description = character(),
                      pfam_domain = character(), panther_description = character(),
                      go_terms = character(), kegg_ko = character())
  if (!file.exists(FUNC_ANNOT_FILE)) {
    msg("  WARNING: functional annotation file not found — eggNOG columns will be NA")
    return(empty)
  }
  msg("Loading eggNOG-mapper functional annotation...")
  # Select only needed columns to keep memory use low (~800 MB TSV)
  fa <- fread(FUNC_ANNOT_FILE, header = TRUE,
              select = c(1L, 5L, 8L, 9L, 15L, 20L),
              col.names = c("mrna_id", "eggnog_description", "go_terms",
                            "kegg_ko", "pfam_domain", "panther_description"))
  # Derive gene_id by stripping .mRNA.<n> suffix
  fa[, gene_id := sub("\\.mRNA\\..*$", "", mrna_id)]
  # Keep first row per gene (sorted by mrna_id, so .mRNA.1 comes first)
  fa <- fa[order(mrna_id)][!duplicated(gene_id)]
  fa[, mrna_id := NULL]
  # Replace literal "NA" strings with proper NA
  for (col in c("eggnog_description","go_terms","kegg_ko","pfam_domain","panther_description"))
    set(fa, which(fa[[col]] == "NA"), col, NA_character_)
  setkey(fa, gene_id)
  msg("  Loaded: ", nrow(fa), " gene functional annotations")
  fa
}

# Load TE-derived gene IDs from TE_IDS_REMOVED.txt.
# Returns a character vector of gene IDs (without mRNA suffix).
load_te_gene_ids <- function() {
  if (!file.exists(TE_IDS_FILE)) {
    msg("  WARNING: TE IDs file not found — is_te_gene will be FALSE for all markers")
    return(character(0))
  }
  mrna_ids <- fread(TE_IDS_FILE, header = FALSE, col.names = "mrna_id")$mrna_id
  gene_ids <- unique(sub("\\.mRNA\\..*$", "", mrna_ids))
  msg("  TE-derived gene IDs loaded: ", length(gene_ids))
  gene_ids
}

# Join functional annotation and TE flag onto any marker/gene table.
# Requires a gene_id column. Adds eggnog_description, pfam_domain,
# panther_description, go_terms, kegg_ko, is_te_gene.
add_functional_annotation <- function(dt, func_annot, te_genes) {
  if (!"gene_id" %in% names(dt)) return(dt)
  if (nrow(func_annot) > 0) {
    dt <- merge(dt, func_annot, by = "gene_id", all.x = TRUE)
  } else {
    dt[, eggnog_description  := NA_character_]
    dt[, pfam_domain         := NA_character_]
    dt[, panther_description := NA_character_]
    dt[, go_terms            := NA_character_]
    dt[, kegg_ko             := NA_character_]
  }
  # Flag gene models with known TE origin — markers in TE genes are noted separately
  dt[, is_te_gene := gene_id %in% te_genes]
  dt
}

# Classify genomic context based on distance to nearest gene.
# Proximal intergenic is defined as ≤2 kb — a conservative promoter-proximal window.
add_annotation_class <- function(dt) {
  dt[, annotation_class := fcase(
    distance_bp == 0L,                           "genic",
    distance_bp > 0L & distance_bp <= 2000L,     "proximal_intergenic",
    distance_bp > 2000L,                          "distal_intergenic",
    default =                                     "unassigned"
  )]
  dt
}

# Load mean beta per site for one cohort × context.
# Returns data.table(site, context, cohort, mean_beta) or NULL if file absent.
load_site_betas <- function(cohort, ctx) {
  if (is.null(BETA_ROOT)) return(NULL)
  f <- file.path(BETA_ROOT,
    sprintf("mean_beta_%s_%s.tsv", tolower(cohort), ctx))
  if (!file.exists(f)) {
    msg("  NOTE: beta file not found: ", basename(f), " — meth_status_absolute = NA")
    return(NULL)
  }
  dt <- fread(f, header = TRUE)
  if (!all(c("site", "mean_beta") %in% names(dt))) {
    msg("  WARNING: ", basename(f), " needs columns: site, mean_beta — skipping")
    return(NULL)
  }
  dt[, site    := as.character(site)]
  dt[, context := ctx]
  dt[, cohort  := cohort]
  dt[, .(site, context, cohort, mean_beta)]
}

# Classify absolute methylation status using context-specific thresholds.
classify_meth_status <- function(beta, ctx) {
  hypo  <- BETA_HYPO [ctx]
  hyper <- BETA_HYPER[ctx]
  ifelse(is.na(beta) | is.na(hypo), NA_character_,
    ifelse(beta < hypo,  "hypomethylated",
      ifelse(beta >= hyper, "hypermethylated", "intermediate")))
}

# Join mean beta values + Option 1 absolute classification onto a sites table.
# site_col: name of the column holding site IDs (default "site").
# The table must already have columns: context, cohort.
add_meth_status <- function(sites_dt, site_col = "site") {
  beta_rows <- list()
  for (co in COHORTS)
    for (ctx in CONTEXTS) {
      b <- load_site_betas(co, ctx)
      if (!is.null(b)) beta_rows[[length(beta_rows) + 1L]] <- b
    }
  beta_all <- if (length(beta_rows)) rbindlist(beta_rows) else data.table()

  if (nrow(beta_all) > 0) {
    dt <- copy(sites_dt)
    if (site_col != "site") setnames(dt, site_col, "site")
    dt <- merge(dt, beta_all, by = c("site", "context", "cohort"), all.x = TRUE)
    if (site_col != "site") setnames(dt, "site", site_col)
    # mapply calls classify_meth_status row-wise over the vectorised table
    dt[, meth_status_absolute := mapply(classify_meth_status, mean_beta, context)]
    return(dt)
  }
  # No beta files available — fill columns with NA so schema remains consistent
  sites_dt[, mean_beta            := NA_real_]
  sites_dt[, meth_status_absolute := NA_character_]
  sites_dt
}

# Back-transform M-values to methylation levels (beta) and return the mean
# beta per site across samples. M = log2[(100*beta + 0.5) / (100 - 100*beta + 0.5)]
# (Step 12ab0), so beta = ((100 + EPS) * 2^M - EPS) / (100 * (1 + 2^M)).
# M-value matrix: rows = samples, columns = site IDs (chr:start-end).
EPSILON_M <- 0.5

compute_mean_beta <- function(cohort, ctx) {
  mat_file <- file.path(MEQTL_INPUTDIR,
    toupper(cohort), ctx, "methylation_mvalues_matrix.rds")
  if (!file.exists(mat_file)) {
    msg("  NOTE: M-value matrix not found: ", mat_file)
    return(NULL)
  }
  msg("  Computing mean betas from M-value matrix: ", cohort, " / ", ctx)
  mmat <- readRDS(mat_file)
  if (is.data.frame(mmat)) mmat <- as.matrix(mmat)
  beta_mat   <- ((100 + EPSILON_M) * 2^mmat - EPSILON_M) / (100 * (1 + 2^mmat))
  beta_means <- colMeans(beta_mat, na.rm = TRUE)
  data.table(site = colnames(mmat), mean_beta = as.numeric(beta_means))
}

# Build gene-centric summary from the combined annotation table.
# func_annot and te_genes are joined to enrich the summary with functional info.
make_gene_summary <- function(dt, func_annot, te_genes) {
  # Rank annotation classes: genic > proximal > distal; report best class per gene
  ann_order <- c(genic = 1L, proximal_intergenic = 2L, distal_intergenic = 3L)
  dt_g <- dt[!is.na(gene_id) & gene_id != ""]   # no gene: NA, or empty when read from file
  if (nrow(dt_g) == 0L) return(data.table())

  # Collapse all markers per gene to one summary row
  summ <- dt_g[, {
    ann_u  <- unique(annotation_class[!is.na(annotation_class)])
    ranked <- ann_u[ann_u %in% names(ann_order)]
    best   <- if (length(ranked)) names(sort(ann_order[ranked]))[1L] else ann_u[1L]
    .(
      gene_chr              = gene_chr[1L],
      gene_start            = gene_start[1L],
      gene_end              = gene_end[1L],
      gene_strand           = gene_strand[1L],
      best_annotation_class = best,
      min_distance_bp       = min(distance_bp, na.rm = TRUE),
      n_markers             = uniqueN(marker_id),
      n_snps                = uniqueN(marker_id[marker_type == "SNP"]),
      n_sites               = uniqueN(marker_id[marker_type == "methylation_site"]),
      marker_ids            = paste(sort(unique(marker_id)),          collapse = ";"),
      sources               = paste(sort(unique(na.omit(source))),    collapse = ";"),
      tools                 = paste(sort(unique(na.omit(tool))),      collapse = ";"),
      cohorts               = paste(sort(unique(na.omit(cohort))),    collapse = ";"),
      contexts              = paste(sort(unique(na.omit(context))),   collapse = ";")
    )
  }, by = gene_id]

  # Join functional annotation and TE flag
  summ <- add_functional_annotation(summ, func_annot, te_genes)

  setorder(summ, gene_chr, gene_start)
  summ
}

############################################################
# 3b) PRE-COMPUTE MEAN BETA FILES (from M-value matrices if not present)
############################################################

msg("Checking / computing mean beta files from M-value matrices...")
dir.create(BETA_ROOT, recursive = TRUE, showWarnings = FALSE)

# Pre-compute beta files up front so all annotation sections can use load_site_betas()
for (.co in COHORTS) {
  for (.ctx in CONTEXTS) {
    out_beta <- file.path(BETA_ROOT,
      sprintf("mean_beta_%s_%s.tsv", tolower(.co), .ctx))
    if (file.exists(out_beta)) {
      msg("  Found: ", basename(out_beta), " — skipping recompute")
      next
    }
    mb <- compute_mean_beta(.co, .ctx)
    if (!is.null(mb)) {
      fwrite(mb, out_beta, sep = "\t")
      msg("  Written: ", basename(out_beta), " (", nrow(mb), " sites)")
    }
  }
}
rm(.co, .ctx)

############################################################
# 4) PREPARE GFF3 GENE BED
############################################################

msg("Preparing GFF3 gene BED (filtering gene features)...")

genes_bed_path <- file.path(TMP_DIR, "genes.bed")

# Extract gene-level features only; create 0-based BED
# Attributes field example: ID=PA_chr01_G000001;
gff_cmd <- sprintf(
  "awk 'BEGIN{OFS=\"\\t\"} !/^#/ && $3==\"gene\" { \
     match($9, /ID=([^;]+)/, a); \
     print $1, $4-1, $5, a[1], \".\", $7 \
   }' %s | sort -k1,1 -k2,2n > %s",
  GFF3_FILE, genes_bed_path)

system(gff_cmd)
n_genes <- as.integer(system(sprintf("wc -l < %s", genes_bed_path), intern = TRUE))
msg("  ", n_genes, " gene features written to genes.bed")

# Load functional annotation and TE gene IDs (used in all subsequent sections)
func_annot <- load_functional_annotation()
te_genes   <- load_te_gene_ids()

############################################################
# 5) TMS KW SVMPs — SELECTED MARKERS FROM STEP 8b
#
# Breeding: top-150 balanced (50 per context, all formal SVMPs)
# Natural:  6 formal SVMPs only (3 CpG + 3 CHH; padj < 0.05)
############################################################

msg("======================================================")
msg("TMS KW SVMPs — annotating selected markers from step 8b")

svmp_files <- list(
  BREEDING = file.path(SVMP_ROOT, "TMS_8B_selected_markers_breeding_top150_BALANCED.tsv"),
  NATURAL  = file.path(SVMP_ROOT, "TMS_8B_selected_markers_natural_formal_SVMPs.tsv")
)

svmp_all_list <- list()

for (cohort in COHORTS) {
  f <- svmp_files[[cohort]]
  if (!file.exists(f)) {
    msg("  Missing: ", basename(f), " — skipping ", cohort)
    next
  }
  dt <- fread(f, header = TRUE)
  msg("  ", cohort, ": ", nrow(dt), " SVMPs loaded")

  # Standardise column names
  if ("Context" %in% names(dt) && !"context" %in% names(dt))
    setnames(dt, "Context", "context")
  dt[, cohort := cohort]
  dt[, pos    := as.integer(start)]   # 1-based genomic position alias

  # BEDTools annotation
  bed_path <- file.path(TMP_DIR,
    sprintf("tms_svmp_%s.bed", tolower(cohort)))
  write_marker_bed(dt, "chr", "pos", "marker_id", bed_path)
  ann <- run_bedtools_closest(bed_path, genes_bed_path)

  if (nrow(ann) > 0)
    dt <- merge(dt,
                ann[, .(marker_id, gene_id, gene_chr, gene_start,
                         gene_end, gene_strand, distance_bp)],
                by = "marker_id", all.x = TRUE)

  dt <- add_annotation_class(dt)

  # Mean beta + hypo/hyper classification from pre-computed beta files.
  # add_meth_status joins on (site = loc, context, cohort).
  dt <- add_meth_status(dt, site_col = "loc")

  # Standard marker metadata columns; SNP-specific fields are NA for methylation sites
  dt[, source                    := "TMS_KW_SVMP"]
  dt[, marker_type               := "methylation_site"]
  dt[, tool                      := NA_character_]
  dt[, DF                        := NA_character_]
  dt[, loading                   := NA_real_]
  dt[, methylation_direction     := NA_character_]
  dt[, higher_methylation_cohort := NA_character_]
  dt[, c("ref","alt","af","dr2") := list(NA_character_, NA_character_,
                                          NA_real_, NA_real_)]
  # meQTL association columns left NA — SVMPs are not from the meQTL analysis
  dt[, n_associated_sites        := NA_integer_]
  dt[, associated_sites          := NA_character_]
  dt[, associated_contexts       := NA_character_]
  dt[, associated_site_chrs      := NA_character_]
  dt[, associated_site_pos       := NA_character_]

  dt <- unique(add_functional_annotation(dt, func_annot, te_genes))

  out_f <- file.path(OUT_ROOT,
    sprintf("tms_svmp_%s_annotated.tsv", tolower(cohort)))
  write_marker_table(dt, out_f)
  msg("  Saved: ", basename(out_f), " (", nrow(dt), " markers)")
  cat("  Annotation class breakdown:\n")
  print(dt[, .N, by = .(context, annotation_class)][order(context)])
  cat("  Methylation status:\n")
  print(dt[, .N, by = .(context, meth_status_absolute)][order(context)])

  svmp_all_list[[cohort]] <- dt
}

if (length(svmp_all_list) == 0)
  msg("  No SVMP files found — section 6b produced no output.")

############################################################
# 6) meQTL ROBUST MARKERS (p_FDR < 1e-10 in BOTH GENESIS AND MatrixEQTL)
############################################################

msg("======================================================")
msg("meQTL robust markers (p_FDR < 1e-10 in GENESIS AND MatrixEQTL)")

meqtl_results <- list()

for (cohort in COHORTS) {
  cohort_key <- tolower(cohort)

  # Robust pairs (pair-level replication in both tools), Step 13ab
  rob_file <- file.path(ROBUST_ROOT,
    sprintf("robust_pairs_%s.tsv", cohort_key))

  if (!file.exists(rob_file)) {
    msg("  WARNING: ", basename(rob_file),
        " not found — run Step 13ab first")
    next
  }

  dt <- fread(rob_file, header = TRUE)

  if ("snp"  %in% names(dt)) dt[, snp  := as.character(snp)]
  if ("site" %in% names(dt)) dt[, site := as.character(site)]
  msg("  ", cohort, ": ", nrow(dt), " robust pairs | ",
      uniqueN(dt$context), " contexts | ",
      uniqueN(dt$snp), " SNPs | ", uniqueN(dt$site), " sites")

  # ---- Associated epimarkers per SNP (aggregated before deduplication) ----
  # Each robust SNP may associate with multiple methylation sites; capture all
  # of them before deduplicating SNP positions for the per-SNP annotation table
  snp_site_map <- dt[!is.na(snp) & !is.na(site), .(
    n_associated_sites   = uniqueN(site),
    associated_sites     = paste(sort(unique(site)),             collapse = ";"),
    associated_contexts  = paste(sort(unique(context)),          collapse = ";"),
    associated_site_chrs = paste(sort(unique(as.character(site_chr))), collapse = ";"),
    associated_site_pos  = paste(sort(unique(as.character(site_pos))), collapse = ";")
  ), by = snp]

  # ---- SNP positions (unique by snp_chr + snp_pos) ----
  snps <- unique(dt[!is.na(snp_chr) & !is.na(snp_pos),
                    .(chr = snp_chr,
                      pos = as.numeric(snp_pos),
                      snp,
                      context)])
  snps[, marker_id := paste0("snp_", snp, "_", chr, ":", as.integer(pos))]

  if (nrow(snps) > 0) {
    bed_path <- file.path(TMP_DIR,
      sprintf("robust_%s_snps.bed", cohort_key))
    write_marker_bed(snps, "chr", "pos", "marker_id", bed_path)
    ann <- run_bedtools_closest(bed_path, genes_bed_path)

    if (nrow(ann) > 0)
      snps <- merge(snps,
                    ann[, .(marker_id, gene_id, gene_chr, gene_start,
                             gene_end, gene_strand, distance_bp)],
                    by = "marker_id", all.x = TRUE)

    # Primary VCF query (imputed — has AF + DR2)
    vi <- query_vcf(
      data.table(chr = snps$chr, pos = as.integer(snps$pos)),
      VCF_FILES[[cohort]],
      file.path(TMP_DIR, sprintf("robust_%s_snps_vcf", cohort_key)))
    if (nrow(vi) > 0)
      snps <- merge(snps, vi, by.x = c("chr","pos"), by.y = c("chr","pos"),
                    all.x = TRUE)
    else
      snps[, c("ref","alt","af","dr2") := list(NA_character_, NA_character_,
                                                NA_real_, NA_real_)]

    # Fill NA ref/alt/af from the cohort VCF before imputation (fallback only)
    snps <- fill_from_raw_vcf(snps, cohort,
      file.path(TMP_DIR, sprintf("robust_%s_snps_rawvcf", cohort_key)))

    # Join associated epimarker information (sites, contexts) aggregated above
    snps <- merge(snps, snp_site_map, by = "snp", all.x = TRUE)

    snps <- add_annotation_class(snps)
    snps[, cohort                    := cohort]
    snps[, tool                      := "GENESIS+MatrixEQTL"]
    snps[, source                    := "robust_meQTL"]
    snps[, marker_type               := "SNP"]
    snps[, DF                        := NA_character_]
    snps[, loading                   := NA_real_]
    snps[, methylation_direction     := NA_character_]
    snps[, mean_beta                 := NA_real_]
    snps[, meth_status_absolute      := NA_character_]
    snps[, higher_methylation_cohort := NA_character_]

    snps <- unique(add_functional_annotation(snps, func_annot, te_genes))
    out_f <- file.path(OUT_ROOT,
      sprintf("robust_%s_snps_annotated.tsv", cohort_key))
    write_marker_table(snps, out_f)
    msg("  Saved: ", basename(out_f), " (", nrow(snps), " unique SNP positions)")
    meqtl_results[[paste(cohort, "snp")]] <- snps
  }

  # ---- Methylation site positions (unique by site_chr + site_pos) ----
  sites <- unique(dt[!is.na(site_chr) & !is.na(site_pos),
                     .(chr = site_chr,
                       pos = as.integer(site_pos),
                       site,
                       context)])
  sites[, marker_id := site]

  if (nrow(sites) > 0) {
    bed_path <- file.path(TMP_DIR,
      sprintf("robust_%s_sites.bed", cohort_key))
    write_marker_bed(sites, "chr", "pos", "marker_id", bed_path)
    ann <- run_bedtools_closest(bed_path, genes_bed_path)

    if (nrow(ann) > 0)
      sites <- merge(sites,
                     ann[, .(marker_id, gene_id, gene_chr, gene_start,
                              gene_end, gene_strand, distance_bp)],
                     by = "marker_id", all.x = TRUE)

    # Methylation sites have no VCF allele information
    sites[, c("ref","alt","af","dr2") := list(NA_character_, NA_character_,
                                               NA_real_, NA_real_)]
    sites <- add_annotation_class(sites)
    sites[, cohort                   := cohort]
    sites[, tool                     := "GENESIS+MatrixEQTL"]
    sites[, source                   := "robust_meQTL"]
    sites[, marker_type              := "methylation_site"]
    sites[, DF                       := NA_character_]
    sites[, loading                  := NA_real_]
    sites[, methylation_direction    := NA_character_]

    # Option 1 — absolute methylation status per site per cohort
    sites <- add_meth_status(sites)

    # higher_methylation_cohort filled in cross-cohort step below
    sites[, higher_methylation_cohort := NA_character_]

    sites <- unique(add_functional_annotation(sites, func_annot, te_genes))
    out_f <- file.path(OUT_ROOT,
      sprintf("robust_%s_sites_annotated.tsv", cohort_key))
    write_marker_table(sites, out_f)
    msg("  Saved: ", basename(out_f), " (", nrow(sites), " unique site positions)")
    meqtl_results[[paste(cohort, "site")]] <- sites
  }

  rm(dt); gc()
}

# Option 2 — cross-cohort relative methylation comparison.
# For each site present in both cohorts, flag which has higher mean beta.
msg("--- Cross-cohort methylation direction (Option 2) ---")
{
  breed_s <- meqtl_results[["BREEDING site"]]
  nat_s   <- meqtl_results[["NATURAL site"]]
  has_beta_b <- !is.null(breed_s) && "mean_beta" %in% names(breed_s)
  has_beta_n <- !is.null(nat_s)   && "mean_beta" %in% names(nat_s)

  if (has_beta_b && has_beta_n) {
    # one value per site and context (a site in two gene models has two rows)
    breed_b <- unique(breed_s[!is.na(mean_beta), .(site, context, beta_b = mean_beta)])
    nat_b   <- unique(nat_s  [!is.na(mean_beta), .(site, context, beta_n = mean_beta)])
    # Full outer join: include sites unique to one cohort as well as shared sites
    cross   <- merge(breed_b, nat_b, by = c("site", "context"), all = TRUE)
    cross[, higher_methylation_cohort := fcase(
      !is.na(beta_b) & !is.na(beta_n) & beta_b > beta_n, "BREEDING",
      !is.na(beta_b) & !is.na(beta_n) & beta_b < beta_n, "NATURAL",
      !is.na(beta_b) & !is.na(beta_n),                   "equal",
      !is.na(beta_b) & is.na(beta_n),                    "BREEDING_only",
      is.na(beta_b)  & !is.na(beta_n),                   "NATURAL_only",
      default = NA_character_
    )]
    cross[, c("beta_b", "beta_n") := NULL]

    # Update per-cohort site tables with the cross-cohort direction column
    for (co in COHORTS) {
      key <- paste(co, "site")
      if (!is.null(meqtl_results[[key]])) {
        meqtl_results[[key]] <- merge(
          meqtl_results[[key]][, higher_methylation_cohort := NULL],
          cross, by = c("site", "context"), all.x = TRUE)
        out_f <- file.path(OUT_ROOT,
          sprintf("robust_%s_sites_annotated.tsv", tolower(co)))
        write_marker_table(meqtl_results[[key]], out_f)
        msg("  Updated: ", basename(out_f), " (higher_methylation_cohort added)")
      }
    }
    msg("  Sites in both cohorts: ",
        cross[!is.na(higher_methylation_cohort) &
              higher_methylation_cohort %in% c("BREEDING","NATURAL","equal"), .N])
  } else {
    msg("  Skipping: beta files absent for one or both cohorts")
  }
}

############################################################
# 7) COMBINED TABLE
############################################################

msg("======================================================")
msg("Building combined annotation table...")

# Common schema across all marker sources; missing columns are added as NA
std_cols <- c("marker_id","marker_type","source","tool","cohort","context","DF",
              "chr","pos","loading",
              "methylation_direction","mean_beta","meth_status_absolute",
              "higher_methylation_cohort",
              "n_associated_sites","associated_sites","associated_contexts",
              "associated_site_chrs","associated_site_pos",
              "ref","alt","af","dr2",
              "gene_id","gene_chr","gene_start","gene_end","gene_strand",
              "distance_bp","annotation_class",
              "eggnog_description","pfam_domain","panther_description",
              "go_terms","kegg_ko","is_te_gene")

# Coerce each source table to the shared schema before row-binding
prep_for_combine <- function(dt, extra_id_col = NULL, pos_col = "pos",
                              chr_col = "chr", loading_col = "loading") {
  out <- copy(dt)
  # Ensure all standard columns exist
  for (col in std_cols)
    if (!col %in% names(out)) out[, (col) := NA]
  if (!is.null(extra_id_col) && extra_id_col %in% names(out) &&
      !"marker_id" %in% names(out))
    out[, marker_id := get(extra_id_col)]
  out[, pos := as.integer(get(pos_col))]
  out[, .SD, .SDcols = intersect(std_cols, names(out))]
}

all_list <- list()

for (key in names(meqtl_results))
  all_list[[key]] <- prep_for_combine(meqtl_results[[key]])

for (cohort in names(svmp_all_list))
  all_list[[paste("SVMP", cohort)]] <- prep_for_combine(svmp_all_list[[cohort]])

if (length(all_list) > 0) {
  combined <- rbindlist(all_list, fill = TRUE, use.names = TRUE)
  setorder(combined, source, cohort, context, chr, pos, na.last = TRUE)

  fwrite(combined, file.path(OUT_ROOT, "all_markers_annotated.tsv"), sep = "\t")
  msg("Combined table: ", nrow(combined), " rows -> all_markers_annotated.tsv")
  cat("\nAnnotation class breakdown:\n")
  print(combined[, .N, by = .(source, annotation_class)][order(source)])
} else {
  msg("WARNING: no markers were annotated — check that input files exist.")
}

############################################################
# 8) GENE-CENTRIC SUMMARY TABLE
############################################################

msg("======================================================")
msg("Building gene-centric summary table...")

# Re-read the combined table to avoid holding all marker data in memory simultaneously
comb_path <- file.path(OUT_ROOT, "all_markers_annotated.tsv")
if (file.exists(comb_path)) {
  combined_for_genes <- fread(comb_path, header = TRUE)
  gene_summ <- make_gene_summary(combined_for_genes, func_annot, te_genes)

  if (nrow(gene_summ) > 0) {
    out_gs <- file.path(OUT_ROOT, "gene_summary_annotated.tsv")
    fwrite(gene_summ, out_gs, sep = "\t")
    msg("Gene summary: ", nrow(gene_summ), " genes → ", basename(out_gs))
    cat("\nGene summary — annotation class breakdown:\n")
    print(gene_summ[, .N, by = best_annotation_class])
    cat("\nGene summary — source breakdown:\n")
    print(gene_summ[, .(n_genes = .N), by = sources][order(-n_genes)])
  } else {
    msg("  No annotated genes found — gene summary not written.")
  }
} else {
  msg("  Combined table not found — skipping gene summary (run section 8 first).")
}

############################################################
# 9) GENES BY METHYLATION-CONTEXT COMBINATION
############################################################
# Source genes: the combination of contexts in which the gene's robust SNPs were
# found, pooled over both cohorts, and per cohort the numbers of unique SNPs and of
# unique methylation sites associated with them. Target genes: the combination of
# contexts of the gene's robust methylation sites and the numbers of unique sites.

msg("======================================================")
msg("Building gene tables by methylation-context combination...")
ctx_combination <- function(x) paste(intersect(CONTEXTS, x), collapse = ";")
snp_all  <- rbindlist(lapply(COHORTS, function(co) meqtl_results[[paste(co, "snp")]]),  fill = TRUE)
site_all <- rbindlist(lapply(COHORTS, function(co) meqtl_results[[paste(co, "site")]]), fill = TRUE)
if (nrow(snp_all) > 0) {
  src <- snp_all[!is.na(gene_id), .(
    contexts           = ctx_combination(unique(context)),
    snps_breeding      = uniqueN(snp[cohort == "BREEDING"]),
    snps_natural       = uniqueN(snp[cohort == "NATURAL"]),
    sites_breeding     = uniqueN(unlist(strsplit(associated_sites[cohort == "BREEDING"], ";"))),
    sites_natural      = uniqueN(unlist(strsplit(associated_sites[cohort == "NATURAL"], ";"))),
    eggnog_description = eggnog_description[1]), by = gene_id]
  fwrite(src, file.path(OUT_ROOT, "source_genes_by_context.tsv"), sep = "\t")
  msg("Source genes: ", nrow(src))
}
if (nrow(site_all) > 0) {
  tgt <- site_all[!is.na(gene_id), .(
    contexts           = ctx_combination(unique(context)),
    sites_breeding     = uniqueN(site[cohort == "BREEDING"]),
    sites_natural      = uniqueN(site[cohort == "NATURAL"]),
    annotation_class   = annotation_class[1],
    eggnog_description = eggnog_description[1]), by = gene_id]
  fwrite(tgt, file.path(OUT_ROOT, "target_genes_by_context.tsv"), sep = "\t")
  msg("Target genes: ", nrow(tgt))
}

# Clean up temp files
unlink(TMP_DIR, recursive = TRUE)

msg("======================================================")
msg("Step 17ab finished. Outputs: ", OUT_ROOT)
msg("======================================================")
```

---
