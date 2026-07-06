# ECS — Exome Capture Sequencing Pipeline

SNP genotyping from exome capture data for the breeding and natural *Picea abies* cohorts.

## Steps

| Step | Script | Input | Output | Tool |
|------|--------|-------|--------|------|
| 1a | `1a.tgc.ecs.fastqc.rawdata.sh` | raw fastq.gz | FastQC/MultiQC reports | FastQC, MultiQC |
| 2a | `2a.tgc.ecs.trimmomatic.and.fastqc.trimmed.sh` | raw fastq.gz | trimmed fastq.gz + QC reports | Trimmomatic, FastQC |
| 3a | `3a.tgc.ecs.bowtie.trimmed.sh` | trimmed fastq.gz | BAM (sorted, deduplicated) | Bowtie2, SAMtools, Picard |
| 4a1 | `4a1.tgc.ecs.snv.calling.bcftools.sh` | BAM files | raw VCF | bcftools |
| 4a2 | `4a2.tgc.ecs.snv.filtering.bcftools.vcftools.sh` | raw VCF | filtered VCF (MAF ≥ 0.05) | bcftools, vcftools |
| 5a | `5a.tgc.ecs.plink.admixture.split.gwasprep.sh` | filtered VCF | PLINK bed/bim/fam, ADMIXTURE Q files | PLINK 1.9, ADMIXTURE |
| 6a | `6a.tgc.ecs.beagle.imputation.sh` | filtered VCF (per cohort) | imputed VCF | BEAGLE |
| 7a | `7a.tgc.ecs.plink.ibd.sh` | PLINK bed (LD-pruned) | IBD pairwise estimates (.genome) | PLINK 1.9 |
| 8a | `8a.tgc.ecs.orchestrator.R` | imputed VCF, IBD, PCA, ADMIXTURE | GDS, GRM, kinship matrices (RDS + TSV) | SNPRelate, gdsfmt |
| 9a | `9a.tgc.ecs.pca.R` | GDS | PCA figures (TIFF/PDF) | SNPRelate, ggplot2 |
| 10a | `10a.tgc.ecs.grm.ibd.dapc.biplot.R` | GRM, IBD, ADMIXTURE | GRM/IBD/DAPC figures (TIFF/PDF) | adegenet, ggplot2 |

## Compute requirements

| Step | Cores | Walltime | Notes |
|------|-------|----------|-------|
| 1a | 48 | ~4 h | MultiQC is sequential |
| 2a | 48 | ~12 h | SLURM array recommended |
| 3a | 48 | ~24 h | One job per sample |
| 4a1 | 48 | ~48 h | Region-split calling |
| 4a2 | 8 | ~2 h | |
| 5a | 16 | ~6 h | ADMIXTURE K=2–10 |
| 6a | 32 | ~24 h | Per-cohort |
| 7a | 8 | ~2 h | |
| 8a–10a | 8 | ~1–2 h | Interactive or single node |

## Dependencies

```
module load gcc/14.2.0
module load r/4.5.2
module load fastqc trimmomatic bowtie2 samtools bcftools vcftools plink/1.9 beagle
```

R packages: `SNPRelate`, `gdsfmt`, `dplyr`, `tibble`, `readr`, `ggplot2`, `adegenet`, `MASS`
