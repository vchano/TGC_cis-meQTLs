# TGC_cis-meQTLs
Cis-meQTL mapping in Norway spruce (Picea abies) integrating exome capture sequencing and targeted bisulfite sequencing across breeding and natural cohorts. Scripts for the TreeGeneClimate (TGC) project.

# Cis-meQTL Mapping in Norway Spruce (*Picea abies*)

## Overview

This repository contains all scripts used for cis-methylation quantitative trait
locus (cis-meQTL) mapping in *Picea abies* (Norway spruce), integrating:

- **ECS** — Exome Capture Sequencing (SNP genotyping)
- **TBS** — Targeted Bisulfite Sequencing (DNA methylation)

Two cohorts were analyzed: a **breeding cohort** (17 full-sib families) and a
**natural cohort** (25 Finnish natural stands). Methylation was quantified in
three cytosine contexts (CpG, CHG, CHH). Cis-meQTLs were mapped within a
100 kb window using two complementary approaches: a GRM-based mixed model
(GENESIS) and a linear model (MatrixEQTL).

This is part of the **TreeGeneClimate (TGC)** project.

## Repository Structure

```
scripts/
├── ECS/     # Exome capture sequencing: QC → trimming → mapping →
             #   variant calling → filtering → imputation → pop. structure
├── TBS/     # Bisulfite sequencing: QC → trimming → mapping →
             #   methylation extraction → filtering → methylation analysis
└── JOINT/   # Joint analyses: meQTL mapping, Manhattan plots,
             #   overlap analysis, figure/table export
```

## Pipeline Steps

### ECS (SNP genotyping)

| Step | Script | Description |
|------|--------|-------------|
| 1a | `1a.tgc.ecs.fastqc.rawdata.sh` | FastQC on raw reads |
| 2a | `2a.tgc.ecs.trimmomatic.and.fastqc.trimmed.sh` | Trimmomatic + FastQC |
| 3a | `3a.tgc.ecs.bowtie.trimmed.sh` | Bowtie2 alignment |
| 4a1 | `4a1.tgc.ecs.snv.calling.bcftools.sh` | SNV calling (bcftools) |
| 4a2 | `4a2.tgc.ecs.snv.filtering.bcftools.vcftools.sh` | SNV filtering |
| 5a | `5a.tgc.ecs.plink.admixture.split.gwasprep.sh` | PLINK + ADMIXTURE |
| 6a | `6a.tgc.ecs.beagle.imputation.sh` | BEAGLE imputation |
| 7a | `7a.tgc.ecs.plink.ibd.sh` | IBD estimation |
| 8a | `8a.tgc.ecs.orchestrator.R` | GDS construction |
| 9a | `9a.tgc.ecs.pca.R` | Principal component analysis |
| 10a | `10a.tgc.ecs.grm.ibd.dapc.biplot.R` | GRM, IBD, DAPC |

### TBS (DNA methylation)

| Step | Script | Description |
|------|--------|-------------|
| 1b | `1b.tgc.tbs.fastqc.rawdata.sh` | FastQC on raw reads |
| 2b | `2b.tgc.tbs.trimmomatic.and.fastqc.trimmed.sh` | Trimmomatic + FastQC |
| 3b | `3b.tgc.tbs.bismark.trimmed.sh` | Bismark alignment |
| 4b1 | `4b1.tgc.tbs.meth.extractor.array.sh` | Methylation extraction |
| 4b2 | `4b2.tgc.tbs.meth.filtering.array.sh` | Coverage filtering |
| 5b | `5b.tgc.methylkit.filtering.sh` / `5b.tgc.tbs.import.data.methylkit.R` | methylKit import |
| 6b | `6b.tgc.tbs.methylation.levels.anova.R` | Methylation level ANOVAs |
| 7b | `7b.tgc.tbs.pcs.dapc.heatmaps.R` | PCA, DAPC, heatmaps |
| 8b | `8b.tgc.tbs.heatmaps.R` | Heatmap figures |

### Joint analyses

| Step | Script | Description |
|------|--------|-------------|
| 11ab | `11ab.tgc.joint.correlation.analysis.R` | ECS–TBS correlation |
| 16ab0–2 | `16ab0.R` – `16ab2.R` | meQTL input preparation |
| 16ab3 | `16ab3.R` | meQTL mapping (GENESIS + MatrixEQTL) |
| 16ab4 | `16ab4.R` | Circular Manhattan plots |
| 16ab5 | `16ab5.R` / `16ab5.sh` | Venn overlap analysis |
| 16ab6 | `16ab6.R` / `16ab6.sh` | Panel figure assembly |
| 17ab | `17ab.R` / `17ab.sh` | Additional summaries |
| 19ab | `19ab_export_figures_tables.R` / `19ab_export_figures_tables.sh` | Figure and table export |
| 20ab | `20ab_summary_tables.R` | ECS/TBS summary tables |

## Dependencies

**Tools (via SLURM modules):**
- `gcc/14.2.0`, `r/4.5.2`, `imagemagick/7.1.1-39`
- FastQC, Trimmomatic, Bowtie2, Bismark, bcftools, vcftools, PLINK, BEAGLE

**R packages:**
- `data.table`, `GENESIS`, `MatrixEQTL`, `SeqArray`, `SeqVarTools`
- `ggplot2`, `circlize`, `methylKit`, `SNPRelate`, `MASS`, `adegenet`

## Data Availability

| Data type | Repository | Accession |
|-----------|-----------|-----------|
| Raw ECS fastq (FASTQ) | NCBI SRA | [pending] |
| Raw TBS fastq (FASTQ) | NCBI SRA / GEO | [pending] |
| Methylation matrices | NCBI GEO | [pending] |
| Filtered VCF, GDS, meQTL inputs/results | GRO.data | [pending] |

## Citation

> [Authors]. [Year]. [Title]. [Journal]. doi:[pending]

## Authors

[Name and affiliations]

## License

MIT License — see [LICENSE](LICENSE) for details.
