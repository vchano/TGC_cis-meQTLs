# Cis-meQTL Mapping in Norway Spruce (*Picea abies*)

## Overview

This repository contains all scripts used for cis-methylation quantitative trait
locus (cis-meQTL) mapping in *Picea abies* (Norway spruce), integrating:

- **ECS** — Exome Capture Sequencing (SNP genotyping)
- **TMS** — Target-enriched Enzymatic Methylation Sequencing (DNA methylation)

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
├── TMS/     # Enzymatic Methylation sequencing: QC → trimming → mapping →
             #   methylation extraction → filtering → methylation analysis
└── JOINT/   # Joint analyses: meQTL mapping, Manhattan plots,
             #   overlap analysis, annotation, summary tables
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
| 8a | `8a.tgc.ecs.orchestrator.R` | GDS construction, GRM, kinship |
| 9a | `9a.tgc.ecs.pca.R` | Figure 1 — genetic PCA panel |
| 10a | `10a.tgc.ecs.grm.ibd.R` | Extended Data Fig. 1 — GRM + IBD panel |

### TMS (DNA methylation)

| Step | Script | Description |
|------|--------|-------------|
| 1b | `1b.tgc.tms.fastqc.rawdata.sh` | FastQC on raw reads |
| 2b | `2b.tgc.tms.trimmomatic.and.fastqc.trimmed.sh` | Trimmomatic + FastQC |
| 3b | `3b.tgc.tms.bismark.trimmed.sh` | Bismark alignment |
| 4b1 | `4b1.tgc.tms.meth.extractor.array.sh` | Methylation extraction |
| 4b2 | `4b2.tgc.tms.meth.filtering.array.sh` | Coverage filtering |
| 5b | `5b.tgc.methylkit.filtering.sh` / `5b.tgc.tms.import.data.methylkit.R` | methylKit import and filtering |
| 6b | `6b.tgc.tms.methylation.levels.R` | Extended Data Fig. 2 — methylation level panel + ANOVAs |
| 7b | `7b.tgc.tms.pca.R` | Figure 2 — epigenomic PCA panel |
| 8b | `8b.tgc.tms.heatmaps.R` | Extended Data Fig. 3 — SVMP heatmap panel |

### Joint analyses

| Step | Script | Description |
|------|--------|-------------|
| 11ab | `11ab.tgc.joint.correlation.analysis.R` | ECS–TMS correlation |
| 12ab0 | `12ab0.tgc.joint.meqtl.input.prep.R` | meQTL input preparation (M-values, GRM, PCs) |
| 12ab1 | `12ab1.tgc.joint.matrixeqtl.mapping.R` | MatrixEQTL cis-meQTL mapping (linear model) |
| 12ab2 | `12ab2.tgc.joint.genesis.mapping.R` | GENESIS cis-meQTL mapping (LMM + GRM) |
| 13ab | `13ab.tgc.joint.meqtl.combined.results.R` | Table 1, Supp. Figs. S2–S3, Supp. Table S5 — combined results, QQ plots, robust pairs |
| 14ab | `14ab.tgc.joint.manhattan.plots.R` / `.sh` | Figure 3, Extended Data Fig. 4 — circular Manhattan plots + combined panels |
| 15ab | `15ab.tgc.joint.venn.overlap.R` / `.sh` | Figure 4, Supp. Fig. S4 — overlap analysis across tools and contexts |
| 16ab | `16ab.tgc.joint.meth.heritability.R` / `.sh` | Figure 5 — SNP-based methylation heritability |
| 17ab | `17ab.tgc.joint.marker.annotation.R` / `.sh` | Table 2, Table 3, Supp. Table S3, S6, S7 — marker annotation against reference GFF3 |

## Dependencies

**Tools (via HPC modules):**
- FastQC, MultiQC, Trimmomatic
- Bowtie2, Bismark, SAMtools
- bcftools, vcftools, PLINK 1.9, BEAGLE
- R 4.5.2, ImageMagick

**R packages:**
- `data.table`, `SNPRelate`, `gdsfmt`, `SeqArray`, `SeqVarTools`
- `GENESIS`, `MatrixEQTL`
- `methylKit`
- `ggplot2`, `circlize`, `ComplexHeatmap`, `patchwork`, `MASS`

## Quick start

1. Set `PROJECT_ROOT` at the top of each script to your local project directory.
2. Adapt SLURM directives (`--account`, `--partition`, `--mail-user`) for your cluster.
3. Run scripts in the order listed above. ECS and TMS steps (1a–10a, 1b–8b) can
   be run in parallel; joint steps (11ab–17ab) require both to be complete first.

## Data availability

| Data type | Repository | Accession |
|-----------|-----------|-----------|
| Raw ECS fastq (FASTQ) | NCBI SRA | [pending] |
| Raw TMS fastq (FASTQ) | NCBI SRA | [pending] |
| Methylation matrices | GRO.data | [pending] |
| Filtered VCF, GDS, meQTL inputs/results | GRO.data | [pending] |

## Citation

> [Authors]. [Year]. [Title]. [Journal]. doi:[pending]

## Authors

[Your name and affiliations here]

## License

MIT License — see [LICENSE](LICENSE) for details.
