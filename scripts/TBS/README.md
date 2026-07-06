# TBS — Targeted Bisulfite Sequencing Pipeline

DNA methylation quantification from bisulfite sequencing data for the breeding and
natural *Picea abies* cohorts, across three cytosine contexts (CpG, CHG, CHH).

## Steps

| Step | Script | Input | Output | Tool |
|------|--------|-------|--------|------|
| 1b | `1b.tgc.tbs.fastqc.rawdata.sh` | raw fastq.gz | FastQC/MultiQC reports | FastQC, MultiQC |
| 2b | `2b.tgc.tbs.trimmomatic.and.fastqc.trimmed.sh` | raw fastq.gz | trimmed fastq.gz + QC reports | Trimmomatic, FastQC |
| 3b | `3b.tgc.tbs.bismark.trimmed.sh` | trimmed fastq.gz | BAM (bisulfite-aligned) | Bismark, Bowtie2 |
| 4b1 | `4b1.tgc.tbs.meth.extractor.array.sh` | BAM | CpG/CHG/CHH coverage files | Bismark methylation extractor |
| 4b2 | `4b2.tgc.tbs.meth.filtering.array.sh` | coverage files | filtered coverage files (≥ 5×) | awk/bash |
| 5b | `5b.tgc.methylkit.filtering.sh` / `5b.tgc.tbs.import.data.methylkit.R` | filtered coverage | methylKit objects per context | methylKit |
| 6b | `6b.tgc.tbs.methylation.levels.anova.R` | methylKit objects | methylation level ANOVA results + figures | R base, ggplot2 |
| 7b | `7b.tgc.tbs.pcs.dapc.heatmaps.R` | methylKit objects | PCA, DAPC, heatmap figures | adegenet, ggplot2 |
| 8b | `8b.tgc.tbs.heatmaps.R` | methylKit objects, DAPC | heatmap figures (TIFF/PDF) | ggplot2, pheatmap |

## Compute requirements

| Step | Cores | Walltime | Notes |
|------|-------|----------|-------|
| 1b | 48 | ~4 h | |
| 2b | 48 | ~12 h | SLURM array recommended |
| 3b | 96 | ~48 h | Exclusive node; plates run sequentially or as array |
| 4b1 | 48 | ~24 h | SLURM array per sample |
| 4b2 | 16 | ~6 h | SLURM array per sample |
| 5b | 16 | ~8 h | |
| 6b–8b | 8 | ~1–4 h | CHH datasets are large; allow extra memory |

## Notes on large datasets

- CHH context produces very large methylation matrices (tens of millions of sites).
  Steps 6b–8b may require 64–128 GB RAM for natural cohort CHH.
- QQ plots for CHH panels should be written as TIFF only (not EPS/PNG) to avoid
  out-of-memory errors with millions of scatter points.

## Dependencies

```
module load gcc/14.2.0
module load r/4.5.2
module load miniforge3 bowtie2 samtools
# Bismark must be installed and accessible in PATH
```

R packages: `methylKit`, `data.table`, `ggplot2`, `adegenet`, `pheatmap`
