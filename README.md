# TGC_cis-meQTLs

Analysis code for:

> Chano V, Krutovsky KV, Wang K, Haapanen M, Asiegbu FO, Gailing O. *Cis-meQTL mapping reveals genetic
> architecture of DNA methylation in a giga-genome species.*

The pipeline takes exome capture sequencing (ECS) and targeted enzymatic methylation sequencing (TMS)
reads from 602 Norway spruce (*Picea abies*) trees (breeding cohort, 17 full-sib families, n = 209;
natural cohort, 25 stands, n = 393) to genotypes, methylation calls, cis-meQTL maps, heritability
estimates and marker annotation.

## Layout

```
scripts/
  ECS/     1a-10a    reads -> SNP calls -> filtering, SNP mask, cohort split, imputation, IBD, GRM, PCA
  TMS/     1b-8b     reads -> Bismark -> deduplication -> methylKit datasets -> methylation levels,
                     PCA, SVMP heatmaps
  JOINT/   11ab-17ab concordance, cis-meQTL mapping (MatrixEQTL, GENESIS), robust pairs, plots,
                     heritability, annotation
```

Every script states its inputs and outputs in its header. R scripts are run with
`Rscript --vanilla <script>` from a project root with `DATA/`, `REFERENCE/`, `RESULTS/` and `LOGS/`,
with `TGC_PROJECT_ROOT` set to it; the `.sh` scripts are SLURM job scripts, whose module names and
resource requests need adapting to other clusters. Results are written as tab-separated tables; each
figure is written with a table of the values it shows (`*_values.tsv`, same folder).

## Order of execution

| Step | Script | What it does |
|---|---|---|
| 1a–3a | `ECS/1a…`, `2a…`, `3a…` | FastQC/MultiQC, Trimmomatic, Bowtie2 + SAMtools (fixmate, sort, markdup) |
| 4a1 | `ECS/4a1…` | Joint SNV calling (BCFtools mpileup/call), all samples |
| 4a2 | `ECS/4a2…` | Sample renaming; biallelic SNPs, MAF ≥ 0.05, missing ≤ 20%, site depth ≥ 10 |
| 4a3 | `ECS/4a3…` | SNP positions used to mask methylation calls |
| 5a | `ECS/5a…` | Cohort split, LD pruning (PLINK 50 10 0.2), PCA |
| 6a | `ECS/6a…` | Within-cohort polymorphic sites, BEAGLE imputation (nthreads fixed at 48: BEAGLE output depends on the number of threads) |
| 7a | `ECS/7a…` | PLINK IBD |
| 8a–10a | `ECS/8a…`, `9a…`, `10a…` | Genotype GDS (called and imputed), GCTA GRM from called genotypes of SNPs called in ≥ 99% of trees (`ECS/tgc.ecs.functions.R`), IBD/PCA import; genetic PCA panel; GRM heatmaps and IBD violins |
| 1b–3b | `TMS/1b…`, `2b…`, `3b…` | FastQC/MultiQC, Trimmomatic, Bismark alignment |
| 3b2 | `TMS/3b2…` | PCR-duplicate removal (SAMtools fixmate/markdup), one sample per array task; then `TGC_MODE=summary` for the duplicate-read summary |
| 4b1–4b2 | `TMS/4b1…`, `4b2…` | Bismark methylation extractor and coverage files (not used downstream) |
| 5b1–5b3 | `TMS/5b1…`, `5b2…`, `5b3…` | methylKit import from deduplicated BAMs; SNP mask; population-level (coverage 10 reads – 99.9th percentile) and cis-meQTL datasets (`TMS/tgc.tms.functions.R`) |
| 6b–8b | `TMS/6b…`, `7b…`, `8b…` | Methylation levels with Kruskal–Wallis and Dunn's tests; epigenomic PCA; SVMPs and heatmaps |
| 11ab | `JOINT/11ab…` | Genomic–epigenomic concordance (Mantel, Procrustes, RV), distances without imputed values |
| 12ab0–12ab2 | `JOINT/12ab0…`, `12ab1…`, `12ab2…` | cis-meQTL inputs; MatrixEQTL; GENESIS |
| 13ab | `JOINT/13ab…` | Panel summaries, robust pairs, QQ panel |
| 14ab, 15ab | `JOINT/14ab…`, `15ab…` | Circular Manhattan plots; Venn diagrams |
| 16ab | `JOINT/16ab…` | Per-site heritability, whole-methylome Mantel test |
| 17ab | `JOINT/17ab…` | Marker annotation; genes by methylation-context combination |

## Software

R 4.5.2 with data.table, methylKit, SNPRelate, gdsfmt, GENESIS, GWASTools, MatrixEQTL, lme4, vegan,
ade4, ggplot2, patchwork, ComplexHeatmap, circlize, ggvenn, dunn.test, multcompView.
Command-line tools: FastQC, MultiQC, Trimmomatic 0.36, Bowtie2, SAMtools 1.21, BCFtools 1.19,
PLINK 1.9, BEAGLE 5, Bismark, BEDTools 2.31, ImageMagick 7. Steps 8a, 9a, 10a, 7b, 8b, 11ab, 12ab0,
14ab and 15ab write the R session information (package versions) with their outputs.

## Data

Raw reads: NCBI SRA, BioProject PRJNA1226161. Processed genotypes, methylation matrices and cis-meQTL
results: GRO.data (University of Göttingen), https://doi.org/10.25625/NF4CWR,
https://doi.org/10.25625/BIQMWU and https://doi.org/10.25625/IPOUJT.

## Licence

MIT (see `LICENSE`).
