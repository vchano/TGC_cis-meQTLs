# JOINT — Integrated ECS + TBS Analyses

Joint analyses combining SNP genotype (ECS) and DNA methylation (TBS) data for
cis-meQTL mapping in *Picea abies*. Requires ECS steps 1a–10a and TBS steps
1b–8b to be complete before running.

## Steps

| Step | Script(s) | Input | Output | Tool |
|------|-----------|-------|--------|------|
| 11ab | `11ab.tgc.joint.correlation.analysis.R` | ECS PCs, TBS PCs | Procrustes / correlation figures | vegan, ggplot2 |
| 12ab0 | `12ab0.tgc.joint.meqtl.input.prep.R` | GDS, GRM, methylKit objects | M-value matrices, SNP annotation, PCs, GRM (RDS + TSV) | data.table, SeqArray |
| 12ab1 | `12ab1.tgc.joint.matrixeqtl.mapping.R` | outputs of 12ab0 | cis-meQTL results per panel (RDS) | MatrixEQTL |
| 12ab2 | `12ab2.tgc.joint.genesis.mapping.R` | outputs of 12ab0 | cis-meQTL results per panel (RDS) | GENESIS |
| 13ab | `13ab.tgc.joint.meqtl.combined.results.R` | outputs of 12ab1 + 12ab2 | comprehensive summary TSV, QQ plots (TIFF), significant site lists | data.table |
| 14ab | `14ab.tgc.joint.manhattan.plots.R` | outputs of 12ab1 + 12ab2 | circular Manhattan plots (TIFF/PDF) | circlize |
| 15ab | `15ab.tgc.joint.venn.overlap.R` / `.sh` | significant site lists from 13ab | Venn diagrams, overlap tables | VennDiagram, data.table |
| 16ab | `16ab.tgc.joint.panel.assembly.R` / `.sh` | TIFF figures | multi-panel figures (TIFF/PDF) | ImageMagick |
| 17ab | `17ab.tgc.joint.marker.annotation.R` / `.sh` | significant markers, reference GFF3 | annotated marker tables (TSV) | data.table |
| 18ab | `18ab.tgc.joint.summary.tables.R` | pipeline outputs, MultiQC stats | manuscript summary tables (TSV + markdown) | data.table |

## Key parameters

- **Cis window:** 100 kb
- **FDR thresholds:** 5×10⁻⁸ (loose) and 1×10⁻¹⁰ (strict), BH-adjusted
- **Panels:** 2 cohorts (BREEDING, NATURAL) × 3 contexts (CpG, CHG, CHH) = 6 panels per tool
- **Models:**
  - MatrixEQTL: `M-value ~ SNP + PC1..PC10` (linear)
  - GENESIS: `M-value ~ PC1..PC10 + random(GRM)` (LMM, AIREML)
    - BREEDING: LOCO-GRM; NATURAL: full GRM

## Compute requirements

| Step | Cores | Walltime | Notes |
|------|-------|----------|-------|
| 11ab | 8 | ~1 h | |
| 12ab0 | 8 | ~2 h | |
| 12ab1 | 48 | ~6 h per panel | Run as SLURM array (6 panels) |
| 12ab2 | 48 | ~24 h per panel | AIREML is slower; BREEDING CHH ~48 h |
| 13ab | 16 | ~2 h | |
| 14ab | 8 | ~2 h | |
| 15ab–16ab | 8 | ~1 h | |
| 17ab | 8 | ~2 h | Requires reference GFF3 |
| 18ab | 8 | ~30 min | |

## Dependencies

```
module load gcc/14.2.0
module load r/4.5.2
module load imagemagick/7.1.1-39
```

R packages: `data.table`, `GENESIS`, `MatrixEQTL`, `SeqArray`, `SeqVarTools`,
`SNPRelate`, `circlize`, `ggplot2`, `VennDiagram`
