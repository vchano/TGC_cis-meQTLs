# JOINT — Integrated ECS + TMS Analyses

Joint analyses combining SNP genotype (ECS) and DNA methylation (TMS) data for
cis-meQTL mapping in *Picea abies*. Requires ECS steps 1a–10a and TMS steps
1b–8b to be complete before running.

## Steps

| Step | Script(s) | Input | Output | Tool |
|------|-----------|-------|--------|------|
| 11ab | `11ab.tgc.joint.correlation.analysis.R` | ECS IBS distance, TMS Euclidean distance | Supp. Fig. S1 — Procrustes panel; Mantel/Procrustes/RV tables | vegan, ade4, ggplot2 |
| 12ab0 | `12ab0.tgc.joint.meqtl.input.prep.R` | GDS, GRM, methylKit objects | M-value matrices, SNP annotation, PCs, GRM (RDS + TSV) | data.table, SeqArray |
| 12ab1 | `12ab1.tgc.joint.matrixeqtl.mapping.R` | outputs of 12ab0 | cis-meQTL results per panel (RDS) | MatrixEQTL |
| 12ab2 | `12ab2.tgc.joint.genesis.mapping.R` | outputs of 12ab0 | cis-meQTL results per panel (RDS) | GENESIS |
| 13ab | `13ab.tgc.joint.meqtl.combined.results.R` | outputs of 12ab1 + 12ab2 | Table 1, Supp. Figs. S2–S3 (QQ plots), Supp. Table S5 (robust pairs) | data.table, openxlsx2 |
| 14ab | `14ab.tgc.joint.manhattan.plots.R` / `.sh` | outputs of 12ab1 + 12ab2 | Figure 3 (GENESIS), Extended Data Fig. 4 (MatrixEQTL) — circular Manhattan + combined panels | circlize, magick |
| 15ab | `15ab.tgc.joint.venn.overlap.R` / `.sh` | significant site lists from 13ab | Figure 4, Supp. Fig. S4 — Venn diagrams, overlap tables | ggvenn, data.table |
| 16ab | `16ab.tgc.joint.meth.heritability.R` / `.sh` | GRM, M-value matrices, robust meQTL sites from 15ab | Figure 5 — SNP-based (and pedigree-based) methylation heritability | data.table, ggplot2 |
| 17ab | `17ab.tgc.joint.marker.annotation.R` / `.sh` | SVMPs (8b), robust markers (15ab), reference GFF3 | Table 2, Table 3, Supp. Table S3, S6, S7 — annotated marker tables | data.table, openxlsx2 |

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
| 14ab | 4 | ~6 h | FDR_AXIS=TRUE on NATURAL/CHH needs ~96 GB RAM |
| 15ab | 8 | ~1 h | |
| 16ab | 8 | ~1–several h | Distribute across nodes for >100k sites (see .sh) |
| 17ab | 8 | ~2 h | Requires reference GFF3 |

## Dependencies

```
module load gcc/14.2.0
module load r/4.5.2
module load imagemagick/7.1.1-39
```

R packages: `data.table`, `GENESIS`, `MatrixEQTL`, `SeqArray`, `SeqVarTools`,
`SNPRelate`, `circlize`, `ggplot2`, `patchwork`, `ggvenn`, `openxlsx2`, `vegan`, `ade4`
