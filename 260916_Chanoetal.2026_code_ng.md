# Code Archive — Cis-meQTL mapping reveals genetically hardwired methylation at immune and chromatin regulatory loci in a giga-genome species

**Victor Chano¹·², Konstantin V. Krutovsky¹·², Kai Wang³·⁴, Matti Haapanen⁵, Fred O. Asiegbu³, Oliver Gailing¹·²**

¹ Department of Forest Genetics and Forest Tree Breeding, University of Göttingen. Büsgenweg 2, 37077 Göttingen, Germany.  
² Center for Integrated Breeding Research (CiBreed), University of Göttingen. Albrecht-Thaer-Weg 3, 37075 Göttingen, Germany.  
³ Department of Forest Science, University of Helsinki. FIN-00014 Helsinki, Finland.  
⁴ College of Forestry, Fujian Agriculture and Forestry University, Fuzhou 350002, China.  
⁵ Natural Resources Institute Finland (LUKE). Latokartanonkaari 9, Helsinki 00790, Finland.  

\* Corresponding authors: Victor Chano, Oliver Gailing

**ORCIDs**  
Victor Chano: 0000-0002-3423-9090  
Konstantin V. Krutovsky: 0000-0002-8819-7084  
Kai Wang: 0000-0002-5205-7346  
Matti Haapanen: 0000-0003-3294-501X  
Fred O. Asiegbu: 0000-0003-0223-7194  
Oliver Gailing: 0000-0002-4572-2408  

---

**Journal:** Nature Genetics (2026)  
**DOI:** *to be assigned*  
**GitHub repository:** https://github.com/vchano/TGC_cis-meQTLs  

---

## Description

This archive contains all analysis scripts used to generate the results, figures, and supplementary tables in Chano et al. (2026). Scripts are organised into three sequential modules corresponding to the data-processing workflow:

- **ECS** (Exome Capture Sequencing): raw data QC, adapter trimming, alignment, SNP calling, filtering, imputation, IBD estimation, PCA, and GRM/DAPC (Steps 1a–10a).
- **TMS** (Targeted Methylation Sequencing): raw data QC, trimming, bisulfite alignment, CpG/CHG/CHH methylation extraction, filtering, and population-level methylation analysis (Steps 1b–8b).
- **JOINT**: epigenomic–genomic concordance analysis, cis-meQTL mapping (MatrixEQTL and GENESIS), combined result tables, circular Manhattan plots, Venn overlaps, SNP heritability estimation, LD-based validation, SNP–gene annotation, and supplementary tables (Steps 11ab–17ab).

All scripts require R ≥ 4.5.2 and/or a SLURM-based HPC cluster. Set `PROJECT_ROOT` at the top of each script before running.

**R packages:** data.table, openxlsx2, methylKit, ggplot2, ComplexHeatmap, circlize, SNPRelate, gdsfmt, GENESIS, GWASTools, MatrixEQTL, vegan, ade4, RColorBrewer, magick, patchwork, car, multcompView  
**Bioinformatics tools:** FastQC, Trimmomatic, Bowtie2, BCFtools, SAMtools, VCFtools, BEAGLE, PLINK2, Bismark, MultiQC, ImageMagick

---


---

## ECS — Exome Capture Sequencing

---

### `10a.tgc.ecs.grm.ibd.dapc.biplot.R`

```r
#!/usr/bin/env Rscript
############################################################
# TreeGeneClimate (TGC) — ECS
# Step 10a: Figures 1 & 2 (a–c) — GRM heatmaps, IBD violins (PI_HAT), DAPC biplots
#
# OUTPUT (FIXED; DO NOT CHANGE):
#   /path/to/your/project/RESULTS/ECS/RANALYSIS/FIGURES/FIG1_FIG2
#   /path/to/your/project/RESULTS/ECS/RANALYSIS/TABLES/dapc_loadings
#
# FIGURE FILES (saved individually):
#   Figure1a_breeding.tiff   (GRM heatmap)
#   Figure1b_breeding.tiff   (IBD violin, PI_HAT)
#   Figure1c_breeding.tiff   (DAPC biplot)
#   Figure2a_natural.tiff
#   Figure2b_natural.tiff
#   Figure2c_natural.tiff
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

  library(adegenet)   # DAPC (Discriminant Analysis of Principal Components)
  library(ggrepel)    # non-overlapping text labels in ggplot biplots
})

options(stringsAsFactors = FALSE)
set.seed(1)  # reproducible DAPC cross-validation and random jitter

# ==============================================================================
# 1) FIXED OUTPUT LOCATIONS (DO NOT CHANGE)
# ==============================================================================
# === USER CONFIGURATION ===
# Set TGC_PROJECT_ROOT as an environment variable, or edit the fallback path below
PROJECT_ROOT <- Sys.getenv("TGC_PROJECT_ROOT",
  unset = "/path/to/your/project")
# ===========================
RANA_DIR     <- file.path(PROJECT_ROOT, "RESULTS/ECS/RANALYSIS")
RDATA_DIR    <- file.path(RANA_DIR, "RDATA")   # RDS objects written by step 8a

FIG_DIR      <- file.path(RANA_DIR, "FIGURES", "FIG1_FIG2")      # <-- fixed
TAB_DAPC_LD  <- file.path(RANA_DIR, "TABLES", "dapc_loadings")   # <-- fixed
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DAPC_LD, recursive = TRUE, showWarnings = FALSE)

# Helper: stop with an informative message if a required file is absent
ensure_file <- function(p) if (!file.exists(p)) stop("Missing file: ", p, call. = FALSE)

# ==============================================================================
# 2) INPUTS (from 8a)
# ==============================================================================
# VanRaden GRM matrices (sample × sample, symmetric) — one per cohort
grm_b_file   <- file.path(RDATA_DIR, "breeding_grm_vanraden.rds")
grm_n_file   <- file.path(RDATA_DIR, "natural_grm_vanraden.rds")

# Sample annotation tables: IID + Family (breeding) or Population (natural)
ann_b_file   <- file.path(RDATA_DIR, "breeding_sample_annotation.rds")
ann_n_file   <- file.path(RDATA_DIR, "natural_sample_annotation.rds")

# PLINK pairwise IBD estimates (long format: IID1, IID2, PI_HAT)
ibd_b_long   <- file.path(RDATA_DIR, "breeding_ibd_plink_long.rds")
ibd_n_long   <- file.path(RDATA_DIR, "natural_ibd_plink_long.rds")

# DAPC input objects: genotype matrix (X), group vector, SNP map (loc/chr/pos)
dapc_b_in    <- file.path(RDATA_DIR, "breeding_dapc_input_maf0.05_miss0.10.rds")
dapc_n_in    <- file.path(RDATA_DIR, "natural_dapc_input_maf0.05_miss0.10.rds")

# Verify all required files exist before running any analysis
ensure_file(grm_b_file); ensure_file(grm_n_file)
ensure_file(ann_b_file); ensure_file(ann_n_file)
ensure_file(ibd_b_long); ensure_file(ibd_n_long)
ensure_file(dapc_b_in);  ensure_file(dapc_n_in)

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
FS_LABEL <- 16  # panel label font size (a), b), c) …)

# ==============================================================================
# 5) HELPERS
# ==============================================================================

# ---- GRM heatmap with "nice" key labels (few decimals) ----
# Produces a ComplexHeatmap with hierarchical clustering and cohort-colour
# annotation bars on both axes. The legend shows evenly spaced "pretty" breaks
# rather than exact quantiles so tick labels are readable.
plot_grm_heatmap <- function(G, groups, palette_named, out_file, title_label = "a)",
                             dpi = 600, w_cm = 16, h_cm = 16,
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
    cluster_rows = TRUE,       # hierarchical clustering reveals population structure
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
         merge_legend = TRUE, show_annotation_legend = FALSE)
    # Panel label (e.g. "a)") placed in the top-left corner using grid coordinates
    grid.text(title_label, x = unit(0.35, "cm"),
              y = unit(1, "npc") - unit(0.35, "cm"),
              just = c("left","top"), gp = gpar(fontsize = FS_LABEL))
  }
  # Write TIFF at 600 dpi (publication quality) with lossless LZW compression
  tiff(out_file, width = w_cm, height = h_cm, units = "cm", res = dpi, compression = "lzw")
  draw_heatmap()
  dev.off()
  pdf(sub("\\.tiff$", ".pdf", out_file), width = w_cm / 2.54, height = h_cm / 2.54)
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
# and 0.50 (full-sib) help interpret relatedness levels.
plot_pihat_violin <- function(df, out_file, title_label = "b)",
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
              gp = gpar(fontsize = FS_LABEL))
  }
  tiff(out_file, width = w_cm, height = h_cm, units = "cm", res = dpi, compression = "lzw")
  draw_labeled()
  dev.off()
  pdf(sub("\\.tiff$", ".pdf", out_file), width = w_cm / 2.54, height = h_cm / 2.54)
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

# ---- DAPC runner ----
# Runs DAPC using adegenet. n.pca is capped to avoid retaining more PCs than
# samples allow; n.da is capped at 2 so biplots can always show DF1 vs DF2.
run_dapc_simple <- function(X, grp) {
  grp <- factor(grp)
  n.pca <- max(5L, min(50L, nrow(X) - nlevels(grp)))  # safe upper bound for PC retention
  n.da  <- max(1L, min(2L, nlevels(grp) - 1L))        # retain at most 2 discriminant axes
  fit <- adegenet::dapc(x = X, grp = grp, n.pca = n.pca, n.da = n.da, var.contrib = TRUE)
  # Reproject training samples onto the discriminant axes for the biplot
  pr <- try(adegenet::predict(fit, newdata = X)$x, silent = TRUE)
  coords <- if (inherits(pr, "try-error") || is.null(pr)) fit$ind.coord else pr
  coords <- as.data.frame(coords)
  colnames(coords) <- paste0("DF", seq_len(ncol(coords)))  # rename to DF1, DF2, …
  list(model = fit, coords = coords)
}

# ---- Extract ALL loadings for DF1 + DF2 and attach chromosome/position ----
# Returns a long-format table (one row per SNP × discriminant function)
# so readers can identify which loci drive each axis of separation.
get_all_loadings_df1_df2 <- function(dapc_model, snp_map) {
  vc <- as.data.frame(dapc_model$var.contr)  # SNP contribution scores from DAPC
  if (is.null(vc) || nrow(vc) == 0) stop("DAPC var.contr is empty.")
  vc$loc <- rownames(vc)

  # Identify the two numeric loading columns (DF1, DF2)
  num_cols <- which(vapply(vc, is.numeric, logical(1)))
  if (length(num_cols) < 2) stop("Could not find 2 numeric loading columns in dapc$var.contr.")
  c1 <- num_cols[1]; c2 <- num_cols[2]

  # Pivot to long format: one row per SNP × DF combination
  out <- vc %>%
    transmute(loc = loc, DF = "DF1", loading = .data[[names(vc)[c1]]]) %>%
    bind_rows(
      vc %>% transmute(loc = loc, DF = "DF2", loading = .data[[names(vc)[c2]]])
    )

  if (is.null(snp_map) || !all(c("loc","chr","pos") %in% colnames(snp_map))) {
    stop("SNP map missing loc/chr/pos. Re-run 8a with SNP map included in DAPC input RDS.")
  }

  # Join chromosome and position for each SNP so the table is directly informative
  out %>%
    left_join(snp_map %>% dplyr::select(loc, chr, pos), by = "loc") %>%
    relocate(chr, pos, DF, loading, .after = loc)
}

# ---- DAPC biplot with top-5 DF1 + top-5 DF2 loading arrows ----
# Plots individual sample scores on DF1 × DF2 coloured by group. Arrows show
# the top-5 SNPs by absolute loading on each axis, scaled so the longest arrow
# reaches 90% of the maximum individual radius (keeps arrows readable).
build_dapc_biplot <- function(dapc_res, snp_map, group_vec, palette_named, legend_title,
                              title_label = "c)") {
  coords <- dapc_res$coords
  coords$Group <- as.character(group_vec)

  # Extract loadings for DF1 and DF2 from the DAPC model
  vc <- as.data.frame(dapc_res$model$var.contr)
  vc$loc <- rownames(vc)
  num_cols <- which(vapply(vc, is.numeric, logical(1)))
  c1 <- num_cols[1]; c2 <- num_cols[2]

  tmp <- vc %>%
    transmute(loc = loc,
              DF1_load = .data[[names(vc)[c1]]],
              DF2_load = .data[[names(vc)[c2]]]) %>%
    left_join(snp_map %>% dplyr::select(loc, chr, pos), by = "loc")

  # Select the top 5 SNPs by absolute loading on DF1 and DF2 independently;
  # deduplicate in case a SNP ranks highly on both axes
  top1 <- tmp %>% arrange(desc(abs(DF1_load))) %>% slice_head(n = 5)
  top2 <- tmp %>% arrange(desc(abs(DF2_load))) %>% slice_head(n = 5)
  topm <- bind_rows(top1, top2) %>%
    distinct(loc, .keep_all = TRUE) %>%
    # Label as chr:pos if coordinates are available, otherwise use the SNP ID
    mutate(label = ifelse(!is.na(chr) & !is.na(pos), paste0(chr, ":", pos), loc))

  # Scale arrows so the longest arrow reaches 90% of the individual score radius,
  # keeping arrows visually proportional to the scatter of sample points
  ind_radius <- sqrt(coords$DF1^2 + coords$DF2^2)
  max_ind_r  <- max(ind_radius, na.rm = TRUE)
  load_radius <- sqrt(tmp$DF1_load^2 + tmp$DF2_load^2)
  max_load_r <- max(load_radius, na.rm = TRUE)
  scale_factor <- ifelse(is.finite(max_load_r) && max_load_r > 0, (0.9 * max_ind_r) / max_load_r, 1)

  topm <- topm %>%
    mutate(DF1_scaled = DF1_load * scale_factor,
           DF2_scaled = DF2_load * scale_factor)

  pal <- subset_palette(palette_named, coords$Group)

  ggplot(coords, aes(x = DF1, y = DF2, color = Group)) +
    geom_point(size = 1.8, alpha = 0.9) +
    scale_color_manual(values = pal) +
    labs(x = "DF1", y = "DF2", color = legend_title) +
    theme_minimal(base_size = FS_BASE) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = FS_BASE),
      legend.text  = element_text(size = FS_BASE - 1)
    ) +
    # Arrows from origin to scaled SNP loading position
    geom_segment(
      data = topm,
      aes(x = 0, y = 0, xend = DF1_scaled, yend = DF2_scaled),
      inherit.aes = FALSE,
      arrow = grid::arrow(length = unit(0.18, "cm")),
      color = "black",
      size = 0.6
    ) +
    # Non-overlapping SNP labels at arrow tips
    ggrepel::geom_text_repel(
      data = topm,
      aes(x = DF1_scaled, y = DF2_scaled, label = label),
      inherit.aes = FALSE,
      size = 3.2,
      min.segment.length = 0.05,
      box.padding = 0.25,
      point.padding = 0.2,
      segment.size = 0.3
    )
}

# Generic save helper: prints a ggplot object and overlays the panel label,
# then writes TIFF, PDF, EPS and PNG at the specified dimensions and resolution
save_tiff_with_label <- function(plot_obj, out_file, title_label,
                                 dpi = 600, w_cm = 16, h_cm = 16) {
  draw_obj <- function() {
    print(plot_obj)
    grid.text(title_label,
              x = unit(0.35, "cm"),
              y = unit(1, "npc") - unit(0.35, "cm"),
              just = c("left","top"),
              gp = gpar(fontsize = FS_LABEL))
  }
  tiff(out_file, width = w_cm, height = h_cm, units = "cm", res = dpi, compression = "lzw")
  draw_obj(); dev.off()
  pdf(sub("\\.tiff$", ".pdf", out_file), width = w_cm / 2.54, height = h_cm / 2.54)
  draw_obj(); dev.off()
  cairo_ps(sub("\\.tiff$", ".eps", out_file), width = w_cm / 2.54, height = h_cm / 2.54)
  draw_obj(); dev.off()
  png(sub("\\.tiff$", ".png", out_file), width = w_cm, height = h_cm, units = "cm", res = 150)
  draw_obj(); dev.off()
  invisible(out_file)
}

# ---- Combined 2×2 panel: a) GRM breeding, b) IBD breeding,
#                          c) GRM natural,  d) IBD natural ----
# Uses grid viewports to tile four sub-figures into a single 32×32 cm canvas.
save_fig2_panel <- function(ht_b, ht_n, p_ibd_b, p_ibd_n, out_file,
                             dpi = 600, w_cm = 32, h_cm = 32) {
  draw_panel <- function() {
    grid.newpage()

    # a) GRM breeding — top-left quadrant (x: 0–0.5, y: 0.5–1.0)
    pushViewport(viewport(x = 0, y = 0.5, width = 0.5, height = 0.5,
                          just = c("left", "bottom")))
    draw(ht_b, newpage = FALSE,
         heatmap_legend_side = "right", annotation_legend_side = "right",
         merge_legend = TRUE, show_annotation_legend = FALSE)
    grid.text("a)", x = unit(0.35, "cm"), y = unit(1, "npc") - unit(0.35, "cm"),
              just = c("left", "top"), gp = gpar(fontsize = FS_LABEL))
    popViewport()

    # b) IBD breeding — top-right quadrant
    pushViewport(viewport(x = 0.5, y = 0.5, width = 0.5, height = 0.5,
                          just = c("left", "bottom")))
    grid.draw(ggplotGrob(p_ibd_b))  # convert ggplot to grob for viewport drawing
    grid.text("b)", x = unit(0.35, "cm"), y = unit(1, "npc") - unit(0.35, "cm"),
              just = c("left", "top"), gp = gpar(fontsize = FS_LABEL))
    popViewport()

    # c) GRM natural — bottom-left quadrant
    pushViewport(viewport(x = 0, y = 0, width = 0.5, height = 0.5,
                          just = c("left", "bottom")))
    draw(ht_n, newpage = FALSE,
         heatmap_legend_side = "right", annotation_legend_side = "right",
         merge_legend = TRUE, show_annotation_legend = FALSE)
    grid.text("c)", x = unit(0.35, "cm"), y = unit(1, "npc") - unit(0.35, "cm"),
              just = c("left", "top"), gp = gpar(fontsize = FS_LABEL))
    popViewport()

    # d) IBD natural — bottom-right quadrant
    pushViewport(viewport(x = 0.5, y = 0, width = 0.5, height = 0.5,
                          just = c("left", "bottom")))
    grid.draw(ggplotGrob(p_ibd_n))
    grid.text("d)", x = unit(0.35, "cm"), y = unit(1, "npc") - unit(0.35, "cm"),
              just = c("left", "top"), gp = gpar(fontsize = FS_LABEL))
    popViewport()
  }

  tiff(out_file, width = w_cm, height = h_cm, units = "cm", res = dpi, compression = "lzw")
  draw_panel(); dev.off()
  pdf(sub("\\.tiff$", ".pdf", out_file), width = w_cm / 2.54, height = h_cm / 2.54)
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
# 7) FIGURE 1a / 2a : GRM heatmaps (16x16 cm), "nice" key labels, no group legend
# ==============================================================================
out_1a <- file.path(FIG_DIR, "Figure1a_breeding.tiff")
out_2a <- file.path(FIG_DIR, "Figure2a_natural.tiff")
res_1a <- plot_grm_heatmap(G_b, grp_b_vec, colors.17, out_1a, title_label = "a)", key_digits = 2)
res_2a <- plot_grm_heatmap(G_n, grp_n_vec, colors.25, out_2a, title_label = "a)", key_digits = 2)

# ==============================================================================
# 8) FIGURE 1b / 2b : IBD violins using PLINK PI_HAT (16x16 cm)
# ==============================================================================
# Classify each pairwise comparison as Within or Between family/population
pihat_b <- prep_pairs_pihat(ibd_b, ann_b, "Family")
pihat_n <- prep_pairs_pihat(ibd_n, ann_n, "Population")

out_1b <- file.path(FIG_DIR, "Figure1b_breeding.tiff")
out_2b <- file.path(FIG_DIR, "Figure2b_natural.tiff")

# y_limits include slight negative slack to prevent clipping of "Between" violin tails
res_1b <- plot_pihat_violin(pihat_b, out_1b, title_label = "b)", y_limits = c(-0.05, 0.70))
res_2b <- plot_pihat_violin(pihat_n, out_2b, title_label = "b)", y_limits = c(-0.05, 0.70))

# ==============================================================================
# 9) FIGURE 1c / 2c : DAPC biplots (16x16 cm) + ALL loadings tables (DF1 + DF2)
# ==============================================================================
b_in <- readRDS(dapc_b_in)  # list: X (genotype matrix), group (IID/Group), snp (loc/chr/pos)
n_in <- readRDS(dapc_n_in)

# Confirm SNP map is present before running DAPC — required for loading arrows
ensure_snp_map <- function(x) {
  if (is.null(x$snp) || !all(c("loc","chr","pos") %in% colnames(x$snp))) {
    stop("DAPC input RDS is missing SNP map with loc/chr/pos. Re-run 8a with SNP map included.")
  }
}
ensure_snp_map(b_in); ensure_snp_map(n_in)

# Align group labels to the row order of the genotype matrix
grp_b <- b_in$group$Group[match(rownames(b_in$X), b_in$group$IID)]
grp_b[is.na(grp_b)] <- "Unknown"
grp_n <- n_in$group$Group[match(rownames(n_in$X), n_in$group$IID)]
grp_n[is.na(grp_n)] <- "Unknown"

# Run DAPC; returns model object and DF1/DF2 coordinates for each sample
dapc_b <- run_dapc_simple(b_in$X, factor(grp_b))
dapc_n <- run_dapc_simple(n_in$X, factor(grp_n))

# Extract all SNP loadings for DF1 and DF2 (full tables, not just top SNPs)
ld_b_all <- get_all_loadings_df1_df2(dapc_b$model, b_in$snp)
ld_n_all <- get_all_loadings_df1_df2(dapc_n$model, n_in$snp)

# Write loading tables as TSV for supplementary data / further inspection
out_ld_b <- file.path(TAB_DAPC_LD, "Figure1c_breeding_DAPC_loadings_all_DF1_DF2.tsv")
out_ld_n <- file.path(TAB_DAPC_LD, "Figure2c_natural_DAPC_loadings_all_DF1_DF2.tsv")
write.table(ld_b_all, out_ld_b, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(ld_n_all, out_ld_n, sep = "\t", quote = FALSE, row.names = FALSE)

# Build biplot objects (ggplot) — arrows show top-5 loading SNPs per DF axis
p1c <- build_dapc_biplot(dapc_b, b_in$snp, grp_b, colors.17, legend_title = "Family",        title_label = "c)")
p2c <- build_dapc_biplot(dapc_n, n_in$snp, grp_n, colors.25, legend_title = "Natural Stand", title_label = "c)")

out_1c <- file.path(FIG_DIR, "Figure1c_breeding.tiff")
out_2c <- file.path(FIG_DIR, "Figure2c_natural.tiff")
save_tiff_with_label(p1c, out_1c, "c)", w_cm = 16, h_cm = 16)
save_tiff_with_label(p2c, out_2c, "c)", w_cm = 16, h_cm = 16)

# ==============================================================================
# 10) FIGURE 2: combined 2×2 panel (GRM + IBD, breeding + natural, 32×32 cm)
# ==============================================================================
out_fig2 <- file.path(FIG_DIR, "Figure2_PANEL_GRM_IBD_breeding_natural.tiff")
save_fig2_panel(res_1a$ht, res_2a$ht, res_1b$plot, res_2b$plot, out_fig2,
                dpi = 600, w_cm = 32, h_cm = 32)

cat(
  "Saved figures in:\n  ", FIG_DIR, "\n\n",
  "  - ", basename(out_1a), "\n",
  "  - ", basename(out_1b), "\n",
  "  - ", basename(out_1c), "\n",
  "  - ", basename(out_2a), "\n",
  "  - ", basename(out_2b), "\n",
  "  - ", basename(out_2c), "\n",
  "  - ", basename(out_fig2), "  [combined panel]\n\n",
  "Saved DAPC loadings tables in:\n  ", TAB_DAPC_LD, "\n\n",
  "  - ", basename(out_ld_b), "\n",
  "  - ", basename(out_ld_n), "\n",
  sep = ""
)

sessionInfo()

```

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

#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH -n 48
#SBATCH -N 1
#SBATCH --job-name=ECS.FQC1
#SBATCH --output=/path/to/your/project/LOGS/%x_%j.out
#SBATCH --error=/path/to/your/project/LOGS/%x_%j.err
#SBATCH --ntasks-per-socket 24
#SBATCH --time=48:00:00
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=YOUR_EMAIL

set -euo pipefail

echo "[$(date)] SLURM job started: ${SLURM_JOB_NAME:-no_slurm}"

module purge
module load fastqc/0.11.4
module load multiqc/1.27.1

# === USER CONFIGURATION ===
PROJECT_ROOT="/path/to/your/project"  # <-- set this
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

#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH -n 48
#SBATCH -N 1
#SBATCH --job-name=ECS.TRIM
#SBATCH --output=/path/to/your/project/LOGS/%x_%j.out
#SBATCH --error=/path/to/your/project/LOGS/%x_%j.err
#SBATCH --ntasks-per-socket 24
#SBATCH --time=48:00:00
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=YOUR_EMAIL

set -euo pipefail

echo "[$(date)] SLURM job started: ${SLURM_JOB_NAME:-no_slurm}"

module purge
module load fastqc/0.11.4
module load anaconda3/2020.11
module load trimmomatic/0.36

# === USER CONFIGURATION ===
PROJECT_ROOT="/path/to/your/project"  # <-- set this
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
# Hardcoded ECS sample list
list='P001_WA02 P001_WA03 P001_WA04 P001_WA05 P001_WA06 P001_WA07 P001_WA08 P001_WA09 P001_WA10 P001_WA11 P001_WB02 P001_WB03 P001_WB04 P001_WB05 P001_WB06 P001_WB07 P001_WB08 P001_WB09 P001_WB10 P001_WB11 P001_WB12 P001_WC01 P001_WC02 P001_WC03 P001_WC05 P001_WC06 P001_WC07 P001_WC09 P001_WC10 P001_WC11 P001_WC12 P001_WD01 P001_WD02 P001_WD03 P001_WD04 P001_WD05 P001_WD06 P001_WD07 P001_WD08 P001_WD09 P001_WD10 P001_WD11 P001_WD12 P001_WE01 P001_WE02 P001_WE03 P001_WE04 P001_WE05 P001_WE06 P001_WE07 P001_WE08 P001_WE09 P001_WE11 P001_WE12 P001_WF01 P001_WF02 P001_WF03 P001_WF04 P001_WF05 P001_WF06 P001_WF07 P001_WF08 P001_WF09 P001_WF10 P001_WF11 P001_WF12 P001_WG01 P001_WG03 P001_WG04 P001_WG05 P001_WG06 P001_WG07 P001_WG08 P001_WG09 P001_WG10 P001_WG11 P001_WG12 P001_WH01 P001_WH02 P001_WH03 P001_WH04 P001_WH05 P001_WH06 P001_WH07 P001_WH08 P001_WH09 P001_WH10 P001_WH11 P001_WH12 P002_WA01 P002_WA02 P002_WA03 P002_WA04 P002_WA05 P002_WA06 P002_WA07 P002_WA08 P002_WA09 P002_WA10 P002_WA11 P002_WA12 P002_WB01 P002_WB02 P002_WB03 P002_WB04 P002_WB05 P002_WB06 P002_WB07 P002_WB08 P002_WB09 P002_WB10 P002_WB11 P002_WB12 P002_WC01 P002_WC02 P002_WC03 P002_WC04 P002_WC05 P002_WC06 P002_WC07 P002_WC08 P002_WC09 P002_WC10 P002_WC11 P002_WC12 P002_WD01 P002_WD02 P002_WD03 P002_WD04 P002_WD05 P002_WD06 P002_WD07 P002_WD08 P002_WD09 P002_WD10 P002_WD11 P002_WD12 P002_WE01 P002_WE02 P002_WE03 P002_WE04 P002_WE05 P002_WE06 P002_WE07 P002_WE08 P002_WE09 P002_WE10 P002_WE11 P002_WE12 P002_WF01 P002_WF02 P002_WF03 P002_WF04 P002_WF05 P002_WF06 P002_WF07 P002_WF08 P002_WF09 P002_WF10 P002_WF11 P002_WF12 P002_WG01 P002_WG02 P002_WG03 P002_WG04 P002_WG05 P002_WG06 P002_WG07 P002_WG08 P002_WG09 P002_WG10 P002_WG11 P002_WG12 P002_WH01 P002_WH02 P002_WH03 P002_WH04 P002_WH05 P002_WH06 P002_WH07 P002_WH08 P002_WH09 P002_WH10 P002_WH11 P002_WH12 P003_WA01 P003_WA02 P003_WA03 P003_WA04 P003_WA05 P003_WA06 P003_WA07 P003_WA08 P003_WA09 P003_WA10 P003_WA11 P003_WA12 P003_WB01 P003_WB03 P003_WB04 P003_WB05 P003_WB06 P003_WB07 P003_WB08 P003_WB09 P003_WB10 P003_WB11 P003_WB12 P003_WC01 P004_WA01 P004_WA02 P004_WA03 P004_WA04 P004_WA05 P004_WA06 P004_WA07 P004_WA08 P004_WA09 P004_WA10 P004_WA11 P004_WA12 P004_WB01 P004_WB02 P004_WB03 P004_WB04 P004_WB05 P004_WB06 P004_WB07 P004_WB08 P004_WB09 P004_WB10 P004_WB11 P004_WB12 P004_WC01 P004_WC02 P004_WC03 P004_WC04 P004_WC05 P004_WC06 P004_WC07 P004_WC08 P004_WC09 P004_WC10 P004_WC11 P004_WC12 P004_WD01 P004_WD02 P004_WD03 P004_WD04 P004_WD05 P004_WD06 P004_WD07 P004_WD08 P004_WD09 P004_WD10 P004_WD11 P004_WD12 P004_WE01 P004_WE02 P004_WE03 P004_WE04 P004_WE05 P004_WE06 P004_WE07 P004_WE08 P004_WE09 P004_WE10 P004_WE11 P004_WE12 P004_WF01 P004_WF02 P004_WF03 P004_WF04 P004_WF05 P004_WF06 P004_WF07 P004_WF08 P004_WF09 P004_WF10 P004_WF11 P004_WF12 P004_WG01 P004_WG02 P004_WG03 P004_WG04 P004_WG05 P004_WG06 P004_WG07 P004_WG08 P004_WG09 P004_WG10 P004_WG11 P004_WG12 P004_WH01 P004_WH02 P004_WH03 P004_WH04 P004_WH05 P004_WH06 P004_WH07 P004_WH08 P004_WH09 P004_WH10 P004_WH11 P004_WH12 P005_WA01 P005_WA02 P005_WA03 P005_WA04 P005_WA05 P005_WA06 P005_WA07 P005_WA08 P005_WA09 P005_WA10 P005_WA11 P005_WA12 P005_WB01 P005_WB02 P005_WB03 P005_WB04 P005_WB05 P005_WB06 P005_WB07 P005_WB08 P005_WB09 P005_WB10 P005_WB11 P005_WB12 P005_WC01 P005_WC02 P005_WC03 P005_WC04 P005_WC05 P005_WC06 P005_WC07 P005_WC08 P005_WC09 P005_WC10 P005_WC11 P005_WC12 P005_WD01 P005_WD02 P005_WD03 P005_WD04 P005_WD05 P005_WD06 P005_WD07 P005_WD08 P005_WD09 P005_WD10 P005_WD11 P005_WD12 P005_WE01 P005_WE02 P005_WE03 P005_WE04 P005_WE05 P005_WE06 P005_WE07 P005_WE08 P005_WE09 P005_WE10 P005_WE11 P005_WE12 P005_WF01 P005_WF02 P005_WF03 P005_WF04 P005_WF05 P005_WF06 P005_WF07 P005_WF08 P005_WF09 P005_WF10 P005_WF11 P005_WF12 P005_WG01 P005_WG02 P005_WG03 P005_WG04 P005_WG05 P005_WG06 P005_WG07 P005_WG08 P005_WG09 P005_WG10 P005_WG11 P005_WG12 P005_WH01 P005_WH02 P005_WH03 P005_WH04 P005_WH05 P005_WH06 P005_WH07 P005_WH08 P005_WH09 P005_WH10 P005_WH11 P005_WH12 P006_WA01 P006_WA02 P006_WA03 P006_WA04 P006_WA05 P006_WA06 P006_WA07 P006_WA08 P006_WA09 P006_WA10 P006_WA11 P006_WB01 P006_WB02 P006_WB03 P006_WB04 P006_WB05 P006_WB06 P006_WB07 P006_WB08 P006_WB09 P006_WB10 P006_WB11 P006_WB12 P006_WC01 P006_WC02 P006_WC03 P006_WC04 P006_WC05 P006_WC06 P006_WC07 P006_WC08 P006_WC09 P006_WC10 P006_WC12 P006_WD01 P006_WD02 P006_WD03 P006_WD04 P006_WD05 P006_WD06 P006_WD07 P006_WD08 P006_WD09 P006_WD10 P006_WD11 P006_WD12 P006_WE01 P006_WE02 P006_WE03 P006_WE04 P006_WE05 P006_WE06 P006_WE07 P006_WE08 P006_WE09 P006_WE10 P006_WE11 P006_WE12 P006_WF01 P006_WF02 P006_WF03 P006_WF04 P006_WF05 P006_WF06 P006_WF07 P006_WF08 P006_WF09 P006_WF10 P006_WF11 P006_WF12 P006_WG01 P006_WG02 P006_WG03 P006_WG04 P006_WG05 P006_WG06 P006_WG07 P006_WG08 P006_WG09 P006_WG10 P006_WG11 P006_WG12 P006_WH01 P006_WH02 P006_WH03 P006_WH04 P006_WH05 P006_WH06 P006_WH07 P006_WH08 P006_WH09 P006_WH10 P006_WH11 P006_WH12 P007_WA01 P007_WA03 P007_WA04 P007_WA05 P007_WA08 P007_WA09 P007_WA10 P007_WA11 P007_WA12 P007_WB01 P007_WB02 P007_WB03 P007_WB04 P007_WB05 P007_WB06 P007_WB07 P007_WB08 P007_WB09 P007_WB10 P007_WB11 P007_WB12 P007_WC01 P007_WC02 P007_WC03 P007_WC04 P007_WC05 P007_WC06 P007_WC07 P007_WC08 P007_WC09 P007_WC10 P007_WC11 P007_WC12 P007_WD01 P007_WD02 P007_WD03 P007_WD04 P007_WD05 P007_WD06 P007_WD09 P007_WD10 P007_WD11 P007_WD12 P007_WE01 P007_WE02 P007_WE03 P007_WE04 P007_WE05 P007_WE06 P007_WE07 P007_WE08 P007_WE09 P007_WE10 P007_WE11 P007_WE12 P007_WF01 P007_WF02 P007_WF03 P007_WF04 P007_WF05 P007_WF06 P007_WF07 P007_WF08 P007_WF09 P007_WF10 P007_WF11 P007_WF12 P007_WG01 P007_WG02 P007_WG03 P007_WG04 P007_WG05 P007_WG06 P007_WG07 P007_WG08 P007_WG09 P007_WG10 P007_WG11 P007_WG12 P007_WH01 P007_WH02 P007_WH04 P007_WH05 P007_WH06 P007_WH07 P007_WH08 P007_WH09 P007_WH10 P007_WH11 P007_WH12 P008_WA01 P008_WA02 P008_WA03 P008_WA04 P008_WA05 P008_WA06 P008_WA07 P008_WA08 P008_WA09 P008_WA10 P008_WA11 P008_WA12 P008_WB02 P008_WB03 P008_WB04 P008_WB05 P008_WB06'

for sample in ${list}; do
  r1="${INPUT}/${sample}_R1.fastq.gz"
  r2="${INPUT}/${sample}_R2.fastq.gz"

  if [[ ! -s "${r1}" || ! -s "${r2}" ]]; then
    echo "WARNING: missing pair for ${sample}; skipping"
    continue
  fi

  java -jar /usr/product/bioinfo/SL_7.0/BIOINFORMATICS/TRIMMOMATIC/0.36/trimmomatic-0.36.jar PE \
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

#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH -n 96
#SBATCH -N 1
#SBATCH --job-name=ECS.BT2
#SBATCH --output=/path/to/your/project/LOGS/%x_%j.out
#SBATCH --error=/path/to/your/project/LOGS/%x_%j.err
#SBATCH --ntasks-per-socket 24
#SBATCH --time=48:00:00
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=YOUR_EMAIL

set -euo pipefail

echo "[$(date)] SLURM job started: ${SLURM_JOB_NAME:-no_slurm}"

module purge
module load gcc/14.2.0
module load bowtie2/2.5.4
module load samtools/1.21

# === USER CONFIGURATION ===
PROJECT_ROOT="/path/to/your/project"  # <-- set this
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
# Hardcoded ECS sample list
list='P001_WA02 P001_WA03 P001_WA04 P001_WA05 P001_WA06 P001_WA07 P001_WA08 P001_WA09 P001_WA10 P001_WA11 P001_WB02 P001_WB03 P001_WB04 P001_WB05 P001_WB06 P001_WB07 P001_WB08 P001_WB09 P001_WB10 P001_WB11 P001_WB12 P001_WC01 P001_WC02 P001_WC03 P001_WC05 P001_WC06 P001_WC07 P001_WC09 P001_WC10 P001_WC11 P001_WC12 P001_WD01 P001_WD02 P001_WD03 P001_WD04 P001_WD05 P001_WD06 P001_WD07 P001_WD08 P001_WD09 P001_WD10 P001_WD11 P001_WD12 P001_WE01 P001_WE02 P001_WE03 P001_WE04 P001_WE05 P001_WE06 P001_WE07 P001_WE08 P001_WE09 P001_WE11 P001_WE12 P001_WF01 P001_WF02 P001_WF03 P001_WF04 P001_WF05 P001_WF06 P001_WF07 P001_WF08 P001_WF09 P001_WF10 P001_WF11 P001_WF12 P001_WG01 P001_WG03 P001_WG04 P001_WG05 P001_WG06 P001_WG07 P001_WG08 P001_WG09 P001_WG10 P001_WG11 P001_WG12 P001_WH01 P001_WH02 P001_WH03 P001_WH04 P001_WH05 P001_WH06 P001_WH07 P001_WH08 P001_WH09 P001_WH10 P001_WH11 P001_WH12 P002_WA01 P002_WA02 P002_WA03 P002_WA04 P002_WA05 P002_WA06 P002_WA07 P002_WA08 P002_WA09 P002_WA10 P002_WA11 P002_WA12 P002_WB01 P002_WB02 P002_WB03 P002_WB04 P002_WB05 P002_WB06 P002_WB07 P002_WB08 P002_WB09 P002_WB10 P002_WB11 P002_WB12 P002_WC01 P002_WC02 P002_WC03 P002_WC04 P002_WC05 P002_WC06 P002_WC07 P002_WC08 P002_WC09 P002_WC10 P002_WC11 P002_WC12 P002_WD01 P002_WD02 P002_WD03 P002_WD04 P002_WD05 P002_WD06 P002_WD07 P002_WD08 P002_WD09 P002_WD10 P002_WD11 P002_WD12 P002_WE01 P002_WE02 P002_WE03 P002_WE04 P002_WE05 P002_WE06 P002_WE07 P002_WE08 P002_WE09 P002_WE10 P002_WE11 P002_WE12 P002_WF01 P002_WF02 P002_WF03 P002_WF04 P002_WF05 P002_WF06 P002_WF07 P002_WF08 P002_WF09 P002_WF10 P002_WF11 P002_WF12 P002_WG01 P002_WG02 P002_WG03 P002_WG04 P002_WG05 P002_WG06 P002_WG07 P002_WG08 P002_WG09 P002_WG10 P002_WG11 P002_WG12 P002_WH01 P002_WH02 P002_WH03 P002_WH04 P002_WH05 P002_WH06 P002_WH07 P002_WH08 P002_WH09 P002_WH10 P002_WH11 P002_WH12 P003_WA01 P003_WA02 P003_WA03 P003_WA04 P003_WA05 P003_WA06 P003_WA07 P003_WA08 P003_WA09 P003_WA10 P003_WA11 P003_WA12 P003_WB01 P003_WB03 P003_WB04 P003_WB05 P003_WB06 P003_WB07 P003_WB08 P003_WB09 P003_WB10 P003_WB11 P003_WB12 P003_WC01 P004_WA01 P004_WA02 P004_WA03 P004_WA04 P004_WA05 P004_WA06 P004_WA07 P004_WA08 P004_WA09 P004_WA10 P004_WA11 P004_WA12 P004_WB01 P004_WB02 P004_WB03 P004_WB04 P004_WB05 P004_WB06 P004_WB07 P004_WB08 P004_WB09 P004_WB10 P004_WB11 P004_WB12 P004_WC01 P004_WC02 P004_WC03 P004_WC04 P004_WC05 P004_WC06 P004_WC07 P004_WC08 P004_WC09 P004_WC10 P004_WC11 P004_WC12 P004_WD01 P004_WD02 P004_WD03 P004_WD04 P004_WD05 P004_WD06 P004_WD07 P004_WD08 P004_WD09 P004_WD10 P004_WD11 P004_WD12 P004_WE01 P004_WE02 P004_WE03 P004_WE04 P004_WE05 P004_WE06 P004_WE07 P004_WE08 P004_WE09 P004_WE10 P004_WE11 P004_WE12 P004_WF01 P004_WF02 P004_WF03 P004_WF04 P004_WF05 P004_WF06 P004_WF07 P004_WF08 P004_WF09 P004_WF10 P004_WF11 P004_WF12 P004_WG01 P004_WG02 P004_WG03 P004_WG04 P004_WG05 P004_WG06 P004_WG07 P004_WG08 P004_WG09 P004_WG10 P004_WG11 P004_WG12 P004_WH01 P004_WH02 P004_WH03 P004_WH04 P004_WH05 P004_WH06 P004_WH07 P004_WH08 P004_WH09 P004_WH10 P004_WH11 P004_WH12 P005_WA01 P005_WA02 P005_WA03 P005_WA04 P005_WA05 P005_WA06 P005_WA07 P005_WA08 P005_WA09 P005_WA10 P005_WA11 P005_WA12 P005_WB01 P005_WB02 P005_WB03 P005_WB04 P005_WB05 P005_WB06 P005_WB07 P005_WB08 P005_WB09 P005_WB10 P005_WB11 P005_WB12 P005_WC01 P005_WC02 P005_WC03 P005_WC04 P005_WC05 P005_WC06 P005_WC07 P005_WC08 P005_WC09 P005_WC10 P005_WC11 P005_WC12 P005_WD01 P005_WD02 P005_WD03 P005_WD04 P005_WD05 P005_WD06 P005_WD07 P005_WD08 P005_WD09 P005_WD10 P005_WD11 P005_WD12 P005_WE01 P005_WE02 P005_WE03 P005_WE04 P005_WE05 P005_WE06 P005_WE07 P005_WE08 P005_WE09 P005_WE10 P005_WE11 P005_WE12 P005_WF01 P005_WF02 P005_WF03 P005_WF04 P005_WF05 P005_WF06 P005_WF07 P005_WF08 P005_WF09 P005_WF10 P005_WF11 P005_WF12 P005_WG01 P005_WG02 P005_WG03 P005_WG04 P005_WG05 P005_WG06 P005_WG07 P005_WG08 P005_WG09 P005_WG10 P005_WG11 P005_WG12 P005_WH01 P005_WH02 P005_WH03 P005_WH04 P005_WH05 P005_WH06 P005_WH07 P005_WH08 P005_WH09 P005_WH10 P005_WH11 P005_WH12 P006_WA01 P006_WA02 P006_WA03 P006_WA04 P006_WA05 P006_WA06 P006_WA07 P006_WA08 P006_WA09 P006_WA10 P006_WA11 P006_WB01 P006_WB02 P006_WB03 P006_WB04 P006_WB05 P006_WB06 P006_WB07 P006_WB08 P006_WB09 P006_WB10 P006_WB11 P006_WB12 P006_WC01 P006_WC02 P006_WC03 P006_WC04 P006_WC05 P006_WC06 P006_WC07 P006_WC08 P006_WC09 P006_WC10 P006_WC12 P006_WD01 P006_WD02 P006_WD03 P006_WD04 P006_WD05 P006_WD06 P006_WD07 P006_WD08 P006_WD09 P006_WD10 P006_WD11 P006_WD12 P006_WE01 P006_WE02 P006_WE03 P006_WE04 P006_WE05 P006_WE06 P006_WE07 P006_WE08 P006_WE09 P006_WE10 P006_WE11 P006_WE12 P006_WF01 P006_WF02 P006_WF03 P006_WF04 P006_WF05 P006_WF06 P006_WF07 P006_WF08 P006_WF09 P006_WF10 P006_WF11 P006_WF12 P006_WG01 P006_WG02 P006_WG03 P006_WG04 P006_WG05 P006_WG06 P006_WG07 P006_WG08 P006_WG09 P006_WG10 P006_WG11 P006_WG12 P006_WH01 P006_WH02 P006_WH03 P006_WH04 P006_WH05 P006_WH06 P006_WH07 P006_WH08 P006_WH09 P006_WH10 P006_WH11 P006_WH12 P007_WA01 P007_WA03 P007_WA04 P007_WA05 P007_WA08 P007_WA09 P007_WA10 P007_WA11 P007_WA12 P007_WB01 P007_WB02 P007_WB03 P007_WB04 P007_WB05 P007_WB06 P007_WB07 P007_WB08 P007_WB09 P007_WB10 P007_WB11 P007_WB12 P007_WC01 P007_WC02 P007_WC03 P007_WC04 P007_WC05 P007_WC06 P007_WC07 P007_WC08 P007_WC09 P007_WC10 P007_WC11 P007_WC12 P007_WD01 P007_WD02 P007_WD03 P007_WD04 P007_WD05 P007_WD06 P007_WD09 P007_WD10 P007_WD11 P007_WD12 P007_WE01 P007_WE02 P007_WE03 P007_WE04 P007_WE05 P007_WE06 P007_WE07 P007_WE08 P007_WE09 P007_WE10 P007_WE11 P007_WE12 P007_WF01 P007_WF02 P007_WF03 P007_WF04 P007_WF05 P007_WF06 P007_WF07 P007_WF08 P007_WF09 P007_WF10 P007_WF11 P007_WF12 P007_WG01 P007_WG02 P007_WG03 P007_WG04 P007_WG05 P007_WG06 P007_WG07 P007_WG08 P007_WG09 P007_WG10 P007_WG11 P007_WG12 P007_WH01 P007_WH02 P007_WH04 P007_WH05 P007_WH06 P007_WH07 P007_WH08 P007_WH09 P007_WH10 P007_WH11 P007_WH12 P008_WA01 P008_WA02 P008_WA03 P008_WA04 P008_WA05 P008_WA06 P008_WA07 P008_WA08 P008_WA09 P008_WA10 P008_WA11 P008_WA12 P008_WB02 P008_WB03 P008_WB04 P008_WB05 P008_WB06'

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

#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --job-name=ECS.CALL
#SBATCH --output=/path/to/your/project/LOGS/%x_%j.out
#SBATCH --error=/path/to/your/project/LOGS/%x_%j.err
#SBATCH --time=48:00:00
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=YOUR_EMAIL

set -euo pipefail

echo "[$(date)] SLURM job started: ${SLURM_JOB_NAME:-no_slurm}"

module purge
module load gcc/14.2.0
module load bcftools/1.19
module load samtools/1.21

# === USER CONFIGURATION ===
PROJECT_ROOT="/path/to/your/project"  # <-- set this
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

### `4a2.tgc.ecs.snv.filtering.bcftools.vcftools.sh`

```bash
#!/bin/bash
#-------------------------------------------------------------------------------
# TreeGeneClimate (TGC) — ECS
# Step 4a2: SNV filtering + basic stats (bcftools) from merged unfiltered VCF
#
# Project root:
#   /path/to/your/project
#
# Input (results):
#   RESULTS/ECS/VARIANT.CALLING/tgc.ecs.allsamples.call.unfilt.snvs.vcf.gz
#
# Output (results):
#   RESULTS/ECS/VARIANT.CALLING/tgc.ecs.allsamples.call.filt.maf10.snvs.vcf.gz
#   RESULTS/ECS/VARIANT.CALLING/tgc.ecs.allsamples.call.filt.maf01.snvs.vcf.gz
#   RESULTS/ECS/VARIANT.CALLING/STATS_ALL/
#-------------------------------------------------------------------------------

#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH -n 96
#SBATCH -N 1
#SBATCH --job-name=ECS.FILT
#SBATCH --output=/path/to/your/project/LOGS/%x_%j.out
#SBATCH --error=/path/to/your/project/LOGS/%x_%j.err
#SBATCH --time=48:00:00
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=YOUR_EMAIL

set -euo pipefail

echo "[$(date)] SLURM job started: ${SLURM_JOB_NAME:-no_slurm}"

module purge
module load gcc/14.2.0
module load bcftools/1.19
module load samtools/1.21

# === USER CONFIGURATION ===
PROJECT_ROOT="/path/to/your/project"  # <-- set this
# ===========================

VCF_BASE="${PROJECT_ROOT}/RESULTS/ECS/VARIANT.CALLING"
STATS_DIR="${VCF_BASE}/STATS_ALL"
LOGS_DIR="${PROJECT_ROOT}/LOGS"

UNFILT_VCF="${VCF_BASE}/tgc.ecs.allsamples.call.unfilt.snvs.renamed.vcf.gz"
FILT_MAF10_VCF="${VCF_BASE}/tgc.ecs.allsamples.call.filt.maf10.snvs.vcf.gz"
FILT_MAF01_VCF="${VCF_BASE}/tgc.ecs.allsamples.call.filt.maf01.snvs.vcf.gz"

mkdir -p "${STATS_DIR}" "${LOGS_DIR}"

if [[ ! -s "${UNFILT_VCF}" ]]; then
  echo "ERROR: Input VCF not found: ${UNFILT_VCF}"
  exit 1
fi

echo "[$(date)] Counting unfiltered SNVs..."
bcftools view -H "${UNFILT_VCF}" | wc -l

# ------------------------------------------------------------------------------
echo "[$(date)] Filtering (MAF >= 0.10)..."
bcftools +fill-tags "${UNFILT_VCF}" -Ou -- -t MAF,F_MISSING,AF \
| bcftools view \
    --types snps \
    --min-alleles 2 \
    --max-alleles 2 \
    --include 'F_MISSING<=0.2 && MAF>=0.10 && INFO/DP>=10' \
    --threads 48 \
    -Oz -o "${FILT_MAF10_VCF}"

# CSI index (safe for large contigs/positions)
bcftools index -c "${FILT_MAF10_VCF}"

echo "[$(date)] Filtering (MAF >= 0.01)..."
bcftools +fill-tags "${UNFILT_VCF}" -Ou -- -t MAF,F_MISSING,AF \
| bcftools view \
    --types snps \
    --min-alleles 2 \
    --max-alleles 2 \
    --include 'F_MISSING<=0.2 && MAF>=0.01 && INFO/DP>=10' \
    --threads 48 \
    -Oz -o "${FILT_MAF01_VCF}"

# CSI index (safe for large contigs/positions)
bcftools index -c "${FILT_MAF01_VCF}"

echo "[$(date)] Counting filtered SNVs (MAF >= 0.10)..."
bcftools view -H "${FILT_MAF10_VCF}" | wc -l

echo "[$(date)] Counting filtered SNVs (MAF >= 0.01)..."
bcftools view -H "${FILT_MAF01_VCF}" | wc -l

# ------------------------------------------------------------------------------
echo "[$(date)] Writing stats tables (unfiltered vs MAF01 filtered)..."

bcftools +fill-tags "${UNFILT_VCF}" -- -t AF,MAF,F_MISSING \
| bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\t%MAF\t%F_MISSING\n' \
> "${STATS_DIR}/unfiltered.af_maf_missing.txt"

bcftools +fill-tags "${FILT_MAF01_VCF}" -- -t AF,MAF,F_MISSING \
| bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\t%MAF\t%F_MISSING\n' \
> "${STATS_DIR}/filtered_maf01.af_maf_missing.txt"

bcftools query -f '%CHROM\t%POS\t%QUAL\n' "${UNFILT_VCF}" > "${STATS_DIR}/unfiltered.qual.txt"
bcftools query -f '%CHROM\t%POS\t%QUAL\n' "${FILT_MAF01_VCF}" > "${STATS_DIR}/filtered_maf01.qual.txt"

bcftools query -f '%CHROM\t%POS\t%INFO/DP\n' "${UNFILT_VCF}" > "${STATS_DIR}/unfiltered.site_dp.txt"
bcftools query -f '%CHROM\t%POS\t%INFO/DP\n' "${FILT_MAF1_VCF}" > "${STATS_DIR}/filtered_maf01.site_dp.txt"

bcftools stats -s - "${UNFILT_VCF}" > "${STATS_DIR}/unfiltered.bcftools.stats.txt"
bcftools stats -s - "${FILT_MAF01_VCF}" > "${STATS_DIR}/filtered_maf01.bcftools.stats.txt"

echo "[$(date)] Outputs:"
echo "  Unfiltered: ${UNFILT_VCF}"
echo "  Filtered MAF10: ${FILT_MAF10_VCF}"
echo "  Filtered MAF01: ${FILT_MAF01_VCF}"
echo "  Stats: ${STATS_DIR}"

echo "[$(date)] Done."
exit 0

```

---

### `5a.tgc.ecs.plink.admixture.split.gwasprep.sh`

```bash
#!/bin/bash
#-------------------------------------------------------------------------------
# TreeGeneClimate (TGC) — ECS
# Step 5a: Split final filtered VCF into BREEDING vs NATURAL cohorts, then:
#          (i) GWAS input BED sets (unpruned)
#          (ii) LD pruning (for structure analyses)
#          (iii) PCA on pruned SNPs
#          (iv) ADMIXTURE on pruned SNPs
#
# Project root:
#   /path/to/your/project
#
# Input (results):
#   RESULTS/ECS/VARIANT.CALLING/tgc.ecs.allsamples.call.filt.maf05.snvs.vcf.gz
#
# Cohort sample lists (frozen metadata):
#   DATA/METADATA/ECS/all_breeding_samples.txt
#   DATA/METADATA/ECS/all_natural_samples.txt
#
# Output (results):
#   RESULTS/ECS/VCF_SPLIT/                     (cohort VCFs)
#   RESULTS/ECS/GWAS/PLINK_INPUT/             (unpruned BED sets)
#   RESULTS/ECS/POPGEN/STRUCTURE/             (prune + PCA + ADMIXTURE inputs/outputs)
#-------------------------------------------------------------------------------

#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --job-name=ECS.PLINK.ADMIX
#SBATCH --output=/path/to/your/project/LOGS/%x_%j.out
#SBATCH --error=/path/to/your/project/LOGS/%x_%j.err
#SBATCH --time=48:00:00
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=YOUR_EMAIL

set -euo pipefail

echo "[$(date)] SLURM job started: ${SLURM_JOB_NAME:-no_slurm}"

# ------------------------------------------------------------------------------
# Load modules
module purge
module load gcc/14.2.0
module load bcftools/1.19
module load plink/1.9
module load miniforge3/24.3.0-0
source activate admixture

# ------------------------------------------------------------------------------
# Paths
# === USER CONFIGURATION ===
PROJECT_ROOT="/path/to/your/project"  # <-- set this
# ===========================
LOGS="${PROJECT_ROOT}/LOGS"

VCF_CALL_DIR="${PROJECT_ROOT}/RESULTS/ECS/VARIANT.CALLING"
VCF_ALL="${VCF_CALL_DIR}/tgc.ecs.allsamples.call.filt.maf05.snvs.vcf.gz"

META_DIR="${PROJECT_ROOT}/DATA/METADATA"
BREEDING_SAMPLES="${META_DIR}/all_breeding_samples.txt"
NATURAL_SAMPLES="${META_DIR}/all_natural_samples.txt"

SPLIT_DIR="${PROJECT_ROOT}/RESULTS/ECS/VCF_SPLIT"
VCF_BREEDING="${SPLIT_DIR}/tgc.ecs.breeding.call.filt.maf05.snvs.vcf.gz"
VCF_NATURAL="${SPLIT_DIR}/tgc.ecs.natural.call.filt.maf05.snvs.vcf.gz"

GWAS_DIR="${PROJECT_ROOT}/RESULTS/ECS/GWAS/PLINK_INPUT"
GWAS_BREEDING_DIR="${GWAS_DIR}/BREEDING"
GWAS_NATURAL_DIR="${GWAS_DIR}/NATURAL"

STRUCT_DIR="${PROJECT_ROOT}/RESULTS/ECS/POPGEN/STRUCTURE"
BREEDING_DIR="${STRUCT_DIR}/BREEDING"
NATURAL_DIR="${STRUCT_DIR}/NATURAL"

BREEDING_PCA_DIR="${BREEDING_DIR}/PCA"
NATURAL_PCA_DIR="${NATURAL_DIR}/PCA"

BREEDING_ADMIX_DIR="${BREEDING_DIR}/ADMIXTURE"
NATURAL_ADMIX_DIR="${NATURAL_DIR}/ADMIXTURE"

mkdir -p \
  "${LOGS}" \
  "${SPLIT_DIR}" \
  "${GWAS_BREEDING_DIR}" "${GWAS_NATURAL_DIR}" \
  "${BREEDING_PCA_DIR}" "${NATURAL_PCA_DIR}" \
  "${BREEDING_ADMIX_DIR}" "${NATURAL_ADMIX_DIR}"

# ------------------------------------------------------------------------------
# Sanity checks (fail fast)
if [[ ! -s "${VCF_ALL}" ]]; then
  echo "ERROR: missing input VCF: ${VCF_ALL}"
  exit 1
fi
if [[ ! -s "${BREEDING_SAMPLES}" ]]; then
  echo "ERROR: missing breeding sample list: ${BREEDING_SAMPLES}"
  exit 1
fi
if [[ ! -s "${NATURAL_SAMPLES}" ]]; then
  echo "ERROR: missing natural sample list: ${NATURAL_SAMPLES}"
  exit 1
fi

# ------------------------------------------------------------------------------
# Index VCF (CSI; robust for large coordinates)
if [[ ! -s "${VCF_ALL}.csi" ]]; then
  echo "[$(date)] Indexing input VCF (CSI)"
  bcftools index -c --threads 48 "${VCF_ALL}"
fi

# ------------------------------------------------------------------------------
# Split cohorts
echo "[$(date)] Splitting cohorts"
bcftools view --threads 48 -S "${BREEDING_SAMPLES}" -Oz -o "${VCF_BREEDING}" "${VCF_ALL}"
bcftools view --threads 48 -S "${NATURAL_SAMPLES}"  -Oz -o "${VCF_NATURAL}"  "${VCF_ALL}"

bcftools index -c --threads 48 "${VCF_BREEDING}"
bcftools index -c --threads 48 "${VCF_NATURAL}"

# ------------------------------------------------------------------------------
# GWAS input BED sets (UNPRUNED)
echo "[$(date)] Creating GWAS BED sets (unpruned)"
plink --vcf "${VCF_BREEDING}" --double-id --allow-extra-chr \
  --set-missing-var-ids @:# \
  --threads 48 \
  --make-bed --out "${GWAS_BREEDING_DIR}/tgc.ecs.breeding.gwas.unpruned"

plink --vcf "${VCF_NATURAL}" --double-id --allow-extra-chr \
  --set-missing-var-ids @:# \
  --threads 48 \
  --make-bed --out "${GWAS_NATURAL_DIR}/tgc.ecs.natural.gwas.unpruned"

# ------------------------------------------------------------------------------
# LD pruning (STRUCTURE analyses)
echo "[$(date)] LD pruning (STRUCTURE analyses)"
plink --vcf "${VCF_BREEDING}" --double-id --allow-extra-chr \
  --set-missing-var-ids @:# \
  --threads 48 \
  --indep-pairwise 50 10 0.2 \
  --out "${BREEDING_DIR}/tgc.ecs.breeding.prune"

plink --vcf "${VCF_NATURAL}" --double-id --allow-extra-chr \
  --set-missing-var-ids @:# \
  --threads 48 \
  --indep-pairwise 50 10 0.2 \
  --out "${NATURAL_DIR}/tgc.ecs.natural.prune"

# ------------------------------------------------------------------------------
# PCA on PRUNED SNP sets
echo "[$(date)] PCA on pruned SNP sets"
plink --vcf "${VCF_BREEDING}" --double-id --allow-extra-chr \
  --set-missing-var-ids @:# \
  --threads 48 \
  --extract "${BREEDING_DIR}/tgc.ecs.breeding.prune.prune.in" \
  --pca --out "${BREEDING_PCA_DIR}/tgc.ecs.breeding.pca.pruned"

plink --vcf "${VCF_NATURAL}" --double-id --allow-extra-chr \
  --set-missing-var-ids @:# \
  --threads 48 \
  --extract "${NATURAL_DIR}/tgc.ecs.natural.prune.prune.in" \
  --pca --out "${NATURAL_PCA_DIR}/tgc.ecs.natural.pca.pruned"

# ------------------------------------------------------------------------------
# ADMIXTURE inputs (PRUNED SNP sets) + runs
echo "[$(date)] Preparing ADMIXTURE inputs (pruned SNP sets)"
plink --vcf "${VCF_BREEDING}" --double-id --allow-extra-chr \
  --set-missing-var-ids @:# \
  --threads 48 \
  --extract "${BREEDING_DIR}/tgc.ecs.breeding.prune.prune.in" \
  --make-bed --out "${BREEDING_DIR}/tgc.ecs.breeding.admix.pruned"

plink --vcf "${VCF_NATURAL}" --double-id --allow-extra-chr \
  --set-missing-var-ids @:# \
  --threads 48 \
  --extract "${NATURAL_DIR}/tgc.ecs.natural.prune.prune.in" \
  --make-bed --out "${NATURAL_DIR}/tgc.ecs.natural.admix.pruned"

# ADMIXTURE expects numeric family IDs; set FID=0 in .fam
awk 'BEGIN{OFS="\t"}{$1="0"; print}' "${BREEDING_DIR}/tgc.ecs.breeding.admix.pruned.fam" > "${BREEDING_DIR}/tgc.ecs.breeding.admix.pruned.fam.tmp" \
  && mv "${BREEDING_DIR}/tgc.ecs.breeding.admix.pruned.fam.tmp" "${BREEDING_DIR}/tgc.ecs.breeding.admix.pruned.fam"

awk 'BEGIN{OFS="\t"}{$1="0"; print}' "${NATURAL_DIR}/tgc.ecs.natural.admix.pruned.fam" > "${NATURAL_DIR}/tgc.ecs.natural.admix.pruned.fam.tmp" \
  && mv "${NATURAL_DIR}/tgc.ecs.natural.admix.pruned.fam.tmp" "${NATURAL_DIR}/tgc.ecs.natural.admix.pruned.fam"

# Breeding: K=1..20
echo "[$(date)] ADMIXTURE (breeding) K=1..20"
for K in $(seq 1 20); do
  admixture -j96 --seed=123 "${BREEDING_DIR}/tgc.ecs.breeding.admix.pruned.bed" "${K}" \
    > "${BREEDING_ADMIX_DIR}/K${K}.out" 2> "${BREEDING_ADMIX_DIR}/K${K}.err"
done
grep -h "CV error" "${BREEDING_ADMIX_DIR}"/K*.out > "${BREEDING_ADMIX_DIR}/breeding_cv_errors.txt" || true

# Natural: K=1..10
echo "[$(date)] ADMIXTURE (natural) K=1..10"
for K in $(seq 1 10); do
  admixture -j96 --seed=123 "${NATURAL_DIR}/tgc.ecs.natural.admix.pruned.bed" "${K}" \
    > "${NATURAL_ADMIX_DIR}/K${K}.out" 2> "${NATURAL_ADMIX_DIR}/K${K}.err"
done
grep -h "CV error" "${NATURAL_ADMIX_DIR}"/K*.out > "${NATURAL_ADMIX_DIR}/natural_cv_errors.txt" || true

# ------------------------------------------------------------------------------
conda deactivate || true

echo "[$(date)] DONE."
echo '**************************************************'
exit 0

```

---

### `6a.tgc.ecs.beagle.imputation.sh`

```bash
#!/bin/bash
#-------------------------------------------------------------------------------
# TreeGeneClimate (TGC) — ECS
# Step 6a: (i) Remove monomorphic sites within each cohort (poly-within-cohort)
#          (ii) BEAGLE imputation per cohort
#               - robust to contigs with only 1 variant (skip those for Beagle,
#                 then append them back unchanged)
#
# Project root:
#   /path/to/your/project
#
# Input (from Step 5a split):
#   RESULTS/ECS/VCF_SPLIT/tgc.ecs.breeding.call.filt.maf05.snvs.vcf.gz
#   RESULTS/ECS/VCF_SPLIT/tgc.ecs.natural.call.filt.maf05.snvs.vcf.gz
#
# Output:
#   RESULTS/ECS/VCF_SPLIT/*.poly.vcf.gz
#   RESULTS/ECS/VCF_SPLIT/*.poly.imputed.vcf.gz         (final, merged)
#   RESULTS/ECS/VCF_SPLIT/BEAGLE_TMP/<COHORT>/...        (intermediates)
#-------------------------------------------------------------------------------

#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --job-name=ECS.POLY.BEAGLE
#SBATCH --output=/path/to/your/project/LOGS/%x_%j.out
#SBATCH --error=/path/to/your/project/LOGS/%x_%j.err
#SBATCH --time=48:00:00
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=YOUR_EMAIL

set -euo pipefail

# ----------------------------- helpers ----------------------------------------
ts() { date +"%a %b %d %H:%M:%S %Z %Y"; }

die() { echo "[$(ts)] ERROR: $*" >&2; exit 1; }

need_file() { [[ -s "$1" ]] || die "Missing/empty file: $1"; }

need_cmd() { command -v "$1" >/dev/null 2>&1 || die "Command not found in PATH: $1"; }

count_sites() {
  # fast if indexed; otherwise falls back to counting lines (slower)
  local vcf="$1"
  if [[ -s "${vcf}.csi" || -s "${vcf}.tbi" ]]; then
    bcftools index -n "$vcf"
  else
    bcftools view -H "$vcf" | wc -l
  fi
}

missing_genotypes_sn() {
  # Extract "number of missing genotypes" from bcftools stats (SN line).
  # Returns integer, or "NA" if not found.
  local vcf="$1"
  local val
  val=$(bcftools stats -s - "$vcf" 2>/dev/null | awk -F'\t' '$1=="SN" && $3=="number of missing genotypes:" {print $4; exit}')
  [[ -n "${val:-}" ]] && echo "$val" || echo "NA"
}

make_contig_lists_by_variant_count() {
  # Writes two files:
  #   contigs_ge2.txt : contigs with >=2 variants
  #   contigs_eq1.txt : contigs with exactly 1 variant
  local vcf="$1"
  local out_ge2="$2"
  local out_eq1="$3"

  bcftools query -f '%CHROM\n' "$vcf" \
    | sort \
    | uniq -c \
    | awk '
        $1==1 {print $2 > eq1}
        $1>=2 {print $2 > ge2}
      ' ge2="$out_ge2" eq1="$out_eq1"
}

# ----------------------------- environment ------------------------------------
echo "================================================================================"
echo "JobID = ${SLURM_JOB_ID:-NA}"
echo "User = ${SLURM_JOB_USER:-$USER}, Account = ${SLURM_JOB_ACCOUNT:-NA}"
echo "Partition = ${SLURM_JOB_PARTITION:-NA}, Nodelist = ${SLURM_JOB_NODELIST:-NA}"
echo "================================================================================"

echo "[$(ts)] SLURM job started: ${SLURM_JOB_NAME:-no_slurm}"
echo "[$(ts)] Host: $(hostname)"
echo "[$(ts)] CPUs: ${SLURM_CPUS_PER_TASK:-1}"

module purge
module load gcc/14.2.0
module load bcftools/1.19

module load miniforge3/24.3.0-0
# robust conda init in batch
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate beagle

need_cmd bcftools
need_cmd beagle

echo "[$(ts)] bcftools: $(bcftools --version | head -n 1)"
echo "[$(ts)] beagle:   $(which beagle)"
echo "[$(ts)] beagle bin: $(file -b "$(which beagle)" || true)"

# ----------------------------- config -----------------------------------------
# === USER CONFIGURATION ===
PROJECT_ROOT="/path/to/your/project"  # <-- set this
# ===========================
LOGS="${PROJECT_ROOT}/LOGS"
SPLIT_DIR="${PROJECT_ROOT}/RESULTS/ECS/VCF_SPLIT"
TMP_BASE="${SPLIT_DIR}/BEAGLE_TMP"

mkdir -p "${LOGS}" "${TMP_BASE}"

BREEDING_IN="${SPLIT_DIR}/tgc.ecs.breeding.call.filt.maf05.snvs.vcf.gz"
NATURAL_IN="${SPLIT_DIR}/tgc.ecs.natural.call.filt.maf05.snvs.vcf.gz"

need_file "${BREEDING_IN}"
need_file "${NATURAL_IN}"

THREADS="${SLURM_CPUS_PER_TASK:-1}"

# Java memory: Beagle wrapper uses Java; _JAVA_OPTIONS is honored by Java.
JAVA_MEM="700g"
export _JAVA_OPTIONS="-Xmx${JAVA_MEM}"

echo "[$(ts)] Java opts: ${_JAVA_OPTIONS}"

# ----------------------------- main -------------------------------------------
run_cohort() {
  local cohort="$1"          # BREEDING / NATURAL
  local vcf_in="$2"          # input VCF (split, cohort)
  local out_prefix="$3"      # output prefix (full path, without extensions)

  local tmp_dir="${TMP_BASE}/${cohort}"
  mkdir -p "${tmp_dir}"

  echo "--------------------------------------------------------------------------------"
  echo "[$(ts)] COHORT: ${cohort}"
  echo "[$(ts)] Input:  ${vcf_in}"
  echo "--------------------------------------------------------------------------------"

  echo "[$(ts)] Indexing input VCF (CSI if missing)"
  if [[ ! -s "${vcf_in}.csi" ]]; then
    bcftools index -c --threads "${THREADS}" "${vcf_in}"
  fi

  echo "[$(ts)] Counting sites BEFORE within-cohort polymorphic filtering"
  local n_before
  n_before=$(count_sites "${vcf_in}")
  echo "  ${cohort} (all sites in split VCF): ${n_before}"

  # 1) poly-within-cohort filter
  local vcf_poly="${out_prefix}.poly.vcf.gz"
  echo "[$(ts)] Filtering to polymorphic-within-cohort sites (AC>0 && AC<AN)"
  bcftools view --threads "${THREADS}" -Oz \
    --include 'AC>0 && AC<AN' \
    -o "${vcf_poly}" \
    "${vcf_in}"

  bcftools index -c --threads "${THREADS}" "${vcf_poly}"

  echo "[$(ts)] Counting sites AFTER within-cohort polymorphic filtering"
  local n_poly
  n_poly=$(count_sites "${vcf_poly}")
  echo "  ${cohort} (poly): ${n_poly}"

  # Missingness before imputation (on poly set)
  echo "[$(ts)] Missing genotypes BEFORE imputation (poly set)"
  local miss_pre
  miss_pre=$(missing_genotypes_sn "${vcf_poly}")
  echo "  ${cohort} missing GT count (poly, pre-impute): ${miss_pre}"

  # 2) prevent Beagle crash on contigs with only 1 variant
  local contigs_ge2="${tmp_dir}/${cohort}.contigs_ge2.txt"
  local contigs_eq1="${tmp_dir}/${cohort}.contigs_eq1.txt"

  echo "[$(ts)] Identifying contigs with >=2 variants vs exactly 1 variant (poly set)"
  : > "${contigs_ge2}"
  : > "${contigs_eq1}"
  make_contig_lists_by_variant_count "${vcf_poly}" "${contigs_ge2}" "${contigs_eq1}"

  local n_contig_ge2 n_contig_eq1
  n_contig_ge2=$(wc -l < "${contigs_ge2}" || echo 0)
  n_contig_eq1=$(wc -l < "${contigs_eq1}" || echo 0)

  echo "  ${cohort} contigs with >=2 variants: ${n_contig_ge2}"
  echo "  ${cohort} contigs with  1 variant : ${n_contig_eq1}"

  # Build subset VCFs
  local vcf_ge2="${tmp_dir}/${cohort}.poly.ge2.vcf.gz"
  local vcf_eq1="${tmp_dir}/${cohort}.poly.eq1.vcf.gz"

  if [[ "${n_contig_ge2}" -gt 0 ]]; then
    echo "[$(ts)] Subsetting to contigs with >=2 variants (for Beagle)"
    # IMPORTANT: contigs_ge2 is a 1-column contig list, so use -r with comma-separated contigs
    local regions_ge2
    regions_ge2=$(paste -sd, "${contigs_ge2}")
    bcftools view --threads "${THREADS}" -Oz \
      -r "${regions_ge2}" \
      -o "${vcf_ge2}" \
      "${vcf_poly}"
    bcftools index -c --threads "${THREADS}" "${vcf_ge2}"
  else
    die "${cohort}: No contigs with >=2 variants found; nothing to impute."
  fi

  if [[ "${n_contig_eq1}" -gt 0 ]]; then
    echo "[$(ts)] Subsetting to contigs with exactly 1 variant (will NOT be imputed; appended back later)"
    # IMPORTANT: contigs_eq1 is a 1-column contig list, so use -r with comma-separated contigs
    local regions_eq1
    regions_eq1=$(paste -sd, "${contigs_eq1}")
    bcftools view --threads "${THREADS}" -Oz \
      -r "${regions_eq1}" \
      -o "${vcf_eq1}" \
      "${vcf_poly}"
    bcftools index -c --threads "${THREADS}" "${vcf_eq1}"
  else
    echo "[$(ts)] No 1-variant contigs for ${cohort} (good)."
  fi

  # missingness specifically in the set that will be imputed
  echo "[$(ts)] Missing genotypes BEFORE imputation (subset sent to Beagle, >=2 variants/contig)"
  local miss_pre_ge2
  miss_pre_ge2=$(missing_genotypes_sn "${vcf_ge2}")
  echo "  ${cohort} missing GT count (ge2 subset, pre-impute): ${miss_pre_ge2}"

  # 3) Beagle on ge2 subset
  local beagle_out_prefix="${tmp_dir}/${cohort}.poly.ge2.imputed"
  local vcf_ge2_imputed="${beagle_out_prefix}.vcf.gz"

  echo "[$(ts)] Running BEAGLE on >=2-variant contigs"
  echo "  Threads: ${THREADS}"
  echo "  Input:   ${vcf_ge2}"
  echo "  Output:  ${vcf_ge2_imputed}"

  beagle \
    gt="${vcf_ge2}" \
    out="${beagle_out_prefix}" \
    nthreads="${THREADS}"

  [[ -s "${vcf_ge2_imputed}" ]] || die "${cohort}: Beagle did not produce output VCF: ${vcf_ge2_imputed}"

  bcftools index -c --threads "${THREADS}" "${vcf_ge2_imputed}"

  echo "[$(ts)] Missing genotypes AFTER imputation (ge2 subset)"
  local miss_post_ge2
  miss_post_ge2=$(missing_genotypes_sn "${vcf_ge2_imputed}")
  echo "  ${cohort} missing GT count (ge2 subset, post-impute): ${miss_post_ge2}"

  # 4) Merge imputed ge2 subset + untouched eq1 subset back into a final poly.imputed VCF
  local vcf_final="${out_prefix}.poly.imputed.vcf.gz"

  if [[ "${n_contig_eq1}" -gt 0 ]]; then
    echo "[$(ts)] Merging imputed (ge2) + untouched (eq1) and sorting"
    bcftools concat -a -Oz \
      "${vcf_ge2_imputed}" \
      "${vcf_eq1}" \
      | bcftools sort -Oz -o "${vcf_final}" -
  else
    echo "[$(ts)] No eq1 subset; final = imputed ge2 (sorted anyway)"
    bcftools sort -Oz -o "${vcf_final}" "${vcf_ge2_imputed}"
  fi

  bcftools index -c --threads "${THREADS}" "${vcf_final}"

  echo "[$(ts)] Final counts and missingness (poly.imputed)"
  local n_final miss_post_all
  n_final=$(count_sites "${vcf_final}")
  miss_post_all=$(missing_genotypes_sn "${vcf_final}")
  echo "  ${cohort} sites (final poly.imputed): ${n_final}"
  echo "  ${cohort} missing GT count (final poly.imputed): ${miss_post_all}"

  if [[ "${n_contig_eq1}" -gt 0 ]]; then
    local n_eq1_sites miss_eq1
    n_eq1_sites=$(count_sites "${vcf_eq1}")
    miss_eq1=$(missing_genotypes_sn "${vcf_eq1}")
    echo "  ${cohort} sites on 1-variant contigs (not imputed): ${n_eq1_sites}"
    echo "  ${cohort} missing GT count on 1-variant contigs:      ${miss_eq1}"
  fi

  echo "[$(ts)] COHORT ${cohort} DONE: ${vcf_final}"
}

echo "[$(ts)] Starting Step 6a: within-cohort poly filtering + Beagle imputation"

BREEDING_PREFIX="${SPLIT_DIR}/tgc.ecs.breeding.call.filt.maf05.snvs"
NATURAL_PREFIX="${SPLIT_DIR}/tgc.ecs.natural.call.filt.maf05.snvs"

run_cohort "BREEDING" "${BREEDING_IN}" "${BREEDING_PREFIX}"
run_cohort "NATURAL"  "${NATURAL_IN}"  "${NATURAL_PREFIX}"

conda deactivate || true

echo "[$(ts)] ALL DONE."
echo "**************************************************"
exit 0

```

---

### `7a.tgc.ecs.plink.ibd.sh`

```bash
#!/bin/bash
#-------------------------------------------------------------------------------
# TreeGeneClimate (TGC) — ECS
# Step 7a: Pairwise relatedness (IBD/PI_HAT) per cohort using PLINK --genome
#          Inputs: LD-pruned BED sets created in Step 5a (admix.pruned)
#          Outputs: .genome files for downstream R (orchestrator becomes Step 8a)
#
# Project root:
#   /path/to/your/project
#
# Inputs (from Step 5a):
#   RESULTS/ECS/POPGEN/STRUCTURE/BREEDING/tgc.ecs.breeding.admix.pruned.{bed,bim,fam}
#   RESULTS/ECS/POPGEN/STRUCTURE/NATURAL/tgc.ecs.natural.admix.pruned.{bed,bim,fam}
#
# Outputs:
#   RESULTS/ECS/POPGEN/RELATEDNESS/IBD/
#     - tgc.ecs.breeding.pruned.ibd.genome
#     - tgc.ecs.natural.pruned.ibd.genome
#     - (plus .log/.nosex)
#-------------------------------------------------------------------------------

#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --job-name=ECS.IBD.PLINK
#SBATCH --output=/path/to/your/project/LOGS/%x_%j.out
#SBATCH --error=/path/to/your/project/LOGS/%x_%j.err
#SBATCH --time=24:00:00
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=YOUR_EMAIL

set -euo pipefail

ts() { date +"%a %b %d %H:%M:%S %Z %Y"; }
die() { echo "[$(ts)] ERROR: $*" >&2; exit 1; }
need_cmd() { command -v "$1" >/dev/null 2>&1 || die "Command not found in PATH: $1"; }
need_file() { [[ -s "$1" ]] || die "Missing/empty file: $1"; }

echo "================================================================================"
echo "JobID = ${SLURM_JOB_ID:-NA}"
echo "User = ${SLURM_JOB_USER:-$USER}, Account = ${SLURM_JOB_ACCOUNT:-NA}"
echo "Partition = ${SLURM_JOB_PARTITION:-NA}, Nodelist = ${SLURM_JOB_NODELIST:-NA}"
echo "================================================================================"
echo "[$(ts)] SLURM job started: ${SLURM_JOB_NAME:-no_slurm}"
echo "[$(ts)] Host: $(hostname)"
echo "[$(ts)] CPUs: ${SLURM_CPUS_PER_TASK:-1}"

# --- modules ---
module purge
module load gcc/14.2.0
module load plink/1.9

need_cmd plink

echo "[$(ts)] plink: $(plink --version 2>&1 | head -n 1)"

# --- paths ---
# === USER CONFIGURATION ===
PROJECT_ROOT="/path/to/your/project"  # <-- set this
# ===========================
LOGS="${PROJECT_ROOT}/LOGS"

STRUCT_DIR="${PROJECT_ROOT}/RESULTS/ECS/POPGEN/STRUCTURE"
BREED_BFILE="${STRUCT_DIR}/BREEDING/tgc.ecs.breeding.admix.pruned"
NATUR_BFILE="${STRUCT_DIR}/NATURAL/tgc.ecs.natural.admix.pruned"

OUT_DIR="${PROJECT_ROOT}/RESULTS/ECS/POPGEN/RELATEDNESS/IBD"
mkdir -p "${LOGS}" "${OUT_DIR}"

# Sanity checks (BED trio)
need_file "${BREED_BFILE}.bed"; need_file "${BREED_BFILE}.bim"; need_file "${BREED_BFILE}.fam"
need_file "${NATUR_BFILE}.bed"; need_file "${NATUR_BFILE}.bim"; need_file "${NATUR_BFILE}.fam"

THREADS="${SLURM_CPUS_PER_TASK:-1}"

run_ibd() {
  local cohort="$1"
  local bfile="$2"
  local outprefix="$3"

  echo "--------------------------------------------------------------------------------"
  echo "[$(ts)] COHORT: ${cohort}"
  echo "[$(ts)] Input bfile: ${bfile}"
  echo "[$(ts)] Output: ${outprefix}.genome"
  echo "--------------------------------------------------------------------------------"

  # --genome full computes PI_HAT + Z0/Z1/Z2 etc.
  # --allow-extra-chr because contig/chrom names are non-standard (PA_cUP..., etc.)
  plink \
    --bfile "${bfile}" \
    --allow-extra-chr \
    --threads "${THREADS}" \
    --genome full \
    --out "${outprefix}"

  need_file "${outprefix}.genome"

  # quick summary
  echo "[$(ts)] ${cohort} .genome rows: $(($(wc -l < "${outprefix}.genome") - 1))"
  echo "[$(ts)] ${cohort} PI_HAT quick peek:"
  awk 'NR==1{print;next} NR<=6{print}' "${outprefix}.genome" | column -t || true
}

run_ibd "BREEDING" "${BREED_BFILE}" "${OUT_DIR}/tgc.ecs.breeding.pruned.ibd"
run_ibd "NATURAL"  "${NATUR_BFILE}" "${OUT_DIR}/tgc.ecs.natural.pruned.ibd"

echo "[$(ts)] ALL DONE."
echo "**************************************************"
exit 0

```

---

### `8a.tgc.ecs.orchestrator.R`

```r
############################################################
# TreeGeneClimate (TGC) — ECS
# Step 8a: R Master Orchestrator (CEPH-HDD layout)
#
# PURPOSE
# - One entry-point to:
#   1) Read cohort genotypes (Step 6a imputed VCFs) and metadata
#   2) Convert VCF -> GDS (once) for fast downstream access
#   3) Build GRM + Kinship (VanRaden; K = GRM/2) from imputed genotypes
#   4) Load PLINK IBD (.genome) produced in Step 7a
#   5) Import PCA + ADMIXTURE results from Step 5a (pruned sets)
#   6) Save reusable R objects + tidy tables under RESULTS/ECS/RANALYSIS
#
# NOTES
# - This orchestrator uses GT from the Beagle-imputed VCFs (hard calls).
#   You *can* use DS (dosage) later for GWAS/GRM, but SNPRelate's VCF->GDS
#   path reads GT. With near-zero missingness after imputation,
#   GT-based GRM is acceptable.
# - No plotting here. Plotting lives in separate scripts.
############################################################

suppressPackageStartupMessages({
  library(SNPRelate)  # GDS I/O and SNP utilities
  library(gdsfmt)     # low-level GDS file handle management
  library(dplyr)      # data wrangling
  library(tibble)     # tidy data frames
  library(readr)      # fast TSV output
})

options(stringsAsFactors = FALSE)
set.seed(1)

############################################################
# 0) TOGGLES
# Each block can be run independently once its dependencies exist.
# Set to FALSE to skip a step when rerunning a partial analysis.
############################################################
RUN_GDS_AND_GRM      <- TRUE   # VCF->GDS + GRM/Kinship
RUN_LOAD_IBD         <- TRUE   # load PLINK IBD from step 7a
RUN_IMPORT_PCA       <- TRUE   # read PLINK PCA outputs from step 5a
RUN_IMPORT_ADMIXTURE <- TRUE   # read ADMIXTURE Q/P + CV summaries from step 5a
RUN_SAVE_DAPC_INPUT  <- TRUE   # save a QC-filtered imputed dosage matrix for DAPC (no DAPC run)

# MAF and missingness thresholds applied when preparing the DAPC dosage matrix.
# These are intentionally lenient: DAPC is a descriptive ordination, not a test.
DAPC_MAF_MIN  <- 0.05
DAPC_MISS_MAX <- 0.10

############################################################
# 1) PATHS (CEPH-HDD)
############################################################
# === USER CONFIGURATION ===
# Set TGC_PROJECT_ROOT as an environment variable, or edit the fallback path below
PROJECT_ROOT <- Sys.getenv("TGC_PROJECT_ROOT",
  unset = "/path/to/your/project")
# ===========================

# Genotypes and population-genetics outputs from the pipeline
VCF_SPLIT_DIR <- file.path(PROJECT_ROOT, "RESULTS/ECS/VCF_SPLIT")
STRUCT_DIR    <- file.path(PROJECT_ROOT, "RESULTS/ECS/POPGEN/STRUCTURE")
IBD_DIR       <- file.path(PROJECT_ROOT, "RESULTS/ECS/POPGEN/RELATEDNESS/IBD")

# R analysis outputs
RANA_DIR   <- file.path(PROJECT_ROOT, "RESULTS/ECS/RANALYSIS")
RDATA_DIR  <- file.path(RANA_DIR, "RDATA")   # RDS objects consumed by downstream scripts
TABLES_DIR <- file.path(RANA_DIR, "TABLES")  # human-readable TSV exports

dir.create(RDATA_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(TABLES_DIR, recursive = TRUE, showWarnings = FALSE)

# Subdirectories under TABLES/ mirror the analysis blocks below
TAB_GRM   <- file.path(TABLES_DIR, "grm_kinship")
TAB_IBD   <- file.path(TABLES_DIR, "ibd")
TAB_PCA   <- file.path(TABLES_DIR, "pca")
TAB_ADMIX <- file.path(TABLES_DIR, "admixture")
TAB_DAPC  <- file.path(TABLES_DIR, "dapc_inputs")

dir.create(TAB_GRM,   recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_IBD,   recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_PCA,   recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_ADMIX, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DAPC,  recursive = TRUE, showWarnings = FALSE)

# Step 6a outputs (imputed)
VCF_BREED_IMP <- file.path(VCF_SPLIT_DIR, "tgc.ecs.breeding.call.filt.maf05.snvs.poly.imputed.vcf.gz")
VCF_NATUR_IMP <- file.path(VCF_SPLIT_DIR, "tgc.ecs.natural.call.filt.maf05.snvs.poly.imputed.vcf.gz")

# Step 7a outputs (IBD on pruned sets)
IBD_BREED <- file.path(IBD_DIR, "tgc.ecs.breeding.pruned.ibd.genome")
IBD_NATUR <- file.path(IBD_DIR, "tgc.ecs.natural.pruned.ibd.genome")

# Step 5a PCA outputs (pruned)
PCA_BREED_EIGENVEC <- file.path(STRUCT_DIR, "BREEDING", "PCA", "tgc.ecs.breeding.pca.pruned.eigenvec")
PCA_BREED_EIGENVAL <- file.path(STRUCT_DIR, "BREEDING", "PCA", "tgc.ecs.breeding.pca.pruned.eigenval")
PCA_NATUR_EIGENVEC <- file.path(STRUCT_DIR, "NATURAL",  "PCA", "tgc.ecs.natural.pca.pruned.eigenvec")
PCA_NATUR_EIGENVAL <- file.path(STRUCT_DIR, "NATURAL",  "PCA", "tgc.ecs.natural.pca.pruned.eigenval")

# Step 5a ADMIXTURE outputs (pruned bed/bim/fam live in STRUCTURE/*)
ADMIX_BREED_DIR <- file.path(STRUCT_DIR, "BREEDING", "ADMIXTURE")
ADMIX_NATUR_DIR <- file.path(STRUCT_DIR, "NATURAL",  "ADMIXTURE")

# Metadata (expected under PROJECT_ROOT/DATA/METADATA)
META_DIR <- file.path(PROJECT_ROOT, "DATA/METADATA")

# Two-column flat files mapping sample IDs to family/population group labels
BREED_MAP_FILE <- file.path(META_DIR, "breeding_sample2family.txt")
NATUR_MAP_FILE <- file.path(META_DIR, "natural_sample2pop.txt")

PHENO_BREED_FILE <- file.path(META_DIR, "tgc.breeding.phenotypes.txt")
PHENO_NATUR_FILE <- file.path(META_DIR, "tgc.natural.phenotypes.txt")

############################################################
# 2) HELPERS
############################################################
# Halt with a clear message if an expected input is absent
ensure_file <- function(x) if (!file.exists(x)) stop("Missing file: ", x)

# SNPRelate keeps an internal registry of open GDS handles; close any
# dangling ones before opening new files to avoid handle conflicts.
close_open_gds <- function() {
  try({
    lst <- gdsfmt::showfile.gds()
    if (!is.null(lst) && length(lst)) {
      for (i in seq_along(lst)) try(gdsfmt::closefn.gds(lst[[i]]), silent = TRUE)
    }
  }, silent = TRUE)
}

# Convert VCF to GDS only when the GDS does not yet exist (idempotent).
# GDS stores genotypes in a compressed columnar format enabling fast random
# access by sample or SNP ID — essential for large spruce datasets.
vcf_to_gds_if_needed <- function(vcf, gds) {
  if (!file.exists(gds)) {
    message("VCF->GDS: ", basename(vcf), " -> ", basename(gds))
    snpgdsVCF2GDS(vcf.fn = vcf, out.fn = gds, method = "biallelic.only", snpfirstdim = TRUE)
  }
  gds
}

# Read a two-column whitespace-delimited map file (no header) into a
# data frame with standardised column names IID and <col2_name>.
read_map2 <- function(file, col2_name) {
  ensure_file(file)
  df <- read.table(file, header = FALSE, stringsAsFactors = FALSE)
  if (ncol(df) < 2) stop("Map file must have at least 2 columns: ", file)
  df <- df[, 1:2]
  colnames(df) <- c("IID", col2_name)
  df$IID <- as.character(df$IID)
  df[[col2_name]] <- as.character(df[[col2_name]])
  df
}

# Compute the VanRaden (2008) genomic relationship matrix.
# X is an n-samples x p-SNPs matrix of allele counts (0/1/2).
# Each SNP is mean-centred by 2*p_j (expected count under HWE);
# the result is scaled by the sum of heterozygosities so that
# diagonal elements approximate 1 for outbred individuals.
grm_vanraden <- function(X) {
  p <- colMeans(X, na.rm = TRUE) / 2          # allele frequency per SNP
  denom <- 2 * sum(p * (1 - p))               # total expected heterozygosity (scaling factor)
  if (!is.finite(denom) || denom <= 0) stop("Non-positive VanRaden denominator.")
  Xc <- sweep(X, 2, 2 * p, "-")              # centre each SNP column
  tcrossprod(Xc) / denom                      # GRM = Z Z' / denom
}

# Parse a PLINK .genome file (all-pairs IBD estimates) and return
# both a long-format tibble (for plotting) and a square matrix (for mixed models).
read_plink_genome <- function(genome_file, ids = NULL) {
  ensure_file(genome_file)
  df <- read.table(genome_file, header = TRUE, stringsAsFactors = FALSE)
  stopifnot(all(c("IID1", "IID2", "PI_HAT") %in% names(df)))
  long <- tibble(IID1 = df$IID1, IID2 = df$IID2, PI_HAT = df$PI_HAT)
  all_ids <- unique(c(long$IID1, long$IID2))
  if (!is.null(ids)) all_ids <- ids
  # Initialise an n x n matrix; PLINK outputs only the upper triangle,
  # so fill both [a,b] and [b,a] to produce a symmetric matrix.
  mat <- matrix(NA_real_, nrow = length(all_ids), ncol = length(all_ids),
                dimnames = list(all_ids, all_ids))
  diag(mat) <- 1
  for (i in seq_len(nrow(long))) {
    a <- long$IID1[i]; b <- long$IID2[i]; v <- long$PI_HAT[i]
    if (a %in% all_ids && b %in% all_ids) { mat[a, b] <- v; mat[b, a] <- v }
  }
  list(long = long, mat = mat)
}

# Read PLINK PCA output, assigning PC1..PCk column names regardless of
# whether the eigenvec file was written with or without a header row.
read_plink_pca <- function(eigenvec_file, eigenval_file) {
  ensure_file(eigenvec_file); ensure_file(eigenval_file)
  ev <- read.table(eigenvec_file, header = FALSE, stringsAsFactors = FALSE)
  colnames(ev)[1:2] <- c("FID", "IID")
  pcs <- paste0("PC", seq_len(ncol(ev) - 2))
  colnames(ev)[3:ncol(ev)] <- pcs
  eval <- scan(eigenval_file, quiet = TRUE)
  list(scores = as_tibble(ev), eigenval = eval)
}

# Parse the ADMIXTURE cross-validation error log produced by step 5a.
# Each line has the form "CV error (K=<k>): <value>"; regex extracts both fields.
read_cv_errors <- function(file) {
  ensure_file(file)
  x <- readLines(file, warn = FALSE)
  tibble(line = x) |>
    mutate(
      K = as.integer(sub(".*\\(K=([0-9]+)\\).*", "\\1", line)),
      CV = as.numeric(sub(".*:\\s*", "", line))
    ) |>
    select(K, CV) |>
    arrange(K)
}

# Collect the paths and K values of all *.Q files produced by ADMIXTURE,
# capping at max_k to exclude exploratory high-K runs if present.
collect_q_files <- function(structure_dir, admix_dir, prefix, max_k = 30) {
  q <- list.files(structure_dir, pattern = paste0("^", prefix, "\\.[0-9]+\\.Q$"), full.names = TRUE)
  if (!length(q)) return(tibble())
  tibble(Q_file = q) |>
    mutate(K = as.integer(sub(".*\\.(\\d+)\\.Q$", "\\1", Q_file))) |>
    filter(K <= max_k) |>
    arrange(K)
}

# Read a single ADMIXTURE Q file (one row per individual, one column per cluster)
# and attach standard column names Q1..QK.
read_q_matrix <- function(q_file) {
  Q <- as.matrix(read.table(q_file, header = FALSE))
  colnames(Q) <- paste0("Q", seq_len(ncol(Q)))
  Q
}

############################################################
# 3) SANITY CHECKS
# Verify critical inputs before any computation starts so failures
# are caught immediately with informative messages.
############################################################
ensure_file(VCF_BREED_IMP)
ensure_file(VCF_NATUR_IMP)
ensure_file(file.path(STRUCT_DIR, "BREEDING", "tgc.ecs.breeding.admix.pruned.fam"))
ensure_file(file.path(STRUCT_DIR, "NATURAL",  "tgc.ecs.natural.admix.pruned.fam"))

# Group-label maps are optional; the GRM is computed regardless.
# Missing labels appear as "Unknown" in downstream figures.
if (!file.exists(BREED_MAP_FILE)) message("NOTE: breeding map not found at ", BREED_MAP_FILE, " (GRM still runs; group labels will be 'Unknown').")
if (!file.exists(NATUR_MAP_FILE)) message("NOTE: natural map not found at ", NATUR_MAP_FILE, " (GRM still runs; group labels will be 'Unknown').")

############################################################
# 4) GDS + GRM/KINSHIP (from imputed VCFs)
############################################################
if (RUN_GDS_AND_GRM) {
  close_open_gds()

  message("=== BREEDING: GDS + GRM/Kinship (from imputed VCF) ===")
  gds_breed <- file.path(RDATA_DIR, "breeding.imputed.snp.gds")
  vcf_to_gds_if_needed(VCF_BREED_IMP, gds_breed)

  breed_map <- if (file.exists(BREED_MAP_FILE)) read_map2(BREED_MAP_FILE, "Family") else tibble(IID = character(), Family = character())

  # local() confines the GDS file handle and large intermediate matrices to a
  # temporary environment, releasing memory automatically when the block exits.
  local({
    gf <- snpgdsOpen(gds_breed, allow.duplicate = FALSE)
    on.exit(try(snpgdsClose(gf), silent = TRUE), add = TRUE)

    sample_ids <- read.gdsn(index.gdsn(gf, "sample.id"))
    snp_ids    <- read.gdsn(index.gdsn(gf, "snp.id"))

    # Extract hard-call genotypes as a SNP x sample integer matrix (0/1/2).
    # snpfirstdim=TRUE keeps SNPs in rows for efficient column-wise allele
    # frequency calculation; mode() coercion converts to numeric for arithmetic.
    geno <- snpgdsGetGeno(gf, sample.id = sample_ids, snp.id = snp_ids, snpfirstdim = TRUE, with.id = FALSE)
    mode(geno) <- "numeric"
    geno[geno > 2 | geno < 0] <- NA_real_   # guard against unexpected encoding values

    # Assign each sample to its family group; samples absent from the map
    # receive "Unknown" so they remain in the GRM without inflating any group.
    grp <- rep("Unknown", length(sample_ids))
    if (nrow(breed_map) > 0) {
      grp2 <- breed_map$Family[match(sample_ids, breed_map$IID)]
      grp2[is.na(grp2)] <- "Unknown"
      grp <- grp2
    }
    grp <- factor(grp)

    # Drop SNPs with undefined allele frequency (e.g. all-missing loci) before
    # imputing residual missing values.
    p_all <- rowMeans(geno, na.rm = TRUE) / 2
    keep_snp <- is.finite(p_all)
    geno2 <- geno[keep_snp, , drop = FALSE]
    p_all <- p_all[keep_snp]

    # Within-group mean imputation: replace residual missing genotypes with 2*p_g
    # where p_g is the group-specific allele frequency.  This is equivalent to
    # setting the missing call to the within-group expectation and avoids
    # artificially inflating between-family relatedness.
    lev <- levels(grp)
    for (g in lev) {
      idx <- which(grp == g)
      if (!length(idx)) next
      pg <- rowMeans(geno2[, idx, drop = FALSE], na.rm = TRUE) / 2
      pg[!is.finite(pg)] <- p_all[!is.finite(pg)]  # fall back to global freq if group is monomorphic
      miss <- is.na(geno2[, idx, drop = FALSE])
      if (any(miss)) {
        fill <- matrix(2 * pg, nrow = nrow(geno2), ncol = length(idx))
        geno2[, idx][miss] <- fill[miss]
      }
    }
    # Global fallback for any remaining NAs (e.g. samples with no group assignment)
    if (anyNA(geno2)) {
      fill2 <- matrix(2 * p_all, nrow = nrow(geno2), ncol = ncol(geno2))
      geno2[is.na(geno2)] <- fill2[is.na(geno2)]
    }
    geno2[geno2 < 0] <- 0; geno2[geno2 > 2] <- 2  # clip numerical noise

    # Transpose to samples x SNPs for VanRaden function, then drop invariant
    # SNPs (sd == 0) which contribute nothing to the GRM and could cause
    # numerical instability in downstream matrix inversions.
    X <- t(geno2)
    rownames(X) <- sample_ids
    sdv <- apply(X, 2, sd)
    keep_var <- is.finite(sdv) & sdv > 0
    if (!all(keep_var)) X <- X[, keep_var, drop = FALSE]

    G <- grm_vanraden(X); rownames(G) <- colnames(G) <- sample_ids
    K <- G / 2  # kinship matrix: GRM/2 converts to the probability-of-IBD scale

    saveRDS(G, file.path(RDATA_DIR, "breeding_grm_vanraden.rds"))
    saveRDS(K, file.path(RDATA_DIR, "breeding_kinship_vanraden_half.rds"))
    write_tsv(as.data.frame(G) |> rownames_to_column("IID"), file.path(TAB_GRM, "breeding_grm_vanraden.tsv"))
    write_tsv(as.data.frame(K) |> rownames_to_column("IID"), file.path(TAB_GRM, "breeding_kinship_vanraden_half.tsv"))

    annot <- tibble(IID = sample_ids, Family = as.character(grp))
    saveRDS(annot, file.path(RDATA_DIR, "breeding_sample_annotation.rds"))
    write_tsv(annot, file.path(TAB_GRM, "breeding_sample_annotation.tsv"))

    message("BREEDING GRM dims: ", nrow(G), " x ", ncol(G))
  })

  message("=== NATURAL: GDS + GRM/Kinship (from imputed VCF) ===")
  gds_natur <- file.path(RDATA_DIR, "natural.imputed.snp.gds")
  vcf_to_gds_if_needed(VCF_NATUR_IMP, gds_natur)

  natur_map <- if (file.exists(NATUR_MAP_FILE)) read_map2(NATUR_MAP_FILE, "Population") else tibble(IID = character(), Population = character())

  # Identical workflow to BREEDING above; local() again scopes memory to this block.
  local({
    gf <- snpgdsOpen(gds_natur, allow.duplicate = FALSE)
    on.exit(try(snpgdsClose(gf), silent = TRUE), add = TRUE)

    sample_ids <- read.gdsn(index.gdsn(gf, "sample.id"))
    snp_ids    <- read.gdsn(index.gdsn(gf, "snp.id"))

    geno <- snpgdsGetGeno(gf, sample.id = sample_ids, snp.id = snp_ids, snpfirstdim = TRUE, with.id = FALSE)
    mode(geno) <- "numeric"
    geno[geno > 2 | geno < 0] <- NA_real_

    # Map sample IDs to Finnish stand (population) labels
    grp <- rep("Unknown", length(sample_ids))
    if (nrow(natur_map) > 0) {
      grp2 <- natur_map$Population[match(sample_ids, natur_map$IID)]
      grp2[is.na(grp2)] <- "Unknown"
      grp <- grp2
    }
    grp <- factor(grp)

    p_all <- rowMeans(geno, na.rm = TRUE) / 2
    keep_snp <- is.finite(p_all)
    geno2 <- geno[keep_snp, , drop = FALSE]
    p_all <- p_all[keep_snp]

    # Within-population mean imputation (same rationale as BREEDING block above)
    lev <- levels(grp)
    for (g in lev) {
      idx <- which(grp == g)
      if (!length(idx)) next
      pg <- rowMeans(geno2[, idx, drop = FALSE], na.rm = TRUE) / 2
      pg[!is.finite(pg)] <- p_all[!is.finite(pg)]
      miss <- is.na(geno2[, idx, drop = FALSE])
      if (any(miss)) {
        fill <- matrix(2 * pg, nrow = nrow(geno2), ncol = length(idx))
        geno2[, idx][miss] <- fill[miss]
      }
    }
    if (anyNA(geno2)) {
      fill2 <- matrix(2 * p_all, nrow = nrow(geno2), ncol = ncol(geno2))
      geno2[is.na(geno2)] <- fill2[is.na(geno2)]
    }
    geno2[geno2 < 0] <- 0; geno2[geno2 > 2] <- 2

    X <- t(geno2)
    rownames(X) <- sample_ids
    sdv <- apply(X, 2, sd)
    keep_var <- is.finite(sdv) & sdv > 0
    if (!all(keep_var)) X <- X[, keep_var, drop = FALSE]

    G <- grm_vanraden(X); rownames(G) <- colnames(G) <- sample_ids
    K <- G / 2

    saveRDS(G, file.path(RDATA_DIR, "natural_grm_vanraden.rds"))
    saveRDS(K, file.path(RDATA_DIR, "natural_kinship_vanraden_half.rds"))
    write_tsv(as.data.frame(G) |> rownames_to_column("IID"), file.path(TAB_GRM, "natural_grm_vanraden.tsv"))
    write_tsv(as.data.frame(K) |> rownames_to_column("IID"), file.path(TAB_GRM, "natural_kinship_vanraden_half.tsv"))

    annot <- tibble(IID = sample_ids, Population = as.character(grp))
    saveRDS(annot, file.path(RDATA_DIR, "natural_sample_annotation.rds"))
    write_tsv(annot, file.path(TAB_GRM, "natural_sample_annotation.tsv"))

    message("NATURAL GRM dims: ", nrow(G), " x ", ncol(G))
  })
}

############################################################
# 5) IBD (PLINK .genome from step 7a)
# The PLINK MoM estimator (PI_HAT) was computed on LD-pruned SNPs
# (step 7a) and is loaded here for relatedness QC and plotting.
# Using the GRM sample order as the canonical ID set ensures the
# square matrix has consistent dimensions with downstream objects.
############################################################
if (RUN_LOAD_IBD) {
  message("=== IBD: load PLINK .genome (step 7a) ===")
  ensure_file(IBD_BREED)
  ensure_file(IBD_NATUR)

  # Derive the canonical sample order from the already-saved GRM
  breed_ids <- readRDS(file.path(RDATA_DIR, "breeding_grm_vanraden.rds")) |> rownames()
  ibd_b <- read_plink_genome(IBD_BREED, ids = breed_ids)
  saveRDS(ibd_b$long, file.path(RDATA_DIR, "breeding_ibd_plink_long.rds"))
  saveRDS(ibd_b$mat,  file.path(RDATA_DIR, "breeding_ibd_plink_mat.rds"))
  write_tsv(ibd_b$long, file.path(TAB_IBD, "breeding_ibd_plink_long.tsv"))

  natur_ids <- readRDS(file.path(RDATA_DIR, "natural_grm_vanraden.rds")) |> rownames()
  ibd_n <- read_plink_genome(IBD_NATUR, ids = natur_ids)
  saveRDS(ibd_n$long, file.path(RDATA_DIR, "natural_ibd_plink_long.rds"))
  saveRDS(ibd_n$mat,  file.path(RDATA_DIR, "natural_ibd_plink_mat.rds"))
  write_tsv(ibd_n$long, file.path(TAB_IBD, "natural_ibd_plink_long.tsv"))

  message("IBD imported: breeding pairs=", nrow(ibd_b$long), " natural pairs=", nrow(ibd_n$long))
}

############################################################
# 6) PCA (PLINK pruned PCA from step 5a)
# Eigenvectors computed by PLINK on LD-pruned SNPs are imported
# here rather than recomputed, to keep PCA consistent with the
# ADMIXTURE analysis which used the same pruned SNP set.
############################################################
if (RUN_IMPORT_PCA) {
  message("=== PCA: import PLINK eigenvec/eigenval (step 5a) ===")

  p_b <- read_plink_pca(PCA_BREED_EIGENVEC, PCA_BREED_EIGENVAL)
  saveRDS(p_b, file.path(RDATA_DIR, "breeding_pca_pruned.rds"))
  write_tsv(p_b$scores, file.path(TAB_PCA, "breeding_pca_pruned_scores.tsv"))
  # Store eigenvalues with their PC index so variance-explained plots can be
  # computed without reloading the full GDS.
  write_tsv(tibble(PC = seq_along(p_b$eigenval), eigenval = p_b$eigenval),
            file.path(TAB_PCA, "breeding_pca_pruned_eigenval.tsv"))

  p_n <- read_plink_pca(PCA_NATUR_EIGENVEC, PCA_NATUR_EIGENVAL)
  saveRDS(p_n, file.path(RDATA_DIR, "natural_pca_pruned.rds"))
  write_tsv(p_n$scores, file.path(TAB_PCA, "natural_pca_pruned_scores.tsv"))
  write_tsv(tibble(PC = seq_along(p_n$eigenval), eigenval = p_n$eigenval),
            file.path(TAB_PCA, "natural_pca_pruned_eigenval.tsv"))

  message("PCA imported.")
}

############################################################
# 7) ADMIXTURE (step 5a outputs)
# Only file paths and CV error summaries are stored here;
# the Q matrices themselves are large and read on demand by
# plotting scripts via the paths in the saved tibbles.
############################################################
if (RUN_IMPORT_ADMIXTURE) {
  message("=== ADMIXTURE: import Q matrices + CV summaries (step 5a) ===")

  cv_b_file <- file.path(ADMIX_BREED_DIR, "breeding_cv_errors.txt")
  cv_n_file <- file.path(ADMIX_NATUR_DIR, "natural_cv_errors.txt")
  # CV errors are generated by ADMIXTURE's --cv flag and are used to
  # select the optimal number of clusters K (minimum CV = best K).
  if (file.exists(cv_b_file)) {
    cv_b <- read_cv_errors(cv_b_file)
    saveRDS(cv_b, file.path(RDATA_DIR, "breeding_admixture_cv.rds"))
    write_tsv(cv_b, file.path(TAB_ADMIX, "breeding_admixture_cv.tsv"))
  } else {
    message("NOTE: missing breeding CV file: ", cv_b_file)
  }
  if (file.exists(cv_n_file)) {
    cv_n <- read_cv_errors(cv_n_file)
    saveRDS(cv_n, file.path(RDATA_DIR, "natural_admixture_cv.rds"))
    write_tsv(cv_n, file.path(TAB_ADMIX, "natural_admixture_cv.tsv"))
  } else {
    message("NOTE: missing natural CV file: ", cv_n_file)
  }

  # Prefix used when naming ADMIXTURE output files (must match step 5a convention)
  breed_prefix <- "tgc.ecs.breeding.admix.pruned"
  natur_prefix <- "tgc.ecs.natural.admix.pruned"

  # Collect paths and K values for all Q files; downstream plotting scripts
  # iterate over this tibble to build structure bar charts for each K.
  q_b <- collect_q_files(file.path(STRUCT_DIR, "BREEDING"), ADMIX_BREED_DIR, breed_prefix, max_k = 30)
  q_n <- collect_q_files(file.path(STRUCT_DIR, "NATURAL"),  ADMIX_NATUR_DIR, natur_prefix, max_k = 30)

  saveRDS(q_b, file.path(RDATA_DIR, "breeding_admixture_q_files.rds"))
  saveRDS(q_n, file.path(RDATA_DIR, "natural_admixture_q_files.rds"))
  write_tsv(q_b, file.path(TAB_ADMIX, "breeding_admixture_q_files.tsv"))
  write_tsv(q_n, file.path(TAB_ADMIX, "natural_admixture_q_files.tsv"))

  message("ADMIXTURE imported (file indices + CV if available).")
}

############################################################
# 8) DAPC INPUT (QC-filtered imputed dosages; no DAPC run)
#    UPDATED: now also saves SNP map (loc/chr/pos) aligned to X
#
# The DAPC dosage matrix (X) uses a relaxed MAF/missingness filter
# relative to the GRM because DAPC is a descriptive clustering method
# rather than a linear mixed model requiring a positive-definite matrix.
# Invariant SNPs (sd == 0) are still removed because they carry no
# discriminant information and can destabilise the DA step.
############################################################
if (RUN_SAVE_DAPC_INPUT) {
  message("=== DAPC INPUT: save QC-filtered imputed dosage matrices ===")

  # -------------------- BREEDING --------------------
  gds_breed <- file.path(RDATA_DIR, "breeding.imputed.snp.gds")
  ensure_file(gds_breed)
  breed_map <- if (file.exists(BREED_MAP_FILE)) read_map2(BREED_MAP_FILE, "Family") else tibble(IID = character(), Family = character())

  local({
    gf <- snpgdsOpen(gds_breed, allow.duplicate = FALSE)
    on.exit(try(snpgdsClose(gf), silent = TRUE), add = TRUE)

    samp <- read.gdsn(index.gdsn(gf, "sample.id"))

    # Read the full SNP annotation vectors once; used later to build the
    # aligned SNP map after filtering, avoiding multiple GDS traversals.
    # cache full SNP vectors once (minimal overhead; avoids ambiguous mapping)
    snp_id_all <- read.gdsn(index.gdsn(gf, "snp.id"))
    chr_all    <- read.gdsn(index.gdsn(gf, "snp.chromosome"))
    pos_all    <- read.gdsn(index.gdsn(gf, "snp.position"))

    # Compute per-SNP MAF and missing-call rate across all samples in one pass
    stat <- snpgdsSNPRateFreq(gf, with.id = TRUE, sample.id = samp)
    maf  <- stat$MinorFreq
    miss <- stat$MissingRate
    snp_id <- stat$snp.id

    keep <- is.finite(maf) & is.finite(miss) & maf >= DAPC_MAF_MIN & miss <= DAPC_MISS_MAX
    if (sum(keep) == 0) stop("BREEDING DAPC: 0 SNPs pass QC (maf/miss thresholds).")

    snp_id_keep <- snp_id[keep]

    geno <- snpgdsGetGeno(gf, sample.id = samp, snp.id = snp_id_keep, snpfirstdim = TRUE, with.id = FALSE)
    mode(geno) <- "numeric"
    geno[geno > 2 | geno < 0] <- NA_real_

    grp <- rep("Unknown", length(samp))
    if (nrow(breed_map) > 0) {
      g2 <- breed_map$Family[match(samp, breed_map$IID)]
      g2[is.na(g2)] <- "Unknown"
      grp <- g2
    }
    grp <- factor(grp)

    # Remove SNPs that remain with undefined allele frequency after MAF filter
    p_all <- rowMeans(geno, na.rm = TRUE) / 2
    keep2 <- is.finite(p_all)
    geno2 <- geno[keep2, , drop = FALSE]
    p_all <- p_all[keep2]

    snp_id_keep2 <- snp_id_keep[keep2]  # track SNP IDs through each filtering step

    # Within-family mean imputation (same rationale as GRM block in Section 4)
    for (g in levels(grp)) {
      idx <- which(grp == g)
      if (!length(idx)) next
      pg <- rowMeans(geno2[, idx, drop = FALSE], na.rm = TRUE) / 2
      pg[!is.finite(pg)] <- p_all[!is.finite(pg)]
      missm <- is.na(geno2[, idx, drop = FALSE])
      if (any(missm)) {
        fill <- matrix(2 * pg, nrow = nrow(geno2), ncol = length(idx))
        geno2[, idx][missm] <- fill[missm]
      }
    }
    if (anyNA(geno2)) {
      fill2 <- matrix(2 * p_all, nrow = nrow(geno2), ncol = ncol(geno2))
      geno2[is.na(geno2)] <- fill2[is.na(geno2)]
    }
    geno2[geno2 < 0] <- 0; geno2[geno2 > 2] <- 2

    X <- t(geno2); rownames(X) <- samp
    sdv <- apply(X, 2, sd)
    keep_var <- is.finite(sdv) & sdv > 0
    if (!all(keep_var)) X <- X[, keep_var, drop = FALSE]

    snp_id_final <- snp_id_keep2[keep_var]  # final post-QC SNP set

    # Build a SNP annotation table (locus ID, chromosome, position) aligned
    # column-for-column with X, so downstream DAPC loadings can be mapped
    # back to genomic coordinates.
    # build snp map aligned to X columns
    idx_map <- match(snp_id_final, snp_id_all)
    snp_df <- tibble(
      loc = as.character(snp_id_final),
      chr = as.character(chr_all[idx_map]),
      pos = as.integer(pos_all[idx_map])
    )

    # ensure loadings rownames == loc
    colnames(X) <- snp_df$loc

    # Bundle all components needed by downstream DAPC scripts into one RDS,
    # including QC parameters for reproducibility documentation.
    out <- list(
      X = X,
      snp = snp_df,
      group = tibble(IID = samp, Group = as.character(grp)),
      qc = list(maf_min = DAPC_MAF_MIN, miss_max = DAPC_MISS_MAX),
      note = "Imputed dosages derived from GT in Beagle-imputed VCF via SNPRelate GDS."
    )
    out_file <- file.path(RDATA_DIR, sprintf("breeding_dapc_input_maf%.2f_miss%.2f.rds", DAPC_MAF_MIN, DAPC_MISS_MAX))
    saveRDS(out, out_file)
    write_tsv(out$group, file.path(TAB_DAPC, "breeding_groups.tsv"))
    write_tsv(out$snp,   file.path(TAB_DAPC, "breeding_snps_loc_chr_pos.tsv"))
    message("Saved: ", out_file, " [", nrow(X), " x ", ncol(X), "]")
  })

  # -------------------- NATURAL --------------------
  gds_natur <- file.path(RDATA_DIR, "natural.imputed.snp.gds")
  ensure_file(gds_natur)
  natur_map <- if (file.exists(NATUR_MAP_FILE)) read_map2(NATUR_MAP_FILE, "Population") else tibble(IID = character(), Population = character())

  # Identical workflow to BREEDING above
  local({
    gf <- snpgdsOpen(gds_natur, allow.duplicate = FALSE)
    on.exit(try(snpgdsClose(gf), silent = TRUE), add = TRUE)

    samp <- read.gdsn(index.gdsn(gf, "sample.id"))

    snp_id_all <- read.gdsn(index.gdsn(gf, "snp.id"))
    chr_all    <- read.gdsn(index.gdsn(gf, "snp.chromosome"))
    pos_all    <- read.gdsn(index.gdsn(gf, "snp.position"))

    stat <- snpgdsSNPRateFreq(gf, with.id = TRUE, sample.id = samp)
    maf  <- stat$MinorFreq
    miss <- stat$MissingRate
    snp_id <- stat$snp.id

    keep <- is.finite(maf) & is.finite(miss) & maf >= DAPC_MAF_MIN & miss <= DAPC_MISS_MAX
    if (sum(keep) == 0) stop("NATURAL DAPC: 0 SNPs pass QC (maf/miss thresholds).")

    snp_id_keep <- snp_id[keep]

    geno <- snpgdsGetGeno(gf, sample.id = samp, snp.id = snp_id_keep, snpfirstdim = TRUE, with.id = FALSE)
    mode(geno) <- "numeric"
    geno[geno > 2 | geno < 0] <- NA_real_

    grp <- rep("Unknown", length(samp))
    if (nrow(natur_map) > 0) {
      g2 <- natur_map$Population[match(samp, natur_map$IID)]
      g2[is.na(g2)] <- "Unknown"
      grp <- g2
    }
    grp <- factor(grp)

    p_all <- rowMeans(geno, na.rm = TRUE) / 2
    keep2 <- is.finite(p_all)
    geno2 <- geno[keep2, , drop = FALSE]
    p_all <- p_all[keep2]

    snp_id_keep2 <- snp_id_keep[keep2]

    for (g in levels(grp)) {
      idx <- which(grp == g)
      if (!length(idx)) next
      pg <- rowMeans(geno2[, idx, drop = FALSE], na.rm = TRUE) / 2
      pg[!is.finite(pg)] <- p_all[!is.finite(pg)]
      missm <- is.na(geno2[, idx, drop = FALSE])
      if (any(missm)) {
        fill <- matrix(2 * pg, nrow = nrow(geno2), ncol = length(idx))
        geno2[, idx][missm] <- fill[missm]
      }
    }
    if (anyNA(geno2)) {
      fill2 <- matrix(2 * p_all, nrow = nrow(geno2), ncol = ncol(geno2))
      geno2[is.na(geno2)] <- fill2[is.na(geno2)]
    }
    geno2[geno2 < 0] <- 0; geno2[geno2 > 2] <- 2

    X <- t(geno2); rownames(X) <- samp
    sdv <- apply(X, 2, sd)
    keep_var <- is.finite(sdv) & sdv > 0
    if (!all(keep_var)) X <- X[, keep_var, drop = FALSE]

    snp_id_final <- snp_id_keep2[keep_var]

    idx_map <- match(snp_id_final, snp_id_all)
    snp_df <- tibble(
      loc = as.character(snp_id_final),
      chr = as.character(chr_all[idx_map]),
      pos = as.integer(pos_all[idx_map])
    )

    colnames(X) <- snp_df$loc

    out <- list(
      X = X,
      snp = snp_df,
      group = tibble(IID = samp, Group = as.character(grp)),
      qc = list(maf_min = DAPC_MAF_MIN, miss_max = DAPC_MISS_MAX),
      note = "Imputed dosages derived from GT in Beagle-imputed VCF via SNPRelate GDS."
    )
    out_file <- file.path(RDATA_DIR, sprintf("natural_dapc_input_maf%.2f_miss%.2f.rds", DAPC_MAF_MIN, DAPC_MISS_MAX))
    saveRDS(out, out_file)
    write_tsv(out$group, file.path(TAB_DAPC, "natural_groups.tsv"))
    write_tsv(out$snp,   file.path(TAB_DAPC, "natural_snps_loc_chr_pos.tsv"))
    message("Saved: ", out_file, " [", nrow(X), " x ", ncol(X), "]")
  })
}

message("=== 8a DONE. Outputs under: ", RANA_DIR)
sessionInfo()

```

---

### `9a.tgc.ecs.pca.R`

```r
suppressPackageStartupMessages({
  library(dplyr)      # data wrangling
  library(tibble)     # tidy data frames
  library(ggplot2)    # plotting
  library(patchwork)  # multi-panel layout
  library(grid)       # low-level graphics (used implicitly by patchwork)
  library(scales)     # axis formatting (percent_format, pretty_breaks)
})

# ------------------------------------------------------------------------------
# Paths (UPDATED to the new project layout)
# ------------------------------------------------------------------------------
# === USER CONFIGURATION ===
rdata_dir <- "/path/to/your/project/RESULTS/ECS/RANALYSIS/RDATA"
fig_dir   <- "/path/to/your/project/RESULTS/ECS/RANALYSIS/PCA"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# PCA outputs from Step 5a (PLINK PCA on pruned SNP sets)
pca_dir <- "/path/to/your/project/RESULTS/ECS/POPGEN/STRUCTURE"

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
# Plot builders
# ------------------------------------------------------------------------------
# Scree plot: bar + line chart showing proportion of variance explained (PVE)
# per PC. Used to judge how many PCs capture meaningful population structure.
make_scree_plot <- function(eigval, title_label = "a)") {
  pve <- eigval / sum(eigval)  # convert raw eigenvalues to proportions
  df <- tibble(PC = seq_along(eigval), PVE = pve)

  ggplot(df, aes(x = PC, y = PVE)) +
    geom_col(fill = "grey40", width = 0.8) +
    geom_line(aes(y = PVE), color = "grey20") +
    geom_point(color = "grey20", size = 1.4) +
    scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.05))) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    labs(title = title_label, x = "Principal component", y = "Variance explained") +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0),
      panel.grid.minor = element_blank()
    )
}

# PCA scatter for a user-specified pair of PCs (x, y).
# Axis labels include the percentage of variance explained by each PC,
# computed from the full eigenvalue vector so the proportion is global.
make_pca_scatter <- function(df_scores, eigval, x = 1, y = 2, group_col = "Group",
                             palette_named, title_label = "a)") {
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
      plot.title = element_text(hjust = 0),
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
# Scree panel (two plots side-by-side), TIFF 24 cm x 12 cm @ 600 dpi
# Both cohorts on the same figure for direct comparison of explained variance.
# ------------------------------------------------------------------------------
scree_b <- make_scree_plot(eigval_b, title_label = "a)")
scree_n <- make_scree_plot(eigval_n, title_label = "b)")

scree_panel <- scree_b + scree_n + plot_layout(ncol = 2, guides = "collect")

# 600 dpi / LZW compression meets most journal submission requirements.
# PDF and EPS are also written for vector-format submission; PNG at 150 dpi
# provides a lightweight preview.
scree_tiff <- file.path(fig_dir, "pca_scree_breeding_natural.tiff")
tiff(scree_tiff, width = 24, height = 12, units = "cm", res = 600, compression = "lzw")
print(scree_panel)
dev.off()
ggsave(sub("\\.tiff$", ".pdf", scree_tiff), plot = scree_panel, width = 24, height = 12, units = "cm")
ggsave(sub("\\.tiff$", ".eps", scree_tiff), plot = scree_panel, device = cairo_ps, width = 24, height = 12, units = "cm")
ggsave(sub("\\.tiff$", ".png", scree_tiff), plot = scree_panel, device = "png",      width = 24, height = 12, units = "cm", dpi = 150)

# ------------------------------------------------------------------------------
# PCA panel — single 6-plot panel (a–f), TIFF 34 cm x 26 cm @ 600 dpi
#   Row 1: breeding  a) PC1v2  b) PC1v3  c) PC2v3
#   Row 2: natural   d) PC1v2  e) PC1v3  f) PC2v3
#
# Showing PC1vs2, PC1vs3, and PC2vs3 captures the three leading axes of
# genetic differentiation without repeating information; this is a standard
# layout for reporting population structure in forest-tree studies.
# ------------------------------------------------------------------------------
p_b_12 <- make_pca_scatter(scores_b, eigval_b, x = 1, y = 2, group_col = "Group",
                           palette_named = colors.17, title_label = "a)")
p_b_13 <- make_pca_scatter(scores_b, eigval_b, x = 1, y = 3, group_col = "Group",
                           palette_named = colors.17, title_label = "b)")
p_b_23 <- make_pca_scatter(scores_b, eigval_b, x = 2, y = 3, group_col = "Group",
                           palette_named = colors.17, title_label = "c)")

p_n_12 <- make_pca_scatter(scores_n, eigval_n, x = 1, y = 2, group_col = "Group",
                           palette_named = colors.25, title_label = "d)")
p_n_13 <- make_pca_scatter(scores_n, eigval_n, x = 1, y = 3, group_col = "Group",
                           palette_named = colors.25, title_label = "e)")
p_n_23 <- make_pca_scatter(scores_n, eigval_n, x = 2, y = 3, group_col = "Group",
                           palette_named = colors.25, title_label = "f)")

# patchwork: '|' composes panels horizontally, '/' stacks rows, '&' applies
# theme modifications to all panels in the assembled layout simultaneously.
# Guides are collected per row so breeding and natural legends remain separate
# and align to the top of their respective row's legend area.
pca_row_b <- (p_b_12 | p_b_13 | p_b_23) + plot_layout(guides = "collect") &
  theme(legend.position = "right", legend.justification = "top")
pca_row_n <- (p_n_12 | p_n_13 | p_n_23) + plot_layout(guides = "collect") &
  theme(legend.position = "right", legend.justification = "top")
pca_panel <- pca_row_b / pca_row_n

pca_tiff <- file.path(fig_dir, "pca_panel_a-f.tiff")
tiff(pca_tiff, width = 34, height = 26, units = "cm", res = 600, compression = "lzw")
print(pca_panel)
dev.off()
ggsave(sub("\\.tiff$", ".pdf", pca_tiff), plot = pca_panel, width = 34, height = 26, units = "cm")
ggsave(sub("\\.tiff$", ".eps", pca_tiff), plot = pca_panel, device = cairo_ps, width = 34, height = 26, units = "cm")
ggsave(sub("\\.tiff$", ".png", pca_tiff), plot = pca_panel, device = "png",      width = 34, height = 26, units = "cm", dpi = 150)

cat("Saved:\n", scree_tiff, "\n", pca_tiff, "\n", sep = "")

sessionInfo()

```

---

## TMS — Targeted Methylation Sequencing

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

#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH -n 48
#SBATCH -N 1
#SBATCH --job-name=TMS.FQC1
#SBATCH --output=/path/to/your/project/LOGS/%x_%j.out
#SBATCH --error=/path/to/your/project/LOGS/%x_%j.err
#SBATCH --ntasks-per-socket 24
#SBATCH --time=48:00:00
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=YOUR_EMAIL

set -euo pipefail

echo "[$(date)] SLURM job started: ${SLURM_JOB_NAME:-no_slurm}"

module load fastqc/0.11.4
module load anaconda3/2020.11

# === USER CONFIGURATION ===
PROJECT_ROOT="/path/to/your/project"  # <-- set this
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

#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH -n 48
#SBATCH -N 1
#SBATCH --job-name=TMS.TRIM
#SBATCH --output=/path/to/your/project/LOGS/%x_%j.out
#SBATCH --error=/path/to/your/project/LOGS/%x_%j.err
#SBATCH --ntasks-per-socket 24
#SBATCH --time=48:00:00
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=YOUR_EMAIL

set -euo pipefail

echo "[$(date)] SLURM job started: ${SLURM_JOB_NAME:-no_slurm}"

module purge
module load fastqc/0.11.4
module load anaconda3/2020.11
module load trimmomatic/0.36

# === USER CONFIGURATION ===
PROJECT_ROOT="/path/to/your/project"  # <-- set this
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

list='
P001_WA01 P001_WA02 P001_WA03 P001_WA04 P001_WA05 P001_WA06 P001_WA07 P001_WA08 P001_WA09 P001_WA10 P001_WA11 P001_WA12
P001_WB01 P001_WB02 P001_WB03 P001_WB04 P001_WB05 P001_WB06 P001_WB07 P001_WB08 P001_WB09 P001_WB10 P001_WB11 P001_WB12
P001_WC01 P001_WC02 P001_WC03 P001_WC04 P001_WC05 P001_WC06 P001_WC07 P001_WC08 P001_WC09 P001_WC10 P001_WC11 P001_WC12
P001_WD01 P001_WD02 P001_WD03 P001_WD04 P001_WD05 P001_WD06 P001_WD07 P001_WD08 P001_WD09 P001_WD10 P001_WD11 P001_WD12
P001_WE01 P001_WE02 P001_WE03 P001_WE04 P001_WE05 P001_WE06 P001_WE07 P001_WE08 P001_WE09 P001_WE10 P001_WE11 P001_WE12
P001_WF01 P001_WF02 P001_WF03 P001_WF04 P001_WF05 P001_WF06 P001_WF07 P001_WF08 P001_WF09 P001_WF10 P001_WF11 P001_WF12
P001_WG01 P001_WG02 P001_WG03 P001_WG04 P001_WG05 P001_WG06 P001_WG07 P001_WG08 P001_WG09 P001_WG10 P001_WG11 P001_WG12
P001_WH01 P001_WH02 P001_WH03 P001_WH04 P001_WH05 P001_WH06 P001_WH07 P001_WH08 P001_WH09 P001_WH10 P001_WH11 P001_WH12

P002_WA01 P002_WA02 P002_WA03 P002_WA04 P002_WA05 P002_WA06 P002_WA07 P002_WA08 P002_WA09 P002_WA10 P002_WA11 P002_WA12
P002_WB01 P002_WB02 P002_WB03 P002_WB04 P002_WB05 P002_WB06 P002_WB07 P002_WB08 P002_WB09 P002_WB10 P002_WB11 P002_WB12
P002_WC01 P002_WC02 P002_WC03 P002_WC04 P002_WC05 P002_WC06 P002_WC07 P002_WC08 P002_WC09 P002_WC10 P002_WC11 P002_WC12
P002_WD01 P002_WD02 P002_WD03 P002_WD04 P002_WD05 P002_WD06 P002_WD07 P002_WD08 P002_WD09 P002_WD10 P002_WD11 P002_WD12
P002_WE01 P002_WE02 P002_WE03 P002_WE04 P002_WE05 P002_WE06 P002_WE07 P002_WE08 P002_WE09 P002_WE10 P002_WE11 P002_WE12
P002_WF01 P002_WF02 P002_WF03 P002_WF04 P002_WF05 P002_WF06 P002_WF07 P002_WF08 P002_WF09 P002_WF10 P002_WF11 P002_WF12
P002_WG01 P002_WG02 P002_WG03 P002_WG04 P002_WG05 P002_WG06 P002_WG07 P002_WG08 P002_WG09 P002_WG10 P002_WG11 P002_WG12
P002_WH01 P002_WH02 P002_WH03 P002_WH04 P002_WH05 P002_WH06 P002_WH07 P002_WH08 P002_WH09 P002_WH10 P002_WH11 P002_WH12

P003_WA01 P003_WA02 P003_WA03 P003_WA04 P003_WA05 P003_WA06 P003_WA07 P003_WA08 P003_WA09 P003_WA10 P003_WA11 P003_WA12
P003_WB01 P003_WB02 P003_WB03 P003_WB04 P003_WB05 P003_WB06 P003_WB07 P003_WB08 P003_WB09 P003_WB10 P003_WB11 P003_WB12
P003_WC01 P003_WC02

P004_WA01 P004_WA02 P004_WA03 P004_WA04 P004_WA05 P004_WA06 P004_WA07 P004_WA08 P004_WA09 P004_WA10 P004_WA11 P004_WA12
P004_WB01 P004_WB02 P004_WB03 P004_WB04 P004_WB05 P004_WB06 P004_WB07 P004_WB08 P004_WB09 P004_WB10 P004_WB11 P004_WB12
P004_WC01 P004_WC02 P004_WC03 P004_WC04 P004_WC05 P004_WC06 P004_WC07 P004_WC08 P004_WC09 P004_WC10 P004_WC11 P004_WC12
P004_WD01 P004_WD02 P004_WD03 P004_WD04 P004_WD05 P004_WD06 P004_WD07 P004_WD08 P004_WD09 P004_WD10 P004_WD11 P004_WD12
P004_WE01 P004_WE02 P004_WE03 P004_WE04 P004_WE05 P004_WE06 P004_WE07 P004_WE08 P004_WE09 P004_WE10 P004_WE11 P004_WE12
P004_WF01 P004_WF02 P004_WF03 P004_WF04 P004_WF05 P004_WF06 P004_WF07 P004_WF08 P004_WF09 P004_WF10 P004_WF11 P004_WF12
P004_WG01 P004_WG02 P004_WG03 P004_WG04 P004_WG05 P004_WG06 P004_WG07 P004_WG08 P004_WG09 P004_WG10 P004_WG11 P004_WG12
P004_WH01 P004_WH02 P004_WH03 P004_WH04 P004_WH05 P004_WH06 P004_WH07 P004_WH08 P004_WH09 P004_WH10 P004_WH11 P004_WH12

P005_WA01 P005_WA02 P005_WA03 P005_WA04 P005_WA05 P005_WA06 P005_WA07 P005_WA08 P005_WA09 P005_WA10 P005_WA11 P005_WA12
P005_WB01 P005_WB02 P005_WB03 P005_WB04 P005_WB05 P005_WB06 P005_WB07 P005_WB08 P005_WB09 P005_WB10 P005_WB11 P005_WB12
P005_WC01 P005_WC02 P005_WC03 P005_WC04 P005_WC05 P005_WC06 P005_WC07 P005_WC08 P005_WC09 P005_WC10 P005_WC11 P005_WC12
P005_WD01 P005_WD02 P005_WD03 P005_WD04 P005_WD05 P005_WD06 P005_WD07 P005_WD08 P005_WD09 P005_WD10 P005_WD11 P005_WD12
P005_WE01 P005_WE02 P005_WE03 P005_WE04 P005_WE05 P005_WE06 P005_WE07 P005_WE08 P005_WE09 P005_WE10 P005_WE11 P005_WE12
P005_WF01 P005_WF02 P005_WF03 P005_WF04 P005_WF05 P005_WF06 P005_WF07 P005_WF08 P005_WF09 P005_WF10 P005_WF11 P005_WF12
P005_WG01 P005_WG02 P005_WG03 P005_WG04 P005_WG05 P005_WG06 P005_WG07 P005_WG08 P005_WG09 P005_WG10 P005_WG11 P005_WG12
P005_WH01 P005_WH02 P005_WH03 P005_WH04 P005_WH05 P005_WH06 P005_WH07 P005_WH08 P005_WH09 P005_WH10 P005_WH11 P005_WH12

P006_WA01 P006_WA02 P006_WA03 P006_WA04 P006_WA05 P006_WA06 P006_WA07 P006_WA08 P006_WA09 P006_WA10 P006_WA11 P006_WA12
P006_WB01 P006_WB02 P006_WB03 P006_WB04 P006_WB05 P006_WB06 P006_WB07 P006_WB08 P006_WB09 P006_WB10 P006_WB11 P006_WB12
P006_WC01 P006_WC02 P006_WC03 P006_WC04 P006_WC05 P006_WC06 P006_WC07 P006_WC08 P006_WC09 P006_WC10 P006_WC11 P006_WC12
P006_WD01 P006_WD02 P006_WD03 P006_WD04 P006_WD05 P006_WD06 P006_WD07 P006_WD08 P006_WD09 P006_WD10 P006_WD11 P006_WD12
P006_WE01 P006_WE02 P006_WE03 P006_WE04 P006_WE05 P006_WE06 P006_WE07 P006_WE08 P006_WE09 P006_WE10 P006_WE11 P006_WE12
P006_WF01 P006_WF02 P006_WF03 P006_WF04 P006_WF05 P006_WF06 P006_WF07 P006_WF08 P006_WF09 P006_WF10 P006_WF11 P006_WF12
P006_WG01 P006_WG02 P006_WG03 P006_WG04 P006_WG05 P006_WG06 P006_WG07 P006_WG08 P006_WG09 P006_WG10 P006_WG11 P006_WG12
P006_WH01 P006_WH02 P006_WH03 P006_WH04 P006_WH05 P006_WH06 P006_WH07 P006_WH08 P006_WH09 P006_WH10 P006_WH11 P006_WH12

P007_WA01 P007_WA02 P007_WA03 P007_WA04 P007_WA05 P007_WA06 P007_WA07 P007_WA08 P007_WA09 P007_WA10 P007_WA11 P007_WA12
P007_WB01 P007_WB02 P007_WB03 P007_WB04 P007_WB05 P007_WB06 P007_WB07 P007_WB08 P007_WB09 P007_WB10 P007_WB11 P007_WB12
P007_WC01 P007_WC02 P007_WC03 P007_WC04 P007_WC05 P007_WC06 P007_WC07 P007_WC08 P007_WC09 P007_WC10 P007_WC11 P007_WC12
P007_WD01 P007_WD02 P007_WD03 P007_WD04 P007_WD05 P007_WD06 P007_WD07 P007_WD08 P007_WD09 P007_WD10 P007_WD11 P007_WD12
P007_WE01 P007_WE02 P007_WE03 P007_WE04 P007_WE05 P007_WE06 P007_WE07 P007_WE08 P007_WE09 P007_WE10 P007_WE11 P007_WE12
P007_WF01 P007_WF02 P007_WF03 P007_WF04 P007_WF05 P007_WF06 P007_WF07 P007_WF08 P007_WF09 P007_WF10 P007_WF11 P007_WF12
P007_WG01 P007_WG02 P007_WG03 P007_WG04 P007_WG05 P007_WG06 P007_WG07 P007_WG08 P007_WG09 P007_WG10 P007_WG11 P007_WG12
P007_WH01 P007_WH02 P007_WH03 P007_WH04 P007_WH05 P007_WH06 P007_WH07 P007_WH08 P007_WH09 P007_WH10 P007_WH11 P007_WH12

P008_WA01 P008_WA02 P008_WA03 P008_WA04 P008_WA05 P008_WA06 P008_WA07 P008_WA08 P008_WA09 P008_WA10 P008_WA11 P008_WA12
P008_WB01 P008_WB02 P008_WB03 P008_WB04 P008_WB05 P008_WB06
'

for sample in ${list}; do
  r1="${INPUT}/HEL_101202_${sample}_R1.fastq.gz"
  r2="${INPUT}/HEL_101202_${sample}_R2.fastq.gz"

  if [[ ! -s "${r1}" || ! -s "${r2}" ]]; then
    echo "WARNING: missing pair for ${sample}; skipping"
    continue
  fi

  java -jar /usr/product/bioinfo/SL_7.0/BIOINFORMATICS/TRIMMOMATIC/0.36/trimmomatic-0.36.jar PE \
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

#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH -n 96
#SBATCH -N 1
#SBATCH --job-name=TMS.BISMARK
#SBATCH --output=/path/to/your/project/LOGS/%x_%j.out
#SBATCH --error=/path/to/your/project/LOGS/%x_%j.err
#SBATCH --exclusive
#SBATCH --time=48:00:00
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=YOUR_EMAIL

set -euo pipefail

echo "[$(date)] SLURM job started: ${SLURM_JOB_NAME:-no_slurm}"

module purge
module load miniforge3/24.3.0-0
module load gcc/14.2.0
module load bowtie2/2.5.4
module load samtools/1.21
source activate bismark

# === USER CONFIGURATION ===
PROJECT_ROOT="/path/to/your/project"  # <-- set this
# ===========================

INPUT="${PROJECT_ROOT}/DATA/TMS/TRIMMED.FASTQ.TMS"
OUTBASE="${PROJECT_ROOT}/DATA/TMS/MAPPED.FILES.TMS"
REF="${PROJECT_ROOT}/REFERENCE/Picab02_genome.assembly/bismark"

LOGS="${PROJECT_ROOT}/LOGS"
mkdir -p "${OUTBASE}" "${LOGS}"

list='
P001_WA01 P001_WA02 P001_WA03 P001_WA04 P001_WA05 P001_WA06 P001_WA07 P001_WA08 P001_WA09 P001_WA10 P001_WA11 P001_WA12
P001_WB01 P001_WB02 P001_WB03 P001_WB04 P001_WB05 P001_WB06 P001_WB07 P001_WB08 P001_WB09 P001_WB10 P001_WB11 P001_WB12
P001_WC01 P001_WC02 P001_WC03 P001_WC04 P001_WC05 P001_WC06 P001_WC07 P001_WC08 P001_WC09 P001_WC10 P001_WC11 P001_WC12
P001_WD01 P001_WD02 P001_WD03 P001_WD04 P001_WD05 P001_WD06 P001_WD07 P001_WD08 P001_WD09 P001_WD10 P001_WD11 P001_WD12
P001_WE01 P001_WE02 P001_WE03 P001_WE04 P001_WE05 P001_WE06 P001_WE07 P001_WE08 P001_WE09 P001_WE10 P001_WE11 P001_WE12
P001_WF01 P001_WF02 P001_WF03 P001_WF04 P001_WF05 P001_WF06 P001_WF07 P001_WF08 P001_WF09 P001_WF10 P001_WF11 P001_WF12
P001_WG01 P001_WG02 P001_WG03 P001_WG04 P001_WG05 P001_WG06 P001_WG07 P001_WG08 P001_WG09 P001_WG10 P001_WG11 P001_WG12
P001_WH01 P001_WH02 P001_WH03 P001_WH04 P001_WH05 P001_WH06 P001_WH07 P001_WH08 P001_WH09 P001_WH10 P001_WH11 P001_WH12

P002_WA01 P002_WA02 P002_WA03 P002_WA04 P002_WA05 P002_WA06 P002_WA07 P002_WA08 P002_WA09 P002_WA10 P002_WA11 P002_WA12
P002_WB01 P002_WB02 P002_WB03 P002_WB04 P002_WB05 P002_WB06 P002_WB07 P002_WB08 P002_WB09 P002_WB10 P002_WB11 P002_WB12
P002_WC01 P002_WC02 P002_WC03 P002_WC04 P002_WC05 P002_WC06 P002_WC07 P002_WC08 P002_WC09 P002_WC10 P002_WC11 P002_WC12
P002_WD01 P002_WD02 P002_WD03 P002_WD04 P002_WD05 P002_WD06 P002_WD07 P002_WD08 P002_WD09 P002_WD10 P002_WD11 P002_WD12
P002_WE01 P002_WE02 P002_WE03 P002_WE04 P002_WE05 P002_WE06 P002_WE07 P002_WE08 P002_WE09 P002_WE10 P002_WE11 P002_WE12
P002_WF01 P002_WF02 P002_WF03 P002_WF04 P002_WF05 P002_WF06 P002_WF07 P002_WF08 P002_WF09 P002_WF10 P002_WF11 P002_WF12
P002_WG01 P002_WG02 P002_WG03 P002_WG04 P002_WG05 P002_WG06 P002_WG07 P002_WG08 P002_WG09 P002_WG10 P002_WG11 P002_WG12
P002_WH01 P002_WH02 P002_WH03 P002_WH04 P002_WH05 P002_WH06 P002_WH07 P002_WH08 P002_WH09 P002_WH10 P002_WH11 P002_WH12

P003_WA01 P003_WA02 P003_WA03 P003_WA04 P003_WA05 P003_WA06 P003_WA07 P003_WA08 P003_WA09 P003_WA10 P003_WA11 P003_WA12
P003_WB01 P003_WB02 P003_WB03 P003_WB04 P003_WB05 P003_WB06 P003_WB07 P003_WB08 P003_WB09 P003_WB10 P003_WB11 P003_WB12
P003_WC01 P003_WC02

P004_WA01 P004_WA02 P004_WA03 P004_WA04 P004_WA05 P004_WA06 P004_WA07 P004_WA08 P004_WA09 P004_WA10 P004_WA11 P004_WA12
P004_WB01 P004_WB02 P004_WB03 P004_WB04 P004_WB05 P004_WB06 P004_WB07 P004_WB08 P004_WB09 P004_WB10 P004_WB11 P004_WB12
P004_WC01 P004_WC02 P004_WC03 P004_WC04 P004_WC05 P004_WC06 P004_WC07 P004_WC08 P004_WC09 P004_WC10 P004_WC11 P004_WC12
P004_WD01 P004_WD02 P004_WD03 P004_WD04 P004_WD05 P004_WD06 P004_WD07 P004_WD08 P004_WD09 P004_WD10 P004_WD11 P004_WD12
P004_WE01 P004_WE02 P004_WE03 P004_WE04 P004_WE05 P004_WE06 P004_WE07 P004_WE08 P004_WE09 P004_WE10 P004_WE11 P004_WE12
P004_WF01 P004_WF02 P004_WF03 P004_WF04 P004_WF05 P004_WF06 P004_WF07 P004_WF08 P004_WF09 P004_WF10 P004_WF11 P004_WF12
P004_WG01 P004_WG02 P004_WG03 P004_WG04 P004_WG05 P004_WG06 P004_WG07 P004_WG08 P004_WG09 P004_WG10 P004_WG11 P004_WG12
P004_WH01 P004_WH02 P004_WH03 P004_WH04 P004_WH05 P004_WH06 P004_WH07 P004_WH08 P004_WH09 P004_WH10 P004_WH11 P004_WH12

P005_WA01 P005_WA02 P005_WA03 P005_WA04 P005_WA05 P005_WA06 P005_WA07 P005_WA08 P005_WA09 P005_WA10 P005_WA11 P005_WA12
P005_WB01 P005_WB02 P005_WB03 P005_WB04 P005_WB05 P005_WB06 P005_WB07 P005_WB08 P005_WB09 P005_WB10 P005_WB11 P005_WB12
P005_WC01 P005_WC02 P005_WC03 P005_WC04 P005_WC05 P005_WC06 P005_WC07 P005_WC08 P005_WC09 P005_WC10 P005_WC11 P005_WC12
P005_WD01 P005_WD02 P005_WD03 P005_WD04 P005_WD05 P005_WD06 P005_WD07 P005_WD08 P005_WD09 P005_WD10 P005_WD11 P005_WD12
P005_WE01 P005_WE02 P005_WE03 P005_WE04 P005_WE05 P005_WE06 P005_WE07 P005_WE08 P005_WE09 P005_WE10 P005_WE11 P005_WE12
P005_WF01 P005_WF02 P005_WF03 P005_WF04 P005_WF05 P005_WF06 P005_WF07 P005_WF08 P005_WF09 P005_WF10 P005_WF11 P005_WF12
P005_WG01 P005_WG02 P005_WG03 P005_WG04 P005_WG05 P005_WG06 P005_WG07 P005_WG08 P005_WG09 P005_WG10 P005_WG11 P005_WG12
P005_WH01 P005_WH02 P005_WH03 P005_WH04 P005_WH05 P005_WH06 P005_WH07 P005_WH08 P005_WH09 P005_WH10 P005_WH11 P005_WH12

P006_WA01 P006_WA02 P006_WA03 P006_WA04 P006_WA05 P006_WA06 P006_WA07 P006_WA08 P006_WA09 P006_WA10 P006_WA11 P006_WA12
P006_WB01 P006_WB02 P006_WB03 P006_WB04 P006_WB05 P006_WB06 P006_WB07 P006_WB08 P006_WB09 P006_WB10 P006_WB11 P006_WB12
P006_WC01 P006_WC02 P006_WC03 P006_WC04 P006_WC05 P006_WC06 P006_WC07 P006_WC08 P006_WC09 P006_WC10 P006_WC11 P006_WC12
P006_WD01 P006_WD02 P006_WD03 P006_WD04 P006_WD05 P006_WD06 P006_WD07 P006_WD08 P006_WD09 P006_WD10 P006_WD11 P006_WD12
P006_WE01 P006_WE02 P006_WE03 P006_WE04 P006_WE05 P006_WE06 P006_WE07 P006_WE08 P006_WE09 P006_WE10 P006_WE11 P006_WE12
P006_WF01 P006_WF02 P006_WF03 P006_WF04 P006_WF05 P006_WF06 P006_WF07 P006_WF08 P006_WF09 P006_WF10 P006_WF11 P006_WF12
P006_WG01 P006_WG02 P006_WG03 P006_WG04 P006_WG05 P006_WG06 P006_WG07 P006_WG08 P006_WG09 P006_WG10 P006_WG11 P006_WG12
P006_WH01 P006_WH02 P006_WH03 P006_WH04 P006_WH05 P006_WH06 P006_WH07 P006_WH08 P006_WH09 P006_WH10 P006_WH11 P006_WH12

P007_WA01 P007_WA02 P007_WA03 P007_WA04 P007_WA05 P007_WA06 P007_WA07 P007_WA08 P007_WA09 P007_WA10 P007_WA11 P007_WA12
P007_WB01 P007_WB02 P007_WB03 P007_WB04 P007_WB05 P007_WB06 P007_WB07 P007_WB08 P007_WB09 P007_WB10 P007_WB11 P007_WB12
P007_WC01 P007_WC02 P007_WC03 P007_WC04 P007_WC05 P007_WC06 P007_WC07 P007_WC08 P007_WC09 P007_WC10 P007_WC11 P007_WC12
P007_WD01 P007_WD02 P007_WD03 P007_WD04 P007_WD05 P007_WD06 P007_WD07 P007_WD08 P007_WD09 P007_WD10 P007_WD11 P007_WD12
P007_WE01 P007_WE02 P007_WE03 P007_WE04 P007_WE05 P007_WE06 P007_WE07 P007_WE08 P007_WE09 P007_WE10 P007_WE11 P007_WE12
P007_WF01 P007_WF02 P007_WF03 P007_WF04 P007_WF05 P007_WF06 P007_WF07 P007_WF08 P007_WF09 P007_WF10 P007_WF11 P007_WF12
P007_WG01 P007_WG02 P007_WG03 P007_WG04 P007_WG05 P007_WG06 P007_WG07 P007_WG08 P007_WG09 P007_WG10 P007_WG11 P007_WG12
P007_WH01 P007_WH02 P007_WH03 P007_WH04 P007_WH05 P007_WH06 P007_WH07 P007_WH08 P007_WH09 P007_WH10 P007_WH11 P007_WH12

P008_WA01 P008_WA02 P008_WA03 P008_WA04 P008_WA05 P008_WA06 P008_WA07 P008_WA08 P008_WA09 P008_WA10 P008_WA11 P008_WA12
P008_WB01 P008_WB02 P008_WB03 P008_WB04 P008_WB05 P008_WB06
'

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

#SBATCH -p YOUR_PARTITION
#SBATCH -t 00-48:00:00
#SBATCH -N 1
#SBATCH -c 24
#SBATCH --mem=100G
#SBATC -C ssd
#SBATCH --job-name=TMS.MEX
#SBATCH --output=/path/to/your/project/LOGS/%x_%A_%a.out
#SBATCH --error=/path/to/your/project/LOGS/%x_%A_%a.err
#SBATCH --array=20-601%4
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=YOUR_EMAIL

set -euo pipefail

module purge
module load gcc/14.2.0
module load bismark/0.25.1
module load bowtie2/2.5.4
module load samtools/1.21

# === USER CONFIGURATION ===
PROJECT_ROOT="/path/to/your/project"  # <-- set this
# ===========================

INPUT="${PROJECT_ROOT}/DATA/TMS/MAPPED.FILES.TMS"
OUTPUT_BASE="${PROJECT_ROOT}/RESULTS/TMS/METHYLATION_CALLS"
REF="${PROJECT_ROOT}/REFERENCE/Picab02_genome.assembly/bismark"

mkdir -p "${OUTPUT_BASE}"

# -----------------------------
# Full 602-sample list
# -----------------------------
mapfile -t samples << 'EOF'
P001_WA02
P001_WA03
P001_WA04
P001_WA05
P001_WA06
P001_WA07
P001_WA08
P001_WA09
P001_WA10
P001_WA11
P001_WB02
P001_WB03
P001_WB04
P001_WB05
P001_WB06
P001_WB07
P001_WB08
P001_WB09
P001_WB10
P001_WB11
P001_WB12
P001_WC01
P001_WC02
P001_WC03
P001_WC05
P001_WC06
P001_WC07
P001_WC09
P001_WC10
P001_WC11
P001_WC12
P001_WD01
P001_WD02
P001_WD03
P001_WD04
P001_WD05
P001_WD06
P001_WD07
P001_WD08
P001_WD09
P001_WD10
P001_WD11
P001_WD12
P001_WE01
P001_WE02
P001_WE03
P001_WE04
P001_WE05
P001_WE06
P001_WE07
P001_WE08
P001_WE09
P001_WE11
P001_WE12
P001_WF01
P001_WF02
P001_WF03
P001_WF04
P001_WF05
P001_WF06
P001_WF07
P001_WF08
P001_WF09
P001_WF10
P001_WF11
P001_WF12
P001_WG01
P001_WG03
P001_WG04
P001_WG05
P001_WG06
P001_WG07
P001_WG08
P001_WG09
P001_WG10
P001_WG11
P001_WG12
P001_WH01
P001_WH02
P001_WH03
P001_WH04
P001_WH05
P001_WH06
P001_WH07
P001_WH08
P001_WH09
P001_WH10
P001_WH11
P001_WH12
P002_WA01
P002_WA02
P002_WA03
P002_WA04
P002_WA05
P002_WA06
P002_WA07
P002_WA08
P002_WA09
P002_WA10
P002_WA11
P002_WA12
P002_WB01
P002_WB02
P002_WB03
P002_WB04
P002_WB05
P002_WB06
P002_WB07
P002_WB08
P002_WB09
P002_WB10
P002_WB11
P002_WB12
P002_WC01
P002_WC02
P002_WC03
P002_WC04
P002_WC05
P002_WC06
P002_WC07
P002_WC08
P002_WC09
P002_WC10
P002_WC11
P002_WC12
P002_WD01
P002_WD02
P002_WD03
P002_WD04
P002_WD05
P002_WD06
P002_WD07
P002_WD08
P002_WD09
P002_WD10
P002_WD11
P002_WD12
P002_WE01
P002_WE02
P002_WE03
P002_WE04
P002_WE05
P002_WE06
P002_WE07
P002_WE08
P002_WE09
P002_WE10
P002_WE11
P002_WE12
P002_WF01
P002_WF02
P002_WF03
P002_WF04
P002_WF05
P002_WF06
P002_WF07
P002_WF08
P002_WF09
P002_WF10
P002_WF11
P002_WF12
P002_WG01
P002_WG02
P002_WG03
P002_WG04
P002_WG05
P002_WG06
P002_WG07
P002_WG08
P002_WG09
P002_WG10
P002_WG11
P002_WG12
P002_WH01
P002_WH02
P002_WH03
P002_WH04
P002_WH05
P002_WH06
P002_WH07
P002_WH08
P002_WH09
P002_WH10
P002_WH11
P002_WH12
P003_WA01
P003_WA02
P003_WA03
P003_WA04
P003_WA05
P003_WA06
P003_WA07
P003_WA08
P003_WA09
P003_WA10
P003_WA11
P003_WA12
P003_WB01
P003_WB03
P003_WB04
P003_WB05
P003_WB06
P003_WB07
P003_WB08
P003_WB09
P003_WB10
P003_WB11
P003_WB12
P003_WC01
P004_WA01
P004_WA02
P004_WA03
P004_WA04
P004_WA05
P004_WA06
P004_WA07
P004_WA08
P004_WA09
P004_WA10
P004_WA11
P004_WA12
P004_WB01
P004_WB02
P004_WB03
P004_WB04
P004_WB05
P004_WB06
P004_WB07
P004_WB08
P004_WB09
P004_WB10
P004_WB11
P004_WB12
P004_WC01
P004_WC02
P004_WC03
P004_WC04
P004_WC05
P004_WC06
P004_WC07
P004_WC08
P004_WC09
P004_WC10
P004_WC11
P004_WC12
P004_WD01
P004_WD02
P004_WD03
P004_WD04
P004_WD05
P004_WD06
P004_WD07
P004_WD08
P004_WD09
P004_WD10
P004_WD11
P004_WD12
P004_WE01
P004_WE02
P004_WE03
P004_WE04
P004_WE05
P004_WE06
P004_WE07
P004_WE08
P004_WE09
P004_WE10
P004_WE11
P004_WE12
P004_WF01
P004_WF02
P004_WF03
P004_WF04
P004_WF05
P004_WF06
P004_WF07
P004_WF08
P004_WF09
P004_WF10
P004_WF11
P004_WF12
P004_WG01
P004_WG02
P004_WG03
P004_WG04
P004_WG05
P004_WG06
P004_WG07
P004_WG08
P004_WG09
P004_WG10
P004_WG11
P004_WG12
P004_WH01
P004_WH02
P004_WH03
P004_WH04
P004_WH05
P004_WH06
P004_WH07
P004_WH08
P004_WH09
P004_WH10
P004_WH11
P004_WH12
P005_WA01
P005_WA02
P005_WA03
P005_WA04
P005_WA05
P005_WA06
P005_WA07
P005_WA08
P005_WA09
P005_WA10
P005_WA11
P005_WA12
P005_WB01
P005_WB02
P005_WB03
P005_WB04
P005_WB05
P005_WB06
P005_WB07
P005_WB08
P005_WB09
P005_WB10
P005_WB11
P005_WB12
P005_WC01
P005_WC02
P005_WC03
P005_WC04
P005_WC05
P005_WC06
P005_WC07
P005_WC08
P005_WC09
P005_WC10
P005_WC11
P005_WC12
P005_WD01
P005_WD02
P005_WD03
P005_WD04
P005_WD05
P005_WD06
P005_WD07
P005_WD08
P005_WD09
P005_WD10
P005_WD11
P005_WD12
P005_WE01
P005_WE02
P005_WE03
P005_WE04
P005_WE05
P005_WE06
P005_WE07
P005_WE08
P005_WE09
P005_WE10
P005_WE11
P005_WE12
P005_WF01
P005_WF02
P005_WF03
P005_WF04
P005_WF05
P005_WF06
P005_WF07
P005_WF08
P005_WF09
P005_WF10
P005_WF11
P005_WF12
P005_WG01
P005_WG02
P005_WG03
P005_WG04
P005_WG05
P005_WG06
P005_WG07
P005_WG08
P005_WG09
P005_WG10
P005_WG11
P005_WG12
P005_WH01
P005_WH02
P005_WH03
P005_WH04
P005_WH05
P005_WH06
P005_WH07
P005_WH08
P005_WH09
P005_WH10
P005_WH11
P005_WH12
P006_WA01
P006_WA02
P006_WA03
P006_WA04
P006_WA05
P006_WA06
P006_WA07
P006_WA08
P006_WA09
P006_WA10
P006_WA11
P006_WB01
P006_WB02
P006_WB03
P006_WB04
P006_WB05
P006_WB06
P006_WB07
P006_WB08
P006_WB09
P006_WB10
P006_WB11
P006_WB12
P006_WC01
P006_WC02
P006_WC03
P006_WC04
P006_WC05
P006_WC06
P006_WC07
P006_WC08
P006_WC09
P006_WC10
P006_WC12
P006_WD01
P006_WD02
P006_WD03
P006_WD04
P006_WD05
P006_WD06
P006_WD07
P006_WD08
P006_WD09
P006_WD10
P006_WD11
P006_WD12
P006_WE01
P006_WE02
P006_WE03
P006_WE04
P006_WE05
P006_WE06
P006_WE07
P006_WE08
P006_WE09
P006_WE10
P006_WE11
P006_WE12
P006_WF01
P006_WF02
P006_WF03
P006_WF04
P006_WF05
P006_WF06
P006_WF07
P006_WF08
P006_WF09
P006_WF10
P006_WF11
P006_WF12
P006_WG01
P006_WG02
P006_WG03
P006_WG04
P006_WG05
P006_WG06
P006_WG07
P006_WG08
P006_WG09
P006_WG10
P006_WG11
P006_WG12
P006_WH01
P006_WH02
P006_WH03
P006_WH04
P006_WH05
P006_WH06
P006_WH07
P006_WH08
P006_WH09
P006_WH10
P006_WH11
P006_WH12
P007_WA01
P007_WA03
P007_WA04
P007_WA05
P007_WA08
P007_WA09
P007_WA10
P007_WA11
P007_WA12
P007_WB01
P007_WB02
P007_WB03
P007_WB04
P007_WB05
P007_WB06
P007_WB07
P007_WB08
P007_WB09
P007_WB10
P007_WB11
P007_WB12
P007_WC01
P007_WC02
P007_WC03
P007_WC04
P007_WC05
P007_WC06
P007_WC07
P007_WC08
P007_WC09
P007_WC10
P007_WC11
P007_WC12
P007_WD01
P007_WD02
P007_WD03
P007_WD04
P007_WD05
P007_WD06
P007_WD09
P007_WD10
P007_WD11
P007_WD12
P007_WE01
P007_WE02
P007_WE03
P007_WE04
P007_WE05
P007_WE06
P007_WE07
P007_WE08
P007_WE09
P007_WE10
P007_WE11
P007_WE12
P007_WF01
P007_WF02
P007_WF03
P007_WF04
P007_WF05
P007_WF06
P007_WF07
P007_WF08
P007_WF09
P007_WF10
P007_WF11
P007_WF12
P007_WG01
P007_WG02
P007_WG03
P007_WG04
P007_WG05
P007_WG06
P007_WG07
P007_WG08
P007_WG09
P007_WG10
P007_WG11
P007_WG12
P007_WH01
P007_WH02
P007_WH04
P007_WH05
P007_WH06
P007_WH07
P007_WH08
P007_WH09
P007_WH10
P007_WH11
P007_WH12
P008_WA01
P008_WA02
P008_WA03
P008_WA04
P008_WA05
P008_WA06
P008_WA07
P008_WA08
P008_WA09
P008_WA10
P008_WA11
P008_WA12
P008_WB02
P008_WB03
P008_WB04
P008_WB05
P008_WB06
EOF

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

#SBATCH -p YOUR_PARTITION
#SBATCH -t 00-12:00:00
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH --job-name=TMS.COVFILTER
#SBATCH --output=/path/to/your/project/LOGS/%x_%j.out
#SBATCH --error=/path/to/your/project/LOGS/%x_%j.err
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=YOUR_EMAIL

set -euo pipefail

# === USER CONFIGURATION ===
PROJECT_ROOT="/path/to/your/project"  # <-- set this
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

### `5b.tgc.methylkit.filtering.sh`

```bash
#!/bin/bash
#-------------------------------------------------------------------------------
# TreeGeneClimate (TGC) — TMS
# Step 5b: methylKit import from sorted Bismark BAMs (by cohort + context)
#   - Ensure BAM indexes exist (CSI/BAI)
#   - Import methylation calls from BAMs (processBismarkAln)
#   - Build methylRawList (proper constructor + treatment)
#   - Filter by coverage (5 <= cov <= 50)
#   - Unite into methylBase (context-aware min.per.group; integer!)
#   - Save TWO RDS per condition:
#       1) methylBase after unite (cov-filtered only)
#       2) methylBase after MEF filter (cov + MEF)
#   - Write summary CSV with site counts at each stage
#
# INPUT:
#   DATA/TMS/MAPPED.FILES.TMS/*_R1_p_bismark_bt2_pe.sorted.bam
#
# OUTPUT:
#   RESULTS/TMS/RANALYSIS/METHYLKIT_OBJECTS/
#     methylBase_<cohort>_<context>_cov5_50_mpg<mpg>.rds
#     methylBase_<cohort>_<context>_cov5_50_mpg<mpg>_mef0.05.rds
#     summary_<cohort>_<context>_cov5_50_mpg<mpg>_mef0.05.csv
#
# JOB ARRAY MAP (6 tasks):
#   0: BREEDING / CpG
#   1: BREEDING / CHG
#   2: BREEDING / CHH
#   3: NATURAL  / CpG
#   4: NATURAL  / CHG
#   5: NATURAL  / CHH
#-------------------------------------------------------------------------------

#SBATCH -p YOUR_PARTITION
#SBATCH -t 02-00:00:00
#SBATCH -N 1
#SBATCH -c 24
#SBATCH --mem=120G
#SBATCH --job-name=TMS.METHYLKIT
#SBATCH --output=/path/to/your/project/LOGS/%x_%A_%a.out
#SBATCH --error=/path/to/your/project/LOGS/%x_%A_%a.err
#SBATCH --array=0-5%2
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=YOUR_EMAIL

set -euo pipefail

# -----------------------------
# Modules / environment
# -----------------------------
module purge
module load gcc/14.2.0
module load samtools/1.21
module load r/4.5.2

export LC_ALL=C.UTF-8
export LANG=C.UTF-8

# IMPORTANT: point to the library you want SLURM jobs to use
export R_LIBS_USER="/path/to/your/Rlibs"

# -----------------------------
# Config
# -----------------------------
# === USER CONFIGURATION ===
PROJECT_ROOT="/path/to/your/project"  # <-- set this
# ===========================
BAM_DIR="${PROJECT_ROOT}/DATA/TMS/MAPPED.FILES.TMS"
OUT_BASE="${PROJECT_ROOT}/RESULTS/TMS/RANALYSIS/METHYLKIT_OBJECTS"
TMPDIR="${OUT_BASE}/_tmp_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"

mkdir -p "${OUT_BASE}" "${TMPDIR}" "${R_LIBS_USER}"

# -----------------------------
# Map array task -> cohort/context
# -----------------------------
contexts=(CpG CHG CHH)
cohorts=(BREEDING NATURAL)

task="${SLURM_ARRAY_TASK_ID}"
cohort="${cohorts[$(( task / 3 ))]}"
ctx="${contexts[$(( task % 3 ))]}"

echo "============================================================"
echo "TGC — TMS — Step 5b"
echo "Task:   ${task}"
echo "Cohort: ${cohort}"
echo "Ctx:    ${ctx}"
echo "Node:   $(hostname)"
echo "Start:  $(date)"
echo "R:      $(which R)"
echo "R_LIBS_USER: ${R_LIBS_USER}"
echo "============================================================"

# -----------------------------
# Select BAMs by cohort
# -----------------------------
shopt -s nullglob
if [[ "${cohort}" == "BREEDING" ]]; then
  bams=( "${BAM_DIR}"/P00[1-3]_*_R1_p_bismark_bt2_pe.sorted.bam )
else
  bams=( "${BAM_DIR}"/P00[4-8]_*_R1_p_bismark_bt2_pe.sorted.bam )
fi

if [[ ${#bams[@]} -eq 0 ]]; then
  echo "ERROR: No BAMs found for ${cohort} under: ${BAM_DIR}"
  exit 1
fi

echo "Found BAMs: ${#bams[@]}"

# -----------------------------
# Ensure BAM indexes exist (CSI or BAI)
# -----------------------------
echo "Checking/creating BAM indexes (.csi OR .bai) ..."
missing_idx=()
for b in "${bams[@]}"; do
  if [[ ! -f "${b}.csi" && ! -f "${b%.bam}.csi" && ! -f "${b}.bai" && ! -f "${b%.bam}.bai" ]]; then
    missing_idx+=( "$b" )
  fi
done

if [[ ${#missing_idx[@]} -gt 0 ]]; then
  echo "Missing indexes: ${#missing_idx[@]} (creating CSI with samtools index -c)"
  printf "%s\n" "${missing_idx[@]}" \
    | xargs -n 1 -P 6 -I {} samtools index -c -@ 2 "{}"
else
  echo "All indexes present."
fi

# -----------------------------
# Install/check methylKit once (avoid array races)
# -----------------------------
PKG_LOCK="${OUT_BASE}/.pkg_install_lock"
(
  flock -n 9 || exit 0
  Rscript --vanilla - <<'RS'
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
if (!requireNamespace("methylKit", quietly = TRUE)) {
  message("methylKit not found -> installing via BiocManager...")
  install.packages("BiocManager", repos="https://cloud.r-project.org")
  BiocManager::install("methylKit", ask=FALSE, update=FALSE)
}
suppressPackageStartupMessages(library(methylKit))
message("methylKit OK: ", as.character(packageVersion("methylKit")))
RS
) 9>"${PKG_LOCK}"

# -----------------------------
# Write and run R script
# -----------------------------
R_SCRIPT="${TMPDIR}/step5b_${cohort}_${ctx}.R"

cat > "${R_SCRIPT}" << 'RSCRIPT'
#!/usr/bin/env Rscript
############################################################
# TGC — TMS — Step 5b
# Save TWO methylBase RDS per cohort/context:
#   1) after cov + unite
#   2) after cov + unite + MEF
############################################################

.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressPackageStartupMessages(library(methylKit))

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT")
BAM_DIR      <- Sys.getenv("BAM_DIR")
OUT_BASE     <- Sys.getenv("OUT_BASE")
COHORT       <- Sys.getenv("COHORT")
CTX          <- Sys.getenv("CTX")

if (any(c(PROJECT_ROOT, BAM_DIR, OUT_BASE, COHORT, CTX) == "")) {
  stop("Missing env vars: PROJECT_ROOT/BAM_DIR/OUT_BASE/COHORT/CTX")
}

assembly_label <- "Picab02" # label only
mincov <- 5L
maxcov <- 50L
mef_threshold <- 0.05

# Context-aware unite threshold (fraction of samples required at a site)
mpg_frac <- switch(
  CTX,
  "CpG" = 0.80,
  "CHG" = 0.70,
  "CHH" = 0.50,
  0.70
)

# --- MEF filter on methylBase
mef_filter_methylBase <- function(mb, thr = 0.05) {
  df <- getData(mb)
  numCs_cols <- grep("^numCs[0-9]+$", colnames(df), value = TRUE)
  numTs_cols <- grep("^numTs[0-9]+$", colnames(df), value = TRUE)
  if (length(numCs_cols) == 0 || length(numTs_cols) == 0 || length(numCs_cols) != length(numTs_cols)) {
    stop("Could not find matching numCs/numTs columns in methylBase.")
  }
  totalCs <- rowSums(df[, numCs_cols, drop = FALSE], na.rm = TRUE)
  totalTs <- rowSums(df[, numTs_cols, drop = FALSE], na.rm = TRUE)
  total   <- totalCs + totalTs
  mef <- ifelse(total > 0, totalCs / total, NA_real_)
  keep <- !is.na(mef) & (mef > thr)
  list(filtered = mb[keep, ], n_before = nrow(mb), n_after = sum(keep))
}

# --- BAM discovery for cohort
pattern <- if (COHORT == "BREEDING") "^P00[1-3]_" else "^P00[4-8]_"
bam_files <- list.files(
  BAM_DIR,
  pattern = paste0(pattern, ".*_R1_p_bismark_bt2_pe\\.sorted\\.bam$"),
  full.names = TRUE
)
if (length(bam_files) == 0) stop("No BAMs found for cohort: ", COHORT)

bam_files <- sort(bam_files)
sample_ids <- sub("_R1_p_bismark_bt2_pe\\.sorted\\.bam$", "", basename(bam_files))

cat("\n============================================================\n")
cat("Step 5b — ", COHORT, " / ", CTX, "\n", sep="")
cat("Samples: ", length(bam_files), "\n", sep="")
cat("============================================================\n")

# --- Import each BAM
mrl_list <- lapply(seq_along(bam_files), function(i) {
  processBismarkAln(
    location     = bam_files[i],
    sample.id    = sample_ids[i],
    assembly     = assembly_label,
    read.context = CTX
  )
})

# methylRawList (dummy treatment; grouping is used later in 6b)
treat <- rep(0L, length(mrl_list))
mrl <- methylRawList(mrl_list, treatment = treat)

# Counts: raw sites per sample
n_sites_raw <- vapply(mrl_list, nrow, integer(1))
raw_min    <- min(n_sites_raw)
raw_median <- as.integer(stats::median(n_sites_raw))
raw_max    <- max(n_sites_raw)
cat("Raw sites/sample (min/median/max): ", raw_min, "/", raw_median, "/", raw_max, "\n", sep="")

# --- Coverage filter per sample
mrl_cov <- filterByCoverage(mrl, lo.count = mincov, hi.count = maxcov)

# Counts: cov sites per sample
mrl_cov_list <- as.list(mrl_cov)
n_sites_cov <- vapply(mrl_cov_list, nrow, integer(1))
cov_min    <- min(n_sites_cov)
cov_median <- as.integer(stats::median(n_sites_cov))
cov_max    <- max(n_sites_cov)
cat("Cov-filter sites/sample (min/median/max): ", cov_min, "/", cov_median, "/", cov_max, "\n", sep="")

# --- Unite (must be integer)
n <- length(mrl_cov)
min_per_group <- as.integer(max(1L, ceiling(mpg_frac * n)))
cat("Using min.per.group = ", min_per_group, " (", mpg_frac, " of n=", n, ")\n", sep="")

mb_unite <- unite(mrl_cov, destrand = FALSE, min.per.group = min_per_group)
n_after_unite <- nrow(mb_unite)
cat("Sites after unite(): ", format(n_after_unite, big.mark=","), "\n", sep="")

if (n_after_unite == 0) {
  stop(
    "unite() returned 0 sites.\n",
    "Lower mpg_frac (especially for CHH) or relax coverage.\n",
    "Current: mincov=", mincov, " maxcov=", maxcov,
    " mpg_frac=", mpg_frac, " min.per.group=", min_per_group, " n=", n
  )
}

# --- Save RDS #1: after unite (cov only)
out_unite <- file.path(
  OUT_BASE,
  sprintf("methylBase_%s_%s_cov%d_%d_mpg%d.rds",
          tolower(COHORT), tolower(CTX), mincov, maxcov, min_per_group)
)
saveRDS(mb_unite, out_unite)
cat("Saved methylBase (cov+unite): ", out_unite, "\n", sep="")

# --- MEF filter
mf <- mef_filter_methylBase(mb_unite, mef_threshold)
mb_mef <- mf$filtered
n_after_mef <- nrow(mb_mef)
cat("Sites after MEF >", mef_threshold, ": ", format(n_after_mef, big.mark=","), "\n", sep="")

# --- Save RDS #2: after MEF
out_mef <- file.path(
  OUT_BASE,
  sprintf("methylBase_%s_%s_cov%d_%d_mpg%d_mef%.2f.rds",
          tolower(COHORT), tolower(CTX), mincov, maxcov, min_per_group, mef_threshold)
)
saveRDS(mb_mef, out_mef)
cat("Saved methylBase (cov+unite+MEF): ", out_mef, "\n", sep="")

# --- Summary CSV
summary_csv <- file.path(
  OUT_BASE,
  sprintf("summary_%s_%s_cov%d_%d_mpg%d_mef%.2f.csv",
          tolower(COHORT), tolower(CTX), mincov, maxcov, min_per_group, mef_threshold)
)

summary_tab <- data.frame(
  cohort = COHORT,
  context = CTX,
  n_samples = length(bam_files),

  raw_sites_min = raw_min,
  raw_sites_median = raw_median,
  raw_sites_max = raw_max,

  cov_sites_min = cov_min,
  cov_sites_median = cov_median,
  cov_sites_max = cov_max,

  mpg_fraction = mpg_frac,
  min_per_group_used = min_per_group,

  n_sites_after_unite = n_after_unite,
  n_sites_after_mef = n_after_mef,

  rds_cov_unite = out_unite,
  rds_cov_unite_mef = out_mef,
  stringsAsFactors = FALSE
)

write.csv(summary_tab, summary_csv, row.names = FALSE)
cat("Saved summary CSV: ", summary_csv, "\n", sep="")
RSCRIPT

chmod +x "${R_SCRIPT}"

export PROJECT_ROOT BAM_DIR OUT_BASE
export COHORT="${cohort}"
export CTX="${ctx}"

Rscript --vanilla "${R_SCRIPT}"

echo "============================================================"
echo "Done:  $(date)"
echo "Cohort/context finished: ${cohort} / ${ctx}"
echo "Outputs in: ${OUT_BASE}"
echo "============================================================"
```

---

### `5b.tgc.tms.import.data.methylkit.R`

```r
#!/usr/bin/env Rscript
############################################################
# TreeGeneClimate (TGC) — TMS
# Step 5b: methylKit import from Bismark BAMs by cohort + context,
#          filter coverage (5–50x) and MEF (>0.05),
#          save methylBase objects + site-count summary table
#
# INPUT:
#   /path/to/your/project/DATA/TMS/MAPPED.FILES.TMS/
#     <SAMPLE>_R1_p_bismark_bt2_pe.bam
#
# OUTPUT:
#   /path/to/your/project/RESULTS/TMS/RANALYSIS/METHYLKIT_OBJECTS/
#     methylBase_breeding_cpg_cov5_50_mef0.05.rds
#     methylBase_breeding_chg_cov5_50_mef0.05.rds
#     methylBase_breeding_chh_cov5_50_mef0.05.rds
#     methylBase_natural_cpg_cov5_50_mef0.05.rds
#     methylBase_natural_chg_cov5_50_mef0.05.rds
#     methylBase_natural_chh_cov5_50_mef0.05.rds
#     methylKit_step5b_sitecounts_cov5_50_mef0.05.csv
#
# NOTES
# - This avoids CX_report files by importing directly from BAMs.
# - Runs sequentially (one cohort × one context at a time) to keep RAM manageable.
############################################################

suppressPackageStartupMessages({
  library(methylKit)
})

# -----------------------------
# Parameters (EDIT ONLY IF NEEDED)
# -----------------------------
# === USER CONFIGURATION ===
# Set TGC_PROJECT_ROOT as an environment variable, or edit the fallback path below
PROJECT_ROOT <- Sys.getenv("TGC_PROJECT_ROOT",
  unset = "/path/to/your/project")
# ===========================
BAM_DIR      <- file.path(PROJECT_ROOT, "DATA/TMS/MAPPED.FILES.TMS")
OUT_DIR      <- file.path(PROJECT_ROOT, "RESULTS/TMS/RANALYSIS/METHYLKIT_OBJECTS")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

assembly_label <- "Picab02"   # label stored inside methylKit objects (not used for mapping)
mincov <- 5
maxcov <- 50    # upper cap removes likely PCR duplicates / pile-ups
mef_threshold <- 0.05   # minimum fraction of samples covered at a site (see MEF filter below)

contexts <- c("CpG", "CHG", "CHH")

# -----------------------------
# Full sample list (copied from Step 4b1)
# -----------------------------
# P001–P003 = BREEDING cohort; P004–P008 = NATURAL cohort
samples <- c(
  "P001_WA02","P001_WA03","P001_WA04","P001_WA05","P001_WA06","P001_WA07","P001_WA08","P001_WA09","P001_WA10","P001_WA11",
  "P001_WB02","P001_WB03","P001_WB04","P001_WB05","P001_WB06","P001_WB07","P001_WB08","P001_WB09","P001_WB10","P001_WB11","P001_WB12",
  "P001_WC01","P001_WC02","P001_WC03","P001_WC05","P001_WC06","P001_WC07","P001_WC09","P001_WC10","P001_WC11","P001_WC12",
  "P001_WD01","P001_WD02","P001_WD03","P001_WD04","P001_WD05","P001_WD06","P001_WD07","P001_WD08","P001_WD09","P001_WD10","P001_WD11","P001_WD12",
  "P001_WE01","P001_WE02","P001_WE03","P001_WE04","P001_WE05","P001_WE06","P001_WE07","P001_WE08","P001_WE09","P001_WE11","P001_WE12",
  "P001_WF01","P001_WF02","P001_WF03","P001_WF04","P001_WF05","P001_WF06","P001_WF07","P001_WF08","P001_WF09","P001_WF10","P001_WF11","P001_WF12",
  "P001_WG01","P001_WG03","P001_WG04","P001_WG05","P001_WG06","P001_WG07","P001_WG08","P001_WG09","P001_WG10","P001_WG11","P001_WG12",
  "P001_WH01","P001_WH02","P001_WH03","P001_WH04","P001_WH05","P001_WH06","P001_WH07","P001_WH08","P001_WH09","P001_WH10","P001_WH11","P001_WH12",
  "P002_WA01","P002_WA02","P002_WA03","P002_WA04","P002_WA05","P002_WA06","P002_WA07","P002_WA08","P002_WA09","P002_WA10","P002_WA11","P002_WA12",
  "P002_WB01","P002_WB02","P002_WB03","P002_WB04","P002_WB05","P002_WB06","P002_WB07","P002_WB08","P002_WB09","P002_WB10","P002_WB11","P002_WB12",
  "P002_WC01","P002_WC02","P002_WC03","P002_WC04","P002_WC05","P002_WC06","P002_WC07","P002_WC08","P002_WC09","P002_WC10","P002_WC11","P002_WC12",
  "P002_WD01","P002_WD02","P002_WD03","P002_WD04","P002_WD05","P002_WD06","P002_WD07","P002_WD08","P002_WD09","P002_WD10","P002_WD11","P002_WD12",
  "P002_WE01","P002_WE02","P002_WE03","P002_WE04","P002_WE05","P002_WE06","P002_WE07","P002_WE08","P002_WE09","P002_WE10","P002_WE11","P002_WE12",
  "P002_WF01","P002_WF02","P002_WF03","P002_WF04","P002_WF05","P002_WF06","P002_WF07","P002_WF08","P002_WF09","P002_WF10","P002_WF11","P002_WF12",
  "P002_WG01","P002_WG02","P002_WG03","P002_WG04","P002_WG05","P002_WG06","P002_WG07","P002_WG08","P002_WG09","P002_WG10","P002_WG11","P002_WG12",
  "P002_WH01","P002_WH02","P002_WH03","P002_WH04","P002_WH05","P002_WH06","P002_WH07","P002_WH08","P002_WH09","P002_WH10","P002_WH11","P002_WH12",
  "P003_WA01","P003_WA02","P003_WA03","P003_WA04","P003_WA05","P003_WA06","P003_WA07","P003_WA08","P003_WA09","P003_WA10","P003_WA11","P003_WA12",
  "P003_WB01","P003_WB03","P003_WB04","P003_WB05","P003_WB06","P003_WB07","P003_WB08","P003_WB09","P003_WB10","P003_WB11","P003_WB12",
  "P003_WC01",
  "P004_WA01","P004_WA02","P004_WA03","P004_WA04","P004_WA05","P004_WA06","P004_WA07","P004_WA08","P004_WA09","P004_WA10","P004_WA11","P004_WA12",
  "P004_WB01","P004_WB02","P004_WB03","P004_WB04","P004_WB05","P004_WB06","P004_WB07","P004_WB08","P004_WB09","P004_WB10","P004_WB11","P004_WB12",
  "P004_WC01","P004_WC02","P004_WC03","P004_WC04","P004_WC05","P004_WC06","P004_WC07","P004_WC08","P004_WC09","P004_WC10","P004_WC11","P004_WC12",
  "P004_WD01","P004_WD02","P004_WD03","P004_WD04","P004_WD05","P004_WD06","P004_WD07","P004_WD08","P004_WD09","P004_WD10","P004_WD11","P004_WD12",
  "P004_WE01","P004_WE02","P004_WE03","P004_WE04","P004_WE05","P004_WE06","P004_WE07","P004_WE08","P004_WE09","P004_WE10","P004_WE11","P004_WE12",
  "P004_WF01","P004_WF02","P004_WF03","P004_WF04","P004_WF05","P004_WF06","P004_WF07","P004_WF08","P004_WF09","P004_WF10","P004_WF11","P004_WF12",
  "P004_WG01","P004_WG02","P004_WG03","P004_WG04","P004_WG05","P004_WG06","P004_WG07","P004_WG08","P004_WG09","P004_WG10","P004_WG11","P004_WG12",
  "P004_WH01","P004_WH02","P004_WH03","P004_WH04","P004_WH05","P004_WH06","P004_WH07","P004_WH08","P004_WH09","P004_WH10","P004_WH11","P004_WH12",
  "P005_WA01","P005_WA02","P005_WA03","P005_WA04","P005_WA05","P005_WA06","P005_WA07","P005_WA08","P005_WA09","P005_WA10","P005_WA11","P005_WA12",
  "P005_WB01","P005_WB02","P005_WB03","P005_WB04","P005_WB05","P005_WB06","P005_WB07","P005_WB08","P005_WB09","P005_WB10","P005_WB11","P005_WB12",
  "P005_WC01","P005_WC02","P005_WC03","P005_WC04","P005_WC05","P005_WC06","P005_WC07","P005_WC08","P005_WC09","P005_WC10","P005_WC11","P005_WC12",
  "P005_WD01","P005_WD02","P005_WD03","P005_WD04","P005_WD05","P005_WD06","P005_WD07","P005_WD08","P005_WD09","P005_WD10","P005_WD11","P005_WD12",
  "P005_WE01","P005_WE02","P005_WE03","P005_WE04","P005_WE05","P005_WE06","P005_WE07","P005_WE08","P005_WE09","P005_WE10","P005_WE11","P005_WE12",
  "P005_WF01","P005_WF02","P005_WF03","P005_WF04","P005_WF05","P005_WF06","P005_WF07","P005_WF08","P005_WF09","P005_WF10","P005_WF11","P005_WF12",
  "P005_WG01","P005_WG02","P005_WG03","P005_WG04","P005_WG05","P005_WG06","P005_WG07","P005_WG08","P005_WG09","P005_WG10","P005_WG11","P005_WG12",
  "P005_WH01","P005_WH02","P005_WH03","P005_WH04","P005_WH05","P005_WH06","P005_WH07","P005_WH08","P005_WH09","P005_WH10","P005_WH11","P005_WH12",
  "P006_WA01","P006_WA02","P006_WA03","P006_WA04","P006_WA05","P006_WA06","P006_WA07","P006_WA08","P006_WA09","P006_WA10","P006_WA11",
  "P006_WB01","P006_WB02","P006_WB03","P006_WB04","P006_WB05","P006_WB06","P006_WB07","P006_WB08","P006_WB09","P006_WB10","P006_WB11","P006_WB12",
  "P006_WC01","P006_WC02","P006_WC03","P006_WC04","P006_WC05","P006_WC06","P006_WC07","P006_WC08","P006_WC09","P006_WC10","P006_WC12",
  "P006_WD01","P006_WD02","P006_WD03","P006_WD04","P006_WD05","P006_WD06","P006_WD07","P006_WD08","P006_WD09","P006_WD10","P006_WD11","P006_WD12",
  "P006_WE01","P006_WE02","P006_WE03","P006_WE04","P006_WE05","P006_WE06","P006_WE07","P006_WE08","P006_WE09","P006_WE10","P006_WE11","P006_WE12",
  "P006_WF01","P006_WF02","P006_WF03","P006_WF04","P006_WF05","P006_WF06","P006_WF07","P006_WF08","P006_WF09","P006_WF10","P006_WF11","P006_WF12",
  "P006_WG01","P006_WG02","P006_WG03","P006_WG04","P006_WG05","P006_WG06","P006_WG07","P006_WG08","P006_WG09","P006_WG10","P006_WG11","P006_WG12",
  "P006_WH01","P006_WH02","P006_WH03","P006_WH04","P006_WH05","P006_WH06","P006_WH07","P006_WH08","P006_WH09","P006_WH10","P006_WH11","P006_WH12",
  "P007_WA01","P007_WA03","P007_WA04","P007_WA05","P007_WA08","P007_WA09","P007_WA10","P007_WA11","P007_WA12",
  "P007_WB01","P007_WB02","P007_WB03","P007_WB04","P007_WB05","P007_WB06","P007_WB07","P007_WB08","P007_WB09","P007_WB10","P007_WB11","P007_WB12",
  "P007_WC01","P007_WC02","P007_WC03","P007_WC04","P007_WC05","P007_WC06","P007_WC07","P007_WC08","P007_WC09","P007_WC10","P007_WC11","P007_WC12",
  "P007_WD01","P007_WD02","P007_WD03","P007_WD04","P007_WD05","P007_WD06","P007_WD09","P007_WD10","P007_WD11","P007_WD12",
  "P007_WE01","P007_WE02","P007_WE03","P007_WE04","P007_WE05","P007_WE06","P007_WE07","P007_WE08","P007_WE09","P007_WE10","P007_WE11","P007_WE12",
  "P007_WF01","P007_WF02","P007_WF03","P007_WF04","P007_WF05","P007_WF06","P007_WF07","P007_WF08","P007_WF09","P007_WF10","P007_WF11","P007_WF12",
  "P007_WG01","P007_WG02","P007_WG03","P007_WG04","P007_WG05","P007_WG06","P007_WG07","P007_WG08","P007_WG09","P007_WG10","P007_WG11","P007_WG12",
  "P007_WH01","P007_WH02","P007_WH04","P007_WH05","P007_WH06","P007_WH07","P007_WH08","P007_WH09","P007_WH10","P007_WH11","P007_WH12",
  "P008_WA01","P008_WA02","P008_WA03","P008_WA04","P008_WA05","P008_WA06","P008_WA07","P008_WA08","P008_WA09","P008_WA10","P008_WA11","P008_WA12",
  "P008_WB02","P008_WB03","P008_WB04","P008_WB05","P008_WB06"
)

# -----------------------------
# Cohort split (by P001–P003 vs P004–P008)
# -----------------------------
# Regex-based split: plate prefix determines cohort membership
is_breeding <- grepl("^P00[1-3]_", samples)
is_natural  <- grepl("^P00[4-8]_", samples)

breeding_ids <- samples[is_breeding]
natural_ids  <- samples[is_natural]

# -----------------------------
# MEF filter (on a methylBase object)
# MEF = fraction of samples with coverage >0 at the site
# -----------------------------
# Sites covered in very few samples inflate missing-data rates downstream;
# the MEF threshold (>5% of samples) removes such low-prevalence sites.
mef_filter_methylBase <- function(mb, thr = 0.05) {
  cov_mat <- getCoverage(mb)           # coverage matrix: sites x samples
  mef <- rowMeans(cov_mat > 0)         # proportion of samples with any coverage
  keep <- mef > thr
  list(filtered = mb[keep, ], mef = mef, keep = keep)
}

# -----------------------------
# Run sequentially (min RAM)
# -----------------------------
# Accumulate site counts per cohort/context for a summary CSV
summary_tab <- data.frame(
  cohort = character(),
  context = character(),
  n_samples = integer(),
  n_sites_after_cov = numeric(),
  n_sites_after_mef = numeric(),
  stringsAsFactors = FALSE
)

process_one <- function(sample_ids, cohort_label, ctx) {
  message("\n============================================================")
  message("Step 5b — Processing: ", cohort_label, " / ", ctx)
  message("============================================================")

  # coordinate-sorted BAMs for methylKit
  bam_files <- file.path(
    BAM_DIR,
    paste0(sample_ids, "_R1_p_bismark_bt2_pe.sorted.bam")
  )

  exists <- file.exists(bam_files)
  if (!all(exists)) {
    missing <- bam_files[!exists]
    stop("Missing sorted BAMs (showing up to 10):\n", paste(head(missing, 10), collapse = "\n"))
  }

  # Optional: BAM index (.bai). Do NOT stop if missing.
  bai1 <- paste0(bam_files, ".bai")
  bai2 <- sub("\\.bam$", ".bai", bam_files)
  has_bai <- file.exists(bai1) | file.exists(bai2)
  if (!all(has_bai)) {
    missing_bai_n <- sum(!has_bai)
    warning(
      "BAM index (.bai) not found for ", missing_bai_n, " file(s). Continuing without index.\n",
      "If you later need random access or speedups, create with: samtools index <file>.bam"
    )
  }

  message("Samples: ", length(bam_files))

  # processBismarkAln is NOT vectorized -> one BAM at a time
  # mincov=1 here to retain all covered sites; the 5x floor is applied after unite()
  mlist <- vector("list", length(bam_files))
  for (i in seq_along(bam_files)) {
    mlist[[i]] <- processBismarkAln(
      location     = bam_files[i],
      sample.id    = sample_ids[i],
      assembly     = assembly_label,
      read.context = ctx,   # specifies which cytosine context to extract: CpG, CHG, or CHH
      mincov       = 1      # minimal pre-filter; main coverage filter applied post-unite()
    )
    if (i %% 10 == 0) message("  imported ", i, "/", length(bam_files))
  }

  # Combine individual methylRaw objects into a list object required by unite()
  mrl <- new("methylRawList", mlist)

  # unite() intersects sites across all samples (only sites present in every sample are kept)
  mb <- unite(mrl, destrand = FALSE)

  # Coverage filter 5–50x applied after uniting to work on the intersected site set
  mb_cov <- filterByCoverage(mb, lo.count = mincov, hi.count = maxcov)
  n_after_cov <- nrow(mb_cov)

  # MEF filter: remove sites covered in fewer than 5% of samples
  mf <- mef_filter_methylBase(mb_cov, thr = mef_threshold)
  mb_mef <- mf$filtered
  n_after_mef <- nrow(mb_mef)

  message("Sites after coverage (", mincov, "–", maxcov, "x): ", format(n_after_cov, big.mark = ","))
  message("Sites after MEF > ", mef_threshold, ": ", format(n_after_mef, big.mark = ","))

  # Save the filtered methylBase object for use in downstream steps (6b–8b)
  out_rds <- file.path(
    OUT_DIR,
    paste0("methylBase_", tolower(cohort_label), "_", tolower(ctx),
           "_cov", mincov, "_", maxcov, "_mef", mef_threshold, ".rds")
  )
  saveRDS(mb_mef, out_rds)

  # Explicitly free large objects before moving to the next context
  rm(mlist, mrl, mb, mb_cov, mf, mb_mef)
  gc()

  c(n_after_cov = n_after_cov, n_after_mef = n_after_mef)
}

# Process BREEDING cohort across all three contexts
for (ctx in contexts) {
  res <- process_one(breeding_ids, "BREEDING", ctx)
  summary_tab <- rbind(summary_tab, data.frame(
    cohort = "BREEDING",
    context = ctx,
    n_samples = length(breeding_ids),
    n_sites_after_cov = as.numeric(res["n_after_cov"]),
    n_sites_after_mef = as.numeric(res["n_after_mef"]),
    stringsAsFactors = FALSE
  ))
}

# Process NATURAL cohort across all three contexts
for (ctx in contexts) {
  res <- process_one(natural_ids, "NATURAL", ctx)
  summary_tab <- rbind(summary_tab, data.frame(
    cohort = "NATURAL",
    context = ctx,
    n_samples = length(natural_ids),
    n_sites_after_cov = as.numeric(res["n_after_cov"]),
    n_sites_after_mef = as.numeric(res["n_after_mef"]),
    stringsAsFactors = FALSE
  ))
}

# Write site-count summary table for QC reporting
out_csv <- file.path(OUT_DIR, paste0("methylKit_step5b_sitecounts_cov",
                                     mincov, "_", maxcov, "_mef", mef_threshold, ".csv"))
write.csv(summary_tab, out_csv, row.names = FALSE)

message("\n============================================================")
message("Step 5b completed.")
message("Objects saved to: ", OUT_DIR)
message("Summary table:   ", out_csv)
message("============================================================\n")

print(summary_tab)
sessionInfo()

```

---

### `6b.tgc.tms.methylation.levels.R`

```r
#!/usr/bin/env Rscript
############################################################
# TreeGeneClimate (TGC) — TMS
# Step 6b: Plot per-sample mean methylation (%) by group
#   - Breeding cohort: per family (17 families) -> panels A, C, E
#   - Natural cohort:  per stand  (25 stands)   -> panels B, D, F
#   - One plot per context (CpG/CHG/CHH)
#   - Statistical tests + posthoc + logs
#
# STYLE:
#   - Panel titles: single bold capital letter (A–F), no extra text
#   - No legend
#   - Annotation: only p-value (no method label)
#   - Posthoc letters angled; spaced above boxes to reduce overlap
#
# PANEL LAYOUT (FIG34_PANEL_COMBINED):
#   Top row:    A (breeding CpG)  |  B (natural CpG)
#   Middle row: C (breeding CHG)  |  D (natural CHG)
#   Bottom row: E (breeding CHH)  |  F (natural CHH)
#
# INPUT:
#   RDS from Step 5b (after unite) in:
#     .../RESULTS/TBS/RANALYSIS/METHYLKIT_OBJECTS/
#   Expected naming (examples):
#     methylBase_breeding_cpg_cov5_50_mpg168.rds
#     methylBase_breeding_cpg_cov5_50_mpg168_mef0.05.rds
#
# OUTPUT:
#   Figures:
#     .../RESULTS/TBS/RANALYSIS/FIGURES/FIG3_FIG4/
#   Logs + posthoc tables:
#     .../RESULTS/TBS/RANALYSIS/ANOVA.METHYL.LEVEL/
############################################################

## ---------------------------
## Packages (auto-install CRAN if missing)
## ---------------------------
cran_if_missing <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss) > 0) install.packages(miss, repos = "https://cloud.r-project.org")
}

cran_if_missing(c("ggplot2", "multcompView", "patchwork", "car"))
suppressPackageStartupMessages({
  library(methylKit)
  library(ggplot2)
  library(multcompView)
  library(patchwork)
  library(car)          # leveneTest
})

## ---------------------------
## Paths
## ---------------------------
# === USER CONFIGURATION ===
# Set TGC_PROJECT_ROOT as an environment variable, or edit the fallback path below
PROJECT_ROOT <- Sys.getenv("TGC_PROJECT_ROOT",
  unset = "/path/to/your/project")
# ===========================

rds_dir  <- file.path(PROJECT_ROOT, "RESULTS/TBS/RANALYSIS/METHYLKIT_OBJECTS")
fig_dir  <- file.path(PROJECT_ROOT, "RESULTS/TBS/RANALYSIS/FIGURES/FIG3_FIG4")
log_dir  <- file.path(PROJECT_ROOT, "RESULTS/TBS/RANALYSIS/ANOVA.METHYL.LEVEL")

map_file_breeding <- file.path(PROJECT_ROOT, "DATA/METADATA/breeding_sample2family.txt")
map_file_natural  <- file.path(PROJECT_ROOT, "DATA/METADATA/natural_sample2pop.txt")

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)

log_file <- file.path(log_dir, "Step6b_methylation_level_stats_posthoc.log")

## ---------------------------
## Color codes (yours)
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
  
  # with MEF
  pat_mef <- sprintf("^methylBase_%s_%s_cov5_50_mpg[0-9]+_mef0\\.05\\.rds$", cohort, ctx)
  # cov+unite only
  pat_cov <- sprintf("^methylBase_%s_%s_cov5_50_mpg[0-9]+\\.rds$", cohort, ctx)
  
  if (prefer_mef) {
    hits <- list.files(rds_dir, pattern = pat_mef, full.names = TRUE)
    if (length(hits) == 0) {
      # fall back to cov-only
      hits <- list.files(rds_dir, pattern = pat_cov, full.names = TRUE)
    }
  } else {
    hits <- list.files(rds_dir, pattern = pat_cov, full.names = TRUE)
    if (length(hits) == 0) {
      # fall back to MEF
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

# Global test: always ANOVA + TukeyHSD (per-sample means are CLT-justified).
# Shapiro-Wilk and Levene tests are run and reported for transparency but do
# not affect the choice of test.
run_group_tests <- function(df) {
  df <- df[!is.na(df$Group) & !is.na(df$Methylation), ]
  df$Group <- droplevels(df$Group)

  aov_fit <- aov(Methylation ~ Group, data = df)

  sh_p  <- tryCatch(shapiro.test(residuals(aov_fit))$p.value, error = function(e) NA_real_)
  lev_p <- tryCatch(car::leveneTest(Methylation ~ Group, data = df)[["Pr(>F)"]][1], error = function(e) NA_real_)

  an_p    <- summary(aov_fit)[[1]][["Pr(>F)"]][1]
  tuk     <- TukeyHSD(aov_fit)
  letters <- multcompView::multcompLetters4(aov_fit, tuk)$Group$Letters

  list(
    method    = "ANOVA + TukeyHSD",
    shapiro_p = sh_p,
    levene_p  = lev_p,
    global_p  = an_p,
    model_obj = aov_fit,
    posthoc   = tuk,
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

# Plot function (style you requested)
make_boxplot_clean <- function(df, color_vec, title_letter,
                               ylab = "Percent methylation",
                               letter_angle = 45,
                               letter_yfactor = 0.12) {
  
  present <- levels(df$Group)
  color_vec <- color_vec[names(color_vec) %in% present]
  
  ggplot(df, aes(x = Group, y = Methylation, fill = Group)) +
    geom_boxplot(outlier.shape = NA, width = 0.7) +
    geom_jitter(width = 0.2, size = 1, color = "black", alpha = 0.7) +
    scale_fill_manual(values = color_vec, drop = FALSE) +
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
                              letter_yfactor = 0.12) {
  
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
  
  p <- p + annotate(
    "text",
    x = 1,
    y = max(df$Methylation, na.rm = TRUE) + 0.22 * y_span,
    label = fmt_p(tests$global_p),
    hjust = 0, vjust = 1.2, size = 4, fontface = "italic"
  )
  
  # Posthoc letters angled + spaced
  letters_vec <- tests$letters
  letters_vec <- letters_vec[levels(df$Group)]
  
  box_max <- aggregate(Methylation ~ Group, data = df, max, na.rm = TRUE)
  box_max$Letter <- unname(letters_vec[as.character(box_max$Group)])
  box_max$y <- box_max$Methylation + letter_yfactor * y_span
  
  p <- p + geom_text(
    data = box_max,
    aes(x = Group, y = y, label = Letter),
    size = 5, fontface = 1, vjust = 0,
    angle = letter_angle
  )
  
  # Save figure (all four formats independently)
  out_fig <- file.path(fig_dir, sprintf("%s_%s_boxplot.tiff", fig_prefix, ctx))
  save_all_formats(p, out_fig, w_cm = 26, h_cm = 14)
  
  # Save posthoc table (always TukeyHSD)
  posthoc_tsv <- file.path(log_dir, sprintf("%s_%s_posthoc.tsv", fig_prefix, ctx))
  tuk_df <- as.data.frame(tests$posthoc$Group)
  tuk_df$comparison <- rownames(tuk_df)
  rownames(tuk_df) <- NULL
  write.table(tuk_df, posthoc_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
  
  list(df = df, plot = p, tests = tests, fig = out_fig, posthoc = posthoc_tsv, rds = rds_path)
}

## ---------------------------
## RUN ALL 6
## ---------------------------
sink(log_file)
cat("TGC — TBS — Step 6b (updated style)\n")
cat("Timestamp: ", format(Sys.time()), "\n\n", sep = "")

cat("RDS dir:  ", rds_dir, "\n", sep = "")
cat("Fig dir:  ", fig_dir, "\n", sep = "")
cat("Log dir:  ", log_dir, "\n\n", sep = "")

# Global methylation levels must be estimated from all covered loci (cov+unite),
# NOT from the MEF-filtered subset, which selects variable sites and inflates means.
# MEF-filtered data is used downstream for SVMP analysis (Step 8b) only.
prefer_mef <- FALSE

# BREEDING -> panels a/c/e (left column, per context)
res_b_CpG <- analyze_one_clean("BREEDING", "CpG", map_file_breeding, colors.17, "FIG3a_BREEDING", "A",
                               prefer_mef = prefer_mef, letter_angle = 45, letter_yfactor = 0.12)
res_b_CHG <- analyze_one_clean("BREEDING", "CHG", map_file_breeding, colors.17, "FIG3c_BREEDING", "C",
                               prefer_mef = prefer_mef, letter_angle = 45, letter_yfactor = 0.12)
res_b_CHH <- analyze_one_clean("BREEDING", "CHH", map_file_breeding, colors.17, "FIG3e_BREEDING", "E",
                               prefer_mef = prefer_mef, letter_angle = 45, letter_yfactor = 0.12)

# NATURAL -> panels b/d/f (right column, per context)
res_n_CpG <- analyze_one_clean("NATURAL", "CpG", map_file_natural, colors.25, "FIG4b_NATURAL", "B",
                               prefer_mef = prefer_mef, letter_angle = 45, letter_yfactor = 0.12)
res_n_CHG <- analyze_one_clean("NATURAL", "CHG", map_file_natural, colors.25, "FIG4d_NATURAL", "D",
                               prefer_mef = prefer_mef, letter_angle = 45, letter_yfactor = 0.12)
res_n_CHH <- analyze_one_clean("NATURAL", "CHH", map_file_natural, colors.25, "FIG4f_NATURAL", "F",
                               prefer_mef = prefer_mef, letter_angle = 45, letter_yfactor = 0.12)

cat("\n\n====================\nSUMMARY\n====================\n")
summ_line <- function(res, cohort, ctx) {
  cat(cohort, " / ", ctx, "\n", sep = "")
  cat("  Method: ", res$tests$method, "\n", sep = "")
  cat("  Shapiro-Wilk p (residuals): ", format(res$tests$shapiro_p, digits = 4), "\n", sep = "")
  cat("  Levene p:                   ", format(res$tests$levene_p,  digits = 4), "\n", sep = "")
  cat("  ANOVA global p: ", format(res$tests$global_p, digits = 8), "\n", sep = "")
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
  file.path(fig_dir, "FIG34_PANEL_COMBINED_CpG_CHG_CHH.tiff"),
  w_cm = 48, h_cm = 42
)

corrected_dir <- file.path(PROJECT_ROOT, "RESULTS/CORRECTED/FIGURES/NEW")
dir.create(corrected_dir, showWarnings = FALSE, recursive = TRUE)
for (ext in c("tiff", "pdf", "eps", "png")) {
  src <- file.path(fig_dir, paste0("FIG34_PANEL_COMBINED_CpG_CHG_CHH.", ext))
  if (file.exists(src)) file.copy(src, corrected_dir, overwrite = TRUE)
}

message("DONE.\nLog: ", log_file, "\nFigures: ", fig_dir, "\nPosthoc tables: ", log_dir)

```

### `7b.tgc.tms.pcs.dapc.heatmaps.R`

```r
#!/usr/bin/env Rscript
############################################################
# TreeGeneClimate (TGC) — TMS
# Step 7b (UPDATED): PCA panel (supp) + DAPC panel (Figure 5) + DAPC loadings TSVs
#
# CHANGES (requested):
# - OMIT heatmaps entirely (to be done later as Step 8b).
# - Save PCA as ONE panel (a–f) with shared legend per row:
#     top row = breeding (CpG, CHG, CHH) + legend at right of c)
#     bottom  = natural  (CpG, CHG, CHH) + legend at right of f)
# - Save DAPC biplots as ONE panel (a–f) with same layout/legend behavior.
# - Remove the bold ggplot titles (no duplicated titles). Panel letters are drawn inside plots.
# - Fix epimarker IDs: use genomic labels chr:start (chr:pos) instead of V123 / V1 etc.
# - Fix loadings tables: output loc/chr/pos (pos=start) + DF + loading.
#
# INPUT:
#   RESULTS/TMS/RANALYSIS/METHYLKIT_OBJECTS/
#     methylBase_<cohort>_<context>_cov5_50_mpg*_mef0.05.rds
#
# OUTPUT:
#   RESULTS/TMS/RANALYSIS/FIGURES/FIG5_FIG6/
#     SUPP_PCA_panel_a-f.tiff
#     Figure5_DAPC_panel_a-f.tiff
#
#   RESULTS/TMS/RANALYSIS/TABLES/dapc_loadings/
#     TMS_DAPC_loadings_<cohort>_<context>_ALL.tsv
############################################################

suppressPackageStartupMessages({
  library(methylKit)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(scales)
  library(grid)
  library(adegenet)   # dapc()
  library(ggrepel)    # non-overlapping text labels for loading vectors
  library(tidyr)
  library(patchwork)
})

options(stringsAsFactors = FALSE)
set.seed(1)   # ensures reproducible DAPC (which uses random SVD internally)

# ==============================================================================
# 1) PATHS
# ==============================================================================
# === USER CONFIGURATION ===
# Set TGC_PROJECT_ROOT as an environment variable, or edit the fallback path below
PROJECT_ROOT <- Sys.getenv("TGC_PROJECT_ROOT",
  unset = "/path/to/your/project")
# ===========================
rds_dir <- file.path(PROJECT_ROOT, "RESULTS/TMS/RANALYSIS/METHYLKIT_OBJECTS")
fig_dir <- file.path(PROJECT_ROOT, "RESULTS/TMS/RANALYSIS/FIGURES/FIG5_FIG6")
tab_dir <- file.path(PROJECT_ROOT, "RESULTS/TMS/RANALYSIS/TABLES/dapc_loadings")

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
# Named vectors: consistent group-to-color mapping used across PCA and DAPC panels
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

# Select the MEF-filtered RDS with the highest mpg (most sites) to maximise power
pick_rds_mef <- function(cohort, ctx) {
  pat <- sprintf("^methylBase_%s_%s_cov5_50_mpg[0-9]+_mef0\\.05\\.rds$",
                 tolower(cohort), tolower(ctx))
  files <- list.files(rds_dir, pattern = pat, full.names = TRUE)
  if (length(files) == 0) stop("No MEF RDS found for ", cohort, " / ", ctx, " under ", rds_dir)

  # extract the mpg number and pick the highest (largest site set)
  mpg_num <- as.integer(sub(".*_mpg([0-9]+)_mef0\\.05\\.rds$", "\\1", basename(files)))
  files[which.max(mpg_num)]
}

# Convert methylBase to methylation % matrix (sites x samples) + site info table
# Row names of the returned matrix are set to genomic "loc" IDs (chr:start-end)
# so that DAPC variable names are interpretable rather than V1, V2, ...
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

  # IMPORTANT FIX:
  # Set rownames to loc so DAPC loadings are loc (not V1, V2, ...)
  rownames(perc_mat) <- site_df$loc
  colnames(perc_mat) <- sample_ids

  list(mat_sites_x_samples = perc_mat, site_df = site_df)
}

# panel label inside plot — uses ggplot title to match ECS PCA format
add_panel_label <- function(p, lab) {
  p + labs(title = lab) +
    theme(plot.title = element_text(hjust = 0))
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
# inflate artificial PCs and destabilise DAPC).
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

# DAPC (Discriminant Analysis of Principal Components) via adegenet.
# n.pca is chosen conservatively to avoid over-fitting (retaining too many PCs
# inflates classification accuracy without biological meaning).
run_dapc <- function(X_samples_x_sites, group_factor, max_na_frac = 0.20) {
  group_factor <- factor(group_factor)
  X2 <- filter_and_impute_sites(X_samples_x_sites, max_na_frac = max_na_frac)$X

  # Cap n.pca below n_samples - n_groups to avoid rank deficiency
  n.pca <- max(10L, min(80L, nrow(X2) - nlevels(group_factor)))
  # Number of discriminant functions = min(2, n_groups - 1) for biplot display
  n.da  <- max(1L, min(2L, nlevels(group_factor) - 1L))

  fit <- adegenet::dapc(x = X2, grp = group_factor, n.pca = n.pca, n.da = n.da, var.contrib = TRUE)

  coords <- as.data.frame(fit$ind.coord)
  colnames(coords) <- paste0("DF", seq_len(ncol(coords)))

  # var.contr: contribution of each locus to each discriminant function
  vc <- as.data.frame(fit$var.contr)
  if (nrow(vc) == 0) stop("DAPC var.contr empty.")
  vc$loc <- rownames(vc)  # should now be chr:start-end

  num_cols <- which(vapply(vc, is.numeric, logical(1)))
  if (length(num_cols) < 2) stop("DAPC loadings missing 2 numeric DF columns.")
  df1_col <- names(vc)[num_cols[1]]
  df2_col <- names(vc)[num_cols[2]]

  # Reshape loadings to long format for easier downstream handling
  load_long <- bind_rows(
    vc %>% transmute(loc = loc, DF = "DF1", loading = .data[[df1_col]]),
    vc %>% transmute(loc = loc, DF = "DF2", loading = .data[[df2_col]])
  )

  list(fit = fit, coords = coords, load_long = load_long, X_used = X2)
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

# Build DAPC biplot: individual scores + top loading vectors for DF1 and DF2.
# Vectors are scaled into individual score space for readability.
build_dapc_biplot <- function(dapc_res, site_df_used, group_vec, palette_named) {
  coords <- dapc_res$coords
  coords$Group <- as.character(group_vec)

  # key: loc -> chrpos (chr:start) (requested)
  key <- site_df_used %>%
    dplyr::select(loc, chr, pos, chrpos)

  # Pivot loadings wide to get DF1 and DF2 as columns for each locus
  ld_w <- dapc_res$load_long %>%
    filter(DF %in% c("DF1","DF2")) %>%
    tidyr::pivot_wider(names_from = DF, values_from = loading)

  # Top 5 loci by |loading| on DF1 and DF2 respectively; deduplicate overlaps
  top1 <- ld_w %>% arrange(desc(abs(DF1))) %>% slice_head(n = 5)
  top2 <- ld_w %>% arrange(desc(abs(DF2))) %>% slice_head(n = 5)

  topm <- bind_rows(top1, top2) %>%
    distinct(loc, .keep_all = TRUE) %>%
    left_join(key, by = "loc") %>%
    mutate(label = ifelse(is.na(chrpos), loc, chrpos))

  # scale loading vectors into individual space
  ind_r <- sqrt(coords$DF1^2 + coords$DF2^2)
  max_ind_r <- max(ind_r, na.rm = TRUE)
  load_r <- sqrt(topm$DF1^2 + topm$DF2^2)
  max_load_r <- max(load_r, na.rm = TRUE)
  # Scale factor: vectors fill 85% of individual score radius
  scale_factor <- ifelse(is.finite(max_load_r) && max_load_r > 0, (0.85 * max_ind_r) / max_load_r, 1)
  topm <- topm %>% mutate(DF1s = DF1 * scale_factor, DF2s = DF2 * scale_factor)

  pal <- subset_palette(palette_named, coords$Group)

  ggplot(coords, aes(x = DF1, y = DF2, color = Group)) +
    geom_point(size = 1.8, alpha = 0.9) +
    scale_color_manual(values = pal) +
    labs(x = "DF1", y = "DF2") +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(),
      legend.title = element_blank()
    ) +
    # Loading vectors as arrows from origin
    geom_segment(
      data = topm,
      aes(x = 0, y = 0, xend = DF1s, yend = DF2s),
      inherit.aes = FALSE,
      arrow = grid::arrow(length = unit(0.18, "cm")),
      color = "black",
      linewidth = 0.6
    ) +
    # Non-overlapping genomic labels (chr:start) at arrow tips
    ggrepel::geom_text_repel(
      data = topm,
      aes(x = DF1s, y = DF2s, label = label),
      inherit.aes = FALSE,
      size = 3.0,
      min.segment.length = 0.05,
      box.padding = 0.25,
      point.padding = 0.2,
      segment.size = 0.3
    )
}

# Save panel in four publication formats; TIFF uses lzw compression for journal submission
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

  # site_df_used must align to columns of X (loc)
  site_df_used0 <- site_df[match(colnames(X), site_df$loc), ]
  if (any(is.na(site_df_used0$loc))) stop("Could not align site_df with matrix columns (loc).")

  max_na_frac <- max_na_frac_by_ctx(ctx)

  # PCA on the full filtered site set (unsupervised; used as a supplementary figure)
  pca <- run_pca(X, max_na_frac = max_na_frac)
  pve <- pca$pve
  scores <- as.data.frame(pca$pr$x[, 1:2, drop = FALSE])
  colnames(scores) <- c("PC1","PC2")
  scores$Group <- as.character(group_factor)
  p_pca <- make_pca_scatter(scores, pve, palette_named)

  # DAPC: supervised ordination maximising between-group variance
  dapc <- run_dapc(X, group_factor, max_na_frac = max_na_frac)

  # Align site_df to DAPC-used columns (sites can differ from PCA after SD filter)
  site_df_used <- site_df[match(colnames(dapc$X_used), site_df$loc), ]
  if (any(is.na(site_df_used$loc))) stop("Could not align site_df with DAPC markers (loc).")

  p_dapc <- build_dapc_biplot(dapc, site_df_used, group_factor, palette_named)

  # Write per-locus loading table with genomic coordinates (chr, pos=start)
  load_all <- dapc$load_long %>%
    left_join(site_df_used %>% dplyr::select(loc, chr, pos), by = "loc") %>%
    dplyr::relocate(chr, pos, DF, loading, .after = loc)

  out_tsv <- file.path(tab_dir, sprintf("TMS_DAPC_loadings_%s_%s_ALL.tsv",
                                        tolower(cohort), tolower(ctx)))
  write.table(load_all, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

  list(
    cohort = cohort,
    ctx = ctx,
    pca_plot = p_pca,
    dapc_plot = p_dapc,
    loadings_tsv = out_tsv
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
# 6) PCA PANEL (a–f) with shared legend per row
# ==============================================================================
# Top row: breeding (a=CpG, b=CHG, c=CHH); bottom row: natural (d=CpG, e=CHG, f=CHH)
pca_a <- add_panel_label(res_b_cpg$pca_plot, "a)")
pca_b <- add_panel_label(res_b_chg$pca_plot, "b)")
pca_c <- add_panel_label(res_b_chh$pca_plot, "c)")

pca_d <- add_panel_label(res_n_cpg$pca_plot, "d)")
pca_e <- add_panel_label(res_n_chg$pca_plot, "e)")
pca_f <- add_panel_label(res_n_chh$pca_plot, "f)")

# Collect legend per row so breeding and natural families/populations have separate legends
pca_row1 <- (pca_a | pca_b | pca_c) + plot_layout(guides = "collect") & theme(legend.position = "right", legend.justification = "top")
pca_row2 <- (pca_d | pca_e | pca_f) + plot_layout(guides = "collect") & theme(legend.position = "right", legend.justification = "top")

pca_panel <- pca_row1 / pca_row2

out_pca_panel <- file.path(fig_dir, "SUPP_PCA_panel_a-f.tiff")
save_panel_tiff(pca_panel, out_pca_panel, w_cm = 34, h_cm = 26, dpi = 600)

# ==============================================================================
# 7) DAPC PANEL (Figure 5 a–f) with shared legend per row
# ==============================================================================
dapc_a <- add_panel_label(res_b_cpg$dapc_plot, "a)")
dapc_b <- add_panel_label(res_b_chg$dapc_plot, "b)")
dapc_c <- add_panel_label(res_b_chh$dapc_plot, "c)")

dapc_d <- add_panel_label(res_n_cpg$dapc_plot, "d)")
dapc_e <- add_panel_label(res_n_chg$dapc_plot, "e)")
dapc_f <- add_panel_label(res_n_chh$dapc_plot, "f)")

dapc_row1 <- (dapc_a | dapc_b | dapc_c) + plot_layout(guides = "collect") & theme(legend.position = "right", legend.justification = "top")
dapc_row2 <- (dapc_d | dapc_e | dapc_f) + plot_layout(guides = "collect") & theme(legend.position = "right", legend.justification = "top")

dapc_panel <- dapc_row1 / dapc_row2

out_dapc_panel <- file.path(fig_dir, "Figure5_DAPC_panel_a-f.tiff")
save_panel_tiff(dapc_panel, out_dapc_panel, w_cm = 34, h_cm = 26, dpi = 600)

# ==============================================================================
# 8) FINISH
# ==============================================================================
cat("\nDONE Step 7b (PCA + DAPC only).\n\n")
cat("Panels saved in:\n  ", fig_dir, "\n\n", sep = "")
cat("  - ", basename(out_pca_panel), "\n", sep = "")
cat("  - ", basename(out_dapc_panel), "\n\n", sep = "")
cat("Loadings TSVs saved in:\n  ", tab_dir, "\n\n", sep = "")
cat("  - ", res_b_cpg$loadings_tsv, "\n", sep = "")
cat("  - ", res_b_chg$loadings_tsv, "\n", sep = "")
cat("  - ", res_b_chh$loadings_tsv, "\n", sep = "")
cat("  - ", res_n_cpg$loadings_tsv, "\n", sep = "")
cat("  - ", res_n_chg$loadings_tsv, "\n", sep = "")
cat("  - ", res_n_chh$loadings_tsv, "\n", sep = "")
sessionInfo()

```

---

### `8b.tgc.tms.heatmaps.R`

```r
#!/usr/bin/env Rscript
############################################################
# TreeGeneClimate (TGC) — TMS
# Step 8b: Heatmaps (per cohort) using cov+MEF filtered epiloci
#
# INPUT (same as 7b):
#   RESULTS/TMS/RANALYSIS/METHYLKIT_OBJECTS/
#     methylBase_<cohort>_<context>_cov5_50_mpg*_mef0.05.rds
#
# METHOD:
#   - For each cohort and each context (CpG/CHG/CHH):
#       * Convert methylBase to methylation % matrix (samples x loci)
#       * Filter loci by missingness (context-aware)
#       * Locus-by-locus Kruskal–Wallis test: %meth ~ Group (Family / Stand)
#       * BH adjust within context
#   - Select loci for heatmap by BH-adjusted p-value:
#       * Default: top 200 loci per cohort, approximately balanced across contexts
#   - Heatmap values are methylation % (0–100)
#   - Clustering requires a numeric matrix without too many NAs:
#       * Remaining NAs are imputed by locus mean (ONLY for heatmap/clustering)
#
# OUTPUT:
#   RESULTS/TMS/RANALYSIS/FIGURES/HEATMAPS_8B/
#     Figure7a_HEATMAP_breeding_top200_KW_BH.tiff
#     Figure7b_HEATMAP_natural_top200_KW_BH.tiff
#
#   RESULTS/TMS/RANALYSIS/TABLES/heatmap_markers_8B/
#     TMS_8B_locus_tests_breeding_<ctx>.tsv
#     TMS_8B_locus_tests_natural_<ctx>.tsv
#     TMS_8B_selected_markers_breeding_top200.tsv
#     TMS_8B_selected_markers_natural_top200.tsv
############################################################

suppressPackageStartupMessages({
  library(methylKit)
  library(dplyr)
  library(tidyr)
  library(tibble)

  library(ComplexHeatmap)   # publication-quality heatmaps with annotations
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

fig_dir <- file.path(PROJECT_ROOT, "RESULTS/TMS/RANALYSIS/FIGURES/HEATMAPS_8B")
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

# Okabe–Ito (colorblind-friendly) for context
context_cols <- c(
  CpG = "#009E73",  # bluish green
  CHG = "#D55E00",  # vermillion
  CHH = "#0072B2"   # blue
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

# Select the MEF-filtered RDS with the highest mpg (most sites) to maximise power
pick_rds_mef <- function(cohort, ctx) {
  pat <- sprintf("^methylBase_%s_%s_cov5_50_mpg[0-9]+_mef0\\.05\\.rds$",
                 tolower(cohort), tolower(ctx))
  files <- list.files(rds_dir, pattern = pat, full.names = TRUE)
  if (length(files) == 0) stop("No MEF RDS found for ", cohort, " / ", ctx, " under ", rds_dir)

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
# KW is used instead of ANOVA because per-locus methylation % is bounded (0–100)
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
# Complementary to balanced selection; can reveal context dominance in the signal.
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

# ---- Save a single heatmap TIFF (20x20 cm journal style by default)
# ComplexHeatmap requires a draw() call inside the graphics device;
# ggsave is not compatible, so each format is opened and closed manually.
save_heatmap_tiff <- function(mat_samples_x_markers,
                              group_vec, group_palette,
                              context_vec,
                              group_legend_title = "Group",
                              out_file, title_label = "a)",
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
              gp = gpar(fontsize = 16))
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

    res_list[[ctx]] <- list(
      X_heat = Xh,
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
  df_all <- NULL

  for (ctx in names(res_list)) {
    Xc <- res_list[[ctx]]$X_heat[samp_common, , drop = FALSE]
    X_all <- if (is.null(X_all)) Xc else cbind(X_all, Xc)

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
    X_all     = X_all
  )
}

# ==============================================================================
# 5) RUN + SAVE HEATMAPS (2 per cohort; NO PANELS)
# ==============================================================================
topN <- 150L   # number of loci shown per heatmap (50 per context when balanced)

# ---- BREEDING
res_b <- process_cohort_for_heatmap("BREEDING", map_file_breeding, colors.17, top_n = topN)

out_b_bal <- file.path(fig_dir, sprintf("Figure7a_HEATMAP_breeding_top%d_BALANCED_KW_BH.tiff", topN))
out_b_all <- file.path(fig_dir, sprintf("Figure7a2_HEATMAP_breeding_top%d_OVERALL_KW_BH.tiff", topN))

save_heatmap_tiff(res_b$X_sel_bal, res_b$group_vec, colors.17, res_b$ctx_vec_bal,
                  group_legend_title = "Family",
                  out_file = out_b_bal, title_label = "a)",
                  w_cm = 20, h_cm = 20, dpi = 600)

save_heatmap_tiff(res_b$X_sel_all, res_b$group_vec, colors.17, res_b$ctx_vec_all,
                  group_legend_title = "Family",
                  out_file = out_b_all, title_label = "a)",
                  w_cm = 20, h_cm = 20, dpi = 600)

# ---- NATURAL
res_n <- process_cohort_for_heatmap("NATURAL", map_file_natural, colors.25, top_n = topN)

out_n_bal <- file.path(fig_dir, sprintf("Figure7b_HEATMAP_natural_top%d_BALANCED_KW_BH.tiff", topN))
out_n_all <- file.path(fig_dir, sprintf("Figure7b2_HEATMAP_natural_top%d_OVERALL_KW_BH.tiff", topN))

save_heatmap_tiff(res_n$X_sel_bal, res_n$group_vec, colors.25, res_n$ctx_vec_bal,
                  group_legend_title = "Natural stand",
                  out_file = out_n_bal, title_label = "b)",
                  w_cm = 20, h_cm = 20, dpi = 600)

save_heatmap_tiff(res_n$X_sel_all, res_n$group_vec, colors.25, res_n$ctx_vec_all,
                  group_legend_title = "Natural stand",
                  out_file = out_n_all, title_label = "b)",
                  w_cm = 20, h_cm = 20, dpi = 600)

# ==============================================================================
# 6) NATURAL: FORMAL SVMPs ONLY (padj < 0.05) + COMBINED PANEL (a + b)
# ==============================================================================
# SVMPs (Stochastic Variation in Methylation Patterns) = loci with significant
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

  # ---- Figure 7b: natural formal SVMPs standalone
  out_n_sig <- file.path(fig_dir, "Figure7b_HEATMAP_natural_formal_SVMPs_KW_BH.tiff")
  save_heatmap_tiff(X_sig, res_n$group_vec, colors.25, ctx_vec_sig,
                    group_legend_title = "Natural stand",
                    out_file         = out_n_sig,
                    title_label      = "b)",
                    cluster_cols     = TRUE,
                    show_annot_names = FALSE,   # annotation names omitted for the narrow panel
                    w_cm = 12, h_cm = 20, dpi = 600)

  # ---- Combined panel (a + b)
  # Breeding balanced (150 markers) gets 70% of width;
  # Natural formal SVMPs (6 markers) gets 30% — not proportional, readable.
  panel_w_cm <- 28
  panel_h_cm <- 20
  frac_a     <- 0.70
  label_gp   <- gpar(fontsize = 16)   # no bold

  # Build heatmap objects for the combined panel (legends merged per sub-panel)
  ht_a_panel <- build_heatmap_object(
    mat_samples_x_markers = res_b$X_sel_bal,
    group_vec             = res_b$group_vec,
    group_palette         = colors.17,
    context_vec           = res_b$ctx_vec_bal,
    group_legend_title    = "Breeding family",
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
    grid.text("a)",
              x    = unit(0.35, "cm"),
              y    = unit(1, "npc") - unit(0.35, "cm"),
              just = c("left", "top"),
              gp   = label_gp)
    popViewport()

    # Panel b — natural formal SVMPs
    pushViewport(viewport(x = frac_a, y = 0, width = 1 - frac_a, height = 1,
                          just = c("left", "bottom")))
    draw(ht_b_panel,
         heatmap_legend_side    = "right",
         annotation_legend_side = "right",
         merge_legend           = TRUE,
         newpage                = FALSE)
    grid.text("b)",
              x    = unit(0.35, "cm"),
              y    = unit(1, "npc") - unit(0.35, "cm"),
              just = c("left", "top"),
              gp   = label_gp)
    popViewport()
  }

  panel_tiff <- file.path(fig_dir, "Figure7_PANEL_ab_BALANCED_SVMPs.tiff")
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

  cat("Figure 7b standalone (formal SVMPs): ", out_n_sig, "\n")
  cat("Figure 7 combined panel (a + b):     ", panel_tiff, "\n")

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
#   1) Build a genomic distance matrix from ECS imputed GDS genotypes (IBS distance = 1 - IBS)
#   2) Build an epigenomic distance matrix from TMS methylBase (Euclidean on filtered+imputed % methylation)
#   3) Align samples robustly (trimws + tolower), save diagnostics for mismatches
#   4) Compare distance matrices:
#       - Mantel test
#       - Procrustes / protest (PCoA ordinations)
#       - RV coefficient
#   5) Save per-comparison tables + combined 3x2 Procrustes panel + aggregated summary
#
# INPUTS
# - ECS:
#   RESULTS/ECS/RANALYSIS/RDATA/breeding.imputed.snp.gds
#   RESULTS/ECS/RANALYSIS/RDATA/natural.imputed.snp.gds
# - TMS:
#   RESULTS/TMS/RANALYSIS/METHYLKIT_OBJECTS/
#     methylBase_<cohort>_<context>_cov5_50_mpg*_mef0.05.rds
# - Metadata:
#   DATA/METADATA/breeding_sample2family.txt
#   DATA/METADATA/natural_sample2pop.txt
#
# OUTPUTS
# - /path/to/your/project
#     FIGURES/11ab/
#     TABLES/11ab/
#
# NOTES
# - PERMUTATIONS fixed at 9999.
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
PROJECT_ROOT <- Sys.getenv("TGC_PROJECT_ROOT",
  unset = "/path/to/your/project")
# ===========================

# ECS (GDS produced by 8a)
ECS_RDATA_DIR <- file.path(PROJECT_ROOT, "RESULTS/ECS/RANALYSIS/RDATA")
GDS_BREED <- file.path(ECS_RDATA_DIR, "breeding.imputed.snp.gds")
GDS_NATUR <- file.path(ECS_RDATA_DIR, "natural.imputed.snp.gds")

# TMS methylKit objects
TMS_RDS_DIR <- file.path(PROJECT_ROOT, "RESULTS/TMS/RANALYSIS/METHYLKIT_OBJECTS")

# Metadata maps
META_DIR <- file.path(PROJECT_ROOT, "DATA/METADATA")
MAP_BREED <- file.path(META_DIR, "breeding_sample2family.txt")
MAP_NATUR <- file.path(META_DIR, "natural_sample2pop.txt")

# JOINT outputs
JOINT_ROOT <- "/path/to/your/project"
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
# Standard SNP QC thresholds used throughout the ECS pipeline
ECS_MAF_MIN  <- 0.05
ECS_MISS_MAX <- 0.10

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
  pat <- sprintf("^methylBase_%s_%s_cov5_50_mpg[0-9]+_mef0\\.05\\.rds$",
                 tolower(cohort), tolower(ctx))
  files <- list.files(TMS_RDS_DIR, pattern = pat, full.names = TRUE)
  if (!length(files)) stop("No TMS MEF RDS found for ", cohort, " / ", ctx, " under ", TMS_RDS_DIR)
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
  keep <- is.finite(maf) & is.finite(miss) & maf >= 0.05 & miss <= 0.10

  if (sum(keep) == 0) stop("No ECS SNPs passed MAF/missingness filters.")

  ibs <- snpgdsIBS(
    g,
    snp.id = stat$snp.id[keep],
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

# ---- TMS epigenomic distance: Euclidean on %methylation (sites filtered+imputed) ----
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

filter_and_impute_tms <- function(X_samples_x_sites, max_na_frac) {
  # Remove sites that are missing in too many samples
  na_frac <- colMeans(is.na(X_samples_x_sites))
  keep1 <- na_frac <= max_na_frac
  X2 <- X_samples_x_sites[, keep1, drop = FALSE]
  if (ncol(X2) < 10) stop("Too few TMS sites after NA filter (n=", ncol(X2), ").")

  # Mean-impute remaining NAs within each site (column-wise)
  for (j in seq_len(ncol(X2))) {
    idx <- is.na(X2[, j])
    if (any(idx)) {
      mu <- mean(X2[, j], na.rm = TRUE)
      if (!is.finite(mu)) mu <- 0
      X2[idx, j] <- mu
    }
  }

  # Drop invariant sites — they contribute nothing to Euclidean distance
  sds <- apply(X2, 2, sd)
  keep2 <- is.finite(sds) & sds > 0
  X3 <- X2[, keep2, drop = FALSE]
  if (ncol(X3) < 10) stop("Too few variable TMS sites after SD filter (n=", ncol(X3), ").")

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
  X2 <- filter_and_impute_tms(X, max_na_frac = max_na)
  # Euclidean distance across all filtered CpG/CHG/CHH sites (high-dimensional space)
  D <- dist(X2, method = "euclidean")
  list(dist = D, ids = rownames(X2), n_sites = ncol(X2))
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
make_procrustes_plot <- function(X, Y, title_label = "a)") {
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
    labs(x = "PC1", y = "PC2") +
    annotate("text", x = -Inf, y = Inf, label = title_label,
             hjust = -0.2, vjust = 1.3, size = 7) +
    theme_minimal(base_size = 18) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.direction = "horizontal",
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 12),
      legend.text = element_text(size = 16),
      legend.key.size = unit(4, "lines")
    )

  list(plot = p, pro = pro)
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

# Panel labels follow publication convention (a–f, top-left to bottom-right)
panel_labels <- c("a)", "b)", "c)", "d)", "e)", "f)")
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
# Constrain y-axis for NATURAL panels (d–f) which have tighter ordination spread
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

out_panel <- file.path(OUT_FIG, "Figure11_ECS_TMS_Procrustes_panel.tiff")

# Save in all required publication formats
tiff(out_panel, width = 34, height = 24, units = "cm", res = 600, compression = "lzw")
print(combined_panel)
dev.off()
ggsave(sub("\\.tiff$", ".pdf", out_panel), plot = combined_panel, width = 34, height = 24, units = "cm")
ggsave(sub("\\.tiff$", ".eps", out_panel), plot = combined_panel, width = 34, height = 24, units = "cm", device = cairo_ps)
ggsave(sub("\\.tiff$", ".png", out_panel), plot = combined_panel, width = 34, height = 24, units = "cm", dpi = 150, device = "png")

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
# Step 12ab0: Prepare meQTL5 inputs (OLD-style inputs)
#
# KEY DIFFERENCES vs 15ab0 (meQTL4):
#   - GDS     : non-imputed  breeding.snp.gds / natural.snp.gds
#   - GRM     : GCTA GRM     breeding_grm_gcta.rds / natural_grm_gcta.rds
#   - Methylation: unfiltered {cohort}_{ctx}_methylkit.rds (no MEF filter)
#   - PCs     : 10 for BOTH cohorts, derived from GCTA GRM eigendecomposition
#   - NO quantile normalisation of M-values
#   - Group assignment (Family / Population) saved per panel for HWE imputation
#   - LOCO GRMs for BREEDING computed from non-imputed GDS (VanRaden formula)
#   - Full GCTA GRM saved for NATURAL (used by GENESIS5)
#
# OUTPUTS (under RESULTS/JOINT/MQTL5/INPUTS/<COHORT>/<CONTEXT>/)
#   methylation_mvalues_matrix.rds / .tsv
#   methylation_site_annot.rds / .tsv
#   shared_sample_ids.rds / .tsv
#   pcs_shared.rds / .tsv            (PC1..PC10, from GRM eigen)
#   grm_shared.rds                   (GCTA GRM subset to shared samples)
#   sample_groups.rds / .tsv         (group assignments for HWE imputation)
#   snp_gds_sample_ids.rds / .tsv
#   snp_variant_annot.rds / .tsv
# Per-cohort LOCO GRMs (BREEDING only):
#   INPUTS/BREEDING/loco_grm/loco_grm_excl_<chr>.rds
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(SNPRelate)
  library(gdsfmt)
  library(methylKit)
})

options(stringsAsFactors = FALSE)
options(datatable.fread.datatable = TRUE)

############################################################
# 1) SETTINGS
############################################################

# 10 PCs used as fixed-effect covariates in both MatrixEQTL and GENESIS models
N_PCS    <- 10L

# Pseudo-count for logit transformation: M = log((pct + eps) / (100 - pct + eps))
# eps = 0.5 avoids log(0) at boundary values without over-shrinking
EPSILON_M <- 0.5

COHORTS  <- c("BREEDING", "NATURAL")
CONTEXTS <- c("CpG", "CHG", "CHH")

############################################################
# 2) PATHS
############################################################

# === USER CONFIGURATION ===
# Set TGC_PROJECT_ROOT as an environment variable, or edit the fallback path below
PROJECT_ROOT <- Sys.getenv("TGC_PROJECT_ROOT",
  unset = "/path/to/your/project")
# ===========================
OUTROOT      <- file.path(PROJECT_ROOT, "RESULTS", "JOINT", "MQTL5")

LOGDIR   <- file.path(OUTROOT, "LOGS")
SUMDIR   <- file.path(OUTROOT, "SUMMARIES")
DBGDIR   <- file.path(OUTROOT, "DEBUG")
INPUTDIR <- file.path(OUTROOT, "INPUTS")

for (d in c(OUTROOT, LOGDIR, SUMDIR, DBGDIR, INPUTDIR))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

LOGFILE      <- file.path(LOGDIR, "step12ab0.log")
SUMMARY_FILE <- file.path(SUMDIR, "step12ab0_summary.tsv")
SESSION_FILE <- file.path(SUMDIR, "step12ab0_sessionInfo.txt")
OVERLAP_FILE <- file.path(DBGDIR, "step12ab0_overlap_summary.tsv")

# Clear any logs from a previous run before starting
for (f in c(LOGFILE, SUMMARY_FILE, OVERLAP_FILE))
  if (file.exists(f)) file.remove(f)

# --- OLD-STYLE input paths ---
RDATA_DIR   <- "/path/to/your/project"
TMS_RDS_DIR <- "/path/to/your/project"

# Non-imputed GDS files — imputed GDS is used only in step 11ab
GDS_FILES <- list(
  BREEDING = file.path(RDATA_DIR, "breeding.snp.gds"),
  NATURAL  = file.path(RDATA_DIR, "natural.snp.gds")
)
# GCTA GRMs pre-computed from the full SNP set
GRM_FILES <- list(
  BREEDING = file.path(RDATA_DIR, "breeding_grm_gcta.rds"),
  NATURAL  = file.path(RDATA_DIR, "natural_grm_gcta.rds")
)
# Sample metadata (IID, Family/Population) for group-aware HWE imputation
ANNOT_FILES <- list(
  BREEDING = file.path(RDATA_DIR, "breeding_sample_annotation.rds"),
  NATURAL  = file.path(RDATA_DIR, "natural_sample_annotation.rds")
)

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

# Canonical sample ID normalisation: remove trailing ".0" (Excel artefact),
# leading "X" (R's make.names prefix), replace "-" with "_"
normalize_id <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("\\.0$", "", x)
  x <- gsub("^X",   "", x)
  x <- gsub("-",    "_", x)
  x
}

# Convert a methylKit methylBase object to a logit M-value matrix (samples x sites).
# Uses non-quantile-normalised values; sites with zero variance are dropped.
methylkit_to_mvalues <- function(meth_path, eps = 0.5) {
  if (!file.exists(meth_path)) stop("Missing methylKit RDS: ", meth_path)
  mb   <- readRDS(meth_path)
  d    <- methylKit::getData(mb)
  cs_cols <- grep("^numCs[0-9]+$", colnames(d), value = TRUE)
  ts_cols <- grep("^numTs[0-9]+$", colnames(d), value = TRUE)
  stopifnot(length(cs_cols) > 0, length(cs_cols) == length(ts_cols))
  sample_ids <- mb@sample.ids
  stopifnot(length(sample_ids) == length(cs_cols))
  site_annot <- data.table(
    chr   = as.character(d$chr),
    start = as.integer(d$start),
    end   = as.integer(d$end)
  )
  site_annot[, pos       := start]
  site_annot[, methyl_loc := paste0(chr, ":", start, "-", end)]
  # Compute percentage methylation: 100 * C / (C + T)
  perc_mat <- vapply(seq_along(cs_cols), function(i) {
    cs  <- d[[cs_cols[i]]]; ts <- d[[ts_cols[i]]]
    pct <- 100 * cs / (cs + ts)
    pct[is.nan(pct)] <- NA_real_
    pct
  }, numeric(nrow(d)))
  colnames(perc_mat) <- sample_ids
  # Remove sites with zero variance or all-NA
  keep_var <- apply(perc_mat, 1, function(x) sum(!is.na(x)) > 1 && sd(x, na.rm = TRUE) > 0)
  perc_mat  <- perc_mat[keep_var, , drop = FALSE]
  site_annot <- site_annot[keep_var]
  # Logit transform: M-values are approximately normally distributed and
  # better suited for linear models than beta/percentage values
  m_mat <- log((perc_mat + eps) / (100 - perc_mat + eps))
  rownames(m_mat) <- site_annot$methyl_loc
  m_mat <- t(m_mat)   # samples x sites
  rownames(m_mat) <- sample_ids
  list(m_matrix = m_mat, site_annot = site_annot)
}

# Derive population structure PCs by eigendecomposition of the GRM.
# This avoids re-running GCTA --pca and keeps PCs on the same scale as the GRM.
grm_pcs <- function(G, sample_norm, n_pcs) {
  # G: matrix with rownames (raw or normalised IDs)
  # Returns data.table: sample_id_raw, sample_id, PC1..PCn_pcs
  ev  <- eigen(G, symmetric = TRUE)
  pcs <- ev$vectors[, seq_len(min(n_pcs, ncol(ev$vectors))), drop = FALSE]
  colnames(pcs) <- paste0("PC", seq_len(ncol(pcs)))
  rownames(pcs) <- sample_norm
  dt  <- as.data.table(pcs, keep.rownames = "sample_id")
  dt
}

############################################################
# 4) LOCO GRM — BREEDING (VanRaden formula, non-imputed GDS)
############################################################

# LOCO (leave-one-chromosome-out) GRMs are required by the GENESIS LMM to avoid
# proximal contamination: the chromosome being tested is excluded from the kinship
# matrix, preventing the GRM from inadvertently capturing the tested signal.
# Computed only for BREEDING because relatedness structure is stronger in family material.
compute_loco_grms <- function(gds_path, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  log_msg("Computing LOCO GRMs from: ", gds_path)
  gds <- snpgdsOpen(gds_path, readonly = TRUE)
  on.exit(snpgdsClose(gds), add = TRUE)

  sample_ids <- read.gdsn(index.gdsn(gds, "sample.id"))
  all_chr    <- as.character(read.gdsn(index.gdsn(gds, "snp.chromosome")))
  chrs       <- sort(unique(all_chr))
  log_msg("  Samples: ", length(sample_ids), " | Chromosomes: ", length(chrs))

  log_msg("  Reading all genotypes...")
  gl   <- snpgdsGetGeno(gds, snpfirstdim = TRUE, with.id = TRUE)
  geno <- gl$genotype   # SNPs x samples
  sample_norm <- normalize_id(gl$sample.id)
  log_msg("  Genotype matrix: ", nrow(geno), " SNPs x ", ncol(geno), " samples")

  # VanRaden (2008) GRM: G = Z Z' / (2 sum p_i (1-p_i))
  # where Z = (g - 2p) is the mean-centered genotype matrix
  vanraden_grm <- function(geno_sub) {
    p    <- rowSums(geno_sub, na.rm = TRUE) / (2 * rowSums(!is.na(geno_sub)))
    poly <- p > 0 & p < 1 & !is.na(p)  # exclude fixed/monomorphic SNPs
    if (sum(poly) < 2L) return(NULL)
    g2 <- geno_sub[poly, , drop = FALSE]
    p2 <- p[poly]
    Xc <- g2 - 2 * p2   # center by 2*allele frequency
    Xc[is.na(Xc)] <- 0  # set missing genotypes to zero after centering
    denom <- 2 * sum(p2 * (1 - p2))
    if (!is.finite(denom) || denom <= 0) return(NULL)
    G <- (t(Xc) %*% Xc) / denom
    rownames(G) <- colnames(G) <- sample_norm
    G
  }

  n_done <- 0L
  for (excl_chr in chrs) {
    out_path <- file.path(out_dir, paste0("loco_grm_excl_", excl_chr, ".rds"))
    if (file.exists(out_path)) {
      log_msg("  LOCO chr ", excl_chr, ": already exists — skipping")
      n_done <- n_done + 1L; next
    }
    # Build GRM from all SNPs except those on the chromosome being tested
    mask        <- all_chr != excl_chr
    n_snps_used <- sum(mask)
    if (n_snps_used < 10L) {
      log_msg("  LOCO chr ", excl_chr, ": only ", n_snps_used, " SNPs — skipping"); next
    }
    G <- vanraden_grm(geno[mask, , drop = FALSE])
    if (is.null(G)) {
      log_msg("  LOCO chr ", excl_chr, ": GRM failed — skipping"); next
    }
    saveRDS(G, out_path)
    log_msg("  LOCO chr ", excl_chr, ": ", nrow(G), "x", ncol(G),
            " (", n_snps_used, " SNPs) — saved")
    n_done <- n_done + 1L
  }
  log_msg("LOCO GRMs: ", n_done, " / ", length(chrs), " chromosomes saved")
  invisible(out_dir)
}

############################################################
# 5) LOAD COHORT-LEVEL OBJECTS
############################################################

sep_line()
log_msg("Step 12ab0 — meQTL5 input preparation (OLD-style inputs)")
log_msg("Output root: ", OUTROOT)

# Pre-load per-cohort objects once (GDS metadata, GRM, PCs, annotations)
# to avoid re-reading for each of the three contexts
gds_meta <- list()
grm_pcs_meta <- list()
grm_meta <- list()
annot_meta <- list()

for (cohort in COHORTS) {
  sep_line()
  log_msg("Loading cohort-level objects: ", cohort)

  # Read SNP and sample metadata from GDS (genotypes are not loaded here)
  gds_path <- GDS_FILES[[cohort]]
  stopifnot(file.exists(gds_path))
  gds      <- snpgdsOpen(gds_path, readonly = TRUE)
  samp_ids <- read.gdsn(index.gdsn(gds, "sample.id"))
  snp_ids  <- read.gdsn(index.gdsn(gds, "snp.id"))
  chr_ids  <- as.character(read.gdsn(index.gdsn(gds, "snp.chromosome")))
  pos_ids  <- as.integer(read.gdsn(index.gdsn(gds, "snp.position")))
  snpgdsClose(gds)

  gds_meta[[cohort]] <- list(
    sample_ids    = data.table(sample_id_raw = as.character(samp_ids),
                               sample_id     = normalize_id(samp_ids)),
    variant_annot = data.table(snp_id = as.character(snp_ids),
                               chr    = chr_ids,
                               pos    = pos_ids)
  )
  log_msg("  GDS: ", nrow(gds_meta[[cohort]]$sample_ids), " samples, ",
          nrow(gds_meta[[cohort]]$variant_annot), " SNPs")

  # GCTA GRM → derive 10 PCs via eigendecomposition
  grm_path <- GRM_FILES[[cohort]]
  stopifnot(file.exists(grm_path))
  G <- readRDS(grm_path)
  stopifnot(is.matrix(G) || inherits(G, "Matrix"))
  stopifnot(!is.null(rownames(G)))
  grm_meta[[cohort]] <- G

  grm_ids_norm <- normalize_id(rownames(G))
  pcs_dt <- grm_pcs(G, grm_ids_norm, N_PCS)
  grm_pcs_meta[[cohort]] <- pcs_dt
  log_msg("  GCTA GRM: ", nrow(G), "x", ncol(G),
          " | PCs derived: PC1..PC", N_PCS)

  # Group column differs by cohort: Family (BREEDING) or Population (NATURAL)
  # used later for within-group HWE imputation of missing genotypes
  annot_path <- ANNOT_FILES[[cohort]]
  stopifnot(file.exists(annot_path))
  annot <- readRDS(annot_path)
  group_col <- if (cohort == "BREEDING") "Family" else "Population"
  stopifnot("IID" %in% names(annot), group_col %in% names(annot))
  annot_dt <- data.table(
    sample_id = normalize_id(as.character(annot$IID)),
    group     = as.character(annot[[group_col]])
  )
  annot_meta[[cohort]] <- annot_dt
  log_msg("  Annotation: ", nrow(annot_dt), " samples, group col = ", group_col)
}

############################################################
# 6) LOCO GRM — BREEDING ONLY
############################################################

sep_line()
log_msg("Computing LOCO GRMs for BREEDING cohort (non-imputed GDS)...")
loco_dir <- file.path(INPUTDIR, "BREEDING", "loco_grm")
compute_loco_grms(GDS_FILES[["BREEDING"]], loco_dir)

############################################################
# 7) PANEL-SPECIFIC INPUTS
############################################################

summary_rows <- list()
overlap_rows <- list()

for (cohort in COHORTS) {
  for (context in CONTEXTS) {
    sep_line()
    log_msg("Panel: ", cohort, " / ", context)
    panel_dir <- file.path(INPUTDIR, cohort, context)
    dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)

    # --- Load unfiltered methylKit RDS → M-values (NO QN) ---
    ctx_lower <- tolower(context)
    coh_lower <- tolower(cohort)
    meth_path <- file.path(TMS_RDS_DIR,
                           paste0(coh_lower, "_", context, "_methylkit.rds"))
    if (!file.exists(meth_path)) {
      log_msg("  Missing methylKit RDS — skipping: ", meth_path)
      summary_rows[[paste(cohort, context)]] <- data.table(
        cohort = cohort, context = context, status = "missing_methylbase")
      next
    }
    log_msg("  methylKit file: ", meth_path)
    conv       <- methylkit_to_mvalues(meth_path, eps = EPSILON_M)
    m_mat      <- conv$m_matrix
    site_annot <- conv$site_annot
    log_msg("  M-value matrix: ", nrow(m_mat), " samples x ", ncol(m_mat), " sites")

    # --- Sample harmonisation ---
    # Find the intersection of samples present in all four data sources
    meth_ids_norm <- normalize_id(rownames(m_mat))
    gds_ids  <- gds_meta[[cohort]]$sample_ids$sample_id
    pcs_ids  <- grm_pcs_meta[[cohort]]$sample_id
    grm_ids  <- normalize_id(rownames(grm_meta[[cohort]]))
    shared_ids <- sort(Reduce(intersect, list(gds_ids, pcs_ids, grm_ids, meth_ids_norm)))

    overlap_rows[[paste(cohort, context)]] <- data.table(
      cohort   = cohort, context = context,
      gds_n    = length(gds_ids),  pcs_n   = length(pcs_ids),
      grm_n    = length(grm_ids),  methyl_n = length(meth_ids_norm),
      shared_n = length(shared_ids)
    )
    log_msg("  Shared samples: ", length(shared_ids),
            " (gds=", length(gds_ids), " pcs=", length(pcs_ids),
            " grm=", length(grm_ids), " meth=", length(meth_ids_norm), ")")

    if (length(shared_ids) < 10L) {
      log_msg("  Too few shared samples — skipping panel")
      summary_rows[[paste(cohort, context)]] <- data.table(
        cohort = cohort, context = context, status = "too_few_samples")
      next
    }

    # --- Restrict matrices ---
    meth_raw  <- rownames(m_mat)
    met_idx   <- match(shared_ids, meth_ids_norm)
    m_shared  <- m_mat[met_idx, , drop = FALSE]
    rownames(m_shared) <- meth_raw[met_idx]

    # Subset PCs to shared samples in the same order
    pcs_shared <- copy(grm_pcs_meta[[cohort]])[
      match(shared_ids, grm_pcs_meta[[cohort]]$sample_id)]

    # Subset GRM to shared samples (symmetric submatrix)
    G      <- grm_meta[[cohort]]
    grm_idx <- match(shared_ids, normalize_id(rownames(G)))
    grm_shared <- G[grm_idx, grm_idx, drop = FALSE]
    rownames(grm_shared) <- colnames(grm_shared) <- shared_ids

    # Group assignments for HWE imputation (used in steps 12ab1/12ab2)
    groups_dt <- annot_meta[[cohort]][match(shared_ids, sample_id)]
    groups_dt[is.na(group), group := "Unknown"]
    groups_dt[, sample_id := shared_ids]

    shared_ids_dt <- data.table(sample_id = shared_ids)

    # --- Write RDS outputs ---
    saveRDS(m_shared,      file.path(panel_dir, "methylation_mvalues_matrix.rds"))
    saveRDS(site_annot,    file.path(panel_dir, "methylation_site_annot.rds"))
    saveRDS(shared_ids_dt, file.path(panel_dir, "shared_sample_ids.rds"))
    saveRDS(pcs_shared,    file.path(panel_dir, "pcs_shared.rds"))
    saveRDS(grm_shared,    file.path(panel_dir, "grm_shared.rds"))
    saveRDS(groups_dt,     file.path(panel_dir, "sample_groups.rds"))
    saveRDS(gds_meta[[cohort]]$sample_ids,    file.path(panel_dir, "snp_gds_sample_ids.rds"))
    saveRDS(gds_meta[[cohort]]$variant_annot, file.path(panel_dir, "snp_variant_annot.rds"))

    # --- Write TSV outputs (human-readable, for inspection) ---
    fwrite(as.data.table(m_shared, keep.rownames = "sample_id_raw"),
           file.path(panel_dir, "methylation_mvalues_matrix.tsv"), sep = "\t")
    fwrite(site_annot,    file.path(panel_dir, "methylation_site_annot.tsv"),  sep = "\t")
    fwrite(shared_ids_dt, file.path(panel_dir, "shared_sample_ids.tsv"),       sep = "\t")
    fwrite(pcs_shared,    file.path(panel_dir, "pcs_shared.tsv"),              sep = "\t")
    fwrite(groups_dt,     file.path(panel_dir, "sample_groups.tsv"),           sep = "\t")
    fwrite(gds_meta[[cohort]]$sample_ids,
           file.path(panel_dir, "snp_gds_sample_ids.tsv"), sep = "\t")
    fwrite(gds_meta[[cohort]]$variant_annot,
           file.path(panel_dir, "snp_variant_annot.tsv"),  sep = "\t")

    summary_rows[[paste(cohort, context)]] <- data.table(
      cohort         = cohort, context = context, status = "ok",
      shared_samples = length(shared_ids),
      methyl_sites   = ncol(m_shared),
      n_pcs          = N_PCS
    )
    log_msg("  Panel inputs written | sites=", ncol(m_shared), " | pcs=", N_PCS)
    rm(m_mat, m_shared, grm_shared, pcs_shared, site_annot, conv, groups_dt); gc()
  }
}

############################################################
# 8) WRITE SUMMARIES
############################################################

sep_line()
log_msg("Writing summaries...")
summary_dt <- rbindlist(summary_rows, fill = TRUE)
overlap_dt <- rbindlist(overlap_rows, fill = TRUE)
fwrite(summary_dt, SUMMARY_FILE, sep = "\t")
fwrite(overlap_dt, OVERLAP_FILE, sep = "\t")
writeLines(capture.output(sessionInfo()), SESSION_FILE)

sep_line()
log_msg("Step 12ab0 finished")
log_msg("Inputs: ",    INPUTDIR)
log_msg("LOCO GRMs: ", loco_dir)
log_msg("Summary: ",   SUMMARY_FILE)
print(summary_dt)
sep_line()

```

---

### `12ab1.tgc.joint.matrixeqtl.mapping.R`

```r
#!/usr/bin/env Rscript
############################################################
# TreeGeneClimate (TGC) — JOINT ECS + TMS
# Step 12ab1: MatrixEQTL5 cis-meQTL mapping
#
# KEY DIFFERENCES vs 15ab2 (MatrixEQTL4):
#   - Inputs from 12ab0 (non-imputed GDS, GCTA GRM PCs, unfiltered methylation)
#   - PC1..PC10 for BOTH cohorts (vs 10/3 in meQTL4)
#   - No quantile-normalised M-values
#   - Cis window: 100 kb (1e5) — same as original ECS pipeline
#   - HWE group imputation of missing genotypes
#   - Lambda (genomic inflation) reported per panel in summary
#
# MODEL
#   M-value ~ SNP + PC1..PC10   (linear, no GRM random effect)
#
# USAGE
# SLURM single panel:
#   TGC_COHORT=BREEDING TGC_CONTEXT=CpG Rscript --vanilla 12ab1.R
# Interactive (all 6 panels):
#   source this script in RStudio (no env vars set)
#
# INPUTS (from 12ab0 — RESULTS/JOINT/MQTL5/INPUTS/<COHORT>/<CONTEXT>/)
#   methylation_mvalues_matrix.rds
#   methylation_site_annot.rds
#   shared_sample_ids.rds
#   pcs_shared.rds
#   sample_groups.rds
#   snp_variant_annot.rds
#   snp_gds_sample_ids.rds
#
# OUTPUTS (RESULTS/JOINT/MATRIXEQTL5/<COHORT>/<CONTEXT>/)
#   cis_meqtl_all_results.rds / .tsv.gz
#   cis_meqtl_significant.tsv        (FDR < 0.05)
#   cis_meqtl_top_per_site.tsv
#   cis_meqtl_panel_summary.tsv      (includes lambda)
#   qq_<cohort>_<context>.tiff
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(MatrixEQTL)
  library(SNPRelate)
  library(gdsfmt)
})

options(stringsAsFactors = FALSE)

############################################################
# 1) SETTINGS
############################################################

CIS_WINDOW   <- 1e5       # 100 kb — same as original ECS pipeline
MIN_SAMPLES  <- 10L
FDR_THRESH   <- 0.05
# Keep all cis pairs (p <= 1) so that genome-wide BH correction is applied
# post-hoc over the complete set rather than a pre-filtered subset
PVAL_CIS_OUT <- 1

############################################################
# 2) PATHS
############################################################

# === USER CONFIGURATION ===
# Set TGC_PROJECT_ROOT as an environment variable, or edit the fallback path below
PROJECT_ROOT <- Sys.getenv("TGC_PROJECT_ROOT",
  unset = "/path/to/your/project")
# ===========================

MQTL5_ROOT  <- file.path(PROJECT_ROOT, "RESULTS", "JOINT", "MQTL5")
INPUTDIR    <- file.path(MQTL5_ROOT, "INPUTS")

MEQTL_ROOT  <- file.path(PROJECT_ROOT, "RESULTS", "JOINT", "MATRIXEQTL5")
OUTROOT     <- MEQTL_ROOT
LOGDIR      <- file.path(MEQTL_ROOT, "LOGS")
SUMDIR      <- file.path(MEQTL_ROOT, "SUMMARIES")

for (d in c(OUTROOT, LOGDIR, SUMDIR))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Non-imputed GDS files (same source as used in 12ab0)
RDATA_DIR <- "/path/to/your/project"
GDS_FILES <- list(
  BREEDING = file.path(RDATA_DIR, "breeding.snp.gds"),
  NATURAL  = file.path(RDATA_DIR, "natural.snp.gds")
)

############################################################
# 3) PHASE DETECTION
############################################################

# When run via SLURM with env vars set, process a single panel.
# When run interactively (no env vars), iterate over all 6 panels.
env_cohort  <- Sys.getenv("TGC_COHORT",  unset = "")
env_context <- Sys.getenv("TGC_CONTEXT", unset = "")
RUN_MAPPING <- nzchar(env_cohort) && nzchar(env_context)
SLURM_MODE  <- RUN_MAPPING

if (RUN_MAPPING) {
  stopifnot(env_cohort  %in% c("BREEDING", "NATURAL"))
  stopifnot(env_context %in% c("CpG", "CHG", "CHH"))
  PANELS <- data.frame(cohort = env_cohort, context = env_context,
                       stringsAsFactors = FALSE)
} else {
  PANELS <- expand.grid(cohort  = c("BREEDING", "NATURAL"),
                        context = c("CpG", "CHG", "CHH"),
                        stringsAsFactors = FALSE)
}

############################################################
# 4) HELPERS
############################################################

CURRENT_LOGFILE <- NULL

log_msg <- function(...) {
  txt <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ",
                paste0(..., collapse = ""))
  cat(txt, "\n")
  if (!is.null(CURRENT_LOGFILE))
    cat(txt, "\n", file = CURRENT_LOGFILE, append = TRUE)
}
sep_line <- function() log_msg(paste(rep("=", 70), collapse = ""))

normalize_id <- function(x) {
  x <- as.character(x); x <- trimws(x)
  x <- gsub("\\.0$", "", x); x <- gsub("^X", "", x); x <- gsub("-", "_", x); x
}

# HWE group imputation: fill missing genotypes with 2*p within each group
# (family for BREEDING, population for NATURAL). Falls back to cohort-wide
# allele frequency if a group has no coverage for a given SNP.
impute_hwe <- function(geno_snp_x_samp, groups) {
  mode(geno_snp_x_samp) <- "numeric"
  geno_snp_x_samp[geno_snp_x_samp > 2 | geno_snp_x_samp < 0] <- NA_real_
  # Cohort-wide allele frequency (fallback for groups with no data)
  p_all <- rowSums(geno_snp_x_samp, na.rm = TRUE) /
           (2 * rowSums(!is.na(geno_snp_x_samp)))
  # Remove SNPs with no valid genotype calls
  keep  <- is.finite(p_all)
  g     <- geno_snp_x_samp[keep, , drop = FALSE]
  p_all <- p_all[keep]
  for (grp in unique(groups)) {
    idx <- which(groups == grp)
    # Within-group allele frequency
    pg  <- rowSums(g[, idx, drop = FALSE], na.rm = TRUE) /
           (2 * rowSums(!is.na(g[, idx, drop = FALSE])))
    pg[!is.finite(pg)] <- p_all[!is.finite(pg)]  # fallback
    miss <- which(is.na(g[, idx, drop = FALSE]), arr.ind = TRUE)
    if (nrow(miss)) g[, idx][miss] <- 2 * pg[miss[, 1]]
  }
  # Final sweep: any remaining NAs filled with cohort-wide expectation
  if (anyNA(g)) {
    fill2 <- matrix(2 * p_all, nrow = nrow(g), ncol = ncol(g))
    g[is.na(g)] <- fill2[is.na(g)]
  }
  list(geno = g, keep = keep)
}

# Genomic inflation factor: ratio of observed to expected median chi-squared statistic
compute_lambda <- function(pvals) {
  pvals <- pvals[is.finite(pvals) & pvals > 0 & pvals <= 1]
  if (length(pvals) < 2) return(NA_real_)
  median(qchisq(1 - pvals, df = 1), na.rm = TRUE) / qchisq(0.5, df = 1)
}

############################################################
# 5) MAPPING LOOP
############################################################

panel_summaries <- list()

for (p in seq_len(nrow(PANELS))) {
  cohort  <- PANELS$cohort[p]
  context <- PANELS$context[p]

  panel_tag    <- paste0(tolower(cohort), "_", tolower(context))
  panel_input  <- file.path(INPUTDIR, cohort, context)
  panel_output <- file.path(OUTROOT, cohort, context)
  dir.create(panel_output, recursive = TRUE, showWarnings = FALSE)

  # Skip completed panels to allow safe re-runs / partial SLURM restarts
  if (file.exists(file.path(panel_output, "cis_meqtl_all_results.rds"))) {
    cat("[", format(Sys.time(), "%H:%M:%S"), "] Panel ", cohort, "/", context,
        " already done — skipping\n", sep = "")
    next
  }

  CURRENT_LOGFILE <- file.path(LOGDIR, paste0("step12ab1_", panel_tag, ".log"))
  SUMMARY_FILE    <- file.path(SUMDIR, paste0("step12ab1_", panel_tag, "_summary.tsv"))
  SESSION_FILE    <- file.path(SUMDIR, paste0("step12ab1_", panel_tag, "_sessionInfo.txt"))
  if (file.exists(CURRENT_LOGFILE)) file.remove(CURRENT_LOGFILE)

  sep_line()
  log_msg("Step 12ab1 — MatrixEQTL5 cis-meQTL mapping")
  log_msg("Panel: ", cohort, " / ", context)
  log_msg("Cis window: +/-", CIS_WINDOW / 1e3, " kb | Model: M ~ SNP + PC1..10")

  # --- Check inputs ---
  req <- c("methylation_mvalues_matrix.rds", "methylation_site_annot.rds",
           "shared_sample_ids.rds", "pcs_shared.rds", "sample_groups.rds",
           "snp_variant_annot.rds")
  missing_f <- req[!file.exists(file.path(panel_input, req))]
  if (length(missing_f)) {
    log_msg("Missing inputs: ", paste(missing_f, collapse = ", "))
    fwrite(data.table(cohort = cohort, context = context,
                      status = "missing_inputs", lambda = NA_real_),
           SUMMARY_FILE, sep = "\t")
    if (SLURM_MODE) quit(save = "no", status = 1L); next
  }
  gds_path <- GDS_FILES[[cohort]]
  if (!file.exists(gds_path)) {
    log_msg("Missing GDS: ", gds_path)
    fwrite(data.table(cohort = cohort, context = context,
                      status = "missing_gds", lambda = NA_real_),
           SUMMARY_FILE, sep = "\t")
    if (SLURM_MODE) quit(save = "no", status = 1L); next
  }

  # --- Load panel inputs ---
  log_msg("Loading panel inputs...")
  m_mat      <- readRDS(file.path(panel_input, "methylation_mvalues_matrix.rds"))
  site_annot <- as.data.table(readRDS(file.path(panel_input, "methylation_site_annot.rds")))
  shared_dt  <- readRDS(file.path(panel_input, "shared_sample_ids.rds"))
  pcs_dt     <- as.data.table(readRDS(file.path(panel_input, "pcs_shared.rds")))
  groups_dt  <- as.data.table(readRDS(file.path(panel_input, "sample_groups.rds")))
  var_annot  <- as.data.table(readRDS(file.path(panel_input, "snp_variant_annot.rds")))

  shared_ids <- shared_dt$sample_id
  n_shared   <- length(shared_ids)
  n_sites    <- ncol(m_mat)
  n_snps     <- nrow(var_annot)
  log_msg("Shared samples: ", n_shared, " | Sites: ", n_sites, " | SNPs: ", n_snps)

  if (n_shared < MIN_SAMPLES) {
    fwrite(data.table(cohort = cohort, context = context,
                      status = "too_few_samples", lambda = NA_real_,
                      n_samples = n_shared),
           SUMMARY_FILE, sep = "\t")
    if (SLURM_MODE) quit(save = "no", status = 1L); next
  }

  # --- Read genotypes from GDS ---
  log_msg("Reading genotypes from GDS...")
  gds           <- snpgdsOpen(gds_path, readonly = TRUE)
  gds_samp_ids  <- read.gdsn(index.gdsn(gds, "sample.id"))
  gds_samp_norm <- normalize_id(gds_samp_ids)
  # Extract only the raw IDs that map to shared normalised IDs
  raw_shared    <- gds_samp_ids[gds_samp_norm %in% shared_ids]
  geno_list     <- snpgdsGetGeno(gds, sample.id = raw_shared,
                                 snpfirstdim = TRUE, with.id = TRUE)
  snpgdsClose(gds)

  # Re-order columns to match shared_ids ordering for alignment with M-value matrix
  samp_reorder <- match(shared_ids, normalize_id(geno_list$sample.id))
  geno_mat     <- geno_list$genotype[, samp_reorder, drop = FALSE]  # SNPs x samples
  rm(geno_list, samp_reorder); gc()
  log_msg("Genotype matrix: ", nrow(geno_mat), " SNPs x ", ncol(geno_mat), " samples")

  # --- HWE group imputation ---
  log_msg("Imputing missing genotypes by group...")
  groups_ord <- groups_dt$group[match(shared_ids, groups_dt$sample_id)]
  groups_ord[is.na(groups_ord)] <- "Unknown"
  imp      <- impute_hwe(geno_mat, groups_ord)
  geno_imp <- imp$geno
  var_filt <- var_annot[imp$keep]  # variant annotation subset to polymorphic SNPs
  rm(geno_mat, imp); gc()

  # Drop zero-variance SNPs (monomorphic after imputation)
  snp_sd   <- apply(geno_imp, 1, sd, na.rm = TRUE)
  keep_var <- is.finite(snp_sd) & snp_sd > 0
  geno_imp <- geno_imp[keep_var, , drop = FALSE]
  var_filt <- var_filt[keep_var]
  log_msg("SNPs after imputation + variance filter: ", nrow(geno_imp))

  # --- Build MatrixEQTL SlicedData objects ---
  log_msg("Building MatrixEQTL SlicedData objects...")
  m_reorder  <- match(shared_ids, normalize_id(rownames(m_mat)))
  m_ord      <- m_mat[m_reorder, , drop = FALSE]

  rownames(geno_imp) <- as.character(var_filt$snp_id)
  colnames(geno_imp) <- shared_ids
  rownames(m_ord)    <- shared_ids

  # MatrixEQTL expects features x samples orientation
  snpsSD <- SlicedData$new(); snpsSD$CreateFromMatrix(geno_imp)
  geneSD <- SlicedData$new(); geneSD$CreateFromMatrix(t(m_ord))   # sites x samples

  pc_cols <- grep("^PC[0-9]+$", names(pcs_dt), value = TRUE)
  pcs_ord <- pcs_dt[match(shared_ids, pcs_dt$sample_id), ..pc_cols]
  cvrtSD  <- SlicedData$new()
  cvrtSD$CreateFromMatrix(t(as.matrix(pcs_ord)))
  log_msg("Covariates: ", paste(pc_cols, collapse = ", "))

  # Position tables required by MatrixEQTL to define the cis window
  snpspos <- data.frame(
    snp = as.character(var_filt$snp_id),
    chr = as.character(var_filt$chr),
    pos = as.integer(var_filt$pos),
    stringsAsFactors = FALSE
  )
  genepos <- data.frame(
    gene  = colnames(m_ord),
    chr   = as.character(site_annot$chr),
    start = as.integer(site_annot$start),
    end   = as.integer(site_annot$end),
    stringsAsFactors = FALSE
  )

  # --- Run MatrixEQTL ---
  # pvOutputThreshold = 0 suppresses trans output entirely (cis-only analysis)
  cis_out_file <- file.path(panel_output,
                             paste0(panel_tag, "_MatrixEQTL_cis.txt"))
  log_msg("Running Matrix_eQTL_main (cis, window=", CIS_WINDOW, " bp)...")
  me <- tryCatch(
    Matrix_eQTL_main(
      snps                  = snpsSD,
      gene                  = geneSD,
      cvrt                  = cvrtSD,
      output_file_name      = "",
      pvOutputThreshold     = 0,
      useModel              = modelLINEAR,
      errorCovariance       = numeric(),
      verbose               = TRUE,
      output_file_name.cis  = cis_out_file,
      pvOutputThreshold.cis = PVAL_CIS_OUT,
      snpspos               = snpspos,
      genepos               = genepos,
      cisDist               = CIS_WINDOW
    ),
    error = function(e) { log_msg("MatrixEQTL error: ", conditionMessage(e)); NULL }
  )

  if (is.null(me)) {
    fwrite(data.table(cohort = cohort, context = context,
                      status = "matrixeqtl_error", lambda = NA_real_),
           SUMMARY_FILE, sep = "\t")
    if (SLURM_MODE) quit(save = "no", status = 1L); next
  }

  cis_res <- me$cis$eqtls
  if (is.null(cis_res) || nrow(cis_res) == 0L) {
    log_msg("No cis results for ", cohort, " / ", context)
    fwrite(data.table(cohort = cohort, context = context,
                      status = "no_cis_results", lambda = NA_real_,
                      n_samples = n_shared, n_sites = n_sites,
                      n_snps = nrow(geno_imp), n_pairs = 0L),
           SUMMARY_FILE, sep = "\t")
    if (SLURM_MODE) quit(save = "no", status = 1L); next
  }

  # --- Post-process ---
  cis_dt <- as.data.table(cis_res)
  setnames(cis_dt, c("snps", "gene"), c("snp", "site"))
  # BH correction applied across all cis pairs for this panel
  cis_dt[, p_FDR := p.adjust(pvalue, method = "BH")]

  # Annotate each pair with genomic positions and physical distance
  snp_pos_map  <- setNames(snpspos$pos,   snpspos$snp)
  site_pos_map <- setNames(genepos$start, genepos$gene)
  cis_dt[, snp_pos  := snp_pos_map[snp]]
  cis_dt[, site_pos := site_pos_map[site]]
  cis_dt[, distance := abs(snp_pos - site_pos)]

  lambda <- compute_lambda(cis_dt$pvalue)
  log_msg("Lambda (genomic inflation): ", round(lambda, 4))

  sig_dt <- cis_dt[p_FDR < FDR_THRESH][order(p_FDR)]
  # Best SNP per methylation site (minimum p-value) for downstream summarisation
  top_dt <- cis_dt[, .SD[which.min(pvalue)], by = site]

  log_msg("Total cis pairs: ", nrow(cis_dt),
          " | FDR<0.05: ", nrow(sig_dt),
          " | Top per site: ", nrow(top_dt))

  # --- Save outputs ---
  saveRDS(cis_dt, file.path(panel_output, "cis_meqtl_all_results.rds"))
  fwrite(cis_dt[, .(snp, site, statistic, pvalue, FDR = p_FDR,
                    snp_pos, site_pos, distance)],
         file.path(panel_output, "cis_meqtl_all_results.tsv.gz"),
         sep = "\t", compress = "gzip")
  fwrite(sig_dt, file.path(panel_output, "cis_meqtl_significant.tsv"),  sep = "\t")
  fwrite(top_dt, file.path(panel_output, "cis_meqtl_top_per_site.tsv"), sep = "\t")

  # --- QQ plot (TIFF) ---
  pv_qq <- cis_dt$pvalue
  pv_qq <- pv_qq[is.finite(pv_qq) & pv_qq > 0]
  if (length(pv_qq) >= 2) {
    qq_file <- file.path(panel_output, paste0("qq_", panel_tag, ".tiff"))
    n_qq <- length(pv_qq)
    qq_title <- sprintf("MatrixEQTL5 QQ — %s %s\nn=%d, lambda=%.3f",
                        cohort, context, n_qq, lambda)
    draw_qq <- function() {
      # Compare observed -log10(p) distribution against uniform expectation
      plot(-log10(ppoints(n_qq)), -log10(sort(pv_qq)),
           pch = 20, cex = 0.4, col = "#2b8cbe",
           xlab = "Expected -log10(p)", ylab = "Observed -log10(p)",
           main = qq_title)
      abline(0, 1, col = "red", lty = 2)
    }
    tiff(qq_file, width = 2400, height = 2400, res = 300, compression = "lzw")
    draw_qq(); dev.off()
    pdf(sub("\\.tiff$", ".pdf", qq_file), width = 8, height = 8)
    draw_qq(); dev.off()
    cairo_ps(sub("\\.tiff$", ".eps", qq_file), width = 8, height = 8)
    draw_qq(); dev.off()
    png(sub("\\.tiff$", ".png", qq_file), width = 800, height = 800)
    draw_qq(); dev.off()
    log_msg("QQ plot: ", qq_file)
  }

  # --- Summary with lambda ---
  sum_dt <- data.table(
    cohort        = cohort,
    context       = context,
    status        = "ok",
    n_samples     = n_shared,
    n_sites       = n_sites,
    n_snps_tested = nrow(geno_imp),
    n_cis_pairs   = nrow(cis_dt),
    n_top_sites   = nrow(top_dt),
    n_sig_fdr05   = nrow(sig_dt),
    lambda        = round(lambda, 4)
  )
  fwrite(sum_dt, SUMMARY_FILE, sep = "\t")
  panel_summaries[[panel_tag]] <- sum_dt
  writeLines(capture.output(sessionInfo()), SESSION_FILE)

  rm(geno_imp, snpsSD, geneSD, cvrtSD, cis_dt, cis_res, me); gc()
  sep_line()
}

# --- Master summary ---
if (length(panel_summaries) > 0) {
  master_dt <- rbindlist(panel_summaries, fill = TRUE)
  fwrite(master_dt,
         file.path(SUMDIR, "step12ab1_all_panels_summary.tsv"), sep = "\t")
  cat("\nMatrixEQTL5 — all panels summary:\n")
  print(master_dt)
}
log_msg("Step 12ab1 finished. Outputs: ", MEQTL_ROOT)

```

---

### `12ab2.tgc.joint.genesis.mapping.R`

```r
#!/usr/bin/env Rscript
############################################################
# TreeGeneClimate (TGC) — JOINT ECS + TMS
# Step 12ab2: GENESIS5 cis-meQTL mapping
#
# KEY DIFFERENCES vs 15ab1 (GENESIS4):
#   - Inputs from 12ab0 (non-imputed GDS, GCTA GRM PCs, unfiltered methylation)
#   - PC1..PC10 for BOTH cohorts
#   - GRM random effect for BOTH cohorts (not just BREEDING)
#     BREEDING : LOCO GRM (from non-imputed GDS, VanRaden formula)
#     NATURAL  : full GCTA GRM (grm_shared.rds)
#   - Cis window: 100 kb (1e5)
#   - HWE group imputation of missing genotypes
#   - Lambda reported per panel
#
# MODEL
#   BREEDING: M-value ~ PC1..PC10 + random(LOCO-GRM)   [AIREML LMM]
#   NATURAL:  M-value ~ PC1..PC10 + random(full GRM)    [AIREML LMM]
#
# USAGE
# SLURM single panel:
#   TGC_COHORT=BREEDING TGC_CONTEXT=CpG Rscript --vanilla 12ab2.R
# Interactive (all 6 panels):
#   source this script in RStudio (no env vars set)
#
# INPUTS (from 12ab0 — RESULTS/JOINT/MQTL5/INPUTS/<COHORT>/<CONTEXT>/)
#   methylation_mvalues_matrix.rds
#   methylation_site_annot.rds
#   shared_sample_ids.rds
#   pcs_shared.rds
#   grm_shared.rds
#   sample_groups.rds
#   snp_variant_annot.rds
#   snp_gds_sample_ids.rds
# BREEDING LOCO GRMs (INPUTS/BREEDING/loco_grm/loco_grm_excl_<chr>.rds)
#
# OUTPUTS (RESULTS/JOINT/GENESIS5/<COHORT>/<CONTEXT>/)
#   cis_meqtl_all_results.rds / .tsv.gz
#   cis_meqtl_significant.tsv
#   cis_meqtl_top_per_site.tsv
#   cis_meqtl_panel_summary.tsv      (includes lambda)
#   qq_<cohort>_<context>.tiff
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(SNPRelate)
  library(gdsfmt)
  library(GENESIS)      # fitNullModel + assocTestSingle
  library(GWASTools)    # GenotypeData infrastructure
  library(Biobase)
  library(parallel)
})

options(stringsAsFactors = FALSE)

############################################################
# 1) SETTINGS
############################################################

CIS_WINDOW  <- 1e5        # 100 kb
MIN_SAMPLES <- 10L
FDR_THRESH  <- 0.05
# Loose AIREML convergence tolerance: variance components are re-estimated
# per site, so exact convergence is less critical than throughput
AIREML_FAST <- 1e-2
N_CORES     <- min(4L, parallel::detectCores(logical = FALSE))

############################################################
# 2) PATHS
############################################################

# === USER CONFIGURATION ===
# Set TGC_PROJECT_ROOT as an environment variable, or edit the fallback path below
PROJECT_ROOT <- Sys.getenv("TGC_PROJECT_ROOT",
  unset = "/path/to/your/project")
# ===========================

MQTL5_ROOT   <- file.path(PROJECT_ROOT, "RESULTS", "JOINT", "MQTL5")
INPUTDIR     <- file.path(MQTL5_ROOT, "INPUTS")

GENESIS5_ROOT <- file.path(PROJECT_ROOT, "RESULTS", "JOINT", "GENESIS5")
OUTROOT       <- file.path(GENESIS5_ROOT, "MEQTL")
LOGDIR        <- file.path(GENESIS5_ROOT, "LOGS")
SUMDIR        <- file.path(GENESIS5_ROOT, "SUMMARIES")

for (d in c(OUTROOT, LOGDIR, SUMDIR))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Non-imputed GDS files (genotypes loaded fresh per panel)
RDATA_DIR <- "/path/to/your/project"
GDS_FILES <- list(
  BREEDING = file.path(RDATA_DIR, "breeding.snp.gds"),
  NATURAL  = file.path(RDATA_DIR, "natural.snp.gds")
)

############################################################
# 3) PHASE DETECTION
############################################################

# Single-panel SLURM mode when env vars are set; all-panels interactive mode otherwise
env_cohort  <- Sys.getenv("TGC_COHORT",  unset = "")
env_context <- Sys.getenv("TGC_CONTEXT", unset = "")
RUN_MAPPING <- nzchar(env_cohort) && nzchar(env_context)
SLURM_MODE  <- RUN_MAPPING

if (RUN_MAPPING) {
  stopifnot(env_cohort  %in% c("BREEDING", "NATURAL"))
  stopifnot(env_context %in% c("CpG", "CHG", "CHH"))
  PANELS <- data.frame(cohort = env_cohort, context = env_context,
                       stringsAsFactors = FALSE)
} else {
  PANELS <- expand.grid(cohort  = c("BREEDING", "NATURAL"),
                        context = c("CpG", "CHG", "CHH"),
                        stringsAsFactors = FALSE)
}

############################################################
# 4) HELPERS
############################################################

CURRENT_LOGFILE   <- NULL
CURRENT_DEBUGFILE <- NULL

log_msg <- function(...) {
  txt <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ",
                paste0(..., collapse = ""))
  cat(txt, "\n")
  if (!is.null(CURRENT_LOGFILE))
    cat(txt, "\n", file = CURRENT_LOGFILE, append = TRUE)
}
# Verbose diagnostics written only to a separate debug log to keep the main log clean
debug_msg <- function(...) {
  txt <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ",
                paste0(..., collapse = ""))
  if (!is.null(CURRENT_DEBUGFILE))
    cat(txt, "\n", file = CURRENT_DEBUGFILE, append = TRUE)
}
sep_line <- function() log_msg(paste(rep("=", 70), collapse = ""))

normalize_id <- function(x) {
  x <- as.character(x); x <- trimws(x)
  x <- gsub("\\.0$", "", x); x <- gsub("^X", "", x); x <- gsub("-", "_", x); x
}

# HWE group imputation (SNPs x samples):
# fill missing calls with 2*p estimated within the same family/population group
impute_hwe <- function(geno_snp_x_samp, groups) {
  mode(geno_snp_x_samp) <- "numeric"
  geno_snp_x_samp[geno_snp_x_samp > 2 | geno_snp_x_samp < 0] <- NA_real_
  # Cohort-wide allele frequency as fallback
  p_all <- rowSums(geno_snp_x_samp, na.rm = TRUE) /
           (2 * rowSums(!is.na(geno_snp_x_samp)))
  keep  <- is.finite(p_all)
  g     <- geno_snp_x_samp[keep, , drop = FALSE]
  p_all <- p_all[keep]
  for (grp in unique(groups)) {
    idx <- which(groups == grp)
    pg  <- rowSums(g[, idx, drop = FALSE], na.rm = TRUE) /
           (2 * rowSums(!is.na(g[, idx, drop = FALSE])))
    pg[!is.finite(pg)] <- p_all[!is.finite(pg)]
    miss <- which(is.na(g[, idx, drop = FALSE]), arr.ind = TRUE)
    if (nrow(miss)) g[, idx][miss] <- 2 * pg[miss[, 1]]
  }
  # Global fallback for any remaining NAs
  if (anyNA(g)) {
    fill2 <- matrix(2 * p_all, nrow = nrow(g), ncol = ncol(g))
    g[is.na(g)] <- fill2[is.na(g)]
  }
  list(geno = g, keep = keep)
}

# Load the LOCO GRM for a given chromosome and subset/reindex to shared samples.
# The GRM row/col names are set to integer scan IDs required by GWASTools.
load_loco_grm <- function(loco_dir, chr_name, shared_ids, scan_int) {
  path <- file.path(loco_dir, paste0("loco_grm_excl_", chr_name, ".rds"))
  if (!file.exists(path)) return(NULL)
  G        <- readRDS(path)
  grm_norm <- normalize_id(rownames(G))
  idx      <- match(shared_ids, grm_norm)
  bad      <- is.na(idx)
  if (any(bad)) {
    debug_msg("LOCO GRM missing ", sum(bad), " shared IDs for chr ", chr_name)
    # Substitute a valid row index for unmatched samples (affects only diagonal)
    idx[bad] <- which(!is.na(grm_norm))[1L]
  }
  G_sub <- G[idx, idx, drop = FALSE]
  rownames(G_sub) <- colnames(G_sub) <- as.character(scan_int)
  diag(G_sub) <- diag(G_sub) + 1e-4  # small ridge to ensure positive definiteness
  G_sub
}

# Genomic inflation factor
compute_lambda <- function(pvals) {
  pvals <- pvals[is.finite(pvals) & pvals > 0 & pvals <= 1]
  if (length(pvals) < 2) return(NA_real_)
  median(qchisq(1 - pvals, df = 1), na.rm = TRUE) / qchisq(0.5, df = 1)
}

############################################################
# 5) MAPPING LOOP
############################################################

panel_summaries <- list()

for (p in seq_len(nrow(PANELS))) {
  cohort  <- PANELS$cohort[p]
  context <- PANELS$context[p]

  panel_tag    <- paste0(tolower(cohort), "_", tolower(context))
  panel_input  <- file.path(INPUTDIR, cohort, context)
  panel_output <- file.path(OUTROOT, cohort, context)
  dir.create(panel_output, recursive = TRUE, showWarnings = FALSE)

  # Skip completed panels — allows safe SLURM restarts
  done_rds <- file.path(panel_output, "cis_meqtl_all_results.rds")
  if (file.exists(done_rds)) {
    cat("[", format(Sys.time(), "%H:%M:%S"), "] Panel ", cohort, "/", context,
        " already done — skipping\n", sep = "")
    next
  }

  CURRENT_LOGFILE   <- file.path(LOGDIR, paste0("step12ab2_", panel_tag, ".log"))
  CURRENT_DEBUGFILE <- file.path(LOGDIR, paste0("step12ab2_", panel_tag, "_debug.log"))
  SUMMARY_FILE      <- file.path(SUMDIR, paste0("step12ab2_", panel_tag, "_summary.tsv"))
  SESSION_FILE      <- file.path(SUMDIR, paste0("step12ab2_", panel_tag, "_sessionInfo.txt"))
  if (file.exists(CURRENT_LOGFILE)) file.remove(CURRENT_LOGFILE)

  loco_dir <- file.path(INPUTDIR, "BREEDING", "loco_grm")
  # Both cohorts use a GRM random effect; BREEDING uses per-chromosome LOCO GRMs
  USE_GRM   <- TRUE
  USE_LOCO  <- (cohort == "BREEDING")

  sep_line()
  log_msg("Step 12ab2 — GENESIS5 cis-meQTL mapping")
  log_msg("Panel: ", cohort, " / ", context)
  log_msg("Cis window: +/-", CIS_WINDOW / 1e3, " kb")
  log_msg("GRM: ", if (USE_LOCO) "LOCO (BREEDING)" else "full GCTA GRM (NATURAL)")

  # --- Check inputs ---
  req <- c("methylation_mvalues_matrix.rds", "methylation_site_annot.rds",
           "shared_sample_ids.rds", "pcs_shared.rds", "grm_shared.rds",
           "sample_groups.rds", "snp_variant_annot.rds")
  missing_f <- req[!file.exists(file.path(panel_input, req))]
  if (length(missing_f)) {
    log_msg("Missing inputs: ", paste(missing_f, collapse = ", "))
    fwrite(data.table(cohort = cohort, context = context,
                      status = "missing_inputs", lambda = NA_real_),
           SUMMARY_FILE, sep = "\t")
    if (SLURM_MODE) quit(save = "no", status = 1L); next
  }
  gds_path <- GDS_FILES[[cohort]]
  if (!file.exists(gds_path)) {
    log_msg("Missing GDS: ", gds_path)
    fwrite(data.table(cohort = cohort, context = context,
                      status = "missing_gds", lambda = NA_real_),
           SUMMARY_FILE, sep = "\t")
    if (SLURM_MODE) quit(save = "no", status = 1L); next
  }

  # --- Load panel inputs ---
  log_msg("Loading panel inputs...")
  m_mat      <- readRDS(file.path(panel_input, "methylation_mvalues_matrix.rds"))
  site_annot <- as.data.table(readRDS(file.path(panel_input, "methylation_site_annot.rds")))
  shared_dt  <- readRDS(file.path(panel_input, "shared_sample_ids.rds"))
  pcs_dt     <- as.data.table(readRDS(file.path(panel_input, "pcs_shared.rds")))
  groups_dt  <- as.data.table(readRDS(file.path(panel_input, "sample_groups.rds")))
  var_annot  <- as.data.table(readRDS(file.path(panel_input, "snp_variant_annot.rds")))

  shared_ids <- shared_dt$sample_id
  n_shared   <- length(shared_ids)
  n_sites    <- ncol(m_mat)
  n_snps     <- nrow(var_annot)
  log_msg("Shared samples: ", n_shared, " | Sites: ", n_sites, " | SNPs: ", n_snps)

  if (n_shared < MIN_SAMPLES) {
    fwrite(data.table(cohort = cohort, context = context,
                      status = "too_few_samples", lambda = NA_real_,
                      n_samples = n_shared),
           SUMMARY_FILE, sep = "\t")
    if (SLURM_MODE) quit(save = "no", status = 1L); next
  }

  # --- Read genotypes from GDS ---
  log_msg("Reading genotypes from GDS...")
  gds           <- snpgdsOpen(gds_path, readonly = TRUE)
  gds_samp_ids  <- read.gdsn(index.gdsn(gds, "sample.id"))
  gds_samp_norm <- normalize_id(gds_samp_ids)
  raw_shared    <- gds_samp_ids[gds_samp_norm %in% shared_ids]
  # snpfirstdim=FALSE returns samples x SNPs (GWASTools convention)
  geno_list     <- snpgdsGetGeno(gds, sample.id = raw_shared,
                                 snpfirstdim = FALSE, with.id = TRUE)
  snpgdsClose(gds)

  samp_reorder <- match(shared_ids, normalize_id(geno_list$sample.id))
  geno_mat     <- geno_list$genotype[samp_reorder, , drop = FALSE]  # samples x SNPs
  rm(geno_list, samp_reorder); gc()

  # Transpose to SNPs x samples for imputation (impute_hwe expects this orientation)
  geno_snp_x_samp <- t(geno_mat); rm(geno_mat); gc()
  log_msg("Genotype matrix: ", nrow(geno_snp_x_samp), " SNPs x ",
          ncol(geno_snp_x_samp), " samples")

  # --- HWE group imputation ---
  log_msg("Imputing missing genotypes by group...")
  groups_ord <- groups_dt$group[match(shared_ids, groups_dt$sample_id)]
  groups_ord[is.na(groups_ord)] <- "Unknown"
  imp             <- impute_hwe(geno_snp_x_samp, groups_ord)
  geno_imp        <- imp$geno     # SNPs x samples (filtered + imputed)
  var_filt        <- var_annot[imp$keep]
  rm(geno_snp_x_samp, imp); gc()

  # Drop zero-variance SNPs (monomorphic after imputation)
  snp_sd   <- apply(geno_imp, 1, sd, na.rm = TRUE)
  keep_var <- is.finite(snp_sd) & snp_sd > 0
  geno_imp <- geno_imp[keep_var, , drop = FALSE]
  var_filt <- var_filt[keep_var]
  log_msg("SNPs after imputation + variance filter: ", nrow(geno_imp))

  # Transpose back to samples x SNPs (integer storage required by MatrixGenotypeReader)
  geno_imp_t <- t(geno_imp); storage.mode(geno_imp_t) <- "integer"
  rm(geno_imp); gc()

  # --- Build GWASTools objects ---
  # GWASTools requires integer scan IDs and integer-coded chromosomes
  log_msg("Building GWASTools GenotypeData object...")
  scan_int   <- seq_len(n_shared)
  snp_int    <- seq_len(nrow(var_filt))
  chr_levels <- sort(unique(var_filt$chr))
  chr_map    <- setNames(seq_along(chr_levels), chr_levels)
  chr_int    <- as.integer(chr_map[var_filt$chr])

  reader    <- MatrixGenotypeReader(
    genotype   = t(geno_imp_t),   # SNPs x samples
    snpID      = snp_int,
    chromosome = chr_int,
    position   = as.integer(var_filt$pos),
    scanID     = scan_int
  )
  scanAnnot <- ScanAnnotationDataFrame(data.frame(scanID = scan_int))
  genoData  <- GenotypeData(reader, scanAnnot = scanAnnot)
  rm(geno_imp_t); gc()

  # --- Full GRM (NATURAL: primary; BREEDING: fallback when LOCO unavailable) ---
  log_msg("Loading full GRM...")
  grm_full <- readRDS(file.path(panel_input, "grm_shared.rds"))
  grm_idx  <- match(shared_ids, normalize_id(rownames(grm_full)))
  grm_full_ord <- grm_full[grm_idx, grm_idx, drop = FALSE]
  # Use integer scan IDs as row/col names to match GWASTools scanID
  rownames(grm_full_ord) <- colnames(grm_full_ord) <- as.character(scan_int)
  # Ensure positive definiteness with small ridge
  diag(grm_full_ord) <- diag(grm_full_ord) + 1e-4
  rm(grm_full); gc()

  # --- PC covariates ---
  pc_cols  <- grep("^PC[0-9]+$", names(pcs_dt), value = TRUE)
  pcs_ord  <- as.data.frame(pcs_dt[match(shared_ids, pcs_dt$sample_id)])
  # pheno_base is the covariate data frame passed to fitNullModel for every site
  pheno_base <- data.frame(sample.id = scan_int)
  for (pc in pc_cols) pheno_base[[pc]] <- pcs_ord[[pc]]
  log_msg("Covariates: ", paste(pc_cols, collapse = ", "))

  # --- Reorder M-values to match shared_ids ---
  m_reorder <- match(shared_ids, normalize_id(rownames(m_mat)))
  m_mat     <- m_mat[m_reorder, , drop = FALSE]

  # --- Index by chromosome for the cis loop ---
  site_annot[, site_idx := .I]  # row index for fast column access in m_mat
  var_filt[,   snp_int  := seq_len(.N)]
  site_by_chr <- split(site_annot, site_annot$chr)
  snp_by_chr  <- split(var_filt,   var_filt$chr)
  chrs_shared <- intersect(names(site_by_chr), names(snp_by_chr))
  log_msg("Chromosomes with sites + SNPs: ", length(chrs_shared))

  # --- cis-meQTL loop (parallel over sites within each chromosome) ---
  # One null model is fitted per site (AIREML estimates variance components
  # from M ~ PCs + GRM). SNP association is then scored against this null.
  log_msg("Starting GENESIS cis-meQTL mapping (", N_CORES, " cores)...")
  all_results       <- list()
  n_tested          <- 0L
  n_skip_no_cis     <- 0L
  n_skip_null_fail  <- 0L
  n_skip_assoc_fail <- 0L
  n_tests_total     <- 0L
  t_start           <- proc.time()

  for (chr_name in chrs_shared) {

    sites_chr   <- site_by_chr[[chr_name]]
    snps_chr    <- snp_by_chr[[chr_name]]
    snp_pos_chr <- snps_chr$pos

    # For BREEDING: swap in the LOCO GRM (chromosome excluded) to avoid
    # proximal contamination of the variance-component estimates.
    # For NATURAL: always use the full GRM.
    if (USE_LOCO) {
      loco_candidate <- load_loco_grm(loco_dir, chr_name, shared_ids, scan_int)
      grm_chr <- if (!is.null(loco_candidate)) loco_candidate else grm_full_ord
    } else {
      grm_chr <- grm_full_ord
    }

    # Parallelise over sites within the chromosome
    chr_res <- mclapply(seq_len(nrow(sites_chr)), function(j) {
      site <- sites_chr[j, ]

      # Identify cis-SNPs within the 100 kb window
      cis_mask <- abs(snp_pos_chr - site[["pos"]]) <= CIS_WINDOW
      if (!any(cis_mask)) return(list(status = "no_cis"))

      cis_snps  <- snps_chr[cis_mask]
      cis_snpID <- cis_snps[["snp_int"]]

      # Phenotype: M-value vector for this site
      y_vec <- m_mat[, site[["site_idx"]]]
      if (sum(!is.na(y_vec)) < MIN_SAMPLES) return(list(status = "no_cis"))

      pheno_null   <- pheno_base
      pheno_null$y <- y_vec

      # Fit null LMM: M ~ PC1..PC10 + random(GRM) via AIREML
      null_mod <- tryCatch(
        fitNullModel(pheno_null, outcome = "y", covars = pc_cols,
                     cov.mat = grm_chr, family = "gaussian",
                     AIREML.tol = AIREML_FAST, verbose = FALSE),
        error = function(e) NULL
      )
      if (is.null(null_mod)) return(list(status = "null_fail"))

      # Score test of each cis-SNP against the null model residuals
      iterator <- GenotypeBlockIterator(genoData, snpInclude = cis_snpID)
      assoc <- tryCatch(
        assocTestSingle(iterator, null.model = null_mod, verbose = FALSE),
        error = function(e) NULL
      )
      if (is.null(assoc) || nrow(assoc) == 0L) return(list(status = "assoc_fail"))

      assoc_dt <- as.data.table(assoc)
      .methyl_loc <- site[["methyl_loc"]]
      .site_pos   <- as.integer(site[["pos"]])
      .snp_pos    <- cis_snps[["pos"]][match(assoc_dt$variant.id, cis_snps[["snp_int"]])]
      assoc_dt[, `:=`(
        site     = .methyl_loc,
        site_chr = chr_name,
        site_pos = .site_pos,
        snp_pos  = .snp_pos,
        distance = abs(.site_pos - .snp_pos)
      )]
      list(status = "ok", dt = assoc_dt)
    }, mc.cores = N_CORES, mc.preschedule = TRUE)

    # Collect results and update counters
    for (res in chr_res) {
      if (is.null(res) || !is.list(res)) { n_skip_assoc_fail <- n_skip_assoc_fail + 1L; next }
      if (res$status == "ok") {
        all_results[[length(all_results) + 1L]] <- res$dt
        n_tested      <- n_tested + 1L
        n_tests_total <- n_tests_total + nrow(res$dt)
      } else if (res$status == "no_cis") {
        n_skip_no_cis <- n_skip_no_cis + 1L
      } else if (res$status == "null_fail") {
        n_skip_null_fail <- n_skip_null_fail + 1L
      } else {
        n_skip_assoc_fail <- n_skip_assoc_fail + 1L
      }
    }

    elapsed <- (proc.time() - t_start)[3]
    log_msg("Chr ", chr_name, " | tested=", n_tested,
            " | null-fail=", n_skip_null_fail,
            " | ", round(elapsed / 60, 1), " min elapsed")
  }

  # --- Combine results ---
  log_msg("Combining results...")
  log_msg("  Sites tested: ", n_tested, " | Pairs: ", n_tests_total)
  log_msg("  Skipped no-cis: ", n_skip_no_cis,
          " | null-fail: ", n_skip_null_fail,
          " | assoc-fail: ", n_skip_assoc_fail)

  if (length(all_results) == 0L) {
    log_msg("No results produced for ", cohort, " / ", context)
    fwrite(data.table(cohort = cohort, context = context,
                      status = "no_results", lambda = NA_real_),
           SUMMARY_FILE, sep = "\t")
    if (SLURM_MODE) quit(save = "no", status = 1L); next
  }

  cis_dt <- rbindlist(all_results, fill = TRUE)
  rm(all_results); gc()

  # Standardise column names across GENESIS versions (column names vary by version)
  p_col <- intersect(c("Score.pval", "Score.Stat.p", "pval"), names(cis_dt))[1]
  t_col <- intersect(c("Score.Stat", "Score"), names(cis_dt))[1]
  snp_col <- intersect(c("snpID", "variant.id"), names(cis_dt))[1]
  setnames(cis_dt, c(snp_col, p_col, t_col), c("snp_int", "pvalue", "statistic"),
           skip_absent = TRUE)

  # Map integer SNP index back to the original SNP identifier string
  snp_name_map <- setNames(var_filt$snp_id, var_filt$snp_int)
  cis_dt[, snp := snp_name_map[as.character(snp_int)]]

  # BH correction across all cis pairs for this panel
  cis_dt[, p_FDR := p.adjust(pvalue, method = "BH")]

  lambda <- compute_lambda(cis_dt$pvalue)
  log_msg("Lambda (genomic inflation): ", round(lambda, 4))

  sig_dt <- cis_dt[p_FDR < FDR_THRESH][order(p_FDR)]
  # Best SNP per methylation site
  top_dt <- cis_dt[, .SD[which.min(pvalue)], by = site]

  log_msg("Total cis pairs: ", nrow(cis_dt),
          " | FDR<0.05: ", nrow(sig_dt),
          " | Top per site: ", nrow(top_dt))

  # --- Save outputs ---
  out_cols <- c("snp", "site", "statistic", "pvalue", "p_FDR",
                "snp_pos", "site_pos", "distance")
  out_cols <- intersect(out_cols, names(cis_dt))

  saveRDS(cis_dt, file.path(panel_output, "cis_meqtl_all_results.rds"))
  fwrite(cis_dt[, ..out_cols],
         file.path(panel_output, "cis_meqtl_all_results.tsv.gz"),
         sep = "\t", compress = "gzip")
  fwrite(sig_dt[, ..out_cols],
         file.path(panel_output, "cis_meqtl_significant.tsv"), sep = "\t")
  fwrite(top_dt[, ..out_cols],
         file.path(panel_output, "cis_meqtl_top_per_site.tsv"), sep = "\t")

  # --- QQ plot (TIFF) ---
  pv_qq <- cis_dt$pvalue
  pv_qq <- pv_qq[is.finite(pv_qq) & pv_qq > 0]
  if (length(pv_qq) >= 2) {
    qq_file <- file.path(panel_output, paste0("qq_", panel_tag, ".tiff"))
    n_qq <- length(pv_qq)
    qq_title <- sprintf("GENESIS5 QQ — %s %s\nn=%d, lambda=%.3f",
                        cohort, context, n_qq, lambda)
    draw_qq <- function() {
      plot(-log10(ppoints(n_qq)), -log10(sort(pv_qq)),
           pch = 20, cex = 0.4, col = "#e34a33",
           xlab = "Expected -log10(p)", ylab = "Observed -log10(p)",
           main = qq_title)
      abline(0, 1, col = "red", lty = 2)
    }
    tiff(qq_file, width = 2400, height = 2400, res = 300, compression = "lzw")
    draw_qq(); dev.off()
    pdf(sub("\\.tiff$", ".pdf", qq_file), width = 8, height = 8)
    draw_qq(); dev.off()
    cairo_ps(sub("\\.tiff$", ".eps", qq_file), width = 8, height = 8)
    draw_qq(); dev.off()
    png(sub("\\.tiff$", ".png", qq_file), width = 800, height = 800)
    draw_qq(); dev.off()
    log_msg("QQ plot: ", qq_file)
  }

  # --- Summary ---
  sum_dt <- data.table(
    cohort        = cohort,
    context       = context,
    status        = "ok",
    grm_type      = if (USE_LOCO) "LOCO" else "full_GCTA",
    n_samples     = n_shared,
    n_sites       = n_sites,
    n_snps_tested = nrow(var_filt),
    n_cis_pairs   = nrow(cis_dt),
    n_top_sites   = nrow(top_dt),
    n_sig_fdr05   = nrow(sig_dt),
    lambda        = round(lambda, 4)
  )
  fwrite(sum_dt, SUMMARY_FILE, sep = "\t")
  panel_summaries[[panel_tag]] <- sum_dt
  writeLines(capture.output(sessionInfo()), SESSION_FILE)

  rm(cis_dt, genoData, grm_full_ord, m_mat); gc()
  sep_line()
}

# --- Master summary ---
if (length(panel_summaries) > 0) {
  master_dt <- rbindlist(panel_summaries, fill = TRUE)
  fwrite(master_dt,
         file.path(SUMDIR, "step12ab2_all_panels_summary.tsv"), sep = "\t")
  cat("\nGENESIS5 — all panels summary:\n")
  print(master_dt)
}
log_msg("Step 12ab2 finished. Outputs: ", GENESIS5_ROOT)

```

---

### `13ab.tgc.joint.meqtl.combined.results.R`

```r
#!/usr/bin/env Rscript
############################################################
# TreeGeneClimate (TGC) — JOINT ECS + TMS
# Step 13ab: Combined results — comprehensive summary, QQ plots, sig-site tables,
#             and robust SNP-methylation pair identification
#
# PURPOSE
# Reads cis-meQTL results from GENESIS5 and MatrixEQTL5, produces:
#   - Comprehensive summary (lambda, n_sites, n_sig, cis window — all panels)
#   - QQ plots (one per panel x tool, 12 total, TIFF)
#   - Long-format site FDR table (site, site_chr, context, min_FDR)
#   - Significant site lists (p_FDR < 5e-8 and p_FDR < 1e-10)
#   - Robust pairs (FDR < 1e-10 in BOTH tools) with true genomic coordinates
#   - Supplementary Table S5 xlsx
#
# SIGNIFICANCE THRESHOLDS (BH-adjusted p-values):
#   FDR_LOOSE  = 5e-8
#   FDR_STRICT = 1e-10
#
# USAGE
#   Rscript --vanilla 13ab.R
#
# INPUTS
#   RESULTS/JOINT/GENESIS5/MEQTL/<COHORT>/<CONTEXT>/cis_meqtl_all_results.rds
#   RESULTS/JOINT/MATRIXEQTL5/<COHORT>/<CONTEXT>/cis_meqtl_all_results.rds
#   RESULTS/JOINT/MQTL5/INPUTS/<COHORT>/<CONTEXT>/snp_variant_annot.rds
#
# OUTPUTS (RESULTS/JOINT/COMBINED5/)
#   comprehensive_summary.tsv
#   qq/qqplot_<tool>_<cohort>_<context>.tiff
#   tables/all_sites_<tool>_<cohort>.tsv   — long format: site, site_chr, context, min_FDR
#   sig_sites/sig_p5e8_<tool>_<cohort>.tsv
#   sig_sites/sig_p1e10_<tool>_<cohort>.tsv
#   overlap/tables/robust_markers_<cohort>.tsv
#   overlap/tables/robust_context_summary.tsv
#   tables/supplementary_table_s5.xlsx
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(openxlsx2)
})

options(stringsAsFactors = FALSE)

############################################################
# 1) SETTINGS
############################################################

COHORTS       <- c("BREEDING", "NATURAL")
CONTEXTS      <- c("CpG", "CHG", "CHH")
TOOLS         <- c("GENESIS5", "MATRIXEQTL5")
FDR_LOOSE     <- 5e-8   # genome-wide FDR threshold (loose; used for sig-site tables)
FDR_STRICT    <- 1e-10  # strict threshold used in all downstream analyses
CIS_WINDOW_KB <- 100L   # cis window radius applied during meQTL mapping

############################################################
# 2) PATHS
############################################################

# === USER CONFIGURATION ===
# Set TGC_PROJECT_ROOT as an environment variable, or edit the fallback path below
PROJECT_ROOT <- Sys.getenv("TGC_PROJECT_ROOT",
  unset = "/path/to/your/project")
# ===========================

# GENESIS5 results are stored under MEQTL/; MatrixEQTL5 directly under MATRIXEQTL5/
RESULT_ROOTS <- list(
  GENESIS5    = file.path(PROJECT_ROOT, "RESULTS", "JOINT", "GENESIS5",    "MEQTL"),
  MATRIXEQTL5 = file.path(PROJECT_ROOT, "RESULTS", "JOINT", "MATRIXEQTL5")
)

# Variant annotation RDS files produced by the MQTL5 input-prep step (12ab0.R)
ANNOT_ROOT <- file.path(PROJECT_ROOT, "RESULTS", "JOINT", "MQTL5", "INPUTS")

OUT_ROOT <- file.path(PROJECT_ROOT, "RESULTS", "JOINT", "COMBINED5")
QQ_DIR   <- file.path(OUT_ROOT, "qq")
TAB_DIR  <- file.path(OUT_ROOT, "tables")
SIG_DIR  <- file.path(OUT_ROOT, "sig_sites")
ROBUST_DIR <- file.path(OUT_ROOT, "overlap", "tables")
S5_PATH    <- file.path(TAB_DIR, "supplementary_table_s5.xlsx")
for (d in c(OUT_ROOT, QQ_DIR, TAB_DIR, SIG_DIR, ROBUST_DIR))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

############################################################
# 3) HELPERS
############################################################

msg <- function(...) cat(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n"))

# Genomic inflation factor lambda — ratio of observed to expected chi-squared median.
# Values > 1 indicate p-value inflation relative to null; used as QC metric.
compute_lambda <- function(pvals) {
  pvals <- pvals[is.finite(pvals) & pvals > 0 & pvals <= 1]
  if (length(pvals) < 2) return(NA_real_)
  median(qchisq(1 - pvals, df = 1), na.rm = TRUE) / qchisq(0.5, df = 1)
}

# Standardise columns from either tool's RDS output.
# Returns a data.table with columns:
#   snp, snp_chr, snp_pos, site, site_chr, site_pos,
#   context, statistic, beta, pvalue, p_FDR, distance
clean_result <- function(dt, tool, context, va) {
  dt <- as.data.table(dt)

  if (tool == "GENESIS5") {
    # snp_int is the raw GENESIS variant.id; snp (added by 12ab2.R) is the
    # canonical integer matching snp_variant_annot$snp_id — keep snp, drop snp_int
    drop_cols <- intersect(names(dt),
                           c("snp_int", "chr", "pos", "n.obs", "freq", "MAC",
                             "Score", "Score.SE", "Est.SE", "PVE"))
    if (length(drop_cols)) dt[, (drop_cols) := NULL]
    # Est -> beta
    if ("Est" %in% names(dt) && !"beta" %in% names(dt))
      setnames(dt, "Est", "beta")
    # Ensure p_FDR is present; compute from raw p-values if missing
    if (!"p_FDR" %in% names(dt) && "FDR" %in% names(dt))
      setnames(dt, "FDR", "p_FDR")
    if (!"p_FDR" %in% names(dt))
      dt[, p_FDR := p.adjust(pvalue, method = "BH")]
  }

  if (tool == "MATRIXEQTL5") {
    # FDR (MatrixEQTL native) and p_FDR (added by 12ab1.R) are identical — keep p_FDR
    if ("FDR" %in% names(dt) && "p_FDR" %in% names(dt))
      dt[, FDR := NULL]
    else if ("FDR" %in% names(dt) && !"p_FDR" %in% names(dt))
      setnames(dt, "FDR", "p_FDR")
    if (!"p_FDR" %in% names(dt))
      dt[, p_FDR := p.adjust(pvalue, method = "BH")]
    if (!"beta" %in% names(dt))
      dt[, beta := NA_real_]
  }

  # site_chr: parse from site name if not already present
  # e.g. "PA_chr03:790910176-790910176" -> "PA_chr03"
  if (!"site_chr" %in% names(dt))
    dt[, site_chr := sub(":.*", "", site)]

  # snp_chr: join from variant annotation on snp = snp_id
  if (!is.null(va) && "snp" %in% names(dt)) {
    if ("chr" %in% names(dt)) dt[, chr := NULL]   # drop any stale chr column
    va_sub <- va[, .(snp_id, snp_chr = chr)]
    dt <- merge(dt, va_sub, by.x = "snp", by.y = "snp_id", all.x = TRUE)
  } else {
    dt[, snp_chr := NA_character_]
  }

  # Add context column
  dt[, context := context]

  # Ensure beta exists
  if (!"beta" %in% names(dt)) dt[, beta := NA_real_]

  # Select and reorder final columns
  keep <- c("snp", "snp_chr", "snp_pos", "site", "site_chr", "site_pos",
            "context", "statistic", "beta", "pvalue", "p_FDR", "distance")
  keep <- intersect(keep, names(dt))
  dt[, .SD, .SDcols = keep]
}

############################################################
# 4) COLLECT RESULTS
############################################################

summary_rows    <- list()
sig_loose_list  <- list()   # p_FDR < 5e-8
sig_strict_list <- list()   # p_FDR < 1e-10

# Outer loops: tool → cohort → context (6 panels per tool)
for (tool in TOOLS) {
  rroot <- RESULT_ROOTS[[tool]]

  for (cohort in COHORTS) {
    long_list <- list()   # per-context long-format entries for all_sites table

    for (context in CONTEXTS) {
      panel_tag <- paste0(tolower(cohort), "_", tolower(context))
      rds_path  <- file.path(rroot, cohort, context, "cis_meqtl_all_results.rds")

      # Record missing panels in the summary so the output table is always complete
      if (!file.exists(rds_path)) {
        msg("Missing: ", rds_path, " — skipping")
        summary_rows[[paste(tool, panel_tag)]] <- data.table(
          tool = tool, cohort = cohort, context = context,
          cis_window_kb = CIS_WINDOW_KB,
          n_pairs = NA_integer_, n_sites_tested = NA_integer_,
          n_snps_tested = NA_integer_,
          n_sig_pairs_p1e10 = NA_integer_, n_sig_sites_p1e10 = NA_integer_,
          n_sig_snps_p1e10 = NA_integer_,
          lambda = NA_real_, status = "missing")
        next
      }

      msg("Reading: ", tool, " ", cohort, " ", context)
      raw_dt <- readRDS(rds_path)

      # Load SNP variant annotation (snp_id, chr, pos)
      va_path <- file.path(ANNOT_ROOT, cohort, context, "snp_variant_annot.rds")
      va <- if (file.exists(va_path)) as.data.table(readRDS(va_path)) else NULL

      # Harmonise column names and types across both tools
      dt <- clean_result(raw_dt, tool, context, va)
      rm(raw_dt); gc()

      # Per-panel summary statistics
      lambda         <- compute_lambda(dt$pvalue)
      n_sites        <- dt[, uniqueN(site)]
      n_snps         <- dt[, uniqueN(snp)]
      n_sig_l        <- dt[p_FDR < FDR_LOOSE,  uniqueN(site)]
      n_sig_s_sites  <- dt[p_FDR < FDR_STRICT, uniqueN(site)]
      n_sig_s_snps   <- dt[p_FDR < FDR_STRICT, uniqueN(snp)]
      n_sig_s_pairs  <- dt[p_FDR < FDR_STRICT, .N]

      msg("  ", panel_tag,
          " | pairs=", nrow(dt),
          " | sites=", n_sites,
          " | snps=", n_snps,
          " | lambda=", round(lambda, 4),
          " | sig_1e-10: pairs=", n_sig_s_pairs,
          "  sites=", n_sig_s_sites,
          "  snps=", n_sig_s_snps)

      # One row per panel in the comprehensive summary
      summary_rows[[paste(tool, panel_tag)]] <- data.table(
        tool              = tool,
        cohort            = cohort,
        context           = context,
        cis_window_kb     = CIS_WINDOW_KB,
        n_pairs           = nrow(dt),
        n_sites_tested    = n_sites,
        n_snps_tested     = n_snps,
        n_sig_pairs_p1e10 = n_sig_s_pairs,
        n_sig_sites_p1e10 = n_sig_s_sites,
        n_sig_snps_p1e10  = n_sig_s_snps,
        lambda            = round(lambda, 4),
        status            = "ok"
      )

      # QQ plots: TIFF only for render performance; EPS/PNG skipped (see note below).
      # Each plot uses all raw p-values (not just significant ones) to show the
      # full null distribution and quantify inflation via lambda.
      pv_qq <- dt$pvalue
      pv_qq <- pv_qq[is.finite(pv_qq) & pv_qq > 0]
      if (length(pv_qq) >= 2) {
        qq_file <- file.path(QQ_DIR,
          sprintf("qqplot_%s_%s_%s.tiff",
                  tolower(tool), tolower(cohort), tolower(context)))
        n_qq    <- length(pv_qq)
        # Colour distinguishes tools: red = GENESIS5 (LMM), blue = MatrixEQTL5 (LM)
        col_pt  <- if (grepl("GENESIS", tool)) "#e34a33" else "#2b8cbe"
        panel_labels <- c(
          BREEDING_CpG = "a)", BREEDING_CHG = "b)", BREEDING_CHH = "c)",
          NATURAL_CpG  = "d)", NATURAL_CHG  = "e)", NATURAL_CHH  = "f)"
        )
        qq_label <- panel_labels[paste0(cohort, "_", context)]
        draw_qq <- function() {
          # ppoints() generates uniform quantiles for expected distribution;
          # sort(pv_qq) orders observed values for the quantile-quantile comparison
          plot(-log10(ppoints(n_qq)), -log10(sort(pv_qq)),
               pch = 20, cex = 0.4, col = col_pt,
               xlab = "Expected -log10(p)", ylab = "Observed -log10(p)",
               main = "",
               cex.lab = 2.0, cex.axis = 2.0)
          mtext(qq_label, side = 3, adj = 0, line = 0.5, cex = 3.0, font = 1)
          abline(0, 1, col = "red", lty = 2)  # null expectation: observed = expected
        }
        # 2400×2400 px at 300 dpi gives an 8×8 inch print-ready panel
        tiff(qq_file, width = 2400, height = 2400, res = 300, compression = "lzw")
        draw_qq(); dev.off()
        pdf(sub("\\.tiff$", ".pdf", qq_file), width = 8, height = 8)
        draw_qq(); dev.off()
        # EPS and PNG skipped: EPS OOMs on large natural-context datasets (>1 GB vector
        # file); PNG render of millions of scatter points is prohibitively slow.
        # 14ab only needs the TIFF; 19ab exports PDF from native source.
      }

      # Significant-site tables (all SNP-site pairs passing threshold)
      sig_l <- dt[p_FDR < FDR_LOOSE][order(p_FDR)]
      sig_s <- dt[p_FDR < FDR_STRICT][order(p_FDR)]
      sig_loose_list[[paste(tool, panel_tag)]]  <- sig_l
      sig_strict_list[[paste(tool, panel_tag)]] <- sig_s

      # Long-format all-sites entry: minimum FDR per methylation site
      # (collapsing all SNP associations per site to its best p-value)
      min_fdr_dt <- dt[, .(
        site_chr = site_chr[1L],
        min_FDR  = min(p_FDR, na.rm = TRUE)
      ), by = site]
      min_fdr_dt[, context := context]
      long_list[[context]] <- min_fdr_dt

      rm(dt); gc()
    }  # context loop

    # Write long-format all-sites table once all three contexts are processed
    if (length(long_list) > 0) {
      long_dt <- rbindlist(long_list, fill = TRUE)
      setcolorder(long_dt, c("site", "site_chr", "context", "min_FDR"))
      setorder(long_dt, context, min_FDR)
      fwrite(long_dt,
             file.path(TAB_DIR,
               sprintf("all_sites_%s_%s.tsv", tolower(tool), tolower(cohort))),
             sep = "\t")
    }

    # Write significant-pair tables for each FDR threshold, combining all contexts
    cohort_key <- paste0(tool, " ", tolower(cohort))
    for (thresh_tag in c("p5e8", "p1e10")) {
      src  <- if (thresh_tag == "p5e8") sig_loose_list else sig_strict_list
      keys <- names(src)[grepl(cohort_key, names(src), ignore.case = TRUE)]
      if (!length(keys)) next
      combined <- rbindlist(src[keys], fill = TRUE)
      if (nrow(combined) > 0)
        fwrite(combined,
               file.path(SIG_DIR,
                 sprintf("sig_%s_%s_%s.tsv",
                         thresh_tag, tolower(tool), tolower(cohort))),
               sep = "\t")
    }
  }  # cohort loop
}  # tool loop

############################################################
# 5) WRITE COMPREHENSIVE SUMMARY
############################################################

# Bind all 12 panel summary rows (2 tools × 2 cohorts × 3 contexts)
summary_dt <- rbindlist(summary_rows, fill = TRUE)

fwrite(summary_dt, file.path(OUT_ROOT, "comprehensive_summary.tsv"), sep = "\t")

msg("Comprehensive summary saved: ", file.path(OUT_ROOT, "comprehensive_summary.tsv"))
cat("\nComprehensive summary (ok panels):\n")
print(summary_dt[status == "ok",
                 .(tool, cohort, context, cis_window_kb,
                   n_pairs, n_sites_tested, n_snps_tested,
                   n_sig_pairs_p1e10, n_sig_sites_p1e10, n_sig_snps_p1e10, lambda)])

msg("Comprehensive summary and sig-site tables done. Starting robust pair identification...")

############################################################
# 6) ROBUST PAIRS
#    A pair is robust when the same snp+site combination passes
#    FDR < 1e-10 independently in BOTH GENESIS5 and MatrixEQTL5.
#    SNP positions are taken from snp_variant_annot.rds (true genomic
#    bp coordinates), not from GDS-internal sequential indices.
#    Site positions are parsed from the site string (PA_chrXX:start-end).
############################################################

msg("Loading SNP position annotation maps for robust pair identification...")

snp_annot <- list()
for (cohort in COHORTS) {
  snp_annot[[cohort]] <- list()
  for (ctx in CONTEXTS) {
    va_path <- file.path(ANNOT_ROOT, cohort, ctx, "snp_variant_annot.rds")
    if (file.exists(va_path)) {
      va <- as.data.table(readRDS(va_path))
      snp_annot[[cohort]][[ctx]] <- va[, .(snp_id = as.integer(snp_id),
                                           snp_chr = as.character(chr),
                                           snp_pos = as.integer(pos))]
    } else {
      msg("  MISSING SNP annot: ", va_path)
      snp_annot[[cohort]][[ctx]] <- data.table()
    }
  }
}

robust_list    <- list()
robust_summary <- list()

for (cohort in COHORTS) {
  for (ctx in CONTEXTS) {
    key_g5  <- paste("GENESIS5",    paste0(tolower(cohort), "_", tolower(ctx)))
    key_me5 <- paste("MATRIXEQTL5", paste0(tolower(cohort), "_", tolower(ctx)))
    g5_dt   <- sig_strict_list[[key_g5]]
    me5_dt  <- sig_strict_list[[key_me5]]

    if (is.null(g5_dt) || !nrow(g5_dt) || is.null(me5_dt) || !nrow(me5_dt)) {
      msg("  ", cohort, "/", ctx, ": one tool has 0 pairs — skipping robust")
      next
    }

    keep_cols <- c("snp", "site", "statistic", "beta", "pvalue", "p_FDR")
    g5_keep  <- intersect(keep_cols, names(g5_dt))
    me5_keep <- intersect(keep_cols, names(me5_dt))
    g5_sub   <- g5_dt [, ..g5_keep]
    me5_sub  <- me5_dt[, ..me5_keep]

    setnames(g5_sub,  setdiff(g5_keep,  c("snp","site")),
             paste0(setdiff(g5_keep,  c("snp","site")), "_GENESIS5"))
    setnames(me5_sub, setdiff(me5_keep, c("snp","site")),
             paste0(setdiff(me5_keep, c("snp","site")), "_MATRIXEQTL5"))

    both <- merge(g5_sub, me5_sub, by = c("snp", "site"))
    both[, context := ctx]
    both[, cohort  := cohort]

    sa <- snp_annot[[cohort]][[ctx]]
    if (nrow(sa))
      both <- merge(both, sa, by.x = "snp", by.y = "snp_id", all.x = TRUE)
    else
      both[, `:=`(snp_chr = NA_character_, snp_pos = NA_integer_)]

    both[, site_chr := sub(":.*", "", site)]
    both[, site_pos := as.integer(sub("^[^:]+:(\\d+)-.*$", "\\1", site))]

    n_pairs <- nrow(both)
    n_snps  <- uniqueN(both$snp)
    n_sites <- uniqueN(both$site)
    msg("  ", cohort, "/", ctx, ": ", n_pairs, " robust pairs | ",
        n_snps, " unique SNPs | ", n_sites, " unique sites")

    robust_list[[paste0(cohort, "_", ctx)]] <- both
    robust_summary[[length(robust_summary) + 1]] <- data.table(
      cohort = cohort, context = ctx,
      robust_pairs        = n_pairs,
      robust_unique_snps  = n_snps,
      robust_unique_sites = n_sites
    )
  }
}

robust_dt         <- rbindlist(robust_list, fill = TRUE)
robust_summary_dt <- rbindlist(robust_summary)

col_order <- c("cohort", "context",
               "snp", "snp_chr", "snp_pos",
               "site", "site_chr", "site_pos",
               grep("_GENESIS5$",    names(robust_dt), value = TRUE),
               grep("_MATRIXEQTL5$", names(robust_dt), value = TRUE))
setcolorder(robust_dt, intersect(col_order, names(robust_dt)))

for (coh in COHORTS) {
  sub <- robust_dt[cohort == coh]
  out <- file.path(ROBUST_DIR, paste0("robust_markers_", tolower(coh), ".tsv"))
  fwrite(sub, out, sep = "\t")
  msg("  Saved: ", basename(out), " (", nrow(sub), " rows)")
}
fwrite(robust_summary_dt,
       file.path(ROBUST_DIR, "robust_context_summary.tsv"), sep = "\t")
msg("  Robust summary saved.")

############################################################
# 7) WRITE SUPPLEMENTARY TABLE S5 XLSX
############################################################

msg("Writing Supplementary Table S5...")

s5_export <- copy(robust_dt)
setnames(s5_export,
  old = c("p_FDR_GENESIS5", "p_FDR_MATRIXEQTL5"),
  new = c("FDR_GENESIS5",   "FDR_MATRIXEQTL5"),
  skip_absent = TRUE)

wb <- wb_workbook()
wb <- wb_add_worksheet(wb, "Supplementary Table S5")
wb <- wb_add_data(wb, sheet = "Supplementary Table S5",
  x = "Supplementary Table S5. Robust cis-meQTL pairs (FDR < 1e-10 in both GENESIS5 and MatrixEQTL5).",
  start_row = 1, start_col = 1, col_names = FALSE)
wb <- wb_add_data(wb, sheet = "Supplementary Table S5",
  x = paste0("snp_pos: true genomic bp position from snp_variant_annot.rds. ",
             "site_chr and site_pos: parsed from site string (chr:start-end). ",
             "Pair-level robust: same snp+site significant in both tools."),
  start_row = 2, start_col = 1, col_names = FALSE)
wb <- wb_add_data(wb, sheet = "Supplementary Table S5",
  x = as.data.frame(s5_export), start_row = 3, start_col = 1, col_names = TRUE)
wb_save(wb, S5_PATH)
msg("  Table S5 saved: ", S5_PATH, " (", nrow(s5_export), " rows)")

msg("Step 13ab finished. Outputs: ", OUT_ROOT)

sessionInfo()
```

---

### `14ab.tgc.joint.manhattan.plots.R`

```r
#!/usr/bin/env Rscript
############################################################
# TreeGeneClimate (TGC) — JOINT ECS + TBS
# Step 14ab: Circular Manhattan plots (GENESIS5 + MatrixEQTL5)
#            Version 5 — 18 cm canvas, BH-FDR axis option
#            + Combined 2-panel figures (BREEDING | NATURAL + legend)
#
# PURPOSE
# Produces circular Manhattan plots from cis-meQTL results (meQTL5 run).
# One plot per tool x cohort x axis mode (8 total):
#   GENESIS5     BREEDING/NATURAL  x  raw-p / BH-FDR
#   MATRIXEQTL5  BREEDING/NATURAL  x  raw-p / BH-FDR
# Three rings: CpG (outer), CHG (middle), CHH (inner).
# Dashed threshold line at BH-FDR < 1e-10.
#
# v5 changes vs original:
#   Canvas: 18 cm x 18 cm (was 24)     CHR_LABEL_CEX: 1.20 (was 0.975)
#   Dot cex: 0.22 (was 0.15)           Margins: c(0.5, 2, 1.0, 2)
#   Legend height: 2.0 cm (was 3.5)    FDR_AXIS flag (new)
#   Y-axis tick/title cex: 0.55/0.65   (larger caused overlap on 18 cm canvas)
#   TRACK_HEIGHT: 0.25 (was 0.18)      track.margin: c(0.005, 0.005)
#
# USAGE
#   Rscript --vanilla 14ab.tgc.joint.manhattan.plots.R
#   # or source from another script: SOURCED_14AB <- TRUE; source("14ab.R")
#
# INPUTS
#   RESULTS/JOINT/GENESIS5/MEQTL/<COHORT>/<CONTEXT>/cis_meqtl_all_results.rds
#   RESULTS/JOINT/MATRIXEQTL5/<COHORT>/<CONTEXT>/cis_meqtl_all_results.rds
#   RESULTS/JOINT/MQTL5/INPUTS/<COHORT>/<CONTEXT>/snp_variant_annot.rds
#   GDS files (chromosome layout)
#   DATA/METADATA/chrom.list  (optional; auto-built if absent)
#
# OUTPUTS (RESULTS/JOINT/COMBINED5/manhattan/)
#   manhattan_circular_<tool>_<cohort>_5.{tiff,pdf,svg,eps}      -- raw p-value
#   manhattan_circular_<tool>_<cohort>_5fdr.{tiff,pdf,svg,eps}   -- BH-FDR axis
#   legend_circular_5.tiff      (36 cm, for 2-panel combined figure)
#   legend_circular_5_wide.tiff (72 cm, for 4-panel comparison)
# OUTPUTS (RESULTS/JOINT/COMBINED5/panels/)
#   Figure3_circular_panel_genesis5_5.{tiff,pdf,svg,png,eps}
#   EDF4_circular_panel_matrixeqtl5_5.{tiff,pdf,svg,png,eps}
#   Figure3_circular_panel_genesis5_5_FDR.{tiff,pdf,svg,png,eps}
#   EDF4_circular_panel_matrixeqtl5_5_FDR.{tiff,pdf,svg,png,eps}
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(SNPRelate)
  library(gdsfmt)
  library(circlize)
  library(grid)
  library(RColorBrewer)
})

options(stringsAsFactors = FALSE)

############################################################
# 1) SETTINGS
############################################################

COHORTS        <- c("BREEDING", "NATURAL")
CTX_ORDER      <- c("CpG", "CHG", "CHH")   # outer -> inner rings
TOOLS          <- c("GENESIS5", "MATRIXEQTL5")
FDR_STRICT     <- 1e-10

# FDR_AXIS: when TRUE, y-axis shows -log10(BH-FDR); threshold constant at 10.
# Set before sourcing, or it defaults to FALSE.
if (!exists("FDR_AXIS")) FDR_AXIS <- FALSE

# v5 rendering settings
OUT_W_CM          <- 18; OUT_H_CM <- 18; OUT_RES <- 600
AXIS_GAP_FRAC     <- 0.80    # wide gap for y-axis label
AXIS_SECTOR_SCALE <- 1.50    # sector 14 (axis_panel2) width multiplier vs v5
TRACK_HEIGHT      <- 0.25
POINT_PALETTE     <- "Set1"
CHR_LABEL_CEX     <- 1.02    # 1.20 * 0.85 — reduced 15% to prevent "Un" clipping

# Combined panel dimensions (BREEDING | NATURAL side-by-side + legend strip)
LEG_H_CM  <- 2.0                          # legend strip height
COMB_W_CM <- OUT_W_CM * 2                 # 36 cm total width
COMB_H_CM <- OUT_W_CM + LEG_H_CM          # 20 cm total height

############################################################
# 2) PATHS  -- set PROJECT_ROOT and RDATA_DIR to your paths
############################################################

# === USER CONFIGURATION ===
# Set TGC_PROJECT_ROOT as an environment variable, or edit the fallback path below
PROJECT_ROOT <- Sys.getenv("TGC_PROJECT_ROOT",
  unset = "/path/to/your/project")
# ===========================

RESULT_ROOTS <- list(
  GENESIS5    = file.path(PROJECT_ROOT, "RESULTS", "JOINT", "GENESIS5",    "MEQTL"),
  MATRIXEQTL5 = file.path(PROJECT_ROOT, "RESULTS", "JOINT", "MATRIXEQTL5")
)

ANNOT_ROOT <- file.path(PROJECT_ROOT, "RESULTS", "JOINT", "MQTL5", "INPUTS")

OUT_MAN       <- file.path(PROJECT_ROOT, "RESULTS", "JOINT", "COMBINED5", "manhattan")
PANEL_DIR     <- file.path(PROJECT_ROOT, "RESULTS", "JOINT", "COMBINED5", "panels")
CORRECTED_DIR <- file.path(PROJECT_ROOT, "RESULTS", "CORRECTED", "FIGURES", "NEW")
LOG_DIR       <- file.path(PROJECT_ROOT, "RESULTS", "JOINT", "COMBINED5", "LOGS")
for (d in c(OUT_MAN, PANEL_DIR, CORRECTED_DIR, LOG_DIR))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

LOGFILE <- file.path(LOG_DIR, "step14ab.log")
if (file.exists(LOGFILE)) file.remove(LOGFILE)

# GDS files are in the ECS RDATA directory produced by step 12ab0
RDATA_DIR <- file.path(PROJECT_ROOT, "RESULTS", "ECS", "RANALYSIS", "RDATA")
GDS_FILES <- list(
  BREEDING = file.path(RDATA_DIR, "breeding.snp.gds"),
  NATURAL  = file.path(RDATA_DIR, "natural.snp.gds")
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
  ctx_colors <- setNames(
    brewer.pal(max(3, length(CTX_ORDER)), POINT_PALETTE)[seq_along(CTX_ORDER)],
    CTX_ORDER
  )

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

    var_path <- file.path(ANNOT_ROOT, cohort, ctx, "snp_variant_annot.rds")
    if (file.exists(var_path) && "snp" %in% names(res_dt)) {
      va <- as.data.table(readRDS(var_path))
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
        log_msg("  WARNING: snp_variant_annot.rds not found: ", var_path)
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
                                    open_device = TRUE, label_cex = 2.80) {

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
                 bg.border = "grey55", bg.col = "white",
      panel.fun = function(x, y) {
        sn <- CELL_META$sector.index

        if (sn == "axis_panel2") {
          circos.rect(CELL_META$xlim[1], CELL_META$ylim[1],
                      CELL_META$xlim[2], CELL_META$ylim[2],
                      col = "white", border = "white", lwd = 3)

          ax_x     <- CELL_META$xlim[2]
          at_vals  <- pretty(c(0, ylim_top), n = 3)
          at_vals  <- at_vals[at_vals >= 0 & at_vals <= ylim_top]
          tick_len <- CELL_META$xrange * 0.12

          circos.lines(c(ax_x, ax_x), c(0, ylim_top), lwd = 0.7, straight = TRUE)
          for (v in at_vals) {
            circos.lines(c(ax_x - tick_len, ax_x), c(v, v), lwd = 0.6, straight = TRUE)
            circos.text(ax_x - tick_len * 2.0, v,
                        labels = as.character(v),
                        adj = c(1, 0.5), cex = 0.55, font = 1, col = "black",
                        facing = "bending.inside", niceFacing = TRUE)
          }
          circos.text(ax_x - CELL_META$xrange * 0.80, ylim_top * 0.5,
                      labels = yaxis_label,
                      facing = "clockwise", niceFacing = FALSE,
                      adj = c(0.5, 0.5), cex = 0.65, col = "black")
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

  # Panel label: A = BREEDING, B = NATURAL
  fig_region  <- par("fig")
  fig_x1 <- fig_region[1]; fig_x2 <- fig_region[2]
  fig_y1 <- fig_region[3]; fig_y2 <- fig_region[4]
  label_x_npc <- fig_x1 + 0.02 * (fig_x2 - fig_x1)
  label_y_npc <- fig_y1 + 0.98 * (fig_y2 - fig_y1)
  upViewport(0)

  panel_label <- if (cohort == "BREEDING") "A" else "B"
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
log_msg("Step 14ab v5 -- Circular Manhattan plots (meQTL5)")
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
  suffix   <- if (fdr_mode) "_5fdr" else "_5"
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

ctx_colors_global <- setNames(
  brewer.pal(max(3, length(CTX_ORDER)), POINT_PALETTE)[seq_along(CTX_ORDER)],
  CTX_ORDER
)

draw_legend <- function(path, width_cm, height_cm = 2.0,
                        pt_cex = 2.0, text_cex = 1.5) {
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

leg_36 <- file.path(OUT_MAN, "legend_circular_5.tiff")
draw_legend(leg_36, width_cm = OUT_W_CM * 2)
log_msg("Legend (36 cm) saved: ", leg_36)

leg_72 <- file.path(OUT_MAN, "legend_circular_5_wide.tiff")
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
  axis_tag  <- if (fdr_mode) "_5fdr" else "_5"

  log_msg("  Loading data: ", tool, " BREEDING + NATURAL [FDR_AXIS=", fdr_mode, "]")
  dat_b <- load_panel_data(tool, "BREEDING")
  dat_n <- load_panel_data(tool, "NATURAL")

  draw_both <- function(dev_open_fn) {
    dev_open_fn()
    on.exit(dev.off(), add = TRUE)
    layout(matrix(c(1L, 2L, 3L, 3L), nrow = 2L, byrow = TRUE),
           widths  = lcm(c(OUT_W_CM, OUT_W_CM)),
           heights = lcm(c(OUT_W_CM, LEG_H_CM)))
    draw_circular_manhattan(dat_b, tool, "BREEDING",
                            open_device = FALSE, label_cex = 2.80)
    draw_circular_manhattan(dat_n, tool, "NATURAL",
                            open_device = FALSE, label_cex = 2.80)
    par(mar = c(0, 0, 0, 0))
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1))
    legend("center", legend = CTX_ORDER, pch = 16, col = ctx_colors_combined,
           pt.cex = 2.0, cex = 1.5, horiz = TRUE, bty = "n", x.intersp = 1.2)
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
    leg_tif <- file.path(OUT_MAN, "legend_circular_5.tiff")

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
      run_magick(paste("-density 300", shQuote(tiff_out), shQuote(eps_out)))
      log_msg("  EPS: ", basename(eps_out))
      unlink(TMP, recursive = TRUE)
    } else {
      log_msg("  WARNING: individual TIFFs missing — TIFF/PNG/EPS skipped")
    }
  }

  invisible(out_base)
}

figure_specs <- list(
  list(tool = "GENESIS5",    fdr = FALSE,
       stem = "Figure3_circular_panel_genesis5_5"),
  list(tool = "MATRIXEQTL5", fdr = FALSE,
       stem = "EDF4_circular_panel_matrixeqtl5_5"),
  list(tool = "GENESIS5",    fdr = TRUE,
       stem = "Figure3_circular_panel_genesis5_5_FDR"),
  list(tool = "MATRIXEQTL5", fdr = TRUE,
       stem = "EDF4_circular_panel_matrixeqtl5_5_FDR")
)

for (spec in figure_specs) {
  sep_line()
  log_msg("Panel: ", spec$stem)
  draw_combined_panel(spec$tool,
                      file.path(PANEL_DIR, spec$stem),
                      spec$fdr)
}

############################################################
# 10) COPY COMBINED PANELS TO CORRECTED/FIGURES/NEW
############################################################

sep_line()
log_msg("Copying combined panels to CORRECTED/FIGURES/NEW...")
copy_fmts <- c("tiff", "pdf", "svg", "png", "eps")
for (spec in figure_specs) {
  for (fmt in copy_fmts) {
    src_f <- file.path(PANEL_DIR, paste0(spec$stem, ".", fmt))
    dst_f <- file.path(CORRECTED_DIR, paste0(fmt, "_", spec$stem, ".", fmt))
    if (file.exists(src_f) && file.info(src_f)$size > 0) {
      file.copy(src_f, dst_f, overwrite = TRUE)
      log_msg("  ", basename(dst_f))
    } else {
      log_msg("  MISSING/EMPTY (skip): ", basename(src_f))
    }
  }
}

sep_line()
log_msg("Step 14ab finished")
log_msg("  Individual plots : ", OUT_MAN)
log_msg("  Combined panels  : ", PANEL_DIR)
log_msg("  Copies           : ", CORRECTED_DIR)
sep_line()

} # end SOURCED_14AB guard
```

---

### `14ab.tgc.joint.manhattan.plots.sh`

```bash
#!/bin/bash
#SBATCH -p scc-cpu
#SBATCH -t 06:00:00
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=96G
#SBATCH --job-name=TGC.14ab
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=vmchano.gaug@gmail.com

# 14ab.tgc.joint.manhattan.plots.sh
#
# Runs step 14ab in a single pass:
#   Part 1 — 8 individual circular Manhattan plots
#             (GENESIS5 + MatrixEQTL5) x (BREEDING + NATURAL) x (raw-p + BH-FDR)
#             Formats: TIFF, PDF, SVG, EPS + standalone legend TIFFs
#   Part 2 — 4 combined two-panel figures (BREEDING | NATURAL + legend strip)
#             Figure3_circular_panel_genesis5_5
#             EDF4_circular_panel_matrixeqtl5_5
#             Figure3_circular_panel_genesis5_5_FDR
#             EDF4_circular_panel_matrixeqtl5_5_FDR
#             Formats: TIFF (ImageMagick), PDF, SVG, PNG, EPS
#             Copies all to RESULTS/CORRECTED/FIGURES/NEW/
#
# Memory note: FDR_AXIS=TRUE on NATURAL/CHH (~500 M tests) requires ~96 GB RAM.

set -euo pipefail

module purge
module load gcc/14.2.0
module load r/4.5.2
module load imagemagick/7.1.1-39

export LC_ALL=C.UTF-8
export LANG=C.UTF-8
# export R_LIBS_USER="/path/to/your/Rlibs"  # uncomment and set if needed

PROJECT_ROOT="/path/to/your/project"  # <-- set this
SCRIPT="${PROJECT_ROOT}/SCRIPTS/JOINT/14ab.tgc.joint.manhattan.plots.R"

mkdir -p "${PROJECT_ROOT}/LOGS"
export TGC_PROJECT_ROOT="${PROJECT_ROOT}"

echo "============================================================"
echo "TGC — Circular Manhattan plots + Combined panels (step 14ab)"
echo "Node:  $(hostname)"
echo "Start: $(date)"
echo "RAM:   $(free -h | awk '/^Mem:/{print $2}')"
echo "============================================================"

Rscript --vanilla "${SCRIPT}"

echo "============================================================"
echo "All done: $(date)"

PANEL_DIR="${PROJECT_ROOT}/RESULTS/JOINT/COMBINED5/panels"
echo ""
echo "Combined panels:"
for stem in \
  "Figure3_circular_panel_genesis5_5" \
  "EDF4_circular_panel_matrixeqtl5_5" \
  "Figure3_circular_panel_genesis5_5_FDR" \
  "EDF4_circular_panel_matrixeqtl5_5_FDR"; do
  for fmt in tiff pdf svg png eps; do
    f="${PANEL_DIR}/${stem}.${fmt}"
    if [[ -f "${f}" ]]; then
      echo "  OK: ${stem}.${fmt}  ($(stat -c%s "${f}") bytes)"
    else
      echo "  MISSING: ${stem}.${fmt}"
    fi
  done
done

echo ""
echo "Copies in: ${PROJECT_ROOT}/RESULTS/CORRECTED/FIGURES/NEW/"
echo "============================================================"
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
# Supplementary (6 Venns, labelled a)–f)):
#   supp/venn_tools_<cohort>_<ctx>.{tiff,pdf}
#   Order: BREEDING/CpG (a), BREEDING/CHG (b), BREEDING/CHH (c),
#          NATURAL/CpG  (d), NATURAL/CHG  (e), NATURAL/CHH  (f)
#   Sets: GENESIS5 vs MatrixEQTL5 — unit: methylation SITES
#
# Main A — cohort comparison by context (3 Venns, a)–c)):
#   main/venn_cohorts_<ctx>.{tiff,pdf}
#   Sets: BREEDING vs NATURAL — unit: methylation SITES (robust)
#
# Main B — context comparison by cohort (2 Venns, a)–b)):
#   main/venn_contexts_<cohort>.{tiff,pdf}
#   Sets: CpG / CHG / CHH — unit: SNPs (robust)
#
# Tables (one per cohort):
#   tables/robust_markers_<cohort>.tsv
#
# USAGE
#   Rscript --vanilla 15ab.R   (or source() in RStudio)
#
# INPUTS
#   RESULTS/JOINT/COMBINED5/sig_sites/sig_p1e10_*.tsv  (from 13ab.R)
############################################################

# ---------------------------------------------------------------------------
# Auto-install missing packages (needed when running interactively)
# ---------------------------------------------------------------------------
local({
  lib <- Sys.getenv("R_LIBS_USER", unset = .libPaths()[1])
  dir.create(lib, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(lib, .libPaths()))
  for (pkg in c("ggvenn", "ggplot2", "patchwork", "data.table")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message("Installing missing package: ", pkg)
      install.packages(pkg, lib = lib, repos = "https://cloud.r-project.org",
                       quiet = TRUE)
    }
  }
})

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
TOOLS    <- c("GENESIS5", "MATRIXEQTL5")

# Per-plot fill colours — each plot within a group gets a distinct colour pair/triple.
# Names match original set identifiers (before display-renaming).

# All three Venn types use SNPs as the unit:
#   Supplementary : GENESIS5 SNPs vs MatrixEQTL5 SNPs  → overlap = robust SNPs
#   Main A        : BREEDING robust SNPs vs NATURAL robust SNPs (per context)
#   Main B        : CpG vs CHG vs CHH robust SNPs (per cohort)
# "Robust SNPs" = SNPs significant in both tools for the same cohort × context
# (SNP-level intersection; independent of which methylation site is associated).
# The separate robust MARKER TABLE (section 9) uses the stricter pair-level criterion
# (same SNP–site pair in both tools) as required for downstream annotation.

# Supplementary (6 plots, GENESIS5 vs MatrixEQTL5): Paired palette
SUPP_FILL_COLORS <- list(
  BREEDING_CpG = c(GENESIS5 = "#A6CEE3", MATRIXEQTL5 = "#1F78B4"),  # light/dark blue
  BREEDING_CHG = c(GENESIS5 = "#B2DF8A", MATRIXEQTL5 = "#33A02C"),  # light/dark green
  BREEDING_CHH = c(GENESIS5 = "#FB9A99", MATRIXEQTL5 = "#E31A1C"),  # light/dark red
  NATURAL_CpG  = c(GENESIS5 = "#FDBF6F", MATRIXEQTL5 = "#FF7F00"),  # light/dark orange
  NATURAL_CHG  = c(GENESIS5 = "#CAB2D6", MATRIXEQTL5 = "#6A3D9A"),  # light/dark purple
  NATURAL_CHH  = c(GENESIS5 = "#FFFF99", MATRIXEQTL5 = "#B15928")   # yellow/brown
)

# Main A — cohort comparison (3 plots, BREEDING vs NATURAL): distinct pairs per context
COHORT_FILL_COLORS <- list(
  CpG = c(BREEDING = "#FEE08B", NATURAL = "#D73027"),  # yellow / red
  CHG = c(BREEDING = "#91BFDB", NATURAL = "#4575B4"),  # light / dark blue
  CHH = c(BREEDING = "#D9F0D3", NATURAL = "#1B7837")   # light / dark green
)

# Main B — context comparison (2 plots, CpG/CHG/CHH): distinct triples per cohort
CTX_FILL_COLORS <- list(
  BREEDING = c(CpG = "#FC8D59", CHG = "#91CF60", CHH = "#91BFDB"),  # orange/green/blue
  NATURAL  = c(CpG = "#D7191C", CHG = "#1A9641", CHH = "#2C7BB6")   # dark: red/green/blue
)

# Supplementary: BREEDING first (a–c), NATURAL second (d–f)
SUPP_LABELS <- setNames(
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
PROJECT_ROOT <- Sys.getenv("TGC_PROJECT_ROOT",
  unset = "/path/to/your/project")
# ===========================

SIG_DIR  <- file.path(PROJECT_ROOT, "RESULTS", "JOINT", "COMBINED5", "sig_sites")
OUT_ROOT <- file.path(PROJECT_ROOT, "RESULTS", "JOINT", "COMBINED5", "overlap")
SUPP_DIR      <- file.path(OUT_ROOT, "supp")
MAIN_DIR      <- file.path(OUT_ROOT, "main")
TAB_DIR       <- file.path(OUT_ROOT, "tables")
PANEL_DIR     <- file.path(OUT_ROOT, "panels")
LOG_DIR       <- file.path(PROJECT_ROOT, "RESULTS", "JOINT", "COMBINED5", "LOGS")
CORRECTED_DIR <- file.path(PROJECT_ROOT, "RESULTS", "CORRECTED", "FIGURES", "NEW")

for (d in c(SUPP_DIR, MAIN_DIR, TAB_DIR, PANEL_DIR, LOG_DIR, CORRECTED_DIR))
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
  for (fmt in c("tiff", "pdf", "eps", "png", "svg")) {
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

# Human-readable labels for set names
SET_DISPLAY <- c(
  BREEDING    = "Breeding cohort",
  NATURAL     = "Natural cohort",
  GENESIS5    = "GENESIS",
  MATRIXEQTL5 = "MatrixEQTL"
)

# Draw a Venn diagram using ggvenn (ggplot2-native, true circles).
# Returns a ggplot object.
# set_list   : named list of character vectors (original names used for logic)
# fill_colors: named character vector of fill colours, keyed by original names
# panel_label: panel letter string, e.g. "a)"
draw_venn <- function(set_list, fill_colors, panel_label = NULL, sname_size = 7) {

  set_list <- set_list[lengths(set_list) > 0]
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
    lbl_up <- toupper(sub("[)]$", "", panel_label))
    p <- p + labs(tag = lbl_up) +
      theme(plot.tag = element_text(size = 22, face = "bold", hjust = 0))
  }

  p
}

############################################################
# 3b) OUTPUT DIMENSIONS (updated)
############################################################

# Supplementary (6 Venns, a)–f)):
#   supp/venn_tools_<cohort>_<ctx>.{tiff,pdf}
#   Sets: GENESIS5 vs MatrixEQTL5 — unit: SNPs
#
# Main A — cohort comparison by context (3 Venns, a)–c)):
#   main/venn_cohorts_<ctx>.{tiff,pdf}
#   Sets: BREEDING vs NATURAL — unit: SNPs (robust)
#
# Main B — context comparison by cohort (2 Venns, a)–b)):
#   main/venn_contexts_<cohort>.{tiff,pdf}
#   Sets: CpG / CHG / CHH — unit: SNPs (robust)

############################################################
# 4) LOAD SIGNIFICANT DATA
############################################################

log_msg("Step 15ab — Venn overlap analysis")
log_msg("ggvenn version: ", as.character(packageVersion("ggvenn")))
log_msg("Reading sig_p1e10 files from: ", SIG_DIR)

sig_raw <- list()
for (tool in TOOLS) {
  sig_raw[[tool]] <- list()
  for (cohort in COHORTS) {
    fname <- file.path(SIG_DIR,
      sprintf("sig_p1e10_%s_%s.tsv", tolower(tool), tolower(cohort)))
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
#     snp_variant_annot.rds has the ground-truth genomic positions.
#     GENESIS snp_pos in the sig files uses a different coordinate system
#     (GDS-internal), so we always look up positions from va.
############################################################

ANNOT_ROOT <- file.path(PROJECT_ROOT, "RESULTS", "JOINT", "MQTL5", "INPUTS")

va_map <- list()   # va_map[[cohort]][[ctx]]: named char vector  pos_string[snp_id]
for (cohort in COHORTS) {
  va_map[[cohort]] <- list()
  for (ctx in CONTEXTS) {
    va_path <- file.path(ANNOT_ROOT, cohort, ctx, "snp_variant_annot.rds")
    if (file.exists(va_path)) {
      va <- as.data.table(readRDS(va_path))
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
    g5  <- sig_raw[["GENESIS5"   ]][[cohort]]
    me5 <- sig_raw[["MATRIXEQTL5"]][[cohort]]

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

    g5_keys  <- g5_sub [, .(snp, site, FDR_GENESIS5    = p_FDR)]
    me5_keys <- me5_sub[, .(snp, site, FDR_MATRIXEQTL5 = p_FDR)]

    both <- merge(g5_keys, me5_keys, by = c("snp", "site"))
    both[, context := ctx]

    # Attach position columns from GENESIS5 if available
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
# 6) SUPPLEMENTARY — 6 Venns: GENESIS5 vs MatrixEQTL5 per cohort × context
#    Panel labels a)–f): BREEDING (a–c), NATURAL (d–f)
#    Unit: unique genomic positions (SNP chr:pos)
############################################################

# Storage for panel assembly (populated in sections 6-8)
supp_plots   <- list()   # 6 tool-comparison Venns  (key = "COHORT_ctx")
cohort_plots <- list()   # 3 cohort-comparison Venns (key = context)
ctx_plots    <- list()   # 2 context-comparison Venns (key = cohort)

log_msg("--- Supplementary: tool comparison (6 Venns) ---")

get_sites <- function(dt, ctx) {
  if (!nrow(dt) || !"context" %in% names(dt)) return(character(0))
  unique(dt[context == ctx, site])
}

for (cohort in COHORTS) {
  for (ctx in CONTEXTS) {
    key <- paste0(cohort, "_", ctx)
    lbl <- SUPP_LABELS[[key]]
    log_msg("  Panel ", lbl, "  [", cohort, " / ", ctx, "]")

    # Supplementary Venns: pair-level SNP counts.
    # The intersection shown = SNPs targeting the SAME methylation site in both tools
    # (robust pair-level), NOT just the SNP-level set intersection.
    # We construct synthetic integer sets so ggvenn computes the correct counts:
    #   set_A (GENESIS5)    : |A| = n_g5,  |A ∩ B| = n_rob
    #   set_B (MatrixEQTL5) : |B| = n_me5, |A ∩ B| = n_rob
    # robust IDs = 1..n_rob; G5-only = (n_rob+1)..(n_rob+n_g5_only)
    # ME5-only starts after G5-only to guarantee G5-only ∩ ME5-only = ∅
    g5_sub_ctx  <- sig_raw[["GENESIS5"   ]][[cohort]]
    me5_sub_ctx <- sig_raw[["MATRIXEQTL5"]][[cohort]]
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
    if (n_g5  > 0L) set_list[["GENESIS5"]]    <- set_a
    if (n_me5 > 0L) set_list[["MATRIXEQTL5"]] <- set_b

    log_msg("  Pair-level counts — G5: ", n_g5, " | ME5: ", n_me5,
            " | Robust: ", n_rob)

    base      <- file.path(SUPP_DIR,
      paste0("venn_tools_", tolower(cohort), "_", tolower(ctx)))
    base_corr <- file.path(CORRECTED_DIR,
      paste0("venn_tools_", tolower(cohort), "_", tolower(ctx)))

    gg <- draw_venn(set_list, SUPP_FILL_COLORS[[key]], lbl, sname_size = 5.5)
    supp_plots[[key]] <- gg
    save_plot(gg, base,      width_cm = 18, height_cm = 12)
    save_plot(gg, base_corr, width_cm = 18, height_cm = 12)
  }
}

############################################################
# 7) MAIN A — 3 Venns: BREEDING vs NATURAL per context (c)–e))
#    Unit: unique genomic positions (SNP chr:pos)
############################################################

log_msg("--- Main A: cohort comparison by context (3 Venns) ---")

for (ctx in CONTEXTS) {
  lbl <- CTX_PANEL_LABELS[[ctx]]
  log_msg("  Panel ", lbl, "  [", ctx, "]")

  breed_snps   <- robust_snp_pos[["BREEDING"]][[ctx]]
  natural_snps <- robust_snp_pos[["NATURAL" ]][[ctx]]

  set_list <- list()
  if (length(breed_snps))   set_list[["BREEDING"]] <- breed_snps
  if (length(natural_snps)) set_list[["NATURAL"]]  <- natural_snps

  base      <- file.path(MAIN_DIR,      paste0("venn_cohorts_", tolower(ctx)))
  base_corr <- file.path(CORRECTED_DIR, paste0("venn_cohorts_", tolower(ctx)))

  gg <- draw_venn(set_list, COHORT_FILL_COLORS[[ctx]], lbl, sname_size = 5.5)
  cohort_plots[[ctx]] <- gg
  save_plot(gg, base,      width_cm = 14, height_cm = 10)
  save_plot(gg, base_corr, width_cm = 14, height_cm = 10)
}

############################################################
# 8) MAIN B — 2 Venns: CpG vs CHG vs CHH per cohort (a)–b))
#    Unit: unique genomic positions (SNP chr:pos)
############################################################

log_msg("--- Main B: context comparison by cohort (2 Venns) ---")

for (cohort in COHORTS) {
  lbl <- COHORT_PANEL_LABELS[[cohort]]
  log_msg("  Panel ", lbl, "  [", cohort, "]")

  snp_sets <- list()
  for (ctx in CONTEXTS) {
    snps <- robust_snp_pos[[cohort]][[ctx]]
    if (length(snps)) snp_sets[[ctx]] <- snps
  }

  base      <- file.path(MAIN_DIR,      paste0("venn_contexts_", tolower(cohort)))
  base_corr <- file.path(CORRECTED_DIR, paste0("venn_contexts_", tolower(cohort)))

  gg <- draw_venn(snp_sets, CTX_FILL_COLORS[[cohort]], lbl)
  ctx_plots[[cohort]] <- gg
  save_plot(gg, base,      width_cm = 21, height_cm = 14)
  save_plot(gg, base_corr, width_cm = 21, height_cm = 14)
}

############################################################
# 8b) PANEL ASSEMBLY
#   Figure7 : row1 = ctx_plots (A,B — context comparison per cohort)
#             row2 = cohort_plots (C,D,E — cohort comparison per context)
#   SuppFig5: 2×3 grid, GENESIS5 vs MatrixEQTL5 per cohort×context
############################################################

log_msg("--- Panel assembly ---")

row1_valid <- Filter(Negate(is.null), ctx_plots)
row2_valid <- Filter(Negate(is.null), cohort_plots)

if (length(row1_valid) >= 1 && length(row2_valid) >= 1) {
  row1 <- Reduce(`+`, row1_valid) + plot_layout(ncol = length(row1_valid))
  row2 <- Reduce(`+`, row2_valid) + plot_layout(ncol = length(row2_valid))
  fig7 <- row1 / row2 + plot_layout(heights = c(1, 1))
  save_plot(fig7,
            file.path(PANEL_DIR, "Figure7_Venn_contexts_cohorts_panel"),
            width_cm = 42, height_cm = 28)
  save_plot(fig7,
            file.path(CORRECTED_DIR, "Figure7_Venn_contexts_cohorts_panel"),
            width_cm = 42, height_cm = 28)
  log_msg("  Figure7 panel saved → ", PANEL_DIR, " and ", CORRECTED_DIR)
} else {
  log_msg("  WARNING: insufficient Venn plots for Figure7 panel")
}

row1s <- Filter(Negate(is.null), supp_plots[paste0("BREEDING_", CONTEXTS)])
row2s <- Filter(Negate(is.null), supp_plots[paste0("NATURAL_",  CONTEXTS)])

if (length(row1s) >= 1 && length(row2s) >= 1) {
  row1 <- Reduce(`+`, row1s) + plot_layout(ncol = 3)
  row2 <- Reduce(`+`, row2s) + plot_layout(ncol = 3)
  fs5  <- row1 / row2 + plot_layout(heights = c(1, 1))
  save_plot(fs5,
            file.path(PANEL_DIR, "SuppFig5_Venn_tools_panel"),
            width_cm = 54, height_cm = 24)
  save_plot(fs5,
            file.path(CORRECTED_DIR, "SuppFig5_Venn_tools_panel"),
            width_cm = 54, height_cm = 24)
  log_msg("  SuppFig5 panel saved → ", PANEL_DIR, " and ", CORRECTED_DIR)
} else {
  log_msg("  WARNING: insufficient Venn plots for SuppFig5 panel")
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
      "context", "FDR_GENESIS5", "FDR_MATRIXEQTL5"),
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

############################################################
# 10) SESSION INFO
############################################################

writeLines(capture.output(sessionInfo()),
           file.path(LOG_DIR, "step15ab_sessionInfo.txt"))

log_msg("Step 15ab finished — outputs in: ", OUT_ROOT)
```

---

### `15ab.tgc.joint.venn.overlap.sh`

```bash
#!/bin/bash
#-------------------------------------------------------------------------------
# TreeGeneClimate (TGC) — JOINT ECS + TMS
# Step 15ab: meQTL overlap analysis — Venn diagrams + panel assembly
#
# Reads significant meQTL results (p_FDR < 1e-10) from 13ab.R and produces:
#   - 6 Venn diagrams (tool comparison per cohort×context; cohort vs context comparisons)
#   - SuppFig5 panel: 3×2 tool-comparison Venns
#   - Figure 6 panel: 1×5 contexts + cohorts Venns
#   - overlap_summary.tsv
#
# REQUIRES
#   R packages: data.table, ggvenn, ggplot2, patchwork
#   (installed automatically if absent)
#
# USAGE
#   sbatch 15ab.sh
#
# INPUTS
#   RESULTS/JOINT/COMBINED5/sig_sites/sig_p1e10_*.tsv  (from 13ab.R)
#
# OUTPUTS
#   RESULTS/JOINT/COMBINED5/overlap/
#-------------------------------------------------------------------------------

#SBATCH -p YOUR_PARTITION
#SBATCH -t 00:30:00
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH --job-name=TGC.15ab
#SBATCH --output=/path/to/your/project/LOGS/%x_%j.out
#SBATCH --error=/path/to/your/project/LOGS/%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=YOUR_EMAIL

set -euo pipefail

# ---------------------------------------------------------------------------
# Environment
# ---------------------------------------------------------------------------
module purge
module load gcc/14.2.0
module load r/4.5.2

export LC_ALL=C.UTF-8
export LANG=C.UTF-8
# export R_LIBS_USER="/path/to/your/Rlibs"  # uncomment and set if needed
[ -n "${R_LIBS_USER:-}" ] && mkdir -p "${R_LIBS_USER}"

# ---------------------------------------------------------------------------
# Install missing R packages (idempotent)
# ---------------------------------------------------------------------------
Rscript --vanilla -e '
  lib <- Sys.getenv("R_LIBS_USER")
  dir.create(lib, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(lib, .libPaths()))
  cran_needed <- c("data.table", "ggvenn", "ggplot2", "patchwork")
  for (pkg in cran_needed) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message("Installing: ", pkg)
      install.packages(pkg, lib = lib, repos = "https://cloud.r-project.org")
    }
  }
'

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
# === USER CONFIGURATION ===
PROJECT_ROOT="/path/to/your/project"  # <-- set this
# ===========================
SCRIPTS="${PROJECT_ROOT}/SCRIPTS/JOINT"

mkdir -p "${PROJECT_ROOT}/LOGS"

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
echo "============================================================"
echo "TGC — JOINT — Step 15ab — meQTL overlap analysis"
echo "Node:   $(hostname)"
echo "Start:  $(date)"
echo "============================================================"

Rscript --vanilla "${SCRIPTS}/15ab.R"

echo "============================================================"
echo "Step 15ab finished: $(date)"
echo "============================================================"
```

---

### `16ab.tgc.joint.meth.heritability.R`

```r
#!/usr/bin/env Rscript
############################################################
# TreeGeneClimate (TGC) — JOINT ECS + TMS
# Step 16ab: Per-site methylation heritability
#
# Estimates SNP-based heritability (h²_SNP) at every
# MEF-retained cytosine site using per-site GREML via
# eigendecomposition of the GRM (EMMA; Kang et al. 2008,
# Genetics 178:1709-1723). The GRM eigendecomposition is
# performed once per cohort; per-site REML is a fast
# one-dimensional optimisation over the variance ratio.
#
# For the breeding cohort, pedigree heritability (h²_ped)
# is additionally estimated as 2 × ICC from a full-sib
# family random-effects model (lme4 REML).
#
# INPUTS
#   <RDATA_DIR>/<cohort>_grm_gcta.rds         — dense GRM (n × n)
#   RESULTS/JOINT/MQTL5/INPUTS/<C>/<X>/
#     methylation_mvalues_matrix.rds            — M-value matrix
#   ECS/breeding_sample2family.txt             — sample → family map
#   RESULTS/JOINT/COMBINED5/overlap/tables/
#     robust_markers_{breeding,natural}.tsv    — meQTL site lists
#
# OUTPUTS (RESULTS/JOINT/HERITABILITY/)
#   h2_all_sites.tsv.gz      — per-site h²_SNP [and h²_ped] results
#   h2_summary.tsv           — aggregate stats per cohort × context
#   h2_meqtl_vs_all.*        — figure: meQTL sites vs genome background
#   h2_grm_vs_ped.*          — figure: h²_SNP vs h²_ped (breeding)
#   h2_combined_panel.*      — combined A/B panel figure
#
# NOTE: This script runs sequentially (all sites in memory per
# cohort × context). For very large datasets (>100k sites) use
# the companion SLURM script to distribute chunks across nodes;
# combine chunk TSVs, then skip to section 2 (aggregation).
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

options(stringsAsFactors = FALSE)

# === USER CONFIGURATION ===
# Set TGC_PROJECT_ROOT as an environment variable, or edit the fallback path below
PROJECT_ROOT <- Sys.getenv("TGC_PROJECT_ROOT",
  unset = "/path/to/your/project")
# ===========================
RDATA_DIR    <- file.path(PROJECT_ROOT, "RESULTS/ECS/RANALYSIS/RDATA")  # directory with <cohort>_grm_gcta.rds
# ===========================

COHORTS  <- c("BREEDING", "NATURAL")
CONTEXTS <- c("CpG", "CHG", "CHH")

H2_ROOT     <- file.path(PROJECT_ROOT, "RESULTS/JOINT/HERITABILITY")
ROBUST_FILE <- file.path(PROJECT_ROOT,
  "RESULTS/JOINT/COMBINED5/overlap/tables/robust_markers_breeding.tsv")
ROBUST_NAT  <- file.path(PROJECT_ROOT,
  "RESULTS/JOINT/COMBINED5/overlap/tables/robust_markers_natural.tsv")
FIG_DIR     <- file.path(PROJECT_ROOT, "RESULTS/DRAFT")

dir.create(H2_ROOT, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

msg <- function(...) {
  cat(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n"))
  flush(stdout())
}

############################################################
# CORE FUNCTIONS
############################################################

# Per-site REML via EMMA eigendecomposition (Kang et al. 2008)
# U, D_pos: pre-computed eigenvectors / positive eigenvalues of GRM
# Returns: h²_SNP = Vg / (Vg + Ve)
emma_reml <- function(y, U, D_pos) {
  na_idx <- is.na(y)
  if (any(na_idx)) y[na_idx] <- mean(y, na.rm = TRUE)
  y  <- y - mean(y)
  Uy <- drop(crossprod(U, y))
  n  <- length(D_pos)
  obj <- function(ld) {
    d  <- D_pos + exp(ld)
    vg <- sum(Uy^2 / d) / n
    if (vg <= 0) return(1e15)
    n * log(vg) + sum(log(d))
  }
  opt <- tryCatch(optimize(obj, c(-20, 20), tol = 1e-8), error = function(e) NULL)
  if (is.null(opt)) return(c(h2 = NA_real_, Vg = NA_real_, Ve = NA_real_))
  delta <- exp(opt$minimum)
  d     <- D_pos + delta
  Vg    <- sum(Uy^2 / d) / n
  Ve    <- Vg * delta
  Vp    <- Vg + Ve
  c(h2 = if (Vp > 1e-15) Vg / Vp else NA_real_, Vg = Vg, Ve = Ve)
}

# Pedigree h² via full-sib family random effect
# h²_ped ≈ 2 × ICC (assumes dominance and shared environment ≈ 0)
pedigree_h2 <- function(y, family_vec) {
  na_idx <- is.na(y)
  if (any(na_idx)) y[na_idx] <- mean(y, na.rm = TRUE)
  df  <- data.frame(y = y, fam = factor(family_vec))
  fit <- tryCatch(
    suppressWarnings(suppressMessages(
      lme4::lmer(y ~ 1 + (1 | fam), data = df, REML = TRUE,
                 control = lme4::lmerControl(
                   optimizer = "nloptwrap",
                   optCtrl   = list(maxfun = 2e5),
                   check.conv.singular = "ignore")))),
    error = function(e) NULL)
  if (is.null(fit))
    return(c(h2_ped = NA_real_, ICC = NA_real_, Vf = NA_real_, Ve_ped = NA_real_))
  vc  <- as.data.frame(lme4::VarCorr(fit))
  Vf  <- vc$vcov[vc$grp == "fam"]
  Ve2 <- vc$vcov[vc$grp == "Residual"]
  if (length(Vf) == 0) Vf <- 0
  Vp  <- Vf + Ve2
  ICC <- if (Vp > 1e-15) Vf / Vp else NA_real_
  c(h2_ped = if (!is.na(ICC)) 2 * ICC else NA_real_,
    ICC = ICC, Vf = Vf, Ve_ped = Ve2)
}

############################################################
# 1) PER-SITE HERITABILITY
############################################################

all_list <- list()

for (cohort in COHORTS) {
  DO_PEDIGREE <- (cohort == "BREEDING")

  msg("Loading GRM: ", cohort)
  K <- as.matrix(readRDS(
    file.path(RDATA_DIR, sprintf("%s_grm_gcta.rds", tolower(cohort)))))

  fam_dt <- NULL
  if (DO_PEDIGREE) {
    fam_file <- file.path(PROJECT_ROOT, "ECS/breeding_sample2family.txt")
    fam_dt   <- fread(fam_file, header = FALSE, col.names = c("sample", "family"))
  }

  for (ctx in CONTEXTS) {
    mval_file <- file.path(PROJECT_ROOT, "RESULTS/JOINT/MQTL5/INPUTS",
                           cohort, ctx, "methylation_mvalues_matrix.rds")
    if (!file.exists(mval_file)) {
      msg("  MISSING: ", mval_file, " — skipping"); next
    }

    Mmat   <- as.matrix(readRDS(mval_file))
    common <- intersect(rownames(K), rownames(Mmat))
    K_sub  <- K[common, common, drop = FALSE]
    M_sub  <- Mmat[common, , drop = FALSE]
    n      <- nrow(K_sub)
    msg("  ", cohort, "/", ctx, ": ", n, " samples | ", ncol(M_sub), " sites")

    eig   <- eigen(K_sub, symmetric = TRUE)
    pos   <- eig$values > 1e-10
    U     <- eig$vectors[, pos, drop = FALSE]
    D_pos <- eig$values[pos]
    rm(eig); gc()

    fam_vec <- if (DO_PEDIGREE)
      factor(fam_dt$family[match(common, fam_dt$sample)]) else NULL

    n_sites <- ncol(M_sub)
    res_lst <- vector("list", n_sites)
    for (i in seq_len(n_sites)) {
      y          <- M_sub[, i]
      gr_res     <- emma_reml(y, U, D_pos)
      res_lst[[i]] <- if (DO_PEDIGREE) c(gr_res, pedigree_h2(y, fam_vec)) else gr_res
    }

    dt <- as.data.table(do.call(rbind, res_lst))
    dt[, `:=`(site = colnames(M_sub), cohort = cohort, context = ctx)]
    all_list[[paste(cohort, ctx)]] <- dt
    rm(M_sub, K_sub, U, D_pos, res_lst); gc()
  }
  rm(K); gc()
}

h2_all <- rbindlist(all_list, use.names = TRUE, fill = TRUE)
setorder(h2_all, cohort, context, site)
out_gz <- file.path(H2_ROOT, "h2_all_sites.tsv.gz")
fwrite(h2_all, out_gz, sep = "\t", compress = "gzip")
msg("Saved: h2_all_sites.tsv.gz (", nrow(h2_all), " rows)")

############################################################
# 2) SUMMARY TABLE
############################################################

summary_dt <- h2_all[, .(
  n_sites      = .N,
  n_valid      = sum(!is.na(h2)),
  mean_h2      = mean(h2, na.rm = TRUE),
  median_h2    = median(h2, na.rm = TRUE),
  sd_h2        = sd(h2, na.rm = TRUE),
  pct_h2_gt0.3 = mean(h2 > 0.30, na.rm = TRUE) * 100,
  pct_h2_gt0.5 = mean(h2 > 0.50, na.rm = TRUE) * 100,
  mean_h2_ped  = if ("h2_ped" %in% names(.SD))
                   mean(h2_ped, na.rm = TRUE) else NA_real_
), by = .(cohort, context)]
setorder(summary_dt, cohort, context)
fwrite(summary_dt, file.path(H2_ROOT, "h2_summary.tsv"), sep = "\t")
print(summary_dt[, .(cohort, context, n_valid, mean_h2, pct_h2_gt0.3, mean_h2_ped)])

############################################################
# 3) FIGURES
############################################################

meqtl_sites <- character(0)
for (rf in c(ROBUST_FILE, ROBUST_NAT)) {
  if (file.exists(rf))
    meqtl_sites <- union(meqtl_sites, unique(fread(rf, select = "site")$site))
}
msg("Robust meQTL sites: ", length(meqtl_sites))
h2_all[, is_meqtl := site %in% meqtl_sites]

save_fig <- function(base, plot, width, height) {
  ggsave(paste0(base, ".tiff"), plot, width = width, height = height,
         dpi = 300, compression = "lzw")
  ggsave(paste0(base, ".pdf"),  plot, width = width, height = height)
  ggsave(paste0(base, ".png"),  plot, width = width, height = height, dpi = 300)
  tryCatch(ggsave(paste0(base, ".eps"), plot, width = width, height = height),
           error = function(e) msg("EPS skipped: ", conditionMessage(e)))
}

ctx_cols   <- c(CpG = "#E63946", CHG = "#457B9D", CHH = "#2A9D8F")
coh_labels <- c(BREEDING = "Breeding (n=209)", NATURAL = "Natural (n=393)")

h2_plot <- h2_all[!is.na(h2)]
h2_plot[, context_f := factor(context, levels = c("CpG", "CHG", "CHH"))]
h2_plot[, cohort_f  := factor(cohort, levels = c("BREEDING", "NATURAL"),
                               labels = coh_labels)]

# Panel A: h²_SNP — meQTL sites vs genome background
h2_meqtl <- h2_plot[, .(
  group = ifelse(is_meqtl, "Robust meQTL sites", "All sites"),
  h2, context_f, cohort_f)]

p_meqtl <- ggplot(h2_meqtl, aes(x = context_f, y = h2, fill = group, colour = group)) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.5,
              position = position_dodge(0.8), linewidth = 0.3) +
  geom_boxplot(width = 0.06, outlier.shape = NA, fill = NA,
               position = position_dodge(0.8), linewidth = 0.4, coef = 0) +
  facet_wrap(~cohort_f, nrow = 1) +
  scale_fill_manual(values   = c("All sites" = "grey60",
                                  "Robust meQTL sites" = "#E63946"), name = NULL) +
  scale_colour_manual(values = c("All sites" = "grey40",
                                  "Robust meQTL sites" = "#9B1A22"), name = NULL) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(x = "Methylation context", y = expression(h[SNP]^2)) +
  theme_bw(base_size = 11) +
  theme(strip.background = element_rect(fill = "grey90"),
        legend.position  = "bottom",
        panel.grid.minor = element_blank())
save_fig(file.path(H2_ROOT, "h2_meqtl_vs_all"), p_meqtl, width = 8, height = 4.5)

# Panel B: h²_SNP vs h²_ped (breeding cohort only)
breed_h2 <- h2_all[cohort == "BREEDING" & !is.na(h2) & !is.na(h2_ped)]
if (nrow(breed_h2) > 0) {
  breed_h2[, context_f := factor(context, levels = c("CpG", "CHG", "CHH"))]
  set.seed(42)
  breed_plot <- breed_h2[, .SD[sample(.N, min(.N, 50000L))], by = context_f]

  p_ped <- ggplot(breed_plot, aes(x = h2, y = h2_ped, colour = context_f)) +
    geom_point(alpha = 0.08, size = 0.4, stroke = 0) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                colour = "black", linewidth = 0.5) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
    facet_wrap(~context_f, nrow = 1) +
    scale_colour_manual(values = ctx_cols, guide = "none") +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(x = expression(h[SNP]^2 * " (GREML / GRM)"),
         y = expression(h[Ped]^2 * " (2 × ICC, full-sib families)")) +
    theme_bw(base_size = 11) +
    theme(strip.background = element_rect(fill = "grey90"),
          panel.grid.minor = element_blank())
  save_fig(file.path(H2_ROOT, "h2_grm_vs_ped"), p_ped, width = 8, height = 3.5)

  # Combined two-panel figure
  combined <- (p_meqtl + labs(tag = "A") +
                 theme(plot.tag = element_text(face = "bold", size = 14))) /
              (p_ped   + labs(tag = "B") +
                 theme(plot.tag = element_text(face = "bold", size = 14))) +
    plot_layout(heights = c(4.5, 3.5))
  save_fig(file.path(FIG_DIR, "h2_combined_panel"), combined, width = 8, height = 8)
  msg("Saved: h2_combined_panel.*")
} else {
  msg("Skipping GRM-vs-pedigree figure (no breeding pedigree data)")
}

msg("Step 16ab complete. Outputs in: ", H2_ROOT)
```

---

### `16ab.tgc.joint.meth.heritability.sh`

```bash
#!/bin/bash
############################################################
# TreeGeneClimate (TGC) — JOINT ECS + TMS
# Step 16ab: Methylation heritability — SLURM submission
#
# For small datasets: run the R script directly.
#   Rscript 16ab.tgc.joint.meth.heritability.R
#
# For large datasets (>100k sites per context): distribute
# across cluster nodes using the array below, then collect.
#
# STAGE 1 — submit array (one task per chunk × cohort × context):
#   bash 16ab.tgc.joint.meth.heritability.sh
#
# STAGE 2 — collect results after all array tasks finish:
#   sbatch --job-name=TGC.16ab.collect \
#          --wrap="Rscript 16ab.tgc.joint.meth.heritability.R" \
#          16ab.tgc.joint.meth.heritability.sh collect
############################################################

set -euo pipefail

# === USER CONFIGURATION ===
PROJECT_ROOT="/path/to/your/project"   # <-- set this
R_LIBS_USER="/path/to/your/Rlibs"     # <-- set this
PARTITION="YOUR_PARTITION"
ACCOUNT="YOUR_ACCOUNT"
EMAIL="YOUR_EMAIL"
# ===========================

SCRIPTS="${PROJECT_ROOT}/SCRIPTS/JOINT"
OUTDIR="${PROJECT_ROOT}/RESULTS/JOINT/HERITABILITY"
LOGS="${PROJECT_ROOT}/LOGS"
mkdir -p "${OUTDIR}" "${LOGS}"

CHUNK_SIZE=1000
MODE="${1:-submit}"

# ---- Collect mode: aggregate chunks and generate figures ----
if [ "${MODE}" = "collect" ]; then
  #SBATCH -p YOUR_PARTITION
  #SBATCH -t 01:00:00
  #SBATCH -N 1 -c 2
  #SBATCH --mem=32G
  #SBATCH --job-name=TGC.16ab.collect
  #SBATCH --output=YOUR_LOGS_DIR/h2_collect_%j.out
  #SBATCH --error=YOUR_LOGS_DIR/h2_collect_%j.err
  #SBATCH --mail-type=END,FAIL
  #SBATCH --mail-user=YOUR_EMAIL
  module purge
  module load r/4.5.2
  export R_LIBS_USER="${R_LIBS_USER}"
  export LC_ALL=C.UTF-8
  Rscript --vanilla "${SCRIPTS}/16ab.tgc.joint.meth.heritability.R"
  exit 0
fi

# ---- Submit mode: launch one SLURM array per cohort × context ----
# Site counts (adjust to match your dataset)
declare -A N_SITES
N_SITES["BREEDING_CpG"]=36662
N_SITES["BREEDING_CHG"]=89636
N_SITES["BREEDING_CHH"]=273738
N_SITES["NATURAL_CpG"]=66872
N_SITES["NATURAL_CHG"]=156957
N_SITES["NATURAL_CHH"]=499222

echo "TGC Step 16ab — heritability array submission"
echo "Chunk size: ${CHUNK_SIZE} sites"

for COHORT in BREEDING NATURAL; do
  for CONTEXT in CpG CHG CHH; do
    KEY="${COHORT}_${CONTEXT}"
    N="${N_SITES[${KEY}]}"
    N_CHUNKS=$(( (N + CHUNK_SIZE - 1) / CHUNK_SIZE ))
    OUTDIR_CTX="${OUTDIR}/${COHORT}/${CONTEXT}"
    mkdir -p "${OUTDIR_CTX}"
    echo "${COHORT}/${CONTEXT}: ${N} sites → ${N_CHUNKS} chunks"

    JOBID=$(sbatch --parsable \
      -p "${PARTITION}" \
      -t 00:30:00 -N 1 -c 1 --mem=8G \
      --job-name="TGC.h2.${COHORT:0:1}.${CONTEXT}" \
      --output="${LOGS}/h2_${COHORT}_${CONTEXT}_%A_%a.out" \
      --error="${LOGS}/h2_${COHORT}_${CONTEXT}_%A_%a.err" \
      --mail-type=END,FAIL --mail-user="${EMAIL}" \
      --array="1-${N_CHUNKS}%50" \
      --wrap="
        set -euo pipefail
        module purge && module load r/4.5.2
        export R_LIBS_USER='${R_LIBS_USER}'
        export LC_ALL=C.UTF-8
        CHUNK_START=\$(( (SLURM_ARRAY_TASK_ID - 1) * ${CHUNK_SIZE} + 1 ))
        CHUNK_END=\$(( SLURM_ARRAY_TASK_ID * ${CHUNK_SIZE} ))
        OUT_FILE='${OUTDIR_CTX}/chunk_'\${SLURM_ARRAY_TASK_ID}'.tsv'
        Rscript --vanilla '${SCRIPTS}/16ab.tgc.joint.meth.heritability.R' \
          ${COHORT} ${CONTEXT} \${CHUNK_START} \${CHUNK_END} \${OUT_FILE}
      ")
    echo "  Submitted: job ${JOBID} (${N_CHUNKS} tasks)"
  done
done

echo ""
echo "When all tasks complete, run:"
echo "  bash $0 collect"
```

---

### `17ab.tgc.joint.marker.annotation.R`

```r
#!/usr/bin/env Rscript
############################################################
# TreeGeneClimate (TGC) — JOINT ECS + TMS
# Step 17ab: Marker annotation against reference genome GFF3
#
# Annotates:
#   - ECS DAPC top-10 SNPs (per cohort x DF)
#   - TMS DAPC top-10 methylation sites (per cohort x context x DF)
#   - TMS KW SVMPs: breeding top-150 balanced + natural 6 formal SVMPs (step 8b)
#   - meQTL ROBUST markers: SNP + site positions significant at
#     p_FDR < 1e-10 in BOTH GENESIS5 AND MatrixEQTL5
#     (reads robust_markers_<cohort>.tsv from 15ab.R)
#
# Each marker is annotated with:
#   - gene_id, gene_start, gene_end, gene_strand (GFF3 gene feature)
#   - distance_bp (0 = overlaps gene; positive = nearest gene distance)
#   - annotation_class (genic / proximal_intergenic / distal_intergenic)
#   - ref, alt, af, dr2 (from imputed VCF; SNP markers only)
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
TOOLS    <- c("GENESIS5", "MATRIXEQTL5")
TOP_N    <- 20L   # number of top markers to retain per DAPC discriminant function

# Methylation classification thresholds (context-specific).
# CHG/CHH methylation is generally much lower than CpG in conifers.
BETA_HYPO  <- c(CpG = 0.30, CHG = 0.20, CHH = 0.10)
BETA_HYPER <- c(CpG = 0.70, CHG = 0.50, CHH = 0.30)

############################################################
# 2) PATHS
############################################################

# === USER CONFIGURATION ===
# Set TGC_PROJECT_ROOT as an environment variable, or edit the fallback path below
PROJECT_ROOT <- Sys.getenv("TGC_PROJECT_ROOT",
  unset = "/path/to/your/project")
# ===========================

# Picea abies v2.0 reference annotation (Nystedt et al. 2013 updated build)
GFF3_FILE <- file.path(PROJECT_ROOT,
  "REFERENCE/Pabies2.0/Picab02_230926_at01_all_sorted.gff3")

# Beagle-imputed VCFs (AF + DR2 available); used for allele-frequency annotation
VCF_FILES <- list(
  BREEDING = file.path(PROJECT_ROOT,
    "RESULTS/ECS/VCF_SPLIT/tgc.ecs.breeding.call.filt.maf05.snvs.poly.imputed.vcf.gz"),
  NATURAL  = file.path(PROJECT_ROOT,
    "RESULTS/ECS/VCF_SPLIT/tgc.ecs.natural.call.filt.maf05.snvs.poly.imputed.vcf.gz")
)

DAPC_ECS_ROOT <- file.path(PROJECT_ROOT, "RESULTS/ECS/RANALYSIS/TABLES/dapc_loadings")
DAPC_TMS_ROOT <- file.path(PROJECT_ROOT, "RESULTS/TMS/RANALYSIS/TABLES/dapc_loadings")
# Robust marker tables produced by 15ab.R
ROBUST_ROOT   <- file.path(PROJECT_ROOT, "RESULTS/JOINT/COMBINED5/overlap/tables")
SVMP_ROOT     <- file.path(PROJECT_ROOT, "RESULTS/TMS/RANALYSIS/TABLES/heatmap_markers_8B")

# Directory with mean methylation beta files (one per cohort × context).
# Expected filename: mean_beta_<cohort>_<context>.tsv  (columns: site, mean_beta)
# Computed automatically from M-value matrices in MQTL5_INPUTDIR if not present.
BETA_ROOT <- file.path(PROJECT_ROOT, "RESULTS/TMS/METHYLATION/mean_betas")

# M-value matrix inputs from 12ab0.R — used to compute mean betas on the fly.
# Path: MQTL5_INPUTDIR/<COHORT>/<CONTEXT>/methylation_mvalues_matrix.rds
# Rows = samples, columns = site IDs (chr:start-end).
MQTL5_INPUTDIR <- file.path(PROJECT_ROOT, "RESULTS/JOINT/MQTL5/INPUTS")

# eggNOG-mapper functional annotation (isoform-level; best mRNA selected per gene).
# Columns used: id (1), eggnog_description (5), eggnog_go (8), eggnog_KEGG_ko (9),
#               interpro_description/pfam (15), interpro_panther_description (20).
FUNC_ANNOT_FILE <- file.path(PROJECT_ROOT,
  "REFERENCE/Pabies2.0/Picab02_230926_at01_all_isoform_annotations_merged_sorted_non_redundant_panthers.tsv")

# List of mRNA IDs classified as TE-derived gene models (to build is_te_gene flag).
TE_IDS_FILE <- file.path(PROJECT_ROOT,
  "REFERENCE/Pabies2.0/TE_IDs_REMOVED.txt")

# Non-imputed (raw) VCFs — fallback for SNPs absent from the imputed VCF.
# GENESIS5 uses non-imputed GDS, so some SNPs may not appear in the imputed VCF.
RAW_VCF_FILES <- list(
  BREEDING = file.path(PROJECT_ROOT,
    "RESULTS/ECS/VCF_SPLIT/tgc.ecs.breeding.call.filt.maf05.snvs.poly.vcf.gz"),
  NATURAL  = file.path(PROJECT_ROOT,
    "RESULTS/ECS/VCF_SPLIT/tgc.ecs.natural.call.filt.maf05.snvs.poly.vcf.gz")
)

OUT_ROOT <- file.path(PROJECT_ROOT, "RESULTS/JOINT/ANNOTATION17")
TMP_DIR  <- file.path(OUT_ROOT, "tmp")
dir.create(TMP_DIR, recursive = TRUE, showWarnings = FALSE)

############################################################
# 3) HELPERS
############################################################

msg <- function(...) cat(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n"))

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
  dt[distance_bp < 0, `:=`(gene_id     = "no_gene_on_chrom",
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
    default =                                     "no_gene_on_chrom"
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

# Back-transform M-values to beta values and return per-site mean betas.
# Uses EPSILON_M = 0.5 matching 12ab0.R methylkit_to_mvalues().
# Inverse: beta = (100.5 * exp(M) - 0.5) / (100 * (1 + exp(M)))
# M-value matrix: rows = samples, columns = site IDs (chr:start-end).
EPSILON_M <- 0.5

compute_mean_beta <- function(cohort, ctx) {
  mat_file <- file.path(MQTL5_INPUTDIR,
    toupper(cohort), ctx, "methylation_mvalues_matrix.rds")
  if (!file.exists(mat_file)) {
    msg("  NOTE: M-value matrix not found: ", mat_file)
    return(NULL)
  }
  msg("  Computing mean betas from M-value matrix: ", cohort, " / ", ctx)
  mmat <- readRDS(mat_file)
  if (is.data.frame(mmat)) mmat <- as.matrix(mmat)
  # colMeans across samples gives the average M-value per site
  m_means    <- colMeans(mmat, na.rm = TRUE)
  beta_means <- (100.5 * exp(m_means) - EPSILON_M) /
                (100   * (1 + exp(m_means)))
  data.table(site = names(m_means), mean_beta = as.numeric(beta_means))
}

# Build gene-centric summary from the combined annotation table.
# func_annot and te_genes are joined to enrich the summary with functional info.
make_gene_summary <- function(dt, func_annot, te_genes) {
  # Rank annotation classes: genic > proximal > distal; report best class per gene
  ann_order <- c(genic = 1L, proximal_intergenic = 2L, distal_intergenic = 3L)
  dt_g <- dt[!is.na(gene_id) & gene_id != "no_gene_on_chrom"]
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
# 5) ECS DAPC — TOP-10 SNPs PER COHORT x DF
############################################################

msg("======================================================")
msg("ECS DAPC — top-10 SNPs per cohort x DF")

ecs_files <- list(
  BREEDING = file.path(DAPC_ECS_ROOT, "breeding_dapc_loadings_all_DF1_DF2.csv"),
  NATURAL  = file.path(DAPC_ECS_ROOT, "natural_dapc_loadings_all_DF1_DF2.csv")
)

ecs_list <- list()
for (cohort in COHORTS) {
  f <- ecs_files[[cohort]]
  if (!file.exists(f)) { msg("  Missing: ", f); next }
  dt <- fread(f, header = TRUE)
  # Standardise column names (CSV header: chromosome,position,DF,loading)
  setnames(dt, names(dt), c("chr", "pos", "DF", "loading"))
  dt[, abs_loading := abs(loading)]
  dt[, cohort := cohort]
  # Top N by |loading| within each DF; loading magnitude indicates discriminatory power
  top <- dt[order(-abs_loading), .SD[seq_len(min(TOP_N, .N))], by = DF]
  top[, abs_loading := NULL]
  ecs_list[[cohort]] <- top
}

ecs_top <- rbindlist(ecs_list, fill = TRUE)

if (nrow(ecs_top) > 0) {
  ecs_top[, marker_id := paste0(chr, ":", pos)]

  # GFF3 annotation via bedtools
  bed_path <- file.path(TMP_DIR, "ecs_dapc.bed")
  write_marker_bed(ecs_top, "chr", "pos", "marker_id", bed_path)
  ann <- run_bedtools_closest(bed_path, genes_bed_path)

  if (nrow(ann) > 0) {
    ecs_top <- merge(ecs_top,
                     ann[, .(marker_id, gene_id, gene_chr, gene_start,
                              gene_end, gene_strand, distance_bp)],
                     by = "marker_id", all.x = TRUE)
  }

  # VCF allele info (per cohort) — imputed VCF first, raw VCF fallback for NA
  vcf_parts <- lapply(COHORTS, function(co) {
    sub <- ecs_top[cohort == co]
    vi  <- query_vcf(data.table(chr = sub$chr, pos = as.integer(sub$pos)),
                     VCF_FILES[[co]],
                     file.path(TMP_DIR, paste0("ecs_vcf_", co)))
    if (nrow(vi) > 0)
      sub <- merge(sub, vi, by.x = c("chr","pos"), by.y = c("chr","pos"), all.x = TRUE)
    else
      sub[, c("ref","alt","af","dr2") := list(NA_character_, NA_character_,
                                               NA_real_, NA_real_)]
    sub <- fill_from_raw_vcf(sub, co,
             file.path(TMP_DIR, paste0("ecs_raw_vcf_", co)))
    sub
  })
  ecs_top <- rbindlist(vcf_parts, fill = TRUE)
  ecs_top <- add_annotation_class(ecs_top)

  # SNP-specific metadata; methylation fields set to NA for non-site markers
  ecs_top[, source                   := "ECS_DAPC"]
  ecs_top[, marker_type              := "SNP"]
  ecs_top[, context                  := NA_character_]
  ecs_top[, methylation_direction    := NA_character_]
  ecs_top[, mean_beta                := NA_real_]
  ecs_top[, meth_status_absolute     := NA_character_]
  ecs_top[, higher_methylation_cohort := NA_character_]

  ecs_top <- add_functional_annotation(ecs_top, func_annot, te_genes)
  fwrite(ecs_top, file.path(OUT_ROOT, "ecs_dapc_top20_annotated.tsv"), sep = "\t")
  msg("  Saved: ecs_dapc_top20_annotated.tsv (", nrow(ecs_top), " markers)")
} else {
  msg("  No ECS DAPC markers found.")
  ecs_top <- data.table()
}

############################################################
# 6) TMS DAPC — TOP-10 SITES PER COHORT x CONTEXT x DF
############################################################

msg("======================================================")
msg("TMS DAPC — top-10 methylation sites per cohort x context x DF")

tms_list <- list()
for (cohort in COHORTS) {
  for (ctx in CONTEXTS) {
    f <- file.path(DAPC_TMS_ROOT,
                   sprintf("TMS_DAPC_loadings_%s_%s_ALL.tsv",
                           tolower(cohort), tolower(ctx)))
    if (!file.exists(f)) { msg("  Missing: ", f); next }
    dt <- fread(f, header = TRUE)
    # TSV header: loc  chr  pos  DF  loading
    setnames(dt, names(dt), c("loc", "chr", "pos", "DF", "loading"))
    dt[, abs_loading := abs(loading)]
    dt[, cohort  := cohort]
    dt[, context := ctx]
    top <- dt[order(-abs_loading), .SD[seq_len(min(TOP_N, .N))], by = DF]
    top[, abs_loading := NULL]
    tms_list[[paste(cohort, ctx)]] <- top
  }
}

tms_top <- rbindlist(tms_list, fill = TRUE)

if (nrow(tms_top) > 0) {
  # GFF3 annotation via bedtools
  bed_path <- file.path(TMP_DIR, "tms_dapc.bed")
  write_marker_bed(tms_top, "chr", "pos", "loc", bed_path)
  ann <- run_bedtools_closest(bed_path, genes_bed_path)

  if (nrow(ann) > 0) {
    setnames(ann, "marker_id", "loc")
    tms_top <- merge(tms_top,
                     ann[, .(loc, gene_id, gene_chr, gene_start,
                              gene_end, gene_strand, distance_bp)],
                     by = "loc", all.x = TRUE)
  }

  # Methylation sites are not in the VCF — set SNP-specific columns to NA
  tms_top[, c("ref","alt","af","dr2") := list(NA_character_, NA_character_,
                                               NA_real_, NA_real_)]
  tms_top <- add_annotation_class(tms_top)
  tms_top[, source      := "TMS_DAPC"]
  tms_top[, marker_type := "methylation_site"]

  # Option 1 — absolute methylation status (via loading sign + optional beta files)
  # methylation_direction: sign of DAPC loading (positive/negative relative to DF axis)
  tms_top[, methylation_direction := ifelse(loading > 0, "positive", "negative")]
  tms_top <- add_meth_status(tms_top, site_col = "loc")
  # Cross-cohort comparison not applicable for DAPC (per-cohort analysis)
  tms_top[, higher_methylation_cohort := NA_character_]

  setnames(tms_top, "loc", "marker_id")

  tms_top <- add_functional_annotation(tms_top, func_annot, te_genes)
  fwrite(tms_top, file.path(OUT_ROOT, "tms_dapc_top20_annotated.tsv"), sep = "\t")
  msg("  Saved: tms_dapc_top20_annotated.tsv (", nrow(tms_top), " markers)")
} else {
  msg("  No TMS DAPC markers found.")
  tms_top <- data.table()
}

############################################################
# 6b) TMS KW SVMPs — SELECTED MARKERS FROM STEP 8b
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

  dt <- add_functional_annotation(dt, func_annot, te_genes)

  out_f <- file.path(OUT_ROOT,
    sprintf("tms_svmp_%s_annotated.tsv", tolower(cohort)))
  fwrite(dt, out_f, sep = "\t")
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
# 7) meQTL ROBUST MARKERS (p_FDR < 1e-10 in BOTH GENESIS5 AND MatrixEQTL5)
############################################################

msg("======================================================")
msg("meQTL robust markers (p_FDR < 1e-10 in GENESIS5 AND MatrixEQTL5)")

meqtl_results <- list()

for (cohort in COHORTS) {
  cohort_key <- tolower(cohort)

  # Robust marker tables were produced by 15ab.R (pair-level intersect of both tools)
  rob_file <- file.path(ROBUST_ROOT,
    sprintf("robust_markers_%s.tsv", cohort_key))

  if (!file.exists(rob_file)) {
    msg("  WARNING: ", basename(rob_file),
        " not found — run 15ab.R first")
    next
  }

  dt <- fread(rob_file, header = TRUE)

  # Correct snp_pos from variant annotation (the robust_markers TSV inherited
  # wrong positions from a metadata indexing bug in 12ab2.R, now fixed in that
  # script). Override here so all downstream BED/VCF/bedtools steps use the
  # correct coordinates.
  va_path <- file.path(MQTL5_INPUTDIR, cohort, "CpG", "snp_variant_annot.tsv")
  if (file.exists(va_path)) {
    va <- fread(va_path, header = TRUE, select = c("snp_id", "pos"))
    dt <- merge(dt,
                va[, .(snp_id, snp_pos_correct = pos)],
                by.x = "snp", by.y = "snp_id", all.x = TRUE)
    dt[!is.na(snp_pos_correct), snp_pos := snp_pos_correct]
    dt[, snp_pos_correct := NULL]
    msg("  snp_pos overridden from variant annotation for ", dt[!is.na(snp_pos), .N], " rows")
  } else {
    msg("  WARNING: variant annotation not found — snp_pos not corrected: ", va_path)
  }

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

    # Fill NA ref/alt/af from raw VCF (covers GENESIS5 SNPs absent from imputed
    # VCF, and genotyped markers whose INFO was stripped by Beagle imputation).
    snps <- fill_from_raw_vcf(snps, cohort,
      file.path(TMP_DIR, sprintf("robust_%s_snps_rawvcf", cohort_key)))

    # Join associated epimarker information (sites, contexts) aggregated above
    snps <- merge(snps, snp_site_map, by = "snp", all.x = TRUE)

    snps <- add_annotation_class(snps)
    snps[, cohort                    := cohort]
    snps[, tool                      := "GENESIS5+MATRIXEQTL5"]
    snps[, source                    := "robust_meQTL"]
    snps[, marker_type               := "SNP"]
    snps[, DF                        := NA_character_]
    snps[, loading                   := NA_real_]
    snps[, methylation_direction     := NA_character_]
    snps[, mean_beta                 := NA_real_]
    snps[, meth_status_absolute      := NA_character_]
    snps[, higher_methylation_cohort := NA_character_]

    snps <- add_functional_annotation(snps, func_annot, te_genes)
    out_f <- file.path(OUT_ROOT,
      sprintf("robust_%s_snps_annotated.tsv", cohort_key))
    fwrite(snps, out_f, sep = "\t")
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
    sites[, tool                     := "GENESIS5+MATRIXEQTL5"]
    sites[, source                   := "robust_meQTL"]
    sites[, marker_type              := "methylation_site"]
    sites[, DF                       := NA_character_]
    sites[, loading                  := NA_real_]
    sites[, methylation_direction    := NA_character_]

    # Option 1 — absolute methylation status per site per cohort
    sites <- add_meth_status(sites)

    # higher_methylation_cohort filled in cross-cohort step below
    sites[, higher_methylation_cohort := NA_character_]

    sites <- add_functional_annotation(sites, func_annot, te_genes)
    out_f <- file.path(OUT_ROOT,
      sprintf("robust_%s_sites_annotated.tsv", cohort_key))
    fwrite(sites, out_f, sep = "\t")
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
    breed_b <- breed_s[!is.na(mean_beta), .(site, context, beta_b = mean_beta)]
    nat_b   <- nat_s  [!is.na(mean_beta), .(site, context, beta_n = mean_beta)]
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
        fwrite(meqtl_results[[key]], out_f, sep = "\t")
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
# 7b) CROSS-REFERENCE: ECS DAPC top SNPs × robust meQTL results
#
# For each ECS DAPC top SNP, checks whether it also appears as a robust meQTL
# SNP (same cohort, matched by chr+pos). If so, fills in the associated
# epimarker columns (n_associated_sites, associated_sites, associated_contexts,
# associated_site_chrs, associated_site_pos) from the meQTL snp_site_map.
# These columns are part of std_cols and will carry through to the combined table.
############################################################

msg("======================================================")
msg("Cross-referencing ECS DAPC top SNPs with robust meQTL results...")

if (nrow(ecs_top) > 0 && length(meqtl_results) > 0) {
  ecs_top[, pos := as.integer(pos)]

  for (co in COHORTS) {
    key  <- paste(co, "snp")
    rob  <- meqtl_results[[key]]
    if (is.null(rob) || nrow(rob) == 0) next

    assoc_cols <- intersect(
      c("chr","pos","n_associated_sites","associated_sites",
        "associated_contexts","associated_site_chrs","associated_site_pos"),
      names(rob))
    lk <- unique(rob[, .SD, .SDcols = assoc_cols])
    lk[, pos := as.integer(pos)]

    ecs_co <- ecs_top[cohort == co, .(marker_id, chr, pos)]
    joined <- merge(ecs_co, lk, by = c("chr","pos"), all.x = TRUE)
    matched <- joined[!is.na(n_associated_sites)]

    if (nrow(matched) > 0) {
      msg("  ", co, ": ", nrow(matched), " / ", nrow(ecs_co),
          " DAPC top SNPs are also robust meQTL SNPs")
      # Use set() for in-place assignment by row index (avoids copy overhead)
      for (i in seq_len(nrow(matched))) {
        mid  <- matched$marker_id[i]
        idx  <- which(ecs_top$marker_id == mid & ecs_top$cohort == co)
        if (!length(idx)) next
        set(ecs_top, idx, "n_associated_sites",   matched$n_associated_sites[i])
        set(ecs_top, idx, "associated_sites",      matched$associated_sites[i])
        set(ecs_top, idx, "associated_contexts",   matched$associated_contexts[i])
        if ("associated_site_chrs" %in% names(matched))
          set(ecs_top, idx, "associated_site_chrs", matched$associated_site_chrs[i])
        if ("associated_site_pos"  %in% names(matched))
          set(ecs_top, idx, "associated_site_pos",  matched$associated_site_pos[i])
      }
    } else {
      msg("  ", co, ": no DAPC top SNPs match robust meQTL SNPs by position")
    }
  }

  # Overwrite the saved ECS DAPC annotation file with cross-reference added
  fwrite(ecs_top, file.path(OUT_ROOT, "ecs_dapc_top20_annotated.tsv"), sep = "\t")
  msg("  ECS DAPC annotation updated with meQTL epimarker cross-reference.")

  in_mqtl_n <- if ("n_associated_sites" %in% names(ecs_top))
    ecs_top[!is.na(n_associated_sites), .N] else 0L
  msg("  ECS DAPC SNPs with associated epimarkers: ", in_mqtl_n, " / ", nrow(ecs_top))
} else {
  msg("  Skipping: ecs_top or meqtl_results is empty.")
}

############################################################
# 8) COMBINED TABLE
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

if (nrow(ecs_top) > 0)
  all_list[["ECS_DAPC"]] <- prep_for_combine(
    ecs_top[, tool := NA_character_][, loading := loading])

if (nrow(tms_top) > 0)
  all_list[["TMS_DAPC"]] <- prep_for_combine(tms_top[, tool := NA_character_])

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
# 9) GENE-CENTRIC SUMMARY TABLE
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

# Clean up temp files
unlink(TMP_DIR, recursive = TRUE)

msg("======================================================")
msg("Step 17ab finished. Outputs: ", OUT_ROOT)
msg("======================================================")

############################################################
# SECTION 10) SUPPLEMENTARY TABLE 6 — FILL MISSING REF/ALT/AF
#
# Strategy:
#   1. Read Supplementary Table 6 xlsx (headers on row 3)
#   2. Collect unique (chr, pos, cohort) where ref is NA
#   3. Query per-cohort imputed VCF (has AF in INFO)
#   4. Fallback to raw per-cohort VCF (no AF)
#   5. Fallback to unfiltered all-samples VCF with bcftools fill-tags
#      (for positions filtered out of per-cohort VCFs due to MAF < 0.05)
#   6. Write corrected xlsx preserving all other sheets and formatting
#
# NOTE: CSI indexing is required for large Norway spruce chromosomes (>512 Mb)
############################################################

suppressPackageStartupMessages(library(openxlsx2))

msg("======================================================")
msg("Section 10: Supplementary Table 6 — fill missing ref/alt/af")
msg("======================================================")

NATGEN_DIR <- file.path(PROJECT_ROOT, "RESULTS/DRAFT/NATURE.GENETICS")

XLSX_S6_IN  <- file.path(NATGEN_DIR,
  "260819_Chano.etal.2026_tgc_supp.tables.xlsx")
XLSX_S6_OUT <- file.path(NATGEN_DIR,
  "260819_Chano.etal.2026_tgc_supp.tables_corrected_s8.xlsx")

VCF_S6 <- list(
  BREEDING = list(
    imputed = file.path(PROJECT_ROOT,
      "RESULTS/ECS/VCF_SPLIT/tgc.ecs.breeding.call.filt.maf05.snvs.poly.imputed.vcf.gz"),
    raw     = file.path(PROJECT_ROOT,
      "RESULTS/ECS/VCF_SPLIT/tgc.ecs.breeding.call.filt.maf05.snvs.poly.vcf.gz")
  ),
  NATURAL = list(
    imputed = file.path(PROJECT_ROOT,
      "RESULTS/ECS/VCF_SPLIT/tgc.ecs.natural.call.filt.maf05.snvs.poly.imputed.vcf.gz"),
    raw     = file.path(PROJECT_ROOT,
      "RESULTS/ECS/VCF_SPLIT/tgc.ecs.natural.call.filt.maf05.snvs.poly.vcf.gz")
  )
)
VCF_UNFILT  <- file.path(PROJECT_ROOT,
  "RESULTS/ECS/VARIANT.CALLING/tgc.ecs.allsamples.call.unfilt.snvs.renamed.vcf.gz")
SAMPS_B     <- file.path(PROJECT_ROOT,
  "RESULTS/ECS/VARIANT.CALLING2/SPLIT/all_breeding_samples.txt")
SAMPS_N     <- file.path(PROJECT_ROOT,
  "RESULTS/ECS/VARIANT.CALLING2/SPLIT/all_natural_samples.txt")
COHORT_SAMPS_S6 <- list(BREEDING = SAMPS_B, NATURAL = SAMPS_N)

TMP_S6 <- file.path(PROJECT_ROOT, "RESULTS/JOINT/COMBINED5/tmp_s8fix")
dir.create(TMP_S6, recursive = TRUE, showWarnings = FALSE)

BCFTOOLS_S6 <- Sys.which("bcftools")
if (!nzchar(BCFTOOLS_S6)) {
  msg("WARNING: bcftools not found — skipping Section 10 (load module bcftools)")
} else if (!file.exists(XLSX_S6_IN)) {
  msg("WARNING: Supplementary Table xlsx not found: ", XLSX_S6_IN)
} else {

  msg("bcftools: ", BCFTOOLS_S6)
  msg("Checking/creating CSI indices...")
  for (coh in names(VCF_S6)) {
    for (tp in names(VCF_S6[[coh]])) {
      vcf <- VCF_S6[[coh]][[tp]]
      csi <- paste0(vcf, ".csi")
      if (file.exists(vcf) && !file.exists(csi)) {
        msg("  Indexing (", coh, " ", tp, "): ", basename(vcf))
        system(paste(shQuote(BCFTOOLS_S6), "index -c", shQuote(vcf)))
      }
    }
  }

  msg("Reading Supplementary Table 6...")
  st6 <- as.data.table(read_xlsx(XLSX_S6_IN, sheet = "Supplementary Table 6",
                                  start_row = 3))
  st6 <- st6[, names(st6)[!is.na(names(st6)) & names(st6) != "NA_"], with = FALSE]
  st6[, .row_idx := .I]
  st6[, pos := as.integer(pos)]
  COL_REF <- which(names(st6) == "ref")
  COL_ALT <- which(names(st6) == "alt")
  COL_AF  <- which(names(st6) == "af")
  msg("  Rows: ", nrow(st6), " | NA ref: ", sum(is.na(st6$ref)))

  # Query a VCF that already has AF in INFO (imputed)
  query_vcf_s6 <- function(markers_dt, vcf_file, tmp_prefix, has_af = TRUE) {
    empty <- data.table(chr = character(), pos = integer(),
                        ref = character(), alt = character(), af = numeric())
    if (!file.exists(vcf_file)) return(empty)
    markers_dt <- unique(markers_dt[!is.na(chr) & !is.na(pos), .(chr, pos)])
    if (!nrow(markers_dt)) return(empty)
    reg_file <- paste0(tmp_prefix, ".bed")
    fwrite(markers_dt[, .(chr, start = pos - 1L, end = pos)],
           reg_file, sep = "\t", col.names = FALSE)
    out_file <- paste0(tmp_prefix, ".tsv")
    fmt <- if (has_af) "'%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/AF\\n'"
           else        "'%CHROM\\t%POS\\t%REF\\t%ALT\\t.\\n'"
    system(sprintf("%s query -R %s -f %s %s > %s 2>/dev/null",
                   shQuote(BCFTOOLS_S6), shQuote(reg_file), fmt,
                   shQuote(vcf_file), shQuote(out_file)))
    file.remove(reg_file)
    if (!file.exists(out_file) || file.info(out_file)$size == 0) return(empty)
    res <- fread(out_file, header = FALSE, sep = "\t", fill = TRUE,
                 col.names = c("chr","pos","ref","alt","af"), na.strings = ".")
    file.remove(out_file)
    res[, pos := as.integer(pos)]; res[, af := suppressWarnings(as.numeric(af))]
    unique(res, by = c("chr","pos"))
  }

  # Query the unfiltered all-samples VCF, computing AF per cohort via fill-tags
  query_vcf_unfilt_s6 <- function(markers_dt, cohort, tmp_prefix) {
    empty <- data.table(chr = character(), pos = integer(),
                        ref = character(), alt = character(), af = numeric())
    if (!file.exists(VCF_UNFILT)) return(empty)
    samps_file <- COHORT_SAMPS_S6[[cohort]]
    if (!file.exists(samps_file)) return(empty)
    markers_dt <- unique(markers_dt[!is.na(chr) & !is.na(pos), .(chr, pos)])
    if (!nrow(markers_dt)) return(empty)
    reg_file <- paste0(tmp_prefix, "_unfilt.bed")
    fwrite(markers_dt[, .(chr, start = pos - 1L, end = pos)],
           reg_file, sep = "\t", col.names = FALSE)
    out_file <- paste0(tmp_prefix, "_unfilt.tsv")
    cmd <- sprintf(
      "%s view --samples-file %s -R %s %s | %s plugin fill-tags -- -t AF | %s query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT\\t%%INFO/AF\\n' > %s 2>/dev/null",
      shQuote(BCFTOOLS_S6), shQuote(samps_file), shQuote(reg_file),
      shQuote(VCF_UNFILT),
      shQuote(BCFTOOLS_S6), shQuote(BCFTOOLS_S6),
      shQuote(out_file)
    )
    system(cmd); file.remove(reg_file)
    if (!file.exists(out_file) || file.info(out_file)$size == 0) return(empty)
    res <- fread(out_file, header = FALSE, sep = "\t", fill = TRUE,
                 col.names = c("chr","pos","ref","alt","af"), na.strings = ".")
    file.remove(out_file)
    res[, pos := as.integer(pos)]; res[, af := suppressWarnings(as.numeric(af))]
    unique(res, by = c("chr","pos"))
  }

  fill_cohort_s6 <- function(sub, cohort) {
    missing <- unique(sub[is.na(ref), .(chr, pos)])
    if (!nrow(missing)) { msg("  No missing in ", cohort); return(sub) }
    msg("  ", cohort, ": ", nrow(missing), " chr:pos to fill")

    # 1. Imputed VCF
    vi <- query_vcf_s6(missing, VCF_S6[[cohort]]$imputed,
                       file.path(TMP_S6, paste0(tolower(cohort), "_imp")), TRUE)
    msg("  Imputed: ", nrow(vi), " hits")
    if (nrow(vi)) {
      sub <- merge(sub, vi[, .(chr, pos, ref_q=ref, alt_q=alt, af_q=af)],
                   by=c("chr","pos"), all.x=TRUE, sort=FALSE)
      sub[is.na(ref) & !is.na(ref_q), `:=`(ref=ref_q, alt=alt_q)]
      sub[is.na(af)  & !is.na(af_q),  af := af_q]
      sub[, c("ref_q","alt_q","af_q") := NULL]
    }

    # 2. Raw per-cohort VCF
    still <- unique(sub[is.na(ref), .(chr, pos)])
    if (nrow(still)) {
      msg("  Still missing: ", nrow(still), " — trying raw VCF")
      vr <- query_vcf_s6(still, VCF_S6[[cohort]]$raw,
                         file.path(TMP_S6, paste0(tolower(cohort), "_raw")), FALSE)
      msg("  Raw: ", nrow(vr), " hits")
      if (nrow(vr)) {
        sub <- merge(sub, vr[, .(chr, pos, ref_q=ref, alt_q=alt)],
                     by=c("chr","pos"), all.x=TRUE, sort=FALSE)
        sub[is.na(ref) & !is.na(ref_q), `:=`(ref=ref_q, alt=alt_q)]
        sub[, c("ref_q","alt_q") := NULL]
      }
    }

    # 3. Unfiltered all-samples VCF with fill-tags (MAF-filtered-out positions)
    still2 <- unique(sub[is.na(ref), .(chr, pos)])
    if (nrow(still2)) {
      msg("  Still missing: ", nrow(still2), " — trying unfiltered VCF (", cohort, " fill-tags)")
      vu <- query_vcf_unfilt_s6(still2, cohort,
                                file.path(TMP_S6, paste0(tolower(cohort), "_unfilt")))
      msg("  Unfiltered: ", nrow(vu), " hits")
      if (nrow(vu)) {
        sub <- merge(sub, vu[, .(chr, pos, ref_q=ref, alt_q=alt, af_q=af)],
                     by=c("chr","pos"), all.x=TRUE, sort=FALSE)
        sub[is.na(ref) & !is.na(ref_q), `:=`(ref=ref_q, alt=alt_q)]
        sub[is.na(af)  & !is.na(af_q),  af := af_q]
        sub[, c("ref_q","alt_q","af_q") := NULL]
      }
    }

    msg("  After fill — NA ref: ", sum(is.na(sub$ref)))
    sub
  }

  st6_b  <- fill_cohort_s6(st6[cohort == "BREEDING"], "BREEDING")
  st6_n  <- fill_cohort_s6(st6[cohort == "NATURAL"],  "NATURAL")
  st6_o  <- st6[!cohort %in% c("BREEDING","NATURAL")]
  st6_fx <- rbind(st6_b, st6_n, st6_o, fill = TRUE)
  setorder(st6_fx, .row_idx); st6_fx[, .row_idx := NULL]

  msg("Summary — NA ref: ", sum(is.na(st6_fx$ref)),
      " | NA alt: ", sum(is.na(st6_fx$alt)),
      " | NA af: ", sum(is.na(st6_fx$af)))

  # Write corrected xlsx (preserve all other sheets/formatting)
  msg("Loading original workbook...")
  wb_s6 <- wb_load(XLSX_S6_IN)
  XLSX_DATA_START <- 4L
  wb_s6 <- wb_add_data(wb_s6, sheet = "Supplementary Table 6",
                        x = data.frame(ref = st6_fx$ref),
                        start_row = XLSX_DATA_START, start_col = COL_REF,
                        col_names = FALSE)
  wb_s6 <- wb_add_data(wb_s6, sheet = "Supplementary Table 6",
                        x = data.frame(alt = st6_fx$alt),
                        start_row = XLSX_DATA_START, start_col = COL_ALT,
                        col_names = FALSE)
  wb_s6 <- wb_add_data(wb_s6, sheet = "Supplementary Table 6",
                        x = data.frame(af = st6_fx$af),
                        start_row = XLSX_DATA_START, start_col = COL_AF,
                        col_names = FALSE)
  msg("Saving: ", basename(XLSX_S6_OUT))
  wb_save(wb_s6, XLSX_S6_OUT)
  msg("Saved: ", XLSX_S6_OUT)

  # Save corrected TSVs
  ANN17 <- file.path(PROJECT_ROOT, "RESULTS/JOINT/ANNOTATION17")
  fwrite(st6_fx[cohort == "BREEDING"],
         file.path(ANN17, "robust_breeding_snps_annotated_corrected.tsv"), sep="\t")
  fwrite(st6_fx[cohort == "NATURAL"],
         file.path(ANN17, "robust_natural_snps_annotated_corrected.tsv"),  sep="\t")
  msg("Corrected TSVs saved to ", ANN17)

  unlink(TMP_S6, recursive = TRUE)
}

sessionInfo()
```

---

### `17ab.tgc.joint.marker.annotation.sh`

```bash
#!/bin/bash
#-------------------------------------------------------------------------------
# TreeGeneClimate (TGC) — JOINT ECS + TMS
# Step 17ab: Marker annotation against reference genome GFF3
#
# Annotates three sets of markers with gene features from the GFF3:
#   1. ECS DAPC   — top-20 SNPs per cohort x DF (DF1, DF2)
#   2. TMS DAPC   — top-20 methylation sites per cohort x context x DF
#   3. meQTL      — significant SNP + methylation site positions
#                   (p_FDR < 1e-10, both tools, both cohorts, all contexts)
#
# For each marker the script reports:
#   - Overlapping gene (distance_bp = 0) or nearest gene with distance
#   - annotation_class: genic | proximal_intergenic | distal_intergenic
#   - REF, ALT, AF, DR2 from the imputed VCF (SNP markers only)
#
# REQUIRES (adjust module names for your HPC environment):
#   bedtools >= 2.27   (loaded below — CHECK MODULE NAME)
#   bcftools >= 1.19   (confirmed: module bcftools/1.19)
#   R >= 4.5           (confirmed: module r/4.5.2)
#
# USAGE
#   sbatch 17ab.sh
#
# INPUTS
#   RESULTS/ECS/RANALYSIS/TABLES/dapc_loadings/
#   RESULTS/TMS/RANALYSIS/TABLES/dapc_loadings/
#   RESULTS/JOINT/COMBINED5/overlap/tables/robust_markers_*.tsv  (needs 15ab.R first)
#   RESULTS/ECS/VCF_SPLIT/tgc.ecs.*.imputed.vcf.gz
#   REFERENCE/Pabies2.0/Picab02_230926_at01_all_sorted.gff3
#
# OUTPUTS (RESULTS/JOINT/ANNOTATION17/)
#   ecs_dapc_top20_annotated.tsv
#   tms_dapc_top20_annotated.tsv
#   robust_{breeding,natural}_{snps,sites}_annotated.tsv
#   all_markers_annotated.tsv
#-------------------------------------------------------------------------------

#SBATCH -p YOUR_PARTITION
#SBATCH -t 02:00:00
#SBATCH -N 1
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH --job-name=TGC.annot17
#SBATCH --output=/path/to/your/project/LOGS/%x_%j.out
#SBATCH --error=/path/to/your/project/LOGS/%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=YOUR_EMAIL

set -euo pipefail

# ---------------------------------------------------------------------------
# Environment
# ---------------------------------------------------------------------------
module purge
module load gcc/14.2.0
module load r/4.5.2
module load bcftools/1.19

# bedtools2/2.31.1 confirmed available (gcc/14.2.0 must be loaded first)
module load bedtools2/2.31.1

export LC_ALL=C.UTF-8
export LANG=C.UTF-8
# export R_LIBS_USER="/path/to/your/Rlibs"  # uncomment and set if needed

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
# === USER CONFIGURATION ===
PROJECT_ROOT="/path/to/your/project"  # <-- set this
# ===========================
SCRIPTS="${PROJECT_ROOT}/SCRIPTS/JOINT"

mkdir -p "${PROJECT_ROOT}/LOGS"

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
echo "============================================================"
echo "TGC — JOINT — Step 17ab — Marker GFF3 annotation"
echo "Node:   $(hostname)"
echo "CPUs:   ${SLURM_CPUS_ON_NODE:-2}"
echo "Memory: ${SLURM_MEM_PER_NODE:-?} MB"
echo "Start:  $(date)"
echo "============================================================"

# Verify bedtools and bcftools are reachable
if ! command -v bedtools &>/dev/null; then
  echo "ERROR: bedtools not found in PATH. Adjust the module name above." >&2
  exit 1
fi
if ! command -v bcftools &>/dev/null; then
  echo "ERROR: bcftools not found in PATH." >&2
  exit 1
fi

echo "bedtools: $(bedtools --version | head -1)"
echo "bcftools: $(bcftools --version | head -1)"

Rscript --vanilla "${SCRIPTS}/17ab.R"

echo "============================================================"
echo "Step 17ab finished: $(date)"
echo "============================================================"
```

