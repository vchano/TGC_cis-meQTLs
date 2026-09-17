#!/usr/bin/env Rscript
############################################################
# TreeGeneClimate (TGC) — ECS
# Step 10a: Extended Data Figure 1 — GRM heatmaps + IBD violins
#           (breeding and natural cohorts)
#
# OUTPUT (FIXED; DO NOT CHANGE):
#   /path/to/your/project/RESULTS/ECS/RANALYSIS/FIGURES/EDF1
#
# FIGURE FILES (saved individually):
#   EDF1a_breeding_GRM.tiff   (GRM heatmap, breeding)
#   EDF1b_breeding_IBD.tiff   (IBD violin, PI_HAT, breeding)
#   EDF1c_natural_GRM.tiff    (GRM heatmap, natural)
#   EDF1d_natural_IBD.tiff    (IBD violin, PI_HAT, natural)
#   ExtendedDataFig1_GRM_IBD_panel.tiff  (combined 2x2 panel, A-D)
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
PROJECT_ROOT <- Sys.getenv("TGC_PROJECT_ROOT",
  unset = "/path/to/your/project")
# ===========================
RANA_DIR     <- file.path(PROJECT_ROOT, "RESULTS/ECS/RANALYSIS")
RDATA_DIR    <- file.path(RANA_DIR, "RDATA")   # RDS objects written by step 8a

FIG_DIR      <- file.path(RANA_DIR, "FIGURES", "EDF1")      # <-- fixed
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

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
plot_grm_heatmap <- function(G, groups, palette_named, out_file, title_label = "A)",
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
         merge_legend = TRUE, show_annotation_legend = TRUE)
    # Panel label (e.g. "A)") placed in the top-left corner using grid coordinates
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
plot_pihat_violin <- function(df, out_file, title_label = "B)",
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

# ---- Combined 2×2 panel: A) GRM breeding, B) IBD breeding,
#                          C) GRM natural,  D) IBD natural ----
# Uses grid viewports to tile four sub-figures into a single 32×32 cm canvas.
# Panel layout matches Extended Data Fig. 1 exactly.
save_edf1_panel <- function(ht_b, ht_n, p_ibd_b, p_ibd_n, out_file,
                             dpi = 600, w_cm = 38, h_cm = 32) {
  draw_panel <- function() {
    grid.newpage()

    # A) GRM breeding — top-left quadrant (x: 0–0.5, y: 0.5–1.0)
    pushViewport(viewport(x = 0, y = 0.5, width = 0.5, height = 0.5,
                          just = c("left", "bottom")))
    draw(ht_b, newpage = FALSE,
         heatmap_legend_side = "right", annotation_legend_side = "right",
         merge_legend = TRUE, show_annotation_legend = TRUE)
    grid.text("A)", x = unit(0.35, "cm"), y = unit(1, "npc") - unit(0.35, "cm"),
              just = c("left", "top"), gp = gpar(fontsize = FS_LABEL))
    popViewport()

    # B) IBD breeding — top-right quadrant
    pushViewport(viewport(x = 0.5, y = 0.5, width = 0.5, height = 0.5,
                          just = c("left", "bottom")))
    grid.draw(ggplotGrob(p_ibd_b))  # convert ggplot to grob for viewport drawing
    grid.text("B)", x = unit(0.35, "cm"), y = unit(1, "npc") - unit(0.35, "cm"),
              just = c("left", "top"), gp = gpar(fontsize = FS_LABEL))
    popViewport()

    # C) GRM natural — bottom-left quadrant
    pushViewport(viewport(x = 0, y = 0, width = 0.5, height = 0.5,
                          just = c("left", "bottom")))
    draw(ht_n, newpage = FALSE,
         heatmap_legend_side = "right", annotation_legend_side = "right",
         merge_legend = TRUE, show_annotation_legend = TRUE)
    grid.text("C)", x = unit(0.35, "cm"), y = unit(1, "npc") - unit(0.35, "cm"),
              just = c("left", "top"), gp = gpar(fontsize = FS_LABEL))
    popViewport()

    # D) IBD natural — bottom-right quadrant
    pushViewport(viewport(x = 0.5, y = 0, width = 0.5, height = 0.5,
                          just = c("left", "bottom")))
    grid.draw(ggplotGrob(p_ibd_n))
    grid.text("D)", x = unit(0.35, "cm"), y = unit(1, "npc") - unit(0.35, "cm"),
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
# 7) EDF1 A / C : GRM heatmaps with family/stand annotation legends
# ==============================================================================
out_a <- file.path(FIG_DIR, "EDF1a_breeding_GRM.tiff")
out_c <- file.path(FIG_DIR, "EDF1c_natural_GRM.tiff")
res_a <- plot_grm_heatmap(G_b, grp_b_vec, colors.17, out_a, title_label = "A)", key_digits = 2)
res_c <- plot_grm_heatmap(G_n, grp_n_vec, colors.25, out_c, title_label = "C)", key_digits = 2,
                            w_cm = 24)

# ==============================================================================
# 8) EDF1 B / D : IBD violins using PLINK PI_HAT (16x16 cm)
# ==============================================================================
# Classify each pairwise comparison as Within or Between family/population
pihat_b <- prep_pairs_pihat(ibd_b, ann_b, "Family")
pihat_n <- prep_pairs_pihat(ibd_n, ann_n, "Population")

out_b <- file.path(FIG_DIR, "EDF1b_breeding_IBD.tiff")
out_d <- file.path(FIG_DIR, "EDF1d_natural_IBD.tiff")

# y_limits include slight negative slack to prevent clipping of "Between" violin tails
res_b <- plot_pihat_violin(pihat_b, out_b, title_label = "B)", y_limits = c(-0.05, 0.70))
res_d <- plot_pihat_violin(pihat_n, out_d, title_label = "D)", y_limits = c(-0.05, 0.70))

# ==============================================================================
# 9) EXTENDED DATA FIG. 1: combined 2×2 panel (GRM + IBD, breeding + natural)
# ==============================================================================
out_panel <- file.path(FIG_DIR, "ExtendedDataFig1_GRM_IBD_panel.tiff")
save_edf1_panel(res_a$ht, res_c$ht, res_b$plot, res_d$plot, out_panel,
                dpi = 600, w_cm = 32, h_cm = 32)

cat(
  "Saved figures in:\n  ", FIG_DIR, "\n\n",
  "  - ", basename(out_a), "\n",
  "  - ", basename(out_b), "\n",
  "  - ", basename(out_c), "\n",
  "  - ", basename(out_d), "\n",
  "  - ", basename(out_panel), "  [Extended Data Fig. 1]\n",
  sep = ""
)

sessionInfo()
