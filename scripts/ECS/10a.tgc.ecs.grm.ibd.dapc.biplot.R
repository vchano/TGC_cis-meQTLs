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
PROJECT_ROOT <- "/path/to/your/project"  # <-- set this
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
