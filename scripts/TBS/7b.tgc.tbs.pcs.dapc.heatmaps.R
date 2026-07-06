#!/usr/bin/env Rscript
############################################################
# TreeGeneClimate (TGC) — TBS
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
#   RESULTS/TBS/RANALYSIS/METHYLKIT_OBJECTS/
#     methylBase_<cohort>_<context>_cov5_50_mpg*_mef0.05.rds
#
# OUTPUT:
#   RESULTS/TBS/RANALYSIS/FIGURES/FIG5_FIG6/
#     SUPP_PCA_panel_a-f.tiff
#     Figure5_DAPC_panel_a-f.tiff
#
#   RESULTS/TBS/RANALYSIS/TABLES/dapc_loadings/
#     TBS_DAPC_loadings_<cohort>_<context>_ALL.tsv
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
PROJECT_ROOT <- "/path/to/your/project"  # <-- set this
# ===========================
rds_dir <- file.path(PROJECT_ROOT, "RESULTS/TBS/RANALYSIS/METHYLKIT_OBJECTS")
fig_dir <- file.path(PROJECT_ROOT, "RESULTS/TBS/RANALYSIS/FIGURES/FIG5_FIG6")
tab_dir <- file.path(PROJECT_ROOT, "RESULTS/TBS/RANALYSIS/TABLES/dapc_loadings")

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

  out_tsv <- file.path(tab_dir, sprintf("TBS_DAPC_loadings_%s_%s_ALL.tsv",
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
