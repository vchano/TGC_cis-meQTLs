#!/usr/bin/env Rscript
############################################################
# TreeGeneClimate (TGC) — TMS
# Step 7b: Figure 2 — Epigenomic PCA panel (breeding + natural cohorts)
#
# Save PCA as ONE panel (A-F) with shared legend per row:
#   top row = breeding (CpG, CHG, CHH) + legend at right of C)
#   bottom  = natural  (CpG, CHG, CHH) + legend at right of F)
# Panel letters are drawn inside plots (ggplot title).
# Epimarker IDs use genomic labels chr:start (chr:pos) instead of V123 / V1 etc.
#
# INPUT:
#   RESULTS/TMS/RANALYSIS/METHYLKIT_OBJECTS/
#     methylBase_<cohort>_<context>_cov5_50_mpg*_mef0.05.rds
#
# OUTPUT:
#   RESULTS/TMS/RANALYSIS/FIGURES/FIG2/
#     Figure2_PCA_panel_A-F.tiff (+ pdf, eps, png)
############################################################

suppressPackageStartupMessages({
  library(methylKit)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(scales)
  library(grid)
  library(patchwork)
})

options(stringsAsFactors = FALSE)

# ==============================================================================
# 1) PATHS
# ==============================================================================
# === USER CONFIGURATION ===
# Set TGC_PROJECT_ROOT as an environment variable, or edit the fallback path below
PROJECT_ROOT <- Sys.getenv("TGC_PROJECT_ROOT",
  unset = "/path/to/your/project")
# ===========================
rds_dir <- file.path(PROJECT_ROOT, "RESULTS/TMS/RANALYSIS/METHYLKIT_OBJECTS")
fig_dir <- file.path(PROJECT_ROOT, "RESULTS/TMS/RANALYSIS/FIGURES/FIG2")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

map_file_breeding <- file.path(PROJECT_ROOT, "DATA/METADATA/breeding_sample2family.txt")
map_file_natural  <- file.path(PROJECT_ROOT, "DATA/METADATA/natural_sample2pop.txt")

ensure_file <- function(p) if (!file.exists(p)) stop("Missing file: ", p, call. = FALSE)
ensure_file(map_file_breeding)
ensure_file(map_file_natural)

# ==============================================================================
# 2) PALETTES
# ==============================================================================
# Named vectors: consistent group-to-color mapping used across all PCA panels
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
# so that PCA loadings are interpretable rather than V1, V2, ...
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
# inflate artificial PCs).
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

  max_na_frac <- max_na_frac_by_ctx(ctx)

  # PCA on the full filtered site set
  pca <- run_pca(X, max_na_frac = max_na_frac)
  pve <- pca$pve
  scores <- as.data.frame(pca$pr$x[, 1:2, drop = FALSE])
  colnames(scores) <- c("PC1","PC2")
  scores$Group <- as.character(group_factor)
  p_pca <- make_pca_scatter(scores, pve, palette_named)

  list(
    cohort = cohort,
    ctx = ctx,
    pca_plot = p_pca
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
# 6) FIGURE 2 — PCA PANEL (A-F) with shared legend per row
# ==============================================================================
# Top row: breeding (A=CpG, B=CHG, C=CHH); bottom row: natural (D=CpG, E=CHG, F=CHH)
pca_a <- add_panel_label(res_b_cpg$pca_plot, "A)")
pca_b <- add_panel_label(res_b_chg$pca_plot, "B)")
pca_c <- add_panel_label(res_b_chh$pca_plot, "C)")

pca_d <- add_panel_label(res_n_cpg$pca_plot, "D)")
pca_e <- add_panel_label(res_n_chg$pca_plot, "E)")
pca_f <- add_panel_label(res_n_chh$pca_plot, "F)")

# Collect legend per row so breeding and natural families/populations have separate legends
pca_row1 <- (pca_a | pca_b | pca_c) + plot_layout(guides = "collect") & theme(legend.position = "right", legend.justification = "top")
pca_row2 <- (pca_d | pca_e | pca_f) + plot_layout(guides = "collect") & theme(legend.position = "right", legend.justification = "top")

pca_panel <- pca_row1 / pca_row2

out_pca_panel <- file.path(fig_dir, "Figure2_PCA_panel_A-F.tiff")
save_panel_tiff(pca_panel, out_pca_panel, w_cm = 34, h_cm = 26, dpi = 600)

# ==============================================================================
# 7) FINISH
# ==============================================================================
cat("\nDONE Step 7b.\n\n")
cat("Saved: ", out_pca_panel, "  [Figure 2]\n", sep = "")
sessionInfo()
