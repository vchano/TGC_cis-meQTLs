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
