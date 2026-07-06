#!/usr/bin/env Rscript
############################################################
# TreeGeneClimate (TGC) — TBS
# Step 8b: Heatmaps (per cohort) using cov+MEF filtered epiloci
#
# INPUT (same as 7b):
#   RESULTS/TBS/RANALYSIS/METHYLKIT_OBJECTS/
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
#   RESULTS/TBS/RANALYSIS/FIGURES/HEATMAPS_8B/
#     Figure7a_HEATMAP_breeding_top200_KW_BH.tiff
#     Figure7b_HEATMAP_natural_top200_KW_BH.tiff
#
#   RESULTS/TBS/RANALYSIS/TABLES/heatmap_markers_8B/
#     TBS_8B_locus_tests_breeding_<ctx>.tsv
#     TBS_8B_locus_tests_natural_<ctx>.tsv
#     TBS_8B_selected_markers_breeding_top200.tsv
#     TBS_8B_selected_markers_natural_top200.tsv
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
PROJECT_ROOT <- "/path/to/your/project"  # <-- set this
# ===========================

rds_dir <- file.path(PROJECT_ROOT, "RESULTS/TBS/RANALYSIS/METHYLKIT_OBJECTS")

fig_dir <- file.path(PROJECT_ROOT, "RESULTS/TBS/RANALYSIS/FIGURES/HEATMAPS_8B")
tab_dir <- file.path(PROJECT_ROOT, "RESULTS/TBS/RANALYSIS/TABLES/heatmap_markers_8B")

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

  log_file <- file.path(tab_dir, sprintf("TBS_8B_LOG_%s.txt", tolower(cohort)))
  log_con <- file(log_file, open = "wt")
  on.exit(close(log_con), add = TRUE)

  writeLines(paste0("TBS Step 8b log — cohort: ", cohort), con = log_con)
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
    out_ctx_tsv <- file.path(tab_dir, sprintf("TBS_8B_locus_tests_%s_%s.tsv",
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
  out_sel_bal_tsv <- file.path(tab_dir, sprintf("TBS_8B_selected_markers_%s_top%d_BALANCED.tsv",
                                                tolower(cohort), top_n))
  write.table(selected_bal, out_sel_bal_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

  # ---- Selection B: OVERALL top_n by padj regardless of context
  selected_all <- select_top_loci_overall(df_all, top_n = top_n)
  out_sel_all_tsv <- file.path(tab_dir, sprintf("TBS_8B_selected_markers_%s_top%d_OVERALL.tsv",
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
              file.path(tab_dir, "TBS_8B_selected_markers_natural_formal_SVMPs.tsv"),
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
