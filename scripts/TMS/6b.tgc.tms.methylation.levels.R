#!/usr/bin/env Rscript
############################################################
# TreeGeneClimate (TGC) — TMS
# Step 6b: Plot per-sample mean methylation (%) by group
#   - Breeding cohort: per family (17 families) -> panels a, c, e
#   - Natural cohort:  per stand  (25 stands)   -> panels b, d, f
#   - One plot per context (CpG/CHG/CHH)
#   - Kruskal-Wallis test per cohort x context and Dunn's post-hoc pairwise
#     comparisons (two-sided, Benjamini-Hochberg adjusted); letters from the
#     adjusted p-values (groups sharing a letter do not differ)
#
# STYLE:
#   - Panel titles: single bold lowercase letter (a–f), no extra text
#   - No legend
#   - Kruskal-Wallis p-value as subtitle; post-hoc letters above the boxes
#   - Posthoc letters angled; spaced above boxes to reduce overlap
#
# PANEL LAYOUT:
#   Top row:    a (breeding CpG)  |  b (natural CpG)
#   Middle row: c (breeding CHG)  |  d (natural CHG)
#   Bottom row: e (breeding CHH)  |  f (natural CHH)
#
# INPUT:
#   RDS from Step 5b2 (after unite, before the methylated-fraction filter):
#     RESULTS/TMS/RANALYSIS/METHYLKIT_OBJECTS/methylBase_<cohort>_<ctx>_cov10_p99.9_mpg<k>.rds
#   DATA/METADATA/breeding_sample2family.txt, natural_sample2pop.txt
#
# OUTPUT:
#   Figures:
#     RESULTS/TMS/RANALYSIS/FIGURES/METHYLATION_LEVELS/
#   Logs, posthoc tables and summary of the tests:
#     RESULTS/TMS/RANALYSIS/METHYLATION.LEVELS/
#       methylation_level_tests_summary.tsv, <cohort>_<ctx>_posthoc.tsv
#   Values shown in the panel:
#     RESULTS/TMS/RANALYSIS/FIGURES/METHYLATION_LEVELS/methylation_levels_panel_values.tsv
#
# R packages: methylKit, ggplot2, dunn.test, multcompView, patchwork
############################################################

## ---------------------------
## Packages
## ---------------------------
suppressPackageStartupMessages({
  library(methylKit)
  library(ggplot2)
  library(multcompView)
  library(patchwork)
  library(dunn.test)    # Dunn's post-hoc test
})

## ---------------------------
## Paths
## ---------------------------
# === USER CONFIGURATION ===
# Set TGC_PROJECT_ROOT as an environment variable, or edit the fallback path below
PROJECT_ROOT <- Sys.getenv("TGC_PROJECT_ROOT", unset = "/path/to/your/project")
# ===========================

rds_dir  <- file.path(PROJECT_ROOT, "RESULTS/TMS/RANALYSIS/METHYLKIT_OBJECTS")
fig_dir  <- file.path(PROJECT_ROOT, "RESULTS/TMS/RANALYSIS/FIGURES/METHYLATION_LEVELS")
log_dir  <- file.path(PROJECT_ROOT, "RESULTS/TMS/RANALYSIS/METHYLATION.LEVELS")

map_file_breeding <- file.path(PROJECT_ROOT, "DATA/METADATA/breeding_sample2family.txt")
map_file_natural  <- file.path(PROJECT_ROOT, "DATA/METADATA/natural_sample2pop.txt")

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)

log_file <- file.path(log_dir, "Step6b_methylation_level_stats_posthoc.log")

## ---------------------------
## Colour codes (families / stands)
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
  
  # with the methylated-fraction filter
  pat_mef <- sprintf("^methylBase_%s_%s_cov10_p99\\.9_mpg[0-9]+_mef0\\.05\\.rds$", cohort, ctx)
  # cov+unite only
  pat_cov <- sprintf("^methylBase_%s_%s_cov10_p99\\.9_mpg[0-9]+\\.rds$", cohort, ctx)
  
  if (prefer_mef) {
    hits <- list.files(rds_dir, pattern = pat_mef, full.names = TRUE)
    if (length(hits) == 0) {
      # fall back to cov-only
      hits <- list.files(rds_dir, pattern = pat_cov, full.names = TRUE)
    }
  } else {
    hits <- list.files(rds_dir, pattern = pat_cov, full.names = TRUE)
    if (length(hits) == 0) {
      # fall back to the methylated-fraction-filtered file
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

# Groups are compared with rank-based tests because per-sample mean methylation is
# not normally distributed within groups (Shapiro-Wilk test on the residuals from
# the group means, reported): Kruskal-Wallis test, then Dunn's test for all pairs
# of groups (two-sided p-values, Benjamini-Hochberg adjustment).
run_group_tests <- function(df) {
  df <- df[!is.na(df$Group) & !is.na(df$Methylation), ]
  df$Group <- droplevels(df$Group)

  resid <- df$Methylation - ave(df$Methylation, df$Group)
  sh_p  <- tryCatch(shapiro.test(resid)$p.value, error = function(e) NA_real_)
  kw    <- kruskal.test(Methylation ~ Group, data = df)
  kw_p  <- kw$p.value

  invisible(capture.output(
    dn <- dunn.test::dunn.test(df$Methylation, df$Group, method = "bh", altp = TRUE,
                               kw = FALSE, table = FALSE, list = FALSE)
  ))
  pairs <- data.frame(comparison = dn$comparisons, Z = dn$Z,
                      p = dn$altP, p_adj = dn$altP.adjusted)
  p_adj <- setNames(pairs$p_adj, gsub(" - ", "-", pairs$comparison, fixed = TRUE))
  letters <- multcompView::multcompLetters(p_adj, threshold = 0.05)$Letters

  list(
    method    = "Kruskal-Wallis + Dunn (BH)",
    shapiro_p = sh_p,
    global_p  = kw_p,
    kw_H      = unname(kw$statistic),
    kw_df     = unname(kw$parameter),
    posthoc   = pairs,
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

# Box plot of per-sample means by group
make_boxplot_clean <- function(df, color_vec, title_letter,
                               ylab = "Percent methylation",
                               letter_angle = 45,
                               letter_yfactor = 0.12) {
  
  present <- levels(df$Group)
  color_vec <- color_vec[names(color_vec) %in% present]
  
  ggplot(df, aes(x = Group, y = Methylation, fill = Group)) +
    geom_boxplot(outlier.shape = NA, width = 0.7) +
    geom_jitter(position = position_jitter(width = 0.2, seed = 1), size = 1, color = "black", alpha = 0.7) +
    scale_fill_manual(values = color_vec, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.03, 0.18))) +   # room for the letters
    coord_cartesian(clip = "off") +
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
                              letter_yfactor = 0.12,
                              letter_size = 5) {
  
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
  
  p <- p + labs(subtitle = fmt_p(tests$global_p)) +
    theme(plot.subtitle = element_text(size = 13, face = "italic", hjust = 0))
  
  # Posthoc letters angled + spaced
  letters_vec <- tests$letters
  letters_vec <- letters_vec[levels(df$Group)]
  
  box_max <- aggregate(Methylation ~ Group, data = df, max, na.rm = TRUE)
  box_max$Letter <- unname(letters_vec[as.character(box_max$Group)])
  box_max$y <- box_max$Methylation + letter_yfactor * y_span
  
  p <- p + geom_text(
    data = box_max,
    aes(x = Group, y = y, label = Letter),
    size = letter_size, fontface = 1, vjust = 0,
    angle = letter_angle
  )
  
  # Save figure (all four formats independently)
  out_fig <- file.path(fig_dir, sprintf("%s_%s_boxplot.tiff", fig_prefix, ctx))
  save_all_formats(p, out_fig, w_cm = 26, h_cm = 14)
  
  # Save the post-hoc table (Dunn's test, BH-adjusted)
  posthoc_tsv <- file.path(log_dir, sprintf("%s_%s_posthoc.tsv", fig_prefix, ctx))
  write.table(tests$posthoc, posthoc_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
  
  list(df = df, plot = p, tests = tests, fig = out_fig, posthoc = posthoc_tsv, rds = rds_path)
}

## ---------------------------
## RUN ALL 6
## ---------------------------
sink(log_file)
cat("TGC — TMS — Step 6b\n")
cat("Timestamp: ", format(Sys.time()), "\n\n", sep = "")

cat("RDS dir:  ", rds_dir, "\n", sep = "")
cat("Fig dir:  ", fig_dir, "\n", sep = "")
cat("Log dir:  ", log_dir, "\n\n", sep = "")

# Global methylation levels must be estimated from all covered loci (cov+unite),
# NOT from the methylated-fraction-filtered subset, which drops sites unmethylated in all trees
# and would inflate the means; that subset is used for PCA (7b), SVMPs (8b) and 11ab.
prefer_mef <- FALSE

# BREEDING -> panels a/c/e (left column, per context)
res_b_CpG <- analyze_one_clean("BREEDING", "CpG", map_file_breeding, colors.17, "BREEDING", "a",
                               prefer_mef = prefer_mef, letter_angle = 45, letter_yfactor = 0.12)
res_b_CHG <- analyze_one_clean("BREEDING", "CHG", map_file_breeding, colors.17, "BREEDING", "c",
                               prefer_mef = prefer_mef, letter_angle = 45, letter_yfactor = 0.12)
res_b_CHH <- analyze_one_clean("BREEDING", "CHH", map_file_breeding, colors.17, "BREEDING", "e",
                               prefer_mef = prefer_mef, letter_angle = 45, letter_yfactor = 0.12)

# NATURAL -> panels b/d/f (right column, per context)
res_n_CpG <- analyze_one_clean("NATURAL", "CpG", map_file_natural, colors.25, "NATURAL", "b",
                               prefer_mef = prefer_mef, letter_angle = 45, letter_yfactor = 0.12, letter_size = 4)
res_n_CHG <- analyze_one_clean("NATURAL", "CHG", map_file_natural, colors.25, "NATURAL", "d",
                               prefer_mef = prefer_mef, letter_angle = 45, letter_yfactor = 0.12, letter_size = 4)
res_n_CHH <- analyze_one_clean("NATURAL", "CHH", map_file_natural, colors.25, "NATURAL", "f",
                               prefer_mef = prefer_mef, letter_angle = 45, letter_yfactor = 0.12, letter_size = 4)

cat("\n\n====================\nSUMMARY\n====================\n")
summ_line <- function(res, cohort, ctx) {
  cat(cohort, " / ", ctx, "\n", sep = "")
  cat("  Method: ", res$tests$method, "\n", sep = "")
  cat("  Shapiro-Wilk p (residuals): ", format(res$tests$shapiro_p, digits = 4), "\n", sep = "")
  cat("  Kruskal-Wallis p: ", format(res$tests$global_p, digits = 8), "\n", sep = "")
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
  file.path(fig_dir, "methylation_levels_panel.tiff"),
  w_cm = 48, h_cm = 42
)

# Summary table of the tests per cohort and context
summ_row <- function(res, cohort, ctx) {
  ph <- res$tests$posthoc
  data.frame(cohort = cohort, context = ctx, n_samples = nrow(res$df),
             n_groups = nlevels(res$df$Group),
             shapiro_p = res$tests$shapiro_p, kruskal_H = res$tests$kw_H,
             kruskal_df = res$tests$kw_df, kruskal_p = res$tests$global_p,
             dunn_sig_pairs = sum(ph$p_adj < 0.05), dunn_total_pairs = nrow(ph),
             mean_methylation = mean(res$df$Methylation))
}
summ <- rbind(summ_row(res_b_CpG, "BREEDING", "CpG"), summ_row(res_b_CHG, "BREEDING", "CHG"),
              summ_row(res_b_CHH, "BREEDING", "CHH"), summ_row(res_n_CpG, "NATURAL", "CpG"),
              summ_row(res_n_CHG, "NATURAL", "CHG"), summ_row(res_n_CHH, "NATURAL", "CHH"))
write.table(summ, file.path(log_dir, "methylation_level_tests_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(summ)

# Values shown in the panel: per-tree mean methylation, with the Dunn's test letters
sd_rows <- function(res, cohort, ctx, panel) {
  d <- res$df
  data.frame(panel = panel, cohort = cohort, context = ctx, tree = d$Sample,
             group = as.character(d$Group), percent_methylation = d$Methylation,
             dunn_letter = unname(res$tests$letters[as.character(d$Group)]),
             kruskal_p = res$tests$global_p)
}
write.table(rbind(sd_rows(res_b_CpG, "BREEDING", "CpG", "a"), sd_rows(res_n_CpG, "NATURAL", "CpG", "b"),
                  sd_rows(res_b_CHG, "BREEDING", "CHG", "c"), sd_rows(res_n_CHG, "NATURAL", "CHG", "d"),
                  sd_rows(res_b_CHH, "BREEDING", "CHH", "e"), sd_rows(res_n_CHH, "NATURAL", "CHH", "f")),
            file.path(fig_dir, "methylation_levels_panel_values.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

message("DONE.\nLog: ", log_file, "\nFigures: ", fig_dir, "\nPosthoc tables: ", log_dir)
