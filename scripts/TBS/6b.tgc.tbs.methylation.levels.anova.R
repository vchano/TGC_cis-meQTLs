#!/usr/bin/env Rscript
############################################################
# TreeGeneClimate (TGC) — TBS
# Step 6b (UPDATED): Plot methylation level (%) per group
#   - Breeding cohort: per family (17 families) -> panels a/b/c
#   - Natural cohort:  per stand  (25 stands)   -> panels d/e/f
#   - One plot per context (CpG/CHG/CHH)
#   - Statistical tests + posthoc + logs
#
# STYLE:
#   - Titles: only "a)" ... "f)" (no extra title text)
#   - No legend
#   - Annotation: only p-value (no method label)
#   - Posthoc letters angled; spaced above boxes to reduce overlap
#
# PANEL LAYOUT (FIG34_PANEL_COMBINED):
#   Top row:    a) breeding CpG  |  b) natural CpG
#   Middle row: c) breeding CHG  |  d) natural CHG
#   Bottom row: e) breeding CHH  |  f) natural CHH
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
PROJECT_ROOT <- "/path/to/your/project"  # <-- set this
# ===========================

rds_dir  <- file.path(PROJECT_ROOT, "RESULTS/TBS/RANALYSIS/METHYLKIT_OBJECTS")
fig_dir  <- file.path(PROJECT_ROOT, "RESULTS/TBS/RANALYSIS/FIGURES/FIG3_FIG4")
log_dir  <- file.path(PROJECT_ROOT, "RESULTS/TBS/RANALYSIS/ANOVA.METHYL.LEVEL")

# Tab-delimited files mapping sample IDs to family (breeding) or population (natural)
map_file_breeding <- file.path(PROJECT_ROOT, "DATA/METADATA/breeding_sample2family.txt")
map_file_natural  <- file.path(PROJECT_ROOT, "DATA/METADATA/natural_sample2pop.txt")

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)

log_file <- file.path(log_dir, "Step6b_methylation_level_stats_posthoc.log")

## ---------------------------
## Color codes (yours)
## ---------------------------
# Named color vectors ensure consistent group-to-color mapping across all panels
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
# The unit of observation for the ANOVA is one value per sample (the genome-wide
# mean % methylation), which satisfies independence assumptions better than
# treating individual sites as replicates.
methylbase_to_sample_means <- function(mb) {
  d <- getData(mb)

  # Identify paired count columns by the methylKit naming convention numCsN / numTsN
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

  # Compute % methylation per site per sample; NaN arises when coverage = 0
  perc_mat <- vapply(seq_len(n), function(i) {
    numCs <- d[[numCs_cols[i]]]
    numTs <- d[[numTs_cols[i]]]
    perc <- 100 * (numCs / (numCs + numTs))
    perc[is.nan(perc)] <- NA_real_
    perc
  }, numeric(nrow(d)))

  colnames(perc_mat) <- sample_ids
  colMeans(perc_mat, na.rm = TRUE)   # one mean value per sample
}

# Global test: always ANOVA + TukeyHSD (per-sample means are CLT-justified).
# Shapiro-Wilk and Levene tests are run and reported for transparency but do
# not affect the choice of test.
run_group_tests <- function(df) {
  df <- df[!is.na(df$Group) & !is.na(df$Methylation), ]
  df$Group <- droplevels(df$Group)

  aov_fit <- aov(Methylation ~ Group, data = df)

  # Normality and homoscedasticity checks reported in log but not used as gatekeepers
  sh_p  <- tryCatch(shapiro.test(residuals(aov_fit))$p.value, error = function(e) NA_real_)
  lev_p <- tryCatch(car::leveneTest(Methylation ~ Group, data = df)[["Pr(>F)"]][1], error = function(e) NA_real_)

  an_p    <- summary(aov_fit)[[1]][["Pr(>F)"]][1]
  tuk     <- TukeyHSD(aov_fit)
  # Convert pairwise Tukey contrasts to compact letter display (CLD) for plotting
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

  # Restrict palette to groups actually present in this cohort/context
  present <- levels(df$Group)
  color_vec <- color_vec[names(color_vec) %in% present]

  ggplot(df, aes(x = Group, y = Methylation, fill = Group)) +
    geom_boxplot(outlier.shape = NA, width = 0.7) +   # outliers suppressed; shown via jitter
    geom_jitter(width = 0.2, size = 1, color = "black", alpha = 0.7) +
    scale_fill_manual(values = color_vec, drop = FALSE) +
    labs(title = title_letter, x = NULL, y = ylab) +
    theme_minimal(base_size = 18) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 13),
      axis.title.y = element_text(size = 16),
      plot.title  = element_text(size = 22, hjust = 0, vjust = 1),
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

  # Align group labels to sample order in the methylBase object
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

  # Place global ANOVA p-value above the first group's box
  p <- p + annotate(
    "text",
    x = 1,
    y = max(df$Methylation, na.rm = TRUE) + 0.22 * y_span,
    label = fmt_p(tests$global_p),
    hjust = 0, vjust = 1.2, size = 4, fontface = "italic"
  )

  # Posthoc letters angled + spaced; letters_vec must be ordered to match group factor levels
  letters_vec <- tests$letters
  letters_vec <- letters_vec[levels(df$Group)]

  # Compute per-group maximum to position CLD letters above the boxes
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
res_b_CpG <- analyze_one_clean("BREEDING", "CpG", map_file_breeding, colors.17, "FIG3a_BREEDING", "a)",
                               prefer_mef = prefer_mef, letter_angle = 45, letter_yfactor = 0.12)
res_b_CHG <- analyze_one_clean("BREEDING", "CHG", map_file_breeding, colors.17, "FIG3c_BREEDING", "c)",
                               prefer_mef = prefer_mef, letter_angle = 45, letter_yfactor = 0.12)
res_b_CHH <- analyze_one_clean("BREEDING", "CHH", map_file_breeding, colors.17, "FIG3e_BREEDING", "e)",
                               prefer_mef = prefer_mef, letter_angle = 45, letter_yfactor = 0.12)

# NATURAL -> panels b/d/f (right column, per context)
res_n_CpG <- analyze_one_clean("NATURAL", "CpG", map_file_natural, colors.25, "FIG4b_NATURAL", "b)",
                               prefer_mef = prefer_mef, letter_angle = 45, letter_yfactor = 0.12)
res_n_CHG <- analyze_one_clean("NATURAL", "CHG", map_file_natural, colors.25, "FIG4d_NATURAL", "d)",
                               prefer_mef = prefer_mef, letter_angle = 45, letter_yfactor = 0.12)
res_n_CHH <- analyze_one_clean("NATURAL", "CHH", map_file_natural, colors.25, "FIG4f_NATURAL", "f)",
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

# Assemble 2x3 combined panel: rows = context (CpG/CHG/CHH), columns = cohort (breeding | natural)
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

message("DONE.\nLog: ", log_file, "\nFigures: ", fig_dir, "\nPosthoc tables: ", log_dir)
sessionInfo()
