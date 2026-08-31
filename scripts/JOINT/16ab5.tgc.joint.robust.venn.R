#!/usr/bin/env Rscript
############################################################
# TreeGeneClimate (TGC) — JOINT ECS + TBS
# Step 16ab5: Venn-diagram overlap analysis of significant meQTL markers
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
#   Rscript --vanilla 16ab5.R   (or source() in RStudio)
#
# INPUTS
#   RESULTS/JOINT/COMBINED5/sig_sites/sig_p1e10_*.tsv  (from 16ab3.R)
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

PROJECT_ROOT <- "/mnt/ceph-hdd/projects/scc_ufff_gailing/chano_TGC"

SIG_DIR  <- file.path(PROJECT_ROOT, "RESULTS", "JOINT", "COMBINED5", "sig_sites")
OUT_ROOT <- file.path(PROJECT_ROOT, "RESULTS", "JOINT", "COMBINED5", "overlap")
SUPP_DIR      <- file.path(OUT_ROOT, "supp")
MAIN_DIR      <- file.path(OUT_ROOT, "main")
TAB_DIR       <- file.path(OUT_ROOT, "tables")
PANEL_DIR     <- file.path(OUT_ROOT, "panels")
LOG_DIR       <- file.path(PROJECT_ROOT, "RESULTS", "JOINT", "COMBINED5", "LOGS")
CORRECTED_DIR <- file.path(PROJECT_ROOT, "RESULTS", "CORRECTED", "FIGURES", "NEW")
REV1_DIR      <- file.path(PROJECT_ROOT, "RESULTS", "DRAFT", "NATURE.GENETICS", "REV1")

for (d in c(SUPP_DIR, MAIN_DIR, TAB_DIR, PANEL_DIR, LOG_DIR, CORRECTED_DIR, REV1_DIR))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

LOGFILE <- file.path(LOG_DIR, "step16ab5.log")
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

log_msg("Step 16ab5 — Venn overlap analysis")
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
      paste0("venn_tools_", tolower(cohort), "_", tolower(ctx), "_fixed"))
    base_corr <- file.path(CORRECTED_DIR,
      paste0("venn_tools_", tolower(cohort), "_", tolower(ctx), "_fixed"))

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

  base      <- file.path(MAIN_DIR,      paste0("venn_cohorts_", tolower(ctx), "_fixed"))
  base_corr <- file.path(CORRECTED_DIR, paste0("venn_cohorts_", tolower(ctx), "_fixed"))

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

  base      <- file.path(MAIN_DIR,      paste0("venn_contexts_", tolower(cohort), "_fixed"))
  base_corr <- file.path(CORRECTED_DIR, paste0("venn_contexts_", tolower(cohort), "_fixed"))

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
            file.path(PANEL_DIR, "Figure7_Venn_contexts_cohorts_panel_fixed"),
            width_cm = 42, height_cm = 28)
  save_plot(fig7,
            file.path(CORRECTED_DIR, "Figure6_Venn_contexts_cohorts_panel_fixed"),
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
            file.path(PANEL_DIR, "SuppFig5_Venn_tools_panel_fixed"),
            width_cm = 54, height_cm = 24)
  save_plot(fs5,
            file.path(CORRECTED_DIR, "SuppFig5_Venn_tools_panel_fixed"),
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
           file.path(LOG_DIR, "step16ab5_sessionInfo.txt"))

log_msg("Step 16ab5 finished — outputs in: ", OUT_ROOT)
