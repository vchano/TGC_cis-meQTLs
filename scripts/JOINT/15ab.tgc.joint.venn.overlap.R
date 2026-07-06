#!/usr/bin/env Rscript
############################################################
# TreeGeneClimate (TGC) — JOINT ECS + TBS
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

# Panel labels follow figure conventions: supplementary a)–f), main A c)–e), main B a)–b)
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
PROJECT_ROOT <- "/path/to/your/project"  # <-- set this
# ===========================

# Input: significant pairs at the strict FDR threshold, produced by 13ab.R
SIG_DIR  <- file.path(PROJECT_ROOT, "RESULTS", "JOINT", "COMBINED5", "sig_sites")
OUT_ROOT <- file.path(PROJECT_ROOT, "RESULTS", "JOINT", "COMBINED5", "overlap")
SUPP_DIR  <- file.path(OUT_ROOT, "supp")
MAIN_DIR  <- file.path(OUT_ROOT, "main")
TAB_DIR   <- file.path(OUT_ROOT, "tables")
PANEL_DIR <- file.path(OUT_ROOT, "panels")
LOG_DIR   <- file.path(PROJECT_ROOT, "RESULTS", "JOINT", "COMBINED5", "LOGS")

for (d in c(SUPP_DIR, MAIN_DIR, TAB_DIR, PANEL_DIR, LOG_DIR))
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

# Save a ggplot in all four formats; eps via cairo_ps for vector-safe rendering
save_plot <- function(gg, base_path,
                      width_cm = OUT_W, height_cm = OUT_H) {
  if (is.null(gg)) return(invisible(NULL))
  for (fmt in c("tiff", "pdf", "eps", "png")) {
    out    <- paste0(base_path, ".", fmt)
    dev    <- if (fmt == "eps") cairo_ps else fmt
    dpi    <- if (fmt %in% c("tiff", "png")) OUT_RES else 150
    tryCatch(
      ggsave(out, plot = gg,
             width = width_cm, height = height_cm, units = "cm",
             device = dev, dpi = dpi),
      error = function(e) log_msg("  ERROR (", fmt, "): ", conditionMessage(e))
    )
    log_msg("  Saved: ", out)
  }
}

# Human-readable labels for set names shown inside the Venn diagrams
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

  # Return an informative placeholder when fewer than 2 sets have any markers
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

  # Apply human-readable display names after the set-ordering logic above
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

  # Crop horizontal whitespace for 2-set Venns (circles at ±0.85, r=1.5,
  # so data spans ±2.35; ggplot default expansion pushes well beyond that)
  if (n == 2)
    p <- p + coord_fixed(ratio = 1, clip = "off")

  if (!is.null(panel_label))
    p <- p + labs(tag = panel_label) +
      theme(plot.tag = element_text(size = 22, face = "plain", hjust = 0))

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

# Read FDR-strict significant pairs for both tools and cohorts
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
    # Ensure snp and site IDs are character to prevent integer-vs-character merge issues
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
# 5) BUILD ROBUST SETS  (intersection of both tools, per cohort × context)
############################################################

log_msg("Building robust sets (intersection of both tools)...")

robust      <- list()   # pair-level robust (same snp+site in both tools) — used for tables
robust_snps <- list()   # SNP-level robust  (same SNP in both tools)     — used for Venns

for (cohort in COHORTS) {
  robust[[cohort]]      <- list()
  robust_snps[[cohort]] <- list()
  for (ctx in CONTEXTS) {
    g5  <- sig_raw[["GENESIS5"   ]][[cohort]]
    me5 <- sig_raw[["MATRIXEQTL5"]][[cohort]]

    empty <- data.table()
    # Skip if either tool produced no results or lacks required columns
    if (!nrow(g5) || !nrow(me5) ||
        !"context" %in% names(g5) || !"context" %in% names(me5) ||
        !"p_FDR"   %in% names(g5) || !"p_FDR"   %in% names(me5)) {
      robust[[cohort]][[ctx]]      <- empty
      robust_snps[[cohort]][[ctx]] <- character(0)
      next
    }

    g5_sub  <- g5 [context == ctx]
    me5_sub <- me5[context == ctx]

    if (!nrow(g5_sub) || !nrow(me5_sub)) {
      log_msg("  ", cohort, "/", ctx, ": one or both tools have 0 pairs")
      robust[[cohort]][[ctx]]      <- empty
      robust_snps[[cohort]][[ctx]] <- character(0)
      next
    }

    g5_keys  <- g5_sub [, .(snp, site, FDR_GENESIS5    = p_FDR)]
    me5_keys <- me5_sub[, .(snp, site, FDR_MATRIXEQTL5 = p_FDR)]

    # Inner join on snp+site: keeps only pairs significant in both tools
    both <- merge(g5_keys, me5_keys, by = c("snp", "site"))
    both[, context := ctx]

    # Attach position columns from GENESIS5 if available (used in downstream tables)
    pos_cols <- intersect(c("snp_chr", "snp_pos", "site_chr", "site_pos"),
                          names(g5_sub))
    if (length(pos_cols)) {
      g5_pos <- g5_sub[, c("snp", "site", pos_cols), with = FALSE]
      both   <- merge(both, g5_pos, by = c("snp", "site"), all.x = TRUE)
    }

    robust[[cohort]][[ctx]] <- both
    log_msg("  ", cohort, "/", ctx, ": ", nrow(both), " robust pairs | ",
            uniqueN(both$site), " sites | ", uniqueN(both$snp), " SNPs (pair-level)")

    # SNP-level robust: SNPs significant in both tools (regardless of which site)
    robust_snps[[cohort]][[ctx]] <- intersect(unique(g5_sub$snp), unique(me5_sub$snp))
    log_msg("  ", cohort, "/", ctx, ": ",
            length(robust_snps[[cohort]][[ctx]]), " SNPs (SNP-level robust)")
  }
}

############################################################
# 6) SUPPLEMENTARY — 6 Venns: GENESIS5 vs MatrixEQTL5 per cohort × context
#    Panel labels a)–f): BREEDING (a–c), NATURAL (d–f)
#    Unit: methylation sites
############################################################

# Storage for panel assembly (populated in sections 6-8)
supp_plots   <- list()   # 6 tool-comparison Venns  (key = "COHORT_ctx")
cohort_plots <- list()   # 3 cohort-comparison Venns (key = context)
ctx_plots    <- list()   # 2 context-comparison Venns (key = cohort)

log_msg("--- Supplementary: tool comparison (6 Venns) ---")

# Helper to extract unique sites for a given context from a results table
get_sites <- function(dt, ctx) {
  if (!nrow(dt) || !"context" %in% names(dt)) return(character(0))
  unique(dt[context == ctx, site])
}

# Helper to extract unique SNPs for a given context from a results table
get_snps <- function(dt, ctx) {
  if (!nrow(dt) || !"context" %in% names(dt)) return(character(0))
  unique(dt[context == ctx, snp])
}

for (cohort in COHORTS) {
  for (ctx in CONTEXTS) {
    key <- paste0(cohort, "_", ctx)
    lbl <- SUPP_LABELS[[key]]
    log_msg("  Panel ", lbl, "  [", cohort, " / ", ctx, "]")

    # Supplementary Venns compare tool agreement — unit: SNPs with significant cis-meQTL
    g5_snps  <- get_snps(sig_raw[["GENESIS5"   ]][[cohort]], ctx)
    me5_snps <- get_snps(sig_raw[["MATRIXEQTL5"]][[cohort]], ctx)

    set_list <- list()
    if (length(g5_snps))  set_list[["GENESIS5"]]    <- g5_snps
    if (length(me5_snps)) set_list[["MATRIXEQTL5"]] <- me5_snps

    base <- file.path(SUPP_DIR,
      paste0("venn_tools_", tolower(cohort), "_", tolower(ctx)))

    gg <- draw_venn(set_list, SUPP_FILL_COLORS[[key]], lbl, sname_size = 5.5)
    supp_plots[[key]] <- gg
    save_plot(gg, base, width_cm = 18, height_cm = 12)
  }
}

############################################################
# 7) MAIN A — 3 Venns: BREEDING vs NATURAL per context (a)–c))
#    Unit: methylation sites (robust markers only)
############################################################

log_msg("--- Main A: cohort comparison by context (3 Venns) ---")

for (ctx in CONTEXTS) {
  lbl <- CTX_PANEL_LABELS[[ctx]]
  log_msg("  Panel ", lbl, "  [", ctx, "]")

  # Robust SNPs per cohort: markers detected by both tools in the same cohort
  breed_snps   <- robust_snps[["BREEDING"]][[ctx]]
  natural_snps <- robust_snps[["NATURAL" ]][[ctx]]

  set_list <- list()
  if (length(breed_snps))   set_list[["BREEDING"]] <- breed_snps
  if (length(natural_snps)) set_list[["NATURAL"]]  <- natural_snps

  base <- file.path(MAIN_DIR, paste0("venn_cohorts_", tolower(ctx)))

  gg <- draw_venn(set_list, COHORT_FILL_COLORS[[ctx]], lbl, sname_size = 5.5)
  cohort_plots[[ctx]] <- gg
  save_plot(gg, base, width_cm = 14, height_cm = 10)
}

############################################################
# 8) MAIN B — 2 Venns: CpG vs CHG vs CHH per cohort (a)–b))
#    Unit: SNPs (robust markers only)
############################################################

log_msg("--- Main B: context comparison by cohort (2 Venns) ---")

for (cohort in COHORTS) {
  lbl <- COHORT_PANEL_LABELS[[cohort]]
  log_msg("  Panel ", lbl, "  [", cohort, "]")

  # Compare robust SNP sets across the three cytosine contexts within each cohort
  snp_sets <- list()
  for (ctx in CONTEXTS) {
    snps <- robust_snps[[cohort]][[ctx]]
    if (length(snps)) snp_sets[[ctx]] <- snps
  }

  base <- file.path(MAIN_DIR, paste0("venn_contexts_", tolower(cohort)))

  gg <- draw_venn(snp_sets, CTX_FILL_COLORS[[cohort]], lbl)
  ctx_plots[[cohort]] <- gg
  save_plot(gg, base, width_cm = 21, height_cm = 14)
}

############################################################
# 8b) PANEL ASSEMBLY — removed; all 11 Venns saved individually above
############################################################

############################################################
# 9) TABLES — one per cohort, robust pairs, all contexts combined
############################################################

log_msg("--- Saving robust marker tables ---")

for (cohort in COHORTS) {
  # Combine all three context-level robust pair tables for this cohort
  rows <- lapply(CONTEXTS, function(ctx) robust[[cohort]][[ctx]])
  rows <- rows[sapply(rows, nrow) > 0]

  if (!length(rows)) {
    log_msg("  No robust pairs for ", cohort, " — no table written")
    next
  }

  tab <- rbindlist(rows, fill = TRUE)
  # Standardise column order: identifiers, then positions, then FDR values
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

# Distinguish pair-level SNPs (same snp+site in both tools) from SNP-level
# (same SNP in either tool, regardless of associated site) for reporting
rob_ctx_rows <- list()
for (cohort in COHORTS) {
  for (ctx in CONTEXTS) {
    rob <- robust[[cohort]][[ctx]]
    n_pairs      <- if (is.data.table(rob) && nrow(rob) > 0) nrow(rob)          else 0L
    n_sites      <- if (is.data.table(rob) && nrow(rob) > 0) uniqueN(rob$site)  else 0L
    n_snps_pair  <- if (is.data.table(rob) && nrow(rob) > 0) uniqueN(rob$snp)   else 0L
    n_snps_snp   <- length(robust_snps[[cohort]][[ctx]])
    rob_ctx_rows[[paste(cohort, ctx)]] <- data.table(
      cohort                    = cohort,
      context                   = ctx,
      robust_pairs              = n_pairs,
      robust_unique_sites       = n_sites,
      robust_unique_snps_pairlevel = n_snps_pair,
      robust_unique_snps_snplevel  = n_snps_snp
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
