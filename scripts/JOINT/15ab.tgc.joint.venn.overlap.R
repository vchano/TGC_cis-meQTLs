#!/usr/bin/env Rscript
############################################################
# TreeGeneClimate (TGC) — JOINT ECS + TMS
# Step 15ab: Venn-diagram overlap analysis of significant meQTL markers
#
# "Robust" = SNP-site pair significant at p_FDR < 1e-10 in BOTH tools.
#
# OUTPUTS
# -------
# Tool comparison (6 Venns, labelled a)–f)):
#   tools/venn_tools_<cohort>_<ctx>.{tiff,pdf}
#   Order: BREEDING/CpG (a), BREEDING/CHG (b), BREEDING/CHH (c),
#          NATURAL/CpG  (d), NATURAL/CHG  (e), NATURAL/CHH  (f)
#   Sets: GENESIS vs MatrixEQTL — unit: SNPs
#
# Cohort comparison by context (3 Venns, c)–e)):
#   comparisons/venn_cohorts_<ctx>.{tiff,pdf}
#   Sets: BREEDING vs NATURAL — unit: SNP positions (robust)
#
# Context comparison by cohort (2 Venns, a)–b)):
#   comparisons/venn_contexts_<cohort>.{tiff,pdf}
#   Sets: CpG / CHG / CHH — unit: SNPs (robust)
#
# Tables (one per cohort):
#   tables/robust_markers_<cohort>.tsv
#
# Counts shown in the Venn diagrams:
#   panels/venn_panels_values.tsv
#
# USAGE
#   Rscript --vanilla 15ab.R   (or source() in RStudio)
#
# INPUTS
#   RESULTS/JOINT/MEQTL/COMBINED/sig_fdr1e10_<tool>_<cohort>.tsv  (Step 13ab)
#   RESULTS/JOINT/MEQTL/INPUTS/<COHORT>/snp_variant_annot.tsv     (Step 12ab0)
#
# R packages: data.table, ggvenn, ggplot2, patchwork
############################################################

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
TOOLS    <- c("GENESIS", "MATRIXEQTL")

# Per-plot fill colours — each plot within a group gets a distinct colour pair/triple.
# Names match original set identifiers (before display-renaming).

# All three Venn types use SNPs as the unit:
#   Tool comparison   : GENESIS SNPs vs MatrixEQTL SNPs  → overlap = robust SNPs
#   Cohort comparison : BREEDING robust SNPs vs NATURAL robust SNPs (per context)
#   Context comparison: CpG vs CHG vs CHH robust SNPs (per cohort)
# "Robust SNPs" = SNPs significant in both tools for the same cohort × context
# (SNP-level intersection; independent of which methylation site is associated).
# The separate robust MARKER TABLE (section 9) uses the stricter pair-level criterion
# (same SNP–site pair in both tools) as required for downstream annotation.

# Methylation-context colours, shared by all figures (magenta, blue, turquoise;
# no red-green pair), and a light tint of each for the first set of a two-set Venn
CTX_COLORS <- c(CpG = "#D81B60", CHG = "#1F4E9C", CHH = "#17BECF")
CTX_TINTS  <- c(CpG = "#EC8DB0", CHG = "#8FA7CE", CHH = "#8BDFE7")

# Tool comparison (6 plots, GENESIS vs MatrixEQTL): tint / full colour of the context
TOOL_FILL_COLORS <- setNames(
  lapply(rep(CONTEXTS, times = 2), function(ctx) c(GENESIS = CTX_TINTS[[ctx]], MATRIXEQTL = CTX_COLORS[[ctx]])),
  paste0(rep(COHORTS, each = 3), "_", CONTEXTS))

# Cohort comparison (3 plots, BREEDING vs NATURAL): tint / full colour of the context
COHORT_FILL_COLORS <- setNames(
  lapply(CONTEXTS, function(ctx) c(BREEDING = CTX_TINTS[[ctx]], NATURAL = CTX_COLORS[[ctx]])), CONTEXTS)

# Context comparison (2 plots, CpG/CHG/CHH): the context colours
CTX_FILL_COLORS <- list(BREEDING = CTX_COLORS, NATURAL = CTX_COLORS)

# Tool comparison: BREEDING first (a–c), NATURAL second (d–f)
TOOL_LABELS <- setNames(
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
# Set TGC_PROJECT_ROOT as an environment variable, or edit the fallback path below
PROJECT_ROOT <- Sys.getenv("TGC_PROJECT_ROOT", unset = "/path/to/your/project")
# ===========================

SIG_DIR  <- file.path(PROJECT_ROOT, "RESULTS", "JOINT", "MEQTL", "COMBINED")
OUT_ROOT <- file.path(PROJECT_ROOT, "RESULTS", "JOINT", "MEQTL", "OVERLAP")
TOOL_DIR      <- file.path(OUT_ROOT, "tools")
CMP_DIR       <- file.path(OUT_ROOT, "comparisons")
TAB_DIR       <- file.path(OUT_ROOT, "tables")
PANEL_DIR     <- file.path(OUT_ROOT, "panels")
LOG_DIR       <- file.path(PROJECT_ROOT, "RESULTS", "JOINT", "MEQTL", "LOGS")

for (d in c(TOOL_DIR, CMP_DIR, TAB_DIR, PANEL_DIR, LOG_DIR))
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

save_plot <- function(gg, base_path,
                      width_cm = OUT_W, height_cm = OUT_H) {
  if (is.null(gg)) return(invisible(NULL))
  w_in <- width_cm  / 2.54
  h_in <- height_cm / 2.54
  for (fmt in c("tiff", "pdf", "eps", "png")) {
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

# Counts of the Venn regions; a region is the combination of sets
# that contain an element and no other set
venn_counts <- list()
venn_regions <- function(set_list, figure, panel, cohort, context) {
  nms <- names(set_list)
  els <- unique(unlist(set_list))
  memb <- vapply(set_list, function(s) els %in% s, logical(length(els)))
  if (!is.matrix(memb)) memb <- matrix(memb, ncol = length(nms), dimnames = list(NULL, nms))
  region <- apply(memb, 1, function(r) paste(nms[r], collapse = " & "))
  tab <- table(region)
  data.table(figure = figure, panel = sub("[)]$", "", panel), cohort = cohort, context = context,
             sets = paste(nms, collapse = " | "), region = names(tab), n = as.integer(tab))
}

# Human-readable labels for set names
SET_DISPLAY <- c(
  BREEDING    = "Breeding cohort",
  NATURAL     = "Natural cohort",
  GENESIS    = "GENESIS",
  MATRIXEQTL = "MatrixEQTL"
)

# Draw a Venn diagram using ggvenn (ggplot2-native, true circles).
# Returns a ggplot object.
# set_list   : named list of character vectors (original names used for logic)
# fill_colors: named character vector of fill colours, keyed by original names
# panel_label: panel letter string, e.g. "a)"
draw_venn <- function(set_list, fill_colors, panel_label = NULL, sname_size = 7, sd = NULL) {

  set_list <- set_list[lengths(set_list) > 0]
  if (!is.null(sd) && length(set_list) >= 2)
    venn_counts[[length(venn_counts) + 1L]] <<- venn_regions(set_list, sd$figure, panel_label, sd$cohort, sd$context)
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
    lbl_up <- sub("[)]$", "", panel_label)
    p <- p + labs(tag = lbl_up) +
      theme(plot.tag = element_text(size = 22, face = "bold", hjust = 0))
  }

  p
}

############################################################
# 3b) OUTPUT DIMENSIONS (updated)
############################################################

# Tool comparison (6 Venns, a)–f)):
#   tools/venn_tools_<cohort>_<ctx>.{tiff,pdf}
#   Sets: GENESIS vs MatrixEQTL — unit: SNPs
#
# Cohort comparison by context (3 Venns, c)–e)):
#   comparisons/venn_cohorts_<ctx>.{tiff,pdf}
#   Sets: BREEDING vs NATURAL — unit: SNPs (robust)
#
# Context comparison by cohort (2 Venns, a)–b)):
#   comparisons/venn_contexts_<cohort>.{tiff,pdf}
#   Sets: CpG / CHG / CHH — unit: SNPs (robust)

############################################################
# 4) LOAD SIGNIFICANT DATA
############################################################

log_msg("Step 15ab — Venn overlap analysis")
log_msg("ggvenn version: ", as.character(packageVersion("ggvenn")))
log_msg("Reading FDR < 1e-10 tables from: ", SIG_DIR)

sig_raw <- list()
for (tool in TOOLS) {
  sig_raw[[tool]] <- list()
  for (cohort in COHORTS) {
    fname <- file.path(SIG_DIR,
      sprintf("sig_fdr1e10_%s_%s.tsv", tolower(tool), tolower(cohort)))
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
#     Positions of SNP identifiers from the cohort variant annotation (Step 12ab0);
#     cross-cohort comparisons use chr:pos, never integer SNP identifiers.
############################################################

ANNOT_ROOT <- file.path(PROJECT_ROOT, "RESULTS", "JOINT", "MEQTL", "INPUTS")

va_map <- list()   # va_map[[cohort]][[ctx]]: named char vector  pos_string[snp_id]
for (cohort in COHORTS) {
  va_map[[cohort]] <- list()
  for (ctx in CONTEXTS) {
    va_path <- file.path(ANNOT_ROOT, cohort, "snp_variant_annot.tsv")
    if (file.exists(va_path)) {
      va <- fread(va_path)
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
    g5  <- sig_raw[["GENESIS"   ]][[cohort]]
    me5 <- sig_raw[["MATRIXEQTL"]][[cohort]]

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

    g5_keys  <- g5_sub [, .(snp, site, FDR_GENESIS    = p_FDR)]
    me5_keys <- me5_sub[, .(snp, site, FDR_MATRIXEQTL = p_FDR)]

    both <- merge(g5_keys, me5_keys, by = c("snp", "site"))
    both[, context := ctx]

    # Attach position columns from GENESIS if available
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
# 6) TOOL COMPARISON — 6 Venns: GENESIS vs MatrixEQTL per cohort × context
#    Panel labels a)–f): BREEDING (a–c), NATURAL (d–f)
#    Unit: unique genomic positions (SNP chr:pos)
############################################################

# Storage for panel assembly (populated in sections 6-8)
tool_plots   <- list()   # 6 tool-comparison Venns  (key = "COHORT_ctx")
cohort_plots <- list()   # 3 cohort-comparison Venns (key = context)
ctx_plots    <- list()   # 2 context-comparison Venns (key = cohort)

log_msg("--- Tool comparison (6 Venns) ---")

get_sites <- function(dt, ctx) {
  if (!nrow(dt) || !"context" %in% names(dt)) return(character(0))
  unique(dt[context == ctx, site])
}

for (cohort in COHORTS) {
  for (ctx in CONTEXTS) {
    key <- paste0(cohort, "_", ctx)
    lbl <- TOOL_LABELS[[key]]
    log_msg("  Panel ", lbl, "  [", cohort, " / ", ctx, "]")

    # Tool-comparison Venns: pair-level SNP counts.
    # The intersection shown = SNPs targeting the SAME methylation site in both tools
    # (robust pair-level), NOT just the SNP-level set intersection.
    # Synthetic integer sets are built so ggvenn computes the correct counts:
    #   set_A (GENESIS)    : |A| = n_g5,  |A ∩ B| = n_rob
    #   set_B (MatrixEQTL) : |B| = n_me5, |A ∩ B| = n_rob
    # robust IDs = 1..n_rob; G5-only = (n_rob+1)..(n_rob+n_g5_only)
    # ME5-only starts after G5-only to guarantee G5-only ∩ ME5-only = ∅
    g5_sub_ctx  <- sig_raw[["GENESIS"   ]][[cohort]]
    me5_sub_ctx <- sig_raw[["MATRIXEQTL"]][[cohort]]
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
    if (n_g5  > 0L) set_list[["GENESIS"]]    <- set_a
    if (n_me5 > 0L) set_list[["MATRIXEQTL"]] <- set_b

    log_msg("  Pair-level counts — G5: ", n_g5, " | ME5: ", n_me5,
            " | Robust: ", n_rob)

    base      <- file.path(TOOL_DIR,
      paste0("venn_tools_", tolower(cohort), "_", tolower(ctx)))

    gg <- draw_venn(set_list, TOOL_FILL_COLORS[[key]], lbl, sname_size = 5.5,
                    sd = list(figure = "tools", cohort = cohort, context = ctx))
    tool_plots[[key]] <- gg
    save_plot(gg, base,      width_cm = 18, height_cm = 12)
  }
}

############################################################
# 7) COHORT COMPARISON — 3 Venns: BREEDING vs NATURAL per context (c)–e))
#    Unit: unique genomic positions (SNP chr:pos)
############################################################

log_msg("--- Cohort comparison by context (3 Venns) ---")

for (ctx in CONTEXTS) {
  lbl <- CTX_PANEL_LABELS[[ctx]]
  log_msg("  Panel ", lbl, "  [", ctx, "]")

  breed_snps   <- robust_snp_pos[["BREEDING"]][[ctx]]
  natural_snps <- robust_snp_pos[["NATURAL" ]][[ctx]]

  set_list <- list()
  if (length(breed_snps))   set_list[["BREEDING"]] <- breed_snps
  if (length(natural_snps)) set_list[["NATURAL"]]  <- natural_snps

  base      <- file.path(CMP_DIR,      paste0("venn_cohorts_", tolower(ctx)))

  gg <- draw_venn(set_list, COHORT_FILL_COLORS[[ctx]], lbl, sname_size = 5.5,
                  sd = list(figure = "cohorts_by_context", cohort = "BREEDING | NATURAL", context = ctx))
  cohort_plots[[ctx]] <- gg
  save_plot(gg, base,      width_cm = 14, height_cm = 10)
}

############################################################
# 8) CONTEXT COMPARISON — 2 Venns: CpG vs CHG vs CHH per cohort (a)–b))
#    Unit: unique genomic positions (SNP chr:pos)
############################################################

log_msg("--- Context comparison by cohort (2 Venns) ---")

for (cohort in COHORTS) {
  lbl <- COHORT_PANEL_LABELS[[cohort]]
  log_msg("  Panel ", lbl, "  [", cohort, "]")

  snp_sets <- list()
  for (ctx in CONTEXTS) {
    snps <- robust_snp_pos[[cohort]][[ctx]]
    if (length(snps)) snp_sets[[ctx]] <- snps
  }

  base      <- file.path(CMP_DIR,      paste0("venn_contexts_", tolower(cohort)))

  gg <- draw_venn(snp_sets, CTX_FILL_COLORS[[cohort]], lbl,
                  sd = list(figure = "contexts_by_cohort", cohort = cohort, context = "CpG | CHG | CHH"))
  ctx_plots[[cohort]] <- gg
  save_plot(gg, base,      width_cm = 21, height_cm = 14)
}

############################################################
# 8b) PANEL ASSEMBLY
#   Comparison panel: row1 = ctx_plots (a, b — context overlap per cohort)
#                     row2 = cohort_plots (c-e — cross-cohort position overlap per context)
#   Tool panel      : 2×3 grid, GENESIS vs MatrixEQTL per cohort×context
############################################################

log_msg("--- Panel assembly ---")

row1_valid <- Filter(Negate(is.null), ctx_plots)
row2_valid <- Filter(Negate(is.null), cohort_plots)

if (length(row1_valid) >= 1 && length(row2_valid) >= 1) {
  row1 <- Reduce(`+`, row1_valid) + plot_layout(ncol = length(row1_valid))
  row2 <- Reduce(`+`, row2_valid) + plot_layout(ncol = length(row2_valid))
  panel_plot <- row1 / row2 + plot_layout(heights = c(1, 1))
  save_plot(panel_plot,
            file.path(PANEL_DIR, "venn_contexts_cohorts_panel"),
            width_cm = 42, height_cm = 28)
  log_msg("  comparison panel saved → ", PANEL_DIR)
} else {
  log_msg("  WARNING: insufficient Venn plots for comparison panel")
}

row1s <- Filter(Negate(is.null), tool_plots[paste0("BREEDING_", CONTEXTS)])
row2s <- Filter(Negate(is.null), tool_plots[paste0("NATURAL_",  CONTEXTS)])

if (length(row1s) >= 1 && length(row2s) >= 1) {
  row1 <- Reduce(`+`, row1s) + plot_layout(ncol = 3)
  row2 <- Reduce(`+`, row2s) + plot_layout(ncol = 3)
  fs4  <- row1 / row2 + plot_layout(heights = c(1, 1))
  save_plot(fs4,
            file.path(PANEL_DIR, "venn_tools_panel"),
            width_cm = 54, height_cm = 24)
  log_msg("  tool panel saved → ", PANEL_DIR)
} else {
  log_msg("  WARNING: insufficient Venn plots for tool panel")
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
      "context", "FDR_GENESIS", "FDR_MATRIXEQTL"),
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

# Counts of every Venn region
fwrite(rbindlist(venn_counts), file.path(PANEL_DIR, "venn_panels_values.tsv"), sep = "\t")

############################################################
# 10) SESSION INFO
############################################################

writeLines(capture.output(sessionInfo()),
           file.path(LOG_DIR, "step15ab_sessionInfo.txt"))

log_msg("Step 15ab finished — outputs in: ", OUT_ROOT)
