#!/usr/bin/env Rscript
############################################################
# TreeGeneClimate (TGC) — JOINT ECS + TBS
# Step 19ab: Combined circular Manhattan panels
#            Version 5 — 18 cm canvas, v5 rendering specs
#
# Produces 4 combined panels (each = BREEDING A + NATURAL B + legend):
#   Fig. 3      : GENESIS5,    raw p-value axis   (_5)
#   EDF4        : MatrixEQTL5, raw p-value axis   (_5)
#   Fig. 3 FDR  : GENESIS5,    BH-FDR axis        (_5fdr)
#   EDF4 FDR    : MatrixEQTL5, BH-FDR axis        (_5fdr)
#
# Each panel saved as: tiff, pdf, svg, png, eps
#   tiff -- assembled from individual per-cohort TIFFs via ImageMagick
#   pdf  -- drawn directly via cairo_pdf (36 cm x 20 cm)
#   svg  -- drawn directly via svg()
#   png  -- derived from tiff via ImageMagick (150 dpi)
#   eps  -- derived from tiff via ImageMagick (300 dpi)
#
# OUTPUTS
#   RESULTS/JOINT/COMBINED5/panels/
#     Figure3_circular_panel_genesis5_5.{tiff,pdf,svg,png,eps}
#     EDF4_circular_panel_matrixeqtl5_5.{tiff,pdf,svg,png,eps}
#     Figure3_circular_panel_genesis5_5_FDR.{tiff,pdf,svg,png,eps}
#     EDF4_circular_panel_matrixeqtl5_5_FDR.{tiff,pdf,svg,png,eps}
#   RESULTS/CORRECTED/FIGURES/NEW/  (copies of all above)
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(SNPRelate)
  library(gdsfmt)
  library(circlize)
  library(grid)
  library(RColorBrewer)
})

options(stringsAsFactors = FALSE)

############################################################
# 1) PATHS + SOURCE 16ab.tgc.joint.panel.assembly.R
#    (functions + settings; guard skips main rendering loop)
############################################################

PROJECT_ROOT  <- "/path/to/your/project"   # <-- set this (one place for both scripts)
SOURCED_16AB4 <- TRUE
FDR_AXIS      <- FALSE   # default; overridden per panel below
source(file.path(PROJECT_ROOT, "SCRIPTS", "JOINT",
                 "16ab.tgc.joint.panel.assembly.R"))

############################################################
# 2) PATHS (continued)
############################################################

# PROJECT_ROOT already set above
MAN_DIR      <- file.path(PROJECT_ROOT, "RESULTS", "JOINT", "COMBINED5", "manhattan")
PANEL_DIR    <- file.path(PROJECT_ROOT, "RESULTS", "JOINT", "COMBINED5", "panels")
NEW_DIR      <- file.path(PROJECT_ROOT, "RESULTS", "CORRECTED", "FIGURES", "NEW")

dir.create(PANEL_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(NEW_DIR,   recursive = TRUE, showWarnings = FALSE)

LOG_DIR <- file.path(PROJECT_ROOT, "RESULTS", "JOINT", "COMBINED5", "LOGS")
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)
LOGFILE <- file.path(LOG_DIR, "step19ab_combined.log")
if (file.exists(LOGFILE)) file.remove(LOGFILE)

log_msg <- function(...) {
  txt <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ",
                paste0(..., collapse = ""))
  cat(txt, "\n")
  cat(txt, "\n", file = LOGFILE, append = TRUE)
}

############################################################
# 3) CHROMOSOME LAYOUT
############################################################

log_msg("Building chromosome layout...")
contig_len_dt <- build_contig_lengths(unlist(GDS_FILES))
chrom_map_df  <- read_chrom_map(CHROM_MAP_FILE)
if (is.null(chrom_map_df)) {
  chrom_map_df <- build_default_chrom_map(contig_len_dt)
}
layout_info    <- prepare_layout(contig_len_dt, chrom_map_df)
contig_info_dt <- layout_info$contig_info
sector_lens    <- layout_info$sector_lens
log_msg("Sectors: ", paste(names(sector_lens), collapse = ", "))

############################################################
# 4) COMBINED DEVICE DIMENSIONS
############################################################

PNL_W_CM  <- OUT_W_CM          # 18 cm (from 16ab script)
PNL_H_CM  <- OUT_H_CM          # 18 cm
LEG_H_CM  <- 2.0               # reduced from 3.5 in v4
COMB_W_CM <- PNL_W_CM * 2      # 36 cm
COMB_H_CM <- PNL_H_CM + LEG_H_CM  # 20 cm

ctx_colors_global <- setNames(
  brewer.pal(max(3L, length(CTX_ORDER)), POINT_PALETTE)[seq_along(CTX_ORDER)],
  CTX_ORDER
)

MAGICK_BIN <- Sys.which("magick")
if (!nzchar(MAGICK_BIN)) MAGICK_BIN <- Sys.which("convert")
HAS_MAGICK <- nzchar(MAGICK_BIN)
if (HAS_MAGICK) log_msg("ImageMagick: ", MAGICK_BIN) else
  log_msg("WARNING: ImageMagick not found -- raster formats skipped")

run_magick <- function(args_str) {
  cmd <- paste(shQuote(MAGICK_BIN), args_str)
  ret <- system(cmd, ignore.stdout = TRUE, ignore.stderr = FALSE)
  if (ret != 0L) stop("magick failed:\n  ", cmd)
}

check_nonempty <- function(path, min_bytes = 5000L) {
  file.exists(path) && file.info(path)$size >= min_bytes
}

############################################################
# 5) DRAW-COMBINED FUNCTION
############################################################

draw_combined_panel <- function(tool, out_base, fdr_mode) {

  FDR_AXIS <<- fdr_mode   # set global for load_panel_data and draw_circular_manhattan
  axis_tag  <- if (fdr_mode) "_5fdr" else "_5"

  log_msg("  Loading data: ", tool, " BREEDING + NATURAL [FDR_AXIS=", fdr_mode, "]")
  dat_b <- load_panel_data(tool, "BREEDING")
  dat_n <- load_panel_data(tool, "NATURAL")

  draw_both <- function(dev_open_fn) {
    dev_open_fn()
    on.exit(dev.off(), add = TRUE)

    layout(matrix(c(1L, 2L, 3L, 3L), nrow = 2L, byrow = TRUE),
           widths  = lcm(c(PNL_W_CM, PNL_W_CM)),
           heights = lcm(c(PNL_H_CM, LEG_H_CM)))

    draw_circular_manhattan(dat_b, tool, "BREEDING",
                            open_device = FALSE, label_cex = 2.80)
    draw_circular_manhattan(dat_n, tool, "NATURAL",
                            open_device = FALSE, label_cex = 2.80)

    par(mar = c(0, 0, 0, 0))
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1))
    legend("center",
           legend    = CTX_ORDER,
           pch       = 16,
           col       = ctx_colors_global,
           pt.cex    = 2.0,
           cex       = 1.5,
           horiz     = TRUE,
           bty       = "n",
           x.intersp = 1.2)
  }

  # PDF (cairo_pdf -- vectorized)
  pdf_out <- paste0(out_base, ".pdf")
  log_msg("  Drawing PDF: ", basename(pdf_out))
  tryCatch(
    draw_both(function()
      cairo_pdf(pdf_out, width = COMB_W_CM / 2.54, height = COMB_H_CM / 2.54)),
    error = function(e) log_msg("  ERROR (pdf): ", conditionMessage(e))
  )
  if (check_nonempty(pdf_out)) {
    log_msg("  PDF OK (", file.info(pdf_out)$size, " bytes)")
  } else {
    log_msg("  WARNING: PDF may be empty -- check file size")
  }

  # SVG (vectorized)
  svg_out <- paste0(out_base, ".svg")
  log_msg("  Drawing SVG: ", basename(svg_out))
  tryCatch(
    draw_both(function()
      svg(svg_out, width = COMB_W_CM / 2.54, height = COMB_H_CM / 2.54)),
    error = function(e) log_msg("  ERROR (svg): ", conditionMessage(e))
  )
  if (check_nonempty(svg_out)) {
    log_msg("  SVG OK (", file.info(svg_out)$size, " bytes)")
  } else {
    log_msg("  WARNING: SVG may be empty")
  }

  # TIFF (ImageMagick assembly from individual per-cohort TIFFs)
  if (HAS_MAGICK) {
    tool_l  <- tolower(tool)
    img_b   <- file.path(MAN_DIR,
                  sprintf("manhattan_circular_%s_breeding%s.tiff", tool_l, axis_tag))
    img_n   <- file.path(MAN_DIR,
                  sprintf("manhattan_circular_%s_natural%s.tiff",  tool_l, axis_tag))
    leg_tif <- file.path(MAN_DIR, "legend_circular_5.tiff")

    if (all(file.exists(img_b, img_n, leg_tif))) {
      TMP <- file.path(PANEL_DIR, ".tmp19combined")
      dir.create(TMP, recursive = TRUE, showWarnings = FALSE)

      pair_tmp <- file.path(TMP, sprintf("%s%s_circs.tiff", tool_l, axis_tag))
      run_magick(paste(shQuote(img_b), shQuote(img_n),
                       "+append -compress lzw", shQuote(pair_tmp)))

      tiff_out <- paste0(out_base, ".tiff")
      run_magick(paste(shQuote(pair_tmp), shQuote(leg_tif),
                       "-append -compress lzw", shQuote(tiff_out)))
      log_msg("  TIFF assembled: ", file.info(tiff_out)$size, " bytes")

      png_out <- paste0(out_base, ".png")
      run_magick(paste("-density 150", shQuote(tiff_out), shQuote(png_out)))
      log_msg("  PNG: ", basename(png_out))

      eps_out <- paste0(out_base, ".eps")
      run_magick(paste("-density 300", shQuote(tiff_out), shQuote(eps_out)))
      log_msg("  EPS: ", basename(eps_out))

      unlink(TMP, recursive = TRUE)
    } else {
      log_msg("  WARNING: individual TIFFs missing -- TIFF/PNG/EPS skipped")
      log_msg("    Breeding: ", img_b)
      log_msg("    Natural:  ", img_n)
      log_msg("    Legend:   ", leg_tif)
    }
  }

  invisible(out_base)
}

############################################################
# 6) RUN -- 2 tools x 2 axis modes = 4 panels
############################################################

figure_specs <- list(
  list(tool = "GENESIS5",    fdr = FALSE,
       stem = "Figure3_circular_panel_genesis5_5"),
  list(tool = "MATRIXEQTL5", fdr = FALSE,
       stem = "EDF4_circular_panel_matrixeqtl5_5"),
  list(tool = "GENESIS5",    fdr = TRUE,
       stem = "Figure3_circular_panel_genesis5_5_FDR"),
  list(tool = "MATRIXEQTL5", fdr = TRUE,
       stem = "EDF4_circular_panel_matrixeqtl5_5_FDR")
)

for (spec in figure_specs) {
  log_msg(paste(rep("=", 65), collapse = ""))
  log_msg("Panel: ", spec$stem)
  out_base <- file.path(PANEL_DIR, spec$stem)
  draw_combined_panel(spec$tool, out_base, spec$fdr)
}

############################################################
# 7) COPY ALL FORMATS TO CORRECTED/FIGURES/NEW
############################################################

log_msg(paste(rep("=", 65), collapse = ""))
log_msg("Copying to CORRECTED/FIGURES/NEW...")

copy_fmts <- c("tiff", "pdf", "svg", "png", "eps")

for (spec in figure_specs) {
  for (fmt in copy_fmts) {
    src <- file.path(PANEL_DIR, paste0(spec$stem, ".", fmt))
    dst <- file.path(NEW_DIR,   paste0(fmt, "_", spec$stem, ".", fmt))
    if (file.exists(src) && file.info(src)$size > 0) {
      file.copy(src, dst, overwrite = TRUE)
      log_msg("  ", basename(dst))
    } else {
      log_msg("  MISSING/EMPTY (skip): ", basename(src))
    }
  }
}

writeLines(capture.output(sessionInfo()),
           file.path(LOG_DIR, "step19ab_combined_sessionInfo.txt"))

log_msg(paste(rep("=", 65), collapse = ""))
log_msg("Step 19ab combined finished")
log_msg("Panels in: ", PANEL_DIR)
log_msg("Copies in: ", NEW_DIR)
