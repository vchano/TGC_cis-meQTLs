#!/usr/bin/env Rscript
############################################################
# TreeGeneClimate (TGC) — JOINT ECS + TMS
# Step 19ab v2: Combined vectorized circular Manhattan panels
#               Version 4 — 24cm canvas, wider tracks (TRACK_HEIGHT 0.25), tighter margins,
#                            white bg, black labels. Ring fills more canvas without scaling elements.
#
# Reads *_4 individual TIFFs from 16ab4.R, assembles combined panels,
# saves with "_4" suffix so prior outputs are not overwritten.
#
# OUTPUTS (RESULTS/JOINT/COMBINED5/panels/)
#   Figure6_circular_panel_genesis5_4.{pdf,svg,tiff,png,eps}
#   SuppFig2_circular_panel_matrixeqtl5_4.{pdf,svg,tiff,png,eps}
#
# OUTPUTS (RESULTS/CORRECTED/FIGURES/NEW/)
#   {fmt}_corrected_Figure6_GENESIS_Manhattan_panel_4.{fmt}
#   {fmt}_corrected_FigureS2_MatrixEQTL_Manhattan_panel_4.{fmt}
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
# 1) SOURCE 16ab4.R (functions + updated settings; guard skips main loop)
############################################################

SOURCED_16AB4 <- TRUE
source("/path/to/your/project/SCRIPTS/JOINT/16ab4.R")

############################################################
# 2) PATHS
############################################################

PROJECT_ROOT <- "/path/to/your/project"
MAN_DIR      <- file.path(PROJECT_ROOT, "RESULTS", "JOINT", "COMBINED5", "manhattan")
PANEL_DIR    <- file.path(PROJECT_ROOT, "RESULTS", "JOINT", "COMBINED5", "panels")
NEW_DIR      <- file.path(PROJECT_ROOT, "RESULTS", "CORRECTED", "FIGURES", "NEW")

dir.create(PANEL_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(NEW_DIR,   recursive = TRUE, showWarnings = FALSE)

LOG_DIR <- file.path(PROJECT_ROOT, "RESULTS", "JOINT", "COMBINED5", "LOGS")
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)
LOGFILE <- file.path(LOG_DIR, "step19ab_v4.log")
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
# 4) COMBINED DEVICE DIMENSIONS  (uses OUT_W_CM = 24 from 16ab4.R, wider tracks)
############################################################

PNL_W_CM  <- OUT_W_CM
PNL_H_CM  <- OUT_H_CM
LEG_H_CM  <- 3.5
COMB_W_CM <- PNL_W_CM * 2
COMB_H_CM <- PNL_H_CM + LEG_H_CM

MAGICK_BIN <- Sys.which("magick")
if (!nzchar(MAGICK_BIN)) MAGICK_BIN <- Sys.which("convert")
HAS_MAGICK <- nzchar(MAGICK_BIN)
if (HAS_MAGICK) log_msg("ImageMagick: ", MAGICK_BIN) else
  log_msg("WARNING: ImageMagick not found — raster assembly skipped")

run_magick <- function(args_str) {
  cmd <- paste(shQuote(MAGICK_BIN), args_str)
  ret <- system(cmd, ignore.stdout = TRUE, ignore.stderr = FALSE)
  if (ret != 0L) stop("magick failed:\n  ", cmd)
}

############################################################
# 5) DRAW-COMBINED FUNCTION
############################################################

draw_combined_panel <- function(tool, out_base) {

  dat_b <- load_panel_data(tool, "BREEDING")
  dat_n <- load_panel_data(tool, "NATURAL")
  ctx_colors_global <- setNames(
    brewer.pal(max(3L, length(CTX_ORDER)), POINT_PALETTE)[seq_along(CTX_ORDER)],
    CTX_ORDER
  )

  draw_both <- function(dev_open_fn) {
    dev_open_fn()
    on.exit(dev.off(), add = TRUE)

    layout(matrix(c(1L, 2L, 3L, 3L), nrow = 2L, byrow = TRUE),
           widths  = lcm(c(PNL_W_CM, PNL_W_CM)),
           heights = lcm(c(PNL_H_CM, LEG_H_CM)))

    draw_circular_manhattan(dat_b, tool, "BREEDING",
                            open_device = FALSE, label_cex = 3.15)
    draw_circular_manhattan(dat_n, tool, "NATURAL",
                            open_device = FALSE, label_cex = 3.15)

    par(mar = c(0, 0, 0, 0))
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1))
    legend("center",
           legend    = CTX_ORDER,
           pch       = 16,
           col       = ctx_colors_global,
           pt.cex    = 3.0,
           cex       = 2.0,
           horiz     = TRUE,
           bty       = "n",
           x.intersp = 1.5)
  }

  # Vectorized PDF
  pdf_out <- paste0(out_base, ".pdf")
  log_msg("  Drawing vectorized PDF: ", basename(pdf_out))
  tryCatch(
    draw_both(function()
      cairo_pdf(pdf_out, width = COMB_W_CM / 2.54, height = COMB_H_CM / 2.54)),
    error = function(e) log_msg("  ERROR (pdf): ", conditionMessage(e))
  )
  log_msg("  Saved: ", pdf_out)

  # Vectorized SVG
  svg_out <- paste0(out_base, ".svg")
  log_msg("  Drawing vectorized SVG: ", basename(svg_out))
  tryCatch(
    draw_both(function()
      svg(svg_out, width = COMB_W_CM / 2.54, height = COMB_H_CM / 2.54)),
    error = function(e) log_msg("  ERROR (svg): ", conditionMessage(e))
  )
  log_msg("  Saved: ", svg_out)

  # Raster assembly from _4 individual TIFFs
  if (HAS_MAGICK) {
    tool_l   <- tolower(tool)
    img_b    <- file.path(MAN_DIR, sprintf("manhattan_circular_%s_breeding_4.tiff", tool_l))
    img_n    <- file.path(MAN_DIR, sprintf("manhattan_circular_%s_natural_4.tiff",  tool_l))
    leg_tiff <- file.path(MAN_DIR, "legend_circular_4.tiff")

    if (all(file.exists(img_b, img_n, leg_tiff))) {
      TMP <- file.path(PANEL_DIR, ".tmp19v4")
      dir.create(TMP, recursive = TRUE, showWarnings = FALSE)

      panel_tmp <- file.path(TMP, sprintf("%s_circs_4.tiff", tool_l))
      run_magick(paste(shQuote(img_b), shQuote(img_n),
                       "+append -compress lzw", shQuote(panel_tmp)))

      tiff_out <- paste0(out_base, ".tiff")
      run_magick(paste(shQuote(panel_tmp), shQuote(leg_tiff),
                       "-append -compress lzw", shQuote(tiff_out)))
      log_msg("  Saved: ", tiff_out)

      for (fmt in c("png", "eps")) {
        out_fmt <- paste0(out_base, ".", fmt)
        run_magick(paste("-density 300", shQuote(tiff_out), shQuote(out_fmt)))
        log_msg("  Saved: ", out_fmt)
      }
      unlink(TMP, recursive = TRUE)
    } else {
      log_msg("  WARNING: pre-generated _4 TIFFs missing — raster formats skipped")
    }
  }

  invisible(out_base)
}

############################################################
# 6) RUN FOR EACH TOOL
############################################################

figure_specs <- list(
  GENESIS5    = list(out_stem = "Figure6_circular_panel_genesis5_4"),
  MATRIXEQTL5 = list(out_stem = "SuppFig2_circular_panel_matrixeqtl5_4")
)

for (tool in names(figure_specs)) {
  log_msg("=== ", tool, " ===")
  spec     <- figure_specs[[tool]]
  out_base <- file.path(PANEL_DIR, spec$out_stem)
  draw_combined_panel(tool, out_base)
}

############################################################
# 7) COPY TO CORRECTED/FIGURES/NEW WITH "_4" SUFFIX
############################################################

log_msg("Copying to CORRECTED/FIGURES/NEW/...")

copy_formats <- c("pdf", "svg", "tiff", "png", "eps")

copies <- list(
  list(
    stem = file.path(PANEL_DIR, "Figure6_circular_panel_genesis5_4"),
    dest = "corrected_Figure6_GENESIS_Manhattan_panel_4"
  ),
  list(
    stem = file.path(PANEL_DIR, "SuppFig2_circular_panel_matrixeqtl5_4"),
    dest = "corrected_FigureS2_MatrixEQTL_Manhattan_panel_4"
  )
)

for (cp in copies) {
  for (fmt in copy_formats) {
    src <- paste0(cp$stem, ".", fmt)
    dst <- file.path(NEW_DIR, paste0(fmt, "_", cp$dest, ".", fmt))
    if (file.exists(src)) {
      file.copy(src, dst, overwrite = TRUE)
      log_msg("  ", basename(dst))
    } else {
      log_msg("  MISSING (skip): ", basename(src))
    }
  }
}

############################################################
# 8) DONE
############################################################

log_msg("Step 19ab v2 finished")
log_msg("Outputs in: ", PANEL_DIR)
log_msg("Corrected copies in: ", NEW_DIR)
