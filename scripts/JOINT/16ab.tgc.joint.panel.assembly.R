#!/usr/bin/env Rscript
############################################################
# TreeGeneClimate (TGC) — JOINT ECS + TBS
# Step 16ab: Panel assembly — QQ plots and circular Manhattan plots
#
# PURPOSE
# Reads individually-saved TIFF files from steps 13ab and 14ab,
# assembles them into two-panel figures using the ImageMagick system binary
# (no R magick package required).
#
#   Figure 5           : circular Manhattan — GENESIS5    (breeding | natural)
#   Supplementary Fig 2: QQ panel          — GENESIS5    (2×3 grid)
#   Supplementary Fig 3: circular Manhattan — MatrixEQTL5 (breeding | natural)
#   Supplementary Fig 4: QQ panel          — MatrixEQTL5  (2×3 grid)
#
# QQ grid order (row 1 = breeding, row 2 = natural; cols = CpG, CHG, CHH):
#   a) breeding CpG   b) breeding CHG   c) breeding CHH
#   d) natural  CpG   e) natural  CHG   f) natural  CHH
#
# Circular panel order: a) breeding (left)   b) natural (right)
# Panel labels "a)" / "b)" are already embedded in the individual plots.
#
# USAGE
#   Rscript --vanilla 16ab.R
#
# INPUTS
#   RESULTS/JOINT/COMBINED5/qq/qqplot_<tool>_<cohort>_<context>.tiff
#   RESULTS/JOINT/COMBINED5/manhattan/manhattan_circular_<tool>_<cohort>.tiff
#
# OUTPUTS (RESULTS/JOINT/COMBINED5/panels/)
#   Figure5_circular_panel_genesis5.{tiff,png,pdf,eps}
#   SuppFig2_QQ_panel_genesis5.{tiff,png,pdf,eps}
#   SuppFig3_circular_panel_matrixeqtl5.{tiff,png,pdf,eps}
#   SuppFig4_QQ_panel_matrixeqtl5.{tiff,png,pdf,eps}
############################################################

options(stringsAsFactors = FALSE)

############################################################
# 1) PATHS
############################################################

# === USER CONFIGURATION ===
PROJECT_ROOT <- "/path/to/your/project"  # <-- set this
# ===========================
COMBINED5    <- file.path(PROJECT_ROOT, "RESULTS", "JOINT", "COMBINED5")
QQ_DIR       <- file.path(COMBINED5, "qq")
MAN_DIR      <- file.path(COMBINED5, "manhattan")
PANEL_DIR    <- file.path(COMBINED5, "panels")
dir.create(PANEL_DIR, recursive = TRUE, showWarnings = FALSE)

############################################################
# 2) HELPERS
############################################################

msg <- function(...) cat(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n"))

# Detect ImageMagick binary (prefer 'magick' IMv7, fall back to 'convert' IMv6)
MAGICK_BIN <- Sys.which("magick")
if (!nzchar(MAGICK_BIN)) {
  MAGICK_BIN <- Sys.which("convert")
  if (!nzchar(MAGICK_BIN))
    stop("ImageMagick not found in PATH. Load with:  module load imagemagick")
}
msg("ImageMagick: ", MAGICK_BIN)

# All ImageMagick calls go through run_magick to centralise error handling
run_magick <- function(args_str) {
  cmd <- paste(shQuote(MAGICK_BIN), args_str)
  ret <- system(cmd, ignore.stdout = TRUE, ignore.stderr = FALSE)
  if (ret != 0L) stop("magick failed:\n  ", cmd)
  invisible(ret)
}

# Create a white blank placeholder TIFF at the specified pixel dimensions
make_blank <- function(w, h, out) {
  run_magick(sprintf("-size %dx%d xc:white -compress lzw %s", w, h, shQuote(out)))
  out
}

# Horizontal append of a vector of TIFF paths → out
# ImageMagick +append places images side by side at the same height
hcat <- function(imgs, out) {
  run_magick(paste(paste(shQuote(imgs), collapse = " "),
                   "+append -compress lzw", shQuote(out)))
  out
}

# Vertical append of a vector of TIFF paths → out
# ImageMagick -append stacks images top-to-bottom at the same width
vcat <- function(imgs, out) {
  run_magick(paste(paste(shQuote(imgs), collapse = " "),
                   "-append -compress lzw", shQuote(out)))
  out
}

# Convert assembled TIFF to png / pdf / eps for journal submission flexibility
save_panel_formats <- function(tiff_path, out_base) {
  specs <- list(
    png = "",
    pdf = "-density 300",
    eps = "-density 300"
  )
  for (fmt in names(specs)) {
    out <- paste0(out_base, ".", fmt)
    tryCatch({
      run_magick(paste(specs[[fmt]], shQuote(tiff_path), shQuote(out)))
      msg("  Saved: ", basename(out))
    }, error = function(e) msg("  ERROR (", fmt, "): ", conditionMessage(e)))
  }
}

# Return existing TIFF path, or create and return a white placeholder.
# Placeholders preserve grid geometry when an individual panel file is missing.
safe_tiff <- function(path, blank_w, blank_h, tmp_dir) {
  if (file.exists(path)) return(path)
  msg("  MISSING (placeholder used): ", basename(path))
  blank <- file.path(tmp_dir, paste0("blank_", basename(path)))
  make_blank(blank_w, blank_h, blank)
  blank
}

############################################################
# 3) QQ PANELS  (2 rows × 3 cols per tool)
############################################################

msg("=== QQ panels ===")

# Temporary directory for intermediate row TIFFs; cleaned up in section 5
TMP <- file.path(PANEL_DIR, ".tmp")
dir.create(TMP, recursive = TRUE, showWarnings = FALSE)

qq_specs <- list(
  genesis5    = list(fig_tag = "SuppFig2", tool_l = "genesis5"),
  matrixeqtl5 = list(fig_tag = "SuppFig4", tool_l = "matrixeqtl5")
)

# Individual QQ TIFFs are 2400×2400 px (8 in × 8 in @ 300 dpi)
QQ_PX <- 2400L

for (nm in names(qq_specs)) {
  spec    <- qq_specs[[nm]]
  tool_l  <- spec$tool_l
  fig_tag <- spec$fig_tag

  # Enumerate all 6 panel TIFFs in reading order (row-major: breeding row, natural row)
  img_paths <- file.path(QQ_DIR,
    sprintf("qqplot_%s_%s_%s.tiff",
            tool_l,
            rep(c("breeding", "natural"), each = 3L),
            rep(c("cpg", "chg", "chh"), 2L)))

  # Replace any missing individual panels with same-sized white blanks
  imgs <- vapply(img_paths, safe_tiff,
                 blank_w = QQ_PX, blank_h = QQ_PX, tmp_dir = TMP,
                 FUN.VALUE = character(1L))

  # Build 2×3 grid: concatenate three panels per row, then stack the two rows
  row1_tmp <- file.path(TMP, sprintf("%s_row1.tiff", nm))
  row2_tmp <- file.path(TMP, sprintf("%s_row2.tiff", nm))
  hcat(imgs[1:3], row1_tmp)
  hcat(imgs[4:6], row2_tmp)

  out_tiff <- file.path(PANEL_DIR, sprintf("%s_QQ_panel_%s.tiff", fig_tag, tool_l))
  msg("Assembling: ", basename(out_tiff))
  vcat(c(row1_tmp, row2_tmp), out_tiff)
  msg("  Saved: ", basename(out_tiff))

  # Convert the master TIFF to additional formats for journal submission
  save_panel_formats(out_tiff, sub("\\.tiff$", "", out_tiff))
}

############################################################
# 4) CIRCULAR MANHATTAN PANELS  (1 row × 2 cols per tool)
############################################################

msg("=== Circular Manhattan panels ===")

man_specs <- list(
  genesis5    = list(fig_tag = "Figure5",  tool_l = "genesis5"),
  matrixeqtl5 = list(fig_tag = "SuppFig3", tool_l = "matrixeqtl5")
)

# Individual circular TIFFs are 24 cm × 24 cm @ 600 dpi ≈ 5669 × 5669 px
MAN_PX <- 5669L

for (nm in names(man_specs)) {
  spec    <- man_specs[[nm]]
  tool_l  <- spec$tool_l
  fig_tag <- spec$fig_tag

  img_b <- file.path(MAN_DIR, sprintf("manhattan_circular_%s_breeding.tiff", tool_l))
  img_n <- file.path(MAN_DIR, sprintf("manhattan_circular_%s_natural.tiff",  tool_l))

  imgs <- vapply(c(img_b, img_n), safe_tiff,
                 blank_w = MAN_PX, blank_h = MAN_PX, tmp_dir = TMP,
                 FUN.VALUE = character(1L))

  # Breeding and natural plots side by side, then legend strip appended below
  panel_tmp <- file.path(TMP, sprintf("%s_circs.tiff", nm))
  hcat(imgs, panel_tmp)
  leg_tiff <- file.path(MAN_DIR, "legend_circular.tiff")
  # Legend height ~827 px matches the 3.5 cm strip at 600 dpi saved by 14ab.R
  leg <- safe_tiff(leg_tiff, blank_w = MAN_PX * 2L, blank_h = 827L, tmp_dir = TMP)
  out_tiff <- file.path(PANEL_DIR,
    sprintf("%s_circular_panel_%s.tiff", fig_tag, tool_l))
  msg("Assembling: ", basename(out_tiff))
  vcat(c(panel_tmp, leg), out_tiff)
  msg("  Saved: ", basename(out_tiff))

  save_panel_formats(out_tiff, sub("\\.tiff$", "", out_tiff))
}

############################################################
# 5) CLEANUP
############################################################

# Remove intermediate row/column TIFFs; final panels in PANEL_DIR are kept
unlink(TMP, recursive = TRUE)
msg("Step 16ab finished — panels saved to: ", PANEL_DIR)

sessionInfo()
