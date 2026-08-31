#!/usr/bin/env Rscript
############################################################
# TreeGeneClimate (TGC) — JOINT ECS + TBS
# Step 16ab: Circular Manhattan plots (GENESIS5 + MatrixEQTL5)
#            Version 5 — 18 cm canvas, BH-FDR axis option
#
# PURPOSE
# Produces circular Manhattan plots from cis-meQTL results (meQTL5 run).
# One plot per tool × cohort × axis mode (8 total):
#   GENESIS5     BREEDING/NATURAL  x  raw-p / BH-FDR
#   MATRIXEQTL5  BREEDING/NATURAL  x  raw-p / BH-FDR
# Three rings: CpG (outer), CHG (middle), CHH (inner).
# Dashed threshold line at BH-FDR < 1e-10.
#
# v5 changes vs v4:
#   Canvas: 18 cm x 18 cm (was 24)     CHR_LABEL_CEX: 1.20 (was 0.975)
#   Dot cex: 0.22 (was 0.15)           Margins: c(0.5, 2, 1.0, 2)
#   Legend height: 2.0 cm (was 3.5)    FDR_AXIS flag (new)
#   Y-axis tick/title cex: 0.55/0.65   (larger caused overlap on 18 cm canvas)
#
# USAGE
#   Rscript --vanilla 16ab.tgc.joint.panel.assembly.R
#   # or source from 19ab with: SOURCED_16AB4 <- TRUE; FDR_AXIS <- FALSE
#
# INPUTS
#   RESULTS/JOINT/GENESIS5/MEQTL/<COHORT>/<CONTEXT>/cis_meqtl_all_results.rds
#   RESULTS/JOINT/MATRIXEQTL5/<COHORT>/<CONTEXT>/cis_meqtl_all_results.rds
#   RESULTS/JOINT/MQTL5/INPUTS/<COHORT>/<CONTEXT>/snp_variant_annot.rds
#   GDS files (chromosome layout)
#   DATA/METADATA/chrom.list  (optional; auto-built if absent)
#
# OUTPUTS (RESULTS/JOINT/COMBINED5/manhattan/)
#   manhattan_circular_<tool>_<cohort>_5.{tiff,pdf,svg,eps}      -- raw p-value
#   manhattan_circular_<tool>_<cohort>_5fdr.{tiff,pdf,svg,eps}   -- BH-FDR axis
#   legend_circular_5.tiff      (36 cm, for 2-panel combined figure)
#   legend_circular_5_wide.tiff (72 cm, for 4-panel comparison)
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
# 1) SETTINGS
############################################################

COHORTS        <- c("BREEDING", "NATURAL")
CTX_ORDER      <- c("CpG", "CHG", "CHH")   # outer -> inner rings
TOOLS          <- c("GENESIS5", "MATRIXEQTL5")
FDR_STRICT     <- 1e-10

# FDR_AXIS: when TRUE, y-axis shows -log10(BH-FDR); threshold constant at 10.
# Set before sourcing, or it defaults to FALSE.
if (!exists("FDR_AXIS")) FDR_AXIS <- FALSE

# v5 rendering settings
OUT_W_CM      <- 18; OUT_H_CM <- 18; OUT_RES <- 600
AXIS_GAP_FRAC <- 0.80    # wide gap for y-axis label
TRACK_HEIGHT  <- 0.25
POINT_PALETTE <- "Set1"
CHR_LABEL_CEX <- 1.20

############################################################
# 2) PATHS  -- set PROJECT_ROOT and RDATA_DIR to your paths
############################################################

PROJECT_ROOT <- "/path/to/your/project"   # <-- set this

RESULT_ROOTS <- list(
  GENESIS5    = file.path(PROJECT_ROOT, "RESULTS", "JOINT", "GENESIS5",    "MEQTL"),
  MATRIXEQTL5 = file.path(PROJECT_ROOT, "RESULTS", "JOINT", "MATRIXEQTL5")
)

ANNOT_ROOT <- file.path(PROJECT_ROOT, "RESULTS", "JOINT", "MQTL5", "INPUTS")

OUT_MAN <- file.path(PROJECT_ROOT, "RESULTS", "JOINT", "COMBINED5", "manhattan")
LOG_DIR <- file.path(PROJECT_ROOT, "RESULTS", "JOINT", "COMBINED5", "LOGS")
dir.create(OUT_MAN, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

LOGFILE <- file.path(LOG_DIR, "step16ab_panel_assembly.log")
if (file.exists(LOGFILE)) file.remove(LOGFILE)

# Directory containing breeding.snp.gds and natural.snp.gds
RDATA_DIR <- "/path/to/your/RDATA"        # <-- set this
GDS_FILES <- list(
  BREEDING = file.path(RDATA_DIR, "breeding.snp.gds"),
  NATURAL  = file.path(RDATA_DIR, "natural.snp.gds")
)
CHROM_MAP_FILE <- file.path(PROJECT_ROOT, "DATA", "METADATA", "chrom.list")

############################################################
# 3) HELPERS
############################################################

log_msg <- function(...) {
  txt <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ",
                paste0(..., collapse = ""))
  cat(txt, "\n")
  cat(txt, "\n", file = LOGFILE, append = TRUE)
}
sep_line <- function() log_msg(paste(rep("=", 70), collapse = ""))

pval_column <- function(res_dt) {
  candidates <- c("pvalue", "pval", "Score.pval",
                  grep("pval$", names(res_dt), value = TRUE))
  candidates <- candidates[candidates %in% names(res_dt)]
  if (!length(candidates)) stop("No p-value column found")
  candidates[1]
}

find_pval_at_fdr <- function(pvals, fdr_thresh) {
  pvals <- sort(pvals[is.finite(pvals) & !is.na(pvals) & pvals > 0])
  if (!length(pvals)) return(NA_real_)
  padj <- p.adjust(pvals, method = "BH")
  idx  <- which(padj < fdr_thresh)
  if (!length(idx)) return(NA_real_)
  pvals[max(idx)]
}

read_gds_snppos <- function(gds_path) {
  if (!file.exists(gds_path)) return(NULL)
  gf <- snpgdsOpen(gds_path); on.exit(snpgdsClose(gf))
  data.table(
    snp = as.character(read.gdsn(index.gdsn(gf, "snp.id"))),
    chr = as.character(read.gdsn(index.gdsn(gf, "snp.chromosome"))),
    pos = as.integer (read.gdsn(index.gdsn(gf, "snp.position")))
  )
}

build_contig_lengths <- function(gds_paths) {
  all_dt <- rbindlist(lapply(gds_paths, read_gds_snppos), fill = TRUE)
  if (!nrow(all_dt)) stop("No contig positions from GDS files")
  cl <- all_dt[, .(len = max(pos, na.rm = TRUE)), by = chr]
  cl[!is.na(chr)]
}

read_chrom_map <- function(path) {
  if (!file.exists(path)) return(NULL)
  cm <- fread(path, data.table = FALSE)
  if (!all(c("Real_chrom", "ChromN") %in% colnames(cm))) return(NULL)
  cm$Real_chrom <- as.character(cm$Real_chrom)
  cm$ChromN     <- as.character(cm$ChromN)
  cm$ChromN[cm$ChromN == "Chr13"] <- "ChrUn"
  cm
}

build_default_chrom_map <- function(contig_len_dt) {
  dt <- copy(contig_len_dt); setorder(dt, -len)
  n  <- min(12L, nrow(dt))
  dt[, ChromN := ifelse(seq_len(.N) <= n, sprintf("Chr%02d", seq_len(.N)), "ChrUn")]
  data.frame(Real_chrom = dt$chr, ChromN = dt$ChromN, stringsAsFactors = FALSE)
}

prepare_layout <- function(contig_len_dt, chrom_map_df, gap_frac = AXIS_GAP_FRAC) {
  cm_dt <- as.data.table(chrom_map_df)
  dt    <- merge(contig_len_dt, cm_dt, by.x = "chr", by.y = "Real_chrom",
                 all.x = TRUE, sort = FALSE)
  dt[is.na(ChromN), ChromN := "ChrUn"]
  desired <- c("ChrUn", sprintf("Chr%02d", 1:12))
  lvls    <- c(intersect(desired, unique(dt$ChromN)),
               sort(setdiff(unique(dt$ChromN), desired)))
  dt[, ChromN := factor(ChromN, levels = lvls)]; setorder(dt, ChromN, chr)
  dt[, offset_within_chromN := {
    tmp <- cumsum(len); if (.N == 0) numeric(0) else c(0, tmp[-.N])
  }, by = ChromN]
  clen <- dt[, .(chromN_len = sum(len)), by = ChromN]
  clen[, ChromN := factor(ChromN, levels = lvls)]; setorder(clen, ChromN)
  gap_len <- max(clen$chromN_len) * gap_frac * 0.5
  idx_un  <- which(lvls == "ChrUn")
  if (length(idx_un)) {
    lvls_g <- append(lvls, "axis_panel2", after = idx_un)
    sl_lst <- as.list(setNames(clen$chromN_len, as.character(clen$ChromN)))
    sl_lst <- append(sl_lst, list(axis_panel2 = gap_len), after = idx_un)
  } else {
    lvls_g <- c("axis_panel2", lvls)
    sl_lst <- c(list(axis_panel2 = gap_len),
                as.list(setNames(clen$chromN_len, as.character(clen$ChromN))))
  }
  sl <- unlist(sl_lst[lvls_g])
  if ("Chr01" %in% names(sl)) {
    idx <- which(names(sl) == "Chr01")
    nms <- names(sl); sl <- sl[c(nms[idx:length(nms)], nms[seq_len(idx - 1)])]
  }
  list(contig_info = dt, chrom_levels = lvls, sector_lens = sl)
}

map_to_chromN <- function(res_dt, contig_info, chrom_map_df, pval_col) {
  dt <- copy(res_dt)
  if (pval_col != "pval" && pval_col %in% names(dt)) dt[, pval := get(pval_col)]
  if (!"snp_pos" %in% names(dt)) return(data.table())
  if ("snp_chr" %in% names(dt)) {
    dt[, chr := snp_chr]
  } else {
    dt[, chr := NA_character_]
  }
  dt[, pos := as.numeric(snp_pos)]
  dt <- dt[!is.na(pos) & is.finite(pval) & pval > 0 & !is.na(chr)]
  if (!nrow(dt)) return(data.table())
  cm_dt <- as.data.table(chrom_map_df)
  dt[, Real_chrom := chr]
  dt <- merge(dt, cm_dt, by.x = "Real_chrom", by.y = "Real_chrom",
              all.x = TRUE, sort = FALSE)
  dt[is.na(ChromN), ChromN := "ChrUn"]
  setkey(contig_info, chr)
  out <- contig_info[dt, on = .(chr)]
  out[is.na(offset_within_chromN), offset_within_chromN := 0]
  out[, pos_within_chromN := offset_within_chromN + pos]
  # FDR_AXIS: use pre-computed BH-FDR column if available, else raw p
  if (FDR_AXIS && "fdr__" %in% names(out)) {
    out[, logp := -log10(fdr__)]
  } else {
    out[, logp := -log10(pval)]
  }
  out
}

############################################################
# 4) DATA LOADING (separated from plotting to allow multiple device types)
############################################################

load_panel_data <- function(tool, cohort) {
  ctx_colors <- setNames(
    brewer.pal(max(3, length(CTX_ORDER)), POINT_PALETTE)[seq_along(CTX_ORDER)],
    CTX_ORDER
  )

  panel_data       <- list()
  threshold_logp   <- list()
  max_logp_per_ctx <- list()

  for (ctx in CTX_ORDER) {
    rpath <- file.path(RESULT_ROOTS[[tool]], cohort, ctx, "cis_meqtl_all_results.rds")
    if (!file.exists(rpath)) {
      log_msg("  Missing RDS: ", rpath)
      panel_data[[ctx]] <- data.table()
      threshold_logp[[ctx]] <- NA_real_
      max_logp_per_ctx[[ctx]] <- 10
      next
    }

    log_msg("  Loading: ", tool, " / ", cohort, " / ", ctx,
            " [FDR_AXIS=", FDR_AXIS, "]")
    res_dt   <- as.data.table(readRDS(rpath))
    pval_col <- tryCatch(pval_column(res_dt), error = function(e) NA_character_)
    if (is.na(pval_col)) {
      log_msg("  No p-value column in ", rpath)
      panel_data[[ctx]] <- data.table()
      threshold_logp[[ctx]] <- NA_real_
      max_logp_per_ctx[[ctx]] <- 10
      next
    }
    log_msg("  p-value column: '", pval_col, "' | nrow=", nrow(res_dt))

    pvals <- res_dt[[pval_col]]

    # Compute BH-FDR on the full context test set before position filtering.
    # NOTE: FDR_AXIS=TRUE on NATURAL/CHH (~500M tests) requires ~96 GB RAM.
    if (FDR_AXIS) {
      log_msg("  Computing BH-FDR on ", length(pvals), " tests...")
      fdrs <- p.adjust(pvals, method = "BH")
      res_dt[, fdr__ := fdrs]
      rm(fdrs); gc()
      log_msg("  BH-FDR done")
    }

    # Rebuild snp_chr and snp_pos from snp_variant_annot.rds (authoritative coords).
    # Drops any stale position/chr columns first.
    for (col_drop in intersect(c("chr", "snp_chr", "pos", "snp_pos"), names(res_dt)))
      res_dt[, (col_drop) := NULL]

    var_path <- file.path(ANNOT_ROOT, cohort, ctx, "snp_variant_annot.rds")
    if (file.exists(var_path) && "snp" %in% names(res_dt)) {
      va <- as.data.table(readRDS(var_path))
      log_msg("  snp_variant_annot cols: ", paste(names(va), collapse = ", "))
      res_dt[, snp   := as.character(snp)]
      va[,   snp_id := as.character(snp_id)]
      if ("pos" %in% names(va)) {
        res_dt <- merge(res_dt,
                        va[, .(snp_id, snp_chr = chr, snp_pos = as.numeric(pos))],
                        by.x = "snp", by.y = "snp_id", all.x = TRUE)
        log_msg("  snp_pos matched: ", sum(!is.na(res_dt$snp_pos)), " / ", nrow(res_dt))
      } else {
        log_msg("  WARNING: snp_variant_annot has no 'pos' column")
        res_dt <- merge(res_dt, va[, .(snp_id, snp_chr = chr)],
                        by.x = "snp", by.y = "snp_id", all.x = TRUE)
      }
      n_matched <- sum(!is.na(res_dt$snp_chr))
      log_msg("  snp_chr matched: ", n_matched, " / ", nrow(res_dt),
              " (", round(100 * n_matched / max(1, nrow(res_dt)), 1), "%)")
    } else {
      if (!file.exists(var_path))
        log_msg("  WARNING: snp_variant_annot.rds not found: ", var_path)
      if (!"snp" %in% names(res_dt))
        log_msg("  WARNING: no 'snp' column -- cannot add snp_chr")
      res_dt[, snp_chr := NA_character_]
    }

    # Threshold logp
    if (FDR_AXIS) {
      threshold_logp[[ctx]] <- -log10(FDR_STRICT)   # constant = 10 for all contexts
    } else {
      p_thr <- find_pval_at_fdr(pvals, FDR_STRICT)
      threshold_logp[[ctx]] <- if (is.na(p_thr)) NA_real_ else -log10(p_thr)
    }
    log_msg("  Threshold logp: ",
            if (is.na(threshold_logp[[ctx]])) "none" else round(threshold_logp[[ctx]], 2))

    mapped <- tryCatch(
      map_to_chromN(res_dt, contig_info_dt, chrom_map_df, pval_col),
      error = function(e) {
        log_msg("  map_to_chromN ERROR: ", conditionMessage(e))
        data.table()
      }
    )
    log_msg("  Mapped points: ", nrow(mapped))
    if (nrow(mapped) == 0 && "snp_chr" %in% names(res_dt)) {
      log_msg("  DIAG snp_chr sample  : ",
              paste(head(unique(na.omit(res_dt$snp_chr)), 5), collapse = ", "))
      log_msg("  DIAG contig_info chr : ",
              paste(head(unique(contig_info_dt$chr), 5), collapse = ", "))
    }

    max_logp_per_ctx[[ctx]] <- if (nrow(mapped) > 0 && any(is.finite(mapped$logp))) {
      max(mapped$logp[is.finite(mapped$logp)])
    } else { 10 }

    panel_data[[ctx]] <- mapped
    rm(res_dt, pvals); gc()
  }

  list(panel_data       = panel_data,
       threshold_logp   = threshold_logp,
       max_logp_per_ctx = max_logp_per_ctx,
       ctx_colors       = ctx_colors)
}

############################################################
# 5) CIRCULAR PLOT FUNCTION
############################################################

draw_circular_manhattan <- function(dat, tool, cohort, out_path = NULL,
                                    open_device = TRUE, label_cex = 2.80) {

  panel_data       <- dat$panel_data
  threshold_logp   <- dat$threshold_logp
  max_logp_per_ctx <- dat$max_logp_per_ctx
  ctx_colors       <- dat$ctx_colors

  if (open_device) {
    if (grepl("\\.pdf$", out_path, ignore.case = TRUE)) {
      pdf(out_path, width = OUT_W_CM / 2.54, height = OUT_H_CM / 2.54)
    } else if (grepl("\\.svg$", out_path, ignore.case = TRUE)) {
      svg(out_path, width = OUT_W_CM / 2.54, height = OUT_H_CM / 2.54)
    } else if (grepl("\\.eps$", out_path, ignore.case = TRUE)) {
      cairo_ps(out_path, width = OUT_W_CM / 2.54, height = OUT_H_CM / 2.54)
    } else if (grepl("\\.png$", out_path, ignore.case = TRUE)) {
      png(out_path, width = OUT_W_CM, height = OUT_H_CM, units = "cm", res = 150)
    } else {
      tiff(out_path, width = OUT_W_CM, height = OUT_H_CM,
           units = "cm", res = OUT_RES, compression = "lzw")
    }
    on.exit({ circos.clear(); dev.off() }, add = TRUE)
  } else {
    on.exit(circos.clear(), add = TRUE)
  }

  par(mar = c(0.5, 2, 1.0, 2))   # reduced: less space above A/B label, below legend

  circos.par(
    cell.padding            = c(0, 0, 0, 0),
    track.margin            = c(0.005, 0.005),
    gap.after               = rep(0.5, length(sector_lens)),
    start.degree            = 90,
    clock.wise              = TRUE,
    points.overflow.warning = FALSE
  )
  circos.initialize(factors = names(sector_lens),
                    xlim    = cbind(rep(0, length(sector_lens)),
                                   unname(sector_lens)))

  # Chromosome label track
  circos.track(ylim = c(0, 1), track.height = 0.05,
               bg.border = NA, bg.col = "white",
    panel.fun = function(x, y) {
      sn <- CELL_META$sector.index
      if (sn == "axis_panel2") return(invisible())
      lbl <- if (sn == "ChrUn") "Un" else sn
      circos.text(CELL_META$xcenter, 0.5, lbl,
                  facing = "bending.inside", niceFacing = TRUE,
                  cex = CHR_LABEL_CEX, col = "black")
    })

  yaxis_label <- if (FDR_AXIS) "-log10(FDR)" else "-log10(p)"

  # Explicit per-track calls (not a loop) to avoid R closure-in-loop issues.
  add_data_track <- function(dat_ctx, col_ctx, thr_logp, ylim_top) {
    circos.track(ylim = c(0, ylim_top), track.height = TRACK_HEIGHT,
                 bg.border = "grey55", bg.col = "white",
      panel.fun = function(x, y) {
        sn <- CELL_META$sector.index

        if (sn == "axis_panel2") {
          circos.rect(CELL_META$xlim[1], CELL_META$ylim[1],
                      CELL_META$xlim[2], CELL_META$ylim[2],
                      col = "white", border = "white", lwd = 3)

          ax_x     <- CELL_META$xlim[2]
          at_vals  <- pretty(c(0, ylim_top), n = 3)
          at_vals  <- at_vals[at_vals >= 0 & at_vals <= ylim_top]
          tick_len <- CELL_META$xrange * 0.12

          circos.lines(c(ax_x, ax_x), c(0, ylim_top), lwd = 0.7, straight = TRUE)
          for (v in at_vals) {
            circos.lines(c(ax_x - tick_len, ax_x), c(v, v), lwd = 0.6, straight = TRUE)
            circos.text(ax_x - tick_len * 2.0, v,
                        labels = as.character(v),
                        adj = c(1, 0.5), cex = 0.55, font = 1, col = "black",
                        facing = "bending.inside", niceFacing = TRUE)
          }
          circos.text(ax_x - CELL_META$xrange * 0.80, ylim_top * 0.5,
                      labels = yaxis_label,
                      facing = "clockwise", niceFacing = FALSE,
                      adj = c(0.5, 0.5), cex = 0.65, col = "black")
          return(invisible())
        }

        if (!nrow(dat_ctx)) return(invisible())
        sub <- dat_ctx[ChromN == sn & !is.na(pos_within_chromN) & !is.na(logp)]
        if (!nrow(sub)) return(invisible())

        circos.points(sub$pos_within_chromN, sub$logp,
                      pch = 16, cex = 0.22,
                      col = adjustcolor(col_ctx, alpha.f = 0.4))

        if (is.finite(thr_logp) && thr_logp > 0 && thr_logp <= ylim_top)
          circos.lines(c(0, CELL_META$xlim[2]),
                       c(thr_logp, thr_logp),
                       lty = 2, lwd = 0.8, col = "black")
      })
  }

  ymax_cpg <- max_logp_per_ctx[["CpG"]]; if (!is.finite(ymax_cpg) || ymax_cpg <= 0) ymax_cpg <- 10
  ymax_chg <- max_logp_per_ctx[["CHG"]]; if (!is.finite(ymax_chg) || ymax_chg <= 0) ymax_chg <- 10
  ymax_chh <- max_logp_per_ctx[["CHH"]]; if (!is.finite(ymax_chh) || ymax_chh <= 0) ymax_chh <- 10

  add_data_track(panel_data[["CpG"]], ctx_colors[["CpG"]], threshold_logp[["CpG"]], ymax_cpg * 1.05)
  add_data_track(panel_data[["CHG"]], ctx_colors[["CHG"]], threshold_logp[["CHG"]], ymax_chg * 1.05)
  add_data_track(panel_data[["CHH"]], ctx_colors[["CHH"]], threshold_logp[["CHH"]], ymax_chh * 1.05)

  # Panel label: A = BREEDING, B = NATURAL
  fig_region  <- par("fig")
  fig_x1 <- fig_region[1]; fig_x2 <- fig_region[2]
  fig_y1 <- fig_region[3]; fig_y2 <- fig_region[4]
  label_x_npc <- fig_x1 + 0.02 * (fig_x2 - fig_x1)
  label_y_npc <- fig_y1 + 0.98 * (fig_y2 - fig_y1)
  upViewport(0)

  panel_label <- if (cohort == "BREEDING") "A" else "B"
  grid.text(panel_label,
            x = unit(label_x_npc, "npc"), y = unit(label_y_npc, "npc"),
            just = c("left", "top"),
            gp = gpar(cex = label_cex, fontface = "bold"))

  invisible(out_path)
}

############################################################
# 6) BUILD CHROMOSOME LAYOUT
############################################################

# Guard: skip main loop when sourced by 19ab.tgc.joint.circular.manhattan.combined.R
if (!exists("SOURCED_16AB4") || !isTRUE(SOURCED_16AB4)) {

sep_line()
log_msg("Step 16ab v5 -- Circular Manhattan plots (meQTL5)")
log_msg("Canvas: ", OUT_W_CM, " cm x ", OUT_H_CM, " cm @ ", OUT_RES, " dpi")
log_msg("Building chromosome layout from GDS files...")

contig_len_dt <- build_contig_lengths(unlist(GDS_FILES))
log_msg("Contigs found: ", nrow(contig_len_dt))

chrom_map_df <- read_chrom_map(CHROM_MAP_FILE)
if (is.null(chrom_map_df)) {
  log_msg("No chrom.list -- building default map from GDS (sorted by length)")
  chrom_map_df <- build_default_chrom_map(contig_len_dt)
  fwrite(as.data.table(chrom_map_df),
         file.path(OUT_MAN, "chrom_map_auto.tsv"), sep = "\t")
}

layout_info    <- prepare_layout(contig_len_dt, chrom_map_df)
contig_info_dt <- layout_info$contig_info
sector_lens    <- layout_info$sector_lens
log_msg("Sector order: ", paste(names(sector_lens), collapse = ", "))

############################################################
# 7) DRAW ALL PLOTS -- both p-value and FDR axis modes
############################################################

for (fdr_mode in c(FALSE, TRUE)) {
  FDR_AXIS <- fdr_mode
  suffix   <- if (fdr_mode) "_5fdr" else "_5"
  sep_line()
  log_msg("=== Axis mode: ", if (fdr_mode) "BH-FDR" else "raw p-value",
          " | suffix: ", suffix, " ===")

  for (tool in TOOLS) {
    for (cohort in COHORTS) {
      sep_line()
      log_msg("Processing: ", tool, " | ", cohort)

      dat <- load_panel_data(tool, cohort)

      out_base <- file.path(OUT_MAN,
        paste0("manhattan_circular_", tolower(tool), "_", tolower(cohort), suffix))

      for (fmt in c("tiff", "pdf", "svg", "eps")) {
        out_path <- paste0(out_base, ".", fmt)
        log_msg("  Drawing (", fmt, "): ", basename(out_path))
        tryCatch(
          draw_circular_manhattan(dat, tool, cohort, out_path),
          error = function(e) log_msg("  ERROR (", fmt, "): ", conditionMessage(e))
        )
        log_msg("  Saved: ", out_path)
      }
    }
  }
}

############################################################
# 8) LEGENDS
############################################################

ctx_colors_global <- setNames(
  brewer.pal(max(3, length(CTX_ORDER)), POINT_PALETTE)[seq_along(CTX_ORDER)],
  CTX_ORDER
)

draw_legend <- function(path, width_cm, height_cm = 2.0,
                        pt_cex = 2.0, text_cex = 1.5) {
  tiff(path, width = width_cm, height = height_cm,
       units = "cm", res = OUT_RES, compression = "lzw")
  par(mar = c(0, 0, 0, 0))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  legend("center", legend = CTX_ORDER, pch = 16, col = ctx_colors_global,
         pt.cex = pt_cex, cex = text_cex, horiz = TRUE, bty = "n",
         x.intersp = 1.2)
  dev.off()
}

leg_36 <- file.path(OUT_MAN, "legend_circular_5.tiff")
draw_legend(leg_36, width_cm = OUT_W_CM * 2)
log_msg("Legend (36 cm) saved: ", leg_36)

leg_72 <- file.path(OUT_MAN, "legend_circular_5_wide.tiff")
draw_legend(leg_72, width_cm = OUT_W_CM * 4)
log_msg("Legend (72 cm) saved: ", leg_72)

writeLines(capture.output(sessionInfo()),
           file.path(LOG_DIR, "step16ab_panel_assembly_sessionInfo.txt"))

sep_line()
log_msg("Step 16ab v5 finished -- plots saved to: ", OUT_MAN)
sep_line()

} # end SOURCED_16AB4 guard
