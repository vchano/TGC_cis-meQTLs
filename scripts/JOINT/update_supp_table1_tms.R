#!/usr/bin/env Rscript
############################################################
# Update Supplementary Table 1 — TMS + raw ECS note
#
# Changes applied to the _s2 xlsx (which already has the ECS
# trimmed read 602-sample update and integer SNP counts):
#
#   (a) Row 5 label: clarify raw ECS reads are n=18 only
#       (raw data was not retained for remaining samples)
#   (b) New TMS section appended at rows 17–20:
#       Row 17: section header for TMS
#       Row 18: TMS — Mean trimmed read pairs per sample, M (range) [n=N]
#       Row 19: TMS — Mean Bismark mapping efficiency, % (range) [n=N]
#       Row 20: TMS — updated abbreviation note (replaces row 15)
#
# Input:  _corrected_s8_s2.xlsx
# Output: _corrected_s8_s2_tms.xlsx
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(openxlsx2)
})

msg <- function(...) cat(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n"))

PROJECT_ROOT <- "/path/to/your/project"
NATGEN_DIR   <- file.path(PROJECT_ROOT, "RESULTS/DRAFT/NATURE.GENETICS")

XLSX_IN  <- file.path(NATGEN_DIR,
  "260819_Chano.etal.2026_tgc_supp.tables_corrected_s8_s2.xlsx")
XLSX_OUT <- file.path(NATGEN_DIR,
  "260819_Chano.etal.2026_tgc_supp.tables_corrected_s8_s2_tms.xlsx")

TSV <- file.path(PROJECT_ROOT,
  "RESULTS/JOINT/COMBINED5/tables/ecs_read_counts.tsv")

TBS_MAP <- file.path(PROJECT_ROOT, "DATA/TMS/MAPPED.FILES.TMS")

SHEET <- "Supplementary Table 1"

############################################################
# 1) LOAD ECS STUDY SAMPLE LIST (602 P-prefix samples)
############################################################

msg("Reading ECS sample list from read count TSV...")
rc <- fread(TSV, sep = "\t")
ecs_samples <- rc[substr(sample, 1, 1) == "P", sample]
msg("  ECS study samples: ", length(ecs_samples))
stopifnot(length(ecs_samples) == 602)

############################################################
# 2) PARSE BISMARK PE REPORTS FOR TMS STATS
############################################################

msg("Finding Bismark PE reports...")
report_files <- list.files(TBS_MAP, pattern = "bismark_bt2_PE_report\\.txt$",
                           recursive = TRUE, full.names = TRUE)
msg("  Found: ", length(report_files), " reports")

parse_bismark_report <- function(f) {
  lines <- readLines(f, warn = FALSE)
  # Sample name from filename: P001_WA02_R1_p_bismark_bt2_PE_report.txt -> P001_WA02
  sname <- sub("_R1_p_bismark_bt2_PE_report\\.txt$", "",
               basename(f))
  # Trimmed read pairs
  pairs_line <- grep("Sequence pairs analysed in total", lines, value = TRUE)
  pairs <- if (length(pairs_line)) {
    as.integer(trimws(sub(".*:\t", "", pairs_line[1])))
  } else NA_integer_

  # Mapping efficiency
  eff_line <- grep("Mapping efficiency:", lines, value = TRUE)
  eff <- if (length(eff_line)) {
    as.numeric(sub("%.*", "", trimws(sub("Mapping efficiency:\t", "", eff_line[1]))))
  } else NA_real_

  data.table(sample = sname, tbs_pairs = pairs, tbs_map_pct = eff)
}

msg("  Parsing reports...")
tms_dt <- rbindlist(lapply(report_files, parse_bismark_report))
msg("  Parsed: ", nrow(tms_dt))

# Filter to ECS study samples
tms_study <- tms_dt[sample %in% ecs_samples]
msg("  Matched to ECS study samples: ", nrow(tms_study))

n_tms <- nrow(tms_study)

# TMS trimmed read pairs
tr  <- tms_study$tbs_pairs
mn_tr <- mean(tr, na.rm = TRUE)
min_tr <- min(tr, na.rm = TRUE)
max_tr <- max(tr, na.rm = TRUE)

fmt_M <- function(x) sprintf("%.1f", x / 1e6)
tms_pairs_val <- sprintf("%s (%s–%s)",
                          fmt_M(mn_tr), fmt_M(min_tr), fmt_M(max_tr))

# TMS mapping efficiency
me  <- tms_study$tbs_map_pct
mn_me <- mean(me, na.rm = TRUE)
min_me <- min(me, na.rm = TRUE)
max_me <- max(me, na.rm = TRUE)
tms_map_val <- sprintf("%.1f (%.1f–%.1f)", mn_me, min_me, max_me)

msg("  TMS trimmed read pairs: mean=", round(mn_tr), " min=", round(min_tr),
    " max=", round(max_tr))
msg("  TMS mapping efficiency: mean=", round(mn_me, 1), "% min=", round(min_me, 1),
    "% max=", round(max_me, 1), "%")
msg("  Formatted pairs: ", tms_pairs_val)
msg("  Formatted mapping: ", tms_map_val)

############################################################
# 3) LOAD WORKBOOK AND APPLY CHANGES
############################################################

msg("Loading workbook: ", basename(XLSX_IN))
wb <- wb_load(XLSX_IN)

write_cell <- function(wb, sheet, row, col, val) {
  wb_add_data(wb, sheet = sheet,
              x = val, start_row = row, start_col = col,
              col_names = FALSE)
}

# --- Row 5 label: clarify raw read limitation ---
new_r5_label <- paste0(
  "Mean raw reads per sample (range) ",
  "[n = 18; raw data not retained after processing for remaining samples]")
msg("  Updating row 5 label")
wb <- write_cell(wb, SHEET, row = 5, col = 1, val = new_r5_label)

# --- New TMS section (rows 16–19) ---
# Row 15 is currently the abbreviation note; we write TMS after a blank at 16

# Row 16: blank separator (leave empty — just skip)

# Row 17: TMS section header (Step col only)
wb <- write_cell(wb, SHEET, row = 17, col = 1,
                 val = "Targeted methylation sequencing (TMS) data processing summary")
wb <- write_cell(wb, SHEET, row = 17, col = 2, val = "")

# Row 18: TMS trimmed read pairs
tms_pairs_label <- sprintf(
  "Mean trimmed read pairs per sample, M (range) [n = %d]", n_tms)
wb <- write_cell(wb, SHEET, row = 18, col = 1, val = tms_pairs_label)
wb <- write_cell(wb, SHEET, row = 18, col = 2, val = tms_pairs_val)

# Row 19: TMS mapping efficiency
tms_map_label <- sprintf(
  "Mean Bismark PE mapping efficiency, %% (range) [n = %d]", n_tms)
wb <- write_cell(wb, SHEET, row = 19, col = 1, val = tms_map_label)
wb <- write_cell(wb, SHEET, row = 19, col = 2, val = tms_map_val)

# Row 21: updated abbreviation note (after blank row 20)
wb <- write_cell(wb, SHEET, row = 21, col = 1,
  val = paste0("SNP: single nucleotide polymorphism; MAF: minimum allele frequency; ",
               "LD: linkage disequilibrium; M: millions of reads; ",
               "TMS: targeted methylation sequencing; PE: paired-end; ",
               "Bismark: bisulfite aligner v0.23.0"))

msg("  TMS rows written at rows 17–19 and abbreviation note at row 21")

############################################################
# 4) SAVE
############################################################

msg("Saving: ", basename(XLSX_OUT))
wb_save(wb, XLSX_OUT)
msg("Done: ", XLSX_OUT)
