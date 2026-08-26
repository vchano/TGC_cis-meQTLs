#!/usr/bin/env Rscript
############################################################
# Update Supplementary Table 1 — S2 fix
#
# Changes:
#   (a) Row 6 label: "[n = 18 samples available]" → "[n = 602]"
#   (b) Row 6 value: recomputed from 602-sample trimmed read counts
#   (c) Rows 9–13 values: correct thousands-separator corruption
#       (488.826 → 488826, etc.)
#
# Input:  XLSX_IN  (corrected_s8 version as base)
#         TSV      (ecs_read_counts.tsv, 602 study samples)
# Output: XLSX_OUT (new _s2 suffix version)
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(openxlsx2)
})

msg <- function(...) cat(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n"))

PROJECT_ROOT <- "/path/to/your/project"
NATGEN_DIR   <- file.path(PROJECT_ROOT, "RESULTS/DRAFT/NATURE.GENETICS")

# Use the corrected_s8 xlsx as input base (preserves S8 fix)
XLSX_IN  <- file.path(NATGEN_DIR, "260819_Chano.etal.2026_tgc_supp.tables_corrected_s8.xlsx")
XLSX_OUT <- file.path(NATGEN_DIR, "260819_Chano.etal.2026_tgc_supp.tables_corrected_s8_s2.xlsx")

TSV <- file.path(PROJECT_ROOT, "RESULTS/JOINT/COMBINED5/tables/ecs_read_counts.tsv")

SHEET <- "Supplementary Table 1"

############################################################
# 1) COMPUTE 602-SAMPLE TRIMMED READ STATS
############################################################

msg("Reading read count TSV: ", basename(TSV))
rc <- fread(TSV, sep = "\t")
# Study samples only (P-prefix; x-prefix are excluded/failed)
rc_study <- rc[substr(sample, 1, 1) == "P"]
msg("  Study samples: ", nrow(rc_study))
stopifnot(nrow(rc_study) == 602)

# trimmed_reads column = read pairs counted from R1 FASTQ
tr <- rc_study$trimmed_reads
n_tr <- sum(!is.na(tr))
mn_tr <- mean(tr, na.rm = TRUE)
min_tr <- min(tr, na.rm = TRUE)
max_tr <- max(tr, na.rm = TRUE)

# Format as "X.X (Y.Y–Z.Z)" in millions, 1 decimal place
fmt_M <- function(x) sprintf("%.1f", x / 1e6)
trimmed_label_val <- sprintf("%s (%s–%s)",
                             fmt_M(mn_tr), fmt_M(min_tr), fmt_M(max_tr))

msg("  n = ", n_tr)
msg("  mean trimmed read pairs = ", round(mn_tr), " (", round(min_tr), "–", round(max_tr), ")")
msg("  formatted: ", trimmed_label_val)

############################################################
# 2) KNOWN CORRECT INTEGER VALUES FOR SNP COUNTS
#    (currently stored as European-locale decimals, e.g. 488.826 = 488,826)
############################################################

snp_vals <- list(
  row9  = 488826L,   # SNPs after quality + MAF > 0.05
  row10 = 450947L,   # SNPs after cohort-split BREEDING
  row11 = 462617L,   # SNPs after cohort-split NATURAL
  row12 = 120923L,   # SNPs after LD pruning BREEDING
  row13 = 132708L    # SNPs after LD pruning NATURAL
)

############################################################
# 3) LOAD WORKBOOK AND APPLY CHANGES
############################################################

msg("Loading workbook: ", basename(XLSX_IN))
wb <- wb_load(XLSX_IN)

# Helper: overwrite a single cell (1-indexed row, col within the sheet)
write_cell <- function(wb, sheet, row, col, val) {
  wb_add_data(wb, sheet = sheet,
              x = val, start_row = row, start_col = col,
              col_names = FALSE)
}

# --- Row 6: label (col A = col 1) ---
new_label_r6 <- "Mean trimmed reads per sample, M (range) [n = 602]"
msg("  Updating row 6 label -> ", new_label_r6)
wb <- write_cell(wb, SHEET, row = 6, col = 1, val = new_label_r6)

# --- Row 6: value (col B = col 2) ---
msg("  Updating row 6 value -> ", trimmed_label_val)
wb <- write_cell(wb, SHEET, row = 6, col = 2, val = trimmed_label_val)

# --- Rows 9–13: integer SNP counts ---
for (nm in names(snp_vals)) {
  rnum <- as.integer(sub("row", "", nm))
  val  <- snp_vals[[nm]]
  msg("  Updating row ", rnum, " value -> ", val)
  wb <- write_cell(wb, SHEET, row = rnum, col = 2, val = val)
}

############################################################
# 4) SAVE
############################################################

msg("Saving: ", basename(XLSX_OUT))
wb_save(wb, XLSX_OUT)
msg("Done: ", XLSX_OUT)
