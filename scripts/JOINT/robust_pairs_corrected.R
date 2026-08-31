#!/usr/bin/env Rscript
############################################################
# Robust meQTL pair correction
#
# Recomputes robust pairs (FDR < 1e-10 in BOTH GENESIS5 and MatrixEQTL5)
# with true genomic coordinates for SNP positions.
#
# Fixes:
#   - snp_pos: was GDS-internal sequential index from GENESIS5 output;
#     replaced with true bp position from snp_variant_annot.rds
#   - site_chr / site_pos: parsed from the site string identifier
#     (format: "PA_chrXX:start-end") — no separate annotation file needed
#   - No effect-direction filter (positive and negative beta both included)
#   - Matching by integer snp_id + site string (consistent across tools
#     since both used the same input GDS/methylation files)
#
# Outputs:
#   overlap/tables/robust_markers_breeding_corrected.tsv
#   overlap/tables/robust_markers_natural_corrected.tsv
#   overlap/tables/robust_context_summary_corrected.tsv
#   overlap/tables/overlap_summary_corrected.tsv
#   supplementary_information_table_s5.xlsx
#   LOGS/robust_pairs_corrected.log
#   LOGS/robust_pairs_corrected_report.md
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(openxlsx2)
})

options(stringsAsFactors = FALSE)

PROJECT_ROOT <- "/mnt/ceph-hdd/projects/scc_ufff_gailing/chano_TGC"

COHORTS  <- c("BREEDING", "NATURAL")
CONTEXTS <- c("CpG", "CHG", "CHH")
TOOLS    <- c("GENESIS5", "MATRIXEQTL5")

SIG_DIR   <- file.path(PROJECT_ROOT, "RESULTS/JOINT/COMBINED5/sig_sites")
ANNOT_DIR <- file.path(PROJECT_ROOT, "RESULTS/JOINT/MQTL5/INPUTS")
OUT_DIR   <- file.path(PROJECT_ROOT, "RESULTS/JOINT/COMBINED5/overlap/tables")
LOG_DIR   <- file.path(PROJECT_ROOT, "RESULTS/JOINT/COMBINED5/LOGS")
NATGEN_DIR <- file.path(PROJECT_ROOT, "RESULTS/DRAFT/NATURE.GENETICS")
REV1_DIR   <- file.path(PROJECT_ROOT, "RESULTS/DRAFT/NATURE.GENETICS/REV1")

dir.create(OUT_DIR,    recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR,    recursive = TRUE, showWarnings = FALSE)
dir.create(NATGEN_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(REV1_DIR,   recursive = TRUE, showWarnings = FALSE)

LOGFILE <- file.path(LOG_DIR, "robust_pairs_corrected.log")
if (file.exists(LOGFILE)) file.remove(LOGFILE)

msg <- function(...) {
  txt <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n")
  cat(txt)
  cat(txt, file = LOGFILE, append = TRUE)
}

sep <- function() msg(paste(rep("=", 65), collapse = ""))

sep()
msg("Robust meQTL pair correction with true genomic positions")
sep()

############################################################
# 1) LOAD SNP ANNOTATION MAPS
#    snp_variant_annot.rds : snp_id (char "1","2",...) -> chr, pos (true genomic)
#
#    NOTE: site positions are NOT needed from an annotation file.
#    The 'site' column in sig files is a string like "PA_chr03:790910176-790910176"
#    which already encodes the true genomic coordinates; site_chr and site_pos
#    are parsed directly from this string.
############################################################

msg("Loading SNP position annotation maps...")

snp_annot <- list()

for (cohort in COHORTS) {
  snp_annot[[cohort]] <- list()
  for (ctx in CONTEXTS) {
    va_path <- file.path(ANNOT_DIR, cohort, ctx, "snp_variant_annot.rds")
    if (file.exists(va_path)) {
      va <- as.data.table(readRDS(va_path))
      # snp_id is character ("1","2",...); cast to integer to match sig file snp column
      snp_annot[[cohort]][[ctx]] <- va[, .(snp_id = as.integer(snp_id),
                                           snp_chr = as.character(chr),
                                           snp_pos = as.integer(pos))]
      msg("  SNP annot: ", cohort, "/", ctx, " — ",
          nrow(snp_annot[[cohort]][[ctx]]), " SNPs")
    } else {
      msg("  MISSING SNP annot: ", va_path)
      snp_annot[[cohort]][[ctx]] <- data.table()
    }
  }
}

############################################################
# 2) LOAD SIG_P1E10 FILES
#    snp  : integer ID (consistent across tools — same GDS file)
#    site : character position string "PA_chrXX:start-end" (consistent across tools)
############################################################

sep()
msg("Loading sig_p1e10 files...")

sig <- list()
for (tool in TOOLS) {
  sig[[tool]] <- list()
  for (cohort in COHORTS) {
    fname <- file.path(SIG_DIR,
      sprintf("sig_p1e10_%s_%s.tsv", tolower(tool), tolower(cohort)))
    if (!file.exists(fname)) {
      msg("  MISSING: ", fname)
      sig[[tool]][[cohort]] <- data.table()
      next
    }
    dt <- fread(fname, showProgress = FALSE)
    dt[, snp  := as.integer(snp)]   # integer SNP ID
    # site stays as character: "PA_chr03:790910176-790910176"
    dt[, site := as.character(site)]
    sig[[tool]][[cohort]] <- dt
    msg("  ", tool, "/", cohort, ": ", nrow(dt), " pairs | contexts: ",
        paste(sort(unique(dt$context)), collapse = ", "))
  }
}

############################################################
# 3) PER-TOOL SUMMARY (before robust filter)
############################################################

sep()
msg("Computing per-tool per-context counts...")

tool_summary <- list()

for (tool in TOOLS) {
  for (cohort in COHORTS) {
    dt <- sig[[tool]][[cohort]]
    if (!nrow(dt)) next
    for (ctx in CONTEXTS) {
      sub <- dt[context == ctx]
      if (!nrow(sub)) next
      n_pairs <- nrow(sub)
      n_snps  <- uniqueN(sub$snp)
      n_sites <- uniqueN(sub$site)
      tool_summary[[length(tool_summary) + 1]] <- data.table(
        tool    = tool,
        cohort  = cohort,
        context = ctx,
        n_pairs = n_pairs,
        n_unique_snps  = n_snps,
        n_unique_sites = n_sites
      )
      msg("  ", tool, "/", cohort, "/", ctx, ": ",
          n_pairs, " pairs | ", n_snps, " SNPs | ", n_sites, " sites")
    }
  }
}

tool_summary_dt <- rbindlist(tool_summary)

############################################################
# 4) COMPUTE ROBUST PAIRS
#    Match on integer snp + character site string — both consistent
#    across GENESIS5 and MatrixEQTL5 (same GDS and methylation inputs).
#    Replace snp_pos with annotation-verified true genomic coordinates.
#    Parse site_chr and site_pos from the site string.
############################################################

sep()
msg("Computing robust pairs (same snp+site in BOTH tools)...")

robust_list    <- list()
robust_summary <- list()

for (cohort in COHORTS) {
  g5  <- sig[["GENESIS5"   ]][[cohort]]
  me5 <- sig[["MATRIXEQTL5"]][[cohort]]
  if (!nrow(g5) || !nrow(me5)) next

  for (ctx in CONTEXTS) {
    g5_ctx  <- g5 [context == ctx]
    me5_ctx <- me5[context == ctx]
    if (!nrow(g5_ctx) || !nrow(me5_ctx)) {
      msg("  ", cohort, "/", ctx, ": one tool has 0 pairs — skipping")
      next
    }

    # Select columns to keep from each tool
    keep_cols <- c("snp", "site", "statistic", "beta", "pvalue", "p_FDR")
    g5_keep  <- intersect(keep_cols, names(g5_ctx))
    me5_keep <- intersect(keep_cols, names(me5_ctx))

    g5_dt  <- g5_ctx [, ..g5_keep]
    me5_dt <- me5_ctx[, ..me5_keep]

    # Rename non-key columns to tool-suffix form
    g5_rename  <- setdiff(g5_keep,  c("snp", "site"))
    me5_rename <- setdiff(me5_keep, c("snp", "site"))
    setnames(g5_dt,  g5_rename,  paste0(g5_rename,  "_GENESIS5"))
    setnames(me5_dt, me5_rename, paste0(me5_rename, "_MATRIXEQTL5"))

    # Pair-level intersection: same snp AND same site in both tools
    both <- merge(g5_dt, me5_dt, by = c("snp", "site"))
    both[, context := ctx]
    both[, cohort  := cohort]

    # Add true SNP positions from snp_variant_annot.rds
    sa <- snp_annot[[cohort]][[ctx]]
    if (nrow(sa)) {
      both <- merge(both, sa, by.x = "snp", by.y = "snp_id", all.x = TRUE)
    } else {
      both[, `:=`(snp_chr = NA_character_, snp_pos = NA_integer_)]
    }

    # Parse site_chr and site_pos from site string "PA_chrXX:start-end"
    both[, site_chr := sub(":.*", "", site)]
    both[, site_pos := as.integer(sub("^[^:]+:(\\d+)-.*$", "\\1", site))]

    n_pairs <- nrow(both)
    n_snps  <- uniqueN(both$snp)
    n_sites <- uniqueN(both$site)

    msg("  ", cohort, "/", ctx, ": ", n_pairs, " robust pairs | ",
        n_snps, " unique SNPs | ", n_sites, " unique sites")

    na_snp  <- sum(is.na(both$snp_pos))
    na_site <- sum(is.na(both$site_pos))
    if (na_snp  > 0) msg("    WARNING: ", na_snp,  " pairs have NA snp_pos")
    if (na_site > 0) msg("    WARNING: ", na_site, " pairs have NA site_pos")

    robust_list[[paste0(cohort, "_", ctx)]] <- both
    robust_summary[[length(robust_summary) + 1]] <- data.table(
      cohort  = cohort, context = ctx,
      robust_pairs        = n_pairs,
      robust_unique_snps  = n_snps,
      robust_unique_sites = n_sites
    )
  }
}

robust_dt         <- rbindlist(robust_list, fill = TRUE)
robust_summary_dt <- rbindlist(robust_summary)

############################################################
# 5) COLUMN ORDER FOR EXPORT
############################################################

col_order <- c("cohort", "context",
               "snp", "snp_chr", "snp_pos",
               "site", "site_chr", "site_pos",
               grep("_GENESIS5$",    names(robust_dt), value = TRUE),
               grep("_MATRIXEQTL5$", names(robust_dt), value = TRUE))
col_order <- intersect(col_order, names(robust_dt))
setcolorder(robust_dt, col_order)

############################################################
# 6) WRITE CORRECTED ROBUST MARKER TABLES
############################################################

sep()
msg("Writing corrected robust marker tables...")

for (coh in COHORTS) {
  sub <- robust_dt[cohort == coh]
  out <- file.path(OUT_DIR, paste0("robust_markers_", tolower(coh), "_corrected.tsv"))
  fwrite(sub, out, sep = "\t")
  msg("  Saved: ", basename(out), " (", nrow(sub), " rows)")
}

fwrite(robust_summary_dt,
       file.path(OUT_DIR, "robust_context_summary_corrected.tsv"), sep = "\t")
fwrite(tool_summary_dt,
       file.path(OUT_DIR, "overlap_summary_corrected.tsv"), sep = "\t")
msg("  Summary tables saved.")

############################################################
# 7) BUILD TABLE 1 — COMBINED TOOL + ROBUST SUMMARY
############################################################

sep()
msg("Building Table 1...")

g5_wide  <- tool_summary_dt[tool == "GENESIS5",
  .(cohort, context,
    G5_pairs = n_pairs, G5_unique_snps = n_unique_snps, G5_unique_sites = n_unique_sites)]
me5_wide <- tool_summary_dt[tool == "MATRIXEQTL5",
  .(cohort, context,
    ME5_pairs = n_pairs, ME5_unique_snps = n_unique_snps, ME5_unique_sites = n_unique_sites)]
rob_wide <- robust_summary_dt[,
  .(cohort, context,
    Robust_pairs = robust_pairs, Robust_unique_snps = robust_unique_snps,
    Robust_unique_sites = robust_unique_sites)]

table1 <- Reduce(function(a, b) merge(a, b, by = c("cohort","context"), all = TRUE),
                 list(g5_wide, me5_wide, rob_wide))

ctx_order <- c("CpG","CHG","CHH")
coh_order <- c("BREEDING","NATURAL")
table1[, context := factor(context, levels = ctx_order)]
table1[, cohort  := factor(cohort,  levels = coh_order)]
setorder(table1, cohort, context)

msg("  Table 1 (", nrow(table1), " rows):")
print(table1)

############################################################
# 8) WRITE SUPPLEMENTARY TABLE S5 XLSX
############################################################

sep()
msg("Writing Supplementary Table S5 xlsx...")

S5_PATH      <- file.path(NATGEN_DIR, "supplementary_information_table_s5.xlsx")
S5_PATH_REV1 <- file.path(REV1_DIR,   "supplementary_information_table_s5.xlsx")

wb <- wb_workbook()
wb <- wb_add_worksheet(wb, "Supplementary Table S5")

wb <- wb_add_data(wb, sheet = "Supplementary Table S5",
                  x = "Supplementary Table S5. Robust cis-meQTL pairs (FDR < 1×10⁻¹⁰ in both GENESIS5 and MatrixEQTL5).",
                  start_row = 1, start_col = 1, col_names = FALSE)

note <- paste0(
  "snp_pos: true genomic bp position from snp_variant_annot.rds (MQTL5/INPUTS). ",
  "site_chr and site_pos: parsed from the site string identifier (format: chr:start-end). ",
  "Matching performed on integer snp_id and site string (consistent across tools; ",
  "both tools used the same GDS and methylation input files). ",
  "No effect-direction filter applied; both positive and negative beta values included."
)
wb <- wb_add_data(wb, sheet = "Supplementary Table S5",
                  x = note, start_row = 2, start_col = 1, col_names = FALSE)

s5_export <- copy(robust_dt)
setnames(s5_export,
  old = c("p_FDR_GENESIS5", "p_FDR_MATRIXEQTL5"),
  new = c("FDR_GENESIS5",   "FDR_MATRIXEQTL5"),
  skip_absent = TRUE)

wb <- wb_add_data(wb, sheet = "Supplementary Table S5",
                  x = as.data.frame(s5_export),
                  start_row = 3, start_col = 1, col_names = TRUE)

wb_save(wb, S5_PATH)
file.copy(S5_PATH, S5_PATH_REV1, overwrite = TRUE)
msg("  S5 saved: ", S5_PATH, " (", nrow(s5_export), " rows)")
msg("  S5 REV1:  ", S5_PATH_REV1)

############################################################
# 9) WRITE MD REPORT
############################################################

sep()
msg("Writing MD report...")

MD_PATH      <- file.path(LOG_DIR,  "robust_pairs_corrected_report.md")
MD_PATH_REV1 <- file.path(REV1_DIR, "robust_pairs_corrected_report.md")

# Pre-compute totals for MD
tot_b_pairs <- sum(table1[as.character(cohort) == "BREEDING"]$Robust_pairs,  na.rm = TRUE)
tot_n_pairs <- sum(table1[as.character(cohort) == "NATURAL" ]$Robust_pairs,  na.rm = TRUE)
tot_b_snps  <- uniqueN(robust_dt[cohort == "BREEDING"]$snp)
tot_n_snps  <- uniqueN(robust_dt[cohort == "NATURAL" ]$snp)
tot_b_sites <- uniqueN(robust_dt[cohort == "BREEDING"]$site)
tot_n_sites <- uniqueN(robust_dt[cohort == "NATURAL" ]$site)

tool_rows <- character(0)
for (i in seq_len(nrow(table1))) {
  r <- table1[i]
  tool_rows <- c(tool_rows, sprintf("| %s | %s | %s | %s | %s | %s | %s | %s |",
    r$cohort, r$context,
    r$G5_pairs,  r$G5_unique_snps,  r$G5_unique_sites,
    r$ME5_pairs, r$ME5_unique_snps, r$ME5_unique_sites))
}

rob_rows <- character(0)
for (i in seq_len(nrow(table1))) {
  r <- table1[i]
  rob_rows <- c(rob_rows, sprintf("| %s | %s | %s | %s | %s |",
    r$cohort, r$context,
    r$Robust_pairs, r$Robust_unique_snps, r$Robust_unique_sites))
}

  # Combined Table 1: all fields per cohort x context in one table
  t1_rows <- character(0)
  for (i in seq_len(nrow(table1))) {
    r <- table1[i]
    t1_rows <- c(t1_rows, sprintf(
      "| %s | %s | %d | %d | %d | %d | %d | %d | %d | %d | %d |",
      as.character(r$cohort), as.character(r$context),
      r$G5_pairs,  r$G5_unique_snps,  r$G5_unique_sites,
      r$ME5_pairs, r$ME5_unique_snps, r$ME5_unique_sites,
      r$Robust_pairs, r$Robust_unique_snps, r$Robust_unique_sites))
  }

md <- c(
  "# Robust meQTL pairs — Corrected analysis report",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M")),
  "",
  "## What was fixed",
  "",
  "- **snp_pos**: was GDS-internal sequential index from GENESIS5 output.",
  "  Replaced with true genomic bp position from `snp_variant_annot.rds`.",
  "- **site_chr / site_pos**: parsed directly from the `site` identifier string",
  "  (format: `PA_chrXX:start-end`). No separate annotation file needed.",
  "- **Matching**: integer `snp` ID + character `site` string — both are consistent",
  "  across GENESIS5 and MatrixEQTL5 (same GDS and methylation input files).",
  "- **No effect-direction filter**: positive and negative beta values included.",
  "- **Pair-level robust**: same snp+site pair must be significant (FDR < 1e-10)",
  "  in **both** GENESIS5 **and** MatrixEQTL5.",
  "",
  "---",
  "",
  "## Table 1. Summary of cis-meQTL pairs and robust pairs per cohort and methylation context",
  "",
  "Columns: pairs = SNP–site pairs at FDR < 1×10⁻¹⁰; unique SNPs = distinct SNP IDs in that set;",
  "unique sites = distinct methylation sites. Robust = same SNP+site pair significant in both tools.",
  "Unique SNPs in robust is identical for both tools (all robust SNPs are by definition in both).",
  "",
  "| Cohort | Context | GENESIS5 pairs | GENESIS5 unique SNPs | GENESIS5 unique sites | MatrixEQTL5 pairs | MatrixEQTL5 unique SNPs | MatrixEQTL5 unique sites | Robust pairs | Robust unique SNPs | Robust unique sites |",
  paste(rep("|---", 11), collapse = ""),
  t1_rows,
  "",
  paste0("**BREEDING total**: ", tot_b_pairs, " robust pairs | ",
         tot_b_snps, " unique SNPs | ", tot_b_sites, " unique sites"),
  "",
  paste0("**NATURAL total**: ", tot_n_pairs, " robust pairs | ",
         tot_n_snps, " unique SNPs | ", tot_n_sites, " unique sites"),
  "",
  "---",
  "",
  "## Output files",
  "",
  paste0("- `overlap/tables/robust_markers_breeding_corrected.tsv` (",
         nrow(robust_dt[cohort == "BREEDING"]), " rows)"),
  paste0("- `overlap/tables/robust_markers_natural_corrected.tsv` (",
         nrow(robust_dt[cohort == "NATURAL"]), " rows)"),
  "- `overlap/tables/robust_context_summary_corrected.tsv`",
  "- `overlap/tables/overlap_summary_corrected.tsv`",
  paste0("- `supplementary_information_table_s5.xlsx` (",
         nrow(s5_export), " robust pairs across all cohorts and contexts)"),
  "",
  "---",
  "",
  "## Next steps (pending agreement)",
  "",
  "1. Update `16ab5.R` Venn intersection to use pair-level positions → regenerate Fig. S4 and Fig. 4",
  "2. Replace `robust_markers_*.tsv` with corrected versions → rerun Step 17ab (annotation, S6)",
  "3. Update Table 1, Table 2, Table 3 counts in manuscript",
  "4. Update Abstract, Results, Methods, Discussion, Supplementary Note 5",
  "5. Grep manuscript DOCX for all old counts and replace systematically",
  "6. Push corrected scripts to GitHub",
  ""
)

writeLines(md, MD_PATH)
file.copy(MD_PATH, MD_PATH_REV1, overwrite = TRUE)
msg("  MD report: ", MD_PATH)
msg("  MD REV1:   ", MD_PATH_REV1)

sep()
msg("Done.")
sep()

writeLines(capture.output(sessionInfo()),
           file.path(LOG_DIR, "robust_pairs_corrected_sessionInfo.txt"))
