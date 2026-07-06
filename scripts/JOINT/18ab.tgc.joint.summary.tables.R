#!/usr/bin/env Rscript
############################################################
# TreeGeneClimate (TGC) — JOINT ECS + TBS
# Step 20ab: Compile summary tables for manuscript
#
# OUTPUTS (RESULTS/JOINT/COMBINED5/tables/)
#   table_ecs_summary.tsv / .md   — ECS genotyping table
#   table_tbs_summary.tsv / .md   — TBS methylation table (replaces Table 1 + 1b)
#
# Manuscript draft:
#   Table 1 + Table 1b blocks replaced with new unified TBS table
#   ECS table inserted as Table S-ECS before Table 1
############################################################

options(stringsAsFactors = FALSE)

############################################################
# 1) PATHS
############################################################

# === USER CONFIGURATION ===
PROJECT_ROOT <- "/path/to/your/project"  # <-- set this
# ===========================
COMBINED5    <- file.path(PROJECT_ROOT, "RESULTS", "JOINT", "COMBINED5")
OUT_DIR      <- file.path(COMBINED5, "tables")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

DRAFT_MD   <- file.path(PROJECT_ROOT, "RESULTS", "DRAFT", "manuscript_draft.md")
ECS_DIR    <- file.path(PROJECT_ROOT, "RESULTS", "ECS", "VARIANT.CALLING")
# MultiQC general stats files from WES raw and trimmed QC runs
MQC_RAW    <- file.path(PROJECT_ROOT, "ECO_WES", "1.QC",  "RAWDATA",  "multiqc_data", "multiqc_general_stats.txt")
MQC_TRIM   <- file.path(PROJECT_ROOT, "ECO_WES", "2.TRIMMED",          "multiqc_data", "multiqc_general_stats.txt")
# bcftools stats output files for raw and MAF-filtered VCF
STATS_DIR  <- file.path(ECS_DIR, "STATS_ALL")
# PLINK BIM files for LD-pruned SNP sets (one per cohort)
BIM_B      <- file.path(PROJECT_ROOT, "RESULTS", "ECS", "POPGEN", "STRUCTURE",
                        "BREEDING", "tgc.ecs.breeding.admix.pruned.bim")
BIM_N      <- file.path(PROJECT_ROOT, "RESULTS", "ECS", "POPGEN", "STRUCTURE",
                        "NATURAL",  "tgc.ecs.natural.admix.pruned.bim")

# methylKit summary CSVs from the TBS input-prep step
TBS_MK_DIR <- file.path(PROJECT_ROOT, "RESULTS", "TBS", "RANALYSIS", "METHYLKIT_OBJECTS")
# ANOVA log from step 6b (methylation level comparisons between populations)
ANOVA_DIR  <- file.path(PROJECT_ROOT, "RESULTS", "TBS", "RANALYSIS", "ANOVA.METHYL.LEVEL")
# Kruskal-Wallis SVMP result files from step 8b
KW_DIR     <- file.path(PROJECT_ROOT, "RESULTS", "TBS", "RANALYSIS", "TABLES", "heatmap_markers_8B")
# Robust cis-meQTL summary produced by 15ab.R
ROBUST_TSV <- file.path(COMBINED5, "overlap", "tables", "robust_context_summary.tsv")

msg <- function(...) cat(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n"))

############################################################
# 2) ECS TABLE
############################################################

msg("=== ECS table ===")

## 2a) MultiQC read counts — R1 rows only to avoid double-counting paired reads
parse_mqc_r1 <- function(path) {
  if (!file.exists(path)) { msg("  MISSING: ", basename(path)); return(NULL) }
  d   <- read.table(path, sep = "\t", header = TRUE, check.names = FALSE)
  col <- "FastQC_mqc-generalstats-fastqc-total_sequences"
  if (!col %in% names(d)) { msg("  Column missing: ", col); return(NULL) }
  # Keep only R1 reads to avoid double-counting paired-end reads
  r1  <- d[grepl("_R1", d$Sample, fixed = TRUE), ]
  if (nrow(r1) == 0) r1 <- d
  setNames(as.numeric(r1[[col]]), r1$Sample)
}

raw_reads  <- parse_mqc_r1(MQC_RAW)
trim_reads <- parse_mqc_r1(MQC_TRIM)

# Format as "mean (min–max)" in millions, suitable for a table cell
fmt_reads <- function(v) {
  if (is.null(v) || length(v) == 0) return("N/A")
  m <- v / 1e6
  sprintf("%.1f (%.1f–%.1f)", mean(m), min(m), max(m))
}

# Compute overall retention rate from total read sums (avoids sample-name mismatch
# between raw and trimmed MultiQC files — sample names may differ slightly)
pct_ret <- if (!is.null(raw_reads) && !is.null(trim_reads) &&
               length(raw_reads) > 0 && length(trim_reads) > 0) {
  # Names differ between raw/trimmed MultiQC; use aggregate sum for retention rate
  sprintf("%.1f", sum(trim_reads) / sum(raw_reads) * 100)
} else "N/A"

n_mqc <- if (!is.null(raw_reads)) length(raw_reads) else 0L

## 2b) bcftools stats SNP counts
parse_bcf_snps <- function(path) {
  if (!file.exists(path)) { msg("  MISSING: ", basename(path)); return(NA_integer_) }
  ln <- readLines(path)
  # SN lines in bcftools stats have format: SN\t0\tnumber of SNPs:\t<n>
  sn <- grep("^SN\t0\tnumber of SNPs:", ln, value = TRUE)
  if (length(sn) == 0) return(NA_integer_)
  as.integer(trimws(sub(".*SNPs:\t", "", sn[1])))
}

snp_raw      <- parse_bcf_snps(file.path(STATS_DIR, "unfiltered.bcftools.stats.txt"))
snp_filtered <- parse_bcf_snps(file.path(STATS_DIR, "filtered_maf05.bcftools.stats.txt"))

## 2c) LD-pruned BIM row counts
count_bim <- function(path) {
  if (!file.exists(path)) { msg("  MISSING: ", basename(path)); return(NA_integer_) }
  # Each BIM row is one variant; wc -l gives the count directly
  as.integer(trimws(system(paste("wc -l <", shQuote(path)), intern = TRUE)))
}

snp_ld_b <- count_bim(BIM_B)
snp_ld_n <- count_bim(BIM_N)

## 2d) Format helper — integers with thousands separator for readability
fnum <- function(x) {
  if (is.na(x) || is.null(x)) return("N/A")
  format(as.integer(x), big.mark = ",", scientific = FALSE)
}

# Note if fewer than all 602 samples appear in MultiQC (e.g. partial run)
mqc_note <- if (n_mqc < 600) paste0(" [n = ", n_mqc, " samples available]") else ""

ecs_tbl <- data.frame(
  Step = c(
    "Samples (total)",
    paste0("Mean raw reads per sample, M (range)", mqc_note),
    paste0("Mean trimmed reads per sample, M (range)", mqc_note),
    "Reads retained after trimming (%)",
    "SNPs called (unfiltered)",
    "SNPs after quality + MAF > 0.05 filter",
    "SNPs after LD pruning — Breeding cohort (n = 209)",
    "SNPs after LD pruning — Natural cohort (n = 393)"
  ),
  Value = c(
    "602",
    fmt_reads(raw_reads),
    fmt_reads(trim_reads),
    pct_ret,
    fnum(snp_raw),
    fnum(snp_filtered),
    fnum(snp_ld_b),
    fnum(snp_ld_n)
  ),
  stringsAsFactors = FALSE
)

write.table(ecs_tbl, file.path(OUT_DIR, "table_ecs_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
msg("  Saved: table_ecs_summary.tsv")

############################################################
# 3) TBS TABLE
############################################################

msg("=== TBS table ===")

## 3a) methylKit summary CSVs
# Files named: summary_<cohort>_cov5_50_mpg*_mef0.05.csv — one per cohort × context
mk_files <- Sys.glob(file.path(TBS_MK_DIR, "summary_*_cov5_50_mpg*_mef0.05.csv"))
if (length(mk_files) == 0) stop("No methylKit summary CSVs found in: ", TBS_MK_DIR)
mk <- do.call(rbind, lapply(mk_files, read.csv))
# Normalise context to canonical case (CpG, CHG, CHH)
ctx_norm   <- function(x) ifelse(tolower(x) == "cpg", "CpG", toupper(x))
mk$context <- ctx_norm(mk$context)
# Title-case cohort names to match the rest of the table
mk$cohort  <- paste0(toupper(substring(tolower(mk$cohort), 1L, 1L)),
                     tolower(substring(mk$cohort, 2L)))

## 3b) Parse ANOVA log
anova_log <- readLines(file.path(ANOVA_DIR, "Step6b_methylation_level_stats_posthoc.log"))

# Extract per-cohort×context blocks from the log file.
# Each block starts with a "COHORT / context" header line and continues
# until the next header or end-of-file.
parse_anova_log <- function(lines) {
  res <- list()
  n   <- length(lines)
  i   <- 1L
  while (i <= n) {
    m <- regmatches(lines[i],
                    regexpr("^(BREEDING|NATURAL) / (CpG|CHG|CHH)$", lines[i]))
    if (length(m) == 1L) {
      parts  <- strsplit(m, " / ")[[1L]]
      cohort <- paste0(toupper(substring(tolower(parts[1L]), 1L, 1L)),
                       tolower(substring(parts[1L], 2L)))
      ctx    <- parts[2L]
      blk    <- list(cohort = cohort, context = ctx)
      j      <- i + 1L
      while (j <= n && !grepl("^(BREEDING|NATURAL) /", lines[j])) {
        if (grepl("Method:",          lines[j]))
          blk$method    <- trimws(sub(".*Method:", "",          lines[j]))
        if (grepl("Shapiro-Wilk p",   lines[j]))
          blk$sw_p      <- as.numeric(trimws(sub(".*Shapiro-Wilk p.*:", "", lines[j])))
        if (grepl("Levene p",         lines[j]))
          blk$levene_p  <- as.numeric(trimws(sub(".*Levene p.*:",      "", lines[j])))
        if (grepl("ANOVA global p:",  lines[j]))
          blk$anova_p   <- as.numeric(trimws(sub(".*ANOVA global p:",  "", lines[j])))
        if (grepl("Posthoc TSV:",     lines[j]))
          blk$posthoc   <- trimws(sub(".*Posthoc TSV:", "", lines[j]))
        j <- j + 1L
      }
      res[[paste0(cohort, "_", ctx)]] <- blk
      i <- j
    } else {
      i <- i + 1L
    }
  }
  res
}

anova_res <- parse_anova_log(anova_log)
msg(sprintf("  Parsed ANOVA log: %d blocks", length(anova_res)))

## 3c) Count significant TukeyHSD pairs
count_tukey <- function(tsv) {
  if (is.null(tsv) || !nzchar(tsv) || !file.exists(tsv)) return(NA_character_)
  d <- read.table(tsv, sep = "\t", header = TRUE, check.names = TRUE)
  # column is "p adj" → R reads as "p.adj"
  padj_col <- grep("p.adj|p_adj", names(d), ignore.case = TRUE, value = TRUE)[1L]
  if (is.na(padj_col)) return(NA_character_)
  total <- nrow(d)
  sig   <- sum(d[[padj_col]] < 0.05, na.rm = TRUE)
  # Report as "significant / total" pairs for transparency
  paste0(sig, " / ", total)
}

## 3d) Count significant KW markers (q < 0.05)
count_svmps <- function(cohort, context) {
  path <- file.path(KW_DIR,
    sprintf("TBS_8B_locus_tests_%s_%s.tsv", tolower(cohort), tolower(context)))
  if (!file.exists(path)) { msg("  MISSING KW: ", basename(path)); return(NA_integer_) }
  d <- read.table(path, sep = "\t", header = TRUE)
  # Count loci whose BH-adjusted KW p-value is below the significance threshold
  sum(d$padj < 0.05, na.rm = TRUE)
}

## 3e) Robust cis-meQTL sites
# From 15ab.R: pair-level intersection of GENESIS5 and MatrixEQTL5 results
robust <- read.table(ROBUST_TSV, sep = "\t", header = TRUE)
robust$cohort <- paste0(toupper(substring(tolower(robust$cohort), 1L, 1L)),
                        tolower(substring(robust$cohort, 2L)))

## 3f) Assemble per-cohort×context rows
COHORTS  <- c("Breeding", "Natural")
CONTEXTS <- c("CpG", "CHG", "CHH")
COMBOS   <- expand.grid(cohort = COHORTS, context = CONTEXTS,
                        stringsAsFactors = FALSE)[, c("cohort", "context")]

tbs_rows <- lapply(seq_len(nrow(COMBOS)), function(i) {
  coh <- COMBOS$cohort[i]
  ctx <- COMBOS$context[i]
  key <- paste0(coh, "_", ctx)

  mk_row <- mk[mk$cohort == coh & mk$context == ctx, , drop = FALSE]
  an     <- anova_res[[key]]
  rob    <- robust[robust$cohort == coh & robust$context == ctx, , drop = FALSE]

  # Safe accessor — returns NA if column absent (handles missing methylKit fields)
  get_mk <- function(col) if (nrow(mk_row) > 0 && col %in% names(mk_row))
                              mk_row[[col]][1L] else NA

  fmt_p <- function(x) if (is.null(x) || is.na(x)) "N/A"
              else formatC(x, format = "e", digits = 2)

  data.frame(
    Cohort        = coh,
    Context       = ctx,
    N             = get_mk("n_samples"),
    Raw_sites_med = get_mk("raw_sites_median"),      # median cytosine sites before coverage filter
    Cov_sites_med = get_mk("cov_sites_median"),      # median sites after 5–50x coverage filter
    Sites_unite   = get_mk("n_sites_after_unite"),   # sites in ≥min.per.group samples (unite())
    Sites_MEF     = get_mk("n_sites_after_mef"),     # sites with methylation effect fraction > 5%
    SW_p          = fmt_p(if (!is.null(an)) an$sw_p    else NA),
    Levene_p      = fmt_p(if (!is.null(an)) an$levene_p else NA),
    ANOVA_p       = fmt_p(if (!is.null(an)) an$anova_p  else NA),
    TukeyHSD_sig  = count_tukey(if (!is.null(an)) an$posthoc else NULL),
    SVMPs_q005    = count_svmps(coh, ctx),
    Robust_sites  = if (nrow(rob) > 0) rob$robust_unique_sites[1L] else NA_integer_,
    stringsAsFactors = FALSE
  )
})

tbs_tbl <- do.call(rbind, tbs_rows)

write.table(tbs_tbl, file.path(OUT_DIR, "table_tbs_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
msg("  Saved: table_tbs_summary.tsv")

############################################################
# 4) FORMAT MARKDOWN TABLES
############################################################

# Integer formatter with thousands separator for markdown table cells
fmt_int <- function(x) {
  v <- suppressWarnings(as.integer(x))
  if (is.na(v)) "N/A" else format(v, big.mark = ",", scientific = FALSE)
}

## ECS markdown
ecs_md_lines <- c(
  "**Table S-ECS. Exome capture sequencing (ECS) data processing summary.**",
  "",
  "| Step | Value |",
  "|:-----|:------|",
  sprintf("| %s | %s |", ecs_tbl$Step, ecs_tbl$Value),
  ""
)
ecs_md <- paste(ecs_md_lines, collapse = "\n")

## TBS unified Table 1
tbs_header <- paste0(
  "| Cohort | Context | N | Raw sites^a^ | Cov. sites^b^ |",
  " After unite^c^ | MEF sites^d^ |",
  " SW p | Levene p | ANOVA p | Sig. pairs^e^ | SVMPs^f^ | Robust sites^g^ |"
)
tbs_sep <- paste0(
  "|--------|---------|--:|-------------:|--------------:|",
  "---------------:|-------------:|",
  "------:|----------:|--------:|:-------------|--------:|-----------------:|"
)

# Construct one markdown row per cohort × context combination
tbs_body <- apply(tbs_tbl, 1L, function(r) {
  paste0(
    "| ", r["Cohort"],
    " | ", r["Context"],
    " | ", fmt_int(r["N"]),
    " | ", fmt_int(r["Raw_sites_med"]),
    " | ", fmt_int(r["Cov_sites_med"]),
    " | ", fmt_int(r["Sites_unite"]),
    " | ", fmt_int(r["Sites_MEF"]),
    " | ", r["SW_p"],
    " | ", r["Levene_p"],
    " | ", r["ANOVA_p"],
    " | ", if (is.na(r["TukeyHSD_sig"])) "N/A" else r["TukeyHSD_sig"],
    " | ", fmt_int(r["SVMPs_q005"]),
    " | ", fmt_int(r["Robust_sites"]),
    " |"
  )
})

# Footnotes explain abbreviations and filtering criteria used for each column
tbs_footnotes <- paste(
  "^a^ Median cytosine sites per sample before coverage filtering.",
  "^b^ Median sites per sample after 5–50× coverage filter.",
  "^c^ Sites present in at least min.per.group samples (after `unite()`).",
  "^d^ Sites with methylation effect fraction (MEF) > 5%.",
  "^e^ Significant TukeyHSD pairwise comparisons (p adj < 0.05) / total pairs.",
  "^f^ Loci with significant Kruskal-Wallis test (BH-adjusted q < 0.05).",
  "^g^ Robust cis-meQTL sites identified by both GENESIS and MatrixEQTL.",
  sep = "  \n"
)

tbs_md_lines <- c(
  "**Table 1. Targeted bisulfite sequencing (TBS) data summary by cohort and cytosine context.**",
  "",
  tbs_header,
  tbs_sep,
  tbs_body,
  "",
  tbs_footnotes,
  ""
)
tbs_md <- paste(tbs_md_lines, collapse = "\n")

writeLines(ecs_md, file.path(OUT_DIR, "table_ecs_summary.md"))
writeLines(tbs_md, file.path(OUT_DIR, "table_tbs_summary.md"))
msg("  Saved: table_ecs_summary.md, table_tbs_summary.md")

############################################################
# 5) UPDATE MANUSCRIPT DRAFT
############################################################

msg("=== Updating manuscript draft ===")

if (!file.exists(DRAFT_MD)) {
  msg("  WARNING: manuscript_draft.md not found — skipping")
} else {
  draft <- readLines(DRAFT_MD)

  # Locate start of Table 1 block (first "**Table 1" line)
  t1_start <- grep("^\\*\\*Table 1[. ]", draft)[1L]
  # Locate Table 1b start (separate sub-table that may follow)
  t1b_start <- grep("^\\*\\*Table 1b\\.", draft)[1L]

  if (is.na(t1_start)) {
    msg("  WARNING: Could not find '**Table 1' in draft — appending")
    writeLines(c(draft, "", tbs_md_lines), DRAFT_MD)
  } else {
    # End of the block = last "|"-prefixed line at or after the last sub-table start
    block_end_search_from <- if (!is.na(t1b_start)) t1b_start else t1_start
    after <- seq(block_end_search_from, length(draft))
    pipe_rows <- after[grepl("^\\|", draft[after])]
    t1_end <- if (length(pipe_rows) > 0) max(pipe_rows) else t1_start

    msg(sprintf("  Replacing Table 1 block (lines %d–%d)", t1_start, t1_end))

    new_draft <- c(
      draft[seq_len(t1_start - 1L)],
      tbs_md_lines,
      draft[seq(t1_end + 1L, length(draft))]
    )
    writeLines(new_draft, DRAFT_MD)
    msg("  Table 1 + Table 1b replaced with unified TBS table")

    # Insert ECS table just before "**Table 1."
    draft2   <- readLines(DRAFT_MD)
    t1_pos2  <- grep("^\\*\\*Table 1\\.", draft2)[1L]
    if (!is.na(t1_pos2) && t1_pos2 > 1L) {
      new_draft2 <- c(
        draft2[seq_len(t1_pos2 - 1L)],
        ecs_md_lines,
        "",
        draft2[seq(t1_pos2, length(draft2))]
      )
      writeLines(new_draft2, DRAFT_MD)
      msg("  ECS table (Table S-ECS) inserted before Table 1")
    }
  }
}

msg("Step 20ab finished.")
msg("  Outputs in: ", OUT_DIR)

sessionInfo()
