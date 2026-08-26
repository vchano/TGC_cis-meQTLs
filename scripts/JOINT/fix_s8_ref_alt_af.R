#!/usr/bin/env Rscript
############################################################
# Fix S8: populate missing ref/alt/af in Supplementary Table 6
#
# Strategy:
#   1. Read Supplementary Table 6 (headers on row 3), preserving original order
#   2. Collect unique (chr, pos, cohort) where ref is NA
#   3. Query per-cohort imputed VCF with bcftools (has AF)
#   4. Fallback to raw (non-imputed) VCF for anything still NA
#   5. Join results back, preserving original row order
#   6. Write corrected xlsx: load original wb, overwrite only the
#      ref / alt / af columns in Sheet 6 (preserves all other formatting)
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(openxlsx2)
})

msg <- function(...) cat(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n"))

PROJECT_ROOT <- "/path/to/your/project"
NATGEN_DIR   <- file.path(PROJECT_ROOT, "RESULTS/DRAFT/NATURE.GENETICS")

XLSX_IN  <- file.path(NATGEN_DIR, "260819_Chano.etal.2026_tgc_supp.tables.xlsx")
XLSX_OUT <- file.path(NATGEN_DIR, "260819_Chano.etal.2026_tgc_supp.tables_corrected_s8.xlsx")

VCF <- list(
  BREEDING = list(
    imputed = file.path(PROJECT_ROOT,
      "RESULTS/ECS/VCF_SPLIT/tgc.ecs.breeding.call.filt.maf05.snvs.poly.imputed.vcf.gz"),
    raw     = file.path(PROJECT_ROOT,
      "RESULTS/ECS/VCF_SPLIT/tgc.ecs.breeding.call.filt.maf05.snvs.poly.vcf.gz")
  ),
  NATURAL = list(
    imputed = file.path(PROJECT_ROOT,
      "RESULTS/ECS/VCF_SPLIT/tgc.ecs.natural.call.filt.maf05.snvs.poly.imputed.vcf.gz"),
    raw     = file.path(PROJECT_ROOT,
      "RESULTS/ECS/VCF_SPLIT/tgc.ecs.natural.call.filt.maf05.snvs.poly.vcf.gz")
  )
)

# Unfiltered all-samples VCF: used as fallback for positions filtered out per cohort
# (e.g., MAF < 0.05 in the cohort-specific subset). AF is computed on the fly
# per cohort using bcftools fill-tags with cohort-specific sample lists.
VCF_UNFILT  <- file.path(PROJECT_ROOT,
  "RESULTS/ECS/VARIANT.CALLING/tgc.ecs.allsamples.call.unfilt.snvs.renamed.vcf.gz")
SAMPS_B <- file.path(PROJECT_ROOT,
  "RESULTS/ECS/VARIANT.CALLING2/SPLIT/all_breeding_samples.txt")
SAMPS_N <- file.path(PROJECT_ROOT,
  "RESULTS/ECS/VARIANT.CALLING2/SPLIT/all_natural_samples.txt")
COHORT_SAMPS <- list(BREEDING = SAMPS_B, NATURAL = SAMPS_N)

TMP_DIR <- file.path(PROJECT_ROOT, "RESULTS/JOINT/COMBINED5/tmp_s8fix")
dir.create(TMP_DIR, recursive = TRUE, showWarnings = FALSE)

BCFTOOLS <- Sys.which("bcftools")
if (!nzchar(BCFTOOLS)) stop("bcftools not found — load module bcftools")
msg("bcftools: ", BCFTOOLS)

############################################################
# 1) INDEX VCF FILES (if TBI absent)
############################################################

msg("Checking/creating VCF indices (CSI — required for large Norway spruce chromosomes >512Mb)...")
for (cohort in names(VCF)) {
  for (type in names(VCF[[cohort]])) {
    vcf <- VCF[[cohort]][[type]]
    csi <- paste0(vcf, ".csi")
    if (!file.exists(csi)) {
      msg("  Indexing CSI (", cohort, " ", type, "): ", basename(vcf))
      ret <- system(paste(shQuote(BCFTOOLS), "index -c", shQuote(vcf)))
      if (ret != 0) msg("  WARNING: index failed for ", basename(vcf))
    } else {
      msg("  Index OK (CSI): ", basename(vcf))
    }
  }
}

############################################################
# 2) READ SUPPLEMENTARY TABLE 6  (headers on row 3, preserving row order)
############################################################

msg("Reading Supplementary Table 6...")
st6 <- as.data.table(read_xlsx(XLSX_IN, sheet = "Supplementary Table 6",
                                start_row = 3))
# Drop trailing all-NA column if present
st6 <- st6[, names(st6)[!is.na(names(st6)) & names(st6) != "NA_"], with = FALSE]

# Preserve original row order
st6[, .row_idx := .I]
st6[, pos := as.integer(pos)]

msg("  Rows: ", nrow(st6), " | NA ref: ", sum(is.na(st6$ref)))

# Column positions in the xlsx (1-based): row 3 = header, data from row 4
# gene_id=1, snp=2, chr=3, pos=4, marker_id=5, context=6, gene_start=7,
# gene_end=8, gene_strand=9, distance_bp=10, ref=11, alt=12, af=13, ...
COL_REF <- which(names(st6) == "ref")
COL_ALT <- which(names(st6) == "alt")
COL_AF  <- which(names(st6) == "af")
msg("  ref col=", COL_REF, " alt col=", COL_ALT, " af col=", COL_AF)

############################################################
# 3) VCF QUERY HELPERS
############################################################

# Query a filtered VCF that already has AF in INFO
query_vcf <- function(markers_dt, vcf_file, tmp_prefix, has_af = TRUE) {
  empty <- data.table(chr = character(), pos = integer(),
                      ref = character(), alt = character(), af = numeric())
  if (!file.exists(vcf_file)) {
    msg("  VCF not found: ", basename(vcf_file)); return(empty)
  }
  markers_dt <- unique(markers_dt[!is.na(chr) & !is.na(pos), .(chr, pos)])
  if (!nrow(markers_dt)) return(empty)

  reg_file <- paste0(tmp_prefix, ".bed")
  reg <- markers_dt[, .(chr, start = pos - 1L, end = pos)]
  setorder(reg, chr, start)
  fwrite(reg, reg_file, sep = "\t", col.names = FALSE)

  out_file <- paste0(tmp_prefix, ".tsv")
  fmt <- if (has_af) "'%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/AF\\n'"
         else        "'%CHROM\\t%POS\\t%REF\\t%ALT\\t.\\n'"
  cmd <- sprintf("%s query -R %s -f %s %s > %s 2>/dev/null",
                 shQuote(BCFTOOLS), shQuote(reg_file), fmt,
                 shQuote(vcf_file), shQuote(out_file))
  system(cmd)
  file.remove(reg_file)

  if (!file.exists(out_file) || file.info(out_file)$size == 0) return(empty)
  res <- fread(out_file, header = FALSE, sep = "\t", fill = TRUE,
               col.names = c("chr","pos","ref","alt","af"), na.strings = ".")
  file.remove(out_file)
  res[, pos := as.integer(pos)]
  res[, af  := suppressWarnings(as.numeric(af))]
  unique(res, by = c("chr","pos"))
}

# Query the unfiltered all-samples VCF, subsetting to cohort samples,
# and compute AF on the fly with bcftools fill-tags.
# Used for positions filtered out of the per-cohort filtered VCFs.
query_vcf_unfilt <- function(markers_dt, cohort, tmp_prefix) {
  empty <- data.table(chr = character(), pos = integer(),
                      ref = character(), alt = character(), af = numeric())
  if (!file.exists(VCF_UNFILT)) {
    msg("  Unfiltered VCF not found: ", basename(VCF_UNFILT)); return(empty)
  }
  samps_file <- COHORT_SAMPS[[cohort]]
  if (!file.exists(samps_file)) {
    msg("  Sample list not found: ", basename(samps_file)); return(empty)
  }
  markers_dt <- unique(markers_dt[!is.na(chr) & !is.na(pos), .(chr, pos)])
  if (!nrow(markers_dt)) return(empty)

  reg_file <- paste0(tmp_prefix, "_unfilt.bed")
  reg <- markers_dt[, .(chr, start = pos - 1L, end = pos)]
  setorder(reg, chr, start)
  fwrite(reg, reg_file, sep = "\t", col.names = FALSE)

  out_file <- paste0(tmp_prefix, "_unfilt.tsv")
  # Subset to cohort samples, compute AF via fill-tags, then extract fields
  cmd <- sprintf(
    "%s view --samples-file %s -R %s %s | %s plugin fill-tags -- -t AF | %s query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT\\t%%INFO/AF\\n' > %s 2>/dev/null",
    shQuote(BCFTOOLS), shQuote(samps_file), shQuote(reg_file),
    shQuote(VCF_UNFILT),
    shQuote(BCFTOOLS), shQuote(BCFTOOLS),
    shQuote(out_file)
  )
  system(cmd)
  file.remove(reg_file)

  if (!file.exists(out_file) || file.info(out_file)$size == 0) return(empty)
  res <- fread(out_file, header = FALSE, sep = "\t", fill = TRUE,
               col.names = c("chr","pos","ref","alt","af"), na.strings = ".")
  file.remove(out_file)
  res[, pos := as.integer(pos)]
  res[, af  := suppressWarnings(as.numeric(af))]
  unique(res, by = c("chr","pos"))
}

############################################################
# 4) FILL MISSING VALUES PER COHORT
############################################################

fill_cohort <- function(sub, cohort) {
  missing <- unique(sub[is.na(ref), .(chr, pos)])
  if (!nrow(missing)) { msg("  No missing in ", cohort); return(sub) }
  msg("  ", cohort, ": ", nrow(missing), " unique chr:pos to query")

  # Step 1: Imputed per-cohort VCF (has AF, fastest path)
  vi <- query_vcf(missing, VCF[[cohort]]$imputed,
                  file.path(TMP_DIR, paste0(tolower(cohort), "_imp")),
                  has_af = TRUE)
  msg("  Imputed VCF: ", nrow(vi), " hits")
  if (nrow(vi)) {
    sub <- merge(sub, vi[, .(chr, pos, ref_q=ref, alt_q=alt, af_q=af)],
                 by=c("chr","pos"), all.x=TRUE, sort=FALSE)
    sub[is.na(ref) & !is.na(ref_q), `:=`(ref=ref_q, alt=alt_q)]
    sub[is.na(af)  & !is.na(af_q),  af := af_q]
    sub[, c("ref_q","alt_q","af_q") := NULL]
  }

  # Step 2: Raw per-cohort VCF (no AF INFO field — positions that aren't imputed)
  still <- unique(sub[is.na(ref), .(chr, pos)])
  if (nrow(still)) {
    msg("  Still missing: ", nrow(still), " — trying raw per-cohort VCF")
    vr <- query_vcf(still, VCF[[cohort]]$raw,
                    file.path(TMP_DIR, paste0(tolower(cohort), "_raw")),
                    has_af = FALSE)
    msg("  Raw per-cohort VCF: ", nrow(vr), " hits")
    if (nrow(vr)) {
      sub <- merge(sub, vr[, .(chr, pos, ref_q=ref, alt_q=alt)],
                   by=c("chr","pos"), all.x=TRUE, sort=FALSE)
      sub[is.na(ref) & !is.na(ref_q), `:=`(ref=ref_q, alt=alt_q)]
      sub[, c("ref_q","alt_q") := NULL]
    }
  }

  # Step 3: Unfiltered all-samples VCF (fallback for positions that were filtered
  # out of per-cohort VCFs due to MAF < 0.05 or other filters). AF is computed
  # per cohort on the fly via bcftools fill-tags with cohort-specific sample lists.
  still2 <- unique(sub[is.na(ref), .(chr, pos)])
  if (nrow(still2)) {
    msg("  Still missing: ", nrow(still2), " — trying unfiltered VCF with fill-tags (", cohort, " samples)")
    vu <- query_vcf_unfilt(still2, cohort,
                           file.path(TMP_DIR, paste0(tolower(cohort), "_unfilt")))
    msg("  Unfiltered VCF: ", nrow(vu), " hits")
    if (nrow(vu)) {
      sub <- merge(sub, vu[, .(chr, pos, ref_q=ref, alt_q=alt, af_q=af)],
                   by=c("chr","pos"), all.x=TRUE, sort=FALSE)
      sub[is.na(ref) & !is.na(ref_q), `:=`(ref=ref_q, alt=alt_q)]
      sub[is.na(af)  & !is.na(af_q),  af := af_q]
      sub[, c("ref_q","alt_q","af_q") := NULL]
    }
  }

  msg("  After fill — NA ref: ", sum(is.na(sub$ref)))
  sub
}

st6_b <- fill_cohort(st6[cohort == "BREEDING"], "BREEDING")
st6_n <- fill_cohort(st6[cohort == "NATURAL"],  "NATURAL")
st6_o <- st6[!cohort %in% c("BREEDING","NATURAL")]

# Reassemble in original row order
st6_fixed <- rbind(st6_b, st6_n, st6_o, fill = TRUE)
setorder(st6_fixed, .row_idx)
st6_fixed[, .row_idx := NULL]

msg("Summary after fill:")
msg("  NA ref: ",  sum(is.na(st6_fixed$ref)))
msg("  NA alt: ",  sum(is.na(st6_fixed$alt)))
msg("  NA af:  ",  sum(is.na(st6_fixed$af)))

############################################################
# 5) WRITE CORRECTED XLSX
#    Load original workbook and overwrite only ref/alt/af columns
#    in Supplementary Table 6 so all other sheets/formatting are kept.
############################################################

msg("Loading original workbook...")
wb <- wb_load(XLSX_IN)

SHEET <- "Supplementary Table 6"
# Data rows in the xlsx start at row 4 (row 1=title, row 2=blank, row 3=headers)
XLSX_DATA_START <- 4L

msg("Writing ref column (col ", COL_REF, ") ...")
wb <- wb_add_data(wb, sheet = SHEET,
                  x = data.frame(ref = st6_fixed$ref),
                  start_row  = XLSX_DATA_START,
                  start_col  = COL_REF,
                  col_names  = FALSE)

msg("Writing alt column (col ", COL_ALT, ") ...")
wb <- wb_add_data(wb, sheet = SHEET,
                  x = data.frame(alt = st6_fixed$alt),
                  start_row  = XLSX_DATA_START,
                  start_col  = COL_ALT,
                  col_names  = FALSE)

msg("Writing af column (col ", COL_AF, ") ...")
wb <- wb_add_data(wb, sheet = SHEET,
                  x = data.frame(af = st6_fixed$af),
                  start_row  = XLSX_DATA_START,
                  start_col  = COL_AF,
                  col_names  = FALSE)

msg("Saving: ", XLSX_OUT)
wb_save(wb, XLSX_OUT)
msg("Saved: ", basename(XLSX_OUT))

############################################################
# 6) ALSO SAVE CORRECTED TSVs (for record)
############################################################

ANN17 <- file.path(PROJECT_ROOT, "RESULTS/JOINT/ANNOTATION17")
fwrite(st6_fixed[cohort == "BREEDING"],
       file.path(ANN17, "robust_breeding_snps_annotated_corrected.tsv"), sep="\t")
fwrite(st6_fixed[cohort == "NATURAL"],
       file.path(ANN17, "robust_natural_snps_annotated_corrected.tsv"),  sep="\t")
msg("Corrected TSVs saved to ", ANN17)

unlink(TMP_DIR, recursive = TRUE)
msg("Done.")
