#!/usr/bin/env Rscript
############################################################
# TreeGeneClimate (TGC) — JOINT ECS + TBS
# Step 17ab: Marker annotation against reference genome GFF3
#
# Annotates:
#   - ECS DAPC top-10 SNPs (per cohort x DF)
#   - TBS DAPC top-10 methylation sites (per cohort x context x DF)
#   - TBS KW SVMPs: breeding top-150 balanced + natural 6 formal SVMPs (step 8b)
#   - meQTL ROBUST markers: SNP + site positions significant at
#     p_FDR < 1e-10 in BOTH GENESIS5 AND MatrixEQTL5
#     (reads robust_markers_<cohort>.tsv from 15ab.R)
#
# Each marker is annotated with:
#   - gene_id, gene_start, gene_end, gene_strand (GFF3 gene feature)
#   - distance_bp (0 = overlaps gene; positive = nearest gene distance)
#   - annotation_class (genic / proximal_intergenic / distal_intergenic)
#   - ref, alt, af, dr2 (from imputed VCF; SNP markers only)
############################################################

suppressPackageStartupMessages({
  library(data.table)
})

options(stringsAsFactors = FALSE)

############################################################
# 1) SETTINGS
############################################################

COHORTS  <- c("BREEDING", "NATURAL")
CONTEXTS <- c("CpG", "CHG", "CHH")
TOOLS    <- c("GENESIS5", "MATRIXEQTL5")
TOP_N    <- 20L   # number of top markers to retain per DAPC discriminant function

# Methylation classification thresholds (context-specific).
# CHG/CHH methylation is generally much lower than CpG in conifers.
BETA_HYPO  <- c(CpG = 0.30, CHG = 0.20, CHH = 0.10)
BETA_HYPER <- c(CpG = 0.70, CHG = 0.50, CHH = 0.30)

############################################################
# 2) PATHS
############################################################

# === USER CONFIGURATION ===
PROJECT_ROOT <- "/path/to/your/project"  # <-- set this
# ===========================

# Picea abies v2.0 reference annotation (Nystedt et al. 2013 updated build)
GFF3_FILE <- file.path(PROJECT_ROOT,
  "REFERENCE/Pabies2.0/Picab02_230926_at01_all_sorted.gff3")

# Beagle-imputed VCFs (AF + DR2 available); used for allele-frequency annotation
VCF_FILES <- list(
  BREEDING = file.path(PROJECT_ROOT,
    "RESULTS/ECS/VCF_SPLIT/tgc.ecs.breeding.call.filt.maf05.snvs.poly.imputed.vcf.gz"),
  NATURAL  = file.path(PROJECT_ROOT,
    "RESULTS/ECS/VCF_SPLIT/tgc.ecs.natural.call.filt.maf05.snvs.poly.imputed.vcf.gz")
)

DAPC_ECS_ROOT <- file.path(PROJECT_ROOT, "RESULTS/ECS/RANALYSIS/TABLES/dapc_loadings")
DAPC_TBS_ROOT <- file.path(PROJECT_ROOT, "RESULTS/TBS/RANALYSIS/TABLES/dapc_loadings")
# Robust marker tables produced by 15ab.R
ROBUST_ROOT   <- file.path(PROJECT_ROOT, "RESULTS/JOINT/COMBINED5/overlap/tables")
SVMP_ROOT     <- file.path(PROJECT_ROOT, "RESULTS/TBS/RANALYSIS/TABLES/heatmap_markers_8B")

# Directory with mean methylation beta files (one per cohort × context).
# Expected filename: mean_beta_<cohort>_<context>.tsv  (columns: site, mean_beta)
# Computed automatically from M-value matrices in MQTL5_INPUTDIR if not present.
BETA_ROOT <- file.path(PROJECT_ROOT, "RESULTS/TBS/METHYLATION/mean_betas")

# M-value matrix inputs from 12ab0.R — used to compute mean betas on the fly.
# Path: MQTL5_INPUTDIR/<COHORT>/<CONTEXT>/methylation_mvalues_matrix.rds
# Rows = samples, columns = site IDs (chr:start-end).
MQTL5_INPUTDIR <- file.path(PROJECT_ROOT, "RESULTS/JOINT/MQTL5/INPUTS")

# eggNOG-mapper functional annotation (isoform-level; best mRNA selected per gene).
# Columns used: id (1), eggnog_description (5), eggnog_go (8), eggnog_KEGG_ko (9),
#               interpro_description/pfam (15), interpro_panther_description (20).
FUNC_ANNOT_FILE <- file.path(PROJECT_ROOT,
  "REFERENCE/Pabies2.0/Picab02_230926_at01_all_isoform_annotations_merged_sorted_non_redundant_panthers.tsv")

# List of mRNA IDs classified as TE-derived gene models (to build is_te_gene flag).
TE_IDS_FILE <- file.path(PROJECT_ROOT,
  "REFERENCE/Pabies2.0/TE_IDs_REMOVED.txt")

# Non-imputed (raw) VCFs — fallback for SNPs absent from the imputed VCF.
# GENESIS5 uses non-imputed GDS, so some SNPs may not appear in the imputed VCF.
RAW_VCF_FILES <- list(
  BREEDING = file.path(PROJECT_ROOT,
    "RESULTS/ECS/VCF_SPLIT/tgc.ecs.breeding.call.filt.maf05.snvs.poly.vcf.gz"),
  NATURAL  = file.path(PROJECT_ROOT,
    "RESULTS/ECS/VCF_SPLIT/tgc.ecs.natural.call.filt.maf05.snvs.poly.vcf.gz")
)

OUT_ROOT <- file.path(PROJECT_ROOT, "RESULTS/JOINT/ANNOTATION17")
TMP_DIR  <- file.path(OUT_ROOT, "tmp")
dir.create(TMP_DIR, recursive = TRUE, showWarnings = FALSE)

############################################################
# 3) HELPERS
############################################################

msg <- function(...) cat(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n"))

# Write a sorted 0-based BED file for a set of markers.
# chr_col, pos_col (1-based genomic), id_col are column names in dt.
write_marker_bed <- function(dt, chr_col, pos_col, id_col, path) {
  bed <- data.table(
    chr   = dt[[chr_col]],
    start = as.numeric(dt[[pos_col]]) - 1,  # BED is 0-based half-open
    end   = as.numeric(dt[[pos_col]]),
    id    = dt[[id_col]]
  )
  bed <- bed[!is.na(chr) & !is.na(start)]
  setorder(bed, chr, start)
  fwrite(bed, path, sep = "\t", col.names = FALSE)
  invisible(path)
}

# Run bedtools closest -d against genes BED and return annotated data.table.
# Marker BED has 4 columns: chr, start, end, id.
# Gene BED has 6 columns: chr, start, end, gene_id, ., strand.
run_bedtools_closest <- function(marker_bed, genes_bed) {
  if (!file.exists(marker_bed) || file.info(marker_bed)$size == 0)
    return(data.table())

  out_file <- tempfile(tmpdir = TMP_DIR, fileext = ".closest.tsv")
  cmd <- sprintf("bedtools closest -a %s -b %s -d > %s",
                 marker_bed, genes_bed, out_file)
  ret <- system(cmd)
  if (ret != 0 || !file.exists(out_file) || file.info(out_file)$size == 0) {
    msg("  WARNING: bedtools closest returned no results for ", marker_bed)
    return(data.table())
  }

  # Cols 1-4: marker; cols 5-10: gene (chr,start,end,id,.,strand); col 11: distance
  dt <- fread(out_file, header = FALSE, sep = "\t", fill = TRUE,
              col.names = c("chr","pos_start","pos_end","marker_id",
                            "gene_chr","gene_start","gene_end","gene_id",
                            "gene_score","gene_strand","distance_bp"))
  file.remove(out_file)

  dt[, gene_score := NULL]
  # distance = -1 means no gene on that chromosome at all
  dt[distance_bp < 0, `:=`(gene_id     = "no_gene_on_chrom",
                             gene_chr    = NA_character_,
                             gene_start  = NA_integer_,
                             gene_end    = NA_integer_,
                             gene_strand = NA_character_,
                             distance_bp = NA_integer_)]
  # bedtools returns "." when there is no overlapping feature
  dt[gene_chr == ".", gene_chr := NA_character_]
  dt
}

# Query a VCF for REF, ALT, AF and optionally DR2 at given marker positions.
# has_dr2 = FALSE for raw (non-imputed) VCFs that lack the Beagle DR2 INFO field.
query_vcf <- function(markers_dt, vcf_file, tmp_prefix, has_dr2 = TRUE) {
  empty <- data.table(chr = character(), pos = integer(),
                      ref = character(), alt = character(),
                      af  = numeric(),  dr2 = numeric())
  if (!file.exists(vcf_file)) {
    msg("  WARNING: VCF not found: ", vcf_file)
    return(empty)
  }
  markers_dt <- unique(markers_dt[!is.na(chr) & !is.na(pos)])
  if (nrow(markers_dt) == 0) return(empty)

  # Write 0-based BED regions for bcftools -R
  reg_file <- paste0(tmp_prefix, ".regions.bed")
  reg <- markers_dt[, .(chr, start = pos - 1L, end = pos)]
  setorder(reg, chr, start)
  fwrite(reg, reg_file, sep = "\t", col.names = FALSE)

  out_file <- paste0(tmp_prefix, ".vcf_query.tsv")
  # vcf_fmt is substituted as a value (not a sprintf format), so use single %
  vcf_fmt  <- if (has_dr2)
    "'%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/AF\\t%INFO/DR2\\n'"
  else
    "'%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/AF\\n'"
  vcf_cols <- if (has_dr2)
    c("chr","pos","ref","alt","af","dr2")
  else
    c("chr","pos","ref","alt","af")
  cmd <- sprintf("bcftools query -R %s -f %s %s > %s 2>/dev/null",
    reg_file, vcf_fmt, vcf_file, out_file)
  system(cmd)
  file.remove(reg_file)

  if (!file.exists(out_file) || file.info(out_file)$size == 0) {
    return(empty)
  }
  # na.strings="." converts bcftools missing-value dots to NA.
  # Explicit numeric cast prevents fread from typing all-NA columns as logical,
  # which would coerce numeric AF values to TRUE/FALSE when assigned later.
  res <- fread(out_file, header = FALSE, sep = "\t", fill = TRUE,
               col.names = vcf_cols, na.strings = ".")
  res[, pos := as.integer(pos)]
  if ("af"  %in% names(res)) res[, af  := as.numeric(af)]
  if ("dr2" %in% names(res)) res[, dr2 := as.numeric(dr2)]
  if (!"dr2" %in% names(res)) res[, dr2 := NA_real_]
  res
}

# For SNPs with NA ref/alt or NA af, query the raw (non-imputed) VCF.
# The imputed VCF strips INFO for directly genotyped markers, so AF is only
# available in the raw VCF for those sites. DR2 is never in the raw VCF.
# Returns dt with ref/alt/af filled where previously NA.
fill_from_raw_vcf <- function(dt, cohort, tmp_prefix) {
  raw_vcf <- RAW_VCF_FILES[[cohort]]
  if (is.null(raw_vcf) || !file.exists(raw_vcf)) return(dt)
  if (!all(c("ref","alt","af","chr","pos") %in% names(dt))) return(dt)
  missing <- unique(dt[is.na(ref) | is.na(af), .(chr, pos = as.integer(pos))])
  if (nrow(missing) == 0) return(dt)
  msg("  Raw-VCF fill for ", nrow(missing),
      " SNPs with NA REF/ALT or NA AF (imputed VCF strips INFO for typed markers)...")
  vi <- query_vcf(missing, raw_vcf, tmp_prefix, has_dr2 = FALSE)
  if (nrow(vi) == 0) return(dt)
  dt <- merge(dt,
    vi[, .(chr, pos, ref_raw = ref, alt_raw = alt, af_raw = af)],
    by = c("chr", "pos"), all.x = TRUE)
  dt[is.na(ref), `:=`(ref = ref_raw, alt = alt_raw)]
  dt[is.na(af) & !is.na(af_raw), af := af_raw]
  dt[, c("ref_raw","alt_raw","af_raw") := NULL]
  msg("    After fill — REF known: ", dt[!is.na(ref), .N],
      " | AF known: ", dt[!is.na(af), .N])
  dt
}

# Load eggNOG-mapper functional annotation.
# Returns a data.table keyed by gene_id with columns:
#   eggnog_description, pfam_domain, panther_description, go_terms, kegg_ko.
# One row per gene (best/first mRNA isoform).
load_functional_annotation <- function() {
  empty <- data.table(gene_id = character(), eggnog_description = character(),
                      pfam_domain = character(), panther_description = character(),
                      go_terms = character(), kegg_ko = character())
  if (!file.exists(FUNC_ANNOT_FILE)) {
    msg("  WARNING: functional annotation file not found — eggNOG columns will be NA")
    return(empty)
  }
  msg("Loading eggNOG-mapper functional annotation...")
  # Select only needed columns to keep memory use low (~800 MB TSV)
  fa <- fread(FUNC_ANNOT_FILE, header = TRUE,
              select = c(1L, 5L, 8L, 9L, 15L, 20L),
              col.names = c("mrna_id", "eggnog_description", "go_terms",
                            "kegg_ko", "pfam_domain", "panther_description"))
  # Derive gene_id by stripping .mRNA.<n> suffix
  fa[, gene_id := sub("\\.mRNA\\..*$", "", mrna_id)]
  # Keep first row per gene (sorted by mrna_id, so .mRNA.1 comes first)
  fa <- fa[order(mrna_id)][!duplicated(gene_id)]
  fa[, mrna_id := NULL]
  # Replace literal "NA" strings with proper NA
  for (col in c("eggnog_description","go_terms","kegg_ko","pfam_domain","panther_description"))
    set(fa, which(fa[[col]] == "NA"), col, NA_character_)
  setkey(fa, gene_id)
  msg("  Loaded: ", nrow(fa), " gene functional annotations")
  fa
}

# Load TE-derived gene IDs from TE_IDS_REMOVED.txt.
# Returns a character vector of gene IDs (without mRNA suffix).
load_te_gene_ids <- function() {
  if (!file.exists(TE_IDS_FILE)) {
    msg("  WARNING: TE IDs file not found — is_te_gene will be FALSE for all markers")
    return(character(0))
  }
  mrna_ids <- fread(TE_IDS_FILE, header = FALSE, col.names = "mrna_id")$mrna_id
  gene_ids <- unique(sub("\\.mRNA\\..*$", "", mrna_ids))
  msg("  TE-derived gene IDs loaded: ", length(gene_ids))
  gene_ids
}

# Join functional annotation and TE flag onto any marker/gene table.
# Requires a gene_id column. Adds eggnog_description, pfam_domain,
# panther_description, go_terms, kegg_ko, is_te_gene.
add_functional_annotation <- function(dt, func_annot, te_genes) {
  if (!"gene_id" %in% names(dt)) return(dt)
  if (nrow(func_annot) > 0) {
    dt <- merge(dt, func_annot, by = "gene_id", all.x = TRUE)
  } else {
    dt[, eggnog_description  := NA_character_]
    dt[, pfam_domain         := NA_character_]
    dt[, panther_description := NA_character_]
    dt[, go_terms            := NA_character_]
    dt[, kegg_ko             := NA_character_]
  }
  # Flag gene models with known TE origin — markers in TE genes are noted separately
  dt[, is_te_gene := gene_id %in% te_genes]
  dt
}

# Classify genomic context based on distance to nearest gene.
# Proximal intergenic is defined as ≤2 kb — a conservative promoter-proximal window.
add_annotation_class <- function(dt) {
  dt[, annotation_class := fcase(
    distance_bp == 0L,                           "genic",
    distance_bp > 0L & distance_bp <= 2000L,     "proximal_intergenic",
    distance_bp > 2000L,                          "distal_intergenic",
    default =                                     "no_gene_on_chrom"
  )]
  dt
}

# Load mean beta per site for one cohort × context.
# Returns data.table(site, context, cohort, mean_beta) or NULL if file absent.
load_site_betas <- function(cohort, ctx) {
  if (is.null(BETA_ROOT)) return(NULL)
  f <- file.path(BETA_ROOT,
    sprintf("mean_beta_%s_%s.tsv", tolower(cohort), ctx))
  if (!file.exists(f)) {
    msg("  NOTE: beta file not found: ", basename(f), " — meth_status_absolute = NA")
    return(NULL)
  }
  dt <- fread(f, header = TRUE)
  if (!all(c("site", "mean_beta") %in% names(dt))) {
    msg("  WARNING: ", basename(f), " needs columns: site, mean_beta — skipping")
    return(NULL)
  }
  dt[, site    := as.character(site)]
  dt[, context := ctx]
  dt[, cohort  := cohort]
  dt[, .(site, context, cohort, mean_beta)]
}

# Classify absolute methylation status using context-specific thresholds.
classify_meth_status <- function(beta, ctx) {
  hypo  <- BETA_HYPO [ctx]
  hyper <- BETA_HYPER[ctx]
  ifelse(is.na(beta) | is.na(hypo), NA_character_,
    ifelse(beta < hypo,  "hypomethylated",
      ifelse(beta >= hyper, "hypermethylated", "intermediate")))
}

# Join mean beta values + Option 1 absolute classification onto a sites table.
# site_col: name of the column holding site IDs (default "site").
# The table must already have columns: context, cohort.
add_meth_status <- function(sites_dt, site_col = "site") {
  beta_rows <- list()
  for (co in COHORTS)
    for (ctx in CONTEXTS) {
      b <- load_site_betas(co, ctx)
      if (!is.null(b)) beta_rows[[length(beta_rows) + 1L]] <- b
    }
  beta_all <- if (length(beta_rows)) rbindlist(beta_rows) else data.table()

  if (nrow(beta_all) > 0) {
    dt <- copy(sites_dt)
    if (site_col != "site") setnames(dt, site_col, "site")
    dt <- merge(dt, beta_all, by = c("site", "context", "cohort"), all.x = TRUE)
    if (site_col != "site") setnames(dt, "site", site_col)
    # mapply calls classify_meth_status row-wise over the vectorised table
    dt[, meth_status_absolute := mapply(classify_meth_status, mean_beta, context)]
    return(dt)
  }
  # No beta files available — fill columns with NA so schema remains consistent
  sites_dt[, mean_beta            := NA_real_]
  sites_dt[, meth_status_absolute := NA_character_]
  sites_dt
}

# Back-transform M-values to beta values and return per-site mean betas.
# Uses EPSILON_M = 0.5 matching 12ab0.R methylkit_to_mvalues().
# Inverse: beta = (100.5 * exp(M) - 0.5) / (100 * (1 + exp(M)))
# M-value matrix: rows = samples, columns = site IDs (chr:start-end).
EPSILON_M <- 0.5

compute_mean_beta <- function(cohort, ctx) {
  mat_file <- file.path(MQTL5_INPUTDIR,
    toupper(cohort), ctx, "methylation_mvalues_matrix.rds")
  if (!file.exists(mat_file)) {
    msg("  NOTE: M-value matrix not found: ", mat_file)
    return(NULL)
  }
  msg("  Computing mean betas from M-value matrix: ", cohort, " / ", ctx)
  mmat <- readRDS(mat_file)
  if (is.data.frame(mmat)) mmat <- as.matrix(mmat)
  # colMeans across samples gives the average M-value per site
  m_means    <- colMeans(mmat, na.rm = TRUE)
  beta_means <- (100.5 * exp(m_means) - EPSILON_M) /
                (100   * (1 + exp(m_means)))
  data.table(site = names(m_means), mean_beta = as.numeric(beta_means))
}

# Build gene-centric summary from the combined annotation table.
# func_annot and te_genes are joined to enrich the summary with functional info.
make_gene_summary <- function(dt, func_annot, te_genes) {
  # Rank annotation classes: genic > proximal > distal; report best class per gene
  ann_order <- c(genic = 1L, proximal_intergenic = 2L, distal_intergenic = 3L)
  dt_g <- dt[!is.na(gene_id) & gene_id != "no_gene_on_chrom"]
  if (nrow(dt_g) == 0L) return(data.table())

  # Collapse all markers per gene to one summary row
  summ <- dt_g[, {
    ann_u  <- unique(annotation_class[!is.na(annotation_class)])
    ranked <- ann_u[ann_u %in% names(ann_order)]
    best   <- if (length(ranked)) names(sort(ann_order[ranked]))[1L] else ann_u[1L]
    .(
      gene_chr              = gene_chr[1L],
      gene_start            = gene_start[1L],
      gene_end              = gene_end[1L],
      gene_strand           = gene_strand[1L],
      best_annotation_class = best,
      min_distance_bp       = min(distance_bp, na.rm = TRUE),
      n_markers             = uniqueN(marker_id),
      n_snps                = uniqueN(marker_id[marker_type == "SNP"]),
      n_sites               = uniqueN(marker_id[marker_type == "methylation_site"]),
      marker_ids            = paste(sort(unique(marker_id)),          collapse = ";"),
      sources               = paste(sort(unique(na.omit(source))),    collapse = ";"),
      tools                 = paste(sort(unique(na.omit(tool))),      collapse = ";"),
      cohorts               = paste(sort(unique(na.omit(cohort))),    collapse = ";"),
      contexts              = paste(sort(unique(na.omit(context))),   collapse = ";")
    )
  }, by = gene_id]

  # Join functional annotation and TE flag
  summ <- add_functional_annotation(summ, func_annot, te_genes)

  setorder(summ, gene_chr, gene_start)
  summ
}

############################################################
# 3b) PRE-COMPUTE MEAN BETA FILES (from M-value matrices if not present)
############################################################

msg("Checking / computing mean beta files from M-value matrices...")
dir.create(BETA_ROOT, recursive = TRUE, showWarnings = FALSE)

# Pre-compute beta files up front so all annotation sections can use load_site_betas()
for (.co in COHORTS) {
  for (.ctx in CONTEXTS) {
    out_beta <- file.path(BETA_ROOT,
      sprintf("mean_beta_%s_%s.tsv", tolower(.co), .ctx))
    if (file.exists(out_beta)) {
      msg("  Found: ", basename(out_beta), " — skipping recompute")
      next
    }
    mb <- compute_mean_beta(.co, .ctx)
    if (!is.null(mb)) {
      fwrite(mb, out_beta, sep = "\t")
      msg("  Written: ", basename(out_beta), " (", nrow(mb), " sites)")
    }
  }
}
rm(.co, .ctx)

############################################################
# 4) PREPARE GFF3 GENE BED
############################################################

msg("Preparing GFF3 gene BED (filtering gene features)...")

genes_bed_path <- file.path(TMP_DIR, "genes.bed")

# Extract gene-level features only; create 0-based BED
# Attributes field example: ID=PA_chr01_G000001;
gff_cmd <- sprintf(
  "awk 'BEGIN{OFS=\"\\t\"} !/^#/ && $3==\"gene\" { \
     match($9, /ID=([^;]+)/, a); \
     print $1, $4-1, $5, a[1], \".\", $7 \
   }' %s | sort -k1,1 -k2,2n > %s",
  GFF3_FILE, genes_bed_path)

system(gff_cmd)
n_genes <- as.integer(system(sprintf("wc -l < %s", genes_bed_path), intern = TRUE))
msg("  ", n_genes, " gene features written to genes.bed")

# Load functional annotation and TE gene IDs (used in all subsequent sections)
func_annot <- load_functional_annotation()
te_genes   <- load_te_gene_ids()

############################################################
# 5) ECS DAPC — TOP-10 SNPs PER COHORT x DF
############################################################

msg("======================================================")
msg("ECS DAPC — top-10 SNPs per cohort x DF")

ecs_files <- list(
  BREEDING = file.path(DAPC_ECS_ROOT, "breeding_dapc_loadings_all_DF1_DF2.csv"),
  NATURAL  = file.path(DAPC_ECS_ROOT, "natural_dapc_loadings_all_DF1_DF2.csv")
)

ecs_list <- list()
for (cohort in COHORTS) {
  f <- ecs_files[[cohort]]
  if (!file.exists(f)) { msg("  Missing: ", f); next }
  dt <- fread(f, header = TRUE)
  # Standardise column names (CSV header: chromosome,position,DF,loading)
  setnames(dt, names(dt), c("chr", "pos", "DF", "loading"))
  dt[, abs_loading := abs(loading)]
  dt[, cohort := cohort]
  # Top N by |loading| within each DF; loading magnitude indicates discriminatory power
  top <- dt[order(-abs_loading), .SD[seq_len(min(TOP_N, .N))], by = DF]
  top[, abs_loading := NULL]
  ecs_list[[cohort]] <- top
}

ecs_top <- rbindlist(ecs_list, fill = TRUE)

if (nrow(ecs_top) > 0) {
  ecs_top[, marker_id := paste0(chr, ":", pos)]

  # GFF3 annotation via bedtools
  bed_path <- file.path(TMP_DIR, "ecs_dapc.bed")
  write_marker_bed(ecs_top, "chr", "pos", "marker_id", bed_path)
  ann <- run_bedtools_closest(bed_path, genes_bed_path)

  if (nrow(ann) > 0) {
    ecs_top <- merge(ecs_top,
                     ann[, .(marker_id, gene_id, gene_chr, gene_start,
                              gene_end, gene_strand, distance_bp)],
                     by = "marker_id", all.x = TRUE)
  }

  # VCF allele info (per cohort) — imputed VCF first, raw VCF fallback for NA
  vcf_parts <- lapply(COHORTS, function(co) {
    sub <- ecs_top[cohort == co]
    vi  <- query_vcf(data.table(chr = sub$chr, pos = as.integer(sub$pos)),
                     VCF_FILES[[co]],
                     file.path(TMP_DIR, paste0("ecs_vcf_", co)))
    if (nrow(vi) > 0)
      sub <- merge(sub, vi, by.x = c("chr","pos"), by.y = c("chr","pos"), all.x = TRUE)
    else
      sub[, c("ref","alt","af","dr2") := list(NA_character_, NA_character_,
                                               NA_real_, NA_real_)]
    sub <- fill_from_raw_vcf(sub, co,
             file.path(TMP_DIR, paste0("ecs_raw_vcf_", co)))
    sub
  })
  ecs_top <- rbindlist(vcf_parts, fill = TRUE)
  ecs_top <- add_annotation_class(ecs_top)

  # SNP-specific metadata; methylation fields set to NA for non-site markers
  ecs_top[, source                   := "ECS_DAPC"]
  ecs_top[, marker_type              := "SNP"]
  ecs_top[, context                  := NA_character_]
  ecs_top[, methylation_direction    := NA_character_]
  ecs_top[, mean_beta                := NA_real_]
  ecs_top[, meth_status_absolute     := NA_character_]
  ecs_top[, higher_methylation_cohort := NA_character_]

  ecs_top <- add_functional_annotation(ecs_top, func_annot, te_genes)
  fwrite(ecs_top, file.path(OUT_ROOT, "ecs_dapc_top20_annotated.tsv"), sep = "\t")
  msg("  Saved: ecs_dapc_top20_annotated.tsv (", nrow(ecs_top), " markers)")
} else {
  msg("  No ECS DAPC markers found.")
  ecs_top <- data.table()
}

############################################################
# 6) TBS DAPC — TOP-10 SITES PER COHORT x CONTEXT x DF
############################################################

msg("======================================================")
msg("TBS DAPC — top-10 methylation sites per cohort x context x DF")

tbs_list <- list()
for (cohort in COHORTS) {
  for (ctx in CONTEXTS) {
    f <- file.path(DAPC_TBS_ROOT,
                   sprintf("TBS_DAPC_loadings_%s_%s_ALL.tsv",
                           tolower(cohort), tolower(ctx)))
    if (!file.exists(f)) { msg("  Missing: ", f); next }
    dt <- fread(f, header = TRUE)
    # TSV header: loc  chr  pos  DF  loading
    setnames(dt, names(dt), c("loc", "chr", "pos", "DF", "loading"))
    dt[, abs_loading := abs(loading)]
    dt[, cohort  := cohort]
    dt[, context := ctx]
    top <- dt[order(-abs_loading), .SD[seq_len(min(TOP_N, .N))], by = DF]
    top[, abs_loading := NULL]
    tbs_list[[paste(cohort, ctx)]] <- top
  }
}

tbs_top <- rbindlist(tbs_list, fill = TRUE)

if (nrow(tbs_top) > 0) {
  # GFF3 annotation via bedtools
  bed_path <- file.path(TMP_DIR, "tbs_dapc.bed")
  write_marker_bed(tbs_top, "chr", "pos", "loc", bed_path)
  ann <- run_bedtools_closest(bed_path, genes_bed_path)

  if (nrow(ann) > 0) {
    setnames(ann, "marker_id", "loc")
    tbs_top <- merge(tbs_top,
                     ann[, .(loc, gene_id, gene_chr, gene_start,
                              gene_end, gene_strand, distance_bp)],
                     by = "loc", all.x = TRUE)
  }

  # Methylation sites are not in the VCF — set SNP-specific columns to NA
  tbs_top[, c("ref","alt","af","dr2") := list(NA_character_, NA_character_,
                                               NA_real_, NA_real_)]
  tbs_top <- add_annotation_class(tbs_top)
  tbs_top[, source      := "TBS_DAPC"]
  tbs_top[, marker_type := "methylation_site"]

  # Option 1 — absolute methylation status (via loading sign + optional beta files)
  # methylation_direction: sign of DAPC loading (positive/negative relative to DF axis)
  tbs_top[, methylation_direction := ifelse(loading > 0, "positive", "negative")]
  tbs_top <- add_meth_status(tbs_top, site_col = "loc")
  # Cross-cohort comparison not applicable for DAPC (per-cohort analysis)
  tbs_top[, higher_methylation_cohort := NA_character_]

  setnames(tbs_top, "loc", "marker_id")

  tbs_top <- add_functional_annotation(tbs_top, func_annot, te_genes)
  fwrite(tbs_top, file.path(OUT_ROOT, "tbs_dapc_top20_annotated.tsv"), sep = "\t")
  msg("  Saved: tbs_dapc_top20_annotated.tsv (", nrow(tbs_top), " markers)")
} else {
  msg("  No TBS DAPC markers found.")
  tbs_top <- data.table()
}

############################################################
# 6b) TBS KW SVMPs — SELECTED MARKERS FROM STEP 8b
#
# Breeding: top-150 balanced (50 per context, all formal SVMPs)
# Natural:  6 formal SVMPs only (3 CpG + 3 CHH; padj < 0.05)
############################################################

msg("======================================================")
msg("TBS KW SVMPs — annotating selected markers from step 8b")

svmp_files <- list(
  BREEDING = file.path(SVMP_ROOT, "TBS_8B_selected_markers_breeding_top150_BALANCED.tsv"),
  NATURAL  = file.path(SVMP_ROOT, "TBS_8B_selected_markers_natural_formal_SVMPs.tsv")
)

svmp_all_list <- list()

for (cohort in COHORTS) {
  f <- svmp_files[[cohort]]
  if (!file.exists(f)) {
    msg("  Missing: ", basename(f), " — skipping ", cohort)
    next
  }
  dt <- fread(f, header = TRUE)
  msg("  ", cohort, ": ", nrow(dt), " SVMPs loaded")

  # Standardise column names
  if ("Context" %in% names(dt) && !"context" %in% names(dt))
    setnames(dt, "Context", "context")
  dt[, cohort := cohort]
  dt[, pos    := as.integer(start)]   # 1-based genomic position alias

  # BEDTools annotation
  bed_path <- file.path(TMP_DIR,
    sprintf("tbs_svmp_%s.bed", tolower(cohort)))
  write_marker_bed(dt, "chr", "pos", "marker_id", bed_path)
  ann <- run_bedtools_closest(bed_path, genes_bed_path)

  if (nrow(ann) > 0)
    dt <- merge(dt,
                ann[, .(marker_id, gene_id, gene_chr, gene_start,
                         gene_end, gene_strand, distance_bp)],
                by = "marker_id", all.x = TRUE)

  dt <- add_annotation_class(dt)

  # Mean beta + hypo/hyper classification from pre-computed beta files.
  # add_meth_status joins on (site = loc, context, cohort).
  dt <- add_meth_status(dt, site_col = "loc")

  # Standard marker metadata columns; SNP-specific fields are NA for methylation sites
  dt[, source                    := "TBS_KW_SVMP"]
  dt[, marker_type               := "methylation_site"]
  dt[, tool                      := NA_character_]
  dt[, DF                        := NA_character_]
  dt[, loading                   := NA_real_]
  dt[, methylation_direction     := NA_character_]
  dt[, higher_methylation_cohort := NA_character_]
  dt[, c("ref","alt","af","dr2") := list(NA_character_, NA_character_,
                                          NA_real_, NA_real_)]
  # meQTL association columns left NA — SVMPs are not from the meQTL analysis
  dt[, n_associated_sites        := NA_integer_]
  dt[, associated_sites          := NA_character_]
  dt[, associated_contexts       := NA_character_]
  dt[, associated_site_chrs      := NA_character_]
  dt[, associated_site_pos       := NA_character_]

  dt <- add_functional_annotation(dt, func_annot, te_genes)

  out_f <- file.path(OUT_ROOT,
    sprintf("tbs_svmp_%s_annotated.tsv", tolower(cohort)))
  fwrite(dt, out_f, sep = "\t")
  msg("  Saved: ", basename(out_f), " (", nrow(dt), " markers)")
  cat("  Annotation class breakdown:\n")
  print(dt[, .N, by = .(context, annotation_class)][order(context)])
  cat("  Methylation status:\n")
  print(dt[, .N, by = .(context, meth_status_absolute)][order(context)])

  svmp_all_list[[cohort]] <- dt
}

if (length(svmp_all_list) == 0)
  msg("  No SVMP files found — section 6b produced no output.")

############################################################
# 7) meQTL ROBUST MARKERS (p_FDR < 1e-10 in BOTH GENESIS5 AND MatrixEQTL5)
############################################################

msg("======================================================")
msg("meQTL robust markers (p_FDR < 1e-10 in GENESIS5 AND MatrixEQTL5)")

meqtl_results <- list()

for (cohort in COHORTS) {
  cohort_key <- tolower(cohort)

  # Robust marker tables were produced by 15ab.R (pair-level intersect of both tools)
  rob_file <- file.path(ROBUST_ROOT,
    sprintf("robust_markers_%s.tsv", cohort_key))

  if (!file.exists(rob_file)) {
    msg("  WARNING: ", basename(rob_file),
        " not found — run 15ab.R first")
    next
  }

  dt <- fread(rob_file, header = TRUE)
  if ("snp"  %in% names(dt)) dt[, snp  := as.character(snp)]
  if ("site" %in% names(dt)) dt[, site := as.character(site)]
  msg("  ", cohort, ": ", nrow(dt), " robust pairs | ",
      uniqueN(dt$context), " contexts | ",
      uniqueN(dt$snp), " SNPs | ", uniqueN(dt$site), " sites")

  # ---- Associated epimarkers per SNP (aggregated before deduplication) ----
  # Each robust SNP may associate with multiple methylation sites; capture all
  # of them before deduplicating SNP positions for the per-SNP annotation table
  snp_site_map <- dt[!is.na(snp) & !is.na(site), .(
    n_associated_sites   = uniqueN(site),
    associated_sites     = paste(sort(unique(site)),             collapse = ";"),
    associated_contexts  = paste(sort(unique(context)),          collapse = ";"),
    associated_site_chrs = paste(sort(unique(as.character(site_chr))), collapse = ";"),
    associated_site_pos  = paste(sort(unique(as.character(site_pos))), collapse = ";")
  ), by = snp]

  # ---- SNP positions (unique by snp_chr + snp_pos) ----
  snps <- unique(dt[!is.na(snp_chr) & !is.na(snp_pos),
                    .(chr = snp_chr,
                      pos = as.numeric(snp_pos),
                      snp,
                      context)])
  snps[, marker_id := paste0("snp_", snp, "_", chr, ":", as.integer(pos))]

  if (nrow(snps) > 0) {
    bed_path <- file.path(TMP_DIR,
      sprintf("robust_%s_snps.bed", cohort_key))
    write_marker_bed(snps, "chr", "pos", "marker_id", bed_path)
    ann <- run_bedtools_closest(bed_path, genes_bed_path)

    if (nrow(ann) > 0)
      snps <- merge(snps,
                    ann[, .(marker_id, gene_id, gene_chr, gene_start,
                             gene_end, gene_strand, distance_bp)],
                    by = "marker_id", all.x = TRUE)

    # Primary VCF query (imputed — has AF + DR2)
    vi <- query_vcf(
      data.table(chr = snps$chr, pos = as.integer(snps$pos)),
      VCF_FILES[[cohort]],
      file.path(TMP_DIR, sprintf("robust_%s_snps_vcf", cohort_key)))
    if (nrow(vi) > 0)
      snps <- merge(snps, vi, by.x = c("chr","pos"), by.y = c("chr","pos"),
                    all.x = TRUE)
    else
      snps[, c("ref","alt","af","dr2") := list(NA_character_, NA_character_,
                                                NA_real_, NA_real_)]

    # Fill NA ref/alt/af from raw VCF (covers GENESIS5 SNPs absent from imputed
    # VCF, and genotyped markers whose INFO was stripped by Beagle imputation).
    snps <- fill_from_raw_vcf(snps, cohort,
      file.path(TMP_DIR, sprintf("robust_%s_snps_rawvcf", cohort_key)))

    # Join associated epimarker information (sites, contexts) aggregated above
    snps <- merge(snps, snp_site_map, by = "snp", all.x = TRUE)

    snps <- add_annotation_class(snps)
    snps[, cohort                    := cohort]
    snps[, tool                      := "GENESIS5+MATRIXEQTL5"]
    snps[, source                    := "robust_meQTL"]
    snps[, marker_type               := "SNP"]
    snps[, DF                        := NA_character_]
    snps[, loading                   := NA_real_]
    snps[, methylation_direction     := NA_character_]
    snps[, mean_beta                 := NA_real_]
    snps[, meth_status_absolute      := NA_character_]
    snps[, higher_methylation_cohort := NA_character_]

    snps <- add_functional_annotation(snps, func_annot, te_genes)
    out_f <- file.path(OUT_ROOT,
      sprintf("robust_%s_snps_annotated.tsv", cohort_key))
    fwrite(snps, out_f, sep = "\t")
    msg("  Saved: ", basename(out_f), " (", nrow(snps), " unique SNP positions)")
    meqtl_results[[paste(cohort, "snp")]] <- snps
  }

  # ---- Methylation site positions (unique by site_chr + site_pos) ----
  sites <- unique(dt[!is.na(site_chr) & !is.na(site_pos),
                     .(chr = site_chr,
                       pos = as.integer(site_pos),
                       site,
                       context)])
  sites[, marker_id := site]

  if (nrow(sites) > 0) {
    bed_path <- file.path(TMP_DIR,
      sprintf("robust_%s_sites.bed", cohort_key))
    write_marker_bed(sites, "chr", "pos", "marker_id", bed_path)
    ann <- run_bedtools_closest(bed_path, genes_bed_path)

    if (nrow(ann) > 0)
      sites <- merge(sites,
                     ann[, .(marker_id, gene_id, gene_chr, gene_start,
                              gene_end, gene_strand, distance_bp)],
                     by = "marker_id", all.x = TRUE)

    # Methylation sites have no VCF allele information
    sites[, c("ref","alt","af","dr2") := list(NA_character_, NA_character_,
                                               NA_real_, NA_real_)]
    sites <- add_annotation_class(sites)
    sites[, cohort                   := cohort]
    sites[, tool                     := "GENESIS5+MATRIXEQTL5"]
    sites[, source                   := "robust_meQTL"]
    sites[, marker_type              := "methylation_site"]
    sites[, DF                       := NA_character_]
    sites[, loading                  := NA_real_]
    sites[, methylation_direction    := NA_character_]

    # Option 1 — absolute methylation status per site per cohort
    sites <- add_meth_status(sites)

    # higher_methylation_cohort filled in cross-cohort step below
    sites[, higher_methylation_cohort := NA_character_]

    sites <- add_functional_annotation(sites, func_annot, te_genes)
    out_f <- file.path(OUT_ROOT,
      sprintf("robust_%s_sites_annotated.tsv", cohort_key))
    fwrite(sites, out_f, sep = "\t")
    msg("  Saved: ", basename(out_f), " (", nrow(sites), " unique site positions)")
    meqtl_results[[paste(cohort, "site")]] <- sites
  }

  rm(dt); gc()
}

# Option 2 — cross-cohort relative methylation comparison.
# For each site present in both cohorts, flag which has higher mean beta.
msg("--- Cross-cohort methylation direction (Option 2) ---")
{
  breed_s <- meqtl_results[["BREEDING site"]]
  nat_s   <- meqtl_results[["NATURAL site"]]
  has_beta_b <- !is.null(breed_s) && "mean_beta" %in% names(breed_s)
  has_beta_n <- !is.null(nat_s)   && "mean_beta" %in% names(nat_s)

  if (has_beta_b && has_beta_n) {
    breed_b <- breed_s[!is.na(mean_beta), .(site, context, beta_b = mean_beta)]
    nat_b   <- nat_s  [!is.na(mean_beta), .(site, context, beta_n = mean_beta)]
    # Full outer join: include sites unique to one cohort as well as shared sites
    cross   <- merge(breed_b, nat_b, by = c("site", "context"), all = TRUE)
    cross[, higher_methylation_cohort := fcase(
      !is.na(beta_b) & !is.na(beta_n) & beta_b > beta_n, "BREEDING",
      !is.na(beta_b) & !is.na(beta_n) & beta_b < beta_n, "NATURAL",
      !is.na(beta_b) & !is.na(beta_n),                   "equal",
      !is.na(beta_b) & is.na(beta_n),                    "BREEDING_only",
      is.na(beta_b)  & !is.na(beta_n),                   "NATURAL_only",
      default = NA_character_
    )]
    cross[, c("beta_b", "beta_n") := NULL]

    # Update per-cohort site tables with the cross-cohort direction column
    for (co in COHORTS) {
      key <- paste(co, "site")
      if (!is.null(meqtl_results[[key]])) {
        meqtl_results[[key]] <- merge(
          meqtl_results[[key]][, higher_methylation_cohort := NULL],
          cross, by = c("site", "context"), all.x = TRUE)
        out_f <- file.path(OUT_ROOT,
          sprintf("robust_%s_sites_annotated.tsv", tolower(co)))
        fwrite(meqtl_results[[key]], out_f, sep = "\t")
        msg("  Updated: ", basename(out_f), " (higher_methylation_cohort added)")
      }
    }
    msg("  Sites in both cohorts: ",
        cross[!is.na(higher_methylation_cohort) &
              higher_methylation_cohort %in% c("BREEDING","NATURAL","equal"), .N])
  } else {
    msg("  Skipping: beta files absent for one or both cohorts")
  }
}

############################################################
# 7b) CROSS-REFERENCE: ECS DAPC top SNPs × robust meQTL results
#
# For each ECS DAPC top SNP, checks whether it also appears as a robust meQTL
# SNP (same cohort, matched by chr+pos). If so, fills in the associated
# epimarker columns (n_associated_sites, associated_sites, associated_contexts,
# associated_site_chrs, associated_site_pos) from the meQTL snp_site_map.
# These columns are part of std_cols and will carry through to the combined table.
############################################################

msg("======================================================")
msg("Cross-referencing ECS DAPC top SNPs with robust meQTL results...")

if (nrow(ecs_top) > 0 && length(meqtl_results) > 0) {
  ecs_top[, pos := as.integer(pos)]

  for (co in COHORTS) {
    key  <- paste(co, "snp")
    rob  <- meqtl_results[[key]]
    if (is.null(rob) || nrow(rob) == 0) next

    assoc_cols <- intersect(
      c("chr","pos","n_associated_sites","associated_sites",
        "associated_contexts","associated_site_chrs","associated_site_pos"),
      names(rob))
    lk <- unique(rob[, .SD, .SDcols = assoc_cols])
    lk[, pos := as.integer(pos)]

    ecs_co <- ecs_top[cohort == co, .(marker_id, chr, pos)]
    joined <- merge(ecs_co, lk, by = c("chr","pos"), all.x = TRUE)
    matched <- joined[!is.na(n_associated_sites)]

    if (nrow(matched) > 0) {
      msg("  ", co, ": ", nrow(matched), " / ", nrow(ecs_co),
          " DAPC top SNPs are also robust meQTL SNPs")
      # Use set() for in-place assignment by row index (avoids copy overhead)
      for (i in seq_len(nrow(matched))) {
        mid  <- matched$marker_id[i]
        idx  <- which(ecs_top$marker_id == mid & ecs_top$cohort == co)
        if (!length(idx)) next
        set(ecs_top, idx, "n_associated_sites",   matched$n_associated_sites[i])
        set(ecs_top, idx, "associated_sites",      matched$associated_sites[i])
        set(ecs_top, idx, "associated_contexts",   matched$associated_contexts[i])
        if ("associated_site_chrs" %in% names(matched))
          set(ecs_top, idx, "associated_site_chrs", matched$associated_site_chrs[i])
        if ("associated_site_pos"  %in% names(matched))
          set(ecs_top, idx, "associated_site_pos",  matched$associated_site_pos[i])
      }
    } else {
      msg("  ", co, ": no DAPC top SNPs match robust meQTL SNPs by position")
    }
  }

  # Overwrite the saved ECS DAPC annotation file with cross-reference added
  fwrite(ecs_top, file.path(OUT_ROOT, "ecs_dapc_top20_annotated.tsv"), sep = "\t")
  msg("  ECS DAPC annotation updated with meQTL epimarker cross-reference.")

  in_mqtl_n <- if ("n_associated_sites" %in% names(ecs_top))
    ecs_top[!is.na(n_associated_sites), .N] else 0L
  msg("  ECS DAPC SNPs with associated epimarkers: ", in_mqtl_n, " / ", nrow(ecs_top))
} else {
  msg("  Skipping: ecs_top or meqtl_results is empty.")
}

############################################################
# 8) COMBINED TABLE
############################################################

msg("======================================================")
msg("Building combined annotation table...")

# Common schema across all marker sources; missing columns are added as NA
std_cols <- c("marker_id","marker_type","source","tool","cohort","context","DF",
              "chr","pos","loading",
              "methylation_direction","mean_beta","meth_status_absolute",
              "higher_methylation_cohort",
              "n_associated_sites","associated_sites","associated_contexts",
              "associated_site_chrs","associated_site_pos",
              "ref","alt","af","dr2",
              "gene_id","gene_chr","gene_start","gene_end","gene_strand",
              "distance_bp","annotation_class",
              "eggnog_description","pfam_domain","panther_description",
              "go_terms","kegg_ko","is_te_gene")

# Coerce each source table to the shared schema before row-binding
prep_for_combine <- function(dt, extra_id_col = NULL, pos_col = "pos",
                              chr_col = "chr", loading_col = "loading") {
  out <- copy(dt)
  # Ensure all standard columns exist
  for (col in std_cols)
    if (!col %in% names(out)) out[, (col) := NA]
  if (!is.null(extra_id_col) && extra_id_col %in% names(out) &&
      !"marker_id" %in% names(out))
    out[, marker_id := get(extra_id_col)]
  out[, pos := as.integer(get(pos_col))]
  out[, .SD, .SDcols = intersect(std_cols, names(out))]
}

all_list <- list()

if (nrow(ecs_top) > 0)
  all_list[["ECS_DAPC"]] <- prep_for_combine(
    ecs_top[, tool := NA_character_][, loading := loading])

if (nrow(tbs_top) > 0)
  all_list[["TBS_DAPC"]] <- prep_for_combine(tbs_top[, tool := NA_character_])

for (key in names(meqtl_results))
  all_list[[key]] <- prep_for_combine(meqtl_results[[key]])

for (cohort in names(svmp_all_list))
  all_list[[paste("SVMP", cohort)]] <- prep_for_combine(svmp_all_list[[cohort]])

if (length(all_list) > 0) {
  combined <- rbindlist(all_list, fill = TRUE, use.names = TRUE)
  setorder(combined, source, cohort, context, chr, pos, na.last = TRUE)

  fwrite(combined, file.path(OUT_ROOT, "all_markers_annotated.tsv"), sep = "\t")
  msg("Combined table: ", nrow(combined), " rows -> all_markers_annotated.tsv")
  cat("\nAnnotation class breakdown:\n")
  print(combined[, .N, by = .(source, annotation_class)][order(source)])
} else {
  msg("WARNING: no markers were annotated — check that input files exist.")
}

############################################################
# 9) GENE-CENTRIC SUMMARY TABLE
############################################################

msg("======================================================")
msg("Building gene-centric summary table...")

# Re-read the combined table to avoid holding all marker data in memory simultaneously
comb_path <- file.path(OUT_ROOT, "all_markers_annotated.tsv")
if (file.exists(comb_path)) {
  combined_for_genes <- fread(comb_path, header = TRUE)
  gene_summ <- make_gene_summary(combined_for_genes, func_annot, te_genes)

  if (nrow(gene_summ) > 0) {
    out_gs <- file.path(OUT_ROOT, "gene_summary_annotated.tsv")
    fwrite(gene_summ, out_gs, sep = "\t")
    msg("Gene summary: ", nrow(gene_summ), " genes → ", basename(out_gs))
    cat("\nGene summary — annotation class breakdown:\n")
    print(gene_summ[, .N, by = best_annotation_class])
    cat("\nGene summary — source breakdown:\n")
    print(gene_summ[, .(n_genes = .N), by = sources][order(-n_genes)])
  } else {
    msg("  No annotated genes found — gene summary not written.")
  }
} else {
  msg("  Combined table not found — skipping gene summary (run section 8 first).")
}

# Clean up temp files
unlink(TMP_DIR, recursive = TRUE)

msg("======================================================")
msg("Step 17ab finished. Outputs: ", OUT_ROOT)
msg("======================================================")

sessionInfo()
