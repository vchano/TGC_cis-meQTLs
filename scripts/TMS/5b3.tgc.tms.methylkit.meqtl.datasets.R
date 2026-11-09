#!/usr/bin/env Rscript
############################################################
# TreeGeneClimate (TGC) — TMS
# Step 5b3: Methylation datasets for cis-meQTL mapping (one cohort x context)
#
# Larger site sets than Step 5b2: no upper coverage limit and no
# methylated-fraction filter, and sites are united family by family
# (breeding) or stand by stand (natural), so that every group is represented.
#
#   1) Remove sites with an ECS SNP at the cytosine or at the bases defining
#      its context (as in Step 5b2)
#   2) Per-sample coverage >= 10 reads (applied at import, Step 5b1)
#   3) unite() with families (breeding) or stands (natural) as groups:
#      a site is kept if it is covered in at least m samples of EVERY group,
#      where m is the size of the smallest group (6 trees for breeding
#      families, 2 for natural stands); strand-specific (destrand = FALSE)
#
# USAGE
#   Rscript 5b3.tgc.tms.methylkit.meqtl.datasets.R <BREEDING|NATURAL> <CpG|CHG|CHH>
#
# INPUTS
#   RESULTS/TMS/METHYLKIT_RAW/<SAMPLE>_<CONTEXT>.rds              (Step 5b1)
#   RESULTS/ECS/VARIANT.CALLING/MASK/tgc.ecs.cytosine.mask.snps.tsv.gz (Step 4a3)
#   DATA/METADATA/breeding_sample2family.txt, natural_sample2pop.txt
#
# OUTPUTS (RESULTS/TMS/RANALYSIS/METHYLKIT_OBJECTS_MEQTL/)
#   methylBase_<cohort>_<ctx>_meqtl.rds   (attribute "groups": sample, group)
#   summary_<cohort>_<ctx>_meqtl.tsv
############################################################

# === USER CONFIGURATION ===
PROJECT_ROOT <- Sys.getenv("TGC_PROJECT_ROOT", unset = "/path/to/your/project")
# ===========================

source(file.path(PROJECT_ROOT, "SCRIPTS/TMS/tgc.tms.functions.R"))

args <- commandArgs(trailingOnly = TRUE)
COHORT <- toupper(args[1]); CTX <- args[2]
stopifnot(COHORT %in% c("BREEDING", "NATURAL"), CTX %in% c("CpG", "CHG", "CHH"))

RAW_DIR   <- file.path(PROJECT_ROOT, "RESULTS/TMS/METHYLKIT_RAW")
MASK_FILE <- file.path(PROJECT_ROOT, "RESULTS/ECS/VARIANT.CALLING/MASK/tgc.ecs.cytosine.mask.snps.tsv.gz")
MAP_FILE  <- file.path(PROJECT_ROOT, "DATA/METADATA",
                       if (COHORT == "BREEDING") "breeding_sample2family.txt" else "natural_sample2pop.txt")
OUT_DIR   <- file.path(PROJECT_ROOT, "RESULTS/TMS/RANALYSIS/METHYLKIT_OBJECTS_MEQTL")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
msg <- function(...) cat("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n", sep = "")

files <- raw_files(RAW_DIR, COHORT, CTX)
groups <- fread(MAP_FILE, header = FALSE, col.names = c("sample", "group"))
groups <- groups[match(names(files), sample)]
if (anyNA(groups$group)) stop("Samples without group: ", paste(names(files)[is.na(groups$group)], collapse = ", "))
group_index <- match(groups$group, sort(unique(groups$group)))
min_group <- min(table(group_index))
msg(COHORT, " / ", CTX, ": ", length(files), " samples in ", length(unique(group_index)),
    " groups; smallest group = ", min_group)

mask <- read_snp_mask(MASK_FILE)
n_raw <- n_masked <- integer(length(files))
objs <- vector("list", length(files))
for (i in seq_along(files)) {
  mr <- readRDS(files[i])
  n_raw[i] <- nrow(mr)
  mk <- mask_methylRaw(mr, mask, CTX)
  objs[[i]] <- mk$object
  n_masked[i] <- mk$n_removed
}
mrl <- methylRawList(objs, treatment = group_index)
rm(objs); gc()

mb <- unite(mrl, destrand = FALSE, min.per.group = as.integer(min_group))
rm(mrl); gc()
attr(mb, "groups") <- data.frame(sample = names(files), group = groups$group)
saveRDS(mb, file.path(OUT_DIR, sprintf("methylBase_%s_%s_meqtl.rds", tolower(COHORT), tolower(CTX))))
msg("Sites after unite: ", nrow(mb))

summ <- data.frame(
  cohort = COHORT, context = CTX, n_samples = length(files),
  n_groups = length(unique(group_index)), min_per_group = min_group,
  raw_sites_median = median(n_raw), masked_sites_median = median(n_masked),
  n_sites_after_unite = nrow(mb)
)
write.table(summ, file.path(OUT_DIR, sprintf("summary_%s_%s_meqtl.tsv", tolower(COHORT), tolower(CTX))),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(summ)
msg("Step 5b3 done")
