# =============================================================================
# 05_hotspot_pseudobulk_CIBERSORTx.R
#
# Purpose
# -------
# Generate section-specific Visium pseudo-bulk expression profiles for external
# CIBERSORTx deconvolution with the LM7 immune signature matrix.
#
# Design
# ------
# For each retained section, raw Visium counts are summed separately across:
#   1) Semla-defined lymphocyte hotspot spots
#   2) all remaining in-tissue non-hotspot spots (background)
#
# This produces paired hotspot/background pseudo-bulks per section. The script
# does NOT perform deconvolution itself; CIBERSORTx is run externally.
#
# Inputs
# ------
#   stage_input("hotspots", "vis/vis_section_<sec>.rds")   for each GOOD_SECTION
#     Written by 02c_add_hotspots_to_coordinate.R. Carries raw Spatial counts
#     plus is_hotspot_semla / local_lymph_score / in_tissue in meta.data.
#
# Outputs
# -------
#   stage_dir("pseudobulk", "pseudobulk_counts.tsv")
#       Raw summed counts, genes x pseudo-bulk samples. Traceability, and for
#       any downstream method requiring count input.
#   stage_dir("pseudobulk", "pseudobulk_cpm_CIBERSORTx.tsv")
#       CPM-normalized, non-log mixture matrix for upload to CIBERSORTx.
#   stage_dir("pseudobulk", "pseudobulk_sample_metadata.tsv")
#       Per pseudo-bulk: section, group, contributing spots, raw library size,
#       CPM column sum, and the number of LM7 genes detected.
#
# Stochastic
# ----------
#   none
#
# Runtime
# -------
#   ~2 minutes
#
# Note on stage 04
# ----------------
# This script previously read 04_lm7_deconvolution/vis_good_lm7_ready.rds. That
# stage implemented the per-spot UCell approach, which has been superseded by
# the CIBERSORTx route and no longer exists. Everything stage 04 contributed
# was UCell scoring, which is unused here, so this script now reads the 02c
# outputs directly and stage 04 is removed from the pipeline.
#
# Important interpretation
# ------------------------
# CIBERSORTx/LM7 relative fractions describe the composition of the modeled
# leukocyte compartment. They are NOT estimates of the fraction of all cells in
# the spatial region that are immune cells.
#
# CIBERSORTx run used for the manuscript (see CIBERSORTX in config.R):
#   signature matrix LM7; B-mode batch correction; QN off; relative; 500 perms.
# =============================================================================

source(here::here("config.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(tidyverse)
})

MIN_SPOTS <- THRESH$min_spots_pseudobulk

# Optional: path to the LM7 matrix, used only to report how many signature
# genes are detected per pseudo-bulk. Not redistributed with this repository;
# see docs/data_sources.md. Set to NA to skip the check.
LM7_PATH <- Sys.getenv("EPN_LM7_PATH", unset = NA_character_)

# -----------------------------------------------------------------------------
# 1. Helpers
# -----------------------------------------------------------------------------
get_counts <- function(obj, assay = "Spatial") {
  a <- obj[[assay]]

  if (inherits(a, "Assay5")) {
    lyrs <- SeuratObject::Layers(a, search = "counts")
    if (length(lyrs) == 0) stop("No counts layer in assay: ", assay)
    if (length(lyrs) > 1) {
      message("  joining ", length(lyrs), " count layers")
      obj[[assay]] <- SeuratObject::JoinLayers(a)
      a <- obj[[assay]]
    }
    return(SeuratObject::LayerData(a, layer = "counts"))
  }

  SeuratObject::GetAssayData(a, slot = "counts")
}

# Sum a section's counts into hotspot and background column vectors.
section_pseudobulk <- function(sec) {
  f <- stage_input("hotspots", file.path("vis", paste0("vis_section_", sec, ".rds")))
  message("Section ", sec)
  obj <- readRDS(f)
  DefaultAssay(obj) <- "Spatial"

  md <- obj@meta.data
  required <- "is_hotspot_semla"
  if (!required %in% colnames(md)) {
    stop("Section ", sec, ": missing metadata column 'is_hotspot_semla'. ",
         "Was 02c_add_hotspots_to_coordinate.R run for this section?")
  }

  in_tiss <- if ("in_tissue" %in% colnames(md)) as.integer(md$in_tissue) else 1L
  hotspot <- as.character(md$is_hotspot_semla) == "TRUE"
  hotspot[is.na(hotspot)] <- FALSE

  keep <- which(in_tiss == 1L)
  if (length(keep) == 0) stop("Section ", sec, ": no in-tissue spots.")

  grp <- ifelse(hotspot[keep], "hotspot", "background")
  n_hot <- sum(grp == "hotspot")
  n_bg  <- sum(grp == "background")

  if (n_hot < MIN_SPOTS || n_bg < MIN_SPOTS) {
    stop("Section ", sec, " has ", n_hot, " hotspot and ", n_bg,
         " background in-tissue spots; at least ", MIN_SPOTS,
         " are required in BOTH groups. The paired design would be incomplete, ",
         "so no pseudo-bulk files were written.")
  }

  counts <- get_counts(obj, "Spatial")
  if (is.null(counts) || nrow(counts) == 0) {
    stop("Section ", sec, ": no usable raw counts in the Spatial assay.")
  }
  counts <- counts[, keep, drop = FALSE]

  ind <- Matrix::sparse.model.matrix(~ 0 + factor(grp, levels = c("hotspot", "background")))
  pb  <- as.matrix(counts %*% ind)
  colnames(pb) <- paste0(sec, "_", c("hotspot", "background"))

  list(
    pb = pb,
    n  = tibble::tibble(
      section = sec,
      group   = c("hotspot", "background"),
      sample  = colnames(pb),
      n_spots = c(n_hot, n_bg)
    )
  )
}

# -----------------------------------------------------------------------------
# 2. Build pseudo-bulks for every retained section
# -----------------------------------------------------------------------------
parts <- lapply(GOOD_SECTIONS, section_pseudobulk)

n_by_sample <- dplyr::bind_rows(lapply(parts, `[[`, "n"))

cat("\n=== Spots per pseudo-bulk ===\n")
print(n_by_sample, n = Inf)

# Align on the union of genes. Sections share a reference so this is normally
# a no-op, but aligning explicitly means a mismatch cannot silently truncate.
all_genes <- sort(unique(unlist(lapply(parts, function(p) rownames(p$pb)))))
gene_sets <- lapply(parts, function(p) rownames(p$pb))
if (length(unique(vapply(gene_sets, length, integer(1)))) > 1) {
  warning("Sections do not share an identical gene set; aligning on the union ",
          "and filling absent genes with zero.")
}

pb <- matrix(
  0,
  nrow = length(all_genes),
  ncol = sum(vapply(parts, function(p) ncol(p$pb), integer(1))),
  dimnames = list(all_genes, unlist(lapply(parts, function(p) colnames(p$pb))))
)
for (p in parts) pb[rownames(p$pb), colnames(p$pb)] <- p$pb

# Order columns section by section, hotspot before background.
sample_levels <- unlist(lapply(GOOD_SECTIONS, function(s) {
  c(paste0(s, "_hotspot"), paste0(s, "_background"))
}))
pb <- pb[, sample_levels[sample_levels %in% colnames(pb)], drop = FALSE]

# -----------------------------------------------------------------------------
# 3. Gene harmonization
# -----------------------------------------------------------------------------
pb <- pb[rowSums(pb) > 0, , drop = FALSE]

# LM7 gene symbols are uppercase throughout, so upper-casing here is the
# correct harmonization step. It can create duplicates, which are summed.
genes <- toupper(trimws(rownames(pb)))
valid_gene <- !is.na(genes) & genes != ""
pb <- pb[valid_gene, , drop = FALSE]
rownames(pb) <- genes[valid_gene]

if (anyDuplicated(rownames(pb)) > 0) {
  message("Collapsing duplicated gene symbols created by harmonization")
  pb <- rowsum(pb, group = rownames(pb), reorder = FALSE)
}

if (any(colSums(pb) == 0)) {
  stop("At least one pseudo-bulk has zero total counts after processing.")
}

message("Pseudo-bulk matrix: ", nrow(pb), " genes x ", ncol(pb), " samples")
cat("\n=== Raw pseudo-bulk library sizes (millions of counts) ===\n")
print(round(colSums(pb) / 1e6, 3))

# -----------------------------------------------------------------------------
# 4. CPM mixture for CIBERSORTx
# -----------------------------------------------------------------------------
# CPM preserves linear (non-log) expression while correcting for the very
# different total depths of hotspot and background pseudo-bulks.
cpm <- sweep(pb, 2, colSums(pb), "/") * 1e6

counts_out <- stage_dir("pseudobulk", "pseudobulk_counts.tsv")
write.table(
  data.frame(Gene = rownames(pb), pb, check.names = FALSE),
  counts_out, sep = "\t", quote = FALSE, row.names = FALSE
)

cpm_out <- stage_dir("pseudobulk", "pseudobulk_cpm_CIBERSORTx.tsv")
write.table(
  data.frame(Gene = rownames(cpm), cpm, check.names = FALSE),
  cpm_out, sep = "\t", quote = FALSE, row.names = FALSE
)

# -----------------------------------------------------------------------------
# 5. Sample metadata and LM7 coverage
# -----------------------------------------------------------------------------
# How many LM7 signature genes are actually present, and non-zero per sample.
# A low count is a better explanation for populations deconvolving to exactly
# zero than any biological reading, so it is recorded alongside the mixture.
lm7_genes <- NULL
if (!is.na(LM7_PATH) && file.exists(LM7_PATH)) {
  lm7 <- read.delim(LM7_PATH, header = TRUE, row.names = 1, check.names = FALSE)
  lm7_genes <- toupper(trimws(rownames(lm7)))
  message("LM7 genes: ", length(lm7_genes),
          "; present in mixture: ", sum(lm7_genes %in% rownames(pb)))
} else {
  message("LM7 matrix not provided (set EPN_LM7_PATH); skipping coverage check.")
}

sample_qc <- n_by_sample %>%
  filter(sample %in% colnames(pb)) %>%
  mutate(
    sample                   = factor(sample, levels = colnames(pb)),
    raw_library_size         = as.numeric(colSums(pb)[as.character(sample)]),
    raw_library_size_million = raw_library_size / 1e6,
    cpm_sum                  = as.numeric(colSums(cpm)[as.character(sample)]),
    n_lm7_genes_nonzero      = if (is.null(lm7_genes)) NA_integer_ else {
      shared <- intersect(lm7_genes, rownames(pb))
      as.integer(colSums(pb[shared, as.character(sample), drop = FALSE] > 0))
    }
  ) %>%
  arrange(sample) %>%
  mutate(sample = as.character(sample))

metadata_out <- stage_dir("pseudobulk", "pseudobulk_sample_metadata.tsv")
readr::write_tsv(sample_qc, metadata_out)

cat("\n=== Pseudo-bulk metadata ===\n")
print(sample_qc, n = Inf)

message("\nDone.")
message("Raw counts:       ", counts_out)
message("CIBERSORTx input: ", cpm_out)
message("Sample metadata:  ", metadata_out)
message("\nRun CIBERSORTx externally with the LM7 signature matrix.")
message("Parameters: ", paste(names(CIBERSORTX), unlist(CIBERSORTX),
                              sep = "=", collapse = "; "))
