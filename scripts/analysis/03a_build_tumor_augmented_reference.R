# =============================================================================
# 03a_build_tumor_augmented_reference.R
#
# Purpose
# -------
# Build a PFA-specific tumor-augmented reference for the Semla NNLS
# sensitivity analysis, leaving the original immune-focused reference and all
# previous augmented-reference outputs untouched.
#
# DESIGN
#  1) Start from the original immune-focused reference (exprMatrix + meta).
#  2) Use the tumor-only single-cell dataset (meta_tumor + paired matrix).
#  3) Restrict tumor cells FIRST to patients annotated PFA1 or PFA2.
#  4) Within that pool, treat each tumor `cluster` as a separate NNLS
#     competitor and target N_PER_CLUSTER cells per cluster.
#  5) Sampling is balanced across PFA1/PFA2 when both are available, then
#     across biological samples within subtype.
#  6) Clusters with fewer than MIN_CELLS_PER_GROUP eligible cells are excluded;
#     very small NNLS reference groups are unstable.
#  7) Matrices are gene-matched BY NAME from the start, using only the
#     intersection of genes; original gene order is preserved.
#
# IMPORTANT INTERPRETATION
# The PFA filter uses the patient/tumor `subtype` field, NOT the text of the
# tumor-cluster name. A RELA-like, YAP-like or PFB-like transcriptional cluster
# is therefore retained if those cells come from a PFA1/PFA2 tumor. This keeps
# within-PFA heterogeneity while excluding non-PFA patients.
#
# The script refuses to overwrite existing outputs.
#
# Inputs
# ------
#   ref_file("scrna_expr"), ref_file("scrna_meta")   original immune reference
#   ref_file("tumor_expr"), ref_file("tumor_meta")   tumor-only dataset
#   All four are READ ONLY.
#
# Outputs
# -------
#   REF$augmented_ref_dir/
#     exprMatrix_tumor_augmented_gene_matched.tsv   <- consumed by 03b
#     meta_tumor_augmented_gene_matched.tsv        <- consumed by 03b
#     tumor_subset_meta.tsv, tumor_subset_exprMatrix.tsv,
#     tumor_sampling_by_cluster{,_subtype,_subtype_sample}.tsv,
#     tumor_clusters_excluded.tsv, gene_matching_summary.tsv,
#     reference_build_summary.txt
#
# Stochastic
# ----------
#   YES -- use_seed("tumor_augmented_reference"), currently 1L.
#   Cell sampling per cluster. DO NOT change this seed: the published
#   augmented reference, and every hotspot call derived from it, depends on it.
#   The built reference is also deposited on Zenodo for this reason.
#
# Runtime
# -------
#   ~45 minutes
# =============================================================================

source(here::here("config.R"))

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
  library(tibble)
})


# =============================================================================
# 0. CONFIG
# =============================================================================

# ---- Original immune-focused reference: READ ONLY ----------------------------
original_expr <- ref_file("scrna_expr")
original_meta <- ref_file("scrna_meta")

# ---- Tumor-only dataset: READ ONLY -------------------------------------------
tumor_meta <- ref_file("tumor_meta")
tumor_expr <- ref_file("tumor_expr")

# ---- Output location ----------------------------------------------------------
out_dir <- REF$augmented_ref_dir
OUT_TAG <- basename(out_dir)

out_tumor_meta    <- file.path(out_dir, "tumor_subset_meta.tsv")
out_tumor_expr    <- file.path(out_dir, "tumor_subset_exprMatrix.tsv")
out_sampling      <- file.path(out_dir, "tumor_sampling_by_cluster.tsv")
out_sampling_sub  <- file.path(out_dir, "tumor_sampling_by_cluster_subtype.tsv")
out_sampling_samp <- file.path(out_dir, "tumor_sampling_by_cluster_subtype_sample.tsv")
out_excluded      <- file.path(out_dir, "tumor_clusters_excluded.tsv")

out_aug_expr <- file.path(out_dir, "exprMatrix_tumor_augmented_gene_matched.tsv")
out_aug_meta <- file.path(out_dir, "meta_tumor_augmented_gene_matched.tsv")
out_gene_summary <- file.path(out_dir, "gene_matching_summary.tsv")
out_summary      <- file.path(out_dir, "reference_build_summary.txt")

# ---- Tumor sampling -----------------------------------------------------------
PFA_SUBTYPES        <- c("PFA1", "PFA2")
N_PER_CLUSTER       <- THRESH$cells_per_cluster
MIN_CELLS_PER_GROUP <- 10L
SEED                <- SEEDS$tumor_augmented_reference

# Optional manual exclusions. Default = none.
# Example: EXCLUDE_CLUSTERS <- c("YAP", "PFB")
EXCLUDE_CLUSTERS <- character(0)

# Gene-matching / streaming settings.
BUCKET_SIZE  <- 500L
REPORT_EVERY <- 1000L

# =============================================================================
# 1. SAFETY CHECKS
# =============================================================================

inputs <- c(original_expr, original_meta, tumor_meta, tumor_expr)
missing_inputs <- inputs[!file.exists(inputs)]

if (length(missing_inputs) > 0) {
  stop(
    "Missing required input file(s):\n",
    paste0("  ", missing_inputs, collapse = "\n"),
    "\n\nIf the tumor expression filename differs, update `tumor_expr` in CONFIG."
  )
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

outputs <- c(
  out_tumor_meta, out_tumor_expr,
  out_sampling, out_sampling_sub, out_sampling_samp, out_excluded,
  out_aug_expr, out_aug_meta, out_gene_summary, out_summary
)

existing_outputs <- outputs[file.exists(outputs)]
if (length(existing_outputs) > 0) {
  stop(
    "\nSAFETY STOP: output file(s) already exist. Nothing was overwritten.\n",
    "Use a new OUT_TAG or manually move/delete the previous PFA reference.\n\n",
    paste0("  ", existing_outputs, collapse = "\n")
  )
}

norm <- function(x) normalizePath(x, winslash = "/", mustWork = FALSE)
if (any(norm(outputs) %in% norm(inputs))) {
  stop("SAFETY STOP: an output path resolves to an input path.")
}


# =============================================================================
# 2. HELPERS
# =============================================================================

first_field <- function(x) {
  pos <- regexpr("\t", x, fixed = TRUE)[1]
  if (pos < 0) x else substr(x, 1, pos - 1)
}

after_first_field <- function(x) {
  pos <- regexpr("\t", x, fixed = TRUE)[1]
  if (pos < 0) "" else substr(x, pos + 1, nchar(x))
}

sanitize_group <- function(x) {
  y <- gsub("[^A-Za-z0-9]+", "_", x)
  y <- gsub("^_+|_+$", "", y)
  paste0("Tumor_", y)
}

# Sample n_target rows approximately evenly across a grouping column.
# Used within one subtype to avoid one patient dominating a tumor-state profile.
balanced_sample_rows <- function(df, n_target, group_col = "sample") {
  
  if (nrow(df) <= n_target) return(df)
  
  if (!group_col %in% colnames(df) || all(is.na(df[[group_col]]))) {
    return(dplyr::slice_sample(df, n = n_target))
  }
  
  gid <- as.character(df[[group_col]])
  gid[is.na(gid) | gid == ""] <- "UNKNOWN"
  
  split_idx <- split(seq_len(nrow(df)), gid)
  avail <- lengths(split_idx)
  
  quota <- rep(floor(n_target / length(split_idx)), length(split_idx))
  names(quota) <- names(split_idx)
  quota <- pmin(quota, avail)
  
  remaining <- n_target - sum(quota)
  
  while (remaining > 0L) {
    spare <- avail - quota
    eligible <- which(spare > 0)
    if (length(eligible) == 0L) break
    
    eligible <- eligible[sample.int(length(eligible))]
    for (j in eligible) {
      if (remaining <= 0L) break
      quota[j] <- quota[j] + 1L
      remaining <- remaining - 1L
    }
  }
  
  picked <- unlist(
    Map(
      function(ids, q) {
        q <- as.integer(q)
        if (q <= 0L) integer(0) else ids[sample.int(length(ids), q)]
      },
      split_idx,
      quota
    ),
    use.names = FALSE
  )
  
  df[picked, , drop = FALSE]
}

# Sample one tumor cluster with two-stage balancing:
#   1) balance PFA1/PFA2 as evenly as availability permits;
#   2) within each subtype, balance across patient/sample.
sample_one_cluster <- function(df, n_target) {
  
  if (nrow(df) <= n_target) return(df)
  
  subtype_levels <- PFA_SUBTYPES[PFA_SUBTYPES %in% unique(df$subtype)]
  if (length(subtype_levels) == 0L) return(df[0, , drop = FALSE])
  
  avail <- vapply(
    subtype_levels,
    function(st) sum(df$subtype == st),
    integer(1)
  )
  
  quota <- rep(floor(n_target / length(subtype_levels)), length(subtype_levels))
  names(quota) <- subtype_levels
  quota <- pmin(quota, avail)
  
  remaining <- n_target - sum(quota)
  
  while (remaining > 0L) {
    spare <- avail - quota
    eligible <- which(spare > 0)
    if (length(eligible) == 0L) break
    
    eligible <- eligible[sample.int(length(eligible))]
    for (j in eligible) {
      if (remaining <= 0L) break
      quota[j] <- quota[j] + 1L
      remaining <- remaining - 1L
    }
  }
  
  picked <- lapply(subtype_levels, function(st) {
    d <- df %>% filter(subtype == st)
    q <- as.integer(quota[st])
    balanced_sample_rows(d, q, group_col = "sample")
  })
  
  bind_rows(picked)
}

read_gene_names <- function(path, report_every = REPORT_EVERY) {
  con <- file(path, open = "r")
  on.exit(close(con), add = TRUE)
  
  header <- readLines(con, n = 1)
  if (length(header) != 1L) stop("Could not read header from: ", path)
  
  genes <- character()
  n <- 0L
  
  repeat {
    ln <- readLines(con, n = 1)
    if (length(ln) == 0L) break
    n <- n + 1L
    genes[n] <- first_field(ln)
    if (n %% report_every == 0L) message("  indexed ", n, " genes from ", basename(path))
  }
  
  genes
}


# =============================================================================
# 3. READ NEW TUMOR METADATA AND RESTRICT TO PFA1/PFA2
# =============================================================================

use_seed("tumor_augmented_reference")

message("\nReading NEW tumor metadata...")
tm <- data.table::fread(tumor_meta, data.table = FALSE)

required_tm <- c("Cell", "sample", "subtype", "cluster")
miss <- setdiff(required_tm, colnames(tm))
if (length(miss) > 0L) {
  stop("meta_tumor.tsv is missing required column(s): ", paste(miss, collapse = ", "))
}

tm <- tm %>%
  mutate(
    Cell = as.character(Cell),
    sample = as.character(sample),
    subtype = trimws(as.character(subtype)),
    cluster = trimws(as.character(cluster))
  )

message("\nTumor cells by patient subtype BEFORE PFA restriction:")
print(sort(table(tm$subtype), decreasing = TRUE))

pfa_pool <- tm %>%
  filter(
    subtype %in% PFA_SUBTYPES,
    !is.na(cluster),
    cluster != "",
    !cluster %in% EXCLUDE_CLUSTERS
  )

if (nrow(pfa_pool) == 0L) {
  stop("No tumor cells remain after filtering to PFA1/PFA2.")
}

message("\nPFA-only cells retained: ", nrow(pfa_pool))
message("PFA-only cells by subtype:")
print(sort(table(pfa_pool$subtype), decreasing = TRUE))

cluster_availability <- pfa_pool %>%
  count(cluster, name = "n_PFA_available") %>%
  arrange(desc(n_PFA_available), cluster)

message("\nPFA-only tumor cells by cluster:")
print(as.data.frame(cluster_availability), row.names = FALSE)

excluded_clusters <- cluster_availability %>%
  filter(n_PFA_available < MIN_CELLS_PER_GROUP) %>%
  mutate(reason = paste0("< ", MIN_CELLS_PER_GROUP, " eligible PFA cells"))

included_clusters <- cluster_availability %>%
  filter(n_PFA_available >= MIN_CELLS_PER_GROUP) %>%
  pull(cluster)

if (length(included_clusters) == 0L) {
  stop("No tumor cluster has enough PFA1/PFA2 cells for the augmented reference.")
}

if (nrow(excluded_clusters) > 0L) {
  message("\nExcluded very small tumor clusters:")
  print(as.data.frame(excluded_clusters), row.names = FALSE)
}

write_tsv(excluded_clusters, out_excluded)


# =============================================================================
# 4. SAMPLE UP TO 250 CELLS PER PFA TUMOR CLUSTER
# =============================================================================

message("\nSampling PFA tumor cells...")

selected_list <- lapply(included_clusters, function(cl) {
  
  d <- pfa_pool %>% filter(cluster == cl)
  n_take <- min(N_PER_CLUSTER, nrow(d))
  
  picked <- sample_one_cluster(d, n_take)
  
  if (nrow(picked) != n_take) {
    stop(
      "Sampling error for cluster ", cl,
      ": expected ", n_take, " cells, obtained ", nrow(picked), "."
    )
  }
  
  picked
})

tumor_selected <- bind_rows(selected_list) %>%
  mutate(
    source_tumor_cluster = cluster,
    source_patient_subtype = subtype,
    cell_type = sanitize_group(cluster)
  )

# Safety: sanitised group names must remain one-to-one with source clusters.
group_map <- tumor_selected %>%
  distinct(source_tumor_cluster, cell_type)

if (anyDuplicated(group_map$cell_type)) {
  stop("Two tumor clusters collapsed to the same sanitized NNLS group name.")
}

if (anyDuplicated(tumor_selected$Cell)) {
  stop("Selected tumor cell IDs are not unique.")
}

sampling_cluster <- tumor_selected %>%
  count(source_tumor_cluster, cell_type, name = "n_selected") %>%
  left_join(cluster_availability, by = c("source_tumor_cluster" = "cluster")) %>%
  arrange(source_tumor_cluster)

sampling_subtype <- tumor_selected %>%
  count(source_tumor_cluster, cell_type, source_patient_subtype, name = "n_selected") %>%
  arrange(source_tumor_cluster, source_patient_subtype)

sampling_sample <- tumor_selected %>%
  count(source_tumor_cluster, cell_type, source_patient_subtype, sample, name = "n_selected") %>%
  arrange(source_tumor_cluster, source_patient_subtype, sample)

cat("\n=== SELECTED PFA TUMOR CELLS BY CLUSTER ===\n")
print(as.data.frame(sampling_cluster), row.names = FALSE)

cat("\n=== SELECTED CELLS BY CLUSTER AND PFA SUBTYPE ===\n")
print(as.data.frame(sampling_subtype), row.names = FALSE)

write_tsv(tumor_selected, out_tumor_meta)
write_tsv(sampling_cluster, out_sampling)
write_tsv(sampling_subtype, out_sampling_sub)
write_tsv(sampling_sample, out_sampling_samp)


# =============================================================================
# 5. BUILD AUGMENTED METADATA
# =============================================================================

message("\nBuilding augmented metadata...")
om <- read_tsv(original_meta, show_col_types = FALSE)

if (!all(c("Cell", "cell_type") %in% colnames(om))) {
  stop("Original meta.tsv must contain both 'Cell' and 'cell_type'.")
}

om$Cell <- as.character(om$Cell)

overlap_cells <- intersect(om$Cell, tumor_selected$Cell)
if (length(overlap_cells) > 0L) {
  stop(
    "Selected tumor cells overlap original reference cell IDs. Example: ",
    paste(head(overlap_cells, 5), collapse = ", ")
  )
}

# Match the ORIGINAL metadata column structure exactly.
tumor_for_original_meta <- as.data.frame(
  matrix(
    NA,
    nrow = nrow(tumor_selected),
    ncol = ncol(om),
    dimnames = list(NULL, colnames(om))
  ),
  stringsAsFactors = FALSE
)

shared_cols <- intersect(colnames(om), colnames(tumor_selected))
for (nm in shared_cols) {
  tumor_for_original_meta[[nm]] <- tumor_selected[[nm]]
}

tumor_for_original_meta$Cell <- tumor_selected$Cell
tumor_for_original_meta$cell_type <- tumor_selected$cell_type

aug_meta <- bind_rows(om, tumor_for_original_meta)

if (anyDuplicated(aug_meta$Cell)) {
  stop("Augmented metadata contains duplicated Cell IDs.")
}


# =============================================================================
# 6. EXTRACT ONLY THE SELECTED PFA TUMOR CELLS FROM THE TUMOR MATRIX
# =============================================================================

message("\nLocating selected PFA tumor cells in tumor expression matrix...")

con_tumor_in <- file(tumor_expr, open = "r")
tumor_header_line <- readLines(con_tumor_in, n = 1)

if (length(tumor_header_line) != 1L) {
  close(con_tumor_in)
  stop("Could not read header from tumor expression matrix.")
}

tumor_header <- strsplit(tumor_header_line, "\t", fixed = TRUE)[[1]]
if (length(tumor_header) < 2L) {
  close(con_tumor_in)
  stop("Tumor expression matrix has fewer than 2 columns.")
}

tumor_cells_all <- tumor_header[-1]
selected_cells <- as.character(tumor_selected$Cell)
idx <- match(selected_cells, tumor_cells_all)

if (anyNA(idx)) {
  missing_cells <- selected_cells[is.na(idx)]
  close(con_tumor_in)
  stop(
    length(missing_cells),
    " selected PFA tumor cells are absent from the tumor expression matrix.\nExamples: ",
    paste(head(missing_cells, 10), collapse = ", ")
  )
}

wanted_idx <- c(1L, idx + 1L)
con_tumor_subset <- file(out_tumor_expr, open = "w")

writeLines(
  paste(c(tumor_header[1], selected_cells), collapse = "\t"),
  con_tumor_subset
)

n_tumor_genes <- 0L
repeat {
  ln <- readLines(con_tumor_in, n = 1)
  if (length(ln) == 0L) break
  
  fields <- strsplit(ln, "\t", fixed = TRUE)[[1]]
  if (length(fields) != length(tumor_header)) {
    close(con_tumor_in)
    close(con_tumor_subset)
    stop(
      "Column-count mismatch in tumor expression matrix at gene row ",
      n_tumor_genes + 1L, "."
    )
  }
  
  writeLines(paste(fields[wanted_idx], collapse = "\t"), con_tumor_subset)
  
  n_tumor_genes <- n_tumor_genes + 1L
  if (n_tumor_genes %% REPORT_EVERY == 0L) {
    message("  extracted ", n_tumor_genes, " tumor gene rows...")
  }
}

close(con_tumor_in)
close(con_tumor_subset)

message(
  "Tumor subset matrix written: ", n_tumor_genes,
  " genes x ", length(selected_cells), " selected PFA tumor cells."
)


# =============================================================================
# 7. GENE-MATCH ORIGINAL IMMUNE REFERENCE AND PFA TUMOR SUBSET
# =============================================================================

message("\nReading expression headers...")

con_orig_header <- file(original_expr, open = "r")
orig_header_line <- readLines(con_orig_header, n = 1)
close(con_orig_header)

con_sub_header <- file(out_tumor_expr, open = "r")
sub_header_line <- readLines(con_sub_header, n = 1)
close(con_sub_header)

orig_header <- strsplit(orig_header_line, "\t", fixed = TRUE)[[1]]
sub_header  <- strsplit(sub_header_line, "\t", fixed = TRUE)[[1]]

orig_cells <- orig_header[-1]
sub_cells  <- sub_header[-1]

if (!identical(sub_cells, selected_cells)) {
  stop("Tumor subset expression cell order does not match selected metadata.")
}

if (!setequal(orig_cells, as.character(om$Cell))) {
  stop("Original exprMatrix.tsv cell IDs do not match original meta.tsv cell IDs.")
}

message("\nIndexing original gene names...")
orig_genes <- read_gene_names(original_expr)

message("\nIndexing selected PFA tumor gene names...")
sub_genes <- read_gene_names(out_tumor_expr)

if (anyDuplicated(orig_genes)) {
  dup <- unique(orig_genes[duplicated(orig_genes)])
  stop("Original expression matrix contains duplicated gene names: ",
       paste(head(dup, 10), collapse = ", "))
}

if (anyDuplicated(sub_genes)) {
  dup <- unique(sub_genes[duplicated(sub_genes)])
  stop("Tumor subset expression matrix contains duplicated gene names: ",
       paste(head(dup, 10), collapse = ", "))
}

sub_gene_set <- unique(sub_genes)
common_genes <- orig_genes[orig_genes %in% sub_gene_set]
orig_only <- setdiff(orig_genes, sub_genes)
tumor_only <- setdiff(sub_genes, orig_genes)

if (length(common_genes) == 0L) {
  stop("The original immune and PFA tumor matrices have zero genes in common.")
}

message("\nGene matching:")
message("  original genes : ", length(orig_genes))
message("  tumor genes    : ", length(sub_genes))
message("  common genes   : ", length(common_genes))
message("  original-only  : ", length(orig_only))
message("  tumor-only     : ", length(tumor_only))

orig_pos <- seq_along(orig_genes)
names(orig_pos) <- orig_genes


# =============================================================================
# 8. INDEX TUMOR ROWS INTO SMALL ON-DISK BUCKETS
# =============================================================================

message("\nIndexing tumor-subset rows on disk for gene-matched merge...")

tmp_dir <- file.path(out_dir, paste0(".gene_match_tmp_", Sys.getpid()))
if (dir.exists(tmp_dir)) stop("Temporary directory already exists: ", tmp_dir)
dir.create(tmp_dir, recursive = TRUE)

cleanup_tmp <- TRUE
on.exit({
  if (cleanup_tmp && dir.exists(tmp_dir)) {
    unlink(tmp_dir, recursive = TRUE, force = TRUE)
  }
}, add = TRUE)

n_buckets <- ceiling(length(orig_genes) / BUCKET_SIZE)
bucket_paths <- file.path(tmp_dir, sprintf("bucket_%04d.tsv", seq_len(n_buckets)))
bucket_cons <- lapply(bucket_paths, function(p) file(p, open = "w"))

close_bucket_connections <- function() {
  for (cc in bucket_cons) try(close(cc), silent = TRUE)
}

con_sub <- file(out_tumor_expr, open = "r")
readLines(con_sub, n = 1)  # header

n_sub_rows <- 0L
n_sub_common <- 0L

repeat {
  ln <- readLines(con_sub, n = 1)
  if (length(ln) == 0L) break
  
  n_sub_rows <- n_sub_rows + 1L
  g <- first_field(ln)
  pos <- unname(orig_pos[g])
  
  if (!is.na(pos)) {
    b <- ((pos - 1L) %/% BUCKET_SIZE) + 1L
    writeLines(paste0(pos, "\t", ln), bucket_cons[[b]])
    n_sub_common <- n_sub_common + 1L
  }
  
  if (n_sub_rows %% REPORT_EVERY == 0L) {
    message("  indexed ", n_sub_rows, " tumor rows...")
  }
}

close(con_sub)
close_bucket_connections()

if (n_sub_common != length(common_genes)) {
  stop(
    "Internal matching error: indexed ", n_sub_common,
    " common tumor rows but expected ", length(common_genes), "."
  )
}


# =============================================================================
# 9. STREAM ORIGINAL MATRIX AND WRITE FINAL GENE-MATCHED AUGMENTED MATRIX
# =============================================================================

message("\nWriting final PFA-specific gene-matched augmented expression matrix...")

tmp_out <- paste0(out_aug_expr, ".tmp")
if (file.exists(tmp_out)) unlink(tmp_out)

con_orig <- file(original_expr, open = "r")
con_out  <- file(tmp_out, open = "w")

# Original header + selected PFA tumor cells.
writeLines(paste(c(orig_header, sub_cells), collapse = "\t"), con_out)
readLines(con_orig, n = 1)  # skip original header

current_bucket <- NA_integer_
bucket_payload_by_pos <- NULL

load_bucket <- function(bucket_id) {
  rows <- readLines(bucket_paths[bucket_id])
  if (length(rows) == 0L) return(character())
  
  first_tabs <- regexpr("\t", rows, fixed = TRUE)
  ord <- as.integer(substr(rows, 1, first_tabs - 1L))
  payload <- substr(rows, first_tabs + 1L, nchar(rows))
  
  if (anyDuplicated(ord)) {
    stop("Temporary bucket contains duplicated original row positions.")
  }
  
  names(payload) <- as.character(ord)
  payload
}

n_orig_rows <- 0L
n_written <- 0L

repeat {
  lo <- readLines(con_orig, n = 1)
  if (length(lo) == 0L) break
  
  n_orig_rows <- n_orig_rows + 1L
  g <- first_field(lo)
  
  if (!g %in% sub_gene_set) next
  
  bucket_id <- ((n_orig_rows - 1L) %/% BUCKET_SIZE) + 1L
  if (!identical(bucket_id, current_bucket)) {
    current_bucket <- bucket_id
    bucket_payload_by_pos <- load_bucket(bucket_id)
  }
  
  payload <- unname(bucket_payload_by_pos[as.character(n_orig_rows)])
  if (length(payload) != 1L || is.na(payload)) {
    close(con_orig)
    close(con_out)
    unlink(tmp_out)
    stop("Could not retrieve tumor row for common gene: ", g)
  }
  
  tumor_values <- after_first_field(payload)
  writeLines(paste0(lo, "\t", tumor_values), con_out)
  n_written <- n_written + 1L
  
  if (n_written %% REPORT_EVERY == 0L) {
    message("  wrote ", n_written, " common genes...")
  }
}

close(con_orig)
close(con_out)

if (n_written != length(common_genes)) {
  unlink(tmp_out)
  stop(
    "Final augmented matrix contains ", n_written,
    " genes but expected ", length(common_genes), "."
  )
}

if (!file.rename(tmp_out, out_aug_expr)) {
  unlink(tmp_out)
  stop("Could not rename temporary augmented matrix to final output.")
}

cleanup_tmp <- FALSE
unlink(tmp_dir, recursive = TRUE, force = TRUE)


# =============================================================================
# 10. WRITE FINAL METADATA IN EXPRESSION-MATRIX CELL ORDER
# =============================================================================

final_cells <- c(orig_cells, sub_cells)

if (!setequal(final_cells, aug_meta$Cell)) {
  stop("Final expression and augmented metadata cell IDs do not match.")
}

aug_meta_out <- aug_meta[match(final_cells, aug_meta$Cell), , drop = FALSE]

if (anyNA(aug_meta_out$Cell)) {
  stop("Unexpected NA while reordering augmented metadata.")
}

write_tsv(aug_meta_out, out_aug_meta)


# =============================================================================
# 11. SUMMARIES / QC
# =============================================================================

gene_summary <- tibble(
  metric = c(
    "original_genes",
    "tumor_genes",
    "common_genes_retained",
    "original_only_removed",
    "tumor_only_removed",
    "percent_original_genes_retained"
  ),
  value = c(
    length(orig_genes),
    length(sub_genes),
    length(common_genes),
    length(orig_only),
    length(tumor_only),
    100 * length(common_genes) / length(orig_genes)
  )
)

write_tsv(gene_summary, out_gene_summary)

summary_lines <- c(
  "PFA-specific tumor-augmented Semla reference",
  "=============================================",
  "",
  paste0("Created: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("PFA patient subtypes retained: ", paste(PFA_SUBTYPES, collapse = ", ")),
  paste0("Target cells per tumor cluster: ", N_PER_CLUSTER),
  paste0("Minimum cells required per tumor NNLS group: ", MIN_CELLS_PER_GROUP),
  paste0("Random seed: ", SEED),
  "",
  "Sampling design:",
  "  - filter to PFA1/PFA2 patients first",
  "  - sample separately within each tumor cluster",
  "  - balance PFA1/PFA2 when both are represented",
  "  - balance across patient/sample within subtype",
  "  - use all eligible cells when fewer than target are available",
  "",
  paste0("Eligible PFA tumor cells before sampling: ", nrow(pfa_pool)),
  paste0("Added PFA tumor cells: ", nrow(tumor_selected)),
  paste0("Included tumor clusters: ", length(included_clusters)),
  paste0("Excluded tumor clusters (<", MIN_CELLS_PER_GROUP, " cells): ", nrow(excluded_clusters)),
  "",
  paste0("Original immune reference cells: ", length(orig_cells)),
  paste0("Final augmented reference cells: ", length(final_cells)),
  paste0("Common genes retained: ", length(common_genes), " / ", length(orig_genes),
         " (", sprintf("%.2f", 100 * length(common_genes) / length(orig_genes)), "%)"),
  "",
  "Added tumor NNLS groups:",
  paste0("  ", sampling_cluster$cell_type, ": ", sampling_cluster$n_selected,
         " cells (", sampling_cluster$n_PFA_available, " PFA cells available)"),
  "",
  "Files for 03b:",
  paste0("  sn_exp_path  <- \"", gsub("\\\\", "/", out_aug_expr), "\""),
  paste0("  sn_meta_path <- \"", gsub("\\\\", "/", out_aug_meta), "\""),
  "",
  "Recommended NEW Semla output root:",
  "  Analysis/02-4_semla_nnls_tumor_augmented_PFA"
)

writeLines(summary_lines, out_summary)

cat("\n============================================================\n")
cat("DONE — PFA-SPECIFIC AUGMENTED REFERENCE CREATED\n")
cat("============================================================\n")
cat("Output folder:\n  ", out_dir, "\n\n", sep = "")
cat("PFA tumor cells added: ", nrow(tumor_selected), "\n", sep = "")
cat("Tumor groups added: ", length(included_clusters), "\n", sep = "")
cat("Common genes retained: ", length(common_genes), " / ", length(orig_genes), "\n", sep = "")
cat("\nUse in 03b:\n")
cat('  sn_exp_path  <- "', out_aug_expr, '"\n', sep = "")
cat('  sn_meta_path <- "', out_aug_meta, '"\n', sep = "")
cat("============================================================\n")
