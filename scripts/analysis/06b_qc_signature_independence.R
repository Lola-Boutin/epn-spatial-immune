# =============================================================================
# 06b_qc_signature_independence.R
#
# Purpose
# -------
# QC: are the zone signatures independent of the Lymphocyte NNLS signal?
#
# Compares the gene sets used for Mesenchymal and Myeloid zone scoring in
# stage 06 (AddModuleScore, fixed curated lists) against the genes that best
# distinguish Lymphocytes from all other cell types in the reference used for
# NNLS in stage 02 (FindMarkers, Lymphocytes vs all).
#
# Why it matters: zone scores use a fixed gene list, while Lymphocyte NNLS fits
# the full transcriptome as a mixture of reference pseudo-bulks. If the two
# gene sets are largely non-overlapping, spatial co-localisation of high
# Lymphocyte NNLS with Mesenchymal / Myeloid zones reflects genuine biological
# co-occurrence rather than a shared gene-measurement artefact.
#
# The settings below MUST match stage 06 exactly.
#
# Inputs
# ------
#   ref_file("zone_genes_xlsx")
#   ref_file("scrna_expr"), ref_file("scrna_meta")
#
# Outputs
# -------
#   stage_dir("signature_overlap")
#     signature_gene_lists.tsv, overlap_summary.tsv,
#     shared_genes_meso_lymph.txt, shared_genes_myeloid_lymph.txt,
#     shared_genes_meso_myeloid.txt,
#     barplot_signature_sizes.png, upset_overlap.png,
#     lymphocyte_markers_full.tsv
#
# Stochastic
# ----------
#   none
#
# Runtime
# -------
#   ~20 minutes
# =============================================================================

source(here::here("config.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(readxl)
  library(dplyr)
  library(tibble)
  library(readr)
  library(purrr)
  library(ggplot2)
  library(tidyr)
})


# -------------------------
# Configuration
# -------------------------
xlsx         <- ref_file("zone_genes_xlsx")
sn_exp_path  <- ref_file("scrna_expr")
sn_meta_path <- ref_file("scrna_meta")

out_dir <- stage_dir("signature_overlap")

# -------------------------
# STAGE 06 SETTINGS  (must match 06 exactly)
# -------------------------
TOP_N_PER_CLUSTER <- 100
MIN_DELTA_PCT     <- 0.15
MIN_PCT1          <- 0.20
MIN_LOG2FC        <- 0.40
MAX_PCT2          <- 0.30
USE_METRIC        <- "score"

# Number of top Lymphocyte marker genes used for comparison. NNLS uses the full
# transcriptome, so the top N most discriminative genes stand in as a
# representative set.
TOP_LYMPH_MARKERS <- 200

# =========================================================
# PART 1 — Build Mesenchymal + Myeloid gene sets
#          (identical logic to build_zone_gene_sets in phase06)
# =========================================================
message("Building zone gene sets from Excel...")

zone_sheets <- list(
  Mesenchymal = c("MEC-A", "MEC-B", "MEC-C", "MEC-D", "UEC-B"),
  Myeloid     = c("classic-M", "hypoxia-M", "chemokine-M")
)

endo_markers <- toupper(c(
  "PECAM1","VWF","KDR","FLT1","TEK","ENG","EMCN","RAMP2","PLVAP",
  "CD34","ESAM","ROBO4","PGF","KLF2","KLF4","CLEC14A","FABP4",
  "MCAM","SOX17","EDNRB","ADGRL4","PROM1"
))

available_sheets <- excel_sheets(xlsx)
missing_sheets   <- setdiff(unlist(zone_sheets), available_sheets)
if (length(missing_sheets) > 0) {
  stop("Missing Excel sheets: ", paste(missing_sheets, collapse = ", "))
}

read_one_sheet <- function(sheet, zone_name) {
  read_excel(xlsx, sheet = sheet) %>%
    transmute(
      zone       = zone_name,
      cluster    = sheet,
      gene       = as.character(gene),
      p_val_adj  = suppressWarnings(as.numeric(p_val_adj)),
      avg_log2fc = suppressWarnings(as.numeric(avg_log2FC)),
      pct_1      = suppressWarnings(as.numeric(`pct.1`)),
      pct_2      = suppressWarnings(as.numeric(`pct.2`)),
      delta_pct  = pct_1 - pct_2,
      score      = avg_log2fc * (pct_1 - pct_2)
    ) %>%
    filter(!is.na(gene), gene != "")
}

markers_raw <- imap_dfr(zone_sheets, function(sheets, zone_name) {
  map_dfr(sheets, function(sh) read_one_sheet(sh, zone_name))
})

markers_clean <- markers_raw %>%
  mutate(
    p_val_adj  = ifelse(is.na(p_val_adj),  1, p_val_adj),
    avg_log2fc = ifelse(is.na(avg_log2fc), 0, avg_log2fc),
    pct_1      = ifelse(is.na(pct_1),      0, pct_1),
    pct_2      = ifelse(is.na(pct_2),      0, pct_2),
    delta_pct  = pct_1 - pct_2,
    score      = avg_log2fc * delta_pct
  ) %>%
  filter(delta_pct >= MIN_DELTA_PCT) %>%
  filter(zone != "Vascular" |
           (pct_1 >= MIN_PCT1 & avg_log2fc >= MIN_LOG2FC & pct_2 <= MAX_PCT2))

metric_sym <- rlang::sym(USE_METRIC)

top_per_cluster <- markers_clean %>%
  group_by(zone, cluster) %>%
  arrange(p_val_adj, desc(!!metric_sym), desc(avg_log2fc), desc(delta_pct)) %>%
  slice_head(n = TOP_N_PER_CLUSTER) %>%
  ungroup()

zone_gene_sets <- top_per_cluster %>%
  group_by(zone) %>%
  summarise(genes = list(sort(unique(toupper(gene)))), .groups = "drop")

zone_genes <- setNames(zone_gene_sets$genes, zone_gene_sets$zone)

genes_mesenchymal <- zone_genes[["Mesenchymal"]]
genes_myeloid     <- zone_genes[["Myeloid"]]

message("  Mesenchymal signature: ", length(genes_mesenchymal), " genes")
message("  Myeloid signature:     ", length(genes_myeloid), " genes")

# =========================================================
# PART 2 — Build Lymphocyte marker gene set from snRNA-seq
#          reference (FindMarkers: Lymphocytes vs all others)
#
# This represents the genes most informative for the NNLS
# Lymphocyte coefficient; the full reference transcriptome
# is used by RunNNLS, so we take the top discriminating genes
# as the effective "signature".
# =========================================================
message("Loading snRNA-seq reference...")

exp  <- read.delim(sn_exp_path, header = TRUE, row.names = 1, check.names = FALSE)
sn   <- CreateSeuratObject(counts = as.matrix(exp), min.cells = 3, min.features = 200)
meta <- read.delim(sn_meta_path, header = TRUE, stringsAsFactors = FALSE)
rownames(meta) <- meta$Cell
meta <- meta[colnames(sn), , drop = FALSE]
sn   <- AddMetaData(sn, metadata = meta)

if (!"cell_type" %in% colnames(sn@meta.data)) {
  stop("snRNA metadata does not contain 'cell_type'. Check sn_meta_path.")
}

if (!"Lymphocytes" %in% unique(sn$cell_type)) {
  available_types <- sort(unique(sn$cell_type))
  stop(
    "No cells labelled 'Lymphocytes' found in cell_type column.\n",
    "Available types: ", paste(available_types, collapse = ", ")
  )
}

# Normalize before FindMarkers
sn <- NormalizeData(sn, verbose = FALSE)
Idents(sn) <- sn$cell_type

message("Running FindMarkers: Lymphocytes vs all other cell types...")
lymph_markers <- FindMarkers(
  sn,
  ident.1  = "Lymphocytes",
  min.pct  = 0.10,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  verbose  = FALSE
) %>%
  rownames_to_column("gene") %>%
  mutate(gene = toupper(gene)) %>%
  arrange(p_val_adj, desc(avg_log2FC))

write_tsv(lymph_markers, file.path(out_dir, "lymphocyte_markers_full.tsv"))
message("  Total Lymphocyte markers (p_adj < 0.05, positive): ",
        sum(lymph_markers$p_val_adj < 0.05 & lymph_markers$avg_log2FC > 0, na.rm = TRUE))

# Use only upregulated markers (positive log2FC) as the lymphocyte "signature"
# — these are the genes that drive a high Lymphocyte NNLS score
genes_lymphocyte <- lymph_markers %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0) %>%
  slice_head(n = TOP_LYMPH_MARKERS) %>%
  pull(gene)

message("  Lymphocyte signature (top ", TOP_LYMPH_MARKERS,
        " upregulated markers): ", length(genes_lymphocyte), " genes")

# =========================================================
# PART 3 — Pairwise overlap analysis
# =========================================================
message("Computing pairwise overlaps...")

jaccard <- function(a, b) {
  a <- unique(a); b <- unique(b)
  n_inter <- length(intersect(a, b))
  n_union <- length(union(a, b))
  if (n_union == 0) return(NA_real_)
  n_inter / n_union
}

gene_lists <- list(
  Mesenchymal = genes_mesenchymal,
  Myeloid     = genes_myeloid,
  Lymphocyte  = genes_lymphocyte
)

pair_names <- combn(names(gene_lists), 2, simplify = FALSE)

overlap_summary <- map_dfr(pair_names, function(pair) {
  a_name <- pair[1]; b_name <- pair[2]
  a <- gene_lists[[a_name]]; b <- gene_lists[[b_name]]
  shared <- intersect(a, b)
  tibble(
    set_A          = a_name,
    set_B          = b_name,
    n_A            = length(a),
    n_B            = length(b),
    n_shared       = length(shared),
    pct_of_A       = round(100 * length(shared) / length(a), 1),
    pct_of_B       = round(100 * length(shared) / length(b), 1),
    jaccard_index  = round(jaccard(a, b), 4),
    shared_genes   = paste(sort(shared), collapse = ", ")
  )
})

print(overlap_summary %>% select(-shared_genes))
write_tsv(overlap_summary, file.path(out_dir, "overlap_summary.tsv"))

# Save shared gene lists as separate files for easy inspection
for (i in seq_len(nrow(overlap_summary))) {
  row    <- overlap_summary[i, ]
  shared <- strsplit(row$shared_genes, ", ")[[1]]
  shared <- shared[nzchar(shared)]
  fname  <- paste0("shared_genes_",
                   tolower(row$set_A), "_",
                   tolower(row$set_B), ".txt")
  writeLines(
    c(paste0("# Shared genes between ", row$set_A, " and ", row$set_B,
             "  (n=", length(shared), ", Jaccard=", row$jaccard_index, ")"),
      sort(shared)),
    file.path(out_dir, fname)
  )
  message("  ", row$set_A, " ∩ ", row$set_B, " = ", length(shared),
          " genes (Jaccard = ", row$jaccard_index, ")")
}

# Save combined gene list table (one row per gene, membership columns)
all_genes <- sort(unique(c(genes_mesenchymal, genes_myeloid, genes_lymphocyte)))
gene_table <- tibble(
  gene         = all_genes,
  Mesenchymal  = gene %in% genes_mesenchymal,
  Myeloid      = gene %in% genes_myeloid,
  Lymphocyte   = gene %in% genes_lymphocyte,
  n_sets       = Mesenchymal + Myeloid + Lymphocyte
) %>%
  arrange(desc(n_sets), gene)

write_tsv(gene_table, file.path(out_dir, "signature_gene_lists.tsv"))
message("  Total unique genes across all three sets: ", nrow(gene_table))
message("  Genes in 2+ sets: ", sum(gene_table$n_sets >= 2))
message("  Genes in all 3:   ", sum(gene_table$n_sets == 3))

# =========================================================
# PART 4 — Plots
# =========================================================

# 4a. Bar chart: signature sizes
size_df <- tibble(
  Signature = names(gene_lists),
  n_genes   = lengths(gene_lists)
) %>%
  mutate(Signature = factor(Signature, levels = c("Mesenchymal", "Myeloid", "Lymphocyte")))

fill_cols <- c(
  "Mesenchymal" = "#E41A1C",
  "Myeloid"     = "#984EA3",
  "Lymphocyte"  = "#377EB8"
)

p_bar <- ggplot(size_df, aes(x = Signature, y = n_genes, fill = Signature)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = n_genes), vjust = -0.4, size = 4.5, fontface = "bold") +
  scale_fill_manual(values = fill_cols) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "Signature gene set sizes",
    x     = NULL,
    y     = "Number of genes"
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = "none",
        axis.title.y    = element_text(face = "bold"),
        plot.title      = element_text(face = "bold", hjust = 0.5))

ggsave(file.path(out_dir, "barplot_signature_sizes.png"),
       p_bar, width = 5.5, height = 4.5, dpi = 300, bg = "white")

# 4b. Overlap dot-matrix (UpSet-style without extra packages)
# Shows per-combination membership counts
combos <- list(
  "Mesenchymal only"        = c(TRUE,  FALSE, FALSE),
  "Myeloid only"            = c(FALSE, TRUE,  FALSE),
  "Lymphocyte only"         = c(FALSE, FALSE, TRUE),
  "Meso ∩ Myeloid"          = c(TRUE,  TRUE,  FALSE),
  "Meso ∩ Lymphocyte"       = c(TRUE,  FALSE, TRUE),
  "Myeloid ∩ Lymphocyte"    = c(FALSE, TRUE,  TRUE),
  "All three"               = c(TRUE,  TRUE,  TRUE)
)

combo_counts <- map_dfr(names(combos), function(nm) {
  mask <- combos[[nm]]
  n <- sum(
    gene_table$Mesenchymal == mask[1] &
      gene_table$Myeloid    == mask[2] &
      gene_table$Lymphocyte == mask[3]
  )
  tibble(combo = nm, n = n,
         in_meso  = mask[1],
         in_myelo = mask[2],
         in_lymph = mask[3])
}) %>%
  mutate(combo = factor(combo, levels = rev(names(combos))))

# dot membership panel
dot_data <- combo_counts %>%
  pivot_longer(cols = c(in_meso, in_myelo, in_lymph),
               names_to = "set", values_to = "member") %>%
  mutate(set = recode(set,
                      in_meso  = "Mesenchymal",
                      in_myelo = "Myeloid",
                      in_lymph = "Lymphocyte"),
         set = factor(set, levels = c("Mesenchymal", "Myeloid", "Lymphocyte")))

p_dots <- ggplot(dot_data, aes(x = combo, y = set)) +
  geom_point(aes(color = member, size = member)) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "grey85")) +
  scale_size_manual(values  = c("TRUE" = 4,       "FALSE" = 2)) +
  geom_text(
    data = combo_counts,
    aes(x = combo, y = 3.75, label = n),
    inherit.aes = FALSE,
    size = 4, fontface = "bold"
  ) +
  coord_flip() +
  labs(title = "Gene-set overlap (UpSet matrix)",
       x = NULL, y = NULL) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.text  = element_text(size = 10))

ggsave(file.path(out_dir, "upset_overlap.png"),
       p_dots, width = 8.5, height = 4.5, dpi = 300, bg = "white")

# =========================================================
# PART 5 — Console interpretation summary
# =========================================================
cat("\n")
cat("=======================================================\n")
cat("  SIGNATURE INDEPENDENCE QC SUMMARY\n")
cat("=======================================================\n")
cat(sprintf("  Mesenchymal  : %d genes (AddModuleScore, tumor MEC/UEC clusters)\n",
            length(genes_mesenchymal)))
cat(sprintf("  Myeloid      : %d genes (AddModuleScore, classic/hypoxia/chemokine-M)\n",
            length(genes_myeloid)))
cat(sprintf("  Lymphocyte   : top %d upregulated FindMarkers genes\n",
            length(genes_lymphocyte)))
cat(sprintf("                 (from snRNA-seq Lymphocytes vs all — drives NNLS coeff)\n"))
cat("-------------------------------------------------------\n")

for (i in seq_len(nrow(overlap_summary))) {
  row    <- overlap_summary[i, ]
  shared <- strsplit(row$shared_genes, ", ")[[1]]
  shared <- shared[nzchar(shared)]
  cat(sprintf("  %s ∩ %s : %d shared genes (Jaccard = %.4f)\n",
              row$set_A, row$set_B, length(shared), row$jaccard_index))
  if (length(shared) > 0 && length(shared) <= 20) {
    cat("    Shared: ", paste(shared, collapse = ", "), "\n")
  } else if (length(shared) > 20) {
    cat("    First 20 shared: ",
        paste(head(shared, 20), collapse = ", "), "...\n")
  }
}

cat("-------------------------------------------------------\n")
interpretation <- dplyr::case_when(
  max(overlap_summary$jaccard_index[
    overlap_summary$set_A == "Mesenchymal" |
      overlap_summary$set_B == "Mesenchymal" |
      overlap_summary$set_A == "Myeloid" |
      overlap_summary$set_B == "Myeloid"],
    na.rm = TRUE) < 0.05 ~
    "GOOD: Lymphocyte NNLS genes are largely independent of both zone signatures.",
  max(overlap_summary$jaccard_index[
    overlap_summary$set_A == "Mesenchymal" |
      overlap_summary$set_B == "Mesenchymal" |
      overlap_summary$set_A == "Myeloid" |
      overlap_summary$set_B == "Myeloid"],
    na.rm = TRUE) < 0.15 ~
    "ACCEPTABLE: Minor gene overlap — spatial co-localisation is not a scoring artefact.",
  TRUE ~
    "CAUTION: Meaningful gene overlap detected — review shared genes carefully."
)
cat("  Interpretation: ", interpretation, "\n")
cat("=======================================================\n")

message("\nDone. QC outputs saved in: ", out_dir)
