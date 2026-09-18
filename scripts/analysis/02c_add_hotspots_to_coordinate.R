# =============================================================================
# 02c_add_hotspots_to_coordinate.R
#
# Purpose
# -------
# Rebuild lymphocyte hotspots from the curated NNLS tables, pulling the tissue
# mask and geometry directly from the stage 00 spot tables.
#
# HOTSPOT DEFINITION
# local_lymph_score = mean Lymphocytes score over the k nearest neighbours.
# A spot is a hotspot when Lymphocytes >= the section's own quantile AND
# local_lymph_score >= the section's own quantile. Both thresholds are
# WITHIN-SECTION, so hotspots exist in every section by construction and are
# capped near 10% of spots. The hotspot percentage is therefore a property of
# the definition, not a measurement, and is not comparable across sections.
#
# These fields are guaranteed to come from the corrected stage 00 backbone:
#   in_tissue, img_file / json_file,
#   pxl_col_in_fullres / pxl_row_in_fullres,
#   pxl_col_in_hires / pxl_row_in_hires
#
# Inputs
# ------
#   stage_dir("raw_spot_tables")   spot_table_<sec>_raw.tsv
#   stage_dir("semla_nnls")        good_sections/{vis,nnls}
#
# Outputs
# -------
#   stage_dir("hotspots")   vis/vis_section_<sec>.rds  (adds is_hotspot_semla,
#                           local_lymph_score to meta.data),
#                           nnls/NNLS_section_<sec>.{rds,tsv},
#                           qc/rebuilt_hotspot_summary_tissue.tsv, plots/
#
# Stochastic
# ----------
#   none
#
# Runtime
# -------
#   ~15 minutes
# =============================================================================

source(here::here("config.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(tibble)
  library(dbscan)
  library(ggplot2)
})


# -------------------------
# Configuration
# -------------------------
phase0_spot_dir <- stage_dir("raw_spot_tables")
phase2_root     <- stage_dir("semla_nnls")

vis_dir  <- file.path(phase2_root, "good_sections", "vis")
nnls_dir <- file.path(phase2_root, "good_sections", "nnls")

out_root <- stage_dir("hotspots")
out_vis  <- ensure_dir(file.path(out_root, "vis"))
out_nnls <- ensure_dir(file.path(out_root, "nnls"))
out_qc   <- ensure_dir(file.path(out_root, "qc"))
out_plot <- ensure_dir(file.path(out_root, "plots"))

section_keep <- GOOD_SECTIONS

k_neighbors    <- THRESH$hotspot_knn
expr_quantile  <- THRESH$hotspot_quantile
local_quantile <- THRESH$hotspot_quantile

extract_barcode16 <- function(x) {
  x <- toupper(as.character(x))
  m <- regexpr("[ACGT]{16}", x)
  out <- rep(NA_character_, length(x))
  ok <- m > 0
  out[ok] <- regmatches(x, m)
  out
}

local_knn_mean <- function(xy, values, k = 6) {
  if (nrow(xy) < 2) return(rep(NA_real_, nrow(xy)))
  k_use <- min(k, nrow(xy) - 1L)
  if (k_use < 1) return(rep(NA_real_, nrow(xy)))

  kn <- dbscan::kNN(xy, k = k_use)
  idx <- kn$id

  out <- numeric(nrow(idx))
  for (i in seq_len(nrow(idx))) {
    nbrs <- idx[i, ]
    out[i] <- mean(values[c(i, nbrs)], na.rm = TRUE)
  }
  out
}

summary_list <- list()

for (sec in section_keep) {
  message("Hotspot reconstruction from corrected Phase 0 backbone for section ", sec)

  nnls_file <- file.path(nnls_dir, paste0("NNLS_section_", sec, ".rds"))
  vis_file  <- file.path(vis_dir,  paste0("vis_section_", sec, ".rds"))
  spot_file <- file.path(phase0_spot_dir, paste0("spot_table_", sec, "_raw.tsv"))

  if (!file.exists(nnls_file)) stop("Missing NNLS file: ", nnls_file)
  if (!file.exists(vis_file))  stop("Missing vis file: ", vis_file)
  if (!file.exists(spot_file)) stop("Missing Phase 0 spot table: ", spot_file)

  df <- readRDS(nnls_file)
  if (!"cell" %in% colnames(df)) stop("NNLS table missing 'cell' for section ", sec)
  if (!"Lymphocytes" %in% colnames(df)) stop("NNLS table missing 'Lymphocytes' for section ", sec)

  df$barcode16_join <- if ("barcode16" %in% colnames(df)) {
    dplyr::coalesce(as.character(df$barcode16), extract_barcode16(df$cell))
  } else {
    extract_barcode16(df$cell)
  }

  spot_tbl <- readr::read_tsv(spot_file, show_col_types = FALSE) %>%
    mutate(
      barcode16_join = dplyr::coalesce(as.character(barcode16), extract_barcode16(cell)),
      in_tissue = as.integer(in_tissue),
      pxl_col_in_fullres = as.numeric(pxl_col_in_fullres),
      pxl_row_in_fullres = as.numeric(pxl_row_in_fullres),
      pxl_col_in_hires = as.numeric(pxl_col_in_hires),
      pxl_row_in_hires = as.numeric(pxl_row_in_hires)
    ) %>%
    distinct(barcode16_join, .keep_all = TRUE) %>%
    select(
      barcode16_join,
      img_file, json_file,
      in_tissue, array_row, array_col,
      pxl_col_in_fullres, pxl_row_in_fullres,
      pxl_col_in_hires, pxl_row_in_hires,
      tissue_hires_scalef, img_width_hires, img_height_hires
    )

  df2 <- df %>%
    select(-any_of(c(
      "in_tissue", "array_row", "array_col",
      "img_file", "json_file",
      "pxl_col_in_fullres.x", "pxl_row_in_fullres.x",
      "pxl_col_in_fullres.y", "pxl_row_in_fullres.y",
      "pxl_col_in_fullres", "pxl_row_in_fullres",
      "pxl_col_in_hires", "pxl_row_in_hires",
      "is_hotspot_semla", "local_lymph_score"
    ))) %>%
    left_join(spot_tbl, by = "barcode16_join") %>%
    mutate(
      section_id = sec,
      Lymphocytes = as.numeric(Lymphocytes),
      in_tissue = as.integer(coalesce(in_tissue, 0L))
    )

  valid <- with(
    df2,
    in_tissue == 1L &
      !is.na(pxl_col_in_hires) &
      !is.na(pxl_row_in_hires) &
      !is.na(Lymphocytes)
  )

  df2$local_lymph_score <- NA_real_
  df2$is_hotspot_semla  <- FALSE

  n_valid <- sum(valid, na.rm = TRUE)

  if (n_valid >= 2) {
    xy <- as.matrix(df2[valid, c("pxl_col_in_hires", "pxl_row_in_hires")])
    df2$local_lymph_score[valid] <- local_knn_mean(xy, df2$Lymphocytes[valid], k = k_neighbors)

    expr_thr  <- stats::quantile(df2$Lymphocytes[valid], expr_quantile, na.rm = TRUE)
    local_thr <- stats::quantile(df2$local_lymph_score[valid], local_quantile, na.rm = TRUE)

    df2$is_hotspot_semla[valid] <-
      (df2$Lymphocytes[valid] >= expr_thr) &
      (df2$local_lymph_score[valid] >= local_thr)
  } else {
    expr_thr  <- NA_real_
    local_thr <- NA_real_
  }

  saveRDS(df2, file.path(out_nnls, paste0("NNLS_section_", sec, ".rds")))
  write_tsv(df2, file.path(out_nnls, paste0("NNLS_section_", sec, ".tsv")))

  obj <- readRDS(vis_file)
  md <- obj@meta.data %>%
    rownames_to_column("cell") %>%
    select(-any_of(c(
      "barcode16_join",
      "in_tissue", "array_row", "array_col",
      "img_file", "json_file",
      "pxl_col_in_fullres", "pxl_row_in_fullres",
      "pxl_col_in_hires", "pxl_row_in_hires",
      "is_hotspot_semla", "local_lymph_score"
    ))) %>%
    left_join(
      df2 %>%
        select(
          cell, barcode16_join,
          in_tissue, array_row, array_col,
          img_file, json_file,
          pxl_col_in_fullres, pxl_row_in_fullres,
          pxl_col_in_hires, pxl_row_in_hires,
          local_lymph_score, is_hotspot_semla
        ),
      by = "cell"
    ) %>%
    column_to_rownames("cell")

  obj@meta.data <- md
  saveRDS(obj, file.path(out_vis, paste0("vis_section_", sec, ".rds")))

  p_mask <- ggplot(df2, aes(pxl_col_in_hires, pxl_row_in_hires, color = factor(in_tissue))) +
    geom_point(size = 0.8) +
    scale_color_manual(values = c("0" = "grey80", "1" = "dodgerblue3")) +
    scale_y_reverse() + coord_equal() + theme_void(base_size = 11) +
    ggtitle(paste0("Tissue mask — ", sec))

  p_hot <- ggplot(df2, aes(pxl_col_in_hires, pxl_row_in_hires)) +
    geom_point(data = df2 %>% filter(in_tissue == 0), color = "grey85", size = 0.8) +
    geom_point(data = df2 %>% filter(in_tissue == 1), color = "grey35", size = 0.35) +
    geom_point(data = df2 %>% filter(is_hotspot_semla), color = "red3", size = 1.1) +
    scale_y_reverse() + coord_equal() + theme_void(base_size = 11) +
    ggtitle(paste0("Hotspots (tissue-only) — ", sec))

  p_local <- ggplot(df2 %>% filter(in_tissue == 1), aes(pxl_col_in_hires, pxl_row_in_hires, color = local_lymph_score)) +
    geom_point(size = 0.8) +
    scale_color_viridis_c(option = "magma", oob = scales::squish) +
    scale_y_reverse() + coord_equal() + theme_void(base_size = 11) +
    ggtitle(paste0("Local lymph score (tissue-only) — ", sec))

  ggsave(file.path(out_plot, paste0("QC_TissueMask_", sec, ".png")), p_mask, width = 6, height = 7, dpi = 300, bg = "white")
  ggsave(file.path(out_plot, paste0("QC_Hotspots_", sec, ".png")), p_hot, width = 6, height = 7, dpi = 300, bg = "white")
  ggsave(file.path(out_plot, paste0("QC_LocalLymphScore_", sec, ".png")), p_local, width = 6, height = 7, dpi = 300, bg = "white")

  summary_list[[sec]] <- tibble(
    section_id = sec,
    n_spots = nrow(df2),
    n_tissue = sum(df2$in_tissue == 1, na.rm = TRUE),
    n_hotspots = sum(df2$is_hotspot_semla, na.rm = TRUE),
    pct_hotspots_of_tissue = 100 * mean(df2$is_hotspot_semla[df2$in_tissue == 1], na.rm = TRUE),
    n_valid_for_hotspots = n_valid,
    expr_threshold = expr_thr,
    local_threshold = local_thr
  )
}

summary_df <- bind_rows(summary_list) %>% arrange(section_id)
write_tsv(summary_df, file.path(out_qc, "rebuilt_hotspot_summary_tissue.tsv"))
write_lines(section_keep, file.path(out_qc, "section_keep_used.txt"))

message("Phase 2b v5 done. Outputs saved under: ", out_root)
