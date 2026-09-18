# =============================================================================
# 02_run_semla_nnls_qc_and_save_lymph_only.R
#
# Purpose
# -------
# Rerun Semla NNLS from per-section raw objects, write coordinate fields into
# metadata, screen sections on QC, save all sections and good sections, and
# build the merged vis_good object.
#
# Only Lymphocytes are retained downstream. Deconvolution still runs against
# the FULL reference (all nnls_group cell types in meta.tsv), because NNLS
# needs competing reference programs; other outputs are simply not kept.
#
# The three QC criteria below are all LOWER bounds and therefore cannot detect
# a uniformly inflated section. See docs/pipeline.md and 02b for the manual
# inspection step this necessitates.
#
# Inputs
# ------
#   stage_input("manifest", "visium_manifest.tsv")
#   stage_input("coord_maps", "barcode_coordinate_map_all.tsv")
#   stage_dir("raw_vis")   vis_section_<sec>_raw.rds
#   ref_file("scrna_expr"), ref_file("scrna_meta")   GEO GSE125969
#
# Outputs
# -------
#   stage_dir("semla_nnls")   all_sections/{vis,nnls}, good_sections/{vis,nnls},
#                             qc/section_qc.tsv, merged/vis_good_semla_ready.rds,
#                             global_maps/, global_summary/
#
# Stochastic
# ----------
#   none
#
# Runtime
# -------
#   ~2 hours
# =============================================================================

source(here::here("config.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(semla)
  library(tidyverse)
  library(readr)
  library(tibble)
  library(patchwork)
})


# -------------------------
# Configuration
# -------------------------
raw_vis_dir   <- stage_dir("raw_vis")
phase1_map    <- stage_input("coord_maps", "barcode_coordinate_map_all.tsv")
manifest_file <- stage_input("manifest", "visium_manifest.tsv")

sn_exp_path  <- ref_file("scrna_expr")
sn_meta_path <- ref_file("scrna_meta")

# Optional: process only a subset of sections in this run
section_keep <- NULL

# QC thresholds for "good" sections (all LOWER bounds -- see header)
min_nonNA_spots  <- THRESH$min_nonNA_spots
min_prop_nonzero <- THRESH$min_prop_nonzero
min_max_ly       <- THRESH$min_max_lymphocyte

out_root            <- stage_dir("semla_nnls")
out_all_vis         <- ensure_dir(file.path(out_root, "all_sections", "vis"))
out_all_nnls        <- ensure_dir(file.path(out_root, "all_sections", "nnls"))
out_good_vis        <- ensure_dir(file.path(out_root, "good_sections", "vis"))
out_good_nnls       <- ensure_dir(file.path(out_root, "good_sections", "nnls"))
out_global_maps     <- ensure_dir(file.path(out_root, "global_maps"))
out_global_summary  <- ensure_dir(file.path(out_root, "global_summary"))
out_qc              <- ensure_dir(file.path(out_root, "qc"))
out_merged          <- ensure_dir(file.path(out_root, "merged"))

# -------------------------
# Helpers
# -------------------------
extract_barcode16 <- function(x) {
  x <- toupper(as.character(x))
  m <- regexpr("[ACGT]{16}", x)
  out <- rep(NA_character_, length(x))
  ok <- m > 0
  out[ok] <- regmatches(x, m)
  out
}

coalesce_from_candidates <- function(df, candidates, default = NULL) {
  existing <- candidates[candidates %in% colnames(df)]
  if (length(existing) == 0) {
    if (!is.null(default)) return(default)
    return(rep(NA, nrow(df)))
  }
  out <- df[[existing[1]]]
  if (length(existing) > 1) {
    for (nm in existing[-1]) {
      out <- dplyr::coalesce(out, df[[nm]])
    }
  }
  out
}

get_staffli_meta <- function(vis_sec) {
  out <- tryCatch({
    st_obj <- GetStaffli(vis_sec)
    st_meta <- st_obj@meta_data
    
    if (is.null(st_meta)) {
      tibble(cell = colnames(vis_sec))
    } else {
      st_meta <- as_tibble(st_meta)
      colnames(st_meta) <- make.unique(colnames(st_meta), sep = "_dup")
      
      if ("cell" %in% colnames(st_meta)) {
        st_meta$cell <- as.character(st_meta$cell)
      } else if ("barcode" %in% colnames(st_meta)) {
        st_meta$cell <- as.character(st_meta$barcode)
      } else if ("barcode_dup" %in% colnames(st_meta)) {
        st_meta$cell <- as.character(st_meta$barcode_dup)
      } else {
        st_meta$cell <- colnames(vis_sec)
      }
      
      st_meta %>% distinct(cell, .keep_all = TRUE)
    }
  }, error = function(e) {
    message("GetStaffli() unavailable or failed; continuing without Staffli metadata.")
    tibble(cell = colnames(vis_sec))
  })
  out
}

prepare_coord_map_section <- function(coord_map_all, sec) {
  coord_sec <- coord_map_all %>%
    filter(section_id == sec) %>%
    as_tibble()
  
  colnames(coord_sec) <- make.unique(colnames(coord_sec), sep = "_dup")
  
  # Safely standardize one barcode column to `cell`
  if ("cell" %in% colnames(coord_sec)) {
    coord_sec$cell <- as.character(coord_sec$cell)
  } else if ("barcode_raw" %in% colnames(coord_sec)) {
    coord_sec$cell <- as.character(coord_sec$barcode_raw)
  } else if ("barcode" %in% colnames(coord_sec)) {
    coord_sec$cell <- as.character(coord_sec$barcode)
  } else {
    stop("No usable raw barcode column found in phase1 map for section ", sec,
         ". Available columns: ", paste(colnames(coord_sec), collapse = ", "))
  }
  
  coord_sec %>%
    distinct(cell, .keep_all = TRUE)
}

choose_plot_coords <- function(df) {
  if (all(c("pxl_col_in_fullres", "pxl_row_in_fullres") %in% colnames(df)) &&
      sum(!is.na(df$pxl_col_in_fullres) & !is.na(df$pxl_row_in_fullres)) > 0) {
    return(list(x = "pxl_col_in_fullres", y = "pxl_row_in_fullres", reverse_y = TRUE))
  }
  
  if (all(c("x_nnls", "y_nnls") %in% colnames(df)) &&
      sum(!is.na(df$x_nnls) & !is.na(df$y_nnls)) > 0) {
    return(list(x = "x_nnls", y = "y_nnls", reverse_y = FALSE))
  }
  
  stop("No usable plot coordinates found. Available columns: ",
       paste(colnames(df), collapse = ", "))
}

# -------------------------
# Load inputs
# -------------------------
manifest <- read_tsv(manifest_file, show_col_types = FALSE)
coord_map_all <- read_tsv(phase1_map, show_col_types = FALSE)

if (!is.null(section_keep)) {
  manifest <- manifest %>% filter(section_id %in% section_keep)
  if (nrow(manifest) == 0) {
    stop("After filtering by section_keep, no sections remain in the manifest.")
  }
}

# snRNA reference
exp <- read.delim(sn_exp_path, header = TRUE, row.names = 1, check.names = FALSE)
sn  <- CreateSeuratObject(counts = as.matrix(exp), min.cells = 3, min.features = 200)

meta <- read.delim(sn_meta_path, header = TRUE, stringsAsFactors = FALSE)
rownames(meta) <- meta$Cell
meta <- meta[colnames(sn), , drop = FALSE]
stopifnot(all(rownames(meta) == colnames(sn)))
sn <- AddMetaData(sn, metadata = meta)

DefaultAssay(sn) <- "RNA"
sn <- NormalizeData(sn)
sn <- FindVariableFeatures(sn)
sn <- ScaleData(sn)
sn <- RunPCA(sn)

if (!"cell_type" %in% colnames(sn@meta.data)) {
  stop("snRNA metadata does not contain 'cell_type'.")
}
sn$nnls_group <- sn$cell_type

# -------------------------
# Process each section independently
# -------------------------
qc_list <- list()
vis_list <- list()
log_lines <- c()

for (sec in manifest$section_id) {
  message("Processing section ", sec)
  f <- file.path(raw_vis_dir, paste0("vis_section_", sec, "_raw.rds"))
  if (!file.exists(f)) stop("Missing raw section object: ", f)
  
  vis_sec <- readRDS(f)
  
  if (!"Spatial" %in% names(vis_sec@assays)) {
    stop("Section ", sec, " does not contain a 'Spatial' assay.")
  }
  
  DefaultAssay(vis_sec) <- "Spatial"
  vis_sec$section_id <- sec
  vis_sec$barcode16  <- extract_barcode16(colnames(vis_sec))
  
  vis_sec <- NormalizeData(vis_sec, normalization.method = "RC", scale.factor = 10000)
  vis_sec <- FindVariableFeatures(vis_sec, nfeatures = 3000)
  vis_sec <- ScaleData(vis_sec)
  vis_sec <- RunPCA(vis_sec)
  
  vis_sec <- RunNNLS(
    object = vis_sec,
    singlecell_object = sn,
    groups = "nnls_group"
  )
  
  DefaultAssay(vis_sec) <- "celltypeprops"
  nnls_mat <- GetAssayData(vis_sec, layer = "data")
  
  vis_sec$Lymphocytes <- if ("Lymphocytes" %in% rownames(nnls_mat)) {
    as.numeric(nnls_mat["Lymphocytes", ])
  } else {
    NA_real_
  }
  
  st_meta <- get_staffli_meta(vis_sec)
  coord_sec <- prepare_coord_map_section(coord_map_all, sec)
  
  md <- vis_sec@meta.data %>%
    rownames_to_column("cell") %>%
    left_join(st_meta, by = "cell") %>%
    left_join(coord_sec, by = "cell")
  
  md$barcode16 <- coalesce_from_candidates(
    md,
    c("barcode16.x", "barcode16.y", "barcode16"),
    default = extract_barcode16(md$cell)
  )
  
  md$section_id <- sec
  
  md$x_native <- coalesce_from_candidates(
    md,
    c("x_native", "x", "xpos", "imagecol", "col")
  )
  md$y_native <- coalesce_from_candidates(
    md,
    c("y_native", "y", "ypos", "imagerow", "row")
  )
  
  md$pxl_col_in_fullres <- coalesce_from_candidates(
    md,
    c("pxl_col_in_fullres", "pxl_col_in_fullres.x", "pxl_col_in_fullres.y", "pxl_col")
  )
  md$pxl_row_in_fullres <- coalesce_from_candidates(
    md,
    c("pxl_row_in_fullres", "pxl_row_in_fullres.x", "pxl_row_in_fullres.y", "pxl_row")
  )
  
  md$x_nnls <- dplyr::coalesce(md$pxl_col_in_fullres, md$x_native)
  md$y_nnls <- dplyr::coalesce(md$pxl_row_in_fullres, md$y_native)
  
  if (!"img_file" %in% colnames(md)) md$img_file <- NA_character_
  if (!"sample_dir" %in% colnames(md)) md$sample_dir <- NA_character_
  
  md <- md %>%
    select(-any_of(c(
      "barcode16.x", "barcode16.y",
      "section_id.x", "section_id.y"
    ))) %>%
    column_to_rownames("cell")
  
  vis_sec@meta.data <- md
  
  saveRDS(vis_sec, file.path(out_all_vis, paste0("vis_section_", sec, "_semla.rds")))
  
  df_sec <- vis_sec@meta.data %>%
    rownames_to_column("cell") %>%
    mutate(section_id = sec) %>%
    select(any_of(c(
      "cell", "section_id", "barcode16",
      "x_native", "y_native",
      "x_nnls", "y_nnls",
      "pxl_col_in_fullres", "pxl_row_in_fullres",
      "Lymphocytes",
      "img_file", "sample_dir"
    )))
  
  saveRDS(df_sec, file.path(out_all_nnls, paste0("NNLS_section_", sec, ".rds")))
  write_tsv(df_sec, file.path(out_all_nnls, paste0("NNLS_section_", sec, ".tsv")))
  
  qc_sec <- tibble(
    section_id    = sec,
    n_spots       = nrow(df_sec),
    n_nonNA       = sum(!is.na(df_sec$Lymphocytes)),
    prop_nonNA    = sum(!is.na(df_sec$Lymphocytes)) / nrow(df_sec),
    mean_ly       = mean(df_sec$Lymphocytes, na.rm = TRUE),
    sd_ly         = sd(df_sec$Lymphocytes, na.rm = TRUE),
    max_ly        = max(df_sec$Lymphocytes, na.rm = TRUE),
    prop_nonzero  = mean(df_sec$Lymphocytes > 0, na.rm = TRUE),
    n_x_native    = sum(!is.na(df_sec$x_native)),
    n_y_native    = sum(!is.na(df_sec$y_native)),
    n_x_nnls      = sum(!is.na(df_sec$x_nnls)),
    n_y_nnls      = sum(!is.na(df_sec$y_nnls)),
    n_pxl         = sum(!is.na(df_sec$pxl_col_in_fullres) & !is.na(df_sec$pxl_row_in_fullres))
  )
  
  qc_list[[sec]] <- qc_sec
  vis_list[[sec]] <- vis_sec
  log_lines <- c(
    log_lines,
    paste0(
      sec,
      " | n_spots=", ncol(vis_sec),
      " | max_ly=", signif(qc_sec$max_ly, 4),
      " | prop_nonzero=", signif(qc_sec$prop_nonzero, 4),
      " | n_pxl=", qc_sec$n_pxl
    )
  )
}

qc <- bind_rows(qc_list) %>% arrange(section_id)

qc <- qc %>%
  mutate(
    keep_good = !is.na(max_ly) &
      n_nonNA >= min_nonNA_spots &
      prop_nonzero >= min_prop_nonzero &
      max_ly >= min_max_ly
  )

write_tsv(qc, file.path(out_qc, "section_qc.tsv"))
write_lines(log_lines, file.path(out_qc, "phase2_log.txt"))

good_sections <- qc %>% filter(keep_good) %>% pull(section_id)
write_lines(good_sections, file.path(out_qc, "good_sections.txt"))

if (length(good_sections) == 0) {
  stop("No good sections passed QC. Check section_qc.tsv before continuing.")
}

for (sec in good_sections) {
  vis_sec <- vis_list[[sec]]
  df_sec  <- readRDS(file.path(out_all_nnls, paste0("NNLS_section_", sec, ".rds")))
  
  saveRDS(vis_sec, file.path(out_good_vis, paste0("vis_section_", sec, ".rds")))
  saveRDS(df_sec,  file.path(out_good_nnls, paste0("NNLS_section_", sec, ".rds")))
  write_tsv(df_sec, file.path(out_good_nnls, paste0("NNLS_section_", sec, ".tsv")))
}

vis_good <- vis_list[[good_sections[1]]]
if (length(good_sections) > 1) {
  for (sec in good_sections[-1]) {
    vis_good <- merge(vis_good, vis_list[[sec]], merge.data = TRUE)
  }
}

saveRDS(vis_good, file.path(out_merged, "vis_good_semla_ready.rds"))
write_tsv(
  vis_good@meta.data %>% rownames_to_column("cell"),
  file.path(out_merged, "vis_good_semla_ready_metadata.tsv")
)

df_all_good <- purrr::map_dfr(
  good_sections,
  ~ readRDS(file.path(out_good_nnls, paste0("NNLS_section_", .x, ".rds")))
)

df_all_good$section_id <- factor(df_all_good$section_id, levels = good_sections)

df_all_good <- df_all_good %>%
  mutate(
    Lymphocytes_cp10k_global = Lymphocytes / sum(Lymphocytes, na.rm = TRUE) * 10000,
    Lymphocytes_z_global = as.numeric(scale(Lymphocytes))
  )

write_tsv(df_all_good, file.path(out_merged, "df_all_good_sections.tsv"))

ly_lim <- quantile(df_all_good$Lymphocytes_cp10k_global, 0.99, na.rm = TRUE)

for (sec in good_sections) {
  dfp <- df_all_good %>% filter(section_id == sec)
  coord_choice <- choose_plot_coords(dfp)
  
  p <- ggplot(dfp, aes(
    x = .data[[coord_choice$x]],
    y = .data[[coord_choice$y]],
    color = Lymphocytes_cp10k_global
  )) +
    geom_point(size = 0.9) +
    scale_color_viridis_c(
      option = "plasma",
      limits = c(0, ly_lim),
      oob = scales::squish
    ) +
    coord_equal() +
    theme_void() +
    ggtitle(paste("Lymphocytes (global cp10k) –", sec))
  
  if (isTRUE(coord_choice$reverse_y)) {
    p <- p + scale_y_reverse()
  }
  
  ggsave(
    file.path(out_global_maps, paste0("Lymphocytes_global_", sec, ".png")),
    p, width = 6, height = 7, dpi = 300
  )
}

summary_bar <- df_all_good %>%
  group_by(section_id) %>%
  summarise(mean_Ly = mean(Lymphocytes_cp10k_global, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_Ly))

write_tsv(summary_bar, file.path(out_global_summary, "mean_lymphocyte_cp10k_by_section.tsv"))

p_bar <- ggplot(summary_bar, aes(x = reorder(section_id, -mean_Ly), y = mean_Ly)) +
  geom_col(fill = "purple") +
  theme_bw() +
  ylab("Mean Lymphocyte cp10k (global)") +
  xlab("Tumor Section") +
  ggtitle("Lymphocyte Infiltration Across Tumors")

p_box <- ggplot(df_all_good, aes(section_id, Lymphocytes_cp10k_global)) +
  geom_boxplot(outlier.size = 0.5, fill = "gold") +
  theme_bw() +
  ylab("Lymphocyte cp10k (global)") +
  xlab("Tumor Section") +
  ggtitle("Distribution of T-cell Infiltration Across Tumors")

p_violin <- ggplot(df_all_good, aes(section_id, Lymphocytes_cp10k_global)) +
  geom_violin(fill = "purple", alpha = 0.7) +
  theme_bw() +
  ylab("Lymphocyte cp10k (global)") +
  xlab("Tumor Section") +
  ggtitle("Lymphocyte Infiltration Across Tumors")

ggsave(
  file.path(out_global_summary, "mean_lymphocyte_cp10k_barplot.png"),
  p_bar, width = 7, height = 4.5, dpi = 300
)
ggsave(
  file.path(out_global_summary, "lymphocyte_cp10k_boxplot.png"),
  p_box, width = 7, height = 4.5, dpi = 300
)
ggsave(
  file.path(out_global_summary, "lymphocyte_cp10k_violin.png"),
  p_violin, width = 7, height = 4.5, dpi = 300
)

message("Phase 2 done. vis_good_semla_ready.rds and clean good_sections outputs are ready.")
