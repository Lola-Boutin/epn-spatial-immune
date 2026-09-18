# =============================================================================
# 01_build_raw_barcode_coordinate_maps.R
#
# Purpose
# -------
# Build clean barcode / coordinate / H&E maps from the per-section raw
# Visium objects created in stage 00.
#
# Inputs
# ------
#   stage_input("manifest", "visium_manifest.tsv")
#   stage_dir("raw_vis")   vis_section_<sec>_raw.rds
#
# Outputs
# -------
#   stage_dir("coord_maps")   barcode_coordinate_map_<sec>.tsv,
#                             barcode_coordinate_map_all.tsv,
#                             barcode_coordinate_map_qc.tsv
#
# Stochastic
# ----------
#   none
#
# Runtime
# -------
#   ~5 minutes
# =============================================================================

source(here::here("config.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(semla)
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
})


# -------------------------
# Configuration
# -------------------------
manifest_file <- stage_input("manifest", "visium_manifest.tsv")
raw_vis_dir   <- stage_dir("raw_vis")

out_dir <- stage_dir("coord_maps")

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

manifest <- read_tsv(manifest_file, show_col_types = FALSE)

section_logs <- list()
map_list <- list()
for (sec in manifest$section_id) {
  f <- file.path(raw_vis_dir, paste0("vis_section_", sec, "_raw.rds"))
  if (!file.exists(f)) stop("Missing raw object for section ", sec, ": ", f)

  obj <- readRDS(f)
  md <- obj@meta.data %>%
    rownames_to_column("cell") %>%
    mutate(
      section_id = sec,
      barcode_raw = cell,
      barcode16   = extract_barcode16(cell)
    )

  st_obj <- tryCatch(GetStaffli(obj), error = function(e) NULL)
  if (is.null(st_obj)) {
    stop("No Staffli object found in raw section ", sec, ". Raw image-linked coordinates are required.")
  }

  st_meta <- st_obj@meta_data %>%
    as_tibble() %>%
    rename(barcode_raw = barcode)

  join_map <- md %>%
    left_join(st_meta, by = "barcode_raw") %>%
    mutate(
      sample_dir  = manifest$sample_dir[match(sec, manifest$section_id)],
      h5_file     = manifest$h5_file[match(sec, manifest$section_id)],
      img_file    = manifest$img_file[match(sec, manifest$section_id)],
      spot_file   = manifest$spot_file[match(sec, manifest$section_id)],
      json_file   = manifest$json_file[match(sec, manifest$section_id)]
    )

  # Canonical fields for downstream joins
  join_map <- join_map %>%
    mutate(
      x_native = if ("x" %in% colnames(join_map)) x else NA_real_,
      y_native = if ("y" %in% colnames(join_map)) y else NA_real_,
      pxl_col_in_fullres = if ("pxl_col_in_fullres" %in% colnames(join_map)) pxl_col_in_fullres else NA_real_,
      pxl_row_in_fullres = if ("pxl_row_in_fullres" %in% colnames(join_map)) pxl_row_in_fullres else NA_real_
    )

  exact_rate <- mean(!is.na(join_map$x_native) | !is.na(join_map$pxl_col_in_fullres))
  dup_barcodes <- sum(duplicated(join_map$barcode_raw))
  dup_bc16     <- sum(duplicated(join_map$barcode16[!is.na(join_map$barcode16)]))

  readr::write_tsv(join_map, file.path(out_dir, paste0("barcode_coordinate_map_", sec, ".tsv")))
  map_list[[sec]] <- join_map

  section_logs[[sec]] <- tibble(
    section_id = sec,
    n_spots = nrow(join_map),
    exact_join_rate = exact_rate,
    n_missing_native = sum(is.na(join_map$x_native) & is.na(join_map$pxl_col_in_fullres)),
    duplicated_barcode_raw = dup_barcodes,
    duplicated_barcode16 = dup_bc16
  )
}

map_all <- bind_rows(map_list)
log_df  <- bind_rows(section_logs) %>% arrange(section_id)

write_tsv(map_all, file.path(out_dir, "barcode_coordinate_map_all.tsv"))
write_tsv(log_df, file.path(out_dir, "barcode_coordinate_map_qc.tsv"))

if (any(log_df$exact_join_rate < 0.95)) {
  warning("Some sections have coordinate join rate < 95%. Check barcode_coordinate_map_qc.tsv")
}
if (any(log_df$duplicated_barcode_raw > 0)) {
  warning("Some sections have duplicated raw barcodes. Check barcode_coordinate_map_qc.tsv")
}

message("Phase 1 done. Master barcode/coordinate maps saved under: ", out_dir)
