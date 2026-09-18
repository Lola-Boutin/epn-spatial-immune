# =============================================================================
# 00_build_manifest_and_raw_sections.R
#
# Purpose
# -------
# Discover raw Visium inputs, build the section manifest, save per-section
# raw Seurat objects, and export a barcode / coordinate / tissue backbone
# with H&E QC overlays.
#
# Reads tissue_positions_list.csv and scalefactors_json.json explicitly,
# converts fullres spot coordinates to hires PNG coordinates, and stores
# both. QC overlays use base-R rasterImage to avoid y-axis flip issues.
#
# Inputs
# ------
#   RAW_ROOT/<section>/spatial/   (from GEO GSE195661)
#
# Outputs
# -------
#   stage_dir("manifest")          visium_manifest.tsv, visium_manifest_batched.tsv,
#                                  visium_manifest_incomplete.tsv, phase0_qc_summary.tsv
#   stage_dir("raw_vis")           vis_section_<sec>_raw.rds
#   stage_dir("raw_spot_tables")   spot_table_<sec>_raw.tsv
#   stage_dir("raw_qc_png")        H&E overlay QC images
#
# Stochastic
# ----------
#   none
#
# Runtime
# -------
#   ~30 minutes
# =============================================================================

source(here::here("config.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(semla)
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
  library(stringr)
  library(jsonlite)
  library(png)
})


# -------------------------
# Configuration
# -------------------------
base_spatial <- RAW_ROOT

# Set to NULL to scan everything automatically
section_keep <- NULL
batch_size   <- 14

# -------------------------
# Output dirs
# -------------------------
manifest_dir <- stage_dir("manifest")
raw_vis_dir  <- stage_dir("raw_vis")
raw_tbl_dir  <- stage_dir("raw_spot_tables")
raw_qc_dir   <- stage_dir("raw_qc_png")
log_dir      <- manifest_dir

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

read_spot_positions <- function(path) {
  x <- tryCatch(
    readr::read_csv(path, show_col_types = FALSE, col_names = FALSE),
    error = function(e) NULL
  )

  if (!is.null(x) && ncol(x) >= 6) {
    colnames(x)[1:6] <- c(
      "barcode", "in_tissue", "array_row", "array_col",
      "pxl_row_in_fullres", "pxl_col_in_fullres"
    )
    x <- x[, 1:6]
  } else {
    x <- readr::read_csv(path, show_col_types = FALSE)
    nms <- colnames(x)
    if (all(c("barcode", "in_tissue", "array_row", "array_col",
              "pxl_row_in_fullres", "pxl_col_in_fullres") %in% nms)) {
      x <- x %>%
        dplyr::select(
          barcode, in_tissue, array_row, array_col,
          pxl_row_in_fullres, pxl_col_in_fullres
        )
    } else {
      stop("Unsupported tissue_positions format: ", path)
    }
  }

  x$barcode <- as.character(x$barcode)
  x$barcode16 <- extract_barcode16(x$barcode)
  x$in_tissue <- as.integer(x$in_tissue)
  x$array_row <- as.integer(x$array_row)
  x$array_col <- as.integer(x$array_col)
  x$pxl_row_in_fullres <- as.numeric(x$pxl_row_in_fullres)
  x$pxl_col_in_fullres <- as.numeric(x$pxl_col_in_fullres)
  x %>% distinct(barcode16, .keep_all = TRUE)
}

read_scalefactors <- function(path) {
  sf <- jsonlite::fromJSON(path)
  tibble(
    tissue_hires_scalef = as.numeric(sf$tissue_hires_scalef),
    tissue_lowres_scalef = as.numeric(sf$tissue_lowres_scalef),
    fiducial_diameter_fullres = as.numeric(sf$fiducial_diameter_fullres),
    spot_diameter_fullres = as.numeric(sf$spot_diameter_fullres)
  )
}

make_raw_he_overlay <- function(df, sec, outfile) {
  img_path <- unique(df$img_file)
  if (length(img_path) != 1 || is.na(img_path) || !file.exists(img_path)) return(FALSE)

  req <- c("pxl_col_in_hires", "pxl_row_in_hires", "in_tissue")
  if (!all(req %in% colnames(df))) return(FALSE)

  img <- png::readPNG(img_path)
  img_h <- nrow(img)
  img_w <- ncol(img)

  # convert top-left image coordinates to base-R bottom-left plotting coordinates
  df_plot <- df %>%
    mutate(
      x_plot = pxl_col_in_hires,
      y_plot = img_h - pxl_row_in_hires
    )

  png(filename = outfile, width = 2100, height = 2100, res = 300)
  op <- par(no.readonly = TRUE)
  on.exit({
    par(op)
    dev.off()
  }, add = TRUE)

  par(mar = c(1, 1, 3, 1))
  plot(
    NA,
    xlim = c(0, img_w),
    ylim = c(0, img_h),
    xaxs = "i", yaxs = "i",
    asp = 1,
    xaxt = "n", yaxt = "n",
    xlab = "", ylab = "",
    bty = "n",
    main = paste0("Raw tissue mask on H&E — ", sec)
  )
  rasterImage(as.raster(img), 0, 0, img_w, img_h)

  # off-tissue rings
  d0 <- df_plot %>% filter(in_tissue == 0)
  if (nrow(d0) > 0) {
    points(d0$x_plot, d0$y_plot, pch = 1, cex = 0.85, col = "grey60", lwd = 1)
  }

  # tissue spots
  d1 <- df_plot %>% filter(in_tissue == 1)
  if (nrow(d1) > 0) {
    points(d1$x_plot, d1$y_plot, pch = 16, cex = 0.18, col = "deepskyblue3")
  }

  TRUE
}

# -------------------------
# Discover raw sample dirs
# -------------------------
h5_files <- list.files(
  base_spatial, recursive = TRUE, full.names = TRUE,
  pattern = "filtered_feature_bc_matrix\\.h5$"
)

if (length(h5_files) == 0) {
  stop("No filtered_feature_bc_matrix.h5 files found under: ", base_spatial)
}

sample_dirs <- dirname(h5_files)

manifest <- tibble(
  section_id = basename(sample_dirs),
  sample_dir = sample_dirs,
  h5_file    = h5_files,
  img_file   = file.path(sample_dirs, "spatial", "tissue_hires_image.png"),
  spot_file  = file.path(sample_dirs, "spatial", "tissue_positions_list.csv"),
  json_file  = file.path(sample_dirs, "spatial", "scalefactors_json.json")
) %>%
  distinct(section_id, .keep_all = TRUE) %>%
  mutate(
    has_h5   = file.exists(h5_file),
    has_img  = file.exists(img_file),
    has_spot = file.exists(spot_file),
    has_json = file.exists(json_file),
    complete = has_h5 & has_img & has_spot & has_json
  ) %>%
  arrange(section_id)

if (!is.null(section_keep)) {
  manifest <- manifest %>% filter(section_id %in% section_keep)
}

write_tsv(manifest, file.path(manifest_dir, "visium_manifest.tsv"))

bad <- manifest %>% filter(!complete)
if (nrow(bad) > 0) {
  write_tsv(bad, file.path(log_dir, "visium_manifest_incomplete.tsv"))
  stop("Some sections are missing raw Visium files. Check visium_manifest_incomplete.tsv")
}

manifest <- manifest %>% mutate(batch_id = ceiling(row_number() / batch_size))
write_tsv(manifest, file.path(manifest_dir, "visium_manifest_batched.tsv"))

# -------------------------
# Load in batches and save per-section raw outputs
# -------------------------
save_log <- c()
qc_rows <- list()

for (bid in sort(unique(manifest$batch_id))) {
  man_b <- manifest %>% filter(batch_id == bid)
  message("Loading raw batch ", bid, ": ", paste(man_b$section_id, collapse = ", "))

  infoTable <- man_b %>%
    transmute(
      samples    = h5_file,
      imgs       = img_file,
      spotfiles  = spot_file,
      json       = json_file,
      section_id = section_id
    )

  vis_b <- ReadVisiumData(infoTable)
  vis_b <- LoadImages(vis_b)
  DefaultAssay(vis_b) <- "Spatial"

  for (sec in man_b$section_id) {
    cells_sec <- colnames(vis_b)[vis_b$section_id == sec]
    vis_sec <- subset(vis_b, cells = cells_sec)
    vis_sec$section_id <- sec
    vis_sec$barcode16  <- extract_barcode16(colnames(vis_sec))

    saveRDS(vis_sec, file.path(raw_vis_dir, paste0("vis_section_", sec, "_raw.rds")))

    man_row <- man_b %>% filter(section_id == sec)

    st_obj <- tryCatch(GetStaffli(vis_sec), error = function(e) NULL)
    if (!is.null(st_obj) && !is.null(st_obj@meta_data)) {
      st_meta <- st_obj@meta_data %>% as_tibble()
      colnames(st_meta) <- make.unique(colnames(st_meta), sep = "_dup")
      if ("barcode" %in% colnames(st_meta)) {
        st_meta$cell <- as.character(st_meta$barcode)
      } else if (!"cell" %in% colnames(st_meta)) {
        st_meta$cell <- colnames(vis_sec)
      }
      st_meta$barcode16 <- extract_barcode16(st_meta$cell)
    } else {
      st_meta <- tibble(
        cell = colnames(vis_sec),
        barcode16 = extract_barcode16(colnames(vis_sec))
      )
    }

    spot_tbl <- read_spot_positions(man_row$spot_file[[1]])
    sf_tbl   <- read_scalefactors(man_row$json_file[[1]])

    img_arr <- png::readPNG(man_row$img_file[[1]])
    img_h <- nrow(img_arr)
    img_w <- ncol(img_arr)

    raw_tbl <- vis_sec@meta.data %>%
      rownames_to_column("cell") %>%
      mutate(
        section_id = sec,
        barcode16 = extract_barcode16(cell),
        sample_dir = man_row$sample_dir[[1]],
        h5_file = man_row$h5_file[[1]],
        img_file = man_row$img_file[[1]],
        spot_file = man_row$spot_file[[1]],
        json_file = man_row$json_file[[1]]
      ) %>%
      left_join(
        st_meta %>% select(any_of(c(
          "cell", "barcode16", "x", "y", "sampleID",
          "pxl_col_in_fullres", "pxl_row_in_fullres"
        ))),
        by = c("cell", "barcode16")
      ) %>%
      left_join(
        spot_tbl %>%
          rename(
            barcode_raw_from_spot = barcode,
            in_tissue = in_tissue,
            array_row = array_row,
            array_col = array_col,
            pxl_row_in_fullres_spot = pxl_row_in_fullres,
            pxl_col_in_fullres_spot = pxl_col_in_fullres
          ) %>%
          select(
            barcode16,
            barcode_raw_from_spot,
            in_tissue,
            array_row,
            array_col,
            pxl_row_in_fullres_spot,
            pxl_col_in_fullres_spot
          ),
        by = "barcode16"
      ) %>%
      mutate(
        tissue_hires_scalef = sf_tbl$tissue_hires_scalef[[1]],
        tissue_lowres_scalef = sf_tbl$tissue_lowres_scalef[[1]],
        fiducial_diameter_fullres = sf_tbl$fiducial_diameter_fullres[[1]],
        spot_diameter_fullres = sf_tbl$spot_diameter_fullres[[1]],
        # canonical fullres coords from spot file
        pxl_col_in_fullres = as.numeric(pxl_col_in_fullres_spot),
        pxl_row_in_fullres = as.numeric(pxl_row_in_fullres_spot),
        # hires PNG coords
        pxl_col_in_hires = pxl_col_in_fullres * tissue_hires_scalef,
        pxl_row_in_hires = pxl_row_in_fullres * tissue_hires_scalef,
        # image canvas info
        img_width_hires = img_w,
        img_height_hires = img_h,
        # optional Staffli/native coords
        x_native = if ("x" %in% colnames(.)) x else NA_real_,
        y_native = if ("y" %in% colnames(.)) y else NA_real_,
        # convenience bottom-left plotting coords
        x_plot_hires = pxl_col_in_hires,
        y_plot_hires = img_h - pxl_row_in_hires
      )

    write_tsv(raw_tbl, file.path(raw_tbl_dir, paste0("spot_table_", sec, "_raw.tsv")))

    qc_ok <- make_raw_he_overlay(
      raw_tbl,
      sec = sec,
      outfile = file.path(raw_qc_dir, paste0("QC_raw_spots_on_HE_", sec, ".png"))
    )

    n_barcode16_matched <- sum(!is.na(raw_tbl$in_tissue), na.rm = TRUE)
    n_tissue <- sum(raw_tbl$in_tissue == 1, na.rm = TRUE)
    n_on_canvas <- sum(
      !is.na(raw_tbl$pxl_col_in_hires) &
        !is.na(raw_tbl$pxl_row_in_hires) &
        raw_tbl$pxl_col_in_hires >= 0 &
        raw_tbl$pxl_col_in_hires <= img_w &
        raw_tbl$pxl_row_in_hires >= 0 &
        raw_tbl$pxl_row_in_hires <= img_h,
      na.rm = TRUE
    )

    qc_rows[[sec]] <- tibble(
      section_id = sec,
      n_spots = ncol(vis_sec),
      n_barcode16_matched = n_barcode16_matched,
      barcode16_match_rate = n_barcode16_matched / ncol(vis_sec),
      n_tissue = n_tissue,
      hires_scale = sf_tbl$tissue_hires_scalef[[1]],
      img_width_hires = img_w,
      img_height_hires = img_h,
      n_on_canvas = n_on_canvas,
      on_canvas_rate = n_on_canvas / ncol(vis_sec),
      qc_png_written = qc_ok
    )

    save_log <- c(
      save_log,
      paste0(
        sec,
        "\tn_spots=", ncol(vis_sec),
        "\tn_barcode16_matched=", n_barcode16_matched,
        "\tn_tissue=", n_tissue,
        "\thires_scale=", signif(sf_tbl$tissue_hires_scalef[[1]], 5),
        "\tn_on_canvas=", n_on_canvas,
        "\tqc_png=", qc_ok
      )
    )
  }

  rm(vis_b)
  invisible(gc())
}

qc_df <- bind_rows(qc_rows) %>% arrange(section_id)
write_tsv(qc_df, file.path(log_dir, "phase0_qc_summary.tsv"))
write_lines(save_log, file.path(log_dir, "phase0_save_log.txt"))

message("Phase 0 v2 fixed done. Raw per-section objects and tissue-aware spot tables saved under: ", rerun_root)
