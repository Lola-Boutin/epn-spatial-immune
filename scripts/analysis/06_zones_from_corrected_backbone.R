# =============================================================================
# 06_zones_from_corrected_backbone.R
#
# Purpose
# -------
# Assign zone identities (Epithelial / Mesenchymal / Vascular / Myeloid /
# Uncertain) per spot and relate them to hotspot context.
#
# Geometry backbone comes from stage 00; hotspot calls from 02c. Zone gene
# programs are built from the source publication's supplementary workbook.
# Normal-brain / granule analysis has been removed; there is no visHE
# dependency.
#
# The gene-set construction here MUST stay in step with 06b, which re-derives
# the same sets to test signature independence.
#
# RUN THIS SCRIPT TWICE -- see ALL_SECTIONS_RUN in the config block. The
# six-section run feeds 07 and 08; the fourteen-section run feeds 09 and 10.
#
# Inputs
# ------
#   stage_dir("raw_vis")           vis_section_<sec>_raw.rds
#   stage_dir("raw_spot_tables")   spot tables
#   stage_dir("hotspots")          nnls/NNLS_section_<sec>.rds
#   ref_file("zone_genes_xlsx")    Donson et al. supplementary data 1
#
# Outputs
# -------
#   stage_dir("zones")   tables/zones_section_<sec>.tsv,
#                        tables/zones_all_sections.tsv,
#                        tables/zone_hotspot_stats.tsv,
#                        tables/barrier_neighbor_stats.tsv,
#                        plots/, rds/zones_section_<sec>.rds
#
# Stochastic
# ----------
#   none
#
# Runtime
# -------
#   ~40 minutes
# =============================================================================

source(here::here("config.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(purrr)
  library(readr)
  library(readxl)
  library(tibble)
  library(ggplot2)
  library(dbscan)
  library(png)
})


# -------------------------
# Configuration
# -------------------------
raw_vis_dir  <- stage_dir("raw_vis")
raw_spot_dir <- stage_dir("raw_spot_tables")
hotspot_dir  <- file.path(stage_dir("hotspots"), "nnls")

xlsx <- ref_file("zone_genes_xlsx")

# -----------------------------------------------------------------------------
# RUN MODE
#
# This stage is run TWICE, producing two separate zone tables:
#
#   FALSE  ->  GOOD_SECTIONS (6)  ->  stage "zones"      ->  consumed by 07, 08
#   TRUE   ->  ALL_SECTIONS (14)  ->  stage "zones_all"  ->  consumed by 09, 10
#
# Run it once each way. Hotspot calls exist only for GOOD_SECTIONS, so in the
# all-sections run the eight remaining sections carry empty hotspot context
# (is_hotspot_semla = FALSE) -- which is what 09 and 10 expect, since they use
# zone_call and not hotspot status.
# -----------------------------------------------------------------------------
ALL_SECTIONS_RUN <- FALSE

section_keep <- if (ALL_SECTIONS_RUN) ALL_SECTIONS else GOOD_SECTIONS

out_root   <- stage_dir(if (ALL_SECTIONS_RUN) "zones_all" else "zones")
out_tables <- ensure_dir(file.path(out_root, "tables"))
out_plots  <- ensure_dir(file.path(out_root, "plots"))
out_rds    <- ensure_dir(file.path(out_root, "rds"))

message("Zone run mode: ", if (ALL_SECTIONS_RUN) "ALL sections" else "GOOD sections only",
        " (", length(section_keep), " sections) -> ", out_root)

# Zone gene-set construction -- 06b re-derives these and MUST match exactly.
TOP_N_PER_CLUSTER <- 100
MIN_DELTA_PCT     <- 0.15
USE_METRIC        <- "score"
MIN_PCT1   <- 0.20
MIN_LOG2FC <- 0.40
MAX_PCT2   <- 0.30
min_max_z <- 0.5
k_neighbors <- THRESH$hotspot_knn

zone_cols <- c(
  "Epithelial"  = "#4DAF4A",
  "Mesenchymal" = "#E41A1C",
  "Vascular"    = "#377EB8",
  "Myeloid"     = "#984EA3",
  "Uncertain"   = "grey70"
)

extract_barcode16 <- function(x) {
  x <- toupper(as.character(x))
  m <- regexpr("[ACGT]{16}", x)
  out <- rep(NA_character_, length(x))
  ok <- m > 0
  out[ok] <- regmatches(x, m)
  out
}

read_expr_matrix <- function(obj) {
  DefaultAssay(obj) <- "Spatial"
  mat <- tryCatch(GetAssayData(obj, assay = "Spatial", layer = "data"), error = function(e) NULL)
  if (is.null(mat) || nrow(mat) == 0 || ncol(mat) == 0) {
    message("  no usable 'data' layer found; normalising from counts")
    obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)
    mat <- tryCatch(GetAssayData(obj, assay = "Spatial", layer = "data"), error = function(e) NULL)
  }
  if (is.null(mat) || nrow(mat) == 0 || ncol(mat) == 0) {
    stop("Could not obtain normalized expression matrix from Spatial assay")
  }
  list(obj = obj, mat = mat)
}

safe_module_score <- function(obj, features, name) {
  features <- features[lengths(features) > 0]
  if (length(features) == 0) stop("No non-empty feature sets supplied to AddModuleScore")
  AddModuleScore(obj, features = features, name = name)
}

z_by_section <- function(v, sec) {
  out <- rep(NA_real_, length(v))
  for (s in unique(sec)) {
    idx <- which(sec == s)
    vv <- v[idx]
    vv_ok <- vv[is.finite(vv)]
    if (length(vv_ok) < 2 || length(unique(vv_ok)) < 2) {
      out[idx] <- NA_real_
    } else {
      out[idx] <- as.numeric(scale(vv))
    }
  }
  out
}

read_one_sheet <- function(sheet, zone_name, xlsx) {
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

build_zone_gene_sets <- function(xlsx) {
  zone_sheets <- list(
    Epithelial  = c("TEC-A","TEC-B","TEC-C","TEC-D","UEC-A","CEC"),
    Mesenchymal = c("MEC-A","MEC-B","MEC-C","MEC-D","UEC-B"),
    Vascular    = c("VE"),
    Myeloid     = c("classic-M","hypoxia-M","chemokine-M")
  )
  endo_markers <- toupper(c(
    "PECAM1","VWF","KDR","FLT1","TEK","ENG","EMCN","RAMP2","PLVAP",
    "CD34","ESAM","ROBO4","PGF","KLF2","KLF4","CLEC14A","FABP4",
    "MCAM","SOX17","EDNRB","ADGRL4","PROM1"
  ))
  available <- excel_sheets(xlsx)
  wanted <- unlist(zone_sheets, use.names = FALSE)
  missing <- setdiff(wanted, available)
  if (length(missing) > 0) stop("Missing Excel sheets: ", paste(missing, collapse = ", "))

  markers_all <- imap_dfr(zone_sheets, \(sheets, zone_name) {
    map_dfr(sheets, \(sh) read_one_sheet(sh, zone_name, xlsx))
  })

  markers_clean <- markers_all %>%
    mutate(
      p_val_adj  = ifelse(is.na(p_val_adj), 1, p_val_adj),
      avg_log2fc = ifelse(is.na(avg_log2fc), 0, avg_log2fc),
      pct_1      = ifelse(is.na(pct_1), 0, pct_1),
      pct_2      = ifelse(is.na(pct_2), 0, pct_2),
      delta_pct  = pct_1 - pct_2,
      score      = avg_log2fc * delta_pct
    ) %>%
    filter(delta_pct >= MIN_DELTA_PCT) %>%
    filter(zone != "Vascular" | (pct_1 >= MIN_PCT1 & avg_log2fc >= MIN_LOG2FC & pct_2 <= MAX_PCT2))

  metric_sym <- rlang::sym(USE_METRIC)
  top_per_cluster <- markers_clean %>%
    group_by(zone, cluster) %>%
    arrange(p_val_adj, desc(!!metric_sym), desc(avg_log2fc), desc(delta_pct)) %>%
    slice_head(n = TOP_N_PER_CLUSTER) %>%
    ungroup()

  zone_gene_sets <- top_per_cluster %>%
    group_by(zone) %>%
    summarise(genes = list(sort(unique(gene))), .groups = "drop")

  zone_genes <- setNames(zone_gene_sets$genes, zone_gene_sets$zone)
  zone_genes <- lapply(zone_genes, toupper)
  zone_genes$Vascular <- intersect(zone_genes$Vascular, endo_markers)

  list(zone_genes = zone_genes)
}

map_genes_to_vis <- function(vis_genes, gene_set_upper) {
  vis_genes_upper <- toupper(vis_genes)
  hit_upper <- intersect(gene_set_upper, vis_genes_upper)
  vis_genes[match(hit_upper, vis_genes_upper)]
}

compute_neighbor_fractions <- function(df, zone_col = "zone_call", k = 6) {
  keep <- is.finite(df$pxl_col_in_hires) & is.finite(df$pxl_row_in_hires)
  out <- df %>% mutate(
    frac_epithelial_nbr = NA_real_,
    frac_mesenchymal_nbr = NA_real_,
    frac_vascular_nbr = NA_real_,
    frac_myeloid_nbr = NA_real_
  )
  if (sum(keep) < 2) return(out)

  xy <- as.matrix(df[keep, c("pxl_col_in_hires", "pxl_row_in_hires")])
  kk <- min(k, nrow(xy) - 1)
  if (kk < 1) return(out)
  kn <- dbscan::kNN(xy, k = kk)$id
  sub <- df[keep, , drop = FALSE]

  get_frac <- function(z) {
    vapply(seq_len(nrow(sub)), function(i) mean(sub[[zone_col]][kn[i, ]] == z, na.rm = TRUE), numeric(1))
  }

  out$frac_epithelial_nbr[keep] <- get_frac("Epithelial")
  out$frac_mesenchymal_nbr[keep] <- get_frac("Mesenchymal")
  out$frac_vascular_nbr[keep] <- get_frac("Vascular")
  out$frac_myeloid_nbr[keep] <- get_frac("Myeloid")
  out
}

plot_he_background <- function(df, title_txt) {
  img_path <- unique(df$img_file)
  if (length(img_path) != 1 || is.na(img_path) || !file.exists(img_path)) {
    stop("Missing or invalid img_file for section ", unique(df$section_id))
  }
  img <- png::readPNG(img_path)
  img_h <- nrow(img)
  img_w <- ncol(img)
  plot(NA, xlim = c(0, img_w), ylim = c(0, img_h), xaxs = "i", yaxs = "i",
       asp = 1, xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n",
       main = title_txt)
  rasterImage(as.raster(img), 0, 0, img_w, img_h)
  invisible(list(img_w = img_w, img_h = img_h))
}

save_zone_overlay_base <- function(df, sec, outfile, mode = c("zones","zones_hotspots","mesenchymal_barrier","myeloid_barrier","epithelial_exclusion")) {
  mode <- match.arg(mode)
  png(outfile, width = 2200, height = 2200, res = 300, bg = "white")
  op <- par(no.readonly = TRUE)
  on.exit({par(op); dev.off()}, add = TRUE)
  par(mar = c(1,1,3,1))

  dims <- plot_he_background(df, paste0(sec, " — ", mode))

  d <- df %>% mutate(
    x_plot = pxl_col_in_hires,
    y_plot = dims$img_h - pxl_row_in_hires
  )
  d0 <- d %>% filter(is.finite(x_plot), is.finite(y_plot), in_tissue == 0)
  d1 <- d %>% filter(is.finite(x_plot), is.finite(y_plot), in_tissue == 1)

  # off-tissue background rings
  if (nrow(d0) > 0) {
    points(d0$x_plot, d0$y_plot, pch = 1, cex = 0.90, col = "grey70", lwd = 0.8)
  }

  if (mode == "zones") {
    cols_now <- zone_cols[as.character(d1$zone_call)]
    cols_now[is.na(cols_now)] <- zone_cols[["Uncertain"]]
    if (nrow(d1) > 0) {
      points(d1$x_plot, d1$y_plot, pch = 16, cex = 0.72, col = adjustcolor(cols_now, alpha.f = 0.9))
    }
  }

  # zones + hotspot rings (no normal-brain layer)
  if (mode == "zones_hotspots") {
    cols_now <- zone_cols[as.character(d1$zone_call)]
    cols_now[is.na(cols_now)] <- zone_cols[["Uncertain"]]
    if (nrow(d1) > 0) {
      points(d1$x_plot, d1$y_plot, pch = 16, cex = 0.68, col = adjustcolor(cols_now, alpha.f = 0.88))
    }
    dh <- d1 %>% filter(is_hotspot_semla %in% TRUE)
    if (nrow(dh) > 0) points(dh$x_plot, dh$y_plot, pch = 1, cex = 1.00, col = "white", lwd = 1.6)
  }

  if (mode %in% c("mesenchymal_barrier","myeloid_barrier","epithelial_exclusion")) {
    val <- switch(mode,
      mesenchymal_barrier = d1$frac_mesenchymal_nbr,
      myeloid_barrier     = d1$frac_myeloid_nbr,
      epithelial_exclusion  = d1$frac_epithelial_nbr
    )
    pal <- colorRampPalette(c("#440154", "#31688e", "#35b779", "#fde725"))(100)
    if (all(!is.finite(val)) || length(unique(stats::na.omit(val))) < 2) {
      cols_val <- rep(pal[1], length(val))
    } else {
      brk <- seq(min(val, na.rm = TRUE), max(val, na.rm = TRUE), length.out = 101)
      idx <- cut(val, breaks = brk, include.lowest = TRUE, labels = FALSE)
      cols_val <- pal[pmax(1, idx)]
    }
    if (nrow(d1) > 0) points(d1$x_plot, d1$y_plot, pch = 16, cex = 0.72, col = cols_val)
    dh <- d1 %>% filter(is_hotspot_semla %in% TRUE)
    if (nrow(dh) > 0) points(dh$x_plot, dh$y_plot, pch = 1, cex = 1.00, col = "white", lwd = 1.6)
  }
}

sig <- build_zone_gene_sets(xlsx)
zone_genes <- sig$zone_genes
all_tables <- list()

for (sec in section_keep) {
  message("\n==============================")
  message("Processing section ", sec)
  message("==============================")

  raw_vis_file <- file.path(raw_vis_dir, paste0("vis_section_", sec, "_raw.rds"))
  raw_spot_file <- file.path(raw_spot_dir, paste0("spot_table_", sec, "_raw.tsv"))
  hotspot_file <- file.path(hotspot_dir, paste0("NNLS_section_", sec, ".rds"))

  message("  phase 0 raw vis file: ", raw_vis_file)
  if (!file.exists(raw_vis_file)) stop("Missing raw vis file: ", raw_vis_file)
  if (!file.exists(raw_spot_file)) stop("Missing raw spot table: ", raw_spot_file)
  has_hotspots <- file.exists(hotspot_file)
  if (!has_hotspots) {
    if (!ALL_SECTIONS_RUN) {
      stop("Missing corrected hotspot table: ", hotspot_file)
    }
    message("  no hotspot table for section ", sec,
            " - hotspot context left empty (expected outside GOOD_SECTIONS)")
  }

  vis_raw <- readRDS(raw_vis_file)
  expr <- read_expr_matrix(vis_raw)
  vis_raw <- expr$obj

  spot_tbl <- read_tsv(raw_spot_file, show_col_types = FALSE) %>%
    mutate(barcode16_join = coalesce(as.character(barcode16), extract_barcode16(cell)),
           in_tissue = as.integer(in_tissue))

  hotspot_tbl <- if (has_hotspots) {
    readRDS(hotspot_file) %>%
      mutate(barcode16_join = coalesce(as.character(barcode16), extract_barcode16(cell)),
             is_hotspot_semla = as.logical(is_hotspot_semla)) %>%
      select(any_of(c(
        "barcode16_join", "cell", "Lymphocytes", "local_lymph_score", "is_hotspot_semla"
      )))
  } else {
    tibble(
      barcode16_join   = character(0),
      Lymphocytes      = numeric(0),
      local_lymph_score = numeric(0),
      is_hotspot_semla = logical(0)
    )
  }

  cells_md <- tibble(
    cell = colnames(vis_raw),
    barcode16_join = coalesce(as.character(vis_raw@meta.data$barcode16), extract_barcode16(colnames(vis_raw)))
  )

  tissue_cells <- cells_md %>%
    left_join(spot_tbl %>% select(barcode16_join, in_tissue), by = "barcode16_join") %>%
    filter(in_tissue == 1) %>%
    pull(cell)

  vis_raw <- subset(vis_raw, cells = tissue_cells)
  DefaultAssay(vis_raw) <- "Spatial"

  vis_genes <- rownames(vis_raw)
  zone_genes_vis <- lapply(zone_genes, function(gs) map_genes_to_vis(vis_genes, gs))

  vis_raw <- safe_module_score(vis_raw, list(zone_genes_vis$Epithelial),  "Epithelial")
  vis_raw <- safe_module_score(vis_raw, list(zone_genes_vis$Mesenchymal), "Mesenchymal")
  vis_raw <- safe_module_score(vis_raw, list(zone_genes_vis$Vascular),    "Vascular")
  vis_raw <- safe_module_score(vis_raw, list(zone_genes_vis$Myeloid),     "Myeloid")
  vis_raw$Epithelial_score  <- vis_raw$Epithelial1
  vis_raw$Mesenchymal_score <- vis_raw$Mesenchymal1
  vis_raw$Vascular_score    <- vis_raw$Vascular1
  vis_raw$Myeloid_score     <- vis_raw$Myeloid1

  vis_raw@meta.data <- vis_raw@meta.data %>% select(-any_of(c("Epithelial1","Mesenchymal1","Vascular1","Myeloid1")))

  md <- vis_raw@meta.data %>%
    rownames_to_column("cell") %>%
    mutate(section_id = sec,
           barcode16_join = coalesce(as.character(barcode16), extract_barcode16(cell))) %>%
    select(cell, section_id, barcode16_join,
           Epithelial_score, Mesenchymal_score, Vascular_score, Myeloid_score) %>%
    left_join(spot_tbl %>%
                distinct(barcode16_join, .keep_all = TRUE) %>%
                select(any_of(c(
                  "barcode16_join", "in_tissue", "img_file", "json_file",
                  "pxl_col_in_fullres", "pxl_row_in_fullres",
                  "pxl_col_in_hires", "pxl_row_in_hires",
                  "tissue_hires_scalef", "img_width_hires", "img_height_hires"
                ))),
              by = "barcode16_join") %>%
    left_join(hotspot_tbl %>% distinct(barcode16_join, .keep_all = TRUE), by = "barcode16_join") %>%
    mutate(
      in_tissue = as.integer(in_tissue),
      is_hotspot_semla = ifelse(is.na(is_hotspot_semla), FALSE, is_hotspot_semla)
    )

  md <- md %>%
    mutate(
      Epithelial_score_z  = z_by_section(Epithelial_score,  section_id),
      Mesenchymal_score_z = z_by_section(Mesenchymal_score, section_id),
      Vascular_score_z    = z_by_section(Vascular_score,    section_id),
      Myeloid_score_z     = z_by_section(Myeloid_score,     section_id)
    )

  score_cols_z <- c("Epithelial_score_z","Mesenchymal_score_z","Vascular_score_z","Myeloid_score_z")
  score_mat <- as.matrix(md[, score_cols_z, drop = FALSE])
  winner <- apply(score_mat, 1, function(v) if (all(is.na(v))) NA_character_ else score_cols_z[which.max(v)])
  zone_call <- dplyr::recode(winner,
    "Epithelial_score_z"  = "Epithelial",
    "Mesenchymal_score_z" = "Mesenchymal",
    "Vascular_score_z"    = "Vascular",
    "Myeloid_score_z"     = "Myeloid",
    .default = NA_character_
  )
  maxz <- apply(score_mat, 1, function(v) suppressWarnings(max(v, na.rm = TRUE)))
  zone_call[maxz < min_max_z] <- "Uncertain"
  md$zone_call <- zone_call
  md <- compute_neighbor_fractions(md, zone_col = "zone_call", k = k_neighbors)

  write_tsv(md, file.path(out_tables, paste0("zones_section_", sec, ".tsv")))

  md2 <- md
  rownames(md2) <- md2$cell
  add_cols <- setdiff(colnames(md2), colnames(vis_raw@meta.data))
  vis_raw@meta.data <- cbind(vis_raw@meta.data[rownames(md2), , drop = FALSE], md2[rownames(vis_raw@meta.data), add_cols, drop = FALSE])
  saveRDS(vis_raw, file.path(out_rds, paste0("zones_section_", sec, ".rds")))

  save_zone_overlay_base(md, sec, file.path(out_plots, paste0("B_", sec, "_zones.png")), "zones")
  save_zone_overlay_base(md, sec, file.path(out_plots, paste0("D_", sec, "_zones_hotspots.png")), "zones_hotspots")
  save_zone_overlay_base(md, sec, file.path(out_plots, paste0("E_", sec, "_mesenchymal_barrier.png")), "mesenchymal_barrier")
  save_zone_overlay_base(md, sec, file.path(out_plots, paste0("F_", sec, "_myeloid_barrier.png")), "myeloid_barrier")
  save_zone_overlay_base(md, sec, file.path(out_plots, paste0("G_", sec, "_epithelial_exclusion.png")), "epithelial_exclusion")

  all_tables[[sec]] <- md
}

all_df <- bind_rows(all_tables)
write_tsv(all_df, file.path(out_tables, "zones_all_sections.tsv"))

zone_stats <- all_df %>%
  group_by(section_id, zone_call) %>%
  summarise(n = n(), n_hotspot = sum(is_hotspot_semla %in% TRUE, na.rm = TRUE), pct_hotspot = 100 * n_hotspot / n(), .groups = "drop")
write_tsv(zone_stats, file.path(out_tables, "zone_hotspot_stats.tsv"))

barrier_stats <- all_df %>%
  mutate(hotspot = is_hotspot_semla %in% TRUE) %>%
  group_by(section_id, hotspot) %>%
  summarise(
    mean_frac_epithelial_nbr = mean(frac_epithelial_nbr, na.rm = TRUE),
    mean_frac_mesenchymal_nbr = mean(frac_mesenchymal_nbr, na.rm = TRUE),
    mean_frac_vascular_nbr = mean(frac_vascular_nbr, na.rm = TRUE),
    mean_frac_myeloid_nbr = mean(frac_myeloid_nbr, na.rm = TRUE),
    .groups = "drop"
  )
write_tsv(barrier_stats, file.path(out_tables, "barrier_neighbor_stats.tsv"))

message("Done. Outputs written to: ", out_root)
