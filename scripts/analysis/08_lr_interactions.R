# =============================================================================
# 08_lr_interactions.R
#
# Purpose
# -------
# Ligand-receptor analysis across CellChatDB, CellPhoneDB and NicheNet.
#
# Uses only the corrected backbone: expression from the stage 00 raw
# per-section objects, hotspot / zone / geometry from the stage 06 per-section
# tables. No vis_good, visHE or merged-object dependency.
#
# SECTION SET
# This stage uses FIVE sections, not the usual six: 723 is absent. Confirm and
# record the reason -- every other stage uses GOOD_SECTIONS. The literal is kept
# explicit below rather than substituting GOOD_SECTIONS, so the difference stays
# visible instead of being silently unified.
#
# Zones are additionally excluded per section where coverage is inadequate
# (see zone_exclude_sections).
#
# Inputs
# ------
#   stage_dir("raw_vis")   vis_section_<sec>_raw.rds
#   stage_dir("zones")     tables/
#   CellChatDB, CellPhoneDB (from packages)
#   NicheNet ligand_target_matrix / lr_network / weighted_networks (downloaded)
#     -- record the NicheNet release version; see docs/data_sources.md
#
# Outputs
# -------
#   stage_dir("lr_interactions")
#     qc/, resources/, observed_scores/, top_pairs/,
#     permutations/, nichenet/, final_pair_sets/,
#     rds/phase08_lr_results_bundle.rds
#
# Stochastic
# ----------
#   YES -- use_seed("lr_permutation") at the start of run_resource_permutations.
#   The hotspot-label permutation null (n_perm draws per zone per section) was
#   previously UNSEEDED, so published significance calls could not be
#   reproduced. Re-run this stage after seeding and check that the calls are
#   unchanged before relying on them.
#
# Runtime
# -------
#   ~3 hours
# =============================================================================

source(here::here("config.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(FNN)
  library(dplyr)
  library(readr)
  library(tibble)
  library(purrr)
  library(tidyr)
})


# -------------------------
# Configuration
# -------------------------
raw_vis_dir    <- stage_dir("raw_vis")
phase06_tables <- file.path(stage_dir("zones"), "tables")

out_root   <- stage_dir("lr_interactions")
out_qc     <- ensure_dir(file.path(out_root, "qc"))
out_res    <- ensure_dir(file.path(out_root, "resources"))
out_obs    <- ensure_dir(file.path(out_root, "observed_scores"))
out_top    <- ensure_dir(file.path(out_root, "top_pairs"))
out_perm   <- ensure_dir(file.path(out_root, "permutations"))
out_nichen <- ensure_dir(file.path(out_root, "nichenet"))
out_final  <- ensure_dir(file.path(out_root, "final_pair_sets"))
out_rds    <- ensure_dir(file.path(out_root, "rds"))

# NOTE: five sections, not GOOD_SECTIONS. 723 is deliberately absent here.
# See header. Do not substitute GOOD_SECTIONS without checking why.
good_sections <- c("459", "812", "821", "928", "1239")

zones_to_test <- c("Myeloid", "Mesenchymal", "Vascular")
zone_score_col <- c(
  Myeloid = "Myeloid_score",
  Mesenchymal = "Mesenchymal_score",
  Vascular = "Vascular_score"
)

zone_exclude_sections <- list(
  Myeloid = c("1239"),
  Mesenchymal = character(0),
  Vascular = c("928")
)

# Build hotspot -> kNN neighbors in corrected pixel space
k_neighbors <- THRESH$hotspot_knn

# Top interactions retained per resource / zone after cross-section aggregation
n_top_pairs <- 30

# Permutation depth for top-pair validation
n_perm <- 200

# Run toggles
run_cellchatdb  <- TRUE
run_cellphonedb <- TRUE
run_nichenet    <- TRUE

# NicheNet settings
nichenet_min_hotspots_per_group <- 15
nichenet_consensus_min_sections <- 1
nichenet_top_genes <- 300

# Fallback curated ligand set if no top-pair ligands are available
fallback_ligands <- c(
  "SPP1","MIF","FN1","COL6A1","COL6A2","COL4A1",
  "APP","PPIA","C3","ANXA1","APOE"
)

# -------------------------
# HELPERS
# -------------------------
`%||%` <- function(x, y) if (is.null(x)) y else x

split_entity <- function(x) {
  x <- toupper(as.character(x))
  if (is.na(x) || x == "") return(character(0))
  unlist(strsplit(x, "[_&|:+]", perl = TRUE), use.names = FALSE)
}

read_expr_matrix <- function(obj) {
  DefaultAssay(obj) <- "Spatial"
  mat <- tryCatch(
    GetAssayData(obj, assay = "Spatial", layer = "data"),
    error = function(e) NULL
  )

  if (is.null(mat) || nrow(mat) == 0 || ncol(mat) == 0) {
    message("  no usable 'data' layer found; normalising from counts")
    obj <- NormalizeData(
      obj,
      normalization.method = "LogNormalize",
      scale.factor = 1e4,
      verbose = FALSE
    )
    mat <- tryCatch(
      GetAssayData(obj, assay = "Spatial", layer = "data"),
      error = function(e) NULL
    )
  }

  if (is.null(mat) || nrow(mat) == 0 || ncol(mat) == 0) {
    stop("Could not recover a usable Spatial data layer.")
  }

  list(obj = obj, mat = mat)
}

safe_read_phase06 <- function(sec) {
  f <- file.path(phase06_tables, paste0("zones_section_", sec, ".tsv"))
  if (!file.exists(f)) stop("Missing phase06 table: ", f)
  read_tsv(f, show_col_types = FALSE)
}

safe_read_vis <- function(sec) {
  f <- file.path(raw_vis_dir, paste0("vis_section_", sec, "_raw.rds"))
  if (!file.exists(f)) stop("Missing raw vis object: ", f)
  readRDS(f)
}

get_cell_names <- function(obj) colnames(obj)

find_cell_col <- function(df) {
  cand <- intersect(c("cell","barcode","barcode_raw","spot","spot_id","Cell","Barcode"), colnames(df))
  if (length(cand) == 0) return(NA_character_)
  cand[1]
}

prepare_phase06_table <- function(df, obj, sec) {
  cell_col <- find_cell_col(df)
  obj_cells <- get_cell_names(obj)

  if (!is.na(cell_col)) {
    df2 <- df %>% mutate(cell = as.character(.data[[cell_col]]))
    n_overlap <- sum(df2$cell %in% obj_cells, na.rm = TRUE)
    if (n_overlap > 0) {
      return(df2 %>% filter(cell %in% obj_cells) %>% distinct(cell, .keep_all = TRUE))
    }
  }

  if (nrow(df) != length(obj_cells)) {
    stop(
      "Phase06 table for section ", sec,
      " has no usable cell column and row count does not match object cells. ",
      "nrow(table)=", nrow(df), " ; ncol(object)=", length(obj_cells)
    )
  }

  message("  falling back to row-order matching for section ", sec)
  df %>% mutate(cell = obj_cells)
}

ensure_required_phase06_cols <- function(df, sec) {
  req <- c("cell","section_id","in_tissue","is_hotspot_semla","zone_call",
           "pxl_col_in_hires","pxl_row_in_hires")
  miss <- setdiff(req, colnames(df))
  if (length(miss) > 0) stop("Section ", sec, " phase06 table missing columns: ", paste(miss, collapse = ", "))
  df
}

scale_01 <- function(x) {
  x <- as.numeric(x)
  if (all(!is.finite(x))) return(rep(0, length(x)))
  xr <- range(x, na.rm = TRUE)
  if (!is.finite(diff(xr)) || diff(xr) <= 0) return(rep(1, length(x)))
  (x - xr[1]) / (xr[2] - xr[1] + 1e-9)
}

entity_expr <- function(entity, expr_lin, cells) {
  g <- split_entity(entity)
  if (length(g) == 0 || !all(g %in% rownames(expr_lin))) return(rep(NA_real_, length(cells)))
  as.numeric(Matrix::colMeans(expr_lin[g, cells, drop = FALSE]))
}

present_entity <- function(entity, genes) {
  g <- split_entity(entity)
  length(g) > 0 && all(g %in% genes)
}

prefilter_lr_pairs <- function(lr_df, genes_union) {
  keep <- vapply(seq_len(nrow(lr_df)), function(i) {
    present_entity(lr_df$ligand[i], genes_union) && present_entity(lr_df$receptor[i], genes_union)
  }, logical(1))
  lr_df[keep, , drop = FALSE]
}

aggregate_top_pairs <- function(score_df, zone_target, topN = 30) {
  if (is.null(score_df) || nrow(score_df) == 0) return(NULL)

  rank_mode <- if (zone_target == "Vascular") "max" else "median"

  agg <- score_df %>%
    group_by(ligand, receptor) %>%
    summarise(
      n_sections = n_distinct(section_id),
      median_score = median(score, na.rm = TRUE),
      max_score = max(score, na.rm = TRUE),
      mean_score = mean(score, na.rm = TRUE),
      .groups = "drop"
    )

  if (rank_mode == "median") {
    agg <- agg %>% arrange(desc(n_sections), desc(median_score), desc(max_score))
  } else {
    agg <- agg %>% arrange(desc(max_score), desc(n_sections), desc(median_score))
  }

  head(agg, topN)
}

# -------------------------
# RESOURCE LOADERS
# -------------------------
load_cellchat_resource <- function() {
  tmp_rda <- tempfile(fileext = ".rda")
  url <- "https://github.com/jinworks/CellChat/raw/main/data/CellChatDB.human.rda"
  utils::download.file(url, tmp_rda, mode = "wb")
  e <- new.env(parent = emptyenv())
  load(tmp_rda, envir = e)

  if (!exists("CellChatDB.human", envir = e)) {
    stop("Failed to load CellChatDB.human from downloaded .rda")
  }

  lr <- unique(e$CellChatDB.human$interaction[, c("ligand","receptor","annotation","pathway_name")])
  lr <- lr %>%
    transmute(
      ligand = toupper(as.character(ligand)),
      receptor = toupper(as.character(receptor)),
      class = as.character(annotation),
      pathway = as.character(pathway_name)
    ) %>%
    filter(class %in% c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact")) %>%
    distinct()

  lr
}

load_cellphonedb_resource <- function() {
  if (!requireNamespace("liana", quietly = TRUE)) {
    warning("Package 'liana' not installed; skipping CellPhoneDB resource.")
    return(NULL)
  }

  cpdb <- liana::select_resource("CellPhoneDB")[[1]]
  col_lig <- if ("ligand" %in% names(cpdb)) "ligand" else if ("source_genesymbol" %in% names(cpdb)) "source_genesymbol" else NA_character_
  col_rec <- if ("receptor" %in% names(cpdb)) "receptor" else if ("target_genesymbol" %in% names(cpdb)) "target_genesymbol" else NA_character_

  if (is.na(col_lig) || is.na(col_rec)) {
    warning("Could not find ligand/receptor columns in CellPhoneDB resource; skipping.")
    return(NULL)
  }

  cpdb %>%
    transmute(
      ligand = toupper(as.character(.data[[col_lig]])),
      receptor = toupper(as.character(.data[[col_rec]])),
      class = "CellPhoneDB",
      pathway = NA_character_
    ) %>%
    filter(!is.na(ligand), !is.na(receptor), ligand != "", receptor != "") %>%
    distinct()
}

# -------------------------
# NicheNet resources
# -------------------------
organism_nichenet <- "human"
nichenet_resource_dir <- file.path(out_res, "nichenet")
dir.create(nichenet_resource_dir, recursive = TRUE, showWarnings = FALSE)

load_nichenet_resources <- function(
  organism = organism_nichenet,
  resource_dir = nichenet_resource_dir,
  download_if_missing = FALSE
) {
  if (!requireNamespace("nichenetr", quietly = TRUE)) {
    warning("Package 'nichenetr' not installed; skipping NicheNet.")
    return(NULL)
  }

  files <- switch(
    organism,
    human = list(
      ligand_target_matrix = file.path(resource_dir, "ligand_target_matrix_nsga2r_final.rds"),
      lr_network           = file.path(resource_dir, "lr_network_human_21122021.rds"),
      weighted_networks    = file.path(resource_dir, "weighted_networks_nsga2r_final.rds")
    ),
    mouse = list(
      ligand_target_matrix = file.path(resource_dir, "ligand_target_matrix_nsga2r_final_mouse.rds"),
      lr_network           = file.path(resource_dir, "lr_network_mouse_21122021.rds"),
      weighted_networks    = file.path(resource_dir, "weighted_networks_nsga2r_final_mouse.rds")
    ),
    stop("organism must be 'human' or 'mouse'")
  )

  missing_local <- names(files)[!file.exists(unlist(files))]
  if (length(missing_local) > 0) {
    stop(
      "Missing NicheNet resource files: ",
      paste(missing_local, collapse = ", "),
      "\nExpected in: ", resource_dir
    )
  }

  ligand_target_matrix <- readRDS(files$ligand_target_matrix)
  lr_network <- readRDS(files$lr_network)
  weighted_networks <- readRDS(files$weighted_networks)

  message(
    "Loaded NicheNet resources: ",
    nrow(ligand_target_matrix), " targets x ",
    ncol(ligand_target_matrix), " ligands"
  )

  list(
    ligand_target_matrix = ligand_target_matrix,
    lr_network = lr_network,
    weighted_networks = weighted_networks,
    files = files,
    organism = organism
  )
}

# -------------------------
# BUILD PER-SECTION BUNDLES
# -------------------------
section_bundles <- list()
section_qc <- list()
all_genes_union <- character()

for (sec in good_sections) {
  message("\n==============================")
  message("Preparing section ", sec)
  message("==============================")

  vis_sec <- safe_read_vis(sec)
  phase06 <- safe_read_phase06(sec)

  expr_res <- read_expr_matrix(vis_sec)
  vis_sec <- expr_res$obj
  expr_log <- expr_res$mat
  phase06 <- prepare_phase06_table(phase06, vis_sec, sec)
  phase06 <- ensure_required_phase06_cols(phase06, sec)

  md <- phase06 %>%
    mutate(
      cell = as.character(cell),
      section_id = as.character(section_id),
      in_tissue = as.integer(in_tissue),
      is_hotspot_semla = ifelse(is.na(is_hotspot_semla), FALSE, as.logical(is_hotspot_semla)),
      cell_uid = paste(section_id, cell, sep = "__")
    )

  # keep only cells shared with expression object
  md <- md %>% filter(cell %in% colnames(expr_log))
  expr_log <- expr_log[, md$cell, drop = FALSE]

  # tumor / tissue spots only
  if ("normal_brain_candidate" %in% colnames(md)) {
    md <- md %>% filter(in_tissue == 1, !(normal_brain_candidate %in% TRUE))
  } else {
    md <- md %>% filter(in_tissue == 1)
  }
  expr_log <- expr_log[, md$cell, drop = FALSE]

  # corrected geometry required for neighbors
  keep_geom <- is.finite(md$pxl_col_in_hires) & is.finite(md$pxl_row_in_hires)
  md <- md[keep_geom, , drop = FALSE]
  expr_log <- expr_log[, md$cell, drop = FALSE]

  # rename expression columns to unique ids for safe cross-section work
  colnames(expr_log) <- md$cell_uid
  md$cell_uid <- colnames(expr_log)

  expr_lin <- expm1(expr_log)

  section_bundles[[sec]] <- list(
    section_id = sec,
    md = md,
    expr_log = expr_log,
    expr_lin = expr_lin
  )

  all_genes_union <- union(all_genes_union, rownames(expr_log))

  section_qc[[sec]] <- tibble(
    section_id = sec,
    n_spots = nrow(md),
    n_hotspots = sum(md$is_hotspot_semla, na.rm = TRUE),
    n_myeloid = sum(md$zone_call == "Myeloid", na.rm = TRUE),
    n_mesenchymal = sum(md$zone_call == "Mesenchymal", na.rm = TRUE),
    n_vascular = sum(md$zone_call == "Vascular", na.rm = TRUE)
  )
}

section_qc_df <- bind_rows(section_qc)
write_tsv(section_qc_df, file.path(out_qc, "section_input_qc.tsv"))
saveRDS(section_bundles, file.path(out_rds, "section_bundles.rds"))

# -------------------------
# EDGE BUILDING
# -------------------------
make_edges_one_section <- function(bundle, k = 6) {
  md <- bundle$md
  coords <- as.matrix(md[, c("pxl_col_in_hires", "pxl_row_in_hires")])
  if (nrow(coords) <= k) return(NULL)

  hs_idx <- which(md$is_hotspot_semla %in% TRUE)
  if (length(hs_idx) == 0) return(NULL)

  nn <- FNN::get.knn(coords, k = k)
  edges <- do.call(rbind, lapply(hs_idx, function(i) {
    nei <- nn$nn.index[i, ]
    data.frame(
      from = md$cell_uid[i],
      to = md$cell_uid[nei],
      section_id = bundle$section_id,
      stringsAsFactors = FALSE
    )
  }))

  edges$to_zone <- md$zone_call[match(edges$to, md$cell_uid)]
  edges
}

edges_all <- bind_rows(lapply(section_bundles, make_edges_one_section, k = k_neighbors))
write_tsv(edges_all, file.path(out_qc, "hotspot_neighbor_edges_all.tsv"))

hotspots_touch_zone <- function(zone_name, edges_df = edges_all) {
  unique(edges_df$from[edges_df$to_zone == zone_name])
}

# -------------------------
# OBSERVED SCORE ENGINE
# -------------------------
score_zone_to_hotspot_one_section <- function(bundle, zone_target, lr_pairs, k = 6) {
  edges <- make_edges_one_section(bundle, k = k)
  if (is.null(edges) || nrow(edges) == 0) return(NULL)

  edges <- edges[edges$to_zone == zone_target, , drop = FALSE]
  if (nrow(edges) == 0) return(NULL)

  # reverse: zone sender -> hotspot receiver
  edges_rev <- transform(edges, from = to, to = from)

  md <- bundle$md
  cells <- md$cell_uid
  expr_lin <- bundle$expr_lin

  Lw <- if ("Lymphocytes" %in% colnames(md)) scale_01(md$Lymphocytes) else rep(1, nrow(md))
  names(Lw) <- md$cell_uid

  zcol <- zone_score_col[[zone_target]]
  z_w <- if (!is.null(zcol) && zcol %in% colnames(md)) pmax(md[[zcol]], 0) else as.numeric(md$zone_call == zone_target)
  names(z_w) <- md$cell_uid

  ligands_u <- unique(lr_pairs$ligand)
  receptors_u <- unique(lr_pairs$receptor)

  Lmat <- vapply(ligands_u, entity_expr, numeric(length(cells)), expr_lin = expr_lin, cells = cells)
  Rmat <- vapply(receptors_u, entity_expr, numeric(length(cells)), expr_lin = expr_lin, cells = cells)
  Lmat <- t(Lmat); rownames(Lmat) <- ligands_u
  Rmat <- t(Rmat); rownames(Rmat) <- receptors_u

  idx_from <- match(edges_rev$from, cells)
  idx_to   <- match(edges_rev$to, cells)

  lig_idx <- match(lr_pairs$ligand, ligands_u)
  rec_idx <- match(lr_pairs$receptor, receptors_u)

  tmp <- Lmat[lig_idx, idx_from, drop = FALSE] * Rmat[rec_idx, idx_to, drop = FALSE]
  tmp <- sweep(tmp, 2, z_w[edges_rev$from] * Lw[edges_rev$to], "*")
  score <- rowMeans(tmp, na.rm = TRUE)

  tibble(
    section_id = bundle$section_id,
    zone = zone_target,
    ligand = lr_pairs$ligand,
    receptor = lr_pairs$receptor,
    class = lr_pairs$class %||% NA_character_,
    pathway = lr_pairs$pathway %||% NA_character_,
    score = score
  ) %>%
    filter(is.finite(score), score > 0)
}

run_resource_observed <- function(resource_name, lr_pairs) {
  lr_pairs <- prefilter_lr_pairs(lr_pairs, all_genes_union)
  write_tsv(lr_pairs, file.path(out_res, paste0(resource_name, "_pairs_filtered.tsv")))

  message("\n", resource_name, " detectable pairs: ", nrow(lr_pairs))
  if (nrow(lr_pairs) == 0) return(NULL)

  obs_out <- list()
  top_out <- list()

  for (zone in zones_to_test) {
    secs_use <- setdiff(good_sections, zone_exclude_sections[[zone]] %||% character(0))

    res_zone <- lapply(secs_use, function(sec) {
      message("Observed scoring | ", resource_name, " | ", zone, " | section ", sec)
      score_zone_to_hotspot_one_section(section_bundles[[sec]], zone, lr_pairs, k = k_neighbors)
    }) %>% bind_rows()

    write_tsv(res_zone, file.path(out_obs, paste0(resource_name, "_", zone, "_observed.tsv")))
    obs_out[[zone]] <- res_zone

    top_zone <- aggregate_top_pairs(res_zone, zone_target = zone, topN = n_top_pairs)
    if (!is.null(top_zone)) {
      top_zone <- top_zone %>% mutate(resource = resource_name, zone = zone)
      write_tsv(top_zone, file.path(out_top, paste0(resource_name, "_", zone, "_top", n_top_pairs, ".tsv")))
    }
    top_out[[zone]] <- top_zone
  }

  list(observed = obs_out, top = top_out)
}

# -------------------------
# PERMUTATION ENGINE
# -------------------------
permute_scores_zone_one_section <- function(bundle, zone_target, lr_pairs, nperm = 200, k = 6, verbose = TRUE) {
  md <- bundle$md
  cells <- md$cell_uid
  coords <- as.matrix(md[, c("pxl_col_in_hires", "pxl_row_in_hires")])
  if (nrow(coords) <= k) return(NULL)

  hot_obs <- md$is_hotspot_semla %in% TRUE
  n_hs <- sum(hot_obs)
  if (n_hs == 0) return(NULL)

  nn <- FNN::get.knn(coords, k = k)

  Lw <- if ("Lymphocytes" %in% colnames(md)) scale_01(md$Lymphocytes) else rep(1, nrow(md))
  zcol <- zone_score_col[[zone_target]]
  z_w <- if (!is.null(zcol) && zcol %in% colnames(md)) pmax(md[[zcol]], 0) else as.numeric(md$zone_call == zone_target)
  zcall <- as.character(md$zone_call)

  expr_lin <- bundle$expr_lin

  ligands_u <- unique(lr_pairs$ligand)
  receptors_u <- unique(lr_pairs$receptor)
  Lmat <- vapply(ligands_u, entity_expr, numeric(length(cells)), expr_lin = expr_lin, cells = cells)
  Rmat <- vapply(receptors_u, entity_expr, numeric(length(cells)), expr_lin = expr_lin, cells = cells)
  Lmat <- t(Lmat); rownames(Lmat) <- ligands_u
  Rmat <- t(Rmat); rownames(Rmat) <- receptors_u
  lig_idx <- match(lr_pairs$ligand, ligands_u)
  rec_idx <- match(lr_pairs$receptor, receptors_u)

  score_given_hot_idx <- function(hs_idx) {
    from_hs <- rep(hs_idx, each = k)
    to_nbr <- as.vector(t(nn$nn.index[hs_idx, , drop = FALSE]))
    from_zone <- to_nbr
    to_hot <- from_hs

    keep_z <- which(zcall[from_zone] == zone_target)
    if (length(keep_z) == 0) return(rep(NA_real_, nrow(lr_pairs)))

    from_zone <- from_zone[keep_z]
    to_hot <- to_hot[keep_z]
    w <- z_w[from_zone] * Lw[to_hot]

    tmp <- Lmat[lig_idx, from_zone, drop = FALSE] * Rmat[rec_idx, to_hot, drop = FALSE]
    tmp <- sweep(tmp, 2, w, "*")
    rowMeans(tmp, na.rm = TRUE)
  }

  obs <- score_given_hot_idx(which(hot_obs))
  null <- matrix(NA_real_, nrow = nperm, ncol = nrow(lr_pairs))

  for (p in seq_len(nperm)) {
    if (verbose && p %% 50 == 0) message("perm ", p, "/", nperm, " | ", bundle$section_id, " | ", zone_target)
    hs_rand <- sample(seq_len(nrow(md)), n_hs, replace = FALSE)
    null[p, ] <- score_given_hot_idx(hs_rand)
  }

  pval <- vapply(seq_along(obs), function(j) {
    nu <- null[, j]
    nu <- nu[is.finite(nu)]
    if (!is.finite(obs[j]) || length(nu) < 10) return(NA_real_)
    (sum(nu >= obs[j]) + 1) / (length(nu) + 1)
  }, numeric(1))

  tibble(
    section_id = bundle$section_id,
    zone = zone_target,
    ligand = lr_pairs$ligand,
    receptor = lr_pairs$receptor,
    obs_score = obs,
    pval = pval
  )
}

run_resource_permutations <- function(resource_name, top_tables) {
  if (length(top_tables) == 0) return(NULL)

  # Reseed at the start of each resource so permutation nulls are
  # reproducible and do not depend on what else consumed the RNG first.
  use_seed("lr_permutation")

  perm_out <- list()
  for (zone in names(top_tables)) {
    top_zone <- top_tables[[zone]]
    if (is.null(top_zone) || nrow(top_zone) == 0) next

    lr_pairs <- top_zone %>% select(ligand, receptor)
    secs_use <- setdiff(good_sections, zone_exclude_sections[[zone]] %||% character(0))

    res <- lapply(secs_use, function(sec) {
      message("Permutation scoring | ", resource_name, " | ", zone, " | section ", sec)
      permute_scores_zone_one_section(section_bundles[[sec]], zone, lr_pairs, nperm = n_perm, k = k_neighbors, verbose = TRUE)
    }) %>% bind_rows()

    write_tsv(res, file.path(out_perm, paste0(resource_name, "_", zone, "_perm_top", n_top_pairs, ".tsv")))
    perm_out[[zone]] <- res
  }

  perm_out
}

# -------------------------
# NICHENET CONSENSUS DE
# -------------------------
run_hotspot_touch_de_one_section <- function(bundle, zone_name, min_cells = 15) {
  edges <- make_edges_one_section(bundle, k = k_neighbors)
  if (is.null(edges) || nrow(edges) == 0) return(NULL)

  hs_touch <- unique(edges$from[edges$to_zone == zone_name])
  hs_all <- bundle$md$cell_uid[bundle$md$is_hotspot_semla %in% TRUE]
  hs_not <- setdiff(hs_all, hs_touch)

  if (length(hs_touch) < min_cells || length(hs_not) < min_cells) return(NULL)

  # use the section object directly; one layer only, so FindMarkers is safe here
  obj <- CreateSeuratObject(counts = bundle$expr_log)
  DefaultAssay(obj) <- "RNA"
  obj <- NormalizeData(obj, verbose = FALSE)
  obj$grp <- ifelse(colnames(obj) %in% hs_touch, "touch", ifelse(colnames(obj) %in% hs_not, "not_touch", NA_character_))
  obj <- subset(obj, cells = colnames(obj)[!is.na(obj$grp)])
  Idents(obj) <- obj$grp

  de <- FindMarkers(
    obj,
    ident.1 = "touch",
    ident.2 = "not_touch",
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0.05,
    assay = "RNA"
  )

  de %>% rownames_to_column("gene") %>% mutate(section_id = bundle$section_id, zone = zone_name)
}

make_consensus_geneset <- function(de_list, min_sections = 2, top_genes = 300) {
  de_all <- bind_rows(de_list)
  if (nrow(de_all) == 0) return(list(genes = character(0), table = de_all))

  de_up <- de_all %>%
    filter(p_val_adj < 0.20, avg_log2FC > 0.10) %>%
    group_by(gene) %>%
    summarise(
      n_sections = n_distinct(section_id),
      mean_log2FC = mean(avg_log2FC, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(n_sections >= min_sections) %>%
    arrange(desc(n_sections), desc(mean_log2FC))

  genes <- head(de_up$gene, top_genes)
  list(genes = genes, table = de_up)
}

auc_from_scores <- function(scores, is_pos) {
  ok <- is.finite(scores)
  scores <- scores[ok]
  is_pos <- is_pos[ok]
  n_pos <- sum(is_pos)
  n_neg <- sum(!is_pos)
  if (n_pos < 10 || n_neg < 10) return(NA_real_)
  r <- rank(scores, ties.method = "average")
  (sum(r[is_pos]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

score_nichenet_ligands <- function(geneset, ligand_candidates, lt_mat, bg_genes) {
  geneset_use <- intersect(geneset, bg_genes)
  if (length(geneset_use) < 15) return(NULL)
  
  out <- lapply(ligand_candidates, function(L) {
    if (!(L %in% colnames(lt_mat))) return(NULL)
    
    s <- as.numeric(lt_mat[bg_genes, L])
    is_pos <- bg_genes %in% geneset_use
    
    tibble(
      ligand = L,
      auc = auc_from_scores(s, is_pos),
      mean_pos = mean(s[is_pos], na.rm = TRUE),
      mean_neg = mean(s[!is_pos], na.rm = TRUE)
    )
  }) %>%
    bind_rows() %>%
    arrange(desc(auc))
  
  out
}



make_joint_lr_shortlists <- function(results_all, zones, out_dir, topN = n_top_pairs) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  if (is.null(results_all$CellChatDB$top) || is.null(results_all$CellPhoneDB$top)) {
    warning("Both CellChatDB and CellPhoneDB top tables are required to build joint shortlists.")
    return(invisible(NULL))
  }

  both_all <- list()
  relaxed_all <- list()

  for (zone in zones) {
    top_cc <- results_all$CellChatDB$top[[zone]]
    top_cp <- results_all$CellPhoneDB$top[[zone]]

    if (is.null(top_cc)) top_cc <- tibble()
    if (is.null(top_cp)) top_cp <- tibble()

    if (nrow(top_cc) > 0) {
      top_cc <- top_cc %>%
        transmute(
          zone = zone,
          ligand,
          receptor,
          score_cellchat = mean_score,
          n_sections_cellchat = n_sections,
          median_cellchat = median_score,
          max_cellchat = max_score
        )
    }
    if (nrow(top_cp) > 0) {
      top_cp <- top_cp %>%
        transmute(
          zone = zone,
          ligand,
          receptor,
          score_cellphonedb = mean_score,
          n_sections_cellphonedb = n_sections,
          median_cellphonedb = median_score,
          max_cellphonedb = max_score
        )
    }

    both <- full_join(top_cc, top_cp, by = c("zone", "ligand", "receptor")) %>%
      mutate(
        in_cellchat = !is.na(score_cellchat),
        in_cellphonedb = !is.na(score_cellphonedb),
        in_both = in_cellchat & in_cellphonedb,
        mean_score_both = rowMeans(cbind(score_cellchat, score_cellphonedb), na.rm = TRUE),
        min_score_both = pmin(score_cellchat, score_cellphonedb, na.rm = TRUE)
      ) %>%
      arrange(desc(in_both), desc(min_score_both), desc(mean_score_both))

    write_tsv(both, file.path(out_dir, paste0("LR_pairs_joint_", zone, ".tsv")))

    map_pairs <- both %>%
      filter(in_both) %>%
      arrange(desc(min_score_both), desc(mean_score_both)) %>%
      slice_head(n = topN)

    write_tsv(map_pairs, file.path(out_dir, paste0("LR_pairs_for_mapping_", zone, "_top", topN, ".tsv")))

    relaxed_pairs <- both %>%
      filter(in_cellchat | in_cellphonedb) %>%
      arrange(desc(in_both), desc(mean_score_both)) %>%
      slice_head(n = topN)

    write_tsv(relaxed_pairs, file.path(out_dir, paste0("LR_pairs_relaxed_", zone, "_top", topN, ".tsv")))

    both_all[[zone]] <- map_pairs
    relaxed_all[[zone]] <- relaxed_pairs
  }

  write_tsv(bind_rows(both_all), file.path(out_dir, paste0("LR_pairs_for_mapping_all_zones_top", topN, ".tsv")))
  write_tsv(bind_rows(relaxed_all), file.path(out_dir, paste0("LR_pairs_relaxed_all_zones_top", topN, ".tsv")))

  invisible(list(mapping = both_all, relaxed = relaxed_all))
}

# -------------------------
# RUN RESOURCES
# -------------------------
results_all <- list()
ligands_for_nichenet <- list()

if (run_cellchatdb) {
  lr_cellchat <- load_cellchat_resource()
  write_tsv(lr_cellchat, file.path(out_res, "CellChatDB_pairs_raw.tsv"))
  results_all$CellChatDB <- run_resource_observed("CellChatDB", lr_cellchat)
  results_all$CellChatDB$perm <- run_resource_permutations("CellChatDB", results_all$CellChatDB$top)
}

if (run_cellphonedb) {
  lr_cpdb <- load_cellphonedb_resource()
  if (!is.null(lr_cpdb)) {
    write_tsv(lr_cpdb, file.path(out_res, "CellPhoneDB_pairs_raw.tsv"))
    results_all$CellPhoneDB <- run_resource_observed("CellPhoneDB", lr_cpdb)
    results_all$CellPhoneDB$perm <- run_resource_permutations("CellPhoneDB", results_all$CellPhoneDB$top)
  }
}

# gather top ligands for optional NicheNet scoring
for (zone in zones_to_test) {
  ligands_for_nichenet[[zone]] <- unique(c(
    if (!is.null(results_all$CellChatDB$top[[zone]])) results_all$CellChatDB$top[[zone]]$ligand else character(0),
    if (!is.null(results_all$CellPhoneDB$top[[zone]])) results_all$CellPhoneDB$top[[zone]]$ligand else character(0)
  ))
  if (length(ligands_for_nichenet[[zone]]) == 0) ligands_for_nichenet[[zone]] <- fallback_ligands
}
write_tsv(bind_rows(lapply(names(ligands_for_nichenet), function(z) {
  tibble(zone = z, ligand = ligands_for_nichenet[[z]])
})), file.path(out_nichen, "candidate_ligands_by_zone.tsv"))

# -------------------------
# RUN NICHENET (optional)
# -------------------------
if (run_nichenet) {
  nn_res <- load_nichenet_resources()
  if (!is.null(nn_res)) {
    lt_mat <- nn_res$ligand_target_matrix

    # ligand_target_matrix:
    # rows = target genes
    # cols = ligands
    bg_genes <- intersect(rownames(lt_mat), all_genes_union)

    for (zone in zones_to_test) {
      secs_use <- setdiff(good_sections, zone_exclude_sections[[zone]] %||% character(0))

      de_list <- setNames(vector("list", length(secs_use)), secs_use)

      for (sec in secs_use) {
        message("NicheNet DE | ", zone, " | section ", sec)

        de_sec <- run_hotspot_touch_de_one_section(
          section_bundles[[sec]],
          zone,
          min_cells = nichenet_min_hotspots_per_group
        )

        de_list[[sec]] <- de_sec

        if (!is.null(de_sec) && nrow(de_sec) > 0) {
          write_tsv(
            de_sec,
            file.path(out_nichen, paste0("DE_touch_vs_not_", zone, "_section_", sec, ".tsv"))
          )
        }
      }

      de_list_nonnull <- Filter(Negate(is.null), de_list)
      de_all <- bind_rows(de_list_nonnull, .id = "section_id_list")

      if (nrow(de_all) > 0) {
        write_tsv(
          de_all,
          file.path(out_nichen, paste0("DE_touch_vs_not_", zone, "_all_sections.tsv"))
        )
      }

      cons <- make_consensus_geneset(
        de_list_nonnull,
        min_sections = nichenet_consensus_min_sections,
        top_genes = nichenet_top_genes
      )

      write_tsv(
        cons$table,
        file.path(out_nichen, paste0("consensus_upgenes_", zone, ".tsv"))
      )
      write_tsv(
        tibble(zone = zone, gene = cons$genes),
        file.path(out_nichen, paste0("geneset_", zone, ".tsv"))
      )

      ligands_use <- unique(intersect(ligands_for_nichenet[[zone]], colnames(lt_mat)))

      message("Zone: ", zone)
      message("Consensus genes: ", length(cons$genes))
      message("Ligands tested: ", length(ligands_use))
      print(head(ligands_use, 10))

      # zone-level NicheNet
      if (length(cons$genes) >= 15 && length(ligands_use) > 0) {
        act <- score_nichenet_ligands(cons$genes, ligands_use, lt_mat, bg_genes)

        if (!is.null(act) && nrow(act) > 0) {
          print(head(act, 10))
          act <- act %>% mutate(zone = zone)
          write_tsv(
            act,
            file.path(out_nichen, paste0("ligand_activity_", zone, ".tsv"))
          )
        }
      } else {
        message("Skipping zone-level NicheNet for ", zone, ": too few consensus genes or ligands.")
      }

      # section-level NicheNet
      section_act_list <- lapply(names(de_list_nonnull), function(sec) {
        de_sec <- de_list_nonnull[[sec]]

        genes_sec <- de_sec %>%
          filter(p_val_adj < 0.20, avg_log2FC > 0.10) %>%
          arrange(p_val_adj, desc(avg_log2FC)) %>%
          pull(gene) %>%
          unique()

        genes_sec <- intersect(genes_sec, bg_genes)

        message(
          "Section-level NicheNet | zone=", zone,
          " | section=", sec,
          " | genes=", length(genes_sec),
          " | ligands=", length(ligands_use)
        )

        if (length(genes_sec) < 15 || length(ligands_use) == 0) {
          return(NULL)
        }

        act_sec <- score_nichenet_ligands(genes_sec, ligands_use, lt_mat, bg_genes)
        if (is.null(act_sec) || nrow(act_sec) == 0) {
          return(NULL)
        }

        act_sec %>%
          mutate(
            zone = zone,
            section_id = sec,
            n_geneset = length(genes_sec)
          )
      })

      section_act <- bind_rows(section_act_list)

      if (nrow(section_act) > 0) {
        write_tsv(
          section_act,
          file.path(out_nichen, paste0("ligand_activity_", zone, "_by_section.tsv"))
        )
      }
    }
  }
}

# -------------------------
# SAVE HIGH-CONFIDENCE LR SHORTLISTS
# -------------------------
if (!is.null(results_all$CellChatDB) && !is.null(results_all$CellPhoneDB)) {
  shortlist_res <- make_joint_lr_shortlists(
    results_all = results_all,
    zones = zones_to_test,
    out_dir = out_final,
    topN = n_top_pairs
  )
  results_all$joint_shortlists <- shortlist_res
}

saveRDS(results_all, file.path(out_rds, "phase08_lr_results_bundle.rds"))
message("\nPhase 08 LR analysis complete. Outputs written to: ", out_root)
