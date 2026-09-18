# =============================================================================
# 07_functional_states.R
#
# Purpose
# -------
# Compute hotspot functional-state scores (UCell) and produce plot-ready
# per-section tables for the Figure 4 scripts.
#
# Geometry, hotspot and zone context come from the stage 06 tables; expression
# comes from the stage 00 raw per-section objects. No visHE dependency and no
# in-memory vis_good dependency.
#
# STATE CALLING
# UCell scores are continuous signature-activity scores. To avoid forcing a
# label when all programs are weak, state calls are made only within lymphocyte
# hotspots and calibrated separately per signature and section. A state is
# eligible when its score is in the top STATE_HOTSPOT_PERCENTILE of hotspot
# scores for that signature in that section AND at least
# MIN_SIGNATURE_GENES_DETECTED signature genes are detected in that spot. If
# more than one state is eligible, the top two percentile scores must differ by
# at least STATE_MIN_MARGIN, otherwise the hotspot is left Unassigned.
#
# A sensitivity grid over percentile and margin is written for QC; the primary
# calls are unaffected by it.
#
# Inputs
# ------
#   stage_dir("zones")     tables/
#   stage_dir("raw_vis")   vis_section_<sec>_raw.rds
#
# Outputs
# -------
#   stage_dir("functional_states")
#     tables/phase07_hotspot_states_all_sections.tsv
#     tables/phase07_hotspot_states_hotspots_only.tsv
#     tables/phase07_state_frequency_by_section.tsv
#     tables/phase07_state_frequency_qc.tsv
#     tables/phase07_signature_gene_summary.tsv
#     sections/phase07_section_<sec>_plotready.{tsv,rds}
#     sections/phase07_section_<sec>_retention_scaled.{tsv,rds}
#     rds/phase07_analysis_bundle.rds
#
# Stochastic
# ----------
#   none -- UCell scoring and the sensitivity grid are deterministic.
#   If a future change introduces sampling here, add a SEEDS entry first.
#
# Runtime
# -------
#   ~30 minutes
# =============================================================================

source(here::here("config.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(UCell)
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
  library(tidyr)
  library(stringr)
})


# -------------------------
# Configuration
# -------------------------
raw_vis_dir    <- stage_dir("raw_vis")
phase06_tables <- file.path(stage_dir("zones"), "tables")

out_root     <- stage_dir("functional_states")
out_tables   <- ensure_dir(file.path(out_root, "tables"))
out_sections <- ensure_dir(file.path(out_root, "sections"))
out_rds      <- ensure_dir(file.path(out_root, "rds"))
out_qc       <- ensure_dir(file.path(out_root, "qc"))

good_sections   <- GOOD_SECTIONS
panelA_sections <- c("459", "928")
panelC_section  <- "459"

# Dominant-state calling -- see header for the rule these implement.
STATE_HOTSPOT_PERCENTILE     <- 0.50
STATE_MIN_MARGIN             <- 0.05
MIN_SIGNATURE_GENES_DETECTED <- 2L

# Sensitivity grid written for QC; primary calls above remain unchanged.
SENSITIVITY_PERCENTILES <- c(0.70, 0.75, 0.80, 0.85, 0.90)
SENSITIVITY_MARGINS     <- c(0.00, 0.05, 0.10, 0.15)

# -------------------------
# Functional-state signatures
# Keep the 5-state dominant-state layout used previously,
# but also compute proliferation as a saved QC score.
# -------------------------
state_sets_requested <- list(
  cyto = c("NKG7","PRF1","GZMB","GNLY","CTSW","KLRD1","KLRK1","XCL1","XCL2","CCL5","IFNG"),
  inflam = c("IFNG","TNF","CCL4","CCL5","XCL1","XCL2","NFKBIA","IRF1","STAT1"),
  exhaust = c("PDCD1","LAG3","TIGIT","HAVCR2","CTLA4","TOX","ENTPD1","CXCL13"),
  il17like = c("RORC","CCR6","KLRB1","IL7R","IL17A","IL17F","AREG"),
  prolif = c("MKI67","TOP2A","TYMS","STMN1","HMGB2","TUBB"),
  retention = c("CD44","ITGAE","CXCR3","CXCR6","CD69","CCL5")
)

state_score_cols_core <- c(
  "state_cyto",
  "state_inflam",
  "state_exhaust",
  "state_il17like",
  "state_retention"
)

pretty_state <- c(
  state_cyto      = "Cytotoxic",
  state_inflam    = "Inflammatory",
  state_exhaust   = "Exhausted",
  state_il17like  = "IL17-like",
  state_retention = "Retention"
)

state_levels <- c("Cytotoxic", "Inflammatory", "Exhausted", "IL17-like", "Retention", "Unassigned")

retention_genes_requested <- c(
  "TNC", "ICAM1", "VCAM1", "ITGA3",
  "PTX3", "CXCL3", "SERPINE1",
  "PPP1R15A", "PDK1", "UPP1", "SLC2A1", "NAMPT", "GPNMB"
)

# -------------------------
# Helpers
# -------------------------
`%||%` <- function(x, y) if (is.null(x)) y else x

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

match_requested_genes <- function(requested, present) {
  present_up <- toupper(present)
  req_up <- unique(toupper(requested))
  idx <- match(req_up, present_up)
  idx <- idx[!is.na(idx)]
  unique(present[idx])
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

get_cell_names <- function(obj) {
  colnames(obj)
}

find_cell_col <- function(df) {
  cand <- intersect(
    c("cell", "barcode", "barcode_raw", "spot", "spot_id", "Cell", "Barcode"),
    colnames(df)
  )
  if (length(cand) == 0) return(NA_character_)
  cand[1]
}

prepare_phase06_table <- function(df, obj, sec) {
  cell_col <- find_cell_col(df)
  obj_cells <- get_cell_names(obj)
  
  if (!is.na(cell_col)) {
    df2 <- df %>%
      mutate(cell = as.character(.data[[cell_col]]))
    n_overlap <- sum(df2$cell %in% obj_cells, na.rm = TRUE)
    
    if (n_overlap > 0) {
      df2 <- df2 %>%
        filter(cell %in% obj_cells) %>%
        distinct(cell, .keep_all = TRUE)
      return(df2)
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
  req <- c(
    "cell", "section_id", "img_file", "in_tissue", "is_hotspot_semla",
    "pxl_col_in_hires", "pxl_row_in_hires"
  )
  miss <- setdiff(req, colnames(df))
  if (length(miss) > 0) {
    stop("Section ", sec, " phase06 table missing columns: ", paste(miss, collapse = ", "))
  }
  df
}

compute_scaled_expr <- function(mat, genes_present) {
  if (length(genes_present) == 0) return(NULL)
  
  expr_mat <- mat[genes_present, , drop = FALSE]
  q99 <- apply(expr_mat, 1, function(x) {
    as.numeric(quantile(x, probs = 0.99, na.rm = TRUE, names = FALSE))
  })
  q99[!is.finite(q99) | q99 <= 0] <- 1
  
  expr_scaled <- expr_mat
  for (g in rownames(expr_scaled)) {
    expr_scaled[g, ] <- pmin(expr_scaled[g, ], q99[g]) / q99[g]
  }
  
  expr_scaled
}

hotspot_percentile <- function(x) {
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  n <- sum(ok)
  if (n == 0) return(out)
  # Average ranks make tied zero values sit together rather than being treated
  # as highly enriched simply because the reference distribution is sparse.
  out[ok] <- (rank(x[ok], ties.method = "average") - 0.5) / n
  out
}

compute_signature_detection <- function(expr_mat, state_sets_used) {
  # Number of detected genes from each core signature in every spot.
  # Detection is equivalent on log-normalised data because zero counts remain 0.
  key_map <- c(
    state_cyto = "cyto",
    state_inflam = "inflam",
    state_exhaust = "exhaust",
    state_il17like = "il17like",
    state_retention = "retention"
  )
  
  out <- tibble(cell = colnames(expr_mat))
  for (sc in names(key_map)) {
    genes <- state_sets_used[[key_map[[sc]]]]
    nd_col <- paste0(sc, "_n_detected")
    if (length(genes) == 0) {
      out[[nd_col]] <- 0L
    } else {
      out[[nd_col]] <- as.integer(Matrix::colSums(expr_mat[genes, , drop = FALSE] > 0))
    }
  }
  out
}

call_hotspot_states <- function(sec_tbl,
                                pct_cut = STATE_HOTSPOT_PERCENTILE,
                                margin_cut = STATE_MIN_MARGIN,
                                min_detect = MIN_SIGNATURE_GENES_DETECTED) {
  hs_idx <- which(sec_tbl$in_tissue == 1 & sec_tbl$is_hotspot_semla)
  
  # Initialise output columns for all spots. Calls are intentionally restricted
  # to Semla-defined in-tissue hotspots.
  pct_cols <- paste0(state_score_cols_core, "_hpct")
  for (pc in pct_cols) sec_tbl[[pc]] <- NA_real_
  
  sec_tbl$dominant_state <- NA_character_
  sec_tbl$dominant_state_pretty <- NA_character_
  sec_tbl$dominant_state_percentile <- NA_real_
  sec_tbl$dominant_state_margin <- NA_real_
  sec_tbl$dominant_state_n_detected <- NA_integer_
  sec_tbl$state_call_reason <- ifelse(
    sec_tbl$in_tissue == 1 & sec_tbl$is_hotspot_semla,
    "not_evaluated", "not_hotspot"
  )
  
  if (length(hs_idx) == 0) return(sec_tbl)
  
  # Signature-specific percentile calibration is performed only among hotspots
  # from the same section. This removes between-signature scale differences and
  # section-wide shifts while preserving relative within-section enrichment.
  for (sc in state_score_cols_core) {
    pc <- paste0(sc, "_hpct")
    sec_tbl[[pc]][hs_idx] <- hotspot_percentile(sec_tbl[[sc]][hs_idx])
  }
  
  P <- as.matrix(sec_tbl[hs_idx, pct_cols, drop = FALSE])
  D <- as.matrix(sec_tbl[hs_idx, paste0(state_score_cols_core, "_n_detected"), drop = FALSE])
  eligible <- is.finite(P) & P >= pct_cut & D >= min_detect
  
  for (ii in seq_along(hs_idx)) {
    row_i <- hs_idx[ii]
    elig_i <- which(eligible[ii, ])
    
    if (length(elig_i) == 0) {
      sec_tbl$dominant_state[row_i] <- "state_unassigned"
      sec_tbl$dominant_state_pretty[row_i] <- "Unassigned"
      sec_tbl$state_call_reason[row_i] <- "no_state_passed_thresholds"
      next
    }
    
    # Rank only states that met both the score and gene-detection criteria.
    ord <- elig_i[order(P[ii, elig_i], decreasing = TRUE)]
    top <- ord[1]
    top_pct <- P[ii, top]
    second_pct <- if (length(ord) >= 2) P[ii, ord[2]] else NA_real_
    margin <- if (length(ord) >= 2) top_pct - second_pct else Inf
    
    sec_tbl$dominant_state_percentile[row_i] <- top_pct
    sec_tbl$dominant_state_margin[row_i] <- ifelse(is.finite(margin), margin, NA_real_)
    sec_tbl$dominant_state_n_detected[row_i] <- D[ii, top]
    
    if (length(ord) >= 2 && margin < margin_cut) {
      sec_tbl$dominant_state[row_i] <- "state_unassigned"
      sec_tbl$dominant_state_pretty[row_i] <- "Unassigned"
      sec_tbl$state_call_reason[row_i] <- "ambiguous_top_states"
    } else {
      sc <- state_score_cols_core[top]
      sec_tbl$dominant_state[row_i] <- sc
      sec_tbl$dominant_state_pretty[row_i] <- unname(pretty_state[sc])
      sec_tbl$state_call_reason[row_i] <- "assigned"
    }
  }
  
  sec_tbl
}

summarise_call_sensitivity <- function(sec_tbl,
                                       pct_grid = SENSITIVITY_PERCENTILES,
                                       margin_grid = SENSITIVITY_MARGINS,
                                       min_detect = MIN_SIGNATURE_GENES_DETECTED) {
  hs <- sec_tbl %>% filter(in_tissue == 1, is_hotspot_semla)
  if (nrow(hs) == 0) return(tibble())
  
  pct_cols <- paste0(state_score_cols_core, "_hpct")
  det_cols <- paste0(state_score_cols_core, "_n_detected")
  P <- as.matrix(hs[, pct_cols, drop = FALSE])
  D <- as.matrix(hs[, det_cols, drop = FALSE])
  
  bind_rows(lapply(pct_grid, function(pc) {
    bind_rows(lapply(margin_grid, function(mc) {
      eligible <- is.finite(P) & P >= pc & D >= min_detect
      calls <- rep("Unassigned", nrow(hs))
      reasons <- rep("no_state_passed_thresholds", nrow(hs))
      
      for (i in seq_len(nrow(hs))) {
        ei <- which(eligible[i, ])
        if (length(ei) == 0) next
        ord <- ei[order(P[i, ei], decreasing = TRUE)]
        if (length(ord) >= 2 && (P[i, ord[1]] - P[i, ord[2]]) < mc) {
          reasons[i] <- "ambiguous_top_states"
          next
        }
        calls[i] <- unname(pretty_state[state_score_cols_core[ord[1]]])
        reasons[i] <- "assigned"
      }
      
      tibble(
        section_id = unique(as.character(hs$section_id))[1],
        hotspot_percentile_cutoff = pc,
        top2_margin_cutoff = mc,
        min_signature_genes_detected = min_detect,
        n_hotspots = nrow(hs),
        n_assigned = sum(calls != "Unassigned"),
        n_unassigned = sum(calls == "Unassigned"),
        n_ambiguous = sum(reasons == "ambiguous_top_states"),
        pct_assigned = 100 * mean(calls != "Unassigned"),
        pct_unassigned = 100 * mean(calls == "Unassigned"),
        Cytotoxic = sum(calls == "Cytotoxic"),
        Inflammatory = sum(calls == "Inflammatory"),
        Exhausted = sum(calls == "Exhausted"),
        IL17_like = sum(calls == "IL17-like"),
        Retention = sum(calls == "Retention")
      )
    }))
  }))
}

# -------------------------
# Main per-section loop
# -------------------------
section_state_summaries <- list()
section_tables <- list()
section_retention <- list()
signature_gene_summary <- list()
state_detection_summary <- list()
state_sensitivity_summary <- list()

for (sec in good_sections) {
  message("\n==============================")
  message("Processing section ", sec)
  message("==============================")
  
  vis_sec  <- safe_read_vis(sec)
  phase06  <- safe_read_phase06(sec)
  
  expr_res <- read_expr_matrix(vis_sec)
  vis_sec  <- expr_res$obj
  expr_mat <- expr_res$mat
  
  phase06 <- prepare_phase06_table(phase06, vis_sec, sec)
  phase06 <- ensure_required_phase06_cols(phase06, sec)
  
  present_genes <- rownames(expr_mat)
  
  state_sets_used <- lapply(state_sets_requested, match_requested_genes, present = present_genes)
  retention_genes_present <- match_requested_genes(retention_genes_requested, present_genes)
  
  signature_gene_summary[[sec]] <- tibble(
    section_id = sec,
    signature = c(names(state_sets_used), "retention_panel"),
    n_genes_kept = c(sapply(state_sets_used, length), length(retention_genes_present)),
    genes_kept = c(
      vapply(state_sets_used, function(x) paste(x, collapse = ","), character(1)),
      paste(retention_genes_present, collapse = ",")
    )
  )
  
  vis_sec <- AddModuleScore_UCell(
    vis_sec,
    features = state_sets_used,
    name = "gdstates",
    ncores = 1
  )
  
  ucell_cols <- grep("gdstates$", colnames(vis_sec@meta.data), value = TRUE)
  
  col_cyto      <- grep("^cyto",      ucell_cols, value = TRUE)[1]
  col_inflam    <- grep("^inflam",    ucell_cols, value = TRUE)[1]
  col_exhaust   <- grep("^exhaust",   ucell_cols, value = TRUE)[1]
  col_il17      <- grep("^il17like",  ucell_cols, value = TRUE)[1]
  col_prolif    <- grep("^prolif",    ucell_cols, value = TRUE)[1]
  col_retention <- grep("^retention", ucell_cols, value = TRUE)[1]
  
  md_state <- vis_sec@meta.data %>%
    rownames_to_column("cell") %>%
    transmute(
      cell,
      state_cyto      = .data[[col_cyto]],
      state_inflam    = .data[[col_inflam]],
      state_exhaust   = .data[[col_exhaust]],
      state_il17like  = .data[[col_il17]],
      state_prolif    = .data[[col_prolif]],
      state_retention = .data[[col_retention]]
    )
  
  sig_detect_df <- compute_signature_detection(expr_mat, state_sets_used)
  
  expr_scaled <- compute_scaled_expr(expr_mat, retention_genes_present)
  if (!is.null(expr_scaled)) {
    scaled_df <- as.data.frame(t(expr_scaled))
    colnames(scaled_df) <- paste0("scaled_", colnames(t(expr_scaled)))
    scaled_df <- scaled_df %>%
      rownames_to_column("cell")
  } else {
    scaled_df <- tibble(cell = get_cell_names(vis_sec))
  }
  
  sec_tbl <- phase06 %>%
    left_join(md_state, by = "cell") %>%
    left_join(sig_detect_df, by = "cell") %>%
    left_join(scaled_df, by = "cell") %>%
    mutate(
      section_id = as.character(section_id),
      in_tissue = as.integer(in_tissue),
      is_hotspot_semla = ifelse(is.na(is_hotspot_semla), FALSE, as.logical(is_hotspot_semla))
    )
  
  # High-confidence dominant-state calls among hotspots only.
  sec_tbl <- call_hotspot_states(sec_tbl)
  
  # Signature-gene detection QC in the hotspot compartment.
  hs_tmp <- sec_tbl %>% filter(in_tissue == 1, is_hotspot_semla)
  state_detection_summary[[sec]] <- bind_rows(lapply(state_score_cols_core, function(sc) {
    nd <- hs_tmp[[paste0(sc, "_n_detected")]]
    tibble(
      section_id = sec,
      signature = unname(pretty_state[sc]),
      score_column = sc,
      n_hotspots = length(nd),
      median_genes_detected = median(nd, na.rm = TRUE),
      q25_genes_detected = as.numeric(quantile(nd, 0.25, na.rm = TRUE, names = FALSE)),
      q75_genes_detected = as.numeric(quantile(nd, 0.75, na.rm = TRUE, names = FALSE)),
      pct_with_0_genes = 100 * mean(nd == 0, na.rm = TRUE),
      pct_with_1_gene = 100 * mean(nd == 1, na.rm = TRUE),
      pct_with_2plus_genes = 100 * mean(nd >= 2, na.rm = TRUE)
    )
  }))
  
  # Sensitivity analysis across nearby percentile and top-two margin cutoffs.
  state_sensitivity_summary[[sec]] <- summarise_call_sensitivity(sec_tbl)
  
  front_cols <- c(
    "cell", "section_id", "img_file",
    "in_tissue", "is_hotspot_semla",
    "zone_call",
    "pxl_col_in_hires", "pxl_row_in_hires",
    "Lymphocytes", "local_lymph_score",
    "frac_epithelial_nbr", "frac_myeloid_nbr", "frac_mesenchymal_nbr",
    "state_cyto", "state_inflam", "state_exhaust", "state_il17like", "state_retention", "state_prolif",
    "state_cyto_hpct", "state_inflam_hpct", "state_exhaust_hpct", "state_il17like_hpct", "state_retention_hpct",
    "state_cyto_n_detected", "state_inflam_n_detected", "state_exhaust_n_detected",
    "state_il17like_n_detected", "state_retention_n_detected",
    "dominant_state", "dominant_state_pretty", "dominant_state_percentile",
    "dominant_state_margin", "dominant_state_n_detected", "state_call_reason"
  )
  keep_front <- intersect(front_cols, colnames(sec_tbl))
  other_cols <- setdiff(colnames(sec_tbl), keep_front)
  sec_tbl <- sec_tbl[, c(keep_front, other_cols), drop = FALSE]
  
  section_tables[[sec]] <- sec_tbl
  
  sec_ret <- sec_tbl %>%
    select(
      any_of(c(
        "cell", "section_id", "img_file", "in_tissue", "is_hotspot_semla",
        "pxl_col_in_hires", "pxl_row_in_hires"
      )),
      starts_with("scaled_")
    )
  section_retention[[sec]] <- sec_ret
  
  state_summary <- sec_tbl %>%
    filter(in_tissue == 1, is_hotspot_semla) %>%
    count(dominant_state_pretty, name = "n_hotspots_state") %>%
    complete(
      dominant_state_pretty = factor(state_levels, levels = state_levels),
      fill = list(n_hotspots_state = 0)
    ) %>%
    mutate(
      section_id = sec,
      n_hotspots_total = sum(n_hotspots_state),
      pct_hotspots = ifelse(n_hotspots_total > 0, 100 * n_hotspots_state / n_hotspots_total, 0)
    ) %>%
    select(section_id, dominant_state_pretty, n_hotspots_state, n_hotspots_total, pct_hotspots)
  
  section_state_summaries[[sec]] <- state_summary
  
  write_tsv(sec_tbl, file.path(out_sections, paste0("phase07_section_", sec, "_plotready.tsv")))
  saveRDS(sec_tbl, file.path(out_sections, paste0("phase07_section_", sec, "_plotready.rds")))
  
  write_tsv(sec_ret, file.path(out_sections, paste0("phase07_section_", sec, "_retention_scaled.tsv")))
  saveRDS(sec_ret, file.path(out_sections, paste0("phase07_section_", sec, "_retention_scaled.rds")))
  
  qc_one <- tibble(
    section_id = sec,
    n_spots_total = nrow(sec_tbl),
    n_tissue = sum(sec_tbl$in_tissue == 1, na.rm = TRUE),
    n_hotspots = sum(sec_tbl$is_hotspot_semla, na.rm = TRUE),
    n_hotspots_in_tissue = sum(sec_tbl$is_hotspot_semla & sec_tbl$in_tissue == 1, na.rm = TRUE),
    n_background_nonhotspot = sum(sec_tbl$in_tissue == 1 & !sec_tbl$is_hotspot_semla, na.rm = TRUE),
    n_hotspots_assigned = sum(sec_tbl$in_tissue == 1 & sec_tbl$is_hotspot_semla &
                                sec_tbl$dominant_state_pretty != "Unassigned", na.rm = TRUE),
    n_hotspots_unassigned = sum(sec_tbl$in_tissue == 1 & sec_tbl$is_hotspot_semla &
                                  sec_tbl$dominant_state_pretty == "Unassigned", na.rm = TRUE),
    pct_hotspots_unassigned = ifelse(
      sum(sec_tbl$in_tissue == 1 & sec_tbl$is_hotspot_semla, na.rm = TRUE) > 0,
      100 * sum(sec_tbl$in_tissue == 1 & sec_tbl$is_hotspot_semla &
                  sec_tbl$dominant_state_pretty == "Unassigned", na.rm = TRUE) /
        sum(sec_tbl$in_tissue == 1 & sec_tbl$is_hotspot_semla, na.rm = TRUE),
      NA_real_
    ),
    state_hotspot_percentile = STATE_HOTSPOT_PERCENTILE,
    state_min_margin = STATE_MIN_MARGIN,
    min_signature_genes_detected = MIN_SIGNATURE_GENES_DETECTED,
    n_hotspots_ambiguous = sum(sec_tbl$in_tissue == 1 & sec_tbl$is_hotspot_semla &
                                 sec_tbl$state_call_reason == "ambiguous_top_states", na.rm = TRUE),
    img_file = unique(sec_tbl$img_file)[1]
  )
  write_tsv(qc_one, file.path(out_qc, paste0("phase07_qc_section_", sec, ".tsv")))
}

# -------------------------
# Global outputs
# -------------------------
all_sections_tbl <- bind_rows(section_tables)
all_hotspots_tbl <- all_sections_tbl %>%
  filter(in_tissue == 1, is_hotspot_semla)

freq_df <- bind_rows(section_state_summaries) %>%
  mutate(
    section_id = factor(section_id, levels = good_sections),
    dominant_state_pretty = factor(dominant_state_pretty, levels = state_levels)
  ) %>%
  arrange(section_id, dominant_state_pretty)

qc_freq <- freq_df %>%
  group_by(section_id) %>%
  summarise(
    n_hotspots_total = unique(n_hotspots_total),
    pct_sum = sum(pct_hotspots),
    .groups = "drop"
  )

sig_gene_df <- bind_rows(signature_gene_summary) %>%
  arrange(section_id, signature)

state_detection_df <- bind_rows(state_detection_summary) %>%
  arrange(section_id, signature)

state_sensitivity_df <- bind_rows(state_sensitivity_summary) %>%
  arrange(section_id, hotspot_percentile_cutoff, top2_margin_cutoff)

write_tsv(all_sections_tbl, file.path(out_tables, "phase07_hotspot_states_all_sections.tsv"))
write_tsv(all_hotspots_tbl, file.path(out_tables, "phase07_hotspot_states_hotspots_only.tsv"))
write_tsv(freq_df, file.path(out_tables, "phase07_state_frequency_by_section.tsv"))
write_tsv(qc_freq, file.path(out_tables, "phase07_state_frequency_qc.tsv"))
write_tsv(sig_gene_df, file.path(out_tables, "phase07_signature_gene_summary.tsv"))
write_tsv(state_detection_df, file.path(out_tables, "phase07_state_signature_detection_qc.tsv"))
write_tsv(state_sensitivity_df, file.path(out_tables, "phase07_state_calling_sensitivity.tsv"))

write_tsv(
  tibble(
    parameter = c("state_hotspot_percentile", "state_min_margin",
                  "minimum_signature_genes_detected"),
    value = c(STATE_HOTSPOT_PERCENTILE, STATE_MIN_MARGIN,
              MIN_SIGNATURE_GENES_DETECTED)
  ),
  file.path(out_tables, "phase07_state_calling_parameters.tsv")
)

bundle <- list(
  good_sections = good_sections,
  panelA_sections = panelA_sections,
  panelC_section = panelC_section,
  state_sets_requested = state_sets_requested,
  retention_genes_requested = retention_genes_requested,
  state_levels = state_levels,
  state_hotspot_percentile = STATE_HOTSPOT_PERCENTILE,
  state_min_margin = STATE_MIN_MARGIN,
  min_signature_genes_detected = MIN_SIGNATURE_GENES_DETECTED,
  all_sections = all_sections_tbl,
  hotspots_only = all_hotspots_tbl,
  frequency = freq_df,
  frequency_qc = qc_freq,
  signature_gene_summary = sig_gene_df,
  state_detection_qc = state_detection_df,
  state_calling_sensitivity = state_sensitivity_df,
  section_tables = section_tables,
  section_retention = section_retention
)
saveRDS(bundle, file.path(out_rds, "phase07_analysis_bundle.rds"))

message("\nPhase 07 analysis complete.")
message("Main outputs written to: ", out_root)
