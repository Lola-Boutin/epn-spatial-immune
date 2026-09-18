# =============================================================================
# 03c_compare_original_vs_tumor_augmented_hotspots.R
#
# Purpose
# -------
# Quantify whether the original lymphocyte hotspot pattern survives after
# adding representative real PFA1/PFA2 tumor cells to the Semla NNLS reference.
#
# Reads only; modifies no original or augmented NNLS files.
#
# The augmented hotspot call is rebuilt using the SAME rule as 02c:
#   within-section 90th percentile of the raw Lymphocytes score, AND
#   within-section 90th percentile of the 6-nearest-neighbour mean.
# Because both thresholds are within-section, hotspots exist in both arms by
# construction. What this script tests is whether the SAME SPOTS are called,
# not whether hotspots exist.
#
# Inputs
# ------
#   stage_dir("semla_nnls")       good_sections/nnls,
#                                 good_sections_hotspots_tissue/nnls
#   stage_dir("semla_nnls_aug")   all_sections/nnls (preferred) or good_sections/nnls
#
# Outputs
# -------
#   stage_dir("hotspot_robustness")
#     score_agreement.tsv, hotspot_overlap.tsv, hotspot_thresholds.tsv,
#     augmented_score_in_original_hotspots.tsv, per_spot_comparison.tsv,
#     Hotspot_tumor_augmented_PFA_robustness.{png,pdf}, comparison_summary.txt
#
# Stochastic
# ----------
#   none
#
# Runtime
# -------
#   ~10 minutes
# =============================================================================

source(here::here("config.R"))

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})


# =============================================================================
# 0. PATHS / SETTINGS
# =============================================================================

original_root  <- stage_dir("semla_nnls")
augmented_root <- stage_dir("semla_nnls_aug")

original_nnls <- file.path(original_root, "good_sections", "nnls")
original_hotspot_nnls <- file.path(stage_dir("hotspots"), "nnls")

# For the sensitivity comparison, use the exact same six sections whether or not
# the augmented run's own QC would label them "good". Prefer all_sections.
augmented_nnls_candidates <- c(
  file.path(augmented_root, "all_sections", "nnls"),
  file.path(augmented_root, "good_sections", "nnls")
)
augmented_nnls <- augmented_nnls_candidates[
  dir.exists(augmented_nnls_candidates)
][1]

if (length(augmented_nnls) == 0 || is.na(augmented_nnls)) {
  stop(
    "Could not find PFA tumor-augmented NNLS output under:\n  ",
    augmented_root,
    "\nExpected all_sections/nnls or good_sections/nnls."
  )
}

out_dir <- stage_dir("hotspot_robustness")

good_sections <- GOOD_SECTIONS

# MUST match the original hotspot definition in 02c
K_NEIGHBORS    <- as.integer(THRESH$hotspot_knn)
EXPR_QUANTILE  <- THRESH$hotspot_quantile
LOCAL_QUANTILE <- THRESH$hotspot_quantile

# =============================================================================
# 1. HELPERS
# =============================================================================

read_one <- function(dir, sec, label) {
  
  f_rds <- file.path(dir, paste0("NNLS_section_", sec, ".rds"))
  f_tsv <- file.path(dir, paste0("NNLS_section_", sec, ".tsv"))
  
  if (file.exists(f_rds)) {
    df <- readRDS(f_rds)
  } else if (file.exists(f_tsv)) {
    df <- readr::read_tsv(f_tsv, show_col_types = FALSE)
  } else {
    stop(
      "Missing NNLS table for section ", sec, " in:\n  ", dir,
      "\nExpected .rds or .tsv."
    )
  }
  
  if (!"cell" %in% names(df) && "cell.x" %in% names(df)) {
    df <- dplyr::rename(df, cell = cell.x)
  }
  
  if (!"cell" %in% names(df)) {
    stop("No 'cell' column in ", label, " section ", sec)
  }
  
  tibble::as_tibble(df)
}


safe_cor <- function(x, y, method = "spearman") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 20L) return(NA_real_)
  if (sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = method))
}


local_knn_mean <- function(xy, values, k = 6L) {
  if (!requireNamespace("dbscan", quietly = TRUE)) {
    stop("Package 'dbscan' is required.")
  }
  
  n <- nrow(xy)
  k_use <- min(k, n - 1L)
  
  if (k_use < 1L) return(values)
  
  kn <- dbscan::kNN(xy, k = k_use)
  
  vapply(
    seq_len(n),
    function(i) mean(values[kn$id[i, ]], na.rm = TRUE),
    numeric(1)
  )
}


jaccard_cells <- function(a, b) {
  a <- unique(as.character(a))
  b <- unique(as.character(b))
  u <- union(a, b)
  if (length(u) == 0L) return(NA_real_)
  length(intersect(a, b)) / length(u)
}


dice_cells <- function(a, b) {
  a <- unique(as.character(a))
  b <- unique(as.character(b))
  den <- length(a) + length(b)
  if (den == 0L) return(NA_real_)
  2 * length(intersect(a, b)) / den
}


# =============================================================================
# 2. LOAD ORIGINAL + AUGMENTED PER-SPOT SCORES
# =============================================================================

message("Original NNLS:       ", original_nnls)
message("Original hotspots:   ", original_hotspot_nnls)
message("PFA tumor-augmented NNLS: ", augmented_nnls)

dat <- lapply(good_sections, function(sec) {
  
  hs  <- read_one(original_hotspot_nnls, sec, "original hotspot")
  org <- read_one(original_nnls, sec, "original NNLS")
  aug <- read_one(augmented_nnls, sec, "PFA tumor-augmented NNLS")
  
  required_hs <- c(
    "cell", "in_tissue", "pxl_col_in_hires", "pxl_row_in_hires",
    "Lymphocytes", "local_lymph_score", "is_hotspot_semla"
  )
  
  missing_hs <- setdiff(required_hs, names(hs))
  if (length(missing_hs) > 0L) {
    stop(
      "Original hotspot table for section ", sec,
      " is missing: ", paste(missing_hs, collapse = ", ")
    )
  }
  
  if (!"Lymphocytes" %in% names(org)) {
    stop("Original NNLS table has no Lymphocytes column for section ", sec)
  }
  
  if (!"Lymphocytes" %in% names(aug)) {
    stop("Augmented NNLS table has no Lymphocytes column for section ", sec)
  }
  
  base <- hs %>%
    transmute(
      cell,
      section = as.character(sec),
      in_tissue = as.integer(in_tissue),
      x = as.numeric(pxl_col_in_hires),
      y = as.numeric(pxl_row_in_hires),
      Ly_original_from_hotspot_table = as.numeric(Lymphocytes),
      local_original = as.numeric(local_lymph_score),
      hotspot_original = as.logical(is_hotspot_semla) %in% TRUE
    )
  
  joined <- base %>%
    left_join(
      org %>% transmute(cell, Ly_original = as.numeric(Lymphocytes)),
      by = "cell"
    ) %>%
    left_join(
      aug %>% transmute(cell, Ly_augmented = as.numeric(Lymphocytes)),
      by = "cell"
    )
  
  n_missing_aug <- sum(!is.finite(joined$Ly_augmented))
  if (n_missing_aug > 0L) {
    warning(
      "Section ", sec, ": ", n_missing_aug,
      " original in-tissue spots lack an augmented score."
    )
  }
  
  joined
}) %>%
  bind_rows() %>%
  filter(
    in_tissue == 1L,
    is.finite(x),
    is.finite(y)
  ) %>%
  mutate(section = factor(section, levels = good_sections))


# Sanity check: original main score should equal original hotspot-table score.
sanity <- dat %>%
  filter(
    is.finite(Ly_original_from_hotspot_table),
    is.finite(Ly_original)
  ) %>%
  summarise(
    max_abs_diff = max(
      abs(Ly_original_from_hotspot_table - Ly_original),
      na.rm = TRUE
    )
  )

message(
  "Sanity: max |original hotspot-table score - original main score| = ",
  signif(sanity$max_abs_diff, 4),
  " (should be ~0)"
)

cat("\nSpots loaded per section:\n")
dat %>%
  count(section, name = "n_spots") %>%
  as.data.frame() %>%
  print(row.names = FALSE)


# =============================================================================
# 3. SCORE AGREEMENT
# =============================================================================

score_agreement <- dat %>%
  group_by(section) %>%
  summarise(
    n_spots = n(),
    
    rho_raw_original_augmented =
      safe_cor(Ly_original, Ly_augmented, method = "spearman"),
    
    pct_nonzero_original =
      100 * mean(Ly_original > 0, na.rm = TRUE),
    
    pct_nonzero_augmented =
      100 * mean(Ly_augmented > 0, na.rm = TRUE),
    
    mean_original =
      mean(Ly_original, na.rm = TRUE),
    
    mean_augmented =
      mean(Ly_augmented, na.rm = TRUE),
    
    median_original =
      median(Ly_original, na.rm = TRUE),
    
    median_augmented =
      median(Ly_augmented, na.rm = TRUE),
    
    .groups = "drop"
  )

cat("\n=== SCORE AGREEMENT ===\n")
print(as.data.frame(score_agreement), row.names = FALSE)

write_tsv(
  score_agreement,
  file.path(out_dir, "score_agreement.tsv")
)


# =============================================================================
# 4. REBUILD HOTSPOTS ON AUGMENTED SCORE
#    SAME 90th raw + 90th 6-NN rule as original
# =============================================================================

redefined_list <- list()
threshold_list <- list()

for (sec in good_sections) {
  
  d <- dat %>%
    filter(section == sec) %>%
    as_tibble()
  
  ok <- is.finite(d$Ly_augmented) & is.finite(d$x) & is.finite(d$y)
  
  d$local_augmented <- NA_real_
  d$hotspot_augmented <- FALSE
  
  if (sum(ok) < 10L) {
    warning("Section ", sec, ": fewer than 10 valid augmented spots.")
    redefined_list[[sec]] <- d
    threshold_list[[sec]] <- tibble(
      section = sec,
      raw_threshold_original = NA_real_,
      local_threshold_original = NA_real_,
      raw_threshold_augmented = NA_real_,
      local_threshold_augmented = NA_real_
    )
    next
  }
  
  xy <- as.matrix(d[ok, c("x", "y")])
  
  d$local_augmented[ok] <- local_knn_mean(
    xy,
    d$Ly_augmented[ok],
    k = K_NEIGHBORS
  )
  
  raw_thr_aug <- quantile(
    d$Ly_augmented[ok],
    EXPR_QUANTILE,
    na.rm = TRUE,
    names = FALSE
  )
  
  local_thr_aug <- quantile(
    d$local_augmented[ok],
    LOCAL_QUANTILE,
    na.rm = TRUE,
    names = FALSE
  )
  
  d$hotspot_augmented[ok] <-
    d$Ly_augmented[ok] >= raw_thr_aug &
    d$local_augmented[ok] >= local_thr_aug
  
  # Reconstruct thresholds from original columns for documentation.
  ok_orig <- is.finite(d$Ly_original) & is.finite(d$local_original)
  
  raw_thr_orig <- if (sum(ok_orig) > 0L) {
    quantile(
      d$Ly_original[ok_orig],
      EXPR_QUANTILE,
      na.rm = TRUE,
      names = FALSE
    )
  } else NA_real_
  
  local_thr_orig <- if (sum(ok_orig) > 0L) {
    quantile(
      d$local_original[ok_orig],
      LOCAL_QUANTILE,
      na.rm = TRUE,
      names = FALSE
    )
  } else NA_real_
  
  threshold_list[[sec]] <- tibble(
    section = sec,
    raw_threshold_original = raw_thr_orig,
    local_threshold_original = local_thr_orig,
    raw_threshold_augmented = raw_thr_aug,
    local_threshold_augmented = local_thr_aug
  )
  
  redefined_list[[sec]] <- d
}

redefined <- bind_rows(redefined_list) %>%
  mutate(section = factor(section, levels = good_sections))

thresholds <- bind_rows(threshold_list)

# Correlation of spatially smoothed scores too.
local_agreement <- redefined %>%
  group_by(section) %>%
  summarise(
    rho_local_original_augmented =
      safe_cor(local_original, local_augmented, method = "spearman"),
    .groups = "drop"
  )

score_agreement <- score_agreement %>%
  left_join(local_agreement, by = "section")

# Rewrite score_agreement with local rho included.
write_tsv(
  score_agreement,
  file.path(out_dir, "score_agreement.tsv")
)

write_tsv(
  thresholds,
  file.path(out_dir, "hotspot_thresholds.tsv")
)


# =============================================================================
# 5. HOTSPOT OVERLAP
# =============================================================================

hotspot_overlap <- redefined %>%
  group_by(section) %>%
  summarise(
    n_hot_original = sum(hotspot_original, na.rm = TRUE),
    n_hot_augmented = sum(hotspot_augmented, na.rm = TRUE),
    n_shared = sum(
      hotspot_original & hotspot_augmented,
      na.rm = TRUE
    ),
    
    jaccard = jaccard_cells(
      cell[hotspot_original],
      cell[hotspot_augmented]
    ),
    
    dice = dice_cells(
      cell[hotspot_original],
      cell[hotspot_augmented]
    ),
    
    # Sensitivity / recovery of the original hotspot set
    pct_original_recovered =
      100 * sum(hotspot_original & hotspot_augmented, na.rm = TRUE) /
      pmax(sum(hotspot_original, na.rm = TRUE), 1),
    
    # Fraction of augmented calls that were already original hotspots
    pct_augmented_already_original =
      100 * sum(hotspot_original & hotspot_augmented, na.rm = TRUE) /
      pmax(sum(hotspot_augmented, na.rm = TRUE), 1),
    
    .groups = "drop"
  )

cat("\n=== HOTSPOT OVERLAP ===\n")
print(as.data.frame(hotspot_overlap), row.names = FALSE)

write_tsv(
  hotspot_overlap,
  file.path(out_dir, "hotspot_overlap.tsv")
)


# =============================================================================
# 6. AUGMENTED SCORE INSIDE ORIGINAL HOTSPOTS
#
# Descriptive only. We deliberately do NOT perform a spot-level Wilcoxon test:
# spots are spatially autocorrelated and are not independent biological samples.
# =============================================================================

score_shift <- redefined %>%
  group_by(section) %>%
  summarise(
    n_hot_original = sum(hotspot_original, na.rm = TRUE),
    
    augmented_mean_original_hotspot =
      mean(Ly_augmented[hotspot_original], na.rm = TRUE),
    
    augmented_mean_original_nonhotspot =
      mean(Ly_augmented[!hotspot_original], na.rm = TRUE),
    
    augmented_median_original_hotspot =
      median(Ly_augmented[hotspot_original], na.rm = TRUE),
    
    augmented_median_original_nonhotspot =
      median(Ly_augmented[!hotspot_original], na.rm = TRUE),
    
    mean_ratio_hot_vs_nonhot =
      augmented_mean_original_hotspot /
      pmax(augmented_mean_original_nonhotspot, 1e-12),
    
    mean_difference_hot_minus_nonhot =
      augmented_mean_original_hotspot -
      augmented_mean_original_nonhotspot,
    
    .groups = "drop"
  )

cat("\n=== AUGMENTED SCORE IN ORIGINAL HOTSPOTS ===\n")
print(as.data.frame(score_shift), row.names = FALSE)

write_tsv(
  score_shift,
  file.path(out_dir, "augmented_score_in_original_hotspots.tsv")
)


# =============================================================================
# 7. SAVE PER-SPOT TABLE
# =============================================================================

per_spot <- redefined %>%
  mutate(
    overlap_class = case_when(
      hotspot_original & hotspot_augmented ~ "Both",
      hotspot_original & !hotspot_augmented ~ "Original only",
      !hotspot_original & hotspot_augmented ~ "PFA tumor-augmented only",
      TRUE ~ "Neither"
    )
  ) %>%
  select(
    cell, section, x, y,
    Ly_original, Ly_augmented,
    local_original, local_augmented,
    hotspot_original, hotspot_augmented,
    overlap_class
  )

write_tsv(
  per_spot,
  file.path(out_dir, "per_spot_comparison.tsv")
)


# =============================================================================
# 8. FIGURE
# =============================================================================

# A. Score rank agreement. Do NOT draw y=x because coefficient scales need not
# be identical after changing the reference.
p_score <- ggplot(
  redefined,
  aes(x = Ly_original, y = Ly_augmented)
) +
  geom_point(size = 0.25, alpha = 0.25, colour = "grey30") +
  facet_wrap(~section, nrow = 2, scales = "free") +
  labs(
    x = "Original Lymphocyte NNLS score",
    y = "PFA tumor-augmented Lymphocyte NNLS score",
    title = "Per-spot score agreement"
  ) +
  theme_bw(base_size = 9) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey93", colour = NA)
  )


# B. Spatial overlap maps
map_cols <- c(
  "Both" = "#4DAF4A",
  "Original only" = "#E41A1C",
  "PFA tumor-augmented only" = "#377EB8"
)

map_dat <- per_spot %>%
  mutate(
    overlap_class = factor(
      overlap_class,
      levels = c("Both", "Original only", "PFA tumor-augmented only", "Neither")
    )
  )

make_map <- function(sec) {
  
  d <- map_dat %>% filter(section == sec)
  dh <- d %>% filter(overlap_class != "Neither")
  
  ggplot() +
    geom_point(
      data = d,
      aes(x, -y),
      colour = "grey88",
      size = 0.25
    ) +
    geom_point(
      data = dh,
      aes(x, -y, colour = overlap_class),
      size = 0.75
    ) +
    scale_colour_manual(
      values = map_cols,
      name = NULL,
      drop = FALSE
    ) +
    coord_equal() +
    ggtitle(sec) +
    theme_void(base_size = 9) +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5,
        size = 9
      )
    )
}

p_maps <- wrap_plots(
  lapply(good_sections, make_map),
  nrow = 2,
  guides = "collect"
) &
  theme(legend.position = "bottom")


# C. Jaccard + recovery overview
metric_plot <- hotspot_overlap %>%
  transmute(
    section,
    Jaccard = jaccard,
    `Original hotspots recovered` = pct_original_recovered / 100
  ) %>%
  pivot_longer(
    -section,
    names_to = "metric",
    values_to = "value"
  ) %>%
  ggplot(
    aes(x = section, y = value, shape = metric, group = metric)
  ) +
  geom_line(position = position_dodge(width = 0.12), linewidth = 0.5) +
  geom_point(position = position_dodge(width = 0.12), size = 2.4) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2)
  ) +
  labs(
    x = "Section",
    y = "Overlap",
    shape = NULL,
    title = "Hotspot-set robustness"
  ) +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank())


pfig <- (p_score / p_maps / metric_plot) +
  plot_layout(heights = c(1.1, 1.2, 0.9))

ggsave(
  file.path(out_dir, "Hotspot_tumor_augmented_PFA_robustness.png"),
  pfig,
  width = 11,
  height = 11,
  dpi = 600,
  bg = "white"
)

pdf_device <- if (isTRUE(capabilities("cairo"))) cairo_pdf else "pdf"

ggsave(
  file.path(out_dir, "Hotspot_tumor_augmented_PFA_robustness.pdf"),
  pfig,
  width = 11,
  height = 11,
  device = pdf_device,
  bg = "white"
)


# =============================================================================
# 9. SUMMARY
# =============================================================================

summary_lines <- c(
  "Original vs PFA tumor-augmented Semla hotspot robustness",
  "=========================================================",
  "",
  "Comparison:",
  "  Original immune-focused single-cell NNLS reference",
  "  vs.",
  "  Same reference augmented with representative tumor cells sampled only from PFA1/PFA2 tumors",
  "",
  paste0(
    "Hotspot rule: raw >= within-section ",
    EXPR_QUANTILE * 100,
    "th percentile AND 6-NN mean >= within-section ",
    LOCAL_QUANTILE * 100,
    "th percentile"
  ),
  "",
  "Per-section raw-score Spearman:",
  paste0(
    "  ", score_agreement$section, ": ",
    sprintf("%.3f", score_agreement$rho_raw_original_augmented)
  ),
  "",
  "Per-section local-score Spearman:",
  paste0(
    "  ", score_agreement$section, ": ",
    sprintf("%.3f", score_agreement$rho_local_original_augmented)
  ),
  "",
  "Hotspot Jaccard:",
  paste0(
    "  ", hotspot_overlap$section, ": ",
    sprintf("%.3f", hotspot_overlap$jaccard)
  ),
  "",
  "Original hotspot recovery:",
  paste0(
    "  ", hotspot_overlap$section, ": ",
    sprintf("%.1f%%", hotspot_overlap$pct_original_recovered)
  ),
  "",
  "NOTE:",
  "  Similar spatial patterns + high score-rank agreement + substantial hotspot",
  "  recovery support robustness to adding explicit tumor competitors.",
  "  Do not use a single arbitrary Jaccard cutoff as a pass/fail criterion;",
  "  interpret the six sections together with the spatial maps.",
  "",
  paste0("Output directory: ", out_dir)
)

writeLines(
  summary_lines,
  file.path(out_dir, "comparison_summary.txt")
)


cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
cat("No original or PFA tumor-augmented NNLS files were modified.\n\n")

cat("Raw-score Spearman range: ",
    sprintf(
      "%.3f to %.3f\n",
      min(score_agreement$rho_raw_original_augmented, na.rm = TRUE),
      max(score_agreement$rho_raw_original_augmented, na.rm = TRUE)
    ),
    sep = "")

cat("Hotspot Jaccard range: ",
    sprintf(
      "%.3f to %.3f\n",
      min(hotspot_overlap$jaccard, na.rm = TRUE),
      max(hotspot_overlap$jaccard, na.rm = TRUE)
    ),
    sep = "")

cat("Original hotspot recovery range: ",
    sprintf(
      "%.1f%% to %.1f%%\n",
      min(hotspot_overlap$pct_original_recovered, na.rm = TRUE),
      max(hotspot_overlap$pct_original_recovered, na.rm = TRUE)
    ),
    sep = "")

cat("\nOutputs written to:\n  ", out_dir, "\n", sep = "")
cat("============================================================\n")
