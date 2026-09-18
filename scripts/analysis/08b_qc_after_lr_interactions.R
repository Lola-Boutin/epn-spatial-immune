# =============================================================================
# 08b_qc_after_lr_interactions.R
#
# Purpose
# -------
# Cross-section QC of the CellChatDB and CellPhoneDB results from stage 08,
# to detect sections that dominate or distort the aggregated pair rankings.
#
# NOTE: this script lists all six sections while stage 08 itself runs on five
# (723 absent). That is intentional here -- a section with no stage 08 output
# simply yields nothing to read -- but it means the QC set and the analysis set
# differ. Confirm this is what you want before publication.
#
# Inputs
# ------
#   stage_dir("lr_interactions")   observed_scores/, top_pairs/, permutations/
#
# Outputs
# -------
#   stage_dir("lr_qc")   burden, consistency and outlier tables
#
# Stochastic
# ----------
#   none
#
# Runtime
# -------
#   ~2 minutes
# =============================================================================

source(here::here("config.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(purrr)
})


# -------------------------
# Configuration
# -------------------------
phase08_root <- stage_dir("lr_interactions")

out_obs  <- file.path(phase08_root, "observed_scores")
out_top  <- file.path(phase08_root, "top_pairs")
out_perm <- file.path(phase08_root, "permutations")

out_qc <- stage_dir("lr_qc")

resources <- c("CellChatDB", "CellPhoneDB")
zones <- c("Myeloid", "Mesenchymal", "Vascular")

# All six -- stage 08 runs on five (see header).
good_sections <- GOOD_SECTIONS

# -------------------------
# helper
# -------------------------
read_one_observed <- function(resource, zone) {
  f <- file.path(out_obs, paste0(resource, "_", zone, "_observed.tsv"))
  if (!file.exists(f)) return(NULL)
  
  read_tsv(f, show_col_types = FALSE) %>%
    mutate(
      resource = resource,
      zone = zone,
      pair = paste(ligand, receptor, sep = " | ")
    )
}

read_one_top <- function(resource, zone, topN = 30) {
  f <- file.path(out_top, paste0(resource, "_", zone, "_top", topN, ".tsv"))
  if (!file.exists(f)) return(NULL)
  
  read_tsv(f, show_col_types = FALSE) %>%
    mutate(
      resource = resource,
      zone = zone,
      pair = paste(ligand, receptor, sep = " | ")
    )
}

read_one_perm <- function(resource, zone, topN = 30) {
  f <- file.path(out_perm, paste0(resource, "_", zone, "_perm_top", topN, ".tsv"))
  if (!file.exists(f)) return(NULL)
  
  read_tsv(f, show_col_types = FALSE) %>%
    mutate(
      resource = resource,
      zone = zone,
      pair = paste(ligand, receptor, sep = " | ")
    )
}

# -------------------------
# load all outputs
# -------------------------
obs_all <- bind_rows(lapply(resources, function(r) {
  bind_rows(lapply(zones, function(z) read_one_observed(r, z)))
}))

top_all <- bind_rows(lapply(resources, function(r) {
  bind_rows(lapply(zones, function(z) read_one_top(r, z, topN = 30)))
}))

perm_all <- bind_rows(lapply(resources, function(r) {
  bind_rows(lapply(zones, function(z) read_one_perm(r, z, topN = 30)))
}))

if (nrow(obs_all) == 0) stop("No observed LR files found.")
if (nrow(top_all) == 0) warning("No top-pair files found.")
if (nrow(perm_all) == 0) warning("No permutation files found.")

# -------------------------
# 1) broad section-level burden
# -------------------------
section_burden <- obs_all %>%
  group_by(resource, zone, section_id) %>%
  summarise(
    n_pairs = n(),
    mean_score = mean(score, na.rm = TRUE),
    median_score = median(score, na.rm = TRUE),
    max_score = max(score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    z_mean_score = ave(mean_score, resource, zone,
                       FUN = function(x) as.numeric(scale(x))),
    z_n_pairs = ave(n_pairs, resource, zone,
                    FUN = function(x) as.numeric(scale(x)))
  ) %>%
  arrange(resource, zone, desc(mean_score))

write_tsv(section_burden, file.path(out_qc, "section_burden_summary.tsv"))

# flag possible outliers
section_outliers <- section_burden %>%
  mutate(
    outlier_mean = abs(z_mean_score) >= 2,
    outlier_n_pairs = abs(z_n_pairs) >= 2
  ) %>%
  filter(outlier_mean | outlier_n_pairs)

write_tsv(section_outliers, file.path(out_qc, "section_burden_outliers.tsv"))

# -------------------------
# 2) top-pair consistency by section
# -------------------------
if (nrow(top_all) > 0) {
  top_pair_section <- obs_all %>%
    semi_join(top_all %>% select(resource, zone, ligand, receptor), 
              by = c("resource", "zone", "ligand", "receptor")) %>%
    group_by(resource, zone, pair, section_id) %>%
    summarise(score = mean(score, na.rm = TRUE), .groups = "drop")
  
  write_tsv(top_pair_section, file.path(out_qc, "top_pairs_by_section.tsv"))
  
  top_pair_consistency <- top_pair_section %>%
    group_by(resource, zone, pair) %>%
    summarise(
      n_sections = n(),
      mean_score = mean(score, na.rm = TRUE),
      sd_score = sd(score, na.rm = TRUE),
      cv_score = ifelse(mean_score > 0, sd_score / mean_score, NA_real_),
      min_score = min(score, na.rm = TRUE),
      max_score = max(score, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(resource, zone, cv_score)
  
  write_tsv(top_pair_consistency, file.path(out_qc, "top_pair_consistency.tsv"))
  
  # section-specific z score per pair
  top_pair_z <- top_pair_section %>%
    group_by(resource, zone, pair) %>%
    mutate(pair_z = as.numeric(scale(score))) %>%
    ungroup()
  
  write_tsv(top_pair_z, file.path(out_qc, "top_pairs_by_section_zscores.tsv"))
  
  pair_outlier_sections <- top_pair_z %>%
    filter(abs(pair_z) >= 2) %>%
    arrange(resource, zone, pair, desc(abs(pair_z)))
  
  write_tsv(pair_outlier_sections, file.path(out_qc, "pair_outlier_sections.tsv"))
}

# -------------------------
# 3) permutation significance per section
# -------------------------
if (nrow(perm_all) > 0) {
  perm_summary <- perm_all %>%
    group_by(resource, zone, section_id) %>%
    summarise(
      n_pairs = n(),
      n_sig_05 = sum(pval < 0.05, na.rm = TRUE),
      n_sig_01 = sum(pval < 0.01, na.rm = TRUE),
      frac_sig_05 = mean(pval < 0.05, na.rm = TRUE),
      mean_obs = mean(obs_score, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      z_frac_sig_05 = ave(frac_sig_05, resource, zone,
                          FUN = function(x) as.numeric(scale(x)))
    ) %>%
    arrange(resource, zone, desc(frac_sig_05))
  
  write_tsv(perm_summary, file.path(out_qc, "permutation_section_summary.tsv"))
  
  perm_outliers <- perm_summary %>%
    filter(abs(z_frac_sig_05) >= 2)
  
  write_tsv(perm_outliers, file.path(out_qc, "permutation_section_outliers.tsv"))
}

# -------------------------
# 4) quick plots
# -------------------------
p1 <- ggplot(section_burden, aes(x = factor(section_id, levels = good_sections), y = mean_score, fill = resource)) +
  geom_col(position = "dodge") +
  facet_wrap(~ zone, scales = "free_y") +
  theme_classic(base_size = 12) +
  labs(x = "Section", y = "Mean observed LR score")

ggsave(file.path(out_qc, "QC_mean_observed_score_by_section.png"),
       p1, width = 9, height = 4.8, dpi = 300, bg = "white")

if (exists("perm_summary")) {
  p2 <- ggplot(perm_summary, aes(x = factor(section_id, levels = good_sections), y = frac_sig_05, fill = resource)) +
    geom_col(position = "dodge") +
    facet_wrap(~ zone, scales = "free_y") +
    theme_classic(base_size = 12) +
    labs(x = "Section", y = "Fraction significant pairs (p < 0.05)")
  
  ggsave(file.path(out_qc, "QC_fraction_significant_pairs_by_section.png"),
         p2, width = 9, height = 4.8, dpi = 300, bg = "white")
}

message("Done. Cross-section LR QC written to: ", out_qc)