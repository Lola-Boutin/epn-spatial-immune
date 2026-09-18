# =========================================================
# Figure 1 — unified script
#
# Panels:
#   a  Mean inferred immune abundance by histology (all populations)
#   b  Mean inferred lymphoid abundance by histology (lymphoid only,
#        independently ordered by lymphoid total)
#   c  Cumulated total immune abundance by Aggressiveness
#   d  Cumulated lymphoid abundance by Aggressiveness
#   e  EPN T-cell abundance by refined anatomical location
#        (CD4T, CD8T, γδT; Posterior Fossa / Cortex / Ventricles /
#         Spinal cord)
#   f  EPN T-cell abundance by tumor type
#        (CD4T, CD8T, γδT; Primary / Progression / Recurrence)
#
# Data pipeline: reads raw CIBERSORTx outputs + metadata TSV,
# merges them, then produces all six panels.
# The merged TSV is also saved for downstream use.
# =========================================================

source(here::here("config.R"))

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
  library(rstatix)
})

# =========================================================
# 0. Paths
# =========================================================
input_dir <- figure_input_path("figure1")

meta_file <- file.path(
  input_dir,
  "merged_abundance_filtered_allimmune_with_tumor_type.tsv"
)
deepTIL_dir <- file.path(input_dir, "deeptil_results")

out_dir <- ensure_dir(figure_path("figure1"))

if (!file.exists(meta_file)) {
  stop(
    "Missing Figure 1 metadata: ", meta_file,
    "\nPlace the file under $EPN_FIGURE_INPUT_ROOT/figure1/ ",
    "or set EPN_FIGURE_INPUT_ROOT."
  )
}
if (!dir.exists(deepTIL_dir)) {
  stop(
    "Missing Figure 1 DeepTIL directory: ", deepTIL_dir,
    "\nExpected SES_CIBERSORTx_*.txt files under this directory."
  )
}

# =========================================================
# 1. Settings
# =========================================================
immune_cols   <- c("Bcells", "TCD4", "TCD8", "Tgd", "NK", "MoMaDC", "granulocytes")
lymphoid_cols <- c("Bcells", "TCD4", "TCD8", "Tgd", "NK")
epn_cols      <- c("TCD4", "TCD8", "Tgd")

histology_short <- c(
  "Atypical_Teratoid_Rhabdoid_Tumor"                      = "ATRT",
  "Atypical_choroid_plexus_papilloma"                     = "ACPP",
  "Choroid_plexus_carcinoma"                              = "CPC",
  "Choroid_plexus_papilloma"                              = "CPP",
  "Desmoplastic_infantile_astrocytoma_and_ganglioglioma"  = "DIA/DIG",
  "Diffuse_fibrillary_astrocytoma"                        = "DFA",
  "Diffuse_intrinsic_pontine_glioma"                      = "DIPG",
  "Diffuse_midline_glioma"                                = "DMG",
  "Embryonal_tumor_with_multilayer_rosettes"              = "ETMR",
  "Ependymoma"                                            = "EPN",
  "Ganglioglioma"                                         = "GG",
  "High-grade_glioma_astrocytoma"                         = "HGG",
  "Low-grade_glioma_astrocytoma"                          = "LGG",
  "Medulloblastoma"                                       = "MB",
  "Meningioma"                                            = "MNG",
  "Oligodendroglioma"                                     = "ODG",
  "Pilocytic_astrocytoma"                                 = "PA",
  "Pineoblastoma"                                         = "PB",
  "Pleomorphic_xanthoastrocytoma"                         = "PXA"
)

cell_labels <- c(
  "Bcells"       = "B",
  "TCD4"         = "CD4T",
  "TCD8"         = "CD8T",
  "Tgd"          = "gdT",
  "NK"           = "NK",
  "MoMaDC"       = "MoMaDC",
  "granulocytes" = "Granulocytes"
)

# Shared colour palette (panels a, b, e, f)
cell_cols <- c(
  "B"            = "#FF7F00",
  "CD4T"         = "#4DAF4A",
  "CD8T"         = "#377EB8",
  "gdT"          = "#E41A1C",
  "NK"           = "#984EA3",
  "MoMaDC"       = "#A65628",
  "Granulocytes" = "#999999"
)

# Aggressiveness palette (panels c, d)
aggr_levels <- c("Low", "Intermediate", "High", "Very High")
aggr_cols   <- c(
  "Low"          = "#74C476",
  "Intermediate" = "#FD8D3C",
  "High"         = "#E31A1C",
  "Very High"    = "#800026"
)
aggr_pairs <- list(
  c("Low", "Intermediate"), c("Low", "High"), c("Low", "Very High"),
  c("Intermediate", "High"), c("Intermediate", "Very High"),
  c("High", "Very High")
)

lymphoid_levels <- c("B", "CD4T", "CD8T", "gdT", "NK")

# =========================================================
# 2. Load and merge data
# =========================================================
message("Loading data...")
meta <- read.delim(meta_file, sep = "\t", check.names = FALSE)

files <- list.files(deepTIL_dir, pattern = "^SES_CIBERSORTx_.*\\.txt$",
                    full.names = TRUE)
if (!length(files)) stop("No SES_CIBERSORTx_*.txt files found in: ", deepTIL_dir)

deconv_all <- purrr::map_dfr(files, function(f) {
  histology_from_file <- gsub("\\.txt$", "",
                              gsub("^SES_CIBERSORTx_", "", basename(f)))
  read.delim(f, sep = "\t", check.names = FALSE) %>%
    dplyr::rename(SAMPLE_ID = Mixture, P_value = `P-value`) %>%
    dplyr::mutate(HISTOLOGY_from_file = histology_from_file) %>%
    dplyr::select(SAMPLE_ID, HISTOLOGY_from_file,
                  dplyr::all_of(immune_cols), P_value, Correlation, RMSE)
})

# Optional QC: check gdT correlation between old and new scores
if ("Tgd" %in% colnames(meta)) {
  qc_tgd <- meta %>%
    dplyr::select(SAMPLE_ID, Tgd_old = Tgd) %>%
    dplyr::left_join(deconv_all %>% dplyr::select(SAMPLE_ID, Tgd_new = Tgd),
                     by = "SAMPLE_ID") %>%
    dplyr::summarise(n_overlap    = sum(!is.na(Tgd_old) & !is.na(Tgd_new)),
                     cor          = cor(Tgd_old, Tgd_new, use = "complete.obs"),
                     max_abs_diff = max(abs(Tgd_old - Tgd_new), na.rm = TRUE))
  message("QC gdT correlation: r=", round(qc_tgd$cor, 3),
          "  max_diff=", round(qc_tgd$max_abs_diff, 4))
}

merged_allimmune <- meta %>%
  dplyr::select(-dplyr::any_of(c(immune_cols, "P_value", "Correlation",
                                  "RMSE", "HISTOLOGY_from_file"))) %>%
  dplyr::left_join(deconv_all, by = "SAMPLE_ID")

# Histology filename vs metadata QC
hist_mismatch <- merged_allimmune %>%
  dplyr::filter(!is.na(HISTOLOGY_from_file), HISTOLOGY != HISTOLOGY_from_file)
if (nrow(hist_mismatch) > 0)
  warning(nrow(hist_mismatch),
          " histology mismatches between metadata and CIBERSORTx filenames.")

write_tsv(merged_allimmune,
          file.path(out_dir,
                    "merged_abundance_filtered_allimmune_with_tumor_type.tsv"))

# =========================================================
# PANELS A + B — Mean abundance stacked bar charts
# =========================================================
message("Panels a/b")

plot_df <- merged_allimmune %>%
  dplyr::filter(dplyr::if_all(dplyr::all_of(immune_cols), ~ !is.na(.x))) %>%
  dplyr::mutate(HISTOLOGY_SHORT = dplyr::recode(HISTOLOGY, !!!histology_short,
                                                 .default = HISTOLOGY)) %>%
  dplyr::group_by(HISTOLOGY, HISTOLOGY_SHORT) %>%
  dplyr::summarise(n = dplyr::n(),
                   dplyr::across(dplyr::all_of(immune_cols),
                                 ~ mean(.x, na.rm = TRUE)),
                   .groups = "drop") %>%
  tidyr::pivot_longer(cols = dplyr::all_of(immune_cols),
                      names_to = "Immune_pop", values_to = "Mean_abundance") %>%
  dplyr::mutate(Immune_pop = dplyr::recode(Immune_pop, !!!cell_labels))

# Panel a: order by total immune abundance
hist_order_a <- plot_df %>%
  dplyr::group_by(HISTOLOGY_SHORT) %>%
  dplyr::summarise(total = sum(Mean_abundance), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(total)) %>%
  dplyr::pull(HISTOLOGY_SHORT)

plot_df_a <- plot_df %>%
  dplyr::mutate(HISTOLOGY_SHORT = factor(HISTOLOGY_SHORT, levels = hist_order_a),
                Immune_pop      = factor(Immune_pop, levels = names(cell_cols)))

write_tsv(plot_df_a, file.path(out_dir,
                                "Figure1a_mean_abundance_by_histology_long.tsv"))

p_a <- ggplot(plot_df_a, aes(HISTOLOGY_SHORT, Mean_abundance, fill = Immune_pop)) +
  geom_col(width = 0.85, color = "black", linewidth = 0.15) +
  scale_fill_manual(values = cell_cols, drop = FALSE) +
  labs(x = NULL, y = "Mean inferred abundance", fill = NULL) +
  theme_classic(base_size = 12) +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.y = element_text(face = "bold"),
        legend.position = "right")

ggsave(file.path(out_dir, "Figure1a_all_immune_by_histology.png"),
       p_a, width = 10, height = 5, dpi = 600, bg = "white")
ggsave(file.path(out_dir, "Figure1a_all_immune_by_histology.pdf"),
       p_a, width = 10, height = 5, device = cairo_pdf, bg = "white")

# Panel b: lymphoid only, order by lymphoid abundance
plot_df_b <- plot_df %>%
  dplyr::filter(Immune_pop %in% lymphoid_levels) %>%
  dplyr::mutate(Immune_pop = factor(Immune_pop, levels = lymphoid_levels))

hist_order_b <- plot_df_b %>%
  dplyr::group_by(HISTOLOGY_SHORT) %>%
  dplyr::summarise(total = sum(Mean_abundance), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(total)) %>%
  dplyr::pull(HISTOLOGY_SHORT)

plot_df_b <- plot_df_b %>%
  dplyr::mutate(HISTOLOGY_SHORT = factor(HISTOLOGY_SHORT, levels = hist_order_b))

write_tsv(plot_df_b, file.path(out_dir,
                                "Figure1b_lymphoid_by_histology_long.tsv"))

p_b <- ggplot(plot_df_b, aes(HISTOLOGY_SHORT, Mean_abundance, fill = Immune_pop)) +
  geom_col(width = 0.85, color = "black", linewidth = 0.15) +
  scale_fill_manual(values = cell_cols[lymphoid_levels], drop = FALSE) +
  labs(x = NULL, y = "Mean inferred abundance", fill = NULL) +
  theme_classic(base_size = 12) +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.y = element_text(face = "bold"),
        legend.position = "right")

ggsave(file.path(out_dir, "Figure1b_lymphoid_by_histology.png"),
       p_b, width = 10, height = 5, dpi = 600, bg = "white")
ggsave(file.path(out_dir, "Figure1b_lymphoid_by_histology.pdf"),
       p_b, width = 10, height = 5, device = cairo_pdf, bg = "white")

# =========================================================
# PANELS C + D — Cumulated abundance by Aggressiveness
# =========================================================
message("Panels c/d")

df_scores <- merged_allimmune %>%
  dplyr::filter(dplyr::if_all(dplyr::all_of(immune_cols), ~ !is.na(.x)),
                !is.na(Aggressiveness)) %>%
  dplyr::mutate(
    Aggressiveness = factor(Aggressiveness, levels = aggr_levels),
    total_immune   = rowSums(dplyr::across(dplyr::all_of(immune_cols))),
    total_lymphoid = rowSums(dplyr::across(dplyr::all_of(lymphoid_cols)))
  )

# Kruskal-Wallis global tests
kw_immune   <- df_scores %>% rstatix::kruskal_test(total_immune   ~ Aggressiveness)
kw_lymphoid <- df_scores %>% rstatix::kruskal_test(total_lymphoid ~ Aggressiveness)
message("Kruskal-Wallis total immune:   p=", kw_immune$p)
message("Kruskal-Wallis total lymphoid: p=", kw_lymphoid$p)

# Pairwise Wilcoxon with BH correction
pw_immune <- df_scores %>%
  rstatix::wilcox_test(total_immune ~ Aggressiveness,
                       p.adjust.method = "BH", comparisons = aggr_pairs) %>%
  rstatix::add_significance() %>%
  rstatix::add_xy_position(x = "Aggressiveness")

pw_lymphoid <- df_scores %>%
  rstatix::wilcox_test(total_lymphoid ~ Aggressiveness,
                       p.adjust.method = "BH", comparisons = aggr_pairs) %>%
  rstatix::add_significance() %>%
  rstatix::add_xy_position(x = "Aggressiveness")

# Sample counts for x-axis labels
n_labels <- df_scores %>%
  dplyr::count(Aggressiveness) %>%
  dplyr::mutate(label = paste0(Aggressiveness, "\n(n=", n, ")"))
aggr_label_map <- setNames(n_labels$label, as.character(n_labels$Aggressiveness))

make_aggr_plot <- function(score_col, y_label, pw_result) {
  ggplot(df_scores, aes(x = Aggressiveness, y = .data[[score_col]],
                        fill = Aggressiveness)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.85, width = 0.55,
                 color = "black", linewidth = 0.35) +
    geom_jitter(aes(color = Aggressiveness), width = 0.15,
                size = 1.2, alpha = 0.5) +
    ggpubr::stat_pvalue_manual(pw_result, label = "p.adj.signif",
                               tip.length = 0.01, step.increase = 0.06,
                               hide.ns = TRUE) +
    scale_x_discrete(labels = aggr_label_map) +
    scale_fill_manual(values  = aggr_cols) +
    scale_color_manual(values = aggr_cols) +
    labs(x = NULL, y = y_label) +
    theme_classic(base_size = 13) +
    theme(axis.title.y    = element_text(face = "bold"),
          axis.text.x     = element_text(size = 11),
          legend.position = "none")
}

p_c <- make_aggr_plot("total_immune",   "Cumulated inferred immune abundance",   pw_immune)
p_d <- make_aggr_plot("total_lymphoid", "Cumulated inferred lymphoid abundance", pw_lymphoid)

ggsave(file.path(out_dir, "Figure1c_total_immune_by_aggressiveness.png"),
       p_c, width = 7, height = 5, dpi = 600, bg = "white")
ggsave(file.path(out_dir, "Figure1c_total_immune_by_aggressiveness.pdf"),
       p_c, width = 7, height = 5, device = cairo_pdf, bg = "white")

ggsave(file.path(out_dir, "Figure1d_lymphoid_by_aggressiveness.png"),
       p_d, width = 7, height = 5, dpi = 600, bg = "white")
ggsave(file.path(out_dir, "Figure1d_lymphoid_by_aggressiveness.pdf"),
       p_d, width = 7, height = 5, device = cairo_pdf, bg = "white")

# =========================================================
# PANELS E + F — EPN T-cell abundance (shared data prep)
# =========================================================
message("Panels e/f")

epn_df <- merged_allimmune %>%
  dplyr::filter(HISTOLOGY == "Ependymoma") %>%
  dplyr::filter(dplyr::if_all(dplyr::all_of(epn_cols), ~ !is.na(.x))) %>%
  tidyr::pivot_longer(cols = dplyr::all_of(epn_cols),
                      names_to = "Immune_pop", values_to = "Abundance") %>%
  dplyr::mutate(
    Immune_pop = dplyr::recode(Immune_pop, !!!cell_labels),
    Immune_pop = factor(Immune_pop, levels = c("CD4T", "CD8T", "gdT"))
  )

# Panel e — by REFINED_LOCATION (4 locations)
epn_df_e <- epn_df %>%
  dplyr::filter(REFINED_LOCATION %in% c("Cerebellum/Posterior Fossa",
                                         "Cortex", "Ventricles", "Spinal cord")) %>%
  dplyr::mutate(REFINED_LOCATION = factor(
    REFINED_LOCATION,
    levels = c("Cerebellum/Posterior Fossa", "Cortex", "Ventricles", "Spinal cord")
  ))

stat_e <- epn_df_e %>%
  dplyr::group_by(Immune_pop) %>%
  rstatix::kruskal_test(Abundance ~ REFINED_LOCATION) %>%
  rstatix::adjust_pvalue(method = "BH") %>%
  rstatix::add_significance()

pairwise_e <- epn_df_e %>%
  dplyr::group_by(Immune_pop) %>%
  rstatix::wilcox_test(Abundance ~ REFINED_LOCATION, p.adjust.method = "BH") %>%
  rstatix::add_significance() %>%
  rstatix::add_xy_position(x = "REFINED_LOCATION", dodge = 0.8)

message("\n--- Panel e: Kruskal-Wallis (EPN by REFINED_LOCATION) ---")
print(stat_e)

p_e <- ggplot(epn_df_e, aes(REFINED_LOCATION, Abundance, fill = Immune_pop)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.6,
               color = "black", linewidth = 0.35) +
  geom_jitter(aes(color = Immune_pop),
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.6),
              size = 1, alpha = 0.6) +
  ggpubr::stat_pvalue_manual(pairwise_e, label = "p.adj.signif",
                             tip.length = 0.01, step.increase = 0.05,
                             hide.ns = TRUE) +
  facet_wrap(~ Immune_pop, scales = "free_y", nrow = 1) +
  scale_fill_manual(values  = cell_cols[c("CD4T", "CD8T", "gdT")]) +
  scale_color_manual(values = cell_cols[c("CD4T", "CD8T", "gdT")]) +
  labs(x = "Refined location", y = "Inferred abundance",
       fill = NULL, color = NULL) +
  theme_classic(base_size = 12) +
  theme(axis.text.x      = element_text(angle = 45, hjust = 1, vjust = 1, size = 9),
        axis.title       = element_text(face = "bold"),
        strip.text       = element_text(face = "bold", size = 11),
        strip.background = element_blank(),
        legend.position  = "none")

ggsave(file.path(out_dir, "Figure1e_EPN_Tcells_by_location.png"),
       p_e, width = 14, height = 5, dpi = 600, bg = "white")
ggsave(file.path(out_dir, "Figure1e_EPN_Tcells_by_location.pdf"),
       p_e, width = 14, height = 5, device = cairo_pdf, bg = "white")

# Panel f — by TUMOR_TYPE
epn_df_f <- epn_df %>%
  dplyr::mutate(
    TUMOR_TYPE = stringr::str_to_title(TUMOR_TYPE),
    TUMOR_TYPE = factor(TUMOR_TYPE,
                        levels = c("Primary", "Progression", "Recurrence"))
  )

pairwise_f <- epn_df_f %>%
  dplyr::group_by(Immune_pop) %>%
  rstatix::wilcox_test(Abundance ~ TUMOR_TYPE, p.adjust.method = "BH") %>%
  rstatix::add_significance() %>%
  rstatix::add_xy_position(x = "TUMOR_TYPE", dodge = 0.8)

message("\n--- Panel f: Wilcoxon (EPN by TUMOR_TYPE) ---")
print(epn_df_f %>%
        dplyr::group_by(Immune_pop) %>%
        rstatix::wilcox_test(Abundance ~ TUMOR_TYPE, p.adjust.method = "BH") %>%
        rstatix::add_significance())

p_f <- ggplot(epn_df_f, aes(TUMOR_TYPE, Abundance, fill = Immune_pop)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.6,
               color = "black", linewidth = 0.35) +
  geom_jitter(aes(color = Immune_pop),
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.6),
              size = 1.5, alpha = 0.6) +
  ggpubr::stat_pvalue_manual(pairwise_f, label = "p.adj.signif",
                             tip.length = 0.01, step.increase = 0.08,
                             hide.ns = TRUE) +
  facet_wrap(~ Immune_pop, scales = "free_y", nrow = 1) +
  scale_fill_manual(values  = cell_cols[c("CD4T", "CD8T", "gdT")]) +
  scale_color_manual(values = cell_cols[c("CD4T", "CD8T", "gdT")]) +
  labs(x = "Tumor type", y = "Inferred abundance",
       fill = NULL, color = NULL) +
  theme_classic(base_size = 12) +
  theme(axis.text.x      = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
        axis.title       = element_text(face = "bold"),
        strip.text       = element_text(face = "bold", size = 11),
        strip.background = element_blank(),
        legend.position  = "none")

ggsave(file.path(out_dir, "Figure1f_EPN_Tcells_by_tumor_type.png"),
       p_f, width = 9, height = 5, dpi = 600, bg = "white")
ggsave(file.path(out_dir, "Figure1f_EPN_Tcells_by_tumor_type.pdf"),
       p_f, width = 9, height = 5, device = cairo_pdf, bg = "white")

message("\nFigure 1 done. Outputs saved to: ", out_dir)
