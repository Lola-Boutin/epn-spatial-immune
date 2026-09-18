# =============================================================================
# 11b_ranking_plots.R
#
# Purpose
# -------
# Supplementary ranking figures from the stage 11 outputs: volcano plots,
# effect concordance, and program barcode plots per target.
#
# Reads only; recomputes no rankings. Writes alongside the stage 11 tables so
# the panels sit with the data they describe.
#
# Inputs
# ------
#   stage_dir("gsea")   tables/11_gene_ranking_<target>.tsv and 11_gsea_* tables
#
# Outputs
# -------
#   stage_dir("gsea")
#     tables/11b_program_ranks_<target>.tsv, 11b_top_nonprogram_<target>.tsv,
#     tables/11b_program_rank_test.tsv
#     plots/11b_volcano_<target>.{png,pdf}
#     plots/11b_effect_concordance_<target>.{png,pdf}
#     plots/11b_program_barcode_<target>.{png,pdf}
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
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

# =========================================================
# 11b — Gene-level views of the opportunity-score rankings
#
# Reads the ranking tables written by script 11 and produces three views per
# target. It recomputes nothing, so it reruns in seconds.
#
# CHANGE FROM THE FIRST VERSION
# -----------------------------
# Program membership is now defined PER TARGET. Script 11 flags any of the 47
# declared genes as a program gene regardless of which program it belongs to,
# so its rank test asked whether ANY program gene ranked highly - and because
# the alphabeta core and the inhibitory program are collinear, their genes
# inflate each other's result. Here each target is tested against its OWN
# members only, which is the question that matters: does a score recover the
# genes it was built from?
#
# THE THREE VIEWS
# ---------------
# 1. VOLCANO. x = median within-section Spearman rho, y = -log10(p) from the
#    combined Stouffer z. Familiar, but the p tests whether the weighted
#    combination of within-section correlations differs from zero: it summarises
#    14 sections, not spots, and does not account for the score being derived
#    from expression. Effect size and concordance are more trustworthy.
#
# 2. EFFECT VERSUS CONCORDANCE. x = combined z, y = sections agreeing in
#    direction. At n = 14 this is the honest display: a gene at z = 0.4
#    concordant in 14/14 sections is better supported than one at z = 0.8
#    concordant in 9/14, and a volcano hides that.
#
# 3. PROGRAM-GENE BARCODE. Where the target's own program genes fall in the
#    ranking, with a rank-sum test against all other genes. A positive control:
#    a score should rank its own members at the top. The shape also shows how
#    steeply the signal falls away once they are passed, distinguishing a score
#    with a broad transcriptional neighbourhood from one that is essentially its
#    own gene list.
#
# OUTPUTS
#   tables/11b_program_gene_ranks_<target>.tsv
#   tables/11b_top_nonprogram_<target>.tsv
#   tables/11b_program_rank_test.tsv
#   plots/11b_volcano_<target>.png | .pdf
#   plots/11b_effect_concordance_<target>.png | .pdf
#   plots/11b_program_barcode_<target>.png | .pdf
# =========================================================


# -------------------------
# CONFIG
# -------------------------
p11_root   <- stage_dir("gsea")
p11_tables <- file.path(p11_root, "tables")

out_tables <- ensure_dir(file.path(p11_root, "tables"))
out_plots  <- ensure_dir(file.path(p11_root, "plots"))

TARGETS <- c("gdT_core", "abT_core", "inhibitory_program")

N_LABEL_TOP     <- 12
N_LABEL_NONPROG <- 10
N_TOP_TABLE     <- 40

TARGET_LABELS <- c(
  gdT_core           = "\u03b3\u03b4T ligand availability",
  abT_core           = "\u03b1\u03b2T ligand availability",
  inhibitory_program = "Immunosuppressive tone"
)

# Per-target membership. Both alias variants listed so the flag is complete
# whichever symbol the reference carries.
PROGRAM_MEMBERS <- list(
  abT_core = c("HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2",
               "PSMB8", "PSMB9", "ERAP1", "ERAP2", "NLRC5", "IRF1"),
  gdT_core = c("BTN2A1", "BTN3A1", "MVK", "PMVK", "MVD", "IDI1",
               "MICA", "MICB", "ULBP2", "ULBP3",
               "ULBP4", "RAET1E", "ULBP5", "RAET1G", "ULBP6", "RAET1L",
               "EPHA2", "PVR", "NECTIN2", "ICAM1"),
  inhibitory_program = c("CD274", "PDCD1LG2", "HLA-E", "LGALS9", "IDO1",
                         "TGFB1", "TGFB2", "TGFB3", "IL10", "VEGFA",
                         "NT5E", "ENTPD1", "PTGS2", "PTGES", "PTGER2", "PTGER4")
)
PROGRAM_MEMBERS <- lapply(PROGRAM_MEMBERS, toupper)

# Genes in the other two programs, flagged separately so cross-program
# collinearity is visible rather than hidden
ALL_PROGRAM_GENES <- toupper(unique(c(unlist(PROGRAM_MEMBERS), "FDPS")))

COL_OWN   <- "#F07167"
COL_OTHER_PROG <- "#A3A500"
COL_REST  <- "#9ecae1"

CLASS_LEVELS <- c("own program", "other program", "other gene")
CLASS_COLS <- c("own program" = COL_OWN,
                "other program" = COL_OTHER_PROG,
                "other gene" = COL_REST)

target_label <- function(x) ifelse(x %in% names(TARGET_LABELS), TARGET_LABELS[x], x)

# =========================================================
# LOOP OVER TARGETS
# =========================================================
rank_test_rows <- list()

for (tgt in TARGETS) {
  
  f <- file.path(p11_tables, paste0("11_gene_ranking_", tgt, ".tsv"))
  if (!file.exists(f)) { warning("Missing ranking table: ", f); next }
  
  message("\n=====================================================")
  message("Target: ", tgt)
  message("=====================================================")
  
  own <- PROGRAM_MEMBERS[[tgt]]
  if (is.null(own)) { warning("No program membership defined for ", tgt); next }
  
  rk <- read_tsv(f, show_col_types = FALSE) %>%
    filter(is.finite(z)) %>%
    mutate(
      gene = toupper(gene),
      p_stouffer = 2 * stats::pnorm(-abs(z)),
      neglog10p  = -log10(pmax(p_stouffer, .Machine$double.xmin)),
      concordance = k_concordant / n_tested,
      is_own = gene %in% own,
      is_any_program = gene %in% ALL_PROGRAM_GENES,
      class = factor(dplyr::case_when(
        is_own ~ "own program",
        is_any_program ~ "other program",
        TRUE ~ "other gene"), levels = CLASS_LEVELS)
    ) %>%
    arrange(desc(z)) %>%
    mutate(rank_pos = row_number())
  
  rk$p_adj <- p.adjust(rk$p_stouffer, method = "BH")
  
  message("Genes ranked: ", nrow(rk),
          " | own-program genes present: ", sum(rk$is_own),
          " | genes from other programs: ", sum(rk$is_any_program & !rk$is_own))
  
  # -------------------------
  # TABLES
  # -------------------------
  prog_ranks <- rk %>%
    filter(is_own) %>%
    select(gene, rank_pos, z, median_rho, k_concordant, n_tested, p_adj) %>%
    arrange(rank_pos)
  write_tsv(prog_ranks, file.path(out_tables,
                                  paste0("11b_program_gene_ranks_", tgt, ".tsv")))
  message("\n--- Own-program genes and their rank positions ---")
  print(prog_ranks, n = Inf)
  
  top_nonprog <- rk %>%
    filter(!is_any_program) %>%
    arrange(desc(abs(z))) %>%
    slice_head(n = N_TOP_TABLE) %>%
    select(gene, rank_pos, z, median_rho, k_concordant, n_tested, p_adj)
  write_tsv(top_nonprog, file.path(out_tables,
                                   paste0("11b_top_nonprogram_", tgt, ".tsv")))
  message("\n--- Top genes outside every scoring program ---")
  print(top_nonprog %>% slice_head(n = 15), n = Inf)
  
  # rank-sum test: own-program genes against everything else
  if (sum(rk$is_own) >= 3) {
    wt <- suppressWarnings(wilcox.test(rank_pos ~ is_own, data = rk))
    rank_test_rows[[tgt]] <- tibble(
      target = tgt,
      n_own_program = sum(rk$is_own),
      n_other = sum(!rk$is_own),
      median_rank_own = median(rk$rank_pos[rk$is_own]),
      median_rank_other = median(rk$rank_pos[!rk$is_own]),
      best_rank_own = min(rk$rank_pos[rk$is_own]),
      p_value = wt$p.value)
  }
  
  # -------------------------
  # 1. VOLCANO
  # -------------------------
  lab_v <- bind_rows(
    rk %>% arrange(desc(abs(z))) %>% slice_head(n = N_LABEL_TOP),
    rk %>% filter(!is_any_program) %>% arrange(desc(abs(z))) %>%
      slice_head(n = N_LABEL_NONPROG)
  ) %>% distinct(gene, .keep_all = TRUE)
  
  p_volc <- ggplot(rk, aes(x = median_rho, y = neglog10p, colour = class)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey70") +
    geom_point(data = rk %>% filter(class == "other gene"),
               size = 0.7, alpha = 0.45) +
    geom_point(data = rk %>% filter(class != "other gene"), size = 2.1) +
    ggrepel::geom_text_repel(data = lab_v, aes(label = gene), size = 2.8,
                             max.overlaps = Inf, min.segment.length = 0,
                             show.legend = FALSE) +
    scale_colour_manual(values = CLASS_COLS, name = NULL, drop = FALSE) +
    labs(x = "Median within-section Spearman \u03c1",
         y = expression(-log[10]~"p (combined across sections)"),
         title = paste0("Genes associated with ", target_label(tgt)),
         subtitle = paste0("p from the weighted combination of within-section ",
                           "correlations; it summarises sections, not spots.\n",
                           "This score's own program genes are expected to rank ",
                           "highly - that is a positive control, not a finding.")) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(size = 8),
          legend.position = "bottom")
  
  ggsave(file.path(out_plots, paste0("11b_volcano_", tgt, ".png")), p_volc,
         width = 8.5, height = 7, dpi = 600, bg = "white")
  ggsave(file.path(out_plots, paste0("11b_volcano_", tgt, ".pdf")), p_volc,
         width = 8.5, height = 7, bg = "white", device = cairo_pdf)
  
  # -------------------------
  # 2. EFFECT VERSUS CONCORDANCE
  # -------------------------
  n_max <- max(rk$n_tested, na.rm = TRUE)
  
  lab_c <- bind_rows(
    rk %>% arrange(desc(abs(z))) %>% slice_head(n = N_LABEL_TOP),
    rk %>% filter(!is_any_program, k_concordant == n_tested) %>%
      arrange(desc(abs(z))) %>% slice_head(n = N_LABEL_NONPROG)
  ) %>% distinct(gene, .keep_all = TRUE)
  
  p_conc <- ggplot(rk, aes(x = z, y = k_concordant, colour = class)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey70") +
    geom_hline(yintercept = n_max, linetype = 3, colour = "grey50") +
    geom_jitter(data = rk %>% filter(class == "other gene"),
                width = 0, height = 0.18, size = 0.7, alpha = 0.4) +
    geom_point(data = rk %>% filter(class != "other gene"), size = 2.1) +
    ggrepel::geom_text_repel(data = lab_c, aes(label = gene), size = 2.8,
                             max.overlaps = Inf, min.segment.length = 0,
                             show.legend = FALSE) +
    scale_colour_manual(values = CLASS_COLS, name = NULL, drop = FALSE) +
    scale_y_continuous(breaks = seq(0, n_max, by = 2)) +
    labs(x = "Combined z (effect across sections)",
         y = "Sections with concordant direction",
         title = paste0("Effect size versus cross-section consistency: ",
                        target_label(tgt)),
         subtitle = paste0("Dotted line: concordant in every section. A modest ",
                           "effect seen in all sections is better supported than ",
                           "a larger effect seen in half of them.")) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(size = 8),
          legend.position = "bottom")
  
  ggsave(file.path(out_plots, paste0("11b_effect_concordance_", tgt, ".png")), p_conc,
         width = 8.5, height = 6.5, dpi = 600, bg = "white")
  ggsave(file.path(out_plots, paste0("11b_effect_concordance_", tgt, ".pdf")), p_conc,
         width = 8.5, height = 6.5, bg = "white", device = cairo_pdf)
  
  # -------------------------
  # 3. PROGRAM-GENE BARCODE
  # -------------------------
  n_genes <- nrow(rk)
  
  bar_top <- ggplot(rk, aes(x = rank_pos, y = z)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey70") +
    geom_line(colour = "grey40", linewidth = 0.4) +
    geom_point(data = rk %>% filter(is_own), colour = COL_OWN, size = 1.8) +
    ggrepel::geom_text_repel(
      data = rk %>% filter(is_own) %>% slice_head(n = 10),
      aes(label = gene), size = 2.6, max.overlaps = Inf,
      min.segment.length = 0, direction = "y", nudge_y = 0.15) +
    scale_x_continuous(limits = c(1, n_genes), expand = c(0.01, 0)) +
    labs(x = NULL, y = "Combined z",
         title = paste0("Where this score's own program genes sit: ",
                        target_label(tgt)),
         subtitle = paste0(n_genes, " genes ranked by association with the score")) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(size = 8),
          axis.text.x = element_blank(), axis.ticks.x = element_blank())
  
  bar_ticks <- ggplot(rk %>% filter(is_own),
                      aes(x = rank_pos, xend = rank_pos, y = 0, yend = 1)) +
    geom_segment(colour = COL_OWN, linewidth = 0.5) +
    scale_x_continuous(limits = c(1, n_genes), expand = c(0.01, 0)) +
    labs(x = "Rank (most positively to most negatively associated)", y = NULL) +
    theme_bw(base_size = 10) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          panel.grid = element_blank())
  
  p_bar <- bar_top / bar_ticks + plot_layout(heights = c(5, 1))
  
  ggsave(file.path(out_plots, paste0("11b_program_barcode_", tgt, ".png")), p_bar,
         width = 9, height = 5.4, dpi = 600, bg = "white")
  ggsave(file.path(out_plots, paste0("11b_program_barcode_", tgt, ".pdf")), p_bar,
         width = 9, height = 5.4, bg = "white", device = cairo_pdf)
}

# =========================================================
# RANK TEST ACROSS TARGETS
# =========================================================
if (length(rank_test_rows) > 0) {
  rank_test <- bind_rows(rank_test_rows) %>%
    mutate(p_adj = p.adjust(p_value, "BH"))
  write_tsv(rank_test, file.path(out_tables, "11b_program_rank_test.tsv"))
  message("\n--- Does each score rank its OWN program genes above the rest? ---")
  print(rank_test, n = Inf)
  message("A much lower median rank for own-program genes confirms the score ",
          "recovers its own members. This is a positive control, not a finding.")
}

message("\nDone. Outputs written to: ", p11_root)