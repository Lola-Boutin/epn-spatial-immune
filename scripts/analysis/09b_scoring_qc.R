# =============================================================================
# 09b_scoring_qc.R
#
# Purpose
# -------
# Methodological QC for the opportunity scores. Reads the cache written by
# 09_opportunity_map.R, so it never redefines a gene set and never repeats the
# expression pass.
#
# Run when the scoring method, the detection floor, or the gene programs change
# -- not on every rerun of 09.
#
# QUESTIONS ANSWERED
#  1. Gene dominance     Is a composite driven by two abundant genes, or do all
#                       members contribute? (why z-mean over raw-mean)
#  2. Depth coupling     Do scores track sequencing depth rather than biology?
#                       (why z-mean over UCell)
#  3. Scheme sensitivity How much does classification depend on the aggregation
#                       choice? Reported as a RESULT, not a diagnostic --
#                       agreement is not complete.
#  4. Axis independence  Are the mechanism axes uncorrelated? They are expected
#                       to be: the composites are formative indices. Low
#                       correlation supports the design; it is not a coherence
#                       failure.
#  5. Split-half         Is a section's score a stable property of the tissue?
#                       If this fails, nothing downstream survives.
#
# Inputs
# ------
#   stage_input("scoring", "rds/09_scoring_inputs.rds")
#   stage_input("scoring", "rds/09_section_summary.rds")
#   stage_dir("raw_vis")   only when RUN_UCELL is TRUE
#
# Outputs
# -------
#   stage_dir("scoring_qc")
#     tables/09b_gene_dominance.tsv, 09b_depth_correlation.tsv,
#     tables/09b_scheme_scores.tsv, 09b_scheme_category_agreement.tsv,
#     tables/09b_axis_independence.tsv, 09b_split_half.tsv
#     plots/S_QC1_gene_dominance.{png,pdf} ... S_QC4_split_half.{png,pdf}
#
# Stochastic
# ----------
#   YES -- use_seed("split_half"), 20260812 as published.
#   N_SPLIT_REP random half-splits of each section's spots.
#
# Runtime
# -------
#   ~25 minutes
# =============================================================================

source(here::here("config.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})


# -------------------------
# CONFIG
# -------------------------
raw_vis_dir  <- stage_dir("raw_vis")
phase09_root <- stage_dir("scoring")
phase09_rds  <- file.path(phase09_root, "rds")

out_root   <- stage_dir("scoring_qc")
out_tables <- ensure_dir(file.path(out_root, "tables"))
out_plots  <- ensure_dir(file.path(out_root, "plots"))

# UCell needs the FULL transcriptome to rank genes within a spot; the cached
# matrix holds only program genes, so enabling this reloads the raw objects.
# Off by default -- the depth result below is what settled the choice.
RUN_UCELL      <- TRUE
UCELL_MAX_RANK <- 1500

N_SPLIT_REP <- 100     # split-half repetitions
SEED        <- SEEDS$split_half

AXIS_LAB_AB <- "\u03b1\u03b2T net"
AXIS_LAB_GD <- "\u03b3\u03b4T net"

use_seed("split_half")

# -------------------------
# LOAD CACHE FROM 09
# -------------------------
cache_file <- file.path(phase09_rds, "09_scoring_inputs.rds")
if (!file.exists(cache_file)) {
  stop("Missing cache from script 09: ", cache_file,
       "\nRerun 09_opportunity_map.R - it writes this file.")
}
cache <- readRDS(cache_file)

expr_prog    <- cache$expr
meta_all     <- cache$meta
gene_sets    <- cache$gene_sets
PROGRAM_AXES <- cache$program_axes
program_names <- cache$programs
axis_names    <- cache$axes
detection     <- cache$detection
cfg           <- cache$config

section_summary <- readRDS(file.path(phase09_rds, "09_section_summary.rds"))

message("Loaded ", nrow(expr_prog), " spots x ", ncol(expr_prog), " genes")
message("Config: ", cfg$scoring_method, " / ", cfg$scoring_aggregation,
        " / floor ", cfg$detection_min_pct, "%")

# -------------------------
# SCORING (mirrors 09; kept here only so alternatives can be computed)
# -------------------------
z_vec <- function(x) {
  ok <- is.finite(x)
  if (sum(ok) < 2) return(rep(NA_real_, length(x)))
  s <- stats::sd(x[ok]); if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mean(x[ok])) / s
}

safe_scale <- function(X) {
  mu  <- colMeans(X, na.rm = TRUE)
  sdv <- apply(X, 2, stats::sd, na.rm = TRUE)
  sdv[!is.finite(sdv) | sdv == 0] <- NA_real_
  Z <- sweep(sweep(X, 2, mu, "-"), 2, sdv, "/")
  Z[!is.finite(Z)] <- 0
  Z
}

score_set_scheme <- function(X, scheme) {
  if (is.null(X) || ncol(X) == 0) return(rep(NA_real_, nrow(X)))
  X[is.na(X)] <- 0
  if (scheme == "rawmean") return(rowMeans(X))
  Z <- safe_scale(X)
  if (scheme == "zmean") return(rowMeans(Z))
  if (scheme == "pc1") {
    if (ncol(Z) < 2) return(as.numeric(Z))
    pc <- stats::prcomp(Z, center = FALSE, scale. = FALSE)
    s  <- as.numeric(Z %*% pc$rotation[, 1])
    if (suppressWarnings(stats::cor(s, rowMeans(Z))) < 0) s <- -s
    return(s)
  }
  stop("Unknown scheme: ", scheme)
}

score_program_scheme <- function(program, scheme, aggregation = cfg$scoring_aggregation) {
  ax <- PROGRAM_AXES[[program]]
  if (aggregation == "gene" || is.null(ax)) {
    g <- gene_sets[[program]]
    if (length(g) == 0) return(rep(NA_real_, nrow(expr_prog)))
    return(score_set_scheme(expr_prog[, g, drop = FALSE], scheme))
  }
  ax_use <- ax[sapply(ax, function(a) length(gene_sets[[a]]) > 0)]
  if (length(ax_use) == 0) return(rep(NA_real_, nrow(expr_prog)))
  M <- sapply(ax_use, function(a)
    z_vec(score_set_scheme(expr_prog[, gene_sets[[a]], drop = FALSE], scheme)))
  if (is.null(dim(M))) M <- matrix(M, ncol = 1)
  rowMeans(M, na.rm = TRUE)
}

classify_opportunity <- function(ab, gd, inh, threshold = cfg$suppressed_inhib_z) {
  out <- dplyr::case_when(
    ab >  0 & gd >  0 ~ "Dual-opportunity",
    ab >  0 & gd <= 0 ~ "abT-favorable",
    ab <= 0 & gd >  0 ~ "gdT-favorable",
    ab <= 0 & gd <= 0 & inh >= threshold ~ "Suppressed",
    TRUE ~ "Immuno-cold")
  factor(out, levels = c("abT-favorable", "Dual-opportunity", "gdT-favorable",
                         "Immuno-cold", "Suppressed"))
}

cor_long <- function(M, label = NA_character_) {
  cm <- suppressWarnings(stats::cor(M, method = "spearman", use = "pairwise.complete.obs"))
  tibble(group = label,
         var_a = rep(rownames(cm), times = ncol(cm)),
         var_b = rep(colnames(cm), each  = nrow(cm)),
         rho   = as.numeric(cm))
}

# =========================================================
# 1. GENE DOMINANCE
#
# Share of a raw-mean score contributed by each gene. Under raw-mean scoring an
# abundant gene dominates; under z-mean every gene contributes equally. This is
# the evidence for the scoring choice.
# =========================================================
dominance <- bind_rows(lapply(program_names, function(nm) {
  g <- gene_sets[[nm]]
  if (length(g) == 0) return(NULL)
  detection %>% filter(gene %in% g) %>%
    mutate(program = nm,
           share_rawmean = 100 * mean_lognorm / sum(mean_lognorm),
           share_zmean = 100 / length(g)) %>%
    select(program, gene, pct_detected, mean_lognorm, share_rawmean, share_zmean)
})) %>% arrange(program, desc(share_rawmean))

write_tsv(dominance, file.path(out_tables, "09b_gene_dominance.tsv"))

message("\n--- Top raw-mean contributors per program ---")
print(dominance %>% group_by(program) %>% slice_head(n = 3) %>%
        select(program, gene, pct_detected, share_rawmean, share_zmean), n = Inf)

p_dom <- dominance %>%
  pivot_longer(c(share_rawmean, share_zmean), names_to = "scheme", values_to = "share") %>%
  mutate(scheme = recode(scheme, share_rawmean = "raw mean", share_zmean = "z-mean")) %>%
  ggplot(aes(x = reorder(gene, share), y = share, fill = scheme)) +
  geom_col(position = "dodge", width = 0.75) +
  coord_flip() +
  facet_wrap(~program, scales = "free_y") +
  scale_fill_manual(values = c(`raw mean` = "#F07167", `z-mean` = "#20B3B7"), name = NULL) +
  labs(x = NULL, y = "Share of composite score (%)",
       title = "Gene contribution under each scoring scheme",
       subtitle = "Raw-mean scoring is dominated by the most abundant members; z-mean weights genes equally") +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

ggsave(file.path(out_plots, "S_QC1_gene_dominance.png"), p_dom,
       width = 12, height = 7, dpi = 600, bg = "white")
ggsave(file.path(out_plots, "S_QC1_gene_dominance.pdf"), p_dom,
       width = 12, height = 7, bg = "white", device = cairo_pdf)

# =========================================================
# 2. SEQUENCING DEPTH
#
# Needs nCount_Spatial, which is not in the cached program-gene matrix, so the
# raw objects are reopened. UCell scoring, if requested, happens in the same
# pass because it also needs the full transcriptome.
# =========================================================
depth_list <- list(); ucell_list <- list()

for (sec in unique(meta_all$section_id)) {
  vis_file <- file.path(raw_vis_dir, paste0("vis_section_", sec, "_raw.rds"))
  if (!file.exists(vis_file)) { warning("Missing vis file: ", sec); next }
  message("Depth pass: section ", sec)

  obj <- readRDS(vis_file)
  keys <- paste(sec, colnames(obj), sep = "|")
  depth_list[[sec]] <- tibble(key = keys,
                              nCount = obj$nCount_Spatial,
                              nFeature = obj$nFeature_Spatial)


  if (RUN_UCELL) {
    if (!requireNamespace("UCell", quietly = TRUE))
      stop("UCell not installed; set RUN_UCELL <- FALSE or install it")
    
    DefaultAssay(obj) <- "Spatial"
    m <- tryCatch(GetAssayData(obj, assay = "Spatial", layer = "data"),
                  error = function(e) NULL)
    if (is.null(m) || nrow(m) == 0 || ncol(m) == 0) {
      message("  normalising from counts")
      obj <- NormalizeData(obj, normalization.method = "LogNormalize",
                           scale.factor = 1e4, verbose = FALSE)
      m <- GetAssayData(obj, assay = "Spatial", layer = "data")
    }
    rn <- toupper(rownames(m)); keep <- !duplicated(rn)
    m <- m[keep, , drop = FALSE]; rownames(m) <- rn[keep]
    
    feats <- lapply(gene_sets[c(program_names, axis_names)],
                    function(g) intersect(g, rownames(m)))
    feats <- feats[lengths(feats) > 0]
    
    u <- UCell::ScoreSignatures_UCell(m, features = feats,
                                      maxRank = UCELL_MAX_RANK, ncores = 1)
    u <- as.data.frame(u); colnames(u) <- sub("_UCell$", "", colnames(u))
    if (nrow(u) == 0) stop("UCell returned no rows for section ", sec)
    u$key <- paste(sec, rownames(u), sep = "|")
    ucell_list[[sec]] <- as_tibble(u)
  }
  rm(obj); gc(verbose = FALSE)
}

depth <- bind_rows(depth_list)
meta_depth <- meta_all %>% left_join(depth, by = "key")

# =========================================================
# 3. SCHEME SENSITIVITY
# =========================================================
schemes <- c("zmean", "rawmean", "pc1")
scheme_scores <- meta_all %>% select(key, section_id)

for (s in schemes) {
  for (nm in program_names) {
    scheme_scores[[paste0(nm, "__", s)]] <- score_program_scheme(nm, s)
  }
}

if (RUN_UCELL && length(ucell_list) > 0) {
  uc <- bind_rows(ucell_list)
  # rebuild composites from UCell axis scores using the same aggregation
  for (nm in program_names) {
    ax <- PROGRAM_AXES[[nm]]
    if (is.null(ax)) {
      v <- uc[[nm]][match(scheme_scores$key, uc$key)]
    } else {
      ax_use <- intersect(ax, colnames(uc))
      ax_use <- ax_use[sapply(ax_use, function(a) any(is.finite(uc[[a]])))]
      M <- sapply(ax_use, function(a) z_vec(uc[[a]][match(scheme_scores$key, uc$key)]))
      v <- if (length(ax_use) == 0) NA_real_ else rowMeans(M, na.rm = TRUE)
    }
    scheme_scores[[paste0(nm, "__ucell")]] <- v
  }
  schemes <- c(schemes, "ucell")
}

for (s in schemes) {
  for (nm in program_names) {
    scheme_scores[[paste0(nm, "__", s, "_z")]] <-
      z_vec(scheme_scores[[paste0(nm, "__", s)]])
  }
  scheme_scores[[paste0("abT_net__", s)]] <-
    scheme_scores[[paste0("abT_core__", s, "_z")]] -
    scheme_scores[[paste0("inhibitory_program__", s, "_z")]]
  scheme_scores[[paste0("gdT_net__", s)]] <-
    scheme_scores[[paste0("gdT_core__", s, "_z")]] -
    scheme_scores[[paste0("inhibitory_program__", s, "_z")]]
}

write_tsv(scheme_scores, file.path(out_tables, "09b_scheme_scores.tsv"))

# --- depth correlation, per scheme ---
depth_cor <- bind_rows(lapply(schemes, function(s) {
  tibble(scheme = s,
         rho_abT_core = cor(scheme_scores[[paste0("abT_core__", s)]],
                            meta_depth$nCount, method = "spearman", use = "complete.obs"),
         rho_gdT_core = cor(scheme_scores[[paste0("gdT_core__", s)]],
                            meta_depth$nCount, method = "spearman", use = "complete.obs"),
         rho_inhibitory = cor(scheme_scores[[paste0("inhibitory_program__", s)]],
                              meta_depth$nCount, method = "spearman", use = "complete.obs"),
         rho_abT_net = cor(scheme_scores[[paste0("abT_net__", s)]],
                           meta_depth$nCount, method = "spearman", use = "complete.obs"),
         rho_gdT_net = cor(scheme_scores[[paste0("gdT_net__", s)]],
                           meta_depth$nCount, method = "spearman", use = "complete.obs"))
}))

# section-level depth check: does depth explain the Figure 5a positions?
sec_depth <- meta_depth %>%
  group_by(section_id) %>%
  summarise(med_depth = median(nCount, na.rm = TRUE), .groups = "drop") %>%
  left_join(section_summary %>% select(section_id, abT_net_median, gdT_net_median),
            by = "section_id")

depth_section <- tibble(
  level = "section",
  rho_abT_net = cor(sec_depth$med_depth, sec_depth$abT_net_median,
                    method = "spearman", use = "complete.obs"),
  rho_gdT_net = cor(sec_depth$med_depth, sec_depth$gdT_net_median,
                    method = "spearman", use = "complete.obs"),
  n = nrow(sec_depth))

write_tsv(depth_cor, file.path(out_tables, "09b_depth_correlation.tsv"))
write_tsv(depth_section, file.path(out_tables, "09b_depth_correlation_section.tsv"))

message("\n--- Depth coupling (spot level) ---")
print(depth_cor, n = Inf)
message("\n--- Depth coupling (section level) ---")
print(depth_section)

# --- classification agreement between schemes ---
sec_scheme <- scheme_scores %>%
  group_by(section_id) %>%
  summarise(across(where(is.numeric), ~ median(.x, na.rm = TRUE)), .groups = "drop")

for (s in schemes) {
  sec_scheme[[paste0("category__", s)]] <- classify_opportunity(
    sec_scheme[[paste0("abT_net__", s)]],
    sec_scheme[[paste0("gdT_net__", s)]],
    sec_scheme[[paste0("inhibitory_program__", s, "_z")]])
}

cat_agree <- expand_grid(a = schemes, b = schemes) %>%
  rowwise() %>%
  mutate(n_same = sum(sec_scheme[[paste0("category__", a)]] ==
                        sec_scheme[[paste0("category__", b)]], na.rm = TRUE),
         n_total = nrow(sec_scheme),
         pct_agreement = 100 * n_same / n_total) %>%
  ungroup()

write_tsv(cat_agree, file.path(out_tables, "09b_scheme_category_agreement.tsv"))
message("\n--- Classification agreement between aggregation schemes ---")
print(cat_agree %>% select(a, b, n_same, n_total, pct_agreement), n = Inf)

p_depth <- depth_cor %>%
  pivot_longer(starts_with("rho_"), names_to = "score", values_to = "rho") %>%
  mutate(score = sub("^rho_", "", score)) %>%
  ggplot(aes(x = score, y = rho, fill = scheme)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey70") +
  geom_col(position = "dodge", width = 0.75) +
  labs(x = NULL, y = "Spearman rho vs nCount",
       title = "Sequencing-depth coupling by scoring scheme",
       subtitle = "Uniform mild coupling is tolerable; opposite signs between modalities are not") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        plot.title = element_text(face = "bold"))

p_agree <- ggplot(cat_agree, aes(x = a, y = b, fill = pct_agreement)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.0f%%\n(%d/%d)", pct_agreement, n_same, n_total)),
            colour = "white", size = 3, fontface = "bold") +
  scale_fill_gradient(low = "#9ecae1", high = "#08519c", name = "% agreement") +
  labs(x = NULL, y = NULL, title = "Classification agreement between schemes") +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(), plot.title = element_text(face = "bold"))

p_qc2 <- (p_depth | p_agree) + plot_annotation(tag_levels = "a")
ggsave(file.path(out_plots, "S_QC2_depth_and_schemes.png"), p_qc2,
       width = 13, height = 5.5, dpi = 600, bg = "white")
ggsave(file.path(out_plots, "S_QC2_depth_and_schemes.pdf"), p_qc2,
       width = 13, height = 5.5, bg = "white", device = cairo_pdf)

# =========================================================
# 4. AXIS INDEPENDENCE
#
# Expected to be low. The composites are formative indices: a tumor may be
# targetable through butyrophilins OR NKG2D ligands OR ephrin independently, so
# the axes define the construct rather than reflecting one latent cause. Low
# correlation confirms the design; it is not a coherence failure, and Cronbach's
# alpha is not an appropriate statistic here.
# =========================================================
axis_scores <- sapply(axis_names, function(a) {
  g <- gene_sets[[a]]
  if (length(g) == 0) return(rep(NA_real_, nrow(expr_prog)))
  z_vec(score_set_scheme(expr_prog[, g, drop = FALSE], "zmean"))
})
axis_scores <- axis_scores[, colSums(is.finite(axis_scores)) > 0, drop = FALSE]

gd_ax  <- intersect(PROGRAM_AXES$gdT_core, colnames(axis_scores))
inh_ax <- intersect(PROGRAM_AXES$inhibitory_program, colnames(axis_scores))

axis_ind <- bind_rows(
  if (length(gd_ax)  >= 2) cor_long(axis_scores[, gd_ax,  drop = FALSE], "gdT axes"),
  if (length(inh_ax) >= 2) cor_long(axis_scores[, inh_ax, drop = FALSE], "inhibitory axes")
)
write_tsv(axis_ind, file.path(out_tables, "09b_axis_independence.tsv"))

message("\n--- Axis independence (spot level, off-diagonal) ---")
print(axis_ind %>% filter(var_a != var_b) %>% arrange(desc(abs(rho))), n = Inf)

p_axind <- axis_ind %>%
  ggplot(aes(x = var_a, y = var_b, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 2.8) +
  facet_wrap(~group, scales = "free") +
  scale_fill_gradient2(low = "#20B3B7", mid = "white", high = "#F07167",
                       midpoint = 0, limits = c(-1, 1), name = "Spearman") +
  labs(x = NULL, y = NULL,
       title = "Mechanism axes are independent by design",
       subtitle = "Low correlation is expected for a formative index and supports averaging axes equally") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(), plot.title = element_text(face = "bold"))

ggsave(file.path(out_plots, "S_QC3_axis_independence.png"), p_axind,
       width = 12, height = 5.5, dpi = 600, bg = "white")
ggsave(file.path(out_plots, "S_QC3_axis_independence.pdf"), p_axind,
       width = 12, height = 5.5, bg = "white", device = cairo_pdf)

# =========================================================
# 5. SPLIT-HALF RELIABILITY
#
# Spots within each section are split at random; section medians are computed
# from each half and compared across sections. This asks whether a section HAS
# an opportunity score, independently of what that score means.
# =========================================================
net_ab <- scheme_scores[[paste0("abT_net__", cfg$scoring_method)]]
net_gd <- scheme_scores[[paste0("gdT_net__", cfg$scoring_method)]]
inh_z  <- scheme_scores[[paste0("inhibitory_program__", cfg$scoring_method, "_z")]]

split_df <- tibble(section_id = scheme_scores$section_id,
                   ab = net_ab, gd = net_gd, inh = inh_z)

split_rows <- vector("list", N_SPLIT_REP)
for (r in seq_len(N_SPLIT_REP)) {
  d <- split_df %>%
    group_by(section_id) %>%
    mutate(half = sample(rep(c("A", "B"), length.out = n()))) %>%
    group_by(section_id, half) %>%
    summarise(ab = median(ab, na.rm = TRUE),
              gd = median(gd, na.rm = TRUE),
              inh = median(inh, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = half, values_from = c(ab, gd, inh))

  cat_a <- classify_opportunity(d$ab_A, d$gd_A, d$inh_A)
  cat_b <- classify_opportunity(d$ab_B, d$gd_B, d$inh_B)

  split_rows[[r]] <- tibble(
    rep = r,
    rho_ab = cor(d$ab_A, d$ab_B, method = "spearman", use = "complete.obs"),
    rho_gd = cor(d$gd_A, d$gd_B, method = "spearman", use = "complete.obs"),
    pct_category_same = 100 * mean(as.character(cat_a) == as.character(cat_b), na.rm = TRUE))
}
split_half <- bind_rows(split_rows)

write_tsv(split_half, file.path(out_tables, "09b_split_half.tsv"))

split_summary <- split_half %>%
  summarise(across(c(rho_ab, rho_gd, pct_category_same),
                   list(mean = ~mean(.x, na.rm = TRUE),
                        q05  = ~quantile(.x, 0.05, na.rm = TRUE),
                        q95  = ~quantile(.x, 0.95, na.rm = TRUE))))
write_tsv(split_summary, file.path(out_tables, "09b_split_half_summary.tsv"))

message("\n--- Split-half reliability (", N_SPLIT_REP, " repetitions) ---")
print(split_summary)
message("High rho and high category agreement mean the section-level score is a ",
        "stable property of the tissue rather than sampling noise.")

p_split <- split_half %>%
  pivot_longer(c(rho_ab, rho_gd, pct_category_same),
               names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric,
                         rho_ab = paste0("rho, ", AXIS_LAB_AB),
                         rho_gd = paste0("rho, ", AXIS_LAB_GD),
                         pct_category_same = "% category agreement")) %>%
  ggplot(aes(x = value)) +
  geom_histogram(bins = 25, fill = "#9ecae1", colour = "white", linewidth = 0.2) +
  facet_wrap(~metric, scales = "free") +
  labs(x = NULL, y = "Repetitions",
       title = "Split-half reliability of section-level scores",
       subtitle = paste0("Spots split at random within each section, ",
                         N_SPLIT_REP, " repetitions")) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(out_plots, "S_QC4_split_half.png"), p_split,
       width = 11, height = 4.2, dpi = 600, bg = "white")
ggsave(file.path(out_plots, "S_QC4_split_half.pdf"), p_split,
       width = 11, height = 4.2, bg = "white", device = cairo_pdf)

writeLines(c(
  paste0("scoring_method_primary: ", cfg$scoring_method),
  paste0("scoring_aggregation: ", cfg$scoring_aggregation),
  paste0("detection_min_pct: ", cfg$detection_min_pct),
  paste0("schemes_compared: ", paste(schemes, collapse = ",")),
  paste0("ucell_run: ", RUN_UCELL),
  paste0("split_half_repetitions: ", N_SPLIT_REP),
  paste0("seed: ", SEED),
  paste0("run_date: ", as.character(Sys.Date()))
), file.path(out_tables, "09b_run_parameters.txt"))

message("\nDone. Outputs written to: ", out_root)
