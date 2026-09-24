# =============================================================================
# 10_zone_and_spatial.R
#
# Purpose
# -------
# Zone-level and spatial characterisation of the opportunity scores:
# Figure 5c (compartment components) and Figure 5d (spatial maps), plus
# depth-residualised, patient-independence and per-axis sensitivity analyses.
#
# Reads the scores written by 09 and never recomputes them, so the two scripts
# cannot drift apart. Requires the hires tissue images to place spots.
#
# Inputs
# ------
#   stage_input("scoring", "rds/09_spot_scores.rds")
#   stage_input("scoring", "rds/09_section_summary.rds")
#   stage_dir("zones_all")   tables/zones_all_sections.tsv
#   stage_dir("raw_vis")     vis_section_<sec>_raw.rds
#   RAW_ROOT                 tissue_hires_image.png|jpg per section
#
# Outputs
# -------
#   stage_dir("zone_spatial")
#     tables/10_depth_coupling.tsv, 10_zone_summaries.tsv, 10_zone_tests.tsv,
#     tables/10_zone_tests_patient_sensitivity.tsv,
#     tables/10_image_lookup.tsv, 10_spot_scores_with_coordinates.tsv,
#     tables/10_section_balance_summary.tsv
#     rds/10_spot_scores_with_coordinates.rds
#     plots/Fig5C_compartment_components.{png,pdf}
#     plots/S_compartment_depth_residualised.{png,pdf}
#     plots/S_zone_gdT_axes.{png,pdf}, S_zone_inhibitory_axes.{png,pdf}
#     plots/S_spatial_balance_all_sections.{png,pdf}
#
# Stochastic
# ----------
#   none -- no sampling or permutation in this stage.
#
# Runtime
# -------
#   ~20 minutes
# =============================================================================

source(here::here("config.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})

# =========================================================
# 10 — Compartment comparison and spatial opportunity maps
#
# Figure 5c (component scores by tumor compartment) and Figure 5d (spatial
# gammadelta-minus-alphabeta balance), plus supplementary panels for the
# gammadelta recognition modalities and the immunosuppressive mechanisms.
#
# Reads script 09's output. It never redefines a gene set or recomputes a
# score, so it cannot drift from 09.
#
# WHY 5c USES COMPONENT SCORES, NOT NET SCORES
# --------------------------------------------
# Net scores (core - inhibitory) are appropriate for ranking whole sections
# (Figure 5a), but not for comparing compartments within a section. The
# immunosuppressive program carries a strong mesenchymal gradient, so
# subtracting it from a compartment-neutral gammadelta core manufactures an
# epithelial-high net: gdT_net separates the compartments in 14/14 sections
# while gdT_core alone shows no difference (6/14 after depth control). Panel 5c
# therefore reports the three components separately.
#
# GENE OVERLAP WITH THE ZONE PROGRAMS
# -----------------------------------
# Zone programs come from published snRNA-seq cluster markers; opportunity
# programs from literature-defined recognition and suppression genes. They were
# assembled independently, but four genes appear in both: ICAM1 and IRF1 (in the
# cores) and PTGS2 and VEGFA (in the inhibitory program). The mesenchymal zone
# program also shares 91 genes with the myeloid program (Jaccard 0.229),
# reflecting the inflammatory character of that compartment. This is expected -
# the mesenchymal compartment IS the inflammatory one - but it means the
# association between immunosuppressive tone and mesenchymal identity is partly
# definitional and should not be presented as an independent observation.
# The per-axis supplementary panel addresses this directly: HLA-E, TGF-beta,
# adenosine and galectin-9 contain no zone-program genes, so if they track
# mesenchymal too, the overlap is not what is driving the result.
#
# STATISTICS
# ----------
# Section-level summaries with a paired Wilcoxon signed-rank test across
# sections. Spots within a section are not independent, so a spot-level test
# would treat thousands of correlated observations as replicates. n = sections.
# The count of sections favouring each compartment is reported alongside the
# p-value; at n = 14 it is the more interpretable statistic.
#
# DEPTH CONTROL
# -------------
# Each score is also residualised on log sequencing depth WITHIN section and the
# tests repeated. NOTE: this works for the composites and for well-detected
# axes, but overcorrects sparsely detected ones (a linear fit cannot represent a
# zero-inflated score, so at high depth the still-zero spots receive large
# negative residuals). Residual depth coupling is reported per score so the
# affected axes can be identified; report raw values for those.
#
# BALANCE (5d)
# ------------
# balance = gdT_net - abT_net = gdT_core_z - abT_core_z. The inhibitory term
# cancels exactly, so this contrast is unaffected by the issue above.
#
# SUMMARY STATISTIC
# -----------------
# Components use the median. Mechanism axes use the MEAN: several retain only
# one sparsely detected gene, and for those the section median collapses to the
# value of zero on the z scale.
#
# OUTPUTS
# -------
#   tables/10_zone_summaries.tsv
#   tables/10_zone_tests.tsv
#   tables/10_zone_tests_patient_sensitivity.tsv
#   tables/10_depth_coupling.tsv
#   tables/10_image_lookup.tsv
#   tables/10_spot_scores_with_coordinates.tsv | rds
#   tables/10_section_balance_summary.tsv
#   plots/Fig5C_compartment_components.png | .pdf
#   plots/Fig5D_spatial_balance.png | .pdf
#   plots/S_compartment_depth_residualised.png | .pdf
#   plots/S_zone_gdT_axes.png | .pdf
#   plots/S_zone_inhibitory_axes.png | .pdf
#   plots/S_spatial_balance_all_sections.png | .pdf
# =========================================================


# -------------------------
# CONFIG
# -------------------------
SPATIAL_ROOT <- RAW_ROOT

phase06_tables <- file.path(stage_dir("zones_all"), "tables")
phase09_rds    <- file.path(stage_dir("scoring"), "rds")
raw_vis_dir    <- stage_dir("raw_vis")

out_root   <- stage_dir("zone_spatial")
out_tables <- ensure_dir(file.path(out_root, "tables"))
out_plots  <- ensure_dir(file.path(out_root, "plots"))
out_rds    <- ensure_dir(file.path(out_root, "rds"))

zones_keep <- c("Epithelial", "Mesenchymal")
main_figure_sections <- c("459", "459_2")

RESIDUALISE_DEPTH <- TRUE
IMAGE_FILENAMES   <- c("tissue_hires_image.png", "tissue_hires_image.jpg")

COMPONENTS <- c("abT_core_z", "gdT_core_z", "inhibitory_program_z")
NETS       <- c("abT_net", "gdT_net")

GD_AXES  <- c("gdT_btn", "gdT_mevalonate", "gdT_nkg2d",
              "gdT_ephrin", "gdT_adhesion_dnam", "gdT_ipp_drain")
INH_AXES <- c("inh_mhc_checkpoint", "inh_classical_ckpt", "inh_tgfb",
              "inh_adenosine", "inh_pge2", "inh_angiogenic")

# Genes shared with the mesenchymal zone program, flagged in the axis panel
ZONE_OVERLAP_AXES <- c("inh_pge2", "inh_angiogenic", "gdT_adhesion_dnam")

# Axis labels name what each axis actually contains after detection filtering,
# so a single-gene axis is not presented as a multi-gene mechanism.
SCORE_LABELS <- c(
  abT_core_z           = "\u03b1\u03b2T ligand availability",
  gdT_core_z           = "\u03b3\u03b4T ligand availability",
  inhibitory_program_z = "Immunosuppressive tone",
  abT_net              = "\u03b1\u03b2T net",
  gdT_net              = "\u03b3\u03b4T net",
  gdT_btn              = "Butyrophilin (BTN2A1/3A1)",
  gdT_mevalonate       = "Mevalonate / IPP synthesis",
  gdT_nkg2d            = "NKG2D ligand (MICA)",
  gdT_ephrin           = "EPHA2",
  gdT_adhesion_dnam    = "DNAM-1 / adhesion",
  gdT_ipp_drain        = "FDPS (IPP consumption)",
  inh_mhc_checkpoint   = "HLA-E / NKG2A",
  inh_classical_ckpt   = "Galectin-9 (LGALS9)",
  inh_tgfb             = "TGF-\u03b2 (TGFB1/2/3)",
  inh_adenosine        = "Adenosine (NT5E/ENTPD1)",
  inh_pge2             = "COX-2 (PTGS2)",
  inh_angiogenic       = "VEGF (VEGFA)"
)

q_clip     <- 0.02
point_size <- 1.35
he_alpha   <- 0.42
dim_tol    <- 2

AXIS_LAB_AB <- "\u03b1\u03b2T net"
AXIS_LAB_GD <- "\u03b3\u03b4T net"
BAL_LAB     <- "\u03b3\u03b4T \u2212 \u03b1\u03b2T"

ZONE_COLORS <- c(Epithelial = "#00BF7D", Mesenchymal = "#F07167")

# -------------------------
# HELPERS
# -------------------------
score_label <- function(x) ifelse(x %in% names(SCORE_LABELS), SCORE_LABELS[x], x)

extract_barcode16 <- function(x) {
  x <- toupper(as.character(x))
  m <- regexpr("[ACGT]{16}", x)
  out <- rep(NA_character_, length(x)); ok <- m > 0
  out[ok] <- regmatches(x, m); out
}

read_image_any <- function(img_path) {
  ext <- tolower(tools::file_ext(img_path))
  if (ext == "png") {
    if (!requireNamespace("png", quietly = TRUE)) stop("Package 'png' required")
    return(png::readPNG(img_path))
  }
  if (ext %in% c("jpg", "jpeg")) {
    if (!requireNamespace("jpeg", quietly = TRUE)) stop("Package 'jpeg' required")
    return(jpeg::readJPEG(img_path))
  }
  stop("Unsupported image type: ", img_path)
}

sym_limits <- function(x, q = q_clip) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(c(-1, 1))
  lim <- max(abs(quantile(x, c(q, 1 - q), na.rm = TRUE)))
  if (!is.finite(lim) || lim == 0) lim <- max(abs(x), na.rm = TRUE)
  if (!is.finite(lim) || lim == 0) lim <- 1
  c(-lim, lim)
}

first_existing <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) NA_character_ else hit[1]
}

paired_zone_test <- function(df, col) {
  w <- df %>%
    select(section_id, zone_call, value = all_of(col)) %>%
    pivot_wider(names_from = zone_call, values_from = value) %>%
    filter(is.finite(Epithelial), is.finite(Mesenchymal))
  if (nrow(w) < 3) {
    return(tibble(score = col, n_sections = nrow(w), median_diff = NA_real_,
                  n_higher_epithelial = NA_integer_, p_value = NA_real_))
  }
  d <- w$Epithelial - w$Mesenchymal
  tt <- suppressWarnings(wilcox.test(w$Epithelial, w$Mesenchymal,
                                     paired = TRUE, exact = FALSE))
  tibble(score = col, n_sections = nrow(w),
         median_diff = median(d, na.rm = TRUE),
         n_higher_epithelial = sum(d > 0, na.rm = TRUE),
         p_value = tt$p.value)
}

resid_on_depth <- function(y, logd) {
  ok <- is.finite(y) & is.finite(logd)
  out <- rep(NA_real_, length(y))
  if (sum(ok) < 10 || length(unique(logd[ok])) < 3) return(out)
  fit <- stats::lm(y[ok] ~ logd[ok])
  out[ok] <- as.numeric(stats::residuals(fit))
  out
}

# -------------------------
# LOAD
# -------------------------
spot_rds <- file.path(phase09_rds, "09_spot_scores.rds")
if (!file.exists(spot_rds)) stop("Missing script 09 output: ", spot_rds)
spot <- readRDS(spot_rds)

section_summary <- readRDS(file.path(phase09_rds, "09_section_summary.rds"))
zones_all <- read_tsv(file.path(phase06_tables, "zones_all_sections.tsv"),
                      show_col_types = FALSE)

req <- c("section_id", "cell", "zone_call", COMPONENTS, NETS)
missing <- setdiff(req, colnames(spot))
if (length(missing) > 0) stop("Missing columns in 09 output: ", paste(missing, collapse = ", "))

spot <- spot %>% mutate(balance = gdT_net - abT_net)

GD_AXES  <- intersect(GD_AXES,  colnames(spot))
INH_AXES <- intersect(INH_AXES, colnames(spot))
ALL_AXES <- c(GD_AXES, INH_AXES)
message("gdT axes available:        ", paste(GD_AXES, collapse = ", "))
message("Inhibitory axes available: ", paste(INH_AXES, collapse = ", "))

section_labels <- section_summary %>%
  select(section_id, section_label, sample_status, opportunity_category)

# -------------------------
# DEPTH AND RESIDUALS
# -------------------------
if (RESIDUALISE_DEPTH) {
  message("\nLoading per-spot sequencing depth...")
  depth <- bind_rows(lapply(unique(spot$section_id), function(s) {
    f <- file.path(raw_vis_dir, paste0("vis_section_", s, "_raw.rds"))
    if (!file.exists(f)) { warning("Missing vis file: ", s); return(NULL) }
    o <- readRDS(f)
    tibble(key = paste(s, colnames(o), sep = "|"), nCount = o$nCount_Spatial)
  }))
  
  spot <- spot %>% left_join(depth, by = "key")
  message("Spots with depth: ", sum(is.finite(spot$nCount)), "/", nrow(spot))
  
  resid_cols <- c(COMPONENTS, NETS, ALL_AXES)
  spot <- spot %>%
    mutate(log_depth = log1p(nCount)) %>%
    group_by(section_id) %>%
    mutate(across(all_of(resid_cols), ~ resid_on_depth(.x, log_depth),
                  .names = "{.col}__resid")) %>%
    ungroup()
  
  depth_tbl <- bind_rows(lapply(resid_cols, function(c) {
    tibble(score = c,
           label = score_label(c),
           rho_spot_raw = suppressWarnings(cor(spot$nCount, spot[[c]],
                                               method = "spearman", use = "complete.obs")),
           rho_spot_resid = suppressWarnings(cor(spot$nCount, spot[[paste0(c, "__resid")]],
                                                 method = "spearman", use = "complete.obs")))
  })) %>%
    mutate(residual_ok = abs(rho_spot_resid) < 0.10)
  
  write_tsv(depth_tbl, file.path(out_tables, "10_depth_coupling.tsv"))
  message("\n--- Depth coupling, spot level ---")
  print(depth_tbl %>% select(score, rho_spot_raw, rho_spot_resid, residual_ok), n = Inf)
  if (any(!depth_tbl$residual_ok)) {
    message("Residualisation overcorrected (|rho| >= 0.10) for: ",
            paste(depth_tbl$score[!depth_tbl$residual_ok], collapse = ", "),
            "\n  These are zero-inflated scores a linear fit cannot represent. ",
            "Report their RAW values.")
  }
}

# =========================================================
# SECTION-LEVEL SUMMARIES BY COMPARTMENT
# =========================================================
med_cols <- c(COMPONENTS, NETS, "balance")
if (RESIDUALISE_DEPTH) med_cols <- c(med_cols, paste0(c(COMPONENTS, NETS), "__resid"))

mean_cols <- ALL_AXES
if (RESIDUALISE_DEPTH) mean_cols <- c(mean_cols, paste0(ALL_AXES, "__resid"))

zone_med <- spot %>%
  filter(zone_call %in% zones_keep) %>%
  group_by(section_id, zone_call) %>%
  summarise(n_spots = n(),
            across(all_of(med_cols), ~ median(.x, na.rm = TRUE)),
            across(all_of(mean_cols), ~ mean(.x, na.rm = TRUE)),
            frac_gd_favoured = mean(balance > 0, na.rm = TRUE),
            .groups = "drop") %>%
  left_join(section_labels, by = "section_id")

write_tsv(zone_med, file.path(out_tables, "10_zone_summaries.tsv"))

run_tests <- function(cols, level, version, df = zone_med) {
  bind_rows(lapply(cols, function(c) paired_zone_test(df, c))) %>%
    mutate(
      level = level,
      version = version,
      score_base = sub("__resid$", "", score),
      p_adj = p.adjust(p_value, method = "BH")
    )
}

zone_tests <- bind_rows(
  run_tests(COMPONENTS, "component", "raw"),
  run_tests(NETS, "net", "raw"),
  run_tests(GD_AXES, "gdT axis", "raw"),
  run_tests(INH_AXES, "inhibitory axis", "raw"),
  if (RESIDUALISE_DEPTH) run_tests(paste0(COMPONENTS, "__resid"), "component", "depth-residualised"),
  if (RESIDUALISE_DEPTH) run_tests(paste0(NETS, "__resid"), "net", "depth-residualised"),
  if (RESIDUALISE_DEPTH) run_tests(paste0(GD_AXES, "__resid"), "gdT axis", "depth-residualised"),
  if (RESIDUALISE_DEPTH) run_tests(paste0(INH_AXES, "__resid"), "inhibitory axis", "depth-residualised")
) %>%
  mutate(label = score_label(score_base),
         zone_program_overlap = score_base %in% ZONE_OVERLAP_AXES) %>%
  arrange(level, score_base, version)

write_tsv(zone_tests, file.path(out_tables, "10_zone_tests.tsv"))

message("\n--- Compartment tests (paired Wilcoxon on section summaries) ---")
print(zone_tests %>%
        select(level, label, version, n_higher_epithelial, n_sections,
               median_diff, p_value, p_adj, zone_program_overlap), n = Inf)

# =========================================================
# PATIENT-INDEPENDENCE SENSITIVITY
#
# The 14 sections represent 11 patients because three patients have matched
# primary-relapse sections. Two raw-score sensitivity analyses are therefore
# performed:
#   1) one section per patient (n = 11), retaining the primary section;
#   2) leave-one-relapse-out (n = 13), repeated for each matched relapse.
#
# These analyses assess whether the section-level compartment comparisons are
# dependent on the repeated sampling of the three matched patients.
# =========================================================
message("\nRunning patient-independence sensitivity analysis...")

relapse_sections <- vapply(
  MATCHED_PAIRS,
  `[[`,
  character(1),
  "relapse"
)

sens_cols <- list(
  component = COMPONENTS,
  net = NETS,
  `gdT axis` = GD_AXES,
  `inhibitory axis` = INH_AXES
)

run_tests_labelled <- function(zone_med_sub, version_tag) {
  bind_rows(lapply(names(sens_cols), function(lvl) {
    run_tests(
      sens_cols[[lvl]],
      level = lvl,
      version = version_tag,
      df = zone_med_sub
    )
  })) %>%
    mutate(
      label = score_label(score_base),
      zone_program_overlap = score_base %in% ZONE_OVERLAP_AXES
    )
}

# 1. One section per patient: retain primary sections, drop matched relapses.
zone_med_one_per_patient <- zone_med %>%
  filter(!section_id %in% relapse_sections)

sens_one_per_patient <- run_tests_labelled(
  zone_med_one_per_patient,
  "one-per-patient"
)

# 2. Leave one relapse section out at a time.
sens_leave_one_relapse_out <- bind_rows(lapply(relapse_sections, function(sec) {
  run_tests_labelled(
    zone_med %>% filter(section_id != sec),
    paste0("excl_", sec)
  )
}))

zone_tests_sensitivity <- bind_rows(
  sens_one_per_patient,
  sens_leave_one_relapse_out
) %>%
  arrange(level, score_base, version)

write_tsv(
  zone_tests_sensitivity,
  file.path(out_tables, "10_zone_tests_patient_sensitivity.tsv")
)

message("\n--- One section per patient (n = 11): components ---")
print(
  sens_one_per_patient %>%
    filter(level == "component") %>%
    select(
      label,
      n_sections,
      n_higher_epithelial,
      median_diff,
      p_value,
      p_adj
    ),
  n = Inf
)

# Directly compare the main component calls with the one-per-patient analysis.
main_components <- zone_tests %>%
  filter(level == "component", version == "raw") %>%
  select(
    score_base,
    label,
    median_diff_main = median_diff,
    p_adj_main = p_adj
  )

sensitivity_components <- sens_one_per_patient %>%
  filter(level == "component") %>%
  select(
    score_base,
    median_diff_sensitivity = median_diff,
    p_adj_sensitivity = p_adj
  )

component_sensitivity_comparison <- main_components %>%
  left_join(sensitivity_components, by = "score_base") %>%
  mutate(
    direction_changed =
      sign(median_diff_main) != sign(median_diff_sensitivity),
    significance_changed =
      (p_adj_main < 0.05) != (p_adj_sensitivity < 0.05)
  )

message("\n--- Components: n = 14 vs n = 11 ---")
print(component_sensitivity_comparison, n = Inf)

if (any(component_sensitivity_comparison$direction_changed, na.rm = TRUE)) {
  warning(
    "Direction changed for at least one opportunity component in the ",
    "one-section-per-patient sensitivity analysis.",
    call. = FALSE
  )
}

if (any(component_sensitivity_comparison$significance_changed, na.rm = TRUE)) {
  warning(
    "BH-adjusted significance changed for at least one opportunity component ",
    "in the one-section-per-patient sensitivity analysis.",
    call. = FALSE
  )
}

# =========================================================
# PANEL BUILDER
# =========================================================
make_compartment_panel <- function(cols, suffix = "", title_txt, subtitle_txt,
                                   ncol = NULL) {
  use <- if (nzchar(suffix)) paste0(cols, suffix) else cols
  use <- intersect(use, colnames(zone_med))
  if (length(use) == 0) return(NULL)
  
  long <- zone_med %>%
    select(section_id, zone_call, all_of(use)) %>%
    pivot_longer(all_of(use), names_to = "score", values_to = "value") %>%
    mutate(score_base = sub("__resid$", "", score),
           score_lab = factor(score_label(score_base), levels = score_label(cols)),
           zone_call = factor(zone_call, levels = zones_keep)) %>%
    filter(is.finite(value))
  
  ann <- zone_tests %>%
    filter(score %in% use) %>%
    mutate(score_lab = factor(score_label(score_base), levels = score_label(cols)),
           label_txt = paste0("p = ", format.pval(p_value, digits = 2, eps = 1e-4),
                              "\nepithelial higher: ", n_higher_epithelial,
                              "/", n_sections))
  
  ggplot(long, aes(x = zone_call, y = value)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey70") +
    geom_boxplot(aes(fill = zone_call), width = 0.55,
                 outlier.shape = NA, alpha = 0.7) +
    geom_line(aes(group = section_id), colour = "grey55", linewidth = 0.35) +
    geom_point(size = 1.7, colour = "grey20") +
    facet_wrap(~score_lab, scales = "free_y", ncol = ncol) +
    geom_text(data = ann, aes(x = 1.5, y = Inf, label = label_txt),
              inherit.aes = FALSE, vjust = 1.15, size = 2.8) +
    scale_fill_manual(values = ZONE_COLORS, guide = "none") +
    labs(x = NULL, y = "Section summary score",
         title = title_txt, subtitle = subtitle_txt) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(size = 8))
}

# =========================================================
# FIGURE 5C — components
# =========================================================
p_5c <- make_compartment_panel(
  COMPONENTS, "",
  "Opportunity components by tumor compartment",
  paste0("Components are shown separately rather than as net scores: the ",
         "immunosuppressive term carries a mesenchymal gradient, so subtracting ",
         "it would impose a compartment difference on the \u03b3\u03b4T score.\n",
         "Zone and opportunity programs were assembled independently; four genes ",
         "overlap (ICAM1, IRF1, PTGS2, VEGFA \u2014 see Methods). Each line is one section."))

if (!is.null(p_5c)) {
  ggsave(file.path(out_plots, "Fig5C_compartment_components.png"), p_5c,
         width = 10.5, height = 5.6, dpi = 600, bg = "white")
  ggsave(file.path(out_plots, "Fig5C_compartment_components.pdf"), p_5c,
         width = 10.5, height = 5.6, bg = "white", device = cairo_pdf)
}

if (RESIDUALISE_DEPTH) {
  p_5c_res <- make_compartment_panel(
    COMPONENTS, "__resid",
    "Opportunity components by compartment, depth-residualised",
    paste0("Each score regressed on log sequencing depth within section before ",
           "summarising. Between-section depth differences are retained ",
           "deliberately, being confounded with patient."))
  if (!is.null(p_5c_res)) {
    ggsave(file.path(out_plots, "S_compartment_depth_residualised.png"), p_5c_res,
           width = 10.5, height = 5.6, dpi = 600, bg = "white")
    ggsave(file.path(out_plots, "S_compartment_depth_residualised.pdf"), p_5c_res,
           width = 10.5, height = 5.6, bg = "white", device = cairo_pdf)
  }
}

# =========================================================
# SUPPLEMENTARY — gdT recognition modalities
# =========================================================
p_gd <- make_compartment_panel(
  GD_AXES, "",
  "\u03b3\u03b4T recognition modalities by tumor compartment",
  paste0("Axis means; y scales differ between panels. No score is subtracted ",
         "from another. FDPS is a modifier and is not part of the composite."),
  ncol = 3)

if (!is.null(p_gd)) {
  ggsave(file.path(out_plots, "S_zone_gdT_axes.png"), p_gd,
         width = 11, height = 7, dpi = 600, bg = "white")
  ggsave(file.path(out_plots, "S_zone_gdT_axes.pdf"), p_gd,
         width = 11, height = 7, bg = "white", device = cairo_pdf)
}

# =========================================================
# SUPPLEMENTARY — immunosuppressive mechanisms
#
# The point of this panel: VEGF and COX-2 are the two axes whose genes also
# appear in the mesenchymal zone program. HLA-E, galectin-9, TGF-beta and
# adenosine do not. If those four track mesenchymal as well, the compartment
# association is not an artifact of the shared genes.
# =========================================================
p_inh <- make_compartment_panel(
  INH_AXES, "",
  "Immunosuppressive mechanisms by tumor compartment",
  paste0("Axis means; y scales differ between panels. VEGF (VEGFA) and COX-2 ",
         "(PTGS2) also appear in the mesenchymal zone program; HLA-E, ",
         "galectin-9, TGF-\u03b2 and adenosine do not."),
  ncol = 3)

if (!is.null(p_inh)) {
  ggsave(file.path(out_plots, "S_zone_inhibitory_axes.png"), p_inh,
         width = 11, height = 7, dpi = 600, bg = "white")
  ggsave(file.path(out_plots, "S_zone_inhibitory_axes.pdf"), p_inh,
         width = 11, height = 7, bg = "white", device = cairo_pdf)
}

# quick console read-out of the overlap question
inh_raw <- zone_tests %>% filter(level == "inhibitory axis", version == "raw")
if (nrow(inh_raw) > 0) {
  n_mes <- sum(inh_raw$n_higher_epithelial < inh_raw$n_sections / 2, na.rm = TRUE)
  message("\nInhibitory axes leaning mesenchymal: ", n_mes, "/", nrow(inh_raw))
  message("Of the axes with NO zone-program gene overlap: ",
          sum(!inh_raw$zone_program_overlap &
                inh_raw$n_higher_epithelial < inh_raw$n_sections / 2, na.rm = TRUE),
          "/", sum(!inh_raw$zone_program_overlap))
}

# =========================================================
# ATTACH COORDINATES
# =========================================================
coord_cols <- c("section_id", "cell", "in_tissue",
                "pxl_col_in_hires", "pxl_row_in_hires",
                "img_width_hires", "img_height_hires")

zones_join <- zones_all %>%
  mutate(barcode16_join = coalesce(
    if ("barcode16_join" %in% names(.)) as.character(barcode16_join) else NA_character_,
    if ("barcode16" %in% names(.)) as.character(barcode16) else NA_character_,
    extract_barcode16(cell))) %>%
  select(any_of(c(coord_cols, "barcode16_join"))) %>%
  distinct(section_id, cell, .keep_all = TRUE)

have_dims <- all(c("img_width_hires", "img_height_hires") %in% colnames(zones_join))
if (!have_dims) warning("Recorded hires dimensions not found - image/coordinate ",
                        "agreement cannot be verified.")

join_cols <- intersect(coord_cols, colnames(zones_join))

spot_xy <- spot %>%
  left_join(zones_join %>% select(all_of(join_cols)), by = c("section_id", "cell"))

if (any(is.na(spot_xy$pxl_col_in_hires))) {
  message("Rescuing coordinates via barcode for ",
          sum(is.na(spot_xy$pxl_col_in_hires)), " spots")
  rescue <- spot %>%
    mutate(barcode16_join = extract_barcode16(cell)) %>%
    select(section_id, cell, barcode16_join) %>%
    left_join(zones_join %>%
                select(any_of(c("section_id", "barcode16_join",
                                setdiff(join_cols, c("section_id", "cell"))))) %>%
                distinct(),
              by = c("section_id", "barcode16_join"))
  spot_xy <- spot_xy %>%
    select(-any_of(setdiff(join_cols, c("section_id", "cell")))) %>%
    left_join(rescue %>% select(-barcode16_join), by = c("section_id", "cell"))
}

spot_xy <- spot_xy %>% left_join(section_labels, by = "section_id")

# -------------------------
# IMAGE PATHS
# -------------------------
img_lookup <- tibble(section_id = sort(unique(spot_xy$section_id))) %>%
  rowwise() %>%
  mutate(img_file = first_existing(file.path(SPATIAL_ROOT, section_id,
                                             "spatial", IMAGE_FILENAMES))) %>%
  ungroup()

if (any(is.na(img_lookup$img_file))) {
  message("Direct path missed ", sum(is.na(img_lookup$img_file)),
          " section(s); searching ", SPATIAL_ROOT)
  found <- list.files(SPATIAL_ROOT, pattern = "tissue_hires_image\\.(png|jpe?g)$",
                      recursive = TRUE, full.names = TRUE)
  if (length(found) > 0) {
    found_tbl <- tibble(path = found, folder = basename(dirname(dirname(found))))
    img_lookup <- img_lookup %>%
      left_join(found_tbl %>% distinct(folder, .keep_all = TRUE),
                by = c("section_id" = "folder")) %>%
      mutate(img_file = coalesce(img_file, path)) %>% select(-path)
  }
}

img_lookup <- img_lookup %>% mutate(found = !is.na(img_file))
write_tsv(img_lookup, file.path(out_tables, "10_image_lookup.tsv"))
message("\n--- Image lookup ---"); print(img_lookup, n = Inf)
if (any(!img_lookup$found)) {
  warning("No image found for: ",
          paste(img_lookup$section_id[!img_lookup$found], collapse = ", "),
          " - plotted without the H&E layer.")
}

spot_xy <- spot_xy %>% left_join(img_lookup %>% select(section_id, img_file),
                                 by = "section_id")

write_tsv(spot_xy, file.path(out_tables, "10_spot_scores_with_coordinates.tsv"))
saveRDS(spot_xy,   file.path(out_rds, "10_spot_scores_with_coordinates.rds"))

balance_summary <- spot_xy %>%
  group_by(section_id, section_label, opportunity_category) %>%
  summarise(n_spots = n(),
            balance_median = median(balance, na.rm = TRUE),
            frac_gd_favoured = mean(balance > 0, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(desc(balance_median))

write_tsv(balance_summary, file.path(out_tables, "10_section_balance_summary.tsv"))
message("\n--- Spatial balance by section ---"); print(balance_summary, n = Inf)

# =========================================================
# FIGURE 5D — spatial balance maps
# =========================================================
bal_limits <- sym_limits(spot_xy$balance)
message("\nShared colour limits: ", paste(round(bal_limits, 3), collapse = " to "))

make_map <- function(df_sec, title_txt, limits = bal_limits) {
  sec <- unique(df_sec$section_id)
  
  dd <- df_sec %>%
    filter(in_tissue == 1, is.finite(balance),
           is.finite(pxl_col_in_hires), is.finite(pxl_row_in_hires))
  if (nrow(dd) == 0) { warning("No plottable spots for section ", sec); return(NULL) }
  
  rec_w <- if (have_dims) unique(na.omit(dd$img_width_hires))  else numeric(0)
  rec_h <- if (have_dims) unique(na.omit(dd$img_height_hires)) else numeric(0)
  
  img <- NULL
  img_path <- unique(na.omit(df_sec$img_file))
  if (length(img_path) == 1 && file.exists(img_path)) {
    img <- tryCatch(read_image_any(img_path), error = function(e) {
      warning("Could not read image for section ", sec, ": ", conditionMessage(e)); NULL })
  }
  
  if (!is.null(img) && length(rec_w) == 1 && length(rec_h) == 1 &&
      (abs(ncol(img) - rec_w) > dim_tol || abs(nrow(img) - rec_h) > dim_tol)) {
    warning("Section ", sec, ": image is ", ncol(img), "x", nrow(img),
            " but coordinates assume ", rec_w, "x", rec_h,
            " - dropping the H&E layer.")
    img <- NULL
  }
  
  if (length(rec_w) == 1 && length(rec_h) == 1) {
    canvas_w <- rec_w; canvas_h <- rec_h
  } else if (!is.null(img)) {
    canvas_w <- ncol(img); canvas_h <- nrow(img)
  } else {
    canvas_w <- max(dd$pxl_col_in_hires, na.rm = TRUE)
    canvas_h <- max(dd$pxl_row_in_hires, na.rm = TRUE)
  }
  
  dd <- dd %>% mutate(x_plot = pxl_col_in_hires, y_plot = canvas_h - pxl_row_in_hires)
  
  p <- ggplot()
  if (!is.null(img)) {
    p <- p + annotation_raster(matrix(rgb(img[, , 1], img[, , 2], img[, , 3],
                                          alpha = he_alpha), nrow = nrow(img)),
                               xmin = 0, xmax = ncol(img), ymin = 0, ymax = nrow(img))
  }
  
  p +
    geom_point(data = dd, aes(x = x_plot, y = y_plot, colour = balance),
               size = point_size) +
    coord_fixed(xlim = c(0, canvas_w), ylim = c(0, canvas_h), expand = FALSE) +
    scale_colour_gradient2(low = "#F07167", mid = "grey92", high = "#00BF7D",
                           midpoint = 0, limits = limits, oob = scales::squish,
                           name = paste0(BAL_LAB, "\nbalance")) +
    labs(title = title_txt) +
    theme_void(base_size = 11) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          legend.position = "right")
}

main_secs <- intersect(main_figure_sections, unique(spot_xy$section_id))
if (length(main_secs) == 0) {
  warning("None of main_figure_sections found; using the first two available")
  main_secs <- head(sort(unique(spot_xy$section_id)), 2)
}

main_maps <- Filter(Negate(is.null), lapply(main_secs, function(s) {
  df <- spot_xy %>% filter(section_id == s)
  lab <- unique(df$section_label); if (length(lab) != 1) lab <- s
  make_map(df, lab)
}))

if (length(main_maps) > 0) {
  p_main <- wrap_plots(main_maps, nrow = 1, guides = "collect") +
    plot_annotation(
      title = "Spatial distribution of recognition-modality balance",
      subtitle = paste0("Green: ", AXIS_LAB_GD, " favoured; red: ", AXIS_LAB_AB,
                        " favoured. The immunosuppressive term cancels in this contrast."),
      theme = theme(plot.title = element_text(face = "bold", size = 13)))
  
  ggsave(file.path(out_plots, "Fig5D_spatial_balance.png"), p_main,
         width = 5.2 * length(main_maps) + 1.4, height = 5.4, dpi = 600, bg = "white")
  ggsave(file.path(out_plots, "Fig5D_spatial_balance.pdf"), p_main,
         width = 5.2 * length(main_maps) + 1.4, height = 5.4,
         bg = "white", device = cairo_pdf)
}

all_maps <- Filter(Negate(is.null), lapply(balance_summary$section_id, function(s) {
  df <- spot_xy %>% filter(section_id == s)
  lab <- unique(df$section_label); if (length(lab) != 1) lab <- s
  make_map(df, lab)
}))

if (length(all_maps) > 0) {
  ncol_s <- 4
  p_all <- wrap_plots(all_maps, ncol = ncol_s, guides = "collect") +
    plot_annotation(title = "Recognition-modality balance across all sections",
                    theme = theme(plot.title = element_text(face = "bold", size = 13)))
  ggsave(file.path(out_plots, "S_spatial_balance_all_sections.png"), p_all,
         width = 4.2 * ncol_s + 1.4,
         height = 4.2 * ceiling(length(all_maps) / ncol_s),
         dpi = 600, bg = "white", limitsize = FALSE)
  ggsave(file.path(out_plots, "S_spatial_balance_all_sections.pdf"), p_all,
         width = 4.2 * ncol_s + 1.4,
         height = 4.2 * ceiling(length(all_maps) / ncol_s),
         bg = "white", device = cairo_pdf, limitsize = FALSE)
}

writeLines(c(
  paste0("input: ", spot_rds),
  paste0("spatial_root: ", SPATIAL_ROOT),
  paste0("images_found: ", sum(img_lookup$found), "/", nrow(img_lookup)),
  paste0("zones_scored: ", paste(zones_keep, collapse = ",")),
  "Fig5C uses component scores, not nets (inhibitory term carries a mesenchymal gradient)",
  "zone_test: paired Wilcoxon signed-rank on section-level summaries",
  paste0("patient_sensitivity_relapses: ", paste(relapse_sections, collapse = ",")),
  "patient_sensitivity: one primary section per patient (n=11) plus leave-one-relapse-out",
  "component summary: median; mechanism axis summary: mean",
  paste0("depth_residualised: ", RESIDUALISE_DEPTH),
  paste0("gdT_axes_tested: ", paste(GD_AXES, collapse = ",")),
  paste0("inhibitory_axes_tested: ", paste(INH_AXES, collapse = ",")),
  "zone/opportunity gene overlap: ICAM1, IRF1, PTGS2, VEGFA",
  paste0("main_figure_sections: ", paste(main_secs, collapse = ",")),
  paste0("colour_limits_balance: ", paste(round(bal_limits, 3), collapse = " to ")),
  paste0("run_date: ", as.character(Sys.Date()))
), file.path(out_tables, "10_run_parameters.txt"))

message("\nDone. Outputs written to: ", out_root)