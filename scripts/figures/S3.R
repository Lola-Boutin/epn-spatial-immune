# =========================================================
# S3.R — supplementary figure for the opportunity analysis
#
# PANELS, ordered by first citation in the Results
#   a  Immunosuppressive tone by section                (phase 09)
#   b  Immunosuppressive mechanisms by section          (phase 09)
#   c  Matched primary vs relapse, BY COMPONENT         (phase 09)
#   d  Opportunity components by compartment, depth-residualised  (phase 10)
#   e  Immunosuppressive mechanisms by compartment      (phase 10)
#
# The all-sections spatial balance map is no longer produced here. It stays in
# the analysis output (10_zone_and_spatial) - fourteen thumbnails are not
# legible at supplementary size, and the claim is already carried by Fig 5e.
#
# CONVENTIONS, MATCHED TO THE ANALYSIS SCRIPTS
# --------------------------------------------
#   * Compartment annotation comes from 10_zone_tests.tsv, not recomputed here.
#     p_adj is Benjamini-Hochberg WITHIN each family (phase 10 calls p.adjust
#     once per level x version); n_higher_epithelial is the direction count.
#   * Component panels use section MEDIANS; mechanism-axis panels use section
#     MEANS, because for sparsely detected axes the section median reduces to
#     zero on the z scale. The y labels say which.
#
# PANEL c - TWO DELIBERATE CHOICES
# --------------------------------
# 1. Components, not net scores: net = ligand availability - immunosuppressive
#    tone, so a rise at relapse can come from more ligand OR less suppression.
# 2. No p value. n = 3 patients; the smallest attainable p from a paired sign
#    or signed-rank test is 0.25. Direction and n are reported instead, with a
#    depth-coupling check printed to the console.
#
# Reads pre-computed outputs only. Run order: 09 -> 10 -> this script.
# =========================================================

source(here::here("config.R"))

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(tibble)
  library(ggplot2); library(patchwork); library(scales)
})

# =========================================================
# 0. Paths
# =========================================================
phase09_root   <- stage_dir("scoring")
phase09_tables <- file.path(phase09_root, "tables")
phase09_rds    <- file.path(phase09_root, "rds")

phase10_root   <- stage_dir("zone_spatial")
phase10_tables <- file.path(phase10_root, "tables")

out_dir <- ensure_dir(figure_path("S3"))

# =========================================================
# 1. Settings
# =========================================================
pairs_map <- tibble(
  patient = paste0("P", vapply(MATCHED_PAIRS, `[[`, character(1), "primary")),
  primary = vapply(MATCHED_PAIRS, `[[`, character(1), "primary"),
  relapse = vapply(MATCHED_PAIRS, `[[`, character(1), "relapse")
)
relapse_sections <- pairs_map$relapse
primary_of_pair  <- pairs_map$primary

zone_levels <- c("Epithelial", "Mesenchymal")
zone_cols   <- c(Epithelial = "#4DAF4A", Mesenchymal = "#E41A1C")

ELEVATED_TONE_Z <- 0.75
RESID_SUFFIX    <- "__resid"
USE_DEPTH_RESIDUALISED_FOR_PAIRS <- TRUE

COMPONENTS <- c(
  abT_core_z           = "\u03b1\u03b2T ligand availability",
  gdT_core_z           = "\u03b3\u03b4T ligand availability",
  inhibitory_program_z = "Immunosuppressive tone"
)

# Labels taken from the `label` column of 10_zone_tests.tsv so figure and table
# use the same names.
INH_AXES <- c(
  inh_mhc_checkpoint = "HLA-E / NKG2A",
  inh_classical_ckpt = "Galectin-9 (LGALS9)",
  inh_tgfb           = "TGF-\u03b2 (TGFB1/2/3)",
  inh_adenosine      = "Adenosine (NT5E/ENTPD1)",
  inh_pge2           = "COX-2 (PTGS2)",
  inh_angiogenic     = "VEGF (VEGFA)"
)

# =========================================================
# 2. Helpers
# =========================================================
fmt_p <- function(p, digits = 2)
  ifelse(is.na(p), "p = NA",
         paste0("p = ", format.pval(p, digits = digits, eps = 1e-300)))

# Errors on ambiguity rather than silently taking the first hit: several
# phase-09 tables carry both _mean and _median variants of the same score.
resolve_col <- function(df, wanted, what = "column", required = TRUE) {
  if (wanted %in% names(df)) return(wanted)
  hit <- grep(wanted, names(df), ignore.case = TRUE, value = TRUE)
  if (length(hit) == 1) {
    message("  resolved ", what, ": '", wanted, "' -> '", hit, "'")
    return(hit)
  }
  if (length(hit) > 1)
    stop("Ambiguous ", what, " '", wanted, "': ", paste(hit, collapse = ", "),
         "\nName it explicitly in the script.")
  if (required)
    stop("Could not find ", what, " '", wanted, "'. Available: ",
         paste(names(df), collapse = ", "))
  NA_character_
}

require_cols <- function(df, need, where) {
  miss <- setdiff(need, names(df))
  if (length(miss) > 0)
    stop("Missing from ", where, ": ", paste(miss, collapse = ", "),
         "\nAvailable: ", paste(names(df), collapse = ", "))
  invisible(TRUE)
}

label_section <- function(s)
  ifelse(s %in% relapse_sections, paste0(s, " (Relapse)"),
         ifelse(s %in% primary_of_pair, paste0(s, " (Primary)"), s))

pdf_device <- if (isTRUE(capabilities("cairo"))) cairo_pdf else "pdf"

save_panel <- function(p, name, width, height) {
  ggsave(file.path(out_dir, paste0(name, ".png")), p,
         width = width, height = height, dpi = PLOT$dpi, bg = "white")
  ggsave(file.path(out_dir, paste0(name, ".pdf")), p,
         width = width, height = height, device = pdf_device, bg = "white")
}

# =========================================================
# 3. Load
# =========================================================
message("Loading phase 09 / 10 outputs...")
spot <- readRDS(file.path(phase09_rds, "09_spot_scores.rds")) %>% as_tibble()
sec_col <- resolve_col(spot, "^section_id$", "section id")
spot <- spot %>% mutate(section = as.character(.data[[sec_col]]))

zone_summ  <- read_tsv(file.path(phase10_tables, "10_zone_summaries.tsv"),
                       show_col_types = FALSE)
zone_tests <- read_tsv(file.path(phase10_tables, "10_zone_tests.tsv"),
                       show_col_types = FALSE)

require_cols(zone_summ, c("section_id", "zone_call"), "10_zone_summaries.tsv")
require_cols(zone_tests,
             c("score_base", "version", "n_sections", "n_higher_epithelial",
               "p_value", "p_adj"),
             "10_zone_tests.tsv")

# =========================================================
# PANEL A — immunosuppressive tone by section (section MEDIAN)
# First citation: the sentence identifying 1239 as the only Suppressed section.
# =========================================================
message("Panel a")
inh_prog <- resolve_col(spot, "^inhibitory_program_z$", "inhibitory composite",
                        required = FALSE)
if (is.na(inh_prog))
  inh_prog <- resolve_col(spot, "^inhibitory_program$", "inhibitory composite")

tone <- spot %>%
  group_by(section) %>%
  summarise(tone = median(.data[[inh_prog]], na.rm = TRUE), .groups = "drop") %>%
  mutate(label    = label_section(section),
         elevated = ifelse(tone >= ELEVATED_TONE_Z, "Yes", "No")) %>%
  arrange(tone) %>%
  mutate(label = factor(label, levels = label))

pS_a <- ggplot(tone, aes(tone, label, fill = elevated)) +
  geom_col(width = 0.72) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey55", linewidth = 0.35) +
  geom_vline(xintercept = ELEVATED_TONE_Z, linetype = 3, colour = "grey25") +
  scale_fill_manual(values = c(No = "#A6CEE3", Yes = "#E76BF3"), name = "Elevated") +
  labs(x = "Median immunosuppressive tone (z)", y = NULL,
       title = "Immunosuppressive tone across sections",
       subtitle = paste0("Dotted line: elevated-tone threshold (z = ",
                         ELEVATED_TONE_Z, "), the criterion separating the ",
                         "Suppressed category.\nScores are z-scaled within the ",
                         "cohort, so this is relative to these ", nrow(tone),
                         " sections.")) +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey35"),
        panel.grid.minor = element_blank())

save_panel(pS_a, "S3a_inhibitory_tone_by_section", 7.5, 5.5)

# =========================================================
# PANEL B — immunosuppressive mechanisms by section (section MEAN)
# Shows that 1239's tone reflects concurrent elevation across mechanisms
# rather than one dominant axis.
# =========================================================
message("Panel b")
inh_cols <- names(INH_AXES)[names(INH_AXES) %in% names(spot)]
if (length(inh_cols) == 0)
  stop("No inhibitory axis columns in the spot table. Available: ",
       paste(names(spot), collapse = ", "))
if (length(inh_cols) < length(INH_AXES))
  warning("Inhibitory axes missing and skipped: ",
          paste(setdiff(names(INH_AXES), inh_cols), collapse = ", "), call. = FALSE)

inh_by_sec <- spot %>%
  select(section, all_of(inh_cols)) %>%
  pivot_longer(-section, names_to = "axis", values_to = "score") %>%
  group_by(section, axis) %>%
  summarise(value = mean(score, na.rm = TRUE), .groups = "drop") %>%
  mutate(axis  = factor(unname(INH_AXES[axis]), levels = unname(INH_AXES[inh_cols])),
         label = label_section(section),
         pos   = value >= 0)

ord <- inh_by_sec %>% filter(axis == levels(axis)[1]) %>%
  arrange(value) %>% pull(label)
inh_by_sec$label <- factor(inh_by_sec$label, levels = ord)

pS_b <- ggplot(inh_by_sec, aes(value, label, fill = pos)) +
  geom_col(width = 0.72) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey55", linewidth = 0.35) +
  facet_wrap(~ axis, nrow = 2, scales = "free_x") +
  scale_fill_manual(values = c(`TRUE` = "#00BF7D", `FALSE` = "#A6CEE3"),
                    guide = "none") +
  labs(x = "Section mean score", y = NULL,
       title = "Immunosuppressive mechanisms by section",
       subtitle = paste0("Mechanisms are independent by design; the composite ",
                         "weights each equally. Axes are summarised as means ",
                         "because for sparsely detected axes\nthe section ",
                         "median reduces to zero on the z scale.")) +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey35"),
        strip.background = element_rect(fill = "grey93", colour = NA),
        panel.grid.minor = element_blank())

save_panel(pS_b, "S3b_inhibitory_mechanism_axes", 11, 6)

# Ranking of 1239 per mechanism, for the Results sentence
cat("\n--- Rank of each section per inhibitory mechanism (1 = highest) ---\n")
inh_by_sec %>%
  group_by(axis) %>%
  mutate(rank = rank(-value, ties.method = "min")) %>%
  filter(section == "1239") %>%
  select(axis, value, rank) %>%
  mutate(value = round(value, 3)) %>%
  as.data.frame() %>% print()

# =========================================================
# PANEL C — matched primary vs relapse, by component
# =========================================================
message("Panel c")

comp_cols   <- names(COMPONENTS)[names(COMPONENTS) %in% names(spot)]
comp_labels <- COMPONENTS
if (length(comp_cols) == 0) {
  alt <- sub("_z$", "", names(COMPONENTS))
  hit <- alt[alt %in% names(spot)]
  if (length(hit) > 0) {
    comp_labels <- setNames(unname(COMPONENTS), alt)
    comp_cols   <- hit
    message("  using un-suffixed component names from the spot table")
  }
}
if (length(comp_cols) == 0)
  stop("No component columns in the spot table. Available: ",
       paste(names(spot), collapse = ", "))

depth_col <- intersect(c("nCount_Spatial", "nCount_RNA", "lib_size", "n_umi"),
                       names(spot))
if (length(depth_col) > 0) {
  dc <- depth_col[1]
  depth_by_sec <- spot %>%
    group_by(section) %>%
    summarise(med_log_depth = median(log(.data[[dc]] + 1), na.rm = TRUE),
              across(all_of(comp_cols), ~ median(.x, na.rm = TRUE)),
              .groups = "drop")
  cat("\n--- Depth coupling across sections (", dc, ") ---\n", sep = "")
  for (cc in comp_cols) {
    rho <- suppressWarnings(cor(depth_by_sec$med_log_depth, depth_by_sec[[cc]],
                                method = "spearman", use = "complete.obs"))
    cat(sprintf("  %-22s Spearman vs median log depth: %+.3f\n", cc, rho))
  }
  cat("  A strong positive value for abT ligand availability would mean the\n",
      "  primary-to-relapse increase could reflect depth rather than biology.\n", sep = "")
} else {
  message("  no depth column in the spot table - depth check skipped")
}

resid_cols <- paste0(comp_cols, RESID_SUFFIX)
use_resid  <- USE_DEPTH_RESIDUALISED_FOR_PAIRS && all(resid_cols %in% names(spot))
score_cols <- if (use_resid) resid_cols else comp_cols
names(score_cols) <- comp_cols
message(if (use_resid) "  using depth-residualised component scores"
        else "  using raw component scores (no residualised columns in the spot table)")

pair_long <- pairs_map %>%
  pivot_longer(c(primary, relapse), names_to = "timepoint", values_to = "section") %>%
  inner_join(
    spot %>%
      select(section, all_of(unname(score_cols))) %>%
      pivot_longer(-section, names_to = "raw_score", values_to = "v") %>%
      group_by(section, raw_score) %>%
      summarise(value = median(v, na.rm = TRUE), .groups = "drop"),
    by = "section") %>%
  mutate(base      = names(score_cols)[match(raw_score, unname(score_cols))],
         component = factor(unname(comp_labels[base]), levels = unname(comp_labels)),
         timepoint = factor(timepoint, levels = c("primary", "relapse"),
                            labels = c("Primary", "Relapse")))

dir_summary <- pair_long %>%
  select(patient, component, timepoint, value) %>%
  pivot_wider(names_from = timepoint, values_from = value) %>%
  mutate(up = Relapse > Primary) %>%
  group_by(component) %>%
  summarise(n_up = sum(up), n = n(), .groups = "drop") %>%
  mutate(lab = paste0("increased in ", n_up, "/", n, " patients"))

pS_c <- ggplot(pair_long, aes(timepoint, value, group = patient, colour = patient)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2.6) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60", linewidth = 0.35) +
  geom_text(data = dir_summary, aes(x = 1.5, y = Inf, label = lab),
            inherit.aes = FALSE, vjust = 1.4, size = 2.7, colour = "grey25") +
  facet_wrap(~ component, nrow = 1, scales = "free_y") +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.2))) +
  labs(x = NULL, colour = "Patient",
       y = paste0("Section median score",
                  if (use_resid) " (depth-residualised)" else ""),
       title = "Matched primary versus relapse sections",
       subtitle = paste0("Within-patient comparison, n = 3. Components are shown ",
                         "separately rather than as net scores, because a change in ",
                         "a net score\ncannot distinguish altered ligand availability ",
                         "from altered immunosuppressive tone. No significance test ",
                         "is reported:\nwith three pairs the smallest attainable p is 0.25.")) +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 7.5, colour = "grey35"),
        strip.background = element_rect(fill = "grey93", colour = NA),
        panel.grid.minor = element_blank())

save_panel(pS_c, "S3c_primary_relapse_components", 10, 4.4)

# =========================================================
# PANELS D and E — by tumor compartment
# 10_zone_summaries.tsv is WIDE: one row per section x zone, scores as columns.
# Annotation is read from 10_zone_tests.tsv (BH within family).
# =========================================================
message("Panels d/e")

zone_long <- function(keep_map, use_resid = FALSE) {
  cols <- names(keep_map)
  if (use_resid) cols <- paste0(cols, RESID_SUFFIX)
  have <- cols[cols %in% names(zone_summ)]
  if (length(have) == 0)
    stop("None of these columns are in 10_zone_summaries.tsv: ",
         paste(cols, collapse = ", "), "\nAvailable: ",
         paste(names(zone_summ), collapse = ", "))
  if (length(have) < length(cols))
    warning("Missing from the zone table and skipped: ",
            paste(setdiff(cols, have), collapse = ", "), call. = FALSE)
  
  zone_summ %>%
    filter(zone_call %in% zone_levels) %>%
    select(section_id, zone_call, all_of(have)) %>%
    pivot_longer(all_of(have), names_to = "col", values_to = "value") %>%
    mutate(base    = sub(paste0(RESID_SUFFIX, "$"), "", col),
           score   = factor(unname(keep_map[base]), levels = unname(keep_map)),
           section = as.character(section_id),
           zone    = factor(zone_call, levels = zone_levels)) %>%
    select(section, zone, score, value)
}

zone_ann <- function(keep_map, version_keep = "raw", p_col = "p_adj") {
  a <- zone_tests %>%
    filter(version == version_keep, score_base %in% names(keep_map)) %>%
    transmute(score = factor(unname(keep_map[score_base]),
                             levels = unname(keep_map)),
              lab = paste0(fmt_p(.data[[p_col]]),
                           "\nepithelial higher: ", n_higher_epithelial,
                           "/", n_sections))
  if (nrow(a) == 0)
    stop("No rows in 10_zone_tests.tsv for version '", version_keep,
         "' and scores: ", paste(names(keep_map), collapse = ", "))
  a
}

compartment_panel <- function(keep_map, title_txt, subtitle_txt, stat_lab,
                              use_resid = FALSE) {
  d   <- zone_long(keep_map, use_resid)
  ann <- zone_ann(keep_map,
                  version_keep = if (use_resid) "depth-residualised" else "raw")
  
  # Placed just inside the top of each facet; at Inf it gets clipped.
  ypos <- d %>% group_by(score) %>%
    summarise(y = max(value, na.rm = TRUE) +
                0.16 * diff(range(value, na.rm = TRUE)), .groups = "drop")
  ann <- left_join(ann, ypos, by = "score")
  
  ggplot(d, aes(zone, value)) +
    geom_boxplot(aes(fill = zone), outlier.shape = NA, alpha = 0.85,
                 width = 0.6, linewidth = 0.35) +
    geom_line(aes(group = section), colour = "grey55", linewidth = 0.3) +
    geom_point(size = 1.2, colour = "grey15") +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey60", linewidth = 0.35) +
    geom_text(data = ann, aes(x = 1.5, y = y, label = lab),
              inherit.aes = FALSE, vjust = 1, size = 2.6, lineheight = 0.95) +
    facet_wrap(~ score, nrow = 2, scales = "free_y") +
    scale_fill_manual(values = zone_cols, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.24))) +
    labs(x = NULL, y = stat_lab, title = title_txt, subtitle = subtitle_txt) +
    theme_bw(base_size = 9) +
    theme(plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(size = 7.5, colour = "grey35"),
          strip.background = element_rect(fill = "grey93", colour = NA),
          panel.grid.minor = element_blank())
}

resid_present <- all(paste0(names(COMPONENTS), RESID_SUFFIX) %in% names(zone_summ))
if (resid_present) {
  pS_d <- compartment_panel(
    COMPONENTS,
    "Opportunity components by compartment, depth-residualised",
    paste0("Each score regressed on log sequencing depth within section before ",
           "summarising. Between-section depth differences are retained ",
           "deliberately, being confounded with patient.\nPaired Wilcoxon ",
           "signed-rank across the 14 sections; p values are ",
           "Benjamini-Hochberg adjusted within the three components."),
    "Section median score", use_resid = TRUE)
  save_panel(pS_d, "S3d_compartment_depth_residualised", 10, 5.6)
} else {
  message("No ", RESID_SUFFIX, " component columns in the zone table - skipping panel d")
}

pS_e <- compartment_panel(
  INH_AXES,
  "Immunosuppressive mechanisms by tumor compartment",
  paste0("Paired Wilcoxon signed-rank across the 14 sections; p values are ",
         "Benjamini-Hochberg adjusted within the six inhibitory mechanisms.\n",
         "y scales differ between panels. VEGF (VEGFA) and COX-2 (PTGS2) also ",
         "appear in the mesenchymal zone program, so those two comparisons are ",
         "not fully independent of the compartment definition."),
  "Section mean score")

save_panel(pS_e, "S3e_inhibitory_axes_by_compartment", 10, 6)

# =========================================================
# 4. Reporting aids
# =========================================================
cat("\n--- Primary vs relapse, direction only (n = 3) ---\n")
pair_long %>%
  select(patient, component, timepoint, value) %>%
  pivot_wider(names_from = timepoint, values_from = value) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)),
         direction = ifelse(Relapse > Primary, "up", "down")) %>%
  as.data.frame() %>% print()

cat("\n--- Compartment tests used in panels d and e ---\n")
zone_tests %>%
  filter((version == "depth-residualised" & score_base %in% names(COMPONENTS)) |
           (version == "raw" & score_base %in% names(INH_AXES))) %>%
  select(any_of(c("level", "version", "label", "n_sections",
                  "n_higher_epithelial", "p_value", "p_adj"))) %>%
  as.data.frame() %>% print()

message("\nDone. Panels written to: ", normalizePath(out_dir))
