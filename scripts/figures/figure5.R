# =========================================================
# figure5.R — main figure assembly 
#
# PANELS
#   a  Opportunity quadrant across sections            (phase 09)
#   b  gdT recognition modalities by section           (phase 09)
#   c  Opportunity components by tumor compartment     (phase 10)
#   d  gdT recognition modalities by compartment       (phase 10)  <- PROMOTED
#   e  Spatial balance, 459 primary vs relapse         (phase 10)
#
# WHY d IS IN THE MAIN FIGURE
# ---------------------------
# Panel c shows gdT ligand availability as flat across compartments, which
# reads as "no spatial structure to gdT opportunity". Panel d shows why that is
# wrong: the modalities are strongly compartmentalised but in OPPOSITE
# directions - mevalonate and NKG2D ligand epithelial-high, EPHA2 and DNAM-1
# mesenchymal-high, butyrophilin flat - so the composite cancels. Without d a
# null is reported where there is a mechanism.
#
# SUMMARY STATISTICS DIFFER BY PANEL, DELIBERATELY
# ------------------------------------------------
#   a  section MEDIAN of the net scores. This matters: phase 09 calls
#      classify_opportunity() on abT_net_median / gdT_net_median, so plotting
#      means against median-derived categories could put a point on the wrong
#      side of a boundary. The columns are named explicitly below for that
#      reason - 09_section_summary.tsv contains BOTH _mean and _median.
#   b  section MEAN. Several axes retain one sparsely detected gene, for which
#      the section median collapses to zero on the z scale.
#   c/d  whatever phase 10 wrote into 10_zone_summaries.tsv 
# NOTHING IS SUBTRACTED IN c OR d
# -------------------------------
# Components are shown separately rather than as net scores because the
# immunosuppressive term carries a mesenchymal gradient; subtracting it would
# impose a compartment difference on the gdT score. Panel e is exempt - in the
# balance contrast the inhibitory term cancels exactly.
#
# Reads pre-computed outputs only. Run order: 09 -> 10 -> this script.
# =========================================================

source(here::here("config.R"))

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(tibble)
  library(ggplot2); library(ggrepel); library(patchwork)
  library(scales); library(grid); library(png)
})

# =========================================================
# 0. Paths
# =========================================================
phase09_root   <- stage_dir("scoring")
phase09_tables <- file.path(phase09_root, "tables")
phase09_rds    <- file.path(phase09_root, "rds")

phase10_root   <- stage_dir("zone_spatial")
phase10_tables <- file.path(phase10_root, "tables")
phase10_rds    <- file.path(phase10_root, "rds")

out_dir <- ensure_dir(figure_path("figure5"))

# =========================================================
# 1. Settings
# =========================================================
all_sections <- c(
  "1101", "1239", "1269", "1513", "459", "459_2", "723", "723_2",
  "727", "812", "821", "848", "928", "928_2"
)
stopifnot(setequal(all_sections, ALL_SECTIONS))

BALANCE_SECTIONS <- c("459", "459_2")   # panel e

relapse_sections <- vapply(MATCHED_PAIRS, `[[`, character(1), "relapse")
primary_of_pair  <- vapply(MATCHED_PAIRS, `[[`, character(1), "primary")

zone_levels <- c("Epithelial", "Mesenchymal")
zone_cols   <- c(Epithelial = "#4DAF4A", Mesenchymal = "#E41A1C")

# Keys must match the strings produced by display_category() in phase 09.
# The script checks this and tells you what it actually found.
category_cols <- c(
  "\u03b1\u03b2T-favorable" = "#F8766D",
  "Dual-opportunity"       = "#A3A500",
  "\u03b3\u03b4T-favorable" = "#00BF7D",
  "Immuno-cold"            = "#00B0F6",
  "Suppressed"             = "#E76BF3"
)

AXIS_LAB_AB <- "\u03b1\u03b2T net"
AXIS_LAB_GD <- "\u03b3\u03b4T net"

RESID_SUFFIX <- "__resid"     # double underscore, as written by phase 10

COMPONENTS <- c(
  abT_core_z           = "\u03b1\u03b2T ligand availability",
  gdT_core_z           = "\u03b3\u03b4T ligand availability",
  inhibitory_program_z = "Immunosuppressive tone"
)

# IPP drain (FDPS) is a MODIFIER and not part of the composite.
GD_AXES <- c(
  gdT_btn           = "Butyrophilin (BTN2A1/3A1)",
  gdT_mevalonate    = "Mevalonate / IPP synthesis",
  gdT_nkg2d         = "NKG2D ligand (MICA)",
  gdT_ephrin        = "EPHA2",
  gdT_adhesion_dnam = "DNAM-1 / adhesion",
  gdT_ipp_drain     = "FDPS (IPP consumption)"
)

HE_WASH_ALPHA <- 0.42
PT_SIZE       <- 1.35
HE_PAD_FRAC   <- 0.015

# =========================================================
# 2. Helpers
# =========================================================
fmt_p <- function(p, digits = 2)
  ifelse(is.na(p), "p = NA",
         paste0("p = ", format.pval(p, digits = digits, eps = 1e-300)))

# Errors on ambiguity rather than silently taking the first hit. Several phase-09
# tables carry both _mean and _median variants of the same score.
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
         width = width, height = height, dpi = 600, bg = "white")
  ggsave(file.path(out_dir, paste0(name, ".pdf")), p,
         width = width, height = height, device = pdf_device, bg = "white")
}

# =========================================================
# 3. Load
# =========================================================
message("Loading phase 09 / 10 outputs...")

section_summary <- read_tsv(file.path(phase09_tables, "09_section_summary.tsv"),
                            show_col_types = FALSE)
spot <- readRDS(file.path(phase09_rds, "09_spot_scores.rds")) %>% as_tibble()

zone_summ <- read_tsv(file.path(phase10_tables, "10_zone_summaries.tsv"),
                      show_col_types = FALSE)

sec_col_spot <- resolve_col(spot, "^section_id$", "section id")
spot <- spot %>% mutate(section = as.character(.data[[sec_col_spot]]))

# =========================================================
# PANEL A — opportunity quadrant
# MEDIANS, named explicitly: the categories were assigned from the medians.
# =========================================================
message("Panel a")

require_cols(section_summary,
             c("section_id", "abT_net_median", "gdT_net_median",
               "opportunity_label"),
             "09_section_summary.tsv")

quad <- section_summary %>%
  transmute(section  = as.character(section_id),
            ab       = abT_net_median,
            gd       = gdT_net_median,
            category = as.character(opportunity_label))

unknown_cat <- setdiff(unique(quad$category), names(category_cols))
if (length(unknown_cat) > 0)
  stop("Category labels not in category_cols: ",
       paste(unknown_cat, collapse = ", "),
       "\nUpdate category_cols to match display_category() in phase 09.")

p_a <- ggplot(quad, aes(ab, gd, colour = category)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(aes(label = section), size = 3, show.legend = FALSE,
                           min.segment.length = 0.2, box.padding = 0.45,
                           point.padding = 0.35, max.overlaps = Inf, seed = 1) +
  scale_colour_manual(values = category_cols, name = "Opportunity\ncategory") +
  labs(x = AXIS_LAB_AB, y = AXIS_LAB_GD,
       title = "Immunotherapy opportunity across ependymoma sections",
       subtitle = paste0("Section medians. Scores are z-scaled within the ",
                         "cohort, so categories are relative to these ",
                         nrow(quad), " sections.")) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey35"),
        panel.grid.minor = element_blank())

save_panel(p_a, "Fig5a_opportunity_quadrant", 7.5, 5.6)

# =========================================================
# PANEL B — gdT modalities by section (section MEAN)
# =========================================================
message("Panel b")

gd_cols <- names(GD_AXES)[names(GD_AXES) %in% names(spot)]
if (length(gd_cols) == 0)
  stop("No gdT axis columns in the spot table. Available: ",
       paste(names(spot), collapse = ", "))
if (length(gd_cols) < length(GD_AXES))
  warning("gdT axes missing from the spot table and skipped: ",
          paste(setdiff(names(GD_AXES), gd_cols), collapse = ", "), call. = FALSE)

axis_by_section <- spot %>%
  select(section, all_of(gd_cols)) %>%
  pivot_longer(-section, names_to = "axis", values_to = "score") %>%
  group_by(section, axis) %>%
  summarise(value = mean(score, na.rm = TRUE), .groups = "drop") %>%
  mutate(axis  = factor(unname(GD_AXES[axis]), levels = unname(GD_AXES[gd_cols])),
         label = label_section(section),
         pos   = value >= 0)

lab_order <- axis_by_section %>% filter(axis == levels(axis)[1]) %>%
  arrange(value) %>% pull(label)
axis_by_section$label <- factor(axis_by_section$label, levels = lab_order)

p_b <- ggplot(axis_by_section, aes(value, label, fill = pos)) +
  geom_col(width = 0.72) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey55", linewidth = 0.35) +
  facet_wrap(~ axis, nrow = 2, scales = "free_x") +
  scale_fill_manual(values = c(`TRUE` = "#00BF7D", `FALSE` = "#A6CEE3"),
                    guide = "none") +
  labs(x = "Section mean score", y = NULL,
       title = paste0(AXIS_LAB_GD, " recognition modalities by section"),
       subtitle = paste0("Modalities are independent by design; the composite ",
                         "weights each equally. FDPS is a modifier and is not ",
                         "part of the composite.")) +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey35"),
        strip.background = element_rect(fill = "grey93", colour = NA),
        panel.grid.minor = element_blank())

save_panel(p_b, "Fig5b_gdT_mechanism_axes", 11, 6)

# =========================================================
# PANELS C and D — by tumor compartment
# 10_zone_summaries.tsv is WIDE: one row per section x zone, scores as columns.
# =========================================================
message("Panels c/d")

require_cols(zone_summ, c("section_id", "zone_call"), "10_zone_summaries.tsv")

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

# Tests computed from the plotted values, so annotation and data cannot diverge.
# Compare the printed table against 10_zone_tests.tsv on first run; a uniform
# shift means a different variant of the test, not different data.
zone_tests <- read_tsv(file.path(phase10_tables, "10_zone_tests.tsv"),
                       show_col_types = FALSE)

zone_ann <- function(keep_map, version_keep = "raw", p_col = "p_adj") {
  zone_tests %>%
    filter(version == version_keep, score_base %in% names(keep_map)) %>%
    transmute(score = factor(unname(keep_map[score_base]),
                             levels = unname(keep_map)),
              lab = paste0(fmt_p(.data[[p_col]]),
                           "\nepithelial higher: ", n_higher_epithelial,
                           "/", n_sections))
}

compartment_panel <- function(keep_map, title_txt, subtitle_txt, stat_lab,
                              use_resid = FALSE) {
  d   <- zone_long(keep_map, use_resid)
  ann <- zone_ann(keep_map,
                  version_keep = if (use_resid) "depth-residualised" else "raw")
  
  
  # Place the annotation just inside the top of each facet rather than at Inf,
  # which gets clipped when the panel is tall.
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

p_c <- compartment_panel(
  COMPONENTS,
  "Opportunity components by tumor compartment",
  paste0("Components are shown separately rather than as net scores: the ",
         "immunosuppressive term carries a mesenchymal gradient, so subtracting ",
         "it would impose a compartment difference on the ", AXIS_LAB_GD,
         " score.\nZone and opportunity programs were assembled independently; ",
         "three scores share a gene with the mesenchymal zone program (DNAM-1/",
         "adhesion: ICAM1; VEGF: VEGFA; COX-2: PTGS2). Each line is one section."),
  "Section median score")

save_panel(p_c, "Fig5c_compartment_components", 10, 5.6)

p_d <- compartment_panel(
  GD_AXES,
  paste0(AXIS_LAB_GD, " recognition modalities by tumor compartment"),
  paste0("The composite in panel c is flat across compartments because its ",
         "modalities are compartmentalised in opposite directions.\n",
         "y scales differ between panels. No score is subtracted from another. ",
         "FDPS is a modifier and is not part of the composite."),
  "Section mean score")

save_panel(p_d, "Fig5d_gdT_axes_by_compartment", 10, 6)

# =========================================================
# PANEL E — spatial balance, 459 primary vs relapse
# balance = gdT_net - abT_net; the inhibitory term cancels exactly.
#
# Colour limits follow phase 10 exactly: symmetric, from the 2nd and 98th
# percentiles of the balance across ALL sections, so the main-figure map and
# the all-sections supplementary map share one scale.
# =========================================================
message("Panel e")

coord_rds <- file.path(phase10_rds,    "10_spot_scores_with_coordinates.rds")
coord_tsv <- file.path(phase10_tables, "10_spot_scores_with_coordinates.tsv")
coords <- if (file.exists(coord_rds)) as_tibble(readRDS(coord_rds)) else
  read_tsv(coord_tsv, show_col_types = FALSE)

c_sec <- resolve_col(coords, "^section_id$",       "section id")
c_x   <- resolve_col(coords, "^pxl_col_in_hires$", "x coordinate")
c_y   <- resolve_col(coords, "^pxl_row_in_hires$", "y coordinate")
c_bal <- resolve_col(coords, "^balance$",          "balance", required = FALSE)
if (is.na(c_bal)) {
  gd_n <- resolve_col(coords, "^gdT_net$", "gdT net (spot level)")
  ab_n <- resolve_col(coords, "^abT_net$", "abT net (spot level)")
  coords$balance <- coords[[gd_n]] - coords[[ab_n]]
  c_bal <- "balance"
}

# Verbatim from 10_zone_and_spatial.R (q_clip = 0.02)
Q_CLIP <- 0.02
sym_limits <- function(x, q = Q_CLIP) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(c(-1, 1))
  lim <- max(abs(quantile(x, c(q, 1 - q), na.rm = TRUE)))
  if (!is.finite(lim) || lim == 0) lim <- max(abs(x), na.rm = TRUE)
  if (!is.finite(lim) || lim == 0) lim <- 1
  c(-lim, lim)
}

bal_limits <- sym_limits(coords[[c_bal]])
message("  shared colour limits: ",
        paste(round(bal_limits, 3), collapse = " to "))

make_balance_map <- function(sec, limits = bal_limits) {
  d <- coords %>% filter(as.character(.data[[c_sec]]) == sec)
  if (nrow(d) == 0) stop("No spots for section ", sec)
  
  ipath <- file.path(RAW_ROOT, sec, "spatial", "tissue_hires_image.png")
  use_he <- file.exists(ipath)
  
  if (use_he) {
    img <- png::readPNG(ipath); ih <- nrow(img); iw <- ncol(img)
    d   <- d %>% mutate(xp = .data[[c_x]], yp = ih - .data[[c_y]])
    pad <- ceiling(max(ih, iw) * HE_PAD_FRAC)
    p <- ggplot() +
      annotation_custom(grid::rasterGrob(img, interpolate = TRUE),
                        xmin = 0, xmax = iw, ymin = 0, ymax = ih) +
      annotate("rect", xmin = 0, xmax = iw, ymin = 0, ymax = ih,
               fill = "white", alpha = HE_WASH_ALPHA) +
      coord_fixed(xlim = c(-pad, iw + pad), ylim = c(-pad, ih + pad), expand = FALSE)
  } else {
    message("  no H&E for ", sec, " - coordinate-only fallback")
    d <- d %>% mutate(xp = .data[[c_x]], yp = -.data[[c_y]])
    p <- ggplot() + coord_fixed()
  }
  
  p +
    geom_point(data = d, aes(xp, yp, colour = .data[[c_bal]]),
               size = PT_SIZE, alpha = 0.95) +
    scale_colour_gradient2(low = "#E41A1C", mid = "grey96", high = "#00BF7D",
                           midpoint = 0, limits = limits,
                           oob = scales::squish,
                           name = paste0(AXIS_LAB_GD, " \u2212 ", AXIS_LAB_AB,
                                         "\nbalance")) +
    ggtitle(label_section(sec)) +
    theme_void(base_size = 10) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 10))
}

p_e <- wrap_plots(lapply(BALANCE_SECTIONS, make_balance_map),
                  nrow = 1, guides = "collect") +
  plot_annotation(
    title = "Spatial distribution of recognition-modality balance",
    subtitle = paste0("Green: ", AXIS_LAB_GD, " favoured; red: ", AXIS_LAB_AB,
                      " favoured. The immunosuppressive term cancels exactly ",
                      "in this contrast. Colour limits are shared across all ",
                      "sections and clipped at the 2nd and 98th percentiles."),
    theme = theme(plot.title = element_text(face = "bold", size = 11),
                  plot.subtitle = element_text(size = 8, colour = "grey35"))) &
  theme(legend.position = "right")

save_panel(p_e, "Fig5e_spatial_balance_459", 11, 5)

# =========================================================
# 4. Reporting aids
# =========================================================
cat("\n--- Panel a values (section medians, as classified) ---\n")
quad %>% arrange(desc(gd)) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>% print()

cat("\n--- Compartment tests as plotted (from 10_zone_tests.tsv, BH within family) ---\n")
zone_tests %>%
  filter(version == "raw",
         score_base %in% c(names(COMPONENTS), names(GD_AXES))) %>%
  select(level, label, n_sections, n_higher_epithelial, p_value, p_adj) %>%
  as.data.frame() %>% print()

cat("\nNOTE: panel a uses section MEDIANS (the statistic the categories were",
    "\nassigned from); panel b uses section MEANS. Phase 10 uses the median for",
    "\ncomponents and the mean for mechanism axes. Say so in the legend.\n")

cat("\nNOTE: Figure 5 uses all", n_distinct(quad$section),
    "sections including 848, which Figure 2 excludes.\nThe 848 exclusion",
    "concerned Lymphocyte score inflation; these are tumor-side scores.\n")

message("\nDone. Panels written to: ", normalizePath(out_dir))