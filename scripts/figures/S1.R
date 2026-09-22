# =========================================================
# Supplementary Figure 1 — unified script
#
# Panels:
#   a  Epithelial neighbourhood spatial maps (459 & 928)
#   b  Epithelial neighbourhood permutation bar chart
#   c  Vascular neighbourhood spatial maps (459 & 928)
#   d  Vascular neighbourhood permutation bar chart
#   e  Section 723 — all four zone neighbourhood maps (1 row)
#
# Reads permutation results from Figure3_all_permutation_results.tsv
# (written by figure3.R). If that file is not yet present, the
# script recomputes the four zone-neighbourhood permutation families independently.
# =========================================================

source(here::here("config.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(readr)
  library(tidyr)
  library(grid)
  library(png)
})

# =========================================================
# 0. Paths
# =========================================================
phase06_tables <- file.path(stage_dir("zones"), "tables")
fig3_out_dir   <- figure_path("figure3")
out_dir        <- ensure_dir(figure_path("S1"))

# =========================================================
# 1. Settings
# =========================================================
good_sections    <- GOOD_SECTIONS
example_sections <- c("459", "928")

zone_levels <- c("Epithelial", "Mesenchymal", "Vascular", "Myeloid", "Uncertain")
zone_cols   <- c(
  "Epithelial"  = "#4DAF4A",
  "Mesenchymal" = "#E41A1C",
  "Vascular"    = "#377EB8",
  "Myeloid"     = "#984EA3",
  "Uncertain"   = "grey70"
)

PT_SIZE_TISSUE  <- 0.55
HOT_RING_SIZE   <- 1.00
HOT_RING_STROKE <- 0.80
HE_PAD_FRAC     <- 0.015
N_PERM          <- 1000L
use_seed("neighborhood_permutation")

# =========================================================
# 2. Helpers (identical to figure3.R)
# =========================================================
stars <- function(p) {
  dplyr::case_when(
    is.na(p)  ~ "",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ "ns"
  )
}

fmt_p <- function(p, digits = 2) {
  ifelse(is.na(p), "p=NA",
         paste0("p=", format.pval(p, digits = digits, eps = 1e-300)))
}

IMAGE_FILENAMES <- c(
  "tissue_hires_image.png",
  "tissue_hires_image.jpg",
  "tissue_hires_image.jpeg"
)

find_hires_image <- function(sec) {
  direct <- file.path(RAW_ROOT, as.character(sec), "spatial", IMAGE_FILENAMES)
  hit <- direct[file.exists(direct)]
  if (length(hit) >= 1) return(hit[1])

  found <- list.files(
    RAW_ROOT,
    pattern = "^tissue_hires_image\\.(png|jpe?g)$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )

  if (length(found)) {
    section_folder <- basename(dirname(dirname(found)))
    hit <- found[section_folder == as.character(sec)]
  } else {
    hit <- character(0)
  }

  if (length(hit) == 1) return(hit)
  if (length(hit) > 1) {
    stop(
      "Multiple H&E images found for section ", sec, ":\n",
      paste(hit, collapse = "\n")
    )
  }
  NA_character_
}

load_sec_df <- function(sec) {
  path <- file.path(phase06_tables, paste0("zones_section_", sec, ".tsv"))
  if (!file.exists(path)) stop("Missing stage-06 table for section ", sec, ": ", path)

  df <- read_tsv(path, show_col_types = FALSE)
  if ("cell.x" %in% colnames(df) && !"cell" %in% colnames(df)) {
    df <- df %>% rename(cell = cell.x)
  }

  required <- c(
    "section_id", "in_tissue", "is_hotspot_semla", "zone_call",
    "pxl_col_in_hires", "pxl_row_in_hires",
    "frac_epithelial_nbr", "frac_mesenchymal_nbr",
    "frac_vascular_nbr", "frac_myeloid_nbr"
  )
  missing <- setdiff(required, colnames(df))
  if (length(missing)) {
    stop(
      "Stage-06 table for section ", sec, " is missing: ",
      paste(missing, collapse = ", ")
    )
  }

  df %>%
    mutate(
      section_id       = as.character(section_id),
      in_tissue        = as.integer(in_tissue),
      is_hotspot_semla = ifelse(
        is.na(is_hotspot_semla), FALSE, as.logical(is_hotspot_semla)
      ),
      zone_call = factor(as.character(zone_call), levels = zone_levels)
    )
}

prepare_plot_df <- function(df) {
  sec <- unique(as.character(df$section_id))
  if (length(sec) != 1) stop("prepare_plot_df() requires exactly one section.")

  img_path <- find_hires_image(sec)
  if (is.na(img_path) || !file.exists(img_path)) {
    stop("Missing H&E image for section ", sec, " under RAW_ROOT: ", RAW_ROOT)
  }

  img   <- png::readPNG(img_path)
  img_h <- nrow(img)
  img_w <- ncol(img)
  d <- df %>%
    mutate(
      x_plot = pxl_col_in_hires,
      y_plot = img_h - pxl_row_in_hires
    )

  list(df = d, img = img, img_h = img_h, img_w = img_w, img_path = img_path)
}

make_he_canvas <- function(prep, title_txt) {
  g   <- grid::rasterGrob(prep$img, interpolate = TRUE)
  pad <- ceiling(max(prep$img_h, prep$img_w) * HE_PAD_FRAC)
  ggplot() +
    annotation_custom(g, xmin = 0, xmax = prep$img_w,
                      ymin = 0, ymax = prep$img_h) +
    coord_fixed(xlim = c(-pad, prep$img_w + pad),
                ylim = c(-pad, prep$img_h + pad), expand = FALSE) +
    theme_void(base_size = 11) +
    ggtitle(title_txt) +
    theme(plot.title   = element_text(face = "bold", hjust = 0.5, size = 12),
          legend.title = element_text(size = 10, face = "bold"),
          legend.text  = element_text(size = 9))
}

plot_barrier_map <- function(sec, frac_col, high_col, legend_title,
                              shared_max = NULL) {
  df <- load_sec_df(sec); prep <- prepare_plot_df(df); d <- prep$df
  p  <- make_he_canvas(prep, sec)
  d1 <- d %>% filter(in_tissue == 1, is.finite(x_plot), is.finite(y_plot),
                     is.finite(.data[[frac_col]]))
  if (nrow(d1) > 0) {
    max_sc <- if (!is.null(shared_max)) shared_max
              else quantile(d1[[frac_col]], 0.99, na.rm = TRUE)
    if (!is.finite(max_sc) || max_sc <= 0) max_sc <- 1
    p <- p +
      geom_point(data = d1, aes(x_plot, y_plot, colour = .data[[frac_col]]),
                 size = PT_SIZE_TISSUE, alpha = 0.50) +
      scale_colour_gradient(low = "grey90", high = high_col,
                            limits = c(0, max_sc), oob = scales::squish,
                            name = legend_title)
    dh <- d1 %>% filter(is_hotspot_semla)
    if (nrow(dh) > 0)
      p <- p +
        geom_point(data = dh, aes(x_plot, y_plot),
                   shape = 21, fill = NA, color = "black",
                   stroke = HOT_RING_STROKE, size = HOT_RING_SIZE,
                   inherit.aes = FALSE) +
        geom_point(data = dh, aes(x_plot, y_plot),
                   shape = 21, fill = NA, color = "white",
                   stroke = HOT_RING_STROKE - 0.20, size = HOT_RING_SIZE - 0.10,
                   inherit.aes = FALSE)
  }
  p
}

perm_test_barrier <- function(sec_df, frac_col,
                               direction = c("higher", "lower"),
                               n_perm = N_PERM) {
  direction <- match.arg(direction)
  df    <- sec_df %>% filter(in_tissue == 1, is.finite(.data[[frac_col]]))
  n_hot <- sum(df$is_hotspot_semla, na.rm = TRUE)
  if (n_hot < 5) return(NULL)
  obs   <- mean(df[[frac_col]][df$is_hotspot_semla %in% TRUE], na.rm = TRUE)
  nulls <- vapply(seq_len(n_perm), function(i)
    mean(df[[frac_col]][sample(nrow(df), n_hot)], na.rm = TRUE), numeric(1))
  p_val <- if (direction == "higher") (sum(nulls >= obs) + 1) / (n_perm + 1)
           else                       (sum(nulls <= obs) + 1) / (n_perm + 1)
  tibble(obs = obs, null_mean = mean(nulls), p = p_val)
}

make_quant_panel <- function(tab, obs_col, null_col, p_col, bar_color, y_label) {
  plot_tab <- tab %>%
    mutate(section_id = factor(section_id, levels = good_sections),
           x_base = seq_along(good_sections),
           x_obs  = x_base - 0.20, x_null = x_base + 0.20)
  bars <- bind_rows(
    plot_tab %>% transmute(section_id, x = x_obs, x_base,
                           group = "Observed (hotspot neighborhoods)",
                           value = .data[[obs_col]]),
    plot_tab %>% transmute(section_id, x = x_null, x_base,
                           group = "Null (permuted hotspots)",
                           value = .data[[null_col]])
  ) %>%
    mutate(group = factor(group, levels = c("Observed (hotspot neighborhoods)",
                                            "Null (permuted hotspots)")))
  ymax <- max(c(plot_tab[[obs_col]], plot_tab[[null_col]]), na.rm = TRUE)
  if (!is.finite(ymax) || ymax <= 0) ymax <- 1
  br <- plot_tab %>%
    transmute(section_id, x1 = x_obs, x2 = x_null, xm = x_base,
              y  = pmax(.data[[obs_col]], .data[[null_col]], na.rm = TRUE) + 0.06 * ymax,
              y0 = pmax(.data[[obs_col]], .data[[null_col]], na.rm = TRUE) + 0.03 * ymax,
              label = paste0(fmt_p(.data[[p_col]]), " ", stars(.data[[p_col]])))
  ggplot(bars, aes(x = x, y = value, fill = group)) +
    geom_col(width = 0.32, color = NA) +
    geom_segment(data = br, aes(x = x1, xend = x1, y = y0, yend = y),
                 inherit.aes = FALSE, linewidth = 0.7) +
    geom_segment(data = br, aes(x = x2, xend = x2, y = y0, yend = y),
                 inherit.aes = FALSE, linewidth = 0.7) +
    geom_segment(data = br, aes(x = x1, xend = x2, y = y, yend = y),
                 inherit.aes = FALSE, linewidth = 0.7) +
    geom_text(data = br, aes(x = xm, y = y + 0.015 * ymax, label = label),
              inherit.aes = FALSE, vjust = 0, size = 4) +
    scale_x_continuous(breaks = plot_tab$x_base, labels = good_sections) +
    scale_fill_manual(values = c("Observed (hotspot neighborhoods)" = bar_color,
                                 "Null (permuted hotspots)" = alpha(bar_color, 0.25)),
                      name = NULL) +
    labs(x = "Section", y = y_label) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    theme_classic(base_size = 12) +
    theme(axis.text.x  = element_text(angle = 25, hjust = 1),
          axis.title.y = element_text(face = "bold"),
          legend.position = "right")
}

# =========================================================
# Load or recompute permutation results
# =========================================================
perm_tsv <- file.path(fig3_out_dir, "Figure3_all_permutation_results.tsv")

if (file.exists(perm_tsv)) {
  message("Loading permutation results from figure3.R output: ", perm_tsv)
  perm_tab <- read_tsv(perm_tsv, show_col_types = FALSE) %>%
    mutate(section_id = factor(section_id, levels = good_sections))
} else {
  message("figure3.R permutation TSV not found — recomputing four zone families...")
  perm_results <- list()

  for (sec in good_sections) {
    df_sec  <- load_sec_df(sec)
    res_epi <- perm_test_barrier(df_sec, "frac_epithelial_nbr", direction = "lower")
    res_my  <- perm_test_barrier(df_sec, "frac_myeloid_nbr", direction = "higher")
    res_mes <- perm_test_barrier(df_sec, "frac_mesenchymal_nbr", direction = "higher")
    res_vasc <- perm_test_barrier(df_sec, "frac_vascular_nbr", direction = "higher")

    perm_results[[sec]] <- tibble(
      section_id                = sec,
      obs_mean_epi_nbr          = if (!is.null(res_epi)) res_epi$obs else NA_real_,
      null_mean_epi_nbr         = if (!is.null(res_epi)) res_epi$null_mean else NA_real_,
      p_less_epi                = if (!is.null(res_epi)) res_epi$p else NA_real_,
      obs_mean_myeloid_nbr      = if (!is.null(res_my)) res_my$obs else NA_real_,
      null_mean_myeloid_nbr     = if (!is.null(res_my)) res_my$null_mean else NA_real_,
      p_more_myeloid            = if (!is.null(res_my)) res_my$p else NA_real_,
      obs_mean_mesenchymal_nbr  = if (!is.null(res_mes)) res_mes$obs else NA_real_,
      null_mean_mesenchymal_nbr = if (!is.null(res_mes)) res_mes$null_mean else NA_real_,
      p_more_mesenchymal        = if (!is.null(res_mes)) res_mes$p else NA_real_,
      obs_mean_vascular_nbr     = if (!is.null(res_vasc)) res_vasc$obs else NA_real_,
      null_mean_vascular_nbr    = if (!is.null(res_vasc)) res_vasc$null_mean else NA_real_,
      p_more_vascular           = if (!is.null(res_vasc)) res_vasc$p else NA_real_
    )
  }

  perm_tab <- bind_rows(perm_results) %>%
    mutate(section_id = factor(section_id, levels = good_sections)) %>%
    arrange(section_id)

  write_tsv(perm_tab, file.path(out_dir, "S1_permutation_results.tsv"))
}

pdf_device <- if (isTRUE(capabilities("cairo"))) cairo_pdf else "pdf"

# =========================================================
# PANELS A + B — Epithelial
# =========================================================
message("Panels a/b — Epithelial")
figS1A <- wrap_plots(
  lapply(example_sections, function(s)
    plot_barrier_map(s, "frac_epithelial_nbr", unname(zone_cols["Epithelial"]),
                     "Mean fraction of\nEpithelial neighbors")),
  nrow = 1, guides = "collect"
) & theme(legend.position = "right")

ggsave(file.path(out_dir, "S1a_epithelial_maps_459_928.png"),
       figS1A, width = 8.8, height = 4.4, dpi = PLOT$dpi, bg = "white")
ggsave(file.path(out_dir, "S1a_epithelial_maps_459_928.pdf"),
       figS1A, width = 8.8, height = 4.4, device = pdf_device, bg = "white")

pS1B <- make_quant_panel(perm_tab, "obs_mean_epi_nbr", "null_mean_epi_nbr",
                          "p_less_epi", unname(zone_cols["Epithelial"]),
                          "Mean fraction of Epithelial neighbors (6-NN)")

ggsave(file.path(out_dir, "S1b_epithelial_permutation.png"),
       pS1B, width = 9.0, height = 4.8, dpi = PLOT$dpi, bg = "white")
ggsave(file.path(out_dir, "S1b_epithelial_permutation.pdf"),
       pS1B, width = 9.0, height = 4.8, device = pdf_device, bg = "white")

# =========================================================
# PANELS C + D — Vascular
# =========================================================
message("Panels c/d — Vascular")

has_vasc <- "p_more_vascular" %in% colnames(perm_tab) &&
            any(!is.na(perm_tab$p_more_vascular))

if (!has_vasc) {
  message("WARNING: no vascular permutation results available. Skipping panels c/d.")
} else {
  # Shared scale across example sections
  vasc_max <- quantile(
    unlist(lapply(example_sections, function(s) {
      d <- load_sec_df(s)
      d$frac_vascular_nbr[d$in_tissue == 1 & is.finite(d$frac_vascular_nbr)]
    })), 0.99, na.rm = TRUE)
  if (!is.finite(vasc_max) || vasc_max <= 0) vasc_max <- 1

  figS1C <- wrap_plots(
    lapply(example_sections, function(s)
      plot_barrier_map(s, "frac_vascular_nbr", unname(zone_cols["Vascular"]),
                       "Mean fraction\nVascular neighbors", shared_max = vasc_max)),
    nrow = 1, guides = "collect"
  ) & theme(legend.position = "right")

  ggsave(file.path(out_dir, "S1c_vascular_maps_459_928.png"),
         figS1C, width = 8.8, height = 4.4, dpi = PLOT$dpi, bg = "white")
  ggsave(file.path(out_dir, "S1c_vascular_maps_459_928.pdf"),
         figS1C, width = 8.8, height = 4.4, device = pdf_device, bg = "white")

  pS1D <- make_quant_panel(perm_tab, "obs_mean_vascular_nbr", "null_mean_vascular_nbr",
                            "p_more_vascular", unname(zone_cols["Vascular"]),
                            "Mean fraction of Vascular neighbors (6-NN)")

  ggsave(file.path(out_dir, "S1d_vascular_permutation.png"),
         pS1D, width = 9.0, height = 4.8, dpi = PLOT$dpi, bg = "white")
  ggsave(file.path(out_dir, "S1d_vascular_permutation.pdf"),
         pS1D, width = 9.0, height = 4.8, device = pdf_device, bg = "white")
}

# =========================================================
# PANEL E — Section 723: all four zones in one row
# =========================================================
message("Panel e — Section 723")

SEC_723 <- "723"
df_723  <- load_sec_df(SEC_723)
prep_723 <- prepare_plot_df(df_723)

img_723   <- prep_723$img
img_h_723 <- prep_723$img_h
img_w_723 <- prep_723$img_w

lighten_png <- function(img, alpha = 0.55) {
  out <- img
  for (k in seq_len(min(3, dim(out)[3])))
    out[,,k] <- 1 - alpha * (1 - out[,,k])
  out
}
img_lite_723 <- lighten_png(img_723)
g_he_723     <- grid::rasterGrob(img_lite_723, interpolate = TRUE)
pad_723      <- ceiling(max(img_h_723, img_w_723) * HE_PAD_FRAC)

df_723 <- df_723 %>%
  mutate(x_plot = pxl_col_in_hires, y_plot = img_h_723 - pxl_row_in_hires)

d1_723 <- df_723 %>% filter(in_tissue == 1, is.finite(x_plot), is.finite(y_plot))
dh_723 <- d1_723 %>% filter(is_hotspot_semla)

zone_fracs_723 <- c(
  "Epithelial"  = "frac_epithelial_nbr",
  "Mesenchymal" = "frac_mesenchymal_nbr",
  "Vascular"    = "frac_vascular_nbr",
  "Myeloid"     = "frac_myeloid_nbr"
)

plot_zone_nbr_723 <- function(zone_name) {
  frac_col <- zone_fracs_723[zone_name]
  col_high <- zone_cols[zone_name]
  d1       <- d1_723 %>% filter(is.finite(.data[[frac_col]]))
  lim_max  <- quantile(d1[[frac_col]], 0.99, na.rm = TRUE)
  if (!is.finite(lim_max) || lim_max <= 0) lim_max <- 1

  p <- ggplot() +
    annotation_custom(g_he_723, xmin = 0, xmax = img_w_723,
                      ymin = 0, ymax = img_h_723) +
    coord_fixed(xlim = c(-pad_723, img_w_723 + pad_723),
                ylim = c(-pad_723, img_h_723 + pad_723), expand = FALSE) +
    theme_void(base_size = 11) +
    labs(title = paste0("Mean fraction ", zone_name, " neighbours")) +
    theme(plot.title   = element_text(face = "bold", hjust = 0.5, size = 10),
          legend.title = element_text(size = 9, face = "bold"),
          legend.text  = element_text(size = 8))

  if (nrow(d1) > 0)
    p <- p +
      geom_point(data = d1, aes(x_plot, y_plot, colour = .data[[frac_col]]),
                 size = PT_SIZE_TISSUE, alpha = 0.75, inherit.aes = FALSE) +
      scale_colour_gradient(low = "grey95", high = col_high,
                            limits = c(0, lim_max), oob = scales::squish,
                            name = "Mean\nfraction\n(6-NN)")

  if (nrow(dh_723) > 0)
    p <- p +
      geom_point(data = dh_723, aes(x_plot, y_plot),
                 shape = 21, fill = NA, color = "black",
                 stroke = 0.50, size = 0.70, inherit.aes = FALSE) +
      geom_point(data = dh_723, aes(x_plot, y_plot),
                 shape = 21, fill = NA, color = "white",
                 stroke = 0.30, size = 0.60, inherit.aes = FALSE)
  p
}

# Retrieve the same raw permutation p-values shown by the original S1 panel.
p_723 <- perm_tab %>%
  filter(as.character(section_id) == SEC_723) %>%
  slice(1)

get_p723 <- function(col) {
  if (!col %in% names(p_723) || nrow(p_723) == 0) return(NA_real_)
  as.numeric(p_723[[col]][1])
}

epi_p_723  <- get_p723("p_less_epi")
my_p_723   <- get_p723("p_more_myeloid")
mes_p_723  <- get_p723("p_more_mesenchymal")
vasc_p_723 <- get_p723("p_more_vascular")

p_text <- function(p) {
  if (is.na(p)) "p=NA" else paste0(fmt_p(p), " ", stars(p))
}

figS1E <- wrap_plots(lapply(names(zone_fracs_723), plot_zone_nbr_723), nrow = 1) +
  plot_annotation(
    title    = "Section 723 — neighbourhood zone fractions",
    subtitle = paste0(
      "Myeloid: ", p_text(my_p_723),
      "  |  Mesenchymal: ", p_text(mes_p_723),
      "  |  Epithelial: ", p_text(epi_p_723),
      "  |  Vascular: ", p_text(vasc_p_723)
    ),
    theme = theme(
      plot.title    = element_text(face = "bold", hjust = 0.5, size = 13),
      plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey40")
    )
  )

ggsave(file.path(out_dir, "S1e_section723_neighbourhood_4panels.png"),
       figS1E, width = 16, height = 5, dpi = PLOT$dpi, bg = "white")
ggsave(file.path(out_dir, "S1e_section723_neighbourhood_4panels.pdf"),
       figS1E, width = 16, height = 5, device = pdf_device, bg = "white")

message("Supplementary Figure 1 done. Outputs saved to: ", out_dir)
