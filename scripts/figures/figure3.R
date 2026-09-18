# =========================================================
# Figure 3 — unified script
#
# Panels:
#   a  H&E / zone maps / hotspot overlay (sections 459 & 928)
#   b  Myeloid neighbourhood spatial maps (459 & 928)
#   c  Myeloid neighbourhood permutation bar chart
#   d  Mesenchymal neighbourhood spatial maps (459 & 928)
#   e  Mesenchymal neighbourhood permutation bar chart
#   f  Myeloid–Mesenchymal interface spatial maps (459 & 928)
#   g  Interface permutation bar chart
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
out_dir <- ensure_dir(figure_path("figure3"))

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
PT_SIZE_OFF     <- 0.55
HOT_RING_SIZE   <- 1.00
HOT_RING_STROKE <- 0.80
HE_PAD_FRAC     <- 0.015
N_PERM          <- 1000
use_seed("neighborhood_permutation")

# =========================================================
# 2. Helpers
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

load_sec_df <- function(sec) {
  path <- file.path(phase06_tables, paste0("zones_section_", sec, ".tsv"))
  if (!file.exists(path)) stop("Missing phase06 table for section ", sec, ": ", path)

  df <- read_tsv(path, show_col_types = FALSE)
  if ("cell.x" %in% colnames(df) && !"cell" %in% colnames(df)) {
    df <- df %>% rename(cell = cell.x)
  }

  # Do not trust absolute image paths embedded in old intermediate tables.
  # Reconstruct the H&E location from the configured RAW_ROOT instead.
  df %>%
    mutate(
      section_id       = as.character(section_id),
      in_tissue        = as.integer(in_tissue),
      is_hotspot_semla = ifelse(
        is.na(is_hotspot_semla), FALSE, as.logical(is_hotspot_semla)
      ),
      zone_call        = factor(as.character(zone_call), levels = zone_levels),
      img_file         = file.path(
        RAW_ROOT, as.character(section_id), "spatial", "tissue_hires_image.png"
      )
    )
}

zones_all <- read_tsv(file.path(phase06_tables, "zones_all_sections.tsv"),
                      show_col_types = FALSE) %>%
  mutate(
    section_id       = as.character(section_id),
    in_tissue        = as.integer(in_tissue),
    is_hotspot_semla = ifelse(is.na(is_hotspot_semla), FALSE,
                              as.logical(is_hotspot_semla)),
    zone_call        = factor(as.character(zone_call), levels = zone_levels)
  ) %>%
  filter(section_id %in% good_sections)

prepare_plot_df <- function(df) {
  img_path <- unique(df$img_file)
  if (length(img_path) != 1 || is.na(img_path) || !file.exists(img_path))
    stop("Missing/invalid img_file for section ", unique(df$section_id)[1])
  img   <- png::readPNG(img_path)
  img_h <- nrow(img); img_w <- ncol(img)
  d <- df %>% mutate(x_plot = pxl_col_in_hires, y_plot = img_h - pxl_row_in_hires)
  list(df = d, img = img, img_h = img_h, img_w = img_w)
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

plot_he_only <- function(sec) {
  df <- load_sec_df(sec); prep <- prepare_plot_df(df); d <- prep$df
  p  <- make_he_canvas(prep, sec)
  d0 <- d %>% filter(in_tissue == 0, is.finite(x_plot), is.finite(y_plot))
  if (nrow(d0) > 0)
    p <- p + geom_point(data = d0, aes(x_plot, y_plot),
                        shape = 1, stroke = 0.35, color = "grey65",
                        size = PT_SIZE_OFF, inherit.aes = FALSE)
  p
}

plot_zones <- function(sec, show_hotspots = FALSE) {
  df <- load_sec_df(sec); prep <- prepare_plot_df(df); d <- prep$df
  p  <- make_he_canvas(prep, sec)
  d0 <- d %>% filter(in_tissue == 0, is.finite(x_plot), is.finite(y_plot))
  d1 <- d %>% filter(in_tissue == 1, is.finite(x_plot), is.finite(y_plot))
  if (nrow(d0) > 0)
    p <- p + geom_point(data = d0, aes(x_plot, y_plot),
                        shape = 1, stroke = 0.35, color = "grey65",
                        size = PT_SIZE_OFF, inherit.aes = FALSE)
  if (nrow(d1) > 0)
    p <- p +
    geom_point(data = d1, aes(x_plot, y_plot, colour = zone_call),
               size = PT_SIZE_TISSUE, alpha = 0.90) +
    scale_colour_manual(values = zone_cols, drop = FALSE,
                        na.value = "grey80", name = "zone")
  if (show_hotspots) {
    dh <- d1 %>% filter(is_hotspot_semla)
    if (nrow(dh) > 0)
      p <- p + geom_point(data = dh, aes(x_plot, y_plot),
                          shape = 21, fill = NA, color = "white",
                          stroke = HOT_RING_STROKE, size = HOT_RING_SIZE,
                          inherit.aes = FALSE)
  }
  p
}

plot_barrier_map <- function(sec, frac_col, high_col, legend_title) {
  df <- load_sec_df(sec); prep <- prepare_plot_df(df); d <- prep$df
  p  <- make_he_canvas(prep, sec)
  d1 <- d %>% filter(in_tissue == 1, is.finite(x_plot), is.finite(y_plot),
                     is.finite(.data[[frac_col]]))
  if (nrow(d1) > 0) {
    max_sc <- quantile(d1[[frac_col]], 0.99, na.rm = TRUE)
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

perm_test_interface <- function(sec_df, n_perm = N_PERM) {
  df <- sec_df %>%
    filter(in_tissue == 1,
           is.finite(frac_myeloid_nbr), is.finite(frac_mesenchymal_nbr)) %>%
    mutate(is_interface = frac_myeloid_nbr > 0 & frac_mesenchymal_nbr > 0)
  n_hot <- sum(df$is_hotspot_semla, na.rm = TRUE)
  if (n_hot < 5) return(NULL)
  obs   <- mean(df$is_interface[df$is_hotspot_semla], na.rm = TRUE)
  nulls <- vapply(seq_len(n_perm), function(i)
    mean(df$is_interface[sample(nrow(df), n_hot)], na.rm = TRUE), numeric(1))
  p_val <- (sum(nulls >= obs) + 1) / (n_perm + 1)
  tibble(obs = obs, null_mean = mean(nulls), p = p_val,
         n_interface_hotspot = sum(df$is_interface[df$is_hotspot_semla]),
         n_hotspot = n_hot)
}

make_quant_panel <- function(tab, obs_col, null_col, p_col, bar_color, y_label) {
  plot_tab <- tab %>%
    mutate(section_id = factor(section_id, levels = good_sections),
           x_base = seq_along(good_sections),
           x_obs  = x_base - 0.20,
           x_null = x_base + 0.20)
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
    geom_segment(data = br, aes(x = x1, xend = x2, y = y,  yend = y),
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
# PANEL A
# =========================================================
message("Panel a")
fig3A <- (
  plot_he_only("459") +
    plot_zones("459", show_hotspots = FALSE) + theme(legend.position = "none") +
    plot_zones("459", show_hotspots = TRUE)
) / (
  plot_he_only("928") +
    plot_zones("928", show_hotspots = FALSE) + theme(legend.position = "none") +
    plot_zones("928", show_hotspots = TRUE)
) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right", legend.direction = "vertical")

ggsave(file.path(out_dir, "Figure3a_HE_zones_hotspots_459_928.png"),
       fig3A, width = 13.5, height = 8.8, dpi = 600, bg = "white")
ggsave(file.path(out_dir, "Figure3a_HE_zones_hotspots_459_928.pdf"),
       fig3A, width = 13.5, height = 8.8, bg = "white")

# =========================================================
# PERMUTATION RESULTS — all zones in one loop
# (Myeloid + Mesenchymal for main fig; Epithelial + Vascular saved for S1)
# =========================================================
message("Running permutations...")

perm_results <- list()
for (sec in good_sections) {
  df_sec <- load_sec_df(sec)
  res_epi  <- perm_test_barrier(df_sec, "frac_epithelial_nbr",  direction = "lower")
  res_my   <- perm_test_barrier(df_sec, "frac_myeloid_nbr",     direction = "higher")
  res_mes  <- perm_test_barrier(df_sec, "frac_mesenchymal_nbr", direction = "higher")
  res_vasc <- if ("frac_vascular_nbr" %in% colnames(df_sec))
    perm_test_barrier(df_sec, "frac_vascular_nbr", direction = "higher")
  else NULL
  res_iface <- perm_test_interface(df_sec)
  
  if (is.null(res_epi) || is.null(res_my) || is.null(res_mes)) next
  
  row <- tibble(
    section_id                = sec,
    obs_mean_epi_nbr          = res_epi$obs,
    null_mean_epi_nbr         = res_epi$null_mean,
    p_less_epi                = res_epi$p,
    obs_mean_myeloid_nbr      = res_my$obs,
    null_mean_myeloid_nbr     = res_my$null_mean,
    p_more_myeloid            = res_my$p,
    obs_mean_mesenchymal_nbr  = res_mes$obs,
    null_mean_mesenchymal_nbr = res_mes$null_mean,
    p_more_mesenchymal        = res_mes$p,
    obs_mean_vascular_nbr     = if (!is.null(res_vasc)) res_vasc$obs      else NA_real_,
    null_mean_vascular_nbr    = if (!is.null(res_vasc)) res_vasc$null_mean else NA_real_,
    p_more_vascular           = if (!is.null(res_vasc)) res_vasc$p         else NA_real_,
    obs_interface             = if (!is.null(res_iface)) res_iface$obs      else NA_real_,
    null_mean_interface       = if (!is.null(res_iface)) res_iface$null_mean else NA_real_,
    p_interface               = if (!is.null(res_iface)) res_iface$p         else NA_real_
  )
  perm_results[[sec]] <- row
}

perm_tab <- bind_rows(perm_results) %>%
  mutate(section_id = factor(section_id, levels = good_sections)) %>%
  arrange(section_id)

# BH correction, WITHIN each test family (one zone/interface test, across the
# 6 hotspot-positive sections) rather than pooled across all five families.
# Each p_* column is a distinct biological question (epithelial depletion,
# myeloid enrichment, mesenchymal enrichment, vascular enrichment, interface
# enrichment); correcting them together would be over-conservative for some
# and under-conservative for others. NAs (sections excluded for n_hot < 5, or
# the vascular test when unavailable) are excluded from each correction and
# preserved as NA in the adjusted column.
p_cols <- c("p_less_epi", "p_more_myeloid", "p_more_mesenchymal",
            "p_more_vascular", "p_interface")
p_cols <- intersect(p_cols, names(perm_tab))

for (pc in p_cols) {
  adj_col <- paste0(pc, "_adj")
  vals <- perm_tab[[pc]]
  adj  <- rep(NA_real_, length(vals))
  ok   <- !is.na(vals)
  adj[ok] <- p.adjust(vals[ok], method = "BH")
  perm_tab[[adj_col]] <- adj
}

write_tsv(perm_tab, file.path(out_dir, "Figure3_all_permutation_results.tsv"))

# =========================================================
# PANELS B + C — Myeloid
# =========================================================
message("Panels b/c")
fig3C <- wrap_plots(
  lapply(example_sections, function(s)
    plot_barrier_map(s, "frac_myeloid_nbr", unname(zone_cols["Myeloid"]),
                     "Mean fraction of\nMyeloid neighbors")),
  nrow = 1, guides = "collect"
) & theme(legend.position = "right")

ggsave(file.path(out_dir, "Figure3b_myeloid_maps_459_928.png"),
       fig3C, width = 8.8, height = 4.4, dpi = 600, bg = "white")
ggsave(file.path(out_dir, "Figure3b_myeloid_maps_459_928.pdf"),
       fig3C, width = 8.8, height = 4.4, bg = "white")

pC <- make_quant_panel(perm_tab, "obs_mean_myeloid_nbr", "null_mean_myeloid_nbr",
                       "p_more_myeloid_adj", unname(zone_cols["Myeloid"]),
                       "Mean fraction of Myeloid neighbors (6-NN)")

ggsave(file.path(out_dir, "Figure3c_myeloid_permutation.png"),
       pC, width = 9.0, height = 4.8, dpi = 600, bg = "white")
ggsave(file.path(out_dir, "Figure3c_myeloid_permutation.pdf"),
       pC, width = 9.0, height = 4.8, bg = "white")

# =========================================================
# PANELS D + E — Mesenchymal
# =========================================================
message("Panels d/e")
fig3E <- wrap_plots(
  lapply(example_sections, function(s)
    plot_barrier_map(s, "frac_mesenchymal_nbr", unname(zone_cols["Mesenchymal"]),
                     "Mean fraction of\nMesenchymal neighbors")),
  nrow = 1, guides = "collect"
) & theme(legend.position = "right")

ggsave(file.path(out_dir, "Figure3d_mesenchymal_maps_459_928.png"),
       fig3E, width = 8.8, height = 4.4, dpi = 600, bg = "white")
ggsave(file.path(out_dir, "Figure3d_mesenchymal_maps_459_928.pdf"),
       fig3E, width = 8.8, height = 4.4, bg = "white")

pE <- make_quant_panel(perm_tab, "obs_mean_mesenchymal_nbr", "null_mean_mesenchymal_nbr",
                       "p_more_mesenchymal_adj", unname(zone_cols["Mesenchymal"]),
                       "Mean fraction of Mesenchymal neighbors (6-NN)")

ggsave(file.path(out_dir, "Figure3e_mesenchymal_permutation.png"),
       pE, width = 9.0, height = 4.8, dpi = 600, bg = "white")
ggsave(file.path(out_dir, "Figure3e_mesenchymal_permutation.pdf"),
       pE, width = 9.0, height = 4.8, bg = "white")

# =========================================================
# PANELS F + G — Myeloid–Mesenchymal interface
# =========================================================
message("Panels f/g")

plot_interface_map <- function(sec) {
  df       <- load_sec_df(sec)
  img_path <- unique(df$img_file)
  use_he   <- length(img_path) == 1 && !is.na(img_path) && file.exists(img_path)
  
  if (use_he) {
    prep <- list()
    prep$img  <- png::readPNG(img_path)
    prep$img_h <- nrow(prep$img); prep$img_w <- ncol(prep$img)
    df <- df %>% mutate(x_plot = pxl_col_in_hires,
                        y_plot = prep$img_h - pxl_row_in_hires)
    d  <- df
  } else {
    message("  img_file not found for section ", sec, " — coord-only fallback")
    df <- df %>% mutate(x_plot = pxl_col_in_hires, y_plot = -pxl_row_in_hires)
    d  <- df
  }
  
  d1 <- d %>%
    filter(in_tissue == 1, is.finite(x_plot), is.finite(y_plot),
           is.finite(frac_myeloid_nbr), is.finite(frac_mesenchymal_nbr)) %>%
    mutate(interface_score = frac_myeloid_nbr * frac_mesenchymal_nbr)
  
  if (use_he) {
    prep_obj <- list(img = prep$img, img_h = prep$img_h, img_w = prep$img_w)
    p <- make_he_canvas(prep_obj, sec)
  } else {
    p <- ggplot() + coord_equal() + theme_void(base_size = 11) + ggtitle(sec) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12))
  }
  
  if (nrow(d1) > 0) {
    max_sc <- quantile(d1$interface_score, 0.99, na.rm = TRUE)
    if (!is.finite(max_sc) || max_sc <= 0) max_sc <- 1
    p <- p +
      geom_point(data = d1, aes(x_plot, y_plot, colour = interface_score),
                 size = PT_SIZE_TISSUE, alpha = 0.70, inherit.aes = FALSE) +
      scale_colour_gradient(low = "grey90", high = "#8B4513",
                            limits = c(0, max_sc), oob = scales::squish,
                            name  = "Interface\nscore\n(Mye × Mes)")
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

fig3G <- wrap_plots(lapply(example_sections, plot_interface_map),
                    nrow = 1, guides = "collect") &
  theme(legend.position = "right")

ggsave(file.path(out_dir, "Figure3f_interface_maps_459_928.png"),
       fig3G, width = 8.8, height = 4.4, dpi = 600, bg = "white")
ggsave(file.path(out_dir, "Figure3f_interface_maps_459_928.pdf"),
       fig3G, width = 8.8, height = 4.4, bg = "white")

pG <- make_quant_panel(perm_tab, "obs_interface", "null_mean_interface",
                       "p_interface_adj", "#8B4513",
                       "Proportion of spots with both Myeloid & Mesenchymal neighbours")

ggsave(file.path(out_dir, "Figure3g_interface_permutation.png"),
       pG, width = 9.0, height = 4.8, dpi = 600, bg = "white")
ggsave(file.path(out_dir, "Figure3g_interface_permutation.pdf"),
       pG, width = 9.0, height = 4.8, bg = "white")

message("Figure 3 done. Outputs saved to: ", out_dir)
message("Permutation TSV saved — use it as input for S1.R.")