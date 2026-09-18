# =============================================================================
# Supplementary Figure S2 — LR scoring QC (4 panels)
#
# Panel a: Fraction significant LR pairs (p < 0.05) — existing PNG from phase 08
# Panel b: Mean observed LR score — existing PNG from phase 08
# Panel c: Full top-30 LR pair dot plot — CellChatDB
# Panel d: Full top-30 LR pair dot plot — CellPhoneDB
#
# Panels a+b read from: 08_lr_interactions/qc/ (existing PNGs)
# Panels c+d read from: 08_lr_interactions/top_pairs/{resource}_{zone}_top30.tsv
# Output: figureS2/FigureS2_LR_QC_4panel.pdf / .png
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(cowplot)
  library(png)
  library(grid)
})

# =============================================================================
# 0. Paths
# =============================================================================
analysis_root <- "D:/Ped-CNS_KBH/Spatial Transcriptomic/GSE195661/Analysis"
qc_dir        <- file.path(analysis_root, "08_lr_interactions", "qc_cross_section")
top_dir       <- file.path(analysis_root, "08_lr_interactions", "top_pairs")
out_dir       <- "D:/Ped-CNS_KBH/Manuscript_2_EPN/S2"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. Settings
# =============================================================================
resources   <- c("CellChatDB", "CellPhoneDB")
zones       <- c("Mesenchymal", "Myeloid", "Vascular")
zone_levels <- c("Mesenchymal", "Myeloid", "Vascular")
n_top_pairs <- 30

db_cols <- c(
  "CellChatDB"  = "#E8846A",
  "CellPhoneDB" = "#45B8AC"
)

# =============================================================================
# 2. Panels a + b — read existing PNGs from phase 08 QC output
# =============================================================================
read_png_panel <- function(path, label) {
  img  <- readPNG(path)
  grob <- rasterGrob(img, interpolate = TRUE)
  ggdraw() +
    draw_grob(grob) +
    draw_label(label, x = 0.01, y = 0.98,
               hjust = 0, vjust = 1,
               fontface = "bold", size = 14)
}

panel_a <- read_png_panel(
  file.path(qc_dir, "QC_fraction_significant_pairs_by_section.png"), "a"
)
panel_b <- read_png_panel(
  file.path(qc_dir, "QC_mean_observed_score_by_section.png"), "b"
)

# =============================================================================
# 3. Panels c + d — generate full top-30 dot plots per database
# =============================================================================
load_top_pairs <- function(resource) {
  top_list <- lapply(zones, function(z) {
    fname <- file.path(top_dir,
                       paste0(resource, "_", z, "_top", n_top_pairs, ".tsv"))
    if (!file.exists(fname)) {
      message("Missing: ", fname); return(NULL)
    }
    read_tsv(fname, show_col_types = FALSE) %>%
      mutate(
        resource   = resource,
        zone       = z,
        pair_label = paste0(toupper(ligand), "\u2013", toupper(receptor)),
        score_plot = if (z == "Vascular") max_score else median_score
      )
  })
  bind_rows(top_list)
}

make_pair_order <- function(df) {
  df %>%
    group_by(pair_label) %>%
    summarise(best = max(score_plot, na.rm = TRUE), .groups = "drop") %>%
    arrange(best) %>%
    pull(pair_label)
}

plot_dotplot <- function(resource) {
  df         <- load_top_pairs(resource)
  pair_order <- make_pair_order(df)
  fill_col   <- db_cols[resource]
  
  df <- df %>%
    mutate(
      zone       = factor(zone, levels = zone_levels),
      pair_label = factor(pair_label, levels = pair_order)
    )
  
  ggplot(df, aes(x = zone, y = pair_label)) +
    geom_point(
      aes(size = n_sections, fill = score_plot),
      shape = 21, colour = "black", stroke = 0.25
    ) +
    scale_fill_gradient(
      low  = "white",
      high = fill_col,
      name = "Score\n(median/max)"
    ) +
    scale_size_continuous(
      range  = c(2, 8),
      breaks = c(1, 2, 3, 4, 5),
      name   = "N sections"
    ) +
    labs(
      title = resource,
      x     = NULL,
      y     = "Ligand\u2013receptor pair"
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title       = element_text(face = "bold", size = 11,
                                      colour = fill_col),
      axis.text.y      = element_text(size = 7),
      axis.text.x      = element_text(face = "bold", size = 10),
      panel.grid.major = element_line(colour = "grey92", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      legend.box       = "vertical",
      legend.key.size  = unit(0.4, "cm")
    )
}

panel_c_gg <- plot_dotplot("CellChatDB")
panel_d_gg <- plot_dotplot("CellPhoneDB")

# Wrap ggplot panels with bold panel labels
panel_c <- ggdraw() +
  draw_plot(panel_c_gg) +
  draw_label("c", x = 0.01, y = 0.98,
             hjust = 0, vjust = 1,
             fontface = "bold", size = 14)

panel_d <- ggdraw() +
  draw_plot(panel_d_gg) +
  draw_label("d", x = 0.01, y = 0.98,
             hjust = 0, vjust = 1,
             fontface = "bold", size = 14)

# =============================================================================
# 4. Assemble 2x2 layout
# Top row: a + b (QC bar charts, shorter)
# Bottom row: c + d (dot plots, taller to fit pair labels)
# =============================================================================
top_row    <- plot_grid(panel_a, panel_b, ncol = 2, rel_widths = c(1, 1))
bottom_row <- plot_grid(panel_c, panel_d, ncol = 2, rel_widths = c(1, 1))

p_combined <- plot_grid(
  top_row,
  bottom_row,
  ncol        = 1,
  rel_heights = c(0.9, 1.8)   # dot plots get more vertical space
)

# =============================================================================
# 5. Save
# =============================================================================
ggsave(
  filename = file.path(out_dir, "FigureS2_LR_QC_4panel.pdf"),
  plot     = p_combined,
  width    = 14,
  height   = 18,
  device   = cairo_pdf
)

ggsave(
  filename = file.path(out_dir, "FigureS2_LR_QC_4panel.png"),
  plot     = p_combined,
  width    = 14,
  height   = 18,
  dpi      = 300
)

message("Supplementary Figure S2 (4 panels) saved to: ", out_dir)