# =========================================================
# Figure 4 — unified script
#
# Panels:
#   a  Hotspot functional state spatial maps (sections 459 & 928)
#   b  % hotspot functional state frequency (6 sections, stacked bar)
#   c  Retention gene spatial expression grid (section 459, 6×2)
#   d  GPNMB spatial expression (sections 459 & 928, Fig4C style)
#   e  GPNMB by zone — hotspot vs non-hotspot boxplots
#   f  LR consensus dot plot (Myeloid / Mesenchymal / Vascular)
#   g  NicheNet ligand-target heatmap
#   h  S100A9 & CCL2 spatial maps (section 928)
#
# Panels f-h consume stage 08 ligand-receptor/NicheNet outputs. Panels g/h are
# generated only when the required NicheNet files are present; otherwise the
# script skips them with an informative message.
# =========================================================

source(here::here("config.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(scales)
  library(grid)
  library(png)
  library(Matrix)
  library(FNN)
  library(dbscan)
  library(stringr)
})

# =========================================================
# 0. Paths
# =========================================================
phase07_root     <- stage_dir("functional_states")
phase07_tables   <- file.path(phase07_root, "tables")
phase07_sections <- file.path(phase07_root, "sections")
phase06_tables   <- file.path(stage_dir("zones"), "tables")

out_dir <- ensure_dir(figure_path("figure4"))

# =========================================================
# 1. Settings
# =========================================================
# All hotspot-positive sections — functional-state analysis
good_sections <- GOOD_SECTIONS

# Sections displaying the recurrent hotspot zone-confinement pattern
# — used for downstream spatial niche analyses
niche_sections <- c("459", "812", "821", "928", "1239")

panelA_sections <- c("459", "928")
panelC_section  <- "459"
panelD_sections <- c("459", "928")

state_levels <- c("Cytotoxic", "Inflammatory", "Exhausted", "IL17-like", "Retention", "Unassigned")
state_cols   <- c(
  "Cytotoxic"    = "#D73027",
  "Inflammatory" = "#FC8D59",
  "Exhausted"    = "#4575B4",
  "IL17-like"    = "#7B3294",
  "Retention"    = "#1A9850",
  "Unassigned"   = "grey55"
)

retention_genes <- c(
  "TNC", "ICAM1", "VCAM1", "ITGA3", "PTX3", "CXCL3",
  "SERPINE1", "PPP1R15A", "PDK1", "UPP1", "SLC2A1", "NAMPT",
  "GPNMB"    # panel d also reads from the same scaled file
)

zone_levels <- c("Epithelial", "Mesenchymal", "Vascular", "Myeloid", "Uncertain")
zone_cols   <- c(
  "Epithelial"  = "#4DAF4A",
  "Mesenchymal" = "#E41A1C",
  "Vascular"    = "#377EB8",
  "Myeloid"     = "#984EA3",
  "Uncertain"   = "grey70"
)

# Plotting constants
PT_SIZE_OFF      <- 0.55
PT_SIZE_TISSUE   <- 0.52
PT_SIZE_STATE    <- 0.62
PT_SIZE_EXPR     <- 0.58
HE_PAD_FRAC      <- 0.015
HE_ALPHA_PANELC  <- 0.55
HOT_RING_SIZE    <- 0.70
HOT_RING_STR1    <- 0.50
HOT_RING_STR2    <- 0.30
EXPR_THRESHOLD_QC <- 0.20
N_TOP_PANELD     <- 20

# =========================================================
# 2. Shared helpers
# =========================================================
lighten_png <- function(img, alpha = 1) {
  alpha <- max(0, min(1, alpha))
  out   <- img
  if (length(dim(out)) < 3) return(out)
  for (k in seq_len(min(3, dim(out)[3])))
    out[,,k] <- 1 - alpha * (1 - out[,,k])
  out
}

section_img_path <- function(sec) {
  candidates <- file.path(
    RAW_ROOT,
    as.character(sec),
    "spatial",
    c(
      "tissue_hires_image.png",
      "tissue_hires_image.jpg",
      "tissue_hires_image.jpeg"
    )
  )

  hit <- candidates[file.exists(candidates)]
  if (length(hit) >= 1) return(hit[1])

  # Fallback for a non-standard extracted GEO directory layout.
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

prepare_he_plot <- function(df) {
  sec <- unique(as.character(df$section_id))
  if (length(sec) != 1 || is.na(sec)) stop("Expected one section in prepare_he_plot().")
  img_path <- section_img_path(sec)
  if (is.na(img_path) || !file.exists(img_path))
    stop("Missing H&E image for section ", sec, " under RAW_ROOT: ", RAW_ROOT)
  img   <- png::readPNG(img_path)
  img_h <- nrow(img); img_w <- ncol(img)
  df2   <- df %>% mutate(x_plot = pxl_col_in_hires,
                         y_plot = img_h - pxl_row_in_hires)
  list(df = df2, img = img, img_h = img_h, img_w = img_w)
}

make_he_canvas <- function(prep, title_txt = NULL, he_alpha = 1) {
  img_use <- lighten_png(prep$img, alpha = he_alpha)
  g   <- grid::rasterGrob(img_use, interpolate = TRUE)
  pad <- ceiling(max(prep$img_h, prep$img_w) * HE_PAD_FRAC)
  ggplot() +
    annotation_custom(g, xmin = 0, xmax = prep$img_w,
                      ymin = 0, ymax = prep$img_h) +
    coord_fixed(xlim = c(-pad, prep$img_w + pad),
                ylim = c(-pad, prep$img_h + pad), expand = FALSE) +
    theme_void(base_size = 11) +
    labs(title = title_txt) +
    theme(plot.title   = element_text(face = "bold", hjust = 0.5, size = 12),
          legend.title = element_text(size = 10, face = "bold"),
          legend.text  = element_text(size = 9))
}

add_hotspot_rings <- function(p, dh) {
  if (nrow(dh) == 0) return(p)
  p +
    geom_point(data = dh, aes(x_plot, y_plot),
               shape = 21, fill = NA, color = "black",
               stroke = HOT_RING_STR1, size = HOT_RING_SIZE,
               inherit.aes = FALSE) +
    geom_point(data = dh, aes(x_plot, y_plot),
               shape = 21, fill = NA, color = "white",
               stroke = HOT_RING_STR2, size = HOT_RING_SIZE,
               inherit.aes = FALSE)
}

read_plotready <- function(sec) {
  f_rds <- file.path(phase07_sections, paste0("phase07_section_", sec, "_plotready.rds"))
  f_tsv <- file.path(phase07_sections, paste0("phase07_section_", sec, "_plotready.tsv"))
  if (file.exists(f_rds)) return(readRDS(f_rds))
  if (file.exists(f_tsv)) return(read_tsv(f_tsv, show_col_types = FALSE))
  stop("Missing plotready file for section ", sec)
}

read_retention <- function(sec) {
  f_rds <- file.path(phase07_sections, paste0("phase07_section_", sec, "_retention_scaled.rds"))
  f_tsv <- file.path(phase07_sections, paste0("phase07_section_", sec, "_retention_scaled.tsv"))
  if (file.exists(f_rds)) return(readRDS(f_rds))
  if (file.exists(f_tsv)) return(read_tsv(f_tsv, show_col_types = FALSE))
  stop("Missing retention_scaled file for section ", sec)
}

load_sec_zone_df <- function(sec) {
  path <- file.path(phase06_tables, paste0("zones_section_", sec, ".tsv"))
  if (!file.exists(path)) stop("Missing zone table for section ", sec)
  df <- read_tsv(path, show_col_types = FALSE)
  if ("cell.x" %in% colnames(df) && !"cell" %in% colnames(df))
    df <- df %>% rename(cell = cell.x)
  df %>% mutate(
    section_id       = as.character(section_id),
    in_tissue        = as.integer(in_tissue),
    is_hotspot_semla = ifelse(is.na(is_hotspot_semla), FALSE, as.logical(is_hotspot_semla)),
    zone_call        = factor(as.character(zone_call), levels = zone_levels)
  )
}

# =========================================================
# PANEL A
# =========================================================
message("Panel a")

# States actually represented across the sections displayed in panel A
panelA_states_present <- panelA_sections %>%
  lapply(function(sec) {
    read_plotready(sec) %>%
      filter(
        in_tissue == 1,
        is_hotspot_semla,
        !is.na(dominant_state_pretty)
      ) %>%
      pull(dominant_state_pretty) %>%
      as.character() %>%
      unique()
  }) %>%
  unlist() %>%
  unique()

panelA_state_levels <- state_levels[state_levels %in% panelA_states_present]


plot_panel_A <- function(sec, show_legend = FALSE) {
  
  df <- read_plotready(sec) %>%
    mutate(
      section_id = as.character(section_id),
      in_tissue = as.integer(in_tissue),
      is_hotspot_semla = ifelse(
        is.na(is_hotspot_semla),
        FALSE,
        as.logical(is_hotspot_semla)
      ),
      dominant_state_pretty = factor(
        dominant_state_pretty,
        levels = panelA_state_levels
      )
    )
  
  prep <- prepare_he_plot(df)
  d <- prep$df
  
  d0 <- d %>%
    filter(in_tissue == 0, is.finite(x_plot), is.finite(y_plot))
  
  d1 <- d %>%
    filter(in_tissue == 1, is.finite(x_plot), is.finite(y_plot))
  
  dh <- d1 %>%
    filter(is_hotspot_semla, !is.na(dominant_state_pretty))
  
  p <- make_he_canvas(prep, sec, he_alpha = 1)
  
  if (nrow(d0) > 0) {
    p <- p +
      geom_point(
        data = d0,
        aes(x_plot, y_plot),
        shape = 1,
        stroke = 0.35,
        color = "grey70",
        size = PT_SIZE_OFF,
        inherit.aes = FALSE
      )
  }
  
  if (nrow(d1) > 0) {
    p <- p +
      geom_point(
        data = d1,
        aes(x_plot, y_plot),
        color = alpha("grey85", 0.55),
        size = PT_SIZE_TISSUE,
        inherit.aes = FALSE,
        show.legend = FALSE
      )
  }
  
  if (nrow(dh) > 0) {
    p <- p +
      geom_point(
        data = dh,
        aes(x_plot, y_plot, colour = dominant_state_pretty),
        size = PT_SIZE_STATE,
        inherit.aes = FALSE,
        show.legend = show_legend
      ) +
      scale_colour_manual(
        values = state_cols,
        limits = panelA_state_levels,
        breaks = panelA_state_levels,
        drop = FALSE,
        name = NULL,
        guide = guide_legend(
          override.aes = list(size = 3)
        )
      )
  }
  
  p
}

# Only the second plot generates the legend
plots_A <- list(
  plot_panel_A(panelA_sections[1], show_legend = FALSE),
  plot_panel_A(panelA_sections[2], show_legend = TRUE)
)

fig4A <- wrap_plots(
  plots_A,
  nrow = 1
) &
  theme(legend.position = "right")

ggsave(file.path(out_dir, "Figure4a_hotspot_states_459_928.png"),
       fig4A, width = 8.4, height = 4.4, dpi = PLOT$dpi, bg = "white")
ggsave(file.path(out_dir, "Figure4a_hotspot_states_459_928.pdf"),
       fig4A, width = 8.4, height = 4.4, bg = "white")

# =========================================================
# PANEL B
# =========================================================
message("Panel b")

freq_path <- file.path(phase07_tables, "phase07_state_frequency_by_section.tsv")
if (!file.exists(freq_path)) stop("Missing frequency table: ", freq_path)

freq_df <- read_tsv(freq_path, show_col_types = FALSE) %>%
  mutate(
    section_id = factor(as.character(section_id), levels = good_sections),
    dominant_state_pretty = as.character(dominant_state_pretty)
  )

# Keep the complete table, including Unassigned, for QC/reporting
write_tsv(
  freq_df,
  file.path(out_dir, "Figure4b_hotspot_state_frequency_all.tsv")
)

# Figure 4b: composition of confidently classified hotspots only
freq_plot <- freq_df %>%
  filter(dominant_state_pretty != "Unassigned") %>%
  group_by(section_id) %>%
  mutate(
    n_classified = sum(n_hotspots_state),
    pct_classified = ifelse(
      n_classified > 0,
      100 * n_hotspots_state / n_classified,
      0
    )
  ) %>%
  ungroup() %>%
  mutate(
    dominant_state_pretty = factor(
      dominant_state_pretty,
      levels = state_levels[state_levels != "Unassigned"]
    )
  ) %>%
  arrange(section_id, dominant_state_pretty)

write_tsv(
  freq_plot,
  file.path(out_dir, "Figure4b_hotspot_state_frequency_classified.tsv")
)

p4B <- ggplot(
  freq_plot,
  aes(section_id, pct_classified, fill = dominant_state_pretty)
) +
  geom_col(width = 0.82, color = "black", linewidth = 0.2) +
  scale_fill_manual(
    values = state_cols,
    drop = FALSE,
    name = NULL
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100),
    expand = c(0, 0)
  ) +
  labs(
    x = NULL,
    y = "% classified hotspots"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title.y = element_text(face = "bold"),
    legend.position = "right",
    plot.margin = margin(8, 8, 8, 8)
  )

ggsave(file.path(out_dir, "Figure4b_hotspot_state_frequency.png"),
       p4B, width = 7.6, height = 4.3, dpi = PLOT$dpi, bg = "white")
ggsave(file.path(out_dir, "Figure4b_hotspot_state_frequency.pdf"),
       p4B, width = 7.6, height = 4.3, bg = "white")

# =========================================================
# PANEL C — Retention gene grid (section 459)
# =========================================================
message("Panel c")

panelC_df <- read_retention(panelC_section) %>%
  mutate(section_id       = as.character(section_id),
         in_tissue        = as.integer(in_tissue),
         is_hotspot_semla = ifelse(is.na(is_hotspot_semla), FALSE,
                                   as.logical(is_hotspot_semla)))

# Retention genes only (GPNMB handled separately in panel d)
retention_genes_c <- setdiff(retention_genes, "GPNMB")

plot_retention_gene <- function(df, gene, show_legend = FALSE) {
  feat <- paste0("scaled_", gene)
  if (!feat %in% colnames(df)) { warning("Missing: ", feat); return(NULL) }
  prep <- prepare_he_plot(df); d <- prep$df
  d1   <- d %>% filter(in_tissue == 1, is.finite(x_plot), is.finite(.data[[feat]]))
  dh   <- d1 %>% filter(is_hotspot_semla)
  p    <- make_he_canvas(prep, NULL, he_alpha = HE_ALPHA_PANELC)
  if (nrow(d1) > 0)
    p <- p +
    geom_point(data = d1, aes(x_plot, y_plot, colour = .data[[feat]]),
               size = PT_SIZE_EXPR - 0.40, alpha = 0.70,
               inherit.aes = FALSE, show.legend = show_legend) +
    scale_colour_gradientn(colours = viridisLite::plasma(256),
                           limits = c(0, 1), oob = scales::squish,
                           name = "Scaled\nexpression")
  add_hotspot_rings(p, dh) +
    labs(title = gene) +
    theme(plot.title = element_text(size = 11, face = "italic", hjust = 0.5))
}

plots_c <- lapply(seq_along(retention_genes_c), function(i)
  plot_retention_gene(panelC_df, retention_genes_c[i],
                      show_legend = (i == length(retention_genes_c))))
plots_c <- Filter(Negate(is.null), plots_c)

fig4C <- wrap_plots(plots_c, ncol = 6, guides = "collect") &
  theme(legend.position = "right")

ggsave(file.path(out_dir, "Figure4c_retention_genes_section459.png"),
       fig4C, width = 15.0, height = 6.2, dpi = PLOT$dpi, bg = "white")
ggsave(file.path(out_dir, "Figure4c_retention_genes_section459.pdf"),
       fig4C, width = 15.0, height = 6.2, bg = "white")

# =========================================================
# PANELS D + E — GPNMB
# =========================================================
message("Panels d/e — GPNMB")

# ---------------------------------------------------------
# Settings
# ---------------------------------------------------------

# Sections retaining the recurrent hotspot zone-confinement pattern.
# If already defined near the top of figure4.R, this line can be removed.
# niche_sections defined centrally above

# Only zones with sufficient hotspot representation and relevant
# to the recurrent hotspot niche are tested/displayed in panel E.
gpnmb_display_zones <- c("Mesenchymal", "Myeloid")

# Minimum number of spots required in BOTH groups
# within a section × zone comparison.
MIN_SPOTS_PER_GROUP <- 10L


# ---------------------------------------------------------
# Check GPNMB is present in phase07 retention files
# ---------------------------------------------------------

test_df   <- read_retention(good_sections[1])
has_gpnmb <- "scaled_GPNMB" %in% colnames(test_df)

if (!has_gpnmb) {
  
  message(
    "WARNING: scaled_GPNMB not found — ",
    "add GPNMB to retention_genes in phase07 and rerun."
  )
  
} else {
  
  
  # =======================================================
  # PANEL D — Spatial GPNMB maps
  # =======================================================
  
  plot_gpnmb_map <- function(sec, show_legend = TRUE) {
    
    df <- read_retention(sec) %>%
      mutate(
        section_id = as.character(section_id),
        in_tissue  = as.integer(in_tissue),
        is_hotspot_semla = ifelse(
          is.na(is_hotspot_semla),
          FALSE,
          as.logical(is_hotspot_semla)
        )
      )
    
    prep <- prepare_he_plot(df)
    d    <- prep$df
    
    d1 <- d %>%
      filter(
        in_tissue == 1,
        is.finite(x_plot),
        is.finite(y_plot),
        is.finite(scaled_GPNMB)
      )
    
    dh <- d1 %>%
      filter(is_hotspot_semla)
    
    p <- make_he_canvas(
      prep,
      NULL,
      he_alpha = HE_ALPHA_PANELC
    )
    
    if (nrow(d1) > 0) {
      
      p <- p +
        geom_point(
          data = d1,
          aes(
            x = x_plot,
            y = y_plot,
            colour = scaled_GPNMB
          ),
          size = PT_SIZE_EXPR - 0.40,
          alpha = 0.70,
          inherit.aes = FALSE,
          show.legend = show_legend
        ) +
        scale_colour_gradientn(
          colours = viridisLite::plasma(256),
          limits = c(0, 1),
          oob = scales::squish,
          name = "Scaled\nexpression"
        )
    }
    
    p <- add_hotspot_rings(p, dh)
    
    p +
      labs(title = "GPNMB") +
      theme(
        plot.title = element_text(
          size = 11,
          face = "italic",
          hjust = 0.5
        )
      )
  }
  
  
  fig4D <- wrap_plots(
    list(
      plot_gpnmb_map(
        "459",
        show_legend = FALSE
      ),
      plot_gpnmb_map(
        "928",
        show_legend = TRUE
      )
    ),
    nrow = 1,
    guides = "collect"
  ) &
    theme(
      legend.position = "right"
    )
  
  
  ggsave(
    file.path(
      out_dir,
      "Figure4d_GPNMB_maps_459_928.png"
    ),
    fig4D,
    width = 8.8,
    height = 4.4,
    dpi = PLOT$dpi,
    bg = "white"
  )
  
  ggsave(
    file.path(
      out_dir,
      "Figure4d_GPNMB_maps_459_928.pdf"
    ),
    fig4D,
    width = 8.8,
    height = 4.4,
    bg = "white"
  )
  
  
  
  # =======================================================
  # PANEL E — Section-wise GPNMB comparison
  # =======================================================
  #
  # Statistical unit = spatial section.
  #
  # Analysis restricted to:
  #   - sections displaying the recurrent hotspot
  #     zone-confinement pattern
  #   - Mesenchymal and Myeloid zones
  #   - section × zone combinations with >=10 spots in
  #     BOTH hotspot and non-hotspot groups
  #
  # =======================================================
  
  message(
    "Panel e — section-wise GPNMB comparison: ",
    paste(gpnmb_display_zones, collapse = ", ")
  )
  
  
  # -------------------------------------------------------
  # E1. Assemble spot-level GPNMB data
  # -------------------------------------------------------
  
  gpnmb_all <- purrr::map_dfr(
    niche_sections,
    function(sec) {
      
      df <- read_retention(sec) %>%
        mutate(
          section_id = as.character(section_id),
          in_tissue  = as.integer(in_tissue),
          is_hotspot_semla = ifelse(
            is.na(is_hotspot_semla),
            FALSE,
            as.logical(is_hotspot_semla)
          )
        )
      
      zone_df <- load_sec_zone_df(sec) %>%
        select(
          cell,
          zone_call
        )
      
      join_col <- intersect(
        c(
          "cell.x",
          "cell",
          "barcode",
          "spot_id"
        ),
        colnames(df)
      )[1]
      
      if (is.na(join_col)) {
        df$cell <- rownames(df)
        join_col <- "cell"
      }
      
      zone_df <- zone_df %>%
        rename(
          !!join_col := cell
        )
      
      df %>%
        filter(
          in_tissue == 1,
          is.finite(scaled_GPNMB)
        ) %>%
        left_join(
          zone_df,
          by = join_col
        ) %>%
        select(
          section_id,
          is_hotspot_semla,
          scaled_GPNMB,
          zone_call
        )
    }
  ) %>%
    mutate(
      group = ifelse(
        is_hotspot_semla,
        "Hotspot",
        "Non-hotspot"
      ),
      group = factor(
        group,
        levels = c(
          "Non-hotspot",
          "Hotspot"
        )
      ),
      zone_call = as.character(zone_call)
    ) %>%
    filter(
      !is.na(zone_call),
      zone_call != "Uncertain"
    )
  
  
  # Save complete niche-section spot-level table for QC
  write_tsv(
    gpnmb_all,
    file.path(
      out_dir,
      "Figure4e_GPNMB_spotlevel_input.tsv"
    )
  )
  
  
  # -------------------------------------------------------
  # E2. Restrict panel E to Mesenchymal + Myeloid
  # -------------------------------------------------------
  
  gpnmb_focus <- gpnmb_all %>%
    filter(
      zone_call %in% gpnmb_display_zones
    ) %>%
    mutate(
      zone_call = factor(
        zone_call,
        levels = gpnmb_display_zones
      )
    )
  
  
  # -------------------------------------------------------
  # E3. Section-level summaries
  # -------------------------------------------------------
  
  gpnmb_section <- gpnmb_focus %>%
    group_by(
      section_id,
      zone_call,
      group
    ) %>%
    summarise(
      mean_GPNMB = mean(
        scaled_GPNMB,
        na.rm = TRUE
      ),
      
      median_GPNMB = median(
        scaled_GPNMB,
        na.rm = TRUE
      ),
      
      n_spots = n(),
      
      .groups = "drop"
    )
  
  
  write_tsv(
    gpnmb_section,
    file.path(
      out_dir,
      "Figure4e_GPNMB_section_summary.tsv"
    )
  )
  
  
  # -------------------------------------------------------
  # E4. Minimum representation filter
  #
  # Require >=10 hotspot AND >=10 non-hotspot spots
  # in the same section × zone.
  # -------------------------------------------------------
  
  gpnmb_section_use <- gpnmb_section %>%
    group_by(
      section_id,
      zone_call
    ) %>%
    filter(
      any(
        group == "Hotspot" &
          n_spots >= MIN_SPOTS_PER_GROUP
      ),
      any(
        group == "Non-hotspot" &
          n_spots >= MIN_SPOTS_PER_GROUP
      )
    ) %>%
    ungroup()
  
  
  # Save representation QC
  write_tsv(
    gpnmb_section_use,
    file.path(
      out_dir,
      "Figure4e_GPNMB_section_summary_filtered.tsv"
    )
  )
  
  
  # -------------------------------------------------------
  # E5. Build paired section-level table
  # -------------------------------------------------------
  
  gpnmb_paired <- gpnmb_section_use %>%
    select(
      section_id,
      zone_call,
      group,
      mean_GPNMB
    ) %>%
    pivot_wider(
      names_from = group,
      values_from = mean_GPNMB
    ) %>%
    filter(
      is.finite(`Non-hotspot`),
      is.finite(`Hotspot`)
    ) %>%
    mutate(
      delta = Hotspot - `Non-hotspot`
    )
  
  
  write_tsv(
    gpnmb_paired,
    file.path(
      out_dir,
      "Figure4e_GPNMB_section_pairs.tsv"
    )
  )
  
  
  # -------------------------------------------------------
  # E6. Paired section-wise statistics
  #
  # Paired Wilcoxon signed-rank test.
  # Minimum 3 paired sections required.
  # -------------------------------------------------------
  
  wt_zone <- gpnmb_paired %>%
    group_by(
      zone_call
    ) %>%
    summarise(
      n_sections = n(),
      
      n_increased = sum(
        delta > 0,
        na.rm = TRUE
      ),
      
      n_decreased = sum(
        delta < 0,
        na.rm = TRUE
      ),
      
      mean_nonhotspot = mean(
        `Non-hotspot`,
        na.rm = TRUE
      ),
      
      mean_hotspot = mean(
        Hotspot,
        na.rm = TRUE
      ),
      
      mean_delta = mean(
        delta,
        na.rm = TRUE
      ),
      
      median_delta = median(
        delta,
        na.rm = TRUE
      ),
      
      p_value = if (n() >= 3) {
        
        tryCatch(
          wilcox.test(
            Hotspot,
            `Non-hotspot`,
            paired = TRUE,
            exact = FALSE
          )$p.value,
          error = function(e) NA_real_
        )
        
      } else {
        
        NA_real_
      },
      
      .groups = "drop"
    ) %>%
    
    # BH correction across the two interpretable zones
    mutate(
      p_adj = p.adjust(
        p_value,
        method = "BH"
      ),
      
      label = case_when(
        
        n_sections < 3 ~ paste0(
          "n = ",
          n_sections,
          "\ninsufficient representation"
        ),
        
        TRUE ~ paste0(
          n_increased,
          "/",
          n_sections,
          " \u2191",
          "\nadj. p = ",
          format.pval(
            p_adj,
            digits = 2,
            eps = 1e-300
          )
        )
      )
    )
  
  
  write_tsv(
    wt_zone,
    file.path(
      out_dir,
      "Figure4e_GPNMB_by_zone_stats.tsv"
    )
  )
  
  
  # -------------------------------------------------------
  # E7. Plot
  #
  # Each point = section-level mean
  # Each line  = paired values from the same section
  # -------------------------------------------------------
  
  plot_df_E <- gpnmb_section_use %>%
    semi_join(
      gpnmb_paired %>%
        select(
          section_id,
          zone_call
        ),
      by = c(
        "section_id",
        "zone_call"
      )
    ) %>%
    mutate(
      group = factor(
        group,
        levels = c(
          "Non-hotspot",
          "Hotspot"
        )
      ),
      
      zone_call = factor(
        zone_call,
        levels = gpnmb_display_zones
      )
    )
  
  
  p4E <- ggplot(
    plot_df_E,
    aes(
      x = group,
      y = mean_GPNMB,
      group = section_id
    )
  ) +
    
    # Paired section lines
    geom_line(
      colour = "grey65",
      linewidth = 0.65,
      alpha = 0.85
    ) +
    
    # Section-level means
    geom_point(
      aes(
        fill = zone_call
      ),
      shape = 21,
      colour = "black",
      stroke = 0.40,
      size = 3.4
    ) +
    
    facet_wrap(
      ~ zone_call,
      nrow = 1,
      drop = TRUE
    ) +
    
    # Directional consistency + BH-adjusted p-value
    geom_text(
      data = wt_zone %>%
        mutate(
          zone_call = factor(
            zone_call,
            levels = gpnmb_display_zones
          )
        ),
      aes(
        x = 1.5,
        y = Inf,
        label = label
      ),
      inherit.aes = FALSE,
      vjust = 1.10,
      size = 3.7,
      lineheight = 0.95
    ) +
    
    scale_fill_manual(
      values = zone_cols,
      drop = FALSE
    ) +
    
    scale_y_continuous(
      expand = expansion(
        mult = c(
          0.05,
          0.18
        )
      )
    ) +
    
    labs(
      x = NULL,
      y = "Mean scaled GPNMB expression"
    ) +
    
    theme_classic(
      base_size = 12
    ) +
    
    theme(
      axis.text.x = element_text(
        angle = 30,
        hjust = 1
      ),
      
      axis.title.y = element_text(
        face = "bold"
      ),
      
      strip.text = element_text(
        face = "bold",
        size = 12
      ),
      
      strip.background = element_blank(),
      
      legend.position = "none",
      
      panel.spacing.x = grid::unit(
        12,
        "pt"
      )
    )
  
  
  write_tsv(
    plot_df_E,
    file.path(
      out_dir,
      "Figure4e_GPNMB_section_plotdata.tsv"
    )
  )
  
  
  ggsave(
    file.path(
      out_dir,
      "Figure4e_GPNMB_by_zone.png"
    ),
    p4E,
    width = 7.2,
    height = 5,
    dpi = PLOT$dpi,
    bg = "white"
  )
  
  ggsave(
    file.path(
      out_dir,
      "Figure4e_GPNMB_by_zone.pdf"
    ),
    p4E,
    width = 7.2,
    height = 5,
    bg = "white"
  )
  
  
  # -------------------------------------------------------
  # Console summary
  # -------------------------------------------------------
  
  message(
    "\nPanel E — paired section-level GPNMB statistics ",
    "(niche sections; minimum ",
    MIN_SPOTS_PER_GROUP,
    " spots/group):"
  )
  
  print(
    wt_zone %>%
      select(
        zone_call,
        n_sections,
        n_increased,
        n_decreased,
        mean_nonhotspot,
        mean_hotspot,
        mean_delta,
        p_value,
        p_adj
      )
  )
}

# =========================================================
# PANEL F — LR consensus dot plot (v7)
# Reads from 08_lr_interactions/final_pair_sets/LR_pairs_joint_*.tsv
# No in-memory objects required.
# =========================================================
message("Panel f — LR dot plot")

out_lr_root <- stage_dir("lr_interactions")
out_final   <- file.path(out_lr_root, "final_pair_sets")
out_nichen  <- file.path(out_lr_root, "nichenet")
out_res     <- file.path(out_lr_root, "resources")
raw_vis_dir <- stage_dir("raw_vis")

all_sections_kept <- niche_sections
zones_to_plot     <- c("Myeloid", "Mesenchymal", "Vascular")
zone_exclude_sections <- list(
  Myeloid     = c("1239"),
  Mesenchymal = character(0),
  Vascular    = c("928")
)
n_top_pairs_per_zone <- 10

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

stop_if_missing <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path)
  path
}

zone_n_tested_tbl <- function(zones, section_keep, zone_exclude) {
  tibble(
    zone     = zones,
    n_tested = vapply(zones, function(z)
      length(setdiff(section_keep, zone_exclude[[z]] %||% character(0))),
      integer(1))
  )
}

read_joint_zone <- function(zone) {
  f <- stop_if_missing(file.path(out_final, paste0("LR_pairs_joint_", zone, ".tsv")))
  read_tsv(f, show_col_types = FALSE) %>%
    mutate(
      zone          = zone,
      ligand        = toupper(ligand),
      receptor      = toupper(receptor),
      pair_label    = paste0(ligand, "–", receptor),
      min_n_sections = pmin(n_sections_cellchat, n_sections_cellphonedb, na.rm = TRUE)
    )
}

tryCatch({
  zone_n_tested <- zone_n_tested_tbl(zones_to_plot, all_sections_kept,
                                     zone_exclude_sections)
  joint_list  <- lapply(zones_to_plot, read_joint_zone)
  pair_overlap <- bind_rows(joint_list) %>% filter(in_both %in% TRUE)
  
  panelD_keep <- pair_overlap %>%
    group_by(zone) %>%
    arrange(desc(min_score_both), desc(mean_score_both), desc(min_n_sections)) %>%
    slice_head(n = n_top_pairs_per_zone) %>%
    ungroup() %>%
    left_join(zone_n_tested, by = "zone") %>%
    mutate(support_pct   = 100 * min_n_sections / n_tested,
           support_label = paste0(min_n_sections, "/", n_tested))
  
  pair_order <- panelD_keep %>%
    group_by(pair_label) %>%
    summarise(order_score = max(min_score_both, na.rm = TRUE), .groups = "drop") %>%
    arrange(order_score) %>% pull(pair_label)
  
  plot_df_F <- panelD_keep %>%
    mutate(zone       = factor(zone,       levels = zones_to_plot),
           pair_label = factor(pair_label, levels = pair_order)) %>%
    select(zone, ligand, receptor, pair_label, min_score_both, mean_score_both,
           min_n_sections, n_tested, support_pct, support_label)
  
  p4F <- ggplot(plot_df_F, aes(x = zone, y = pair_label)) +
    geom_point(aes(size = support_pct, fill = min_score_both),
               shape = 21, colour = "black", stroke = 0.25) +
    scale_fill_viridis_c(option = "plasma", direction = 1,
                         name = "Shared-pair\nstrength\n(min score)") +
    scale_size_continuous(range = c(2.5, 10),
                          breaks = c(25, 50, 75, 100),
                          labels = function(x) paste0(x, "%"),
                          limits = c(0, 100),
                          name = "Supporting sections\n(% of included sections)") +
    labs(x = NULL, y = "Ligand\u2013receptor pair") +
    theme_bw(base_size = 12) +
    theme(axis.text.y = element_text(size = 9),
          axis.text.x = element_text(face = "bold"),
          panel.grid.major = element_line(colour = "grey92", linewidth = 0.25),
          panel.grid.minor = element_blank(),
          legend.box = "vertical")
  
  write_tsv(plot_df_F, file.path(out_dir, "Figure4f_LR_dotplot.tsv"))
  ggsave(file.path(out_dir, "Figure4f_LR_dotplot.png"),
         p4F, width = 8.8, height = 7.8, dpi = PLOT$dpi, bg = "white")
  ggsave(file.path(out_dir, "Figure4f_LR_dotplot.pdf"),
         p4F, width = 8.8, height = 7.8, device = cairo_pdf, bg = "white")
  message("Panel f done.")
  
}, error = function(e) {
  message("Panel f skipped — LR joint pair TSVs not found.\n  ", conditionMessage(e))
})

# =========================================================
# PANELS G + H — NicheNet heatmap + S100A9/CCL2 spatial maps
# Reads from 08_lr_interactions/nichenet/ and 00_raw_sections/vis/
# =========================================================
message("Panels g/h — NicheNet")

panelE_zone          <- "Myeloid"
panelE_n_ligands     <- 6
panelE_n_targets     <- 18
manual_panelE_ligands <- NULL
manual_best_section  <- "928"
manual_rep_ligand    <- "APP"
manual_rep_targets   <- c("S100A9", "CCL2")

PT_SIZE_EXPR_GH      <- 0.55
PT_SIZE_OFF_GH       <- 0.45
HOT_RING_SIZE_GH     <- 1.00
HOT_RING_STR_DARK    <- 0.80
HOT_RING_STR_LIGHT   <- 0.60

read_expr_matrix <- function(obj) {
  DefaultAssay(obj) <- "Spatial"
  mat <- tryCatch(GetAssayData(obj, assay = "Spatial", layer = "data"),
                  error = function(e) NULL)
  if (is.null(mat) || nrow(mat) == 0 || ncol(mat) == 0) {
    message("  normalising from counts")
    obj <- NormalizeData(obj, normalization.method = "LogNormalize",
                         scale.factor = 1e4, verbose = FALSE)
    mat <- tryCatch(GetAssayData(obj, assay = "Spatial", layer = "data"),
                    error = function(e) NULL)
  }
  if (is.null(mat) || nrow(mat) == 0) stop("Could not recover Spatial data layer.")
  list(obj = obj, mat = mat)
}

find_cell_col <- function(df) {
  cand <- intersect(c("cell","barcode","barcode_raw","spot","spot_id","Cell","Barcode"),
                    colnames(df))
  if (length(cand) == 0) NA_character_ else cand[1]
}

prepare_phase06_table <- function(df, obj, sec) {
  df <- df
  if ("cell.x" %in% colnames(df) && !"cell" %in% colnames(df))
    df <- df %>% rename(cell = cell.x)
  df$img_file <- section_img_path(sec)
  cell_col <- find_cell_col(df)
  obj_cells <- colnames(obj)
  if (!is.na(cell_col)) {
    df2 <- df %>% mutate(cell = as.character(.data[[cell_col]]))
    if (sum(df2$cell %in% obj_cells) > 0)
      return(df2 %>% filter(cell %in% obj_cells) %>% distinct(cell, .keep_all = TRUE))
  }
  if (nrow(df) != length(obj_cells))
    stop("Phase06 table row count mismatch for section ", sec)
  df %>% mutate(cell = obj_cells)
}

feature_expr_vec <- function(mat, feature, cells) {
  feature <- toupper(feature)
  rn_up   <- toupper(rownames(mat))
  if (!(feature %in% rn_up)) return(rep(NA_real_, length(cells)))
  as.numeric(mat[match(feature, rn_up), cells, drop = TRUE])
}

prepare_he_plot_gh <- function(df_sec) {
  sec <- unique(as.character(df_sec$section_id))
  if (length(sec) != 1 || is.na(sec)) stop("Expected one section in prepare_he_plot_gh().")
  img_path <- section_img_path(sec)
  if (is.na(img_path) || !file.exists(img_path))
    stop("Missing H&E image for section ", sec, " under RAW_ROOT: ", RAW_ROOT)
  img   <- png::readPNG(img_path)
  img_h <- nrow(img); img_w <- ncol(img)
  df2   <- df_sec %>% mutate(x_plot = pxl_col_in_hires,
                             y_plot = img_h - pxl_row_in_hires)
  list(df = df2, img = img, img_h = img_h, img_w = img_w)
}

plot_feature_on_he <- function(sec, feature, title_txt = NULL) {
  vis_sec  <- readRDS(file.path(raw_vis_dir, paste0("vis_section_", sec, "_raw.rds")))
  expr_res <- read_expr_matrix(vis_sec)
  expr_mat <- expr_res$mat
  ph  <- read_tsv(file.path(phase06_tables, paste0("zones_section_", sec, ".tsv")),
                  show_col_types = FALSE)
  ph  <- prepare_phase06_table(ph, expr_res$obj, sec) %>%
    mutate(section_id       = as.character(section_id),
           in_tissue        = as.integer(in_tissue),
           is_hotspot_semla = ifelse(is.na(is_hotspot_semla), FALSE,
                                     as.logical(is_hotspot_semla)),
           cell             = as.character(cell))
  ph$expr <- feature_expr_vec(expr_mat, feature, ph$cell)
  prep <- prepare_he_plot_gh(ph); d <- prep$df
  d0 <- d %>% filter(in_tissue == 0, is.finite(x_plot), is.finite(y_plot))
  d1 <- d %>% filter(in_tissue == 1, is.finite(x_plot), is.finite(y_plot),
                     is.finite(expr))
  dh <- d1 %>% filter(is_hotspot_semla)
  if (nrow(d1) > 0) {
    q99 <- quantile(d1$expr, 0.99, na.rm = TRUE)
    if (!is.finite(q99) || q99 <= 0) q99 <- max(d1$expr, na.rm = TRUE)
    if (!is.finite(q99) || q99 <= 0) q99 <- 1
    d1 <- d1 %>% mutate(expr_plot = pmin(expr, q99) / q99)
  }
  pal <- scales::alpha(viridisLite::plasma(256),
                       seq(0.18, 1, length.out = 256))
  g   <- grid::rasterGrob(prep$img, interpolate = TRUE)
  pad <- ceiling(max(prep$img_h, prep$img_w) * HE_PAD_FRAC)
  p   <- ggplot() +
    annotation_custom(g, xmin = 0, xmax = prep$img_w,
                      ymin = 0, ymax = prep$img_h) +
    coord_fixed(xlim = c(-pad, prep$img_w + pad),
                ylim = c(-pad, prep$img_h + pad),
                expand = FALSE, clip = "off") +
    theme_void(base_size = 11) +
    ggtitle(title_txt) +
    theme(plot.title   = element_text(face = "italic", hjust = 0.5, size = 11),
          legend.title = element_text(size = 9.5, face = "bold"),
          legend.text  = element_text(size = 8.5))
  if (nrow(d0) > 0)
    p <- p + geom_point(data = d0, aes(x_plot, y_plot),
                        shape = 1, stroke = 0.30,
                        color = scales::alpha("grey65", 0.45),
                        size = PT_SIZE_OFF_GH, inherit.aes = FALSE)
  if (nrow(d1) > 0)
    p <- p +
    geom_point(data = d1, aes(x_plot, y_plot, colour = expr_plot),
               size = PT_SIZE_EXPR_GH, alpha = 0.95, inherit.aes = FALSE) +
    scale_colour_gradientn(colours = pal, limits = c(0, 1),
                           oob = scales::squish, name = "Scaled\nexpression")
  if (nrow(dh) > 0)
    p <- p +
    geom_point(data = dh, aes(x_plot, y_plot),
               shape = 21, fill = NA, colour = "black",
               stroke = HOT_RING_STR_DARK, size = HOT_RING_SIZE_GH,
               inherit.aes = FALSE) +
    geom_point(data = dh, aes(x_plot, y_plot),
               shape = 21, fill = NA, colour = "white",
               stroke = HOT_RING_STR_LIGHT, size = HOT_RING_SIZE_GH - 0.10,
               inherit.aes = FALSE)
  p
}

tryCatch({
  act_zone <- read_tsv(stop_if_missing(
    file.path(out_nichen, paste0("ligand_activity_", panelE_zone, ".tsv"))),
    show_col_types = FALSE) %>% mutate(ligand = toupper(ligand))
  act_sec  <- read_tsv(stop_if_missing(
    file.path(out_nichen, paste0("ligand_activity_", panelE_zone, "_by_section.tsv"))),
    show_col_types = FALSE) %>% mutate(ligand = toupper(ligand),
                                       section_id = as.character(section_id))
  cons_df  <- read_tsv(stop_if_missing(
    file.path(out_nichen, paste0("consensus_upgenes_", panelE_zone, ".tsv"))),
    show_col_types = FALSE) %>% mutate(gene = toupper(gene))
  lt_mat <- readRDS(stop_if_missing(
    file.path(out_res, "nichenet", "ligand_target_matrix_nsga2r_final.rds")))
  
  lt_ligands <- toupper(colnames(lt_mat))
  lt_targets <- toupper(rownames(lt_mat))
  
  joint_mye  <- read_joint_zone(panelE_zone) %>%
    filter(in_both %in% TRUE) %>%
    arrange(desc(min_score_both), desc(mean_score_both), desc(min_n_sections))
  candidate_ligands <- unique(joint_mye$ligand)
  
  ligands_selected <- if (!is.null(manual_panelE_ligands)) {
    toupper(manual_panelE_ligands)
  } else {
    act_zone %>% filter(ligand %in% candidate_ligands) %>%
      arrange(desc(auc), desc(mean_pos - mean_neg)) %>%
      pull(ligand) %>% unique() %>% head(panelE_n_ligands)
  }
  ligands_selected <- ligands_selected[ligands_selected %in% lt_ligands]
  
  cons_use <- cons_df %>%
    mutate(gene_weight = mean_log2FC * pmax(n_sections, 1L),
           lt_gene = rownames(lt_mat)[match(gene, lt_targets)]) %>%
    filter(!is.na(lt_gene))
  
  weighted_long <- bind_rows(lapply(ligands_selected, function(L) {
    col_idx <- match(L, lt_ligands)
    vals    <- as.numeric(lt_mat[cons_use$lt_gene, col_idx])
    tibble(ligand = L, gene = cons_use$gene,
           n_sections = cons_use$n_sections,
           mean_log2FC = cons_use$mean_log2FC,
           regulatory_potential = vals,
           weighted_contribution = vals * cons_use$gene_weight)
  }))
  
  selected_targets <- weighted_long %>%
    group_by(gene) %>%
    summarise(max_weighted = max(weighted_contribution, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(max_weighted)) %>% slice_head(n = panelE_n_targets) %>% pull(gene)
  if (!is.null(manual_rep_targets))
    selected_targets <- unique(c(selected_targets, toupper(manual_rep_targets)))
  
  heat_df <- weighted_long %>% filter(gene %in% selected_targets)
  max_wc  <- max(heat_df$weighted_contribution, na.rm = TRUE)
  if (!is.finite(max_wc) || max_wc <= 0) max_wc <- 1
  heat_df <- heat_df %>% mutate(weighted_scaled = weighted_contribution / max_wc)
  
  ligand_labels <- joint_mye %>%
    filter(ligand %in% ligands_selected) %>%
    group_by(ligand) %>%
    arrange(desc(min_score_both)) %>% slice_head(n = 1) %>% ungroup() %>%
    transmute(ligand, ligand_label = paste0(ligand, " [", ligand, "\u2013", receptor, "]"))
  
  heat_df <- heat_df %>%
    left_join(ligand_labels, by = "ligand") %>%
    mutate(ligand_label = ifelse(is.na(ligand_label), ligand, ligand_label))
  
  target_order <- heat_df %>%
    group_by(gene) %>%
    summarise(mx = max(weighted_contribution, na.rm = TRUE), .groups = "drop") %>%
    arrange(mx) %>% pull(gene)
  
  heat_df$gene         <- factor(heat_df$gene, levels = target_order)
  heat_df$ligand_label <- factor(heat_df$ligand_label)
  
  # Panel g — NicheNet heatmap
  p4G <- ggplot(heat_df, aes(x = gene, y = ligand_label, fill = weighted_scaled)) +
    geom_tile(colour = "white", linewidth = 0.25) +
    scale_fill_viridis_c(option = "plasma", direction = 1,
                         limits = c(0, 1), oob = scales::squish,
                         name = "Scaled weighted\ncontribution") +
    labs(x = "Predicted downstream target genes",
         y = "Selected ligands [best overlapping LR pair]") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1,
                                     face = "italic", size = 8.5),
          axis.text.y = element_text(size = 9.2),
          axis.title  = element_text(face = "bold"),
          panel.grid  = element_blank())
  
  write_tsv(heat_df, file.path(out_dir, "Figure4g_NicheNet_heatmap.tsv"))
  ggsave(file.path(out_dir, "Figure4g_NicheNet_heatmap.png"),
         p4G, width = 10.0, height = 5.5, dpi = PLOT$dpi, bg = "white")
  ggsave(file.path(out_dir, "Figure4g_NicheNet_heatmap.pdf"),
         p4G, width = 10.0, height = 5.5, device = cairo_pdf, bg = "white")
  message("Panel g done.")
  
  # Panel h — S100A9 & CCL2 spatial maps (side by side)
  rep_ligand  <- toupper(manual_rep_ligand %||%
                           as.character(heat_df %>% arrange(desc(weighted_contribution)) %>%
                                          slice_head(n = 1) %>% pull(ligand)))
  rep_targets <- toupper(manual_rep_targets %||%
                           (heat_df %>% filter(ligand == rep_ligand) %>%
                              arrange(desc(weighted_contribution)) %>%
                              slice_head(n = 2) %>% pull(gene)))
  
  best_sec <- manual_best_section %||%
    (act_sec %>% filter(ligand == rep_ligand) %>%
       arrange(desc(auc)) %>% slice_head(n = 1) %>% pull(section_id))
  
  maps_h <- wrap_plots(
    lapply(rep_targets, function(tg)
      plot_feature_on_he(best_sec, tg,
                         paste0(best_sec, " | ", tg, " (target)"))),
    nrow = 1, guides = "collect"
  ) & theme(legend.position = "right")
  
  write_tsv(tibble(section = best_sec, ligand = rep_ligand,
                   targets = paste(rep_targets, collapse = ";")),
            file.path(out_dir, "Figure4h_spatial_targets_meta.tsv"))
  ggsave(file.path(out_dir, "Figure4h_spatial_targets.png"),
         maps_h, width = 9.0, height = 4.4, dpi = PLOT$dpi, bg = "white")
  ggsave(file.path(out_dir, "Figure4h_spatial_targets.pdf"),
         maps_h, width = 9.0, height = 4.4, device = cairo_pdf, bg = "white")
  message("Panel h done.")
  
}, error = function(e) {
  message("Panels g/h skipped — NicheNet inputs not found.\n  ", conditionMessage(e))
})

message("\nFigure 4 complete. Outputs saved to: ", out_dir)
