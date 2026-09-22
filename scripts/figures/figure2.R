# =========================================================
# Figure 2 — unified script
#
# Panels:
#   a  scRNA-seq UMAP coloured by manuscript cell group
#   b  % spots with non-zero lymphocyte NNLS coefficient across all 14 Visium sections
#   c  Lymphocyte hotspot spatial maps (6 retained sections, pink)
#   d  DeepTIL-inferred lymphoid composition of hotspot pseudobulks
#   e  Paired DeepTIL lymphoid abundance: hotspot vs non-hotspot pseudobulks
#
# Optional bonus:
#   b_plasma  Hotspot spatial maps coloured by cp10k (not in main fig)
# =========================================================

source(here::here("config.R"))

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(ggrepel)
  library(dbscan)
  library(patchwork)
  library(scales)
  library(jsonlite)
})

# =========================================================
# 0. Paths & global settings
# =========================================================

# Figure-specific scRNA-seq display inputs.
# `meta.tsv` under figure_inputs/figure2 is optional; if absent, use the
# single-cell metadata configured as REF$scrna_meta.
figure2_input_dir <- figure_input_path("figure2")
meta_candidate    <- file.path(figure2_input_dir, "meta.tsv")
meta_file         <- if (file.exists(meta_candidate)) meta_candidate else ref_file("scrna_meta")
umap_file         <- file.path(figure2_input_dir, "harmony_umap.coords.tsv.gz")

if (!file.exists(umap_file)) {
  stop(
    "Missing Figure 2 UMAP coordinates: ", umap_file,
    "\nPlace harmony_umap.coords.tsv.gz under EPN_FIGURE_INPUT_ROOT/figure2/. ",
    "See docs/data_sources.md."
  )
}

# Stage-derived inputs.
section_qc_file <- stage_input("semla_nnls", "qc/section_qc.tsv")
nnls_dir        <- file.path(stage_dir("hotspots"), "nnls")

# DeepTIL inputs (panels d, e).
deepTIL_results_dir <- stage_dir("deeptil")
deepTIL_file <- file.path(deepTIL_results_dir, DEEPTIL$output)

# Fallback: if the expected filename differs, accept a single SES_CIBERSORTx_*.txt.
if (!file.exists(deepTIL_file)) {
  deepTIL_hits <- list.files(
    deepTIL_results_dir,
    pattern = "^SES_CIBERSORTx_.*\\.txt$",
    full.names = TRUE
  )
  if (length(deepTIL_hits) == 0) {
    stop("No DeepTIL SES_CIBERSORTx_*.txt file found in: ", deepTIL_results_dir)
  }
  if (length(deepTIL_hits) > 1) {
    stop(
      "More than one SES_CIBERSORTx_*.txt file found. Set DEEPTIL$output ",
      "or the expected file explicitly.\n",
      paste(deepTIL_hits, collapse = "\n")
    )
  }
  deepTIL_file <- deepTIL_hits[1]
}

# Manuscript output directory.
out_dir <- ensure_dir(figure_path("figure2"))

# =========================================================
# 1. Section lists & shared settings
# =========================================================
good_sections <- GOOD_SECTIONS

# Preserve the historical lexicographic display order used in Figure 2
# without duplicating the cohort definition from config.R.
all_sections <- sort(ALL_SECTIONS)

# =========================================================
# 2. Colour palettes
# =========================================================
group_cols <- c(
  "Lymphocytes"               = "#F781BF",  # keep
  "Myeloid"                   = "#984EA3",  # keep
  "B cells"                   = "#FF7F00",
  "DC"                        = "#377EB8",
  "Neutrophils"               = "#A6CEE3",
  "Contaminating tumor cells" = "#4DAF4A"
)

lm7_cols <- c(
  "B"    = "#FF7F00",
  "CD4T" = "#4DAF4A",
  "CD8T" = "#377EB8",
  "gdT"  = "#E41A1C",
  "NK"   = "#984EA3"
)

# =========================================================
# 3. Helpers
# =========================================================
guess_id_col <- function(df) {
  hit <- c("cell_id","Cell","cell","barcode","Barcode","NAME","name")
  hit <- hit[hit %in% colnames(df)]
  if (length(hit)) hit[1] else colnames(df)[1]
}

guess_umap_cols <- function(df) {
  cn  <- colnames(df)
  xh  <- c("UMAP1","umap_1","UMAP_1","x","X")
  yh  <- c("UMAP2","umap_2","UMAP_2","y","Y")
  x   <- xh[xh %in% cn]; y <- yh[yh %in% cn]
  if (length(x) && length(y)) return(c(x[1], y[1]))
  num <- cn[sapply(df, is.numeric)]
  num[1:2]
}

read_hires_scale <- function(json_file) {
  as.numeric(jsonlite::fromJSON(json_file)$tissue_hires_scalef)
}

prep_section_coords <- function(df_sec, use_hires = TRUE) {
  if (use_hires) {
    if (!all(c("pxl_col_in_hires","pxl_row_in_hires") %in% colnames(df_sec))) {
      jf <- unique(df_sec$json_file)
      if (length(jf) != 1 || is.na(jf) || !file.exists(jf))
        stop("Cannot build hires coords for section ", unique(df_sec$section_id))
      sc <- read_hires_scale(jf)
      df_sec <- df_sec %>%
        mutate(pxl_col_in_hires = pxl_col_in_fullres * sc,
               pxl_row_in_hires = pxl_row_in_fullres * sc)
    }
    x_col <- "pxl_col_in_hires"; y_col <- "pxl_row_in_hires"
  } else {
    x_col <- "pxl_col_in_fullres"; y_col <- "pxl_row_in_fullres"
  }
  coords <- as.matrix(df_sec[, c(x_col, y_col)])
  pitch  <- median(dbscan::kNN(coords, k = 1)$dist[, 1], na.rm = TRUE)
  x_min  <- min(df_sec[[x_col]], na.rm = TRUE)
  y_max  <- max(df_sec[[y_col]], na.rm = TRUE)
  df_sec %>%
    mutate(
      spot_pitch = pitch,
      x_plot = (.data[[x_col]] - x_min) / pitch,
      y_plot = (y_max - .data[[y_col]]) / pitch
    )
}

# =========================================================
# PANEL A — scRNA-seq UMAP
# =========================================================
message("Panel a — UMAP")

meta_raw <- fread(meta_file)
umap_raw <- fread(umap_file)

meta_id  <- guess_id_col(meta_raw)
umap_id  <- guess_id_col(umap_raw)
umap_xy  <- guess_umap_cols(umap_raw)

meta_u <- meta_raw %>%
  rename(cell_id = all_of(meta_id))

umap_u <- umap_raw %>%
  rename(
    cell_id = all_of(umap_id),
    UMAP1   = all_of(umap_xy[1]),
    UMAP2   = all_of(umap_xy[2])
  )

df_umap <- inner_join(umap_u, meta_u, by = "cell_id") %>%
  mutate(
    cell_annot = trimws(as.character(cell_type)),
    annot_low  = str_to_lower(cell_annot),
    
    manuscript_group = case_when(
      
      # Lymphoid
      cell_annot %in% c(
        "Lymphocytes",
        "TRegs"
      ) ~ "Lymphocytes",
      
      # Myeloid continuum
      cell_annot %in% c(
        "Microglial",
        "Classical M1",
        "Alternative M2",
        "Hypoxia",
        "Mitotic",
        "Undefined M"
      ) ~ "Myeloid",
      
      # Dendritic cells
      cell_annot == "DC" ~ "DC",
      
      # Neutrophils
      cell_annot == "Neutrophil" ~ "Neutrophils",
      
      # B lineage — same manuscript group / same colour
      cell_annot %in% c(
        "B-cells",
        "Plasma-B-cells"
      ) ~ "B cells",
      
      # Tumor contamination
      cell_annot == "Contaminating Tumor cells" ~
        "Contaminating tumor cells",
      
      TRUE ~ "Other"
    ),
    
    # Separate DISPLAY labels for the UMAP only
    umap_label = case_when(
      cell_annot == "B-cells"        ~ "B cells",
      cell_annot == "Plasma-B-cells" ~ "Plasma cells",
      TRUE                           ~ manuscript_group
    )
  ) %>%
  
  filter(manuscript_group != "Other") %>%
  
  mutate(
    manuscript_group = factor(
      manuscript_group,
      levels = c(
        "Myeloid",
        "Lymphocytes",
        "B cells",
        "DC",
        "Neutrophils",
        "Contaminating tumor cells"
      )
    )
  )


# ---------------------------------------------------------
# Label positions
# B cells and Plasma cells are labelled separately,
# although both retain the same manuscript_group colour.
# ---------------------------------------------------------
label_df <- df_umap %>%
  group_by(umap_label) %>%
  summarise(
    UMAP1 = median(UMAP1),
    UMAP2 = median(UMAP2),
    .groups = "drop"
  ) %>%
  
  mutate(
    nudge_x = case_when(
      umap_label == "Lymphocytes"               ~ -0.8,
      umap_label == "Myeloid"                   ~ -0.5,
      umap_label == "B cells"                   ~  0.3,
      umap_label == "Plasma cells"              ~  0.2,
      umap_label == "DC"                        ~  0.3,
      umap_label == "Neutrophils"               ~  0.3,
      umap_label == "Contaminating tumor cells" ~  0.3,
      TRUE ~ 0
    ),
    
    nudge_y = case_when(
      umap_label == "Lymphocytes"               ~ -0.3,
      umap_label == "Myeloid"                   ~  0.4,
      umap_label == "B cells"                   ~ -0.2,
      umap_label == "Plasma cells"              ~ -0.2,
      umap_label == "DC"                        ~  0.3,
      umap_label == "Neutrophils"               ~  0.2,
      umap_label == "Contaminating tumor cells" ~ -0.2,
      TRUE ~ 0
    ),
    
    label_size = if_else(
      umap_label == "Lymphocytes",
      6,
      4
    )
  )


# ---------------------------------------------------------
# Plot
# ---------------------------------------------------------
p_a <- ggplot(df_umap, aes(UMAP1, UMAP2)) +
  
  # All non-lymphocyte populations
  geom_point(
    data = df_umap %>%
      filter(manuscript_group != "Lymphocytes"),
    aes(color = manuscript_group),
    size = 0.38,
    alpha = 0.85,
    stroke = 0
  ) +
  
  # Lymphocytes plotted last, slightly larger and fully opaque
  geom_point(
    data = df_umap %>%
      filter(manuscript_group == "Lymphocytes"),
    aes(color = manuscript_group),
    size = 0.55,
    alpha = 1,
    stroke = 0
  ) +
  
  scale_color_manual(
    values = group_cols,
    drop = FALSE
  ) +
  
  geom_text_repel(
    data = label_df,
    aes(
      UMAP1,
      UMAP2,
      label = umap_label,
      size = label_size
    ),
    inherit.aes = FALSE,
    nudge_x = label_df$nudge_x,
    nudge_y = label_df$nudge_y,
    fontface = "bold",
    color = "black",
    bg.color = NA,
    segment.color = NA,
    seed = 1
  ) +
  
  # Lymphocytes = 6, all other labels = 4
  scale_size_identity() +
  
  coord_equal() +
  
  labs(
    title = "scRNA-seq UMAP",
    x = "UMAP1",
    y = "UMAP2"
  ) +
  
  theme_classic(base_size = 15) +
  
  theme(
    legend.position = "none",
    
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 22
    ),
    
    axis.title = element_text(
      face = "bold",
      size = 16
    ),
    
    axis.text = element_blank(),
    
    axis.line = element_line(
      linewidth = 0.6,
      color = "black"
    ),
    
    axis.ticks = element_blank()
  )


# ---------------------------------------------------------
# Save
# ---------------------------------------------------------
ggsave(
  file.path(out_dir, "Figure2a_scRNAseq_UMAP.png"),
  p_a,
  width = 7.2,
  height = 5.2,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(out_dir, "Figure2a_scRNAseq_UMAP.pdf"),
  p_a,
  width = 7.2,
  height = 5.2,
  bg = "white"
)

# =========================================================
# PANELS B, C — cohort-wide NNLS QC + hotspot spatial maps
# =========================================================
message("Panels b/c — cohort-wide NNLS QC and Visium hotspot maps")

rds_files <- list.files(nnls_dir, pattern = "^NNLS_section_.*\\.rds$",
                        full.names = TRUE)
if (!length(rds_files)) stop("No NNLS section RDS files found in: ", nnls_dir)

df_all <- purrr::map_dfr(rds_files, readRDS)

miss <- setdiff(c("section_id","Lymphocytes","in_tissue","is_hotspot_semla"),
                colnames(df_all))
if (length(miss)) stop("Missing columns in NNLS tables: ", paste(miss, collapse = ", "))

df_good <- df_all %>%
  filter(section_id %in% good_sections, in_tissue == 1) %>%
  mutate(
    section_id = factor(section_id, levels = good_sections),
    is_hotspot = as.logical(is_hotspot_semla),
    Lymphocytes_cp10k_global = Lymphocytes / sum(Lymphocytes, na.rm = TRUE) * 10000
  )

df_vis <- df_good %>%
  split(.$section_id) %>%
  purrr::map_dfr(prep_section_coords)

x_lim_vis <- max(df_vis$x_plot, na.rm = TRUE)
y_lim_vis <- max(df_vis$y_plot, na.rm = TRUE)
ly_lim    <- quantile(df_vis$Lymphocytes_cp10k_global, 0.99, na.rm = TRUE)

BG_SZ  <- 0.38
HOT_SZ <- 0.50

# =========================================================
# PANEL B — % spots with non-zero lymphocyte NNLS coefficient
# =========================================================
# Stage 02 calculates prop_nonzero for all 14 sections as:
#   mean(Lymphocytes > 0, na.rm = TRUE)
#
# This panel therefore shows a measured cohort-wide QC quantity.
# It is NOT the percentage of hotspot spots.
#
# The dashed line indicates the automatic stage-02 QC threshold:
#   prop_nonzero >= 0.10
#
# Section 848 passed this numerical criterion but was subsequently excluded
# after spatial inspection because its lymphocyte NNLS signal was elevated
# near-uniformly across the capture area rather than showing a focal pattern.

message("Panel b — non-zero lymphocyte NNLS coefficient across all sections")

section_qc <- readr::read_tsv(
  section_qc_file,
  show_col_types = FALSE
)

required_qc_cols <- c("section_id", "prop_nonzero")

missing_qc_cols <- setdiff(
  required_qc_cols,
  colnames(section_qc)
)

if (length(missing_qc_cols)) {
  stop(
    "section_qc.tsv is missing required columns: ",
    paste(missing_qc_cols, collapse = ", ")
  )
}

panel_b_df <- section_qc %>%
  transmute(
    section_id = as.character(section_id),
    pct_nonzero = 100 * prop_nonzero
  ) %>%
  filter(section_id %in% all_sections) %>%
  tidyr::complete(section_id = all_sections) %>%
  mutate(
    status = if_else(
      section_id %in% good_sections,
      "Retained for hotspot analysis",
      "Not retained for hotspot analysis"
    ),
    status = factor(
      status,
      levels = c(
        "Retained for hotspot analysis",
        "Not retained for hotspot analysis"
      )
    ),
    section_id = factor(
      section_id,
      levels = all_sections
    )
  )

if (anyNA(panel_b_df$pct_nonzero)) {
  warning(
    "Missing stage-02 QC values for: ",
    paste(
      as.character(
        panel_b_df$section_id[is.na(panel_b_df$pct_nonzero)]
      ),
      collapse = ", "
    )
  )
}

# Automatic stage-02 threshold from config.R.
qc_prop_nonzero_threshold <- 100 * THRESH$min_prop_nonzero

p_b <- ggplot(
  panel_b_df,
  aes(
    x = section_id,
    y = pct_nonzero,
    fill = status
  )
) +
  geom_col(
    width = 0.82,
    color = "black",
    linewidth = 0.2
  ) +
  
  # Automatic QC threshold
  geom_hline(
    yintercept = qc_prop_nonzero_threshold,
    linetype = "dashed",
    linewidth = 0.5,
    colour = "grey40"
  ) +
  
  # Threshold label
  annotate(
    "text",
    x = 1.2,
    y = qc_prop_nonzero_threshold + 2.2,
    label = paste0("Automatic QC threshold (",
                   format(qc_prop_nonzero_threshold, trim = TRUE), "%)"),
    hjust = 0,
    vjust = 0,
    size = 3.5,
    colour = "grey35"
  ) +
  
  scale_fill_manual(
    values = c(
      "Retained for hotspot analysis"     = "#F781BF",
      "Not retained for hotspot analysis" = "grey75"
    ),
    name = NULL
  ) +
  
  scale_y_continuous(
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100),
    expand = expansion(mult = c(0, 0.02))
  ) +
  
  labs(
    x = NULL,
    y = "% spots with non-zero\nlymphocyte NNLS coefficient"
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.title.y = element_text(
      face = "bold"
    ),
    legend.position = "top",
    legend.justification = "left",
    plot.margin = margin(8, 8, 8, 8)
  )


# ---------------------------------------------------------
# Save
# ---------------------------------------------------------

ggsave(
  file.path(
    out_dir,
    "Figure2b_nonzero_lymphocyte_NNLS_all_sections.png"
  ),
  p_b,
  width = 7.8,
  height = 4.3,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(
    out_dir,
    "Figure2b_nonzero_lymphocyte_NNLS_all_sections.pdf"
  ),
  p_b,
  width = 7.8,
  height = 4.3,
  device = cairo_pdf,
  bg = "white"
)

write_tsv(
  panel_b_df,
  file.path(
    out_dir,
    "Figure2b_nonzero_lymphocyte_NNLS_all_sections.tsv"
  )
)

# Panel c — pink hotspot spatial maps (6 retained sections)
plot_panel_c <- function(sec) {
  df_sec <- df_vis %>% filter(section_id == sec)
  ggplot() +
    geom_point(data = df_sec, aes(x_plot, y_plot),
               color = "grey82", size = BG_SZ) +
    geom_point(data = df_sec %>% filter(is_hotspot),
               aes(x_plot, y_plot), color = "#F781BF", size = HOT_SZ) +
    coord_equal(xlim = c(0, x_lim_vis), ylim = c(0, y_lim_vis), expand = FALSE) +
    ggtitle(as.character(sec)) +
    theme_void(base_size = 12) +
    theme(plot.title  = element_text(hjust = 0.5, face = "bold", size = 15),
          plot.margin = margin(6, 6, 6, 6))
}

fig_c <- wrap_plots(purrr::map(good_sections, plot_panel_c), nrow = 1)

ggsave(file.path(out_dir, "Figure2c_hotspot_maps_6sections.png"),
       fig_c, width = 15.0, height = 4.2, dpi = 600, bg = "white")
ggsave(file.path(out_dir, "Figure2c_hotspot_maps_6sections.pdf"),
       fig_c, width = 15.0, height = 4.2, device = cairo_pdf, bg = "white")

# Bonus: panel c coloured by cp10k (not in main figure, useful for QC)
plot_plasma <- function(sec, show_legend = FALSE) {
  df_sec <- df_vis %>% filter(section_id == sec)
  p <- ggplot() +
    geom_point(data = df_sec, aes(x_plot, y_plot),
               color = "grey82", size = BG_SZ) +
    geom_point(data = df_sec %>% filter(is_hotspot),
               aes(x_plot, y_plot, color = Lymphocytes_cp10k_global),
               size = HOT_SZ) +
    scale_color_viridis_c(option = "plasma", limits = c(0, ly_lim),
                          oob  = scales::squish, name = "Lymphocyte\ncp10k",
                          guide = guide_colorbar(barheight = unit(35,"mm"),
                                                 barwidth  = unit(6,"mm"))) +
    coord_equal(xlim = c(0, x_lim_vis), ylim = c(0, y_lim_vis), expand = FALSE) +
    ggtitle(as.character(sec)) +
    theme_void(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 15))
  if (!show_legend) p <- p + theme(legend.position = "none")
  p
}

fig_plasma <- wrap_plots(
  purrr::map(good_sections,
             ~ plot_plasma(.x, show_legend = (.x == tail(good_sections, 1)))),
  nrow = 1, guides = "collect"
) & theme(legend.position = "right")

ggsave(file.path(out_dir, "Figure2_bonus_hotspot_maps_plasma_cp10k.png"),
       fig_plasma, width = 15.0, height = 4.2, dpi = 600, bg = "white")
ggsave(file.path(out_dir, "Figure2_bonus_hotspot_maps_plasma_cp10k.pdf"),
       fig_plasma, width = 15.0, height = 4.2, device = cairo_pdf, bg = "white")

# =========================================================
# PANELS D, E — DeepTIL lymphoid composition and abundance
# =========================================================
message("Panels d/e — DeepTIL lymphoid composition and paired abundance")

# DeepTIL output contains inferred abundance values (arbitrary units) for the
# seven LM7 populations. Here we retain the five lymphoid populations.
deepTIL_raw <- read.delim(
  deepTIL_file,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

deepTIL_required <- c("Mixture", "Bcells", "TCD4", "TCD8", "Tgd", "NK")
deepTIL_missing  <- setdiff(deepTIL_required, colnames(deepTIL_raw))
if (length(deepTIL_missing)) {
  stop(
    "DeepTIL file is missing required columns: ",
    paste(deepTIL_missing, collapse = ", ")
  )
}

deepTIL_df <- deepTIL_raw %>%
  mutate(
    section = sub("_(background|hotspot)$", "", Mixture),
    group   = sub("^.*_(background|hotspot)$", "\\1", Mixture)
  ) %>%
  filter(section %in% good_sections) %>%
  mutate(
    section = factor(section, levels = good_sections),
    group   = factor(group, levels = c("background", "hotspot"))
  )

# Require one matched hotspot/non-hotspot pseudobulk per section.
deepTIL_pair_qc <- deepTIL_df %>%
  count(section, group, name = "n") %>%
  tidyr::complete(section, group, fill = list(n = 0))

if (any(deepTIL_pair_qc$n != 1)) {
  print(deepTIL_pair_qc)
  stop("Expected exactly one background and one hotspot DeepTIL sample per section.")
}

deepTIL_long <- deepTIL_df %>%
  select(Mixture, section, group, Bcells, TCD4, TCD8, Tgd, NK) %>%
  pivot_longer(
    cols = c(Bcells, TCD4, TCD8, Tgd, NK),
    names_to = "population_raw",
    values_to = "Abundance"
  ) %>%
  mutate(
    Immune_pop = recode(
      population_raw,
      "Bcells" = "B",
      "TCD4"   = "CD4T",
      "TCD8"   = "CD8T",
      "Tgd"    = "gdT",
      "NK"     = "NK"
    ),
    Immune_pop = factor(Immune_pop, levels = names(lm7_cols))
  )

write_tsv(
  deepTIL_long,
  file.path(out_dir, "Figure2de_DeepTIL_lymphoid_abundance_long.tsv")
)

# ---------------------------------------------------------
# PANEL D — hotspot-only lymphoid composition
# ---------------------------------------------------------
# DeepTIL abundances are normalised within the five-population lymphoid
# compartment separately for each hotspot pseudobulk. These percentages are
# descriptive composition of the LM7/DeepTIL-inferred lymphoid signal, not
# absolute cell fractions.
panel_d_df <- deepTIL_long %>%
  filter(group == "hotspot") %>%
  group_by(section) %>%
  mutate(
    total_lymphoid = sum(Abundance, na.rm = TRUE),
    pct_lymphoid   = if_else(
      total_lymphoid > 0,
      100 * Abundance / total_lymphoid,
      NA_real_
    )
  ) %>%
  ungroup()

write_tsv(
  panel_d_df,
  file.path(out_dir, "Figure2d_DeepTIL_hotspot_lymphoid_composition.tsv")
)

p_d <- ggplot(
  panel_d_df,
  aes(x = section, y = pct_lymphoid, fill = Immune_pop)
) +
  geom_col(width = 0.78, color = "black", linewidth = 0.20) +
  scale_fill_manual(values = lm7_cols, drop = FALSE, name = NULL) +
  scale_y_continuous(
    breaks = c(0, 25, 50, 75, 100),
    expand = c(0, 0),
    name = "% hotspot lymphoid abundance"
  ) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(x = NULL) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x     = element_text(size = 11),
    axis.title.y    = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave(
  file.path(out_dir, "Figure2d_DeepTIL_hotspot_lymphoid_composition.png"),
  p_d, width = 7.6, height = 4.4, dpi = 600, bg = "white"
)
ggsave(
  file.path(out_dir, "Figure2d_DeepTIL_hotspot_lymphoid_composition.pdf"),
  p_d, width = 7.6, height = 4.4, device = cairo_pdf, bg = "white"
)

# ---------------------------------------------------------
# PANEL E — paired hotspot vs non-hotspot abundance
# ---------------------------------------------------------
panel_e_df <- deepTIL_long %>%
  mutate(
    group_plot = factor(
      as.character(group),
      levels = c("background", "hotspot"),
      labels = c("Non-hotspot", "Hotspot")
    )
  )

# Paired section-level statistics.
# The displayed p-value is the unadjusted paired Wilcoxon p-value.
# BH-adjusted values across the five lymphoid populations are also exported.
fmt_p <- function(p) {
  dplyr::case_when(
    is.na(p) ~ NA_character_,
    p < 0.001 ~ format(p, scientific = TRUE, digits = 2),
    TRUE ~ formatC(p, format = "f", digits = 2)
  )
}

panel_e_stats <- panel_e_df %>%
  select(section, group_plot, Immune_pop, Abundance) %>%
  pivot_wider(names_from = group_plot, values_from = Abundance) %>%
  group_by(Immune_pop) %>%
  summarise(
    n_sections = sum(complete.cases(`Non-hotspot`, Hotspot)),
    n_up = sum(Hotspot > `Non-hotspot`, na.rm = TRUE),
    max_abundance = max(c(`Non-hotspot`, Hotspot), na.rm = TRUE),
    p_value = {
      ok <- complete.cases(`Non-hotspot`, Hotspot)
      if (sum(ok) >= 2 && any(Hotspot[ok] != `Non-hotspot`[ok])) {
        suppressWarnings(
          wilcox.test(
            Hotspot[ok],
            `Non-hotspot`[ok],
            paired = TRUE,
            exact = FALSE
          )$p.value
        )
      } else {
        NA_real_
      }
    },
    .groups = "drop"
  ) %>%
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH"),
    y_annot = if_else(max_abundance > 0, max_abundance * 1.13, 0.01),
    label = if_else(
      is.na(p_value),
      paste0(n_up, "/", n_sections, " \u2191"),
      paste0(n_up, "/", n_sections, " \u2191\np = ", fmt_p(p_value))
    )
  )

write_tsv(
  panel_e_stats,
  file.path(out_dir, "Figure2e_DeepTIL_hotspot_vs_nonhotspot_stats.tsv")
)

p_e <- ggplot(
  panel_e_df,
  aes(x = group_plot, y = Abundance, group = section)
) +
  geom_line(
    colour = "grey70",
    linewidth = 0.55,
    alpha = 0.9
  ) +
  geom_point(
    aes(fill = Immune_pop),
    shape = 21,
    size = 3.0,
    stroke = 0.35,
    colour = "black"
  ) +
  stat_summary(
    aes(group = group_plot),
    fun = median,
    geom = "crossbar",
    width = 0.18,
    linewidth = 0.65,
    colour = "black"
  ) +
  geom_text(
    data = panel_e_stats,
    aes(x = 1.5, y = y_annot, label = label),
    inherit.aes = FALSE,
    size = 3.5,
    lineheight = 1.0
  ) +
  facet_wrap(
    ~ Immune_pop,
    scales = "free_y",
    nrow = 1
  ) +
  scale_fill_manual(values = lm7_cols, drop = FALSE) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.20))
  ) +
  labs(
    x = NULL,
    y = "DeepTIL inferred abundance (a.u.)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.title.y = element_text(
      face = "bold",
      size = 12,
      margin = margin(r = 7)
    ),
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 8.5),
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_blank(),
    panel.spacing.x = grid::unit(1.0, "lines"),
    legend.position = "none"
  )

ggsave(
  file.path(out_dir, "Figure2e_DeepTIL_lymphoid_hotspot_vs_nonhotspot.png"),
  p_e, width = 13.0, height = 3.8, dpi = 600, bg = "white"
)
ggsave(
  file.path(out_dir, "Figure2e_DeepTIL_lymphoid_hotspot_vs_nonhotspot.pdf"),
  p_e, width = 13.0, height = 3.8, device = cairo_pdf, bg = "white"
)

message("Figure 2 — all panels done. Outputs saved to: ", out_dir)

