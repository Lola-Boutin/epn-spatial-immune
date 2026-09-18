# =============================================================================
# 09_opportunity_map.R
#
# Purpose
# -------
# Immunotherapy opportunity map: Figure 5a (opportunity quadrant) and Figure 5b
# (mechanism axes), plus the supplementary tables documenting what was scored.
#
# Self-contained: gene programs, scoring and classification are all defined
# here. Script 10 reads this script's output and never redefines a gene set, so
# the two cannot drift apart.
#
# DESIGN
# Gene programs state the biological hypothesis and are NOT edited per dataset.
# Genes are removed only by the detection filter, which is data-driven and
# reported, so the same lists apply unchanged to another cohort.
#
# Composites average MECHANISM AXES rather than genes, so axis gene count does
# not set its weight. This is the right construction for a formative index,
# where members define the construct rather than reflecting a latent cause.
#
# Scores are z-transformed across the POOLED tumor spots, not within section:
# between-section differences are the quantity of interest. A consequence is
# that they reflect tumor biology together with inter-patient and technical
# variation, which this design cannot separate.
#
# Inputs
# ------
#   stage_input("manifest", "visium_manifest.tsv")
#   stage_dir("raw_vis")     vis_section_<sec>_raw.rds
#   stage_dir("zones_all")   tables/  -- ALL 14 sections (see docs/pipeline.md)
#
# Outputs
# -------
#   stage_dir("scoring")
#     tables/09_spot_scores.tsv | rds/09_spot_scores.rds
#     tables/09_section_summary.tsv | rds/09_section_summary.rds
#     tables/09_section_by_zone_summary.tsv
#     tables/09_gene_detection.tsv, 09_gene_detection_by_section.tsv
#     tables/09_gene_provenance.tsv, 09_effective_gene_sets.tsv
#     tables/09_primary_relapse_pairs.tsv, 09_core_revision_comparison.tsv
#     rds/09_scoring_inputs.rds   <- cache consumed by 09b
#     plots/Fig5A_opportunity_quadrant.{png,pdf}
#     plots/Fig5B_gdT_mechanism_axes.{png,pdf}
#     plots/S_inhibitory_mechanism_axes.{png,pdf}
#     plots/S_inhibitory_tone_by_section.{png,pdf}
#     plots/S_primary_relapse_paired.{png,pdf}
#
# Stochastic
# ----------
#   none
#
# Runtime
# -------
#   ~40 minutes
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
  library(ggrepel)
  library(patchwork)
})


# -------------------------
# CONFIG
# -------------------------
manifest_file  <- stage_input("manifest", "visium_manifest.tsv")
raw_vis_dir    <- stage_dir("raw_vis")
phase06_tables <- file.path(stage_dir("zones_all"), "tables")

out_root   <- stage_dir("scoring")
out_tables <- ensure_dir(file.path(out_root, "tables"))
out_plots  <- ensure_dir(file.path(out_root, "plots"))
out_rds    <- ensure_dir(file.path(out_root, "rds"))

zones_keep <- c("Epithelial", "Mesenchymal")

DETECTION_MIN_PCT   <- THRESH$detection_floor * 100  # % of tumor spots
SUPPRESSED_INHIB_Z  <- 0.75     # median inhibitory z defining elevated tone
SCORING_METHOD      <- "zmean"  # "zmean" (primary) or "rawmean" (sensitivity)
SCORING_AGGREGATION <- "axis"   # "axis" (equal weight per mechanism) or "gene"

manual_section_metadata <- tibble::tribble(
  ~section_id, ~patient_id, ~sample_status,
  "459",   "P459", "Primary",
  "459_2", "P459", "Relapse",
  "723",   "P723", "Primary",
  "723_2", "P723", "Relapse",
  "928",   "P928", "Primary",
  "928_2", "P928", "Relapse"
)

# =========================================================
# GENE PROGRAMS
#
# "A|B" marks alias alternatives for one gene; whichever symbol the reference
# carries is used (the better-detected one if both are present).
# =========================================================
PROGRAMS <- list(

  # MHC class I antigen presentation machinery and immunoproteasome:
  # what MHC-I-restricted alphabeta T cell recognition requires.
  abT_core = c(
    "HLA-A", "HLA-B", "HLA-C", "B2M",
    "TAP1", "TAP2", "PSMB8", "PSMB9",
    "ERAP1", "ERAP2", "NLRC5", "IRF1"
  ),

  # gammadelta LIGAND AVAILABILITY across four recognition modalities plus the
  # phosphoantigen biosynthetic pathway. This is not a phosphoantigen-only
  # score - name it accordingly in the manuscript.
  #
  # Mevalonate is included because the core otherwise contained the presenting
  # molecules (BTN2A1/BTN3A1) but no measure of the antigen they present:
  # mevalonate flux produces IPP, IPP binding the BTN3A1 B30.2 domain drives
  # the conformational change Vgamma9Vdelta2 TCRs read, and BTN2A1 is required
  # for that recognition.
  #
  # IDI1, not IDI2 - IDI2 is largely skeletal-muscle restricted.
  # FDPS is deliberately absent: it CONSUMES IPP, so it runs opposite in
  # direction to the rest of the pathway. It is scored separately (MODIFIERS).
  gdT_core = c(
    "BTN2A1", "BTN3A1",
    "MVK", "PMVK", "MVD", "IDI1",
    "MICA", "MICB", "ULBP2", "ULBP3",
    "ULBP4|RAET1E", "ULBP5|RAET1G", "ULBP6|RAET1L",
    "EPHA2",
    "PVR", "NECTIN2", "ICAM1"
  ),

  # Immunosuppressive tone. NOT a checkpoint score - report it as what is
  # actually detected, which in this cohort excludes CD274 / PDCD1LG2 / IDO1.
  inhibitory_program = c(
    "CD274", "PDCD1LG2",
    "HLA-E",
    "LGALS9",
    "IDO1",
    "TGFB1", "TGFB2", "TGFB3",
    "IL10",
    "VEGFA",
    "NT5E", "ENTPD1",
    "PTGS2", "PTGES", "PTGER2", "PTGER4"
  )
)

# Pre-revision gdT core (no mevalonate), scored alongside so the effect of the
# revision can be reported rather than asserted.
GDT_CORE_LEGACY <- c(
  "BTN2A1", "BTN3A1", "EPHA2",
  "MICA", "MICB", "ULBP2", "ULBP3",
  "ULBP4|RAET1E", "ULBP5|RAET1G", "ULBP6|RAET1L",
  "PVR", "NECTIN2", "ICAM1"
)

# Mechanism axes partition the composites above. With SCORING_AGGREGATION =
# "axis" they also set the composite weighting.
AXES <- list(
  gdT_btn            = c("BTN2A1", "BTN3A1"),
  gdT_mevalonate     = c("MVK", "PMVK", "MVD", "IDI1"),
  gdT_nkg2d          = c("MICA", "MICB", "ULBP2", "ULBP3",
                         "ULBP4|RAET1E", "ULBP5|RAET1G", "ULBP6|RAET1L"),
  gdT_ephrin         = c("EPHA2"),
  gdT_adhesion_dnam  = c("PVR", "NECTIN2", "ICAM1"),

  inh_mhc_checkpoint = c("HLA-E"),
  inh_classical_ckpt = c("CD274", "PDCD1LG2", "IDO1", "LGALS9"),
  inh_tgfb           = c("TGFB1", "TGFB2", "TGFB3"),
  inh_adenosine      = c("NT5E", "ENTPD1"),
  inh_pge2           = c("PTGS2", "PTGES", "PTGER2", "PTGER4"),
  inh_angiogenic     = c("VEGFA"),
  inh_il10           = c("IL10")
)

PROGRAM_AXES <- list(
  abT_core           = NULL,   # single coherent construct, not partitioned
  gdT_core           = c("gdT_btn", "gdT_mevalonate", "gdT_nkg2d",
                         "gdT_ephrin", "gdT_adhesion_dnam"),
  inhibitory_program = c("inh_mhc_checkpoint", "inh_classical_ckpt", "inh_tgfb",
                         "inh_adenosine", "inh_pge2", "inh_angiogenic", "inh_il10")
)

# Scored and reported, never added to a composite. FDPS consumes IPP, so a
# section high in mevalonate flux AND high in FDPS is where zoledronate
# (an FDPS inhibitor) would have the most to act on.
MODIFIERS <- list(
  gdT_ipp_drain = c("FDPS")
)

# -------------------------
# LABELS
# -------------------------
OPPORTUNITY_LEVELS <- c("abT-favorable", "Dual-opportunity",
                        "gdT-favorable", "Immuno-cold", "Suppressed")

OPPORTUNITY_LABELS <- c(
  "abT-favorable"    = "\u03b1\u03b2T-favorable",
  "Dual-opportunity" = "Dual-opportunity",
  "gdT-favorable"    = "\u03b3\u03b4T-favorable",
  "Immuno-cold"      = "Immuno-cold",
  "Suppressed"       = "Suppressed"
)

OPPORTUNITY_COLORS <- c(
  "abT-favorable"    = "#F8766D",
  "Dual-opportunity" = "#A3A500",
  "gdT-favorable"    = "#00BF7D",
  "Immuno-cold"      = "#00B0F6",
  "Suppressed"       = "#E76BF3"
)

AXIS_LAB_AB <- "\u03b1\u03b2T net"
AXIS_LAB_GD <- "\u03b3\u03b4T net"

AXIS_LABELS <- c(
  gdT_btn            = "Butyrophilin (BTN2A1/3A1)",
  gdT_mevalonate     = "Mevalonate / IPP synthesis",
  gdT_nkg2d          = "NKG2D ligands",
  gdT_ephrin         = "EPHA2",
  gdT_adhesion_dnam  = "DNAM-1 / adhesion",
  gdT_ipp_drain      = "FDPS (IPP consumption)",
  inh_mhc_checkpoint = "HLA-E / NKG2A",
  inh_classical_ckpt = "Classical checkpoints",
  inh_tgfb           = "TGF-beta",
  inh_adenosine      = "Adenosine",
  inh_pge2           = "PGE2 / COX",
  inh_angiogenic     = "VEGF",
  inh_il10           = "IL-10"
)

axis_label <- function(x) ifelse(x %in% names(AXIS_LABELS), AXIS_LABELS[x], x)

# =========================================================
# HELPERS
# =========================================================
z_vec <- function(x) {
  ok <- is.finite(x)
  if (sum(ok) < 2) return(rep(NA_real_, length(x)))
  s <- stats::sd(x[ok])
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mean(x[ok])) / s
}

mean_or_na   <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
median_or_na <- function(x) if (all(is.na(x))) NA_real_ else stats::median(x, na.rm = TRUE)

safe_scale <- function(X) {
  mu  <- colMeans(X, na.rm = TRUE)
  sdv <- apply(X, 2, stats::sd, na.rm = TRUE)
  sdv[!is.finite(sdv) | sdv == 0] <- NA_real_
  Z <- sweep(sweep(X, 2, mu, "-"), 2, sdv, "/")
  Z[!is.finite(Z)] <- 0
  Z
}

extract_barcode16 <- function(x) {
  x <- toupper(as.character(x))
  m <- regexpr("[ACGT]{16}", x)
  out <- rep(NA_character_, length(x)); ok <- m > 0
  out[ok] <- regmatches(x, m); out
}

read_expr_matrix <- function(obj) {
  DefaultAssay(obj) <- "Spatial"
  mat <- tryCatch(GetAssayData(obj, assay = "Spatial", layer = "data"),
                  error = function(e) NULL)
  if (is.null(mat) || nrow(mat) == 0 || ncol(mat) == 0) {
    message("  no usable 'data' layer; normalising from counts")
    obj <- NormalizeData(obj, normalization.method = "LogNormalize",
                         scale.factor = 1e4, verbose = FALSE)
    mat <- GetAssayData(obj, assay = "Spatial", layer = "data")
  }
  if (is.null(mat) || nrow(mat) == 0) stop("Could not obtain normalized matrix")
  rn <- toupper(rownames(mat)); keep <- !duplicated(rn)
  mat <- mat[keep, , drop = FALSE]; rownames(mat) <- rn[keep]
  list(obj = obj, mat = mat)
}

# Alias resolution: "A|B" -> whichever symbol exists (better detected if both)
resolve_symbol <- function(entry, present_symbols, detection_pct = NULL) {
  opts  <- toupper(trimws(strsplit(entry, "|", fixed = TRUE)[[1]]))
  avail <- opts[opts %in% present_symbols]
  if (length(avail) == 0) return(NA_character_)
  if (length(avail) == 1) return(avail)
  if (is.null(detection_pct)) return(avail[1])
  avail[which.max(detection_pct[match(avail, names(detection_pct))])]
}

build_effective_sets <- function(all_sets, present_symbols, detection_pct,
                                 min_pct = DETECTION_MIN_PCT) {
  prov <- list(); eff <- list()
  for (set_name in names(all_sets)) {
    tb <- bind_rows(lapply(all_sets[[set_name]], function(entry) {
      sym <- resolve_symbol(entry, present_symbols, detection_pct)
      pct <- if (!is.na(sym)) unname(detection_pct[sym]) else NA_real_
      tibble(set = set_name, declared = entry, resolved_symbol = sym,
             pct_detected = pct,
             status = dplyr::case_when(
               is.na(sym)      ~ "absent from reference",
               !is.finite(pct) ~ "absent from reference",
               pct >= min_pct  ~ "scored",
               TRUE            ~ paste0("below ", min_pct, "% detection")))
    }))
    prov[[set_name]] <- tb
    eff[[set_name]]  <- tb$resolved_symbol[tb$status == "scored"]
  }
  list(sets = eff, provenance = bind_rows(prov))
}

score_set <- function(X, method = SCORING_METHOD) {
  if (is.null(X) || ncol(X) == 0) return(rep(NA_real_, nrow(X)))
  X[is.na(X)] <- 0
  switch(method,
         zmean   = rowMeans(safe_scale(X)),
         rawmean = rowMeans(X),
         stop("Unknown scoring method: ", method))
}

# Composite score honouring SCORING_AGGREGATION
score_program <- function(expr_mat, program, sets,
                          aggregation = SCORING_AGGREGATION,
                          method = SCORING_METHOD) {
  ax <- PROGRAM_AXES[[program]]
  if (aggregation == "gene" || is.null(ax)) {
    g <- sets[[program]]
    if (length(g) == 0) return(rep(NA_real_, nrow(expr_mat)))
    return(score_set(expr_mat[, g, drop = FALSE], method))
  }
  if (aggregation != "axis") stop("Unknown aggregation: ", aggregation)
  ax_use <- ax[sapply(ax, function(a) length(sets[[a]]) > 0)]
  if (length(ax_use) == 0) return(rep(NA_real_, nrow(expr_mat)))
  M <- sapply(ax_use, function(a) z_vec(score_set(expr_mat[, sets[[a]], drop = FALSE], method)))
  if (is.null(dim(M))) M <- matrix(M, ncol = 1)
  rowMeans(M, na.rm = TRUE)
}

classify_opportunity <- function(ab_net, gd_net, inhibitory_z,
                                 threshold = SUPPRESSED_INHIB_Z) {
  if (any(is.na(inhibitory_z))) {
    stop("inhibitory_z contains NA - sections would fall through to Immuno-cold.")
  }
  out <- dplyr::case_when(
    ab_net >  0 & gd_net >  0 ~ "Dual-opportunity",
    ab_net >  0 & gd_net <= 0 ~ "abT-favorable",
    ab_net <= 0 & gd_net >  0 ~ "gdT-favorable",
    ab_net <= 0 & gd_net <= 0 & inhibitory_z >= threshold ~ "Suppressed",
    TRUE ~ "Immuno-cold"
  )
  factor(out, levels = OPPORTUNITY_LEVELS)
}

display_category <- function(x) {
  factor(unname(OPPORTUNITY_LABELS[as.character(x)]),
         levels = unname(OPPORTUNITY_LABELS[OPPORTUNITY_LEVELS]))
}

run_spearman_safe <- function(df, x, y) {
  dd <- df %>% filter(is.finite(.data[[x]]), is.finite(.data[[y]]))
  if (nrow(dd) < 4) return(tibble(x = x, y = y, rho = NA_real_, p_value = NA_real_, n = nrow(dd)))
  ct <- suppressWarnings(cor.test(dd[[x]], dd[[y]], method = "spearman", exact = FALSE))
  tibble(x = x, y = y, rho = unname(ct$estimate), p_value = ct$p.value, n = nrow(dd))
}

# =========================================================
# LOAD INPUTS
# =========================================================
manifest <- read_tsv(manifest_file, show_col_types = FALSE) %>%
  filter(complete %in% TRUE) %>% arrange(section_id)

zones_all <- read_tsv(file.path(phase06_tables, "zones_all_sections.tsv"),
                      show_col_types = FALSE)

ALL_SETS <- c(PROGRAMS, AXES, MODIFIERS, list(gdT_core_legacy = GDT_CORE_LEGACY))
all_symbols <- unique(toupper(unlist(strsplit(unlist(ALL_SETS), "|", fixed = TRUE))))

# =========================================================
# PASS 1 — collect program-gene expression across tumor spots
#
# Detection rates must be known before scoring, because the detection floor
# decides which genes enter each score. Expression is therefore gathered first
# and scored afterwards.
# =========================================================
mat_list <- list(); meta_list <- list(); present_list <- list()

for (sec in manifest$section_id) {
  message("\n=== Section ", sec, " ===")

  vis_file <- file.path(raw_vis_dir, paste0("vis_section_", sec, "_raw.rds"))
  if (!file.exists(vis_file)) { warning("Missing vis file: ", sec); next }

  zone_sec <- zones_all %>%
    filter(section_id == sec, in_tissue == 1, zone_call %in% zones_keep)
  if (nrow(zone_sec) == 0) { warning("No tumor spots: ", sec); next }

  vis_raw  <- readRDS(vis_file)
  expr     <- read_expr_matrix(vis_raw)
  expr_mat <- expr$mat

  keep_cells <- intersect(zone_sec$cell, colnames(expr_mat))
  if (length(keep_cells) == 0) {
    md_vis <- tibble(cell = colnames(expr_mat),
                     barcode16_join = coalesce(as.character(expr$obj@meta.data$barcode16),
                                               extract_barcode16(colnames(expr_mat))))
    zone_match <- zone_sec %>%
      mutate(barcode16_join = coalesce(as.character(barcode16_join),
                                       as.character(barcode16),
                                       extract_barcode16(cell))) %>%
      left_join(md_vis, by = "barcode16_join", suffix = c(".zone", ".vis"))
    keep_cells <- unique(zone_match$cell.vis[!is.na(zone_match$cell.vis)])
  }
  keep_cells <- intersect(keep_cells, colnames(expr_mat))
  if (length(keep_cells) == 0) { warning("No matched cells: ", sec); next }

  expr_use <- expr_mat[, keep_cells, drop = FALSE]
  zone_use <- zone_sec %>%
    filter(cell %in% keep_cells) %>%
    distinct(cell, .keep_all = TRUE) %>%
    slice(match(colnames(expr_use), cell))

  hit <- intersect(all_symbols, rownames(expr_use))
  present_list[[sec]] <- tibble(section_id = sec, gene = hit)

  sub  <- as.matrix(t(expr_use[hit, , drop = FALSE]))
  full <- matrix(NA_real_, nrow = nrow(sub), ncol = length(all_symbols),
                 dimnames = list(paste(sec, rownames(sub), sep = "|"), all_symbols))
  full[, hit] <- sub

  mat_list[[sec]]  <- full
  meta_list[[sec]] <- tibble(key = rownames(full), section_id = sec,
                             cell = rownames(sub), zone_call = zone_use$zone_call)
  rm(vis_raw, expr, expr_mat, expr_use, sub); gc(verbose = FALSE)
}

expr_prog <- do.call(rbind, mat_list)
meta_all  <- bind_rows(meta_list)
present_genes <- bind_rows(present_list)
if (nrow(expr_prog) == 0) stop("No tumor spots collected.")
stopifnot(nrow(expr_prog) == nrow(meta_all))

message("\nCollected ", nrow(expr_prog), " tumor spots across ",
        n_distinct(meta_all$section_id), " sections")

# =========================================================
# DETECTION, ALIAS RESOLUTION, EFFECTIVE GENE SETS
# =========================================================
detection <- tibble(
  gene = colnames(expr_prog),
  pct_detected = 100 * colMeans(expr_prog > 0, na.rm = TRUE),
  mean_lognorm = colMeans(expr_prog, na.rm = TRUE)
) %>% mutate(across(c(pct_detected, mean_lognorm), ~ ifelse(is.finite(.x), .x, 0)))

detection_by_section <- bind_rows(lapply(unique(meta_all$section_id), function(sec) {
  idx <- meta_all$section_id == sec
  tibble(section_id = sec, gene = colnames(expr_prog),
         pct_detected = 100 * colMeans(expr_prog[idx, , drop = FALSE] > 0, na.rm = TRUE))
})) %>% mutate(pct_detected = ifelse(is.finite(pct_detected), pct_detected, 0))

det_vec <- setNames(detection$pct_detected, detection$gene)
eff <- build_effective_sets(ALL_SETS, unique(present_genes$gene), det_vec)
gene_sets <- eff$sets

write_tsv(detection,            file.path(out_tables, "09_gene_detection.tsv"))
write_tsv(detection_by_section, file.path(out_tables, "09_gene_detection_by_section.tsv"))
write_tsv(eff$provenance,       file.path(out_tables, "09_gene_provenance.tsv"))

message("\n--- Genes not scored ---")
print(eff$provenance %>% filter(status != "scored") %>%
        select(set, declared, resolved_symbol, pct_detected, status), n = Inf)

set_sizes <- tibble(
  set = names(gene_sets),
  n_declared = lengths(ALL_SETS[names(gene_sets)]),
  n_scored   = lengths(gene_sets),
  genes_scored = sapply(gene_sets, paste, collapse = ";")
)
write_tsv(set_sizes, file.path(out_tables, "09_effective_gene_sets.tsv"))
message("\n--- Effective gene sets ---")
print(set_sizes %>% select(-genes_scored), n = Inf)

for (nm in names(PROGRAMS)) {
  if (length(gene_sets[[nm]]) == 0) stop("No genes passed the detection filter for: ", nm)
}

# Cache everything 09b_scoring_qc.R needs, so the QC script never redefines a
# gene set and never repeats the expression pass.
saveRDS(list(expr = expr_prog,
             meta = meta_all,
             gene_sets = gene_sets,
             program_axes = PROGRAM_AXES,
             programs = names(PROGRAMS),
             axes = names(AXES),
             detection = detection,
             config = list(scoring_method = SCORING_METHOD,
                           scoring_aggregation = SCORING_AGGREGATION,
                           detection_min_pct = DETECTION_MIN_PCT,
                           suppressed_inhib_z = SUPPRESSED_INHIB_Z,
                           zones_keep = zones_keep)),
        file.path(out_rds, "09_scoring_inputs.rds"))

# =========================================================
# PASS 2 — score
# =========================================================
spot <- meta_all

# composites (axis-weighted)
for (nm in names(PROGRAMS)) {
  spot[[nm]] <- score_program(expr_prog, nm, gene_sets)
}

# individual mechanism axes and modifiers (always gene-level)
for (nm in c(names(AXES), names(MODIFIERS))) {
  g <- gene_sets[[nm]]
  spot[[nm]] <- if (length(g) == 0) NA_real_ else score_set(expr_prog[, g, drop = FALSE])
}

# legacy core, for the revision comparison
spot$gdT_core_legacy <- if (length(gene_sets$gdT_core_legacy) == 0) NA_real_ else
  score_set(expr_prog[, gene_sets$gdT_core_legacy, drop = FALSE])

score_cols <- c(names(PROGRAMS), names(AXES), names(MODIFIERS), "gdT_core_legacy")
score_cols <- score_cols[sapply(score_cols, function(c) !all(is.na(spot[[c]])))]

# global z across pooled tumor spots (deliberate - see header)
for (nm in score_cols) spot[[paste0(nm, "_z")]] <- z_vec(spot[[nm]])

spot <- spot %>%
  mutate(
    abT_net    = abT_core_z - inhibitory_program_z,
    gdT_net    = gdT_core_z - inhibitory_program_z,
    gdT_vs_abT = gdT_core_z - abT_core_z,
    gdT_net_legacy = if ("gdT_core_legacy_z" %in% names(.))
      gdT_core_legacy_z - inhibitory_program_z else NA_real_
  )

write_tsv(spot, file.path(out_tables, "09_spot_scores.tsv"))
saveRDS(spot,   file.path(out_rds, "09_spot_scores.rds"))

# =========================================================
# SECTION SUMMARIES
# =========================================================
summary_cols <- c(score_cols, paste0(score_cols, "_z"),
                  "abT_net", "gdT_net", "gdT_vs_abT", "gdT_net_legacy")
summary_cols <- intersect(summary_cols, colnames(spot))

section_by_zone <- spot %>%
  group_by(section_id, zone_call) %>%
  summarise(n_spots = n(),
            across(all_of(summary_cols), list(mean = mean_or_na, median = median_or_na),
                   .names = "{.col}_{.fn}"),
            .groups = "drop")

mes_fraction <- zones_all %>%
  filter(in_tissue == 1, zone_call %in% zones_keep) %>%
  count(section_id, zone_call, name = "n_spots") %>%
  group_by(section_id) %>% mutate(frac = n_spots / sum(n_spots)) %>%
  filter(zone_call == "Mesenchymal") %>% ungroup() %>%
  transmute(section_id, mesenchymal_fraction = frac)

section_summary <- spot %>%
  group_by(section_id) %>%
  summarise(n_spots = n(),
            across(all_of(summary_cols), list(mean = mean_or_na, median = median_or_na),
                   .names = "{.col}_{.fn}"),
            .groups = "drop") %>%
  left_join(manual_section_metadata, by = "section_id") %>%
  left_join(mes_fraction, by = "section_id")

section_summary <- section_summary %>%
  mutate(
    opportunity_category = classify_opportunity(abT_net_median, gdT_net_median,
                                                inhibitory_program_z_median),
    opportunity_label = display_category(opportunity_category),
    section_label = ifelse(is.na(sample_status) | sample_status == "",
                           section_id, paste0(section_id, " (", sample_status, ")"))
  )

write_tsv(section_by_zone, file.path(out_tables, "09_section_by_zone_summary.tsv"))
write_tsv(section_summary, file.path(out_tables, "09_section_summary.tsv"))
saveRDS(section_summary,   file.path(out_rds, "09_section_summary.rds"))

message("\n--- Section classification ---")
print(section_summary %>%
        select(section_id, sample_status, abT_net_median, gdT_net_median,
               inhibitory_program_z_median, opportunity_category), n = Inf)

high_inhib <- section_summary %>%
  filter(is.finite(inhibitory_program_z_median),
         inhibitory_program_z_median >= SUPPRESSED_INHIB_Z)
message("\nElevated inhibitory tone (z >= ", SUPPRESSED_INHIB_Z, "): ",
        if (nrow(high_inhib) == 0) "none" else paste(high_inhib$section_id, collapse = ", "))

# =========================================================
# CORE REVISION COMPARISON (mevalonate added vs not)
# =========================================================
if ("gdT_net_legacy_median" %in% colnames(section_summary)) {
  rev_cmp <- section_summary %>%
    mutate(category_legacy = classify_opportunity(abT_net_median, gdT_net_legacy_median,
                                                  inhibitory_program_z_median)) %>%
    select(section_id, gdT_net_median, gdT_net_legacy_median,
           opportunity_category, category_legacy) %>%
    mutate(changed = as.character(opportunity_category) != as.character(category_legacy))

  write_tsv(rev_cmp, file.path(out_tables, "09_core_revision_comparison.tsv"))

  rho_rev <- suppressWarnings(cor(rev_cmp$gdT_net_median, rev_cmp$gdT_net_legacy_median,
                                  method = "spearman", use = "complete.obs"))
  message("\n--- gdT core revision (mevalonate added) ---")
  message("Section-level Spearman vs legacy core: ", sprintf("%.3f", rho_rev))
  message("Sections changing category: ", sum(rev_cmp$changed, na.rm = TRUE), "/", nrow(rev_cmp))
  if (any(rev_cmp$changed, na.rm = TRUE)) {
    print(rev_cmp %>% filter(changed) %>%
            select(section_id, category_legacy, opportunity_category), n = Inf)
  }
}

# =========================================================
# MESENCHYMAL BURDEN CORRELATIONS
# =========================================================
corr_tbl <- bind_rows(
  run_spearman_safe(section_summary, "mesenchymal_fraction", "abT_net_median"),
  run_spearman_safe(section_summary, "mesenchymal_fraction", "gdT_net_median"),
  run_spearman_safe(section_summary, "mesenchymal_fraction", "inhibitory_program_z_median")
) %>% mutate(p_adj = p.adjust(p_value, method = "BH"))
write_tsv(corr_tbl, file.path(out_tables, "09_mesenchymal_correlations.tsv"))

# =========================================================
# FIGURE 5A — opportunity quadrant
# =========================================================
p_quad <- ggplot(section_summary,
                 aes(x = abT_net_median, y = gdT_net_median,
                     label = section_id, color = opportunity_label)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey70") +
  geom_vline(xintercept = 0, linetype = 2, color = "grey70") +
  geom_point(size = 3.5) +
  ggrepel::geom_text_repel(size = 3, max.overlaps = Inf, min.segment.length = 0) +
  scale_color_manual(values = setNames(unname(OPPORTUNITY_COLORS[OPPORTUNITY_LEVELS]),
                                       unname(OPPORTUNITY_LABELS[OPPORTUNITY_LEVELS])),
                     drop = FALSE) +
  theme_bw(base_size = 11) +
  labs(x = AXIS_LAB_AB, y = AXIS_LAB_GD, color = "Opportunity\ncategory",
       title = paste0("Immunotherapy opportunity across ependymoma sections")) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(out_plots, "Fig5A_opportunity_quadrant.png"), p_quad,
       width = 8.8, height = 6.2, dpi = 600, bg = "white")
ggsave(file.path(out_plots, "Fig5A_opportunity_quadrant.pdf"), p_quad,
       width = 8.8, height = 6.2, bg = "white", device = cairo_pdf)

# =========================================================
# FIGURE 5B — gdT mechanism axes
#
# Bar panels with a free x scale per axis: one outlier section cannot flatten
# the contrast in the others, and the axes are independent by design so a
# shared scale would imply a comparability they do not have.
# =========================================================
make_axis_panel <- function(axis_names, title_txt, subtitle_txt, ncol = 3) {
  use <- axis_names[sapply(axis_names, function(a) length(gene_sets[[a]]) > 0)]
  dropped <- setdiff(axis_names, use)
  if (length(dropped) > 0) {
    message("Axes with no gene above the detection floor (not plotted): ",
            paste(dropped, collapse = ", "))
  }
  if (length(use) == 0) return(NULL)

  cols <- paste0(use, "_mean")
  cols <- intersect(cols, colnames(section_summary))
  if (length(cols) == 0) return(NULL)

  ord <- section_summary %>% arrange(desc(gdT_net_median)) %>% pull(section_label)

  df <- section_summary %>%
    select(section_label, all_of(cols)) %>%
    pivot_longer(all_of(cols), names_to = "axis", values_to = "value") %>%
    mutate(axis = sub("_mean$", "", axis),
           axis_lab = factor(axis_label(axis), levels = axis_label(use)),
           section_label = factor(section_label, levels = rev(ord))) %>%
    filter(is.finite(value))

  ggplot(df, aes(x = section_label, y = value, fill = value > 0)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey70") +
    geom_col(width = 0.7) +
    coord_flip() +
    facet_wrap(~axis_lab, scales = "free_x", ncol = ncol) +
    scale_fill_manual(values = c(`TRUE` = "#00BF7D", `FALSE` = "#9ecae1"),
                      guide = "none") +
    labs(x = NULL, y = "Section mean score",
         title = title_txt, subtitle = subtitle_txt) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold"),
          panel.grid.major.y = element_blank())
}

p_gd_axes <- make_axis_panel(
  c(PROGRAM_AXES$gdT_core, names(MODIFIERS)),
  paste0(AXIS_LAB_GD, " recognition modalities by section"),
  "Modalities are independent by design; the composite weights each equally"
)

if (!is.null(p_gd_axes)) {
  h <- max(6, 0.22 * nrow(section_summary) * 2)
  ggsave(file.path(out_plots, "Fig5B_gdT_mechanism_axes.png"), p_gd_axes,
         width = 11, height = h, dpi = 600, bg = "white")
  ggsave(file.path(out_plots, "Fig5B_gdT_mechanism_axes.pdf"), p_gd_axes,
         width = 11, height = h, bg = "white", device = cairo_pdf)
}

p_inh_axes <- make_axis_panel(
  PROGRAM_AXES$inhibitory_program,
  "Immunosuppressive mechanisms by section",
  "Mechanisms are independent by design; the composite weights each equally"
)

if (!is.null(p_inh_axes)) {
  h <- max(6, 0.22 * nrow(section_summary) * 3)
  ggsave(file.path(out_plots, "S_inhibitory_mechanism_axes.png"), p_inh_axes,
         width = 11, height = h, dpi = 600, bg = "white")
  ggsave(file.path(out_plots, "S_inhibitory_mechanism_axes.pdf"), p_inh_axes,
         width = 11, height = h, bg = "white", device = cairo_pdf)
}

# =========================================================
# SUPPLEMENTARY — inhibitory tone, matched pairs
# =========================================================
p_tone <- section_summary %>%
  mutate(elevated = inhibitory_program_z_median >= SUPPRESSED_INHIB_Z) %>%
  ggplot(aes(x = reorder(section_label, inhibitory_program_z_median),
             y = inhibitory_program_z_median, fill = elevated)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey70") +
  geom_hline(yintercept = SUPPRESSED_INHIB_Z, linetype = 3, colour = "grey30") +
  geom_col(width = 0.7) + coord_flip() +
  scale_fill_manual(values = c(`FALSE` = "#9ecae1", `TRUE` = "#E76BF3"),
                    labels = c("No", "Yes"), name = "Elevated") +
  labs(x = NULL, y = "Median immunosuppressive tone (z)",
       title = "Immunosuppressive tone across sections",
       subtitle = paste0("Dotted line: elevated-tone threshold (z = ", SUPPRESSED_INHIB_Z, ")")) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), panel.grid.major.y = element_blank())

ggsave(file.path(out_plots, "S_inhibitory_tone_by_section.png"), p_tone,
       width = 7.5, height = max(4.5, 0.35 * nrow(section_summary)), dpi = 600, bg = "white")
ggsave(file.path(out_plots, "S_inhibitory_tone_by_section.pdf"), p_tone,
       width = 7.5, height = max(4.5, 0.35 * nrow(section_summary)),
       bg = "white", device = cairo_pdf)

paired <- section_summary %>%
  filter(!is.na(sample_status), !is.na(patient_id)) %>%
  select(patient_id, sample_status, section_id,
         abT_net_median, gdT_net_median, inhibitory_program_z_median,
         opportunity_category)

if (nrow(paired) > 0) {
  write_tsv(paired, file.path(out_tables, "09_primary_relapse_pairs.tsv"))

  paired_long <- paired %>%
    select(patient_id, sample_status, abT_net_median, gdT_net_median) %>%
    pivot_longer(c(abT_net_median, gdT_net_median),
                 names_to = "score", values_to = "value") %>%
    mutate(score = recode(score, abT_net_median = AXIS_LAB_AB,
                          gdT_net_median = AXIS_LAB_GD),
           sample_status = factor(sample_status, levels = c("Primary", "Relapse")))

  p_paired <- ggplot(paired_long,
                     aes(x = sample_status, y = value, group = patient_id, colour = patient_id)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey70") +
    geom_line(linewidth = 0.7) + geom_point(size = 3) +
    facet_wrap(~score) +
    labs(x = NULL, y = "Section median net score", colour = "Patient",
         title = "Matched primary versus relapse sections",
         subtitle = "Within-patient comparison; each line is one patient") +
    theme_bw(base_size = 11) + theme(plot.title = element_text(face = "bold"))

  ggsave(file.path(out_plots, "S_primary_relapse_paired.png"), p_paired,
         width = 8.5, height = 4.8, dpi = 600, bg = "white")
  ggsave(file.path(out_plots, "S_primary_relapse_paired.pdf"), p_paired,
         width = 8.5, height = 4.8, bg = "white", device = cairo_pdf)

  message("\n--- Matched primary / relapse pairs ---")
  print(paired %>% arrange(patient_id, sample_status), n = Inf)
}

# =========================================================
# RUN PROVENANCE
# =========================================================
writeLines(c(
  paste0("scoring_method: ", SCORING_METHOD),
  paste0("scoring_aggregation: ", SCORING_AGGREGATION),
  paste0("detection_min_pct: ", DETECTION_MIN_PCT),
  paste0("suppressed_inhibitory_threshold_z: ", SUPPRESSED_INHIB_Z),
  "z_scope: global across pooled tumor spots",
  paste0("zones_scored: ", paste(zones_keep, collapse = ",")),
  paste0("n_spots: ", nrow(spot)),
  paste0("n_sections: ", n_distinct(spot$section_id)),
  paste0("run_date: ", as.character(Sys.Date()))
), file.path(out_tables, "09_run_parameters.txt"))

message("\nDone. Outputs written to: ", out_root)
