# =============================================================================
# config.R — single source of truth for paths, constants and seeds.
#
# Every script sources this file and NEVER hardcodes a path or a threshold.
#
#   source(here::here("config.R"))
#
# The STAGE list below maps a STABLE KEY to the folder name CURRENTLY on disk.
# Adopt this file first without moving anything; rename folders later by
# editing the values here only. Scripts never see folder names.
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Roots
# -----------------------------------------------------------------------------
# Set in .Renviron, e.g.
#   EPN_DATA_ROOT=D:/Ped-CNS_KBH/Spatial Transcriptomic/GSE195661
DATA_ROOT <- Sys.getenv("EPN_DATA_ROOT", unset = file.path(getwd(), "data"))

ANALYSIS_ROOT <- file.path(DATA_ROOT, "Analysis")
RAW_ROOT      <- file.path(DATA_ROOT, "neuro_onc_spatial_files")

FIGURE_ROOT <- Sys.getenv("EPN_FIGURE_ROOT",
                          unset = file.path(getwd(), "outputs", "figures"))

if (!dir.exists(DATA_ROOT)) {
  stop("DATA_ROOT does not exist: ", DATA_ROOT,
       "\nSet EPN_DATA_ROOT in .Renviron, or run scripts/00_build_manifest.R.")
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

# -----------------------------------------------------------------------------
# 2. External reference files (live at DATA_ROOT, not under Analysis/)
#
# Provenance for each is in docs/data_sources.md. None are redistributed.
# -----------------------------------------------------------------------------
REF <- list(
  # GSE125969 single-cell reference (via pneuroonccellatlas.org repackaging)
  scrna_expr        = file.path(DATA_ROOT, "exprMatrix.tsv"),
  scrna_meta        = file.path(DATA_ROOT, "meta.tsv"),

  # Neoplastic-population subset used to build the augmented reference (03a)
  tumor_expr        = file.path(DATA_ROOT, "exprMatrix_tumor.tsv"),
  tumor_meta        = file.path(DATA_ROOT, "meta_tumor.tsv"),

  # Zone gene programs — Donson et al. Neuro-Oncology supplementary data 1
  zone_genes_xlsx   = file.path(DATA_ROOT, "noac219_suppl_supplementary_data1.xlsx"),

  # Built by 03a
  augmented_ref_dir = file.path(DATA_ROOT, "tumor_augmented_reference_PFA")
)

ref_file <- function(key) {
  p <- REF[[key]]
  if (is.null(p)) stop("Unknown reference key: '", key, "'")
  if (!file.exists(p) && !dir.exists(p)) {
    stop("Missing reference file: ", p,
         "\nSee docs/data_sources.md for where to obtain it.")
  }
  p
}

# -----------------------------------------------------------------------------
# 3. Stage directories (values = folder names as they exist TODAY)
#
# Where the current folder name disagrees with the script number, the mismatch
# is noted. Fix by editing the value here and moving the folder once.
# -----------------------------------------------------------------------------
STAGE <- list(
  manifest           = "00_manifest",
  raw_vis            = "00_raw_sections/vis",
  raw_spot_tables    = "00_raw_sections/spot_tables",
  raw_qc_png         = "00_raw_sections/qc_png",

  coord_maps         = "01_barcode_coordinate_maps",

  semla_nnls         = "02_semla_nnls",
  hotspots           = "02_semla_nnls/good_sections_hotspots_tissue",

  # written by 03b — folder is numbered 02-4
  semla_nnls_aug     = "02-4_semla_nnls_tumor_augmented_PFA",
  # written by 03c — folder is unnumbered
  hotspot_robustness = "Hotspot_tumor_augmented_PFA_robustness",

  # written by 04 — folder is numbered 03c, which COLLIDES with script 03c
  he_overlays        = "03c_he_overlays_from_png_tissue_v2",

  # NOTE: there is no stage 04. The former 04_lm7_deconvolution (per-spot
  # UCell) was superseded by the CIBERSORTx route and its script no longer
  # exists; stage 05 now reads the 02c hotspot objects directly.
  pseudobulk         = "04_hotspot_pseudobulk_lm7",
  deeptil            = "04_hotspot_pseudobulk_lm7/deeptil/results",

  # 06 is run TWICE with different section_keep, producing two zone tables:
  #   zones      -- GOOD_SECTIONS only; consumed by 07 and 08
  #   zones_all  -- all 14 sections;    consumed by 09 and 10
  # Only the first run is scripted. See docs/pipeline.md.
  zones              = "06_zones",
  zones_all          = "06_zones_all_sections",
  # written by 06b — folder is unnumbered
  signature_overlap  = "qc_signature_overlap",

  functional_states  = "07_functional_states",
  lr_interactions    = "08_lr_interactions",
  lr_qc              = "08_lr_interactions/qc_cross_section",

  scoring            = "09_opportunity_map",
  scoring_qc         = "09b_scoring_qc",
  zone_spatial       = "10_zone_and_spatial",
  gsea               = "11_gsea_characterization"
)

stage_dir <- function(key, ...) {
  if (!key %in% names(STAGE)) {
    stop("Unknown stage key: '", key, "'. Known: ",
         paste(names(STAGE), collapse = ", "))
  }
  d <- ensure_dir(file.path(ANALYSIS_ROOT, STAGE[[key]]))
  if (length(list(...)) == 0) d else file.path(d, ...)
}

# Read-only accessor: fails loudly if an upstream stage has not been run.
stage_input <- function(key, file) {
  p <- file.path(ANALYSIS_ROOT, STAGE[[key]], file)
  if (!file.exists(p)) {
    stop("Missing input: ", p,
         "\nProduced by an earlier stage. See docs/pipeline.md for run order.")
  }
  p
}

figure_path <- function(...) file.path(ensure_dir(FIGURE_ROOT), ...)

# -----------------------------------------------------------------------------
# 4. Cohort constants
# -----------------------------------------------------------------------------
ALL_SECTIONS <- c("459", "459_2", "723", "723_2", "727", "812", "821",
                  "848", "928", "928_2", "1101", "1239", "1269", "1513")

GOOD_SECTIONS <- c("459", "723", "812", "821", "928", "1239")

MATCHED_PAIRS <- list(
  c(primary = "459", relapse = "459_2"),
  c(primary = "723", relapse = "723_2"),
  c(primary = "928", relapse = "928_2")
)

# -----------------------------------------------------------------------------
# 5. Analysis thresholds
# -----------------------------------------------------------------------------
THRESH <- list(
  # 02 automatic section screen — all LOWER bounds; see docs/pipeline.md
  min_nonNA_spots      = 50L,
  min_prop_nonzero     = 0.10,
  min_max_lymphocyte   = 0.05,

  # 02c hotspot definition: within-section quantiles
  hotspot_quantile     = 0.90,
  hotspot_knn          = 6L,

  # 03a augmented reference construction
  cells_per_cluster    = 250L,

  # 05 pseudobulk pairing: minimum spots in BOTH groups per section
  min_spots_pseudobulk = 20L,

  # 09 opportunity scoring gene detection floor (fraction of tumor spots)
  detection_floor      = 0.05,

  # 11 GSEA — separate, stricter floor
  gsea_detection_floor = 0.10,
  gsea_min_sections    = 10L,
  gsea_min_set_size    = 10L,
  gsea_max_set_size    = 500L
)

# -----------------------------------------------------------------------------
# 6. External tool parameters (runs happen outside R; recorded for Methods)
# -----------------------------------------------------------------------------
CIBERSORTX <- list(
  signature_matrix   = "LM7 (Tosolini et al. 2017, OncoImmunology 6:e1284723)",
  batch_correction   = "B-mode",
  quantile_normalise = FALSE,
  mode               = "relative",
  permutations       = 500L
)

DEEPTIL <- list(
  input  = "pseudobulk_cpm_CIBERSORTx.tsv via CIBERSORTx fractions",
  output = "SES_CIBERSORTx_EPN_pseudobulk.txt",
  version = NA_character_   # TODO: record tool version used
)

# -----------------------------------------------------------------------------
# 7. Seeds
#
# 03a, 09b and 11 were already seeded and their published values are preserved
# here verbatim. 08 was NOT seeded and its value is new, so stage 08 must be
# re-run. Do not change any value once results are published.
# -----------------------------------------------------------------------------
SEEDS <- list(
  tumor_augmented_reference = 1L,          # 03a -- as published, do not change
  lr_permutation            = 20260801L,   # 08  -- newly added; re-run required
  split_half                = 20260812L,   # 09b -- as published, do not change
  gsea                      = 20260812L    # 11  -- as published, do not change
)
# 10_zone_and_spatial.R contains no stochastic step and needs no seed.

use_seed <- function(key) {
  if (!key %in% names(SEEDS)) stop("Unknown seed key: '", key, "'")
  set.seed(SEEDS[[key]])
  message("Seed set: ", key, " = ", SEEDS[[key]])
  invisible(SEEDS[[key]])
}

# -----------------------------------------------------------------------------
# 8. Plot conventions
# -----------------------------------------------------------------------------
PLOT <- list(
  dpi           = 600,
  spatial_scale = "plasma",
  device_vector = "cairo_pdf"
)

invisible(TRUE)
