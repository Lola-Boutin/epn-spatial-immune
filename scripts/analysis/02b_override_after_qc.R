# =============================================================================
# 02b_override_after_qc.R
#
# Purpose
# -------
# Apply the manual section-selection decision on top of the automatic screen
# from stage 02, then rebuild good_sections/ and the merged vis_good object.
#
# WHY A MANUAL STEP IS NEEDED
# The three criteria in stage 02 are all LOWER bounds, so a section whose
# Lymphocyte score is uniformly inflated across the whole capture area passes
# by construction. Section 848 did exactly that: it cleared the automatic
# screen but its per-spot map showed near-uniform elevation with no focal
# structure, against the sparse focal pattern of the six retained sections
# (compare Lymphocytes_global_848.png with Lymphocytes_global_459.png).
# It is quantifiable on the screen's own metrics: prop_nonzero 0.978 versus
# 0.144-0.484 for the retained six, and mean_ly 0.298 versus 0.016-0.076.
#
# The control deconvolution against the full reference came AFTERWARDS and
# confirmed the call. It is independent confirmation, not the selection
# criterion, and Methods must state the order that way.
#
# Inputs
# ------
#   stage_dir("semla_nnls")   all_sections/{vis,nnls}
#
# Outputs
# -------
#   stage_dir("semla_nnls")   good_sections/{vis,nnls} rebuilt,
#                             qc/good_sections_manual_override.txt,
#                             merged/vis_good_semla_ready_manual_override.rds
#
# Stochastic
# ----------
#   none
#
# Runtime
# -------
#   ~10 minutes
# =============================================================================

source(here::here("config.R"))

library(Seurat)
library(readr)
library(dplyr)
library(purrr)
library(tibble)


# -------------------------
# Configuration
# -------------------------
out_root      <- stage_dir("semla_nnls")
out_all_vis   <- ensure_dir(file.path(out_root, "all_sections", "vis"))
out_all_nnls  <- ensure_dir(file.path(out_root, "all_sections", "nnls"))
out_good_vis  <- ensure_dir(file.path(out_root, "good_sections", "vis"))
out_good_nnls <- ensure_dir(file.path(out_root, "good_sections", "nnls"))
out_qc        <- ensure_dir(file.path(out_root, "qc"))
out_merged    <- ensure_dir(file.path(out_root, "merged"))

# -------------------------
# Manual override
# -------------------------
good_sections_override <- GOOD_SECTIONS

# Optional: save the decision
write_lines(good_sections_override, file.path(out_qc, "good_sections_manual_override.txt"))

# -------------------------
# Clean old good_sections outputs
# -------------------------
old_good_vis  <- list.files(out_good_vis,  full.names = TRUE)
old_good_nnls <- list.files(out_good_nnls, full.names = TRUE)

if (length(old_good_vis) > 0) file.remove(old_good_vis)
if (length(old_good_nnls) > 0) file.remove(old_good_nnls)

# -------------------------
# Rebuild good_sections folders
# -------------------------
for (sec in good_sections_override) {
  vis_file <- file.path(out_all_vis,  paste0("vis_section_", sec, "_semla.rds"))
  nnls_rds <- file.path(out_all_nnls, paste0("NNLS_section_", sec, ".rds"))
  nnls_tsv <- file.path(out_all_nnls, paste0("NNLS_section_", sec, ".tsv"))
  
  if (!file.exists(vis_file)) stop("Missing vis file: ", vis_file)
  if (!file.exists(nnls_rds)) stop("Missing NNLS rds: ", nnls_rds)
  
  file.copy(
    from = vis_file,
    to   = file.path(out_good_vis, paste0("vis_section_", sec, ".rds")),
    overwrite = TRUE
  )
  
  file.copy(
    from = nnls_rds,
    to   = file.path(out_good_nnls, paste0("NNLS_section_", sec, ".rds")),
    overwrite = TRUE
  )
  
  if (file.exists(nnls_tsv)) {
    file.copy(
      from = nnls_tsv,
      to   = file.path(out_good_nnls, paste0("NNLS_section_", sec, ".tsv")),
      overwrite = TRUE
    )
  }
}

# -------------------------
# Rebuild vis_good from manual set
# -------------------------
vis_list <- lapply(good_sections_override, function(sec) {
  readRDS(file.path(out_good_vis, paste0("vis_section_", sec, ".rds")))
})
names(vis_list) <- good_sections_override

vis_good <- vis_list[[1]]
if (length(vis_list) > 1) {
  for (sec in names(vis_list)[-1]) {
    vis_good <- merge(vis_good, vis_list[[sec]], merge.data = TRUE)
  }
}

saveRDS(vis_good, file.path(out_merged, "vis_good_semla_ready_manual_override.rds"))
write_tsv(
  vis_good@meta.data %>% rownames_to_column("cell"),
  file.path(out_merged, "vis_good_semla_ready_manual_override_metadata.tsv")
)

# -------------------------
# Rebuild merged flat table
# -------------------------
df_all_good <- map_dfr(
  good_sections_override,
  ~ readRDS(file.path(out_good_nnls, paste0("NNLS_section_", .x, ".rds")))
)

write_tsv(
  df_all_good,
  file.path(out_merged, "df_all_good_sections_manual_override.tsv")
)

message("Manual override complete.")