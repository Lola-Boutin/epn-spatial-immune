# =============================================================================
# 11_gsea_characterization.R
#
# Purpose
# -------
# Pathway characterisation of the opportunity scores.
#
# For each target score, genes are ranked by their cross-section correlation
# with that score, then fgsea is run on the combined z ranking with redundant
# pathways collapsed. fgsea is also run per section, and the number of sections
# in which a pathway is independently significant is reported.
#
# Collapsing matters: without it a single signal can appear as many separate
# "findings" because member pathways overlap heavily.
#
# Compositional-artifact genes (ribosomal, mitochondrial, haemoglobin) are
# excluded from the ranking but still diagnosed and reported.
#
# Inputs
# ------
#   stage_input("manifest", "visium_manifest.tsv")
#   stage_input("scoring", "rds/09_spot_scores.rds")
#   stage_dir("raw_vis")   vis_section_<sec>_raw.rds
#   MSigDB Hallmark and C2:CP:REACTOME via msigdbr (downloaded at runtime;
#     not redistributed -- see docs/data_sources.md)
#
# Outputs
# -------
#   stage_dir("gsea")
#     tables/11_gene_ranking_<target>.tsv
#     tables/11_gsea_full_<target>.tsv, 11_gsea_main_<target>.tsv,
#     tables/11_gsea_discovery_<target>.tsv
#     tables/11_pathway_program_overlap.tsv
#     tables/11_artifact_diagnostic_<target>.tsv
#     rds/11_per_section_correlations.rds   (cache; see reuse_cached)
#
# Stochastic
# ----------
#   YES -- use_seed("gsea"), 20260812 as published.
#   fgsea::fgseaMultilevel uses randomised multilevel sampling.
#
# Runtime
# -------
#   ~1 hour
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
  library(patchwork)
  library(fgsea)
})

# =========================================================
# 11 — GSEA characterisation of the opportunity scores
#
# DESIGN
# ------
# Genes are ranked by their association with the opportunity score WITHIN each
# section, and the 14 per-section rankings are then combined. Spots provide the
# resolution to estimate each section's effect; sections provide the
# replication for inference. Nothing is compared across patients, so patient
# identity cannot drive the result - which is the flaw in a
# one-section-versus-one-section contrast.
#
# Per section: Spearman rho between each gene and the spot-level score.
# Across sections: Fisher z-transform, weighted Stouffer combination
# (weight = sqrt(n_spots - 3)), giving one combined z per gene.
# Enrichment: fgsea on the combined z ranking, with redundant pathways
# collapsed. fgsea is also run per section, and the number of sections in which
# each pathway's NES has the same sign as the combined result is reported. At
# n = 14, "concordant in 13/14" is stronger evidence than any single p-value.
#
# TWO ARTIFACTS THIS SCRIPT CONTROLS FOR
# --------------------------------------
# 1. COMPOSITIONAL. Ribosomal protein transcripts are a large share of counts
#    per spot, so when any other program rises their normalised fraction falls.
#    Library-size normalisation turns that into an apparent negative
#    correlation with every score, and because RP genes move as a coordinated
#    family of several hundred it produces extreme NES values. RP, MT and
#    haemoglobin genes are therefore excluded from the ranking by default. The
#    excluded genes are still scored and reported separately (see the
#    artifact-diagnostic table) so the effect is documented rather than hidden.
#
# 2. REDUNDANCY. The Reactome hierarchy contains many near-duplicate sets, so a
#    single signal can appear as seven "findings". fgsea::collapsePathways is
#    applied and only main pathways are plotted; the full table is still
#    written.
#
# CONSISTENCY VERSUS DISCOVERY
# ----------------------------
# Program genes are NOT removed from the ranking. Pathways containing them are
# expected to enrich by construction, so each pathway is annotated with the
# number of program genes it contains and the fraction of its leading edge they
# represent. This separates CONSISTENCY (the score behaves as designed) from
# DISCOVERY (associations independent of the scoring genes).
#
# CAVEAT
# ------
# The score is derived from expression, so anything co-regulated with the
# program genes will correlate with it. These are associations, not causes; the
# cross-section concordance is what makes them findings rather than
# descriptions.
#
# OUTPUTS
# -------
#   tables/11_gene_ranking_<target>.tsv
#   tables/11_gsea_full_<target>.tsv            all pathways
#   tables/11_gsea_main_<target>.tsv            after collapsing
#   tables/11_gsea_discovery_<target>.tsv       no program-gene overlap
#   tables/11_gsea_section_concordance_<target>.tsv
#   tables/11_artifact_diagnostic_<target>.tsv  RP/MT/HB association
#   tables/11_pathway_program_overlap.tsv
#   plots/11_gsea_main_<target>.png | .pdf
#   plots/11_gsea_discovery_<target>.png | .pdf
#   plots/11_artifact_diagnostic_<target>.png | .pdf
# =========================================================


# -------------------------
# CONFIG
# -------------------------
manifest_file <- stage_input("manifest", "visium_manifest.tsv")
raw_vis_dir   <- stage_dir("raw_vis")
phase09_rds   <- file.path(stage_dir("scoring"), "rds")

out_root   <- stage_dir("gsea")
out_tables <- ensure_dir(file.path(out_root, "tables"))
out_plots  <- ensure_dir(file.path(out_root, "plots"))
out_rds    <- ensure_dir(file.path(out_root, "rds"))

# Scores to characterise. Each must be a column in 09_spot_scores.rds.
TARGETS <- c("gdT_core", "abT_core", "inhibitory_program",
             "gdT_btn", "gdT_mevalonate", "gdT_nkg2d",
             "gdT_ephrin", "gdT_adhesion_dnam")

# Gene must be detected in at least this fraction of a section's tumor spots
GENE_MIN_DETECT <- THRESH$gsea_detection_floor

# A gene must be testable in at least this many sections to enter the ranking
MIN_SECTIONS <- THRESH$gsea_min_sections

# Compositional-artifact genes. Excluded from the ranking; still diagnosed.
EXCLUDE_ARTIFACT_GENES <- TRUE
ARTIFACT_PATTERN <- "^(RPL|RPS|MRPL|MRPS|MT-|MTRNR|HBA|HBB|HBD|HBG|HBM|HBQ|HBZ)"

# Gene set collections
USE_HALLMARK <- TRUE
USE_REACTOME <- TRUE
GS_MIN <- THRESH$gsea_min_set_size
GS_MAX <- THRESH$gsea_max_set_size

COLLAPSE_PATHWAYS <- TRUE
PADJ_CUT   <- 0.05
N_TOP_PLOT <- 20
SEED <- SEEDS$gsea

reuse_cached <- FALSE   # TRUE only if the ranking config is unchanged

use_seed("gsea")

# =========================================================
# PROGRAM GENES (annotation only - NOT removed from the ranking)
# Must match 09_opportunity_map.R; both alias variants listed.
# =========================================================
PROGRAM_GENES <- toupper(unique(c(
  "HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "PSMB8", "PSMB9",
  "ERAP1", "ERAP2", "NLRC5", "IRF1",
  "BTN2A1", "BTN3A1", "MVK", "PMVK", "MVD", "IDI1",
  "MICA", "MICB", "ULBP2", "ULBP3",
  "ULBP4", "RAET1E", "ULBP5", "RAET1G", "ULBP6", "RAET1L",
  "EPHA2", "PVR", "NECTIN2", "ICAM1",
  "CD274", "PDCD1LG2", "HLA-E", "LGALS9", "IDO1",
  "TGFB1", "TGFB2", "TGFB3", "IL10", "VEGFA",
  "NT5E", "ENTPD1", "PTGS2", "PTGES", "PTGER2", "PTGER4",
  "FDPS"
)))

# =========================================================
# HELPERS
# =========================================================
read_expr_matrix <- function(obj) {
  DefaultAssay(obj) <- "Spatial"
  mat <- tryCatch(GetAssayData(obj, assay = "Spatial", layer = "data"),
                  error = function(e) NULL)
  if (is.null(mat) || nrow(mat) == 0 || ncol(mat) == 0) {
    message("  normalising from counts")
    obj <- NormalizeData(obj, normalization.method = "LogNormalize",
                         scale.factor = 1e4, verbose = FALSE)
    mat <- GetAssayData(obj, assay = "Spatial", layer = "data")
  }
  if (is.null(mat) || nrow(mat) == 0) stop("Could not obtain normalized matrix")
  rn <- toupper(rownames(mat)); keep <- !duplicated(rn)
  mat <- mat[keep, , drop = FALSE]; rownames(mat) <- rn[keep]
  mat
}

# Weighted Stouffer combination of Fisher-z transformed correlations
combine_rho <- function(rho, n) {
  ok <- is.finite(rho) & is.finite(n) & n > 3
  if (!any(ok)) return(NA_real_)
  z <- atanh(pmax(pmin(rho[ok], 0.999999), -0.999999))
  w <- sqrt(n[ok] - 3)
  sum(w * z) / sqrt(sum(w^2))
}

load_gene_sets <- function() {
  if (!requireNamespace("msigdbr", quietly = TRUE))
    stop("msigdbr not installed. Run: install.packages('msigdbr')")
  
  # msigdbr >= 10 renamed 'category' to 'collection'
  get_sets <- function(...) {
    args <- list(species = "Homo sapiens", ...)
    tryCatch(do.call(msigdbr::msigdbr, args),
             error = function(e) {
               names(args)[names(args) == "collection"]    <- "category"
               names(args)[names(args) == "subcollection"] <- "subcategory"
               do.call(msigdbr::msigdbr, args)
             })
  }
  
  out <- list()
  if (USE_HALLMARK) {
    h <- get_sets(collection = "H")
    out <- c(out, split(toupper(h$gene_symbol), h$gs_name))
    message("Hallmark sets: ", length(unique(h$gs_name)))
  }
  if (USE_REACTOME) {
    r <- get_sets(collection = "C2", subcollection = "CP:REACTOME")
    out <- c(out, split(toupper(r$gene_symbol), r$gs_name))
    message("Reactome sets: ", length(unique(r$gs_name)))
  }
  lapply(out, unique)
}

pretty_path <- function(x, width = 55) {
  substr(gsub("_", " ", sub("^(HALLMARK|REACTOME)_", "", x)), 1, width)
}

# =========================================================
# PER-SECTION CORRELATIONS
#
# Artifact genes are correlated too, so their behaviour can be reported, but
# they are flagged and dropped before ranking.
# =========================================================
cache_file <- file.path(out_rds, "11_per_section_correlations.rds")

if (reuse_cached && file.exists(cache_file)) {
  message("Reusing cached per-section correlations")
  cor_long <- readRDS(cache_file)
} else {
  
  spot <- readRDS(file.path(phase09_rds, "09_spot_scores.rds"))
  manifest <- read_tsv(manifest_file, show_col_types = FALSE) %>%
    filter(complete %in% TRUE) %>% arrange(section_id)
  
  missing_targets <- setdiff(TARGETS, colnames(spot))
  if (length(missing_targets) > 0)
    stop("Target(s) not in 09 output: ", paste(missing_targets, collapse = ", "))
  
  cor_rows <- list()
  
  for (sec in intersect(manifest$section_id, unique(spot$section_id))) {
    message("\n=== Section ", sec, " ===")
    
    vis_file <- file.path(raw_vis_dir, paste0("vis_section_", sec, "_raw.rds"))
    if (!file.exists(vis_file)) { warning("Missing vis file: ", sec); next }
    
    sp  <- spot %>% filter(section_id == sec)
    obj <- readRDS(vis_file)
    m   <- read_expr_matrix(obj)
    
    cells <- intersect(sp$cell, colnames(m))
    if (length(cells) < 50) {
      warning("Too few matched spots in section ", sec, " (", length(cells), ")"); next
    }
    sp <- sp %>% filter(cell %in% cells) %>% slice(match(cells, cell))
    m  <- m[, cells, drop = FALSE]
    
    m <- m[Matrix::rowMeans(m > 0) >= GENE_MIN_DETECT, , drop = FALSE]
    message("  ", nrow(m), " genes x ", ncol(m), " spots")
    
    mt      <- as.matrix(Matrix::t(m))
    rank_mt <- apply(mt, 2, rank)
    
    for (tgt in TARGETS) {
      y  <- sp[[tgt]]
      ok <- is.finite(y)
      if (sum(ok) < 50) next
      r <- suppressWarnings(as.numeric(cor(rank_mt[ok, , drop = FALSE], rank(y[ok]))))
      cor_rows[[paste(sec, tgt, sep = "__")]] <- tibble(
        section_id = sec, target = tgt, gene = colnames(mt),
        rho = r, n_spots = sum(ok))
    }
    
    rm(obj, m, mt, rank_mt); gc(verbose = FALSE)
  }
  
  cor_long <- bind_rows(cor_rows)
  if (nrow(cor_long) == 0) stop("No correlations computed.")
  saveRDS(cor_long, cache_file)
}

cor_long <- cor_long %>%
  mutate(is_artifact_gene = grepl(ARTIFACT_PATTERN, gene),
         is_program_gene  = gene %in% PROGRAM_GENES)

message("\nArtifact-class genes detected: ",
        n_distinct(cor_long$gene[cor_long$is_artifact_gene]))

# =========================================================
# GENE SETS
# =========================================================
gene_sets <- load_gene_sets()
gene_sets <- gene_sets[lengths(gene_sets) >= GS_MIN & lengths(gene_sets) <= GS_MAX]
message("Gene sets retained: ", length(gene_sets))

set_overlap <- tibble(
  pathway = names(gene_sets),
  set_size = lengths(gene_sets),
  n_program_genes = sapply(gene_sets, function(g) sum(g %in% PROGRAM_GENES)),
  program_genes = sapply(gene_sets, function(g)
    paste(intersect(g, PROGRAM_GENES), collapse = ";")),
  frac_artifact_genes = sapply(gene_sets, function(g)
    mean(grepl(ARTIFACT_PATTERN, g)))
)
write_tsv(set_overlap, file.path(out_tables, "11_pathway_program_overlap.tsv"))

# Sets dominated by artifact genes are flagged so a reader can see which
# results the exclusion was protecting against.
message("Gene sets >50% ribosomal/mitochondrial/haemoglobin: ",
        sum(set_overlap$frac_artifact_genes > 0.5))

# =========================================================
# PER TARGET
# =========================================================
for (tgt in TARGETS) {
  message("\n=====================================================")
  message("Target: ", tgt)
  message("=====================================================")
  
  ct <- cor_long %>% filter(target == tgt)
  
  # ---- artifact diagnostic (before exclusion) ----
  art_diag <- ct %>%
    group_by(gene_class = ifelse(is_artifact_gene,
                                 "ribosomal / mito / haemoglobin", "other")) %>%
    summarise(n_genes = n_distinct(gene),
              median_rho = median(rho, na.rm = TRUE),
              q25 = quantile(rho, 0.25, na.rm = TRUE),
              q75 = quantile(rho, 0.75, na.rm = TRUE),
              .groups = "drop")
  write_tsv(art_diag, file.path(out_tables, paste0("11_artifact_diagnostic_", tgt, ".tsv")))
  
  message("\n--- Compositional artifact check ---")
  print(art_diag)
  message("A clearly more negative median for the ribosomal class indicates the ",
          "compositional effect that motivates excluding these genes.")
  
  p_art <- ct %>%
    mutate(gene_class = ifelse(is_artifact_gene,
                               "ribosomal /\nmito / haemoglobin", "other")) %>%
    ggplot(aes(x = gene_class, y = rho, fill = gene_class)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey70") +
    geom_violin(alpha = 0.7, linewidth = 0.3) +
    geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.9) +
    scale_fill_manual(values = c("#F07167", "#9ecae1"), guide = "none") +
    labs(x = NULL, y = paste0("Spearman rho with ", tgt),
         title = paste0("Compositional artifact check: ", tgt),
         subtitle = "Per-gene, per-section correlations before exclusion") +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold"))
  
  ggsave(file.path(out_plots, paste0("11_artifact_diagnostic_", tgt, ".png")), p_art,
         width = 6.5, height = 4.5, dpi = 600, bg = "white")
  ggsave(file.path(out_plots, paste0("11_artifact_diagnostic_", tgt, ".pdf")), p_art,
         width = 6.5, height = 4.5, bg = "white", device = cairo_pdf)
  
  # ---- exclusion ----
  if (EXCLUDE_ARTIFACT_GENES) {
    n_before <- n_distinct(ct$gene)
    ct <- ct %>% filter(!is_artifact_gene)
    message("Excluded ", n_before - n_distinct(ct$gene),
            " artifact-class genes from the ranking")
  }
  
  # ---- combined ranking ----
  ranking <- ct %>%
    group_by(gene) %>%
    summarise(n_sections = n(),
              median_rho = median(rho, na.rm = TRUE),
              z = combine_rho(rho, n_spots),
              .groups = "drop") %>%
    filter(n_sections >= MIN_SECTIONS, is.finite(z))
  
  conc <- ct %>%
    inner_join(ranking %>% select(gene, z), by = "gene") %>%
    group_by(gene) %>%
    summarise(k_concordant = sum(sign(rho) == sign(z), na.rm = TRUE),
              n_tested = n(), .groups = "drop")
  
  ranking <- ranking %>%
    left_join(conc, by = "gene") %>%
    mutate(is_program_gene = gene %in% PROGRAM_GENES) %>%
    arrange(desc(z))
  
  write_tsv(ranking, file.path(out_tables, paste0("11_gene_ranking_", tgt, ".tsv")))
  
  message("\nGenes ranked: ", nrow(ranking))
  message("Top 8 positive:")
  print(ranking %>% slice_head(n = 8) %>%
          select(gene, z, median_rho, k_concordant, n_tested, is_program_gene))
  message("Top 8 negative:")
  print(ranking %>% slice_tail(n = 8) %>%
          select(gene, z, median_rho, k_concordant, n_tested, is_program_gene))
  
  # ---- enrichment ----
  stats_vec <- setNames(ranking$z, ranking$gene)
  fg_raw <- fgsea::fgseaMultilevel(pathways = gene_sets, stats = stats_vec,
                                   minSize = GS_MIN, maxSize = GS_MAX)
  
  fg <- fg_raw %>%
    as_tibble() %>%
    mutate(n_le = lengths(leadingEdge),
           n_le_program = sapply(leadingEdge, function(g) sum(g %in% PROGRAM_GENES)),
           frac_le_program = ifelse(n_le > 0, n_le_program / n_le, NA_real_),
           leadingEdge_chr = sapply(leadingEdge, paste, collapse = ";")) %>%
    select(-leadingEdge) %>%
    left_join(set_overlap %>%
                select(pathway, n_program_genes, program_genes, frac_artifact_genes),
              by = "pathway") %>%
    mutate(evidence_type = ifelse(n_program_genes == 0, "discovery", "consistency")) %>%
    arrange(padj, desc(abs(NES)))
  
  # ---- collapse redundant pathways ----
  sig <- fg %>% filter(padj < PADJ_CUT)
  main_paths <- sig$pathway
  if (COLLAPSE_PATHWAYS && nrow(sig) > 1) {
    coll <- tryCatch(
      fgsea::collapsePathways(fg_raw[fg_raw$padj < PADJ_CUT, ],
                              gene_sets, stats_vec),
      error = function(e) { warning("collapsePathways failed: ",
                                    conditionMessage(e)); NULL })
    if (!is.null(coll)) {
      main_paths <- coll$mainPathways
      message("\nCollapsed ", nrow(sig), " significant pathways to ",
              length(main_paths), " independent ones")
    }
  }
  fg <- fg %>% mutate(is_main = pathway %in% main_paths)
  
  # ---- per-section enrichment, for concordance ----
  sec_fg <- bind_rows(lapply(unique(ct$section_id), function(s) {
    sr <- ct %>% filter(section_id == s, is.finite(rho))
    sv <- setNames(atanh(pmax(pmin(sr$rho, 0.999999), -0.999999)), sr$gene)
    r <- tryCatch(fgsea::fgseaMultilevel(gene_sets, sv,
                                         minSize = GS_MIN, maxSize = GS_MAX),
                  error = function(e) NULL)
    if (is.null(r)) return(NULL)
    as_tibble(r) %>% transmute(section_id = s, pathway, NES_sec = NES)
  }))
  
  if (nrow(sec_fg) > 0) {
    path_conc <- sec_fg %>%
      inner_join(fg %>% select(pathway, NES), by = "pathway") %>%
      group_by(pathway) %>%
      summarise(k_sections_same_sign = sum(sign(NES_sec) == sign(NES), na.rm = TRUE),
                n_sections_tested = sum(is.finite(NES_sec)), .groups = "drop")
    fg <- fg %>% left_join(path_conc, by = "pathway")
    write_tsv(sec_fg, file.path(out_tables,
                                paste0("11_gsea_section_concordance_", tgt, ".tsv")))
  }
  
  write_tsv(fg, file.path(out_tables, paste0("11_gsea_full_", tgt, ".tsv")))
  
  fg_main <- fg %>% filter(padj < PADJ_CUT, is_main) %>% arrange(desc(abs(NES)))
  write_tsv(fg_main, file.path(out_tables, paste0("11_gsea_main_", tgt, ".tsv")))
  
  discovery <- fg_main %>% filter(evidence_type == "discovery")
  write_tsv(discovery, file.path(out_tables, paste0("11_gsea_discovery_", tgt, ".tsv")))
  
  n_pos <- sum(fg_main$NES > 0, na.rm = TRUE)
  n_neg <- sum(fg_main$NES < 0, na.rm = TRUE)
  message("\nMain pathways (padj < ", PADJ_CUT, "): ", nrow(fg_main),
          "  [", n_pos, " positive, ", n_neg, " negative]")
  message("Of these, discovery (no program-gene overlap): ", nrow(discovery))
  if (n_pos == 0 || n_neg == 0) {
    message("NOTE: all enrichment runs in one direction. Check the artifact ",
            "diagnostic above - a one-sided result usually means a global ",
            "composition effect rather than a biological program.")
  }
  if (nrow(fg_main) > 0) {
    print(fg_main %>% slice_head(n = 12) %>%
            select(pathway, NES, padj, size, evidence_type,
                   k_sections_same_sign, n_sections_tested))
  }
  
  # ---- plots ----
  plot_fg <- function(df, title_txt, subtitle_txt, file_stem) {
    if (nrow(df) == 0) { message("Nothing to plot for ", file_stem); return(invisible(NULL)) }
    d <- df %>% slice_head(n = N_TOP_PLOT) %>%
      mutate(pathway_lab = pretty_path(pathway),
             conc_lab = if ("k_sections_same_sign" %in% names(.))
               paste0(k_sections_same_sign, "/", n_sections_tested) else "")
    
    p <- ggplot(d, aes(x = reorder(pathway_lab, NES), y = NES, fill = evidence_type)) +
      geom_col(width = 0.72) +
      coord_flip() +
      geom_text(aes(label = conc_lab, hjust = ifelse(NES > 0, -0.15, 1.15)), size = 2.6) +
      scale_fill_manual(values = c(consistency = "#9ecae1", discovery = "#F07167"),
                        name = NULL, drop = FALSE) +
      labs(x = NULL, y = "Normalised enrichment score",
           title = title_txt, subtitle = subtitle_txt) +
      theme_bw(base_size = 10) +
      theme(plot.title = element_text(face = "bold"),
            plot.subtitle = element_text(size = 8),
            panel.grid.major.y = element_blank(),
            legend.position = "bottom")
    
    ggsave(file.path(out_plots, paste0(file_stem, ".png")), p,
           width = 10, height = max(4, 0.32 * nrow(d) + 2), dpi = 600, bg = "white")
    ggsave(file.path(out_plots, paste0(file_stem, ".pdf")), p,
           width = 10, height = max(4, 0.32 * nrow(d) + 2),
           bg = "white", device = cairo_pdf)
  }
  
  plot_fg(fg_main,
          paste0("Pathways associated with ", tgt),
          paste0("Redundant pathways collapsed; ribosomal/mito/haemoglobin genes excluded. ",
                 "Blue: contains scoring-program genes. Red: no overlap. ",
                 "Labels: sections with concordant sign."),
          paste0("11_gsea_main_", tgt))
  
  plot_fg(discovery,
          paste0("Associations independent of the scoring genes: ", tgt),
          "Main pathways containing no gene from the scoring programs",
          paste0("11_gsea_discovery_", tgt))
}

writeLines(c(
  paste0("targets: ", paste(TARGETS, collapse = ",")),
  paste0("gene_min_detect_per_section: ", GENE_MIN_DETECT),
  paste0("min_sections_per_gene: ", MIN_SECTIONS),
  paste0("exclude_artifact_genes: ", EXCLUDE_ARTIFACT_GENES),
  paste0("artifact_pattern: ", ARTIFACT_PATTERN),
  paste0("collections: ",
         paste(c(if (USE_HALLMARK) "H", if (USE_REACTOME) "C2:CP:REACTOME"),
               collapse = ",")),
  paste0("gene_set_size: ", GS_MIN, "-", GS_MAX),
  paste0("collapse_pathways: ", COLLAPSE_PATHWAYS),
  paste0("padj_cut: ", PADJ_CUT),
  "ranking: within-section Spearman, weighted Stouffer across sections",
  "program genes retained in the ranking and annotated, not removed",
  paste0("seed: ", SEED),
  paste0("run_date: ", as.character(Sys.Date()))
), file.path(out_tables, "11_run_parameters.txt"))

message("\nDone. Outputs written to: ", out_root)