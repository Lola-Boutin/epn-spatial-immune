# =============================================================================
# figure6_curated.R — Pathway associations of the opportunity scores
#
#   Figure 6a : curated pathways x the 3 composite opportunity scores
#   Figure 6b : curated pathways x the 5 gamma-delta mechanism axes
#   Figure 6c : running enrichment, IFN-alpha response under the butyrophilin
#               axis vs under the mevalonate / IPP axis
#
# Main-figure heatmaps are deliberately curated rather than displaying the full
# union of significant GSEA terms. The complete GSEA results remain in the TSVs.
#
# Reads the per-target tables written by 11_gsea_characterization.R.
# Writes to Manuscript_2_EPN. Nothing in phase 11 needs to be rerun.
#
# NOTE ON THE OVERLAP ANNOTATION
# ------------------------------
# Phase 11 computes evidence_type against the POOLED PROGRAM_GENES vector, so
# it is identical for all eight targets. That can mislabel, for example,
# IFN-alpha response under the butyrophilin axis as "consistency" because the
# pathway contains abT_core genes (B2M, PSMB8/9, IRF1), none of which belong to
# the BTN axis. This script therefore recomputes overlap per (score, pathway)
# from SCORE_GENES below.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(forcats)
  library(ggplot2)
  library(patchwork)
  library(fgsea)
})

# -----------------------------------------------------------------------------
# 0. Paths
# -----------------------------------------------------------------------------
rerun_root  <- "D:/Ped-CNS_KBH/Spatial Transcriptomic/GSE195661/Analysis"
gsea_tables <- file.path(rerun_root, "11_gsea_characterization", "tables")
figure_dir  <- "D:/Ped-CNS_KBH/Manuscript_2_EPN/figure6"
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 1. Config
# -----------------------------------------------------------------------------
PADJ_CUT  <- 0.05
MAIN_ONLY <- TRUE    # only collapsePathways main pathways are coloured
CONC_CUT  <- 0.75    # sign concordance below this is drawn at low opacity
GS_MIN    <- 10      # must match phase 11
GS_MAX    <- 500

composite_scores <- c(
  abT_core           = "αβT ligand availability",
  gdT_core           = "γδT ligand availability",
  inhibitory_program = "Immunosuppressive tone"
)

# NKG2D is intentionally placed last. Its only selected significant main pathway
# is gluconeogenesis and is poorly concordant across sections (8/14), so it acts
# as a useful visual contrast rather than interrupting the stronger axes.
axis_scores <- c(
  gdT_btn           = "Butyrophilin",
  gdT_mevalonate    = "Mevalonate / IPP",
  gdT_ephrin        = "EPHA2",
  gdT_adhesion_dnam = "DNAM-1 / adhesion",
  gdT_nkg2d         = "NKG2D ligands"
)

all_scores <- c(composite_scores, axis_scores)

# Panel c ---------------------------------------------------------------------
PANEL_C_PATHWAY <- "HALLMARK_INTERFERON_ALPHA_RESPONSE"
PANEL_C_SCORES  <- c("gdT_btn", "gdT_mevalonate")
PANEL_C_YLIM    <- c(-0.15, 0.65)  # shared scale makes the BTN/mevalonate contrast directly comparable

# -----------------------------------------------------------------------------
# 2. Curated pathway sets for the MAIN figure
# -----------------------------------------------------------------------------
# These vectors define BOTH inclusion and top-to-bottom row order.
# A pathway is only coloured in a cell when it passes PADJ_CUT and MAIN_ONLY for
# that score; non-significant score/pathway combinations stay grey.
#
# Full GSEA results remain in the phase-11 TSV files, so the main figure can be
# selective and mechanistic rather than exhaustive.

MAIN_PATHWAYS_A <- c(
  # Interferon & antigen presentation
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "REACTOME_ANTIGEN_PRESENTATION_FOLDING_ASSEMBLY_AND_PEPTIDE_LOADING_OF_CLASS_I_MHC",
  "REACTOME_IMMUNOREGULATORY_INTERACTIONS_BETWEEN_A_LYMPHOID_AND_A_NON_LYMPHOID_CELL",
  
  # Inflammation & immunoregulation
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "REACTOME_TGF_BETA_RECEPTOR_SIGNALING_ACTIVATES_SMADS",
  
  # Isoprenoid & glycosylation
  "REACTOME_CHOLESTEROL_BIOSYNTHESIS",
  "REACTOME_LANOSTEROL_BIOSYNTHESIS",
  "REACTOME_BIOSYNTHESIS_OF_THE_N_GLYCAN_PRECURSOR_DOLICHOL_LIPID_LINKED_OLIGOSACCHARIDE_LLO_AND_TRANSFER_TO_A_NASCENT_PROTEIN",
  
  # Hypoxia, glycolysis & stress
  "HALLMARK_HYPOXIA",
  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_APOPTOSIS",
  
  # Oxidative metabolism
  "REACTOME_FORMATION_OF_ATP_BY_CHEMIOSMOTIC_COUPLING",
  
  # ECM & adhesion
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",
  "REACTOME_CELL_SURFACE_INTERACTIONS_AT_THE_VASCULAR_WALL"
)

MAIN_PATHWAYS_B <- c(
  # Interferon & antigen presentation
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "REACTOME_ANTIGEN_PRESENTATION_FOLDING_ASSEMBLY_AND_PEPTIDE_LOADING_OF_CLASS_I_MHC",
  "REACTOME_IMMUNOREGULATORY_INTERACTIONS_BETWEEN_A_LYMPHOID_AND_A_NON_LYMPHOID_CELL",
  
  # Inflammation & cytokine signalling
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "REACTOME_DAP12_INTERACTIONS",
  
  # Isoprenoid & glycosylation
  "REACTOME_CHOLESTEROL_BIOSYNTHESIS",
  "REACTOME_LANOSTEROL_BIOSYNTHESIS",
  "REACTOME_DAG1_GLYCOSYLATIONS",
  
  # Hypoxia, glycolysis & stress
  "HALLMARK_HYPOXIA",
  "REACTOME_GLUCOSE_METABOLISM",
  "HALLMARK_APOPTOSIS",
  
  # Oxidative metabolism
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "REACTOME_FORMATION_OF_ATP_BY_CHEMIOSMOTIC_COUPLING",
  
  # ECM & adhesion
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",
  "REACTOME_CELL_SURFACE_INTERACTIONS_AT_THE_VASCULAR_WALL",
  
  # Sole significant NKG2D main-pathway hit; intentionally retained despite
  # low cross-section concordance so the weak/inconsistent NKG2D signal is shown.
  "REACTOME_GLUCONEOGENESIS"
)

# -----------------------------------------------------------------------------
# 3. Per-score gene programs
#    Genes as declared in 09_opportunity_map.R, split by the score they build.
#    Used ONLY for the overlap annotation; they were never removed from ranking.
#    Alias variants are listed so the annotation matches phase 11.
# -----------------------------------------------------------------------------
SCORE_GENES <- list(
  abT_core           = c("HLA-A","HLA-B","HLA-C","B2M","TAP1","TAP2",
                         "PSMB8","PSMB9","ERAP1","ERAP2","NLRC5","IRF1"),
  gdT_core           = c("BTN2A1","BTN3A1","MVK","PMVK","MVD","IDI1",
                         "MICA","MICB","ULBP2","ULBP3","ULBP4","RAET1E",
                         "ULBP5","RAET1G","ULBP6","RAET1L",
                         "EPHA2","PVR","NECTIN2","ICAM1"),
  inhibitory_program = c("CD274","PDCD1LG2","HLA-E","LGALS9","IDO1",
                         "TGFB1","TGFB2","TGFB3","IL10","VEGFA",
                         "NT5E","ENTPD1","PTGS2","PTGES","PTGER2","PTGER4"),
  gdT_btn            = c("BTN2A1","BTN3A1"),
  gdT_mevalonate     = c("MVK","PMVK","MVD","IDI1"),
  gdT_nkg2d          = c("MICA","MICB","ULBP2","ULBP3","ULBP4","RAET1E",
                         "ULBP5","RAET1G","ULBP6","RAET1L"),
  gdT_ephrin         = c("EPHA2"),
  gdT_adhesion_dnam  = c("PVR","NECTIN2","ICAM1")
)
SCORE_GENES <- lapply(SCORE_GENES, toupper)

# Pathways with NO shared genes but a direct metabolic dependence on the score's
# program. Manual, mechanistic, and declared as such in the legend.
# Dolichol is synthesised from farnesyl-PP, downstream of the mevalonate enzymes.
DOWNSTREAM_OF <- list(
  gdT_mevalonate = c(
    "REACTOME_SYNTHESIS_OF_DOLICHYL_PHOSPHATE_MANNOSE",
    "REACTOME_SYNTHESIS_OF_SUBSTRATES_IN_N_GLYCAN_BIOSYTHESIS",
    "REACTOME_ASPARAGINE_N_LINKED_GLYCOSYLATION",
    "REACTOME_DISEASES_OF_GLYCOSYLATION",
    "REACTOME_GLYCOSAMINOGLYCAN_METABOLISM"
  )
)
DOWNSTREAM_OF$gdT_core <- DOWNSTREAM_OF$gdT_mevalonate

# -----------------------------------------------------------------------------
# 4. Row grouping and publication labels
# -----------------------------------------------------------------------------
# Only selected pathways need to be grouped for the main figure. Unlisted terms
# still fall through to "Other" as a safety check.
PATHWAY_GROUP <- c(
  # Interferon & antigen presentation
  HALLMARK_INTERFERON_GAMMA_RESPONSE       = "Interferon & antigen presentation",
  HALLMARK_INTERFERON_ALPHA_RESPONSE       = "Interferon & antigen presentation",
  REACTOME_ANTIGEN_PRESENTATION_FOLDING_ASSEMBLY_AND_PEPTIDE_LOADING_OF_CLASS_I_MHC =
    "Interferon & antigen presentation",
  REACTOME_IMMUNOREGULATORY_INTERACTIONS_BETWEEN_A_LYMPHOID_AND_A_NON_LYMPHOID_CELL =
    "Interferon & antigen presentation",
  
  # Inflammation & cytokine signalling
  HALLMARK_TNFA_SIGNALING_VIA_NFKB                    = "Inflammation & cytokine signaling",
  HALLMARK_INFLAMMATORY_RESPONSE                      = "Inflammation & cytokine signaling",
  HALLMARK_IL6_JAK_STAT3_SIGNALING                     = "Inflammation & cytokine signaling",
  REACTOME_TGF_BETA_RECEPTOR_SIGNALING_ACTIVATES_SMADS = "Inflammation & cytokine signaling",
  REACTOME_DAP12_INTERACTIONS                          = "Inflammation & cytokine signaling",
  
  # Isoprenoid & glycosylation
  REACTOME_CHOLESTEROL_BIOSYNTHESIS = "Isoprenoid & glycosylation",
  REACTOME_LANOSTEROL_BIOSYNTHESIS   = "Isoprenoid & glycosylation",
  REACTOME_BIOSYNTHESIS_OF_THE_N_GLYCAN_PRECURSOR_DOLICHOL_LIPID_LINKED_OLIGOSACCHARIDE_LLO_AND_TRANSFER_TO_A_NASCENT_PROTEIN =
    "Isoprenoid & glycosylation",
  REACTOME_DAG1_GLYCOSYLATIONS = "Isoprenoid & glycosylation",
  
  # Hypoxia / glycolysis / stress
  HALLMARK_HYPOXIA            = "Hypoxia, glycolysis & stress",
  HALLMARK_GLYCOLYSIS         = "Hypoxia, glycolysis & stress",
  HALLMARK_APOPTOSIS          = "Hypoxia, glycolysis & stress",
  REACTOME_GLUCOSE_METABOLISM = "Hypoxia, glycolysis & stress",
  REACTOME_GLUCONEOGENESIS    = "Hypoxia, glycolysis & stress",
  
  # Oxidative metabolism
  HALLMARK_OXIDATIVE_PHOSPHORYLATION                 = "Oxidative metabolism",
  REACTOME_FORMATION_OF_ATP_BY_CHEMIOSMOTIC_COUPLING = "Oxidative metabolism",
  
  # ECM / adhesion
  HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION          = "ECM & adhesion",
  REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION          = "ECM & adhesion",
  REACTOME_CELL_SURFACE_INTERACTIONS_AT_THE_VASCULAR_WALL = "ECM & adhesion"
)

GROUP_LEVELS <- c(
  "Interferon & antigen presentation",
  "Inflammation & cytokine signaling",
  "Isoprenoid & glycosylation",
  "Hypoxia, glycolysis & stress",
  "Oxidative metabolism",
  "ECM & adhesion",
  "Other"
)

# Override table keyed on RAW pathway IDs, so gene symbols/complexes are not
# mangled by sentence-casing. Long Reactome names are shortened where useful.
LABEL_FIXES <- c(
  HALLMARK_TNFA_SIGNALING_VIA_NFKB          = "TNF-α signaling via NF-κB",
  HALLMARK_IL6_JAK_STAT3_SIGNALING          = "IL6-JAK-STAT3 signaling",
  HALLMARK_INTERFERON_GAMMA_RESPONSE        = "Interferon-γ response",
  HALLMARK_INTERFERON_ALPHA_RESPONSE        = "Interferon-α response",
  HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION = "Epithelial-mesenchymal transition",
  HALLMARK_OXIDATIVE_PHOSPHORYLATION        = "Oxidative phosphorylation",
  
  REACTOME_ANTIGEN_PRESENTATION_FOLDING_ASSEMBLY_AND_PEPTIDE_LOADING_OF_CLASS_I_MHC =
    "Antigen presentation: folding, assembly and peptide loading of class I MHC",
  REACTOME_IMMUNOREGULATORY_INTERACTIONS_BETWEEN_A_LYMPHOID_AND_A_NON_LYMPHOID_CELL =
    "Immunoregulatory interactions between lymphoid and non-lymphoid cells",
  REACTOME_TGF_BETA_RECEPTOR_SIGNALING_ACTIVATES_SMADS =
    "TGF-β receptor signaling activates SMADs",
  REACTOME_CHOLESTEROL_BIOSYNTHESIS = "Cholesterol biosynthesis",
  REACTOME_LANOSTEROL_BIOSYNTHESIS   = "Lanosterol biosynthesis",
  REACTOME_BIOSYNTHESIS_OF_THE_N_GLYCAN_PRECURSOR_DOLICHOL_LIPID_LINKED_OLIGOSACCHARIDE_LLO_AND_TRANSFER_TO_A_NASCENT_PROTEIN =
    "N-glycan precursor biosynthesis (dolichol-linked)",
  REACTOME_DAG1_GLYCOSYLATIONS = "DAG1 glycosylation",
  REACTOME_GLUCOSE_METABOLISM  = "Glucose metabolism",
  REACTOME_GLUCONEOGENESIS     = "Gluconeogenesis",
  REACTOME_FORMATION_OF_ATP_BY_CHEMIOSMOTIC_COUPLING =
    "Formation of ATP by chemiosmotic coupling",
  REACTOME_DAP12_INTERACTIONS = "DAP12 interactions",
  REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION = "Extracellular matrix organization",
  REACTOME_CELL_SURFACE_INTERACTIONS_AT_THE_VASCULAR_WALL =
    "Cell-surface interactions at the vascular wall"
)

pretty_pathway <- function(x, width = 38) {
  base <- gsub("_", " ", sub("^(HALLMARK|REACTOME)_", "", x))
  base <- paste0(toupper(substr(base, 1, 1)), tolower(substr(base, 2, nchar(base))))
  out  <- ifelse(x %in% names(LABEL_FIXES), unname(LABEL_FIXES[x]), base)
  str_wrap(out, width = width)          # wrap, never truncate
}

# Staggered (quincunx) column headers -----------------------------------------
stagger_labels <- function(x, width = 18) {
  x <- str_wrap(x, width = width)
  ifelse(seq_along(x) %% 2 == 1, paste0(x, "\n"), paste0("\n", x))
}

# -----------------------------------------------------------------------------
# 5. Load the per-target GSEA tables
# -----------------------------------------------------------------------------
read_target <- function(tgt) {
  f <- file.path(gsea_tables, paste0("11_gsea_full_", tgt, ".tsv"))
  if (!file.exists(f)) stop("Missing phase-11 table: ", f)
  read_tsv(f, show_col_types = FALSE) %>% mutate(score = tgt)
}

res <- bind_rows(lapply(names(all_scores), read_target))

stopifnot(all(c("pathway", "NES", "pval", "padj", "size", "is_main",
                "k_sections_same_sign", "n_sections_tested") %in% names(res)))

# Validate curated pathway IDs before plotting. This catches typos or changes in
# MSigDB naming rather than silently producing missing rows.
curated_all <- unique(c(MAIN_PATHWAYS_A, MAIN_PATHWAYS_B))
missing_curated <- setdiff(curated_all, unique(res$pathway))
if (length(missing_curated) > 0) {
  stop("Curated pathway(s) absent from GSEA tables:\n  ",
       paste(missing_curated, collapse = "\n  "))
}

# Pathway -> gene mapping, rebuilt with the same filters as phase 11 so set
# membership is identical to what fgsea saw.
load_gene_sets <- function() {
  if (!requireNamespace("msigdbr", quietly = TRUE))
    stop("msigdbr not installed.")
  
  get_sets <- function(...) {
    args <- list(species = "Homo sapiens", ...)
    tryCatch(
      do.call(msigdbr::msigdbr, args),
      error = function(e) {
        names(args)[names(args) == "collection"]    <- "category"
        names(args)[names(args) == "subcollection"] <- "subcategory"
        do.call(msigdbr::msigdbr, args)
      }
    )
  }
  
  h  <- get_sets(collection = "H")
  r  <- get_sets(collection = "C2", subcollection = "CP:REACTOME")
  gs <- c(split(toupper(h$gene_symbol), h$gs_name),
          split(toupper(r$gene_symbol), r$gs_name))
  gs <- lapply(gs, unique)
  gs[lengths(gs) >= GS_MIN & lengths(gs) <= GS_MAX]
}

gene_sets <- load_gene_sets()

classify_overlap <- function(score_id, pathway_id) {
  gs <- gene_sets[[pathway_id]]
  if (is.null(gs)) return(NA_character_)
  if (any(gs %in% SCORE_GENES[[score_id]])) return("contains")
  dn <- DOWNSTREAM_OF[[score_id]]
  if (!is.null(dn) && pathway_id %in% dn) return("downstream")
  "independent"
}

res <- res %>%
  rowwise() %>%
  mutate(overlap = classify_overlap(score, pathway)) %>%
  ungroup()

# How often does the corrected annotation disagree with phase 11? -------------
if ("evidence_type" %in% names(res)) {
  d <- res %>%
    filter(padj < PADJ_CUT, is_main) %>%
    mutate(phase11 = ifelse(evidence_type == "discovery", "independent", "contains"))
  message("Overlap annotation corrected in ",
          sum(d$overlap != d$phase11, na.rm = TRUE), " of ", nrow(d),
          " significant cells")
}

# -----------------------------------------------------------------------------
# 6. Prepare a CURATED panel
# -----------------------------------------------------------------------------
prep_panel <- function(res, score_map, pathway_keep) {
  
  # Safety: each selected pathway should be a significant main pathway in at
  # least one score of the panel. Grey cells in the other columns are retained.
  sig_union <- res %>%
    filter(score %in% names(score_map),
           padj < PADJ_CUT,
           !MAIN_ONLY | is_main) %>%
    distinct(pathway) %>%
    pull(pathway)
  
  selected_not_sig <- setdiff(pathway_keep, sig_union)
  if (length(selected_not_sig) > 0) {
    warning("Selected pathway(s) are not significant main pathways in this panel:\n  ",
            paste(selected_not_sig, collapse = "\n  "))
  }
  
  label_levels <- rev(pretty_pathway(pathway_keep))
  
  d <- res %>%
    filter(score %in% names(score_map), pathway %in% pathway_keep) %>%
    mutate(
      sig       = padj < PADJ_CUT & (!MAIN_ONLY | is_main),
      conc_frac = k_sections_same_sign / n_sections_tested,
      score_lab = factor(unname(score_map[score]), levels = unname(score_map)),
      NES_plot  = ifelse(sig, NES, NA_real_),
      conc_bin  = factor(
        ifelse(!sig, NA_character_,
               ifelse(conc_frac >= CONC_CUT, "≥75%", "<75%")),
        levels = c("≥75%", "<75%")
      ),
      cell_lab = ifelse(sig,
                        paste0(k_sections_same_sign, "/", n_sections_tested),
                        NA_character_),
      ov_plot = factor(
        ifelse(sig, overlap, NA_character_),
        levels = c("contains", "downstream", "independent")
      ),
      group = factor(
        ifelse(pathway %in% names(PATHWAY_GROUP),
               unname(PATHWAY_GROUP[pathway]), "Other"),
        levels = GROUP_LEVELS
      ),
      lab = factor(pretty_pathway(pathway), levels = label_levels),
      pathway_order = match(pathway, pathway_keep)
    ) %>%
    arrange(pathway_order, score_lab)
  
  if (any(d$group == "Other")) {
    warning("Selected pathways falling through to 'Other': ",
            paste(unique(d$pathway[d$group == "Other"]), collapse = ", "))
  }
  
  d
}

dat_a <- prep_panel(res, composite_scores, MAIN_PATHWAYS_A)
dat_b <- prep_panel(res, axis_scores,      MAIN_PATHWAYS_B)

# One common NES scale across panels a and b.
nes_lim <- max(abs(c(dat_a$NES_plot, dat_b$NES_plot)), na.rm = TRUE)

# -----------------------------------------------------------------------------
# 7. Heatmap
#    fill    = NES
#    opacity = cross-section sign concordance
#    border  = overlap of the pathway with THAT score's gene program
# -----------------------------------------------------------------------------
OVERLAP_COLS <- c(
  contains    = "#08519C",
  downstream  = "#E08214",
  independent = "grey75"
)

build_heatmap <- function(d, lim, title_txt) {
  ggplot(d, aes(x = score_lab, y = lab)) +
    geom_tile(fill = "grey90", colour = "white", linewidth = 0.6) +
    geom_tile(
      data = filter(d, sig),
      aes(fill = NES_plot, alpha = conc_bin, colour = ov_plot),
      linewidth = 0.7
    ) +
    geom_text(aes(label = cell_lab), size = 2.5, colour = "grey15") +
    scale_fill_gradient2(
      low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
      midpoint = 0, limits = c(-lim, lim), name = "NES",
      na.value = "grey90"
    ) +
    scale_alpha_manual(
      values = c("≥75%" = 1, "<75%" = 0.45),
      name = "Sign concordance\nacross sections",
      na.translate = FALSE, drop = FALSE
    ) +
    scale_colour_manual(
      values = OVERLAP_COLS,
      labels = c(
        contains    = "contains program genes",
        downstream  = "downstream of program",
        independent = "independent"
      ),
      name = "Overlap with that\nscore's gene program",
      na.translate = FALSE, drop = TRUE
    ) +
    scale_x_discrete(
      position = "top",
      labels = function(x) stagger_labels(x),
      drop = FALSE
    ) +
    scale_y_discrete(expand = c(0, 0), drop = TRUE) +
    facet_grid(
      group ~ .,
      scales = "free_y",
      space = "free_y",
      drop = TRUE,
      labeller = labeller(group = label_wrap_gen(18))
    ) +
    guides(
      colour = guide_legend(
        override.aes = list(fill = "white", linewidth = 1)
      )
    ) +
    labs(x = NULL, y = NULL, title = title_txt) +
    theme_minimal(base_size = 9) +
    theme(
      plot.title        = element_text(face = "bold", size = 10),
      panel.grid        = element_blank(),
      axis.text.x.top   = element_text(vjust = 0.5, lineheight = 0.95),
      axis.text.y       = element_text(lineheight = 0.95),
      strip.text.y      = element_text(angle = 0, face = "bold", hjust = 0),
      strip.background  = element_rect(fill = "grey93", colour = NA),
      panel.spacing.y   = unit(4, "pt"),
      legend.key.height = unit(11, "pt")
    )
}

panel_a <- build_heatmap(dat_a, nes_lim, "Composite opportunity scores")
panel_b <- build_heatmap(dat_b, nes_lim, "γδT mechanism axes")

# -----------------------------------------------------------------------------
# 8. Panel c — running enrichment
#    Stats rebuilt from the phase-11 ranking tables; identical to what fgsea saw.
# -----------------------------------------------------------------------------
load_stats <- function(tgt) {
  f <- file.path(gsea_tables, paste0("11_gene_ranking_", tgt, ".tsv"))
  if (!file.exists(f)) stop("Missing ranking table: ", f)
  rk <- read_tsv(f, show_col_types = FALSE) %>% arrange(desc(z))
  setNames(rk$z, rk$gene)
}

enrichment_panel <- function(score_id, pathway_id) {
  stats <- load_stats(score_id)
  gs    <- gene_sets[[pathway_id]]
  if (is.null(gs)) stop("Pathway not found in gene sets: ", pathway_id)
  
  pd <- plotEnrichmentData(pathway = gs, stats = stats)
  row <- res %>%
    filter(score == score_id, pathway == pathway_id) %>%
    slice(1)
  
  has_stats <- nrow(row) == 1 && is.finite(row$NES) && is.finite(row$padj)
  
  sub <- if (nrow(row) == 0) {
    "pathway not present in this score's table"
  } else if (!has_stats) {
    if (is.finite(row$n_sections_tested)) {
      sprintf(
        "No reproducible enrichment across sections (%d/%d concordant)",
        row$k_sections_same_sign, row$n_sections_tested
      )
    } else {
      "No reproducible enrichment across sections"
    }
  } else {
    p_txt <- if (row$padj < 1e-4) {
      "adjusted p < 1e-4"
    } else {
      sprintf("adjusted p = %.2g", row$padj)
    }
    
    sprintf(
      "NES = %.2f, %s; concordant in %d/%d sections%s",
      row$NES, p_txt,
      row$k_sections_same_sign, row$n_sections_tested,
      if (row$padj >= PADJ_CUT) " — not significant" else ""
    )
  }
  
  tick_y <- -0.025
  peak_ES <- pd$curve$ES[which.max(abs(pd$curve$ES))]
  
  ggplot() +
    geom_hline(yintercept = 0, colour = "grey20", linewidth = 0.4) +
    geom_hline(
      yintercept = peak_ES,
      linetype = "dashed", colour = "grey45", linewidth = 0.4
    ) +
    geom_line(
      data = pd$curve,
      aes(x = rank, y = ES),
      colour = "#3B7A57", linewidth = 0.6
    ) +
    geom_segment(
      data = pd$ticks,
      aes(x = rank, xend = rank, y = tick_y, yend = tick_y * 0.2),
      colour = "grey25", linewidth = 0.3
    ) +
    labs(
      title = paste0(
        pretty_pathway(pathway_id, width = 60), " — ",
        unname(all_scores[score_id])
      ),
      subtitle = sub,
      x = "Gene rank (combined across sections)",
      y = "Running enrichment score"
    ) +
    coord_cartesian(ylim = PANEL_C_YLIM) +
    theme_minimal(base_size = 9) +
    theme(
      plot.title       = element_text(face = "bold", size = 9),
      plot.subtitle    = element_text(colour = "grey35", size = 8),
      panel.grid.minor = element_blank()
    )
}

panel_c <- enrichment_panel(PANEL_C_SCORES[1], PANEL_C_PATHWAY) +
  enrichment_panel(PANEL_C_SCORES[2], PANEL_C_PATHWAY)

# -----------------------------------------------------------------------------
# 9. Save
# -----------------------------------------------------------------------------
save_panel <- function(plot, name, width, height) {
  ggsave(
    file.path(figure_dir, paste0(name, ".png")), plot,
    width = width, height = height, dpi = 600, units = "in", bg = "white"
  )
  ggsave(
    file.path(figure_dir, paste0(name, ".pdf")), plot,
    width = width, height = height, units = "in", device = cairo_pdf
  )
}

n_rows <- function(d) n_distinct(d$pathway)

# More compact than the original 0.26 inch/row because panels a/b are now curated.
heatmap_height <- function(d) max(4.2, 0.20 * n_rows(d) + 2.2)

save_panel(
  panel_a, "Figure6a_composite_pathway_heatmap",
  7.5, heatmap_height(dat_a)
)
save_panel(
  panel_b, "Figure6b_gdT_axis_pathway_heatmap",
  9.0, heatmap_height(dat_b)
)
save_panel(panel_c, "Figure6c_running_enrichment", 9.0, 3.2)

# -----------------------------------------------------------------------------
# 10. Values / QC to print
# -----------------------------------------------------------------------------
cat("\n--- Curated main-figure pathway counts ---\n")
cat("Figure 6a:", n_rows(dat_a), "pathways\n")
cat("Figure 6b:", n_rows(dat_b), "pathways\n")

cat("\n--- BH family check (whichever difference is ~0 is the family) ---\n")
res %>%
  group_by(score) %>%
  mutate(padj_within = p.adjust(pval, "BH")) %>%
  ungroup() %>%
  mutate(padj_global = p.adjust(pval, "BH")) %>%
  summarise(
    diff_within = max(abs(padj - padj_within), na.rm = TRUE),
    diff_global = max(abs(padj - padj_global), na.rm = TRUE)
  ) %>%
  print()

cat("\n--- Significant cells in CURATED panels by concordance and overlap ---\n")
bind_rows(
  mutate(dat_a, panel = "a"),
  mutate(dat_b, panel = "b")
) %>%
  filter(sig) %>%
  count(panel, conc_bin, ov_plot) %>%
  print(n = Inf)

cat("\n--- Selected NKG2D cells in panel b ---\n")
dat_b %>%
  filter(score == "gdT_nkg2d", sig) %>%
  select(pathway, NES, padj, k_sections_same_sign, n_sections_tested,
         conc_frac, overlap) %>%
  print(n = Inf)

cat("\n--- Selected pathways falling through to 'Other' ---\n")
bind_rows(dat_a, dat_b) %>%
  filter(group == "Other") %>%
  distinct(pathway) %>%
  print(n = Inf)

message("\nDone. Figures written to: ", figure_dir)
