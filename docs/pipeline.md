# Pipeline: run order and dependencies

Reconstructed from the scripts as they currently stand. Rows marked ⚠ need
confirming against `Analysis/` on disk.

Every script sources `config.R` and resolves paths through `stage_dir()`,
`stage_input()` and `ref_file()`. No script contains a literal path.

---

## Issues to resolve before publishing

**-1. Stage 06 must be run TWICE — RESOLVED.**
`06_zones_from_corrected_backbone.R` now carries an `ALL_SECTIONS_RUN` switch:

| Setting | Sections | Stage key | Consumed by |
|---|---|---|---|
| `FALSE` (default) | `GOOD_SECTIONS` (6) | `zones` | 07, 08 |
| `TRUE` | `ALL_SECTIONS` (14) | `zones_all` | 09, 10 |

Run it once each way. Hotspot calls exist only for `GOOD_SECTIONS`, so the
all-sections run leaves hotspot context empty (`is_hotspot_semla = FALSE`) for
the other eight — which is what 09 and 10 expect, since they use `zone_call`
rather than hotspot status. The hard `stop()` on a missing hotspot table now
applies only when `ALL_SECTIONS_RUN` is `FALSE`, so the six-section run is
unchanged.

**0. Stage 04 no longer exists — RESOLVED.**
`04_lm7_deconvolution.R` implemented the per-spot UCell approach and has been
deleted. Its only output, `vis_good_lm7_ready.rds`, was the merged good-sections
object plus UCell scores; stage 05 used none of those scores. Stage 05 now reads
the `02c` per-section objects directly, which already carry `is_hotspot_semla`,
`local_lymph_score` and `in_tissue`. There is no stage 04 in the pipeline and no
orphan RDS to deposit. ⚠ Delete any leftover `04_lm7_deconvolution/` folder.

**1. `05` writes to the wrong folder — RESOLVED, fix pending.**
The correct folder is `Analysis/04_hotspot_pseudobulk_lm7`, which is what
`05b_plot_DeepTIL_stackedbars.R` and `05c_plot_DeepTIL_hotspot_vs_nonhotspot.R`
already read. `05_hotspot_pseudobulk_CIBERSORTx.R` hardcodes
`04b_hotspot_pseudobulk_lm7` and is the file that must change:

```r
# line 57 — wrong
out_dir <- file.path(analysis_root, "04b_hotspot_pseudobulk_lm7")
# replace with
out_dir <- stage_dir("pseudobulk")
```

Run fresh today, `05b`/`05c` would fail on a missing input. ⚠ Also check for a
stale `04b_` folder on disk holding an older copy of the pseudobulks; delete it
rather than keep two.

**2. Script number vs folder number collisions.**

| Script | Writes to folder | Problem |
|---|---|---|
| `04_he_overlays_from_png_tissue_hires.R` | `03c_he_overlays_from_png_tissue_v2` | `03c` is also a script name |
| `03b_run_semla_nnls_from_augmented.R` | `02-4_semla_nnls_tumor_augmented_PFA` | number mismatch |
| `03c_compare_original_vs_tumor_augmented_hotspots.R` | `Hotspot_tumor_augmented_PFA_robustness` | unnumbered |
| `06b_qc_signature_independence.R` | `qc_signature_overlap` | unnumbered |
| `05_hotspot_pseudobulk_CIBERSORTx.R` | `04b_hotspot_pseudobulk_lm7` | number mismatch |

Fix by editing values in `STAGE` and moving each folder once. Scripts do not
change, because they refer to keys.

**3. `08_lr_interactions.R` permutations -- SEEDED, re-run needed.**
`use_seed("lr_permutation")` is now called at the start of
`run_resource_permutations`, so the null is reproducible and independent of what
else consumed the RNG first. ⚠ The published significance calls were produced
UNSEEDED. Re-run stage 08 and confirm the calls are unchanged before relying on
them.

**6. Stage 08 uses five sections, not six.**
`08_lr_interactions.R` runs on 459, 812, 821, 928, 1239 -- `723` is absent,
unlike every other stage. `08b` lists all six. ⚠ Confirm the reason and record it
in the script header; the literal has been kept explicit rather than unified to
`GOOD_SECTIONS` so the difference stays visible.

**4. Filename style.** `02b_override-after-qc.R` and `08b_qc-after-phase08.R`
use hyphens; everything else uses underscores. Output file prefixes mix
`phase0`, `phase06_`, `phase07_`, `phase3c_` with numbered directories.

**5. Line endings.** Most scripts are CRLF. Add `.gitattributes`:

```
* text=auto eol=lf
*.R text eol=lf
```

---

## Dependency table

| # | Script | Reads | Writes | Random |
|---|---|---|---|---|
| 00 | `00_phase0_build_manifest_and_raw_sections.R` | `RAW_ROOT` Visium files | `manifest` (`visium_manifest.tsv`, `visium_manifest_batched.tsv`, `phase0_qc_summary.tsv`), `raw_vis` (`vis_section_<sec>_raw.rds`), `raw_spot_tables`, `raw_qc_png` | — |
| 01 | `01_build_raw_barcode_coordinate_maps.R` | manifest, `raw_vis` | `coord_maps` (`barcode_coordinate_map_all.tsv`, per-section, QC) | — |
| 02 | `02_run_semla_nnls_qc_and_save_lymph_only.R` | manifest, coord maps, `raw_vis`, `ref_file("scrna_expr")`, `ref_file("scrna_meta")` | `semla_nnls` — `all_sections/`, `good_sections/`, `qc/section_qc.tsv`, `merged/vis_good_semla_ready.rds`, global maps | — |
| 02b | `02b_override-after-qc.R` | `semla_nnls/all_sections` | `merged/vis_good_semla_ready_manual_override.rds`, revised good_sections | — |
| 02c | `02c_add_hotspots_to_coordinate.R` | `semla_nnls/good_sections`, `raw_spot_tables` | `hotspots` — `vis/`, `nnls/`, `qc/rebuilt_hotspot_summary_tissue.tsv`, QC plots | — |
| 03a | `03a_build_tumor_augmented_reference.R` | `ref_file("tumor_expr")`, `ref_file("tumor_meta")`, `ref_file("scrna_meta")` | `REF$augmented_ref_dir` — `exprMatrix_tumor_augmented_gene_matched.tsv`, `meta_..._gene_matched.tsv`, sampling tables | **Yes** — `SEED <- 1L`, 250 cells/cluster |
| 03b | `03b_run_semla_nnls_from_augmented.R` | augmented reference, manifest, coord maps, `raw_vis` | `semla_nnls_aug` | — |
| 03c | `03c_compare_original_vs_tumor_augmented_hotspots.R` | `semla_nnls`, `hotspots`, `semla_nnls_aug` | `hotspot_robustness` | — |
| 04 | `04_he_overlays_from_png_tissue_hires.R` | `hotspots/nnls`, hires PNGs | `he_overlays` — plots, `phase3c_he_overlay_summary_tissue.tsv` | — |
| 05 | `05_hotspot_pseudobulk_CIBERSORTx.R` | `hotspots/vis/vis_section_<sec>.rds` | `pseudobulk` — `pseudobulk_cpm_CIBERSORTx.tsv`, `pseudobulk_counts.tsv`, `pseudobulk_sample_metadata.tsv` | — |
| — | **external: CIBERSORTx** | mixture TSV + LM7 | fraction tables | not seedable |
| — | **external: DeepTIL / SES** | CIBERSORTx fractions | `deeptil/SES_CIBERSORTx_EPN_pseudobulk.txt` | ⚠ version unrecorded |
| 05b | `05b_plot_DeepTIL_stackedbars.R` | `deeptil` | Fig 2d stacked bars — ⚠ currently writes into `pseudobulk`; should use `figure_path()` | — |
| 05c | `05c_plot_DeepTIL_hotspot_vs_nonhotspot.R` | `deeptil` | Fig 2e paired comparison — ⚠ same, should use `figure_path()` | — |
| 06 | `06_zones_from_corrected_backbone.R` (`ALL_SECTIONS_RUN=FALSE`) | `raw_vis`, `raw_spot_tables`, `hotspots/nnls`, `ref_file("zone_genes_xlsx")` | `zones` — `zones_section_<sec>.tsv/.rds`, `zones_all_sections.tsv`, `zone_hotspot_stats.tsv`, `barrier_neighbor_stats.tsv` | none |
| 06 | `06_zones_from_corrected_backbone.R` (`ALL_SECTIONS_RUN=TRUE`) | same, hotspots optional | `zones_all` — same file set, 14 sections | none |
| 06b | `06b_qc_signature_independence.R` | zone xlsx, scRNA reference | `signature_overlap` — overlap summary, gene lists, UpSet | — |
| 07 | `07_functional_states.R` | `zones/tables`, `raw_vis` | `functional_states` — per-section plotready/retention tables, `phase07_analysis_bundle.rds`, sensitivity QC | none — UCell and the sensitivity grid are deterministic |
| 08 | `08_lr_interactions.R` | `zones/tables`, `raw_vis`, CellChatDB, CellPhoneDB, NicheNet | `lr_interactions` — resources, observed scores, top pairs, permutations, NicheNet, `phase08_lr_results_bundle.rds` | **Yes** — `use_seed("lr_permutation")` |
| 08b | `08b_qc_after_lr_interactions.R` | `lr_interactions` observed/top/perm | `lr_qc` — burden, consistency, outliers | — |
| 09 | `09_opportunity_map.R` | manifest, `raw_vis`, `zones_all/tables` | `scoring` — spot/section scores, gene detection and provenance tables, `rds/09_scoring_inputs.rds`, Fig 5a/5b | none |
| 09b | `09b_scoring_qc.R` | `scoring/rds/09_scoring_inputs.rds`, `raw_vis` (UCell only) | `scoring_qc` — dominance, depth, scheme, axis independence, split-half | **Yes** — `use_seed("split_half")`, 20260812 |
| 10 | `10_zone_and_spatial.R` | `scoring/rds`, `zones_all/tables`, `raw_vis`, `RAW_ROOT` hires images | `zone_spatial` — Fig 5c/5d, depth-residualised and per-axis panels | none |
| 11 | `11_gsea_characterization.R` | manifest, `scoring/rds`, `raw_vis`, msigdbr | `gsea` — per-target rankings, full/main/discovery GSEA tables, artifact diagnostics | **Yes** — `use_seed("gsea")`, 20260812 |
| 11b | `11b_ranking_plots.R` | `gsea/tables` | `gsea` — volcano, concordance, program barcode panels | none |

---

## The two external breaks

The pipeline cannot run unattended end to end. Two steps happen outside R:

1. **CIBERSORTx** — stage 05 writes the mixture file; the run is manual, with
   parameters recorded in `config.R` under `CIBERSORTX`.
2. **DeepTIL / SES** — converts CIBERSORTx fractions into the abundance
   estimates that `05b` and `05c` plot. ⚠ Tool version is not currently
   recorded anywhere; add it to `DEEPTIL$version`.

For tier-2 reproduction, download both sets of results from Zenodo and skip
straight to `05b`.

---

## Script header template

```r
# =============================================================================
# <NN>_<name>.R
#
# Purpose
# -------
# <one paragraph>
#
# Inputs
# ------
#   stage_input("<key>", "<file>")
#   ref_file("<key>")
#
# Outputs
# -------
#   stage_dir("<key>", "<file>")   <one line each>
#
# Stochastic
# ----------
#   use_seed("<key>")   OR   none
#
# Runtime
# -------
#   ~<n> minutes
# =============================================================================

source(here::here("config.R"))
```

---

## Section selection: the manual step

The automatic screen in stage 02 applies `THRESH$min_nonNA_spots`,
`THRESH$min_prop_nonzero` and `THRESH$min_max_lymphocyte`. All three are
**lower** bounds, so a section with a uniformly inflated Lymphocyte score
passes by construction.

Section 848 passed the automatic screen and was excluded on inspection of the
per-spot maps: its Lymphocyte score was elevated near-uniformly across the
capture area with no focal structure, against the sparse focal pattern of the
retained sections. The exclusion is applied by `02b_override-after-qc.R`, which
currently hardcodes the six retained IDs with no reason recorded — add a
comment there pointing to this section.

The control deconvolution against the full reference came afterwards and
confirmed the call. It is independent confirmation, **not** the selection
criterion, and Methods must state the order that way.
