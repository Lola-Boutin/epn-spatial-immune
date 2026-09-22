# Pipeline: run order and dependencies

This document describes the current analysis workflow for the pediatric ependymoma
spatial immune profiling study.

All analysis scripts are under `scripts/analysis/` and source the repository-level
`config.R`:

```r
source(here::here("config.R"))
```

`config.R` is the single source of truth for data roots, stage directories,
analysis thresholds, cohort definitions and random seeds.

---

## 1. Required local configuration

Create a `.Renviron` file in the repository root, or otherwise define the
following environment variables before running the analysis:

```text
EPN_DATA_ROOT=/absolute/path/to/data
EPN_FIGURE_ROOT=/absolute/path/to/figures
```

Optional:

```text
EPN_LM7_PATH=/absolute/path/to/LM7_matrix.txt
```

`EPN_LM7_PATH` is used only by stage 05 to report LM7 signature-gene coverage in
the pseudobulk mixtures. The LM7 matrix itself is not distributed with this
repository.

The expected spatial-data layout under `EPN_DATA_ROOT` is:

```text
neuro_onc_spatial_files/
├── 459/
│   ├── filtered_feature_bc_matrix.h5
│   └── spatial/
│       ├── tissue_positions_list.csv
│       ├── scalefactors_json.json
│       └── tissue_hires_image.png
├── 459_2/
├── 723/
└── ...
```

Stage 00 derives `section_id` directly from these directory names.

---

## 2. Analysis run order

| Stage | Script | Main purpose | Approx. runtime | Stochastic? |
|---|---|---|---:|---|
| 00 | `scripts/analysis/00_build_manifest_and_raw_sections.R` | Build the Visium manifest, raw per-section Seurat objects, spot-coordinate tables and raw H&E QC overlays | ~30 min | No |
| 01 | `scripts/analysis/01_build_raw_barcode_coordinate_maps.R` | Build barcode/coordinate maps used to align spatial outputs | — | No |
| 02 | `scripts/analysis/02_run_semla_nnls_qc_and_save_lymph_only.R` | Run Semla NNLS with the full scRNA reference and retain the lymphocyte signal; perform automatic section QC | ~2 h | No |
| 02b | `scripts/analysis/02b_override_after_qc.R` | Apply the documented manual section-selection override after QC review | — | No |
| 02c | `scripts/analysis/02c_add_hotspots_to_coordinate.R` | Define spatial lymphocyte hotspots from the retained sections | — | No |
| 03a | `scripts/analysis/03a_build_tumor_augmented_reference.R` | Build the PFA tumor-augmented reference used for NNLS sensitivity analysis | ~45 min | Yes |
| 03b | `scripts/analysis/03b_run_semla_nnls_from_augmented.R` | Re-run Semla NNLS with the tumor-augmented reference | — | No |
| 03c | `scripts/analysis/03c_compare_original_vs_tumor_augmented_hotspots.R` | Compare original and tumor-augmented lymphocyte scores/hotspots | — | No |
| 05 | `scripts/analysis/05_hotspot_pseudobulk_CIBERSORTx.R` | Build hotspot/background Visium pseudobulks for external CIBERSORTx deconvolution | ~2 min | No |
| 06 | `scripts/analysis/06_zones_from_corrected_backbone.R` | Assign spatial zones. Run twice: retained sections and all 14 sections | — | No |
| 06b | `scripts/analysis/06b_qc_signature_independence.R` | Test overlap between zone signatures and the lymphocyte NNLS marker signal | ~20 min | No |
| 07 | `scripts/analysis/07_functional_states.R` | Score lymphocyte functional states and retention-related programs | — | No |
| 08 | `scripts/analysis/08_lr_interactions.R` | Ligand–receptor analysis using CellChatDB, CellPhoneDB and optional NicheNet | ~3 h | Yes |
| 08b | `scripts/analysis/08b_qc_after_lr_interactions.R` | Cross-section QC of ligand–receptor results | ~2 min | No |
| 09 | `scripts/analysis/09_opportunity_map.R` | Build the immunotherapy opportunity scores and mechanism axes | ~40 min | No |
| 09b | `scripts/analysis/09b_scoring_qc.R` | Methodological QC and sensitivity analyses for the opportunity scores | ~25 min | Yes |
| 10 | `scripts/analysis/10_zone_and_spatial.R` | Zone-level and spatial characterization of opportunity scores | ~20 min | No |
| 11 | `scripts/analysis/11_gsea_characterization.R` | Pathway characterization of opportunity-score-associated genes | ~1 h | Yes |
| 11b | `scripts/analysis/11b_ranking_plots.R` | Supplementary ranking, concordance and program-gene plots | ~5 min | No |

A missing stage number is intentional: the former per-spot LM7/UCell stage was
removed after the workflow was replaced by the CIBERSORTx pseudobulk route.

---

## 3. Stage 06 must be run twice

`06_zones_from_corrected_backbone.R` contains:

```r
ALL_SECTIONS_RUN <- FALSE
```

Run the script once with:

```r
ALL_SECTIONS_RUN <- FALSE
```

to generate the six retained-section zone tables under the `zones` stage. These
are consumed by stages 07 and 08.

Then run it again with:

```r
ALL_SECTIONS_RUN <- TRUE
```

to generate the 14-section zone tables under `zones_all`. These are consumed by
stages 09 and 10.

Hotspot calls exist only for the retained sections. In the all-section run, the
other sections carry no hotspot context; stages 09 and 10 use the zone calls
rather than hotspot status.

---

## 4. External/manual steps

### CIBERSORTx

Stage 05 writes:

```text
pseudobulk_counts.tsv
pseudobulk_cpm_CIBERSORTx.tsv
pseudobulk_sample_metadata.tsv
```

The CPM mixture is run through CIBERSORTx with the LM7 signature matrix using
the parameters recorded in `config.R`:

```text
Batch correction: B-mode
Quantile normalisation: off
Mode: relative
Permutations: 500
```

The LM7 signature matrix is obtained from Tosolini et al. (OncoImmunology,
2017) and is not distributed in this repository.

### DeepTIL / SES

The CIBERSORTx fractions are subsequently processed externally with DeepTIL/SES.
The manuscript Figure 2 script expects:

```text
SES_CIBERSORTx_EPN_pseudobulk.txt
```

under the configured DeepTIL results directory.

This file is used for Figure 2d and 2e: lymphoid composition of hotspot
pseudobulks and paired hotspot-versus-background lymphoid abundance.

The manuscript analysis used
   `AutoCompare_SES_windows_011.pl` (calling `SES-Fred-007.r`), followed by
   `deeptil-004.jl`. Exact script checksums and provenance are recorded in
   `docs/data_sources.md`.

---

## 5. Random seeds

Seeds are centralized in `config.R`.

| Analysis | Seed |
|---|---:|
| Tumor-augmented reference sampling (`03a`) | 1 |
| Ligand–receptor permutation analysis (`08`) | 20260801 |
| Split-half scoring QC (`09b`) | 20260812 |
| GSEA (`11`) | 20260812 |

Do not change these values once the corresponding results are finalized.

Stage 08 uses 1000 hotspot-label permutations with the fixed
`lr_permutation` seed. The final seeded run was validated against the earlier
unseeded analysis: observed ligand–receptor scores were unchanged, and the
downstream LR shortlist is based on observed scores rather than permutation
p-values.

---

## 6. Ligand–receptor resources

Stage 08 uses:

- CellChatDB human interactions, downloaded by the script.
- CellPhoneDB interactions accessed through `liana`.
- NicheNet resources supplied as local RDS files when NicheNet analysis is
  enabled.

The NicheNet files expected by the current script are:

```text
ligand_target_matrix_nsga2r_final.rds
lr_network_human_21122021.rds
weighted_networks_nsga2r_final.rds
```

These third-party resource files are not part of this repository.

---

## 7. Figure generation

Figure scripts are stored under:

```text
scripts/figures/
```

The current repository contains:

```text
figure1.R
figure2.R
figure3.R
figure4.R
figure5.R
figure6.R
S1.R
S2.R
S3.R
```

The analysis pipeline and figure-generation scripts are intentionally separated:
analysis scripts generate reusable tables/RDS objects, while figure scripts turn
those outputs into manuscript panels.

### Figure 2

`figure2.R` currently generates:

- Figure 2a: scRNA-seq UMAP.
- Figure 2b: percentage of spots with a non-zero lymphocyte NNLS coefficient
  across all 14 Visium sections.
- Figure 2c: spatial lymphocyte hotspot maps for the retained sections.
- Figure 2d: DeepTIL-inferred hotspot lymphoid composition.
- Figure 2e: paired DeepTIL lymphoid abundance in hotspot versus background
  pseudobulks.

Figure-specific inputs resolve through `EPN_FIGURE_INPUT_ROOT`, while analysis
outputs resolve through `config.R` and the configured stage directories.

---

## 8. Section selection

Stage 02 applies the automatic lower-bound QC criteria defined in `config.R`:

```text
min_nonNA_spots
min_prop_nonzero
min_max_lymphocyte
```

Because all three are lower bounds, they cannot identify a section with a
uniformly elevated lymphocyte score. Stage 02b therefore records the subsequent
manual inspection/override used to define the retained section set.

The final retained set in `config.R` is:

```text
459, 723, 812, 821, 928, 1239
```

The reason for every manual exclusion should be described in the manuscript
Methods and/or stage-02b comments rather than being inferred from the code alone.

---

## 9. Remaining release tasks

Before public release:

1. Deposit the derived reproducibility artifacts described in
   `docs/data_sources.md`.
2. Ensure the final software environment is documented in
   `docs/software_environment.md`.
3. Run the cleaned analysis and figure scripts once from the final repository
   layout to confirm that all configured paths and dependencies resolve.
4. Verify that any clinical annotations not originating from the documented
   public datasets are appropriately described and handled.

