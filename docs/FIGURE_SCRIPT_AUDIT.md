# Figure-script cleanup audit

This audit accompanies the cleaned manuscript figure scripts.

## Changes applied across all figure scripts

- Every figure script now sources the repository-level `config.R`.
- Absolute local paths such as `D:/Ped-CNS_KBH/...` were removed.
- Analysis-stage inputs now resolve through `stage_dir()`, `RAW_ROOT`, `REF`,
  `GOOD_SECTIONS`, `ALL_SECTIONS`, and `MATCHED_PAIRS` where appropriate.
- Figure outputs now resolve through `figure_path()` under `EPN_FIGURE_ROOT`.
- H&E image paths are reconstructed from `RAW_ROOT/<section>/spatial/` rather
  than trusting absolute paths embedded in older intermediate tables.
- Line endings were normalized to LF.

## config.R changes required by the cleaned figure scripts

The supplied cleaned `config.R` adds:

- `FIGURE_INPUT_ROOT`, controlled by `EPN_FIGURE_INPUT_ROOT` and defaulting to
  `$EPN_DATA_ROOT/figure_inputs`;
- `figure_input_path()`;
- `SEEDS$neighborhood_permutation = 1L`, preserving the original Figure 3 / S1
  permutation seed;
- the previously discussed corrected `DATA_ROOT` setup error message;
- removal of the obsolete `he_overlays` stage key.

## Figure-specific changes

### Figure 1

Portable inputs are expected under:

```text
$EPN_FIGURE_INPUT_ROOT/figure1/
├── merged_abundance_filtered_allimmune_with_tumor_type.tsv
└── deeptil_results/
    └── SES_CIBERSORTx_*.txt
```

The script now fails with an informative message if these inputs are absent.

### Figure 2

Portable inputs are expected under:

```text
$EPN_FIGURE_INPUT_ROOT/figure2/
├── harmony_umap.coords.tsv.gz
└── meta.tsv   # optional; falls back to REF$scrna_meta
```

Spatial hotspot inputs resolve from the `hotspots` stage and DeepTIL results
resolve from the configured `deeptil` stage. The panel-A label was corrected
from `snRNA-seq UMAP` to `scRNA-seq UMAP`.

**Important unresolved manuscript issue:** panel 2b currently computes hotspot
percentages only for `GOOD_SECTIONS` and then assigns 0% to the other sections
when building the all-section bar chart. Because hotspot calls do not exist for
those other sections, those zeros are not measured values. In addition, the
hotspot definition is based on within-section quantiles, so hotspot percentage
is largely a property of the definition rather than a directly comparable
biological quantity. The cleaned script preserves the existing panel for
compatibility but adds an explicit warning comment. This panel should be
reconsidered before final preprint release.

### Figure 3 and Supplementary Figure S1

- Inputs resolve from the `zones` stage.
- H&E images resolve from `RAW_ROOT`.
- The original permutation seed of 1 is now centralized as
  `SEEDS$neighborhood_permutation`.
- S1 still reads Figure 3 permutation output when present and otherwise
  recomputes its required permutations, preserving existing behavior.

The permutation calculations still live in the figure scripts. They are now
deterministic, but moving them into a dedicated analysis/QC stage would make the
analysis/figure separation cleaner in a future refactor.

### Figure 4

- Functional-state inputs resolve from `functional_states`.
- Zone inputs resolve from `zones`.
- LR/NicheNet inputs resolve from `lr_interactions`.
- Raw Visium expression and H&E images resolve through `raw_vis` / `RAW_ROOT`.
- The obsolete comment referring to a separate `figure4_nicheNet.R` was removed;
  the current unified script already contains conditional panels g/h.
- The Vascular section-exclusion denominator was corrected to match stage 08
  (`Vascular = c("928")`), rather than treating all five LR sections as tested.

### Figure 5 and Supplementary Figure S3

- Stage 09 inputs resolve from `scoring`.
- Stage 10 inputs resolve from `zone_spatial`.
- Matched primary/relapse samples are derived from `MATCHED_PAIRS`.
- The manuscript display order is retained but checked against `ALL_SECTIONS`.
- Figure 5 reconstructs H&E paths from `RAW_ROOT` instead of using old absolute
  paths stored in `10_image_lookup.tsv`.

### Figure 6

- GSEA inputs resolve from the `gsea` stage.
- Outputs resolve to `FIGURE_ROOT/figure6`.
- The stale manuscript-local output-path comment was removed.

### Supplementary Figure S2

- Cross-section QC inputs resolve from `lr_qc`.
- top-pair inputs resolve from `lr_interactions/top_pairs`.
- outputs resolve to `FIGURE_ROOT/S2`.

## Validation performed

All cleaned scripts were checked for:

- remaining absolute Windows paths;
- stale `Semla integration` path rewrites;
- remaining references to removed path variables;
- balanced parentheses, brackets, braces and quoted strings.

An R interpreter is not installed in the current execution environment, so the
scripts could not be executed or parsed with `Rscript` here. They should still
be run once locally against the final data layout before the repository is made
public.
