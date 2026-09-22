# Data sources and provenance

This repository contains analysis code for secondary analyses of publicly
available pediatric ependymoma datasets. Large primary datasets and
third-party resources are not redistributed here.

This document records the source and expected local representation of each
analysis input.

---

## 1. Spatial transcriptomics

| Item | Details |
|---|---|
| Accession | GEO **GSE195661** |
| Data type | 10x Genomics Visium spatial transcriptomics |
| Cohort used here | 14 PFA ependymoma sections |
| Source publication | Donson et al., *Neuro-Oncology* (2023) |
| Repository handling | Primary data are not committed; users obtain them from GEO |

The analysis expects one local directory per internal section ID under:

```text
$EPN_DATA_ROOT/neuro_onc_spatial_files/
```

For example:

```text
459/
459_2/
723/
723_2/
727/
812/
821/
848/
928/
928_2/
1101/
1239/
1269/
1513/
```

Each section directory must contain the Visium matrix plus the associated
`spatial/` files required by stage 00.

### Important mapping requirement

Stage 00 derives `section_id` directly from the directory name. Therefore an
external user needs an explicit mapping from GEO sample identifiers to the
internal section IDs above.

Before public release, provide either:

- a small mapping table in the repository, or
- a download/arrangement script that creates the expected directory names from
  the GEO files.

The mapping should contain no unnecessary clinical information.

---

## 2. Single-cell RNA-seq reference

| Item | Details |
|---|---|
| Accession | GEO **GSE125969** |
| Data type | scRNA-seq |
| Role | Reference for Semla NNLS and source material for the tumor-augmented reference |
| Repository handling | Cell-level expression matrices are not committed |

The current pipeline expects the following files under `$EPN_DATA_ROOT`:

```text
exprMatrix.tsv
meta.tsv
exprMatrix_tumor.tsv
meta_tumor.tsv
```

These paths are defined in `config.R`.

### Use in stage 02

`02_run_semla_nnls_qc_and_save_lymph_only.R` loads the complete
`exprMatrix.tsv` and `meta.tsv`, creates a Seurat reference object, assigns
`cell_type` as the NNLS group, and runs Semla NNLS against that full reference.

The code does **not** subset the reference separately for each spatial patient.
The safest wording for the manuscript/repository is therefore that the
single-cell reference cohort overlaps with the spatial cohort, rather than
describing the implemented deconvolution as patient-matched per section.

### Use in stage 03a

`03a_build_tumor_augmented_reference.R` uses the tumor expression/metadata,
restricts cells to PFA1/PFA2, and samples up to 250 cells per tumor cluster
(seed = 1) before merging those tumor populations with the original reference.

The exact tumor-augmented reference used for the final results should be
deposited as a derived artifact because the sampling step is stochastic.

---

## 3. Figure 2 scRNA-seq UMAP

`figure2.R` currently expects:

```text
meta.tsv
harmony_umap.coords.tsv.gz
```

The metadata are used to group cells into the manuscript display categories.

The provenance of `harmony_umap.coords.tsv.gz` is not documented in the current
repository. Before public release, either:

1. document exactly how these coordinates were generated and provide the
   generating code, or
2. deposit the coordinate table as a derived artifact with clear provenance.

The figure script should also be refactored to remove hard-coded local paths.

---

## 4. Zone gene programs

The spatial-zone analysis uses:

```text
noac219_suppl_supplementary_data1.xlsx
```

from the supplementary material of the GSE195661 source publication
(Donson et al., *Neuro-Oncology*, 2023).

The file is referenced through `ref_file("zone_genes_xlsx")` in `config.R`.

It should not be copied into this repository unless its redistribution terms
explicitly allow that. Users should instead obtain the source supplementary file
from the original publication.

---

## 5. LM7 signature matrix

The LM7 immune signature matrix is from:

Tosolini et al., *OncoImmunology* (2017), 6(3):e1284723.

It is used with CIBERSORTx for seven immune subsets:

```text
B cells
CD4 T cells
CD8 T cells
gamma-delta T cells
NK cells
Mo-Ma-DC
granulocytes
```

The LM7 matrix is not redistributed in this repository.

Stage 05 optionally reads a local LM7 file through:

```text
EPN_LM7_PATH
```

only to report signature-gene coverage in the generated pseudobulk mixtures.

---

## 6. CIBERSORTx

Stage 05 generates the CIBERSORTx input files:

```text
pseudobulk_counts.tsv
pseudobulk_cpm_CIBERSORTx.tsv
pseudobulk_sample_metadata.tsv
```

The manuscript run uses the parameters recorded in `config.R`:

```text
LM7 signature matrix
B-mode batch correction
quantile normalisation off
relative mode
500 permutations
```

Because the CIBERSORTx run occurs outside the R pipeline, the exact input and
output files used for the manuscript should be deposited as derived artifacts.

---

## 7. DeepTIL / SES

Figure 2d and Figure 2e use the external DeepTIL/SES output:

```text
SES_CIBERSORTx_EPN_pseudobulk.txt
```

The figure script expects this table to contain the CIBERSORTx-derived abundance
estimates for the LM7 populations and uses the five lymphoid populations:

```text
B cells
CD4 T cells
CD8 T cells
gamma-delta T cells
NK cells
```

SES scores were generated using AutoCompare_SES_windows_011.pl
(calling SES-Fred-007.r), followed by deeptil-004.jl.
Exact script SHA-256 checksums:
AutoCompare_SES_windows_011.pl:
4f60a52cfee887830c934cce3921d63eed0862ec313861271360f0a51587a739

deeptil-004.jl:
d32771b5eecaf6fc08a0255958ce37f535b25eae84373347c79108a2d87c53ce

Because this is an external, non-scripted step, the final output table used for
the paper should be included in the derived-data deposit.

---

## 8. Ligand–receptor resources

Stage 08 uses three external ligand–receptor resources.

### CellChatDB

The current script downloads `CellChatDB.human.rda` from the CellChat GitHub
repository at runtime and extracts ligand–receptor interactions belonging to
secreted signaling, ECM–receptor and cell–cell contact classes.

Before public release, pin the exact CellChat release or commit used rather than
depending on a moving `main` branch.

### CellPhoneDB

CellPhoneDB interactions are obtained through the R package `liana`:

```r
liana::select_resource("CellPhoneDB")
```

The `liana` package version should be captured in the reproducible R environment.

### NicheNet

When NicheNet is enabled, stage 08 expects local copies of:

```text
ligand_target_matrix_nsga2r_final.rds
lr_network_human_21122021.rds
weighted_networks_nsga2r_final.rds
```

These files are not downloaded by the current script and are not distributed in
this repository.

Before public release, record the exact source URL/release/version used for each
resource.

---

## 9. MSigDB

Stage 11 obtains Hallmark and Reactome gene sets at runtime through `msigdbr`.

The gene sets themselves should not be copied into this repository. The
`msigdbr` package version should be captured in the final reproducible software
environment.

---

## 10. Derived artifacts recommended for deposition

A derived-data archive should contain the files necessary to reproduce the
manuscript figures without repeating slow or external steps.

At minimum, consider depositing:

- the exact tumor-augmented reference generated by stage 03a;
- CIBERSORTx input tables from stage 05;
- the exact CIBERSORTx output used downstream;
- `SES_CIBERSORTx_EPN_pseudobulk.txt`;
- the Figure 2 UMAP coordinate table if it is not regenerated by repository
  code;
- selected analysis tables/RDS objects required by the figure scripts.

Do not use the derived-data archive as a mirror of GEO primary data.

---

## 11. What should not be committed to GitHub

Do not commit:

- full GSE195661 Visium data;
- full GSE125969 cell-level expression matrices;
- large Seurat/RDS intermediates unless specifically needed for a small
  reproducibility example;
- the LM7 signature matrix;
- Donson et al. supplementary files unless redistribution is clearly permitted;
- NicheNet resource RDS files unless their redistribution terms explicitly allow
  it;
- institutional proxy URLs or absolute local filesystem paths.

The GitHub repository should contain code and lightweight study-generated
documentation, while primary data remain at their original repositories and
derived reproducibility artifacts are deposited separately.

---

## 12. Provenance items to finalize before release

Before making the repository public, confirm and document:

1. GEO sample → internal section-ID mapping for GSE195661.
2. Provenance/generation of `harmony_umap.coords.tsv.gz`.
3. Exact DeepTIL/SES version.
4. Exact NicheNet resource source/release.
5. Exact CellChatDB release/commit.
6. Final software/package versions for the R environment.
7. Whether any additional clinical annotations used in the manuscript came
   from published supplementary material or from collaborators, and ensure the
   appropriate publication/consent handling is documented.
