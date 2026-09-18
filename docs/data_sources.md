# Data sources and provenance

All primary data in this study are secondary analyses of publicly deposited
datasets. This document records, for every input file, where it came from,
whether it is redistributable, and where a user of this repository should
obtain it.

Items marked **⚠ VERIFY** are unconfirmed and must be checked before the
repository is made public.

---

## 1. Primary spatial transcriptomics data

| Item | Value |
|---|---|
| Accession | GEO **GSE195661** |
| Contents | Visium spatial gene expression, 14 PFA ependymoma sections (11 primary, 3 matched recurrence) |
| Includes | Count matrices, `tissue_positions_list.csv`, `scalefactors_json.json`, hires/lores tissue images |
| Source publication | Donson et al., *Neuro-Oncology* 2023 (spatial transcriptomic analysis of epithelial/mesenchymal subpopulations in childhood ependymoma) |
| Consent | COMIRB #95-500 (University of Colorado) |
| Redistribute in this repo? | **No** — fetch from GEO |
| Redistribute on Zenodo? | **No** |

**Required repository asset:** `resources/gsm_section_map.tsv` — mapping from
GSM accession to the internal section IDs used throughout this codebase
(459, 459_2, 723, 723_2, 727, 812, 821, 848, 928, 928_2, 1101, 1239, 1269,
1513). Without this table no external user can connect any figure to GEO.
⚠ VERIFY: extract from the Visium manifest used in phase 00.

**Required script:** `scripts/00_fetch_and_arrange_geo.R` — downloads GEO
supplementary files and reshapes them into the
`neuro_onc_spatial_files/<section>/spatial/` layout the pipeline expects.
⚠ VERIFY: confirm all 14 sections have complete `spatial/` folders in GEO,
including both hires and lores images, and that image dimensions match those
asserted in `10_zone_and_spatial.R`.

---

## 2. Single-cell reference (deconvolution)

| Item | Value |
|---|---|
| Accession | GEO **GSE125969** |
| Contents | scRNA-seq of posterior fossa ependymoma; defines CEC / TEC / UEC / MEC neoplastic subpopulations |
| Processed mirror | `full_exprMatrix.tsv` from https://www.pneuroonccellatlas.org |
| Note | **scRNA-seq**, not snRNA-seq (resolved) |
| Cohort relationship | **Same patients as GSE195661, plus additional cases.** Deconvolution is patient-matched — state this in Methods (resolved) |
| Redistribute in this repo? | No (size) |
| Redistribute on Zenodo? | Derived products only — see §3 |

**RESOLVED:** the cell-type annotations used by
`03a_build_tumor_augmented_reference.R` (TEC-A/CEC, TEC-D, MEC, UEC-A,
UEC-A/profil, CEC, TEC-D low) are present in the GSE125969 cell metadata on
GEO. No labels originate with this study.

Consequence: the atlas download is a **convenience repackaging**, not a
dependency. Both counts and annotation are obtainable from GEO, so no
mirroring of atlas files is required and no redistribution permission is
needed for them.

Practical note: GEO likely serves this as per-sample matrices plus a
metadata table, while `03a` streams from a single concatenated
`full_exprMatrix.tsv`. Rewriting `03a` against the GEO layout is not
worthwhile — instead, document the equivalence and deposit the resulting
augmented reference (see §3), which is required regardless for the sampling
reason in §6.

**Required repository asset:** `resources/scrna_section_overlap.tsv` — which
GSE125969 samples correspond to which GSE195661 sections, and which scRNA
samples have no spatial counterpart. `03a` samples cells balanced by sample,
so this table documents what actually entered the augmented reference.

Framing consequence: reference-dependence is closed *for this cohort*.
Fit quality is not evidence of generalization to independent cohorts.

Licence / terms of use for pneuroonccellatlas.org downloads: no longer
blocking, since nothing atlas-derived is redistributed. Cite the atlas as a
resource, not as a data source.

---

## 3. Derived references (deposit on Zenodo)

| File | Origin | Why it must be deposited |
|---|---|---|
| `expMatrix_pseudobulk.tsv` (30,279 genes × 14 cell types) | Aggregated from GSE125969 / atlas matrix | Input to Semla NNLS; small; avoids redistributing the full cell-level matrix |
| `tumor_augmented_reference/` | `03a` — 250 cells sampled per neoplastic population, appended to the immune-focused reference | **Sampling step; see §6.** The exact reference used cannot otherwise be rebuilt |
| CIBERSORTx output fraction tables | CIBERSORTx web/Docker service | Not reproducible by script — see §4 |
| `SES_CIBERSORTx_EPN_pseudobulk.txt` | DeepTIL/SES, external | Second non-scriptable step; sole input to stages 05b and 05c |
| CIBERSORTx input files | `05_hotspot_pseudobulk_CIBERSORTx.R` — `pseudobulk_cpm_CIBERSORTx.tsv` (CPM, linear space), `pseudobulk_counts.tsv`, `pseudobulk_sample_metadata.tsv` | Allows an independent re-run by anyone with an account |

---

## 4. Third-party tools and matrices with restricted redistribution

| Item | Status | Handling |
|---|---|---|
| **LM7 signature matrix** | Tosolini et al., *OncoImmunology* 2017;6(3):e1284723 — supplementary Table S3. Built from n=265 purified-leucocyte Affymetrix HGU133 Plus 2.0 transcriptomes. Seven subsets: B, CD4 T, CD8 T, γδ T, NK, Mo-Ma-DC, granulocytes | **Do not commit.** Used unmodified as a CIBERSORTx input; no LM7-derived artifact is produced by this pipeline, so nothing needs redistributing. Cite the paper and link the supplementary |
| **CIBERSORTx** | Registration-gated web/Docker service | Not scriptable. Document run parameters in `docs/cibersortx.md`: B-mode batch correction, quantile normalisation off, relative mode, 500 permutations. Deposit inputs and outputs |
| **Zone gene programs** (`noac219_suppl_supplementary_data1.xlsx`) | Donson et al., *Neuro-Oncology* 2023, supplementary data 1 — the GSE195661 source publication. Used by stages 06 and 06b | ⚠ VERIFY the article licence. Do not commit unless CC-BY. Cite and link |
| **NicheNet** (`ligand_target_matrix`, `lr_network`, `weighted_networks`) | Downloaded at runtime by stage 08 | Do not commit. Pin the download URL and release version in `docs/pipeline.md` |
| **CellChatDB / CellPhoneDB** | Loaded in stage 08; extracted to `CellChatDB_pairs_raw.tsv` and `CellPhoneDB_pairs_raw.tsv` | ⚠ VERIFY each licence. These extractions redistribute database content — regenerate at runtime rather than committing |
| **DeepTIL / SES** | External tool converting CIBERSORTx fractions to abundance estimates | Not scriptable. ⚠ Record the tool version. Deposit `SES_CIBERSORTx_EPN_pseudobulk.txt` on Zenodo |
| **MSigDB** (Hallmark, C2:CP:REACTOME via `msigdbr`) | Registration / redistribution terms apply | Do not commit gene sets. Download at runtime in `11_gsea_characterization.R`; pin the `msigdbr` version in `renv.lock` |
| **CellChatDB / CellPhoneDB** | ⚠ VERIFY licence for each | Load from package at runtime rather than committing extracted tables |

---

## 5. Original to this study (commit to GitHub)

| File | Notes |
|---|---|
| `SupplementaryTableS2_ZoneGenes.xlsx` | Zone programs: Epithelial, Mesenchymal (306), Vascular (≤22), Myeloid (182); sheet S2d holds zone × lymphocyte Jaccard indices |
| Opportunity-score gene programs | γδ core, αβ core, inhibitory program, and the five γδ mechanism axes, with `GENE_A\|GENE_B` alias syntax |
| `section_qc.tsv` | All 14 sections, raw NNLS metrics, automatic-screen outcome |
| `qc_sensitivity_n40_n60_n80.tsv` | LM7 gene-set-size sensitivity |
| `11_gsea_*_<target>.tsv` | Per-target GSEA results and gene rankings |
| Cohort clinical table | Age at surgery, sex, WHO grade, methylation subtype, classifier score, 1q status. Derived from the GSE195661 source publication's supplementary material — cite, do not present as new |

⚠ **VERIFY:** treatment exposure of the matched pairs (which relapse sections
were irradiated). If this was supplied directly by a co-author rather than
taken from the published supplement, it is new clinical annotation on a
consented paediatric cohort and needs explicit sign-off before publication.

---

## 6. Stochastic steps requiring seeds

Public inputs only reproduce published numbers if every random step is
seeded. Record seeds in `docs/pipeline.md`.

| Script | Random step | Status |
|---|---|---|
| `03a_build_tumor_augmented_reference.R` | Samples 250 cells per neoplastic population | ⚠ **Highest priority.** Seed it, and deposit the reference actually used |
| `09b_scoring_qc.R` | 100-replicate split-half | ⚠ Seed |
| Zone-enrichment null | Rigid rotation/translation of hotspot mask | ⚠ Seed |
| `11_gsea_characterization.R` | `fgseaMultilevel` | ⚠ Seed |
| CIBERSORTx | 500 permutations | Not seedable — deposit outputs |

---

## 7. Two-tier reproducibility

**Tier 1 — full pipeline from public data.** Fully feasible: GSE195661 and
GSE125969 are both open, including the cell-type annotations. Runtime hours.
The only external dependency is the CIBERSORTx step, which must be run
manually or skipped using the deposited outputs.

**Tier 2 — figures from deposited intermediates.** Runtime minutes. This is
the tier reviewers will actually attempt and the one to test on a clean
machine before submission.

---

## 8. Deposition plan

| Record | Type | DOI timing |
|---|---|---|
| Zenodo data record | Manual deposit: derived references + CIBERSORTx inputs/outputs | Reserve DOI before submission; cite in Data Availability |
| Zenodo code record | Created automatically by the GitHub release webhook | Minted at release; cannot be pre-reserved |

Note: DOI pre-reservation is not available through the GitHub–Zenodo
integration, which is why the data record must be a manual deposit. Test the
webhook on sandbox.zenodo.org first — published Zenodo records cannot be
deleted.

Add `CITATION.cff` (or `.zenodo.json`) to the repository root so the code
record inherits proper author metadata and ORCIDs rather than GitHub
usernames.

---

## 9. Pre-publication hygiene

- **Strip institutional proxy URLs.** Grep the manuscript, scripts, and all
  repo files for `proxy.kib.ki.se` and replace with canonical URLs. Proxied
  links are gated and broken for every external reader.
- Grep for absolute local paths (`D:/Ped-CNS_KBH/...`) — all path handling
  belongs in `config.R`.
- Confirm no superseded scripts reach the public branch.
- Platform-mismatch note for Limitations: LM7 is microarray-derived and is
  applied here to Visium UMI counts.
