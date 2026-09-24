# Spatial immune profiling in pediatric ependymoma

Analysis code for:

> Truong P, Mirzazadeh R, Lundeberg J, Blomgren K, Boutin L.
> Spatial immune profiling reveals lymphocyte confinement to myeloid–mesenchymal niches and heterogeneous therapeutic T-cell opportunity in pediatric ependymoma


Correspondence: lola.boutin@ki.se

---

## What this is

A secondary analysis of publicly deposited data. This repository contains
analysis code, documentation and small reproducibility files. Primary data are
obtained from GEO and are not redistributed here. Derived reproducibility
artifacts required to reproduce selected analyses will be deposited separately
for the public release.

See [`docs/data_sources.md`](docs/data_sources.md) for the provenance and
expected local representation of each input.

| | |
|---|---|
| Spatial transcriptomics | GEO [GSE195661](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE195661) — 14 PFA sections |
| Single-cell reference | GEO [GSE125969](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125969) — overlapping pediatric ependymoma reference cohort |
| Derived artifacts | To be deposited for public release |
| R | 4.5.3 — see `docs/software_environment.md` |

---

## Setup

```bash
git clone https://github.com/Lola-Boutin/epn-spatial-immune.git
cd epn-spatial-immune
```

Create `.Renviron` in the repository root:

```
EPN_DATA_ROOT=/absolute/path/to/data
EPN_FIGURE_ROOT=/absolute/path/to/figures
```

---

## Running the analysis

Obtain the GSE195661 Visium data from GEO and arrange the section directories
under `$EPN_DATA_ROOT/neuro_onc_spatial_files/` using the accession-to-section
mapping provided in:

`docs/geo_section_mapping.tsv`

Then run the analysis pipeline beginning with
`scripts/analysis/00_build_manifest_and_raw_sections.R`.

Then run stages in the order given in [`docs/pipeline.md`](docs/pipeline.md).

**One stage cannot be automated.** Stage 05 writes a mixture file that must be
uploaded to [CIBERSORTx](https://cibersortx.stanford.edu/) and run manually
with the LM7 signature matrix (B-mode batch correction, quantile normalisation
off, relative mode, 500 permutations). Results are placed in the stage-05 results directory. Parameters are recorded
in `config.R`. For the public release, the corresponding derived output used by
the manuscript will be included among the deposited reproducibility artifacts.

---

## Layout

```text
config.R                 central paths, thresholds, cohort definitions and seeds
scripts/
├── analysis/            numbered analysis and QC stages
└── figures/             manuscript and supplementary figure scripts
docs/
├── pipeline.md          run order and stage dependencies
├── data_sources.md      input and resource provenance
├── geo_section_mapping.tsv
├── software_environment.md
└── FIGURE_SCRIPT_AUDIT.md

No script contains a literal filesystem path. Everything resolves through
`stage_dir()` and `stage_input()` in `config.R`.

---

## Third-party resources

The **LM7** signature matrix is not redistributed here. Obtain it from the
supplementary material of Tosolini et al., *OncoImmunology* 2017;6(3):e1284723.
It is used unmodified.

**MSigDB** gene sets (Hallmark, C2:CP:REACTOME) are obtained at runtime via
`msigdbr`. The manuscript analysis used `msigdbr` 26.1.0; software versions
are recorded in `docs/software_environment.md`.

---

## License

Code in this repository is distributed under the terms specified in `LICENSE`.

Third-party datasets and resources retain their original licenses and are not
redistributed here unless explicitly stated.

## Citing

Citation information for the manuscript and the archived code release will be
added when the preprint and repository archive are public.
