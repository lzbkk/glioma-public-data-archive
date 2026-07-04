# LAP3 Glioma State Public Multi-Omics Study

Status: public repository for manuscript source-data and analysis-code support.

Repository URL: https://github.com/lzbkk/glioma-public-data-archive

Git SSH remote for maintainers:

```text
git@github.com:lzbkk/glioma-public-data-archive.git
```

## Overview

This repository contains analysis scripts and reproducibility support files
for a public-data study of LAP3 as a benchmark-supported anchor of a
malignant-microenvironmental state in glioma. The study integrates bulk glioma
cohorts, Core GBmap single-cell profiles, GBM-Space Visium spatial
transcriptomics, CPTAC proteogenomic/metabolomic data and GLASS paired
longitudinal tumors.

The repository is intended to support the manuscript's derived source data,
figures, supplementary tables and evidence-boundary analyses. It does not
redistribute controlled-access or third-party raw data.

## Manuscript Artifacts

Generated source-data and supplementary-table packages are included here under
`archive/` and are archived through the linked Zenodo dataset record. Analysis
scripts are kept in this GitHub repository rather than in a separate Zenodo code
tarball:

```text
Zenodo version DOI: 10.5281/zenodo.21161701
Zenodo concept DOI: 10.5281/zenodo.21161700
Release version: v20260703
```

Current archive files:

```text
archive/Source_Data_Package_20260702.tar.gz
archive/Journal_Supplementary_Tables_20260702.tar.gz
archive/Supplementary_Tables_1_11_JournalFacing.xlsx
```

## Repository Contents

Public repository contents:

```text
archive/
metadata/
scripts/Data_Bulk_TCGA/Data_Merged/scripts/
scripts/Data_scRNA_GEO/scripts/
scripts/Data_Spatial_Public/scripts/
scripts/Data_Protein_Public/scripts/
scripts/Data_Longitudinal_Public/GLASS/scripts/
scripts/Data_Perturbation_Public/LINCS_GSE70138/scripts/
scripts/Project_Management/Operations/
ARCHIVE_MANIFEST.md
CHECKSUMS.sha256
LICENSE_NOTICE.md
```

Do not include:

```text
controlled-access raw data
third-party raw matrices unless redistribution is explicitly allowed
GBM-Space raw H5AD/tarball files
GLASS raw Synapse release files
large local caches
tmux logs or private operational logs
API keys, tokens, cookies, `.Renviron`, proxy settings or credentials
```

## Third-Party Data Access

Readers should obtain third-party inputs from their original sources:

```text
TCGA GBM: DOI 10.1038/nature07385
TCGA LGG: DOI 10.1056/NEJMoa1402121
CGGA: DOI 10.1016/j.gpb.2020.10.005; PMID 33662628
Core GBmap: DOI 10.1093/neuonc/noaf113
GBM-Space preprint: DOI 10.1101/2025.05.13.653495
GBM-Space EGA dataset: EGAD00001015527
GBM-Space EGA study: EGAS00001005801
GBM-Space DAC: EGAC00001000205
EGA archive paper: DOI 10.1093/nar/gkab1059
CPTAC GBM: DOI 10.1016/j.ccell.2021.01.006
GLASS: DOI 10.1038/s41586-019-1775-1; Synapse resources as listed in the manuscript
```

## Rebuild Commands

Run from the repository root after required third-party inputs are present:

```bash
Rscript --vanilla Project_Management/Operations/build_submission_source_data_package.R
Rscript --vanilla Project_Management/Operations/build_journal_supplementary_tables.R
```

In this public archive, those scripts are available under
`scripts/Project_Management/Operations/`. The original project-relative paths are
preserved in the scripts and source-data manifests; users should adapt paths
after downloading the required third-party inputs.

Individual analysis scripts are organized by data layer. Some scripts require
large third-party inputs that are not redistributed in this repository.

## Licences

```text
Generated derived source data and supplementary tables: CC BY 4.0
```

This licence applies only to this study's generated derived outputs and
documentation in the Zenodo dataset record. It does not change access terms for
third-party datasets. Code licensing for the associated GitHub repository is
handled separately from the Zenodo dataset record.

## Controlled-Access Boundary

The public source-data package does not include GBM-Space raw Visium files,
complete expression matrices, H5AD/tarball files or controlled-access
spot-level spatial coordinate/score tables. The Supplementary Fig. S4A public
source-data file is a section-level provenance/summary placeholder; users with
appropriate EGA access can regenerate the full panel source table from the
analysis scripts and controlled source data.

## Citation

Use the Zenodo version DOI for the source-data and supplementary-table package:

```text
10.5281/zenodo.21161701
```

For the associated analysis scripts, cite the GitHub repository:

```text
https://github.com/lzbkk/glioma-public-data-archive
```
