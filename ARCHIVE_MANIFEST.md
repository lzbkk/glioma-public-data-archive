# Archive Manifest

Date: 2026-07-03

Repository URL: https://github.com/lzbkk/glioma-public-data-archive

## Scope

This public archive contains analysis scripts, derived source-data packages,
journal-facing supplementary tables and metadata for a public multi-omic glioma
study centered on LAP3 as a malignant-microenvironmental state anchor.

The archive does not redistribute controlled-access raw data, third-party raw
omics files or controlled-access spot-level spatial coordinate/score tables.

## Top-Level Contents

| Path | Contents |
|---|---|
| `README.md` | Repository overview, third-party data access routes and rebuild notes |
| `archive/` | Derived source-data and supplementary-table packages for manuscript support |
| `metadata/` | Data/code availability notes, supplementary manifests and reference metadata |
| `scripts/` | Analysis and package-building scripts organized by data layer |
| `LICENSE_NOTICE.md` | Licence and third-party data-use boundary |
| `.gitignore` | Guardrails against raw data, caches and credentials |

## Main Archived Files

| File | SHA256 |
|---|---|
| `archive/Source_Data_Package_20260702.tar.gz` | `4040417bfc65863b3698d9d7ca5bc2884c7d90ab92f24370f1a98a9b0e57021c` |
| `archive/Journal_Supplementary_Tables_20260702.tar.gz` | `061ad7e295f4fba4c046f4efd46abd68ba34f9bd13c119328773d49b89bd9b43` |
| `archive/Supplementary_Tables_1_11_JournalFacing.xlsx` | `52f0b9e830cef7ef0f364a0acf11145973f9b5355ee4134f4736ae85c4a562e8` |

## Do Not Upload / Do Not Redistribute

The following file types or data classes should not be added to this repository:

- GBM-Space raw H5AD/tarball files.
- GLASS raw Synapse files.
- Raw TCGA, CGGA, Core GBmap or CPTAC matrices unless redistribution is
  explicitly permitted by the original source.
- Large local caches and partial downloads.
- `.Renviron`, API keys, tokens, cookies, proxy settings or private logs.

## Citation and Version

Archive version: `v20260703`

Zenodo version DOI for this dataset record:

```text
10.5281/zenodo.21161701
```

Zenodo concept DOI for all versions:

```text
10.5281/zenodo.21161700
```

Associated public code repository:

```text
https://github.com/lzbkk/glioma-public-data-archive
```

This archive record is intended for the generated source-data and
supplementary-table artifacts. Analysis scripts are provided through the
associated GitHub repository rather than as a separate code tarball in this
dataset record.
