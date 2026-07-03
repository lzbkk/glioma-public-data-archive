# Archive Manifest

Date: 2026-07-03

Repository URL: https://github.com/lzbkk/glioma-public-data-archive

## Scope

This public archive contains analysis scripts, derived source-data packages,
journal-facing supplementary tables and metadata for a public multi-omic glioma
study centered on LAP3 as a malignant-microenvironmental state anchor.

The archive does not redistribute controlled-access raw data or third-party raw
omics files.

## Top-Level Contents

| Path | Contents |
|---|---|
| `README.md` | Repository overview, third-party data access routes and rebuild notes |
| `archive/` | Derived source-data and supplementary-table packages for manuscript support |
| `metadata/` | Data/code availability draft, supplementary manifests and reference metadata |
| `scripts/` | Analysis and package-building scripts organized by data layer |
| `LICENSE_NOTICE.md` | Licence boundary and pending author confirmation |
| `.gitignore` | Guardrails against raw data, caches and credentials |

## Main Archived Files

| File | SHA256 |
|---|---|
| `archive/Source_Data_Package_20260702.tar.gz` | `b4311b25255167e675efc9bd03ef83dd3b23af6d3bce8c58b524a074518e4a8d` |
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

## Remaining Before DOI Freeze

- Confirm final licence choices.
- Create a GitHub release.
- Create or link a Zenodo record and obtain DOI/record URL.
- Backfill Zenodo DOI, release version and licence into manuscript
  Data/Code Availability.
