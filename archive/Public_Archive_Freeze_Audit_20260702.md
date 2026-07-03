# Public Archive Freeze Audit 20260702

Date: 2026-07-02

## Routing

```text
Task: GitHub + Zenodo DOI / licence / release-version freeze audit
Skill: nature-data
Target journal: Neuro-Oncology Advances
Runtime: foreground documentation and package audit; no tmux required
Status: route frozen; GitHub repository created; Zenodo DOI/record URL pending author-side upload
```

## Official Sources Rechecked

| Source | Relevant rule captured | Use in this project |
|---|---|---|
| GitHub Docs, referencing and citing content | Zenodo can archive a public GitHub repository and issue a DOI; Zenodo issues a new DOI for each GitHub release; GitHub recommends including a licence | Supports GitHub + Zenodo route and confirms repository must be public before Zenodo GitHub archiving |
| GitHub Docs, licensing a repository | Without a licence, default copyright law applies; GitHub strongly encourages an open-source licence for open-source projects | Supports adding MIT licence for code before public repository release |
| Zenodo Docs, linking GitHub account | GitHub linking enables automatic repository archiving through Zenodo GitHub integration | Supports author-side GitHub-Zenodo workflow |
| Zenodo Docs, describe records | Deposit metadata includes DOI, resource type, title, creators, description, licences/rights, keywords and funding fields | Supports current Zenodo metadata draft fields |
| Zenodo Docs, manage records/versions | Published metadata can be edited; files generally require a new version or support contact after publication; versions have separate persistent identifiers linked across versions | Supports freezing file set and checksums before DOI publication |

Source URLs:

```text
https://docs.github.com/en/repositories/archiving-a-github-repository/referencing-and-citing-content
https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/licensing-a-repository
https://help.zenodo.org/docs/profile/linking-accounts/
https://help.zenodo.org/docs/deposit/describe-records/
https://help.zenodo.org/docs/deposit/manage-records/
https://help.zenodo.org/docs/deposit/manage-versions/
```

## Executive Verdict

```text
GO: Freeze the public archive route as GitHub repository + Zenodo archived release DOI.

GO: Use one primary Zenodo record for generated source-data/supplementary packages,
    linked to a public GitHub release for code. If the author team prefers strict
    separation, create a linked software record for the GitHub release later.

GO WITH AUTHOR CONFIRMATION: Recommended licences are CC BY 4.0 for generated
    derived data/tables and MIT for code.

NO-GO: Do not publish DOI or Data/Code Availability as final until the Zenodo DOI,
    Zenodo record URL, release version, licence and creator metadata are actually
    created and author-confirmed. The GitHub repository URL is fixed as
    https://github.com/lzbkk/glioma-public-data-archive.
```

## Package Files Verified

| File | Size | SHA256 | Archive status |
|---|---:|---|---|
| `Project_Management/Submission_Package/Source_Data_Package_20260702.tar.gz` | 2.3M | `b4311b25255167e675efc9bd03ef83dd3b23af6d3bce8c58b524a074518e4a8d` | ready after GitHub URL backfill rebuild |
| `Project_Management/Submission_Package/Journal_Supplementary_Tables_20260702.tar.gz` | 8.6M | `061ad7e295f4fba4c046f4efd46abd68ba34f9bd13c119328773d49b89bd9b43` | ready after GitHub URL backfill rebuild |
| `Project_Management/Submission_Package/Journal_Supplementary_Tables_20260702/Supplementary_Tables_1_11_JournalFacing.xlsx` | 4.6M | `52f0b9e830cef7ef0f364a0acf11145973f9b5355ee4134f4736ae85c4a562e8` | ready after GitHub URL backfill rebuild |
| `Project_Management/Submission_Package/PUBLIC_ARCHIVE_UPLOAD_CHECKLIST_20260702.md` | text | not checksum-critical | include in record |
| `Project_Management/Submission_Package/ZENODO_METADATA_DRAFT_20260702.md` | text | not checksum-critical | include in record |
| `Project_Management/Submission_Package/GITHUB_REPOSITORY_README_DRAFT_20260702.md` | text | not checksum-critical | use as GitHub README starting point and include in record if useful |
| `Project_Management/Audits/20260702/Public_Archive_Freeze_Audit_20260702.md` | text | not checksum-critical | include in record or repository governance folder |

Note: checksums above were verified after rebuilding the source-data package,
journal-facing supplementary tables and `.tar.gz` archives on 2026-07-03
following GitHub repository URL backfill.
If packages are rebuilt again after this audit, regenerate this section.

## Upload / Do-Not-Upload Classification

| Category | Route | Upload to Zenodo/GitHub? | Notes |
|---|---|---|---|
| Figure 1-5 panel source data and Supplementary Fig. S4A source data | generated derived source data | Yes | In `Source_Data_Package_20260702.tar.gz` |
| Supporting figure key-results, panel maps and export-QC tables | generated derived support tables | Yes | In source-data package |
| Supplementary Tables 1-11 journal-facing workbook | generated derived supplement | Yes | Upload combined workbook and/or table tarball |
| Analysis scripts for Figure 1-5 and supporting analyses | generated code | Yes, via public GitHub repository | Add MIT licence and README; do not include raw data |
| TCGA/CGGA source expression/clinical matrices | reused third-party public sources | No, cite/access through original sources | Share only derived summaries/source data |
| Core GBmap H5AD | reused public/third-party source | No, cite/access through original resource | Share donor-state summaries only |
| GBM-Space raw Visium H5AD/tarball | controlled-access EGA / third-party restricted | No | EGA dataset `EGAD00001015527`, study `EGAS00001005801`, DAC `EGAC00001000205` |
| GLASS raw Synapse files | third-party restricted/reused source | No | Readers use Synapse access route |
| CPTAC publication Table S2/raw workbook | reused public source | No by default | Share derived tables and cite Wang et al. |
| Local caches, temporary downloads, tmux logs, `.Renviron`, tokens, credentials | local/private operational data | No | Explicitly excluded |

## GitHub Repository Freeze Plan

Recommended repository scope:

```text
README.md
LICENSE
CITATION.cff or citation note after Zenodo DOI exists
Data_Bulk_TCGA/Data_Merged/scripts/
Data_scRNA_GEO/scripts/
Data_Spatial_Public/scripts/
Data_Protein_Public/scripts/
Data_Longitudinal_Public/GLASS/scripts/
Data_Perturbation_Public/LINCS_GSE70138/scripts/
Project_Management/Operations/build_submission_source_data_package.R
Project_Management/Operations/build_journal_supplementary_tables.R
Project_Management/Plans/Data_Code_Availability_Draft_20260702.md
Project_Management/References/NOA_Reference_Final_Pass_20260702.md
Project_Management/Submission_Package/PUBLIC_ARCHIVE_UPLOAD_CHECKLIST_20260702.md
Project_Management/Submission_Package/ZENODO_METADATA_DRAFT_20260702.md
```

Do not include:

```text
controlled-access raw data
third-party raw matrices unless redistribution is explicitly allowed
large H5AD/tarball/cache objects
private logs
API keys, cookies, tokens, `.Renviron`, proxy settings
personal author contact details not intended for publication
```

Recommended release tag:

```text
v20260702-noa-pre-submission
```

If the author team wants a more conservative name before actual submission:

```text
v20260702-preprint-review
```

## Zenodo Record Freeze Plan

Recommended resource type:

```text
Dataset
```

Rationale: the primary Zenodo files are derived source data and supplementary
tables. The code should be browsable through GitHub and DOI-archived through a
GitHub release. If the author team wants code and data in one Zenodo record,
use `Software and dataset` or explain the mixed record in the description.

Recommended title:

```text
Source data, supplementary tables and analysis code for a LAP3-centered malignant-microenvironment state study in glioma
```

Recommended version:

```text
v20260702
```

Recommended licences pending author confirmation:

```text
Generated derived data and tables: CC BY 4.0
Analysis code: MIT
```

Boundary:

```text
These licences apply only to this study's generated derived outputs and code.
They do not alter third-party source-data access terms for TCGA, CGGA, Core
GBmap, GBM-Space, CPTAC or GLASS.
```

## Data Availability Draft for Final Backfill

Use after Zenodo DOI is created:

```text
Processed source data generated by this study for Figures 1-5 and Supplementary Fig. S4A, together with journal-facing supplementary tables and archive metadata, are available in Zenodo at [Zenodo DOI / record URL], version [release version]. Analysis scripts are available in a public GitHub repository at https://github.com/lzbkk/glioma-public-data-archive and are archived with the corresponding release [GitHub release URL / Zenodo software DOI if separate]. The deposited files contain derived source data, panel maps, key-results tables, export-quality-control tables, Supplementary Tables 1-11 and metadata/checksum records. The archive does not redistribute controlled-access or third-party raw data.
```

Keep the third-party source paragraph from
`Project_Management/Plans/Data_Code_Availability_Draft_20260702.md`, including:

```text
GBM-Space: EGA dataset EGAD00001015527, study EGAS00001005801, DAC EGAC00001000205; stable preprint DOI 10.1101/2025.05.13.653495; upload-day journal-version recheck required.
CGGA: Zhao et al. GPB 2021, DOI 10.1016/j.gpb.2020.10.005, PMID 33662628.
```

## FAIR / Metadata Audit

| FAIR item | Status | Action |
|---|---|---|
| Persistent identifier | partial | GitHub repository created; create Zenodo DOI and GitHub release URL |
| File list and package contents | ready | Source/supplement package manifests exist |
| Checksums | ready, but regenerate after rebuild | Current three main artifact checksums recorded |
| README/data dictionary | partial | Source and supplement README exist; GitHub README still needed |
| Provenance | strong | Source-data index, panel maps, key-results and build scripts exist |
| Licence | pending author confirmation | Recommended CC BY 4.0 for generated derived data/tables; MIT for code |
| Third-party restrictions | explicit | GBM-Space/EGA and GLASS/Synapse routes retained; raw redistribution blocked |
| Versioning | proposed | `v20260702` for Zenodo; `v20260702-noa-pre-submission` for GitHub release |
| Creator metadata | partial | Zibo Li and Weiqing Wan entered provisionally; additional authors/ORCID/departments pending |
| Reviewer accessibility | pending | Public release before submission is simplest; otherwise test private reviewer link |

## Blocking Items Before DOI Publication

```text
P0.1 Confirm final author/creator list and order.
P0.2 Confirm official English names, departments, affiliations and ORCID IDs.
P0.3 Confirm licence choices with author team/institution.
P0.4 Populate public GitHub repository and create release.
P0.5 Upload Zenodo record and obtain DOI/record URL.
P0.6 Backfill Data/Code Availability, Source Data Package README, Supplementary Tables README,
     NOA working manuscript, cover letter and project management files.
P0.7 Recheck GBM-Space journal-version status on upload day.
```

## 中文核对

```text
1. 路线可以冻结：GitHub + Zenodo DOI。
2. GitHub repository 已创建：`https://github.com/lzbkk/glioma-public-data-archive`。
3. 现在不能伪造 DOI，也不能把 Zenodo pending 字段写成已完成。
4. 能上传的是本研究生成的 derived source data、补表、脚本和元数据。
5. 不能上传 GBM-Space raw Visium、GLASS raw Synapse、H5AD/tarball/cache、原始第三方矩阵和任何凭据。
6. 许可证推荐：数据/补表 CC BY 4.0，代码 MIT；但必须作者/单位确认。
7. CGGA 引用已经解决；GBM-Space 仍需正式上传当天再查一次有没有 journal version。
```
