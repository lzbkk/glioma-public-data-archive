# Public Archive Upload Checklist

Date: 2026-07-02

Decision: GitHub repository + Zenodo archived release DOI

## 1. Pre-upload Checks

```text
[ ] Confirm final manuscript title or working archive title.
[ ] Confirm creator names, ORCID IDs and affiliations.
[ ] Confirm licence: recommended CC BY 4.0 for generated derived data; recommended MIT for code, unless author team chooses otherwise.
[x] Confirm GitHub repository visibility and URL: https://github.com/lzbkk/glioma-public-data-archive
[x] Confirm whether code and data use one Zenodo release or separate code/data records: default is one Zenodo release for generated source-data/supplement package plus GitHub release/archive for code, unless author team chooses split records.
[x] Recheck GBM-Space journal-version status on 2026-07-02: no peer-reviewed article found; repeat on upload day.
[x] Recheck CGGA portal-preferred citation on 2026-07-02: use Zhao et al. GPB 2021, DOI 10.1016/j.gpb.2020.10.005, PMID 33662628.
```

## 2. Files to Upload to Zenodo

```text
[ ] Project_Management/Submission_Package/Source_Data_Package_20260702.tar.gz
[ ] Project_Management/Submission_Package/Journal_Supplementary_Tables_20260702.tar.gz
[ ] Project_Management/Submission_Package/Journal_Supplementary_Tables_20260702/Supplementary_Tables_1_11_JournalFacing.xlsx
[ ] Project_Management/Submission_Package/PUBLIC_ARCHIVE_UPLOAD_CHECKLIST_20260702.md
[ ] Project_Management/Submission_Package/ZENODO_METADATA_DRAFT_20260702.md
[ ] Project_Management/Submission_Package/GITHUB_REPOSITORY_README_DRAFT_20260702.md
[ ] Project_Management/Audits/20260702/Public_Archive_Freeze_Audit_20260702.md
```

## 3. Expected Checksums

```text
b4311b25255167e675efc9bd03ef83dd3b23af6d3bce8c58b524a074518e4a8d  Source_Data_Package_20260702.tar.gz
061ad7e295f4fba4c046f4efd46abd68ba34f9bd13c119328773d49b89bd9b43  Journal_Supplementary_Tables_20260702.tar.gz
52f0b9e830cef7ef0f364a0acf11145973f9b5355ee4134f4736ae85c4a562e8  Supplementary_Tables_1_11_JournalFacing.xlsx
```

## 4. GitHub Repository Minimum Contents

```text
[ ] README with manuscript overview, third-party data access routes and reproduction boundaries.
[ ] Use `Project_Management/Submission_Package/GITHUB_REPOSITORY_README_DRAFT_20260702.md` as the starting README.
[ ] LICENSE for code.
[ ] Data_Bulk_TCGA/Data_Merged/scripts/
[ ] Data_scRNA_GEO/scripts/
[ ] Data_Spatial_Public/scripts/
[ ] Data_Protein_Public/scripts/
[ ] Data_Longitudinal_Public/GLASS/scripts/
[ ] Data_Perturbation_Public/LINCS_GSE70138/scripts/
[ ] Project_Management/Operations/build_submission_source_data_package.R
[ ] Project_Management/Operations/build_journal_supplementary_tables.R
[ ] Project_Management/Plans/Data_Code_Availability_Draft_20260702.md
[ ] Link to Zenodo DOI after release.
[ ] Add GitHub remote for local publishing if needed: `git@github.com:lzbkk/glioma-public-data-archive.git`.
```

Do not include:

```text
controlled-access raw data
third-party raw matrices that cannot be redistributed
local credentials, API keys, cookies or private `.Renviron`
large local caches unless explicitly intended for archive
```

## 5. Zenodo Metadata Fields

```text
Title: see ZENODO_METADATA_DRAFT_20260702.md
Creators: [TO ADD]
Description: see ZENODO_METADATA_DRAFT_20260702.md
Resource type: Dataset; if GitHub release is also archived through Zenodo, use Software and dataset or create a linked software record.
Version: v20260702
Licence: [TO CONFIRM; recommended CC BY 4.0 for generated derived data/tables and MIT for code]
Keywords: glioma, glioblastoma, LAP3, single-cell RNA-seq, spatial transcriptomics, public multi-omics
Related identifiers: GitHub URL `https://github.com/lzbkk/glioma-public-data-archive`, manuscript DOI if available, third-party data DOI/accessions
```

## 6. After DOI Creation

Update:

```text
[ ] Project_Management/Plans/Data_Code_Availability_Draft_20260702.md
[ ] Project_Management/Submission_Package/Source_Data_Package_20260702/README.md
[ ] Project_Management/Submission_Package/Journal_Supplementary_Tables_20260702/README.md
[ ] Project_Management/Audits/20260702/Public_Archive_Freeze_Audit_20260702.md
[ ] Project_Management/Core/项目执行看板.md
[ ] 实验记录本.md
[ ] Final manuscript Data Availability
[ ] Final manuscript Code Availability
```

Insert:

```text
GitHub URL: https://github.com/lzbkk/glioma-public-data-archive
Zenodo DOI: [TO ADD]
Zenodo record URL: [TO ADD]
Release version: [TO ADD]
Archive publication date: [TO ADD]
```

## 7. Current Verdict

```text
Ready for author-side repository creation and Zenodo upload preparation.
GitHub repository is created. Not DOI-frozen until Zenodo DOI/record URL, release version and licence are confirmed.
```
