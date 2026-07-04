# Data and Code Availability Draft 20260702

Date: 2026-07-02

Status: repository route fixed; GitHub repository pushed; Zenodo draft created and validated with reserved DOI; public record URL to be inserted after Zenodo publication

## 2026-07-02 Public Archive Decision

```text
Decision: GitHub repository for analysis scripts + Zenodo archived release DOI for generated source-data and supplementary-table files.
Public archive staging repository: https://github.com/lzbkk/glioma-public-data-archive.
Zenodo draft deposition: https://zenodo.org/deposit/21161701.
Reserved Zenodo version DOI: 10.5281/zenodo.21161701.
Zenodo concept DOI: 10.5281/zenodo.21161700.
```

The generated source-data package and journal-facing supplementary tables are archived through the Zenodo route, with analysis scripts maintained in the public GitHub repository. The public GitHub repository has been created and pushed at `https://github.com/lzbkk/glioma-public-data-archive`. A Zenodo draft has been created at `https://zenodo.org/deposit/21161701` with reserved version DOI `10.5281/zenodo.21161701` and concept DOI `10.5281/zenodo.21161700`; the draft is currently unsubmitted and should be cited as a final public record only after publication. The public Zenodo package is restricted to generated source-data/supplementary-table artifacts and public-facing documentation, and excludes controlled-access raw data, third-party raw omics files, controlled-access GBM-Space spot-level coordinate/score tables and internal project-governance checklists.

## Routing

```text
Task: Data/Code Availability draft for the LAP3 glioma manuscript
Skill: nature-data
Target journal working assumption: Neuro-Oncology Advances
Note: this is a formatting/policy working assumption, not a locked submission target.
Runtime: foreground text/provenance task; no tmux required
Local inputs:
  Project_Management/Plans/Supplementary_Source_Data_Index_20260701.md
  Project_Management/Plans/supplementary_source_data_index_20260701.csv
  Project_Management/Plans/supplementary_table_manifest_20260701.csv
  Project_Management/References/Methods_Provenance_Fill_20260701.md
  Project_Management/References/Claim_to_Reference_Citation_Pass_20260701.md
  Project_Management/Plans/Manuscript_Methods_Draft_20260701_polished.md
```

## Journal Policy Check

Neuro-Oncology Advances author instructions, checked 2026-07-02, state that authors are strongly encouraged to make the deidentified data and software code supporting paper conclusions available, and that underlying deidentified data may be requested for journal inspection during peer review. Therefore, this draft prioritizes formal dataset accessions, especially EGA accessions for GBM-Space, and keeps the GBM-Space preprint status as a submission-time caveat if no peer-reviewed version is available.

## Data Availability

Ready-to-paste working draft:

```text
Data Availability

This study did not generate new primary human sequencing, proteomic, metabolomic or spatial transcriptomic data. All primary datasets analysed in this study were obtained from previously generated public or controlled-access resources under their original access conditions.

TCGA glioma bulk expression and clinical data were derived from TCGA GBM and TCGA lower-grade glioma resources. The TCGA GBM and LGG source studies are identified by DOI 10.1038/nature07385 and DOI 10.1056/NEJMoa1402121, respectively. CGGA validation data were obtained from the CGGA mRNAseq_693 and mRNAseq_325 RNA-seq cohorts; the CGGA portal-preferred citation was rechecked on 2026-07-02 and is the journal article by Zhao et al. in Genomics, Proteomics & Bioinformatics (DOI 10.1016/j.gpb.2020.10.005; PMID 33662628).

Single-cell analyses used the Core GBmap CELLxGENE resource, cited by DOI 10.1093/neuonc/noaf113. The local Core GBmap H5AD file used for analysis had SHA256 checksum 459444da84efc45565ac9ccb68c4873dac033c07abfe395bc814b3a33f60c56c. GBM-Space Visium spatial transcriptomic data are available through the European Genome-phenome Archive under dataset EGAD00001015527 and study EGAS00001005801, subject to the access conditions of the repository and associated Data Access Committee. The EGA Data Access Committee accession is EGAC00001000205. The EGA study page recommends citing study accession EGAS00001005801 together with the EGA archive paper (DOI 10.1093/nar/gkab1059). A 2026-07-02 live recheck found no peer-reviewed journal article for the accompanying GBM-Space preprint; the stable preprint DOI is 10.1101/2025.05.13.653495. The EGA dataset page notes two new samples added on 2026-05-19; this manuscript reports analyses of the local downloaded/cache version described in Methods.

CPTAC GBM proteogenomic and metabolomic analyses used the public CPTAC GBM study by Wang et al. (DOI 10.1016/j.ccell.2021.01.006), including the published Supplementary Table S2 workbook used for mRNA, protein, phosphosite and metabolite extraction. GLASS longitudinal analyses used the current-release GLASS Synapse resources, including project root syn17038081, current release syn26465623, gene TPM matrix syn69961520, clinical cases syn69931132, clinical surgeries syn69931127, tumor pairs syn69931139 and RNA silver set syn69930927, subject to Synapse access terms.

Processed source data generated by this study for Figures 1-5 and Supplementary Fig. S4A are provided in the Source Data files and mapped in Supplementary Table 5. The frozen LAP3-state gene set, LAP3-state submodule assignments, malignant-state M1/M2/M3 cluster assignments, Core GBmap donor-state projection summaries, GBM-Space topology summaries, CPTAC boundary analyses, GLASS paired analyses, software/package inventory and pre-specified mechanism landscape evidence matrix are provided in Supplementary Tables 1-11. The source data and supplementary tables are derived summaries and do not redistribute controlled-access raw single-cell, spatial, clinical or omics files from third-party repositories.

The final manuscript submission should include the public GitHub repository at https://github.com/lzbkk/glioma-public-data-archive for analysis scripts and a versioned Zenodo record for the generated source-data package and journal-facing supplementary tables. A Zenodo draft has been created at https://zenodo.org/deposit/21161701 with reserved version DOI 10.5281/zenodo.21161701 and concept DOI 10.5281/zenodo.21161700. The draft is currently unsubmitted; after publication, replace the draft URL with the public Zenodo record URL. The public Zenodo package does not redistribute controlled-access raw data, third-party raw omics files or controlled-access spot-level spatial coordinate/score tables. The generated derived source-data/supplementary-table artifacts and public archive documentation are released under CC BY 4.0; software/code licensing is handled through the GitHub repository.
```

## Code Availability

Ready-to-paste working draft:

```text
Code Availability

Analysis scripts used to generate the reported tables, source data and figures are archived in the public GitHub repository at https://github.com/lzbkk/glioma-public-data-archive. The repository includes R scripts for TCGA/CGGA bulk analyses, LAP3-state construction, alternative-anchor benchmarking, submodule decomposition, Core GBmap projection, GBM-Space topology analysis, CPTAC/GLASS projection, figure generation and source-data indexing, together with a README describing expected inputs, access-controlled inputs that must be obtained from third-party repositories, software versions and the order of execution. The code repository does not include controlled-access raw data or third-party datasets that the authors are not permitted to redistribute. Generated source-data and supplementary-table artifacts are separately staged in the Zenodo draft at https://zenodo.org/deposit/21161701 with reserved DOI 10.5281/zenodo.21161701; the Zenodo record should be cited as a final public record only after publication.
```

## Dataset Access Route Table

| Dataset/resource | Role | Access route | Identifier/status | Redistribution decision |
|---|---|---|---|---|
| TCGA GBM/LGG | Bulk discovery and clinical-molecular context | Reused public source | TCGA GBM DOI `10.1038/nature07385`; TCGA LGG DOI `10.1056/NEJMoa1402121` | Do not redistribute full third-party matrices unless licence allows; share derived source data and scripts |
| CGGA mRNAseq_693/325 | External bulk validation | Reused public source | CGGA portal-preferred journal article DOI `10.1016/j.gpb.2020.10.005`; PMID `33662628`; checked 2026-07-02 | Do not redistribute full third-party matrices unless portal terms allow; share derived summaries |
| Core GBmap | Single-cell malignant-state projection | Reused public source | DOI `10.1093/neuonc/noaf113`; local SHA256 recorded | Do not redistribute original H5AD by default; share donor-state summaries and scripts |
| GBM-Space Visium | Spatial topology | Controlled-access repository / third-party restricted | EGA dataset `EGAD00001015527`; study `EGAS00001005801`; DAC `EGAC00001000205`; stable preprint DOI `10.1101/2025.05.13.653495`; no peer-reviewed article found in 2026-07-02 live recheck; EGA page notes two new samples added on 2026-05-19 | Do not redistribute raw H5AD/tarball; readers request access through EGA/DAC; analyses refer to local downloaded/cache version |
| CPTAC GBM | Protein, metabolite and phospho boundary | Reused public source | DOI `10.1016/j.ccell.2021.01.006`; Supplementary Table S2 workbook | Share derived source data; cite original resource |
| GLASS | Longitudinal paired tumor analysis | Reused source with Synapse access terms | `syn17038081`, `syn26465623`, `syn69961520`, `syn69931132`, `syn69931127`, `syn69931139`, `syn69930927` | Share derived paired summaries; readers obtain original data through Synapse |
| Generated Figure 1-5 source data | Direct support for manuscript figures | Within paper/supplement plus public archive recommended | 30 source-data rows in `supplementary_source_data_index_20260701.csv`; Zenodo draft DOI reserved `10.5281/zenodo.21161701`, draft/unsubmitted; S4A spot-level controlled spatial table replaced by public provenance/summary placeholder | Deposit with manuscript Source Data and final repository |
| Supplementary Tables 1-11 | Gene sets, projection tables, software, provenance and mechanism landscape audit | Within paper/supplement plus public archive recommended | `supplementary_table_manifest_20260701.csv`; Zenodo draft DOI reserved `10.5281/zenodo.21161701`, draft/unsubmitted | Deposit as Supplementary Tables and final repository |
| Analysis scripts | Reproducibility | Public code repository plus Zenodo release DOI recommended | GitHub repository `https://github.com/lzbkk/glioma-public-data-archive`; Zenodo draft DOI reserved `10.5281/zenodo.21161701`, draft/unsubmitted | Archive release before submission; exclude credentials and controlled raw data |

## Repository and Citation Actions

1. Create a public release for generated manuscript artifacts before submission.
   - Fixed route: GitHub for code plus Zenodo DOI for a versioned code/source-data release.
   - Minimum contents: `source_data/`, final Supplementary Tables 1-11, figure panel maps, key-results tables, export QC tables, README/data dictionary, script manifest, mechanism landscape audit tables and software/package inventory.

2. Do not upload third-party raw data unless redistribution is explicitly permitted.
   - GBM-Space should stay EGA-routed.
   - GLASS should stay Synapse-routed.
   - TCGA/CGGA/GBmap/CPTAC should be cited and accessed through original resources unless their terms clearly permit redistribution.

3. Add a README to the generated data/code archive.
   - Explain which files reproduce each Figure 1-5 panel and Supplementary Fig. S4A.
   - Explain which third-party inputs must be downloaded separately.
   - State that `LAP3_STATE_UNION` excludes LAP3.
   - Preserve inference-unit definitions: patient/tumor, donor-state, tumor-level spatial summary, CPTAC tumor/patient and GLASS patient pair.

4. Submission-time citation checks still required.
   - GBM-Space was rechecked on 2026-07-02 and still had no peer-reviewed journal article; recheck again immediately before final upload.
   - CGGA portal-preferred citation was rechecked on 2026-07-02 and resolved to the GPB 2021 journal article, DOI `10.1016/j.gpb.2020.10.005`.
   - Confirm whether the final target journal will allow a bioRxiv item parenthetically if no journal article exists; otherwise use EGA accessions and avoid a formal preprint reference.

## FAIR / Metadata Audit

| Check | Current status | Action |
|---|---|---|
| Persistent identifier for generated artifacts | GitHub repository pushed; Zenodo DOI reserved; Zenodo draft validated but unsubmitted | Preview and publish Zenodo draft before final submission |
| Source data mapped to figures | Present | Use `supplementary_source_data_index_20260701.csv` as Source Data map |
| Supplementary table manifest | Present | Convert 10-table manifest into final journal files |
| Third-party data provenance | Mostly present | Keep DOI/accession/Synapse/EGA identifiers in statement |
| Access restrictions explicit | Partial | Make GBM-Space EGA/DAC route explicit; confirm GLASS Synapse terms |
| README/data dictionary | Needed for final archive | Add variable definitions, file descriptions and script-to-output map |
| Licence | Public Zenodo artifacts use CC BY 4.0; software/code licensing handled separately in GitHub | Do not apply open licence to third-party raw data |
| Checksums | Partial | Core GBmap SHA256 recorded; consider checksums for final archive zip |
| Reviewer access | Not set | Use repository private-review link if archive supports it |

## Missing Information / Risk Flags

1. **Blocking before submission:** generated source-data/code archive has a Zenodo draft and reserved DOI, but the public Zenodo record URL and release version are not yet submission-frozen.
2. **Resolved as of 2026-07-03:** GitHub repository URL is fixed as `https://github.com/lzbkk/glioma-public-data-archive`.
3. **Important:** Zenodo generated derived artifacts use CC BY 4.0; confirm whether the GitHub code repository should carry a separate MIT or institutional software licence.
4. **Important:** GBM-Space journal-version status was checked on 2026-07-02 and remains preprint-only; recheck immediately before submission.
5. **Resolved as of 2026-07-02:** CGGA portal-preferred citation is the GPB 2021 journal article DOI 10.1016/j.gpb.2020.10.005, PMID 33662628.
6. **Important:** GLASS Synapse data-use terms should be checked before deciding how much derived patient-pair summary data can be redistributed.
7. **Moderate:** exact data retrieval dates are not fully recoverable from this draft; include them only if target journal requires them and logs/portal records support them.
8. **Moderate:** Neuro-Oncology has a specific AI-use policy. A separate AI-use disclosure may be needed in Methods/Acknowledgements; do not bury that disclosure inside Data Availability.

## 中文核对

1. 这篇纯干实验稿没有产生新的原始人类组学数据，所以 Data Availability 的核心不是“上传原始数据”，而是说明每个公共/受限数据集的来源、登录号和访问路线。
2. 本研究自己产生的内容主要是 derived source data、Supplementary Tables、figure source data、panel maps、key results、QC tables 和分析脚本。这些应该在投稿前打包成一个有 DOI 的公开归档。
3. GBM-Space 原始 Visium 数据不能由我们随意二次分发，应让读者通过 EGA dataset `EGAD00001015527` / study `EGAS00001005801` / DAC `EGAC00001000205` 申请或访问。
4. GLASS 原始数据也不要直接打包进我们的仓库，应保留 Synapse 访问路线。
5. 不能写“数据可向通讯作者合理索取”作为主要方案。这对 Neuro-Oncology 风险较高，除非有明确限制、审核主体和数据使用协议。
6. 归档路线已固定为 `GitHub + Zenodo DOI`。GitHub repository 已创建并推送：`https://github.com/lzbkk/glioma-public-data-archive`。Zenodo draft 已创建：`https://zenodo.org/deposit/21161701`，reserved DOI 为 `10.5281/zenodo.21161701`。下一步是 preview、确认 licence、publish Zenodo record，并把 public record URL/release version 回填到本文档和最终稿。
