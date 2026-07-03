# NOA Reference Final Pass 20260702

Date: 2026-07-02

## Routing

```text
Task: NOA-facing reference final pass for active manuscript package
Skill: nature-academic-search
Runtime: foreground citation/provenance task; no tmux required
Target journal: Neuro-Oncology Advances
Inputs:
  Project_Management/Plans/Manuscript_Continuous_Draft_20260702_polished_v3_NOA_mechanism_landscape.md
  Project_Management/Plans/Data_Code_Availability_Draft_20260702.md
  Project_Management/References/Claim_to_Reference_Citation_Pass_20260701.md
  Project_Management/References/GBM_Space_Manuscript_Preprint_Check_20260701.md
Outputs:
  Project_Management/References/NOA_Reference_Final_Pass_20260702.md
  Project_Management/References/references_noa_core_20260702.ris
```

## Executive Verdict

```text
GO: CGGA citation is now corrected from the 2020 bioRxiv preprint to the
     portal-preferred GPB 2021 journal article.

GO WITH CAVEAT: GBM-Space remains preprint + EGA accession anchored.
     No peer-reviewed journal article was found in the 2026-07-02 live recheck.

BOUNDARY: The EGA dataset page reports two new samples added on 2026-05-19.
     The manuscript must state that the analysis uses the local downloaded/cache
     version: 97 sections, 12 tumors, 343,799 in-tissue spots.
```

This pass supersedes the CGGA citation entry in
`Claim_to_Reference_Citation_Pass_20260701.md`. Historical notes may remain in
older files, but the active manuscript and NOA package should use the records
below.

## Live Checks Completed

| Item | Source checked | Finding | Manuscript action |
|---|---|---|---|
| CGGA preferred citation | CGGA portal, CrossRef DOI, PubMed PMID | Portal-preferred record is Zhao et al., *Genomics, Proteomics & Bioinformatics* 2021, DOI `10.1016/j.gpb.2020.10.005`, PMID `33662628` | Replace the CGGA 2020 bioRxiv citation in active manuscript/Data Availability |
| GBM-Space article status | CrossRef DOI/title searches, PubMed/title searches, EGA dataset/study pages | No peer-reviewed journal article found on 2026-07-02; DOI `10.1101/2025.05.13.653495` remains bioRxiv posted content | Cite stable preprint only if journal permits; always retain EGA dataset/study/DAC accessions |
| GBM-Space EGA accessions | EGA dataset `EGAD00001015527`, study `EGAS00001005801`, DAC `EGAC00001000205` | Dataset is official Visium GBM-Space record; EGA page notes two new samples added on 2026-05-19 | Tie reported analyses to local downloaded/cache version rather than implying all current EGA files were analyzed |
| EGA archive citation | DOI `10.1093/nar/gkab1059` | Formal EGA archive paper verified | Use as archive citation when a formal database citation is needed |

## Corrected Core Reference Decisions

| Resource/claim | Active citation decision | Status |
|---|---|---|
| TCGA GBM | TCGA Research Network, *Nature* 2008, DOI `10.1038/nature07385` | keep |
| TCGA LGG | TCGA Research Network, *NEJM* 2015, DOI `10.1056/NEJMoa1402121` | keep |
| CGGA | Zhao et al., *Genomics, Proteomics & Bioinformatics* 2021, DOI `10.1016/j.gpb.2020.10.005`, PMID `33662628` | corrected; replaces `10.1101/2020.01.20.911982` as primary citation |
| Core/Extended GBmap | Ruiz-Moreno et al., *Neuro-Oncology* 2025, DOI `10.1093/neuonc/noaf113` | keep |
| GBM malignant-cell states | Neftel et al., *Cell* 2019, DOI `10.1016/j.cell.2019.06.024` | keep |
| GBM-Space spatial/trajectory manuscript | de Jong et al., bioRxiv 2025, DOI `10.1101/2025.05.13.653495` | preprint-only as of 2026-07-02 |
| GBM-Space data access | EGA dataset `EGAD00001015527`, study `EGAS00001005801`, DAC `EGAC00001000205`; EGA archive DOI `10.1093/nar/gkab1059` | stable data anchor |
| CPTAC GBM | Wang et al., *Cancer Cell* 2021, DOI `10.1016/j.ccell.2021.01.006` | keep |
| GLASS | Barthel et al., *Nature* 2019, DOI `10.1038/s41586-019-1775-1` | keep |
| LAP3 glioma background | He et al., *International Journal of Biological Macromolecules* 2015, DOI `10.1016/j.ijbiomac.2014.10.021` | background only |
| Amino-acid/mTORC1 rationale | Sancak et al. 2008, DOI `10.1126/science.1157535`; Wolfson et al. 2016, DOI `10.1126/science.aab2674` | rationale only; not evidence that this manuscript proves the mechanism |

## Active Manuscript Edits Completed

```text
Updated:
  Project_Management/Plans/Manuscript_Continuous_Draft_20260702_polished_v3_NOA_mechanism_landscape.md
  Project_Management/Plans/Data_Code_Availability_Draft_20260702.md
  Project_Management/Operations/build_noa_submission_working_package.R

Expected rebuilt outputs:
  Project_Management/Submission_Package/NOA_Submission_Working_Package_20260702/NOA_Submission_Working_Manuscript_20260702.md
  Project_Management/Submission_Package/NOA_Submission_Working_Package_20260702/NOA_Submission_Checklist_20260702.md
  Project_Management/Submission_Package/NOA_Submission_Working_Package_20260702/README.md
```

## Submission-Day Checks Still Required

1. Recheck DOI `10.1101/2025.05.13.653495`, PubMed, CrossRef and the GBM-Space/EGA pages for a journal version immediately before final upload.
2. Confirm NOA's current policy on preprints in the reference list. If NOA does not want a preprint reference, keep the GBM-Space citation in text/Data Availability as EGA accession plus parenthetical preprint status, and cite the EGA archive paper formally.
3. Backfill GitHub URL, Zenodo DOI, author statements, funding, conflicts and final reference numbering.
4. Do not cite BioStudies DOI `10.6019/s-biad2898` for the current Visium analysis; it is platform-mismatched to the current GBM-Space Visium workflow.

## Source Links Used

```text
CGGA portal:
https://www.cgga.org.cn/

CGGA formal journal DOI:
https://doi.org/10.1016/j.gpb.2020.10.005

GBM-Space stable preprint DOI:
https://doi.org/10.1101/2025.05.13.653495

GBM-Space EGA dataset:
https://ega-archive.org/datasets/EGAD00001015527

GBM-Space EGA study:
https://ega-archive.org/studies/EGAS00001005801

EGA archive paper DOI:
https://doi.org/10.1093/nar/gkab1059
```

