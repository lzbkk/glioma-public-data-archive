#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

src <- file.path(
  "Project_Management", "Plans",
  "Manuscript_Continuous_Draft_20260702_polished_v3_NOA_mechanism_landscape.md"
)
out_dir <- file.path(
  "Project_Management", "Submission_Package",
  "NOA_Submission_Working_Package_20260702"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_lines <- function(path) readLines(path, warn = FALSE)

lines <- read_lines(src)

section_bounds <- function(heading) {
  start <- grep(paste0("^## ", heading, "$"), lines)
  if (length(start) != 1) {
    stop("Expected exactly one section: ", heading)
  }
  next_heads <- grep("^## ", lines)
  end_candidates <- next_heads[next_heads > start]
  end <- if (length(end_candidates) == 0) length(lines) + 1 else end_candidates[1]
  c(start + 1, end - 1)
}

get_section <- function(heading) {
  b <- section_bounds(heading)
  if (b[1] > b[2]) character() else lines[b[1]:b[2]]
}

shift_headings <- function(x) {
  x <- sub("^### ", "#### ", x)
  x <- sub("^## ", "### ", x)
  x
}

word_count <- function(x) {
  text <- paste(x, collapse = " ")
  text <- gsub("`[^`]+`", " ", text)
  text <- gsub("[^A-Za-z0-9_/-]+", " ", text)
  words <- unlist(strsplit(trimws(text), "[[:space:]]+"))
  sum(nchar(words) > 0)
}

char_count <- function(x) nchar(x, type = "chars", allowNA = FALSE)

title <- "LAP3 marks a reproducible malignant-microenvironmental state in glioma"
running_title <- "LAP3-state in glioma"
keywords <- c(
  "glioma",
  "LAP3",
  "single-cell RNA-seq",
  "spatial transcriptomics",
  "public multi-omics"
)

structured_abstract <- c(
  "Background: Glioma biomarkers often conflate malignant-cell programs, microenvironmental composition, and spatial organization. We asked whether leucine aminopeptidase 3 (LAP3) marks a reproducible glioma state, rather than a single-gene signal, and how public multi-omic data constrain its mechanistic interpretation.",
  "",
  "Methods: We integrated TCGA and CGGA bulk cohorts, Core GBmap single-cell profiles, GBM-Space Visium spatial transcriptomics, CPTAC proteogenomic and metabolomic data, and GLASS paired longitudinal tumors. We constructed a frozen 207-gene LAP3-state score excluding LAP3, decomposed it into submodules, and projected it across cellular, spatial, protein, metabolite, phospho-readout and longitudinal layers.",
  "",
  "Results: The LAP3-state was strongly coupled to LAP3 across TCGA and CGGA and resolved into malignant-state, myeloid/TAM, anabolic/translation, proteostasis/stress and hypoxia/perinecrotic modules. Core GBmap localized the state across malignant-cell contexts, while GBM-Space linked it most strongly to TAM/myeloid and gliosis-associated topology. CPTAC supported LAP3 mRNA-protein concordance and state-level cross-modal consistency, but direct BCAA metabolites and phospho-mTORC1 readouts did not validate a direct LAP3-BCAA-mTORC1 mechanism. GLASS supported coordinated longitudinal state variation without uniform recurrence activation. A pre-specified mechanism landscape audit prioritized malignant-state, TAM/myeloid, proteostasis/stress and hypoxia-gliosis evidence directions.",
  "",
  "Conclusions: LAP3 is a benchmark-supported anchor of a reproducible malignant-microenvironmental state in glioma. Public data nominate ecosystem and metabolic hypotheses while treating LAP3-BCAA-mTORC1 biology as a future perturbational hypothesis."
)

key_points <- c(
  "LAP3 anchors a reproducible malignant-microenvironmental glioma state.",
  "The state maps to malignant contexts and TAM/myeloid spatial topology.",
  "CPTAC and GLASS constrain direct BCAA/mTORC1 and recurrence claims."
)

importance <- paste(
  "Single-gene glioma biomarkers are difficult to interpret because they may reflect tumor subtype, malignant-cell state, immune composition, spatial niche or pathway activity.",
  "This public-data-only study uses LAP3 as a biologically motivated and experimentally tractable anchor, but does not treat it as a proven causal driver.",
  "By integrating bulk, single-cell, spatial, proteogenomic and longitudinal datasets, we define a frozen LAP3-state excluding LAP3 itself, decompose it into interpretable malignant and microenvironmental modules, and show that its strongest spatial topology is TAM/myeloid and gliosis-linked.",
  "Cross-modal CPTAC and GLASS analyses support protein and longitudinal state coordination while limiting direct BCAA/mTORC1 claims.",
  "The work provides a reusable evidence framework and sharper hypotheses for future perturbational experiments."
)

intro <- get_section("Introduction")
methods <- get_section("Methods")
results <- get_section("Results")
discussion <- get_section("Discussion")
data_availability <- get_section("Data Availability")
code_availability <- get_section("Code Availability")
figure_legends <- get_section("Figure Legends")

body_lines <- c(
  "### Introduction", "",
  intro, "",
  "### Materials and Methods", "",
  shift_headings(methods), "",
  "### Results", "",
  shift_headings(results), "",
  "### Discussion", "",
  discussion
)

required_statements <- c(
  "## Required Statements",
  "",
  "### Funding",
  "",
  "[To be completed by authors.]",
  "",
  "### Conflict of Interest",
  "",
  "[To be completed by authors. All authors should disclose conflicts according to NOA/OUP and ICMJE requirements.]",
  "",
  "### Authorship Contributions",
  "",
  "[To be completed by authors in the format requested by NOA: Contribution or task: author names or initials.]",
  "",
  "### Data Availability",
  "",
  data_availability,
  "",
  "### Code Availability",
  "",
  code_availability,
  "",
  "### AI-Assisted Work Disclosure",
  "",
  "[Author verification required.] AI-assisted tools were used to support project organization, code drafting/debugging, analysis documentation and manuscript language/structure editing. The authors verified all analyses, outputs, claims and wording. No AI tool is listed as an author, and no AI-generated data, images or citations were used as evidence.",
  "",
  "### Acknowledgments",
  "",
  "[To be completed by authors.]"
)

manuscript_lines <- c(
  "# NOA Submission Working Manuscript 20260702",
  "",
  "Status: NOA-formatted working manuscript generated from the v3 NOA mechanism-landscape draft; not submission-frozen.",
  "",
  "## Title Page",
  "",
  paste0("Title: ", title),
  "",
  paste0("Running title: ", running_title),
  "",
  "Authors: [To be completed by authors]",
  "",
  "Affiliations: [To be completed by authors]",
  "",
  "Corresponding author: [Name, full address, telephone and email to be completed by authors]",
  "",
  "Prior publication / concurrent submission: [To be completed in cover letter and title page if applicable.]",
  "",
  paste0("Manuscript text body word count, Introduction through Discussion: ", word_count(body_lines)),
  "",
  "## Abstract and Keywords",
  "",
  structured_abstract,
  "",
  paste0("Keywords: ", paste(keywords, collapse = "; ")),
  "",
  "## Key Points",
  "",
  paste0("- ", key_points),
  "",
  "## Importance of the Study",
  "",
  importance,
  "",
  "## Text",
  "",
  body_lines,
  "",
  required_statements,
  "",
  "## References",
  "",
  "[To be inserted from the citation pass / reference manager export after final journal-style formatting.]",
  "",
  "## Captions for All Illustrations",
  "",
  shift_headings(figure_legends),
  ""
)

manuscript_path <- file.path(out_dir, "NOA_Submission_Working_Manuscript_20260702.md")
writeLines(manuscript_lines, manuscript_path)

cover_letter <- c(
  "# NOA Cover Letter Draft 20260702",
  "",
  "Status: working draft; author names, affiliations, preprint/concurrent-submission status, conflicts and repository DOI must be verified before use.",
  "",
  "Dear Editor,",
  "",
  "We are pleased to submit our manuscript entitled \"LAP3 marks a reproducible malignant-microenvironmental state in glioma\" for consideration as an original article in Neuro-Oncology Advances.",
  "",
  "This public-data-only study addresses a common challenge in glioma biomarker interpretation: single-gene signals may reflect malignant-cell state, immune composition, spatial niche and pathway context rather than a direct causal mechanism. We use LAP3 as a biologically motivated and experimentally tractable anchor, but deliberately avoid presenting it as a proven LAP3-BCAA-mTORC1 causal driver.",
  "",
  "The manuscript integrates TCGA and CGGA bulk cohorts, Core GBmap single-cell profiles, GBM-Space Visium spatial transcriptomics, CPTAC proteogenomic and metabolomic data, and GLASS paired longitudinal tumors. We define a frozen 207-gene LAP3-state score excluding LAP3 itself, decompose it into malignant-state and microenvironmental modules, project it across cellular and spatial contexts, and use CPTAC/GLASS to define cross-modal support and public-data boundaries. A pre-specified mechanism landscape audit prioritizes malignant-state, TAM/myeloid, proteostasis/stress and hypoxia-gliosis evidence directions, while retaining amino-acid/mTORC1 biology as an experimentally testable hypothesis.",
  "",
  "We believe the manuscript fits Neuro-Oncology Advances because it is a CNS tumor-focused, integrative computational study using large public data resources, and because it explicitly distinguishes robust state/ecosystem organization from unsupported causal mechanism claims.",
  "",
  "The manuscript, or any part of it, [has not been previously published and is not under consideration elsewhere / AUTHOR TO VERIFY]. All authors have seen and approved the submitted version of the manuscript [AUTHOR TO VERIFY]. Conflicts of interest are disclosed in the manuscript [AUTHOR TO VERIFY].",
  "",
  "AI-assisted tools were used to support project organization, code drafting/debugging, analysis documentation and manuscript language/structure editing. The authors verified all analyses, outputs, claims and wording. No AI tool is listed as an author, and no AI-generated data, images or citations were used as evidence.",
  "",
  "Generated source data, supplementary tables and analysis code will be archived through a public GitHub repository linked to a versioned Zenodo DOI before final submission [GitHub URL and Zenodo DOI pending]. Controlled-access and third-party raw data are not redistributed and remain available through their original repositories and access conditions.",
  "",
  "Sincerely,",
  "",
  "[Corresponding author name]",
  "[Institution]",
  "[Email]"
)
writeLines(cover_letter, file.path(out_dir, "NOA_Cover_Letter_Draft_20260702.md"))

checklist <- c(
  "# NOA Submission Checklist 20260702",
  "",
  "Status: working checklist based on NOA/OUP author instructions checked on 2026-07-02; recheck the live submission interface before upload.",
  "",
  "## Format Constraints",
  "",
  paste0("- [x] Title <=160 characters: ", char_count(title), " characters."),
  paste0("- [x] Running title <=50 characters: ", char_count(running_title), " characters."),
  paste0("- [x] Structured abstract <=250 words: ", word_count(structured_abstract), " words."),
  "- [x] Up to five keywords: 5 keywords.",
  paste0("- [x] 2-3 key points: ", length(key_points), " key points."),
  paste0("- [x] Each key point <=85 characters: ", paste(char_count(key_points), collapse = ", "), " characters."),
  paste0("- [x] Key points total <=260 characters: ", sum(char_count(key_points)), " characters."),
  paste0("- [x] Importance of the Study <=150 words: ", word_count(importance), " words."),
  paste0("- [x] Text body Introduction through Discussion <=6000 words: ", word_count(body_lines), " words."),
  "- [x] Display items <=6: five main figures planned.",
  "- [x] Core references <=50: 19-item NOA core RIS generated; final journal-style numbering/formatting still pending.",
  "",
  "Official basis checked on 2026-07-02:",
  "",
  "- NOA manuscript preparation instructions: https://academic.oup.com/noa/pages/General_Instructions",
  "- NOA aims/scope page: https://academic.oup.com/noa/pages/About",
  "- NOA online submission system: https://www.editorialmanager.com/noa/",
  "",
  "## Files",
  "",
  "- [x] NOA-formatted working manuscript generated.",
  "- [x] Cover letter draft generated.",
  "- [x] Author submission information checklist generated as a companion file.",
  "- [x] Source Data Package rebuilt with Supplementary Table 11.",
  "- [x] Journal-facing Supplementary Tables rebuilt as 1-11.",
  "- [ ] Final Word `.docx` manuscript file created by authors or later conversion step.",
  "- [ ] Final figure upload files checked against NOA accepted formats.",
  "- [ ] Supplement uploaded separately; decide combined workbook, individual workbooks, or both.",
  "",
  "## Submission Gates",
  "",
  "- [ ] Author list, affiliations and corresponding author details completed.",
  "- [ ] Funding statement completed.",
  "- [ ] Conflict-of-interest statement completed.",
  "- [ ] Author contributions completed in NOA format.",
  "- [ ] Concurrent submission / prior publication statement confirmed.",
  "- [ ] AI-assisted work disclosure approved by authors and included in cover letter plus manuscript.",
  "- [ ] GitHub URL and Zenodo DOI created and backfilled.",
  "- [x] GBM-Space journal-version status rechecked on 2026-07-02; no peer-reviewed article found, final upload-day recheck still required.",
  "- [x] CGGA portal-preferred citation rechecked on 2026-07-02; use Zhao et al. GPB 2021, DOI 10.1016/j.gpb.2020.10.005, PMID 33662628.",
  "- [ ] Final reference list formatted and checked.",
  "- [ ] Final PDF/Word preview checked for figure/table/caption order.",
  "",
  "## Boundary Locks",
  "",
  "- [x] Do not write LAP3 as a proven causal driver.",
  "- [x] Do not claim direct LAP3-BCAA-mTORC1 validation from public data.",
  "- [x] Keep Figure 4 as spatial topology / neighborhood association.",
  "- [x] Keep Figure 5 as cross-modal support and evidence boundary."
)
writeLines(checklist, file.path(out_dir, "NOA_Submission_Checklist_20260702.md"))

qa <- data.frame(
  item = c(
    "title_characters",
    "running_title_characters",
    "abstract_words",
    "importance_words",
    "key_points_count",
    "key_points_total_characters",
    "body_words_intro_to_discussion",
    "display_items",
    "supplementary_table_items",
    "supplementary_expanded_inputs"
  ),
  value = c(
    char_count(title),
    char_count(running_title),
    word_count(structured_abstract),
    word_count(importance),
    length(key_points),
    sum(char_count(key_points)),
    word_count(body_lines),
    5,
    11,
    61
  ),
  limit_or_note = c(
    "<=160",
    "<=50",
    "<=250",
    "<=150",
    "2-3",
    "<=260",
    "<=6000",
    "<=6",
    "current package",
    "current package"
  ),
  stringsAsFactors = FALSE
)
write.csv(qa, file.path(out_dir, "NOA_Submission_Format_QA_20260702.csv"), row.names = FALSE)

readme <- c(
  "# NOA Submission Working Package 20260702",
  "",
  "This directory contains the working Neuro-Oncology Advances submission package generated from the v3 NOA mechanism-landscape manuscript draft.",
  "",
  "## Files",
  "",
  "- `NOA_Submission_Working_Manuscript_20260702.md`: NOA-ordered manuscript scaffold with structured abstract, keywords, Key Points, Importance of the Study, Materials and Methods before Results, required statements and figure captions.",
  "- `NOA_Cover_Letter_Draft_20260702.md`: cover letter draft with author-verification placeholders.",
  "- `NOA_Submission_Checklist_20260702.md`: format and submission-gate checklist.",
  "- `NOA_Submission_Format_QA_20260702.csv`: machine-readable word/character/display checks.",
  "- `../Author_Submission_Info_Checklist_20260702.md`: companion author-fillable checklist for author roster, affiliations, corresponding author details, funding, conflicts, contributions, ethics, AI disclosure, cover-letter declarations and upload-day sign-off.",
  "",
  "## Not Submission-Frozen",
  "",
  "- Author list, affiliations, corresponding-author details, funding, conflicts, author contributions and acknowledgments require author input.",
  "- GitHub URL and Zenodo DOI are pending.",
  "- GBM-Space journal-version status was rechecked on 2026-07-02 and remained preprint-only; final upload-day recheck is still required.",
  "- CGGA portal-preferred citation was rechecked on 2026-07-02 and resolved to Zhao et al. GPB 2021, DOI 10.1016/j.gpb.2020.10.005, PMID 33662628.",
  "- The final Word `.docx` and journal-upload file set have not yet been created.",
  "",
  "## Rebuild",
  "",
  "Run from the repository root:",
  "",
  "```bash",
  "Rscript --vanilla Project_Management/Operations/build_noa_submission_working_package.R",
  "```"
)
writeLines(readme, file.path(out_dir, "README.md"))

cat("NOA submission working package built at:", out_dir, "\n")
cat("Structured abstract words:", word_count(structured_abstract), "\n")
cat("Importance words:", word_count(importance), "\n")
cat("Body words:", word_count(body_lines), "\n")
cat("Key point characters:", paste(char_count(key_points), collapse = ", "), "\n")
