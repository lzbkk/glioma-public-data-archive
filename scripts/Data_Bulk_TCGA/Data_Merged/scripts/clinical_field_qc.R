#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

setwd("/home/lzb/glioma/Data_Bulk_TCGA/Data_Merged")

out_dir <- file.path("results", "Clinical_Field_QC")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

clean_chr <- function(x) {
  x <- as.character(x)
  x[x == "" | x == "NA" | is.na(x)] <- NA_character_
  x
}

mode_value <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[1]
}

clinical <- readRDS(file.path("data_analysis", "clinical_glioma_uni.rds"))

column_qc <- data.frame(
  column = names(clinical),
  class = sapply(clinical, function(x) paste(class(x), collapse = "|")),
  is_paper = grepl("^paper_", names(clinical)),
  missing = sapply(clinical, function(x) sum(is.na(x))),
  missing_pct = round(sapply(clinical, function(x) mean(is.na(x)) * 100), 2),
  n_unique = sapply(clinical, function(x) {
    if (is.list(x)) NA_integer_ else length(unique(x[!is.na(x)]))
  }),
  example_values = sapply(clinical, function(x) {
    if (is.list(x)) return("<list>")
    paste(head(unique(as.character(x[!is.na(x)])), 8), collapse = " | ")
  }),
  stringsAsFactors = FALSE
) %>%
  arrange(is_paper, missing_pct, column)
write.csv(column_qc, file.path(out_dir, "clinical_column_qc_summary.csv"), row.names = FALSE)

key_fields <- c(
  "cohort", "barcode", "patient", "shortLetterCode",
  "age_at_diagnosis", "paper_Age..years.at.diagnosis.",
  "vital_status", "days_to_death", "days_to_last_follow_up",
  "paper_Survival..months.", "paper_Vital.status..1.dead.",
  "tumor_grade", "paper_Grade", "primary_diagnosis", "morphology",
  "paper_IDH.status", "paper_X1p.19q.codeletion", "paper_IDH.codel.subtype",
  "paper_MGMT.promoter.status", "paper_Transcriptome.Subtype",
  "paper_Pan.Glioma.RNA.Expression.Cluster",
  "paper_ESTIMATE.stromal.score", "paper_ESTIMATE.immune.score",
  "paper_ESTIMATE.combined.score"
)

missing_by_cohort <- lapply(intersect(key_fields, names(clinical)), function(f) {
  clinical %>%
    group_by(cohort) %>%
    summarise(
      variable = f,
      missing = sum(is.na(.data[[f]])),
      total = n(),
      missing_pct = round(mean(is.na(.data[[f]])) * 100, 2),
      .groups = "drop"
    )
}) %>%
  bind_rows() %>%
  select(variable, cohort, missing, total, missing_pct)
write.csv(missing_by_cohort, file.path(out_dir, "clinical_key_fields_missing_by_cohort.csv"), row.names = FALSE)

raw_age_years <- clinical$age_at_diagnosis / 365.25
paper_age_years <- suppressWarnings(as.numeric(clinical$paper_Age..years.at.diagnosis.))

raw_event <- case_when(
  clinical$vital_status == "Dead" ~ 1L,
  clinical$vital_status == "Alive" ~ 0L,
  TRUE ~ NA_integer_
)
paper_event <- suppressWarnings(as.integer(clinical$paper_Vital.status..1.dead.))

raw_os_days <- ifelse(raw_event == 1L, clinical$days_to_death, clinical$days_to_last_follow_up)
raw_os_months <- raw_os_days / 30.4375
raw_os_months[!is.na(raw_os_months) & raw_os_months <= 0] <- NA_real_
paper_os_months <- suppressWarnings(as.numeric(clinical$paper_Survival..months.))
paper_os_months[!is.na(paper_os_months) & paper_os_months <= 0] <- NA_real_

grade_raw <- clean_chr(clinical$tumor_grade)
grade_paper <- clean_chr(clinical$paper_Grade)
grade_from_diagnosis <- case_when(
  clinical$cohort == "GBM" ~ "G4",
  grepl("glioblastoma", clinical$primary_diagnosis, ignore.case = TRUE) ~ "G4",
  grepl("anaplastic", clinical$primary_diagnosis, ignore.case = TRUE) ~ "G3",
  TRUE ~ NA_character_
)
grade_unified <- coalesce(grade_paper, grade_raw, grade_from_diagnosis)
grade_source <- case_when(
  !is.na(grade_paper) ~ "paper_Grade",
  is.na(grade_paper) & !is.na(grade_raw) ~ "tumor_grade",
  is.na(grade_paper) & is.na(grade_raw) & !is.na(grade_from_diagnosis) ~ "primary_diagnosis/cohort",
  TRUE ~ NA_character_
)

age_unified <- coalesce(raw_age_years, paper_age_years)
age_source <- case_when(
  !is.na(raw_age_years) ~ "age_at_diagnosis",
  is.na(raw_age_years) & !is.na(paper_age_years) ~ "paper_Age",
  TRUE ~ NA_character_
)

os_months_unified <- coalesce(raw_os_months, paper_os_months)
os_time_source <- case_when(
  !is.na(raw_os_months) ~ "GDC days_to_death/last_follow_up",
  is.na(raw_os_months) & !is.na(paper_os_months) ~ "paper_Survival_months",
  TRUE ~ NA_character_
)

os_event_unified <- coalesce(raw_event, paper_event)
os_event_source <- case_when(
  !is.na(raw_event) ~ "GDC vital_status",
  is.na(raw_event) & !is.na(paper_event) ~ "paper_Vital_status",
  TRUE ~ NA_character_
)

analysis_clinical <- clinical %>%
  transmute(
    barcode,
    patient,
    sample,
    shortLetterCode,
    cohort,
    primary_diagnosis,
    morphology,
    tissue_or_organ_of_origin,
    age_years = age_unified,
    age_source,
    os_months = os_months_unified,
    os_event = os_event_unified,
    os_time_source,
    os_event_source,
    vital_status_raw = vital_status,
    grade = factor(grade_unified, levels = c("G2", "G3", "G4"), ordered = TRUE),
    grade_source,
    tumor_grade_raw = tumor_grade,
    paper_Grade,
    idh_status = factor(clean_chr(paper_IDH.status), levels = c("Mutant", "WT")),
    codel_1p19q = factor(clean_chr(paper_X1p.19q.codeletion), levels = c("codel", "non-codel")),
    idh_codel_subtype = factor(
      clean_chr(paper_IDH.codel.subtype),
      levels = c("IDHmut-codel", "IDHmut-non-codel", "IDHwt")
    ),
    mgmt_status = factor(clean_chr(paper_MGMT.promoter.status), levels = c("Methylated", "Unmethylated")),
    transcriptome_subtype = clean_chr(paper_Transcriptome.Subtype),
    pan_glioma_rna_cluster = clean_chr(paper_Pan.Glioma.RNA.Expression.Cluster),
    estimate_stromal_score = suppressWarnings(as.numeric(paper_ESTIMATE.stromal.score)),
    estimate_immune_score = suppressWarnings(as.numeric(paper_ESTIMATE.immune.score)),
    estimate_combined_score = suppressWarnings(as.numeric(paper_ESTIMATE.combined.score))
  )

saveRDS(analysis_clinical, file.path(out_dir, "clinical_glioma_analysis_fields.rds"))
write.csv(analysis_clinical, file.path(out_dir, "clinical_glioma_analysis_fields.csv"), row.names = FALSE)

comparison <- list(
  grade_cross_tab = table(raw = clinical$tumor_grade, paper = clinical$paper_Grade, useNA = "ifany"),
  primary_diagnosis_grade = table(primary = clinical$primary_diagnosis, paper = clinical$paper_Grade, useNA = "ifany"),
  vital_cross_tab = table(raw = clinical$vital_status, paper = clinical$paper_Vital.status..1.dead., useNA = "ifany")
)

age_idx <- !is.na(raw_age_years) & !is.na(paper_age_years)
surv_idx <- !is.na(raw_os_months) & !is.na(paper_os_months)

qc_metrics <- data.frame(
  metric = c(
    "raw_age_available",
    "paper_age_available",
    "age_raw_paper_cor",
    "raw_survival_available",
    "paper_survival_available",
    "survival_raw_paper_cor",
    "grade_unified_available",
    "idh_available",
    "1p19q_available",
    "mgmt_available"
  ),
  value = c(
    sum(!is.na(raw_age_years)),
    sum(!is.na(paper_age_years)),
    ifelse(sum(age_idx) > 2, cor(raw_age_years[age_idx], paper_age_years[age_idx]), NA),
    sum(!is.na(raw_os_months) & !is.na(raw_event)),
    sum(!is.na(paper_os_months) & !is.na(paper_event)),
    ifelse(sum(surv_idx) > 2, cor(raw_os_months[surv_idx], paper_os_months[surv_idx]), NA),
    sum(!is.na(analysis_clinical$grade)),
    sum(!is.na(analysis_clinical$idh_status)),
    sum(!is.na(analysis_clinical$codel_1p19q)),
    sum(!is.na(analysis_clinical$mgmt_status))
  )
)
write.csv(qc_metrics, file.path(out_dir, "clinical_qc_key_metrics.csv"), row.names = FALSE)

field_recommendations <- data.frame(
  concept = c(
    "sample identity",
    "cohort",
    "age",
    "overall survival time",
    "overall survival event",
    "grade",
    "primary diagnosis / histology",
    "IDH status",
    "1p/19q codeletion",
    "IDH/1p19q molecular subtype",
    "MGMT promoter status",
    "transcriptome subtype",
    "immune/stromal score",
    "treatment",
    "tumor size / sample dimensions"
  ),
  recommended_field = c(
    "barcode, patient, sample, shortLetterCode",
    "cohort",
    "age_years",
    "os_months",
    "os_event",
    "grade",
    "primary_diagnosis, morphology",
    "idh_status",
    "codel_1p19q",
    "idh_codel_subtype",
    "mgmt_status",
    "transcriptome_subtype / pan_glioma_rna_cluster",
    "estimate_stromal_score / estimate_immune_score",
    "not recommended for first-pass analysis",
    "not recommended for first-pass analysis"
  ),
  source_priority = c(
    "GDC/raw identifiers",
    "project-derived TCGA-GBM/TCGA-LGG",
    "GDC age_at_diagnosis first; paper age fallback",
    "GDC days_to_death/days_to_last_follow_up first; paper survival fallback",
    "GDC vital_status first; paper vital fallback",
    "paper_Grade first; tumor_grade and primary_diagnosis/cohort fallback",
    "GDC/raw diagnosis fields",
    "paper molecular annotation; no complete raw alternative in current table",
    "paper molecular annotation; no complete raw alternative in current table",
    "paper molecular annotation",
    "paper molecular annotation",
    "paper annotation; use as exploratory/sensitivity",
    "paper ESTIMATE; use as exploratory/sensitivity",
    "list-column treatments is heterogeneous and requires separate parsing",
    "many missing values and unclear direct relevance"
  ),
  role = c(
    "mandatory",
    "mandatory",
    "main covariate",
    "main survival endpoint",
    "main survival endpoint",
    "main covariate / stratification",
    "descriptive / sensitivity",
    "main molecular covariate",
    "main molecular covariate",
    "stratification",
    "main molecular covariate",
    "exploratory",
    "exploratory",
    "defer",
    "defer"
  ),
  stringsAsFactors = FALSE
)
write.csv(field_recommendations, file.path(out_dir, "clinical_field_recommendations.csv"), row.names = FALSE)

sink(file.path(out_dir, "clinical_field_qc.log"))
cat("Clinical field QC for glioma analysis\n")
cat("Run time:", format(Sys.time()), "\n\n")
cat("Input dimensions:", dim(clinical), "\n\n")
cat("Key QC metrics:\n")
print(qc_metrics)
cat("\nGrade raw vs paper:\n")
print(comparison$grade_cross_tab)
cat("\nPrimary diagnosis vs paper grade:\n")
print(comparison$primary_diagnosis_grade)
cat("\nVital status raw vs paper event:\n")
print(comparison$vital_cross_tab)
cat("\nUnified analysis field missingness:\n")
print(data.frame(
  variable = names(analysis_clinical),
  missing = sapply(analysis_clinical, function(x) sum(is.na(x))),
  total = nrow(analysis_clinical)
))
cat("\nField recommendations:\n")
print(field_recommendations)
sink()

message("Done. Outputs written to: ", normalizePath(out_dir))
