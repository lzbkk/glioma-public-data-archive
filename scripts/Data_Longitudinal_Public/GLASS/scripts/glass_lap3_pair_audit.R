#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

set.seed(20260630)

root_dir <- "Data_Longitudinal_Public/GLASS"
data_dir <- file.path(root_dir, "data")
result_dir <- file.path(root_dir, "results", "LAP3_Longitudinal_Feasibility")
table_dir <- file.path(result_dir, "tables")
log_dir <- file.path(result_dir, "logs")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

expr_file <- file.path(data_dir, "expression", "gene_tpm_matrix_all_samples.tsv")
case_file <- file.path(data_dir, "tables", "clinical_cases.csv")
surgery_file <- file.path(data_dir, "tables", "clinical_surgeries.csv")
pair_file <- file.path(data_dir, "tables", "analysis_rna_silver_set.csv")
stopifnot(file.exists(expr_file), file.exists(case_file))
stopifnot(file.exists(surgery_file), file.exists(pair_file))

log_file <- file.path(log_dir, "glass_lap3_pair_audit.log")
log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste(...))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

normalize_barcode <- function(x) {
  gsub("\\.", "-", x)
}

to_sample_barcode <- function(x) {
  sub("-[0-9]{2}R-RNA-[A-Z0-9]+$", "", x)
}

sample_type_from_barcode <- function(x) {
  parts <- strsplit(x, "-", fixed = TRUE)
  vapply(parts, function(z) {
    hit <- grep("^(TP|R[0-9]+)$", z, value = TRUE)
    if (length(hit) == 0L) NA_character_ else hit[1]
  }, character(1))
}

clean_text <- function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  x
}

log_msg("Reading GLASS gene TPM matrix")
expr <- fread(expr_file, showProgress = FALSE)
stopifnot(names(expr)[1] == "Gene_symbol")
lap3_rows <- expr[Gene_symbol == "LAP3"]
if (nrow(lap3_rows) != 1L) {
  stop("Expected exactly one LAP3 row; found ", nrow(lap3_rows))
}

expr_samples <- data.table(
  expression_barcode = names(lap3_rows)[-1],
  lap3_tpm = as.numeric(lap3_rows[1, -1])
)
expr_samples[, normalized_barcode := normalize_barcode(expression_barcode)]
expr_samples[, sample_barcode := to_sample_barcode(normalized_barcode)]
expr_samples[, sample_type := sample_type_from_barcode(normalized_barcode)]
stopifnot(!anyDuplicated(expr_samples$normalized_barcode))

cases <- fread(case_file)
surgeries <- fread(surgery_file)
pairs <- fread(pair_file)

for (col in c(
  "case_barcode", "sample_barcode", "histology", "grade", "idh_status",
  "codel_status", "who_classification", "mgmt_methylation",
  "idh_codel_subtype"
)) {
  if (col %in% names(surgeries)) {
    set(surgeries, j = col, value = clean_text(surgeries[[col]]))
  }
}

pairs[, `:=`(
  pair_row = .I,
  sample_barcode_a = to_sample_barcode(tumor_barcode_a),
  sample_barcode_b = to_sample_barcode(tumor_barcode_b),
  sample_type_a = sample_type_from_barcode(tumor_barcode_a),
  sample_type_b = sample_type_from_barcode(tumor_barcode_b)
)]

expr_a <- copy(expr_samples)[, .(
  tumor_barcode_a = normalized_barcode,
  expression_barcode_a = expression_barcode,
  lap3_tpm_a = lap3_tpm
)]
expr_b <- copy(expr_samples)[, .(
  tumor_barcode_b = normalized_barcode,
  expression_barcode_b = expression_barcode,
  lap3_tpm_b = lap3_tpm
)]

pair_x <- merge(pairs, expr_a, by = "tumor_barcode_a", all.x = TRUE)
pair_x <- merge(pair_x, expr_b, by = "tumor_barcode_b", all.x = TRUE)

surgery_keep <- c(
  "case_barcode", "sample_barcode", "surgery_number", "surgical_interval_mo",
  "histology", "grade", "idh_status", "codel_status", "who_classification",
  "mgmt_methylation", "idh_codel_subtype", "treatment_tmz",
  "treatment_radiotherapy", "treatment_alkylating_agent"
)
surgery_meta <- surgeries[, ..surgery_keep]

surgery_a <- copy(surgery_meta)
setnames(
  surgery_a,
  setdiff(names(surgery_a), "case_barcode"),
  paste0(setdiff(names(surgery_a), "case_barcode"), "_a")
)
surgery_b <- copy(surgery_meta)
setnames(
  surgery_b,
  setdiff(names(surgery_b), "case_barcode"),
  paste0(setdiff(names(surgery_b), "case_barcode"), "_b")
)

pair_x <- merge(
  pair_x,
  surgery_a,
  by = c("case_barcode", "sample_barcode_a"),
  all.x = TRUE,
  allow.cartesian = TRUE
)
pair_x <- merge(
  pair_x,
  surgery_b,
  by = c("case_barcode", "sample_barcode_b"),
  all.x = TRUE,
  allow.cartesian = TRUE
)

pair_x[, `:=`(
  expression_match_a = is.finite(lap3_tpm_a),
  expression_match_b = is.finite(lap3_tpm_b),
  surgery_match_a = !is.na(surgery_number_a),
  surgery_match_b = !is.na(surgery_number_b),
  is_tp_recurrence = sample_type_a == "TP" & grepl("^R[0-9]+$", sample_type_b),
  log2_lap3_a = log2(lap3_tpm_a + 1),
  log2_lap3_b = log2(lap3_tpm_b + 1)
)]
pair_x[, delta_log2_lap3 := log2_lap3_b - log2_lap3_a]

pair_x[, molecular_group := fcase(
  idh_status_a == "IDHwt", "IDH-wildtype",
  idh_status_a == "IDHmut" & codel_status_a == "codel",
  "IDH-mutant, 1p/19q-codeleted",
  idh_status_a == "IDHmut" & codel_status_a == "noncodel",
  "IDH-mutant astrocytoma",
  idh_codel_subtype_a == "IDHwt", "IDH-wildtype",
  idh_codel_subtype_a == "IDHmut-codel",
  "IDH-mutant, 1p/19q-codeleted",
  idh_codel_subtype_a == "IDHmut-noncodel",
  "IDH-mutant astrocytoma",
  default = "Unknown"
)]

pair_x[, strict_stratum := fcase(
  molecular_group == "IDH-wildtype" &
    (
      histology_a %chin% c("Glioblastoma", "Gliosarcoma") |
        grepl("Glioblastoma", who_classification_a)
    ),
  "IDH-wildtype GBM",
  molecular_group == "IDH-mutant, 1p/19q-codeleted",
  "IDH-mutant, 1p/19q-codeleted oligodendroglioma",
  molecular_group == "IDH-mutant astrocytoma" &
    (
      grade_a == "IV" |
        grepl("Glioblastoma", who_classification_a)
    ),
  "IDH-mutant astrocytoma, grade 4",
  molecular_group == "IDH-mutant astrocytoma",
  "IDH-mutant astrocytoma, grade 2/3",
  molecular_group == "IDH-wildtype",
  "IDH-wildtype non-GBM/uncertain",
  default = "Unknown"
)]

pair_x[, `:=`(
  idh_status_consistent = (
    !is.na(idh_status_a) & !is.na(idh_status_b) &
      idh_status_a == idh_status_b
  ),
  codel_status_consistent = (
    !is.na(codel_status_a) & !is.na(codel_status_b) &
      codel_status_a == codel_status_b
  )
)]

pair_x[, strict_evaluable := (
  is_tp_recurrence &
    expression_match_a & expression_match_b &
    surgery_match_a & surgery_match_b
)]

strict_candidates <- pair_x[strict_evaluable == TRUE]
setorder(
  strict_candidates,
  case_barcode,
  surgery_number_b,
  surgical_interval_mo_b,
  pair_row
)
strict_primary <- strict_candidates[, .SD[1], by = case_barcode]
strict_primary[, primary_set := TRUE]

pair_x[, primary_set := FALSE]
if (nrow(strict_primary) > 0L) {
  pair_x[pair_row %in% strict_primary$pair_row, primary_set := TRUE]
}

audit_summary <- data.table(
  metric = c(
    "expression_samples",
    "lap3_rows",
    "rna_silver_pair_rows",
    "rna_silver_unique_cases",
    "pair_rows_both_expression_matched",
    "pair_rows_both_surgeries_matched",
    "tp_recurrence_pair_rows",
    "strict_evaluable_pair_rows",
    "strict_primary_unique_cases"
  ),
  value = c(
    nrow(expr_samples),
    nrow(lap3_rows),
    nrow(pairs),
    uniqueN(pairs$case_barcode),
    pair_x[expression_match_a & expression_match_b, .N],
    pair_x[surgery_match_a & surgery_match_b, .N],
    pair_x[is_tp_recurrence == TRUE, .N],
    pair_x[strict_evaluable == TRUE, .N],
    nrow(strict_primary)
  )
)

group_counts <- strict_primary[, .(
  n_pairs = .N,
  median_baseline_tpm = median(lap3_tpm_a, na.rm = TRUE),
  median_recurrence_tpm = median(lap3_tpm_b, na.rm = TRUE),
  median_delta_log2 = median(delta_log2_lap3, na.rm = TRUE),
  n_increased = sum(delta_log2_lap3 > 0, na.rm = TRUE),
  n_decreased = sum(delta_log2_lap3 < 0, na.rm = TRUE),
  n_unchanged = sum(delta_log2_lap3 == 0, na.rm = TRUE),
  n_idh_consistent = sum(idh_status_consistent, na.rm = TRUE),
  n_codel_consistent = sum(codel_status_consistent, na.rm = TRUE)
), by = strict_stratum][order(strict_stratum)]

paired_tests <- rbindlist(lapply(
  c("All", sort(unique(strict_primary$strict_stratum))),
  function(group) {
    dat <- if (group == "All") {
      strict_primary
    } else {
      strict_primary[strict_stratum == group]
    }
    if (nrow(dat) < 3L) {
      return(NULL)
    }
    test <- suppressWarnings(wilcox.test(
      dat$log2_lap3_b,
      dat$log2_lap3_a,
      paired = TRUE,
      exact = FALSE
    ))
    data.table(
      strict_stratum = group,
      n_pairs = nrow(dat),
      median_baseline_tpm = median(dat$lap3_tpm_a),
      median_recurrence_tpm = median(dat$lap3_tpm_b),
      median_delta_log2 = median(dat$delta_log2_lap3),
      p_value = test$p.value
    )
  }
), fill = TRUE)
if (nrow(paired_tests) > 0L) {
  paired_tests[, fdr := p.adjust(p_value, method = "BH")]
}

fwrite(expr_samples, file.path(table_dir, "glass_lap3_expression_samples.csv"))
fwrite(pair_x, file.path(table_dir, "glass_rna_pair_crosswalk.csv"))
fwrite(strict_primary, file.path(table_dir, "glass_strict_primary_recurrence_pairs.csv"))
fwrite(audit_summary, file.path(table_dir, "glass_pair_audit_summary.csv"))
fwrite(group_counts, file.path(table_dir, "glass_pair_counts_by_molecular_group.csv"))
fwrite(paired_tests, file.path(table_dir, "glass_lap3_paired_tests.csv"))

log_msg(
  "Completed:",
  " expression_samples=", nrow(expr_samples),
  " pair_rows=", nrow(pairs),
  " strict_primary_cases=", nrow(strict_primary)
)
