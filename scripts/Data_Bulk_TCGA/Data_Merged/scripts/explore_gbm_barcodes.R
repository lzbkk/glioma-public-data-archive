## Explore TCGA GBM barcode structure and duplicate-sample retention evidence.
## This script is intentionally read-only for source RDS files.

options(stringsAsFactors = FALSE)

work_dir <- "~/glioma/Data_Bulk_TCGA/Data_Merged"
setwd(path.expand(work_dir))

log_dir <- "logs"
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)
run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")

summary_log <- file.path(log_dir, paste0("gbm_barcode_summary_", run_id, ".log"))
sink(summary_log, split = TRUE)
message_log_con <- file(summary_log, open = "at")
sink(message_log_con, type = "message")
on.exit({
  while (sink.number(type = "message") > 0) sink(type = "message")
  while (sink.number() > 0) sink()
  close(message_log_con)
}, add = TRUE)

cat("GBM barcode exploration\n")
cat("Run time:", as.character(Sys.time()), "\n")
cat("Working directory:", getwd(), "\n\n")

parse_tcga_barcode <- function(barcode) {
  data.frame(
    barcode = barcode,
    participant = substr(barcode, 1, 12),
    sample_type_code = substr(barcode, 14, 15),
    vial = substr(barcode, 16, 16),
    portion = substr(barcode, 18, 19),
    analyte = substr(barcode, 20, 20),
    plate = substr(barcode, 22, 25),
    center = substr(barcode, 27, 28),
    stringsAsFactors = FALSE
  )
}

sample_type_priority <- function(short_code) {
  ## Main patient-level tumor analyses usually prefer primary tumor.
  ifelse(short_code == "TP", 1L,
    ifelse(short_code == "TR", 2L,
      ifelse(short_code == "NT", 3L, 9L)
    )
  )
}

analyte_priority <- function(analyte) {
  ifelse(analyte == "R", 1L, ifelse(analyte == "T", 2L, 9L))
}

first_non_na <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) NA else x[1]
}

make_pattern <- function(x) {
  paste(sort(unique(x)), collapse = "|")
}

write_csv_without_list_cols <- function(x, file) {
  list_cols <- names(x)[vapply(x, is.list, logical(1))]
  if (length(list_cols) > 0) {
    message("Dropping list columns for CSV: ", paste(list_cols, collapse = ", "))
    x <- x[setdiff(names(x), list_cols)]
  }
  write.csv(x, file, row.names = FALSE)
}

cat("Reading clinical data...\n")
cli <- readRDS(file.path("data_raw", "clinical_glioma.rds"))
gbm <- cli[cli$cohort == "GBM", , drop = FALSE]
barcode_parts <- parse_tcga_barcode(gbm$barcode)
gbm_audit <- cbind(gbm, barcode_parts[setdiff(names(barcode_parts), "barcode")])

cat("Clinical rows:", nrow(cli), "\n")
cat("GBM rows:", nrow(gbm_audit), "\n")
cat("GBM patients:", length(unique(gbm_audit$patient)), "\n\n")

cat("GBM sample type counts:\n")
print(table(gbm_audit$shortLetterCode, useNA = "ifany"))
cat("\nSample type code vs shortLetterCode:\n")
print(table(gbm_audit$sample_type_code, gbm_audit$shortLetterCode, useNA = "ifany"))
cat("\nAnalyte counts:\n")
print(table(gbm_audit$analyte, useNA = "ifany"))
cat("\nCenter counts:\n")
print(sort(table(gbm_audit$center), decreasing = TRUE))
cat("\nTop plate counts:\n")
print(head(sort(table(gbm_audit$plate), decreasing = TRUE), 30))

cat("\nReading TPM expression matrix for GBM sample QC...\n")
expr_tpm <- readRDS(file.path("data_annotated", "expr_tpm_glioma_anno.rds"))
sample_cols <- setdiff(colnames(expr_tpm), "gene_type")
missing_in_expr <- setdiff(gbm_audit$barcode, sample_cols)
extra_expr <- setdiff(sample_cols, cli$barcode)

cat("Expression sample columns:", length(sample_cols), "\n")
cat("GBM barcodes missing in expression:", length(missing_in_expr), "\n")
if (length(missing_in_expr) > 0) print(missing_in_expr)
cat("Expression samples not in clinical:", length(extra_expr), "\n")
if (length(extra_expr) > 0) print(head(extra_expr, 20))

gbm_barcodes <- intersect(gbm_audit$barcode, sample_cols)
expr_gbm <- as.matrix(expr_tpm[, gbm_barcodes, drop = FALSE])
storage.mode(expr_gbm) <- "double"
gene_type <- expr_tpm$gene_type
protein_coding <- !is.na(gene_type) & gene_type == "protein_coding"

qc <- data.frame(
  barcode = gbm_barcodes,
  total_tpm = as.numeric(colSums(expr_gbm, na.rm = TRUE)),
  detected_tpm_gt0 = as.integer(colSums(expr_gbm > 0, na.rm = TRUE)),
  detected_tpm_gt1 = as.integer(colSums(expr_gbm > 1, na.rm = TRUE)),
  pc_detected_tpm_gt0 = as.integer(colSums(expr_gbm[protein_coding, , drop = FALSE] > 0, na.rm = TRUE)),
  pc_detected_tpm_gt1 = as.integer(colSums(expr_gbm[protein_coding, , drop = FALSE] > 1, na.rm = TRUE)),
  stringsAsFactors = FALSE
)
if ("LAP3" %in% rownames(expr_tpm)) {
  qc$LAP3_tpm <- as.numeric(expr_gbm["LAP3", qc$barcode])
}

gbm_audit <- merge(gbm_audit, qc, by = "barcode", all.x = TRUE, sort = FALSE)
gbm_audit$patient_n_barcode <- ave(gbm_audit$barcode, gbm_audit$patient, FUN = length)
gbm_audit$sample_priority <- sample_type_priority(gbm_audit$shortLetterCode)
gbm_audit$analyte_priority <- analyte_priority(gbm_audit$analyte)

## QC-driven preliminary candidate. This is evidence, not a final biological rule.
ord <- order(
  gbm_audit$patient,
  gbm_audit$sample_priority,
  gbm_audit$analyte_priority,
  -gbm_audit$pc_detected_tpm_gt1,
  -gbm_audit$detected_tpm_gt1,
  gbm_audit$vial,
  gbm_audit$portion,
  gbm_audit$plate,
  gbm_audit$center,
  na.last = TRUE
)
gbm_ranked <- gbm_audit[ord, ]
gbm_ranked$retention_rank_qc <- ave(
  seq_len(nrow(gbm_ranked)),
  gbm_ranked$patient,
  FUN = seq_along
)
gbm_audit <- merge(
  gbm_audit,
  gbm_ranked[, c("barcode", "retention_rank_qc")],
  by = "barcode",
  all.x = TRUE,
  sort = FALSE
)
gbm_audit$keep_candidate_qc <- gbm_audit$retention_rank_qc == 1L
gbm_audit$retention_note <- ifelse(
  gbm_audit$shortLetterCode == "NT",
  "normal_sample_not_for_main_tumor_analysis",
  ifelse(
    gbm_audit$shortLetterCode == "TR" &
      ave(gbm_audit$shortLetterCode == "TP", gbm_audit$patient, FUN = any),
    "recurrent_sample_when_primary_exists",
    ifelse(gbm_audit$keep_candidate_qc, "candidate_by_type_analyte_expression_qc", "lower_rank_within_patient")
  )
)

patients <- split(gbm_audit, gbm_audit$patient)
patient_summary <- do.call(rbind, lapply(patients, function(d) {
  data.frame(
    patient = d$patient[1],
    n_barcode = nrow(d),
    n_tp = sum(d$shortLetterCode == "TP", na.rm = TRUE),
    n_tr = sum(d$shortLetterCode == "TR", na.rm = TRUE),
    n_nt = sum(d$shortLetterCode == "NT", na.rm = TRUE),
    has_tp = any(d$shortLetterCode == "TP", na.rm = TRUE),
    has_tr = any(d$shortLetterCode == "TR", na.rm = TRUE),
    has_nt = any(d$shortLetterCode == "NT", na.rm = TRUE),
    sample_type_pattern = make_pattern(d$shortLetterCode),
    sample_code_pattern = make_pattern(d$sample_type_code),
    vial_pattern = make_pattern(d$vial),
    portion_pattern = make_pattern(d$portion),
    analyte_pattern = make_pattern(d$analyte),
    plate_pattern = make_pattern(d$plate),
    center_pattern = make_pattern(d$center),
    keep_barcode_qc = first_non_na(d$barcode[d$keep_candidate_qc]),
    keep_shortLetterCode_qc = first_non_na(d$shortLetterCode[d$keep_candidate_qc]),
    keep_plate_qc = first_non_na(d$plate[d$keep_candidate_qc]),
    keep_center_qc = first_non_na(d$center[d$keep_candidate_qc]),
    stringsAsFactors = FALSE
  )
}))
row.names(patient_summary) <- NULL
patient_summary <- patient_summary[order(-patient_summary$n_barcode, patient_summary$patient), ]

patient_summary$duplicate_class <- with(patient_summary, ifelse(
  n_barcode == 1 & n_tp == 1, "single_TP",
  ifelse(n_barcode == 1 & n_tr == 1, "single_TR_only",
    ifelse(n_barcode == 1 & n_nt == 1, "single_NT_only",
      ifelse(n_tp > 1 & n_tr == 0 & n_nt == 0, "multi_TP_only",
        ifelse(n_tp >= 1 & n_tr >= 1, "TP_plus_TR",
          ifelse(n_tp == 0 & n_tr > 1, "multi_TR_only",
            ifelse(n_tp == 0 & n_nt > 0, "NT_only", "other")
          )
        )
      )
    )
  )
))

dup_detail <- gbm_audit[gbm_audit$patient_n_barcode > 1, ]
dup_detail <- dup_detail[order(
  dup_detail$patient,
  dup_detail$sample_priority,
  dup_detail$vial,
  dup_detail$portion,
  dup_detail$analyte,
  dup_detail$plate,
  dup_detail$center
), ]

cat("\nPatient duplicate class counts:\n")
print(sort(table(patient_summary$duplicate_class), decreasing = TRUE))

cat("\nPatient n_barcode distribution:\n")
print(table(patient_summary$n_barcode))

cat("\nPatients with only recurrent tumor samples:\n")
print(patient_summary$patient[patient_summary$n_tp == 0 & patient_summary$n_tr > 0])

cat("\nPatients with only normal samples:\n")
print(patient_summary$patient[patient_summary$n_tp == 0 & patient_summary$n_tr == 0 & patient_summary$n_nt > 0])

cat("\nTop duplicate patient summary:\n")
print(head(patient_summary, 30), row.names = FALSE)

cat("\nComputing within-patient pairwise expression correlations for duplicate patients...\n")
pc_expr_log <- log2(expr_gbm[protein_coding, , drop = FALSE] + 1)
pair_rows <- list()
pair_i <- 0L
dup_patients <- names(patients)[vapply(patients, nrow, integer(1)) > 1]

for (patient_id in dup_patients) {
  d <- patients[[patient_id]]
  b <- intersect(d$barcode, colnames(pc_expr_log))
  if (length(b) < 2) next
  combos <- combn(b, 2, simplify = FALSE)
  for (combo in combos) {
    d1 <- d[d$barcode == combo[1], ][1, ]
    d2 <- d[d$barcode == combo[2], ][1, ]
    pair_i <- pair_i + 1L
    pair_rows[[pair_i]] <- data.frame(
      patient = patient_id,
      barcode_1 = combo[1],
      barcode_2 = combo[2],
      short_1 = d1$shortLetterCode,
      short_2 = d2$shortLetterCode,
      sample_code_1 = d1$sample_type_code,
      sample_code_2 = d2$sample_type_code,
      vial_1 = d1$vial,
      vial_2 = d2$vial,
      portion_1 = d1$portion,
      portion_2 = d2$portion,
      analyte_1 = d1$analyte,
      analyte_2 = d2$analyte,
      plate_1 = d1$plate,
      plate_2 = d2$plate,
      center_1 = d1$center,
      center_2 = d2$center,
      pearson_log2_tpm_pc = suppressWarnings(cor(
        pc_expr_log[, combo[1]],
        pc_expr_log[, combo[2]],
        use = "pairwise.complete.obs",
        method = "pearson"
      )),
      stringsAsFactors = FALSE
    )
  }
}

pairwise_cor <- if (length(pair_rows) > 0) do.call(rbind, pair_rows) else data.frame()
if (nrow(pairwise_cor) > 0) {
  pairwise_cor$same_short <- pairwise_cor$short_1 == pairwise_cor$short_2
  pairwise_cor$same_sample_code <- pairwise_cor$sample_code_1 == pairwise_cor$sample_code_2
  pairwise_cor$same_vial <- pairwise_cor$vial_1 == pairwise_cor$vial_2
  pairwise_cor$same_portion <- pairwise_cor$portion_1 == pairwise_cor$portion_2
  pairwise_cor$same_analyte <- pairwise_cor$analyte_1 == pairwise_cor$analyte_2
  pairwise_cor$same_center <- pairwise_cor$center_1 == pairwise_cor$center_2
  pairwise_cor <- pairwise_cor[order(pairwise_cor$patient, -pairwise_cor$pearson_log2_tpm_pc), ]
}

cat("\nPairwise duplicate correlation summary:\n")
if (nrow(pairwise_cor) > 0) {
  print(summary(pairwise_cor$pearson_log2_tpm_pc))
  cat("\nCorrelation by shortLetterCode pair:\n")
  pairwise_cor$short_pair <- paste(pairwise_cor$short_1, pairwise_cor$short_2, sep = "_")
  print(aggregate(pearson_log2_tpm_pc ~ short_pair, pairwise_cor, summary))
  cat("\nCorrelation by center pair:\n")
  pairwise_cor$center_pair <- paste(pairwise_cor$center_1, pairwise_cor$center_2, sep = "_")
  print(aggregate(pearson_log2_tpm_pc ~ center_pair, pairwise_cor, summary))
} else {
  cat("No duplicate pair correlations were computed.\n")
}

candidate <- gbm_audit[gbm_audit$keep_candidate_qc, ]
candidate <- candidate[order(candidate$patient), ]

output_files <- c(
  audit = file.path(log_dir, paste0("gbm_barcode_audit_", run_id, ".csv")),
  audit_rds = file.path(log_dir, paste0("gbm_barcode_audit_", run_id, ".rds")),
  patient_summary = file.path(log_dir, paste0("gbm_patient_summary_", run_id, ".csv")),
  duplicate_detail = file.path(log_dir, paste0("gbm_duplicate_detail_", run_id, ".csv")),
  duplicate_detail_rds = file.path(log_dir, paste0("gbm_duplicate_detail_", run_id, ".rds")),
  pairwise_cor = file.path(log_dir, paste0("gbm_duplicate_pairwise_cor_", run_id, ".csv")),
  candidate = file.path(log_dir, paste0("gbm_retention_candidate_qc_", run_id, ".csv"))
)

write_csv_without_list_cols(gbm_audit, output_files[["audit"]])
saveRDS(gbm_audit, output_files[["audit_rds"]])
write.csv(patient_summary, output_files[["patient_summary"]], row.names = FALSE)
write_csv_without_list_cols(dup_detail, output_files[["duplicate_detail"]])
saveRDS(dup_detail, output_files[["duplicate_detail_rds"]])
write.csv(pairwise_cor, output_files[["pairwise_cor"]], row.names = FALSE)
write_csv_without_list_cols(candidate, output_files[["candidate"]])

cat("\nOutput files:\n")
print(output_files)
cat("\nDone.\n")
