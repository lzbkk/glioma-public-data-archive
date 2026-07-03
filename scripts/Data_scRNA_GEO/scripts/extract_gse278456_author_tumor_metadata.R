#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(data.table)
  library(dplyr)
})

setwd("/home/lzb/glioma")
data.table::setDTthreads(8)

author_file <- "Data_scRNA_GEO/GSE278456/GSE278456_TumorSeurat.RDS"
local_file <- "Data_scRNA_GEO/GSE278456/Tumor_Integrated_SeuratV5.rds"
out_dir <- "Data_scRNA_GEO/results/GSE278456_Tumor_Object_Audit"
table_dir <- file.path(out_dir, "tables")
export_dir <- file.path(out_dir, "exports")
log_dir <- file.path(out_dir, "logs")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(log_dir, "extract_gse278456_author_tumor_metadata.log")
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message")
on.exit({
  sink(type = "message")
  sink(type = "output")
  close(log_con)
}, add = TRUE)

cat("Started:", format(Sys.time()), "\n")
cat("R:", R.version.string, "\n")
cat("Seurat:", as.character(packageVersion("Seurat")), "\n")
cat("Available memory before load:\n")
print(system("free -h", intern = TRUE))

extract_barcode <- function(cell_id) {
  match <- regexpr("[ACGT]{16}", cell_id)
  out <- regmatches(cell_id, match)
  out[match < 0] <- NA_character_
  out
}

cat("\nLoading author object:", author_file, "\n")
load_start <- Sys.time()
author_obj <- readRDS(author_file)
cat("Author ReadRDS minutes:", round(difftime(Sys.time(), load_start, units = "mins"), 2), "\n")
cat("Author class:", paste(class(author_obj), collapse = ";"), "\n")
cat("Author version:", as.character(author_obj@version), "\n")
cat("Author dimensions:", nrow(author_obj), "features x", ncol(author_obj), "cells\n")
cat("Author object size GB:", round(as.numeric(object.size(author_obj)) / 1024^3, 2), "\n")

author_meta <- author_obj[[]]
author_meta$cell_id_author <- rownames(author_meta)
rownames(author_meta) <- NULL

patient_field <- intersect(c("Patient", "Patient_ID", "patient", "sample"), names(author_meta))
if (length(patient_field) == 0) {
  stop("No patient field found in author metadata")
}
patient_field <- patient_field[1]
author_meta$patient_key <- as.character(author_meta[[patient_field]])
author_meta$cell_barcode <- extract_barcode(author_meta$cell_id_author)
author_meta$universal_key <- paste(author_meta$cell_barcode, author_meta$patient_key, sep = "||")

author_field_audit <- data.frame(
  field = names(author_meta),
  class = vapply(author_meta, function(x) paste(class(x), collapse = ";"), character(1)),
  n_non_missing = vapply(author_meta, function(x) sum(!is.na(x) & as.character(x) != ""), integer(1)),
  n_unique = vapply(author_meta, function(x) length(unique(x[!is.na(x)])), integer(1)),
  example_values = vapply(author_meta, function(x) {
    vals <- unique(as.character(x[!is.na(x) & as.character(x) != ""]))
    paste(head(vals, 8), collapse = ";")
  }, character(1))
)
fwrite(author_field_audit, file.path(table_dir, "gse278456_author_metadata_field_audit.csv"))

author_assay_rows <- list()
for (assay_name in Assays(author_obj)) {
  assay <- author_obj[[assay_name]]
  for (layer_name in Layers(assay)) {
    layer <- LayerData(assay, layer = layer_name, fast = TRUE)
    author_assay_rows[[length(author_assay_rows) + 1L]] <- data.frame(
      assay = assay_name,
      layer = layer_name,
      n_features = nrow(layer),
      n_cells = ncol(layer),
      matrix_class = paste(class(layer), collapse = ";")
    )
  }
}
fwrite(
  bind_rows(author_assay_rows),
  file.path(table_dir, "gse278456_author_assay_layer_audit.csv")
)

fwrite(
  author_meta,
  file.path(export_dir, "gse278456_author_tumor_metadata.csv.gz"),
  compress = "gzip"
)
saveRDS(
  author_meta,
  file.path(export_dir, "gse278456_author_tumor_metadata.rds"),
  compress = "xz"
)

cat("\nAuthor metadata fields:\n")
print(author_field_audit)
cat("\nPatient field:", patient_field, "\n")
cat("Author universal-key duplicates:", sum(duplicated(author_meta$universal_key)), "\n")

author_keys <- author_meta$universal_key
author_meta_small <- author_meta
rm(author_obj, author_meta)
invisible(gc())
cat("Released author object. Available memory:\n")
print(system("free -h", intern = TRUE))

cat("\nLoading local integrated object for crosswalk...\n")
local_start <- Sys.time()
local_obj <- readRDS(local_file)
cat("Local ReadRDS minutes:", round(difftime(Sys.time(), local_start, units = "mins"), 2), "\n")
local_meta <- local_obj[[]]
local_meta$cell_id_local <- rownames(local_meta)
rownames(local_meta) <- NULL
if (!"Patient_ID" %in% names(local_meta)) {
  stop("Patient_ID absent from local object")
}
local_meta$cell_barcode <- extract_barcode(local_meta$cell_id_local)
local_meta$patient_key <- as.character(local_meta$Patient_ID)
local_meta$universal_key <- paste(local_meta$cell_barcode, local_meta$patient_key, sep = "||")

match_idx <- match(local_meta$universal_key, author_keys)
matched <- !is.na(match_idx)

transfer_candidates <- intersect(
  c(
    "Patient", "Phase", "Tumor.Grade", "Pathology", "IDH.status",
    "EGFR.amplification", "EGFR.mutation", "PTEN.mutation",
    "TP53.mutation", "CDKN2A..2G.loss", "Tumor_State",
    "Cell.Type", "cell_type", "celltype", "seurat_clusters"
  ),
  names(author_meta_small)
)

crosswalk <- local_meta %>%
  transmute(
    cell_id_local,
    Patient_ID,
    File_Source,
    cell_barcode,
    universal_key,
    matched_author = matched,
    cell_id_author = ifelse(matched, author_meta_small$cell_id_author[match_idx], NA_character_)
  )
for (field in transfer_candidates) {
  crosswalk[[paste0("author_", field)]] <- author_meta_small[[field]][match_idx]
}

fwrite(
  crosswalk,
  file.path(export_dir, "gse278456_local_author_metadata_crosswalk.csv.gz"),
  compress = "gzip"
)
saveRDS(
  crosswalk,
  file.path(export_dir, "gse278456_local_author_metadata_crosswalk.rds"),
  compress = "xz"
)

match_summary <- data.frame(
  metric = c(
    "author_cells", "local_cells", "matched_local_cells",
    "match_rate", "author_key_duplicates", "local_key_duplicates",
    "author_metadata_fields", "transfer_candidate_fields"
  ),
  value = c(
    nrow(author_meta_small),
    nrow(local_meta),
    sum(matched),
    mean(matched),
    sum(duplicated(author_keys)),
    sum(duplicated(local_meta$universal_key)),
    ncol(author_meta_small),
    length(transfer_candidates)
  ),
  detail = c(
    "", "", "", "", "", "",
    paste(names(author_meta_small), collapse = ";"),
    paste(transfer_candidates, collapse = ";")
  )
)
fwrite(match_summary, file.path(table_dir, "gse278456_author_local_match_summary.csv"))

field_match_summary <- bind_rows(lapply(transfer_candidates, function(field) {
  values <- author_meta_small[[field]][match_idx]
  data.frame(
    field = field,
    matched_non_missing = sum(matched & !is.na(values) & as.character(values) != ""),
    matched_unique = length(unique(values[matched & !is.na(values)])),
    example_values = paste(head(unique(as.character(values[matched & !is.na(values)])), 12), collapse = ";")
  )
}))
fwrite(
  field_match_summary,
  file.path(table_dir, "gse278456_transfer_field_match_summary.csv")
)

cat("\nCrosswalk summary:\n")
print(match_summary)
cat("\nTransfer fields:\n")
print(field_match_summary)
cat("Completed:", format(Sys.time()), "\n")
