#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(data.table)
  library(dplyr)
})

setwd("/home/lzb/glioma")

object_file <- "Data_scRNA_GEO/GSE278456/Tumor_Integrated_SeuratV5.rds"
out_dir <- "Data_scRNA_GEO/results/GSE278456_Tumor_Object_Audit"
table_dir <- file.path(out_dir, "tables")
log_dir <- file.path(out_dir, "logs")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(log_dir, "audit_gse278456_tumor_object.log")
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
cat("SeuratObject:", as.character(packageVersion("SeuratObject")), "\n")
cat("Object file:", object_file, "\n")
cat("Object bytes:", file.info(object_file)$size, "\n")

load_start <- Sys.time()
obj <- readRDS(object_file)
cat("ReadRDS minutes:", round(difftime(Sys.time(), load_start, units = "mins"), 2), "\n")
cat("Class:", paste(class(obj), collapse = ";"), "\n")
cat("Object version:", as.character(obj@version), "\n")
cat("Object size GB:", round(as.numeric(object.size(obj)) / 1024^3, 2), "\n")
cat("Dimensions:", nrow(obj), "features x", ncol(obj), "cells\n")

metadata <- obj[[]]
metadata$cell_id <- rownames(metadata)
meta_summary <- data.frame(
  field = names(metadata),
  class = vapply(metadata, function(x) paste(class(x), collapse = ";"), character(1)),
  n_non_missing = vapply(metadata, function(x) sum(!is.na(x) & as.character(x) != ""), integer(1)),
  n_unique = vapply(metadata, function(x) length(unique(x[!is.na(x)])), integer(1)),
  example_values = vapply(metadata, function(x) {
    vals <- unique(as.character(x[!is.na(x) & as.character(x) != ""]))
    paste(head(vals, 8), collapse = ";")
  }, character(1))
)
fwrite(meta_summary, file.path(table_dir, "gse278456_tumor_metadata_field_audit.csv"))

assay_rows <- list()
for (assay_name in Assays(obj)) {
  assay <- obj[[assay_name]]
  assay_layers <- Layers(assay)
  for (layer_name in assay_layers) {
    layer <- LayerData(assay, layer = layer_name, fast = TRUE)
    assay_rows[[length(assay_rows) + 1L]] <- data.frame(
      assay = assay_name,
      layer = layer_name,
      n_features = nrow(layer),
      n_cells = ncol(layer),
      matrix_class = paste(class(layer), collapse = ";")
    )
  }
}
assay_audit <- bind_rows(assay_rows)
fwrite(assay_audit, file.path(table_dir, "gse278456_tumor_assay_layer_audit.csv"))

reduction_audit <- bind_rows(lapply(Reductions(obj), function(reduction_name) {
  embedding <- Embeddings(obj, reduction = reduction_name)
  data.frame(
    reduction = reduction_name,
    n_cells = nrow(embedding),
    n_dimensions = ncol(embedding),
    key = Key(obj[[reduction_name]])
  )
}))
fwrite(reduction_audit, file.path(table_dir, "gse278456_tumor_reduction_audit.csv"))

command_audit <- bind_rows(lapply(names(obj@commands), function(command_name) {
  command <- obj@commands[[command_name]]
  data.frame(
    command = command_name,
    assay_used = paste(command@assay.used, collapse = ";"),
    call_string = paste(deparse(command@call.string), collapse = " "),
    time_stamp = as.character(command@time.stamp)
  )
}))
fwrite(command_audit, file.path(table_dir, "gse278456_tumor_command_audit.csv"))

candidate_fields <- intersect(
  c(
    "Patient_ID", "Patient", "File_Source", "orig.ident", "seurat_clusters",
    "Tumor_State", "Phase", "Tumor.Grade", "Pathology", "IDH.status",
    "EGFR.amplification", "EGFR.mutation", "PTEN.mutation", "TP53.mutation",
    "CDKN2A..2G.loss"
  ),
  names(metadata)
)

metadata_counts <- bind_rows(lapply(candidate_fields, function(field) {
  metadata %>%
    transmute(value = as.character(.data[[field]])) %>%
    mutate(value = ifelse(is.na(value) | value == "", "NA", value)) %>%
    count(value, name = "n_cells") %>%
    mutate(field = field, .before = 1)
}))
fwrite(metadata_counts, file.path(table_dir, "gse278456_tumor_metadata_value_counts.csv"))

rna_features <- rownames(obj[["RNA"]])
gene_check <- data.frame(
  gene = c(
    "LAP3", "BCAT1", "BCAT2", "SLC7A5", "SLC3A2", "MTOR", "RPTOR",
    "RPS6KB1", "EIF4EBP1", "RPS6", "CHI3L1", "CD44", "GFAP", "AQP4",
    "PDGFRA", "OLIG2", "SOX4", "SOX11"
  ),
  present = c(
    "LAP3", "BCAT1", "BCAT2", "SLC7A5", "SLC3A2", "MTOR", "RPTOR",
    "RPS6KB1", "EIF4EBP1", "RPS6", "CHI3L1", "CD44", "GFAP", "AQP4",
    "PDGFRA", "OLIG2", "SOX4", "SOX11"
  ) %in% rna_features
)
fwrite(gene_check, file.path(table_dir, "gse278456_tumor_key_gene_coverage.csv"))

required_direct_fields <- c("Patient_ID")
state_fields <- intersect(c("Tumor_State", "Phase"), names(metadata))
counts_available <- "counts" %in% Layers(obj[["RNA"]])
data_available <- "data" %in% Layers(obj[["RNA"]])
patient_available <- all(required_direct_fields %in% names(metadata))
state_available <- length(state_fields) > 0

readiness <- data.frame(
  check = c(
    "patient_field",
    "author_or_validated_tumor_state",
    "RNA_counts_layer",
    "RNA_data_layer",
    "LAP3_present",
    "clinical_annotation",
    "direct_patient_aware_analysis_ready"
  ),
  pass = c(
    patient_available,
    state_available,
    counts_available,
    data_available,
    "LAP3" %in% rna_features,
    any(c("IDH.status", "Tumor.Grade", "Pathology") %in% names(metadata)),
    patient_available && state_available && counts_available && ("LAP3" %in% rna_features)
  ),
  detail = c(
    paste(intersect(required_direct_fields, names(metadata)), collapse = ";"),
    paste(state_fields, collapse = ";"),
    paste(Layers(obj[["RNA"]]), collapse = ";"),
    paste(Layers(obj[["RNA"]]), collapse = ";"),
    ifelse("LAP3" %in% rna_features, "LAP3 found", "LAP3 absent"),
    paste(intersect(c("IDH.status", "Tumor.Grade", "Pathology"), names(metadata)), collapse = ";"),
    "Requires patient field, validated state, counts layer, and LAP3"
  )
)
fwrite(readiness, file.path(table_dir, "gse278456_tumor_analysis_readiness.csv"))

cat("\nMetadata fields:\n")
print(meta_summary)
cat("\nAssays/layers:\n")
print(assay_audit)
cat("\nReductions:\n")
print(reduction_audit)
cat("\nReadiness:\n")
print(readiness)
cat("\nCandidate metadata counts:\n")
print(metadata_counts)
cat("Completed:", format(Sys.time()), "\n")
