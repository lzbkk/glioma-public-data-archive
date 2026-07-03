#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)

input_file <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(rhdf5)
})

message("Reading Extended GBmap H5AD structure without loading expression matrix")
h5_index <- h5ls(input_file, recursive = TRUE, all = TRUE)
write.csv(
  h5_index,
  file.path(output_dir, "extended_gbmap_h5_structure.csv"),
  row.names = FALSE
)

obs_index <- h5_index[
  h5_index$group == "/obs" | startsWith(h5_index$group, "/obs/"),
  ,
  drop = FALSE
]
write.csv(
  obs_index,
  file.path(output_dir, "extended_gbmap_obs_structure.csv"),
  row.names = FALSE
)

obs_top_level <- unique(c(
  obs_index$name[obs_index$group == "/obs"],
  sub("^/obs/([^/]+).*$", "\\1", obs_index$group[obs_index$group != "/obs"])
))
obs_top_level <- sort(obs_top_level[nzchar(obs_top_level)])

required_patterns <- c(
  "author", "study", "source", "dataset", "donor", "patient", "sample",
  "cell_type", "annotation", "malignant", "neoplastic", "state", "assay",
  "suspension", "disease", "tissue", "iCNV"
)
candidate_fields <- obs_top_level[
  grepl(paste(required_patterns, collapse = "|"), obs_top_level, ignore.case = TRUE)
]

field_table <- data.frame(
  field = obs_top_level,
  candidate_for_analysis = obs_top_level %in% candidate_fields,
  stringsAsFactors = FALSE
)
write.csv(
  field_table,
  file.path(output_dir, "extended_gbmap_obs_fields.csv"),
  row.names = FALSE
)

required_checks <- data.frame(
  requirement = c(
    "study_or_source",
    "patient_or_donor",
    "cell_type_or_annotation",
    "malignant_or_neoplastic",
    "state",
    "assay_or_platform",
    "cnv_or_malignancy_support"
  ),
  pattern = c(
    "author|study|source|dataset",
    "patient|donor",
    "cell_type|annotation",
    "malignant|neoplastic",
    "state",
    "assay|platform",
    "iCNV|cnv|malignant|neoplastic"
  ),
  stringsAsFactors = FALSE
)
required_checks$matched_fields <- vapply(
  required_checks$pattern,
  function(pattern) {
    hits <- obs_top_level[grepl(pattern, obs_top_level, ignore.case = TRUE)]
    paste(hits, collapse = ";")
  },
  character(1)
)
required_checks$available <- nzchar(required_checks$matched_fields)
write.csv(
  required_checks,
  file.path(output_dir, "extended_gbmap_entry_requirements.csv"),
  row.names = FALSE
)

summary_lines <- c(
  "# Extended GBmap H5AD Entry Audit",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Input",
  "",
  paste0("- File: `", input_file, "`"),
  paste0("- Size: ", file.info(input_file)$size, " bytes"),
  "",
  "## Metadata fields",
  "",
  paste0("- Total top-level `obs` fields: ", length(obs_top_level)),
  paste0("- Candidate source/patient/annotation fields: ",
         if (length(candidate_fields)) paste(candidate_fields, collapse = ", ") else "none"),
  "",
  "## Entry decision",
  "",
  if (all(required_checks$available[c(1, 2, 3)])) {
    "PASS for detailed metadata audit: source, patient/donor, and annotation-like fields are present."
  } else {
    "CONDITIONAL/FAIL: one or more source, patient/donor, or annotation field classes are absent."
  },
  "",
  "This structure audit does not load the expression matrix. Before using Extended",
  "GBmap for biology, decode categorical metadata, quantify author/donor/sample",
  "overlap with Core and GSE131928, and freeze a separate entry rule."
)
writeLines(summary_lines, file.path(output_dir, "README.md"))

message("Extended GBmap H5AD structure audit complete")
