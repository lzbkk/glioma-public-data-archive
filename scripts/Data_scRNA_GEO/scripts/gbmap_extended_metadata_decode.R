#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)

input_file <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(rhdf5)
})
data.table::setDTthreads(8)

read_categorical_path <- function(file, path) {
  categories <- h5read(file, paste0(path, "/categories"))
  codes <- as.integer(h5read(file, paste0(path, "/codes")))
  values <- rep(NA_character_, length(codes))
  keep <- codes >= 0L
  values[keep] <- categories[codes[keep] + 1L]
  values
}

read_categorical <- function(file, field) {
  read_categorical_path(file, paste0("/obs/", field))
}

h5_index <- h5ls(input_file, recursive = TRUE, all = TRUE)
obs_index <- h5_index[h5_index$group == "/obs" | startsWith(h5_index$group, "/obs/"), ]
obs_top <- unique(c(
  obs_index$name[obs_index$group == "/obs"],
  sub("^/obs/([^/]+).*$", "\\1", obs_index$group[obs_index$group != "/obs"])
))
obs_top <- sort(obs_top[nzchar(obs_top)])

candidate_fields <- c(
  "author", "donor_id", "sample", "assay", "suspension_type",
  "annotation_level_1", "annotation_level_2", "annotation_level_3",
  "annotation_level_4", "cell_type", "celltype_original", "iCNV"
)
fields <- intersect(candidate_fields, obs_top)
message("Decoding fields: ", paste(fields, collapse = ", "))
obs <- as.data.table(setNames(lapply(fields, function(field) {
  read_categorical(input_file, field)
}), fields))

expected_cells <- 1135677L
if (nrow(obs) != expected_cells) {
  warning("Expected ", expected_cells, " cells, observed ", nrow(obs))
}

summary_by_author <- obs[, .(
  cells = .N,
  donors = if ("donor_id" %in% names(obs)) uniqueN(donor_id) else NA_integer_,
  samples = if ("sample" %in% names(obs)) uniqueN(sample) else NA_integer_,
  assays = if ("assay" %in% names(obs)) paste(sort(unique(assay)), collapse = ";") else NA_character_,
  suspensions = if ("suspension_type" %in% names(obs)) paste(sort(unique(suspension_type)), collapse = ";") else NA_character_,
  level1 = if ("annotation_level_1" %in% names(obs)) paste(sort(unique(annotation_level_1)), collapse = ";") else NA_character_,
  icnv_values = if ("iCNV" %in% names(obs)) paste(sort(unique(na.omit(iCNV))), collapse = ";") else NA_character_
), by = author][order(-cells)]
fwrite(summary_by_author, file.path(output_dir, "extended_gbmap_author_summary.csv"))

if (all(c("author", "annotation_level_1") %in% names(obs))) {
  fwrite(
    obs[, .(cells = .N, donors = uniqueN(donor_id)), by = .(author, annotation_level_1)][order(author, -cells)],
    file.path(output_dir, "extended_gbmap_author_level1_summary.csv")
  )
}

if (all(c("author", "annotation_level_4") %in% names(obs))) {
  fwrite(
    obs[, .(cells = .N, donors = uniqueN(donor_id)), by = .(author, annotation_level_4)][order(author, -cells)],
    file.path(output_dir, "extended_gbmap_author_level4_summary.csv")
  )
}

if (all(c("annotation_level_1", "annotation_level_2", "annotation_level_3", "annotation_level_4") %in% names(obs))) {
  hierarchy <- obs[, .(
    cells = .N,
    donors = uniqueN(donor_id),
    authors = uniqueN(author)
  ), by = .(annotation_level_1, annotation_level_2, annotation_level_3, annotation_level_4)][order(annotation_level_1, -cells)]
  fwrite(hierarchy, file.path(output_dir, "extended_gbmap_annotation_hierarchy_summary.csv"))
}

field_levels <- rbindlist(lapply(fields, function(field) {
  values <- obs[[field]]
  levels <- sort(unique(na.omit(values)))
  data.table(
    field = field,
    level = levels,
    cells = vapply(levels, function(level) sum(values == level, na.rm = TRUE), integer(1))
  )
}), fill = TRUE)
fwrite(field_levels, file.path(output_dir, "extended_gbmap_selected_field_levels.csv"))

entry_count_list <- list(
  data.table(population = "full_extended", cells = nrow(obs), donors = uniqueN(obs$donor_id), authors = uniqueN(obs$author)),
  obs[author != "Neftel2019", .(population = "exclude_Neftel2019", cells = .N, donors = uniqueN(donor_id), authors = uniqueN(author))],
  obs[author != "Neftel2019" & annotation_level_1 == "Neoplastic",
      .(population = "exclude_Neftel2019_neoplastic", cells = .N, donors = uniqueN(donor_id), authors = uniqueN(author))]
)
if ("iCNV" %in% names(obs)) {
  entry_count_list[[length(entry_count_list) + 1L]] <- obs[
    author != "Neftel2019" & annotation_level_1 == "Neoplastic" & iCNV == "aneuploid",
    .(population = "exclude_Neftel2019_neoplastic_aneuploid", cells = .N, donors = uniqueN(donor_id), authors = uniqueN(author))
  ]
} else {
  entry_count_list[[length(entry_count_list) + 1L]] <- data.table(
    population = "exclude_Neftel2019_neoplastic_aneuploid",
    cells = NA_integer_,
    donors = NA_integer_,
    authors = NA_integer_,
    note = "iCNV field unavailable in Extended object"
  )
}
entry_counts <- rbindlist(entry_count_list, fill = TRUE)
fwrite(entry_counts, file.path(output_dir, "extended_gbmap_analysis_entry_counts.csv"))

immune_levels <- if (all(c("annotation_level_1", "annotation_level_2", "annotation_level_3") %in% names(obs))) {
  by_cols <- c("annotation_level_2", "annotation_level_3")
  if ("annotation_level_4" %in% names(obs)) by_cols <- c(by_cols, "annotation_level_4")
  obs[annotation_level_1 == "Non-neoplastic", .(
    cells = .N,
    donors = uniqueN(donor_id),
    authors = uniqueN(author)
  ), by = by_cols][order(-cells)]
} else {
  data.table()
}
fwrite(immune_levels, file.path(output_dir, "extended_gbmap_non_neoplastic_immune_niche_summary.csv"))

summary_lines <- c(
  "# Extended GBmap Decoded Metadata Audit",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  paste0("- Cells: ", format(nrow(obs), big.mark = ",")),
  paste0("- Authors/studies: ", uniqueN(obs$author)),
  paste0("- Donors: ", uniqueN(obs$donor_id)),
  paste0("- Samples: ", if ("sample" %in% names(obs)) uniqueN(obs$sample) else "NA"),
  paste0("- Assays: ", if ("assay" %in% names(obs)) paste(sort(unique(obs$assay)), collapse = "; ") else "NA"),
  paste0("- Suspension types: ", if ("suspension_type" %in% names(obs)) paste(sort(unique(obs$suspension_type)), collapse = "; ") else "NA"),
  "",
  "## Entry counts",
  "",
  paste(capture.output(print(entry_counts)), collapse = "\n"),
  "",
  "## Use",
  "",
  "Extended GBmap is appropriate as a second-stage atlas asset for rare immune/niche,",
  "cell-cell communication, deconvolution, and broader state-signature work. It should",
  "not be merged into Core-first LAP3-mTORC1 inference until source overlap and",
  "author/suspension/assay sensitivity rules are frozen."
)
writeLines(summary_lines, file.path(output_dir, "README_Decoded_Metadata.md"))

message("Extended metadata decode complete")
