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

fields <- c(
  "author", "donor_id", "sample", "assay", "annotation_level_1",
  "annotation_level_2", "annotation_level_3", "annotation_level_4",
  "celltype_original", "iCNV"
)
message("Decoding selected obs fields")
obs <- as.data.table(setNames(
  lapply(fields, function(field) read_categorical(input_file, field)),
  fields
))

stopifnot(nrow(obs) == 338564L)

author_summary <- obs[, .(
  cells = .N,
  donors = uniqueN(donor_id),
  samples = uniqueN(sample),
  assays = paste(sort(unique(assay)), collapse = ";"),
  icnv_values = paste(sort(unique(na.omit(iCNV))), collapse = ";")
), by = author][order(-cells)]
fwrite(author_summary, file.path(output_dir, "core_gbmap_author_summary.csv"))

author_icnv <- obs[, .(cells = .N, donors = uniqueN(donor_id)),
                   by = .(author, iCNV)][order(author, iCNV)]
fwrite(author_icnv, file.path(output_dir, "core_gbmap_author_icnv_summary.csv"))

author_level1 <- obs[, .(cells = .N, donors = uniqueN(donor_id)),
                     by = .(author, annotation_level_1)][order(author, -cells)]
fwrite(author_level1, file.path(output_dir, "core_gbmap_author_level1_summary.csv"))

author_level4 <- obs[, .(cells = .N, donors = uniqueN(donor_id)),
                     by = .(author, annotation_level_4)][order(author, -cells)]
fwrite(author_level4, file.path(output_dir, "core_gbmap_author_level4_summary.csv"))

donor_crosswalk <- unique(
  obs[, .(author, donor_id, sample, assay)]
)[order(author, donor_id, sample)]
fwrite(donor_crosswalk, file.path(output_dir, "core_gbmap_donor_sample_crosswalk.csv"))

field_levels <- rbindlist(lapply(fields, function(field) {
  values <- obs[[field]]
  levels <- sort(unique(na.omit(values)))
  out <- data.table(
    field = field,
    level = levels,
    cells = vapply(levels, function(level) sum(values == level, na.rm = TRUE), integer(1))
  )
  if (anyNA(values)) {
    out <- rbind(out, data.table(field = field, level = NA_character_,
                                 cells = sum(is.na(values))))
  }
  out
}), use.names = TRUE)
fwrite(field_levels, file.path(output_dir, "core_gbmap_selected_field_levels.csv"))

lap3_checks <- list()
for (var_group in c("/var", "/raw/var")) {
  ids <- h5read(input_file, paste0(var_group, "/_index"))
  names <- read_categorical_path(input_file, paste0(var_group, "/feature_name"))
  lap3_checks[[var_group]] <- data.table(
    var_group = var_group,
    feature_id = ids,
    feature_name = names
  )[feature_name == "LAP3" | grepl("^LAP3([._-]|$)", feature_id)]
}
lap3_table <- rbindlist(lap3_checks, fill = TRUE)
fwrite(lap3_table, file.path(output_dir, "core_gbmap_lap3_feature_check.csv"))

gse131928_candidates <- author_summary[
  grepl("Neftel|131928", author, ignore.case = TRUE)
]
fwrite(
  gse131928_candidates,
  file.path(output_dir, "core_gbmap_gse131928_source_candidates.csv")
)

analysis_entry <- rbindlist(list(
  data.table(
    population = "full_core",
    cells = nrow(obs),
    donors = uniqueN(obs$donor_id),
    authors = uniqueN(obs$author)
  ),
  obs[author != "Neftel2019", .(
    population = "exclude_Neftel2019",
    cells = .N,
    donors = uniqueN(donor_id),
    authors = uniqueN(author)
  )],
  obs[author != "Neftel2019" & annotation_level_1 == "Neoplastic", .(
    population = "exclude_Neftel2019_neoplastic",
    cells = .N,
    donors = uniqueN(donor_id),
    authors = uniqueN(author)
  )],
  obs[author != "Neftel2019" & annotation_level_1 == "Neoplastic" &
        iCNV == "aneuploid", .(
    population = "exclude_Neftel2019_neoplastic_aneuploid",
    cells = .N,
    donors = uniqueN(donor_id),
    authors = uniqueN(author)
  )]
), fill = TRUE)
fwrite(analysis_entry, file.path(output_dir, "core_gbmap_analysis_entry_counts.csv"))

summary_lines <- c(
  "# Core GBmap Decoded Metadata Audit",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  paste0("- Cells: ", format(nrow(obs), big.mark = ",")),
  paste0("- Authors/studies: ", uniqueN(obs$author)),
  paste0("- Donors: ", uniqueN(obs$donor_id)),
  paste0("- Samples: ", uniqueN(obs$sample)),
  paste0("- iCNV levels: ", paste(sort(unique(na.omit(obs$iCNV))), collapse = ", ")),
  paste0("- LAP3 feature entries: ", nrow(lap3_table)),
  paste0("- GSE131928/Neftel source candidates: ",
         if (nrow(gse131928_candidates)) paste(gse131928_candidates$author, collapse = ", ") else "none"),
  paste0("- Cells after excluding Neftel2019: ",
         format(analysis_entry[population == "exclude_Neftel2019", cells], big.mark = ",")),
  paste0("- Donors after excluding Neftel2019: ",
         analysis_entry[population == "exclude_Neftel2019", donors]),
  paste0("- Neoplastic cells after exclusion: ",
         format(analysis_entry[population == "exclude_Neftel2019_neoplastic", cells], big.mark = ",")),
  "",
  "The source candidate must be reconciled with the paper's study-accession table",
  "before exclusion rules are frozen."
)
writeLines(summary_lines, file.path(output_dir, "README_Decoded_Metadata.md"))

message("Decoded metadata audit complete")
