#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

root <- getwd()
package_date <- "20260702"
out_dir <- file.path("Project_Management", "Submission_Package",
                     paste0("Source_Data_Package_", package_date))

source_index_path <- file.path("Project_Management", "Plans",
                               "supplementary_source_data_index_20260701.csv")
supp_manifest_path <- file.path("Project_Management", "Plans",
                                "supplementary_table_manifest_20260701.csv")

if (dir.exists(out_dir)) {
  unlink(out_dir, recursive = TRUE)
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "source_data"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "supporting_figure_tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "supplementary_tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "manifest"), recursive = TRUE, showWarnings = FALSE)

sanitize <- function(x) {
  x <- gsub("[^A-Za-z0-9_.-]+", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

copy_one <- function(src, dest) {
  if (!file.exists(src)) {
    return(FALSE)
  }
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  file.copy(src, dest, overwrite = TRUE)
}

read_csv <- function(path) {
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

stop_if_missing_cols <- function(x, cols, label) {
  miss <- setdiff(cols, names(x))
  if (length(miss) > 0) {
    stop(label, " missing columns: ", paste(miss, collapse = ", "))
  }
}

source_index <- read_csv(source_index_path)
supp_manifest <- read_csv(supp_manifest_path)
supplementary_table_count <- nrow(supp_manifest)

stop_if_missing_cols(
  source_index,
  c("item_id", "figure", "panel", "manuscript_role", "conclusion",
    "source_data_path", "key_results_path", "panel_map_path", "export_qc_path"),
  "source index"
)

stop_if_missing_cols(
  supp_manifest,
  c("supplementary_item", "proposed_title", "primary_files", "status", "caveat"),
  "supplementary manifest"
)

source_rows <- list()
for (i in seq_len(nrow(source_index))) {
  row <- source_index[i, ]
  fig <- sanitize(row$figure)
  panel <- sanitize(row$panel)
  src <- row$source_data_path
  dest_rel <- file.path(
    "source_data",
    fig,
    paste0(row$item_id, "_", fig, "_Panel_", panel, "_", basename(src))
  )
  dest <- file.path(out_dir, dest_rel)
  copied <- copy_one(src, dest)
  source_rows[[length(source_rows) + 1]] <- data.frame(
    item_id = row$item_id,
    figure = row$figure,
    panel = row$panel,
    role = row$manuscript_role,
    conclusion = row$conclusion,
    original_path = src,
    package_path = dest_rel,
    copied = copied,
    stringsAsFactors = FALSE
  )
}

source_manifest <- do.call(rbind, source_rows)
write.csv(source_manifest,
          file.path(out_dir, "manifest", "source_data_manifest.csv"),
          row.names = FALSE, quote = TRUE)

figure_table_paths <- unique(c(
  source_index$key_results_path,
  source_index$panel_map_path,
  source_index$export_qc_path
))

support_rows <- list()
for (src in figure_table_paths) {
  fig_guess <- "all_figures"
  m <- regmatches(src, regexpr("Figure[0-9][^/]*", src))
  if (length(m) == 1 && nchar(m) > 0) {
    fig_guess <- sanitize(m)
  }
  dest_rel <- file.path("supporting_figure_tables", fig_guess, basename(src))
  copied <- copy_one(src, file.path(out_dir, dest_rel))
  support_rows[[length(support_rows) + 1]] <- data.frame(
    original_path = src,
    package_path = dest_rel,
    copied = copied,
    stringsAsFactors = FALSE
  )
}
support_manifest <- do.call(rbind, support_rows)
write.csv(support_manifest,
          file.path(out_dir, "manifest", "supporting_figure_tables_manifest.csv"),
          row.names = FALSE, quote = TRUE)

expand_primary_file <- function(path_text) {
  path_text <- trimws(path_text)
  if (path_text == "Figure 1-5 tables/figure*_panel_map.csv") {
    return(unique(source_index$panel_map_path))
  }
  if (path_text == "Figure 1-5 tables/figure*_key_results.csv") {
    return(unique(source_index$key_results_path))
  }
  if (grepl("[*?]", path_text)) {
    hits <- Sys.glob(path_text)
    return(hits)
  }
  path_text
}

supp_rows <- list()
for (i in seq_len(nrow(supp_manifest))) {
  item <- supp_manifest[i, ]
  item_dir <- file.path(
    "supplementary_tables",
    paste0(sanitize(gsub("Supplementary Table ", "Table_", item$supplementary_item)),
           "__", sanitize(item$proposed_title))
  )
  parts <- unlist(strsplit(item$primary_files, ";", fixed = TRUE))
  expanded <- unique(unlist(lapply(parts, expand_primary_file)))
  if (length(expanded) == 0) {
    expanded <- character()
  }
  for (src in expanded) {
    dest_rel <- file.path(item_dir, basename(src))
    copied <- copy_one(src, file.path(out_dir, dest_rel))
    supp_rows[[length(supp_rows) + 1]] <- data.frame(
      supplementary_item = item$supplementary_item,
      proposed_title = item$proposed_title,
      original_path = src,
      package_path = dest_rel,
      copied = copied,
      status = item$status,
      caveat = item$caveat,
      stringsAsFactors = FALSE
    )
  }
}

supp_expanded_manifest <- if (length(supp_rows) > 0) {
  do.call(rbind, supp_rows)
} else {
  data.frame()
}
write.csv(supp_expanded_manifest,
          file.path(out_dir, "manifest", "supplementary_tables_expanded_manifest.csv"),
          row.names = FALSE, quote = TRUE)

write.csv(supp_manifest,
          file.path(out_dir, "manifest", "supplementary_table_manifest_original.csv"),
          row.names = FALSE, quote = TRUE)

all_manifest <- rbind(
  data.frame(category = "source_data", source_manifest[, c("original_path", "package_path", "copied")]),
  data.frame(category = "supporting_figure_tables", support_manifest),
  data.frame(category = "supplementary_tables", supp_expanded_manifest[, c("original_path", "package_path", "copied")])
)
write.csv(all_manifest,
          file.path(out_dir, "manifest", "package_file_manifest.csv"),
          row.names = FALSE, quote = TRUE)

missing_files <- all_manifest[!all_manifest$copied, , drop = FALSE]
write.csv(missing_files,
          file.path(out_dir, "manifest", "missing_files.csv"),
          row.names = FALSE, quote = TRUE)

summary_lines <- c(
  "# LAP3 Glioma Source Data Package",
  "",
  paste0("Date: ", package_date),
  "",
  "## Scope",
  "",
  "This package organizes derived source data and supplementary-table inputs for the LAP3 glioma manuscript.",
  "It is intended as a submission-preparation package and not yet a public repository record.",
  "",
  "## Contents",
  "",
  "- `source_data/`: panel-level source data for Figures 1-5 and Supplementary Fig. S4A.",
  "- `supporting_figure_tables/`: figure key-results, panel-map and export-QC tables.",
  paste0("- `supplementary_tables/`: expanded input files for the proposed Supplementary Tables 1-", supplementary_table_count, "."),
  "- `manifest/`: machine-readable file manifests, original paths and copy status.",
  "",
  "## Important Boundaries",
  "",
  "- This package contains derived tables generated by the project. It does not redistribute controlled-access raw data.",
  "- GBM-Space raw Visium data should be accessed through EGA dataset EGAD00001015527 and study EGAS00001005801.",
  "- GLASS raw/current-release data should be accessed through the original Synapse resources.",
  "- TCGA, CGGA, Core GBmap and CPTAC source data should be cited and accessed through their original resources unless redistribution is explicitly allowed.",
  "- `LAP3_STATE_UNION` excludes LAP3 and should be treated as this-study source data.",
  "",
  "## Submission Notes",
  "",
  "- Public archive route is fixed as GitHub repository plus Zenodo archived release DOI.",
  "- Public GitHub repository: https://github.com/lzbkk/glioma-public-data-archive.",
  "- Add a final licence for generated source data and code. Current recommended default is CC BY 4.0 for generated derived data/tables and MIT for code, pending author confirmation.",
  "- CGGA portal-preferred citation was rechecked on 2026-07-02 and fixed to DOI 10.1016/j.gpb.2020.10.005. GBM-Space journal-version status was rechecked on 2026-07-02 and remained preprint-only; recheck again on upload day.",
  "- Neuro-Oncology-family journals, including the current NOA target route, encourage depositing deidentified data and software code in public repositories when possible.",
  "",
  "## Build Summary",
  "",
  paste0("- Source-data rows: ", nrow(source_manifest)),
  paste0("- Supporting figure-table files: ", nrow(support_manifest)),
  paste0("- Supplementary-table expanded files: ", nrow(supp_expanded_manifest)),
  paste0("- Missing files: ", nrow(missing_files)),
  "",
  "## Rebuild",
  "",
  "Run from the repository root:",
  "",
  "```bash",
  "Rscript --vanilla Project_Management/Operations/build_submission_source_data_package.R",
  "```"
)

writeLines(summary_lines, file.path(out_dir, "README.md"))

data_integrity_lines <- c(
  "# Data Integrity Supplement Draft",
  "",
  "Date: 20260702",
  "",
  "Status: draft scaffold for possible R1-stage data/figure responsibility reporting; not submission-frozen.",
  "",
  "## Scope",
  "",
  "This draft supports journal or reviewer requests for a data integrity supplement mapping manuscript figures, supplementary figures and source-data packages to responsible authors.",
  "",
  "## Figure Responsibility Table",
  "",
  "| Display item | Primary data/source package | Responsible author(s) | Verification status | Notes |",
  "|---|---|---|---|---|",
  "| Figure 1 | Source Data Package / Figure_1 | [TBD] | [TBD] | bulk context, anchor benchmark, composition boundary |",
  "| Figure 2 | Source Data Package / Figure_2 | [TBD] | [TBD] | frozen state and submodule decomposition |",
  "| Figure 3 | Source Data Package / Figure_3 | [TBD] | [TBD] | Core GBmap donor-state decomposition |",
  "| Figure 4 | Source Data Package / Figure_4 | [TBD] | [TBD] | GBM-Space spatial topology; tumor-level inference |",
  "| Figure 5 | Source Data Package / Figure_5 | [TBD] | [TBD] | CPTAC/GLASS evidence boundary |",
  "| Supplementary Fig. S4A | Source Data Package / Supplementary_Figure_S4A | [TBD] | [TBD] | representative GBM-Space spatial map; visualization only |",
  "",
  "## Package Responsibility Table",
  "",
  "| Package/file | Responsible author(s) | Verification status | Notes |",
  "|---|---|---|---|",
  "| Source_Data_Package_20260702.tar.gz | [TBD] | [TBD] | derived source data and supplementary inputs only |",
  "| Journal_Supplementary_Tables_20260702.tar.gz | [TBD] | [TBD] | journal-facing supplementary workbooks |",
  "| Supplementary_Tables_1_11_JournalFacing.xlsx | [TBD] | [TBD] | combined supplementary workbook |",
  "",
  "## Boundary",
  "",
  "- This draft does not replace author contribution statements.",
  "- Do not list AI tools as responsible authors.",
  "- Fill only after final author roster and contribution roles are confirmed.",
  "- The package does not redistribute controlled-access raw data or third-party raw matrices."
)
writeLines(data_integrity_lines, file.path(out_dir, "DATA_INTEGRITY_SUPPLEMENT_DRAFT.md"))

audit_lines <- c(
  "# Source Data Package Build Audit",
  "",
  paste0("Date: ", package_date),
  "",
  "## Counts",
  "",
  paste0("- Source-data rows copied: ", sum(source_manifest$copied), "/", nrow(source_manifest)),
  paste0("- Supporting figure tables copied: ", sum(support_manifest$copied), "/", nrow(support_manifest)),
  paste0("- Supplementary expanded files copied: ", sum(supp_expanded_manifest$copied), "/", nrow(supp_expanded_manifest)),
  paste0("- Missing files: ", nrow(missing_files)),
  "",
  "## Blocking Issues",
  "",
  if (nrow(missing_files) == 0) {
    "- None at the file-copy level."
  } else {
    "- See `manifest/missing_files.csv`."
  },
  "",
  "## Not Yet Submission-Frozen",
  "",
  "- Zenodo DOI/record URL is not created.",
  "- Public GitHub repository is fixed as https://github.com/lzbkk/glioma-public-data-archive.",
  "- Final generated-data and code licences require author confirmation.",
  "- GBM-Space journal-version status was rechecked on 2026-07-02 and must be rechecked again on upload day.",
  "- CGGA portal-preferred citation was rechecked on 2026-07-02 and fixed to DOI 10.1016/j.gpb.2020.10.005."
)
writeLines(audit_lines, file.path(out_dir, "BUILD_AUDIT.md"))

cat("Source data package built at:", out_dir, "\n")
cat("Source-data rows:", nrow(source_manifest), "\n")
cat("Supporting figure tables:", nrow(support_manifest), "\n")
cat("Supplementary expanded files:", nrow(supp_expanded_manifest), "\n")
cat("Missing files:", nrow(missing_files), "\n")

if (nrow(missing_files) > 0) {
  quit(status = 2)
}
