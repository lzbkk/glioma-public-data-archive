#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

pkg_dir <- file.path(
  "Project_Management", "Submission_Package", "Source_Data_Package_20260702"
)
out_dir <- file.path(
  "Project_Management", "Submission_Package", "Journal_Supplementary_Tables_20260702"
)
manifest_path <- file.path(pkg_dir, "manifest", "supplementary_tables_expanded_manifest.csv")
source_manifest_path <- file.path(pkg_dir, "manifest", "source_data_manifest.csv")

if (dir.exists(out_dir)) {
  unlink(out_dir, recursive = TRUE)
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "individual_tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "manifest"), recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  stop("The openxlsx R package is required.")
}

read_manifest <- function(path) {
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

safe_sheet_name <- function(x, used = character()) {
  x <- gsub("\\.[A-Za-z0-9]+$", "", basename(x))
  x <- gsub("[\\[\\]\\*\\?/\\\\:]", "_", x)
  x <- gsub("[^A-Za-z0-9_]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  if (nchar(x) == 0) {
    x <- "sheet"
  }
  x <- substr(x, 1, 28)
  base <- x
  idx <- 1
  while (tolower(x) %in% tolower(used)) {
    suffix <- paste0("_", idx)
    x <- paste0(substr(base, 1, 31 - nchar(suffix)), suffix)
    idx <- idx + 1
  }
  x
}

read_table_file <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "csv") {
    return(read.csv(path, check.names = FALSE, stringsAsFactors = FALSE))
  }
  if (ext %in% c("tsv", "txt")) {
    return(read.delim(path, check.names = FALSE, stringsAsFactors = FALSE))
  }
  if (ext == "md") {
    lines <- readLines(path, warn = FALSE)
    return(data.frame(line_number = seq_along(lines), text = lines, stringsAsFactors = FALSE))
  }
  data.frame(path = path, note = paste("Unsupported file extension:", ext), stringsAsFactors = FALSE)
}

write_sheet <- function(wb, sheet, dat) {
  openxlsx::addWorksheet(wb, sheet)
  openxlsx::writeData(wb, sheet, dat, withFilter = TRUE)
  openxlsx::freezePane(wb, sheet, firstRow = TRUE)
  widths <- pmin(pmax(nchar(names(dat)) + 2, 12), 45)
  openxlsx::setColWidths(wb, sheet, cols = seq_along(widths), widths = widths)
}

manifest <- read_manifest(manifest_path)
source_manifest <- read_manifest(source_manifest_path)

manifest$package_full_path <- file.path(pkg_dir, manifest$package_path)
manifest$file_exists <- file.exists(manifest$package_full_path)
manifest$file_ext <- tolower(tools::file_ext(manifest$package_full_path))
manifest$file_size_bytes <- ifelse(
  manifest$file_exists,
  file.info(manifest$package_full_path)$size,
  NA_real_
)

table_items <- unique(manifest$supplementary_item)
combined_workbook_name <- paste0(
  "Supplementary_Tables_1_", length(table_items), "_JournalFacing.xlsx"
)
individual_table_range <- paste0(
  "individual_tables/Supplementary_Table_1.xlsx` through `Supplementary_Table_",
  length(table_items), ".xlsx"
)

table_overview <- data.frame(
  supplementary_item = table_items,
  proposed_title = vapply(table_items, function(item) {
    unique(manifest$proposed_title[manifest$supplementary_item == item])[1]
  }, character(1)),
  input_files = vapply(table_items, function(item) {
    sum(manifest$supplementary_item == item)
  }, integer(1)),
  missing_files = vapply(table_items, function(item) {
    sum(!manifest$file_exists[manifest$supplementary_item == item])
  }, integer(1)),
  caveat = vapply(table_items, function(item) {
    unique(manifest$caveat[manifest$supplementary_item == item])[1]
  }, character(1)),
  stringsAsFactors = FALSE
)

global_wb <- openxlsx::createWorkbook()
write_sheet(global_wb, "Overview", table_overview)
write_sheet(global_wb, "SourceDataMap", source_manifest)
write_sheet(global_wb, "SuppInputManifest", manifest)

sheet_manifest_rows <- list()

for (item in table_items) {
  item_rows <- manifest[manifest$supplementary_item == item, , drop = FALSE]
  item_num <- gsub("[^0-9]+", "", item)
  item_title <- unique(item_rows$proposed_title)[1]
  out_name <- paste0("Supplementary_Table_", item_num, ".xlsx")
  item_wb <- openxlsx::createWorkbook()

  item_overview <- data.frame(
    field = c("supplementary_item", "proposed_title", "input_files", "caveat", "status"),
    value = c(item, item_title, nrow(item_rows), unique(item_rows$caveat)[1], unique(item_rows$status)[1]),
    stringsAsFactors = FALSE
  )
  write_sheet(item_wb, "Overview", item_overview)
  write_sheet(item_wb, "InputManifest", item_rows)

  used_item_sheets <- c("Overview", "InputManifest")
  used_global_sheets <- names(global_wb)

  for (i in seq_len(nrow(item_rows))) {
    file_path <- item_rows$package_full_path[i]
    dat <- if (file.exists(file_path)) {
      read_table_file(file_path)
    } else {
      data.frame(path = file_path, error = "missing file", stringsAsFactors = FALSE)
    }
    sheet_base <- safe_sheet_name(file_path, used_item_sheets)
    used_item_sheets <- c(used_item_sheets, sheet_base)
    write_sheet(item_wb, sheet_base, dat)

    global_sheet <- safe_sheet_name(paste0("T", item_num, "_", sheet_base), used_global_sheets)
    used_global_sheets <- c(used_global_sheets, global_sheet)
    write_sheet(global_wb, global_sheet, dat)

    sheet_manifest_rows[[length(sheet_manifest_rows) + 1]] <- data.frame(
      supplementary_item = item,
      proposed_title = item_title,
      original_path = item_rows$original_path[i],
      package_path = item_rows$package_path[i],
      individual_workbook = file.path("individual_tables", out_name),
      individual_sheet = sheet_base,
      combined_workbook = combined_workbook_name,
      combined_sheet = global_sheet,
      rows = nrow(dat),
      columns = ncol(dat),
      file_exists = file.exists(file_path),
      stringsAsFactors = FALSE
    )
  }

  openxlsx::saveWorkbook(
    item_wb,
    file.path(out_dir, "individual_tables", out_name),
    overwrite = TRUE
  )
}

sheet_manifest <- do.call(rbind, sheet_manifest_rows)
write.csv(
  sheet_manifest,
  file.path(out_dir, "manifest", "journal_supplementary_tables_sheet_manifest.csv"),
  row.names = FALSE,
  quote = TRUE
)

write.csv(
  table_overview,
  file.path(out_dir, "manifest", "journal_supplementary_tables_overview.csv"),
  row.names = FALSE,
  quote = TRUE
)

write.csv(
  manifest,
  file.path(out_dir, "manifest", "journal_supplementary_tables_input_manifest.csv"),
  row.names = FALSE,
  quote = TRUE
)

combined_workbook <- file.path(out_dir, combined_workbook_name)
openxlsx::saveWorkbook(global_wb, combined_workbook, overwrite = TRUE)

missing_files <- manifest[!manifest$file_exists, , drop = FALSE]
write.csv(
  missing_files,
  file.path(out_dir, "manifest", "missing_files.csv"),
  row.names = FALSE,
  quote = TRUE
)

qa_lines <- c(
  "# Journal-Facing Supplementary Tables QA",
  "",
  "Date: 2026-07-02",
  "",
  "## Outputs",
  "",
  paste0("- Combined workbook: `", combined_workbook_name, "`"),
  paste0("- Individual workbooks: `", individual_table_range, "`"),
  "- Sheet manifest: `manifest/journal_supplementary_tables_sheet_manifest.csv`",
  "",
  "## Counts",
  "",
  paste0("- Supplementary table items: ", length(table_items)),
  paste0("- Expanded input files: ", nrow(manifest)),
  paste0("- Missing files: ", nrow(missing_files)),
  paste0("- Sheets mapped from source files: ", nrow(sheet_manifest)),
  "",
  "## Boundaries",
  "",
  "- These workbooks organize derived source-data and supplementary-table inputs.",
  "- They do not include GBM-Space raw Visium H5AD/tarball files.",
  "- They do not include GLASS raw Synapse release files.",
  "- Public archive route is fixed as GitHub repository plus Zenodo archived release DOI.",
  "- Public GitHub repository is fixed as https://github.com/lzbkk/glioma-public-data-archive.",
  "- Zenodo DOI/record URL, author-confirmed licence and journal-specific upload naming remain pending.",
  "",
  "## Manual Checks Before Submission",
  "",
  "- Confirm whether the journal wants one combined workbook or separate Supplementary Table files.",
  "- Confirm final Zenodo DOI/record URL and GitHub release URL.",
  "- GBM-Space journal-version status was rechecked on 2026-07-02 and remained preprint-only; recheck again on upload day.",
  "- CGGA portal-preferred citation was rechecked on 2026-07-02 and fixed to DOI 10.1016/j.gpb.2020.10.005.",
  "- Replace provisional table titles if the final manuscript changes figure/table numbering."
)
writeLines(qa_lines, file.path(out_dir, "JOURNAL_SUPPLEMENTARY_TABLES_QA.md"))

readme_lines <- c(
  "# Journal-Facing Supplementary Tables",
  "",
  "Date: 2026-07-02",
  "",
  "This directory converts the Source Data Package inputs into Excel workbooks for manuscript review and journal upload preparation.",
  "",
  "## Files",
  "",
  paste0("- `", combined_workbook_name, "`: one combined workbook containing overview, source-data map, input manifest and all expanded input sheets."),
  "- `individual_tables/`: one workbook per proposed Supplementary Table.",
  "- `manifest/`: sheet-level and input-file-level mapping tables.",
  "- `JOURNAL_SUPPLEMENTARY_TABLES_QA.md`: QA and remaining submission checks.",
  "",
  "## Rebuild",
  "",
  "Run from the repository root:",
  "",
  "```bash",
  "Rscript --vanilla Project_Management/Operations/build_journal_supplementary_tables.R",
  "```"
)
writeLines(readme_lines, file.path(out_dir, "README.md"))

cat("Journal-facing supplementary tables built at:", out_dir, "\n")
cat("Supplementary table items:", length(table_items), "\n")
cat("Expanded input files:", nrow(manifest), "\n")
cat("Missing files:", nrow(missing_files), "\n")
cat("Combined workbook:", combined_workbook, "\n")

if (nrow(missing_files) > 0) {
  quit(status = 2)
}
