#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(rhdf5)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: audit_gbmspace_h5ad_fields.R <extract_dir> <out_dir>", call. = FALSE)
}

extract_dir <- normalizePath(args[[1]], mustWork = TRUE)
out_dir <- args[[2]]
tables_dir <- file.path(out_dir, "tables")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

h5ad_files <- list.files(
  file.path(extract_dir, "anndata"),
  pattern = "\\.h5ad$",
  full.names = TRUE
)
if (length(h5ad_files) == 0) {
  stop("No .h5ad files found under: ", file.path(extract_dir, "anndata"), call. = FALSE)
}

safe_h5ls <- function(file) {
  out <- tryCatch(
    as.data.table(h5ls(file, all = TRUE)),
    error = function(e) data.table(error = conditionMessage(e))
  )
  out[, source_file := basename(file)]
  out
}

safe_read <- function(file, path) {
  tryCatch(h5read(file, path), error = function(e) NULL)
}

read_first_available <- function(file, paths) {
  for (p in paths) {
    value <- safe_read(file, p)
    if (!is.null(value)) {
      return(list(path = p, value = value))
    }
  }
  list(path = NA_character_, value = NULL)
}

as_character_vector <- function(x) {
  if (is.null(x)) return(character())
  if (is.data.frame(x)) return(as.character(unlist(x, use.names = FALSE)))
  if (is.list(x) && !is.atomic(x)) return(as.character(unlist(x, use.names = FALSE)))
  as.character(x)
}

path_exists <- function(ls_dt, path) {
  any(file.path(ls_dt$group, ls_dt$name) == paste0("/", path))
}

get_group_children <- function(ls_dt, group_name) {
  group_name <- paste0("/", sub("^/", "", group_name))
  unique(ls_dt[group == group_name, name])
}

summaries <- list()
structures <- list()
obs_cols <- list()
var_cols <- list()
obsm_keys <- list()
uns_keys <- list()
lap3_checks <- list()
feature_type_counts <- list()

for (file in h5ad_files) {
  sample_id <- sub("\\.h5ad$", "", basename(file))
  ls_dt <- safe_h5ls(file)
  structures[[sample_id]] <- copy(ls_dt)

  if ("error" %in% names(ls_dt)) {
    summaries[[sample_id]] <- data.table(
      sample_id = sample_id,
      file = file,
      file_bytes = file.info(file)$size,
      readable = FALSE,
      n_h5_entries = NA_integer_,
      x_shape = NA_character_,
      obs_n = NA_character_,
      var_n = NA_character_,
      error = ls_dt$error[1]
    )
    next
  }

  obs_children <- get_group_children(ls_dt, "obs")
  var_children <- get_group_children(ls_dt, "var")
  obsm_children <- get_group_children(ls_dt, "obsm")
  uns_children <- get_group_children(ls_dt, "uns")

  obs_cols[[sample_id]] <- data.table(sample_id = sample_id, obs_column = obs_children)
  var_cols[[sample_id]] <- data.table(sample_id = sample_id, var_column = var_children)
  obsm_keys[[sample_id]] <- data.table(sample_id = sample_id, obsm_key = obsm_children)
  uns_keys[[sample_id]] <- data.table(sample_id = sample_id, uns_key = uns_children)

  x_shape <- read_first_available(file, c("X/shape", "raw/X/shape"))
  x_shape_text <- if (length(x_shape$value) > 0) paste(as_character_vector(x_shape$value), collapse = " x ") else NA_character_

  var_index <- read_first_available(
    file,
    c("var/_index", "var/index", "var/gene_symbols", "var/gene_symbol", "var/gene_ids", "var/features", "var/feature_name")
  )
  var_index_vec <- as_character_vector(var_index$value)
  var_index_vec <- var_index_vec[nzchar(var_index_vec)]

  feature_types <- character(length(var_index_vec))
  feature_type_path <- NA_character_
  ft_categories <- safe_read(file, "var/feature_types/categories")
  ft_codes <- safe_read(file, "var/feature_types/codes")
  if (!is.null(ft_categories) && !is.null(ft_codes)) {
    feature_type_path <- "var/feature_types"
    ft_categories <- as_character_vector(ft_categories)
    ft_codes <- as.integer(ft_codes)
    feature_types <- ft_categories[ft_codes + 1L]
  } else {
    ft_direct <- read_first_available(file, c("var/feature_type", "var/feature_types"))
    if (!is.null(ft_direct$value)) {
      feature_type_path <- ft_direct$path
      feature_types <- as_character_vector(ft_direct$value)
    }
  }
  if (length(feature_types) != length(var_index_vec)) {
    feature_types <- rep(NA_character_, length(var_index_vec))
  }

  ft_tab <- as.data.table(table(feature_types, useNA = "ifany"))
  setnames(ft_tab, c("feature_type", "n"))
  ft_tab[, sample_id := sample_id]
  feature_type_counts[[sample_id]] <- ft_tab[, .(sample_id, feature_type, n)]

  lap3_exact <- any(toupper(var_index_vec) == "LAP3")
  lap3_feature_type <- feature_types[match("LAP3", toupper(var_index_vec))]
  if (length(lap3_feature_type) == 0 || is.na(lap3_feature_type)) lap3_feature_type <- NA_character_
  lap3_contains <- unique(var_index_vec[grepl("(^|[^A-Z0-9])LAP3([^A-Z0-9]|$)", toupper(var_index_vec))])
  if (length(lap3_contains) == 0) lap3_contains <- NA_character_

  obs_index <- read_first_available(file, c("obs/_index", "obs/index", "obs/barcode", "obs/barcodes"))
  obs_index_vec <- as_character_vector(obs_index$value)

  lap3_checks[[sample_id]] <- data.table(
    sample_id = sample_id,
    gene_source_path = var_index$path,
    feature_type_path = feature_type_path,
    n_gene_names_read = length(var_index_vec),
    lap3_exact = lap3_exact,
    lap3_feature_type = lap3_feature_type,
    lap3_matching_values = paste(head(lap3_contains, 20), collapse = ";")
  )

  summaries[[sample_id]] <- data.table(
    sample_id = sample_id,
    file = file,
    file_bytes = file.info(file)$size,
    readable = TRUE,
    n_h5_entries = nrow(ls_dt),
    x_shape_path = x_shape$path,
    x_shape = x_shape_text,
    obs_index_path = obs_index$path,
    n_obs_index_read = length(obs_index_vec),
    var_index_path = var_index$path,
    n_var_index_read = length(var_index_vec),
    feature_type_path = feature_type_path,
    n_obs_children = length(obs_children),
    n_var_children = length(var_children),
    n_obsm_children = length(obsm_children),
    n_uns_children = length(uns_children),
    has_spatial_obsm = any(grepl("spatial|coord|position|array", obsm_children, ignore.case = TRUE)),
    has_cell_state_feature_type = any(grepl("cell state", feature_types, ignore.case = TRUE)),
    has_spatial_niche_feature_type = any(grepl("spatial niche", feature_types, ignore.case = TRUE)),
    has_histology_feature_type = any(grepl("histopath", feature_types, ignore.case = TRUE)),
    has_patient_section_field = any(grepl("patient|donor|tumou?r|section|slide|sample", obs_children, ignore.case = TRUE)),
    error = NA_character_
  )
}

structure_dt <- rbindlist(structures, fill = TRUE)
summary_dt <- rbindlist(summaries, fill = TRUE)
obs_dt <- rbindlist(obs_cols, fill = TRUE)
var_dt <- rbindlist(var_cols, fill = TRUE)
obsm_dt <- rbindlist(obsm_keys, fill = TRUE)
uns_dt <- rbindlist(uns_keys, fill = TRUE)
lap3_dt <- rbindlist(lap3_checks, fill = TRUE)
feature_type_dt <- rbindlist(feature_type_counts, fill = TRUE)

fwrite(summary_dt, file.path(tables_dir, "gbmspace_h5ad_summary.tsv"), sep = "\t")
fwrite(lap3_dt, file.path(tables_dir, "gbmspace_h5ad_lap3_gene_check.tsv"), sep = "\t")
fwrite(feature_type_dt, file.path(tables_dir, "gbmspace_h5ad_feature_type_counts.tsv"), sep = "\t")
fwrite(obs_dt, file.path(tables_dir, "gbmspace_h5ad_obs_columns.tsv"), sep = "\t")
fwrite(var_dt, file.path(tables_dir, "gbmspace_h5ad_var_columns.tsv"), sep = "\t")
fwrite(obsm_dt, file.path(tables_dir, "gbmspace_h5ad_obsm_keys.tsv"), sep = "\t")
fwrite(uns_dt, file.path(tables_dir, "gbmspace_h5ad_uns_keys.tsv"), sep = "\t")
fwrite(structure_dt, file.path(tables_dir, "gbmspace_h5ad_structure_long.tsv"), sep = "\t")

readme_path <- file.path(extract_dir, "anndata", "README.md")
readme_note <- if (file.exists(readme_path)) {
  paste(readLines(readme_path, warn = FALSE), collapse = "\n")
} else {
  "README.md was not extracted."
}

cat(
  "# GBM-Space H5AD Field Audit\n\n",
  "生成时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n",
  "## 输入\n\n",
  "- selective extraction directory: `", extract_dir, "`\n",
  "- inspected h5ad files: `", length(h5ad_files), "`\n\n",
  "## 关键输出\n\n",
  "- `tables/gbmspace_h5ad_summary.tsv`\n",
  "- `tables/gbmspace_h5ad_lap3_gene_check.tsv`\n",
  "- `tables/gbmspace_h5ad_feature_type_counts.tsv`\n",
  "- `tables/gbmspace_h5ad_obs_columns.tsv`\n",
  "- `tables/gbmspace_h5ad_var_columns.tsv`\n",
  "- `tables/gbmspace_h5ad_obsm_keys.tsv`\n",
  "- `tables/gbmspace_h5ad_uns_keys.tsv`\n",
  "- `tables/gbmspace_h5ad_structure_long.tsv`\n\n",
  "## 第一版结果摘要\n\n",
  "- readable files: `", sum(summary_dt$readable, na.rm = TRUE), "/", nrow(summary_dt), "`\n",
  "- LAP3 exact gene hits: `", sum(lap3_dt$lap3_exact, na.rm = TRUE), "/", nrow(lap3_dt), "`\n",
  "- files with spatial-like `obsm` keys: `", sum(summary_dt$has_spatial_obsm, na.rm = TRUE), "/", nrow(summary_dt), "`\n",
  "- files with patient/section-like `obs` fields: `", sum(summary_dt$has_patient_section_field, na.rm = TRUE), "/", nrow(summary_dt), "`\n",
  "- files with cell-state feature types: `", sum(summary_dt$has_cell_state_feature_type, na.rm = TRUE), "/", nrow(summary_dt), "`\n",
  "- files with spatial-niche feature types: `", sum(summary_dt$has_spatial_niche_feature_type, na.rm = TRUE), "/", nrow(summary_dt), "`\n",
  "- files with histopath feature types: `", sum(summary_dt$has_histology_feature_type, na.rm = TRUE), "/", nrow(summary_dt), "`\n\n",
  "## 解释边界\n\n",
  "本步骤只审计少量代表性 H5AD 的字段结构和 LAP3 gene coverage，不进行空间邻域统计。\n",
  "若字段结构一致，下一步应批量抽取所有 97 个 H5AD 的 lightweight cache；若字段不一致，需先按文件类型分层。\n\n",
  "## Extracted README\n\n",
  "```text\n", readme_note, "\n```\n",
  file = file.path(out_dir, "README_H5AD_Field_Audit.md")
)
