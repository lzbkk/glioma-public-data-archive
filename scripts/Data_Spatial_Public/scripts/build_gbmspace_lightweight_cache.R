#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(rhdf5)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop(
    "Usage: build_gbmspace_lightweight_cache.R <h5ad_root_dir> <gene_set_csv> <out_dir> [max_files]",
    call. = FALSE
  )
}

h5ad_root <- normalizePath(args[[1]], mustWork = TRUE)
gene_set_csv <- normalizePath(args[[2]], mustWork = TRUE)
out_dir <- args[[3]]
max_files <- if (length(args) >= 4) as.integer(args[[4]]) else NA_integer_

tables_dir <- file.path(out_dir, "tables")
per_section_dir <- file.path(out_dir, "per_section_non_gene_features")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(per_section_dir, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n")
  flush.console()
}

safe_read <- function(file, path) {
  tryCatch(h5read(file, path), error = function(e) NULL)
}

read_categorical <- function(file, path, n_expected = NULL) {
  categories <- safe_read(file, file.path(path, "categories"))
  codes <- safe_read(file, file.path(path, "codes"))
  if (is.null(categories) || is.null(codes)) {
    value <- safe_read(file, path)
    if (is.null(value)) return(rep(NA_character_, n_expected %||% 0L))
    return(as.character(value))
  }
  categories <- as.character(categories)
  codes <- as.integer(codes)
  out <- categories[codes + 1L]
  out[is.na(codes) | codes < 0L] <- NA_character_
  out
}

`%||%` <- function(x, y) if (is.null(x)) y else x

parse_sample_name <- function(x) {
  tumor_id <- sub("^((AT[0-9]+).*)$", "\\2", x)
  data.table(
    sample_name = x,
    tumor_id = tumor_id,
    site_token = sub("^.*-FO-([^_]+).*$", "\\1", x),
    replicate_token = fifelse(grepl("_", x), sub("^.*_([^_]+)$", "\\1", x), NA_character_)
  )
}

h5ad_files <- list.files(file.path(h5ad_root, "anndata"), pattern = "\\.h5ad$", full.names = TRUE)
h5ad_files <- sort(h5ad_files)
if (!is.na(max_files)) {
  h5ad_files <- head(h5ad_files, max_files)
}
if (length(h5ad_files) == 0) {
  stop("No H5AD files found under: ", file.path(h5ad_root, "anndata"), call. = FALSE)
}

gene_sets <- fread(gene_set_csv)
stopifnot(all(c("signature", "gene") %in% names(gene_sets)))
gene_sets[, gene := toupper(gene)]
gene_sets <- unique(gene_sets[gene != "LAP3", .(signature, gene, fdr_family)])

spot_meta_list <- list()
spot_score_list <- list()
feature_catalog_list <- list()
coverage_list <- list()
section_summary_list <- list()

for (file_i in seq_along(h5ad_files)) {
  file <- h5ad_files[[file_i]]
  sample_id <- sub("\\.h5ad$", "", basename(file))
  log_msg(sprintf("[%d/%d] %s", file_i, length(h5ad_files), sample_id))

  attrs <- h5readAttributes(file, "X")
  shape <- as.integer(attrs$shape)
  n_obs <- shape[[1]]
  n_var <- shape[[2]]

  var_names <- as.character(h5read(file, "var/_index"))
  feature_types <- read_categorical(file, "var/feature_types", n_expected = length(var_names))
  stopifnot(length(var_names) == n_var, length(feature_types) == n_var)
  var_upper <- toupper(var_names)

  obs_index <- as.character(h5read(file, "obs/_index"))
  array_col <- as.integer(h5read(file, "obs/array_col"))
  array_row <- as.integer(h5read(file, "obs/array_row"))
  in_tissue <- as.integer(h5read(file, "obs/in_tissue"))
  sample <- read_categorical(file, "obs/sample", n_expected = n_obs)
  sample_name <- read_categorical(file, "obs/sample_name", n_expected = n_obs)
  spatial <- h5read(file, "obsm/spatial")
  if (nrow(spatial) == 2L && ncol(spatial) == n_obs) {
    spatial_x <- as.numeric(spatial[1, ])
    spatial_y <- as.numeric(spatial[2, ])
  } else if (ncol(spatial) == 2L && nrow(spatial) == n_obs) {
    spatial_x <- as.numeric(spatial[, 1])
    spatial_y <- as.numeric(spatial[, 2])
  } else {
    spatial_x <- rep(NA_real_, n_obs)
    spatial_y <- rep(NA_real_, n_obs)
  }

  data <- as.numeric(h5read(file, "X/data"))
  indices <- as.integer(h5read(file, "X/indices")) + 1L
  indptr <- as.integer(h5read(file, "X/indptr"))

  gene_feature <- feature_types == "Gene Expression"
  non_gene_ids <- which(!gene_feature)
  non_gene_names <- var_names[non_gene_ids]
  non_gene_types <- feature_types[non_gene_ids]
  non_gene_mat <- if (length(non_gene_ids) > 0) {
    matrix(0, nrow = n_obs, ncol = length(non_gene_ids), dimnames = list(NULL, make.names(non_gene_names, unique = TRUE)))
  } else {
    matrix(numeric(), nrow = n_obs, ncol = 0)
  }
  non_gene_pos <- integer(n_var)
  if (length(non_gene_ids) > 0) non_gene_pos[non_gene_ids] <- seq_along(non_gene_ids)

  lap3_idx <- match("LAP3", var_upper)
  signature_genes <- split(gene_sets$gene, gene_sets$signature)
  signature_indices <- lapply(signature_genes, function(g) which(var_upper %in% g & gene_feature))
  signature_sizes <- vapply(signature_indices, length, integer(1))

  lap3_raw <- numeric(n_obs)
  lap3_log1p_cp10k <- numeric(n_obs)
  gene_library_size <- numeric(n_obs)
  detected_gene_features <- integer(n_obs)
  score_matrix <- matrix(
    0,
    nrow = n_obs,
    ncol = length(signature_indices),
    dimnames = list(NULL, names(signature_indices))
  )

  sig_bool <- lapply(signature_indices, function(idx) {
    x <- logical(n_var)
    x[idx] <- TRUE
    x
  })

  for (i in seq_len(n_obs)) {
    start <- indptr[[i]] + 1L
    end <- indptr[[i + 1L]]
    if (end < start) next
    idx <- indices[start:end]
    vals <- data[start:end]

    gene_mask <- gene_feature[idx]
    lib <- sum(vals[gene_mask])
    gene_library_size[[i]] <- lib
    detected_gene_features[[i]] <- sum(gene_mask)

    lap3_hit <- which(idx == lap3_idx)
    if (length(lap3_hit) > 0) {
      lap3_raw[[i]] <- vals[lap3_hit[[1]]]
    }
    if (lib > 0) {
      if (lap3_raw[[i]] > 0) {
        lap3_log1p_cp10k[[i]] <- log1p(lap3_raw[[i]] / lib * 10000)
      }
      log_norm <- log1p(vals / lib * 10000)
      for (sig in names(sig_bool)) {
        hit <- sig_bool[[sig]][idx]
        if (any(hit)) {
          score_matrix[i, sig] <- sum(log_norm[hit]) / signature_sizes[[sig]]
        }
      }
    }

    ng_pos <- non_gene_pos[idx]
    keep_ng <- ng_pos > 0L
    if (any(keep_ng)) {
      non_gene_mat[i, ng_pos[keep_ng]] <- vals[keep_ng]
    }
  }

  parsed <- parse_sample_name(sample_name)
  spot_id <- paste(sample_id, obs_index, sep = "__")
  spot_meta <- data.table(
    spot_id = spot_id,
    barcode = obs_index,
    h5ad_sample_id = sample_id,
    sample = sample,
    sample_name = sample_name,
    tumor_id = parsed$tumor_id,
    site_token = parsed$site_token,
    replicate_token = parsed$replicate_token,
    array_row = array_row,
    array_col = array_col,
    in_tissue = in_tissue,
    spatial_x = spatial_x,
    spatial_y = spatial_y,
    gene_library_size = gene_library_size,
    detected_gene_features = detected_gene_features
  )

  score_dt <- data.table(
    spot_id = spot_id,
    h5ad_sample_id = sample_id,
    sample_name = sample_name,
    tumor_id = parsed$tumor_id,
    LAP3_raw = lap3_raw,
    LAP3_detected = lap3_raw > 0,
    LAP3_log1p_cp10k = lap3_log1p_cp10k,
    gene_library_size = gene_library_size,
    detected_gene_features = detected_gene_features
  )
  score_dt <- cbind(score_dt, as.data.table(score_matrix))

  non_gene_dt <- cbind(data.table(spot_id = spot_id, h5ad_sample_id = sample_id), as.data.table(non_gene_mat))
  non_gene_file <- file.path(per_section_dir, paste0(sample_id, ".non_gene_features.tsv"))
  fwrite(non_gene_dt, non_gene_file, sep = "\t")

  feature_catalog <- data.table(
    h5ad_sample_id = sample_id,
    feature_index_1based = seq_along(var_names),
    feature = var_names,
    feature_upper = var_upper,
    feature_type = feature_types
  )
  feature_catalog_list[[sample_id]] <- feature_catalog
  spot_meta_list[[sample_id]] <- spot_meta
  spot_score_list[[sample_id]] <- score_dt

  coverage_list[[sample_id]] <- rbindlist(lapply(names(signature_indices), function(sig) {
    data.table(
      h5ad_sample_id = sample_id,
      signature = sig,
      n_genes_requested = length(signature_genes[[sig]]),
      n_genes_present = signature_sizes[[sig]],
      present_fraction = signature_sizes[[sig]] / length(signature_genes[[sig]])
    )
  }))

  section_summary_list[[sample_id]] <- data.table(
    h5ad_sample_id = sample_id,
    file = file,
    n_spots = n_obs,
    n_features = n_var,
    n_gene_expression_features = sum(gene_feature),
    n_cell_state_features = sum(feature_types == "Cell state abundances"),
    n_spatial_niche_features = sum(feature_types == "Spatial niche abundances"),
    n_histopath_features = sum(feature_types == "Histopath annotation overlap"),
    has_lap3 = !is.na(lap3_idx),
    n_lap3_detected_spots = sum(lap3_raw > 0),
    lap3_detection_rate = mean(lap3_raw > 0),
    median_gene_library_size = median(gene_library_size),
    non_gene_feature_file = non_gene_file
  )
}

spot_meta <- rbindlist(spot_meta_list, use.names = TRUE, fill = TRUE)
spot_scores <- rbindlist(spot_score_list, use.names = TRUE, fill = TRUE)
feature_catalog <- rbindlist(feature_catalog_list, use.names = TRUE, fill = TRUE)
coverage <- rbindlist(coverage_list, use.names = TRUE, fill = TRUE)
section_summary <- rbindlist(section_summary_list, use.names = TRUE, fill = TRUE)

fwrite(spot_meta, file.path(tables_dir, "gbmspace_spot_metadata.tsv"), sep = "\t")
fwrite(spot_scores, file.path(tables_dir, "gbmspace_spot_lap3_pathway_scores.tsv"), sep = "\t")
fwrite(feature_catalog, file.path(tables_dir, "gbmspace_feature_catalog.tsv"), sep = "\t")
fwrite(coverage, file.path(tables_dir, "gbmspace_pathway_gene_coverage.tsv"), sep = "\t")
fwrite(section_summary, file.path(tables_dir, "gbmspace_section_summary.tsv"), sep = "\t")

saveRDS(
  list(
    generated_at = Sys.time(),
    h5ad_root = h5ad_root,
    gene_set_csv = gene_set_csv,
    section_summary = section_summary,
    pathway_gene_coverage = coverage
  ),
  file.path(out_dir, "gbmspace_lightweight_cache_manifest.rds")
)

cat(
  "# GBM-Space Lightweight Cache\n\n",
  "生成时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n",
  "## 输入\n\n",
  "- H5AD root: `", h5ad_root, "`\n",
  "- frozen pathway genes: `", gene_set_csv, "`\n",
  "- inspected files: `", length(h5ad_files), "`\n\n",
  "## 关键输出\n\n",
  "- `tables/gbmspace_spot_metadata.tsv`\n",
  "- `tables/gbmspace_spot_lap3_pathway_scores.tsv`\n",
  "- `tables/gbmspace_feature_catalog.tsv`\n",
  "- `tables/gbmspace_pathway_gene_coverage.tsv`\n",
  "- `tables/gbmspace_section_summary.tsv`\n",
  "- `per_section_non_gene_features/*.non_gene_features.tsv`\n",
  "- `gbmspace_lightweight_cache_manifest.rds`\n\n",
  "## 第一版结果摘要\n\n",
  "- sections processed: `", nrow(section_summary), "`\n",
  "- total spots: `", nrow(spot_meta), "`\n",
  "- sections with LAP3: `", sum(section_summary$has_lap3), "/", nrow(section_summary), "`\n",
  "- total LAP3-detected spots: `", sum(section_summary$n_lap3_detected_spots), "`\n",
  "- sections with cell-state features: `", sum(section_summary$n_cell_state_features > 0), "/", nrow(section_summary), "`\n",
  "- sections with spatial-niche features: `", sum(section_summary$n_spatial_niche_features > 0), "/", nrow(section_summary), "`\n",
  "- sections with histopath features: `", sum(section_summary$n_histopath_features > 0), "/", nrow(section_summary), "`\n\n",
  "## 方法要点\n\n",
  "H5AD `X` 为 CSR sparse matrix。基因表达、cell-state abundance、spatial niche abundance\n",
  "和 histopath overlap 均通过 `var/feature_types` 区分。Pathway scores 只使用\n",
  "`Gene Expression` features，且从 gene sets 中排除 LAP3，避免自相关。\n\n",
  "## 解释边界\n\n",
  "本缓存只为后续空间统计提供输入，不构成 LAP3 spatial neighborhood 结果。\n",
  "正式分析仍需以 tumor/section 为推断层级，避免把 spot 当独立样本。\n",
  file = file.path(out_dir, "README_Lightweight_Cache.md")
)

log_msg("DONE")
