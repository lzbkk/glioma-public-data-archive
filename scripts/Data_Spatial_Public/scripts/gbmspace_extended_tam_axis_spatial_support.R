#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(rhdf5)
  library(FNN)
})

project_root <- "/home/lzb/glioma"
args <- commandArgs(trailingOnly = TRUE)
max_files <- if (length(args) >= 1L && nzchar(args[[1]])) as.integer(args[[1]]) else NA_integer_

cache_dir <- file.path(project_root, "Data_Spatial_Public/GBM_Space/results/Lightweight_Cache")
manifest_file <- file.path(cache_dir, "gbmspace_lightweight_cache_manifest.rds")
state_score_file <- file.path(
  project_root,
  "Data_Spatial_Public/GBM_Space/results/LAP3_State_Spatial_Topology/source_data/gbmspace_spot_lap3_state_scores.tsv.gz"
)
metadata_file <- file.path(cache_dir, "tables/gbmspace_spot_metadata.tsv")
nomination_file <- file.path(
  project_root,
  "Data_scRNA_GEO/GBmap_Extended/results/Focused_Malignant_TAM_Communication/tables/extended_focused_tam_nomination_table.csv"
)

out_dir <- file.path(project_root, "Data_Spatial_Public/GBM_Space/results/Extended_TAM_Axis_Spatial_Support")
tables_dir <- file.path(out_dir, "tables")
source_dir <- file.path(out_dir, "source_data")
logs_dir <- file.path(out_dir, "logs")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

data.table::setDTthreads(8)

message_ts <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n")
  flush.console()
}

`%||%` <- function(x, y) if (is.null(x)) y else x

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

parse_tumor_id <- function(x) {
  out <- sub("^((AT[0-9]+).*)$", "\\2", x)
  out[!grepl("^AT[0-9]+$", out)] <- NA_character_
  out
}

safe_cor <- function(x, y, method = "spearman") {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 20 || uniqueN(x[keep]) < 3 || uniqueN(y[keep]) < 3) return(NA_real_)
  suppressWarnings(cor(x[keep], y[keep], method = method))
}

residual_rank_cor <- function(dt, xcol, ycol) {
  keep <- is.finite(dt[[xcol]]) &
    is.finite(dt[[ycol]]) &
    is.finite(dt$gene_library_size) &
    is.finite(dt$detected_gene_features)

  if (sum(keep) < 20 ||
      uniqueN(dt[[xcol]][keep]) < 3 ||
      uniqueN(dt[[ycol]][keep]) < 3) {
    return(NA_real_)
  }

  sub <- dt[keep]
  rx <- rank(sub[[xcol]], ties.method = "average")
  ry <- rank(sub[[ycol]], ties.method = "average")
  rlib <- rank(log1p(sub$gene_library_size), ties.method = "average")
  rfeat <- rank(sub$detected_gene_features, ties.method = "average")

  xres <- residuals(lm(rx ~ rlib + rfeat))
  yres <- residuals(lm(ry ~ rlib + rfeat))
  suppressWarnings(cor(xres, yres, method = "pearson"))
}

wilcox_vs_zero <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 4 || all(x == 0)) return(NA_real_)
  suppressWarnings(wilcox.test(x, mu = 0, exact = FALSE)$p.value)
}

z_within <- function(x) {
  if (sum(is.finite(x)) < 3) return(rep(NA_real_, length(x)))
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  as.numeric((x - mean(x, na.rm = TRUE)) / s)
}

make_neighbor_means <- function(sec, readout_cols, k = 6L) {
  coords <- as.matrix(sec[, .(spatial_x, spatial_y)])
  ok <- is.finite(coords[, 1]) & is.finite(coords[, 2])
  out <- matrix(NA_real_, nrow = nrow(sec), ncol = length(readout_cols))
  colnames(out) <- paste0("neighbor_", readout_cols)
  if (sum(ok) <= k + 1L || length(readout_cols) == 0L) return(as.data.table(out))

  nn <- get.knn(coords[ok, , drop = FALSE], k = k)$nn.index
  mat <- as.matrix(sec[ok, ..readout_cols])
  storage.mode(mat) <- "double"
  neigh <- matrix(NA_real_, nrow = nrow(mat), ncol = ncol(mat))
  for (j in seq_len(ncol(mat))) {
    vals <- mat[, j]
    neigh[, j] <- rowMeans(matrix(vals[nn], nrow = nrow(nn)), na.rm = TRUE)
  }
  out[which(ok), ] <- neigh
  as.data.table(out)
}

extract_axis_genes_for_h5ad <- function(file, target_genes) {
  sample_id <- sub("\\.h5ad$", "", basename(file))
  attrs <- h5readAttributes(file, "X")
  shape <- as.integer(attrs$shape)
  n_obs <- shape[[1]]
  n_var <- shape[[2]]

  var_names <- as.character(h5read(file, "var/_index"))
  feature_types <- read_categorical(file, "var/feature_types", n_expected = length(var_names))
  stopifnot(length(var_names) == n_var, length(feature_types) == n_var)
  var_upper <- toupper(var_names)
  gene_feature <- feature_types == "Gene Expression"

  obs_index <- as.character(h5read(file, "obs/_index"))
  sample_name <- read_categorical(file, "obs/sample_name", n_expected = n_obs)
  if (all(is.na(sample_name))) sample_name <- rep(sample_id, n_obs)

  data <- as.numeric(h5read(file, "X/data"))
  indices <- as.integer(h5read(file, "X/indices")) + 1L
  indptr <- as.integer(h5read(file, "X/indptr"))

  target_idx <- match(toupper(target_genes), var_upper)
  target_present <- !is.na(target_idx) & gene_feature[target_idx]
  target_idx_present <- target_idx[target_present]
  gene_names_present <- toupper(target_genes)[target_present]

  raw_mat <- matrix(0, nrow = n_obs, ncol = length(target_genes))
  log_mat <- matrix(0, nrow = n_obs, ncol = length(target_genes))
  colnames(raw_mat) <- paste0("raw_", toupper(target_genes))
  colnames(log_mat) <- paste0("log1p_cp10k_", toupper(target_genes))

  for (i in seq_len(n_obs)) {
    start <- indptr[[i]] + 1L
    end <- indptr[[i + 1L]]
    if (end < start) next
    idx <- indices[start:end]
    vals <- data[start:end]
    gene_mask <- gene_feature[idx]
    lib <- sum(vals[gene_mask])
    if (lib <= 0 || length(target_idx_present) == 0L) next

    hit_pos <- match(target_idx_present, idx)
    found <- !is.na(hit_pos)
    if (!any(found)) next
    target_col <- match(gene_names_present[found], toupper(target_genes))
    target_vals <- vals[hit_pos[found]]
    raw_mat[i, target_col] <- target_vals
    log_mat[i, target_col] <- log1p(target_vals / lib * 10000)
  }

  spot_id <- paste(sample_id, obs_index, sep = "__")
  scores <- data.table(
    spot_id = spot_id,
    h5ad_sample_id = sample_id,
    sample_name = sample_name,
    tumor_id = parse_tumor_id(sample_name)
  )
  scores <- cbind(scores, as.data.table(raw_mat), as.data.table(log_mat))

  coverage <- data.table(
    h5ad_sample_id = sample_id,
    gene = toupper(target_genes),
    present = target_present,
    var_index = target_idx
  )

  list(scores = scores, coverage = coverage)
}

axis_definitions <- data.table(
  axis = c("MIF_CD74_CXCR4_CD44", "SPP1_CD44_ITGAV_ITGB1", "LGALS3_CD44_ITGB1"),
  ligand = c("MIF", "SPP1", "LGALS3"),
  receptors = c("CD74;CXCR4;CD44", "CD44;ITGAV;ITGB1", "CD44;ITGB1"),
  nominated_direction = c("malignant_to_tam_or_tam_to_malignant", "tam_to_malignant", "tam_to_malignant")
)
target_genes <- sort(unique(c(axis_definitions$ligand, unlist(strsplit(axis_definitions$receptors, ";")))))

stopifnot(file.exists(manifest_file), file.exists(state_score_file), file.exists(metadata_file))
manifest <- readRDS(manifest_file)
section_summary <- as.data.table(manifest$section_summary)
if (!is.na(max_files)) section_summary <- head(section_summary, max_files)
if (nrow(section_summary) == 0L) stop("No sections available in manifest.", call. = FALSE)

message_ts("Loading spot metadata and LAP3-state scores")
metadata <- fread(metadata_file, select = c(
  "spot_id", "h5ad_sample_id", "tumor_id", "spatial_x", "spatial_y",
  "gene_library_size", "detected_gene_features"
))
state_scores <- fread(state_score_file, select = c(
  "spot_id", "h5ad_sample_id", "tumor_id",
  "LAP3_STATE_UNION", "LAP3_STATE_UNION_NO_TRANSLATION_PROTEOSTASIS",
  "gene_library_size", "detected_gene_features"
))

message_ts("Extracting candidate-axis genes from H5AD files")
gene_score_rows <- list()
coverage_rows <- list()
for (i in seq_len(nrow(section_summary))) {
  file <- section_summary$file[[i]]
  sample_id <- section_summary$h5ad_sample_id[[i]]
  message_ts(sprintf("[%d/%d] %s", i, nrow(section_summary), sample_id))
  obj <- extract_axis_genes_for_h5ad(file, target_genes)
  gene_score_rows[[sample_id]] <- obj$scores
  coverage_rows[[sample_id]] <- obj$coverage
}

gene_scores <- rbindlist(gene_score_rows, use.names = TRUE, fill = TRUE)
gene_coverage <- rbindlist(coverage_rows, use.names = TRUE, fill = TRUE)

fwrite(gene_scores, file.path(source_dir, "gbmspace_extended_tam_axis_gene_scores.tsv.gz"), sep = "\t")
fwrite(gene_coverage, file.path(tables_dir, "gbmspace_extended_tam_axis_gene_coverage_by_section.tsv"), sep = "\t")

coverage_summary <- gene_coverage[, .(
  n_sections = .N,
  n_present_sections = sum(present),
  present_fraction = mean(present)
), by = gene][order(gene)]
fwrite(coverage_summary, file.path(tables_dir, "gbmspace_extended_tam_axis_gene_coverage_summary.tsv"), sep = "\t")

message_ts("Building spot-level axis scores")
spot_dt <- merge(metadata, state_scores[
  ,
  .(
    spot_id, h5ad_sample_id,
    LAP3_STATE_UNION, LAP3_STATE_UNION_NO_TRANSLATION_PROTEOSTASIS
  )
], by = c("spot_id", "h5ad_sample_id"), all.x = TRUE)
spot_dt <- merge(spot_dt, gene_scores, by = c("spot_id", "h5ad_sample_id"), all.x = TRUE, suffixes = c("", "_gene"))
spot_dt[, tumor_id := fifelse(!is.na(tumor_id), tumor_id, tumor_id_gene)]
spot_dt[, tumor_id_gene := NULL]
spot_dt[, sample_name := NULL]

axis_score_cols <- character()
for (j in seq_len(nrow(axis_definitions))) {
  axis <- axis_definitions$axis[[j]]
  ligand <- toupper(axis_definitions$ligand[[j]])
  receptors <- toupper(unlist(strsplit(axis_definitions$receptors[[j]], ";")))
  ligand_col <- paste0("log1p_cp10k_", ligand)
  receptor_cols <- paste0("log1p_cp10k_", receptors)
  receptor_score_col <- paste0(axis, "_receptor_mean")
  ligand_z_col <- paste0(axis, "_ligand_z")
  receptor_z_col <- paste0(axis, "_receptor_z")
  axis_score_col <- paste0(axis, "_axis_score")

  spot_dt[, (receptor_score_col) := rowMeans(.SD, na.rm = TRUE), .SDcols = receptor_cols]
  spot_dt[!is.finite(get(receptor_score_col)), (receptor_score_col) := NA_real_]
  spot_dt[, (ligand_z_col) := z_within(get(ligand_col)), by = h5ad_sample_id]
  spot_dt[, (receptor_z_col) := z_within(get(receptor_score_col)), by = h5ad_sample_id]
  spot_dt[, (axis_score_col) := rowMeans(.SD, na.rm = TRUE), .SDcols = c(ligand_z_col, receptor_z_col)]
  spot_dt[!is.finite(get(axis_score_col)), (axis_score_col) := NA_real_]
  axis_score_cols <- c(axis_score_cols, axis_score_col)
}

fwrite(
  spot_dt[, c(
    "spot_id", "h5ad_sample_id", "tumor_id",
    "LAP3_STATE_UNION", "LAP3_STATE_UNION_NO_TRANSLATION_PROTEOSTASIS",
    axis_score_cols
  ), with = FALSE],
  file.path(source_dir, "gbmspace_extended_tam_axis_spot_scores.tsv.gz"),
  sep = "\t"
)

priority_tam_readouts <- c(
  "Proliferative.TAMs", "Resident.TAMs", "Resident.BAM.TAMs", "RTN1..TAMs",
  "Immune..TAMs.", "Monocytes", "Pro.inflammatory.TAMs", "Anti.inflammatory.TAMs",
  "Interferon.TAMs", "Angiogenic.TAMs", "Dendritic.cells"
)

message_ts("Testing section-level spatial support")
section_effect_rows <- list()
for (i in seq_len(nrow(section_summary))) {
  sample_id <- section_summary$h5ad_sample_id[[i]]
  ng_file <- section_summary$non_gene_feature_file[[i]]
  if (is.na(ng_file) || !file.exists(ng_file)) next
  sec <- spot_dt[h5ad_sample_id == sample_id]
  ng <- fread(ng_file)
  sec <- merge(sec, ng, by = c("spot_id", "h5ad_sample_id"), all.x = TRUE)
  tam_cols <- intersect(priority_tam_readouts, names(sec))
  if (nrow(sec) < 20 || length(tam_cols) == 0L) next

  neighbor_dt <- make_neighbor_means(sec, tam_cols, k = 6L)
  sec <- cbind(sec, neighbor_dt)

  for (j in seq_len(nrow(axis_definitions))) {
    axis <- axis_definitions$axis[[j]]
    axis_score_col <- paste0(axis, "_axis_score")
    target_cols <- c("LAP3_STATE_UNION", "LAP3_STATE_UNION_NO_TRANSLATION_PROTEOSTASIS")
    for (target in target_cols) {
      section_effect_rows[[length(section_effect_rows) + 1L]] <- data.table(
        h5ad_sample_id = sample_id,
        tumor_id = sec$tumor_id[[1]],
        axis = axis,
        readout = target,
        readout_class = "lap3_state",
        analysis_type = "same_spot_colocalization",
        n_spots = sum(is.finite(sec[[axis_score_col]]) & is.finite(sec[[target]])),
        raw_rho = safe_cor(sec[[axis_score_col]], sec[[target]]),
        depth_adjusted_rho = residual_rank_cor(sec, axis_score_col, target)
      )
    }
    for (target in tam_cols) {
      section_effect_rows[[length(section_effect_rows) + 1L]] <- data.table(
        h5ad_sample_id = sample_id,
        tumor_id = sec$tumor_id[[1]],
        axis = axis,
        readout = target,
        readout_class = "tam_same_spot",
        analysis_type = "same_spot_colocalization",
        n_spots = sum(is.finite(sec[[axis_score_col]]) & is.finite(sec[[target]])),
        raw_rho = safe_cor(sec[[axis_score_col]], sec[[target]]),
        depth_adjusted_rho = residual_rank_cor(sec, axis_score_col, target)
      )
      neighbor_col <- paste0("neighbor_", target)
      section_effect_rows[[length(section_effect_rows) + 1L]] <- data.table(
        h5ad_sample_id = sample_id,
        tumor_id = sec$tumor_id[[1]],
        axis = axis,
        readout = target,
        readout_class = "tam_neighbor_k6",
        analysis_type = "neighbor_topology_k6",
        n_spots = sum(is.finite(sec[[axis_score_col]]) & is.finite(sec[[neighbor_col]])),
        raw_rho = safe_cor(sec[[axis_score_col]], sec[[neighbor_col]]),
        depth_adjusted_rho = residual_rank_cor(sec, axis_score_col, neighbor_col)
      )
    }
  }
}

section_effects <- rbindlist(section_effect_rows, use.names = TRUE, fill = TRUE)
fwrite(section_effects, file.path(tables_dir, "gbmspace_extended_tam_axis_section_effects.tsv"), sep = "\t")

tumor_effects <- section_effects[
  is.finite(depth_adjusted_rho),
  .(
    n_sections = .N,
    median_raw_rho = median(raw_rho, na.rm = TRUE),
    median_depth_adjusted_rho = median(depth_adjusted_rho, na.rm = TRUE)
  ),
  by = .(tumor_id, axis, readout, readout_class, analysis_type)
]
fwrite(tumor_effects, file.path(tables_dir, "gbmspace_extended_tam_axis_tumor_effects.tsv"), sep = "\t")

summary_dt <- tumor_effects[
  ,
  .(
    n_tumors = .N,
    median_raw_rho = median(median_raw_rho, na.rm = TRUE),
    median_depth_adjusted_rho = median(median_depth_adjusted_rho, na.rm = TRUE),
    n_positive_depth_adjusted = sum(median_depth_adjusted_rho > 0, na.rm = TRUE),
    p_depth_adjusted = wilcox_vs_zero(median_depth_adjusted_rho)
  ),
  by = .(axis, readout, readout_class, analysis_type)
]
summary_dt[, fdr_depth_adjusted := p.adjust(p_depth_adjusted, method = "BH"), by = .(readout_class, analysis_type)]
summary_dt[, abs_median_depth_adjusted_rho := abs(median_depth_adjusted_rho)]
setorder(summary_dt, axis, readout_class, analysis_type, fdr_depth_adjusted, -abs_median_depth_adjusted_rho)
fwrite(summary_dt, file.path(tables_dir, "gbmspace_extended_tam_axis_spatial_support_summary.tsv"), sep = "\t")

loto_rows <- list()
for (ax in unique(tumor_effects$axis)) {
  for (ro in unique(tumor_effects$readout)) {
    for (cls in unique(tumor_effects[axis == ax & readout == ro]$readout_class)) {
      te <- tumor_effects[axis == ax & readout == ro & readout_class == cls]
      if (nrow(te) < 4) next
      for (drop_tumor in unique(te$tumor_id)) {
        sub <- te[tumor_id != drop_tumor]
        loto_rows[[length(loto_rows) + 1L]] <- data.table(
          axis = ax,
          readout = ro,
          readout_class = cls,
          dropped_tumor_id = drop_tumor,
          n_tumors = nrow(sub),
          median_depth_adjusted_rho = median(sub$median_depth_adjusted_rho, na.rm = TRUE),
          n_positive_depth_adjusted = sum(sub$median_depth_adjusted_rho > 0, na.rm = TRUE),
          p_depth_adjusted = wilcox_vs_zero(sub$median_depth_adjusted_rho)
        )
      }
    }
  }
}
loto <- rbindlist(loto_rows, use.names = TRUE, fill = TRUE)
fwrite(loto, file.path(tables_dir, "gbmspace_extended_tam_axis_leave_one_tumor_out.tsv"), sep = "\t")

state_support <- summary_dt[
  readout == "LAP3_STATE_UNION_NO_TRANSLATION_PROTEOSTASIS",
  .(
    state_median_rho = median_depth_adjusted_rho,
    state_positive_tumors = n_positive_depth_adjusted,
    state_fdr = fdr_depth_adjusted
  ),
  by = axis
]
tam_neighbor_best <- summary_dt[
  readout_class == "tam_neighbor_k6",
  .SD[which.max(median_depth_adjusted_rho)],
  by = axis
][
  ,
  .(
    axis,
    best_neighbor_tam_readout = readout,
    neighbor_median_rho = median_depth_adjusted_rho,
    neighbor_positive_tumors = n_positive_depth_adjusted,
    neighbor_fdr = fdr_depth_adjusted
  )
]
same_spot_best <- summary_dt[
  readout_class == "tam_same_spot",
  .SD[which.max(median_depth_adjusted_rho)],
  by = axis
][
  ,
  .(
    axis,
    best_same_spot_tam_readout = readout,
    same_spot_median_rho = median_depth_adjusted_rho,
    same_spot_positive_tumors = n_positive_depth_adjusted,
    same_spot_fdr = fdr_depth_adjusted
  )
]

verdict <- merge(axis_definitions, state_support, by = "axis", all.x = TRUE)
verdict <- merge(verdict, tam_neighbor_best, by = "axis", all.x = TRUE)
verdict <- merge(verdict, same_spot_best, by = "axis", all.x = TRUE)
verdict[, state_supportive := state_median_rho > 0 & state_positive_tumors >= 8]
verdict[, tam_neighbor_supportive := neighbor_median_rho > 0 & neighbor_positive_tumors >= 8]
verdict[, same_spot_tam_supportive := same_spot_median_rho > 0 & same_spot_positive_tumors >= 8]
verdict[, spatial_support_verdict := fifelse(
  state_supportive & (tam_neighbor_supportive | same_spot_tam_supportive),
  "supportive_spatial_colocalization",
  fifelse(state_supportive, "state_only_support", "not_supported")
)]
verdict[, interpretation_boundary := "Visium spot-level co-expression/topology support only; not directional cell-cell communication or causal mechanism."]
fwrite(verdict, file.path(tables_dir, "gbmspace_extended_tam_axis_closure_verdict.tsv"), sep = "\t")

nomination <- if (file.exists(nomination_file)) fread(nomination_file) else data.table()
if (nrow(nomination) > 0L) {
  merged_nomination <- merge(nomination, verdict, by = "axis", all.x = TRUE)
  fwrite(merged_nomination, file.path(tables_dir, "extended_nomination_with_gbmspace_spatial_support.tsv"), sep = "\t")
}

top_summary <- summary_dt[
  readout_class %in% c("lap3_state", "tam_same_spot", "tam_neighbor_k6")
][order(axis, readout_class, fdr_depth_adjusted, -abs_median_depth_adjusted_rho)]

cat(
  "# GBM-Space Extended TAM Axis Spatial Support\n\n",
  "Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n",
  "## Scope\n\n",
  "This mini-check projects the three Extended GBmap nominated malignant/TAM axes into GBM-Space Visium sections.\n",
  "It tests spot-level co-expression and k=6 spatial-neighbor topology only. It does not infer directional ligand-receptor communication.\n\n",
  "## Inputs\n\n",
  "- H5AD sections: `", nrow(section_summary), "`\n",
  "- Tumors: `", uniqueN(metadata$tumor_id), "`\n",
  "- Candidate axes: `", paste(axis_definitions$axis, collapse = "`, `"), "`\n",
  "- Target genes: `", paste(target_genes, collapse = "`, `"), "`\n\n",
  "## Gene coverage\n\n",
  "```text\n",
  paste(capture.output(print(coverage_summary)), collapse = "\n"),
  "\n```\n\n",
  "## Closure verdict\n\n",
  "```text\n",
  paste(capture.output(print(verdict[, .(
    axis, ligand, receptors, state_median_rho, state_positive_tumors,
    best_neighbor_tam_readout, neighbor_median_rho, neighbor_positive_tumors,
    best_same_spot_tam_readout, same_spot_median_rho, same_spot_positive_tumors,
    spatial_support_verdict
  )])), collapse = "\n"),
  "\n```\n\n",
  "## Top statistics\n\n",
  "```text\n",
  paste(capture.output(print(head(top_summary[, .(
    axis, readout_class, analysis_type, readout, n_tumors,
    median_depth_adjusted_rho, n_positive_depth_adjusted, p_depth_adjusted, fdr_depth_adjusted
  )], 40))), collapse = "\n"),
  "\n```\n\n",
  "## Outputs\n\n",
  "- `source_data/gbmspace_extended_tam_axis_gene_scores.tsv.gz`\n",
  "- `source_data/gbmspace_extended_tam_axis_spot_scores.tsv.gz`\n",
  "- `tables/gbmspace_extended_tam_axis_gene_coverage_summary.tsv`\n",
  "- `tables/gbmspace_extended_tam_axis_section_effects.tsv`\n",
  "- `tables/gbmspace_extended_tam_axis_tumor_effects.tsv`\n",
  "- `tables/gbmspace_extended_tam_axis_spatial_support_summary.tsv`\n",
  "- `tables/gbmspace_extended_tam_axis_leave_one_tumor_out.tsv`\n",
  "- `tables/gbmspace_extended_tam_axis_closure_verdict.tsv`\n\n",
  file = file.path(out_dir, "README.md")
)

message_ts("DONE")
