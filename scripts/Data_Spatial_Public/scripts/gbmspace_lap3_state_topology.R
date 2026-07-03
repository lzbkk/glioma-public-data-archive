#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(rhdf5)
  library(FNN)
})

project_root <- "/home/lzb/glioma"
args <- commandArgs(trailingOnly = TRUE)
max_files <- if (length(args) >= 1L) as.integer(args[[1]]) else NA_integer_
h5ad_root <- file.path(project_root, "Data_Spatial_Public/GBM_Space/cache/all_h5ad")
cache_dir <- file.path(project_root, "Data_Spatial_Public/GBM_Space/results/Lightweight_Cache")
state_gene_file <- file.path(project_root, "Data_Bulk_TCGA/Data_Merged/results/LAP3_State_Module/tables/lap3_state_frozen_gene_sets.csv")
out_dir <- file.path(project_root, "Data_Spatial_Public/GBM_Space/results/LAP3_State_Spatial_Topology")
tables_dir <- file.path(out_dir, "tables")
source_dir <- file.path(out_dir, "source_data")
logs_dir <- file.path(out_dir, "logs")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

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

parse_sample_name <- function(x) {
  data.table(
    sample_name = x,
    tumor_id = sub("^((AT[0-9]+).*)$", "\\2", x)
  )
}

wilcox_vs_zero <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 4 || all(x == 0)) return(NA_real_)
  suppressWarnings(wilcox.test(x, mu = 0, exact = FALSE)$p.value)
}

safe_cor <- function(x, y, method = "spearman") {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 20 || uniqueN(x[keep]) < 3 || uniqueN(y[keep]) < 3) {
    return(NA_real_)
  }
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

classify_readout <- function(x) {
  lx <- tolower(x)
  fifelse(
    grepl("HALLMARK_MTORC1|MTORC1_READOUT|LEUCINE_BCAA|REACTOME_TRANSLATION", x),
    "pathway",
    fifelse(
      grepl("tam|monocyte|dendritic|immune", lx),
      "myeloid_tam",
      fifelse(
        grepl("ac|opc|npc|proliferative|hypoxic|dev.like", lx),
        "malignant_state",
        fifelse(
          grepl("niche|gliosis|vasculature|grey.matter|white.matter", lx),
          "spatial_niche",
          fifelse(
            grepl("necrosis|pseudopalisading|cellular.tumor|leading.edge|microvascular|infiltrating|blood.vessels|perinecrotic", lx),
            "histopath",
            "other"
          )
        )
      )
    )
  )
}

priority_readout <- function(x) {
  lx <- tolower(x)
  grepl(
    paste(c(
      "HALLMARK_MTORC1", "LEUCINE_BCAA", "MTORC1_READOUT", "REACTOME_TRANSLATION",
      "immune..tams", "immune..resident", "resident.tams", "pro.inflammatory.tams",
      "anti.inflammatory.tams", "angiogenic.tams", "proliferative.tams", "monocytes",
      "hypoxic.1", "dev.like..ac", "dev.like..opc", "proliferative",
      "cellular.tumor", "necrosis", "pseudopalisading", "microvascular",
      "vasculature", "gliosis", "leading.edge", "infiltrating.tumor"
    ), collapse = "|"),
    lx
  )
}

score_state_sets_for_h5ad <- function(file, state_sets) {
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
  in_tissue <- as.integer(h5read(file, "obs/in_tissue"))
  array_col <- as.integer(h5read(file, "obs/array_col"))
  array_row <- as.integer(h5read(file, "obs/array_row"))
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

  lap3_idx <- match("LAP3", var_upper)
  state_indices <- lapply(state_sets, function(g) which(var_upper %in% g & gene_feature))
  state_sizes <- vapply(state_indices, length, integer(1))
  state_bool <- lapply(state_indices, function(idx) {
    x <- logical(n_var)
    x[idx] <- TRUE
    x
  })

  lap3_raw <- numeric(n_obs)
  lap3_log1p_cp10k <- numeric(n_obs)
  gene_library_size <- numeric(n_obs)
  detected_gene_features <- integer(n_obs)
  score_matrix <- matrix(
    0,
    nrow = n_obs,
    ncol = length(state_indices),
    dimnames = list(NULL, names(state_indices))
  )

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
    if (length(lap3_hit) > 0) lap3_raw[[i]] <- vals[lap3_hit[[1]]]

    if (lib > 0) {
      if (lap3_raw[[i]] > 0) {
        lap3_log1p_cp10k[[i]] <- log1p(lap3_raw[[i]] / lib * 10000)
      }
      log_norm <- log1p(vals / lib * 10000)
      for (sig in names(state_bool)) {
        hit <- state_bool[[sig]][idx]
        if (any(hit) && state_sizes[[sig]] > 0L) {
          score_matrix[i, sig] <- sum(log_norm[hit]) / state_sizes[[sig]]
        }
      }
    }
  }

  parsed <- parse_sample_name(sample_name)
  spot_id <- paste(sample_id, obs_index, sep = "__")
  score_dt <- data.table(
    spot_id = spot_id,
    barcode = obs_index,
    h5ad_sample_id = sample_id,
    sample_name = sample_name,
    tumor_id = parsed$tumor_id,
    in_tissue = in_tissue,
    array_row = array_row,
    array_col = array_col,
    spatial_x = spatial_x,
    spatial_y = spatial_y,
    LAP3_raw = lap3_raw,
    LAP3_detected = lap3_raw > 0,
    LAP3_log1p_cp10k = lap3_log1p_cp10k,
    gene_library_size = gene_library_size,
    detected_gene_features = detected_gene_features
  )
  score_dt <- cbind(score_dt, as.data.table(score_matrix))

  coverage <- rbindlist(lapply(names(state_indices), function(sig) {
    data.table(
      h5ad_sample_id = sample_id,
      state_set = sig,
      n_genes_requested = length(state_sets[[sig]]),
      n_genes_present = state_sizes[[sig]],
      present_fraction = state_sizes[[sig]] / length(state_sets[[sig]])
    )
  }))

  list(scores = score_dt, coverage = coverage)
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

state_gene_sets <- fread(state_gene_file)
state_gene_sets <- unique(state_gene_sets[
  state_set %in% c("LAP3_STATE_UNION", "LAP3_STATE_UNION_NO_TRANSLATION_PROTEOSTASIS") &
    toupper(gene) != "LAP3",
  .(state_set, gene = toupper(gene))
])
state_sets <- split(state_gene_sets$gene, state_gene_sets$state_set)
stopifnot(length(state_sets) == 2L)

section_summary <- fread(file.path(cache_dir, "tables/gbmspace_section_summary.tsv"))
pathway_scores <- fread(file.path(cache_dir, "tables/gbmspace_spot_lap3_pathway_scores.tsv"))
pathway_cols <- c("HALLMARK_MTORC1_SIGNALING", "LEUCINE_BCAA_CORE", "MTORC1_READOUT_CORE", "REACTOME_TRANSLATION")
pathway_scores <- pathway_scores[, c("spot_id", "h5ad_sample_id", pathway_cols), with = FALSE]

h5ad_files <- list.files(file.path(h5ad_root, "anndata"), pattern = "\\.h5ad$", full.names = TRUE)
h5ad_files <- sort(h5ad_files)
if (!is.na(max_files)) h5ad_files <- head(h5ad_files, max_files)
if (length(h5ad_files) == 0L) stop("No H5AD files found under ", file.path(h5ad_root, "anndata"), call. = FALSE)

section_effect_rows <- list()
section_contrast_rows <- list()
state_score_rows <- list()
coverage_rows <- list()
readout_catalog_rows <- list()

for (file_i in seq_along(h5ad_files)) {
  file <- h5ad_files[[file_i]]
  sample_id <- sub("\\.h5ad$", "", basename(file))
  message_ts(sprintf("[%d/%d] %s", file_i, length(h5ad_files), sample_id))

  score_obj <- score_state_sets_for_h5ad(file, state_sets)
  sec_score <- score_obj$scores[in_tissue == 1 | is.na(in_tissue)]
  coverage_rows[[sample_id]] <- score_obj$coverage
  state_score_rows[[sample_id]] <- sec_score[, .(
    spot_id, h5ad_sample_id, tumor_id,
    LAP3_log1p_cp10k, LAP3_STATE_UNION, LAP3_STATE_UNION_NO_TRANSLATION_PROTEOSTASIS,
    gene_library_size, detected_gene_features
  )]

  ng_file <- section_summary[h5ad_sample_id == sample_id, non_gene_feature_file][1]
  if (is.na(ng_file) || !file.exists(ng_file)) next

  ng <- fread(ng_file)
  sec <- merge(sec_score, pathway_scores, by = c("spot_id", "h5ad_sample_id"), all.x = TRUE)
  sec <- merge(sec, ng, by = c("spot_id", "h5ad_sample_id"), all.x = TRUE)
  if (nrow(sec) < 20) next

  readout_cols <- setdiff(names(sec), c(
    "spot_id", "barcode", "h5ad_sample_id", "sample_name", "tumor_id", "in_tissue",
    "array_row", "array_col", "spatial_x", "spatial_y",
    "LAP3_raw", "LAP3_detected", "LAP3_log1p_cp10k",
    "gene_library_size", "detected_gene_features",
    names(state_sets)
  ))
  readout_cols <- readout_cols[vapply(sec[, ..readout_cols], is.numeric, logical(1))]
  if (length(readout_cols) == 0L) next

  readout_catalog <- data.table(readout = readout_cols)
  readout_catalog[, readout_class := classify_readout(readout)]
  readout_catalog[, priority := priority_readout(readout) | readout_class %in% c("pathway", "myeloid_tam", "spatial_niche", "histopath")]
  readout_catalog_rows[[sample_id]] <- readout_catalog

  neighbor_cols <- readout_catalog[priority == TRUE, readout]
  neighbor_dt <- make_neighbor_means(sec, neighbor_cols, k = 6L)
  sec <- cbind(sec, neighbor_dt)

  for (state in names(state_sets)) {
    q25 <- as.numeric(quantile(sec[[state]], 0.25, na.rm = TRUE))
    q75 <- as.numeric(quantile(sec[[state]], 0.75, na.rm = TRUE))
    high <- sec[[state]] >= q75
    low <- sec[[state]] <= q25

    for (ro in readout_cols) {
      cls <- readout_catalog[readout == ro, readout_class][1]
      pr <- readout_catalog[readout == ro, priority][1]
      section_effect_rows[[length(section_effect_rows) + 1L]] <- data.table(
        h5ad_sample_id = sample_id,
        tumor_id = sec$tumor_id[[1]],
        state_set = state,
        readout = ro,
        readout_class = cls,
        priority = pr,
        analysis_type = "same_spot_colocalization",
        n_spots = sum(is.finite(sec[[state]]) & is.finite(sec[[ro]])),
        raw_rho = safe_cor(sec[[state]], sec[[ro]]),
        depth_adjusted_rho = residual_rank_cor(sec, state, ro)
      )

      section_contrast_rows[[length(section_contrast_rows) + 1L]] <- data.table(
        h5ad_sample_id = sample_id,
        tumor_id = sec$tumor_id[[1]],
        state_set = state,
        readout = ro,
        readout_class = cls,
        priority = pr,
        contrast_type = "same_spot_top25_vs_bottom25",
        n_high = sum(high, na.rm = TRUE),
        n_low = sum(low, na.rm = TRUE),
        mean_high = mean(sec[[ro]][high], na.rm = TRUE),
        mean_low = mean(sec[[ro]][low], na.rm = TRUE),
        delta_high_minus_low = mean(sec[[ro]][high], na.rm = TRUE) - mean(sec[[ro]][low], na.rm = TRUE)
      )
    }

    for (ro in neighbor_cols) {
      ncol <- paste0("neighbor_", ro)
      cls <- readout_catalog[readout == ro, readout_class][1]
      pr <- readout_catalog[readout == ro, priority][1]
      section_effect_rows[[length(section_effect_rows) + 1L]] <- data.table(
        h5ad_sample_id = sample_id,
        tumor_id = sec$tumor_id[[1]],
        state_set = state,
        readout = ro,
        readout_class = cls,
        priority = pr,
        analysis_type = "neighbor_topology_k6",
        n_spots = sum(is.finite(sec[[state]]) & is.finite(sec[[ncol]])),
        raw_rho = safe_cor(sec[[state]], sec[[ncol]]),
        depth_adjusted_rho = residual_rank_cor(sec, state, ncol)
      )

      section_contrast_rows[[length(section_contrast_rows) + 1L]] <- data.table(
        h5ad_sample_id = sample_id,
        tumor_id = sec$tumor_id[[1]],
        state_set = state,
        readout = ro,
        readout_class = cls,
        priority = pr,
        contrast_type = "neighbor_top25_vs_bottom25_k6",
        n_high = sum(high, na.rm = TRUE),
        n_low = sum(low, na.rm = TRUE),
        mean_high = mean(sec[[ncol]][high], na.rm = TRUE),
        mean_low = mean(sec[[ncol]][low], na.rm = TRUE),
        delta_high_minus_low = mean(sec[[ncol]][high], na.rm = TRUE) - mean(sec[[ncol]][low], na.rm = TRUE)
      )
    }
  }
}

section_effects <- rbindlist(section_effect_rows, use.names = TRUE, fill = TRUE)
section_contrasts <- rbindlist(section_contrast_rows, use.names = TRUE, fill = TRUE)
state_scores <- rbindlist(state_score_rows, use.names = TRUE, fill = TRUE)
coverage <- rbindlist(coverage_rows, use.names = TRUE, fill = TRUE)
readout_catalog <- unique(rbindlist(readout_catalog_rows, use.names = TRUE, fill = TRUE))

tumor_effects <- section_effects[
  is.finite(depth_adjusted_rho),
  .(
    n_sections = .N,
    median_raw_rho = median(raw_rho, na.rm = TRUE),
    median_depth_adjusted_rho = median(depth_adjusted_rho, na.rm = TRUE)
  ),
  by = .(tumor_id, state_set, readout, readout_class, priority, analysis_type)
]

summary <- tumor_effects[
  ,
  .(
    n_tumors = .N,
    median_raw_rho = median(median_raw_rho, na.rm = TRUE),
    median_depth_adjusted_rho = median(median_depth_adjusted_rho, na.rm = TRUE),
    n_positive_depth_adjusted = sum(median_depth_adjusted_rho > 0, na.rm = TRUE),
    p_depth_adjusted = wilcox_vs_zero(median_depth_adjusted_rho)
  ),
  by = .(state_set, readout, readout_class, priority, analysis_type)
]
summary[, fdr_depth_adjusted_all := p.adjust(p_depth_adjusted, method = "BH"), by = .(state_set, analysis_type)]
summary[, fdr_depth_adjusted_priority := {
  if (isTRUE(priority[1])) p.adjust(p_depth_adjusted, method = "BH") else rep(NA_real_, .N)
}, by = .(state_set, analysis_type, priority)]
summary[, abs_median_depth_adjusted_rho := abs(median_depth_adjusted_rho)]
setorder(summary, state_set, analysis_type, fdr_depth_adjusted_priority, fdr_depth_adjusted_all, -abs_median_depth_adjusted_rho)

tumor_contrasts <- section_contrasts[
  is.finite(delta_high_minus_low),
  .(
    n_sections = .N,
    median_delta = median(delta_high_minus_low, na.rm = TRUE)
  ),
  by = .(tumor_id, state_set, readout, readout_class, priority, contrast_type)
]
contrast_summary <- tumor_contrasts[
  ,
  .(
    n_tumors = .N,
    median_tumor_delta = median(median_delta, na.rm = TRUE),
    n_positive_delta = sum(median_delta > 0, na.rm = TRUE),
    p_delta = wilcox_vs_zero(median_delta)
  ),
  by = .(state_set, readout, readout_class, priority, contrast_type)
]
contrast_summary[, fdr_delta_all := p.adjust(p_delta, method = "BH"), by = .(state_set, contrast_type)]
contrast_summary[, fdr_delta_priority := {
  if (isTRUE(priority[1])) p.adjust(p_delta, method = "BH") else rep(NA_real_, .N)
}, by = .(state_set, contrast_type, priority)]
contrast_summary[, abs_median_tumor_delta := abs(median_tumor_delta)]
setorder(contrast_summary, state_set, contrast_type, fdr_delta_priority, fdr_delta_all, -abs_median_tumor_delta)

loto_rows <- list()
for (state in unique(summary$state_set)) {
  for (atype in unique(summary$analysis_type)) {
    keep_summary <- summary[state_set == state & analysis_type == atype & priority == TRUE]
    for (ro in keep_summary$readout) {
      te <- tumor_effects[state_set == state & analysis_type == atype & readout == ro]
      for (drop_tumor in unique(te$tumor_id)) {
        sub <- te[tumor_id != drop_tumor]
        loto_rows[[length(loto_rows) + 1L]] <- data.table(
          state_set = state,
          analysis_type = atype,
          readout = ro,
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

coverage_summary <- coverage[
  ,
  .(
    n_sections = .N,
    min_present_fraction = min(present_fraction),
    median_present_fraction = median(present_fraction),
    n_sections_below_0_8 = sum(present_fraction < 0.8)
  ),
  by = state_set
]

fwrite(section_effects, file.path(tables_dir, "gbmspace_lap3_state_section_effects.tsv"), sep = "\t")
fwrite(tumor_effects, file.path(tables_dir, "gbmspace_lap3_state_tumor_effects.tsv"), sep = "\t")
fwrite(summary, file.path(tables_dir, "gbmspace_lap3_state_topology_summary.tsv"), sep = "\t")
fwrite(section_contrasts, file.path(tables_dir, "gbmspace_lap3_state_section_contrasts.tsv"), sep = "\t")
fwrite(tumor_contrasts, file.path(tables_dir, "gbmspace_lap3_state_tumor_contrasts.tsv"), sep = "\t")
fwrite(contrast_summary, file.path(tables_dir, "gbmspace_lap3_state_contrast_summary.tsv"), sep = "\t")
fwrite(loto, file.path(tables_dir, "gbmspace_lap3_state_leave_one_tumor_out.tsv"), sep = "\t")
fwrite(coverage, file.path(tables_dir, "gbmspace_lap3_state_gene_coverage_by_section.tsv"), sep = "\t")
fwrite(coverage_summary, file.path(tables_dir, "gbmspace_lap3_state_gene_coverage_summary.tsv"), sep = "\t")
fwrite(readout_catalog, file.path(tables_dir, "gbmspace_lap3_state_readout_catalog.tsv"), sep = "\t")
fwrite(state_scores, file.path(source_dir, "gbmspace_spot_lap3_state_scores.tsv.gz"), sep = "\t")

top_priority <- summary[
  state_set == "LAP3_STATE_UNION" &
    priority == TRUE &
    is.finite(fdr_depth_adjusted_priority)
][order(analysis_type, fdr_depth_adjusted_priority, -abs(median_depth_adjusted_rho))]
top_priority <- head(top_priority, 30)

top_contrast <- contrast_summary[
  state_set == "LAP3_STATE_UNION" &
    priority == TRUE &
    is.finite(fdr_delta_priority)
][order(contrast_type, fdr_delta_priority, -abs(median_tumor_delta))]
top_contrast <- head(top_contrast, 30)

cat(
  "# GBM-Space LAP3-State Spatial Topology\n\n",
  "生成时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n",
  "## 输入\n\n",
  "- H5AD root: `", h5ad_root, "`\n",
  "- Lightweight cache: `", cache_dir, "`\n",
  "- LAP3-state gene sets: `", state_gene_file, "`\n",
  "- Sections processed: `", uniqueN(state_scores$h5ad_sample_id), "`\n",
  "- Tumors processed: `", uniqueN(state_scores$tumor_id), "`\n",
  "- In-tissue spots scored: `", nrow(state_scores), "`\n\n",
  "## 方法概述\n\n",
  "`LAP3_STATE_UNION` 和 `LAP3_STATE_UNION_NO_TRANSLATION_PROTEOSTASIS` 以每个 spot 的 mean log1p(CP10K) 计算，LAP3 本身已从 gene set 排除。\n",
  "同一 spot co-localization 使用 section 内 Spearman 和 rank-residual correlation；残差化变量为 log1p(gene_library_size) 与 detected_gene_features。\n",
  "Spatial topology 使用每个 spot 的 6-nearest spatial neighbors 的 readout 平均值，再执行同样的 depth/complexity-adjusted residual correlation。\n",
  "所有主推断先在 section 内得到效应量，再按 tumor 取 median，并在 tumor-level effects 上进行 Wilcoxon signed-rank test。\n\n",
  "## 关键输出\n\n",
  "- `tables/gbmspace_lap3_state_topology_summary.tsv`\n",
  "- `tables/gbmspace_lap3_state_contrast_summary.tsv`\n",
  "- `tables/gbmspace_lap3_state_leave_one_tumor_out.tsv`\n",
  "- `tables/gbmspace_lap3_state_gene_coverage_summary.tsv`\n",
  "- `source_data/gbmspace_spot_lap3_state_scores.tsv.gz`\n\n",
  "## Gene coverage\n\n",
  "```text\n",
  paste(capture.output(print(coverage_summary)), collapse = "\n"),
  "\n```\n\n",
  "## Top priority depth-adjusted effects\n\n",
  "```text\n",
  paste(capture.output(print(top_priority[, .(
    state_set, analysis_type, readout, readout_class, n_tumors,
    median_depth_adjusted_rho, n_positive_depth_adjusted,
    p_depth_adjusted, fdr_depth_adjusted_priority
  )])), collapse = "\n"),
  "\n```\n\n",
  "## Top priority high-vs-low contrasts\n\n",
  "```text\n",
  paste(capture.output(print(top_contrast[, .(
    state_set, contrast_type, readout, readout_class, n_tumors,
    median_tumor_delta, n_positive_delta, p_delta, fdr_delta_priority
  )])), collapse = "\n"),
  "\n```\n\n",
  "## 解释边界\n\n",
  "本分析检验的是 LAP3-state spatial co-localization / neighbor topology，而不是因果机制。\n",
  "若同一 spot 和 neighbor topology 均通过 depth-adjusted tumor-level sensitivity，可支持 Figure 4 spatial niche/state panel；\n",
  "若仅 raw 或少数 readout 显著，应降级为 spatial context 或补充材料。\n",
  file = file.path(out_dir, "README.md")
)

saveRDS(
  list(
    generated_at = Sys.time(),
    state_gene_file = state_gene_file,
    h5ad_root = h5ad_root,
    coverage_summary = coverage_summary,
    summary = summary,
    contrast_summary = contrast_summary
  ),
  file.path(out_dir, "gbmspace_lap3_state_topology_results.rds")
)

message_ts("DONE")
