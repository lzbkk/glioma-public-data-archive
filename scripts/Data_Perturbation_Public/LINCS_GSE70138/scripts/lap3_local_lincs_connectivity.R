#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(cmapR)
  library(data.table)
  library(fgsea)
  library(rhdf5)
})

root <- "Data_Perturbation_Public/LINCS_GSE70138"
gctx_path <- file.path(
  root,
  "cache/GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx"
)
metadata_dir <- file.path(root, "data")
query_dir <- "Data_Bulk_TCGA/Data_Merged/results/LAP3_CMap_CLUE_Adjusted/inputs"
out_dir <- Sys.getenv(
  "LAP3_LINCS_OUT_DIR",
  file.path(root, "results/LAP3_Local_Connectivity")
)
status_path <- Sys.getenv(
  "LAP3_LINCS_STATUS",
  file.path(root, "logs/lap3_local_connectivity.status")
)
workers <- as.integer(Sys.getenv("LAP3_LINCS_WORKERS", "16"))
chunk_size <- as.integer(Sys.getenv("LAP3_LINCS_CHUNK_SIZE", "512"))
max_profiles <- as.integer(Sys.getenv("LAP3_LINCS_MAX_PROFILES", "0"))

stopifnot(
  file.exists(gctx_path),
  dir.exists(query_dir),
  workers >= 1L,
  chunk_size >= 1L,
  max_profiles >= 0L
)

table_dir <- file.path(out_dir, "tables")
source_dir <- file.path(out_dir, "source_data")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

read_gmt_genes <- function(path) {
  fields <- strsplit(readLines(path, warn = FALSE)[1L], "\t", fixed = TRUE)[[1L]]
  stopifnot(length(fields) >= 3L)
  unique(fields[-c(1L, 2L)])
}

query_genes <- list(
  full = list(
    up = read_gmt_genes(file.path(query_dir, "lap3_adjusted_full_up.gmt")),
    down = read_gmt_genes(file.path(query_dir, "lap3_adjusted_full_down.gmt"))
  ),
  no_lap3 = list(
    up = read_gmt_genes(file.path(query_dir, "lap3_adjusted_no_lap3_up.gmt")),
    down = read_gmt_genes(file.path(query_dir, "lap3_adjusted_no_lap3_down.gmt"))
  )
)

fast_weighted_es <- function(stats, hit_idx) {
  n_genes <- length(stats)
  n_hits <- length(hit_idx)
  stopifnot(n_hits > 0L, n_hits < n_genes)

  ord <- order(stats, decreasing = TRUE)
  hit_positions <- sort(match(hit_idx, ord))
  hit_weights <- abs(stats[ord[hit_positions]])
  weight_sum <- sum(hit_weights)
  if (!is.finite(weight_sum) || weight_sum == 0) {
    hit_weights <- rep(1 / n_hits, n_hits)
  } else {
    hit_weights <- hit_weights / weight_sum
  }

  miss_fraction <- (hit_positions - seq_len(n_hits)) / (n_genes - n_hits)
  running_before <- c(0, head(cumsum(hit_weights), -1L)) - miss_fraction
  running_after <- cumsum(hit_weights) - miss_fraction
  max_positive <- max(c(0, running_after))
  min_negative <- min(c(0, running_before))

  if (max_positive >= -min_negative) max_positive else min_negative
}

validate_fast_es <- function() {
  set.seed(20260630)
  for (iteration in seq_len(100L)) {
    stats <- rnorm(257L)
    hit_idx <- sample(
      seq_along(stats),
      size = sample(10:60, size = 1L),
      replace = FALSE
    )
    ord <- order(stats, decreasing = TRUE)
    hit_positions <- sort(match(hit_idx, ord))
    observed <- fast_weighted_es(stats, hit_idx)
    expected <- fgsea::calcGseaStat(
      stats[ord],
      selectedStats = hit_positions,
      gseaParam = 1,
      scoreType = "std"
    )
    if (!isTRUE(all.equal(observed, expected, tolerance = 1e-12))) {
      stop("fast_weighted_es validation failed at iteration ", iteration)
    }
  }
}

validate_fast_es()

row_ids <- cmapR::read_gctx_ids(gctx_path, "row")
col_ids <- cmapR::read_gctx_ids(gctx_path, "col")
stopifnot(!anyDuplicated(row_ids), !anyDuplicated(col_ids))

query_index <- lapply(query_genes, function(query) {
  list(
    up = match(query$up, row_ids, nomatch = 0L),
    down = match(query$down, row_ids, nomatch = 0L)
  )
})
query_index <- lapply(query_index, function(query) {
  list(up = query$up[query$up > 0L], down = query$down[query$down > 0L])
})

query_mapping <- rbindlist(lapply(names(query_genes), function(query_name) {
  rbindlist(lapply(c("up", "down"), function(direction) {
    genes <- query_genes[[query_name]][[direction]]
    data.table(
      query = query_name,
      direction = direction,
      gene_id = genes,
      in_gctx = genes %chin% row_ids
    )
  }))
}))
fwrite(query_mapping, file.path(table_dir, "lap3_lincs_query_gene_mapping.csv"))

query_qc <- query_mapping[, .(
  input_genes = .N,
  matched_genes = sum(in_gctx),
  match_fraction = mean(in_gctx)
), by = .(query, direction)]
if (query_qc[, min(matched_genes)] < 10L) {
  stop("Fewer than 10 query genes matched the GCTX matrix")
}

sig_info_path <- file.path(
  metadata_dir,
  "GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz"
)
sig_metrics_path <- file.path(
  metadata_dir,
  "GSE70138_Broad_LINCS_sig_metrics_2017-03-06.txt.gz"
)
sig_info <- fread(cmd = paste("gzip -dc", shQuote(sig_info_path)))
sig_metrics <- fread(cmd = paste("gzip -dc", shQuote(sig_metrics_path)))
stopifnot(
  nrow(sig_info) == length(col_ids),
  nrow(sig_metrics) == length(col_ids),
  !anyDuplicated(sig_info$sig_id),
  !anyDuplicated(sig_metrics$sig_id),
  all(col_ids %chin% sig_info$sig_id),
  all(col_ids %chin% sig_metrics$sig_id)
)

metadata <- sig_info[match(col_ids, sig_id)]
metric_columns <- c(
  "sig_id", "distil_cc_q75", "distil_ss", "tas",
  "ngenes_modulated_up_lm", "ngenes_modulated_dn_lm", "distil_nsample"
)
metadata <- sig_metrics[, ..metric_columns][metadata, on = "sig_id"]
stopifnot(identical(metadata$sig_id, col_ids))

if (max_profiles > 0L) {
  profile_indices <- seq_len(min(max_profiles, length(col_ids)))
} else {
  profile_indices <- seq_along(col_ids)
}

score_profile <- function(stats) {
  output <- numeric(length(query_index) * 3L)
  names(output) <- unlist(lapply(names(query_index), function(query_name) {
    paste0(query_name, c("_es_up", "_es_down", "_wtcs"))
  }))

  for (query_name in names(query_index)) {
    es_up <- fast_weighted_es(stats, query_index[[query_name]]$up)
    es_down <- fast_weighted_es(stats, query_index[[query_name]]$down)
    wtcs <- if (sign(es_up) == sign(es_down)) 0 else (es_up - es_down) / 2
    output[paste0(query_name, c("_es_up", "_es_down", "_wtcs"))] <- c(
      es_up, es_down, wtcs
    )
  }
  output
}

score_chunks <- vector("list", ceiling(length(profile_indices) / chunk_size))
chunk_starts <- seq.int(1L, length(profile_indices), by = chunk_size)
started_at <- Sys.time()

for (chunk_number in seq_along(chunk_starts)) {
  from <- chunk_starts[chunk_number]
  to <- min(from + chunk_size - 1L, length(profile_indices))
  selected <- profile_indices[from:to]

  matrix_chunk <- rhdf5::h5read(
    gctx_path,
    "/0/DATA/0/matrix",
    index = list(NULL, selected)
  )
  if (is.null(dim(matrix_chunk))) {
    matrix_chunk <- matrix(matrix_chunk, ncol = 1L)
  }
  stopifnot(nrow(matrix_chunk) == length(row_ids))

  profile_scores <- parallel::mclapply(
    seq_len(ncol(matrix_chunk)),
    function(column_index) score_profile(matrix_chunk[, column_index]),
    mc.cores = min(workers, ncol(matrix_chunk)),
    mc.preschedule = TRUE
  )
  score_chunks[[chunk_number]] <- as.data.table(do.call(rbind, profile_scores))
  score_chunks[[chunk_number]][, profile_index := selected]

  status_line <- sprintf(
    "%s RUNNING pid=%d stage=score profiles=%d/%d chunks=%d/%d workers=%d",
    format(Sys.time(), "%F %T"),
    Sys.getpid(),
    to,
    length(profile_indices),
    chunk_number,
    length(chunk_starts),
    workers
  )
  writeLines(status_line, status_path)
  cat(status_line, "\n")
  rm(matrix_chunk, profile_scores)
  gc(verbose = FALSE)
}

scores <- rbindlist(score_chunks)
setorder(scores, profile_index)
scores <- cbind(
  metadata[scores$profile_index],
  scores[, !"profile_index"]
)

normalize_wtcs <- function(values) {
  positive_mean <- mean(values[values > 0], na.rm = TRUE)
  negative_mean <- abs(mean(values[values < 0], na.rm = TRUE))
  output <- rep(NA_real_, length(values))
  if (is.finite(positive_mean) && positive_mean > 0) {
    output[values > 0] <- values[values > 0] / positive_mean
  }
  if (is.finite(negative_mean) && negative_mean > 0) {
    output[values < 0] <- values[values < 0] / negative_mean
  }
  output[values == 0] <- 0
  output
}

for (query_name in names(query_index)) {
  wtcs_column <- paste0(query_name, "_wtcs")
  ncs_column <- paste0(query_name, "_ncs")
  scores[, (ncs_column) := normalize_wtcs(get(wtcs_column)), by = .(cell_id, pert_type)]
}

fwrite(
  scores,
  file.path(source_dir, "lap3_lincs_signature_level_scores.csv.gz"),
  compress = "gzip"
)

summarize_compounds <- function(query_name) {
  wtcs_column <- paste0(query_name, "_wtcs")
  ncs_column <- paste0(query_name, "_ncs")
  scores[pert_type == "trt_cp", .(
    n_signatures = .N,
    n_cells = uniqueN(cell_id),
    n_doses = uniqueN(pert_idose),
    median_wtcs = median(get(wtcs_column), na.rm = TRUE),
    median_ncs = median(get(ncs_column), na.rm = TRUE),
    mean_ncs = mean(get(ncs_column), na.rm = TRUE),
    ncs_q25 = quantile(get(ncs_column), 0.25, na.rm = TRUE),
    ncs_q75 = quantile(get(ncs_column), 0.75, na.rm = TRUE),
    fraction_negative = mean(get(ncs_column) < 0, na.rm = TRUE),
    fraction_positive = mean(get(ncs_column) > 0, na.rm = TRUE),
    median_tas = median(tas, na.rm = TRUE),
    max_tas = max(tas, na.rm = TRUE)
  ), by = .(pert_id, pert_iname)][
    , query := query_name
  ][
    order(median_ncs, -n_cells, -n_signatures)
  ][
    , reverse_rank := seq_len(.N)
  ]
}

compound_summary <- rbindlist(lapply(names(query_index), summarize_compounds))
fwrite(
  compound_summary,
  file.path(table_dir, "lap3_lincs_compound_summary_all.csv")
)

robust_candidates <- compound_summary[
  n_signatures >= 3L & n_cells >= 2L
][
  order(query, median_ncs, -fraction_negative, -n_cells)
][
  , robust_reverse_rank := seq_len(.N), by = query
]
fwrite(
  robust_candidates,
  file.path(table_dir, "lap3_lincs_robust_reverse_candidates.csv")
)

full_scores <- scores$full_ncs
no_lap3_scores <- scores$no_lap3_ncs
full_compounds <- compound_summary[query == "full"]
no_lap3_compounds <- compound_summary[query == "no_lap3"]
compound_compare <- merge(
  full_compounds,
  no_lap3_compounds,
  by = c("pert_id", "pert_iname"),
  suffixes = c("_full", "_no_lap3")
)
compound_compare[, median_ncs_delta := median_ncs_full - median_ncs_no_lap3]
setorder(compound_compare, median_ncs_full)
fwrite(
  compound_compare,
  file.path(table_dir, "lap3_lincs_full_vs_no_lap3_compounds.csv")
)

run_qc <- rbindlist(list(
  query_qc[, .(
    metric = paste0("genes_", query, "_", direction),
    value = as.character(matched_genes),
    detail = sprintf("%d/%d matched", matched_genes, input_genes)
  )],
  data.table(
    metric = c(
      "gctx_genes", "gctx_signatures", "profiles_scored", "workers",
      "chunk_size", "signature_ncs_spearman", "compound_median_ncs_spearman",
      "runtime_minutes"
    ),
    value = as.character(c(
      length(row_ids),
      length(col_ids),
      length(profile_indices),
      workers,
      chunk_size,
      cor(full_scores, no_lap3_scores, method = "spearman", use = "complete.obs"),
      cor(
        compound_compare$median_ncs_full,
        compound_compare$median_ncs_no_lap3,
        method = "spearman",
        use = "complete.obs"
      ),
      as.numeric(difftime(Sys.time(), started_at, units = "mins"))
    )),
    detail = ""
  )
))
fwrite(run_qc, file.path(table_dir, "lap3_lincs_run_qc.csv"))

cat(
  sprintf(
    paste0(
      "DONE profiles=%d/%d compounds=%d workers=%d runtime_minutes=%.2f ",
      "signature_spearman=%.6f compound_spearman=%.6f\n"
    ),
    length(profile_indices),
    length(col_ids),
    uniqueN(compound_summary$pert_id),
    workers,
    as.numeric(difftime(Sys.time(), started_at, units = "mins")),
    cor(full_scores, no_lap3_scores, method = "spearman", use = "complete.obs"),
    cor(
      compound_compare$median_ncs_full,
      compound_compare$median_ncs_no_lap3,
      method = "spearman",
      use = "complete.obs"
    )
  )
)

if (max_profiles > 0L) {
  writeLines(
    sprintf(
      "%s DONE pid=%d profiles=%d/%d workers=%d mode=validation",
      format(Sys.time(), "%F %T"),
      Sys.getpid(),
      length(profile_indices),
      length(col_ids),
      workers
    ),
    status_path
  )
}
