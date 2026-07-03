#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3L)

input_file <- normalizePath(args[[1]], mustWork = TRUE)
gene_set_file <- normalizePath(args[[2]], mustWork = TRUE)
output_file <- args[[3]]
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
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

decode_obs <- function(file, fields, h5_index) {
  as.data.table(setNames(
    lapply(fields, function(field) {
      path <- paste0("/obs/", field)
      if (!field %in% h5_index$name[h5_index$group == "/obs"]) {
        stop("Missing obs field: ", field)
      }
      if (any(h5_index$group == path & h5_index$name == "categories")) {
        read_categorical_path(file, path)
      } else {
        h5read(file, path)
      }
    }),
    fields
  ))
}

message("Reading frozen gene sets")
gene_sets <- readRDS(gene_set_file)
neftel_sets <- gene_sets$neftel[c("AC", "OPC", "NPC1", "NPC2", "MES1", "MES2")]
pathway_sets <- gene_sets$pathways
analysis_sets <- c(neftel_sets, pathway_sets)
target_gene_names <- unique(c("LAP3", unlist(analysis_sets, use.names = FALSE)))

feature_ids <- h5read(input_file, "/var/_index")
feature_names <- read_categorical_path(input_file, "/var/feature_name")
feature_table <- data.table(
  feature_position = seq_along(feature_ids),
  feature_id = feature_ids,
  feature_name = feature_names
)
target_features <- feature_table[feature_name %chin% target_gene_names]
setorder(target_features, feature_position)
if (!"LAP3" %chin% target_features$feature_name) {
  stop("LAP3 is absent from Core GBmap")
}

coverage <- rbindlist(lapply(names(analysis_sets), function(set_name) {
  requested <- unique(analysis_sets[[set_name]])
  available <- intersect(requested, target_features$feature_name)
  data.table(
    gene_set = set_name,
    requested = length(requested),
    available = length(available),
    coverage = length(available) / length(requested)
  )
}))
fwrite(coverage, sub("\\.rds$", "_gene_coverage.csv", output_file))

message("Decoding selected cell metadata")
h5_index <- h5ls(input_file, recursive = TRUE, datasetinfo = FALSE)
obs_fields <- c(
  "_index", "author", "donor_id", "sample", "assay",
  "annotation_level_1", "annotation_level_2", "annotation_level_3",
  "annotation_level_4", "celltype_original", "iCNV"
)
obs <- decode_obs(input_file, obs_fields, h5_index)
obs[, cell_index := .I]
stopifnot(nrow(obs) == 338564L)

indptr <- as.numeric(h5read(input_file, "/X/indptr"))
nnz <- tail(indptr, 1L)
chunk_size <- 10000000
starts <- seq(1, nnz, by = chunk_size)
max_chunks <- suppressWarnings(as.integer(Sys.getenv("GBMAP_MAX_CHUNKS", "0")))
if (is.finite(max_chunks) && max_chunks > 0L) {
  starts <- head(starts, max_chunks)
  message("Dry-run chunk limit active: ", max_chunks)
}
feature_lookup <- integer(length(feature_ids))
feature_lookup[target_features$feature_position] <- seq_len(nrow(target_features))

i_parts <- vector("list", length(starts))
j_parts <- vector("list", length(starts))
x_parts <- vector("list", length(starts))
raw_parts <- vector("list", length(starts))

message("Scanning ", format(nnz, big.mark = ","), " sparse entries in ",
        length(starts), " chunks")
for (chunk_id in seq_along(starts)) {
  start <- starts[[chunk_id]]
  count <- min(chunk_size, nnz - start + 1)
  indices <- as.integer(h5read(
    input_file, "/X/indices", start = start, count = count
  )) + 1L
  target_column <- feature_lookup[indices]
  keep <- target_column > 0L

  if (any(keep)) {
    local_positions <- which(keep)
    global_positions_zero <- (start - 1) + local_positions - 1
    i_parts[[chunk_id]] <- as.integer(findInterval(global_positions_zero, indptr))
    j_parts[[chunk_id]] <- target_column[keep]
    x_chunk <- h5read(input_file, "/X/data", start = start, count = count)
    raw_chunk <- h5read(input_file, "/raw/X/data", start = start, count = count)
    x_parts[[chunk_id]] <- as.numeric(x_chunk[keep])
    raw_parts[[chunk_id]] <- as.numeric(raw_chunk[keep])
  }
  message(
    sprintf(
      "chunk=%d/%d scanned=%s selected=%s",
      chunk_id, length(starts),
      format(min(start + count - 1, nnz), big.mark = ","),
      format(sum(lengths(i_parts)), big.mark = ",")
    )
  )
}

i <- unlist(i_parts, use.names = FALSE)
j <- unlist(j_parts, use.names = FALSE)
x <- unlist(x_parts, use.names = FALSE)
raw_x <- unlist(raw_parts, use.names = FALSE)
stopifnot(length(i) == length(j), length(j) == length(x), length(x) == length(raw_x))

normalized <- sparseMatrix(
  i = i, j = j, x = x,
  dims = c(nrow(obs), nrow(target_features)),
  dimnames = list(obs[["_index"]], target_features$feature_name),
  giveCsparse = TRUE
)
raw <- sparseMatrix(
  i = i, j = j, x = raw_x,
  dims = c(nrow(obs), nrow(target_features)),
  dimnames = list(obs[["_index"]], target_features$feature_name),
  giveCsparse = TRUE
)

cache <- list(
  source_file = input_file,
  generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  obs = obs,
  normalized = normalized,
  raw = raw,
  target_features = target_features,
  gene_sets = analysis_sets,
  coverage = coverage,
  exclusion_rule = "author != 'Neftel2019'"
)
saveRDS(cache, output_file, compress = "gzip")

message("Cache written: ", output_file)
message("Dimensions: ", nrow(normalized), " cells x ", ncol(normalized), " genes")
message("Selected nonzero entries: ", length(x))
