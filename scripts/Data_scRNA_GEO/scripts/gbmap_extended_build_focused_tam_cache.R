#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 4L)

input_file <- normalizePath(args[[1]], mustWork = TRUE)
gene_set_file <- normalizePath(args[[2]], mustWork = TRUE)
submodule_file <- normalizePath(args[[3]], mustWork = TRUE)
output_file <- args[[4]]
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(rhdf5)
})
data.table::setDTthreads(16)

read_categorical_path <- function(file, path) {
  categories <- h5read(file, paste0(path, "/categories"))
  codes <- as.integer(h5read(file, paste0(path, "/codes")))
  values <- rep(NA_character_, length(codes))
  keep <- codes >= 0L
  values[keep] <- categories[codes[keep] + 1L]
  values
}

read_obs_field <- function(file, field) {
  path <- paste0("/obs/", field)
  has_categories <- !is.null(tryCatch(
    h5read(file, paste0(path, "/categories")),
    error = function(e) NULL
  ))
  if (has_categories) {
    read_categorical_path(file, path)
  } else {
    tryCatch(
      h5read(file, path),
      error = function(e) {
        stop("Missing or unreadable obs field: ", field, call. = FALSE)
      }
    )
  }
}

decode_obs <- function(file, fields) {
  out <- as.data.table(setNames(
    lapply(fields, function(field) read_obs_field(file, field)),
    fields
  ))
  missing <- setdiff(c("_index", "author", "donor_id", "annotation_level_1",
                      "annotation_level_2", "annotation_level_3"), names(out))
  if (length(missing) > 0L) {
    stop("Missing required obs fields: ", paste(missing, collapse = ", "))
  }
  out
}

axis_table <- data.table(
  direction = c(
    rep("tam_to_malignant", 10L),
    rep("malignant_to_tam", 9L)
  ),
  axis = c(
    "SPP1_CD44_ITGAV_ITGB1", "LGALS3_CD44_ITGB1",
    "TGFB1_TGFBR1_TGFBR2", "IL1B_IL1R1_IL1RAP",
    "TNF_TNFRSF1A_TNFRSF1B", "OSM_OSMR_LIFR",
    "EGF_EGFR", "HBEGF_EGFR_ERBB4",
    "VEGFA_KDR_FLT1_NRP1", "MIF_CD74_CXCR4_CD44",
    "CSF1_CSF1R", "CCL2_CCR2", "CX3CL1_CX3CR1",
    "MIF_CD74_CXCR4_CD44", "CD47_SIRPA",
    "LGALS9_HAVCR2", "SPP1_CD44_ITGAV_ITGB1",
    "TGFB1_TGFBR1_TGFBR2", "VEGFA_KDR_FLT1_NRP1"
  ),
  ligand = c(
    "SPP1", "LGALS3", "TGFB1", "IL1B", "TNF", "OSM", "EGF", "HBEGF",
    "VEGFA", "MIF", "CSF1", "CCL2", "CX3CL1", "MIF", "CD47",
    "LGALS9", "SPP1", "TGFB1", "VEGFA"
  ),
  receptors = c(
    "CD44;ITGAV;ITGB1", "CD44;ITGB1", "TGFBR1;TGFBR2",
    "IL1R1;IL1RAP", "TNFRSF1A;TNFRSF1B", "OSMR;LIFR",
    "EGFR", "EGFR;ERBB4", "KDR;FLT1;NRP1", "CD74;CXCR4;CD44",
    "CSF1R", "CCR2", "CX3CR1", "CD74;CXCR4;CD44", "SIRPA",
    "HAVCR2", "CD44;ITGAV;ITGB1", "TGFBR1;TGFBR2", "FLT1;KDR;NRP1"
  )
)

tam_marker_genes <- c(
  "PTPRC", "AIF1", "TYROBP", "C1QA", "C1QB", "C1QC", "CSF1R",
  "CD68", "CD163", "TREM2", "APOE", "LILRB4", "FCGR3A",
  "ITGAM", "CX3CR1", "TMEM119", "P2RY12"
)

message("Reading frozen LAP3-state gene sets")
gene_table <- fread(gene_set_file)
required_cols <- c("gene", "state_set")
stopifnot(all(required_cols %in% names(gene_table)))
state_sets <- lapply(
  c("LAP3_STATE_UNION", "LAP3_STATE_UNION_NO_TRANSLATION_PROTEOSTASIS"),
  function(set_name) sort(unique(gene_table[state_set == set_name, gene]))
)
names(state_sets) <- c(
  "LAP3_STATE_UNION", "LAP3_STATE_UNION_NO_TRANSLATION_PROTEOSTASIS"
)
state_sets <- lapply(state_sets, function(x) sort(setdiff(x, "LAP3")))

message("Reading LAP3-state submodule gene sets")
submodule_table <- fread(submodule_file)
if (!all(c("gene", "primary_submodule") %in% names(submodule_table))) {
  stop("Submodule file must contain columns: gene, primary_submodule")
}
submodule_sets <- split(
  submodule_table$gene,
  paste0("SUBMODULE_", make.names(submodule_table$primary_submodule))
)
submodule_sets <- lapply(submodule_sets, function(x) sort(setdiff(unique(x), "LAP3")))

gene_sets <- c(
  state_sets,
  submodule_sets,
  list(TAM_MARKER_SCORE = sort(unique(tam_marker_genes)))
)
axis_genes <- unique(c(axis_table$ligand, unlist(strsplit(axis_table$receptors, ";"))))
target_gene_names <- sort(unique(c("LAP3", axis_genes, unlist(gene_sets, use.names = FALSE))))

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
  stop("LAP3 is absent from Extended GBmap")
}

coverage <- rbindlist(lapply(names(gene_sets), function(set_name) {
  requested <- unique(gene_sets[[set_name]])
  available <- intersect(requested, target_features$feature_name)
  data.table(
    gene_set = set_name,
    requested = length(requested),
    available = length(available),
    coverage = length(available) / length(requested),
    missing_genes = paste(setdiff(requested, available), collapse = ";")
  )
}), use.names = TRUE)

axis_coverage <- rbindlist(lapply(seq_len(nrow(axis_table)), function(i) {
  receptors <- unlist(strsplit(axis_table$receptors[[i]], ";"))
  available_receptors <- intersect(receptors, target_features$feature_name)
  data.table(
    direction = axis_table$direction[[i]],
    axis = axis_table$axis[[i]],
    ligand = axis_table$ligand[[i]],
    ligand_available = axis_table$ligand[[i]] %chin% target_features$feature_name,
    receptors = axis_table$receptors[[i]],
    receptor_requested = length(receptors),
    receptor_available = length(available_receptors),
    available_receptors = paste(available_receptors, collapse = ";"),
    missing_receptors = paste(setdiff(receptors, available_receptors), collapse = ";")
  )
}), use.names = TRUE)

coverage_file <- sub("\\.rds$", "_gene_coverage.csv", output_file)
axis_coverage_file <- sub("\\.rds$", "_axis_coverage.csv", output_file)
fwrite(coverage, coverage_file)
fwrite(axis_coverage, axis_coverage_file)

message("Decoding selected cell metadata")
obs_fields <- c(
  "_index", "author", "donor_id", "sample", "assay", "suspension_type",
  "annotation_level_1", "annotation_level_2", "annotation_level_3",
  "celltype_original", "organism_ontology_term_id"
)
optional_fields <- c("sample", "assay", "suspension_type", "celltype_original",
                     "organism_ontology_term_id")
required_fields <- setdiff(obs_fields, optional_fields)
obs <- decode_obs(input_file, required_fields)
for (field in optional_fields) {
  obs[[field]] <- tryCatch(
    read_obs_field(input_file, field),
    error = function(e) {
      message("Optional obs field unavailable: ", field)
      rep(NA_character_, nrow(obs))
    }
  )
}
obs[, cell_index := .I]

indptr <- as.numeric(h5read(input_file, "/X/indptr", bit64conversion = "double"))
nnz <- tail(indptr, 1L)
chunk_size <- as.integer(Sys.getenv("GBMAP_CHUNK_SIZE", "10000000"))
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
  indices <- as.integer(h5read(input_file, "/X/indices", start = start, count = count)) + 1L
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

  selected <- sum(lengths(i_parts))
  message(sprintf(
    "[%d/%d] scanned=%s selected=%s",
    chunk_id, length(starts),
    format(min(start + count - 1, nnz), big.mark = ","),
    format(selected, big.mark = ",")
  ))
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
  gene_sets = gene_sets,
  axis_table = axis_table,
  coverage = coverage,
  axis_coverage = axis_coverage,
  exclusion_rule = "author != 'Neftel2019'",
  primary_unit = "author::donor_id",
  note = "Focused lightweight cache for Extended GBmap malignant-TAM feasibility analysis."
)
saveRDS(cache, output_file, compress = "gzip")

message("Cache written: ", output_file)
message("Coverage written: ", coverage_file)
message("Axis coverage written: ", axis_coverage_file)
message("Dimensions: ", nrow(normalized), " cells x ", ncol(normalized), " genes")
message("Selected nonzero entries: ", length(x))
