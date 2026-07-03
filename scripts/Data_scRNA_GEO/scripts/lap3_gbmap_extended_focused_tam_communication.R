#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)

cache_file <- normalizePath(args[[1]], mustWork = TRUE)
out_dir <- args[[2]]
table_dir <- file.path(out_dir, "tables")
source_dir <- file.path(out_dir, "source_data")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})
data.table::setDTthreads(16)

mean_sparse_rows <- function(mat, genes, min_genes = 3L) {
  genes <- intersect(genes, colnames(mat))
  if (length(genes) < min_genes) {
    return(rep(NA_real_, nrow(mat)))
  }
  Matrix::rowMeans(mat[, genes, drop = FALSE])
}

zscore <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) {
    return(rep(NA_real_, length(x)))
  }
  (x - mean(x, na.rm = TRUE)) / s
}

safe_spearman <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 8L || length(unique(x[ok])) < 3L || length(unique(y[ok])) < 3L) {
    return(list(rho = NA_real_, p = NA_real_, n = sum(ok)))
  }
  test <- suppressWarnings(stats::cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  list(rho = unname(test$estimate), p = test$p.value, n = sum(ok))
}

prefix_columns <- function(dt, prefix, keep) {
  old <- setdiff(names(dt), keep)
  setnames(dt, old, paste0(prefix, old))
  dt
}

message("Reading cache: ", cache_file)
cache <- readRDS(cache_file)
obs <- copy(cache$obs)
expr <- cache$normalized
raw <- cache$raw
axis_table <- as.data.table(cache$axis_table)
gene_sets <- cache$gene_sets

obs[, author_donor := paste(author, donor_id, sep = "::")]
obs[, compartment := fifelse(
  annotation_level_1 == "Neoplastic", "malignant",
  fifelse(
    annotation_level_1 == "Non-neoplastic" &
      annotation_level_2 == "Myeloid" &
      annotation_level_3 %chin% c("TAM-BDM", "TAM-MG"),
    "tam",
    NA_character_
  )
)]
obs[, included_for_primary := !is.na(compartment) & author != "Neftel2019"]

entry_counts <- obs[, .(
  cells = .N,
  donors = uniqueN(author_donor),
  authors = uniqueN(author)
), by = .(author, annotation_level_1, annotation_level_2, annotation_level_3, compartment)]
setorder(entry_counts, author, -cells)
fwrite(entry_counts, file.path(table_dir, "extended_focused_tam_entry_counts.csv"))

keep_cells <- which(obs$included_for_primary)
if (length(keep_cells) == 0L) {
  stop("No cells passed focused malignant/TAM inclusion rules")
}
obs_keep <- obs[keep_cells]
expr_keep <- expr[keep_cells, , drop = FALSE]
raw_keep <- raw[keep_cells, , drop = FALSE]

score_sets <- c(
  "LAP3_STATE_UNION",
  "LAP3_STATE_UNION_NO_TRANSLATION_PROTEOSTASIS",
  grep("^SUBMODULE_", names(gene_sets), value = TRUE),
  "TAM_MARKER_SCORE"
)
score_sets <- intersect(score_sets, names(gene_sets))
score_dt <- data.table(cell_index = seq_len(nrow(obs_keep)))
for (set_name in score_sets) {
  score_dt[[set_name]] <- mean_sparse_rows(expr_keep, gene_sets[[set_name]], min_genes = 3L)
}
if ("LAP3" %chin% colnames(expr_keep)) {
  score_dt[, LAP3_expr := as.numeric(expr_keep[, "LAP3"])]
  score_dt[, LAP3_detected := as.numeric(raw_keep[, "LAP3"] > 0)]
}

comm_genes <- sort(unique(c(axis_table$ligand, unlist(strsplit(axis_table$receptors, ";")))))
comm_genes <- intersect(comm_genes, colnames(expr_keep))
message("Summarising ", length(comm_genes), " curated ligand/receptor genes")
comm_expr <- as.data.table(as.matrix(expr_keep[, comm_genes, drop = FALSE]))
setnames(comm_expr, comm_genes, paste0(comm_genes, "_expr"))
comm_detect <- as.data.table(as.matrix(raw_keep[, comm_genes, drop = FALSE] > 0))
setnames(comm_detect, comm_genes, paste0(comm_genes, "_detect"))

cell_dt <- cbind(
  obs_keep[, .(author, donor_id, author_donor, assay, suspension_type,
               annotation_level_1, annotation_level_2, annotation_level_3,
               compartment)],
  score_dt[, !"cell_index"],
  comm_expr,
  comm_detect
)
rm(comm_expr, comm_detect, score_dt)
gc(verbose = FALSE)

score_cols <- c(
  score_sets,
  "LAP3_expr", "LAP3_detected",
  paste0(comm_genes, "_expr"), paste0(comm_genes, "_detect")
)
score_cols <- intersect(score_cols, names(cell_dt))

compartment_scores <- cell_dt[, c(
  .(
    cells = .N,
    assays = paste(sort(unique(assay)), collapse = ";"),
    suspension_types = paste(sort(unique(suspension_type)), collapse = ";")
  ),
  lapply(.SD, function(x) mean(x, na.rm = TRUE))
), by = .(author, donor_id, author_donor, compartment), .SDcols = score_cols]
setorder(compartment_scores, author, donor_id, compartment)
fwrite(compartment_scores, file.path(table_dir, "extended_focused_tam_compartment_scores.csv"))
fwrite(compartment_scores, file.path(source_dir, "extended_focused_tam_compartment_scores.csv.gz"))

mal <- compartment_scores[compartment == "malignant"]
tam <- compartment_scores[compartment == "tam"]
mal <- prefix_columns(copy(mal), "mal__", keep = c("author", "donor_id", "author_donor"))
tam <- prefix_columns(copy(tam), "tam__", keep = c("author", "donor_id", "author_donor"))
paired <- merge(mal, tam, by = c("author", "donor_id", "author_donor"))
paired[, tam_fraction_within_pair := tam__cells / (tam__cells + mal__cells)]

eligible <- paired[mal__cells >= 20 & tam__cells >= 20]
fwrite(eligible, file.path(table_dir, "extended_focused_tam_donor_compartment_summary.csv"))

axis_rows <- list()
response_map <- list(
  tam_to_malignant = c(
    "mal__LAP3_STATE_UNION",
    "mal__LAP3_STATE_UNION_NO_TRANSLATION_PROTEOSTASIS",
    "mal__LAP3_expr"
  ),
  malignant_to_tam = c(
    "tam_fraction_within_pair",
    "tam__TAM_MARKER_SCORE"
  )
)

for (i in seq_len(nrow(axis_table))) {
  direction <- axis_table$direction[[i]]
  axis <- axis_table$axis[[i]]
  ligand <- axis_table$ligand[[i]]
  receptors <- unlist(strsplit(axis_table$receptors[[i]], ";"))

  if (direction == "tam_to_malignant") {
    sender_prefix <- "tam__"
    receiver_prefix <- "mal__"
  } else {
    sender_prefix <- "mal__"
    receiver_prefix <- "tam__"
  }

  ligand_expr_col <- paste0(sender_prefix, ligand, "_expr")
  ligand_detect_col <- paste0(sender_prefix, ligand, "_detect")
  receptor_expr_cols <- paste0(receiver_prefix, receptors, "_expr")
  receptor_detect_cols <- paste0(receiver_prefix, receptors, "_detect")
  receptor_expr_cols <- intersect(receptor_expr_cols, names(eligible))
  receptor_detect_cols <- intersect(receptor_detect_cols, names(eligible))

  if (!ligand_expr_col %chin% names(eligible) || length(receptor_expr_cols) == 0L) {
    next
  }

  sender_ligand <- eligible[[ligand_expr_col]]
  receiver_receptor <- rowMeans(as.matrix(eligible[, ..receptor_expr_cols]), na.rm = TRUE)
  axis_score <- (zscore(sender_ligand) + zscore(receiver_receptor)) / 2
  ligand_detect_mean <- if (ligand_detect_col %chin% names(eligible)) {
    mean(eligible[[ligand_detect_col]], na.rm = TRUE)
  } else {
    NA_real_
  }
  receptor_detect_mean <- if (length(receptor_detect_cols) > 0L) {
    mean(rowMeans(as.matrix(eligible[, ..receptor_detect_cols]), na.rm = TRUE), na.rm = TRUE)
  } else {
    NA_real_
  }

  for (response in intersect(response_map[[direction]], names(eligible))) {
    cor_result <- safe_spearman(axis_score, eligible[[response]])
    axis_rows[[length(axis_rows) + 1L]] <- data.table(
      direction = direction,
      axis = axis,
      ligand = ligand,
      receptors_requested = paste(receptors, collapse = ";"),
      receptors_available = paste(sub(paste0("^", receiver_prefix), "", sub("_expr$", "", receptor_expr_cols)), collapse = ";"),
      response = response,
      n_author_donors = cor_result$n,
      n_authors = uniqueN(eligible$author[is.finite(axis_score) & is.finite(eligible[[response]])]),
      ligand_detect_mean = ligand_detect_mean,
      receptor_detect_mean = receptor_detect_mean,
      spearman_rho = cor_result$rho,
      p_value = cor_result$p
    )
  }
}

axis_results <- rbindlist(axis_rows, use.names = TRUE, fill = TRUE)
if (nrow(axis_results) > 0L) {
  axis_results[, fdr_bh := p.adjust(p_value, method = "BH"), by = .(direction, response)]
  axis_results[, abs_rho := abs(spearman_rho)]
  setorder(axis_results, direction, response, fdr_bh, -abs_rho)
}
tam_to_malignant <- axis_results[direction == "tam_to_malignant"]
malignant_to_tam <- axis_results[direction == "malignant_to_tam"]
fwrite(tam_to_malignant, file.path(table_dir, "extended_focused_tam_axes_tam_to_malignant.csv"))
fwrite(malignant_to_tam, file.path(table_dir, "extended_focused_tam_axes_malignant_to_tam.csv"))

primary_summary <- data.table(
  metric = c(
    "source_cache", "generated_at", "cells_in_cache", "primary_cells_used",
    "paired_author_donors_ge20_ge20", "paired_authors_ge20_ge20",
    "paired_unique_donors_ge20_ge20",
    "tam_to_malignant_axes_tested", "malignant_to_tam_axes_tested",
    "best_tam_to_malignant_axis",
    "best_tam_to_malignant_response",
    "best_tam_to_malignant_rho",
    "best_tam_to_malignant_fdr",
    "best_malignant_to_tam_axis",
    "best_malignant_to_tam_response",
    "best_malignant_to_tam_rho",
    "best_malignant_to_tam_fdr"
  ),
  value = NA_character_
)
best_t2m <- tam_to_malignant[is.finite(fdr_bh)][order(fdr_bh, -abs_rho)][1]
best_m2t <- malignant_to_tam[is.finite(fdr_bh)][order(fdr_bh, -abs_rho)][1]
summary_values <- list(
  cache_file,
  format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  as.character(nrow(cache$obs)),
  as.character(nrow(obs_keep)),
  as.character(nrow(eligible)),
  as.character(uniqueN(eligible$author)),
  as.character(uniqueN(eligible$donor_id)),
  as.character(uniqueN(tam_to_malignant$axis)),
  as.character(uniqueN(malignant_to_tam$axis)),
  if (nrow(best_t2m)) best_t2m$axis else NA_character_,
  if (nrow(best_t2m)) best_t2m$response else NA_character_,
  if (nrow(best_t2m)) sprintf("%.4f", best_t2m$spearman_rho) else NA_character_,
  if (nrow(best_t2m)) sprintf("%.4g", best_t2m$fdr_bh) else NA_character_,
  if (nrow(best_m2t)) best_m2t$axis else NA_character_,
  if (nrow(best_m2t)) best_m2t$response else NA_character_,
  if (nrow(best_m2t)) sprintf("%.4f", best_m2t$spearman_rho) else NA_character_,
  if (nrow(best_m2t)) sprintf("%.4g", best_m2t$fdr_bh) else NA_character_
)
primary_summary[, value := unlist(summary_values)]
fwrite(primary_summary, file.path(table_dir, "extended_focused_tam_primary_summary.csv"))

readme <- c(
  "# Extended GBmap Focused Malignant-TAM Communication Feasibility",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Purpose",
  "",
  "This module tests a focused, conservative malignant-cell/TAM co-variation layer in Extended GBmap. It is not treated as causal cell-cell communication. It asks whether curated TAM-to-malignant or malignant-to-TAM ligand/receptor axes are jointly detectable and co-vary with donor-level LAP3-state summaries.",
  "",
  "## Primary Design",
  "",
  "- Exclude Neftel2019 to avoid overlap with GSE131928-derived analyses.",
  "- Use author::donor_id as the inference unit.",
  "- Restrict primary cells to Neoplastic malignant cells and TAM-BDM/TAM-MG myeloid cells.",
  "- Require at least 20 malignant cells and 20 TAM cells per author::donor_id pair.",
  "- Summarise expression and detection at donor-compartment level before testing correlations.",
  "- Control FDR separately by direction and response.",
  "",
  "## Outputs",
  "",
  "- `tables/extended_focused_tam_entry_counts.csv`: entry and compartment counts.",
  "- `tables/extended_focused_tam_compartment_scores.csv`: donor-compartment source table.",
  "- `tables/extended_focused_tam_donor_compartment_summary.csv`: paired malignant/TAM donor summaries.",
  "- `tables/extended_focused_tam_axes_tam_to_malignant.csv`: TAM-to-malignant axis tests.",
  "- `tables/extended_focused_tam_axes_malignant_to_tam.csv`: malignant-to-TAM axis tests.",
  "- `tables/extended_focused_tam_primary_summary.csv`: compact summary for project tracking.",
  "- `source_data/extended_focused_tam_compartment_scores.csv.gz`: source data export.",
  "",
  "## Interpretation Boundary",
  "",
  "A positive result supports a spatial/ecological compatibility hypothesis for LAP3-state tumors and TAM-rich microenvironments. A negative or weak result does not refute the main LAP3-state story because this module is limited by cross-study integration, cell annotation heterogeneity, donor-level aggregation, and lack of spatial adjacency."
)
writeLines(readme, file.path(out_dir, "README.md"))

message("Analysis complete: ", out_dir)
