#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(data.table)
  library(dplyr)
})

setwd("/home/lzb/glioma")
data.table::setDTthreads(8)
set.seed(20260629)

object_file <- "Data_scRNA_GEO/GSE278456/Tumor_Integrated_SeuratV5.rds"
crosswalk_file <- "Data_scRNA_GEO/results/GSE278456_Tumor_Object_Audit/exports/gse278456_local_author_metadata_crosswalk.rds"
gene_set_file <- "Data_scRNA_GEO/results/LAP3_CellState_Phase0/source_data/frozen_cellstate_gene_sets.rds"
cell_cache_file <- "Data_scRNA_GEO/results/GSE278456_LAP3_CellState/source_data/gse278456_cell_neftel_scores_qc.csv.gz"

out_dir <- "Data_scRNA_GEO/results/GSE278456_LAP3_CellState/source_data"
log_dir <- "Data_scRNA_GEO/results/GSE278456_LAP3_CellState/logs"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(log_dir, "lap3_gse278456_extend_pathway_cache.log")
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message")
on.exit({
  sink(type = "message")
  sink(type = "output")
  close(log_con)
}, add = TRUE)

score_gene_sets <- function(log_expr, gene_sets) {
  scores <- matrix(
    NA_real_,
    nrow = ncol(log_expr),
    ncol = length(gene_sets),
    dimnames = list(colnames(log_expr), names(gene_sets))
  )
  for (signature in names(gene_sets)) {
    genes <- intersect(gene_sets[[signature]], rownames(log_expr))
    if (length(genes) < 5L) {
      next
    }
    x <- log_expr[genes, , drop = FALSE]
    z <- t(scale(t(as.matrix(x))))
    z[!is.finite(z)] <- 0
    scores[, signature] <- colMeans(z)
  }
  scores
}

cat("Started:", format(Sys.time()), "\n")
cat("Reading existing cell cache header...\n")
cell_cache <- fread(cell_cache_file, nThread = 8)
stopifnot("cell_id" %in% names(cell_cache))

cat("Loading GSE278456 local tumor object...\n")
obj <- readRDS(object_file)
crosswalk <- readRDS(crosswalk_file)
stopifnot(ncol(obj) == nrow(crosswalk), identical(colnames(obj), crosswalk$cell_id_local))

matched_cells <- which(crosswalk$matched_author)
stopifnot(length(matched_cells) == nrow(cell_cache))
stopifnot(identical(colnames(obj)[matched_cells], cell_cache$cell_id))

gene_sets <- readRDS(gene_set_file)
pathway_sets <- gene_sets$pathways
target_genes <- unique(unlist(pathway_sets, use.names = FALSE))

cat("Extracting normalized data layer for pathway genes...\n")
data_layer <- LayerData(obj[["RNA"]], layer = "data", fast = FALSE)
genes_present <- rownames(data_layer)
target_present <- intersect(target_genes, genes_present)
expr <- data_layer[target_present, matched_cells, drop = FALSE]
rm(data_layer, obj)
invisible(gc())

coverage <- bind_rows(lapply(names(pathway_sets), function(signature) {
  requested <- pathway_sets[[signature]]
  present <- intersect(requested, genes_present)
  data.frame(
    signature = signature,
    genes_requested = length(requested),
    genes_present = length(present),
    coverage = length(present) / length(requested),
    missing_genes = paste(setdiff(requested, present), collapse = ";"),
    stringsAsFactors = FALSE
  )
}))
fwrite(coverage, "Data_scRNA_GEO/results/GSE278456_LAP3_CellState/tables/gse278456_pathway_cache_gene_coverage.csv")
stopifnot(all(coverage$genes_present >= 10))

cat("Scoring pathway sets at cell level...\n")
pathway_scores <- as.data.frame(score_gene_sets(expr, pathway_sets), check.names = FALSE)
pathway_cache <- data.frame(
  cell_id = colnames(expr),
  pathway_scores,
  check.names = FALSE
)

out_csv <- file.path(out_dir, "gse278456_cell_pathway_scores.csv.gz")
out_rds <- file.path(out_dir, "gse278456_cell_pathway_scores.rds")
fwrite(pathway_cache, out_csv)
saveRDS(pathway_cache, out_rds)

cat("Cells scored:", nrow(pathway_cache), "\n")
cat("Pathway scores:", paste(names(pathway_sets), collapse = ", "), "\n")
cat("Wrote:", out_csv, "\n")
cat("Completed:", format(Sys.time()), "\n")
