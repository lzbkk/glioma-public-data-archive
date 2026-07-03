#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(ggplot2)
  library(rhdf5)
  library(cluster)
})

project_dir <- "/home/lzb/glioma"
setwd(project_dir)
data.table::setDTthreads(12)
set.seed(20260701)

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(prefix, default = NA_character_) {
  hit <- args[startsWith(args, prefix)]
  if (!length(hit)) return(default)
  sub(prefix, "", hit[[1]])
}
max_spatial_files <- suppressWarnings(as.integer(get_arg("--max-spatial-files=", NA_character_)))
skip_spatial <- "--skip-spatial" %in% args

out_dir <- "Data_Bulk_TCGA/Data_Merged/results/LAP3_Malignant_State_Module_Audit"
table_dir <- file.path(out_dir, "tables")
source_dir <- file.path(out_dir, "source_data")
plot_dir <- file.path(out_dir, "plots")
export_dir <- file.path(out_dir, "exports")
log_dir <- file.path(out_dir, "logs")
for (d in c(table_dir, source_dir, plot_dir, export_dir, log_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

log_file <- file.path(log_dir, "lap3_malignant_state_module_audit.log")
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message")
on.exit({
  sink(type = "message")
  sink(type = "output")
  close(log_con)
}, add = TRUE)

message_ts <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n")
  flush.console()
}

write_table <- function(x, filename, dir = table_dir) {
  fwrite(as.data.table(x), file.path(dir, filename))
}

save_plot <- function(plot, filename, width = 7, height = 5) {
  ggsave(file.path(plot_dir, paste0(filename, ".pdf")), plot, width = width, height = height)
  ggsave(file.path(plot_dir, paste0(filename, ".png")), plot, width = width, height = height, dpi = 300)
  ggsave(file.path(plot_dir, paste0(filename, ".tiff")), plot, width = width, height = height, dpi = 300)
  if (requireNamespace("svglite", quietly = TRUE)) {
    ggsave(file.path(plot_dir, paste0(filename, ".svg")), plot, width = width, height = height)
  }
}

collapse_expr <- function(mat) {
  mat <- as.matrix(mat)
  storage.mode(mat) <- "double"
  rn <- rownames(mat)
  if (is.null(rn)) stop("Matrix rownames are required")
  keep <- !is.na(rn) & nzchar(rn)
  mat <- mat[keep, , drop = FALSE]
  rn <- rn[keep]
  if (!anyDuplicated(rn)) return(mat)
  dt <- as.data.table(mat)
  dt[, gene := rn]
  dt <- dt[, lapply(.SD, mean, na.rm = TRUE), by = gene]
  out <- as.matrix(dt[, -"gene"])
  rownames(out) <- dt$gene
  storage.mode(out) <- "double"
  out
}

z_rows <- function(mat) {
  z <- t(scale(t(mat)))
  z[!is.finite(z)] <- NA_real_
  z
}

z_num <- function(x) {
  x <- as.numeric(x)
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  as.numeric((x - mean(x, na.rm = TRUE)) / s)
}

module_score_dense <- function(expr_log2, genes, min_genes = 3L) {
  genes <- intersect(unique(genes), rownames(expr_log2))
  if (length(genes) < min_genes) return(rep(NA_real_, ncol(expr_log2)))
  z <- z_rows(expr_log2[genes, , drop = FALSE])
  colMeans(z, na.rm = TRUE)
}

module_score_sparse_cells <- function(expr, genes, min_genes = 3L) {
  genes <- intersect(unique(genes), colnames(expr))
  if (length(genes) < min_genes) return(rep(NA_real_, nrow(expr)))
  Matrix::rowMeans(expr[, genes, drop = FALSE])
}

spearman_safe <- function(x, y, min_n = 25L) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) < min_n || length(unique(x)) < 3L || length(unique(y)) < 3L) {
    return(list(n = length(x), rho = NA_real_, p_value = NA_real_))
  }
  test <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  list(n = length(x), rho = unname(test$estimate), p_value = test$p.value)
}

safe_cor <- function(x, y, method = "spearman", min_n = 20L) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < min_n || uniqueN(x[keep]) < 3L || uniqueN(y[keep]) < 3L) return(NA_real_)
  suppressWarnings(cor(x[keep], y[keep], method = method))
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

residual_rank_cor <- function(dt, xcol, ycol) {
  keep <- is.finite(dt[[xcol]]) &
    is.finite(dt[[ycol]]) &
    is.finite(dt$gene_library_size) &
    is.finite(dt$detected_gene_features)
  if (sum(keep) < 20 || uniqueN(dt[[xcol]][keep]) < 3L || uniqueN(dt[[ycol]][keep]) < 3L) {
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
  if (length(x) < 4L || all(x == 0)) return(NA_real_)
  suppressWarnings(wilcox.test(x, mu = 0, exact = FALSE)$p.value)
}

map_state <- function(annotation_level_3, annotation_level_4) {
  state <- rep(NA_character_, length(annotation_level_3))
  state[annotation_level_3 == "AC-like" | grepl("^AC-like", annotation_level_4)] <- "AC"
  state[annotation_level_3 == "OPC-like" | grepl("^OPC-like", annotation_level_4)] <- "OPC"
  state[annotation_level_3 == "NPC-like" | grepl("^NPC-like", annotation_level_4)] <- "NPC"
  state[annotation_level_3 == "MES-like" | grepl("^MES-like", annotation_level_4)] <- "MES"
  state
}

reference_sets <- list(
  MES_LIKE = c("CHI3L1", "CD44", "VIM", "ANXA1", "ANXA2", "LGALS3", "EMP3", "TNC", "SERPINE1", "TIMP1", "ITGA5", "FN1"),
  AC_LIKE = c("GFAP", "AQP4", "SLC1A3", "ALDH1L1", "CLU", "FABP7", "SPARCL1", "S100B", "APOE"),
  OPC_LIKE = c("PDGFRA", "OLIG1", "OLIG2", "CSPG4", "BCAN", "SOX10", "PTPRZ1", "SOX6"),
  NPC_LIKE = c("SOX4", "SOX11", "DLL3", "HES6", "ASCL1", "DCX", "STMN2", "TUBB3"),
  PROLIFERATION = c("MKI67", "TOP2A", "PCNA", "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "CDK1", "CCNB1", "CCNB2", "AURKA", "AURKB", "BUB1", "BUB1B", "UBE2C"),
  HYPOXIA_STRESS = c("CA9", "VEGFA", "SLC2A1", "LDHA", "ENO1", "PGK1", "BNIP3", "NDRG1", "P4HA1", "EGLN3", "ADM", "ANGPTL4", "HSPA5", "CALR"),
  MYELOID_TAM = c("AIF1", "CD68", "LST1", "CSF1R", "TYROBP", "FCER1G", "C1QA", "C1QB", "C1QC", "CD163", "MRC1", "TREM2", "APOE")
)

message_ts("Started LAP3 malignant-state module audit")
message_ts("Args:", paste(args, collapse = " "))

message_ts("Reading frozen submodule object")
sub_obj <- readRDS("Data_Bulk_TCGA/Data_Merged/results/LAP3_State_Submodules/exports/lap3_state_submodules_projection.rds")
gene_assignment <- as.data.table(sub_obj$gene_assignment)
malignant_genes <- sort(unique(sub_obj$submodule_sets$LAP3_MALIGNANT_STATE_MODULE))
if (length(malignant_genes) != 144L) {
  stop("Expected 144 malignant-state genes, got ", length(malignant_genes))
}

message_ts("Reading bulk expression matrices")
tcga_tpm_raw <- readRDS("Data_Bulk_TCGA/Data_Merged/data_analysis/expr_tpm_glioma_uni.rds")
tcga_tpm <- as.matrix(tcga_tpm_raw[, -1, drop = FALSE])
rownames(tcga_tpm) <- rownames(tcga_tpm_raw)
tcga_log2 <- log2(collapse_expr(tcga_tpm) + 1)

cgga_validation <- as.data.table(readRDS("Data_Bulk_CGGA/results/LAP3_CGGA/exports/cgga_lap3_validation_dataset.rds"))
cgga693_mat <- as.matrix(readRDS("Data_Bulk_CGGA/mRNAseq_693/tpm_data.rds"))
cgga693_log2 <- log2(collapse_expr(cgga693_mat) + 1)
cgga325_raw <- fread("Data_Bulk_CGGA/mRNAseq_325/CGGA.mRNAseq_325.RSEM-genes.20200506.txt", data.table = FALSE, check.names = FALSE)
cgga325_mat <- as.matrix(cgga325_raw[, -1, drop = FALSE])
rownames(cgga325_mat) <- cgga325_raw$Gene_Name
cgga325_log2 <- log2(collapse_expr(cgga325_mat) + 1)

bulk_mats <- list(
  TCGA = tcga_log2,
  CGGA_mRNAseq_693 = cgga693_log2,
  CGGA_mRNAseq_325 = cgga325_log2
)

message_ts("Clustering malignant-state genes from bulk expression")
z_blocks <- lapply(names(bulk_mats), function(ds) {
  mat <- bulk_mats[[ds]]
  out <- matrix(NA_real_, nrow = length(malignant_genes), ncol = ncol(mat),
                dimnames = list(malignant_genes, paste(ds, colnames(mat), sep = "::")))
  present <- intersect(malignant_genes, rownames(mat))
  out[present, ] <- z_rows(mat[present, , drop = FALSE])
  out
})
combined_z <- do.call(cbind, z_blocks)
gene_cor <- suppressWarnings(cor(t(combined_z), use = "pairwise.complete.obs", method = "pearson"))
gene_cor[!is.finite(gene_cor)] <- 0
diag(gene_cor) <- 1
dist_obj <- as.dist(1 - pmax(pmin(gene_cor, 1), -1))
hc <- hclust(dist_obj, method = "average")

silhouette_table <- rbindlist(lapply(2:4, function(k) {
  cl <- cutree(hc, k = k)
  size <- as.integer(table(cl))
  sil <- tryCatch(cluster::silhouette(cl, dist_obj), error = function(e) NULL)
  data.table(
    k = k,
    min_cluster_size = min(size),
    max_cluster_size = max(size),
    mean_silhouette = if (is.null(sil)) NA_real_ else mean(sil[, "sil_width"])
  )
}))
eligible <- silhouette_table[min_cluster_size >= 5L & is.finite(mean_silhouette)]
chosen_k <- if (nrow(eligible)) eligible[which.max(mean_silhouette), k] else 4L
clusters_raw <- cutree(hc, k = chosen_k)
names(clusters_raw) <- rownames(gene_cor)

merge_small_clusters <- function(cl, cor_mat, min_size = 5L) {
  nms <- names(cl)
  cl <- as.character(cl)
  names(cl) <- nms
  repeat {
    sizes <- sort(table(cl))
    if (!length(sizes) || min(sizes) >= min_size || length(sizes) <= 1L) break
    small <- names(sizes)[[1]]
    small_genes <- names(cl)[cl == small]
    candidates <- setdiff(unique(cl), small)
    if (!length(small_genes) || !length(candidates)) break
    target_scores <- vapply(candidates, function(candidate) {
      target_genes <- names(cl)[cl == candidate]
      mean(cor_mat[small_genes, target_genes, drop = FALSE], na.rm = TRUE)
    }, numeric(1))
    target_scores[!is.finite(target_scores)] <- -Inf
    if (all(!is.finite(target_scores)) || all(target_scores == -Inf)) {
      target <- names(sort(table(cl), decreasing = TRUE))[[1]]
    } else {
      target <- candidates[which.max(target_scores)]
    }
    cl[small_genes] <- target
  }
  old_levels <- sort(unique(cl))
  new_levels <- paste0("M", seq_along(old_levels))
  setNames(new_levels[match(cl, old_levels)], names(cl))
}

clusters <- merge_small_clusters(setNames(clusters_raw, names(clusters_raw)), gene_cor, min_size = 5L)
chosen_k <- uniqueN(clusters)

cluster_sets <- split(names(clusters), clusters)
cluster_sets <- lapply(cluster_sets, sort)
write_table(silhouette_table, "malignant_module_cluster_silhouette_selection.csv")

cluster_gene_table <- data.table(
  gene = names(clusters),
  malignant_cluster = unname(clusters)
)
cluster_gene_table <- merge(cluster_gene_table, gene_assignment, by = "gene", all.x = TRUE)
cluster_gene_table[, cluster_size := .N, by = malignant_cluster]

cluster_ref_overlap <- rbindlist(lapply(names(cluster_sets), function(cl) {
  genes <- cluster_sets[[cl]]
  rbindlist(lapply(names(reference_sets), function(ref) {
    ref_genes <- reference_sets[[ref]]
    overlap <- intersect(genes, ref_genes)
    data.table(
      malignant_cluster = cl,
      reference_set = ref,
      cluster_genes = length(genes),
      reference_genes = length(ref_genes),
      overlap_n = length(overlap),
      overlap_genes = paste(overlap, collapse = ";"),
      overlap_fraction_cluster = length(overlap) / length(genes)
    )
  }))
}))
write_table(cluster_ref_overlap, "malignant_module_cluster_reference_overlap.csv")

cluster_labels <- cluster_ref_overlap[
  order(malignant_cluster, -overlap_n, -overlap_fraction_cluster)
][, .SD[1], by = malignant_cluster][
  , .(malignant_cluster, cluster_label = paste0(malignant_cluster, "_", reference_set))
]
cluster_gene_table <- merge(cluster_gene_table, cluster_labels, by = "malignant_cluster", all.x = TRUE)
write_table(cluster_gene_table[order(malignant_cluster, gene)], "malignant_module_gene_clusters.csv")

message_ts("Scoring bulk cluster modules")
bulk_score_list <- lapply(names(bulk_mats), function(ds) {
  mat <- bulk_mats[[ds]]
  projection <- if (ds == "TCGA") {
    fread("Data_Bulk_TCGA/Data_Merged/results/LAP3_State_Submodules/tables/tcga_lap3_state_submodule_projection.csv")
  } else {
    fread("Data_Bulk_TCGA/Data_Merged/results/LAP3_State_Submodules/tables/cgga_lap3_state_submodule_projection.csv")[dataset == ds]
  }
  sample_col <- if (ds == "TCGA") "barcode" else "sample_id"
  out <- copy(projection)
  out[, analysis_dataset := ds]
  out[, sample_key := get(sample_col)]
  for (cl in names(cluster_sets)) {
    out[, (paste0("LAP3_MALIGNANT_", cl)) := module_score_dense(mat, cluster_sets[[cl]])[match(sample_key, colnames(mat))]]
  }
  out
})
bulk_scores <- rbindlist(bulk_score_list, fill = TRUE)

comp_file <- "Data_Bulk_TCGA/Data_Merged/results/LAP3_Bulk_Composition_Audit/source_data/bulk_composition_audit_analysis_dataset.csv"
if (file.exists(comp_file)) {
  comp <- fread(comp_file)
  comp_keep <- c(
    "analysis_dataset", "sample_key", "IMMUNE_PANLEUKOCYTE_CORE", "STROMAL_FIBROVASCULAR_CORE",
    "PURITY_LIKE_LOW_IMMUNE_STROMAL_PROXY", "TAM_MYELOID_EXTENDED_CORE",
    "HYPOXIA_EXTENDED_CORE", "PROLIFERATION_EXTENDED_CORE",
    "NEURAL_GLIAL_TUMOR_CONTENT_PROXY"
  )
  comp_keep <- comp_keep[comp_keep %in% names(comp)]
  bulk_scores <- merge(bulk_scores, comp[, ..comp_keep], by = c("analysis_dataset", "sample_key"), all.x = TRUE)
}
fwrite(bulk_scores, file.path(source_dir, "bulk_malignant_cluster_scores.csv"))

cluster_score_cols <- paste0("LAP3_MALIGNANT_", names(cluster_sets))
bulk_variables <- c(
  "LAP3_log2_expr", "LAP3_STATE_UNION", "LAP3_MALIGNANT_STATE_MODULE",
  "LAP3_MYELOID_TAM_CONTEXT_MODULE", "LAP3_ANABOLIC_TRANSLATION_MODULE",
  "LAP3_PROTEOSTASIS_STRESS_MODULE", "LAP3_HYPOXIA_PERINECROTIC_MODULE",
  "IMMUNE_PANLEUKOCYTE_CORE", "STROMAL_FIBROVASCULAR_CORE",
  "PURITY_LIKE_LOW_IMMUNE_STROMAL_PROXY", "TAM_MYELOID_EXTENDED_CORE",
  "HYPOXIA_EXTENDED_CORE", "PROLIFERATION_EXTENDED_CORE"
)
bulk_variables <- bulk_variables[bulk_variables %in% names(bulk_scores)]

make_bulk_strata <- function(d) {
  out <- list(data.table(stratum = "all", row_id = seq_len(nrow(d))))
  if ("cohort" %in% names(d) && any(as.character(d$cohort) == "GBM", na.rm = TRUE)) {
    out <- c(out, list(data.table(stratum = "cohort_GBM", row_id = which(as.character(d$cohort) == "GBM"))))
    out <- c(out, list(data.table(stratum = "cohort_LGG", row_id = which(as.character(d$cohort) == "LGG"))))
    idh <- as.character(d$idh_status)
    out <- c(out, list(data.table(stratum = "GBM_IDH_WT_or_Wildtype", row_id = which(as.character(d$cohort) == "GBM" & idh %in% c("WT", "Wildtype")))))
  }
  if ("tumor_class" %in% names(d)) {
    vals <- sort(unique(na.omit(as.character(d$tumor_class))))
    out <- c(out, lapply(vals, function(v) data.table(stratum = paste0("tumor_class_", v), row_id = which(as.character(d$tumor_class) == v))))
    idh <- as.character(d$idh_status)
    out <- c(out, list(data.table(stratum = "GBM_grade4_IDH_Wildtype", row_id = which(as.character(d$tumor_class) == "GBM_grade4" & idh %in% c("WT", "Wildtype")))))
  }
  rbindlist(out, fill = TRUE)
}

bulk_cor <- rbindlist(lapply(split(bulk_scores, bulk_scores$analysis_dataset), function(d0) {
  d0 <- as.data.table(d0)
  ds <- unique(d0$analysis_dataset)[1]
  strata <- make_bulk_strata(d0)
  rbindlist(lapply(unique(strata$stratum), function(st) {
    d <- d0[strata[stratum == st, row_id]]
    rbindlist(lapply(cluster_score_cols, function(cl_col) {
      rbindlist(lapply(bulk_variables, function(v) {
        res <- spearman_safe(d[[cl_col]], d[[v]], min_n = 25L)
        data.table(dataset = ds, stratum = st, cluster_score = cl_col, variable = v, n = res$n, rho = res$rho, p_value = res$p_value)
      }))
    }))
  }))
}), fill = TRUE)
bulk_cor[, p_adj_BH := p.adjust(p_value, method = "BH")]
write_table(bulk_cor, "bulk_malignant_cluster_score_correlations.csv")

message_ts("Scoring Core GBmap cluster modules from lightweight cache")
gbmap_cache <- readRDS("Data_scRNA_GEO/GBmap_Core/cache/core_gbmap_lap3_state_union_lightweight.rds")
obs <- as.data.table(gbmap_cache$obs)
expr <- gbmap_cache$normalized
raw <- gbmap_cache$raw
obs[, author_donor := paste(author, donor_id, sep = "::")]
obs[, author_state := map_state(annotation_level_3, annotation_level_4)]
obs[, main_entry := author != "Neftel2019" & annotation_level_1 == "Neoplastic" & !is.na(author_state)]
obs[, strict_entry := main_entry & iCNV == "aneuploid"]
main_idx <- which(obs$main_entry)
expr_main <- expr[main_idx, , drop = FALSE]
raw_main <- raw[main_idx, , drop = FALSE]
obs_main <- copy(obs[main_idx])

gbmap_scores <- obs_main[, .(
  cell_id = `_index`,
  author,
  donor_id,
  sample,
  assay,
  annotation_level_3,
  annotation_level_4,
  iCNV,
  author_donor,
  author_state,
  strict_entry
)]
gbmap_scores[, lap3_log_norm := as.numeric(expr_main[, "LAP3"])]
gbmap_scores[, lap3_raw := as.numeric(raw_main[, "LAP3"])]
gbmap_scores[, lap3_detected := lap3_raw > 0]
gbmap_scores[, LAP3_STATE_UNION := module_score_sparse_cells(expr_main, setdiff(gbmap_cache$gene_sets$LAP3_STATE_UNION, "LAP3"), min_genes = 20L)]
gbmap_scores[, LAP3_MALIGNANT_STATE_MODULE := module_score_sparse_cells(expr_main, malignant_genes, min_genes = 20L)]
for (cl in names(cluster_sets)) {
  gbmap_scores[, (paste0("LAP3_MALIGNANT_", cl)) := module_score_sparse_cells(expr_main, cluster_sets[[cl]])]
}
fwrite(gbmap_scores, file.path(source_dir, "gbmap_core_malignant_cluster_cell_scores.csv.gz"))

gbmap_donor_state <- gbmap_scores[
  , c(
    .(n_cells = .N, lap3_mean = mean(lap3_log_norm, na.rm = TRUE), lap3_detected_fraction = mean(lap3_detected, na.rm = TRUE)),
    lapply(.SD, mean, na.rm = TRUE)
  ),
  by = .(author, donor_id, author_donor, author_state, strict_entry),
  .SDcols = c("LAP3_STATE_UNION", "LAP3_MALIGNANT_STATE_MODULE", cluster_score_cols)
]
fwrite(gbmap_donor_state, file.path(source_dir, "gbmap_core_malignant_cluster_donor_state_scores.csv"))

gbmap_assoc <- rbindlist(lapply(c(FALSE, TRUE), function(strict_flag) {
  d0 <- gbmap_donor_state[strict_entry == strict_flag & n_cells >= 20]
  entry <- if (strict_flag) "strict_aneuploid" else "main_neoplastic"
  rbindlist(lapply(sort(unique(d0$author_state)), function(st) {
    d <- d0[author_state == st]
    rbindlist(lapply(cluster_score_cols, function(cl_col) {
      res_lap3 <- spearman_safe(d[[cl_col]], d$lap3_mean, min_n = 6L)
      res_union <- spearman_safe(d[[cl_col]], d$LAP3_STATE_UNION, min_n = 6L)
      rbindlist(list(
        data.table(entry_variant = entry, author_state = st, cluster_score = cl_col, target = "lap3_mean", n = res_lap3$n, rho = res_lap3$rho, p_value = res_lap3$p_value),
        data.table(entry_variant = entry, author_state = st, cluster_score = cl_col, target = "LAP3_STATE_UNION", n = res_union$n, rho = res_union$rho, p_value = res_union$p_value)
      ))
    }))
  }))
}), fill = TRUE)
gbmap_assoc[, p_adj_BH := p.adjust(p_value, method = "BH")]
write_table(gbmap_assoc, "gbmap_core_malignant_cluster_donor_state_associations.csv")

score_state_sets_for_h5ad <- function(file, gene_sets) {
  sample_id <- sub("\\.h5ad$", "", basename(file))
  attrs <- h5readAttributes(file, "X")
  shape <- as.integer(attrs$shape)
  n_obs <- shape[[1]]
  n_var <- shape[[2]]
  var_names <- as.character(h5read(file, "var/_index"))
  feature_types <- read_categorical(file, "var/feature_types", n_expected = length(var_names))
  var_upper <- toupper(var_names)
  gene_feature <- feature_types == "Gene Expression"

  obs_index <- as.character(h5read(file, "obs/_index"))
  sample_name <- read_categorical(file, "obs/sample_name", n_expected = n_obs)
  tumor_id <- sub("^((AT[0-9]+).*)$", "\\2", sample_name)

  data <- as.numeric(h5read(file, "X/data"))
  indices <- as.integer(h5read(file, "X/indices")) + 1L
  indptr <- as.integer(h5read(file, "X/indptr"))

  lap3_idx <- match("LAP3", var_upper)
  set_indices <- lapply(gene_sets, function(g) which(var_upper %in% toupper(g) & gene_feature))
  set_sizes <- vapply(set_indices, length, integer(1))
  set_bool <- lapply(set_indices, function(idx) {
    x <- logical(n_var)
    x[idx] <- TRUE
    x
  })

  lap3_raw <- numeric(n_obs)
  lap3_log1p_cp10k <- numeric(n_obs)
  gene_library_size <- numeric(n_obs)
  detected_gene_features <- integer(n_obs)
  score_matrix <- matrix(0, nrow = n_obs, ncol = length(set_indices), dimnames = list(NULL, names(set_indices)))

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
    if (length(lap3_hit) > 0L) lap3_raw[[i]] <- vals[lap3_hit[[1]]]
    if (lib > 0) {
      if (lap3_raw[[i]] > 0) lap3_log1p_cp10k[[i]] <- log1p(lap3_raw[[i]] / lib * 10000)
      log_norm <- log1p(vals / lib * 10000)
      for (sig in names(set_bool)) {
        hit <- set_bool[[sig]][idx]
        if (any(hit) && set_sizes[[sig]] > 0L) {
          score_matrix[i, sig] <- sum(log_norm[hit]) / set_sizes[[sig]]
        }
      }
    }
  }

  score_dt <- data.table(
    spot_id = paste(sample_id, obs_index, sep = "__"),
    h5ad_sample_id = sample_id,
    tumor_id = tumor_id,
    LAP3_raw = lap3_raw,
    LAP3_detected = lap3_raw > 0,
    LAP3_log1p_cp10k = lap3_log1p_cp10k,
    gene_library_size = gene_library_size,
    detected_gene_features = detected_gene_features
  )
  score_dt <- cbind(score_dt, as.data.table(score_matrix))
  coverage <- rbindlist(lapply(names(set_indices), function(sig) {
    data.table(
      h5ad_sample_id = sample_id,
      cluster_score = sig,
      n_genes_requested = length(gene_sets[[sig]]),
      n_genes_present = set_sizes[[sig]],
      present_fraction = set_sizes[[sig]] / length(gene_sets[[sig]])
    )
  }))
  list(scores = score_dt, coverage = coverage)
}

gbmspace_summary <- data.table()
gbmspace_coverage <- data.table()
if (!skip_spatial) {
  message_ts("Scoring GBM-Space H5AD sections for malignant clusters")
  spatial_manifest <- readRDS("Data_Spatial_Public/GBM_Space/results/Lightweight_Cache/gbmspace_lightweight_cache_manifest.rds")
  files <- spatial_manifest$section_summary$file
  files <- files[file.exists(files)]
  if (is.finite(max_spatial_files)) files <- head(files, max_spatial_files)
  spatial_gene_sets <- c(list(LAP3_MALIGNANT_STATE_MODULE = malignant_genes), setNames(cluster_sets, cluster_score_cols))
  spatial_scores <- vector("list", length(files))
  spatial_cov <- vector("list", length(files))
  for (i in seq_along(files)) {
    message_ts("Spatial", i, "/", length(files), basename(files[[i]]))
    tmp <- score_state_sets_for_h5ad(files[[i]], spatial_gene_sets)
    spatial_scores[[i]] <- tmp$scores
    spatial_cov[[i]] <- tmp$coverage
  }
  gbmspace_spot <- rbindlist(spatial_scores, fill = TRUE)
  gbmspace_coverage <- rbindlist(spatial_cov, fill = TRUE)
  state_spot <- fread("Data_Spatial_Public/GBM_Space/results/LAP3_State_Spatial_Topology/source_data/gbmspace_spot_lap3_state_scores.tsv.gz",
                      select = c("spot_id", "LAP3_STATE_UNION", "LAP3_STATE_UNION_NO_TRANSLATION_PROTEOSTASIS"))
  gbmspace_spot <- merge(gbmspace_spot, state_spot, by = "spot_id", all.x = TRUE)
  fwrite(gbmspace_spot, file.path(source_dir, "gbmspace_spot_malignant_cluster_scores.tsv.gz"), sep = "\t")
  write_table(gbmspace_coverage, "gbmspace_malignant_cluster_gene_coverage_by_section.tsv")

  section_cor <- rbindlist(lapply(split(gbmspace_spot, gbmspace_spot$h5ad_sample_id), function(sec) {
    sec <- as.data.table(sec)
    rbindlist(lapply(c("LAP3_log1p_cp10k", "LAP3_STATE_UNION"), function(target) {
      rbindlist(lapply(c("LAP3_MALIGNANT_STATE_MODULE", cluster_score_cols), function(cl_col) {
        data.table(
          h5ad_sample_id = unique(sec$h5ad_sample_id),
          tumor_id = unique(sec$tumor_id)[1],
          target = target,
          cluster_score = cl_col,
          n_spots = nrow(sec),
          raw_rho = safe_cor(sec[[cl_col]], sec[[target]], method = "spearman"),
          depth_adjusted_rho = residual_rank_cor(sec, cl_col, target)
        )
      }))
    }))
  }), fill = TRUE)
  write_table(section_cor, "gbmspace_malignant_cluster_section_correlations.tsv")
  gbmspace_summary <- section_cor[, .(
    n_sections = .N,
    n_tumors = uniqueN(tumor_id),
    median_raw_rho = median(raw_rho, na.rm = TRUE),
    median_depth_adjusted_rho = median(depth_adjusted_rho, na.rm = TRUE),
    n_positive_depth_adjusted = sum(depth_adjusted_rho > 0, na.rm = TRUE),
    p_depth_adjusted = wilcox_vs_zero(depth_adjusted_rho)
  ), by = .(target, cluster_score)]
  gbmspace_summary[, p_adj_BH := p.adjust(p_depth_adjusted, method = "BH")]
  write_table(gbmspace_summary, "gbmspace_malignant_cluster_spatial_summary.tsv")
}

message_ts("Building summary plots")
gene_order <- cluster_gene_table[order(malignant_cluster, -tcga_t, gene), gene]
tcga_projection <- fread("Data_Bulk_TCGA/Data_Merged/results/LAP3_State_Submodules/tables/tcga_lap3_state_submodule_projection.csv")
gene_feature_cor <- rbindlist(lapply(malignant_genes, function(g) {
  x <- if (g %in% rownames(tcga_log2)) as.numeric(tcga_log2[g, match(tcga_projection$barcode, colnames(tcga_log2))]) else rep(NA_real_, nrow(tcga_projection))
  rbindlist(lapply(c("LAP3_log2_expr", "LAP3_STATE_UNION", "TAM_MYELOID_CORE", "HALLMARK_HYPOXIA", "PROLIFERATION_CORE"), function(v) {
    res <- spearman_safe(x, tcga_projection[[v]], min_n = 50L)
    data.table(gene = g, feature = v, rho = res$rho, p_value = res$p_value)
  }))
}))
gene_feature_cor <- merge(gene_feature_cor, cluster_gene_table[, .(gene, malignant_cluster, cluster_label)], by = "gene", all.x = TRUE)
write_table(gene_feature_cor, "tcga_malignant_gene_feature_correlations.csv")

plot_gene_cor <- copy(gene_feature_cor)
plot_gene_cor[, gene := factor(gene, levels = gene_order)]
plot_gene_cor[, feature := factor(feature, levels = c("LAP3_log2_expr", "LAP3_STATE_UNION", "TAM_MYELOID_CORE", "HALLMARK_HYPOXIA", "PROLIFERATION_CORE"))]
p_gene <- ggplot(plot_gene_cor, aes(feature, gene, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.1) +
  facet_grid(malignant_cluster ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C", midpoint = 0, limits = c(-1, 1), na.value = "grey90") +
  labs(x = NULL, y = NULL, fill = "rho") +
  theme_bw(base_size = 7) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.grid = element_blank(), strip.text.y = element_text(angle = 0))
save_plot(p_gene, "tcga_malignant_gene_feature_correlation_heatmap", width = 7.2, height = 10)

bulk_plot <- bulk_cor[
  stratum == "all" &
    variable %in% c("LAP3_log2_expr", "LAP3_STATE_UNION", "TAM_MYELOID_EXTENDED_CORE", "HYPOXIA_EXTENDED_CORE", "PROLIFERATION_EXTENDED_CORE")
]
bulk_plot[, cluster_score := factor(cluster_score, levels = cluster_score_cols)]
bulk_plot[, variable := factor(variable, levels = c("LAP3_log2_expr", "LAP3_STATE_UNION", "TAM_MYELOID_EXTENDED_CORE", "HYPOXIA_EXTENDED_CORE", "PROLIFERATION_EXTENDED_CORE"))]
p_bulk <- ggplot(bulk_plot, aes(variable, cluster_score, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.2) +
  geom_text(aes(label = ifelse(is.finite(rho), sprintf("%.2f", rho), "")), size = 2.4) +
  facet_wrap(~ dataset, nrow = 1) +
  scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C", midpoint = 0, limits = c(-1, 1), na.value = "grey90") +
  labs(x = NULL, y = NULL, fill = "rho") +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.grid = element_blank())
save_plot(p_bulk, "bulk_malignant_cluster_score_correlation_heatmap", width = 10.5, height = 4.8)

gbmap_plot <- gbmap_assoc[entry_variant == "main_neoplastic" & target == "lap3_mean"]
gbmap_plot[, cluster_score := factor(cluster_score, levels = cluster_score_cols)]
p_gbmap <- ggplot(gbmap_plot, aes(author_state, cluster_score, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.2) +
  geom_text(aes(label = ifelse(is.finite(rho), sprintf("%.2f", rho), "")), size = 2.8) +
  scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C", midpoint = 0, limits = c(-1, 1), na.value = "grey90") +
  labs(x = "Core GBmap malignant state", y = NULL, fill = "rho") +
  theme_bw(base_size = 8) +
  theme(panel.grid = element_blank())
save_plot(p_gbmap, "gbmap_malignant_cluster_lap3_association_heatmap", width = 6.2, height = 3.8)

if (nrow(gbmspace_summary)) {
  sp_plot <- gbmspace_summary[target == "LAP3_STATE_UNION"]
  sp_plot[, cluster_score := factor(cluster_score, levels = c("LAP3_MALIGNANT_STATE_MODULE", cluster_score_cols))]
  p_sp <- ggplot(sp_plot, aes(cluster_score, median_depth_adjusted_rho)) +
    geom_hline(yintercept = 0, color = "grey45", linewidth = 0.25) +
    geom_col(fill = "#5B8DB8", width = 0.68) +
    coord_flip() +
    labs(x = NULL, y = "Median section depth-adjusted rho with LAP3_STATE_UNION") +
    theme_bw(base_size = 8) +
    theme(panel.grid.minor = element_blank())
  save_plot(p_sp, "gbmspace_malignant_cluster_spatial_summary", width = 6.2, height = 3.8)
}

coverage_table <- rbindlist(lapply(names(cluster_sets), function(cl) {
  genes <- cluster_sets[[cl]]
  rbindlist(lapply(names(bulk_mats), function(ds) {
    data.table(dataset = ds, cluster_score = paste0("LAP3_MALIGNANT_", cl), n_genes = length(genes), n_present = sum(genes %in% rownames(bulk_mats[[ds]])))
  }))
}))
coverage_table[, coverage := n_present / n_genes]
if (nrow(gbmspace_coverage)) {
  sp_cov <- gbmspace_coverage[, .(
    dataset = "GBM_Space",
    n_sections = .N,
    min_coverage = min(present_fraction, na.rm = TRUE),
    median_coverage = median(present_fraction, na.rm = TRUE)
  ), by = cluster_score]
} else {
  sp_cov <- data.table()
}
write_table(coverage_table, "bulk_malignant_cluster_gene_coverage.csv")
if (nrow(sp_cov)) write_table(sp_cov, "gbmspace_malignant_cluster_gene_coverage_summary.tsv")

primary_bulk <- bulk_cor[
  stratum == "all" & variable %in% c("LAP3_log2_expr", "LAP3_STATE_UNION"),
  .(dataset, cluster_score, variable, n, rho = round(rho, 3), p_value = signif(p_value, 3))
]
primary_gbmap <- gbmap_assoc[
  entry_variant == "main_neoplastic" & target == "lap3_mean",
  .(author_state, cluster_score, n, rho = round(rho, 3), p_value = signif(p_value, 3))
]
primary_spatial <- if (nrow(gbmspace_summary)) {
  gbmspace_summary[target == "LAP3_STATE_UNION", .(
    cluster_score,
    n_sections,
    n_tumors,
    median_depth_adjusted_rho = round(median_depth_adjusted_rho, 3),
    p_depth_adjusted = signif(p_depth_adjusted, 3),
    p_adj_BH = signif(p_adj_BH, 3)
  )]
} else {
  data.table()
}

verdict <- data.table(
  item = c(
    "malignant_module_gene_count",
    "chosen_cluster_count",
    "cluster_sizes",
    "bulk_primary_pattern",
    "gbmap_primary_pattern",
    "gbmspace_primary_pattern",
    "interpretation"
  ),
  value = c(
    as.character(length(malignant_genes)),
    as.character(chosen_k),
    paste(names(cluster_sets), lengths(cluster_sets), sep = "=", collapse = "; "),
    "Cluster scores are tested against LAP3 and LAP3_STATE_UNION in TCGA/CGGA with GBM/LGG strata.",
    "Cluster scores are tested against LAP3 donor-state mean within AC/MES/NPC/OPC in Core GBmap.",
    if (nrow(primary_spatial)) "Cluster scores are tested against LAP3_STATE_UNION in GBM-Space using section-level depth-adjusted correlations." else "Spatial scoring skipped or not available.",
    "The 144-gene malignant-state module is decomposed into reproducible bulk co-expression clusters and projected across bulk, Core GBmap and GBM-Space. Positive clusters support interpretable malignant-state structure, not causal LAP3 biology."
  )
)
write_table(verdict, "malignant_state_module_audit_verdict.csv")

saveRDS(
  list(
    generated_at = format(Sys.time()),
    malignant_genes = malignant_genes,
    cluster_sets = cluster_sets,
    cluster_gene_table = cluster_gene_table,
    silhouette_table = silhouette_table,
    bulk_cor = bulk_cor,
    gbmap_assoc = gbmap_assoc,
    gbmspace_summary = gbmspace_summary,
    verdict = verdict
  ),
  file.path(export_dir, "lap3_malignant_state_module_audit_results.rds")
)

readme <- c(
  "# LAP3 Malignant-State Module Audit",
  "",
  paste0("Generated: ", format(Sys.time())),
  "",
  "## Purpose",
  "",
  "Audit the 144-gene `LAP3_MALIGNANT_STATE_MODULE` so it is not treated as a black-box composite score.",
  "",
  "## Inputs",
  "",
  "- `Data_Bulk_TCGA/Data_Merged/results/LAP3_State_Submodules/exports/lap3_state_submodules_projection.rds`",
  "- TCGA TPM, CGGA693 TPM, CGGA325 RSEM expression matrices",
  "- `Data_scRNA_GEO/GBmap_Core/cache/core_gbmap_lap3_state_union_lightweight.rds`",
  "- GBM-Space H5AD lightweight cache manifest when spatial scoring is not skipped",
  "",
  "## Methods",
  "",
  paste0("- Clustered 144 malignant-state genes using combined TCGA/CGGA gene-wise z-scored bulk expression; evaluated k=2:4, then merged sub-5-gene outlier clusters into the nearest larger co-expression cluster, yielding ", chosen_k, " scoreable clusters."),
  "- Scored each malignant gene cluster in TCGA, CGGA693 and CGGA325.",
  "- Projected cluster scores to Core GBmap malignant cells using the existing lightweight sparse cache and summarized at donor-state level.",
  "- Projected cluster scores to GBM-Space spots from H5AD when enabled, then summarized section-level depth-adjusted correlations.",
  "",
  "## Key Outputs",
  "",
  "- `tables/malignant_module_gene_clusters.csv`",
  "- `tables/bulk_malignant_cluster_score_correlations.csv`",
  "- `tables/gbmap_core_malignant_cluster_donor_state_associations.csv`",
  "- `tables/gbmspace_malignant_cluster_spatial_summary.tsv` if spatial scoring was enabled",
  "- `source_data/bulk_malignant_cluster_scores.csv`",
  "- `source_data/gbmap_core_malignant_cluster_cell_scores.csv.gz`",
  "- `source_data/gbmspace_spot_malignant_cluster_scores.tsv.gz` if spatial scoring was enabled",
  "- `plots/*.{svg,pdf,tiff,png}`",
  "- `exports/lap3_malignant_state_module_audit_results.rds`",
  "",
  "## Verdict",
  "",
  paste(capture.output(print(verdict)), collapse = "\n"),
  "",
  "## Bulk Primary Pattern",
  "",
  paste(capture.output(print(primary_bulk)), collapse = "\n"),
  "",
  "## Core GBmap Primary Pattern",
  "",
  paste(capture.output(print(primary_gbmap)), collapse = "\n"),
  "",
  "## GBM-Space Primary Pattern",
  "",
  if (nrow(primary_spatial)) paste(capture.output(print(primary_spatial)), collapse = "\n") else "Spatial scoring skipped.",
  "",
  "## Interpretation Boundary",
  "",
  "This audit decomposes the malignant-state component for interpretability and figure defense. It does not prove a LAP3-driven malignant-cell-intrinsic mechanism, enzymatic activity, leucine flux, or causal mTORC1 activation."
)
writeLines(readme, file.path(out_dir, "README.md"))

message_ts("Completed malignant-state module audit")
