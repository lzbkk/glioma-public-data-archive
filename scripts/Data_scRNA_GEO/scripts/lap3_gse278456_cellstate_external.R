#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
  library(dplyr)
  library(tidyr)
})

setwd("/home/lzb/glioma")
data.table::setDTthreads(8)
set.seed(20260629)

source("Data_scRNA_GEO/scripts/helpers/scRNA_inference_helpers.R")

object_file <- "Data_scRNA_GEO/GSE278456/Tumor_Integrated_SeuratV5.rds"
crosswalk_file <- "Data_scRNA_GEO/results/GSE278456_Tumor_Object_Audit/exports/gse278456_local_author_metadata_crosswalk.rds"
gene_set_file <- "Data_scRNA_GEO/results/LAP3_CellState_Phase0/source_data/frozen_cellstate_gene_sets.rds"
threshold_file <- "Data_scRNA_GEO/results/GSE131928_LAP3_CellState/tables/gse131928_frozen_state_thresholds.csv"

out_dir <- "Data_scRNA_GEO/results/GSE278456_LAP3_CellState"
table_dir <- file.path(out_dir, "tables")
source_dir <- file.path(out_dir, "source_data")
log_dir <- file.path(out_dir, "logs")
plot_dir <- file.path(out_dir, "plots")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(log_dir, "lap3_gse278456_cellstate_external.log")
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message")
on.exit({
  sink(type = "message")
  sink(type = "output")
  close(log_con)
}, add = TRUE)

write_table <- function(x, filename) {
  fwrite(x, file.path(table_dir, filename))
}

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

mean_score <- function(expr, genes) {
  genes <- intersect(genes, rownames(expr))
  if (length(genes) < 2L) {
    return(rep(NA_real_, ncol(expr)))
  }
  x <- expr[genes, , drop = FALSE]
  z <- t(scale(t(as.matrix(x))))
  z[!is.finite(z)] <- 0
  colMeans(z)
}

spearman_safe <- function(x, y, min_n = 6) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) < min_n || length(unique(x)) < 3 || length(unique(y)) < 3) {
    return(list(n = length(x), rho = NA_real_, p_value = NA_real_))
  }
  test <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  list(n = length(x), rho = unname(test$estimate), p_value = test$p.value)
}

finite_range <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(c(min = NA_real_, max = NA_real_))
  }
  c(min = min(x), max = max(x))
}

patient_consensus <- function(data, patient_col = "Patient_ID") {
  data %>%
    group_by(.data[[patient_col]]) %>%
    summarise(
      across(
        c(pathology, tumor_grade, idh_status, analysis_class),
        ~ {
          values <- unique(na.omit(as.character(.x)))
          if (!length(values)) NA_character_ else paste(values, collapse = ";")
        }
      ),
      n_cells = n(),
      .groups = "drop"
    ) %>%
    rename(patient = all_of(patient_col))
}

add_state_calls <- function(df, margin_main, margin_strict) {
  state_matrix <- as.matrix(df[, c("state_AC", "state_OPC", "state_NPC", "state_MES")])
  colnames(state_matrix) <- c("AC", "OPC", "NPC", "MES")
  top_index <- max.col(state_matrix, ties.method = "first")
  sorted <- t(apply(state_matrix, 1, sort, decreasing = TRUE))
  df$dominant_state <- colnames(state_matrix)[top_index]
  df$dominant_score <- sorted[, 1]
  df$second_score <- sorted[, 2]
  df$dominant_margin <- df$dominant_score - df$second_score
  df$high_conf_gse131928_10x_main <- df$dominant_score > 0 & df$dominant_margin >= margin_main
  df$high_conf_gse131928_10x_strict <- df$dominant_score > 0 & df$dominant_margin >= margin_strict
  df
}

cat("Started:", format(Sys.time()), "\n")
cat("Loading local GSE278456 tumor object...\n")
obj <- readRDS(object_file)
crosswalk <- readRDS(crosswalk_file)
stopifnot(ncol(obj) == nrow(crosswalk), identical(colnames(obj), crosswalk$cell_id_local))

crosswalk <- crosswalk %>%
  mutate(
    pathology = author_Pathology,
    tumor_grade = author_Tumor.Grade,
    idh_status = author_IDH.status,
    analysis_class = case_when(
      pathology == "GBM" & tumor_grade == "IV" & idh_status == "wild type" ~ "GBM_grade4_IDHwt",
      pathology %in% c("Astrocytoma", "Anaplastic Astrocytoma", "Oligodendroglioma") &
        tumor_grade %in% c("II", "III") & idh_status == "mutant" ~ "LGG_grade2_3_IDHmut",
      pathology == "Normal" | tumor_grade == "Normal" ~ "Normal",
      TRUE ~ "Other_tumor"
    )
  )

matched_cells <- which(crosswalk$matched_author)
stopifnot(length(matched_cells) == 193730)

gene_sets <- readRDS(gene_set_file)
neftel_sets <- gene_sets$neftel[c("AC", "OPC", "NPC1", "NPC2", "MES1", "MES2")]
pathway_sets <- gene_sets$pathways
all_gene_sets <- c(neftel_sets, pathway_sets)
qc_sets <- list(
  immune = c("PTPRC", "LST1", "TYROBP", "C1QA", "C1QB", "CD3D", "CD3E", "NKG7"),
  endothelial = c("PECAM1", "VWF", "KDR", "CLDN5", "FLT1"),
  oligodendrocyte = c("MBP", "PLP1", "MOBP", "MOG", "MAG"),
  malignant_glial = c("SOX2", "OLIG1", "OLIG2", "PDGFRA", "EGFR", "NES", "PROM1")
)
target_genes <- unique(c("LAP3", unlist(all_gene_sets, use.names = FALSE), unlist(qc_sets, use.names = FALSE)))

cat("Extracting normalized data layer for target genes...\n")
data_layer <- LayerData(obj[["RNA"]], layer = "data", fast = FALSE)
genes_present <- rownames(data_layer)
target_present <- intersect(target_genes, genes_present)
expr <- data_layer[target_present, matched_cells, drop = FALSE]
rm(data_layer)
invisible(gc())

coverage <- bind_rows(lapply(names(all_gene_sets), function(signature) {
  requested <- all_gene_sets[[signature]]
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
write_table(coverage, "gse278456_cellstate_gene_coverage.csv")
stopifnot("LAP3" %in% rownames(expr), all(coverage$genes_present >= 10))

cat("Scoring QC marker sets and Neftel states...\n")
qc_scores <- as.data.frame(lapply(qc_sets, function(genes) mean_score(expr, genes)))
names(qc_scores) <- paste0("qc_", names(qc_scores))
neftel_scores <- as.data.frame(score_gene_sets(expr, neftel_sets), check.names = FALSE)
pathway_scores <- as.data.frame(score_gene_sets(expr, pathway_sets), check.names = FALSE)

thresholds <- fread(threshold_file)
tenx_main_margin <- thresholds[platform == "10X" & threshold_set == "main_median_margin", margin_threshold][[1]]
tenx_strict_margin <- thresholds[platform == "10X" & threshold_set == "strict_q75_margin", margin_threshold][[1]]

cell_scores <- data.frame(
  cell_id = colnames(expr),
  patient = crosswalk$Patient_ID[matched_cells],
  pathology = crosswalk$pathology[matched_cells],
  tumor_grade = crosswalk$tumor_grade[matched_cells],
  idh_status = crosswalk$idh_status[matched_cells],
  analysis_class = crosswalk$analysis_class[matched_cells],
  seurat_cluster = as.character(obj$seurat_clusters[matched_cells]),
  nCount_RNA = obj$nCount_RNA[matched_cells],
  nFeature_RNA = obj$nFeature_RNA[matched_cells],
  percent_mt = obj$percent.mt[matched_cells],
  lap3_log_norm = as.numeric(expr["LAP3", ]),
  lap3_detected = as.numeric(expr["LAP3", ]) > 0,
  qc_scores,
  neftel_scores,
  pathway_scores,
  check.names = FALSE
) %>%
  mutate(
    state_AC = AC,
    state_OPC = OPC,
    state_NPC = rowMeans(across(c(NPC1, NPC2)), na.rm = TRUE),
    state_MES = rowMeans(across(c(MES1, MES2)), na.rm = TRUE),
    contamination_score = pmax(qc_immune, qc_endothelial, qc_oligodendrocyte, na.rm = TRUE),
    trusted_malignant = analysis_class != "Normal" &
      qc_malignant_glial > -0.25 &
      contamination_score < 0.75
  ) %>%
  add_state_calls(margin_main = tenx_main_margin, margin_strict = tenx_strict_margin)

cell_qc_summary <- cell_scores %>%
  group_by(analysis_class) %>%
  summarise(
    n_cells = n(),
    n_patients = n_distinct(patient),
    trusted_malignant_cells = sum(trusted_malignant),
    trusted_fraction = trusted_malignant_cells / n_cells,
    median_immune_score = median(qc_immune),
    median_endothelial_score = median(qc_endothelial),
    median_oligodendrocyte_score = median(qc_oligodendrocyte),
    median_malignant_glial_score = median(qc_malignant_glial),
    .groups = "drop"
  )
write_table(cell_qc_summary, "gse278456_malignant_qc_summary.csv")

state_distribution <- cell_scores %>%
  filter(trusted_malignant) %>%
  group_by(analysis_class, patient, dominant_state) %>%
  summarise(
    n_cells = n(),
    high_conf_main_cells = sum(high_conf_gse131928_10x_main),
    high_conf_strict_cells = sum(high_conf_gse131928_10x_strict),
    mean_lap3_log_norm = mean(lap3_log_norm),
    lap3_detection_rate = mean(lap3_detected),
    mean_state_AC = mean(state_AC),
    mean_state_OPC = mean(state_OPC),
    mean_state_NPC = mean(state_NPC),
    mean_state_MES = mean(state_MES),
    .groups = "drop"
  )
write_table(state_distribution, "gse278456_patient_state_distribution_all_dominant.csv")

highconf_patient_state <- cell_scores %>%
  filter(trusted_malignant, high_conf_gse131928_10x_main) %>%
  group_by(analysis_class, patient, dominant_state) %>%
  summarise(
    n_cells = n(),
    lap3_detection_rate = mean(lap3_detected),
    mean_lap3_log_norm = mean(lap3_log_norm),
    across(all_of(names(pathway_sets)), ~ mean(.x, na.rm = TRUE), .names = "{.col}"),
    .groups = "drop"
  ) %>%
  mutate(eligible_n20 = n_cells >= 20)
write_table(highconf_patient_state, "gse278456_highconf_patient_state_scores.csv")

continuous_patient <- cell_scores %>%
  filter(trusted_malignant) %>%
  group_by(analysis_class, patient) %>%
  summarise(
    n_cells = n(),
    lap3_detection_rate = mean(lap3_detected),
    mean_lap3_log_norm = mean(lap3_log_norm),
    mean_AC = mean(state_AC),
    mean_OPC = mean(state_OPC),
    mean_NPC = mean(state_NPC),
    mean_MES = mean(state_MES),
    mean_margin = mean(dominant_margin),
    high_conf_main_fraction = mean(high_conf_gse131928_10x_main),
    .groups = "drop"
  ) %>%
  mutate(
    eligible_n20 = n_cells >= 20
  )
write_table(continuous_patient, "gse278456_patient_continuous_state_scores.csv")

state_migration <- cell_scores %>%
  filter(trusted_malignant) %>%
  group_by(analysis_class) %>%
  summarise(
    n_cells = n(),
    n_patients = n_distinct(patient),
    median_dominant_score = median(dominant_score),
    q25_dominant_score = quantile(dominant_score, 0.25),
    q75_dominant_score = quantile(dominant_score, 0.75),
    median_margin = median(dominant_margin),
    q25_margin = quantile(dominant_margin, 0.25),
    q75_margin = quantile(dominant_margin, 0.75),
    high_conf_main_fraction = mean(high_conf_gse131928_10x_main),
    high_conf_strict_fraction = mean(high_conf_gse131928_10x_strict),
    .groups = "drop"
  ) %>%
  mutate(
    gse131928_10x_main_margin = tenx_main_margin,
    gse131928_10x_strict_margin = tenx_strict_margin
  )
write_table(state_migration, "gse278456_state_threshold_migration_check.csv")

continuous_assoc <- bind_rows(lapply(sort(unique(continuous_patient$analysis_class)), function(subset_name) {
  d <- continuous_patient %>% filter(analysis_class == subset_name, eligible_n20)
  bind_rows(lapply(c("mean_AC", "mean_OPC", "mean_NPC", "mean_MES"), function(signature) {
    res <- spearman_safe(d$mean_lap3_log_norm, d[[signature]], min_n = 6)
    stat_fun <- function(x) {
      suppressWarnings(cor(x$mean_lap3_log_norm, x[[signature]], method = "spearman", use = "complete.obs"))
    }
    boot <- if (is.finite(res$rho) && res$n >= 6) {
      cluster_bootstrap(d, cluster = "patient", statistic = stat_fun, replicates = 2000)
    } else {
      list(ci_low = NA_real_, ci_high = NA_real_)
    }
    data.frame(
      analysis_class = subset_name,
      estimand = "continuous_state_score_across_patients",
      signature = signature,
      n_patients = res$n,
      spearman_rho = res$rho,
      ci_low = boot$ci_low,
      ci_high = boot$ci_high,
      p_value = res$p_value,
      stringsAsFactors = FALSE
    )
  }))
})) %>%
  group_by(analysis_class) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
  ungroup()
write_table(continuous_assoc, "gse278456_lap3_continuous_state_associations.csv")

pathway_assoc <- highconf_patient_state %>%
  filter(eligible_n20) %>%
  group_by(analysis_class, dominant_state) %>%
  group_modify(function(.x, .y) {
    bind_rows(lapply(names(pathway_sets), function(pathway) {
      res <- spearman_safe(.x$mean_lap3_log_norm, .x[[pathway]], min_n = 6)
      stat_fun <- function(d) {
        suppressWarnings(cor(d$mean_lap3_log_norm, d[[pathway]], method = "spearman", use = "complete.obs"))
      }
      boot <- if (is.finite(res$rho) && res$n >= 6) {
        cluster_bootstrap(.x, cluster = "patient", statistic = stat_fun, replicates = 2000)
      } else {
        list(ci_low = NA_real_, ci_high = NA_real_)
      }
      lopo <- if (is.finite(res$rho) && res$n >= 6) {
        leave_one_cluster_out(.x, cluster = "patient", statistic = stat_fun)
      } else {
        data.frame(estimate = NA_real_)
      }
      lopo_range <- finite_range(lopo$estimate)
      data.frame(
        pathway = pathway,
        fdr_family = ifelse(pathway %in% c("HALLMARK_MTORC1_SIGNALING", "LEUCINE_BCAA_CORE"), "primary", "secondary"),
        n_patients = res$n,
        spearman_rho = res$rho,
        ci_low = boot$ci_low,
        ci_high = boot$ci_high,
        p_value = res$p_value,
        lopo_min_rho = lopo_range[["min"]],
        lopo_max_rho = lopo_range[["max"]],
        stringsAsFactors = FALSE
      )
    }))
  }) %>%
  ungroup() %>%
  adjust_fdr_by_family(p_column = "p_value", family_column = "fdr_family")
write_table(pathway_assoc, "gse278456_highconf_state_pathway_associations.csv")

cell_score_export <- cell_scores %>%
  select(
    cell_id, patient, pathology, tumor_grade, idh_status, analysis_class,
    seurat_cluster, nCount_RNA, nFeature_RNA, percent_mt,
    lap3_log_norm, lap3_detected, starts_with("qc_"),
    AC, OPC, NPC1, NPC2, MES1, MES2,
    state_AC, state_OPC, state_NPC, state_MES,
    dominant_state, dominant_score, second_score, dominant_margin,
    trusted_malignant, high_conf_gse131928_10x_main, high_conf_gse131928_10x_strict
  )
fwrite(cell_score_export, file.path(source_dir, "gse278456_cell_neftel_scores_qc.csv.gz"))
saveRDS(
  list(
    thresholds = thresholds,
    qc_summary = cell_qc_summary,
    migration = state_migration,
    continuous_patient = continuous_patient,
    highconf_patient_state = highconf_patient_state
  ),
  file.path(source_dir, "gse278456_cellstate_external_projection.rds")
)

cat("Matched cells:", length(matched_cells), "\n")
cat("Trusted malignant cells:", sum(cell_scores$trusted_malignant), "\n")
cat("Trusted malignant patients:", n_distinct(cell_scores$patient[cell_scores$trusted_malignant]), "\n")
cat("GSE131928 10X main margin:", tenx_main_margin, "\n")
cat("GSE131928 10X strict margin:", tenx_strict_margin, "\n")
cat("Completed:", format(Sys.time()), "\n")
