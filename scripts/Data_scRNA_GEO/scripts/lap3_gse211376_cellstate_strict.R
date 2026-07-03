#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(Matrix)
})

setwd("/home/lzb/glioma")
data.table::setDTthreads(8)
set.seed(20260629)
bootstrap_replicates_main <- 500L

source("Data_scRNA_GEO/scripts/helpers/scRNA_inference_helpers.R")

matrix_file <- "Data_scRNA_GEO/GSE211376/GSE211376_raw_counts_Ruiz2022_all_samples_filtered_cells.tsv.gz"
metadata_file <- "Data_scRNA_GEO/GSE211376/GSE211376_metadata_Ruiz2022_all_samples_filtered_cells.csv.gz"
gene_set_file <- "Data_scRNA_GEO/results/LAP3_CellState_Phase0/source_data/frozen_cellstate_gene_sets.rds"

out_dir <- "Data_scRNA_GEO/results/GSE211376_LAP3_CellState"
table_dir <- file.path(out_dir, "tables")
log_dir <- file.path(out_dir, "logs")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(log_dir, "lap3_gse211376_cellstate_strict.log")
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message")
on.exit({
  sink(type = "message")
  sink(type = "output")
  close(log_con)
}, add = TRUE)

write_table <- function(x, filename) {
  data.table::fwrite(x, file.path(table_dir, filename))
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

residual_spearman_test <- function(data, x, y, patient, min_n = 6) {
  required <- c(x, y, patient)
  keep <- complete.cases(data[, required, drop = FALSE])
  analysis_data <- data[keep, , drop = FALSE]
  if (
    nrow(analysis_data) < min_n ||
    n_distinct(analysis_data[[patient]]) < 3 ||
    length(unique(analysis_data[[x]])) < 3 ||
    length(unique(analysis_data[[y]])) < 3
  ) {
    return(list(n = nrow(analysis_data), rho = NA_real_, p_value = NA_real_))
  }
  fit_x <- lm(reformulate(patient, response = x), data = analysis_data)
  fit_y <- lm(reformulate(patient, response = y), data = analysis_data)
  test <- suppressWarnings(cor.test(
    residuals(fit_x),
    residuals(fit_y),
    method = "spearman",
    exact = FALSE
  ))
  list(n = nrow(analysis_data), rho = unname(test$estimate), p_value = test$p.value)
}

score_gene_sets <- function(log_cpm, gene_sets, columns) {
  scores <- matrix(
    NA_real_,
    nrow = length(columns),
    ncol = length(gene_sets),
    dimnames = list(colnames(log_cpm)[columns], names(gene_sets))
  )
  for (signature in names(gene_sets)) {
    genes <- intersect(gene_sets[[signature]], rownames(log_cpm))
    if (length(genes) < 5L) {
      next
    }
    x <- log_cpm[genes, columns, drop = FALSE]
    z <- t(scale(t(x)))
    z[!is.finite(z)] <- 0
    scores[, signature] <- colMeans(z)
  }
  scores
}

score_cell_rank_gene_sets <- function(log_cpm, gene_sets) {
  gene_universe <- rownames(log_cpm)
  rank_percentiles <- apply(log_cpm, 2, function(x) {
    rank(x, ties.method = "average", na.last = "keep") / length(x)
  })
  rownames(rank_percentiles) <- gene_universe
  scores <- matrix(
    NA_real_,
    nrow = ncol(log_cpm),
    ncol = length(gene_sets),
    dimnames = list(colnames(log_cpm), names(gene_sets))
  )
  for (signature in names(gene_sets)) {
    genes <- intersect(gene_sets[[signature]], gene_universe)
    if (length(genes) < 5L) {
      next
    }
    scores[, signature] <- colMeans(rank_percentiles[genes, , drop = FALSE], na.rm = TRUE)
  }
  scores
}

run_pathway_within_state <- function(
    data,
    pathway_names,
    threshold_label,
    value_suffix = "",
    bootstrap_replicates = bootstrap_replicates_main) {
  bind_rows(lapply(sort(unique(data$author_state)), function(state) {
    d_state <- data %>% filter(author_state == state)
    bind_rows(lapply(pathway_names, function(pathway) {
      pathway_column <- paste0(pathway, value_suffix)
      res <- spearman_safe(d_state$lap3_log1p_cpm, d_state[[pathway_column]], min_n = 6)
      stat_fun <- function(d) {
        suppressWarnings(cor(
          d$lap3_log1p_cpm,
          d[[pathway_column]],
          method = "spearman",
          use = "complete.obs"
        ))
      }
      boot <- if (bootstrap_replicates > 0L) {
        cluster_bootstrap(d_state, cluster = "patient", statistic = stat_fun, replicates = bootstrap_replicates)
      } else {
        list(ci_low = NA_real_, ci_high = NA_real_)
      }
      lopo <- leave_one_cluster_out(d_state, cluster = "patient", statistic = stat_fun)
      data.frame(
        estimand = "within_author_state_across_patients",
        threshold = threshold_label,
        score_type = ifelse(nzchar(value_suffix), sub("^_", "", value_suffix), "pseudobulk_logcpm_zmean"),
        author_state = state,
        pathway = pathway,
        fdr_family = ifelse(pathway %in% c("HALLMARK_MTORC1_SIGNALING", "LEUCINE_BCAA_CORE"), "primary", "secondary"),
        n_patients = res$n,
        spearman_rho = res$rho,
        ci_low = boot$ci_low,
        ci_high = boot$ci_high,
        p_value = res$p_value,
        lopo_min_rho = suppressWarnings(min(lopo$estimate, na.rm = TRUE)),
        lopo_max_rho = suppressWarnings(max(lopo$estimate, na.rm = TRUE)),
        stringsAsFactors = FALSE
      )
    }))
  })) %>%
    adjust_fdr_by_family(p_column = "p_value", family_column = "fdr_family")
}

run_lopo_detail <- function(data, pathway_names, threshold_label, value_suffix = "") {
  bind_rows(lapply(sort(unique(data$author_state)), function(state) {
    d_state <- data %>% filter(author_state == state)
    bind_rows(lapply(pathway_names, function(pathway) {
      pathway_column <- paste0(pathway, value_suffix)
      stat_fun <- function(d) {
        suppressWarnings(cor(
          d$lap3_log1p_cpm,
          d[[pathway_column]],
          method = "spearman",
          use = "complete.obs"
        ))
      }
      leave_one_cluster_out(d_state, cluster = "patient", statistic = stat_fun) %>%
        transmute(
          threshold = threshold_label,
          score_type = ifelse(nzchar(value_suffix), sub("^_", "", value_suffix), "pseudobulk_logcpm_zmean"),
          author_state = state,
          pathway = pathway,
          omitted_patient = omitted_cluster,
          lopo_rho = estimate
        )
    }))
  }))
}

run_depth_adjusted_pathway <- function(data, pathway_names, threshold_label, value_suffix = "") {
  bind_rows(lapply(sort(unique(data$author_state)), function(state) {
    d_state <- data %>% filter(author_state == state)
    bind_rows(lapply(pathway_names, function(pathway) {
      pathway_column <- paste0(pathway, value_suffix)
      keep <- is.finite(d_state$lap3_log1p_cpm) &
        is.finite(d_state[[pathway_column]]) &
        is.finite(d_state$library_size) &
        is.finite(d_state$n_cells)
      analysis_data <- d_state[keep, , drop = FALSE] %>%
        mutate(log10_library_size = log10(library_size + 1))
      if (
        nrow(analysis_data) < 6L ||
        length(unique(analysis_data$lap3_log1p_cpm)) < 3L ||
        length(unique(analysis_data[[pathway_column]])) < 3L
      ) {
        return(data.frame(
          threshold = threshold_label,
          score_type = ifelse(nzchar(value_suffix), sub("^_", "", value_suffix), "pseudobulk_logcpm_zmean"),
          author_state = state,
          pathway = pathway,
          fdr_family = ifelse(pathway %in% c("HALLMARK_MTORC1_SIGNALING", "LEUCINE_BCAA_CORE"), "primary", "secondary"),
          n_patients = nrow(analysis_data),
          depth_adjusted_spearman_rho = NA_real_,
          p_value = NA_real_,
          stringsAsFactors = FALSE
        ))
      }
      lap3_resid <- residuals(lm(lap3_log1p_cpm ~ log10_library_size + n_cells, data = analysis_data))
      pathway_resid <- residuals(lm(reformulate(c("log10_library_size", "n_cells"), response = pathway_column), data = analysis_data))
      test <- suppressWarnings(cor.test(lap3_resid, pathway_resid, method = "spearman", exact = FALSE))
      data.frame(
        threshold = threshold_label,
        score_type = ifelse(nzchar(value_suffix), sub("^_", "", value_suffix), "pseudobulk_logcpm_zmean"),
        author_state = state,
        pathway = pathway,
        fdr_family = ifelse(pathway %in% c("HALLMARK_MTORC1_SIGNALING", "LEUCINE_BCAA_CORE"), "primary", "secondary"),
        n_patients = nrow(analysis_data),
        depth_adjusted_spearman_rho = unname(test$estimate),
        p_value = test$p.value,
        stringsAsFactors = FALSE
      )
    }))
  })) %>%
    adjust_fdr_by_family(p_column = "p_value", family_column = "fdr_family")
}

cat("Started:", format(Sys.time()), "\n")
cat("R:", R.version.string, "\n")
cat("data.table threads:", data.table::getDTthreads(), "\n")

gene_sets <- readRDS(gene_set_file)
neftel_sets <- gene_sets$neftel[c("AC", "OPC", "NPC1", "NPC2", "MES1", "MES2")]
pathway_sets <- gene_sets$pathways
all_gene_sets <- c(neftel_sets, pathway_sets)

metadata <- read.csv(gzfile(metadata_file), stringsAsFactors = FALSE, check.names = FALSE)
metadata$cell_id <- rownames(metadata)
rownames(metadata) <- NULL
names(metadata)[names(metadata) == "predicted.high_hierarchy"] <- "cell_state_raw"

malignant_map <- c(
  "AC-like" = "AC",
  "MES-like" = "MES",
  "NPC-like" = "NPC",
  "OPC-like" = "OPC"
)
tam_states <- c("TAM-BDM", "TAM-MG")
metadata <- metadata %>%
  mutate(
    author_state = unname(malignant_map[cell_state_raw]),
    compartment = case_when(
      !is.na(author_state) ~ "Malignant",
      cell_state_raw %in% tam_states ~ "TAM",
      TRUE ~ "Other"
    )
  )

stopifnot(
  nrow(metadata) == 39355,
  n_distinct(metadata$patient) == 11,
  all(names(malignant_map) %in% metadata$cell_state_raw)
)

cell_ids <- strsplit(readLines(gzfile(matrix_file), n = 1), "\t", fixed = TRUE)[[1]]
stopifnot(length(cell_ids) == nrow(metadata), identical(cell_ids, metadata$cell_id))

target_genes <- unique(c("LAP3", unlist(all_gene_sets, use.names = FALSE)))
cat("Target genes requested:", length(target_genes), "\n")
cat("Loading count matrix...\n")
counts_dt <- data.table::fread(
  matrix_file,
  skip = 1,
  header = FALSE,
  sep = "\t",
  data.table = TRUE,
  showProgress = TRUE,
  nThread = 8
)
setnames(counts_dt, 1, "gene")
stopifnot(nrow(counts_dt) == 27102, ncol(counts_dt) == length(cell_ids) + 1)
genes_present <- counts_dt$gene

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
write_table(coverage, "gse211376_cellstate_gene_coverage.csv")
stopifnot("LAP3" %in% genes_present)

library_sizes <- vapply(counts_dt[, -1], sum, numeric(1))
selected_dt <- counts_dt[gene %in% target_genes]
selected_genes <- selected_dt$gene
selected_mat <- as.matrix(selected_dt[, -1])
storage.mode(selected_mat) <- "numeric"
rownames(selected_mat) <- selected_genes
colnames(selected_mat) <- cell_ids
rm(counts_dt)
invisible(gc())

metadata$library_size <- library_sizes
metadata$lap3_count <- selected_mat["LAP3", ]
metadata$lap3_detected <- metadata$lap3_count > 0
metadata$lap3_log1p_cpm <- log1p(1e6 * metadata$lap3_count / metadata$library_size)

cat("Computing cell-level rank-based pathway scores...\n")
selected_log_cpm <- log1p(sweep(selected_mat, 2, metadata$library_size, "/") * 1e6)
cell_rank_scores <- score_cell_rank_gene_sets(selected_log_cpm, pathway_sets)
cell_rank_score_df <- data.frame(
  cell_id = rownames(cell_rank_scores),
  cell_rank_scores,
  check.names = FALSE
)
names(cell_rank_score_df)[-1] <- paste0(names(cell_rank_score_df)[-1], "_rank_score")
metadata <- metadata %>%
  left_join(cell_rank_score_df, by = "cell_id")
rm(selected_log_cpm, cell_rank_scores, cell_rank_score_df)
invisible(gc())

cell_detection <- metadata %>%
  filter(compartment %in% c("Malignant", "TAM")) %>%
  group_by(compartment, cell_state_raw, author_state, patient) %>%
  summarise(
    n_cells = n(),
    library_size = sum(library_size),
    lap3_positive_cells = sum(lap3_detected),
    lap3_detection_rate = mean(lap3_detected),
    mean_lap3_log1p_cpm = mean(lap3_log1p_cpm),
    median_lap3_log1p_cpm = median(lap3_log1p_cpm),
    .groups = "drop"
  )
write_table(cell_detection, "gse211376_lap3_detection_by_patient_state.csv")

state_detection <- cell_detection %>%
  group_by(compartment, cell_state_raw, author_state) %>%
  summarise(
    n_patients = n_distinct(patient),
    n_cells = sum(n_cells),
    lap3_positive_cells = sum(lap3_positive_cells),
    overall_lap3_detection_rate = lap3_positive_cells / n_cells,
    median_patient_detection_rate = median(lap3_detection_rate),
    median_patient_lap3_log1p_cpm = median(mean_lap3_log1p_cpm),
    .groups = "drop"
  )
write_table(state_detection, "gse211376_lap3_detection_state_summary.csv")

group_index <- metadata %>%
  filter(compartment == "Malignant") %>%
  mutate(group_id = paste(patient, author_state, sep = "||")) %>%
  select(cell_id, patient, cell_state_raw, author_state, group_id)
group_factor <- factor(group_index$group_id, levels = unique(group_index$group_id))
group_cells <- match(group_index$cell_id, cell_ids)
design <- Matrix::sparse.model.matrix(~ 0 + group_factor)
colnames(design) <- sub("^group_factor", "", colnames(design))

cat("Aggregating malignant patient-state pseudobulk...\n")
pseudobulk_counts <- selected_mat[, group_cells, drop = FALSE] %*% design
pseudobulk_counts <- as.matrix(pseudobulk_counts)
pseudobulk_lib <- as.numeric(rowsum(metadata$library_size[group_cells], group_factor, reorder = FALSE))
group_meta <- group_index %>%
  distinct(group_id, patient, cell_state_raw, author_state) %>%
  mutate(
    n_cells = as.integer(table(group_factor)[group_id]),
    library_size = pseudobulk_lib[match(group_id, levels(group_factor))]
  )
stopifnot(identical(colnames(pseudobulk_counts), group_meta$group_id))
log_cpm <- log1p(sweep(pseudobulk_counts, 2, group_meta$library_size, "/") * 1e6)

eligible <- group_meta$n_cells >= 20
neftel_scores <- score_gene_sets(log_cpm, neftel_sets, which(eligible))
pathway_scores <- score_gene_sets(log_cpm, pathway_sets, which(eligible))
score_df <- data.frame(
  group_id = rownames(neftel_scores),
  neftel_scores,
  pathway_scores,
  check.names = FALSE
)

rank_score_columns <- paste0(names(pathway_sets), "_rank_score")
patient_state_rank_scores <- metadata %>%
  filter(compartment == "Malignant") %>%
  mutate(group_id = paste(patient, author_state, sep = "||")) %>%
  group_by(group_id) %>%
  summarise(
    across(all_of(rank_score_columns), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

pseudobulk <- group_meta %>%
  mutate(
    eligible_n20 = n_cells >= 20,
    lap3_count = as.numeric(pseudobulk_counts["LAP3", ]),
    lap3_log1p_cpm = as.numeric(log_cpm["LAP3", ])
  ) %>%
  left_join(score_df, by = "group_id") %>%
  left_join(patient_state_rank_scores, by = "group_id") %>%
  mutate(
    neftel_mes_score = rowMeans(across(c(MES1, MES2)), na.rm = TRUE),
    neftel_npc_score = rowMeans(across(c(NPC1, NPC2)), na.rm = TRUE),
    dominant_continuous_state = case_when(
      !eligible_n20 ~ NA_character_,
      TRUE ~ c("AC", "OPC", "NPC", "MES")[max.col(cbind(AC, OPC, neftel_npc_score, neftel_mes_score), ties.method = "first")]
    ),
    dominant_margin = apply(
      cbind(AC, OPC, neftel_npc_score, neftel_mes_score),
      1,
      function(x) {
        if (any(!is.finite(x))) return(NA_real_)
        sx <- sort(x, decreasing = TRUE)
        sx[1] - sx[2]
      }
    )
  )
write_table(pseudobulk, "gse211376_patient_state_pseudobulk_cellstate_scores.csv")

rank_score_state_summary <- metadata %>%
  filter(compartment == "Malignant") %>%
  group_by(author_state, patient) %>%
  summarise(
    n_cells = n(),
    across(all_of(rank_score_columns), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = all_of(rank_score_columns),
    names_to = "pathway_rank_score",
    values_to = "mean_cell_rank_score"
  ) %>%
  group_by(author_state, pathway_rank_score) %>%
  summarise(
    n_patients = n_distinct(patient),
    n_cells = sum(n_cells),
    median_patient_mean_rank_score = median(mean_cell_rank_score, na.rm = TRUE),
    iqr_patient_mean_rank_score = IQR(mean_cell_rank_score, na.rm = TRUE),
    .groups = "drop"
  )
write_table(rank_score_state_summary, "gse211376_cell_rank_pathway_state_summary.csv")

analysis_pb <- pseudobulk %>% filter(eligible_n20)
state_pairs <- combn(sort(unique(analysis_pb$author_state)), 2, simplify = FALSE)
state_preference <- bind_rows(lapply(state_pairs, function(pair) {
  wide <- analysis_pb %>%
    filter(author_state %in% pair) %>%
    select(patient, author_state, lap3_log1p_cpm) %>%
    pivot_wider(names_from = author_state, values_from = lap3_log1p_cpm) %>%
    filter(is.finite(.data[[pair[1]]]), is.finite(.data[[pair[2]]]))
  wide$delta <- wide[[pair[2]]] - wide[[pair[1]]]
  stat_fun <- function(d) mean(d$delta, na.rm = TRUE)
  boot <- cluster_bootstrap(wide, cluster = "patient", statistic = stat_fun, replicates = bootstrap_replicates_main)
  test <- if (nrow(wide) >= 4 && length(unique(wide$delta)) >= 3) {
    suppressWarnings(wilcox.test(wide$delta, mu = 0, exact = FALSE))
  } else {
    list(p.value = NA_real_)
  }
  data.frame(
    contrast = paste(pair[2], "minus", pair[1]),
    state_low = pair[1],
    state_high = pair[2],
    n_patients = nrow(wide),
    mean_delta = mean(wide$delta),
    median_delta = median(wide$delta),
    ci_low = boot$ci_low,
    ci_high = boot$ci_high,
    p_value = test$p.value,
    stringsAsFactors = FALSE
  )
})) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH"))
write_table(state_preference, "gse211376_lap3_author_state_preference_paired.csv")

continuous_state_assoc <- bind_rows(lapply(c("AC", "OPC", "NPC1", "NPC2", "MES1", "MES2", "neftel_npc_score", "neftel_mes_score"), function(signature) {
  stat_fun <- function(d) {
    suppressWarnings(cor(d$lap3_log1p_cpm, d[[signature]], method = "spearman", use = "complete.obs"))
  }
  boot <- cluster_bootstrap(analysis_pb, cluster = "patient", statistic = stat_fun, replicates = bootstrap_replicates_main)
  res <- residual_spearman_test(analysis_pb, x = "lap3_log1p_cpm", y = signature, patient = "patient")
  data.frame(
    estimand = "continuous_state_score_patient_residual",
    signature = signature,
    n_patient_state = nrow(analysis_pb),
    n_patients = n_distinct(analysis_pb$patient),
    spearman_rho = res$rho,
    ci_low = boot$ci_low,
    ci_high = boot$ci_high,
    p_value = res$p_value,
    stringsAsFactors = FALSE
  )
})) %>%
  mutate(p_adj_BH_exploratory_state = p.adjust(p_value, method = "BH"))
write_table(continuous_state_assoc, "gse211376_lap3_continuous_state_score_associations.csv")

pathway_names <- names(pathway_sets)
pathway_assoc <- run_pathway_within_state(analysis_pb, pathway_names, threshold_label = "n_cells_ge_20")
write_table(pathway_assoc, "gse211376_lap3_pathway_within_state_associations.csv")

rank_pathway_assoc <- run_pathway_within_state(
  analysis_pb,
  pathway_names,
  threshold_label = "n_cells_ge_20",
  value_suffix = "_rank_score",
  bootstrap_replicates = 0L
)
write_table(rank_pathway_assoc, "gse211376_lap3_rank_pathway_within_state_associations.csv")

depth_adjusted_pathway <- run_depth_adjusted_pathway(
  analysis_pb,
  pathway_names,
  threshold_label = "n_cells_ge_20"
)
write_table(depth_adjusted_pathway, "gse211376_lap3_pathway_depth_adjusted_associations.csv")

lopo_detail <- run_lopo_detail(analysis_pb, pathway_names, threshold_label = "n_cells_ge_20")
write_table(lopo_detail, "gse211376_lap3_pathway_lopo_detail.csv")

threshold_retention <- bind_rows(lapply(c(20L, 50L, 100L), function(threshold) {
  pseudobulk %>%
    filter(n_cells >= threshold) %>%
    group_by(author_state) %>%
    summarise(
      threshold_n_cells = threshold,
      n_patient_state = n(),
      n_patients = n_distinct(patient),
      .groups = "drop"
    )
}))
write_table(threshold_retention, "gse211376_patient_state_threshold_retention.csv")

threshold_pathway_sensitivity <- bind_rows(lapply(c(20L, 50L, 100L), function(threshold) {
  d_threshold <- pseudobulk %>% filter(n_cells >= threshold)
  run_pathway_within_state(
    d_threshold,
    pathway_names,
    threshold_label = paste0("n_cells_ge_", threshold),
    bootstrap_replicates = 0L
  )
}))
write_table(threshold_pathway_sensitivity, "gse211376_lap3_pathway_threshold_sensitivity.csv")

cat("Eligible malignant patient-state rows n>=20:", nrow(analysis_pb), "\n")
cat("Eligible patients:", n_distinct(analysis_pb$patient), "\n")
cat("Completed:", format(Sys.time()), "\n")
