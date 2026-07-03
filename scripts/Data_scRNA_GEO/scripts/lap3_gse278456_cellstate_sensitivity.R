#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
})

setwd("/home/lzb/glioma")
data.table::setDTthreads(8)
set.seed(20260629)

source("Data_scRNA_GEO/scripts/helpers/scRNA_inference_helpers.R")

cache_file <- "Data_scRNA_GEO/results/GSE278456_LAP3_CellState/source_data/gse278456_cell_neftel_scores_qc.csv.gz"
pathway_cache_file <- "Data_scRNA_GEO/results/GSE278456_LAP3_CellState/source_data/gse278456_cell_pathway_scores.csv.gz"
projection_rds <- "Data_scRNA_GEO/results/GSE278456_LAP3_CellState/source_data/gse278456_cellstate_external_projection.rds"
out_dir <- "Data_scRNA_GEO/results/GSE278456_LAP3_CellState_Sensitivity"
table_dir <- file.path(out_dir, "tables")
source_dir <- file.path(out_dir, "source_data")
plot_dir <- file.path(out_dir, "plots")
log_dir <- file.path(out_dir, "logs")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(log_dir, "lap3_gse278456_cellstate_sensitivity.log")
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

make_rule <- function(data, rule_name) {
  non_normal <- data$analysis_class != "Normal"
  contamination <- pmax(
    data$qc_immune,
    data$qc_endothelial,
    data$qc_oligodendrocyte,
    na.rm = TRUE
  )
  if (rule_name == "relaxed") {
    return(non_normal & data$qc_malignant_glial > -0.50 & contamination < 1.00)
  }
  if (rule_name == "main") {
    return(non_normal & data$qc_malignant_glial > -0.25 & contamination < 0.75)
  }
  if (rule_name == "strict") {
    return(non_normal & data$qc_malignant_glial > 0.00 & contamination < 0.50)
  }
  stop("Unknown rule: ", rule_name)
}

summarise_continuous <- function(data, rule_name) {
  d <- data %>%
    filter(.data[[paste0("trusted_", rule_name)]]) %>%
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
      high_conf_strict_fraction = mean(high_conf_gse131928_10x_strict),
      .groups = "drop"
    ) %>%
    mutate(rule = rule_name, eligible_n20 = n_cells >= 20)
  d
}

continuous_assoc_by_rule <- function(patient_data) {
  bind_rows(lapply(sort(unique(patient_data$rule)), function(rule_name) {
    bind_rows(lapply(sort(unique(patient_data$analysis_class)), function(class_name) {
      d <- patient_data %>%
        filter(rule == rule_name, analysis_class == class_name, eligible_n20)
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
          rule = rule_name,
          analysis_class = class_name,
          signature = signature,
          n_patients = res$n,
          spearman_rho = res$rho,
          ci_low = boot$ci_low,
          ci_high = boot$ci_high,
          p_value = res$p_value,
          stringsAsFactors = FALSE
        )
      })) %>%
        mutate(p_adj_BH = p.adjust(p_value, method = "BH"))
    }))
  }))
}

highconf_pathway_from_patient_state <- function(patient_state, rule_name = "main", highconf_rule = "high_conf_gse131928_10x_main") {
  pathways <- c(
    "HALLMARK_MTORC1_SIGNALING",
    "LEUCINE_BCAA_CORE",
    "MTORC1_READOUT_CORE",
    "REACTOME_TRANSLATION"
  )

  patient_state %>%
    filter(eligible_n20) %>%
    group_by(analysis_class, dominant_state) %>%
    group_modify(function(.x, .y) {
      bind_rows(lapply(pathways, function(pathway) {
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
        rule = rule_name,
        highconf_rule = highconf_rule,
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
}

highconf_pathway_from_cells <- function(data, rule_name, highconf_col) {
  pathways <- c(
    "HALLMARK_MTORC1_SIGNALING",
    "LEUCINE_BCAA_CORE",
    "MTORC1_READOUT_CORE",
    "REACTOME_TRANSLATION"
  )

  d_patient_state <- data %>%
    filter(.data[[paste0("trusted_", rule_name)]], .data[[highconf_col]]) %>%
    group_by(analysis_class, patient, dominant_state) %>%
    summarise(
      n_cells = n(),
      mean_lap3_log_norm = mean(lap3_log_norm),
      across(all_of(pathways), ~ mean(.x, na.rm = TRUE), .names = "{.col}"),
      .groups = "drop"
    ) %>%
    mutate(eligible_n20 = n_cells >= 20)

  highconf_pathway_from_patient_state(
    d_patient_state,
    rule_name = rule_name,
    highconf_rule = highconf_col
  )
}

cat("Started:", format(Sys.time()), "\n")
cat("Reading cache:", cache_file, "\n")
cells <- fread(cache_file, nThread = 8)
projection <- readRDS(projection_rds)
stopifnot(
  all(c(
    "patient", "analysis_class", "seurat_cluster", "lap3_log_norm",
    "state_AC", "state_OPC", "state_NPC", "state_MES",
    "qc_immune", "qc_endothelial", "qc_oligodendrocyte", "qc_malignant_glial"
  ) %in% names(cells))
)

for (rule_name in c("relaxed", "main", "strict")) {
  cells[[paste0("trusted_", rule_name)]] <- make_rule(cells, rule_name)
}

retention <- bind_rows(lapply(c("relaxed", "main", "strict"), function(rule_name) {
  cells %>%
    group_by(analysis_class) %>%
    summarise(
      rule = rule_name,
      n_cells_total = n(),
      n_patients_total = n_distinct(patient),
      n_cells_retained = sum(.data[[paste0("trusted_", rule_name)]]),
      n_patients_retained = n_distinct(patient[.data[[paste0("trusted_", rule_name)]]]),
      retained_fraction = n_cells_retained / n_cells_total,
      .groups = "drop"
    )
}))
write_table(retention, "gse278456_malignant_definition_retention.csv")

cluster_summary <- cells %>%
  filter(trusted_main) %>%
  group_by(seurat_cluster) %>%
  summarise(
    n_cells = n(),
    n_patients = n_distinct(patient),
    dominant_analysis_class = names(sort(table(analysis_class), decreasing = TRUE))[1],
    lap3_detection_rate = mean(lap3_detected),
    mean_lap3_log_norm = mean(lap3_log_norm),
    mean_AC = mean(state_AC),
    mean_OPC = mean(state_OPC),
    mean_NPC = mean(state_NPC),
    mean_MES = mean(state_MES),
    dominant_state_by_mean = c("AC", "OPC", "NPC", "MES")[
      which.max(c(mean(state_AC), mean(state_OPC), mean(state_NPC), mean(state_MES)))
    ],
    high_conf_main_fraction = mean(high_conf_gse131928_10x_main),
    median_margin = median(dominant_margin),
    .groups = "drop"
  ) %>%
  arrange(desc(n_cells))
write_table(cluster_summary, "gse278456_cluster_state_summary_main_qc.csv")

cluster_class_state <- cells %>%
  filter(trusted_main) %>%
  group_by(analysis_class, seurat_cluster, dominant_state) %>%
  summarise(
    n_cells = n(),
    n_patients = n_distinct(patient),
    mean_lap3_log_norm = mean(lap3_log_norm),
    lap3_detection_rate = mean(lap3_detected),
    high_conf_main_fraction = mean(high_conf_gse131928_10x_main),
    .groups = "drop"
  )
write_table(cluster_class_state, "gse278456_cluster_class_state_crosswalk.csv")

patient_cluster_influence <- cells %>%
  filter(trusted_main) %>%
  group_by(analysis_class, patient, seurat_cluster) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  group_by(analysis_class, patient) %>%
  mutate(patient_fraction = n_cells / sum(n_cells)) %>%
  ungroup()
write_table(patient_cluster_influence, "gse278456_patient_cluster_composition_main_qc.csv")

patient_by_rule <- bind_rows(lapply(c("relaxed", "main", "strict"), function(rule_name) {
  summarise_continuous(cells, rule_name)
}))
write_table(patient_by_rule, "gse278456_patient_continuous_state_by_malignant_rule.csv")

continuous_sensitivity <- continuous_assoc_by_rule(patient_by_rule)
write_table(continuous_sensitivity, "gse278456_continuous_state_malignant_rule_sensitivity.csv")

if (file.exists(pathway_cache_file)) {
  cat("Reading pathway score cache:", pathway_cache_file, "\n")
  pathway_cache <- fread(pathway_cache_file, nThread = 8)
  stopifnot("cell_id" %in% names(pathway_cache), identical(cells$cell_id, pathway_cache$cell_id))
  cells <- bind_cols(cells, pathway_cache %>% select(-cell_id))
  hard_sensitivity <- bind_rows(lapply(c("relaxed", "main", "strict"), function(rule_name) {
    bind_rows(
      highconf_pathway_from_cells(cells, rule_name, "high_conf_gse131928_10x_main"),
      highconf_pathway_from_cells(cells, rule_name, "high_conf_gse131928_10x_strict")
    )
  }))
} else {
  cat("Pathway score cache not found; using main-rule patient-state projection only.\n")
  hard_sensitivity <- highconf_pathway_from_patient_state(
    projection$highconf_patient_state,
    rule_name = "main",
    highconf_rule = "high_conf_gse131928_10x_main"
  )
}
write_table(hard_sensitivity, "gse278456_highconf_pathway_malignant_rule_sensitivity.csv")

key_summary <- bind_rows(
  continuous_sensitivity %>%
    filter(
      analysis_class %in% c("GBM_grade4_IDHwt", "LGG_grade2_3_IDHmut"),
      signature %in% c("mean_AC", "mean_MES")
    ) %>%
    transmute(
      evidence_type = "continuous_state",
      rule,
      highconf_rule = NA_character_,
      analysis_class,
      state_or_signature = signature,
      pathway = NA_character_,
      n_patients,
      spearman_rho,
      ci_low,
      ci_high,
      p_value,
      p_adj_BH
    ),
  hard_sensitivity %>%
    filter(
      analysis_class == "GBM_grade4_IDHwt",
      dominant_state == "AC",
      pathway %in% c("HALLMARK_MTORC1_SIGNALING", "LEUCINE_BCAA_CORE")
    ) %>%
    transmute(
      evidence_type = "highconf_state_pathway",
      rule,
      highconf_rule,
      analysis_class,
      state_or_signature = dominant_state,
      pathway,
      n_patients,
      spearman_rho,
      ci_low,
      ci_high,
      p_value,
      p_adj_BH
    )
)
write_table(key_summary, "gse278456_key_sensitivity_summary.csv")

saveRDS(
  list(
    retention = retention,
    cluster_summary = cluster_summary,
    patient_by_rule = patient_by_rule,
    continuous_sensitivity = continuous_sensitivity,
    hard_sensitivity = hard_sensitivity,
    key_summary = key_summary
  ),
  file.path(source_dir, "gse278456_cellstate_sensitivity.rds")
)

cat("Cells:", nrow(cells), "\n")
cat("Patients:", n_distinct(cells$patient), "\n")
cat("Main trusted cells:", sum(cells$trusted_main), "\n")
cat("Completed:", format(Sys.time()), "\n")
