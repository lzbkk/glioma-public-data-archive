#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(readxl)
})

setwd("/home/lzb/glioma")
data.table::setDTthreads(16)
set.seed(20260629)

source("Data_scRNA_GEO/scripts/helpers/scRNA_inference_helpers.R")

bootstrap_replicates_main <- 500L
smartseq_file <- "Data_scRNA_GEO/GSE131928/GSE131928_RAW/GSM3828672_Smartseq2_GBM_IDHwt_processed_TPM.tsv.gz"
tenx_file <- "Data_scRNA_GEO/GSE131928/GSE131928_RAW/GSM3828673_10X_GBM_IDHwt_processed_TPM.tsv.gz"
metadata_file <- "Data_scRNA_GEO/GSE131928/GSE131928_single_cells_tumor_name_and_adult_or_peidatric.xlsx"
gene_set_file <- "Data_scRNA_GEO/results/LAP3_CellState_Phase0/source_data/frozen_cellstate_gene_sets.rds"

out_dir <- "Data_scRNA_GEO/results/GSE131928_LAP3_CellState"
table_dir <- file.path(out_dir, "tables")
source_dir <- file.path(out_dir, "source_data")
log_dir <- file.path(out_dir, "logs")
plot_dir <- file.path(out_dir, "plots")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(log_dir, "lap3_gse131928_cellstate_calibration.log")
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
    z <- t(scale(t(x)))
    z[!is.finite(z)] <- 0
    scores[, signature] <- colMeans(z)
  }
  scores
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

summarise_patient_state <- function(cell_scores, confidence_column, minimum_cells) {
  cell_scores %>%
    filter(.data[[confidence_column]]) %>%
    group_by(platform, patient, age_group, dominant_state) %>%
    summarise(
      n_cells = n(),
      lap3_detection_rate = mean(lap3_detected),
      mean_lap3_log1p_tpm = mean(lap3_log1p_tpm),
      median_lap3_log1p_tpm = median(lap3_log1p_tpm),
      mean_log1p_detected_genes = mean(log1p(n_genes_detected)),
      median_detected_genes = median(n_genes_detected),
      mean_total_tpm = mean(total_tpm),
      across(
        all_of(names(pathway_sets)),
        ~ mean(.x, na.rm = TRUE),
        .names = "{.col}"
      ),
      .groups = "drop"
    ) %>%
    filter(n_cells >= minimum_cells) %>%
    mutate(
      confidence_rule = confidence_column,
      minimum_cells = minimum_cells
    )
}

run_pathway_associations <- function(patient_state_data, adjusted_for_depth = FALSE) {
  patient_state_data %>%
    group_by(confidence_rule, minimum_cells, platform, dominant_state) %>%
    group_modify(function(.x, .y) {
      bind_rows(lapply(names(pathway_sets), function(pathway) {
        analysis_data <- .x %>%
          select(patient, mean_lap3_log1p_tpm, mean_log1p_detected_genes, all_of(pathway)) %>%
          filter(if_all(everything(), ~ !is.na(.x)))
        n_patients <- n_distinct(analysis_data$patient)
        if (adjusted_for_depth && n_patients >= 8L &&
            n_distinct(analysis_data$mean_log1p_detected_genes) >= 3L) {
          lap3_value <- residuals(lm(
            mean_lap3_log1p_tpm ~ mean_log1p_detected_genes,
            data = analysis_data
          ))
          pathway_value <- residuals(lm(
            reformulate("mean_log1p_detected_genes", response = pathway),
            data = analysis_data
          ))
        } else if (!adjusted_for_depth) {
          lap3_value <- analysis_data$mean_lap3_log1p_tpm
          pathway_value <- analysis_data[[pathway]]
        } else {
          lap3_value <- rep(NA_real_, nrow(analysis_data))
          pathway_value <- rep(NA_real_, nrow(analysis_data))
        }
        res <- spearman_safe(lap3_value, pathway_value, min_n = ifelse(adjusted_for_depth, 8, 6))
        data.frame(
          pathway = pathway,
          fdr_family = ifelse(
            pathway %in% c("HALLMARK_MTORC1_SIGNALING", "LEUCINE_BCAA_CORE"),
            "primary",
            "secondary"
          ),
          n_patients = n_patients,
          spearman_rho = res$rho,
          p_value = res$p_value,
          depth_adjusted = adjusted_for_depth,
          stringsAsFactors = FALSE
        )
      }))
    }) %>%
    ungroup() %>%
    adjust_fdr_by_family(p_column = "p_value", family_column = "fdr_family")
}

collapse_state_scores <- function(score_df) {
  score_df %>%
    mutate(
      state_AC = AC,
      state_OPC = OPC,
      state_NPC = rowMeans(across(c(NPC1, NPC2)), na.rm = TRUE),
      state_MES = rowMeans(across(c(MES1, MES2)), na.rm = TRUE)
    )
}

add_dominant_state <- function(df, margin_main, margin_strict) {
  state_matrix <- as.matrix(df[, c("state_AC", "state_OPC", "state_NPC", "state_MES")])
  colnames(state_matrix) <- c("AC", "OPC", "NPC", "MES")
  top_index <- max.col(state_matrix, ties.method = "first")
  sorted <- t(apply(state_matrix, 1, sort, decreasing = TRUE))
  df$dominant_state <- colnames(state_matrix)[top_index]
  df$dominant_score <- sorted[, 1]
  df$second_score <- sorted[, 2]
  df$dominant_margin <- df$dominant_score - df$second_score
  df$high_conf_state_main <- df$dominant_score > 0 & df$dominant_margin >= margin_main
  df$high_conf_state_strict <- df$dominant_score > 0 & df$dominant_margin >= margin_strict
  df
}

cat("Started:", format(Sys.time()), "\n")
cat("R:", R.version.string, "\n")
cat("data.table threads:", data.table::getDTthreads(), "\n")

gene_sets <- readRDS(gene_set_file)
neftel_sets <- gene_sets$neftel[c("AC", "OPC", "NPC1", "NPC2", "MES1", "MES2")]
pathway_sets <- gene_sets$pathways
all_gene_sets <- c(neftel_sets, pathway_sets)
target_genes <- unique(c("LAP3", unlist(all_gene_sets, use.names = FALSE)))

metadata <- read_excel(metadata_file, skip = 43) %>%
  rename(
    cell_id = `Sample name`,
    patient = `tumour name`,
    age_group = `adult/pediatric`,
    processed_file = `processed data file`
  ) %>%
  mutate(
    platform = ifelse(grepl("Smartseq2", processed_file), "Smartseq2", "10X")
  )
stopifnot(nrow(metadata) == 24131, !anyDuplicated(metadata$cell_id))
write_table(metadata, "gse131928_cell_metadata.csv")

process_platform <- function(matrix_file, platform_name) {
  cat("Loading", platform_name, "processed TPM matrix...\n")
  dt <- fread(
    cmd = paste("zcat", shQuote(matrix_file)),
    header = TRUE,
    sep = "\t",
    data.table = TRUE,
    showProgress = TRUE,
    nThread = 16
  )
  gene_col <- names(dt)[1]
  cell_ids <- names(dt)[-1]
  meta <- metadata %>% filter(platform == platform_name)
  stopifnot(length(cell_ids) == nrow(meta), setequal(cell_ids, meta$cell_id))
  meta <- meta[match(cell_ids, meta$cell_id), ]

  genes_present <- dt[[gene_col]]
  expression_dt <- dt[, -1]
  total_tpm <- colSums(expression_dt)
  n_genes_detected <- colSums(expression_dt > 0)
  rm(expression_dt)
  coverage <- bind_rows(lapply(names(all_gene_sets), function(signature) {
    requested <- all_gene_sets[[signature]]
    present <- intersect(requested, genes_present)
    data.frame(
      platform = platform_name,
      signature = signature,
      genes_requested = length(requested),
      genes_present = length(present),
      coverage = length(present) / length(requested),
      missing_genes = paste(setdiff(requested, present), collapse = ";"),
      stringsAsFactors = FALSE
    )
  }))

  selected <- dt[get(gene_col) %in% target_genes]
  expr <- as.matrix(selected[, -1])
  storage.mode(expr) <- "numeric"
  rownames(expr) <- selected[[gene_col]]
  colnames(expr) <- cell_ids
  rm(dt, selected)
  invisible(gc())

  log_expr <- log1p(expr)
  neftel_scores <- score_gene_sets(log_expr, neftel_sets)
  pathway_scores <- score_gene_sets(log_expr, pathway_sets)
  score_df <- data.frame(
    cell_id = cell_ids,
    neftel_scores,
    pathway_scores,
    check.names = FALSE
  ) %>%
    collapse_state_scores()

  margin_distribution <- score_df %>%
    mutate(
      dominant_score_pre = apply(
        as.matrix(across(c(state_AC, state_OPC, state_NPC, state_MES))),
        1,
        max,
        na.rm = TRUE
      ),
      dominant_margin_pre = {
        state_matrix <- as.matrix(across(c(state_AC, state_OPC, state_NPC, state_MES)))
        sorted <- t(apply(state_matrix, 1, sort, decreasing = TRUE))
        sorted[, 1] - sorted[, 2]
      }
    )

  margin_main <- as.numeric(quantile(
    margin_distribution$dominant_margin_pre[margin_distribution$dominant_score_pre > 0],
    probs = 0.50,
    na.rm = TRUE
  ))
  margin_strict <- as.numeric(quantile(
    margin_distribution$dominant_margin_pre[margin_distribution$dominant_score_pre > 0],
    probs = 0.75,
    na.rm = TRUE
  ))

  score_df <- score_df %>%
    add_dominant_state(margin_main = margin_main, margin_strict = margin_strict) %>%
    mutate(
      platform = platform_name,
      patient = meta$patient,
      age_group = meta$age_group,
      lap3_processed_tpm = as.numeric(expr["LAP3", cell_id]),
      lap3_log1p_tpm = log1p(lap3_processed_tpm),
      lap3_detected = lap3_processed_tpm > 0,
      total_tpm = as.numeric(total_tpm[cell_id]),
      n_genes_detected = as.integer(n_genes_detected[cell_id])
    )

  threshold <- data.frame(
    platform = platform_name,
    threshold_set = c("main_median_margin", "strict_q75_margin"),
    dominant_score_rule = "dominant_score > 0 after platform-wise gene z-scoring",
    margin_threshold = c(margin_main, margin_strict),
    source = "GSE131928 platform-specific state-score margin distribution; no LAP3/pathway tuning",
    stringsAsFactors = FALSE
  )

  patient_state <- summarise_patient_state(
    score_df,
    confidence_column = "high_conf_state_main",
    minimum_cells = 1L
  ) %>%
    mutate(eligible_n20 = n_cells >= 20)

  list(
    cell_scores = score_df,
    patient_state = patient_state,
    coverage = coverage,
    threshold = threshold
  )
}

smartseq <- process_platform(smartseq_file, "Smartseq2")
tenx <- process_platform(tenx_file, "10X")

cell_scores <- bind_rows(smartseq$cell_scores, tenx$cell_scores)
patient_state <- bind_rows(smartseq$patient_state, tenx$patient_state)
coverage <- bind_rows(smartseq$coverage, tenx$coverage)
thresholds <- bind_rows(smartseq$threshold, tenx$threshold)

stopifnot(all(coverage$genes_present >= 10), !("LAP3" %in% unlist(all_gene_sets)))
write_table(coverage, "gse131928_cellstate_gene_coverage.csv")
write_table(thresholds, "gse131928_frozen_state_thresholds.csv")
write_table(patient_state, "gse131928_highconf_patient_state_scores.csv")

patient_state_sensitivity <- bind_rows(lapply(
  c("high_conf_state_main", "high_conf_state_strict"),
  function(confidence_column) {
    bind_rows(lapply(c(20L, 50L, 100L), function(minimum_cells) {
      summarise_patient_state(cell_scores, confidence_column, minimum_cells)
    }))
  }
))
write_table(
  patient_state_sensitivity,
  "gse131928_patient_state_threshold_sensitivity.csv"
)

cell_state_summary <- cell_scores %>%
  group_by(platform, age_group, dominant_state) %>%
  summarise(
    n_cells = n(),
    n_patients = n_distinct(patient),
    high_conf_main_cells = sum(high_conf_state_main),
    high_conf_strict_cells = sum(high_conf_state_strict),
    lap3_detection_rate = mean(lap3_detected),
    mean_lap3_log1p_tpm = mean(lap3_log1p_tpm),
    median_margin = median(dominant_margin),
    .groups = "drop"
  )
write_table(cell_state_summary, "gse131928_cell_state_summary.csv")

patient_platform_overlap <- cell_scores %>%
  distinct(platform, patient, age_group) %>%
  group_by(patient) %>%
  summarise(
    n_platforms = n_distinct(platform),
    platforms = paste(sort(unique(platform)), collapse = ";"),
    age_group = paste(sort(unique(age_group)), collapse = ";"),
    .groups = "drop"
  ) %>%
  mutate(overlap_between_platforms = n_platforms > 1)
write_table(patient_platform_overlap, "gse131928_patient_platform_overlap.csv")

pathway_assoc <- patient_state %>%
  filter(age_group == "adult", eligible_n20) %>%
  group_by(platform, dominant_state) %>%
  group_modify(function(.x, .y) {
    bind_rows(lapply(names(pathway_sets), function(pathway) {
      res <- spearman_safe(.x$mean_lap3_log1p_tpm, .x[[pathway]], min_n = 6)
      stat_fun <- function(d) {
        suppressWarnings(cor(d$mean_lap3_log1p_tpm, d[[pathway]], method = "spearman", use = "complete.obs"))
      }
      boot <- if (is.finite(res$rho) && res$n >= 6 && n_distinct(.x$patient) >= 6) {
        cluster_bootstrap(
          .x,
          cluster = "patient",
          statistic = stat_fun,
          replicates = bootstrap_replicates_main
        )
      } else {
        list(ci_low = NA_real_, ci_high = NA_real_)
      }
      lopo <- if (is.finite(res$rho) && res$n >= 6 && n_distinct(.x$patient) >= 6) {
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
write_table(pathway_assoc, "gse131928_highconf_state_pathway_associations.csv")

threshold_pathway_assoc <- run_pathway_associations(
  patient_state_sensitivity %>% filter(age_group == "adult"),
  adjusted_for_depth = FALSE
)
write_table(
  threshold_pathway_assoc,
  "gse131928_pathway_threshold_sensitivity.csv"
)

depth_adjusted_assoc <- run_pathway_associations(
  patient_state_sensitivity %>%
    filter(
      age_group == "adult",
      confidence_rule == "high_conf_state_main",
      minimum_cells == 20L
    ),
  adjusted_for_depth = TRUE
)
write_table(
  depth_adjusted_assoc,
  "gse131928_pathway_depth_adjusted_associations.csv"
)

lopo_detail <- patient_state_sensitivity %>%
  filter(
    age_group == "adult",
    confidence_rule == "high_conf_state_main",
    minimum_cells == 20L
  ) %>%
  group_by(confidence_rule, minimum_cells, platform, dominant_state) %>%
  group_modify(function(.x, .y) {
    bind_rows(lapply(names(pathway_sets), function(pathway) {
      if (n_distinct(.x$patient) < 6L) {
        return(data.frame(
          pathway = pathway,
          omitted_patient = NA_character_,
          lopo_rho = NA_real_
        ))
      }
      statistic <- function(d) {
        suppressWarnings(cor(
          d$mean_lap3_log1p_tpm,
          d[[pathway]],
          method = "spearman",
          use = "complete.obs"
        ))
      }
      leave_one_cluster_out(.x, cluster = "patient", statistic = statistic) %>%
        transmute(
          pathway = pathway,
          omitted_patient = omitted_cluster,
          lopo_rho = estimate
        )
    }))
  }) %>%
  ungroup()
write_table(lopo_detail, "gse131928_pathway_lopo_detail.csv")

depth_detection_summary <- cell_scores %>%
  group_by(platform, age_group, dominant_state) %>%
  summarise(
    n_cells = n(),
    n_patients = n_distinct(patient),
    lap3_detection_rate = mean(lap3_detected),
    mean_nonzero_lap3_log1p_tpm = mean(
      lap3_log1p_tpm[lap3_detected],
      na.rm = TRUE
    ),
    median_detected_genes = median(n_genes_detected),
    median_total_tpm = median(total_tpm),
    .groups = "drop"
  )
write_table(depth_detection_summary, "gse131928_platform_detection_depth_summary.csv")

cell_score_export <- cell_scores %>%
  select(
    platform, cell_id, patient, age_group, lap3_processed_tpm, lap3_log1p_tpm,
    lap3_detected, total_tpm, n_genes_detected, AC, OPC, NPC1, NPC2, MES1, MES2,
    state_AC, state_OPC, state_NPC, state_MES,
    dominant_state, dominant_score, second_score, dominant_margin,
    high_conf_state_main, high_conf_state_strict
  )
fwrite(cell_score_export, file.path(source_dir, "gse131928_cell_neftel_scores.csv.gz"))
saveRDS(
  list(
    thresholds = thresholds,
    cell_scores = cell_score_export,
    patient_state = patient_state
  ),
  file.path(source_dir, "gse131928_cellstate_calibration.rds")
)

cat("Cells scored:", nrow(cell_scores), "\n")
cat("Patients:", n_distinct(cell_scores$patient), "\n")
cat("Cross-platform patient overlaps:", sum(patient_platform_overlap$overlap_between_platforms), "\n")
print(thresholds)
cat("Completed:", format(Sys.time()), "\n")
