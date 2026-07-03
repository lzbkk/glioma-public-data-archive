#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})

setwd("/home/lzb/glioma")
data.table::setDTthreads(16)
set.seed(20260630)

source("Data_scRNA_GEO/scripts/helpers/scRNA_inference_helpers.R")

cache_file <- "Data_scRNA_GEO/GBmap_Core/cache/core_gbmap_lap3_state_union_lightweight.rds"
out_dir <- "Data_scRNA_GEO/GBmap_Core/results/LAP3_State_Union_Projection"
table_dir <- file.path(out_dir, "tables")
source_dir <- file.path(out_dir, "source_data")
log_dir <- file.path(out_dir, "logs")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(log_dir, "lap3_gbmap_core_state_union_projection.log")
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

map_state <- function(annotation_level_3, annotation_level_4) {
  state <- rep(NA_character_, length(annotation_level_3))
  state[annotation_level_3 == "AC-like" | grepl("^AC-like", annotation_level_4)] <- "AC"
  state[annotation_level_3 == "OPC-like" | grepl("^OPC-like", annotation_level_4)] <- "OPC"
  state[annotation_level_3 == "NPC-like" | grepl("^NPC-like", annotation_level_4)] <- "NPC"
  state[annotation_level_3 == "MES-like" | grepl("^MES-like", annotation_level_4)] <- "MES"
  state
}

spearman_safe <- function(x, y, min_n = 6L) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) < min_n || length(unique(x)) < 3L || length(unique(y)) < 3L) {
    return(list(n = length(x), rho = NA_real_, p_value = NA_real_))
  }
  test <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  list(n = length(x), rho = unname(test$estimate), p_value = test$p.value)
}

finite_min <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else min(x)
}

finite_max <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else max(x)
}

score_gene_sets <- function(expr, gene_sets) {
  scores <- matrix(
    NA_real_,
    nrow = nrow(expr),
    ncol = length(gene_sets),
    dimnames = list(rownames(expr), names(gene_sets))
  )
  for (set_name in names(gene_sets)) {
    genes <- intersect(gene_sets[[set_name]], colnames(expr))
    if (length(genes) < 5L) {
      next
    }
    message("Scoring ", set_name, " with ", length(genes), " genes")
    scores[, set_name] <- Matrix::rowMeans(expr[, genes, drop = FALSE])
  }
  as.data.table(scores, keep.rownames = "cell_id")
}

run_within_state_lap3_association <- function(ps) {
  state_sets <- c("LAP3_STATE_UNION", "LAP3_STATE_UNION_NO_TRANSLATION_PROTEOSTASIS")
  out <- rbindlist(lapply(sort(unique(ps$entry_variant)), function(variant) {
    rbindlist(lapply(sort(unique(ps$threshold)), function(threshold_value) {
      rbindlist(lapply(sort(unique(ps$author_state)), function(state) {
        d_state <- ps[
          entry_variant == variant &
            threshold == threshold_value &
            author_state == state
        ]
        rbindlist(lapply(state_sets, function(state_set) {
          res <- spearman_safe(d_state$lap3_mean, d_state[[state_set]], min_n = 6L)
          stat_fun <- function(d) {
            suppressWarnings(cor(
              d$lap3_mean,
              d[[state_set]],
              method = "spearman",
              use = "complete.obs"
            ))
          }
          boot <- if (res$n >= 6L && length(unique(d_state$author_donor)) >= 3L) {
            cluster_bootstrap(d_state, cluster = "author_donor", statistic = stat_fun, replicates = 500L)
          } else {
            list(ci_low = NA_real_, ci_high = NA_real_)
          }
          lopo <- if (res$n >= 6L && length(unique(d_state$author)) >= 3L) {
            leave_one_cluster_out(d_state, cluster = "author", statistic = stat_fun)
          } else {
            data.frame(estimate = NA_real_)
          }
          data.table(
            entry_variant = variant,
            threshold = threshold_value,
            author_state = state,
            state_set = state_set,
            estimand = "within_state_donor_lap3_vs_state_score",
            n_donor_states = res$n,
            n_authors = uniqueN(d_state$author),
            n_donors = uniqueN(d_state$author_donor),
            spearman_rho = res$rho,
            ci_low = boot$ci_low,
            ci_high = boot$ci_high,
            p_value = res$p_value,
            leave_one_author_min_rho = finite_min(lopo$estimate),
            leave_one_author_max_rho = finite_max(lopo$estimate),
            fdr_family = "state_union_lap3_association"
          )
        }))
      }))
    }))
  }), fill = TRUE)
  as.data.frame(out) |>
    adjust_fdr_by_family(p_column = "p_value", family_column = "fdr_family") |>
    as.data.table()
}

run_state_preference_lm <- function(ps) {
  state_sets <- c("LAP3_STATE_UNION", "LAP3_STATE_UNION_NO_TRANSLATION_PROTEOSTASIS")
  out <- rbindlist(lapply(sort(unique(ps$entry_variant)), function(variant) {
    rbindlist(lapply(sort(unique(ps$threshold)), function(threshold_value) {
      d0 <- ps[entry_variant == variant & threshold == threshold_value]
      d0 <- d0[, if (.N >= 2L) .SD, by = author_donor]
      rbindlist(lapply(state_sets, function(state_set) {
        d <- copy(d0)
        if (nrow(d) < 12L || uniqueN(d$author_donor) < 6L || uniqueN(d$author_state) < 2L) {
          return(data.table(
            entry_variant = variant,
            threshold = threshold_value,
            state_set = state_set,
            model = "state_score ~ author_state + author_donor_fixed_effect",
            term = NA_character_,
            estimate = NA_real_,
            p_value = NA_real_,
            n_rows = nrow(d),
            n_donors = uniqueN(d$author_donor)
          ))
        }
        d[, author_state := factor(author_state, levels = c("AC", "OPC", "NPC", "MES"))]
        d[, state_score := get(state_set)]
        fit <- lm(state_score ~ author_state + author_donor, data = d)
        coef_table <- as.data.table(coef(summary(fit)), keep.rownames = "term")
        setnames(coef_table, c("Estimate", "Std. Error", "t value", "Pr(>|t|)"),
                 c("estimate", "std_error", "t_value", "p_value"))
        coef_table[grepl("^author_state", term),
          .(
            entry_variant = variant,
            threshold = threshold_value,
            state_set = state_set,
            model = "state_score ~ author_state + author_donor_fixed_effect",
            term,
            estimate,
            std_error,
            t_value,
            p_value,
            n_rows = nrow(d),
            n_donors = uniqueN(d$author_donor)
          )
        ]
      }))
    }))
  }), fill = TRUE)
  out[, p_adj_BH := p.adjust(p_value, method = "BH")]
  out[]
}

cat("Started:", format(Sys.time()), "\n")
cache <- readRDS(cache_file)
stopifnot(all(c("obs", "normalized", "raw", "gene_sets", "coverage") %in% names(cache)))
obs <- as.data.table(cache$obs)
expr <- cache$normalized
raw <- cache$raw
stopifnot(nrow(obs) == nrow(expr), identical(rownames(expr), obs$`_index`))
stopifnot("LAP3" %in% colnames(expr), "LAP3" %in% colnames(raw))

gene_sets <- lapply(cache$gene_sets, function(x) sort(setdiff(x, "LAP3")))
obs[, author_donor := paste(author, donor_id, sep = "::")]
obs[, author_state := map_state(annotation_level_3, annotation_level_4)]
obs[, main_entry := author != "Neftel2019" & annotation_level_1 == "Neoplastic" & !is.na(author_state)]
obs[, strict_entry := main_entry & iCNV == "aneuploid"]

entry_counts <- rbindlist(list(
  obs[, .(entry_variant = "all_core", n_cells = .N, n_authors = uniqueN(author), n_donors = uniqueN(author_donor))],
  obs[author != "Neftel2019", .(entry_variant = "exclude_neftel2019", n_cells = .N, n_authors = uniqueN(author), n_donors = uniqueN(author_donor))],
  obs[main_entry == TRUE, .(entry_variant = "main_neoplastic_exclude_neftel2019", n_cells = .N, n_authors = uniqueN(author), n_donors = uniqueN(author_donor))],
  obs[strict_entry == TRUE, .(entry_variant = "strict_neoplastic_aneuploid_exclude_neftel2019", n_cells = .N, n_authors = uniqueN(author), n_donors = uniqueN(author_donor))]
), fill = TRUE)
write_table(entry_counts, "gbmap_core_lap3_state_entry_counts.csv")

analysis_index <- which(obs$main_entry)
cat("Main analysis cells:", length(analysis_index), "\n")
expr_main <- expr[analysis_index, , drop = FALSE]
raw_main <- raw[analysis_index, , drop = FALSE]
obs_main <- copy(obs[analysis_index])

coverage <- rbindlist(lapply(names(gene_sets), function(state_set) {
  requested <- gene_sets[[state_set]]
  present <- intersect(requested, colnames(expr_main))
  data.table(
    state_set = state_set,
    requested = length(requested),
    available = length(present),
    coverage = length(present) / length(requested),
    missing_genes = paste(setdiff(requested, present), collapse = ";")
  )
}))
write_table(coverage, "gbmap_core_lap3_state_gene_coverage.csv")

scores <- score_gene_sets(expr_main, gene_sets)
score_cols <- setdiff(names(scores), "cell_id")
cell_dt <- cbind(
  obs_main[, .(
    cell_id = `_index`,
    author,
    donor_id,
    sample,
    assay,
    annotation_level_1,
    annotation_level_3,
    annotation_level_4,
    iCNV,
    author_donor,
    author_state,
    strict_entry
  )],
  scores[, ..score_cols]
)
cell_dt[, lap3_log_norm := as.numeric(expr_main[, "LAP3"])]
cell_dt[, lap3_raw := as.numeric(raw_main[, "LAP3"])]
cell_dt[, lap3_detected := lap3_raw > 0]
cell_dt[, state_target_raw_sum := Matrix::rowSums(raw_main)]
cell_dt[, entry_variant := "main_neoplastic_exclude_neftel2019"]
strict_dt <- copy(cell_dt[strict_entry == TRUE])
strict_dt[, entry_variant := "strict_neoplastic_aneuploid_exclude_neftel2019"]
cell_dt <- rbindlist(list(cell_dt, strict_dt), use.names = TRUE, fill = TRUE)

fwrite(
  cell_dt[, .(
    cell_id, entry_variant, author, donor_id, sample, assay,
    annotation_level_3, annotation_level_4, iCNV, author_state,
    lap3_log_norm, lap3_raw, lap3_detected,
    LAP3_STATE_UNION, LAP3_STATE_UNION_NO_TRANSLATION_PROTEOSTASIS,
    state_target_raw_sum
  )],
  file.path(source_dir, "gbmap_core_lap3_state_cell_scores.csv.gz")
)

state_cell_summary <- cell_dt[
  ,
  .(
    n_cells = .N,
    n_authors = uniqueN(author),
    n_donors = uniqueN(author_donor),
    lap3_mean = mean(lap3_log_norm),
    lap3_detection_rate = mean(lap3_detected),
    LAP3_STATE_UNION = mean(LAP3_STATE_UNION),
    LAP3_STATE_UNION_NO_TRANSLATION_PROTEOSTASIS = mean(LAP3_STATE_UNION_NO_TRANSLATION_PROTEOSTASIS)
  ),
  by = .(entry_variant, author_state)
][order(entry_variant, author_state)]
write_table(state_cell_summary, "gbmap_core_lap3_state_cell_summary.csv")

thresholds <- c(20L, 50L, 100L)
patient_state <- rbindlist(lapply(thresholds, function(threshold_value) {
  cell_dt[
    ,
    .(
      n_cells = .N,
      lap3_mean = mean(lap3_log_norm),
      lap3_median = median(lap3_log_norm),
      lap3_detection_rate = mean(lap3_detected),
      lap3_raw_mean = mean(lap3_raw),
      state_target_raw_sum_mean = mean(state_target_raw_sum),
      LAP3_STATE_UNION = mean(LAP3_STATE_UNION),
      LAP3_STATE_UNION_NO_TRANSLATION_PROTEOSTASIS = mean(LAP3_STATE_UNION_NO_TRANSLATION_PROTEOSTASIS)
    ),
    by = .(entry_variant, author, donor_id, author_donor, author_state)
  ][n_cells >= threshold_value][, threshold := threshold_value][]
}), fill = TRUE)
write_table(patient_state, "gbmap_core_lap3_state_patient_state_summary.csv")

threshold_retention <- patient_state[
  ,
  .(
    n_donor_states = .N,
    n_authors = uniqueN(author),
    n_donors = uniqueN(author_donor),
    median_cells = as.numeric(median(n_cells)),
    min_cells = as.numeric(min(n_cells))
  ),
  by = .(entry_variant, threshold, author_state)
][order(entry_variant, threshold, author_state)]
write_table(threshold_retention, "gbmap_core_lap3_state_threshold_retention.csv")

lap3_associations <- run_within_state_lap3_association(patient_state)
write_table(lap3_associations, "gbmap_core_lap3_state_within_state_lap3_associations.csv")

state_preference <- run_state_preference_lm(patient_state)
write_table(state_preference, "gbmap_core_lap3_state_preference_fixed_effect.csv")

primary_summary <- lap3_associations[
  entry_variant == "main_neoplastic_exclude_neftel2019" &
    threshold == 20L &
    state_set == "LAP3_STATE_UNION"
][order(author_state)]
write_table(primary_summary, "gbmap_core_lap3_state_primary_summary.csv")

readme <- c(
  "# Core GBmap LAP3-State Union Projection",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Input",
  "",
  paste0("- Cache: `", cache_file, "`"),
  "- Frozen gene sets: `Data_Bulk_TCGA/Data_Merged/results/LAP3_State_Module/tables/lap3_state_frozen_gene_sets.csv`",
  "- Main entry: `author != Neftel2019`, `annotation_level_1 == Neoplastic`, mapped AC/OPC/NPC/MES state.",
  "- Strict sensitivity: main entry plus `iCNV == aneuploid`.",
  "",
  "## Methods",
  "",
  "- `LAP3_STATE_UNION` and `LAP3_STATE_UNION_NO_TRANSLATION_PROTEOSTASIS` are scored as mean normalized expression across available state genes.",
  "- LAP3 itself is excluded from both state gene sets and retained only as an anchor/readout.",
  "- Cells are summarized to `author x donor x state`; donor-state summaries are the inference units.",
  "- Main threshold is 20 cells per donor-state; 50 and 100 are sensitivity thresholds.",
  "- State preference uses donor fixed effects: `state_score ~ author_state + author_donor`.",
  "",
  "## Key Outputs",
  "",
  "- `tables/gbmap_core_lap3_state_gene_coverage.csv`",
  "- `tables/gbmap_core_lap3_state_patient_state_summary.csv`",
  "- `tables/gbmap_core_lap3_state_within_state_lap3_associations.csv`",
  "- `tables/gbmap_core_lap3_state_preference_fixed_effect.csv`",
  "- `source_data/gbmap_core_lap3_state_cell_scores.csv.gz`",
  "",
  "## Interpretation Boundary",
  "",
  "This module projects the frozen bulk-derived LAP3-state into Core GBmap malignant cells.",
  "It is a cell-state decomposition and context analysis, not a causal LAP3 mechanism test."
)
writeLines(readme, file.path(out_dir, "README.md"))

cat("Finished:", format(Sys.time()), "\n")
