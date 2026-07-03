#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})

setwd("/home/lzb/glioma")
data.table::setDTthreads(16)
set.seed(20260630)

source("Data_scRNA_GEO/scripts/helpers/scRNA_inference_helpers.R")

cache_file <- "Data_scRNA_GEO/GBmap_Core/cache/core_gbmap_lap3_cellstate_lightweight.rds"
out_dir <- "Data_scRNA_GEO/GBmap_Core/results/LAP3_CellState_Strict"
table_dir <- file.path(out_dir, "tables")
source_dir <- file.path(out_dir, "source_data")
log_dir <- file.path(out_dir, "logs")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(log_dir, "lap3_gbmap_core_cellstate_strict.log")
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

map_state <- function(annotation_level_3, annotation_level_4) {
  state <- rep(NA_character_, length(annotation_level_3))
  state[annotation_level_3 == "AC-like" | grepl("^AC-like", annotation_level_4)] <- "AC"
  state[annotation_level_3 == "OPC-like" | grepl("^OPC-like", annotation_level_4)] <- "OPC"
  state[annotation_level_3 == "NPC-like" | grepl("^NPC-like", annotation_level_4)] <- "NPC"
  state[annotation_level_3 == "MES-like" | grepl("^MES-like", annotation_level_4)] <- "MES"
  state
}

score_gene_sets <- function(expr, gene_sets) {
  scores <- matrix(
    NA_real_,
    nrow = nrow(expr),
    ncol = length(gene_sets),
    dimnames = list(rownames(expr), names(gene_sets))
  )
  for (signature in names(gene_sets)) {
    genes <- intersect(gene_sets[[signature]], colnames(expr))
    if (length(genes) < 5L) {
      next
    }
    message("Scoring ", signature, " with ", length(genes), " genes")
    x <- as.matrix(expr[, genes, drop = FALSE])
    z <- scale(x)
    z[!is.finite(z)] <- 0
    scores[, signature] <- rowMeans(z)
    rm(x, z)
    invisible(gc())
  }
  as.data.table(scores, keep.rownames = "cell_id")
}

summarise_patient_state <- function(cell_dt, variant, threshold) {
  cell_dt[
    entry_variant == variant & !is.na(author_state),
    .(
      n_cells = .N,
      lap3_mean = mean(lap3_log_norm),
      lap3_median = median(lap3_log_norm),
      lap3_detection_rate = mean(lap3_detected),
      lap3_raw_mean = mean(lap3_raw),
      target_raw_sum_mean = mean(target_raw_sum),
      state_AC = mean(state_AC),
      state_OPC = mean(state_OPC),
      state_NPC = mean(state_NPC),
      state_MES = mean(state_MES),
      HALLMARK_MTORC1_SIGNALING = mean(HALLMARK_MTORC1_SIGNALING),
      LEUCINE_BCAA_CORE = mean(LEUCINE_BCAA_CORE),
      MTORC1_READOUT_CORE = mean(MTORC1_READOUT_CORE),
      REACTOME_TRANSLATION = mean(REACTOME_TRANSLATION)
    ),
    by = .(entry_variant, author, donor_id, author_donor, author_state)
  ][n_cells >= threshold][, threshold := threshold][]
}

run_within_state_pathway <- function(ps, threshold) {
  pathways <- c(
    "HALLMARK_MTORC1_SIGNALING",
    "LEUCINE_BCAA_CORE",
    "MTORC1_READOUT_CORE",
    "REACTOME_TRANSLATION"
  )
  out <- rbindlist(lapply(sort(unique(ps$entry_variant)), function(variant) {
    rbindlist(lapply(sort(unique(ps$author_state)), function(state) {
      d_state <- ps[entry_variant == variant & author_state == state]
      rbindlist(lapply(pathways, function(pathway) {
        res <- spearman_safe(d_state$lap3_mean, d_state[[pathway]], min_n = 6L)
        stat_fun <- function(d) {
          suppressWarnings(cor(
            d$lap3_mean,
            d[[pathway]],
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
          threshold = threshold,
          estimand = "within_state_across_donors",
          author_state = state,
          pathway = pathway,
          fdr_family = ifelse(pathway %in% c("HALLMARK_MTORC1_SIGNALING", "LEUCINE_BCAA_CORE"), "primary", "secondary"),
          n_donor_states = res$n,
          n_authors = uniqueN(d_state$author),
          n_donors = uniqueN(d_state$author_donor),
          spearman_rho = res$rho,
          ci_low = boot$ci_low,
          ci_high = boot$ci_high,
          p_value = res$p_value,
          leave_one_author_min_rho = finite_min(lopo$estimate),
          leave_one_author_max_rho = finite_max(lopo$estimate)
        )
      }))
    }))
  }), fill = TRUE)
  as.data.frame(out) |>
    adjust_fdr_by_family(p_column = "p_value", family_column = "fdr_family") |>
    as.data.table()
}

run_continuous_state_association <- function(donor_dt) {
  state_scores <- c("state_AC", "state_OPC", "state_NPC", "state_MES")
  out <- rbindlist(lapply(sort(unique(donor_dt$entry_variant)), function(variant) {
    d_variant <- donor_dt[entry_variant == variant]
    rbindlist(lapply(state_scores, function(score) {
      res <- spearman_safe(d_variant$lap3_mean, d_variant[[score]], min_n = 8L)
      data.table(
        entry_variant = variant,
        estimand = "donor_level_continuous_state_score",
        state_score = score,
        n_donors = res$n,
        n_authors = uniqueN(d_variant$author),
        spearman_rho = res$rho,
        p_value = res$p_value
      )
    }))
  }), fill = TRUE)
  out[, fdr_family := "state_continuum"]
  as.data.frame(out) |>
    adjust_fdr_by_family(p_column = "p_value", family_column = "fdr_family") |>
    as.data.table()
}

run_state_preference_lm <- function(ps, threshold) {
  rbindlist(lapply(sort(unique(ps$entry_variant)), function(variant) {
    d <- copy(ps[entry_variant == variant])
    d <- d[, if (.N >= 2L) .SD, by = author_donor]
    if (nrow(d) < 12L || uniqueN(d$author_donor) < 6L || uniqueN(d$author_state) < 2L) {
      return(data.table(
        entry_variant = variant,
        threshold = threshold,
        model = "lap3_mean ~ author_state + author_donor_fixed_effect",
        term = NA_character_,
        estimate = NA_real_,
        p_value = NA_real_,
        n_rows = nrow(d),
        n_donors = uniqueN(d$author_donor)
      ))
    }
    d[, author_state := factor(author_state, levels = c("AC", "OPC", "NPC", "MES"))]
    fit <- lm(lap3_mean ~ author_state + author_donor, data = d)
    coef_table <- as.data.table(coef(summary(fit)), keep.rownames = "term")
    setnames(coef_table, c("Estimate", "Std. Error", "t value", "Pr(>|t|)"),
             c("estimate", "std_error", "t_value", "p_value"))
    coef_table[grepl("^author_state", term),
      .(
        entry_variant = variant,
        threshold = threshold,
        model = "lap3_mean ~ author_state + author_donor_fixed_effect",
        term,
        estimate,
        std_error,
        t_value,
        p_value,
        n_rows = nrow(d),
        n_donors = uniqueN(d$author_donor)
      )
    ]
  }), fill = TRUE)
}

cat("Started:", format(Sys.time()), "\n")
cat("R:", R.version.string, "\n")
cat("data.table threads:", data.table::getDTthreads(), "\n")

cache <- readRDS(cache_file)
stopifnot(all(c("obs", "normalized", "raw", "gene_sets", "coverage") %in% names(cache)))
obs <- as.data.table(cache$obs)
expr <- cache$normalized
raw <- cache$raw
stopifnot(nrow(obs) == nrow(expr), identical(rownames(expr), obs$`_index`))
stopifnot("LAP3" %in% colnames(expr), "LAP3" %in% colnames(raw))

gene_sets <- cache$gene_sets
gene_sets <- gene_sets[names(gene_sets) %in% c(
  "AC", "OPC", "NPC1", "NPC2", "MES1", "MES2",
  "HALLMARK_MTORC1_SIGNALING", "LEUCINE_BCAA_CORE",
  "MTORC1_READOUT_CORE", "REACTOME_TRANSLATION"
)]
gene_sets <- lapply(gene_sets, function(x) setdiff(x, "LAP3"))

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
write_table(entry_counts, "gbmap_core_entry_counts.csv")

analysis_index <- which(obs$main_entry)
cat("Main analysis cells:", length(analysis_index), "\n")
expr_main <- expr[analysis_index, , drop = FALSE]
raw_main <- raw[analysis_index, , drop = FALSE]
obs_main <- copy(obs[analysis_index])

coverage <- rbindlist(lapply(names(gene_sets), function(signature) {
  requested <- gene_sets[[signature]]
  present <- intersect(requested, colnames(expr_main))
  data.table(
    gene_set = signature,
    requested = length(requested),
    available = length(present),
    coverage = length(present) / length(requested),
    missing_genes = paste(setdiff(requested, present), collapse = ";")
  )
}))
write_table(coverage, "gbmap_core_cellstate_gene_coverage.csv")
stopifnot(all(coverage$available >= 10L))

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
    annotation_level_2,
    annotation_level_3,
    annotation_level_4,
    iCNV,
    author_donor,
    author_state,
    strict_entry
  )],
  scores[, ..score_cols]
)

cell_dt[, state_AC := AC]
cell_dt[, state_OPC := OPC]
cell_dt[, state_NPC := rowMeans(.SD, na.rm = TRUE), .SDcols = c("NPC1", "NPC2")]
cell_dt[, state_MES := rowMeans(.SD, na.rm = TRUE), .SDcols = c("MES1", "MES2")]
cell_dt[, lap3_log_norm := as.numeric(expr_main[, "LAP3"])]
cell_dt[, lap3_raw := as.numeric(raw_main[, "LAP3"])]
cell_dt[, lap3_detected := lap3_raw > 0]
cell_dt[, target_raw_sum := Matrix::rowSums(raw_main)]
cell_dt[, entry_variant := "main_neoplastic_exclude_neftel2019"]
strict_dt <- copy(cell_dt[strict_entry == TRUE])
strict_dt[, entry_variant := "strict_neoplastic_aneuploid_exclude_neftel2019"]
cell_dt <- rbindlist(list(cell_dt, strict_dt), use.names = TRUE, fill = TRUE)

fwrite(
  cell_dt[, .(
    cell_id, entry_variant, author, donor_id, sample, assay,
    annotation_level_3, annotation_level_4, iCNV, author_state,
    lap3_log_norm, lap3_raw, lap3_detected,
    state_AC, state_OPC, state_NPC, state_MES,
    HALLMARK_MTORC1_SIGNALING, LEUCINE_BCAA_CORE,
    MTORC1_READOUT_CORE, REACTOME_TRANSLATION
  )],
  file.path(source_dir, "gbmap_core_cell_scores.csv.gz")
)

state_cell_summary <- cell_dt[
  ,
  .(
    n_cells = .N,
    n_authors = uniqueN(author),
    n_donors = uniqueN(author_donor),
    lap3_mean = mean(lap3_log_norm),
    lap3_detection_rate = mean(lap3_detected)
  ),
  by = .(entry_variant, author_state)
][order(entry_variant, author_state)]
write_table(state_cell_summary, "gbmap_core_state_cell_summary.csv")

thresholds <- c(20L, 50L, 100L)
patient_state_all <- rbindlist(lapply(thresholds, function(th) {
  rbindlist(lapply(c(
    "main_neoplastic_exclude_neftel2019",
    "strict_neoplastic_aneuploid_exclude_neftel2019"
  ), function(variant) summarise_patient_state(cell_dt, variant, th)), fill = TRUE)
}), fill = TRUE)
write_table(patient_state_all, "gbmap_core_patient_state_pseudobulk.csv")

threshold_retention <- patient_state_all[
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
write_table(threshold_retention, "gbmap_core_patient_state_threshold_retention.csv")

pathway_results <- rbindlist(lapply(thresholds, function(th) {
  run_within_state_pathway(patient_state_all[threshold == th], th)
}), fill = TRUE)
write_table(pathway_results, "gbmap_core_within_state_pathway_associations.csv")

state_preference <- rbindlist(lapply(thresholds, function(th) {
  run_state_preference_lm(patient_state_all[threshold == th], th)
}), fill = TRUE)
if ("p_value" %in% names(state_preference)) {
  state_preference[, p_adj_BH := p.adjust(p_value, method = "BH")]
}
write_table(state_preference, "gbmap_core_lap3_state_preference_fixed_effect.csv")

donor_summary <- cell_dt[
  ,
  .(
    n_cells = .N,
    lap3_mean = mean(lap3_log_norm),
    lap3_detection_rate = mean(lap3_detected),
    state_AC = mean(state_AC),
    state_OPC = mean(state_OPC),
    state_NPC = mean(state_NPC),
    state_MES = mean(state_MES)
  ),
  by = .(entry_variant, author, donor_id, author_donor)
][n_cells >= 50L]
continuous_results <- run_continuous_state_association(donor_summary)
write_table(donor_summary, "gbmap_core_donor_continuous_state_scores.csv")
write_table(continuous_results, "gbmap_core_lap3_continuous_state_associations.csv")

primary_summary <- pathway_results[
  threshold == 20L &
    entry_variant == "main_neoplastic_exclude_neftel2019" &
    pathway %in% c("HALLMARK_MTORC1_SIGNALING", "LEUCINE_BCAA_CORE")
][order(author_state, pathway)]
write_table(primary_summary, "gbmap_core_primary_result_summary.csv")

readme <- c(
  "# Core GBmap LAP3 Cell-State Strict Analysis",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Input",
  "",
  paste0("- Cache: `", cache_file, "`"),
  "- Main entry: `author != Neftel2019`, `annotation_level_1 == Neoplastic`, mapped AC/OPC/NPC/MES state.",
  "- Strict sensitivity: main entry plus `iCNV == aneuploid`.",
  "",
  "## Methods",
  "",
  "- Scores use frozen Neftel and pathway gene sets from Phase 0; LAP3 is excluded from pathway sets.",
  "- Cells are first summarized to `author × donor × state`; donor/state summaries are the inference units.",
  "- Main within-state pathway association uses Spearman correlation across donor-state summaries with donor-cluster bootstrap CIs and leave-one-author sensitivity.",
  "- LAP3 state preference is estimated with donor fixed effects: `lap3_mean ~ author_state + author_donor`.",
  "- Thresholds of 20, 50, and 100 cells per donor-state are reported.",
  "",
  "## Key Outputs",
  "",
  "- `tables/gbmap_core_entry_counts.csv`",
  "- `tables/gbmap_core_patient_state_pseudobulk.csv`",
  "- `tables/gbmap_core_within_state_pathway_associations.csv`",
  "- `tables/gbmap_core_lap3_state_preference_fixed_effect.csv`",
  "- `tables/gbmap_core_lap3_continuous_state_associations.csv`",
  "- `source_data/gbmap_core_cell_scores.csv.gz`",
  "",
  "## Interpretation Boundary",
  "",
  "This module is a patient-level association analysis in a harmonized public atlas.",
  "It does not prove LAP3 enzymatic causality, leucine availability, or direct mTORC1 activation."
)
writeLines(readme, file.path(out_dir, "README.md"))

cat("Finished:", format(Sys.time()), "\n")
