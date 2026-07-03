#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

setwd("/home/lzb/glioma")
data.table::setDTthreads(16)

input_file <- "Data_scRNA_GEO/GBmap_Core/results/LAP3_CellState_Strict/source_data/gbmap_core_cell_scores.csv.gz"
out_dir <- "Data_scRNA_GEO/GBmap_Core/results/LAP3_Substate_Profiling"
table_dir <- file.path(out_dir, "tables")
source_dir <- file.path(out_dir, "source_data")
log_dir <- file.path(out_dir, "logs")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(log_dir, "lap3_gbmap_core_substate_profiling.log")
cat("Started: ", format(Sys.time()), "\n", file = log_file)

score_cols <- c(
  "HALLMARK_MTORC1_SIGNALING",
  "LEUCINE_BCAA_CORE",
  "MTORC1_READOUT_CORE",
  "REACTOME_TRANSLATION",
  "state_AC",
  "state_OPC",
  "state_NPC",
  "state_MES"
)

wilcox_safe <- function(x, y) {
  if (sum(is.finite(x)) < 3L || sum(is.finite(y)) < 3L) return(NA_real_)
  suppressWarnings(wilcox.test(x, y, exact = FALSE)$p.value)
}

signed_rank_safe <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 4L || length(unique(x)) < 2L) return(NA_real_)
  suppressWarnings(wilcox.test(x, mu = 0, exact = FALSE)$p.value)
}

message("Reading cell score table")
cell <- fread(input_file)
cell <- cell[entry_variant == "main_neoplastic_exclude_neftel2019"]
cell[, author_donor := paste(author, donor_id, sep = "::")]
cell[, donor_state := paste(author, donor_id, author_state, sep = "::")]

message("Calling LAP3-high and LAP3-detected substates within state")
state_cutoffs <- cell[
  ,
  .(
    lap3_state_median = median(lap3_log_norm, na.rm = TRUE),
    lap3_state_q75 = quantile(lap3_log_norm, 0.75, na.rm = TRUE),
    cells = .N,
    donors = uniqueN(author_donor),
    authors = uniqueN(author)
  ),
  by = author_state
]
fwrite(state_cutoffs, file.path(table_dir, "gbmap_core_state_lap3_cutoffs.csv"))

cell <- merge(cell, state_cutoffs[, .(author_state, lap3_state_median, lap3_state_q75)],
              by = "author_state", all.x = TRUE, sort = FALSE)
cell[, lap3_group_state_median := fifelse(lap3_log_norm > lap3_state_median, "LAP3_high", "LAP3_low")]
cell[, lap3_group_state_q75 := fifelse(lap3_log_norm >= lap3_state_q75, "LAP3_top_q75", "LAP3_non_top_q75")]
cell[, lap3_group_detected := fifelse(lap3_detected, "LAP3_detected", "LAP3_not_detected")]

cell_group_summary <- rbindlist(lapply(c(
  "lap3_group_state_median",
  "lap3_group_state_q75",
  "lap3_group_detected"
), function(group_col) {
  cell[, c(
    .(
      grouping = group_col,
      cells = .N,
      donors = uniqueN(author_donor),
      authors = uniqueN(author),
      lap3_mean = mean(lap3_log_norm),
      lap3_raw_mean = mean(lap3_raw)
    ),
    lapply(.SD, mean)
  ), by = .(author_state, group = get(group_col)), .SDcols = score_cols]
}), fill = TRUE)
fwrite(cell_group_summary, file.path(table_dir, "gbmap_core_lap3_substate_cell_group_summary.csv"))

message("Building donor-state paired high-low summaries")
paired <- rbindlist(lapply(c(
  "lap3_group_state_median",
  "lap3_group_state_q75",
  "lap3_group_detected"
), function(group_col) {
  group_dt <- cell[
    ,
    c(
      .(cells = .N, lap3_mean = mean(lap3_log_norm), lap3_raw_mean = mean(lap3_raw)),
      lapply(.SD, mean)
    ),
    by = .(author, donor_id, author_donor, author_state, group = get(group_col)),
    .SDcols = score_cols
  ]
  dcast(
    group_dt,
    author + donor_id + author_donor + author_state ~ group,
    value.var = c("cells", "lap3_mean", "lap3_raw_mean", score_cols)
  )[, grouping := group_col][]
}), fill = TRUE)
fwrite(paired, file.path(source_dir, "gbmap_core_lap3_substate_donor_state_group_means.csv"))

paired_effect_rows <- list()
k <- 0L
group_pairs <- list(
  lap3_group_state_median = c("LAP3_high", "LAP3_low"),
  lap3_group_state_q75 = c("LAP3_top_q75", "LAP3_non_top_q75"),
  lap3_group_detected = c("LAP3_detected", "LAP3_not_detected")
)
for (grp in names(group_pairs)) {
  high <- group_pairs[[grp]][1]
  low <- group_pairs[[grp]][2]
  for (state in sort(unique(paired$author_state))) {
    d <- paired[grouping == grp & author_state == state]
    for (score in score_cols) {
      high_col <- paste(score, high, sep = "_")
      low_col <- paste(score, low, sep = "_")
      if (!all(c(high_col, low_col) %in% names(d))) next
      delta <- d[[high_col]] - d[[low_col]]
      keep <- is.finite(delta)
      k <- k + 1L
      paired_effect_rows[[k]] <- data.table(
        grouping = grp,
        author_state = state,
        score = score,
        n_donor_states = sum(keep),
        n_authors = uniqueN(d$author[keep]),
        median_delta = median(delta[keep], na.rm = TRUE),
        mean_delta = mean(delta[keep], na.rm = TRUE),
        p_value = signed_rank_safe(delta),
        positive_fraction = mean(delta[keep] > 0)
      )
    }
  }
}
paired_effect <- rbindlist(paired_effect_rows, fill = TRUE)
paired_effect[, fdr_family := fifelse(
  score %in% c("HALLMARK_MTORC1_SIGNALING", "LEUCINE_BCAA_CORE"),
  "primary_pathway_substate",
  "secondary_or_state_substate"
)]
paired_effect[, p_adj_BH := p.adjust(p_value, method = "BH"), by = .(grouping, fdr_family)]
fwrite(paired_effect, file.path(table_dir, "gbmap_core_lap3_substate_paired_effects.csv"))

message("Author-level contribution table")
author_effect <- rbindlist(lapply(names(group_pairs), function(grp) {
  high <- group_pairs[[grp]][1]
  low <- group_pairs[[grp]][2]
  rbindlist(lapply(sort(unique(paired$author_state)), function(state) {
    d <- paired[grouping == grp & author_state == state]
    rbindlist(lapply(score_cols, function(score) {
      high_col <- paste(score, high, sep = "_")
      low_col <- paste(score, low, sep = "_")
      if (!all(c(high_col, low_col) %in% names(d))) return(NULL)
      d[, .(
        n_donor_states = .N,
        median_delta = median(get(high_col) - get(low_col), na.rm = TRUE),
        mean_delta = mean(get(high_col) - get(low_col), na.rm = TRUE)
      ), by = author][
        ,
        `:=`(grouping = grp, author_state = state, score = score)
      ]
    }), fill = TRUE)
  }), fill = TRUE)
}), fill = TRUE)
setcolorder(author_effect, c("grouping", "author_state", "score", "author"))
fwrite(author_effect, file.path(table_dir, "gbmap_core_lap3_substate_author_effects.csv"))

primary_matrix <- paired_effect[
  grouping == "lap3_group_state_median" &
    score %in% c("HALLMARK_MTORC1_SIGNALING", "LEUCINE_BCAA_CORE", "MTORC1_READOUT_CORE", "REACTOME_TRANSLATION")
][order(author_state, score)]
fwrite(primary_matrix, file.path(table_dir, "gbmap_core_lap3_substate_primary_matrix.csv"))

readme <- c(
  "# Core GBmap LAP3 Substate Profiling",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Purpose",
  "",
  "Explore whether LAP3-high cells within each malignant state carry deeper anabolic/translation substate information.",
  "",
  "## Input",
  "",
  paste0("- `", input_file, "`"),
  "",
  "## Methods",
  "",
  "- Uses Core GBmap cells after excluding Neftel2019 and retaining author-neoplastic malignant states.",
  "- Defines LAP3-high/low within each author_state by state-level median, top quartile, and detection.",
  "- Main inference table uses donor-state paired high-low score differences rather than cell-level p-values.",
  "- Outputs author-level contribution tables to detect study-structured effects.",
  "",
  "## Key Outputs",
  "",
  "- `tables/gbmap_core_lap3_substate_paired_effects.csv`",
  "- `tables/gbmap_core_lap3_substate_primary_matrix.csv`",
  "- `tables/gbmap_core_lap3_substate_author_effects.csv`",
  "- `source_data/gbmap_core_lap3_substate_donor_state_group_means.csv`",
  "",
  "## Boundary",
  "",
  "This is exploratory substate profiling. It nominates LAP3-high malignant substates but does not prove LAP3 enzymatic causality."
)
writeLines(readme, file.path(out_dir, "README.md"))

cat("Finished: ", format(Sys.time()), "\n", file = log_file, append = TRUE)
message("Done")
