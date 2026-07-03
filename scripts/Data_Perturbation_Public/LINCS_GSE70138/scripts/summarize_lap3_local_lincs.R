#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

result_dir <- "Data_Perturbation_Public/LINCS_GSE70138/results/LAP3_Local_Connectivity"
table_dir <- file.path(result_dir, "tables")
comparison_path <- file.path(
  table_dir,
  "lap3_lincs_full_vs_no_lap3_compounds.csv"
)
stopifnot(file.exists(comparison_path))

comparison <- fread(comparison_path)
comparison[, `:=`(
  worst_case_median_ncs = pmax(median_ncs_full, median_ncs_no_lap3),
  best_case_median_ncs = pmin(median_ncs_full, median_ncs_no_lap3),
  min_fraction_negative = pmin(
    fraction_negative_full,
    fraction_negative_no_lap3
  ),
  min_median_tas = pmin(median_tas_full, median_tas_no_lap3),
  min_cells = pmin(n_cells_full, n_cells_no_lap3),
  min_signatures = pmin(n_signatures_full, n_signatures_no_lap3)
)]

comparison[, stable_reverse := (
  min_signatures >= 10L &
    min_cells >= 5L &
    worst_case_median_ncs <= -0.5 &
    min_fraction_negative >= 0.6 &
    min_median_tas >= 0.2
)]
setorder(
  comparison,
  -stable_reverse,
  worst_case_median_ncs,
  -min_fraction_negative,
  -min_cells
)
comparison[, consensus_reverse_rank := seq_len(.N)]

fwrite(
  comparison,
  file.path(table_dir, "lap3_lincs_consensus_compound_ranking.csv")
)
fwrite(
  comparison[stable_reverse == TRUE],
  file.path(table_dir, "lap3_lincs_stable_reverse_candidates.csv")
)

mechanism_patterns <- list(
  "PI3K-mTOR" = paste(
    c(
      "torin", "BGT-226", "GDC-0980", "NVP-BEZ235", "MLN-0128",
      "PI-103", "sirolimus", "everolimus", "temsirolimus"
    ),
    collapse = "|"
  ),
  "translation" = "homoharringtonine|bruceantin",
  "proteostasis" = "bortezomib|MG-132",
  "BET-epigenetic" = "JQ1|romidepsin|dacinostat|JNJ-26481585",
  "HSP90" = "BIIB-021|XL-888"
)

mechanism_focused <- rbindlist(lapply(names(mechanism_patterns), function(group) {
  comparison[
    grepl(mechanism_patterns[[group]], pert_iname, ignore.case = TRUE)
  ][
    , mechanism_group := group
  ]
}), use.names = TRUE)
setorder(
  mechanism_focused,
  mechanism_group,
  worst_case_median_ncs,
  -min_fraction_negative
)
fwrite(
  mechanism_focused,
  file.path(table_dir, "lap3_lincs_mechanism_focused_candidates.csv")
)

top_overlap <- rbindlist(lapply(c(20L, 50L, 100L, 200L), function(top_n) {
  full_ids <- comparison[order(median_ncs_full)][seq_len(top_n), pert_id]
  no_lap3_ids <- comparison[order(median_ncs_no_lap3)][seq_len(top_n), pert_id]
  data.table(
    top_n = top_n,
    overlap_n = length(intersect(full_ids, no_lap3_ids)),
    overlap_fraction = length(intersect(full_ids, no_lap3_ids)) / top_n
  )
}))
fwrite(
  top_overlap,
  file.path(table_dir, "lap3_lincs_top_rank_overlap.csv")
)

summary_qc <- data.table(
  metric = c(
    "compounds_total",
    "stable_reverse_candidates",
    "top20_overlap",
    "top50_overlap",
    "top100_overlap",
    "top200_overlap"
  ),
  value = c(
    nrow(comparison),
    comparison[stable_reverse == TRUE, .N],
    top_overlap[top_n == 20L, overlap_fraction],
    top_overlap[top_n == 50L, overlap_fraction],
    top_overlap[top_n == 100L, overlap_fraction],
    top_overlap[top_n == 200L, overlap_fraction]
  )
)
fwrite(
  summary_qc,
  file.path(table_dir, "lap3_lincs_consensus_qc.csv")
)

cat(
  sprintf(
    "compounds=%d stable_reverse=%d top100_overlap=%.2f\n",
    nrow(comparison),
    comparison[stable_reverse == TRUE, .N],
    top_overlap[top_n == 100L, overlap_fraction]
  )
)
