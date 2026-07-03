#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})
data.table::setDTthreads(8)

root <- "/home/lzb/glioma"
result_dir <- file.path(root, "Data_scRNA_GEO/GBmap_Extended/results/Focused_Malignant_TAM_Communication")
table_dir <- file.path(result_dir, "tables")

paired <- fread(file.path(table_dir, "extended_focused_tam_donor_compartment_summary.csv"))
t2m <- fread(file.path(table_dir, "extended_focused_tam_axes_tam_to_malignant.csv"))
m2t <- fread(file.path(table_dir, "extended_focused_tam_axes_malignant_to_tam.csv"))
axis_results <- rbindlist(list(t2m, m2t), use.names = TRUE, fill = TRUE)

zscore <- function(x) {
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) {
    return(rep(NA_real_, length(x)))
  }
  (x - mean(x, na.rm = TRUE)) / s
}

safe_spearman <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 8L || length(unique(x[ok])) < 3L || length(unique(y[ok])) < 3L) {
    return(c(rho = NA_real_, p = NA_real_, n = sum(ok)))
  }
  test <- suppressWarnings(stats::cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  c(rho = unname(test$estimate), p = test$p.value, n = sum(ok))
}

axis_score_from_row <- function(row, dt) {
  direction <- row[["direction"]]
  ligand <- row[["ligand"]]
  receptors <- unlist(strsplit(row[["receptors_available"]], ";", fixed = TRUE))
  receptors <- receptors[nzchar(receptors)]
  if (direction == "tam_to_malignant") {
    sender_prefix <- "tam__"
    receiver_prefix <- "mal__"
  } else {
    sender_prefix <- "mal__"
    receiver_prefix <- "tam__"
  }
  ligand_col <- paste0(sender_prefix, ligand, "_expr")
  receptor_cols <- paste0(receiver_prefix, receptors, "_expr")
  receptor_cols <- intersect(receptor_cols, names(dt))
  if (!ligand_col %chin% names(dt) || length(receptor_cols) == 0L) {
    return(rep(NA_real_, nrow(dt)))
  }
  ligand_score <- dt[[ligand_col]]
  receptor_score <- rowMeans(as.matrix(dt[, ..receptor_cols]), na.rm = TRUE)
  (zscore(ligand_score) + zscore(receptor_score)) / 2
}

author_summary <- paired[, .(
  paired_author_donors = .N,
  malignant_cells = sum(mal__cells),
  tam_cells = sum(tam__cells),
  median_malignant_cells = median(mal__cells),
  median_tam_cells = median(tam__cells),
  mean_tam_fraction = mean(tam_fraction_within_pair)
), by = author]
author_summary[, donor_fraction := paired_author_donors / sum(paired_author_donors)]
setorder(author_summary, -paired_author_donors, author)
fwrite(author_summary, file.path(table_dir, "extended_focused_tam_author_retention_summary.csv"))

loa_rows <- list()
for (i in seq_len(nrow(axis_results))) {
  row <- axis_results[i]
  response <- row$response
  if (!response %chin% names(paired)) {
    next
  }
  axis_score <- axis_score_from_row(row, paired)
  for (drop_author in sort(unique(paired$author))) {
    keep <- paired$author != drop_author
    cor_result <- safe_spearman(axis_score[keep], paired[[response]][keep])
    loa_rows[[length(loa_rows) + 1L]] <- data.table(
      direction = row$direction,
      axis = row$axis,
      response = response,
      dropped_author = drop_author,
      n_author_donors = as.integer(cor_result[["n"]]),
      spearman_rho = cor_result[["rho"]],
      p_value = cor_result[["p"]]
    )
  }
}
loa <- rbindlist(loa_rows, use.names = TRUE, fill = TRUE)
fwrite(loa, file.path(table_dir, "extended_focused_tam_leave_one_author.csv"))

loa_summary <- loa[, .(
  leave_one_author_min_rho = min(spearman_rho, na.rm = TRUE),
  leave_one_author_max_rho = max(spearman_rho, na.rm = TRUE),
  leave_one_author_median_rho = median(spearman_rho, na.rm = TRUE),
  leave_one_author_sign_flips = sum(sign(spearman_rho) != sign(median(spearman_rho, na.rm = TRUE)), na.rm = TRUE),
  leave_one_author_max_p = max(p_value, na.rm = TRUE)
), by = .(direction, axis, response)]

axis_audit <- merge(
  axis_results,
  loa_summary,
  by = c("direction", "axis", "response"),
  all.x = TRUE
)
axis_audit[, robust_positive := (
  is.finite(spearman_rho) &
    spearman_rho > 0 &
    fdr_bh < 0.05 &
    ligand_detect_mean >= 0.05 &
    receptor_detect_mean >= 0.05 &
    leave_one_author_min_rho > 0
)]
axis_audit[, robust_negative := (
  is.finite(spearman_rho) &
    spearman_rho < 0 &
    fdr_bh < 0.05 &
    ligand_detect_mean >= 0.05 &
    receptor_detect_mean >= 0.05 &
    leave_one_author_max_rho < 0
)]
axis_audit[, nomination_tier := fifelse(
  robust_positive == TRUE & abs_rho >= 0.5, "Tier1_positive_candidate",
  fifelse(
    robust_positive == TRUE, "Tier2_positive_candidate",
    fifelse(
      robust_negative == TRUE, "negative_or_inverse_context",
      fifelse(fdr_bh < 0.05, "statistical_but_context_sensitive", "not_prioritized")
    )
  )
)]
setorder(axis_audit, direction, response, nomination_tier, fdr_bh, -abs_rho)
fwrite(axis_audit, file.path(table_dir, "extended_focused_tam_axis_closure_audit.csv"))

nominations <- axis_audit[
  nomination_tier %chin% c("Tier1_positive_candidate", "Tier2_positive_candidate"),
  .(
    direction, axis, ligand, receptors_available, response,
    n_author_donors, n_authors, ligand_detect_mean, receptor_detect_mean,
    spearman_rho, fdr_bh, leave_one_author_min_rho,
    leave_one_author_max_p, nomination_tier
  )
]
if (nrow(nominations) > 0L) {
  nominations[, abs_rho_nomination := abs(spearman_rho)]
  setorder(nominations, nomination_tier, fdr_bh, -abs_rho_nomination)
}
fwrite(nominations, file.path(table_dir, "extended_focused_tam_nomination_table.csv"))

top_author_fraction <- author_summary$donor_fraction[1]
verdict <- data.table(
  item = c(
    "analysis_status",
    "inference_units",
    "author_balance",
    "tier1_positive_candidates",
    "recommended_use",
    "main_limitation",
    "manuscript_boundary"
  ),
  verdict = c(
    "complete",
    paste0(nrow(paired), " paired author::donor units from ", uniqueN(paired$author), " authors"),
    paste0("largest author contributes ", sprintf("%.1f%%", 100 * top_author_fraction), " of paired units"),
    as.character(nrow(nominations[nomination_tier == "Tier1_positive_candidate"])),
    "Supplementary or Figure 6 reserve; use for focused malignant-TAM hypothesis nomination, not core Figure 1-5.",
    "Donor-level co-variation across integrated sc/snRNA studies lacks spatial adjacency and perturbational directionality.",
    "Use ecological compatibility / candidate axis language; avoid causal cell-cell communication claims."
  )
)
fwrite(verdict, file.path(table_dir, "extended_focused_tam_closure_verdict.csv"))

readme_append <- c(
  "",
  "## Closure Audit",
  "",
  paste0("Closure audit generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "- `tables/extended_focused_tam_author_retention_summary.csv`: author-level retention and balance.",
  "- `tables/extended_focused_tam_leave_one_author.csv`: leave-one-author correlations for tested axes.",
  "- `tables/extended_focused_tam_axis_closure_audit.csv`: all-axis robustness and nomination tiers.",
  "- `tables/extended_focused_tam_nomination_table.csv`: positive candidate axes after detection and leave-one-author filters.",
  "- `tables/extended_focused_tam_closure_verdict.csv`: project-facing decision summary.",
  "",
  "Recommended use: Supplementary or Figure 6 reserve. These results nominate focused malignant-TAM ecological axes but do not prove directional cell-cell communication."
)
write(readme_append, file.path(result_dir, "README.md"), append = TRUE)

message("Closure audit complete: ", result_dir)
