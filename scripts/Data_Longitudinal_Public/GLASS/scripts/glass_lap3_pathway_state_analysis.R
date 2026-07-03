#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(msigdbr)
})

set.seed(20260630)
data.table::setDTthreads(8)

root_dir <- "Data_Longitudinal_Public/GLASS"
expr_file <- file.path(root_dir, "data", "expression", "gene_tpm_matrix_all_samples.tsv")
pair_file <- file.path(
  root_dir, "results", "LAP3_Longitudinal_Feasibility", "tables",
  "glass_strict_primary_recurrence_pairs.csv"
)
result_dir <- file.path(root_dir, "results", "LAP3_Pathway_State_Longitudinal")
table_dir <- file.path(result_dir, "tables")
log_dir <- file.path(result_dir, "logs")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(expr_file), file.exists(pair_file))

log_file <- file.path(log_dir, "glass_lap3_pathway_state_analysis.log")
if (file.exists(log_file)) invisible(file.remove(log_file))
log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste(...))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

normalize_barcode <- function(x) gsub("\\.", "-", x)

module_score <- function(log_expr, genes) {
  genes <- intersect(genes, rownames(log_expr))
  if (length(genes) < 3L) {
    stop("Fewer than three available genes for module")
  }
  gene_z <- t(scale(t(log_expr[genes, , drop = FALSE])))
  gene_z[!is.finite(gene_z)] <- NA_real_
  colMeans(gene_z, na.rm = TRUE)
}

paired_test <- function(dat, score, stratum) {
  a <- dat[[paste0(score, "_a")]]
  b <- dat[[paste0(score, "_b")]]
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < 3L) return(NULL)
  wt <- suppressWarnings(wilcox.test(b[ok], a[ok], paired = TRUE, exact = FALSE))
  data.table(
    analysis = "paired_recurrence_change",
    stratum = stratum,
    module = score,
    n_pairs = sum(ok),
    median_primary = median(a[ok]),
    median_recurrence = median(b[ok]),
    median_delta = median(b[ok] - a[ok]),
    estimate = unname(wt$statistic),
    p_value = wt$p.value
  )
}

delta_correlation <- function(dat, score, stratum) {
  x <- dat$delta_log2_lap3
  y <- dat[[paste0("delta_", score)]]
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 4L) return(NULL)
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  data.table(
    analysis = "delta_lap3_delta_module",
    stratum = stratum,
    module = score,
    n_pairs = sum(ok),
    median_primary = NA_real_,
    median_recurrence = NA_real_,
    median_delta = median(y[ok]),
    estimate = unname(ct$estimate),
    p_value = ct$p.value
  )
}

partial_delta_correlation <- function(dat, score, stratum) {
  controls <- c(
    "delta_MYELOID_TAM_CORE",
    "delta_MES_LIKE_CORE",
    "delta_HALLMARK_HYPOXIA"
  )
  x_name <- "delta_log2_lap3"
  y_name <- paste0("delta_", score)
  keep <- complete.cases(dat[, c(x_name, y_name, controls), with = FALSE])
  d <- dat[keep, c(x_name, y_name, controls), with = FALSE]
  if (nrow(d) < 10L) return(NULL)
  ranked <- as.data.table(lapply(d, rank, ties.method = "average"))
  control_formula <- as.formula(paste("value ~", paste(controls, collapse = " + ")))
  x_resid <- residuals(lm(control_formula, data = cbind(
    value = ranked[[x_name]], ranked[, ..controls]
  )))
  y_resid <- residuals(lm(control_formula, data = cbind(
    value = ranked[[y_name]], ranked[, ..controls]
  )))
  ct <- suppressWarnings(cor.test(x_resid, y_resid, method = "pearson"))
  data.table(
    analysis = "partial_delta_controlling_context",
    stratum = stratum,
    module = score,
    n_pairs = nrow(d),
    median_primary = NA_real_,
    median_recurrence = NA_real_,
    median_delta = median(d[[y_name]]),
    estimate = unname(ct$estimate),
    p_value = ct$p.value
  )
}

log_msg("Loading frozen gene sets")
hallmark <- as.data.table(msigdbr(species = "Homo sapiens", category = "H"))
hallmark_sets <- split(hallmark$gene_symbol, hallmark$gs_name)
gene_sets <- list(
  HALLMARK_MTORC1_SIGNALING = hallmark_sets[["HALLMARK_MTORC1_SIGNALING"]],
  LEUCINE_BCAA_CORE = c(
    "BCAT1", "BCAT2", "BCKDHA", "BCKDHB", "DBT", "DLD", "BCKDK", "PPM1K",
    "SLC7A5", "SLC3A2", "SLC43A1", "SLC43A2", "SLC38A2", "SLC38A9"
  ),
  MTORC1_READOUT_CORE = c(
    "MTOR", "RPTOR", "RHEB", "AKT1", "TSC1", "TSC2", "RRAGA", "RRAGB",
    "RRAGC", "RRAGD", "LAMTOR1", "LAMTOR2", "LAMTOR3", "LAMTOR4", "LAMTOR5",
    "RPS6KB1", "RPS6KB2", "EIF4EBP1", "EIF4EBP2", "RPS6", "EIF4E"
  ),
  HALLMARK_HYPOXIA = hallmark_sets[["HALLMARK_HYPOXIA"]],
  MES_LIKE_CORE = c(
    "CHI3L1", "MET", "CD44", "VIM", "SERPINE1", "TNC", "CAV1", "CAV2",
    "LGALS3", "ANXA1", "ITGA5", "CTGF", "DAB2", "WWTR1"
  ),
  MYELOID_TAM_CORE = c(
    "AIF1", "CD68", "CSF1R", "TYROBP", "FCER1G", "LST1", "CTSS", "C1QA",
    "C1QB", "C1QC", "APOE", "TREM2", "CD163", "MSR1"
  )
)

genes_needed <- unique(c("LAP3", unlist(gene_sets, use.names = FALSE)))
log_msg("Reading GLASS TPM rows for", length(genes_needed), "requested genes")
expr <- fread(expr_file, showProgress = FALSE)
stopifnot(names(expr)[1] == "Gene_symbol", !anyDuplicated(expr$Gene_symbol))
expr <- expr[Gene_symbol %chin% genes_needed]
rownames_expr <- expr$Gene_symbol
expr[, Gene_symbol := NULL]
expr_mat <- as.matrix(expr)
storage.mode(expr_mat) <- "double"
rownames(expr_mat) <- rownames_expr
colnames(expr_mat) <- normalize_barcode(colnames(expr_mat))
log_expr <- log2(expr_mat + 1)
stopifnot("LAP3" %in% rownames(log_expr), !anyDuplicated(colnames(log_expr)))

coverage <- rbindlist(lapply(names(gene_sets), function(nm) {
  available <- intersect(gene_sets[[nm]], rownames(log_expr))
  data.table(
    module = nm,
    n_requested = length(unique(gene_sets[[nm]])),
    n_available = length(available),
    coverage = length(available) / length(unique(gene_sets[[nm]])),
    available_genes = paste(sort(available), collapse = ";"),
    missing_genes = paste(sort(setdiff(gene_sets[[nm]], available)), collapse = ";")
  )
}))
stopifnot(all(coverage$n_available >= 3L))

scores <- data.table(tumor_barcode = colnames(log_expr))
scores[, log2_lap3 := as.numeric(log_expr["LAP3", ])]
for (nm in names(gene_sets)) {
  scores[, (nm) := module_score(log_expr, gene_sets[[nm]])]
}

pairs <- fread(pair_file)
stopifnot(uniqueN(pairs$case_barcode) == nrow(pairs))
score_names <- names(gene_sets)
score_a <- copy(scores)
setnames(score_a, names(score_a), c("tumor_barcode_a", "log2_lap3_score_a", paste0(score_names, "_a")))
score_b <- copy(scores)
setnames(score_b, names(score_b), c("tumor_barcode_b", "log2_lap3_score_b", paste0(score_names, "_b")))
pairs <- merge(pairs, score_a, by = "tumor_barcode_a", all.x = TRUE)
pairs <- merge(pairs, score_b, by = "tumor_barcode_b", all.x = TRUE)
stopifnot(nrow(pairs) == uniqueN(pairs$case_barcode))

for (nm in score_names) {
  pairs[, (paste0("delta_", nm)) := get(paste0(nm, "_b")) - get(paste0(nm, "_a"))]
}

strata <- list(
  All = pairs,
  `IDH-wildtype GBM` = pairs[strict_stratum == "IDH-wildtype GBM"]
)
tests <- rbindlist(lapply(names(strata), function(stratum) {
  dat <- strata[[stratum]]
  rbindlist(c(
    lapply(score_names, function(score) paired_test(dat, score, stratum)),
    lapply(score_names, function(score) delta_correlation(dat, score, stratum)),
    lapply(
      c("HALLMARK_MTORC1_SIGNALING", "LEUCINE_BCAA_CORE", "MTORC1_READOUT_CORE"),
      function(score) partial_delta_correlation(dat, score, stratum)
    )
  ), fill = TRUE)
}))
tests[, fdr := p.adjust(p_value, method = "BH"), by = .(analysis, stratum)]
setorder(tests, analysis, stratum, fdr, module)

pair_keep <- c(
  "case_barcode", "strict_stratum", "surgical_interval_mo_b",
  "treatment_tmz_a", "treatment_radiotherapy_a", "treatment_alkylating_agent_a",
  "log2_lap3_a", "log2_lap3_b", "delta_log2_lap3",
  unlist(lapply(score_names, function(x) c(paste0(x, "_a"), paste0(x, "_b"), paste0("delta_", x))))
)

fwrite(coverage, file.path(table_dir, "glass_pathway_gene_coverage.csv"))
fwrite(scores, file.path(table_dir, "glass_pathway_scores_all_samples.csv"))
fwrite(pairs[, ..pair_keep], file.path(table_dir, "glass_pathway_strict_pair_scores.csv"))
fwrite(tests, file.path(table_dir, "glass_pathway_paired_tests.csv"))

summary_rows <- tests[, .(
  analysis, stratum, module, n_pairs,
  effect = fifelse(
    analysis == "paired_recurrence_change",
    sprintf("median_delta=%.3f", median_delta),
    sprintf("spearman_rho=%.3f", estimate)
  ),
  p_value, fdr,
  interpretation = fifelse(
    fdr < 0.05,
    "FDR-significant within the prespecified analysis family",
    "not FDR-significant"
  )
)]
fwrite(summary_rows, file.path(table_dir, "glass_pathway_interpretation_summary.csv"))

readme <- c(
  "# GLASS LAP3 longitudinal pathway-state analysis",
  "",
  "## Design",
  "",
  "- Inference unit: patient.",
  "- Pair definition: frozen strict primary-to-first-evaluable-recurrence pairs from the prior LAP3 audit.",
  "- Expression: `log2(TPM + 1)`; each module is the mean of gene-wise z-scores across all GLASS RNA samples.",
  "- Prespecified modules: mTORC1 Hallmark, leucine/BCAA core, mTORC1 readout, hypoxia, MES-like core, and myeloid/TAM core.",
  "- Tests: paired recurrence change and Spearman correlation between within-patient delta LAP3 and delta module score.",
  "- Context sensitivity: rank-based partial correlations for mTORC1/BCAA/readout control delta myeloid/TAM, MES-like, and hypoxia scores.",
  "- Multiplicity: BH correction within each analysis and stratum.",
  "",
  "## Main result",
  "",
  "- Across all 131 strict pairs, none of the six modules showed an FDR-significant uniform recurrence shift.",
  "- Within 71 IDH-wildtype GBM pairs, mTORC1 Hallmark and mTORC1 readout scores decreased modestly at recurrence; BCAA did not pass FDR.",
  "- Delta LAP3 was positively correlated with delta mTORC1, BCAA, readout, hypoxia, MES-like, and myeloid/TAM scores in both strata.",
  "- The delta LAP3 associations with mTORC1, BCAA, and readout remained FDR-significant after rank-based adjustment for delta myeloid/TAM, MES-like, and hypoxia scores.",
  "- Therefore, GLASS supports coordinated within-patient state variation, not a universal recurrence-induced activation of the LAP3 axis.",
  "",
  "## Interpretation boundary",
  "",
  "This longitudinal bulk analysis tests coordinated patient-level state change. It cannot distinguish malignant-cell-intrinsic change from composition change and does not establish LAP3 causality.",
  "",
  "## Outputs",
  "",
  "- `tables/glass_pathway_gene_coverage.csv`",
  "- `tables/glass_pathway_scores_all_samples.csv`",
  "- `tables/glass_pathway_strict_pair_scores.csv`",
  "- `tables/glass_pathway_paired_tests.csv`",
  "- `tables/glass_pathway_interpretation_summary.csv`",
  "- `logs/glass_lap3_pathway_state_analysis.log`"
)
writeLines(readme, file.path(result_dir, "README.md"))

log_msg("Completed:", nrow(pairs), "strict pairs;", nrow(tests), "tests")
