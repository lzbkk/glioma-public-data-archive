#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(openxlsx)
})

project_dir <- "/home/lzb/glioma"
setwd(project_dir)
data.table::setDTthreads(16)
set.seed(20260701)

out_dir <- "Data_Bulk_TCGA/Data_Merged/results/LAP3_State_CPTAC_GLASS_Projection"
table_dir <- file.path(out_dir, "tables")
source_dir <- file.path(out_dir, "source_data")
export_dir <- file.path(out_dir, "exports")
log_dir <- file.path(out_dir, "logs")
for (d in c(table_dir, source_dir, export_dir, log_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

log_file <- file.path(log_dir, "lap3_state_cptac_glass_projection.log")
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message")
on.exit({
  sink(type = "message")
  sink(type = "output")
  close(log_con)
}, add = TRUE)

log_msg <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n")
write_table <- function(x, filename) fwrite(as.data.table(x), file.path(table_dir, filename))
`%||%` <- function(x, y) if (is.null(x)) y else x

collapse_expr <- function(mat) {
  mat <- as.matrix(mat)
  storage.mode(mat) <- "double"
  rn <- rownames(mat)
  if (is.null(rn)) stop("Matrix rownames are required")
  if (!anyDuplicated(rn)) return(mat)
  dt <- as.data.table(mat)
  dt[, gene := rn]
  dt <- dt[, lapply(.SD, mean, na.rm = TRUE), by = gene]
  out <- as.matrix(dt[, -"gene"])
  rownames(out) <- dt$gene
  storage.mode(out) <- "double"
  out
}

module_score <- function(expr, genes, min_genes = 3L) {
  genes <- intersect(unique(genes), rownames(expr))
  if (length(genes) < min_genes) return(rep(NA_real_, ncol(expr)))
  z <- t(scale(t(expr[genes, , drop = FALSE])))
  z[!is.finite(z)] <- NA_real_
  colMeans(z, na.rm = TRUE)
}

spearman_safe <- function(x, y, min_n = 10L) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) < min_n || length(unique(x)) < 3L || length(unique(y)) < 3L) {
    return(list(n = length(x), rho = NA_real_, p_value = NA_real_))
  }
  test <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  list(n = length(x), rho = unname(test$estimate), p_value = test$p.value)
}

paired_wilcox <- function(a, b, min_n = 10L) {
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < min_n) return(list(n = sum(ok), p_value = NA_real_))
  wt <- suppressWarnings(wilcox.test(b[ok], a[ok], paired = TRUE, exact = FALSE))
  list(n = sum(ok), p_value = wt$p.value)
}

coverage_table <- function(dataset, assay, feature_genes, gene_sets, row_genes) {
  rbindlist(lapply(names(gene_sets), function(nm) {
    genes <- unique(gene_sets[[nm]])
    present <- intersect(genes, row_genes)
    data.table(
      dataset = dataset,
      assay = assay,
      module = nm,
      requested = length(genes),
      available = length(present),
      coverage = length(present) / length(genes),
      missing_genes = paste(sort(setdiff(genes, row_genes)), collapse = ";")
    )
  }))
}

score_matrix <- function(expr, prefix, gene_sets) {
  out <- data.table(sample_id = colnames(expr))
  for (nm in names(gene_sets)) {
    out[, (paste0(prefix, nm)) := module_score(expr, gene_sets[[nm]])]
  }
  out
}

log_msg("Loading LAP3 state submodule sets")
submodule_export <- readRDS("Data_Bulk_TCGA/Data_Merged/results/LAP3_State_Submodules/exports/lap3_state_submodules_projection.rds")
submodule_sets <- submodule_export$submodule_sets
gene_sets <- c(list(LAP3_STATE_UNION = unique(submodule_export$gene_assignment$gene)), submodule_sets)
genes_needed <- sort(unique(c("LAP3", unlist(gene_sets, use.names = FALSE))))

log_msg("CPTAC: reading Table S2 transcriptome and proteome rows")
table_s2 <- "Data_Protein_Public/data/CPTAC_GBM_Publication/Wang2021_CancerCell_Table_S2_mmc3.xlsx"
stopifnot(file.exists(table_s2))

cptac_mrna_dt <- as.data.table(read.xlsx(table_s2, sheet = "gene_expression_fpkm_uq", check.names = FALSE))
stopifnot(all(c("gene_name", "gene_type") %in% names(cptac_mrna_dt)))
cptac_mrna_dt <- cptac_mrna_dt[gene_name %chin% genes_needed]
cptac_mrna_sample_cols <- grep("^C3[NL]-", names(cptac_mrna_dt), value = TRUE)
cptac_mrna <- as.matrix(cptac_mrna_dt[, ..cptac_mrna_sample_cols])
rownames(cptac_mrna) <- cptac_mrna_dt$gene_name
cptac_mrna <- collapse_expr(cptac_mrna)
cptac_mrna_log <- log2(cptac_mrna + 1)

cptac_prot_dt <- as.data.table(read.xlsx(table_s2, sheet = "proteome_normalized", check.names = FALSE))
stopifnot("symbol" %in% names(cptac_prot_dt))
cptac_prot_dt <- cptac_prot_dt[symbol %chin% genes_needed]
cptac_prot_sample_cols <- grep("^C3[NL]-", names(cptac_prot_dt), value = TRUE)
cptac_prot <- as.matrix(cptac_prot_dt[, ..cptac_prot_sample_cols])
rownames(cptac_prot) <- cptac_prot_dt$symbol
cptac_prot <- collapse_expr(cptac_prot)

cptac_mrna_scores <- score_matrix(cptac_mrna_log, "mrna_", gene_sets)
cptac_prot_scores <- score_matrix(cptac_prot, "protein_", gene_sets)

cptac_samples <- Reduce(intersect, list(cptac_mrna_scores$sample_id, cptac_prot_scores$sample_id))
cptac_projection <- merge(
  cptac_mrna_scores[sample_id %chin% cptac_samples],
  cptac_prot_scores[sample_id %chin% cptac_samples],
  by = "sample_id",
  all = TRUE
)
cptac_projection[, lap3_mrna_log2 := as.numeric(cptac_mrna_log["LAP3", sample_id])]
cptac_projection[, lap3_protein := as.numeric(cptac_prot["LAP3", sample_id])]
setcolorder(cptac_projection, c("sample_id", "lap3_mrna_log2", "lap3_protein"))

log_msg("CPTAC: adding direct BCAA and phosphosite readouts")
bcaa_file <- "Data_Protein_Public/results/LAP3_CPTAC_BCAA_Metabolomics/tables/cptac_lap3_bcaa_matched_analysis_data.csv"
if (file.exists(bcaa_file)) {
  bcaa <- fread(bcaa_file)
  bcaa_keep <- intersect(
    c("sampleId", "l_leucine", "dl_isoleucine", "l_valine", "bcaa_composite", "idh_mutation_status"),
    names(bcaa)
  )
  bcaa <- bcaa[, ..bcaa_keep]
  setnames(bcaa, "sampleId", "sample_id")
  cptac_projection <- merge(cptac_projection, bcaa, by = "sample_id", all.x = TRUE)
}

phospho_file <- "Data_Protein_Public/results/LAP3_Protein_Readout/tables/cptac_gbm_target_phosphosite_values.csv"
if (file.exists(phospho_file)) {
  phospho <- fread(phospho_file)
  phospho <- phospho[gene_symbol %chin% c("RPS6", "RPS6KB1", "EIF4EBP1", "MTOR")]
  phospho[, z_value := as.numeric(scale(value)), by = stableId]
  phospho_score <- phospho[, .(
    phospho_mtorc1_target_score = mean(z_value, na.rm = TRUE),
    phospho_n_sites = uniqueN(stableId),
    phospho_n_values = sum(is.finite(value))
  ), by = .(sample_id = sampleId)]
  cptac_projection <- merge(cptac_projection, phospho_score, by = "sample_id", all.x = TRUE)
}

if (!"idh_mutation_status" %in% names(cptac_projection)) {
  cptac_projection[, idh_mutation_status := NA_character_]
}
idh_mut_file <- "Data_Protein_Public/results/LAP3_CPTAC_BCAA_Metabolomics/tables/cptac_idh1_idh2_mutation_samples.csv"
if (file.exists(idh_mut_file)) {
  idh_mut <- fread(idh_mut_file)
  cptac_projection[, idh_mutation_status := fifelse(
    sample_id %chin% idh_mut$sampleId,
    "IDH1_or_IDH2_mutant",
    "no_IDH1_or_IDH2_mutation_detected"
  )]
}

coverage <- rbindlist(list(
  coverage_table("CPTAC_GBM", "mRNA_FPKM_UQ_log2", genes_needed, gene_sets, rownames(cptac_mrna_log)),
  coverage_table("CPTAC_GBM", "total_protein_normalized", genes_needed, gene_sets, rownames(cptac_prot))
), fill = TRUE)

cptac_exposures <- c(
  "lap3_mrna_log2", "lap3_protein",
  paste0("mrna_", names(gene_sets)),
  paste0("protein_", names(gene_sets))
)
cptac_exposures <- intersect(cptac_exposures, names(cptac_projection))
cptac_outcomes <- intersect(c(
  "lap3_protein", "l_leucine", "dl_isoleucine", "l_valine",
  "bcaa_composite", "phospho_mtorc1_target_score"
), names(cptac_projection))
cptac_strata <- list(
  All = cptac_projection,
  `IDH1/2-wildtype sensitivity` = cptac_projection[
    idh_mutation_status == "no_IDH1_or_IDH2_mutation_detected"
  ]
)
cptac_cor <- rbindlist(lapply(names(cptac_strata), function(stratum) {
  dat <- cptac_strata[[stratum]]
  rbindlist(lapply(cptac_exposures, function(exposure) {
    rbindlist(lapply(setdiff(cptac_outcomes, exposure), function(outcome) {
      res <- spearman_safe(dat[[exposure]], dat[[outcome]], min_n = 20L)
      data.table(
        dataset = "CPTAC_GBM",
        stratum = stratum,
        exposure = exposure,
        outcome = outcome,
        n_complete = res$n,
        spearman_rho = res$rho,
        p_value = res$p_value
      )
    }))
  }))
}), fill = TRUE)
cptac_cor[, fdr := p.adjust(p_value, method = "BH"), by = .(stratum, outcome)]

log_msg("GLASS: reading TPM matrix and scoring submodules")
glass_expr_file <- "Data_Longitudinal_Public/GLASS/data/expression/gene_tpm_matrix_all_samples.tsv"
glass_pair_file <- "Data_Longitudinal_Public/GLASS/results/LAP3_Longitudinal_Feasibility/tables/glass_strict_primary_recurrence_pairs.csv"
stopifnot(file.exists(glass_expr_file), file.exists(glass_pair_file))
glass_expr <- fread(glass_expr_file, showProgress = FALSE)
stopifnot(names(glass_expr)[1] == "Gene_symbol")
glass_expr <- glass_expr[Gene_symbol %chin% genes_needed]
glass_genes <- glass_expr$Gene_symbol
glass_expr[, Gene_symbol := NULL]
glass_mat <- as.matrix(glass_expr)
storage.mode(glass_mat) <- "double"
rownames(glass_mat) <- glass_genes
glass_mat <- collapse_expr(glass_mat)
colnames(glass_mat) <- gsub("\\.", "-", colnames(glass_mat))
glass_log <- log2(glass_mat + 1)

coverage <- rbindlist(list(
  coverage,
  coverage_table("GLASS", "mRNA_TPM_log2", genes_needed, gene_sets, rownames(glass_log))
), fill = TRUE)

glass_scores <- score_matrix(glass_log, "", gene_sets)
setnames(glass_scores, "sample_id", "tumor_barcode")
glass_scores[, log2_lap3 := as.numeric(glass_log["LAP3", tumor_barcode])]
setcolorder(glass_scores, c("tumor_barcode", "log2_lap3"))

pairs <- fread(glass_pair_file)
score_cols <- names(gene_sets)
score_a <- copy(glass_scores)
setnames(score_a, names(score_a), c("tumor_barcode_a", "log2_lap3_a_scored", paste0(score_cols, "_a")))
score_b <- copy(glass_scores)
setnames(score_b, names(score_b), c("tumor_barcode_b", "log2_lap3_b_scored", paste0(score_cols, "_b")))
pairs_scored <- merge(pairs, score_a, by = "tumor_barcode_a", all.x = TRUE)
pairs_scored <- merge(pairs_scored, score_b, by = "tumor_barcode_b", all.x = TRUE)
stopifnot(nrow(pairs_scored) == nrow(pairs))
for (nm in score_cols) {
  pairs_scored[, (paste0("delta_", nm)) := get(paste0(nm, "_b")) - get(paste0(nm, "_a"))]
}
pairs_scored[, delta_log2_lap3_scored := log2_lap3_b_scored - log2_lap3_a_scored]

glass_strata <- c("All", sort(unique(na.omit(pairs_scored$strict_stratum))))
glass_tests <- rbindlist(lapply(glass_strata, function(stratum) {
  dat <- if (stratum == "All") pairs_scored else pairs_scored[strict_stratum == stratum]
  rbindlist(c(
    lapply(score_cols, function(score) {
      a <- dat[[paste0(score, "_a")]]
      b <- dat[[paste0(score, "_b")]]
      wt <- paired_wilcox(a, b, min_n = 10L)
      data.table(
        dataset = "GLASS",
        analysis = "paired_recurrence_change",
        stratum = stratum,
        module = score,
        n_pairs = wt$n,
        median_primary = median(a, na.rm = TRUE),
        median_recurrence = median(b, na.rm = TRUE),
        median_delta = median(b - a, na.rm = TRUE),
        estimate = NA_real_,
        p_value = wt$p_value
      )
    }),
    lapply(score_cols, function(score) {
      res <- spearman_safe(dat$delta_log2_lap3_scored, dat[[paste0("delta_", score)]], min_n = 10L)
      data.table(
        dataset = "GLASS",
        analysis = "delta_lap3_delta_module",
        stratum = stratum,
        module = score,
        n_pairs = res$n,
        median_primary = NA_real_,
        median_recurrence = NA_real_,
        median_delta = median(dat[[paste0("delta_", score)]], na.rm = TRUE),
        estimate = res$rho,
        p_value = res$p_value
      )
    })
  ), fill = TRUE)
}), fill = TRUE)
glass_tests[, fdr := p.adjust(p_value, method = "BH"), by = .(analysis, stratum)]

log_msg("Writing outputs")
write_table(coverage, "projection_gene_coverage.csv")
write_table(cptac_projection, "cptac_lap3_state_submodule_projection.csv")
write_table(cptac_cor, "cptac_lap3_state_submodule_correlations.csv")
write_table(glass_scores, "glass_lap3_state_submodule_scores_all_samples.csv")
write_table(pairs_scored, "glass_lap3_state_submodule_strict_pair_scores.csv")
write_table(glass_tests, "glass_lap3_state_submodule_paired_tests.csv")

primary_summary <- rbindlist(list(
  cptac_cor[
    stratum == "All" &
      exposure %chin% c("mrna_LAP3_STATE_UNION", "protein_LAP3_STATE_UNION",
                        paste0("mrna_", names(submodule_sets)), paste0("protein_", names(submodule_sets))) &
      outcome %chin% c("lap3_protein", "bcaa_composite", "phospho_mtorc1_target_score")
  ][order(outcome, fdr, -abs(spearman_rho))][
    , .(dataset, analysis = "cross_sectional_correlation", stratum,
        feature = exposure, outcome, n = n_complete, effect = spearman_rho,
        p_value, fdr)
  ],
  glass_tests[
    stratum %chin% c("All", "IDH-wildtype GBM") &
      module %chin% names(gene_sets)
  ][order(analysis, stratum, fdr, module)][
    , .(dataset, analysis, stratum, feature = module, outcome = analysis,
        n = n_pairs, effect = fifelse(analysis == "paired_recurrence_change", median_delta, estimate),
        p_value, fdr)
  ]
), fill = TRUE)
write_table(primary_summary, "lap3_state_cptac_glass_primary_summary.csv")

saveRDS(
  list(
    generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    gene_sets = gene_sets,
    coverage = coverage,
    cptac_projection = cptac_projection,
    cptac_correlations = cptac_cor,
    glass_scores = glass_scores,
    glass_pairs = pairs_scored,
    glass_tests = glass_tests,
    primary_summary = primary_summary
  ),
  file.path(export_dir, "lap3_state_cptac_glass_projection.rds")
)

fmt_top <- function(dt, n = 10L) {
  if (!nrow(dt)) return("No rows.")
  paste(capture.output(print(head(dt, n))), collapse = "\n")
}
cptac_top <- cptac_cor[
  stratum == "All" &
    outcome %chin% c("lap3_protein", "bcaa_composite", "phospho_mtorc1_target_score") &
    grepl("LAP3_STATE_UNION|LAP3_MALIGNANT_STATE_MODULE|LAP3_MYELOID_TAM_CONTEXT_MODULE|LAP3_ANABOLIC_TRANSLATION_MODULE|LAP3_PROTEOSTASIS_STRESS_MODULE|LAP3_HYPOXIA_PERINECROTIC_MODULE", exposure)
][order(outcome, fdr, -abs(spearman_rho))]
glass_top <- glass_tests[
  stratum %chin% c("All", "IDH-wildtype GBM")
][order(analysis, stratum, fdr, module)]

readme <- c(
  "# LAP3 State CPTAC/GLASS Projection",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Purpose",
  "",
  "Project the frozen `LAP3_STATE_UNION` and five interpretable submodules into CPTAC GBM multi-omics and GLASS longitudinal RNA-seq.",
  "",
  "## Inputs",
  "",
  "- Submodules: `Data_Bulk_TCGA/Data_Merged/results/LAP3_State_Submodules/exports/lap3_state_submodules_projection.rds`",
  "- CPTAC: Wang et al. Table S2 transcriptome/proteome/metabolome plus existing phosphosite tables.",
  "- GLASS: current-release gene TPM matrix and strict primary-to-recurrence patient pairs.",
  "",
  "## Design",
  "",
  "- CPTAC inference unit: tumor sample/patient; analyses are cross-sectional correlations.",
  "- GLASS inference unit: patient pair; analyses are paired recurrence change and within-patient delta correlations.",
  "- Module score: mean of gene-wise z-scores using available genes in each assay.",
  "- Multiplicity: BH correction within CPTAC stratum/outcome and GLASS analysis/stratum.",
  "",
  "## Key Output Tables",
  "",
  "- `tables/projection_gene_coverage.csv`",
  "- `tables/cptac_lap3_state_submodule_projection.csv`",
  "- `tables/cptac_lap3_state_submodule_correlations.csv`",
  "- `tables/glass_lap3_state_submodule_strict_pair_scores.csv`",
  "- `tables/glass_lap3_state_submodule_paired_tests.csv`",
  "- `tables/lap3_state_cptac_glass_primary_summary.csv`",
  "- `exports/lap3_state_cptac_glass_projection.rds`",
  "",
  "## Main Result",
  "",
  "- CPTAC direct BCAA composite does not support a submodule association, preserving the direct metabolite negative boundary.",
  "- CPTAC LAP3 protein is strongly associated with protein/mRNA state components, especially malignant-state, myeloid/TAM context and the union score.",
  "- CPTAC phospho mTORC1 target aggregate shows cross-sectional state/phospho-readout associations, but this is not causal LAP3 phosphorylation evidence.",
  "- GLASS strict pairs show strong coordinated delta relationships between LAP3 and malignant-state/union/myeloid/proteostasis modules, without universal recurrence activation.",
  "",
  "## CPTAC Snapshot",
  "",
  fmt_top(cptac_top, 12L),
  "",
  "## GLASS Snapshot",
  "",
  fmt_top(glass_top, 16L),
  "",
  "## Interpretation Boundary",
  "",
  "CPTAC can test concordance with total protein, phosphosite readouts and direct BCAA metabolites, but cannot establish LAP3 enzymatic causality. GLASS tests coordinated longitudinal state change, not malignant-cell-intrinsic mechanism."
)
writeLines(readme, file.path(out_dir, "README.md"))

log_msg("Completed LAP3 state CPTAC/GLASS projection")
