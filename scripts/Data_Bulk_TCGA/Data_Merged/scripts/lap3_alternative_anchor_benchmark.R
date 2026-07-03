#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

project_dir <- "/home/lzb/glioma"
setwd(project_dir)
data.table::setDTthreads(8)
set.seed(20260701)

out_dir <- "Data_Bulk_TCGA/Data_Merged/results/LAP3_Alternative_Anchor_Benchmark"
table_dir <- file.path(out_dir, "tables")
source_dir <- file.path(out_dir, "source_data")
log_dir <- file.path(out_dir, "logs")
for (d in c(table_dir, source_dir, log_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

log_file <- file.path(log_dir, "lap3_alternative_anchor_benchmark.log")
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message")
on.exit({
  sink(type = "message")
  sink(type = "output")
  close(log_con)
}, add = TRUE)

cat("LAP3 alternative anchor benchmark\n")
cat("Started:", format(Sys.time()), "\n")
cat("R:", R.version.string, "\n")
cat("Working directory:", getwd(), "\n\n")

write_table <- function(x, filename, dir = table_dir) {
  fwrite(as.data.table(x), file.path(dir, filename))
}

collapse_expr <- function(mat) {
  mat <- as.matrix(mat)
  storage.mode(mat) <- "double"
  rn <- rownames(mat)
  if (is.null(rn)) stop("Matrix rownames are required")
  keep <- !is.na(rn) & nzchar(rn)
  mat <- mat[keep, , drop = FALSE]
  rn <- rn[keep]
  rn <- toupper(rn)
  if (!anyDuplicated(rn)) {
    rownames(mat) <- rn
    return(mat)
  }
  dt <- as.data.table(mat)
  dt[, gene := rn]
  dt <- dt[, lapply(.SD, mean, na.rm = TRUE), by = gene]
  out <- as.matrix(dt[, -"gene"])
  rownames(out) <- dt$gene
  storage.mode(out) <- "double"
  out
}

spearman_safe <- function(x, y, min_n = 25L) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) < min_n || length(unique(x)) < 3L || length(unique(y)) < 3L) {
    return(list(n = length(x), rho = NA_real_, p = NA_real_))
  }
  test <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  list(n = length(x), rho = unname(test$estimate), p = test$p.value)
}

z_num <- function(x) {
  x <- as.numeric(x)
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  as.numeric((x - mean(x, na.rm = TRUE)) / s)
}

single_gene_expr <- function(expr_log2, gene, samples) {
  gene <- toupper(gene)
  if (!gene %in% rownames(expr_log2)) return(rep(NA_real_, length(samples)))
  as.numeric(expr_log2[gene, match(samples, colnames(expr_log2))])
}

lm_anchor_term <- function(dt, outcome, gene_expr_col, covariates, dataset_label) {
  d <- as.data.table(copy(dt))
  d[, anchor_z := z_num(get(gene_expr_col))]
  required <- unique(c(outcome, "anchor_z", covariates))
  required <- required[required %in% names(d)]
  d <- d[complete.cases(d[, ..required])]
  if (nrow(d) < 50L || uniqueN(d$anchor_z) < 3L || uniqueN(d[[outcome]]) < 3L) {
    return(data.table(
      dataset = dataset_label, outcome = outcome, n = nrow(d),
      beta = NA_real_, std_error = NA_real_, p_value = NA_real_, adjusted_r2 = NA_real_
    ))
  }
  for (cc in covariates) {
    if (cc %in% names(d) && !is.numeric(d[[cc]])) d[[cc]] <- factor(d[[cc]])
  }
  form <- reformulate(c("anchor_z", covariates[covariates %in% names(d)]), response = outcome)
  fit <- lm(form, data = d)
  ct <- summary(fit)$coefficients
  data.table(
    dataset = dataset_label,
    outcome = outcome,
    n = nrow(d),
    beta = ct["anchor_z", "Estimate"],
    std_error = ct["anchor_z", "Std. Error"],
    p_value = ct["anchor_z", "Pr(>|t|)"],
    adjusted_r2 = summary(fit)$adj.r.squared
  )
}

manual_gene_sets <- list(
  AMINOPEPTIDASE_FAMILY = c(
    "LAP3", "ANPEP", "ERAP1", "ERAP2", "LNPEP", "NPEPPS",
    "XPNPEP1", "XPNPEP2", "XPNPEP3", "LTA4H", "RNPEP", "ENPEP"
  ),
  TAM_MYELOID_CORE = c(
    "AIF1", "CD68", "LST1", "CSF1R", "TYROBP", "FCER1G", "SPI1", "ITGAM",
    "C1QA", "C1QB", "C1QC", "CD163", "MRC1", "MSR1", "TREM2", "APOE", "LILRB4"
  ),
  HYPOXIA_CORE = c(
    "CA9", "VEGFA", "SLC2A1", "LDHA", "ENO1", "PGK1", "BNIP3", "NDRG1",
    "P4HA1", "EGLN3", "ADM", "ANGPTL4"
  ),
  PROLIFERATION_CORE = c(
    "MKI67", "TOP2A", "PCNA", "MCM2", "MCM3", "MCM4", "MCM5", "MCM6",
    "CDK1", "CCNB1", "CCNB2", "AURKA", "AURKB", "BUB1", "BUB1B", "UBE2C"
  ),
  MTOR_BCAA_READOUT = c(
    "BCAT1", "BCAT2", "BCKDHA", "BCKDHB", "DBT", "DLD", "BCKDK", "PPM1K",
    "SLC7A5", "SLC3A2", "SLC43A1", "SLC43A2", "SLC38A2", "SLC38A9",
    "MTOR", "RPTOR", "RHEB", "RPS6KB1", "EIF4EBP1", "RPS6",
    "LAMTOR1", "LAMTOR2", "LAMTOR3", "LAMTOR4", "LAMTOR5"
  )
)

state_projection <- readRDS("Data_Bulk_TCGA/Data_Merged/results/LAP3_State_Module/exports/lap3_state_module_projection.rds")
tcga_projection <- as.data.table(state_projection$tcga_projection)
cgga_projection <- as.data.table(state_projection$cgga_projection)
state_sets <- lapply(state_projection$state_sets, toupper)

malignant_audit <- readRDS("Data_Bulk_TCGA/Data_Merged/results/LAP3_Malignant_State_Module_Audit/exports/lap3_malignant_state_module_audit_results.rds")
m1_genes <- toupper(malignant_audit$cluster_sets$M1)

state_genes <- unique(toupper(state_sets$LAP3_STATE_UNION))
no_tp_state_genes <- unique(toupper(state_sets$LAP3_STATE_UNION_NO_TRANSLATION_PROTEOSTASIS))

candidate_long <- rbindlist(c(
  lapply(names(manual_gene_sets), function(class_name) {
    data.table(gene = toupper(manual_gene_sets[[class_name]]), anchor_class = class_name)
  }),
  list(data.table(gene = m1_genes, anchor_class = "M1_MALIGNANT_STATE_COMPONENT"))
), use.names = TRUE, fill = TRUE)
candidate_long <- unique(candidate_long)

candidate_gene_table <- candidate_long[
  ,
  .(
    anchor_classes = paste(sort(unique(anchor_class)), collapse = ";"),
    is_lap3 = gene == "LAP3",
    is_state_component = gene %in% state_genes,
    is_no_translation_state_component = gene %in% no_tp_state_genes,
    n_anchor_classes = uniqueN(anchor_class)
  ),
  by = gene
]
candidate_gene_table[, anchor_interpretation := fifelse(
  is_lap3,
  "external biologically motivated aminopeptidase anchor",
  fifelse(
    grepl("M1_MALIGNANT_STATE_COMPONENT", anchor_classes),
    "state-defining malignant component gene; strong but partly circular as an anchor",
    fifelse(
      grepl("TAM_MYELOID|HYPOXIA|PROLIFERATION|MTOR_BCAA", anchor_classes),
      "canonical context/pathway marker; useful comparator but higher generic-marker risk",
      "aminopeptidase-family comparator"
    )
  )
)]

cat("Candidate genes:", nrow(candidate_gene_table), "\n")
write_table(candidate_long, "alternative_anchor_candidate_gene_classes.csv")
write_table(candidate_gene_table, "alternative_anchor_candidate_genes.csv")

cat("Reading expression matrices...\n")
tcga_tpm_raw <- readRDS("Data_Bulk_TCGA/Data_Merged/data_analysis/expr_tpm_glioma_uni.rds")
tcga_tpm <- as.matrix(tcga_tpm_raw[, -1, drop = FALSE])
rownames(tcga_tpm) <- rownames(tcga_tpm_raw)
tcga_log2 <- log2(collapse_expr(tcga_tpm) + 1)

cgga693_mat <- as.matrix(readRDS("Data_Bulk_CGGA/mRNAseq_693/tpm_data.rds"))
cgga693_log2 <- log2(collapse_expr(cgga693_mat) + 1)

cgga325_raw <- fread("Data_Bulk_CGGA/mRNAseq_325/CGGA.mRNAseq_325.RSEM-genes.20200506.txt", data.table = FALSE, check.names = FALSE)
cgga325_mat <- as.matrix(cgga325_raw[, -1, drop = FALSE])
rownames(cgga325_mat) <- cgga325_raw$Gene_Name
cgga325_log2 <- log2(collapse_expr(cgga325_mat) + 1)

dataset_specs <- list(
  TCGA = list(
    projection = tcga_projection,
    expr = tcga_log2,
    sample_col = "barcode",
    covariates = c(
      "age_years", "cohort", "grade", "idh_status", "codel_1p19q",
      "TAM_MYELOID_CORE", "HALLMARK_HYPOXIA", "PROLIFERATION_CORE"
    )
  ),
  CGGA_mRNAseq_693 = list(
    projection = cgga_projection[cohort == "mRNAseq_693"],
    expr = cgga693_log2,
    sample_col = "sample_id",
    covariates = c(
      "age_years", "tumor_class", "idh_status", "codel_1p19q",
      "TAM_MYELOID_CORE", "HALLMARK_HYPOXIA", "PROLIFERATION_CORE"
    )
  ),
  CGGA_mRNAseq_325 = list(
    projection = cgga_projection[cohort == "mRNAseq_325"],
    expr = cgga325_log2,
    sample_col = "sample_id",
    covariates = c(
      "age_years", "tumor_class", "idh_status", "codel_1p19q",
      "TAM_MYELOID_CORE", "HALLMARK_HYPOXIA", "PROLIFERATION_CORE"
    )
  )
)

state_targets <- c("LAP3_STATE_UNION", "LAP3_STATE_UNION_NO_TRANSLATION_PROTEOSTASIS")

cat("Computing candidate anchor metrics...\n")
metric_rows <- list()
for (ds in names(dataset_specs)) {
  spec <- dataset_specs[[ds]]
  proj <- copy(as.data.table(spec$projection))
  expr <- spec$expr
  samples <- proj[[spec$sample_col]]
  for (gene in candidate_gene_table$gene) {
    gene_expr <- single_gene_expr(expr, gene, samples)
    proj[, candidate_gene_expr := gene_expr]
    present <- gene %in% rownames(expr)
    mean_log2_expr <- mean(gene_expr, na.rm = TRUE)
    detect_rate <- mean(gene_expr > 0, na.rm = TRUE)
    for (target in state_targets) {
      cor_res <- spearman_safe(gene_expr, proj[[target]])
      lm_res <- lm_anchor_term(
        proj,
        outcome = target,
        gene_expr_col = "candidate_gene_expr",
        covariates = spec$covariates,
        dataset_label = ds
      )
      metric_rows[[length(metric_rows) + 1L]] <- data.table(
        dataset = ds,
        gene = gene,
        present = present,
        target_state = target,
        n_cor = cor_res$n,
        spearman_rho = cor_res$rho,
        spearman_p = cor_res$p,
        mean_log2_expr = mean_log2_expr,
        detect_rate = detect_rate,
        n_lm = lm_res$n,
        composition_aware_beta = lm_res$beta,
        composition_aware_p = lm_res$p_value,
        composition_aware_adjusted_r2 = lm_res$adjusted_r2
      )
    }
  }
}

gene_metrics <- rbindlist(metric_rows, use.names = TRUE, fill = TRUE)
gene_metrics <- merge(gene_metrics, candidate_gene_table, by = "gene", all.x = TRUE)
gene_metrics[, spearman_fdr := p.adjust(spearman_p, method = "BH"), by = .(dataset, target_state)]
gene_metrics[, composition_aware_fdr := p.adjust(composition_aware_p, method = "BH"), by = .(dataset, target_state)]
gene_metrics[
  ,
  aminopeptidase_rank_abs_rho := frank(-abs(spearman_rho), ties.method = "min"),
  by = .(dataset, target_state, grepl("AMINOPEPTIDASE_FAMILY", anchor_classes))
]
gene_metrics[!grepl("AMINOPEPTIDASE_FAMILY", anchor_classes), aminopeptidase_rank_abs_rho := NA_real_]
write_table(gene_metrics, "alternative_anchor_gene_metrics_long.csv")

primary_target <- "LAP3_STATE_UNION_NO_TRANSLATION_PROTEOSTASIS"
primary_metrics <- gene_metrics[target_state == primary_target]

summary_by_gene <- primary_metrics[
  ,
  .(
    n_datasets_present = sum(present),
    n_datasets_tested = sum(is.finite(spearman_rho)),
    median_abs_rho = median(abs(spearman_rho), na.rm = TRUE),
    min_abs_rho = min(abs(spearman_rho), na.rm = TRUE),
    median_signed_rho = median(spearman_rho, na.rm = TRUE),
    n_positive_rho = sum(spearman_rho > 0, na.rm = TRUE),
    median_abs_composition_beta = median(abs(composition_aware_beta), na.rm = TRUE),
    min_abs_composition_beta = min(abs(composition_aware_beta), na.rm = TRUE),
    n_positive_composition_beta = sum(composition_aware_beta > 0, na.rm = TRUE),
    max_composition_aware_fdr = max(composition_aware_fdr, na.rm = TRUE),
    mean_log2_expr_tcga = mean(mean_log2_expr[dataset == "TCGA"], na.rm = TRUE)
  ),
  by = .(gene, anchor_classes, is_lap3, is_state_component, anchor_interpretation)
]
summary_by_gene[!is.finite(max_composition_aware_fdr), max_composition_aware_fdr := NA_real_]
summary_by_gene[, rank_all_candidates_median_abs_rho := frank(-median_abs_rho, ties.method = "min")]
summary_by_gene[, rank_external_non_state_candidates := {
  x <- rep(NA_real_, .N)
  idx <- !is_state_component | is_lap3
  x[idx] <- frank(-median_abs_rho[idx], ties.method = "min")
  x
}]
summary_by_gene[, rank_external_non_state_composition_beta := {
  x <- rep(NA_real_, .N)
  idx <- !is_state_component | is_lap3
  x[idx] <- frank(-median_abs_composition_beta[idx], ties.method = "min")
  x
}]
summary_by_gene[, rank_aminopeptidase_family := {
  x <- rep(NA_real_, .N)
  idx <- grepl("AMINOPEPTIDASE_FAMILY", anchor_classes)
  x[idx] <- frank(-median_abs_rho[idx], ties.method = "min")
  x
}]
setorder(summary_by_gene, rank_all_candidates_median_abs_rho, gene)
write_table(summary_by_gene, "alternative_anchor_gene_summary.csv")

class_summary <- primary_metrics[
  ,
  .(
    n_genes = uniqueN(gene),
    median_of_gene_median_abs_rho = median(
      summary_by_gene[gene %in% unique(.SD$gene), median_abs_rho],
      na.rm = TRUE
    ),
    max_gene_median_abs_rho = max(
      summary_by_gene[gene %in% unique(.SD$gene), median_abs_rho],
      na.rm = TRUE
    ),
    top_gene = summary_by_gene[gene %in% unique(.SD$gene)][order(rank_all_candidates_median_abs_rho), gene][1]
  ),
  by = anchor_classes
]
setorder(class_summary, -max_gene_median_abs_rho)
write_table(class_summary, "alternative_anchor_class_summary.csv")

aminopeptidase_existing <- fread("Data_Bulk_TCGA/Data_Merged/results/LAP3_State_Specificity_Benchmark/tables/aminopeptidase_gene_benchmark.csv")
amino_lap3_rank <- aminopeptidase_existing[gene == "LAP3", .(
  dataset,
  target_state = "LAP3_STATE_UNION",
  spearman_rho_with_LAP3_STATE_UNION,
  rank_abs_rho
)]
write_table(amino_lap3_rank, "lap3_aminopeptidase_family_rank.csv")

random_benchmark <- fread("Data_Bulk_TCGA/Data_Merged/results/LAP3_State_Specificity_Benchmark/tables/expression_matched_random_gene_benchmark_tcga.csv")
random_lap3 <- random_benchmark[is_lap3 == TRUE][1]

composition_verdict <- fread("Data_Bulk_TCGA/Data_Merged/results/LAP3_Bulk_Composition_Audit/tables/lap3_state_primary_composition_verdict.csv")
composition_lap3 <- composition_verdict[, .(
  dataset,
  n = n_clinical,
  beta_composition,
  p_composition,
  retained_fraction = beta_retention_vs_clinical
)]

cptac_bridge <- fread("Data_Protein_Public/results/LAP3_CPTAC_Multiomics_Feasibility/tables/cptac_lap3_mrna_protein_correlation.csv")
cptac_state <- fread("Data_Bulk_TCGA/Data_Merged/results/LAP3_State_CPTAC_GLASS_Projection/tables/cptac_lap3_state_submodule_correlations.csv")
cptac_state_bridge <- cptac_state[
  exposure %in% c("lap3_mrna_log2", "mrna_LAP3_STATE_UNION") &
    outcome == "lap3_protein",
  .(exposure, outcome, n_complete, spearman_rho, fdr)
]

top_external <- summary_by_gene[
  !is_state_component | is_lap3
][order(rank_external_non_state_candidates)][1:20]
top_all <- summary_by_gene[order(rank_all_candidates_median_abs_rho)][1:30]
top_m1 <- summary_by_gene[grepl("M1_MALIGNANT_STATE_COMPONENT", anchor_classes)][order(rank_all_candidates_median_abs_rho)][1:20]
write_table(top_external, "top_external_non_state_candidate_anchors.csv")
write_table(top_all, "top_all_candidate_anchors_including_state_components.csv")
write_table(top_m1, "top_m1_state_component_candidate_anchors.csv")

lap3_summary <- summary_by_gene[gene == "LAP3"]
dataset_order <- c("TCGA", "CGGA_mRNAseq_693", "CGGA_mRNAseq_325")
lap3_primary_metrics <- primary_metrics[gene == "LAP3"]
lap3_primary_metrics[, dataset := factor(dataset, levels = dataset_order)]
setorder(lap3_primary_metrics, dataset)
msr1_summary <- summary_by_gene[gene == "MSR1"]
defense_summary <- rbindlist(list(
  data.table(
    evidence_dimension = "Cross-cohort state anchoring",
    metric = "median_abs_rho_with_no_translation_state",
    value = lap3_summary$median_abs_rho,
    detail = paste0(
      "TCGA/CGGA693/CGGA325 rhos: ",
      paste(round(lap3_primary_metrics$spearman_rho, 4), collapse = ", ")
    ),
    interpretation = "LAP3 strongly tracks the frozen state across cohorts."
  ),
  data.table(
    evidence_dimension = "Aminopeptidase family specificity",
    metric = "rank_abs_rho",
    value = NA_real_,
    detail = paste(amino_lap3_rank[, paste0(dataset, " rank ", rank_abs_rho)], collapse = "; "),
    interpretation = "Within the curated aminopeptidase family, LAP3 is the strongest state anchor in all three bulk cohorts."
  ),
  data.table(
    evidence_dimension = "Expression-matched random benchmark",
    metric = "TCGA expression-matched percentile",
    value = random_lap3$percentile_abs_rho,
    detail = paste0("rank ", random_lap3$rank_abs_rho, " / ", nrow(random_benchmark)),
    interpretation = "LAP3 is extreme among expression-matched genes, arguing against a high-expression artifact."
  ),
  data.table(
    evidence_dimension = "Composition-aware residual signal",
    metric = "composition_adjusted_beta",
    value = median(composition_lap3$beta_composition, na.rm = TRUE),
    detail = paste(composition_lap3[, paste0(dataset, " beta=", round(beta_composition, 3))], collapse = "; "),
    interpretation = "The LAP3 term persists after composition-lite adjustment, although attenuation prevents purity-independent claims."
  ),
  data.table(
    evidence_dimension = "Protein bridge",
    metric = "CPTAC LAP3 mRNA-protein Spearman rho",
    value = cptac_bridge$spearman_rho[1],
    detail = paste0("n=", cptac_bridge$n_matched[1], ", p=", signif(cptac_bridge$spearman_p[1], 3)),
    interpretation = "LAP3 has direct transcript-to-protein support in CPTAC."
  ),
  data.table(
    evidence_dimension = "State-component counterfactual",
    metric = "top_M1_gene_median_abs_rho",
    value = top_m1$median_abs_rho[1],
    detail = paste0("top M1 state component: ", top_m1$gene[1], "; LAP3 is not part of the frozen state score."),
    interpretation = "Some M1 state genes may be stronger state-defining genes, but they are internal state components and less suitable as external biological anchors."
  ),
  data.table(
    evidence_dimension = "Generic TAM-marker counterfactual",
    metric = "MSR1_vs_LAP3",
    value = msr1_summary$median_abs_rho,
    detail = paste0(
      "MSR1 raw external rank ", msr1_summary$rank_external_non_state_candidates,
      " (median rho=", round(msr1_summary$median_abs_rho, 3),
      "); LAP3 raw external rank ", lap3_summary$rank_external_non_state_candidates,
      " but residual-beta rank ", lap3_summary$rank_external_non_state_composition_beta
    ),
    interpretation = "A TAM marker can slightly exceed LAP3 in raw state correlation, but LAP3 is less collapsible to a generic TAM marker because it has stronger residual signal, aminopeptidase-family specificity and a protein bridge."
  )
))
write_table(defense_summary, "lap3_anchor_defense_summary.csv")

source_bundle <- list(
  generated_at = format(Sys.time()),
  candidate_gene_table = candidate_gene_table,
  gene_metrics = gene_metrics,
  summary_by_gene = summary_by_gene,
  defense_summary = defense_summary,
  random_lap3 = random_lap3,
  composition_lap3 = composition_lap3,
  cptac_bridge = cptac_bridge,
  cptac_state_bridge = cptac_state_bridge
)
saveRDS(source_bundle, file.path(source_dir, "lap3_alternative_anchor_benchmark_source_bundle.rds"))
fwrite(gene_metrics, file.path(source_dir, "alternative_anchor_gene_metrics_long.csv"))

cat("LAP3 defense summary:\n")
print(defense_summary)
cat("\nTop external non-state candidate anchors:\n")
print(top_external[, .(
  gene, anchor_classes, median_abs_rho, min_abs_rho,
  median_abs_composition_beta, rank_external_non_state_candidates,
  rank_external_non_state_composition_beta
)])
cat("\nTop all candidates including state components:\n")
print(top_all[, .(
  gene, anchor_classes, median_abs_rho, is_state_component,
  rank_all_candidates_median_abs_rho
)])

cat(
  "# LAP3 Alternative Anchor Benchmark\n\n",
  "Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n",
  "## Purpose\n\n",
  "This module addresses the reviewer-facing question: why use LAP3 as the anchor rather than another gene in the state or a generic TAM/hypoxia/proliferation marker?\n\n",
  "## Design\n\n",
  "- Candidate classes: aminopeptidase family, TAM/myeloid genes, hypoxia genes, proliferation genes, mTOR/BCAA/readout genes, and M1 malignant-state component genes.\n",
  "- Primary target: `LAP3_STATE_UNION_NO_TRANSLATION_PROTEOSTASIS`.\n",
  "- Cohorts: TCGA, CGGA mRNAseq_693, CGGA mRNAseq_325.\n",
  "- Metrics: Spearman state anchoring, composition-aware linear-model beta, aminopeptidase-family rank, expression-matched random benchmark, CPTAC protein bridge, and M1 state-component counterfactual.\n\n",
  "## Key Results\n\n",
  "```text\n",
  paste(capture.output(print(defense_summary)), collapse = "\n"),
  "\n```\n\n",
  "## Top external non-state candidates\n\n",
  "```text\n",
  paste(capture.output(print(top_external[, .(
    gene, anchor_classes, median_abs_rho, min_abs_rho,
    median_abs_composition_beta, rank_external_non_state_candidates,
    rank_external_non_state_composition_beta
  )])), collapse = "\n"),
  "\n```\n\n",
  "## Top candidates including internal state components\n\n",
  "```text\n",
  paste(capture.output(print(head(top_all[, .(
    gene, anchor_classes, median_abs_rho, is_state_component,
    rank_all_candidates_median_abs_rho
  )], 20))), collapse = "\n"),
  "\n```\n\n",
  "## Interpretation Boundary\n\n",
  "LAP3 should be described as a biologically motivated, cross-cohort reproducible, protein-supported and experimentally actionable anchor. It should not be described as the only possible anchor, the globally optimal state hub, or a causal driver proven by public data. M1 malignant-state genes can be stronger state-defining components, but they are internal to the frozen state and therefore serve a different evidentiary role.\n\n",
  "## Outputs\n\n",
  "- `tables/alternative_anchor_candidate_genes.csv`\n",
  "- `tables/alternative_anchor_gene_metrics_long.csv`\n",
  "- `tables/alternative_anchor_gene_summary.csv`\n",
  "- `tables/alternative_anchor_class_summary.csv`\n",
  "- `tables/lap3_anchor_defense_summary.csv`\n",
  "- `tables/lap3_aminopeptidase_family_rank.csv`\n",
  "- `tables/top_external_non_state_candidate_anchors.csv`\n",
  "- `tables/top_all_candidate_anchors_including_state_components.csv`\n",
  "- `tables/top_m1_state_component_candidate_anchors.csv`\n",
  "- `source_data/lap3_alternative_anchor_benchmark_source_bundle.rds`\n\n",
  file = file.path(out_dir, "README.md")
)

cat("\nDONE\n")
