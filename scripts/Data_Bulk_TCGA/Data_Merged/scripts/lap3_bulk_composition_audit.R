#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

project_dir <- "/home/lzb/glioma"
setwd(project_dir)
data.table::setDTthreads(8)
set.seed(20260701)

out_dir <- "Data_Bulk_TCGA/Data_Merged/results/LAP3_Bulk_Composition_Audit"
table_dir <- file.path(out_dir, "tables")
source_dir <- file.path(out_dir, "source_data")
plot_dir <- file.path(out_dir, "plots")
export_dir <- file.path(out_dir, "exports")
log_dir <- file.path(out_dir, "logs")
for (d in c(table_dir, source_dir, plot_dir, export_dir, log_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

log_file <- file.path(log_dir, "lap3_bulk_composition_audit.log")
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message")
on.exit({
  sink(type = "message")
  sink(type = "output")
  close(log_con)
}, add = TRUE)

cat("LAP3 bulk composition / purity-like audit\n")
cat("Started:", format(Sys.time()), "\n")
cat("R:", R.version.string, "\n")
cat("Working directory:", getwd(), "\n\n")

write_table <- function(x, filename, dir = table_dir) {
  fwrite(as.data.table(x), file.path(dir, filename))
}

save_plot <- function(plot, filename, width = 7, height = 5) {
  ggsave(file.path(plot_dir, paste0(filename, ".pdf")), plot, width = width, height = height)
  ggsave(file.path(plot_dir, paste0(filename, ".png")), plot, width = width, height = height, dpi = 300)
  ggsave(file.path(plot_dir, paste0(filename, ".tiff")), plot, width = width, height = height, dpi = 300)
  if (requireNamespace("svglite", quietly = TRUE)) {
    ggsave(file.path(plot_dir, paste0(filename, ".svg")), plot, width = width, height = height)
  }
}

collapse_expr <- function(mat) {
  mat <- as.matrix(mat)
  storage.mode(mat) <- "double"
  rn <- rownames(mat)
  if (is.null(rn)) stop("Matrix rownames are required")
  keep <- !is.na(rn) & nzchar(rn)
  mat <- mat[keep, , drop = FALSE]
  rn <- rn[keep]
  if (!anyDuplicated(rn)) return(mat)
  dt <- as.data.table(mat)
  dt[, gene := rn]
  dt <- dt[, lapply(.SD, mean, na.rm = TRUE), by = gene]
  out <- as.matrix(dt[, -"gene"])
  rownames(out) <- dt$gene
  storage.mode(out) <- "double"
  out
}

module_score <- function(expr_log2, genes, min_genes = 3L) {
  genes <- intersect(unique(genes), rownames(expr_log2))
  if (length(genes) < min_genes) return(rep(NA_real_, ncol(expr_log2)))
  z <- t(scale(t(expr_log2[genes, , drop = FALSE])))
  z[!is.finite(z)] <- NA_real_
  colMeans(z, na.rm = TRUE)
}

single_gene_expr <- function(expr_log2, gene) {
  if (!gene %in% rownames(expr_log2)) return(rep(NA_real_, ncol(expr_log2)))
  as.numeric(expr_log2[gene, ])
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

composition_sets <- list(
  STROMAL_FIBROVASCULAR_CORE = c(
    "COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL6A1", "COL6A2", "FN1", "DCN",
    "LUM", "THY1", "ACTA2", "TAGLN", "PDGFRB", "RGS5", "MCAM", "VWF",
    "PECAM1", "ENG", "KDR", "CLDN5"
  ),
  IMMUNE_PANLEUKOCYTE_CORE = c(
    "PTPRC", "CD3D", "CD3E", "CD2", "CD8A", "CD4", "TRAC", "MS4A1", "CD79A",
    "NKG7", "GNLY", "LST1", "AIF1", "TYROBP", "FCER1G", "ITGAM", "CSF1R",
    "C1QA", "C1QB", "C1QC"
  ),
  TAM_MYELOID_EXTENDED_CORE = c(
    "AIF1", "CD68", "LST1", "CSF1R", "TYROBP", "FCER1G", "SPI1", "ITGAM",
    "C1QA", "C1QB", "C1QC", "CD163", "MRC1", "MSR1", "TREM2", "APOE",
    "LILRB4", "FCGR3A", "FCGR2A", "FCGR2B", "MYD88", "LILRB1", "LILRB2"
  ),
  HYPOXIA_EXTENDED_CORE = c(
    "CA9", "VEGFA", "SLC2A1", "LDHA", "ENO1", "PGK1", "BNIP3", "NDRG1",
    "P4HA1", "EGLN3", "ADM", "ANGPTL4", "BHLHE40", "HK2", "ALDOA", "GAPDH",
    "TPI1"
  ),
  PROLIFERATION_EXTENDED_CORE = c(
    "MKI67", "TOP2A", "PCNA", "MCM2", "MCM3", "MCM4", "MCM5", "MCM6",
    "CDK1", "CCNB1", "CCNB2", "AURKA", "AURKB", "BUB1", "BUB1B", "UBE2C"
  ),
  NEURAL_GLIAL_TUMOR_CONTENT_PROXY = c(
    "OLIG2", "SOX2", "SOX9", "EGFR", "PDGFRA", "GFAP", "ALDH1L1", "SLC1A3",
    "BCAN", "NES", "PTPRZ1", "SOX10"
  )
)

cat("Reading frozen projection tables...\n")
tcga <- fread("Data_Bulk_TCGA/Data_Merged/results/LAP3_State_Submodules/tables/tcga_lap3_state_submodule_projection.csv")
cgga <- fread("Data_Bulk_TCGA/Data_Merged/results/LAP3_State_Submodules/tables/cgga_lap3_state_submodule_projection.csv")

cat("Reading expression matrices for composition proxies...\n")
tcga_tpm_raw <- readRDS("Data_Bulk_TCGA/Data_Merged/data_analysis/expr_tpm_glioma_uni.rds")
tcga_tpm <- as.matrix(tcga_tpm_raw[, -1, drop = FALSE])
rownames(tcga_tpm) <- rownames(tcga_tpm_raw)
tcga_log2 <- log2(collapse_expr(tcga_tpm) + 1)

cgga_validation <- as.data.table(readRDS("Data_Bulk_CGGA/results/LAP3_CGGA/exports/cgga_lap3_validation_dataset.rds"))
cgga693_mat <- as.matrix(readRDS("Data_Bulk_CGGA/mRNAseq_693/tpm_data.rds"))
cgga693_log2 <- log2(collapse_expr(cgga693_mat) + 1)
cgga325_raw <- fread("Data_Bulk_CGGA/mRNAseq_325/CGGA.mRNAseq_325.RSEM-genes.20200506.txt", data.table = FALSE, check.names = FALSE)
cgga325_mat <- as.matrix(cgga325_raw[, -1, drop = FALSE])
rownames(cgga325_mat) <- cgga325_raw$Gene_Name
cgga325_log2 <- log2(collapse_expr(cgga325_mat) + 1)

score_composition <- function(projection, expr_log2, sample_col) {
  out <- copy(as.data.table(projection))
  out[, sample_key := get(sample_col)]
  for (nm in names(composition_sets)) {
    out[, (nm) := module_score(expr_log2, composition_sets[[nm]])[match(sample_key, colnames(expr_log2))]]
  }
  out[, LAP3_log2_expr_recomputed := single_gene_expr(expr_log2, "LAP3")[match(sample_key, colnames(expr_log2))]]
  out[, PURITY_LIKE_LOW_IMMUNE_STROMAL_PROXY := -0.5 * (
    z_num(IMMUNE_PANLEUKOCYTE_CORE) + z_num(STROMAL_FIBROVASCULAR_CORE)
  ), by = dataset]
  out
}

tcga_scored <- score_composition(tcga, tcga_log2, "barcode")
cgga693_scored <- score_composition(cgga[cohort == "mRNAseq_693"], cgga693_log2, "sample_id")
cgga325_scored <- score_composition(cgga[cohort == "mRNAseq_325"], cgga325_log2, "sample_id")
cgga_scored <- rbindlist(list(cgga693_scored, cgga325_scored), fill = TRUE)

tcga_scored[, analysis_dataset := "TCGA"]
cgga_scored[, analysis_dataset := dataset]
all_data <- rbindlist(list(tcga_scored, cgga_scored), fill = TRUE, use.names = TRUE)

coverage_table <- rbindlist(lapply(names(composition_sets), function(nm) {
  genes <- unique(composition_sets[[nm]])
  rbindlist(list(
    data.table(dataset = "TCGA", gene_set = nm, n_genes = length(genes), n_present = sum(genes %in% rownames(tcga_log2))),
    data.table(dataset = "CGGA_mRNAseq_693", gene_set = nm, n_genes = length(genes), n_present = sum(genes %in% rownames(cgga693_log2))),
    data.table(dataset = "CGGA_mRNAseq_325", gene_set = nm, n_genes = length(genes), n_present = sum(genes %in% rownames(cgga325_log2)))
  ))
}))
coverage_table[, coverage := n_present / n_genes]
write_table(coverage_table, "composition_proxy_gene_set_coverage.csv")

fwrite(all_data, file.path(source_dir, "bulk_composition_audit_analysis_dataset.csv"))
saveRDS(
  list(
    generated_at = format(Sys.time()),
    composition_sets = composition_sets,
    analysis_dataset = all_data,
    coverage = coverage_table
  ),
  file.path(export_dir, "lap3_bulk_composition_audit_dataset.rds")
)

score_cols <- c(
  "LAP3_STATE_UNION",
  "LAP3_STATE_UNION_NO_TRANSLATION_PROTEOSTASIS",
  "LAP3_MALIGNANT_STATE_MODULE",
  "LAP3_MYELOID_TAM_CONTEXT_MODULE",
  "LAP3_ANABOLIC_TRANSLATION_MODULE",
  "LAP3_PROTEOSTASIS_STRESS_MODULE",
  "LAP3_HYPOXIA_PERINECROTIC_MODULE"
)
proxy_cols <- c(
  "LAP3_log2_expr",
  "IMMUNE_PANLEUKOCYTE_CORE",
  "STROMAL_FIBROVASCULAR_CORE",
  "PURITY_LIKE_LOW_IMMUNE_STROMAL_PROXY",
  "TAM_MYELOID_EXTENDED_CORE",
  "TAM_MYELOID_CORE",
  "HYPOXIA_EXTENDED_CORE",
  "HYPOXIA_CORE",
  "PROLIFERATION_EXTENDED_CORE",
  "PROLIFERATION_CORE",
  "NEURAL_GLIAL_TUMOR_CONTENT_PROXY"
)
proxy_cols <- proxy_cols[proxy_cols %in% names(all_data)]

make_strata <- function(dt) {
  out <- list(data.table(stratum = "all", row_id = seq_len(nrow(dt))))
  if ("cohort" %in% names(dt)) {
    vals <- sort(unique(na.omit(as.character(dt$cohort))))
    out <- c(out, lapply(vals, function(v) data.table(stratum = paste0("cohort_", v), row_id = which(as.character(dt$cohort) == v))))
  }
  if ("tumor_class" %in% names(dt)) {
    vals <- sort(unique(na.omit(as.character(dt$tumor_class))))
    out <- c(out, lapply(vals, function(v) data.table(stratum = paste0("tumor_class_", v), row_id = which(as.character(dt$tumor_class) == v))))
  }
  if (all(c("cohort", "idh_status") %in% names(dt)) && any(as.character(dt$cohort) == "GBM", na.rm = TRUE)) {
    idh <- as.character(dt$idh_status)
    out <- c(out, list(data.table(stratum = "GBM_IDH_WT_or_Wildtype", row_id = which(as.character(dt$cohort) == "GBM" & idh %in% c("WT", "Wildtype")))))
  }
  if (all(c("tumor_class", "idh_status") %in% names(dt)) && any(as.character(dt$tumor_class) == "GBM_grade4", na.rm = TRUE)) {
    idh <- as.character(dt$idh_status)
    out <- c(out, list(data.table(stratum = "GBM_grade4_IDH_Wildtype", row_id = which(as.character(dt$tumor_class) == "GBM_grade4" & idh %in% c("WT", "Wildtype")))))
  }
  rbindlist(out, fill = TRUE)
}

split_by_dataset <- split(all_data, all_data$analysis_dataset)

cor_table <- rbindlist(lapply(names(split_by_dataset), function(ds) {
  d0 <- as.data.table(split_by_dataset[[ds]])
  strata <- make_strata(d0)
  rbindlist(lapply(unique(strata$stratum), function(st) {
    ids <- strata[stratum == st, row_id]
    d <- d0[ids]
    rbindlist(lapply(score_cols[score_cols %in% names(d)], function(score) {
      rbindlist(lapply(proxy_cols, function(proxy) {
        res <- spearman_safe(d[[score]], d[[proxy]], min_n = 25L)
        data.table(
          dataset = ds,
          stratum = st,
          score = score,
          proxy = proxy,
          n = res$n,
          spearman_rho = res$rho,
          p_value = res$p
        )
      }))
    }))
  }))
}), fill = TRUE)
cor_table[, p_adj_BH := p.adjust(p_value, method = "BH")]
write_table(cor_table, "lap3_state_composition_proxy_correlations.csv")

lm_one <- function(d, outcome, model_name, terms, term_of_interest = "LAP3_log2_expr", min_n = 50L) {
  terms <- terms[terms %in% names(d)]
  if (!term_of_interest %in% terms) {
    return(data.table(
      model = model_name, outcome = outcome, term = term_of_interest,
      n = nrow(d), beta = NA_real_, p_value = NA_real_, adjusted_r2 = NA_real_
    ))
  }
  terms <- terms[vapply(terms, function(nm) {
    x <- d[[nm]]
    if (nm == term_of_interest) return(TRUE)
    if (is.numeric(x)) return(sum(is.finite(x)) >= min_n && uniqueN(x[is.finite(x)]) >= 3L)
    uniqueN(na.omit(as.character(x))) >= 2L
  }, logical(1))]
  required <- unique(c(outcome, terms))
  d <- copy(as.data.table(d))
  d <- d[complete.cases(d[, ..required])]
  if (nrow(d) < min_n || !term_of_interest %in% terms || uniqueN(d[[term_of_interest]]) < 5L) {
    return(data.table(
      model = model_name, outcome = outcome, term = term_of_interest,
      n = nrow(d), beta = NA_real_, p_value = NA_real_, adjusted_r2 = NA_real_
    ))
  }

  numeric_terms <- terms[vapply(d[, ..terms], is.numeric, logical(1))]
  d[, (outcome) := z_num(get(outcome))]
  for (nm in numeric_terms) d[, (nm) := z_num(get(nm))]
  for (nm in setdiff(terms, numeric_terms)) d[, (nm) := factor(get(nm))]
  d <- d[complete.cases(d[, ..required])]
  if (nrow(d) < min_n) {
    return(data.table(
      model = model_name, outcome = outcome, term = term_of_interest,
      n = nrow(d), beta = NA_real_, p_value = NA_real_, adjusted_r2 = NA_real_
    ))
  }

  form <- reformulate(terms, response = outcome)
  fit <- tryCatch(lm(form, data = d), error = function(e) NULL)
  if (is.null(fit)) {
    return(data.table(
      model = model_name, outcome = outcome, term = term_of_interest,
      n = nrow(d), beta = NA_real_, p_value = NA_real_, adjusted_r2 = NA_real_
    ))
  }
  sm <- summary(fit)
  ct <- as.data.table(sm$coefficients, keep.rownames = "term")
  setnames(ct, c("Estimate", "Std. Error", "t value", "Pr(>|t|)"), c("beta", "std_error", "t_value", "p_value"))
  ct <- ct[term != "(Intercept)"]
  ct[, `:=`(model = model_name, outcome = outcome, n = nrow(d), adjusted_r2 = sm$adj.r.squared)]
  setcolorder(ct, c("model", "outcome", "term", "n", "beta", "std_error", "t_value", "p_value", "adjusted_r2"))
  ct
}

clinical_terms <- function(ds, d) {
  if (ds == "TCGA") {
    c("cohort", "grade", "idh_status", "codel_1p19q", "age_years")
  } else {
    c("tumor_class", "grade", "idh_status", "codel_1p19q", "age_years")
  }
}

model_table <- rbindlist(lapply(names(split_by_dataset), function(ds) {
  d0 <- as.data.table(split_by_dataset[[ds]])
  strata <- make_strata(d0)
  rbindlist(lapply(unique(strata$stratum), function(st) {
    ids <- strata[stratum == st, row_id]
    d <- d0[ids]
    base_clin <- clinical_terms(ds, d)
    model_defs <- list(
      LAP3_only = c("LAP3_log2_expr"),
      clinical = c("LAP3_log2_expr", base_clin),
      clinical_plus_immune_stromal = c("LAP3_log2_expr", base_clin, "IMMUNE_PANLEUKOCYTE_CORE", "STROMAL_FIBROVASCULAR_CORE"),
      clinical_plus_purity_like = c("LAP3_log2_expr", base_clin, "PURITY_LIKE_LOW_IMMUNE_STROMAL_PROXY"),
      clinical_plus_composition_lite = c(
        "LAP3_log2_expr", base_clin, "IMMUNE_PANLEUKOCYTE_CORE", "STROMAL_FIBROVASCULAR_CORE",
        "TAM_MYELOID_EXTENDED_CORE", "HYPOXIA_EXTENDED_CORE", "PROLIFERATION_EXTENDED_CORE"
      )
    )
    if (ds == "TCGA" && all(c("estimate_immune_score", "estimate_stromal_score") %in% names(d))) {
      model_defs$clinical_plus_tcga_estimate <- c("LAP3_log2_expr", base_clin, "estimate_immune_score", "estimate_stromal_score")
    }
    rbindlist(lapply(score_cols[score_cols %in% names(d)], function(outcome) {
      rbindlist(lapply(names(model_defs), function(mn) {
        res <- lm_one(d, outcome, mn, model_defs[[mn]])
        res[, `:=`(dataset = ds, stratum = st)]
        res
      }), fill = TRUE)
    }), fill = TRUE)
  }), fill = TRUE)
}), fill = TRUE)
model_table[, p_adj_BH := p.adjust(p_value, method = "BH")]
setcolorder(model_table, c("dataset", "stratum", "model", "outcome", "term", "n", "beta", "std_error", "t_value", "p_value", "p_adj_BH", "adjusted_r2"))
write_table(model_table, "lap3_state_composition_adjusted_lm_models.csv")

lap3_models <- model_table[term == "LAP3_log2_expr"]
lap3_wide <- dcast(
  lap3_models,
  dataset + stratum + outcome ~ model,
  value.var = c("beta", "p_value", "p_adj_BH", "n", "adjusted_r2")
)
if (!"beta_clinical" %in% names(lap3_wide)) lap3_wide[, beta_clinical := NA_real_]
if (!"beta_clinical_plus_composition_lite" %in% names(lap3_wide)) lap3_wide[, beta_clinical_plus_composition_lite := NA_real_]
lap3_wide[, beta_retention_vs_clinical := beta_clinical_plus_composition_lite / beta_clinical]
lap3_wide[, audit_call := fifelse(
  is.na(beta_clinical_plus_composition_lite),
  "insufficient_model",
  fifelse(beta_clinical_plus_composition_lite > 0 & p_value_clinical_plus_composition_lite < 0.05,
          "retained_after_composition_adjustment",
          fifelse(beta_clinical_plus_composition_lite > 0 & beta_retention_vs_clinical >= 0.25,
                  "attenuated_but_positive",
                  "composition_sensitive_or_not_retained"))
)]
write_table(lap3_wide, "lap3_state_lap3_term_model_verdict.csv")

primary_verdict <- lap3_wide[outcome == "LAP3_STATE_UNION" & stratum == "all", .(
  dataset,
  n_clinical = n_clinical,
  beta_clinical = round(beta_clinical, 3),
  p_clinical = signif(p_value_clinical, 3),
  beta_composition = round(beta_clinical_plus_composition_lite, 3),
  p_composition = signif(p_value_clinical_plus_composition_lite, 3),
  beta_retention_vs_clinical = round(beta_retention_vs_clinical, 3),
  audit_call
)]
write_table(primary_verdict, "lap3_state_primary_composition_verdict.csv")

pdat <- cor_table[
  stratum == "all" &
    score %in% c("LAP3_STATE_UNION", "LAP3_MALIGNANT_STATE_MODULE", "LAP3_MYELOID_TAM_CONTEXT_MODULE") &
    proxy %in% c("LAP3_log2_expr", "IMMUNE_PANLEUKOCYTE_CORE", "STROMAL_FIBROVASCULAR_CORE", "PURITY_LIKE_LOW_IMMUNE_STROMAL_PROXY", "TAM_MYELOID_EXTENDED_CORE", "HYPOXIA_EXTENDED_CORE", "PROLIFERATION_EXTENDED_CORE", "NEURAL_GLIAL_TUMOR_CONTENT_PROXY")
]
pdat[, score := factor(score, levels = c("LAP3_STATE_UNION", "LAP3_MALIGNANT_STATE_MODULE", "LAP3_MYELOID_TAM_CONTEXT_MODULE"))]
pdat[, proxy := factor(proxy, levels = rev(c("LAP3_log2_expr", "IMMUNE_PANLEUKOCYTE_CORE", "STROMAL_FIBROVASCULAR_CORE", "PURITY_LIKE_LOW_IMMUNE_STROMAL_PROXY", "TAM_MYELOID_EXTENDED_CORE", "HYPOXIA_EXTENDED_CORE", "PROLIFERATION_EXTENDED_CORE", "NEURAL_GLIAL_TUMOR_CONTENT_PROXY")))]
p_heat <- ggplot(pdat, aes(score, proxy, fill = spearman_rho)) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(aes(label = ifelse(is.finite(spearman_rho), sprintf("%.2f", spearman_rho), "")), size = 2.6) +
  facet_wrap(~ dataset, nrow = 1) +
  scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C", midpoint = 0, limits = c(-1, 1), na.value = "grey90") +
  labs(x = NULL, y = NULL, fill = "rho") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.grid = element_blank())
save_plot(p_heat, "lap3_state_composition_proxy_correlation_heatmap", width = 11, height = 4.8)

coef_dat <- lap3_models[
  stratum == "all" &
    outcome %in% c("LAP3_STATE_UNION", "LAP3_MALIGNANT_STATE_MODULE", "LAP3_MYELOID_TAM_CONTEXT_MODULE") &
    model %in% c("clinical", "clinical_plus_immune_stromal", "clinical_plus_purity_like", "clinical_plus_composition_lite")
]
coef_dat[, model := factor(model, levels = c("clinical", "clinical_plus_immune_stromal", "clinical_plus_purity_like", "clinical_plus_composition_lite"))]
p_coef <- ggplot(coef_dat, aes(model, beta, color = dataset)) +
  geom_hline(yintercept = 0, linewidth = 0.25, color = "grey40") +
  geom_point(position = position_dodge(width = 0.55), size = 1.8) +
  geom_errorbar(aes(ymin = beta - 1.96 * std_error, ymax = beta + 1.96 * std_error),
                position = position_dodge(width = 0.55), width = 0.2, linewidth = 0.35) +
  facet_wrap(~ outcome, ncol = 1, scales = "free_y") +
  coord_flip() +
  labs(x = NULL, y = "Standardized beta for LAP3 expression", color = NULL) +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank())
save_plot(p_coef, "lap3_term_composition_adjusted_model_forest", width = 8.5, height = 7)

readme <- c(
  "# LAP3 Bulk Composition Audit",
  "",
  paste0("Generated: ", format(Sys.time())),
  "",
  "## Purpose",
  "",
  "This module audits whether the frozen `LAP3_STATE_UNION` bulk signal can be reduced to a trivial immune/stromal/TAM or purity-like artifact.",
  "",
  "## Inputs",
  "",
  "- `Data_Bulk_TCGA/Data_Merged/results/LAP3_State_Submodules/tables/tcga_lap3_state_submodule_projection.csv`",
  "- `Data_Bulk_TCGA/Data_Merged/results/LAP3_State_Submodules/tables/cgga_lap3_state_submodule_projection.csv`",
  "- `Data_Bulk_TCGA/Data_Merged/data_analysis/expr_tpm_glioma_uni.rds`",
  "- `Data_Bulk_CGGA/mRNAseq_693/tpm_data.rds`",
  "- `Data_Bulk_CGGA/mRNAseq_325/CGGA.mRNAseq_325.RSEM-genes.20200506.txt`",
  "",
  "## Methods",
  "",
  "- Computed marker-based immune, stromal/fibrovascular, TAM/myeloid, hypoxia, proliferation and neural/glial tumor-content proxy scores using mean gene-wise z-scored log2 expression.",
  "- Defined `PURITY_LIKE_LOW_IMMUNE_STROMAL_PROXY` as the negative mean of standardized immune and stromal scores within each dataset.",
  "- Tested Spearman correlations between LAP3-state/submodules and composition proxies across TCGA, CGGA693 and CGGA325 strata.",
  "- Fitted standardized linear models for state outcomes, comparing clinical-only models with composition-adjusted models.",
  "- Used TCGA built-in `estimate_immune_score` and `estimate_stromal_score` only as an optional sensitivity model where complete cases were available.",
  "",
  "## Key Outputs",
  "",
  "- `source_data/bulk_composition_audit_analysis_dataset.csv`",
  "- `tables/composition_proxy_gene_set_coverage.csv`",
  "- `tables/lap3_state_composition_proxy_correlations.csv`",
  "- `tables/lap3_state_composition_adjusted_lm_models.csv`",
  "- `tables/lap3_state_lap3_term_model_verdict.csv`",
  "- `tables/lap3_state_primary_composition_verdict.csv`",
  "- `plots/lap3_state_composition_proxy_correlation_heatmap.*`",
  "- `plots/lap3_term_composition_adjusted_model_forest.*`",
  "- `exports/lap3_bulk_composition_audit_dataset.rds`",
  "",
  "## Primary Verdict",
  "",
  paste(capture.output(print(primary_verdict)), collapse = "\n"),
  "",
  "## Interpretation Boundary",
  "",
  "This is a lightweight composition audit, not a full cell-fraction deconvolution benchmark. It should be used to test whether the LAP3-state is obviously reducible to broad immune/stromal abundance. It should not be used to claim purity-independent causal biology. A positive retained LAP3 term supports a non-trivial LAP3-centered state within a malignant-microenvironment framework; attenuation supports the current ecosystem interpretation."
)
writeLines(readme, file.path(out_dir, "README.md"))

cat("\nPrimary verdict:\n")
print(primary_verdict)
cat("\nCompleted:", format(Sys.time()), "\n")
