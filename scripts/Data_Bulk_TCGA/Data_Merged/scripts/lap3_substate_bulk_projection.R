#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(survival)
})

project_dir <- "/home/lzb/glioma"
setwd(project_dir)
data.table::setDTthreads(16)
set.seed(20260630)

out_dir <- "Data_Bulk_TCGA/Data_Merged/results/LAP3_Substate_Bulk_Projection"
table_dir <- file.path(out_dir, "tables")
export_dir <- file.path(out_dir, "exports")
log_dir <- file.path(out_dir, "logs")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(log_dir, "lap3_substate_bulk_projection.log")
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message")
on.exit({
  sink(type = "message")
  sink(type = "output")
  close(log_con)
}, add = TRUE)

cat("Started:", format(Sys.time()), "\n")
cat("R:", R.version.string, "\n")

cache_file <- "Data_scRNA_GEO/GBmap_Core/cache/core_gbmap_lap3_cellstate_lightweight.rds"
tcga_clinical_file <- "Data_Bulk_TCGA/Data_Merged/results/Clinical_Field_QC/clinical_glioma_analysis_fields.rds"
tcga_tpm_file <- "Data_Bulk_TCGA/Data_Merged/data_analysis/expr_tpm_glioma_uni.rds"
cgga_693_tpm_file <- "Data_Bulk_CGGA/mRNAseq_693/tpm_data.rds"
cgga_693_clinical_file <- "Data_Bulk_CGGA/mRNAseq_693/CGGA.mRNAseq_693_clinical.20200506.txt"
cgga_325_rsem_file <- "Data_Bulk_CGGA/mRNAseq_325/CGGA.mRNAseq_325.RSEM-genes.20200506.txt"
cgga_325_clinical_file <- "Data_Bulk_CGGA/mRNAseq_325/CGGA.mRNAseq_325_clinical.20200506.txt"
tcga_pathway_file <- "Data_Bulk_TCGA/Data_Merged/results/LAP3_Pathway/tables/lap3_pathway_module_scores.csv"
cgga_validation_file <- "Data_Bulk_CGGA/results/LAP3_CGGA/exports/cgga_lap3_validation_dataset.rds"

write_table <- function(x, filename) {
  fwrite(x, file.path(table_dir, filename))
}

clean_na <- function(x) {
  x <- trimws(as.character(x))
  x[x == "" | x == "NA" | x == "NaN" | is.na(x)] <- NA_character_
  x
}

map_cgga_grade <- function(x) {
  x <- clean_na(x)
  fifelse(x == "WHO II", "G2",
    fifelse(x == "WHO III", "G3",
      fifelse(x == "WHO IV", "G4", NA_character_)
    )
  )
}

score_signature <- function(expr_log2, up_genes, down_genes = character()) {
  up <- intersect(up_genes, rownames(expr_log2))
  down <- intersect(down_genes, rownames(expr_log2))
  if (length(up) < 3L) {
    return(rep(NA_real_, ncol(expr_log2)))
  }
  z_up <- t(scale(t(expr_log2[up, , drop = FALSE])))
  z_up[!is.finite(z_up)] <- 0
  up_score <- colMeans(z_up, na.rm = TRUE)
  if (length(down) >= 3L) {
    z_down <- t(scale(t(expr_log2[down, , drop = FALSE])))
    z_down[!is.finite(z_down)] <- 0
    down_score <- colMeans(z_down, na.rm = TRUE)
    return(up_score - down_score)
  }
  up_score
}

spearman_safe <- function(x, y, min_n = 20L) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) < min_n || length(unique(x)) < 3L || length(unique(y)) < 3L) {
    return(list(n = length(x), rho = NA_real_, p = NA_real_))
  }
  test <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  list(n = length(x), rho = unname(test$estimate), p = test$p.value)
}

lm_term <- function(data, y, x, covariates, model_name) {
  required <- c(y, x, covariates)
  d <- data[complete.cases(data[, ..required])]
  if (nrow(d) < 30L || length(unique(d[[x]])) < 3L) {
    return(data.table(model = model_name, outcome = y, term = x, n = nrow(d), beta = NA_real_, p_value = NA_real_))
  }
  for (cc in covariates) {
    if (!is.numeric(d[[cc]])) d[[cc]] <- factor(d[[cc]])
  }
  fit <- lm(reformulate(c(x, covariates), response = y), data = d)
  coef_table <- summary(fit)$coefficients
  if (!x %in% rownames(coef_table)) {
    return(data.table(model = model_name, outcome = y, term = x, n = nrow(d), beta = NA_real_, p_value = NA_real_))
  }
  data.table(
    model = model_name,
    outcome = y,
    term = x,
    n = nrow(d),
    beta = coef_table[x, "Estimate"],
    std_error = coef_table[x, "Std. Error"],
    p_value = coef_table[x, "Pr(>|t|)"]
  )
}

cat("Reading Core GBmap lightweight cache...\n")
cache <- readRDS(cache_file)
obs <- as.data.table(cache$obs)
expr <- cache$normalized
raw <- cache$raw
obs[, author_donor := paste(author, donor_id, sep = "::")]
obs[, author_state := fifelse(annotation_level_3 == "AC-like" | grepl("^AC-like", annotation_level_4), "AC",
  fifelse(annotation_level_3 == "OPC-like" | grepl("^OPC-like", annotation_level_4), "OPC",
    fifelse(annotation_level_3 == "NPC-like" | grepl("^NPC-like", annotation_level_4), "NPC",
      fifelse(annotation_level_3 == "MES-like" | grepl("^MES-like", annotation_level_4), "MES", NA_character_)
    )
  )
)]
main_index <- which(obs$author != "Neftel2019" & obs$annotation_level_1 == "Neoplastic" & !is.na(obs$author_state))
obs_main <- copy(obs[main_index])
expr_main <- expr[main_index, , drop = FALSE]
raw_main <- raw[main_index, , drop = FALSE]
lap3_raw <- as.numeric(raw_main[, "LAP3"])
obs_main[, lap3_detected := lap3_raw > 0]
obs_main[, target_raw_sum := Matrix::rowSums(raw_main)]

genes <- setdiff(colnames(expr_main), "LAP3")
cat("Main cells:", nrow(obs_main), "Genes for signature:", length(genes), "\n")

cat("Deriving donor-state LAP3-detected signature deltas...\n")
delta_rows <- list()
k <- 0L
for (state in sort(unique(obs_main$author_state))) {
  state_indices <- which(obs_main$author_state == state)
  donor_states <- unique(obs_main$author_donor[state_indices])
  for (ds in donor_states) {
    idx <- state_indices[obs_main$author_donor[state_indices] == ds]
    high <- idx[obs_main$lap3_detected[idx]]
    low <- idx[!obs_main$lap3_detected[idx]]
    if (length(high) < 10L || length(low) < 10L) next
    high_mean <- Matrix::colMeans(expr_main[high, genes, drop = FALSE])
    low_mean <- Matrix::colMeans(expr_main[low, genes, drop = FALSE])
    meta <- obs_main[idx[1L], .(author, donor_id, author_donor, author_state)]
    k <- k + 1L
    delta_rows[[k]] <- data.table(
      meta,
      gene = genes,
      delta = as.numeric(high_mean - low_mean),
      high_cells = length(high),
      low_cells = length(low),
      target_raw_sum_delta = mean(obs_main$target_raw_sum[high]) - mean(obs_main$target_raw_sum[low])
    )
  }
}
delta_dt <- rbindlist(delta_rows, use.names = TRUE)
fwrite(delta_dt, file.path(export_dir, "core_lap3_detected_donor_state_gene_deltas.csv.gz"))

signature_stats <- delta_dt[
  ,
  .(
    n_donor_states = .N,
    n_authors = uniqueN(author),
    median_delta = median(delta, na.rm = TRUE),
    mean_delta = mean(delta, na.rm = TRUE),
    positive_fraction = mean(delta > 0, na.rm = TRUE),
    p_value = suppressWarnings(wilcox.test(delta, mu = 0, exact = FALSE)$p.value),
    median_target_raw_sum_delta = median(target_raw_sum_delta, na.rm = TRUE)
  ),
  by = gene
][n_donor_states >= 20L]
signature_stats[, p_adj_BH := p.adjust(p_value, method = "BH")]
signature_stats[, abs_median_delta := abs(median_delta)]
setorder(signature_stats, p_adj_BH, -abs_median_delta)
write_table(signature_stats, "core_lap3_detected_signature_gene_stats.csv")

up_genes <- signature_stats[
  median_delta > 0 & positive_fraction >= 0.65 & p_adj_BH < 0.05
][order(-median_delta)][1:min(.N, 60), gene]
down_genes <- signature_stats[
  median_delta < 0 & positive_fraction <= 0.35 & p_adj_BH < 0.05
][order(median_delta)][1:min(.N, 60), gene]
if (length(up_genes) < 10L) {
  up_genes <- signature_stats[median_delta > 0][order(-median_delta)][1:min(.N, 40), gene]
}
if (length(down_genes) < 10L) {
  down_genes <- signature_stats[median_delta < 0][order(median_delta)][1:min(.N, 40), gene]
}

signature_gene_table <- rbindlist(list(
  data.table(signature = "LAP3_detected_substate_up", gene = up_genes, direction = "up"),
  data.table(signature = "LAP3_detected_substate_down", gene = down_genes, direction = "down")
))
signature_gene_table <- merge(signature_gene_table, signature_stats, by = "gene", all.x = TRUE)
write_table(signature_gene_table, "core_lap3_detected_substate_signature_genes.csv")

cat("Signature up genes:", length(up_genes), "down genes:", length(down_genes), "\n")

read_tcga_expr <- function() {
  clinical <- as.data.table(readRDS(tcga_clinical_file))
  expr_tpm <- readRDS(tcga_tpm_file)
  stopifnot(identical(clinical$barcode, colnames(expr_tpm)[-1]))
  mat <- as.matrix(expr_tpm[, -1, drop = FALSE])
  storage.mode(mat) <- "double"
  rownames(mat) <- rownames(expr_tpm)
  list(clinical = clinical, expr = mat)
}

read_cgga_693 <- function() {
  validation <- as.data.table(readRDS(cgga_validation_file))
  validation <- validation[cohort == "mRNAseq_693"]
  mat <- as.matrix(readRDS(cgga_693_tpm_file))
  storage.mode(mat) <- "double"
  common <- intersect(validation$sample_id, colnames(mat))
  validation <- validation[match(common, sample_id)]
  mat <- mat[, validation$sample_id, drop = FALSE]
  list(clinical = validation, expr = mat)
}

read_cgga_325 <- function() {
  validation <- as.data.table(readRDS(cgga_validation_file))
  validation <- validation[cohort == "mRNAseq_325"]
  rsem <- fread(cgga_325_rsem_file, data.table = FALSE, check.names = FALSE)
  gene_names <- rsem$Gene_Name
  mat <- as.matrix(rsem[, -1, drop = FALSE])
  storage.mode(mat) <- "double"
  rownames(mat) <- gene_names
  common <- intersect(validation$sample_id, colnames(mat))
  validation <- validation[match(common, sample_id)]
  mat <- mat[, validation$sample_id, drop = FALSE]
  list(clinical = validation, expr = mat)
}

cat("Projecting signature to TCGA and CGGA...\n")
tcga <- read_tcga_expr()
cgga693 <- read_cgga_693()
cgga325 <- read_cgga_325()

make_projection <- function(dataset_name, clinical, expr) {
  expr_log2 <- log2(expr + 1)
  signature_score <- score_signature(expr_log2, up_genes, down_genes)
  up_score <- score_signature(expr_log2, up_genes, character())
  down_score <- score_signature(expr_log2, down_genes, character())
  lap3_expr <- if ("LAP3" %in% rownames(expr_log2)) as.numeric(expr_log2["LAP3", ]) else rep(NA_real_, ncol(expr_log2))
  out <- copy(clinical)
  sample_col <- if ("barcode" %in% names(out)) "barcode" else "sample_id"
  out[, dataset := dataset_name]
  out[, sample_key := get(sample_col)]
  out[, LAP3_log2_expr_projection := lap3_expr[match(sample_key, colnames(expr_log2))]]
  out[, LAP3_detected_substate_signature := signature_score[match(sample_key, colnames(expr_log2))]]
  out[, LAP3_detected_substate_up_score := up_score[match(sample_key, colnames(expr_log2))]]
  out[, LAP3_detected_substate_down_score := down_score[match(sample_key, colnames(expr_log2))]]
  out
}

tcga_proj <- make_projection("TCGA", tcga$clinical, tcga$expr)
cgga_proj <- rbindlist(list(
  make_projection("CGGA_mRNAseq_693", cgga693$clinical, cgga693$expr),
  make_projection("CGGA_mRNAseq_325", cgga325$clinical, cgga325$expr)
), fill = TRUE)

tcga_path <- fread(tcga_pathway_file)
tcga_path[, sample_key := barcode]
path_cols_tcga <- intersect(c("HALLMARK_MTORC1_SIGNALING", "LEUCINE_BCAA_CORE", "MTORC1_READOUT_CORE", "REACTOME_TRANSLATION"), names(tcga_path))
tcga_proj <- merge(
  tcga_proj,
  tcga_path[, c("sample_key", path_cols_tcga), with = FALSE],
  by = "sample_key",
  all.x = TRUE,
  suffixes = c("", "_pathway")
)

write_table(tcga_proj, "tcga_lap3_substate_signature_projection.csv")
write_table(cgga_proj, "cgga_lap3_substate_signature_projection.csv")
saveRDS(list(tcga = tcga_proj, cgga = cgga_proj, up_genes = up_genes, down_genes = down_genes),
        file.path(export_dir, "lap3_substate_bulk_projection_dataset.rds"))

score_col <- "LAP3_detected_substate_signature"
pathways <- c("HALLMARK_MTORC1_SIGNALING", "LEUCINE_BCAA_CORE", "MTORC1_READOUT_CORE", "REACTOME_TRANSLATION")

tcga_cor <- rbindlist(lapply(c("all", sort(unique(tcga_proj$cohort))), function(group) {
  d <- if (group == "all") tcga_proj else tcga_proj[cohort == group]
  rbindlist(lapply(c("LAP3_log2_expr_projection", pathways, "estimate_immune_score", "estimate_stromal_score"), function(v) {
    if (!v %in% names(d)) return(NULL)
    res <- spearman_safe(d[[score_col]], d[[v]], min_n = 30L)
    data.table(dataset = "TCGA", group = group, variable = v, n = res$n, spearman_rho = res$rho, p_value = res$p)
  }), fill = TRUE)
}), fill = TRUE)
tcga_cor[, p_adj_BH := p.adjust(p_value, method = "BH")]
write_table(tcga_cor, "tcga_lap3_substate_signature_correlations.csv")

cgga_cor <- rbindlist(lapply(sort(unique(cgga_proj$dataset)), function(ds) {
  d0 <- cgga_proj[dataset == ds]
  rbindlist(lapply(c("all", sort(unique(as.character(d0$tumor_class)))), function(group) {
    d <- if (group == "all") d0 else d0[as.character(tumor_class) == group]
    rbindlist(lapply(c("LAP3_log2_expr_projection", pathways), function(v) {
      if (!v %in% names(d)) return(NULL)
      res <- spearman_safe(d[[score_col]], d[[v]], min_n = 25L)
      data.table(dataset = ds, group = group, variable = v, n = res$n, spearman_rho = res$rho, p_value = res$p)
    }), fill = TRUE)
  }), fill = TRUE)
}), fill = TRUE)
cgga_cor[, p_adj_BH := p.adjust(p_value, method = "BH")]
write_table(cgga_cor, "cgga_lap3_substate_signature_correlations.csv")

tcga_lm <- rbindlist(lapply(c("HALLMARK_MTORC1_SIGNALING", "LEUCINE_BCAA_CORE", "MTORC1_READOUT_CORE", "REACTOME_TRANSLATION"), function(y) {
  if (!y %in% names(tcga_proj)) return(NULL)
  lm_term(tcga_proj, y = y, x = score_col, covariates = c("cohort", "age_years", "idh_status", "codel_1p19q"), model_name = "TCGA_adjusted_tumor_class_molecular")
}), fill = TRUE)
tcga_lm[, p_adj_BH := p.adjust(p_value, method = "BH")]
write_table(tcga_lm, "tcga_lap3_substate_signature_adjusted_lm.csv")

cgga_lm <- rbindlist(lapply(sort(unique(cgga_proj$dataset)), function(ds) {
  d <- cgga_proj[dataset == ds]
  rbindlist(lapply(c("HALLMARK_MTORC1_SIGNALING", "LEUCINE_BCAA_CORE", "MTORC1_READOUT_CORE", "REACTOME_TRANSLATION"), function(y) {
    if (!y %in% names(d)) return(NULL)
    lm_term(d, y = y, x = score_col, covariates = c("tumor_class", "age_years", "idh_status", "codel_1p19q"), model_name = paste0(ds, "_adjusted_tumor_class_molecular"))
  }), fill = TRUE)
}), fill = TRUE)
cgga_lm[, p_adj_BH := p.adjust(p_value, method = "BH")]
write_table(cgga_lm, "cgga_lap3_substate_signature_adjusted_lm.csv")

survival_rows <- list()
if (all(c("os_months", "os_event") %in% names(tcga_proj))) {
  d <- tcga_proj[complete.cases(tcga_proj[, .(os_months, os_event, LAP3_detected_substate_signature, cohort, age_years, idh_status, codel_1p19q)])]
  if (nrow(d) >= 100L) {
    fit <- coxph(Surv(os_months, os_event) ~ LAP3_detected_substate_signature + cohort + age_years + idh_status + codel_1p19q, data = d)
    s <- summary(fit)
    survival_rows[["TCGA"]] <- data.table(
      dataset = "TCGA",
      n = s$n,
      events = s$nevent,
      term = rownames(s$coefficients),
      HR = s$conf.int[, "exp(coef)"],
      conf_low = s$conf.int[, "lower .95"],
      conf_high = s$conf.int[, "upper .95"],
      p_value = s$coefficients[, "Pr(>|z|)"]
    )[term == "LAP3_detected_substate_signature"]
  }
}
for (ds in sort(unique(cgga_proj$dataset))) {
  d <- cgga_proj[dataset == ds]
  d <- d[complete.cases(d[, .(os_months, os_event, LAP3_detected_substate_signature, tumor_class, age_years, idh_status, codel_1p19q)])]
  if (nrow(d) >= 100L && length(unique(d$os_event)) == 2L) {
    fit <- coxph(Surv(os_months, os_event) ~ LAP3_detected_substate_signature + tumor_class + age_years + idh_status + codel_1p19q, data = d)
    s <- summary(fit)
    survival_rows[[ds]] <- data.table(
      dataset = ds,
      n = s$n,
      events = s$nevent,
      term = rownames(s$coefficients),
      HR = s$conf.int[, "exp(coef)"],
      conf_low = s$conf.int[, "lower .95"],
      conf_high = s$conf.int[, "upper .95"],
      p_value = s$coefficients[, "Pr(>|z|)"]
    )[term == "LAP3_detected_substate_signature"]
  }
}
survival_table <- rbindlist(survival_rows, fill = TRUE)
if (nrow(survival_table)) survival_table[, p_adj_BH := p.adjust(p_value, method = "BH")]
write_table(survival_table, "lap3_substate_signature_survival_cox.csv")

coverage <- rbindlist(list(
  data.table(dataset = "TCGA", direction = "up", requested = length(up_genes), available = length(intersect(up_genes, rownames(tcga$expr)))),
  data.table(dataset = "TCGA", direction = "down", requested = length(down_genes), available = length(intersect(down_genes, rownames(tcga$expr)))),
  data.table(dataset = "CGGA_mRNAseq_693", direction = "up", requested = length(up_genes), available = length(intersect(up_genes, rownames(cgga693$expr)))),
  data.table(dataset = "CGGA_mRNAseq_693", direction = "down", requested = length(down_genes), available = length(intersect(down_genes, rownames(cgga693$expr)))),
  data.table(dataset = "CGGA_mRNAseq_325", direction = "up", requested = length(up_genes), available = length(intersect(up_genes, rownames(cgga325$expr)))),
  data.table(dataset = "CGGA_mRNAseq_325", direction = "down", requested = length(down_genes), available = length(intersect(down_genes, rownames(cgga325$expr))))
))
coverage[, coverage := available / requested]
write_table(coverage, "lap3_substate_signature_bulk_gene_coverage.csv")

readme <- c(
  "# LAP3 Substate Bulk Projection",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Purpose",
  "",
  "Bridge Core GBmap LAP3-detected/high malignant substate to bulk TCGA/CGGA dry-lab evidence.",
  "",
  "## Inputs",
  "",
  paste0("- Core cache: `", cache_file, "`"),
  paste0("- TCGA clinical/expression: `", tcga_clinical_file, "`, `", tcga_tpm_file, "`"),
  "- CGGA mRNAseq_693 and mRNAseq_325 RSEM/clinical inputs",
  "",
  "## Methods",
  "",
  "- Derive gene-level LAP3-detected vs not-detected deltas within donor-state units in Core GBmap.",
  "- Select consensus up/down genes from 725 cached target genes, excluding LAP3.",
  "- Project an up-minus-down z-score signature into TCGA and CGGA bulk TPM/RSEM data.",
  "- Test associations with LAP3 expression, mTORC1/BCAA/translation module scores, clinical-molecular covariates, and survival.",
  "",
  "## Boundary",
  "",
  "The signature is limited to genes available in the Core lightweight cache and can carry LAP3 detection/RNA-content effects. It supports dry-lab mechanistic continuity but does not establish causal LAP3-leucine-mTORC1 signaling.",
  "",
  "## Key outputs",
  "",
  "- `tables/core_lap3_detected_substate_signature_genes.csv`",
  "- `tables/tcga_lap3_substate_signature_correlations.csv`",
  "- `tables/cgga_lap3_substate_signature_correlations.csv`",
  "- `tables/tcga_lap3_substate_signature_adjusted_lm.csv`",
  "- `tables/cgga_lap3_substate_signature_adjusted_lm.csv`",
  "- `tables/lap3_substate_signature_survival_cox.csv`"
)
writeLines(readme, file.path(out_dir, "README.md"))

cat("Finished:", format(Sys.time()), "\n")
