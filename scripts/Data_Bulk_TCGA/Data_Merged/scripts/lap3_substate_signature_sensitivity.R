#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

project_dir <- "/home/lzb/glioma"
setwd(project_dir)
data.table::setDTthreads(16)
set.seed(20260630)

out_dir <- "Data_Bulk_TCGA/Data_Merged/results/LAP3_Substate_Bulk_Projection"
table_dir <- file.path(out_dir, "tables")
log_dir <- file.path(out_dir, "logs")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(log_dir, "lap3_substate_signature_sensitivity.log")
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

tcga_clinical_file <- "Data_Bulk_TCGA/Data_Merged/results/Clinical_Field_QC/clinical_glioma_analysis_fields.rds"
tcga_tpm_file <- "Data_Bulk_TCGA/Data_Merged/data_analysis/expr_tpm_glioma_uni.rds"
tcga_pathway_file <- "Data_Bulk_TCGA/Data_Merged/results/LAP3_Pathway/tables/lap3_pathway_module_scores.csv"
cgga_693_tpm_file <- "Data_Bulk_CGGA/mRNAseq_693/tpm_data.rds"
cgga_325_rsem_file <- "Data_Bulk_CGGA/mRNAseq_325/CGGA.mRNAseq_325.RSEM-genes.20200506.txt"
cgga_validation_file <- "Data_Bulk_CGGA/results/LAP3_CGGA/exports/cgga_lap3_validation_dataset.rds"
signature_file <- file.path(table_dir, "core_lap3_detected_substate_signature_genes.csv")

write_table <- function(x, filename) {
  fwrite(x, file.path(table_dir, filename))
}

score_signature <- function(expr_log2, up_genes, down_genes = character()) {
  up <- intersect(up_genes, rownames(expr_log2))
  down <- intersect(down_genes, rownames(expr_log2))
  if (length(up) < 3L) return(rep(NA_real_, ncol(expr_log2)))
  z_up <- t(scale(t(expr_log2[up, , drop = FALSE])))
  z_up[!is.finite(z_up)] <- 0
  up_score <- colMeans(z_up, na.rm = TRUE)
  if (length(down) >= 3L) {
    z_down <- t(scale(t(expr_log2[down, , drop = FALSE])))
    z_down[!is.finite(z_down)] <- 0
    return(up_score - colMeans(z_down, na.rm = TRUE))
  }
  up_score
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
  mat <- as.matrix(rsem[, -1, drop = FALSE])
  storage.mode(mat) <- "double"
  rownames(mat) <- rsem$Gene_Name
  common <- intersect(validation$sample_id, colnames(mat))
  validation <- validation[match(common, sample_id)]
  mat <- mat[, validation$sample_id, drop = FALSE]
  list(clinical = validation, expr = mat)
}

get_translation_proteostasis_genes <- function(signature_genes) {
  msig_genes <- character()
  if (requireNamespace("msigdbr", quietly = TRUE)) {
    msig <- as.data.table(msigdbr::msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME"))
    reactome_keep <- c(
      "REACTOME_TRANSLATION",
      "REACTOME_EUKARYOTIC_TRANSLATION_INITIATION",
      "REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION",
      "REACTOME_RRNA_PROCESSING",
      "REACTOME_RIBOSOMAL_SCANNING_AND_START_CODON_RECOGNITION",
      "REACTOME_FORMATION_OF_A_POOL_OF_FREE_40S_SUBUNITS",
      "REACTOME_NONSENSE_MEDIATED_DECAY_NMD",
      "REACTOME_PROTEIN_FOLDING",
      "REACTOME_CHAPERONIN_MEDIATED_PROTEIN_FOLDING",
      "REACTOME_ACTIVATION_OF_CHAPERONE_GENES_BY_XBP1S",
      "REACTOME_REGULATION_OF_HSF1_MEDIATED_HEAT_SHOCK_RESPONSE"
    )
    msig_genes <- unique(msig[gs_name %in% reactome_keep, gene_symbol])
  }
  rule_genes <- signature_genes[
    grepl("^(RPL|RPS|MRPL|MRPS|EEF|EIF|HSP|DNAJ|CCT|PSM|RPN)", signature_genes) |
      signature_genes %in% c("CALR", "PPIA", "NACA", "TPT1", "UBC", "UBB", "RHEB", "RPS6")
  ]
  unique(c(msig_genes, rule_genes))
}

cat("Reading signature and bulk matrices...\n")
signature <- fread(signature_file)
up_genes <- signature[direction == "up", gene]
down_genes <- signature[direction == "down", gene]
all_signature_genes <- unique(c(up_genes, down_genes))
translation_proteostasis_genes <- get_translation_proteostasis_genes(all_signature_genes)

variants <- list(
  original_up_minus_down = list(up = up_genes, down = down_genes),
  up_only = list(up = up_genes, down = character()),
  no_translation_proteostasis = list(
    up = setdiff(up_genes, translation_proteostasis_genes),
    down = setdiff(down_genes, translation_proteostasis_genes)
  ),
  no_translation_proteostasis_up_only = list(
    up = setdiff(up_genes, translation_proteostasis_genes),
    down = character()
  )
)

variant_gene_counts <- rbindlist(lapply(names(variants), function(v) {
  data.table(
    variant = v,
    up_genes = length(variants[[v]]$up),
    down_genes = length(variants[[v]]$down),
    removed_translation_proteostasis = length(intersect(all_signature_genes, translation_proteostasis_genes))
  )
}))
write_table(variant_gene_counts, "lap3_substate_signature_sensitivity_gene_counts.csv")

variant_gene_table <- rbindlist(lapply(names(variants), function(v) {
  rbindlist(list(
    data.table(variant = v, direction = "up", gene = variants[[v]]$up),
    data.table(variant = v, direction = "down", gene = variants[[v]]$down)
  ))
}))
variant_gene_table[, removed_translation_proteostasis := gene %in% translation_proteostasis_genes]
write_table(variant_gene_table, "lap3_substate_signature_sensitivity_gene_lists.csv")

project_variants <- function(dataset_name, clinical, expr) {
  expr_log2 <- log2(expr + 1)
  out <- copy(clinical)
  sample_col <- if ("barcode" %in% names(out)) "barcode" else "sample_id"
  out[, dataset := dataset_name]
  out[, sample_key := get(sample_col)]
  out[, LAP3_log2_expr_projection := if ("LAP3" %in% rownames(expr_log2)) as.numeric(expr_log2["LAP3", match(sample_key, colnames(expr_log2))]) else NA_real_]
  for (v in names(variants)) {
    score <- score_signature(expr_log2, variants[[v]]$up, variants[[v]]$down)
    out[, paste0("signature_", v) := score[match(sample_key, colnames(expr_log2))]]
  }
  out
}

tcga <- read_tcga_expr()
cgga693 <- read_cgga_693()
cgga325 <- read_cgga_325()

tcga_proj <- project_variants("TCGA", tcga$clinical, tcga$expr)
tcga_path <- fread(tcga_pathway_file)
tcga_path[, sample_key := barcode]
pathways <- c("HALLMARK_MTORC1_SIGNALING", "LEUCINE_BCAA_CORE", "MTORC1_READOUT_CORE", "REACTOME_TRANSLATION")
tcga_proj <- merge(
  tcga_proj,
  tcga_path[, c("sample_key", intersect(pathways, names(tcga_path))), with = FALSE],
  by = "sample_key",
  all.x = TRUE
)
cgga_proj <- rbindlist(list(
  project_variants("CGGA_mRNAseq_693", cgga693$clinical, cgga693$expr),
  project_variants("CGGA_mRNAseq_325", cgga325$clinical, cgga325$expr)
), fill = TRUE)

write_table(tcga_proj, "tcga_lap3_substate_signature_sensitivity_projection.csv")
write_table(cgga_proj, "cgga_lap3_substate_signature_sensitivity_projection.csv")

variant_cols <- paste0("signature_", names(variants))

tcga_cor <- rbindlist(lapply(c("all", sort(unique(tcga_proj$cohort))), function(group) {
  d <- if (group == "all") tcga_proj else tcga_proj[cohort == group]
  rbindlist(lapply(variant_cols, function(score_col) {
    rbindlist(lapply(c("LAP3_log2_expr_projection", pathways), function(v) {
      if (!v %in% names(d)) return(NULL)
      res <- spearman_safe(d[[score_col]], d[[v]], min_n = 30L)
      data.table(dataset = "TCGA", group = group, variant = sub("^signature_", "", score_col), variable = v, n = res$n, spearman_rho = res$rho, p_value = res$p)
    }), fill = TRUE)
  }), fill = TRUE)
}), fill = TRUE)
tcga_cor[, p_adj_BH := p.adjust(p_value, method = "BH")]
write_table(tcga_cor, "tcga_lap3_substate_signature_sensitivity_correlations.csv")

cgga_cor <- rbindlist(lapply(sort(unique(cgga_proj$dataset)), function(ds) {
  d0 <- cgga_proj[dataset == ds]
  rbindlist(lapply(c("all", sort(unique(as.character(d0$tumor_class)))), function(group) {
    d <- if (group == "all") d0 else d0[as.character(tumor_class) == group]
    rbindlist(lapply(variant_cols, function(score_col) {
      rbindlist(lapply(c("LAP3_log2_expr_projection", pathways), function(v) {
        if (!v %in% names(d)) return(NULL)
        res <- spearman_safe(d[[score_col]], d[[v]], min_n = 25L)
        data.table(dataset = ds, group = group, variant = sub("^signature_", "", score_col), variable = v, n = res$n, spearman_rho = res$rho, p_value = res$p)
      }), fill = TRUE)
    }), fill = TRUE)
  }), fill = TRUE)
}), fill = TRUE)
cgga_cor[, p_adj_BH := p.adjust(p_value, method = "BH")]
write_table(cgga_cor, "cgga_lap3_substate_signature_sensitivity_correlations.csv")

tcga_lm <- rbindlist(lapply(variant_cols, function(score_col) {
  rbindlist(lapply(pathways, function(y) {
    if (!y %in% names(tcga_proj)) return(NULL)
    lm_term(tcga_proj, y = y, x = score_col, covariates = c("cohort", "age_years", "idh_status", "codel_1p19q"), model_name = "TCGA_adjusted_tumor_class_molecular")
  }), fill = TRUE)[, variant := sub("^signature_", "", score_col)]
}), fill = TRUE)
tcga_lm[, p_adj_BH := p.adjust(p_value, method = "BH")]
setcolorder(tcga_lm, c("model", "variant", "outcome", "term", "n", "beta", "std_error", "p_value", "p_adj_BH"))
write_table(tcga_lm, "tcga_lap3_substate_signature_sensitivity_adjusted_lm.csv")

cgga_lm <- rbindlist(lapply(sort(unique(cgga_proj$dataset)), function(ds) {
  d <- cgga_proj[dataset == ds]
  rbindlist(lapply(variant_cols, function(score_col) {
    rbindlist(lapply(pathways, function(y) {
      if (!y %in% names(d)) return(NULL)
      lm_term(d, y = y, x = score_col, covariates = c("tumor_class", "age_years", "idh_status", "codel_1p19q"), model_name = paste0(ds, "_adjusted_tumor_class_molecular"))
    }), fill = TRUE)[, variant := sub("^signature_", "", score_col)]
  }), fill = TRUE)
}), fill = TRUE)
cgga_lm[, p_adj_BH := p.adjust(p_value, method = "BH")]
setcolorder(cgga_lm, c("model", "variant", "outcome", "term", "n", "beta", "std_error", "p_value", "p_adj_BH"))
write_table(cgga_lm, "cgga_lap3_substate_signature_sensitivity_adjusted_lm.csv")

summary_rows <- rbindlist(list(
  tcga_lm[, .(dataset = "TCGA", variant, outcome, n, beta, p_adj_BH)],
  cgga_lm[, .(dataset = sub("_adjusted_tumor_class_molecular$", "", model), variant, outcome, n, beta, p_adj_BH)]
), fill = TRUE)
summary_rows[, direction := fifelse(is.na(beta), "NA", fifelse(beta > 0, "positive", "negative"))]
write_table(summary_rows, "lap3_substate_signature_sensitivity_summary.csv")

readme_add <- c(
  "",
  "## Signature Sensitivity",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "Sensitivity variants test whether the bulk bridge is driven mainly by the translation/proteostasis-heavy up signature.",
  "",
  "- `original_up_minus_down`: original 60-up/20-down projection.",
  "- `up_only`: original up arm only.",
  "- `no_translation_proteostasis`: remove Reactome/rule-based translation, ribosome, chaperone, proteasome, and proteostasis genes from up/down arms.",
  "- `no_translation_proteostasis_up_only`: deconfounded up arm only.",
  "",
  "Key outputs:",
  "",
  "- `tables/lap3_substate_signature_sensitivity_gene_counts.csv`",
  "- `tables/lap3_substate_signature_sensitivity_gene_lists.csv`",
  "- `tables/lap3_substate_signature_sensitivity_summary.csv`",
  "- `tables/tcga_lap3_substate_signature_sensitivity_adjusted_lm.csv`",
  "- `tables/cgga_lap3_substate_signature_sensitivity_adjusted_lm.csv`"
)
cat(readme_add, sep = "\n", file = file.path(out_dir, "README.md"), append = TRUE)

cat("Finished:", format(Sys.time()), "\n")
