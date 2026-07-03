#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

project_dir <- "/home/lzb/glioma"
setwd(project_dir)
data.table::setDTthreads(8)
set.seed(20260630)

state_dir <- "Data_Bulk_TCGA/Data_Merged/results/LAP3_State_Module"
bench_dir <- "Data_Bulk_TCGA/Data_Merged/results/LAP3_State_Specificity_Benchmark"
state_table_dir <- file.path(state_dir, "tables")
state_export_dir <- file.path(state_dir, "exports")
state_plot_dir <- file.path(state_dir, "plots")
bench_table_dir <- file.path(bench_dir, "tables")
bench_plot_dir <- file.path(bench_dir, "plots")
log_dir <- file.path(state_dir, "logs")

for (d in c(state_table_dir, state_export_dir, state_plot_dir, bench_table_dir, bench_plot_dir, log_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

log_file <- file.path(log_dir, "lap3_state_module_benchmark.log")
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message")
on.exit({
  sink(type = "message")
  sink(type = "output")
  close(log_con)
}, add = TRUE)

cat("LAP3 state/module construction and specificity benchmark\n")
cat("Started:", format(Sys.time()), "\n")
cat("R:", R.version.string, "\n")
cat("Working directory:", getwd(), "\n\n")

write_table <- function(x, dir, filename) {
  data.table::fwrite(as.data.table(x), file.path(dir, filename))
}

save_plot <- function(plot, dir, filename, width = 7, height = 5) {
  ggsave(file.path(dir, paste0(filename, ".pdf")), plot, width = width, height = height)
  ggsave(file.path(dir, paste0(filename, ".png")), plot, width = width, height = height, dpi = 300)
}

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

module_score <- function(expr_log2, genes, min_genes = 3L) {
  genes <- intersect(unique(genes), rownames(expr_log2))
  if (length(genes) < min_genes) {
    return(rep(NA_real_, ncol(expr_log2)))
  }
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

group_test <- function(data, score_col, group_col, dataset, group_scope) {
  d <- as.data.table(data)
  d <- d[is.finite(get(score_col)) & !is.na(get(group_col))]
  if (nrow(d) < 20L || uniqueN(d[[group_col]]) < 2L) {
    return(data.table(dataset = dataset, group_scope = group_scope, score = score_col, variable = group_col, n = nrow(d), method = NA_character_, p_value = NA_real_))
  }
  if (uniqueN(d[[group_col]]) == 2L) {
    p <- wilcox.test(d[[score_col]] ~ d[[group_col]])$p.value
    method <- "Wilcoxon"
  } else {
    p <- kruskal.test(d[[score_col]] ~ d[[group_col]])$p.value
    method <- "Kruskal-Wallis"
  }
  data.table(
    dataset = dataset,
    group_scope = group_scope,
    score = score_col,
    variable = group_col,
    n = nrow(d),
    n_groups = uniqueN(d[[group_col]]),
    method = method,
    p_value = p
  )
}

lm_term <- function(data, outcome, term, covariates, model_name, dataset) {
  d <- as.data.table(data)
  required <- unique(c(outcome, term, covariates))
  required <- required[required %in% names(d)]
  d <- d[complete.cases(d[, ..required])]
  if (nrow(d) < 50L || uniqueN(d[[term]]) < 3L) {
    return(data.table(dataset = dataset, model = model_name, outcome = outcome, term = term, n = nrow(d), beta = NA_real_, p_value = NA_real_))
  }
  for (cc in covariates) {
    if (cc %in% names(d) && !is.numeric(d[[cc]])) d[[cc]] <- factor(d[[cc]])
  }
  form <- reformulate(c(term, covariates[covariates %in% names(d)]), response = outcome)
  fit <- lm(form, data = d)
  ct <- summary(fit)$coefficients
  if (!term %in% rownames(ct)) {
    return(data.table(dataset = dataset, model = model_name, outcome = outcome, term = term, n = nrow(d), beta = NA_real_, p_value = NA_real_))
  }
  data.table(
    dataset = dataset,
    model = model_name,
    outcome = outcome,
    term = term,
    n = nrow(d),
    beta = ct[term, "Estimate"],
    std_error = ct[term, "Std. Error"],
    p_value = ct[term, "Pr(>|t|)"],
    adjusted_r2 = summary(fit)$adj.r.squared
  )
}

get_msig <- function(category, subcategory = NULL) {
  if (!requireNamespace("msigdbr", quietly = TRUE)) return(data.table(gs_name = character(), gene_symbol = character()))
  if (is.null(subcategory)) {
    as.data.table(msigdbr::msigdbr(species = "Homo sapiens", category = category))[, .(gs_name, gene_symbol)]
  } else {
    as.data.table(msigdbr::msigdbr(species = "Homo sapiens", category = category, subcategory = subcategory))[, .(gs_name, gene_symbol)]
  }
}

manual_gene_sets <- function() {
  hallmark <- get_msig("H")
  reactome <- get_msig("C2", "CP:REACTOME")
  msig <- rbind(hallmark, reactome, fill = TRUE)
  split_msig <- split(msig$gene_symbol, msig$gs_name)
  list(
    HALLMARK_MTORC1_SIGNALING = split_msig[["HALLMARK_MTORC1_SIGNALING"]],
    HALLMARK_HYPOXIA = split_msig[["HALLMARK_HYPOXIA"]],
    HALLMARK_E2F_TARGETS = split_msig[["HALLMARK_E2F_TARGETS"]],
    HALLMARK_G2M_CHECKPOINT = split_msig[["HALLMARK_G2M_CHECKPOINT"]],
    REACTOME_TRANSLATION = split_msig[["REACTOME_TRANSLATION"]],
    LEUCINE_BCAA_CORE = c(
      "BCAT1", "BCAT2", "BCKDHA", "BCKDHB", "DBT", "DLD", "BCKDK", "PPM1K",
      "SLC7A5", "SLC3A2", "SLC43A1", "SLC43A2", "SLC38A2", "SLC38A9",
      "MTOR", "RPTOR", "RRAGA", "RRAGB", "RRAGC", "RRAGD", "LAMTOR1", "LAMTOR2",
      "LAMTOR3", "LAMTOR4", "LAMTOR5", "RHEB", "RPS6KB1", "EIF4EBP1", "RPS6"
    ),
    MTORC1_READOUT_CORE = c(
      "MTOR", "RPTOR", "RHEB", "AKT1", "TSC1", "TSC2", "RRAGA", "RRAGB", "RRAGC",
      "RRAGD", "LAMTOR1", "LAMTOR2", "LAMTOR3", "LAMTOR4", "LAMTOR5",
      "RPS6KB1", "RPS6KB2", "EIF4EBP1", "EIF4EBP2", "RPS6", "EIF4E"
    ),
    TAM_MYELOID_CORE = c("AIF1", "CD68", "LST1", "CSF1R", "TYROBP", "FCER1G", "SPI1", "ITGAM", "C1QA", "C1QB", "C1QC", "CD163", "MRC1", "MSR1", "TREM2", "APOE", "LILRB4"),
    PROLIFERATION_CORE = c("MKI67", "TOP2A", "PCNA", "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "CDK1", "CCNB1", "CCNB2", "AURKA", "AURKB", "BUB1", "BUB1B", "UBE2C"),
    HYPOXIA_CORE = c("CA9", "VEGFA", "SLC2A1", "LDHA", "ENO1", "PGK1", "BNIP3", "NDRG1", "P4HA1", "EGLN3", "ADM", "ANGPTL4"),
    AMINOPEPTIDASE_FAMILY = c("LAP3", "ANPEP", "ERAP1", "ERAP2", "LNPEP", "NPEPPS", "XPNPEP1", "XPNPEP2", "XPNPEP3", "LTA4H", "RNPEP", "ENPEP")
  )
}

translation_proteostasis_genes <- function(genes) {
  sets <- manual_gene_sets()
  reactome_translation <- sets$REACTOME_TRANSLATION
  prefix <- genes[grepl("^(RPL|RPS|MRPL|MRPS|EEF|EIF|HSP|DNAJ|CCT|PSM|RPN)", genes)]
  unique(c(reactome_translation, prefix, "CALR", "PPIA", "NACA", "TPT1", "UBC", "UBB", "RHEB", "RPS6"))
}

cat("Reading TCGA inputs...\n")
tcga_clin <- as.data.table(readRDS("Data_Bulk_TCGA/Data_Merged/results/Clinical_Field_QC/clinical_glioma_analysis_fields.rds"))
tcga_tpm_raw <- readRDS("Data_Bulk_TCGA/Data_Merged/data_analysis/expr_tpm_glioma_uni.rds")
stopifnot(identical(tcga_clin$barcode, colnames(tcga_tpm_raw)[-1]))
tcga_tpm <- as.matrix(tcga_tpm_raw[, -1, drop = FALSE])
rownames(tcga_tpm) <- rownames(tcga_tpm_raw)
tcga_tpm <- collapse_expr(tcga_tpm)
tcga_log2 <- log2(tcga_tpm + 1)
stopifnot("LAP3" %in% rownames(tcga_log2))

cat("Reading TCGA discovery signature...\n")
deg <- fread("Data_Bulk_TCGA/Data_Merged/results/LAP3_Pathway/tables/deg_cohort_balanced_tertile_adjusted.csv")
deg <- deg[order(-rank_t)]
tcga_up <- deg[gene != "LAP3" & adj.P.Val < 0.05 & logFC > 0][1:min(.N, 150), gene]
if (length(tcga_up) < 50L) stop("Too few TCGA discovery genes")

cat("Reading Core GBmap substate genes...\n")
substate <- fread("Data_Bulk_TCGA/Data_Merged/results/LAP3_Substate_Bulk_Projection/tables/core_lap3_detected_substate_signature_genes.csv")
gbmap_up <- substate[direction == "up" & p_adj_BH < 0.05, gene]
gbmap_up <- setdiff(unique(gbmap_up), "LAP3")

all_candidate <- unique(c(tcga_up, gbmap_up))
tp_genes <- translation_proteostasis_genes(all_candidate)

state_sets <- list(
  LAP3_STATE_TCGA_TOP150 = tcga_up,
  LAP3_STATE_GBMAP_UP = gbmap_up,
  LAP3_STATE_UNION = all_candidate,
  LAP3_STATE_UNION_NO_TRANSLATION_PROTEOSTASIS = setdiff(all_candidate, tp_genes)
)
state_sets <- lapply(state_sets, function(x) sort(unique(setdiff(x, "LAP3"))))

state_gene_table <- rbindlist(lapply(names(state_sets), function(set_name) {
  data.table(
    state_set = set_name,
    gene = state_sets[[set_name]],
    in_tcga_top150 = state_sets[[set_name]] %in% tcga_up,
    in_gbmap_up = state_sets[[set_name]] %in% gbmap_up,
    translation_proteostasis_flag = state_sets[[set_name]] %in% tp_genes
  )
}))
state_gene_table <- merge(
  state_gene_table,
  deg[, .(gene, tcga_logFC = logFC, tcga_t = t, tcga_p = P.Value, tcga_fdr = adj.P.Val, tcga_rank_t = rank_t)],
  by = "gene",
  all.x = TRUE
)
state_gene_table <- merge(
  state_gene_table,
  substate[, .(gene, gbmap_direction = direction, gbmap_median_delta = median_delta, gbmap_p_adj_BH = p_adj_BH)],
  by = "gene",
  all.x = TRUE
)
write_table(state_gene_table, state_table_dir, "lap3_state_frozen_gene_sets.csv")

state_counts <- state_gene_table[, .(
  n_genes = .N,
  n_tcga_top150 = sum(in_tcga_top150),
  n_gbmap_up = sum(in_gbmap_up),
  n_translation_proteostasis = sum(translation_proteostasis_flag)
), by = state_set]
write_table(state_counts, state_table_dir, "lap3_state_gene_set_counts.csv")

cat("State gene set counts:\n")
print(state_counts)

score_dataset <- function(dataset, clinical, expr_log2, sample_col, expression_col_name) {
  out <- as.data.table(clinical)
  out[, dataset := dataset]
  out[, sample_key := get(sample_col)]
  out[, LAP3_log2_expr := single_gene_expr(expr_log2, "LAP3")[match(sample_key, colnames(expr_log2))]]
  setnames(out, "LAP3_log2_expr", expression_col_name)
  for (nm in names(state_sets)) {
    out[, (nm) := module_score(expr_log2, state_sets[[nm]])[match(sample_key, colnames(expr_log2))]]
  }
  sets <- manual_gene_sets()
  benchmark_sets <- sets[c(
    "HALLMARK_MTORC1_SIGNALING", "HALLMARK_HYPOXIA", "HALLMARK_E2F_TARGETS",
    "HALLMARK_G2M_CHECKPOINT", "REACTOME_TRANSLATION", "LEUCINE_BCAA_CORE",
    "MTORC1_READOUT_CORE", "TAM_MYELOID_CORE", "PROLIFERATION_CORE",
    "HYPOXIA_CORE", "AMINOPEPTIDASE_FAMILY"
  )]
  for (nm in names(benchmark_sets)) {
    out[, (nm) := module_score(expr_log2, benchmark_sets[[nm]])[match(sample_key, colnames(expr_log2))]]
  }
  out
}

tcga_projection <- score_dataset("TCGA", tcga_clin, tcga_log2, "barcode", "LAP3_log2_expr")
write_table(tcga_projection, state_table_dir, "tcga_lap3_state_score_projection.csv")

cat("Reading CGGA validation dataset and matrices...\n")
cgga_validation <- as.data.table(readRDS("Data_Bulk_CGGA/results/LAP3_CGGA/exports/cgga_lap3_validation_dataset.rds"))

read_cgga_693 <- function() {
  clinical <- cgga_validation[cohort == "mRNAseq_693"]
  mat <- as.matrix(readRDS("Data_Bulk_CGGA/mRNAseq_693/tpm_data.rds"))
  mat <- collapse_expr(mat)
  common <- intersect(clinical$sample_id, colnames(mat))
  clinical <- clinical[match(common, sample_id)]
  mat <- mat[, clinical$sample_id, drop = FALSE]
  list(clinical = clinical, log2 = log2(mat + 1))
}

read_cgga_325 <- function() {
  clinical <- cgga_validation[cohort == "mRNAseq_325"]
  rsem <- fread("Data_Bulk_CGGA/mRNAseq_325/CGGA.mRNAseq_325.RSEM-genes.20200506.txt", data.table = FALSE, check.names = FALSE)
  mat <- as.matrix(rsem[, -1, drop = FALSE])
  rownames(mat) <- rsem$Gene_Name
  mat <- collapse_expr(mat)
  common <- intersect(clinical$sample_id, colnames(mat))
  clinical <- clinical[match(common, sample_id)]
  mat <- mat[, clinical$sample_id, drop = FALSE]
  list(clinical = clinical, log2 = log2(mat + 1))
}

cgga693 <- read_cgga_693()
cgga325 <- read_cgga_325()
cgga_projection <- rbindlist(list(
  score_dataset("CGGA_mRNAseq_693", cgga693$clinical, cgga693$log2, "sample_id", "LAP3_log2_expr"),
  score_dataset("CGGA_mRNAseq_325", cgga325$clinical, cgga325$log2, "sample_id", "LAP3_log2_expr")
), fill = TRUE)
write_table(cgga_projection, state_table_dir, "cgga_lap3_state_score_projection.csv")

saveRDS(
  list(
    state_sets = state_sets,
    tcga_projection = tcga_projection,
    cgga_projection = cgga_projection
  ),
  file.path(state_export_dir, "lap3_state_module_projection.rds")
)

score_cols <- names(state_sets)
benchmark_cols <- c(
  "HALLMARK_MTORC1_SIGNALING", "HALLMARK_HYPOXIA", "HALLMARK_E2F_TARGETS",
  "HALLMARK_G2M_CHECKPOINT", "REACTOME_TRANSLATION", "LEUCINE_BCAA_CORE",
  "MTORC1_READOUT_CORE", "TAM_MYELOID_CORE", "PROLIFERATION_CORE",
  "HYPOXIA_CORE", "AMINOPEPTIDASE_FAMILY"
)

cor_table <- rbindlist(list(
  rbindlist(lapply(c("all", sort(unique(tcga_projection$cohort))), function(group) {
    d <- if (group == "all") tcga_projection else tcga_projection[cohort == group]
    rbindlist(lapply(score_cols, function(score_col) {
      rbindlist(lapply(c("LAP3_log2_expr", benchmark_cols), function(v) {
        if (!v %in% names(d)) return(NULL)
        res <- spearman_safe(d[[score_col]], d[[v]], min_n = 30L)
        data.table(dataset = "TCGA", group = group, score = score_col, variable = v, n = res$n, spearman_rho = res$rho, p_value = res$p)
      }))
    }))
  })),
  rbindlist(lapply(c("all", sort(unique(cgga_projection$dataset))), function(ds) {
    d0 <- if (ds == "all") cgga_projection else cgga_projection[dataset == ds]
    rbindlist(lapply(c("all", sort(unique(as.character(d0$tumor_class)))), function(group) {
      d <- if (group == "all") d0 else d0[as.character(tumor_class) == group]
      rbindlist(lapply(score_cols, function(score_col) {
        rbindlist(lapply(c("LAP3_log2_expr", benchmark_cols), function(v) {
          if (!v %in% names(d)) return(NULL)
          res <- spearman_safe(d[[score_col]], d[[v]], min_n = 25L)
          data.table(dataset = ds, group = group, score = score_col, variable = v, n = res$n, spearman_rho = res$rho, p_value = res$p)
        }))
      }))
    }))
  }))
), fill = TRUE)
cor_table[, p_adj_BH := p.adjust(p_value, method = "BH")]
write_table(cor_table, bench_table_dir, "lap3_state_benchmark_score_correlations.csv")

clinical_tests <- rbindlist(list(
  rbindlist(lapply(score_cols, function(score) {
    rbindlist(lapply(c("cohort", "grade", "idh_status", "codel_1p19q", "mgmt_status"), function(v) {
      if (!v %in% names(tcga_projection)) return(NULL)
      group_test(tcga_projection, score, v, "TCGA", "all")
    }), fill = TRUE)
  }), fill = TRUE),
  rbindlist(lapply(score_cols, function(score) {
    rbindlist(lapply(c("dataset", "tumor_class", "grade", "idh_status", "codel_1p19q", "mgmt_status"), function(v) {
      if (!v %in% names(cgga_projection)) return(NULL)
      group_test(cgga_projection, score, v, "CGGA", "all")
    }), fill = TRUE)
  }), fill = TRUE)
), fill = TRUE)
clinical_tests[, p_adj_BH := p.adjust(p_value, method = "BH")]
write_table(clinical_tests, bench_table_dir, "lap3_state_clinical_group_tests.csv")

model_table <- rbindlist(list(
  rbindlist(lapply(score_cols, function(score) {
    rbindlist(list(
      lm_term(tcga_projection, score, "LAP3_log2_expr", c("cohort", "grade", "idh_status", "codel_1p19q", "age_years"), "clinical_covariates", "TCGA"),
      lm_term(tcga_projection, score, "LAP3_log2_expr", c("cohort", "grade", "idh_status", "codel_1p19q", "HALLMARK_HYPOXIA", "TAM_MYELOID_CORE", "PROLIFERATION_CORE", "HALLMARK_MTORC1_SIGNALING"), "clinical_plus_benchmark_scores", "TCGA")
    ), fill = TRUE)
  }), fill = TRUE),
  rbindlist(lapply(score_cols, function(score) {
    rbindlist(list(
      lm_term(cgga_projection, score, "LAP3_log2_expr", c("dataset", "tumor_class", "grade", "idh_status", "codel_1p19q", "age_years"), "clinical_covariates", "CGGA"),
      lm_term(cgga_projection, score, "LAP3_log2_expr", c("dataset", "tumor_class", "grade", "idh_status", "codel_1p19q", "HALLMARK_HYPOXIA", "TAM_MYELOID_CORE", "PROLIFERATION_CORE", "HALLMARK_MTORC1_SIGNALING"), "clinical_plus_benchmark_scores", "CGGA")
    ), fill = TRUE)
  }), fill = TRUE)
), fill = TRUE)
model_table[, p_adj_BH := p.adjust(p_value, method = "BH")]
write_table(model_table, bench_table_dir, "lap3_state_adjusted_lm_models.csv")

sets <- manual_gene_sets()
amino_genes <- unique(sets$AMINOPEPTIDASE_FAMILY)
amino_table <- rbindlist(lapply(c("TCGA", "CGGA_mRNAseq_693", "CGGA_mRNAseq_325"), function(ds) {
  if (ds == "TCGA") {
    expr <- tcga_log2
    state <- tcga_projection$LAP3_STATE_UNION[match(colnames(expr), tcga_projection$sample_key)]
  } else if (ds == "CGGA_mRNAseq_693") {
    expr <- cgga693$log2
    state <- cgga_projection[dataset == ds]$LAP3_STATE_UNION[match(colnames(expr), cgga_projection[dataset == ds]$sample_key)]
  } else {
    expr <- cgga325$log2
    state <- cgga_projection[dataset == ds]$LAP3_STATE_UNION[match(colnames(expr), cgga_projection[dataset == ds]$sample_key)]
  }
  rbindlist(lapply(amino_genes, function(g) {
    x <- single_gene_expr(expr, g)
    res <- spearman_safe(x, state, min_n = 25L)
    data.table(dataset = ds, gene = g, present = g %in% rownames(expr), n = res$n, spearman_rho_with_LAP3_STATE_UNION = res$rho, p_value = res$p)
  }))
}), fill = TRUE)
amino_table[, p_adj_BH := p.adjust(p_value, method = "BH"), by = dataset]
amino_table[, rank_abs_rho := frank(-abs(spearman_rho_with_LAP3_STATE_UNION), ties.method = "min"), by = dataset]
write_table(amino_table, bench_table_dir, "aminopeptidase_gene_benchmark.csv")

cat("Running expression-matched random gene benchmark...\n")
mean_expr <- rowMeans(tcga_log2, na.rm = TRUE)
lap3_mean <- mean_expr[["LAP3"]]
pool <- names(mean_expr)[is.finite(mean_expr) & names(mean_expr) != "LAP3"]
pool <- pool[!pool %in% unlist(state_sets)]
pool <- pool[order(abs(mean_expr[pool] - lap3_mean))]
random_pool <- head(pool, min(1000L, length(pool)))
state_tcga <- tcga_projection$LAP3_STATE_UNION[match(colnames(tcga_log2), tcga_projection$sample_key)]
random_gene_cor <- rbindlist(lapply(c("LAP3", random_pool), function(g) {
  res <- spearman_safe(single_gene_expr(tcga_log2, g), state_tcga, min_n = 50L)
  data.table(gene = g, mean_log2_expr = mean_expr[[g]], n = res$n, spearman_rho = res$rho, p_value = res$p, is_lap3 = g == "LAP3")
}))
random_gene_cor[, rank_abs_rho := frank(-abs(spearman_rho), ties.method = "min")]
random_gene_cor[, percentile_abs_rho := 1 - (rank_abs_rho - 1) / .N]
write_table(random_gene_cor, bench_table_dir, "expression_matched_random_gene_benchmark_tcga.csv")

summary_table <- data.table(
  metric = c(
    "frozen_primary_score",
    "primary_score_n_genes",
    "tcga_lap3_state_lap3_rho",
    "cgga693_lap3_state_lap3_rho",
    "cgga325_lap3_state_lap3_rho",
    "tcga_lap3_random_rank_abs_rho",
    "tcga_lap3_random_percentile_abs_rho"
  ),
  value = c(
    "LAP3_STATE_UNION",
    as.character(length(state_sets$LAP3_STATE_UNION)),
    as.character(cor_table[dataset == "TCGA" & group == "all" & score == "LAP3_STATE_UNION" & variable == "LAP3_log2_expr", spearman_rho][1]),
    as.character(cor_table[dataset == "CGGA_mRNAseq_693" & group == "all" & score == "LAP3_STATE_UNION" & variable == "LAP3_log2_expr", spearman_rho][1]),
    as.character(cor_table[dataset == "CGGA_mRNAseq_325" & group == "all" & score == "LAP3_STATE_UNION" & variable == "LAP3_log2_expr", spearman_rho][1]),
    as.character(random_gene_cor[is_lap3 == TRUE, rank_abs_rho][1]),
    as.character(random_gene_cor[is_lap3 == TRUE, percentile_abs_rho][1])
  )
)
write_table(summary_table, state_table_dir, "lap3_state_first_pass_summary.csv")

get_rho <- function(dataset_name, variable_name, score_name = "LAP3_STATE_UNION", group_name = "all") {
  x <- cor_table[dataset == dataset_name & group == group_name & score == score_name & variable == variable_name, spearman_rho]
  if (length(x) < 1L) return(NA_real_)
  x[1]
}

interpretation_summary <- data.table(
  item = c(
    "primary_score",
    "primary_score_genes",
    "tcga_lap3_anchor_rho",
    "cgga693_lap3_anchor_rho",
    "cgga325_lap3_anchor_rho",
    "tcga_tam_rho",
    "tcga_hypoxia_rho",
    "tcga_mtorc1_rho",
    "tcga_bcaa_rho",
    "tcga_translation_rho",
    "tcga_adjusted_lap3_beta_after_benchmarks",
    "cgga_adjusted_lap3_beta_after_benchmarks",
    "random_gene_benchmark_rank",
    "interpretation"
  ),
  value = c(
    "LAP3_STATE_UNION",
    as.character(length(state_sets$LAP3_STATE_UNION)),
    sprintf("%.3f", get_rho("TCGA", "LAP3_log2_expr")),
    sprintf("%.3f", get_rho("CGGA_mRNAseq_693", "LAP3_log2_expr")),
    sprintf("%.3f", get_rho("CGGA_mRNAseq_325", "LAP3_log2_expr")),
    sprintf("%.3f", get_rho("TCGA", "TAM_MYELOID_CORE")),
    sprintf("%.3f", get_rho("TCGA", "HALLMARK_HYPOXIA")),
    sprintf("%.3f", get_rho("TCGA", "HALLMARK_MTORC1_SIGNALING")),
    sprintf("%.3f", get_rho("TCGA", "LEUCINE_BCAA_CORE")),
    sprintf("%.3f", get_rho("TCGA", "REACTOME_TRANSLATION")),
    sprintf("%.3f", model_table[dataset == "TCGA" & model == "clinical_plus_benchmark_scores" & outcome == "LAP3_STATE_UNION", beta][1]),
    sprintf("%.3f", model_table[dataset == "CGGA" & model == "clinical_plus_benchmark_scores" & outcome == "LAP3_STATE_UNION", beta][1]),
    paste0(random_gene_cor[is_lap3 == TRUE, rank_abs_rho][1], "/", nrow(random_gene_cor)),
    "Reproducible LAP3-centered malignant-microenvironment/anabolic state; strongly coupled to TAM/hypoxia/mTORC1/BCAA axes, therefore not a clean malignant-intrinsic or causal mTORC1 mechanism."
  )
)
write_table(interpretation_summary, state_table_dir, "lap3_state_interpretation_summary.csv")

p_cor <- cor_table[
  score == "LAP3_STATE_UNION" &
    dataset %in% c("TCGA", "CGGA_mRNAseq_693", "CGGA_mRNAseq_325") &
    group == "all" &
    variable %in% c("LAP3_log2_expr", "HALLMARK_MTORC1_SIGNALING", "LEUCINE_BCAA_CORE", "REACTOME_TRANSLATION", "TAM_MYELOID_CORE", "HALLMARK_HYPOXIA", "PROLIFERATION_CORE")
]
p_cor[, variable := factor(variable, levels = unique(variable))]
p <- ggplot(p_cor, aes(variable, spearman_rho, fill = dataset)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  coord_flip() +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  labs(x = NULL, y = "Spearman rho with LAP3_STATE_UNION", fill = NULL) +
  theme_bw(base_size = 10)
save_plot(p, bench_plot_dir, "lap3_state_benchmark_correlations", width = 8, height = 4.8)

p_rand <- random_gene_cor[order(rank_abs_rho)][1:min(.N, 80)]
p_rand[, label := ifelse(is_lap3, "LAP3", "expression-matched genes")]
p2 <- ggplot(p_rand, aes(rank_abs_rho, abs(spearman_rho), color = is_lap3)) +
  geom_point(alpha = 0.75, size = 1.3) +
  scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "#C23B23"), labels = c("Matched genes", "LAP3")) +
  labs(x = "Rank by absolute correlation", y = "|Spearman rho| with LAP3_STATE_UNION", color = NULL) +
  theme_bw(base_size = 10)
save_plot(p2, bench_plot_dir, "lap3_expression_matched_random_gene_rank", width = 6, height = 4.3)

readme_state <- c(
  "# LAP3 State Module",
  "",
  paste0("Generated: ", format(Sys.time())),
  "",
  "## Purpose",
  "",
  "Build a first-pass candidate LAP3-centered state/module score for the upgraded manuscript route.",
  "",
  "## Primary Score",
  "",
  "`LAP3_STATE_UNION` is the primary first-pass candidate score. It combines:",
  "",
  "- TCGA cohort-balanced LAP3-high vs low top 150 up-regulated genes;",
  "- Core GBmap LAP3-detected/high substate up genes passing BH FDR < 0.05;",
  "- LAP3 itself is excluded.",
  "",
  "## Key Outputs",
  "",
  "- `tables/lap3_state_frozen_gene_sets.csv`",
  "- `tables/lap3_state_gene_set_counts.csv`",
  "- `tables/lap3_state_interpretation_summary.csv`",
  "- `tables/tcga_lap3_state_score_projection.csv`",
  "- `tables/cgga_lap3_state_score_projection.csv`",
  "- `exports/lap3_state_module_projection.rds`",
  "",
  "## First-Pass Summary",
  "",
  paste(capture.output(print(summary_table)), collapse = "\n"),
  "",
  "## Interpretation Summary",
  "",
  paste(capture.output(print(interpretation_summary)), collapse = "\n"),
  "",
  "## Interpretation Boundary",
  "",
  "This score is a transcriptional state score. It supports a LAP3-centered malignant-microenvironment/anabolic state, not LAP3 enzymatic causality, intracellular leucine flux, or phospho-mTORC1 activation."
)
writeLines(readme_state, file.path(state_dir, "README.md"))

readme_bench <- c(
  "# LAP3 State Specificity Benchmark",
  "",
  paste0("Generated: ", format(Sys.time())),
  "",
  "## Purpose",
  "",
  "Test whether the first-pass LAP3-state score behaves as a reproducible state rather than a simple unexamined marker.",
  "",
  "## Benchmarks",
  "",
  "- aminopeptidase genes;",
  "- hypoxia scores;",
  "- TAM/myeloid scores;",
  "- proliferation scores;",
  "- mTORC1/BCAA/translation/readout scores;",
  "- grade/IDH/1p19q/MGMT clinical group tests;",
  "- expression-matched random genes in TCGA.",
  "",
  "## Key Outputs",
  "",
  "- `tables/lap3_state_benchmark_score_correlations.csv`",
  "- `tables/lap3_state_clinical_group_tests.csv`",
  "- `tables/lap3_state_adjusted_lm_models.csv`",
  "- `tables/aminopeptidase_gene_benchmark.csv`",
  "- `tables/expression_matched_random_gene_benchmark_tcga.csv`",
  "- `plots/lap3_state_benchmark_correlations.*`",
  "- `plots/lap3_expression_matched_random_gene_rank.*`",
  "",
  "## Boundary",
  "",
  "The benchmark is diagnostic and hypothesis-generating. Because the state is discovered from LAP3-associated transcriptional structure, the expression-matched random benchmark should be read as an anchoring/specificity check, not as independent causal validation."
  ,
  "",
  "First-pass result: the state is highly reproducible across TCGA and CGGA, but it is also strongly coupled to TAM/myeloid, hypoxia, mTORC1, BCAA and translation scores. This favors a state/ecosystem framing over a single clean pathway or malignant-cell-intrinsic mechanism."
)
writeLines(readme_bench, file.path(bench_dir, "README.md"))

cat("\nCompleted:", format(Sys.time()), "\n")
