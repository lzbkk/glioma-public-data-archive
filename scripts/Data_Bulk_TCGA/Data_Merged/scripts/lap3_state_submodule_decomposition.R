#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})

project_dir <- "/home/lzb/glioma"
setwd(project_dir)
data.table::setDTthreads(16)
set.seed(20260701)

out_dir <- "Data_Bulk_TCGA/Data_Merged/results/LAP3_State_Submodules"
table_dir <- file.path(out_dir, "tables")
source_dir <- file.path(out_dir, "source_data")
export_dir <- file.path(out_dir, "exports")
log_dir <- file.path(out_dir, "logs")
for (d in c(table_dir, source_dir, export_dir, log_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

log_file <- file.path(log_dir, "lap3_state_submodule_decomposition.log")
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message")
on.exit({
  sink(type = "message")
  sink(type = "output")
  close(log_con)
}, add = TRUE)

write_table <- function(x, filename) {
  fwrite(as.data.table(x), file.path(table_dir, filename))
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

spearman_safe <- function(x, y, min_n = 20L) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) < min_n || length(unique(x)) < 3L || length(unique(y)) < 3L) {
    return(list(n = length(x), rho = NA_real_, p_value = NA_real_))
  }
  test <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  list(n = length(x), rho = unname(test$estimate), p_value = test$p.value)
}

finite_min <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else min(x)
}

finite_max <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else max(x)
}

cluster_bootstrap <- function(data, cluster, statistic, replicates = 500L, seed = 20260701L) {
  clusters <- unique(as.character(data[[cluster]]))
  clusters <- clusters[!is.na(clusters)]
  observed <- statistic(data)
  if (length(clusters) < 2L || !is.finite(observed)) {
    return(list(estimate = observed, ci_low = NA_real_, ci_high = NA_real_))
  }
  cluster_index <- split(seq_len(nrow(data)), as.character(data[[cluster]]))
  cluster_index <- cluster_index[names(cluster_index) %in% clusters]
  set.seed(seed)
  vals <- replicate(replicates, {
    sampled <- sample(clusters, length(clusters), replace = TRUE)
    idx <- unlist(cluster_index[sampled], use.names = FALSE)
    suppressWarnings(statistic(data[idx]))
  })
  vals <- vals[is.finite(vals)]
  if (length(vals) < 50L) {
    return(list(estimate = observed, ci_low = NA_real_, ci_high = NA_real_))
  }
  ci <- quantile(vals, c(0.025, 0.975), names = FALSE, na.rm = TRUE)
  list(estimate = observed, ci_low = ci[[1]], ci_high = ci[[2]])
}

leave_one_cluster_out <- function(data, cluster, statistic) {
  clusters <- unique(as.character(data[[cluster]]))
  clusters <- clusters[!is.na(clusters)]
  data.table(
    omitted_cluster = clusters,
    estimate = vapply(clusters, function(cl) {
      suppressWarnings(statistic(data[as.character(data[[cluster]]) != cl]))
    }, numeric(1))
  )
}

manual_gene_sets <- function() {
  list(
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
    TAM_MYELOID_CORE = c(
      "AIF1", "CD68", "LST1", "CSF1R", "TYROBP", "FCER1G", "SPI1", "ITGAM",
      "C1QA", "C1QB", "C1QC", "CD163", "MRC1", "MSR1", "TREM2", "APOE",
      "LILRB4", "FCGR3A", "FCGR2A", "FCGR2B", "MYD88", "LILRB1", "LILRB2"
    ),
    HYPOXIA_CORE = c(
      "CA9", "VEGFA", "SLC2A1", "LDHA", "ENO1", "PGK1", "BNIP3", "NDRG1",
      "P4HA1", "EGLN3", "ADM", "ANGPTL4", "BHLHE40", "HK2", "ALDOA", "GAPDH",
      "TPI1"
    ),
    PROLIFERATION_CORE = c(
      "MKI67", "TOP2A", "PCNA", "MCM2", "MCM3", "MCM4", "MCM5", "MCM6",
      "CDK1", "CCNB1", "CCNB2", "AURKA", "AURKB", "BUB1", "BUB1B", "UBE2C"
    )
  )
}

assign_submodule <- function(gene_table) {
  sets <- manual_gene_sets()
  dt <- copy(gene_table)
  dt <- dt[state_set == "LAP3_STATE_UNION"]
  dt[, gene_upper := toupper(gene)]

  translation_regex <- "^(RPL|RPS|MRPL|MRPS|EEF|EIF)"
  proteostasis_regex <- "^(HSP|DNAJ|CCT|PSM|RPN|UBB|UBC)"
  immune_regex <- "^(HLA|IFI|IFIT|IFITM|ISG|OAS|GBP|IRF|STAT|CXCL|CCL|LILR|FCGR|C1Q|C1R|C1S)"

  dt[, flag_translation := grepl(translation_regex, gene_upper) |
       gene_upper %in% c("NACA", "TPT1", sets$MTORC1_READOUT_CORE, sets$LEUCINE_BCAA_CORE)]
  dt[, flag_proteostasis := grepl(proteostasis_regex, gene_upper) |
       gene_upper %in% c("CALR", "PPIA", "HSPD1", "HSP90B1", "SELENON")]
  dt[, flag_myeloid_tam := gene_upper %in% sets$TAM_MYELOID_CORE |
       grepl(immune_regex, gene_upper) |
       gene_upper %in% c("B2M", "PARP9", "PLSCR1", "NMI", "APOBEC3C", "IFI16", "DTX3L", "CMTM6")]
  dt[, flag_hypoxia_perinecrotic := gene_upper %in% sets$HYPOXIA_CORE]
  dt[, flag_malignant_state := in_gbmap_up | gene_upper %in% c(
    "FABP7", "S100B", "GPM6B", "BCAS1", "DLL1", "TUBA1A", "TIMP1", "CLU",
    "ANXA2", "ANXA5", "DBI", "PRDX1", "CST3"
  )]

  dt[, primary_submodule := fifelse(
    flag_myeloid_tam,
    "LAP3_MYELOID_TAM_CONTEXT_MODULE",
    fifelse(
      flag_hypoxia_perinecrotic,
      "LAP3_HYPOXIA_PERINECROTIC_MODULE",
      fifelse(
        flag_proteostasis,
        "LAP3_PROTEOSTASIS_STRESS_MODULE",
        fifelse(
          flag_translation,
          "LAP3_ANABOLIC_TRANSLATION_MODULE",
          "LAP3_MALIGNANT_STATE_MODULE"
        )
      )
    )
  )]

  dt[, assignment_rule := fifelse(
    primary_submodule == "LAP3_MYELOID_TAM_CONTEXT_MODULE", "immune/myeloid marker or interferon/HLA/complement pattern",
    fifelse(
      primary_submodule == "LAP3_HYPOXIA_PERINECROTIC_MODULE", "hypoxia/glycolysis/perinecrotic marker",
      fifelse(
        primary_submodule == "LAP3_PROTEOSTASIS_STRESS_MODULE", "proteostasis/chaperone/proteasome/ER-stress marker",
        fifelse(
          primary_submodule == "LAP3_ANABOLIC_TRANSLATION_MODULE", "translation/mTORC1/amino-acid anabolic marker",
          "GBmap malignant-state/substate or residual state-associated gene"
        )
      )
    )
  )]

  dt[, c("gene_upper") := NULL]
  setcolorder(dt, c(
    "gene", "primary_submodule", "assignment_rule",
    "flag_translation", "flag_proteostasis", "flag_myeloid_tam",
    "flag_hypoxia_perinecrotic", "flag_malignant_state"
  ))
  dt[order(primary_submodule, gene)]
}

score_bulk_projection <- function(dataset_label, clinical, expr_log2, sample_col, submodule_sets) {
  out <- copy(as.data.table(clinical))
  out[, dataset := dataset_label]
  out[, sample_key := get(sample_col)]
  out[, LAP3_log2_expr := single_gene_expr(expr_log2, "LAP3")[match(sample_key, colnames(expr_log2))]]
  for (nm in names(submodule_sets)) {
    out[, (nm) := module_score(expr_log2, submodule_sets[[nm]])[match(sample_key, colnames(expr_log2))]]
  }
  out
}

summarise_bulk_correlations <- function(proj, score_cols, dataset_label) {
  benchmark_cols <- c(
    "LAP3_log2_expr", "LAP3_STATE_UNION", "LAP3_STATE_UNION_NO_TRANSLATION_PROTEOSTASIS",
    "HALLMARK_MTORC1_SIGNALING", "LEUCINE_BCAA_CORE", "MTORC1_READOUT_CORE",
    "REACTOME_TRANSLATION", "TAM_MYELOID_CORE", "HALLMARK_HYPOXIA",
    "PROLIFERATION_CORE", "HYPOXIA_CORE"
  )
  benchmark_cols <- intersect(benchmark_cols, names(proj))
  rbindlist(lapply(c("all", sort(unique(as.character(proj$cohort %||% proj$tumor_class)))), function(group) {
    d <- if (group == "all") proj else {
      if ("cohort" %in% names(proj) && group %in% proj$cohort) proj[as.character(cohort) == group] else proj[as.character(tumor_class) == group]
    }
    rbindlist(lapply(score_cols, function(score_col) {
      rbindlist(lapply(benchmark_cols, function(v) {
        res <- spearman_safe(d[[score_col]], d[[v]], min_n = 25L)
        data.table(dataset = dataset_label, group = group, submodule = score_col, variable = v, n = res$n, spearman_rho = res$rho, p_value = res$p_value)
      }))
    }))
  }), fill = TRUE)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

run_gbmap_projection <- function(submodule_sets) {
  cache_file <- "Data_scRNA_GEO/GBmap_Core/cache/core_gbmap_lap3_state_union_lightweight.rds"
  if (!file.exists(cache_file)) {
    warning("Missing Core GBmap LAP3-state cache; skipping Core GBmap projection")
    return(NULL)
  }
  cache <- readRDS(cache_file)
  obs <- as.data.table(cache$obs)
  expr <- cache$normalized
  raw <- cache$raw
  stopifnot(nrow(obs) == nrow(expr), "LAP3" %in% colnames(expr))

  map_state <- function(annotation_level_3, annotation_level_4) {
    state <- rep(NA_character_, length(annotation_level_3))
    state[annotation_level_3 == "AC-like" | grepl("^AC-like", annotation_level_4)] <- "AC"
    state[annotation_level_3 == "OPC-like" | grepl("^OPC-like", annotation_level_4)] <- "OPC"
    state[annotation_level_3 == "NPC-like" | grepl("^NPC-like", annotation_level_4)] <- "NPC"
    state[annotation_level_3 == "MES-like" | grepl("^MES-like", annotation_level_4)] <- "MES"
    state
  }

  obs[, author_donor := paste(author, donor_id, sep = "::")]
  obs[, author_state := map_state(annotation_level_3, annotation_level_4)]
  obs[, main_entry := author != "Neftel2019" & annotation_level_1 == "Neoplastic" & !is.na(author_state)]
  obs[, strict_entry := main_entry & iCNV == "aneuploid"]

  analysis_index <- which(obs$main_entry)
  expr_main <- expr[analysis_index, , drop = FALSE]
  raw_main <- raw[analysis_index, , drop = FALSE]
  obs_main <- copy(obs[analysis_index])

  score_dt <- obs_main[, .(
    cell_id = `_index`, author, donor_id, sample, assay, annotation_level_3,
    annotation_level_4, iCNV, author_donor, author_state, strict_entry
  )]
  for (nm in names(submodule_sets)) {
    genes <- intersect(submodule_sets[[nm]], colnames(expr_main))
    score_dt[, (nm) := if (length(genes) >= 3L) Matrix::rowMeans(expr_main[, genes, drop = FALSE]) else NA_real_]
  }
  score_dt[, lap3_log_norm := as.numeric(expr_main[, "LAP3"])]
  score_dt[, lap3_raw := as.numeric(raw_main[, "LAP3"])]
  score_dt[, lap3_detected := lap3_raw > 0]
  score_dt[, entry_variant := "main_neoplastic_exclude_neftel2019"]
  strict_dt <- copy(score_dt[strict_entry == TRUE])
  strict_dt[, entry_variant := "strict_neoplastic_aneuploid_exclude_neftel2019"]
  score_dt <- rbindlist(list(score_dt, strict_dt), use.names = TRUE, fill = TRUE)

  fwrite(score_dt, file.path(source_dir, "gbmap_core_lap3_state_submodule_cell_scores.csv.gz"))

  thresholds <- c(20L, 50L, 100L)
  patient_state <- rbindlist(lapply(thresholds, function(threshold_value) {
    score_dt[
      ,
      c(
        list(
          n_cells = .N,
          lap3_mean = mean(lap3_log_norm),
          lap3_detection_rate = mean(lap3_detected)
        ),
        lapply(.SD, mean, na.rm = TRUE)
      ),
      by = .(entry_variant, author, donor_id, author_donor, author_state),
      .SDcols = names(submodule_sets)
    ][n_cells >= threshold_value][, threshold := threshold_value][]
  }), fill = TRUE)
  write_table(patient_state, "gbmap_core_lap3_state_submodule_patient_state_summary.csv")

  cat("Running Core GBmap donor-state associations with cluster bootstrap...\n")
  association <- rbindlist(lapply(sort(unique(patient_state$entry_variant)), function(variant) {
    rbindlist(lapply(sort(unique(patient_state$threshold)), function(threshold_value) {
      rbindlist(lapply(sort(unique(patient_state$author_state)), function(state) {
        d_state <- patient_state[entry_variant == variant & threshold == threshold_value & author_state == state]
        rbindlist(lapply(names(submodule_sets), function(module) {
          res <- spearman_safe(d_state$lap3_mean, d_state[[module]], min_n = 6L)
          stat_fun <- function(d) suppressWarnings(cor(d$lap3_mean, d[[module]], method = "spearman", use = "complete.obs"))
          boot <- if (res$n >= 6L && uniqueN(d_state$author_donor) >= 3L) cluster_bootstrap(d_state, "author_donor", stat_fun, 500L) else list(ci_low = NA_real_, ci_high = NA_real_)
          lopo <- if (res$n >= 6L && uniqueN(d_state$author) >= 3L) leave_one_cluster_out(d_state, "author", stat_fun) else data.table(estimate = NA_real_)
          data.table(
            entry_variant = variant,
            threshold = threshold_value,
            author_state = state,
            submodule = module,
            n_donor_states = res$n,
            n_authors = uniqueN(d_state$author),
            n_donors = uniqueN(d_state$author_donor),
            spearman_rho = res$rho,
            ci_low = boot$ci_low,
            ci_high = boot$ci_high,
            p_value = res$p_value,
            leave_one_author_min_rho = finite_min(lopo$estimate),
            leave_one_author_max_rho = finite_max(lopo$estimate)
          )
        }))
      }))
    }))
  }), fill = TRUE)
  association[, p_adj_BH := p.adjust(p_value, method = "BH")]
  write_table(association, "gbmap_core_lap3_state_submodule_lap3_associations.csv")

  list(cell_scores = score_dt, patient_state = patient_state, association = association)
}

cat("Started:", format(Sys.time()), "\n")

gene_file <- "Data_Bulk_TCGA/Data_Merged/results/LAP3_State_Module/tables/lap3_state_frozen_gene_sets.csv"
state_gene_table <- fread(gene_file)
assignment <- assign_submodule(state_gene_table)
write_table(assignment, "lap3_state_submodule_gene_assignment.csv")

submodule_sets <- split(assignment$gene, assignment$primary_submodule)
submodule_counts <- assignment[, .(
  n_genes = .N,
  n_tcga_top150 = sum(in_tcga_top150),
  n_gbmap_up = sum(in_gbmap_up),
  n_translation_flag = sum(translation_proteostasis_flag),
  median_tcga_t = median(tcga_t, na.rm = TRUE),
  median_gbmap_delta = median(gbmap_median_delta, na.rm = TRUE)
), by = primary_submodule][order(primary_submodule)]
write_table(submodule_counts, "lap3_state_submodule_gene_counts.csv")

coverage_for_projection <- function(expr_genes, dataset) {
  rbindlist(lapply(names(submodule_sets), function(module) {
    genes <- submodule_sets[[module]]
    present <- intersect(genes, expr_genes)
    data.table(
      dataset = dataset,
      submodule = module,
      requested = length(genes),
      available = length(present),
      coverage = length(present) / length(genes),
      missing_genes = paste(setdiff(genes, present), collapse = ";")
    )
  }))
}

cat("Reading existing bulk projections and raw expression matrices...\n")
state_projection <- readRDS("Data_Bulk_TCGA/Data_Merged/results/LAP3_State_Module/exports/lap3_state_module_projection.rds")
tcga_base <- as.data.table(state_projection$tcga_projection)
cgga_base <- as.data.table(state_projection$cgga_projection)

tcga_tpm_raw <- readRDS("Data_Bulk_TCGA/Data_Merged/data_analysis/expr_tpm_glioma_uni.rds")
tcga_tpm <- as.matrix(tcga_tpm_raw[, -1, drop = FALSE])
rownames(tcga_tpm) <- rownames(tcga_tpm_raw)
tcga_log2 <- log2(collapse_expr(tcga_tpm) + 1)

cgga_validation <- as.data.table(readRDS("Data_Bulk_CGGA/results/LAP3_CGGA/exports/cgga_lap3_validation_dataset.rds"))
clinical693 <- cgga_validation[cohort == "mRNAseq_693"]
mat693 <- as.matrix(readRDS("Data_Bulk_CGGA/mRNAseq_693/tpm_data.rds"))
mat693 <- collapse_expr(mat693)
common693 <- intersect(clinical693$sample_id, colnames(mat693))
clinical693 <- clinical693[match(common693, sample_id)]
mat693 <- mat693[, clinical693$sample_id, drop = FALSE]
log693 <- log2(mat693 + 1)

clinical325 <- cgga_validation[cohort == "mRNAseq_325"]
rsem325 <- fread("Data_Bulk_CGGA/mRNAseq_325/CGGA.mRNAseq_325.RSEM-genes.20200506.txt", data.table = FALSE, check.names = FALSE)
mat325 <- as.matrix(rsem325[, -1, drop = FALSE])
rownames(mat325) <- rsem325$Gene_Name
mat325 <- collapse_expr(mat325)
common325 <- intersect(clinical325$sample_id, colnames(mat325))
clinical325 <- clinical325[match(common325, sample_id)]
mat325 <- mat325[, clinical325$sample_id, drop = FALSE]
log325 <- log2(mat325 + 1)

tcga_sub <- score_bulk_projection("TCGA", tcga_base, tcga_log2, "sample_key", submodule_sets)
cgga_sub <- rbindlist(list(
  score_bulk_projection("CGGA_mRNAseq_693", cgga_base[dataset == "CGGA_mRNAseq_693"], log693, "sample_key", submodule_sets),
  score_bulk_projection("CGGA_mRNAseq_325", cgga_base[dataset == "CGGA_mRNAseq_325"], log325, "sample_key", submodule_sets)
), fill = TRUE)

write_table(tcga_sub, "tcga_lap3_state_submodule_projection.csv")
write_table(cgga_sub, "cgga_lap3_state_submodule_projection.csv")

coverage <- rbindlist(list(
  coverage_for_projection(rownames(tcga_log2), "TCGA"),
  coverage_for_projection(rownames(log693), "CGGA_mRNAseq_693"),
  coverage_for_projection(rownames(log325), "CGGA_mRNAseq_325")
), fill = TRUE)
write_table(coverage, "lap3_state_submodule_bulk_gene_coverage.csv")

score_cols <- names(submodule_sets)
bulk_cor <- rbindlist(list(
  summarise_bulk_correlations(tcga_sub, score_cols, "TCGA"),
  summarise_bulk_correlations(cgga_sub[dataset == "CGGA_mRNAseq_693"], score_cols, "CGGA_mRNAseq_693"),
  summarise_bulk_correlations(cgga_sub[dataset == "CGGA_mRNAseq_325"], score_cols, "CGGA_mRNAseq_325")
), fill = TRUE)
bulk_cor[, p_adj_BH := p.adjust(p_value, method = "BH")]
write_table(bulk_cor, "lap3_state_submodule_bulk_correlations.csv")

cat("Projecting submodules into Core GBmap cache...\n")
gbmap <- run_gbmap_projection(submodule_sets)
gbmap_primary <- if (!is.null(gbmap)) {
  gbmap$association[
    entry_variant == "main_neoplastic_exclude_neftel2019" &
      threshold == 20L
  ][order(author_state, submodule)]
} else {
  data.table()
}
write_table(gbmap_primary, "gbmap_core_lap3_state_submodule_primary_summary.csv")

interpretation <- data.table(
  item = c(
    "n_submodules",
    "primary_next_step",
    "boundary",
    "figure_use"
  ),
  value = c(
    as.character(length(submodule_sets)),
    "Use submodules before CPTAC/GLASS projection and figure freeze",
    "Submodules explain a composite transcriptional state; they do not prove LAP3 enzymatic causality or leucine-mTORC1 activation",
    "Figure 2 internal structure; Figure 3 malignant-state decomposition; Figure 4 spatial topology sensitivity; Figure 5 conditional nomination"
  )
)
write_table(interpretation, "lap3_state_submodule_interpretation_summary.csv")

saveRDS(
  list(
    generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    submodule_sets = submodule_sets,
    gene_assignment = assignment,
    tcga_projection = tcga_sub,
    cgga_projection = cgga_sub,
    bulk_correlations = bulk_cor,
    gbmap_primary = gbmap_primary
  ),
  file.path(export_dir, "lap3_state_submodules_projection.rds")
)

readme <- c(
  "# LAP3 State Submodules",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Purpose",
  "",
  "Decompose the frozen 207-gene `LAP3_STATE_UNION` into interpretable submodules before CPTAC/GLASS projection and figure freezing.",
  "",
  "## Submodules",
  "",
  "- `LAP3_ANABOLIC_TRANSLATION_MODULE`",
  "- `LAP3_PROTEOSTASIS_STRESS_MODULE`",
  "- `LAP3_MYELOID_TAM_CONTEXT_MODULE`",
  "- `LAP3_HYPOXIA_PERINECROTIC_MODULE`",
  "- `LAP3_MALIGNANT_STATE_MODULE`",
  "",
  "## Assignment Rule",
  "",
  "Each union gene receives transparent flags for translation/anabolic, proteostasis/stress, myeloid/TAM, hypoxia/perinecrotic and malignant-state evidence. The primary submodule is assigned by deterministic priority: myeloid/TAM, hypoxia, proteostasis, anabolic/translation, then malignant-state/residual.",
  "",
  "## Key Outputs",
  "",
  "- `tables/lap3_state_submodule_gene_assignment.csv`",
  "- `tables/lap3_state_submodule_gene_counts.csv`",
  "- `tables/tcga_lap3_state_submodule_projection.csv`",
  "- `tables/cgga_lap3_state_submodule_projection.csv`",
  "- `tables/lap3_state_submodule_bulk_correlations.csv`",
  "- `tables/gbmap_core_lap3_state_submodule_lap3_associations.csv`",
  "- `source_data/gbmap_core_lap3_state_submodule_cell_scores.csv.gz`",
  "- `exports/lap3_state_submodules_projection.rds`",
  "",
  "## Gene Counts",
  "",
  paste(capture.output(print(submodule_counts)), collapse = "\n"),
  "",
  "## Interpretation Boundary",
  "",
  "These submodules explain the internal structure of a composite transcriptional state. They are not causal evidence for LAP3 enzymatic activity, intracellular leucine availability, or direct mTORC1 phosphorylation."
)
writeLines(readme, file.path(out_dir, "README.md"))

cat("Completed:", format(Sys.time()), "\n")
