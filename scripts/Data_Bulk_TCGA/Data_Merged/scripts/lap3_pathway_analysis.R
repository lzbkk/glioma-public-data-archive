#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(limma)
  library(edgeR)
  library(fgsea)
  library(BiocParallel)
  library(msigdbr)
  library(pheatmap)
  library(ggrepel)
})

setwd("/home/lzb/glioma/Data_Bulk_TCGA/Data_Merged")
set.seed(123)

out_dir <- file.path("results", "LAP3_Pathway")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

plot_dir <- file.path(out_dir, "plots")
table_dir <- file.path(out_dir, "tables")
cmap_dir <- file.path(out_dir, "cmap_inputs")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cmap_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "lap3_pathway_analysis.log")
sink(log_file, split = TRUE)
on.exit({
  while (sink.number() > 0) sink()
}, add = TRUE)

cat("LAP3 pathway analysis\n")
cat("Run time:", as.character(Sys.time()), "\n")
cat("Working directory:", getwd(), "\n\n")

save_plot <- function(plot, filename, width = 7, height = 5) {
  ggsave(file.path(plot_dir, paste0(filename, ".pdf")), plot, width = width, height = height)
  ggsave(file.path(plot_dir, paste0(filename, ".png")), plot, width = width, height = height, dpi = 300)
}

write_table <- function(x, filename) {
  write.csv(x, file.path(table_dir, filename), row.names = FALSE)
}

clean_factor <- function(x) {
  x <- as.character(x)
  x[x == "" | x == "NA" | is.na(x)] <- NA_character_
  x
}

module_score <- function(expr_log2, genes) {
  genes <- intersect(genes, rownames(expr_log2))
  if (length(genes) < 3) {
    return(rep(NA_real_, ncol(expr_log2)))
  }
  z <- t(scale(t(expr_log2[genes, , drop = FALSE])))
  colMeans(z, na.rm = TRUE)
}

make_pathway_list <- function(msig) {
  split(msig$gene_symbol, msig$gs_name)
}

run_de <- function(counts, meta, group_col, covariates = character(), label) {
  cat("\nDE analysis:", label, "\n")
  keep_meta <- !is.na(meta[[group_col]])
  if (length(covariates) > 0) {
    for (cc in covariates) keep_meta <- keep_meta & !is.na(meta[[cc]])
  }
  meta2 <- meta[keep_meta, , drop = FALSE]
  counts2 <- counts[, meta2$barcode, drop = FALSE]
  meta2$group <- relevel(factor(meta2[[group_col]]), ref = "Low")
  if (nlevels(meta2$group) < 2) stop("Group has <2 levels for ", label)

  usable_covariates <- covariates[vapply(covariates, function(cc) {
    x <- meta2[[cc]]
    if (is.numeric(x)) {
      return(sum(!is.na(x)) > 5 && stats::sd(x, na.rm = TRUE) > 0)
    }
    length(unique(x[!is.na(x)])) >= 2
  }, logical(1))]

  for (cc in usable_covariates) {
    if (!is.numeric(meta2[[cc]])) meta2[[cc]] <- factor(meta2[[cc]])
  }

  y <- DGEList(counts = counts2)
  keep <- filterByExpr(y, group = meta2$group)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)

  rhs <- c("group", usable_covariates)
  design <- model.matrix(as.formula(paste("~", paste(rhs, collapse = " + "))), data = meta2)
  non_estimable <- limma::nonEstimable(design)
  if (length(non_estimable) > 0) {
    if ("groupHigh" %in% non_estimable) {
      stop("groupHigh is not estimable for ", label)
    }
    cat("Dropping non-estimable design columns:", paste(non_estimable, collapse = ", "), "\n")
    design <- design[, !colnames(design) %in% non_estimable, drop = FALSE]
  }
  v <- voom(y, design, plot = FALSE)
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  coef_name <- grep("^groupHigh$", colnames(design), value = TRUE)
  if (length(coef_name) != 1) stop("Cannot find groupHigh coefficient for ", label)

  res <- topTable(fit, coef = coef_name, number = Inf, sort.by = "P")
  res$gene <- rownames(res)
  res <- res %>%
    relocate(gene) %>%
    mutate(
      comparison = label,
      n_samples = nrow(meta2),
      n_high = sum(meta2$group == "High"),
      n_low = sum(meta2$group == "Low"),
      rank_t = t,
      rank_signed_logp = sign(logFC) * -log10(P.Value)
    )

  cat("Samples:", nrow(meta2), "High:", sum(meta2$group == "High"), "Low:", sum(meta2$group == "Low"), "\n")
  cat("Genes after filterByExpr:", nrow(res), "\n")
  cat("Covariates:", ifelse(length(usable_covariates), paste(usable_covariates, collapse = ", "), "none"), "\n")
  res
}

run_fgsea <- function(de_table, pathways, label) {
  ranks <- de_table$rank_t
  names(ranks) <- de_table$gene
  ranks <- ranks[!is.na(ranks)]
  ranks <- sort(ranks, decreasing = TRUE)
  fg <- fgsea(
    pathways = pathways,
    stats = ranks,
    minSize = 10,
    maxSize = 500,
    eps = 0,
    BPPARAM = BiocParallel::SerialParam()
  )
  fg %>%
    arrange(padj, desc(abs(NES))) %>%
    mutate(comparison = label, leadingEdge = vapply(leadingEdge, paste, character(1), collapse = ";")) %>%
    relocate(comparison)
}

plot_volcano <- function(de_table, filename, title) {
  d <- de_table %>%
    mutate(
      neg_log10_p = -log10(P.Value),
      sig = case_when(
        adj.P.Val < 0.05 & logFC > 0 ~ "High in LAP3-high",
        adj.P.Val < 0.05 & logFC < 0 ~ "Low in LAP3-high",
        TRUE ~ "NS"
      )
    )
  label_genes <- d %>%
    filter(adj.P.Val < 0.05) %>%
    arrange(P.Value) %>%
    slice_head(n = 18)
  p <- ggplot(d, aes(logFC, neg_log10_p, color = sig)) +
    geom_point(alpha = 0.55, size = 0.9) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey70") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey70") +
    ggrepel::geom_text_repel(data = label_genes, aes(label = gene), size = 2.6, max.overlaps = 25) +
    scale_color_manual(values = c("High in LAP3-high" = "#C23B23", "Low in LAP3-high" = "#2878B5", "NS" = "grey75")) +
    labs(title = title, x = "log2 fold-change, LAP3-high vs low", y = "-log10(P value)", color = NULL) +
    theme_bw(base_size = 11) +
    theme(legend.position = "top")
  save_plot(p, filename, width = 7, height = 5.5)
}

plot_fgsea_dot <- function(fg_table, filename, title, patterns = NULL, n = 25) {
  d <- fg_table
  if (!is.null(patterns)) {
    d <- d %>% filter(grepl(patterns, pathway, ignore.case = TRUE))
  }
  d <- d %>% arrange(padj, desc(abs(NES))) %>% slice_head(n = n)
  if (!nrow(d)) return(invisible(NULL))
  d$pathway_label <- gsub("^HALLMARK_|^REACTOME_", "", d$pathway)
  d$pathway_label <- gsub("_", " ", d$pathway_label)
  d$pathway_label <- factor(d$pathway_label, levels = rev(d$pathway_label))
  p <- ggplot(d, aes(NES, pathway_label, size = -log10(padj), color = NES)) +
    geom_point() +
    scale_color_gradient2(low = "#2878B5", mid = "grey85", high = "#C23B23", midpoint = 0) +
    labs(title = title, x = "Normalized enrichment score", y = NULL, size = "-log10(FDR)", color = "NES") +
    theme_bw(base_size = 10)
  save_plot(p, filename, width = 7.5, height = max(4, 0.22 * nrow(d) + 1.8))
}

plot_enrichment_if_present <- function(de_table, pathways, pathway_name, filename, title) {
  if (!pathway_name %in% names(pathways)) return(invisible(NULL))
  ranks <- de_table$rank_t
  names(ranks) <- de_table$gene
  ranks <- sort(ranks[!is.na(ranks)], decreasing = TRUE)
  p <- fgsea::plotEnrichment(pathways[[pathway_name]], ranks) +
    labs(title = title, x = "Rank", y = "Enrichment score") +
    theme_bw(base_size = 11)
  save_plot(p, filename, width = 6.8, height = 4.6)
}

cat("Reading inputs...\n")
clinical <- readRDS(file.path("results", "Clinical_Field_QC", "clinical_glioma_analysis_fields.rds"))
expr_tpm <- readRDS(file.path("data_analysis", "expr_tpm_glioma_uni.rds"))
expr_count <- readRDS(file.path("data_analysis", "expr_count_glioma_uni.rds"))

stopifnot(identical(clinical$barcode, colnames(expr_tpm)[-1]))
stopifnot(identical(clinical$barcode, colnames(expr_count)[-1]))
stopifnot("LAP3" %in% rownames(expr_tpm))
stopifnot("LAP3" %in% rownames(expr_count))

gene_type <- expr_tpm$gene_type
protein_coding <- !is.na(gene_type) & gene_type == "protein_coding"
cat("Clinical samples:", nrow(clinical), "\n")
cat("Expression genes:", nrow(expr_tpm), "\n")
cat("Protein-coding genes:", sum(protein_coding), "\n")

count_mat <- as.matrix(expr_count[protein_coding, -1, drop = FALSE])
storage.mode(count_mat) <- "integer"
tpm_mat <- as.matrix(expr_tpm[protein_coding, -1, drop = FALSE])
storage.mode(tpm_mat) <- "double"
expr_log2 <- log2(tpm_mat + 1)

lap3_tpm <- as.numeric(expr_tpm["LAP3", -1, drop = TRUE])
names(lap3_tpm) <- colnames(expr_tpm)[-1]
clinical <- clinical %>%
  mutate(
    LAP3_tpm = lap3_tpm[barcode],
    LAP3_log2_tpm = log2(LAP3_tpm + 1),
    LAP3_group_all = ifelse(LAP3_log2_tpm >= median(LAP3_log2_tpm, na.rm = TRUE), "High", "Low")
  ) %>%
  group_by(cohort) %>%
  mutate(LAP3_group_by_cohort = ifelse(LAP3_log2_tpm >= median(LAP3_log2_tpm, na.rm = TRUE), "High", "Low")) %>%
  ungroup() %>%
  mutate(
    LAP3_tertile_all = case_when(
      LAP3_log2_tpm <= quantile(LAP3_log2_tpm, 1 / 3, na.rm = TRUE) ~ "Low",
      LAP3_log2_tpm >= quantile(LAP3_log2_tpm, 2 / 3, na.rm = TRUE) ~ "High",
      TRUE ~ NA_character_
    )
  ) %>%
  group_by(cohort) %>%
  mutate(
    LAP3_tertile_by_cohort = case_when(
      LAP3_log2_tpm <= quantile(LAP3_log2_tpm, 1 / 3, na.rm = TRUE) ~ "Low",
      LAP3_log2_tpm >= quantile(LAP3_log2_tpm, 2 / 3, na.rm = TRUE) ~ "High",
      TRUE ~ NA_character_
    )
  ) %>%
  ungroup() %>%
  mutate(
    cohort = factor(cohort),
    grade = factor(clean_factor(grade), levels = c("G2", "G3", "G4")),
    idh_status = factor(clean_factor(idh_status)),
    codel_1p19q = factor(clean_factor(codel_1p19q)),
    mgmt_status = factor(clean_factor(mgmt_status))
  )

write_table(
  clinical %>%
    select(
      barcode, patient, cohort, grade, idh_status, codel_1p19q, mgmt_status,
      LAP3_tpm, LAP3_log2_tpm, LAP3_group_all, LAP3_group_by_cohort,
      LAP3_tertile_all, LAP3_tertile_by_cohort
    ),
  "lap3_pathway_sample_groups.csv"
)

stratification_audit <- clinical %>%
  pivot_longer(
    cols = c(LAP3_group_all, LAP3_group_by_cohort, LAP3_tertile_all, LAP3_tertile_by_cohort),
    names_to = "grouping_rule",
    values_to = "LAP3_group"
  ) %>%
  filter(!is.na(LAP3_group)) %>%
  count(grouping_rule, cohort, grade, idh_status, LAP3_group, name = "n") %>%
  arrange(grouping_rule, cohort, grade, idh_status, LAP3_group)
write_table(stratification_audit, "lap3_grouping_stratification_audit.csv")

cat("\nLoading MSigDB gene sets...\n")
hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>% select(gs_name, gene_symbol)
reactome <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>% select(gs_name, gene_symbol)

pathways_h_reactome <- make_pathway_list(bind_rows(hallmark, reactome))
pathways_all <- pathways_h_reactome

focus_pathways <- c(
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_G2M_CHECKPOINT",
  "REACTOME_AMINO_ACIDS_REGULATE_MTORC1",
  "REACTOME_MTORC1_MEDIATED_SIGNALLING",
  "REACTOME_MTOR_SIGNALLING",
  "REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES",
  "REACTOME_AMINO_ACID_TRANSPORT_ACROSS_THE_PLASMA_MEMBRANE",
  "REACTOME_BRANCHED_CHAIN_AMINO_ACID_CATABOLISM",
  "REACTOME_EUKARYOTIC_TRANSLATION_INITIATION",
  "REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION",
  "REACTOME_TRANSLATION",
  "REACTOME_RRNA_PROCESSING",
  "REACTOME_RIBOSOMAL_SCANNING_AND_START_CODON_RECOGNITION"
)
focus_pathways <- intersect(focus_pathways, names(pathways_all))

cat("Hallmark + Reactome pathways:", length(pathways_h_reactome), "\n")
cat("Main GSEA pathways, Hallmark + Reactome:", length(pathways_all), "\n")
cat("Focus pathways found:", paste(focus_pathways, collapse = ", "), "\n")

manual_sets <- list(
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
  )
)
score_sets <- c(pathways_all[focus_pathways], manual_sets)

cat("\nRunning differential expression...\n")
de_results <- list()
de_results$all_median_unadjusted <- run_de(count_mat, clinical, "LAP3_group_all", character(), "all_median_unadjusted")
de_results$all_median_adjusted <- run_de(count_mat, clinical, "LAP3_group_all", c("grade", "idh_status", "codel_1p19q", "age_years"), "all_median_adjusted")
de_results$all_tertile_unadjusted <- run_de(count_mat, clinical, "LAP3_tertile_all", character(), "all_tertile_unadjusted")
de_results$all_tertile_adjusted <- run_de(count_mat, clinical, "LAP3_tertile_all", c("grade", "idh_status", "codel_1p19q", "age_years"), "all_tertile_adjusted")
de_results$cohort_balanced_tertile_adjusted <- run_de(count_mat, clinical, "LAP3_tertile_by_cohort", c("cohort", "grade", "idh_status", "codel_1p19q", "age_years"), "cohort_balanced_tertile_adjusted")
de_results$GBM_median_unadjusted <- run_de(count_mat, clinical %>% filter(cohort == "GBM"), "LAP3_group_by_cohort", character(), "GBM_median_unadjusted")
de_results$LGG_median_unadjusted <- run_de(count_mat, clinical %>% filter(cohort == "LGG"), "LAP3_group_by_cohort", character(), "LGG_median_unadjusted")
de_results$GBM_tertile_unadjusted <- run_de(count_mat, clinical %>% filter(cohort == "GBM"), "LAP3_tertile_by_cohort", character(), "GBM_tertile_unadjusted")
de_results$LGG_tertile_unadjusted <- run_de(count_mat, clinical %>% filter(cohort == "LGG"), "LAP3_tertile_by_cohort", character(), "LGG_tertile_unadjusted")

for (nm in names(de_results)) {
  write_table(de_results[[nm]], paste0("deg_", nm, ".csv"))
}
write_table(bind_rows(de_results), "deg_all_comparisons_combined.csv")

plot_volcano(de_results$all_median_unadjusted, "Fig3A_volcano_all_median_unadjusted", "All glioma: LAP3-high vs low")
plot_volcano(de_results$all_tertile_unadjusted, "Fig3A_volcano_all_tertile_unadjusted", "All glioma: LAP3 upper vs lower tertile")

cat("\nRunning fgsea...\n")
fgsea_results <- lapply(names(de_results), function(nm) {
  run_fgsea(de_results[[nm]], pathways_all, nm)
})
names(fgsea_results) <- names(de_results)

for (nm in names(fgsea_results)) {
  write_table(fgsea_results[[nm]], paste0("fgsea_", nm, ".csv"))
}
write_table(bind_rows(fgsea_results), "fgsea_all_comparisons_combined.csv")

focus_fgsea <- bind_rows(fgsea_results) %>% filter(pathway %in% focus_pathways)
write_table(focus_fgsea, "fgsea_focus_pathways.csv")

plot_fgsea_dot(
  fgsea_results$all_median_unadjusted,
  "Fig3B_fgsea_focus_all_median_unadjusted",
  "All glioma: focused GSEA pathways",
  patterns = "MTOR|AMINO|BRANCHED|TRANSLATION|RIBOSOM|MYC|E2F|G2M",
  n = 30
)
plot_fgsea_dot(
  fgsea_results$all_tertile_unadjusted,
  "Fig3B_fgsea_focus_all_tertile_unadjusted",
  "All glioma tertiles: focused GSEA pathways",
  patterns = "MTOR|AMINO|BRANCHED|TRANSLATION|RIBOSOM|MYC|E2F|G2M",
  n = 30
)

plot_enrichment_if_present(
  de_results$all_median_unadjusted, pathways_all, "HALLMARK_MTORC1_SIGNALING",
  "Fig3C_enrichment_HALLMARK_MTORC1_all_median",
  "HALLMARK_MTORC1_SIGNALING"
)
plot_enrichment_if_present(
  de_results$all_median_unadjusted, pathways_all, "REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES",
  "Fig3C_enrichment_REACTOME_amino_acid_metabolism_all_median",
  "REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES"
)
plot_enrichment_if_present(
  de_results$all_median_unadjusted, pathways_all, "REACTOME_AMINO_ACIDS_REGULATE_MTORC1",
  "Fig3C_enrichment_REACTOME_amino_acids_regulate_mtorc1_all_median",
  "REACTOME_AMINO_ACIDS_REGULATE_MTORC1"
)
plot_enrichment_if_present(
  de_results$all_median_unadjusted, pathways_all, "REACTOME_BRANCHED_CHAIN_AMINO_ACID_CATABOLISM",
  "Fig3C_enrichment_REACTOME_BCAA_catabolism_all_median",
  "REACTOME_BRANCHED_CHAIN_AMINO_ACID_CATABOLISM"
)

cat("\nComputing module scores...\n")
score_df <- clinical %>%
  select(
    barcode, patient, cohort, grade, idh_status, codel_1p19q, mgmt_status,
    LAP3_log2_tpm, LAP3_group_all, LAP3_group_by_cohort,
    LAP3_tertile_all, LAP3_tertile_by_cohort
  )

for (set_name in names(score_sets)) {
  score_df[[set_name]] <- module_score(expr_log2, score_sets[[set_name]])
}
write_table(score_df, "lap3_pathway_module_scores.csv")

score_long <- score_df %>%
  pivot_longer(cols = all_of(names(score_sets)), names_to = "pathway", values_to = "score")

group_tests <- score_long %>%
  filter(!is.na(score), !is.na(LAP3_group_all)) %>%
  group_by(pathway) %>%
  summarise(
    n = n(),
    median_high = median(score[LAP3_group_all == "High"], na.rm = TRUE),
    median_low = median(score[LAP3_group_all == "Low"], na.rm = TRUE),
    p_value = wilcox.test(score ~ LAP3_group_all)$p.value,
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)
write_table(group_tests, "pathway_score_high_low_tests.csv")

cohort_tertile_group_tests <- score_long %>%
  filter(!is.na(score), !is.na(LAP3_tertile_by_cohort)) %>%
  group_by(cohort, pathway) %>%
  summarise(
    n = n(),
    median_high = median(score[LAP3_tertile_by_cohort == "High"], na.rm = TRUE),
    median_low = median(score[LAP3_tertile_by_cohort == "Low"], na.rm = TRUE),
    p_value = wilcox.test(score ~ LAP3_tertile_by_cohort)$p.value,
    .groups = "drop"
  ) %>%
  group_by(cohort) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  arrange(cohort, p_adj)
write_table(cohort_tertile_group_tests, "pathway_score_high_low_tests_by_cohort_tertile.csv")

cor_tests <- score_long %>%
  filter(!is.na(score), !is.na(LAP3_log2_tpm)) %>%
  group_by(pathway) %>%
  summarise(
    n = n(),
    spearman_rho = suppressWarnings(cor(LAP3_log2_tpm, score, method = "spearman", use = "complete.obs")),
    p_value = suppressWarnings(cor.test(LAP3_log2_tpm, score, method = "spearman")$p.value),
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)
write_table(cor_tests, "pathway_score_lap3_correlations.csv")

cohort_cor_tests <- score_long %>%
  filter(!is.na(score), !is.na(LAP3_log2_tpm)) %>%
  group_by(cohort, pathway) %>%
  summarise(
    n = n(),
    spearman_rho = suppressWarnings(cor(LAP3_log2_tpm, score, method = "spearman", use = "complete.obs")),
    p_value = suppressWarnings(cor.test(LAP3_log2_tpm, score, method = "spearman")$p.value),
    .groups = "drop"
  ) %>%
  group_by(cohort) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  arrange(cohort, p_adj)
write_table(cohort_cor_tests, "pathway_score_lap3_correlations_by_cohort.csv")

score_plot_sets <- intersect(c(
  "HALLMARK_MTORC1_SIGNALING",
  "REACTOME_AMINO_ACIDS_REGULATE_MTORC1",
  "REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES",
  "REACTOME_BRANCHED_CHAIN_AMINO_ACID_CATABOLISM",
  "REACTOME_EUKARYOTIC_TRANSLATION_INITIATION",
  "LEUCINE_BCAA_CORE",
  "MTORC1_READOUT_CORE"
), names(score_sets))

plot_score_df <- score_long %>%
  filter(pathway %in% score_plot_sets) %>%
  mutate(
    pathway_label = gsub("^HALLMARK_|^REACTOME_", "", pathway),
    pathway_label = gsub("_", " ", pathway_label)
  )

p_box <- ggplot(plot_score_df, aes(LAP3_group_all, score, fill = LAP3_group_all)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.25, size = 0.45) +
  facet_wrap(~ pathway_label, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c("Low" = "#2878B5", "High" = "#C23B23")) +
  labs(title = "Pathway module scores by LAP3 group", x = NULL, y = "Module score", fill = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position = "top")
save_plot(p_box, "Fig3D_pathway_scores_by_LAP3_group", width = 8.5, height = 8)

scatter_df <- score_long %>%
  filter(pathway %in% score_plot_sets) %>%
  mutate(
    pathway_label = gsub("^HALLMARK_|^REACTOME_", "", pathway),
    pathway_label = gsub("_", " ", pathway_label)
  )
p_scatter <- ggplot(scatter_df, aes(LAP3_log2_tpm, score, color = cohort)) +
  geom_point(alpha = 0.55, size = 0.9) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.55, color = "black") +
  facet_wrap(~ pathway_label, scales = "free_y", ncol = 2) +
  scale_color_manual(values = c("GBM" = "#C23B23", "LGG" = "#2878B5")) +
  labs(title = "LAP3 expression vs pathway module scores", x = "LAP3 log2(TPM + 1)", y = "Module score") +
  theme_bw(base_size = 10) +
  theme(legend.position = "top")
save_plot(p_scatter, "Fig3E_LAP3_vs_pathway_scores", width = 8.5, height = 8)

key_genes <- unique(c(
  "LAP3", "MTOR", "RPTOR", "RHEB", "RPS6KB1", "EIF4EBP1", "RPS6",
  "SLC7A5", "SLC3A2", "SLC38A2", "SLC38A9", "BCAT1", "BCAT2",
  "BCKDHA", "BCKDHB", "DBT", "DLD", "BCKDK", "PPM1K",
  "MKI67", "TOP2A", "MCM2", "CDK1", "SOX2", "NES", "PROM1", "VIM", "CD44"
))
key_genes <- intersect(key_genes, rownames(expr_log2))
ann <- clinical %>%
  select(barcode, cohort, grade, idh_status, LAP3_tertile_by_cohort) %>%
  as.data.frame()
rownames(ann) <- ann$barcode
ann$barcode <- NULL
heat_mat <- expr_log2[key_genes, clinical$barcode, drop = FALSE]
heat_mat_z <- t(scale(t(heat_mat)))
pdf(file.path(plot_dir, "Fig3F_key_gene_heatmap.pdf"), width = 10, height = 7)
pheatmap(
  heat_mat_z,
  annotation_col = ann,
  show_colnames = FALSE,
  fontsize_row = 8,
  clustering_method = "ward.D2",
  main = "LAP3, mTORC1, BCAA/leucine, proliferation and stemness genes"
)
dev.off()
png(file.path(plot_dir, "Fig3F_key_gene_heatmap.png"), width = 3000, height = 2100, res = 300)
pheatmap(
  heat_mat_z,
  annotation_col = ann,
  show_colnames = FALSE,
  fontsize_row = 8,
  clustering_method = "ward.D2",
  main = "LAP3, mTORC1, BCAA/leucine, proliferation and stemness genes"
)
dev.off()

cat("\nPreparing CMap input gene lists...\n")
cmap_source <- de_results$cohort_balanced_tertile_adjusted %>%
  filter(!is.na(logFC), !is.na(adj.P.Val)) %>%
  arrange(desc(logFC))
cmap_up <- cmap_source %>%
  filter(logFC > 0) %>%
  arrange(adj.P.Val, desc(logFC)) %>%
  slice_head(n = 150) %>%
  pull(gene)
cmap_down <- cmap_source %>%
  filter(logFC < 0) %>%
  arrange(adj.P.Val, logFC) %>%
  slice_head(n = 150) %>%
  pull(gene)
writeLines(cmap_up, file.path(cmap_dir, "LAP3_high_vs_low_cohort_balanced_tertile_top150_up.txt"))
writeLines(cmap_down, file.path(cmap_dir, "LAP3_high_vs_low_cohort_balanced_tertile_top150_down.txt"))
write.csv(
  data.frame(direction = c(rep("up", length(cmap_up)), rep("down", length(cmap_down))), gene = c(cmap_up, cmap_down)),
  file.path(cmap_dir, "LAP3_high_vs_low_cohort_balanced_tertile_CMap_genes.csv"),
  row.names = FALSE
)

summary_tables <- list(
  sample_counts = clinical %>% count(cohort, LAP3_group_all, LAP3_tertile_by_cohort),
  focus_fgsea_all_median = focus_fgsea %>% filter(comparison == "all_median_unadjusted") %>% arrange(padj),
  focus_fgsea_cohort_balanced = focus_fgsea %>% filter(comparison == "cohort_balanced_tertile_adjusted") %>% arrange(padj),
  score_tests_top = group_tests %>% filter(pathway %in% score_plot_sets) %>% arrange(p_adj),
  score_tests_by_cohort_tertile_top = cohort_tertile_group_tests %>% filter(pathway %in% score_plot_sets) %>% arrange(cohort, p_adj),
  score_correlations_top = cor_tests %>% filter(pathway %in% score_plot_sets) %>% arrange(p_adj)
)
saveRDS(summary_tables, file.path(out_dir, "lap3_pathway_summary_tables.rds"))

cat("\nWriting README...\n")
fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "fg"))
}
fmt_sci <- function(x) {
  ifelse(is.na(x), "NA", formatC(x, digits = 3, format = "e"))
}
focus_readme <- focus_fgsea %>%
  filter(comparison == "all_median_unadjusted", pathway %in% c(
    "HALLMARK_E2F_TARGETS",
    "HALLMARK_G2M_CHECKPOINT",
    "HALLMARK_MTORC1_SIGNALING",
    "HALLMARK_MYC_TARGETS_V1"
  )) %>%
  arrange(padj) %>%
  mutate(line = paste0("| ", pathway, " | ", fmt_num(NES), " | ", fmt_sci(padj), " |"))
focus_readme_tertile <- focus_fgsea %>%
  filter(comparison == "all_tertile_unadjusted", pathway %in% c(
    "HALLMARK_E2F_TARGETS",
    "HALLMARK_G2M_CHECKPOINT",
    "HALLMARK_MTORC1_SIGNALING"
  )) %>%
  arrange(padj) %>%
  mutate(line = paste0("| ", pathway, " | ", fmt_num(NES), " | ", fmt_sci(padj), " |"))
focus_readme_balanced <- focus_fgsea %>%
  filter(comparison == "cohort_balanced_tertile_adjusted", pathway %in% c(
    "HALLMARK_E2F_TARGETS",
    "HALLMARK_G2M_CHECKPOINT",
    "HALLMARK_MTORC1_SIGNALING",
    "HALLMARK_MYC_TARGETS_V1"
  )) %>%
  arrange(padj) %>%
  mutate(line = paste0("| ", pathway, " | ", fmt_num(NES), " | ", fmt_sci(padj), " |"))
focus_readme_gbm <- focus_fgsea %>%
  filter(comparison == "GBM_median_unadjusted", pathway %in% c(
    "HALLMARK_MTORC1_SIGNALING",
    "HALLMARK_MYC_TARGETS_V1",
    "REACTOME_TRANSLATION",
    "REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES"
  )) %>%
  arrange(padj) %>%
  mutate(line = paste0("| ", pathway, " | ", fmt_num(NES), " | ", fmt_sci(padj), " |"))
score_readme <- cor_tests %>%
  filter(pathway %in% c(
    "HALLMARK_MTORC1_SIGNALING",
    "LEUCINE_BCAA_CORE",
    "REACTOME_MTORC1_MEDIATED_SIGNALLING",
    "REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES",
    "MTORC1_READOUT_CORE",
    "REACTOME_TRANSLATION",
    "REACTOME_BRANCHED_CHAIN_AMINO_ACID_CATABOLISM"
  )) %>%
  arrange(p_adj) %>%
  mutate(line = paste0("| ", pathway, " | ", fmt_num(spearman_rho), " |"))
cohort_readme <- cohort_cor_tests %>%
  filter(pathway %in% c("HALLMARK_MTORC1_SIGNALING", "LEUCINE_BCAA_CORE")) %>%
  arrange(cohort, pathway) %>%
  mutate(line = paste0("- ", cohort, "：LAP3 与 ", pathway, " rho = ", fmt_num(spearman_rho), "。"))
cohort_tertile_readme <- cohort_tertile_group_tests %>%
  filter(pathway %in% c("HALLMARK_MTORC1_SIGNALING", "LEUCINE_BCAA_CORE")) %>%
  arrange(cohort, pathway) %>%
  mutate(line = paste0("- ", cohort, "：", pathway, " high-low median difference = ", fmt_num(median_high - median_low), "。"))

readme <- c(
  "# LAP3 Pathway Analysis",
  "",
  paste0("生成时间：", as.character(Sys.time())),
  "",
  "## 输入",
  "",
  "- `results/Clinical_Field_QC/clinical_glioma_analysis_fields.rds`",
  "- `data_analysis/expr_count_glioma_uni.rds`",
  "- `data_analysis/expr_tpm_glioma_uni.rds`",
  "",
  "## 方法概述",
  "",
  "- 使用 protein-coding gene 的 count 矩阵进行 `limma-voom` 差异分析。",
  "- 严格框架把 GBM 与 LGG 分开处理：主判断优先使用 cohort 内 median/tertile；all glioma 比较只作为背景敏感性。",
  "- 主要比较包括 all glioma median/tertile、cohort-balanced tertile adjusted、GBM median/tertile、LGG median/tertile。",
  "- 使用 `fgsea` 和 MSigDB Hallmark/Reactome 做第一版主 GSEA。",
  "- 使用 log2(TPM + 1) 计算透明的 z-score module score。",
  "- 导出 cohort 内 upper/lower tertile 并经 cohort、grade、IDH、1p/19q 和 age 调整后的 top 150 up/down genes 作为 CMap/LINCS 输入。",
  "",
  "## 关键输出",
  "",
  "- `tables/deg_*.csv`：差异分析结果。",
  "- `tables/fgsea_*.csv`：GSEA 结果。",
  "- `tables/fgsea_focus_pathways.csv`：mTORC1、氨基酸代谢、BCAA、translation、cell-cycle 等重点通路。",
  "- `tables/lap3_pathway_module_scores.csv`：样本级 pathway/module score。",
  "- `tables/pathway_score_lap3_correlations.csv`：LAP3 与 pathway score 的相关性。",
  "- `tables/lap3_grouping_stratification_audit.csv`：不同 LAP3 分组规则下 GBM/LGG、grade 和 IDH 组成审计。",
  "- `tables/pathway_score_high_low_tests_by_cohort_tertile.csv`：GBM/LGG 内部 upper/lower tertile 的 module score 检验。",
  "- `cmap_inputs/`：CMap/LINCS 输入基因列表。",
  "- `plots/`：Figure 3 第一版图表。",
  "",
  "## 第一版结果摘要",
  "",
  "### GSEA",
  "",
  "All glioma median split 中，LAP3-high 样本显著富集：",
  "",
  "| Pathway | NES | FDR |",
  "|---|---:|---:|",
  focus_readme$line,
  "",
  "All glioma upper/lower tertile 比较中趋势一致：",
  "",
  "| Pathway | NES | FDR |",
  "|---|---:|---:|",
  focus_readme_tertile$line,
  "",
  "严格框架下的 cohort-balanced upper/lower tertile adjusted 比较：",
  "",
  "| Pathway | NES | FDR |",
  "|---|---:|---:|",
  focus_readme_balanced$line,
  "",
  "GBM 内部 median split 中，LAP3-high 显示：",
  "",
  "| Pathway | NES | FDR |",
  "|---|---:|---:|",
  focus_readme_gbm$line,
  "",
  "LGG 内部 median split 中，LAP3-high 主要富集 E2F/G2M/mTORC1；translation 相关 Reactome 项在 ranked GSEA 中呈负向富集。这提示 LGG 内部的 LAP3 相关生物学可能更偏 cell-cycle/mTORC1，而 translation Reactome 结果需要谨慎解释。",
  "",
  "### Module score",
  "",
  "样本级 z-score module score 显示，LAP3 表达与多种机制相关 score 呈显著正相关：",
  "",
  "| Score | Spearman rho, all glioma |",
  "|---|---:|",
  score_readme$line,
  "",
  "分队列内相关性仍然明显：",
  "",
  cohort_readme$line,
  "",
  "分队列 upper/lower tertile 的 module score 差异用于替代 pan-glioma high/low 作为主分组判断：",
  "",
  cohort_tertile_readme$line,
  "",
  "### 阶段解释",
  "",
  "第一版结果支持以下克制结论：",
  "",
  "```text",
  "LAP3-high gliomas are associated with an mTORC1-active, cell-cycle/proliferative,",
  "and amino-acid/BCAA metabolism-related transcriptional state.",
  "```",
  "",
  "但 Reactome amino-acid/translation GSEA 在不同分组中方向并不完全一致，因此更适合作为 module score、cohort-balanced signature 和分队列结果共同支持的机制线索，而不是单独依赖某一个 pan-glioma GSEA 条目。",
  "",
  "## 解释边界",
  "",
  "本分析用于证明 LAP3-high glioma 与 mTORC1/氨基酸代谢/translation 等程序相关，属于计算机制线索和 hypothesis-generating evidence。不能单独证明 LAP3 直接调控亮氨酸或 mTORC1。",
  "主文结论应优先引用 GBM/LGG 内部分组、cohort-balanced adjusted signature 和协变量敏感性结果；all glioma global median/tertile 只保留为背景描述。"
)
writeLines(readme, file.path(out_dir, "README_LAP3_Pathway.md"))

cat("\nDone. Outputs written to:", normalizePath(out_dir), "\n")
