#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(msigdbr)
})

setwd("/home/lzb/glioma/Data_Bulk_TCGA/Data_Merged")
set.seed(123)

depmap_release <- "22Q2"
raw_dir <- file.path("data_raw", "depmap_22Q2")
out_dir <- file.path("results", "LAP3_DepMap")
plot_dir <- file.path(out_dir, "plots")
table_dir <- file.path(out_dir, "tables")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "lap3_depmap_analysis.log")
sink(log_file, split = TRUE)
on.exit({
  while (sink.number() > 0) sink()
}, add = TRUE)

cat("LAP3 DepMap initial screen\n")
cat("Run time:", as.character(Sys.time()), "\n")
cat("Working directory:", getwd(), "\n")
cat("DepMap release:", depmap_release, "\n\n")

metadata_file <- file.path(raw_dir, "sample_info.csv")
expression_file <- file.path(raw_dir, "CCLE_expression.csv")
crispr_file <- file.path(raw_dir, "CRISPR_gene_effect.csv")
stopifnot(file.exists(metadata_file), file.exists(expression_file), file.exists(crispr_file))

save_plot <- function(plot, filename, width = 7, height = 5) {
  ggsave(file.path(plot_dir, paste0(filename, ".pdf")), plot, width = width, height = height)
  ggsave(file.path(plot_dir, paste0(filename, ".png")), plot, width = width, height = height, dpi = 300)
}

write_table <- function(x, filename) {
  write.csv(x, file.path(table_dir, filename), row.names = FALSE)
}

fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "fg"))
}

fmt_sci <- function(x) {
  ifelse(is.na(x), "NA", formatC(x, digits = 3, format = "e"))
}

gene_symbol_from_col <- function(x) {
  sub(" \\([0-9]+\\)$", "", x)
}

find_gene_col <- function(headers, gene) {
  hits <- grep(paste0("^", gene, "( \\(|$)"), headers, value = TRUE)
  if (length(hits) != 1) {
    stop("Expected exactly one column for gene ", gene, ", found: ", paste(hits, collapse = ", "))
  }
  hits
}

z <- function(x) {
  x <- as.numeric(x)
  if (all(is.na(x)) || stats::sd(x, na.rm = TRUE) == 0) return(rep(NA_real_, length(x)))
  as.numeric(scale(x))
}

module_score <- function(expr_df, genes) {
  genes <- unique(genes)
  cols <- intersect(genes, names(expr_df))
  if (length(cols) < 3) {
    return(rep(NA_real_, nrow(expr_df)))
  }
  mat <- as.matrix(expr_df[, cols, drop = FALSE])
  storage.mode(mat) <- "double"
  zmat <- scale(mat)
  rowMeans(zmat, na.rm = TRUE)
}

cat("Reading metadata...\n")
metadata <- fread(metadata_file) %>%
  as.data.frame() %>%
  mutate(
    cell_line_name = na_if(cell_line_name, ""),
    stripped_cell_line_name = na_if(stripped_cell_line_name, ""),
    CCLE_Name = na_if(CCLE_Name, ""),
    lineage = na_if(lineage, ""),
    lineage_subtype = na_if(lineage_subtype, ""),
    lineage_sub_subtype = na_if(lineage_sub_subtype, ""),
    primary_disease = na_if(primary_disease, ""),
    Subtype = na_if(Subtype, ""),
    is_cns = lineage == "central_nervous_system",
    is_glioma = lineage == "central_nervous_system" & lineage_subtype == "glioma",
    is_gbm = is_glioma & lineage_sub_subtype == "glioblastoma"
  )

cat("Metadata rows:", nrow(metadata), "\n")
cat("CNS rows:", sum(metadata$is_cns, na.rm = TRUE), "\n")
cat("Glioma rows:", sum(metadata$is_glioma, na.rm = TRUE), "\n")
cat("GBM rows:", sum(metadata$is_gbm, na.rm = TRUE), "\n\n")

cat("Preparing gene sets...\n")
hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>%
  select(gs_name, gene_symbol) %>%
  distinct()

hallmark_mtorc1 <- hallmark %>%
  filter(gs_name == "HALLMARK_MTORC1_SIGNALING") %>%
  pull(gene_symbol) %>%
  unique()

manual_sets <- list(
  BCAA_LEUCINE_TRANSPORT_METABOLISM = c(
    "BCAT1", "BCAT2", "BCKDHA", "BCKDHB", "DBT", "DLD", "BCKDK", "PPM1K",
    "SLC7A5", "SLC3A2", "SLC43A1", "SLC43A2", "SLC38A2", "SLC38A9"
  ),
  MTORC1_READOUT_CORE = c(
    "MTOR", "RPTOR", "RHEB", "AKT1", "TSC1", "TSC2", "RRAGA", "RRAGB", "RRAGC",
    "RRAGD", "LAMTOR1", "LAMTOR2", "LAMTOR3", "LAMTOR4", "LAMTOR5",
    "RPS6KB1", "RPS6KB2", "EIF4EBP1", "EIF4EBP2", "RPS6", "EIF4E"
  )
)

score_sets <- c(
  list(HALLMARK_MTORC1_SIGNALING = hallmark_mtorc1),
  manual_sets
)
genes_needed <- unique(c("LAP3", unlist(score_sets, use.names = FALSE)))

cat("Inspecting expression header...\n")
expr_header <- names(fread(expression_file, nrows = 0))
expr_id_col <- expr_header[1]
expr_lap3_col <- find_gene_col(expr_header, "LAP3")
expr_gene_cols <- expr_header[gene_symbol_from_col(expr_header) %in% genes_needed]
expr_select <- unique(c(expr_id_col, expr_gene_cols))
cat("Expression ID column:", expr_id_col, "\n")
cat("Expression LAP3 column:", expr_lap3_col, "\n")
cat("Expression selected gene columns:", length(expr_gene_cols), "\n\n")

cat("Reading selected expression columns...\n")
expr <- fread(expression_file, select = expr_select)
setnames(expr, expr_id_col, "DepMap_ID")
expr_names <- names(expr)
expr_names[-1] <- gene_symbol_from_col(expr_names[-1])
setnames(expr, names(expr), make.unique(expr_names))
expr <- as.data.frame(expr)

cat("Inspecting CRISPR header...\n")
crispr_header <- names(fread(crispr_file, nrows = 0))
crispr_lap3_col <- find_gene_col(crispr_header, "LAP3")
cat("CRISPR LAP3 column:", crispr_lap3_col, "\n\n")

cat("Reading LAP3 CRISPR gene effect...\n")
crispr_lap3 <- fread(crispr_file, select = c("DepMap_ID", crispr_lap3_col))
setnames(crispr_lap3, crispr_lap3_col, "LAP3_gene_effect")
crispr_lap3 <- as.data.frame(crispr_lap3)

cat("Computing module scores...\n")
expr_scores <- expr %>%
  transmute(
    DepMap_ID,
    LAP3_expression_log2_tpm1 = LAP3,
    HALLMARK_MTORC1_SIGNALING_score = module_score(expr, score_sets$HALLMARK_MTORC1_SIGNALING),
    BCAA_LEUCINE_TRANSPORT_METABOLISM_score = module_score(expr, score_sets$BCAA_LEUCINE_TRANSPORT_METABOLISM),
    MTORC1_READOUT_CORE_score = module_score(expr, score_sets$MTORC1_READOUT_CORE)
  )

available_score_genes <- lapply(score_sets, function(gs) intersect(gs, names(expr)))
score_gene_coverage <- data.frame(
  score = names(available_score_genes),
  n_genes_in_set = vapply(score_sets, length, integer(1)),
  n_genes_available = vapply(available_score_genes, length, integer(1)),
  genes_available = vapply(available_score_genes, paste, character(1), collapse = ";")
)
write_table(score_gene_coverage, "depmap_score_gene_coverage.csv")
print(score_gene_coverage[, c("score", "n_genes_in_set", "n_genes_available")])

depmap_df <- metadata %>%
  left_join(expr_scores, by = "DepMap_ID") %>%
  left_join(crispr_lap3, by = "DepMap_ID") %>%
  mutate(
    display_cell_line = dplyr::coalesce(cell_line_name, stripped_cell_line_name, CCLE_Name, DepMap_ID),
    LAP3_dependency = -LAP3_gene_effect,
    has_expression = !is.na(LAP3_expression_log2_tpm1),
    has_crispr = !is.na(LAP3_gene_effect),
    has_both = has_expression & has_crispr,
    lineage_group = case_when(
      is_gbm ~ "CNS glioma: GBM",
      is_glioma ~ "CNS glioma: non-GBM",
      is_cns ~ "CNS non-glioma",
      TRUE ~ "Other"
    )
  )

write_table(depmap_df, "depmap_lap3_cell_line_dataset.csv")
write_table(depmap_df %>% filter(is_cns), "depmap_lap3_cns_cell_lines.csv")
write_table(depmap_df %>% filter(is_glioma), "depmap_lap3_glioma_cell_lines.csv")

coverage <- depmap_df %>%
  summarise(
    metadata_n = n(),
    expression_n = sum(has_expression),
    crispr_n = sum(has_crispr),
    both_n = sum(has_both),
    cns_metadata_n = sum(is_cns, na.rm = TRUE),
    cns_expression_n = sum(is_cns & has_expression, na.rm = TRUE),
    cns_crispr_n = sum(is_cns & has_crispr, na.rm = TRUE),
    cns_both_n = sum(is_cns & has_both, na.rm = TRUE),
    glioma_metadata_n = sum(is_glioma, na.rm = TRUE),
    glioma_expression_n = sum(is_glioma & has_expression, na.rm = TRUE),
    glioma_crispr_n = sum(is_glioma & has_crispr, na.rm = TRUE),
    glioma_both_n = sum(is_glioma & has_both, na.rm = TRUE),
    gbm_metadata_n = sum(is_gbm, na.rm = TRUE),
    gbm_expression_n = sum(is_gbm & has_expression, na.rm = TRUE),
    gbm_crispr_n = sum(is_gbm & has_crispr, na.rm = TRUE),
    gbm_both_n = sum(is_gbm & has_both, na.rm = TRUE)
  )
write_table(coverage, "depmap_lap3_coverage_summary.csv")
print(coverage)

lineage_summary <- depmap_df %>%
  filter(has_both) %>%
  group_by(lineage) %>%
  summarise(
    n = n(),
    median_LAP3_expression = median(LAP3_expression_log2_tpm1, na.rm = TRUE),
    median_LAP3_gene_effect = median(LAP3_gene_effect, na.rm = TRUE),
    median_LAP3_dependency = median(LAP3_dependency, na.rm = TRUE),
    median_HALLMARK_MTORC1 = median(HALLMARK_MTORC1_SIGNALING_score, na.rm = TRUE),
    median_BCAA = median(BCAA_LEUCINE_TRANSPORT_METABOLISM_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n))
write_table(lineage_summary, "depmap_lap3_lineage_summary.csv")

group_summary <- depmap_df %>%
  filter(has_both) %>%
  group_by(lineage_group) %>%
  summarise(
    n = n(),
    median_LAP3_expression = median(LAP3_expression_log2_tpm1, na.rm = TRUE),
    median_LAP3_gene_effect = median(LAP3_gene_effect, na.rm = TRUE),
    median_LAP3_dependency = median(LAP3_dependency, na.rm = TRUE),
    median_HALLMARK_MTORC1 = median(HALLMARK_MTORC1_SIGNALING_score, na.rm = TRUE),
    median_BCAA = median(BCAA_LEUCINE_TRANSPORT_METABOLISM_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(lineage_group)
write_table(group_summary, "depmap_lap3_group_summary.csv")

score_cols <- c(
  "LAP3_expression_log2_tpm1",
  "HALLMARK_MTORC1_SIGNALING_score",
  "BCAA_LEUCINE_TRANSPORT_METABOLISM_score",
  "MTORC1_READOUT_CORE_score"
)

cor_one <- function(df, subset_label, x_col, y_col = "LAP3_dependency") {
  d <- df %>% filter(!is.na(.data[[x_col]]), !is.na(.data[[y_col]]))
  if (nrow(d) < 5) {
    return(data.frame(subset = subset_label, x = x_col, y = y_col, n = nrow(d), spearman_rho = NA_real_, p_value = NA_real_))
  }
  ct <- suppressWarnings(cor.test(d[[x_col]], d[[y_col]], method = "spearman"))
  data.frame(
    subset = subset_label,
    x = x_col,
    y = y_col,
    n = nrow(d),
    spearman_rho = unname(ct$estimate),
    p_value = ct$p.value
  )
}

cor_results <- bind_rows(
  lapply(score_cols, function(x) cor_one(depmap_df %>% filter(has_both), "pan_cancer", x)),
  lapply(score_cols, function(x) cor_one(depmap_df %>% filter(is_cns, has_both), "CNS", x)),
  lapply(score_cols, function(x) cor_one(depmap_df %>% filter(is_glioma, has_both), "CNS_glioma", x)),
  lapply(score_cols, function(x) cor_one(depmap_df %>% filter(is_gbm, has_both), "CNS_glioma_GBM", x))
) %>%
  group_by(subset) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  arrange(subset, p_adj)
write_table(cor_results, "depmap_lap3_dependency_correlations.csv")

candidate <- depmap_df %>%
  filter(is_glioma, has_both) %>%
  mutate(
    LAP3_expression_z = z(LAP3_expression_log2_tpm1),
    LAP3_dependency_z = z(LAP3_dependency),
    HALLMARK_MTORC1_z = z(HALLMARK_MTORC1_SIGNALING_score),
    BCAA_z = z(BCAA_LEUCINE_TRANSPORT_METABOLISM_score),
    candidate_score = rowMeans(
      cbind(LAP3_expression_z, LAP3_dependency_z, HALLMARK_MTORC1_z, BCAA_z),
      na.rm = TRUE
    )
  ) %>%
  arrange(desc(candidate_score)) %>%
  select(
    DepMap_ID, display_cell_line, cell_line_name, stripped_cell_line_name, CCLE_Name,
    primary_disease, Subtype, lineage, lineage_subtype, lineage_sub_subtype,
    default_growth_pattern, source,
    LAP3_expression_log2_tpm1, LAP3_gene_effect, LAP3_dependency,
    HALLMARK_MTORC1_SIGNALING_score, BCAA_LEUCINE_TRANSPORT_METABOLISM_score,
    MTORC1_READOUT_CORE_score, candidate_score
  )
write_table(candidate, "depmap_candidate_glioma_cell_lines.csv")

common_glioma_cell_lines <- depmap_df %>%
  filter(
    is_glioma,
    grepl(
      "U87|U-87|U251|U-251|T98G|LN18|LN-18|A172|A-172|U118|U-118",
      paste(cell_line_name, stripped_cell_line_name, CCLE_Name),
      ignore.case = TRUE
    )
  ) %>%
  arrange(display_cell_line, DepMap_ID) %>%
  select(
    DepMap_ID, display_cell_line, cell_line_name, stripped_cell_line_name, CCLE_Name,
    lineage_sub_subtype, has_expression, has_crispr,
    LAP3_expression_log2_tpm1, LAP3_gene_effect, LAP3_dependency,
    HALLMARK_MTORC1_SIGNALING_score, BCAA_LEUCINE_TRANSPORT_METABOLISM_score,
    MTORC1_READOUT_CORE_score
  )
write_table(common_glioma_cell_lines, "depmap_common_glioma_cell_lines.csv")

cat("\nTop candidate glioma cell lines:\n")
print(head(candidate, 15))

top_lineages <- lineage_summary %>%
  filter(!is.na(lineage), n >= 20) %>%
  pull(lineage)
plot_lineages <- unique(c(top_lineages, "central_nervous_system"))

plot_df <- depmap_df %>%
  filter(has_both, lineage %in% plot_lineages) %>%
  mutate(lineage = factor(lineage, levels = lineage_summary$lineage[lineage_summary$lineage %in% plot_lineages]))

p_expr <- ggplot(plot_df, aes(lineage, LAP3_expression_log2_tpm1, fill = lineage == "central_nervous_system")) +
  geom_boxplot(outlier.shape = NA, width = 0.65) +
  geom_jitter(width = 0.18, alpha = 0.35, size = 0.6) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#C23B23", "FALSE" = "grey75")) +
  labs(title = "DepMap 22Q2: LAP3 expression by lineage", x = NULL, y = "LAP3 log2(TPM + 1)") +
  theme_bw(base_size = 10) +
  theme(legend.position = "none")
save_plot(p_expr, "Fig_DepMap_A_LAP3_expression_by_lineage", width = 7.5, height = 6.5)

p_effect <- ggplot(plot_df, aes(lineage, LAP3_gene_effect, fill = lineage == "central_nervous_system")) +
  geom_hline(yintercept = 0, color = "grey70") +
  geom_boxplot(outlier.shape = NA, width = 0.65) +
  geom_jitter(width = 0.18, alpha = 0.35, size = 0.6) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#C23B23", "FALSE" = "grey75")) +
  labs(title = "DepMap 22Q2: LAP3 CRISPR gene effect by lineage", x = NULL, y = "LAP3 gene effect; more negative = more dependent") +
  theme_bw(base_size = 10) +
  theme(legend.position = "none")
save_plot(p_effect, "Fig_DepMap_B_LAP3_gene_effect_by_lineage", width = 7.5, height = 6.5)

scatter_dep <- depmap_df %>%
  filter(has_both) %>%
  mutate(plot_group = factor(lineage_group, levels = c("CNS glioma: GBM", "CNS glioma: non-GBM", "CNS non-glioma", "Other")))

p_mtor <- ggplot(scatter_dep, aes(HALLMARK_MTORC1_SIGNALING_score, LAP3_dependency, color = plot_group)) +
  geom_point(alpha = 0.55, size = 1.1) +
  geom_smooth(data = scatter_dep %>% filter(is_glioma), method = "lm", se = FALSE, color = "black", linewidth = 0.55) +
  scale_color_manual(values = c(
    "CNS glioma: GBM" = "#C23B23",
    "CNS glioma: non-GBM" = "#E28A2B",
    "CNS non-glioma" = "#2878B5",
    "Other" = "grey75"
  )) +
  labs(title = "LAP3 dependency vs mTORC1 score", x = "HALLMARK_MTORC1_SIGNALING module score", y = "LAP3 dependency (-gene effect)", color = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position = "top")
save_plot(p_mtor, "Fig_DepMap_C_LAP3_dependency_vs_MTORC1_score", width = 7, height = 5.4)

p_bcaa <- ggplot(scatter_dep, aes(BCAA_LEUCINE_TRANSPORT_METABOLISM_score, LAP3_dependency, color = plot_group)) +
  geom_point(alpha = 0.55, size = 1.1) +
  geom_smooth(data = scatter_dep %>% filter(is_glioma), method = "lm", se = FALSE, color = "black", linewidth = 0.55) +
  scale_color_manual(values = c(
    "CNS glioma: GBM" = "#C23B23",
    "CNS glioma: non-GBM" = "#E28A2B",
    "CNS non-glioma" = "#2878B5",
    "Other" = "grey75"
  )) +
  labs(title = "LAP3 dependency vs BCAA/leucine score", x = "BCAA/leucine transport-metabolism module score", y = "LAP3 dependency (-gene effect)", color = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position = "top")
save_plot(p_bcaa, "Fig_DepMap_D_LAP3_dependency_vs_BCAA_score", width = 7, height = 5.4)

candidate_plot <- candidate %>%
  slice_head(n = 25) %>%
  mutate(display_cell_line = factor(display_cell_line, levels = rev(display_cell_line)))

p_candidate <- ggplot(candidate_plot, aes(candidate_score, display_cell_line, color = lineage_sub_subtype)) +
  geom_point(size = 2.2) +
  labs(title = "Candidate glioma cell lines for LAP3 axis follow-up", x = "Composite score: LAP3 expression + dependency + mTORC1 + BCAA", y = NULL, color = "Subtype") +
  theme_bw(base_size = 10) +
  theme(legend.position = "top")
save_plot(p_candidate, "Fig_DepMap_E_candidate_glioma_cell_lines", width = 7, height = 6)

cat("\nWriting README...\n")
cov <- coverage[1, ]
top_candidates <- candidate %>%
  slice_head(n = 10) %>%
  mutate(line = paste0(
    "| ", display_cell_line, " | ", DepMap_ID, " | ", lineage_sub_subtype, " | ",
    fmt_num(LAP3_expression_log2_tpm1), " | ", fmt_num(LAP3_gene_effect), " | ",
    fmt_num(HALLMARK_MTORC1_SIGNALING_score), " | ", fmt_num(BCAA_LEUCINE_TRANSPORT_METABOLISM_score), " |"
  ))
common_readme <- common_glioma_cell_lines %>%
  mutate(line = paste0(
    "| ", display_cell_line, " | ", DepMap_ID, " | ", lineage_sub_subtype, " | ",
    fmt_num(LAP3_expression_log2_tpm1), " | ", fmt_num(LAP3_gene_effect), " | ",
    fmt_num(HALLMARK_MTORC1_SIGNALING_score), " | ", fmt_num(BCAA_LEUCINE_TRANSPORT_METABOLISM_score), " |"
  ))
key_cors <- cor_results %>%
  filter(x %in% c("HALLMARK_MTORC1_SIGNALING_score", "BCAA_LEUCINE_TRANSPORT_METABOLISM_score", "LAP3_expression_log2_tpm1")) %>%
  mutate(line = paste0(
    "| ", subset, " | ", x, " | ", n, " | ", fmt_num(spearman_rho), " | ", fmt_sci(p_value), " |"
  ))

readme <- c(
  "# LAP3 DepMap Initial Screen",
  "",
  paste0("生成时间：", as.character(Sys.time())),
  "",
  "## 输入数据",
  "",
  paste0("- DepMap release：`", depmap_release, "`。"),
  "- `data_raw/depmap_22Q2/sample_info.csv`",
  "- `data_raw/depmap_22Q2/CCLE_expression.csv`",
  "- `data_raw/depmap_22Q2/CRISPR_gene_effect.csv`",
  "",
  "说明：本分析使用同一 DepMap release 的 metadata、expression 和 CRISPR gene effect，适合作为第一轮细胞系与依赖性初筛。它不是最新 release 的终版结论，后续如进入投稿图，可再用最新版 DepMap 复核。",
  "",
  "## 方法概述",
  "",
  "- 提取 `LAP3 (51056)` 的 log2(TPM + 1) expression。",
  "- 提取 `LAP3 (51056)` 的 CRISPR gene effect；gene effect 越负，表示敲除后细胞增殖/适应度下降越明显。",
  "- 为了便于解释，额外定义 `LAP3_dependency = -LAP3_gene_effect`，数值越高表示越依赖 LAP3。",
  "- CNS/glioma 细胞系根据 `lineage == central_nervous_system` 且 `lineage_subtype == glioma` 定义。",
  "- mTORC1 score 使用 MSigDB `HALLMARK_MTORC1_SIGNALING` 的 DepMap 表达 z-score 均值。",
  "- BCAA/leucine score 使用 BCAA 代谢与亮氨酸转运核心基因，且不包含 LAP3 本身，避免循环相关。",
  "",
  "## 覆盖情况",
  "",
  paste0("- metadata 总细胞系：", cov$metadata_n, "。"),
  paste0("- 同时有 LAP3 expression 和 CRISPR gene effect 的细胞系：", cov$both_n, "。"),
  paste0("- CNS metadata：", cov$cns_metadata_n, "；CNS 同时有 expression + CRISPR：", cov$cns_both_n, "。"),
  paste0("- CNS/glioma metadata：", cov$glioma_metadata_n, "；CNS/glioma 同时有 expression + CRISPR：", cov$glioma_both_n, "。"),
  paste0("- CNS/glioma/GBM metadata：", cov$gbm_metadata_n, "；GBM 同时有 expression + CRISPR：", cov$gbm_both_n, "。"),
  "",
  "## LAP3 dependency 与表达/通路 score 的相关性",
  "",
  "| Subset | Variable | n | Spearman rho with LAP3 dependency | P value |",
  "|---|---|---:|---:|---:|",
  key_cors$line,
  "",
  "解释边界：CNS/glioma 子集样本量远小于 pan-cancer，因此相关性应主要作为细胞系选择和假说生成依据，而不是机制证明。",
  "",
  "## 候选 glioma 细胞系",
  "",
  "候选排序综合考虑：LAP3 expression、LAP3 dependency、mTORC1 score、BCAA/leucine score。",
  "",
  "| Cell line | DepMap ID | Subtype | LAP3 expr | LAP3 gene effect | mTORC1 score | BCAA score |",
  "|---|---|---|---:|---:|---:|---:|",
  top_candidates$line,
  "",
  "## 常用胶质瘤细胞系对照",
  "",
  "| Cell line | DepMap ID | Subtype | LAP3 expr | LAP3 gene effect | mTORC1 score | BCAA score |",
  "|---|---|---|---:|---:|---:|---:|",
  common_readme$line,
  "",
  "## 关键输出",
  "",
  "- `tables/depmap_lap3_cell_line_dataset.csv`：整合后的细胞系级数据。",
  "- `tables/depmap_lap3_glioma_cell_lines.csv`：CNS/glioma 细胞系清单。",
  "- `tables/depmap_lap3_dependency_correlations.csv`：LAP3 dependency 与 expression、mTORC1、BCAA score 的相关性。",
  "- `tables/depmap_candidate_glioma_cell_lines.csv`：候选 glioma 细胞系排序。",
  "- `tables/depmap_common_glioma_cell_lines.csv`：常用 glioma 细胞系对照表。",
  "- `plots/Fig_DepMap_*.pdf/png`：第一版图表。",
  "",
  "## 第一版结论",
  "",
  "DepMap 模块的主要用途不是证明 LAP3 机制，而是回答两个实用问题：哪些 CNS/glioma 细胞系具有较高 LAP3 表达和较强 LAP3 dependency，以及这种 dependency 是否与 mTORC1/BCAA-leucine 转录状态同向。第一版结果显示，CNS/glioma 内 LAP3 dependency 与 mTORC1/BCAA score 的相关性较弱，不支持把 DepMap 依赖性作为 LAP3-mTORC1 机制的强证据。更合适的用法是筛选可操作细胞系，并把 DepMap 结果作为湿实验选型和风险提示。最终选型应结合实验室现有细胞系、培养稳定性、转染/感染效率和文献常用性共同决定。"
)
writeLines(readme, file.path(out_dir, "README_LAP3_DepMap.md"))

summary_obj <- list(
  coverage = coverage,
  lineage_summary = lineage_summary,
  group_summary = group_summary,
  correlations = cor_results,
  candidate_glioma_cell_lines = candidate,
  score_gene_coverage = score_gene_coverage
)
saveRDS(summary_obj, file.path(out_dir, "lap3_depmap_summary_tables.rds"))

cat("\nDone.\n")
