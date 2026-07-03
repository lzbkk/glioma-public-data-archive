#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

setwd("/home/lzb/glioma/Data_Bulk_TCGA/Data_Merged")

depmap_dir <- file.path("results", "LAP3_DepMap")
table_dir <- file.path(depmap_dir, "tables")
export_dir <- file.path(depmap_dir, "exports")
dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)

candidate_file <- file.path(table_dir, "depmap_candidate_glioma_cell_lines.csv")
common_file <- file.path(table_dir, "depmap_common_glioma_cell_lines.csv")
cor_file <- file.path(table_dir, "depmap_lap3_dependency_correlations.csv")
all_file <- file.path(table_dir, "depmap_lap3_cell_line_dataset.csv")
stopifnot(file.exists(candidate_file), file.exists(common_file), file.exists(cor_file), file.exists(all_file))

candidate <- read_csv(candidate_file, show_col_types = FALSE)
common <- read_csv(common_file, show_col_types = FALSE)
correlations <- read_csv(cor_file, show_col_types = FALSE)
all_lines <- read_csv(all_file, show_col_types = FALSE)

fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), NA_character_, formatC(x, digits = digits, format = "fg"))
}

classify_use <- function(display_cell_line, expr, dependency, mtorc1, bcaa, is_common) {
  if (is.na(dependency)) return("insufficient DepMap data")
  if (display_cell_line == "U-87 MG") return("expression-positive comparator; not KD-priority")
  if (dependency >= 0.15 && expr >= 5.8 && (mtorc1 >= 0.2 || bcaa >= 0.2)) {
    return("KD priority")
  }
  if (dependency >= 0.08 && expr >= 5.8 && (mtorc1 >= 0.15 || bcaa >= 0.15)) {
    return("KD candidate")
  }
  if (expr < 5.8 && dependency >= 0.08) {
    return("OE or KD pilot candidate")
  }
  if (dependency < 0.05) {
    return("comparator or caution")
  }
  if (is_common) return("common-line pilot candidate")
  "exploratory candidate"
}

priority_tier <- function(recommended_use, is_common, candidate_score) {
  case_when(
    recommended_use == "KD priority" & is_common ~ "Tier 1: wet-lab priority",
    recommended_use == "KD priority" ~ "Tier 2: data-driven, availability check",
    recommended_use %in% c("KD candidate", "OE or KD pilot candidate", "common-line pilot candidate") ~ "Tier 2: pilot/backup",
    TRUE ~ "Tier 3: comparator/caution"
  )
}

candidate_short <- candidate %>%
  mutate(source_candidate_rank = row_number()) %>%
  slice_head(n = 25) %>%
  mutate(in_top_candidate_table = TRUE)

common_short <- common %>%
  mutate(in_common_table = TRUE)

delivery <- full_join(
  candidate_short,
  common_short %>%
    select(
      DepMap_ID, in_common_table,
      common_display_cell_line = display_cell_line,
      common_Subtype = lineage_sub_subtype,
      common_lineage_sub_subtype = lineage_sub_subtype,
      common_LAP3_expression_log2_tpm1 = LAP3_expression_log2_tpm1,
      common_LAP3_gene_effect = LAP3_gene_effect,
      common_LAP3_dependency = LAP3_dependency,
      common_HALLMARK_MTORC1_SIGNALING_score = HALLMARK_MTORC1_SIGNALING_score,
      common_BCAA_LEUCINE_TRANSPORT_METABOLISM_score = BCAA_LEUCINE_TRANSPORT_METABOLISM_score,
      common_MTORC1_READOUT_CORE_score = MTORC1_READOUT_CORE_score
    ),
  by = "DepMap_ID"
) %>%
  left_join(
    all_lines %>%
      select(
        DepMap_ID,
        all_default_growth_pattern = default_growth_pattern,
        all_source = source,
        all_Subtype = Subtype,
        all_lineage_sub_subtype = lineage_sub_subtype
      ),
    by = "DepMap_ID"
  ) %>%
  mutate(
    in_top_candidate_table = ifelse(is.na(in_top_candidate_table), FALSE, in_top_candidate_table),
    in_common_table = ifelse(is.na(in_common_table), FALSE, in_common_table),
    display_cell_line = coalesce(display_cell_line, common_display_cell_line),
    Subtype = coalesce(Subtype, common_Subtype, all_Subtype),
    lineage_sub_subtype = coalesce(lineage_sub_subtype, common_lineage_sub_subtype, all_lineage_sub_subtype),
    default_growth_pattern = coalesce(default_growth_pattern, all_default_growth_pattern),
    source = coalesce(source, all_source),
    LAP3_expression_log2_tpm1 = coalesce(LAP3_expression_log2_tpm1, common_LAP3_expression_log2_tpm1),
    LAP3_gene_effect = coalesce(LAP3_gene_effect, common_LAP3_gene_effect),
    LAP3_dependency = coalesce(LAP3_dependency, common_LAP3_dependency),
    HALLMARK_MTORC1_SIGNALING_score = coalesce(HALLMARK_MTORC1_SIGNALING_score, common_HALLMARK_MTORC1_SIGNALING_score),
    BCAA_LEUCINE_TRANSPORT_METABOLISM_score = coalesce(BCAA_LEUCINE_TRANSPORT_METABOLISM_score, common_BCAA_LEUCINE_TRANSPORT_METABOLISM_score),
    MTORC1_READOUT_CORE_score = coalesce(MTORC1_READOUT_CORE_score, common_MTORC1_READOUT_CORE_score),
    wet_lab_common = display_cell_line %in% c("LN-18", "T98G", "U-251 MG", "A-172", "U-87 MG", "U-118 MG"),
    recommended_use = mapply(
      classify_use,
      display_cell_line,
      LAP3_expression_log2_tpm1,
      LAP3_dependency,
      HALLMARK_MTORC1_SIGNALING_score,
      BCAA_LEUCINE_TRANSPORT_METABOLISM_score,
      wet_lab_common
    ),
    priority = priority_tier(recommended_use, wet_lab_common, candidate_score),
    evidence_summary = case_when(
      display_cell_line == "LN-18" ~ "Common GBM line; high LAP3 expression, high mTORC1/BCAA scores, modest LAP3 dependency.",
      display_cell_line == "T98G" ~ "Common GBM line; high LAP3 expression, stronger dependency among common lines, high BCAA score.",
      display_cell_line == "U-251 MG" ~ "Common astrocytoma/GBM-use line; two DepMap entries, one shows stronger LAP3 dependency and higher BCAA score.",
      display_cell_line == "A-172" ~ "Common GBM line; moderate dependency but lower LAP3 expression and negative BCAA score.",
      display_cell_line == "U-87 MG" ~ "Common line; high LAP3 expression but LAP3 gene effect near zero and negative BCAA score.",
      display_cell_line == "U-118 MG" ~ "Common line; LAP3 gene effect near zero, better as comparator than KD-priority line.",
      in_top_candidate_table ~ "Data-driven candidate from composite LAP3 expression/dependency/mTORC1/BCAA ranking.",
      TRUE ~ "Common glioma comparator line."
    ),
    caution = case_when(
      display_cell_line == "U-87 MG" ~ "Do not expect strong LAP3 KD proliferation phenotype from DepMap alone.",
      display_cell_line == "U-118 MG" ~ "Low/negative dependency; use as comparator or only if lab availability dictates.",
      !wet_lab_common & source %in% c("Academic lab", "HSRRB", "RIKEN", "JCRB") ~ "Availability and handling should be confirmed before recommendation.",
      TRUE ~ "Final choice should consider lab availability, culture stability, transduction efficiency, and authentication."
    )
  ) %>%
  arrange(
    factor(priority, levels = c("Tier 1: wet-lab priority", "Tier 2: pilot/backup", "Tier 2: data-driven, availability check", "Tier 3: comparator/caution")),
    desc(wet_lab_common),
    desc(coalesce(candidate_score, -Inf))
  ) %>%
  transmute(
    priority,
    recommended_use,
    display_cell_line,
    DepMap_ID,
    Subtype,
    lineage_sub_subtype,
    default_growth_pattern,
    source,
    wet_lab_common,
    in_top_candidate_table,
    source_candidate_rank,
    candidate_score,
    LAP3_expression_log2_tpm1,
    LAP3_gene_effect,
    LAP3_dependency,
    HALLMARK_MTORC1_SIGNALING_score,
    BCAA_LEUCINE_TRANSPORT_METABOLISM_score,
    MTORC1_READOUT_CORE_score,
    evidence_summary,
    caution
  )

write_csv(delivery, file.path(export_dir, "candidate_cell_lines.csv"))

top_for_readme <- bind_rows(
  delivery %>% filter(priority == "Tier 1: wet-lab priority"),
  delivery %>% filter(wet_lab_common, priority != "Tier 1: wet-lab priority"),
  delivery %>% filter(priority == "Tier 2: data-driven, availability check") %>% arrange(desc(coalesce(candidate_score, -Inf))) %>% slice_head(n = 5)
) %>%
  distinct(DepMap_ID, .keep_all = TRUE) %>%
  slice_head(n = 12)

cns_glioma_cor <- correlations %>%
  filter(subset %in% c("CNS_glioma", "CNS_glioma_GBM"))

readme <- c(
  "# LAP3 DepMap Candidate Cell Line Delivery Package",
  "",
  paste("生成时间：", as.character(Sys.time())),
  "",
  "## 目的",
  "",
  "本交付包把 `LAP3_DepMap` 第一版分析结果整理成湿实验团队可直接讨论的候选细胞系清单。它用于实验选型和风险提示，不用于证明 LAP3-亮氨酸-mTORC1 机制因果。",
  "",
  "## 输入",
  "",
  "```text",
  "results/LAP3_DepMap/tables/depmap_candidate_glioma_cell_lines.csv",
  "results/LAP3_DepMap/tables/depmap_common_glioma_cell_lines.csv",
  "results/LAP3_DepMap/tables/depmap_lap3_dependency_correlations.csv",
  "```",
  "",
  "## 输出",
  "",
  "```text",
  "results/LAP3_DepMap/exports/candidate_cell_lines.csv",
  "results/LAP3_DepMap/exports/candidate_cell_lines_README.md",
  "```",
  "",
  "## 字段解释",
  "",
  "- `priority`：综合湿实验常用性、LAP3 expression、LAP3 dependency、mTORC1/BCAA score 后给出的交付等级。",
  "- `recommended_use`：建议用途，例如 KD priority、KD candidate、OE or KD pilot candidate、comparator or caution。",
  "- `LAP3_gene_effect`：DepMap CRISPR gene effect，越负代表敲除后适应度下降越明显。",
  "- `LAP3_dependency`：`-LAP3_gene_effect`，越高代表越依赖 LAP3。",
  "- `wet_lab_common`：是否属于当前湿实验讨论中常见的 glioma 细胞系。",
  "- `evidence_summary` 和 `caution`：给湿实验讨论用的简短解释。",
  "",
  "## 推荐优先级",
  "",
  "| Priority | Use | Cell line | DepMap ID | LAP3 expr | Gene effect | mTORC1 | BCAA | Summary |",
  "|---|---|---|---|---:|---:|---:|---:|---|",
  paste0(
    "| ", top_for_readme$priority,
    " | ", top_for_readme$recommended_use,
    " | ", top_for_readme$display_cell_line,
    " | ", top_for_readme$DepMap_ID,
    " | ", fmt_num(top_for_readme$LAP3_expression_log2_tpm1),
    " | ", fmt_num(top_for_readme$LAP3_gene_effect),
    " | ", fmt_num(top_for_readme$HALLMARK_MTORC1_SIGNALING_score),
    " | ", fmt_num(top_for_readme$BCAA_LEUCINE_TRANSPORT_METABOLISM_score),
    " | ", top_for_readme$evidence_summary,
    " |"
  ),
  "",
  "## 当前建议",
  "",
  "- 优先讨论 `T98G`、`LN-18` 和 `U-251 MG`：三者兼具常用性和较可解释的 LAP3/mTORC1/BCAA 背景，适合优先设计 LAP3 knockdown 或 knockdown + mTORC1 readout。",
  "- `A-172` 可作为 pilot/backup：LAP3 dependency 不弱，但 BCAA score 为负，机制 readout 预期需要谨慎。",
  "- `U-87 MG` 可作为表达阳性但 dependency 较弱的比较对象，不建议作为首选 LAP3 knockdown 表型验证细胞系。",
  "- `8-MG-BA`、`H4`、`SF126`、`LNZ308` 等是数据驱动候选，若实验室可获得，可作为备选或后续扩展。",
  "",
  "## 解释边界",
  "",
  "- CNS/glioma 子集中 LAP3 dependency 与 mTORC1/BCAA score 相关性较弱，因此不能把 DepMap 结果解释为机制证明。",
  "- DepMap 主要反映二维培养中的细胞增殖/适应度依赖，不能覆盖侵袭、干性、代谢适应和体内微环境。",
  "- 最终细胞系选择必须结合细胞系来源、STR 鉴定、污染状态、培养稳定性、转染/感染效率和实验室现有条件。",
  "",
  "## 相关性提示",
  "",
  "| Subset | Variable | n | Spearman rho | P value |",
  "|---|---|---:|---:|---:|",
  paste0(
    "| ", cns_glioma_cor$subset,
    " | ", cns_glioma_cor$x,
    " | ", cns_glioma_cor$n,
    " | ", fmt_num(cns_glioma_cor$spearman_rho),
    " | ", fmt_num(cns_glioma_cor$p_value),
    " |"
  )
)

writeLines(readme, file.path(export_dir, "candidate_cell_lines_README.md"))

cat("Wrote delivery package to:", export_dir, "\n")
