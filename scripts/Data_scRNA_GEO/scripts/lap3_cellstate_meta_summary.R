#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

set.seed(20260629)
data.table::setDTthreads(8)

project_dir <- "/home/lzb/glioma"
setwd(project_dir)

out_dir <- file.path("Data_scRNA_GEO", "results", "LAP3_CellState_MetaSummary")
tables_dir <- file.path(out_dir, "tables")
plots_dir <- file.path(out_dir, "plots")
source_dir <- file.path(out_dir, "source_data")
logs_dir <- file.path(out_dir, "logs")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(logs_dir, "lap3_cellstate_meta_summary.log")
cat("", file = log_file)
log_msg <- function(...) {
  text <- paste0(...)
  message(text)
  cat(text, "\n", file = log_file, append = TRUE)
}

log_msg("LAP3 CellState MetaSummary started: ", Sys.time())

read_csv <- function(path) {
  if (!file.exists(path)) stop("Missing input: ", path)
  data.table::fread(path)
}

clean_state <- function(x) {
  x <- gsub("^mean_", "", x)
  x <- fifelse(x == "neftel_npc_score", "NPC", x)
  x <- fifelse(x == "neftel_mes_score", "MES", x)
  x
}

classify_signal <- function(rho, p_adj, p_value, n_patients,
                            sensitivity_flag = NA_character_,
                            fdr_cut = 0.05,
                            trend_fdr_cut = 0.15) {
  out <- rep("not_evaluable", length(rho))
  ok <- !is.na(rho) & !is.na(n_patients) & n_patients >= 6
  out[ok] <- "not_reproduced"
  out[ok & !is.na(p_adj) & p_adj < fdr_cut] <- "supported_in_this_analysis"
  out[ok & !is.na(p_adj) & p_adj >= fdr_cut & p_adj < trend_fdr_cut] <- "exploratory_trend"
  out[ok & is.na(p_adj) & !is.na(p_value) & p_value < 0.05] <- "nominal_only"
  out[ok & !is.na(sensitivity_flag) & sensitivity_flag != "main"] <-
    paste0(out[ok & !is.na(sensitivity_flag) & sensitivity_flag != "main"], "_sensitivity")
  out
}

input_manifest <- data.table(
  role = c(
    "phase0_rules",
    "gse211376_continuous_state",
    "gse211376_within_state_pathway",
    "gse211376_rank_pathway",
    "gse211376_depth_adjusted_pathway",
    "gse211376_threshold_sensitivity",
    "gse278456_continuous_state",
    "gse278456_highconf_pathway",
    "gse278456_sensitivity_key",
    "gse278456_highconf_sensitivity",
    "gse131928_highconf_pathway",
    "gse131928_depth_adjusted_pathway",
    "gse131928_threshold_sensitivity",
    "gse131928_detection_depth_summary"
  ),
  path = c(
    "Data_scRNA_GEO/results/LAP3_CellState_Phase0/README.md",
    "Data_scRNA_GEO/results/GSE211376_LAP3_CellState/tables/gse211376_lap3_continuous_state_score_associations.csv",
    "Data_scRNA_GEO/results/GSE211376_LAP3_CellState/tables/gse211376_lap3_pathway_within_state_associations.csv",
    "Data_scRNA_GEO/results/GSE211376_LAP3_CellState/tables/gse211376_lap3_rank_pathway_within_state_associations.csv",
    "Data_scRNA_GEO/results/GSE211376_LAP3_CellState/tables/gse211376_lap3_pathway_depth_adjusted_associations.csv",
    "Data_scRNA_GEO/results/GSE211376_LAP3_CellState/tables/gse211376_lap3_pathway_threshold_sensitivity.csv",
    "Data_scRNA_GEO/results/GSE278456_LAP3_CellState/tables/gse278456_lap3_continuous_state_associations.csv",
    "Data_scRNA_GEO/results/GSE278456_LAP3_CellState/tables/gse278456_highconf_state_pathway_associations.csv",
    "Data_scRNA_GEO/results/GSE278456_LAP3_CellState_Sensitivity/tables/gse278456_key_sensitivity_summary.csv",
    "Data_scRNA_GEO/results/GSE278456_LAP3_CellState_Sensitivity/tables/gse278456_highconf_pathway_malignant_rule_sensitivity.csv",
    "Data_scRNA_GEO/results/GSE131928_LAP3_CellState/tables/gse131928_highconf_state_pathway_associations.csv",
    "Data_scRNA_GEO/results/GSE131928_LAP3_CellState/tables/gse131928_pathway_depth_adjusted_associations.csv",
    "Data_scRNA_GEO/results/GSE131928_LAP3_CellState/tables/gse131928_pathway_threshold_sensitivity.csv",
    "Data_scRNA_GEO/results/GSE131928_LAP3_CellState/tables/gse131928_platform_detection_depth_summary.csv"
  )
)
input_manifest[, exists := file.exists(path)]
if (any(!input_manifest$exists)) {
  print(input_manifest[exists == FALSE])
  stop("One or more MetaSummary inputs are missing.")
}
fwrite(input_manifest, file.path(source_dir, "metasummary_input_manifest.csv"))

g211_cont <- read_csv(input_manifest[role == "gse211376_continuous_state", path])
g211_path <- read_csv(input_manifest[role == "gse211376_within_state_pathway", path])
g211_rank <- read_csv(input_manifest[role == "gse211376_rank_pathway", path])
g211_depth <- read_csv(input_manifest[role == "gse211376_depth_adjusted_pathway", path])
g211_thresh <- read_csv(input_manifest[role == "gse211376_threshold_sensitivity", path])

g278_cont <- read_csv(input_manifest[role == "gse278456_continuous_state", path])
g278_path <- read_csv(input_manifest[role == "gse278456_highconf_pathway", path])
g278_key <- read_csv(input_manifest[role == "gse278456_sensitivity_key", path])
g278_sens_path <- read_csv(input_manifest[role == "gse278456_highconf_sensitivity", path])

g131_path <- read_csv(input_manifest[role == "gse131928_highconf_pathway", path])
g131_depth <- read_csv(input_manifest[role == "gse131928_depth_adjusted_pathway", path])
g131_thresh <- read_csv(input_manifest[role == "gse131928_threshold_sensitivity", path])
g131_detect <- read_csv(input_manifest[role == "gse131928_detection_depth_summary", path])

state_continuum <- rbindlist(list(
  g211_cont[signature %in% c("AC", "OPC", "neftel_npc_score", "neftel_mes_score"),
            .(dataset = "GSE211376",
              analysis_class = "IDHwt_GBM_Ruiz_NH",
              malignant_rule = "author_malignant_state",
              estimand,
              state = clean_state(signature),
              n_patients,
              n_patient_state,
              spearman_rho,
              ci_low,
              ci_high,
              p_value,
              p_adj = p_adj_BH_exploratory_state,
              source_result = "continuous_state_patient_residual")],
  g278_cont[analysis_class %in% c("GBM_grade4_IDHwt", "LGG_grade2_3_IDHmut"),
            .(dataset = "GSE278456",
              analysis_class,
              malignant_rule = "main_marker_qc",
              estimand,
              state = clean_state(signature),
              n_patients,
              n_patient_state = NA_integer_,
              spearman_rho,
              ci_low,
              ci_high,
              p_value,
              p_adj = p_adj_BH,
              source_result = "continuous_state_across_patients")]
), fill = TRUE)
state_continuum[, fdr_family := "state_continuum_exploratory"]
state_continuum[, evidence_level := classify_signal(
  spearman_rho, p_adj, p_value, n_patients
)]
setcolorder(state_continuum, c(
  "dataset", "analysis_class", "malignant_rule", "estimand", "state",
  "n_patients", "n_patient_state", "spearman_rho", "ci_low", "ci_high",
  "p_value", "p_adj", "fdr_family", "evidence_level", "source_result"
))
fwrite(state_continuum, file.path(tables_dir, "state_continuum_summary.csv"))

within_state_pathway <- rbindlist(list(
  g211_path[, .(dataset = "GSE211376",
                analysis_class = "IDHwt_GBM_Ruiz_NH",
                platform = "snRNAseq",
                malignant_rule = "author_malignant_state",
                confidence_rule = threshold,
                score_type,
                state = author_state,
                pathway,
                fdr_family,
                n_patients,
                spearman_rho,
                ci_low,
                ci_high,
                p_value,
                p_adj = p_adj_BH,
                lopo_min_rho,
                lopo_max_rho,
                source_result = "main_pseudobulk_state_pathway")],
  g278_path[, .(dataset = "GSE278456",
                analysis_class,
                platform = "10X_like_snRNAseq",
                malignant_rule = "main_marker_qc",
                confidence_rule = "gse131928_10x_main_highconf",
                score_type = "cell_cache_zmean_patient_state",
                state = dominant_state,
                pathway,
                fdr_family,
                n_patients,
                spearman_rho,
                ci_low,
                ci_high,
                p_value,
                p_adj = p_adj_BH,
                lopo_min_rho,
                lopo_max_rho,
                source_result = "external_threshold_main_pathway")],
  g131_path[, .(dataset = "GSE131928",
                analysis_class = "IDHwt_GBM_adult",
                platform,
                malignant_rule = "tumor_cells_processed_TPM",
                confidence_rule = "platform_main_highconf_n20",
                score_type = "log1p_TPM_platform_zmean",
                state = dominant_state,
                pathway,
                fdr_family,
                n_patients,
                spearman_rho,
                ci_low,
                ci_high,
                p_value,
                p_adj = p_adj_BH,
                lopo_min_rho,
                lopo_max_rho,
                source_result = "platform_main_pathway")]
), fill = TRUE)
within_state_pathway[, evidence_level := classify_signal(
  spearman_rho, p_adj, p_value, n_patients
)]
within_state_pathway[, lopo_direction_stable := fifelse(
  is.na(lopo_min_rho) | is.na(lopo_max_rho), NA,
  (lopo_min_rho > 0 & lopo_max_rho > 0) | (lopo_min_rho < 0 & lopo_max_rho < 0)
)]
within_state_pathway[, primary_readout := fdr_family == "primary"]
fwrite(within_state_pathway, file.path(tables_dir, "within_state_pathway_summary.csv"))

rank_sensitivity <- g211_rank[, .(
  dataset = "GSE211376",
  analysis_class = "IDHwt_GBM_Ruiz_NH",
  sensitivity_type = "rank_based_pathway_score",
  state = author_state,
  pathway,
  fdr_family,
  n_patients,
  spearman_rho,
  p_value,
  p_adj = p_adj_BH
)]

g211_depth_sensitivity <- g211_depth[, .(
  dataset = "GSE211376",
  analysis_class = "IDHwt_GBM_Ruiz_NH",
  sensitivity_type = "library_size_ncells_depth_adjusted",
  state = author_state,
  pathway,
  fdr_family,
  n_patients,
  spearman_rho = depth_adjusted_spearman_rho,
  p_value,
  p_adj = p_adj_BH
)]

g131_depth_sensitivity <- g131_depth[, .(
  dataset = "GSE131928",
  analysis_class = paste("IDHwt_GBM_adult", platform, sep = "_"),
  sensitivity_type = "detected_genes_depth_adjusted",
  state = dominant_state,
  pathway,
  fdr_family,
  n_patients,
  spearman_rho,
  p_value,
  p_adj = p_adj_BH
)]

g278_cont_sensitivity <- g278_key[evidence_type == "continuous_state", .(
  dataset = "GSE278456",
  analysis_class,
  sensitivity_type = paste0("malignant_rule_", rule),
  state = clean_state(state_or_signature),
  pathway = NA_character_,
  fdr_family = "state_continuum_exploratory",
  n_patients,
  spearman_rho,
  p_value,
  p_adj = p_adj_BH
)]

g278_path_sensitivity <- g278_sens_path[, .(
  dataset = "GSE278456",
  analysis_class,
  sensitivity_type = paste(rule, highconf_rule, sep = "_"),
  state = dominant_state,
  pathway,
  fdr_family,
  n_patients,
  spearman_rho,
  p_value,
  p_adj = p_adj_BH
)]

g131_threshold_sensitivity <- g131_thresh[!is.na(spearman_rho), .(
  dataset = "GSE131928",
  analysis_class = paste("IDHwt_GBM_adult", platform, sep = "_"),
  sensitivity_type = paste(confidence_rule, paste0("n", minimum_cells), sep = "_"),
  state = dominant_state,
  pathway,
  fdr_family,
  n_patients,
  spearman_rho,
  p_value,
  p_adj = p_adj_BH
)]

sensitivity_summary <- rbindlist(list(
  rank_sensitivity,
  g211_depth_sensitivity,
  g131_depth_sensitivity,
  g278_cont_sensitivity,
  g278_path_sensitivity,
  g131_threshold_sensitivity
), fill = TRUE)
sensitivity_summary[, evidence_level := classify_signal(
  spearman_rho, p_adj, p_value, n_patients, sensitivity_type
)]
fwrite(sensitivity_summary, file.path(tables_dir, "sensitivity_evidence_summary.csv"))

primary_pathway <- within_state_pathway[fdr_family == "primary"]
primary_pathway[, positive_primary := !is.na(spearman_rho) & spearman_rho > 0]
primary_pathway[, fdr_supported := !is.na(p_adj) & p_adj < 0.05]
primary_pathway[, trend_or_better := !is.na(p_adj) & p_adj < 0.15]

pathway_recurrence <- primary_pathway[, .(
  n_analyses = .N,
  n_evaluable = sum(!is.na(spearman_rho) & n_patients >= 6),
  n_positive = sum(positive_primary, na.rm = TRUE),
  n_fdr_supported = sum(fdr_supported, na.rm = TRUE),
  n_trend_or_better = sum(trend_or_better, na.rm = TRUE),
  median_rho = median(spearman_rho, na.rm = TRUE),
  max_rho = suppressWarnings(max(spearman_rho, na.rm = TRUE)),
  min_rho = suppressWarnings(min(spearman_rho, na.rm = TRUE))
), by = .(pathway)]
pathway_recurrence[is.infinite(max_rho), max_rho := NA_real_]
pathway_recurrence[is.infinite(min_rho), min_rho := NA_real_]
pathway_recurrence[, cross_dataset_interpretation := fifelse(
  n_fdr_supported >= 2,
  "cross_dataset_supported",
  fifelse(
    n_trend_or_better >= 2 & n_positive >= 2,
    "recurrent_exploratory_trend",
    fifelse(n_positive >= 2, "directional_but_not_significant", "not_reproduced")
  )
)]
fwrite(pathway_recurrence, file.path(tables_dir, "primary_pathway_recurrence_summary.csv"))

interpretation_matrix <- data.table(
  question = c(
    "Does LAP3 associate with malignant-state continuum?",
    "Does LAP3 show robust within-state mTORC1 coupling?",
    "Does LAP3 show robust within-state BCAA/leucine coupling?",
    "Is GSE278456 AC/mTORC1 signal enough for a main claim?",
    "Can scRNA prove LAP3-leucine-mTORC1 causality?"
  ),
  answer = c(
    "Yes, as an exploratory state/composition signal. GSE211376 supports AC/MES-like and low NPC-like direction; GSE278456 supports GBM AC and LGG MES directions with rule-sensitive FDR.",
    "No. Main within-state analyses are not FDR-robust across datasets. GSE278456 GBM AC mTORC1 becomes significant only in strict high-confidence sensitivity subsets.",
    "No. BCAA/leucine readout is not consistently supported in state-specific scRNA analyses.",
    "No. It is a threshold-sensitive exploratory signal, useful for Figure 5 supplement or a cautious panel, not a standalone mechanism claim.",
    "No. These are cross-sectional expression associations and state summaries. Wet-lab phospho/metabolite/rescue experiments remain necessary."
  ),
  evidence_grade = c(
    "exploratory_supported",
    "not_robust",
    "not_reproduced",
    "exploratory_threshold_sensitive",
    "not_addressable_by_current_data"
  ),
  recommended_wording = c(
    "LAP3 was associated with AC/MES-like malignant-state continua.",
    "Within-state LAP3-mTORC1 transcriptional coupling was not robustly reproduced.",
    "Within-state LAP3-BCAA/leucine transcriptional coupling was not supported.",
    "A GBM AC-like mTORC1 trend was observed in sensitivity analyses but remained threshold-sensitive.",
    "Single-cell data nominate a cellular context but do not establish LAP3-dependent mTORC1 activation."
  )
)
fwrite(interpretation_matrix, file.path(tables_dir, "meta_interpretation_matrix.csv"))

dataset_status <- data.table(
  dataset = c("GSE211376", "GSE278456", "GSE131928"),
  main_strength = c(
    "LAP3-state continuum association in author malignant states",
    "Independent malignant-rule sensitivity for AC/MES continuum; AC/mTORC1 sensitivity signal",
    "Platform calibration and adult-only state/pathway sensitivity"
  ),
  main_limitation = c(
    "Only 11 patients; within-state pathway n is small; rank/depth sensitivity weakens pathway inference",
    "No author malignant state labels; hard state threshold migrates imperfectly from GSE131928",
    "Processed TPM, platform repeats, depth adjustment undermines Smart-seq2 trends"
  ),
  meta_role = c(
    "discovery of state preference",
    "external projection and malignant-rule sensitivity",
    "threshold calibration and platform sensitivity"
  ),
  final_evidence_level = c(
    "state_signal_supported_pathway_not_robust",
    "state_signal_exploratory_pathway_threshold_sensitive",
    "platform_trend_not_depth_robust"
  )
)
fwrite(dataset_status, file.path(tables_dir, "dataset_role_and_limitation_summary.csv"))

plot_state <- state_continuum[!is.na(spearman_rho) & state %in% c("AC", "MES", "NPC", "OPC")]
plot_state[, label := paste(dataset, analysis_class, sep = "\n")]
plot_state[, state := factor(state, levels = c("AC", "MES", "NPC", "OPC"))]

p1 <- ggplot(plot_state, aes(x = state, y = spearman_rho, color = dataset)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey55") +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.12, linewidth = 0.45,
                position = position_dodge(width = 0.55), na.rm = TRUE) +
  geom_point(position = position_dodge(width = 0.55), size = 2.2) +
  facet_wrap(~analysis_class, ncol = 1) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  labs(x = NULL, y = "Spearman rho", color = "Dataset",
       title = "LAP3 association with malignant-state continuum") +
  theme_classic(base_size = 10) +
  theme(legend.position = "bottom", strip.background = element_blank())

ggsave(file.path(plots_dir, "state_continuum_effects.pdf"), p1, width = 7, height = 7)
ggsave(file.path(plots_dir, "state_continuum_effects.svg"), p1, width = 7, height = 7)
ggsave(file.path(plots_dir, "state_continuum_effects.png"), p1, width = 7, height = 7, dpi = 300)

plot_path <- primary_pathway[!is.na(spearman_rho) & n_patients >= 6]
plot_path[, analysis_label := paste(dataset, analysis_class, state, sep = " | ")]
plot_path[, pathway := factor(pathway, levels = c("HALLMARK_MTORC1_SIGNALING", "LEUCINE_BCAA_CORE"))]
plot_path[, evidence_level := factor(evidence_level, levels = c(
  "supported_in_this_analysis", "exploratory_trend", "nominal_only",
  "not_reproduced", "not_evaluable"
))]

p2 <- ggplot(plot_path, aes(x = spearman_rho, y = reorder(analysis_label, spearman_rho),
                            shape = pathway, color = evidence_level)) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey55") +
  geom_point(size = 2.3, alpha = 0.9) +
  scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  scale_color_manual(values = c(
    supported_in_this_analysis = "#0072B2",
    exploratory_trend = "#D55E00",
    nominal_only = "#CC79A7",
    not_reproduced = "#666666",
    not_evaluable = "#BBBBBB"
  ), drop = FALSE) +
  labs(x = "Spearman rho", y = NULL, shape = "Primary pathway",
       color = "Evidence", title = "Within-state LAP3-pathway associations") +
  theme_classic(base_size = 9) +
  theme(legend.position = "bottom")

ggsave(file.path(plots_dir, "within_state_primary_pathway_effects.pdf"), p2, width = 8, height = 7)
ggsave(file.path(plots_dir, "within_state_primary_pathway_effects.svg"), p2, width = 8, height = 7)
ggsave(file.path(plots_dir, "within_state_primary_pathway_effects.png"), p2, width = 8, height = 7, dpi = 300)

readme <- c(
  "# LAP3 CellState MetaSummary",
  "",
  paste0("生成时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## 目的",
  "",
  "本结果包整合 GSE211376、GSE278456 和 GSE131928 的 LAP3 malignant cell-state 严格分析结果。MetaSummary 不重新读取大对象，也不重新估计单细胞分数；它只消费已经完成的 CSV 结果表，统一 estimand、effect size、敏感性风险和解释边界。",
  "",
  "核心区分：",
  "",
  "- 状态偏好 / state continuum：LAP3 是否更偏 AC/MES/NPC/OPC 等恶性状态。",
  "- 状态内通路关联：同一状态内，不同患者的 LAP3 是否与 mTORC1/BCAA/translation readout 相关。",
  "",
  "## 输入",
  "",
  "`source_data/metasummary_input_manifest.csv` 记录了所有输入表。三套数据集的原始脚本和 README 保持不变。",
  "",
  "## 关键输出",
  "",
  "- `tables/state_continuum_summary.csv`：LAP3 与连续恶性状态分数的统一 effect-size 表。",
  "- `tables/within_state_pathway_summary.csv`：状态内 LAP3-pathway 患者级关联统一表。",
  "- `tables/sensitivity_evidence_summary.csv`：rank/depth/阈值/恶性定义敏感性证据。",
  "- `tables/primary_pathway_recurrence_summary.csv`：mTORC1 与 BCAA 主假设的跨分析复现概览。",
  "- `tables/meta_interpretation_matrix.csv`：可直接写入项目结论的判断矩阵。",
  "- `tables/dataset_role_and_limitation_summary.csv`：每个数据集在综合证据链中的角色和限制。",
  "- `plots/state_continuum_effects.{pdf,svg,png}`。",
  "- `plots/within_state_primary_pathway_effects.{pdf,svg,png}`。",
  "",
  "## 综合结论",
  "",
  "1. 三数据集综合后，LAP3 与恶性细胞状态 continuum 的关系有较一致的探索性支持：GSE211376 偏 AC/MES-like、低 NPC-like；GSE278456 中 GBM 偏 AC、LGG 偏 MES，且方向对恶性定义有一定稳定性。",
  "2. 状态内 LAP3-mTORC1/BCAA 关联没有形成跨数据集、跨阈值、跨深度调整均稳健的证据。GSE278456 GBM AC-like mTORC1 信号在 strict high-conf 敏感性子集中可达 FDR<0.05，但 main high-conf 不显著，且 BCAA 不支持。",
  "3. GSE131928 Smart-seq2 中 AC/NPC 的 mTORC1 或 BCAA 有名义趋势，但成人限定、平台内、深度调整后不稳定；10X 可估计状态少。",
  "4. 因此 Figure 5 应定位为 cell-state/context/composition 证据，而不是 malignant-cell-intrinsic LAP3-mTORC1/BCAA 机制证明。",
  "",
  "## 推荐写法",
  "",
  "可以写：",
  "",
  "> LAP3 was associated with AC/MES-like malignant-state continua, while within-state LAP3-mTORC1/BCAA transcriptional coupling was not robustly reproduced across datasets and sensitivity analyses.",
  "",
  "不应写：",
  "",
  "> Single-cell data demonstrate that LAP3 activates mTORC1/BCAA metabolism within malignant glioma cells.",
  "",
  "## 下一步",
  "",
  "1. 用本目录输出冻结 Figure 5 的主图/补图结构。",
  "2. 若 Figure 5 作主图，优先展示 state continuum 与 cell-source/context；把状态内 pathway 结果放为谨慎补图或结果边界。",
  "3. 后续如补 CNV/inferCNV 或新数据集，应追加到本 MetaSummary，而不是替换已有负结果。",
  "",
  "## 复现命令",
  "",
  "```bash",
  "Rscript --vanilla Data_scRNA_GEO/scripts/lap3_cellstate_meta_summary.R",
  "```"
)
writeLines(readme, file.path(out_dir, "README.md"))

log_msg("Rows written:")
log_msg("state_continuum_summary: ", nrow(state_continuum))
log_msg("within_state_pathway_summary: ", nrow(within_state_pathway))
log_msg("sensitivity_evidence_summary: ", nrow(sensitivity_summary))
log_msg("MetaSummary completed: ", Sys.time())
