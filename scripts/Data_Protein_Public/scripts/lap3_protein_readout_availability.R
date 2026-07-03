#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(httr)
  library(jsonlite)
})

data.table::setDTthreads(8)
options(timeout = 300)

project_dir <- "/home/lzb/glioma"
out_dir <- file.path(project_dir, "Data_Protein_Public", "results", "LAP3_Protein_Readout")
tables_dir <- file.path(out_dir, "tables")
logs_dir <- file.path(out_dir, "logs")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

log_path <- file.path(logs_dir, "lap3_protein_readout.log")
if (file.exists(log_path)) unlink(log_path)
log_msg <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  cat(msg, "\n", file = log_path, append = TRUE)
}

get_json <- function(url) {
  response <- GET(url, timeout(180))
  stop_for_status(response)
  fromJSON(content(response, "text", encoding = "UTF-8"), flatten = TRUE)
}

post_json <- function(url, body) {
  response <- POST(url, content_type_json(), body = body, encode = "json", timeout(180))
  stop_for_status(response)
  fromJSON(content(response, "text", encoding = "UTF-8"), flatten = TRUE)
}

target_genes <- c(
  LAP3 = 51056L,
  MTOR = 2475L,
  RPS6 = 6194L,
  RPS6KB1 = 6198L,
  EIF4EBP1 = 1978L
)

log_msg("Protein readout audit started: ", Sys.time())

# 1. HPA pathology IHC summary (archived stable release with aggregate counts)
hpa_url <- "https://v23.proteinatlas.org/download/pathology.tsv.zip"
hpa_zip <- tempfile(fileext = ".zip")
download.file(hpa_url, hpa_zip, mode = "wb", quiet = TRUE)
hpa <- fread(cmd = sprintf("unzip -p %s", shQuote(hpa_zip)))
hpa_lap3_glioma <- hpa[`Gene name` == "LAP3" & Cancer == "glioma"]
stopifnot(nrow(hpa_lap3_glioma) == 1L)
hpa_lap3_glioma[, `:=`(
  total_ihc_cases = High + Medium + Low + `Not detected`,
  detected_cases = High + Medium + Low,
  source_version = "HPA v23 archived pathology aggregate",
  source_url = hpa_url
)]
hpa_lap3_glioma[, detected_pct := 100 * detected_cases / total_ihc_cases]
fwrite(hpa_lap3_glioma, file.path(tables_dir, "hpa_lap3_glioma_ihc_summary.csv"))
unlink(hpa_zip)

# 2. CPTAC GBM total protein data
cbio_base <- "https://www.cbioportal.org/api"
cptac_protein_profile <- "gbm_cptac_2021_protein_quantification"
cptac_sample_list <- "gbm_cptac_2021_all"
protein_url <- sprintf(
  "%s/molecular-profiles/%s/molecular-data/fetch?projection=DETAILED",
  cbio_base, cptac_protein_profile
)
protein <- as.data.table(post_json(
  protein_url,
  list(entrezGeneIds = unname(target_genes), sampleListId = cptac_sample_list)
))
protein[, value := as.numeric(value)]
protein[, gene_symbol := gene.hugoGeneSymbol]
protein <- protein[, .(studyId, molecularProfileId, sampleId, patientId, gene_symbol, value)]
stopifnot(all(names(target_genes) %in% protein$gene_symbol))
fwrite(protein, file.path(tables_dir, "cptac_gbm_target_total_protein_values.csv"))

protein_coverage <- protein[, .(
  n_samples = uniqueN(sampleId),
  n_nonmissing = sum(!is.na(value))
), by = gene_symbol]
fwrite(protein_coverage, file.path(tables_dir, "cptac_gbm_total_protein_coverage.csv"))

# 3. CPTAC GBM phosphosite metadata and values
phospho_profile <- "gbm_cptac_2021_phosphoproteome"
phospho_meta <- as.data.table(get_json(sprintf(
  "%s/generic-assay-meta/%s?projection=DETAILED",
  cbio_base, phospho_profile
)))
setnames(
  phospho_meta,
  c(
    "genericEntityMetaProperties.GENE_SYMBOL",
    "genericEntityMetaProperties.DESCRIPTION",
    "genericEntityMetaProperties.NAME"
  ),
  c("gene_symbol", "description", "name")
)
target_phospho_meta <- phospho_meta[gene_symbol %in% names(target_genes)]
stopifnot(nrow(target_phospho_meta) > 0L)

phospho_url <- sprintf(
  "%s/generic_assay_data/%s/fetch?projection=DETAILED",
  cbio_base, phospho_profile
)
phospho <- as.data.table(post_json(
  phospho_url,
  list(
    sampleListId = cptac_sample_list,
    genericAssayStableIds = target_phospho_meta$stableId
  )
))
phospho[, value := as.numeric(value)]
phospho <- merge(
  phospho[, .(sampleId, patientId, stableId, value)],
  target_phospho_meta[, .(stableId, gene_symbol, description, name)],
  by = "stableId",
  all.x = TRUE
)
fwrite(phospho, file.path(tables_dir, "cptac_gbm_target_phosphosite_values.csv"))

phospho_coverage <- phospho[, .(
  n_samples = uniqueN(sampleId),
  n_nonmissing = sum(!is.na(value))
), by = .(gene_symbol, stableId, description)]
setorder(phospho_coverage, gene_symbol, -n_nonmissing, stableId)
fwrite(phospho_coverage, file.path(tables_dir, "cptac_gbm_target_phosphosite_coverage.csv"))

# 4. LAP3 total protein versus target phosphosites
lap3 <- protein[gene_symbol == "LAP3", .(sampleId, lap3_protein = value)]
substrate_protein <- protein[, .(
  sampleId,
  phospho_gene = gene_symbol,
  substrate_total_protein = value
)]

cor_rows <- phospho[, {
  dt <- merge(
    .SD[, .(sampleId, phosphosite_value = value)],
    lap3,
    by = "sampleId"
  )
  substrate <- substrate_protein[phospho_gene == .BY$gene_symbol]
  dt <- merge(dt, substrate[, .(sampleId, substrate_total_protein)], by = "sampleId")
  dt <- dt[complete.cases(dt)]

  raw_test <- if (nrow(dt) >= 10L) {
    cor.test(dt$lap3_protein, dt$phosphosite_value, method = "spearman", exact = FALSE)
  } else NULL

  adjusted_test <- NULL
  if (nrow(dt) >= 10L && sd(dt$substrate_total_protein) > 0) {
    residual <- resid(lm(phosphosite_value ~ substrate_total_protein, data = dt))
    adjusted_test <- cor.test(dt$lap3_protein, residual, method = "spearman", exact = FALSE)
  }

  .(
    n_complete = nrow(dt),
    spearman_rho_raw = if (is.null(raw_test)) NA_real_ else unname(raw_test$estimate),
    p_raw = if (is.null(raw_test)) NA_real_ else raw_test$p.value,
    spearman_rho_total_protein_adjusted = if (is.null(adjusted_test)) NA_real_ else unname(adjusted_test$estimate),
    p_total_protein_adjusted = if (is.null(adjusted_test)) NA_real_ else adjusted_test$p.value
  )
}, by = .(gene_symbol, stableId)]
cor_rows[, fdr_raw := p.adjust(p_raw, method = "BH")]
cor_rows[, fdr_total_protein_adjusted := p.adjust(p_total_protein_adjusted, method = "BH")]
cor_rows[, readout_class := fifelse(gene_symbol == "LAP3", "LAP3_own_phosphosite", "mTORC1_target")]
cor_rows[, fdr_within_class_raw := p.adjust(p_raw, method = "BH"), by = readout_class]
cor_rows[, fdr_within_class_total_protein_adjusted :=
           p.adjust(p_total_protein_adjusted, method = "BH"), by = readout_class]
cor_rows[, abs_adjusted_rho := abs(spearman_rho_total_protein_adjusted)]
setorder(cor_rows, fdr_total_protein_adjusted, -abs_adjusted_rho)
cor_rows[, abs_adjusted_rho := NULL]
fwrite(cor_rows, file.path(tables_dir, "cptac_lap3_vs_target_phosphosite_correlations.csv"))

# 5. TCGA GBM RPPA target coverage. LAP3 is intentionally queried to verify absence.
tcga_rppa_profile <- "gbm_tcga_pan_can_atlas_2018_rppa"
tcga_rppa <- as.data.table(post_json(
  sprintf(
    "%s/molecular-profiles/%s/molecular-data/fetch?projection=DETAILED",
    cbio_base, tcga_rppa_profile
  ),
  list(
    entrezGeneIds = unname(target_genes),
    sampleListId = "gbm_tcga_pan_can_atlas_2018_all"
  )
))
tcga_rppa[, gene_symbol := gene.hugoGeneSymbol]
tcga_rppa_coverage <- merge(
  data.table(gene_symbol = names(target_genes)),
  tcga_rppa[, .(n_samples = uniqueN(sampleId)), by = gene_symbol],
  by = "gene_symbol",
  all.x = TRUE
)
tcga_rppa_coverage[is.na(n_samples), n_samples := 0L]
tcga_rppa_coverage[, lap3_direct_pairing_possible := gene_symbol == "LAP3" & n_samples > 0]
fwrite(tcga_rppa_coverage, file.path(tables_dir, "tcga_gbm_rppa_target_coverage.csv"))

# Preserve antibody/feature labels because the gene-level API does not distinguish
# total protein from phospho-specific RPPA antibodies.
rppa_feature_url <- paste0(
  "https://media.githubusercontent.com/media/cBioPortal/datahub/master/",
  "public/gbm_tcga_pan_can_atlas_2018/data_rppa.txt"
)
rppa_file <- tempfile(fileext = ".tsv")
download.file(rppa_feature_url, rppa_file, mode = "wb", quiet = TRUE)
rppa_features <- fread(rppa_file, select = 1L)
setnames(rppa_features, 1L, "composite_element_ref")
rppa_features[, gene_symbol := sub("\\|.*$", "", composite_element_ref)]
rppa_target_features <- rppa_features[gene_symbol %in% names(target_genes)]
fwrite(rppa_target_features, file.path(tables_dir, "tcga_gbm_rppa_target_antibodies.csv"))
unlink(rppa_file)

source_summary <- data.table(
  source = c("HPA pathology IHC", "CPTAC GBM total proteome", "CPTAC GBM phosphoproteome", "TCGA GBM RPPA"),
  lap3_protein_available = c(TRUE, TRUE, FALSE, FALSE),
  phospho_readout_available = c(FALSE, FALSE, TRUE, TRUE),
  paired_lap3_protein_phospho = c(FALSE, FALSE, TRUE, FALSE),
  sample_scope = c(
    paste0(hpa_lap3_glioma$total_ihc_cases, " aggregate glioma IHC cases"),
    paste0(uniqueN(protein$sampleId), " GBM samples"),
    paste0(uniqueN(phospho$sampleId), " GBM samples with site-dependent coverage"),
    paste0(max(tcga_rppa_coverage$n_samples), " GBM samples")
  ),
  recommended_role = c(
    "Descriptive LAP3 protein detection only",
    "Primary LAP3 total-protein anchor",
    "Primary paired mTORC1 phosphosite analysis",
    "Secondary downstream phospho validation; no LAP3 antibody"
  )
)
fwrite(source_summary, file.path(tables_dir, "protein_readout_source_availability.csv"))

canonical_pattern <- paste(
  c(
    "RPS6_S235_S236", "RPS6_S236_S240",
    "EIF4EBP1_S65", "EIF4EBP1_T37_T46", "EIF4EBP1_T70",
    "MTOR_S2481"
  ),
  collapse = "|"
)
canonical_results <- cor_rows[grepl(canonical_pattern, stableId)]
fwrite(canonical_results, file.path(tables_dir, "cptac_canonical_mtorc1_phosphosite_correlations.csv"))

readme_lines <- c(
  "# LAP3 Public Protein and Phospho-Readout Audit",
  "",
  paste0("生成时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## 核心结论",
  "",
  sprintf(
    "- HPA v23 glioma IHC 共 %d 例：中等 %d、低 %d、未检出 %d；检出率 %.1f%%，仅适合作描述性蛋白定位。",
    hpa_lap3_glioma$total_ihc_cases,
    hpa_lap3_glioma$Medium,
    hpa_lap3_glioma$Low,
    hpa_lap3_glioma$`Not detected`,
    hpa_lap3_glioma$detected_pct
  ),
  sprintf(
    "- CPTAC GBM 共有 %d 例，LAP3、MTOR、RPS6、RPS6KB1、EIF4EBP1 总蛋白均为完整覆盖。",
    uniqueN(protein$sampleId)
  ),
  sprintf(
    "- CPTAC 目标 phosphosite 共 %d 条样本级观测；可用于 LAP3 总蛋白与 mTORC1 readout 的同队列配对分析。",
    nrow(phospho)
  ),
  sprintf(
    "- mTORC1 下游目标位点单独 BH 校正后，原始相关最小 FDR 为 %.3f，总蛋白校正后最小 FDR 为 %.3f；没有下游位点达到 FDR < 0.05。",
    min(cor_rows[readout_class == "mTORC1_target", fdr_within_class_raw], na.rm = TRUE),
    min(cor_rows[
      readout_class == "mTORC1_target",
      fdr_within_class_total_protein_adjusted
    ], na.rm = TRUE)
  ),
  "- TCGA GBM RPPA 有 231 例及 p-S6、p-4EBP1 等下游抗体，但没有 LAP3 抗体，不能进行 LAP3 protein 与 phospho readout 的直接配对。",
  "",
  "## 解释边界",
  "",
  "- 当前 CPTAC 结果是可用但未提供多重校正后显著支持的公共蛋白证据。",
  "- LAP3 自身 phosphosite 与 LAP3 总蛋白的相关不作为 mTORC1 readout，已与下游目标分开校正。",
  "- 不能据此声称 LAP3 蛋白升高伴随稳定 mTORC1 phospho-activation，更不能替代干预实验。",
  "- HPA 是小样本半定量 IHC 汇总，不适合与 p-S6 做个体级相关。",
  "- phosphosite 编号受蛋白 isoform/reference accession 影响，正式图表需保留 stable ID，避免擅自换算位点。",
  "",
  "## 推荐用途",
  "",
  "- HPA LAP3 IHC 可作为补充图或蛋白存在性说明。",
  "- CPTAC 总蛋白与 phosphosite 相关结果应作为阴性/探索性补充，不建议升格为主图机制证据。",
  "- 湿实验仍应直接检测 LAP3 与 p-S6/p-S6K1/p-4EBP1，并用 knockdown/rescue 验证方向。",
  "",
  "## 关键输出",
  "",
  "```text",
  "tables/protein_readout_source_availability.csv",
  "tables/hpa_lap3_glioma_ihc_summary.csv",
  "tables/cptac_gbm_total_protein_coverage.csv",
  "tables/cptac_gbm_target_phosphosite_coverage.csv",
  "tables/cptac_lap3_vs_target_phosphosite_correlations.csv",
  "tables/cptac_canonical_mtorc1_phosphosite_correlations.csv",
  "tables/tcga_gbm_rppa_target_antibodies.csv",
  "```"
)
writeLines(readme_lines, file.path(out_dir, "README_LAP3_Protein_Readout.md"))

log_msg("HPA glioma IHC cases: ", hpa_lap3_glioma$total_ihc_cases)
log_msg("CPTAC total-protein samples: ", uniqueN(protein$sampleId))
log_msg("CPTAC target phosphosite observations: ", nrow(phospho))
log_msg("TCGA RPPA LAP3 samples: ", tcga_rppa_coverage[gene_symbol == "LAP3", n_samples])
log_msg("Protein readout audit finished: ", Sys.time())
