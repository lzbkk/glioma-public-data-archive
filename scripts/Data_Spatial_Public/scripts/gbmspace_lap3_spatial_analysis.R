#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

project_root <- "/home/lzb/glioma"
cache_dir <- file.path(project_root, "Data_Spatial_Public/GBM_Space/results/Lightweight_Cache")
out_dir <- file.path(project_root, "Data_Spatial_Public/GBM_Space/results/LAP3_Spatial_Association")
tables_dir <- file.path(out_dir, "tables")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

score_file <- file.path(cache_dir, "tables/gbmspace_spot_lap3_pathway_scores.tsv")
meta_file <- file.path(cache_dir, "tables/gbmspace_spot_metadata.tsv")
section_file <- file.path(cache_dir, "tables/gbmspace_section_summary.tsv")
coverage_file <- file.path(cache_dir, "tables/gbmspace_pathway_gene_coverage.tsv")

message_ts <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n")
  flush.console()
}

safe_cor <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 20 || uniqueN(x[keep]) < 3 || uniqueN(y[keep]) < 3) {
    return(list(rho = NA_real_, p = NA_real_, n = sum(keep)))
  }
  ct <- suppressWarnings(cor.test(x[keep], y[keep], method = "spearman", exact = FALSE))
  list(rho = unname(ct$estimate), p = ct$p.value, n = sum(keep))
}

wilcox_greater_zero <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 4 || all(x == 0)) return(NA_real_)
  suppressWarnings(wilcox.test(x, mu = 0, exact = FALSE)$p.value)
}

classify_readout <- function(x) {
  lx <- tolower(x)
  fifelse(
    grepl("HALLMARK_MTORC1|MTORC1_READOUT|LEUCINE_BCAA|REACTOME_TRANSLATION", x),
    "pathway",
    fifelse(
      grepl("tam|monocyte|dendritic|immune", lx),
      "myeloid_tam",
      fifelse(
        grepl("ac|opc|npc|proliferative|hypoxic|dev.like", lx),
        "malignant_state",
        fifelse(
          grepl("niche|gliosis|vasculature|grey.matter|white.matter|pial", lx),
          "spatial_niche",
          fifelse(
            grepl("necrosis|pseudopalisading|cellular.tumor|leading.edge|microvascular|infiltrating|blood.vessels|perinecrotic", lx),
            "histopath",
            "other"
          )
        )
      )
    )
  )
}

priority_readout <- function(x) {
  lx <- tolower(x)
  grepl(
    paste(
      c(
        "HALLMARK_MTORC1", "LEUCINE_BCAA", "MTORC1_READOUT", "REACTOME_TRANSLATION",
        "immune..tams", "resident.tams", "pro.inflammatory.tams", "anti.inflammatory.tams",
        "angiogenic.tams", "proliferative.tams", "monocytes",
        "dev.like..ac", "dev.like..opc", "proliferative", "hypoxic.1",
        "ac.progenitor", "opc.like", "npc.neuronal",
        "cellular.tumor", "necrosis", "pseudopalisading", "microvascular",
        "vasculature", "gliosis"
      ),
      collapse = "|"
    ),
    lx
  )
}

scores <- fread(score_file)
meta <- fread(meta_file, select = c("spot_id", "h5ad_sample_id", "tumor_id", "in_tissue"))
sections <- fread(section_file)
coverage <- fread(coverage_file)

scores <- merge(scores, meta, by = c("spot_id", "h5ad_sample_id", "tumor_id"), all.x = TRUE)
scores <- scores[in_tissue == 1 | is.na(in_tissue)]

pathway_cols <- c(
  "HALLMARK_MTORC1_SIGNALING",
  "LEUCINE_BCAA_CORE",
  "MTORC1_READOUT_CORE",
  "REACTOME_TRANSLATION"
)

cor_rows <- list()
contrast_rows <- list()
readout_catalog_rows <- list()

non_gene_files <- sections$non_gene_feature_file
names(non_gene_files) <- sections$h5ad_sample_id

for (i in seq_along(non_gene_files)) {
  sample_id <- names(non_gene_files)[[i]]
  file <- non_gene_files[[i]]
  message_ts(sprintf("[%d/%d] %s", i, length(non_gene_files), sample_id))

  sec_score <- scores[h5ad_sample_id == sample_id]
  if (nrow(sec_score) < 20 || !file.exists(file)) next

  ng <- fread(file)
  sec <- merge(
    sec_score,
    ng,
    by = c("spot_id", "h5ad_sample_id"),
    all.x = TRUE
  )

  readout_cols <- setdiff(names(sec), c(
    "spot_id", "h5ad_sample_id", "sample_name", "tumor_id",
    "LAP3_raw", "LAP3_detected", "LAP3_log1p_cp10k",
    "gene_library_size", "detected_gene_features", "in_tissue"
  ))
  readout_cols <- readout_cols[sapply(sec[, ..readout_cols], is.numeric)]

  q25 <- as.numeric(quantile(sec$LAP3_log1p_cp10k, 0.25, na.rm = TRUE))
  q75 <- as.numeric(quantile(sec$LAP3_log1p_cp10k, 0.75, na.rm = TRUE))
  high <- sec$LAP3_log1p_cp10k >= q75
  low <- sec$LAP3_log1p_cp10k <= q25

  for (readout in readout_cols) {
    cls <- classify_readout(readout)
    readout_catalog_rows[[length(readout_catalog_rows) + 1L]] <- data.table(
      readout = readout,
      readout_class = cls
    )

    cv <- safe_cor(sec$LAP3_log1p_cp10k, sec[[readout]])
    cor_rows[[length(cor_rows) + 1L]] <- data.table(
      h5ad_sample_id = sample_id,
      tumor_id = sec$tumor_id[[1]],
      readout = readout,
      readout_class = cls,
      n_spots = cv$n,
      spearman_rho = cv$rho,
      spearman_p = cv$p
    )

    contrast_rows[[length(contrast_rows) + 1L]] <- data.table(
      h5ad_sample_id = sample_id,
      tumor_id = sec$tumor_id[[1]],
      readout = readout,
      readout_class = cls,
      contrast = "top25_vs_bottom25",
      n_high = sum(high, na.rm = TRUE),
      n_low = sum(low, na.rm = TRUE),
      mean_high = mean(sec[[readout]][high], na.rm = TRUE),
      mean_low = mean(sec[[readout]][low], na.rm = TRUE),
      delta_high_minus_low = mean(sec[[readout]][high], na.rm = TRUE) - mean(sec[[readout]][low], na.rm = TRUE)
    )

    det <- sec$LAP3_detected == TRUE
    undet <- sec$LAP3_detected == FALSE
    if (sum(det, na.rm = TRUE) >= 20 && sum(undet, na.rm = TRUE) >= 20) {
      contrast_rows[[length(contrast_rows) + 1L]] <- data.table(
        h5ad_sample_id = sample_id,
        tumor_id = sec$tumor_id[[1]],
        readout = readout,
        readout_class = cls,
        contrast = "detected_vs_undetected",
        n_high = sum(det, na.rm = TRUE),
        n_low = sum(undet, na.rm = TRUE),
        mean_high = mean(sec[[readout]][det], na.rm = TRUE),
        mean_low = mean(sec[[readout]][undet], na.rm = TRUE),
        delta_high_minus_low = mean(sec[[readout]][det], na.rm = TRUE) - mean(sec[[readout]][undet], na.rm = TRUE)
      )
    }
  }
}

section_cor <- rbindlist(cor_rows, use.names = TRUE, fill = TRUE)
section_contrast <- rbindlist(contrast_rows, use.names = TRUE, fill = TRUE)
readout_catalog <- unique(rbindlist(readout_catalog_rows, use.names = TRUE, fill = TRUE))
readout_catalog[, priority := priority_readout(readout)]

section_cor <- merge(section_cor, readout_catalog, by = c("readout", "readout_class"), all.x = TRUE)
section_contrast <- merge(section_contrast, readout_catalog, by = c("readout", "readout_class"), all.x = TRUE)

tumor_cor <- section_cor[
  is.finite(spearman_rho),
  .(
    n_sections = .N,
    median_section_rho = median(spearman_rho),
    mean_section_rho = mean(spearman_rho)
  ),
  by = .(tumor_id, readout, readout_class, priority)
]

tumor_contrast <- section_contrast[
  is.finite(delta_high_minus_low),
  .(
    n_sections = .N,
    median_section_delta = median(delta_high_minus_low),
    mean_section_delta = mean(delta_high_minus_low)
  ),
  by = .(tumor_id, readout, readout_class, priority, contrast)
]

cor_summary <- tumor_cor[
  ,
  .(
    n_tumors = .N,
    median_tumor_rho = median(median_section_rho),
    mean_tumor_rho = mean(median_section_rho),
    n_positive_tumors = sum(median_section_rho > 0),
    p_wilcox_vs_zero = wilcox_greater_zero(median_section_rho)
  ),
  by = .(readout, readout_class, priority)
]
cor_summary[, fdr := p.adjust(p_wilcox_vs_zero, method = "BH")]
cor_summary[, abs_median_tumor_rho := abs(median_tumor_rho)]
setorder(cor_summary, fdr, -abs_median_tumor_rho)

contrast_summary <- tumor_contrast[
  ,
  .(
    n_tumors = .N,
    median_tumor_delta = median(median_section_delta),
    mean_tumor_delta = mean(median_section_delta),
    n_positive_tumors = sum(median_section_delta > 0),
    p_wilcox_vs_zero = wilcox_greater_zero(median_section_delta)
  ),
  by = .(readout, readout_class, priority, contrast)
]
contrast_summary[, fdr := p.adjust(p_wilcox_vs_zero, method = "BH"), by = contrast]
contrast_summary[, abs_median_tumor_delta := abs(median_tumor_delta)]
setorder(contrast_summary, contrast, fdr, -abs_median_tumor_delta)

pathway_coverage_summary <- coverage[
  ,
  .(
    n_sections = .N,
    min_present_fraction = min(present_fraction),
    median_present_fraction = median(present_fraction),
    n_sections_below_0_8 = sum(present_fraction < 0.8)
  ),
  by = signature
]

priority_cor <- cor_summary[priority == TRUE | readout_class == "pathway"]
priority_contrast <- contrast_summary[priority == TRUE | readout_class == "pathway"]

fwrite(section_cor, file.path(tables_dir, "gbmspace_lap3_section_spearman_all_readouts.tsv"), sep = "\t")
fwrite(section_contrast, file.path(tables_dir, "gbmspace_lap3_section_high_low_contrasts_all_readouts.tsv"), sep = "\t")
fwrite(tumor_cor, file.path(tables_dir, "gbmspace_lap3_tumor_median_correlations.tsv"), sep = "\t")
fwrite(tumor_contrast, file.path(tables_dir, "gbmspace_lap3_tumor_median_high_low_contrasts.tsv"), sep = "\t")
fwrite(cor_summary, file.path(tables_dir, "gbmspace_lap3_tumor_level_correlation_summary.tsv"), sep = "\t")
fwrite(contrast_summary, file.path(tables_dir, "gbmspace_lap3_tumor_level_contrast_summary.tsv"), sep = "\t")
fwrite(priority_cor, file.path(tables_dir, "gbmspace_lap3_priority_correlation_summary.tsv"), sep = "\t")
fwrite(priority_contrast, file.path(tables_dir, "gbmspace_lap3_priority_contrast_summary.tsv"), sep = "\t")
fwrite(pathway_coverage_summary, file.path(tables_dir, "gbmspace_pathway_coverage_analysis_summary.tsv"), sep = "\t")
fwrite(readout_catalog, file.path(tables_dir, "gbmspace_spatial_readout_catalog.tsv"), sep = "\t")

top_cor <- head(priority_cor[is.finite(fdr)], 20)
top_delta <- head(priority_contrast[contrast == "top25_vs_bottom25" & is.finite(fdr)], 20)

cat(
  "# GBM-Space LAP3 Spatial Association\n\n",
  "生成时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n",
  "## 输入\n\n",
  "- Lightweight cache: `", cache_dir, "`\n",
  "- sections: `", uniqueN(scores$h5ad_sample_id), "`\n",
  "- tumors: `", uniqueN(scores$tumor_id), "`\n",
  "- in-tissue spots used: `", nrow(scores), "`\n\n",
  "## 方法概述\n\n",
  "每个 section 内计算 LAP3 log1p(CP10K) 与 pathway、cell-state、spatial-niche、histopath readout 的 Spearman 相关；\n",
  "同时计算 LAP3 top quartile vs bottom quartile、LAP3 detected vs undetected 的 section 内 readout 差异。\n",
  "所有统计检验均先在 section 内得到效应量，再按 tumor 汇总，以 tumor-level median effect 作为推断单位。\n",
  "Wilcoxon signed-rank test 只作用在 12 个 tumor-level effects 上；spot-level p 值不作为主推断。\n\n",
  "## 关键输出\n\n",
  "- `tables/gbmspace_lap3_section_spearman_all_readouts.tsv`\n",
  "- `tables/gbmspace_lap3_section_high_low_contrasts_all_readouts.tsv`\n",
  "- `tables/gbmspace_lap3_tumor_level_correlation_summary.tsv`\n",
  "- `tables/gbmspace_lap3_tumor_level_contrast_summary.tsv`\n",
  "- `tables/gbmspace_lap3_priority_correlation_summary.tsv`\n",
  "- `tables/gbmspace_lap3_priority_contrast_summary.tsv`\n",
  "- `tables/gbmspace_pathway_coverage_analysis_summary.tsv`\n\n",
  "## 第一版结果摘要\n\n",
  "- readouts tested: `", uniqueN(readout_catalog$readout), "`\n",
  "- priority readouts: `", sum(readout_catalog$priority | readout_catalog$readout_class == "pathway"), "`\n",
  "- pathway coverage caveat: REACTOME_TRANSLATION has `",
  pathway_coverage_summary[signature == "REACTOME_TRANSLATION", n_sections_below_0_8],
  "` sections below 0.8 gene coverage.\n\n",
  "Top priority correlations by tumor-level FDR:\n\n",
  "```text\n",
  paste(capture.output(print(top_cor[, .(readout, readout_class, n_tumors, median_tumor_rho, n_positive_tumors, fdr)])), collapse = "\n"),
  "\n```\n\n",
  "Top priority top25-vs-bottom25 contrasts by tumor-level FDR:\n\n",
  "```text\n",
  paste(capture.output(print(top_delta[, .(readout, readout_class, n_tumors, median_tumor_delta, n_positive_tumors, fdr)])), collapse = "\n"),
  "\n```\n\n",
  "## 解释边界\n\n",
  "本分析支持或反驳的是 LAP3-high spatial readout co-localization / neighborhood association。\n",
  "Visium spots 是多细胞混合单位；即便结果显著，也不能单独证明 LAP3 的 malignant-cell-intrinsic 因果机制。\n",
  "正式论文中应使用 association、co-localization、spatial niche 等措辞，并保留 wet-lab perturbation requirement。\n",
  file = file.path(out_dir, "README_LAP3_Spatial_Association.md")
)

message_ts("DONE")
