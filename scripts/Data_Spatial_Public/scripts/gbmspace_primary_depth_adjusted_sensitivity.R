#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

project_root <- "/home/lzb/glioma"
cache_dir <- file.path(project_root, "Data_Spatial_Public/GBM_Space/results/Lightweight_Cache")
out_dir <- file.path(project_root, "Data_Spatial_Public/GBM_Space/results/LAP3_Spatial_Depth_Adjusted_Sensitivity")
tables_dir <- file.path(out_dir, "tables")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

wilcox_vs_zero <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 4 || all(x == 0)) return(NA_real_)
  suppressWarnings(wilcox.test(x, mu = 0, exact = FALSE)$p.value)
}

residual_spearman <- function(sec, ycol) {
  keep <- is.finite(sec$LAP3_log1p_cp10k) &
    is.finite(sec[[ycol]]) &
    is.finite(sec$gene_library_size) &
    is.finite(sec$detected_gene_features)

  if (sum(keep) < 20 ||
      uniqueN(sec$LAP3_log1p_cp10k[keep]) < 3 ||
      uniqueN(sec[[ycol]][keep]) < 3) {
    return(NA_real_)
  }

  sub <- sec[keep]
  rx <- rank(sub$LAP3_log1p_cp10k, ties.method = "average")
  ry <- rank(sub[[ycol]], ties.method = "average")
  rlib <- rank(log1p(sub$gene_library_size), ties.method = "average")
  rfeat <- rank(sub$detected_gene_features, ties.method = "average")

  xres <- residuals(lm(rx ~ rlib + rfeat))
  yres <- residuals(lm(ry ~ rlib + rfeat))
  suppressWarnings(cor(xres, yres, method = "pearson"))
}

raw_spearman <- function(sec, ycol) {
  suppressWarnings(cor(
    sec$LAP3_log1p_cp10k,
    sec[[ycol]],
    method = "spearman",
    use = "complete.obs"
  ))
}

score_file <- file.path(cache_dir, "tables/gbmspace_spot_lap3_pathway_scores.tsv")
meta_file <- file.path(cache_dir, "tables/gbmspace_spot_metadata.tsv")

primary_readouts <- c(
  "HALLMARK_MTORC1_SIGNALING",
  "LEUCINE_BCAA_CORE",
  "MTORC1_READOUT_CORE"
)

scores <- fread(
  score_file,
  select = c(
    "spot_id", "h5ad_sample_id", "tumor_id",
    "LAP3_log1p_cp10k", "gene_library_size", "detected_gene_features",
    primary_readouts
  )
)
meta <- fread(meta_file, select = c("spot_id", "h5ad_sample_id", "tumor_id", "in_tissue"))
dt <- merge(scores, meta, by = c("spot_id", "h5ad_sample_id", "tumor_id"), all.x = TRUE)
dt <- dt[in_tissue == 1 | is.na(in_tissue)]

section_rows <- list()
for (sample_id in unique(dt$h5ad_sample_id)) {
  sec <- dt[h5ad_sample_id == sample_id]
  for (readout in primary_readouts) {
    section_rows[[length(section_rows) + 1L]] <- data.table(
      h5ad_sample_id = sample_id,
      tumor_id = sec$tumor_id[[1]],
      readout = readout,
      n_spots = nrow(sec),
      raw_rho = raw_spearman(sec, readout),
      depth_adjusted_rho = residual_spearman(sec, readout),
      lap3_library_rho = suppressWarnings(cor(
        sec$LAP3_log1p_cp10k,
        log1p(sec$gene_library_size),
        method = "spearman",
        use = "complete.obs"
      )),
      lap3_features_rho = suppressWarnings(cor(
        sec$LAP3_log1p_cp10k,
        sec$detected_gene_features,
        method = "spearman",
        use = "complete.obs"
      ))
    )
  }
}

section_effects <- rbindlist(section_rows, use.names = TRUE)

tumor_effects <- section_effects[
  ,
  .(
    n_sections = .N,
    median_raw_rho = median(raw_rho, na.rm = TRUE),
    median_depth_adjusted_rho = median(depth_adjusted_rho, na.rm = TRUE),
    median_lap3_library_rho = median(lap3_library_rho, na.rm = TRUE),
    median_lap3_features_rho = median(lap3_features_rho, na.rm = TRUE)
  ),
  by = .(tumor_id, readout)
]

summary <- tumor_effects[
  ,
  .(
    n_tumors = .N,
    median_raw_rho = median(median_raw_rho, na.rm = TRUE),
    median_depth_adjusted_rho = median(median_depth_adjusted_rho, na.rm = TRUE),
    n_positive_depth_adjusted = sum(median_depth_adjusted_rho > 0, na.rm = TRUE),
    p_depth_adjusted = wilcox_vs_zero(median_depth_adjusted_rho),
    median_lap3_library_rho = median(median_lap3_library_rho, na.rm = TRUE),
    median_lap3_features_rho = median(median_lap3_features_rho, na.rm = TRUE)
  ),
  by = readout
]
summary[, fdr_depth_adjusted := p.adjust(p_depth_adjusted, method = "BH")]

fwrite(section_effects, file.path(tables_dir, "primary_depth_adjusted_section_effects.tsv"), sep = "\t")
fwrite(tumor_effects, file.path(tables_dir, "primary_depth_adjusted_tumor_effects.tsv"), sep = "\t")
fwrite(summary, file.path(tables_dir, "primary_depth_adjusted_summary.tsv"), sep = "\t")

cat(
  "# GBM-Space Primary Pathway Depth-Adjusted Sensitivity\n\n",
  "生成时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n",
  "## 输入\n\n",
  "- Lightweight cache: `", cache_dir, "`\n",
  "- Primary readouts: `", paste(primary_readouts, collapse = ", "), "`\n",
  "- Spots after in-tissue filter: `", nrow(dt), "`\n",
  "- Sections: `", uniqueN(dt$h5ad_sample_id), "`\n",
  "- Tumors: `", uniqueN(dt$tumor_id), "`\n\n",
  "## 方法\n\n",
  "每个 section 内对 LAP3 和 primary pathway score 分别取 rank，\n",
  "再用 rank(log1p(gene_library_size)) 与 rank(detected_gene_features) 回归残差化，\n",
  "最后计算两个残差的 Pearson correlation，作为 depth/complexity-adjusted Spearman 近似。\n",
  "section effect 再按 tumor 取 median，并在 12 个 tumor-level effects 上做 Wilcoxon signed-rank test。\n\n",
  "## 结果\n\n",
  "```text\n",
  paste(capture.output(print(summary)), collapse = "\n"),
  "\n```\n\n",
  "## 解释边界\n\n",
  "该敏感性专门评估 primary pathway spatial co-expression 是否受 spot-level sequencing depth / gene complexity 影响。\n",
  "若 depth-adjusted rho 接近 0 且 FDR 不显著，GBM-Space pathway 结果应降级为未校正表达共变/探索性证据，\n",
  "不能作为强机制或强 spatial pathway closure 证据。\n",
  file = file.path(out_dir, "README_Primary_Depth_Adjusted_Sensitivity.md")
)
