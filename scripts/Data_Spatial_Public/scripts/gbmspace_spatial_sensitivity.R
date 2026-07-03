#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

project_root <- "/home/lzb/glioma"
assoc_dir <- file.path(project_root, "Data_Spatial_Public/GBM_Space/results/LAP3_Spatial_Association")
cache_dir <- file.path(project_root, "Data_Spatial_Public/GBM_Space/results/Lightweight_Cache")
out_dir <- file.path(project_root, "Data_Spatial_Public/GBM_Space/results/LAP3_Spatial_Sensitivity")
tables_dir <- file.path(out_dir, "tables")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

wilcox_vs_zero <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 4 || all(x == 0)) return(NA_real_)
  suppressWarnings(wilcox.test(x, mu = 0, exact = FALSE)$p.value)
}

summarise_cor <- function(section_cor, readouts, fdr_family_name) {
  tumor_eff <- section_cor[
    readout %in% readouts & is.finite(spearman_rho),
    .(
      n_sections = .N,
      median_section_rho = median(spearman_rho),
      mean_section_rho = mean(spearman_rho)
    ),
    by = .(tumor_id, readout)
  ]
  out <- tumor_eff[
    ,
    .(
      fdr_family = fdr_family_name,
      n_tumors = .N,
      median_tumor_rho = median(median_section_rho),
      mean_tumor_rho = mean(median_section_rho),
      n_positive_tumors = sum(median_section_rho > 0),
      p_wilcox_vs_zero = wilcox_vs_zero(median_section_rho)
    ),
    by = readout
  ]
  out[, fdr := p.adjust(p_wilcox_vs_zero, method = "BH")]
  setorder(out, fdr, readout)
  out[]
}

summarise_contrast <- function(section_contrast, readouts, contrast_name, fdr_family_name) {
  tumor_eff <- section_contrast[
    readout %in% readouts & contrast == contrast_name & is.finite(delta_high_minus_low),
    .(
      n_sections = .N,
      median_section_delta = median(delta_high_minus_low),
      mean_section_delta = mean(delta_high_minus_low)
    ),
    by = .(tumor_id, readout)
  ]
  out <- tumor_eff[
    ,
    .(
      fdr_family = fdr_family_name,
      contrast = contrast_name,
      n_tumors = .N,
      median_tumor_delta = median(median_section_delta),
      mean_tumor_delta = mean(median_section_delta),
      n_positive_tumors = sum(median_section_delta > 0),
      p_wilcox_vs_zero = wilcox_vs_zero(median_section_delta)
    ),
    by = readout
  ]
  out[, fdr := p.adjust(p_wilcox_vs_zero, method = "BH")]
  setorder(out, fdr, readout)
  out[]
}

section_cor <- fread(file.path(assoc_dir, "tables/gbmspace_lap3_section_spearman_all_readouts.tsv"))
section_contrast <- fread(file.path(assoc_dir, "tables/gbmspace_lap3_section_high_low_contrasts_all_readouts.tsv"))
coverage <- fread(file.path(cache_dir, "tables/gbmspace_pathway_gene_coverage.tsv"))
readout_catalog <- fread(file.path(assoc_dir, "tables/gbmspace_spatial_readout_catalog.tsv"))

primary_readouts <- c("HALLMARK_MTORC1_SIGNALING", "LEUCINE_BCAA_CORE", "MTORC1_READOUT_CORE")
translation_readout <- "REACTOME_TRANSLATION"
supportive_readouts <- readout_catalog[
  priority == TRUE | readout_class %in% c("myeloid_tam", "malignant_state", "spatial_niche", "histopath"),
  unique(readout)
]

all_tumors <- sort(unique(section_cor$tumor_id))

primary_base_cor <- summarise_cor(section_cor, primary_readouts, "primary")
primary_base_contrast <- rbindlist(lapply(
  c("top25_vs_bottom25", "detected_vs_undetected"),
  function(x) summarise_contrast(section_contrast, primary_readouts, x, "primary")
))

loto_cor <- rbindlist(lapply(all_tumors, function(drop_tumor) {
  out <- summarise_cor(section_cor[tumor_id != drop_tumor], primary_readouts, "primary")
  out[, dropped_tumor := drop_tumor]
  out
}))
setcolorder(loto_cor, c("dropped_tumor", setdiff(names(loto_cor), "dropped_tumor")))

loto_contrast <- rbindlist(lapply(all_tumors, function(drop_tumor) {
  out <- rbindlist(lapply(
    c("top25_vs_bottom25", "detected_vs_undetected"),
    function(x) summarise_contrast(section_contrast[tumor_id != drop_tumor], primary_readouts, x, "primary")
  ))
  out[, dropped_tumor := drop_tumor]
  out
}))
setcolorder(loto_contrast, c("dropped_tumor", setdiff(names(loto_contrast), "dropped_tumor")))

loto_stability <- loto_cor[
  ,
  .(
    n_loto = .N,
    min_n_tumors = min(n_tumors),
    min_median_tumor_rho = min(median_tumor_rho),
    max_median_tumor_rho = max(median_tumor_rho),
    min_positive_tumors = min(n_positive_tumors),
    max_fdr = max(fdr, na.rm = TRUE),
    all_positive_after_loto = all(n_positive_tumors == n_tumors),
    all_fdr_lt_0_05 = all(fdr < 0.05, na.rm = TRUE)
  ),
  by = readout
]

good_translation_sections <- coverage[signature == translation_readout & present_fraction >= 0.8, h5ad_sample_id]
translation_cor_filtered <- summarise_cor(
  section_cor[readout == translation_readout & h5ad_sample_id %in% good_translation_sections],
  translation_readout,
  "supportive_translation_coverage_filtered"
)
translation_contrast_filtered <- rbindlist(lapply(
  c("top25_vs_bottom25", "detected_vs_undetected"),
  function(x) summarise_contrast(
    section_contrast[readout == translation_readout & h5ad_sample_id %in% good_translation_sections],
    translation_readout,
    x,
    "supportive_translation_coverage_filtered"
  )
))
translation_filter_summary <- data.table(
  readout = translation_readout,
  sections_total = coverage[signature == translation_readout, .N],
  sections_kept = length(good_translation_sections),
  sections_excluded = coverage[signature == translation_readout & present_fraction < 0.8, .N],
  tumors_kept = uniqueN(section_cor[readout == translation_readout & h5ad_sample_id %in% good_translation_sections, tumor_id])
)

histopath_readouts <- readout_catalog[readout_class == "histopath", unique(readout)]
histopath_cor <- summarise_cor(section_cor, histopath_readouts, "supportive_histopath_only")
histopath_contrast <- rbindlist(lapply(
  c("top25_vs_bottom25", "detected_vs_undetected"),
  function(x) summarise_contrast(section_contrast, histopath_readouts, x, "supportive_histopath_only")
))
histopath_scope <- section_cor[
  readout %in% histopath_readouts,
  .(
    n_sections = uniqueN(h5ad_sample_id),
    n_tumors = uniqueN(tumor_id),
    n_readouts = uniqueN(readout)
  )
]

supportive_cor <- summarise_cor(section_cor, supportive_readouts, "secondary_supportive")
supportive_contrast <- rbindlist(lapply(
  c("top25_vs_bottom25", "detected_vs_undetected"),
  function(x) summarise_contrast(section_contrast, supportive_readouts, x, "secondary_supportive")
))

fwrite(primary_base_cor, file.path(tables_dir, "primary_base_correlation_summary.tsv"), sep = "\t")
fwrite(primary_base_contrast, file.path(tables_dir, "primary_base_contrast_summary.tsv"), sep = "\t")
fwrite(loto_cor, file.path(tables_dir, "primary_leave_one_tumor_out_correlations.tsv"), sep = "\t")
fwrite(loto_contrast, file.path(tables_dir, "primary_leave_one_tumor_out_contrasts.tsv"), sep = "\t")
fwrite(loto_stability, file.path(tables_dir, "primary_leave_one_tumor_out_stability.tsv"), sep = "\t")
fwrite(translation_filter_summary, file.path(tables_dir, "translation_coverage_filter_scope.tsv"), sep = "\t")
fwrite(translation_cor_filtered, file.path(tables_dir, "translation_coverage_filtered_correlation.tsv"), sep = "\t")
fwrite(translation_contrast_filtered, file.path(tables_dir, "translation_coverage_filtered_contrast.tsv"), sep = "\t")
fwrite(histopath_scope, file.path(tables_dir, "histopath_only_scope.tsv"), sep = "\t")
fwrite(histopath_cor, file.path(tables_dir, "histopath_only_correlation_summary.tsv"), sep = "\t")
fwrite(histopath_contrast, file.path(tables_dir, "histopath_only_contrast_summary.tsv"), sep = "\t")
fwrite(supportive_cor, file.path(tables_dir, "secondary_supportive_correlation_summary.tsv"), sep = "\t")
fwrite(supportive_contrast, file.path(tables_dir, "secondary_supportive_contrast_summary.tsv"), sep = "\t")

closure_call <- if (
  all(loto_stability$all_positive_after_loto) &&
    all(loto_stability$all_fdr_lt_0_05) &&
    all(primary_base_cor$fdr < 0.05)
) {
  "GO / SUPPORTIVE for spatial co-localization and niche association; NOT causal mechanism."
} else {
  "PARTIAL / SUPPORTIVE, sensitivity caveats remain."
}

cat(
  "# GBM-Space LAP3 Spatial Sensitivity\n\n",
  "生成时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n",
  "## Fixed FDR Families\n\n",
  "- Primary: HALLMARK_MTORC1_SIGNALING, LEUCINE_BCAA_CORE, MTORC1_READOUT_CORE.\n",
  "- Secondary/supportive: REACTOME_TRANSLATION, cell-state, TAM/myeloid, spatial niche, histopath readouts.\n\n",
  "## Sensitivity Summary\n\n",
  "Closure call:\n\n",
  "```text\n", closure_call, "\n```\n\n",
  "Primary base correlations:\n\n",
  "```text\n",
  paste(capture.output(print(primary_base_cor)), collapse = "\n"),
  "\n```\n\n",
  "Leave-one-tumor-out stability:\n\n",
  "```text\n",
  paste(capture.output(print(loto_stability)), collapse = "\n"),
  "\n```\n\n",
  "Translation coverage filter:\n\n",
  "```text\n",
  paste(capture.output(print(translation_filter_summary)), collapse = "\n"),
  "\n",
  paste(capture.output(print(translation_cor_filtered)), collapse = "\n"),
  "\n```\n\n",
  "Histopath-only scope and results:\n\n",
  "```text\n",
  paste(capture.output(print(histopath_scope)), collapse = "\n"),
  "\n",
  paste(capture.output(print(histopath_cor)), collapse = "\n"),
  "\n```\n\n",
  "## Interpretation Boundary\n\n",
  "The sensitivity analyses test robustness of spatial co-localization / niche association.\n",
  "They do not establish LAP3 causality, enzymatic leucine flux, or malignant-cell-intrinsic mTORC1 activation.\n",
  file = file.path(out_dir, "README_LAP3_Spatial_Sensitivity.md")
)

