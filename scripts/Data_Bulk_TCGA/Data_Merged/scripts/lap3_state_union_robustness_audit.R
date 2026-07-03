#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})
data.table::setDTthreads(8)

root <- "/home/lzb/glioma"
out_dir <- file.path(root, "Data_Bulk_TCGA/Data_Merged/results/LAP3_State_Union_Robustness_Audit")
table_dir <- file.path(out_dir, "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

read_dt <- function(path, sep = ",") {
  fread(file.path(root, path), sep = sep)
}

submodule_genes <- read_dt("Data_Bulk_TCGA/Data_Merged/results/LAP3_State_Submodules/tables/lap3_state_submodule_gene_assignment.csv")
bulk_cor <- read_dt("Data_Bulk_TCGA/Data_Merged/results/LAP3_State_Submodules/tables/lap3_state_submodule_bulk_correlations.csv")
gbmap_cor <- read_dt("Data_Bulk_TCGA/Data_Merged/results/LAP3_State_Submodules/tables/gbmap_core_lap3_state_submodule_lap3_associations.csv")
cptac_cor <- read_dt("Data_Bulk_TCGA/Data_Merged/results/LAP3_State_CPTAC_GLASS_Projection/tables/cptac_lap3_state_submodule_correlations.csv")
glass_tests <- read_dt("Data_Bulk_TCGA/Data_Merged/results/LAP3_State_CPTAC_GLASS_Projection/tables/glass_lap3_state_submodule_paired_tests.csv")
gene_set_counts <- read_dt("Data_Bulk_TCGA/Data_Merged/results/LAP3_State_Module/tables/lap3_state_gene_set_counts.csv")
spatial_contrast <- read_dt("Data_Spatial_Public/GBM_Space/results/LAP3_State_Spatial_Topology/tables/gbmspace_lap3_state_contrast_summary.tsv", sep = "\t")
spatial_loto <- read_dt("Data_Spatial_Public/GBM_Space/results/LAP3_State_Spatial_Topology/tables/gbmspace_lap3_state_leave_one_tumor_out.tsv", sep = "\t")

module_counts <- submodule_genes[, .(
  n_genes = uniqueN(gene),
  n_translation_flag = sum(flag_translation, na.rm = TRUE),
  n_proteostasis_flag = sum(flag_proteostasis, na.rm = TRUE),
  n_myeloid_tam_flag = sum(flag_myeloid_tam, na.rm = TRUE),
  n_hypoxia_flag = sum(flag_hypoxia_perinecrotic, na.rm = TRUE),
  n_malignant_flag = sum(flag_malignant_state, na.rm = TRUE),
  n_tcga_top150 = sum(in_tcga_top150, na.rm = TRUE),
  n_gbmap_up = sum(in_gbmap_up, na.rm = TRUE)
), by = .(submodule = primary_submodule)]

bulk_lap3 <- bulk_cor[
  group == "all" & variable == "LAP3_log2_expr",
  .(dataset, submodule, bulk_n = n, bulk_rho = spearman_rho, bulk_fdr = p_adj_BH)
]
bulk_wide <- dcast(
  bulk_lap3,
  submodule ~ dataset,
  value.var = c("bulk_n", "bulk_rho", "bulk_fdr")
)

gbmap_main <- gbmap_cor[
  entry_variant == "main_neoplastic_exclude_neftel2019" & threshold == 20
]
gbmap_summary <- gbmap_main[, .(
  gbmap_states_tested = .N,
  gbmap_supported_states = sum(p_adj_BH < 0.05 & spearman_rho > 0, na.rm = TRUE),
  gbmap_best_state = author_state[which.min(p_adj_BH)],
  gbmap_best_rho = spearman_rho[which.min(p_adj_BH)],
  gbmap_best_fdr = min(p_adj_BH, na.rm = TRUE),
  gbmap_leave_one_author_min_rho = min(leave_one_author_min_rho, na.rm = TRUE),
  gbmap_leave_one_author_max_rho = max(leave_one_author_max_rho, na.rm = TRUE)
), by = submodule]

cptac_summary <- cptac_cor[
  stratum == "All" & grepl("^mrna_LAP3_", exposure),
  .(
    module = sub("^mrna_", "", exposure),
    outcome,
    cptac_n = n_complete,
    cptac_rho = spearman_rho,
    cptac_fdr = fdr
  )
]
cptac_wide <- dcast(
  cptac_summary[outcome %chin% c("lap3_protein", "phospho_mtorc1_target_score", "bcaa_composite")],
  module ~ outcome,
  value.var = c("cptac_n", "cptac_rho", "cptac_fdr")
)
setnames(cptac_wide, "module", "submodule")

glass_summary <- glass_tests[
  stratum == "All" &
    analysis %chin% c("paired_recurrence_change", "delta_lap3_delta_module"),
  .(
    submodule = module,
    analysis,
    glass_n_pairs = n_pairs,
    glass_median_delta = median_delta,
    glass_estimate = estimate,
    glass_fdr = fdr
  )
]
glass_wide <- dcast(
  glass_summary,
  submodule ~ analysis,
  value.var = c("glass_n_pairs", "glass_median_delta", "glass_estimate", "glass_fdr")
)

component_matrix <- Reduce(
  function(x, y) merge(x, y, by = "submodule", all.x = TRUE),
  list(module_counts, bulk_wide, gbmap_summary, cptac_wide, glass_wide)
)

component_matrix[, support_score := 0L]
if (all(c("bulk_fdr_TCGA", "bulk_rho_TCGA") %in% names(component_matrix))) {
  component_matrix[bulk_fdr_TCGA < 0.05 & bulk_rho_TCGA > 0, support_score := support_score + 1L]
}
cgga_datasets <- sub("^bulk_fdr_", "", grep("^bulk_fdr_CGGA", names(component_matrix), value = TRUE))
if (length(cgga_datasets) > 0L) {
  cgga_supported <- rep(FALSE, nrow(component_matrix))
  for (dataset in cgga_datasets) {
    fdr_col <- paste0("bulk_fdr_", dataset)
    rho_col <- paste0("bulk_rho_", dataset)
    if (all(c(fdr_col, rho_col) %in% names(component_matrix))) {
      cgga_supported <- cgga_supported |
        (component_matrix[[fdr_col]] < 0.05 & component_matrix[[rho_col]] > 0)
    }
  }
  component_matrix[cgga_supported == TRUE, support_score := support_score + 1L]
}
component_matrix[gbmap_supported_states >= 1, support_score := support_score + 1L]
component_matrix[cptac_fdr_lap3_protein < 0.05 & cptac_rho_lap3_protein > 0, support_score := support_score + 1L]
component_matrix[cptac_fdr_phospho_mtorc1_target_score < 0.05 & cptac_rho_phospho_mtorc1_target_score > 0, support_score := support_score + 1L]
component_matrix[glass_fdr_delta_lap3_delta_module < 0.05 & glass_estimate_delta_lap3_delta_module > 0, support_score := support_score + 1L]
component_matrix[, robustness_class := fifelse(
  support_score >= 4, "strong_multimodal_component",
  fifelse(support_score >= 2, "supportive_component", "weak_or_context_specific_component")
)]
setorder(component_matrix, -support_score, submodule)
fwrite(component_matrix, file.path(table_dir, "lap3_state_union_component_evidence_matrix.csv"))

union_rows <- rbindlist(list(
  gene_set_counts[, .(
    state_set,
    metric = "gene_set_composition",
    n = n_genes,
    detail = paste0(
      "tcga_top150=", n_tcga_top150,
      ";gbmap_up=", n_gbmap_up,
      ";translation_proteostasis=", n_translation_proteostasis
    )
  )],
  bulk_cor[
    group == "all" & submodule %chin% c("LAP3_STATE_UNION", "LAP3_STATE_UNION_NO_TRANSLATION_PROTEOSTASIS") &
      variable %chin% c("LAP3_log2_expr", "HALLMARK_MTORC1_SIGNALING"),
    .(
      state_set = submodule,
      metric = paste(dataset, variable, sep = "::"),
      n,
      detail = paste0("rho=", signif(spearman_rho, 4), ";fdr=", signif(p_adj_BH, 4))
    )
  ]
), use.names = TRUE, fill = TRUE)
fwrite(union_rows, file.path(table_dir, "lap3_state_union_variant_support.csv"))

spatial_loto_summary <- spatial_loto[
  state_set == "LAP3_STATE_UNION" & n_tumors >= 8,
  .(
    n_leave_one_runs = .N,
    min_n_tumors = min(n_tumors, na.rm = TRUE),
    min_median_depth_adjusted_rho = min(median_depth_adjusted_rho, na.rm = TRUE),
    max_median_depth_adjusted_rho = max(median_depth_adjusted_rho, na.rm = TRUE),
    all_positive_after_leave_one = all(n_positive_depth_adjusted > 0, na.rm = TRUE),
    max_p_depth_adjusted = max(p_depth_adjusted, na.rm = TRUE)
  ),
  by = .(analysis_type, readout)
]
setorder(spatial_loto_summary, analysis_type, -min_median_depth_adjusted_rho)
fwrite(spatial_loto_summary, file.path(table_dir, "lap3_state_union_spatial_leave_one_tumor_summary.csv"))

spatial_priority <- spatial_contrast[
  state_set == "LAP3_STATE_UNION" & priority == TRUE,
  .(
    readout_class,
    readout,
    n_tumors,
    median_tumor_delta,
    n_positive_delta,
    p_delta,
    fdr_delta_priority
  )
]
setorder(spatial_priority, readout_class, -median_tumor_delta)
fwrite(spatial_priority, file.path(table_dir, "lap3_state_union_spatial_priority_readouts.csv"))

verdict <- data.table(
  question = c(
    "Is LAP3_STATE_UNION dominated by translation/proteostasis genes?",
    "Does the non-translation/proteostasis variant remain usable?",
    "Are multiple biological components supported?",
    "Is the spatial signal driven by one tumor?",
    "Can this audit prove LAP3 causal mechanism?"
  ),
  answer = c(
    "No. Translation/proteostasis contributes to the union, but the frozen no-translation/proteostasis variant retains 170 genes and is tracked separately.",
    "Yes, as a robustness variant. It should be used to support the state/ecology interpretation, not the mTOR/BCAA causal claim.",
    paste0(sum(component_matrix$support_score >= 2, na.rm = TRUE), " submodules have at least two independent support layers; ",
           sum(component_matrix$support_score >= 4, na.rm = TRUE), " reach strong multimodal support."),
    if (nrow(spatial_loto_summary) > 0 && all(spatial_loto_summary$all_positive_after_leave_one, na.rm = TRUE)) {
      "No obvious one-tumor domination among retained priority spatial readouts."
    } else {
      "Some spatial readouts require cautious wording; inspect leave-one-tumor table."
    },
    "No. It strengthens robustness and boundary control, but mechanism still requires wet-lab perturbation."
  ),
  manuscript_use = c(
    "Methods/Results robustness sentence and Supplementary table.",
    "Supplementary sensitivity and Results caveat.",
    "Figure 2/Discussion support for a composite state program.",
    "Figure 4/spatial closure support.",
    "Discussion limitation and wet-lab rationale."
  )
)
fwrite(verdict, file.path(table_dir, "lap3_state_union_robustness_verdict.csv"))

readme <- c(
  "# LAP3_STATE_UNION Robustness Audit",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Purpose",
  "",
  "This audit checks whether the frozen LAP3_STATE_UNION is a single-component artifact, especially whether it is overly driven by translation/proteostasis genes. It uses existing validated result tables and does not rerun matrix-level scoring.",
  "",
  "## Inputs",
  "",
  "- LAP3 state frozen gene-set counts.",
  "- LAP3 submodule gene assignment and bulk correlations.",
  "- GBmap Core patient-state submodule associations.",
  "- CPTAC/GLASS submodule projection results.",
  "- GBM-Space LAP3_STATE_UNION spatial topology and leave-one-tumor-out summaries.",
  "",
  "## Outputs",
  "",
  "- `tables/lap3_state_union_component_evidence_matrix.csv`: component-level multimodal evidence matrix.",
  "- `tables/lap3_state_union_variant_support.csv`: full union vs no-translation/proteostasis support.",
  "- `tables/lap3_state_union_spatial_leave_one_tumor_summary.csv`: spatial robustness against one-tumor domination.",
  "- `tables/lap3_state_union_spatial_priority_readouts.csv`: priority spatial readout summary.",
  "- `tables/lap3_state_union_robustness_verdict.csv`: manuscript-facing verdicts.",
  "",
  "## Interpretation Boundary",
  "",
  "This audit supports robustness of a composite LAP3-state/ecology program. It does not rescue the rejected LAP3-BCAA-mTORC1 causal mechanism and should not be described as causal proof."
)
writeLines(readme, file.path(out_dir, "README.md"))

message("Robustness audit written: ", out_dir)
