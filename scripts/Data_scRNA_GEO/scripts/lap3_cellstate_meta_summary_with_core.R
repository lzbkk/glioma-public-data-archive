#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

setwd("/home/lzb/glioma")
data.table::setDTthreads(8)

out_dir <- "Data_scRNA_GEO/results/LAP3_CellState_MetaSummary_v2_Core"
tables_dir <- file.path(out_dir, "tables")
source_dir <- file.path(out_dir, "source_data")
logs_dir <- file.path(out_dir, "logs")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(logs_dir, "lap3_cellstate_meta_summary_with_core.log")
cat("Started: ", format(Sys.time()), "\n", file = log_file)

read_csv <- function(path) {
  if (!file.exists(path)) stop("Missing input: ", path)
  fread(path)
}

classify_signal <- function(rho, p_adj, n, fdr_cut = 0.05, trend_cut = 0.15) {
  out <- rep("not_evaluable", length(rho))
  ok <- !is.na(rho) & !is.na(n) & n >= 6L
  out[ok] <- "not_reproduced"
  out[ok & !is.na(p_adj) & p_adj < fdr_cut] <- "supported"
  out[ok & !is.na(p_adj) & p_adj >= fdr_cut & p_adj < trend_cut] <- "trend"
  out
}

manifest <- data.table(
  role = c(
    "old_within_state_pathway",
    "old_state_continuum",
    "old_primary_recurrence",
    "core_within_state_pathway",
    "core_author_adjusted_pathway",
    "core_state_preference",
    "core_continuous_state"
  ),
  path = c(
    "Data_scRNA_GEO/results/LAP3_CellState_MetaSummary/tables/within_state_pathway_summary.csv",
    "Data_scRNA_GEO/results/LAP3_CellState_MetaSummary/tables/state_continuum_summary.csv",
    "Data_scRNA_GEO/results/LAP3_CellState_MetaSummary/tables/primary_pathway_recurrence_summary.csv",
    "Data_scRNA_GEO/GBmap_Core/results/LAP3_CellState_Strict/tables/gbmap_core_within_state_pathway_associations.csv",
    "Data_scRNA_GEO/GBmap_Core/results/LAP3_CellState_Strict/tables/gbmap_core_author_adjusted_within_state_pathway_associations.csv",
    "Data_scRNA_GEO/GBmap_Core/results/LAP3_CellState_Strict/tables/gbmap_core_lap3_state_preference_fixed_effect.csv",
    "Data_scRNA_GEO/GBmap_Core/results/LAP3_CellState_Strict/tables/gbmap_core_lap3_continuous_state_associations.csv"
  )
)
manifest[, exists := file.exists(path)]
if (any(!manifest$exists)) {
  print(manifest[exists == FALSE])
  stop("Missing inputs")
}
fwrite(manifest, file.path(source_dir, "metasummary_v2_input_manifest.csv"))

old_pathway <- read_csv(manifest[role == "old_within_state_pathway", path])
old_state <- read_csv(manifest[role == "old_state_continuum", path])
old_primary <- read_csv(manifest[role == "old_primary_recurrence", path])
core_raw <- read_csv(manifest[role == "core_within_state_pathway", path])
core_adj <- read_csv(manifest[role == "core_author_adjusted_pathway", path])
core_pref <- read_csv(manifest[role == "core_state_preference", path])
core_cont <- read_csv(manifest[role == "core_continuous_state", path])

core_adj_primary <- core_adj[
  fdr_family == "primary",
  .(
    dataset = "Core_GBmap_minus_Neftel2019",
    analysis_class = "IDHwt_GBM_harmonized_atlas",
    evidence_role = "atlas_level_supportive",
    entry_variant,
    threshold,
    state = author_state,
    pathway,
    n_donor_states = n,
    n_authors,
    raw_rho,
    raw_p_adj,
    author_adjusted_rho,
    author_adjusted_p_adj,
    author_adjusted_level = classify_signal(author_adjusted_rho, author_adjusted_p_adj, n)
  )
][order(entry_variant, threshold, state, pathway)]
fwrite(core_adj_primary, file.path(tables_dir, "core_author_adjusted_primary_summary.csv"))

old_primary_detail <- old_pathway[
  primary_readout == TRUE,
  .(
    dataset,
    analysis_class,
    evidence_role = "independent_dataset",
    state,
    pathway,
    n_patients,
    rho = spearman_rho,
    p_adj,
    evidence_level
  )
]

core_for_compare <- core_adj_primary[
  entry_variant == "main_neoplastic_exclude_neftel2019" & threshold == 20,
  .(
    dataset,
    analysis_class,
    evidence_role,
    state,
    pathway,
    n_patients = n_donor_states,
    rho = author_adjusted_rho,
    p_adj = author_adjusted_p_adj,
    evidence_level = paste0("core_author_adjusted_", author_adjusted_level)
  )
]

within_state_v2 <- rbindlist(list(old_primary_detail, core_for_compare), fill = TRUE)
fwrite(within_state_v2, file.path(tables_dir, "within_state_primary_pathway_summary_with_core.csv"))

core_cont_v2 <- core_cont[
  ,
  .(
    dataset = "Core_GBmap_minus_Neftel2019",
    analysis_class = "IDHwt_GBM_harmonized_atlas",
    evidence_role = "atlas_level_supportive",
    malignant_rule = entry_variant,
    estimand,
    state = sub("^state_", "", state_score),
    n_patients = n_donors,
    n_patient_state = NA_integer_,
    spearman_rho,
    ci_low = NA_real_,
    ci_high = NA_real_,
    p_value,
    p_adj = p_adj_BH,
    fdr_family,
    evidence_level = classify_signal(spearman_rho, p_adj_BH, n_donors),
    source_result = "core_donor_level_continuous_state_score"
  )
]
old_state[, evidence_role := "independent_dataset"]
state_v2 <- rbindlist(list(old_state, core_cont_v2), fill = TRUE)
fwrite(state_v2, file.path(tables_dir, "state_continuum_summary_with_core.csv"))

interpretation <- data.table(
  claim = c(
    "LAP3 is associated with malignant-state continuum",
    "LAP3 is preferentially AC-like in Core GBmap",
    "LAP3-mTORC1 within-state coupling",
    "LAP3-BCAA within-state coupling",
    "MES-like LAP3-pathway coupling",
    "single-cell evidence proves LAP3-leucine-mTORC1 causality"
  ),
  v1_three_dataset = c(
    "supported_exploratory",
    "partial_AC_MES_context",
    "not_robust",
    "not_robust",
    "not_robust",
    "not_supported"
  ),
  core_gbmap_update = c(
    "strengthens_AC_like_signal",
    "supported_by_donor_fixed_effect_state_preference",
    "atlas_supportive_in_AC_NPC_OPC_after_author_adjustment",
    "limited_support_mainly_NPC_after_author_adjustment",
    "largely_author_structured",
    "still_not_supported"
  ),
  recommended_writing = c(
    "LAP3 marks an AC-like malignant-state continuum with context-dependent MES signals.",
    "Core GBmap supports LAP3 enrichment in AC-like malignant states after donor fixed effects.",
    "Core GBmap provides atlas-level supportive evidence for LAP3-mTORC1 transcriptional coupling, especially in AC/NPC/OPC states.",
    "BCAA coupling is weaker and should remain secondary/exploratory.",
    "MES results should be treated as author-sensitive and not emphasized.",
    "Keep causal wording out of the single-cell section."
  )
)
fwrite(interpretation, file.path(tables_dir, "meta_interpretation_matrix_v2_with_core.csv"))

readme <- c(
  "# LAP3 CellState MetaSummary v2 With Core GBmap",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Purpose",
  "",
  "This v2 result package adds Core GBmap as atlas-level supportive evidence.",
  "It does not overwrite the original three-dataset MetaSummary because Core GBmap is a harmonized multi-study atlas with stronger author/study structure.",
  "",
  "## Key decision",
  "",
  "Core GBmap strengthens the LAP3-AC/mTORC1 story, but does not make single-cell data a causal LAP3-leucine-mTORC1 proof.",
  "",
  "## Key outputs",
  "",
  "- `tables/core_author_adjusted_primary_summary.csv`",
  "- `tables/within_state_primary_pathway_summary_with_core.csv`",
  "- `tables/state_continuum_summary_with_core.csv`",
  "- `tables/meta_interpretation_matrix_v2_with_core.csv`",
  "",
  "## Recommended wording",
  "",
  "LAP3 was associated with an AC-like malignant-state continuum. Core GBmap, after excluding Neftel2019/GSE131928-derived cells, provided atlas-level supportive evidence for LAP3-mTORC1 transcriptional coupling in selected malignant states, whereas BCAA and MES-like signals were less stable and sensitive to author/study adjustment."
)
writeLines(readme, file.path(out_dir, "README.md"))

cat("Finished: ", format(Sys.time()), "\n", file = log_file, append = TRUE)
cat("Wrote:", out_dir, "\n")
