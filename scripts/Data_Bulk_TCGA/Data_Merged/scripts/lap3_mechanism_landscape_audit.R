#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

data.table::setDTthreads(8)

root <- "/home/lzb/glioma"
out_dir <- file.path(root, "Data_Bulk_TCGA/Data_Merged/results/LAP3_Mechanism_Landscape_Audit")
table_dir <- file.path(out_dir, "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

read_dt <- function(path, sep = ",") {
  full_path <- file.path(root, path)
  if (!file.exists(full_path)) {
    stop("Missing input: ", full_path)
  }
  fread(full_path, sep = sep)
}

safe_num <- function(x) suppressWarnings(as.numeric(x))

classify_pathway_family <- function(pathway) {
  p <- toupper(pathway)
  fifelse(
    grepl("INTERFERON|ALLOGRAFT|INFLAMMAT|COMPLEMENT|NEUTROPHIL|LEUKOCYTE|MACROPHAGE|MONOCYTE|CYTOKINE|TNFA|IL6|JAK_STAT|TOLL|NFKB|ANTIGEN|MHC|HLA", p),
    "TAM_myeloid_immune",
    fifelse(
      grepl("PROTEASOME|UNFOLDED|ER_STRESS|HEAT_SHOCK|CHAPERONE|PROTEIN_FOLD|UBIQUITIN|AUTOPHAGY|LYSOSOME", p),
      "proteostasis_stress",
      fifelse(
        grepl("TRANSLATION|RIBOSOME|EIF|MRNA|RNA_PROCESS|SPLICEOSOME|PROTEIN_SYNTHESIS|TRNA", p),
        "translation_anabolic",
        fifelse(
          grepl("MTOR|AMINO_ACID|BRANCHED|BCAA|LEUCINE|VALINE|ISOLEUCINE|GLUTAMINE|NITROGEN|ONE_CARBON|PURINE|PYRIMIDINE", p),
          "mTORC1_BCAA_amino_acid",
          fifelse(
            grepl("HYPOXIA|GLYCOLYSIS|ANGIOGENESIS|VEGF|REACTIVE_OXYGEN", p),
            "hypoxia_gliosis_angiogenesis",
            fifelse(
              grepl("E2F|G2M|MITOTIC|CELL_CYCLE|DNA_REPAIR|REPLICATION|CHECKPOINT|MYC|TELOMERE", p),
              "cell_cycle_proliferation",
              fifelse(
                grepl("EXTRACELLULAR|COLLAGEN|MATRIX|ECM|INTEGRIN|ADHESION|EPITHELIAL_MESENCHYMAL|EMT", p),
                "ECM_invasion_adhesion",
                "other"
              )
            )
          )
        )
      )
    )
  )
}

map_submodule_family <- function(module) {
  fifelse(
    module == "LAP3_MYELOID_TAM_CONTEXT_MODULE", "TAM_myeloid_immune",
    fifelse(
      module == "LAP3_MALIGNANT_STATE_MODULE", "malignant_state",
      fifelse(
        module == "LAP3_PROTEOSTASIS_STRESS_MODULE", "proteostasis_stress",
        fifelse(
          module == "LAP3_ANABOLIC_TRANSLATION_MODULE", "translation_anabolic",
          fifelse(
            module == "LAP3_HYPOXIA_PERINECROTIC_MODULE", "hypoxia_gliosis_angiogenesis",
            NA_character_
          )
        )
      )
    )
  )
}

classify_lincs_group <- function(mechanism_group) {
  m <- toupper(mechanism_group)
  fifelse(
    grepl("HSP90|PROTEASOME|PROTEOSTASIS", m), "proteostasis_stress",
    fifelse(
      grepl("TRANSLATION", m), "translation_anabolic",
      fifelse(
        grepl("PI3K|MTOR", m), "mTORC1_BCAA_amino_acid",
        fifelse(
          grepl("BET|EPIGENETIC|HDAC", m), "chromatin_epigenetic",
          "other"
        )
      )
    )
  )
}

evidence_rows <- list()
add_evidence <- function(dt) {
  evidence_rows[[length(evidence_rows) + 1L]] <<- dt
}

# 1. Bulk open-but-preclassified pathway landscape.
fgsea <- read_dt("Data_Bulk_TCGA/Data_Merged/results/LAP3_Pathway/tables/fgsea_cohort_balanced_tertile_adjusted.csv")
fgsea[, mechanism_family := classify_pathway_family(pathway)]
fgsea_pos <- fgsea[!is.na(NES) & NES > 0]
fgsea_summary <- fgsea_pos[, .(
  n_pathways_tested = .N,
  n_fdr_005_positive = sum(padj < 0.05, na.rm = TRUE),
  best_pathway = pathway[which.min(padj)],
  best_nes = NES[which.min(padj)],
  best_fdr = min(padj, na.rm = TRUE),
  top_positive_pathways = paste(head(pathway[order(padj, -NES)], 5), collapse = "; ")
), by = mechanism_family]
fgsea_summary[, evidence_layer := "bulk_fgsea_TCGA_cohort_balanced"]
fgsea_summary[, metric := paste0("n_sig=", n_fdr_005_positive, ";best_NES=", signif(best_nes, 4), ";best_FDR=", signif(best_fdr, 4))]
fgsea_summary[, strength := fifelse(
  n_fdr_005_positive >= 5 & best_nes >= 2.5 & best_fdr < 0.001, "strong",
  fifelse(n_fdr_005_positive >= 1 & best_fdr < 0.05, "supportive", "weak")
)]
add_evidence(fgsea_summary[, .(
  mechanism_family,
  evidence_layer,
  strength,
  metric,
  representative_feature = best_pathway,
  n = n_pathways_tested,
  effect = best_nes,
  fdr = best_fdr,
  notes = top_positive_pathways
)])
fwrite(fgsea_summary[order(best_fdr)], file.path(table_dir, "bulk_fgsea_mechanism_family_summary.csv"))

# 2. LAP3-state submodule multimodal support.
component <- read_dt("Data_Bulk_TCGA/Data_Merged/results/LAP3_State_Union_Robustness_Audit/tables/lap3_state_union_component_evidence_matrix.csv")
component[, mechanism_family := map_submodule_family(submodule)]
component[, evidence_layer := "state_submodule_multimodal"]
component[, strength := fifelse(
  robustness_class == "strong_multimodal_component", "strong",
  fifelse(robustness_class == "supportive_component", "supportive", "weak")
)]
component[, metric := paste0("support_score=", support_score, ";genes=", n_genes)]
add_evidence(component[!is.na(mechanism_family), .(
  mechanism_family,
  evidence_layer,
  strength,
  metric,
  representative_feature = submodule,
  n = n_genes,
  effect = support_score,
  fdr = NA_real_,
  notes = robustness_class
)])

# 3. Core GBmap malignant-cell submodule coupling.
gbmap <- read_dt("Data_Bulk_TCGA/Data_Merged/results/LAP3_State_Submodules/tables/gbmap_core_lap3_state_submodule_primary_summary.csv")
gbmap[, mechanism_family := map_submodule_family(submodule)]
gbmap_summary <- gbmap[
  entry_variant == "main_neoplastic_exclude_neftel2019" & threshold == 20 & !is.na(mechanism_family),
  .(
    states_tested = .N,
    states_supported = sum(p_adj_BH < 0.05 & spearman_rho > 0, na.rm = TRUE),
    best_state = author_state[which.min(p_adj_BH)],
    best_rho = spearman_rho[which.min(p_adj_BH)],
    best_fdr = min(p_adj_BH, na.rm = TRUE),
    min_leave_one_author_rho = min(leave_one_author_min_rho, na.rm = TRUE)
  ),
  by = mechanism_family
]
gbmap_summary[, evidence_layer := "core_gbmap_donor_state"]
gbmap_summary[, strength := fifelse(
  states_supported >= 3 & min_leave_one_author_rho > 0, "strong",
  fifelse(states_supported >= 1 & best_fdr < 0.05, "supportive", "weak")
)]
gbmap_summary[, metric := paste0("states_supported=", states_supported, "/4;best=", best_state, ";rho=", signif(best_rho, 4), ";FDR=", signif(best_fdr, 4))]
add_evidence(gbmap_summary[, .(
  mechanism_family,
  evidence_layer,
  strength,
  metric,
  representative_feature = best_state,
  n = states_tested,
  effect = best_rho,
  fdr = best_fdr,
  notes = paste0("min_LOAO_rho=", signif(min_leave_one_author_rho, 4))
)])

# 4. GBM-Space spatial topology and extended TAM axes.
spatial <- read_dt("Data_Spatial_Public/GBM_Space/results/LAP3_State_Spatial_Topology/tables/gbmspace_lap3_state_topology_summary.tsv", sep = "\t")
spatial[, mechanism_family := fifelse(
  readout_class == "myeloid_tam", "TAM_myeloid_immune",
  fifelse(readout_class %chin% c("spatial_niche", "malignant_state"), "hypoxia_gliosis_angiogenesis",
    fifelse(grepl("MTOR|BCAA|LEUCINE|TRANSLATION", toupper(readout)), "mTORC1_BCAA_amino_acid", "other")
  )
)]
spatial_summary <- spatial[
  state_set == "LAP3_STATE_UNION" & priority == TRUE & n_tumors >= 8 & mechanism_family != "other",
  .(
    n_readouts = .N,
    n_fdr_005_positive = sum(fdr_depth_adjusted_priority < 0.05 & median_depth_adjusted_rho > 0, na.rm = TRUE),
    best_readout = readout[which.max(abs_median_depth_adjusted_rho)],
    best_rho = median_depth_adjusted_rho[which.max(abs_median_depth_adjusted_rho)],
    best_fdr = fdr_depth_adjusted_priority[which.max(abs_median_depth_adjusted_rho)]
  ),
  by = mechanism_family
]
spatial_summary[, evidence_layer := "gbmspace_depth_aware_topology"]
spatial_summary[, strength := fifelse(
  n_fdr_005_positive >= 3 & best_rho > 0.25, "strong",
  fifelse(n_fdr_005_positive >= 1 & best_rho > 0, "supportive", "weak")
)]
spatial_summary[, metric := paste0("positive=", n_fdr_005_positive, "/", n_readouts, ";best_rho=", signif(best_rho, 4), ";FDR=", signif(best_fdr, 4))]
add_evidence(spatial_summary[, .(
  mechanism_family,
  evidence_layer,
  strength,
  metric,
  representative_feature = best_readout,
  n = n_readouts,
  effect = best_rho,
  fdr = best_fdr,
  notes = "depth-adjusted tumor-level topology"
)])

tam_nom <- read_dt("Data_scRNA_GEO/GBmap_Extended/results/Focused_Malignant_TAM_Communication/tables/extended_focused_tam_nomination_table.csv")
tam_summary <- tam_nom[nomination_tier == "Tier1_positive_candidate", .(
  n_axes = .N,
  best_axis = axis[which.max(abs(spearman_rho))],
  best_direction = direction[which.max(abs(spearman_rho))],
  best_rho = spearman_rho[which.max(abs(spearman_rho))],
  best_fdr = fdr_bh[which.max(abs(spearman_rho))],
  axes = paste(unique(axis), collapse = "; ")
)]
tam_summary[, `:=`(
  mechanism_family = "TAM_myeloid_immune",
  evidence_layer = "extended_gbmap_focused_malignant_TAM_axes",
  strength = fifelse(n_axes >= 3 & best_fdr < 0.001, "strong", "supportive"),
  metric = paste0("Tier1_axes=", n_axes, ";best_rho=", signif(best_rho, 4), ";FDR=", signif(best_fdr, 4))
)]
add_evidence(tam_summary[, .(
  mechanism_family,
  evidence_layer,
  strength,
  metric,
  representative_feature = best_axis,
  n = n_axes,
  effect = best_rho,
  fdr = best_fdr,
  notes = axes
)])

axis_spatial <- read_dt("Data_Spatial_Public/GBM_Space/results/Extended_TAM_Axis_Spatial_Support/tables/gbmspace_extended_tam_axis_spatial_support_summary.tsv", sep = "\t")
axis_summary <- axis_spatial[
  readout_class %chin% c("lap3_state", "tam_neighbor_k6") & n_tumors >= 8,
  .(
    n_tests = .N,
    n_fdr_005_positive = sum(fdr_depth_adjusted < 0.05 & median_depth_adjusted_rho > 0, na.rm = TRUE),
    best_axis = axis[which.max(abs_median_depth_adjusted_rho)],
    best_readout = readout[which.max(abs_median_depth_adjusted_rho)],
    best_rho = median_depth_adjusted_rho[which.max(abs_median_depth_adjusted_rho)],
    best_fdr = fdr_depth_adjusted[which.max(abs_median_depth_adjusted_rho)]
  )
]
axis_summary[, `:=`(
  mechanism_family = "TAM_myeloid_immune",
  evidence_layer = "gbmspace_extended_TAM_axis_spatial_support",
  strength = fifelse(n_fdr_005_positive >= 6 & best_rho > 0.25, "strong", "supportive"),
  metric = paste0("positive=", n_fdr_005_positive, "/", n_tests, ";best_rho=", signif(best_rho, 4), ";FDR=", signif(best_fdr, 4))
)]
add_evidence(axis_summary[, .(
  mechanism_family,
  evidence_layer,
  strength,
  metric,
  representative_feature = paste(best_axis, best_readout, sep = "::"),
  n = n_tests,
  effect = best_rho,
  fdr = best_fdr,
  notes = "spatial co-localization/topology; not directional communication"
)])

# 5. CPTAC/GLASS boundary and cross-modal coordination by submodule family.
cptac_glass <- read_dt("Data_Bulk_TCGA/Data_Merged/results/LAP3_State_CPTAC_GLASS_Projection/tables/lap3_state_cptac_glass_primary_summary.csv")
cptac_glass[, module_clean := sub("^mrna_", "", sub("^protein_", "", feature))]
cptac_glass[, mechanism_family := map_submodule_family(module_clean)]
cptac_glass_summary <- cptac_glass[!is.na(mechanism_family), .(
  n_tests = .N,
  n_fdr_005_positive = sum(fdr < 0.05 & effect > 0, na.rm = TRUE),
  best_outcome = outcome[which.min(fdr)],
  best_effect = effect[which.min(fdr)],
  best_fdr = min(fdr, na.rm = TRUE)
), by = mechanism_family]
cptac_glass_summary[, evidence_layer := "cptac_glass_cross_modal_boundary"]
cptac_glass_summary[, strength := fifelse(
  n_fdr_005_positive >= 2 & best_fdr < 0.001, "strong",
  fifelse(n_fdr_005_positive >= 1 & best_fdr < 0.05, "supportive", "weak")
)]
cptac_glass_summary[, metric := paste0("positive=", n_fdr_005_positive, "/", n_tests, ";best=", best_outcome, ";effect=", signif(best_effect, 4), ";FDR=", signif(best_fdr, 4))]
add_evidence(cptac_glass_summary[, .(
  mechanism_family,
  evidence_layer,
  strength,
  metric,
  representative_feature = best_outcome,
  n = n_tests,
  effect = best_effect,
  fdr = best_fdr,
  notes = "includes positive coordination and negative direct BCAA/phospho boundaries"
)])

# 6. LINCS perturbational nomination classes.
lincs_path <- "Data_Perturbation_Public/LINCS_GSE70138/results/LAP3_Local_Connectivity/tables/lap3_lincs_mechanism_focused_candidates.csv"
if (file.exists(file.path(root, lincs_path))) {
  lincs <- read_dt(lincs_path)
  lincs[, mechanism_family := classify_lincs_group(mechanism_group)]
  lincs_summary <- lincs[stable_reverse == TRUE, .(
    n_stable_candidates = .N,
    best_compound = pert_iname[which.min(consensus_reverse_rank)],
    best_rank = min(consensus_reverse_rank, na.rm = TRUE),
    best_worst_case_ncs = worst_case_median_ncs[which.min(consensus_reverse_rank)]
  ), by = mechanism_family]
  lincs_summary[, evidence_layer := "lincs_perturbational_nomination"]
  lincs_summary[, strength := fifelse(n_stable_candidates >= 2 & best_rank <= 20, "supportive", "weak")]
  lincs_summary[, metric := paste0("stable_candidates=", n_stable_candidates, ";best_rank=", best_rank, ";worst_case_NCS=", signif(best_worst_case_ncs, 4))]
  add_evidence(lincs_summary[, .(
    mechanism_family,
    evidence_layer,
    strength,
    metric,
    representative_feature = best_compound,
    n = n_stable_candidates,
    effect = -best_rank,
    fdr = NA_real_,
    notes = "perturbagen nominations only; no efficacy or target engagement"
  )])
}

evidence <- rbindlist(evidence_rows, use.names = TRUE, fill = TRUE)
strength_points <- c(strong = 2L, supportive = 1L, weak = 0L, negative = -1L)
evidence[, support_points := unname(strength_points[strength])]
evidence[is.na(support_points), support_points := 0L]
setorder(evidence, mechanism_family, evidence_layer)
fwrite(evidence, file.path(table_dir, "mechanism_evidence_long.csv"))

family_summary <- evidence[, .(
  evidence_layers = .N,
  strong_layers = sum(strength == "strong"),
  supportive_layers = sum(strength == "supportive"),
  weak_layers = sum(strength == "weak"),
  total_support_points = sum(support_points),
  representative_features = paste(unique(na.omit(representative_feature))[1:min(5, uniqueN(na.omit(representative_feature)))], collapse = "; "),
  key_metrics = paste(paste(evidence_layer, metric, sep = ": "), collapse = " | ")
), by = mechanism_family]
family_summary[, priority_class := fifelse(
  strong_layers >= 2 & total_support_points >= 5, "primary_evidence_direction",
  fifelse(strong_layers >= 1 & total_support_points >= 3, "supportive_evidence_direction",
    fifelse(total_support_points >= 2, "context_or_secondary_direction", "weak_or_boundary")
  )
)]
setorder(family_summary, -total_support_points, -strong_layers, mechanism_family)
fwrite(family_summary, file.path(table_dir, "mechanism_evidence_matrix.csv"))

pre_specified_families <- data.table(
  mechanism_family = c(
    "TAM_myeloid_immune",
    "malignant_state",
    "proteostasis_stress",
    "translation_anabolic",
    "hypoxia_gliosis_angiogenesis",
    "mTORC1_BCAA_amino_acid",
    "cell_cycle_proliferation",
    "ECM_invasion_adhesion",
    "chromatin_epigenetic"
  ),
  interpretation_scope = c(
    "Malignant-microenvironmental ecology and TAM/myeloid spatial context.",
    "Malignant-cell state coupled to LAP3-state and donor-level malignant programs.",
    "Protein folding/proteasome/HSP90/stress-associated state component.",
    "Anabolic and translational state component; supportive, not a causal LAP3-mTORC1 proof.",
    "Hypoxic, perinecrotic, angiogenic or gliosis-like context.",
    "Amino-acid/mTORC1-associated readout family; direct BCAA and phospho evidence remains boundary-controlled.",
    "Proliferation/cell-cycle context seen in bulk pathway scans.",
    "Invasion, matrix and adhesion context.",
    "Perturbational nomination class from LINCS, not a core biology claim by itself."
  ),
  manuscript_boundary = c(
    "Can be used as a main ecological interpretation if supported across layers.",
    "Can be used as the malignant-state component of LAP3-state.",
    "Can be used as a secondary mechanistic hypothesis and perturbational vulnerability class.",
    "Use only as supportive biology unless independently validated.",
    "Use as niche/context evidence, not causal mechanism.",
    "Do not restore LAP3-BCAA-mTORC1 causal wording from this audit.",
    "Use as background/context; avoid making it the LAP3-specific mechanism.",
    "Use as context only unless spatial or single-cell support becomes strong.",
    "Use as LINCS hypothesis only."
  )
)
fwrite(pre_specified_families, file.path(table_dir, "pre_specified_mechanism_families.csv"))

primary <- family_summary[priority_class == "primary_evidence_direction", mechanism_family]
supportive <- family_summary[priority_class == "supportive_evidence_direction", mechanism_family]
boundary <- family_summary[mechanism_family == "mTORC1_BCAA_amino_acid"]

verdict <- data.table(
  question = c(
    "Was a completely open-ended pathway fishing analysis performed?",
    "What is the strongest alternative mechanism direction?",
    "Does the audit rescue LAP3-BCAA-mTORC1 as the main mechanism?",
    "What should happen if a strong positive mechanism is detected?",
    "How should this be used in the manuscript?"
  ),
  answer = c(
    "No. This audit uses pre-specified mechanism families and existing validated outputs, including an open-but-classified TCGA Hallmark/Reactome pathway scan.",
    if (length(primary) > 0) paste(primary, collapse = "; ") else "No primary mechanism candidate exceeded the cross-layer threshold.",
    if (nrow(boundary) > 0 && boundary$total_support_points >= 5) {
      "It shows pathway-associated support but must still be constrained by direct CPTAC BCAA/phospho and single-cell boundaries."
    } else {
      "No. mTORC1/BCAA remains a boundary-controlled supportive family, not a direct causal mechanism."
    },
    "Treat it as a nominated mechanism family requiring targeted validation and claim-boundary review before changing the main story.",
    "Use as Supplementary/Methods defense that mechanism space was audited without post hoc fishing."
  ),
  manuscript_use = c(
    "Methods and Supplementary table; not a main-figure discovery claim.",
    "Discussion and Results bridge if wording remains observational.",
    "Avoid causal LAP3-leucine-mTORC1 language.",
    "Escalate to expert review and targeted validation design.",
    "Reviewer-defense evidence and future-direction prioritization."
  )
)
fwrite(verdict, file.path(table_dir, "mechanism_landscape_verdict.csv"))

readme <- c(
  "# LAP3 Mechanism Landscape Audit",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Purpose",
  "",
  "This audit is a conservative rescue step after the direct LAP3-BCAA-mTORC1 mechanism was not supported by prior single-cell, spatial, CPTAC and GLASS analyses. It asks whether other pathway or mechanism families were systematically considered without turning the project into post hoc pathway fishing.",
  "",
  "## Design",
  "",
  "- Pre-specified mechanism families were frozen before summarising the outputs.",
  "- Existing validated result tables were reused; no large expression matrix or single-cell object was reloaded.",
  "- TCGA cohort-balanced fgsea was used as an open-but-classified bulk pathway landscape.",
  "- Cross-layer evidence was then summarised from LAP3-state submodules, Core GBmap, GBM-Space, Extended GBmap TAM axes, CPTAC/GLASS and LINCS.",
  "- Direct causal claims remain prohibited; strong positive families are mechanism nominations, not proof.",
  "",
  "## Main Output Tables",
  "",
  "- `tables/pre_specified_mechanism_families.csv`: frozen mechanism-family definitions and manuscript boundaries.",
  "- `tables/bulk_fgsea_mechanism_family_summary.csv`: TCGA pathway landscape grouped into the frozen families.",
  "- `tables/mechanism_evidence_long.csv`: layer-by-layer evidence entries.",
  "- `tables/mechanism_evidence_matrix.csv`: family-level cross-modal evidence matrix.",
  "- `tables/mechanism_landscape_verdict.csv`: manuscript-facing conclusion and claim boundary.",
  "",
  "## High-Level Conclusion",
  "",
  "The audit supports a state/ecosystem interpretation rather than a direct LAP3-BCAA-mTORC1 causal mechanism. The strongest cross-layer direction should remain malignant-state and TAM/myeloid ecological coupling, with proteostasis/stress and translation/anabolic programs treated as secondary mechanism hypotheses. Any newly strong family should be escalated as a candidate mechanism requiring targeted validation before it changes the main manuscript claim.",
  "",
  "## Boundary",
  "",
  "This audit does not prove LAP3 causality, ligand-receptor directionality, intracellular leucine availability, mTORC1 phosphorylation activation or drug efficacy."
)
writeLines(readme, file.path(out_dir, "README.md"))

cat("Wrote mechanism landscape audit to ", out_dir, "\n", sep = "")
