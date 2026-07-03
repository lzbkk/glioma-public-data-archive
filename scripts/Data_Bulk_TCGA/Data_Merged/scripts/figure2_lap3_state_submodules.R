suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

required_pkgs <- c("svglite", "ragg")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))

set.seed(20260701)
data.table::setDTthreads(8)

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)][1])
if (is.na(file_arg) || !nzchar(file_arg)) {
  file_arg <- "Data_Bulk_TCGA/Data_Merged/scripts/figure2_lap3_state_submodules.R"
}
repo <- normalizePath(file.path(dirname(file_arg), "../../.."), mustWork = FALSE)
if (!dir.exists(file.path(repo, "Project_Management"))) repo <- normalizePath(".", mustWork = TRUE)

state_dir <- file.path(repo, "Data_Bulk_TCGA/Data_Merged/results/LAP3_State_Module")
submodule_dir <- file.path(repo, "Data_Bulk_TCGA/Data_Merged/results/LAP3_State_Submodules")
malignant_dir <- file.path(repo, "Data_Bulk_TCGA/Data_Merged/results/LAP3_Malignant_State_Module_Audit")
out_dir <- file.path(repo, "Data_Bulk_TCGA/Data_Merged/results/Figure2_LAP3_State_Submodules")
source_dir <- file.path(out_dir, "source_data")
plot_dir <- file.path(out_dir, "plots")
table_dir <- file.path(out_dir, "tables")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

read_dt <- function(path) {
  if (!file.exists(path)) stop("Missing input: ", path)
  fread(path)
}
write_source <- function(dt, name) fwrite(copy(dt), file.path(source_dir, name))
write_table <- function(dt, name) fwrite(copy(dt), file.path(table_dir, name))

label_dataset <- c(
  "TCGA" = "TCGA",
  "CGGA_mRNAseq_693" = "CGGA693",
  "CGGA_mRNAseq_325" = "CGGA325"
)
clean_dataset <- function(x) {
  out <- as.character(x)
  idx <- out %in% names(label_dataset)
  out[idx] <- unname(label_dataset[out[idx]])
  out
}

submodule_labels <- c(
  "LAP3_MALIGNANT_STATE_MODULE" = "Malignant-state",
  "LAP3_MYELOID_TAM_CONTEXT_MODULE" = "Myeloid/TAM",
  "LAP3_ANABOLIC_TRANSLATION_MODULE" = "Anabolic/translation",
  "LAP3_PROTEOSTASIS_STRESS_MODULE" = "Proteostasis/stress",
  "LAP3_HYPOXIA_PERINECROTIC_MODULE" = "Hypoxia/perinecrotic"
)
submodule_order <- names(submodule_labels)
submodule_colors <- c(
  "Malignant-state" = "#315E8A",
  "Myeloid/TAM" = "#B15A2C",
  "Anabolic/translation" = "#7A9B55",
  "Proteostasis/stress" = "#7B5B8F",
  "Hypoxia/perinecrotic" = "#8A6F3D"
)
source_colors <- c(
  "TCGA only" = "#6E91B8",
  "GBmap only" = "#8DBB77",
  "TCGA + GBmap" = "#2F6F9F"
)

theme_fig <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.3, colour = "#2A2D31"),
      axis.line = element_line(linewidth = 0.25, colour = "#343941"),
      axis.ticks = element_line(linewidth = 0.25, colour = "#343941"),
      strip.background = element_rect(fill = "#F2F4F6", colour = NA),
      strip.text = element_text(size = base_size, face = "bold"),
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size - 0.3),
      plot.title = element_text(size = base_size + 1, face = "bold", hjust = 0),
      plot.margin = margin(4, 4, 4, 4)
    )
}

frozen_genes <- read_dt(file.path(state_dir, "tables/lap3_state_frozen_gene_sets.csv"))
gene_counts <- read_dt(file.path(state_dir, "tables/lap3_state_gene_set_counts.csv"))
submodule_assign <- read_dt(file.path(submodule_dir, "tables/lap3_state_submodule_gene_assignment.csv"))
submodule_counts <- read_dt(file.path(submodule_dir, "tables/lap3_state_submodule_gene_counts.csv"))
submodule_cor <- read_dt(file.path(submodule_dir, "tables/lap3_state_submodule_bulk_correlations.csv"))
submodule_coverage <- read_dt(file.path(submodule_dir, "tables/lap3_state_submodule_bulk_gene_coverage.csv"))
gbmap_submodule <- read_dt(file.path(submodule_dir, "tables/gbmap_core_lap3_state_submodule_primary_summary.csv"))
malignant_clusters <- read_dt(file.path(malignant_dir, "tables/malignant_module_gene_clusters.csv"))
malignant_bulk <- read_dt(file.path(malignant_dir, "tables/bulk_malignant_cluster_score_correlations.csv"))
malignant_gbmap <- read_dt(file.path(malignant_dir, "tables/gbmap_core_malignant_cluster_donor_state_associations.csv"))
malignant_space <- read_dt(file.path(malignant_dir, "tables/gbmspace_malignant_cluster_spatial_summary.tsv"))
malignant_verdict <- read_dt(file.path(malignant_dir, "tables/malignant_state_module_audit_verdict.csv"))

union_genes <- unique(frozen_genes[state_set == "LAP3_STATE_UNION"], by = "gene")
panel_a <- union_genes[, .(
  n_genes = .N
), by = .(
  source = fifelse(in_tcga_top150 & in_gbmap_up, "TCGA + GBmap",
    fifelse(in_tcga_top150, "TCGA only", "GBmap only")
  )
)]
panel_a[, state_set := "LAP3_STATE_UNION"]
panel_a[, source := factor(source, levels = names(source_colors))]
panel_a[, total_union := sum(n_genes)]
panel_a <- panel_a[order(source)]

panel_b <- unique(submodule_assign[primary_submodule %chin% submodule_order], by = "gene")
panel_b[, submodule_label := submodule_labels[primary_submodule]]
panel_b[, source := fifelse(in_tcga_top150 & in_gbmap_up, "TCGA + GBmap",
  fifelse(in_tcga_top150, "TCGA only", "GBmap only")
)]
panel_b <- panel_b[, .(n_genes = .N), by = .(primary_submodule, submodule_label, source)]
panel_b[, submodule_label := factor(submodule_label, levels = rev(submodule_labels[submodule_order]))]
panel_b[, source := factor(source, levels = names(source_colors))]
panel_b <- panel_b[order(submodule_label, source)]

bulk_vars <- c(
  "LAP3_log2_expr",
  "LAP3_STATE_UNION",
  "HALLMARK_MTORC1_SIGNALING",
  "LEUCINE_BCAA_CORE",
  "REACTOME_TRANSLATION",
  "TAM_MYELOID_CORE",
  "HALLMARK_HYPOXIA"
)
bulk_var_labels <- c(
  "LAP3_log2_expr" = "LAP3",
  "LAP3_STATE_UNION" = "Union",
  "HALLMARK_MTORC1_SIGNALING" = "mTORC1",
  "LEUCINE_BCAA_CORE" = "BCAA",
  "REACTOME_TRANSLATION" = "Translation",
  "TAM_MYELOID_CORE" = "TAM/myeloid",
  "HALLMARK_HYPOXIA" = "Hypoxia"
)
panel_c_raw <- submodule_cor[
  group == "all" & submodule %chin% submodule_order & variable %chin% bulk_vars
]
panel_c_raw[, dataset_label := clean_dataset(dataset)]
panel_c_raw[, submodule_label := submodule_labels[submodule]]
panel_c_raw[, variable_label := bulk_var_labels[variable]]
panel_c <- panel_c_raw[, .(
  median_rho = median(spearman_rho, na.rm = TRUE),
  min_rho = min(spearman_rho, na.rm = TRUE),
  max_rho = max(spearman_rho, na.rm = TRUE),
  min_fdr = min(p_adj_BH, na.rm = TRUE),
  n_cohorts = uniqueN(dataset_label)
), by = .(submodule, submodule_label, variable, variable_label)]
panel_c[, submodule_label := factor(submodule_label, levels = rev(submodule_labels[submodule_order]))]
panel_c[, variable_label := factor(variable_label, levels = bulk_var_labels[bulk_vars])]
panel_c[, label := sprintf("%.2f", median_rho)]

panel_d <- gbmap_submodule[
  entry_variant == "main_neoplastic_exclude_neftel2019" &
    threshold == 20 &
    submodule %chin% submodule_order
]
panel_d[, submodule_label := submodule_labels[submodule]]
panel_d[, author_state := factor(author_state, levels = c("AC", "MES", "NPC", "OPC"))]
panel_d[, submodule_label := factor(submodule_label, levels = rev(submodule_labels[submodule_order]))]
panel_d[, label := sprintf("%.2f", spearman_rho)]

panel_e <- malignant_clusters[, .(
  n_genes = .N,
  label = first(cluster_label)
), by = malignant_cluster]
panel_e[, cluster_label := paste0(malignant_cluster, "\n", sub("_", "-", label), "\n", n_genes, " genes")]
panel_e[, display_label := fifelse(malignant_cluster == "M1", "M1\nMES-like",
  fifelse(malignant_cluster == "M2", "M2\nAC-like", "M3\nminor arm")
)]
panel_e[, count_label := paste0(n_genes, " genes")]
panel_e[, malignant_cluster := factor(malignant_cluster, levels = c("M1", "M2", "M3"))]

bulk_malignant_summary <- malignant_bulk[
  stratum == "all" &
    cluster_score %chin% paste0("LAP3_MALIGNANT_M", 1:3) &
    variable %chin% c("LAP3_log2_expr", "LAP3_STATE_UNION")
][, .(
  effect = median(rho, na.rm = TRUE),
  min_effect = min(rho, na.rm = TRUE),
  max_effect = max(rho, na.rm = TRUE),
  min_fdr = min(p_adj_BH, na.rm = TRUE),
  n_units = median(n, na.rm = TRUE)
), by = .(cluster_score, variable)]
bulk_malignant_summary[, modality := fifelse(variable == "LAP3_log2_expr", "Bulk: LAP3", "Bulk: union")]

gbmap_malignant_summary <- malignant_gbmap[
  entry_variant == "main_neoplastic" &
    target == "lap3_mean" &
    cluster_score %chin% paste0("LAP3_MALIGNANT_M", 1:3)
][, .(
  effect = median(rho, na.rm = TRUE),
  min_effect = min(rho, na.rm = TRUE),
  max_effect = max(rho, na.rm = TRUE),
  min_fdr = min(p_adj_BH, na.rm = TRUE),
  n_units = median(n, na.rm = TRUE)
), by = cluster_score]
gbmap_malignant_summary[, modality := "Core GBmap: LAP3"]

space_malignant_summary <- malignant_space[
  target == "LAP3_STATE_UNION" &
    cluster_score %chin% paste0("LAP3_MALIGNANT_M", 1:3)
][, .(
  effect = median_depth_adjusted_rho,
  min_effect = median_depth_adjusted_rho,
  max_effect = median_depth_adjusted_rho,
  min_fdr = p_adj_BH,
  n_units = n_sections
), by = cluster_score]
space_malignant_summary[, modality := "GBM-Space: union"]

panel_f <- rbindlist(list(
  bulk_malignant_summary[, .(cluster_score, modality, effect, min_effect, max_effect, min_fdr, n_units)],
  gbmap_malignant_summary[, .(cluster_score, modality, effect, min_effect, max_effect, min_fdr, n_units)],
  space_malignant_summary[, .(cluster_score, modality, effect, min_effect, max_effect, min_fdr, n_units)]
), fill = TRUE)
panel_f[, cluster_label := sub("LAP3_MALIGNANT_", "", cluster_score)]
panel_f[, cluster_label := factor(cluster_label, levels = c("M1", "M2", "M3"))]
panel_f[, modality := factor(modality, levels = c("Bulk: LAP3", "Bulk: union", "Core GBmap: LAP3", "GBM-Space: union"))]
panel_f[, label := sprintf("%.2f", effect)]

panel_coverage <- copy(submodule_coverage)
panel_coverage[, dataset_label := clean_dataset(dataset)]
panel_coverage[, submodule_label := submodule_labels[submodule]]

panel_map <- data.table(
  panel = LETTERS[1:6],
  source_data = c(
    "figure2_panel_a_frozen_union_sources.csv",
    "figure2_panel_b_submodule_gene_sources.csv",
    "figure2_panel_c_bulk_submodule_context_correlations.csv",
    "figure2_panel_d_gbmap_submodule_lap3_coupling.csv",
    "figure2_panel_e_malignant_module_cluster_sizes.csv",
    "figure2_panel_f_malignant_cluster_cross_modal_support.csv"
  ),
  conclusion = c(
    "The frozen 207-gene union is built from TCGA and Core GBmap evidence with LAP3 excluded.",
    "The union separates into interpretable malignant-state, myeloid/TAM, anabolic/translation, proteostasis/stress and hypoxia/perinecrotic components.",
    "Bulk cohorts show distinct but coordinated submodule associations with LAP3, union and biological context readouts.",
    "Core GBmap supports donor-state-level LAP3 coupling mainly for malignant-state and selected context modules.",
    "The 144-gene malignant-state component is decomposed into three scoreable co-expression clusters.",
    "Malignant clusters show cross-modal support, with M1 as the dominant component and M2/M3 as smaller context-sensitive components."
  )
)

key_results <- rbindlist(list(
  data.table(panel = "A", feature = "LAP3_STATE_UNION", statistic = "n genes", value = sum(panel_a$n_genes), note = "LAP3 excluded"),
  panel_b[, .(panel = "B", feature = as.character(submodule_label), statistic = "n genes", value = sum(n_genes), note = "primary submodule") , by = primary_submodule][, primary_submodule := NULL],
  panel_c[variable %chin% c("LAP3_STATE_UNION", "HALLMARK_MTORC1_SIGNALING", "LEUCINE_BCAA_CORE"),
          .(panel = "C", feature = paste(as.character(submodule_label), variable_label, sep = " vs "), statistic = "median Spearman rho across cohorts", value = median_rho, note = paste0("range ", sprintf("%.2f", min_rho), " to ", sprintf("%.2f", max_rho)))],
  panel_d[submodule == "LAP3_MALIGNANT_STATE_MODULE",
          .(panel = "D", feature = paste("GBmap", author_state, "malignant-state"), statistic = "Spearman rho with LAP3 donor-state mean", value = spearman_rho, note = paste0("FDR=", signif(p_adj_BH, 3)))],
  panel_e[, .(panel = "E", feature = as.character(malignant_cluster), statistic = "n genes", value = n_genes, note = label)],
  panel_f[, .(panel = "F", feature = paste(cluster_label, modality, sep = " / "), statistic = "effect", value = effect, note = paste0("FDR=", signif(min_fdr, 3)))]
), fill = TRUE)

write_source(panel_a, "figure2_panel_a_frozen_union_sources.csv")
write_source(panel_b, "figure2_panel_b_submodule_gene_sources.csv")
write_source(panel_c_raw, "figure2_panel_c_bulk_submodule_context_correlations_raw.csv")
write_source(panel_c, "figure2_panel_c_bulk_submodule_context_correlations.csv")
write_source(panel_d, "figure2_panel_d_gbmap_submodule_lap3_coupling.csv")
write_source(panel_e, "figure2_panel_e_malignant_module_cluster_sizes.csv")
write_source(panel_f, "figure2_panel_f_malignant_cluster_cross_modal_support.csv")
write_table(panel_map, "figure2_panel_map.csv")
write_table(key_results, "figure2_key_results.csv")
write_table(panel_coverage, "figure2_submodule_bulk_gene_coverage_input.csv")
write_table(malignant_verdict, "figure2_malignant_state_module_audit_verdict_input.csv")

p_a <- ggplot(panel_a, aes(x = state_set, y = n_genes, fill = source)) +
  geom_col(width = 0.55, colour = "white", linewidth = 0.25) +
  geom_text(aes(label = paste0(n_genes, "\n", source)), position = position_stack(vjust = 0.5),
            colour = "white", size = 2.05, fontface = "bold", lineheight = 0.82) +
  annotate("text", x = 1, y = sum(panel_a$n_genes) + 12, label = paste0("207 genes\nLAP3 excluded"), size = 2.4, fontface = "bold", lineheight = 0.9) +
  scale_fill_manual(values = source_colors, name = NULL) +
  scale_y_continuous(limits = c(0, 235), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "A  Frozen LAP3-state union", x = NULL, y = "Genes") +
  theme_fig() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none")

p_b <- ggplot(panel_b, aes(x = n_genes, y = submodule_label, fill = source)) +
  geom_col(width = 0.68, colour = "white", linewidth = 0.25) +
  geom_text(
    data = panel_b[, .(n_genes = sum(n_genes)), by = submodule_label],
    aes(x = n_genes + 4, y = submodule_label, label = n_genes),
    inherit.aes = FALSE, hjust = 0, size = 2.2, colour = "#2A2D31"
  ) +
  scale_fill_manual(values = source_colors, name = NULL) +
  scale_x_continuous(limits = c(0, 160), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "B  Interpretable submodules", x = "Genes", y = NULL) +
  theme_fig() +
  theme(legend.position = "none")

p_c <- ggplot(panel_c, aes(x = variable_label, y = submodule_label, fill = median_rho)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_text(aes(label = label), size = 1.85, colour = "#1E252B") +
  scale_fill_gradient2(low = "#B7C7D6", mid = "#F3F4F2", high = "#B15A2C", midpoint = 0.45,
                       limits = c(0, 1), oob = squish, name = "median rho") +
  labs(title = "C  Bulk submodule context", x = NULL, y = NULL) +
  theme_fig() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1), legend.position = "right")

p_d <- ggplot(panel_d, aes(x = author_state, y = submodule_label, fill = spearman_rho)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_text(aes(label = label), size = 1.85, colour = "#1E252B") +
  scale_fill_gradient2(low = "#D8DEE6", mid = "#F4F4F0", high = "#315E8A", midpoint = 0.25,
                       limits = c(-0.3, 0.8), oob = squish, name = "rho") +
  labs(title = "D  Core GBmap coupling", x = "Cell state", y = NULL) +
  theme_fig() +
  theme(legend.position = "right")

p_e <- ggplot(panel_e, aes(x = malignant_cluster, y = n_genes, fill = malignant_cluster)) +
  geom_col(width = 0.58, colour = "white", linewidth = 0.25) +
  geom_text(aes(label = count_label), vjust = -0.25, size = 2.2) +
  scale_fill_manual(values = c("M1" = "#315E8A", "M2" = "#7B5B8F", "M3" = "#8A6F3D"), guide = "none") +
  scale_y_continuous(limits = c(0, 140), expand = expansion(mult = c(0, 0.03))) +
  scale_x_discrete(labels = setNames(as.character(panel_e$display_label), as.character(panel_e$malignant_cluster))) +
  labs(title = "E  Malignant module split", x = NULL, y = "Genes") +
  theme_fig() +
  theme(axis.text.x = element_text(lineheight = 0.9))

p_f <- ggplot(panel_f, aes(x = modality, y = cluster_label, fill = effect)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_text(aes(label = label), size = 1.95, colour = "#1E252B") +
  scale_fill_gradient2(low = "#D9E1EA", mid = "#F4F4F0", high = "#315E8A", midpoint = 0.45,
                       limits = c(-0.2, 1), oob = squish, name = "effect") +
  labs(title = "F  Cross-modal malignant support", x = NULL, y = NULL) +
  theme_fig() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1), legend.position = "right")

caption <- paste(
  "Figure 2 decomposes the frozen LAP3-state into interpretable components and audits the malignant-state arm.\n",
  "These panels support composite interpretability and cross-modal consistency, not a new causal LAP3 mechanism."
)

figure2 <- (p_a | p_b) / (p_c | p_d) / (p_e | p_f) +
  plot_layout(heights = c(0.9, 1.15, 0.95), guides = "keep") +
  plot_annotation(
    title = "Figure 2. The frozen LAP3-state is reproducible and interpretable",
    caption = caption,
    theme = theme(
      plot.title = element_text(face = "bold", size = 10, family = "sans"),
      plot.caption = element_text(size = 6.5, colour = "#3A3D42", hjust = 0, family = "sans"),
      plot.margin = margin(5, 5, 5, 5)
    )
  )

plot_base <- file.path(plot_dir, "Figure2_LAP3_State_Submodules")
svglite::svglite(paste0(plot_base, ".svg"), width = 7.2, height = 7.6)
print(figure2)
dev.off()
ggsave(paste0(plot_base, ".pdf"), figure2, width = 7.2, height = 7.6, units = "in", device = grDevices::pdf)
ragg::agg_tiff(paste0(plot_base, ".tiff"), width = 7.2, height = 7.6, units = "in", res = 600, compression = "lzw")
print(figure2)
dev.off()
ragg::agg_png(paste0(plot_base, ".png"), width = 7.2, height = 7.6, units = "in", res = 300)
print(figure2)
dev.off()

export_files <- c(
  paste0(plot_base, c(".svg", ".pdf", ".tiff", ".png")),
  file.path(source_dir, panel_map$source_data),
  file.path(source_dir, "figure2_panel_c_bulk_submodule_context_correlations_raw.csv"),
  file.path(table_dir, c(
    "figure2_panel_map.csv",
    "figure2_key_results.csv",
    "figure2_submodule_bulk_gene_coverage_input.csv",
    "figure2_malignant_state_module_audit_verdict_input.csv"
  ))
)
export_qc <- data.table(
  file = export_files,
  exists = file.exists(export_files),
  bytes = as.numeric(file.info(export_files)$size)
)
write_table(export_qc, "figure2_export_qc.csv")

readme <- c(
  "# Figure 2 LAP3 State Submodules",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Purpose",
  "",
  "Create the active manuscript Figure 2 for the state/ecosystem route.",
  "",
  "## Core Conclusion",
  "",
  "The frozen `LAP3_STATE_UNION` is a reproducible and interpretable composite transcriptional state. It separates into malignant-state, myeloid/TAM, anabolic/translation, proteostasis/stress and hypoxia/perinecrotic components, and its 144-gene malignant-state arm can be further decomposed into scoreable M1/M2/M3 clusters.",
  "",
  "## Inputs",
  "",
  "- `Data_Bulk_TCGA/Data_Merged/results/LAP3_State_Module/`",
  "- `Data_Bulk_TCGA/Data_Merged/results/LAP3_State_Submodules/`",
  "- `Data_Bulk_TCGA/Data_Merged/results/LAP3_Malignant_State_Module_Audit/`",
  "",
  "## Panels",
  "",
  "- Panel A: frozen `LAP3_STATE_UNION` source evidence.",
  "- Panel B: five interpretable submodule gene-source composition.",
  "- Panel C: bulk cohort median submodule-context correlations.",
  "- Panel D: Core GBmap donor-state LAP3 coupling for submodules.",
  "- Panel E: malignant-state module M1/M2/M3 split.",
  "- Panel F: cross-modal support for malignant-state clusters.",
  "",
  "## Outputs",
  "",
  "- `plots/Figure2_LAP3_State_Submodules.svg`",
  "- `plots/Figure2_LAP3_State_Submodules.pdf`",
  "- `plots/Figure2_LAP3_State_Submodules.tiff`",
  "- `plots/Figure2_LAP3_State_Submodules.png`",
  "- `source_data/figure2_panel_*.csv`",
  "- `tables/figure2_key_results.csv`",
  "- `tables/figure2_panel_map.csv`",
  "- `tables/figure2_export_qc.csv`",
  "",
  "## Interpretation Boundary",
  "",
  "Figure 2 should be written as an interpretability and robustness figure for a composite LAP3-centered transcriptional state. It should not be written as evidence that LAP3 enzymatic activity, intracellular leucine flux or causal mTORC1 phosphorylation has been proven."
)
writeLines(readme, file.path(out_dir, "README.md"))

if (!all(export_qc$exists) || any(is.na(export_qc$bytes) | export_qc$bytes <= 0)) {
  stop("Export QC failed: ", file.path(table_dir, "figure2_export_qc.csv"))
}

message("Figure 2 export completed: ", out_dir)
