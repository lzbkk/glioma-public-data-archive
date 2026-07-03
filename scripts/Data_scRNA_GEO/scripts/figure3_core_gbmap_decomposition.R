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
  file_arg <- "Data_scRNA_GEO/scripts/figure3_core_gbmap_decomposition.R"
}
repo <- normalizePath(file.path(dirname(file_arg), "../.."), mustWork = FALSE)
if (!dir.exists(file.path(repo, "Project_Management"))) repo <- normalizePath(".", mustWork = TRUE)

gbmap_dir <- file.path(repo, "Data_scRNA_GEO/GBmap_Core/results/LAP3_State_Union_Projection")
submodule_dir <- file.path(repo, "Data_Bulk_TCGA/Data_Merged/results/LAP3_State_Submodules")
malignant_dir <- file.path(repo, "Data_Bulk_TCGA/Data_Merged/results/LAP3_Malignant_State_Module_Audit")
out_dir <- file.path(repo, "Data_scRNA_GEO/GBmap_Core/results/Figure3_Core_GBmap_Decomposition")
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

entry_counts <- read_dt(file.path(gbmap_dir, "tables/gbmap_core_lap3_state_entry_counts.csv"))
coverage <- read_dt(file.path(gbmap_dir, "tables/gbmap_core_lap3_state_gene_coverage.csv"))
primary_summary <- read_dt(file.path(gbmap_dir, "tables/gbmap_core_lap3_state_primary_summary.csv"))
within_state <- read_dt(file.path(gbmap_dir, "tables/gbmap_core_lap3_state_within_state_lap3_associations.csv"))
preference <- read_dt(file.path(gbmap_dir, "tables/gbmap_core_lap3_state_preference_fixed_effect.csv"))
retention <- read_dt(file.path(gbmap_dir, "tables/gbmap_core_lap3_state_threshold_retention.csv"))
patient_state <- read_dt(file.path(gbmap_dir, "tables/gbmap_core_lap3_state_patient_state_summary.csv"))
submodule_primary <- read_dt(file.path(submodule_dir, "tables/gbmap_core_lap3_state_submodule_primary_summary.csv"))
malignant_gbmap <- read_dt(file.path(malignant_dir, "tables/gbmap_core_malignant_cluster_donor_state_associations.csv"))

state_levels <- c("AC", "MES", "NPC", "OPC")
state_cols <- c("AC" = "#315E8A", "MES" = "#B15A2C", "NPC" = "#7B5B8F", "OPC" = "#7A9B55")
set_labels <- c(
  "LAP3_STATE_UNION" = "Union",
  "LAP3_STATE_UNION_NO_TRANSLATION_PROTEOSTASIS" = "No T/P"
)
submodule_labels <- c(
  "LAP3_MALIGNANT_STATE_MODULE" = "Malignant-state",
  "LAP3_MYELOID_TAM_CONTEXT_MODULE" = "Myeloid/TAM",
  "LAP3_ANABOLIC_TRANSLATION_MODULE" = "Anabolic/translation",
  "LAP3_PROTEOSTASIS_STRESS_MODULE" = "Proteostasis/stress",
  "LAP3_HYPOXIA_PERINECROTIC_MODULE" = "Hypoxia/perinecrotic"
)
submodule_order <- names(submodule_labels)

panel_a <- copy(entry_counts)
panel_a[, entry_label := factor(
  entry_variant,
  levels = c("all_core", "exclude_neftel2019", "main_neoplastic_exclude_neftel2019", "strict_neoplastic_aneuploid_exclude_neftel2019"),
  labels = c("Core atlas", "Exclude Neftel", "Author-neoplastic", "Strict aneuploid")
)]
panel_a[, label := paste0(comma(n_cells), "\n", n_authors, " authors / ", n_donors, " donors")]
panel_a[, short_label := paste0(round(n_cells / 1000, 1), "k\n", n_authors, "A/", n_donors, "D")]

panel_b <- preference[
  entry_variant == "main_neoplastic_exclude_neftel2019" &
    threshold == 20 &
    state_set %chin% names(set_labels)
]
panel_b[, state := sub("^author_state", "", term)]
panel_b[, state := factor(state, levels = c("MES", "NPC", "OPC"))]
panel_b[, state_set_label := factor(set_labels[state_set], levels = set_labels)]
panel_b[, label := sprintf("%.2f", estimate)]

panel_c <- within_state[
  entry_variant == "main_neoplastic_exclude_neftel2019" &
    threshold == 20 &
    state_set %chin% names(set_labels)
]
panel_c[, author_state := factor(author_state, levels = state_levels)]
panel_c[, state_set_label := factor(set_labels[state_set], levels = set_labels)]
panel_c[, sig_label := fifelse(p_adj_BH < 0.05, "FDR < 0.05", "FDR >= 0.05")]

panel_d <- copy(retention[entry_variant == "main_neoplastic_exclude_neftel2019" & threshold %in% c(20, 50, 100)])
panel_d[, author_state := factor(author_state, levels = state_levels)]
panel_d[, threshold := factor(threshold, levels = c(20, 50, 100))]

panel_e <- submodule_primary[
  entry_variant == "main_neoplastic_exclude_neftel2019" &
    threshold == 20 &
    submodule %chin% submodule_order
]
panel_e[, author_state := factor(author_state, levels = state_levels)]
panel_e[, submodule_label := factor(submodule_labels[submodule], levels = rev(submodule_labels[submodule_order]))]
panel_e[, label := sprintf("%.2f", spearman_rho)]

panel_f <- malignant_gbmap[
  entry_variant == "main_neoplastic" &
    target == "lap3_mean" &
    cluster_score %chin% paste0("LAP3_MALIGNANT_M", 1:3)
]
panel_f[, author_state := factor(author_state, levels = state_levels)]
panel_f[, cluster_label := factor(sub("LAP3_MALIGNANT_", "", cluster_score), levels = rev(c("M1", "M2", "M3")))]
panel_f[, label := sprintf("%.2f", rho)]

panel_coverage <- copy(coverage)
panel_patient_state_summary <- patient_state[
  entry_variant == "main_neoplastic_exclude_neftel2019" & threshold == 20,
  .(
    n_donor_states = .N,
    n_donors = uniqueN(author_donor),
    median_cells = as.numeric(median(n_cells)),
    median_lap3_detection_rate = as.numeric(median(lap3_detection_rate)),
    median_union_score = as.numeric(median(LAP3_STATE_UNION))
  ),
  by = author_state
][order(match(author_state, state_levels))]

panel_map <- data.table(
  panel = LETTERS[1:6],
  source_data = c(
    "figure3_panel_a_core_gbmap_inference_frame.csv",
    "figure3_panel_b_state_preference_fixed_effect.csv",
    "figure3_panel_c_within_state_lap3_state_coupling.csv",
    "figure3_panel_d_threshold_retention.csv",
    "figure3_panel_e_submodule_decomposition.csv",
    "figure3_panel_f_malignant_cluster_decomposition.csv"
  ),
  conclusion = c(
    "Core GBmap provides a main author-neoplastic entry and strict aneuploid sensitivity without reusing Neftel2019.",
    "Donor fixed-effect models place the frozen union highest in AC-like reference state and lower in OPC/NPC/MES.",
    "Within-state donor-level LAP3-state coupling is strongest in AC/MES/OPC; NPC is weaker for union but supported by no-translation/proteostasis sensitivity.",
    "The 20-cell donor-state threshold retains all four malignant states and 50-cell sensitivity remains feasible.",
    "Submodule decomposition shows malignant-state and context modules contribute differently across AC/MES/NPC/OPC.",
    "The malignant-state arm decomposes into M1/M2/M3 with state-dependent LAP3 coupling."
  )
)

key_results <- rbindlist(list(
  panel_a[entry_variant == "main_neoplastic_exclude_neftel2019",
          .(panel = "A", feature = "Main author-neoplastic entry", statistic = "cells/authors/donors",
            value = n_cells, note = paste0(n_authors, " authors; ", n_donors, " donors"))],
  panel_coverage[, .(panel = "A", feature = state_set, statistic = "coverage", value = coverage, note = paste0(available, "/", requested, " genes"))],
  panel_b[, .(panel = "B", feature = paste(state, state_set_label, sep = " / "), statistic = "fixed-effect estimate vs AC", value = estimate, note = paste0("FDR=", signif(p_adj_BH, 3)))],
  panel_c[, .(panel = "C", feature = paste(author_state, state_set_label, sep = " / "), statistic = "Spearman rho with LAP3", value = spearman_rho, note = paste0("FDR=", signif(p_adj_BH, 3)))],
  panel_e[submodule == "LAP3_MALIGNANT_STATE_MODULE",
          .(panel = "E", feature = paste(author_state, "malignant-state submodule"), statistic = "Spearman rho with LAP3", value = spearman_rho, note = paste0("FDR=", signif(p_adj_BH, 3)))],
  panel_f[, .(panel = "F", feature = paste(cluster_label, author_state, sep = " / "), statistic = "Spearman rho with LAP3", value = rho, note = paste0("FDR=", signif(p_adj_BH, 3)))]
), fill = TRUE)

write_source(panel_a, "figure3_panel_a_core_gbmap_inference_frame.csv")
write_source(panel_b, "figure3_panel_b_state_preference_fixed_effect.csv")
write_source(panel_c, "figure3_panel_c_within_state_lap3_state_coupling.csv")
write_source(panel_d, "figure3_panel_d_threshold_retention.csv")
write_source(panel_e, "figure3_panel_e_submodule_decomposition.csv")
write_source(panel_f, "figure3_panel_f_malignant_cluster_decomposition.csv")
write_table(panel_map, "figure3_panel_map.csv")
write_table(key_results, "figure3_key_results.csv")
write_table(panel_coverage, "figure3_gene_coverage_input.csv")
write_table(panel_patient_state_summary, "figure3_patient_state_summary_input.csv")

p_a <- ggplot(panel_a, aes(x = entry_label, y = n_cells)) +
  geom_col(fill = "#315E8A", width = 0.62) +
  geom_text(aes(y = n_cells * 0.56, label = short_label), colour = "white",
            size = 2.05, lineheight = 0.86, fontface = "bold") +
  scale_y_continuous(labels = comma, limits = c(0, max(panel_a$n_cells) * 1.05), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "A  Core GBmap entry frame", x = NULL, y = "Cells") +
  theme_fig() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1))

p_b <- ggplot(panel_b, aes(x = estimate, y = state, colour = state_set_label)) +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "#8B9199") +
  geom_errorbarh(aes(xmin = estimate - 1.96 * std_error, xmax = estimate + 1.96 * std_error),
                 height = 0.14, linewidth = 0.35, position = position_dodge(width = 0.45)) +
  geom_point(size = 2.1, position = position_dodge(width = 0.45)) +
  scale_colour_manual(values = c("Union" = "#315E8A", "No T/P" = "#B15A2C"), name = NULL) +
  labs(title = "B  Donor fixed-effect preference", x = "State-score estimate vs AC", y = NULL) +
  theme_fig() +
  theme(legend.position = "bottom")

p_c <- ggplot(panel_c, aes(x = author_state, y = spearman_rho, colour = state_set_label, shape = sig_label)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "#8B9199") +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.12, linewidth = 0.35, position = position_dodge(width = 0.45)) +
  geom_point(size = 2.25, position = position_dodge(width = 0.45)) +
  scale_colour_manual(values = c("Union" = "#315E8A", "No T/P" = "#B15A2C"), name = NULL) +
  scale_shape_manual(values = c("FDR < 0.05" = 16, "FDR >= 0.05" = 1), name = NULL) +
  scale_y_continuous(limits = c(-0.25, 0.9), breaks = seq(-0.2, 0.8, 0.2)) +
  labs(title = "C  Within-state LAP3 coupling", x = "Cell state", y = "Spearman rho") +
  theme_fig() +
  guides(colour = "none") +
  theme(legend.position = "bottom")

p_d <- ggplot(panel_d, aes(x = author_state, y = n_donor_states, fill = threshold)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.62, colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = c("20" = "#315E8A", "50" = "#6E91B8", "100" = "#BAC8D8"), name = "cells/state") +
  labs(title = "D  Threshold retention", x = "Cell state", y = "Donor-states") +
  theme_fig() +
  theme(legend.position = "bottom")

p_e <- ggplot(panel_e, aes(x = author_state, y = submodule_label, fill = spearman_rho)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_text(aes(label = label), size = 1.85, colour = "#1E252B") +
  scale_fill_gradient2(low = "#D8DEE6", mid = "#F4F4F0", high = "#315E8A", midpoint = 0.25,
                       limits = c(-0.3, 0.8), oob = squish, name = "rho") +
  labs(title = "E  Submodule decomposition", x = "Cell state", y = NULL) +
  theme_fig() +
  theme(legend.position = "right")

p_f <- ggplot(panel_f, aes(x = author_state, y = cluster_label, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_text(aes(label = label), size = 1.95, colour = "#1E252B") +
  scale_fill_gradient2(low = "#D8DEE6", mid = "#F4F4F0", high = "#7B5B8F", midpoint = 0.25,
                       limits = c(-0.3, 1), oob = squish, name = "rho") +
  labs(title = "F  Malignant-cluster decomposition", x = "Cell state", y = NULL) +
  theme_fig() +
  theme(legend.position = "right")

caption <- paste(
  "Figure 3 projects the frozen LAP3-state into Core GBmap malignant cells and decomposes it by state and module.\n",
  "The analysis uses donor-state summaries as inference units and does not test causal malignant-cell-intrinsic biology."
)

figure3 <- (p_a | p_b) / (p_c | p_d) / (p_e | p_f) +
  plot_layout(heights = c(0.95, 1.05, 1.15), guides = "keep") +
  plot_annotation(
    title = "Figure 3. Core GBmap resolves the LAP3-state within malignant-cell states",
    caption = caption,
    theme = theme(
      plot.title = element_text(face = "bold", size = 10, family = "sans"),
      plot.caption = element_text(size = 6.5, colour = "#3A3D42", hjust = 0, family = "sans"),
      plot.margin = margin(5, 5, 5, 5)
    )
  )

plot_base <- file.path(plot_dir, "Figure3_Core_GBmap_Decomposition")
svglite::svglite(paste0(plot_base, ".svg"), width = 7.2, height = 7.6)
print(figure3)
dev.off()
ggsave(paste0(plot_base, ".pdf"), figure3, width = 7.2, height = 7.6, units = "in", device = grDevices::pdf)
ragg::agg_tiff(paste0(plot_base, ".tiff"), width = 7.2, height = 7.6, units = "in", res = 600, compression = "lzw")
print(figure3)
dev.off()
ragg::agg_png(paste0(plot_base, ".png"), width = 7.2, height = 7.6, units = "in", res = 300)
print(figure3)
dev.off()

export_files <- c(
  paste0(plot_base, c(".svg", ".pdf", ".tiff", ".png")),
  file.path(source_dir, panel_map$source_data),
  file.path(table_dir, c(
    "figure3_panel_map.csv",
    "figure3_key_results.csv",
    "figure3_gene_coverage_input.csv",
    "figure3_patient_state_summary_input.csv"
  ))
)
export_qc <- data.table(
  file = export_files,
  exists = file.exists(export_files),
  bytes = as.numeric(file.info(export_files)$size)
)
write_table(export_qc, "figure3_export_qc.csv")

readme <- c(
  "# Figure 3 Core GBmap Decomposition",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Purpose",
  "",
  "Create the active manuscript Figure 3 for the state/ecosystem route.",
  "",
  "## Core Conclusion",
  "",
  "The frozen `LAP3_STATE_UNION` projects into Core GBmap malignant cells and is resolved across AC/MES/NPC/OPC donor-state summaries. The signal is strongest in AC/MES/OPC for the union score, NPC is sensitivity-supported after removing translation/proteostasis genes, and submodule/cluster views show an interpretable malignant-state decomposition.",
  "",
  "## Inputs",
  "",
  "- `Data_scRNA_GEO/GBmap_Core/results/LAP3_State_Union_Projection/`",
  "- `Data_Bulk_TCGA/Data_Merged/results/LAP3_State_Submodules/`",
  "- `Data_Bulk_TCGA/Data_Merged/results/LAP3_Malignant_State_Module_Audit/`",
  "",
  "## Panels",
  "",
  "- Panel A: Core GBmap inference frame.",
  "- Panel B: donor fixed-effect state preference.",
  "- Panel C: within-state LAP3-state coupling.",
  "- Panel D: donor-state threshold retention.",
  "- Panel E: submodule decomposition.",
  "- Panel F: malignant-cluster decomposition.",
  "",
  "## Outputs",
  "",
  "- `plots/Figure3_Core_GBmap_Decomposition.svg`",
  "- `plots/Figure3_Core_GBmap_Decomposition.pdf`",
  "- `plots/Figure3_Core_GBmap_Decomposition.tiff`",
  "- `plots/Figure3_Core_GBmap_Decomposition.png`",
  "- `source_data/figure3_panel_*.csv`",
  "- `tables/figure3_key_results.csv`",
  "- `tables/figure3_panel_map.csv`",
  "- `tables/figure3_export_qc.csv`",
  "",
  "## Interpretation Boundary",
  "",
  "Figure 3 should be written as malignant-cell-state decomposition of a bulk-derived LAP3-centered state. It should not be written as proof of malignant-cell-intrinsic LAP3 causality, intracellular leucine flux, or direct mTORC1 activation."
)
writeLines(readme, file.path(out_dir, "README.md"))

if (!all(export_qc$exists) || any(is.na(export_qc$bytes) | export_qc$bytes <= 0)) {
  stop("Export QC failed: ", file.path(table_dir, "figure3_export_qc.csv"))
}

message("Figure 3 export completed: ", out_dir)
