suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

required_pkgs <- c("svglite", "ragg")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))
}

set.seed(20260701)
data.table::setDTthreads(8)

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)][1])
if (is.na(file_arg) || !nzchar(file_arg)) {
  file_arg <- "Data_Bulk_TCGA/Data_Merged/scripts/figure5_cptac_glass_boundary.R"
}
repo <- normalizePath(file.path(dirname(file_arg), "../../.."), mustWork = FALSE)
if (!dir.exists(file.path(repo, "Project_Management"))) {
  repo <- normalizePath(".", mustWork = TRUE)
}

input_dir <- file.path(repo, "Data_Bulk_TCGA/Data_Merged/results/LAP3_State_CPTAC_GLASS_Projection")
table_dir <- file.path(input_dir, "tables")
out_dir <- file.path(repo, "Data_Bulk_TCGA/Data_Merged/results/Figure5_CPTAC_GLASS_Boundary")
source_dir <- file.path(out_dir, "source_data")
plot_dir <- file.path(out_dir, "plots")
out_table_dir <- file.path(out_dir, "tables")
log_dir <- file.path(out_dir, "logs")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

read_dt <- function(name) {
  path <- file.path(table_dir, name)
  if (!file.exists(path)) stop("Missing input table: ", path)
  fread(path)
}

write_source <- function(dt, name) {
  path <- file.path(source_dir, name)
  fwrite(copy(dt), path)
  invisible(path)
}

write_table <- function(dt, name) {
  path <- file.path(out_table_dir, name)
  fwrite(copy(dt), path)
  invisible(path)
}

module_map <- data.table(
  module = c(
    "LAP3_STATE_UNION",
    "LAP3_MALIGNANT_STATE_MODULE",
    "LAP3_MYELOID_TAM_CONTEXT_MODULE",
    "LAP3_PROTEOSTASIS_STRESS_MODULE",
    "LAP3_ANABOLIC_TRANSLATION_MODULE",
    "LAP3_HYPOXIA_PERINECROTIC_MODULE"
  ),
  module_label = c(
    "State union",
    "Malignant state",
    "Myeloid/TAM context",
    "Proteostasis/stress",
    "Anabolic/translation",
    "Hypoxia/perinecrotic"
  ),
  module_order = 6:1
)

parse_exposure <- function(x) {
  assay <- fifelse(startsWith(x, "mrna_"), "mRNA",
    fifelse(startsWith(x, "protein_"), "Protein",
      fifelse(x == "lap3_mrna_log2", "LAP3 mRNA",
        fifelse(x == "lap3_protein", "LAP3 protein", "Other")
      )
    )
  )
  module <- gsub("^mrna_|^protein_", "", x)
  module <- fifelse(x == "lap3_mrna_log2", "LAP3_mRNA", module)
  module <- fifelse(x == "lap3_protein", "LAP3_protein", module)
  data.table(exposure = x, assay = assay, module = module)
}

label_exposures <- function(dt) {
  parsed <- parse_exposure(dt$exposure)
  out <- cbind(copy(dt), parsed[, .(assay, module)])
  out <- merge(out, module_map, by = "module", all.x = TRUE, sort = FALSE)
  out[module == "LAP3_mRNA", `:=`(module_label = "LAP3 mRNA", module_order = 7)]
  out[module == "LAP3_protein", `:=`(module_label = "LAP3 protein", module_order = 7)]
  out[, sig_label := fifelse(fdr < 0.05, "FDR < 0.05", "FDR >= 0.05")]
  out[]
}

fmt_p <- function(x) {
  ifelse(is.na(x), "NA", ifelse(x < 1e-3, scientific(x, digits = 2), number(x, accuracy = 0.001)))
}

cptac_cor <- read_dt("cptac_lap3_state_submodule_correlations.csv")
glass_tests <- read_dt("glass_lap3_state_submodule_paired_tests.csv")
coverage <- read_dt("projection_gene_coverage.csv")
summary_dt <- read_dt("lap3_state_cptac_glass_primary_summary.csv")

stopifnot(
  nrow(cptac_cor) > 0,
  nrow(glass_tests) > 0,
  all(c("dataset", "stratum", "exposure", "outcome", "n_complete", "spearman_rho", "p_value", "fdr") %in% names(cptac_cor)),
  all(c("dataset", "analysis", "stratum", "module", "n_pairs", "median_delta", "estimate", "p_value", "fdr") %in% names(glass_tests))
)

module_regex <- paste(module_map$module, collapse = "|")

panel_a <- data.table(
  panel = c("CPTAC", "CPTAC", "CPTAC", "GLASS"),
  evidence_layer = c("Protein bridge", "Direct BCAA boundary", "Phospho readout boundary", "Longitudinal coordination"),
  inference_unit = c("Tumor sample/patient", "Tumor sample/patient", "Tumor sample/patient", "Strict patient pair"),
  n = c(
    cptac_cor[stratum == "All" & outcome == "lap3_protein", max(n_complete, na.rm = TRUE)],
    cptac_cor[stratum == "All" & outcome == "bcaa_composite", max(n_complete, na.rm = TRUE)],
    cptac_cor[stratum == "All" & outcome == "phospho_mtorc1_target_score", max(n_complete, na.rm = TRUE)],
    glass_tests[analysis == "delta_lap3_delta_module" & stratum == "All", max(n_pairs, na.rm = TRUE)]
  ),
  result_role = c(
    "Concordance",
    "Negative boundary",
    "Cross-sectional readout",
    "Paired delta signal"
  )
)
panel_a[, y := rev(seq_len(.N))]
panel_a[, xmin := 0]
panel_a[, xmax := 1]
panel_a[, label_left := paste0(panel, "\n", evidence_layer)]
panel_a[, label_mid := paste0("n=", n, "\n", inference_unit)]
panel_a[, label_right := result_role]

panel_b <- cptac_cor[
  stratum == "All" &
    outcome == "lap3_protein" &
    (exposure == "lap3_mrna_log2" | grepl(module_regex, exposure))
]
panel_b <- label_exposures(panel_b)
panel_b[, y_label := factor(module_label, levels = module_map[order(module_order), module_label] |> c("LAP3 mRNA") |> rev())]
panel_b[, assay_plot := fifelse(assay %chin% c("mRNA", "Protein"), assay, "LAP3 mRNA")]

panel_c <- cptac_cor[
  stratum == "All" &
    outcome %chin% c("bcaa_composite", "phospho_mtorc1_target_score") &
    (exposure %chin% c("lap3_mrna_log2", "lap3_protein") | grepl(module_regex, exposure))
]
panel_c <- label_exposures(panel_c)
panel_c[, outcome_label := fifelse(outcome == "bcaa_composite", "Direct BCAA composite", "mTORC1 phospho readout")]
panel_c[, y_label := factor(module_label, levels = c("LAP3 protein", "LAP3 mRNA", module_map[order(module_order), module_label]))]

panel_d <- glass_tests[
  analysis == "delta_lap3_delta_module" &
    stratum %chin% c("All", "IDH-wildtype GBM") &
    module %chin% module_map$module
]
panel_d <- merge(panel_d, module_map, by = "module", all.x = TRUE, sort = FALSE)
panel_d[, y_label := factor(module_label, levels = module_map[order(module_order), module_label])]
panel_d[, sig_label := fifelse(fdr < 0.05, "FDR < 0.05", "FDR >= 0.05")]

panel_e <- glass_tests[
  analysis == "paired_recurrence_change" &
    stratum %chin% c("All", "IDH-wildtype GBM") &
    module %chin% module_map$module
]
panel_e <- merge(panel_e, module_map, by = "module", all.x = TRUE, sort = FALSE)
panel_e[, y_label := factor(module_label, levels = module_map[order(module_order), module_label])]
panel_e[, sig_label := fifelse(fdr < 0.05, "FDR < 0.05", "FDR >= 0.05")]

key_results <- rbindlist(list(
  panel_b[, .(
    panel = "B",
    dataset, stratum, feature = exposure, outcome,
    n = n_complete, effect = spearman_rho, p_value, fdr,
    interpretation = "CPTAC LAP3 protein bridge"
  )],
  panel_c[, .(
    panel = "C",
    dataset, stratum, feature = exposure, outcome,
    n = n_complete, effect = spearman_rho, p_value, fdr,
    interpretation = fifelse(outcome == "bcaa_composite", "direct BCAA negative boundary", "cross-sectional phospho readout boundary")
  )],
  panel_d[, .(
    panel = "D",
    dataset, stratum, feature = module, outcome = analysis,
    n = n_pairs, effect = estimate, p_value, fdr,
    interpretation = "GLASS paired delta coordination"
  )],
  panel_e[, .(
    panel = "E",
    dataset, stratum, feature = module, outcome = analysis,
    n = n_pairs, effect = median_delta, p_value, fdr,
    interpretation = "GLASS recurrence change boundary"
  )]
), fill = TRUE)

panel_map <- data.table(
  panel = c("A", "B", "C", "D", "E"),
  source_data = c(
    "figure5_panel_a_cohort_frame.csv",
    "figure5_panel_b_cptac_lap3_protein_bridge.csv",
    "figure5_panel_c_cptac_bcaa_phospho_boundary.csv",
    "figure5_panel_d_glass_delta_coordination.csv",
    "figure5_panel_e_glass_recurrence_boundary.csv"
  ),
  conclusion = c(
    "CPTAC and GLASS layers are separated by inference unit and evidentiary role.",
    "LAP3 protein tracks LAP3 mRNA and state/submodule protein-mRNA scores.",
    "Direct BCAA remains null whereas phospho mTORC1 target score is only cross-sectional.",
    "Patient-level delta LAP3 is strongly coordinated with delta state modules.",
    "There is no uniform recurrence activation across all paired tumors."
  )
)

write_source(panel_a, "figure5_panel_a_cohort_frame.csv")
write_source(panel_b, "figure5_panel_b_cptac_lap3_protein_bridge.csv")
write_source(panel_c, "figure5_panel_c_cptac_bcaa_phospho_boundary.csv")
write_source(panel_d, "figure5_panel_d_glass_delta_coordination.csv")
write_source(panel_e, "figure5_panel_e_glass_recurrence_boundary.csv")
write_table(key_results, "figure5_key_results.csv")
write_table(panel_map, "figure5_panel_map.csv")
write_table(coverage, "figure5_projection_gene_coverage.csv")
write_table(summary_dt, "figure5_input_primary_summary.csv")

theme_fig <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 1, hjust = 0),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.5, colour = "#222222"),
      strip.background = element_rect(fill = "#F2F3F5", colour = NA),
      strip.text = element_text(face = "bold", size = base_size),
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size - 0.5),
      legend.key.height = unit(3.5, "mm"),
      plot.margin = margin(4, 4, 4, 4)
    )
}

col_layer <- c(
  "Protein bridge" = "#2F6F9F",
  "Direct BCAA boundary" = "#747C87",
  "Phospho readout boundary" = "#B15A2C",
  "Longitudinal coordination" = "#2F7D5B"
)
col_sig <- c("FDR < 0.05" = "#1F6F8B", "FDR >= 0.05" = "#B7BDC5")
col_stratum <- c("All" = "#38598A", "IDH-wildtype GBM" = "#C05A37")

p_a <- ggplot(panel_a) +
  geom_rect(aes(xmin = 0, xmax = 1, ymin = y - 0.42, ymax = y + 0.42, fill = evidence_layer),
            colour = "white", linewidth = 0.35) +
  geom_text(aes(x = 0.03, y = y, label = label_left), hjust = 0, colour = "white",
            size = 2.2, lineheight = 0.9, fontface = "bold") +
  geom_text(aes(x = 0.48, y = y, label = label_mid), hjust = 0.5, colour = "white",
            size = 2.05, lineheight = 0.9) +
  geom_text(aes(x = 0.97, y = y, label = label_right), hjust = 1, colour = "white",
            size = 2.05, lineheight = 0.9) +
  scale_fill_manual(values = col_layer, guide = "none") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0.45, 4.55), expand = FALSE) +
  labs(title = "A  Multi-omic evidence frame") +
  theme_void(base_family = "sans") +
  theme(
    plot.title = element_text(face = "bold", size = 8, hjust = 0),
    plot.margin = margin(4, 4, 4, 4)
  )

p_b <- ggplot(panel_b, aes(x = spearman_rho, y = y_label)) +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "#7C828A") +
  geom_segment(aes(x = 0, xend = spearman_rho, yend = y_label, colour = sig_label),
               linewidth = 0.45, alpha = 0.7) +
  geom_point(aes(shape = assay_plot, colour = sig_label), size = 2.3, stroke = 0.4) +
  scale_colour_manual(values = col_sig, name = NULL) +
  scale_shape_manual(values = c("mRNA" = 16, "Protein" = 17, "LAP3 mRNA" = 15), name = NULL) +
  scale_x_continuous(limits = c(-0.08, 0.82), breaks = seq(0, 0.8, 0.2)) +
  labs(
    title = "B  CPTAC LAP3 protein bridge",
    x = "Spearman rho with LAP3 protein",
    y = NULL
  ) +
  theme_fig()

p_c <- ggplot(panel_c, aes(x = spearman_rho, y = y_label)) +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "#7C828A") +
  geom_segment(aes(x = 0, xend = spearman_rho, yend = y_label, colour = sig_label),
               linewidth = 0.4, alpha = 0.75) +
  geom_point(aes(shape = assay, colour = sig_label), size = 1.9, stroke = 0.4) +
  facet_wrap(~ outcome_label, nrow = 1) +
  scale_colour_manual(values = col_sig, name = NULL) +
  scale_shape_manual(values = c("mRNA" = 16, "Protein" = 17, "LAP3 mRNA" = 15, "LAP3 protein" = 18), name = NULL) +
  scale_x_continuous(limits = c(-0.16, 0.42), breaks = seq(-0.1, 0.4, 0.1)) +
  labs(
    title = "C  CPTAC direct metabolite and phospho-readout boundaries",
    x = "Spearman rho",
    y = NULL
  ) +
  theme_fig()

p_d <- ggplot(panel_d, aes(x = estimate, y = y_label, colour = stratum)) +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "#7C828A") +
  geom_segment(aes(x = 0, xend = estimate, yend = y_label), linewidth = 0.45, alpha = 0.55,
               position = position_dodge(width = 0.45)) +
  geom_point(aes(shape = sig_label), position = position_dodge(width = 0.45), size = 2.2, stroke = 0.4) +
  scale_colour_manual(values = col_stratum, name = NULL) +
  scale_shape_manual(values = c("FDR < 0.05" = 16, "FDR >= 0.05" = 1), name = NULL) +
  scale_x_continuous(limits = c(0, 0.9), breaks = seq(0, 0.8, 0.2)) +
  labs(
    title = "D  GLASS paired delta coordination",
    x = "Spearman rho: delta LAP3 vs delta module",
    y = NULL
  ) +
  theme_fig()

p_e <- ggplot(panel_e, aes(x = median_delta, y = y_label, colour = stratum)) +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "#7C828A") +
  geom_segment(aes(x = 0, xend = median_delta, yend = y_label), linewidth = 0.45, alpha = 0.55,
               position = position_dodge(width = 0.45)) +
  geom_point(aes(shape = sig_label), position = position_dodge(width = 0.45), size = 2.2, stroke = 0.4) +
  scale_colour_manual(values = col_stratum, name = NULL) +
  scale_shape_manual(values = c("FDR < 0.05" = 16, "FDR >= 0.05" = 1), name = NULL) +
  scale_x_continuous(limits = c(-0.23, 0.11), breaks = seq(-0.2, 0.1, 0.1)) +
  labs(
    title = "E  GLASS recurrence-change boundary",
    x = "Median recurrence minus primary score",
    y = NULL
  ) +
  theme_fig()

caption_text <- paste(
  "CPTAC analyses are cross-sectional tumor-sample correlations; GLASS analyses use strict primary-to-recurrence patient pairs.\n",
  "Direct BCAA metabolites are negative; phospho and longitudinal signals are interpreted as state-linked associations."
)

figure5 <- (p_a / p_b / p_c / (p_d | p_e)) +
  plot_layout(heights = c(0.95, 1.45, 1.65, 1.55)) +
  plot_annotation(
    title = "Figure 5. CPTAC/GLASS evidence boundaries for the LAP3-state program",
    caption = caption_text,
    theme = theme(
      plot.title = element_text(face = "bold", size = 10, family = "sans"),
      plot.caption = element_text(size = 6.5, colour = "#3A3D42", hjust = 0, family = "sans"),
      plot.margin = margin(5, 5, 5, 5)
    )
  )

plot_base <- file.path(plot_dir, "Figure5_CPTAC_GLASS_Boundary")
svglite::svglite(paste0(plot_base, ".svg"), width = 7.2, height = 8.4)
print(figure5)
dev.off()
ggsave(paste0(plot_base, ".pdf"), figure5, width = 7.2, height = 8.4, units = "in", device = grDevices::pdf)
ragg::agg_tiff(paste0(plot_base, ".tiff"), width = 7.2, height = 8.4, units = "in", res = 600, compression = "lzw")
print(figure5)
dev.off()
ragg::agg_png(paste0(plot_base, ".png"), width = 7.2, height = 8.4, units = "in", res = 300)
print(figure5)
dev.off()

export_files <- c(
  paste0(plot_base, ".svg"),
  paste0(plot_base, ".pdf"),
  paste0(plot_base, ".tiff"),
  paste0(plot_base, ".png"),
  file.path(source_dir, panel_map$source_data),
  file.path(out_table_dir, c("figure5_key_results.csv", "figure5_panel_map.csv", "figure5_projection_gene_coverage.csv", "figure5_input_primary_summary.csv"))
)
export_qc <- data.table(
  file = export_files,
  exists = file.exists(export_files),
  bytes = as.numeric(file.info(export_files)$size)
)
write_table(export_qc, "figure5_export_qc.csv")

readme <- c(
  "# Figure 5 CPTAC/GLASS Boundary",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Purpose",
  "",
  "Create a manuscript-grade compact Figure 5 for the active LAP3-state/ecosystem route.",
  "",
  "## Input",
  "",
  "- `Data_Bulk_TCGA/Data_Merged/results/LAP3_State_CPTAC_GLASS_Projection/`",
  "- Frozen `LAP3_STATE_UNION` and five submodules from the LAP3-state submodule decomposition.",
  "- CPTAC GBM mRNA/protein/direct BCAA/phospho readout tables and GLASS strict primary-to-recurrence pairs.",
  "",
  "## Panels",
  "",
  "- Panel A: cohort and inference-unit frame.",
  "- Panel B: CPTAC LAP3 protein bridge.",
  "- Panel C: CPTAC direct BCAA negative boundary and phospho-readout cross-sectional boundary.",
  "- Panel D: GLASS paired delta coordination.",
  "- Panel E: GLASS recurrence-change boundary.",
  "",
  "## Key Interpretation",
  "",
  "- CPTAC supports LAP3 mRNA-protein concordance and LAP3-state protein bridge.",
  "- CPTAC direct BCAA metabolites remain a negative boundary.",
  "- CPTAC phospho mTORC1 target score is cross-sectional readout evidence, not causal LAP3 phosphorylation evidence.",
  "- GLASS supports patient-level coordinated state change, but not universal recurrence activation.",
  "",
  "## Outputs",
  "",
  "- `plots/Figure5_CPTAC_GLASS_Boundary.svg`",
  "- `plots/Figure5_CPTAC_GLASS_Boundary.pdf`",
  "- `plots/Figure5_CPTAC_GLASS_Boundary.tiff`",
  "- `plots/Figure5_CPTAC_GLASS_Boundary.png`",
  "- `source_data/figure5_panel_*.csv`",
  "- `tables/figure5_key_results.csv`",
  "- `tables/figure5_panel_map.csv`",
  "- `tables/figure5_export_qc.csv`",
  "",
  "## Method Boundary",
  "",
  "CPTAC analyses are cross-sectional tumor-sample correlations. GLASS analyses use strict patient pairs. The figure should be written as evidence-boundary support for a LAP3-centered state program, not as proof of LAP3-BCAA-mTORC1 causality."
)
writeLines(readme, file.path(out_dir, "README.md"))

if (!all(export_qc$exists) || any(is.na(export_qc$bytes) | export_qc$bytes <= 0)) {
  stop("Export QC failed. See ", file.path(out_table_dir, "figure5_export_qc.csv"))
}

message("Figure 5 export completed: ", out_dir)
