#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(ragg)
})

options(stringsAsFactors = FALSE)
set.seed(20260701)

project_root <- normalizePath(file.path(getwd()), mustWork = TRUE)

topology_dir <- file.path(
  project_root,
  "Data_Spatial_Public/GBM_Space/results/LAP3_State_Spatial_Topology"
)
module_dir <- file.path(
  project_root,
  "Data_Bulk_TCGA/Data_Merged/results/LAP3_Malignant_State_Module_Audit"
)
out_dir <- file.path(
  project_root,
  "Data_Spatial_Public/GBM_Space/results/Figure4_LAP3_State_Topology"
)
source_dir <- file.path(out_dir, "source_data")
plot_dir <- file.path(out_dir, "plots")
table_dir <- file.path(out_dir, "tables")
log_dir <- file.path(out_dir, "logs")

dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(
  log_dir,
  sprintf("figure4_lap3_state_topology.%s.log", format(Sys.time(), "%Y%m%d_%H%M%S"))
)
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)

message("START Figure 4 LAP3 state topology reconstruction")
message("Project root: ", project_root)
message("Output dir: ", out_dir)

need_files <- c(
  topology_summary = file.path(topology_dir, "tables/gbmspace_lap3_state_topology_summary.tsv"),
  contrast_summary = file.path(topology_dir, "tables/gbmspace_lap3_state_contrast_summary.tsv"),
  tumor_effects = file.path(topology_dir, "tables/gbmspace_lap3_state_tumor_effects.tsv"),
  loto = file.path(topology_dir, "tables/gbmspace_lap3_state_leave_one_tumor_out.tsv"),
  malignant_spatial = file.path(module_dir, "tables/gbmspace_malignant_cluster_spatial_summary.tsv")
)
missing_files <- need_files[!file.exists(need_files)]
if (length(missing_files) > 0) {
  stop("Missing required input files:\n", paste(missing_files, collapse = "\n"))
}

topology <- fread(need_files[["topology_summary"]])
contrast <- fread(need_files[["contrast_summary"]])
tumor_effects <- fread(need_files[["tumor_effects"]])
loto <- fread(need_files[["loto"]])
malignant_spatial <- fread(need_files[["malignant_spatial"]])

required_topology_cols <- c(
  "state_set", "readout", "readout_class", "priority", "analysis_type",
  "n_tumors", "median_raw_rho", "median_depth_adjusted_rho",
  "fdr_depth_adjusted_all", "fdr_depth_adjusted_priority"
)
required_contrast_cols <- c(
  "state_set", "readout", "readout_class", "priority", "contrast_type",
  "n_tumors", "median_tumor_delta", "fdr_delta_all", "fdr_delta_priority"
)
required_tumor_cols <- c(
  "tumor_id", "state_set", "readout", "readout_class", "priority",
  "analysis_type", "n_sections", "median_depth_adjusted_rho"
)
required_loto_cols <- c(
  "state_set", "analysis_type", "readout", "dropped_tumor_id",
  "n_tumors", "median_depth_adjusted_rho"
)
required_module_cols <- c(
  "target", "cluster_score", "n_sections", "n_tumors",
  "median_raw_rho", "median_depth_adjusted_rho", "p_adj_BH"
)

check_cols <- function(x, req, label) {
  miss <- setdiff(req, names(x))
  if (length(miss) > 0) {
    stop(label, " missing columns: ", paste(miss, collapse = ", "))
  }
}
check_cols(topology, required_topology_cols, "topology")
check_cols(contrast, required_contrast_cols, "contrast")
check_cols(tumor_effects, required_tumor_cols, "tumor_effects")
check_cols(loto, required_loto_cols, "loto")
check_cols(malignant_spatial, required_module_cols, "malignant_spatial")

class_palette <- c(
  myeloid_tam = "#C44E52",
  malignant_state = "#009E73",
  spatial_niche = "#0072B2",
  pathway = "#7E57C2",
  histopath = "#7A7A7A",
  other = "#999999"
)
class_labels <- c(
  myeloid_tam = "TAM/myeloid",
  malignant_state = "Malignant state",
  spatial_niche = "Spatial niche",
  pathway = "Pathway",
  histopath = "Histopathology",
  other = "Other"
)
target_labels <- c(
  LAP3_STATE_UNION = "LAP3 state union",
  LAP3_log1p_cp10k = "LAP3 expression"
)
cluster_labels <- c(
  LAP3_MALIGNANT_STATE_MODULE = "144-gene module",
  LAP3_MALIGNANT_M1 = "M1 dominant",
  LAP3_MALIGNANT_M2 = "M2 AC/stress",
  LAP3_MALIGNANT_M3 = "M3 boundary"
)

pretty_readout <- function(x) {
  map <- c(
    HALLMARK_MTORC1_SIGNALING = "mTORC1 hallmark",
    LEUCINE_BCAA_CORE = "BCAA core",
    MTORC1_READOUT_CORE = "mTORC1 readout",
    REACTOME_TRANSLATION = "Reactome translation",
    "RTN1..TAMs" = "RTN1/TAMs",
    "Immune..TAMs." = "Immune/TAMs",
    "Immune..resident." = "Immune/resident",
    "Dev.like..AC." = "Dev-like AC",
    "Dev.like..OPC." = "Dev-like OPC",
    "Gliosis..diffuse." = "Gliosis diffuse",
    "Leading.edge..white.matter." = "Leading edge/white matter",
    "Infiltrating.tumor..white.matter." = "Infiltrating tumor/white matter",
    "Infiltrating.tumor..grey.matter." = "Infiltrating tumor/grey matter"
  )
  out <- unname(map[x])
  idx <- is.na(out)
  out[idx] <- x[idx]
  out <- gsub("\\.\\.", "/", out)
  out <- gsub("_", " ", out)
  out <- gsub("\\.", " ", out)
  out <- gsub("\\s+", " ", out)
  trimws(out)
}

sig_label <- function(q) {
  fifelse(
    is.na(q), "not tested",
    fifelse(q < 0.001, "FDR<0.001",
      fifelse(q < 0.01, "FDR<0.01",
        fifelse(q < 0.05, "FDR<0.05", "n.s.")
      )
    )
  )
}

base_theme <- theme_classic(base_family = "sans", base_size = 8.5) +
  theme(
    plot.title = element_text(face = "bold", size = 9.5, margin = margin(b = 4)),
    plot.subtitle = element_text(size = 7.5, color = "#4D4D4D", margin = margin(b = 5)),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7.2, color = "#2B2B2B"),
    legend.title = element_text(size = 7.5),
    legend.text = element_text(size = 7),
    legend.key.height = unit(3.5, "mm"),
    strip.background = element_rect(fill = "#F3F4F6", color = NA),
    strip.text = element_text(face = "bold", size = 7.2, color = "#333333"),
    plot.margin = margin(6, 6, 6, 6)
  )

cohort_summary <- data.table(
  metric = c(
    "GBM-Space tumors",
    "Spatial sections",
    "Primary readout",
    "Inference unit",
    "Adjustment"
  ),
  value = c(
    as.character(unique(topology[state_set == "LAP3_STATE_UNION" &
      analysis_type == "neighbor_topology_k6" &
      readout == "Proliferative.TAMs", n_tumors])),
    as.character(unique(malignant_spatial[target == "LAP3_STATE_UNION" &
      cluster_score == "LAP3_MALIGNANT_STATE_MODULE", n_sections])),
    "Union",
    "Tumor",
    "Depth"
  ),
  note = c(
    "leave-one-tumor-out sensitivity",
    "spatial topology sections",
    "fixed 207-gene state signature",
    "tumor-aware summaries",
    "not causal neighborhood inference"
  )
)
cohort_summary[, y_metric := factor(metric, levels = rev(metric))]
cohort_summary[, x_tile := 0.55]
fwrite(cohort_summary, file.path(source_dir, "figure4_panel_a_cohort_summary.csv"))

main_readouts <- c(
  "Proliferative.TAMs", "RTN1..TAMs", "Monocytes", "Resident.BAM.TAMs",
  "Anti.inflammatory.TAMs", "Dendritic.cells",
  "AC.gliosis.like.1", "AC.gliosis.like.2", "AC.gliosis.like.3",
  "Proliferative.nIPC.like",
  "Gliosis.like", "Gliosis.transition", "Gliosis", "Vasculature",
  "HALLMARK_MTORC1_SIGNALING", "LEUCINE_BCAA_CORE",
  "MTORC1_READOUT_CORE", "REACTOME_TRANSLATION"
)

topology_plot_data <- topology[
  state_set == "LAP3_STATE_UNION" &
    analysis_type == "neighbor_topology_k6" &
    readout %in% main_readouts
]
topology_plot_data[, readout_label := pretty_readout(readout)]
topology_plot_data[, readout_class_label := class_labels[readout_class]]
topology_plot_data[, fdr_label := sig_label(fdr_depth_adjusted_priority)]
setorder(topology_plot_data, readout_class, median_depth_adjusted_rho)
topology_plot_data[, readout_label := factor(readout_label, levels = unique(readout_label))]
fwrite(topology_plot_data, file.path(source_dir, "figure4_panel_b_neighbor_topology.csv"))

contrast_readouts <- c(
  "Gliosis.transition", "Dev.like..AC.", "Dev.like..OPC.",
  "AC.gliosis.like.1", "Gliosis..diffuse.", "Immune..TAMs.",
  "Vasculature", "Resident.TAMs", "Proliferative",
  "HALLMARK_MTORC1_SIGNALING", "LEUCINE_BCAA_CORE", "MTORC1_READOUT_CORE",
  "REACTOME_TRANSLATION"
)
contrast_plot_data <- contrast[
  state_set == "LAP3_STATE_UNION" &
    contrast_type == "neighbor_top25_vs_bottom25_k6" &
    readout %in% contrast_readouts
]
contrast_plot_data[, readout_label := pretty_readout(readout)]
contrast_plot_data[, readout_class_label := class_labels[readout_class]]
contrast_plot_data[, fdr_label := sig_label(fdr_delta_priority)]
setorder(contrast_plot_data, readout_class, median_tumor_delta)
contrast_plot_data[, readout_label := factor(readout_label, levels = unique(readout_label))]
fwrite(contrast_plot_data, file.path(source_dir, "figure4_panel_c_neighbor_contrast.csv"))

tumor_readouts <- c(
  "Proliferative.TAMs", "RTN1..TAMs", "AC.gliosis.like.1",
  "AC.gliosis.like.3", "Gliosis.like", "Gliosis.transition",
  "HALLMARK_MTORC1_SIGNALING", "LEUCINE_BCAA_CORE"
)
tumor_plot_data <- tumor_effects[
  state_set == "LAP3_STATE_UNION" &
    analysis_type == "neighbor_topology_k6" &
    readout %in% tumor_readouts
]
tumor_plot_data[, readout_label := pretty_readout(readout)]
tumor_plot_data[, readout_class_label := class_labels[readout_class]]
setorder(tumor_plot_data, readout_class, median_depth_adjusted_rho)
tumor_levels <- tumor_plot_data[
  ,
  .(median_rho = median(median_depth_adjusted_rho, na.rm = TRUE)),
  by = .(readout, readout_label, readout_class)
][order(readout_class, median_rho), readout_label]
tumor_plot_data[, readout_label := factor(readout_label, levels = unique(tumor_levels))]
fwrite(tumor_plot_data, file.path(source_dir, "figure4_panel_d_tumor_effects.csv"))

loto_plot_data <- loto[
  state_set == "LAP3_STATE_UNION" &
    analysis_type == "neighbor_topology_k6" &
    readout %in% tumor_readouts
]
loto_summary <- loto_plot_data[
  ,
  .(
    n_leave_one_out = .N,
    median_loto_rho = median(median_depth_adjusted_rho, na.rm = TRUE),
    min_loto_rho = min(median_depth_adjusted_rho, na.rm = TRUE),
    max_loto_rho = max(median_depth_adjusted_rho, na.rm = TRUE),
    n_positive_loto = sum(median_depth_adjusted_rho > 0, na.rm = TRUE)
  ),
  by = .(readout)
]
loto_summary <- merge(
  loto_summary,
  unique(topology_plot_data[, .(readout, readout_class, readout_class_label, n_tumors)]),
  by = "readout",
  all.x = TRUE
)
loto_summary[, readout_label := pretty_readout(readout)]
setorder(loto_summary, readout_class, median_loto_rho)
loto_summary[, readout_label := factor(readout_label, levels = unique(readout_label))]
fwrite(loto_summary, file.path(source_dir, "figure4_panel_e_leave_one_tumor_out.csv"))

module_plot_data <- copy(malignant_spatial)
module_plot_data[, target_label := target_labels[target]]
module_plot_data[, cluster_label := cluster_labels[cluster_score]]
module_plot_data[, cluster_label := factor(
  cluster_label,
  levels = c("M3 boundary", "M2 AC/stress", "M1 dominant", "144-gene module")
)]
module_plot_data[, target_label := factor(
  target_label,
  levels = c("LAP3 expression", "LAP3 state union")
)]
fwrite(module_plot_data, file.path(source_dir, "figure4_panel_f_malignant_module_projection.csv"))

panel_map <- data.table(
  panel = LETTERS[1:6],
  source_data = c(
    "source_data/figure4_panel_a_cohort_summary.csv",
    "source_data/figure4_panel_b_neighbor_topology.csv",
    "source_data/figure4_panel_c_neighbor_contrast.csv",
    "source_data/figure4_panel_d_tumor_effects.csv",
    "source_data/figure4_panel_e_leave_one_tumor_out.csv",
    "source_data/figure4_panel_f_malignant_module_projection.csv"
  ),
  conclusion = c(
    "GBM-Space analysis uses tumor-aware, depth-adjusted topology summaries.",
    "LAP3_STATE_UNION is spatially enriched near TAM/myeloid and gliosis/malignant-state readouts.",
    "Top-neighborhood contrasts reproduce the gliosis and malignant-state enrichment pattern.",
    "Tumor-level effects show whether topology signals are broadly distributed across tumors.",
    "Leave-one-tumor-out intervals test whether selected signals depend on a single tumor.",
    "The 144-gene malignant-state module and M1/M2/M3 clusters project back onto GBM-Space topology."
  )
)
fwrite(panel_map, file.path(table_dir, "figure4_panel_map.csv"))

key_results <- rbindlist(list(
  topology_plot_data[
    ,
    .(
      source = "neighbor_topology",
      readout,
      readout_class,
      n_tumors,
      effect = median_depth_adjusted_rho,
      fdr = fdr_depth_adjusted_priority
    )
  ],
  contrast_plot_data[
    ,
    .(
      source = "neighbor_contrast",
      readout,
      readout_class,
      n_tumors,
      effect = median_tumor_delta,
      fdr = fdr_delta_priority
    )
  ],
  module_plot_data[
    target == "LAP3_STATE_UNION",
    .(
      source = "malignant_module_projection",
      readout = cluster_score,
      readout_class = "malignant_module",
      n_tumors,
      effect = median_depth_adjusted_rho,
      fdr = p_adj_BH
    )
  ]
), fill = TRUE)
fwrite(key_results, file.path(table_dir, "figure4_key_results.csv"))

p_a <- ggplot(cohort_summary, aes(y = y_metric)) +
  geom_tile(aes(x = x_tile, fill = metric), width = 0.62, height = 0.78, color = "white", linewidth = 0.4) +
  geom_text(aes(x = x_tile, label = value), fontface = "bold", size = 3.2, color = "white") +
  geom_text(aes(x = 0.94, label = metric), hjust = 0, fontface = "bold", size = 2.7, color = "#252525") +
  geom_text(aes(x = 0.94, label = note), hjust = 0, nudge_y = -0.19, size = 2.25, color = "#555555") +
  scale_fill_manual(values = c(
    "GBM-Space tumors" = "#2A9D8F",
    "Spatial sections" = "#457B9D",
    "Primary readout" = "#E76F51",
    "Inference unit" = "#8D6E63",
    "Adjustment" = "#6A4C93"
  ), guide = "none") +
  coord_cartesian(xlim = c(0.15, 2.45), clip = "off") +
  labs(title = "Cohort and inference frame") +
  theme_void(base_family = "sans", base_size = 8.5) +
  theme(
    plot.title = element_text(face = "bold", size = 9.5, margin = margin(b = 4)),
    plot.margin = margin(6, 18, 6, 6)
  )

p_b <- ggplot(
  topology_plot_data,
  aes(x = median_depth_adjusted_rho, y = readout_label, color = readout_class)
) +
  geom_vline(xintercept = 0, linewidth = 0.25, color = "#B0B0B0") +
  geom_segment(aes(x = 0, xend = median_depth_adjusted_rho, yend = readout_label),
    linewidth = 0.45, alpha = 0.65
  ) +
  geom_point(size = 3.4, shape = 17, stroke = 0.35) +
  scale_color_manual(values = class_palette, labels = class_labels, name = NULL) +
  labs(
    title = "Depth-adjusted spatial-neighborhood topology",
    subtitle = "LAP3_STATE_UNION versus selected GBM-Space readouts",
    x = "Median depth-adjusted rho",
    y = NULL
  ) +
  base_theme +
  theme(legend.position = "right")

p_c <- ggplot(
  contrast_plot_data,
  aes(x = median_tumor_delta, y = readout_label, color = readout_class)
) +
  geom_vline(xintercept = 0, linewidth = 0.25, color = "#B0B0B0") +
  geom_segment(aes(x = 0, xend = median_tumor_delta, yend = readout_label),
    linewidth = 0.45, alpha = 0.65
  ) +
  geom_point(size = 3.0, shape = 16, stroke = 0.35) +
  scale_color_manual(values = class_palette, labels = class_labels, name = NULL) +
  labs(
    title = "Top-neighborhood contrast",
    subtitle = "Top 25% versus bottom 25% LAP3_STATE_UNION neighborhood score",
    x = "Median tumor delta",
    y = NULL
  ) +
  base_theme +
  theme(legend.position = "none")

p_d <- ggplot(
  tumor_plot_data,
  aes(x = median_depth_adjusted_rho, y = readout_label, color = readout_class)
) +
  geom_vline(xintercept = 0, linewidth = 0.25, color = "#B0B0B0") +
  geom_point(
    position = position_jitter(height = 0.11, width = 0),
    alpha = 0.72,
    size = 1.65
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    shape = 95,
    size = 7.5,
    color = "#111111"
  ) +
  scale_color_manual(values = class_palette, labels = class_labels, name = NULL) +
  labs(
    title = "Tumor-level distribution",
    subtitle = "Each point is one tumor-level median effect",
    x = "Depth-adjusted rho",
    y = NULL
  ) +
  base_theme +
  theme(legend.position = "none")

p_e <- ggplot(
  loto_summary,
  aes(y = readout_label, color = readout_class)
) +
  geom_vline(xintercept = 0, linewidth = 0.25, color = "#B0B0B0") +
  geom_errorbarh(aes(xmin = min_loto_rho, xmax = max_loto_rho), height = 0.18, linewidth = 0.55) +
  geom_point(aes(x = median_loto_rho), shape = 16, size = 3.0) +
  scale_color_manual(values = class_palette, labels = class_labels, name = NULL) +
  labs(
    title = "Leave-one-tumor-out robustness",
    subtitle = "Intervals show min-to-max effect after dropping each tumor",
    x = "Median depth-adjusted rho after one-tumor drop",
    y = NULL
  ) +
  base_theme +
  theme(legend.position = "none")

p_f <- ggplot(
  module_plot_data,
  aes(x = median_depth_adjusted_rho, y = cluster_label, shape = target_label)
) +
  geom_vline(xintercept = 0, linewidth = 0.25, color = "#B0B0B0") +
  geom_segment(aes(x = 0, xend = median_depth_adjusted_rho, yend = cluster_label),
    linewidth = 0.45, color = "#666666", alpha = 0.55
  ) +
  geom_point(aes(fill = target_label), size = 3.1, color = "#222222", stroke = 0.4) +
  scale_shape_manual(values = c("LAP3 expression" = 21, "LAP3 state union" = 24), name = NULL) +
  scale_fill_manual(values = c("LAP3 expression" = "#F4A261", "LAP3 state union" = "#2A9D8F"), name = NULL) +
  labs(
    title = "144-gene malignant-state module projection",
    subtitle = "Full module and M1-M3 clusters projected to GBM-Space",
    x = "Median depth-adjusted rho",
    y = NULL
  ) +
  base_theme +
  theme(legend.position = "bottom")

figure4 <- (
  p_a + p_b +
    p_c + p_d +
    p_e + p_f
) +
  plot_layout(widths = c(0.86, 1.16), heights = c(1.0, 1.05, 0.98)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(face = "bold", size = 11),
      plot.margin = margin(8, 8, 8, 8)
    )
  )

fig_base <- file.path(plot_dir, "Figure4_LAP3_State_Topology")
ggsave(paste0(fig_base, ".svg"), figure4, width = 13.8, height = 11.2, units = "in", device = svglite)
ggsave(paste0(fig_base, ".pdf"), figure4, width = 13.8, height = 11.2, units = "in", device = cairo_pdf)
ragg::agg_tiff(paste0(fig_base, ".tiff"), width = 13.8, height = 11.2, units = "in", res = 600, compression = "lzw")
print(figure4)
dev.off()
ragg::agg_png(paste0(fig_base, ".png"), width = 13.8, height = 11.2, units = "in", res = 220)
print(figure4)
dev.off()

plot_files <- list.files(plot_dir, pattern = "^Figure4_LAP3_State_Topology\\.", full.names = TRUE)
plot_qc <- data.table(
  file = basename(plot_files),
  bytes = file.info(plot_files)$size,
  modified = format(file.info(plot_files)$mtime, "%Y-%m-%d %H:%M:%S")
)
if (any(plot_qc$bytes < 10000)) {
  stop("One or more plot files are unexpectedly small:\n", paste(plot_qc$file[plot_qc$bytes < 10000], collapse = "\n"))
}
fwrite(plot_qc, file.path(table_dir, "figure4_export_qc.csv"))

readme <- c(
  "# Figure 4 LAP3 State Topology",
  "",
  sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Figure conclusion",
  "",
  "LAP3_STATE_UNION shows a depth-adjusted GBM-Space topology concentrated in TAM/myeloid and gliosis or malignant-state neighborhoods. The figure supports spatial co-localization/topology, not causal mechanism.",
  "",
  "## Inputs",
  "",
  "- `Data_Spatial_Public/GBM_Space/results/LAP3_State_Spatial_Topology/tables/gbmspace_lap3_state_topology_summary.tsv`",
  "- `Data_Spatial_Public/GBM_Space/results/LAP3_State_Spatial_Topology/tables/gbmspace_lap3_state_contrast_summary.tsv`",
  "- `Data_Spatial_Public/GBM_Space/results/LAP3_State_Spatial_Topology/tables/gbmspace_lap3_state_tumor_effects.tsv`",
  "- `Data_Spatial_Public/GBM_Space/results/LAP3_State_Spatial_Topology/tables/gbmspace_lap3_state_leave_one_tumor_out.tsv`",
  "- `Data_Bulk_TCGA/Data_Merged/results/LAP3_Malignant_State_Module_Audit/tables/gbmspace_malignant_cluster_spatial_summary.tsv`",
  "",
  "## Outputs",
  "",
  "- `plots/Figure4_LAP3_State_Topology.svg`",
  "- `plots/Figure4_LAP3_State_Topology.pdf`",
  "- `plots/Figure4_LAP3_State_Topology.tiff`",
  "- `plots/Figure4_LAP3_State_Topology.png`",
  "- `source_data/figure4_panel_*.csv`",
  "- `tables/figure4_panel_map.csv`",
  "- `tables/figure4_key_results.csv`",
  "- `tables/figure4_export_qc.csv`",
  "",
  "## Interpretation boundary",
  "",
  "This reconstruction is intended as a manuscript-grade spatial topology figure. It intentionally does not claim LAP3-mTORC1 or LAP3-BCAA causal co-localization; mTORC1/BCAA remain primary-pathway boundary checks rather than the main spatial conclusion."
)
writeLines(readme, file.path(out_dir, "README.md"))

message("Wrote source data: ", source_dir)
message("Wrote plots: ", plot_dir)
message("Wrote tables: ", table_dir)
message("DONE Figure 4 LAP3 state topology reconstruction")
