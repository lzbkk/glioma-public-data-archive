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

project_root <- normalizePath(getwd(), mustWork = TRUE)

cache_dir <- file.path(project_root, "Data_Spatial_Public/GBM_Space/results/Lightweight_Cache")
topology_dir <- file.path(project_root, "Data_Spatial_Public/GBM_Space/results/LAP3_State_Spatial_Topology")
out_dir <- file.path(project_root, "Data_Spatial_Public/GBM_Space/results/Figure4_LAP3_State_Topology")
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
  sprintf("figure4_representative_spatial_map.%s.log", format(Sys.time(), "%Y%m%d_%H%M%S"))
)
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)

message("START Figure 4 representative spatial map")
message("Project root: ", project_root)
message("Output dir: ", out_dir)

need_files <- c(
  spot_metadata = file.path(cache_dir, "tables/gbmspace_spot_metadata.tsv"),
  state_scores = file.path(topology_dir, "source_data/gbmspace_spot_lap3_state_scores.tsv.gz"),
  section_effects = file.path(topology_dir, "tables/gbmspace_lap3_state_section_effects.tsv"),
  topology_summary = file.path(topology_dir, "tables/gbmspace_lap3_state_topology_summary.tsv")
)
missing_files <- need_files[!file.exists(need_files)]
if (length(missing_files) > 0) {
  stop("Missing required input files:\n", paste(missing_files, collapse = "\n"))
}

metadata <- fread(need_files[["spot_metadata"]], select = c(
  "spot_id", "h5ad_sample_id", "tumor_id", "spatial_x", "spatial_y", "in_tissue",
  "gene_library_size", "detected_gene_features"
))
state_scores <- fread(need_files[["state_scores"]], select = c(
  "spot_id", "h5ad_sample_id", "tumor_id", "LAP3_log1p_cp10k",
  "LAP3_STATE_UNION", "gene_library_size", "detected_gene_features"
))
section_effects <- fread(need_files[["section_effects"]])
topology_summary <- fread(need_files[["topology_summary"]])

check_cols <- function(x, req, label) {
  miss <- setdiff(req, names(x))
  if (length(miss) > 0) {
    stop(label, " missing columns: ", paste(miss, collapse = ", "))
  }
}
check_cols(metadata, c("spot_id", "h5ad_sample_id", "tumor_id", "spatial_x", "spatial_y"), "metadata")
check_cols(state_scores, c("spot_id", "LAP3_log1p_cp10k", "LAP3_STATE_UNION"), "state_scores")
check_cols(section_effects, c(
  "h5ad_sample_id", "tumor_id", "state_set", "readout", "readout_class",
  "analysis_type", "n_spots", "depth_adjusted_rho"
), "section_effects")
check_cols(topology_summary, c(
  "state_set", "readout", "readout_class", "analysis_type",
  "n_tumors", "median_depth_adjusted_rho", "fdr_depth_adjusted_priority"
), "topology_summary")

context_priority <- c(
  "Proliferative.TAMs",
  "RTN1..TAMs",
  "Monocytes",
  "Resident.BAM.TAMs",
  "Immune..TAMs.",
  "Gliosis.transition",
  "Gliosis",
  "AC.gliosis.like.1"
)
available_context <- intersect(
  context_priority,
  topology_summary[
    state_set == "LAP3_STATE_UNION" &
      analysis_type == "neighbor_topology_k6",
    readout
  ]
)
if (length(available_context) == 0) {
  stop("No preferred context readout available in topology summary.")
}
target_readout <- available_context[1]
target_label <- switch(
  target_readout,
  "Proliferative.TAMs" = "Proliferative TAMs",
  "RTN1..TAMs" = "RTN1+ TAMs",
  "Resident.BAM.TAMs" = "Resident BAM TAMs",
  "Immune..TAMs." = "Immune TAMs",
  gsub("\\.+", " ", target_readout)
)

desired_rho <- topology_summary[
  state_set == "LAP3_STATE_UNION" &
    analysis_type == "neighbor_topology_k6" &
    readout == target_readout,
  median_depth_adjusted_rho
][1]

section_summary <- state_scores[, .(
  n_state_spots = .N,
  lap3_detection_rate = mean(LAP3_log1p_cp10k > 0, na.rm = TRUE),
  state_iqr = as.numeric(IQR(LAP3_STATE_UNION, na.rm = TRUE)),
  lap3_iqr = as.numeric(IQR(LAP3_log1p_cp10k, na.rm = TRUE))
), by = .(h5ad_sample_id, tumor_id)]

candidates <- section_effects[
  state_set == "LAP3_STATE_UNION" &
    analysis_type == "neighbor_topology_k6" &
    readout == target_readout &
    is.finite(depth_adjusted_rho) &
    n_spots >= 1500
]
candidates <- merge(
  candidates,
  section_summary,
  by = c("h5ad_sample_id", "tumor_id"),
  all.x = TRUE
)
candidates <- candidates[
  is.finite(lap3_detection_rate) &
    lap3_detection_rate >= 0.15 &
    lap3_detection_rate <= 0.90 &
    is.finite(state_iqr) &
    state_iqr > 0
]
if (nrow(candidates) == 0) {
  stop("No eligible representative sections after filtering.")
}

candidates[, target_median_rho := desired_rho]
candidates[, rho_distance := abs(depth_adjusted_rho - desired_rho)]
candidates[, spot_distance := abs(n_spots - 3000) / 3000]
candidates[, selection_score := rho_distance + 0.05 * spot_distance - 0.01 * pmin(state_iqr, 1)]
setorder(candidates, selection_score, -n_spots)

non_gene_dir <- file.path(cache_dir, "per_section_non_gene_features")
selected <- NULL
context_features <- NULL
for (i in seq_len(nrow(candidates))) {
  candidate_section <- candidates$h5ad_sample_id[i]
  candidate_file <- file.path(non_gene_dir, paste0(candidate_section, ".non_gene_features.tsv"))
  if (!file.exists(candidate_file)) {
    next
  }
  header <- names(fread(candidate_file, nrows = 0))
  if (!target_readout %in% header) {
    next
  }
  selected <- candidates[i]
  context_features <- fread(candidate_file, select = c("spot_id", "h5ad_sample_id", target_readout))
  setnames(context_features, target_readout, "context_readout")
  break
}
if (is.null(selected)) {
  stop("No eligible section had the target context readout in per-section features.")
}

selected_section <- selected$h5ad_sample_id[1]
selected_tumor <- selected$tumor_id[1]
message("Selected section: ", selected_section)
message("Selected tumor: ", selected_tumor)
message("Target readout: ", target_readout)

map_data <- merge(
  metadata[h5ad_sample_id == selected_section],
  state_scores[h5ad_sample_id == selected_section, .(
    spot_id, LAP3_log1p_cp10k, LAP3_STATE_UNION
  )],
  by = "spot_id",
  all = FALSE
)
map_data <- merge(
  map_data,
  context_features[, .(spot_id, context_readout)],
  by = "spot_id",
  all.x = TRUE
)
if (nrow(map_data) == 0) {
  stop("Selected section produced an empty map data table.")
}
map_data[, y_plot := -spatial_y]
map_data[, selected_context_readout := target_readout]
map_data[, selected_context_label := target_label]

selection_table <- copy(candidates)
selection_table[, selected := h5ad_sample_id == selected_section]
selection_table <- selection_table[, .(
  selected, h5ad_sample_id, tumor_id, readout, readout_class, n_spots,
  depth_adjusted_rho, target_median_rho, rho_distance, lap3_detection_rate,
  state_iqr, lap3_iqr, selection_score
)]

key_results <- data.table(
  metric = c(
    "selected_section",
    "selected_tumor",
    "target_context_readout",
    "section_n_spots",
    "section_depth_adjusted_rho_to_context",
    "target_median_depth_adjusted_rho",
    "lap3_detection_rate",
    "state_iqr"
  ),
  value = c(
    selected_section,
    selected_tumor,
    target_readout,
    as.character(selected$n_spots[1]),
    sprintf("%.6f", selected$depth_adjusted_rho[1]),
    sprintf("%.6f", desired_rho),
    sprintf("%.6f", selected$lap3_detection_rate[1]),
    sprintf("%.6f", selected$state_iqr[1])
  )
)

fwrite(map_data, file.path(source_dir, "figure4_representative_spatial_map.csv"))
fwrite(selection_table, file.path(table_dir, "figure4_representative_spatial_map_selection_candidates.csv"))
fwrite(key_results, file.path(table_dir, "figure4_representative_spatial_map_key_results.csv"))

plot_theme <- theme_void(base_family = "sans", base_size = 7) +
  theme(
    plot.title = element_text(face = "bold", size = 8.5, hjust = 0, margin = margin(b = 3)),
    plot.subtitle = element_text(size = 6.2, color = "#4D4D4D", hjust = 0, margin = margin(b = 3)),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 6.5),
    legend.text = element_text(size = 6),
    legend.key.width = unit(13, "mm"),
    legend.key.height = unit(2.2, "mm"),
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
    legend.box.margin = margin(t = -2, r = 0, b = 0, l = 0),
    plot.margin = margin(4, 4, 4, 4)
  )

clip_and_scale <- function(x, probs = c(0.01, 0.99)) {
  q <- stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE)
  clipped <- pmin(pmax(x, q[1]), q[2])
  denom <- q[2] - q[1]
  if (!is.finite(denom) || denom <= 0) {
    return(rep(0, length(x)))
  }
  (clipped - q[1]) / denom
}

plot_map <- function(dt, value_col, title, subtitle, palette, legend_title) {
  plot_dt <- copy(dt)
  plot_dt[, plot_value := clip_and_scale(get(value_col))]
  ggplot(plot_dt, aes(x = spatial_x, y = y_plot, color = plot_value)) +
    geom_point(size = 0.32, alpha = 0.95, stroke = 0) +
    coord_fixed(expand = FALSE) +
    scale_color_gradientn(
      colors = palette,
      name = paste0(legend_title, "\nscaled intensity"),
      limits = c(0, 1),
      breaks = c(0, 0.5, 1),
      labels = c("low", "mid", "high")
    ) +
    labs(title = title, subtitle = subtitle) +
    guides(color = guide_colorbar(
      direction = "horizontal",
      barwidth = unit(17, "mm"),
      barheight = unit(2.2, "mm"),
      title.position = "top",
      title.hjust = 0.5,
      label.position = "bottom"
    )) +
    plot_theme
}

state_palette <- c("#F4F6F8", "#BBD7EA", "#4D92C6", "#1F4E79", "#7A2F1B")
lap3_palette <- c("#F7F7F7", "#D9EAD3", "#7FB069", "#2D6A4F", "#7A3B12")
context_palette <- c("#F7F7F7", "#F4D6C6", "#DB8A66", "#A9433C", "#5A1E24")

section_note <- sprintf(
  "%s; n=%s spots",
  selected_section,
  format(nrow(map_data), big.mark = ",")
)
rho_note <- sprintf(
  "depth-adjusted rho = %.3f",
  selected$depth_adjusted_rho[1]
)

p_state <- plot_map(
  map_data,
  "LAP3_STATE_UNION",
  "A  LAP3-state",
  section_note,
  state_palette,
  "Union score"
)
p_lap3 <- plot_map(
  map_data,
  "LAP3_log1p_cp10k",
  "B  LAP3 expression",
  "log1p(CP10K)",
  lap3_palette,
  "LAP3"
)
p_context <- plot_map(
  map_data,
  "context_readout",
  paste0("C  ", target_label),
  rho_note,
  context_palette,
  target_label
)

combined <- p_state + p_lap3 + p_context +
  plot_layout(ncol = 3, widths = c(1, 1, 1)) +
  plot_annotation(
    title = "Representative GBM-Space Visium section",
    subtitle = "Panel values are 1%-99% clipped and scaled for visualization; source data preserve original values.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 10, margin = margin(b = 2)),
      plot.subtitle = element_text(size = 7.2, color = "#4D4D4D", margin = margin(b = 4))
    )
  )

save_pub <- function(plot, filename, width_mm = 183, height_mm = 62, dpi = 600) {
  w <- width_mm / 25.4
  h <- height_mm / 25.4
  svglite::svglite(paste0(filename, ".svg"), width = w, height = h)
  print(plot)
  dev.off()
  grDevices::cairo_pdf(paste0(filename, ".pdf"), width = w, height = h, family = "sans")
  print(plot)
  dev.off()
  ragg::agg_tiff(paste0(filename, ".tiff"), width = w, height = h, units = "in", res = dpi, compression = "lzw")
  print(plot)
  dev.off()
  ragg::agg_png(paste0(filename, ".png"), width = w, height = h, units = "in", res = 220)
  print(plot)
  dev.off()
}

plot_base <- file.path(plot_dir, "Figure4_Representative_Spatial_Map")
save_pub(combined, plot_base)

plot_files <- paste0(plot_base, c(".svg", ".pdf", ".tiff", ".png"))
export_qc <- data.table(
  file = basename(plot_files),
  path = plot_files,
  exists = file.exists(plot_files),
  size_bytes = as.numeric(file.info(plot_files)$size)
)
export_qc[, nonempty := exists & size_bytes > 1000]
fwrite(export_qc, file.path(table_dir, "figure4_representative_spatial_map_export_qc.csv"))
if (!all(export_qc$nonempty)) {
  stop("One or more representative spatial map exports are missing or empty.")
}

readme <- c(
  "# Figure 4 Representative Spatial Map",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Figure conclusion",
  "",
  "A representative GBM-Space Visium section visually shows the spatial organization of `LAP3_STATE_UNION` alongside LAP3 expression and a TAM/myeloid context readout. This is a visual companion to the tumor-aware topology statistics, not an independent causal test.",
  "",
  "## Selection rule",
  "",
  paste0("- Target context readout: `", target_readout, "` (", target_label, ")."),
  "- Candidate sections required at least 1,500 spots, finite depth-adjusted topology effect, nonzero LAP3-state heterogeneity, and LAP3 detection rate between 0.15 and 0.90.",
  "- The selected section is closest to the tumor-level median depth-adjusted topology effect for the target context readout, with a small tie-breaker for spot count and state heterogeneity.",
  "",
  "## Selected section",
  "",
  paste0("- Section: `", selected_section, "`."),
  paste0("- Tumor: `", selected_tumor, "`."),
  paste0("- Spots in map: ", format(nrow(map_data), big.mark = ","), "."),
  paste0("- Section depth-adjusted rho to ", target_label, ": ", sprintf("%.3f", selected$depth_adjusted_rho[1]), "."),
  paste0("- Tumor-level median depth-adjusted rho for ", target_label, ": ", sprintf("%.3f", desired_rho), "."),
  "",
  "## Inputs",
  "",
  "- `Data_Spatial_Public/GBM_Space/results/Lightweight_Cache/tables/gbmspace_spot_metadata.tsv`",
  "- `Data_Spatial_Public/GBM_Space/results/LAP3_State_Spatial_Topology/source_data/gbmspace_spot_lap3_state_scores.tsv.gz`",
  "- `Data_Spatial_Public/GBM_Space/results/Lightweight_Cache/per_section_non_gene_features/`",
  "- `Data_Spatial_Public/GBM_Space/results/LAP3_State_Spatial_Topology/tables/gbmspace_lap3_state_section_effects.tsv`",
  "- `Data_Spatial_Public/GBM_Space/results/LAP3_State_Spatial_Topology/tables/gbmspace_lap3_state_topology_summary.tsv`",
  "",
  "## Outputs",
  "",
  "- `plots/Figure4_Representative_Spatial_Map.svg`",
  "- `plots/Figure4_Representative_Spatial_Map.pdf`",
  "- `plots/Figure4_Representative_Spatial_Map.tiff`",
  "- `plots/Figure4_Representative_Spatial_Map.png`",
  "- `source_data/figure4_representative_spatial_map.csv`",
  "- `tables/figure4_representative_spatial_map_selection_candidates.csv`",
  "- `tables/figure4_representative_spatial_map_key_results.csv`",
  "- `tables/figure4_representative_spatial_map_export_qc.csv`",
  "",
  "## Interpretation boundary",
  "",
  "This map should be described as representative spatial visualization of the LAP3-state topology. It should not be used as spot-level causal evidence or as proof of a pathway-only LAP3-mTORC1/BCAA spatial mechanism."
)
writeLines(readme, file.path(out_dir, "README_Representative_Spatial_Map.md"))

message("DONE Figure 4 representative spatial map")
message("Selected section: ", selected_section)
message("Plot base: ", plot_base)
message("Log file: ", log_file)
