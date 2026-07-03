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
  file_arg <- "Data_Bulk_TCGA/Data_Merged/scripts/figure1_lap3_context_benchmark.R"
}
repo <- normalizePath(file.path(dirname(file_arg), "../../.."), mustWork = FALSE)
if (!dir.exists(file.path(repo, "Project_Management"))) repo <- normalizePath(".", mustWork = TRUE)

state_dir <- file.path(repo, "Data_Bulk_TCGA/Data_Merged/results/LAP3_State_Module")
spec_dir <- file.path(repo, "Data_Bulk_TCGA/Data_Merged/results/LAP3_State_Specificity_Benchmark")
composition_dir <- file.path(repo, "Data_Bulk_TCGA/Data_Merged/results/LAP3_Bulk_Composition_Audit")
out_dir <- file.path(repo, "Data_Bulk_TCGA/Data_Merged/results/Figure1_LAP3_Context_Benchmark")
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
write_source <- function(dt, name) {
  fwrite(copy(dt), file.path(source_dir, name))
}
write_table <- function(dt, name) {
  fwrite(copy(dt), file.path(table_dir, name))
}

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

label_dataset <- c(
  "TCGA" = "TCGA",
  "CGGA_mRNAseq_693" = "CGGA693",
  "CGGA_mRNAseq_325" = "CGGA325",
  "mRNAseq_693" = "CGGA693",
  "mRNAseq_325" = "CGGA325",
  "CGGA" = "CGGA"
)
clean_dataset <- function(x) {
  out <- as.character(x)
  idx <- out %in% names(label_dataset)
  out[idx] <- unname(label_dataset[out[idx]])
  out
}

tcga_projection <- read_dt(file.path(state_dir, "tables/tcga_lap3_state_score_projection.csv"))
cgga_projection <- read_dt(file.path(state_dir, "tables/cgga_lap3_state_score_projection.csv"))
first_summary <- read_dt(file.path(state_dir, "tables/lap3_state_first_pass_summary.csv"))
clinical_tests <- read_dt(file.path(spec_dir, "tables/lap3_state_clinical_group_tests.csv"))
benchmark_cor <- read_dt(file.path(spec_dir, "tables/lap3_state_benchmark_score_correlations.csv"))
random_benchmark <- read_dt(file.path(spec_dir, "tables/expression_matched_random_gene_benchmark_tcga.csv"))
adjusted_models <- read_dt(file.path(spec_dir, "tables/lap3_state_adjusted_lm_models.csv"))
composition_cor <- read_dt(file.path(composition_dir, "tables/lap3_state_composition_proxy_correlations.csv"))
composition_verdict <- read_dt(file.path(composition_dir, "tables/lap3_state_primary_composition_verdict.csv"))

stopifnot(nrow(tcga_projection) > 0, nrow(cgga_projection) > 0)

panel_a <- rbindlist(list(
  data.table(
    cohort = "TCGA",
    n = nrow(tcga_projection),
    layer = "Discovery",
    role = "Primary tumors"
  ),
  cgga_projection[, .(n = .N), by = cohort][
    , .(cohort = clean_dataset(cohort), n, layer = "External projection", role = "Primary glioma")
  ],
  data.table(
    cohort = c("Composition audit", "Figure 4", "Figure 5"),
    n = c(3, 12, 2),
    layer = c("Purity/TME", "Spatial topology", "Multi-omic"),
    role = c("TCGA + CGGA cohorts", "GBM-Space tumors", "CPTAC + GLASS")
  )
), fill = TRUE)
panel_a[, y := rev(seq_len(.N))]
panel_a[, label_left := paste0(cohort, "\n", layer)]
panel_a[cohort == "CGGA693", label_left := "CGGA693\nValidation"]
panel_a[cohort == "CGGA325", label_left := "CGGA325\nValidation"]
panel_a[cohort == "Composition audit", label_left := "Composition\nBoundary"]
panel_a[cohort == "Figure 5", label_left := "Figure 5\nMulti-omic"]
panel_a[, label_right := fifelse(cohort == "Composition audit", "3 cohorts",
  fifelse(cohort == "Figure 4", "12 tumors",
    fifelse(cohort == "Figure 5", "CPTAC + GLASS", paste0("n=", n))
  )
)]

clinical_keep <- c("cohort", "dataset", "tumor_class", "grade", "idh_status", "codel_1p19q", "mgmt_status")
clinical_labels <- c(
  cohort = "TCGA GBM/LGG",
  dataset = "CGGA platform",
  tumor_class = "CGGA GBM/LGG",
  grade = "Grade",
  idh_status = "IDH",
  codel_1p19q = "1p/19q",
  mgmt_status = "MGMT"
)
panel_b <- clinical_tests[
  score == "LAP3_STATE_UNION" & variable %chin% clinical_keep
]
panel_b[, dataset_label := clean_dataset(dataset)]
panel_b[, variable_label := clinical_labels[variable]]
panel_b[, neg_log10_fdr := pmin(-log10(p_adj_BH), 30)]
panel_b[, variable_label := factor(variable_label, levels = rev(c("TCGA GBM/LGG", "CGGA platform", "CGGA GBM/LGG", "Grade", "IDH", "1p/19q", "MGMT")))]
panel_b[, dataset_label := factor(dataset_label, levels = c("TCGA", "CGGA"))]

panel_c <- data.table(
  dataset = c("TCGA", "CGGA693", "CGGA325"),
  n = c(nrow(tcga_projection), sum(cgga_projection$cohort == "mRNAseq_693"), sum(cgga_projection$cohort == "mRNAseq_325")),
  spearman_rho = c(
    as.numeric(first_summary[metric == "tcga_lap3_state_lap3_rho", value]),
    as.numeric(first_summary[metric == "cgga693_lap3_state_lap3_rho", value]),
    as.numeric(first_summary[metric == "cgga325_lap3_state_lap3_rho", value])
  )
)
panel_c[, dataset := factor(dataset, levels = dataset)]
panel_c[, label := paste0("n=", n)]

panel_d <- copy(random_benchmark)
panel_d[, group := fifelse(is_lap3, "LAP3", "Expression-matched genes")]
lap3_row <- panel_d[is_lap3 == TRUE][1]

benchmark_keep <- c(
  "LAP3_log2_expr",
  "TAM_MYELOID_CORE",
  "HALLMARK_HYPOXIA",
  "HALLMARK_MTORC1_SIGNALING",
  "LEUCINE_BCAA_CORE",
  "REACTOME_TRANSLATION",
  "PROLIFERATION_CORE",
  "AMINOPEPTIDASE_FAMILY"
)
benchmark_labels <- c(
  LAP3_log2_expr = "LAP3 expression",
  TAM_MYELOID_CORE = "TAM/myeloid",
  HALLMARK_HYPOXIA = "Hypoxia",
  HALLMARK_MTORC1_SIGNALING = "mTORC1",
  LEUCINE_BCAA_CORE = "BCAA module",
  REACTOME_TRANSLATION = "Translation",
  PROLIFERATION_CORE = "Proliferation",
  AMINOPEPTIDASE_FAMILY = "Aminopeptidase family"
)
panel_e <- benchmark_cor[
  score == "LAP3_STATE_UNION" &
    group == "all" &
    dataset %chin% c("TCGA", "CGGA_mRNAseq_693", "CGGA_mRNAseq_325") &
    variable %chin% benchmark_keep
]
panel_e[, dataset_label := clean_dataset(dataset)]
panel_e[, variable_label := benchmark_labels[variable]]
panel_e[, variable_label := factor(variable_label, levels = rev(benchmark_labels[benchmark_keep]))]
panel_e[, dataset_label := factor(dataset_label, levels = c("TCGA", "CGGA693", "CGGA325"))]
panel_e[, sig_label := fifelse(p_adj_BH < 0.05, "FDR < 0.05", "FDR >= 0.05")]

panel_f <- copy(composition_verdict)
panel_f[, dataset_label := clean_dataset(dataset)]
panel_f[, dataset_label := factor(dataset_label, levels = c("TCGA", "CGGA693", "CGGA325"))]
panel_f_long <- melt(
  panel_f,
  id.vars = c("dataset", "dataset_label", "n_clinical", "audit_call"),
  measure.vars = c("beta_clinical", "beta_composition"),
  variable.name = "model",
  value.name = "beta"
)
panel_f_long[, model_label := fifelse(model == "beta_clinical", "Clinical", "Clinical + composition")]
panel_f_long[, model_label := factor(model_label, levels = c("Clinical", "Clinical + composition"))]

key_results <- rbindlist(list(
  panel_c[, .(panel = "C", dataset = as.character(dataset), feature = "LAP3_STATE_UNION vs LAP3", n, effect = spearman_rho, statistic = "Spearman rho")],
  panel_e[, .(panel = "E", dataset = as.character(dataset_label), feature = as.character(variable_label), n, effect = spearman_rho, statistic = "Spearman rho")],
  panel_f[, .(panel = "F", dataset = as.character(dataset_label), feature = "LAP3 term retained after composition adjustment", n = n_clinical, effect = beta_composition, statistic = "standardized beta")]
), fill = TRUE)

panel_map <- data.table(
  panel = LETTERS[1:6],
  source_data = c(
    "figure1_panel_a_evidence_frame.csv",
    "figure1_panel_b_clinical_molecular_tests.csv",
    "figure1_panel_c_lap3_state_anchor_reproducibility.csv",
    "figure1_panel_d_expression_matched_random_benchmark.csv",
    "figure1_panel_e_state_context_benchmark_correlations.csv",
    "figure1_panel_f_composition_adjusted_boundary.csv"
  ),
  conclusion = c(
    "Figure 1 defines the cohorts and evidence roles used to establish LAP3 as a state anchor.",
    "The frozen LAP3-state is associated with major clinical and molecular glioma structure.",
    "LAP3-state is reproducibly anchored to LAP3 expression across TCGA and CGGA cohorts.",
    "LAP3 ranks near the extreme of expression-matched genes for state anchoring.",
    "The state is coupled to TAM/myeloid, hypoxia, mTORC1/BCAA and translation context.",
    "The LAP3 term persists after lightweight composition adjustment but is attenuated."
  )
)

write_source(panel_a, "figure1_panel_a_evidence_frame.csv")
write_source(panel_b, "figure1_panel_b_clinical_molecular_tests.csv")
write_source(panel_c, "figure1_panel_c_lap3_state_anchor_reproducibility.csv")
write_source(panel_d, "figure1_panel_d_expression_matched_random_benchmark.csv")
write_source(panel_e, "figure1_panel_e_state_context_benchmark_correlations.csv")
write_source(panel_f_long, "figure1_panel_f_composition_adjusted_boundary.csv")
write_table(panel_map, "figure1_panel_map.csv")
write_table(key_results, "figure1_key_results.csv")
write_table(adjusted_models, "figure1_adjusted_lm_models_input.csv")
write_table(composition_cor, "figure1_composition_proxy_correlations_input.csv")

col_layer <- c(
  "Discovery" = "#2F6F9F",
  "Discovery and benchmark" = "#2F6F9F",
  "External projection" = "#537A5A",
  "Purity/TME boundary" = "#7D8390",
  "Purity/TME" = "#7D8390",
  "Spatial topology" = "#B15A2C",
  "Multi-omic boundary" = "#7B5B8F",
  "Multi-omic" = "#7B5B8F"
)
col_sig <- c("FDR < 0.05" = "#1F6F8B", "FDR >= 0.05" = "#B7BDC5")
col_dataset <- c("TCGA" = "#315E8A", "CGGA" = "#BA6A3C", "CGGA693" = "#BA6A3C", "CGGA325" = "#7A9B55")

p_a <- ggplot(panel_a) +
  geom_rect(aes(xmin = 0, xmax = 1, ymin = y - 0.42, ymax = y + 0.42, fill = layer),
            colour = "white", linewidth = 0.35) +
  geom_text(aes(x = 0.03, y = y, label = label_left), hjust = 0, colour = "white",
            size = 2.15, lineheight = 0.9, fontface = "bold") +
  geom_text(aes(x = 0.94, y = y, label = label_right), hjust = 1, colour = "white",
            size = 1.85, lineheight = 0.88) +
  scale_fill_manual(values = col_layer, guide = "none") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0.45, nrow(panel_a) + 0.55), expand = FALSE) +
  labs(title = "A  Evidence frame") +
  theme_void(base_family = "sans") +
  theme(plot.title = element_text(face = "bold", size = 8, hjust = 0), plot.margin = margin(4, 4, 4, 4))

p_b <- ggplot(panel_b, aes(x = dataset_label, y = variable_label)) +
  geom_point(aes(size = neg_log10_fdr, fill = neg_log10_fdr), shape = 21, colour = "#3A3D42", stroke = 0.2) +
  scale_fill_gradient(low = "#DDE3EA", high = "#2F6F9F", name = "-log10 FDR", limits = c(0, 30)) +
  scale_size_continuous(range = c(1.5, 5), name = "-log10 FDR", limits = c(0, 30)) +
  guides(size = "none") +
  labs(title = "B  Clinical-molecular structure", x = NULL, y = NULL) +
  theme_fig() +
  theme(legend.position = "right")

p_c <- ggplot(panel_c, aes(x = dataset, y = spearman_rho)) +
  geom_col(fill = "#315E8A", width = 0.58) +
  geom_text(aes(label = sprintf("rho=%.3f\n%s", spearman_rho, label)), vjust = -0.2, size = 2.2, lineheight = 0.9) +
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "C  Cross-cohort state anchor", x = NULL, y = "Spearman rho with LAP3") +
  theme_fig()

p_d <- ggplot(panel_d, aes(x = spearman_rho)) +
  geom_histogram(data = panel_d[is_lap3 == FALSE], bins = 40, fill = "#CBD2DA", colour = "white", linewidth = 0.15) +
  geom_vline(xintercept = lap3_row$spearman_rho, colour = "#B15A2C", linewidth = 0.65) +
  annotate("text", x = lap3_row$spearman_rho, y = Inf,
           label = sprintf("LAP3 rank %d/1001\npercentile %.3f", lap3_row$rank_abs_rho, lap3_row$percentile_abs_rho),
           hjust = 1.02, vjust = 1.15, size = 2.25, colour = "#7F3B1D") +
  labs(title = "D  Expression-matched benchmark", x = "Spearman rho with LAP3-state", y = "Genes") +
  theme_fig()

p_e <- ggplot(panel_e, aes(x = spearman_rho, y = variable_label, colour = dataset_label)) +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "#7C828A") +
  geom_segment(aes(x = 0, xend = spearman_rho, yend = variable_label), linewidth = 0.45, alpha = 0.65,
               position = position_dodge(width = 0.45)) +
  geom_point(size = 2.2, stroke = 0.4, position = position_dodge(width = 0.45)) +
  scale_colour_manual(values = col_dataset, name = NULL) +
  scale_x_continuous(limits = c(-0.05, 1), breaks = seq(0, 1, 0.25)) +
  labs(title = "E  State-context benchmark", x = "Spearman rho with LAP3-state", y = NULL) +
  theme_fig()

p_f <- ggplot(panel_f_long, aes(x = beta, y = dataset_label, colour = model_label)) +
  geom_segment(
    data = panel_f,
    aes(x = beta_clinical, xend = beta_composition, y = dataset_label, yend = dataset_label),
    inherit.aes = FALSE, linewidth = 0.45, colour = "#A8AFB8"
  ) +
  geom_point(size = 2.3, position = position_dodge(width = 0.25)) +
  scale_colour_manual(values = c("Clinical" = "#315E8A", "Clinical + composition" = "#B15A2C"), name = NULL) +
  scale_x_continuous(limits = c(0, 0.92), breaks = seq(0, 0.9, 0.3)) +
  labs(title = "F  Composition-aware boundary", x = "Standardized LAP3 beta for LAP3-state", y = NULL) +
  theme_fig()

caption <- paste(
  "Figure 1 defines LAP3 as a reproducible glioma state anchor rather than a standalone causal marker.\n",
  "Composition adjustment tests broad immune/stromal/TAM confounding; it does not establish purity-independent causality."
)

figure1 <- (p_a | p_b) / (p_c | p_d) / (p_e | p_f) +
  plot_layout(heights = c(1.0, 1.0, 1.25)) +
  plot_annotation(
    title = "Figure 1. LAP3 anchors a non-trivial glioma state context",
    caption = caption,
    theme = theme(
      plot.title = element_text(face = "bold", size = 10, family = "sans"),
      plot.caption = element_text(size = 6.5, colour = "#3A3D42", hjust = 0, family = "sans"),
      plot.margin = margin(5, 5, 5, 5)
    )
  )

plot_base <- file.path(plot_dir, "Figure1_LAP3_Context_Benchmark")
svglite::svglite(paste0(plot_base, ".svg"), width = 7.2, height = 7.6)
print(figure1)
dev.off()
ggsave(paste0(plot_base, ".pdf"), figure1, width = 7.2, height = 7.6, units = "in", device = grDevices::pdf)
ragg::agg_tiff(paste0(plot_base, ".tiff"), width = 7.2, height = 7.6, units = "in", res = 600, compression = "lzw")
print(figure1)
dev.off()
ragg::agg_png(paste0(plot_base, ".png"), width = 7.2, height = 7.6, units = "in", res = 300)
print(figure1)
dev.off()

export_files <- c(
  paste0(plot_base, c(".svg", ".pdf", ".tiff", ".png")),
  file.path(source_dir, panel_map$source_data),
  file.path(table_dir, c(
    "figure1_panel_map.csv",
    "figure1_key_results.csv",
    "figure1_adjusted_lm_models_input.csv",
    "figure1_composition_proxy_correlations_input.csv"
  ))
)
export_qc <- data.table(
  file = export_files,
  exists = file.exists(export_files),
  bytes = as.numeric(file.info(export_files)$size)
)
write_table(export_qc, "figure1_export_qc.csv")

readme <- c(
  "# Figure 1 LAP3 Context Benchmark",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Purpose",
  "",
  "Create the active manuscript Figure 1 for the state/ecosystem route.",
  "",
  "## Core Conclusion",
  "",
  "`LAP3_STATE_UNION` is reproducibly anchored to LAP3 expression and glioma clinical-molecular structure, but the signal is embedded in TAM/myeloid, hypoxia and metabolic context; lightweight composition adjustment supports a non-trivial LAP3-centered state rather than purity-independent causality.",
  "",
  "## Inputs",
  "",
  "- `Data_Bulk_TCGA/Data_Merged/results/LAP3_State_Module/`",
  "- `Data_Bulk_TCGA/Data_Merged/results/LAP3_State_Specificity_Benchmark/`",
  "- `Data_Bulk_TCGA/Data_Merged/results/LAP3_Bulk_Composition_Audit/`",
  "",
  "## Panels",
  "",
  "- Panel A: evidence frame.",
  "- Panel B: clinical-molecular group tests for `LAP3_STATE_UNION`.",
  "- Panel C: cross-cohort LAP3-state/LAP3 anchor correlations.",
  "- Panel D: TCGA expression-matched random gene benchmark.",
  "- Panel E: state-context benchmark correlations.",
  "- Panel F: clinical-only versus composition-adjusted LAP3 term.",
  "",
  "## Outputs",
  "",
  "- `plots/Figure1_LAP3_Context_Benchmark.svg`",
  "- `plots/Figure1_LAP3_Context_Benchmark.pdf`",
  "- `plots/Figure1_LAP3_Context_Benchmark.tiff`",
  "- `plots/Figure1_LAP3_Context_Benchmark.png`",
  "- `source_data/figure1_panel_*.csv`",
  "- `tables/figure1_key_results.csv`",
  "- `tables/figure1_panel_map.csv`",
  "- `tables/figure1_export_qc.csv`",
  "",
  "## Interpretation Boundary",
  "",
  "Figure 1 should be written as a clinical-molecular and specificity benchmark for a LAP3-centered state. It should not be written as proof that LAP3 is a purity-independent causal driver."
)
writeLines(readme, file.path(out_dir, "README.md"))

if (!all(export_qc$exists) || any(is.na(export_qc$bytes) | export_qc$bytes <= 0)) {
  stop("Export QC failed: ", file.path(table_dir, "figure1_export_qc.csv"))
}

message("Figure 1 export completed: ", out_dir)
