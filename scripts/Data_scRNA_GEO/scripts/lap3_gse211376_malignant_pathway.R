#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(msigdbr)
})

setwd("/home/lzb/glioma")
data.table::setDTthreads(8)
set.seed(20260629)

matrix_file <- "Data_scRNA_GEO/GSE211376/GSE211376_raw_counts_Ruiz2022_all_samples_filtered_cells.tsv.gz"
metadata_file <- "Data_scRNA_GEO/GSE211376/GSE211376_metadata_Ruiz2022_all_samples_filtered_cells.csv.gz"
out_dir <- "Data_scRNA_GEO/results/GSE211376_LAP3_Malignant_Pathway"
table_dir <- file.path(out_dir, "tables")
plot_dir <- file.path(out_dir, "plots")
log_dir <- file.path(out_dir, "logs")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(log_dir, "lap3_gse211376_malignant_pathway.log")
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message")
on.exit({
  sink(type = "message")
  sink(type = "output")
  close(log_con)
}, add = TRUE)

cat("Started:", format(Sys.time()), "\n")
cat("R:", R.version.string, "\n")
cat("data.table threads:", data.table::getDTthreads(), "\n")

write_table <- function(x, filename) {
  data.table::fwrite(x, file.path(table_dir, filename))
}

fmt_p <- function(p) {
  ifelse(
    is.na(p), "NA",
    ifelse(p < 1e-4, formatC(p, format = "e", digits = 2), formatC(p, format = "f", digits = 3))
  )
}

cor_test_safe <- function(x, y, min_n = 6) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) < min_n || length(unique(x)) < 3 || length(unique(y)) < 3) {
    return(c(n = length(x), rho = NA_real_, ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_))
  }
  test <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  boot_rho <- replicate(2000, {
    index <- sample.int(length(x), replace = TRUE)
    suppressWarnings(cor(x[index], y[index], method = "spearman"))
  })
  ci <- quantile(boot_rho[is.finite(boot_rho)], c(0.025, 0.975), na.rm = TRUE)
  c(
    n = length(x),
    rho = unname(test$estimate),
    ci_low = unname(ci[1]),
    ci_high = unname(ci[2]),
    p_value = test$p.value
  )
}

make_gene_sets <- function() {
  hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>%
    select(gs_name, gene_symbol)
  reactome <- msigdbr(
    species = "Homo sapiens",
    category = "C2",
    subcategory = "CP:REACTOME"
  ) %>%
    select(gs_name, gene_symbol)
  pathways <- split(
    bind_rows(hallmark, reactome)$gene_symbol,
    bind_rows(hallmark, reactome)$gs_name
  )
  selected <- c(
    "HALLMARK_MTORC1_SIGNALING",
    "HALLMARK_MYC_TARGETS_V1",
    "HALLMARK_E2F_TARGETS",
    "HALLMARK_G2M_CHECKPOINT",
    "REACTOME_TRANSLATION"
  )
  manual <- list(
    LEUCINE_BCAA_CORE = c(
      "BCAT1", "BCAT2", "BCKDHA", "BCKDHB", "DBT", "DLD", "BCKDK", "PPM1K",
      "SLC7A5", "SLC3A2", "SLC43A1", "SLC43A2", "SLC38A2", "SLC38A9",
      "MTOR", "RPTOR", "RRAGA", "RRAGB", "RRAGC", "RRAGD", "LAMTOR1", "LAMTOR2",
      "LAMTOR3", "LAMTOR4", "LAMTOR5", "RHEB", "RPS6KB1", "EIF4EBP1", "RPS6"
    ),
    MTORC1_READOUT_CORE = c(
      "MTOR", "RPTOR", "RHEB", "AKT1", "TSC1", "TSC2", "RRAGA", "RRAGB", "RRAGC",
      "RRAGD", "LAMTOR1", "LAMTOR2", "LAMTOR3", "LAMTOR4", "LAMTOR5",
      "RPS6KB1", "RPS6KB2", "EIF4EBP1", "EIF4EBP2", "RPS6", "EIF4E"
    )
  )
  sets <- c(pathways[intersect(selected, names(pathways))], manual)
  lapply(sets, function(x) setdiff(unique(x), "LAP3"))
}

metadata <- read.csv(gzfile(metadata_file), stringsAsFactors = FALSE, check.names = FALSE)
metadata$cell_id <- rownames(metadata)
rownames(metadata) <- NULL
names(metadata)[names(metadata) == "predicted.high_hierarchy"] <- "cell_state"

malignant_states <- c("AC-like", "MES-like", "NPC-like", "OPC-like")
tam_states <- c("TAM-BDM", "TAM-MG")
metadata <- metadata %>%
  mutate(
    compartment = case_when(
      cell_state %in% malignant_states ~ "Malignant",
      cell_state %in% tam_states ~ "TAM",
      TRUE ~ "Other"
    )
  )

stopifnot(
  nrow(metadata) == 39355,
  length(unique(metadata$patient)) == 11,
  all(malignant_states %in% metadata$cell_state),
  all(tam_states %in% metadata$cell_state)
)

gene_sets <- make_gene_sets()
target_genes <- unique(c("LAP3", unlist(gene_sets, use.names = FALSE)))
cat("Target genes requested:", length(target_genes), "\n")

cell_ids <- strsplit(readLines(gzfile(matrix_file), n = 1), "\t", fixed = TRUE)[[1]]
stopifnot(
  length(cell_ids) == nrow(metadata),
  identical(cell_ids, metadata$cell_id)
)

cat("Loading full count matrix with fread...\n")
load_start <- Sys.time()
counts_dt <- data.table::fread(
  cmd = paste("zcat", shQuote(matrix_file)),
  skip = 1,
  header = FALSE,
  sep = "\t",
  data.table = TRUE,
  showProgress = TRUE,
  nThread = 8
)
cat("Loaded matrix in", round(difftime(Sys.time(), load_start, units = "mins"), 2), "minutes\n")
stopifnot(nrow(counts_dt) == 27102, ncol(counts_dt) == length(cell_ids) + 1)
setnames(counts_dt, 1, "gene")

genes_present <- counts_dt$gene
coverage <- bind_rows(lapply(names(gene_sets), function(pathway) {
  requested <- gene_sets[[pathway]]
  present <- intersect(requested, genes_present)
  data.frame(
    pathway = pathway,
    genes_requested = length(requested),
    genes_present = length(present),
    coverage = length(present) / length(requested),
    present_genes = paste(present, collapse = ";"),
    missing_genes = paste(setdiff(requested, present), collapse = ";")
  )
}))
write_table(coverage, "gse211376_pathway_gene_coverage.csv")
stopifnot(all(coverage$genes_present >= 10), "LAP3" %in% genes_present)

cat("Computing per-cell library sizes...\n")
library_sizes <- vapply(counts_dt[, -1], sum, numeric(1))
stopifnot(length(library_sizes) == nrow(metadata), all(library_sizes > 0))

selected_dt <- counts_dt[gene %in% target_genes]
selected_genes <- selected_dt$gene
selected_mat <- as.matrix(selected_dt[, -1])
storage.mode(selected_mat) <- "numeric"
rownames(selected_mat) <- selected_genes
colnames(selected_mat) <- cell_ids
rm(counts_dt)
invisible(gc())

lap3_counts <- selected_mat["LAP3", ]
metadata$library_size <- library_sizes
metadata$lap3_count <- lap3_counts
metadata$lap3_detected <- lap3_counts > 0
metadata$lap3_log1p_cpm <- log1p(1e6 * lap3_counts / library_sizes)

cell_qc <- metadata %>%
  group_by(compartment, cell_state, patient) %>%
  summarise(
    n_cells = n(),
    lap3_positive_cells = sum(lap3_detected),
    pct_lap3_positive = 100 * mean(lap3_detected),
    mean_lap3_log1p_cpm = mean(lap3_log1p_cpm),
    median_lap3_log1p_cpm = median(lap3_log1p_cpm),
    .groups = "drop"
  )
write_table(cell_qc, "gse211376_lap3_detection_by_patient_state.csv")

state_summary <- metadata %>%
  group_by(compartment, cell_state) %>%
  summarise(
    n_cells = n(),
    n_patients = n_distinct(patient),
    lap3_positive_cells = sum(lap3_detected),
    pct_lap3_positive = 100 * mean(lap3_detected),
    mean_lap3_log1p_cpm = mean(lap3_log1p_cpm),
    median_lap3_log1p_cpm = median(lap3_log1p_cpm),
    .groups = "drop"
  )
write_table(state_summary, "gse211376_lap3_detection_state_summary.csv")

group_index <- metadata %>%
  filter(compartment %in% c("Malignant", "TAM")) %>%
  mutate(group_id = paste(patient, cell_state, sep = "||")) %>%
  select(cell_id, patient, cell_state, compartment, group_id)

group_factor <- factor(group_index$group_id, levels = unique(group_index$group_id))
group_cells <- match(group_index$cell_id, cell_ids)
design <- Matrix::sparse.model.matrix(~ 0 + group_factor)
colnames(design) <- sub("^group_factor", "", colnames(design))

cat("Aggregating target-gene counts to patient x state pseudobulk...\n")
pseudobulk_counts <- selected_mat[, group_cells, drop = FALSE] %*% design
pseudobulk_counts <- as.matrix(pseudobulk_counts)
pseudobulk_lib <- as.numeric(rowsum(
  metadata$library_size[group_cells],
  group = group_factor,
  reorder = FALSE
))
group_meta <- group_index %>%
  distinct(group_id, patient, cell_state, compartment) %>%
  mutate(
    n_cells = as.integer(table(group_factor)[group_id]),
    library_size = pseudobulk_lib[match(group_id, levels(group_factor))]
  )
stopifnot(identical(colnames(pseudobulk_counts), group_meta$group_id))

log_cpm <- log1p(sweep(pseudobulk_counts, 2, group_meta$library_size, "/") * 1e6)

eligible_groups <- group_meta$n_cells >= 20
score_matrix <- matrix(
  NA_real_,
  nrow = ncol(log_cpm),
  ncol = length(gene_sets),
  dimnames = list(colnames(log_cpm), names(gene_sets))
)
for (compartment_name in c("Malignant", "TAM")) {
  eligible_compartment <- eligible_groups & group_meta$compartment == compartment_name
  for (pathway in names(gene_sets)) {
    genes <- intersect(gene_sets[[pathway]], rownames(log_cpm))
    x <- log_cpm[genes, eligible_compartment, drop = FALSE]
    z <- t(scale(t(x)))
    z[!is.finite(z)] <- 0
    score_matrix[eligible_compartment, pathway] <- colMeans(z)
  }
}
if (is.null(dim(score_matrix))) {
  score_matrix <- matrix(score_matrix, ncol = 1, dimnames = list(colnames(log_cpm), names(gene_sets)[1]))
}

pseudobulk <- group_meta %>%
  mutate(
    lap3_count = as.numeric(pseudobulk_counts["LAP3", ]),
    lap3_log1p_cpm = as.numeric(log_cpm["LAP3", ]),
    lap3_detected_fraction = cell_qc$lap3_positive_cells[
      match(group_id, paste(cell_qc$patient, cell_qc$cell_state, sep = "||"))
    ] / n_cells
  ) %>%
  bind_cols(as.data.frame(score_matrix)) %>%
  mutate(eligible_for_state_inference = n_cells >= 20)
write_table(pseudobulk, "gse211376_patient_state_pseudobulk_scores.csv")

analysis_pseudobulk <- pseudobulk %>%
  filter(eligible_for_state_inference)

make_aggregate_pseudobulk <- function(compartment_name) {
  keep_cells <- which(metadata$compartment == compartment_name)
  patient_factor <- factor(metadata$patient[keep_cells], levels = unique(metadata$patient))
  patient_design <- Matrix::sparse.model.matrix(~ 0 + patient_factor)
  colnames(patient_design) <- sub("^patient_factor", "", colnames(patient_design))
  agg_counts <- selected_mat[, keep_cells, drop = FALSE] %*% patient_design
  agg_counts <- as.matrix(agg_counts)
  agg_lib <- as.numeric(rowsum(
    metadata$library_size[keep_cells],
    group = patient_factor,
    reorder = FALSE
  ))
  agg_log_cpm <- log1p(sweep(agg_counts, 2, agg_lib, "/") * 1e6)
  agg_scores <- sapply(names(gene_sets), function(pathway) {
    genes <- intersect(gene_sets[[pathway]], rownames(agg_log_cpm))
    x <- agg_log_cpm[genes, , drop = FALSE]
    z <- t(scale(t(x)))
    z[!is.finite(z)] <- 0
    colMeans(z)
  })
  data.frame(
    compartment = compartment_name,
    patient = colnames(agg_counts),
    n_cells = as.integer(table(patient_factor)[colnames(agg_counts)]),
    library_size = agg_lib,
    lap3_count = as.numeric(agg_counts["LAP3", ]),
    lap3_log1p_cpm = as.numeric(agg_log_cpm["LAP3", ]),
    agg_scores,
    check.names = FALSE
  )
}

patient_aggregate <- bind_rows(
  make_aggregate_pseudobulk("Malignant"),
  make_aggregate_pseudobulk("TAM")
)
write_table(patient_aggregate, "gse211376_patient_compartment_pseudobulk_scores.csv")

correlation_rows <- list()
row_id <- 1L
for (comp in c("Malignant", "TAM")) {
  d_comp <- analysis_pseudobulk %>% filter(compartment == comp)
  d_patient <- patient_aggregate %>% filter(compartment == comp)
  for (pathway in names(gene_sets)) {
    res_patient <- cor_test_safe(d_patient$lap3_log1p_cpm, d_patient[[pathway]])
    correlation_rows[[row_id]] <- data.frame(
      compartment = comp,
      analysis = "patient_aggregate",
      cell_state = "All",
      pathway = pathway,
      n = res_patient["n"],
      spearman_rho = res_patient["rho"],
      ci_low = res_patient["ci_low"],
      ci_high = res_patient["ci_high"],
      p_value = res_patient["p_value"]
    )
    row_id <- row_id + 1L

    for (state in sort(unique(d_comp$cell_state))) {
      d_state <- d_comp %>% filter(cell_state == state)
      res_state <- cor_test_safe(d_state$lap3_log1p_cpm, d_state[[pathway]])
      correlation_rows[[row_id]] <- data.frame(
        compartment = comp,
        analysis = "within_state_across_patients",
        cell_state = state,
        pathway = pathway,
          n = res_state["n"],
          spearman_rho = res_state["rho"],
          ci_low = res_state["ci_low"],
          ci_high = res_state["ci_high"],
          p_value = res_state["p_value"]
      )
      row_id <- row_id + 1L
    }

    complete_patients <- d_comp %>%
      count(patient) %>%
      filter(n >= 2) %>%
      pull(patient)
    d_residual <- d_comp %>% filter(patient %in% complete_patients)
    fit_x <- lm(lap3_log1p_cpm ~ patient + cell_state, data = d_residual)
    fit_y <- lm(reformulate(c("patient", "cell_state"), response = pathway), data = d_residual)
    res_adjusted <- cor_test_safe(residuals(fit_x), residuals(fit_y))
    correlation_rows[[row_id]] <- data.frame(
      compartment = comp,
      analysis = "patient_and_state_residual",
      cell_state = "All",
      pathway = pathway,
      n = res_adjusted["n"],
      spearman_rho = res_adjusted["rho"],
      ci_low = res_adjusted["ci_low"],
      ci_high = res_adjusted["ci_high"],
      p_value = res_adjusted["p_value"]
    )
    row_id <- row_id + 1L
  }
}

correlations <- bind_rows(correlation_rows) %>%
  group_by(compartment, analysis) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
  ungroup()
write_table(correlations, "gse211376_lap3_pathway_correlations.csv")

malignant_high_low <- analysis_pseudobulk %>%
  filter(compartment == "Malignant") %>%
  group_by(cell_state) %>%
  mutate(
    lap3_group_within_state = ifelse(
      lap3_log1p_cpm >= median(lap3_log1p_cpm),
      "High", "Low"
    )
  ) %>%
  ungroup()

high_low_tests <- bind_rows(lapply(names(gene_sets), function(pathway) {
  bind_rows(lapply(malignant_states, function(state) {
    d <- malignant_high_low %>% filter(cell_state == state)
    high <- d[[pathway]][d$lap3_group_within_state == "High"]
    low <- d[[pathway]][d$lap3_group_within_state == "Low"]
    test <- suppressWarnings(wilcox.test(high, low, exact = FALSE))
    data.frame(
      cell_state = state,
      pathway = pathway,
      n_high = length(high),
      n_low = length(low),
      median_high = median(high),
      median_low = median(low),
      median_delta = median(high) - median(low),
      p_value = test$p.value
    )
  }))
})) %>%
  group_by(cell_state) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
  ungroup()
write_table(high_low_tests, "gse211376_malignant_high_low_sensitivity.csv")

focus_pathways <- c(
  "HALLMARK_MTORC1_SIGNALING",
  "LEUCINE_BCAA_CORE",
  "MTORC1_READOUT_CORE",
  "REACTOME_TRANSLATION"
)
plot_cor <- correlations %>%
  filter(
    compartment == "Malignant",
    analysis %in% c("patient_aggregate", "patient_and_state_residual"),
    pathway %in% focus_pathways
  ) %>%
  mutate(
    analysis = recode(
      analysis,
      patient_aggregate = "Patient aggregate",
      patient_and_state_residual = "Patient + state residual"
    ),
    pathway = recode(
      pathway,
      HALLMARK_MTORC1_SIGNALING = "mTORC1",
      LEUCINE_BCAA_CORE = "Leucine/BCAA core",
      MTORC1_READOUT_CORE = "mTORC1 readout",
      REACTOME_TRANSLATION = "Translation"
    ),
    label = paste0(sprintf("%.2f", spearman_rho), "\nP=", fmt_p(p_value))
  )

p_a <- ggplot(plot_cor, aes(pathway, analysis, fill = spearman_rho)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = label), size = 3) +
  scale_fill_gradient2(
    low = "#2F6F9F", mid = "white", high = "#C44E52",
    midpoint = 0, limits = c(-1, 1), name = "Spearman rho"
  ) +
  labs(
    title = "GSE211376 malignant-cell LAP3 pathway associations",
    x = NULL, y = NULL
  ) +
  theme_classic(base_size = 9, base_family = "Arial") +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

plot_detection <- state_summary %>%
  filter(compartment %in% c("Malignant", "TAM")) %>%
  mutate(cell_state = reorder(cell_state, pct_lap3_positive))

p_b <- ggplot(plot_detection, aes(pct_lap3_positive, cell_state, fill = compartment)) +
  geom_col(width = 0.65) +
  scale_fill_manual(values = c(Malignant = "#C44E52", TAM = "#4B9B8F")) +
  labs(
    title = "LAP3 detection by state",
    x = "LAP3-positive cells (%)", y = NULL, fill = NULL
  ) +
  theme_classic(base_size = 9, base_family = "Arial") +
  theme(legend.position = "top")

plot_high_low <- high_low_tests %>%
  filter(pathway %in% focus_pathways) %>%
  mutate(
    pathway = recode(
      pathway,
      HALLMARK_MTORC1_SIGNALING = "mTORC1",
      LEUCINE_BCAA_CORE = "Leucine/BCAA core",
      MTORC1_READOUT_CORE = "mTORC1 readout",
      REACTOME_TRANSLATION = "Translation"
    )
  )

p_c <- ggplot(plot_high_low, aes(pathway, cell_state, fill = median_delta)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", median_delta)), size = 3) +
  scale_fill_gradient2(
    low = "#2F6F9F", mid = "white", high = "#C44E52",
    midpoint = 0, name = "Median delta"
  ) +
  labs(
    title = "Within-state high-low score delta",
    x = NULL, y = NULL
  ) +
  theme_classic(base_size = 9, base_family = "Arial") +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

fig <- p_a / (p_b | p_c) +
  plot_layout(heights = c(0.85, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(face = "bold", size = 11),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

ggsave(
  file.path(plot_dir, "gse211376_lap3_malignant_pathway_summary.pdf"),
  fig, width = 183 / 25.4, height = 135 / 25.4, device = cairo_pdf
)
ggsave(
  file.path(plot_dir, "gse211376_lap3_malignant_pathway_summary.png"),
  fig, width = 183 / 25.4, height = 135 / 25.4, dpi = 300,
  bg = "white"
)

cat("\nKey malignant correlations:\n")
print(
  correlations %>%
    filter(
      compartment == "Malignant",
      analysis %in% c("patient_aggregate", "patient_and_state_residual"),
      pathway %in% focus_pathways
    ) %>%
    select(analysis, pathway, n, spearman_rho, p_value, p_adj_BH)
)
cat("Completed:", format(Sys.time()), "\n")
