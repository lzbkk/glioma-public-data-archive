#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(readxl)
  library(msigdbr)
  library(ggplot2)
  library(patchwork)
})

setwd("/home/lzb/glioma")
data.table::setDTthreads(8)
set.seed(20260629)

smartseq_file <- "Data_scRNA_GEO/GSE131928/GSE131928_RAW/GSM3828672_Smartseq2_GBM_IDHwt_processed_TPM.tsv.gz"
tenx_file <- "Data_scRNA_GEO/GSE131928/GSE131928_RAW/GSM3828673_10X_GBM_IDHwt_processed_TPM.tsv.gz"
metadata_file <- "Data_scRNA_GEO/GSE131928/GSE131928_single_cells_tumor_name_and_adult_or_peidatric.xlsx"
out_dir <- "Data_scRNA_GEO/results/GSE131928_LAP3_Patient_Pathway"
table_dir <- file.path(out_dir, "tables")
plot_dir <- file.path(out_dir, "plots")
log_dir <- file.path(out_dir, "logs")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(log_dir, "lap3_gse131928_patient_pathway.log")
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message")
on.exit({
  sink(type = "message")
  sink(type = "output")
  close(log_con)
}, add = TRUE)

write_table <- function(x, filename) {
  fwrite(x, file.path(table_dir, filename))
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

cat("Started:", format(Sys.time()), "\n")

metadata <- read_excel(metadata_file, skip = 43) %>%
  rename(
    cell_id = `Sample name`,
    patient = `tumour name`,
    age_group = `adult/pediatric`,
    processed_file = `processed data file`
  ) %>%
  mutate(
    platform = ifelse(grepl("Smartseq2", processed_file), "Smartseq2", "10X")
  )

stopifnot(
  nrow(metadata) == 24131,
  !anyDuplicated(metadata$cell_id),
  all(c("adult", "pediatric") %in% unique(metadata$age_group))
)
write_table(metadata, "gse131928_cell_metadata.csv")

gene_sets <- make_gene_sets()
target_genes <- unique(c("LAP3", unlist(gene_sets, use.names = FALSE)))

aggregate_platform <- function(matrix_file, platform_name) {
  cat("Loading", platform_name, "matrix...\n")
  dt <- fread(
    cmd = paste("zcat", shQuote(matrix_file)),
    header = TRUE,
    sep = "\t",
    data.table = TRUE,
    showProgress = TRUE,
    nThread = 8
  )
  gene_col <- names(dt)[1]
  cell_ids <- names(dt)[-1]
  meta <- metadata %>% filter(platform == platform_name)
  stopifnot(
    length(cell_ids) == nrow(meta),
    setequal(cell_ids, meta$cell_id)
  )
  meta <- meta[match(cell_ids, meta$cell_id), ]

  genes_present <- dt[[gene_col]]
  selected <- dt[get(gene_col) %in% target_genes]
  selected_genes <- selected[[gene_col]]
  expr <- as.matrix(selected[, -1])
  storage.mode(expr) <- "numeric"
  rownames(expr) <- selected_genes
  colnames(expr) <- cell_ids

  patient_cells <- split(seq_along(cell_ids), meta$patient)
  patient_mean <- sapply(patient_cells, function(index) {
    rowMeans(expr[, index, drop = FALSE])
  })
  if (is.null(dim(patient_mean))) {
    patient_mean <- matrix(patient_mean, ncol = 1)
  }

  patient_meta <- bind_rows(lapply(names(patient_cells), function(patient_id) {
    index <- patient_cells[[patient_id]]
    data.frame(
      platform = platform_name,
      patient = patient_id,
      age_group = paste(unique(meta$age_group[index]), collapse = ";"),
      n_cells = length(index)
    )
  }))

  eligible <- patient_meta$n_cells >= 20
  scores <- sapply(names(gene_sets), function(pathway) {
    genes <- intersect(gene_sets[[pathway]], rownames(patient_mean))
    x <- patient_mean[genes, eligible, drop = FALSE]
    z <- t(scale(t(x)))
    z[!is.finite(z)] <- 0
    score <- rep(NA_real_, ncol(patient_mean))
    score[eligible] <- colMeans(z)
    score
  })

  patient_scores <- patient_meta %>%
    mutate(
      eligible = eligible,
      lap3_mean_processed_tpm = as.numeric(patient_mean["LAP3", patient]),
      lap3_detection_fraction = vapply(patient_cells, function(index) {
        mean(expr["LAP3", index] > 0)
      }, numeric(1))[patient]
    ) %>%
    bind_cols(as.data.frame(scores))

  coverage <- bind_rows(lapply(names(gene_sets), function(pathway) {
    requested <- gene_sets[[pathway]]
    present <- intersect(requested, genes_present)
    data.frame(
      platform = platform_name,
      pathway = pathway,
      genes_requested = length(requested),
      genes_present = length(present),
      coverage = length(present) / length(requested),
      missing_genes = paste(setdiff(requested, present), collapse = ";")
    )
  }))

  rm(dt, selected, expr)
  invisible(gc())
  list(scores = patient_scores, coverage = coverage, patient_mean = patient_mean)
}

smartseq <- aggregate_platform(smartseq_file, "Smartseq2")
tenx <- aggregate_platform(tenx_file, "10X")
patient_scores <- bind_rows(smartseq$scores, tenx$scores)
coverage <- bind_rows(smartseq$coverage, tenx$coverage)

write_table(patient_scores, "gse131928_patient_mean_scores.csv")
write_table(coverage, "gse131928_pathway_gene_coverage.csv")
stopifnot(all(coverage$genes_present >= 10))

analysis_sets <- list(
  Smartseq2_all = list(
    data = patient_scores %>% filter(platform == "Smartseq2", eligible),
    matrix = smartseq$patient_mean
  ),
  Smartseq2_adult = list(
    data = patient_scores %>% filter(platform == "Smartseq2", eligible, age_group == "adult"),
    matrix = smartseq$patient_mean
  ),
  `10X_all` = list(
    data = patient_scores %>% filter(platform == "10X", eligible),
    matrix = tenx$patient_mean
  ),
  `10X_adult` = list(
    data = patient_scores %>% filter(platform == "10X", eligible, age_group == "adult"),
    matrix = tenx$patient_mean
  )
)

correlations <- bind_rows(lapply(names(analysis_sets), function(analysis_name) {
  d <- analysis_sets[[analysis_name]]$data
  patient_mean <- analysis_sets[[analysis_name]]$matrix
  bind_rows(lapply(names(gene_sets), function(pathway) {
    genes <- intersect(gene_sets[[pathway]], rownames(patient_mean))
    x <- patient_mean[genes, d$patient, drop = FALSE]
    z <- t(scale(t(x)))
    z[!is.finite(z)] <- 0
    subset_score <- colMeans(z)
    result <- cor_test_safe(patient_mean["LAP3", d$patient], subset_score)
    data.frame(
      analysis = analysis_name,
      pathway = pathway,
      n_patients = result["n"],
      spearman_rho = result["rho"],
      ci_low = result["ci_low"],
      ci_high = result["ci_high"],
      p_value = result["p_value"]
    )
  }))
})) %>%
  group_by(analysis) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
  ungroup()

pooled_rows <- bind_rows(lapply(names(gene_sets), function(pathway) {
  pooled_values <- bind_rows(lapply(
    list(Smartseq2 = smartseq, `10X` = tenx),
    function(platform_data) {
      d <- platform_data$scores %>% filter(eligible)
      patient_mean <- platform_data$patient_mean
      genes <- intersect(gene_sets[[pathway]], rownames(patient_mean))
      x <- patient_mean[genes, d$patient, drop = FALSE]
      z <- t(scale(t(x)))
      z[!is.finite(z)] <- 0
      data.frame(
        lap3_platform_z = as.numeric(scale(patient_mean["LAP3", d$patient])),
        pathway_score = colMeans(z)
      )
    }
  ))
  result <- cor_test_safe(pooled_values$lap3_platform_z, pooled_values$pathway_score)
  data.frame(
    analysis = "platform_residual_pooled",
    pathway = pathway,
    n_patients = result["n"],
    spearman_rho = result["rho"],
    ci_low = result["ci_low"],
    ci_high = result["ci_high"],
    p_value = result["p_value"]
  )
})) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH"))

correlations <- bind_rows(correlations, pooled_rows)
write_table(correlations, "gse131928_lap3_pathway_correlations.csv")

focus <- c(
  "HALLMARK_MTORC1_SIGNALING",
  "LEUCINE_BCAA_CORE",
  "MTORC1_READOUT_CORE",
  "REACTOME_TRANSLATION"
)
plot_df <- correlations %>%
  filter(pathway %in% focus) %>%
  mutate(
    pathway_label = recode(
      pathway,
      HALLMARK_MTORC1_SIGNALING = "mTORC1",
      LEUCINE_BCAA_CORE = "Leucine/BCAA core",
      MTORC1_READOUT_CORE = "mTORC1 readout",
      REACTOME_TRANSLATION = "Translation"
    ),
    analysis_label = recode(
      analysis,
      Smartseq2_all = "Smart-seq2 all",
      Smartseq2_adult = "Smart-seq2 adult",
      `10X_all` = "10x all",
      `10X_adult` = "10x adult",
      platform_residual_pooled = "Platform-adjusted pooled"
    ),
    label = paste0(sprintf("%.2f", spearman_rho), "\nP=", formatC(p_value, format = "g", digits = 2))
  )

p_a <- ggplot(plot_df, aes(pathway_label, analysis_label, fill = spearman_rho)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = label), size = 2.5) +
  scale_fill_gradient2(
    low = "#2F6F9F", mid = "white", high = "#C44E52",
    midpoint = 0, limits = c(-1, 1), name = "Spearman rho"
  ) +
  labs(title = "GSE131928 patient-level validation", x = NULL, y = NULL) +
  theme_classic(base_size = 9, base_family = "Arial") +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

p_b <- patient_scores %>%
  filter(eligible) %>%
  count(platform, age_group, name = "n_patients") %>%
  ggplot(aes(platform, n_patients, fill = age_group)) +
  geom_col(position = "dodge", width = 0.65) +
  geom_text(
    aes(label = n_patients),
    position = position_dodge(width = 0.65),
    vjust = -0.3,
    size = 3
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Patient groups", x = NULL, y = "Patients", fill = NULL) +
  theme_classic(base_size = 9, base_family = "Arial") +
  theme(
    legend.position = "top",
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

fig <- (p_a | p_b) +
  plot_layout(widths = c(1.6, 0.7)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 11))

ggsave(
  file.path(plot_dir, "gse131928_patient_pathway_summary.pdf"),
  fig, width = 183 / 25.4, height = 115 / 25.4, device = cairo_pdf
)
ggsave(
  file.path(plot_dir, "gse131928_patient_pathway_summary.png"),
  fig, width = 183 / 25.4, height = 115 / 25.4, dpi = 300, bg = "white"
)

cat("\nPatient counts:\n")
print(patient_scores %>% count(platform, age_group, eligible))
cat("\nFocus correlations:\n")
print(correlations %>% filter(pathway %in% focus))
cat("Completed:", format(Sys.time()), "\n")
