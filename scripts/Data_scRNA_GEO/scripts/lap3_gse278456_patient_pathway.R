#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(msigdbr)
  library(ggplot2)
  library(patchwork)
})

setwd("/home/lzb/glioma")
data.table::setDTthreads(8)
set.seed(20260629)

object_file <- "Data_scRNA_GEO/GSE278456/Tumor_Integrated_SeuratV5.rds"
crosswalk_file <- "Data_scRNA_GEO/results/GSE278456_Tumor_Object_Audit/exports/gse278456_local_author_metadata_crosswalk.rds"
out_dir <- "Data_scRNA_GEO/results/GSE278456_LAP3_Patient_Pathway"
table_dir <- file.path(out_dir, "tables")
plot_dir <- file.path(out_dir, "plots")
log_dir <- file.path(out_dir, "logs")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(log_dir, "lap3_gse278456_patient_pathway.log")
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

fmt_p <- function(p) {
  ifelse(
    is.na(p), "NA",
    ifelse(p < 0.001, "<0.001", formatC(p, format = "f", digits = 2))
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
cat("Loading local object...\n")
obj <- readRDS(object_file)
crosswalk <- readRDS(crosswalk_file)
stopifnot(ncol(obj) == nrow(crosswalk), identical(colnames(obj), crosswalk$cell_id_local))

crosswalk <- crosswalk %>%
  mutate(
    pathology = author_Pathology,
    tumor_grade = author_Tumor.Grade,
    idh_status = author_IDH.status,
    analysis_class = case_when(
      pathology == "GBM" & tumor_grade == "IV" & idh_status == "wild type" ~ "GBM_grade4_IDHwt",
      pathology %in% c("Astrocytoma", "Anaplastic Astrocytoma", "Oligodendroglioma") &
        tumor_grade %in% c("II", "III") & idh_status == "mutant" ~ "LGG_grade2_3_IDHmut",
      pathology == "Normal" | tumor_grade == "Normal" ~ "Normal",
      TRUE ~ "Other_tumor"
    )
  )

matched_cells <- which(crosswalk$matched_author)
stopifnot(length(matched_cells) == 193730)

clinical_fields <- c(
  "Patient_ID", "pathology", "tumor_grade", "idh_status", "analysis_class",
  "author_EGFR.amplification", "author_EGFR.mutation", "author_PTEN.mutation",
  "author_TP53.mutation", "author_CDKN2A..2G.loss"
)
patient_annotation <- crosswalk[matched_cells, clinical_fields] %>%
  group_by(Patient_ID) %>%
  summarise(
    across(
      everything(),
      ~ {
        values <- unique(na.omit(as.character(.x)))
        if (length(values) == 0) NA_character_ else paste(values, collapse = ";")
      }
    ),
    n_matched_cells = n(),
    .groups = "drop"
  )
write_table(patient_annotation, "gse278456_patient_annotation.csv")

ambiguity <- patient_annotation %>%
  pivot_longer(
    cols = -c(Patient_ID, n_matched_cells),
    names_to = "field",
    values_to = "value"
  ) %>%
  mutate(ambiguous = grepl(";", value, fixed = TRUE))
write_table(ambiguity, "gse278456_patient_annotation_ambiguity.csv")
stopifnot(!any(ambiguity$ambiguous, na.rm = TRUE))

gene_sets <- make_gene_sets()
target_genes <- unique(c("LAP3", unlist(gene_sets, use.names = FALSE)))
counts <- LayerData(obj[["RNA"]], layer = "counts", fast = FALSE)
genes_present <- rownames(counts)

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
write_table(coverage, "gse278456_pathway_gene_coverage.csv")
stopifnot(all(coverage$genes_present >= 10), "LAP3" %in% genes_present)

selected_counts <- counts[intersect(target_genes, genes_present), matched_cells, drop = FALSE]
patient_factor <- factor(
  crosswalk$Patient_ID[matched_cells],
  levels = unique(crosswalk$Patient_ID[matched_cells])
)
design <- sparse.model.matrix(~ 0 + patient_factor)
colnames(design) <- sub("^patient_factor", "", colnames(design))
patient_counts <- as.matrix(selected_counts %*% design)

library_sizes <- as.numeric(rowsum(
  obj$nCount_RNA[matched_cells],
  group = patient_factor,
  reorder = FALSE
))
names(library_sizes) <- levels(patient_factor)
stopifnot(identical(colnames(patient_counts), names(library_sizes)))
log_cpm <- log1p(sweep(patient_counts, 2, library_sizes, "/") * 1e6)

score_matrix <- sapply(names(gene_sets), function(pathway) {
  genes <- intersect(gene_sets[[pathway]], rownames(log_cpm))
  x <- log_cpm[genes, , drop = FALSE]
  z <- t(scale(t(x)))
  z[!is.finite(z)] <- 0
  colMeans(z)
})

patient_scores <- patient_annotation %>%
  transmute(
    patient = Patient_ID,
    pathology,
    tumor_grade,
    idh_status,
    analysis_class,
    n_matched_cells,
    library_size = library_sizes[Patient_ID],
    lap3_count = as.numeric(patient_counts["LAP3", Patient_ID]),
    lap3_log1p_cpm = as.numeric(log_cpm["LAP3", Patient_ID])
  ) %>%
  bind_cols(as.data.frame(score_matrix[match(patient_annotation$Patient_ID, rownames(score_matrix)), , drop = FALSE]))

if (is.null(rownames(score_matrix))) {
  patient_scores <- patient_annotation %>%
    transmute(
      patient = Patient_ID,
      pathology,
      tumor_grade,
      idh_status,
      analysis_class,
      n_matched_cells,
      library_size = library_sizes[Patient_ID],
      lap3_count = as.numeric(patient_counts["LAP3", Patient_ID]),
      lap3_log1p_cpm = as.numeric(log_cpm["LAP3", Patient_ID])
    ) %>%
    bind_cols(as.data.frame(score_matrix))
}
write_table(patient_scores, "gse278456_patient_pseudobulk_scores.csv")

analysis_subsets <- list(
  all_non_normal_tumor = patient_scores %>%
    filter(analysis_class != "Normal"),
  GBM_grade4_IDHwt = patient_scores %>%
    filter(analysis_class == "GBM_grade4_IDHwt"),
  LGG_grade2_3_IDHmut = patient_scores %>%
    filter(analysis_class == "LGG_grade2_3_IDHmut")
)

score_subset <- function(data, pathway) {
  patients <- data$patient
  genes <- intersect(gene_sets[[pathway]], rownames(log_cpm))
  x <- log_cpm[genes, patients, drop = FALSE]
  z <- t(scale(t(x)))
  z[!is.finite(z)] <- 0
  colMeans(z)
}

correlations <- bind_rows(lapply(names(analysis_subsets), function(subset_name) {
  d <- analysis_subsets[[subset_name]]
  bind_rows(lapply(names(gene_sets), function(pathway) {
    subset_score <- score_subset(d, pathway)
    result <- cor_test_safe(d$lap3_log1p_cpm, subset_score)
    data.frame(
      analysis_subset = subset_name,
      pathway = pathway,
      n_patients = result["n"],
      spearman_rho = result["rho"],
      ci_low = result["ci_low"],
      ci_high = result["ci_high"],
      p_value = result["p_value"]
    )
  }))
})) %>%
  group_by(analysis_subset) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
  ungroup()
write_table(correlations, "gse278456_lap3_pathway_correlations.csv")

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
    subset_label = recode(
      analysis_subset,
      all_non_normal_tumor = "All tumors",
      GBM_grade4_IDHwt = "GBM G4 / IDH-wt",
      LGG_grade2_3_IDHmut = "LGG G2/3 / IDH-mut"
    ),
    label = paste0(sprintf("%.2f", spearman_rho), "\nP", fmt_p(p_value))
  )

p_a <- ggplot(plot_df, aes(pathway_label, subset_label, fill = spearman_rho)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = label), size = 2.5) +
  scale_fill_gradient2(
    low = "#2F6F9F", mid = "white", high = "#C44E52",
    midpoint = 0, limits = c(-1, 1), name = "Spearman rho"
  ) +
  labs(
    title = "Patient-level pathway associations",
    x = NULL, y = NULL
  ) +
  theme_classic(base_size = 9, base_family = "Arial") +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

p_b <- patient_scores %>%
  count(analysis_class, name = "n_patients") %>%
  ggplot(aes(n_patients, reorder(analysis_class, n_patients), fill = analysis_class)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = n_patients), hjust = -0.15, size = 3) +
  scale_x_continuous(
    limits = c(0, max(patient_scores %>% count(analysis_class) %>% pull(n)) * 1.25),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(title = "Patient strata", x = "Patients", y = NULL) +
  theme_classic(base_size = 9, base_family = "Arial") +
  theme(
    legend.position = "none",
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

fig <- (p_a | p_b) +
  plot_layout(widths = c(1.65, 0.75)) +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(face = "bold", size = 11),
    plot.background = element_rect(fill = "white", colour = NA)
  )

ggsave(
  file.path(plot_dir, "gse278456_patient_pathway_summary.pdf"),
  fig, width = 183 / 25.4, height = 105 / 25.4, device = cairo_pdf
)
ggsave(
  file.path(plot_dir, "gse278456_patient_pathway_summary.png"),
  fig, width = 183 / 25.4, height = 105 / 25.4, dpi = 300, bg = "white"
)

cat("\nPatient strata:\n")
print(patient_scores %>% count(analysis_class))
cat("\nFocus correlations:\n")
print(correlations %>% filter(pathway %in% focus))
cat("Completed:", format(Sys.time()), "\n")
