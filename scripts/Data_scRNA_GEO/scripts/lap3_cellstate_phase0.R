#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readxl)
  library(msigdbr)
})

setwd("/home/lzb/glioma")
data.table::setDTthreads(8)
set.seed(20260629)

out_dir <- "Data_scRNA_GEO/results/LAP3_CellState_Phase0"
table_dir <- file.path(out_dir, "tables")
source_dir <- file.path(out_dir, "source_data")
log_dir <- file.path(out_dir, "logs")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(log_dir, "lap3_cellstate_phase0.log")
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

cat("Started:", format(Sys.time()), "\n")
cat("R:", R.version.string, "\n")
cat("data.table threads:", data.table::getDTthreads(), "\n")

# 1. Dataset and patient provenance
gse211_metadata_file <- paste0(
  "Data_scRNA_GEO/GSE211376/",
  "GSE211376_metadata_Ruiz2022_all_samples_filtered_cells.csv.gz"
)
gse131_metadata_file <- paste0(
  "Data_scRNA_GEO/GSE131928/",
  "GSE131928_single_cells_tumor_name_and_adult_or_peidatric.xlsx"
)

gse211 <- read.csv(
  gzfile(gse211_metadata_file),
  stringsAsFactors = FALSE,
  check.names = FALSE
) %>%
  transmute(
    dataset = "GSE211376",
    patient = as.character(patient),
    platform = "snRNA-seq",
    age_group = NA_character_
  ) %>%
  distinct()

gse131 <- read_excel(gse131_metadata_file, skip = 43) %>%
  transmute(
    dataset = "GSE131928",
    patient = as.character(`tumour name`),
    platform = ifelse(
      grepl("Smartseq2", `processed data file`),
      "Smartseq2",
      "10X"
    ),
    age_group = as.character(`adult/pediatric`)
  ) %>%
  distinct()

patient_platform <- bind_rows(gse211, gse131) %>%
  arrange(dataset, patient, platform)
write_table(patient_platform, "dataset_patient_platform_inventory.csv")

gse131_overlap <- gse131 %>%
  group_by(patient) %>%
  summarise(
    n_platforms = n_distinct(platform),
    platforms = paste(sort(unique(platform)), collapse = ";"),
    age_group = paste(sort(unique(age_group)), collapse = ";"),
    .groups = "drop"
  ) %>%
  mutate(overlap_between_platforms = n_platforms > 1L)
write_table(gse131_overlap, "gse131928_patient_platform_overlap.csv")

cross_dataset_overlap <- tidyr::crossing(
  gse211_patient = sort(unique(gse211$patient)),
  gse131_patient = sort(unique(gse131$patient))
) %>%
  mutate(
    exact_match = gse211_patient == gse131_patient,
    normalized_match = gsub("[^A-Z0-9]", "", toupper(gse211_patient)) ==
      gsub("[^A-Z0-9]", "", toupper(gse131_patient))
  ) %>%
  filter(exact_match | normalized_match)
write_table(cross_dataset_overlap, "gse211376_gse131928_patient_id_overlap.csv")

provenance <- data.frame(
  dataset = c("GSE211376", "GSE131928", "GSE278456"),
  local_input = c(
    "Ruiz2022_all_samples filtered counts and metadata; 11 NH patients",
    "Smart-seq2 and 10X processed TPM plus GEO metadata workbook",
    "Author TumorSeurat metadata crosswalk plus local _Tumor/_Tumor2 count object"
  ),
  phase0_finding = c(
    "Local files are the 11 newly generated NH samples, not the complete core/extended GBmap matrix; no patient ID overlap with GSE131928",
    "Six patients occur on both Smart-seq2 and 10X; platform-pooled inference must retain patient clustering",
    "Tumor-labelled input is enriched/selected but lacks author malignant-state labels; malignancy and contamination require expression QC before state scoring"
  ),
  evidence_status = c(
    "supported_by_filename_patient_ids_and_GBmap_data_availability",
    "confirmed_from_GEO_cell_metadata",
    "confirmed_from_local_object_audit_file_labels; biological purity pending_QC"
  ),
  source = c(
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC12526130/",
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131928",
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE278456"
  ),
  stringsAsFactors = FALSE
)
write_table(provenance, "dataset_provenance_audit.csv")

# 2. Freeze Neftel meta-modules from the original Table S2
neftel_file <- file.path(source_dir, "Neftel_2019_Table_S2_meta_modules.xlsx")
stopifnot(file.exists(neftel_file))
neftel_wide <- read_excel(neftel_file, sheet = "Table S2", skip = 4)
expected_modules <- c("MES2", "MES1", "AC", "OPC", "NPC1", "NPC2", "G1/S", "G2/M")
stopifnot(identical(names(neftel_wide), expected_modules))

neftel_long <- bind_rows(lapply(expected_modules, function(module_name) {
  genes <- na.omit(as.character(neftel_wide[[module_name]]))
  genes <- genes[nzchar(genes)]
  data.frame(
    signature_family = ifelse(module_name %in% c("G1/S", "G2/M"), "cell_cycle", "Neftel_state"),
    signature = module_name,
    rank = seq_along(genes),
    gene = genes,
    source = "Neftel et al. Cell 2019 Table S2",
    stringsAsFactors = FALSE
  )
}))
stopifnot(
  !anyDuplicated(neftel_long[, c("signature", "gene")]),
  !("LAP3" %in% neftel_long$gene),
  all(c("MES2", "MES1", "AC", "OPC", "NPC1", "NPC2") %in% neftel_long$signature)
)
write_table(neftel_long, "frozen_neftel_meta_modules.csv")

# 3. Freeze pathway sets used by the project
hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>%
  select(gs_name, gene_symbol)
reactome <- msigdbr(
  species = "Homo sapiens",
  category = "C2",
  subcategory = "CP:REACTOME"
) %>%
  select(gs_name, gene_symbol)
msig_sets <- split(
  bind_rows(hallmark, reactome)$gene_symbol,
  bind_rows(hallmark, reactome)$gs_name
)

pathway_sets <- list(
  HALLMARK_MTORC1_SIGNALING = msig_sets[["HALLMARK_MTORC1_SIGNALING"]],
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
  ),
  REACTOME_TRANSLATION = msig_sets[["REACTOME_TRANSLATION"]]
)
pathway_sets <- lapply(pathway_sets, function(genes) {
  sort(setdiff(unique(na.omit(genes)), "LAP3"))
})

pathway_long <- bind_rows(lapply(names(pathway_sets), function(signature_name) {
  family <- if (signature_name %in% c(
    "HALLMARK_MTORC1_SIGNALING",
    "LEUCINE_BCAA_CORE"
  )) {
    "primary"
  } else {
    "secondary"
  }
  data.frame(
    fdr_family = family,
    signature = signature_name,
    gene = pathway_sets[[signature_name]],
    source = ifelse(
      grepl("HALLMARK|REACTOME", signature_name),
      paste0("MSigDB via msigdbr ", as.character(packageVersion("msigdbr"))),
      "project_predefined_manual_set_frozen_20260629"
    ),
    stringsAsFactors = FALSE
  )
}))
stopifnot(
  !("LAP3" %in% pathway_long$gene),
  all(table(pathway_long$signature) >= 10L)
)
write_table(pathway_long, "frozen_lap3_pathway_gene_sets.csv")
saveRDS(
  list(neftel = split(neftel_long$gene, neftel_long$signature), pathways = pathway_sets),
  file.path(source_dir, "frozen_cellstate_gene_sets.rds")
)

gene_set_summary <- bind_rows(
  neftel_long %>% count(signature_family, signature, name = "n_genes"),
  pathway_long %>%
    transmute(signature_family = paste0("pathway_", fdr_family), signature) %>%
    count(signature_family, signature, name = "n_genes")
)
write_table(gene_set_summary, "frozen_gene_set_summary.csv")

fdr_manifest <- data.frame(
  fdr_family = c("primary", "secondary", "exploratory_state"),
  included_tests = c(
    "Pre-specified disease stratum x state x HALLMARK_MTORC1_SIGNALING/LEUCINE_BCAA_CORE",
    "Pre-specified disease stratum x state x MTORC1_READOUT_CORE/REACTOME_TRANSLATION",
    "Six continuous Neftel meta-modules and classification-sensitivity tests"
  ),
  adjustment = "Benjamini-Hochberg within family",
  stringsAsFactors = FALSE
)
write_table(fdr_manifest, "fdr_family_manifest.csv")

# 4. Record package readiness
packages <- c("UCell", "AUCell", "lme4", "nlme", "data.table", "dplyr", "msigdbr", "testthat")
package_status <- bind_rows(lapply(packages, function(package_name) {
  installed <- requireNamespace(package_name, quietly = TRUE)
  data.frame(
    package = package_name,
    installed = installed,
    version = if (installed) as.character(packageVersion(package_name)) else NA_character_,
    note = if (package_name == "lme4") {
      "Installed but Matrix ABI warning observed; do not use until compatibility is restored"
    } else {
      ""
    },
    stringsAsFactors = FALSE
  )
}))
write_table(package_status, "r_package_readiness.csv")

cat("GSE211376 patients:", n_distinct(gse211$patient), "\n")
cat("GSE131928 unique patients:", n_distinct(gse131$patient), "\n")
cat("GSE131928 cross-platform overlaps:", sum(gse131_overlap$overlap_between_platforms), "\n")
cat("Cross-dataset patient ID overlaps:", nrow(cross_dataset_overlap), "\n")
print(gene_set_summary)
cat("Completed:", format(Sys.time()), "\n")
