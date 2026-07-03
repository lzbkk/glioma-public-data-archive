#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(httr)
  library(jsonlite)
})

data.table::setDTthreads(8)
options(timeout = 300)

project_dir <- "/home/lzb/glioma"
protein_dir <- file.path(
  project_dir, "Data_Protein_Public", "results", "LAP3_Protein_Readout", "tables"
)
result_dir <- file.path(
  project_dir, "Data_Protein_Public", "results", "LAP3_CPTAC_Multiomics_Feasibility"
)
table_dir <- file.path(result_dir, "tables")
log_dir <- file.path(result_dir, "logs")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

protein_file <- file.path(protein_dir, "cptac_gbm_target_total_protein_values.csv")
stopifnot(file.exists(protein_file))
log_file <- file.path(log_dir, "lap3_cptac_multiomics_feasibility.log")
if (file.exists(log_file)) invisible(file.remove(log_file))
log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste(...))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

get_json <- function(url) {
  response <- GET(url, timeout(180))
  stop_for_status(response)
  fromJSON(content(response, "text", encoding = "UTF-8"), flatten = TRUE)
}

post_json <- function(url, body) {
  response <- POST(url, content_type_json(), body = body, encode = "json", timeout(180))
  stop_for_status(response)
  fromJSON(content(response, "text", encoding = "UTF-8"), flatten = TRUE)
}

cbio_base <- "https://www.cbioportal.org/api"
study_id <- "gbm_cptac_2021"
sample_list_id <- "gbm_cptac_2021_all"
mrna_profile_id <- "gbm_cptac_2021_mrna"
lap3_entrez <- 51056L

log_msg("Querying CPTAC GBM profile inventory")
profiles <- as.data.table(get_json(sprintf(
  "%s/studies/%s/molecular-profiles", cbio_base, study_id
)))
profile_keep <- intersect(
  c(
    "molecularProfileId", "molecularAlterationType", "genericAssayType",
    "datatype", "name", "description", "patientLevel"
  ),
  names(profiles)
)
fwrite(
  profiles[, ..profile_keep],
  file.path(table_dir, "cptac_gbm_molecular_profile_inventory.csv")
)

log_msg("Fetching LAP3 mRNA")
mrna <- as.data.table(post_json(
  sprintf(
    "%s/molecular-profiles/%s/molecular-data/fetch?projection=DETAILED",
    cbio_base, mrna_profile_id
  ),
  list(entrezGeneIds = list(lap3_entrez), sampleListId = sample_list_id)
))
mrna[, lap3_mrna := as.numeric(value)]
mrna <- mrna[, .(sampleId, patientId, lap3_mrna)]

protein <- fread(protein_file)
lap3_protein <- protein[gene_symbol == "LAP3", .(
  sampleId, patientId, lap3_protein = value
)]
stopifnot(!anyDuplicated(mrna$sampleId), !anyDuplicated(lap3_protein$sampleId))

matched <- merge(
  mrna,
  lap3_protein,
  by = c("sampleId", "patientId"),
  all = TRUE
)
matched[, `:=`(
  mrna_available = is.finite(lap3_mrna),
  protein_available = is.finite(lap3_protein),
  matched_available = is.finite(lap3_mrna) & is.finite(lap3_protein)
)]
fwrite(matched, file.path(table_dir, "cptac_lap3_mrna_protein_matched_values.csv"))

complete <- matched[matched_available == TRUE]
stopifnot(nrow(complete) >= 10L)
spearman <- cor.test(
  complete$lap3_mrna,
  complete$lap3_protein,
  method = "spearman",
  exact = FALSE
)
pearson <- cor.test(
  complete$lap3_mrna,
  complete$lap3_protein,
  method = "pearson"
)
cor_result <- data.table(
  gene = "LAP3",
  n_mrna = sum(matched$mrna_available),
  n_protein = sum(matched$protein_available),
  n_matched = nrow(complete),
  spearman_rho = unname(spearman$estimate),
  spearman_p = spearman$p.value,
  pearson_r = unname(pearson$estimate),
  pearson_p = pearson$p.value,
  mrna_profile = mrna_profile_id,
  protein_profile = "gbm_cptac_2021_protein_quantification"
)
fwrite(cor_result, file.path(table_dir, "cptac_lap3_mrna_protein_correlation.csv"))

metabolomics_audit <- data.table(
  source = c(
    "CPTAC GBM study publication",
    "cBioPortal gbm_cptac_2021",
    "LinkedOmics CPTAC-GBM",
    "CPTAC Python current GBM sources",
    "Publication processed-data supplement"
  ),
  polar_metabolome_status = c(
    "generated",
    "not_exposed",
    "not_listed",
    "not_exposed_in_current_loader",
    "candidate_entry"
  ),
  evidence = c(
    "GC-MS polar metabolomics was performed on the 99-treatment-naive-GBM cohort",
    "Profile inventory contains lipidome but no metabolome generic assay",
    "Download inventory lists RNA/protein/phosphoprotein but no metabolome",
    "Current multi-source GBM loader exposes genomics/transcriptomics/proteomics sources",
    "Paper states processed data are available in Table S2; metabolite sheet and sample IDs require direct inspection"
  ),
  leucine_bcaa_direct_analysis = c(
    "potentially_possible",
    "not_possible_from_this_entry",
    "not_possible_from_this_entry",
    "not_confirmed",
    "pending_file_access_and_feature_audit"
  )
)
fwrite(
  metabolomics_audit,
  file.path(table_dir, "cptac_gbm_metabolomics_entry_audit.csv")
)

availability <- data.table(
  data_layer = c(
    "LAP3 mRNA",
    "LAP3 total protein",
    "matched LAP3 mRNA-protein",
    "phosphoproteome",
    "lipidome",
    "polar metabolome",
    "leucine/isoleucine/valine"
  ),
  current_status = c(
    "available",
    "available",
    "available",
    "available",
    "available",
    "exists_in_study_but_matrix_not_yet_located",
    "not_yet_audited_in_feature_matrix"
  ),
  n_samples_or_scope = c(
    as.character(sum(matched$mrna_available)),
    as.character(sum(matched$protein_available)),
    as.character(nrow(complete)),
    "99 with site-dependent coverage",
    "positive and negative mode profiles",
    "reported for the same study cohort",
    "requires metabolite feature table"
  ),
  decision = c(
    "GO",
    "GO",
    "GO",
    "completed_first_pass",
    "out_of_scope_for_direct_BCAA",
    "CONDITIONAL",
    "CONDITIONAL"
  )
)
fwrite(availability, file.path(table_dir, "cptac_multiomics_availability_summary.csv"))

readme <- c(
  "# CPTAC GBM LAP3 multi-omics feasibility audit",
  "",
  "## Matched mRNA-protein result",
  "",
  sprintf(
    "- LAP3 mRNA: %d samples; LAP3 total protein: %d samples; exact matched pairs: %d.",
    cor_result$n_mrna, cor_result$n_protein, cor_result$n_matched
  ),
  sprintf(
    "- LAP3 mRNA-protein Spearman rho = %.3f (P = %.3g); Pearson r = %.3f (P = %.3g).",
    cor_result$spearman_rho, cor_result$spearman_p,
    cor_result$pearson_r, cor_result$pearson_p
  ),
  "",
  "## Metabolomics feasibility",
  "",
  "- The source publication confirms GC-MS polar metabolomics in the 99-treatment-naive-GBM cohort.",
  "- The current cBioPortal study exposes mRNA, total protein, phosphoproteome, acetylproteome, and lipidome, but no polar-metabolome profile.",
  "- LinkedOmics and the current CPTAC Python GBM loader do not expose a metabolome matrix through their standard inventory.",
  "- Follow-up completed: publication Table S2 contains 134 metabolites across 87 samples, including direct L-leucine, DL-isoleucine, and L-valine measurements.",
  "",
  "## Decision",
  "",
  "- Matched LAP3 mRNA-protein: GO and completed.",
  "- Direct leucine/BCAA metabolomics: GO and completed in `LAP3_CPTAC_BCAA_Metabolomics`.",
  "- Do not infer metabolite abundance from transcript or lipidome features.",
  "",
  "## Outputs",
  "",
  "- `tables/cptac_gbm_molecular_profile_inventory.csv`",
  "- `tables/cptac_lap3_mrna_protein_matched_values.csv`",
  "- `tables/cptac_lap3_mrna_protein_correlation.csv`",
  "- `tables/cptac_gbm_metabolomics_entry_audit.csv`",
  "- `tables/cptac_multiomics_availability_summary.csv`"
)
writeLines(readme, file.path(result_dir, "README.md"))
log_msg("Completed matched mRNA-protein and metabolomics entry audit")
