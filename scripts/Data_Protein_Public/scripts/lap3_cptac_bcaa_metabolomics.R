#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(httr)
  library(jsonlite)
  library(openxlsx)
})

data.table::setDTthreads(8)
options(timeout = 300)
set.seed(20260630)

project_dir <- "/home/lzb/glioma"
table_s2 <- file.path(
  project_dir, "Data_Protein_Public", "data", "CPTAC_GBM_Publication",
  "Wang2021_CancerCell_Table_S2_mmc3.xlsx"
)
matched_file <- file.path(
  project_dir, "Data_Protein_Public", "results", "LAP3_CPTAC_Multiomics_Feasibility",
  "tables", "cptac_lap3_mrna_protein_matched_values.csv"
)
table_s2_lap3_cache <- file.path(
  project_dir, "Data_Protein_Public", "data", "CPTAC_GBM_Publication",
  "lap3_table_s2_lightweight.rds"
)
result_dir <- file.path(
  project_dir, "Data_Protein_Public", "results", "LAP3_CPTAC_BCAA_Metabolomics"
)
table_dir <- file.path(result_dir, "tables")
log_dir <- file.path(result_dir, "logs")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(
  file.exists(table_s2),
  file.exists(matched_file),
  file.exists(table_s2_lap3_cache)
)

log_file <- file.path(log_dir, "lap3_cptac_bcaa_metabolomics.log")
if (file.exists(log_file)) invisible(file.remove(log_file))
log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste(...))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

post_json <- function(url, body) {
  response <- POST(url, content_type_json(), body = body, encode = "json", timeout(180))
  stop_for_status(response)
  fromJSON(content(response, "text", encoding = "UTF-8"), flatten = TRUE)
}

run_cor <- function(dat, exposure, outcome, stratum) {
  x <- dat[[exposure]]
  y <- dat[[outcome]]
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 10L) return(NULL)
  test <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  data.table(
    source = unique(dat$analysis_source),
    stratum = stratum,
    exposure = exposure,
    metabolite = outcome,
    n_complete = sum(ok),
    spearman_rho = unname(test$estimate),
    p_value = test$p.value
  )
}

log_msg("Reading metabolome_normalized sheet")
metabolome_wide <- as.data.table(read.xlsx(
  table_s2,
  sheet = "metabolome_normalized",
  check.names = FALSE
))
stopifnot(names(metabolome_wide)[1] == "Metabolite")

target_map <- c(
  dl_isoleucine = "DL-isoleucine",
  l_leucine = "L-leucine",
  l_valine = "L-valine"
)
target_rows <- metabolome_wide[Metabolite %chin% unname(target_map)]
stopifnot(nrow(target_rows) == length(target_map))

target_long <- melt(
  target_rows,
  id.vars = "Metabolite",
  variable.name = "sampleId",
  value.name = "abundance"
)
target_long[, metabolite_key := names(target_map)[match(Metabolite, target_map)]]
target_wide <- dcast(target_long, sampleId ~ metabolite_key, value.var = "abundance")
for (nm in names(target_map)) {
  target_wide[, (paste0(nm, "_z")) := as.numeric(scale(get(nm)))]
}
target_wide[, bcaa_composite := rowMeans(.SD, na.rm = TRUE), .SDcols = paste0(names(target_map), "_z")]

coverage <- target_long[, .(
  n_samples = uniqueN(sampleId),
  n_nonmissing = sum(is.finite(abundance)),
  min_abundance = min(abundance, na.rm = TRUE),
  median_abundance = median(abundance, na.rm = TRUE),
  max_abundance = max(abundance, na.rm = TRUE)
), by = .(metabolite_key, Metabolite)]

log_msg("Fetching IDH1/IDH2 mutation status")
mutation_url <- paste0(
  "https://www.cbioportal.org/api/molecular-profiles/",
  "gbm_cptac_2021_mutations/mutations/fetch?projection=DETAILED"
)
mutations <- as.data.table(post_json(
  mutation_url,
  list(
    entrezGeneIds = list(3417L, 3418L),
    sampleListId = "gbm_cptac_2021_all"
  )
))
idh_mut_samples <- unique(mutations$sampleId)
mutation_summary <- if (nrow(mutations) > 0L) {
  mutations[, .(
    genes = paste(sort(unique(gene.hugoGeneSymbol)), collapse = ";"),
    protein_changes = paste(sort(unique(proteinChange)), collapse = ";")
  ), by = .(sampleId, patientId)]
} else {
  data.table(sampleId = character(), patientId = character())
}

matched <- fread(matched_file)
analysis_data <- merge(matched, target_wide, by = "sampleId", all = FALSE)
analysis_data[, analysis_source := "cBioPortal_matched"]
analysis_data[, idh_mutation_status := fifelse(
  sampleId %chin% idh_mut_samples,
  "IDH1_or_IDH2_mutant",
  "no_IDH1_or_IDH2_mutation_detected"
)]
analysis_data <- merge(
  analysis_data,
  mutation_summary,
  by = c("sampleId", "patientId"),
  all.x = TRUE
)
stopifnot(!anyDuplicated(analysis_data$sampleId))

outcomes <- c(names(target_map), "bcaa_composite")
exposures <- c("lap3_mrna", "lap3_protein")
table_s2_cache <- readRDS(table_s2_lap3_cache)
extract_lap3_row <- function(x, value_name) {
  sample_cols <- grep("^C3[NL]-", names(x), value = TRUE)
  out <- melt(
    as.data.table(x)[, ..sample_cols],
    measure.vars = sample_cols,
    variable.name = "sampleId",
    value.name = value_name
  )
  out[, (value_name) := as.numeric(get(value_name))]
  out
}
table_s2_mrna <- extract_lap3_row(table_s2_cache$mrna, "lap3_mrna")
table_s2_protein <- extract_lap3_row(table_s2_cache$protein, "lap3_protein")
table_s2_analysis <- merge(
  table_s2_mrna,
  table_s2_protein,
  by = "sampleId",
  all = FALSE
)
table_s2_analysis <- merge(
  table_s2_analysis,
  target_wide,
  by = "sampleId",
  all = FALSE
)
table_s2_analysis[, `:=`(
  patientId = sampleId,
  analysis_source = "TableS2_workbook_consistent",
  idh_mutation_status = fifelse(
    sampleId %chin% idh_mut_samples,
    "IDH1_or_IDH2_mutant",
    "not_used_for_IDH_sensitivity"
  )
)]

strata <- list(
  `cBioPortal matched / All` = analysis_data,
  `cBioPortal matched / IDH1/2-wildtype sensitivity` = analysis_data[
    idh_mutation_status == "no_IDH1_or_IDH2_mutation_detected"
  ],
  `Table S2 workbook / All tumors` = table_s2_analysis
)
cor_results <- rbindlist(lapply(names(strata), function(stratum) {
  rbindlist(lapply(exposures, function(exposure) {
    rbindlist(lapply(outcomes, function(outcome) {
      run_cor(strata[[stratum]], exposure, outcome, stratum)
    }))
  }))
}))
cor_results[, fdr := p.adjust(p_value, method = "BH"), by = .(source, stratum, exposure)]
setorder(cor_results, stratum, exposure, fdr, metabolite)

fwrite(coverage, file.path(table_dir, "cptac_bcaa_metabolite_coverage.csv"))
fwrite(
  target_long,
  file.path(table_dir, "cptac_bcaa_metabolite_values_long.csv")
)
fwrite(
  analysis_data,
  file.path(table_dir, "cptac_lap3_bcaa_matched_analysis_data.csv")
)
fwrite(
  table_s2_analysis,
  file.path(table_dir, "cptac_lap3_bcaa_table_s2_analysis_data.csv")
)
fwrite(
  mutation_summary,
  file.path(table_dir, "cptac_idh1_idh2_mutation_samples.csv")
)
fwrite(
  cor_results,
  file.path(table_dir, "cptac_lap3_bcaa_correlations.csv")
)

best_rows <- cor_results[metabolite == "bcaa_composite"]
readme <- c(
  "# CPTAC GBM LAP3-BCAA metabolomics analysis",
  "",
  "## Data and design",
  "",
  "- Source: Wang et al. 2021 Cancer Cell Table S2 (`mmc3.xlsx`).",
  sprintf(
    "- Metabolome sheet: %d metabolites across %d samples; log2 transformed and global-median normalized.",
    nrow(metabolome_wide), ncol(metabolome_wide) - 1L
  ),
  "- Direct BCAA features: L-leucine, DL-isoleucine, and L-valine.",
  "- LAP3 exposures: matched mRNA and total protein from the same CPTAC GBM study.",
  "- Table S2 contains 80 tumor metabolome samples; 75 have complete LAP3 mRNA, LAP3 protein, and BCAA measurements.",
  "- The 75 complete cases are identical through direct Table S2 extraction and the cBioPortal interface.",
  "- IDH sensitivity: 69 cBioPortal-matched samples without any IDH1 or IDH2 mutation in the study mutation profile.",
  "- Multiplicity: BH correction across the three metabolites and the prespecified BCAA composite within each exposure and stratum.",
  "",
  "## BCAA composite result",
  "",
  paste0(
    "- ", best_rows$stratum, " / ", best_rows$exposure,
    ": n=", best_rows$n_complete,
    ", Spearman rho=", sprintf("%.3f", best_rows$spearman_rho),
    ", FDR=", format(best_rows$fdr, digits = 3, scientific = TRUE), "."
  ),
  "",
  "## Interpretation boundary",
  "",
  "- This is direct matched metabolite evidence, stronger than transcript-derived BCAA scores.",
  "- Cross-sectional correlation does not establish that LAP3 enzymatic activity generates intracellular leucine or activates mTORC1.",
  "- Bulk tumor metabolite abundance may still reflect tissue composition, substrate availability, and other metabolic enzymes.",
  "- Wet-lab LAP3 perturbation with intracellular leucine/BCAA measurement and rescue remains required for causality.",
  "",
  "## Outputs",
  "",
  "- `tables/cptac_bcaa_metabolite_coverage.csv`",
  "- `tables/cptac_bcaa_metabolite_values_long.csv`",
  "- `tables/cptac_lap3_bcaa_matched_analysis_data.csv`",
  "- `tables/cptac_lap3_bcaa_table_s2_analysis_data.csv`",
  "- `tables/cptac_idh1_idh2_mutation_samples.csv`",
  "- `tables/cptac_lap3_bcaa_correlations.csv`"
)
writeLines(readme, file.path(result_dir, "README.md"))
log_msg(
  "Completed:", nrow(analysis_data), "matched metabolome samples;",
  length(idh_mut_samples), "IDH1/2-mutant samples"
)
