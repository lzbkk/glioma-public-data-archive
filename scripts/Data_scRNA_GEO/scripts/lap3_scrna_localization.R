#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

set.seed(20260628)
data.table::setDTthreads(8)

project_dir <- "/home/lzb/glioma"
data_dir <- file.path(project_dir, "Data_scRNA_GEO")
out_dir <- file.path(data_dir, "results", "LAP3_scRNA_Localization")
tables_dir <- file.path(out_dir, "tables")
plots_dir <- file.path(out_dir, "plots")
logs_dir <- file.path(out_dir, "logs")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

run_log <- file.path(logs_dir, "lap3_scrna_localization.log")
if (file.exists(run_log)) {
  unlink(run_log)
}

log_msg <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  cat(msg, "\n", file = run_log, append = TRUE)
}

log_msg("LAP3 scRNA localization started: ", Sys.time())
log_msg("R version: ", getRversion())

safe_shell <- function(cmd) {
  out <- system(cmd, intern = TRUE)
  if (!is.null(attr(out, "status"))) {
    stop("Command failed: ", cmd)
  }
  out
}

read_annotation <- function(path, modality) {
  dt <- fread(path, header = FALSE, col.names = c("cell_id", "cell_type"))
  dt[, modality := modality]
  dt[, sample_id := sub("_.*$", "", cell_id)]
  dt
}

tar_member_list <- function(tar_path) {
  safe_shell(sprintf("tar -tf %s", shQuote(tar_path)))
}

extract_tar_member <- function(tar_path, member, out_path) {
  cmd <- sprintf("tar -xOf %s %s > %s", shQuote(tar_path), shQuote(member), shQuote(out_path))
  status <- system(cmd)
  if (status != 0) {
    stop("Failed to extract tar member: ", member)
  }
  out_path
}

sample_id_from_member <- function(member) {
  basename(member) |>
    sub("_features.tsv.gz$", "", x = _) |>
    sub("_genes.tsv.gz$", "", x = _) |>
    sub("_barcodes.tsv.gz$", "", x = _) |>
    sub("_matrix.mtx.gz$", "", x = _) |>
    sub("^GSM[0-9]+_", "", x = _)
}

read_10x_lap3_from_tar <- function(tar_path, annotation) {
  members <- tar_member_list(tar_path)
  feature_members <- members[grepl("_(features|genes)\\.tsv\\.gz$", members)]
  feature_members <- feature_members[!grepl("snATAC|peaks", feature_members)]

  tmp_dir <- tempfile("gse138794_10x_")
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  all_rows <- list()
  qc_rows <- list()

  for (feature_member in feature_members) {
    sample_id <- sample_id_from_member(feature_member)
    barcode_member <- sub("_(features|genes)\\.tsv\\.gz$", "_barcodes.tsv.gz", feature_member)
    matrix_member <- sub("_(features|genes)\\.tsv\\.gz$", "_matrix.mtx.gz", feature_member)

    if (!barcode_member %in% members || !matrix_member %in% members) {
      qc_rows[[sample_id]] <- data.table(
        dataset = "GSE138794",
        sample_id = sample_id,
        status = "missing barcode or matrix member"
      )
      next
    }

    feature_path <- file.path(tmp_dir, basename(feature_member))
    barcode_path <- file.path(tmp_dir, basename(barcode_member))
    matrix_path <- file.path(tmp_dir, basename(matrix_member))
    extract_tar_member(tar_path, feature_member, feature_path)
    extract_tar_member(tar_path, barcode_member, barcode_path)
    extract_tar_member(tar_path, matrix_member, matrix_path)

    features <- fread(feature_path, header = FALSE)
    gene_symbol <- if (ncol(features) >= 2) features[[2]] else features[[1]]
    lap3_idx <- which(gene_symbol == "LAP3")
    barcodes <- fread(barcode_path, header = FALSE)[[1]]

    if (length(lap3_idx) == 0) {
      qc_rows[[sample_id]] <- data.table(
        dataset = "GSE138794",
        sample_id = sample_id,
        status = "LAP3 not found"
      )
      next
    }

    mat <- Matrix::readMM(gzfile(matrix_path))
    if (ncol(mat) != length(barcodes)) {
      stop("Barcode/matrix dimension mismatch for ", sample_id)
    }

    lap3_count <- Matrix::colSums(mat[lap3_idx, , drop = FALSE])
    library_size <- Matrix::colSums(mat)
    cell_id <- paste(sample_id, barcodes, sep = "_")

    sample_dt <- data.table(
      dataset = "GSE138794",
      sample_id = sample_id,
      cell_id = cell_id,
      lap3_count = as.numeric(lap3_count),
      library_size = as.numeric(library_size)
    )
    sample_dt <- merge(sample_dt, annotation, by = "cell_id", all.x = FALSE, all.y = FALSE)
    if ("sample_id.x" %in% names(sample_dt)) {
      sample_dt[, sample_id := sample_id.x]
      sample_dt[, c("sample_id.x", "sample_id.y") := NULL]
    }
    sample_dt[, lap3_cpm := fifelse(library_size > 0, lap3_count / library_size * 1e6, NA_real_)]
    sample_dt[, lap3_log1p_cpm := log1p(lap3_cpm)]
    sample_dt[, lap3_detected := lap3_count > 0]

    all_rows[[sample_id]] <- sample_dt
    status <- if (nrow(sample_dt) > 0) "ok" else "no matched annotation"
    qc_rows[[sample_id]] <- data.table(
      dataset = "GSE138794",
      sample_id = sample_id,
      modality = paste(sort(unique(sample_dt$modality)), collapse = ";"),
      matrix_cells = length(barcodes),
      annotated_cells = nrow(sample_dt),
      lap3_positive_cells = sum(sample_dt$lap3_detected),
      status = status
    )
  }

  list(
    cells = rbindlist(all_rows, fill = TRUE),
    qc = rbindlist(qc_rows, fill = TRUE)
  )
}

read_lap3_from_wide_matrix <- function(matrix_gz, metadata_gz, dataset) {
  metadata <- fread(metadata_gz, header = TRUE)
  if (ncol(metadata) < 3) {
    metadata <- fread(metadata_gz, header = FALSE)
    setnames(metadata, c("cell_id", "cell_type", "patient"))
  } else {
    setnames(metadata, names(metadata)[1:3], c("cell_id", "cell_type", "patient"))
  }
  metadata <- metadata[, .(cell_id, cell_type, patient)]
  metadata[, modality := "scRNA"]
  metadata[, sample_id := patient]

  header_line <- readLines(gzfile(matrix_gz), n = 1)
  lap3_line <- safe_shell(sprintf("zgrep -m 1 '^LAP3\t' %s", shQuote(matrix_gz)))
  if (length(lap3_line) == 0) {
    stop("LAP3 row not found in ", matrix_gz)
  }

  cells <- strsplit(header_line, "\t", fixed = TRUE)[[1]]
  values <- strsplit(lap3_line, "\t", fixed = TRUE)[[1]]
  gene <- values[[1]]
  if (gene != "LAP3") {
    stop("Unexpected gene row while reading ", matrix_gz)
  }
  counts <- suppressWarnings(as.numeric(values[-1]))
  if (length(cells) != length(counts)) {
    stop("Header/value length mismatch for ", dataset)
  }

  dt <- data.table(
    dataset = dataset,
    cell_id = cells,
    lap3_count = counts
  )
  dt <- merge(dt, metadata, by = "cell_id", all.x = FALSE, all.y = FALSE)
  dt[, library_size := NA_real_]
  dt[, lap3_cpm := NA_real_]
  dt[, lap3_log1p_cpm := log1p(lap3_count)]
  dt[, lap3_detected := lap3_count > 0]
  dt
}

marker_sets <- list(
  MYELOID_MACROPHAGE = c(
    "AIF1", "C1QA", "C1QB", "C1QC", "CD14", "CD68", "CSF1R",
    "FCER1G", "FCGR3A", "ITGAM", "LST1", "MS4A7", "TYROBP"
  ),
  FOAM_LIPID_TAF = c(
    "APOE", "APOC1", "LPL", "TREM2", "GPNMB", "SPP1", "FABP5",
    "PLIN2", "LIPA", "CD36", "LGALS3", "CTSB", "CTSL"
  ),
  MALIGNANT_GLIOMA = c(
    "EGFR", "OLIG2", "PDGFRA", "SOX2", "NES", "GFAP", "VIM", "MKI67"
  ),
  T_NK = c("CD3D", "CD3E", "TRAC", "NKG7", "GNLY"),
  ENDOTHELIAL = c("PECAM1", "VWF", "KDR", "CLDN5")
)

read_marker_panel_from_10x_tar <- function(tar_path, dataset) {
  members <- tar_member_list(tar_path)
  feature_members <- members[grepl("_(features|genes)\\.tsv\\.gz$", members)]
  target_genes <- sort(unique(c("LAP3", unlist(marker_sets, use.names = FALSE))))

  tmp_dir <- tempfile(paste0(tolower(dataset), "_10x_"))
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  all_rows <- list()
  qc_rows <- list()
  coverage_rows <- list()

  for (feature_member in feature_members) {
    sample_id <- sample_id_from_member(feature_member)
    gsm_id <- sub("_.*$", "", basename(feature_member))
    patient_id <- sub(".*_Tumour_", "", sample_id)
    patient_id <- sub("_S[0-9]+$", "", patient_id)
    barcode_member <- sub("_(features|genes)\\.tsv\\.gz$", "_barcodes.tsv.gz", feature_member)
    matrix_member <- sub("_(features|genes)\\.tsv\\.gz$", "_matrix.mtx.gz", feature_member)

    if (!barcode_member %in% members || !matrix_member %in% members) {
      qc_rows[[sample_id]] <- data.table(
        dataset = dataset,
        sample_id = sample_id,
        gsm_id = gsm_id,
        patient_id = patient_id,
        status = "missing barcode or matrix member"
      )
      next
    }

    feature_path <- file.path(tmp_dir, basename(feature_member))
    barcode_path <- file.path(tmp_dir, basename(barcode_member))
    matrix_path <- file.path(tmp_dir, basename(matrix_member))
    extract_tar_member(tar_path, feature_member, feature_path)
    extract_tar_member(tar_path, barcode_member, barcode_path)
    extract_tar_member(tar_path, matrix_member, matrix_path)

    features <- fread(feature_path, header = FALSE)
    gene_symbol <- if (ncol(features) >= 2) features[[2]] else features[[1]]
    keep_idx <- which(gene_symbol %in% target_genes)
    barcodes <- fread(barcode_path, header = FALSE)[[1]]
    mat <- Matrix::readMM(gzfile(matrix_path))
    if (ncol(mat) != length(barcodes)) {
      stop("Barcode/matrix dimension mismatch for ", sample_id)
    }

    library_size <- Matrix::colSums(mat)
    cell_id <- paste(sample_id, barcodes, sep = "_")
    sample_dt <- data.table(
      dataset = dataset,
      sample_id = sample_id,
      gsm_id = gsm_id,
      patient_id = patient_id,
      cell_id = cell_id,
      library_size = as.numeric(library_size)
    )

    gene_counts <- matrix(0, nrow = length(target_genes), ncol = ncol(mat))
    rownames(gene_counts) <- target_genes
    if (length(keep_idx) > 0) {
      detected_genes <- gene_symbol[keep_idx]
      count_mat <- as.matrix(mat[keep_idx, , drop = FALSE])
      for (gene in unique(detected_genes)) {
        gene_counts[gene, ] <- Matrix::colSums(count_mat[detected_genes == gene, , drop = FALSE])
      }
    }

    add_gene_metric <- function(dt, gene) {
      counts <- as.numeric(gene_counts[gene, ])
      set(dt, j = paste0(gene, "_count"), value = counts)
      set(dt, j = paste0(gene, "_log1p_cpm"), value = log1p(fifelse(library_size > 0, counts / library_size * 1e6, NA_real_)))
      dt
    }
    sample_dt <- add_gene_metric(sample_dt, "LAP3")
    setnames(sample_dt, c("LAP3_count", "LAP3_log1p_cpm"), c("lap3_count", "lap3_log1p_cpm"))
    set(sample_dt, j = "lap3_detected", value = sample_dt$lap3_count > 0)

    for (set_name in names(marker_sets)) {
      genes <- intersect(marker_sets[[set_name]], rownames(gene_counts))
      present_genes <- intersect(genes, gene_symbol)
      if (length(genes) == 0) {
        set(sample_dt, j = paste0(tolower(set_name), "_score"), value = NA_real_)
        set(sample_dt, j = paste0(tolower(set_name), "_detected_genes"), value = NA_integer_)
        next
      }
      cpm_mat <- t(t(gene_counts[genes, , drop = FALSE]) / library_size * 1e6)
      log_mat <- log1p(cpm_mat)
      score_col <- paste0(tolower(set_name), "_score")
      detected_col <- paste0(tolower(set_name), "_detected_genes")
      set(sample_dt, j = score_col, value = colMeans(log_mat, na.rm = TRUE))
      set(sample_dt, j = detected_col, value = Matrix::colSums(gene_counts[genes, , drop = FALSE] > 0))
      coverage_rows[[paste(sample_id, set_name, sep = "__")]] <- data.table(
        dataset = dataset,
        sample_id = sample_id,
        marker_set = set_name,
        n_marker_genes = length(genes),
        n_present_in_features = length(present_genes),
        present_genes = paste(sort(present_genes), collapse = ";")
      )
    }

    marker_compartment <- fifelse(
      sample_dt$myeloid_macrophage_detected_genes >= 2 & sample_dt$foam_lipid_taf_detected_genes >= 2,
      "myeloid_foam_like",
      fifelse(
        sample_dt$myeloid_macrophage_detected_genes >= 2,
        "myeloid_marker_high",
        fifelse(
          sample_dt$foam_lipid_taf_detected_genes >= 2,
          "foam_lipid_marker_high",
          fifelse(sample_dt$malignant_glioma_detected_genes >= 2, "glioma_marker_high", "other_or_low_marker")
        )
      )
    )
    set(sample_dt, j = "marker_compartment", value = marker_compartment)

    all_rows[[sample_id]] <- sample_dt
    qc_rows[[sample_id]] <- data.table(
      dataset = dataset,
      sample_id = sample_id,
      gsm_id = gsm_id,
      patient_id = patient_id,
      matrix_cells = length(barcodes),
      matrix_genes = nrow(mat),
      lap3_positive_cells = sum(sample_dt$lap3_detected),
      status = "ok"
    )
  }

  list(
    cells = rbindlist(all_rows, fill = TRUE),
    qc = rbindlist(qc_rows, fill = TRUE),
    coverage = rbindlist(coverage_rows, fill = TRUE)
  )
}

summarize_gse237673_marker_compartments <- function(cells) {
  cells |>
    as_tibble() |>
    group_by(dataset, marker_compartment) |>
    summarise(
      n_cells = n(),
      n_samples = n_distinct(sample_id),
      lap3_positive_cells = sum(lap3_detected, na.rm = TRUE),
      pct_lap3_positive = 100 * mean(lap3_detected, na.rm = TRUE),
      mean_lap3_log1p_cpm = mean(lap3_log1p_cpm, na.rm = TRUE),
      median_lap3_log1p_cpm = median(lap3_log1p_cpm, na.rm = TRUE),
      mean_myeloid_score = mean(myeloid_macrophage_score, na.rm = TRUE),
      mean_foam_lipid_score = mean(foam_lipid_taf_score, na.rm = TRUE),
      mean_glioma_score = mean(malignant_glioma_score, na.rm = TRUE),
      .groups = "drop"
    ) |>
    arrange(desc(mean_lap3_log1p_cpm), desc(pct_lap3_positive)) |>
    as.data.table()
}

summarize_gse237673_sample_compartments <- function(cells) {
  cells |>
    as_tibble() |>
    group_by(dataset, sample_id, patient_id, marker_compartment) |>
    summarise(
      n_cells = n(),
      lap3_positive_cells = sum(lap3_detected, na.rm = TRUE),
      pct_lap3_positive = 100 * mean(lap3_detected, na.rm = TRUE),
      mean_lap3_log1p_cpm = mean(lap3_log1p_cpm, na.rm = TRUE),
      mean_myeloid_score = mean(myeloid_macrophage_score, na.rm = TRUE),
      mean_foam_lipid_score = mean(foam_lipid_taf_score, na.rm = TRUE),
      .groups = "drop"
    ) |>
    arrange(sample_id, marker_compartment) |>
    as.data.table()
}

parse_gse318564_sample_metadata <- function(tar_path) {
  members <- tar_member_list(tar_path)
  archives <- members[grepl("_filtered_feature_bc_matrix\\.tar\\.gz$", members)]
  if (length(archives) == 0) {
    return(data.table())
  }
  dt <- data.table(archive_member = archives)
  dt[, archive_file := basename(archive_member)]
  dt[, gsm_id := sub("_.*$", "", archive_file)]
  dt[, sample_label := sub("^GSM[0-9]+_", "", archive_file)]
  dt[, sample_label := sub("_filtered_feature_bc_matrix\\.tar\\.gz$", "", sample_label)]
  dt[, region := fifelse(grepl("core", sample_label, ignore.case = TRUE), "core",
                         fifelse(grepl("edge", sample_label, ignore.case = TRUE), "edge", "unknown"))]
  dt[, cd45_selection := fifelse(grepl("CD45neg", sample_label, ignore.case = TRUE), "CD45neg", "not_specified")]
  dt[, patient_id := fifelse(
    grepl("^(WU-[0-9]+)", sample_label),
    sub("^((WU-[0-9]+)).*$", "\\1", sample_label),
    fifelse(
      grepl("^(GBM[0-9]+)", sample_label),
      sub("^((GBM[0-9]+)).*$", "\\1", sample_label),
      sub("_.*$", "", sample_label)
    )
  )]
  dt[, sample_index := fifelse(grepl("_[0-9]+$", sample_label), sub("^.*_([0-9]+)$", "\\1", sample_label), NA_character_)]
  dt[, geo_sample_type := fifelse(
    patient_id %in% c("B178", "B183"),
    "same-patient_replicate",
    fifelse(patient_id == "B189", "single_sample", "region_sample")
  )]
  dt[, region_confidence := fifelse(region %in% c("core", "edge"), "geo_title_confirmed", "not_provided_by_geo")]
  dt[, include_core_edge := region %in% c("core", "edge")]
  dt[, metadata_source := "NCBI GEO SOFT accessed 2026-06-28"]
  dt[, note := fifelse(
    region == "unknown",
    fifelse(
      patient_id %in% c("B178", "B183"),
      "GEO describes _1/_2/_3 as replicates of one GBM sample; no core/edge label is provided.",
      "GEO provides no core/edge label; retain for overall validation but exclude from regional comparison."
    ),
    "Core/edge label is explicitly present in the GEO sample title."
  )]
  dt[]
}

read_gse318564_marker_panel <- function(tar_path, metadata) {
  target_genes <- sort(unique(c("LAP3", unlist(marker_sets, use.names = FALSE))))
  tmp_dir <- tempfile("gse318564_nested_")
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  all_rows <- list()
  qc_rows <- list()
  coverage_rows <- list()

  for (i in seq_len(nrow(metadata))) {
    meta <- metadata[i]
    nested_path <- file.path(tmp_dir, meta$archive_file)
    extract_tar_member(tar_path, meta$archive_member, nested_path)
    nested_members <- safe_shell(sprintf("tar -tzf %s", shQuote(nested_path)))
    feature_member <- nested_members[grepl("(^|/)(features|genes)\\.tsv\\.gz$", nested_members)][1]
    barcode_member <- nested_members[grepl("(^|/)barcodes\\.tsv\\.gz$", nested_members)][1]
    matrix_member <- nested_members[grepl("(^|/)matrix\\.mtx\\.gz$", nested_members)][1]

    if (any(is.na(c(feature_member, barcode_member, matrix_member)))) {
      qc_rows[[meta$sample_label]] <- data.table(
        dataset = "GSE318564",
        sample_id = meta$sample_label,
        patient_id = meta$patient_id,
        region = meta$region,
        matrix_cells = NA_integer_,
        matrix_genes = NA_integer_,
        lap3_positive_cells = NA_integer_,
        status = "incomplete_archive",
        note = paste0(
          "Missing nested 10x members: ",
          paste(c("features", "barcodes", "matrix")[is.na(c(
            feature_member, barcode_member, matrix_member
          ))], collapse = ";")
        )
      )
      log_msg("GSE318564 skipped incomplete archive: ", meta$sample_label)
      unlink(nested_path)
      next
    }

    sample_dir <- file.path(tmp_dir, meta$sample_label)
    dir.create(sample_dir)
    feature_path <- file.path(sample_dir, "features.tsv.gz")
    barcode_path <- file.path(sample_dir, "barcodes.tsv.gz")
    matrix_path <- file.path(sample_dir, "matrix.mtx.gz")
    extract_nested <- function(member, out_path) {
      status <- system(sprintf(
        "tar -xzOf %s %s > %s",
        shQuote(nested_path), shQuote(member), shQuote(out_path)
      ))
      if (status != 0) stop("Failed to extract nested member: ", member)
    }
    extract_nested(feature_member, feature_path)
    extract_nested(barcode_member, barcode_path)
    extract_nested(matrix_member, matrix_path)

    features <- fread(feature_path, header = FALSE)
    gene_symbol <- if (ncol(features) >= 2) features[[2]] else features[[1]]
    barcodes <- fread(barcode_path, header = FALSE)[[1]]
    mat <- Matrix::readMM(gzfile(matrix_path))
    stopifnot(ncol(mat) == length(barcodes), nrow(mat) == length(gene_symbol))

    keep_idx <- which(gene_symbol %in% target_genes)
    library_size <- Matrix::colSums(mat)
    gene_counts <- matrix(0, nrow = length(target_genes), ncol = ncol(mat),
                          dimnames = list(target_genes, NULL))
    if (length(keep_idx) > 0) {
      count_mat <- mat[keep_idx, , drop = FALSE]
      detected_genes <- gene_symbol[keep_idx]
      for (gene in unique(detected_genes)) {
        gene_counts[gene, ] <- Matrix::colSums(count_mat[detected_genes == gene, , drop = FALSE])
      }
    }

    sample_dt <- data.table(
      dataset = "GSE318564",
      gsm_id = meta$gsm_id,
      sample_id = meta$sample_label,
      patient_id = meta$patient_id,
      region = meta$region,
      cd45_selection = meta$cd45_selection,
      include_core_edge = meta$include_core_edge,
      cell_id = paste(meta$sample_label, barcodes, sep = "_"),
      library_size = as.numeric(library_size)
    )
    lap3_count <- as.numeric(gene_counts["LAP3", ])
    sample_dt[, lap3_count := lap3_count]
    sample_dt[, lap3_log1p_cpm := log1p(fifelse(
      library_size > 0, lap3_count / library_size * 1e6, NA_real_
    ))]
    sample_dt[, lap3_detected := lap3_count > 0]

    for (set_name in names(marker_sets)) {
      genes <- marker_sets[[set_name]]
      cpm_mat <- t(t(gene_counts[genes, , drop = FALSE]) / library_size * 1e6)
      sample_dt[, (paste0(tolower(set_name), "_score")) := colMeans(log1p(cpm_mat), na.rm = TRUE)]
      sample_dt[, (paste0(tolower(set_name), "_detected_genes")) :=
                  Matrix::colSums(gene_counts[genes, , drop = FALSE] > 0)]
      coverage_rows[[paste(meta$sample_label, set_name, sep = "__")]] <- data.table(
        dataset = "GSE318564",
        sample_id = meta$sample_label,
        marker_set = set_name,
        n_marker_genes = length(genes),
        n_present_in_features = sum(genes %in% gene_symbol),
        present_genes = paste(sort(intersect(genes, gene_symbol)), collapse = ";")
      )
    }

    sample_dt[, marker_compartment := fifelse(
      myeloid_macrophage_detected_genes >= 2 & foam_lipid_taf_detected_genes >= 2,
      "myeloid_foam_like",
      fifelse(
        myeloid_macrophage_detected_genes >= 2,
        "myeloid_marker_high",
        fifelse(
          foam_lipid_taf_detected_genes >= 2,
          "foam_lipid_marker_high",
          fifelse(malignant_glioma_detected_genes >= 2, "glioma_marker_high", "other_or_low_marker")
        )
      )
    )]

    all_rows[[meta$sample_label]] <- sample_dt
    qc_rows[[meta$sample_label]] <- data.table(
      dataset = "GSE318564",
      sample_id = meta$sample_label,
      patient_id = meta$patient_id,
      region = meta$region,
      matrix_cells = ncol(mat),
      matrix_genes = nrow(mat),
      lap3_positive_cells = sum(sample_dt$lap3_detected),
      status = "ok"
    )
    rm(mat, gene_counts, sample_dt)
    unlink(c(nested_path, sample_dir), recursive = TRUE)
    gc(verbose = FALSE)
    log_msg("GSE318564 processed sample: ", meta$sample_label)
  }

  list(
    cells = rbindlist(all_rows, fill = TRUE),
    qc = rbindlist(qc_rows, fill = TRUE),
    coverage = rbindlist(coverage_rows, fill = TRUE)
  )
}

summarize_gse318564 <- function(cells) {
  sample_summary <- cells[, .(
    n_cells = .N,
    lap3_positive_cells = sum(lap3_detected),
    pct_lap3_positive = 100 * mean(lap3_detected),
    mean_lap3_log1p_cpm = mean(lap3_log1p_cpm),
    mean_myeloid_score = mean(myeloid_macrophage_score),
    mean_foam_lipid_score = mean(foam_lipid_taf_score),
    mean_glioma_score = mean(malignant_glioma_score)
  ), by = .(dataset, gsm_id, sample_id, patient_id, region, cd45_selection, include_core_edge)]

  compartment_summary <- cells[, .(
    n_cells = .N,
    lap3_positive_cells = sum(lap3_detected),
    pct_lap3_positive = 100 * mean(lap3_detected),
    mean_lap3_log1p_cpm = mean(lap3_log1p_cpm)
  ), by = .(dataset, sample_id, patient_id, region, cd45_selection, marker_compartment)]

  list(sample = sample_summary, compartment = compartment_summary)
}

summarize_by_cell_type <- function(cells) {
  cells |>
    as_tibble() |>
    group_by(dataset, modality, cell_type) |>
    summarise(
      n_cells = n(),
      n_samples = n_distinct(sample_id),
      lap3_positive_cells = sum(lap3_detected, na.rm = TRUE),
      pct_lap3_positive = 100 * mean(lap3_detected, na.rm = TRUE),
      mean_lap3_count = mean(lap3_count, na.rm = TRUE),
      mean_lap3_log1p_cpm = mean(lap3_log1p_cpm, na.rm = TRUE),
      median_lap3_log1p_cpm = median(lap3_log1p_cpm, na.rm = TRUE),
      .groups = "drop"
    ) |>
    arrange(dataset, modality, desc(mean_lap3_log1p_cpm), desc(pct_lap3_positive)) |>
    as.data.table()
}

summarize_by_sample <- function(cells) {
  cells |>
    as_tibble() |>
    group_by(dataset, modality, sample_id, cell_type) |>
    summarise(
      n_cells = n(),
      lap3_positive_cells = sum(lap3_detected, na.rm = TRUE),
      pct_lap3_positive = 100 * mean(lap3_detected, na.rm = TRUE),
      mean_lap3_log1p_cpm = mean(lap3_log1p_cpm, na.rm = TRUE),
      .groups = "drop"
    ) |>
    arrange(dataset, modality, sample_id, cell_type) |>
    as.data.table()
}

write_plot <- function(summary_dt, cell_dt) {
  gse138 <- summary_dt[dataset == "GSE138794"]
  if (nrow(gse138) > 0) {
    cell_type_order <- gse138[
      ,
      .(rank_value = max(mean_lap3_log1p_cpm, na.rm = TRUE)),
      by = cell_type
    ][order(rank_value), cell_type]
    gse138[, cell_type := factor(cell_type, levels = rev(cell_type_order))]
    p <- ggplot(gse138, aes(x = mean_lap3_log1p_cpm, y = cell_type)) +
      geom_point(aes(size = pct_lap3_positive, color = modality), alpha = 0.85) +
      facet_wrap(~ modality, scales = "free_y") +
      scale_size_continuous(range = c(1.5, 7)) +
      labs(
        x = "Mean LAP3 log1p(CPM)",
        y = NULL,
        size = "LAP3+ cells (%)",
        color = "Modality",
        title = "GSE138794 LAP3 cell-type localization"
      ) +
      theme_bw(base_size = 11)
    ggsave(file.path(plots_dir, "gse138794_lap3_cell_type_dotplot.png"), p, width = 8.5, height = 5.5, dpi = 300)
  }

  violin_cells <- cell_dt[dataset == "GSE138794" & !is.na(lap3_log1p_cpm)]
  if (nrow(violin_cells) > 0) {
    top_types <- summary_dt[dataset == "GSE138794"][
      ,
      .(rank_value = max(mean_lap3_log1p_cpm, na.rm = TRUE)),
      by = cell_type
    ][order(-rank_value)][1:min(.N, 12), cell_type]
    plot_dt <- violin_cells[cell_type %in% top_types]
    plot_dt[, cell_type := factor(cell_type, levels = rev(top_types))]
    p2 <- ggplot(plot_dt, aes(x = lap3_log1p_cpm, y = cell_type, fill = modality)) +
      geom_boxplot(outlier.shape = NA, width = 0.55, alpha = 0.72) +
      facet_wrap(~ modality, scales = "free_y") +
      labs(
        x = "LAP3 log1p(CPM)",
        y = NULL,
        fill = "Modality",
        title = "GSE138794 LAP3 expression by annotated cell type"
      ) +
      theme_bw(base_size = 11)
    ggsave(file.path(plots_dir, "gse138794_lap3_cell_type_boxplot.png"), p2, width = 8.5, height = 5.5, dpi = 300)
  }
}

write_gse237673_plot <- function(gse237673_cells, gse237673_summary) {
  if (nrow(gse237673_summary) == 0) {
    return(invisible(NULL))
  }
  plot_dt <- copy(gse237673_summary)
  plot_dt[, marker_compartment := factor(marker_compartment, levels = rev(marker_compartment))]
  p <- ggplot(plot_dt, aes(x = mean_lap3_log1p_cpm, y = marker_compartment)) +
    geom_point(aes(size = pct_lap3_positive, color = mean_foam_lipid_score), alpha = 0.9) +
    scale_size_continuous(range = c(2, 8)) +
    labs(
      x = "Mean LAP3 log1p(CPM)",
      y = NULL,
      size = "LAP3+ cells (%)",
      color = "Foam/lipid score",
      title = "GSE237673 LAP3 in marker-defined myeloid/foam-like compartments"
    ) +
    theme_bw(base_size = 11)
  ggsave(file.path(plots_dir, "gse237673_lap3_marker_compartment_dotplot.png"), p, width = 8.5, height = 4.8, dpi = 300)

  sample_plot <- gse237673_cells[
    ,
    .(
      lap3_log1p_cpm = mean(lap3_log1p_cpm, na.rm = TRUE),
      myeloid_macrophage_score = mean(myeloid_macrophage_score, na.rm = TRUE),
      foam_lipid_taf_score = mean(foam_lipid_taf_score, na.rm = TRUE)
    ),
    by = .(sample_id, marker_compartment)
  ]
  p2 <- ggplot(sample_plot, aes(x = myeloid_macrophage_score, y = foam_lipid_taf_score)) +
    geom_point(aes(size = lap3_log1p_cpm, color = marker_compartment), alpha = 0.85) +
    facet_wrap(~ marker_compartment, scales = "free") +
    labs(
      x = "Mean myeloid/macrophage marker score",
      y = "Mean foam/lipid marker score",
      size = "Mean LAP3 log1p(CPM)",
      color = "Marker compartment",
      title = "GSE237673 sample-level marker scores and LAP3"
    ) +
    theme_bw(base_size = 10)
  ggsave(file.path(plots_dir, "gse237673_marker_scores_vs_lap3.png"), p2, width = 9.5, height = 5.8, dpi = 300)
}

write_readme <- function(summary_dt, dataset_qc, validation_plan) {
  top_gse138 <- summary_dt[dataset == "GSE138794"][order(-mean_lap3_log1p_cpm)][1:min(.N, 8)]
  top_gse211376 <- summary_dt[dataset == "GSE211376"][order(-pct_lap3_positive, -mean_lap3_count)][1:min(.N, 8)]
  gse237673_marker_path <- file.path(tables_dir, "gse237673_lap3_marker_compartment_summary.csv")
  if (file.exists(gse237673_marker_path)) {
    top_gse237673 <- fread(gse237673_marker_path)[order(-mean_lap3_log1p_cpm, -pct_lap3_positive)][1:min(.N, 6)]
  } else {
    top_gse237673 <- data.table()
  }
  if (nrow(top_gse237673) > 0) {
    gse237673_lines <- sprintf(
      "| %s | %s | %s | %s | %.2f | %.3f | %.3f | %.3f |",
      top_gse237673$dataset,
      top_gse237673$marker_compartment,
      top_gse237673$n_cells,
      top_gse237673$n_samples,
      top_gse237673$pct_lap3_positive,
      top_gse237673$mean_lap3_log1p_cpm,
      top_gse237673$mean_myeloid_score,
      top_gse237673$mean_foam_lipid_score
    )
  } else {
    gse237673_lines <- "| NA | NA | 0 | 0 | 0 | 0 | 0 | 0 |"
  }

  lines <- c(
    "# LAP3 scRNA Cell-Type Localization",
    "",
    paste0("生成时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "",
    "## 输入",
    "",
    "```text",
    "Data_scRNA_GEO/GSE138794/GSE138794_RAW.tar",
    "Data_scRNA_GEO/GSE138794/GSE138794_scRNA_Seq_cell_types.txt.gz",
    "Data_scRNA_GEO/GSE138794/GSE138794_snRNA_Seq_cell_types.txt.gz",
    "Data_scRNA_GEO/GSE211376/GSE211376_raw_counts_Ruiz2022_all_samples_filtered_cells.tsv.gz",
    "Data_scRNA_GEO/GSE211376/GSE211376_metadata_Ruiz2022_all_samples_filtered_cells.csv.gz",
    "Data_scRNA_GEO/GSE237673/GSE237673_RAW.tar",
    "Data_scRNA_GEO/GSE318564/GSE318564_RAW.tar",
    "```",
    "",
    "## 方法概述",
    "",
    "- `GSE138794`：从 raw tar 中逐样本抽取 10x matrix，仅计算 `LAP3` raw count、library size、CPM 和 `log1p(CPM)`，再与 scRNA/snRNA cell-type annotation 合并。",
    "- `GSE211376`：作为第一批独立验证资源，读取现成 raw-count wide matrix 的 `LAP3` 行，并与 metadata 中的 cell type 合并。",
    "- `GSE237673`：作为 macrophage/foam-cell focused 验证资源，缺少本地官方 cell annotation，因此读取 marker panel 并输出 marker-defined myeloid/foam-like compartments。",
    "- `GSE318564`：依据 GEO SOFT 校准 sample title，逐个读取嵌套 10x matrix；对三个有明确标签的患者做配对 core-edge 描述，缺少区域标签的样本不纳入区域比较。",
    "- `GSE117891`、`GSE200984`：本轮登记为逐步纳入队列；若缺少 cell annotation，则先不输出正式 cell-type 结论。",
    "",
    "## 关键输出",
    "",
    "```text",
    "tables/lap3_scrna_cell_type_summary.csv",
    "tables/lap3_scrna_sample_cell_type_summary.csv",
    "tables/lap3_scrna_dataset_qc.csv",
    "tables/lap3_scrna_validation_plan.csv",
    "tables/gse237673_lap3_marker_compartment_summary.csv",
    "tables/gse237673_lap3_marker_sample_summary.csv",
    "tables/gse237673_marker_gene_coverage.csv",
    "tables/gse318564_sample_metadata_calibrated.csv",
    "tables/gse318564_lap3_sample_summary.csv",
    "tables/gse318564_lap3_compartment_summary.csv",
    "tables/gse318564_marker_gene_coverage.csv",
    "plots/gse138794_lap3_cell_type_dotplot.png",
    "plots/gse138794_lap3_cell_type_boxplot.png",
    "plots/gse237673_lap3_marker_compartment_dotplot.png",
    "plots/gse237673_marker_scores_vs_lap3.png",
    "plots/gse318564_paired_core_edge_lap3.png",
    "logs/lap3_scrna_localization.log",
    "```",
    "",
    "## 第一版结果摘要",
    "",
    "`GSE138794` 的 LAP3 高表达/高检出 cell type 前几项如下：",
    "",
    "| Dataset | Modality | Cell type | Cells | Samples | LAP3+ % | Mean log1p(CPM) |",
    "|---|---|---|---:|---:|---:|---:|"
  )

  if (nrow(top_gse138) > 0) {
    table_lines <- sprintf(
      "| %s | %s | %s | %s | %s | %.2f | %.3f |",
      top_gse138$dataset,
      top_gse138$modality,
      top_gse138$cell_type,
      top_gse138$n_cells,
      top_gse138$n_samples,
      top_gse138$pct_lap3_positive,
      top_gse138$mean_lap3_log1p_cpm
    )
  } else {
    table_lines <- "| NA | NA | NA | 0 | 0 | 0 | 0 |"
  }

  if (nrow(top_gse211376) > 0) {
    gse211376_lines <- sprintf(
      "| %s | %s | %s | %s | %.2f | %.3f |",
      top_gse211376$dataset,
      top_gse211376$cell_type,
      top_gse211376$n_cells,
      top_gse211376$n_samples,
      top_gse211376$pct_lap3_positive,
      top_gse211376$mean_lap3_count
    )
  } else {
    gse211376_lines <- "| NA | NA | 0 | 0 | 0 | 0 |"
  }

  lines <- c(
    lines,
    table_lines,
    "",
    "`GSE211376` 作为第一批验证数据集，LAP3 检出比例较高的 cell type 如下；该数据集本轮使用 raw count 行，因此不与 `GSE138794` 的 CPM 值直接比较：",
    "",
    "| Dataset | Cell type | Cells | Samples | LAP3+ % | Mean raw count |",
    "|---|---|---:|---:|---:|---:|",
    gse211376_lines,
    "",
    "`GSE237673` 作为 macrophage/foam-cell focused 验证数据集，本轮没有官方 cell annotation，因此只做 marker-defined compartment 验证，不把结果写成正式细胞类型注释：",
    "",
    "| Dataset | Marker compartment | Cells | Samples | LAP3+ % | Mean LAP3 log1p(CPM) | Myeloid score | Foam/lipid score |",
    "|---|---|---:|---:|---:|---:|---:|---:|",
    gse237673_lines,
    "",
    "## 数据集处理状态",
    "",
    "| Dataset | Status | Note |",
    "|---|---|---|",
    sprintf("| %s | %s | %s |", validation_plan$dataset, validation_plan$status, validation_plan$note),
    "",
    "## 解释边界",
    "",
    "- 第一版定位结果以 `GSE138794` 为主，`GSE211376` 为验证性补充；不同数据集的表达单位暂不直接横向比较绝对值。",
    "- `GSE211376` 本轮仅抽取 `LAP3` raw count 行，缺少 library size normalization；其验证重点是 cell-type 检出方向，而不是表达量绝对大小。",
    "- `GSE237673` 缺少独立 cell annotation，本轮采用 myeloid/macrophage 与 foam/lipid marker-defined compartments；这不是作者原始 cell-type annotation。",
    "- `GSE318564` 中只有 WU-1225、WU-1226、GBM034 具有 GEO-confirmed core/edge 配对；B178/B183/B189 没有区域标签，不进入区域比较。",
    "- `GSE318564` 的 B183 三个 GEO archive 缺少 features/barcodes，仅含 matrix，无法可靠解析，已标记为 `incomplete_archive`。",
    "- 三组 core-edge 配对方向不完全一致，且 WU-1226 为 CD45-negative selected；该结果支持 LAP3 跨空间区域可检出，不支持稳定的 core 或 edge 特异性。",
    "- 缺 annotation 的 10x raw matrix 数据集不能直接支持正式 cell-type localization，需先补 metadata / annotation 或保守使用 marker-defined 结果。",
    "",
    "## 下一步",
    "",
    "1. 复核 `GSE138794` top cell types 是否符合论文原始注释和 glioma biology。",
    "2. 若能定位 `GSE237673` 原文/补充材料中的 cell annotation，再替换或校准当前 marker-defined compartments。",
    "3. 若能获得 `GSE318564` 作者 cell annotation 或完整对象，再按作者细胞类型复核 marker-defined 结果。",
    "4. Figure 5 继续写作 malignant-microenvironment spatial context，而不是 malignant-only 或 core/edge-specific。"
  )

  writeLines(lines, file.path(out_dir, "README_LAP3_scRNA_Localization.md"))
}

gse138794_annotation <- rbindlist(list(
  read_annotation(file.path(data_dir, "GSE138794", "GSE138794_scRNA_Seq_cell_types.txt.gz"), "scRNA"),
  read_annotation(file.path(data_dir, "GSE138794", "GSE138794_snRNA_Seq_cell_types.txt.gz"), "snRNA")
), fill = TRUE)

log_msg("GSE138794 annotated cells: ", nrow(gse138794_annotation))
gse138794 <- read_10x_lap3_from_tar(
  file.path(data_dir, "GSE138794", "GSE138794_RAW.tar"),
  gse138794_annotation
)
log_msg("GSE138794 matched cells: ", nrow(gse138794$cells))

gse211376 <- tryCatch(
  read_lap3_from_wide_matrix(
    file.path(data_dir, "GSE211376", "GSE211376_raw_counts_Ruiz2022_all_samples_filtered_cells.tsv.gz"),
    file.path(data_dir, "GSE211376", "GSE211376_metadata_Ruiz2022_all_samples_filtered_cells.csv.gz"),
    "GSE211376"
  ),
  error = function(e) {
    log_msg("GSE211376 validation extraction failed: ", conditionMessage(e))
    data.table()
  }
)
log_msg("GSE211376 matched cells: ", nrow(gse211376))

gse237673 <- tryCatch(
  read_marker_panel_from_10x_tar(
    file.path(data_dir, "GSE237673", "GSE237673_RAW.tar"),
    "GSE237673"
  ),
  error = function(e) {
    log_msg("GSE237673 marker-defined extraction failed: ", conditionMessage(e))
    list(cells = data.table(), qc = data.table(), coverage = data.table())
  }
)
log_msg("GSE237673 marker-defined cells: ", nrow(gse237673$cells))

gse318564_metadata <- tryCatch(
  parse_gse318564_sample_metadata(file.path(data_dir, "GSE318564", "GSE318564_RAW.tar")),
  error = function(e) {
    log_msg("GSE318564 metadata skeleton failed: ", conditionMessage(e))
    data.table()
  }
)
log_msg("GSE318564 metadata skeleton rows: ", nrow(gse318564_metadata))

gse318564 <- tryCatch(
  read_gse318564_marker_panel(
    file.path(data_dir, "GSE318564", "GSE318564_RAW.tar"),
    gse318564_metadata
  ),
  error = function(e) {
    log_msg("GSE318564 marker extraction failed: ", conditionMessage(e))
    list(cells = data.table(), qc = data.table(), coverage = data.table())
  }
)
log_msg("GSE318564 marker-defined cells: ", nrow(gse318564$cells))

all_cells <- rbindlist(list(gse138794$cells, gse211376), fill = TRUE)
if (nrow(all_cells) == 0) {
  stop("No cells were processed.")
}

cell_type_summary <- summarize_by_cell_type(all_cells)
sample_summary <- summarize_by_sample(all_cells)

dataset_qc <- rbindlist(list(
  gse138794$qc,
  data.table(
    dataset = "GSE211376",
    sample_id = if (nrow(gse211376) > 0) "all" else NA_character_,
    modality = if (nrow(gse211376) > 0) "scRNA" else NA_character_,
    matrix_cells = if (nrow(gse211376) > 0) length(unique(gse211376$cell_id)) else 0L,
    annotated_cells = nrow(gse211376),
    lap3_positive_cells = if (nrow(gse211376) > 0) sum(gse211376$lap3_detected) else 0L,
    status = if (nrow(gse211376) > 0) "ok" else "not processed"
  ),
  gse237673$qc,
  gse318564$qc
), fill = TRUE)

validation_plan <- data.table(
  dataset = c("GSE138794", "GSE211376", "GSE117891", "GSE237673", "GSE200984", "GSE318564"),
  status = c(
    "processed_main",
    if (nrow(gse211376) > 0) "processed_validation" else "pending_fix",
    "pending_annotation",
    if (nrow(gse237673$cells) > 0) "processed_marker_defined" else "pending_annotation",
    "pending_metadata_tissue_state",
    if (nrow(gse318564$cells) > 0) "processed_marker_defined_region" else "pending_processing"
  ),
  note = c(
    "Main first-pass localization using scRNA/snRNA cell-type annotations.",
    "Metadata contains cell type; current run extracts LAP3 raw-count row.",
    "Wide matrix and sample metadata available, but cell-type annotation is not yet identified.",
    "No official cell annotation available locally; current run uses myeloid/macrophage and foam/lipid marker-defined compartments.",
    "Raw tar contains matrix files; tissue-state metadata/annotation needs localization before cell-type conclusion.",
    "GEO-confirmed region labels are available for three paired patients; B178/B183 replicates and B189 lack region labels and are excluded from core-edge comparison."
  )
)

if (nrow(gse237673$cells) > 0) {
  gse237673_compartment_summary <- summarize_gse237673_marker_compartments(gse237673$cells)
  gse237673_sample_summary <- summarize_gse237673_sample_compartments(gse237673$cells)
  fwrite(gse237673$cells, file.path(tables_dir, "gse237673_lap3_marker_cell_values.csv"))
  fwrite(gse237673_compartment_summary, file.path(tables_dir, "gse237673_lap3_marker_compartment_summary.csv"))
  fwrite(gse237673_sample_summary, file.path(tables_dir, "gse237673_lap3_marker_sample_summary.csv"))
  fwrite(gse237673$coverage, file.path(tables_dir, "gse237673_marker_gene_coverage.csv"))
  write_gse237673_plot(gse237673$cells, gse237673_compartment_summary)
}
if (nrow(gse318564_metadata) > 0) {
  fwrite(gse318564_metadata, file.path(tables_dir, "gse318564_sample_metadata_calibrated.csv"))
}
if (nrow(gse318564$cells) > 0) {
  gse318564_summary <- summarize_gse318564(gse318564$cells)
  fwrite(gse318564$cells, file.path(tables_dir, "gse318564_lap3_marker_cell_values.csv"))
  fwrite(gse318564_summary$sample, file.path(tables_dir, "gse318564_lap3_sample_summary.csv"))
  fwrite(gse318564_summary$compartment, file.path(tables_dir, "gse318564_lap3_compartment_summary.csv"))
  fwrite(gse318564$coverage, file.path(tables_dir, "gse318564_marker_gene_coverage.csv"))

  paired_plot_dt <- gse318564_summary$sample[include_core_edge == TRUE]
  if (nrow(paired_plot_dt) > 0) {
    p3 <- ggplot(
      paired_plot_dt,
      aes(x = region, y = mean_lap3_log1p_cpm, group = patient_id, color = cd45_selection)
    ) +
      geom_line(alpha = 0.65) +
      geom_point(size = 3) +
      geom_text(aes(label = patient_id), hjust = -0.08, size = 3, show.legend = FALSE) +
      coord_cartesian(clip = "off") +
      labs(
        x = NULL,
        y = "Mean LAP3 log1p(CPM)",
        color = "Selection",
        title = "GSE318564 paired core-edge LAP3 expression"
      ) +
      theme_bw(base_size = 11) +
      theme(plot.margin = margin(5.5, 65, 5.5, 5.5))
    ggsave(
      file.path(plots_dir, "gse318564_paired_core_edge_lap3.png"),
      p3, width = 7.2, height = 4.8, dpi = 300
    )
  }
}

fwrite(all_cells, file.path(tables_dir, "lap3_scrna_cell_level_values.csv"))
fwrite(cell_type_summary, file.path(tables_dir, "lap3_scrna_cell_type_summary.csv"))
fwrite(sample_summary, file.path(tables_dir, "lap3_scrna_sample_cell_type_summary.csv"))
fwrite(dataset_qc, file.path(tables_dir, "lap3_scrna_dataset_qc.csv"))
fwrite(validation_plan, file.path(tables_dir, "lap3_scrna_validation_plan.csv"))

write_plot(cell_type_summary, all_cells)
write_readme(cell_type_summary, dataset_qc, validation_plan)

log_msg("Wrote outputs to: ", out_dir)
log_msg("LAP3 scRNA localization finished: ", Sys.time())
