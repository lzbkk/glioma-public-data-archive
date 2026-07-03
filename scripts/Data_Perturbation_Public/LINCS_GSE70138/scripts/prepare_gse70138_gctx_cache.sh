#!/usr/bin/env bash
set -euo pipefail

cd /home/lzb/glioma

root="Data_Perturbation_Public/LINCS_GSE70138"
data_dir="${root}/data"
cache_dir="${root}/cache"
log_dir="${root}/logs"
source_gz="${data_dir}/GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx.gz"
target_gctx="${cache_dir}/GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx"
structure_tsv="${cache_dir}/gse70138_gctx_h5_structure.tsv"
status="${log_dir}/gse70138_cache.status"
run_id="$(date '+%Y%m%d_%H%M%S')_$$"
run_log="${log_dir}/gse70138_cache.${run_id}.log"
latest_log="${log_dir}/gse70138_cache.log"

mkdir -p "${cache_dir}" "${log_dir}"
ln -sfn "$(basename "${run_log}")" "${latest_log}"
exec >>"${run_log}" 2>&1

on_exit() {
  exit_code=$?
  if [[ $exit_code -ne 0 ]]; then
    date "+%F %T FAILED pid=$$ run_id=$run_id exit_code=$exit_code log=$run_log" >"${status}"
    date "+%F %T | FAILED GSE70138 GCTX cache exit_code=$exit_code"
  fi
}
trap on_exit EXIT

date "+%F %T RUNNING pid=$$ run_id=$run_id stage=decompress log=$run_log" >"${status}"
date "+%F %T | START GSE70138 GCTX cache"

if [[ ! -s "${target_gctx}" ]]; then
  gzip -dc "${source_gz}" >"${target_gctx}.partial"
  mv "${target_gctx}.partial" "${target_gctx}"
else
  date "+%F %T | REUSE existing uncompressed GCTX"
fi

date "+%F %T RUNNING pid=$$ run_id=$run_id stage=inspect log=$run_log" >"${status}"
GCTX_PATH="${target_gctx}" STRUCTURE_TSV="${structure_tsv}" Rscript --vanilla - <<'RS'
suppressPackageStartupMessages(library(rhdf5))

gctx_path <- Sys.getenv("GCTX_PATH")
structure_tsv <- Sys.getenv("STRUCTURE_TSV")

h5 <- h5ls(gctx_path, recursive = TRUE)
write.table(h5, structure_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

matrix_row <- h5[h5$name == "matrix" & h5$otype == "H5I_DATASET", , drop = FALSE]
if (nrow(matrix_row) != 1L) {
  stop("Expected exactly one HDF5 dataset named matrix; found ", nrow(matrix_row))
}

cat("matrix dimensions:", matrix_row$dim[[1]], "\n")
cat("HDF5 objects:", nrow(h5), "\n")
RS

date "+%F %T DONE pid=$$ run_id=$run_id cache=$target_gctx structure=$structure_tsv log=$run_log" >"${status}"
date "+%F %T | DONE GSE70138 GCTX cache"
