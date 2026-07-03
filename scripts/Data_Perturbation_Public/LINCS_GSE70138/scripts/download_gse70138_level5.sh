#!/usr/bin/env bash
set -euo pipefail

cd /home/lzb/glioma

BASE_URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/"
ROOT="Data_Perturbation_Public/LINCS_GSE70138"
DATA_DIR="${ROOT}/data"
LOG_DIR="${ROOT}/logs"
STATUS="${LOG_DIR}/gse70138_download.status"
LATEST_LOG="${LOG_DIR}/gse70138_download.log"
RUN_ID="$(date '+%Y%m%d_%H%M%S')_$$"
RUN_LOG="${LOG_DIR}/gse70138_download.${RUN_ID}.log"

mkdir -p "${DATA_DIR}" "${LOG_DIR}"
ln -sfn "$(basename "${RUN_LOG}")" "${LATEST_LOG}"
exec >>"${RUN_LOG}" 2>&1

files=(
  "GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx.gz"
  "GSE70138_Broad_LINCS_gene_info_2017-03-06.txt.gz"
  "GSE70138_Broad_LINCS_pert_info_2017-03-06.txt.gz"
  "GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz"
  "GSE70138_Broad_LINCS_sig_metrics_2017-03-06.txt.gz"
  "GSE70138_SHA512SUMS.txt.gz"
)

date "+%F %T RUNNING pid=$$ run_id=$RUN_ID stage=download files=${#files[@]} log=$RUN_LOG" >"${STATUS}"
date "+%F %T | START GSE70138 Level 5 download run_id=$RUN_ID"

on_exit() {
  exit_code=$?
  if [[ $exit_code -ne 0 ]]; then
    date "+%F %T FAILED pid=$$ run_id=$RUN_ID exit_code=$exit_code log=$RUN_LOG" >"${STATUS}"
    date "+%F %T | FAILED GSE70138 download exit_code=$exit_code"
  fi
}
trap on_exit EXIT

download_with_retry() {
  local url="$1"
  local destination="$2"
  local max_attempts=12
  local attempt
  local wait_seconds

  for ((attempt = 1; attempt <= max_attempts; attempt++)); do
    date "+%F %T | DOWNLOAD attempt=${attempt}/${max_attempts} destination=${destination}"
    if curl \
      --fail \
      --location \
      --connect-timeout 30 \
      --continue-at - \
      --output "${destination}" \
      "${url}"; then
      return 0
    fi

    if ((attempt < max_attempts)); then
      wait_seconds=$((attempt * 10))
      ((wait_seconds > 120)) && wait_seconds=120
      date "+%F %T | RETRY in ${wait_seconds}s destination=${destination}"
      sleep "${wait_seconds}"
    fi
  done

  date "+%F %T | DOWNLOAD exhausted retries destination=${destination}"
  return 1
}

for i in "${!files[@]}"; do
  file="${files[$i]}"
  index=$((i + 1))
  date "+%F %T RUNNING pid=$$ run_id=$RUN_ID stage=download file=${index}/${#files[@]} name=$file log=$RUN_LOG" >"${STATUS}"
  download_with_retry \
    "${BASE_URL}${file}" \
    "${DATA_DIR}/${file}"
done

gzip -t "${DATA_DIR}/GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx.gz"
gzip -dc "${DATA_DIR}/GSE70138_SHA512SUMS.txt.gz" >"${DATA_DIR}/GSE70138_SHA512SUMS.txt"

(
  cd "${DATA_DIR}"
  sha512sum --check --ignore-missing GSE70138_SHA512SUMS.txt
)

date "+%F %T DONE pid=$$ run_id=$RUN_ID files=${#files[@]} failed=0 log=$RUN_LOG" >"${STATUS}"
date "+%F %T | DONE GSE70138 Level 5 download"
