#!/usr/bin/env bash
set -euo pipefail

cd /home/lzb/glioma

CREDENTIALS="/home/lzb/.config/glioma/synapse.env"
LOG_DIR="Data_Longitudinal_Public/GLASS/logs"
STATUS="${LOG_DIR}/glass_fetch.status"
LATEST_LOG="${LOG_DIR}/glass_fetch.log"
RUN_ID="$(date '+%Y%m%d_%H%M%S')_$$"
RUN_LOG="${LOG_DIR}/glass_fetch.${RUN_ID}.log"

mkdir -p "${LOG_DIR}"
ln -sfn "$(basename "${RUN_LOG}")" "${LATEST_LOG}"
exec >>"${RUN_LOG}" 2>&1

if [[ ! -r "${CREDENTIALS}" ]]; then
  date "+%F %T FAILED pid=$$ run_id=$RUN_ID reason=credentials_unreadable" >"${STATUS}"
  exit 1
fi

source "${CREDENTIALS}"
if [[ -z "${SYNAPSE_AUTH_TOKEN:-}" ]]; then
  date "+%F %T FAILED pid=$$ run_id=$RUN_ID reason=token_unset" >"${STATUS}"
  exit 1
fi

date "+%F %T RUNNING pid=$$ run_id=$RUN_ID stage=download log=$RUN_LOG" >"${STATUS}"
date "+%F %T | START GLASS current-release fetch run_id=$RUN_ID"

on_exit() {
  exit_code=$?
  if [[ $exit_code -ne 0 ]]; then
    date "+%F %T FAILED pid=$$ run_id=$RUN_ID exit_code=$exit_code log=$RUN_LOG" >"${STATUS}"
    date "+%F %T | FAILED GLASS current-release fetch exit_code=$exit_code"
  fi
}
trap on_exit EXIT

Data_Longitudinal_Public/GLASS/.venv/bin/python \
  Data_Longitudinal_Public/GLASS/scripts/glass_fetch_current_release.py

test -s Data_Longitudinal_Public/GLASS/data/tables/clinical_cases.csv
test -s Data_Longitudinal_Public/GLASS/data/tables/clinical_surgeries.csv
test -s Data_Longitudinal_Public/GLASS/data/tables/analysis_tumor_pairs.csv
test -s Data_Longitudinal_Public/GLASS/data/tables/analysis_rna_silver_set.csv
test -s Data_Longitudinal_Public/GLASS/data/expression/gene_tpm_matrix_all_samples.tsv

date "+%F %T DONE pid=$$ run_id=$RUN_ID files=5 failed=0 log=$RUN_LOG" >"${STATUS}"
date "+%F %T | DONE GLASS current-release fetch"
