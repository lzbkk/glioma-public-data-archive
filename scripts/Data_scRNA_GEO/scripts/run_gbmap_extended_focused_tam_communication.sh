#!/usr/bin/env bash
set -euo pipefail

ROOT="/home/lzb/glioma"
RUN_ID="$(date +%Y%m%d_%H%M%S)_$$"
RESULT_DIR="${ROOT}/Data_scRNA_GEO/GBmap_Extended/results/Focused_Malignant_TAM_Communication"
CACHE_DIR="${ROOT}/Data_scRNA_GEO/GBmap_Extended/cache"
LOG_DIR="${RESULT_DIR}/logs"
STATUS_FILE="${LOG_DIR}/gbmap_extended_focused_tam.status"
WRAPPER_LOG="${LOG_DIR}/gbmap_extended_focused_tam.${RUN_ID}.log"
HEARTBEAT_FILE="${LOG_DIR}/gbmap_extended_focused_tam.heartbeat.log"

INPUT_H5AD="${ROOT}/Data_scRNA_GEO/GBmap_Extended/data/extended_gbmap_cellxgene.h5ad"
GENE_SET_FILE="${ROOT}/Data_Bulk_TCGA/Data_Merged/results/LAP3_State_Module/tables/lap3_state_frozen_gene_sets.csv"
SUBMODULE_FILE="${ROOT}/Data_Bulk_TCGA/Data_Merged/results/LAP3_State_Submodules/tables/lap3_state_submodule_gene_assignment.csv"
CACHE_FILE="${CACHE_DIR}/extended_gbmap_focused_tam_lightweight.rds"

mkdir -p "${CACHE_DIR}" "${LOG_DIR}"

status() {
  printf '%s %s\n' "$(date '+%Y-%m-%d %H:%M:%S %Z')" "$*" | tee -a "${STATUS_FILE}"
}

heartbeat() {
  while true; do
    printf '%s RUNNING pid=%s run_id=%s log=%s\n' "$(date '+%Y-%m-%d %H:%M:%S %Z')" "$$" "${RUN_ID}" "${WRAPPER_LOG}" > "${HEARTBEAT_FILE}"
    sleep 60
  done
}

heartbeat &
HEARTBEAT_PID=$!
trap 'kill "${HEARTBEAT_PID}" >/dev/null 2>&1 || true' EXIT

{
  status "RUNNING START Extended GBmap focused malignant-TAM communication pid=$$ run_id=${RUN_ID} log=${WRAPPER_LOG}"
  status "INPUT ${INPUT_H5AD}"

  export OMP_NUM_THREADS="${OMP_NUM_THREADS:-16}"
  export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-16}"
  export MKL_NUM_THREADS="${MKL_NUM_THREADS:-16}"
  export VECLIB_MAXIMUM_THREADS="${VECLIB_MAXIMUM_THREADS:-16}"
  export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-16}"

  if [[ ! -s "${CACHE_FILE}" || "${REBUILD_CACHE:-0}" == "1" ]]; then
    status "RUNNING CACHE_BUILD pid=$$ run_id=${RUN_ID} log=${WRAPPER_LOG} cache=${CACHE_FILE}"
    Rscript --vanilla "${ROOT}/Data_scRNA_GEO/scripts/gbmap_extended_build_focused_tam_cache.R" \
      "${INPUT_H5AD}" \
      "${GENE_SET_FILE}" \
      "${SUBMODULE_FILE}" \
      "${CACHE_FILE}"
  else
    status "RUNNING CACHE_REUSE pid=$$ run_id=${RUN_ID} log=${WRAPPER_LOG} cache=${CACHE_FILE}"
  fi

  status "RUNNING ANALYSIS pid=$$ run_id=${RUN_ID} log=${WRAPPER_LOG} result=${RESULT_DIR}"
  Rscript --vanilla "${ROOT}/Data_scRNA_GEO/scripts/lap3_gbmap_extended_focused_tam_communication.R" \
    "${CACHE_FILE}" \
    "${RESULT_DIR}"

  status "DONE run_id=${RUN_ID} result=${RESULT_DIR} cache=${CACHE_FILE}"
} >> "${WRAPPER_LOG}" 2>&1
