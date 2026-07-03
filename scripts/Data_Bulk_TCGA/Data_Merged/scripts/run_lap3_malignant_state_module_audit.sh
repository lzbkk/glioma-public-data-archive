#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="/home/lzb/glioma"
OUT_DIR="${PROJECT_DIR}/Data_Bulk_TCGA/Data_Merged/results/LAP3_Malignant_State_Module_Audit"
LOG_DIR="${OUT_DIR}/logs"
mkdir -p "${LOG_DIR}"

RUN_ID="$(date +%Y%m%d_%H%M%S)_$$"
STATUS_FILE="${LOG_DIR}/lap3_malignant_state_module_audit.status"
WRAPPER_LOG="${LOG_DIR}/lap3_malignant_state_module_audit_wrapper.${RUN_ID}.log"
HEARTBEAT_LOG="${LOG_DIR}/lap3_malignant_state_module_audit.heartbeat.${RUN_ID}.log"

echo "$(date '+%F %T') RUNNING pid=$$ run_id=${RUN_ID} log=${WRAPPER_LOG}" > "${STATUS_FILE}"

heartbeat() {
  while true; do
    echo "$(date '+%F %T') heartbeat pid=$$ run_id=${RUN_ID}" >> "${HEARTBEAT_LOG}"
    sleep 60
  done
}
heartbeat &
HB_PID=$!
trap 'kill "${HB_PID}" >/dev/null 2>&1 || true' EXIT

cd "${PROJECT_DIR}"
{
  echo "$(date '+%F %T') START LAP3 malignant-state module audit run_id=${RUN_ID}"
  Rscript --vanilla Data_Bulk_TCGA/Data_Merged/scripts/lap3_malignant_state_module_audit.R
  echo "$(date '+%F %T') DONE LAP3 malignant-state module audit run_id=${RUN_ID}"
} > "${WRAPPER_LOG}" 2>&1

echo "$(date '+%F %T') DONE pid=$$ run_id=${RUN_ID} log=${WRAPPER_LOG}" > "${STATUS_FILE}"
