#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="/home/lzb/glioma"
SCRIPT="${PROJECT_ROOT}/Data_Spatial_Public/scripts/gbmspace_lap3_state_topology.R"
OUT_DIR="${PROJECT_ROOT}/Data_Spatial_Public/GBM_Space/results/LAP3_State_Spatial_Topology"
LOG_DIR="${OUT_DIR}/logs"
mkdir -p "${LOG_DIR}"

RUN_ID="$(date +%Y%m%d_%H%M%S)_$$"
LOG_FILE="${LOG_DIR}/gbmspace_lap3_state_topology.${RUN_ID}.log"
STATUS_FILE="${LOG_DIR}/gbmspace_lap3_state_topology.status"
HEARTBEAT_FILE="${LOG_DIR}/gbmspace_lap3_state_topology.heartbeat.log"

echo "$(date '+%F %T %Z') RUNNING pid=$$ run_id=${RUN_ID} log=${LOG_FILE}" > "${STATUS_FILE}"

(
  while true; do
    echo "$(date '+%F %T %Z') heartbeat pid=$$ run_id=${RUN_ID}" >> "${HEARTBEAT_FILE}"
    sleep 60
  done
) &
HB_PID=$!

cleanup() {
  kill "${HB_PID}" >/dev/null 2>&1 || true
}
trap cleanup EXIT

cd "${PROJECT_ROOT}"
echo "$(date '+%F %T %Z') START GBM-Space LAP3-state topology" | tee -a "${LOG_FILE}"

if Rscript --vanilla "${SCRIPT}" >> "${LOG_FILE}" 2>&1; then
  echo "$(date '+%F %T %Z') DONE exit_code=0 run_id=${RUN_ID} log=${LOG_FILE} result=${OUT_DIR}" > "${STATUS_FILE}"
  echo "$(date '+%F %T %Z') DONE" | tee -a "${LOG_FILE}"
else
  ec=$?
  echo "$(date '+%F %T %Z') FAILED exit_code=${ec} run_id=${RUN_ID} log=${LOG_FILE}" > "${STATUS_FILE}"
  echo "$(date '+%F %T %Z') FAILED exit_code=${ec}" | tee -a "${LOG_FILE}"
  exit "${ec}"
fi
