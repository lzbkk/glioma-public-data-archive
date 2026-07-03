#!/usr/bin/env bash
set -uo pipefail

cd /home/lzb/glioma

result_dir="Data_scRNA_GEO/results/GSE131928_LAP3_CellState"
log_dir="${result_dir}/logs"
log_file="${log_dir}/gse131928_cellstate_sensitivity_wrapper.log"
status_file="${log_dir}/gse131928_cellstate_sensitivity.status"
heartbeat_file="${log_dir}/gse131928_cellstate_sensitivity.heartbeat.log"

mkdir -p "${log_dir}"
printf '%s RUNNING pid=%s\n' "$(date)" "$$" > "${status_file}"

heartbeat() {
  while true; do
    printf '%s RUNNING pid=%s\n' "$(date)" "$$" >> "${heartbeat_file}"
    sleep 300
  done
}

heartbeat &
heartbeat_pid=$!
trap 'kill "${heartbeat_pid}" 2>/dev/null || true' EXIT

export OMP_NUM_THREADS=16
export OPENBLAS_NUM_THREADS=16
export MKL_NUM_THREADS=16

if Rscript --vanilla Data_scRNA_GEO/scripts/lap3_gse131928_cellstate_calibration.R \
  > "${log_file}" 2>&1; then
  printf '%s DONE pid=%s\n' "$(date)" "$$" > "${status_file}"
else
  exit_code=$?
  printf '%s FAILED exit=%s pid=%s\n' "$(date)" "${exit_code}" "$$" > "${status_file}"
  exit "${exit_code}"
fi
