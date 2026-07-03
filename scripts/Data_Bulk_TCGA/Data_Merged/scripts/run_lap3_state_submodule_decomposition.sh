#!/usr/bin/env bash
set -uo pipefail

project_dir="/home/lzb/glioma"
result_dir="${project_dir}/Data_Bulk_TCGA/Data_Merged/results/LAP3_State_Submodules"
logs_dir="${result_dir}/logs"
status_file="${logs_dir}/lap3_state_submodule_decomposition.status"
wrapper_log="${logs_dir}/lap3_state_submodule_decomposition_wrapper.log"
heartbeat_file="${logs_dir}/lap3_state_submodule_decomposition.heartbeat.log"
run_id="$(date '+%Y%m%d_%H%M%S')_$$"
terminated_signal=""

mkdir -p "${logs_dir}"
exec > >(tee -a "${wrapper_log}") 2>&1

write_status() {
  printf '%s %s pid=%s run_id=%s %s\n' \
    "$(date '+%F %T %Z')" "$1" "$$" "${run_id}" "${2:-}" > "${status_file}"
}

heartbeat() {
  while kill -0 "$$" 2>/dev/null; do
    printf '%s run_id=%s pid=%s result_bytes=%s\n' \
      "$(date '+%F %T %Z')" "${run_id}" "$$" \
      "$(du -sb "${result_dir}" 2>/dev/null | awk '{print $1}' || printf '0')" \
      > "${heartbeat_file}"
    sleep 60
  done
}

finish() {
  local exit_code=$?
  if [[ -n "${terminated_signal}" && "${exit_code}" -eq 0 ]]; then
    exit_code=143
  fi
  if [[ -n "${heartbeat_pid:-}" ]]; then
    kill "${heartbeat_pid}" 2>/dev/null || true
    wait "${heartbeat_pid}" 2>/dev/null || true
  fi
  if ((exit_code == 0)); then
    write_status "DONE" "exit_code=0 result=${result_dir} log=${wrapper_log}"
  else
    write_status "FAILED" "exit_code=${exit_code} log=${wrapper_log}"
  fi
  exit "${exit_code}"
}
trap finish EXIT
trap 'terminated_signal=INT; exit 130' INT
trap 'terminated_signal=TERM; exit 143' TERM

write_status "RUNNING" "stage=submodule_decomposition log=${wrapper_log}"
heartbeat &
heartbeat_pid=$!

export OMP_NUM_THREADS=16
export OPENBLAS_NUM_THREADS=16
export MKL_NUM_THREADS=16

Rscript --vanilla "${project_dir}/Data_Bulk_TCGA/Data_Merged/scripts/lap3_state_submodule_decomposition.R"
