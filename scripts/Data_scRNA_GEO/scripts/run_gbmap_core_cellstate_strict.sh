#!/usr/bin/env bash
set -uo pipefail

project_dir="/home/lzb/glioma"
status_file="${project_dir}/Data_scRNA_GEO/GBmap_Core/results/LAP3_CellState_Strict/logs/gbmap_core_cellstate_strict.status"
wrapper_log="${project_dir}/Data_scRNA_GEO/GBmap_Core/results/LAP3_CellState_Strict/logs/gbmap_core_cellstate_strict_wrapper.log"
heartbeat_file="${project_dir}/Data_scRNA_GEO/GBmap_Core/results/LAP3_CellState_Strict/logs/gbmap_core_cellstate_strict.heartbeat.log"
run_id="$(date '+%Y%m%d_%H%M%S')_$$"

mkdir -p "$(dirname "${status_file}")"
exec > >(tee -a "${wrapper_log}") 2>&1

write_status() {
  printf '%s %s pid=%s run_id=%s %s\n' \
    "$(date '+%F %T')" "$1" "$$" "${run_id}" "${2:-}" > "${status_file}"
}

heartbeat() {
  while kill -0 "$$" 2>/dev/null; do
    printf '%s run_id=%s pid=%s\n' "$(date '+%F %T')" "${run_id}" "$$" > "${heartbeat_file}"
    sleep 60
  done
}

finish() {
  local exit_code=$?
  if [[ -n "${heartbeat_pid:-}" ]]; then
    kill "${heartbeat_pid}" 2>/dev/null || true
    wait "${heartbeat_pid}" 2>/dev/null || true
  fi
  if ((exit_code == 0)); then
    write_status "DONE" "log=${wrapper_log}"
  else
    write_status "FAILED" "exit_code=${exit_code} log=${wrapper_log}"
  fi
  exit "${exit_code}"
}
trap finish EXIT INT TERM

write_status "RUNNING" "stage=gbmap_core_cellstate_strict"
heartbeat &
heartbeat_pid=$!

cd "${project_dir}" || exit 1
export OMP_NUM_THREADS=16
export OPENBLAS_NUM_THREADS=16
export MKL_NUM_THREADS=16
Rscript --vanilla "${project_dir}/Data_scRNA_GEO/scripts/lap3_gbmap_core_cellstate_strict.R"
