#!/usr/bin/env bash
set -uo pipefail

project_dir="/home/lzb/glioma"
module_dir="${project_dir}/Data_scRNA_GEO/GBmap_Core"
logs_dir="${module_dir}/logs"
status_file="${logs_dir}/gbmap_core_cache.status"
wrapper_log="${logs_dir}/gbmap_core_cache_wrapper.log"
heartbeat_file="${logs_dir}/gbmap_core_cache.heartbeat.log"
input_file="${module_dir}/data/core_gbmap_cellxgene.h5ad"
gene_set_file="${project_dir}/Data_scRNA_GEO/results/LAP3_CellState_Phase0/source_data/frozen_cellstate_gene_sets.rds"
output_file="${module_dir}/cache/core_gbmap_lap3_cellstate_lightweight.rds"
run_id="$(date '+%Y%m%d_%H%M%S')_$$"

mkdir -p "${logs_dir}" "$(dirname "${output_file}")"
exec > >(tee -a "${wrapper_log}") 2>&1

write_status() {
  printf '%s %s pid=%s run_id=%s %s\n' \
    "$(date '+%F %T')" "$1" "$$" "${run_id}" "${2:-}" > "${status_file}"
}

heartbeat() {
  while kill -0 "$$" 2>/dev/null; do
    printf '%s run_id=%s pid=%s cache_bytes=%s\n' \
      "$(date '+%F %T')" "${run_id}" "$$" \
      "$(stat -c %s "${output_file}" 2>/dev/null || printf '0')" \
      > "${heartbeat_file}"
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
    write_status "DONE" "cache=${output_file} bytes=$(stat -c %s "${output_file}")"
  else
    write_status "FAILED" "exit_code=${exit_code}"
  fi
  exit "${exit_code}"
}
trap finish EXIT INT TERM

write_status "RUNNING" "stage=build_lightweight_cache"
heartbeat &
heartbeat_pid=$!

export OMP_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=8
export MKL_NUM_THREADS=8

Rscript --vanilla \
  "${project_dir}/Data_scRNA_GEO/scripts/gbmap_core_build_lightweight_cache.R" \
  "${input_file}" "${gene_set_file}" "${output_file}"
