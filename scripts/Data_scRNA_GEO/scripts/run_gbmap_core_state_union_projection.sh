#!/usr/bin/env bash
set -uo pipefail

project_dir="/home/lzb/glioma"
module_dir="${project_dir}/Data_scRNA_GEO/GBmap_Core"
result_dir="${module_dir}/results/LAP3_State_Union_Projection"
logs_dir="${result_dir}/logs"
status_file="${logs_dir}/gbmap_core_state_union_projection.status"
wrapper_log="${logs_dir}/gbmap_core_state_union_projection_wrapper.log"
heartbeat_file="${logs_dir}/gbmap_core_state_union_projection.heartbeat.log"
input_file="${module_dir}/data/core_gbmap_cellxgene.h5ad"
gene_set_file="${project_dir}/Data_Bulk_TCGA/Data_Merged/results/LAP3_State_Module/tables/lap3_state_frozen_gene_sets.csv"
cache_file="${module_dir}/cache/core_gbmap_lap3_state_union_lightweight.rds"
run_id="$(date '+%Y%m%d_%H%M%S')_$$"

mkdir -p "${logs_dir}" "$(dirname "${cache_file}")"
exec > >(tee -a "${wrapper_log}") 2>&1

write_status() {
  printf '%s %s pid=%s run_id=%s %s\n' \
    "$(date '+%F %T %Z')" "$1" "$$" "${run_id}" "${2:-}" > "${status_file}"
}

heartbeat() {
  while kill -0 "$$" 2>/dev/null; do
    printf '%s run_id=%s pid=%s cache_bytes=%s\n' \
      "$(date '+%F %T %Z')" "${run_id}" "$$" \
      "$(stat -c %s "${cache_file}" 2>/dev/null || printf '0')" \
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
    write_status "DONE" "exit_code=0 cache=${cache_file} result=${result_dir} log=${wrapper_log}"
  else
    write_status "FAILED" "exit_code=${exit_code} log=${wrapper_log}"
  fi
  exit "${exit_code}"
}
trap finish EXIT INT TERM

write_status "RUNNING" "stage=build_state_cache cache=${cache_file} log=${wrapper_log}"
heartbeat &
heartbeat_pid=$!

export OMP_NUM_THREADS=16
export OPENBLAS_NUM_THREADS=16
export MKL_NUM_THREADS=16

Rscript --vanilla \
  "${project_dir}/Data_scRNA_GEO/scripts/gbmap_core_build_lap3_state_cache.R" \
  "${input_file}" "${gene_set_file}" "${cache_file}"

write_status "RUNNING" "stage=project_state_union cache=${cache_file} log=${wrapper_log}"

Rscript --vanilla \
  "${project_dir}/Data_scRNA_GEO/scripts/lap3_gbmap_core_state_union_projection.R"
