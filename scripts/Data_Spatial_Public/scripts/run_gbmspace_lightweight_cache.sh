#!/usr/bin/env bash
set -euo pipefail

project_root="/home/lzb/glioma"
tarball="${project_root}/Data_Spatial_Public/GBM_Space/data/spatial_data_visium.tar.gz"
extract_dir="${project_root}/Data_Spatial_Public/GBM_Space/cache/all_h5ad"
out_dir="${project_root}/Data_Spatial_Public/GBM_Space/results/Lightweight_Cache"
gene_set_csv="${project_root}/Data_scRNA_GEO/results/LAP3_CellState_Phase0/tables/frozen_lap3_pathway_gene_sets.csv"
logs_dir="${project_root}/Data_Spatial_Public/GBM_Space/logs"
run_id="$(date +%Y%m%d_%H%M%S)_$$"
status_file="${logs_dir}/gbmspace_lightweight_cache.status"
log_file="${logs_dir}/gbmspace_lightweight_cache.${run_id}.log"

mkdir -p "${extract_dir}" "${out_dir}/tables" "${logs_dir}"

log() {
  printf '%s | %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"
}

write_status() {
  printf '%s %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" > "${status_file}"
}

{
  log "START GBM-Space lightweight cache run_id=${run_id}"
  write_status "RUNNING pid=$$ run_id=${run_id} step=precheck log=${log_file}"

  if [[ ! -s "${tarball}" ]]; then
    write_status "FAILED pid=$$ run_id=${run_id} reason=tarball_missing log=${log_file}"
    log "FAILED tarball missing: ${tarball}"
    exit 1
  fi
  if [[ ! -s "${gene_set_csv}" ]]; then
    write_status "FAILED pid=$$ run_id=${run_id} reason=gene_set_missing log=${log_file}"
    log "FAILED gene set missing: ${gene_set_csv}"
    exit 1
  fi

  current_h5ad_count="$(find "${extract_dir}/anndata" -maxdepth 1 -name '*.h5ad' 2>/dev/null | wc -l || true)"
  if [[ "${current_h5ad_count}" -lt 97 ]]; then
    write_status "RUNNING pid=$$ run_id=${run_id} step=extract_all_h5ad current=${current_h5ad_count}/97 log=${log_file}"
    log "Extracting all H5AD files to ${extract_dir}; current=${current_h5ad_count}/97"
    tar -xzf "${tarball}" -C "${extract_dir}" anndata
  else
    log "All H5AD files already extracted: ${current_h5ad_count}/97"
  fi

  final_h5ad_count="$(find "${extract_dir}/anndata" -maxdepth 1 -name '*.h5ad' | wc -l)"
  log "Extracted H5AD count=${final_h5ad_count}"
  if [[ "${final_h5ad_count}" -ne 97 ]]; then
    write_status "FAILED pid=$$ run_id=${run_id} reason=h5ad_count_mismatch count=${final_h5ad_count} log=${log_file}"
    log "FAILED expected 97 H5AD files"
    exit 1
  fi

  write_status "RUNNING pid=$$ run_id=${run_id} step=build_cache files=97 log=${log_file}"
  log "Running R lightweight cache builder"
  Rscript --vanilla "${project_root}/Data_Spatial_Public/scripts/build_gbmspace_lightweight_cache.R" \
    "${extract_dir}" \
    "${gene_set_csv}" \
    "${out_dir}"

  write_status "DONE pid=$$ run_id=${run_id} files=97 log=${log_file} result=${out_dir}"
  log "DONE result=${out_dir}"
} >> "${log_file}" 2>&1
