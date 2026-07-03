#!/usr/bin/env bash
set -euo pipefail

project_root="/home/lzb/glioma"
tarball="${project_root}/Data_Spatial_Public/GBM_Space/data/spatial_data_visium.tar.gz"
manifest_verbose="${project_root}/Data_Spatial_Public/GBM_Space/results/Entry_Audit/tables/gbmspace_visium_tar_verbose.txt"
extract_dir="${project_root}/Data_Spatial_Public/GBM_Space/cache/selected_h5ad_field_audit"
out_dir="${project_root}/Data_Spatial_Public/GBM_Space/results/H5AD_Field_Audit"
tables_dir="${out_dir}/tables"
logs_dir="${project_root}/Data_Spatial_Public/GBM_Space/logs"
run_id="$(date +%Y%m%d_%H%M%S)_$$"
status_file="${logs_dir}/gbmspace_h5ad_field_audit.status"
log_file="${logs_dir}/gbmspace_h5ad_field_audit.${run_id}.log"

mkdir -p "${extract_dir}" "${tables_dir}" "${logs_dir}"

log() {
  printf '%s | %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"
}

write_status() {
  printf '%s %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" > "${status_file}"
}

{
  log "START GBM-Space H5AD field audit run_id=${run_id}"
  write_status "RUNNING pid=$$ run_id=${run_id} step=select_files log=${log_file}"

  if [[ ! -s "${tarball}" ]]; then
    write_status "FAILED pid=$$ run_id=${run_id} reason=tarball_missing log=${log_file}"
    log "FAILED tarball missing: ${tarball}"
    exit 1
  fi
  if [[ ! -s "${manifest_verbose}" ]]; then
    write_status "FAILED pid=$$ run_id=${run_id} reason=manifest_missing log=${log_file}"
    log "FAILED manifest missing: ${manifest_verbose}"
    exit 1
  fi

  selected_tsv="${tables_dir}/gbmspace_h5ad_selected_files.tsv"
  small_path="$(awk '$1 ~ /^-/ && $6 ~ /\.h5ad$/ {print $3"\t"$6}' "${manifest_verbose}" | sort -n | awk 'NR==1{print $2}')"
  medium_path="anndata/AT12-BRA-4-FO-2_1.h5ad"
  large_path="anndata/AT15-BRA-4-FO-F1-S31.h5ad"
  readme_path="anndata/README.md"

  {
    printf 'role\ttar_path\n'
    printf 'readme\t%s\n' "${readme_path}"
    printf 'smallest_h5ad\t%s\n' "${small_path}"
    printf 'medium_representative_h5ad\t%s\n' "${medium_path}"
    printf 'large_at15_representative_h5ad\t%s\n' "${large_path}"
  } > "${selected_tsv}"

  log "Selected files:"
  cat "${selected_tsv}"

  write_status "RUNNING pid=$$ run_id=${run_id} step=extract_selected log=${log_file}"
  log "Extracting selected members from tar.gz"
  tar -xzf "${tarball}" -C "${extract_dir}" \
    "${readme_path}" \
    "${small_path}" \
    "${medium_path}" \
    "${large_path}"

  write_status "RUNNING pid=$$ run_id=${run_id} step=inspect_h5ad log=${log_file}"
  log "Inspecting H5AD files with rhdf5"
  Rscript --vanilla "${project_root}/Data_Spatial_Public/scripts/audit_gbmspace_h5ad_fields.R" \
    "${extract_dir}" \
    "${out_dir}"

  write_status "DONE pid=$$ run_id=${run_id} selected=3 log=${log_file} result=${out_dir}"
  log "DONE result=${out_dir}"
} >> "${log_file}" 2>&1
