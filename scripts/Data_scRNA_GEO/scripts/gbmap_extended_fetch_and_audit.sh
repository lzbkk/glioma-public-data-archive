#!/usr/bin/env bash
set -uo pipefail

project_dir="/home/lzb/glioma"
module_dir="${project_dir}/Data_scRNA_GEO/GBmap_Extended"
data_dir="${module_dir}/data"
results_dir="${module_dir}/results/Entry_Audit"
logs_dir="${module_dir}/logs"
status_file="${logs_dir}/gbmap_extended_fetch.status"
wrapper_log="${logs_dir}/gbmap_extended_fetch_wrapper.log"
heartbeat_file="${logs_dir}/gbmap_extended_fetch.heartbeat.log"
target_file="${data_dir}/extended_gbmap_cellxgene.h5ad"
source_url="https://datasets.cellxgene.cziscience.com/230c8701-7291-4dd5-bf36-688490c681ee.h5ad"
expected_size="11224714727"
run_id="$(date '+%Y%m%d_%H%M%S')_$$"

mkdir -p "${data_dir}" "${results_dir}" "${logs_dir}"
exec > >(tee -a "${wrapper_log}") 2>&1

write_status() {
  printf '%s %s pid=%s run_id=%s %s\n' \
    "$(date '+%F %T')" "$1" "$$" "${run_id}" "${2:-}" > "${status_file}"
}

heartbeat() {
  while kill -0 "$$" 2>/dev/null; do
    printf '%s run_id=%s pid=%s file_bytes=%s\n' \
      "$(date '+%F %T')" "${run_id}" "$$" \
      "$(stat -c %s "${target_file}" 2>/dev/null || printf '0')" \
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
    write_status "DONE" "file=${target_file} bytes=$(stat -c %s "${target_file}")"
  else
    write_status "FAILED" "exit_code=${exit_code}"
  fi
  exit "${exit_code}"
}
trap finish EXIT INT TERM

write_status "RUNNING" "stage=download file=${target_file}"
heartbeat &
heartbeat_pid=$!

printf '%s | START GBmap Extended fetch and entry audit\n' "$(date '+%F %T')"
printf '%s | source=%s\n' "$(date '+%F %T')" "${source_url}"

# aria2c on this host cannot parse the inherited all_proxy format.
unset all_proxy ALL_PROXY

aria2c \
  --continue=true \
  --max-connection-per-server=8 \
  --split=8 \
  --min-split-size=16M \
  --max-tries=0 \
  --retry-wait=15 \
  --timeout=60 \
  --connect-timeout=30 \
  --file-allocation=none \
  --allow-overwrite=false \
  --auto-file-renaming=false \
  --dir="${data_dir}" \
  --out="$(basename "${target_file}")" \
  "${source_url}"

actual_size="$(stat -c %s "${target_file}")"
if [[ "${actual_size}" != "${expected_size}" ]]; then
  printf 'ERROR: expected %s bytes, observed %s bytes\n' "${expected_size}" "${actual_size}" >&2
  exit 1
fi

write_status "RUNNING" "stage=hdf5_structure_audit bytes=${actual_size}"
Rscript --vanilla "${project_dir}/Data_scRNA_GEO/scripts/gbmap_extended_h5ad_entry_audit.R" \
  "${target_file}" "${results_dir}"

printf '%s  %s\n' "$(sha256sum "${target_file}" | cut -d' ' -f1)" "$(basename "${target_file}")" \
  > "${results_dir}/extended_gbmap_cellxgene.sha256"

printf '%s | DONE GBmap Extended fetch and entry audit\n' "$(date '+%F %T')"
