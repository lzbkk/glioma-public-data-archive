#!/usr/bin/env bash
set -uo pipefail

ROOT="/home/lzb/glioma/Data_Spatial_Public/GBM_Space"
DATA_DIR="${ROOT}/data"
LOG_DIR="${ROOT}/logs"
OUT="${DATA_DIR}/spatial_data_visium.tar.gz"
PART="${OUT}.part"
LOG="${LOG_DIR}/gbmspace_visium_download.log"
STATUS="${LOG_DIR}/gbmspace_visium_download.status"
URL="https://gbmspace.cog.sanger.ac.uk/spatial_data_visium.tar.gz"
EXPECTED_BYTES=4106130155
RUN_ID="$(date +%Y%m%d_%H%M%S)_$$"

mkdir -p "${DATA_DIR}" "${LOG_DIR}"

write_status() {
  printf '%s %s pid=%s run_id=%s file=%s bytes=%s expected=%s log=%s\n' \
    "$(date '+%Y-%m-%d %H:%M:%S')" "$1" "$$" "${RUN_ID}" "${OUT}" "$2" \
    "${EXPECTED_BYTES}" "${LOG}" > "${STATUS}"
}

log_msg() {
  printf '%s | %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" | tee -a "${LOG}"
}

: > "${LOG}"
current_bytes=0
[[ -f "${PART}" ]] && current_bytes="$(stat -c %s "${PART}")"
write_status "RUNNING" "${current_bytes}"
log_msg "START GBM-Space Visium download run_id=${RUN_ID}"

if [[ -f "${OUT}" ]] && [[ "$(stat -c %s "${OUT}")" -eq "${EXPECTED_BYTES}" ]]; then
  log_msg "Existing validated file found"
  write_status "DONE" "${EXPECTED_BYTES}"
  exit 0
fi

success=0
for attempt in 1 2 3 4 5 6; do
  current_bytes=0
  [[ -f "${PART}" ]] && current_bytes="$(stat -c %s "${PART}")"
  write_status "RUNNING attempt=${attempt}/6" "${current_bytes}"
  log_msg "Attempt ${attempt}/6; current_bytes=${current_bytes}"

  if curl -fL -C - --connect-timeout 30 --max-time 7200 \
    --speed-time 180 --speed-limit 1024 \
    -o "${PART}" "${URL}" >> "${LOG}" 2>&1; then
    downloaded_bytes="$(stat -c %s "${PART}")"
    if [[ "${downloaded_bytes}" -eq "${EXPECTED_BYTES}" ]]; then
      mv "${PART}" "${OUT}"
      success=1
      break
    fi
    log_msg "Byte check failed: expected=${EXPECTED_BYTES} observed=${downloaded_bytes}"
  else
    log_msg "curl attempt ${attempt} failed; retaining partial file"
  fi
  sleep 30
done

if [[ "${success}" -eq 1 ]]; then
  log_msg "DONE bytes=${EXPECTED_BYTES}"
  write_status "DONE" "${EXPECTED_BYTES}"
  exit 0
fi

current_bytes=0
[[ -f "${PART}" ]] && current_bytes="$(stat -c %s "${PART}")"
log_msg "FAILED bytes=${current_bytes}"
write_status "FAILED" "${current_bytes}"
exit 1
