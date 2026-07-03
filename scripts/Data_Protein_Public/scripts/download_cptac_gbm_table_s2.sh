#!/usr/bin/env bash
set -uo pipefail

ROOT="/home/lzb/glioma/Data_Protein_Public"
DATA_DIR="${ROOT}/data/CPTAC_GBM_Publication"
LOG_DIR="${ROOT}/logs"
OUT="${DATA_DIR}/Wang2021_CancerCell_Table_S2_mmc3.xlsx"
PART="${OUT}.part"
LOG="${LOG_DIR}/cptac_gbm_table_s2_download.log"
STATUS="${LOG_DIR}/cptac_gbm_table_s2_download.status"
URL="https://ars.els-cdn.com/content/image/1-s2.0-S1535610821000507-mmc3.xlsx"
EXPECTED_BYTES=129239538
RUN_ID="$(date +%Y%m%d_%H%M%S)_$$"

mkdir -p "${DATA_DIR}" "${LOG_DIR}"

write_status() {
  printf '%s %s pid=%s run_id=%s file=%s bytes=%s log=%s\n' \
    "$(date '+%Y-%m-%d %H:%M:%S')" "$1" "$$" "${RUN_ID}" "${OUT}" "$2" "${LOG}" \
    > "${STATUS}"
}

log_msg() {
  printf '%s | %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" | tee -a "${LOG}"
}

: > "${LOG}"
current_bytes=0
[[ -f "${PART}" ]] && current_bytes="$(stat -c %s "${PART}")"
write_status "RUNNING" "${current_bytes}"
log_msg "START CPTAC GBM Table S2 download run_id=${RUN_ID}"

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

  if curl -fL -C - --connect-timeout 30 --max-time 1800 \
    --speed-time 120 --speed-limit 1024 \
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
  sleep 15
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
