#!/usr/bin/env bash
set -euo pipefail

cd /home/lzb/glioma

root="Data_Perturbation_Public/LINCS_GSE70138"
result_dir="${root}/results/LAP3_Local_Connectivity"
log_dir="${root}/logs"
status="${log_dir}/lap3_local_connectivity.status"
run_id="$(date '+%Y%m%d_%H%M%S')_$$"
run_log="${log_dir}/lap3_local_connectivity.${run_id}.log"
latest_log="${log_dir}/lap3_local_connectivity.log"

mkdir -p "${result_dir}" "${log_dir}"
ln -sfn "$(basename "${run_log}")" "${latest_log}"
exec >>"${run_log}" 2>&1

on_exit() {
  exit_code=$?
  if [[ $exit_code -ne 0 ]]; then
    date "+%F %T FAILED pid=$$ run_id=$run_id exit_code=$exit_code log=$run_log" >"${status}"
    date "+%F %T | FAILED LAP3 local LINCS connectivity exit_code=$exit_code"
  fi
}
trap on_exit EXIT

date "+%F %T RUNNING pid=$$ run_id=$run_id stage=initialize workers=16 log=$run_log" >"${status}"
date "+%F %T | START LAP3 local LINCS connectivity run_id=$run_id"

OMP_NUM_THREADS=1 \
OPENBLAS_NUM_THREADS=1 \
MKL_NUM_THREADS=1 \
LAP3_LINCS_OUT_DIR="${result_dir}" \
LAP3_LINCS_STATUS="${status}" \
LAP3_LINCS_MAX_PROFILES=0 \
LAP3_LINCS_WORKERS=16 \
LAP3_LINCS_CHUNK_SIZE=512 \
Rscript --vanilla "${root}/scripts/lap3_local_lincs_connectivity.R"

test -s "${result_dir}/source_data/lap3_lincs_signature_level_scores.csv.gz"
test -s "${result_dir}/tables/lap3_lincs_compound_summary_all.csv"
test -s "${result_dir}/tables/lap3_lincs_robust_reverse_candidates.csv"
test -s "${result_dir}/tables/lap3_lincs_full_vs_no_lap3_compounds.csv"
test -s "${result_dir}/tables/lap3_lincs_run_qc.csv"
gzip -t "${result_dir}/source_data/lap3_lincs_signature_level_scores.csv.gz"

date "+%F %T DONE pid=$$ run_id=$run_id profiles=118050 workers=16 failed=0 log=$run_log" >"${status}"
date "+%F %T | DONE LAP3 local LINCS connectivity"
