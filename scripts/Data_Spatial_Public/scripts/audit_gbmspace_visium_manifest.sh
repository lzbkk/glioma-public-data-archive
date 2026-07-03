#!/usr/bin/env bash
set -euo pipefail

project_root="/home/lzb/glioma"
tarball="${project_root}/Data_Spatial_Public/GBM_Space/data/spatial_data_visium.tar.gz"
out_dir="${project_root}/Data_Spatial_Public/GBM_Space/results/Entry_Audit"
tables_dir="${out_dir}/tables"
logs_dir="${project_root}/Data_Spatial_Public/GBM_Space/logs"
run_id="$(date +%Y%m%d_%H%M%S)_$$"
status_file="${logs_dir}/gbmspace_manifest_audit.status"
log_file="${logs_dir}/gbmspace_manifest_audit.${run_id}.log"

expected_bytes=4106130155

mkdir -p "${tables_dir}" "${logs_dir}"

log() {
  printf '%s | %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"
}

write_status() {
  printf '%s %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" > "${status_file}"
}

{
  log "START GBM-Space Visium manifest audit run_id=${run_id}"
  write_status "RUNNING pid=$$ run_id=${run_id} step=precheck log=${log_file}"

  if [[ ! -s "${tarball}" ]]; then
    write_status "FAILED pid=$$ run_id=${run_id} reason=tarball_missing_or_empty log=${log_file}"
    log "FAILED tarball missing or empty: ${tarball}"
    exit 1
  fi

  actual_bytes="$(stat -c '%s' "${tarball}")"
  log "Tarball bytes=${actual_bytes}; expected=${expected_bytes}"
  if [[ "${actual_bytes}" != "${expected_bytes}" ]]; then
    write_status "FAILED pid=$$ run_id=${run_id} reason=byte_mismatch bytes=${actual_bytes} expected=${expected_bytes} log=${log_file}"
    log "FAILED byte mismatch"
    exit 1
  fi

  manifest_paths="${tables_dir}/gbmspace_visium_tar_paths.txt"
  manifest_verbose="${tables_dir}/gbmspace_visium_tar_verbose.txt"
  top_level="${tables_dir}/gbmspace_visium_top_level_counts.tsv"
  extension_counts="${tables_dir}/gbmspace_visium_extension_counts.tsv"
  candidate_files="${tables_dir}/gbmspace_visium_candidate_files.tsv"
  summary="${tables_dir}/gbmspace_visium_manifest_summary.tsv"

  write_status "RUNNING pid=$$ run_id=${run_id} step=list_paths log=${log_file}"
  log "Listing tar paths"
  tar -tzf "${tarball}" > "${manifest_paths}"

  write_status "RUNNING pid=$$ run_id=${run_id} step=list_verbose log=${log_file}"
  log "Listing tar verbose records"
  tar -tvzf "${tarball}" > "${manifest_verbose}"

  write_status "RUNNING pid=$$ run_id=${run_id} step=summarize log=${log_file}"
  log "Summarizing manifest"

  awk -F'/' 'BEGIN{OFS="\t"; print "top_level","n"} {k=$1; if(k=="") k="<root>"; n[k]++} END{for(k in n) print k,n[k]}' \
    "${manifest_paths}" | sort -k2,2nr > "${top_level}"

  awk '
    BEGIN{OFS="\t"; print "extension","n"}
    {
      path=$0
      n=split(path, a, "/")
      file=a[n]
      if (file == "" || file !~ /\./) ext="<none>"
      else {
        ext=file
        sub(/^.*\./, "", ext)
        ext=tolower(ext)
      }
      count[ext]++
    }
    END{for(ext in count) print ext,count[ext]}
  ' "${manifest_paths}" | sort -k2,2nr > "${extension_counts}"

  awk '
    BEGIN{
      OFS="\t"
      print "category","path"
    }
    {
      lower=tolower($0)
      category=""
      if (lower ~ /metadata|meta|annotation|annot|clinical|sample|patient|donor|section|slide/) category="metadata_or_annotation"
      if (lower ~ /h5ad|h5|hdf5|mtx|matrix|expression|counts|normalized/) category=(category=="" ? "matrix_or_expression" : category ";matrix_or_expression")
      if (lower ~ /spatial|visium|tissue_positions|scalefactors|image|histology|spot|coordinate|coords/) category=(category=="" ? "spatial_or_histology" : category ";spatial_or_histology")
      if (lower ~ /cell.*state|state.*abundance|neftel|malignant|myeloid|tme/) category=(category=="" ? "cell_state_or_tme" : category ";cell_state_or_tme")
      if (lower ~ /readme|license|manifest|json|yaml|yml/) category=(category=="" ? "docs_or_config" : category ";docs_or_config")
      if (category != "") print category,$0
    }
  ' "${manifest_paths}" > "${candidate_files}"

  {
    printf 'metric\tvalue\n'
    printf 'run_id\t%s\n' "${run_id}"
    printf 'tarball\t%s\n' "${tarball}"
    printf 'tarball_bytes\t%s\n' "${actual_bytes}"
    printf 'total_entries\t%s\n' "$(wc -l < "${manifest_paths}")"
    printf 'candidate_entries\t%s\n' "$(( $(wc -l < "${candidate_files}") - 1 ))"
    printf 'top_level_count\t%s\n' "$(( $(wc -l < "${top_level}") - 1 ))"
    printf 'generated_at\t%s\n' "$(date '+%Y-%m-%d %H:%M:%S %Z')"
    printf 'log_file\t%s\n' "${log_file}"
  } > "${summary}"

  cat > "${out_dir}/README_Entry_Audit.md" <<EOF
# GBM-Space Visium Entry Audit

生成时间：$(date '+%Y-%m-%d %H:%M:%S %Z')

## 输入

- \`${tarball}\`
- expected bytes: \`${expected_bytes}\`
- observed bytes: \`${actual_bytes}\`

## 方法概述

本审计仅扫描官方 \`spatial_data_visium.tar.gz\` 的 tar manifest，不执行完整解压。
输出包括完整路径清单、verbose 清单、顶层目录计数、扩展名计数和候选 metadata/matrix/spatial/cell-state 文件列表。

## 关键输出

- \`tables/gbmspace_visium_tar_paths.txt\`
- \`tables/gbmspace_visium_tar_verbose.txt\`
- \`tables/gbmspace_visium_top_level_counts.tsv\`
- \`tables/gbmspace_visium_extension_counts.tsv\`
- \`tables/gbmspace_visium_candidate_files.tsv\`
- \`tables/gbmspace_visium_manifest_summary.tsv\`

## 第一版结果摘要

见 \`tables/gbmspace_visium_manifest_summary.tsv\`、\`tables/gbmspace_visium_top_level_counts.tsv\`
和 \`tables/gbmspace_visium_candidate_files.tsv\`。

## 解释边界

该步骤只判断包内结构与后续可用入口，不判断 LAP3 空间邻域结果。
正式空间分析必须等待 LAP3 gene coverage、spot 坐标、patient/section 层级、
cell-state abundance 和 histopathology annotation 字段确认后再启动。

## 下一步

基于候选文件列表选择性解压 metadata / feature names / spot annotation，
建立 lightweight cache，再决定是否开展患者/切片分层的 LAP3-neighborhood 分析。
EOF

  write_status "DONE pid=$$ run_id=${run_id} entries=$(wc -l < "${manifest_paths}") log=${log_file} result=${out_dir}"
  log "DONE entries=$(wc -l < "${manifest_paths}")"
} >> "${log_file}" 2>&1
