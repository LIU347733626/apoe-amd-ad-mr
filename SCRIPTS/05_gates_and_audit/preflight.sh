#!/usr/bin/env bash
# preflight.sh — 资源预检 + 线程环境（2026-09-09 建立）
#
# 目的：把"本机有多少核/多少内存、该怎么并行"从**备注**变成**可执行的一步**，
#       避免每个窗口重新摸索、或者忘了用多核而拖慢。
#
# 用法：
#   bash scripts/preflight.sh            # 打印一行摘要 + 建议
#   source scripts/preflight.sh          # 导出线程环境变量（供同 shell 后续命令用）
#   bash scripts/preflight.sh --gate     # 额外跑 llm_gate（开子代理/长任务前）
#   bash scripts/preflight.sh --json     # 机器可读
#
# 导出（source 时）：
#   RES_CORES  物理核数
#   RES_JOBS   建议并行任务数（= 核数）
#   OMP_NUM_THREADS / OPENBLAS_NUM_THREADS / MKL_NUM_THREADS / R_DATATABLE_NUM_THREADS
#
# 纪律（写进 CLAUDE.md 顶部「资源铁律」）：
#   - 重活（>1 min 或 >1 GB 输入）前先 source 本脚本；
#   - R 多任务用 mclapply(..., mc.cores = min(n_tasks, RES_JOBS))；
#   - shell 多任务用 xargs -P "$RES_JOBS"；
#   - R 脚本一律走 `bash run_r.sh`（它会自动 source 本脚本并打印横幅）。

set -u
JSON=0; GATE=0
for a in "$@"; do
  case "$a" in
    --json) JSON=1 ;;
    --gate) GATE=1 ;;
  esac
done

RES_CORES="$(nproc 2>/dev/null || echo 1)"
RES_JOBS="$RES_CORES"
# BLAS/OpenMP 超过 ~32 线程收益很小且易互相抢核
_res_blas="$RES_CORES"; [ "$_res_blas" -gt 32 ] && _res_blas=32
export OMP_NUM_THREADS="$RES_CORES"
export OPENBLAS_NUM_THREADS="$_res_blas"
export MKL_NUM_THREADS="$_res_blas"
export R_DATATABLE_NUM_THREADS="$RES_CORES"
export RES_CORES RES_JOBS

_res_mem_gb="$(awk '/MemTotal/{printf "%d", $2/1048576}' /proc/meminfo 2>/dev/null || echo 0)"
_res_avail_gb="$(awk '/MemAvailable/{printf "%d", $2/1048576}' /proc/meminfo 2>/dev/null || echo 0)"
_res_load="$(awk '{printf "%.1f", $1}' /proc/loadavg 2>/dev/null || echo '?')"

if [ "$JSON" = "1" ]; then
  printf '{"cores":%s,"jobs":%s,"mem_total_gb":%s,"mem_avail_gb":%s,"load1":%s,"omp":%s,"blas":%s,"dt":%s}\n' \
    "$RES_CORES" "$RES_JOBS" "$_res_mem_gb" "$_res_avail_gb" "$_res_load" \
    "$OMP_NUM_THREADS" "$OPENBLAS_NUM_THREADS" "$R_DATATABLE_NUM_THREADS"
else
  echo "[preflight] cores=$RES_CORES jobs=$RES_JOBS mem=${_res_mem_gb}G(avail ${_res_avail_gb}G) load1=$_res_load | threads: omp=$OMP_NUM_THREADS blas=$OPENBLAS_NUM_THREADS dt=$R_DATATABLE_NUM_THREADS"
  echo "[preflight] R 多任务: mclapply(seq_len(n), f, mc.cores=min(n, $RES_JOBS)) | shell 多任务: xargs -P $RES_JOBS | R 一律走: bash run_r.sh < script.R"
fi

if [ "$GATE" = "1" ] && [ -x "$HOME/.agent-tools/llm_gate.sh" ]; then
  bash "$HOME/.agent-tools/llm_gate.sh" 2>&1 | tail -3
fi
