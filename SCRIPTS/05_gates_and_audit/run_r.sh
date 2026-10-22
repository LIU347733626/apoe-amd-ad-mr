#!/usr/bin/env bash
# run_r.sh — R 执行 + 资源预检（线程自动拉满）+ 自动图质检门禁 (figgate)
# 用法: bash run_r.sh < 脚本.R   或   bash run_r.sh -e 'R代码'
# 关闭质检: FIGGATE_ON=0
#
# 2026-09-09：本脚本现在会 source scripts/preflight.sh，自动导出
#   OMP/OPENBLAS/MKL/R_DATATABLE_NUM_THREADS 与 RES_CORES/RES_JOBS，
#   并打印一行资源横幅——这样"用上本机多核"不再依赖备注被记住。
cd "$(dirname "$0")"
if [ -f scripts/preflight.sh ]; then
    # shellcheck disable=SC1091
    source scripts/preflight.sh
fi
export RSCRIPT="${RSCRIPT:-/home/huan/software/miniconda3/envs/biomni_r/bin/Rscript}"
export FIGGATE_DIRS="${FIGGATE_DIRS:-figures manuscript revision}"
exec bash /home/huan/.agent-vision/rgate.sh "$@"
