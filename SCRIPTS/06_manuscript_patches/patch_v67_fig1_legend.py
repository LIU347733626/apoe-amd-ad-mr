#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""patch_v67_fig1_legend.py — Figure 1 重构后的图注联动（2026-09-10）

Figure 1 已重构两处：
  (1) 图例补齐 Second validation，并把「颜色=队列 / 形状=APOE」两个维度拆开
      （原先 APOE 三行统一紫色菱形 → 队列信息丢失，且 Second validation 无图例项）；
  (2) Second validation 两层由 IVW 改为 GSMR（`revision/98_second_validation_gsmr.R` 补跑），
      全图六行估计量统一为 GSMR。
原图注写 "APOE is highlighted (purple diamonds, bold)"，与 (1) 不再相符，且未声明估计量。
"""
import io, os
PROJ = "/home/huan/projects/08_APOE_AMD_AD_MR"
OLD = "APOE is highlighted (purple diamonds, bold)."
NEW = ("APOE is highlighted by marker shape (a wide diamond, whereas all other proteins are filled "
       "squares); marker colour denotes the cohort tier throughout (blue = Discovery, orange = First "
       "validation, green = Second validation), and the tier is also named in the Cohort column. "
       "All estimates in this figure are GSMR estimates (gsmr 1.0.6) from the same clumped APOE "
       "cis-pQTL instrument set; the number of instruments retained after HEIDI filtering is reported "
       "per row in the source data.")
for src, dst in [("final_analysis/00_manuscript_NC_v66_20260910.md",
                  "final_analysis/00_manuscript_NC_v67_20260910.md"),
                 ("NA_package/manuscript/manuscript_NA_v5.md",
                  "NA_package/manuscript/manuscript_NA_v6.md")]:
    p = os.path.join(PROJ, src); t = io.open(p, encoding="utf-8").read()
    n = t.count(OLD)
    assert n == 1, f"{src}: 锚点 {n} 次"
    t = t.replace(OLD, NEW)
    io.open(os.path.join(PROJ, dst), "w", encoding="utf-8").write(t)
    print(f"  ✓ {os.path.basename(dst)}（图注已更新）")
