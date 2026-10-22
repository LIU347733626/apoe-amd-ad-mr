#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""build_na_from_v66.py — 从 NC 定稿 v66 重建 Nature Aging 投稿版（2026-09-10）

背景：NA_package 的手稿（manuscript_NA_v1..v4）派生自工作稿 v47，其后未跟进 v48–v66
      的任何修改（含 ScPagwas 撤回、Figure 5 换机制图、S8 换亚型图、共定位撤销、A1 修复等）。
      PI 决定：NA 版整体重建，手稿与图同步到 v66 定稿。

NA 格式要求（NA_package/admin/NA_requirements.md）：
  · 摘要 ≤ 150 词、不加引用
  · Introduction 不设标题；正文结构 = Introduction / Results / Discussion / Online Methods
  · 正文 ≤ 4,500 词（不含摘要、Methods、refs、图注）
  · 补充图引用写作 "Figure Sx"（本构建统一转换自 v66 的 "Supplementary Fig. x"）

输入 : final_analysis/00_manuscript_NC_v66_20260910.md
输出 : NA_package/manuscript/manuscript_NA_v5.md
"""
import io, os, re, sys

PROJ = "/home/huan/projects/08_APOE_AMD_AD_MR"
SRC  = os.path.join(PROJ, "final_analysis", "00_manuscript_NC_v66_20260910.md")
DST  = os.path.join(PROJ, "NA_package", "manuscript", "manuscript_NA_v5.md")

NA_ABSTRACT = (
"**Background:** Alzheimer's disease (AD) and age-related macular degeneration (AMD) are leading "
"degenerative diseases of ageing, yet the APOE ε4 allele raises AD risk while protecting against AMD. "
"**Methods:** Using cis-pQTL instruments for plasma APOE, we tested by Mendelian randomization whether "
"APOE protein abundance, not isoform function, is the shared causal variable. "
"**Results:** ε4 lowered APOE abundance (β = −1.012 per allele) and ε2 raised it (+0.745). Predicted "
"APOE abundance increased AMD risk (β = +0.182, P = 7.9 × 10⁻¹¹) and decreased AD risk "
"(β = −0.942, P = 4.1 × 10⁻¹⁸), with 90.1% (AMD) and 95.8% (AD) of the ε4 effect level-mediated, "
"replicating across cohorts. Retinal eQTL MR was null, and single-cell analyses localised the switch "
"to microglial TREM2 signalling (AD) and lipid-transport pathways (AMD). Complement was coregulated "
"with, not driven by, APOE. "
"**Conclusions:** APOE abundance is a systemic state that protects the brain while endangering the retina."
)

t = io.open(SRC, encoding="utf-8").read()
orig = t
log = []

# ---- 1) 去掉 NC 的版本注记（NA 版不保留内部 changelog）----
lines = [l for l in t.split("\n") if not l.startswith("> **Version note")]
t = "\n".join(lines)
log.append("移除 NC 版本注记")

# ---- 2) 替换摘要为 NA 版（≤150 词）----
i = t.index("## Abstract")
j = t.index("## Introduction")
kw_start = t.rindex("**Keywords:**", i, j) if "**Keywords:**" in t[i:j] else j
old_kw = t[kw_start:j].strip()
new_block = "## Abstract\n\n" + NA_ABSTRACT + "\n\n" + old_kw + "\n\n"
t = t[:i] + new_block + t[j:]
log.append("摘要替换为 NA 版（≤150 词）")

# ---- 3) Introduction 去标题（NA：引言不设标题）----
t = t.replace("## Introduction\n", "", 1)
log.append("Introduction 去标题")

# ---- 4) Methods → Online Methods ----
t = t.replace("## Methods\n", "## Online Methods\n", 1)
log.append("Methods → Online Methods")

# ---- 5) 补充图引用：Supplementary Fig. x → Figure Sx（NA 风格）----
n1 = len(re.findall(r"Supplementary Fig\. ", t))
t = re.sub(r"Supplementary Fig\. (\d+)", r"Figure S\1", t)
log.append(f"补充图引用转换 Supplementary Fig. → Figure S（{n1} 处）")

# ---- 6) 标题下的 running title / 关键词保持 ----
t = t.rstrip("\n") + "\n"

# ---- 7) 门禁：词数 ----
def wc(s):
    s = re.sub(r"^#+ .*$", "", s, flags=re.M)
    s = re.sub(r"^\|.*$", "", s, flags=re.M)          # 去表格行
    s = re.sub(r"^> .*$", "", s, flags=re.M)           # 去引用块
    return len(re.findall(r"\S+", s))

abs_txt = t[t.index("## Abstract"):t.index("**Keywords:**", t.index("## Abstract"))]
abs_words = wc(abs_txt)
# 正文 = Introduction + Results + Discussion（不含摘要、Online Methods、参考文献、图注）
_a = t.index("**Keywords:**", t.index("## Abstract"))
_b = t.index("## Online Methods")
main_words = wc(t[_a:_b])
assert abs_words <= 150, f"摘要 {abs_words} 词 > 150，需再压缩"
log.append(f"摘要 {abs_words} 词（限 150）✓")
log.append(f"正文（含摘要与图注，不含 Methods/参考文献）≈ {main_words} 词（限 4500）")

io.open(DST, "w", encoding="utf-8").write(t)
print("已写:", DST)
for l in log: print("  ✓", l)
