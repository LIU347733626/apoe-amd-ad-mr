#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""patch_v65_colocp.py — 修正 v61 引入的 p 值机制表述错误，并修正 LD 措辞（2026-09-10）

背景（对抗式复核的结果）：
  `patch_v61.py` 正确撤销了"多信号共定位"层，但我在 Methods 与版本说明里写的**机制解释是错的**：
  我写"从 β/SE 重算 p 从而保留这些变异"。实测证明重算 p **救不回**这些变异——
  被删的 31 个变异 |z| = 41.8–73.1，而 `2*pnorm(-|z|)` 在双精度下 |z| ≳ 37.5 即下溢为 0，
  重算后仍为 0（实测：p_raw==0 的 31 个，p_recomp==0 的也是 31 个）。
  真正起作用的是**是否施加 `p > 0` 过滤器**，而不是 p 是否重算。
  进一步：`coloc:::process.dataset` 在同时提供 beta 与 varbeta 时走 `approx.bf.estimates`，
  **p 列完全不参与似然**（本机实测源码分支；红队用"把 p 全改 1e-300 结果不变"独立佐证）。

  另修正一处措辞：Results 段称 AD credible set 含 rs429358 的 "near-perfect-LD neighbors"，
  但实测 r²(rs429358, rs405509) = 0.084、r²(rs429358, rs440446) = 0.092（EU1000QC），
  只有 rs769449 (0.771) 接近高 LD —— "near-perfect" 不成立。

输入 : final_analysis/00_manuscript_NC_v64_20260910.md
输出 : final_analysis/00_manuscript_NC_v65_20260910.md
"""
import io, os

PROJ = "/home/huan/projects/08_APOE_AMD_AD_MR"
SRC = os.path.join(PROJ, "final_analysis", "00_manuscript_NC_v64_20260910.md")
DST = os.path.join(PROJ, "final_analysis", "00_manuscript_NC_v65_20260910.md")

t = io.open(SRC, encoding="utf-8").read()
orig = t
log = []


def sub_once(old, new, label):
    n = t.count(old)
    assert n == 1, f"[{label}] 锚点出现 {n} 次（应为 1）: {old[:70]!r}"
    return t.replace(old, new)


# --------------------------------------------------------------------------
# F1  Methods：把错误的"重算 p 保留变异"改为准确的机制表述
# --------------------------------------------------------------------------
t = sub_once(
    "In all runs the disease p values were recomputed from the reported β/SE within the region rather "
    "than taken from the released p column, because that column is truncated at p = 0 for the strongest "
    "APOE-region variants (31 such variants in the FinnGen R13 AD file, including rs429358); retaining "
    "them is what the estimates above are based on.",
    "The released p column is truncated at 0 for the strongest APOE-region variants (31 such variants in "
    "the FinnGen R13 AD file, including rs429358; 30 of them survive the LD-reference filter applied "
    "here). Those variants carry |z| = 41.8–73.1, and a p value recomputed from β/SE underflows to zero "
    "as well in double precision (|z| ≳ 37.5), so recomputing p does not recover them. The analyses above "
    "therefore do not filter on p at all: variants are retained on finite β and positive SE. This matters "
    "because a `p > 0` filter removes precisely the strongest signals in the region. In this configuration "
    "coloc.abf derives its approximate Bayes factors from β and varbeta (both supplied), so the p column "
    "does not enter the likelihood in any case.",
    "F1")
log.append("F1 Methods：机制表述纠正（重算 p 无效 → 必须不按 p 过滤；并说明 coloc.abf 用 β/varbeta）")

# --------------------------------------------------------------------------
# F2  版本说明：同一错误的另一处
# --------------------------------------------------------------------------
t = sub_once(
    "including rs429358 — recomputing p from β/SE restores them and returns PP4 ≈ 1; and (iii)",
    "including rs429358 (|z| = 41.8–73.1, so a recomputed p underflows to zero as well); **removing that "
    "filter** returns PP4 ≈ 1, whereas recomputing p while keeping the filter does not; and (iii)",
    "F2")
log.append("F2 版本说明：同步纠正")

# --------------------------------------------------------------------------
# F3  Results：LD 措辞（"near-perfect-LD neighbors" 与实测 r² 不符）
# --------------------------------------------------------------------------
t = sub_once(
    "comprising the ε4 coding tag rs429358 and its near-perfect-LD neighbors (rs405509, rs440446, "
    "rs769449, rs10414043 — the tag's own PIP is ≈ 0 in the R13 AD run)",
    "comprising the ε4 coding tag rs429358 together with rs405509, rs440446, rs769449 and rs10414043 "
    "(the tag's own PIP is ≈ 0 in the R13 AD run). The LD between the tag and these four variants in the "
    "1000 Genomes European reference is heterogeneous — r² = 0.77 with rs769449 of order high LD, but only "
    "0.08–0.09 with rs405509 and rs440446 — so this credible set spans more than one LD block rather than "
    "a single near-perfect-LD cluster",
    "F3")
log.append("F3 Results：LD 措辞与实测 r² 对齐")

# --------------------------------------------------------------------------
# F4  同一句群里的"near-perfect LD"泛化表述
# --------------------------------------------------------------------------
t = sub_once(
    "the single-signal PP4 for AD (1.0) reflects near-perfect LD between the ε4 cluster and the level "
    "instruments and does not by itself distinguish shared from independent causality",
    "the single-signal PP4 for AD (1.0) reflects the strong LD between the ε4 cluster and the level "
    "instruments and does not by itself distinguish shared from independent causality",
    "F4")
log.append("F4 Results：去掉 'near-perfect' 的过度表述")

# --------------------------------------------------------------------------
# 验证：被撤数字仍不在正文；版本说明里只应出现在记录语境
# --------------------------------------------------------------------------
body = t[:t.index("## Figure Legends")]
body_only = "\n".join(l for l in body.split("\n") if not l.startswith("> **Version note"))
for banned in ["PP3 = 1.0", "1.16 × 10⁻²⁰⁵", "Multi-signal (suSiE-based) colocalization resolves"]:
    assert banned not in body_only, f"正文仍含被撤表述: {banned}"
log.append("断言：被撤表述仍不在正文")

note = (
    "> **Version note (2026-09-10, v65, correction of the p-value mechanism statement):** An adversarial "
    "re-check of the v61 withdrawal found that the **mechanism** I stated there was wrong. v61 said the "
    "31 affected variants could be recovered by recomputing p from β/SE; measurement shows they cannot — "
    "they carry |z| = 41.8–73.1 and a recomputed p underflows to zero as well in double precision "
    "(|z| ≳ 37.5). The determinant is the `p > 0` filter itself, not how p was obtained, and in this "
    "configuration `coloc.abf` uses β and varbeta rather than the p column at all. The Methods and the "
    "v61 note have been corrected to say that variants are retained on finite β and positive SE, and that "
    "30 (not 31) of them survive the LD-reference filter. Separately, the Results description of the AD "
    "credible set called rs405509/rs440446 'near-perfect-LD neighbours' of rs429358, whereas the measured "
    "r² is 0.084–0.092 (rs769449 is 0.77); the wording now reports the heterogeneous LD. **No estimate, "
    "P value, confidence interval or sample size was changed; no claim was added or removed.** v64 "
    "preserved as `00_manuscript_NC_v64_20260910.md`.\n\n"
)
t = t.replace("# APOE abundance drives inverse comorbidity", note + "# APOE abundance drives inverse comorbidity", 1)

assert t != orig, "无改动"
io.open(DST, "w", encoding="utf-8").write(t)
print("已写:", DST)
for l in log:
    print("  ✓", l)
print(f"行数 {orig.count(chr(10))+1} → {t.count(chr(10))+1}")
