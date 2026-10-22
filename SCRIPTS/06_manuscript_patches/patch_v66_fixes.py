#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""patch_v66_fixes.py — 把本轮修复的结果写回手稿（2026-09-10）

依据：docs/repro_audit_20260910/FINDINGS_MASTER_20260910.md
  1) A1 修复后 IAMDGC 复制变为 11 工具（原 :110 的 +0.211/+0.212 是 header bug 下的 10 工具值）；
  2) S-1 Bellenguez 发布文件的 per-variant N 上限与论文声称不一致，须如实标注；
  3) §11 新增官方 11 工具集的 leave-one-out 稳健性（直接回应"是否由单个工具驱动"）。

输入 : final_analysis/00_manuscript_NC_v65_20260910.md
输出 : final_analysis/00_manuscript_NC_v66_20260910.md
"""
import io, os

PROJ = "/home/huan/projects/08_APOE_AMD_AD_MR"
SRC = os.path.join(PROJ, "final_analysis", "00_manuscript_NC_v65_20260910.md")
DST = os.path.join(PROJ, "final_analysis", "00_manuscript_NC_v66_20260910.md")

t = io.open(SRC, encoding="utf-8").read()
orig = t
log = []


def sub_once(old, new, label):
    n = t.count(old)
    assert n == 1, f"[{label}] 锚点出现 {n} 次（应为 1）: {old[:70]!r}"
    return t.replace(old, new)


# --------------------------------------------------------------------------
# F1  IAMDGC 复制值：10 工具（header bug）→ 11 工具（修复后实测）
# --------------------------------------------------------------------------
t = sub_once(
    "replicated the APOE→AMD effect with the same clumped instruments: GSMR β = +0.211 (SE = 0.035, "
    "P = 1.3 × 10⁻⁹), IVW β = +0.212 (P = 4.3 × 10⁻⁷), MR-Egger intercept P = 0.67 (no directional "
    "pleiotropy).",
    "replicated the APOE→AMD effect with the same clumped instruments: GSMR β = +0.212 (SE = 0.035, "
    "P = 8.7 × 10⁻¹⁰), IVW β = +0.213 (P = 8.5 × 10⁻⁸), MR-Egger intercept P = 0.59 (no directional "
    "pleiotropy).",
    "F1")
log.append("F1 IAMDGC 复制值更新为 11 工具（修复 header bug 后）")

# --------------------------------------------------------------------------
# F2  Bellenguez：如实标注发布文件的 per-variant N 上限
# --------------------------------------------------------------------------
t = sub_once(
    "and Bellenguez 2022[22] (AD; 111,326 cases + 677,663 controls; European ancestry; contains the "
    "FinnGen component; APOE ε4/ε2 coding variants absent from the released summary statistics).",
    "and Bellenguez 2022[22] (AD; European ancestry; contains the FinnGen component; APOE ε4/ε2 coding "
    "variants absent from the released summary statistics). The published case–control counts of that "
    "meta-analysis (111,326 cases + 677,663 controls) are not recoverable from the released summary "
    "statistics: the GCST90027158 file reports per-variant n = n_cases + n_controls, whose maximum is "
    "487,511 (85,934 cases + 401,577 controls), and 12 classical APOE-region tag variants are absent "
    "from it altogether. This cohort therefore enters only as a directional validation, and no "
    "case-fraction-dependent quantity is derived from it.",
    "F2")
log.append("F2 Bellenguez 样本量：标注文件实际 per-variant N 上限 487,511")

# --------------------------------------------------------------------------
# F3  新增 leave-one-out 稳健性（官方 11 工具集）
# --------------------------------------------------------------------------
anchor = "**Instrument robustness.**"
i = t.find(anchor)
assert i >= 0, "[F3] 未找到 Instrument-definition robustness 段"
seg = ("**Instrument-leverage check (leave-one-out, official clumped set).** The 11-instrument clumped "
       "APOE pQTL set (which by design excludes the ε4 coding tag rs429358) was re-estimated with each "
       "instrument removed in turn. The AMD estimate is not carried by any single instrument: β ranged "
       "0.154–0.214 across the 11 leave-one-out fits, with a maximum P of 1.3 × 10⁻⁶ (vs β = +0.182, "
       "P = 7.9 × 10⁻¹¹ for the full set), and the heterogeneity statistic is non-significant "
       "(Q = 11.0, df = 10, P = 0.36). The AD estimate likewise held (β = −1.008 to −0.862, maximum "
       "P = 1.1 × 10⁻¹¹), although instrument heterogeneity there is substantial (Q P ≈ 3 × 10⁻⁴²), so "
       "the AD interval is reported under the overdispersion-corrected convention that "
       "MendelianRandomization::mr_ivw applies by default; the classical fixed-effect interval is "
       "narrower (SE 0.0230 vs 0.1086) with the same point estimate, and both are given in Methods. "
       "The IAMDGC replication behaves the same way (β 0.181–0.268, maximum P = 9.8 × 10⁻⁶). A separate "
       "leave-one-out on instrument sets that *include* rs429358 has been reported by others to lose "
       "significance for AMD when that variant is removed; that set is a different exposure definition "
       "(the ε4 coding tag is the single strongest level regulator) and is not the instrument set used "
       "here, which is why the two analyses are not directly comparable.\n\n")
t = t[:i] + seg + t[i:]
log.append("F3 新增官方 11 工具集的 leave-one-out 稳健性段")

# --------------------------------------------------------------------------
# 断言
# --------------------------------------------------------------------------
assert "GSMR β = +0.211 (SE = 0.035, P = 1.3 × 10⁻⁹)" not in t, "旧 IAMDGC 值仍在"
body = t[:t.index("## Figure Legends")]
body_only = "\n".join(l for l in body.split("\n") if not l.startswith("> **Version note"))
for banned in ["PP3 = 1.0", "1.16 × 10⁻²⁰⁵"]:
    assert banned not in body_only, f"正文仍含: {banned}"
log.append("断言通过")

note = (
    "> **Version note (2026-09-10, v66, fixes arising from the reproducibility audit):** Three changes. "
    "(i) The IAMDGC replication estimates were recomputed: the script that produced them passed "
    "`header = TRUE` to a `grep` extraction whose output contains no header, so the first matching "
    "variant was silently consumed as column names and **one instrument (rs28399649) was dropped** "
    "(11 → 10). With the fix, GSMR β = +0.212 (SE 0.035, P = 8.7 × 10⁻¹⁰) and IVW β = +0.213 "
    "(P = 8.5 × 10⁻⁸), which now agrees with the independent IAMDGC line already reported in the "
    "Results. (ii) The Bellenguez 2022 cohort description no longer implies that its published "
    "case–control counts are available in the released summary statistics: the GCST90027158 file "
    "reports per-variant n ≤ 487,511 and lacks 12 classical APOE-region tag variants, so the cohort is "
    "a directional validation only. (iii) A leave-one-out leverage check on the official 11-instrument "
    "set was added; it shows the AMD estimate is not driven by any single instrument. **No other "
    "estimate, P value, confidence interval or sample size was changed.** Audit records: "
    "`docs/repro_audit_20260910/FINDINGS_MASTER_20260910.md`. v65 preserved as "
    "`00_manuscript_NC_v65_20260910.md`.\n\n"
)
t = t.replace("# APOE abundance drives inverse comorbidity",
              note + "# APOE abundance drives inverse comorbidity", 1)

assert t != orig, "无改动"
io.open(DST, "w", encoding="utf-8").write(t)
print("已写:", DST)
for l in log:
    print("  ✓", l)
print(f"行数 {orig.count(chr(10))+1} → {t.count(chr(10))+1}")
