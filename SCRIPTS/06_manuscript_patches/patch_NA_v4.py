#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""patch_NA_v4.py — 把 v66 的 A1 修复结果与 leave-one-out 同步到 NA 投稿版（2026-09-10）

与 v66 的差异：NA 版的 Bellenguez 描述未声称 111,326+677,663，故不需要 F2。
输入 : NA_package/manuscript/manuscript_NA_v3.md
输出 : NA_package/manuscript/manuscript_NA_v4.md
"""
import io, os

PROJ = "/home/huan/projects/08_APOE_AMD_AD_MR"
SRC = os.path.join(PROJ, "NA_package", "manuscript", "manuscript_NA_v3.md")
DST = os.path.join(PROJ, "NA_package", "manuscript", "manuscript_NA_v4.md")

t = io.open(SRC, encoding="utf-8").read()
orig = t
log = []


def sub_once(old, new, label):
    n = t.count(old)
    assert n == 1, f"[{label}] 锚点出现 {n} 次（应为 1）: {old[:70]!r}"
    return t.replace(old, new)


# F1: IAMDGC 值（10 工具 header bug → 11 工具）
t = sub_once(
    "replicated the APOE→AMD effect with the same clumped instruments: GSMR β = +0.211 (SE = 0.035, "
    "P = 1.3 × 10⁻⁹), IVW β = +0.212 (P = 4.3 × 10⁻⁷), MR-Egger intercept P = 0.67 (no directional "
    "pleiotropy).",
    "replicated the APOE→AMD effect with the same clumped instruments: GSMR β = +0.212 (SE = 0.035, "
    "P = 8.7 × 10⁻¹⁰), IVW β = +0.213 (P = 8.5 × 10⁻⁸), MR-Egger intercept P = 0.59 (no directional "
    "pleiotropy).",
    "F1")
log.append("F1 IAMDGC 值更新为 11 工具")

# F3: leave-one-out 稳健性
anchor = "**Instrument robustness.**"
i = t.find(anchor)
assert i >= 0, "[F3] 未找到 Instrument robustness"
seg = ("**Instrument-leverage check (leave-one-out, official clumped set).** The 11-instrument clumped "
       "APOE pQTL set (which by design excludes the ε4 coding tag rs429358) was re-estimated with each "
       "instrument removed in turn. The AMD estimate is not carried by any single instrument: β ranged "
       "0.154–0.214 across the 11 leave-one-out fits, with a maximum P of 1.3 × 10⁻⁶ (vs β = +0.182, "
       "P = 7.9 × 10⁻¹¹ for the full set), and the heterogeneity statistic is non-significant "
       "(Q = 11.0, df = 10, P = 0.36). The AD estimate likewise held (β = −1.008 to −0.862, maximum "
       "P = 1.1 × 10⁻¹¹), although instrument heterogeneity there is substantial, so the AD interval is "
       "reported under the overdispersion-corrected convention that MendelianRandomization::mr_ivw "
       "applies by default; the classical fixed-effect interval is narrower (SE 0.0230 vs 0.1086) with "
       "the same point estimate. The IAMDGC replication behaves the same way (β 0.181–0.268, maximum "
       "P = 9.8 × 10⁻⁶). A separate leave-one-out on instrument sets that *include* rs429358 has been "
       "reported by others to lose significance for AMD when that variant is removed; that set is a "
       "different exposure definition (the ε4 coding tag is the single strongest level regulator) and is "
       "not the instrument set used here.\n\n")
t = t[:i] + seg + t[i:]
log.append("F3 新增 leave-one-out 段")

prov = (
    "> **v4 provenance (2026-09-10).** Two changes from the reproducibility audit. (i) The IAMDGC "
    "replication estimates were recomputed: the producing script passed `header = TRUE` to a `grep` "
    "extraction whose output contains no header, so the first matching variant was consumed as column "
    "names and **one instrument (rs28399649) was silently dropped** (11 → 10). With the fix, GSMR "
    "β = +0.212 (SE 0.035, P = 8.7 × 10⁻¹⁰) and IVW β = +0.213 (P = 8.5 × 10⁻⁸). (ii) A leave-one-out "
    "leverage check on the official 11-instrument set was added; the AMD estimate is not driven by any "
    "single instrument. **No other estimate, P value, confidence interval or sample size was changed.** "
    "Audit records: `docs/repro_audit_20260910/FINDINGS_MASTER_20260910.md`.\n\n"
)
t = t.replace("> **v3 provenance (2026-09-10).**", prov + "> **v3 provenance (2026-09-10).**", 1)
log.append("provenance 追加 v4")

assert t != orig, "无改动"
io.open(DST, "w", encoding="utf-8").write(t)
print("已写:", DST)
for l in log:
    print("  ✓", l)
print(f"行数 {orig.count(chr(10))+1} → {t.count(chr(10))+1}")
