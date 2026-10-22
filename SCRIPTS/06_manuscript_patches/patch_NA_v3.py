#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""patch_NA_v3.py — 把 v65 的 p 值机制表述修正同步到 NA 投稿版（2026-09-10）

v2 撤销了"多信号共定位"层，但继承了我在 v61 里写错的一条机制解释：
"从 β/SE 重算 p 从而保留这些变异"。实测证明重算 p 救不回这些变异
（|z| = 41.8–73.1；双精度下 |z| ≳ 37.5 时 2*pnorm(-|z|) 即下溢为 0，重算后仍为 0），
决定性因素是 `p > 0` 过滤器本身；且 coloc.abf 在同时提供 β 与 varbeta 时走
approx.bf.estimates，p 列不参与似然。另修正 "near-perfect-LD neighbours" 措辞
（实测 r²(rs429358, rs405509) = 0.084、rs440446 = 0.092；rs769449 = 0.77）。

输入 : NA_package/manuscript/manuscript_NA_v2.md
输出 : NA_package/manuscript/manuscript_NA_v3.md
"""
import io, os

PROJ = "/home/huan/projects/08_APOE_AMD_AD_MR"
SRC = os.path.join(PROJ, "NA_package", "manuscript", "manuscript_NA_v2.md")
DST = os.path.join(PROJ, "NA_package", "manuscript", "manuscript_NA_v3.md")

t = io.open(SRC, encoding="utf-8").read()
orig = t
log = []


def sub_once(old, new, label):
    n = t.count(old)
    assert n == 1, f"[{label}] 锚点出现 {n} 次（应为 1）: {old[:70]!r}"
    return t.replace(old, new)


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
log.append("F1 Methods：机制表述纠正")

t = sub_once(
    "comprising the ε4 coding tag rs429358 and its near-perfect-LD neighbors (rs405509, rs440446, "
    "rs769449, rs10414043 — the tag's own PIP is ≈ 0 in the R13 AD run)",
    "comprising the ε4 coding tag rs429358 together with rs405509, rs440446, rs769449 and rs10414043 "
    "(the tag's own PIP is ≈ 0 in the R13 AD run). The LD between the tag and these four variants in the "
    "1000 Genomes European reference is heterogeneous — r² = 0.77 with rs769449 of order high LD, but only "
    "0.08–0.09 with rs405509 and rs440446 — so this credible set spans more than one LD block rather than "
    "a single near-perfect-LD cluster",
    "F3")
log.append("F3 LD 措辞与实测 r² 对齐")

t = sub_once(
    "the single-signal PP4 for AD (1.0) reflects near-perfect LD between the ε4 cluster and the level "
    "instruments and does not by itself distinguish shared from independent causality",
    "the single-signal PP4 for AD (1.0) reflects the strong LD between the ε4 cluster and the level "
    "instruments and does not by itself distinguish shared from independent causality",
    "F4")
log.append("F4 去掉 near-perfect 的过度表述")

prov = (
    "> **v3 provenance (2026-09-10).** Corrects a mechanism statement introduced in v2. v2 said the "
    "31 affected APOE-region variants could be recovered by recomputing p from β/SE; measurement shows "
    "they cannot (|z| = 41.8–73.1; a recomputed p underflows to zero as well in double precision, "
    "|z| ≳ 37.5). The determinant is the `p > 0` filter itself, and in this configuration `coloc.abf` "
    "uses β and varbeta rather than the p column at all. The text now states that variants are retained "
    "on finite β and positive SE, and that 30 (not 31) survive the LD-reference filter. The description "
    "of the AD credible set, which called rs405509/rs440446 'near-perfect-LD neighbours' of rs429358 "
    "(measured r² = 0.084–0.092; rs769449 is 0.77), now reports the heterogeneous LD. **No estimate, "
    "P value, confidence interval or sample size was changed.**\n\n"
)
t = t.replace("> **v2 provenance (2026-09-10).**", prov + "> **v2 provenance (2026-09-10).**", 1)
log.append("provenance 追加 v3 说明")

assert t != orig, "无改动"
io.open(DST, "w", encoding="utf-8").write(t)
print("已写:", DST)
for l in log:
    print("  ✓", l)
print(f"行数 {orig.count(chr(10))+1} → {t.count(chr(10))+1}")
