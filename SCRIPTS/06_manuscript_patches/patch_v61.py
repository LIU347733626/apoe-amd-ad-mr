#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""patch_v61.py — 撤销不可复现的"多信号共定位"层（2026-09-10）

依据（可核对）：
  docs/repro_audit_20260910/FINDINGS_L0_L1_20260910.md  —— D-1：PP3=1.0 由 p>0 过滤删除
      AD 区域 31 个 p 钳制变异（含 rs429358）造成，且 55 脚本实为 coloc.abf 非 coloc.susie。
  docs/repro_audit_20260910/AD_COLOC_3LAYER_RERUN_20260910.md —— 真 coloc.susie 跑通但退化：
      credible set 全为 singleton、LD 矩阵非半正定、s_rss = 0.68/0.71、r²(rs7412,rs1065853)=1.0
      → 本数据不支持多信号共定位，第三层应整层撤销，而非改数字。

本次改动原则：
  1) 只删除/改写"多信号（suSiE）共定位"相关声明与其推论，**不新增任何未经本次验证的数字**；
  2) 保留 single-signal coloc（L1 PP4=1.000；L2 clumped PP4=0.9961）与 fine-mapping 的既有数字，
     两者均由重算-p 的干净实现产出，且本次独立复算确认；
  3) 保留全部既有多效性 caveat（MR-Egger 截距、GSMR HEIDI），它们不依赖被撤销的那一层；
  4) 任何锚点未命中 → 断言失败，不写文件。

输入 : final_analysis/00_manuscript_NC_v60_20260910.md
输出 : final_analysis/00_manuscript_NC_v61_20260910.md
"""
import io, os, sys, re

PROJ = "/home/huan/projects/08_APOE_AMD_AD_MR"
SRC = os.path.join(PROJ, "final_analysis", "00_manuscript_NC_v60_20260910.md")
DST = os.path.join(PROJ, "final_analysis", "00_manuscript_NC_v61_20260910.md")

t = io.open(SRC, encoding="utf-8").read()
orig = t
log = []


def sub_once(old, new, label):
    """全局恰好出现一次才替换"""
    n = t.count(old)
    assert n == 1, f"[{label}] 锚点出现 {n} 次（应为 1）: {old[:80]!r}"
    return t.replace(old, new)


def cut_between(a, b, repl, label):
    """删除从 a 到 b（含）之间的片段，替换为 repl"""
    i = t.find(a)
    assert i >= 0, f"[{label}] 起始锚点未找到: {a[:80]!r}"
    j = t.find(b, i)
    assert j >= 0, f"[{label}] 结束锚点未找到: {b[:80]!r}"
    return t[:i] + repl + t[j + len(b):]


# --------------------------------------------------------------------------
# R1  Results：层级标题 三 → 二
# --------------------------------------------------------------------------
t = sub_once("**Colocalization (three tiers).**", "**Colocalization (two tiers).**", "R1")
log.append("R1 Results 层级 三→二")

# --------------------------------------------------------------------------
# R2  Results：整段撤销"多信号共定位"，替换为不可行说明 + 以两层为准
# --------------------------------------------------------------------------
t = cut_between(
    "Multi-signal (suSiE-based) colocalization resolves the within-locus attribution for AD",
    "single-signal + fine-mapping chain. ",
    ("Multi-signal (suSiE-based) colocalization was **not** used to support any claim: a suSiE-based "
     "analysis of the same contrast does not resolve independent signals in this region (every credible "
     "set collapsed to a single variant, the LD matrix derived from the 1000 Genomes European reference "
     "was not positive semidefinite, and the consistency between the summary statistics and that LD "
     "reference was only s_rss = 0.68 for the pQTL and 0.71 for the AD GWAS, with r² = 1.0 between two "
     "of the retained variants) — expected, because the AD outcome is measured in a Finnish founder "
     "population and the exposure in UK Biobank, so neither shares the LD structure of the reference "
     "panel. Within-locus attribution is therefore based on the single-signal colocalization above "
     "together with fine-mapping (below), which is what the credible-set panel in **Figure 3D** shows. "),
    "R2")
log.append("R2 撤销 Results 多信号段")

# --------------------------------------------------------------------------
# R3  Results：L84 段末，去掉对多信号的依赖
# --------------------------------------------------------------------------
t = sub_once(
    "the distinct-variant attribution (PP3 = 1.0) states that the two signals sit on different variants "
    "within the locus, not that the AD mechanism is non-level. The classic single-signal PP4 for AD (1.0) "
    "reflects near-perfect LD between the ε4 cluster and the level instruments and does not by itself "
    "distinguish shared from independent causality; multi-signal colocalization and fine-mapping do.",
    "the single-signal PP4 for AD (1.0) reflects near-perfect LD between the ε4 cluster and the level "
    "instruments and does not by itself distinguish shared from independent causality; fine-mapping, "
    "which localises the two credible sets to different variants within the block, does.",
    "R3")
log.append("R3 去 PP3 依赖（Results 解读段）")

# --------------------------------------------------------------------------
# R4  Results：fine-mapping 段末——把 ε4 簇解读方向纠正为"与水平模型一致"
# --------------------------------------------------------------------------
t = sub_once(
    "The AD signal is therefore concentrated on the ε4 coding cluster, supporting a predominantly "
    "isoform-specific architecture in AD with a superimposed level effect detectable by multi-instrument GSMR.",
    "The AD signal is therefore concentrated on the ε4 coding cluster. Because the ε4 allele is itself the "
    "strongest regulator of APOE abundance (β = −1.012 per allele), this localisation is consistent with a "
    "predominantly level-mediated AD mechanism rather than evidence against it; a superimposed "
    "isoform-specific component of ε4 cannot be separated from the level effect with these data.",
    "R4")
log.append("R4 纠正 ε4 簇的解读方向")

# --------------------------------------------------------------------------
# R5  Results 首段：HEIDI 失败的归因不再引用被撤层
# --------------------------------------------------------------------------
t = sub_once(
    "the latter is expected if the AD signal and the level instruments do not share a single causal variant, "
    "and is consistent with the tiered colocalization below (see Limitations)",
    "the latter is expected if the AD signal sits on variants that the level instruments do not fully "
    "capture, and is consistent with the fine-mapping below (see Limitations)",
    "R5")
log.append("R5 首段归因改写")

# --------------------------------------------------------------------------
# R6  Discussion：去掉两处 "(multi-signal PP3 = 1.0)" 引用
# --------------------------------------------------------------------------
t = sub_once(
    "consistent with the tiered colocalization result that the AD GWAS signal sits on variants distinct "
    "from the level instrument set (multi-signal PP3 = 1.0)",
    "consistent with fine-mapping placing the AD signal on the ε4 coding cluster while the level-regulating "
    "variants carry the AMD signal",
    "R6a")
t = sub_once(
    "both consistent with the tiered colocalization result that the AD signal sits on variants distinct "
    "from the level instrument set (multi-signal PP3 = 1.0)",
    "both consistent with fine-mapping placing the AD signal on the ε4 coding cluster rather than on the "
    "level-regulating variants that carry the AMD signal",
    "R6b")
log.append("R6 去 Discussion 两处 PP3 引用")

# --------------------------------------------------------------------------
# R7  Discussion L184：去掉 PP4=1.16e-205 与"isoform-specific architecture"结论
# --------------------------------------------------------------------------
t = sub_once(
    "while the AD signal, despite single-signal PP4 = 1.0 at classic resolution, concentrates on variants "
    "distinct from the level instruments (multi-signal PP3 = 1.0, PP4 = 1.16 × 10⁻²⁰⁵; the independent-cohort "
    "descriptive Jansen 2019 z-profile shows the same ε4-cluster dominance). Fine-mapping resolves the "
    "within-locus attribution: the AD signal is concentrated on the ε4 coding cluster (PIP = 1.0 at "
    "45.40–45.42 Mb), while the level-regulating variants that drive AMD carry PIP ≈ 0 in AD "
    "(Supplementary Fig. 5). This supports a predominantly isoform-specific AD architecture with a "
    "superimposed, multi-instrument-detectable level effect.",
    "while the AD signal additionally colocalizes with the level instrument set at single-signal resolution "
    "(PP4 = 1.000 with the full set and 0.996 with the ε4-independent clumped set; the independent-cohort "
    "descriptive Jansen 2019 z-profile shows the same ε4-cluster dominance). Fine-mapping resolves the "
    "within-locus attribution: the AD signal is concentrated on the ε4 coding cluster (PIP = 1.0 at "
    "45.40–45.42 Mb), while the level-regulating variants that drive AMD carry PIP ≈ 0 in AD "
    "(Supplementary Fig. 5). Because ε4 is the strongest regulator of APOE abundance (β = −1.012 per "
    "allele), this pattern is consistent with a predominantly level-mediated AD mechanism; a superimposed "
    "isoform-specific component of ε4 remains possible but is not separable from the level effect here.",
    "R7")
log.append("R7 撤销 isoform-specific 结论")

# --------------------------------------------------------------------------
# R8  Limitations 第二条：重写，改为"不可分离"而非引用 PP3
# --------------------------------------------------------------------------
t = sub_once(
    "Second, the AD signal sits on variants distinct from the pQTL (level) instrument set (multi-signal "
    "PP3 = 1.0); we argue — and the fine-mapping (AD signal on the ε4 coding cluster), the robustness grid, "
    "the independent-cohort replications, and the directionally concordant ε4 dose–response support — that "
    "the AD association is a level effect dominated by the ε4 allele itself (the strongest level regulator), "
    "with a possible non-separable isoform-specific component; this remains a limitation to be addressed by "
    "larger or alternative instrumental approaches (e.g., isoform-resolved proteomics).",
    "Second, the ε4 coding variant is simultaneously the strongest regulator of APOE abundance and a "
    "protein-coding variant, so a level effect and an isoform-specific effect of ε4 on AD cannot be "
    "separated in these data. Fine-mapping places the AD signal on the ε4 coding cluster (PIP = 1.0), which "
    "is consistent with a level effect dominated by ε4 (the strongest level regulator), while the separate "
    "single-signal colocalization of the level instrument set with AD (PP4 = 1.000 with the full set and "
    "0.996 with the ε4-independent clumped set) speaks against a purely isoform-specific mechanism. We "
    "therefore report the AD arm as level-mediated with a possible non-separable isoform-specific component. "
    "This remains a limitation to be addressed by larger or alternative instrumental approaches (e.g., "
    "isoform-resolved proteomics). A multi-signal (suSiE-based) colocalization was attempted for this "
    "contrast and is not reported, because it does not resolve independent signals in this region: every "
    "credible set collapsed to a single variant, the 1000 Genomes European LD matrix was not positive "
    "semidefinite, the summary-statistic/LD consistency was only s_rss = 0.68 (pQTL) and 0.71 (AD), and two "
    "retained variants were in perfect LD (r² = 1.0); the AD outcome is measured in a Finnish founder "
    "population and the exposure in UK Biobank, so both depart from the reference LD structure.",
    "R8")
log.append("R8 Limitations 第二条重写 + 新增多信号不可行说明")

# --------------------------------------------------------------------------
# R9  Methods：Colocalization 段落 三 → 二 + 保留实施细节说明
# --------------------------------------------------------------------------
t = sub_once(
    "**Colocalization.** Three tiers. (i) Classic single-signal: coloc R package[44], coloc.abf, default "
    "priors (p1 = p2 = 1 × 10⁻⁴, p12 = 1 × 10⁻⁵), with the full 1,407-variant APOE pQTL set (unclumped) and "
    "the 11-instrument clumped level set (ε4 coding tag excluded by design).",
    "**Colocalization.** Two tiers. (i) Classic single-signal: coloc R package[44], coloc.abf, default "
    "priors (p1 = p2 = 1 × 10⁻⁴, p12 = 1 × 10⁻⁵), with the full 1,407-variant APOE pQTL set (unclumped) and "
    "the 11-instrument clumped level set (ε4 coding tag excluded by design). In all runs the disease "
    "p values were recomputed from the reported β/SE within the region rather than taken from the released "
    "p column, because that column is truncated at p = 0 for the strongest APOE-region variants (31 such "
    "variants in the FinnGen R13 AD file, including rs429358); retaining them is what the estimates above "
    "are based on.",
    "R9")
t = sub_once(
    "(ii) Multi-signal: suSiE-based colocalization (susieR 0.14.2; posterior H0–H4 with H3 = both traits "
    "present via distinct causal variants, H4 = shared). (iii) Fine-mapping: see below.",
    "(ii) Fine-mapping: see below. A multi-signal (suSiE-based) colocalization was attempted and is "
    "deliberately not reported: in this region every credible set collapsed to a single variant, the LD "
    "matrix was not positive semidefinite, and the agreement between the summary statistics and the "
    "1000 Genomes European LD reference was low (s_rss = 0.68 for the pQTL and 0.71 for the AD GWAS), so "
    "the analysis would report LD structure rather than shared or distinct causal variants.",
    "R10")
log.append("R9+ Methods Colocalization 段改两层 + 重算 p 说明 + 多信号不可行说明")

# --------------------------------------------------------------------------
# R11/R12 指代清理：撤销该层后，"tiered"/"distinct-variant attribution" 失去所指
# --------------------------------------------------------------------------
t = sub_once(
    "This tiered colocalization structure is a key mechanistic insight",
    "This two-layer colocalization structure is a key mechanistic insight",
    "R11")
t = sub_once(
    "The distinct-variant attribution for AD does NOT invalidate the GSMR estimate",
    "The two-layer colocalization structure for AD does not invalidate the GSMR estimate",
    "R12")
log.append("R11/R12 清理悬空指代")

# --------------------------------------------------------------------------
# 断言：被撤销的数字不得再出现于正文
# --------------------------------------------------------------------------
body = t[:t.index("## Figure Legends")] if "## Figure Legends" in t else t
body_only = "\n".join(l for l in body.split("\n") if not l.startswith("> **Version note"))
for banned in ["PP3 = 1.0", "1.16 × 10⁻²⁰⁵", "0.88) and shared (0.12", "Multi-signal (suSiE-based) colocalization resolves"]:
    assert banned not in body_only, f"撤销后正文仍含: {banned}"
log.append("断言：被撤数字已从正文清除")

# 版本说明
note = (
    "> **Version note (2026-09-10, v61, withdrawal of the multi-signal colocalization layer):** "
    "A reproducibility audit found that the previously reported multi-signal (suSiE-based) colocalization "
    "for AD (`PP3 = 1.0`, `PP4 = 1.16 × 10⁻²⁰⁵`) was not reproducible and should not have been reported: "
    "(i) the script that produced it called `coloc.abf` (single-signal), not `coloc.susie`, despite the "
    "label; (ii) its input filter `p > 0` silently removed the 31 APOE-region variants whose p value is "
    "truncated at 0 in the FinnGen R13 AD file, including rs429358 — recomputing p from β/SE restores them "
    "and returns PP4 ≈ 1; and (iii) a genuine `coloc.susie` run on the clean input does not resolve "
    "independent signals here (all credible sets collapsed to single variants; the LD matrix was not "
    "positive semidefinite; summary-statistic/LD agreement was s_rss = 0.68–0.71; two retained variants "
    "were in perfect LD). The multi-signal layer is therefore withdrawn in full, together with every claim "
    "that rested on it, including the inference of a predominantly isoform-specific AD architecture. The "
    "single-signal colocalization results (full set PP4 = 1.000; ε4-independent clumped set PP4 = 0.996), "
    "which were always computed from β/SE-recomputed p values, are unchanged and are now the basis of the "
    "two-tier colocalization statement. The fine-mapping results (Figure 3D) are unchanged. Because ε4 is "
    "itself the strongest APOE level regulator, the ε4-cluster localisation of the AD signal is now read as "
    "consistent with a level-mediated mechanism. **No estimate, P value, confidence interval or sample size "
    "was changed; only claims resting on the withdrawn layer were removed.** Audit records: "
    "`docs/repro_audit_20260910/FINDINGS_L0_L1_20260910.md` and "
    "`docs/repro_audit_20260910/AD_COLOC_3LAYER_RERUN_20260910.md`. v60 preserved as "
    "`00_manuscript_NC_v60_20260910.md`.\n\n"
)
t = t.replace("# APOE abundance drives inverse comorbidity", note + "# APOE abundance drives inverse comorbidity", 1)
assert t != orig, "无改动"

io.open(DST, "w", encoding="utf-8").write(t)
print("已写:", DST)
for l in log:
    print("  ✓", l)
print(f"行数 {orig.count(chr(10))+1} → {t.count(chr(10))+1}")
