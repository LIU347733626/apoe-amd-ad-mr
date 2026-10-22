#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""patch_NA_v2.py — 把 v61 的"撤销多信号共定位层"同步到 NA 投稿版（2026-09-10）

背景：`NA_package/manuscript/manuscript_NA_v1.md` 是 2026-09-09 从工作稿 v47 派生的
      独立投稿版（Nature Aging framing），其后未跟进工作稿 v48–v61 的修改，因此仍含
      已被证伪的"多信号（suSiE）共定位"声明。

本脚本只做一件事：把 v61 对该层的撤销**逐字同步**到 NA 版。
  - 依据：docs/repro_audit_20260910/FINDINGS_L0_L1_20260910.md（D-1）
          docs/repro_audit_20260910/AD_COLOC_3LAYER_RERUN_20260910.md（真 suSiE 退化）
  - 图编号按 NA 版风格（`Figure S5` 而非 `Supplementary Fig. 5`）。
  - **不动** NA 版中语义正确的 "three tiers"（Figure 1 图注里的**队列三档**）
    与 "isoform-specific"（研究设计 / 中介路径语境）。
  - 任何锚点未命中或偏离 → 断言失败，不写文件。

输入 : NA_package/manuscript/manuscript_NA_v1.md
输出 : NA_package/manuscript/manuscript_NA_v2.md
"""
import io, os

PROJ = "/home/huan/projects/08_APOE_AMD_AD_MR"
SRC = os.path.join(PROJ, "NA_package", "manuscript", "manuscript_NA_v1.md")
DST = os.path.join(PROJ, "NA_package", "manuscript", "manuscript_NA_v2.md")

t = io.open(SRC, encoding="utf-8").read()
orig = t
log = []


def sub_once(old, new, label):
    n = t.count(old)
    assert n == 1, f"[{label}] 锚点出现 {n} 次（应为 1）: {old[:70]!r}"
    return t.replace(old, new)


def cut_between(a, b, repl, label):
    i = t.find(a)
    assert i >= 0, f"[{label}] 起始锚点未找到: {a[:70]!r}"
    j = t.find(b, i)
    assert j >= 0, f"[{label}] 结束锚点未找到: {b[:70]!r}"
    return t[:i] + repl + t[j + len(b):]


# --- 与 v61 同样的 12 处，文本按 NA 版图编号风格适配 -----------------------
t = sub_once("**Colocalization (three tiers).**", "**Colocalization (two tiers).**", "R1")
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
t = sub_once(
    "the distinct-variant attribution (PP3 = 1.0) states that the two signals sit on different variants "
    "within the locus, not that the AD mechanism is non-level. The classic single-signal PP4 for AD (1.0) "
    "reflects near-perfect LD between the ε4 cluster and the level instruments and does not by itself "
    "distinguish shared from independent causality; multi-signal colocalization and fine-mapping do.",
    "the single-signal PP4 for AD (1.0) reflects near-perfect LD between the ε4 cluster and the level "
    "instruments and does not by itself distinguish shared from independent causality; fine-mapping, "
    "which localises the two credible sets to different variants within the block, does.",
    "R3")
t = sub_once(
    "The AD signal is therefore concentrated on the ε4 coding cluster, supporting a predominantly "
    "isoform-specific architecture in AD with a superimposed level effect detectable by multi-instrument GSMR.",
    "The AD signal is therefore concentrated on the ε4 coding cluster. Because the ε4 allele is itself the "
    "strongest regulator of APOE abundance (β = −1.012 per allele), this localisation is consistent with a "
    "predominantly level-mediated AD mechanism rather than evidence against it; a superimposed "
    "isoform-specific component of ε4 cannot be separated from the level effect with these data.",
    "R4")
t = sub_once(
    "the latter is expected if the AD signal and the level instruments do not share a single causal variant, "
    "and is consistent with the tiered colocalization below (see Limitations)",
    "the latter is expected if the AD signal sits on variants that the level instruments do not fully "
    "capture, and is consistent with the fine-mapping below (see Limitations)",
    "R5")
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
t = sub_once(
    "while the AD signal, despite single-signal PP4 = 1.0 at classic resolution, concentrates on variants "
    "distinct from the level instruments (multi-signal PP3 = 1.0, PP4 = 1.16 × 10⁻²⁰⁵; the independent-cohort "
    "descriptive Jansen 2019 z-profile shows the same ε4-cluster dominance). Fine-mapping resolves the "
    "within-locus attribution: the AD signal is concentrated on the ε4 coding cluster (PIP = 1.0 at "
    "45.40–45.42 Mb), while the level-regulating variants that drive AMD carry PIP ≈ 0 in AD "
    "(Figure S5). This supports a predominantly isoform-specific AD architecture with a "
    "superimposed, multi-instrument-detectable level effect.",
    "while the AD signal additionally colocalizes with the level instrument set at single-signal resolution "
    "(PP4 = 1.000 with the full set and 0.996 with the ε4-independent clumped set; the independent-cohort "
    "descriptive Jansen 2019 z-profile shows the same ε4-cluster dominance). Fine-mapping resolves the "
    "within-locus attribution: the AD signal is concentrated on the ε4 coding cluster (PIP = 1.0 at "
    "45.40–45.42 Mb), while the level-regulating variants that drive AMD carry PIP ≈ 0 in AD "
    "(Figure S5). Because ε4 is the strongest regulator of APOE abundance (β = −1.012 per allele), this "
    "pattern is consistent with a predominantly level-mediated AD mechanism; a superimposed "
    "isoform-specific component of ε4 remains possible but is not separable from the level effect here.",
    "R7")
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
t = sub_once(
    "**Colocalization.** Three tiers. (i) Classic single-signal: coloc R package[16], coloc.abf, default "
    "priors (p1 = p2 = 1 × 10⁻⁴, p12 = 1 × 10⁻⁵), with the full 1,407-variant APOE pQTL set (unclumped) and "
    "the 11-instrument clumped level set (ε4 coding tag excluded by design).",
    "**Colocalization.** Two tiers. (i) Classic single-signal: coloc R package[16], coloc.abf, default "
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
t = sub_once(
    "This tiered colocalization structure is a key mechanistic insight",
    "This two-layer colocalization structure is a key mechanistic insight",
    "R11")
t = sub_once(
    "The distinct-variant attribution for AD does NOT invalidate the GSMR estimate",
    "The two-layer colocalization structure for AD does not invalidate the GSMR estimate",
    "R12")
log.append("12 处替换（与 v61 同源，图编号按 NA 风格）")

# --- 断言：被撤数字不得再出现于正文（排除 A/B 节与 provenance）-------------
body = t[t.index("## Abstract"):]
for banned in ["PP3 = 1.0", "1.16 × 10⁻²⁰⁵", "0.88) and shared (0.12",
               "Multi-signal (suSiE-based) colocalization resolves"]:
    assert banned not in body, f"撤销后正文仍含: {banned}"
log.append("断言：被撤数字已从正文清除")

# --- provenance 更新 -------------------------------------------------------
prov = (
    "> **v2 provenance (2026-09-10).** This version syncs the NA submission manuscript with the "
    "withdrawal of the multi-signal (suSiE-based) colocalization layer made in the working manuscript "
    "v61. A reproducibility audit showed that the previously reported AD multi-signal result "
    "(`PP3 = 1.0`, `PP4 = 1.16 × 10⁻²⁰⁵`) was an artefact of a `p > 0` filter that silently removed the "
    "31 APOE-region variants whose p value is truncated at 0 in the FinnGen R13 AD file (including "
    "rs429358), and that the producing script had called `coloc.abf` (single-signal) rather than "
    "`coloc.susie`; a genuine suSiE run on clean input does not resolve independent signals in this "
    "region (all credible sets collapsed to single variants; LD matrix not positive semidefinite; "
    "summary-statistic/LD agreement s_rss = 0.68–0.71; two variants in perfect LD). Every claim that "
    "rested on that layer has been removed, including the inference of a predominantly isoform-specific "
    "AD architecture; because ε4 is itself the strongest APOE level regulator, the ε4-cluster "
    "localisation of the AD signal is now read as consistent with a level-mediated mechanism. The "
    "single-signal colocalization results and the fine-mapping results are unchanged. **No estimate, "
    "P value, confidence interval or sample size was changed.** Audit records: "
    "`docs/repro_audit_20260910/FINDINGS_L0_L1_20260910.md`, "
    "`docs/repro_audit_20260910/AD_COLOC_3LAYER_RERUN_20260910.md`. Note that `manuscript_NA_v1.md` "
    "derives from working-manuscript v47 and therefore predates working versions v48–v60; a full "
    "re-derivation from the current working manuscript is a separate decision and has **not** been "
    "performed here.\n\n"
)
t = t.replace("> **Provenance.**", prov + "> **Provenance.**", 1)
log.append("provenance 更新（含'仍基于 v47、未全面重建'的如实声明）")

assert t != orig, "无改动"
io.open(DST, "w", encoding="utf-8").write(t)
print("已写:", DST)
for l in log:
    print("  ✓", l)
print(f"行数 {orig.count(chr(10))+1} → {t.count(chr(10))+1}")
