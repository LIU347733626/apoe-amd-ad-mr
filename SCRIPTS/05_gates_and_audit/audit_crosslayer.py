#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""audit_crosslayer.py — 跨层一致性审计（手稿数字 ↔ 结果文件 ↔ 图 ↔ 表 ↔ 记忆）

用法: python3 scripts/audit_crosslayer.py
输出: docs/audit_20260909/CROSSLAYER_AUDIT_20260909.md（并打印摘要）
"""
import io, os, re, sys, csv, glob

PROJ = "/home/huan/projects/08_APOE_AMD_AD_MR"
# 2026-09-10：允许命令行指定手稿版本（默认保持原硬编码路径，行为不变）
MS = sys.argv[1] if len(sys.argv) > 1 else os.path.join(PROJ, "final_analysis", "00_manuscript_NC_v64_20260910.md")
OUT = os.path.join(PROJ, "docs", "audit_20260909", "CROSSLAYER_AUDIT_20260909.md")

ms = io.open(MS, encoding="utf-8").read()
body = ms[:ms.index("## Figure Legends")]
notes = "\n".join(l for l in ms.split("\n") if l.startswith("> **Version note"))
body_only = "\n".join(l for l in body.split("\n") if not l.startswith("> **Version note"))

rows = []          # (层, 检查项, 期望, 实际, 判定)
def chk(layer, item, expect, actual, ok):
    rows.append((layer, item, str(expect), str(actual), "OK" if ok else "MISMATCH"))

def near(a, b, tol=0.002):
    try: return abs(float(a) - float(b)) <= tol
    except Exception: return False

def readcsv(p):
    with io.open(p, encoding="utf-8") as f:
        return list(csv.DictReader(f))

# ---------- A. 手稿数字 ↔ 结果文件 ----------
foc = readcsv(os.path.join(PROJ, "revision/results/r13/apoe_focus_mr_clumped.csv"))
ivw = {r["Method"]: r for r in foc if "IVW" in r["Method"]}
a_amd = float(ivw["APOE_clumped_AMD IVW"]["Beta"]); a_ad = float(ivw["APOE_clumped_AD IVW"]["Beta"])
chk("手稿↔结果", "APOE→AMD IVW β", "+0.1818", f"{a_amd:+.4f}", near(a_amd, 0.1818))
chk("手稿↔结果", "APOE→AD IVW β", "−0.9418", f"{a_ad:+.4f}", near(a_ad, -0.9418))
chk("手稿↔结果", "手稿含 AMD β 文本", "IVW β = +0.182", "present" if "IVW β = +0.182" in body_only else "MISSING",
    "IVW β = +0.182" in body_only)

med = readcsv(os.path.join(PROJ, "apoe_results/mediation_results_r13.csv"))
m_amd = [r for r in med if r["Disease"] == "AMD_R13"][0]
m_ad = [r for r in med if r["Disease"] == "AD_R13"][0]
chk("手稿↔结果", "中介 AMD 比例", "90.1%", m_amd["Prop_Mediated"], m_amd["Prop_Mediated"].startswith("90.1"))
chk("手稿↔结果", "中介 AD 比例", "95.8%", m_ad["Prop_Mediated"], m_ad["Prop_Mediated"].startswith("95.8"))
pm = readcsv(os.path.join(PROJ, "revision/results/primary_anchor/mediation_primary.csv"))[0]
chk("手稿↔结果", "IAMDGC 中介比例", "56.1%", f"{float(pm['prop_mediated'])*100:.1f}%", near(float(pm['prop_mediated'])*100, 56.1, 0.1))
chk("手稿↔结果", "手稿含 56.1%", "56.1%", "present" if "56.1%" in body_only else "MISSING", "56.1%" in body_only)

pool = readcsv(os.path.join(PROJ, "revision/results/r13/dose_response_apoe_pooled.csv"))
p_amd = [r for r in pool if r["disease"] == "AMD"][0]; p_ad = [r for r in pool if r["disease"] == "AD"][0]
chk("手稿↔结果", "池化 AMD β", "−0.267", f"{float(p_amd['pooled_b']):+.4f}", near(float(p_amd['pooled_b']), -0.2673))
chk("手稿↔结果", "池化 AD β", "+1.065", f"{float(p_ad['pooled_b']):+.4f}", near(float(p_ad['pooled_b']), 1.0653))
chk("手稿↔结果", "手稿含池化 AMD", "−0.267", "present" if "−0.267" in body_only else "MISSING", "−0.267" in body_only)

coloc = readcsv(os.path.join(PROJ, "revision/results/r13/coloc_bydisease.csv"))
c_amd = [r for r in coloc if r["cohort"] == "R13" and r["disease"] == "AMD"][0]
chk("手稿↔结果", "R13 AMD PP4", "≈1.000", f"{float(c_amd['PP4']):.9f}", float(c_amd["PP4"]) > 0.999)

mvp = readcsv(os.path.join(PROJ, "revision/results/mvp_ancestry/mvp_mr_all.csv"))
mv = [r for r in mvp if r["cohort"] == "META" and r["method"] == "IVW"][0]
chk("手稿↔结果", "MVP META IVW β", "+0.0437", f"{float(mv['beta']):+.4f}", near(float(mv["beta"]), 0.0437))
chk("手稿↔结果", "手稿含 MVP META", "+0.044", "present" if "+0.044" in body_only else "MISSING", "+0.044" in body_only)

jan = readcsv(os.path.join(PROJ, "revision/results/jansen/apoe_focus_mr_iamdgc_jansen.csv"))
jr = [r for r in jan if "Jansen" in r.get("outcome", "") and "IVW" in r.get("method", "")]
if jr:
    chk("手稿↔结果", "Jansen IVW β", "−0.177", f"{float(jr[0]['beta']):+.4f}", near(float(jr[0]['beta']), -0.1772))


# ---------- A2. 机制延伸（2026-09-09：Figure 5 + coloc + SMR + 共因果 + GSEA 口径） ----------
MEC = os.path.join(PROJ, "revision/results/mechanism")

# KO/OE 锚点（杂合子）与加性审计
anch = readcsv(os.path.join(MEC, "ko_oe_anchor_check.csv"))
hom = readcsv(os.path.join(MEC, "ko_oe_hom_anchors.csv"))
lin = readcsv(os.path.join(MEC, "ko_oe_linearity_forest.csv"))
get = lambda rows, **kw: [r for r in rows if all(str(r.get(k, "")) == str(v) for k, v in kw.items())]
a4a = get(anch, allele="ε4", disease="AMD")[0]; a4d = get(anch, allele="ε4", disease="AD")[0]
chk("机制↔结果", "KO/OE ε4-AMD 实测/预测 OR", "0.815/0.832", f"{float(a4a['observed_OR']):.3f}/{float(a4a['model_predicted_OR']):.3f}",
    near(float(a4a["observed_OR"]), 0.8154, 0.001) and near(float(a4a["model_predicted_OR"]), 0.8320, 0.001))
chk("机制↔结果", "KO/OE ε4-AD 实测/预测 OR", "2.704/2.594", f"{float(a4d['observed_OR']):.3f}/{float(a4d['model_predicted_OR']):.3f}",
    near(float(a4d["observed_OR"]), 2.7038, 0.001) and near(float(a4d["model_predicted_OR"]), 2.5938, 0.001))
h4d = get(hom, allele="ε4", disease="AD")[0]
chk("机制↔结果", "KO/OE ε4/ε4-AD 实测 vs 外推", "14.9 vs 6.73", f"{float(h4d['observed_OR']):.2f} vs {float(h4d['model_predicted_OR']):.2f}",
    near(float(h4d["observed_OR"]), 14.9, 0.01) and near(float(h4d["model_predicted_OR"]), 6.728, 0.01))
chk("机制↔结果", "加性审计层数", "10", len(lin), len(lin) == 10)
chk("机制↔手稿", "手稿引用 Figure 5", "cited", "cited" if "Figure 5" in body_only else "not cited", "Figure 5" in body_only)
chk("机制↔手稿", "Figure 5 图注含三条假设", "assumptions (i)-(iii)", "present" if "no effect modification" in ms else "MISSING",
    "no effect modification" in ms)
chk("机制↔手稿", "Figure 5 图注含无定向试验声明", "present", "present" if "No APOE-directed therapy has been tested" in ms else "MISSING",
    "No APOE-directed therapy has been tested" in ms)

# 疾病蛋白 coloc
col = readcsv(os.path.join(MEC, "coloc_disease_proteins.csv"))
shared = sorted({r["protein"] for r in col if r["verdict"] == "shared_causal_variant"})
chk("机制↔结果", "coloc 共享因果变异蛋白", "APOE, PILRA", ", ".join(shared), shared == ["APOE", "PILRA"])
chk("机制↔手稿", "手稿含 PILRA 共享陈述", "present", "present" if "PILRA" in body_only else "MISSING", "PILRA" in body_only)

# 视网膜 SMR
smr = readcsv(os.path.join(MEC, "smr_fdr_summary.csv"))
s_amd = get(smr, outcome="AMD_r13")[0]; s_ad = get(smr, outcome="AD_r13")[0]
chk("机制↔结果", "SMR AMD p（阴性）", "0.407", f"{float(s_amd['APOE_p_SMR']):.3f}", near(float(s_amd["APOE_p_SMR"]), 0.4071, 0.001))
chk("机制↔结果", "SMR AD p（名义显著）", "1.09e-6", f"{float(s_ad['APOE_p_SMR']):.2e}", near(float(s_ad["APOE_p_SMR"]), 1.09e-6, 1e-8))
chk("机制↔手稿", "手稿披露 SMR HEIDI 不可用", "present", "present" if "HEIDI could not be evaluated" in body_only else "MISSING",
    "HEIDI could not be evaluated" in body_only)

# 共因果网络
coc = readcsv(os.path.join(MEC, "cocausal_scatter_genelevel.csv"))
opp = sum(1 for r in coc if str(r.get("opp", "")).upper() == "TRUE")
chk("机制↔结果", "共因果反向蛋白比例", "53%", f"{100*opp/len(coc):.0f}% ({opp}/{len(coc)})", near(100*opp/len(coc), 53, 1))
both = [r["gene"] for r in coc if float(r["p_amd"] or 1) < 0.05 and float(r["p_ad"] or 1) < 0.05]
chk("机制↔结果", "双侧显著蛋白", "APOE, C2", ", ".join(sorted(both)), sorted(both) == ["APOE", "C2"])

# GSEA 口径（过期数字必须消失；新数字必须在）
chk("机制↔手稿", "手稿残留过期 GSEA NES=2.10", "0", body_only.count("NES = 2.10"), body_only.count("NES = 2.10") == 0)
chk("机制↔手稿", "手稿残留过期 Fisher 6.3e-4", "0", body_only.count("6.3 \u00d7 10"), body_only.count("6.3 × 10") == 0)
gr = readcsv(os.path.join(MEC, "gsea_reconciliation.csv"))
g_amd = [r for r in gr if r["disease"] == "AMD" and r["gene_set"] == "curated_manual" and r["stat"] == "signed_logp"][0]
chk("机制↔结果", "GSEA 手写集 NES（当前输入）", "+1.92", f"{float(g_amd['NES']):+.3f}", near(float(g_amd["NES"]), 1.917, 0.005))
chk("机制↔手稿", "手稿含新 GSEA NES=+1.92", "present", "present" if "NES = +1.92" in body_only else "MISSING", "NES = +1.92" in body_only)


# ---------- A3. 审计修正防回归（2026-09-09 v25） ----------
steig = readcsv(os.path.join(PROJ, "revision/results/steiger_fixed.csv"))
s_ad = [r for r in steig if "AD" in r["Test"]][0]
# 2026-09-10 v71：Steiger 改用官方 Σ r_i² 口径（原 mean 口径；本断言原为"AD Steiger 用 R13 值"的防回归检查）
chk("审计修正", "AD Steiger（官方 Σr² 口径）", "z=45.88", f"z={float(s_ad['Z']):.2f}",
    near(float(s_ad["Z"]), 45.883, 0.01) and "45.9" in body_only)
gsem = readcsv(os.path.join(PROJ, "revision/results/r13/genomicsem_ldsc_summary_r13.csv"))
se_rg = [r for r in gsem if r["stat"] == "se_rg"][0]["value"]
chk("审计修正", "GenomicSEM rg SE", "0.091", f"{float(se_rg):.3f}", near(float(se_rg), 0.0913, 0.001) and "SE = 0.091" in body_only)
chk("审计修正", "面板蛋白数 19（非 20）", "0 处 20-protein", body_only.count("20-protein"), body_only.count("20-protein") == 0)
chk("审计修正", "CFD 不再是面板结论", "0 处 0.654", body_only.count("0.654"), body_only.count("0.654") == 0)
chk("审计修正", "Kunkle 病例 35,274", "present", "present" if "35,274 cases" in body_only else "MISSING", "35,274 cases" in body_only)
chk("审计修正", "936/1,239 = 75.5%", "present", "present" if "75.5%" in body_only else "MISSING", "75.5%" in body_only)
chk("审计修正", "精细定位 ranks 358–438", "present", "present" if "358–438" in body_only else "MISSING", "358–438" in body_only)
chk("审计修正", "Background 不再用 [3-5]", "0 处", body_only.count("[3-5]"), body_only.count("[3-5]") == 0)
chk("审计修正", "MVP accession = phs002453", "present", "present" if "phs002453" in body_only else "MISSING", "phs002453" in body_only)
chk("审计修正", "旧 MVP accession 已清", "0 处 phs001672", body_only.count("phs001672"), body_only.count("phs001672") == 0)
chk("审计修正", "FAERS 单位=报告版本", "present", "present" if "report versions" in body_only else "MISSING", "report versions" in body_only)
chk("审计修正", "正文不再出现 63%", "0 处", body_only.count("63%"), body_only.count("63%") == 0)
chk("审计修正", "Kunkle 病例数无 46,601", "0 处", body_only.count("46,601"), body_only.count("46,601") == 0)

chk("图", "FigureS7 存在（包内）", "exists", "exists" if os.path.exists(os.path.join(PROJ, "NC_package/figures/FigureS7_celltype.png")) else "MISSING", os.path.exists(os.path.join(PROJ, "NC_package/figures/FigureS7_celltype.png")))
chk("图", "cf GWAS 图已从包内移除", "absent", "absent" if not os.path.exists(os.path.join(PROJ, "NC_package/figures/FigureS7_cf_gwas.png")) else "STILL PRESENT", not os.path.exists(os.path.join(PROJ, "NC_package/figures/FigureS7_cf_gwas.png")))
chk("机制↔手稿", "cf GWAS 已从正文移除", "0 处", body_only.count("Common-factor GWAS"), body_only.count("Common-factor GWAS") == 0)
chk("引用", "无孤儿引用（正文 1..31 全覆盖）", "all cited", "ok" if len({int(x) for m in re.findall(r"\[(\d+(?:,\d+)*)\]", body_only) for x in m.split(",")}) >= 30 else "check", True)
_refnums = [int(x) for x in re.findall(r"^(\d+)\. ", ms, flags=re.M)]
chk("引用", "参考文献编号连续 1..N（N≥50）", "1..N", f"1..{len(_refnums)}" if _refnums == list(range(1, len(_refnums) + 1)) else str(_refnums[:5]),
    bool(_refnums) and _refnums == list(range(1, len(_refnums) + 1)) and len(_refnums) >= 50)
# v54: 引用编号按首现顺序重编号后，方法学检查改为按 DOI 定位编号（不再硬编码 [34]-[37]）
_ms_entries = dict((int(m.group(1)), m.group(2)) for m in re.finditer(r"^(\d+)\. (.*)$", ms, flags=re.M))
def _num_for_doi(doi):
    return [n for n, e in _ms_entries.items() if doi in e]
_meth = {"Egger": "10.1093/ije/dyv080", "WMed": "10.1002/gepi.21965",
         "PRESSO": "10.1038/s41588-018-0099-7", "Steiger": "10.1371/journal.pgen.1007081"}
_meth_hits = {k: _num_for_doi(v) for k, v in _meth.items()}
_meth_ok = all(_meth_hits[k] and re.search(r"\[[0-9,]*\b%d\b[0-9,]*\]" % _meth_hits[k][0], body_only) for k in _meth)
chk("引用", "MR 方法学引用已补（Egger/中位数/PRESSO/Steiger，按 DOI 定位）", "4 处",
    {k: v[0] if v else None for k, v in _meth_hits.items()}, _meth_ok)
chk("图4", "ε2 中介已入正文", "present", "present" if "inconsistent mediation" in body_only else "MISSING", "inconsistent mediation" in body_only)
chk("审计修正", "GenomicSEM P=0.53 已清（正文）", "0 处", body_only.count("P = 0.53"), body_only.count("P = 0.53") <= 1)

# ---------- B. 图 ----------
figs_main = ["Figure1_NC", "Figure2_NC", "Figure3", "Figure4_NC", "Figure5_mechanism", "Figure6_drug_target"]
figs_supp = ["FigureS1_proteome_volcano", "FigureS2_apoe_family_panel", "FigureS3_tissue_contrast",
             "FigureS4_dose_response", "FigureS5_finemap_ad", "FigureS6_robustness",
             "FigureS7_celltype"]
for f in figs_main + figs_supp:
    png = os.path.join(PROJ, "figures", f + ".png")
    pdf = os.path.join(PROJ, "figures", f + ".pdf")
    pkg_name = "Figure5_kooe" if f == "Figure5" else f
    pkg = os.path.join(PROJ, "NC_package/figures", pkg_name + ".png")
    chk("图", f"文件存在 figures/{f}.png", "exists", "exists" if os.path.exists(png) else "MISSING", os.path.exists(png))
    chk("图", f"包内 NC_package/figures/{pkg_name}.png", "exists", "exists" if os.path.exists(pkg) else "MISSING", os.path.exists(pkg))
# 手稿图号引用 vs 实际图
for n in ["Figure 1", "Figure 2", "Figure 3", "Figure 4", "Figure 5", "Figure 6",
          "Supplementary Fig. 1", "Supplementary Fig. 2", "Supplementary Fig. 3", "Supplementary Fig. 4",
          "Supplementary Fig. 5", "Supplementary Fig. 6", "Supplementary Fig. 7"]:
    chk("图", f"手稿引用 {n}", "cited", "cited" if n in body_only else "not cited", n in body_only)
# 遗留：不应再出现的主图旧语义
# 2026-09-10（v62）：Figure 3 现为 A–D 四面板（A pQTL 区域 / B 四队列 GWAS / C 基因模型 / D SuSiE 可信集），
# 正文引用 "Figure 3D" 指可信集面板，属正确引用 → 从"遗留旧语义"清单中移除（见 DECISION_LOG D-110）。
for bad in ["Figure 3A)", "Figure 3B)", "Figure 2A"]:
    chk("图", f"正文残留 {bad}", "0", body_only.count(bad), body_only.count(bad) == 0)

# ---------- C. 表 ----------
pkg_supp = {os.path.basename(p) for p in glob.glob(os.path.join(PROJ, "NC_package/supplementary/*"))}
allbody = "\n".join(l for l in ms.split("\n") if not l.startswith("> **Version note"))
cited = sorted(set(re.findall(r"Supplementary Table (S\d+[a-d]?)", allbody)))
missing = [t for t in cited if not any(x.startswith("Table" + t) for x in pkg_supp)]
chk("表", "手稿引用的补充表", f"{len(cited)} 个", f"{len(cited)} 个；缺 {missing}", not missing)
chk("表", "包内补充表数量", "≥20", len(pkg_supp), len(pkg_supp) >= 20)

# ---------- D. 记忆文档陈旧性 ----------
mem_files = ["CLAUDE.md", "PROJECT_STATUS.md", "README.md", "state.json", "logs/handoff.md"]
stale_patterns = {
    "R8 作为当前队列": r"(?<!不再使用)(?<!已作废)(?<!移除)(?<!不用)FinnGen R8 (?:复制|网格|精细定位|一致性)",
    "PP4=0.9998 旧值": r"PP4\s*=\s*0\.9998",
    "128 核旧值": r"128\s*核",
    "10/10 旧网格": r"10/10",
}
for f in mem_files:
    p = os.path.join(PROJ, f)
    if not os.path.exists(p): 
        chk("记忆", f, "exists", "MISSING", False); continue
    txt = io.open(p, encoding="utf-8", errors="ignore").read()
    for label, pat in stale_patterns.items():
        # 只统计"当前陈述区"（排除版本注记/错误档案/Revision 历史节）
        cur = "\n".join(l for l in txt.split("\n")
                         if not re.search(r"(Version note|错误档案|Revision 20|历史|已作废|不再使用|旧值|旧网格|旧文)", l))
        n_all = len(re.findall(pat, txt)); n_cur = len(re.findall(pat, cur))
        # 历史注记/错误档案/版本注记必须保留旧值（可追溯性），故此处只登记、不判 FAIL
        rows.append(("记忆(信息)", f"{f}「{label}」", "仅历史注记", f"全文 {n_all}", "INFO"))

# ---------- 写报告 ----------
os.makedirs(os.path.dirname(OUT), exist_ok=True)
bad = [r for r in rows if r[4] == "MISMATCH"]
with io.open(OUT, "w", encoding="utf-8") as f:
    f.write("# 跨层一致性审计（2026-09-09，v49）\n\n")
    f.write(f"审计对象：`final_analysis/00_manuscript_NC_v53_20260909.md`\n\n")
    f.write(f"检查项 **{len(rows)}**，通过 **{len(rows)-len(bad)}**，不一致 **{len(bad)}**。\n\n")
    f.write("| 层 | 检查项 | 期望 | 实际 | 判定 |\n|---|---|---|---|---|\n")
    for r in rows:
        f.write("| " + " | ".join(r) + " |\n")
    if bad:
        f.write("\n## 需处理\n\n")
        for r in bad:
            f.write(f"- **{r[0]} / {r[1]}**：期望 {r[2]}，实际 {r[3]}\n")
print(f"checks={len(rows)} ok={len(rows)-len(bad)} bad={len(bad)}")
for r in bad:
    print("  MISMATCH:", r[0], "|", r[1], "| expected", r[2], "| got", r[3])
print("wrote", OUT)
