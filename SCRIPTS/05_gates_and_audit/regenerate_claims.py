#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""regenerate_claims.py — 数值一律由结果文件程序化生成后与手稿比对（监督机制 ②，2026-09-10）

规则：手稿里的**每一个派生数值**都必须能由本脚本从结果文件重新算出并格式化，
      格式化结果必须**逐字出现**在手稿中。任何"人肉转录"都会在这里 FAIL。

用法: python3 scripts/regenerate_claims.py [手稿路径]
输出: docs/audit_20260909/REGENERATE_CLAIMS_<date>.md
"""
import io, os, re, sys, csv, json, math, datetime

PROJ = "/home/huan/projects/08_APOE_AMD_AD_MR"
MS = sys.argv[1] if len(sys.argv) > 1 else os.path.join(PROJ, "final_analysis", "00_manuscript_NC_v64_20260910.md")
OUTDIR = os.path.join(PROJ, "docs", "audit_20260909")
OUTD = os.path.join(PROJ, "revision", "results")
MD = os.path.join(OUTD, "mechanism")

_SUP = str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻")
def sig(p):
    m, e = f"{p:.1e}".split("e")
    return f"{m} × 10{str(int(e)).translate(_SUP)}"

rows = []
def claim(cid, desc, expected, source):
    rows.append((cid, desc, expected, source))

# ---------------------------------------------------------------- 1. 主锚点 & 亚型 MR
mr = {(r["label"], r["Method"]): r for r in csv.DictReader(io.open(f"{OUTD}/amd_subtype/amd_subtype_mr.csv", encoding="utf-8"))}
def ivw(label):
    r = mr[(label, "IVW")]; return float(r["Beta"]), float(r["SE"]), float(r["P"])
b_all, se_all, p_all = ivw("APOE_H7")
assert abs(b_all - 0.181777146570588) < 1e-6, "H7 锚点漂移"
claim("A1", "H7_AMD 全体 IVW β", f"β = +{b_all:.3f}", "revision/results/amd_subtype/amd_subtype_mr.csv")
claim("A2", "H7_AMD 全体 IVW P", sig(p_all), "同上")
for tag, lab, name in [("dry", "APOE_DRY", "DRY_AMD"), ("wet", "APOE_WET", "WET_AMD")]:
    b, se, p = ivw(lab)
    claim(f"S-{tag}-beta", f"{name} IVW β", f"β = +{b:.3f}", "amd_subtype_mr.csv")
    claim(f"S-{tag}-se",   f"{name} IVW SE", f"SE {se:.3f}", "amd_subtype_mr.csv")
    claim(f"S-{tag}-p",    f"{name} IVW P",  sig(p), "amd_subtype_mr.csv")
b_dry, se_dry, _ = ivw("APOE_DRY"); b_wet, se_wet, _ = ivw("APOE_WET")
d = b_dry - b_wet; sd = math.sqrt(se_dry**2 + se_wet**2)
pd_ = 2 * (1 - 0.5 * (1 + math.erf(abs(d / sd) / math.sqrt(2))))
claim("S-delta-beta", "干-湿 IVW Δβ", f"Δβ = {d:.3f}", "计算：amd_subtype_mr.csv")
claim("S-delta-se",   "干-湿 Δ SE",   f"SE {sd:.3f}", "计算同上")
claim("S-delta-p",    "干-湿 Δ P",    f"P = {pd_:.2f}", "计算同上")

# 亚型镜像等位
al = {(r["subtype"], r["SNP"]): r for r in csv.DictReader(io.open(f"{OUTD}/amd_subtype/amd_subtype_allele_effects.csv", encoding="utf-8"))}
for sub in ("DRY", "WET"):
    for snp, nm in (("rs429358", "ε4"), ("rs7412", "ε2")):
        r = al[(sub, snp)]
        claim(f"AL-{sub}-{snp}", f"{sub} {nm} per-allele OR", f"OR = {float(r['OR']):.2f}", "amd_subtype_allele_effects.csv")

# MVP 亚型
mvp = {r["cohort"] + "|" + r["method"]: r for r in csv.DictReader(io.open(f"{OUTD}/mvp_ancestry/mvp_mr_all.csv", encoding="utf-8"))}
for key, nm in (("dry_META|IVW", "MVP dry"), ("wet_META|IVW", "MVP wet")):
    r = mvp[key]
    claim(f"MVP-{key}", f"{nm} IVW β", f"β = +{float(r['beta']):.3f}", "mvp_ancestry/mvp_mr_all.csv")

# ---------------------------------------------------------------- 2. 虚拟扰动（两轮）
for ct in ("Microglia", "Muller_Glia", "RPE", "Endothelial", "Rod", "Unassigned"):
    for suf, run in (("", "primary"), ("_nNet10", "higher-depth")):
        p = f"{MD}/sctenifold_KO_APOE_{ct}{suf}.csv"
        if not os.path.exists(p):
            continue
        rows_ct = list(csv.DictReader(io.open(p, encoding="utf-8")))
        n_sig = sum(1 for r in rows_ct if float(r["p.adj"]) < 0.05)
        claim(f"VKO-{ct}-{run}-nsig", f"{ct} {run} KO FDR<0.05 个数", str(n_sig), os.path.relpath(p, PROJ))
# 秩相关（KO vs OE）声明已在 v67 随虚拟过表达（OE）删除（决策 D-113）——
# 手稿不再含该比较，故此处不再生成对应声明；旧件留档
# NC_package/supplementary/_withdrawn_TableS32_APOE_virtual_KO_OE.csv（不参与投递）。
# TableS32 现为 KO-only：NC_package/supplementary/TableS32_APOE_virtual_KO.csv
_t32_ko = f"{PROJ}/NC_package/supplementary/TableS32_APOE_virtual_KO.csv"
if os.path.exists(_t32_ko):
    t32ko = list(csv.DictReader(io.open(_t32_ko, encoding="utf-8")))
    sig_counts = sorted({int(r["KO_sig_FDR05"]) for r in t32ko})
    claim("VKO-ko-nsig-range", "KO（含 APOE）FDR<0.05 个数区间",
          f"{sig_counts[0]}–{sig_counts[-1]}", "TableS32_APOE_virtual_KO.csv")

# ---------------------------------------------------------------- 3. 与手稿比对
ms = io.open(MS, encoding="utf-8").read()
body = "\n".join(l for l in ms.split("\n") if not l.startswith("> **Version note"))
ok, bad = [], []
for cid, desc, expected, source in rows:
    (ok if expected in body else bad).append((cid, desc, expected, source))

date = datetime.date.today().strftime("%Y%m%d")
out = os.path.join(OUTDIR, f"REGENERATE_CLAIMS_{date}.md")
with io.open(out, "w", encoding="utf-8") as f:
    f.write(f"# 数值再生成门禁（{datetime.date.today().isoformat()}）\n\n")
    f.write(f"手稿: `{os.path.relpath(MS, PROJ)}` ｜ 声明数值 {len(rows)} 条 ｜ **一致 {len(ok)} / 不一致 {len(bad)}**\n\n")
    if bad:
        f.write("## ❌ 手稿中找不到（或已被改写）的再生成值\n\n| id | 含义 | 再生成值 | 来源 |\n|---|---|---|---|\n")
        for cid, desc, exp, src in bad:
            f.write(f"| {cid} | {desc} | `{exp}` | `{src}` |\n")
        f.write("\n> 处理方式：要么把该数值**按上面这串字符**写回手稿，要么修改本脚本的格式化规则（须在 DECISION_LOG 记录）。\n")
    f.write("\n## ✅ 一致\n\n| id | 含义 | 值 | 来源 |\n|---|---|---|---|\n")
    for cid, desc, exp, src in ok:
        f.write(f"| {cid} | {desc} | `{exp}` | `{src}` |\n")

print(f"claims={len(rows)} consistent={len(ok)} inconsistent={len(bad)}")
for cid, desc, exp, src in bad:
    print(f"  MISSING {cid}: expected '{exp}'  <- {desc} ({src})")
print("report:", out)
sys.exit(1 if bad else 0)
