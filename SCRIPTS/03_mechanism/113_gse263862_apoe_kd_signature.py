#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""113_gse263862_apoe_kd_signature.py — 人类 APOE 敲低转录组的基因水平签名（监督机制 ④）

数据: GSE263862（Homo sapiens；人星形胶质细胞）
  - siAPOE: 5 个独立 siRNA（2,3,4,5,pool），day 2
  - siAPOE NTC: 非靶向对照 1 个
  - 另有 siAANCR 系列（不同靶点）→ 用作**负对照**（特异性检验）
  数据形式: bigWig 覆盖度（.plus / .minus 各一），hg38
注释: UCSC refGene hg38（基因体 min(start)–max(end)）

产出:
  revision/results/external_validation/GSE263862_gene_coverage.csv   （基因 × 样本，CPM）
  revision/results/external_validation/GSE263862_APOE_KD_signature.csv（基因, log2FC, t, rank）
"""
import gzip, glob, os, re, json
import numpy as np, pandas as pd, pyBigWig

PROJ = "/home/huan/projects/08_APOE_AMD_AD_MR"
D = os.path.join(PROJ, "data/validation/GSE263862")
OUTD = os.path.join(PROJ, "revision/results/external_validation")
os.makedirs(OUTD, exist_ok=True)

# ---------------- 1. 基因体（hg38）
genes = {}
with gzip.open(os.path.join(D, "annot/refGene_hg38.txt.gz"), "rt") as f:
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) < 13: continue
        chrom, strand, txStart, txEnd, sym = p[2], p[3], int(p[4]), int(p[5]), p[12]
        if chrom not in genes.get(sym, ("", 0, 1 << 60, 0))[0:1] or sym not in genes:
            genes[sym] = [chrom, txStart, txEnd, strand]
        else:
            g = genes[sym]
            g[1] = min(g[1], txStart); g[2] = max(g[2], txEnd)
genes = {k: v for k, v in genes.items() if not k.startswith("NR_")}
print(f"[1] genes: {len(genes)}  (hg38)")

# ---------------- 2. 每个 bigWig 的基因覆盖度
bws = sorted(glob.glob(os.path.join(D, "bw", "*.bw")))
samples = {}
for b in bws:
    # 显式按 GSM 编号映射（两个对照的标题都是 "NTC"，按名字会撞车）
    GSM2NAME = {
        "GSM8200301": "siAANCR-2", "GSM8200302": "siAANCR-3", "GSM8200303": "siAANCR+siAP1",
        "GSM8200304": "siAANCR-pool", "GSM8200305": "AANCR-ASO-Gapmer5", "GSM8200306": "siAANCR-NTC",
        "GSM8200307": "siAPOE-2", "GSM8200308": "siAPOE-3", "GSM8200309": "siAPOE-4",
        "GSM8200310": "siAPOE-5", "GSM8200311": "siAPOE-pool", "GSM8200312": "siAPOE-NTC"}
    gsm = re.search(r"(GSM\d+)", os.path.basename(b))
    if not gsm or gsm.group(1) not in GSM2NAME: continue
    samples.setdefault(GSM2NAME[gsm.group(1)], {})["plus" if ".plus." in b else "minus"] = b
print(f"[2] samples: {len(samples)}")

genes_df = pd.DataFrame([(g, v[0], v[1], v[2], v[3]) for g, v in genes.items()],
                        columns=["symbol", "chrom", "start", "end", "strand"])
cov = {}
for s, files in samples.items():
    tot = np.zeros(len(genes_df))
    for st, path in files.items():
        bw = pyBigWig.open(path)
        chroms = bw.chroms()
        vals = np.zeros(len(genes_df))
        for i, (c, st_, en) in enumerate(zip(genes_df.chrom, genes_df.start, genes_df.end)):
            if c not in chroms: continue
            en = min(en, chroms[c])
            try:
                v = bw.stats(c, int(st_), int(en), type="sum")[0]
                if v: vals[i] = v
            except Exception:
                pass
        bw.close()
        tot += vals
    cov[s] = tot
    print(f"    {s}: total coverage {tot.sum():.3e}")

X = pd.DataFrame(cov, index=genes_df.symbol.values)
X = X.groupby(level=0).sum()                      # 同名基因合并
Xcpm = X / X.sum(axis=0) * 1e6
Xcpm.to_csv(os.path.join(OUTD, "GSE263862_gene_coverage.csv"))
print(f"[3] coverage matrix: {Xcpm.shape}")

# ---------------- 3. siAPOE 敲低签名
kd = [c for c in Xcpm.columns if re.match(r"siAPOE-\d+|siAPOE-pool", c)]
ctrl = [c for c in Xcpm.columns if c == "siAPOE-NTC"]
neg = [c for c in Xcpm.columns if re.match(r"siAANCR-\d+|siAANCR-pool", c)]
negc = [c for c in Xcpm.columns if c == "siAANCR-NTC"]
print(f"[4] siAPOE={kd}  NTC={ctrl}  负对照 siAANCR={neg}  NTC={negc}")

def signature(trt, ctl):
    l2 = np.log2(Xcpm[trt] + 1).mean(axis=1) - np.log2(Xcpm[ctl] + 1).mean(axis=1)
    sd = np.log2(Xcpm[trt] + 1).std(axis=1, ddof=1)
    return l2, sd

l2_apoe, sd_apoe = signature(kd, ctrl)
expr = Xcpm[kd + ctrl].mean(axis=1)
res = pd.DataFrame({"symbol": Xcpm.index, "mean_cpm": expr.values,
                    "log2FC_siAPOE_vs_NTC": l2_apoe.values, "sd_siAPOE": sd_apoe.values})
res["z"] = res.log2FC_siAPOE_vs_NTC / res.sd_siAPOE.replace(0, np.nan)
res = res.sort_values("log2FC_siAPOE_vs_NTC")
res.to_csv(os.path.join(OUTD, "GSE263862_APOE_KD_signature.csv"), index=False)
apoe = res[res.symbol == "APOE"]
print("[5] APOE row:", apoe.to_dict("records"))

# 负对照签名（siAANCR）
if neg and negc:
    l2_n, _ = signature(neg, negc)
    pd.DataFrame({"symbol": Xcpm.index, "log2FC_siAANCR_vs_NTC": l2_n.values}) \
      .to_csv(os.path.join(OUTD, "GSE263862_AANCR_KD_signature.csv"), index=False)
    print("[6] negative-control signature written")
print("DONE")
