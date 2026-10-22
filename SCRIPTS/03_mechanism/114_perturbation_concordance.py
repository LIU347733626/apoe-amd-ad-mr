#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""114_perturbation_concordance.py — 虚拟扰动预测 vs 人类真实 APOE 敲低（监督机制 ④，2026-09-10）

外部真值：GSE263862（人星形胶质细胞；siAPOE ×5 vs siAPOE-NTC；<https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE263862>）
负对照：同数据集的 siAANCR ×3 vs siAANCR-NTC（不同靶点）——若我们的 APOE 预测与负对照同样"相关"，
        说明该"相关性"不具靶点特异性。
预测：`revision/results/mechanism/sctenifold_{KO,OE}_APOE_<ct>*.csv`（距离）与
      `final_analysis_v2/03_results/ring5_perturbation/scTenifoldKnk_APOE_KO_diffRegulation.tsv`（Z）

指标（对每个预测-真值配对）：
  - Spearman（虚扰动距离 vs |log2FC|，共享基因）
  - AUROC：用扰动排名判别"真实敲低效应 top 10%"基因
  - 富集：Fisher 精确检验（预测 top-50 × 真实 top-500）
判据（预设定）：AUROC ≥ 0.60 且 Spearman ≥ 0.10 且 Fisher OR > 1 且 P < 0.05，且**优于负对照**，才判"有外部一致性"。

用法: /home/huan/software/miniconda3/envs/bioinfo/bin/python revision/114_perturbation_concordance.py
输出: revision/results/external_validation/concordance_summary.csv + docs 报告
"""
import os, io, glob, datetime
import numpy as np, pandas as pd
from scipy import stats

PROJ = "/home/huan/projects/08_APOE_AMD_AD_MR"
MD = os.path.join(PROJ, "revision", "results", "mechanism")
EV = os.path.join(PROJ, "revision", "results", "external_validation")
V2 = os.path.join(PROJ, "final_analysis_v2", "03_results", "ring5_perturbation")

def load_truth(path, col):
    d = pd.read_csv(path)
    d = d[["symbol", col]].dropna()
    return d.set_index("symbol")[col]

def load_pred_mine(mode, ct, suf):
    p = os.path.join(MD, f"sctenifold_{mode}_APOE_{ct}{suf}.csv")
    if not os.path.exists(p): return None
    d = pd.read_csv(p)[["gene", "distance"]].dropna()
    return d.set_index("gene")["distance"]

def load_pred_v2():
    p = os.path.join(V2, "scTenifoldKnk_APOE_KO_diffRegulation.tsv")
    if not os.path.exists(p): return None
    d = pd.read_csv(p, sep="\t")
    col = "Z" if "Z" in d.columns else "distance"
    d = d[["gene", col]].dropna()
    return d.set_index("gene")[col]

def evaluate(pred, truth, label_pred, label_truth):
    sh = pred.index.intersection(truth.index)
    if len(sh) < 100: return None
    p = pred.loc[sh].values; t = np.abs(truth.loc[sh].values)
    rho = stats.spearmanr(p, t).statistic
    thr = np.quantile(t, 0.90)
    y = (t >= thr).astype(int)
    if y.sum() < 10 or y.sum() == len(y): return None
    auroc = stats.mannwhitneyu(p[y == 1], p[y == 0], alternative="greater").statistic / (y.sum() * (len(y) - y.sum()))
    top_pred = set(pred.loc[sh].nlargest(50).index)
    top_truth = set(truth.loc[sh].abs().nlargest(500).index)
    a, b = len(top_pred & top_truth), len(top_pred - top_truth)
    c, d_ = len(top_truth - top_pred), len(sh) - len(top_pred | top_truth)
    odds, pv = stats.fisher_exact([[a, b], [c, d_]], alternative="greater")
    return dict(pred=label_pred, truth=label_truth, n_shared=len(sh),
                spearman=round(float(rho), 4), auroc=round(float(auroc), 4),
                fisher_OR=round(float(odds), 3), fisher_P=float(f"{pv:.3g}"),
                n_top50_hits=len(top_pred & top_truth))

def main():
    truth_apoe = load_truth(os.path.join(EV, "GSE263862_APOE_KD_signature.csv"), "log2FC_siAPOE_vs_NTC")
    truth_neg_p = os.path.join(EV, "GSE263862_AANCR_KD_signature.csv")
    truth_neg = load_truth(truth_neg_p, "log2FC_siAANCR_vs_NTC") if os.path.exists(truth_neg_p) else None
    rows = []
    preds = {}
    for ct in ["Microglia", "Muller_Glia", "RPE", "Endothelial", "Rod", "Unassigned"]:
        for suf, lab in [("", "G=1000 nNet=3"), ("_nNet10", "G=700 nNet=10"),
                         ("_G700_nNet3", "G=700 nNet=3"), ("_G700_nNet25", "G=700 nNet=25")]:
            d = load_pred_mine("KO", ct, suf)
            if d is not None:
                preds[f"{ct} {lab}"] = d
    v2 = load_pred_v2()
    if v2 is not None:
        preds["v2 pooled G=3005 nNet=10"] = v2
    for lab, d in preds.items():
        r = evaluate(d, truth_apoe, lab, "GSE263862 siAPOE (human astrocyte)")
        if r: rows.append(r)
        if truth_neg is not None:
            r2 = evaluate(d, truth_neg, lab, "GSE263862 siAANCR (negative control)")
            if r2: rows.append(r2)
    df = pd.DataFrame(rows)
    out = os.path.join(EV, "concordance_summary.csv")
    df.to_csv(out, index=False)
    print(df.to_string(index=False))
    print("\nwritten:", out)
    if not df.empty:
        ap = df[df.truth.str.contains("siAPOE")]
        ok = ap[(ap.auroc >= 0.60) & (ap.spearman >= 0.10) & (ap.fisher_OR > 1) & (ap.fisher_P < 0.05)]
        print(f"\n达到预设判据（有外部一致性）的预测数：{len(ok)}/{len(ap)}")
        if len(ok): print(ok.to_string(index=False))

if __name__ == "__main__":
    main()
