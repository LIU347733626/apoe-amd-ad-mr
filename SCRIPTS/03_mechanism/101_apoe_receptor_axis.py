#!/usr/bin/env python3
"""101_apoe_receptor_axis.py — APOE 受体/效应轴在两病相关细胞类型中的表达（GSE230348）

动机：解释"同一分子为何在两病上分岔"——检查 APOE 的已知受体/脂质转运伙伴
(LRP1/LDLR/TREM2/AXL/MERTK/ABCA1/LRP8/SORL1/PLAUR/…) 在各视网膜细胞类型中的表达。

输出：revision/results/mechanism/apoe_receptor_axis.csv（gene × cell_type：% 表达细胞 + 平均 log-norm）
"""
import numpy as np, scipy.io as sio, scipy.sparse as sp, csv, os
D = "/home/huan/projects/08_APOE_AMD_AD_MR/revision/tmp_scpagwas/gse230348_10x"
OUT = "/home/huan/projects/08_APOE_AMD_AD_MR/revision/results/mechanism/apoe_receptor_axis.csv"
genes = [l.rstrip("\n").split("\t")[-1] for l in open(D + "/genes.tsv")]
bc = [l.strip() for l in open(D + "/barcodes.tsv")]
ct_lines = [l.rstrip("\n") for l in open(D + "/celltype.tsv")]
if ct_lines[0].split("\t")[0].lower() in ("barcode", "cell"): ct_lines = ct_lines[1:]
ct = {l.split("\t")[0]: l.split("\t")[-1] for l in ct_lines}
M = sio.mmread(D + "/matrix.mtx").tocsr()
tot = np.asarray(M.sum(axis=0)).ravel(); tot[tot == 0] = 1
Mn = M.tocsc().astype(np.float32).multiply(sp.csr_matrix(1e4 / tot)).tocsr(); Mn.data = np.log1p(Mn.data)
gi = {g: i for i, g in enumerate(genes)}
TARGETS = ["APOE", "LRP1", "LDLR", "TREM2", "AXL", "MERTK", "TYRO3", "ABCA1", "ABCG1",
           "LRP8", "LRP2", "SORL1", "PLAUR", "SCARB1", "VLDLR", "ABCB1"]
cts = sorted(set(ct.values()))
rows = []
for g in TARGETS:
    if g not in gi:
        rows.append({"gene": g, "cell_type": "NOT_IN_MATRIX", "pct_expressing": "", "mean_lognorm": "", "n_cells": ""}); continue
    row = np.asarray(Mn[gi[g]].todense()).ravel()
    for c in cts:
        m = np.array([ct.get(b) == c for b in bc])
        rows.append({"gene": g, "cell_type": c, "pct_expressing": f"{100*(row[m]>0).mean():.2f}",
                     "mean_lognorm": f"{row[m].mean():.4f}", "n_cells": int(m.sum())})
with open(OUT, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["gene", "cell_type", "pct_expressing", "mean_lognorm", "n_cells"])
    w.writeheader()
    for r in rows: w.writerow(r)
print("wrote", OUT, len(rows), "rows")
# 打印矩阵（% 表达）
print(f"\n{'gene':8s} " + " ".join(f"{c[:10]:>10s}" for c in cts))
for g in TARGETS:
    sub = [r for r in rows if r["gene"] == g]
    if not sub or sub[0]["cell_type"] == "NOT_IN_MATRIX":
        print(f"{g:8s}  (not in matrix)"); continue
    d = {r["cell_type"]: r["pct_expressing"] for r in sub}
    print(f"{g:8s} " + " ".join(f"{float(d[c]):10.1f}" for c in cts))
