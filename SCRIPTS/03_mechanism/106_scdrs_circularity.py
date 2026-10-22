#!/usr/bin/env python
"""106_scdrs_circularity.py — scDRS 循环性检验：剔除 APOE 区域基因后重跑

目的：scDRS/seismic 的基因集来自 MAGMA（含 APOE 区域基因），可能循环。
做法：把 chr19:44.0–46.0 Mb 的基因从两个基因集中剔除，重跑 scDRS，比较细胞类型排序是否保持。
输出：revision/results/scdrs_apoe/scdrs_celltype_tests_noAPOEregion.csv
"""
import os, sys, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, scanpy as sc, scdrs
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

PROJ = "/home/huan/projects/08_APOE_AMD_AD_MR"
H5AD = "/home/huan/database/GEO/AMD/GSE230348/processed_01_results/adata_final_with_scores.h5ad"
AMD_Z = "/home/huan/projects/02_AMD_AD_Network_Pharm/results/finn_amd_gene_level_scores.csv"
AD_Z  = f"{PROJ}/revision/results/magma_ad.genes.out"
REGION = f"{PROJ}/revision/results/mechanism/apoe_region_genes_grch37.csv"
OUTD  = f"{PROJ}/revision/results/scdrs_apoe"

reg = pd.read_csv(REGION)
gcol = "gene" if "gene" in reg.columns else reg.columns[0]
region_genes = set(reg[gcol].astype(str))
print(f"APOE 区域基因: {len(region_genes)}", flush=True)

adata = sc.read_h5ad(H5AD)
scdrs.preprocess(adata)

def make_set(path, sep, gcol_, zcol, pcol, drop_region=True):
    df = pd.read_csv(path, sep=sep)
    df = df[(df[pcol] < 5e-8) & df[zcol].notna()].drop_duplicates(subset=[gcol_])
    df = df[~df[gcol_].astype(str).isin(region_genes)] if drop_region else df
    return dict(zip(df[gcol_].astype(str), df[zcol].astype(float)))

for tag, drop in [("noAPOEregion", True), ("withAPOEregion", False)]:
    amd_set = make_set(AMD_Z, ",", "GENE", "TOP_Z", "TOP_PVAL", drop)
    ad_set  = make_set(AD_Z, r"\s+", "GENE", "ZSTAT", "P", drop)
    print(f"[{tag}] AMD={len(amd_set)} AD={len(ad_set)} genes | APOE in AMD: {'APOE' in amd_set}", flush=True)
    rows = []
    for name, gset in [("AMD", amd_set), ("AD", ad_set)]:
        genes = [g for g in gset if g in adata.var_names]
        w = np.array([gset[g] for g in genes], dtype=np.float64)
        res = scdrs.score_cell(adata, gene_list=genes, gene_weight=w, n_ctrl=200)
        z = np.asarray(res["zscore"]).ravel()
        ct = adata.obs["cell_type"].values
        for c in sorted(set(ct)):
            m = ct == c
            u, p = mannwhitneyu(z[m], z[~m], alternative="two-sided")
            rb = 2*u/(m.sum()*(~m).sum()) - 1
            rows.append({"set": tag, "disease": name, "cell_type": c, "n_cells": int(m.sum()),
                         "mean_z": float(z[m].mean()), "MW_P": float(p), "rank_biserial": float(rb)})
    d = pd.DataFrame(rows)
    d["FDR"] = np.nan
    for dis in ["AMD", "AD"]:
        sel = d["disease"] == dis
        d.loc[sel, "FDR"] = multipletests(d.loc[sel, "MW_P"], method="fdr_bh")[1]
    out = f"{OUTD}/scdrs_celltype_tests_{tag}.csv"
    d.to_csv(out, index=False)
    print(d.sort_values(["disease", "FDR"]).to_string(index=False), flush=True)
    print(f"wrote {out}\n", flush=True)
