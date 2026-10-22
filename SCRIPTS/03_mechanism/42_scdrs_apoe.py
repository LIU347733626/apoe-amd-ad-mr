#!/usr/bin/env python
# 42: scDRS on GSE230348 (92,439 scored cells, cell_type) — AMD & AD gene sets (MAGMA Z).
# Pattern: 09 project run_scdrs_fixed.py (scdrs.score_cell, n_ctrl=200).
# Gene sets: lead SNP p < 5e-8 (defensible); weights = MAGMA TOP_Z / ZSTAT.
# Output: per-cell z + per-cell-type summary + APOE expression context.
import os, sys, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, scanpy as sc, scdrs

H5AD  = "/home/huan/database/GEO/AMD/GSE230348/processed_01_results/adata_final_with_scores.h5ad"
AMD_Z = "/home/huan/projects/02_AMD_AD_Network_Pharm/results/finn_amd_gene_level_scores.csv"
AD_Z  = "/home/huan/projects/08_APOE_AMD_AD_MR/revision/results/magma_ad.genes.out"
OUTD  = "/home/huan/projects/08_APOE_AMD_AD_MR/revision/results/scdrs_apoe"
os.makedirs(OUTD, exist_ok=True)

print("loading adata...", flush=True)
adata = sc.read_h5ad(H5AD)
print(f"cells={adata.n_obs} genes={adata.n_vars}", flush=True)
if "cell_type" not in adata.obs.columns:
    print("FATAL: no cell_type column", flush=True); sys.exit(1)
scdrs.preprocess(adata)

def make_set(path, sep, gcol, zcol, pcol):
    df = pd.read_csv(path, sep=sep)
    df = df[(df[pcol] < 5e-8) & df[zcol].notna()]
    df = df.drop_duplicates(subset=[gcol])
    return dict(zip(df[gcol].astype(str), df[zcol].astype(float)))

amd_set = make_set(AMD_Z, ",", "GENE", "TOP_Z", "TOP_PVAL")
ad_set  = make_set(AD_Z, r"\s+", "GENE", "ZSTAT", "P")
print(f"gene sets: AMD={len(amd_set)} AD={len(ad_set)} genes", flush=True)
print(f"APOE in AMD set: {amd_set.get('APOE', 'ABSENT')} | in AD set: {ad_set.get('APOE', 'ABSENT')}", flush=True)

def run(name, gset):
    genes = [g for g in gset if g in adata.var_names]
    print(f"{name}: {len(genes)}/{len(gset)} genes found in adata", flush=True)
    if len(genes) < 10:
        print(f"{name}: too few genes, skipping", flush=True); return None
    w = np.array([gset[g] for g in genes], dtype=np.float64)
    res = scdrs.score_cell(adata, gene_list=genes, gene_weight=w, n_ctrl=200)
    z = np.asarray(res["zscore"]).ravel()
    adata.obs[f"scdrs_{name}_z"] = z
    print(f"{name}: z mean={z.mean():.4f} sd={z.std():.4f}", flush=True)
    return z

zamd = run("AMD", amd_set)
zad  = run("AD", ad_set)

apoe_expr = None
if "APOE" in adata.var_names:
    X = adata[:, "APOE"].X
    apoe_expr = X.toarray().ravel() if hasattr(X, "toarray") else np.asarray(X).ravel()

out = pd.DataFrame({
    "cell_type": adata.obs["cell_type"].values,
    "z_AMD": zamd if zamd is not None else np.nan,
    "z_AD": zad if zad is not None else np.nan,
    "APOE_expr": apoe_expr,
})
out.to_csv(f"{OUTD}/scdrs_cell_scores.csv", index=False)
g = out.groupby("cell_type").agg(
    n_cells=("z_AMD", "size"),
    mean_z_AMD=("z_AMD", "mean"), sd_z_AMD=("z_AMD", "std"),
    mean_z_AD=("z_AD", "mean"), sd_z_AD=("z_AD", "std"),
    mean_APOE_expr=("APOE_expr", "mean"),
).reset_index().sort_values("mean_z_AMD", ascending=False)
g.to_csv(f"{OUTD}/scdrs_celltype_summary.csv", index=False)
print(g.to_string(index=False), flush=True)
print("DONE 42 scdrs", flush=True)
