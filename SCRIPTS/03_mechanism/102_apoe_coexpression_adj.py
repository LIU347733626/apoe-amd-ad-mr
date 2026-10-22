#!/usr/bin/env python3
"""102_apoe_coexpression_adj.py — APOE 共表达程序（**去技术混杂**版）

问题：原始共表达 top 命中全是 MALAT1/MT-/核糖体基因 = 文库大小与细胞质量伪迹。
做法：偏相关——控制每细胞的 log 总 UMI 与 log 检出基因数后，再算 APOE 与各基因的相关。
     用偏协方差公式在稀疏矩阵上直接算，不稠密化。
输出：revision/results/mechanism/apoe_coexpression_adj_<celltype>.csv
"""
import numpy as np, scipy.io as sio, scipy.sparse as sp, os, sys
D="/home/huan/projects/08_APOE_AMD_AD_MR/revision/tmp_scpagwas/gse230348_10x"
OUT="/home/huan/projects/08_APOE_AMD_AD_MR/revision/results/mechanism"
genes=[l.rstrip("\n").split("\t")[-1] for l in open(D+"/genes.tsv")]
bc=[l.strip() for l in open(D+"/barcodes.tsv")]
ct_lines=[l.rstrip("\n") for l in open(D+"/celltype.tsv")]
if ct_lines[0].split("\t")[0].lower() in ("barcode","cell"): ct_lines=ct_lines[1:]
ct={l.split("\t")[0]: l.split("\t")[-1] for l in ct_lines}
M=sio.mmread(D+"/matrix.mtx").tocsr()
tot=np.asarray(M.sum(axis=0)).ravel(); tot[tot==0]=1
ng=np.asarray((M>0).sum(axis=0)).ravel()
Mn=M.tocsc().astype(np.float32).multiply(sp.csr_matrix(1e4/tot)).tocsr(); Mn.data=np.log1p(Mn.data)
gi={g:i for i,g in enumerate(genes)}
ia=gi["APOE"]
Y=Mn  # 基因×细胞
for ctn in sorted(set(ct.values())):
    keep=np.array([i for i,b in enumerate(bc) if ct.get(b)==ctn])
    if len(keep)<200: continue
    y=np.asarray(Y[ia][:,keep].todense()).ravel()
    X=np.column_stack([np.ones(len(keep)), np.log1p(tot[keep]), np.log1p(ng[keep])])
    XtX_inv=np.linalg.pinv(X.T@X)
    # 所有基因与 X 的叉积
    XtY=np.asarray((Y[:,keep]@X))          # 基因×3
    YtY=np.asarray(Y[:,keep].multiply(Y[:,keep]).sum(axis=1)).ravel()
    Xty=X.T@y
    yty=y@y
    n=len(keep)
    # 偏协方差：cov(g,y|X) = cov(g,y) - cov(g,X) Sigma_X^-1 cov(X,y)
    cov_gy=(Y[:,keep]@y)/n - (Y[:,keep].sum(axis=1).A.ravel()/n)*(y.mean())
    cov_gX=XtY/n - (Y[:,keep].sum(axis=1).A.ravel()/n)[:,None]*X.mean(axis=0)[None,:]
    cov_Xy=Xty/n - X.mean(axis=0)*y.mean()
    beta_X=X.mean(axis=0)@XtX_inv@cov_Xy
    num=cov_gy - cov_gX@XtX_inv@cov_Xy
    var_g=YtY/n - (Y[:,keep].sum(axis=1).A.ravel()/n)**2
    var_y=yty/n - y.mean()**2
    adj_g=var_g - np.einsum('ij,jk,ik->i', cov_gX, XtX_inv, cov_gX)
    adj_y=var_y - X.mean(axis=0)@XtX_inv@X.mean(axis=0)*0 - cov_Xy@XtX_inv@cov_Xy
    with np.errstate(divide="ignore", invalid="ignore"):
        r=num/np.sqrt(adj_g*adj_y)
    r[~np.isfinite(r)]=0.0
    order=np.argsort(-np.abs(r))
    out=os.path.join(OUT, f"apoe_coexpression_adj_{ctn.replace(' ','_')}.csv")
    with open(out,"w") as fh:
        fh.write("gene,partial_r,n_cells\n")
        for i in order:
            if genes[i]=="APOE": continue
            fh.write(f"{genes[i]},{r[i]:.6f},{n}\n")
    top=[(genes[i],r[i]) for i in order[:12] if genes[i]!="APOE"]
    print(f"{ctn:14s} n={n:6d} | " + ", ".join(f"{g}({v:+.2f})" for g,v in top), flush=True)
print("done")
