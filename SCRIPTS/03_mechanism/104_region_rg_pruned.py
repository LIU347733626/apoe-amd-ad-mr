#!/usr/bin/env python3
"""104_region_rg_pruned.py — APOE 区域 vs 全基因组：AMD/AD 效应一致性（LD 剪枝版）
一次性读入两个 .ma，提取区域与 HM3 剪枝集，算 z 相关 + 符号一致性。"""
import numpy as np, math, os
REG = set(l.strip() for l in open("revision/tmp_region_rg/region_prune.prune.in"))
GEN = set(l.strip() for l in open("revision/tmp_region_rg/genome_prune.prune.in"))
WANT = REG | GEN
def load(path, header):
    out = {}
    with open(path) as fh:
        if header: next(fh)
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 8 or p[0] not in WANT: continue
            try: b = float(p[4]); se = float(p[5])
            except ValueError: continue
            if se <= 0: continue
            out[p[0]] = (p[1].upper(), p[2].upper(), b, se)
    return out
print("读 AMD ...", flush=True); amd = load("smr_data/AMD_r13.ma", True)
print("读 AD ...", flush=True);  ad  = load("smr_data/AD_r13.ma", False)
print(f"载入 AMD {len(amd)} | AD {len(ad)}", flush=True)
def harm(d1, d2):
    rows = []
    for s, (a1, a2, b1, se1) in d1.items():
        if s not in d2: continue
        c1, c2, b2, se2 = d2[s]
        if a1 == c1 and a2 == c2: rows.append((s, b1/se1, b2/se2))
        elif a1 == c2 and a2 == c1: rows.append((s, b1/se1, -b2/se2))
    return rows
def report(rows, label):
    z1 = np.array([r[1] for r in rows]); z2 = np.array([r[2] for r in rows])
    r = float(np.corrcoef(z1, z2)[0, 1]); n = len(rows)
    zf = 0.5*math.log((1+r)/(1-r)); se = 1/math.sqrt(n-3)
    lo, hi = math.tanh(zf-1.96*se), math.tanh(zf+1.96*se)
    opp = int(((z1*z2) < 0).sum())
    # 二项检验（正态近似，n 很大时精确组合数不可行）
    from scipy.stats import binomtest
    p = binomtest(opp, n, 0.5, alternative="two-sided").pvalue
    print(f"{label}: n={n}  r={r:+.4f}  Fisher95%CI {lo:+.3f}..{hi:+.3f}  |  反向符号 {opp}/{n} ({100*opp/n:.1f}%)  二项P={p:.2e}")
    return r, n
r1, n1 = report(harm({k: v for k, v in amd.items() if k in REG}, {k: v for k, v in ad.items() if k in REG}), "APOE 区域（剪枝 r²<0.2）")
r2, n2 = report(harm({k: v for k, v in amd.items() if k in GEN}, {k: v for k, v in ad.items() if k in GEN}), "全基因组 HM3（剪枝 r²<0.2）")
print(f"\n|r_region| / |r_genome| = {abs(r1)/abs(r2):.1f}x")
print("done")
