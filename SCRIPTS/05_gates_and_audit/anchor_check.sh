#!/usr/bin/env bash
# anchor_check.sh — 黄金锚点预检。任何新 MR 脚本（或 mr_lib 改动）上线前必跑。
# 断言：APOE(clumped 10Mb) → AMD IVW = +0.18177715 (p=7.89243e-11)  [FinnGen R13]
#        APOE(clumped 10Mb) → AD  IVW = -0.94183915 (p=4.08201e-18) [FinnGen R13]
# 容差：|Δ| < 1e-6（确定性管线应逐位一致）
set -uo pipefail
/home/huan/software/miniconda3/envs/biomni_r/bin/Rscript -e '
source("/home/huan/projects/08_APOE_AMD_AD_MR/revision/mr_lib.R")
pq <- fread("/home/huan/database/pQTL/pqtl/ukb-ppp/summary/ukb-ppp-2023-5e-08-clumped-10000kb-0.001r2.csv")
iv <- pq[gene.exposure == "APOE"]
iv <- iv[, .(SNP, qtl_effect = effect_allele.exposure, qtl_other = other_allele.exposure,
             beta = beta.exposure, se = se.exposure, p = pval.exposure)]
amd <- read_ma("/home/huan/projects/08_APOE_AMD_AD_MR/smr_data/AMD_r13.ma", TRUE)
ad  <- read_ma("/home/huan/projects/08_APOE_AMD_AD_MR/smr_data/AD_r13.ma", TRUE)
r_amd <- run_mr_pair(iv, amd, "APOE_AMD")
r_ad  <- run_mr_pair(iv, ad,  "APOE_AD")
g1 <- r_amd[Method == "IVW"]; g2 <- r_ad[Method == "IVW"]
stopifnot(nrow(g1) == 1, nrow(g2) == 1)
ok1 <- abs(g1$Beta - 0.181777146570588) < 1e-6 && abs(g1$P - 7.89243476286554e-11) < 1e-6
ok2 <- abs(g2$Beta - (-0.941839148400983)) < 1e-6 && abs(g2$P - 4.08201490071021e-18) < 1e-6
cat(sprintf("ANCHOR AMD: beta=%s p=%s -> %s\n", format(g1$Beta, digits=10), format(g1$P, digits=6), ifelse(ok1, "PASS", "FAIL")))
cat(sprintf("ANCHOR AD:  beta=%s p=%s -> %s\n", format(g2$Beta, digits=10), format(g2$P, digits=6), ifelse(ok2, "PASS", "FAIL")))
if (!(ok1 && ok2)) quit(status = 1)
cat("ALL ANCHORS PASS\n")
'
