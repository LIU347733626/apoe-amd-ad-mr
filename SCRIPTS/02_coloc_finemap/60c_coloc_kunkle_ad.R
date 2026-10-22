#!/usr/bin/env Rscript
# 60c: Kunkle 2019 AD APOE-region coloc (rsID-complete) vs UKB-PPP APOE cis-pQTL.
# Rationale: Bellenguez 2022 sumstats lacks rs429358 / level-regulating key variants
# (rsID intersection with the pQTL set is partial), so its coloc PP3=1.0 is an
# alignment artifact, not a biological "distinct". Kunkle 2019 has the full APOE
# region (rs429358 / rs6859 / rs157580 / rs7412 present); this run tests whether
# Kunkle gives a reliable PP4 (shared) for AD.
# D1 = UKB-PPP APOE cis-pQTL (identical variant set to revision/60).
# D2 = Kunkle 2019 AD, APOE-region SNPs, case fraction s = 46601/94437.
# NOTE: Kunkle rs429358 has the known effect-allele sign anomaly (Effect=T, beta<0);
# coloc uses PP (sign-agnostic) so this does NOT affect the shared/distinct call.
# Output: revision/results/r13/coloc_kunkle_ad.csv
suppressPackageStartupMessages({ library(data.table); library(coloc) })

PROJ <- "/home/huan/projects/08_APOE_AMD_AD_MR"
OUT  <- file.path(PROJ, "revision/results/r13")
PQTL <- "/home/huan/database/pQTL/pqtl/ukb-ppp/summary/ukb-ppp-2023-5e-08-unclumped.csv"
KDAT <- file.path(PROJ, "smr_data/Kunkle_2019_AD.txt")
N_K  <- 94437; CA_K <- 46601; s <- CA_K/N_K

pq <- fread(PQTL, select=c("SNP","effect_allele","other_allele","eaf","beta","se",
                           "pval","pos_hg19","gene","samplesize"))
apoe <- pq[gene=="APOE"]
setnames(apoe, c("beta","se","pval","eaf","pos_hg19","samplesize"),
         c("beta_pqtl","se_pqtl","pval_pqtl","eaf_pqtl","pos_hg19_pqtl","samplesize_pqtl"))
cat("APOE cis-pQTL SNPs:", nrow(apoe), "\n")
snp_file <- file.path(OUT, "apoe_coloc_snps_60.txt"); writeLines(apoe$SNP, snp_file)

# Kunkle: space-delimited, has header
# FIX(2026-09-10): grep 只输出匹配行，表头行不匹配 rsID → 输出中本无表头；
# 原 header=TRUE 会把首个匹配变异（rs62116612）当列名吃掉。
k <- fread(cmd=sprintf("grep -w -F -f %s %s", snp_file, KDAT), header=FALSE,
           colClasses=c("integer","integer","character","character","character",
                        "numeric","numeric","numeric"))
setnames(k, c("Chromosome","Position","MarkerName","Effect_allele","Non_Effect_allele","b","se","p"))
k <- k[!is.na(b) & !is.na(se) & se>0]
setnames(k, "MarkerName", "SNP")
cat("Kunkle matched SNPs:", nrow(k), "\n")

m <- merge(apoe, k, by="SNP")
m[, beta_harm := fifelse(toupper(effect_allele)==toupper(Effect_allele) &
                           toupper(other_allele)==toupper(Non_Effect_allele), beta_pqtl,
                         fifelse(toupper(effect_allele)==toupper(Non_Effect_allele) &
                                   toupper(other_allele)==toupper(Effect_allele), -beta_pqtl, NA_real_))]
# Kunkle file has no allele-frequency column; use the pQTL (UKB-PPP) EAF as the
# population MAF proxy for the same SNP (coloc MAF is a filter/computation input).
m[, maf := ifelse(eaf_pqtl<=0.5, eaf_pqtl, 1-eaf_pqtl)]
m <- m[!is.na(beta_harm) & se_pqtl>0 & se>0 & !duplicated(SNP) & maf>0 & maf<1]
cat("Harmonized SNPs:", nrow(m), "\n")
if (nrow(m) < 10) stop("too few harmonized SNPs", call.=FALSE)

p_out <- 2*pnorm(-abs(m$b/m$se))
D1 <- list(type="quant", N=unique(m$samplesize_pqtl)[1], sdY=1,
           snp=m$SNP, pvalues=m$pval_pqtl, beta=m$beta_harm,
           varbeta=m$se_pqtl^2, MAF=m$maf, position=m$pos_hg19_pqtl)
D2 <- list(type="cc", N=N_K, s=s, snp=m$SNP, pvalues=p_out,
           beta=m$b, varbeta=m$se^2, MAF=m$maf, position=m$pos_hg19_pqtl)

fit <- coloc.abf(D1, D2); ssum <- fit$summary
out <- data.table(cohort="Kunkle", disease="AD", N_SNP=nrow(m),
                  PP3=as.numeric(ssum["PP.H3.abf"]), PP4=as.numeric(ssum["PP.H4.abf"]),
                  PP.H0=as.numeric(ssum["PP.H0.abf"]), PP.H1=as.numeric(ssum["PP.H1.abf"]),
                  PP.H2=as.numeric(ssum["PP.H2.abf"]), s_outcome=s)
fwrite(out, file.path(OUT, "coloc_kunkle_ad.csv"), col.names=TRUE)
print(out)
cat("\nKunkle AD coloc done. PP3=%.4f PP4=%.6f N=%d\n", out$PP3, out$PP4, out$N_SNP)
