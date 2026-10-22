#!/usr/bin/env Rscript
# 60: APOE pQTL x {disease x cohort} colocalization (coloc.abf), full-cis variant set
#
# D1 = UKB-PPP APOE cis-pQTL (unclumped; identical variant set to revision/23 / 23b / 59)
# D2 = outcome GWAS APOE-region SNPs (harmonized, per-cohort case fraction)
#
# Cohorts (disease x cohort), 3 per disease:
#   AMD: IAMDGC (already in 59), FinnGen R13 (23_r13 area), MVP Phe_362_2 (NEW)
#   AD : Bellenguez 2022 (NEW), FinnGen R13 (already in 23b), Jansen 2019 (already in 59)
#
# Case fractions (documented values; reused from prior revision scripts where computed):
#   IAMDGC   : s = 16144 / 32339   (Fritsche 2016 late-AMD cases; 59)
#   AMD R13  : s = 5500  / 472750  (project-convention R8-era value; robust band runs
#                                   in 59b to ~16k gave identical PP4 ~ 0.974)
#   MVP      : s = 70936 / (70936+536119)  (MVP R4 META Phe_362_2 header case/control counts)
#   Bellenguez: s = 111326 / (111326+677663)  (Bellenguez 2022 published counts)
#   AD R13   : s = 16652 / 500186  (23b value)
#   Jansen   : s = 71880 / 436498  (59 value)
#
# p recomputed from b/se (robust to p=0 clamp in converted files).
#
# Output: revision/results/r13/coloc_bydisease.csv
##############################################################################
suppressPackageStartupMessages({
  library(data.table); library(coloc)
})

PROJ  <- "/home/huan/projects/08_APOE_AMD_AD_MR"
PQTL  <- "/home/huan/database/pQTL/pqtl/ukb-ppp/summary/ukb-ppp-2023-5e-08-unclumped.csv"
OUT   <- file.path(PROJ, "revision/results/r13"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

cat("===== 60: APOE pQTL coloc x {disease x cohort} =====\n")

# ---------- D1: APOE cis-pQTL (identical to 23/23b/59) ----------
pq <- fread(PQTL, select = c("SNP", "effect_allele", "other_allele",
                              "eaf", "beta", "se", "pval", "pos_hg19",
                              "gene", "samplesize"))
apoe <- pq[gene == "APOE"]
setnames(apoe, c("beta", "se", "pval", "eaf", "pos_hg19", "samplesize"),
         c("beta_pqtl", "se_pqtl", "pval_pqtl", "eaf_pqtl",
           "pos_hg19_pqtl", "samplesize_pqtl"))
cat(sprintf("APOE cis-pQTL SNPs: %d\n", nrow(apoe)))
snp_file <- file.path(OUT, "apoe_coloc_snps_60.txt"); writeLines(apoe$SNP, snp_file)

# ---------- outcome readers -> .(SNP,A1,A2,freq,b,se) ----------
read_ma <- function(p, header) {
  # FIX(2026-09-10): grep 输出恒无表头（表头行不匹配 rsID）；原传 header=TRUE 会吃掉 rs62116612。
  d <- fread(cmd = sprintf("grep -w -F -f %s %s", snp_file, p), header = FALSE,
             colClasses = c("character","character","character","numeric",
                            "numeric","numeric","numeric","integer"))
  setnames(d, c("SNP","A1","A2","freq","b","se","p","n"))
  d[!is.na(b) & !is.na(se) & se > 0]
}
# MVP reader (Phe_362_2 META): ea=effect allele, or=OR per-ea, n from header case/control
# .gz -> must zgrep; grep on binary drops header, so read header=FALSE and set names by position
MVP_COLS <- c("SNP_ID","chrom","pos","ref","alt","ea","af","num_samples","case_af",
              "num_cases","control_af","num_controls","or","ci","pval","r2","q_pval","i2","direction")
read_mvp <- function(p) {
  z <- fread(cmd = sprintf("zcat %s | grep -w -F -f %s", p, snp_file), header = FALSE)
  setnames(z, MVP_COLS)
  z <- z[, .(SNP=SNP_ID, A1=ea, A2=ref, freq=as.numeric(af),
             b=log(as.numeric(or)),
             se=(log(as.numeric(sub(".*,", "", ci))) - log(as.numeric(sub(",.*", "", ci)))) / 3.92,
             p=as.numeric(pval))]
  z[!is.na(b) & !is.na(se) & se > 0 & freq > 0 & freq < 1]
}

# ---------- coloc runner ----------
run_coloc <- function(disease, cohort, gwas, header, N, cases) {
  cat(sprintf("\n##### %s | %s  (N=%d, cases=%d, s=%.5f)\n",
              disease, cohort, N, cases, cases / N))
  g <- if (cohort == "MVP") read_mvp(gwas) else read_ma(gwas, header)
  cat(sprintf("Matched outcome SNPs: %d\n", nrow(g)))
  m <- merge(apoe, g, by = "SNP")
  m[, beta_harm := fifelse(toupper(effect_allele) == toupper(A1) &
                             toupper(other_allele) == toupper(A2), beta_pqtl,
                           fifelse(toupper(effect_allele) == toupper(A2) &
                                     toupper(other_allele) == toupper(A1),
                                   -beta_pqtl, NA_real_))]
  m <- m[!is.na(beta_harm) & se_pqtl > 0 & se > 0 & !duplicated(SNP) &
           !is.na(freq) & freq > 0 & freq < 1 & !is.na(eaf_pqtl) &
           eaf_pqtl > 0 & eaf_pqtl < 1]
  cat(sprintf("Harmonized SNPs (with MAF): %d\n", nrow(m)))
  if (nrow(m) < 5) stop(disease, "|", cohort, ": too few harmonized SNPs", call. = FALSE)

  p_out <- 2 * pnorm(-abs(m$b / m$se))   # recompute p from z (robust to p=0 clamp)

  D1 <- list(type = "quant", N = unique(m$samplesize_pqtl)[1], sdY = 1,
             snp = m$SNP, pvalues = m$pval_pqtl, beta = m$beta_harm,
             varbeta = m$se_pqtl^2,
             MAF = ifelse(m$eaf_pqtl <= 0.5, m$eaf_pqtl, 1 - m$eaf_pqtl),
             position = m$pos_hg19_pqtl)
  D2 <- list(type = "cc", N = N, s = cases / N,
             snp = m$SNP, pvalues = p_out, beta = m$b, varbeta = m$se^2,
             MAF = ifelse(m$freq <= 0.5, m$freq, 1 - m$freq),
             position = m$pos_hg19_pqtl)

  fit <- coloc.abf(D1, D2); s <- fit$summary
  out <- data.table(disease = disease, cohort = cohort,
                    N_SNP = nrow(m),
                    PP3 = as.numeric(s["PP.H3.abf"]), PP4 = as.numeric(s["PP.H4.abf"]),
                    PP.H0 = as.numeric(s["PP.H0.abf"]), PP.H1 = as.numeric(s["PP.H1.abf"]),
                    PP.H2 = as.numeric(s["PP.H2.abf"]),
                    s_outcome = cases / N)
  cat(sprintf("  -> PP3=%.4f  PP4=%.4f  N_SNP=%d\n", out$PP3, out$PP4, out$N_SNP))
  out
}

res <- rbindlist(list(
  ## AMD
  run_coloc("AMD", "IAMDGC", "/home/huan/database/GWAS/AMD/IAMDGC_late.ma", TRUE,
            32339, 16144),
  run_coloc("AMD", "R13",   file.path(PROJ, "smr_data/AMD_r13.ma"), TRUE,
            472750, 5500),
  run_coloc("AMD", "MVP",   "/mnt/ssd4t/GWAS/AMD/MvP/extract/MVP_R4.1000G_AGR.Phe_362_2.META.GIA.dbGaP.txt.gz",
            FALSE, 607055, 70936),
  ## AD
  run_coloc("AD", "Bellenguez", file.path(PROJ, "smr_data/AD_Bellenguez2022.ma"), FALSE,
            788989, 111326),
  run_coloc("AD", "R13",   file.path(PROJ, "smr_data/AD_r13_hdr.ma"), TRUE,
            500186, 16652),
  run_coloc("AD", "Jansen", file.path(PROJ, "smr_data/AD_Jansen2019_converted.ma"), FALSE,
            436498, 71880)
), fill = TRUE)

fwrite(res, file.path(OUT, "coloc_bydisease.csv"), col.names = TRUE)
print(res[])
cat("\n=== Summary ===\n")
cat("AMD: ", paste(sprintf("%s=PP4 %.4f (PP3 %.4f, N=%d)", res[disease=="AMD"]$cohort,
      res[disease=="AMD"]$PP4, res[disease=="AMD"]$PP3, res[disease=="AMD"]$N_SNP), collapse="; "), "\n")
cat("AD:  ", paste(sprintf("%s=PP4 %.4f (PP3 %.4f, N=%d)", res[disease=="AD"]$cohort,
      res[disease=="AD"]$PP4, res[disease=="AD"]$PP3, res[disease=="AD"]$N_SNP), collapse="; "), "\n")
cat("\nWritten:", file.path(OUT, "coloc_bydisease.csv"), "\n")
cat("Done.\n")
