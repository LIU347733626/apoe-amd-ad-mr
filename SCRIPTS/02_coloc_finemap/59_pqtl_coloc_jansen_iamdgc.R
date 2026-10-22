#!/usr/bin/env Rscript
# 59: APOE pQTL x {Jansen 2019 AD, IAMDGC late AMD} colocalization (coloc.abf)
#
# D1 = UKB-PPP APOE cis-pQTL (unclumped, same as revision/23 — the R13-AMD
#      coloc provenance)
# D2 = outcome GWAS APOE-region SNPs (harmonized, per-SNP N)
#
# Expected pattern (per R13 mechanism story):
#   AMD (IAMDGC):  PP4 high  — same causal signal, APOE-level mediated
#   AD  (Jansen):  PP3 high  — pQTL (level variants) separate from the AD
#                         signal concentrated in the epsilon4 cluster
#
# Data caveats (documented in manuscript):
#   - Jansen: AD-by-proxy composition; case fraction s=71,880/436,498
#     (published case count / file Nsum); per-allele scale differs from
#     FinnGen (z-score based inference below is scale-invariant).
#   - IAMDGC: 7 flipped variants in chr1/chr10 complement regions do NOT
#     affect chr19 APOE; case fraction s=16,144/32,339 (Fritsche 2016:
#     16,144 late AMD cases; file per-SNP N=32,339).
#   - pvalues for coloc recomputed from b/se (2*pnorm(|z|)) so extreme
#     p=0 clamp artifacts in the converted Jansen file cannot perturb z.
#
# Outputs: revision/results/jansen/apoe_pqtl_coloc_jansen_ad.csv
#          revision/results/jansen/apoe_pqtl_coloc_iamdgc_amd.csv
##############################################################################
suppressPackageStartupMessages({
  library(data.table); library(coloc)
})

PROJ <- "/home/huan/projects/08_APOE_AMD_AD_MR"
PQTL <- "/home/huan/database/pQTL/pqtl/ukb-ppp/summary/ukb-ppp-2023-5e-08-unclumped.csv"
OUT  <- file.path(PROJ, "revision/results/jansen")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

cat("===== 59: APOE pQTL coloc x {Jansen AD, IAMDGC AMD} =====\n")

# 1. D1: APOE cis-pQTL (identical to revision/23) --------------------------
pq <- fread(PQTL, select = c("SNP", "effect_allele", "other_allele",
                              "eaf", "beta", "se", "pval", "pos_hg19",
                              "gene", "samplesize"))
apoe <- pq[gene == "APOE"]
cat(sprintf("APOE cis-pQTL SNPs: %d\n", nrow(apoe)))
setnames(apoe, c("beta", "se", "pval", "eaf", "pos_hg19", "samplesize"),
         c("beta_pqtl", "se_pqtl", "pval_pqtl", "eaf_pqtl",
           "pos_hg19_pqtl", "samplesize_pqtl"))
snp_file <- file.path(OUT, "apoe_coloc_snps_59.txt")
writeLines(apoe$SNP, snp_file)

# 2. Runner -----------------------------------------------------------------
run_coloc <- function(tag, gwas_file, use_header, cases, n_total) {
  cat(sprintf("\n##### %s (N=%d, cases~%d, s=%.4f)\n", tag, n_total, cases, cases / n_total))
  g <- fread(cmd = sprintf("grep -w -F -f %s %s", snp_file, gwas_file),
  # FIX(2026-09-10): grep 只输出匹配行，表头行不匹配 rsID → 输出中恒无表头。
  # 原按『文件是否有表头』传参，传 TRUE 时会静默吃掉首个匹配变异（rs62116612）。
             header = FALSE,
             colClasses = c("character", "character", "character", "numeric",
                            "numeric", "numeric", "numeric", "integer"))
  setnames(g, c("SNP", "A1", "A2", "freq", "b", "se", "p", "n"))
  g <- g[!is.na(b) & !is.na(se) & se > 0]
  cat(sprintf("Matched outcome SNPs: %d\n", nrow(g)))
  m <- merge(apoe, g, by = "SNP")
  m[, beta_harm := fifelse(toupper(effect_allele) == toupper(A1) &
                             toupper(other_allele) == toupper(A2), beta_pqtl,
                           fifelse(toupper(effect_allele) == toupper(A2) &
                                     toupper(other_allele) == toupper(A1),
                                   -beta_pqtl, NA_real_))]
  m <- m[!is.na(beta_harm) & se_pqtl > 0 & se > 0 & !duplicated(SNP) &
           !is.na(freq) & freq > 0 & freq < 1 &
           !is.na(eaf_pqtl) & eaf_pqtl > 0 & eaf_pqtl < 1]
  cat(sprintf("Harmonized SNPs (with MAF): %d\n", nrow(m)))
  if (nrow(m) < 5) stop(tag, ": too few harmonized SNPs for coloc", call. = FALSE)

  # recomputed p from z (robust to p=0 clamp in converted files)
  p_out <- 2 * pnorm(-abs(m$b / m$se))

  D1 <- list(type = "quant", N = unique(m$samplesize_pqtl)[1], sdY = 1,
             snp = m$SNP, pvalues = m$pval_pqtl, beta = m$beta_harm,
             varbeta = m$se_pqtl^2,
             MAF = ifelse(m$eaf_pqtl <= 0.5, m$eaf_pqtl, 1 - m$eaf_pqtl),
             position = m$pos_hg19_pqtl)
  D2 <- list(type = "cc", N = n_total, s = cases / n_total,
             snp = m$SNP, pvalues = p_out, beta = m$b, varbeta = m$se^2,
             MAF = ifelse(m$freq <= 0.5, m$freq, 1 - m$freq),
             position = m$pos_hg19_pqtl)

  fit <- coloc.abf(D1, D2)
  s <- fit$summary
  out <- data.table(
    Module = tag, Gene = "APOE", N_SNP = nrow(m),
    PP.H0 = s["PP.H0.abf"], PP.H1 = s["PP.H1.abf"], PP.H2 = s["PP.H2.abf"],
    PP.H3 = s["PP.H3.abf"], PP.H4 = s["PP.H4.abf"])
  fwrite(out, file.path(OUT, paste0("apoe_pqtl_coloc_", tag, ".csv")), col.names = TRUE)
  print(out)
  out
}

res_jansen  <- run_coloc("jansen_ad",  file.path(PROJ, "smr_data/AD_Jansen2019_converted.ma"),
                         use_header = FALSE, cases = 71880, n_total = 436498)
res_iamdgc  <- run_coloc("iamdgc_amd", "/home/huan/database/GWAS/AMD/IAMDGC_late.ma",
                         use_header = TRUE,  cases = 16144, n_total = 32339)

cat("\n=== Summary ===\n")
cat(sprintf("Jansen AD  : PP3=%.4f  PP4=%.4f  (expect PP3 high: eps4-driven AD signal separate from level pQTL)\n",
            res_jansen$PP.H3, res_jansen$PP.H4))
cat(sprintf("IAMDGC AMD : PP3=%.4f  PP4=%.4f  (expect PP4 high: same signal, level-mediated; R13 anchor PP4=0.9998)\n",
            res_iamdgc$PP.H3, res_iamdgc$PP.H4))
cat("\nDone.\n")
