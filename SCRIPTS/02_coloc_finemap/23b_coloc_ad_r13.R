#!/usr/bin/env Rscript
# Reproducible APOE pQTL vs AMD coloc using UKB-PPP and FinnGen summary data.

library(data.table)
library(coloc)

DB <- "/home/huan/database/pQTL/pqtl/ukb-ppp/summary"
AD <- "/home/huan/projects/08_APOE_AMD_AD_MR/smr_data/AD_r13.ma"
OUT <- "/home/huan/projects/08_APOE_AMD_AD_MR/revision/results/r13"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

pq <- fread(file.path(DB, "ukb-ppp-2023-5e-08-unclumped.csv"),
            select = c("SNP", "effect_allele", "other_allele",
                       "eaf", "beta", "se",
                       "pval", "pos_hg19", "gene",
                       "samplesize"))
apoe <- pq[gene == "APOE"]
cat(sprintf("APOE cis-pQTL SNPs: %d\n", nrow(apoe)))
setnames(apoe,
         c("beta", "se", "pval", "eaf", "pos_hg19", "samplesize"),
         c("beta_pqtl", "se_pqtl", "pval_pqtl", "eaf_pqtl",
           "pos_hg19_pqtl", "samplesize_pqtl"))

snp_file <- file.path(OUT, "apoe_coloc_snps.txt")
writeLines(apoe$SNP, snp_file)

amd <- fread(cmd = sprintf("grep -w -F -f %s %s", snp_file, AD),
             header = FALSE,
             colClasses = c("character", "character", "character", "numeric",
                            "numeric", "numeric", "numeric", "integer"))
setnames(amd, c("SNP", "A1", "A2", "freq", "b", "se", "p", "n"))
amd <- amd[!is.na(b) & !is.na(se) & se > 0 & !is.na(p)]
cat(sprintf("Matched AD SNPs: %d\n", nrow(amd)))
amd[, c("b_gwas", "se_gwas", "p_gwas", "n_gwas") := .(b, se, p, n)]

m <- merge(apoe, amd, by = "SNP")
m[, beta_harm := fifelse(
  toupper(effect_allele) == toupper(A1) &
    toupper(other_allele) == toupper(A2),
  beta_pqtl,
  fifelse(toupper(effect_allele) == toupper(A2) &
            toupper(other_allele) == toupper(A1),
          -beta_pqtl, NA_real_))]
m <- m[!is.na(beta_harm) & !is.na(se_pqtl) & !is.na(se_gwas) &
         se_pqtl > 0 & se_gwas > 0 & !duplicated(SNP)]
cat(sprintf("Harmonized SNPs: %d\n", nrow(m)))

D1 <- list(
  type = "quant",
  N = unique(m$samplesize_pqtl)[1],
  sdY = 1,
  snp = m$SNP,
  pvalues = m$pval_pqtl,
  beta = m$beta_harm,
  varbeta = m$se_pqtl^2,
  MAF = ifelse(m$eaf_pqtl <= 0.5, m$eaf_pqtl, 1 - m$eaf_pqtl),
  position = m$pos_hg19_pqtl
)

D2 <- list(
  type = "cc",
  N = unique(m$n)[1],
  s = 16652 / unique(m$n)[1],
  snp = m$SNP,
  pvalues = m$p_gwas,
  beta = m$b_gwas,
  varbeta = m$se_gwas^2,
  MAF = ifelse(m$freq <= 0.5, m$freq, 1 - m$freq),
  position = m$pos_hg19_pqtl
)

cat("Running coloc.abf (comparable with R8 & AMD side)...\n")
resab <- coloc.abf(D1, D2)
print(resab$summary)

out <- data.table(
  Module = "pQTL",
  Gene = "APOE",
  N_SNP = nrow(m),
  PP.H0 = resab$summary["PP.H0.abf"],
  PP.H1 = resab$summary["PP.H1.abf"],
  PP.H2 = resab$summary["PP.H2.abf"],
  PP.H3 = resab$summary["PP.H3.abf"],
  PP.H4 = resab$summary["PP.H4.abf"],
  PP.H0.abf = resab$summary["PP.H0.abf"],
  PP.H1.abf = resab$summary["PP.H1.abf"],
  PP.H2.abf = resab$summary["PP.H2.abf"],
  PP.H3.abf = resab$summary["PP.H3.abf"],
  PP.H4.abf = resab$summary["PP.H4.abf"]
)
write.csv(out, file.path(OUT, "apoe_pqtl_coloc_ad.csv"), row.names = FALSE)
print(out)
cat("\nDONE\n")
