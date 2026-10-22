#!/usr/bin/env Rscript
# Revision analysis: APOE focus MR using clumped pQTL instruments only.
# Addresses reviewer concern that the original _03 used unclumped APOE-region SNPs.

library(data.table)
library(MendelianRandomization)
library(MRPRESSO)

PROJ <- "/home/huan/projects/08_APOE_AMD_AD_MR"
DB <- "/home/huan/database/pQTL/pqtl/ukb-ppp/summary"
AMD <- "/home/huan/projects/08_APOE_AMD_AD_MR/smr_data/AMD_finn.ma"
AD <- "/home/huan/projects/08_APOE_AMD_AD_MR/smr_data/AD_Bellenguez2022.ma"
OUT <- file.path(PROJ, "revision", "results")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

pq <- fread(file.path(DB, "ukb-ppp-2023-5e-08-clumped-10000kb-0.001r2.csv"),
            header = TRUE, sep = ",", data.table = TRUE)
apoe <- pq[gene.exposure == "APOE", ]
cat(sprintf("Clumped APOE instruments: %d\n", nrow(apoe)))

read_ma <- function(path, header = TRUE) {
  if (header) {
    d <- fread(path, header = TRUE, data.table = TRUE)
    setnames(d, c("SNP", "A1", "A2", "freq", "b", "se", "p", "n"))
  } else {
    d <- fread(path, header = FALSE, data.table = TRUE,
               colClasses = c("character", "character", "character", "numeric",
                              "numeric", "numeric", "numeric", "integer"))
    setnames(d, c("SNP", "A1", "A2", "freq", "b", "se", "p", "n"))
  }
  d[!is.na(b) & !is.na(se) & se > 0]
}

amd <- read_ma(AMD, header = TRUE)
ad <- read_ma(AD, header = FALSE)

harmonize <- function(qtl_a1, qtl_a2, qtl_beta, gwas_a1, gwas_a2) {
  qa1 <- toupper(qtl_a1); qa2 <- toupper(qtl_a2)
  ga1 <- toupper(gwas_a1); ga2 <- toupper(gwas_a2)
  fifelse(qa1 == ga1 & qa2 == ga2, qtl_beta,
          fifelse(qa1 == ga2 & qa2 == ga1, -qtl_beta, NA_real_))
}

merge_iv <- function(outcome, outcome_name) {
  m <- merge(apoe, outcome[, .(SNP, A1, A2, b, se)], by = "SNP",
             suffixes = c(".pqtl", ".gwas"))
  m[, beta_harm := harmonize(effect_allele.exposure, other_allele.exposure,
                             beta.exposure, A1, A2)]
  m <- m[!is.na(beta_harm) & !is.na(se.exposure) & !is.na(se) &
           se.exposure > 0 & se > 0]
  cat(sprintf("%s: %d instruments after merge/harmonization\n",
              outcome_name, nrow(m)))
  m
}

run_mr <- function(bx, bxse, by, byse, snps, label) {
  obj <- mr_input(bx = bx, bxse = bxse, by = by, byse = byse, snps = snps)
  ivw <- tryCatch(mr_ivw(obj), error = function(e) NULL)
  egger <- tryCatch(mr_egger(obj), error = function(e) NULL)
  median <- tryCatch(mr_median(obj), error = function(e) NULL)
  presso <- tryCatch(
    # FIX(2026-09-10): mr_presso 的 data= 无默认值；原调用缺该参数 → 报错被
    # tryCatch 静默吞掉 → MR-PRESSO 从未运行。且 BetaOutcome 等须为 data 中的
    # 列名（原代码传的是向量）。
    mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
              SdOutcome = "sd.outcome", SdExposure = "sd.exposure",
              data = data.frame(beta.outcome = by, beta.exposure = bx,
                                sd.outcome = byse, sd.exposure = bxse),
              OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
              NbDistribution = 1000, SignifThreshold = 0.05,
              seed = 20260910),
    error = function(e) NULL
  )
  loo <- tryCatch(mr_leaveoneout(obj), error = function(e) NULL)

  rows <- list()
  if (!is.null(ivw)) {
    rows[[length(rows) + 1]] <- data.table(
      Method = paste0(label, " IVW"), N_SNP = length(bx),
      Beta = ivw@Estimate, SE = ivw@StdError, P = ivw@Pvalue,
      Hetero_Q = ivw@Heter.Stat[1], Hetero_p = ivw@Heter.Stat[2])
  }
  if (!is.null(egger)) {
    rows[[length(rows) + 1]] <- data.table(
      Method = paste0(label, " Egger"), N_SNP = length(bx),
      Beta = egger@Estimate, SE = egger@StdError.Est, P = egger@Pvalue.Est,
      Egger_Intercept = egger@Intercept, Egger_InterceptSE = egger@StdError.Int,
      Egger_InterceptP = egger@Pvalue.Int)
  }
  if (!is.null(median)) {
    rows[[length(rows) + 1]] <- data.table(
      Method = paste0(label, " Weighted_median"), N_SNP = length(bx),
      Beta = median@Estimate, SE = median@StdError, P = median@Pvalue)
  }
  if (!is.null(presso) && length(presso[["Main MR results"]]) > 0) {
    mr <- presso[["Main MR results"]]
    rows[[length(rows) + 1]] <- data.table(
      Method = paste0(label, " MR-PRESSO_raw"), N_SNP = length(bx),
      Beta = mr[["Causal Estimate"]][1], SE = mr[["Sd"]][1], P = mr[["P-value"]][1])
    if (nrow(mr) >= 2) {
      rows[[length(rows) + 1]] <- data.table(
        Method = paste0(label, " MR-PRESSO_corrected"), N_SNP = length(bx),
        Beta = mr[["Causal Estimate"]][2], SE = mr[["Sd"]][2], P = mr[["P-value"]][2])
    }
  }
  if (!is.null(loo)) {
    loo_beta <- sapply(loo, function(x) x@Estimate)
    rows[[length(rows) + 1]] <- data.table(
      Method = paste0(label, " Leave-one-out_range"), N_SNP = length(bx),
      Beta = mean(loo_beta), SE = sd(loo_beta),
      P = NA_real_, LOO_min = min(loo_beta), LOO_max = max(loo_beta))
  }
  r <- rbindlist(rows, fill = TRUE)
  # Guard against P=0 underflow in the MendelianRandomization package
  # (2*pnorm(-|z|) underflows to 0 for |z|>~38). Report a computed floor
  # (z-test approximation) instead of a literal 0.
  if (nrow(r) > 0 && "P" %in% names(r)) {
    r[, P := fifelse(P == 0, pmax(2*pnorm(-abs(Beta/SE)), 1e-300), P)]
  }
  r
}

m_amd <- merge_iv(amd, "AMD")
m_ad <- merge_iv(ad, "AD")
r_amd <- run_mr(m_amd$beta_harm, m_amd$se.exposure, m_amd$b, m_amd$se,
                m_amd$SNP, "APOE_clumped_AMD")
r_ad <- run_mr(m_ad$beta_harm, m_ad$se.exposure, m_ad$b, m_ad$se,
               m_ad$SNP, "APOE_clumped_AD")
res <- rbindlist(list(r_amd, r_ad), fill = TRUE)
write.csv(res, file.path(OUT, "apoe_focus_mr_clumped.csv"), row.names = FALSE)
print(res)

cat("\nDONE\n")
