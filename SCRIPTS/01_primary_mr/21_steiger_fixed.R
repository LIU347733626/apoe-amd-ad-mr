#!/usr/bin/env Rscript
# Correct Steiger directionality test using matched clumped APOE instruments.

library(data.table)

DB <- "/home/huan/database/pQTL/pqtl/ukb-ppp/summary"
AMD <- "/home/huan/projects/08_APOE_AMD_AD_MR/smr_data/AMD_finn.ma"
AD <- "/home/huan/projects/08_APOE_AMD_AD_MR/smr_data/AD_Bellenguez2022.ma"
OUT <- "/home/huan/projects/08_APOE_AMD_AD_MR/revision/results"

pq <- fread(file.path(DB, "ukb-ppp-2023-5e-08-clumped-10000kb-0.001r2.csv"))
apoe <- pq[gene.exposure == "APOE"]

read_ma <- function(path, header = TRUE) {
  if (header) {
    d <- fread(path, header = TRUE)
    setnames(d, c("SNP", "A1", "A2", "freq", "b", "se", "p", "n"))
  } else {
    d <- fread(path, header = FALSE,
               colClasses = c("character", "character", "character", "numeric",
                              "numeric", "numeric", "numeric", "integer"))
    setnames(d, c("SNP", "A1", "A2", "freq", "b", "se", "p", "n"))
  }
  d[!is.na(b) & !is.na(se) & se > 0 & !is.na(p)]
}

amd <- read_ma(AMD)
ad <- read_ma(AD, header = FALSE)

harm <- function(qa1, qa2, beta, ga1, ga2) {
  qa1 <- toupper(qa1); qa2 <- toupper(qa2)
  ga1 <- toupper(ga1); ga2 <- toupper(ga2)
  fifelse(qa1 == ga1 & qa2 == ga2, beta,
          fifelse(qa1 == ga2 & qa2 == ga1, -beta, NA_real_))
}

steiger_one <- function(outcome, label, n_outcome) {
  m <- merge(apoe, outcome[, .(SNP, A1, A2, b, se, p, n)], by = "SNP")
  m[, beta_harm := harm(effect_allele.exposure, other_allele.exposure,
                        beta.exposure, A1, A2)]
  m <- m[!is.na(beta_harm) & !is.na(se.exposure) & !is.na(se) &
           se.exposure > 0 & se > 0]
  if (nrow(m) < 3) {
    cat(sprintf("%s: insufficient matched instruments\n", label))
    return(NULL)
  }
  zx <- m$beta_harm / m$se.exposure
  zy <- m$b / m$se
  ## 2026-09-10 修正（审计 F4/F5/F6）：官方 Steiger 的多工具 R² = Σ r_i²（Hemani 2017；
  ## TwoSampleMR::mr_steiger 源码 `r_exp <- sqrt(sum(r_exp^2))`），原用 mean() 即人均 r²，
  ## 导致 z 不在官方尺度上（差 ~k 倍）。证据：docs/audit_20260909/CODE_AUDIT_SELFWRITTEN_20260910.md
  r2x <- sum(zx^2 / (zx^2 + m$samplesize.exposure - 2))
  r2y <- sum(zy^2 / (zy^2 + m$n - 2))
  n_x <- mean(m$samplesize.exposure)
  n_y <- mean(m$n)
  z_steiger <- (atanh(sqrt(r2x)) - atanh(sqrt(r2y))) /
    sqrt(1 / (n_x - 3) + 1 / (n_y - 3))
  p_steiger <- 2 * pnorm(-abs(z_steiger))
  direction <- ifelse(r2x > r2y & p_steiger < 0.05, "Forward",
                      ifelse(r2y > r2x & p_steiger < 0.05, "Reverse", "Inconclusive"))
  cat(sprintf("%s: r2_exposure=%.6f, r2_outcome=%.6f, z=%.3f, p=%.2e, direction=%s\n",
              label, r2x, r2y, z_steiger, p_steiger, direction))
  data.table(Test = label, N_SNP = nrow(m), R2_exposure = r2x,
             R2_outcome = r2y, Z = z_steiger, P = p_steiger,
             Direction = direction)
}

r_amd <- steiger_one(amd, "APOE_pQTL_to_AMD", NA)
r_ad <- steiger_one(ad, "APOE_pQTL_to_AD", NA)
res <- rbindlist(list(r_amd, r_ad), fill = TRUE)
write.csv(res, file.path(OUT, "steiger_fixed.csv"), row.names = FALSE)
print(res)
cat("\nDONE\n")
