#!/usr/bin/env Rscript
# Harmonized APOE safety MR using clumped instruments.

library(data.table)
library(MendelianRandomization)

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
  d[!is.na(b) & !is.na(se) & se > 0]
}

amd <- read_ma(AMD)
ad <- read_ma(AD, header = FALSE)

harm <- function(qa1, qa2, beta, ga1, ga2) {
  qa1 <- toupper(qa1); qa2 <- toupper(qa2)
  ga1 <- toupper(ga1); ga2 <- toupper(ga2)
  fifelse(qa1 == ga1 & qa2 == ga2, beta,
          fifelse(qa1 == ga2 & qa2 == ga1, -beta, NA_real_))
}

run_mr <- function(outcome, label) {
  m <- merge(apoe, outcome[, .(SNP, A1, A2, b, se)], by = "SNP")
  m[, beta_harm := harm(effect_allele.exposure, other_allele.exposure,
                        beta.exposure, A1, A2)]
  m <- m[!is.na(beta_harm) & se.exposure > 0 & se > 0]
  cat(sprintf("\n%s: %d instruments\n", label, nrow(m)))
  if (nrow(m) < 4) return(NULL)
  obj <- mr_input(bx = m$beta_harm, bxse = m$se.exposure,
                  by = m$b, byse = m$se, snps = m$SNP)
  ivw <- mr_ivw(obj)
  egger <- tryCatch(mr_egger(obj), error = function(e) NULL)
  rows <- list(data.table(Method = paste0(label, " IVW"), N_SNP = nrow(m),
                          Beta = ivw@Estimate, SE = ivw@StdError, P = ivw@Pvalue))
  if (!is.null(egger)) {
    rows[[length(rows) + 1]] <- data.table(
      Method = paste0(label, " Egger"), N_SNP = nrow(m),
      Beta = egger@Estimate, SE = egger@StdError.Est, P = egger@Pvalue.Est,
      Egger_Intercept = egger@Intercept, Egger_InterceptSE = egger@StdError.Int,
      Egger_InterceptP = egger@Pvalue.Int)
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

res <- rbindlist(list(run_mr(amd, "APOE_clumped_AMD"),
                      run_mr(ad, "APOE_clumped_AD")), fill = TRUE)
write.csv(res, file.path(OUT, "apoe_safety_mr_fixed.csv"), row.names = FALSE)
print(res)
cat("\nDONE\n")
