#!/usr/bin/env Rscript
# Run MR-PRESSO on the clumped APOE instruments.

library(data.table)
library(MRPRESSO)

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

run_presso <- function(outcome, label) {
  m <- merge(apoe, outcome[, .(SNP, A1, A2, b, se)], by = "SNP")
  m[, beta_harm := harm(effect_allele.exposure, other_allele.exposure,
                        beta.exposure, A1, A2)]
  m <- m[!is.na(beta_harm) & se.exposure > 0 & se > 0]
  cat(sprintf("\n%s: %d instruments\n", label, nrow(m)))
  dat <- data.frame(beta.outcome = m$b, beta.exposure = m$beta_harm,
                    sd.outcome = m$se, sd.exposure = m$se.exposure)
  res <- mr_presso(BetaOutcome = "beta.outcome",
                   BetaExposure = "beta.exposure",
                   SdOutcome = "sd.outcome",
                   SdExposure = "sd.exposure",
                   data = dat,
                   OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
                   NbDistribution = 1000, SignifThreshold = 0.05,
              seed = 20260910)
  print(res[["Main MR results"]])
  print(res[["MR-PRESSO results"]])
  data.table::fwrite(as.data.frame(res[["Main MR results"]]),
                     file.path(OUT, paste0("mr_presso_", label, ".csv")),
                     row.names = FALSE)
}

run_presso(amd, "AMD")
run_presso(ad, "AD")
cat("\nDONE\n")
