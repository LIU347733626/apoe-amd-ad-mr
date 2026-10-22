#!/usr/bin/env Rscript
# Independent AMD replication using IAMDGC late AMD summary statistics.

library(data.table)
library(gsmr)
library(MendelianRandomization)
library(MRPRESSO)

DB <- "/home/huan/database/pQTL/pqtl/ukb-ppp/summary"
IAMDGC <- "/home/huan/database/GWAS/AMD/IAMDGC_late.ma"
BFILE <- "/home/huan/references/1000g/EU1000QC"
PLINK <- "/home/huan/software/configs_G19/bin/plink"
OUT <- "/home/huan/projects/08_APOE_AMD_AD_MR/revision/results"
TMP <- "/home/huan/projects/08_APOE_AMD_AD_MR/revision/tmp_iamdgc"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
dir.create(TMP, showWarnings = FALSE, recursive = TRUE)

pq <- fread(file.path(DB, "ukb-ppp-2023-5e-08-clumped-10000kb-0.001r2.csv"))
apoe <- pq[gene.exposure == "APOE"]
cat(sprintf("APOE clumped instruments: %d\n", nrow(apoe)))

snp_file <- file.path(OUT, "apoe_iamdgc_snps.txt")
writeLines(apoe$SNP, snp_file)

iam <- fread(cmd = sprintf("grep -w -F -f %s %s", snp_file, IAMDGC),
             header = FALSE,   # FIX(2026-09-10): grep 只输出匹配行；IAMDGC 的表头行
                               # 不匹配 rsID，故输出中本就没有表头。原写 header=TRUE
                               # 会把第一行数据当作列名，静默丢掉 1 个工具（rs28399649）。
             colClasses = c("character", "character", "character", "numeric",
                            "numeric", "numeric", "numeric", "integer"))
setnames(iam, c("SNP", "A1", "A2", "freq", "b", "se", "p", "n"))
iam <- iam[!is.na(b) & !is.na(se) & se > 0 & !is.na(p)]
# 断言：每个工具都必须能在结局文件中找到，否则报出缺失（防再次静默丢工具）
miss_tool <- setdiff(apoe$SNP, unique(iam$SNP))
if (length(miss_tool) > 0) {
  cat(sprintf("WARNING: %d/%d instruments absent from IAMDGC: %s\n",
              length(miss_tool), nrow(apoe), paste(miss_tool, collapse = ", ")))
}
stopifnot(nrow(iam) >= nrow(apoe))
cat(sprintf("Matched IAMDGC SNPs: %d (instruments covered: %d/%d)\n",
            nrow(iam), sum(apoe$SNP %in% iam$SNP), nrow(apoe)))

m <- merge(apoe, iam, by = "SNP")
m[, beta_harm := fifelse(
  toupper(effect_allele.exposure) == toupper(A1) &
    toupper(other_allele.exposure) == toupper(A2),
  beta.exposure,
  fifelse(toupper(effect_allele.exposure) == toupper(A2) &
            toupper(other_allele.exposure) == toupper(A1),
          -beta.exposure, NA_real_))]
m <- m[!is.na(beta_harm) & !is.na(se.exposure) & !is.na(se) &
         se.exposure > 0 & se > 0 & !duplicated(SNP)]
cat(sprintf("Harmonized SNPs: %d\n", nrow(m)))

# GSMR with 1000G LD
snp_f <- file.path(TMP, "snps.txt")
ld_f <- file.path(TMP, "ld")
writeLines(m$SNP, snp_f)
system(sprintf("%s --bfile %s --extract %s --r square --write-snplist --out %s --threads 2 --silent",
               PLINK, BFILE, snp_f, ld_f))
ld_snps <- readLines(paste0(ld_f, ".snplist"))
ld <- as.matrix(read.table(paste0(ld_f, ".ld")))
colnames(ld) <- rownames(ld) <- ld_snps
m <- m[match(ld_snps, m$SNP), ]
m <- m[!is.na(m$SNP), ]

gsmr_res <- gsmr(m$beta_harm, m$se.exposure, m$pval.exposure,
                 m$b, m$se, ld, ld_snps, n_ref = 503,
                 gwas_thresh = 5e-8, nsnps_thresh = 2,
                 heidi_outlier_flag = TRUE)
gsmr_out <- data.table(
  Replication = "IAMDGC_late_AMD",
  Gene = "APOE",
  N_SNP = nrow(m),
  N_used = if (!is.null(gsmr_res$used_index)) length(gsmr_res$used_index) else NA,
  Beta = if (length(gsmr_res$bxy) > 0) gsmr_res$bxy[1] else NA,
  SE = if (length(gsmr_res$bxy_se) > 0) gsmr_res$bxy_se[1] else NA,
  P = if (length(gsmr_res$bxy_pval) > 0) gsmr_res$bxy_pval[1] else NA
)

# MR sensitivity
obj <- mr_input(bx = m$beta_harm, bxse = m$se.exposure,
                by = m$b, byse = m$se, snps = m$SNP)
ivw <- mr_ivw(obj)
egger <- tryCatch(mr_egger(obj), error = function(e) NULL)
median <- tryCatch(mr_median(obj), error = function(e) NULL)
dat <- data.frame(beta.outcome = m$b, beta.exposure = m$beta_harm,
                  sd.outcome = m$se, sd.exposure = m$se.exposure)
presso <- tryCatch(
  mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
            SdOutcome = "sd.outcome", SdExposure = "sd.exposure",
            data = dat, OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
            NbDistribution = 1000, SignifThreshold = 0.05,
            seed = 20260910),   # FIX(2026-09-10): 未固定种子时 MR-PRESSO 每次
                                # 结果不同（实测同数据 5 次给 0.122–0.157），必须固定
  error = function(e) NULL
)
mr_out <- rbindlist(list(
  data.table(Method = "IAMDGC IVW", N_SNP = nrow(m),
             Beta = ivw@Estimate, SE = ivw@StdError, P = ivw@Pvalue),
  if (!is.null(egger)) data.table(Method = "IAMDGC Egger", N_SNP = nrow(m),
                                   Beta = egger@Estimate, SE = egger@StdError.Est,
                                   P = egger@Pvalue.Est,
                                   Intercept = egger@Intercept,
                                   Intercept_P = egger@Pvalue.Int),
  if (!is.null(median)) data.table(Method = "IAMDGC Weighted_median", N_SNP = nrow(m),
                                    Beta = median@Estimate, SE = median@StdError,
                                    P = median@Pvalue),
  if (!is.null(presso)) data.table(
    Method = "IAMDGC MR-PRESSO_raw", N_SNP = nrow(m),
    Beta = presso[["Main MR results"]][1, "Causal Estimate"],
    SE = presso[["Main MR results"]][1, "Sd"],
    P = presso[["Main MR results"]][1, "P-value"])
), fill = TRUE)

write.csv(gsmr_out, file.path(OUT, "replication_iamdgc_gsmr.csv"), row.names = FALSE)
write.csv(mr_out, file.path(OUT, "replication_iamdgc_mr.csv"), row.names = FALSE)
print(gsmr_out)
print(mr_out)
cat("\nDONE\n")
