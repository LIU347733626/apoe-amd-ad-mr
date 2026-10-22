#!/usr/bin/env Rscript
##############################################################################
# 74_mvp_ancestry_mr.R (2026-09-09)
#
# 缺口补齐（AMD 第二验证队列 = Million Veteran Program）：
#   分族裔（EUR/AFR/AMR/EAS）+ 跨族裔 META + 干性/湿性亚型 的 APOE 水平 MR。
#   暴露：UKB-PPP 11 个 clumped APOE 工具（与主分析同一套）。
#   方法：IVW / MR-Egger / weighted median / MR-PRESSO / leave-one-out。
#
# 输入：/mnt/ssd4t/GWAS/AMD/MvP/extract/MVP_R4.1000G_AGR.Phe_362_{2,21,22}.{EUR,AFR,AMR,EAS,META}.GIA.dbGaP.txt.gz
# 产出：revision/results/mvp_ancestry/mvp_mr_all.csv
##############################################################################
suppressMessages({library(data.table); library(MendelianRandomization); library(MRPRESSO)})

PROJ <- "/home/huan/projects/08_APOE_AMD_AD_MR"
EXT  <- "/mnt/ssd4t/GWAS/AMD/MvP/extract"
OUTD <- file.path(PROJ, "revision", "results", "mvp_ancestry")
dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)

## ---------- 暴露：UKB-PPP 11 工具 ----------
iv <- fread(file.path(PROJ, "revision/results/r13/supp_S1c_instruments.csv"))[gene == "APOE"]
iv <- unique(iv[, .(SNP, eff = effect_allele, oth = other_allele, beta = beta_exposure, se = se_exposure)])
cat(sprintf("[0] instruments: %d\n", nrow(iv)))
snpfile <- file.path(OUTD, "instruments.txt")
writeLines(iv$SNP, snpfile)

## ---------- 结局：逐个 MVP 文件 ----------
read_mvp <- function(path, tag) {
  pat <- paste(readLines(snpfile), collapse = "|")
  cmd <- sprintf("zcat %s | awk -F'\\t' 'NR==1 || $1 ~ /^(%s)$/'", path, pat)
  d <- tryCatch(fread(cmd = cmd, header = TRUE), error = function(e) NULL)
  if (is.null(d) || !nrow(d)) { cat(sprintf("    %s: 0 matched\n", tag)); return(NULL) }
  # or = odds ratio, ci = "lo,hi" -> beta = ln(or), se = (ln hi - ln lo)/3.92
  orv <- as.numeric(d[["or"]]); civ <- as.character(d[["ci"]])
  lo <- as.numeric(sub(",.*", "", civ)); hi <- as.numeric(sub(".*,", "", civ))
  d[, b  := log(orv)]
  d[, gse := (log(hi) - log(lo)) / 3.92]
  d[, tag := tag]
  d[is.finite(b) & is.finite(gse) & gse > 0, .(SNP = SNP_ID, A1 = ea, A2 = ref, b, gse, tag)]
}

files <- c(
  "EUR"  = "MVP_R4.1000G_AGR.Phe_362_2.EUR.GIA.dbGaP.txt.gz",
  "AFR"  = "MVP_R4.1000G_AGR.Phe_362_2.AFR.GIA.dbGaP.txt.gz",
  "AMR"  = "MVP_R4.1000G_AGR.Phe_362_2.AMR.GIA.dbGaP.txt.gz",
  "EAS"  = "MVP_R4.1000G_AGR.Phe_362_2.EAS.GIA.dbGaP.txt.gz",
  "META" = "MVP_R4.1000G_AGR.Phe_362_2.META.GIA.dbGaP.txt.gz",
  "dry_META" = "MVP_R4.1000G_AGR.Phe_362_21.META.GIA.dbGaP.txt.gz",
  "wet_META" = "MVP_R4.1000G_AGR.Phe_362_22.META.GIA.dbGaP.txt.gz")

harmonize <- function(qa1, qa2, qb, ga1, ga2) {
  qa1 <- toupper(qa1); qa2 <- toupper(qa2); ga1 <- toupper(ga1); ga2 <- toupper(ga2)
  fifelse(qa1 == ga1 & qa2 == ga2, qb,
          fifelse(qa1 == ga2 & qa2 == ga1, -qb, NA_real_))
}

res <- list()
for (tag in names(files)) {
  path <- file.path(EXT, files[[tag]])
  if (!file.exists(path)) { cat(sprintf("[1] %s: FILE MISSING\n", tag)); next }
  oc <- read_mvp(path, tag)
  if (is.null(oc)) next
  m <- merge(iv, oc, by = "SNP")
  if (!nrow(m)) { cat(sprintf("[1] %s: no merged rows\n", tag)); next }
  m[, beta_h := harmonize(eff, oth, beta, A1, A2)]     # A1 = MVP effect allele
  m <- m[!is.na(beta_h) & se > 0 & gse > 0]
  cat(sprintf("[1] %s: matched %d/%d instruments\n", tag, nrow(m), nrow(iv)))
  if (nrow(m) < 3) {
    res[[tag]] <- data.table(cohort = tag, method = "TOO_FEW", n = nrow(m), beta = NA_real_, se = NA_real_, p = NA_real_)
    next
  }
  obj <- mr_input(bx = m$beta_h, bxse = m$se, by = m$b, byse = m$gse, snps = m$SNP)
  ivw <- tryCatch(mr_ivw(obj), error = function(e) NULL)
  egr <- tryCatch(mr_egger(obj), error = function(e) NULL)
  wm  <- tryCatch(mr_median(obj), error = function(e) NULL)
  pr  <- tryCatch(mr_presso(BetaOutcome = "b", BetaExposure = "beta_h", SdOutcome = "gse", SdExposure = "se",
                            data = as.data.frame(m), OUTLIERtest = TRUE, DISTORTIONtest = TRUE, NbDistribution = 1000, SignifThreshold = 0.05,
              seed = 20260910),
                  error = function(e) NULL)
  rows <- list(
    data.table(cohort = tag, method = "IVW", n = nrow(m), beta = ivw@Estimate, se = ivw@StdError, p = ivw@Pvalue),
    data.table(cohort = tag, method = "Egger", n = nrow(m), beta = egr@Estimate, se = egr@StdError.Est, p = egr@Pvalue.Est),
    data.table(cohort = tag, method = "WMed", n = nrow(m), beta = wm@Estimate, se = wm@StdError, p = wm@Pvalue))
  if (!is.null(pr)) {
    g <- pr$`Main MR results`
    rows[[length(rows) + 1]] <- data.table(cohort = tag, method = "PRESSO",
      n = nrow(m), beta = g$`Causal Estimate`[1], se = g$`Sd`[1], p = g$`P-value`[1])
  }
  res[[tag]] <- rbindlist(rows)
  # leave-one-out (IVW)
  loo <- rbindlist(lapply(seq_len(nrow(m)), function(i) {
    mm <- m[-i]
    o <- mr_input(bx = mm$beta_h, bxse = mm$se, by = mm$b, byse = mm$gse, snps = mm$SNP)
    r <- tryCatch(mr_ivw(o), error = function(e) NULL)
    if (is.null(r)) return(NULL)
    data.table(cohort = tag, method = paste0("LOO_drop_", m$SNP[i]), n = nrow(mm),
               beta = r@Estimate, se = r@StdError, p = r@Pvalue)
  }), fill = TRUE)
  if (nrow(loo)) res[[paste0(tag, "_loo")]] <- loo
}
out <- rbindlist(res, fill = TRUE)
fwrite(out, file.path(OUTD, "mvp_mr_all.csv"))
cat(sprintf("[2] wrote %d rows -> %s\n", nrow(out), file.path(OUTD, "mvp_mr_all.csv")))
print(out[!grepl("^LOO", method)])
cat("DONE\n")
