#!/usr/bin/env Rscript
# 96_leaveoneout_official11.R — 项目官方 11 工具集的 leave-one-out 与 FE/RE 双口径
# 目的：直接回答"AMD 结论是否由单个 SNP 驱动"与"固定/随机效应下显著性如何"
# 数据：官方 clumped 文件（项目主分析所用）：ukb-ppp-2023-5e-08-clumped-10000kb-0.001r2.csv
# 输出：revision/results/r13/loo_official11_{AMD_R13,AD_R13,IAMDGC}.csv
#       revision/results/r13/loo_official11_summary.csv
suppressPackageStartupMessages({library(data.table); library(MendelianRandomization)})
PROJ <- "/home/huan/projects/08_APOE_AMD_AD_MR"
OUT  <- file.path(PROJ, "revision", "results", "r13")
PQ   <- "/home/huan/database/pQTL/pqtl/ukb-ppp/summary/ukb-ppp-2023-5e-08-clumped-10000kb-0.001r2.csv"
ENDP <- list(
  AMD_R13   = list(f = file.path(PROJ, "smr_data/AMD_r13.ma"),    hdr = TRUE),
  AD_R13    = list(f = file.path(PROJ, "smr_data/AD_r13_hdr.ma"), hdr = TRUE),
  IAMDGC    = list(f = "/home/huan/database/GWAS/AMD/IAMDGC_late.ma", hdr = TRUE)
)

pq <- fread(PQ)
apoe <- pq[gene.exposure == "APOE"]
cat(sprintf("Official clumped APOE instruments: %d\n", nrow(apoe)))
cat("  ", paste(apoe$SNP, collapse = ", "), "\n")
cat(sprintf("  rs429358 present in instrument set? %s\n", "rs429358" %in% apoe$SNP))

sf <- file.path(OUT, "loo11_snps.txt"); writeLines(apoe$SNP, sf)

read_out <- function(f) {
  d <- fread(cmd = sprintf("grep -w -F -f %s %s", sf, f), header = FALSE)  # grep 输出恒无表头
  setnames(d, c("SNP","A1","A2","freq","b","se","p","n"))
  d[, p := as.numeric(p)][!is.na(b) & !is.na(se) & se > 0]
}

summ <- list()
for (nm in names(ENDP)) {
  o <- read_out(ENDP[[nm]]$f)
  m <- merge(apoe, o, by = "SNP")
  m[, beta_h := fifelse(toupper(effect_allele.exposure) == toupper(A1) &
                          toupper(other_allele.exposure) == toupper(A2), beta.exposure,
                fifelse(toupper(effect_allele.exposure) == toupper(A2) &
                          toupper(other_allele.exposure) == toupper(A1), -beta.exposure, NA_real_))]
  m <- m[!is.na(beta_h) & !duplicated(SNP)]
  cat(sprintf("\n### %s: harmonized k=%d\n", nm, nrow(m)))
  rows <- list()
  add <- function(tag, bx, bxse, by, byse, snps) {
    obj <- mr_input(bx = bx, bxse = bxse, by = by, byse = byse, snps = snps)
    for (mdl in c("fixed", "random")) {
      r <- tryCatch(mr_ivw(obj, model = mdl), error = function(e) NULL)
      if (is.null(r)) next
      rows[[length(rows) + 1]] <<- data.table(
        outcome = nm, excluded = tag, model = mdl, k = length(bx),
        beta = r@Estimate, se = r@StdError, p = r@Pvalue,
        Q = tryCatch(r@Heter.Stat[1], error = function(e) NA_real_),
        Q_p = tryCatch(r@Heter.Stat[2], error = function(e) NA_real_))
    }
  }
  add("none (full set)", m$beta_h, m$se.exposure, m$b, m$se, m$SNP)
  for (i in seq_len(nrow(m))) {
    s <- m[-i]
    add(paste0("minus ", m$SNP[i]), s$beta_h, s$se.exposure, s$b, s$se, s$SNP)
  }
  r <- rbindlist(rows)
  fwrite(r, file.path(OUT, paste0("loo_official11_", nm, ".csv")))
  full <- r[excluded == "none (full set)"]
  cat("  full set : "); print(full[, .(model, k, beta = round(beta,4), se = round(se,4), p = signif(p,3))])
  cat("  LOO range (fixed): beta", sprintf("%.4f ~ %.4f",
      min(r[excluded != "none (full set)" & model == "fixed"]$beta),
      max(r[excluded != "none (full set)" & model == "fixed"]$beta)),
      "| p max", signif(max(r[excluded != "none (full set)" & model == "fixed"]$p), 3), "\n")
  cat("  LOO range (random): beta", sprintf("%.4f ~ %.4f",
      min(r[excluded != "none (full set)" & model == "random"]$beta),
      max(r[excluded != "none (full set)" & model == "random"]$beta)),
      "| p max", signif(max(r[excluded != "none (full set)" & model == "random"]$p), 3), "\n")
  summ[[nm]] <- r
}
S <- rbindlist(summ)
fwrite(S, file.path(OUT, "loo_official11_summary.csv"))
cat("\n=== 汇总（full set, 每结局每口径一行）===\n")
print(S[excluded == "none (full set)", .(outcome, model, k, beta = round(beta,4),
                                          se = round(se,4), p = signif(p,3), Q_p = signif(Q_p,3))])
cat("\nDONE\n")
