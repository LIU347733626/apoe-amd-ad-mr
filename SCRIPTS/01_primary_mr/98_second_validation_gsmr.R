#!/usr/bin/env Rscript
# 98_second_validation_gsmr.R — 为 Figure 1 的 Second validation 补跑 GSMR（2026-09-10）
#
# 背景：PI 要求 Figure 1 全部行统一为 GSMR。实测发现 Discovery / First validation 已是 GSMR，
#       但 Second validation 两行都不是：AD（Bellenguez）用的是 IVW −0.7627，
#       AMD（MVP）用的是 META,IVW 0.0437 —— 且两个队列**都没有 GSMR 结果**
#       （74_mvp_ancestry_mr.R 的方法注释只列 IVW/Egger/WMed/PRESSO/LOO；
#        全库 grep Bellenguez 的 GSMR 零命中，CLAUDE.md 记的"GSMR −0.764"系误记）。
#
# 本脚本：用与主分析同一套 11 个 clumped UKB-PPP APOE cis-pQTL 工具 + EU1000QC LD，
#         对 Bellenguez 2022（AD）与 MVP META（AMD）各跑一次 GSMR（gsmr 1.0.6 的 7 参数 API）。
#
# 输出：revision/results/r13/second_validation_gsmr.csv
suppressPackageStartupMessages({library(data.table); library(gsmr)})
PROJ <- "/home/huan/projects/08_APOE_AMD_AD_MR"
OUT  <- file.path(PROJ, "revision", "results", "r13")
TMP  <- file.path(PROJ, "tmp", "sv_gsmr"); dir.create(TMP, recursive = TRUE, showWarnings = FALSE)
PQ   <- "/home/huan/database/pQTL/pqtl/ukb-ppp/summary/ukb-ppp-2023-5e-08-clumped-10000kb-0.001r2.csv"
BFILE<- "/home/huan/references/1000g/EU1000QC"; PLINK <- "/home/huan/software/configs_G19/bin/plink"
BELL <- file.path(PROJ, "smr_data/AD_Bellenguez2022.ma")          # 无表头，8 列
MVP  <- "/mnt/ssd4t/GWAS/AMD/MvP/extract/MVP_R4.1000G_AGR.Phe_362_2.META.GIA.dbGaP.txt.gz"

pq <- fread(PQ); iv <- pq[gene.exposure == "APOE"]
cat(sprintf("instruments: %d (rs429358 present: %s)\n", nrow(iv), "rs429358" %in% iv$SNP))
sf <- file.path(TMP, "snps.txt"); writeLines(iv$SNP, sf)

## ---- 两个结局的读取（各自格式不同，均只取工具位点）----
read_bell <- function() {
  d <- fread(cmd = sprintf("grep -w -F -f %s %s", sf, BELL), header = FALSE)
  setnames(d, c("SNP","A1","A2","freq","b","se","p","n"))
  d[, .(SNP, A1, A2, freq, b = as.numeric(b), se = as.numeric(se))]
}
MVP_COLS <- c("SNP_ID","chrom","pos","ref","alt","ea","af","num_samples","case_af",
              "num_cases","control_af","num_controls","or","ci","pval","r2","q_pval","i2","direction")
read_mvp <- function() {
  z <- fread(cmd = sprintf("zcat %s | grep -w -F -f %s", MVP, sf), header = FALSE, sep = "\t")
  setnames(z, MVP_COLS[seq_len(ncol(z))])
  z <- z[, .(SNP = SNP_ID, A1 = ea, A2 = ref, freq = as.numeric(af),
             b = log(as.numeric(or)),
             se = (log(as.numeric(sub(".*,", "", ci))) - log(as.numeric(sub(",.*", "", ci)))) / 3.92)]
  z[is.finite(b) & is.finite(se) & se > 0]
}

## ---- LD（1000G EUR，与主分析同一参考）----
prep <- function(ev) {
  m <- merge(iv, ev, by = "SNP")
  m[, beta_h := fifelse(toupper(effect_allele.exposure) == toupper(A1) &
                          toupper(other_allele.exposure) == toupper(A2), beta.exposure,
                fifelse(toupper(effect_allele.exposure) == toupper(A2) &
                          toupper(other_allele.exposure) == toupper(A1), -beta.exposure, NA_real_))]
  m <- m[!is.na(beta_h) & !duplicated(SNP)]
  f <- file.path(TMP, "ld_snps.txt"); writeLines(m$SNP, f)
  pfx <- file.path(TMP, "ld")
  system(sprintf("%s --bfile %s --extract %s --r square --write-snplist --out %s --silent",
                 PLINK, BFILE, f, pfx))
  ids <- readLines(paste0(pfx, ".snplist"))
  ld <- as.matrix(read.table(paste0(pfx, ".ld"))); colnames(ld) <- rownames(ld) <- ids
  m <- m[match(ids, m$SNP), ]
  list(m = m, ld = ld, ids = ids)
}

run_gsmr <- function(tag, ev, n_ref = 503) {
  p <- prep(ev)
  cat(sprintf("\n=== %s: harmonized %d, in-LD %d ===\n", tag, nrow(p$m), length(p$ids)))
  r <- tryCatch(
    gsmr(p$m$beta_h, p$m$se.exposure, p$m$pval.exposure,
         p$m$b, p$m$se, p$ld, p$ids,
         n_ref = n_ref, gwas_thresh = 5e-8, nsnps_thresh = 2, heidi_outlier_flag = TRUE),
    error = function(e) { cat("  GSMR ERROR:", conditionMessage(e), "\n"); NULL })
  if (is.null(r)) return(NULL)
  dt <- data.table(
    cohort = tag, method = "GSMR",
    n_snp = nrow(p$m), n_used = if (!is.null(r$used_index)) length(r$used_index) else NA_integer_,
    beta = if (length(r$bxy)) r$bxy[1] else NA_real_,
    se   = if (length(r$bxy_se)) r$bxy_se[1] else NA_real_,
    p    = if (length(r$bxy_pval)) r$bxy_pval[1] else NA_real_)
  print(dt)
  dt
}

res <- rbindlist(list(run_gsmr("Bellenguez_2022_AD", read_bell()),
                      run_gsmr("MVP_META_AMD",     read_mvp())), fill = TRUE)
fwrite(res, file.path(OUT, "second_validation_gsmr.csv"))
cat("\n=== 汇总 ===\n"); print(res)
cat("\n对照（图1 现值）：Bellenguez IVW lnOR=-0.76266 ; MVP META IVW lnOR=+0.04372\n")
cat("DONE\n")
