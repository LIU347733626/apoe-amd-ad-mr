#!/usr/bin/env Rscript
# 59b: CLUMPED pQTL instrument set (11, epsilon4-tag excluded by design)
#      x {AMD R13, AD R13, AD Jansen 2019, AMD IAMDGC} colocalization
#
# 目的：复现并扩展手稿 "colocalization asymmetry" 叙事的可复现证据链
#   —— clumped（ε4-independent）水平仪器集：
#        AMD 侧预期 PP4（AMD 信号就在水平变异上）
#        AD  侧预期 PP3（AD 信号在 ε4 簇，不在仪器集内）
#   —— 与 59 的 unclumped（1407 SNP，含 ε4 tag）run 对照：
#        unclumped 下 ε4 tag 本身是 APOE 最强水平变异 → AD 也 PP4 高
#   —— 两种仪器集的说法都报，手稿以 clumped 为 asymmetry 主声明，
#      unclumped 作为 "含 ε4 时 AD 亦共定位" 的补充句（level 模型自洽）。
#
# s 值溯源：
#   AD R13   = 16652/500186（23b_coloc_ad_r13.R 已用值，有盘上溯源）
#   AMD R13  = 5500/472750（23_r13.R 已用值，R8 时代沿用 —— 不精确但为历史
#               run 的原始参数；另跑 s=16000/N 敏感性行，若 PP4 不变即稳健）
#   Jansen   = 71880/436498（Jansen 2019 临床+proxy 病例数，前期会话验证）
#   IAMDGC   = 16144/32339（Fritsche 2016 late AMD 病例数，kp4cd/早期会话验证）
#
# 输出：revision/results/jansen/clumped_coloc_{tag}.csv
##############################################################################
suppressPackageStartupMessages({
  library(data.table); library(coloc)
})

PROJ <- "/home/huan/projects/08_APOE_AMD_AD_MR"
IV   <- file.path(PROJ, "revision/results/r13/apoe_pqtl_instruments_r13.csv")
OUT  <- file.path(PROJ, "revision/results/jansen")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

cat("===== 59b: CLUMPED (11-instrument, eps4-excluded) pQTL x 4 cohorts =====\n")

# D1: clumped 11 instruments ----------------------------------------------
iv <- fread(IV)
iv <- iv[order(pos.exposure)]
cat(sprintf("Clumped instruments: %d (range %.4f-%.4f Mb)\n", nrow(iv),
            min(iv$pos.exposure)/1e6, max(iv$pos.exposure)/1e6))

run_coloc <- function(tag, gwas_file, use_header, cases, n_total, note = "") {
  cat(sprintf("\n##### %s (N=%d, s=%d/%d=%.4f)%s\n", tag, n_total, cases, n_total,
              cases/n_total, if (nzchar(note)) sprintf(" [%s]", note) else ""))
  g <- fread(gwas_file, header = use_header)
  setnames(g, c("SNP", "A1", "A2", "freq", "b", "se", "p", "n"))
  g <- g[!is.na(SNP) & SNP %in% iv$SNP & !is.na(b) & !is.na(se) & se > 0]
  cat(sprintf("Matched outcome SNPs: %d / %d instruments\n", nrow(g), nrow(iv)))
  if (nrow(g) < 5) { cat("TOO FEW - skip\n"); return(NULL) }
  m <- merge(iv, g, by = "SNP")
  m[, beta_harm := fifelse(toupper(effect_allele.exposure) == toupper(A1) &
                             toupper(other_allele.exposure) == toupper(A2),
                           beta.exposure,
                           fifelse(toupper(effect_allele.exposure) == toupper(A2) &
                                     toupper(other_allele.exposure) == toupper(A1),
                                   -beta.exposure, NA_real_))]
  m <- m[!is.na(beta_harm) & se.exposure > 0 & se > 0 &
           !is.na(freq) & freq > 0 & freq < 1 &
           !is.na(eaf.exposure) & eaf.exposure > 0 & eaf.exposure < 1 & !duplicated(SNP)]
  cat(sprintf("Harmonized instruments: %d\n", nrow(m)))
  if (nrow(m) < 5) { cat("TOO FEW after harmonization - skip\n"); return(NULL) }

  p_out <- 2 * pnorm(-abs(m$b / m$se))
  p_qtl <- 2 * pnorm(-abs(m$beta.exposure / m$se.exposure))
  D1 <- list(type = "quant", N = unique(m$samplesize.exposure)[1], sdY = 1,
             snp = m$SNP, pvalues = p_qtl, beta = m$beta.exposure,
             varbeta = m$se.exposure^2,
             MAF = ifelse(m$eaf.exposure <= 0.5, m$eaf.exposure, 1 - m$eaf.exposure),
             position = m$pos.exposure)
  D2 <- list(type = "cc", N = n_total, s = cases / n_total,
             snp = m$SNP, pvalues = p_out, beta = m$b, varbeta = m$se^2,
             MAF = ifelse(m$freq <= 0.5, m$freq, 1 - m$freq),
             position = m$pos.exposure)
  fit <- coloc.abf(D1, D2)
  s <- fit$summary
  out <- data.table(Module = tag, Gene = "APOE", D1_set = "clumped11",
                    N_SNP = nrow(m), s_outcome = cases / n_total,
                    PP.H0 = s["PP.H0.abf"], PP.H1 = s["PP.H1.abf"],
                    PP.H2 = s["PP.H2.abf"], PP.H3 = s["PP.H3.abf"],
                    PP.H4 = s["PP.H4.abf"], note = note)
  fwrite(out, file.path(OUT, paste0("clumped_coloc_", tag, ".csv")), col.names = TRUE)
  print(out)
  out
}

r_amd  <- run_coloc("amd_r13",  file.path(PROJ, "smr_data/AMD_r13.ma"),
                    TRUE, 5500, 472750, note = "s from 23_r13 line (5500/N; R8-era value)")
r_amds <- run_coloc("amd_r13_hi", file.path(PROJ, "smr_data/AMD_r13.ma"),
                    TRUE, 16000, 472750, note = "sensitivity: upper plausible case band (~16k/N)")
r_ad   <- run_coloc("ad_r13",   file.path(PROJ, "smr_data/AD_r13.ma"),
                    FALSE, 16652, 500186, note = "s from 23b_coloc_ad_r13 line (16652/N)")
r_jana <- run_coloc("ad_jansen", file.path(PROJ, "smr_data/AD_Jansen2019_converted.ma"),
                    FALSE, 71880, 436498, note = "AD-by-proxy composition; per-allele scale caveat")
r_iam  <- run_coloc("amd_iamdgc", "/home/huan/database/GWAS/AMD/IAMDGC_late.ma",
                    TRUE, 16144, 32339, note = "Fritsche 2016 late AMD; 7 chr1/chr10 flipped variants (non-APOE) caveat")

cat("\n=== Summary (clumped D1) ===\n")
for (r in list(r_amd, r_amds, r_ad, r_jana, r_iam)) {
  if (is.null(r)) next
  cat(sprintf("%-14s: PP3=%.4f  PP4=%.4f  (N_SNP=%d)\n",
              r$Module, r$PP.H3, r$PP.H4, r$N_SNP))
}
cat("\nDone.\n")
