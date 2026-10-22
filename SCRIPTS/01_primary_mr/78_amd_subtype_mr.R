#!/usr/bin/env Rscript
##############################################################################
# 78_amd_subtype_mr.R (2026-09-09)
#
# 干性/湿性 AMD 亚型分层 MR（PI 2026-09-09 要求，补充材料）
#   - 暴露：UKB-PPP clumped APOE cis-pQTL 工具（与主分析同一套 11 个）
#   - 结局：FinnGen R13 DRY_AMD（干性，含地图状萎缩）/ WET_AMD（湿性，新生血管）
#           + H7_AMD（全体，作为管线锚点：必须复现 IVW β=+0.181777146570588）
#   - 另出 ε4(rs429358)/ε2(rs7412) 的 per-allele（alt 等位）效应，用于镜像检查
# 产出：
#   revision/results/amd_subtype/amd_subtype_mr.csv
#   revision/results/amd_subtype/amd_subtype_allele_effects.csv
#   revision/results/amd_subtype/amd_subtype_meta.csv  (R13 + MVP 亚型汇总)
##############################################################################
suppressMessages({library(data.table); library(MendelianRandomization)})
source("/home/huan/projects/08_APOE_AMD_AD_MR/revision/mr_lib.R")

PROJ <- "/home/huan/projects/08_APOE_AMD_AD_MR"
OUTD <- file.path(PROJ, "revision", "results", "amd_subtype")
dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)

## ---------- 1. 暴露：UKB-PPP clumped APOE 工具（与 anchor_check.sh 同源） ----------
pq <- fread("/home/huan/database/pQTL/pqtl/ukb-ppp/summary/ukb-ppp-2023-5e-08-clumped-10000kb-0.001r2.csv")
iv <- pq[gene.exposure == "APOE",
         .(SNP, qtl_effect = effect_allele.exposure, qtl_other = other_allele.exposure,
           beta = beta.exposure, se = se.exposure, p = pval.exposure)]
cat(sprintf("[0] APOE instruments: %d\n", nrow(iv)))

## ---------- 2. 结局：FinnGen R13 亚型文件 ----------
FF <- "/home/huan/database/GWAS/AMD"
files <- c(DRY = "finngen_R13_DRY_AMD.gz",
           WET = "finngen_R13_WET_AMD.gz",
           H7  = "finngen_R13_H7_AMD.gz")

read_finngen_subset <- function(path, snps) {
  d <- fread(cmd = sprintf("zcat %s | awk -F'\\t' 'NR==1 || $5 ~ /^(%s)$/'",
                           path, paste(snps, collapse = "|")), header = TRUE)
  setnames(d, gsub("^#", "", names(d)))
  # FinnGen 约定：effect allele = alt（beta 按 alt 编码）
  d[, .(SNP = rsids, A1 = toupper(alt), A2 = toupper(ref),
        b = as.numeric(beta), se = as.numeric(sebeta), p = as.numeric(pval))]
}

TAGS <- c("rs429358", "rs7412")           # ε4 / ε2 编码变异
res_mr <- list(); res_allele <- list()

for (nm in names(files)) {
  path <- file.path(FF, files[[nm]])
  oc <- read_finngen_subset(path, c(iv$SNP, TAGS))
  oc <- unique(oc, by = "SNP")
  cat(sprintf("[1] %s: %d rows matched (%d instruments + tags)\n", nm, nrow(oc), nrow(oc)))

  ## 2a. 工具 MR（IVW / Egger / WMed / Wald(best)）
  r <- run_mr_pair(iv, oc[SNP %in% iv$SNP], paste0("APOE_", nm))
  r[, cohort := "FinnGen R13"]
  r[, subtype := nm]
  res_mr[[nm]] <- r

  ## 2b. ε4 / ε2 per-allele（alt 等位）——镜像检查
  for (tg in TAGS) {
    row <- oc[SNP == tg]
    if (nrow(row) != 1) { cat(sprintf("    WARN: %s missing in %s\n", tg, nm)); next }
    res_allele[[paste(nm, tg)]] <- data.table(
      cohort = "FinnGen R13", subtype = nm, SNP = tg,
      effect_allele = row$A1, other_allele = row$A2,
      beta = row$b, se = row$se, p = row$p,
      OR = exp(row$b),
      OR_lo = exp(row$b - 1.96 * row$se), OR_hi = exp(row$b + 1.96 * row$se))
  }
}

mr <- rbindlist(res_mr, fill = TRUE)
al <- rbindlist(res_allele, fill = TRUE)

## ---------- 3. 锚点校验：H7 必须复现主分析 IVW ----------
g <- mr[label == "APOE_H7" & Method == "IVW"]
stopifnot(nrow(g) == 1)
cat(sprintf("\n[2] ANCHOR H7: IVW beta=%.12f  p=%.6g  (expected +0.181777146570588 / 7.89243e-11)\n",
            g$Beta, g$P))
if (abs(g$Beta - 0.181777146570588) > 1e-6) {
  stop(sprintf("ANCHOR FAIL: H7 IVW beta=%.12f", g$Beta), call. = FALSE)
}
cat("    ANCHOR PASS\n")

## ---------- 4. MVP 亚型（已在 74 脚本中跑出，此处只并入汇总表） ----------
mvp_path <- file.path(PROJ, "revision", "results", "mvp_ancestry", "mvp_mr_all.csv")
mvp <- fread(mvp_path)
mvp_sub <- mvp[cohort %in% c("dry_META", "wet_META", "META") &
               method %in% c("IVW", "Egger", "WMed", "PRESSO")]
if (nrow(mvp_sub)) {
  mvp_sub[, `:=`(cohort_full = "Million Veteran Program (META)",
                 subtype = fifelse(cohort == "dry_META", "DRY",
                                   fifelse(cohort == "wet_META", "WET", "ALL")))]
  setnames(mvp_sub, c("beta", "se", "p"), c("Beta", "SE", "P"))
  mvp_sub[, Method := fifelse(method == "WMed", "WMed", method)]
}
r13_out <- mr[, .(cohort, subtype, Method, Beta, SE, P, N_matched, F_mean)]
r13_out[, exposure := "UKB-PPP cis-pQTL"]
mvp_out <- mvp_sub[, .(cohort = cohort_full, subtype, Method, Beta, SE, P,
                       N_matched = n, F_mean = NA_real_, exposure = "UKB-PPP cis-pQTL")]
meta <- rbind(r13_out, mvp_out, fill = TRUE)

fwrite(mr, file.path(OUTD, "amd_subtype_mr.csv"))
fwrite(al, file.path(OUTD, "amd_subtype_allele_effects.csv"))
fwrite(meta, file.path(OUTD, "amd_subtype_meta.csv"))

cat("\n[3] FinnGen R13 subtype MR (IVW):\n")
print(mr[Method == "IVW", .(subtype, Beta = round(Beta, 4), SE = round(SE, 4),
                            P = signif(P, 4), N_matched, F_mean = round(F_mean, 1))])
cat("\n[4] per-allele mirror (FinnGen R13):\n")
print(al[, .(subtype, SNP, effect_allele, beta = round(beta, 4), se = round(se, 4),
             P = signif(p, 4), OR = round(OR, 3))])
cat("\nDONE ->", OUTD, "\n")
