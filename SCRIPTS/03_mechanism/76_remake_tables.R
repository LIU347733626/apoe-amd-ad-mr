#!/usr/bin/env Rscript
##############################################################################
# 76_remake_tables.R (2026-09-09)
#
# 按 v14 队列架构重制补充表：
#   TableS7  robustness grid  -> FinnGen R13 网格 + IAMDGC 第一验证网格
#   TableS16 cross-cohort PIP -> 三队列（R13 AD / R13 AMD / IAMDGC late），去掉 R8
#   TableS21 MVP 分族裔/亚型 MR（新增）
#   TableS22 IAMDGC 网格 + 中介（新增）
##############################################################################
suppressMessages(library(data.table))
PROJ <- "/home/huan/projects/08_APOE_AMD_AD_MR"
OUT  <- file.path(PROJ, "NC_package", "supplementary")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## ---------- Table S7: robustness grids (R13 + IAMDGC) ----------
r13 <- fread(file.path(PROJ, "revision/results/r13/robustness_grid_full.csv"))
r13[, cohort := "FinnGen R13 (discovery)"]
iam <- fread(file.path(PROJ, "revision/results/primary_anchor/robustness_grid_iamdgc.csv"))
iam[, cohort := "IAMDGC late AMD (first validation)"]
cols <- intersect(names(r13), names(iam))
s7 <- rbindlist(list(r13[, ..cols], iam[, ..cols]), fill = TRUE)
setcolorder(s7, c("cohort", setdiff(names(s7), "cohort")))
fwrite(s7, file.path(OUT, "TableS7_robustness_grid.csv"))
cat(sprintf("[S7] %d rows\n", nrow(s7)))

## ---------- Table S16: cross-cohort PIP audit (3 cohorts) ----------
ad_r13 <- fread(file.path(PROJ, "revision/results/r13/finemap_ad_finngen_pip.csv")); setnames(ad_r13, c("SNP","PIP","pos"))
ad_kun <- fread(file.path(PROJ, "revision/results/r13/finemap_ad_kunkle_pip.csv")); setnames(ad_kun, c("SNP","PIP","pos"))
amd_r13 <- fread(file.path(PROJ, "revision/results/finemap_amd_r13_pip.csv")); setnames(amd_r13, c("SNP","PIP","pos"))
iamd <- fread(file.path(PROJ, "revision/results/iamdgc_apoe_finemap_pip.csv"))
setnames(iamd, c("SNP","variable","PIP","pos"), skip_absent = TRUE)
iamd <- unique(iamd[, .(SNP, PIP = as.numeric(PIP), pos = as.numeric(pos))], by = "SNP")

posmap <- unique(rbindlist(list(ad_r13[, .(SNP, pos)], amd_r13[, .(SNP, pos)], ad_kun[, .(SNP, pos)], iamd[, .(SNP, pos)]), fill = TRUE))
posmap <- posmap[!is.na(SNP)][order(SNP)][!duplicated(SNP)]
m <- posmap
for (nm in c("ad_R13", "ad_Kunkle", "amd_R13", "IAMDGC_late")) {
  src <- switch(nm, ad_R13 = ad_r13, ad_Kunkle = ad_kun, amd_R13 = amd_r13, IAMDGC_late = iamd)
  tmp <- unique(src[, .(SNP, pip = as.numeric(PIP))], by = "SNP")
  setnames(tmp, "pip", nm)
  m <- merge(m, tmp, by = "SNP", all.x = TRUE)
}
m[, maxPIP := pmax(ad_R13, ad_Kunkle, amd_R13, IAMDGC_late, na.rm = TRUE)]
m <- m[is.finite(maxPIP) & maxPIP >= 1]
setorder(m, pos)
fwrite(m, file.path(OUT, "TableS16_pip_cross_cohort.csv"))
cat(sprintf("[S16] %d rows with maxPIP = 1\n", nrow(m)))

## ---------- Table S21: MVP per-ancestry / subtype MR ----------
mvp <- fread(file.path(PROJ, "revision/results/mvp_ancestry/mvp_mr_all.csv"))
mvp <- mvp[!grepl("^LOO", method)]
fwrite(mvp, file.path(OUT, "TableS21_mvp_ancestry_mr.csv"))
cat(sprintf("[S21] %d rows\n", nrow(mvp)))

## ---------- Table S22: IAMDGC grid + mediation ----------
med <- fread(file.path(PROJ, "revision/results/primary_anchor/mediation_primary.csv"))
s22 <- rbindlist(list(
  iam[, .(analysis = "robustness grid", cohort, Grid, Outcome, Method, n = N_iv, N_matched, beta = Beta, se = SE, p = P)],
  med[, .(analysis = "mediation", cohort, Grid = NA_character_, Outcome = disease, Method = "product of coefficients",
          n = NA_integer_, N_matched = NA_integer_, beta = indirect, se = indirect_se, p = indirect_p)],
  med[, .(analysis = "mediation (direct)", cohort, Grid = NA_character_, Outcome = disease, Method = "residual direct",
          n = NA_integer_, N_matched = NA_integer_, beta = direct, se = direct_se, p = direct_p)]), fill = TRUE)
fwrite(s22, file.path(OUT, "TableS22_iamdgc_grid_mediation.csv"))
cat(sprintf("[S22] %d rows\n", nrow(s22)))
cat("DONE\n")
