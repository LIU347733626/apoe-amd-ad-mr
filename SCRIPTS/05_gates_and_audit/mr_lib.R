#!/usr/bin/env Rscript
# mr_lib.R — shared MR utilities for 08 project. NEW MR SCRIPTS: source this,
# do NOT re-implement read_ma/harmonize/run_mr (they were the bug farm on 2026-09-01).
#
# CONVENTIONS (violations = wrong signs, see CLAUDE.md 错误档案 2026-09-01 #1/#2):
#  - harmonize(qa1,qa2,qb,ga1,ga2): qa1 MUST be the QTL EFFECT allele (beta is per-effect)
#  - mr_input bx MUST be the harmonized beta (beta_h), never raw
library(data.table)
library(MendelianRandomization)

# 2026-09-09: 线程自动拉满（run_r.sh 会导出 R_DATATABLE_NUM_THREADS / RES_CORES）
# 目的：大文件读取/合并不再只用一半核。手动调用时退回本机核数。
setDTthreads(as.integer(Sys.getenv("R_DATATABLE_NUM_THREADS",
              unset = as.character(max(1L, parallel::detectCores())))))

read_ma <- function(path, header) {
  d <- if (isTRUE(header)) fread(path, header = TRUE) else
    fread(path, header = FALSE,
          colClasses = c("character","character","character","numeric",
                         "numeric","numeric","numeric","integer"))
  setnames(d, c("SNP","A1","A2","freq","b","se","p","n"))
  d[!is.na(b) & !is.na(se) & se > 0 & !is.na(p)]
}

# qa1 = QTL effect allele, qa2 = QTL other allele; ga1 = GWAS A1 (effect), ga2 = GWAS A2
harmonize <- function(qa1, qa2, qb, ga1, ga2) {
  qa1 <- toupper(qa1); qa2 <- toupper(qa2); ga1 <- toupper(ga1); ga2 <- toupper(ga2)
  fifelse(qa1 == ga1 & qa2 == ga2, qb,
          fifelse(qa1 == ga2 & qa2 == ga1, -qb, NA_real_))
}

# iv: data.table with SNP + qtl_effect (effect allele) + qtl_other + beta + se
# oc: data.table with SNP + A1 + A2 + b + se (read_ma output)
run_mr_pair <- function(iv, oc, label) {
  oc2 <- oc[, .(SNP, gA1 = A1, gA2 = A2, gb = b, gse = se)]
  m <- merge(iv, oc2, by = "SNP")
  m[, beta_h := harmonize(qtl_effect, qtl_other, beta, gA1, gA2)]
  m <- m[!is.na(beta_h) & se > 0 & gse > 0]
  out <- data.table(label = label, N_iv = nrow(iv), N_matched = nrow(m),
                    F_mean = if (nrow(m) > 0) mean((m$beta / m$se)^2) else NA_real_,
                    Method = NA_character_, Beta = NA_real_, SE = NA_real_, P = NA_real_)
  if (nrow(m) == 0) return(out)
  rows <- list()
  bi  <- which.min(m$p)
  bx  <- m$beta_h[bi]; by <- m$gb[bi]; sbx <- m$se[bi]; sby <- m$gse[bi]
  bz  <- by / bx
  bse <- sqrt(sby^2 / bx^2 + by^2 * sbx^2 / bx^4)
  rows[[1]] <- data.table(Method = "Wald(best)", Beta = bz, SE = bse, P = 2 * pnorm(-abs(bz / bse)))
  if (nrow(m) >= 3) {
    obj  <- mr_input(bx = m$beta_h, bxse = m$se, by = m$gb, byse = m$gse, snps = m$SNP)
    ivw  <- tryCatch(mr_ivw(obj), error = function(e) NULL)
    egr  <- tryCatch(mr_egger(obj), error = function(e) NULL)
    wm   <- tryCatch(mr_median(obj), error = function(e) NULL)
    if (!is.null(ivw)) rows[["IVW"]] <- data.table(Method = "IVW", Beta = ivw@Estimate, SE = ivw@StdError, P = ivw@Pvalue)
    if (!is.null(egr)) rows[["Eg"]]  <- data.table(Method = "Egger", Beta = egr@Estimate, SE = egr@StdError.Est, P = egr@Pvalue.Est)
    if (!is.null(wm))  rows[["WM"]]  <- data.table(Method = "WMed", Beta = wm@Estimate, SE = wm@StdError, P = wm@Pvalue)
  }
  r <- rbindlist(rows, fill = TRUE)
  r[, label := label]
  r[, N_iv := nrow(iv)]
  r[, N_matched := nrow(m)]
  r[, F_mean := out$F_mean[1]]
  r
}
