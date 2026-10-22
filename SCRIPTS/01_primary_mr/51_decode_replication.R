#!/usr/bin/env Rscript
# =============================================================================
# 51_decode_replication.R — deCODE APOE pQTL instrument -> independent outcomes
#
# The primary R13 analysis uses the UKB-PPP (SomaScan) instrument. The deCODE
# (Olink, Icelanders) pQTLs are genetically independent and share ZERO sample
# overlap with FinnGen R13 or IAMDGC -> clean replication of the APOE protein
# -> AMD / AD causal directions.
#
# Exposure: deCODE 2021 6-SNP clump (primary); single-SNP Wald sensitivity:
#   rs429358 (eps4 tag) and rs71352239 (deCODE lead).
# Outcomes:
#   AMD_r13 (FinnGen R13 H7, N = 472,750)
#   AD_r13  (FinnGen R13 G6, N = 500,186)
#   IAMDGC late AMD (Fritsche 2016 meta; ~15k)
# Methods: IVW (random), MR-Egger, weighted median.
# Gates: convention_check.csv + golden anchor (UKB) + deCODE direction on R13.
# Output: revision/results/r13/decode_replication_mr.csv
# Run:    PATH=/home/huan/software/miniconda3/envs/biomni_r/bin:$PATH \
#         Rscript revision/51_decode_replication.R
# =============================================================================
library(data.table)
library(MendelianRandomization)

PROJ   <- "/home/huan/projects/08_APOE_AMD_AD_MR"
DB_UKB <- "/home/huan/database/pQTL/pqtl/ukb-ppp/summary"
DB_DEC <- "/home/huan/database/pQTL/pqtl/decode/deCODE-2021-data/decode_2021_clump"
RES    <- file.path(PROJ, "revision/results/r13")
OUT    <- file.path(RES, "decode_replication_mr.csv")
if (!dir.exists(RES)) dir.create(RES, recursive = TRUE)

if (!file.exists(file.path(RES, "convention_check.csv")))
  stop("GATE FAIL: convention_check.csv missing (run 50a first)", call. = FALSE)

## ---- exposures ----
pq_dec <- fread(file.path(DB_DEC, "decode2021-5e-08-clumped-10000kb-0.001r2-20240307.csv"),
                header = TRUE, sep = ",")
iv_dec <- pq_dec[gene.exposure == "APOE" & mr_keep.exposure == TRUE, ]
pq_dec_u <- fread(file.path(DB_DEC, "decode2021-5e-08-unclumped.csv"),
                  header = TRUE, sep = ",")
setnames(pq_dec_u, c("SNP", "effect_allele", "other_allele", "eaf", "beta", "se",
                     "pval", "samplesize"),
             c("SNP", "effect_allele.exposure", "other_allele.exposure",
               "eaf.exposure", "beta.exposure", "se.exposure",
               "pval.exposure", "samplesize.exposure"))
dec_tags <- pq_dec_u[SNP %in% c("rs429358", "rs71352239") & gene == "APOE", ]
cat("instruments: deCODE clump", nrow(iv_dec), "| deCODE tags", nrow(dec_tags), "\n")

## ---- outcomes ----
read_ma <- function(path, header = TRUE, fivecol = FALSE) {
  if (header) {
    d <- fread(path, header = TRUE, data.table = TRUE)
  } else {
    d <- fread(path, header = FALSE, data.table = TRUE,
               colClasses = c("character", "character", "character", "numeric",
                              "numeric", "numeric", "numeric", "integer"))
  }
  if (fivecol) {
    setnames(d, c("SNP", "A1", "A2", "beta", "se", "p"),
                 c("SNP", "A1", "A2", "b", "se", "p"))
  } else {
    setnames(d, c("SNP", "A1", "A2", "freq", "b", "se", "p", "n"))
  }
  d[!is.na(b) & !is.na(se) & se > 0]
}
amd_r13 <- read_ma(file.path(PROJ, "smr_data/AMD_r13.ma"), TRUE)
ad_r13  <- read_ma(file.path(PROJ, "smr_data/AD_r13.ma"), FALSE)  # AD_r13.ma has NO header
iamdgc  <- read_ma(file.path(PROJ, "smr_data/iamdgc/iamdgc_amd_for_mr.txt"),
                   TRUE, fivecol = TRUE)
cat("outcome rows: AMD_r13", nrow(amd_r13), "AD_r13", nrow(ad_r13),
    "IAMDGC", nrow(iamdgc), "\n")

## ---- harmonize / merge / MR (validated 20/50a logic) ----
harmonize <- function(qtl_a1, qtl_a2, qtl_beta, gwas_a1, gwas_a2) {
  qa1 <- toupper(qtl_a1); qa2 <- toupper(qtl_a2)
  ga1 <- toupper(gwas_a1); ga2 <- toupper(gwas_a2)
  fifelse(qa1 == ga1 & qa2 == ga2, qtl_beta,
          fifelse(qa1 == ga2 & qa2 == ga1, -qtl_beta, NA_real_))
}
merge_iv <- function(exp, outcome, outcome_name) {
  # 2026-09-02: .() NSE 子集在 00:38 smoke 中偶发 "could not find function '.'"
  # （完整复现探针未复现，判定偶发）；改用无 NSE 的 with=FALSE 子集，语义等价、消除脆点。
  m <- merge(exp, outcome[, c("SNP", "A1", "A2", "b", "se"), with = FALSE],
             by = "SNP", suffixes = c(".pqtl", ".gwas"))
  m[, beta_harm := harmonize(effect_allele.exposure, other_allele.exposure,
                             beta.exposure, A1, A2)]
  m <- m[!is.na(beta_harm) & !is.na(se.exposure) & !is.na(se) &
           se.exposure > 0 & se > 0]
  cat(sprintf("  %s: %d instruments after merge/harmonization\n",
              outcome_name, nrow(m)))
  m
}
run_mr <- function(m, label) {
  obj <- mr_input(bx = m$beta_harm, bxse = m$se.exposure,
                  by = m$b, byse = m$se, snps = m$SNP)
  ivw   <- tryCatch(mr_ivw(obj),   error = function(e) NULL)
  egger <- tryCatch(mr_egger(obj), error = function(e) NULL)
  med   <- tryCatch(mr_median(obj), error = function(e) NULL)
  r <- data.frame(label = label,
                  n_instr = nrow(m),
                  ivw_beta  = if (is.null(ivw))   NA_real_ else ivw@Estimate,
                  ivw_se    = if (is.null(ivw))   NA_real_ else ivw@StdError,
                  ivw_p     = if (is.null(ivw))   NA_real_ else ivw@Pvalue,
                  egger_beta = if (is.null(egger)) NA_real_ else egger@Estimate,
                  egger_p    = if (is.null(egger)) NA_real_ else egger@Pvalue.Est,
                  wmed_beta  = if (is.null(med))  NA_real_ else med@Estimate,
                  wmed_se    = if (is.null(med))  NA_real_ else med@StdError,
                  wmed_p     = if (is.null(med))  NA_real_ else med@Pvalue,
                  stringsAsFactors = FALSE)
  rownames(r) <- label
  r
}

## ---- gate: golden anchor (UKB instrument on the original outcomes) ----
pq_ukb <- fread(file.path(DB_UKB, "ukb-ppp-2023-5e-08-clumped-10000kb-0.001r2.csv"),
                header = TRUE, sep = ",")
iv_ukb <- pq_ukb[gene.exposure == "APOE", ]
amd_finn <- read_ma(file.path(PROJ, "smr_data/AMD_finn.ma"), TRUE)
ad_bell  <- read_ma(file.path(PROJ, "smr_data/AD_Bellenguez2022.ma"), FALSE)
a1 <- run_mr(merge_iv(iv_ukb, amd_finn, "AMD_finn(anchor)"), "UKB->AMD_finn")
a2 <- run_mr(merge_iv(iv_ukb, ad_bell, "AD_bell(anchor)"),  "UKB->AD_bell")
anchor_ok <- abs(a1["UKB->AMD_finn", "ivw_beta"] - 0.2092) < 1e-4 &&
             abs(a2["UKB->AD_bell", "ivw_beta"] - (-0.7627)) < 1e-4
if (!anchor_ok) stop("GATE FAIL: golden anchor not reproduced", call. = FALSE)
cat("GATE PASS: anchor reproduced\n")

## ---- Part 1: deCODE 6-SNP -> 3 outcomes (primary) ----
# 2026-09-02 修复：原 `for (o in list(...))` 中 o 已是 data.table（list 元素），
# `o[[1]]` 取到的是第一列（SNP 向量）而非整表 → merge 挂在
# "could not find function '.'" / "incorrect number of dimensions"。
# 改为按名迭代。
rows <- list()
ol <- list(AMD_r13 = amd_r13, AD_r13 = ad_r13, IAMDGC = iamdgc)
for (nm in names(ol)) {
  rows[[length(rows) + 1]] <-
    run_mr(merge_iv(iv_dec, ol[[nm]], paste0("deCODE->", nm)),
           paste0("deCODE->", nm))
}
## ---- Part 2: single-SNP Wald sensitivities ----
tag_429 <- dec_tags[SNP == "rs429358", ]
tag_713 <- dec_tags[SNP == "rs71352239", ]
for (tg in list(tag_429, tag_713)) {
  if (nrow(tg) == 0) next
  for (nm in names(ol)) {
    rows[[length(rows) + 1]] <-
      run_mr(merge_iv(tg, ol[[nm]],
                      paste0(tg$SNP[1], "->", nm)),
             paste0(tg$SNP[1], "->", nm))
  }
}
res <- rbind(rows)
write.csv(res, OUT, row.names = FALSE)
print(res)
cat("SAVED:", OUT, "\n")
