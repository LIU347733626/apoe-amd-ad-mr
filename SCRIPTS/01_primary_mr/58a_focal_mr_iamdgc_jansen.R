##############################################################################
# 58a: APOE focal clumped MR — same 11-instrument set as R13 primary
#       applied to IAMDGC late AMD (N=32,339) and Jansen 2019 AD (N~430k).
#
# Instrument set: revision/results/r13/apoe_pqtl_instruments_r13.csv
#   (11 APOE cis-pQTLs from UKB-PPP official clump; all level-regulating,
#    the ε4 coding tag rs429358 is NOT in the set by design).
#
# Direction gates (run BEFORE any output, per error archive 2026-09-01 #1):
#   IAMDGC: IVW β must be > 0 (↑APOE level → ↑AMD)
#           ε4 tag rs429358 (C→T, T=ε4) must be β < 0 (ε4 protective for AMD)
#   Jansen: IVW β must be < 0 (↑APOE level → ↓AD)
#           ε4 tag rs429358 (C→T, T=ε4) must be β > 0 (ε4 risk for AD)
#
# Outputs: revision/results/jansen/apoe_focus_mr_iamdgc_jansen.csv
#          revision/results/jansen/apoe_single_snp_tags_iamdgc_jansen.csv
#
# TwoSampleMR 0.7.9: FE IVW (not MRE — error archive 2026-08-14 #1),
# MR-Egger, weighted median; MRPRESSO 1.0 for outliers.
##############################################################################
suppressPackageStartupMessages({
  library(data.table); library(TwoSampleMR); library(MRPRESSO)
  library(MendelianRandomization)  # mr_input / mr_* live here in this TwoSampleMR build
})

PROJ <- "/home/huan/projects/08_APOE_AMD_AD_MR"
OUT  <- file.path(PROJ, "revision/results/jansen")
dir.create(OUT, showWarnings=F, recursive=T)

cat("===== 58a: APOE focal MR (11 instruments) x IAMDGC / Jansen =====\n")

# 1. Instrument set (identical to R13 primary) -------------------------
iv <- fread(file.path(PROJ, "revision/results/r13/apoe_pqtl_instruments_r13.csv"))
iv <- iv[, .(SNP, A1 = effect_allele.exposure, A2 = other_allele.exposure,
             beta = beta.exposure, se = se.exposure, pval = pval.exposure)]
cat(sprintf("Instruments: %d\n", nrow(iv)))
if (nrow(iv) != 11) stop("Expected 11 instruments, got ", nrow(iv), call.=FALSE)

# 2. Outcomes -----------------------------------------------------------
load_outcome <- function(path, use_header) {
  d <- fread(path, header = use_header, data.table = T, nThread = 4,
             colClasses = c("character","character","character","numeric",
                            "numeric","numeric","numeric","integer"))
  setnames(d, 1:8, c("SNP","A1","A2","freq","b","se","p","n"))
  d[!is.na(b) & se > 0, .(SNP, A1, A2, b, se, p)]
}
iamdgc <- load_outcome("/home/huan/database/GWAS/AMD/IAMDGC_late.ma", use_header = TRUE)
jansen <- load_outcome(file.path(PROJ, "smr_data/AD_Jansen2019_converted.ma"), use_header = FALSE)
cat(sprintf("IAMDGC: %d SNPs | Jansen: %d SNPs\n", nrow(iamdgc), nrow(jansen)))

# 3. Single-SNP ε4 tag direction gate (MUST pass before MR) ------------
tag_gate <- function(occ, label, expect_sign) {
  t <- occ[occ$SNP == "rs429358"]
  if (nrow(t) != 1) stop("rs429358 not found in ", label, " outcome file", call.=FALSE)
  # rs429358 canonical encoding: C -> T (T = ε4). If A1/A2 are flipped in the
  # file, re-express beta on the C->T direction.
  b <- if (t$A1 == "C" && t$A2 == "T") t$b else if (t$A1 == "T" && t$A2 == "C") -t$b
       else stop("rs429358 unexpected alleles in ", label, call.=FALSE)
  ok <- (expect_sign == "pos" && b > 0) || (expect_sign == "neg" && b < 0)
  if (!ok) stop(label, " ε4 tag direction gate FAILED: C->T beta = ",
                round(b, 4), " (expected sign: ", expect_sign, ")", call.=FALSE)
  cat(sprintf("  GATE PASS %s: rs429358 (C->T, T=ε4) beta = %.4f (p = %.3e)\n",
              label, b, t$p))
  data.table(outcome = label, SNP = "rs429358", A1 = "C", A2 = "T",
             beta = b, se = t$se, p = t$p, n = t$n, tag = "e4", check = "direction_gate")
}
e2_tag <- function(occ, label) {
  r <- occ[occ$SNP == "rs7412"]
  if (nrow(r) != 1)
    return(data.table(outcome = label, SNP = "rs7412", A1 = NA_character_, A2 = NA_character_,
                      beta = NA_real_, se = NA_real_, p = NA_real_, n = NA_integer_,
                      tag = "e2", check = "absent in file"))
  data.table(outcome = label, SNP = "rs7412", A1 = r$A1, A2 = r$A2,
             beta = r$b, se = r$se, p = r$p, n = r$n,
             tag = "e2", check = "info(19q13.32 tag LD with e2; direction not gated)")
}
tags <- rbindlist(list(
  tag_gate(iamdgc, "IAMDGC", "neg"),
  tag_gate(jansen, "Jansen", "pos"),
  e2_tag(iamdgc, "IAMDGC"),
  e2_tag(jansen, "Jansen")
), fill = TRUE)
fwrite(tags, file.path(OUT, "apoe_single_snp_tags_iamdgc_jansen.csv"), col.names = TRUE)

# 4. Harmonized focal MR (identical instrument set per outcome) --------
# TwoSampleMR 0.7.9 proven pattern (revision/20, 40a, 51):
#   mr_input() -> S4 object; mr_*() -> iRef; access via @Estimate/@StdError/@Pvalue
#   (Egger slope: @StdError.Est/@Pvalue.Est; Egger intercept: @Intercept/@StdError.Int/@Pvalue.Int)
#   CI reported as beta ± 1.96*SE (paper statistical-brief convention).
ci196 <- function(b, s) if (is.na(b) || is.na(s)) c(NA_real_, NA_real_) else c(b - 1.959964 * s, b + 1.959964 * s)

run_focal <- function(occ, outcome_name, expected_sign) {
  # iv cols:  SNP, A1, A2, beta, se, pval  (renamed at load, line 34)
  # occ cols: SNP, A1, A2, b, se, p
  # merge suffixes .iv/.oc apply ONLY to duplicated cols: A1, A2, se.
  # 'beta' (iv) does NOT collide with 'b' (occ) -> stays unsuffixed.
  # SUFFIX RULE (learned, 58a v3-v5): data.table merge suffixes apply ONLY to
  # duplicated cols (A1, A2, se). iv-only 'beta'/'pval' and occ-only 'b'/'p'
  # keep their plain names in the merged table.
  m <- merge(iv, occ, by = "SNP", suffixes = c(".iv", ".oc"))
  m[, beta_oc := ifelse(toupper(A1.iv) == toupper(A1.oc) &
                          toupper(A2.iv) == toupper(A2.oc), b,
                ifelse(toupper(A1.iv) == toupper(A2.oc) &
                          toupper(A2.iv) == toupper(A1.oc), -b, NA_real_))]
  n_drop <- sum(is.na(m$beta_oc))
  m <- m[!is.na(beta_oc) & se.oc > 0 & se.iv > 0]
  cat(sprintf("%s: %d/%d instruments matched & harmonized (dropped %d)\n",
              outcome_name, nrow(m), nrow(iv), n_drop))
  if (nrow(m) < 3) stop(outcome_name, ": need >=3 harmonized instruments, got ", nrow(m), call.=FALSE)

  bx   <- m$beta       # iv exposure beta (unique col name, no suffix)
  bxse <- m$se.iv
  by   <- m$beta_oc
  byse <- m$se.oc
  snps <- m$SNP

  obj     <- mr_input(bx = bx, bxse = bxse, by = by, byse = byse, snps = snps)
  ivw     <- tryCatch(mr_ivw(obj),  error = function(e) NULL)
  egger   <- tryCatch(mr_egger(obj), error = function(e) NULL)
  med     <- tryCatch(mr_median(obj), error = function(e) NULL)
  pressed <- tryCatch(
    # FIX(2026-09-10): mr_presso 的 data= 无默认值；原调用缺该参数 → 报错被
    # tryCatch 静默吞掉 → MR-PRESSO 从未运行。且 BetaOutcome 等须为 data 中的
    # 列名（原代码传的是向量）。
    mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
              SdOutcome = "sd.outcome", SdExposure = "sd.exposure",
              data = data.frame(beta.outcome = by, beta.exposure = bx,
                                sd.outcome = byse, sd.exposure = bxse),
              OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
              NbDistribution = 1000, SignifThreshold = 0.05,
              seed = 20260910),
    error = function(e) NULL)
  if (is.null(ivw)) stop(outcome_name, ": FE IVW failed", call.=FALSE)

  rows <- list()
  ci  <- ci196(ivw@Estimate, ivw@StdError)
  rows[[1]] <- data.table(
    outcome = outcome_name, method = "IVW", n_snp = length(bx),
    beta = ivw@Estimate, se = ivw@StdError, p = ivw@Pvalue,
    ci_lo = ci[1], ci_hi = ci[2],
    het_q = ivw@Heter.Stat[1], het_p = ivw@Heter.Stat[2], egger_int_p = NA_real_)
  if (!is.null(egger)) {
    ci  <- ci196(egger@Estimate, egger@StdError.Est)
    rows[[length(rows) + 1]] <- data.table(
      outcome = outcome_name, method = "MR_Egger", n_snp = length(bx),
      beta = egger@Estimate, se = egger@StdError.Est, p = egger@Pvalue.Est,
      ci_lo = ci[1], ci_hi = ci[2],
      het_q = NA_real_, het_p = NA_real_, egger_int_p = egger@Pvalue.Int)
  }
  if (!is.null(med)) {
    ci  <- ci196(med@Estimate, med@StdError)
    rows[[length(rows) + 1]] <- data.table(
      outcome = outcome_name, method = "WMed", n_snp = length(bx),
      beta = med@Estimate, se = med@StdError, p = med@Pvalue,
      ci_lo = ci[1], ci_hi = ci[2],
      het_q = NA_real_, het_p = NA_real_, egger_int_p = NA_real_)
  }
  if (!is.null(pressed) && length(pressed[["Main MR results"]]) > 0) {
    mr  <- pressed[["Main MR results"]]
    ci  <- ci196(mr[["Causal Estimate"]][1], mr[["Sd"]][1])
    rows[[length(rows) + 1]] <- data.table(
      outcome = outcome_name, method = "MR-PRESSO_raw", n_snp = length(bx),
      beta = mr[["Causal Estimate"]][1], se = mr[["Sd"]][1], p = mr[["P-value"]][1],
      ci_lo = ci[1], ci_hi = ci[2],
      het_q = NA_real_, het_p = NA_real_, egger_int_p = NA_real_)
  }
  r <- rbindlist(rows, fill = TRUE)

  # direction gate on IVW
  b <- r[method == "IVW", beta]
  ok <- (expected_sign == "pos" && b > 0) || (expected_sign == "neg" && b < 0)
  if (!ok) stop(outcome_name, " IVW direction gate FAILED: beta = ", round(b, 4), call.=FALSE)
  cat(sprintf("  GATE PASS %s: IVW beta = %.4f (p = %.2e) — expected %s\n",
              outcome_name, b, r[method == "IVW", p], expected_sign))
  cat(sprintf("  MR-PRESSO available: %s\n", !is.null(pressed) && length(pressed[["Main MR results"]]) > 0))

  r
}
res <- rbindlist(list(
  run_focal(iamdgc, "IAMDGC_late_AMD", "pos"),
  run_focal(jansen, "Jansen2019_AD",   "neg")
), fill = TRUE)
fwrite(res, file.path(OUT, "apoe_focus_mr_iamdgc_jansen.csv"), col.names = TRUE)

cat("\n=== Focal MR results (11-instrument set, FE IVW) ===\n")
print(res[, .(outcome, method, n_snp, beta, se, p, ci_lo, ci_hi, egger_int_p)], row.names = FALSE)
cat(sprintf("\nResults dir: %s\n", OUT))
cat("--- Session info ---\n"); print(sessionInfo()$loadedOnly)
cat("\nDone.\n")
