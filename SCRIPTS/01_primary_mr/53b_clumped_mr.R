#!/usr/bin/env Rscript
# =============================================================================
# 53b_clumped_mr.R — deCODE APOE protein pQTL GWAS -> MR with LD-clump grid
# =============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(MendelianRandomization)
})
PROJ <- "/home/huan/projects/08_APOE_AMD_AD_MR"
RES  <- file.path(PROJ, "revision", "results", "r13")
OUT  <- file.path(RES, "decode_protein_mr_clump_grid.csv")
dir.create(RES, showWarnings = FALSE, recursive = TRUE)
set.seed(42)
## ---------- Load plink r2 matrix ----------
LD_REF <- "/tmp/apoe_ld_1kg.ld"
ld_full <- fread(LD_REF)
ld_full <- rbindlist(list(
  ld_full,
  ld_full[, .(CHR_A = CHR_B, BP_A = BP_B, SNP_A = SNP_B,
              CHR_B = CHR_A, BP_B = BP_A, SNP_B = SNP_A, R2 = R2)]
), use.names = TRUE)
## ---------- Instrument positions ----------
instrument_positions <- data.table(
  SNP = c("rs143539524", "rs77332277", "rs60286743", "rs57355367",
          "rs58819137",  "rs148023052", "rs113928991",
          "rs777208606", "rs766506382"),
  chr = 19,
  pos = c(45552198, 45690709, 45698052, 45698085,
          45698174, 45698590, 45700180,
          45218908,  45536359)
)
# Add self-pairs (r2=1) for all instrument positions
self_pairs <- instrument_positions[, .(SNP_A = SNP, CHR_A = chr, 
                                       BP_A = pos, CHR_B = chr, BP_B = pos, 
                                       SNP_B = SNP, R2 = 1.0)]
ld_full <- rbindlist(list(ld_full, self_pairs), use.names = TRUE)
# Zero out any NA R2 (SNPs not overlapping with 1000G -> treat as independent, r2=0)
n_na_before <- sum(is.na(ld_full$R2))
ld_full[is.na(R2), R2 := 0.0]
## ---------- Load deCODE Tier B ----------
ivB <- fread(file.path(RES, "decode_protein_iv_tierB.csv"))
cat(sprintf("[53b] Tier B: %d instruments loaded\n", nrow(ivB)))
## ---------- Load deCODE Tier A ----------
ivA <- fread("/home/huan/database/pQTL/pqtl/decode/deCODE-2021-data/decode_2021_clump/decode2021-5e-08-clumped-10000kb-0.001r2-20240307.csv")
ivA <- ivA[gene.exposure == "APOE" & mr_keep.exposure == TRUE]
ivA <- ivA[, .(SNP, effect_allele.exposure, other_allele.exposure,
               beta.exposure, se.exposure, pval.exposure)]
cat(sprintf("[53b] Tier A: %d instruments (deCODE 6-SNP clumped)\n", nrow(ivA)))
## ---------- Outcome files ----------
AMD_fp   <- file.path(PROJ, "smr_data", "AMD_r13.ma")
AD_fp    <- file.path(PROJ, "smr_data", "AD_r13.ma")
IAMD_fp  <- file.path(PROJ, "smr_data", "iamdgc", "iamdgc_amd_for_mr.txt")
read_ma <- function(path, has_header = TRUE) {
  d <- fread(path, header = has_header, data.table = TRUE)
  if (ncol(d) == 9) setnames(d, c("SNP","A1","A2","freq","b","se","p","n","x"))
  else if (ncol(d) == 8) setnames(d, c("SNP","A1","A2","freq","b","se","p","n"))
  else if (ncol(d) == 6) setnames(d, c("SNP","A1","A2","b","se","p"))
  d[nchar(A1) <= 1 & nchar(A2) <= 1, ]
  d[!is.na(b) & !is.na(se) & se > 0, ]
}
amd_d <- read_ma(AMD_fp, has_header = TRUE)
ad_d  <- read_ma(AD_fp,  has_header = FALSE)
iam_d <- read_ma(IAMD_fp, has_header = TRUE)
cat(sprintf("[53b] AMD_r13=%d  AD_r13=%d  IAMDGC=%d\n",
            nrow(amd_d), nrow(ad_d), nrow(iam_d)))

## Pre-extract outcome SNPs (only the 9 APOE-region rsids) for fast MR
amd_snp <- amd_d[amd_d$SNP %in% c("rs143539524","rs77332277","rs60286743","rs57355367",
                                    "rs58819137","rs148023052","rs113928991","rs777208606","rs766506382"),
                                    .(SNP, A1, A2, b, se)]
ad_snp  <- ad_d [ad_d$SNP %in% c("rs143539524","rs77332277","rs60286743","rs57355367",
                                  "rs58819137","rs148023052","rs113928991","rs777208606","rs766506382"),
                                  .(SNP, A1, A2, b, se)]
iam_snp <- iam_d[iam_d$SNP %in% c("rs143539524","rs77332277","rs60286743","rs57355367",
                                   "rs58819137","rs148023052","rs113928991","rs777208606","rs766506382"),
                                   .(SNP, A1, A2, b, se)]
cat(sprintf("[53b] AMD=%d  AD=%d  IAMDGC=%d pre-extracted SNPs\n",
            nrow(amd_snp), nrow(ad_snp), nrow(iam_snp)))
## ---------- harmonize + merge + MR (following 53) ----------
harmonize <- function(qa1, qa2, qb1, qb2, qbeta) {
  qa1 <- toupper(qa1); qa2 <- toupper(qa2)
  qb1 <- toupper(qb1); qb2 <- toupper(qb2)
  fifelse(qa1 == qb1 & qa2 == qb2, qbeta,
          fifelse(qa1 == qb2 & qa2 == qb1, -qbeta, NA_real_))
}
merge_iv <- function(exposure, outcome, outcome_name) {
  inc <- exposure[, .(SNP, effect_allele.exposure, other_allele.exposure,
                      beta.exposure, se.exposure)]
  m <- merge(inc, outcome, by = "SNP", suffixes = c(".iv", ".gw"))
  m <- m[!is.na(A1) & !is.na(A2), ]
  m[, beta_harm := harmonize(effect_allele.exposure, other_allele.exposure,
                             A1, A2, beta.exposure)]
  m <- m[!is.na(beta_harm) & !is.na(se.exposure) & !is.na(se) &
          se.exposure > 0 & se > 0, ]
  cat(sprintf("[53b] %s: %d instruments after merge\n", outcome_name, nrow(m)))
  m
}
run_mr <- function(m, label, tier) {
  n <- nrow(m)
  if (n < 3) return(data.table(outcome = label, tier = tier,
                               n_instruments = n, method = "TOO_FEW",
                               beta = NA_real_, se = NA_real_, p = NA_real_))
  obj <- tryCatch(mr_input(bx = m$beta_harm, bxse = m$se.exposure,
                           by = m$b, byse = m$se, snps = m$SNP),
                  error = function(e) {
    cat(sprintf("[53b] %s %s mr_input failed: %s\n",
                tier, label, conditionMessage(e)))
    NULL
  })
  if (is.null(obj)) return(data.table(outcome = label, tier = tier,
                                      n_instruments = n, method = "FAIL_INPUT",
                                      beta = NA_real_, se = NA_real_, p = NA_real_))
  rows <- list()
  t0 <- Sys.time()
  err_log <- function(fn) function(e) {
    cat(sprintf("  [53b] %s failed: %s\n", fn, conditionMessage(e)))
    NULL
  }
  ivw   <- tryCatch(mr_ivw(obj),   error = err_log("mr_ivw"))
  egger <- tryCatch(mr_egger(obj), error = err_log("mr_egger"))
  wmed  <- tryCatch(mr_median(obj), error = err_log("mr_median"))
  cat(sprintf("  [53b %s %s] MR 3 methods: %.3fs\n", tier, label, as.numeric(Sys.time() - t0, units = "secs")))
  if (!is.null(ivw))
    rows[["IVW"]]   <- data.table(method = "IVW",   beta = ivw@Estimate,
                                  se = ivw@StdError, p = ivw@Pvalue)
  if (!is.null(egger))
    rows[["Egger"]] <- data.table(method = "Egger", beta = egger@Estimate,
                                  se = egger@StdError.Est, p = egger@Pvalue.Est)
  if (!is.null(wmed))
    rows[["WMed"]]  <- data.table(method = "WMed",  beta = wmed@Estimate,
                                  se = wmed@StdError, p = wmed@Pvalue)
  r <- rbindlist(rows, fill = TRUE)
  r[, outcome := label]
  r[, tier := tier]
  r[, n_instruments := n]
  r
}
## ---------- Clump grid ----------
clump_instruments <- function(r2thr, window_mb) {
  ivB_sorted <- ivB[order(-log10(pval.exposure), decreasing = TRUE), ]
  keep <- logical(nrow(ivB_sorted))
  for (i in seq_len(nrow(ivB_sorted))) {
    s_cur <- ivB_sorted[i, SNP]
    p_cur <- instrument_positions[SNP == s_cur, pos]
    if (is.na(p_cur)) p_cur <- 0
    drop_it <- FALSE
    for (j in seq_len(i - 1)) {
      s_prev <- ivB_sorted[j, SNP]
      r2_idx <- ld_full[SNP_A == s_cur & SNP_B == s_prev][, "R2"]
      if (length(r2_idx) == 0L) {
        r2_val <- 0.0
      } else {
        r2_val <- as.numeric(r2_idx[1L])
        if (is.na(r2_val)) r2_val <- 0.0
      }
      p_prev <- instrument_positions[SNP == s_prev, pos]
      if (is.na(p_prev)) p_prev <- 0
      if (r2_val >= r2thr &&
          abs(p_cur - p_prev) <= window_mb * 1e6) {
        drop_it <- TRUE; break
      }
    }
    keep[i] <- !drop_it
  }
  ivB_sorted[keep]
}
## ---------- Build jobs ----------
jobs <- list()
for (r2thr in c(0.0001, 0.001, 0.01)) {
  for (wm in c(5, 10)) {
    lab <- sprintf("r2=%.4f_w=%dMb", r2thr, wm)
    y <- tryCatch(clump_instruments(r2thr, wm), error = function(e) {
      cat("clump failed for", lab, ":", conditionMessage(e), "\n"); NULL
    })
    if (is.null(y) || nrow(y) == 0) {
      jobs[[length(jobs) + 1]] <- list(iv = NULL, nm = NULL, out = NULL, tier = lab)
      next
    }
    cat(sprintf("[53b] %s -> %d instruments: %s\n",
                lab, nrow(y), paste(y$SNP, collapse = ",")))
    for (nm in c("AMD_r13", "AD_r13", "IAMDGC_late")) {
      # FIX(2026-09-10): 原为未定义的 amd_e，应为 amd_d（见 :61 定义）
      m <- if (nm == "AMD_r13") amd_d
           else if (nm == "AD_r13") ad_d
           else iam_d
      jobs[[length(jobs) + 1]] <- list(iv = y, nm = nm, out = m, tier = lab)
    }
  }
}
for (nm in c("AMD_r13", "AD_r13", "IAMDGC_late")) {
  # FIX(2026-09-10): 原为未定义的 amd_e，应为 amd_d（见 :61 定义）
  m <- if (nm == "AMD_r13") amd_d
       else if (nm == "AD_r13") ad_d
       else iam_d
  jobs[[length(jobs) + 1]] <- list(iv = ivA, nm = nm, out = m, tier = "A_6SNP")
}
chunks <- list()
for (k in seq_along(jobs)) {
  j <- jobs[[k]]
  if (is.null(j$out)) next
  tb <- Sys.time()
  mm <- merge_iv(j$iv, j$out, j$nm)
  cat(sprintf("[53b] job %d/%d merge: %.3fs (n=%d)\n",
              k, length(jobs),
              as.numeric(Sys.time() - tb, units = "secs"), nrow(mm)))
  tr <- Sys.time()
  rows_<- run_mr(mm, j$nm, j$tier)
  cat(sprintf("[53b] job %d/%d MR: %.3fs\n",
              k, length(jobs),
              as.numeric(Sys.time() - tr, units = "secs")))
  chunks <- c(chunks, list(rows_))
}
res <- rbindlist(chunks, fill = TRUE)
res <- res[, .(outcome, tier, n_instruments, method, beta, se, p)]
res <- res[order(outcome, tier, method)]
cat("\n[53b] COMPLETE CLUMP-GRID RESULTS\n")
print(res)
fwrite(res, OUT, col.names = TRUE)
cat(sprintf("\n[53b] DONE -> %s\n", OUT))
