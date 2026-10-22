#!/usr/bin/env Rscript
# 97_iv_qc_and_rerun.R — 工具变量严格质控 + 剔除后重跑 MR（2026-09-10）
# 用户要求：严格检测 MR 每一步，找出需剔除的 SNP，剔除后重跑
# 输出：revision/results/r13/iv_qc_report.csv（逐工具判定）
#       revision/results/r13/mr_after_qc_{AMD_R13,AD_R13,IAMDGC}.csv（剔除前后对比）
suppressPackageStartupMessages({library(data.table); library(MendelianRandomization)})
PROJ <- "/home/huan/projects/08_APOE_AMD_AD_MR"
OUT  <- file.path(PROJ, "revision", "results", "r13")
PQ   <- "/home/huan/database/pQTL/pqtl/ukb-ppp/summary/ukb-ppp-2023-5e-08-clumped-10000kb-0.001r2.csv"
PLINK<- "/home/huan/software/configs_G19/bin/plink"; BFILE <- "/home/huan/references/1000g/EU1000QC"
ENDP <- list(AMD_R13=file.path(PROJ,"smr_data/AMD_r13.ma"),
             AD_R13 =file.path(PROJ,"smr_data/AD_r13_hdr.ma"),
             IAMDGC ="/home/huan/database/GWAS/AMD/IAMDGC_late.ma")

pq <- fread(PQ); iv <- pq[gene.exposure == "APOE"]
iv[, F_stat := (beta.exposure/se.exposure)^2]
iv[, maf := pmin(eaf.exposure, 1 - eaf.exposure)]
iv[, palindromic := toupper(paste0(effect_allele.exposure, other_allele.exposure)) %in% c("AT","TA","CG","GC")]

## ---- 1) 工具间 LD（1000G EUR）----
sf <- file.path(PROJ,"tmp/qc2_snps.txt"); writeLines(iv$SNP, sf)
ldpfx <- file.path(PROJ,"tmp/qc2_ld")
system(sprintf("%s --bfile %s --extract %s --r square --write-snplist --out %s --silent", PLINK, BFILE, sf, ldpfx))
ids <- readLines(paste0(ldpfx,".snplist")); ld <- as.matrix(read.table(paste0(ldpfx,".ld")))
colnames(ld) <- rownames(ld) <- ids
ld2 <- ld^2; diag(ld2) <- 0
iv[, max_r2_with_other := sapply(SNP, function(s) if (s %in% ids) max(ld2[s, ]) else NA_real_)]
iv[, n_partners_gt_0001 := sapply(SNP, function(s) if (s %in% ids) sum(ld2[s, ] > 0.001) else NA_integer_)]

## ---- 2) 回文变异的链方向验证（暴露 vs 结局 EAF 一致性）----
read_out <- function(f) {
  d <- fread(cmd = sprintf("grep -w -F -f %s %s", sf, f), header = FALSE)
  setnames(d, c("SNP","A1","A2","freq","b","se","p","n"))
  d[, p := as.numeric(p)][!is.na(b) & !is.na(se) & se > 0]
}
OUTD <- list(); for (nm in names(ENDP)) OUTD[[nm]] <- read_out(ENDP[[nm]])

pal_check <- rbindlist(lapply(names(OUTD), function(nm) {
  o <- OUTD[[nm]]; m <- merge(iv[palindromic == TRUE, .(SNP, ea=effect_allele.exposure, eaf_exposure=eaf.exposure)],
                              o[, .(SNP, A1, A2, freq, b, se)], by="SNP")
  if (nrow(m)==0) return(NULL)
  m[, `:=`(eaf_outcome = fifelse(toupper(A1)==toupper(ea), freq, 1-freq),
           eaf_diff = abs(eaf_exposure - fifelse(toupper(A1)==toupper(ea), freq, 1-freq)))]
  m[, outcome := nm]
  m[, .(outcome, SNP, ea, eaf_exposure=round(eaf_exposure,5), eaf_outcome=round(eaf_outcome,5),
        eaf_diff=round(eaf_diff,5), b_out=round(b,4))]
}), fill=TRUE)

## ---- 3) 剔除决策 ----
# 规则 R1: 回文且两侧 EAF 差异 > 0.05（同 SNP 同链应几乎相等）→ 链方向存疑，剔除
# 规则 R2: 工具间 r² > 0.001 → 剔除 p 值较大（较弱）的那一个
# 规则 R3: F < 10 或 MAF < 0.01 → 剔除（本集无）
drop <- character(0); reasons <- character(0)
if (nrow(pal_check) > 0) {
  bad <- pal_check[eaf_diff > 0.05]
  for (s in unique(bad$SNP)) { drop <- c(drop, s); reasons <- c(reasons, sprintf("%s: 回文且暴露/结局 EAF 差 >0.05（链方向存疑）", s)) }
}
pairs <- which(ld2 > 0.001 & upper.tri(ld2), arr.ind = TRUE)
for (k in seq_len(nrow(pairs))) {
  a <- ids[pairs[k,1]]; b <- ids[pairs[k,2]]
  pa <- iv$pval.exposure[iv$SNP==a]; pb <- iv$pval.exposure[iv$SNP==b]
  weak <- if (pa > pb) a else b
  drop <- c(drop, weak)
  reasons <- c(reasons, sprintf("%s: 与 %s 的 r²=%.4f 超 0.001 阈值，二者取 p 较大者剔除", weak, if (weak==a) b else a, ld2[pairs[k,1],pairs[k,2]]))
}
drop <- unique(drop)

cat("=== 工具变量质控判定 ===\n")
print(iv[, .(SNP, ea=effect_allele.exposure, oa=other_allele.exposure, maf=round(maf,4),
             F=round(F_stat,1), p=pval.exposure, palindromic,
             max_r2=round(max_r2_with_other,5), n_partners=n_partners_gt_0001)])
cat("\n=== 回文变异 EAF 一致性（暴露 vs 结局）===\n")
if (nrow(pal_check)) print(pal_check) else cat("  （无）\n")
cat("\n=== 剔除清单 ===\n")
if (length(drop)) { for (k in seq_along(drop)) cat(sprintf("  - %s\n", reasons[k])) } else cat("  （无）\n")
cat(sprintf("\n保留 %d / %d 个工具\n", nrow(iv)-length(drop), nrow(iv)))

iv[, verdict := fifelse(SNP %in% drop, "EXCLUDE", "KEEP")]
iv[, reason := sapply(SNP, function(s) { i <- which(drop==s); if (length(i)) reasons[i[1]] else "" })]
fwrite(iv[, .(SNP, effect_allele=effect_allele.exposure, other_allele=other_allele.exposure,
              eaf=eaf.exposure, beta=beta.exposure, se=se.exposure, p=pval.exposure,
              F_stat, maf, palindromic, max_r2_with_other, verdict, reason)],
       file.path(OUT, "iv_qc_report.csv"))

## ---- 4) 剔除前后重跑 MR ----
res <- list()
for (nm in names(ENDP)) {
  o <- OUTD[[nm]]
  m <- merge(iv, o, by="SNP")
  m[, beta_h := fifelse(toupper(effect_allele.exposure)==toupper(A1) &
                          toupper(other_allele.exposure)==toupper(A2), beta.exposure,
                fifelse(toupper(effect_allele.exposure)==toupper(A2) &
                          toupper(other_allele.exposure)==toupper(A1), -beta.exposure, NA_real_))]
  m <- m[!is.na(beta_h) & !duplicated(SNP)]
  for (tag in c("all", "after_qc")) {
    s <- if (tag=="all") m else m[!SNP %in% drop]
    if (nrow(s) < 2) next
    obj <- mr_input(bx=s$beta_h, bxse=s$se.exposure, by=s$b, byse=s$se, snps=s$SNP)
    for (mdl in c("fixed","random")) {
      r <- tryCatch(mr_ivw(obj, model=mdl), error=function(e) NULL); if (is.null(r)) next
      e <- tryCatch(mr_egger(obj), error=function(e) NULL)
      res[[length(res)+1]] <- data.table(outcome=nm, set=tag, model=mdl, k=nrow(s),
        beta=r@Estimate, se=r@StdError, p=r@Pvalue,
        Q=tryCatch(r@Heter.Stat[1],error=function(e) NA_real_), Q_p=tryCatch(r@Heter.Stat[2],error=function(e) NA_real_),
        egger_int_p=if(!is.null(e)) e@Pvalue.Int else NA_real_)
    }
  }
}
R <- rbindlist(res)
fwrite(R, file.path(OUT, "mr_after_qc_summary.csv"))
cat("\n=== 剔除前后对比 ===\n")
print(R[, .(outcome, set, model, k, beta=round(beta,4), se=round(se,4), p=signif(p,3), Q_p=signif(Q_p,3))])
cat("\nDONE\n")
