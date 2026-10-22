##############################################################################
# Mediation analysis, R13 re-estimation (ε4 → APOE level → AMD / AD)
#
# _17_mediation.R 同法 (product of coefficients + delta method)。
# R13 重估只把"总效应"换成 R13,三条路径统一到当前主链:
#   - path a (ε4→APOE):     UKB-PPP rs429358 (unclumped),不变
#   - path b (APOE→disease): R13 focal clumped MR (apoe_focus_mr_clumped.csv)
#   - Total (ε4→disease):   R13 rs429358 per-allele (AMD_r13.ma / AD_r13.ma)
#
# Scale caveat: a (pQTL β) 与 b (MR β) 须同一 APOE 尺度 (均 Z 或均 raw) 时
#   proportion 才是严格比例;否则按 semi-quantitative / 方向性解读。
#   Sobel / indirect 的 p 是方向证据 (与尺度无关)。
##############################################################################
library(data.table)

PROJ <- "/home/huan/projects/08_APOE_AMD_AD_MR"
OUT  <- file.path(PROJ, "apoe_results")
DB   <- "/home/huan/database/01_raw"

cat("===== Mediation: R13 re-estimation (ε4 → APOE → AMD/AD) =====\n\n")

# ---------- path a: ε4 → APOE (UKB-PPP rs429358) ----------
pq <- fread(file.path(DB, "pqtl/ukb-ppp/summary/ukb-ppp-2023-5e-08-unclumped.csv"),
            header = TRUE, sep = ",", data.table = TRUE,
            select = c("SNP", "beta", "se", "pval", "gene"))
e4 <- pq[gene == "APOE" & SNP == "rs429358", ]
if (nrow(e4) == 0) stop("rs429358 not in UKB-PPP pQTL!")
a <- e4$beta; a_se <- e4$se; a_p <- e4$pval
cat(sprintf("  path a (ε4→APOE, UKB-PPP): beta=%+.4f (SE=%.4f, p=%.1e)\n", a, a_se, a_p))

# ---------- path b: APOE → disease (R13 focal clumped MR) ----------
foc <- fread(file.path(PROJ, "revision/results/r13/apoe_focus_mr_clumped.csv"))
# Method col = "APOE_clumped_<disease> <method>"; Beta/SE hold the estimates
b_amd <- foc[grepl("APOE_clumped_AMD IVW", Method), ][1, .(Beta, SE)]
b_ad  <- foc[grepl("APOE_clumped_AD IVW", Method), ][1, .(Beta, SE)]
if (nrow(b_amd) == 0 || nrow(b_ad) == 0) stop("R13 focal IVW (AMD/AD) not found in focal file")
b_amd_v <- b_amd$Beta; b_amd_se <- b_amd$SE
b_ad_v  <- b_ad$Beta;  b_ad_se  <- b_ad$SE
cat(sprintf("  path b (APOE→AMD, R13 focal IVW): beta=%+.4f (SE=%.4f)\n", b_amd_v, b_amd_se))
cat(sprintf("  path b (APOE→AD,  R13 focal IVW): beta=%+.4f (SE=%.4f)\n", b_ad_v, b_ad_se))

# ---------- Total: ε4 → disease (R13 rs429358 per-allele) ----------
read_e4 <- function(f) {
  # AMD_r13.ma has a header line ("SNP A1 ..."); AD_r13.ma does not. Detect.
  first_line <- readLines(f, n = 1, warn = FALSE)
  has_header <- grepl("^SNP", first_line)
  d <- fread(f, header = has_header, data.table = TRUE)
  if (!has_header) setnames(d, 1:8, c("SNP","A1","A2","freq","b","se","p","n"))
  r <- d[SNP == "rs429358", ]
  if (nrow(r) == 0) stop("rs429358 not in ", f)
  b <- as.numeric(r$b[1]); s <- as.numeric(r$se[1])
  if (is.na(b) || is.na(s)) stop("b/se NA in ", f, " -> re-check format")
  list(b = b, se = s)
}
cat("\n  Reading R13 Total effect (rs429358 per-allele)...\n")
t_amd <- read_e4(file.path(PROJ, "smr_data/AMD_r13.ma"))
t_ad  <- read_e4(file.path(PROJ, "smr_data/AD_r13.ma"))
cat(sprintf("  Total (ε4→AMD, R13): beta=%+.4f (SE=%.4f)\n", t_amd$b, t_amd$se))
cat(sprintf("  Total (ε4→AD,  R13): beta=%+.4f (SE=%.4f)\n", t_ad$b, t_ad$se))

# ---------- mediation: product of coefficients + delta ----------
calc_mediation <- function(a, a_se, b, b_se, total, total_se, label) {
  indirect <- a * b
  indirect_se <- sqrt(b^2 * a_se^2 + a^2 * b_se^2)
  indirect_z <- indirect / indirect_se
  indirect_p <- 2 * pnorm(-abs(indirect_z))
  direct <- total - indirect
  direct_se <- sqrt(total_se^2 + indirect_se^2)
  direct_z <- direct / direct_se
  direct_p <- 2 * pnorm(-abs(direct_z))
  prop <- if (abs(total) > 1e-10) indirect / total * 100 else NA_real_
  cat(sprintf("\n  --- %s (R13) ---\n", label))
  cat(sprintf("    Total (ε4→disease, R13):  beta=%+.4f (SE=%.4f)\n", total, total_se))
  cat(sprintf("    path a (ε4→APOE):         beta=%+.4f (SE=%.4f)\n", a, a_se))
  cat(sprintf("    path b (APOE→disease):    beta=%+.4f (SE=%.4f)\n", b, b_se))
  cat(sprintf("    Indirect (a×b):           beta=%+.4f (SE=%.4f, z=%.2f, p=%.2e)\n", indirect, indirect_se, indirect_z, indirect_p))
  cat(sprintf("    Direct (Total-Indirect):  beta=%+.4f (SE=%.4f, z=%.2f, p=%.2e)\n", direct, direct_se, direct_z, direct_p))
  cat(sprintf("    Proportion reproduced:    %.1f%%  (scale caveat: see header)\n", prop))
  if (!is.na(prop) && prop > 100) cat("    [suppression / inconsistent mediation: direct opposes indirect]\n")
  data.table(Disease = label, Total = total, Total_SE = total_se,
             Indirect = indirect, Indirect_SE = indirect_se, Indirect_P = indirect_p,
             Sobel_Z = indirect_z, Sobel_P = indirect_p,
             Direct = direct, Direct_SE = direct_se, Direct_P = direct_p,
             Prop_Mediated = sprintf("%.1f%%", prop),
             Interpretation = ifelse(direct_p < 0.05,
               "Partial mediation via APOE level; significant residual isoform-specific (direct) effect",
               "APOE level accounts for full detectable effect; no significant residual (isoform-specific) direct pathway"))
}

results <- rbind(
  calc_mediation(a, a_se, b_amd_v, b_amd_se, t_amd$b, t_amd$se, "AMD_R13"),
  calc_mediation(a, a_se, b_ad_v, b_ad_se, t_ad$b, t_ad$se, "AD_R13"),
  fill = TRUE
)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
fwrite(results, file.path(OUT, "mediation_results_r13.csv"))

# ---------- compare with R8-era version ----------
cat("\n===== R8-era (old _17) for comparison =====\n")
old_file <- file.path(OUT, "mediation_results.csv")
if (file.exists(old_file)) {
  old <- fread(old_file)
  print(old[, .(Disease, Total, Indirect, Direct, Prop_Mediated, Interpretation)])
} else cat("  (old mediation_results.csv not found)\n")

cat("\n===== R13 re-estimation done =====\n")
cat("  output: apoe_results/mediation_results_r13.csv\n")
cat("  Scale caveat: proportion is strict only if a and b share APOE scale\n")
cat("  (both Z-score or both raw); verify before writing Methods.\n")
cat("\n--- sessionInfo ---\n"); sessionInfo()
cat("\nOK mediation R13.\n")
