#!/usr/bin/env Rscript
# =============================================================================
# build_fig1_retiered.R (v14, 2026-09-09)
# Figure 1 forest data under the finalised cohort architecture:
#   AMD: Discovery = FinnGen R13 | First validation = IAMDGC late | Second validation = MVP (META, APOE only)
#   AD : Discovery = FinnGen R13 | First validation = Jansen 2019  | Second validation = Bellenguez 2022 (APOE only)
# Proteins included = union of Bonferroni-significant hits (p_adj < 0.05) across that disease's cohorts.
# Output: figures/Figure1_forest_data.csv
# =============================================================================
suppressMessages({library(data.table)})

R13_AMD <- "pqtl_gsmr_results/pqtl_gsmr_v2_full.csv"
R13_AD  <- "pqtl_gsmr_results/pqtl_gsmr_v2_AD_full.csv"
IAM_AMD <- "pqtl_gsmr_results/iamdgc/pqtl_gsmr_v2_full.csv"
JAN_AD  <- "pqtl_gsmr_results/jansen/pqtl_gsmr_v2_full.csv"
MVP_MR  <- "revision/results/mvp_ancestry/mvp_mr_all.csv"

pick <- function(path, cohort, alpha = 0.05) {
  d <- fread(path); setnames(d, tolower(names(d)))
  d <- d[is.finite(as.numeric(bxy)) & is.finite(as.numeric(bxy_se))]
  d[, b := as.numeric(bxy)][, s := as.numeric(bxy_se)][, p := as.numeric(bxy_p)]
  d[, padj := as.numeric(p_adj)]
  d <- d[is.finite(padj) & padj < alpha]
  data.table(disease = NA_character_, gene = toupper(d$gene), cohort = cohort,
             or = exp(d$b), ci_lo = exp(d$b - 1.96 * d$s), ci_hi = exp(d$b + 1.96 * d$s),
             p = d$p, nsnp = as.integer(d$n_used))
}

amd <- rbindlist(list(
  pick(R13_AMD, "Discovery")[, disease := "AMD"],
  pick(IAM_AMD, "First validation")[, disease := "AMD"]), fill = TRUE)
ad <- rbindlist(list(
  pick(R13_AD, "Discovery")[, disease := "AD"],
  pick(JAN_AD, "First validation")[, disease := "AD"]), fill = TRUE)

# --- second-validation cohorts: APOE only ---
# FIX(2026-09-10): Second validation 两层原先用 IVW（MVP 取 META,IVW；Bellenguez 取 IVW(RE)），
# 与 Discovery / First validation 的 GSMR 不一致。现改为统一 GSMR —— 由
# revision/98_second_validation_gsmr.R 用同一套 11 个 clumped 工具 + EU1000QC LD 补跑：
#   Bellenguez 2022 (AD) GSMR beta=-0.7634932, se=0.02748117, P=7.08e-170, n_used=7（HEIDI 剔除 1）
#   MVP META      (AMD) GSMR beta=+0.0437596, se=0.01441239, P=2.395e-3,  n_used=10
SV_GSMR <- "revision/results/r13/second_validation_gsmr.csv"
sv <- fread(SV_GSMR)
get_sv <- function(tag) { r <- sv[cohort == tag]; stopifnot(nrow(r) == 1); r }
mk_row <- function(disease, r) data.table(disease = disease, gene = "APOE", cohort = "Second validation",
  or = exp(r$beta), ci_lo = exp(r$beta - 1.96 * r$se), ci_hi = exp(r$beta + 1.96 * r$se),
  p = r$p, nsnp = as.integer(r$n_used))
amd <- rbind(amd, mk_row("AMD", get_sv("MVP_META_AMD")), fill = TRUE)
ad  <- rbind(ad,  mk_row("AD",  get_sv("Bellenguez_2022_AD")), fill = TRUE)

d <- rbindlist(list(amd, ad), fill = TRUE)
# F13B (IAMDGC only, OR = 10.26) compresses the log axis; kept in Supplementary Table S1b,
# reported in the Figure 1 legend instead.
d <- d[!(gene == "F13B")]
setorder(d, disease, gene, cohort)
fwrite(d, "figures/Figure1_forest_data.csv")
cat("Figure 1 forest data (v14 tiering):", nrow(d), "rows\n")
print(d[, .(disease, gene, cohort, or = round(or, 3), lo = round(ci_lo, 3), hi = round(ci_hi, 3),
            p = signif(p, 3), n = nsnp)], digits = 4)
