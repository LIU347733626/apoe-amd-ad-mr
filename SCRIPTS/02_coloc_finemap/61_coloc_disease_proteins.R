#!/usr/bin/env Rscript
# 61: disease-significant-protein cis-pQTL x {AMD, AD} colocalization (coloc.abf)
#
# Purpose: the manuscript so far only ran coloc for APOE (revision/23, PP4=0.9998).
#          Here the 13 IAMDGC late-AMD Bonferroni proteins and the 5 Jansen AD
#          proteins (17 unique) are tested one by one against the two R13
#          discovery GWAS, to ask which of them share a causal variant with the
#          disease signal (like APOE does) and which do not.
#
# Inputs
#   D1: UKB-PPP cis-pQTL, pre-extracted in one awk pass to
#       revision/tmp_coloc_prot/cis_raw.tsv (gene, SNP, other_allele,
#       effect_allele, eaf, beta, se, pval, samplesize, pos_hg38, gene_start,
#       gene_end, cis_trans).  cis window re-applied here: gene +/- 1 Mb.
#   D2: smr_data/AMD_r13.ma (header) / smr_data/AD_r13.ma (no header)
#       columns SNP A1 A2 freq b se p n
#
# Conventions reused from revision/23_pqtl_coloc_reproducible.R and
# revision/60_coloc_bydisease.R:
#   - pQTL beta flipped onto the outcome effect allele (beta_harm); coloc's ABF
#     is sign-agnostic, harmonization is for direction bookkeeping only.
#   - outcome p recomputed from b/se (robust to p=0 clamp in converted files).
#   - MAF = min(eaf, 1-eaf); D1 type="quant", sdY=1, N=samplesize;
#     D2 type="cc", s = cases/N.
#   - priors: p1=p2=1e-4, p12=1e-5 (primary); p12=1e-6 (sensitivity).
#   - case fractions: AMD 13947/472750, AD 16652/500186.
#
# CAVEAT (documented in the output): the local UKB-PPP file is a p<5e-8
#   *unclumped* subset, not full regional summary statistics.  coloc.abf assumes
#   the complete set of regional variants, so PP.H0/PP.H1 are structurally
#   inflated and PP3 (distinct causal variants) is hard to reach.  The verdict
#   column therefore distinguishes "shared", "distinct" and "not evaluable".
#
# Output: revision/results/mechanism/coloc_disease_proteins.csv
##############################################################################
suppressPackageStartupMessages({
  library(data.table); library(coloc); library(parallel)
})

PROJ <- "/home/huan/projects/08_APOE_AMD_AD_MR"
CIS  <- file.path(PROJ, "revision/tmp_coloc_prot/cis_raw.tsv")
OUT  <- file.path(PROJ, "revision/results/mechanism")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

AMD_F <- file.path(PROJ, "smr_data/AMD_r13.ma")
AD_F  <- file.path(PROJ, "smr_data/AD_r13.ma")
AMD_S <- 13947 / 472750     # 0.029501
AD_S  <- 16652 / 500186     # 0.033292

cat("===== 61: disease-significant-protein coloc =====\n")
cat(sprintf("coloc %s | data.table %s | R %s\n",
            packageVersion("coloc"), packageVersion("data.table"),
            paste(R.version$major, R.version$minor, sep = ".")))

## ---------- D1: cis-pQTL, one awk pass output ----------
d <- fread(CIS, header = FALSE)
setnames(d, c("gene","SNP","other_allele","effect_allele","eaf","beta","se",
              "pval","samplesize","pos_hg38","gene_start","gene_end","cis_trans"))
num <- c("eaf","beta","se","pval","samplesize","pos_hg38","gene_start","gene_end")
d[, (num) := lapply(.SD, as.numeric), .SDcols = num]
d <- d[!is.na(SNP) & !is.na(se) & se > 0 & !is.na(eaf) & eaf > 0 & eaf < 1 &
         !is.na(pos_hg38) & !is.na(beta)]
## NOTE 2026-09-10: do NOT filter on pval > 0.  49/1407 APOE cis rows have the
## UKB-PPP p-value underflowed to exactly 0 (the strongest pQTL signals); a
## `pval > 0` filter silently removed them and flipped the APOE golden anchor
## from PP4 ~ 1.0 to PP4 = 3.8e-4.  p is recomputed from beta/se below.
d[, p_z := 2 * pnorm(-abs(beta / se))]
n0 <- nrow(d)
d <- d[pos_hg38 >= gene_start - 1e6 & pos_hg38 <= gene_end + 1e6]
cat(sprintf("cis rows: %d -> %d after gene +/-1Mb window\n", n0, nrow(d)))
d <- unique(d, by = c("gene","SNP"))
cat(sprintf("unique gene-SNP pairs: %d | genes: %d\n", nrow(d), uniqueN(d$gene)))
print(d[, .(n_snp = .N, n_ss = uniqueN(samplesize),
            samplesize = paste(sort(unique(samplesize)), collapse = ";")), by = gene][order(gene)])

## ---------- D2: outcome rows for the union of cis SNPs ----------
read_out <- function(path, label) {
  snps <- unique(d$SNP)
  pat  <- tempfile(); writeLines(snps, pat)
  cat(sprintf("[%s] grepping %d SNP ids ...\n", label, length(snps)))
  x <- fread(cmd = sprintf("grep -w -F -f %s %s", pat, path),
             header = FALSE, sep = "\t",
             colClasses = c("character","character","character","numeric",
                            "numeric","numeric","numeric","integer"),
             col.names = c("SNP","o_A1","o_A2","o_freq","o_b","o_se","o_p","o_n"))
  file.remove(pat)
  x <- x[!is.na(o_b) & !is.na(o_se) & o_se > 0 & !is.na(o_freq) &
           o_freq > 0 & o_freq < 1]
  x <- unique(x, by = "SNP")
  cat(sprintf("[%s] outcome rows: %d\n", label, nrow(x)))
  x
}
amd <- read_out(AMD_F, "AMD_r13")
ad  <- read_out(AD_F,  "AD_r13")

## ---------- coloc runner ----------
run_one <- function(prot, disease, p12 = 1e-5) {
  # NOTE 2026-09-10: argument MUST NOT be called `gene` - inside a data.table
  # `i` expression the column `gene` shadows the function argument, so
  # d[gene == gene] is always TRUE and every protein silently uses all 103k
  # cis rows (caught by the APOE golden anchor gate below).
  g  <- d[gene == prot]
  od <- if (disease == "AMD") amd else ad
  s_case <- if (disease == "AMD") AMD_S else AD_S
  m <- merge(g, od, by = "SNP")   # outcome cols are o_* -> no suffix collisions
  if (nrow(m) < 5) return(NULL)
  # harmonize pQTL beta onto the outcome effect allele
  m[, beta_harm := fifelse(toupper(effect_allele) == toupper(o_A1) &
                             toupper(other_allele) == toupper(o_A2), beta,
                           fifelse(toupper(effect_allele) == toupper(o_A2) &
                                     toupper(other_allele) == toupper(o_A1), -beta,
                                   NA_real_))]
  m <- m[!is.na(beta_harm) & !duplicated(SNP)]
  m[, maf := pmin(eaf, 1 - eaf)]
  m[, maf_out := pmin(o_freq, 1 - o_freq)]
  m <- m[maf > 0 & maf <= 0.5 & maf_out > 0 & maf_out <= 0.5]
  if (nrow(m) < 5) return(NULL)
  N_pqtl <- as.numeric(names(sort(table(m$samplesize), decreasing = TRUE))[1])
  N_out  <- as.numeric(names(sort(table(m$o_n), decreasing = TRUE))[1])
  p_out  <- 2 * pnorm(-abs(m$o_b / m$o_se))

  D1 <- list(type = "quant", N = N_pqtl, sdY = 1, snp = m$SNP,
             pvalues = m$p_z, beta = m$beta_harm, varbeta = m$se^2,
             MAF = m$maf, position = m$pos_hg38)
  D2 <- list(type = "cc", N = N_out, s = s_case, snp = m$SNP, pvalues = p_out,
             beta = m$o_b, varbeta = m$o_se^2, MAF = m$maf_out,
             position = m$pos_hg38)
  fit <- try(coloc.abf(D1, D2, p1 = 1e-4, p2 = 1e-4, p12 = p12), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  s <- fit$summary
  lead <- m[which.max(abs(beta / se))]   # pval can underflow to 0; use the z
  lead_out <- m[which.max(abs(o_b / o_se))]  # same: min(p) is ambiguous at p=0
  data.table(
    protein = prot, disease = disease, n_snp = nrow(m),
    PP.H0 = as.numeric(s["PP.H0.abf"]), PP.H1 = as.numeric(s["PP.H1.abf"]),
    PP.H2 = as.numeric(s["PP.H2.abf"]), PP3 = as.numeric(s["PP.H3.abf"]),
    PP4 = as.numeric(s["PP.H4.abf"]), p12 = p12,
    N_pqtl = N_pqtl, s_outcome = s_case,
    n_pqtl_snp_p05_in_disease = sum(m$o_p < 0.05),
    n_pqtl_snp_gws_in_disease = sum(m$o_p < 5e-8),
    min_p_disease_at_pqtl_snp = min(m$o_p),
    lead_pqtl_snp = lead$SNP, lead_pqtl_p = lead$p_z,
    lead_pqtl_p_disease = lead$o_p,
    lead_pqtl_beta_disease = lead$o_b,
    lead_disease_snp = lead_out$SNP, lead_disease_p = lead_out$o_p,
    same_lead_snp = (lead$SNP == lead_out$SNP))
}

## hard length guard (2026-09-10: caught the merge-suffix bug that silently
## emptied p_out because the outcome se was renamed `se_out` by merge)
chk <- run_one("APOE", "AMD", 1e-5)
if (is.null(chk)) stop("APOE anchor returned NULL - pipeline broken", call. = FALSE)
cat(sprintf("APOE/AMD anchor: n=%d PP3=%.4g PP4=%.6g\n",
            chk$n_snp, chk$PP3, chk$PP4))
## golden anchor gate: revision/23 and revision/60 both give APOE PP4 ~ 1.0
## with ~1.2-1.4k harmonized SNPs.  If n_snp is 10x larger the protein filter
## has failed (column/argument shadowing) - stop rather than emit nonsense.
if (chk$n_snp > 3000 || chk$PP4 < 0.99)
  stop(sprintf("APOE golden anchor FAILED: n_snp=%d PP4=%.4g (expected ~1.2-1.4k SNPs, PP4>0.99)",
               chk$n_snp, chk$PP4), call. = FALSE)

grid <- CJ(gene = unique(d$gene), disease = c("AMD","AD"), sorted = FALSE)
res  <- rbindlist(mclapply(seq_len(nrow(grid)), function(i)
  run_one(grid$gene[i], grid$disease[i], 1e-5), mc.cores = min(nrow(grid), 34)),
  fill = TRUE)
sens <- rbindlist(mclapply(seq_len(nrow(grid)), function(i)
  run_one(grid$gene[i], grid$disease[i], 1e-6), mc.cores = min(nrow(grid), 34)),
  fill = TRUE)
if (nrow(res) == 0) stop("primary coloc pass returned no rows", call. = FALSE)
if (nrow(sens) > 0 && all(c("PP3","PP4") %in% names(sens))) {
  setnames(sens, c("PP3","PP4"), c("PP3_p12_1e6","PP4_p12_1e6"))
  res <- merge(res, sens[, .(protein, disease, PP3_p12_1e6, PP4_p12_1e6)],
               by = c("protein","disease"), all.x = TRUE)
} else {
  res[, `:=`(PP3_p12_1e6 = NA_real_, PP4_p12_1e6 = NA_real_)]
}

## ---------- LD between the two leads (EU1000QC) ----------
## A PP3-high call with the two leads in near-perfect LD cannot be read as
## "distinct causal variants"; a PP4-high call with the two leads identical is
## the strongest form of "shared".  plink --r2 labelled pair output (the
## --r2 square matrix of this plink build has no row/column labels).
LDREF <- "/home/huan/references/1000g/EU1000QC"
ldprefix <- file.path(PROJ, "revision/tmp_coloc_prot/lead_ld_pairs")
snps <- unique(c(res$lead_pqtl_snp, res$lead_disease_snp))
snps <- snps[!is.na(snps)]
writeLines(snps, paste0(ldprefix, ".snplist"))
cmd <- sprintf("plink --bfile %s --extract %s --r2 --ld-window-kb 5000 --ld-window 99999 --ld-window-r2 0 --out %s 2>&1 | tail -1",
               LDREF, paste0(ldprefix, ".snplist"), ldprefix)
cat("\n[LD]", cmd, "\n"); system(cmd)
ldf <- paste0(ldprefix, ".ld")
if (file.exists(ldf) && file.info(ldf)$size > 0) {
  lp <- fread(ldf, header = TRUE, sep = " ", strip.white = TRUE)
  setnames(lp, trimws(names(lp)))
  lp[, `:=`(SNP_A = trimws(SNP_A), SNP_B = trimws(SNP_B))]
  key <- lp[, .(pair = paste(pmin(SNP_A, SNP_B), pmax(SNP_A, SNP_B), sep = "|"), R2)]
  setkey(key, pair)
  get_r2 <- function(a, b) {
    if (is.na(a) || is.na(b)) return(NA_real_)
    if (a == b) return(1)
    key[.(paste(pmin(a, b), pmax(a, b), sep = "|"))]$R2[1]
  }
  res[, r2_lead_pqtl_vs_lead_disease := mapply(get_r2, lead_pqtl_snp, lead_disease_snp)]
  res[, n_leads_in_ldref := uniqueN(c(lead_pqtl_snp, lead_disease_snp)) -
        sum(is.na(match(unique(c(lead_pqtl_snp, lead_disease_snp)),
                        unique(c(lp$SNP_A, lp$SNP_B)))))]
} else {
  res[, `:=`(r2_lead_pqtl_vs_lead_disease = NA_real_, n_leads_in_ldref = 0L)]
}

## ---------- verdict ----------
## PP4 >= 0.8  : shared causal variant (like APOE)
## PP3 >= 0.8  : distinct causal variants
## PP3 >= 0.8 but the two leads are the same SNP or r2 >= 0.8 : the SNP set
##               cannot separate the two signals -> not interpretable
## PP4 and PP3 both < 0.5 and no disease signal at the pQTL SNPs : not evaluable
res[, verdict := fifelse(PP4 >= 0.8, "shared_causal_variant",
                  fifelse(PP3 >= 0.8 & (same_lead_snp | (!is.na(r2_lead_pqtl_vs_lead_disease) &
                                                            r2_lead_pqtl_vs_lead_disease >= 0.8)),
                          "PP3_high_but_leads_in_high_LD_not_separable",
                  fifelse(PP3 >= 0.8, "distinct_causal_variants",
                  fifelse(PP4 < 0.5 & PP3 < 0.5 & min_p_disease_at_pqtl_snp > 0.05,
                          "not_evaluable_no_disease_signal_at_pqtl_snps",
                          "inconclusive"))))]
setorder(res, disease, -PP4)
## how strong is the disease-side evidence at the pQTL lead SNP?
res[, shared_support := fifelse(PP4 < 0.8, NA_character_,
                          fifelse(lead_pqtl_p_disease < 5e-8,
                                  "disease signal at pQTL lead is genome-wide significant",
                                  "disease signal at pQTL lead is only nominally/weak"))]
fwrite(res, file.path(OUT, "coloc_disease_proteins.csv"))
print(res[, .(protein, disease, n_snp, PP3 = signif(PP3, 3), PP4 = signif(PP4, 3),
              n_gws = n_pqtl_snp_gws_in_disease,
              lead = lead_pqtl_snp, lead_disease_snp, same_lead_snp,
              r2 = signif(r2_lead_pqtl_vs_lead_disease, 3), verdict)])

## ---------- golden anchor: APOE must reproduce revision/23 / 60 ----------
cat("\n=== APOE anchor (expect PP4 ~ 1.0 vs AMD_r13 and AD_r13) ===\n")
print(res[protein == "APOE", .(disease, n_snp, PP3, PP4)])
cat("\nWritten:", file.path(OUT, "coloc_disease_proteins.csv"), "\n")
cat("Done.\n")
