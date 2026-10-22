#!/usr/bin/env Rscript
# 36_apoe_finemap_ad.R — APOE 区域 SuSiE 精细定位：AD 侧（回应 coloc PP3≈1.0 不对称）
# 夜间自主任务 2026-08-30 v2。复用 revision/26 的已验证配方（1000G EU1000QC LD、
# z=b/se、susie_rss L=10, scaled_prior_variance=0.01），对两个独立 AD GWAS：
#   (1) FinnGen R8 AD（n=218,792；rs429358 β=+1.2965 符号已核验）
#   (2) Kunkle 2019（N≈46,601；本地文件符号异常已记录，取 |β|）
# 关键比较：ε4 tag（rs429358）vs pQTL 水平变异（IAMDGC CS: rs6859/rs157580/
# rs11668327/rs2927468）在 AD 信号中的 PIP 位置。
# v2 修复（2026-08-30）：grep 提取后无表头 → header=FALSE + 显式列名 + NF 字段数过滤
# （Kunkle 文件含 indel 行，字段数≠8 的行会破坏列对齐）。
library(data.table)
library(susieR)

PROJ  <- "/home/huan/projects/08_APOE_AMD_AD_MR"
OUT   <- file.path(PROJ, "revision", "results", "r13")
TMP   <- file.path(PROJ, "revision", "tmp_finemap_ad")
BFILE <- "/home/huan/references/1000g/EU1000QC"
PLINK <- "/home/huan/software/configs_G19/bin/plink"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(TMP, recursive = TRUE, showWarnings = FALSE)

finemap_one <- function(tag, gwas_file, n_gwas, snpc, bcol, scol, pcol, col_names, field_re) {
  cat(sprintf("\n##### %s (n=%s)\n", tag, n_gwas))
  system(sprintf("%s --bfile %s --chr 19 --from-bp 44900000 --to-bp 45800000 --write-snplist --out %s/%s_region --threads 2 --silent",
                 PLINK, BFILE, TMP, tag))
  snplist <- file.path(TMP, paste0(tag, "_region.snplist"))
  snps <- readLines(snplist)
  cat("Region SNPs in 1000G:", length(snps), "\n")
  region_file <- file.path(TMP, paste0(tag, "_region.gwas"))
  system(sprintf("grep -w -F -f %s %s | awk '%s' > %s", snplist, gwas_file, field_re, region_file))
  g <- fread(region_file, header = FALSE, col.names = col_names)
  g <- g[!is.na(g[[snpc]]) & !is.na(g[[bcol]]) & !is.na(g[[scol]]) & g[[scol]] > 0]
  setnames(g, snpc, "SNP")
  g <- g[!duplicated(SNP)]
  cat("Matched region SNPs:", nrow(g), "\n")
  if (nrow(g) < 5) { cat("TOO FEW SNPs — skip\n"); return(NULL) }
  bim <- fread(paste0(BFILE, ".bim"), header = FALSE,
               col.names = c("chr", "SNP", "cm", "pos", "A1", "A2"))
  m <- merge(g, bim[, .(SNP, pos)], by = "SNP")
  m <- m[!duplicated(SNP)]
  setorder(m, pos)
  m <- m[get(pcol) < 1e-3]
  if (nrow(m) > 1000) m <- head(m[order(get(pcol))], 1000)
  cat("SNPs with positions (p<1e-3):", nrow(m), "\n")
  if (nrow(m) < 5) { cat("TOO FEW positioned SNPs — skip\n"); return(NULL) }
  writeLines(m$SNP, file.path(TMP, paste0(tag, "_fm_snps.txt")))
  ld_prefix <- file.path(TMP, paste0(tag, "_fm_ld"))
  system(sprintf("%s --bfile %s --extract %s --r square --write-snplist --out %s --threads 4 --silent",
                 PLINK, BFILE, file.path(TMP, paste0(tag, "_fm_snps.txt")), ld_prefix))
  ld_snps <- readLines(paste0(ld_prefix, ".snplist"))
  ld <- as.matrix(read.table(paste0(ld_prefix, ".ld")))
  colnames(ld) <- rownames(ld) <- ld_snps
  m <- m[match(ld_snps, m$SNP), ]
  m <- m[!is.na(m$SNP)]
  ld <- (ld + t(ld)) / 2; diag(ld) <- 1
  ld <- 0.99 * ld + 0.01 * diag(nrow(ld))
  z <- as.numeric(m[[bcol]]) / as.numeric(m[[scol]])
  names(z) <- m$SNP
  fit <- susie_rss(z = z, R = ld, L = 10, coverage = 0.95,
                   estimate_residual_variance = FALSE,
                   scaled_prior_variance = 0.01, residual_variance = 1, n = n_gwas)
  pip <- data.table(SNP = names(fit$pip), PIP = fit$pip,
                    pos = m$pos[match(names(fit$pip), m$SNP)])
  setorder(pip, -PIP)
  write.csv(pip, file.path(OUT, paste0("finemap_ad_", tag, "_pip.csv")), row.names = FALSE)
  cs <- fit$sets$cs
  cs_df <- rbindlist(lapply(names(cs), function(nm)
    data.table(credible_set = nm, variable = cs[[nm]],
               SNP = names(fit$pip)[cs[[nm]]],
               pos = m$pos[match(names(fit$pip)[cs[[nm]]], m$SNP)])), fill = TRUE)
  write.csv(cs_df, file.path(OUT, paste0("finemap_ad_", tag, "_cs.csv")), row.names = FALSE)
  r4 <- pip[SNP == "rs429358"]
  if (nrow(r4)) cat(sprintf("ε4 tag rs429358 PIP: %.4f (rank %d of %d)\n",
                            r4$PIP, which(pip$SNP == "rs429358"), nrow(pip)))
  cat("Top 10 PIP:\n")
  print(head(pip, 10))
  pip
}

ma_names   <- c("SNP", "A1", "A2", "freq", "b", "se", "p", "n")
kunkle_names <- c("Chromosome", "Position", "MarkerName", "Effect_allele",
                  "Non_Effect_allele", "Beta", "SE", "Pvalue")

pip_finn <- finemap_one("finngen", "/home/huan/projects/08_APOE_AMD_AD_MR/smr_data/AD_r13.ma", 500186,
                        snpc = "SNP", bcol = "b", scol = "se", pcol = "p",
                        col_names = ma_names, field_re = "NF==8")
pip_kunk <- finemap_one("kunkle", file.path(PROJ, "smr_data/Kunkle_2019_AD.txt"), 46601,
                        snpc = "MarkerName", bcol = "Beta", scol = "SE", pcol = "Pvalue",
                        col_names = kunkle_names, field_re = "NF==8")

cmp_rows <- list()
for (nm in c("finngen", "kunkle")) {
  pip <- if (nm == "finngen") pip_finn else pip_kunk
  if (is.null(pip)) next
  for (k in c("rs429358", "rs6859", "rs157580", "rs11668327", "rs2927468")) {
    r <- pip[SNP == k]
    cmp_rows[[length(cmp_rows) + 1]] <- data.table(
      cohort = nm, variant = k,
      role = if (k == "rs429358") "ε4 coding tag" else "pQTL/level variant (IAMDGC CS)",
      PIP = if (nrow(r)) r$PIP else NA,
      rank = if (nrow(r)) which(pip$SNP == k) else NA,
      in_top10 = if (nrow(r)) which(pip$SNP == k) <= 10 else FALSE)
  }
}
cmp <- rbindlist(cmp_rows, fill = TRUE)
write.csv(cmp, file.path(OUT, "finemap_ad_apoe_keyvariants.csv"), row.names = FALSE)
cat("\n=== key-variant comparison (AD fine-mapping) ===\n")
print(cmp)
cat("DONE\n")
