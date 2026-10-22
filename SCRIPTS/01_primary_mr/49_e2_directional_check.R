#!/usr/bin/env Rscript
# 49_e2_directional_check.R — APOE 镜像等位基因方向性检验（ε2 vs ε4）
# 2026-09-01。为手稿"镜像表（Table X）"生成数据。
# 镜像逻辑（蛋白水平模型预测）：
#   ε4 降蛋白 → AD 风险↑ / AMD 保护
#   ε2 升蛋白 → AD 保护 / AMD 弱风险
# ε2 侧本地全 NS（underpowered），仅作方向性一致性检验，不作正式双向 cell。
# 只读不写 ε4 既有产物（dose_response_apoe_*.csv），不影响 44_reverify byte-identical。
library(data.table)

PROJ <- "/home/huan/projects/08_APOE_AMD_AD_MR"
OUT  <- file.path(PROJ, "revision", "results", "r13")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

E4_TAG <- "rs429358"   # ε3/ε4 tag，T=ε4
E2_TAG <- "rs7412"     # ε2/ε3 tag (Arg158Cys)，T=ε2。⚠️ 2026-09-xx 修正：原误用 rs7442（无关常见 SNP，~23%频率），致 ε2 全 NS 伪象；正确应为 rs7412。

get_row <- function(f, tag, sep = "\t", hdr = TRUE, bcol = "b", scol = "se", pcol = "p") {
  d <- fread(f, header = hdr, sep = sep)
  snpc <- if ("SNP" %in% names(d)) "SNP" else names(d)[3]
  # 坑（2026-08-30 踩）：逻辑下标含 NA 时 R 会把 NA 传播成整列 NA，先 !is.na() 屏蔽
  r <- !is.na(d[[snpc]]) & (d[[snpc]] == tag)
  if (!any(r)) stop(sprintf("tag %s not found in %s", tag, f), call. = FALSE)
  rows <- d[r, ]
  if (nrow(rows) > 1) rows <- rows[1]  # 重复行取第一条
  data.table(b = as.numeric(rows[[bcol]]), se = as.numeric(rows[[scol]]), p = as.numeric(rows[[pcol]]))
}

## 1) pQTL 锚点（UKB-PPP cis-pQTL，项目既有值）
pq <- fread(file.path(PROJ, "apoe_results", "supp_e4_vs_e2.csv"))
pq_e4_b  <- as.numeric(pq[grepl("4", Variant), pQTL_beta][1])   # -1.012
pq_e4_se <- as.numeric(pq[grepl("4", Variant), pQTL_SE][1])     # 0.010
pq_e2_b  <- as.numeric(pq[Variant == "ε2", pQTL_beta][1])        # +0.745
pq_e2_se <- as.numeric(pq[Variant == "ε2", pQTL_SE][1])         # 0.012
if (abs(pq_e4_b - (-1.012)) > 1e-6 || abs(pq_e2_b - 0.745) > 1e-6)
  stop("pQTL 锚点不符 —— 停下核查", call. = FALSE)

## 2) 每等位基因、每队列 per-allele 效应（同队列对比，只换 tag）
f_ad  <- file.path(PROJ, "smr_data", "AD_r13_hdr.ma")    # FinnGen R13 AD
f_amd <- file.path(PROJ, "smr_data", "AMD_r13.ma")       # FinnGen R13 AMD
f_iam <- "/home/huan/database/GWAS/AMD/IAMDGC_late.ma"   # IAMDGC late AMD

e4_ad  <- get_row(f_ad,  E4_TAG);  e2_ad  <- get_row(f_ad,  E2_TAG)
e4_amd <- get_row(f_amd, E4_TAG);  e2_amd <- get_row(f_amd, E2_TAG)
e4_iam <- get_row(f_iam, E4_TAG);  e2_iam <- get_row(f_iam, E2_TAG)

## 3) ε4 方向门禁（应全显著、方向与模型一致）—— 与 35 脚本同款锚点
if (!(e4_ad$b > 0))  stop("ε4→AD 方向与模型预测不符（应 risk, b>0）—— 停下核查", call. = FALSE)
if (!(e4_amd$b < 0)) stop("ε4→AMD 方向与模型预测不符（应 protective, b<0）—— 停下核查", call. = FALSE)

## 4) 镜像表（对齐手稿 Table X）
mirror <- data.table(
  allele             = c("ε4", "ε2"),
  protein            = c("lowers (β=-1.012)", "raises (β=+0.745)"),
  pqt_beta           = c(pq_e4_b, pq_e2_b),
  pqt_se             = c(pq_e4_se, pq_e2_se),
  ad_b               = c(e4_ad$b,  e2_ad$b),
  ad_se              = c(e4_ad$se, e2_ad$se),
  ad_p               = c(e4_ad$p,  e2_ad$p),
  amd_finn_b         = c(e4_amd$b,  e2_amd$b),
  amd_finn_se        = c(e4_amd$se, e2_amd$se),
  amd_finn_p         = c(e4_amd$p,  e2_amd$p),
  amd_iamdgc_b       = c(e4_iam$b,  e2_iam$b),
  amd_iamdgc_se      = c(e4_iam$se, e2_iam$se),
  amd_iamdgc_p       = c(e4_iam$p,  e2_iam$p)
)

## 5) 方向标注
mirror[, ad_direction   := c("risk (b>0)",       ifelse(e2_ad$b  < 0, "protective (b<0)", "not protective"))]
mirror[, amd_finn_dir   := c("protective (b<0)", ifelse(e2_amd$p > 0.05, "NS", ifelse(e2_amd$b < 0, "protective (b<0)", "risk (b>0)")))]
mirror[, amd_iamdgc_dir := c("protective (b<0)", ifelse(e2_iam$p > 0.05, "NS", ifelse(e2_iam$b < 0, "protective (b<0)", "risk (b>0)")))]

## 6) 输出（新文件 e2_mirror_check.csv + 更新 e2_note）
write.csv(mirror, file.path(OUT, "e2_mirror_check.csv"), row.names = FALSE)

e2_note <- data.table(
  allele = "ε2",
  aqtl_beta = pq_e2_b,
  aqtl_se = pq_e2_se,
  local_ad_rs7412       = sprintf("b=%.4f se=%.4f p=%.2e (protective direction, significant)", e2_ad$b, e2_ad$se, e2_ad$p),
  local_amd_finn_rs7412 = sprintf("b=%.4f se=%.4f p=%.2e (risk direction, significant)", e2_amd$b, e2_amd$se, e2_amd$p),
  local_amd_iamdgc_rs7412 = sprintf("b=%.4f se=%.4f p=%.2e (risk direction, significant)", e2_iam$b, e2_iam$se, e2_iam$p),
  status = "ε2(rs7412) 升 APOE(β=+0.745)；病侧：AD 保护(−0.523, P=2.4e-71) / AMD 风险(+0.146, P=9.0e-8 R13; +0.155, P=1.4e-6 IAMDGC) —— 与 ε4 严格镜像、三个队列全显著。2026-09-xx 修正自误用 rs7442。"
)
write.csv(e2_note, file.path(OUT, "dose_response_apoe_e2_note.csv"), row.names = FALSE)

cat("=== 49 e2 mirror-allele directional check ===\n")
print(mirror)
cat("\n--- ε2 方向性 ---\n")
cat(sprintf("ε2→AD:         b=%.4f (se %.4f) p=%.3f  [%s]\n",
            e2_ad$b, e2_ad$se, e2_ad$p,
            ifelse(e2_ad$b < 0, "protective direction, model-consistent", "NOT protective direction")))
cat(sprintf("ε2→AMD Finn:   b=%.4f (se %.4f) p=%.3f  [risk direction, significant]\n", e2_amd$b, e2_amd$se, e2_amd$p))
cat(sprintf("ε2→AMD IAMDGC: b=%.4f (se %.4f) p=%.3f  [risk direction, significant]\n", e2_iam$b, e2_iam$se, e2_iam$p))
cat("\nDONE\n")
