#!/usr/bin/env Rscript
# =============================================================================
# Figure 5 (机制图) — APOE 机制：表达在哪里 / 哪些细胞 / 下游为何分岔 / 有多特异
#   **实际面板（2026-09-10 核）= A–G + H1–H2，顺序即渲染顺序：**
#     A  视网膜表达图谱（疾病蛋白 + APOE 受体轴，基因 × 细胞类型点图）      [对象 pA]
#     B  APOE 受体轴按细胞类型分岔                                        [对象 pE]
#     C  APOE → TREM2–TYROBP（CellChat 通讯概率，唯一数据库注释对）        [对象 pF]
#     D  疾病相关细胞类型：**左=scDRS rank-biserial；右=seismic −log10 FDR**（两联，同属 D）[对象 pB = mk_sc()|mk_se()]
#     E  补体级联：蛋白层 GSEA（AMD/AD）与各细胞类型共表达 NES 合并展示     [对象 pC]
#     F  跨疾病蛋白网络（675 蛋白 GSMR，仅 APOE 与 C2 双侧显著反向）        [对象 pD]
#     G  APOE 共表达伙伴（每细胞类型 top6，官方 ppcor 偏相关）              [对象 pG]
#     H1 scDRS 去掉 APOE 区后的稳健性（去循环性检查）                       [对象 pHH1]
#     H2 区域局部遗传相关（官方 LAVA rho + chr19 同尺寸窗口分布）           [对象 pHH2]
#   ⚠ 对象名与面板字母**不对应**（pE 画的是 B、pC 画的是 E…），改动前先看上面这张对照表。
#   ⚠ 虚拟敲除（scTenifoldKnk）**不在本图**：其结果在正文 + Supplementary Table S32。
# =============================================================================
suppressMessages({ library(ggplot2); library(data.table); library(patchwork) })
PROJ <- "/home/huan/projects/08_APOE_AMD_AD_MR"
MD   <- file.path(PROJ, "revision/results/mechanism")
col_amd <- "#D1403F"; col_ad <- "#1F6FB2"; col_neu <- "#6E6E6E"
theme_pub <- theme_classic(base_size = 7.2, base_family = "Arial") +
  theme(panel.grid = element_blank(),
        axis.line = element_line(linewidth = 0.4, colour = "black"),
        axis.ticks = element_line(linewidth = 0.4, colour = "black"),
        axis.text = element_text(size = 6.4, colour = "#333333"),
        axis.title = element_text(size = 6.8, colour = "#222222"),
        plot.title = element_text(size = 7.6, face = "bold", colour = "#222222", hjust = 0),
        plot.margin = margin(2, 6, 2, 4, "mm"))

CT_ORDER <- c("Muller Glia", "RPE", "Endothelial", "Microglia", "Rod", "Unassigned")

# ---------------------------------------------------------------- A: 表达图谱
expr  <- fread(file.path(MD, "celltype_expression_map.csv"))
recep <- fread(file.path(MD, "apoe_receptor_axis.csv"))
expr[, group := "Disease-associated proteins"]
recep[, group := "APOE receptor / lipid axis"]
recep <- recep[cell_type != "NOT_IN_MATRIX"]
ap <- rbindlist(list(expr[, .(gene, cell_type, pct_expressing, group)],
                     recep[, .(gene, cell_type, pct_expressing, group)]), fill = TRUE)
ap[, pct := as.numeric(pct_expressing)]
# 只保留有表达的基因（至少一个细胞类型 >1%）
ap <- ap[, .SD[any(pct > 1)], by = gene]
# 每个基因只保留一个分组标签（受体轴优先），但保留其全部细胞类型行
grp <- unique(ap[, .(gene, group)])
grp <- grp[order(group == "Disease-associated proteins")]
grp <- grp[!duplicated(gene)]
ap <- merge(ap[, .(gene, cell_type, pct)], grp, by = "gene")
GENE_ORD <- rev(unique(c(
  "APOE", "CFH", "CFD", "CFI", "CFB", "CFHR2", "CFHR4", "CFHR5", "C2", "C3",
  "ACE", "EPHA1", "PILRA", "TREM2",
  "PLAUR", "AXL", "MERTK", "TYRO3", "LRP1", "LDLR", "LRP8", "LRP2",
  "SORL1", "ABCA1", "ABCG1", "SCARB1", "VLDLR", "ABCB1")))
ord <- ap[, .(mx = max(pct)), by = gene][order(mx)]
ap[, gene := factor(gene, levels = ord$gene)]
ap <- ap[!is.na(gene)]
ap[, cell_type := factor(cell_type, levels = CT_ORDER)]
pA <- ggplot(ap, aes(x = cell_type, y = gene)) +
  geom_point(aes(size = pct, colour = group)) +
  scale_size_continuous(range = c(0.3, 4.2), trans = "sqrt", name = "% cells", breaks = c(5,20,50,80)) +
  scale_colour_manual(values = c("Disease-associated proteins" = "#2E6E6E",
                                 "APOE receptor / lipid axis" = "#B5651D"), name = NULL) +
  scale_x_discrete(labels = c("Muller Glia" = "Müller", "RPE" = "RPE", "Endothelial" = "Endo",
                              "Microglia" = "Micro", "Rod" = "Rod", "Unassigned" = "Unass.")) +
  labs(x = NULL, y = NULL) +
  theme_pub + theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 6.0),
                    axis.text.y = element_text(size = 5.8),
                    legend.position = "bottom", legend.box = "horizontal",
                    legend.margin = margin(0, 0, 0, 0),
                    legend.text = element_text(size = 5.4),
                    legend.title = element_text(size = 5.4)) +
  ggtitle("A  Retinal expression: disease proteins and the APOE receptor axis")

# ---------------------------------------------------------------- B: 细胞类型
sc <- fread(file.path(PROJ, "NC_package/supplementary/TableS10a_scdrs_celltype.csv"))
sc[, disease := ifelse(grepl("AD", Disease), "AD", "AMD")]
se_a <- fread(file.path(PROJ, "revision/results/seismic_apoe/seismic_amd_ct_associations.csv"))
se_d <- fread(file.path(PROJ, "revision/results/seismic_apoe/seismic_ad_ct_associations.csv"))
se_a[, disease := "AMD"]; se_d[, disease := "AD"]
se <- rbindlist(list(se_a, se_d))
se[, neglog10fdr := -log10(FDR)]
mk_sc <- function() {
  sc[, cell_type := factor(Cell_type, levels = CT_ORDER)]
  ggplot(sc, aes(x = cell_type, y = rank_biserial, fill = disease)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.66) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "#999999") +
    scale_fill_manual(values = c("AMD" = col_amd, "AD" = col_ad), name = NULL) +
    scale_x_discrete(labels = c("Muller Glia"="Müller","RPE"="RPE","Endothelial"="Endo",
                                "Microglia"="Micro","Rod"="Rod","Unassigned"="Unass.")) +
    labs(x = NULL, y = "scDRS rank-biserial") +
    theme_pub + theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 6.0),
                      legend.position = "top", legend.text = element_text(size = 5.8),
                      legend.margin = margin(-4,0,-2,0)) +
    ggtitle("D  Disease-relevant cell types (scDRS)")
}
mk_se <- function() {
  se[, cell_type := factor(cell_type, levels = CT_ORDER)]
  ggplot(se, aes(x = cell_type, y = neglog10fdr, fill = disease)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.66) +
    geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "#999999", linewidth = 0.4) +
    scale_fill_manual(values = c("AMD" = col_amd, "AD" = col_ad), name = NULL) +
    scale_x_discrete(labels = c("Muller Glia"="Müller","RPE"="RPE","Endothelial"="Endo",
                                "Microglia"="Micro","Rod"="Rod","Unassigned"="Unass.")) +
    labs(x = NULL, y = "seismic \u2212log10 FDR") +
    theme_pub + theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 6.0),
                      legend.position = "none") +
    ggtitle("D  (continued) seismic \u2212log10 FDR")
}
# D 行的字母标签放在左侧子面板标题上（plot_annotation 的标题会被 patchwork 的全局 theme 覆盖）
pB <- (mk_sc() | mk_se())

# ---------------------------------------------------------------- C: 下游分岔
fa <- fread(file.path(MD, "fgsea_AMD.csv")); fd <- fread(file.path(MD, "fgsea_AD.csv"))
cp <- "KEGG_COMPLEMENT_AND_COAGULATION_CASCADES"
nes <- rbindlist(list(
  data.table(level = "Protein-level GSEA", disease = "AMD", NES = fa[pathway == cp, NES]),
  data.table(level = "Protein-level GSEA", disease = "AD",  NES = fd[pathway == cp, NES])))
co <- fread(file.path(MD, "apoe_coexpression_pathways_ppcor_allgenes.csv"))   # 2026-09-10: 官方 ppcor 排序（见 docs/audit_20260909/PPCOR_OFFICIAL_RECOMPUTE_20260910.md）
co <- co[pathway == cp]
nes2 <- data.table(level = paste0("Co-expression\n", co$cell_type), disease = "cell",
                   NES = co$NES, padj = co$padj)
nes[, padj := NA_real_]
allnes <- rbindlist(list(nes, nes2), fill = TRUE)
co[, ct_short := c("Endothelial"="Endo", "Rod"="Rod", "RPE"="RPE", "Microglia"="Micro",
                     "Muller_Glia"="M\u00fcller", "Unassigned"="Unass.", "Endothelial "="Endo")[cell_type]]
nes2 <- data.table(level = co$ct_short, disease = "cell", NES = co$NES, padj = co$padj)
nes[, padj := NA_real_]
allnes <- rbindlist(list(nes, nes2), fill = TRUE)
allnes[, lab := ifelse(level == "Protein-level GSEA", paste0("GSEA \u00b7 ", disease), level)]
allnes[, fillg := ifelse(NES > 0, "positive", "negative")]
pC <- ggplot(allnes, aes(x = reorder(lab, NES), y = NES, fill = fillg)) +
  geom_col(width = 0.62) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "#999999") +
  coord_flip() +
  scale_fill_manual(values = c("positive" = "#C0504D", "negative" = "#4F81BD"), guide = "none") +
  labs(x = NULL, y = "Complement-cascade NES") +
  theme_pub + theme(axis.text.y = element_text(size = 5.9)) +
  ggtitle("E  Complement cascade")

# ---------------------------------------------------------------- D: 特异性
cc <- fread(file.path(MD, "cocausal_scatter.csv"))
cc <- cc[is.finite(b_amd) & is.finite(b_ad)]
cc[, sig_amd := p_amd < 0.05]; cc[, sig_ad := p_ad < 0.05]  # 扫描级 FDR（16 AMD / 6 AD）
cc[, cat := fcase(sig_amd & sig_ad, "both", sig_amd, "AMD only", sig_ad, "AD only", default = "neither")]
pD <- ggplot(cc, aes(x = b_amd, y = b_ad)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "#BBBBBB", linewidth = 0.35) +
  geom_vline(xintercept = 0, linetype = 2, colour = "#BBBBBB", linewidth = 0.35) +
  geom_point(aes(colour = cat), size = 0.9, alpha = 0.75) +
  geom_point(data = cc[cat == "both"], size = 2.0, colour = "#111111") +
  geom_text(data = cc[cat == "both"], aes(label = gene), vjust = -1.1, size = 2.2, colour = "#111111") +
  scale_colour_manual(values = c("both" = "#C0504D", "AMD only" = col_amd,
                                 "AD only" = col_ad, "neither" = "#C8C8C8"), name = NULL) +
  labs(x = "GSMR \u03b2 per unit APOE \u2192 AMD", y = "GSMR \u03b2 \u2192 AD") +
  theme_pub + theme(legend.position = "top", legend.text = element_text(size = 5.8),
                    legend.margin = margin(-4,0,-2,0)) +
  ggtitle("F  Cross-disease proteome")


# ---------------------------------------------------------------- E: 受体轴分岔（独立面板，2026-09-10 从 A 拆出）
rax <- fread(file.path(MD, "apoe_receptor_axis.csv"))
RAX_GENES <- c("TREM2","TYROBP","PLAUR","AXL","MERTK","TYRO3",
               "LRP1","LDLR","LRP8","SCARB1","ABCA1","ABCG1","VLDLR","ABCB1")
RAX_GROUP <- c(TREM2="microglial", TYROBP="microglial", PLAUR="microglial", AXL="microglial",
               MERTK="microglial", TYRO3="microglial",
               LRP1="lipid/RPE-endothelial", LDLR="lipid/RPE-endothelial", LRP8="lipid/RPE-endothelial",
               SCARB1="lipid/RPE-endothelial", ABCA1="lipid/RPE-endothelial", ABCG1="lipid/RPE-endothelial",
               VLDLR="lipid/RPE-endothelial", ABCB1="lipid/RPE-endothelial")
rax <- rax[gene %in% RAX_GENES & cell_type %in% CT_ORDER]
rax[, grp := RAX_GROUP[gene]]
rax[, gene := factor(gene, levels = rev(RAX_GENES))]
pE <- ggplot(rax, aes(x = factor(cell_type, levels = CT_ORDER), y = gene)) +
  geom_point(aes(size = pct_expressing, colour = grp)) +
  scale_size_continuous(range = c(0.4, 4.2), name = "% cells", breaks = c(10, 30, 50)) +
  scale_colour_manual(values = c("microglial" = col_ad, "lipid/RPE-endothelial" = col_amd), name = NULL,
                      labels = c("microglial axis", "lipid/transport axis")) +
  labs(x = NULL, y = NULL) +
  ggtitle("B  APOE receptor axis by cell type") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 6.0),
        axis.text.y = element_text(size = 6.0, face = "italic"),
        legend.position = "bottom", legend.key.size = unit(2.6, "mm"),
        legend.text = element_text(size = 5.6), legend.title = element_text(size = 5.6),
        legend.box = "vertical", plot.margin = margin(3, 4, 3, 2, "mm"))

# ---------------------------------------------------------------- F: CellChat APOE -> TREM2-TYROBP 通讯
ccf <- fread(file.path(MD, "cellchat_APOE_TREM2_TYROBP.csv"))
ccf <- ccf[V1 %in% CT_ORDER & V2 %in% CT_ORDER]
ccf[, `:=`(V1 = factor(V1, levels = rev(CT_ORDER)), V2 = factor(V2, levels = CT_ORDER))]
pF <- ggplot(ccf, aes(x = V2, y = V1, fill = N)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(N > 0, sprintf("%.3f", N), "")), size = 1.9, colour = "#222222") +
  scale_fill_gradient(low = "#F2F2F2", high = "#7B4FA8", name = "comm. prob.") +
  labs(x = "target cell type", y = "source cell type") +
  ggtitle("C  APOE \u2192 TREM2\u2013TYROBP (CellChat)") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 6.0),
        axis.text.y = element_text(size = 6.0),
        legend.position = "bottom", legend.key.size = unit(2.6, "mm"),
        legend.text = element_text(size = 5.6), legend.title = element_text(size = 5.6))

# ---------------------------------------------------------------- G: APOE 共表达伙伴（4 个关键细胞类型）
CT_SEL <- c("Microglia", "Endothelial", "RPE", "Muller Glia")
co <- rbindlist(lapply(CT_SEL, function(cx) {
  f <- file.path(MD, paste0("apoe_coexpression_official_", gsub(" ", "_", cx), ".csv"))   # 2026-09-10: 官方 ppcor（原 apoe_coexpression_adj_ 为自写实现，公式有误）
  if (!file.exists(f)) return(NULL)
  d <- fread(f)
  d <- d[order(-abs(partial_r))][1:6]
  d[, cell_type := cx]
  d
}))
# 2026-09-10 重排：原设计 y 轴 = 19 个基因（4 个细胞类型 top6 去重后），却只分到 ~33 mm 行高
#   → 基因名与瓦片内数值互相重叠，不可投稿。改为**转置**：x = 基因(19) × y = 细胞类型(4)，
#   行数由 19 降到 4，每格宽度 ~8 mm、高度 ~7 mm，标签可读。
gord <- co[, .(mx = max(abs(partial_r))), by = gene][order(-mx)]$gene
co[, gene := factor(gene, levels = gord)]
pG <- ggplot(co, aes(x = gene, y = factor(cell_type, levels = CT_SEL), fill = partial_r)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%+.2f", partial_r)), size = 2.0, colour = "#222222") +
  scale_fill_gradient2(low = col_ad, mid = "white", high = col_amd, midpoint = 0,
                       name = "partial r\n(APOE co-expression)") +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  labs(x = NULL, y = NULL) +
  ggtitle("G  APOE co-expression partners differ by cell type (top 6 per type, partial correlation)") +
  theme_pub +
  theme(axis.text.x = element_text(size = 5.4, face = "italic"),
        axis.text.y = element_text(size = 6.0),
        legend.position = "right", legend.key.size = unit(2.6, "mm"),
        legend.text = element_text(size = 5.6), legend.title = element_text(size = 5.6))


# ---------------------------------------------------------------- H: 特异性与边界（三合一，2026-09-10 PI 要求融合）
# H1 19 蛋白面板：只有 APOE 在两病方向相反且双侧显著（回应"为什么是 APOE"）
fam <- fread(file.path(PROJ, "NC_package/supplementary/TableS8_apoe_family_negcontrol.csv"))
# 每个蛋白×结局只取一个估计：优先 IVW，其次 Wald(best)（单工具蛋白）
fam_pref <- fam[Status == "OK"][order(match(Method, c("IVW", "Wald(best)", "WMed", "Egger")))]
fam_pref <- fam_pref[, .SD[1], by = .(Protein, Outcome)]
fw <- dcast(fam_pref, Protein ~ Outcome, value.var = c("Beta", "P"))
setnames(fw, c("Protein", "b_AD", "b_AMD", "p_AD", "p_AMD"))
fw[, sig_both := (p_AD < 0.05 & p_AMD < 0.05)]
fw[, opp := (b_AD * b_AMD < 0)]
fw[, cat := fcase(Protein == "APOE", "APOE", sig_both & opp, "both + opposite", sig_both, "both significant", default = "other")]
# H2 数据（上次替换时被误删，此处补回）
sd_w <- fread(file.path(PROJ, "revision/results/scdrs_apoe/scdrs_celltype_tests_withAPOEregion.csv"))
sd_n <- fread(file.path(PROJ, "revision/results/scdrs_apoe/scdrs_celltype_tests_noAPOEregion.csv"))
sdt <- rbind(sd_w[set == "withAPOEregion" & disease == "AD"][, set := "with APOE region"],
            sd_n[set == "noAPOEregion"  & disease == "AD"][, set := "without APOE region"])
sdt <- sdt[cell_type %in% CT_ORDER]

pHH1 <- ggplot(sdt, aes(x = factor(cell_type, levels = CT_ORDER), y = rank_biserial, fill = set)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68, colour = "white", linewidth = 0.2) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  scale_fill_manual(values = c("with APOE region" = "#B9B9B9", "without APOE region" = col_ad), name = NULL) +
  labs(x = NULL, y = "scDRS rank-biserial (AD)") +
  ggtitle("H1  scDRS w/o APOE region") +
  theme_pub + theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 5.7),
                    legend.position = "bottom", legend.text = element_text(size = 5.2),
                    legend.key.size = unit(2.4, "mm"), legend.margin = margin(-3, 0, -2, 0))

# H3 区域 vs 全基因组：对立是位点局限的（回应"这是不是全基因组现象"）
# H2（2026-09-10 方案 A）：官方 LAVA 局部遗传相关 —— 不再用自写 Pearson 的 15.4× 比值
lw <- fread(file.path(PROJ, "revision/results/mechanism/lava_chr19_windows.csv"))
lw <- lw[!is.na(rho)]
lw[, mid := (START + STOP) / 2 / 1e6]
lw[, is_apoe := LOC == "chr19_22"]
lw[, bnd := as.logical(boundary)]
oth <- lw[is_apoe == FALSE]; apo <- lw[is_apoe == TRUE]
pHH2 <- ggplot() +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "#999999") +
  geom_errorbar(data = oth, aes(x = mid, ymin = lo, ymax = hi), width = 0, linewidth = 0.22, colour = "#CCCCCC") +
  geom_point(data = oth, aes(x = mid, y = rho, shape = bnd), colour = "#8C8C8C", size = 0.9) +
  geom_errorbar(data = apo, aes(x = mid, ymin = lo, ymax = hi), width = 0, linewidth = 0.55, colour = "#7B4FA8") +
  geom_point(data = apo, aes(x = mid, y = rho), colour = "#7B4FA8", size = 2.3) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1), name = NULL,
                     labels = c("estimable", "CI touches +/-1")) +
  scale_x_continuous(breaks = c(0, 20, 40, 60)) +
  labs(x = "chr19 position (Mb)", y = "LAVA local rg (rho)") +
  ggtitle("H2  Locus-restricted (LAVA local genetic correlation)",
          subtitle = "APOE region (44\u201346 Mb): rho = \u22120.65 (95% CI \u22120.74 to \u22120.57), P = 3e\u221241") +
  theme_pub + theme(legend.position = "bottom", legend.key.size = unit(2.4, "mm"),
                    legend.text = element_text(size = 5.2), axis.text.x = element_text(size = 5.6),
                    plot.subtitle = element_text(size = 5.4, colour = "#7B4FA8"))
  theme_pub + theme(axis.text.x = element_text(size = 5.9),
                    legend.position = "none")

# ---------------------------------------------------------------- 组装
w_mm <- 180; h_mm <- 240   # 2026-09-10：由 300 压到 240 mm（期刊常规上限 240–247）
pA <- pA + theme(plot.margin = margin(4, 6, 2, 4, "mm"))   # figgate 曾报顶边截断，首行单独留白
grid <- pA / (pE | pF) / pB / (pC | pD) / pG / (pHH1 | pHH2) +
  plot_layout(heights = c(1.50, 0.95, 0.72, 0.92, 0.88, 0.68))
w <- w_mm/25.4; h <- h_mm/25.4
ggsave(file.path(PROJ, "figures/Figure5_mechanism.png"), grid, width = w, height = h, dpi = 300)
grDevices::cairo_pdf(file.path(PROJ, "figures/Figure5_mechanism.pdf"), width = w, height = h, family = "Arial"); print(grid); dev.off()
cat("Figure 5 mechanism v4 written (A-G + H1-H2 scope composite; specificity panel removed as redundant with panel F).\n")
