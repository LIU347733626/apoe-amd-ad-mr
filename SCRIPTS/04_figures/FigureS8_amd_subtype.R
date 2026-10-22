#!/usr/bin/env Rscript
# =============================================================================
# Supplementary Fig. 8 — 干性/湿性 AMD 亚型分层（APOE 水平 MR + ε 等位镜像）
#   A: 工具变量 MR（APOE 丰度 → AMD 亚型）：FinnGen R13 (ALL/DRY/WET) + MVP (ALL/DRY/WET)
#      IVW（圆点）+ 加权中位数（三角）；数值写进 y 轴标签（避免画布溢出）
#   B: 镜像等位基因 per-allele 效应（FinnGen R13）：ε4(rs429358) 保护 / ε2(rs7412) 风险
# 数据：revision/results/amd_subtype/{amd_subtype_meta,amd_subtype_allele_effects}.csv
# =============================================================================
suppressMessages({ library(ggplot2); library(data.table); library(patchwork); library(Cairo) })
PROJ <- "/home/huan/projects/08_APOE_AMD_AD_MR"
MD   <- file.path(PROJ, "revision/results/amd_subtype")
OUT  <- file.path(PROJ, "figures", "FigureS8_amd_subtype")

col_r13 <- "#7B4FA8"; col_mvp <- "#D1403F"
col_e4  <- "#1F6FB2"; col_e2  <- "#C4721A"
SUB_ORD <- c("All AMD", "Dry AMD", "Wet AMD")

theme_pub <- theme_classic(base_size = 7.2, base_family = "Arial") +
  theme(panel.grid = element_blank(),
        axis.line = element_line(linewidth = 0.4, colour = "black"),
        axis.ticks = element_line(linewidth = 0.4, colour = "black"),
        axis.text = element_text(size = 6.4, colour = "#333333"),
        axis.title = element_text(size = 6.8, colour = "#222222"),
        plot.title = element_text(size = 7.6, face = "bold", colour = "#222222", hjust = 0),
        legend.position = "top", legend.title = element_blank(),
        legend.text = element_text(size = 6.2),
        legend.key.size = unit(3.2, "mm"),
        plot.margin = margin(4, 8, 4, 4, "mm"))

# ---------------------------------------------------------------- A: subtype MR
m <- fread(file.path(MD, "amd_subtype_meta.csv"))
m <- m[Method %in% c("IVW", "WMed")]
m[, OR := exp(Beta)]
m[, `:=`(OR_lo = exp(Beta - 1.96 * SE), OR_hi = exp(Beta + 1.96 * SE))]
m[, cohort_lab := fifelse(grepl("FinnGen", cohort), "FinnGen R13", "MVP")]
m[, subtype3 := fifelse(subtype %in% c("ALL", "H7"), "ALL", subtype)]
m[, subtype_lab := factor(subtype3, levels = c("ALL", "DRY", "WET"), labels = SUB_ORD)]
stopifnot(!any(is.na(m$subtype_lab)))

iv <- m[Method == "IVW"]
iv[, est_lab := sprintf("%.2f (%.2f\u2013%.2f), P %s", OR, OR_lo, OR_hi,
                        fifelse(P < 1e-4, sprintf("= %.1e", P), sprintf("= %.3f", P)))]
iv[, row_lab := sprintf("%s | %s\n%s", cohort_lab, subtype_lab, est_lab)]
iv_ord <- iv[order(match(cohort_lab, c("FinnGen R13", "MVP")), match(subtype_lab, SUB_ORD))]
m <- merge(m, iv[, .(cohort_lab, subtype_lab, row_lab)], by = c("cohort_lab", "subtype_lab"))
m[, row_lab := factor(row_lab, levels = rev(iv_ord$row_lab))]

pA <- ggplot(m, aes(x = OR, y = row_lab, colour = cohort_lab, shape = Method)) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = 0.35, colour = "grey55") +
  geom_errorbarh(aes(xmin = OR_lo, xmax = OR_hi), height = 0.2, linewidth = 0.5) +
  geom_point(size = 1.9, fill = "white", stroke = 0.7) +
  scale_colour_manual(values = c("FinnGen R13" = col_r13, "MVP" = col_mvp)) +
  scale_shape_manual(values = c(IVW = 21, WMed = 24)) +
  scale_x_continuous(trans = "log", breaks = c(0.8, 0.9, 1.0, 1.1, 1.2, 1.3),
                     labels = c("0.8", "0.9", "1.0", "1.1", "1.2", "1.3")) +
  coord_cartesian(xlim = c(0.82, 1.33)) +
  labs(title = "A  APOE abundance \u2192 AMD subtype",
       x = "OR per unit genetically predicted plasma APOE", y = NULL) +
  theme_pub +
  theme(axis.text.y = element_text(size = 6.0, lineheight = 1.05, hjust = 1))

# ---------------------------------------------------------------- B: mirror alleles
al <- fread(file.path(MD, "amd_subtype_allele_effects.csv"))
al[, allele := fifelse(SNP == "rs429358", "\u03b54", "\u03b52")]
al[, subtype3 := fifelse(subtype == "H7", "ALL", subtype)]
al[, subtype_lab := factor(subtype3, levels = c("ALL", "DRY", "WET"), labels = SUB_ORD)]
stopifnot(!any(is.na(al$subtype_lab)))
al[, est := sprintf("%.2f (%.2f\u2013%.2f)", OR, OR_lo, OR_hi)]
e4 <- al[SNP == "rs429358", .(subtype_lab, e4_est = est)]
e2 <- al[SNP == "rs7412",   .(subtype_lab, e2_est = est)]
lab <- merge(e4, e2, by = "subtype_lab")
lab[, row_lab := sprintf("%s\n\u03b54 %s | \u03b52 %s", subtype_lab, e4_est, e2_est)]
lab_ord <- lab[order(match(subtype_lab, SUB_ORD))]
al <- merge(al, lab[, .(subtype_lab, row_lab)], by = "subtype_lab")
al[, row_lab := factor(row_lab, levels = rev(lab_ord$row_lab))]

pB <- ggplot(al, aes(x = OR, y = row_lab, colour = allele)) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = 0.35, colour = "grey55") +
  geom_errorbarh(aes(xmin = OR_lo, xmax = OR_hi), height = 0.2,
                 linewidth = 0.5, position = position_dodge(width = 0.4)) +
  geom_point(size = 1.9, shape = 21, fill = "white", stroke = 0.7,
             position = position_dodge(width = 0.4)) +
  scale_colour_manual(values = c("\u03b52" = col_e2, "\u03b54" = col_e4)) +
  scale_x_continuous(breaks = c(0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3)) +
  coord_cartesian(xlim = c(0.68, 1.33)) +
  labs(title = "B  Mirror alleles (FinnGen R13)",
       x = "OR per allele", y = NULL) +
  theme_pub +
  theme(axis.text.y = element_text(size = 6.0, lineheight = 1.05, hjust = 1))

fig <- pA / pB + plot_layout(heights = c(1.3, 1))

CairoPDF(paste0(OUT, ".pdf"), width = 180/25.4, height = 104/25.4)
print(fig); dev.off()
CairoPNG(paste0(OUT, ".png"), width = 180/25.4, height = 104/25.4,
         units = "in", dpi = 600, bg = "white")
print(fig); dev.off()
cat("wrote", OUT, "\n")
