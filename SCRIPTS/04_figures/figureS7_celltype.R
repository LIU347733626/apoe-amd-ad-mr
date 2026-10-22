#!/usr/bin/env Rscript
# figureS7_celltype.R — Supplementary Fig. 7（Cell-type resolution）重出至 180 mm
#
# 背景（2026-09-10 可复现性审计）：
#   原 figures/FigureS7_celltype.{pdf,png} 宽 228.6 mm（= 9.0 inch），超出 180 mm 上限。
#   它由 final_analysis/figure_tissue_cellular.R 末尾的 `save2(figC, "FigureS8_celltype", 9.0, 4.2)`
#   生成（故"无生成脚本"的旧档案记载不准确），内容 = 该脚本的 C|D 两面板。
#
# 本脚本：**逐字复制** figure_tissue_cellular.R 的 C / D 面板代码（逻辑、主题、数据源全不改），
#         仅把画布宽度由 9.0 in（228.6 mm）改为 7.087 in（180.0 mm）、输出文件名为 FigureS7_celltype。
#         因此图注无需改动。
# 数据源（mtime 08-31，早于原图 09-09，未改动）：
#   C: revision/results/scdrs_apoe/scdrs_celltype_summary.csv
#   D: revision/results/seismic_apoe/seismic_{amd,ad}_ct_associations.csv
# 输出：figures/FigureS7_celltype.{png,pdf}
source("/home/huan/software/plots/theme_journal.R")
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})

PROJ <- "/home/huan/projects/08_APOE_AMD_AD_MR"
R    <- file.path(PROJ, "revision", "results")
CT5  <- c("Muller Glia", "Microglia", "Rod", "RPE", "Endothelial")

## ---------- Panel C: scDRS cell-type relevance（复制自 figure_tissue_cellular.R:78-100）----------
s4 <- fread(file.path(R, "scdrs_apoe", "scdrs_celltype_summary.csv"))
s4 <- s4[cell_type %in% CT5]
ctd <- rbindlist(list(
  s4[, .(cell_type, mz = mean_z_AMD, sdz = sd_z_AMD, disease = "AMD")],
  s4[, .(cell_type, mz = mean_z_AD,  sdz = sd_z_AD,  disease = "AD")]))
ctd[, y := match(cell_type, CT5)]
ctd[, sig := ifelse(disease == "AMD", "**",
             ifelse(cell_type %in% c("RPE", "Endothelial", "Muller Glia"), "**",
                    ifelse(cell_type == "Microglia", "*", "")))]

pC <- ggplot(ctd, aes(x = mz, y = y, colour = disease)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = COL_NEUTRAL) +
  geom_errorbarh(aes(xmin = mz - sdz, xmax = mz + sdz), height = 0.22, linewidth = 0.4, alpha = 0.75) +
  geom_point(size = 2.6, shape = 16) +
  geom_text(aes(label = sig), hjust = -0.35, size = 3.4) +
  facet_grid(~ disease) +
  scale_colour_manual(values = c(AMD = COL_AMD, AD = COL_AD)) +
  scale_y_continuous(breaks = c("Muller Glia" = 1, "Microglia" = 2, "Rod" = 3, "RPE" = 4, "Endothelial" = 5),
                     limits = c(0.4, 5.6)) +
  labs(x = "scDRS mean z (higher = more disease-relevant)", y = NULL, disease = "Disease") +
  pub_theme() + guides(colour = "none")

## ---------- Panel D: seismic FDR（复制自 figure_tissue_cellular.R:101-118）----------
seis <- rbindlist(list(
  fread(file.path(R, "seismic_apoe", "seismic_amd_ct_associations.csv")),
  fread(file.path(R, "seismic_apoe", "seismic_ad_ct_associations.csv"))))
seis <- seis[cell_type %in% CT5]
seis[, disease := ifelse(rownames(seis) %in% 1:5, "AMD", "AD")]
seis[, ldfdr := -log10(pmax(FDR, 1e-300))]
seis[, y := match(cell_type, CT5)]

pD <- ggplot(seis, aes(x = ldfdr, y = y, colour = disease)) +
  geom_vline(xintercept = -log10(0.05), linetype = "dotted", colour = COL_NEUTRAL) +
  geom_point(size = 3, shape = 16) +
  facet_grid(~ disease) +
  scale_colour_manual(values = c(AMD = COL_AMD, AD = COL_AD)) +
  scale_y_continuous(breaks = c("Muller Glia" = 1, "Microglia" = 2, "Rod" = 3, "RPE" = 4, "Endothelial" = 5),
                     limits = c(0.4, 5.6)) +
  labs(x = "-log10(FDR) (cell type x trait)", y = NULL, disease = "Disease") +
  pub_theme() + guides(colour = "none")

## ---------- 拼接与输出（仅此处与原脚本不同：宽度 9.0 in → 7.087 in = 180 mm）----------
figC <- (pC | pD) + plot_annotation(title = "Cell-type resolution (scDRS and seismic)",
                                    subtitle = "*FDR<0.05, **FDR<0.001",
                                    theme = theme_void(base_size = 7))

W_IN <- 180 / 25.4   # 7.087 in = 180.0 mm
H_IN <- 4.2          # 高度与原图一致
pub_png(figC, file.path(PROJ, "figures", "FigureS7_celltype.png"), width = W_IN, height = H_IN, dpi = 300)
pdf(file.path(PROJ, "figures", "FigureS7_celltype.pdf"), width = W_IN, height = H_IN); print(figC); dev.off()
cat(sprintf("FigureS7_celltype 已重出：%.1f mm x %.1f mm (原 228.6 mm)\n", W_IN * 25.4, H_IN * 25.4))
