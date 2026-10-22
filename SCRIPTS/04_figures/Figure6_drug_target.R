#!/usr/bin/env Rscript
# =============================================================================
# Figure 6 (药物靶点图, v1 — 2026-09-09) — APOE 作为药物靶点：方向决定的双疾病权衡
#   A: APOE 丰度轴上的疾病风险（AMD 红升 / AD 蓝降），叠加在研/假想干预的位置与预测方向
#   B: 干预 × 方向 × 预测风险变化（表格式条图）
# =============================================================================
suppressMessages({ library(ggplot2); library(data.table); library(patchwork) })
PROJ <- "/home/huan/projects/08_APOE_AMD_AD_MR"
col_amd <- "#D1403F"; col_ad <- "#1F6FB2"
theme_pub <- theme_classic(base_size = 7.2, base_family = "Arial") +
  theme(panel.grid = element_blank(),
        axis.line = element_line(linewidth = 0.4, colour = "black"),
        axis.ticks = element_line(linewidth = 0.4, colour = "black"),
        axis.text = element_text(size = 6.6, colour = "#333333"),
        axis.title = element_text(size = 6.8, colour = "#222222"),
        plot.title = element_text(size = 7.8, face = "bold", colour = "#222222", hjust = 0),
        plot.margin = margin(5, 8, 5, 5, "mm"))

B_AMD <- 0.182; B_AD <- -0.942      # R13 IVW per-unit log-OR

## ---------- A: 丰度轴上的风险曲线 + 干预位置 ----------
xs <- seq(-2.2, 2.2, length.out = 300)
curves <- rbind(
  data.table(x = xs, y = B_AMD * xs, disease = "AMD"),
  data.table(x = xs, y = B_AD  * xs, disease = "AD"))

# 干预：名称, Δlog(APOE), 标签, 类型
iv <- data.table(
  label = c("LX1001\n(APOE2 gene therapy,\nintrathecal)",
            "ALZ-801\n(APOE4 conformation,\nlevel unchanged)",
            "RO7812653\n(APOE ASO, intrathecal,\nphase 1, 2026)"),
  x = c(1.5, 0, -1.5),
  type = c("raising", "neutral", "lowering"))

pA <- ggplot() +
  annotate("rect", xmin = -1.012, xmax = 0.745, ymin = -Inf, ymax = Inf,
           fill = "#F0F0F0", alpha = 0.7) +
  geom_line(data = curves, aes(x = x, y = y, colour = disease), linewidth = 0.9) +
  geom_hline(yintercept = 0, linetype = 2, colour = "#999999", linewidth = 0.35) +
  geom_vline(xintercept = 0, linetype = 2, colour = "#999999", linewidth = 0.35) +
  geom_point(data = iv, aes(x = x, y = 0), shape = 21, size = 3.2, stroke = 0.5,
             fill = c("#7A5195", "#8A8A8A", "#2E8B57")) +
  geom_text(data = iv, aes(x = x, y = 0, label = label), vjust = c(-1.1, 1.6, -1.1),
            size = 2.1, colour = "#333333", lineheight = 0.95) +
  annotate("text", x = -1.012, y = 2.05, label = "observed\nheterozygous\nrange",
           size = 1.9, colour = "#777777", hjust = 0.5) +
  scale_colour_manual(values = c("AMD" = col_amd, "AD" = col_ad), name = NULL) +
  scale_x_continuous(limits = c(-2.2, 2.2), breaks = c(-2, -1.012, 0, 0.745, 2),
                     labels = c("-2", "\u03b54/4", "0", "\u03b52", "+2"),
                     name = "\u0394 log(APOE abundance)") +
  scale_y_continuous(limits = c(-2.1, 2.1), name = "Predicted disease log-OR") +
  theme_pub + theme(legend.position = "top", legend.text = element_text(size = 6.2),
                    legend.margin = margin(-4, 0, -2, 0),
                    plot.subtitle = element_text(size = 5.8, colour = "#555555")) +
  ggtitle("A  The APOE abundance axis: opposite risk directions",
          subtitle = "Intervention positions are schematic (direction only; the magnitude of APOE change is not known for any agent)")

## ---------- B: 干预 → 预测风险变化 ----------
pred <- data.table(
  agent = c("LX1001", "ALZ-801", "APOE-lowering"),
  dlog  = c(1.5, 0, -1.5))
pred[, amd_logor := B_AMD * dlog][, ad_logor := B_AD * dlog]
pb <- melt(pred[, .(agent, AMD = amd_logor, AD = ad_logor)], id.vars = "agent",
           variable.name = "disease", value.name = "logor")
pb[, agent := factor(agent, levels = c("APOE-lowering", "ALZ-801", "LX1001"))]
pB <- ggplot(pb, aes(x = logor, y = agent, fill = disease)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "#999999", linewidth = 0.35) +
  geom_col(position = position_dodge(width = 0.7), width = 0.62) +
  geom_text(aes(label = sprintf("%+.2f", logor)),
            position = position_dodge(width = 0.7), hjust = ifelse(pb$logor > 0, -0.25, 1.25),
            size = 2.0, colour = "#333333") +
  scale_fill_manual(values = c("AMD" = col_amd, "AD" = col_ad), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.28, 0.28)),
                     name = "Predicted change in disease log-OR") +
  labs(y = NULL) +
  theme_pub + theme(legend.position = "top", legend.text = element_text(size = 6.2),
                    legend.margin = margin(-4, 0, -2, 0)) +
  ggtitle("B  Predicted risk change by intervention direction")

grid <- pA / pB + plot_layout(heights = c(1.35, 0.8))
w_mm <- 180; h_mm <- 150
ggsave(file.path(PROJ, "figures/Figure6_drug_target.png"), grid,
       width = w_mm/25.4, height = h_mm/25.4, dpi = 300)
grDevices::cairo_pdf(file.path(PROJ, "figures/Figure6_drug_target.pdf"),
                     width = w_mm/25.4, height = h_mm/25.4, family = "Arial"); print(grid); dev.off()
cat("Figure 6 drug-target figure written.\n")
