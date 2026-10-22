#!/usr/bin/env Rscript
# =============================================================================
# Figure 4 (NC, 定量, v8 — 2026-09-09) —— 4 行 × AMD|AD (A-F) + G
#   行: ①丰度轴剂量-反应(+因果模型线+外部纯合锚点+外推带) → ②中介 → ③per-unit 三角验证 → ④加性审计(G)
#   v8 变更: 并入原 Figure 5 的锚点预测与加性审计(Figure 5 改为机制图).
#   每个病独立坐标(free scale), 误差棒用满坐标轴; 配色 NC 克制柔和.
#   A/B: 点位统一中性灰, 线按病着色(AMD红/AD蓝), X轴标签短(ε4/ε3/ε2).
# =============================================================================
suppressMessages({library(ggplot2); library(data.table); library(patchwork)})

col_amd <- "#D1403F"; col_ad <- "#1F6FB2"
col_e4 <- "#3B4CC0"; col_e2 <- "#00A087"; col_mr <- "#E08214"
theme_pub <- theme_classic(base_size = 7.2, base_family = "Arial") +
  theme(panel.grid = element_blank(),
        axis.line = element_line(linewidth = 0.4, colour = "black"),
        axis.ticks = element_line(linewidth = 0.4, colour = "black"),
        axis.text = element_text(size = 6.6, colour = "#333333"),
        axis.title = element_text(size = 6.8, colour = "#222222"),
        plot.title = element_text(size = 7.6, face = "bold", colour = "#222222", hjust = 0),
        plot.margin = margin(5, 10, 5, 5, "mm"))

# =============================================================================
# DATA
# =============================================================================
med <- fread("apoe_results/mediation_results_r13.csv")
dis <- fread("revision/results/r13/s4_doseresponse_disease.csv")
pqt <- fread("revision/results/r13/s4_doseresponse_pqtl.csv")
pq <- pqt[exposure == "UKB-PPP"]; setnames(pq, "pqtl_beta", "pb")
pq[, allele := ifelse(allele == "eps4", "ε4", "ε2")]
dis[, allele := ifelse(allele == "eps4", "ε4", "ε2")]
dis[, cohort := fcase(cohort == "IAMDGC" | cohort == "IAMDGC late", "IAMDGC late",
                      cohort == "IAMDGC all", "IAMDGC all",
                      grepl("Jansen", cohort), "Jansen 2019",
                      grepl("Kunkle", cohort), "Kunkle 2019", default = cohort)]

# --- ① 剂量-反应（含外部纯合锚点）---
hom <- fread("revision/results/mechanism/ko_oe_hom_anchors.csv")
hom[, disease := ifelse(disease == "AD", "AD", "AMD")]
hom[, obs_logOR := log(observed_OR)]
hom[, lo_logOR := log(lo95)]; hom[, hi_logOR := log(hi95)]
hom[, alab := ifelse(allele == "\u03b54", "\u03b54/\u03b54", "\u03b52/\u03b52")]

# --- ① 剂量-反应 ---
dose <- merge(dis, pq[, .(allele, pb)], by = "allele", all.x = TRUE)[cohort == "FinnGen R13"]
dose <- rbindlist(list(dose[, .(disease, allele, gw_beta, gw_se, pb)],
                       data.table(allele = "ε3", disease = c("AMD", "AD"), gw_beta = 0, gw_se = 0, pb = 0)),
                  fill = TRUE)
dose[, al := factor(allele, levels = c("ε2", "ε3", "ε4"))]
dose[, dcol := "#00A087"]
dose1 <- function(d) {
  dd  <- dose[disease == d]
  hh  <- hom[disease == d]
  bmr <- if (d == "AMD") 0.182 else -0.942          # R13 IVW per-unit causal slope
  xr  <- c(-2.35, 1.85)
  cv  <- data.table(x = seq(xr[1], xr[2], length.out = 120))
  cv[, y := bmr * x]
  ggplot(dd, aes(x = pb, y = gw_beta)) +
    annotate("rect", xmin = -Inf, xmax = -1.012, ymin = -Inf, ymax = Inf,
             fill = "#D9D9D9", alpha = 0.40) +
    annotate("rect", xmin = 0.745, xmax = Inf, ymin = -Inf, ymax = Inf,
             fill = "#D9D9D9", alpha = 0.40) +
    geom_hline(yintercept = 0, linetype = 2, colour = "#8A8A8A", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = 2, colour = "#8A8A8A", linewidth = 0.4) +
    geom_line(data = cv, aes(x = x, y = y), inherit.aes = FALSE,
              linetype = 2, colour = "#555555", linewidth = 0.55) +
    geom_line(aes(group = 1), linewidth = 0.9, colour = unique(dd$dcol)) +
    geom_point(size = 2.8, colour = "#4A4A4A") +
    geom_text(aes(label = al), vjust = -1.6, size = 2.4, colour = "#333333") +
    geom_errorbar(data = hh, aes(x = dlog_apoe, ymin = lo_logOR, ymax = hi_logOR),
                  width = 0.07, colour = "#7A5195", linewidth = 0.5, inherit.aes = FALSE) +
    geom_point(data = hh, aes(x = dlog_apoe, y = obs_logOR),
               shape = 22, fill = "white", colour = "#7A5195", size = 2.3, inherit.aes = FALSE) +
    geom_text(data = hh, aes(x = dlog_apoe, y = obs_logOR, label = alab),
              vjust = 1.9, size = 2.0, colour = "#7A5195", inherit.aes = FALSE) +
    scale_x_continuous(breaks = c(-1.012, 0, 0.745), labels = c("\u03b54", "\u03b53", "\u03b52"),
                       limits = xr, name = "Per-copy APOE change (\u0394log)") +
    scale_y_continuous(name = "Per-allele log-OR") +
    labs(y = NULL) + theme_pub + ggtitle(d)
}

# --- ② 中介（两个等位基因 × 间接/直接）---
med2 <- fread("revision/results/r13/mediation_both_alleles_r13.csv")
med2[, allele := ifelse(allele == "e4", "\u03b54", "\u03b52")]
med2[, rowlab := paste0(allele, " \u00b7 ", ifelse(TRUE, "", ""))]
md <- rbindlist(list(
  med2[, .(disease, allele, grp = "indirect", estimate = indirect, se = se_indirect)],
  med2[, .(disease, allele, grp = "direct",   estimate = direct,   se = se_direct)]))
md[, row := factor(paste0(allele, " \u00b7 ", ifelse(grp == "indirect", "Indirect (level)", "Direct (residual)")),
                   levels = c("\u03b54 \u00b7 Indirect (level)", "\u03b54 \u00b7 Direct (residual)",
                              "\u03b52 \u00b7 Indirect (level)", "\u03b52 \u00b7 Direct (residual)"))]
md[, lo := estimate - 1.96*se]; md[, hi := estimate + 1.96*se]
md[, sig := ifelse(abs(estimate/se) >= 1.96, "*", "")]
md[, lab := sprintf("%+.2f%s", estimate, sig)]
med1 <- function(d) {
  dd <- md[disease == d]
  ggplot(dd, aes(x = estimate, y = row, colour = grp)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "#8A8A8A", linewidth = 0.4) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.11, linewidth = 0.55) +
    geom_point(size = 2.4) +
    geom_text(aes(x = hi + 0.03, label = lab), hjust = 0, size = 2.0) +
    scale_colour_manual(values = c("indirect" = "#E64B35", "direct" = "#4DBBD5"), guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0.10, 0.26)), name = "log-OR per allele (null = 0)") +
    labs(y = NULL) + theme_pub + theme(axis.text.y = element_text(size = 5.7, face = "bold")) + ggtitle(d)
}

# --- ③ per-unit 三角验证 ---
disa <- merge(dis, pq[, .(allele, pb)], by = "allele", all.x = TRUE)
disa[, wald := gw_beta / pb]; disa[, se := gw_se / abs(pb)]
disa[, lo := wald - 1.96*se]; disa[, hi := wald + 1.96*se]
disa <- disa[cohort %in% c("FinnGen R13", "IAMDGC late", "Kunkle 2019")]
disa[, est := allele]
mr <- data.table(disease = c("AMD", "AD"), cohort = "FinnGen R13", est = "MR (IVW)",
                 wald = c(0.182, -0.942), se = c(0.028, 0.109))
mr[, lo := wald - 1.96*se]; mr[, hi := wald + 1.96*se]
bD <- rbindlist(list(disa[, .(disease, cohort, est, wald, lo, hi)], mr), fill = TRUE)
bD[, rowlab := paste0(cohort, " \u00b7 ", est)]
bD[, qprio := fcase(cohort == "FinnGen R13", 1, cohort == "IAMDGC late", 2,
                    cohort == "Kunkle 2019", 3, default = 4)]
bD[, estprio := fcase(est == "MR (IVW)", 1, est == "\u03b54", 2, est == "\u03b52", 3, default = 4)]
setorder(bD, disease, qprio, estprio)
tri1 <- function(d) {
  dd <- bD[disease == d]
  dd[, rowlab := factor(rowlab, levels = unique(rowlab))]
  ggplot(dd, aes(x = wald, y = rowlab, colour = est)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "#8A8A8A", linewidth = 0.4) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.13, linewidth = 0.55) +
    geom_point(size = 2.4) +
    geom_text(aes(x = hi + 0.04, label = sprintf("%+.2f", wald)), hjust = 0, size = 2.05) +
    scale_colour_manual(values = c("\u03b54" = col_e4, "\u03b52" = col_e2, "MR (IVW)" = col_mr),
                        name = NULL,
                        labels = c("\u03b54" = "\u03b54 Wald", "\u03b52" = "\u03b52 Wald", "MR (IVW)" = "Instrument MR")) +
    scale_x_continuous(expand = expansion(mult = c(0.10, 0.22)), name = "Per-unit APOE causal effect (log-OR, null = 0)") +
    labs(y = NULL, title = d) +
    theme_pub + theme(legend.position = "bottom", legend.text = element_text(size = 5.6),
                      legend.key.size = unit(0.7, "lines"), legend.margin = margin(-6, 0, 0, 0),
                      axis.text.y = element_text(size = 5.8, face = "bold"),
                      plot.title = element_text(size = 7.4, face = "bold"),
                      plot.subtitle = element_text(size = 5.8, colour = "#555555"))
}

# ---- 3 行 × AMD|AD, 每病独立坐标 ----
grid <- (dose1("AMD") | dose1("AD")) /
        (med1("AMD") | med1("AD")) /
        (tri1("AMD") | tri1("AD"))
grid <- grid + plot_layout(heights = c(1.15, 1.0, 1.0))
grid <- grid + plot_annotation(tag_levels = "A",
                               theme = theme(plot.tag = element_text(size = 9, face = "bold"),
                                             plot.tag.position = c(0, 1)))

w_mm <- 180; h_mm <- 185
w <- w_mm/25.4; h <- h_mm/25.4
svglite::svglite("figures/Figure4_NC.svg", width = w, height = h); print(grid); dev.off()
grDevices::cairo_pdf("figures/Figure4_NC.pdf", width = w, height = h, family = "Arial"); print(grid); dev.off()
ragg::agg_tiff("figures/Figure4_NC.tiff", width = w, height = h, units = "in", res = 600); print(grid); dev.off()
png("figures/Figure4_NC.png", width = w, height = h, units = "in", res = 300); print(grid); dev.off()
cat("Figure 4 v8 (anchored forecast in A/B, both-allele mediation in C/D) done.\n")
