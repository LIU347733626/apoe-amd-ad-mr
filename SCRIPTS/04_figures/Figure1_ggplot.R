#!/usr/bin/env Rscript
# =============================================================================
# Figure 1 (NC, ggplot2 re-implementation of Figure1_NC_v12)
# v14: three cohort tiers (Discovery / First validation / Second validation)
# Per-gene paired (alphabetic), bold header + full-width underline, AMD/AD
# full-width divider, APOE highlighted (bold + purple diamond), right-aligned
# numeric columns. Reads figures/Figure1_forest_data.csv.
#
# Coordinates: single continuous x in [0,100] (hand-mapped; the forest band is
# log-spaced via xpos(or)=X0+k*log10(or)); y = -pos (pos = row index, 1 = top),
# plain y scale so the top row has the largest y.
# =============================================================================
suppressMessages({library(ggplot2); library(data.table)})

d <- fread("figures/Figure1_forest_data.csv", data.table = FALSE)
d$p <- as.numeric(d$p)
d <- d[order(match(d$disease, c("AMD", "AD")), d$gene,
             match(d$cohort, c("Discovery", "First validation", "Second validation"))), ]

# ---- vertical positions (pos: 1 = top; within-gene 0.95, gene gap 1.05) ----
WITHIN <- 0.95; GAP <- 1.05
d$pos <- NA_real_
cur <- 0.9
for (i in seq_len(nrow(d))) {
  d$pos[i] <- cur
  if (i < nrow(d) && d$gene[i] == d$gene[i + 1] && d$disease[i] == d$disease[i + 1]) {
    cur <- cur + WITHIN          # within a gene pair
  } else {
    cur <- cur + GAP             # to next gene
  }
}
d$y <- -d$pos                      # plain scale: top row has largest y

# ---- text shown once ----
d$dis_show <- ifelse(duplicated(d$disease), "", d$disease)
# gene shown once per (disease,gene) — FIRST row only (handles non-paired rows
# like AD PVR/CD177 which have only a Validation row: show the gene there)
d$gene_show <- ifelse(duplicated(paste(d$disease, d$gene)), "", d$gene)
d$or_ci <- sprintf("%.2f (%.2f, %.2f)", d$or, d$ci_lo, d$ci_hi)
# p format: p<1e-3 -> scientific (2 dp mantissa); else 3 significant digits.
# NOTE: format() decides scientific for the WHOLE vector (one tiny value forces
# e-notation everywhere) -> format per element with sapply.
d$p_show <- sapply(d$p, function(x) ifelse(x < 1e-3, sprintf("%.2e", x),
                                           format(signif(x, 3), trim = TRUE)))

# ---- forest x-position (log scale) ----
xpos <- function(or) 47.34 + 38.7 * log10(or)
d$xf  <- xpos(d$or); d$xlo <- xpos(d$ci_lo); d$xhi <- xpos(d$ci_hi)

# column anchors (x in [0,100]); protein pulled left & cohort pushed right so
# the longest gene name (HLA-DRA, ~15.5 mm) never overlaps the Cohort column
dise_x <- 2.5; prot_x <- 11.0; coh_x <- 22.0
or_x <- 86.3;  pval_x <- 96.2

is_apoe <- d$gene == "APOE"
# FIX(2026-09-10): 颜色只编码**队列**；APOE 的高亮由**形状**（宽菱形）承担。
# 原实现把 APOE 三行统一归为 "APOE" 组 → APOE 的 Discovery/First/Second 三行同色，
# 队列信息在颜色上丢失（且 col_v2「Second validation 绿」定义了却用不上）。
d$grp <- factor(d$cohort, levels = c("Discovery", "First validation", "Second validation"))

# ---- APOE wide-diamond polygons (one group per row so they stay separate) ----
polys <- do.call(rbind, lapply(seq_len(nrow(d))[is_apoe], function(i) {
  w <- 1.9; h <- 0.10
  data.frame(grp = i, cohort = d$cohort[i],
             x = c(d$xf[i] - w/2, d$xf[i], d$xf[i] + w/2, d$xf[i]),
             y = c(d$y[i], d$y[i] + h, d$y[i], d$y[i] - h))
}))
d_non <- d[!is_apoe, ]

# ---- CI end caps (short perpendicular ticks at each end; thicker bars) ----
hcap <- 0.09
caps <- do.call(rbind, lapply(seq_len(nrow(d_non)), function(i) {
  r <- d_non[i, ]
  data.frame(x = c(r$xlo, r$xhi), y = c(r$y - hcap, r$y - hcap),
             yend = c(r$y + hcap, r$y + hcap),
             cohort = r$cohort, xend = c(r$xlo, r$xhi))
}))

# ---- layout constants (y = -pos now) ----
y_header <- -min(d$pos) + 1.05                 # header text (above everything)
y_uline  <- -min(d$pos) + 0.45                 # header underline
div_pos  <- mean(c(max(d$pos[d$disease == "AMD"]), min(d$pos[d$disease == "AD"])))
ydiv     <- -div_pos
xdiv     <- c(2.6, 97.4)
y_axis   <- -(max(d$pos) + 0.85)               # x-axis below last row
tick_x   <- c(xpos(0.5), xpos(1.0), xpos(2.0), xpos(4.0))
tick_lab <- c("0.50", "1.0", "2.0", "4.0")
y_ticks  <- 0.35
y_xlab   <- y_axis - 1.4
y_cap    <- y_xlab - 1.3
y_leg    <- y_cap - 1.15                # manual legend row (below the axis caption)

# manual legend (auto legend-label colouring is unsupported in ggplot2, so we
# draw a reliable legend with per-colour labels ourselves)
# FIX(2026-09-10): 图例拆成两个维度，并补上原先缺失的 "Second validation"：
#   (a) 三个方点 = 三个队列（颜色）
#   (b) 一个菱形   = APOE 标记（形状；颜色随队列）
legdf <- data.frame(
  x = c(12, 34, 56), y = y_leg,
  cohort = factor(c("Discovery", "First validation", "Second validation"),
                  levels = c("Discovery", "First validation", "Second validation")),
  txt = c("Discovery", "First validation", "Second validation"))
legdf_apoe <- data.frame(x = 76, y = y_leg, txt = "APOE")

# Okabe-Ito palette (NC / colour-blind-safe): blue, orange, red-purple
col_d <- "#0072B2"; col_v <- "#E69F00"; col_v2 <- "#009E73"; col_a <- "#CC79A7"
c_purp <- col_a; c_null <- "#555555"

p <- ggplot() +
  annotate("segment", x = xdiv[1], xend = xdiv[2], y = y_uline, yend = y_uline,
           linewidth = 0.55, colour = "black") +
  annotate("segment", x = xdiv[1], xend = xdiv[2], y = ydiv, yend = ydiv,
           linewidth = 0.45, colour = "black") +
  annotate("segment", x = xpos(1), xend = xpos(1), y = y_uline, yend = y_axis,
           linewidth = 0.3, colour = c_null) +
  # CI bars (non-APOE only) — thicker, coloured by cohort, with end caps
  geom_segment(data = d_non, aes(x = xlo, xend = xhi, y = y, yend = y,
                                 colour = cohort), linewidth = 0.6) +
  geom_segment(data = caps, aes(x = x, xend = xend, y = y, yend = yend,
                                colour = cohort), linewidth = 0.6) +
  # square markers (non-APOE) — bigger, coloured by cohort
  geom_point(data = d_non, aes(x = xf, y = y, colour = cohort),
             shape = 15, size = 2.2) +
  # APOE wide-diamond (still purple, added to the same colour legend)
  geom_polygon(data = polys, aes(x = x, y = y, group = as.factor(grp),
                                 colour = cohort, fill = cohort),
               linewidth = 0.6) +
  # left text
  geom_text(data = d, aes(x = dise_x, y = y, label = dis_show), size = 3.5,
            hjust = 0, family = "Arial", fontface = "bold", colour = "black") +
  geom_text(data = d, aes(x = prot_x, y = y, label = gene_show), size = 3.5,
            hjust = 0, family = "Arial",
            fontface = ifelse(is_apoe, "bold", "plain"),
            colour = ifelse(is_apoe, c_purp, "black")) +
  geom_text(data = d, aes(x = coh_x, y = y, label = cohort), size = 3.5,
            hjust = 0, family = "Arial", colour = "black") +
  # right numeric (right-aligned)
  geom_text(data = d, aes(x = or_x, y = y, label = or_ci), size = 3.5,
            hjust = 1, family = "Arial",
            fontface = ifelse(is_apoe, "bold", "plain"),
            colour = ifelse(is_apoe, c_purp, "black")) +
  geom_text(data = d, aes(x = pval_x, y = y, label = p_show), size = 3.5,
            hjust = 1, family = "Arial",
            fontface = ifelse(is_apoe, "bold", "plain"),
            colour = ifelse(is_apoe, c_purp, "black")) +
  # header
  geom_text(aes(x = dise_x, y = y_header, label = "Disease"), size = 3.6,
            hjust = 0, family = "Arial", fontface = "bold") +
  geom_text(aes(x = prot_x, y = y_header, label = "Protein"), size = 3.6,
            hjust = 0, family = "Arial", fontface = "bold") +
  geom_text(aes(x = coh_x, y = y_header, label = "Cohort"), size = 3.6,
            hjust = 0, family = "Arial", fontface = "bold") +
  geom_text(aes(x = or_x, y = y_header, label = "OR (95% CI)"), size = 3.6,
            hjust = 1, family = "Arial", fontface = "bold") +
  geom_text(aes(x = pval_x, y = y_header, label = "P value"), size = 3.6,
            hjust = 1, family = "Arial", fontface = "bold") +
  # x-axis
  annotate("segment", x = tick_x[1] - 0.2, xend = tick_x[4] + 0.9,
           y = y_axis, yend = y_axis, linewidth = 0.4, colour = "black") +
  annotate("segment", x = tick_x[1] - 0.2, xend = tick_x[1] - 0.2,
           y = y_axis, yend = y_axis - y_ticks, linewidth = 0.4, colour = "black") +
  annotate("segment", x = tick_x[4] + 0.9, xend = tick_x[4] + 0.9,
           y = y_axis, yend = y_axis - y_ticks, linewidth = 0.4, colour = "black") +
  annotate("segment", x = tick_x, xend = tick_x, y = y_axis, yend = y_axis - y_ticks,
           linewidth = 0.4, colour = "black") +
  annotate("text", x = tick_x, y = y_xlab, label = tick_lab, size = 3.2,
           family = "Arial", vjust = 1) +
  annotate("text", x = 47.3, y = y_cap,
           label = "OR (95% CI) per unit plasma protein (log scale)", size = 3.2,
           family = "Arial", vjust = 1) +
  # manual marker legend (coloured shapes + coloured labels, bottom)
  geom_point(data = legdf, aes(x = x, y = y, colour = cohort), shape = 15, size = 2.2) +
  geom_text(data = legdf, aes(x = x + 1.6, y = y, label = txt, colour = cohort),
            size = 3.2, family = "Arial", hjust = 0) +
  geom_polygon(data = data.frame(x = c(76 - 0.9, 76, 76 + 0.9, 76),
                                 y = c(y_leg, y_leg + 0.10, y_leg, y_leg - 0.10)),
               aes(x = x, y = y), colour = col_a, fill = col_a, linewidth = 0.6) +
  geom_text(data = legdf_apoe, aes(x = x + 1.6, y = y, label = txt),
            size = 3.2, family = "Arial", hjust = 0, colour = c_purp) +
  scale_x_continuous(limits = c(0, 100), expand = expansion(0, 0)) +
  scale_y_continuous(limits = c(y_leg - 0.6, y_header + 0.6), expand = expansion(0, 0)) +
  scale_colour_manual(values = c(Discovery = col_d, `First validation` = col_v,
                                 `Second validation` = col_v2, APOE = col_a), name = NULL,
                      breaks = c("Discovery", "First validation", "Second validation", "APOE")) +
  scale_shape_manual(values = c(Discovery = 15, `First validation` = 15,
                                `Second validation` = 15, APOE = 18),
                     name = NULL, breaks = c("Discovery", "First validation", "Second validation", "APOE")) +
  scale_fill_manual(values = c(Discovery = col_d, `First validation` = col_v,
                               `Second validation` = col_v2), name = NULL, guide = "none") +
  coord_cartesian(clip = "off") +
  theme_void() +
  # match forestplot's vertical centring (content ~19%..86% of height) + side margins
  theme(plot.margin = margin(41, 4, 8, 4, "mm"), legend.position = "none")

w_in <- 180 / 25.4; h_in <- 220 / 25.4
write_fig <- function(f, dev, res = 300) {
  if (dev == "pdf") { grDevices::cairo_pdf(f, width = w_in, height = h_in); print(p); dev.off() }
  else if (dev == "png") { png(f, width = w_in, height = h_in, units = "in", res = res); print(p); dev.off() }
}
write_fig("figures/Figure1_NC.pdf", "pdf")
write_fig("figures/Figure1_NC.png", "png")
write_fig("figures/Figure1_NC.svg", "svg")

cat("ggplot2 Figure 1 written (rows=", nrow(d),
    ", AMD genes=", length(unique(d$gene[d$disease=="AMD"])),
    ", AD genes=", length(unique(d$gene[d$disease=="AD"])),
    ", div_pos=", round(div_pos, 2), ")\n")
