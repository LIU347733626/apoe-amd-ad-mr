#!/usr/bin/env Rscript
# =============================================================================
# Figure 3 (v19) — APOE locus: discovery + validation regional tracks,
#                 exon-level gene track, credible-set attribution
#   A 五条轨道共享 x：pQTL(UKB-PPP) / AMD-R13 / AMD-IAMDGC / AD-R13 / AD-Jansen
#     散点按与 rs429358 的 r² 着色（EU1000QC，plink --recode → 剂量 r²）
#     两条阴影带（水平块 / ε4 簇）+ rs429358 虚线
#   B 外显子级基因轨道（GENCODE v19 = GRCh37；APOE 高亮）
#   C 可信集归因点阵：AMD-R13 | AMD-IAMDGC | AD-R13
# 输出 figures/Figure3_v19.{png,pdf}
# =============================================================================
source("/home/huan/software/plots/theme_journal.R")
suppressPackageStartupMessages({library(ggplot2); library(data.table); library(patchwork)})
PROJ <- "/home/huan/projects/08_APOE_AMD_AD_MR"
R <- file.path(PROJ, "revision/results/r13")
W <- file.path(PROJ, "tmp", "fig3_ld"); dir.create(W, showWarnings = FALSE, recursive = TRUE)
XMIN <- 44.90; XMAX <- 45.85
BAND_LEVEL <- c(45.380, 45.400); BAND_E4 <- c(45.404, 45.418); RS429358 <- 45.411941
## ---- NC 风格配色（低饱和、单色阶为主、色盲安全） ----
R2_RAMP   <- c("#CDE7E7", "#9BD3D3", "#66B8B8", "#2E9B9B", "#15676A")  # r2 / PIP 单色阶（青碧）
COL_LEVEL <- "#5B7FA6"   # 水平调控块（雾蓝）
COL_E4CL  <- "#C0603C"   # ε4 编码簇（陶土）

## ---------- 1. 读五条轨道 ----------
read_trk <- function(f, lab) {
  d <- fread(f); d <- d[pos/1e6 >= XMIN & pos/1e6 <= XMAX]
  data.table(SNP = d$SNP, pos = d$pos/1e6, logp = as.numeric(d$log10p),
             z = abs(as.numeric(d$b) / as.numeric(d$se)), track = lab)
}
trk <- rbindlist(list(
  read_trk(file.path(R, "fig3_pqtl_ukb.csv"),          "APOE pQTL (UKB-PPP)"),
  read_trk(file.path(R, "fig3_gwas_AMD_R13.csv"),      "AMD \u00b7 FinnGen R13"),
  read_trk(file.path(R, "fig3_gwas_AMD_IAMDGC.csv"),   "AMD \u00b7 IAMDGC late"),
  read_trk(file.path(R, "fig3_gwas_AD_R13.csv"),       "AD \u00b7 FinnGen R13"),
  read_trk(file.path(R, "fig3_gwas_AD_Jansen.csv"),    "AD \u00b7 Jansen 2019")))
trk[, disease := fifelse(grepl("^AMD", track), "AMD", fifelse(grepl("^AD", track), "AD", NA_character_))]
trk[, tier := fifelse(grepl("R13", track), "Discovery", fifelse(is.na(disease), "Exposure", "Validation"))]

## ---------- 2. r²(rs429358)（EU1000QC，ped 剂量法；41 号脚本已验证） ----------
r2f <- file.path(W, "r2_rs429358.csv")
if (!file.exists(r2f)) {
  snps <- unique(c("rs429358", trk$SNP))
  writeLines(snps, file.path(W, "snps.txt"))
  PLINK <- "/home/huan/software/configs_G19/bin/plink"
  system(sprintf("%s --bfile /home/huan/references/1000g/EU1000QC --extract %s --recode --out %s/ld >/dev/null 2>&1",
                 PLINK, file.path(W, "snps.txt"), W))
  ped <- fread(file.path(W, "ld.ped"), sep = " ", header = FALSE)
  map <- fread(file.path(W, "ld.map"), header = FALSE)
  ids <- map$V2
  geno <- as.matrix(ped[, -(1:6)])
  # PED 每 SNP 两列；按众数碱基编码剂量
  dose <- matrix(NA_real_, nrow = nrow(geno), ncol = length(ids))
  for (j in seq_along(ids)) {
    a <- geno[, 2 * j - 1]; b <- geno[, 2 * j]
    al <- sort(unique(c(a, b))); ref <- al[1]
    dose[, j] <- (a != ref) + (b != ref)
  }
  colnames(dose) <- ids
  idx <- which(ids == "rs429358")
  if (length(idx) != 1) stop("rs429358 not in LD panel")
  v <- dose[, idx]
  r2 <- sapply(seq_along(ids), function(j) {
    if (j == idx) return(1)
    x <- dose[, j]
    if (sd(x) == 0 || sd(v) == 0) return(NA_real_)
    cor(x, v)^2
  })
  fwrite(data.table(SNP = ids, r2_rs429358 = r2), r2f)
  cat("[2] r2 computed for", length(ids), "SNPs\n")
}
r2 <- fread(r2f)
trk <- merge(trk, r2, by = "SNP", all.x = TRUE)
cat("[2] r2 merged; NA:", sum(is.na(trk$r2_rs429358)), "/", nrow(trk), "\n")

# log10p 在 300 处触顶、且触顶区大量变异 |z| 仅差 0.06 → 不标"每轨 lead"，
# 统一标注 index variant rs429358（LocusZoom 惯例）
idx <- trk[SNP == "rs429358"]
lab_idx <- idx[track == "APOE pQTL (UKB-PPP)"][1]

band_layer <- list(
  annotate("rect", xmin = BAND_LEVEL[1], xmax = BAND_LEVEL[2], ymin = -Inf, ymax = Inf, alpha = 0.16, fill = COL_LEVEL),
  annotate("rect", xmin = BAND_E4[1], xmax = BAND_E4[2], ymin = -Inf, ymax = Inf, alpha = 0.16, fill = COL_E4CL),
  geom_vline(xintercept = RS429358, linetype = "dashed", colour = "grey30", linewidth = 0.35))

base_theme <- theme_bw(base_size = 6.8) +
  theme(strip.placement = "outside", strip.background = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(margin = margin(r = -2.5)),
        strip.text.y.left = element_text(angle = 0, hjust = 0, size = 6.2, face = "bold"),
        strip.text.x = element_text(size = 6.6, face = "bold"),
        panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        legend.position = "right", legend.key.height = unit(0.55, "cm"),
        legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(1, "mm"),
        legend.title = element_text(size = 6), legend.text = element_text(size = 6.0),
        plot.margin = margin(2, 4, 0, 0, "mm"))

# pQTL（全宽，暴露参照）
pq <- trk[tier == "Exposure"]
pQ <- ggplot(pq, aes(pos, logp)) + band_layer +
  geom_point(aes(colour = r2_rs429358), size = 0.6, alpha = 1) +
  geom_point(data = lab_idx, aes(x = pos, y = logp), shape = 18, size = 2.2, colour = "black", inherit.aes = FALSE) +
  ggrepel::geom_text_repel(data = lab_idx, aes(x = pos, y = logp, label = "rs429358 (ε4)"), size = 2.1,
                           nudge_x = 0.16, hjust = 0, direction = "y",
                           min.segment.length = 0, segment.size = 0.22, box.padding = 0.4, seed = 2) +
  scale_colour_gradientn(colours = R2_RAMP, limits = c(0, 1), breaks = c(0, 0.5, 1),
                         name = expression(r^2)) +
  scale_x_continuous(limits = c(XMIN, XMAX), breaks = seq(45.0, 45.8, 0.2)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  guides(colour = guide_colourbar(barwidth = unit(0.16, "cm"), barheight = unit(1.0, "cm"),
                                  ticks.colour = "white", frame.colour = NA)) +
  annotate("text", x = XMIN + 0.015, y = Inf, label = "\u2212log10 P", hjust = 0, vjust = 1.6,
           size = 2.1, colour = "grey25") +
  labs(x = NULL, y = NULL, title = "APOE pQTL (UKB-PPP)") +
  base_theme +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.title = element_text(size = 6.8, face = "bold"))

# 四轨：行 = Discovery/Validation，列 = AMD/AD
d4 <- trk[tier != "Exposure"]
d4[, tier := factor(tier, levels = c("Discovery", "Validation"))]
d4[, disease := factor(disease, levels = c("AMD", "AD"))]
# 分面标题写清"发现/验证 + 疾病 + 队列"
d4[, panel := factor(fifelse(
  tier == "Discovery"  & disease == "AMD", "Discovery \u00b7 AMD \u2014 FinnGen R13",
  fifelse(tier == "Discovery"  & disease == "AD",  "Discovery \u00b7 AD \u2014 FinnGen R13",
  fifelse(tier == "Validation" & disease == "AMD", "Validation \u00b7 AMD \u2014 IAMDGC late",
                                                   "Validation \u00b7 AD \u2014 Jansen 2019"))),
  levels = c("Discovery \u00b7 AMD \u2014 FinnGen R13", "Discovery \u00b7 AD \u2014 FinnGen R13",
             "Validation \u00b7 AMD \u2014 IAMDGC late", "Validation \u00b7 AD \u2014 Jansen 2019"))]
pA <- ggplot(d4, aes(pos, logp)) + band_layer +
  geom_point(aes(colour = r2_rs429358), size = 0.6, alpha = 1) +
  geom_point(data = d4[SNP == "rs429358"], aes(x = pos, y = logp), shape = 18, size = 2.2,
             colour = "black", inherit.aes = FALSE) +
  ggrepel::geom_text_repel(data = d4[SNP == "rs429358"],
             aes(x = pos, y = logp, label = "rs429358 (\u03b54)"), size = 1.9, colour = "black",
             nudge_x = 0.16, hjust = 0, direction = "y", min.segment.length = 0, segment.size = 0.22,
             box.padding = 0.3, seed = 3) +
  facet_wrap(~ panel, ncol = 2, scales = "free_y") +
  scale_colour_gradientn(colours = R2_RAMP, limits = c(0, 1), breaks = c(0, 0.5, 1),
                         name = expression(r^2)) +
  scale_x_continuous(limits = c(XMIN, XMAX), breaks = seq(45.0, 45.8, 0.2)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  guides(colour = guide_colourbar(barwidth = unit(0.16, "cm"), barheight = unit(1.0, "cm"),
                                  ticks.colour = "white", frame.colour = NA)) +
  annotate("text", x = XMIN + 0.015, y = Inf, label = "\u2212log10 P", hjust = 0, vjust = 1.6,
           size = 2.1, colour = "grey25") +
  labs(x = "chr19 position (Mb, GRCh37)", y = NULL) +
  base_theme + theme(axis.text.y = element_text(margin = margin(r = 3)))

## ---------- 3. 外显子级基因轨道（v19：链向箭头 + 关键基因分层配色） ----------
GTF <- "/home/huan/software/magma/ref_data/gencode.v19.gtf.gz"
cmd <- sprintf("zcat %s | awk -F'\\t' '$1==\"chr19\" && ($3==\"gene\" || $3==\"exon\") && $4>44900000 && $5<45850000'", GTF)
gt <- fread(cmd = cmd, header = FALSE, sep = "\t")
gt[, gene := sub('.*gene_name "([^"]+)".*', '\\1', V9)]
gt[, biotype := sub('.*gene_type "([^"]+)".*', '\\1', V9)]
keep_g <- gt[V3 == "gene" & biotype == "protein_coding", unique(gene)]
gn <- gt[V3 == "gene" & biotype == "protein_coding" & gene %in% keep_g,
         .(gene, start = V4/1e6, end = V5/1e6, strand = V7)]
ex <- gt[V3 == "exon" & gene %in% keep_g, .(gene, s = V4/1e6, e = V5/1e6)]
setorder(gn, start)
gn[, row := 0L]
for (i in seq_len(nrow(gn))) for (r in 1:4)
  if (!any(gn[row == r, end] > gn$start[i] - 0.030 & gn[row == r, start] < gn$end[i] + 0.030)) { gn$row[i] <- r; break }
ex <- merge(ex, gn[, .(gene, row)], by = "gene")
KEYG <- c("APOE","TOMM40","PVRL2"); KEY2 <- c("TOMM40","PVRL2")
gn[, `:=`(lab = fifelse(gene %in% KEYG, gene, NA_character_),
          is_key = gene %in% KEYG,
          arr_from = fifelse(strand == "+", end - 0.032, start + 0.032),
          arr_to   = fifelse(strand == "+", end, start))]
ex[, `:=`(is_key = gene %in% KEYG, is_apoe = gene == "APOE")]
arr <- arrow(length = unit(0.5, "mm"), type = "open")

pG <- ggplot() +
  annotate("rect", xmin = BAND_LEVEL[1], xmax = BAND_LEVEL[2], ymin = -Inf, ymax = Inf, alpha = 0.16, fill = COL_LEVEL) +
  annotate("rect", xmin = BAND_E4[1], xmax = BAND_E4[2], ymin = -Inf, ymax = Inf, alpha = 0.16, fill = COL_E4CL) +
  geom_vline(xintercept = RS429358, linetype = "dashed", colour = "grey30", linewidth = 0.35) +
  # 其余基因：浅灰
  geom_segment(data = gn[is_key == FALSE], aes(x = start, xend = end, y = -row, yend = -row),
               colour = "grey58", linewidth = 0.30) +
  geom_segment(data = gn[is_key == FALSE], aes(x = arr_from, xend = arr_to, y = -row, yend = -row),
               colour = "grey52", linewidth = 0.30, arrow = arr) +
  # TOMM40 / PVRL2：中灰
  geom_segment(data = gn[gene %in% KEY2], aes(x = start, xend = end, y = -row, yend = -row),
               colour = "grey30", linewidth = 0.38) +
  geom_segment(data = gn[gene %in% KEY2], aes(x = arr_from, xend = arr_to, y = -row, yend = -row),
               colour = "grey30", linewidth = 0.38, arrow = arr) +
  # APOE：紫色
  geom_segment(data = gn[gene == "APOE"], aes(x = start, xend = end, y = -row, yend = -row),
               colour = COL_APOE, linewidth = 0.42) +
  geom_segment(data = gn[gene == "APOE"], aes(x = arr_from, xend = arr_to, y = -row, yend = -row),
               colour = COL_APOE, linewidth = 0.42, arrow = arr) +
  # 外显子：三级配色
  geom_segment(data = ex[is_key == FALSE], aes(x = s, xend = e, y = -row, yend = -row),
               colour = "grey42", linewidth = 0.78) +
  geom_segment(data = ex[is_key == TRUE & is_apoe == FALSE], aes(x = s, xend = e, y = -row, yend = -row),
               colour = "grey15", linewidth = 1.10) +
  geom_segment(data = ex[is_apoe == TRUE], aes(x = s, xend = e, y = -row, yend = -row),
               colour = COL_APOE, linewidth = 1.55) +
  ggrepel::geom_text_repel(data = gn[!is.na(lab) & gene != "APOE"],
                           aes(x = (start + end)/2, y = -row, label = lab),
                           size = 2.15, fontface = "bold", colour = "grey20",
                           nudge_y = 0.55, direction = "y", box.padding = 0.5, point.padding = 0.25,
                           segment.size = 0.22, segment.colour = "grey55",
                           min.segment.length = 0, max.overlaps = 20, seed = 1) +
  geom_text(data = gn[gene == "APOE"],
            aes(x = (start + end)/2, y = -row + 0.52, label = lab),
            size = 2.15, fontface = "bold", colour = COL_APOE, hjust = 0.5) +
  scale_x_continuous(limits = c(XMIN, XMAX), breaks = seq(45.0, 45.8, 0.2)) +
  scale_y_continuous(limits = c(-3.75, 0.80)) +
  labs(x = "chr19 position (Mb, GRCh37)", y = NULL,
       title = "Gene models \u2014 thick bars = exons, arrows = strand (GENCODE v19 / GRCh37)") +
  theme_bw(base_size = 6.8) +
  theme(panel.grid = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        plot.title = element_text(size = 6.4, face = "bold"), plot.margin = margin(0, 4, 2, 4, "mm"))

## ---------- 4. 归因点阵（v19：按位置排序 + 按功能块分组） ----------
KEY <- c("rs6859","rs157580","rs11668327","rs2927468","rs405509","rs769449","rs429358","rs440446","rs10414043")
LVL <- c("rs2927468","rs6859","rs157580","rs11668327")   # 水平调控块
ad_r13 <- fread(file.path(R, "fig3_finemap_AD.csv")); setnames(ad_r13, c("snp","pip","pos","disease","cohort"), skip_absent = TRUE)
am_iam <- fread(file.path(R, "fig3_finemap_AMD.csv")); setnames(am_iam, c("snp","variable","pip","pos","disease","cohort"), skip_absent = TRUE)
am_r13 <- fread(file.path(PROJ, "revision/results/finemap_amd_r13_pip.csv")); setnames(am_r13, c("snp","pip","pos"))
att <- rbindlist(list(
  am_r13[snp %in% KEY, .(snp, pip = as.numeric(pip), col = "AMD\nFinnGen R13")],
  am_iam[snp %in% KEY, .(snp, pip = as.numeric(pip), col = "AMD\nIAMDGC late")],
  ad_r13[snp %in% KEY, .(snp, pip = as.numeric(pip), col = "AD\nFinnGen R13")]), fill = TRUE)
att <- att[!is.na(pip)]
posmap <- unique(rbind(am_r13[snp %in% KEY, .(snp, pos = as.numeric(pos))],
                       am_iam[snp %in% KEY, .(snp, pos = as.numeric(pos))],
                       ad_r13[snp %in% KEY, .(snp, pos = as.numeric(pos))]), by = "snp")
posmap[, pos := pos/1e6]
posmap[, grp := factor(fifelse(snp %in% LVL, "Level block", "\u03b54 cluster"),
                       levels = c("Level block", "\u03b54 cluster"))]
setorder(posmap, pos)
grid_all <- CJ(snp = posmap$snp, col = c("AMD\nFinnGen R13", "AMD\nIAMDGC late", "AD\nFinnGen R13"))
grid_all <- merge(grid_all, posmap[, .(snp, grp)], by = "snp")
att <- merge(grid_all, att, by = c("snp", "col"), all.x = TRUE)
att[, snp := factor(snp, levels = rev(posmap$snp))]
att[, col := factor(col, levels = c("AMD\nFinnGen R13", "AMD\nIAMDGC late", "AD\nFinnGen R13"))]

pC <- ggplot(att, aes(x = col, y = snp)) +
  geom_point(aes(size = pip, fill = pip), shape = 21, colour = "grey25", stroke = 0.2) +
  geom_point(data = att[is.na(pip)], aes(x = col, y = snp), shape = 95, size = 2.4, colour = "grey55") +
  geom_text(data = unique(att[, .(snp, grp)]), aes(x = 0.52, y = snp, label = snp),
            hjust = 0, size = 2.1, colour = "grey20") +
  facet_grid(grp ~ ., scales = "free_y", space = "free_y") +
  scale_size_continuous(range = c(0.5, 2.6), limits = c(0, 1), guide = "none") +
  scale_fill_gradient(low = R2_RAMP[1], high = R2_RAMP[5], limits = c(0, 1), name = "PIP",
                      guide = guide_colourbar(barwidth = unit(0.16, "cm"), barheight = unit(1.0, "cm"),
                                              ticks.colour = "white", frame.colour = NA)) +
  
  labs(x = NULL, y = NULL, title = "Credible-set attribution (SuSiE PIP)") +
  theme_bw(base_size = 6.8) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text.y = element_text(angle = 90, hjust = 0.5, size = 6.4, face = "bold",
                                    colour = c(COL_LEVEL, COL_E4CL)),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 6.0, lineheight = 0.9),
        plot.title = element_text(size = 6.8, face = "bold"),
        legend.position = "right", legend.title = element_text(size = 6),
        legend.text = element_text(size = 6.0),
        plot.margin = margin(2, 4, 4, 4, "mm"))

fig <- pQ / pA / pG / pC + plot_layout(heights = c(1, 2.2, 0.85, 1.25)) +
  plot_annotation(tag_levels = "A", theme = theme(plot.tag = element_text(size = 8, face = "bold")))
w_in <- 180/25.4; h_in <- 175/25.4
png(file.path(PROJ, "figures", "Figure3_v19.png"), width = w_in, height = h_in, units = "in", res = 600)
print(fig); dev.off()
cairo_pdf(file.path(PROJ, "figures", "Figure3_v19.pdf"), width = w_in, height = h_in)
print(fig); dev.off()
cat("Figure3 v19 written\n")
cat("index variant rows:", nrow(idx), "| tracks:", paste(unique(idx$track), collapse="; "), "\n")
cat("genes:", nrow(gn), "exons:", nrow(ex), "| attribution rows:", nrow(att), "\n")
