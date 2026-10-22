#!/usr/bin/env Rscript
# Figure 3 (NC main): Tissue and cellular resolution of the APOE level effect.
# (A) Tissue contrast: plasma pQTL (reference) vs retina eQTL (EyeGEx/Strunz) -> AMD/AD
# (B) APOE-family negative-control panel (proteome specificity)
# (C) scDRS cell-type disease relevance (mean z +/- sd, FDR stars)
# (D) seismic cell-type-trait FDR (APOE driver annotation)
source("/home/huan/software/plots/theme_journal.R")
library(data.table)
library(ggplot2)
library(patchwork)

PROJ <- "/home/huan/projects/08_APOE_AMD_AD_MR"
R <- file.path(PROJ, "revision", "results")
CT5 <- c("Muller Glia", "Microglia", "Rod", "RPE", "Endothelial")

## ---------- Panel A: tissue contrast ----------
ret <- fread(file.path(R, "retina_tissue_mr.csv"))
ret_ivw <- ret[Method == "IVW" & Outcome %in% c("AMD", "AD")]
## 2026-09-09 修正：血浆参考值改为当前 R13 clumped IVW（原硬编码 0.2092/-0.7627 是 R8 时代值）
foc <- fread(file.path(PROJ, "revision", "results", "r13", "apoe_focus_mr_clumped.csv"))
foc_ivw <- foc[grepl("IVW", Method)]
b_amd <- foc_ivw[grepl("AMD", Method), Beta][1]; se_amd <- foc_ivw[grepl("AMD", Method), SE][1]
b_ad  <- foc_ivw[grepl("AD",  Method), Beta][1]; se_ad  <- foc_ivw[grepl("AD",  Method), SE][1]
cat(sprintf("[panelA] plasma R13 IVW: AMD %+.4f (SE %.4f) | AD %+.4f (SE %.4f)\n", b_amd, se_amd, b_ad, se_ad))
ta <- data.table(
  tissue = c("Plasma (UKB-PPP, N=34k)", "Retina (EyeGEx, n=129)",
             "Retina (Strunz 2020, n=311, F~13)", "Plasma (UKB-PPP, N=34k)",
             "Retina (EyeGEx, n=129)", "Retina (Strunz 2020, n=311, F~13)"),
  disease = c("AMD", "AMD", "AMD", "AD", "AD", "AD"),
  b = c(b_amd, NA, NA, b_ad, NA, NA),
  se = c(se_amd, NA, NA, se_ad, NA, NA))
eg <- ret_ivw[InstrumentSet == "EyeGEx_n129"]
st <- ret_ivw[InstrumentSet == "Strunz2020_n311_top10"]
ta[2, b := eg[Outcome == "AMD", Beta]]
ta[2, se := eg[Outcome == "AMD", SE]]
ta[3, b := st[Outcome == "AMD", Beta]]
ta[3, se := st[Outcome == "AMD", SE]]
ta[5, b := eg[Outcome == "AD", Beta]]
ta[5, se := eg[Outcome == "AD", SE]]
ta[6, b := st[Outcome == "AD", Beta]]
ta[6, se := st[Outcome == "AD", SE]]
ta[, grp := ifelse(grepl("^Plasma", tissue), "Plasma", "Retina")]
ta[, y := seq_len(nrow(ta))]

pA <- ggplot(ta, aes(x = b, y = y, colour = grp)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = COL_NEUTRAL) +
  geom_errorbarh(aes(xmin = b - 1.96 * se, xmax = b + 1.96 * se), height = 0.28, linewidth = 0.45) +
  geom_point(size = 2.4, shape = 16) +
  scale_colour_manual(values = c(Plasma = COL_APOE, Retina = COL_NEUTRAL)) +
  scale_y_continuous(breaks = ta$y, labels = ta$tissue, limits = c(0.4, nrow(ta) + 0.8)) +
  annotate("text", x = -0.60, y = nrow(ta) + 0.65, label = "plasma effect >> retinal (systemic level model)",
           hjust = 0, size = 2.6, colour = COL_NEUTRAL) +
  labs(x = "Causal effect of APOE abundance, MR IVW", y = NULL) +
  pub_theme() + guides(colour = "none")

## ---------- Panel B: APOE-family negative controls ----------
nc <- fread(file.path(PROJ, "revision", "results", "r13", "apoe_family_negcontrol.csv"))  # 2026-09-09: R13 面板（原为 R8 时代 revision/results/apoe_family_negcontrol.csv）
nc <- nc[Status == "OK"]
nc1 <- nc[Method == "IVW"]
nc1w <- nc[Method == "Wald(best)"][!paste(Protein, Outcome) %in% paste(nc1$Protein, nc1$Outcome)]
nc1 <- rbindlist(list(nc1, nc1w), fill = TRUE)
comp <- c("C3", "CFD", "CFH", "C9")
nc1[, role := ifelse(Protein == "APOE", "APOE", ifelse(Protein %in% comp, "Complement", "APO-family"))]
prot_order <- unique(nc1[order(Protein)]$Protein)
nc1[, y := match(Protein, prot_order)]

pB <- ggplot(nc1, aes(x = Beta, y = y, colour = role)) +
  facet_grid(~ Outcome) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = COL_NEUTRAL) +
  geom_errorbarh(aes(xmin = Beta - 1.96 * SE, xmax = Beta + 1.96 * SE), height = 0.25, linewidth = 0.4) +
  geom_point(aes(size = ifelse(role == "APOE", 3.2, 2.0)), shape = 16) +
  scale_colour_manual(values = c(APOE = COL_APOE, Complement = "#D4A017", `APO-family` = COL_NEUTRAL)) +
  scale_y_continuous(breaks = seq_along(prot_order), labels = prot_order, limits = c(0.4, length(prot_order) + 0.6)) +
  scale_size(range = c(1.8, 3.4), guide = "none") +
  labs(x = "Protein-level MR (IVW / Wald)", y = NULL, Outcome = "Disease") +
  pub_theme() + guides(colour = "none")

## ---------- Panel C: scDRS cell-type relevance ----------
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

## ---------- Panel D: seismic FDR ----------
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

# v14 (2026-09-09): split into three submission figures
#   FigureS2 = APOE-family specificity (pB); FigureS3 = tissue contrast (pA); FigureS8 = cell type (pC | pD)
figB <- pB + plot_annotation(title = "APOE-family and complement specificity panel",
  theme = theme_void(base_size = 7))
figA <- pA + plot_annotation(title = "Tissue contrast: plasma vs retinal APOE",
  theme = theme_void(base_size = 7))
figC <- (pC | pD) + plot_annotation(title = "Cell-type resolution (scDRS and seismic)",
  subtitle = "*FDR<0.05, **FDR<0.001", theme = theme_void(base_size = 7))

save2 <- function(fig, base, w, h) {
  pub_png(fig, file.path(PROJ, "figures", paste0(base, ".png")), width = w, height = h, dpi = 300)
  pdf(file.path(PROJ, "figures", paste0(base, ".pdf")), width = w, height = h); print(fig); dev.off()
}
save2(figB, "FigureS2_apoe_family_panel", 6.4, 5.0)
save2(figA, "FigureS3_tissue_contrast", 6.4, 4.0)
save2(figC, "FigureS8_celltype", 9.0, 4.2)
cat("FIGURES DONE (S2 / S3 / S8)\n")
