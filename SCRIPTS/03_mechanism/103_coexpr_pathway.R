## 103_coexpr_pathway.R — APOE 共表达程序的通路富集（按细胞类型）
suppressMessages({ library(data.table); library(fgsea) })
PROJ <- "/home/huan/projects/08_APOE_AMD_AD_MR"
MD   <- file.path(PROJ, "revision/results/mechanism")
GMT  <- file.path(PROJ, "data/pathways/c2.cp.kegg.v2023.1.Hs.symbols.gmt")
pw   <- fgsea::gmtPathways(GMT)
files <- list.files(MD, pattern = "^apoe_coexpression_adj_.*csv$", full.names = TRUE)
res_all <- list()
for (f in files) {
  ct <- sub("^apoe_coexpression_adj_", "", sub("[.]csv$", "", basename(f)))
  d  <- fread(f)
  d  <- d[is.finite(partial_r)]
  stats <- setNames(d$partial_r, d$gene)
  stats <- sort(stats, decreasing = TRUE)
  set.seed(20260909)
  r <- fgsea(pathways = pw, stats = stats, minSize = 10, maxSize = 500, nPermSimple = 10000)
  r[, cell_type := ct]
  res_all[[ct]] <- r
  cat("\n===", ct, "===\n")
  print(r[order(-NES)][1:6, .(pathway, size, NES, pval, padj)])
  print(r[order(NES)][1:4, .(pathway, size, NES, pval, padj)])
}
out <- rbindlist(res_all, fill = TRUE)
fwrite(out, file.path(MD, "apoe_coexpression_pathways.csv"))
cat("\nwrote apoe_coexpression_pathways.csv\n")
# 重点通路在各细胞类型的排名
key <- c("KEGG_COMPLEMENT_AND_COAGULATION_CASCADES", "KEGG_LYSOSOME", "KEGG_OXIDATIVE_PHOSPHORYLATION",
         "KEGG_PEROXISOME", "KEGG_PPAR_SIGNALING_PATHWAY", "KEGG_FOCAL_ADHESION")
cat("\n关键通路 NES（细胞类型 × 通路）:\n")
k <- out[pathway %in% key]
if (nrow(k)) print(dcast(k, pathway ~ cell_type, value.var = "NES"))
