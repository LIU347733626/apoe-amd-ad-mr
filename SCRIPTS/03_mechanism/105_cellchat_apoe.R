## 105_cellchat_apoe.R — 细胞通讯（CellChat）中 APOE 相关信号（GSE230348, 6 细胞类型）
suppressMessages({ library(Matrix); library(CellChat); library(data.table) })
PROJ <- "/home/huan/projects/08_APOE_AMD_AD_MR"
D    <- file.path(PROJ, "revision/tmp_scpagwas/gse230348_10x")
OUT  <- file.path(PROJ, "revision/results/mechanism")
genes <- fread(file.path(D, "genes.tsv"), header = FALSE)$V2
bc    <- fread(file.path(D, "barcodes.tsv"), header = FALSE)$V1
ctl   <- fread(file.path(D, "celltype.tsv"), header = TRUE)
M <- as(readMM(file.path(D, "matrix.mtx")), "CsparseMatrix")
rownames(M) <- make.unique(genes); colnames(M) <- bc
keep <- colnames(M) %in% ctl$barcode
M <- M[, keep, drop = FALSE]
labels <- ctl$cell_type[match(colnames(M), ctl$barcode)]
names(labels) <- colnames(M)
cat("cells:", ncol(M), "| types:", paste(sort(unique(labels)), collapse=", "), "\n")
# 归一化（CellChat 需要 normalized data）
M <- CellChat::normalizeData(M, do.log = TRUE, do.sparse = TRUE)  # 保持稀疏，避免 24 GB 稠密矩阵
cellchat <- createCellChat(object = M, meta = data.frame(labels = labels, row.names = colnames(M)),
                           group.by = "labels")
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cat("\n=== APOE 作为配体/受体出现在哪些 LR 对 ===\n")
lr <- as.data.table(cellchat@LR$LRsig)
ap <- lr[grepl("APOE", ligand) | grepl("APOE", receptor)]
print(ap[, .(interaction_name, ligand, receptor, pathway_name)])
fwrite(ap, file.path(OUT, "cellchat_apoe_LR.csv"))
cat("\n=== 各细胞类型间 APOE 相关信号强度 ===\n")
saveRDS(cellchat, file.path(OUT, "cellchat_gse230348.rds"))
## 交互层面通讯矩阵（net$prob: source x target x interaction）
prob <- cellchat@net$prob
for (i in seq_len(nrow(ap))) {
  iname <- ap$interaction_name[i]
  if (iname %in% dimnames(prob)[[3]]) {
    m <- prob[, , iname]
    cat("\n--", iname, "(", ap$pathway_name[i], ") --\n")
    print(round(m, 5))
    fwrite(as.data.table(as.table(m)), file.path(OUT, paste0("cellchat_", gsub("[^A-Za-z0-9]", "_", iname), ".csv")))
  } else {
    cat("\n--", iname, ": 不在交互矩阵中（可能被 filterCommunication 过滤）--\n")
  }
}
## 通路层面（ApoE 通路）
pp <- cellchat@netP$prob
if ("ApoE" %in% dimnames(pp)[[3]]) {
  cat("\n-- pathway ApoE --\n"); print(round(pp[, , "ApoE"], 5))
  fwrite(as.data.table(as.table(pp[, , "ApoE"])), file.path(OUT, "cellchat_pathway_ApoE.csv"))
}
cat("\nwrote cellchat outputs\n")
