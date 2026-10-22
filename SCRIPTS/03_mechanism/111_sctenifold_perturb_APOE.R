## 111_sctenifold_perturb_APOE.R — APOE 虚拟敲除(KO) + 虚拟过表达(OE)
## 2026-09-09（PI 指出：过表达可用 scTenifold 框架做）
##
## 工具：scTenifoldPerturb 0.1.3（weigenwu/scTenifoldPerturb，GitHub；依赖 CRAN scTenifoldNet 1.4）
##   - KO: perturb_factor = 0  → 目标基因在重建网络中的出边置零
##   - OE: perturb_factor > 1  → 目标基因出边按因子放大（本脚本用 2）
##   同一个 WT 网络复用（scTenifoldPerturbAll），保证 KO/OE 可比。
##   注意：这是**网络层**扰动模型（对重建邻接矩阵的出边做线性缩放），不是表达量扰动模型；
##   过表达因子 2 是模型选择，不是实测的生物学倍数（写进图注/局限）。
##
## 输入：revision/tmp_scpagwas/gse230348_10x（10x mtx，基因×细胞）+ celltype.tsv（原始 counts）
## 基因集：force-include 疾病/APOE 轴基因 + 变异度最高的基因，凑到 VKO_G_TARGET（默认 1000）
##         —— 全基因集(32,849)在 pcNet 上不可行（O(G^2)），见 docs 记录。
## 输出：revision/results/mechanism/sctenifold_{KO,OE}_APOE_<celltype>.csv
##       revision/results/mechanism/sctenifold_APOE_<celltype>_meta.json
suppressMessages({ library(Matrix); library(scTenifoldPerturb); library(data.table) })
PROJ <- "/home/huan/projects/08_APOE_AMD_AD_MR"
D    <- file.path(PROJ, "revision", "tmp_scpagwas", "gse230348_10x")
OUTD <- file.path(PROJ, "revision", "results", "mechanism")
dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)

CT       <- Sys.getenv("VKO_CT", "Microglia")
G_TARGET <- as.integer(Sys.getenv("VKO_G_TARGET", "1000"))
NNET     <- as.integer(Sys.getenv("VKO_NNET", "3"))
NCELLS   <- as.integer(Sys.getenv("VKO_NCELLS", "300"))
NCORES   <- as.integer(Sys.getenv("VKO_CORES", "8"))
OE_FAC   <- as.numeric(Sys.getenv("VKO_OE_FACTOR", "2"))
SEED     <- as.integer(Sys.getenv("VKO_SEED", "20260909"))
SUFFIX   <- Sys.getenv("VKO_OUT_SUFFIX", "")   # 输出文件名后缀（高功率重跑用）
set.seed(SEED)

say <- function(...) { cat(sprintf(...), "\n"); flush.console() }

## ---------- 1. 读原始 counts ----------
genes <- fread(file.path(D, "genes.tsv"), header = FALSE)$V2
bc    <- fread(file.path(D, "barcodes.tsv"), header = FALSE)$V1
ct    <- fread(file.path(D, "celltype.tsv"), header = TRUE)
M     <- as(readMM(file.path(D, "matrix.mtx")), "CsparseMatrix")
rownames(M) <- make.unique(genes); colnames(M) <- bc
say("[1] full matrix: %d genes x %d cells", nrow(M), ncol(M))

keep <- ct$barcode[ct$cell_type == CT]
keep <- intersect(keep, colnames(M))
if (length(keep) < 200) stop("细胞数太少: ", length(keep))
M <- M[, keep, drop = FALSE]
say("[2] %s: %d cells", CT, ncol(M))
stopifnot("APOE" %in% rownames(M))

## ---------- 2. 基因集：force-include + 高变异 ----------
FORCE <- unique(c(
  "APOE", "APOC1", "APOC2", "APOC3", "APOC4", "APOD", "APOH", "APOA1", "APOA2",
  "APOA4", "APOA5", "APOB", "APOL1", "APOL2", "APOL3", "APOL4", "APOL5", "APOL6",
  # 疾病相关蛋白（Figure 5A）
  "CFH", "CFD", "CFI", "CFB", "C2", "C3", "CFHR2", "CFHR4", "CFHR5", "CFHR1", "CFHR3",
  "C1QA", "C1QB", "C1QC", "C4A", "C4B", "C5", "C7", "C9", "CFP", "CD46", "CD55", "CD59",
  "ACE", "EPHA1", "PILRA", "TREM2", "TYROBP", "CR1", "CLU", "PICALM", "BIN1", "MS4A6A",
  "CD33", "ABCA7", "SORL1", "PLCG2", "INPP5D", "SPI1", "MS4A4A", "MS4A4E", "ABI3", "TREM1",
  # APOE 受体 / 脂质转运轴
  "PLAUR", "AXL", "MERTK", "TYRO3", "LRP1", "LRP2", "LRP8", "LDLR", "VLDLR", "SORL1",
  "ABCA1", "ABCG1", "ABCB1", "SCARB1", "AP2A1", "AP2M1", "CLTC", "CUBN", "LRPAP1",
  # 视网膜/AMD 相关
  "ARMS2", "HTRA1", "CETP", "LIPC", "VEGFA", "TIMP3", "COL8A1", "TGFBR1", "VEGFB",
  # 一般参考
  "ACTB", "GAPDH", "B2M", "RPL13A", "TUBB", "GFAP", "RLBP1", "RPE65", "BEST1", "RHO"
))
FORCE <- intersect(FORCE, rownames(M))
say("[3] force-include genes present: %d", length(FORCE))

frac <- Matrix::rowSums(M > 0) / ncol(M)
ok   <- frac >= 0.05                     # 表达 >=5% 细胞（与 scTenifoldNet 的 qc_minCells 精神一致）
cand <- setdiff(rownames(M)[ok], FORCE)
# 变异度（log-normalized 后的方差）
tot <- Matrix::colSums(M); tot[tot == 0] <- 1
Mn  <- M; Mn@x <- log1p(1e4 * Mn@x / rep(tot, diff(Mn@p)))
v   <- sparseMatrixStats::rowVars(Mn[cand, , drop = FALSE])
names(v) <- cand
cand <- cand[order(v, decreasing = TRUE)]
n_fill <- max(0, G_TARGET - length(FORCE))
GENES <- c(FORCE, head(cand, n_fill))
Msub  <- M[GENES, , drop = FALSE]
say("[4] gene set: %d genes (%d forced + %d HVG); cells = %d", length(GENES),
    length(FORCE), n_fill, ncol(Msub))
stopifnot("APOE" %in% rownames(Msub))

## ---------- 3. KO + OE（同一 WT 网络） ----------
t0 <- Sys.time()
res <- scTenifoldPerturbAll(
  countMatrix = Msub, gene = "APOE", modes = c("KO", "OE"),
  perturb_factors = c(KO = 0, OE = OE_FAC),
  preset = "custom", qc = TRUE, qc_minLSize = 500, qc_minCells = 5, qc_mtThreshold = 0.20,
  nc_nNet = NNET, nc_nCells = NCELLS, nc_nComp = 5, nc_lambda = 0.05,
  nc_scaleScores = TRUE, nc_symmetric = FALSE, nc_q = 0.98,
  td_K = 3, td_maxIter = 1000, td_maxError = 1e-5, td_nDecimal = 3, ma_nDim = 3,
  seed = SEED, nCores = NCORES, safe_acceleration = TRUE, memory_efficient = TRUE)
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
say("[5] scTenifoldPerturbAll done in %.1f min", elapsed)

tag <- gsub(" ", "_", CT)
sanity <- list()
for (m in c("KO", "OE")) {
  r <- res[[m]]
  dr <- r$diffRegulation
  dr <- as.data.table(dr)
  dr[, rank := frank(p.value, ties.method = "first")]
  setorder(dr, rank)
  fwrite(dr, file.path(OUTD, sprintf("sctenifold_%s_APOE_%s%s.csv", m, tag, SUFFIX)))
  wt_row <- sum(r$tensorNetworks$WT["APOE", ])
  pe_row <- sum(r$tensorNetworks$Perturbed["APOE", ])
  sanity[[m]] <- list(rows = nrow(dr), wt_rowsum = wt_row, perturbed_rowsum = pe_row,
                      top10 = head(dr$gene, 10),
                      apoe_p = dr[gene == "APOE", p.value])
  say("[6] %s: %d rows; APOE out-edge row sum WT=%.4f -> perturbed=%.4f; APOE p=%s",
      m, nrow(dr), wt_row, pe_row, format(sanity[[m]]$apoe_p, digits = 3))
}

meta <- list(cell_type = CT, n_cells = ncol(Msub), n_genes = nrow(Msub),
             n_forced = length(FORCE), n_hvg = n_fill,
             nNet = NNET, nCells = NCELLS, nComp = 5, td_K = 3, ma_nDim = 3,
             oe_factor = OE_FAC, seed = SEED, elapsed_min = elapsed,
             out_suffix = SUFFIX,
             tool = "scTenifoldPerturb 0.1.3 (weigenwu/scTenifoldPerturb) + scTenifoldNet 1.4",
             sanity = sanity)
jsonlite::write_json(meta, file.path(OUTD, sprintf("sctenifold_APOE_%s%s_meta.json", tag, SUFFIX)),
                     auto_unbox = TRUE, pretty = TRUE)
say("[7] wrote %s and meta json", tag)
say("DONE")
