#!/usr/bin/env Rscript
# 43: seismicGWAS (Lai 2025 Nat Commun) on GSE230348 (92k scored cells) — AMD & AD.
# Load: python 10x export -> readMM -> SingleCellExperiment (no Seurat/SeuratDiscover needed).
# Pipeline: calc_specificity -> get_ct_trait_associations -> find_inf_genes (check APOE).
suppressMessages({
  library(Matrix)
  library(SingleCellExperiment)
  library(seismicGWAS)
  library(data.table)
})
DIR  <- "/home/huan/projects/08_APOE_AMD_AD_MR/revision/tmp_seismic_10x"
AMD_Z <- "/home/huan/projects/02_AMD_AD_Network_Pharm/results/finn_amd_gene_level_scores.csv"
AD_Z  <- "/home/huan/projects/08_APOE_AMD_AD_MR/revision/results/magma_ad.genes.out"
OUTD  <- "/home/huan/projects/08_APOE_AMD_AD_MR/revision/results/seismic_apoe"
dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)

mat <- as(Matrix::readMM(file.path(DIR, "matrix.mtx")), "CsparseMatrix")
genes <- fread(file.path(DIR, "genes.tsv"), sep = "\t", header = FALSE)
setnames(genes, c("id", "symbol"))
bc    <- fread(file.path(DIR, "barcodes.tsv"), header = FALSE)[, V1]
ct <- fread(file.path(DIR, "celltype.tsv"))
setnames(ct, c("barcode", "cell_type"))
setkey(ct, barcode)
stopifnot(ncol(mat) == length(bc))
colnames(mat) <- bc
rownames(mat) <- genes$symbol
ct2 <- ct[bc]
ct2[, cell_type := fifelse(is.na(cell_type) | cell_type == "", "Unassigned", cell_type)]

mx <- max(head(mat@x, 200000))
if (mx > 100) { cat(sprintf("X max=%.1f -> raw counts, log1p applied (sparse-safe)\n", mx)); mat@x <- log1p(mat@x) } else
              { cat(sprintf("X max=%.3f -> treated as log-normalized\n", mx)) }
rd <- genes[, .(symbol)]
rownames(rd) <- genes$symbol
sce <- SingleCellExperiment(assays = list(logcounts = mat),
                            rowData = rd,
                            colData = ct2)
cat(sprintf("sce: %d cells, %d genes; cell types: %s\n", ncol(sce), nrow(sce),
            paste(head(sort(table(ct2$cell_type)), 12), collapse = ", ")))

sscore <- calc_specificity(sce, assay_name = "logcounts", ct_label_col = "cell_type", min_ct_size = 50)
write.csv(as.data.frame(sscore), file.path(OUTD, "seismic_specificity_scores.csv"), row.names = FALSE)
cat("specificity scores:", nrow(sscore), "rows\n")

mag_amd <- fread(AMD_Z, header = TRUE)
mag_amd <- unique(mag_amd[, .(GENE, ZSTAT = TOP_Z)])
mag_ad  <- fread(AD_Z, sep = " ", header = TRUE)
mag_ad  <- unique(mag_ad[, .(GENE, ZSTAT)])
cat(sprintf("magma sets: AMD=%d genes, AD=%d genes\n", nrow(mag_amd), nrow(mag_ad)))

assoc_amd <- tryCatch(get_ct_trait_associations(sscore, mag_amd), error = function(e) { cat("AMD assoc error:", conditionMessage(e), "\n"); NULL })
assoc_ad  <- tryCatch(get_ct_trait_associations(sscore, mag_ad),  error = function(e) { cat("AD assoc error:", conditionMessage(e), "\n"); NULL })
if (!is.null(assoc_amd)) write.csv(assoc_amd, file.path(OUTD, "seismic_amd_ct_associations.csv"), row.names = FALSE)
if (!is.null(assoc_ad))  write.csv(assoc_ad,  file.path(OUTD, "seismic_ad_ct_associations.csv"),  row.names = FALSE)
if (!is.null(assoc_amd)) print(head(as.data.frame(assoc_amd), 15))
if (!is.null(assoc_ad))  print(head(as.data.frame(assoc_ad), 15))

cts <- intersect(c("RPE", "Muller Glia", "Microglia", "Rod", "Endothelial"), unique(ct2$cell_type))
for (ctn in cts) {
  for (nm in c("AMD", "AD")) {
    mag <- if (nm == "AMD") mag_amd else mag_ad
    inf <- tryCatch(find_inf_genes(ctn, sscore, mag), error = function(e) NULL)
    if (!is.null(inf)) {
      infdf <- as.data.frame(inf)
      write.csv(infdf, file.path(OUTD, sprintf("seismic_inf_genes_%s_%s.csv", ctn, nm)), row.names = FALSE)
      apoe_row <- infdf[grepl("^APOE$", infdf[[1]]), ]
      cat(sprintf("inf_genes %s/%s: %d rows; APOE: %s\n", ctn, nm, nrow(infdf),
                  if (nrow(apoe_row) > 0) "PRESENT" else "absent"))
    } else {
      cat(sprintf("inf_genes %s/%s: FAILED\n", ctn, nm))
    }
  }
}
cat("DONE 43 seismic\n")
