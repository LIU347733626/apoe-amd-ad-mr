##############################################################################
# 57a: Proteome-wide GSMR — IAMDGC late AMD (independent validation scan)
# Exposure: UKB-PPP cis-pQTL (official clumped 10000kb/0.001r2, 2923 proteins)
# Outcome:  IAMDGC late AMD, N=32,339 (Iowa cohorts; no UKB / no FinnGen overlap)
# File:     /home/huan/database/GWAS/AMD/IAMDGC_late.ma (header=T, .ma 8 cols)
#
# CAVEAT (write-up): 7 known direction-flipped variants in CFH/ARMS2 regions
#   (chr1/chr10; rs204878/rs635348/rs707937/rs2230199 etc.) are NOT re-oriented
#   here -> interpret non-APOE hits with caution; APOE (chr19) is unaffected.
#
# SMOKE=1  -> APOE-only cell (_pipeline test), N_CORES forced to 2
# N_CORES env -> overrides parallelism (default 128)
#
# Pipeline identical to final_analysis/_01_proteome_gsmr.R (gsmr 1.0.6,
# 7-arg API, n_ref=503, gwas_thresh=5e-8, nsnps_thresh=2, HEIDI internal).
##############################################################################
library(data.table); library(gsmr)

PROJ  <- "/home/huan/projects/08_APOE_AMD_AD_MR"
DB    <- "/home/huan/database/01_raw"
BIN   <- "/home/huan/software/configs_G19/bin"

SMOKE <- Sys.getenv("SMOKE","0") == "1"

AMD_FILE  <- "/home/huan/database/GWAS/AMD/IAMDGC_late.ma"
PQTL_FILE <- file.path(DB, "pqtl/ukb-ppp/summary/ukb-ppp-2023-5e-08-clumped-10000kb-0.001r2.csv")
BFILE  <- file.path(DB, "reference/1000g/EU1000QC")
PQTL_BFILE <- file.path(PROJ, "tmp/pqtl_ref")
BEDFILE <- PQTL_BFILE
PLINK  <- file.path(BIN, "plink")
OUT    <- file.path(PROJ, if (SMOKE) "pqtl_gsmr_results/iamdgc_smoke"
                     else          "pqtl_gsmr_results/iamdgc")
TMP    <- file.path(PROJ, "tmp", if (SMOKE) "gsmr_iamdgc_smoke" else "gsmr_iamdgc")
N_TOTAL_PROTEINS <- 2923
N_CORES <- if (SMOKE) 2L else as.integer(Sys.getenv("N_CORES", "128"))
PREF <- if (SMOKE) "smk" else "iamdgc"

dir.create(OUT, showWarnings=F, recursive=T)
dir.create(TMP, showWarnings=F, recursive=T)

cat(sprintf("===== pQTL -> IAMDGC late AMD GSMR (%s) =====\n",
            if (SMOKE) "SMOKE: APOE only" else "FULL scan"))

# 1. Load
cat("Loading data...\n"); flush(stdout())
amd <- fread(AMD_FILE, header=T, data.table=T,
  colClasses=c("character","character","character","numeric","numeric","numeric","numeric","integer"),
  nThread=4)
setnames(amd, c("SNP","A1","A2","freq","b","se","p","n"))
amd <- amd[!is.na(b) & se>0, .(SNP, A1, A2, b, se, p)]
cat(sprintf("  IAMDGC: %d SNPs\n", nrow(amd)))

pq <- fread(PQTL_FILE, data.table=T,
  select=c("SNP","effect_allele.exposure","other_allele.exposure",
           "beta.exposure","se.exposure","pval.exposure","gene.exposure","exposure"))
setnames(pq, c("SNP","A1","A2","beta","se","pval","gene","exposure_id"))
cat(sprintf("  pQTL: %d instruments, %d proteins\n", nrow(pq), uniqueN(pq$exposure_id)))

if (!file.exists(paste0(PQTL_BFILE, ".bed"))) {
  cat("Creating subsetted pQTL reference...\n")
  all_snps <- unique(pq$SNP)
  writeLines(all_snps, file.path(TMP, "all_pqtl_snps.txt"))
  system(sprintf("%s --bfile %s --extract %s --make-bed --out %s --threads 4 --silent 2>/dev/null",
                 PLINK, BFILE, file.path(TMP, "all_pqtl_snps.txt"), PQTL_BFILE))
}

# 2. Merge + allele harmonization
cat("Merging + harmonizing alleles...\n")
m <- merge(pq, amd, by="SNP", suffixes=c(".pqtl",".amd"))
m[, beta_harm := ifelse(A1.pqtl == A1.amd & A2.pqtl == A2.amd, beta,
              ifelse(A1.pqtl == A2.amd & A2.pqtl == A1.amd, -beta, NA_real_))]
n_unharmonizable <- sum(is.na(m$beta_harm))
m <- m[!is.na(beta_harm) & se.amd>0 & se.pqtl>0 & !is.na(b) & !is.na(pval), ]
cat(sprintf("  Merged: %d, Proteins: %d (unharmonizable: %d)\n",
            nrow(m), uniqueN(m$exposure_id), n_unharmonizable))

# 3. Per-protein LD (plink on pQTL-only reference)
compute_ld <- function(snps, job_id) {
  snp_file <- file.path(TMP, paste0("snps_", job_id, ".txt"))
  writeLines(snps, snp_file)
  ld_file <- file.path(TMP, paste0("ld_", job_id))
  status <- suppressWarnings(system(
    sprintf("%s --bfile %s --extract %s --r square --write-snplist --out %s --threads 2 --silent 2>/dev/null",
            PLINK, BEDFILE, snp_file, ld_file), ignore.stdout=T, ignore.stderr=T))
  if (status != 0) { unlink(snp_file); return(NULL) }
  ld_mat <- paste0(ld_file, ".ld")
  snp_list <- paste0(ld_file, ".snplist")
  if(!file.exists(ld_mat) || !file.exists(snp_list)) {
    unlink(c(snp_file, ld_mat, snp_list, paste0(ld_file, ".log"), paste0(ld_file, ".nosex")))
    return(NULL)
  }
  ld_snps <- tryCatch(readLines(snp_list), error=function(e) NULL)
  if(is.null(ld_snps) || length(ld_snps)<2) {
    unlink(c(snp_file, ld_mat, snp_list, paste0(ld_file, ".log"), paste0(ld_file, ".nosex")))
    return(NULL)
  }
  ld <- as.matrix(read.table(ld_mat))
  colnames(ld) <- ld_snps; rownames(ld) <- ld_snps
  unlink(c(snp_file, ld_mat, snp_list, paste0(ld_file, ".log"), paste0(ld_file, ".nosex")))
  list(ld=ld, snps=ld_snps)
}

# 4. GSMR per protein (gsmr 1.0.6, 7-arg positional API)
run_one <- function(prot, m) {
  sub <- m[m$exposure_id == prot, ]
  if(nrow(sub) < 2) return(NULL)
  sub <- sub[!is.na(pval) & pval > 0, ]
  if(nrow(sub) < 2) return(NULL)

  ld_res <- compute_ld(sub$SNP, paste0(PREF, "_", gsub("[^A-Za-z0-9]","_",prot)))
  if(is.null(ld_res) || nrow(ld_res$ld) < 2) return(NULL)

  sub <- sub[match(ld_res$snps, sub$SNP), ]
  sub <- sub[!is.na(SNP) & !is.na(beta_harm) & !is.na(b), ]
  if(nrow(sub) < 2) return(NULL)

  gsmr_res <- tryCatch({
    gsmr(sub$beta_harm, sub$se.pqtl, sub$pval, sub$b, sub$se.amd,
         ld_res$ld, ld_res$snps, n_ref=503,
         gwas_thresh=5e-8, nsnps_thresh=2, heidi_outlier_flag=T)
  }, error=function(e) NULL)
  if(is.null(gsmr_res)) return(NULL)

  bxy      <- if(length(gsmr_res$bxy)>0)       gsmr_res$bxy[1]       else NA
  bxy_se   <- if(length(gsmr_res$bxy_se)>0)    gsmr_res$bxy_se[1]    else NA
  bxy_p    <- if(length(gsmr_res$bxy_pval)>0)  gsmr_res$bxy_pval[1]  else NA
  n_used   <- if(!is.null(gsmr_res$used_index)) length(gsmr_res$used_index) else nrow(sub)
  heidi_pleio_snps <- if(!is.null(gsmr_res$pleio_snps)) length(gsmr_res$pleio_snps) else 0L
  heidi_pass_flag  <- ifelse(heidi_pleio_snps == 0, "PASS", "FAIL")
  n_heidi <- nrow(sub) - n_used

  data.table(protein=prot, gene=sub$gene[1], n_snp=nrow(sub),
             n_used=as.integer(n_used), bxy=bxy, bxy_se=bxy_se, bxy_p=bxy_p,
             heidi_pleio_snps=heidi_pleio_snps, heidi_pass=heidi_pass_flag,
             n_heidi_removed=n_heidi)
}

proteins <- unique(m$exposure_id)
if (SMOKE) {
  proteins <- proteins[grepl("APOE", proteins)]
  cat(sprintf("SMOKE: running %d APOE protein(s)\n", length(proteins)))
}
N_CORES <- if (SMOKE) 2L else as.integer(Sys.getenv("N_CORES", "128"))
cat(sprintf("Running GSMR for %d proteins on %d cores...\n", length(proteins), N_CORES))
flush(stdout())
raw <- parallel::mclapply(proteins, run_one, m=m, mc.cores=N_CORES, mc.preschedule=TRUE)
results <- Filter(Negate(is.null), raw)
cat(sprintf("  Completed: %d proteins\n", length(results)))
flush(stdout())
if(length(results)==0) stop("No GSMR results. Check data merge and LD reference.", call.=FALSE)

r <- do.call(rbind, results)
r$p_adj <- p.adjust(r$bxy_p, "bonferroni", N_TOTAL_PROTEINS)
r$p_fdr <- p.adjust(r$bxy_p, "fdr")
r <- r[order(r$bxy_p),]
write.csv(r, file.path(OUT,"pqtl_gsmr_v2_full.csv"), row.names=F)
r_pass <- r[r$heidi_pass != "FAIL", ]
write.csv(r_pass, file.path(OUT,"pqtl_gsmr_v2_heidi_pass.csv"), row.names=F)

cat("=== GSMR Results (IAMDGC late AMD) ===\n")
cat(sprintf("Total tested (input): %d | with result: %d | HEIDI pass: %d / FAIL: %d\n",
            N_TOTAL_PROTEINS, nrow(r),
            sum(r$heidi_pass=="PASS"), sum(r$heidi_pass=="FAIL")))
cat(sprintf("Bonferroni < 0.05 (denom=%d): %d | FDR < 0.05: %d\n",
            N_TOTAL_PROTEINS, sum(r$p_adj<0.05,na.rm=T), sum(r$p_fdr<0.05,na.rm=T)))
cat("\nKnown AMD targets:\n")
for(gene in c("APOE","CFH","VEGFA","CFB","C3","TIMP3","CFHR4","CFHR5","ARMS2")){
  hits <- r[grepl(gene, r$gene, ignore.case=T),]
  if(nrow(hits)>0) {
    for(i in 1:nrow(hits)) cat(sprintf("  %s: bxy=%.4f, p=%.2e, p_bonf=%.2e, n_used=%d\n",
        hits$gene[i], hits$bxy[i], hits$bxy_p[i], hits$p_adj[i], hits$n_used[i]))
  }
}
cat("\n--- Session info ---\n"); sessionInfo(); cat("\nDone.\n")
