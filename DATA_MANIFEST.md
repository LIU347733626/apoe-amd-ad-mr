# 数据来源与结果清单（FINAL_RELEASE_20260910）

## 一、原始数据（不在包内，路径为唯一真源）

| 用途 | 路径 | 说明 |
|---|---|---|
| 暴露：血浆 APOE cis-pQTL | `/home/huan/database/pQTL/pqtl/ukb-ppp/summary/ukb-ppp-2023-5e-08-unclumped.csv` | UKB-PPP；`gene=="APOE"` 行 = 1,407 个 cis 变异；samplesize 33,995 |
| 暴露工具（**终稿所用**） | `.../ukb-ppp-2023-5e-08-clumped-10000kb-0.001r2.csv` | **官方给定**的 clumped 文件；APOE 得 11 个工具（**不含 rs429358**，故为 ε4-independent 水平工具集） |
| AMD 发现（FinnGen R13 H7） | `smr_data/AMD_r13.ma` | N = 472,750；13,947 例（病例数经 eaf 反解与表型清单双路径核验一致） |
| AD 发现（FinnGen R13 G6） | `smr_data/AD_r13_hdr.ma` | N = 500,186；16,652 例。与 `AD_r13.ma` 去表头后 md5 相同（同一数据） |
| AMD 第一验证 | `/home/huan/database/GWAS/AMD/IAMDGC_late.ma` | N = 32,339（Fritsche 2016） |
| AD 第一验证 | `smr_data/AD_Jansen2019_converted.ma` | N ≈ 436,498（rs429358 行 Nsum）；**尺度与 FinnGen 不同**，仅作方向验证 |
| AD 第二验证 | `smr_data/AD_Bellenguez2022.ma` | GCST90027158；**APOE 区受限**（见 README 限制 5） |
| AMD 第二验证 | `/mnt/ssd4t/GWAS/AMD/MvP/extract/MVP_R4.1000G_AGR.Phe_362_2.META.GIA.dbGaP.txt.gz` | MVP trans-ancestry meta（dbGaP） |
| ε4 剂量-反应外部锚 | `smr_data/Kunkle_2019_AD.txt` | 与 FinnGen 独立方向检查（10/10 位点方向一致） |
| LD 参考 | `/home/huan/references/1000g/EU1000QC.{bed,bim,fam}` | 1000 Genomes EUR，503 个体，hg19，8,550,156 变异 |
| 单细胞 | GSE230348（92,439 细胞） | 细胞类型层（scDRS / seismic / CellChat / 共表达程序） |

## 二、关键结果文件（在项目 `revision/results/` 下）

| 结果 | 生成脚本 |
|---|---|
| `r13/apoe_focus_mr_clumped.csv`（主 MR） | `01_primary_mr/20_r13.R` |
| `replication_iamdgc_{mr,gsmr}.csv` | `01_primary_mr/24_replication_iamdgc.R` |
| `mvp_ancestry/mvp_mr_all.csv` | `01_primary_mr/74_mvp_ancestry_mr.R` |
| `r13/second_validation_gsmr.csv`（**本轮新增**） | `01_primary_mr/98_second_validation_gsmr.R` |
| `r13/loo_official11_summary.csv`（**本轮新增**） | `01_primary_mr/96_leaveoneout_official11.R` |
| `r13/coloc_bydisease.csv` / `r13/coloc_kunkle_ad.csv` | `02_coloc_finemap/60_*.R` |
| `jansen/clumped_coloc_ad_r13.csv` | `02_coloc_finemap/59b_clumped_coloc.R` |
| `mechanism/*`（受体轴/共表达/CellChat/扰动） | `03_mechanism/10x_*.{R,py}` |
| `figures/Figure1_forest_data.csv` | `04_figures/build_fig1_retiered.R` |

## 三、数据完整性核验（2026-09-10 实测）

- FinnGen R13 病例/对照：`data/provenance/finngen_r13_phenos_20260909.json` 记 H7_AMD 13,947+458,803 = **472,750** ✓；G6_ALZHEIMER 16,652+483,534 = **500,186** ✓
- `.ma` 相对官方 `.gz` 的转换：行数差**全部**为"无 rsID 行"（893,576 / 893,606），转换规则可复现
- 效应方向：rs429358（ε4）在 AMD 侧 β<0、AD 侧 β>0；rs7412（ε2）逐位镜像 ✓（7 个文件全部符合领域共识）
- Kunkle 2019：**ε4 定义等位是 C 不是 T**，其 "Effect_allele=T, β=−1.2017" 方向正确（10/10 位点与 FinnGen 一致）
